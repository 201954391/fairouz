"""Tests for target preparation, ligand preparation and docking.

Fast unit tests run by default. They never run Vina: where a docking score
is needed to test the plumbing, a clearly labelled stub stands in for Vina.

Tests marked ``integration`` run the real Vina executable:

    python -m pytest -m integration
"""

import csv
import json
import math
import sys
from pathlib import Path

import pytest

import config
from backend import docking, ligands
from backend.descriptors import InvalidSmilesError
from backend.docking import DockingError, DockingSetupError, VinaNotFoundError
from backend.ligands import LigandPreparationError, prepare_ligand
from backend.receptor import (
    BindingBox,
    calculate_binding_box,
    file_sha256,
    reference_ligand_coordinates,
    validate_binding_box,
)

ETHANOL = {"compound_id": "T-1", "name": "Ethanol", "smiles": "CCO"}
ASPIRIN = {"compound_id": "T-2", "name": "Aspirin", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"}
INVALID = {"compound_id": "T-3", "name": "Broken", "smiles": "C1CC("}

# Test double only: lets unit tests check how results are structured without
# running Vina. These numbers never reach the app.
STUB_AFFINITIES = [-1.0, -0.5]

VINA_OUTPUT_SAMPLE = """MODEL 1
REMARK VINA RESULT:    -5.085      0.000      0.000
ENDMDL
MODEL 2
REMARK VINA RESULT:    -4.900      2.097      4.399
ENDMDL
"""


# --- Target provenance and binding box -------------------------------------------

def test_raw_target_is_unmodified_download() -> None:
    assert file_sha256(config.RAW_TARGET_FILE) == config.RAW_TARGET_SHA256


def test_prepared_receptor_exists_and_is_separate_from_raw() -> None:
    assert config.PREPARED_RECEPTOR_FILE.is_file()
    assert config.PREPARED_RECEPTOR_FILE != config.RAW_TARGET_FILE


def test_configured_box_is_reproducible_from_reference_ligand() -> None:
    coordinates = reference_ligand_coordinates(config.RAW_TARGET_FILE, config.REFERENCE_LIGAND_CHAIN)
    box = calculate_binding_box(coordinates, config.BINDING_BOX_PADDING)

    assert len(coordinates) == 49  # Heavy atoms of the N3 inhibitor in 6LU7.
    assert box.center == pytest.approx(config.BINDING_BOX_CENTER, abs=1e-3)
    assert box.size == pytest.approx(config.BINDING_BOX_SIZE, abs=1e-3)


def test_binding_box_from_known_coordinates() -> None:
    box = calculate_binding_box([(0.0, 0.0, 0.0), (2.0, 4.0, 6.0)], padding=1.0)

    assert box.center == pytest.approx((1.0, 2.0, 3.0))
    assert box.size == pytest.approx((4.0, 6.0, 8.0))


@pytest.mark.parametrize(
    "box",
    [
        BindingBox(center=(0.0, 0.0, 0.0), size=(10.0, -1.0, 10.0)),
        BindingBox(center=(0.0, 0.0, 0.0), size=(10.0, 0.0, 10.0)),
        BindingBox(center=(0.0, math.nan, 0.0), size=(10.0, 10.0, 10.0)),
        BindingBox(center=(0.0, 0.0), size=(10.0, 10.0, 10.0)),
    ],
)
def test_invalid_binding_box_is_rejected(box: BindingBox) -> None:
    with pytest.raises(ValueError):
        validate_binding_box(box)


def test_binding_box_needs_atoms_and_non_negative_padding() -> None:
    with pytest.raises(ValueError):
        calculate_binding_box([], padding=5.0)
    with pytest.raises(ValueError):
        calculate_binding_box([(0.0, 0.0, 0.0)], padding=-1.0)


# --- Ligand preparation ------------------------------------------------------------

def test_prepare_ligand_returns_pdbqt() -> None:
    pdbqt = prepare_ligand(ASPIRIN["smiles"])

    assert "ROOT" in pdbqt and "TORSDOF" in pdbqt
    assert any(line.startswith(("ATOM", "HETATM")) for line in pdbqt.splitlines())


def test_ligand_preparation_is_deterministic() -> None:
    assert prepare_ligand(ASPIRIN["smiles"]) == prepare_ligand(ASPIRIN["smiles"])


def test_invalid_smiles_is_rejected_before_3d() -> None:
    with pytest.raises(InvalidSmilesError):
        prepare_ligand(INVALID["smiles"])


def test_conformer_failure_is_reported(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(ligands.AllChem, "EmbedMolecule", lambda *args, **kwargs: -1)

    with pytest.raises(LigandPreparationError, match="3D conformer"):
        prepare_ligand(ETHANOL["smiles"])


# --- Vina command and output -------------------------------------------------------

def test_parse_vina_affinities_reads_every_pose() -> None:
    assert docking.parse_vina_affinities(VINA_OUTPUT_SAMPLE) == [-5.085, -4.9]


def test_vina_command_uses_configured_settings() -> None:
    box = docking.configured_binding_box()
    command = docking.build_vina_command(Path("vina"), Path("r.pdbqt"), Path("l.pdbqt"), Path("o.pdbqt"), box)

    assert command[command.index("--exhaustiveness") + 1] == str(config.VINA_EXHAUSTIVENESS)
    assert command[command.index("--seed") + 1] == str(config.RANDOM_SEED)
    assert command[command.index("--size_y") + 1] == str(config.BINDING_BOX_SIZE[1])


def test_failing_vina_process_raises_docking_error(tmp_path: Path) -> None:
    # Python stands in for a broken Vina: it rejects Vina's arguments and exits non-zero.
    with pytest.raises(DockingError, match="exited with code"):
        docking.dock_ligand(
            Path(sys.executable),
            tmp_path / "receptor.pdbqt",
            tmp_path / "ligand.pdbqt",
            tmp_path / "out.pdbqt",
            docking.configured_binding_box(),
        )


# --- Setup errors --------------------------------------------------------------------

def test_missing_vina_gives_install_instructions(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.delenv(config.VINA_EXECUTABLE_ENV_VAR, raising=False)
    monkeypatch.setattr(docking.shutil, "which", lambda name: None)
    monkeypatch.setattr(docking, "_interpreter_dir_candidate", lambda: tmp_path / "vina.exe")

    with pytest.raises(VinaNotFoundError, match="releases"):
        docking.find_vina_executable()


def test_bad_vina_env_var_is_reported(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv(config.VINA_EXECUTABLE_ENV_VAR, str(tmp_path / "nope.exe"))

    with pytest.raises(VinaNotFoundError, match=config.VINA_EXECUTABLE_ENV_VAR):
        docking.find_vina_executable()


def test_missing_receptor_is_a_setup_error(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(config, "PREPARED_RECEPTOR_FILE", tmp_path / "missing.pdbqt")

    with pytest.raises(DockingSetupError, match="Prepared receptor not found"):
        docking.run_docking([ETHANOL], output_dir=tmp_path)


# --- Structured output and batches (Vina stubbed) --------------------------------

@pytest.fixture
def stubbed_vina(monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace the Vina process with a stub so batch logic can be tested fast."""
    monkeypatch.setattr(docking, "find_vina_executable", lambda: Path("stub-vina"))
    monkeypatch.setattr(docking, "vina_version", lambda executable: "stub Vina")
    monkeypatch.setattr(docking, "dock_ligand", lambda *args: list(STUB_AFFINITIES))


def test_successful_result_has_expected_structure(stubbed_vina: None, tmp_path: Path) -> None:
    result = docking.run_docking([ETHANOL], output_dir=tmp_path)["results"][0]

    assert result["compound_id"] == "T-1"
    assert result["docking_success"] is True
    assert result["best_affinity_kcal_mol"] == STUB_AFFINITIES[0]
    assert result["poses_generated"] == len(STUB_AFFINITIES)
    assert result["error"] is None


def test_batch_continues_after_one_ligand_fails(stubbed_vina: None, tmp_path: Path) -> None:
    results = docking.run_docking([ETHANOL, INVALID, ASPIRIN], output_dir=tmp_path)["results"]

    assert [result["compound_id"] for result in results] == ["T-1", "T-3", "T-2"]
    assert [result["docking_success"] for result in results] == [True, False, True]
    failed = results[1]
    assert failed["best_affinity_kcal_mol"] is None
    assert failed["ligand_prepared"] is False
    assert "Ligand preparation failed" in failed["error"]


def test_vina_failure_for_one_ligand_is_recorded(monkeypatch: pytest.MonkeyPatch, stubbed_vina: None, tmp_path: Path) -> None:
    def dock_or_fail(executable, receptor, ligand, output, box):
        if ligand.stem == "T-2":
            raise DockingError("stub failure")
        return list(STUB_AFFINITIES)

    monkeypatch.setattr(docking, "dock_ligand", dock_or_fail)
    results = docking.run_docking([ETHANOL, ASPIRIN], output_dir=tmp_path)["results"]

    assert results[1]["ligand_prepared"] is True
    assert results[1]["docking_success"] is False
    assert results[1]["best_affinity_kcal_mol"] is None
    assert "stub failure" in results[1]["error"]


def test_run_is_saved_with_reproducibility_metadata(stubbed_vina: None, tmp_path: Path) -> None:
    run = docking.run_docking([ETHANOL, INVALID], output_dir=tmp_path)
    run_dir = Path(run["run_dir"])

    with (run_dir / "docking_results.csv").open(encoding="utf-8") as file:
        rows = list(csv.DictReader(file))
    metadata = json.loads((run_dir / "run_metadata.json").read_text(encoding="utf-8"))

    assert [row["compound_id"] for row in rows] == ["T-1", "T-3"]
    assert all(row["target_id"] == config.TARGET_ID for row in rows)
    assert metadata["vina"]["seed"] == config.RANDOM_SEED
    assert metadata["binding_box"]["size"] == list(config.BINDING_BOX_SIZE)
    assert metadata["receptor_sha256"] == file_sha256(config.PREPARED_RECEPTOR_FILE)


def test_runs_never_overwrite_each_other(tmp_path: Path) -> None:
    first = docking.create_run_dir(tmp_path)
    second = docking.create_run_dir(tmp_path)

    assert first != second and first.is_dir() and second.is_dir()


# --- Integration: real AutoDock Vina ---------------------------------------------

def _require_vina() -> Path:
    try:
        return docking.find_vina_executable()
    except VinaNotFoundError as error:
        pytest.skip(str(error))


def _first_pose_coordinates(pose_file: Path) -> list[tuple[float, float, float]]:
    coordinates = []
    for line in pose_file.read_text(encoding="utf-8").splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")):
            coordinates.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return coordinates


@pytest.mark.integration
def test_real_docking_returns_vina_scores_inside_box(tmp_path: Path) -> None:
    _require_vina()
    result = docking.run_docking([ASPIRIN], output_dir=tmp_path)["results"][0]

    assert result["docking_success"] is True
    assert math.isfinite(result["best_affinity_kcal_mol"])
    assert result["best_affinity_kcal_mol"] < 0
    assert result["poses_generated"] >= 1

    box = docking.configured_binding_box()
    for point in _first_pose_coordinates(Path(result["pose_file"])):
        for value, center, size in zip(point, box.center, box.size):
            assert abs(value - center) <= size / 2 + 0.5  # 0.5 Å tolerance for hydrogens.


@pytest.mark.integration
def test_real_docking_is_reproducible_with_fixed_seed(tmp_path: Path) -> None:
    _require_vina()
    first = docking.run_docking([ASPIRIN], output_dir=tmp_path)["results"][0]
    second = docking.run_docking([ASPIRIN], output_dir=tmp_path)["results"][0]

    assert first["best_affinity_kcal_mol"] == second["best_affinity_kcal_mol"]


@pytest.mark.integration
def test_real_batch_with_invalid_compound(tmp_path: Path) -> None:
    _require_vina()
    results = docking.run_docking([ETHANOL, INVALID], output_dir=tmp_path)["results"]

    assert results[0]["docking_success"] is True
    assert results[1]["docking_success"] is False
    assert results[1]["best_affinity_kcal_mol"] is None
