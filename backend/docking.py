"""Molecular docking with AutoDock Vina.

Vina runs as an external program (the official ``vina`` executable from
https://github.com/forlilab/AutoDock-Vina/releases), because the
``vina`` Python package has no prebuilt Windows wheels. Fairouz looks for it
in this order:

1. The path in the ``FAIROUZ_VINA`` environment variable.
2. ``vina`` on the system PATH.
3. ``vina`` next to the running Python interpreter (e.g. ``.venv/Scripts``).

Every affinity Fairouz reports is read directly from Vina's output file.
A compound that fails preparation or docking is marked as failed with a
reason; it never stops the rest of the batch and never gets a score.

Vina scores are *predicted docking affinities* in kcal/mol: estimates from
a simplified scoring function, not measured binding affinities.
"""

import csv
import json
import logging
import os
import re
import shutil
import subprocess
import sys
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Any

import config
from backend.descriptors import InvalidSmilesError
from backend.ligands import LigandPreparationError, prepare_ligand
from backend.receptor import BindingBox, file_sha256, validate_binding_box

logger = logging.getLogger(__name__)

VINA_RESULT_PATTERN = re.compile(r"^REMARK VINA RESULT:\s+(-?\d+\.\d+)", re.MULTILINE)

RESULT_CSV_COLUMNS = (
    "target_id",
    "compound_id",
    "name",
    "smiles",
    "ligand_prepared",
    "docking_success",
    "best_affinity_kcal_mol",
    "poses_generated",
    "error",
)


class DockingSetupError(Exception):
    """Raised when docking cannot run at all (missing Vina, receptor or bad box)."""


class VinaNotFoundError(DockingSetupError):
    """Raised when the AutoDock Vina executable cannot be found or run."""


class DockingError(Exception):
    """Raised when docking a single ligand fails."""


def configured_binding_box() -> BindingBox:
    """Return the binding box defined in ``config``."""
    return BindingBox(center=config.BINDING_BOX_CENTER, size=config.BINDING_BOX_SIZE)


# --- Locating Vina -------------------------------------------------------------

def _interpreter_dir_candidate() -> Path:
    """Path where Vina would sit next to the Python interpreter."""
    name = "vina.exe" if os.name == "nt" else "vina"
    return Path(sys.executable).parent / name


def find_vina_executable() -> Path:
    """Locate the Vina executable.

    Returns:
        Path to the executable.

    Raises:
        VinaNotFoundError: With installation instructions if Vina is missing.
    """
    explicit = os.environ.get(config.VINA_EXECUTABLE_ENV_VAR)
    if explicit:
        path = Path(explicit)
        if not path.is_file():
            raise VinaNotFoundError(
                f"{config.VINA_EXECUTABLE_ENV_VAR} is set to {path}, but no file exists there."
            )
        return path

    on_path = shutil.which("vina")
    if on_path:
        return Path(on_path)

    candidate = _interpreter_dir_candidate()
    if candidate.is_file():
        return candidate

    raise VinaNotFoundError(
        "AutoDock Vina was not found. Download the executable for your system from "
        "https://github.com/forlilab/AutoDock-Vina/releases, then either put it "
        f"on your PATH as 'vina', save it as {candidate}, or set "
        f"{config.VINA_EXECUTABLE_ENV_VAR} to its full path. See README › Installation."
    )


def vina_version(executable: Path) -> str:
    """Run ``vina --version`` to confirm the executable works.

    Raises:
        VinaNotFoundError: If the executable cannot be run.
    """
    try:
        completed = subprocess.run(
            [str(executable), "--version"], capture_output=True, text=True, timeout=30
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise VinaNotFoundError(f"AutoDock Vina at {executable} could not be run: {error}") from error
    output = (completed.stdout + completed.stderr).strip()
    if completed.returncode != 0 or "Vina" not in output:
        raise VinaNotFoundError(f"{executable} did not respond like AutoDock Vina: {output[:200]}")
    return output.splitlines()[0]


# --- Docking one ligand ------------------------------------------------------------

def build_vina_command(
    executable: Path,
    receptor: Path,
    ligand: Path,
    output: Path,
    box: BindingBox,
) -> list[str]:
    """Build the exact Vina command line used for one ligand."""
    (cx, cy, cz), (sx, sy, sz) = box.center, box.size
    return [
        str(executable),
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--out", str(output),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", str(config.VINA_EXHAUSTIVENESS),
        "--num_modes", str(config.VINA_NUM_POSES),
        "--seed", str(config.RANDOM_SEED),
    ]


def parse_vina_affinities(pdbqt_text: str) -> list[float]:
    """Read the affinity of every pose from Vina's output PDBQT, best first."""
    return [float(value) for value in VINA_RESULT_PATTERN.findall(pdbqt_text)]


def dock_ligand(
    executable: Path,
    receptor: Path,
    ligand: Path,
    output: Path,
    box: BindingBox,
) -> list[float]:
    """Dock one prepared ligand and return Vina's pose affinities (kcal/mol).

    Raises:
        DockingError: If Vina fails, times out or produces no poses.
    """
    command = build_vina_command(executable, receptor, ligand, output, box)
    try:
        completed = subprocess.run(
            command, capture_output=True, text=True, timeout=config.VINA_TIMEOUT_SECONDS
        )
    except subprocess.TimeoutExpired as error:
        raise DockingError(f"Vina timed out after {config.VINA_TIMEOUT_SECONDS} s.") from error
    except OSError as error:
        raise DockingError(f"Vina could not be started: {error}") from error

    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()[-300:]
        raise DockingError(f"Vina exited with code {completed.returncode}: {detail}")
    if not output.is_file():
        raise DockingError("Vina finished but wrote no output file.")

    affinities = parse_vina_affinities(output.read_text(encoding="utf-8"))
    if not affinities:
        raise DockingError("Vina output contains no docked poses.")
    return affinities


def _docking_result(compound: Mapping[str, Any], **fields: Any) -> dict[str, Any]:
    """Build a docking result with every field present."""
    result = {
        "compound_id": compound.get("compound_id", ""),
        "name": compound.get("name", ""),
        "ligand_prepared": False,
        "docking_success": False,
        "best_affinity_kcal_mol": None,
        "poses_generated": 0,
        "pose_file": None,
        "error": None,
    }
    result.update(fields)
    return result


def dock_compound(
    compound: Mapping[str, Any],
    run_dir: Path,
    executable: Path,
    receptor: Path,
    box: BindingBox,
) -> dict[str, Any]:
    """Prepare and dock one compound, capturing any failure in the result.

    Args:
        compound: Record with ``compound_id``, ``name`` and ``smiles``.
        run_dir: Run folder; ligand and pose files are written inside it.
        executable: The Vina executable.
        receptor: The prepared receptor PDBQT.
        box: The docking search box.

    Returns:
        A structured docking result (see module docs). Never raises for a
        per-compound problem.
    """
    compound_id = compound.get("compound_id", "")
    try:
        ligand_text = prepare_ligand(compound.get("smiles", ""))
    except (InvalidSmilesError, LigandPreparationError) as error:
        logger.warning("Ligand preparation failed for %s: %s", compound_id, error)
        return _docking_result(compound, error=f"Ligand preparation failed: {error}")

    ligand_path = run_dir / "ligands" / f"{compound_id}.pdbqt"
    pose_path = run_dir / "poses" / f"{compound_id}_poses.pdbqt"
    ligand_path.write_text(ligand_text, encoding="utf-8")

    try:
        affinities = dock_ligand(executable, receptor, ligand_path, pose_path, box)
    except DockingError as error:
        logger.warning("Docking failed for %s: %s", compound_id, error)
        return _docking_result(compound, ligand_prepared=True, error=f"Docking failed: {error}")

    logger.info("Docked %s: best affinity %.3f kcal/mol", compound_id, affinities[0])
    return _docking_result(
        compound,
        ligand_prepared=True,
        docking_success=True,
        best_affinity_kcal_mol=affinities[0],
        poses_generated=len(affinities),
        pose_file=str(pose_path),
    )


# --- Batch runs and storage ----------------------------------------------------

def create_run_dir(base_dir: Path) -> Path:
    """Create a new timestamped run folder; never reuse an existing one."""
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    for attempt in range(1, 100):
        suffix = "" if attempt == 1 else f"-{attempt}"
        run_dir = base_dir / f"{stamp}_{config.TARGET_ID}{suffix}"
        try:
            run_dir.mkdir(parents=True)
        except FileExistsError:
            continue
        (run_dir / "ligands").mkdir()
        (run_dir / "poses").mkdir()
        return run_dir
    raise DockingSetupError(f"Could not create a unique run folder in {base_dir}.")


def run_metadata(executable_version: str, receptor: Path, box: BindingBox) -> dict[str, Any]:
    """Everything needed to reproduce a docking run."""
    return {
        "target_id": config.TARGET_ID,
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "receptor_file": receptor.name,
        "receptor_sha256": file_sha256(receptor),
        "binding_box": {**asdict(box), "units": "angstrom"},
        "vina": {
            "version": executable_version,
            "exhaustiveness": config.VINA_EXHAUSTIVENESS,
            "num_modes": config.VINA_NUM_POSES,
            "seed": config.RANDOM_SEED,
        },
        "ligand_preparation": {
            "conformer": "RDKit ETKDGv3",
            "force_field": "MMFF94 (UFF fallback)",
            "seed": config.RANDOM_SEED,
            "rdkit": version("rdkit"),
            "meeko": version("meeko"),
        },
    }


def save_docking_run(
    run_dir: Path,
    results: Sequence[Mapping[str, Any]],
    compounds: Sequence[Mapping[str, Any]],
    metadata: Mapping[str, Any],
) -> None:
    """Write ``docking_results.csv`` and ``run_metadata.json`` to the run folder."""
    smiles_by_id = {compound.get("compound_id"): compound.get("smiles", "") for compound in compounds}
    with (run_dir / "docking_results.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=RESULT_CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for result in results:
            writer.writerow(
                {**result, "target_id": config.TARGET_ID, "smiles": smiles_by_id.get(result["compound_id"], "")}
            )
    (run_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def run_docking(
    compounds: Sequence[Mapping[str, Any]],
    output_dir: Path | None = None,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Dock a batch of compounds against the configured target.

    Args:
        compounds: Valid compound records (``compound_id``, ``name``, ``smiles``).
        output_dir: Parent folder for the run. Defaults to ``config.DOCKING_RESULTS_DIR``.
        progress: Optional callback receiving a message per compound.

    Returns:
        ``stage``, ``status``, ``results`` (one per compound, in order),
        ``run_dir`` and ``metadata``.

    Raises:
        DockingSetupError: If Vina, the prepared receptor or the box is
            unusable. Per-compound failures never raise.
    """
    receptor = config.PREPARED_RECEPTOR_FILE
    if not receptor.is_file():
        raise DockingSetupError(
            f"Prepared receptor not found: {receptor}. Run 'python -m backend.receptor' first."
        )
    box = configured_binding_box()
    try:
        validate_binding_box(box)
    except ValueError as error:
        raise DockingSetupError(f"Invalid binding box in config: {error}") from error

    executable = find_vina_executable()
    executable_version = vina_version(executable)
    run_dir = create_run_dir(output_dir or config.DOCKING_RESULTS_DIR)
    logger.info("Docking %d compound(s) with %s into %s", len(compounds), executable_version, run_dir)

    results = []
    for number, compound in enumerate(compounds, start=1):
        if progress:
            progress(f"Docking {number} of {len(compounds)}: {compound.get('name', '')}...")
        results.append(dock_compound(compound, run_dir, executable, receptor, box))

    metadata = run_metadata(executable_version, receptor, box)
    save_docking_run(run_dir, results, compounds, metadata)
    docked = sum(result["docking_success"] for result in results)
    logger.info("Docking finished: %d of %d succeeded.", docked, len(results))

    return {
        "stage": "docking",
        "status": config.STATUS_COMPLETED,
        "results": results,
        "run_dir": str(run_dir),
        "metadata": metadata,
    }
