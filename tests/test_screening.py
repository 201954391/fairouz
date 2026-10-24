"""Tests for the screening pipeline coordinator.

Docking is replaced by a stub here so these tests stay fast; real docking is
covered by the integration tests in ``test_docking.py``.
"""

from pathlib import Path

import pytest

import config
from backend import screening
from backend.docking import VinaNotFoundError
from backend.screening import run_screening

TARGET = config.AVAILABLE_TARGETS[0]


def _stub_docking(compounds, progress=None):
    """Stand-in for run_docking: marks every compound as not docked."""
    return {
        "stage": "docking",
        "status": config.STATUS_COMPLETED,
        "results": [
            {
                "compound_id": compound["compound_id"],
                "name": compound["name"],
                "ligand_prepared": True,
                "docking_success": False,
                "best_affinity_kcal_mol": None,
                "poses_generated": 0,
                "pose_file": None,
                "error": "stub: docking not run in unit tests",
            }
            for compound in compounds
        ],
        "run_dir": str(Path("stub")),
        "metadata": {},
    }


@pytest.fixture(autouse=True)
def stub_docking(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(screening, "run_docking", _stub_docking)


def test_requesting_more_than_available_uses_all_without_duplicates() -> None:
    response = run_screening(TARGET, config.MAX_COMPOUNDS)
    summary = response["summary"]

    assert summary["processed"] == summary["loaded"] < config.MAX_COMPOUNDS
    ids = [compound["compound_id"] for compound in response["compounds"]]
    assert len(ids) == len(set(ids))
    assert any("were requested" in warning for warning in response["warnings"])


def test_counts_add_up() -> None:
    summary = run_screening(TARGET, config.MIN_COMPOUNDS)["summary"]

    assert summary["valid"] + summary["invalid"] == summary["processed"]
    assert summary["docked"] + summary["docking_failed"] == summary["valid"]


def test_each_compound_carries_its_docking_result() -> None:
    compounds = run_screening(TARGET, config.MIN_COMPOUNDS)["compounds"]

    for compound in compounds:
        assert compound["docking"]["compound_id"] == compound["compound_id"]


def test_stage_statuses() -> None:
    stages = run_screening(TARGET, config.MIN_COMPOUNDS)["stages"]

    assert stages["compound_loading"] == config.STATUS_COMPLETED
    assert stages["molecular_descriptors"] == config.STATUS_COMPLETED
    assert stages["docking"] == config.STATUS_COMPLETED
    assert stages["ai_ranking"] == config.STATUS_NOT_IMPLEMENTED
    assert stages["report"] == config.STATUS_NOT_IMPLEMENTED


def test_missing_vina_keeps_descriptors_and_reports_unavailable(monkeypatch: pytest.MonkeyPatch) -> None:
    def no_vina(compounds, progress=None):
        raise VinaNotFoundError("Vina missing (test)")

    monkeypatch.setattr(screening, "run_docking", no_vina)
    response = run_screening(TARGET, config.MIN_COMPOUNDS)

    assert response["stages"]["docking"] == config.STATUS_UNAVAILABLE
    assert response["summary"]["valid"] > 0
    assert response["summary"]["docked"] == 0
    assert all(compound["docking"] is None for compound in response["compounds"])
    assert "Vina missing" in response["docking"]["error"]


def test_progress_reports_real_stages_and_no_ai_ranking() -> None:
    messages: list[str] = []

    run_screening(TARGET, config.MIN_COMPOUNDS, progress=messages.append)

    assert messages[0] == "Loading compounds..."
    assert any("descriptors" in message for message in messages)
    assert any("3D ligands" in message for message in messages)
    assert messages[-1] == "Complete."
    assert not any("ranking" in message.lower() for message in messages)


@pytest.mark.parametrize(
    ("target", "n_compounds"),
    [("Unknown Protein", 10), (TARGET, config.MIN_COMPOUNDS - 1), (TARGET, config.MAX_COMPOUNDS + 1)],
)
def test_invalid_requests_raise(target: str, n_compounds: int) -> None:
    with pytest.raises(ValueError):
        run_screening(target, n_compounds)
