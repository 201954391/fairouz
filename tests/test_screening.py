"""Tests for the screening pipeline coordinator.

Docking is replaced by a stub in the unit tests so they stay fast. The test
marked ``integration`` runs the whole pipeline with real Vina.
"""

import json
from pathlib import Path

import pytest

import config
from backend import docking, ranking, screening
from backend.docking import VinaNotFoundError
from backend.screening import run_screening

TARGET = config.AVAILABLE_TARGETS[0]

# Test double only: affinities for the stub below, never shown in the app.
STUB_AFFINITY_START = -8.0
STUB_AFFINITY_STEP = 0.5


def _stub_docking_success(compounds, progress=None):
    """Stand-in for run_docking that "docks" every compound with a distinct stub affinity."""
    return {
        "stage": "docking",
        "status": config.STATUS_COMPLETED,
        "results": [
            {
                "compound_id": compound["compound_id"],
                "name": compound["name"],
                "ligand_prepared": True,
                "docking_success": True,
                "best_affinity_kcal_mol": STUB_AFFINITY_START + index * STUB_AFFINITY_STEP,
                "poses_generated": 9,
                "pose_file": None,
                "error": None,
            }
            for index, compound in enumerate(compounds)
        ],
        "run_dir": str(Path("stub")),
        "metadata": {},
    }


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
def stub_docking(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> None:
    if "integration" not in request.keywords:
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


def test_stage_statuses_when_every_docking_fails() -> None:
    response = run_screening(TARGET, config.MIN_COMPOUNDS)  # Stub: ligands prepared, every docking fails.

    assert list(response["stages"]) == list(config.PIPELINE_STAGE_LABELS)
    assert response["stages"] == {
        "compound_loading": config.STATUS_COMPLETED,
        "molecular_analysis": config.STATUS_COMPLETED,
        "ligand_preparation": config.STATUS_COMPLETED,
        "docking": config.STATUS_FAILED,
        "ranking": config.STATUS_SKIPPED,
        "explainability": config.STATUS_COMPLETED,
        "report": config.STATUS_COMPLETED,
    }
    assert response["message"] == screening.NOTHING_DOCKED_MESSAGE


def test_stage_statuses_for_a_full_run(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(screening, "run_docking", _stub_docking_success)
    response = run_screening(TARGET, config.MIN_COMPOUNDS)

    assert set(response["stages"].values()) == {config.STATUS_COMPLETED}
    assert response["message"] == screening.COMPLETE_MESSAGE
    for stage in response["pipeline"]:
        assert stage["label"] == config.PIPELINE_STAGE_LABELS[stage["key"]]
        assert stage["detail"]


@pytest.mark.parametrize(
    ("succeeded", "attempted", "status"),
    [(0, 0, config.STATUS_SKIPPED), (3, 3, config.STATUS_COMPLETED), (0, 3, config.STATUS_FAILED), (2, 3, config.STATUS_PARTIAL)],
)
def test_ratio_status(succeeded: int, attempted: int, status: str) -> None:
    assert screening.ratio_status(succeeded, attempted) == status


def test_timings_are_recorded_for_every_stage() -> None:
    timings = run_screening(TARGET, config.MIN_COMPOUNDS)["timings"]

    for key in (*config.PIPELINE_STAGE_LABELS, "total"):
        assert timings[key] >= 0.0
    assert timings["total"] >= sum(value for key, value in timings.items() if key != "total") - 1e-6


def test_missing_vina_keeps_descriptors_and_skips_ranking(monkeypatch: pytest.MonkeyPatch) -> None:
    def no_vina(compounds, progress=None):
        raise VinaNotFoundError("Vina missing (test). Install it.")

    monkeypatch.setattr(screening, "run_docking", no_vina)
    response = run_screening(TARGET, config.MIN_COMPOUNDS)

    assert response["stages"]["docking"] == config.STATUS_FAILED
    assert response["stages"]["ligand_preparation"] == config.STATUS_SKIPPED
    assert response["stages"]["ranking"] == config.STATUS_SKIPPED
    assert response["stages"]["molecular_analysis"] == config.STATUS_COMPLETED
    docking_stage = next(stage for stage in response["pipeline"] if stage["key"] == "docking")
    assert docking_stage["detail"] == "Vina missing (test)."
    assert response["summary"]["valid"] > 0
    assert response["summary"]["docked"] == 0
    assert all(compound["docking"] is None for compound in response["compounds"])
    assert "Vina missing" in response["docking"]["error"]
    # Without docking, nothing may receive a priority score.
    assert response["summary"]["ranked"] == 0
    assert all(compound["ranking"]["ranking_status"] == ranking.STATUS_NOT_DOCKED for compound in response["compounds"])


def test_progress_reports_real_stages_in_order() -> None:
    messages: list[str] = []

    run_screening(TARGET, config.MIN_COMPOUNDS, progress=messages.append)

    assert messages[0] == "Loading compounds..."
    stage_positions = [
        next(i for i, message in enumerate(messages) if keyword in message)
        for keyword in ("descriptors", "3D ligands", "Screening Priority Score")
    ]
    assert stage_positions == sorted(stage_positions)
    assert messages[-1] == "Complete."
    assert not any("AI" in message for message in messages)


def test_failed_docking_is_not_ranked_in_pipeline() -> None:
    response = run_screening(TARGET, config.MIN_COMPOUNDS)  # Stub: every docking fails.

    assert response["summary"]["ranked"] == 0
    for compound in response["compounds"]:
        assert compound["ranking"]["rank"] is None
        assert compound["ranking"]["screening_priority_score"] is None


def test_pipeline_ranks_every_docked_compound(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(screening, "run_docking", _stub_docking_success)
    response = run_screening(TARGET, config.MIN_COMPOUNDS)

    entries = response["ranking"]["ranked_compounds"]
    assert response["stages"]["ranking"] == config.STATUS_COMPLETED
    assert response["summary"]["ranked"] == response["summary"]["docked"] == response["summary"]["valid"]
    assert [entry["rank"] for entry in entries] == list(range(1, len(entries) + 1))
    for compound in response["compounds"]:
        assert compound["ranking"]["compound_id"] == compound["compound_id"]


@pytest.mark.integration
def test_full_pipeline_with_real_vina(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    try:
        docking.find_vina_executable()
    except VinaNotFoundError as error:
        pytest.skip(str(error))
    monkeypatch.setattr(config, "DOCKING_RESULTS_DIR", tmp_path)

    first = run_screening(TARGET, config.MIN_COMPOUNDS)
    second = run_screening(TARGET, config.MIN_COMPOUNDS)

    ranked = [entry for entry in first["ranking"]["ranked_compounds"] if entry["rank"] is not None]
    assert len(ranked) == first["summary"]["docked"] > 0
    scores = [entry["screening_priority_score"] for entry in ranked]
    assert scores == sorted(scores, reverse=True)
    assert all(0.0 <= score <= 100.0 for score in scores)
    # Same inputs, same seeds: the same ranking and the same explanations.
    assert [(e["compound_id"], e["screening_priority_score"]) for e in first["ranking"]["ranked_compounds"]] == [
        (e["compound_id"], e["screening_priority_score"]) for e in second["ranking"]["ranked_compounds"]
    ]
    assert first["explanations"] == second["explanations"]
    # The real report is generated and agrees with the pipeline.
    assert first["report"]["status"] == config.REPORT_GENERATED
    report_json = json.loads(Path(first["report"]["files"]["json"]).read_text(encoding="utf-8"))
    assert [item["screening_priority_score"] for item in report_json["ranking"]["ranked"]] == scores
    for entry, explanation in zip(first["ranking"]["ranked_compounds"], first["explanations"]):
        if entry["rank"] is not None:
            assert explanation["docking"]["points"] + explanation["properties"]["points"] == entry["screening_priority_score"]


@pytest.mark.parametrize(
    ("target", "n_compounds"),
    [("Unknown Protein", 10), (TARGET, config.MIN_COMPOUNDS - 1), (TARGET, config.MAX_COMPOUNDS + 1)],
)
def test_invalid_requests_raise(target: str, n_compounds: int) -> None:
    with pytest.raises(ValueError):
        run_screening(target, n_compounds)


def test_message_never_claims_full_success_when_a_stage_is_partial() -> None:
    pipeline = [
        {"label": "Molecular analysis", "status": config.STATUS_PARTIAL},
        {"label": "Report generation", "status": config.STATUS_FAILED},
        {"label": "Docking", "status": config.STATUS_COMPLETED},
    ]

    message = screening.outcome_message(pipeline, docking_error=None, docked=5)

    assert message != screening.COMPLETE_MESSAGE
    assert "Molecular analysis (partial)" in message and "Report generation (failed)" in message
    assert "Docking" not in message


def test_message_after_report_failure(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    blocker = tmp_path / "blocked"
    blocker.write_text("x", encoding="utf-8")
    monkeypatch.setattr(config, "REPORTS_DIR", blocker)
    monkeypatch.setattr(screening, "run_docking", _stub_docking_success)

    response = run_screening(TARGET, config.MIN_COMPOUNDS)

    assert "Report generation (failed)" in response["message"]
