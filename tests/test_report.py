"""Tests for report generation.

Screening responses come from the real pipeline (``run_screening``), with
docking replaced by stubs so the tests stay fast. Stub affinities are test
doubles only and never reach the app.
"""

import copy
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest

import config
from backend import compounds as compounds_module
from backend import explainability, ranking, report, screening
from backend.compounds import CompoundLoadResult
from backend.docking import VinaNotFoundError

TARGET = config.AVAILABLE_TARGETS[0]
FIXED_NOW = datetime(2026, 10, 24, 12, 30, 0, tzinfo=timezone.utc)

SECTIONS = (
    "# Fairouz Virtual Screening Report",
    "## Run Summary",
    "## Target",
    "## Methodology",
    "## Screening Results",
    "## Docking Results",
    "## Molecular Descriptors",
    "## Unranked Compounds",
    "## Explainability Summary",
    "## Limitations",
    "## Responsible Use",
    "## Reproducibility",
)

FORBIDDEN_PHRASES = (
    "proven treatment", "effective drug", "safe drug", "cure", "clinically effective",
    "successful therapy", "guaranteed binding", "confirmed antiviral", "antiviral activity was",
    "effective antiviral", "likely treatment", "promising drug", "safe compound",
    "clinically useful", "likely to work", "therapeutic probability",
)


# --- Building realistic responses -------------------------------------------------

def _docking_result(compound, index, *, success=True, ligand_prepared=True):
    return {
        "compound_id": compound["compound_id"],
        "name": compound["name"],
        "ligand_prepared": ligand_prepared,
        "docking_success": success,
        "best_affinity_kcal_mol": (-9.0 + 0.5 * index) if success else None,
        "poses_generated": 9 if success else 0,
        "pose_file": str(config.PROJECT_ROOT / "results" / "docking" / "stub" / f"{compound['compound_id']}.pdbqt") if success else None,
        "error": None if success else ("Ligand preparation failed: stub" if not ligand_prepared else "Docking failed: stub"),
    }


def _stub_docking(fail_ids=(), ligand_fail_ids=()):
    def run(compounds, progress=None):
        return {
            "stage": "docking",
            "status": config.STATUS_COMPLETED,
            "results": [
                _docking_result(
                    compound,
                    index,
                    success=compound["compound_id"] not in fail_ids and compound["compound_id"] not in ligand_fail_ids,
                    ligand_prepared=compound["compound_id"] not in ligand_fail_ids,
                )
                for index, compound in enumerate(compounds)
            ],
            "run_dir": str(config.PROJECT_ROOT / "results" / "docking" / "stub"),
            "metadata": {
                "vina": {"version": "stub Vina", "exhaustiveness": 8, "num_modes": 9, "seed": 42},
                "ligand_preparation": {"conformer": "RDKit ETKDGv3", "force_field": "MMFF94 (UFF fallback)", "seed": 42},
            },
        }

    return run


@pytest.fixture
def ranked_response(monkeypatch: pytest.MonkeyPatch) -> dict:
    monkeypatch.setattr(screening, "run_docking", _stub_docking())
    return screening.run_screening(TARGET, config.MIN_COMPOUNDS)


@pytest.fixture
def mixed_response(monkeypatch: pytest.MonkeyPatch) -> dict:
    """Real library plus an invalid molecule, one ligand-prep failure and one docking failure."""
    library = compounds_module.load_compounds()
    extra = {"compound_id": "FZ-BAD", "name": "Broken SMILES", "smiles": "C1CC(", "source": "test"}
    monkeypatch.setattr(screening, "load_compounds", lambda: CompoundLoadResult([*library.compounds, extra], []))
    monkeypatch.setattr(screening, "run_docking", _stub_docking(fail_ids=("FZ-0002",), ligand_fail_ids=("FZ-0003",)))
    return screening.run_screening(TARGET, 11)


def generate(response, tmp_path, **kwargs):
    return report.generate_report(response, output_dir=tmp_path / "reports", now=FIXED_NOW, **kwargs)


def read_json(result):
    return json.loads(Path(result["files"]["json"]).read_text(encoding="utf-8"))


def read_md(result):
    return Path(result["files"]["md"]).read_text(encoding="utf-8")


# --- Files and folders -------------------------------------------------------------

def test_creates_markdown_json_and_html(ranked_response: dict, tmp_path: Path) -> None:
    result = generate(ranked_response, tmp_path)

    assert result["status"] == config.REPORT_GENERATED
    assert set(result["files"]) == {"md", "json", "html"}
    assert all(Path(path).is_file() for path in result["files"].values())
    assert result["errors"] == {}


def test_creates_missing_report_directory(ranked_response: dict, tmp_path: Path) -> None:
    target = tmp_path / "deep" / "nested" / "reports"

    result = report.generate_report(ranked_response, output_dir=target, now=FIXED_NOW)

    assert target.is_dir() and result["status"] == config.REPORT_GENERATED


def test_never_overwrites_existing_reports(ranked_response: dict, tmp_path: Path) -> None:
    first = generate(ranked_response, tmp_path)
    second = generate(ranked_response, tmp_path)

    assert first["files"]["md"] != second["files"]["md"]
    assert Path(second["files"]["md"]).stem.endswith("_2")
    assert Path(first["files"]["md"]).is_file()


def test_filenames_are_windows_safe(ranked_response: dict, tmp_path: Path) -> None:
    result = generate(ranked_response, tmp_path)

    for path in result["files"].values():
        assert re.fullmatch(r"[A-Za-z0-9._-]+", Path(path).name)
    assert Path(result["files"]["md"]).name == "fairouz_screening_6LU7_2026-10-24_123000.md"


@pytest.mark.parametrize(
    ("raw", "expected"),
    [("6LU7", "6LU7"), ("a:b<c>d|e?f*g", "a_b_c_d_e_f_g"), ("..\\..\\x", "x"), ("", "report"), ("CON / name", "CON_name")],
)
def test_safe_filename_part(raw: str, expected: str) -> None:
    assert report.safe_filename_part(raw) == expected


def test_json_has_expected_top_level_fields(ranked_response: dict, tmp_path: Path) -> None:
    data = read_json(generate(ranked_response, tmp_path))

    for field in ("metadata", "target", "methodology", "summary", "compounds", "ranking", "explanations",
                  "limitations", "responsible_use", "reproducibility"):
        assert field in data


def test_markdown_has_sections_in_order(ranked_response: dict, tmp_path: Path) -> None:
    text = read_md(generate(ranked_response, tmp_path))

    positions = [text.index(section + "\n") for section in SECTIONS]
    assert positions == sorted(positions)
    assert "### Top Ranked Compounds" in text and "### Full Ranking" in text


def test_html_is_standalone_and_escaped(ranked_response: dict, tmp_path: Path) -> None:
    response = copy.deepcopy(ranked_response)
    response["compounds"][0]["name"] = "<script>alert(1)</script>"
    html_text = Path(generate(response, tmp_path)["files"]["html"]).read_text(encoding="utf-8")

    assert html_text.startswith("<!doctype html>")
    assert "<script>alert(1)</script>" not in html_text
    assert "&lt;script&gt;" in html_text


# --- Metadata and target -------------------------------------------------------------

def test_run_metadata(ranked_response: dict, tmp_path: Path) -> None:
    result = generate(ranked_response, tmp_path)
    data = read_json(result)
    summary = ranked_response["summary"]

    assert re.fullmatch(r"20261024-123000_6LU7_[0-9a-f]{10}-[0-9a-f]{6}", data["metadata"]["run_id"])
    assert result["run_id"] == data["metadata"]["run_id"] == data["reproducibility"]["run_id"]
    assert data["metadata"]["generated_at"] == FIXED_NOW.isoformat()
    assert data["metadata"]["project_version"] == config.PROJECT_VERSION
    assert data["summary"]["pdb_id"] == "6LU7"
    assert data["summary"]["requested"] == summary["requested"]
    assert data["summary"]["processed"] == summary["processed"]
    assert data["summary"]["successfully_analysed"] == summary["valid"]
    assert data["summary"]["docked"] == summary["docked"]
    assert data["summary"]["ranked"] == summary["ranked"]
    assert data["summary"]["unranked"] == summary["processed"] - summary["ranked"]


def test_target_data_comes_from_project_sources(ranked_response: dict, tmp_path: Path) -> None:
    target = read_json(generate(ranked_response, tmp_path))["target"]

    assert target["target_id"] == config.TARGET_ID
    assert target["structure"]["method"] == "X-RAY DIFFRACTION"
    assert target["structure"]["resolution_angstrom"] == 2.16
    assert target["structure"]["citation_doi"] == "10.1038/S41586-020-2223-Y"
    assert target["binding_box"]["center"] == list(config.BINDING_BOX_CENTER)
    assert target["binding_box"]["size"] == list(config.BINDING_BOX_SIZE)
    assert len(target["preparation"]["steps"]) == 7


def test_reproducibility_metadata(ranked_response: dict, tmp_path: Path) -> None:
    repro = read_json(generate(ranked_response, tmp_path))["reproducibility"]

    assert repro["configuration"]["ranking_weights"] == {"docking": 0.7, "properties": 0.3}
    assert repro["configuration"]["random_seed"] == config.RANDOM_SEED
    assert repro["integrity"]["raw_target"]["sha256"] == config.RAW_TARGET_SHA256
    assert repro["integrity"]["prepared_receptor"]["sha256"]
    assert "line endings" in repro["integrity"]["prepared_receptor"]["note"]
    assert repro["software_versions"]["python"]
    assert repro["software_versions"]["autodock_vina"] == "stub Vina"


def test_configuration_hash_is_deterministic() -> None:
    first = report.configuration_hash(report.configuration_snapshot())
    second = report.configuration_hash(report.configuration_snapshot())

    assert first == second and re.fullmatch(r"[0-9a-f]{10}", first)


# --- Values come straight from the pipeline -------------------------------------

def test_ranked_compounds_and_scores_match_pipeline(ranked_response: dict, tmp_path: Path) -> None:
    data = read_json(generate(ranked_response, tmp_path))
    entries = [entry for entry in ranked_response["ranking"]["ranked_compounds"] if entry["rank"] is not None]

    assert [item["compound_id"] for item in data["ranking"]["ranked"]] == [entry["compound_id"] for entry in entries]
    for item, entry in zip(data["ranking"]["ranked"], entries):
        assert item["rank"] == entry["rank"]
        assert item["screening_priority_score"] == entry["screening_priority_score"]
        assert item["docking_points"] == entry["docking_points"]
        assert item["property_points"] == entry["property_points"]


def test_top_three(ranked_response: dict, tmp_path: Path) -> None:
    data = read_json(generate(ranked_response, tmp_path))

    assert data["ranking"]["top"] == data["ranking"]["ranked"][:3]
    assert [item["rank"] for item in data["ranking"]["top"]] == [1, 2, 3]
    text = read_md(generate(ranked_response, tmp_path))
    for item in data["ranking"]["top"]:
        assert f"#### #{item['rank']} · {item['name']} ({item['compound_id']})" in text


def test_docking_values_match_pipeline(ranked_response: dict, tmp_path: Path) -> None:
    data = read_json(generate(ranked_response, tmp_path))

    for record, compound in zip(data["compounds"], ranked_response["compounds"]):
        assert record["docking"]["predicted_affinity_kcal_mol"] == compound["docking"]["best_affinity_kcal_mol"]
        assert record["docking"]["status"] == "Docked"
        assert not Path(record["docking"]["pose_file"]).is_absolute()


def test_descriptor_values_match_pipeline(ranked_response: dict, tmp_path: Path) -> None:
    data = read_json(generate(ranked_response, tmp_path))

    for record, compound in zip(data["compounds"], ranked_response["compounds"]):
        for key in config.DESCRIPTORS:
            assert record["descriptors"][key] == compound[key]


def test_explainability_values_match_pipeline(ranked_response: dict, tmp_path: Path) -> None:
    data = read_json(generate(ranked_response, tmp_path))

    for record, explanation in zip(data["explanations"], ranked_response["explanations"]):
        assert record["docking_points"] == explanation["docking"]["points"]
        assert record["property_points"] == explanation["properties"]["points"]
        assert record["main_positive_factor"] == explanation["main_positive_factor"]
        assert record["main_limiting_factor"] == explanation["main_limiting_factor"]
        assert record["primary_driver"] == explanation["primary_driver"]


def test_report_reads_ranking_instead_of_recomputing(ranked_response: dict, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    def forbidden(*args, **kwargs):
        raise AssertionError("report must not recompute rankings or explanations")

    for name in ("rank_compounds", "score_compound", "docking_component", "property_subscores"):
        monkeypatch.setattr(ranking, name, forbidden)
    monkeypatch.setattr(explainability, "explain_rankings", forbidden)
    response = copy.deepcopy(ranked_response)
    first = next(entry for entry in response["ranking"]["ranked_compounds"] if entry["rank"] == 1)
    first["screening_priority_score"] = 12.345  # A sentinel no formula would produce.

    data = read_json(generate(response, tmp_path))

    assert data["ranking"]["ranked"][0]["screening_priority_score"] == 12.345


def test_pipeline_input_is_not_mutated(ranked_response: dict, tmp_path: Path) -> None:
    before = copy.deepcopy(ranked_response)

    generate(ranked_response, tmp_path)

    assert ranked_response == before


# --- Partial and failed runs -----------------------------------------------------

def test_unranked_compounds_appear_with_stage_and_reason(mixed_response: dict, tmp_path: Path) -> None:
    result = generate(mixed_response, tmp_path)
    data, text = read_json(result), read_md(result)
    unranked = {item["compound_id"]: item for item in data["ranking"]["unranked"]}

    assert unranked["FZ-BAD"]["status"] == ranking.STATUS_INVALID
    assert unranked["FZ-BAD"]["failed_stage"] == "SMILES validation (RDKit)"
    assert unranked["FZ-0002"]["status"] == ranking.STATUS_DOCKING_FAILED
    assert unranked["FZ-0003"]["status"] == ranking.STATUS_LIGAND_FAILED
    assert "none was fabricated" in text
    assert "Broken SMILES" in text
    assert data["summary"]["completeness"] == "partial"
    assert data["summary"]["unranked"] == 3


def test_vina_missing_report(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    def no_vina(compounds, progress=None):
        raise VinaNotFoundError("AutoDock Vina was not found (test).")

    monkeypatch.setattr(screening, "run_docking", no_vina)
    response = screening.run_screening(TARGET, config.MIN_COMPOUNDS)
    result = generate(response, tmp_path)
    data, text = read_json(result), read_md(result)

    assert result["status"] == config.REPORT_GENERATED
    assert data["ranking"]["top"] == [] and data["ranking"]["ranked"] == []
    assert "### Top Ranked Compounds" not in text
    assert "Docking could not run" in text
    assert "did not run in this session" in text
    assert data["summary"]["completeness"] == "partial"


def test_zero_ranked_report(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    library_ids = [compound["compound_id"] for compound in compounds_module.load_compounds().compounds]
    monkeypatch.setattr(screening, "run_docking", _stub_docking(fail_ids=tuple(library_ids)))
    response = screening.run_screening(TARGET, config.MIN_COMPOUNDS)
    result = generate(response, tmp_path)
    text = read_md(result)

    assert result["status"] == config.REPORT_GENERATED
    assert "### Top Ranked Compounds" not in text
    assert "No compound received a Screening Priority Score" in text
    assert read_json(result)["summary"]["ranked"] == 0


def test_empty_compound_list(tmp_path: Path) -> None:
    empty = {
        "summary": {"requested": 10, "loaded": 0, "processed": 0, "valid": 0, "invalid": 0,
                    "ligands_prepared": 0, "docked": 0, "docking_failed": 0, "ranked": 0},
        "stages": {}, "compounds": [], "docking": {"run_dir": None, "metadata": None, "error": None},
        "ranking": {"ranked_compounds": [], "method": None}, "explanations": [], "warnings": [],
    }
    result = generate(empty, tmp_path)
    data = read_json(result)

    assert result["status"] == config.REPORT_GENERATED
    assert data["compounds"] == [] and data["ranking"]["ranked"] == []
    assert data["summary"]["completeness"] == "empty"
    assert "No compounds were processed" in read_md(result)


def test_completely_empty_response_does_not_crash(tmp_path: Path) -> None:
    result = generate({}, tmp_path)

    assert result["status"] == config.REPORT_GENERATED


def test_single_compound(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    only = compounds_module.load_compounds().compounds[:1]
    monkeypatch.setattr(screening, "load_compounds", lambda: CompoundLoadResult(only, []))
    monkeypatch.setattr(screening, "run_docking", _stub_docking())
    data = read_json(generate(screening.run_screening(TARGET, config.MIN_COMPOUNDS), tmp_path))

    assert len(data["ranking"]["ranked"]) == 1
    assert len(data["ranking"]["top"]) == 1
    assert data["summary"]["completeness"] == "complete"


# --- Failure handling ------------------------------------------------------------------

def test_unwritable_folder_gives_not_generated(ranked_response: dict, tmp_path: Path) -> None:
    blocker = tmp_path / "not_a_folder"
    blocker.write_text("x", encoding="utf-8")

    result = report.generate_report(ranked_response, output_dir=blocker, now=FIXED_NOW)

    assert result["status"] == config.REPORT_NOT_GENERATED
    assert result["files"] == {}
    assert "folder" in result["errors"]


def test_one_failing_format_gives_partial(ranked_response: dict, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    real_render = report._render

    def render_or_fail(data, blocks, fmt):
        if fmt == "html":
            raise OSError("disk full (test)")
        return real_render(data, blocks, fmt)

    monkeypatch.setattr(report, "_render", render_or_fail)
    result = generate(ranked_response, tmp_path)

    assert result["status"] == config.REPORT_PARTIAL
    assert set(result["files"]) == {"md", "json"}
    assert "disk full" in result["errors"]["html"]


def test_build_failure_is_contained(ranked_response: dict, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    def broken(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(report, "build_report_data", broken)
    result = generate(ranked_response, tmp_path)

    assert result["status"] == config.REPORT_NOT_GENERATED
    assert "boom" in result["errors"]["build"]


def test_report_failure_does_not_lose_screening_results(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    blocker = tmp_path / "blocked"
    blocker.write_text("x", encoding="utf-8")
    monkeypatch.setattr(config, "REPORTS_DIR", blocker)
    monkeypatch.setattr(screening, "run_docking", _stub_docking())

    response = screening.run_screening(TARGET, config.MIN_COMPOUNDS)

    assert response["stages"]["report"] == config.STATUS_FAILED
    assert response["report"]["status"] == config.REPORT_NOT_GENERATED
    assert response["summary"]["ranked"] > 0
    assert len(response["compounds"]) == response["summary"]["processed"]
    assert any("Report not generated" in warning for warning in response["warnings"])


def test_missing_optional_metadata(ranked_response: dict, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr(report, "git_revision", lambda: None)
    monkeypatch.setattr(report, "_package_version", lambda name: None)
    monkeypatch.setattr(config, "PREPARATION_LOG_FILE", tmp_path / "missing_log.json")
    response = copy.deepcopy(ranked_response)
    response["docking"]["metadata"] = None

    result = generate(response, tmp_path)
    data, text = read_json(result), read_md(result)

    assert result["status"] == config.REPORT_GENERATED
    assert data["metadata"]["git"] is None
    assert data["reproducibility"]["software_versions"]["rdkit"] is None
    assert "not available" in text
    assert "not detected" in text


# --- Serialisation -------------------------------------------------------------------

def test_numpy_values_serialise() -> None:
    converted = report.to_jsonable(
        {"f": np.float64(1.5), "i": np.int64(3), "b": np.bool_(True), "a": np.array([1.0, 2.0]), "n": np.float32("nan")}
    )

    assert converted == {"f": 1.5, "i": 3, "b": True, "a": [1.0, 2.0], "n": None}
    json.dumps(converted, allow_nan=False)


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_nan_and_infinity_become_null(bad: float) -> None:
    assert report.to_jsonable({"x": bad}) == {"x": None}


def test_paths_datetimes_and_tuples_serialise() -> None:
    converted = report.to_jsonable({"p": Path("a") / "b.txt", "d": FIXED_NOW, "t": (1, 2), "s": {"b", "a"}})

    assert converted == {"p": "a/b.txt", "d": FIXED_NOW.isoformat(), "t": [1, 2], "s": ["a", "b"]}


def test_report_with_numpy_and_nan_inputs_writes_valid_json(ranked_response: dict, tmp_path: Path) -> None:
    response = copy.deepcopy(ranked_response)
    response["compounds"][0]["molecular_weight"] = np.float64(response["compounds"][0]["molecular_weight"])
    response["compounds"][1]["logp"] = math.nan

    result = generate(response, tmp_path)
    data = read_json(result)

    assert result["status"] == config.REPORT_GENERATED
    assert data["compounds"][1]["descriptors"]["logp"] is None


# --- Wording and determinism ----------------------------------------------------------

@pytest.mark.parametrize("fixture_name", ["ranked_response", "mixed_response"])
def test_report_wording_avoids_unsupported_claims(fixture_name: str, request: pytest.FixtureRequest, tmp_path: Path) -> None:
    result = generate(request.getfixturevalue(fixture_name), tmp_path)

    for path in result["files"].values():
        text = Path(path).read_text(encoding="utf-8").lower()
        for phrase in FORBIDDEN_PHRASES:
            assert not re.search(rf"\b{re.escape(phrase)}", text), phrase
        assert "experimental validation is required" in text


def test_report_uses_cautious_wording(ranked_response: dict, tmp_path: Path) -> None:
    text = read_md(generate(ranked_response, tmp_path))

    for phrase in ("Predicted AutoDock Vina affinity", "computational estimate", "prioritised for further investigation",
                   "not an experimental binding energy", "research prioritisation only"):
        assert phrase in text


def test_methodology_wording_is_deterministic(ranked_response: dict) -> None:
    first = report.build_report_data(ranked_response, FIXED_NOW)
    second = report.build_report_data(copy.deepcopy(ranked_response), FIXED_NOW)

    assert first["methodology"] == second["methodology"]
    assert report.render_markdown(report.build_document(first)) == report.render_markdown(report.build_document(second))
    assert "Screening Priority Score = docking_points + property_points" in json.dumps(first["methodology"])


def test_one_sided_ranges_are_not_described_as_zero_at_zero(ranked_response: dict, tmp_path: Path) -> None:
    text = read_md(generate(ranked_response, tmp_path))

    assert "TPSA full score up to 140" in text
    assert "zero at ≤0 " not in text


# --- Pipeline integration ----------------------------------------------------------

def test_pipeline_generates_report(ranked_response: dict) -> None:
    result = ranked_response["report"]

    assert ranked_response["stages"]["report"] == config.STATUS_COMPLETED
    assert result["status"] == config.REPORT_GENERATED
    assert all(Path(path).is_file() for path in result["files"].values())
    assert Path(result["files"]["md"]).parent == config.REPORTS_DIR  # Redirected to tmp by conftest.
    assert result["preview"]["ranked"] == ranked_response["summary"]["ranked"]
    assert len(result["preview"]["top"]) == 3
    assert read_json(result)["metadata"]["run_id"] == result["run_id"]


def test_ui_report_panel_renders(ranked_response: dict, mixed_response: dict) -> None:
    from streamlit.testing.v1 import AppTest

    def panel_app(report_result):
        from UI.components.report_panel import render_report

        render_report(report_result)

    for response in (ranked_response, mixed_response):
        app = AppTest.from_function(panel_app, args=(response["report"],)).run()
        assert not app.exception
        assert any("Report status" in element.value for element in app.success)


def test_different_outcomes_in_the_same_second_get_different_run_ids(ranked_response: dict, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    def no_vina(compounds, progress=None):
        raise VinaNotFoundError("AutoDock Vina was not found (test).")

    monkeypatch.setattr(screening, "run_docking", no_vina)
    unranked_response = screening.run_screening(TARGET, config.MIN_COMPOUNDS)

    first = generate(ranked_response, tmp_path)["run_id"]
    second = generate(unranked_response, tmp_path)["run_id"]
    again = generate(ranked_response, tmp_path)["run_id"]

    assert first != second
    assert first == again  # Same configuration, same results, same second: same ID.


# --- Pipeline status model in reports and degraded runs -------------------------

def test_mixed_run_reports_partial_stages(mixed_response: dict) -> None:
    stages = mixed_response["stages"]

    assert stages["molecular_analysis"] == config.STATUS_PARTIAL  # One invalid SMILES.
    assert stages["ligand_preparation"] == config.STATUS_PARTIAL  # One ligand failed.
    assert stages["docking"] == config.STATUS_PARTIAL  # One docking failed.
    assert stages["ranking"] == config.STATUS_COMPLETED  # Every docked compound was ranked.
    assert stages["report"] == config.STATUS_COMPLETED


def test_report_includes_pipeline_and_labelled_timings(mixed_response: dict, tmp_path: Path) -> None:
    result = generate(mixed_response, tmp_path)
    data, text = read_json(result), read_md(result)

    # The report lists exactly the stages in the response it was built from.
    assert [(stage["stage"], stage["status"]) for stage in data["summary"]["pipeline"]] == [
        (stage["label"], stage["status"]) for stage in mixed_response["pipeline"]
    ]
    assert set(data["metadata"]["timings_seconds"]) >= {"compound_loading", "docking", "ranking"}
    assert "### Pipeline stages" in text
    assert "they vary between runs and are not a benchmark" in text.lower()


def test_run_id_is_windows_safe_and_shows_results_hash(ranked_response: dict, tmp_path: Path) -> None:
    run_id = generate(ranked_response, tmp_path)["run_id"]

    assert re.fullmatch(r"[A-Za-z0-9_-]+", run_id)
    assert run_id.endswith(report.results_hash(ranked_response))


def test_results_hash_ignores_timings(ranked_response: dict) -> None:
    changed = copy.deepcopy(ranked_response)
    for compound in changed["compounds"]:
        compound["docking"]["docking_seconds"] = 999.0
    changed["timings"]["total"] = 999.0

    assert report.results_hash(changed) == report.results_hash(ranked_response)


def test_report_affinities_use_display_precision(ranked_response: dict, tmp_path: Path) -> None:
    text = read_md(generate(ranked_response, tmp_path))
    best = ranked_response["ranking"]["ranked_compounds"][0]

    assert f"Predicted AutoDock Vina affinity: {best['docking_affinity']:.1f} kcal/mol" in text
    assert "-0.00" not in text
