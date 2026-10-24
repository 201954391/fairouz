"""Screening pipeline coordinator.

Runs the virtual screening workflow; each stage consumes the previous
stage's actual output and nothing is recalculated downstream:

1. Compound loading (``backend.compounds``)
2. Molecular analysis: RDKit descriptors (``backend.descriptors``)
3. Ligand preparation (``backend.ligands``, run inside docking)
4. Docking with AutoDock Vina (``backend.docking``)
5. Ranking by Screening Priority Score (``backend.ranking``)
6. Explainability (``backend.explainability``)
7. Report generation (``backend.report``)

Each stage reports a status (completed, partial, failed or skipped), a short
detail and its wall-clock time. If docking cannot run (for example Vina is
not installed), descriptors are still returned and ranking is skipped. The
report is generated last; if writing it fails, the results are still
returned. No values are fabricated.
"""

import logging
import time
from collections.abc import Callable
from typing import Any

import config
from backend.compounds import load_compounds
from backend.descriptors import calculate_descriptors
from backend.docking import DockingSetupError, run_docking
from backend.explainability import explain_rankings
from backend.ranking import rank_compounds
from backend.report import generate_report

logger = logging.getLogger(__name__)

COMPLETE_MESSAGE = "Screening complete: molecular analysis, docking, ranking, explanations and report."
DOCKING_UNAVAILABLE_MESSAGE = (
    "Molecular analysis complete, but docking could not run, so ranking was skipped. "
    "The report describes what was completed."
)
NOTHING_DOCKED_MESSAGE = (
    "Molecular analysis complete, but no compound docked successfully, so ranking was skipped."
)
REPORT_STAGE_STATUS = {
    config.REPORT_GENERATED: config.STATUS_COMPLETED,
    config.REPORT_PARTIAL: config.STATUS_PARTIAL,
    config.REPORT_NOT_GENERATED: config.STATUS_FAILED,
}


def validate_screening_request(target: str, n_compounds: int) -> None:
    """Check that a screening request is within the supported options.

    Args:
        target: Name of the viral protein target.
        n_compounds: Number of compounds to screen.

    Raises:
        ValueError: If the target is unknown or the count is out of range.
    """
    if target not in config.AVAILABLE_TARGETS:
        available = ", ".join(config.AVAILABLE_TARGETS)
        raise ValueError(f"Unknown target '{target}'. Available targets: {available}.")

    if not config.MIN_COMPOUNDS <= n_compounds <= config.MAX_COMPOUNDS:
        raise ValueError(
            f"Number of compounds must be between {config.MIN_COMPOUNDS} "
            f"and {config.MAX_COMPOUNDS}, got {n_compounds}."
        )


def _shortfall_warning(requested: int, available: int) -> str | None:
    """Explain when fewer compounds exist than were requested."""
    if available >= requested:
        return None
    return (
        f"{requested} compounds were requested, but the dataset contains only "
        f"{available}. All {available} available compounds were processed."
    )


def _attach_docking(compounds: list[dict[str, Any]], docking_results: list[dict[str, Any]]) -> None:
    """Store each compound's docking result under its ``docking`` key.

    Compounds that were not docked (invalid molecules, or docking
    unavailable) get ``None``.
    """
    by_id = {result["compound_id"]: result for result in docking_results}
    for compound in compounds:
        compound["docking"] = by_id.get(compound["compound_id"])


def _attach_ranking(compounds: list[dict[str, Any]], ranking_entries: list[dict[str, Any]]) -> None:
    """Store each compound's ranking entry (ranked or not) under its ``ranking`` key."""
    by_id = {entry["compound_id"]: entry for entry in ranking_entries}
    for compound in compounds:
        compound["ranking"] = by_id.get(compound["compound_id"])


# --- Stage status model ------------------------------------------------------------

def ratio_status(succeeded: int, attempted: int) -> str:
    """Stage status from how many items succeeded out of those attempted."""
    if attempted == 0:
        return config.STATUS_SKIPPED
    if succeeded == attempted:
        return config.STATUS_COMPLETED
    if succeeded == 0:
        return config.STATUS_FAILED
    return config.STATUS_PARTIAL


def _stage(key: str, status: str, detail: str, seconds: float | None) -> dict[str, Any]:
    return {
        "key": key,
        "label": config.PIPELINE_STAGE_LABELS[key],
        "status": status,
        "detail": detail,
        "seconds": seconds,
    }


def _first_sentence(text: str) -> str:
    return text.split(". ")[0].rstrip(".") + "."


def build_pipeline_stages(
    summary: dict[str, Any],
    skipped_rows: int,
    docking_error: str | None,
    timings: dict[str, float],
) -> list[dict[str, Any]]:
    """Describe stages 1–6 from the run's real counts (report is added afterwards)."""
    loaded, processed, valid = summary["loaded"], summary["processed"], summary["valid"]
    prepared, docked, ranked = summary["ligands_prepared"], summary["docked"], summary["ranked"]

    loading_detail = f"{loaded} usable compounds in the dataset; {processed} selected."
    if skipped_rows:
        loading_detail += f" {skipped_rows} malformed row(s) skipped."

    if docking_error:
        ligand = _stage("ligand_preparation", config.STATUS_SKIPPED, "Skipped because docking could not run.", None)
        docking = _stage("docking", config.STATUS_FAILED, _first_sentence(docking_error), timings.get("docking"))
    else:
        ligand = _stage(
            "ligand_preparation", ratio_status(prepared, valid),
            f"{prepared} of {valid} valid molecules converted to 3D ligands.", timings.get("ligand_preparation"),
        )
        docking = _stage(
            "docking", ratio_status(docked, prepared),
            f"{docked} of {prepared} prepared ligands docked with AutoDock Vina.", timings.get("docking"),
        )

    if docked == 0:
        ranking = _stage("ranking", config.STATUS_SKIPPED, "Skipped: no docking results to rank.", timings.get("ranking"))
    else:
        ranking = _stage(
            "ranking", ratio_status(ranked, docked),
            f"{ranked} of {docked} docked compounds given a Screening Priority Score.", timings.get("ranking"),
        )

    return [
        _stage("compound_loading", config.STATUS_PARTIAL if skipped_rows else config.STATUS_COMPLETED,
               loading_detail, timings.get("compound_loading")),
        _stage("molecular_analysis", ratio_status(valid, processed),
               f"RDKit descriptors calculated for {valid} of {processed} compounds.", timings.get("molecular_analysis")),
        ligand,
        docking,
        ranking,
        _stage(
            "explainability",
            config.STATUS_COMPLETED if processed else config.STATUS_SKIPPED,
            f"{processed} explanations: {ranked} ranked, {processed - ranked} explaining why a compound was not ranked.",
            timings.get("explainability"),
        ),
    ]


def outcome_message(pipeline: list[dict[str, Any]], docking_error: str | None, docked: int) -> str:
    """One-line outcome built from the final stage statuses (so it is never over-optimistic)."""
    if docking_error:
        return DOCKING_UNAVAILABLE_MESSAGE
    if docked == 0:
        return NOTHING_DOCKED_MESSAGE
    issues = [f"{stage['label']} ({stage['status']})" for stage in pipeline if stage["status"] != config.STATUS_COMPLETED]
    if not issues:
        return COMPLETE_MESSAGE
    return f"Screening finished; some stages did not fully succeed: {', '.join(issues)}. Details are listed per stage."


def _report_stage(report: dict[str, Any], seconds: float) -> dict[str, Any]:
    if report["status"] == config.REPORT_GENERATED:
        detail = f"Markdown, JSON and HTML written (run {report['run_id']})."
    else:
        detail = "; ".join(report["errors"].values()) or "No report files were written."
    return _stage("report", REPORT_STAGE_STATUS[report["status"]], detail, seconds)


# --- The pipeline ------------------------------------------------------------------

def run_screening(
    target: str,
    n_compounds: int,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Run the full screening pipeline.

    Args:
        target: Name of the viral protein target to screen against.
        n_compounds: Number of compounds to screen. If the dataset holds
            fewer, all of them are used; compounds are never duplicated.
        progress: Optional callback receiving a short message as each stage
            (and each docked compound) starts.

    Returns:
        A response with the keys:

        * ``message``: overall outcome for the user.
        * ``request``: the target and requested count.
        * ``summary``: requested, loaded, processed, valid, invalid,
          ligands_prepared, docked, docking_failed and ranked counts.
        * ``pipeline``: one entry per stage (key, label, status, detail, seconds).
        * ``stages``: stage key → status (a compact view of ``pipeline``).
        * ``timings``: seconds per stage and in total (time-dependent).
        * ``compounds``: descriptor results, each with ``docking`` and
          ``ranking`` entries.
        * ``docking``: run folder and settings, or the setup error.
        * ``ranking``: entries in rank order (unranked last) and the method used.
        * ``explanations``: one explanation per ranking entry, in the same order.
        * ``report``: report status, run ID, file paths, errors and preview.
        * ``warnings``: skipped rows, count shortfalls, docking and report problems.

    Raises:
        ValueError: If the request fails validation.
        CompoundDataError: If the compound dataset cannot be used.
    """
    validate_screening_request(target, n_compounds)
    report_progress = progress or (lambda message: None)
    timings: dict[str, float] = {}
    run_started = time.perf_counter()
    logger.info("Screening started: target=%r, n_compounds=%d", target, n_compounds)

    report_progress("Loading compounds...")
    started = time.perf_counter()
    library = load_compounds()
    selected = library.compounds[:n_compounds]
    warnings = list(library.issues)
    shortfall = _shortfall_warning(n_compounds, len(library.compounds))
    if shortfall:
        logger.warning(shortfall)
        warnings.append(shortfall)
    timings["compound_loading"] = time.perf_counter() - started

    report_progress(f"Calculating molecular descriptors for {len(selected)} compounds...")
    started = time.perf_counter()
    analysed = calculate_descriptors(selected)
    valid = [compound for compound in analysed if compound["valid"]]
    timings["molecular_analysis"] = time.perf_counter() - started

    report_progress(f"Preparing 3D ligands and running AutoDock Vina for {len(valid)} compounds (a few seconds each)...")
    started = time.perf_counter()
    try:
        docking = run_docking(valid, progress=report_progress)
        docking_results = docking["results"]
        docking_info = {"run_dir": docking["run_dir"], "metadata": docking["metadata"], "error": None}
    except DockingSetupError as error:
        logger.error("Docking unavailable: %s", error)
        warnings.append(f"Docking could not run: {error}")
        docking_results = []
        docking_info = {"run_dir": None, "metadata": None, "error": str(error)}
    docking_wall = time.perf_counter() - started
    timings["ligand_preparation"] = sum(result.get("ligand_seconds", 0.0) for result in docking_results)
    timings["docking"] = max(0.0, docking_wall - timings["ligand_preparation"])
    _attach_docking(analysed, docking_results)

    report_progress("Calculating Screening Priority Scores...")
    started = time.perf_counter()
    ranking = rank_compounds(analysed)
    _attach_ranking(analysed, ranking["ranked_compounds"])
    ranked_count = sum(entry["rank"] is not None for entry in ranking["ranked_compounds"])
    timings["ranking"] = time.perf_counter() - started

    report_progress("Building explanations...")
    started = time.perf_counter()
    explanations = explain_rankings(ranking["ranked_compounds"], analysed)
    timings["explainability"] = time.perf_counter() - started

    docked = sum(result["docking_success"] for result in docking_results)
    summary = {
        "requested": n_compounds,
        "loaded": len(library.compounds),
        "processed": len(analysed),
        "valid": len(valid),
        "invalid": len(analysed) - len(valid),
        "ligands_prepared": sum(result["ligand_prepared"] for result in docking_results),
        "docked": docked,
        "docking_failed": len(docking_results) - docked,
        "ranked": ranked_count,
    }
    pipeline = build_pipeline_stages(summary, len(library.issues), docking_info["error"], timings)
    response = {
        "message": outcome_message(pipeline, docking_info["error"], docked),
        "request": {"target": target, "n_compounds": n_compounds},
        "summary": summary,
        "pipeline": pipeline,
        "stages": {stage["key"]: stage["status"] for stage in pipeline},
        "timings": timings,
        "compounds": analysed,
        "docking": docking_info,
        "ranking": {"ranked_compounds": ranking["ranked_compounds"], "method": ranking["method"]},
        "explanations": explanations,
        "warnings": warnings,
    }

    # The report only reads the finished results. generate_report() never
    # raises, so a failed report cannot lose them.
    report_progress("Generating report...")
    started = time.perf_counter()
    report = generate_report(response)
    timings["report"] = time.perf_counter() - started
    response["report"] = report
    pipeline.append(_report_stage(report, timings["report"]))
    response["stages"]["report"] = pipeline[-1]["status"]
    response["message"] = outcome_message(pipeline, docking_info["error"], docked)
    if report["status"] != config.REPORT_GENERATED:
        warnings.append(f"Report {report['status'].replace('_', ' ')}: {'; '.join(report['errors'].values())}")

    timings["total"] = time.perf_counter() - run_started
    report_progress("Complete.")
    logger.info(
        "Screening finished in %.1f s: %d valid, %d docked, %d ranked.", timings["total"], len(valid), docked, ranked_count
    )
    return response
