"""Screening pipeline coordinator.

Orchestrates the virtual screening workflow:

1. Validate the request.
2. Load the compound library (``backend.compounds``).          Implemented
3. Calculate molecular descriptors (``backend.descriptors``).  Implemented
4. Prepare 3D ligands and dock them (``backend.docking``).     Implemented
5. Rank compounds (``backend.ranking``).                       Not implemented
6. Generate a report (``backend.report``).                     Not implemented

If docking cannot run at all (for example Vina is not installed), the
descriptors are still returned and the docking stage is reported as
unavailable with the reason. Unimplemented stages are still called, so the
response reports their real status. No values are fabricated.
"""

import logging
from collections.abc import Callable
from typing import Any

import config
from backend.compounds import load_compounds
from backend.descriptors import calculate_descriptors
from backend.docking import DockingSetupError, run_docking
from backend.ranking import rank_compounds
from backend.report import generate_report

logger = logging.getLogger(__name__)

DOCKING_COMPLETE_MESSAGE = (
    "Molecular analysis and docking complete. AI ranking and the final report "
    "will be implemented in later stages."
)
DOCKING_UNAVAILABLE_MESSAGE = (
    "Molecular analysis complete, but docking could not run. AI ranking and the "
    "final report will be implemented in later stages."
)


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


def run_screening(
    target: str,
    n_compounds: int,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Run the screening pipeline as far as it is currently implemented.

    Args:
        target: Name of the viral protein target to screen against.
        n_compounds: Number of compounds to screen. If the dataset holds
            fewer, all of them are used; compounds are never duplicated.
        progress: Optional callback receiving a short message as each real
            stage (and each docked compound) starts.

    Returns:
        A response with the keys:

        * ``message``: overall outcome for the user.
        * ``request``: the target and requested count.
        * ``summary``: requested, loaded, processed, valid, invalid,
          ligands_prepared, docked and docking_failed counts.
        * ``stages``: status of each pipeline stage.
        * ``compounds``: descriptor results, each with a ``docking`` entry.
        * ``docking``: run folder and settings, or the setup error.
        * ``warnings``: skipped rows, count shortfalls and docking problems.

    Raises:
        ValueError: If the request fails validation.
        CompoundDataError: If the compound dataset cannot be used.
    """
    validate_screening_request(target, n_compounds)
    report_progress = progress or (lambda message: None)
    logger.info("Screening started: target=%r, n_compounds=%d", target, n_compounds)

    report_progress("Loading compounds...")
    library = load_compounds()
    selected = library.compounds[:n_compounds]
    warnings = list(library.issues)
    shortfall = _shortfall_warning(n_compounds, len(library.compounds))
    if shortfall:
        logger.warning(shortfall)
        warnings.append(shortfall)

    report_progress(f"Calculating molecular descriptors for {len(selected)} compounds...")
    analysed = calculate_descriptors(selected)
    valid = [compound for compound in analysed if compound["valid"]]

    report_progress(f"Preparing 3D ligands and docking {len(valid)} compounds (this takes longer)...")
    try:
        docking = run_docking(valid, progress=report_progress)
        docking_results = docking["results"]
        docking_info = {"run_dir": docking["run_dir"], "metadata": docking["metadata"], "error": None}
        docking_status = docking["status"]
    except DockingSetupError as error:
        logger.error("Docking unavailable: %s", error)
        warnings.append(f"Docking could not run: {error}")
        docking_results, docking_status = [], config.STATUS_UNAVAILABLE
        docking_info = {"run_dir": None, "metadata": None, "error": str(error)}
    _attach_docking(analysed, docking_results)

    ranking = rank_compounds(valid)
    report = generate_report()

    report_progress("Complete.")
    docked = sum(result["docking_success"] for result in docking_results)
    logger.info("Screening finished: %d valid, %d docked.", len(valid), docked)

    return {
        "message": DOCKING_COMPLETE_MESSAGE if docking_info["error"] is None else DOCKING_UNAVAILABLE_MESSAGE,
        "request": {"target": target, "n_compounds": n_compounds},
        "summary": {
            "requested": n_compounds,
            "loaded": len(library.compounds),
            "processed": len(analysed),
            "valid": len(valid),
            "invalid": len(analysed) - len(valid),
            "ligands_prepared": sum(result["ligand_prepared"] for result in docking_results),
            "docked": docked,
            "docking_failed": len(docking_results) - docked,
        },
        "stages": {
            "compound_loading": config.STATUS_COMPLETED,
            "molecular_descriptors": config.STATUS_COMPLETED,
            "docking": docking_status,
            "ai_ranking": ranking["status"],
            "report": report["status"],
        },
        "compounds": analysed,
        "docking": docking_info,
        "warnings": warnings,
    }
