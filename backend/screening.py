"""Screening pipeline coordinator.

This module will eventually orchestrate the full virtual screening workflow:

1. Load the selected viral protein target.
2. Load the compound library.
3. Calculate molecular descriptors (``backend.descriptors``).
4. Run molecular docking (``backend.docking``).
5. Rank compounds (``backend.ranking``).
6. Generate results and a report (``backend.report``).

For now it validates the request and reports that the scientific stages are
not implemented yet. It never fabricates results.
"""

import logging
from typing import Any

import config

logger = logging.getLogger(__name__)

NOT_IMPLEMENTED_MESSAGE = (
    "The screening pipeline is ready for implementation. Scientific screening "
    "modules will be added in the next development stage."
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


def run_screening(target: str, n_compounds: int) -> dict[str, Any]:
    """Run the virtual screening pipeline (placeholder).

    Args:
        target: Name of the viral protein target to screen against.
        n_compounds: Number of compounds to screen.

    Returns:
        A structured status response with the keys ``status``, ``message``,
        ``request``, ``stages`` and ``results``. ``results`` stays empty until
        the scientific modules are implemented.

    Raises:
        ValueError: If the request fails validation.
    """
    validate_screening_request(target, n_compounds)
    logger.info("Screening requested: target=%r, n_compounds=%d", target, n_compounds)
    logger.warning("Screening pipeline not implemented yet; no calculations were run.")

    return {
        "status": config.STATUS_NOT_IMPLEMENTED,
        "message": NOT_IMPLEMENTED_MESSAGE,
        "request": {"target": target, "n_compounds": n_compounds},
        "stages": [name for name, _ in config.PIPELINE_STAGES],
        "results": [],
    }
