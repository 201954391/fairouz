"""Research report generation (placeholder).

Later this module will build a research-style summary from actual screening
results: target, data provenance, top-ranked compounds, the reasons they
ranked highly, uncertainty and limitations, and a clear statement that the
findings are computational predictions requiring laboratory validation.
"""

import logging
from typing import Any

import config

logger = logging.getLogger(__name__)


def generate_report(screening_results: dict[str, Any] | None = None) -> dict[str, Any]:
    """Generate a research report from screening results (placeholder).

    Args:
        screening_results: Output of the screening pipeline, once implemented.

    Returns:
        A status dictionary whose ``report`` is ``None`` for now.
    """
    logger.warning("Report generation requested, but it is not implemented yet.")
    return {
        "stage": "report",
        "status": config.STATUS_NOT_IMPLEMENTED,
        "message": (
            "Report generation is not implemented yet. A report will be "
            "generated from real screening results in a later stage."
        ),
        "report": None,
    }
