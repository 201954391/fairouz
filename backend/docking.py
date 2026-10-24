"""Molecular docking (placeholder).

AutoDock Vina will be integrated here in a later stage. Docking estimates how
strongly a compound may bind to the target protein; the result is a predicted
binding affinity (kcal/mol), which is an estimate, not a measurement.
"""

import logging
from typing import Any

import config

logger = logging.getLogger(__name__)


def run_docking(target: str, compounds: list[str]) -> dict[str, Any]:
    """Dock compounds against a viral protein target (placeholder).

    Planned behaviour: prepare the receptor and ligands, run AutoDock Vina for
    each compound and return its predicted binding affinity and pose. Until
    then no docking scores are produced or simulated.

    Args:
        target: Name of the prepared viral protein target.
        compounds: Compound identifiers (e.g. SMILES strings) to dock.

    Returns:
        A status dictionary whose ``scores`` list is always empty for now.
    """
    logger.warning(
        "Docking requested for %d compound(s) against %r, but it is not implemented yet.",
        len(compounds),
        target,
    )
    return {
        "stage": "docking",
        "status": config.STATUS_NOT_IMPLEMENTED,
        "message": "Molecular docking (AutoDock Vina) will be added in a later stage.",
        "scores": [],
    }
