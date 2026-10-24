"""Molecular descriptor calculation (placeholder).

RDKit will be integrated here in Step 2 to calculate properties such as:

* molecular weight
* LogP (lipophilicity)
* TPSA (topological polar surface area)
* hydrogen bond donors
* hydrogen bond acceptors
"""

import logging
from typing import Any

import config

logger = logging.getLogger(__name__)


def calculate_descriptors(smiles: list[str]) -> dict[str, Any]:
    """Calculate molecular descriptors for each compound (placeholder).

    Args:
        smiles: Compounds as SMILES strings.

    Returns:
        A status dictionary whose ``descriptors`` list is always empty for now.
        No descriptor values are produced until RDKit is integrated.
    """
    logger.warning(
        "Descriptor calculation requested for %d compound(s), but it is not implemented yet.",
        len(smiles),
    )
    return {
        "stage": "descriptors",
        "status": config.STATUS_NOT_IMPLEMENTED,
        "message": "Descriptor calculation (RDKit) will be added in the next stage.",
        "descriptors": [],
    }
