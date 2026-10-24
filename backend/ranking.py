"""Compound ranking (placeholder).

Later this module will combine relevant scientific features (molecular
descriptors, docking scores) with machine-learning predictions to rank
compounds, report the uncertainty of each prediction and explain which
features drove each ranking.
"""

import logging
from typing import Any

import config

logger = logging.getLogger(__name__)


def rank_compounds(compounds: list[dict[str, Any]]) -> dict[str, Any]:
    """Rank compounds by predicted research priority (placeholder).

    Args:
        compounds: One dictionary per compound holding its calculated features.

    Returns:
        A status dictionary whose ``ranked_compounds`` list is always empty for
        now. No AI scores are generated until a model is trained on real data.
    """
    logger.warning(
        "Ranking requested for %d compound(s), but it is not implemented yet.",
        len(compounds),
    )
    return {
        "stage": "ranking",
        "status": config.STATUS_NOT_IMPLEMENTED,
        "message": "AI-assisted ranking will be added in a later stage.",
        "ranked_compounds": [],
    }
