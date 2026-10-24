"""Screening Priority Score: transparent, evidence-weighted compound ranking.

This is a documented heuristic, **not** a trained model. With only 10
compounds and no labelled activity data, a supervised model would learn
nothing meaningful, so Fairouz does not pretend to have one.

For each compound that was docked successfully:

* **Docking component** (0–1): its Vina affinity on a fixed scale.
  ``config.DOCKING_AFFINITY_WEAK`` (−4 kcal/mol) or weaker scores 0;
  ``config.DOCKING_AFFINITY_STRONG`` (−10 kcal/mol) or stronger scores 1;
  linear in between.
* **Property component** (0–1): the mean of five property sub-scores. Each
  is 1 inside a typical screening range and falls linearly to 0 outside it
  (``config.PROPERTY_SCORING_RANGES``).
* **Screening Priority Score** (0–100) =
  100 × (0.7 × docking component + 0.3 × property component).

The score only prioritises compounds for further investigation. It is not
a probability of efficacy, success or clinical benefit.

Compounds that are invalid, were not docked or failed docking are listed
with an explicit status and no score. Identical inputs always give
identical rankings: there is no randomness, and ties are broken by docking
affinity, then compound ID.
"""

import logging
import math
from collections.abc import Mapping, Sequence
from typing import Any

import config

logger = logging.getLogger(__name__)

STATUS_RANKED = "Ranked"
STATUS_INVALID = "Invalid molecule"
STATUS_NOT_DOCKED = "Not docked"
STATUS_LIGAND_FAILED = "Ligand preparation failed"
STATUS_DOCKING_FAILED = "Docking failed"
STATUS_MISSING_DATA = "Not ranked: missing data"

# Short phrases used in explanations when a property is outside its full-score range.
PROPERTY_PHRASES = {
    "molecular_weight": "molecular weight",
    "logp": "LogP",
    "tpsa": "polar surface area",
    "h_bond_donors": "hydrogen-bond donor count",
    "h_bond_acceptors": "hydrogen-bond acceptor count",
}


def _is_number(value: Any) -> bool:
    """True for a finite int or float (bools and NaN excluded)."""
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


# --- Components ----------------------------------------------------------------

def docking_component(affinity: float) -> float:
    """Map a Vina affinity (kcal/mol) onto 0–1, where more negative is higher.

    Raises:
        ValueError: If the affinity is not a finite number.
    """
    if not _is_number(affinity):
        raise ValueError(f"Docking affinity must be a finite number, got {affinity!r}.")
    weak, strong = config.DOCKING_AFFINITY_WEAK, config.DOCKING_AFFINITY_STRONG
    fraction = (weak - affinity) / (weak - strong)
    return min(1.0, max(0.0, fraction))


def property_subscore(value: float, scoring_range: tuple[float, float, float, float]) -> float:
    """Score one property: 1 inside its full range, falling linearly to 0 outside.

    Args:
        value: The property value.
        scoring_range: ``(zero_below, full_from, full_to, zero_above)``.

    Returns:
        A score from 0 to 1. Boundary values (exactly ``full_from`` or
        ``full_to``) score 1; values at or beyond a zero limit score 0.
    """
    zero_below, full_from, full_to, zero_above = scoring_range
    if full_from <= value <= full_to:
        return 1.0
    if value < full_from:
        if value <= zero_below:
            return 0.0
        return (value - zero_below) / (full_from - zero_below)
    if value >= zero_above:
        return 0.0
    return (zero_above - value) / (zero_above - full_to)


def property_subscores(properties: Mapping[str, Any]) -> dict[str, float]:
    """Score every configured property.

    Raises:
        ValueError: If any property is missing or not a finite number.
    """
    missing = [key for key in config.PROPERTY_SCORING_RANGES if not _is_number(properties.get(key))]
    if missing:
        raise ValueError(f"Missing or non-numeric properties: {', '.join(missing)}.")
    return {
        key: property_subscore(float(properties[key]), scoring_range)
        for key, scoring_range in config.PROPERTY_SCORING_RANGES.items()
    }


def _mean(subscores: Mapping[str, float]) -> float:
    return sum(subscores.values()) / len(subscores)


def property_component(properties: Mapping[str, Any]) -> float:
    """Mean of the property sub-scores, from 0 to 1."""
    return _mean(property_subscores(properties))


DOCKING_MAX_POINTS = 100.0 * config.RANKING_DOCKING_WEIGHT
PROPERTY_MAX_POINTS = 100.0 * config.RANKING_PROPERTY_WEIGHT


def contribution_points(docking: float, properties: float) -> tuple[float, float]:
    """Convert the two components into score points (out of 70 and 30 by default)."""
    return DOCKING_MAX_POINTS * docking, PROPERTY_MAX_POINTS * properties


def screening_priority_score(docking: float, properties: float) -> float:
    """Combine the two components into a 0–100 score: the sum of their points."""
    docking_points, property_points = contribution_points(docking, properties)
    return docking_points + property_points


# --- Explanations --------------------------------------------------------------

def _docking_phrase(affinity: float) -> str:
    if affinity <= config.DOCKING_STRONG_THRESHOLD:
        return "Strong predicted docking affinity"
    if affinity <= config.DOCKING_MODERATE_THRESHOLD:
        return "Moderate predicted docking affinity"
    return "Weak predicted docking affinity"


def _outside_range_phrases(properties: Mapping[str, Any], subscores: Mapping[str, float]) -> list[str]:
    """E.g. ``["low molecular weight", "high LogP"]`` for properties below full score."""
    phrases = []
    for key, score in subscores.items():
        if score < 1.0:
            full_from = config.PROPERTY_SCORING_RANGES[key][1]
            direction = "low" if properties[key] < full_from else "high"
            phrases.append(f"{direction} {PROPERTY_PHRASES[key]}")
    return phrases


def explain_ranking(affinity: float, properties: Mapping[str, Any], subscores: Mapping[str, float]) -> str:
    """Write a short rule-based explanation of a compound's score."""
    docking_text = f"{_docking_phrase(affinity)} ({affinity:.1f} kcal/mol)"
    outside = _outside_range_phrases(properties, subscores)
    if not outside:
        return f"{docking_text}, with all molecular properties within typical screening ranges."
    listed = outside[0] if len(outside) == 1 else f"{', '.join(outside[:-1])} and {outside[-1]}"
    verb = "reduces" if len(outside) == 1 else "reduce"
    return f"{docking_text}; {listed} (outside typical screening ranges) {verb} its overall priority."


# --- Ranking ----------------------------------------------------------------------

def _base_entry(compound: Mapping[str, Any]) -> dict[str, Any]:
    """Ranking entry with identity filled in and every score empty."""
    docking = compound.get("docking") or {}
    return {
        "rank": None,
        "compound_id": compound.get("compound_id", ""),
        "name": compound.get("name", ""),
        "docking_affinity": docking.get("best_affinity_kcal_mol"),
        "docking_component": None,
        "property_component": None,
        "property_subscores": None,
        "docking_points": None,
        "property_points": None,
        "screening_priority_score": None,
        "ranking_status": STATUS_RANKED,
        "explanation": None,
    }


def _unranked_status(compound: Mapping[str, Any]) -> str | None:
    """Why a compound cannot be ranked, or None if it can."""
    if not compound.get("valid"):
        return STATUS_INVALID
    docking = compound.get("docking")
    if docking is None:
        return STATUS_NOT_DOCKED
    if not docking.get("ligand_prepared"):
        return STATUS_LIGAND_FAILED
    if not docking.get("docking_success"):
        return STATUS_DOCKING_FAILED
    return None


def score_compound(compound: Mapping[str, Any]) -> dict[str, Any]:
    """Score one compound, or explain why it is not ranked.

    Args:
        compound: A descriptor result with a ``docking`` entry, as produced
            by the screening pipeline. It is not modified.

    Returns:
        A ranking entry; ``rank`` is assigned later by ``rank_compounds``.
    """
    entry = _base_entry(compound)
    status = _unranked_status(compound)
    if status is not None:
        entry["ranking_status"] = status
        entry["explanation"] = f"Not ranked: {status.lower()}."
        return entry

    try:
        dock = docking_component(entry["docking_affinity"])
        subscores = property_subscores(compound)
    except ValueError as error:
        logger.warning("Compound %s not ranked: %s", entry["compound_id"], error)
        entry["ranking_status"] = STATUS_MISSING_DATA
        entry["explanation"] = f"Not ranked: {error}"
        return entry

    properties = _mean(subscores)
    docking_points, property_points = contribution_points(dock, properties)
    entry.update(
        docking_component=dock,
        property_component=properties,
        property_subscores=subscores,
        docking_points=docking_points,
        property_points=property_points,
        # The score is exactly the sum of the points the explanation layer shows.
        screening_priority_score=docking_points + property_points,
        explanation=explain_ranking(entry["docking_affinity"], compound, subscores),
    )
    return entry


def _sort_key(entry: Mapping[str, Any]) -> tuple[float, float, str]:
    """Highest score first; ties go to the more negative affinity, then compound ID."""
    return (-entry["screening_priority_score"], entry["docking_affinity"], entry["compound_id"])


def rank_compounds(compounds: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Rank compounds by Screening Priority Score.

    Args:
        compounds: Pipeline compound records (descriptors plus ``docking``).
            Invalid and undocked compounds may be included; they are listed
            as not ranked. The input is not modified.

    Returns:
        ``stage``, ``status``, ``ranked_compounds`` (ranked entries from rank 1
        down, followed by unranked entries in input order) and ``method``
        (the weights and scale used).
    """
    entries = [score_compound(compound) for compound in compounds]
    ranked = sorted((entry for entry in entries if entry["ranking_status"] == STATUS_RANKED), key=_sort_key)
    for rank, entry in enumerate(ranked, start=1):
        entry["rank"] = rank
    unranked = [entry for entry in entries if entry["ranking_status"] != STATUS_RANKED]

    logger.info("Ranked %d compound(s); %d not ranked.", len(ranked), len(unranked))
    return {
        "stage": "ranking",
        "status": config.STATUS_COMPLETED,
        "ranked_compounds": ranked + unranked,
        "method": {
            "name": "Screening Priority Score (transparent heuristic, not a trained model)",
            "docking_weight": config.RANKING_DOCKING_WEIGHT,
            "property_weight": config.RANKING_PROPERTY_WEIGHT,
            "docking_scale_kcal_mol": {
                "zero_at": config.DOCKING_AFFINITY_WEAK,
                "full_at": config.DOCKING_AFFINITY_STRONG,
            },
            "property_ranges": config.PROPERTY_SCORING_RANGES,
        },
    }
