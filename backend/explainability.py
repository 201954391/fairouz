"""Explainability for the Screening Priority Score.

Explains *how Fairouz calculated* each compound's ranking. It does not
explain biological efficacy, safety, toxicity or mechanism of action.

Every number here is read from the ranking entries produced by
``backend.ranking``: docking points, property points, sub-scores and the
score itself. Nothing is recalculated with a different formula. The only
derived numbers are a property's share of the property points
(``sub-score × 30 / 5``) and points not earned (maximum minus earned), which
are simple rearrangements of the same formula.

There is no machine-learning model, so no ML explainability method (such as
SHAP) applies; the ranking formula itself is the explanation. All wording
comes from fixed rules, so the same ranking always gives the same
explanation.
"""

from collections.abc import Mapping, Sequence
from typing import Any

import config
from backend import ranking

PROPERTY_KEYS = tuple(config.PROPERTY_SCORING_RANGES)
PER_PROPERTY_MAX_POINTS = ranking.PROPERTY_MAX_POINTS / len(PROPERTY_KEYS)

DISPLAY_DECIMALS = {key: descriptor["decimals"] for key, descriptor in config.DESCRIPTORS.items()}

POSITION_WITHIN = "within configured range"
POSITION_BELOW = "below configured range"
POSITION_ABOVE = "above configured range"

DRIVER_DOCKING = "Docking"
DRIVER_PROPERTIES = "Molecular properties"
DRIVER_BOTH = "Both"
DRIVER_NONE = "None (no points earned)"

DOCKING_NOTE = (
    "More negative AutoDock Vina affinities generally correspond to stronger "
    "predicted binding. The affinity is a computational estimate, not "
    "experimental binding evidence."
)
DRIVER_NOTE = (
    "Describes where this compound's points came from in Fairouz's formula. "
    "It is not causal evidence."
)
PRIORITISATION_NOTE = "The score is intended only for computational prioritisation."

# Pipeline stage at which each unranked status arises.
FAILED_STAGES = {
    ranking.STATUS_INVALID: "SMILES validation (RDKit)",
    ranking.STATUS_LIGAND_FAILED: "3D ligand preparation (RDKit and Meeko)",
    ranking.STATUS_DOCKING_FAILED: "Molecular docking (AutoDock Vina)",
    ranking.STATUS_NOT_DOCKED: "Molecular docking (did not run)",
    ranking.STATUS_MISSING_DATA: "Ranking input check",
}


# --- Formatting --------------------------------------------------------------------

def format_score(score: float) -> str:
    """Scores are displayed to one decimal place."""
    return f"{score:.1f}"


def _format_number(value: float) -> str:
    return f"{value:g}"


def _ordinal(number: int) -> str:
    suffix = "th" if 10 <= number % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(number % 10, "th")
    return f"{number}{suffix}"


def format_property_value(key: str, value: float) -> str:
    """E.g. ``"206.3 g/mol"`` or ``"1"``."""
    unit = config.DESCRIPTORS[key]["unit"]
    decimals = DISPLAY_DECIMALS[key]
    rounded = round(value, decimals) + 0.0  # + 0.0 turns -0.0 into 0.0, avoiding "-0.00".
    text = f"{rounded:.{decimals}f}"
    return f"{text} {unit}".strip()


def range_text(key: str) -> str:
    """The configured full-score range, e.g. ``"160–500 g/mol"`` or ``"≤ 140 Å²"``."""
    zero_below, full_from, full_to, _ = config.PROPERTY_SCORING_RANGES[key]
    unit = config.DESCRIPTORS[key]["unit"]
    one_sided = zero_below == full_from == 0
    text = f"≤ {_format_number(full_to)}" if one_sided else f"{_format_number(full_from)} to {_format_number(full_to)}"
    return f"{text} {unit}".strip()


# --- Properties ----------------------------------------------------------------------

def property_position(key: str, value: float) -> str:
    """Whether a value is within, below or above its configured full-score range."""
    _, full_from, full_to, _ = config.PROPERTY_SCORING_RANGES[key]
    if value < full_from:
        return POSITION_BELOW
    if value > full_to:
        return POSITION_ABOVE
    return POSITION_WITHIN


def explain_property(key: str, value: float, subscore: float) -> dict[str, Any]:
    """Describe one descriptor's value, range, sub-score and points.

    Args:
        key: Descriptor key, e.g. ``"molecular_weight"``.
        value: The RDKit value.
        subscore: The sub-score from the ranking entry (0–1).
    """
    label = config.DESCRIPTORS[key]["label"]
    points = PER_PROPERTY_MAX_POINTS * subscore
    points_lost = PER_PROPERTY_MAX_POINTS - points
    position = property_position(key, value)
    display_value = format_property_value(key, value)

    if position == POSITION_WITHIN:
        interpretation = f"{label}: {display_value} — within configured screening range"
    else:
        direction = "below" if position == POSITION_BELOW else "above"
        interpretation = (
            f"{label}: {display_value} — {direction} configured screening range, "
            f"reducing the property score by {points_lost:.1f} points"
        )
    return {
        "key": key,
        "label": label,
        "value": value,
        "display_value": display_value,
        "range_text": range_text(key),
        "position": position,
        "subscore": subscore,
        "points": points,
        "max_points": PER_PROPERTY_MAX_POINTS,
        "points_lost": points_lost,
        "interpretation": interpretation,
    }


# --- Docking and drivers ----------------------------------------------------------

def explain_docking(affinity: float, component: float, points: float) -> dict[str, Any]:
    """Describe the docking contribution exactly as the ranking calculated it."""
    return {
        "affinity": affinity,
        "component": component,
        "points": points,
        "max_points": ranking.DOCKING_MAX_POINTS,
        "points_lost": ranking.DOCKING_MAX_POINTS - points,
        "text": (
            f"Predicted docking affinity: {affinity:.1f} kcal/mol. On Fairouz's fixed "
            f"docking scale, this contributes {points:.1f} points out of a possible "
            f"{ranking.DOCKING_MAX_POINTS:.0f}."
        ),
        "note": DOCKING_NOTE,
    }


def primary_driver(docking_points: float, property_points: float) -> tuple[str, float | None]:
    """Which part of the formula supplied most of the points.

    Returns:
        The driver label and docking's share of the points (None if no points).
    """
    total = docking_points + property_points
    if total <= 0:
        return DRIVER_NONE, None
    docking_share = docking_points / total
    if docking_share >= config.EXPLANATION_DRIVER_SHARE:
        return DRIVER_DOCKING, docking_share
    if docking_share <= 1 - config.EXPLANATION_DRIVER_SHARE:
        return DRIVER_PROPERTIES, docking_share
    return DRIVER_BOTH, docking_share


def contribution_chart_rows(docking_points: float, property_points: float) -> list[dict[str, Any]]:
    """Rows for a stacked bar: points earned and not earned for each component."""
    rows = []
    for component, earned, maximum in (
        ("Docking", docking_points, ranking.DOCKING_MAX_POINTS),
        ("Molecular properties", property_points, ranking.PROPERTY_MAX_POINTS),
    ):
        rows.append({"component": component, "segment": "Earned", "points": earned})
        rows.append({"component": component, "segment": "Not earned", "points": maximum - earned})
    return rows


# --- Strengths and limitations ---------------------------------------------------

def _strengths(docking: Mapping[str, Any], properties: Sequence[Mapping[str, Any]]) -> list[tuple[str, str]]:
    """Strengths as ``(group, text)`` pairs; group is "docking" or "properties"."""
    strengths = []
    if docking["component"] >= config.EXPLANATION_DOCKING_STRENGTH_COMPONENT:
        strengths.append((
            "docking",
            f"Stronger predicted docking affinity ({docking['affinity']:.1f} kcal/mol, "
            f"{docking['points']:.1f} of {docking['max_points']:.0f} points).",
        ))
    within = [item for item in properties if item["position"] == POSITION_WITHIN]
    if len(within) == len(properties):
        strengths.append((
            "properties",
            f"All five molecular-property descriptors fall within the configured screening "
            f"ranges ({ranking.PROPERTY_MAX_POINTS:.0f} of {ranking.PROPERTY_MAX_POINTS:.0f} points).",
        ))
    elif within:
        labels = ", ".join(item["label"] for item in within)
        strengths.append(("properties", f"{len(within)} of {len(properties)} descriptors within configured screening ranges: {labels}."))
    return strengths


def _limiting_factors(docking: Mapping[str, Any], properties: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Every factor that lost points, largest loss first (ties keep formula order)."""
    factors = []
    if docking["points_lost"] > 1e-9:
        weak = docking["component"] < config.EXPLANATION_DOCKING_STRENGTH_COMPONENT
        phrase = (
            "relatively weak predicted docking affinity"
            if weak
            else "predicted docking affinity short of the top of Fairouz's fixed scale"
        )
        factors.append(
            {
                "factor": "docking",
                "points_lost": docking["points_lost"],
                "phrase": f"{phrase} ({docking['points_lost']:.1f} of {docking['max_points']:.0f} docking points not earned)",
                "text": (
                    f"{phrase.capitalize()} ({docking['affinity']:.1f} kcal/mol): "
                    f"{docking['points_lost']:.1f} of {docking['max_points']:.0f} docking points not earned."
                ),
            }
        )
    for item in properties:
        if item["points_lost"] > 1e-9:
            factors.append(
                {
                    "factor": item["key"],
                    "points_lost": item["points_lost"],
                    "phrase": f"{item['label']} {item['position']} ({item['points_lost']:.1f} of {item['max_points']:.0f} points not earned)",
                    "text": (
                        f"{item['label']} {item['position']} ({item['display_value']}; configured range "
                        f"{item['range_text']}): {item['points_lost']:.1f} of {item['max_points']:.0f} points not earned."
                    ),
                }
            )
    # sorted() is stable, so equal losses keep docking-first, then config order.
    return sorted(factors, key=lambda factor: -factor["points_lost"])


def _main_positive_factor(docking: Mapping[str, Any], properties: Sequence[Mapping[str, Any]], property_points: float) -> str:
    if docking["points"] <= 0 and property_points <= 0:
        return "None: the compound earned no points."
    if property_points >= docking["points"]:
        within = sum(item["position"] == POSITION_WITHIN for item in properties)
        if within == len(properties):
            return "All five molecular-property descriptors fall within the configured screening ranges."
        return (
            f"Its molecular-property profile ({within} of {len(properties)} descriptors within "
            f"configured ranges; {property_points:.1f} of {ranking.PROPERTY_MAX_POINTS:.0f} points)."
        )
    return (
        f"Its predicted docking affinity ({docking['affinity']:.1f} kcal/mol; "
        f"{docking['points']:.1f} of {docking['max_points']:.0f} points)."
    )


# --- Comparisons and ties ----------------------------------------------------------

def explain_tie(higher: Mapping[str, Any], lower: Mapping[str, Any]) -> str | None:
    """Explain the tie-break when two neighbours show the same (displayed) score."""
    shown = format_score(higher["screening_priority_score"])
    if shown != format_score(lower["screening_priority_score"]):
        return None
    if higher["screening_priority_score"] != lower["screening_priority_score"]:
        return (
            f"{higher['name']} and {lower['name']} have the same displayed Screening Priority "
            f"Score ({shown}), but {higher['name']}'s unrounded score is higher "
            f"({higher['screening_priority_score']:.3f} vs {lower['screening_priority_score']:.3f})."
        )
    if higher["docking_affinity"] != lower["docking_affinity"]:
        return (
            f"{higher['name']} and {lower['name']} have the same Screening Priority Score ({shown}). "
            f"{higher['name']} ranks higher because its docking affinity is more negative "
            f"({higher['docking_affinity']:.1f} vs {lower['docking_affinity']:.1f} kcal/mol); "
            "compound ID is used only as the final deterministic tie-break."
        )
    return (
        f"{higher['name']} and {lower['name']} have the same Screening Priority Score ({shown}) "
        f"and the same docking affinity, so compound ID ({higher['compound_id']} before "
        f"{lower['compound_id']}) is used as the final deterministic tie-break."
    )


def compare_ranked(higher: Mapping[str, Any], lower: Mapping[str, Any]) -> str:
    """Explain why ``higher`` ranks above ``lower``, using their contributions only."""
    if higher["screening_priority_score"] == lower["screening_priority_score"]:
        return explain_tie(higher, lower)

    upper, other = higher["name"], lower["name"]
    docking_gap = higher["docking_points"] - lower["docking_points"]
    property_gap = higher["property_points"] - lower["property_points"]
    similar = config.EXPLANATION_SIMILAR_POINTS
    docking_detail = (
        f"{higher['docking_affinity']:.1f} vs {lower['docking_affinity']:.1f} kcal/mol; "
        f"{higher['docking_points']:.1f} vs {lower['docking_points']:.1f} docking points"
    )
    property_detail = f"{higher['property_points']:.1f} vs {lower['property_points']:.1f} property points"

    if abs(docking_gap) < similar and abs(property_gap) < similar:
        total_gap = higher["screening_priority_score"] - lower["screening_priority_score"]
        return (
            f"{upper} and {other} have very similar docking and property contributions; "
            f"{upper} scores {total_gap:.2f} points higher in total."
        )
    if abs(docking_gap) < similar:
        if higher["docking_points"] == lower["docking_points"] == 0:
            docking_text = (
                f"Both compounds earn 0 docking points: their affinities "
                f"({higher['docking_affinity']:.1f} and {lower['docking_affinity']:.1f} kcal/mol) are at or "
                f"weaker than the {config.DOCKING_AFFINITY_WEAK:.0f} kcal/mol end of Fairouz's fixed scale"
            )
        else:
            docking_text = f"Both compounds have similar docking contributions ({docking_detail})"
        return f"{docking_text}, so {upper} ranks higher because of its stronger molecular-property profile ({property_detail})."
    if docking_gap >= property_gap:
        text = f"{upper} ranks above {other} mainly because its predicted docking affinity is more favourable ({docking_detail})."
        if property_gap <= -similar:
            text += f" This outweighs {other}'s stronger molecular-property profile ({property_detail})."
        return text
    text = f"{upper} ranks above {other} mainly because it has a stronger molecular-property profile ({property_detail})."
    if docking_gap <= -similar:
        text += f" This outweighs {other}'s more favourable predicted docking affinity ({docking_detail})."
    return text


# --- Whole explanations ------------------------------------------------------------

def _docking_position_phrase(entry: Mapping[str, Any], ranked: Sequence[Mapping[str, Any]]) -> str:
    """Where this compound's affinity sits among the ranked compounds (ties share a place)."""
    position = 1 + sum(other["docking_affinity"] < entry["docking_affinity"] for other in ranked)
    if position == 1:
        return "the most negative predicted docking affinity in the screened set"
    return f"the {_ordinal(position)} most negative predicted docking affinity of {len(ranked)} ranked compounds"


def _summary(entry: Mapping[str, Any], ranked: Sequence[Mapping[str, Any]], properties: Sequence[Mapping[str, Any]], limiting: str) -> str:
    within = sum(item["position"] == POSITION_WITHIN for item in properties)
    property_phrase = (
        "a favourable descriptor profile (all five descriptors within configured ranges)"
        if within == len(properties)
        else f"{within} of {len(properties)} descriptors within configured ranges"
    )
    return (
        f"{entry['name']} ranks #{entry['rank']} of {len(ranked)} with a Screening Priority Score of "
        f"{format_score(entry['screening_priority_score'])}. It has "
        f"{_docking_position_phrase(entry, ranked)} and {property_phrase}. "
        f"Main limiting factor: {limiting}. {PRIORITISATION_NOTE}"
    )


def explain_ranked(
    entry: Mapping[str, Any],
    compound: Mapping[str, Any],
    ranked: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Build the full explanation for one ranked compound.

    Args:
        entry: The compound's ranking entry.
        compound: The compound's pipeline record (for descriptor values).
        ranked: All ranked entries in rank order (for comparisons).
    """
    docking = explain_docking(entry["docking_affinity"], entry["docking_component"], entry["docking_points"])
    properties = [explain_property(key, compound[key], entry["property_subscores"][key]) for key in PROPERTY_KEYS]
    limiting = _limiting_factors(docking, properties)
    main_limiting = limiting[0]["phrase"] if limiting else "none: the compound earned the maximum score"
    driver, docking_share = primary_driver(entry["docking_points"], entry["property_points"])
    strengths = _strengths(docking, properties)
    # Mirrors _main_positive_factor: the larger contribution is the main positive factor.
    main_positive_group = "properties" if entry["property_points"] >= docking["points"] else "docking"

    index = entry["rank"] - 1
    above = ranked[index - 1] if index > 0 else None
    below = ranked[index + 1] if index + 1 < len(ranked) else None
    tie_notes = [note for note in (above and explain_tie(above, entry), below and explain_tie(entry, below)) if note]

    return {
        "compound_id": entry["compound_id"],
        "name": entry["name"],
        "ranked": True,
        "rank": entry["rank"],
        "total_ranked": len(ranked),
        "screening_priority_score": entry["screening_priority_score"],
        "score_display": format_score(entry["screening_priority_score"]),
        "docking": docking,
        "properties": {
            "component": entry["property_component"],
            "points": entry["property_points"],
            "max_points": ranking.PROPERTY_MAX_POINTS,
            "items": properties,
        },
        "chart_rows": contribution_chart_rows(entry["docking_points"], entry["property_points"]),
        "primary_driver": driver,
        "docking_share": docking_share,
        "driver_note": DRIVER_NOTE,
        "strengths": [text for _, text in strengths],
        "limitations": [factor["text"] for factor in limiting],
        # The same lists without the main factor, so displays do not repeat it.
        "other_strengths": [text for group, text in strengths if group != main_positive_group],
        "other_limitations": [factor["text"] for factor in limiting[1:]],
        "main_positive_factor": _main_positive_factor(docking, properties, entry["property_points"]),
        "main_limiting_factor": main_limiting[0].upper() + main_limiting[1:],
        "comparison_above": compare_ranked(above, entry) if above else None,
        "comparison_below": compare_ranked(entry, below) if below else None,
        "tie_notes": tie_notes,
        "summary": _summary(entry, ranked, properties, main_limiting),
    }


def explain_unranked(entry: Mapping[str, Any], compound: Mapping[str, Any] | None) -> dict[str, Any]:
    """Explain why a compound has no score and where the pipeline stopped."""
    status = entry["ranking_status"]
    stage = FAILED_STAGES.get(status, "Unknown stage")
    compound = compound or {}
    docking = compound.get("docking") or {}
    detail = compound.get("error") or docking.get("error") or entry.get("explanation") or "No further detail recorded."
    return {
        "compound_id": entry["compound_id"],
        "name": entry["name"],
        "ranked": False,
        "rank": None,
        "status": status,
        "failed_stage": stage,
        "detail": detail,
        "summary": (
            f"{entry['name']} was not ranked: {status.lower()} at the {stage} stage. "
            "No Screening Priority Score was calculated for it, and none was fabricated."
        ),
    }


def explain_rankings(
    ranking_entries: Sequence[Mapping[str, Any]],
    compounds: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Explain every compound's ranking, in the same order as ``ranking_entries``.

    Args:
        ranking_entries: ``ranked_compounds`` from ``backend.ranking.rank_compounds``.
        compounds: The pipeline compound records (descriptor values, errors).
            Neither argument is modified.

    Returns:
        One explanation dictionary per entry: a full breakdown for ranked
        compounds, and a failure explanation for unranked ones.
    """
    compounds_by_id = {compound.get("compound_id"): compound for compound in compounds}
    ranked = sorted((entry for entry in ranking_entries if entry["rank"] is not None), key=lambda entry: entry["rank"])
    return [
        explain_ranked(entry, compounds_by_id[entry["compound_id"]], ranked)
        if entry["rank"] is not None
        else explain_unranked(entry, compounds_by_id.get(entry["compound_id"]))
        for entry in ranking_entries
    ]
