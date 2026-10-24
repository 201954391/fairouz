"""Tests for the Screening Priority Score ranking.

Compound records are built by hand so each test controls exactly one
variable. The values are inputs to the formula, not claims about any
real molecule.
"""

import copy
import math

import pytest

import config
from backend import ranking

# Properties that sit inside every full-score range, so the property component is 1.
TYPICAL_PROPERTIES = {
    "molecular_weight": 300.0,
    "logp": 2.0,
    "tpsa": 60.0,
    "h_bond_donors": 1,
    "h_bond_acceptors": 4,
}


def make_compound(
    compound_id: str,
    affinity: float | None = -7.0,
    *,
    valid: bool = True,
    ligand_prepared: bool = True,
    docking_success: bool = True,
    **properties: float | None,
) -> dict:
    """A pipeline-shaped compound record with a docking entry."""
    return {
        "compound_id": compound_id,
        "name": f"Compound {compound_id}",
        "valid": valid,
        **TYPICAL_PROPERTIES,
        **properties,
        "docking": {
            "compound_id": compound_id,
            "ligand_prepared": ligand_prepared,
            "docking_success": docking_success,
            "best_affinity_kcal_mol": affinity if docking_success else None,
        },
    }


def ranked(compounds: list[dict]) -> list[dict]:
    return [entry for entry in ranking.rank_compounds(compounds)["ranked_compounds"] if entry["rank"] is not None]


# --- Contract ----------------------------------------------------------------------

def test_module_exposes_rank_compounds() -> None:
    assert callable(ranking.rank_compounds)


def test_result_structure() -> None:
    result = ranking.rank_compounds([make_compound("A")])
    entry = result["ranked_compounds"][0]

    assert result["stage"] == "ranking"
    assert result["status"] == config.STATUS_COMPLETED
    for key in (
        "rank", "compound_id", "name", "docking_affinity", "docking_component",
        "property_component", "screening_priority_score", "ranking_status", "explanation",
    ):
        assert key in entry
    assert entry["ranking_status"] == ranking.STATUS_RANKED


def test_weights_sum_to_one() -> None:
    assert config.RANKING_DOCKING_WEIGHT + config.RANKING_PROPERTY_WEIGHT == pytest.approx(1.0)
    assert config.RANKING_DOCKING_WEIGHT > config.RANKING_PROPERTY_WEIGHT  # Docking is the strongest signal.


def test_does_not_modify_input() -> None:
    compounds = [make_compound("A"), make_compound("B", valid=False)]
    original = copy.deepcopy(compounds)

    ranking.rank_compounds(compounds)

    assert compounds == original


# --- Docking component -------------------------------------------------------

def test_stronger_docking_ranks_higher_when_properties_equal() -> None:
    entries = ranked([make_compound("WEAK", -5.0), make_compound("STRONG", -8.0)])

    assert [entry["compound_id"] for entry in entries] == ["STRONG", "WEAK"]
    assert entries[0]["screening_priority_score"] > entries[1]["screening_priority_score"]


@pytest.mark.parametrize(
    ("affinity", "expected"),
    [
        (config.DOCKING_AFFINITY_WEAK, 0.0),
        (config.DOCKING_AFFINITY_STRONG, 1.0),
        (-1.0, 0.0),  # Weaker than the scale: clamped.
        (-15.0, 1.0),  # Stronger than the scale: clamped.
        ((config.DOCKING_AFFINITY_WEAK + config.DOCKING_AFFINITY_STRONG) / 2, 0.5),
    ],
)
def test_docking_component_scale(affinity: float, expected: float) -> None:
    assert ranking.docking_component(affinity) == pytest.approx(expected)


@pytest.mark.parametrize("bad", [None, math.nan, math.inf, "−7", True])
def test_docking_component_rejects_non_numbers(bad) -> None:
    with pytest.raises(ValueError):
        ranking.docking_component(bad)


def test_equal_docking_affinities_do_not_break_scoring() -> None:
    entries = ranked([make_compound(compound_id, -6.0) for compound_id in ("A", "B", "C")])

    scores = [entry["screening_priority_score"] for entry in entries]
    assert all(math.isfinite(score) for score in scores)
    assert len(set(scores)) == 1


def test_score_does_not_depend_on_other_compounds_in_batch() -> None:
    alone = ranked([make_compound("A", -6.0)])[0]["screening_priority_score"]
    with_others = ranked([make_compound("A", -6.0), make_compound("B", -9.5)])
    score_in_batch = next(e for e in with_others if e["compound_id"] == "A")["screening_priority_score"]

    assert alone == score_in_batch


# --- Property component ------------------------------------------------------

def test_property_component_is_one_inside_all_ranges() -> None:
    assert ranking.property_component(TYPICAL_PROPERTIES) == pytest.approx(1.0)


def test_property_component_is_deterministic() -> None:
    values = {**TYPICAL_PROPERTIES, "molecular_weight": 610.0, "logp": 6.1}

    assert ranking.property_component(values) == ranking.property_component(values)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (160.0, 1.0),  # Lower full-score boundary.
        (500.0, 1.0),  # Upper full-score boundary.
        (80.0, 0.0),  # Lower zero boundary.
        (700.0, 0.0),  # Upper zero boundary.
        (120.0, 0.5),  # Halfway down the lower slope.
        (600.0, 0.5),  # Halfway down the upper slope.
        (20.0, 0.0),
        (900.0, 0.0),
    ],
)
def test_property_subscore_boundaries(value: float, expected: float) -> None:
    molecular_weight_range = config.PROPERTY_SCORING_RANGES["molecular_weight"]

    assert ranking.property_subscore(value, molecular_weight_range) == pytest.approx(expected)


def test_one_sided_range_has_no_division_by_zero() -> None:
    # zero_below == full_from: values below the full range score 0 directly.
    assert ranking.property_subscore(-1.0, (0.0, 0.0, 5.0, 8.0)) == 0.0
    assert ranking.property_subscore(0.0, (0.0, 0.0, 5.0, 8.0)) == 1.0


def test_worse_properties_lower_the_score_when_docking_equal() -> None:
    entries = ranked([make_compound("HEAVY", -7.0, molecular_weight=650.0), make_compound("TYPICAL", -7.0)])

    assert [entry["compound_id"] for entry in entries] == ["TYPICAL", "HEAVY"]


# --- Score range and order -----------------------------------------------------

@pytest.mark.parametrize(
    "compound",
    [
        make_compound("BEST", -20.0),
        make_compound("WORST", 5.0, molecular_weight=10.0, logp=20.0, tpsa=500.0, h_bond_donors=20, h_bond_acceptors=40),
        make_compound("MID", -6.5, logp=5.5),
    ],
)
def test_score_stays_within_0_to_100(compound: dict) -> None:
    score = ranked([compound])[0]["screening_priority_score"]

    assert 0.0 <= score <= 100.0


def test_extremes_hit_the_ends_of_the_range() -> None:
    best = ranked([make_compound("BEST", -20.0)])[0]["screening_priority_score"]
    worst = ranked(
        [make_compound("WORST", 5.0, molecular_weight=10.0, logp=20.0, tpsa=500.0, h_bond_donors=20, h_bond_acceptors=40)]
    )[0]["screening_priority_score"]

    assert best == pytest.approx(100.0)
    assert worst == pytest.approx(0.0)


def test_ranking_order_and_ranks_start_at_one() -> None:
    entries = ranked([make_compound("C", -5.0), make_compound("A", -9.0), make_compound("B", -7.0)])

    assert [entry["compound_id"] for entry in entries] == ["A", "B", "C"]
    assert [entry["rank"] for entry in entries] == [1, 2, 3]


def test_ties_break_by_affinity_then_compound_id() -> None:
    # Both affinities are weaker than the scale, so both docking components are 0
    # and the scores tie exactly.
    entries = ranked([make_compound("B", -2.0), make_compound("C", -3.0), make_compound("A", -2.0)])

    assert entries[0]["screening_priority_score"] == entries[1]["screening_priority_score"]
    assert [entry["compound_id"] for entry in entries] == ["C", "A", "B"]


def test_identical_inputs_give_identical_rankings() -> None:
    compounds = [make_compound("A", -6.2), make_compound("B", -7.7, logp=5.8), make_compound("C", -4.1)]

    assert ranking.rank_compounds(compounds) == ranking.rank_compounds(copy.deepcopy(compounds))


def test_input_order_does_not_change_ranking() -> None:
    compounds = [make_compound("A", -6.2), make_compound("B", -7.7), make_compound("C", -4.1)]

    forward = [entry["compound_id"] for entry in ranked(compounds)]
    backward = [entry["compound_id"] for entry in ranked(list(reversed(compounds)))]
    assert forward == backward


# --- Compounds that must not be ranked ---------------------------------------

@pytest.mark.parametrize(
    ("compound", "status"),
    [
        (make_compound("X", valid=False), ranking.STATUS_INVALID),
        (make_compound("X", ligand_prepared=False, docking_success=False), ranking.STATUS_LIGAND_FAILED),
        (make_compound("X", docking_success=False), ranking.STATUS_DOCKING_FAILED),
        ({**make_compound("X"), "docking": None}, ranking.STATUS_NOT_DOCKED),
        (make_compound("X", affinity=None), ranking.STATUS_MISSING_DATA),
        (make_compound("X", logp=None), ranking.STATUS_MISSING_DATA),
        (make_compound("X", tpsa=math.nan), ranking.STATUS_MISSING_DATA),
    ],
)
def test_unrankable_compounds_get_status_and_no_score(compound: dict, status: str) -> None:
    entry = ranking.rank_compounds([compound])["ranked_compounds"][0]

    assert entry["ranking_status"] == status
    assert entry["rank"] is None
    assert entry["screening_priority_score"] is None
    assert entry["docking_component"] is None
    assert entry["property_component"] is None


def test_unranked_compounds_follow_ranked_ones_and_do_not_take_ranks() -> None:
    compounds = [
        make_compound("BAD", valid=False),
        make_compound("OK2", -6.0),
        make_compound("FAIL", docking_success=False),
        make_compound("OK1", -8.0),
    ]
    entries = ranking.rank_compounds(compounds)["ranked_compounds"]

    assert [entry["compound_id"] for entry in entries] == ["OK1", "OK2", "BAD", "FAIL"]
    assert [entry["rank"] for entry in entries] == [1, 2, None, None]


# --- Edge cases ------------------------------------------------------------------

def test_empty_list() -> None:
    result = ranking.rank_compounds([])

    assert result["ranked_compounds"] == []
    assert result["status"] == config.STATUS_COMPLETED


def test_single_compound() -> None:
    entries = ranked([make_compound("ONLY", -7.0)])

    assert len(entries) == 1
    assert entries[0]["rank"] == 1
    assert 0.0 <= entries[0]["screening_priority_score"] <= 100.0


# --- Explanations ----------------------------------------------------------------

def test_explanations_are_rule_based_and_specific() -> None:
    strong = ranked([make_compound("S", -9.0)])[0]["explanation"]
    weak_heavy = ranked([make_compound("W", -5.0, molecular_weight=650.0)])[0]["explanation"]

    assert strong.startswith("Strong predicted docking affinity (-9.0 kcal/mol)")
    assert "all molecular properties within typical screening ranges" in strong
    assert weak_heavy.startswith("Weak predicted docking affinity")
    assert "high molecular weight" in weak_heavy


def test_explanations_avoid_efficacy_language() -> None:
    compounds = [make_compound("A", -9.0), make_compound("B", -6.5, logp=6.0), make_compound("C", -3.0, tpsa=180.0)]

    for entry in ranked(compounds):
        text = entry["explanation"].lower()
        for word in ("effective", "efficacy", "cure", "treat", "probability", "clinical"):
            assert word not in text
