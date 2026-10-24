"""Tests for the explainability layer.

Compound records are built by hand so each test controls one variable;
their values are inputs to the formula, not claims about real molecules.
"""

import copy
import re
from pathlib import Path

import pytest

import config
from backend import explainability, ranking, screening
from backend.explainability import (
    DRIVER_BOTH,
    DRIVER_DOCKING,
    DRIVER_NONE,
    DRIVER_PROPERTIES,
    POSITION_ABOVE,
    POSITION_BELOW,
    POSITION_WITHIN,
    explain_rankings,
)

TYPICAL_PROPERTIES = {
    "molecular_weight": 300.0,
    "logp": 2.0,
    "tpsa": 60.0,
    "h_bond_donors": 1,
    "h_bond_acceptors": 4,
}

# Wording the explanations must never use (README › Explainability).
FORBIDDEN_PATTERNS = [
    r"\befficac", r"\beffective\b", r"\bsafe\b", r"\bsafety\b", r"\bclinical", r"\bcure",
    r"\btreat", r"\btherapeutic", r"\bpromising\b", r"likely to work", r"\bprobabilit",
    r"\bantiviral", r"\btoxic",
]


def make_compound(compound_id, affinity=-7.0, *, valid=True, ligand_prepared=True, docking_success=True, error=None, **properties):
    return {
        "compound_id": compound_id,
        "name": f"Compound {compound_id}",
        "smiles": "C",
        "valid": valid,
        "error": error,
        **TYPICAL_PROPERTIES,
        **properties,
        "docking": {
            "compound_id": compound_id,
            "ligand_prepared": ligand_prepared,
            "docking_success": docking_success,
            "best_affinity_kcal_mol": affinity if docking_success else None,
            "error": None if docking_success else "Docking failed: Vina exited with code 1",
        },
    }


def explain(compounds):
    entries = ranking.rank_compounds(compounds)["ranked_compounds"]
    return entries, explain_rankings(entries, compounds)


def by_id(explanations):
    return {explanation["compound_id"]: explanation for explanation in explanations}


def all_text(value):
    """Every string inside a nested explanation structure."""
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for item in value.values():
            yield from all_text(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from all_text(item)


# --- Contributions come straight from the ranking ---------------------------------

def test_contribution_values_match_ranking_exactly() -> None:
    entries, explanations = explain([make_compound("A", -6.1, molecular_weight=150.0), make_compound("B", -8.3)])

    for entry in entries:
        explanation = by_id(explanations)[entry["compound_id"]]
        assert explanation["screening_priority_score"] == entry["screening_priority_score"]
        assert explanation["docking"]["points"] == entry["docking_points"]
        assert explanation["docking"]["component"] == entry["docking_component"]
        assert explanation["properties"]["points"] == entry["property_points"]
        assert explanation["properties"]["component"] == entry["property_component"]
        # The score is exactly the sum of the two contributions shown.
        assert entry["docking_points"] + entry["property_points"] == entry["screening_priority_score"]


def test_property_points_add_up_to_property_contribution() -> None:
    entries, explanations = explain([make_compound("A", -6.0, molecular_weight=120.0, logp=6.0)])
    items = explanations[0]["properties"]["items"]

    assert sum(item["points"] for item in items) == pytest.approx(entries[0]["property_points"], abs=1e-9)
    for item in items:
        assert item["subscore"] == entries[0]["property_subscores"][item["key"]]


def test_maximum_points_follow_configured_weights() -> None:
    explanation = explain([make_compound("A")])[1][0]

    assert explanation["docking"]["max_points"] == pytest.approx(100 * config.RANKING_DOCKING_WEIGHT)
    assert explanation["properties"]["max_points"] == pytest.approx(100 * config.RANKING_PROPERTY_WEIGHT)


def test_docking_text_uses_the_real_numbers() -> None:
    explanation = explain([make_compound("A", -6.1)])[1][0]
    points = explanation["docking"]["points"]

    assert explanation["docking"]["text"] == (
        f"Predicted docking affinity: -6.1 kcal/mol. On Fairouz's fixed docking scale, "
        f"this contributes {points:.1f} points out of a possible 70."
    )
    assert "not experimental binding evidence" in explanation["docking"]["note"]


# --- Descriptor interpretations ----------------------------------------------------

@pytest.mark.parametrize(
    ("value", "position"),
    [(159.9, POSITION_BELOW), (160.0, POSITION_WITHIN), (500.0, POSITION_WITHIN), (500.1, POSITION_ABOVE), (20.0, POSITION_BELOW)],
)
def test_descriptor_position_respects_configured_thresholds(value: float, position: str) -> None:
    explanation = explain([make_compound("A", molecular_weight=value)])[1][0]
    item = explanation["properties"]["items"][0]

    assert item["key"] == "molecular_weight"
    assert item["position"] == position
    assert (item["subscore"] == 1.0) == (position == POSITION_WITHIN)


def test_descriptor_interpretation_wording() -> None:
    items = {item["key"]: item for item in explain([make_compound("A", molecular_weight=78.1)])[1][0]["properties"]["items"]}

    assert items["logp"]["interpretation"] == "LogP: 2.00 — within configured screening range"
    assert items["molecular_weight"]["interpretation"].startswith(
        "Molecular Weight: 78.1 g/mol — below configured screening range, reducing the property score by"
    )


def test_every_descriptor_is_explained_in_config_order() -> None:
    items = explain([make_compound("A")])[1][0]["properties"]["items"]

    assert [item["key"] for item in items] == list(config.PROPERTY_SCORING_RANGES)


def test_negative_zero_is_not_displayed() -> None:
    items = explain([make_compound("A", logp=-0.0014)])[1][0]["properties"]["items"]

    assert items[1]["display_value"] == "0.00"


# --- Strengths, limitations and drivers --------------------------------------------

def test_main_factors_for_weak_docking_and_typical_properties() -> None:
    explanation = explain([make_compound("A", -6.0)])[1][0]

    assert explanation["main_positive_factor"] == (
        "All five molecular-property descriptors fall within the configured screening ranges."
    )
    assert explanation["main_limiting_factor"].startswith("Relatively weak predicted docking affinity")


def test_property_can_be_the_main_limiting_factor() -> None:
    explanation = explain([make_compound("A", -10.0, molecular_weight=650.0)])[1][0]

    assert explanation["main_limiting_factor"].startswith("Molecular Weight above configured range")
    assert explanation["main_positive_factor"].startswith("Its predicted docking affinity")
    assert any("Stronger predicted docking affinity" in strength for strength in explanation["strengths"])


def test_maximum_score_has_no_limiting_factor() -> None:
    explanation = explain([make_compound("A", -12.0)])[1][0]

    assert explanation["limitations"] == []
    assert explanation["main_limiting_factor"].startswith("None")


def test_limitations_are_ordered_by_points_lost() -> None:
    explanation = explain([make_compound("A", -9.5, molecular_weight=650.0, logp=5.5)])[1][0]
    lost = [float(re.search(r"([\d.]+) of \d+ (?:docking )?points not earned", text).group(1)) for text in explanation["limitations"]]

    assert lost == sorted(lost, reverse=True)


def test_factors_are_deterministic() -> None:
    compounds = [make_compound("A", -6.5, logp=5.6), make_compound("B", -9.0, tpsa=150.0)]

    first, second = explain(compounds)[1], explain(copy.deepcopy(compounds))[1]
    for one, two in zip(first, second):
        assert one["main_positive_factor"] == two["main_positive_factor"]
        assert one["main_limiting_factor"] == two["main_limiting_factor"]
        assert one["strengths"] == two["strengths"]
        assert one["limitations"] == two["limitations"]


@pytest.mark.parametrize(
    ("docking_points", "property_points", "driver"),
    [(60.0, 10.0, DRIVER_DOCKING), (5.0, 30.0, DRIVER_PROPERTIES), (25.0, 30.0, DRIVER_BOTH), (0.0, 0.0, DRIVER_NONE)],
)
def test_primary_driver_from_contribution_shares(docking_points: float, property_points: float, driver: str) -> None:
    assert explainability.primary_driver(docking_points, property_points)[0] == driver


# --- Comparisons and ties ----------------------------------------------------------

def test_comparison_when_docking_decides() -> None:
    explanations = by_id(explain([make_compound("A", -8.0), make_compound("B", -6.0)])[1])

    assert explanations["A"]["comparison_below"].startswith(
        "Compound A ranks above Compound B mainly because its predicted docking affinity is more favourable"
    )
    assert explanations["B"]["comparison_above"] == explanations["A"]["comparison_below"]


def test_comparison_when_properties_decide() -> None:
    explanations = by_id(explain([make_compound("A", -7.0), make_compound("B", -7.05, molecular_weight=650.0)])[1])

    text = explanations["A"]["comparison_below"]
    assert text.startswith("Both compounds have similar docking contributions")
    assert "stronger molecular-property profile" in text


def test_comparison_when_both_docking_contributions_are_zero() -> None:
    explanations = by_id(explain([make_compound("A", -2.0), make_compound("B", -3.5, molecular_weight=60.0)])[1])

    assert explanations["A"]["comparison_below"].startswith("Both compounds earn 0 docking points")


def test_tie_with_different_affinities_explains_affinity_tiebreak() -> None:
    # Both weaker than the scale: identical scores, different affinities.
    explanations = by_id(explain([make_compound("A", -2.0), make_compound("B", -3.0)])[1])

    assert explanations["B"]["rank"] == 1
    note = explanations["B"]["tie_notes"][0]
    assert "same Screening Priority Score" in note
    assert "docking affinity is more negative" in note
    assert "compound ID is used only as the final deterministic tie-break" in note
    assert explanations["A"]["tie_notes"] == [note]


def test_tie_with_identical_affinities_explains_compound_id_tiebreak() -> None:
    explanations = by_id(explain([make_compound("B", -6.0), make_compound("A", -6.0)])[1])

    note = explanations["A"]["tie_notes"][0]
    assert explanations["A"]["rank"] == 1
    assert "compound ID (A before B)" in note


def test_displayed_tie_with_different_unrounded_scores() -> None:
    higher = {"name": "H", "compound_id": "H", "screening_priority_score": 41.04, "docking_affinity": -5.0}
    lower = {"name": "L", "compound_id": "L", "screening_priority_score": 41.01, "docking_affinity": -5.0}

    note = explainability.explain_tie(higher, lower)
    assert "same displayed Screening Priority Score (41.0)" in note
    assert "unrounded score is higher" in note


def test_no_tie_note_when_scores_differ() -> None:
    explanations = explain([make_compound("A", -8.0), make_compound("B", -6.0)])[1]

    assert all(explanation["tie_notes"] == [] for explanation in explanations)


# --- Unranked compounds --------------------------------------------------------

def test_failed_docking_gets_unranked_explanation() -> None:
    explanation = explain([make_compound("X", docking_success=False)])[1][0]

    assert explanation["ranked"] is False
    assert explanation["status"] == ranking.STATUS_DOCKING_FAILED
    assert explanation["failed_stage"] == "Molecular docking (AutoDock Vina)"
    assert "none was fabricated" in explanation["summary"]
    assert "Vina exited" in explanation["detail"]
    assert "screening_priority_score" not in explanation


def test_invalid_molecule_gets_appropriate_explanation() -> None:
    explanation = explain([make_compound("X", valid=False, error="RDKit could not parse SMILES 'C1CC('.")])[1][0]

    assert explanation["status"] == ranking.STATUS_INVALID
    assert explanation["failed_stage"] == "SMILES validation (RDKit)"
    assert "could not parse" in explanation["detail"]


@pytest.mark.parametrize(
    ("compound", "status", "stage"),
    [
        (make_compound("X", ligand_prepared=False, docking_success=False), ranking.STATUS_LIGAND_FAILED, "3D ligand preparation"),
        ({**make_compound("X"), "docking": None}, ranking.STATUS_NOT_DOCKED, "did not run"),
        (make_compound("X", logp=None), ranking.STATUS_MISSING_DATA, "Ranking input check"),
        (make_compound("X", affinity=None), ranking.STATUS_MISSING_DATA, "Ranking input check"),
    ],
)
def test_other_unranked_statuses(compound: dict, status: str, stage: str) -> None:
    explanation = explain([compound])[1][0]

    assert explanation["status"] == status
    assert stage in explanation["failed_stage"]
    assert explanation["rank"] is None


# --- Language ----------------------------------------------------------------

def test_explanations_never_claim_efficacy_safety_or_clinical_benefit() -> None:
    compounds = [
        make_compound("A", -9.5),
        make_compound("B", -6.5, logp=6.0),
        make_compound("C", -3.0, tpsa=180.0, molecular_weight=60.0),
        make_compound("D", -3.0),
        make_compound("E", valid=False, error="bad SMILES"),
        make_compound("F", docking_success=False),
    ]
    text = " ".join(all_text(explain(compounds)[1])).lower()

    for pattern in FORBIDDEN_PATTERNS:
        assert not re.search(pattern, text), pattern


# --- Determinism, input safety and edge cases ------------------------------------

def test_repeated_runs_give_identical_explanations() -> None:
    compounds = [make_compound("A", -6.2), make_compound("B", -7.7, logp=5.8), make_compound("C", valid=False)]

    assert explain(compounds)[1] == explain(copy.deepcopy(compounds))[1]


def test_inputs_are_not_mutated() -> None:
    compounds = [make_compound("A", -6.2), make_compound("B", docking_success=False)]
    entries = ranking.rank_compounds(compounds)["ranked_compounds"]
    compounds_before, entries_before = copy.deepcopy(compounds), copy.deepcopy(entries)

    explain_rankings(entries, compounds)

    assert compounds == compounds_before
    assert entries == entries_before


def test_empty_list() -> None:
    assert explain_rankings([], []) == []


def test_single_compound() -> None:
    explanation = explain([make_compound("ONLY", -7.0)])[1][0]

    assert explanation["rank"] == 1
    assert explanation["total_ranked"] == 1
    assert explanation["comparison_above"] is None
    assert explanation["comparison_below"] is None
    assert explanation["summary"].startswith("Compound ONLY ranks #1 of 1")
    assert "most negative predicted docking affinity in the screened set" in explanation["summary"]


def test_order_matches_ranking_entries() -> None:
    entries, explanations = explain([make_compound("BAD", valid=False), make_compound("B", -6.0), make_compound("A", -8.0)])

    assert [explanation["compound_id"] for explanation in explanations] == [entry["compound_id"] for entry in entries]


# --- UI data structures ------------------------------------------------------------

def test_chart_rows_are_valid() -> None:
    explanation = explain([make_compound("A", -6.5, molecular_weight=120.0)])[1][0]
    rows = explanation["chart_rows"]

    assert {(row["component"], row["segment"]) for row in rows} == {
        ("Docking", "Earned"), ("Docking", "Not earned"),
        ("Molecular properties", "Earned"), ("Molecular properties", "Not earned"),
    }
    assert all(row["points"] >= 0 for row in rows)
    totals = {}
    for row in rows:
        totals[row["component"]] = totals.get(row["component"], 0.0) + row["points"]
    assert totals["Docking"] == pytest.approx(ranking.DOCKING_MAX_POINTS)
    assert totals["Molecular properties"] == pytest.approx(ranking.PROPERTY_MAX_POINTS)


def test_ui_charts_and_tables_build_from_explanations() -> None:
    from UI.components import explanation_panel

    explanations = explain([make_compound("A", -8.0), make_compound("B", -5.0), make_compound("C", valid=False)])[1]

    assert explanation_panel.overview_chart(explanations).to_dict()["mark"]["type"] == "bar"
    ranked = [explanation for explanation in explanations if explanation["ranked"]]
    assert explanation_panel._contribution_chart(ranked[0]["chart_rows"]).to_dict()
    rows = explanation_panel._descriptor_rows(ranked[0]["properties"]["items"])
    assert len(rows) == 5 and {"Property", "Value", "Configured range", "Position", "Sub-score", "Points"} <= set(rows[0])


# --- Pipeline integration ------------------------------------------------------------

def test_pipeline_returns_explanations(monkeypatch: pytest.MonkeyPatch) -> None:
    def stub_docking(compounds, progress=None):
        # Test double: distinct stub affinities, never shown in the app.
        return {
            "stage": "docking",
            "status": config.STATUS_COMPLETED,
            "results": [
                {
                    "compound_id": compound["compound_id"],
                    "name": compound["name"],
                    "ligand_prepared": True,
                    "docking_success": True,
                    "best_affinity_kcal_mol": -9.0 + 0.5 * index,
                    "poses_generated": 9,
                    "pose_file": None,
                    "error": None,
                }
                for index, compound in enumerate(compounds)
            ],
            "run_dir": str(Path("stub")),
            "metadata": {},
        }

    monkeypatch.setattr(screening, "run_docking", stub_docking)
    response = screening.run_screening(config.AVAILABLE_TARGETS[0], config.MIN_COMPOUNDS)

    entries, explanations = response["ranking"]["ranked_compounds"], response["explanations"]
    assert [e["compound_id"] for e in explanations] == [e["compound_id"] for e in entries]
    for entry, explanation in zip(entries, explanations):
        assert explanation["ranked"] is (entry["rank"] is not None)
        if explanation["ranked"]:
            assert explanation["screening_priority_score"] == entry["screening_priority_score"]


def test_other_strengths_and_limitations_exclude_the_main_factor() -> None:
    # Weak docking with typical properties: the property profile is the main positive factor.
    weak = explain([make_compound("A", -6.0)])[1][0]
    assert weak["other_strengths"] == []
    assert weak["other_limitations"] == weak["limitations"][1:]

    # Strong docking with one property out of range: docking is the main positive factor.
    strong = explain([make_compound("B", -9.5, molecular_weight=650.0)])[1][0]
    assert strong["main_positive_factor"].startswith("Its predicted docking affinity")
    assert all(not text.startswith("Stronger predicted docking") for text in strong["other_strengths"])
    assert strong["other_strengths"] and strong["other_strengths"][0].startswith("4 of 5 descriptors")
