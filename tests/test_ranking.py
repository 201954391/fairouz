"""Tests for the ranking placeholder.

These check the module's contract while no ML model exists. They
deliberately do not test any scores.
"""

import config
from backend import ranking

SAMPLE_COMPOUNDS = [{"name": "compound-a"}, {"name": "compound-b"}]


def test_module_exposes_rank_compounds() -> None:
    assert callable(ranking.rank_compounds)


def test_placeholder_reports_not_implemented() -> None:
    result = ranking.rank_compounds(SAMPLE_COMPOUNDS)

    assert result["stage"] == "ranking"
    assert result["status"] == config.STATUS_NOT_IMPLEMENTED
    assert result["message"]


def test_placeholder_returns_no_ranked_compounds() -> None:
    result = ranking.rank_compounds(SAMPLE_COMPOUNDS)

    assert result["ranked_compounds"] == []


def test_placeholder_does_not_modify_input() -> None:
    compounds = [dict(compound) for compound in SAMPLE_COMPOUNDS]

    ranking.rank_compounds(compounds)

    assert compounds == SAMPLE_COMPOUNDS


def test_placeholder_is_deterministic() -> None:
    first = ranking.rank_compounds(SAMPLE_COMPOUNDS)
    second = ranking.rank_compounds(SAMPLE_COMPOUNDS)

    assert first == second
