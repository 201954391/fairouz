"""Tests for UI display helpers: consistent precision and no "-0.00"."""

import math

import pytest

import config
from UI.utils import helpers


@pytest.mark.parametrize(
    ("value", "decimals", "expected"),
    [(-0.0014, 2, 0.0), (-0.0004, 1, 0.0), (1.2345, 2, 1.23), (206.285, 1, 206.3), (2.6, 0, 3), (None, 1, None), (math.nan, 1, None)],
)
def test_round_for_display(value, decimals, expected) -> None:
    result = helpers.round_for_display(value, decimals)

    assert result == expected
    if result == 0.0:
        assert math.copysign(1.0, result) == 1.0  # Never negative zero.


@pytest.mark.parametrize("value", [-0.0014, -0.004, -0.0])
def test_format_value_never_shows_negative_zero(value: float) -> None:
    assert helpers.format_value(value, decimals=2) == "0.00"


def test_affinity_and_score_use_configured_precision() -> None:
    assert config.AFFINITY_DECIMALS == 1 and config.SCORE_DECIMALS == 1
    assert helpers.format_affinity(-6.073) == "-6.1 kcal/mol"
    assert helpers.format_priority(54.185) == "54.2 / 100"
    assert helpers.format_affinity(None) == helpers.EMPTY_VALUE
    assert helpers.format_priority(None) == helpers.EMPTY_VALUE


def test_descriptors_use_configured_precision_and_units() -> None:
    assert helpers.format_descriptor("molecular_weight", 206.285) == "206.3 g/mol"
    assert helpers.format_descriptor("logp", -0.0014) == "0.00"
    assert helpers.format_descriptor("tpsa", 37.3) == "37.3 Å²"
    assert helpers.format_descriptor("h_bond_donors", 1) == "1"
    assert helpers.descriptor_number_format("h_bond_acceptors") == "%d"
    assert helpers.descriptor_number_format("logp") == "%.2f"


@pytest.mark.parametrize(("seconds", "expected"), [(0.01, "< 0.1 s"), (2.44, "2.4 s"), (None, helpers.EMPTY_VALUE)])
def test_format_seconds(seconds, expected) -> None:
    assert helpers.format_seconds(seconds) == expected


def test_every_stage_status_has_a_distinct_symbol_and_word() -> None:
    statuses = (config.STATUS_COMPLETED, config.STATUS_PARTIAL, config.STATUS_FAILED, config.STATUS_SKIPPED)
    symbols = {helpers.STAGE_STATUS_DISPLAY[status][0] for status in statuses}

    assert len(symbols) == len(statuses)  # Status is never conveyed by colour alone.
    assert all(helpers.STAGE_STATUS_DISPLAY[status][1] for status in statuses)


def test_pipeline_status_html_shows_words_and_escapes_text() -> None:
    pipeline = [
        {"label": "Docking", "status": config.STATUS_FAILED, "detail": "Vina <missing> & \"quoted\""},
        {"label": "Ranking", "status": config.STATUS_SKIPPED, "detail": "No docking results."},
        {"label": "Report generation", "status": config.STATUS_COMPLETED, "detail": "Written."},
    ]

    markup = helpers.pipeline_status_html(pipeline)

    assert "✕ Docking (Failed)" in markup
    assert "⤼ Ranking (Skipped)" in markup
    assert "✓ Report generation<" in markup  # Completed needs no extra word.
    assert "<missing>" not in markup and "&lt;missing&gt;" in markup


def test_format_status() -> None:
    assert helpers.format_status("not_generated") == "Not generated"
