"""Tests for the descriptor placeholder.

These check the module's contract while RDKit is not integrated. They
deliberately do not test any scientific values.
"""

import config
from backend import descriptors

SAMPLE_SMILES = ["CCO"]  # Ethanol, used only as an input string.


def test_module_exposes_calculate_descriptors() -> None:
    assert callable(descriptors.calculate_descriptors)


def test_placeholder_reports_not_implemented() -> None:
    result = descriptors.calculate_descriptors(SAMPLE_SMILES)

    assert result["stage"] == "descriptors"
    assert result["status"] == config.STATUS_NOT_IMPLEMENTED
    assert result["message"]


def test_placeholder_returns_no_descriptor_values() -> None:
    result = descriptors.calculate_descriptors(SAMPLE_SMILES)

    assert result["descriptors"] == []


def test_placeholder_handles_empty_input() -> None:
    result = descriptors.calculate_descriptors([])

    assert result["status"] == config.STATUS_NOT_IMPLEMENTED
    assert result["descriptors"] == []


def test_placeholder_is_deterministic() -> None:
    first = descriptors.calculate_descriptors(SAMPLE_SMILES)
    second = descriptors.calculate_descriptors(SAMPLE_SMILES)

    assert first == second
