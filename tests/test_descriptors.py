"""Tests for RDKit descriptor calculation.

Expected values for ethanol and benzene are standard reference values:
average molecular weight from atomic masses, TPSA from Ertl's fragment
contributions (one OH = 20.23 Å²), and the Wildman–Crippen LogP estimate.
Floats are compared with a tolerance, not exactly.
"""

import copy

import pytest

import config
from backend import descriptors

ETHANOL = {"compound_id": "T-1", "name": "Ethanol", "smiles": "CCO", "source": "test"}
BENZENE = {"compound_id": "T-2", "name": "Benzene", "smiles": "c1ccccc1", "source": "test"}
INVALID = {"compound_id": "T-3", "name": "Broken", "smiles": "C1CC(", "source": "test"}


@pytest.fixture
def ethanol_result() -> dict:
    return descriptors.analyse_compound(ETHANOL)


# --- Parsing ---------------------------------------------------------------

def test_valid_smiles_is_parsed() -> None:
    molecule = descriptors.parse_smiles("CCO")

    assert molecule.GetNumAtoms() == 3  # heavy atoms: C, C, O


@pytest.mark.parametrize("smiles", ["C1CC(", "not a smiles", "C(C)(C)(C)(C)C"])
def test_invalid_smiles_raises(smiles: str) -> None:
    with pytest.raises(descriptors.InvalidSmilesError):
        descriptors.parse_smiles(smiles)


@pytest.mark.parametrize("smiles", ["", "   ", None])
def test_empty_smiles_raises(smiles: str | None) -> None:
    # RDKit turns "" into a zero-atom molecule, so this needs an explicit check.
    with pytest.raises(descriptors.InvalidSmilesError, match="empty"):
        descriptors.parse_smiles(smiles)


def test_invalid_compound_is_marked_not_crashed() -> None:
    result = descriptors.analyse_compound(INVALID)

    assert result["valid"] is False
    assert "could not parse" in result["error"]
    assert all(result[key] is None for key in config.DESCRIPTORS)


# --- Individual descriptors (ethanol) ---------------------------------------

def test_ethanol_is_valid(ethanol_result: dict) -> None:
    assert ethanol_result["valid"] is True
    assert ethanol_result["error"] is None


def test_molecular_weight(ethanol_result: dict) -> None:
    assert ethanol_result["molecular_weight"] == pytest.approx(46.07, abs=0.01)


def test_logp(ethanol_result: dict) -> None:
    assert ethanol_result["logp"] == pytest.approx(0.0, abs=0.01)


def test_tpsa(ethanol_result: dict) -> None:
    assert ethanol_result["tpsa"] == pytest.approx(20.23, abs=0.01)


def test_h_bond_donors(ethanol_result: dict) -> None:
    assert ethanol_result["h_bond_donors"] == 1


def test_h_bond_acceptors(ethanol_result: dict) -> None:
    assert ethanol_result["h_bond_acceptors"] == 1


def test_molecule_without_polar_groups() -> None:
    result = descriptors.analyse_compound(BENZENE)

    assert result["molecular_weight"] == pytest.approx(78.11, abs=0.01)
    assert result["tpsa"] == pytest.approx(0.0, abs=0.01)
    assert result["h_bond_donors"] == 0
    assert result["h_bond_acceptors"] == 0


# --- Batches ------------------------------------------------------------------

def test_multiple_compounds_keep_order_and_continue_after_invalid() -> None:
    results = descriptors.calculate_descriptors([ETHANOL, INVALID, BENZENE])

    assert [result["compound_id"] for result in results] == ["T-1", "T-3", "T-2"]
    assert [result["valid"] for result in results] == [True, False, True]


def test_empty_input_returns_empty_list() -> None:
    assert descriptors.calculate_descriptors([]) == []


def test_input_is_not_modified() -> None:
    compounds = [ETHANOL, INVALID]
    original = copy.deepcopy(compounds)

    descriptors.calculate_descriptors(compounds)

    assert compounds == original


def test_result_carries_compound_identity(ethanol_result: dict) -> None:
    for column in config.REQUIRED_COMPOUND_COLUMNS:
        assert ethanol_result[column] == ETHANOL[column]


# --- Consistency with config --------------------------------------------------

def test_calculators_match_configured_descriptors() -> None:
    assert set(descriptors.DESCRIPTOR_CALCULATORS) == set(config.DESCRIPTORS)
