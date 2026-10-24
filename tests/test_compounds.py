"""Tests for loading the compound dataset."""

from pathlib import Path

import pytest

import config
from backend.compounds import CompoundDataError, load_compounds
from backend.descriptors import parse_smiles

HEADER = "compound_id,name,smiles,source\n"


def write_csv(tmp_path: Path, content: str) -> Path:
    path = tmp_path / "compounds.csv"
    path.write_text(content, encoding="utf-8")
    return path


def test_loads_valid_file(tmp_path: Path) -> None:
    path = write_csv(tmp_path, HEADER + "A,Ethanol,CCO,test\nB,Benzene,c1ccccc1,test\n")

    result = load_compounds(path)

    assert [compound["compound_id"] for compound in result.compounds] == ["A", "B"]
    assert result.issues == []


def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(CompoundDataError, match="not found"):
        load_compounds(tmp_path / "missing.csv")


def test_missing_column_raises(tmp_path: Path) -> None:
    path = write_csv(tmp_path, "compound_id,name,smiles\nA,Ethanol,CCO\n")

    with pytest.raises(CompoundDataError, match="source"):
        load_compounds(path)


def test_malformed_rows_are_skipped_and_reported(tmp_path: Path) -> None:
    content = HEADER + "A,Ethanol,CCO,test\nB,Missing SMILES,,test\nC,Extra,CCO,test,surplus\n"
    path = write_csv(tmp_path, content)

    result = load_compounds(path)

    assert [compound["compound_id"] for compound in result.compounds] == ["A"]
    assert len(result.issues) == 2
    assert "Line 3" in result.issues[0] and "smiles" in result.issues[0]
    assert "Line 4" in result.issues[1]


def test_duplicate_ids_raise(tmp_path: Path) -> None:
    path = write_csv(tmp_path, HEADER + "A,Ethanol,CCO,test\nA,Benzene,c1ccccc1,test\n")

    with pytest.raises(CompoundDataError, match="Duplicate"):
        load_compounds(path)


def test_file_without_usable_rows_raises(tmp_path: Path) -> None:
    path = write_csv(tmp_path, HEADER)

    with pytest.raises(CompoundDataError, match="no usable compounds"):
        load_compounds(path)


def test_values_are_stripped(tmp_path: Path) -> None:
    path = write_csv(tmp_path, HEADER + " A , Ethanol , CCO , test \n")

    compound = load_compounds(path).compounds[0]

    assert compound == {"compound_id": "A", "name": "Ethanol", "smiles": "CCO", "source": "test"}


def test_project_dataset_is_complete_and_parsable() -> None:
    result = load_compounds(config.COMPOUNDS_FILE)

    assert len(result.compounds) == 10
    assert result.issues == []
    for compound in result.compounds:
        assert compound["source"].startswith("PubChem CID ")
        parse_smiles(compound["smiles"])  # raises if any structure is invalid
