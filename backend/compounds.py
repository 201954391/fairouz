"""Compound dataset loading.

Reads the compound library CSV and checks its structure before any science
runs. Kept separate from ``descriptors.py`` because docking and ranking will
load the same library in later stages.

Problems are handled by severity:

* Problems with the whole file (missing file, missing columns, duplicate
  IDs, no usable rows) raise ``CompoundDataError``.
* Problems with a single row (missing values, wrong number of fields) skip
  that row and are reported in ``CompoundLoadResult.issues``, never silently.
"""

import csv
import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import config

logger = logging.getLogger(__name__)


class CompoundDataError(Exception):
    """Raised when the compound dataset cannot be used at all."""


@dataclass(frozen=True)
class CompoundLoadResult:
    """Compounds read from the dataset, plus any rows that were skipped.

    Attributes:
        compounds: One dictionary per usable row, holding the required columns.
        issues: Human-readable descriptions of skipped rows.
    """

    compounds: list[dict[str, str]]
    issues: list[str]


def _check_columns(fieldnames: list[str] | None, path: Path) -> None:
    """Raise if the CSV header lacks any required column."""
    present = {name.strip() for name in fieldnames or []}
    missing = [column for column in config.REQUIRED_COMPOUND_COLUMNS if column not in present]
    if missing:
        raise CompoundDataError(
            f"Compound dataset {path.name} is missing required column(s): {', '.join(missing)}."
        )


def _find_row_problem(row: dict[str | None, str | None]) -> str | None:
    """Describe what is wrong with a CSV row, or return None if it is usable."""
    # DictReader stores surplus fields under the key None.
    if None in row:
        return "row has more fields than the header"
    empty = [
        column
        for column in config.REQUIRED_COMPOUND_COLUMNS
        if not (row.get(column) or "").strip()
    ]
    if empty:
        return f"missing value for {', '.join(empty)}"
    return None


def _check_unique_ids(compounds: list[dict[str, str]]) -> None:
    """Raise if two compounds share an ID, since results are keyed by ID."""
    counts = Counter(compound["compound_id"] for compound in compounds)
    duplicates = sorted(compound_id for compound_id, count in counts.items() if count > 1)
    if duplicates:
        raise CompoundDataError(f"Duplicate compound_id value(s): {', '.join(duplicates)}.")


def load_compounds(path: Path | None = None) -> CompoundLoadResult:
    """Load and check the compound dataset.

    Args:
        path: CSV file to read. Defaults to ``config.COMPOUNDS_FILE``.

    Returns:
        The usable compounds in file order and a list of skipped-row issues.

    Raises:
        CompoundDataError: If the file is missing, unreadable, lacks required
            columns, has duplicate IDs or contains no usable compounds.
    """
    path = path or config.COMPOUNDS_FILE
    if not path.is_file():
        raise CompoundDataError(f"Compound dataset not found: {path}")

    compounds: list[dict[str, str]] = []
    issues: list[str] = []
    try:
        # utf-8-sig also accepts files saved by Excel, which start with a BOM.
        with path.open(newline="", encoding="utf-8-sig") as file:
            reader = csv.DictReader(file)
            reader.fieldnames = [name.strip() for name in reader.fieldnames or []]
            _check_columns(reader.fieldnames, path)
            for row in reader:
                problem = _find_row_problem(row)
                if problem:
                    issue = f"Line {reader.line_num}: {problem}; row skipped."
                    logger.warning("Compound dataset %s: %s", path.name, issue)
                    issues.append(issue)
                    continue
                compounds.append(
                    {column: row[column].strip() for column in config.REQUIRED_COMPOUND_COLUMNS}
                )
    except (UnicodeDecodeError, csv.Error) as error:
        raise CompoundDataError(f"Compound dataset {path.name} could not be read: {error}") from error

    _check_unique_ids(compounds)
    if not compounds:
        raise CompoundDataError(f"Compound dataset {path.name} contains no usable compounds.")

    logger.info(
        "Loaded %d compound(s) from %s (%d row(s) skipped).",
        len(compounds),
        path.name,
        len(issues),
    )
    return CompoundLoadResult(compounds=compounds, issues=issues)
