"""Shared test setup."""

from pathlib import Path

import pytest

import config


@pytest.fixture(autouse=True)
def isolated_outputs(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Send reports and docking runs from every test to a temporary folder.

    Keeps tests from writing into the real results/ folder.
    """
    monkeypatch.setattr(config, "REPORTS_DIR", tmp_path / "reports")
    monkeypatch.setattr(config, "DOCKING_RESULTS_DIR", tmp_path / "docking")
