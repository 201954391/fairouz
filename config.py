"""Central configuration for Fairouz.

Every path is built with ``pathlib`` relative to the project root, so the
project behaves the same on Windows, macOS and Linux. Import this module
instead of hardcoding names, paths or limits elsewhere.
"""

import logging
from pathlib import Path

# --- Project identity --------------------------------------------------------
PROJECT_NAME = "Fairouz"
PROJECT_SUBTITLE = "AI-Assisted Virtual Drug Screening"
PROJECT_DESCRIPTION = (
    "Fairouz is a research proof-of-concept that explores how computational "
    "screening and artificial intelligence could help researchers prioritise "
    "promising compounds for further laboratory investigation."
)
PROJECT_VERSION = "0.1.0"
DEVELOPMENT_STAGE = "Step 1 · Foundation & Interface"

# Shown in the sidebar and on the home page. Update as modules are built.
PROJECT_STATUS = {
    "User interface": "Ready",
    "Screening pipeline": "Planned",
    "Compound data": "Not loaded",
    "Docking & AI ranking": "Planned",
}

RESEARCH_DISCLAIMER = (
    "Fairouz is a research and educational proof-of-concept. It is not intended "
    "for clinical use, medical diagnosis, treatment decisions, or direct drug "
    "development decisions without expert scientific review and experimental "
    "validation."
)

# --- Paths -------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data"
COMPOUNDS_DIR = DATA_DIR / "compounds"
TARGETS_DIR = DATA_DIR / "targets"
MODELS_DIR = PROJECT_ROOT / "models"
RESULTS_DIR = PROJECT_ROOT / "results"

# --- Screening settings ------------------------------------------------------
AVAILABLE_TARGETS = ("Demo Viral Protein",)
MIN_COMPOUNDS = 10
MAX_COMPOUNDS = 50
DEFAULT_COMPOUNDS = 10

# (stage name, short description) in pipeline order.
PIPELINE_STAGES = (
    ("Target", "Select a prepared viral protein structure."),
    ("Compound Library", "Load a curated set of public compounds."),
    ("Molecular Analysis", "Calculate drug-like molecular properties."),
    ("Docking", "Estimate how each compound may bind the target."),
    ("AI Ranking", "Combine features to prioritise compounds."),
    ("Results", "Review ranked candidates and the report."),
)

# Status returned by any pipeline stage that has not been built yet.
STATUS_NOT_IMPLEMENTED = "not_implemented"

# --- Logging -----------------------------------------------------------------
LOG_LEVEL = logging.INFO
LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"


def configure_logging(level: int = LOG_LEVEL) -> None:
    """Configure application-wide logging.

    Safe to call on every Streamlit rerun: ``basicConfig`` does nothing once
    the root logger already has a handler.
    """
    logging.basicConfig(level=level, format=LOG_FORMAT)
