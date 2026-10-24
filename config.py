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
DEVELOPMENT_STAGE = "Step 3 · Molecular Docking"

# Shown in the sidebar and on the home page. Update as modules are built.
PROJECT_STATUS = {
    "Compound data": "Demo set loaded",
    "Molecular descriptors": "Implemented",
    "Target structure": "6LU7 prepared",
    "Docking": "Implemented (Vina)",
    "AI ranking & report": "Planned",
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
COMPOUNDS_FILE = COMPOUNDS_DIR / "compounds.csv"
MODELS_DIR = PROJECT_ROOT / "models"
RESULTS_DIR = PROJECT_ROOT / "results"
DOCKING_RESULTS_DIR = RESULTS_DIR / "docking"

# --- Docking target ----------------------------------------------------------
# Provenance is documented in data/targets/README.md.
TARGET_ID = "6LU7"
TARGET_LABEL = "SARS-CoV-2 Main Protease (PDB 6LU7)"
TARGET_DIR = TARGETS_DIR / TARGET_ID
RAW_TARGET_FILE = TARGET_DIR / "raw" / "6LU7.pdb"
# SHA-256 of the file downloaded from RCSB; detects accidental edits.
RAW_TARGET_SHA256 = "385497151c2c093d764ffbc694589c2c1db35be638e49d715d4d8af38903966e"
PREPARED_TARGET_DIR = TARGET_DIR / "prepared"
CLEAN_RECEPTOR_FILE = PREPARED_TARGET_DIR / "6LU7_receptor_clean.pdb"
PREPARED_RECEPTOR_FILE = PREPARED_TARGET_DIR / "6LU7_receptor.pdbqt"
PREPARATION_LOG_FILE = PREPARED_TARGET_DIR / "preparation_log.json"
RECEPTOR_CHAIN = "A"  # The protease monomer.
REFERENCE_LIGAND_CHAIN = "C"  # The co-crystallised inhibitor N3.

# Binding box, in ångströms, derived from the N3 inhibitor coordinates by
# backend.receptor.calculate_binding_box(): the centre is the midpoint of
# N3's extent on each axis, and each side is N3's extent plus
# BINDING_BOX_PADDING on both ends. A test recomputes these from the raw file.
BINDING_BOX_PADDING = 5.0
BINDING_BOX_CENTER = (-9.7675, 11.436, 68.9035)
BINDING_BOX_SIZE = (19.259, 29.956, 21.343)

# --- Docking settings (AutoDock Vina) ----------------------------------------
VINA_EXECUTABLE_ENV_VAR = "FAIROUZ_VINA"  # Optional explicit path to Vina.
VINA_EXHAUSTIVENESS = 8  # Vina's default search effort.
VINA_NUM_POSES = 9  # Vina's default number of output poses.
VINA_TIMEOUT_SECONDS = 300  # Per ligand.
# Used for RDKit conformer generation and the Vina search, so runs repeat.
RANDOM_SEED = 42

# --- Screening settings ------------------------------------------------------
AVAILABLE_TARGETS = (TARGET_LABEL,)
MIN_COMPOUNDS = 10
MAX_COMPOUNDS = 50
DEFAULT_COMPOUNDS = 10

# (stage name, short description, build status) in pipeline order.
PIPELINE_STAGES = (
    ("Target", "Use a prepared viral protein structure.", "Implemented"),
    ("Compound Library", "Load a curated set of public compounds.", "Implemented"),
    ("Molecular Analysis", "Calculate molecular properties with RDKit.", "Implemented"),
    ("Docking", "Predict binding poses with AutoDock Vina.", "Implemented"),
    ("AI Ranking", "Combine features to prioritise compounds.", "Planned"),
    ("Results", "Review ranked candidates and the report.", "Partial"),
)

# Statuses reported by pipeline stages.
STATUS_COMPLETED = "completed"
STATUS_UNAVAILABLE = "unavailable"  # Implemented, but could not run (e.g. Vina missing).
STATUS_NOT_IMPLEMENTED = "not_implemented"

# --- Compound data -----------------------------------------------------------
REQUIRED_COMPOUND_COLUMNS = ("compound_id", "name", "smiles", "source")

# Descriptor key -> display metadata. Keys must match the calculators in
# backend/descriptors.py (a test checks this). Help text is deliberately
# neutral: a property value says nothing about whether a molecule works.
DESCRIPTORS = {
    "molecular_weight": {
        "label": "Molecular Weight",
        "unit": "g/mol",
        "help": "Approximate mass of the molecule.",
    },
    "logp": {
        "label": "LogP",
        "unit": "",
        "help": "Estimate of how strongly the molecule prefers oily environments over water.",
    },
    "tpsa": {
        "label": "TPSA",
        "unit": "Å²",
        "help": "Topological polar surface area: a measure of the molecule's polar surface.",
    },
    "h_bond_donors": {
        "label": "H-Bond Donors",
        "unit": "",
        "help": "Count of groups that can donate a hydrogen bond.",
    },
    "h_bond_acceptors": {
        "label": "H-Bond Acceptors",
        "unit": "",
        "help": "Count of atoms that can accept a hydrogen bond.",
    },
}

# --- Logging -----------------------------------------------------------------
LOG_LEVEL = logging.INFO
LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"


def configure_logging(level: int = LOG_LEVEL) -> None:
    """Configure application-wide logging.

    Safe to call on every Streamlit rerun: ``basicConfig`` does nothing once
    the root logger already has a handler.
    """
    logging.basicConfig(level=level, format=LOG_FORMAT)
