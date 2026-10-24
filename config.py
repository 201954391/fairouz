"""Central configuration for Fairouz.

Every path is built with ``pathlib`` relative to the project root, so the
project behaves the same on Windows, macOS and Linux. Import this module
instead of hardcoding names, paths or limits elsewhere.
"""

import logging
from pathlib import Path

# --- Project identity --------------------------------------------------------
PROJECT_NAME = "Fairouz"
PROJECT_SUBTITLE = "Transparent Computational Virtual Screening"
PROJECT_DESCRIPTION = (
    "Fairouz is a transparent virtual-screening platform that combines RDKit "
    "molecular descriptors, AutoDock Vina docking, deterministic ranking and "
    "explainability to prioritise compounds for further investigation."
)
PROJECT_VERSION = "0.7.0"
DEVELOPMENT_STAGE = "Step 7 · Demo-ready prototype"

# The components behind each result, shown in the sidebar.
PROJECT_COMPONENTS = {
    "Target": "SARS-CoV-2 Mpro · PDB 6LU7",
    "Compounds": "Public demo set (PubChem)",
    "Descriptors": "RDKit",
    "Docking": "AutoDock Vina 1.2",
    "Ranking": "Transparent heuristic",
    "Reports": "Markdown · JSON · HTML",
}

LIMITATIONS_SUMMARY = (
    "One rigid receptor structure (PDB 6LU7); no waters, fixed protonation.",
    "Small demo compound set of well-known molecules.",
    "Docking affinities are computational estimates, not measured binding.",
    "The ranking is a documented heuristic; no labelled data exist to train a model.",
    "Every result needs experimental follow-up.",
)

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
REPORTS_DIR = RESULTS_DIR / "reports"
REPORT_FORMATS = ("md", "json", "html")

# --- Docking target ----------------------------------------------------------
# Provenance is documented in data/targets/README.md.
TARGET_ID = "6LU7"
TARGET_LABEL = "SARS-CoV-2 Main Protease (PDB 6LU7)"
TARGET_SOURCE = "RCSB Protein Data Bank"
TARGET_SOURCE_URL = "https://www.rcsb.org/structure/6LU7"
TARGET_DOWNLOAD_URL = "https://files.rcsb.org/download/6LU7.pdb"
TARGET_DOWNLOADED_ON = "2026-10-24"
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

# (step name, short description) for the workflow shown on the Home page.
PIPELINE_STAGES = (
    ("Target", "Real X-ray structure of the SARS-CoV-2 main protease (PDB 6LU7)."),
    ("Compounds", "Public molecules with PubChem-verified structures."),
    ("Molecular Analysis", "RDKit calculates five molecular descriptors."),
    ("Docking", "AutoDock Vina predicts binding poses and affinities."),
    ("Ranking", "A transparent Screening Priority Score orders compounds."),
    ("Explainability", "Rule-based breakdown of every score."),
    ("Report", "Reproducible Markdown, JSON and HTML summary."),
)

# Pipeline stages, in order, with the labels used everywhere in the UI.
PIPELINE_STAGE_LABELS = {
    "compound_loading": "Compound loading",
    "molecular_analysis": "Molecular analysis",
    "ligand_preparation": "Ligand preparation",
    "docking": "Docking",
    "ranking": "Ranking",
    "explainability": "Explainability",
    "report": "Report generation",
}

# Statuses reported by pipeline stages.
STATUS_COMPLETED = "completed"  # Every item succeeded.
STATUS_PARTIAL = "partial"  # Some items succeeded, some failed.
STATUS_FAILED = "failed"  # The stage ran, or tried to, and produced nothing.
STATUS_SKIPPED = "skipped"  # Not run because an earlier stage left nothing to do.

# Report-generation statuses.
REPORT_GENERATED = "generated"  # Every requested format was written.
REPORT_PARTIAL = "partial"  # Some formats were written, some failed.
REPORT_NOT_GENERATED = "not_generated"  # Nothing was written.

# --- Display precision -------------------------------------------------------
# Used by the UI, explanations and reports so numbers look the same everywhere.
SCORE_DECIMALS = 1
AFFINITY_DECIMALS = 1  # Vina's precision is not meaningful beyond 0.1 kcal/mol.

# --- Compound data -----------------------------------------------------------
REQUIRED_COMPOUND_COLUMNS = ("compound_id", "name", "smiles", "source")

# Descriptor key -> display metadata (label, unit, display decimals, help). Keys must match the calculators in
# backend/descriptors.py (a test checks this). Help text is deliberately
# neutral: a property value says nothing about whether a molecule works.
DESCRIPTORS = {
    "molecular_weight": {
        "label": "Molecular Weight",
        "unit": "g/mol",
        "decimals": 1,
        "help": "Approximate mass of the molecule.",
    },
    "logp": {
        "label": "LogP",
        "unit": "",
        "decimals": 2,
        "help": "Estimate of how strongly the molecule prefers oily environments over water.",
    },
    "tpsa": {
        "label": "TPSA",
        "unit": "Å²",
        "decimals": 1,
        "help": "Topological polar surface area: a measure of the molecule's polar surface.",
    },
    "h_bond_donors": {
        "label": "H-Bond Donors",
        "unit": "",
        "decimals": 0,
        "help": "Count of groups that can donate a hydrogen bond.",
    },
    "h_bond_acceptors": {
        "label": "H-Bond Acceptors",
        "unit": "",
        "decimals": 0,
        "help": "Count of atoms that can accept a hydrogen bond.",
    },
}

# --- Screening priority ranking ----------------------------------------------
# A transparent heuristic, not a trained model. See README › Screening
# Priority Score. Score = 100 × (docking weight × docking component
#                                + property weight × property component).
RANKING_DOCKING_WEIGHT = 0.7
RANKING_PROPERTY_WEIGHT = 0.3

# Fixed reference scale for the docking component (kcal/mol). Affinities at
# or above WEAK score 0, at or below STRONG score 1, linear in between. A
# fixed scale (rather than rescaling each batch) keeps a compound's score
# independent of which other compounds were screened alongside it.
DOCKING_AFFINITY_WEAK = -4.0
DOCKING_AFFINITY_STRONG = -10.0

# Wording thresholds for explanations (kcal/mol), on the same scale.
DOCKING_STRONG_THRESHOLD = -8.0
DOCKING_MODERATE_THRESHOLD = -6.0

# Property scoring: (zero_below, full_from, full_to, zero_above). A property
# inside [full_from, full_to] scores 1, falls linearly to 0 at the zero_*
# limits, and scores 0 beyond them. Full-score ranges follow widely used
# screening rules of thumb (Lipinski's rule of five, Ghose for the lower
# molecular-weight and LogP limits, Veber for TPSA); the fall-off widths are
# a documented judgement call. These ranges describe typical oral drugs. They
# do not show that a molecule is suitable or effective.
PROPERTY_SCORING_RANGES = {
    "molecular_weight": (80.0, 160.0, 500.0, 700.0),
    "logp": (-2.0, -0.4, 5.0, 7.0),
    "tpsa": (0.0, 0.0, 140.0, 200.0),
    "h_bond_donors": (0.0, 0.0, 5.0, 8.0),
    "h_bond_acceptors": (0.0, 0.0, 10.0, 15.0),
}

# --- Explainability ----------------------------------------------------------
# Thresholds used only for wording explanations; they never change a score.
# "Primary driver" is docking (or properties) when it supplies at least this
# share of a compound's points; otherwise it is "Both".
EXPLANATION_DRIVER_SHARE = 0.6
# Contribution differences smaller than this many points count as "similar".
EXPLANATION_SIMILAR_POINTS = 1.0
# Docking is listed as a strength at or above this component (≈ −7 kcal/mol).
EXPLANATION_DOCKING_STRENGTH_COMPONENT = 0.5

# --- Logging -----------------------------------------------------------------
LOG_LEVEL = logging.INFO
LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"


def configure_logging(level: int = LOG_LEVEL) -> None:
    """Configure application-wide logging.

    Safe to call on every Streamlit rerun: ``basicConfig`` does nothing once
    the root logger already has a handler.
    """
    logging.basicConfig(level=level, format=LOG_FORMAT)
