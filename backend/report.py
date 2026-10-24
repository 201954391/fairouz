"""Deterministic research report generation.

Turns a finished screening response into a reproducible report in Markdown,
JSON and HTML. The report is assembled from fixed templates and the real
pipeline data. There is no language model, and nothing is recalculated:
descriptor values, docking affinities, ranking points and explanations are
read from the response exactly as the pipeline produced them.

Structure:

1. ``build_report_data(response)``: the structured report (pure, no I/O).
   This is also the JSON format.
2. ``build_document(data)``: an ordered list of simple blocks (headings,
   paragraphs, lists, tables).
3. ``render_markdown`` / ``render_html``: render the same blocks, so the two
   human-readable formats always agree.
4. ``generate_report(response)``: writes the files to ``results/reports/``
   and reports whether it succeeded. It never raises, so a failed report
   cannot destroy the screening results.

The report documents computational results only. It states that
experimental validation is required.
"""

import hashlib
import html
import json
import logging
import math
import platform
import re
import subprocess
from collections.abc import Mapping, Sequence
from datetime import date, datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path, PurePath
from typing import Any

import config
from backend import ranking
from backend.receptor import read_structure_provenance

logger = logging.getLogger(__name__)

SCHEMA_VERSION = "1.0"
REPORT_TITLE = "Fairouz Virtual Screening Report"
TOP_COUNT = 3

VALIDATION_STATEMENT = (
    "All results in this report are computational. Experimental validation is "
    "required before drawing any conclusion about binding, biological activity "
    "or suitability for further development."
)

LIMITATIONS = (
    "Docking is a computational method: it predicts likely poses and approximate interaction scores; it does not measure binding.",
    "The receptor is treated as rigid, without waters, with Gasteiger charges and a fixed histidine protonation state (HIE).",
    "One protein structure (PDB 6LU7) is used; it captures a single conformation of the protease.",
    "The predicted AutoDock Vina affinity is a computational estimate, not experimental binding evidence or an experimental binding energy.",
    "The ranking weights (70% docking, 30% properties) and the fixed docking scale are documented heuristic choices, not values learned from data.",
    "The descriptor screening ranges are heuristic rules of thumb for typical oral drugs; falling inside them does not show that a molecule is suitable.",
    "The current compound set is small (a demo set of well-known molecules) and does not represent the chemical space relevant to drug discovery.",
    "No supervised machine-learning model is used, because no appropriate labelled dataset is available.",
    "The explainability layer describes how the scoring formula produced each score; it does not describe biology or mechanism of action.",
    "Experimental validation is required for every compound before any conclusion is drawn.",
)

RESPONSIBLE_USE = (
    "No patient data are used.",
    "No personally identifiable information is processed.",
    "Compound structures come from PubChem and the protein structure from the RCSB Protein Data Bank, both public scientific sources.",
    "Provenance is documented in data/compounds/README.md and data/targets/README.md, and summarised in this report.",
    "These computational predictions must not be used for clinical decision-making.",
    "Fairouz is intended for research prioritisation only: it suggests which compounds may merit further investigation.",
)

DESCRIPTOR_KEYS = tuple(config.DESCRIPTORS)


# --- Serialisation ---------------------------------------------------------------

def to_jsonable(value: Any) -> Any:
    """Convert a value into plain JSON types.

    Paths become strings, datetimes ISO strings, NumPy scalars and arrays
    Python values, and NaN or infinity ``None`` (JSON has no NaN).
    """
    if value is None or isinstance(value, (bool, str)):
        return value
    if type(value).__module__ == "numpy":
        return to_jsonable(value.tolist() if hasattr(value, "tolist") else value.item())
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, PurePath):
        return value.as_posix()
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if isinstance(value, Mapping):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset)):
        return [to_jsonable(item) for item in sorted(value, key=str)]
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    return str(value)


# --- Metadata helpers ------------------------------------------------------------

def relative_path(path: str | Path | None) -> str | None:
    """A path relative to the project root when possible, using forward slashes."""
    if path is None:
        return None
    path = Path(path)
    try:
        return path.resolve().relative_to(config.PROJECT_ROOT).as_posix()
    except (ValueError, OSError):
        return path.as_posix()


def file_sha256(path: Path) -> str | None:
    """SHA-256 of a file, or None if it cannot be read."""
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def git_revision() -> dict[str, Any] | None:
    """The current Git commit and whether the working tree has changes, if available."""
    try:
        commit = subprocess.run(
            ["git", "-C", str(config.PROJECT_ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        status = subprocess.run(
            ["git", "-C", str(config.PROJECT_ROOT), "status", "--porcelain"],
            capture_output=True, text=True, timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if commit.returncode != 0 or status.returncode != 0:
        return None
    return {"commit": commit.stdout.strip(), "uncommitted_changes": bool(status.stdout.strip())}


def _package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def software_versions(docking_metadata: Mapping[str, Any] | None) -> dict[str, str | None]:
    """Versions of the main tools; any that cannot be detected are None."""
    vina = ((docking_metadata or {}).get("vina") or {}).get("version")
    return {
        "python": platform.python_version(),
        "rdkit": _package_version("rdkit"),
        "meeko": _package_version("meeko"),
        "gemmi": _package_version("gemmi"),
        "streamlit": _package_version("streamlit"),
        "autodock_vina": vina,
    }


def configuration_snapshot() -> dict[str, Any]:
    """Every configured input that affects results, for hashing and the report."""
    return to_jsonable(
        {
            "target_id": config.TARGET_ID,
            "receptor_chain": config.RECEPTOR_CHAIN,
            "binding_box": {"center": config.BINDING_BOX_CENTER, "size": config.BINDING_BOX_SIZE, "padding": config.BINDING_BOX_PADDING},
            "vina": {"exhaustiveness": config.VINA_EXHAUSTIVENESS, "num_modes": config.VINA_NUM_POSES},
            "random_seed": config.RANDOM_SEED,
            "ranking_weights": {"docking": config.RANKING_DOCKING_WEIGHT, "properties": config.RANKING_PROPERTY_WEIGHT},
            "docking_scale_kcal_mol": {"zero_at": config.DOCKING_AFFINITY_WEAK, "full_at": config.DOCKING_AFFINITY_STRONG},
            "property_ranges": config.PROPERTY_SCORING_RANGES,
            "compound_dataset_sha256": file_sha256(config.COMPOUNDS_FILE),
            "raw_target_sha256": config.RAW_TARGET_SHA256,
        }
    )


def configuration_hash(snapshot: Mapping[str, Any]) -> str:
    """Short, deterministic hash of the configuration snapshot."""
    canonical = json.dumps(snapshot, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:10]


def results_hash(response: Mapping[str, Any]) -> str:
    """Short, deterministic hash of what the run produced.

    Two runs in the same second with the same configuration but different
    outcomes (e.g. one without Vina) therefore get different run IDs.
    """
    outcome = [
        [
            compound.get("compound_id"),
            bool(compound.get("valid")),
            _docking_status(compound.get("docking")),
            (compound.get("docking") or {}).get("best_affinity_kcal_mol"),
            (compound.get("ranking") or {}).get("screening_priority_score"),
        ]
        for compound in response.get("compounds") or []
    ]
    canonical = json.dumps(to_jsonable(outcome), separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:6]


def make_run_id(now: datetime, config_hash: str, outcome_hash: str) -> str:
    """Traceable run ID: UTC timestamp, target, configuration hash and results hash.

    No randomness: the same inputs and results always give the same ID.
    """
    return f"{now:%Y%m%d-%H%M%S}_{config.TARGET_ID}_{config_hash}-{outcome_hash}"


def safe_filename_part(text: str) -> str:
    """Keep only characters that are safe in filenames on Windows, macOS and Linux."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("._")
    return cleaned or "report"


# --- Building the structured report ----------------------------------------------

def _read_preparation_log() -> dict[str, Any] | None:
    try:
        return json.loads(config.PREPARATION_LOG_FILE.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _target_section() -> dict[str, Any]:
    try:
        provenance = read_structure_provenance(config.RAW_TARGET_FILE)
    except OSError:
        provenance = {}
    preparation = _read_preparation_log() or {}
    return {
        "target_id": config.TARGET_ID,
        "name": config.TARGET_LABEL,
        "source": config.TARGET_SOURCE,
        "source_url": config.TARGET_SOURCE_URL,
        "download_url": config.TARGET_DOWNLOAD_URL,
        "downloaded_on": config.TARGET_DOWNLOADED_ON,
        "raw_file": relative_path(config.RAW_TARGET_FILE),
        "structure": provenance,
        "receptor_chain": config.RECEPTOR_CHAIN,
        "reference_ligand_chain": config.REFERENCE_LIGAND_CHAIN,
        "preparation": {
            "steps": preparation.get("steps", []),
            "clean_receptor": preparation.get("clean_receptor"),
            "tool_versions": preparation.get("tool_versions"),
            "prepared_at": preparation.get("prepared_at"),
            "prepared_receptor_file": relative_path(config.PREPARED_RECEPTOR_FILE),
        },
        "binding_box": {
            "center": config.BINDING_BOX_CENTER,
            "size": config.BINDING_BOX_SIZE,
            "padding": config.BINDING_BOX_PADDING,
            "units": "angstrom",
            "definition": (
                "Centre: midpoint of the co-crystallised inhibitor N3's extent on each axis. "
                f"Size: N3's extent plus {config.BINDING_BOX_PADDING:g} Å padding on each side."
            ),
        },
    }


def _format_triplet(values: Sequence[float]) -> str:
    return "(" + ", ".join(f"{value:g}" for value in values) + ")"


def _has_lower_limit(scoring_range: Sequence[float]) -> bool:
    """False for ranges like TPSA's, where zero_below == full_from (no lower penalty)."""
    return scoring_range[0] != scoring_range[1]


def _range_sentence(key: str, scoring_range: Sequence[float]) -> str:
    zero_below, full_from, full_to, zero_above = scoring_range
    label = config.DESCRIPTORS[key]["label"]
    if _has_lower_limit(scoring_range):
        return f"{label} full score {full_from:g} to {full_to:g}, zero at ≤{zero_below:g} or ≥{zero_above:g}"
    return f"{label} full score up to {full_to:g}, zero at ≥{zero_above:g}"


def _methodology(summary: Mapping[str, Any], docking_info: Mapping[str, Any], target: Mapping[str, Any]) -> list[dict[str, str]]:
    """Methods, worded from what this run actually did."""
    metadata = docking_info.get("metadata") or {}
    vina = metadata.get("vina") or {}
    ligand_prep = metadata.get("ligand_preparation") or {}
    ranges = "; ".join(_range_sentence(key, scoring_range) for key, scoring_range in config.PROPERTY_SCORING_RANGES.items())
    if docking_info.get("error"):
        docking_text = (
            "AutoDock Vina docking is configured but did not run in this session, so no "
            f"docking results exist: {docking_info['error']}"
        )
    else:
        docking_text = (
            f"{vina.get('version', 'AutoDock Vina')} docked each prepared ligand into the rigid receptor "
            f"within a box centred at {_format_triplet(config.BINDING_BOX_CENTER)} Å with edges "
            f"{_format_triplet(config.BINDING_BOX_SIZE)} Å (exhaustiveness {vina.get('exhaustiveness', config.VINA_EXHAUSTIVENESS)}, "
            f"{vina.get('num_modes', config.VINA_NUM_POSES)} poses, seed {vina.get('seed', config.RANDOM_SEED)}). "
            f"{target['binding_box']['definition']}"
        )
    preparation = target["preparation"]
    counts = preparation.get("clean_receptor") or {}
    if preparation["steps"]:
        receptor_text = (
            f"The raw structure was reduced to protease chain {config.RECEPTOR_CHAIN}"
            + (f" ({counts['residues']} residues)" if counts.get("residues") else "")
            + "; the co-crystallised inhibitor N3 and "
            + (f"all {counts['waters_removed']} " if counts.get("waters_removed") is not None else "all ")
            + "crystallographic waters were removed; and Meeko added polar hydrogens, Gasteiger charges and "
            f"AutoDock atom types to produce a rigid PDBQT receptor. The {len(preparation['steps'])} recorded "
            "preparation steps are listed under Target."
        )
    else:
        receptor_text = "The receptor preparation log is not available."
    return [
        {
            "title": "Compound source and provenance",
            "text": (
                f"{summary.get('loaded', 0)} compounds were loaded from {relative_path(config.COMPOUNDS_FILE)}, a demo set "
                "of well-known small molecules whose SMILES strings were retrieved from PubChem; each row records "
                "its PubChem compound ID. Inclusion in the set does not imply any biological activity."
            ),
        },
        {
            "title": "SMILES validation",
            "text": (
                "Each SMILES string was parsed with RDKit. Empty or unparsable strings were marked invalid "
                "with a recorded reason and excluded from later stages."
            ),
        },
        {
            "title": "Molecular descriptors",
            "text": (
                "RDKit calculated average molecular weight, Wildman–Crippen LogP, Ertl topological polar "
                "surface area (N and O), and Lipinski-style hydrogen-bond donor and acceptor counts."
            ),
        },
        {
            "title": "Ligand 3D preparation",
            "text": (
                f"Hydrogens were added and one 3D conformer was embedded with {ligand_prep.get('conformer', 'RDKit ETKDGv3')} "
                f"(seed {ligand_prep.get('seed', config.RANDOM_SEED)}), optimised with "
                f"{ligand_prep.get('force_field', 'MMFF94 (UFF fallback)')}, and converted to PDBQT with Meeko. "
                "Molecules were used in the protonation state written in their SMILES."
            ),
        },
        {"title": "Receptor preparation", "text": receptor_text},
        {"title": "Molecular docking", "text": docking_text},
        {
            "title": "Docking affinity interpretation",
            "text": (
                "The reported value is the predicted AutoDock Vina affinity (kcal/mol) of each compound's best pose. "
                "More negative values indicate a more favourable predicted pose within this docking setup. "
                "It is a computational estimate, not an experimental binding energy."
            ),
        },
        {
            "title": "Ranking (Screening Priority Score)",
            "text": (
                "Screening Priority Score = docking_points + property_points. Docking contributes up to "
                f"{ranking.DOCKING_MAX_POINTS:g} points on a fixed scale ({config.DOCKING_AFFINITY_WEAK:g} kcal/mol or weaker "
                f"scores 0; {config.DOCKING_AFFINITY_STRONG:g} kcal/mol or stronger scores the maximum; linear in between). "
                f"Molecular properties contribute up to {ranking.PROPERTY_MAX_POINTS:g} points: the average of five "
                f"descriptor sub-scores ({ranges}). Ties are ordered by more negative docking affinity, then compound ID. "
                "The score is a heuristic for research prioritisation, not a trained model."
            ),
        },
        {
            "title": "Explainability",
            "text": (
                "Each explanation is generated by fixed rules from the ranking's own numbers (docking points, property "
                "points and sub-scores). No language model and no machine-learning explainability method is used, "
                "because no machine-learning model exists."
            ),
        },
    ]


def _docking_status(docking: Mapping[str, Any] | None) -> str:
    if docking is None:
        return "Not docked"
    if docking.get("docking_success"):
        return "Docked"
    if not docking.get("ligand_prepared"):
        return "Ligand preparation failed"
    return "Docking failed"


def _compound_record(compound: Mapping[str, Any]) -> dict[str, Any]:
    docking = compound.get("docking")
    entry = compound.get("ranking") or {}
    return {
        "compound_id": compound.get("compound_id"),
        "name": compound.get("name"),
        "smiles": compound.get("smiles"),
        "source": compound.get("source"),
        "valid": bool(compound.get("valid")),
        "descriptor_error": compound.get("error"),
        "descriptors": {key: compound.get(key) for key in DESCRIPTOR_KEYS},
        "docking": {
            "status": _docking_status(docking),
            "predicted_affinity_kcal_mol": (docking or {}).get("best_affinity_kcal_mol"),
            "poses_generated": (docking or {}).get("poses_generated"),
            "pose_file": relative_path((docking or {}).get("pose_file")),
            "error": (docking or {}).get("error"),
        },
        "ranking": {
            "rank": entry.get("rank"),
            "status": entry.get("ranking_status"),
            "screening_priority_score": entry.get("screening_priority_score"),
        },
    }


def _ranking_section(response: Mapping[str, Any], explanations_by_id: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    entries = (response.get("ranking") or {}).get("ranked_compounds") or []
    ranked = [
        {
            "rank": entry["rank"],
            "compound_id": entry["compound_id"],
            "name": entry["name"],
            "screening_priority_score": entry["screening_priority_score"],
            "predicted_affinity_kcal_mol": entry["docking_affinity"],
            "docking_points": entry["docking_points"],
            "property_points": entry["property_points"],
            "docking_component": entry["docking_component"],
            "property_component": entry["property_component"],
            "primary_driver": (explanations_by_id.get(entry["compound_id"]) or {}).get("primary_driver"),
            "explanation": entry["explanation"],
        }
        for entry in entries
        if entry.get("rank") is not None
    ]
    ranked.sort(key=lambda item: item["rank"])
    unranked = [
        {
            "compound_id": entry["compound_id"],
            "name": entry["name"],
            "status": entry["ranking_status"],
            "failed_stage": (explanations_by_id.get(entry["compound_id"]) or {}).get("failed_stage"),
            "detail": (explanations_by_id.get(entry["compound_id"]) or {}).get("detail"),
        }
        for entry in entries
        if entry.get("rank") is None
    ]
    return {
        "formula": "Screening Priority Score = docking_points + property_points",
        "docking_max_points": ranking.DOCKING_MAX_POINTS,
        "property_max_points": ranking.PROPERTY_MAX_POINTS,
        "method": (response.get("ranking") or {}).get("method"),
        "top": ranked[:TOP_COUNT],
        "ranked": ranked,
        "unranked": unranked,
    }


def _explanation_record(explanation: Mapping[str, Any]) -> dict[str, Any]:
    if not explanation.get("ranked"):
        return {
            "compound_id": explanation["compound_id"],
            "name": explanation["name"],
            "ranked": False,
            "status": explanation["status"],
            "failed_stage": explanation["failed_stage"],
            "detail": explanation["detail"],
            "summary": explanation["summary"],
        }
    return {
        "compound_id": explanation["compound_id"],
        "name": explanation["name"],
        "ranked": True,
        "rank": explanation["rank"],
        "docking_points": explanation["docking"]["points"],
        "property_points": explanation["properties"]["points"],
        "primary_driver": explanation["primary_driver"],
        "main_positive_factor": explanation["main_positive_factor"],
        "main_limiting_factor": explanation["main_limiting_factor"],
        "comparison_above": explanation["comparison_above"],
        "comparison_below": explanation["comparison_below"],
        "tie_notes": list(explanation["tie_notes"]),
        "summary": explanation["summary"],
        "properties": [
            {
                "label": item["label"],
                "value": item["value"],
                "display_value": item["display_value"],
                "configured_range": item["range_text"],
                "position": item["position"],
                "subscore": item["subscore"],
                "points": item["points"],
            }
            for item in explanation["properties"]["items"]
        ],
    }


def _outcome(summary: Mapping[str, Any], docking_info: Mapping[str, Any]) -> tuple[str, str]:
    """A one-line description of what happened, and the run's completeness."""
    processed, ranked_count = summary.get("processed", 0), summary.get("ranked", 0)
    if processed == 0:
        return "No compounds were processed, so there are no results to report.", "empty"
    if docking_info.get("error"):
        return (
            f"Molecular descriptors were calculated for {summary.get('valid', 0)} of {processed} compounds, but docking "
            "could not run, so no compound was docked or ranked.",
            "partial",
        )
    if ranked_count == 0:
        return f"{processed} compounds were processed, but none was successfully docked, so none was ranked.", "partial"
    completeness = "complete" if ranked_count == processed else "partial"
    return f"{ranked_count} of {processed} processed compounds were docked and ranked.", completeness


def build_report_data(response: Mapping[str, Any], now: datetime | None = None) -> dict[str, Any]:
    """Assemble the structured report from a screening response.

    Args:
        response: The output of ``backend.screening.run_screening`` (not modified).
        now: Report timestamp (UTC). Defaults to the current time.

    Returns:
        A JSON-serialisable dictionary; also the content of the JSON report.
    """
    now = now or datetime.now(timezone.utc)
    summary_in = response.get("summary") or {}
    docking_info = response.get("docking") or {}
    compounds = response.get("compounds") or []
    explanations = response.get("explanations") or []
    explanations_by_id = {explanation["compound_id"]: explanation for explanation in explanations}

    snapshot = configuration_snapshot()
    config_hash = configuration_hash(snapshot)
    run_id = make_run_id(now, config_hash, results_hash(response))
    target = _target_section()
    outcome, completeness = _outcome(summary_in, docking_info)
    processed = summary_in.get("processed", len(compounds))
    ranked_count = summary_in.get("ranked", 0)

    data = {
        "schema_version": SCHEMA_VERSION,
        "metadata": {
            "title": REPORT_TITLE,
            "run_id": run_id,
            "generated_at": now,
            "project": config.PROJECT_NAME,
            "project_version": config.PROJECT_VERSION,
            "development_stage": config.DEVELOPMENT_STAGE,
            "git": git_revision(),
            # Wall-clock time per stage; varies between runs and machines.
            "timings_seconds": dict(response.get("timings") or {}),
        },
        "target": target,
        "methodology": _methodology(summary_in, docking_info, target),
        "summary": {
            "target": config.TARGET_LABEL,
            "pdb_id": config.TARGET_ID,
            "requested": summary_in.get("requested"),
            "loaded": summary_in.get("loaded"),
            "processed": processed,
            "successfully_analysed": summary_in.get("valid", 0),
            "invalid": summary_in.get("invalid", 0),
            "ligands_prepared": summary_in.get("ligands_prepared", 0),
            "docked": summary_in.get("docked", 0),
            "docking_failed": summary_in.get("docking_failed", 0),
            "ranked": ranked_count,
            "unranked": processed - ranked_count,
            "stages": response.get("stages") or {},
            "warnings": list(response.get("warnings") or []),
            "docking_error": docking_info.get("error"),
            "outcome": outcome,
            "completeness": completeness,
            "validation_statement": VALIDATION_STATEMENT,
            "pipeline": [
                {"stage": stage["label"], "status": stage["status"], "detail": stage["detail"]}
                for stage in response.get("pipeline") or []
            ],
        },
        "compounds": [_compound_record(compound) for compound in compounds],
        "ranking": _ranking_section(response, explanations_by_id),
        "explanations": [_explanation_record(explanation) for explanation in explanations],
        "limitations": list(LIMITATIONS),
        "responsible_use": list(RESPONSIBLE_USE),
        "reproducibility": {
            "run_id": run_id,
            "configuration_hash": config_hash,
            "configuration": snapshot,
            "software_versions": software_versions(docking_info.get("metadata")),
            "docking_run_folder": relative_path(docking_info.get("run_dir")),
            "integrity": {
                "compound_dataset": {"path": relative_path(config.COMPOUNDS_FILE), "sha256": snapshot["compound_dataset_sha256"]},
                "raw_target": {
                    "path": relative_path(config.RAW_TARGET_FILE),
                    "sha256": file_sha256(config.RAW_TARGET_FILE),
                    "expected_sha256": config.RAW_TARGET_SHA256,
                },
                "prepared_receptor": {
                    "path": relative_path(config.PREPARED_RECEPTOR_FILE),
                    "sha256": file_sha256(config.PREPARED_RECEPTOR_FILE),
                    "note": (
                        "Prepared files are written with LF line endings, matching how Git stores them, so this "
                        "checksum matches a fresh checkout. A copy with converted line endings would differ."
                    ),
                },
            },
        },
    }
    return to_jsonable(data)


# --- Document blocks -------------------------------------------------------------

def _num(value: Any, decimals: int) -> str:
    if value is None:
        return "—"
    return f"{round(value, decimals) + 0.0:.{decimals}f}"  # + 0.0 avoids "-0.00".


def _summary_blocks(data: Mapping[str, Any]) -> list[tuple]:
    meta, summary = data["metadata"], data["summary"]
    git = meta.get("git")
    git_text = "not available" if not git else f"{git['commit']}{' (with uncommitted changes)' if git['uncommitted_changes'] else ''}"
    rows = [
        ["Run ID", meta["run_id"]],
        ["Generated (UTC)", meta["generated_at"]],
        ["Fairouz version", f"{meta['project_version']} ({meta['development_stage']})"],
        ["Git commit", git_text],
        ["Target", summary["target"]],
        ["Compounds requested", summary["requested"]],
        ["Compounds processed", summary["processed"]],
        ["Successfully analysed (descriptors)", summary["successfully_analysed"]],
        ["Docked", summary["docked"]],
        ["Ranked", summary["ranked"]],
        ["Unranked", summary["unranked"]],
        ["Run completeness", summary["completeness"]],
    ]
    blocks = [("h", 2, "Run Summary"), ("p", summary["outcome"]), ("table", ["Field", "Value"], rows)]
    if summary["pipeline"]:
        blocks += [
            ("h", 3, "Pipeline stages"),
            (
                "table",
                ["Stage", "Status", "Detail"],
                [[stage["stage"], stage["status"].capitalize(), stage["detail"]] for stage in summary["pipeline"]],
            ),
        ]
    if summary["warnings"]:
        blocks += [("p", "Warnings recorded during the run:"), ("ul", summary["warnings"])]
    blocks.append(("note", summary["validation_statement"]))
    return blocks


def _target_blocks(data: Mapping[str, Any]) -> list[tuple]:
    target = data["target"]
    structure = target["structure"] or {}
    box = target["binding_box"]
    resolution = structure.get("resolution_angstrom")
    rows = [
        ["Name", target["name"]],
        ["PDB ID", target["target_id"]],
        ["Source", f"{target['source']} ({target['source_url']})"],
        ["Structure title", structure.get("title") or "—"],
        ["Organism", structure.get("organism") or "—"],
        ["Method", structure.get("method") or "—"],
        ["Resolution", f"{resolution:g} Å" if resolution else "—"],
        ["Citation", f"{structure.get('citation_title') or '—'} doi:{structure.get('citation_doi') or '—'}"],
        ["Downloaded", f"{target['downloaded_on']} from {target['download_url']}"],
        ["Raw file", target["raw_file"]],
        ["Prepared receptor", target["preparation"]["prepared_receptor_file"]],
        ["Receptor chain", target["receptor_chain"]],
        ["Binding box centre (Å)", _format_triplet(box["center"])],
        ["Binding box size (Å)", _format_triplet(box["size"])],
        ["Binding box definition", box["definition"]],
    ]
    blocks = [("h", 2, "Target"), ("table", ["Field", "Value"], rows)]
    if target["preparation"]["steps"]:
        blocks += [("h", 3, "Preparation steps"), ("ol", target["preparation"]["steps"])]
    blocks.append(("p", "The target is used only as a computational demonstration; docking into it does not show biological activity."))
    return blocks


def _results_blocks(data: Mapping[str, Any]) -> list[tuple]:
    rank_data = data["ranking"]
    blocks: list[tuple] = [("h", 2, "Screening Results")]
    if not rank_data["ranked"]:
        reason = data["summary"]["docking_error"]
        blocks.append(
            (
                "p",
                "No compound received a Screening Priority Score in this run, so there are no top-ranked results. "
                + (f"Docking could not run: {reason}" if reason else "No compound was successfully docked."),
            )
        )
        return blocks

    blocks += [
        ("h", 3, "Top Ranked Compounds"),
        ("p", "Highest screening priority in this run, prioritised for further investigation. These are computational results."),
    ]
    for item in rank_data["top"]:
        blocks += [
            ("h", 4, f"#{item['rank']} · {item['name']} ({item['compound_id']})"),
            (
                "ul",
                [
                    f"Screening Priority Score: {_num(item['screening_priority_score'], config.SCORE_DECIMALS)} / 100",
                    f"Predicted AutoDock Vina affinity: {_num(item['predicted_affinity_kcal_mol'], config.AFFINITY_DECIMALS)} kcal/mol",
                    f"Docking points: {_num(item['docking_points'], 1)} / {rank_data['docking_max_points']:g}",
                    f"Property points: {_num(item['property_points'], 1)} / {rank_data['property_max_points']:g}",
                    f"Primary driver: {item['primary_driver']}",
                    f"Explanation: {item['explanation']}",
                ],
            ),
        ]
    blocks += [
        ("h", 3, "Full Ranking"),
        (
            "table",
            ["Rank", "Compound", "ID", "Score", "Predicted affinity (kcal/mol)", "Docking pts", "Property pts", "Primary driver"],
            [
                [
                    item["rank"], item["name"], item["compound_id"], _num(item["screening_priority_score"], config.SCORE_DECIMALS),
                    _num(item["predicted_affinity_kcal_mol"], config.AFFINITY_DECIMALS), _num(item["docking_points"], 1),
                    _num(item["property_points"], 1), item["primary_driver"],
                ]
                for item in rank_data["ranked"]
            ],
        ),
    ]
    return blocks


def _docking_blocks(data: Mapping[str, Any]) -> list[tuple]:
    blocks: list[tuple] = [("h", 2, "Docking Results")]
    if data["summary"]["docking_error"]:
        blocks.append(("p", f"Docking did not run in this session: {data['summary']['docking_error']}"))
        return blocks
    blocks += [
        ("p", "Predicted AutoDock Vina affinity of each compound's best pose. A computational estimate, not an experimental binding energy."),
        (
            "table",
            ["Compound", "ID", "Predicted affinity (kcal/mol)", "Status", "Poses", "Pose file"],
            [
                [
                    compound["name"], compound["compound_id"], _num(compound["docking"]["predicted_affinity_kcal_mol"], config.AFFINITY_DECIMALS),
                    compound["docking"]["status"], compound["docking"]["poses_generated"] or "—", compound["docking"]["pose_file"] or "—",
                ]
                for compound in data["compounds"]
            ],
        ),
    ]
    return blocks


def _descriptor_blocks(data: Mapping[str, Any]) -> list[tuple]:
    positions = {
        explanation["compound_id"]: {item["label"]: item["position"] for item in explanation["properties"]}
        for explanation in data["explanations"]
        if explanation["ranked"]
    }
    labels = [config.DESCRIPTORS[key]["label"] for key in DESCRIPTOR_KEYS]
    decimals = {key: descriptor["decimals"] for key, descriptor in config.DESCRIPTORS.items()}
    rows = []
    for compound in data["compounds"]:
        if not compound["valid"]:
            rows.append([compound["name"], *["—"] * len(labels), f"Invalid: {compound['descriptor_error']}"])
            continue
        outside = [
            f"{label} {position}"
            for label, position in positions.get(compound["compound_id"], {}).items()
            if position != "within configured range"
        ]
        note = ("; ".join(outside) if outside else "All within configured ranges") if compound["compound_id"] in positions else "Not ranked"
        rows.append(
            [compound["name"], *[_num(compound["descriptors"][key], decimals[key]) for key in DESCRIPTOR_KEYS], note]
        )
    return [
        ("h", 2, "Molecular Descriptors"),
        ("p", "Calculated by RDKit. Descriptors are molecular properties, not evidence of activity. Range notes compare against the configured screening heuristics."),
        ("table", ["Compound", *labels, "Configured-range notes"], rows),
    ]


def _unranked_blocks(data: Mapping[str, Any]) -> list[tuple]:
    unranked = data["ranking"]["unranked"]
    blocks: list[tuple] = [("h", 2, "Unranked Compounds")]
    if not unranked:
        return blocks + [("p", "Every processed compound was ranked.")]
    return blocks + [
        ("p", "These compounds have no Screening Priority Score. No score was calculated for them, and none was fabricated."),
        (
            "table",
            ["Compound", "ID", "Status", "Failure stage", "Recorded error"],
            [[item["name"], item["compound_id"], item["status"], item["failed_stage"] or "—", item["detail"] or "—"] for item in unranked],
        ),
    ]


def _explainability_blocks(data: Mapping[str, Any]) -> list[tuple]:
    ranked = [explanation for explanation in data["explanations"] if explanation["ranked"]]
    blocks: list[tuple] = [("h", 2, "Explainability Summary")]
    if not ranked:
        return blocks + [("p", "No compound was ranked, so there are no score explanations.")]
    blocks.append(("p", "How the scoring formula produced each score. This describes the calculation, not biology."))
    for explanation in sorted(ranked, key=lambda item: item["rank"]):
        items = [
            f"Contributions: docking {_num(explanation['docking_points'], 1)} + properties {_num(explanation['property_points'], 1)} points",
            f"Primary driver: {explanation['primary_driver']}",
            f"Main positive factor: {explanation['main_positive_factor']}",
            f"Main limiting factor: {explanation['main_limiting_factor']}",
        ]
        if explanation["comparison_below"]:
            items.append(f"Comparison: {explanation['comparison_below']}")
        # A tie note is shown once (on the higher compound) and not when it repeats the comparison.
        items.extend(
            f"Tie-break: {note}"
            for note in explanation["tie_notes"]
            if note.startswith(explanation["name"]) and note != explanation["comparison_below"]
        )
        blocks += [("h", 4, f"#{explanation['rank']} · {explanation['name']}"), ("ul", items)]
    return blocks


def _reproducibility_blocks(data: Mapping[str, Any]) -> list[tuple]:
    repro = data["reproducibility"]
    integrity = repro["integrity"]
    configuration = repro["configuration"]
    version_rows = [[name, value or "not detected"] for name, value in repro["software_versions"].items()]
    rows = [
        ["Run ID", repro["run_id"]],
        ["Configuration hash", repro["configuration_hash"]],
        ["Random seed", configuration["random_seed"]],
        ["Ranking weights", f"docking {configuration['ranking_weights']['docking']}, properties {configuration['ranking_weights']['properties']}"],
        ["Docking scale (kcal/mol)", f"0 at {configuration['docking_scale_kcal_mol']['zero_at']}, maximum at {configuration['docking_scale_kcal_mol']['full_at']}"],
        ["Vina settings", f"exhaustiveness {configuration['vina']['exhaustiveness']}, {configuration['vina']['num_modes']} poses"],
        ["Docking run folder", repro["docking_run_folder"] or "—"],
        ["Compound dataset", f"{integrity['compound_dataset']['path']} (SHA-256 {integrity['compound_dataset']['sha256']})"],
        ["Raw target", f"{integrity['raw_target']['path']} (SHA-256 {integrity['raw_target']['sha256']})"],
        ["Prepared receptor", f"{integrity['prepared_receptor']['path']} (SHA-256 {integrity['prepared_receptor']['sha256']})"],
    ]
    range_rows = [
        [
            config.DESCRIPTORS[key]["label"],
            f"{zero_below:g}" if _has_lower_limit((zero_below, full_from)) else "— (no lower limit)",
            f"{full_from:g}",
            f"{full_to:g}",
            f"{zero_above:g}",
        ]
        for key, (zero_below, full_from, full_to, zero_above) in configuration["property_ranges"].items()
    ]
    return [
        ("h", 2, "Reproducibility"),
        ("table", ["Field", "Value"], rows),
        ("p", integrity["prepared_receptor"]["note"]),
        ("h", 3, "Descriptor ranges"),
        ("table", ["Property", "0 at or below", "Full from", "Full to", "0 at or above"], range_rows),
        ("h", 3, "Software versions"),
        ("table", ["Software", "Version"], version_rows),
        *_timing_blocks(data["metadata"].get("timings_seconds") or {}),
    ]


def _timing_blocks(timings: Mapping[str, float]) -> list[tuple]:
    """Stage timings, clearly labelled as time-dependent (not part of the results)."""
    if not timings:
        return []
    rows = [
        [config.PIPELINE_STAGE_LABELS.get(key, "Total" if key == "total" else key), f"{seconds:.1f} s"]
        for key, seconds in timings.items()
    ]
    return [
        ("h", 3, "Stage timings"),
        ("p", "Wall-clock times from this run on this machine. They vary between runs and are not a benchmark."),
        ("table", ["Stage", "Time"], rows),
    ]


def build_document(data: Mapping[str, Any]) -> list[tuple]:
    """Order the report into blocks shared by the Markdown and HTML renderers."""
    blocks: list[tuple] = [("h", 1, data["metadata"]["title"])]
    blocks += _summary_blocks(data)
    blocks += _target_blocks(data)
    blocks.append(("h", 2, "Methodology"))
    for method in data["methodology"]:
        blocks += [("h", 3, method["title"]), ("p", method["text"])]
    blocks += _results_blocks(data)
    blocks += _docking_blocks(data)
    blocks += _descriptor_blocks(data)
    blocks += _unranked_blocks(data)
    blocks += _explainability_blocks(data)
    blocks += [("h", 2, "Limitations"), ("ul", data["limitations"])]
    blocks += [("h", 2, "Responsible Use"), ("ul", data["responsible_use"])]
    blocks += _reproducibility_blocks(data)
    blocks.append(("note", VALIDATION_STATEMENT))
    return blocks


# --- Renderers -------------------------------------------------------------------

def _md_cell(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def render_markdown(blocks: Sequence[tuple]) -> str:
    """Render document blocks as GitHub-flavoured Markdown."""
    lines: list[str] = []
    for block in blocks:
        kind = block[0]
        if kind == "h":
            lines += [f"{'#' * block[1]} {block[2]}", ""]
        elif kind == "p":
            lines += [block[1], ""]
        elif kind == "note":
            lines += [f"> **Note:** {block[1]}", ""]
        elif kind in ("ul", "ol"):
            marker = "-" if kind == "ul" else None
            lines += [f"{marker or f'{number}.'} {item}" for number, item in enumerate(block[1], start=1)] + [""]
        elif kind == "table":
            headers, rows = block[1], block[2]
            lines.append("| " + " | ".join(_md_cell(header) for header in headers) + " |")
            lines.append("| " + " | ".join("---" for _ in headers) + " |")
            lines += ["| " + " | ".join(_md_cell(cell) for cell in row) + " |" for row in rows]
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"


HTML_STYLE = """
body { font-family: system-ui, -apple-system, "Segoe UI", sans-serif; color: #1F2933; max-width: 1000px; margin: 2rem auto; padding: 0 1rem; line-height: 1.5; }
h1 { letter-spacing: 0.04em; } h2 { border-bottom: 1px solid #D9E2E7; padding-bottom: 0.2rem; margin-top: 2rem; }
table { border-collapse: collapse; width: 100%; margin: 0.5rem 0 1rem; font-size: 0.9rem; }
th, td { border: 1px solid #D9E2E7; padding: 0.35rem 0.5rem; text-align: left; vertical-align: top; }
th { background: #F3F6F8; } blockquote { border-left: 4px solid #0F766E; margin: 1rem 0; padding: 0.5rem 1rem; background: #F3F6F8; }
td, p, li { overflow-wrap: break-word; } .table-wrap { overflow-x: auto; } .table-wrap table { min-width: 560px; }
@media (max-width: 640px) { body { margin: 1rem auto; } table { font-size: 0.8rem; } }
"""


def render_html(blocks: Sequence[tuple], title: str) -> str:
    """Render document blocks as a standalone HTML page (all text escaped)."""
    escape = html.escape
    parts = [
        "<!doctype html>",
        '<html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{escape(title)}</title><style>{HTML_STYLE}</style></head><body>",
    ]
    for block in blocks:
        kind = block[0]
        if kind == "h":
            parts.append(f"<h{block[1]}>{escape(str(block[2]))}</h{block[1]}>")
        elif kind == "p":
            parts.append(f"<p>{escape(str(block[1]))}</p>")
        elif kind == "note":
            parts.append(f"<blockquote><strong>Note:</strong> {escape(str(block[1]))}</blockquote>")
        elif kind in ("ul", "ol"):
            items = "".join(f"<li>{escape(str(item))}</li>" for item in block[1])
            parts.append(f"<{kind}>{items}</{kind}>")
        elif kind == "table":
            head = "".join(f"<th>{escape(str(header))}</th>" for header in block[1])
            body = "".join("<tr>" + "".join(f"<td>{escape(str(cell))}</td>" for cell in row) + "</tr>" for row in block[2])
            parts.append(f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>')
    parts.append("</body></html>")
    return "\n".join(parts) + "\n"


# --- Writing ---------------------------------------------------------------------

def report_paths(directory: Path, now: datetime, formats: Sequence[str]) -> dict[str, Path]:
    """Windows-safe filenames that do not collide with existing reports.

    All formats of one report share a stem; ``_2``, ``_3``… is added if any
    file with that stem already exists.
    """
    stem = safe_filename_part(f"fairouz_screening_{config.TARGET_ID}_{now:%Y-%m-%d_%H%M%S}")
    for attempt in range(1, 1000):
        candidate = stem if attempt == 1 else f"{stem}_{attempt}"
        paths = {fmt: directory / f"{candidate}.{fmt}" for fmt in formats}
        if not any(path.exists() for path in paths.values()):
            return paths
    raise FileExistsError(f"Could not find a free report filename in {directory}.")


def build_preview(data: Mapping[str, Any]) -> dict[str, Any]:
    """The short summary shown in the app."""
    summary = data["summary"]
    return {
        "run_id": data["metadata"]["run_id"],
        "target": summary["target"],
        "processed": summary["processed"],
        "ranked": summary["ranked"],
        "unranked": summary["unranked"],
        "outcome": summary["outcome"],
        "top": [
            {
                "rank": item["rank"],
                "name": item["name"],
                "score": item["screening_priority_score"],
                "predicted_affinity_kcal_mol": item["predicted_affinity_kcal_mol"],
            }
            for item in data["ranking"]["top"]
        ],
        "key_limitations": list(data["limitations"][:4]),
    }


def _render(data: Mapping[str, Any], blocks: Sequence[tuple], fmt: str) -> str:
    if fmt == "md":
        return render_markdown(blocks)
    if fmt == "json":
        return json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if fmt == "html":
        return render_html(blocks, f"{data['metadata']['title']} · {data['metadata']['run_id']}")
    raise ValueError(f"Unsupported report format: {fmt}")


def generate_report(
    response: Mapping[str, Any],
    output_dir: Path | None = None,
    formats: Sequence[str] | None = None,
    now: datetime | None = None,
) -> dict[str, Any]:
    """Build and save the report for a screening response.

    Never raises: problems are recorded in the result, so the screening
    results already produced are never lost.

    Args:
        response: The screening response (not modified).
        output_dir: Where to write. Defaults to ``config.REPORTS_DIR``.
        formats: Any of ``"md"``, ``"json"``, ``"html"``. Defaults to all.
        now: Report timestamp (UTC), mainly for tests.

    Returns:
        ``stage``, ``status`` (generated / partial / not_generated),
        ``run_id``, ``files`` (format → path), ``errors`` (format → message),
        ``preview`` and ``markdown`` (the Markdown text, for previews).
    """
    formats = tuple(formats or config.REPORT_FORMATS)
    now = now or datetime.now(timezone.utc)
    result: dict[str, Any] = {
        "stage": "report", "status": config.REPORT_NOT_GENERATED, "run_id": None,
        "files": {}, "errors": {}, "preview": None, "markdown": None,
    }
    try:
        data = build_report_data(response, now)
        blocks = build_document(data)
        result.update(run_id=data["metadata"]["run_id"], preview=build_preview(data), markdown=render_markdown(blocks))
    except Exception as error:  # Report problems must never break the pipeline.
        logger.exception("Report could not be built")
        result["errors"]["build"] = f"Report could not be built: {error}"
        return result

    directory = output_dir or config.REPORTS_DIR
    try:
        directory.mkdir(parents=True, exist_ok=True)
        paths = report_paths(directory, now, formats)
    except (OSError, FileExistsError) as error:
        logger.error("Report folder unavailable: %s", error)
        result["errors"]["folder"] = f"Report folder unavailable: {error}"
        return result

    for fmt, path in paths.items():
        try:
            text = result["markdown"] if fmt == "md" else _render(data, blocks, fmt)  # Markdown is rendered once.
            path.write_text(text, encoding="utf-8", newline="\n")
            result["files"][fmt] = str(path)
        except (OSError, ValueError) as error:
            logger.error("Could not write %s report: %s", fmt, error)
            result["errors"][fmt] = str(error)

    if len(result["files"]) == len(formats):
        result["status"] = config.REPORT_GENERATED
    elif result["files"]:
        result["status"] = config.REPORT_PARTIAL
    logger.info("Report %s: %s", result["status"], ", ".join(result["files"].values()) or "no files")
    return result
