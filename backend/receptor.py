"""Receptor (protein target) preparation and binding-box calculation.

Turns the raw 6LU7 file downloaded from the RCSB Protein Data Bank into the
rigid receptor file AutoDock Vina needs. The raw file is only ever read,
never modified. Every transformation is listed in ``PREPARATION_STEPS`` and
recorded in a JSON log next to the prepared files.

Run once (the prepared files are committed to the repository)::

    python -m backend.receptor

Preparation is kept separate from docking: docking only reads the prepared
receptor and never changes it.
"""

import hashlib
import json
import logging
import math
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import gemmi

import config

logger = logging.getLogger(__name__)

# Human-readable record of what prepare_receptor() does, written to the log.
PREPARATION_STEPS = (
    "Read the raw RCSB PDB file without modifying it; verify its SHA-256 checksum.",
    "Keep model 1 and only the protease chain (chain A).",
    "Remove the co-crystallised inhibitor N3 (chain C), which is covalently "
    "bonded to Cys145 in the crystal; Cys145 is kept as a free cysteine.",
    "Remove all 84 crystallographic waters and any other non-protein residues.",
    "Remove any hydrogens present so they are rebuilt consistently.",
    "Write the cleaned heavy-atom receptor as PDB.",
    "Run Meeko mk_prepare_receptor: match each residue to a template, add "
    "polar hydrogens, assign Gasteiger charges and AutoDock atom types, and "
    "write a rigid receptor PDBQT. Histidines use Meeko's default HIE "
    "protonation state.",
)


class ReceptorPreparationError(Exception):
    """Raised when the target structure cannot be read or prepared."""


@dataclass(frozen=True)
class BindingBox:
    """A docking search box in ångströms.

    Attributes:
        center: (x, y, z) coordinates of the box centre.
        size: Edge lengths along x, y and z.
    """

    center: tuple[float, float, float]
    size: tuple[float, float, float]


def validate_binding_box(box: BindingBox) -> None:
    """Check that a binding box is usable by Vina.

    Raises:
        ValueError: If the centre or size is malformed, not finite, or a size
            is not positive.
    """
    if len(box.center) != 3 or len(box.size) != 3:
        raise ValueError("Binding box centre and size must each have three values (x, y, z).")
    if not all(math.isfinite(value) for value in (*box.center, *box.size)):
        raise ValueError("Binding box values must be finite numbers.")
    if any(edge <= 0 for edge in box.size):
        raise ValueError(f"Binding box edges must be positive, got {box.size}.")


def calculate_binding_box(
    coordinates: Sequence[tuple[float, float, float]],
    padding: float,
) -> BindingBox:
    """Derive a docking box that encloses a reference ligand.

    The centre is the midpoint of the ligand's extent on each axis; each
    edge is that extent plus ``padding`` on both sides.

    Args:
        coordinates: Atom coordinates of the reference ligand, in ångströms.
        padding: Space added on each side of the ligand, in ångströms.

    Returns:
        The binding box, rounded to 4 decimal places.

    Raises:
        ValueError: If no coordinates are given or the padding is negative.
    """
    if not coordinates:
        raise ValueError("Cannot derive a binding box from zero atoms.")
    if padding < 0:
        raise ValueError(f"Padding must not be negative, got {padding}.")

    axes = list(zip(*coordinates))
    center = tuple(round((min(axis) + max(axis)) / 2, 4) for axis in axes)
    size = tuple(round(max(axis) - min(axis) + 2 * padding, 4) for axis in axes)
    box = BindingBox(center=center, size=size)
    validate_binding_box(box)
    return box


def file_sha256(path: Path) -> str:
    """Return the SHA-256 hex digest of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_raw_structure(path: Path, expected_sha256: str) -> None:
    """Check that the raw structure exists and is the file that was downloaded.

    Raises:
        ReceptorPreparationError: If the file is missing or its checksum differs.
    """
    if not path.is_file():
        raise ReceptorPreparationError(f"Raw target structure not found: {path}")
    actual = file_sha256(path)
    if actual != expected_sha256:
        raise ReceptorPreparationError(
            f"Raw target {path.name} has changed since download "
            f"(SHA-256 {actual}, expected {expected_sha256})."
        )


def read_structure(path: Path) -> gemmi.Structure:
    """Read a PDB file with gemmi, raising a clear error if it fails."""
    try:
        return gemmi.read_structure(str(path))
    except (OSError, RuntimeError, ValueError) as error:
        raise ReceptorPreparationError(f"Could not read structure {path}: {error}") from error


def reference_ligand_coordinates(path: Path, chain: str) -> list[tuple[float, float, float]]:
    """Return the coordinates of every atom in the reference ligand chain.

    Raises:
        ReceptorPreparationError: If the chain is absent or empty.
    """
    model = read_structure(path)[0]
    if model.find_chain(chain) is None:
        raise ReceptorPreparationError(f"Chain {chain} not found in {path.name}.")
    coordinates = [(atom.pos.x, atom.pos.y, atom.pos.z) for residue in model[chain] for atom in residue]
    if not coordinates:
        raise ReceptorPreparationError(f"Chain {chain} in {path.name} has no atoms.")
    return coordinates


def clean_receptor(raw_path: Path, output_path: Path, chain: str) -> dict[str, int]:
    """Write a heavy-atom, protein-only copy of one chain.

    Args:
        raw_path: The raw structure (read only).
        output_path: Where to write the cleaned PDB.
        chain: The protein chain to keep.

    Returns:
        Counts describing the result: residues, atoms and waters removed.
    """
    structure = read_structure(raw_path)
    del structure[1:]  # Keep model 1 only.
    model = structure[0]
    if model.find_chain(chain) is None:
        raise ReceptorPreparationError(f"Chain {chain} not found in {raw_path.name}.")

    waters_removed = sum(residue.is_water() for residue in model[chain])
    for other in [c.name for c in model if c.name != chain]:
        model.remove_chain(other)
    structure.remove_ligands_and_waters()
    structure.remove_hydrogens()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    structure.write_pdb(str(output_path))
    return {
        "residues": len(model[chain]),
        "heavy_atoms": model.count_atom_sites(),
        "waters_removed": waters_removed,
    }


def convert_to_pdbqt(clean_path: Path, pdbqt_path: Path) -> None:
    """Convert a cleaned receptor PDB to PDBQT with Meeko's receptor tool.

    Raises:
        ReceptorPreparationError: If Meeko fails or writes no output.
    """
    command = [
        sys.executable, "-m", "meeko.cli.mk_prepare_receptor",
        "--read_pdb", str(clean_path),
        "-o", str(pdbqt_path.with_suffix("")),
        "-p",  # Write a PDBQT file.
    ]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=600)
    if completed.returncode != 0 or not pdbqt_path.is_file():
        raise ReceptorPreparationError(
            f"Meeko receptor preparation failed: {completed.stdout[-500:]} {completed.stderr[-500:]}"
        )


def prepare_receptor() -> dict:
    """Prepare the configured target and write a preparation log.

    Returns:
        The preparation log, also saved to ``config.PREPARATION_LOG_FILE``.
    """
    verify_raw_structure(config.RAW_TARGET_FILE, config.RAW_TARGET_SHA256)
    logger.info("Preparing receptor %s from %s", config.TARGET_ID, config.RAW_TARGET_FILE.name)

    counts = clean_receptor(config.RAW_TARGET_FILE, config.CLEAN_RECEPTOR_FILE, config.RECEPTOR_CHAIN)
    convert_to_pdbqt(config.CLEAN_RECEPTOR_FILE, config.PREPARED_RECEPTOR_FILE)
    box = calculate_binding_box(
        reference_ligand_coordinates(config.RAW_TARGET_FILE, config.REFERENCE_LIGAND_CHAIN),
        config.BINDING_BOX_PADDING,
    )

    log = {
        "target_id": config.TARGET_ID,
        "prepared_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "raw_file": config.RAW_TARGET_FILE.name,
        "raw_sha256": config.RAW_TARGET_SHA256,
        "receptor_chain": config.RECEPTOR_CHAIN,
        "steps": list(PREPARATION_STEPS),
        "clean_receptor": counts,
        "prepared_receptor_sha256": file_sha256(config.PREPARED_RECEPTOR_FILE),
        "binding_box": {**asdict(box), "padding": config.BINDING_BOX_PADDING, "units": "angstrom"},
        "tool_versions": {"gemmi": gemmi.__version__, "meeko": version("meeko")},
    }
    config.PREPARATION_LOG_FILE.write_text(json.dumps(log, indent=2), encoding="utf-8")
    logger.info("Receptor prepared: %s", config.PREPARED_RECEPTOR_FILE)
    return log


if __name__ == "__main__":
    config.configure_logging()
    result = prepare_receptor()
    print(json.dumps(result["binding_box"], indent=2))
