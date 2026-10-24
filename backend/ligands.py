"""Ligand preparation: SMILES → 3D structure → PDBQT for AutoDock Vina.

For each compound:

1. Parse and validate the SMILES with RDKit (shared with ``descriptors``).
2. Add explicit hydrogens.
3. Embed a 3D conformer with ETKDGv3 using a fixed random seed.
4. Optimise it with the MMFF94 force field, or UFF if MMFF lacks
   parameters for the molecule.
5. Convert it to PDBQT (atom types, charges, rotatable bonds) with Meeko.

The same SMILES and seed always give the same PDBQT. Molecules are used in
the protonation state written in their SMILES; no pH-dependent protonation
or tautomer enumeration is done.
"""

import logging

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

import config
from backend.descriptors import parse_smiles

logger = logging.getLogger(__name__)

FORCE_FIELD_MAX_ITERATIONS = 2000


class LigandPreparationError(Exception):
    """Raised when a valid molecule cannot be turned into a docking-ready ligand."""


def generate_conformer(molecule: Chem.Mol, seed: int) -> Chem.Mol:
    """Add hydrogens and embed one 3D conformer.

    Args:
        molecule: A parsed RDKit molecule (not modified).
        seed: Random seed, so embedding is reproducible.

    Returns:
        A new molecule with explicit hydrogens and one 3D conformer.

    Raises:
        LigandPreparationError: If RDKit cannot embed the molecule in 3D.
    """
    with_hydrogens = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    if AllChem.EmbedMolecule(with_hydrogens, parameters) != 0:
        raise LigandPreparationError("RDKit could not generate a 3D conformer.")
    return with_hydrogens


def optimise_conformer(molecule: Chem.Mol) -> str:
    """Optimise the conformer in place with MMFF94, falling back to UFF.

    Returns:
        The name of the force field used.

    Raises:
        LigandPreparationError: If neither force field can handle the molecule.
    """
    if AllChem.MMFFHasAllMoleculeParams(molecule):
        force_field, outcome = "MMFF94", AllChem.MMFFOptimizeMolecule(
            molecule, maxIters=FORCE_FIELD_MAX_ITERATIONS
        )
    elif AllChem.UFFHasAllMoleculeParams(molecule):
        force_field, outcome = "UFF", AllChem.UFFOptimizeMolecule(
            molecule, maxIters=FORCE_FIELD_MAX_ITERATIONS
        )
    else:
        raise LigandPreparationError("No force field (MMFF94 or UFF) has parameters for this molecule.")

    # 0 = converged, 1 = more iterations needed, -1 = setup failed.
    if outcome == -1:
        raise LigandPreparationError(f"{force_field} optimisation could not be set up.")
    if outcome == 1:
        logger.warning("%s optimisation did not fully converge; using the last geometry.", force_field)
    return force_field


def molecule_to_pdbqt(molecule: Chem.Mol) -> str:
    """Convert a 3D molecule with hydrogens into a Vina PDBQT string with Meeko.

    Raises:
        LigandPreparationError: If Meeko cannot prepare the molecule.
    """
    try:
        setups = MoleculePreparation().prepare(molecule)
    except Exception as error:  # Meeko raises various errors for unusual chemistry.
        raise LigandPreparationError(f"Meeko could not prepare the ligand: {error}") from error
    if len(setups) != 1:
        raise LigandPreparationError("Expected one molecule; multi-fragment SMILES are not supported.")

    pdbqt, is_ok, error_message = PDBQTWriterLegacy.write_string(setups[0])
    if not is_ok:
        raise LigandPreparationError(f"Meeko could not write PDBQT: {error_message}")
    return pdbqt


def prepare_ligand(smiles: str, seed: int = config.RANDOM_SEED) -> str:
    """Turn a SMILES string into a docking-ready PDBQT string.

    Args:
        smiles: The compound's SMILES.
        seed: Random seed for conformer generation.

    Returns:
        The ligand in PDBQT format.

    Raises:
        InvalidSmilesError: If the SMILES cannot be parsed.
        LigandPreparationError: If 3D generation, optimisation or conversion fails.
    """
    molecule = generate_conformer(parse_smiles(smiles), seed)
    force_field = optimise_conformer(molecule)
    logger.debug("Prepared ligand %s with %s", smiles, force_field)
    return molecule_to_pdbqt(molecule)
