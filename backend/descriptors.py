"""Molecular descriptor calculation with RDKit.

Each compound goes through:

    validate SMILES → parse into an RDKit Mol → calculate descriptors
    → structured result

A molecule that cannot be parsed is marked invalid, logged, and skipped; it
never stops the rest of the batch. Descriptor values are only ever produced
by RDKit.

Definitions used (other tools, such as PubChem, use slightly different
definitions, so their values can differ for the same molecule):

* Molecular weight: average molecular weight (``Descriptors.MolWt``).
* LogP: Wildman–Crippen estimate (``Crippen.MolLogP``).
* TPSA: Ertl topological polar surface area, N and O only (``CalcTPSA``).
* H-bond donors / acceptors: RDKit's Lipinski-style counts
  (``CalcNumHBD`` / ``CalcNumHBA``).
"""

import logging
from collections.abc import Callable, Mapping, Sequence
from typing import Any

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

logger = logging.getLogger(__name__)

# RDKit writes its own parse errors straight to stderr. We log a clearer
# message per compound instead.
RDLogger.DisableLog("rdApp.*")

# Descriptor key -> RDKit function. Keys match config.DESCRIPTORS.
DESCRIPTOR_CALCULATORS: dict[str, Callable[[Chem.Mol], float | int]] = {
    "molecular_weight": Descriptors.MolWt,
    "logp": Crippen.MolLogP,
    "tpsa": rdMolDescriptors.CalcTPSA,
    "h_bond_donors": rdMolDescriptors.CalcNumHBD,
    "h_bond_acceptors": rdMolDescriptors.CalcNumHBA,
}


class InvalidSmilesError(ValueError):
    """Raised when a SMILES string cannot be turned into a molecule."""


def parse_smiles(smiles: str | None) -> Chem.Mol:
    """Validate a SMILES string and convert it into an RDKit molecule.

    Args:
        smiles: The SMILES string to parse.

    Returns:
        The parsed RDKit molecule.

    Raises:
        InvalidSmilesError: If the string is empty or RDKit cannot parse it.
    """
    cleaned = str(smiles).strip() if smiles is not None else ""
    if not cleaned:
        raise InvalidSmilesError("SMILES string is empty.")

    molecule = Chem.MolFromSmiles(cleaned)
    # MolFromSmiles returns None on failure; a zero-atom Mol is not a molecule.
    if molecule is None or molecule.GetNumAtoms() == 0:
        raise InvalidSmilesError(f"RDKit could not parse SMILES '{cleaned}'.")
    return molecule


def calculate_molecule_descriptors(molecule: Chem.Mol) -> dict[str, float | int]:
    """Calculate every configured descriptor for one RDKit molecule."""
    return {name: calculate(molecule) for name, calculate in DESCRIPTOR_CALCULATORS.items()}


def analyse_compound(compound: Mapping[str, str]) -> dict[str, Any]:
    """Calculate descriptors for one compound record.

    Args:
        compound: A record with ``compound_id``, ``name``, ``smiles`` and
            optionally ``source``. It is not modified.

    Returns:
        A new dictionary with the compound's identity, ``valid``, ``error``
        and one key per descriptor. Descriptor values are ``None`` and
        ``error`` explains why when the compound is invalid.
    """
    result: dict[str, Any] = {
        "compound_id": compound.get("compound_id", ""),
        "name": compound.get("name", ""),
        "smiles": compound.get("smiles", ""),
        "source": compound.get("source", ""),
        "valid": False,
        "error": None,
        **dict.fromkeys(DESCRIPTOR_CALCULATORS),
    }

    try:
        molecule = parse_smiles(result["smiles"])
        result.update(calculate_molecule_descriptors(molecule))
    except InvalidSmilesError as error:
        result["error"] = str(error)
        logger.warning("Compound %s is invalid: %s", result["compound_id"], error)
    except Exception:
        # RDKit can raise on unusual structures; one molecule must not stop the batch.
        result["error"] = "Descriptor calculation failed unexpectedly."
        logger.exception("Descriptor calculation failed for compound %s", result["compound_id"])
    else:
        result["valid"] = True

    return result


def calculate_descriptors(compounds: Sequence[Mapping[str, str]]) -> list[dict[str, Any]]:
    """Calculate RDKit descriptors for a batch of compounds.

    Args:
        compounds: Compound records (see ``analyse_compound``).

    Returns:
        One result per compound, in the same order. Invalid compounds are
        included and marked with ``valid=False``.
    """
    results = [analyse_compound(compound) for compound in compounds]
    valid_count = sum(result["valid"] for result in results)
    logger.info(
        "Descriptors calculated for %d compound(s): %d valid, %d invalid.",
        len(results),
        valid_count,
        len(results) - valid_count,
    )
    return results
