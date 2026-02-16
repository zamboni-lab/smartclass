"""Enumerate structural variants of molecules."""

from __future__ import annotations

from rdkit.Chem import Mol, rdMolEnumerator

from smartclass.logging import get_logger


__all__ = [
    "enumerate_structures",
]

logger = get_logger(__name__)


def enumerate_structures(mol: Mol) -> list[Mol]:
    """
    Enumerate structural variants of a molecule.

    Uses RDKit's MolEnumerator to generate structural variants
    (e.g., for handling tautomers or stereoisomers in queries).

    :param mol: RDKit Mol object to enumerate.
    :returns: List of enumerated Mol objects. Falls back to the original
        molecule if enumeration fails or produces no results.
    """
    try:
        enumerated_mols = rdMolEnumerator.Enumerate(mol)
    except Exception as e:
        logger.debug(f"Enumeration failed: {e}")
        return [mol]

    if not enumerated_mols:
        return [mol]

    return list(enumerated_mols)
