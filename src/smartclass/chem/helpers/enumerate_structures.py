"""Enumerate structures."""

from __future__ import annotations

import logging

from rdkit.Chem import Mol, rdMolEnumerator

__all__ = [
    "enumerate_structures",
]


def enumerate_structures(
    mol: Mol,
) -> list[Mol]:
    """
    Enumerate structures.

    :param mol: A structure MOL.
    :type mol: Mol


    :returns: A list of enumerated molecules.
    :rtype: list[Mol]
    """
    try:
        enumerated_mols = rdMolEnumerator.Enumerate(mol)
    except Exception as e:
        logging.error(f"Enumeration failed: {e}")
        enumerated_mols = []

    # Check if enumeration was successful
    if not enumerated_mols:
        # If enumeration failed, use the original mol as a fallback
        enumerated_mols = [mol]

    return enumerated_mols
