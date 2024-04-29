"""Convert a structure MOL to InChI."""

from __future__ import annotations

from rdkit.Chem import MolToInchi

__all__ = [
    "mol_to_inchi",
]


def mol_to_inchi(mol: str) -> str:
    """
    Convert a structure MOL to InChI.

    :param mol: A structure MOL.
    :type mol: str

    :returns: An InChI.
    :rtype: str
    """
    return MolToInchi(mol)
