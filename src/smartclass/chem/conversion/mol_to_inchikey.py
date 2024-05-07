"""Convert a structure MOL to InChIKey."""

from __future__ import annotations

from rdkit.Chem import MolToInchiKey

__all__ = [
    "mol_to_inchikey",
]


def mol_to_inchikey(mol: str) -> str:
    """
    Convert a structure MOL to InChIKey.

    :param mol: A structure MOL.
    :type mol: str

    :returns: An InChIKey.
    :rtype: str
    """
    return MolToInchiKey(mol)
