"""Convert a structure MOL to CXSMILES."""

from __future__ import annotations

from rdkit.Chem import MolToCXSmiles


__all__ = [
    "convert_mol_to_cxsmiles",
]


def convert_mol_to_cxsmiles(mol: str) -> str:
    """
    Convert a structure MOL to CXSMILES.

    :param mol: A structure MOL.
    :type mol: str

    :returns: A CXSMILES.
    :rtype: str
    """
    return MolToCXSmiles(mol)
