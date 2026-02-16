"""Convert a structure MOL to InChI."""

from __future__ import annotations

from rdkit.Chem import MolToInchi


__all__ = [
    "convert_mol_to_inchi",
]


def convert_mol_to_inchi(mol: str) -> str:
    """
    Convert a structure MOL to InChI.

    :param mol: A structure MOL.
    :type mol: str

    :returns: An InChI.
    :rtype: str
    """
    return MolToInchi(mol)
