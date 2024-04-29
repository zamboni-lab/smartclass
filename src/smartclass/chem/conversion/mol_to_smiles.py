"""Convert a structure MOL to SMILES."""

from __future__ import annotations

from rdkit.Chem import MolToSmiles

__all__ = [
    "mol_to_smiles",
]


def mol_to_smiles(mol: str) -> str:
    """
    Convert a structure MOL to SMILES.

    :param mol: A structure MOL.
    :type mol: str

    :returns: A SMILES.
    :rtype: str
    """
    return MolToSmiles(mol)
