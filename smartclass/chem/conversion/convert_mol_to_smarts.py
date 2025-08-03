"""Convert a structure MOL to SMILES."""

from __future__ import annotations

from rdkit.Chem import Mol, MolToSmarts

__all__ = [
    "convert_mol_to_smarts",
]


def convert_mol_to_smarts(mol: Mol) -> str:
    """
    Convert a structure MOL to SMARTS.

    :param mol: A structure MOL.
    :type mol: Mol

    :returns: A SMARTS.
    :rtype: str
    """
    return MolToSmarts(mol)
