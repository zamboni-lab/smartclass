"""Convert a structure MOL to SMILES."""

from __future__ import annotations

from rdkit.Chem import Mol, MolToSmiles


__all__ = [
    "convert_mol_to_smiles",
]


def convert_mol_to_smiles(mol: Mol) -> str:
    """
    Convert a structure MOL to SMILES.

    :param mol: A structure MOL.
    :type mol: Mol

    :returns: A SMILES.
    :rtype: str
    """
    return MolToSmiles(mol)
