"""Convert a structure MOL to SMILES."""

from __future__ import annotations

from rdkit.Chem import Mol, MolToSmiles


__all__ = [
    "convert_mol_to_smiles",
]


def convert_mol_to_smiles(mol: Mol) -> str:
    """Convert a structure MOL to SMILES.

Parameters
----------
mol : Mol
    MOL.

Returns
-------
str
    SMILES.
    """
    return MolToSmiles(mol)
