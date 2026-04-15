"""Convert a structure MOL to SMILES."""

from __future__ import annotations

from rdkit.Chem import Mol, MolToSmarts


__all__ = [
    "convert_mol_to_smarts",
]


def convert_mol_to_smarts(mol: Mol) -> str:
    """Convert a structure MOL to SMARTS.

    Parameters
    ----------
    mol : Mol
        MOL.

    Returns
    -------
    str
        SMARTS.
    """
    return MolToSmarts(mol)
