"""Convert a structure MOL to CXSMILES."""

from __future__ import annotations

from rdkit.Chem import Mol, MolToCXSmiles

__all__ = [
    "convert_mol_to_cxsmiles",
]


def convert_mol_to_cxsmiles(mol: Mol) -> str:
    """Convert a structure MOL to CXSMILES.

    Parameters
    ----------
    mol : Mol
        MOL.

    Returns
    -------
    str
        CXSMILES.
    """
    return MolToCXSmiles(mol)
