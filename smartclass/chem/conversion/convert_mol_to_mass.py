"""Convert a structure MOL to mass."""

from __future__ import annotations

from rdkit.Chem import Mol
from rdkit.Chem.Descriptors import ExactMolWt


__all__ = [
    "convert_mol_to_mass",
]


def convert_mol_to_mass(mol: Mol) -> float:
    """Convert a structure MOL to mass.

    Parameters
    ----------
    mol : Mol
        MOL.

    Returns
    -------
    float
        A mass.
    """
    return ExactMolWt(mol)
