"""Convert a structure MOL to mass."""

from __future__ import annotations

from rdkit.Chem import Mol
from rdkit.Chem.Descriptors import ExactMolWt

__all__ = [
    "convert_mol_to_mass",
]


def convert_mol_to_mass(mol: Mol) -> float:
    """
    Convert a structure MOL to mass.

    :param mol: A structure MOL.
    :type mol: Mol

    :returns: A mass.
    :rtype: float
    """
    return ExactMolWt(mol)
