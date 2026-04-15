"""Convert a structure MOL to formula."""

from __future__ import annotations

from rdkit.Chem import Mol
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


__all__ = [
    "convert_mol_to_formula",
]


def convert_mol_to_formula(mol: Mol) -> str:
    """Convert a structure MOL to formula.

    Parameters
    ----------
    mol : Mol
        MOL.

    Returns
    -------
    str
        A formula.
    """
    return CalcMolFormula(mol)
