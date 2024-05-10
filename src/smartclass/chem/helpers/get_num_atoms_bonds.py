"""Get number of atoms and bonds."""

from __future__ import annotations

from rdkit.Chem import Mol

__all__ = [
    "get_num_atoms_bonds",
]


def get_num_atoms_bonds(
    mol: Mol,
) -> int:
    """
    Get number of atoms and bonds.

    :param mol: A structure MOL.
    :type mol: Mol

    :returns: Number of atoms and bonds.
    :rtype: int

    """
    return mol.GetNumAtoms() + mol.GetNumBonds()
