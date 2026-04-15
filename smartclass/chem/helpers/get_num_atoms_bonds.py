"""Get number of atoms and bonds."""

from __future__ import annotations

from rdkit.Chem import Mol


__all__ = [
    "get_num_atoms_bonds",
]


def get_num_atoms_bonds(
    mol: Mol,
) -> int:
    """Get number of atoms and bonds.

    Parameters
    ----------
    mol : Mol
        MOL.

    Returns
    -------
    int
        Number of atoms and bonds.
    """
    return mol.GetNumAtoms() + mol.GetNumBonds()
