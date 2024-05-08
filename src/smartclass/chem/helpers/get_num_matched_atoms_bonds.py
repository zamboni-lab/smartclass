"""Get number of matched atoms and bonds."""

from __future__ import annotations

from rdkit.Chem import Mol, rdFMCS

__all__ = [
    "get_num_matched_atoms_bonds",
]


def get_num_matched_atoms_bonds(
    mol_1: Mol,
    mol_2: Mol,
) -> int:
    """
    Get number of matched atoms and bonds.

    :param mol_1: A structure MOL.
    :type mol_1: Mol

    :param mol_2: A structure MOL.
    :type mol_2: Mol

    :returns: Number of matched atoms and bonds.
    :rtype: int

    """
    mcs = rdFMCS.FindMCS(
        [mol_1, mol_2],
        atomCompare=rdFMCS.AtomCompare.CompareAny,
        bondCompare=rdFMCS.BondCompare.CompareAny,
    )

    return mcs.numAtoms + mcs.numBonds
