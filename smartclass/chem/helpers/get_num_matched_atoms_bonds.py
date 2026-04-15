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
    """Get number of matched atoms and bonds.

    Parameters
    ----------
    mol_1 : Mol
        MOL.
    mol_2 : Mol
        MOL.

    Returns
    -------
    int
        Number of matched atoms and bonds.
    """
    # TODO See if using the RASCAL MCES
    mcs = rdFMCS.FindMCS(
        [mol_1, mol_2],
        atomCompare=rdFMCS.AtomCompare.CompareAny,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        # ringCompare=rdFMCS.ringCompare.IgnoreRingFusion,
    )

    return mcs.numAtoms + mcs.numBonds
