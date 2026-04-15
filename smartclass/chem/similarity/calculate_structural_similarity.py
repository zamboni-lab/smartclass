"""Calculate structural similarity."""

from __future__ import annotations

from rdkit.Chem import Mol, rdRascalMCES


opts = rdRascalMCES.RascalOptions()
opts.returnEmptyMCES = True
opts.similarityThreshold = -1
opts.timeout = 0


def calculate_structural_similarity(
    mol_1: Mol,
    mol_2: Mol,
    options: rdRascalMCES.RascalOptions = opts,
) -> list[rdRascalMCES.RascalResult]:
    """Calculate structural similarity.

Parameters
----------
mol_1 : Mol
    Mol 1.
mol_2 : Mol
    Mol 2.
options : rdRascalMCES.RascalOptions
    Default is opts.

Returns
-------
list[rdRascalMCES.RascalResult]
    Rascal results with similarity metrics.
    """
    return rdRascalMCES.FindMCES(mol_1, mol_2, options)
