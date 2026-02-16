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
) -> rdRascalMCES.RascalResult:
    """
    Calculate structural similarity.

    :param mol_1: Mol 1.
    :type mol_1: Mol

    :param mol_2: Mol 2.
    :type mol_2: Mol

    :param options: Options.
    :type options: rdRascalMCES.RascalOptions

    :returns: RascalResult with similarity metrics.
    :rtype: rdRascalMCES.RascalResult
    """
    return rdRascalMCES.FindMCES(mol_1, mol_2, options)
