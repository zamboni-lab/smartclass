"""Calculate structural similarity."""

from rdkit.Chem import rdRascalMCES

opts = rdRascalMCES.RascalOptions()
opts.returnEmptyMCES = True
opts.similarityThreshold = -1
opts.timeout = 0


def calculate_structural_similarity(mol_1, mol_2, options=opts):
    """
    Calculate structural similarity.

    :param mol_1: Mol 1.
    :type mol_1: TODO

    :param mol_2: Mol 2.
    :type mol_2: TODO

    :param options: Options.
    :type options: TODO

    :returns: TODO.
    :rtype: TODO
    """
    return rdRascalMCES.FindMCES(mol_1, mol_2, options)
