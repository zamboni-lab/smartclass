"""Calculate myopic MCES."""

from myopic_mces.myopic_mces import MCES


def calculate_myopic_mces(s_1: str, s_2: str, thr: int = -1):
    """
    Calculate myopic MCES.

    :param s_1: SMILES 1.
    :type s_1: str

    :param s_2: SMILES 2.
    :type s_2: str

    :param thr: Threshold.
    :type thr: int

    :returns: TODO.
    :rtype: TODO
    """
    return MCES(ind=0, s1=s_1, s2=s_2, threshold=thr, solver="default")
