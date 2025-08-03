"""Calculate myopic MCES."""

from __future__ import annotations

from myopic_mces.myopic_mces import MCES


def calculate_myopic_mces(
    s_1: str,
    s_2: str,
    the: int = -1,
) -> tuple[int, float, float, int]:
    """
    Calculate myopic MCES.

    :param s_1: SMILES 1.
    :type s_1: str

    :param s_2: SMILES 2.
    :type s_2: str

    :param the: Threshold.
    :type the: int

    :returns: Index, Distance between the molecules, Time taken for the calculation, Type of Distance:
            - 1 : Exact Distance
            - 2 : Lower bound (if the exact distance is above the threshold; bound chosen dynamically)
            - 4 : Lower bound (second lower bound was used)
    :rtype: int, float, float, int
    """
    return MCES(ind=0, s1=s_1, s2=s_2, threshold=the, solver="default")
