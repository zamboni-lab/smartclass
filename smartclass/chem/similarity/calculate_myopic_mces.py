"""Calculate myopic MCES."""

from __future__ import annotations

from myopic_mces.myopic_mces import MCES


def calculate_myopic_mces(
    s_1: str,
    s_2: str,
    the: int = -1,
) -> tuple[int, float, float, int]:
    """Calculate myopic MCES.

    Parameters
    ----------
    s_1 : str
        SMILES 1.
    s_2 : str
        SMILES 2.
    the : int
        Default is -1.

    Returns
    -------
    tuple[int, float, float, int]
        Tuple of (Index, Distance, Time, Distance Type). Distance Type values: - 1: Exact Distance - 2: Lower bound (distance above threshold) - 4: Lower bound (second bound used)
    """
    return MCES(smiles1=s_1, smiles2=s_2, threshold=the, i=0, solver="default")
