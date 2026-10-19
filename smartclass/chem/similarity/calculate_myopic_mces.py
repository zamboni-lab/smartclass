"""Calculate myopic MCES."""

from __future__ import annotations

from typing import cast

from myopic_mces.myopic_mces import MCES


def calculate_myopic_mces(
    s_1: str,
    s_2: str,
    threshold: int = -1,
) -> tuple[int, float, float, int]:
    """Calculate myopic MCES.

    Parameters
    ----------
    s_1 : str
        SMILES 1.
    s_2 : str
        SMILES 2.
    threshold : int
        Threshold for the MCES calculation. Default is -1.

    Returns
    -------
    tuple[int, float, float, int]
        Tuple of (Index, Distance, Time, Distance Type). Distance Type values:
        - 1: Exact Distance
        - 2: Lower bound (distance above threshold)
        - 4: Lower bound (second bound used)
    """
    return cast(
        "tuple[int, float, float, int]",
        MCES(smiles1=s_1, smiles2=s_2, threshold=threshold, i=0, solver="default"),
    )
