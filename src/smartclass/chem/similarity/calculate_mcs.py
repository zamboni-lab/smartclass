"""Calculate MCS."""

from __future__ import annotations

from rdkit.Chem import rdFMCS

__all__ = [
    "calculate_mcs",
]


def calculate_mcs(
    mols: list, threshold: float = 0.7, ring_matches_ring_only: bool = False
) -> rdFMCS.MCSResult:
    """
    Calculate MCS.

    :param mols: List of mols.
    :type mols: list

    :param threshold: Threshold.
    :type threshold: float

    :param ring_matches_ring_only: Flag to indicate if ring matches ring only.
    :type ring_matches_ring_only: bool

    :returns: A MCS result object.
    :rtype: rdFMCS.MCSResult
    """
    # Implement check to see if result ended up too early
    return rdFMCS.FindMCS(
        mols, threshold=threshold, ring_matches_ring_only=ring_matches_ring_only, timeout=60
    )
