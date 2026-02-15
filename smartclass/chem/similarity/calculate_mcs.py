"""Calculate Maximum Common Substructure (MCS) for molecules."""

from __future__ import annotations

from rdkit.Chem import rdFMCS

from smartclass.config import get_config

__all__ = [
    "calculate_mcs",
]


def calculate_mcs(
    mols: list,
    threshold: float | None = None,
    ring_matches_ring_only: bool = False,
    timeout: int | None = None,
) -> rdFMCS.MCSResult:
    """
    Calculate Maximum Common Substructure (MCS) for a list of molecules.

    Uses RDKit's FindMCS algorithm with flexible atom and bond matching.

    :param mols: List of RDKit Mol objects.
    :param threshold: Fraction of molecules that must contain the substructure.
        Uses config default if None.
    :param ring_matches_ring_only: If True, ring atoms only match ring atoms.
    :param timeout: Maximum time in seconds for MCS calculation.
        Uses config default if None.
    :returns: MCS result object containing the common substructure.
    """
    config = get_config()

    threshold = threshold if threshold is not None else config.mcs_threshold
    timeout = timeout if timeout is not None else config.mcs_timeout

    return rdFMCS.FindMCS(
        mols=mols,
        threshold=threshold,
        ringMatchesRingOnly=ring_matches_ring_only,
        timeout=timeout,
        atomCompare=rdFMCS.AtomCompare.CompareAny,
        bondCompare=rdFMCS.BondCompare.CompareAny,
    )
