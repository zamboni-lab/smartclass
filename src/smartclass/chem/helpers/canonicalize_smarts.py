"""Canonicalize SMARTS."""

from __future__ import annotations
import logging
from rdcanon import canon_smarts

__all__ = [
    "canonicalize_smarts",
]


def canonicalize_smarts(smarts: str, return_mapping: bool = False) -> str | None:
    """
    Canonicalize SMARTS.

    :param smarts: SMARTS string to canonicalize.
    :type smarts: str

    :param return_mapping: Whether to return the atom mapping.
    :type return_mapping: bool

    :returns: Canonical SMARTS string, optionally with atom mapping. Returns None on failure.
    :rtype: Union[str, Tuple[str, List[int]], None]
    """
    try:
        return canon_smarts(smarts, return_mapping)
    except Exception as e:
        logging.error(f"Error canonicalizing SMARTS: {smarts}")
        return None


if __name__ == "__main__":
    test_smarts = [
        "[$([NX3H,NX4H2+]),$([NX3](C)(C)(C))]1[CX4H]([CH2][CH2][CH2]1)[CX3](=[OX1])[OX2H,OX1-,N]",
        "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]",
        "[CX3](=O)[OX1H0-,OX2H1]",
        "[CX3](=[OX1])[OX2][CX3](=[OX1])",
        "[N&H2&+0:4]-[C&H1&+0:2](-[C&H2&+0:8])-[O&H1&+0:3]",
        "INVALID_SMARTS",  # For testing error handling
    ]

    for smarts in test_smarts:
        print("Original:", smarts)
        print("Canonical:", canonicalize_smarts(smarts))
        print(
            "Canonical with mapping:", canonicalize_smarts(smarts, return_mapping=True)
        )
        print()
