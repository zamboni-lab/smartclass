"""Checks if SMILES contains no dot."""

from __future__ import annotations

__all__ = [
    "check_smiles_contains_no_dot",
]


def check_smiles_contains_no_dot(
    smiles: str,
) -> bool:
    """
    Checks if SMILES contains no dot.

    :param smiles: The SMILES.

    :type smiles: str

    :returns: A boolean if the SMILES meets the given criteria.

    :rtype: bool
    """
    return "." not in smiles
