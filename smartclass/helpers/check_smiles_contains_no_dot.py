"""Checks if SMILES contains no dot."""

from __future__ import annotations


__all__ = [
    "check_smiles_contains_no_dot",
]


def check_smiles_contains_no_dot(
    smiles: str,
) -> bool:
    """Checks if SMILES contains no dot.

Parameters
----------
smiles : str
    SMILES.

Returns
-------
bool
    SMILES meets the given criteria.
    """
    return "." not in smiles
