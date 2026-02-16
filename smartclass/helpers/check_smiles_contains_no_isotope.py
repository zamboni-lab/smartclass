"""Checks if SMILES contains no isotope."""

from __future__ import annotations

import re


__all__ = [
    "check_smiles_contains_no_isotope",
]


def check_smiles_contains_no_isotope(
    smiles: str,
) -> bool:
    """
    Checks if SMILES contains no isotope.

    :param smiles: The SMILES.

    :type smiles: str

    :returns: A boolean if the SMILES meets the given criteria.

    :rtype: bool
    """
    return not re.search(r"[2-9]\]", smiles) and not re.search(r"\[[1-9]", smiles)
