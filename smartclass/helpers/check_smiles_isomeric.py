"""Check isomeric smiles."""

from __future__ import annotations

import re

from smartclass.chem.conversion.convert_smiles_to_canonical_smiles import (
    convert_smiles_to_canonical_smiles,
)


__all__ = [
    "check_smiles_isomeric",
]


def check_smiles_isomeric(
    smiles_isomeric: str,
    transform_to_canonical: bool = False,
) -> str | None:
    """Check isomeric smiles.

    Parameters
    ----------
    smiles_isomeric : str
        SMILES.
    transform_to_canonical : bool
        False. Default is False.

    Returns
    -------
    str | None
        SMILES.
    """
    # COMMENT this could be updated live using
    # """SELECT * WHERE { wd:P2017 wdt:P1793 ?reg. }"""
    pattern = re.compile(
        r"^[A-Za-z0-9+\-\*=#$:().>\\[\]%]*([@\/]|\[\d)[A-Za-z0-9+\-\*=#$:().>@\/\\[\]%]+",
    )
    # Use re.match to check if the molecular formula matches the pattern
    if re.match(pattern, smiles_isomeric):
        if transform_to_canonical is True:
            return convert_smiles_to_canonical_smiles(smiles_isomeric)
        else:
            pass
    else:
        # Avoid unknown value
        if "http" not in smiles_isomeric:
            return smiles_isomeric
    return None


if __name__ == "__main__":
    convert_smiles_to_test = ["N[C@@H](CCCNC(N)=N)C(O)=O"]

    for smiles in convert_smiles_to_test:
        check_smiles_isomeric(smiles)
