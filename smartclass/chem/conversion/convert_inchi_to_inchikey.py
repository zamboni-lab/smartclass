"""Convert a structure InChI to InChIKey."""

from __future__ import annotations

from rdkit.Chem.inchi import InchiToInchiKey


__all__ = [
    "convert_inchi_to_inchikey",
]


def convert_inchi_to_inchikey(inchi: str) -> str:
    """
    Convert a structure InChI to InChIKey.

    :param inchi: An InChI.
    :type inchi: str

    :returns: An InChIKey.
    :rtype: str
    """
    return InchiToInchiKey(inchi)


if __name__ == "__main__":
    inchis_to_test = ["InChI=1S/C3H8/c1-3-2/h3H2,1-2H3"]

    for inchi in inchis_to_test:
        convert_inchi_to_inchikey(inchi)
