"""Convert a structure InChI to MOL."""

from __future__ import annotations

from rdkit.Chem import Mol, MolFromInchi

__all__ = [
    "convert_inchi_to_mol",
]


def convert_inchi_to_mol(inchi: str) -> Mol | None:
    """
    Convert a structure InChI to MOL.

    :param inchi: An InChI.
    :type inchi: str

    :returns: A MOL.
    :rtype: Union[Mol, None]
    """
    return MolFromInchi(inchi)


if __name__ == "__main__":
    inchis_to_test = ["InChI=1S/C3H8/c1-3-2/h3H2,1-2H3"]

    for inchi in inchis_to_test:
        convert_inchi_to_mol(inchi)
