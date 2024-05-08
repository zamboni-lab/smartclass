"""Convert a structure SMARTS to MOL."""

from __future__ import annotations

from rdkit.Chem import Mol, MolFromSmarts

__all__ = [
    "convert_smarts_to_mol",
]


def convert_smarts_to_mol(smarts: str) -> Mol | None:
    """
    Convert a structure SMARTS to MOL.

    :param smarts: A SMARTS.
    :type smarts: str

    :returns: A MOL.
    :rtype: Union[Mol, None]
    """
    return MolFromSmarts(smarts)


if __name__ == "__main__":
    smarts_to_test = ["TODO"]

    for smarts in smarts_to_test:
        convert_smarts_to_mol(smarts)
