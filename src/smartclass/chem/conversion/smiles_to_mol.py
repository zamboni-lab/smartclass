"""Convert a structure SMILES to MOL."""

from __future__ import annotations

from rdkit.Chem import MolFromSmiles

__all__ = [
    "smiles_to_mol",
]


def smiles_to_mol(smiles: str) -> str | None:
    """
    Convert a structure SMILES to MOL.

    :param smiles: A SMILES.
    :type smiles: str

    :returns: A MOL.
    :rtype: Union[str, None]
    """
    return MolFromSmiles(smiles)


if __name__ == "__main__":
    smiles_to_test = ["N[C@@H](CCCNC(N)=N)C(O)=O"]

    for smiles in smiles_to_test:
        smiles_to_mol(smiles)
