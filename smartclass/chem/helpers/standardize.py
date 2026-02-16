"""Standardize."""

from __future__ import annotations

from rdkit.Chem.MolStandardize.rdMolStandardize import StandardizeSmiles


__all__ = [
    "standardize",
]


def standardize(smiles: str) -> str | None:
    """
    Standardize.

    :param smiles: SMILES.
    :type smiles: str

    :returns: SMILES.
    :rtype: Union[str, None]
    """
    return StandardizeSmiles(smiles)


if __name__ == "__main__":
    smiles_to_test = [
        "CC(C)C1=CC[C@H]2C(=C1)CC[C@@H]3[C@@]2(CCC[C@@]3(C)C(=O)O)C",
        "O=C(O)CCC(N=C(O)C(N=C(O)C(N=C(O)C(N)CC=1C=CC=CC1)CC(C)C)CC2=CN=CN2)C(O)=NC(C(O)=NC3C(O)=NCCCC3)C(C)CC",
        "CC[C@H](C)C(N=C(O)[C@@H](CCC(=O)O)N=C(O)[C@H](Cc1cnc[nH]1)N=C(O)[C@H](CC(C)C)N=C(O)[C@H](N)Cc1ccccc1)C(O)=N[C@H]1CCCCN=C1O",
    ]

    for smiles in smiles_to_test:
        standardize(smiles)
