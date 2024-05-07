"""Calculate descriptors."""

from __future__ import annotations

from rdkit.Chem.Descriptors import CalcMolDescriptors

from smartclass.chem.conversion.convert_smiles_to_mol import convert_smiles_to_mol

__all__ = [
    "calculate_descriptors",
]


def calculate_descriptors(smiles: str) -> str | None:
    """
    Calculate descriptors.

    :param smiles: SMILES.
    :type smiles: str

    :returns: SMILES.
    :rtype: Union[str, None]
    """
    mol = convert_smiles_to_mol(smiles)

    if mol is not None:
        return CalcMolDescriptors(mol)
    else:
        return None


if __name__ == "__main__":
    convert_smiles_to_test = [
        "CC(C)C1=CC[C@H]2C(=C1)CC[C@@H]3[C@@]2(CCC[C@@]3(C)C(=O)O)C",
        "O=C(O)CCC(N=C(O)C(N=C(O)C(N=C(O)C(N)CC=1C=CC=CC1)CC(C)C)CC2=CN=CN2)C(O)=NC(C(O)=NC3C(O)=NCCCC3)C(C)CC",
        "CC[C@H](C)C(N=C(O)[C@@H](CCC(=O)O)N=C(O)[C@H](Cc1cnc[nH]1)N=C(O)[C@H](CC(C)C)N=C(O)[C@H](N)Cc1ccccc1)C(O)=N[C@H]1CCCCN=C1O",  # noqa:E501
    ]

    for smiles in convert_smiles_to_test:
        calculate_descriptors(smiles)
