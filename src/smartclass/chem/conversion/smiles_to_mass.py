"""Convert a structure SMILES to an exact mass."""

from __future__ import annotations

from rdkit.Chem.Descriptors import ExactMolWt

from smartclass.chem.conversion.smiles_to_mol import smiles_to_mol

__all__ = [
    "smiles_to_mass",
]


def smiles_to_mass(smiles: str) -> float | None:
    """
    Convert a structure SMILES to an exact mass.

    :param smiles: A SMILES.
    :type smiles: str

    :returns: An exact mass.
    :rtype: Union[float, None]
    """
    mol = smiles_to_mol(smiles)
    if mol is not None:
        return ExactMolWt(mol)
    else:
        return None


if __name__ == "__main__":
    smiles_to_test = ["N[C@@H](CCCNC(N)=N)C(O)=O"]

    for smiles in smiles_to_test:
        smiles_to_mass(smiles)
