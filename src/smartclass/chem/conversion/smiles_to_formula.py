"""Convert a structure SMILES to InChI."""

from __future__ import annotations

from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from smartclass.chem.conversion.smiles_to_mol import smiles_to_mol

__all__ = [
    "smiles_to_formula",
]


def smiles_to_formula(smiles: str) -> str | None:
    """
    Convert a structure SMILES to InChI.

    :param smiles: A SMILES.
    :type smiles: str

    :returns: An InChI.
    :rtype: Union[str, None]
    """
    mol = smiles_to_mol(smiles)
    if mol is not None:
        return CalcMolFormula(mol)
    else:
        return None


if __name__ == "__main__":
    smiles_to_test = ["N[C@@H](CCCNC(N)=N)C(O)=O"]

    for smiles in smiles_to_test:
        smiles_to_formula(smiles)
