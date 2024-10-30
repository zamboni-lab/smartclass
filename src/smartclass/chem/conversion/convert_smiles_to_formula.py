"""Convert a structure SMILES to a molecular formula."""

from __future__ import annotations

from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from smartclass.chem.conversion.convert_mol_to_formula import convert_mol_to_formula

from smartclass.chem.conversion.convert_smiles_to_mol import convert_smiles_to_mol

__all__ = [
    "convert_smiles_to_formula",
]


def convert_smiles_to_formula(smiles: str) -> str | None:
    """
    Convert a structure SMILES to a molecular formula.

    :param smiles: A SMILES.
    :type smiles: str

    :returns: A molecular formula.
    :rtype: Union[str, None]
    """
    mol = convert_smiles_to_mol(smiles)
    if mol is not None:
        return convert_mol_to_formula(mol)
    else:
        return None


if __name__ == "__main__":
    smiles_to_test = ["N[C@@H](CCCNC(N)=N)C(O)=O"]

    for smiles in smiles_to_test:
        convert_smiles_to_formula(smiles)
