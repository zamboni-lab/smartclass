"""Convert a structure SMILES to InChI."""

from __future__ import annotations

from smartclass.chem.conversion.convert_mol_to_inchi import convert_mol_to_inchi
from smartclass.chem.conversion.convert_smiles_to_mol import convert_smiles_to_mol

convert_smiles_to_mol

__all__ = [
    "convert_smiles_to_inchi",
]


def convert_smiles_to_inchi(smiles: str) -> str | None:
    """
    Convert a structure SMILES to InChI.

    :param smiles: A SMILES.
    :type smiles: str

    :returns: An InChI.
    :rtype: Union[str, None]
    """
    mol = convert_smiles_to_mol(smiles)
    if mol is not None:
        return convert_mol_to_inchi(mol)
    else:
        return None


if __name__ == "__main__":
    smiles_to_test = ["N[C@@H](CCCNC(N)=N)C(O)=O"]

    for smiles in smiles_to_test:
        convert_smiles_to_inchi(smiles)
