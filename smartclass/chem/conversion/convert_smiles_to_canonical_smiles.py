"""Convert a structure SMILES to canonical SMILES."""

from __future__ import annotations

from rdkit.Chem import RemoveStereochemistry

from smartclass.chem.conversion.convert_mol_to_smiles import convert_mol_to_smiles
from smartclass.chem.conversion.convert_smiles_to_mol import convert_smiles_to_mol


__all__ = [
    "convert_smiles_to_canonical_smiles",
]


def convert_smiles_to_canonical_smiles(smiles: str) -> str | None:
    """
    Convert a structure SMILES to canonical SMILES.

    :param smiles: A SMILES.
    :type smiles: str

    :returns: A canonical SMILES.
    :rtype: Union[str, None]
    """
    mol = convert_smiles_to_mol(smiles)
    if mol is not None:
        RemoveStereochemistry(mol)
        return convert_mol_to_smiles(mol)
    else:
        return None


if __name__ == "__main__":
    smiles_to_test = ["N[C@@H](CCCNC(N)=N)C(O)=O"]

    for smiles in smiles_to_test:
        convert_smiles_to_canonical_smiles(smiles)
