"""Convert a structure SMILES to canonical SMILES."""

from __future__ import annotations

from rdkit.Chem import RemoveStereochemistry

from smartclass.chem.conversion.mol_to_smiles import mol_to_smiles
from smartclass.chem.conversion.smiles_to_mol import smiles_to_mol

__all__ = [
    "smiles_to_canonical_smiles",
]


def smiles_to_canonical_smiles(smiles: str) -> str | None:
    """
    Convert a structure SMILES to canonical SMILES.

    :param smiles: A SMILES.
    :type smiles: str

    :returns: A canonical SMILES.
    :rtype: Union[str, None]
    """
    mol = smiles_to_mol(smiles)
    if mol is not None:
        RemoveStereochemistry(mol)
        return mol_to_smiles(mol)
    else:
        return None


if __name__ == "__main__":
    smiles_to_test = ["N[C@@H](CCCNC(N)=N)C(O)=O"]

    for smiles in smiles_to_test:
        smiles_to_canonical_smiles(smiles)
