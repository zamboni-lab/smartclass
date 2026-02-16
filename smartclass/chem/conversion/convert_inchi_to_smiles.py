"""Convert a structure InChI to SMILES."""

from __future__ import annotations

from smartclass.chem.conversion.convert_inchi_to_mol import convert_inchi_to_mol
from smartclass.chem.conversion.convert_mol_to_smiles import convert_mol_to_smiles


__all__ = [
    "convert_inchi_to_smiles",
]


def convert_inchi_to_smiles(inchi: str) -> str | None:
    """
    Convert a structure InChI to SMILES.

    :param inchi: An InChI.
    :type inchi: str

    :returns: A SMILES.
    :rtype: Union[str, None]
    """
    mol = convert_inchi_to_mol(inchi)
    if mol is not None:
        return convert_mol_to_smiles(mol)
    else:
        return None


if __name__ == "__main__":
    inchis_to_test = ["InChI=1S/C3H8/c1-3-2/h3H2,1-2H3"]

    for inchi in inchis_to_test:
        convert_inchi_to_smiles(inchi)
