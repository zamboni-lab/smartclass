"""Convert a structure InChI to mass."""

from __future__ import annotations

from rdkit.Chem.Descriptors import ExactMolWt

from smartclass.chem.conversion.convert_inchi_to_mol import convert_inchi_to_mol

__all__ = [
    "convert_inchi_to_mass",
]


def convert_inchi_to_mass(inchi: str) -> str | None:
    """
    Convert a structure InChI to mass.

    :param inchi: An InChI.
    :type inchi: str

    :returns: A mass.
    :rtype: Union[str, None]
    """
    mol = convert_inchi_to_mol(inchi)
    if mol is not None:
        return ExactMolWt(mol)
    else:
        return None


if __name__ == "__main__":
    inchis_to_test = [
        "InChI=1S/C16H22O9/c1-2-7-8-3-4-22-14(21)9(8)6-23-15(7)25-16-13(20)12(19)11(18)10(5-17)24-16/h2,6-8,10-13,15-20H,1,3-5H2/t7-,8+,10-,11-,12+,13-,15+,16+/m1/s1",  # noqa:E501
    ]

    for inchi in inchis_to_test:
        convert_inchi_to_mass(inchi)
