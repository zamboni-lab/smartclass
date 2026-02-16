"""Convert a structure InChI to mass."""

from __future__ import annotations

from smartclass.chem.conversion.convert_inchi_to_mol import convert_inchi_to_mol
from smartclass.chem.conversion.convert_mol_to_mass import convert_mol_to_mass


__all__ = [
    "convert_inchi_to_mass",
]


def convert_inchi_to_mass(inchi: str) -> float | None:
    """
    Convert a structure InChI to mass.

    :param inchi: An InChI.
    :type inchi: str

    :returns: A mass.
    :rtype: Union[float, None]
    """
    mol = convert_inchi_to_mol(inchi)
    if mol is not None:
        return convert_mol_to_mass(mol)
    else:
        return None


if __name__ == "__main__":
    inchis_to_test = [
        "InChI=1S/C16H22O9/c1-2-7-8-3-4-22-14(21)9(8)6-23-15(7)25-16-13(20)12(19)11(18)10(5-17)24-16/h2,6-8,10-13,15-20H,1,3-5H2/t7-,8+,10-,11-,12+,13-,15+,16+/m1/s1",
    ]

    for inchi in inchis_to_test:
        convert_inchi_to_mass(inchi)
