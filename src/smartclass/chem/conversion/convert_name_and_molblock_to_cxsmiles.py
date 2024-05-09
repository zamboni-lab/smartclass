"""Convert a name and a structure MOLBLOCK to CXSMILES."""

from __future__ import annotations

import logging

from smartclass.chem.conversion.convert_mol_to_cxsmiles import convert_mol_to_cxsmiles
from smartclass.chem.conversion.convert_molblock_to_mol import convert_molblock_to_mol

__all__ = [
    "convert_name_and_molblock_to_cxsmiles",
]


def convert_name_and_molblock_to_cxsmiles(args: tuple[str, str]) -> tuple[str, str | None]:
    """
    Convert a name and a structure MOLBLOCK to CXSMILES.

    :param args: A tuple containing the structure name and structure MOL.
    :type args: tuple[str, str]

    :returns: A tuple containing the structure name and its CXSMILES representation.
    :rtype: tuple[str, Union[str,None]]
    """
    structure_name, structure_mol = args

    structure = convert_molblock_to_mol(structure_mol)

    # Check if structure parsing was successful
    if structure is not None:
        cxsmiles = convert_mol_to_cxsmiles(structure)
    else:
        logging.warning(f"Warning: {structure_name} could not be parsed.")
        cxsmiles = None

    return structure_name, cxsmiles


# Example usage:
if __name__ == "__main__":
    structure_name = "Example structure"
    structure_mol = """
    ExampleMolecule
    Marvin  11280415412D

  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0
    0.5000    0.8660    0.0000 H   0  0  0  0  0  0
    0.5000   -0.8660    0.0000 H   0  0  0  0  0  0
M  END
"""
    result = convert_name_and_molblock_to_cxsmiles((structure_name, structure_mol))
