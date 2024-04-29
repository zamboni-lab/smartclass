"""Convert a structure MOL to CXSMILES."""

from __future__ import annotations

import logging

from rdkit import Chem

__all__ = [
    "mol_to_cxsmiles",
]


def mol_to_cxsmiles(args: tuple[str, str]) -> tuple[str, str]:
    """
    Convert a structure MOL to CXSMILES.

    :param args: A tuple containing the structure name and structure MOL.
    :type args: tuple[str, str]

    :returns: A tuple containing the structure name and its CXSMILES representation.
    :rtype: tuple[str, str]
    """
    # Extract the arguments from the tuple
    structure_name, structure_mol = args

    # Convert the structure MOL to a molecule
    structure = Chem.MolFromMolBlock(structure_mol)

    # Check if structure parsing was successful
    if structure is not None:
        # Convert the structure to CXSMILES
        cxsmiles = Chem.MolToCXSmiles(structure)
    else:
        logging.warning(f"Warning: {structure_name} could not be parsed.")
        cxsmiles = None

    # Return the structure name and its CXSMILES representation
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
    result = mol_to_cxsmiles((structure_name, structure_mol))
