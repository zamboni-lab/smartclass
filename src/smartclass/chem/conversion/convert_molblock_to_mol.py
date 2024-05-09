"""Convert a structure MOLBlock to MOL."""

from __future__ import annotations

from rdkit.Chem import Mol, MolFromMolBlock

__all__ = [
    "convert_molblock_to_mol",
]


def convert_molblock_to_mol(molblock: str) -> Mol | None:
    """
    Convert a structure MOLBlock to MOL.

    :param molblock: A MOLBlock.
    :type molblock: str

    :returns: A MOL.
    :rtype: Union[Mol, None]
    """
    return MolFromMolBlock(molblock)


# Example usage:
if __name__ == "__main__":
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
    result = convert_molblock_to_mol(structure_mol)
