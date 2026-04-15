"""Convert RDKit Mol objects to InChIKey."""

from __future__ import annotations

from rdkit.Chem import Mol, MolToInchiKey


__all__ = [
    "convert_mol_to_inchikey",
]


def convert_mol_to_inchikey(mol: Mol) -> str:
    """Convert an RDKit Mol object to an InChIKey.

    Parameters
    ----------
    mol : Mol
        Mol object.

    Returns
    -------
    str
        InChIKey string.
    """
    return MolToInchiKey(mol)
