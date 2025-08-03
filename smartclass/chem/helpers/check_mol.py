"""Check MOL."""

from __future__ import annotations

from rdkit.Chem import KekulizeException, Mol, SanitizeMol

__all__ = [
    "check_mol",
]


def check_mol(mol: Mol) -> tuple[Mol, list[str]]:
    """
    Perform chemical validity checks on an RDKit Mol object.

    :param mol: Molecule object.
    :type mol: rdkit.Chem.rdchem.Mol

    :returns: Tuple of (mol, list of error strings if any).
    :rtype: Tuple[Mol, List[str]]
    """
    errors = []

    try:
        SanitizeMol(mol)
    except KekulizeException as e:
        errors.append(f"Kekulization error: {e}")
    except ValueError as e:
        errors.append(f"Sanitization error: {e}")
    except Exception as e:
        errors.append(f"Unknown sanitization error: {e}")

    return mol, errors
