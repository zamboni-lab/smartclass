"""Convert SMILES strings to RDKit Mol objects."""

from __future__ import annotations

from rdkit.Chem import Mol, MolFromSmiles

from smartclass.chem.helpers.check_mol import check_mol
from smartclass.logging import get_logger


__all__ = [
    "convert_smiles_to_mol",
]

logger = get_logger(__name__)


def convert_smiles_to_mol(smiles: str, sanitize: bool = True) -> Mol | None:
    """
    Convert a SMILES string to an RDKit Mol object.

    Performs sanitization and validation of the molecule. Returns None
    if the SMILES is invalid or the molecule fails validation.

    :param smiles: SMILES string to convert.
    :param sanitize: Whether to sanitize the molecule. Default True.
    :returns: RDKit Mol object, or None if conversion fails.
    """
    if not smiles or not isinstance(smiles, str):
        return None

    mol = MolFromSmiles(smiles, sanitize=sanitize)

    if mol is None:
        logger.debug(f"Failed to parse SMILES: {smiles}...")
        return None

    mol, errors = check_mol(mol)

    if errors:
        logger.debug(
            f"Molecule validation failed for SMILES {smiles[:30]}...: {errors}"
        )
        return None

    return mol


if __name__ == "__main__":
    test_smiles = ["N[C@@H](CCCNC(N)=N)C(O)=O", "invalid_smiles", ""]

    for smi in test_smiles:
        result = convert_smiles_to_mol(smi)
        print(f"{smi[:30]}: {result is not None}")
