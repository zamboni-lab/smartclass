"""Convert a structure SMARTS to MOL."""

from __future__ import annotations

from rdkit.Chem import Mol, MolFromSmarts, MolFromSmiles

from smartclass.chem.helpers.check_mol import check_mol


__all__ = [
    "convert_smarts_to_mol",
]


def convert_smarts_to_mol(smarts: str) -> Mol | None:
    """Convert a structure SMARTS to MOL.

    Parameters
    ----------
    smarts : str
        SMARTS.

    Returns
    -------
    Mol | None
        MOL.
    """
    # Try SMARTS first; if that fails, fall back to SMILES/CXSMILES parsing.
    mol = MolFromSmarts(smarts)
    if mol is None:
        mol = MolFromSmiles(smarts)

    if mol is None:
        return None

    mol, errors = check_mol(mol)

    if errors:
        return None

    return mol


if __name__ == "__main__":
    smarts_to_test = [
        "[#8]-[#6]-[#6@H]-1-[#8]-[#6@@H](-[#8]-[#6@@H]-2-[#8]-[#6]=[#6]-3-[#6@@H](-[#6]-[#6]-[#8]-[#6]-3=O)-[#6@H]-2-[#6]=[#6])-[#6@H](-[#8])-[#6@@H](-[#8])-[#6@@H]-1-[#8]",
    ]

    for smarts in smarts_to_test:
        convert_smarts_to_mol(smarts)
