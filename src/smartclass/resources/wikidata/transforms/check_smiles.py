"""Check SMILES."""

from __future__ import annotations

from smartclass.chem.conversion.convert_smiles_to_mol import convert_smiles_to_mol


def check_smiles(result: dict) -> dict | None:
    """
    Check if the SMILES string is chemically valid (valence, charges, etc.)

    :param result: A single query result (dictionary).
    :type result: dict

    :returns: Dictionary with full validity information.
    :rtype: Union[dict, None]
    """
    smiles = result.get("smiles", "")
    qid = result.get("structure", "").replace("http://www.wikidata.org/entity/", "")

    mol = convert_smiles_to_mol(smiles)

    if mol is None:
        return {
            "qid": qid,
            "smiles": f'"{smiles}"',
        }
    else:
        return None
