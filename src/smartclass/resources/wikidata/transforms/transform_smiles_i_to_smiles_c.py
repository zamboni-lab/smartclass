"""Transform smiles_i to smiles_c."""

from __future__ import annotations

from smartclass.helpers.check_smiles_isomeric import check_smiles_isomeric


def transform_smiles_i_to_smiles_c(result: dict) -> dict | None:
    """
    Transform smiles_i to smiles_c.

    :param result: A single query result (dictionary).
    :type result: dict

    :returns: Transformed query result.
    :rtype: Union[dict, None]
    """
    s = result.get("smiles_isomeric", "")
    smiles = check_smiles_isomeric(s, transform_to_canonical=True)
    if smiles is not None:
        transformed_result = {
            "qid": result.get("structure", "").replace("http://www.wikidata.org/entity/", ""),
            "P233": '"' + smiles + '"',
            "S887": "Q123282952",
        }
        return transformed_result
    else:
        return None
