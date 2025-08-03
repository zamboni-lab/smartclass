"""Transform smiles_i to smiles_i."""

from __future__ import annotations

from smartclass.helpers.check_smiles_isomeric import check_smiles_isomeric


def transform_smiles_i_to_smiles_i(result: dict) -> dict | None:
    """
    Transform smiles_i to smiles_i.

    :param result: A single query result (dictionary).
    :type result: dict

    :returns: Transformed query result.
    :rtype: Union[dict, None]
    """
    s = result.get("smiles", "")
    smiles = check_smiles_isomeric(s)
    if smiles:
        transformed_result = {
            "qid": result.get("structure", "").replace(
                "http://www.wikidata.org/entity/", ""
            ),
            "-P2017": '"' + smiles + '"',
        }
        return transformed_result
    else:
        return None
