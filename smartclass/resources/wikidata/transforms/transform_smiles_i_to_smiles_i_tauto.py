"""Transform smiles_i to smiles_i_tauto."""

from __future__ import annotations

from smartclass.chem.helpers.fix_inchi_tautomerization import fix_inchi_tautomerization


def transform_smiles_i_to_smiles_i_tauto(result: dict) -> dict | None:
    """
    Transform smiles_i to smiles_i_tauto.

    :param result: A single query result (dictionary).
    :type result: dict

    :returns: Transformed query result.
    :rtype: Union[dict, None]
    """
    s = result.get("smiles", "")
    smiles = fix_inchi_tautomerization(s)
    if smiles is not s and smiles:
        transformed_result = {
            "qid": result.get("structure", "").replace("http://www.wikidata.org/entity/", ""),
            "-P2017": '"' + s + '"',
            "P2017": '"' + smiles + '"',
            # "S887": "Q123137214",
        }
        return transformed_result
    else:
        return None
