"""Transform smiles_c to smiles_c_tauto."""

from __future__ import annotations

from smartclass.chem.helpers.fix_inchi_tautomerization import fix_inchi_tautomerization


def transform_smiles_c_to_smiles_c_tauto(result: dict) -> dict | None:
    """
    Transform smiles_c to smiles_c_tauto.

    :param result: A single query result (dictionary).
    :type result: dict

    :returns: Transformed query result.
    :rtype: Union[dict, None]
    """
    s = result.get("smiles", "")
    smiles = fix_inchi_tautomerization(s)
    if smiles is not s:
        transformed_result = {
            "qid": result.get("structure", "").replace("http://www.wikidata.org/entity/", ""),
            "-P233": '"' + s + '"',
            "P233": '"' + smiles + '"',
            # "S887": "Q123282952",
        }
        return transformed_result
    else:
        return None
