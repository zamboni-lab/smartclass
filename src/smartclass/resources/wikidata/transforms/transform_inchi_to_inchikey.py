"""Transform InChI to InChIKey."""

from __future__ import annotations

from smartclass.chem.conversion.inchi_to_inchikey import inchi_to_inchikey


def transform_inchi_to_inchikey(result: dict) -> dict | None:
    """
    Transform InChI to InChIKey.

    :param result: A single query result (dictionary).
    :type result: dict

    :returns: Transformed query result.
    :rtype: Union[dict, None]
    """
    inchikey = inchi_to_inchikey(result.get("inchi", ""))
    if inchikey is not None:
        transformed_result = {
            "qid": result.get("structure", "").replace("http://www.wikidata.org/entity/", ""),
            "P235": '"' + inchikey + '"',
            "S887": "Q123137214",
        }
        return transformed_result
    else:
        return None
