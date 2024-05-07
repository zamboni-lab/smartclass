"""Transform InChI to mass."""

from __future__ import annotations

from smartclass.chem.conversion.convert_inchi_to_mass import convert_inchi_to_mass


def transform_inchi_to_mass(result: dict) -> dict | None:
    """
    Transform InChI to mass.

    :param result: A single query result (dictionary).
    :type result: dict

    :returns: Transformed query result.
    :rtype: Union[dict, None]
    """
    mass = convert_inchi_to_mass(result.get("inchi", ""))
    if mass is not None:
        transformed_result = {
            "qid": result.get("structure", "").replace("http://www.wikidata.org/entity/", ""),
            "P2067": "+" + str(mass) + "U483261",
            "S887": "Q123137214",
        }
        return transformed_result
    else:
        return None
