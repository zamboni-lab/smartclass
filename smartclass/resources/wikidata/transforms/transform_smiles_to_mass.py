"""Transform SMILES to mass."""

from __future__ import annotations

from smartclass.chem.conversion.convert_smiles_to_mass import convert_smiles_to_mass


def transform_smiles_to_mass(result: dict) -> dict | None:
    """
    Transform SMILES to mass.

    :param result: A single query result (dictionary).
    :type result: dict

    :returns: Transformed query result.
    :rtype: Union[dict, None]
    """
    mass = convert_smiles_to_mass(result.get("smiles", ""))
    if mass:
        transformed_result = {
            "qid": result.get("structure", "").replace(
                "http://www.wikidata.org/entity/",
                "",
            ),
            "P2067": "+" + str(mass) + "U483261",
            "S887": "Q113907573",
        }
        return transformed_result
    else:
        return None
