"""Transform formula to formula."""

from __future__ import annotations

from smartclass.chem.conversion.smiles_to_formula import smiles_to_formula
from smartclass.helpers.convert_chemical_formula import convert_chemical_formula


def transform_formula_to_formula(result: dict) -> dict | None:
    """
    Transform formula to formula.

    :param result: A single query result (dictionary).
    :type result: dict

    :returns: Transformed query result.
    :rtype: Union[dict, None]
    """
    formula = result.get("formula", "")
    nf = smiles_to_formula(result.get("smiles", ""))

    if nf is not None:
        new_formula = convert_chemical_formula(nf)
    else:
        new_formula = None

    if formula is not None and new_formula is not None:
        transformed_result = {
            "qid": result.get("structure", "").replace("http://www.wikidata.org/entity/", ""),
            "-P274": '"' + formula + '"',
            "P274": '"' + new_formula + '"',
            "S887": "Q113907573",
        }
        return transformed_result
    else:
        return None
