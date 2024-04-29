"""Transform SMILES to formula."""

from __future__ import annotations

from smartclass.chem.conversion.smiles_to_formula import smiles_to_formula
from smartclass.helpers.convert_chemical_formula import convert_chemical_formula


def transform_smiles_to_formula(result: dict) -> dict | None:
    """
    Transform SMILES to formula.

    :param result: A single query result (dictionary).
    :type result: dict

    :returns: Transformed query result.
    :rtype: Union[dict, None]
    """
    formula = smiles_to_formula(result.get("smiles", ""))
    if formula is not None:
        transformed_result = {
            "qid": result.get("structure", "").replace("http://www.wikidata.org/entity/", ""),
            "P274": '"' + convert_chemical_formula(formula) + '"',
            "S887": "Q113907573",
        }
        return transformed_result
    else:
        return None
