"""Transform group of stereoisomers to chemical entities."""

from __future__ import annotations

from smartclass.chem.helpers.check_missing_stereochemistry import (
    check_missing_stereochemistry,
)


def transform_stereoisomers_to_entities(result: dict) -> dict | None:
    """Transform group of stereoisomers to chemical entities.

Parameters
----------
result : dict
    A single query result (dictionary).

Returns
-------
dict | None
    Transformed query result.
    """
    flag = check_missing_stereochemistry(result.get("smiles", ""))
    if flag is False:
        transformed_result = {
            "qid": result.get("structure", "").replace(
                "http://www.wikidata.org/entity/",
                "",
            ),
            "P31": "Q113145171",
            "-P31": "Q59199015",
        }
        return transformed_result
    else:
        return None
