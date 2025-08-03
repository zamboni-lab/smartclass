"""Transform group of stereoisomers to chemical entities."""

from __future__ import annotations

from smartclass.chem.helpers.check_missing_stereochemistry import (
    check_missing_stereochemistry,
)


def transform_stereoisomers_to_entities(result: dict) -> dict | None:
    """
    Transform group of stereoisomers to chemical entities.

    :param result: A single query result (dictionary).
    :type result: dict

    :returns: Transformed query result.
    :rtype: Union[dict, None]
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
