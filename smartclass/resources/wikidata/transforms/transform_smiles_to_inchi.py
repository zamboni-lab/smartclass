"""Transform SMILES to InChI."""

from __future__ import annotations

from smartclass.chem.conversion.convert_smiles_to_inchi import convert_smiles_to_inchi


def transform_smiles_to_inchi(result: dict) -> dict | None:
    """Transform SMILES to InChI.

    Parameters
    ----------
    result : dict
        A single query result (dictionary).

    Returns
    -------
    dict | None
        Transformed query result.
    """
    inchi = convert_smiles_to_inchi(result.get("smiles", ""))
    inchi_wd = result.get("inchi")
    if inchi and len(inchi) < 1500:
        if inchi_wd is None:
            transformed_result = {
                "qid": result.get("structure", "").replace(
                    "http://www.wikidata.org/entity/",
                    "",
                ),
                "P234": '"' + inchi + '"',
                "S887": "Q113907573",
            }
        else:
            if inchi != inchi_wd:
                transformed_result = {
                    "qid": result.get("structure", "").replace(
                        "http://www.wikidata.org/entity/",
                        "",
                    ),
                    "P234": '"' + inchi + '"',
                    "S887": "Q113907573",
                    "-P234": '"' + inchi_wd + '"',
                }
            else:
                return None
        return transformed_result
    else:
        return None
