"""Transform InChI to mass."""

from __future__ import annotations

from smartclass.chem.conversion.convert_inchi_to_mass import convert_inchi_to_mass


def transform_inchi_mass_to_inchi_mass(result: dict, tol: float = 0.001) -> dict | None:
    """Transform InChI to mass.

Parameters
----------
result : dict
    A single query result (dictionary).
tol : float
    Tolerance (in dalton). Default is 0.001.

Returns
-------
dict | None
    Transformed query result.
    """
    mass = convert_inchi_to_mass(result.get("inchi", ""))
    mass_wd = float(result.get("mass", "").split("±")[0])
    if mass:
        if abs(float(mass) - mass_wd) > tol:
            transformed_result = {
                "qid": result.get("structure", "").replace(
                    "http://www.wikidata.org/entity/",
                    "",
                ),
                "-P2067": "+" + str(mass_wd) + "U483261",
                "P2067": "+" + str(mass) + "U483261",
                "S887": "Q123137214",
            }
            return transformed_result
        else:
            return None
    else:
        return None
