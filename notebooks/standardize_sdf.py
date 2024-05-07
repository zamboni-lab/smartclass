"""Standardize SDF."""

from __future__ import annotations

from rdkit.Chem import MultithreadedSDMolSupplier

from smartclass.chem.conversion.convert_mol_to_smiles import convert_mol_to_smiles
from smartclass.chem.conversion.convert_smiles_to_formula import convert_smiles_to_formula
from smartclass.chem.conversion.convert_smiles_to_mass import convert_smiles_to_mass
from smartclass.chem.helpers.standardize import standardize
from smartclass.io.export_results import export_results


def standardize_sdf(input: str, output: str) -> None:
    """
    Standardize SDF.

    :param input: The path to the input SDF file.
    :type input: str

    :param output: The path to the output CSV file.
    :type output: str

    :returns: SMILES.
    :rtype: Union[str, None]
    """
    s = [
        (convert_mol_to_smiles(m), m.GetProp("CompoundName"))
        for m in MultithreadedSDMolSupplier(input)
        if m is not None
    ]

    results = []
    for smiles, CompoundName in s:
        result = {
            "name": CompoundName,
            "formula": convert_smiles_to_formula(smiles),
            "exact_mass": convert_smiles_to_mass(smiles),
            "smiles": smiles,
            "new_smiles": standardize(smiles),
        }
        results.append(result)

    export_results(output=output, results=results)


if __name__ == "__main__":
    standardize_sdf(input="../../Downloads/20230818_Pickedcompounds.sdf", output="testMario.csv")
