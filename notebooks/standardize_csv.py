"""Standardize CSV."""

from __future__ import annotations

from smartclass.chem.helpers.standardize import standardize
from smartclass.io.export_results import export_results
from smartclass.io.load_smiles import load_smiles


def standardize_csv(input: str, output: str) -> None:
    """
    Standardize CSV.

    :param input: The path to the input CSV file.
    :type input: str

    :param output: The path to the output CSV file.
    :type output: str

    :returns: SMILES.
    :rtype: Union[str, None]
    """
    s = load_smiles(input=input)

    results = []
    for smiles in s:
        result = {"smiles": smiles, "new_smiles": standardize(smiles)}
        results.append(result)

    export_results(output=output, results=results)


if __name__ == "__main__":
    standardize_csv(input="~/Downloads/nexus_smiles_exact_mass.csv", output="testMario.csv")
