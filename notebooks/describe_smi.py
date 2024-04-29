"""Describe SMILES."""

from __future__ import annotations

from smartclass.chem.helpers.calculate_descriptors import calculate_descriptors
from smartclass.io.export_results import export_results
from smartclass.io.load_smiles import load_smiles


def describe_smi(input: str, output: str) -> None:
    """
    Describe SMILES.

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
        descriptors = calculate_descriptors(smiles)
        if descriptors is not None:
            result = {"smiles": smiles, **descriptors}
            results.append(result)

    export_results(output=output, results=results)


if __name__ == "__main__":
    describe_smi(input="~/Downloads/smiles_monika.csv", output="descriptors_monika.tsv")
