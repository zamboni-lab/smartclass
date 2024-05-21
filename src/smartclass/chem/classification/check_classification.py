"""Check classification."""

from __future__ import annotations

import logging

from polars import DataFrame, col

from smartclass.chem.conversion.convert_smiles_to_mol import convert_smiles_to_mol
from smartclass.chem.similarity.calculate_mcs import calculate_mcs
from smartclass.helpers.sample_list import sample_list
from smartclass.io.load_tsv_from_path import load_tsv_from_path

__all__ = ["check_classification"]


def process_classification(
    classification: str,
    group: DataFrame,
    samples_min: int = 6,
    samples_max: int = 1000,
    threshold: float = 0.7,
) -> tuple:
    """
    Processes a chemical classification group.

    :param classification: Chemical classification name.
    :type classification: str

    :param group: DataFrame group containing structures.
    :type group: DataFrame

    :param samples_min: Minimum number of samples to check. Default to 1.
    :type samples_min: int

    :param samples_max: Minimum number of samples to check. Default to 1000.
    :type samples_max: int

    :param threshold: Threshold. Default to 0.7.
    :type threshold: float

    :returns: A SMARTS.
    :rtype: tuple
    """
    logging.info(f"Processing {classification}")
    smiles_list = group["smiles"].to_list()

    mols = [mol for smiles in smiles_list if (mol := convert_smiles_to_mol(smiles)) is not None]

    if len(mols) > samples_max:
        logging.warning(f"Too many structures of this class. Sampling {samples_max}...")
        mols = sample_list(mols, samples_max=samples_max)

    mcs = calculate_mcs(mols=mols, threshold=threshold)

    return (mcs, len(mols))


def check_classification(
    classes: tuple = ("Quassinoids", "4'-hydroxyflavonoids"),
    classified_mols: str = "scratch/chebi_matched_molecules.tsv",
    threshold: float = 0.7,
) -> list[dict]:
    """
    Check classification.

    :param classes: Tuple of classes to check.
    :type classes: tuple

    :param classified_mols: TSV file containing classified molecules.
    :type classified_mols: str

    :param threshold: Threshold. Default to 0.7.
    :type threshold: float

    :returns: A list (of dictionaries).
    :rtype: list[dict]
    """
    try:
        classifications = DataFrame(load_tsv_from_path(classified_mols))
    except Exception as e:
        logging.error(f"Error loading classified mols: {e}")
        return []

    mcses = []
    # Iterate through the specified chemical classifications
    for c in classes:
        group = classifications.filter(col("ParentName") == c)
        if len(group) > 0:
            results = process_classification(c, group)
            mcs = results[0]
            n = results[1]
            mcses.append(
                {
                    "class": c,
                    "class_structure": mcs.smartsString,
                    "structure_ab": mcs.numAtoms + mcs.numBonds,
                    "threshold": threshold,
                    "n": n,
                }
            )
        else:
            logging.warning(f"Chemical Classification {c} not found in the data.")
            mcses.append(
                {
                    "class": c,
                    "class_structure": "",
                    "structure_ab": 0,
                    "threshold": threshold,
                    "n": 0,
                }
            )
    return mcses


if __name__ == "__main__":
    from smartclass.io.export_results import export_results

    logging.basicConfig(level=logging.INFO)  # Set logging level
    logging.info("This is slow for now...")
    results = check_classification()
    logging.info(results)
    export_results(output="scratch/classes_checked.tsv", results=results)
