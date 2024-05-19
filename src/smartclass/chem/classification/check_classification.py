"""TODO.

TODO.
"""

# MISSING DEF
from __future__ import annotations

import logging
import random

import polars

from smartclass.chem.conversion.convert_smiles_to_mol import convert_smiles_to_mol
from smartclass.chem.similarity.calculate_mcs import calculate_mcs

__all__ = ["sample_structures"]

logging.basicConfig(level=logging.INFO)  # Set logging level


def sample_structures(smiles_list: list, max_samples: int = 1000) -> list:
    """
    Sample structures from a list of SMILES strings.

    :param smiles_list: List of SMILES strings.
    :type smiles_list: list

    :param max_samples: Maximum number of samples to return.
    :type max_samples: int

    :returns: Sampled SMILES strings.
    :rtype: list
    """
    return random.sample(smiles_list, min(len(smiles_list), max_samples))


def process_classification(classification: str, group: polars.DataFrame):
    """
    Processes a chemical classification group.

    :param classification: Chemical classification name.
    :type classification: str

    :param group: DataFrame group containing structures.
    :type group: polars.DataFrame
    """
    smiles_list = group["smiles"].to_list()

    # Sample structures
    sampled_smiles = sample_structures(smiles_list, max_samples=1000)
    num_sampled_structures = len(sampled_smiles)

    if num_sampled_structures < 6:
        logging.warning(f"Not enough structures ({num_sampled_structures}) for {classification}.")
        return

    if num_sampled_structures > 1000:
        logging.warning(f"Too many structures ({num_sampled_structures}) for {classification}.")
        return

    # Convert SMILES to molecular objects
    mols = [
        convert_smiles_to_mol(smiles)
        for smiles in sampled_smiles
        if convert_smiles_to_mol(smiles) is not None
    ]

    if not mols:
        logging.warning(f"No valid molecular objects for {classification}.")
        return

    # Calculate MCS
    mcs = calculate_mcs(mols=mols).smartsString
    logging.info(f"Chemical Classification: {classification}")
    logging.info(f"Maximal Common Substructure: {mcs}")
    logging.info(f"Based on {len(mols)} Substructures")


def main():
    """TODO."""
    # Load data from a CSV file with "ParentName" and "SMILES" columns
    csv_file = "scratch/chebi_matched_molecules.tsv"  # Replace with the path to your CSV file
    try:
        df = polars.read_csv(csv_file, separator="\t")
    except Exception as e:
        logging.error(f"Error loading data from CSV file: {e}")
        return

    # Define the chemical classifications you want to process
    target_classifications = ["Quassinoids", "4'-hydroxyflavonoids"]

    # Iterate through the specified chemical classifications
    for classification in target_classifications:
        filtered_df = df.filter(polars.col("ParentName") == classification)
        if len(filtered_df) > 0:
            process_classification(classification, filtered_df)
        else:
            logging.warning(f"Chemical Classification {classification} not found in the data.")
    else:
        logging.warning(f"Chemical Classification {classification} not found in the data.")


if __name__ == "__main__":
    main()
