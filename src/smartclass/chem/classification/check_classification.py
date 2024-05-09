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

# Load data from a CSV file with "ParentName" and "SMILES" columns
csv_file = "matched_molecules.csv"  # Replace with the path to your CSV file
df = polars.read_csv(csv_file)


def sample_structures(smiles_list: list, max_samples: int = 1000) -> list:
    """
    Description.

    :param smiles_list: Description.
    :type smiles_list: list

    :param max_samples: Description.
    :type max_samples: int

    :returns: Description.
    :rtype: list
    """
    if len(smiles_list) > max_samples:
        return random.sample(smiles_list, max_samples)
    return smiles_list


grouped = df.group_by("ParentName")

# Define the chemical classifications you want to process
target_classifications = ["Quassinoids", "4'-hydroxyflavonoids"]

# Iterate through the specified chemical classifications
for classification in target_classifications:
    if classification in grouped.groups():
        group = grouped.get_group(classification)
        smiles_list = group["SMILES"].to_list()
        num_structures = len(smiles_list)
        sampled_smiles = sample_structures(smiles_list, max_samples=1000)
        num_sampled_structures = len(sampled_smiles)

        # Check if the group has at least 6 structures but no more than 1000
        if 6 <= num_sampled_structures <= 1000:
            # take this out and pre-compute uniquely
            mols = [
                convert_smiles_to_mol(smiles)
                for smiles in smiles_list
                if convert_smiles_to_mol(smiles) is not None
            ]
            if mols:
                # WAYYYYY TOOOOOO SLOOOOOW
                mcs = calculate_mcs(mols).smartsString
                # mcs = rdFMCS.FindMCS(mols).smartsString
                # Note: Implement the following:
                # compute MCS with no threshold,
                # match the structure against all back,
                # if specificity <1, lower the threshold iteratively
                logging.debug(f"Chemical Classification: {classification}")
                logging.debug(f"Maximal Common Substructure: {mcs}")
                # TODO be more precise about discarded ones
                logging.debug(f"Based on {len(mols)} Substructures")
        else:
            logging.debug(
                f"Chemical Classification ,{classification}, does not have the required number of structures."
            )
    else:
        logging.debug(f"Chemical Classification ,{classification}, not found in the data.")
