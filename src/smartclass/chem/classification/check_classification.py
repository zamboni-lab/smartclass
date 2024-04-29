"""TODO.

TODO.
"""

# MISSING DEF
from __future__ import annotations

import logging
import random

import polars
from rdkit import Chem
from rdkit.Chem import rdFMCS

__all__ = ["compute_mcs", "sample_structures"]

# Load data from a CSV file with "ParentName" and "SMILES" columns
csv_file = "matched_molecules.csv"  # Replace with the path to your CSV file
df = polars.read_csv(csv_file)


def compute_mcs(
    mols: Chem.Mol, threshold: float = 0.7, ring_matches_ring_only: bool = False
) -> rdFMCS.MCSResult:
    """
    Description.

    :param mols: Description.
    :type mols: Chem.Mol

    :param threshold: Description.
    :type threshold: float

    :param ring_matches_ring_only: Description.
    :type ring_matches_ring_only: bool

    :returns: Description.
    :rtype: rdFMCS.MCSResult
    """
    # Implement check to see if result ended up too early
    return rdFMCS.FindMCS(
        mols, threshold=threshold, ring_matches_ring_only=ring_matches_ring_only, timeout=60
    )


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
                Chem.MolFromSmiles(smiles)
                for smiles in smiles_list
                if Chem.MolFromSmiles(smiles) is not None
            ]
            if mols:
                # WAYYYYY TOOOOOO SLOOOOOW
                mcs = compute_mcs(mols).smartsString
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
