"""Check and validate chemical classification rules.

This module provides utilities for verifying that chemical class definitions
(SMARTS patterns) correctly capture the expected structural features.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from polars import DataFrame, col

from smartclass.chem.conversion.convert_smiles_to_mol import convert_smiles_to_mol
from smartclass.chem.similarity.calculate_mcs import calculate_mcs
from smartclass.config import get_config
from smartclass.exceptions import DataLoadingError
from smartclass.helpers.sample_list import sample_list
from smartclass.io.load_tsv_from_path import load_tsv_from_path
from smartclass.logging import get_logger

if TYPE_CHECKING:
    from rdkit.Chem import rdFMCS

__all__ = ["check_classification", "process_classification"]

logger = get_logger(__name__)


def process_classification(
    classification: str,
    group: DataFrame,
    max_samples: int = 1000,
    threshold: float | None = None,
) -> tuple[rdFMCS.MCSResult, int]:
    """
    Process a chemical classification group to find common substructure.

    :param classification: Chemical classification name.
    :param group: DataFrame group containing structures with 'smiles' column.
    :param max_samples: Maximum number of molecules to sample. Default 1000.
    :param threshold: MCS threshold. Uses config default if None.
    :returns: Tuple of (MCS result, number of molecules processed).
    """
    config = get_config()
    threshold = threshold if threshold is not None else config.mcs_threshold

    logger.info(f"Processing classification: {classification}")
    smiles_list = group["smiles"].to_list()

    mols = [
        mol
        for smiles in smiles_list
        if (mol := convert_smiles_to_mol(smiles)) is not None
    ]

    if len(mols) > max_samples:
        logger.warning(
            f"Too many structures ({len(mols)}) for {classification}. "
            f"Sampling {max_samples}..."
        )
        mols = sample_list(mols, max_samples=max_samples)

    mcs = calculate_mcs(mols=mols, threshold=threshold)
    return (mcs, len(mols))


def check_classification(
    classes: tuple[str, ...] = ("Quassinoids", "4'-hydroxyflavonoids"),
    classified_mols: str | Path = "scratch/chebi_matched_molecules.tsv",
    threshold: float | None = None,
) -> list[dict]:
    """
    Check classification rules by finding MCS for each class.

    Analyzes molecules belonging to each specified class and finds their
    maximum common substructure, which can be used to validate or generate
    SMARTS patterns for classification.

    :param classes: Tuple of class names to check.
    :param classified_mols: Path to TSV file with classified molecules.
        Must have 'ParentName' and 'smiles' columns.
    :param threshold: MCS threshold. Uses config default if None.
    :returns: List of dictionaries with MCS results for each class.
    :raises DataLoadingError: If the classified molecules file cannot be loaded.
    """
    config = get_config()
    threshold = threshold if threshold is not None else config.mcs_threshold

    try:
        classifications = DataFrame(load_tsv_from_path(str(classified_mols)))
    except Exception as e:
        raise DataLoadingError(str(classified_mols), reason=str(e)) from e

    results = []
    for class_name in classes:
        group = classifications.filter(col("ParentName") == class_name)
        if len(group) > 0:
            mcs, n_mols = process_classification(class_name, group, threshold=threshold)
            results.append({
                "class": class_name,
                "class_structure": mcs.smartsString,
                "structure_ab": mcs.numAtoms + mcs.numBonds,
                "threshold": threshold,
                "n": n_mols,
            })
        else:
            logger.warning(f"Classification '{class_name}' not found in data.")
            results.append({
                "class": class_name,
                "class_structure": "",
                "structure_ab": 0,
                "threshold": threshold,
                "n": 0,
            })
    return results


if __name__ == "__main__":
    from smartclass.io.export_results import export_results

    logger.info("Running classification check (this may be slow)...")
    results = check_classification()
    logger.info(f"Results: {results}")
    export_results(output="scratch/classes_checked.tsv", results=results)
