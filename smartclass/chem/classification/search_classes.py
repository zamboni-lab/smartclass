"""Substructure search for chemical classes.

This module provides the main classification functionality for smartclass,
matching chemical structures against SMARTS-based chemical class definitions.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from rdkit.Chem import SubstructMatchParameters

from smartclass.chem.classification.bfs_search_classes_generator import (
    tqdm_bfs_search_classes_generator,
)
from smartclass.chem.conversion.convert_smiles_to_mol import convert_smiles_to_mol
from smartclass.config import get_config
from smartclass.helpers import convert_list_of_dict
from smartclass.io import (
    export_results,
    load_external_classes_file,
    load_pkg_chemical_hierarchy,
    load_pkg_classes,
    load_smiles,
)
from smartclass.logging import get_logger
from smartclass.resources.chembl import load_latest_chembl

if TYPE_CHECKING:
    from rdkit.Chem import Mol

__all__ = [
    "search_classes",
    "ClassificationResult",
]

logger = get_logger(__name__)

# Result field names for export
RESULT_FIELDS = [
    "class_id",
    "class_structure",
    "structure_inchikey",
    "structure_smarts",
    "structure_ab",
    "matched_ab",
]


def _parse_smiles_to_mols(smiles_set: set[str]) -> list[Mol]:
    """Convert a set of SMILES strings to RDKit Mol objects.

    :param smiles_set: Set of SMILES strings to convert.
    :returns: List of valid Mol objects (invalid SMILES are filtered out).
    """
    structures = []
    for smi in smiles_set:
        # Skip header rows that might be in the data
        if smi.lower() == "smiles":
            continue
        mol = convert_smiles_to_mol(smi)
        if mol is not None:
            structures.append(mol)
        else:
            logger.debug(f"Failed to parse SMILES: {smi[:50]}...")
    return structures


def _load_classes(
    classes_file: str | None,
    classes_name_id: str | None,
    classes_name_smarts: str | None,
) -> list[dict[str, list[str]]]:
    """Load chemical classes from file or package defaults.

    :param classes_file: Optional path to external classes file.
    :param classes_name_id: Name of ID column in external file.
    :param classes_name_smarts: Name of SMARTS column in external file.
    :returns: List containing a dictionary mapping class IDs to SMARTS patterns.
    """
    if classes_file:
        c = load_external_classes_file(
            file=classes_file,
            id_name=classes_name_id,
            smarts_name=classes_name_smarts,
        )
    else:
        logger.info("No classes file provided, loading default package classes.")
        c = load_pkg_classes()

    # Build dictionary mapping class_id to list of SMARTS patterns
    classes_dict: dict[str, list[str]] = {}
    for row in c.iter_rows():
        key = row[0]
        value = row[1]
        if key in classes_dict:
            classes_dict[key].append(value)
        else:
            classes_dict[key] = [value]

    return [classes_dict]


def _filter_closest_matches(results: list[dict]) -> list[dict]:
    """Filter results to keep only the closest class match for each structure.

    :param results: List of classification results.
    :returns: Filtered list with only closest matches per InChIKey.
    """
    max_ab_per_inchikey: dict[str, int] = {}
    for result in results:
        inchikey = result["structure_inchikey"]
        matched_ab = result["matched_ab"]
        if (
            inchikey not in max_ab_per_inchikey
            or matched_ab > max_ab_per_inchikey[inchikey]
        ):
            max_ab_per_inchikey[inchikey] = matched_ab

    return [
        result
        for result in results
        if result["matched_ab"] == max_ab_per_inchikey[result["structure_inchikey"]]
    ]


def _export_classification_results(
    results: list[dict],
    output_dir: Path | None = None,
) -> None:
    """Export classification results to multiple formats.

    :param results: List of classification results to export.
    :param output_dir: Directory for output files. Uses config default if None.
    """
    if not results:
        logger.warning("No results to export.")
        return

    config = get_config()
    output_path = output_dir or config.output_dir

    # Export sorted by matched_ab (descending)
    results_by_match = sorted(results, key=lambda x: x["matched_ab"], reverse=True)
    export_results(
        output=str(output_path / "results_by_class.tsv"),
        results=results_by_match,
    )

    # Export sorted by structure
    results_by_structure = sorted(
        results_by_match, key=lambda x: x["structure_inchikey"]
    )
    export_results(
        output=str(output_path / "results_by_structure.tsv"),
        results=results_by_structure,
    )

    logger.info(f"Results exported to {output_path}")


def search_classes(
    classes_file: str | None = None,
    classes_name_id: str | None = None,
    classes_name_smarts: str | None = None,
    closest_only: bool = True,
    include_hierarchy: bool = False,
    input_smiles: str | None = None,
    smiles: str | list[str] | None = None,
    export: bool = True,
    output_dir: Path | str | None = None,
) -> list[dict]:
    """
    Perform substructure search to classify chemical structures.

    This function matches input structures against a set of chemical class
    definitions using SMARTS patterns. Results include the class ID, matching
    SMARTS pattern, and structural similarity metrics.

    :param classes_file: Path to TSV file with chemical class definitions.
        If None, uses the default package classes.
    :param classes_name_id: Column name for class IDs in the classes file.
        Defaults to "class".
    :param classes_name_smarts: Column name for SMARTS in the classes file.
        Defaults to "structure".
    :param closest_only: If True, return only the closest matching class
        for each structure. Default is True.
    :param include_hierarchy: If True, use chemical hierarchy for faster
        searching. Default is False.
    :param input_smiles: Path to file containing SMILES strings to classify.
    :param smiles: Single SMILES string or list of SMILES to classify.
    :param export: If True, export results to files. Default is True.
    :param output_dir: Directory for output files. Uses config default if None.
    :returns: List of dictionaries with classification results.
    :raises ValueError: If no structures are provided and ChEMBL fallback fails.
    """
    # Set column name defaults
    if not classes_name_id:
        classes_name_id = "class"
    if not classes_name_smarts:
        classes_name_smarts = "structure"

    # Collect SMILES from all sources
    smiles_set: set[str] = set()
    if input_smiles:
        smiles_set.update(load_smiles(input=input_smiles))
    if smiles:
        if isinstance(smiles, str):
            smiles_set.add(smiles)
        else:
            smiles_set.update(smiles)

    # Convert SMILES to Mol objects
    structures = _parse_smiles_to_mols(smiles_set)

    # Fallback to ChEMBL if no structures provided
    if not structures:
        logger.info("No structures provided, loading ChEMBL library as fallback.")
        structures = load_latest_chembl()
        if not structures:
            raise ValueError("No structures available for classification.")

    # Load chemical classes
    classes = _load_classes(classes_file, classes_name_id, classes_name_smarts)

    # Load class hierarchy if requested
    class_hierarchy: dict[str, list[str]] = {}
    if include_hierarchy:
        class_hierarchy = load_pkg_chemical_hierarchy()

    # Configure substructure matching parameters
    params = SubstructMatchParameters()
    params.useGenericMatchers = True

    logger.info(
        f"Classifying {len(structures)} structures against {len(classes[0])} chemical classes..."
    )

    # Perform classification
    results = list(
        tqdm_bfs_search_classes_generator(
            classes=classes,
            class_hierarchy=class_hierarchy,
            structures=structures,
            params=params,
        )
    )

    # Filter to closest matches if requested
    if closest_only and results:
        results = _filter_closest_matches(results)

    # Export results
    if export and results:
        output_path = Path(output_dir) if output_dir else None
        _export_classification_results(results, output_path)

    logger.info(f"Classification complete. Found {len(results)} matches.")
    return results


if __name__ == "__main__":
    search_classes()
