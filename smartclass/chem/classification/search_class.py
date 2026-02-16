"""Search for a single class and return the results as a list of dictionaries.

This module performs substructure matching between chemical structures and
SMARTS-based class definitions using RDKit's FilterCatalog for efficiency.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from rdkit.Chem import FilterCatalog, SubstructMatchParameters

from smartclass.chem.conversion.convert_mol_to_inchikey import convert_mol_to_inchikey
from smartclass.chem.conversion.convert_mol_to_smarts import convert_mol_to_smarts
from smartclass.chem.conversion.convert_smarts_to_mol import convert_smarts_to_mol
from smartclass.chem.helpers.enumerate_structures import enumerate_structures
from smartclass.chem.helpers.get_num_atoms_bonds import get_num_atoms_bonds
from smartclass.chem.helpers.get_num_matched_atoms_bonds import (
    get_num_matched_atoms_bonds,
)
from smartclass.logging import get_logger


if TYPE_CHECKING:
    from rdkit.Chem import Mol

__all__ = [
    "search_class",
    "build_filter_catalog",
]

logger = get_logger(__name__)


def build_filter_catalog(
    class_id: str,
    class_structures: list[str],
) -> tuple[FilterCatalog.FilterCatalog, list[tuple[str, Mol]]]:
    """Build a FilterCatalog from SMARTS patterns.

    Pre-builds the catalog and enumerated patterns for reuse across
    multiple structure queries.

    Args:
        class_id: Identifier for the chemical class.
        class_structures: List of SMARTS patterns defining the class.

    Returns:
        Tuple of (FilterCatalog, list of (smarts, pattern_mol) tuples).
    """
    catalog = FilterCatalog.FilterCatalog()
    pattern_info: list[tuple[str, Mol]] = []

    for class_smarts in class_structures:
        smarts_mol = convert_smarts_to_mol(class_smarts)
        if smarts_mol is None:
            logger.debug(f"Failed to parse SMARTS: {class_smarts[:50]}...")
            continue

        enumerated_patterns = enumerate_structures(mol=smarts_mol)

        for idx, pattern in enumerate(enumerated_patterns):
            entry_name = f"{class_id}_{idx}"
            catalog.AddEntry(
                FilterCatalog.FilterCatalogEntry(
                    name=entry_name,
                    matcher=FilterCatalog.SmartsMatcher(pattern),
                )
            )
            pattern_info.append((class_smarts, pattern))

    return catalog, pattern_info


def search_class(
    class_dict: dict[str, list[str]],
    structures: list[Mol],
    params: SubstructMatchParameters | None = None,
) -> list[dict[str, str | int]]:
    """Search for a single class and return results as dictionaries.

    Performs substructure matching between input structures and class
    SMARTS patterns using RDKit's FilterCatalog for efficient batch
    matching.

    Args:
        class_dict: Dictionary mapping class_id to list of SMARTS patterns.
        structures: List of RDKit Mol objects to classify.
        params: Optional SubstructMatchParameters (currently unused but
            kept for API compatibility).

    Returns:
        List of dictionaries with match information:
        - class_id: Identifier of matched class
        - class_structure: SMARTS pattern that matched
        - structure_inchikey: InChIKey of matched structure
        - structure_smarts: SMARTS representation of structure
        - structure_ab: Total atoms + bonds in structure
        - matched_ab: Number of matched atoms + bonds
    """
    results: list[dict[str, str | int]] = []

    for class_id, class_structures in class_dict.items():
        # Build catalog once per class for efficiency
        catalog, pattern_info = build_filter_catalog(class_id, class_structures)

        if not pattern_info:
            logger.debug(f"No valid patterns for class {class_id}")
            continue

        # Process each structure against the pre-built catalog
        for structure in structures:
            if not catalog.HasMatch(structure):
                continue

            # Get the last pattern for matched_ab calculation
            # TODO: Consider tracking which specific pattern matched
            class_smarts, pattern = pattern_info[-1]

            results.append({
                "class_id": class_id,
                "class_structure": class_smarts,
                "structure_inchikey": convert_mol_to_inchikey(structure),
                "structure_smarts": convert_mol_to_smarts(structure),
                "structure_ab": get_num_atoms_bonds(structure),
                "matched_ab": get_num_matched_atoms_bonds(
                    mol_1=structure,
                    mol_2=pattern,
                ),
            })

    return results
