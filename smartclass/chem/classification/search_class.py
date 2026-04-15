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

Parameters
----------
class_id : str
    Identifier for the chemical class.
class_structures : list[str]
    SMARTS patterns defining the class.

Returns
-------
tuple[FilterCatalog.FilterCatalog, list[tuple[str, Mol]]]
    Filter catalog together with `(original_smarts, compiled_pattern)` entries.
    """
    catalog = FilterCatalog.FilterCatalog()
    pattern_info: list[tuple[str, Mol]] = []

    for class_smarts in class_structures:
        smarts_mol = convert_smarts_to_mol(class_smarts)
        if smarts_mol is None:
            logger.debug(f"Failed to parse SMARTS: {class_smarts}")
            continue

        enumerated_patterns = enumerate_structures(mol=smarts_mol)

        for idx, pattern in enumerate(enumerated_patterns):
            entry_name = f"{class_id}_{idx}"
            catalog.AddEntry(
                FilterCatalog.FilterCatalogEntry(
                    name=entry_name,
                    matcher=FilterCatalog.SmartsMatcher(pattern),
                ),
            )
            pattern_info.append((class_smarts, pattern))

    return catalog, pattern_info


def _find_matching_pattern(
    structure: Mol,
    pattern_info: list[tuple[str, Mol]],
    params: SubstructMatchParameters | None,
) -> tuple[str, Mol] | None:
    """Return the first matching class pattern for a structure.

Parameters
----------
structure : Mol
    Structure to test.
pattern_info : list[tuple[str, Mol]]
    Sequence of `(class_smarts, compiled_pattern)` pairs.
params : SubstructMatchParameters | None
    Optional substructure match parameters.

Returns
-------
tuple[str, Mol] | None
    The first matching `(class_smarts, pattern)` pair, or `None`.
    """
    for class_smarts, pattern in pattern_info:
        if params is None:
            is_match = structure.HasSubstructMatch(pattern)
        else:
            is_match = structure.HasSubstructMatch(pattern, params=params)
        if is_match:
            return class_smarts, pattern
    return None


def search_class(
    class_dict: dict[str, list[str]],
    structures: list[Mol],
    params: SubstructMatchParameters | None = None,
) -> list[dict[str, str | int]]:
    """Search for a single class and return results as dictionaries.

    Performs substructure matching between input structures and class
    SMARTS patterns using RDKit's FilterCatalog for efficient batch
    matching.

Parameters
----------
class_dict : dict[str, list[str]]
    SMARTS patterns.
structures : list[Mol]
    Structures to classify.
params : SubstructMatchParameters | None
    None. Default is None.

Returns
-------
list[dict[str, str | int]]
    Match dictionaries with fields: - class_id: Identifier of matched class - class_structure: SMARTS pattern that matched - structure_inchikey: InChIKey of matched structure - structure_smarts: SMARTS representation of structure - structure_ab: Total atoms + bonds in structure - matched_ab: Number of matched atoms + bonds
    """
    results: list[dict[str, str | int]] = []
    structure_cache: dict[int, dict[str, str | int]] = {}

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

            matched = _find_matching_pattern(structure, pattern_info, params)
            if matched is None:
                # Defensive fallback; catalog matched but no specific pattern found.
                continue

            class_smarts, pattern = matched
            structure_id = id(structure)
            if structure_id not in structure_cache:
                structure_cache[structure_id] = {
                    "structure_inchikey": convert_mol_to_inchikey(structure),
                    "structure_smarts": convert_mol_to_smarts(structure),
                    "structure_ab": get_num_atoms_bonds(structure),
                }

            cached = structure_cache[structure_id]
            results.append({
                "class_id": class_id,
                "class_structure": class_smarts,
                "structure_inchikey": cached["structure_inchikey"],
                "structure_smarts": cached["structure_smarts"],
                "structure_ab": cached["structure_ab"],
                "matched_ab": pattern.GetNumAtoms() + pattern.GetNumBonds(),
            })

    return results
