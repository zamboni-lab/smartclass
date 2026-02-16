"""Perform substructure search for chemical classes and yield results.

This module implements Breadth-First Search (BFS) traversal of chemical
class hierarchies for efficient substructure matching.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Generator, Iterator
from typing import TYPE_CHECKING

from tqdm import tqdm

from smartclass.chem.classification.get_class_structures import get_class_structures
from smartclass.chem.classification.search_class import search_class
from smartclass.logging import get_logger


if TYPE_CHECKING:
    from rdkit.Chem import Mol, SubstructMatchParameters

__all__ = [
    "bfs_search_classes_generator",
    # "dfs_search_classes_generator"
    "tqdm_bfs_search_classes_generator",
]

logger = get_logger(__name__)

# Type aliases for clarity
ClassDict = dict[str, list[str]]
ClassHierarchy = dict[str, list[str]]
MatchResult = dict[str, str | int]


def _search_single_class(
    class_id: str,
    class_structure: list[str],
    structures: list[Mol],
    params: SubstructMatchParameters,
) -> Iterator[MatchResult]:
    """Search a single class against all structures.

    Args:
        class_id: Identifier for the chemical class.
        class_structure: List of SMARTS patterns for this class.
        structures: List of RDKit Mol objects to search.
        params: Substructure matching parameters.

    Yields:
        Match result dictionaries.
    """
    results = search_class(
        class_dict={class_id: class_structure},
        structures=structures,
        params=params,
    )
    yield from results


def _get_hierarchical_classes(
    class_hierarchy: ClassHierarchy,
    class_structures: ClassDict,
) -> Generator[tuple[str, list[str]], None, set[str]]:
    """Traverse hierarchy using BFS and yield classes with their structures.

    Args:
        class_hierarchy: Parent-to-children mapping.
        class_structures: Class ID to SMARTS patterns mapping.

    Yields:
        Tuples of (class_id, class_structure_list).

    Returns:
        Set of searched class IDs.
    """
    searched: set[str] = set()
    queue: deque[str] = deque(class_hierarchy.keys())
    visited: set[str] = set()

    while queue:
        current_class_id = queue.popleft()

        if current_class_id in visited:
            continue
        visited.add(current_class_id)

        class_structure = class_structures.get(current_class_id)
        if class_structure is None:
            logger.debug(f"No SMARTS pattern found for class: {current_class_id}")
            continue

        yield current_class_id, class_structure
        searched.add(current_class_id)

        # Add children to queue
        for child_class in class_hierarchy.get(current_class_id, []):
            if child_class not in visited and child_class not in searched:
                queue.append(child_class)

    return searched


def bfs_search_classes_generator(
    classes: list[ClassDict],
    structures: list[Mol],
    params: SubstructMatchParameters,
    class_hierarchy: ClassHierarchy | None = None,
) -> Generator[MatchResult]:
    """Perform substructure search using Breadth-First Search traversal.

    When a class hierarchy is provided, classes are searched in BFS order
    starting from root classes. This can improve efficiency by allowing
    early termination when parent classes don't match.

    Args:
        classes: List of dictionaries mapping class IDs to SMARTS patterns.
        structures: List of RDKit Mol objects to classify.
        params: Parameters for substructure matching.
        class_hierarchy: Optional parent-to-children mapping for hierarchical
            search. Keys are parent class IDs, values are lists of child IDs.

    Yields:
        Dictionary containing match information for each structure-class pair.
    """
    class_structures = get_class_structures(classes)
    searched_classes: set[str] = set()

    # Phase 1: Search hierarchical classes (BFS order)
    if class_hierarchy:
        queue: deque[str] = deque(class_hierarchy.keys())
        visited: set[str] = set()

        while queue:
            current_class_id = queue.popleft()

            if current_class_id in visited:
                continue
            visited.add(current_class_id)

            class_structure = class_structures.get(current_class_id)
            if class_structure is None:
                logger.debug(f"No SMARTS pattern for class: {current_class_id}")
                continue

            yield from _search_single_class(
                current_class_id,
                class_structure,
                structures,
                params,
            )
            searched_classes.add(current_class_id)

            # Add children to queue
            for child_class in class_hierarchy.get(current_class_id, []):
                if child_class not in visited and child_class not in searched_classes:
                    queue.append(child_class)

    # Phase 2: Search remaining non-hierarchical classes
    for class_dict in classes:
        for class_id, class_structure_list in class_dict.items():
            # Skip if already searched or is a child in hierarchy
            if class_id in searched_classes:
                continue
            if class_hierarchy and class_id in class_hierarchy:
                continue

            yield from _search_single_class(
                class_id,
                class_structure_list,
                structures,
                params,
            )


def tqdm_bfs_search_classes_generator(
    classes: list[ClassDict],
    structures: list[Mol],
    params: SubstructMatchParameters,
    class_hierarchy: ClassHierarchy | None = None,
) -> Generator[MatchResult]:
    """BFS search with tqdm progress bar.

    Wraps bfs_search_classes_generator with a progress indicator.
    See bfs_search_classes_generator for parameter documentation.

    Yields:
        Dictionary containing match information for each match found.
    """
    with tqdm(desc="Searching classes", unit=" matches") as pbar:
        for result in bfs_search_classes_generator(
            classes,
            structures,
            params,
            class_hierarchy,
        ):
            yield result
            pbar.update(1)


if __name__ == "__main__":
    from rdkit.Chem import SubstructMatchParameters

    from smartclass.io import load_pkg_chemical_hierarchy, load_pkg_classes
    from smartclass.resources.chembl import load_latest_chembl

    # Load classes
    c = load_pkg_classes()
    classes: list = [dict(c.iter_rows())]

    # Load class hierarchy
    class_hierarchy = load_pkg_chemical_hierarchy()

    # Load the ChEMBL library
    structures = load_latest_chembl()

    # Use generic matches
    params = SubstructMatchParameters()
    params.maxMatches = 100
    params.useGenericMatchers = True

    # Perform BFS search
    results_bfs = list(
        bfs_search_classes_generator(
            classes=classes,
            class_hierarchy=class_hierarchy,
            structures=structures,
            params=params,
        ),
    )
    print(f"Found {len(results_bfs)} matches")
