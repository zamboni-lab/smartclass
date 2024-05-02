"""Perform substructure search for chemical classes and yield results."""

from __future__ import annotations

import logging
from collections.abc import Generator
from typing import TYPE_CHECKING

from rdkit.Chem import SubstructMatchParameters
from tqdm import tqdm

from smartclass.chem.classification.get_class_structures import get_class_structures
from smartclass.chem.classification.search_class import search_class

if TYPE_CHECKING:
    from rdkit.Chem import rdSubstructLibrary

__all__ = [
    "bfs_search_classes_generator",
    # "dfs_search_classes_generator"
    "tqdm_bfs_search_classes_generator",
]


def bfs_search_classes_generator(
    classes: list[dict[str, list[str]]],
    structures: rdSubstructLibrary,
    params: SubstructMatchParameters,
    tautomer_insensitive: bool,
    class_hierarchy: dict[str, list[str]] | None = None,
) -> Generator:
    """
    Perform substructure search for chemical classes using Breadth-First Search (BFS) and yield results.

    :param classes: List of chemical classes to search.
    :param structures: RDKit substructure library for searching.
    :param params: Parameters for substructure search.
    :param tautomer_insensitive: Flag indicating whether tautomer-insensitive search is required.
    :param class_hierarchy: Dictionary representing hierarchy between classes (optional).
    :yields: A dictionary containing class_id, result_id, and smiles for each match.
    """
    class_structures = get_class_structures(classes)
    searched_classes = set()

    if class_hierarchy:
        queue = list(class_hierarchy.keys())
        visited = set()

        while queue:
            current_class_id = queue.pop(0)
            visited.add(current_class_id)

            class_structure = class_structures.get(current_class_id)
            if class_structure is None:
                logging.error(f"No class structure found for class_id {current_class_id}")
                continue

            for class_s in class_structure:
                results = search_class(
                    current_class_id,
                    class_s,
                    structures,
                    params,
                    tautomer_insensitive,
                )

                for result in results:
                    yield result

            searched_classes.add(current_class_id)

            if current_class_id in class_hierarchy:
                for child_class in class_hierarchy[current_class_id]:
                    if child_class not in visited and child_class not in searched_classes:
                        queue.append(child_class)

    for class_id, class_structure in class_structures.items():
        if class_id not in searched_classes and (
            class_hierarchy is None or class_id not in class_hierarchy
        ):
            for class_s in class_structure:
                results = search_class(
                    class_id,
                    class_s,
                    structures,
                    params,
                    tautomer_insensitive,
                )

                for result in results:
                    yield result

            # results = search_class(
            #     class_id,
            #     class_structure,
            #     structures,
            #     params,
            #     tautomer_insensitive,
            # )

            # for result in results:
            #     yield result


def tqdm_bfs_search_classes_generator(
    classes: list[dict[str, list[str]]],
    structures: rdSubstructLibrary,
    params: SubstructMatchParameters,
    tautomer_insensitive: bool,
    class_hierarchy: dict[str, list[str]] | None = None,
) -> Generator:
    """
    Perform substructure search for chemical classes using Breadth-First Search (BFS) and yield results with tqdm.

    :param classes: List of chemical classes to search.
    :param structures: RDKit substructure library for searching.
    :param params: Parameters for substructure search.
    :param tautomer_insensitive: Flag indicating whether tautomer-insensitive search is required.
    :param class_hierarchy: Dictionary representing hierarchy between classes (optional).
    :yields: A dictionary containing class_id, result_id, and smiles for each match.
    """
    with tqdm(
        desc="Searching classes",
    ) as pbar:
        for result in bfs_search_classes_generator(
            classes, structures, params, tautomer_insensitive, class_hierarchy
        ):
            yield result
            pbar.update(1)


# def dfs_search_classes_generator(
#         classes: List[Dict[str, str]],
#         structures: rdSubstructLibrary,
#         params: SubstructMatchParameters,
#         tautomer_insensitive: bool,
#         class_hierarchy: Optional[Dict[str, List[str]]] = None,
# ) -> Generator:
#     """
#     Perform substructure search for chemical classes using Depth-First Search (DFS) and yield results.

#     :param classes: List of chemical classes to search.
#     :param structures: RDKit substructure library for searching.
#     :param params: Parameters for substructure search.
#     :param tautomer_insensitive: Flag indicating whether tautomer-insensitive search is required.
#     :param class_hierarchy: Dictionary representing hierarchy between classes (optional).
#     :yields: A dictionary containing class_id, result_id, and smiles for each match.
#     """
#     class_structures = get_class_structures(classes)
#     searched_classes = set()

#     def dfs(class_id):
#         if class_id in searched_classes:
#             return

#         class_structure = class_structures.get(class_id)
#         if class_structure is None:
#             logging.error(f"No class structure found for class_id {class_id}")
#             return

#         results = search_class(
#             class_id,
#             class_structure,
#             structures,
#             params,
#             tautomer_insensitive,
#         )

#         for result in results:
#             yield result

#         searched_classes.add(class_id)

#         if class_hierarchy and class_id in class_hierarchy:
#             for child_class in class_hierarchy[class_id]:
#                 dfs(child_class)

#     for class_id in class_structures.keys():
#         if class_id not in searched_classes and (not class_hierarchy or class_id not in class_hierarchy):
#             dfs(class_id)

if __name__ == "__main__":
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

    tautomer_insensitive = True

    # Perform BFS search
    results_bfs = list(
        bfs_search_classes_generator(
            classes=classes,
            class_hierarchy=class_hierarchy,
            structures=structures,
            params=params,
            tautomer_insensitive=tautomer_insensitive,
        )
    )

    # Perform DFS search
    # results_dfs = list(
    #     dfs_search_classes_generator(
    #         classes=classes,
    #         class_hierarchy=class_hierarchy,
    #         structures=structures,
    #         params=params,
    #         tautomer_insensitive=tautomer_insensitive,
    #     )
    # )
