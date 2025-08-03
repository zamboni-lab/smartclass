"""Substructure search for chemical classes."""

from __future__ import annotations

import csv
import json
import logging

from rdkit.Chem import SubstructMatchParameters

from smartclass.chem.classification.bfs_search_classes_generator import (
    tqdm_bfs_search_classes_generator,
)
from smartclass.chem.conversion.convert_smiles_to_mol import convert_smiles_to_mol
from smartclass.helpers import convert_list_of_dict
from smartclass.io import (
    load_external_classes_file,
    load_pkg_chemical_hierarchy,
    load_pkg_classes,
    load_smiles,
)
from smartclass.resources.chembl import load_latest_chembl

__all__ = [
    "search_classes",
]


def search_classes(
    classes_file: None | str = None,
    classes_name_id: None | str = None,
    classes_name_smarts: None | str = None,
    closest_only: bool = True,
    include_hierarchy: bool = False,
    input_smiles: None | str = None,
    smiles: None | (str | list[str]) = None,
) -> list[dict]:
    """
    Substructure search for chemical classes.

    :param classes_file: File providing the chemical classes.
    :type classes_file: Union[None,str]

    :param classes_name_id: Name of the ID column in the classes file.
    :type classes_name_id: Union[None,str]

    :param classes_name_smarts: Name of the SMARTS column in the classes file.
    :type classes_name_smarts: Union[None,str]

    :param closest_only: Flag to return only the closest class.
    :type closest_only: bool

    :param include_hierarchy: Flag to include hierarchy search (default is False).
    :type include_hierarchy: bool

    :param input_smiles: File providing the (list of) structure(s) to classify.
    :type input_smiles: Union[None,str]

    :param smiles: (List of) structure(s) to classify.
    :type smiles: Union[None,str,list[str]]

    :returns: A list of matched classes.
    :rtype: list[dict]
    """
    # Load structures
    s: set[str] = set()
    if input_smiles:
        s.update(load_smiles(input=input_smiles))
    if smiles:
        s.update(smiles)

    # TODO redundant with get_latest_chembl
    structures: list = list()
    for smi in s:
        if smi != "smiles":
            mol = convert_smiles_to_mol(smi)
            # TODO looks important
            # Kekulize(mol)
            if mol is not None:
                structures.append(mol)

    # TODO change this
    if not structures:
        logging.basicConfig(level=logging.INFO)
        logging.info("No structures given, loading ChEMBL library instead.")
        structures = load_latest_chembl()

    # Load classes
    if classes_file:
        c = load_external_classes_file(
            file=classes_file, id_name=classes_name_id, smarts_name=classes_name_smarts
        )
    else:
        logging.basicConfig(level=logging.INFO)
        logging.info("No classes given, loading default package classes instead.")
        c = load_pkg_classes()
    classes = []
    classes_dict: dict[str, list[str]] = {}
    for row in c.iter_rows():
        key = row[0]
        value = row[1]
        if key in classes_dict:
            classes_dict[key].append(value)
        else:
            classes_dict[key] = [value]
    classes.append(classes_dict)

    # Load class hierarchy
    class_hierarchy: dict[str, list[str]] = {}
    if include_hierarchy:
        class_hierarchy = load_pkg_chemical_hierarchy()

    # Use generic matches
    params = SubstructMatchParameters()
    params.useGenericMatchers = True

    logging.basicConfig(level=logging.INFO)
    logging.info(f"Classifying {len(structures)} structures...")
    logging.basicConfig(level=logging.INFO)
    logging.info(f"...against {len(classes[0])} chemical classes...")

    results = list(
        tqdm_bfs_search_classes_generator(
            classes=classes,
            class_hierarchy=class_hierarchy,
            structures=structures,
            params=params,
        ),
    )
    # results = list(
    #     dfs_search_classes_generator(
    #         classes=classes,
    #         class_hierarchy=class_hierarchy,
    #         structures=structures,
    #         params=params,
    #     )
    # )

    # Filter the results to keep only the result with the closest class for each unique InChIKey
    if closest_only:
        max_ab_per_inchikey: dict = {}
        for result in results:
            inchikey = result["structure_inchikey"]
            matched_ab = result["matched_ab"]

            if (
                inchikey not in max_ab_per_inchikey
                or matched_ab > max_ab_per_inchikey[inchikey]
            ):
                max_ab_per_inchikey[inchikey] = matched_ab
        results = [
            result
            for result in results
            if result["matched_ab"] == max_ab_per_inchikey[result["structure_inchikey"]]
        ]

    # Export
    key = "class_id"
    value = "structure_inchikey"
    fields = [
        "class_id",
        "class_structure",
        "structure_inchikey",
        "structure_smarts",
        "structure_ab",
        "matched_ab",
    ]
    results_kv = convert_list_of_dict(results, key, value)
    # Export results to JSON as key_value
    with open("scratch/results_kv.json", "w") as file:
        json.dump(results_kv, file, indent=4)
    # Export results to CSV
    results_sorted = sorted(results, key=lambda x: x["matched_ab"], reverse=True)
    with open("scratch/results_kv.tsv", "w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(results)
    results_sorted = sorted(results_sorted, key=lambda x: x["structure_inchikey"])
    with open("scratch/results_vk.tsv", "w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(results_sorted)

    # Revert dict
    results_vk = convert_list_of_dict(results, key, value, invert=True)
    # Export results to JSON as value_key
    with open("scratch/results_vk.json", "w") as file:
        json.dump(results_vk, file, indent=4)
    logging.basicConfig(level=logging.INFO)
    logging.info("Done")

    return results


if __name__ == "__main__":
    search_classes()
