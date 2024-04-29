"""Substructure search for chemical classes."""

from __future__ import annotations

import csv
import json
import logging

from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary

from smartclass.chem.classification.bfs_search_classes_generator import (  # dfs_search_classes_generator
    bfs_search_classes_generator,
)
from smartclass.helpers import convert_list_of_dict
from smartclass.io import load_external_classes_file, load_pkg_chemical_hierarchy, load_pkg_classes
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
) -> None:
    """
    Substructure search for chemical classes.

    :param classes_file: File providing the chemical classes to classify.
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
    """
    # Load structures
    s: set[str] = set()
    # TODO io load SMILES
    if input_smiles:
        with open(input_smiles) as f:
            s.update(f.read().splitlines())
    if smiles:
        s.update(smiles)

    # TODO redundant with get_latest_chembl
    mols = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
    for smi in s:
        if smi != "smiles":
            mols.AddSmiles(smi)
    structures = rdSubstructLibrary.SubstructLibrary(mols)

    if not structures:
        logging.info("No structures given, loading ChEMBL library instead.")
        structures = load_latest_chembl()

    # Load classes
    if classes_file:
        # TODO
        c = load_external_classes_file(
            file=classes_file, id_name=classes_name_id, smarts_name=classes_name_smarts
        )
    else:
        logging.info("No classes given, loading default package classes instead.")
        c = load_pkg_classes()
    classes: list = [dict(c.iter_rows())]

    # Load class hierarchy
    class_hierarchy: dict = {}
    if include_hierarchy:
        class_hierarchy = load_pkg_chemical_hierarchy()

    # Use generic matches
    params = Chem.SubstructMatchParameters()
    params.useGenericMatchers = True

    tautomer_insensitive = True
    max_results = 100

    results = list(
        bfs_search_classes_generator(
            classes=classes,
            class_hierarchy=class_hierarchy,
            structures=structures,
            params=params,
            tautomer_insensitive=tautomer_insensitive,
            max_results=max_results,
        )
    )
    # results = list(
    #     dfs_search_classes_generator(
    #         classes=classes,
    #         class_hierarchy=class_hierarchy,
    #         structures=structures,
    #         params=params,
    #         tautomer_insensitive=tautomer_insensitive,
    #         max_results=max_results,
    #     )
    # )

    # Filter the results to keep only the result with the closest class for each unique InChIKey
    if closest_only:
        max_ab_per_inchikey = {}
        for result in results:
            inchikey = result["inchikey"]
            matched_ab = result["matched_ab"]
        
            if inchikey not in max_ab_per_inchikey or matched_ab > max_ab_per_inchikey[inchikey]:
                max_ab_per_inchikey[inchikey] = matched_ab
        results = [result for result in results if result["matched_ab"] == max_ab_per_inchikey[result["inchikey"]]]

    # Export
    key = "class_id"
    value = "inchikey"
    results_kv = convert_list_of_dict(results, key, value)
    # Export results to JSON as key_value
    with open("scratch/results_kv.json", "w") as file:
        json.dump(results_kv, file, indent=4)
    # Export results to CSV
    with open("scratch/results.tsv", "w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=["class_id", "class_structure", "inchikey", "matched_ab"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(results)
    # Revert dict
    results_vk = convert_list_of_dict(results, key, value, invert=True)
    # Export results to JSON as value_key
    with open("scratch/results_vk.json", "w") as file:
        json.dump(results_vk, file, indent=4)
    logging.basicConfig(level=logging.INFO)
    logging.info("Done")


if __name__ == "__main__":
    search_classes()
