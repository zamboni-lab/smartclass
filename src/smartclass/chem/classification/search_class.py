"""Search for a single class and return the results as a list of dictionaries."""

from __future__ import annotations

# import logging
from rdkit.Chem import FilterCatalog, MolFromSmarts, SubstructMatchParameters

from smartclass.chem.conversion.convert_mol_to_inchikey import convert_mol_to_inchikey

# from smartclass.chem.helpers.enumerate_structures import enumerate_structures
from smartclass.chem.helpers.get_num_matched_atoms_bonds import get_num_matched_atoms_bonds

# from rdkit.Chem import rdTautomerQuery

# from smartclass.chem.similarity.calculate_myopic_mces import calculate_myopic_mces
# from smartclass.chem.similarity.calculate_structural_similarity import (
#     calculate_structural_similarity,
# )


def search_class(
    class_dict: dict[str, list[str]],
    structures: list,
    params: SubstructMatchParameters,
    tautomer_insensitive: bool,
) -> list[dict[str, str]]:
    """
    Search for a single class and return the results as a list of dictionaries.

    :param class_dict: Dictionary of the chemical classes.
    :type class_dict: dict[str, list[str]]

    :param structures: List of structures.
    :type structures: list

    :param params: Parameters for matching.
    :type params: SubstructMatchParameters

    :param tautomer_insensitive: Flag indicating whether tautomer-insensitive search is required.
    :type tautomer_insensitive: bool

    :returns: A list of dictionaries containing matches.
    :rtype: list[dict[str, str]]
    """
    results: list = []

    # COMMENT: This was the first version with enumeration of tautomers and queries.
    # queries = {}
    # for key, value in class_dict.items():
    #     for s in value:
    #         q = enumerate_structures(structure=s, tautomer_insensitive=tautomer_insensitive)
    #         queries[(key, s)] = q
    # for structure in structures:
    #     for (class_id, class_structure), class_structures in queries.items():
    #         for query in class_structures:
    #             try:
    #                 if type(query) is rdTautomerQuery.TautomerQuery:
    #                     matches = query.GetSubstructMatchesWithTautomers(
    #                         structure,
    #                         params,
    #                         )
    #                     for _, match in matches:
    #                         if match:
    #                             results.append(
    #                                 {
    #                                 "class_id": class_id,
    #                                 "class_structure": class_structure,
    #                                 "inchikey": convert_mol_to_inchikey(structure),
    #                                 "matched_ab": get_num_matched_atoms_bonds(
    #                                     mol_1=structure, mol_2=match
    #                                     ),
    #                                 }
    #                                 )
    #                 else:
    #                     if structure.HasSubstructMatch(
    #                         query,
    #                         params,
    #                     ):
    #                         results.append(
    #                             {
    #                                 "class_id": class_id,
    #                                 "class_structure": class_structure,
    #                                 "inchikey": convert_mol_to_inchikey(structure),
    #                                 "matched_ab": get_num_matched_atoms_bonds(
    #                                     mol_1=structure, mol_2=query
    #                                 ),
    #                             }
    #                         )
    #             except Exception as e:
    #                 logging.error(e)
    #                 logging.error(f"Error while searching for class_id {class_id, class_structure}")

    # COMMENT: This was the second version.
    # for class_id, class_structures in class_dict.items():
    #     for class_structure in class_structures:
    #         pattern = MolFromSmarts(class_structure)
    #         for structure in structures:
    #             if structure.HasSubstructMatch(pattern, params):
    #                 result = {
    #                     "class_id": class_id,
    #                     "class_structure": class_structure,
    #                     "inchikey": convert_mol_to_inchikey(structure),
    #                     "matched_ab": get_num_matched_atoms_bonds(mol_1=structure, mol_2=pattern),
    #                 }
    #                 results.append(result)

    # COMMENT: This is the actual version. TODO params could be removed
    for class_id, class_structures in class_dict.items():
        for class_structure in class_structures:
            catalog = FilterCatalog.FilterCatalog()
            pattern = MolFromSmarts(class_structure)
            catalog.AddEntry(
                FilterCatalog.FilterCatalogEntry(class_id, FilterCatalog.SmartsMatcher(pattern))
            )
            for structure in structures:
                if catalog.HasMatch(structure):
                    result = {
                        "class_id": class_id,
                        "class_structure": class_structure,
                        "inchikey": convert_mol_to_inchikey(structure),
                        "matched_ab": get_num_matched_atoms_bonds(mol_1=structure, mol_2=pattern),
                    }
                    results.append(result)

    return results
