"""Search for a single class and return the results as a list of dictionaries."""

from __future__ import annotations

import logging

from rdkit.Chem import SubstructMatchParameters, rdTautomerQuery

from smartclass.chem.conversion.convert_mol_to_inchikey import convert_mol_to_inchikey
from smartclass.chem.helpers.enumerate_structures import enumerate_structures
from smartclass.chem.helpers.get_num_matched_atoms_bonds import get_num_matched_atoms_bonds

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
    :type class_dict: Union[None,str]

    :param structures: List of structures.
    :type structures: list

    :param params: Parameters for matching.
    :type params: SubstructMatchParameters

    :param tautomer_insensitive: Flag indicating whether tautomer-insensitive search is required.
    :type tautomer_insensitive: bool

    :returns: A list of dictionaries containing matches.
    :rtype: list[dict]
    """
    results: list = []

    queries: list[dict[str, dict[str, list[str]]]] = []
    for k1, v1 in class_dict.items():
        for s1 in v1:
            q = enumerate_structures(structure=s1, tautomer_insensitive=tautomer_insensitive)
            queries.append({k1: q})

    # TODO WIP #5
    # matches = []
    # for k, v in smarts_library.items():
    #     smarts_pattern = Chem.MolFromSmarts(v)
    #     match = m.GetSubstructMatch(smarts_pattern)
    #     if match:
    #         matches.append((k, v))

    for structure in structures:
        for query in queries:
            for k2, v2 in query.items():
                for k3, v3 in v2.items():
                    for s in v3:
                        try:
                            if type(s) is rdTautomerQuery.TautomerQuery:
                                matches = s.GetSubstructMatchesWithTautomers(
                                    structure,
                                    params,
                                )
                                for _, match in matches:
                                    if match:
                                        results.append(
                                            {
                                                "class_id": k2,
                                                "class_structure": k3,
                                                "inchikey": convert_mol_to_inchikey(structure),
                                                "matched_ab": get_num_matched_atoms_bonds(
                                                    mol_1=structure, mol_2=match
                                                ),
                                            }
                                        )
                            else:
                                if structure.HasSubstructMatch(
                                    s,
                                    params,
                                ):
                                    results.append(
                                        {
                                            "class_id": k2,
                                            "class_structure": k3,
                                            "inchikey": convert_mol_to_inchikey(structure),
                                            "matched_ab": get_num_matched_atoms_bonds(
                                                mol_1=structure, mol_2=s
                                            ),
                                        }
                                    )
                        except Exception as e:
                            logging.error(e)
                            logging.error(f"Error while searching for class_id {k2, k3}")

    # TODO this is not optimal, see later
    return [dict(t) for t in set([tuple(sorted(d.items())) for d in results])]
