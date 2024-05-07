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
    class_id: str,
    class_structure: str,
    structures: list,
    params: SubstructMatchParameters,
    tautomer_insensitive: bool,
) -> list[dict[str, str]]:
    """Search for a single class and return the results as a list of dictionaries."""
    results: list = []

    # Create queries list
    queries = []
    q = enumerate_structures(structure=class_structure, tautomer_insensitive=tautomer_insensitive)
    queries.extend(q)
    for structure in structures:
        try:
            for query in queries:
                if type(query) is rdTautomerQuery.TautomerQuery:
                    matches = query.GetSubstructMatchesWithTautomers(
                        structure,
                        params,
                    )
                    for _, match in matches:
                        if match:
                            results.append(
                                {
                                    "class_id": class_id,
                                    "class_structure": class_structure,
                                    "inchikey": convert_mol_to_inchikey(structure),
                                    "matched_ab": get_num_matched_atoms_bonds(
                                        mol_1=structure, mol_2=match
                                    ),
                                }
                            )
                else:
                    if structure.HasSubstructMatch(
                        query,
                        params,
                    ):
                        results.append(
                            {
                                "class_id": class_id,
                                "class_structure": class_structure,
                                "inchikey": convert_mol_to_inchikey(structure),
                                "matched_ab": get_num_matched_atoms_bonds(
                                    mol_1=structure, mol_2=query
                                ),
                            }
                        )
        except Exception as e:
            logging.error(e)
            logging.error(f"Error while searching for class_id {class_id}: {class_structure}")

    # TODO this is not optimal, see later
    return [dict(t) for t in set([tuple(sorted(d.items())) for d in results])]
