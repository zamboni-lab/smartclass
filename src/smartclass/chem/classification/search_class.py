"""Search for a single class and return the results as a list of dictionaries."""

from __future__ import annotations

import logging

from rdkit import Chem
from rdkit.Chem import rdFMCS

from smartclass.chem.helpers.enumerate_structures import enumerate_structures

# from smartclass.chem.similarity.calculate_myopic_mces import calculate_myopic_mces
# from smartclass.chem.similarity.calculate_structural_similarity import (
#     calculate_structural_similarity,
# )


def search_class(
    class_id: str,
    class_structure: str,
    structures: list,
    params: Chem.SubstructMatchParameters,
    tautomer_insensitive: bool,
) -> list[dict[str, str]]:
    """Search for a single class and return the results as a list of dictionaries."""
    results: list = []

    # Create queries list
    queries = enumerate_structures(
        structure=class_structure, tautomer_insensitive=tautomer_insensitive
    )

    try:
        for structure in structures:
            for query in queries:
                if tautomer_insensitive:
                    matches = query.GetSubstructMatchesWithTautomers(
                        structure,
                        params,
                    )
                    for _, match in matches:
                        if match:
                            mols = [structure, match]
                            mcs = rdFMCS.FindMCS(
                                mols,
                                atomCompare=rdFMCS.AtomCompare.CompareAny,
                                bondCompare=rdFMCS.BondCompare.CompareAny,
                            )
                            num_ab = mcs.numAtoms + mcs.numBonds
                            results.append(
                                {
                                    "class_id": class_id,
                                    "class_structure": class_structure,
                                    "inchikey": Chem.inchi.MolToInchiKey(structure),
                                    "matched_ab": num_ab,
                                }
                            )
                else:
                    if structure.HasSubstructMatch(
                    query,
                    params,
                    ):
                        mols = [structure, query]
                        mcs = rdFMCS.FindMCS(
                            mols,
                            atomCompare=rdFMCS.AtomCompare.CompareAny,
                            bondCompare=rdFMCS.BondCompare.CompareAny,
                        )
                        num_ab = mcs.numAtoms + mcs.numBonds
                        results.append(
                            {
                                "class_id": class_id,
                                "class_structure": class_structure,
                                "inchikey": Chem.inchi.MolToInchiKey(structure),
                                "matched_ab": num_ab,
                            }
                        )
    except Exception as e:
        logging.error(e)
        logging.error(f"Error while searching for class_id {class_id}: {class_structure}")

    # TODO this is not optimal, see later
    return [dict(t) for t in set([tuple(sorted(d.items())) for d in results])]
