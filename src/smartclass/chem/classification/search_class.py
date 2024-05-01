"""Search for a single class and return the results as a list of dictionaries."""

from __future__ import annotations

import logging

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdSubstructLibrary

from smartclass.chem.helpers.enumerate_structures import enumerate_structures

# from smartclass.chem.similarity.calculate_myopic_mces import calculate_myopic_mces
# from smartclass.chem.similarity.calculate_structural_similarity import (
#     calculate_structural_similarity,
# )


def search_class(
    class_id: str,
    class_structure: str,
    structures: rdSubstructLibrary,
    params: Chem.SubstructMatchParameters,
    tautomer_insensitive: bool,
    max_results: int,
) -> list[dict[str, str]]:
    """Search for a single class and return the results as a list of dictionaries."""
    results: list = []

    # Create queries list
    queries = enumerate_structures(
        structure=class_structure, tautomer_insensitive=tautomer_insensitive
    )

    try:
        for query in queries:
            sub = query.GetTemplateMolecule()
            matches = structures.GetMatches(
                query,
                params,
                maxResults=max_results,
            )
            for match in matches:
                mol = structures.GetMol(match)
                mols = [mol, sub]
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
                        "inchikey": Chem.inchi.MolToInchiKey(mol),
                        "matched_ab": num_ab,
                    }
                )
                if len(results) >= max_results:
                    break
            if len(results) >= max_results:
                break
    except Exception as e:
        logging.error(e)
        logging.error(f"Error while searching for class_id {class_id}: {class_structure}")

    return results
