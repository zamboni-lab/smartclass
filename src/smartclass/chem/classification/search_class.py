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
        matches = [
            match
            for q in queries
            for match in structures.GetMatches(
                q,
                params,
                maxResults=max_results,
            )
        ]
    except Exception as e:
        logging.error(e)
        logging.error(f"Error while searching for class_id {class_id}: {class_structure}")
        return results

    # Sort matches and yield results
    tmols = [(x, structures.GetMol(x)) for x in matches]
    mols = sorted(tmols, key=lambda x: x[1].GetNumAtoms())

    for _, mol in mols[:max_results]:
        queries = enumerate_structures(
            structure=class_structure, tautomer_insensitive=tautomer_insensitive
        )
        # max_sim = 0
        # min_sim = 666
        num_ab = 0
        for q in queries:
            # sim = calculate_structural_similarity(
            #     mol_1=q.GetTemplateMolecule(), mol_2=Chem.MolFromSmiles(smiles, sanitize=False)
            # )
            # try:
            #     sim = sim[0].similarity
            #     print(f"sim:{sim}")
            # except IndexError:
            #     sim = None
            # if sim > max_sim:
            #     max_sim = sim
            # sim = float(calculate_myopic_mces(
            #     s_1=smiles, s_2=Chem.MolToSmiles(q.GetTemplateMolecule()))[2])
            # print(sim)
            # if sim < min_sim:
            #     min_sim = sim
            mols = [mol, q.GetTemplateMolecule()]
            mcs = rdFMCS.FindMCS(
                mols,
                atomCompare=rdFMCS.AtomCompare.CompareAny,
                bondCompare=rdFMCS.BondCompare.CompareAny,
            )
            sim = mcs.numAtoms + mcs.numBonds
            if sim > num_ab:
                num_ab = sim
        results.append(
            {
                "class_id": class_id,
                "class_structure": class_structure,
                "inchikey": Chem.inchi.MolToInchiKey(mol),
                # "similarity": max_sim,
                # "similarity": min_sim,
                "matched_ab": num_ab,
            }
        )

    return results
