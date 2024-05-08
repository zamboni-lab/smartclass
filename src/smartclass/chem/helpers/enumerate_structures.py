"""Enumerate structures."""

from __future__ import annotations

import logging

from rdkit import Chem
from rdkit.Chem import Mol, rdMolEnumerator, rdTautomerQuery

__all__ = [
    "enumerate_structures",
]


def enumerate_structures(
    structure: str,
    tautomer_insensitive: bool,
) -> list[rdTautomerQuery.TautomerQuery | Mol]:
    """Enumerate structures."""
    query = Chem.MolFromSmarts(structure)
    # TODO looks important
    # Chem.Kekulize(query)
    Chem.SetGenericQueriesFromProperties(query)

    # Enumerate the query molecules
    try:
        bndl = rdMolEnumerator.Enumerate(query)
    except Exception as e:
        bndl = []
        logging.error(e)

    # Check if enumeration was successful
    if len(bndl) == 0:
        # If enumeration failed, use the original query as a fallback
        bndl = [query]

    # Create queries list
    queries = []
    for q in bndl:
        # TODO see if needed
        # Seems more like it is not for now
        q = Chem.AdjustQueryProperties(q)
        if tautomer_insensitive:
            try:
                q = rdTautomerQuery.TautomerQuery(q)
            except Exception as e:
                logging.error(e)
        queries.append(q)

    # TODO See https://github.com/rdkit/rdkit/commit/908e47cc03607b86eb58cd5512555b741356f15e

    return queries
