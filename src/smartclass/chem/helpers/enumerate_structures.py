"""Enumerate structures."""

from __future__ import annotations

import logging

from rdkit import Chem
from rdkit.Chem import rdMolEnumerator, rdTautomerQuery

__all__ = [
    "enumerate_structures",
]


def enumerate_structures(
    structure: str,
    tautomer_insensitive: bool,
) -> rdTautomerQuery:
    """Enumerate structures."""
    query = Chem.MolFromSmarts(structure)
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
        try:
            query = rdTautomerQuery.TautomerQuery(q) if tautomer_insensitive else q
            queries.append(query)
        except Exception as e:
            logging.error(e)

    # TODO See https://github.com/rdkit/rdkit/commit/908e47cc03607b86eb58cd5512555b741356f15e

    return queries
