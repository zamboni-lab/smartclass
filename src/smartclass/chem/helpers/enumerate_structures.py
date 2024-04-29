"""Enumerate structures."""

from __future__ import annotations

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
    bndl = rdMolEnumerator.Enumerate(query)
    if bndl.Size() == 0:
        bndl = [query]

    # Create queries list
    queries = [rdTautomerQuery.TautomerQuery(q) if tautomer_insensitive else q for q in bndl]

    # TODO See https://github.com/rdkit/rdkit/commit/908e47cc03607b86eb58cd5512555b741356f15e

    return queries
