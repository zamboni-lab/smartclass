"""Smartclass classifies structures using SMARTS.resources."""

from __future__ import annotations

from smartclass.resources.chebi.extract_chebi import extract_chebi  # noqa:F401
from smartclass.resources.chebi.get_chebi import get_chebi  # noqa:F401
from smartclass.resources.chembl import get_latest_chembl  # noqa:F401
from smartclass.resources.chembl import latest_chembl_paths  # noqa:F401
from smartclass.resources.chembl import load_latest_chembl  # noqa:F401
from smartclass.resources.get_existing_classes import get_existing_classes  # noqa:F401

# from .lotus import extract_smiles_lotus
from smartclass.resources.wikidata import query_wikidata  # noqa:F401
from smartclass.resources.wikidata import transforms  # noqa:F401
