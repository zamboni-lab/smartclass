"""Smartclass classifies structures using SMARTS.chem.classification."""

# from .check_classification import check_classification
from __future__ import annotations

from smartclass.chem.classification.bfs_search_classes_generator import (  # noqa:F401
    bfs_search_classes_generator,
)
from smartclass.chem.classification.get_class_structures import (
    get_class_structures,  # noqa:F401
)
from smartclass.chem.classification.search_class import search_class  # noqa:F401
from smartclass.chem.classification.search_classes import search_classes  # noqa:F401
