"""Load the package data."""

from __future__ import annotations

from typing import TYPE_CHECKING

from smartclass.io.load_pkg_classes import load_pkg_classes
from smartclass.io.load_pkg_mappings import load_pkg_mappings
from smartclass.io.load_pkg_mia import load_pkg_mia
from smartclass.logging import get_logger

if TYPE_CHECKING:
    from polars import DataFrame

__all__ = ["load_pkg_data"]

logger = get_logger(__name__)


def load_pkg_data() -> tuple[DataFrame, DataFrame, DataFrame]:
    """Load the package data.

    Returns
    -------
    tuple[DataFrame, DataFrame, DataFrame]
        DataFrame containing the package data.
    """
    classes = load_pkg_classes()
    mappings = load_pkg_mappings()
    mia = load_pkg_mia()
    return classes, mappings, mia


if __name__ == "__main__":
    classes, mappings, mia = load_pkg_data()
    logger.info(classes)
    logger.info(mappings)
    logger.info(mia)
