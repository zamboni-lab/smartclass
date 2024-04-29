"""Load the package data."""

from __future__ import annotations

from typing import TYPE_CHECKING

from smartclass.io.load_pkg_classes import load_pkg_classes
from smartclass.io.load_pkg_mappings import load_pkg_mappings
from smartclass.io.load_pkg_mia import load_pkg_mia

if TYPE_CHECKING:
    from polars import DataFrame

__all__ = ["load_pkg_data"]


def load_pkg_data() -> tuple[DataFrame, DataFrame, DataFrame]:
    """
    Load the package data.

    :returns: A tuple of Polars DataFrame containing the package data.
    :rtype: tuple[DataFrame]
    """
    classes = load_pkg_classes()
    mappings = load_pkg_mappings()
    mia = load_pkg_mia()
    return classes, mappings, mia


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)

    classes, mappings, mia = load_pkg_data()
    logging.info(classes)
    logging.info(mappings)
    logging.info(mia)
