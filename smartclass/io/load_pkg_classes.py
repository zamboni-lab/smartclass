"""Load a Polars DataFrame from the package file."""

from __future__ import annotations

from typing import TYPE_CHECKING

from smartclass.io.load_pkg_file import load_pkg_file


if TYPE_CHECKING:
    from polars import DataFrame

__all__ = ["load_pkg_classes"]


def load_pkg_classes() -> DataFrame:
    """
    Load chemical classes data from the package file into a Polars DataFrame.

    :returns: DataFrame containing chemical classes.
    :rtype: DataFrame
    """
    return load_pkg_file(file="classes_smarts.tsv")


if __name__ == "__main__":
    df = load_pkg_classes()
