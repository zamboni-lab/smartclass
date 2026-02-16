"""Load a Polars DataFrame from the package file."""

from __future__ import annotations

from typing import TYPE_CHECKING

from smartclass.io.load_pkg_file import load_pkg_file


if TYPE_CHECKING:
    from polars import DataFrame

__all__ = ["load_pkg_mappings"]


def load_pkg_mappings() -> DataFrame:
    """
    Load chemont__wd mappings data from the package file into a Polars DataFrame.

    :returns: DataFrame containing chemont__wd mappings.
    :rtype: DataFrame
    """
    return load_pkg_file(file="chemontids__wdids.tsv")


if __name__ == "__main__":
    df = load_pkg_mappings()
