"""Load data files bundled with the smartclass package."""

from __future__ import annotations

from typing import TYPE_CHECKING

import importlib_resources

from smartclass.exceptions import DataLoadingError
from smartclass.io.load_csv_from_path import load_csv_from_path
from smartclass.io.load_json_from_path import load_json_from_path
from smartclass.io.load_tsv_from_path import load_tsv_from_path
from smartclass.logging import get_logger


if TYPE_CHECKING:
    from polars import DataFrame

__all__ = ["load_pkg_file"]

logger = get_logger(__name__)


def load_pkg_file(file: str, directory: str = "smartclass.data") -> DataFrame:
    """
    Load data from a package-bundled file into a Polars DataFrame.

    Supports TSV, CSV, and JSON file formats.

    :param file: Name of the file to load (e.g., "classes_smarts.tsv").
    :param directory: Package directory containing the file.
        Defaults to "smartclass.data".
    :returns: Polars DataFrame containing the loaded data.
    :raises DataLoadingError: If the file cannot be loaded or format is unsupported.
    """
    try:
        ref = importlib_resources.files(directory) / file
        with importlib_resources.as_file(ref) as path:
            if file.endswith(".tsv"):
                return load_tsv_from_path(path)
            elif file.endswith(".csv"):
                return load_csv_from_path(path)
            elif file.endswith(".json"):
                from polars import DataFrame

                return DataFrame(load_json_from_path(path))
            else:
                raise DataLoadingError(
                    file,
                    reason=f"Unsupported file format: {file.split('.')[-1]}",
                )
    except DataLoadingError:
        raise
    except Exception as e:
        raise DataLoadingError(file, reason=str(e)) from e


if __name__ == "__main__":
    df1 = load_pkg_file(file="chemontids__wdids.tsv")
    df2 = load_pkg_file(file="classes_smarts.tsv")
    df3 = load_pkg_file(file="mia_smarts.tsv")
