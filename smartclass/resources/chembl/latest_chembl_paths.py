"""Retrieve the latest ChEMBL version and generate file paths for data and library."""

from __future__ import annotations

import logging

from chembl_downloader import latest
from pystow import join


__all__ = [
    "latest_chembl_paths",
]


def latest_chembl_paths() -> tuple:
    """
    Retrieve the latest ChEMBL version and generate file paths for data and library.

    :raises RuntimeError: If the ChEMBL version or file paths cannot be retrieved.

    :returns: A tuple containing the ChEMBL version, data path, and library path.
    :rtype: tuple
    """
    try:
        # Get the latest ChEMBL version
        version = latest()

        # Generate file paths for data and library
        data_path = join("chembl", version, name="sssdata.pkl")
        lib_path = join("chembl", version, name="ssslib.pkl")

        return version, data_path, lib_path
    except Exception as e:
        raise RuntimeError(f"Failed to retrieve ChEMBL paths: {e!s}") from e


if __name__ == "__main__":
    try:
        version, data_path, lib_path = latest_chembl_paths()
        logging.debug(f"ChEMBL version: {version}")
        logging.debug(f"Data path: {data_path}")
        logging.debug(f"Library path: {lib_path}")
    except RuntimeError as e:
        logging.exception(f"Error: {e}")
