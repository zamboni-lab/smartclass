"""Retrieve the latest ChEMBL version and generate file paths for data and library."""

from __future__ import annotations

from typing import TYPE_CHECKING

from chembl_downloader import latest
from pystow import join

from smartclass.logging import get_logger

if TYPE_CHECKING:
    from pathlib import Path

    from chembl_downloader.api import VersionInfo

__all__ = [
    "latest_chembl_paths",
]

logger = get_logger(__name__)


def latest_chembl_paths() -> tuple[str | VersionInfo, Path, Path]:
    """Retrieve the latest ChEMBL version and generate file paths for data and library.

    Returns
    -------
    tuple[str | VersionInfo, Path, Path]
        ChEMBL version, data path, and library path.

    Raises
    ------
    RuntimeError
        If the ChEMBL version or file paths cannot be retrieved.
    """
    try:
        # Get the latest ChEMBL version
        version = latest()
        version_str = str(version)

        # Generate file paths for data and library
        data_path = join("chembl", version_str, name="sssdata.pkl")
        lib_path = join("chembl", version_str, name="ssslib.pkl")

        return version, data_path, lib_path
    except Exception as e:
        raise RuntimeError(f"Failed to retrieve ChEMBL paths: {e!s}") from e


if __name__ == "__main__":
    try:
        version, data_path, lib_path = latest_chembl_paths()
        logger.debug(f"ChEMBL version: {version}")
        logger.debug(f"Data path: {data_path}")
        logger.debug(f"Library path: {lib_path}")
    except RuntimeError:
        logger.exception("Error retrieving ChEMBL paths")
