"""Download file if not exists."""

from __future__ import annotations

import os
import shutil

import pooch

from smartclass.logging import get_logger

__all__ = [
    "download_file_if_not_exists",
]

logger = get_logger(__name__)

# Create a pooch instance for file downloads
FILE_DOWNLOADER = pooch.create(
    path=pooch.os_cache("file_downloads"),
    base_url="",
    registry={},
)


def download_file_if_not_exists(url: str, output: str) -> None:
    """Downloads a file from the specified URL if it does not exist.

    Parameters
    ----------
    url : str
        URL.
    output : str
        Output file path.
    """
    if os.path.exists(output):
        logger.info(f"File already exists at {output}")
        return

    try:
        fetcher = pooch.HTTPDownloader(timeout=60)
        downloaded_file = FILE_DOWNLOADER.fetch(url, downloader=fetcher)

        # Copy the downloaded file to the desired output location
        if downloaded_file != output:
            shutil.copyfile(downloaded_file, output)

        logger.info(f"File downloaded and saved to {output}")
    except (OSError, ValueError):
        logger.exception("Failed to download the file")
