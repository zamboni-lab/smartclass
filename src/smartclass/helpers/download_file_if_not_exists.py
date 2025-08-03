"""Download file if not exists."""

from __future__ import annotations

import logging
import os

import pooch

__all__ = [
    "download_file_if_not_exists",
]

# Create a pooch instance for file downloads
FILE_DOWNLOADER = pooch.create(
    path=pooch.os_cache("file_downloads"),
    base_url="",
    registry={},
)


def download_file_if_not_exists(url: str, output: str) -> None:
    """
    Downloads a file from the specified URL if it does not exist.

    Args:
        url: Input URL.
        output: Output file path.
    """
    if os.path.exists(output):
        logging.info(f"File already exists at {output}")
        return

    try:
        fetcher = pooch.HTTPDownloader(timeout=60)
        downloaded_file = FILE_DOWNLOADER.fetch(url, downloader=fetcher)

        # Copy the downloaded file to the desired output location
        if downloaded_file != output:
            with open(downloaded_file, "rb") as src, open(output, "wb") as dst:
                dst.write(src.read())

        logging.info(f"File downloaded and saved to {output}")
    except Exception as e:
        logging.exception(f"Failed to download the file: {e}")
