"""Download file if not exists."""

from __future__ import annotations

import logging
import os

import requests

__all__ = [
    "download_file_if_not_exists",
]


def download_file_if_not_exists(url: str, output: str):
    """
    Downloads a file from the specified URL if it does not already exist at the output path.

    :param url: Input.
    :type url: str

    :param output: Output.
    :type output: str
    """
    if os.path.exists(output):
        logging.info(f"File already exists at {output}")
        return
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        with open(output, "wb") as file:
            file.write(response.content)
        logging.info(f"File downloaded and saved to {output}")
    except requests.RequestException as e:
        logging.error(f"Failed to download the file: {e}")
