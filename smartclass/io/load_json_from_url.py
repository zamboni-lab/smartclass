"""Load json from url."""

from __future__ import annotations

import json
import logging
from typing import Any

import pooch


__all__ = [
    "load_json_from_url",
]

# Create a pooch instance for JSON loading
JSON_DOWNLOADER = pooch.create(
    path=pooch.os_cache("json_cache"),
    base_url="",
    registry={},
)


def load_json_from_url(url: str) -> dict[str, Any] | None:
    """
    Load JSON from URL.

    Args:
        url: URL of the JSON file.

    Returns:
        Parsed JSON data as a dictionary or None if loading fails.
    """
    try:
        fetcher = pooch.HTTPDownloader(timeout=60)
        data_file = JSON_DOWNLOADER.fetch(url, downloader=fetcher)

        with open(data_file) as f:
            data = json.load(f)
        logging.debug(f"Got {url} successfully.")
        return data
    except Exception as e:
        logging.debug(f"Failed to get {url}: {e}")
        return None
