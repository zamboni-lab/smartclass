"""Load json from url."""

from __future__ import annotations

import json
from typing import Any

import pooch

from smartclass.logging import get_logger

__all__ = [
    "load_json_from_url",
]

logger = get_logger(__name__)

# Create a pooch instance for JSON loading
JSON_DOWNLOADER = pooch.create(
    path=pooch.os_cache("json_cache"),
    base_url="",
    registry={},
)


def load_json_from_url(url: str) -> dict[str, Any] | None:
    """Load JSON from URL.

    Parameters
    ----------
    url : str
        JSON file.

    Returns
    -------
    dict[str, Any] | None
        None if loading fails.
    """
    try:
        fetcher = pooch.HTTPDownloader(timeout=60)
        data_file = JSON_DOWNLOADER.fetch(url, downloader=fetcher)

        with open(data_file) as f:
            data = json.load(f)
        logger.debug(f"Got {url} successfully.")
        return data
    except (OSError, ValueError) as e:
        logger.debug(f"Failed to get {url}: {e}")
        return None
