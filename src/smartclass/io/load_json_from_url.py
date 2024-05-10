"""Load json from url."""

from __future__ import annotations

import json
import logging

import requests

__all__ = [
    "load_json_from_url",
]


def load_json_from_url(url: str) -> dict | None:
    """
    Load json from url.

    :param url: url of the file.
    :type url: str

    :returns: A dictionary.
    :rtype: Union[dict, None]
    """
    response = requests.get(url, timeout=60)
    if response.status_code == 200:
        data = json.loads(response.text)
        logging.debug(f"Got {url} successfully.")
    else:
        logging.debug(f"Failed to get {url}.")
        data = None
    return data
