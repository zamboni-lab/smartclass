"""Load json from URL or path."""

from __future__ import annotations

import os

from smartclass.io.export_dict_to_json import export_dict_to_json
from smartclass.io.load_json_from_path import load_json_from_path
from smartclass.io.load_json_from_url import load_json_from_url
from smartclass.logging import get_logger

__all__ = [
    "load_json_from_url_or_path",
]

logger = get_logger(__name__)


def load_json_from_url_or_path(url: str, name: str) -> dict | None:
    """Load json from URL or path.

    Parameters
    ----------
    url : str
        URL to get the json from.
    name : str
        Name of the file.

    Returns
    -------
    dict | None
        The loaded JSON data, or ``None`` if loading failed.
    """
    file_name = f"{name}.json"
    file_path = os.path.join("scratch", file_name)

    if os.path.exists(file_path):
        logger.debug(f"JSON file {file_path} already exists. Loading from the file.")
        return load_json_from_path(path=file_path)
    else:
        data = load_json_from_url(url=url)
        if data is not None:
            export_dict_to_json(dic=data, output=file_path)
        return data
