"""Load json from path."""

from __future__ import annotations

import json


__all__ = [
    "load_json_from_path",
]


def load_json_from_path(path: str) -> dict:
    """
    Load json from path.

    :param path: Path of the file.
    :type path: str

    :returns: A dictionary.
    :rtype: dict
    """
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return data
