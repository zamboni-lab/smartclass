"""Read SPARQL queries from files or URLs."""

from __future__ import annotations

from pathlib import Path

import requests

from smartclass.exceptions import DataLoadingError
from smartclass.logging import get_logger

__all__ = [
    "read_query",
]

logger = get_logger(__name__)


def read_query(query: str | Path) -> str:
    """
    Read a SPARQL query from a file path or URL.

    :param query: Path to local file or URL to fetch query from.
    :returns: The SPARQL query string.
    :raises DataLoadingError: If the query cannot be loaded.
    """
    query_str = str(query)

    if query_str.startswith(("http://", "https://")):
        try:
            response = requests.get(query_str, timeout=30)
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            raise DataLoadingError(query_str, reason=str(e)) from e
    else:
        query_path = Path(query)
        if not query_path.exists():
            raise DataLoadingError(str(query_path), reason="File not found")
        try:
            return query_path.read_text(encoding="utf-8")
        except OSError as e:
            raise DataLoadingError(str(query_path), reason=str(e)) from e
