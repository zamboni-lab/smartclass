"""Send a GET request to a SPARQL endpoint and retrieve JSON data."""

from __future__ import annotations

import logging
import random
import time

import requests
from requests.exceptions import RequestException

__all__ = ["get_request"]

# Configure module-level logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def get_request(
    url: str,
    query: str,
    max_retries: int = 3,
    base_delay: float = 2.0,
    timeout: int = 60,
) -> list[dict[str, str]]:
    """
    Send a GET request to a SPARQL endpoint and retrieve JSON data.

    :param url: The SPARQL endpoint URL.
    :param query: The SPARQL query string.
    :param max_retries: Maximum number of retry attempts.
    :param base_delay: Base delay (in seconds) for retry backoff.
    :param timeout: Timeout for the request in seconds.
    :return: A list of dictionaries representing the query results.
    """
    headers = {
        "User-Agent": "SmartClassBot/1.0 (mailto:your_email@example.com)",
        "Accept": "application/sparql-results+json",
    }

    params = {
        "format": "json",
        "query": query,
    }

    attempt = 0
    qlever_url = "https://qlever.cs.uni-freiburg.de/api/wikidata"
    while attempt < max_retries:
        try:
            response = requests.get(url, headers=headers, params=params, timeout=timeout)
            response.raise_for_status()

            data = response.json()
            bindings = data.get("results", {}).get("bindings", [])

            results = [{key: value["value"] for key, value in binding.items()} for binding in bindings]

            logger.info(f"Query successful: {len(results)} results retrieved.")
            return results

        except RequestException as e:
            status_code = getattr(response, "status_code", None)
            retriable = status_code in {429, 503}

            if retriable and attempt < max_retries - 1:
                wait_time = base_delay * (2**attempt) + random.uniform(0, 1)
                logger.warning(f"Request failed with status {status_code}. Retrying in {wait_time:.1f} seconds...")
                time.sleep(wait_time)
                attempt += 1
            elif url != qlever_url:
                logger.warning(
                    f"Request failed with status {status_code} on WDQS. Retrying one last time on QLever endpoint...",
                )
                url = qlever_url
            else:
                logger.error(f"Request failed: {e}")
                raise RuntimeError(f"Failed to retrieve data from {url} after {max_retries} attempts.") from e

    # If somehow exits loop without success
    return []
