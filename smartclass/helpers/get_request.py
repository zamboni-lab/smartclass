"""Send a GET request to a SPARQL endpoint and retrieve JSON data."""

from __future__ import annotations

import random
import time

import requests
from requests.exceptions import RequestException

from smartclass.config import get_config
from smartclass.exceptions import NetworkError
from smartclass.logging import get_logger

__all__ = ["get_request"]

logger = get_logger(__name__)


def get_request(
    url: str,
    query: str,
    max_retries: int | None = None,
    base_delay: float | None = None,
    timeout: int | None = None,
) -> list[dict[str, str]]:
    """
    Send a GET request to a SPARQL endpoint and retrieve JSON data.

    Uses exponential backoff with jitter for retries. Falls back to QLever
    endpoint if the primary endpoint fails.

    :param url: The SPARQL endpoint URL.
    :param query: The SPARQL query string.
    :param max_retries: Maximum retry attempts. Uses config default if None.
    :param base_delay: Base delay (seconds) for backoff. Uses config default if None.
    :param timeout: Request timeout in seconds. Uses config default if None.
    :returns: List of dictionaries representing query results.
    :raises NetworkError: If the request fails after all retries.
    """
    config = get_config()

    # Use config defaults if not specified
    max_retries = max_retries if max_retries is not None else config.http_max_retries
    base_delay = base_delay if base_delay is not None else config.http_base_delay
    timeout = timeout if timeout is not None else config.http_timeout

    headers = {
        "User-Agent": config.user_agent,
        "Accept": "application/sparql-results+json",
    }

    params = {
        "format": "json",
        "query": query,
    }

    attempt = 0
    response = None  # Initialize to avoid UnboundLocalError

    while attempt < max_retries:
        try:
            response = requests.get(
                url,
                headers=headers,
                params=params,
                timeout=timeout,
            )
            response.raise_for_status()

            data = response.json()
            bindings = data.get("results", {}).get("bindings", [])

            results = [
                {key: value["value"] for key, value in binding.items()}
                for binding in bindings
            ]

            logger.info(f"Query successful: {len(results)} results retrieved.")
            return results

        except RequestException as e:
            status_code = getattr(response, "status_code", None) if response else None
            retriable_codes = {429, 503}

            if status_code in retriable_codes and attempt < max_retries - 1:
                wait_time = base_delay * (2**attempt) + random.uniform(0, 1)
                logger.warning(
                    f"Request failed with status {status_code}. "
                    f"Retrying in {wait_time:.1f} seconds (attempt {attempt + 1}/{max_retries})..."
                )
                time.sleep(wait_time)
                attempt += 1
            elif url != config.qlever_endpoint:
                logger.warning(
                    f"Request failed with status {status_code}. "
                    "Trying fallback QLever endpoint..."
                )
                url = config.qlever_endpoint
                attempt = 0  # Reset attempts for fallback
            else:
                raise NetworkError(
                    url=url,
                    status_code=status_code,
                    reason=str(e),
                ) from e

    # Should not reach here, but return empty list as safeguard
    return []
