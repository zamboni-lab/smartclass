"""Send a GET request and retrieve JSON data."""

from __future__ import annotations

import logging
import time
from typing import List, Dict

import requests
from requests.exceptions import RequestException

__all__ = [
    "get_request",
]

# Configure logging
logging.basicConfig(level=logging.INFO)

def get_request(url: str, query: str, max_retries: int = 3, retry_delay: int = 60) -> List[Dict[str, str]]:
    """
    Send a GET request and retrieve JSON data with retry logic for handling rate limits.

    :param url: The URL to send the GET request to.
    :type url: str

    :param query: The query string to include in the request.
    :type query: str

    :param max_retries: The maximum number of retry attempts on rate limit errors.
    :type max_retries: int

    :param retry_delay: The delay between retry attempts in seconds.
    :type retry_delay: int

    :raises RequestException: If there is an error after retry attempts.

    :returns: A list of dictionaries containing the retrieved JSON data.
    :rtype: list[dict]
    """
    for attempt in range(max_retries):
        try:
            params = {"format": "json", "query": query}
            response = requests.get(url, timeout=60, params=params)
            response.raise_for_status()  # Raise an exception for non-2xx status codes

            data = response.json()
            bindings = data.get("results", {}).get("bindings", [])

            results = [{key: value["value"] for key, value in binding.items()} for binding in bindings]

            return results

        except RequestException as e:
            if response.status_code == 429 and attempt < max_retries - 1:
                logging.warning(f"Rate limit exceeded. Retrying after {retry_delay} seconds...")
                time.sleep(retry_delay)
            else:
                logging.error(f"Error making the GET request: {e}")
                raise

# Example usage:
if __name__ == "__main__":
    try:
        url = "https://query.wikidata.org/sparql"
        query = "SELECT ?item ?itemLabel WHERE { ?item wdt:P31 wd:Q5. SERVICE wikibase:label { bd:serviceParam wikibase:language '[AUTO_LANGUAGE],en'. } } LIMIT 10"
        results = get_request(url, query)
        for result in results:
            print(result)
    except RequestException as e:
        logging.error(f"Request error: {e}")
