"""Send a GET request and retrieve JSON data."""

from __future__ import annotations

import logging

import requests

__all__ = [
    "get_request",
]


def get_request(url: str, query: str) -> list[dict]:
    """
    Send a GET request and retrieve JSON data.

    :param url: The URL to send the GET request to.
    :type url: str

    :param query: The query string to include in the request.
    :type query: str

    :raises requests.exceptions.RequestException: If there is an error.

    :returns: A list of dictionaries containing the retrieved JSON data.
    :rtype: list[dict]
    """
    try:
        params = {"format": "json", "query": query}
        response = requests.get(url, timeout=60, params=params)
        response.raise_for_status()  # Raise an exception for non-2xx status codes

        data = response.json()
        bindings = data["results"]["bindings"]

        results = []
        for binding in bindings:
            result = {key: value["value"] for key, value in binding.items()}
            results.append(result)

        return results
    except requests.exceptions.RequestException as e:
        logging.error(f"Error making the GET request: {e}")
        raise


# Example usage:
if __name__ == "__main__":
    try:
        url = "https://example.com/sparql-endpoint"
        query = "SELECT * WHERE {?s ?p ?o}"
        results = get_request(url, query)
    except requests.exceptions.RequestException as e:
        logging.error(f"Request error: {e}")
