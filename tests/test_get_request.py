"""Tests for the SPARQL GET request helper."""

from __future__ import annotations

import os
import unittest
from unittest.mock import patch

import requests_mock

from smartclass.config import reset_config
from smartclass.exceptions import NetworkError
from smartclass.helpers.get_request import get_request

URL = "https://example.org/sparql"
QLEVER_URL = "https://qlever.example.org/api/wikidata"
QUERY = "SELECT * WHERE { ?s ?p ?o } LIMIT 1"

_OK_BODY = {"results": {"bindings": [{"x": {"value": "1"}}, {"x": {"value": "2"}}]}}


class TestGetRequest(unittest.TestCase):
    """Test cases for get_request."""

    def tearDown(self):
        """Reset the global config so env overrides do not leak between tests."""
        reset_config()

    def test_successful_request_returns_results(self):
        """A successful request returns the parsed bindings."""
        with requests_mock.Mocker() as mock:
            mock.get(URL, json=_OK_BODY)

            results = get_request(URL, QUERY)

        self.assertEqual(results, [{"x": "1"}, {"x": "2"}])

    def test_retries_on_retriable_status(self):
        """A 503 response is retried and a later success is returned."""
        with (
            requests_mock.Mocker() as mock,
            patch("smartclass.helpers.get_request.random.uniform", return_value=0),
        ):
            mock.get(URL, [{"status_code": 503}, {"json": _OK_BODY}])

            results = get_request(URL, QUERY, max_retries=3, base_delay=0)

        self.assertEqual(results, [{"x": "1"}, {"x": "2"}])

    def test_falls_back_to_qlever_endpoint(self):
        """A non-retriable failure on the primary endpoint falls back to QLever."""
        with (
            requests_mock.Mocker() as mock,
            patch.dict(os.environ, {"SMARTCLASS_QLEVER_ENDPOINT": QLEVER_URL}),
        ):
            reset_config()
            mock.get(URL, status_code=404)
            mock.get(QLEVER_URL, json=_OK_BODY)

            results = get_request(URL, QUERY, max_retries=1)

        self.assertEqual(results, [{"x": "1"}, {"x": "2"}])

    def test_raises_network_error_when_all_endpoints_fail(self):
        """A NetworkError is raised when every endpoint has failed."""
        with (
            requests_mock.Mocker() as mock,
            patch.dict(os.environ, {"SMARTCLASS_QLEVER_ENDPOINT": QLEVER_URL}),
            patch("smartclass.helpers.get_request.random.uniform", return_value=0),
        ):
            reset_config()
            mock.get(URL, status_code=503)
            mock.get(QLEVER_URL, status_code=503)

            with self.assertRaises(NetworkError):
                get_request(URL, QUERY, max_retries=1, base_delay=0)


if __name__ == "__main__":
    unittest.main()
