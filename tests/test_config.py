"""Tests for the configuration module."""

from __future__ import annotations

import os
import unittest
from pathlib import Path
from unittest.mock import patch

from smartclass.config import Config, get_config, reset_config


class TestConfig(unittest.TestCase):
    """Test cases for Config class."""

    def setUp(self):
        """Reset config before each test."""
        reset_config()

    def tearDown(self):
        """Clean up after each test."""
        reset_config()

    def test_default_values(self):
        """Test that Config has sensible defaults."""
        config = Config()

        self.assertEqual(config.wikidata_endpoint, "https://query.wikidata.org/sparql")
        self.assertEqual(config.http_timeout, 60)
        self.assertEqual(config.http_max_retries, 3)
        self.assertEqual(config.chembl_fp_length, 2048)
        self.assertEqual(config.mcs_threshold, 0.7)
        self.assertEqual(config.log_level, "INFO")

    def test_get_config_singleton(self):
        """Test that get_config returns the same instance."""
        config1 = get_config()
        config2 = get_config()

        self.assertIs(config1, config2)

    def test_reset_config(self):
        """Test that reset_config clears the singleton."""
        config1 = get_config()
        reset_config()
        config2 = get_config()

        self.assertIsNot(config1, config2)

    @patch.dict(os.environ, {"SMARTCLASS_HTTP_TIMEOUT": "120"})
    def test_env_var_override(self):
        """Test that environment variables override defaults."""
        config = Config()

        self.assertEqual(config.http_timeout, 120)

    @patch.dict(os.environ, {"SMARTCLASS_LOG_LEVEL": "DEBUG"})
    def test_env_var_log_level(self):
        """Test log level can be set via environment."""
        config = Config()

        self.assertEqual(config.log_level, "DEBUG")

    def test_output_path_helper(self):
        """Test get_output_path returns correct path."""
        config = Config()

        path = config.get_output_path("test.tsv")

        self.assertIsInstance(path, Path)
        self.assertTrue(str(path).endswith("test.tsv"))

    def test_cache_path_helper(self):
        """Test get_cache_path returns correct path."""
        config = Config()

        path = config.get_cache_path("data.pkl")

        self.assertIsInstance(path, Path)
        self.assertTrue(str(path).endswith("data.pkl"))


if __name__ == "__main__":
    unittest.main()
