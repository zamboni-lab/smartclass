"""Tests for the logging module."""

from __future__ import annotations

import logging
import unittest

from smartclass import logging as smartclass_logging
from smartclass.logging import configure_logging, get_logger


class TestGetLogger(unittest.TestCase):
    """Test cases for get_logger."""

    def test_returns_logger(self):
        """get_logger returns a logging.Logger instance."""
        logger = get_logger("smartclass.tests")
        self.assertIsInstance(logger, logging.Logger)

    def test_child_of_smartclass_logger(self):
        """Loggers are children of the 'smartclass' logger."""
        logger = get_logger("smartclass.chem.example")
        self.assertTrue(logger.name.startswith("smartclass"))

    def test_same_name_same_logger(self):
        """Requesting the same name returns the same logger."""
        self.assertIs(get_logger("smartclass.tests"), get_logger("smartclass.tests"))


class TestConfigureLogging(unittest.TestCase):
    """Test cases for configure_logging."""

    def tearDown(self):
        """Reset the package logger handlers after each test."""
        smartclass_logging._logging_configured = False
        root = logging.getLogger("smartclass")
        for handler in root.handlers[:]:
            root.removeHandler(handler)

    def test_configures_smartclass_logger_level(self):
        """configure_logging sets the level on the 'smartclass' logger."""
        configure_logging(level=logging.DEBUG, force=True)

        self.assertEqual(logging.getLogger("smartclass").level, logging.DEBUG)

    def test_accepts_string_level(self):
        """configure_logging accepts a string level."""
        configure_logging(level="INFO", force=True)

        self.assertEqual(logging.getLogger("smartclass").level, logging.INFO)

    def test_does_not_propagate(self):
        """The 'smartclass' logger does not propagate to the root logger."""
        configure_logging(force=True)

        self.assertFalse(logging.getLogger("smartclass").propagate)

    def test_ignored_without_force(self):
        """A second call without force is ignored."""
        configure_logging(level=logging.INFO, force=True)
        configure_logging(level=logging.ERROR, force=False)

        self.assertEqual(logging.getLogger("smartclass").level, logging.INFO)

    def test_force_reconfigures(self):
        """A second call with force reconfigures the level."""
        configure_logging(level=logging.INFO, force=True)
        configure_logging(level=logging.ERROR, force=True)

        self.assertEqual(logging.getLogger("smartclass").level, logging.ERROR)


if __name__ == "__main__":
    unittest.main()
