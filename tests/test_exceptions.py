"""Tests for custom exceptions."""

from __future__ import annotations

import unittest

from smartclass.exceptions import (
    ChemicalConversionError,
    ClassificationError,
    ConfigurationError,
    DataExportError,
    DataLoadingError,
    InChIError,
    InvalidInputError,
    MoleculeParsingError,
    NetworkError,
    SMARTSError,
    SMILESError,
    SmartclassError,
)


class TestExceptionHierarchy(unittest.TestCase):
    """Test exception class hierarchy."""

    def test_all_inherit_from_smartclass_error(self):
        """Test all custom exceptions inherit from SmartclassError."""
        exceptions = [
            ChemicalConversionError("test"),
            ClassificationError("test"),
            ConfigurationError("test"),
            DataLoadingError("source"),
            DataExportError("dest"),
            InvalidInputError("param", "value"),
            NetworkError("url"),
            MoleculeParsingError("input", "type"),
            SMILESError("smiles"),
            SMARTSError("smarts"),
            InChIError("inchi"),
        ]

        for exc in exceptions:
            self.assertIsInstance(exc, SmartclassError)

    def test_molecule_parsing_errors_inherit_from_conversion(self):
        """Test parsing errors inherit from ChemicalConversionError."""
        exceptions = [
            SMILESError("test"),
            SMARTSError("test"),
            InChIError("test"),
        ]

        for exc in exceptions:
            self.assertIsInstance(exc, ChemicalConversionError)


class TestSmartclassError(unittest.TestCase):
    """Test SmartclassError base class."""

    def test_message_attribute(self):
        """Test error message is accessible."""
        exc = SmartclassError("Test error message")

        self.assertEqual(exc.message, "Test error message")
        self.assertEqual(str(exc), "Test error message")


class TestDataLoadingError(unittest.TestCase):
    """Test DataLoadingError."""

    def test_basic_message(self):
        """Test basic error message."""
        exc = DataLoadingError("/path/to/file")

        self.assertIn("/path/to/file", str(exc))

    def test_with_reason(self):
        """Test error message with reason."""
        exc = DataLoadingError("/path/to/file", reason="File not found")

        self.assertIn("/path/to/file", str(exc))
        self.assertIn("File not found", str(exc))


class TestNetworkError(unittest.TestCase):
    """Test NetworkError."""

    def test_basic_message(self):
        """Test basic error message."""
        exc = NetworkError("https://example.com")

        self.assertIn("https://example.com", str(exc))

    def test_with_status_code(self):
        """Test error with status code."""
        exc = NetworkError("https://example.com", status_code=503)

        self.assertIn("503", str(exc))

    def test_with_reason(self):
        """Test error with reason."""
        exc = NetworkError("https://example.com", reason="Connection refused")

        self.assertIn("Connection refused", str(exc))


class TestSMILESError(unittest.TestCase):
    """Test SMILESError."""

    def test_short_smiles(self):
        """Test error with short SMILES."""
        exc = SMILESError("CCO")

        self.assertIn("CCO", str(exc))
        self.assertIn("SMILES", str(exc))

    def test_long_smiles_truncated(self):
        """Test long SMILES are truncated in message."""
        long_smiles = "C" * 100
        exc = SMILESError(long_smiles)

        # Should truncate to 50 chars + "..."
        self.assertIn("...", str(exc))


class TestInvalidInputError(unittest.TestCase):
    """Test InvalidInputError."""

    def test_basic_message(self):
        """Test basic error message."""
        exc = InvalidInputError("threshold", -1)

        self.assertIn("threshold", str(exc))
        self.assertIn("-1", str(exc))

    def test_with_reason(self):
        """Test error with reason."""
        exc = InvalidInputError("threshold", -1, reason="must be positive")

        self.assertIn("must be positive", str(exc))


if __name__ == "__main__":
    unittest.main()
