"""Custom exceptions for smartclass.

This module defines a hierarchy of exceptions specific to the smartclass
package, enabling more precise error handling and clearer error messages.

Exception Hierarchy:
    SmartclassError (base)
    ├── ChemicalConversionError
    │   └── MoleculeParsingError
    │       ├── SMILESError
    │       ├── SMARTSError
    │       └── InChIError
    ├── ClassificationError
    ├── ConfigurationError
    ├── DataLoadingError
    ├── DataExportError
    ├── InvalidInputError
    └── NetworkError
"""

from __future__ import annotations

from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from typing import Any

__all__ = [
    "SmartclassError",
    "ChemicalConversionError",
    "ClassificationError",
    "ConfigurationError",
    "DataLoadingError",
    "DataExportError",
    "InvalidInputError",
    "NetworkError",
    "MoleculeParsingError",
    "SMARTSError",
    "SMILESError",
    "InChIError",
]

# Maximum length for error message snippets
_MAX_SNIPPET_LENGTH = 1500


def _truncate(text: str, max_length: int = _MAX_SNIPPET_LENGTH) -> str:
    """Truncate text with ellipsis if too long."""
    if len(text) <= max_length:
        return text
    return f"{text[:max_length]}..."


class SmartclassError(Exception):
    """Base exception for all smartclass errors.

    All custom exceptions in smartclass inherit from this class,
    making it easy to catch all smartclass-specific errors.

    Attributes:
        message: Human-readable error description.
    """

    def __init__(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Initialize with a descriptive message.

        Args:
            message: Human-readable error description.
        """
        self.message = message
        super().__init__(message, *args, **kwargs)

    def __str__(self) -> str:
        """Return the error message."""
        return self.message


# Chemical Conversion Errors
class ChemicalConversionError(SmartclassError):
    """Error during chemical format conversion.

    Raised when conversion between chemical representations fails
    (e.g., SMILES to Mol, InChI to SMILES).
    """

    pass


class MoleculeParsingError(ChemicalConversionError):
    """Error parsing a molecule representation.

    Raised when RDKit fails to parse a molecular structure.

    Attributes:
        input_string: The string that failed to parse.
        input_type: Type of input (SMILES, InChI, SMARTS, etc.).
    """

    def __init__(
        self,
        input_string: str,
        input_type: str = "unknown",
        *args: Any,
        **kwargs: Any,
    ) -> None:
        """Initialize with details about the failed parsing.

        Args:
            input_string: The string that failed to parse.
            input_type: Type of input (SMILES, InChI, SMARTS, etc.).
        """
        self.input_string = input_string
        self.input_type = input_type
        message = f"Failed to parse {input_type}: '{_truncate(input_string)}'"
        super().__init__(message, *args, **kwargs)


class SMILESError(MoleculeParsingError):
    """Error parsing a SMILES string."""

    def __init__(self, smiles: str, *args, **kwargs) -> None:
        super().__init__(smiles, "SMILES", *args, **kwargs)


class SMARTSError(MoleculeParsingError):
    """Error parsing a SMARTS pattern."""

    def __init__(self, smarts: str, *args, **kwargs) -> None:
        super().__init__(smarts, "SMARTS", *args, **kwargs)


class InChIError(MoleculeParsingError):
    """Error parsing an InChI string."""

    def __init__(self, inchi: str, *args, **kwargs) -> None:
        super().__init__(inchi, "InChI", *args, **kwargs)


# Classification Errors
class ClassificationError(SmartclassError):
    """Error during chemical classification.

    Raised when the classification process encounters an issue.
    """

    pass


# Configuration Errors
class ConfigurationError(SmartclassError):
    """Error in configuration settings.

    Raised when configuration is invalid or missing required values.
    """

    pass


# Data I/O Errors
class DataLoadingError(SmartclassError):
    """Error loading data from a file or URL.

    Raised when data cannot be loaded or parsed correctly.
    """

    def __init__(self, source: str, reason: str | None = None, *args, **kwargs) -> None:
        """Initialize with details about the failed loading.

        :param source: Path or URL that failed to load.
        :param reason: Optional reason for the failure.
        """
        self.source = source
        self.reason = reason
        message = f"Failed to load data from '{source}'"
        if reason:
            message += f": {reason}"
        super().__init__(message, *args, **kwargs)


class DataExportError(SmartclassError):
    """Error exporting data to a file.

    Raised when data cannot be written to the specified destination.
    """

    def __init__(
        self,
        destination: str,
        reason: str | None = None,
        *args,
        **kwargs,
    ) -> None:
        """Initialize with details about the failed export.

        :param destination: Path that failed to write.
        :param reason: Optional reason for the failure.
        """
        self.destination = destination
        self.reason = reason
        message = f"Failed to export data to '{destination}'"
        if reason:
            message += f": {reason}"
        super().__init__(message, *args, **kwargs)


# Input Validation Errors
class InvalidInputError(SmartclassError):
    """Error for invalid user input.

    Raised when input validation fails.
    """

    def __init__(
        self,
        parameter: str,
        value: object,
        reason: str | None = None,
        *args,
        **kwargs,
    ) -> None:
        """Initialize with details about the invalid input.

        :param parameter: Name of the parameter with invalid value.
        :param value: The invalid value provided.
        :param reason: Optional reason why the value is invalid.
        """
        self.parameter = parameter
        self.value = value
        self.reason = reason
        message = f"Invalid value for '{parameter}': {value!r}"
        if reason:
            message += f" ({reason})"
        super().__init__(message, *args, **kwargs)


# Network Errors
class NetworkError(SmartclassError):
    """Error during network operations.

    Raised when HTTP requests fail after retries.
    """

    def __init__(
        self,
        url: str,
        status_code: int | None = None,
        reason: str | None = None,
        *args,
        **kwargs,
    ) -> None:
        """Initialize with details about the network failure.

        :param url: URL that failed.
        :param status_code: HTTP status code if available.
        :param reason: Optional reason for the failure.
        """
        self.url = url
        self.status_code = status_code
        self.reason = reason
        message = f"Network request to '{url}' failed"
        if status_code:
            message += f" (status {status_code})"
        if reason:
            message += f": {reason}"
        super().__init__(message, *args, **kwargs)
