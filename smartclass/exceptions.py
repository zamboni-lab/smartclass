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
    """Truncate text with ellipsis if it exceeds *max_length* characters.

    Parameters
    ----------
    text : str
        Input text to truncate.
    max_length : int
        Maximum allowed length before truncation. Default is _MAX_SNIPPET_LENGTH.

    Returns
    -------
    str
        Original text if within limit, otherwise truncated text ending with ``...``.
    """
    if len(text) <= max_length:
        return text
    return f"{text[:max_length]}..."


class SmartclassError(Exception):
    """Base exception for all smartclass errors.

    All custom exceptions in smartclass inherit from this class,
    making it easy to catch all smartclass-specific errors.
    """

    def __init__(self, message: str, *args: Any, **kwargs: Any) -> None:
        """Initialize with a descriptive message.

        Parameters
        ----------
        message : str
            Human-readable error description.
        *args : Any
            Additional positional arguments forwarded to the base ``Exception``.
        **kwargs : Any
            Additional keyword arguments forwarded to the base ``Exception``.
        """
        self.message = message
        super().__init__(message, *args, **kwargs)

    def __str__(self) -> str:
        """Return the human-readable error message.

        Returns
        -------
        str
            The message passed at construction time.
        """
        return self.message


# Chemical Conversion Errors
class ChemicalConversionError(SmartclassError):
    """Error during chemical format conversion.

    Raised when conversion between chemical representations fails,
    e.g. SMILES to Mol or InChI to SMILES.
    """

    pass


class MoleculeParsingError(ChemicalConversionError):
    """Error parsing a molecule representation.

    Raised when RDKit fails to parse a molecular structure.
    """

    def __init__(
        self,
        input_string: str,
        input_type: str = "unknown",
        *args: Any,
        **kwargs: Any,
    ) -> None:
        """Initialize with details about the failed parsing.

        Parameters
        ----------
        input_string : str
            The string that failed to parse.
        input_type : str
            Format of the input string (e.g. SMILES, InChI, SMARTS). Default is 'unknown'.
        *args : Any
            Additional positional arguments forwarded to the parent exception.
        **kwargs : Any
            Additional keyword arguments forwarded to the parent exception.
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
    """Error during chemical classification."""

    pass


# Configuration Errors
class ConfigurationError(SmartclassError):
    """Error in configuration settings."""

    pass


# Data I/O Errors
class DataLoadingError(SmartclassError):
    """Error loading data from a file or URL."""

    def __init__(self, source: str, reason: str | None = None, *args, **kwargs) -> None:
        """Initialize with details about the failed loading.

        Parameters
        ----------
        source : str
            Path or URL that failed to load.
        reason : str | None
            Optional human-readable reason for the failure. Default is None.
        *args : Any
            Additional positional arguments forwarded to the base ``Exception``.
        **kwargs : Any
            Additional keyword arguments forwarded to the base ``Exception``.
        """
        self.source = source
        self.reason = reason
        message = f"Failed to load data from '{source}'"
        if reason:
            message += f": {reason}"
        super().__init__(message, *args, **kwargs)


class DataExportError(SmartclassError):
    """Error exporting data to a file."""

    def __init__(
        self,
        destination: str,
        reason: str | None = None,
        *args,
        **kwargs,
    ) -> None:
        """Initialize with details about the failed export.

        Parameters
        ----------
        destination : str
            Output path that failed to write.
        reason : str | None
            Optional human-readable reason for the failure. Default is None.
        *args : Any
            Additional positional arguments forwarded to the base ``Exception``.
        **kwargs : Any
            Additional keyword arguments forwarded to the base ``Exception``.
        """
        self.destination = destination
        self.reason = reason
        message = f"Failed to export data to '{destination}'"
        if reason:
            message += f": {reason}"
        super().__init__(message, *args, **kwargs)


# Input Validation Errors
class InvalidInputError(SmartclassError):
    """Error for invalid user input."""

    def __init__(
        self,
        parameter: str,
        value: object,
        reason: str | None = None,
        *args,
        **kwargs,
    ) -> None:
        """Initialize with details about the invalid input.

        Parameters
        ----------
        parameter : str
            Name of the parameter with invalid value.
        value : object
            The invalid value provided.
        reason : str | None
            Optional human-readable reason why the value is invalid. Default is None.
        *args : Any
            Additional positional arguments forwarded to the base ``Exception``.
        **kwargs : Any
            Additional keyword arguments forwarded to the base ``Exception``.
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
    """Error during network operations."""

    def __init__(
        self,
        url: str,
        status_code: int | None = None,
        reason: str | None = None,
        *args,
        **kwargs,
    ) -> None:
        """Initialize with details about the network failure.

        Parameters
        ----------
        url : str
            URL that failed.
        status_code : int | None
            HTTP status code if available. Default is None.
        reason : str | None
            Optional human-readable reason for the failure. Default is None.
        *args : Any
            Additional positional arguments forwarded to the base ``Exception``.
        **kwargs : Any
            Additional keyword arguments forwarded to the base ``Exception``.
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
