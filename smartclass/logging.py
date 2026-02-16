"""Centralized logging configuration for smartclass.

This module provides a single point for configuring logging across the entire
smartclass package. Import and use `get_logger(__name__)` in all modules.
"""

from __future__ import annotations

import logging
import sys
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from pathlib import Path

__all__ = [
    "configure_logging",
    "get_logger",
]

# Track whether logging has been configured
_logging_configured = False

# Default logging format
DEFAULT_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def configure_logging(
    level: int | str = logging.INFO,
    format_string: str = DEFAULT_FORMAT,
    date_format: str = DEFAULT_DATE_FORMAT,
    log_file: Path | str | None = None,
    force: bool = False,
) -> None:
    """
    Configure logging for the smartclass package.

    This should be called once at application startup. Subsequent calls
    will be ignored unless `force=True`.

    :param level: Logging level (e.g., logging.DEBUG, logging.INFO, "DEBUG", "INFO").
    :param format_string: Format string for log messages.
    :param date_format: Format string for timestamps.
    :param log_file: Optional path to a log file. If provided, logs will be
        written to both console and file.
    :param force: If True, reconfigure logging even if already configured.
    """
    global _logging_configured

    if _logging_configured and not force:
        return

    # Convert string level to int if necessary
    if isinstance(level, str):
        level = getattr(logging, level.upper(), logging.INFO)

    # Create handlers
    handlers: list[logging.Handler] = []

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(logging.Formatter(format_string, date_format))
    handlers.append(console_handler)

    # File handler (optional)
    if log_file is not None:
        file_handler = logging.FileHandler(str(log_file))
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter(format_string, date_format))
        handlers.append(file_handler)

    # Configure root logger for smartclass
    root_logger = logging.getLogger("smartclass")
    root_logger.setLevel(level)

    # Remove existing handlers to avoid duplicates
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    for handler in handlers:
        root_logger.addHandler(handler)

    # Prevent propagation to root logger to avoid duplicate messages
    root_logger.propagate = False

    _logging_configured = True


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger for the given module name.

    This is the preferred way to obtain a logger in smartclass modules.
    The logger will be a child of the 'smartclass' logger, inheriting
    its configuration.

    :param name: Module name (typically __name__).
    :returns: Configured logger instance.

    Example:
        >>> from smartclass.logging import get_logger
        >>> logger = get_logger(__name__)
        >>> logger.info("Processing started")
    """
    # Ensure logging is configured with defaults if not already
    if not _logging_configured:
        configure_logging()

    return logging.getLogger(name)
