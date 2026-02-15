"""Centralized configuration management for smartclass.

This module provides configuration settings with sensible defaults that can be
overridden via environment variables or programmatically.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import ClassVar

__all__ = [
    "Config",
    "get_config",
    "reset_config",
]


@dataclass
class Config:
    """Configuration settings for smartclass.

    Settings can be overridden via environment variables prefixed with SMARTCLASS_.
    For example, SMARTCLASS_OUTPUT_DIR=/path/to/output.

    Attributes:
        output_dir: Directory for output files (default: "output").
        cache_dir: Directory for cached data (default: ".smartclass_cache").
        wikidata_endpoint: Default Wikidata SPARQL endpoint URL.
        qlever_endpoint: Fallback QLever SPARQL endpoint URL.
        http_timeout: Default HTTP request timeout in seconds.
        http_max_retries: Maximum number of HTTP retry attempts.
        http_base_delay: Base delay for exponential backoff in seconds.
        chembl_fp_length: Default fingerprint length for ChEMBL processing.
        chembl_max_atoms: Maximum number of atoms for ChEMBL molecules.
        chembl_report_interval: Reporting interval for ChEMBL processing.
        mcs_timeout: Timeout for MCS calculations in seconds.
        mcs_threshold: Default threshold for MCS calculations.
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
    """

    # Output settings
    output_dir: Path = field(default_factory=lambda: Path("output"))
    cache_dir: Path = field(default_factory=lambda: Path(".smartclass_cache"))

    # SPARQL endpoints
    wikidata_endpoint: str = "https://query.wikidata.org/sparql"
    qlever_endpoint: str = "https://qlever.dev/api/wikidata"

    # HTTP settings
    http_timeout: int = 60
    http_max_retries: int = 3
    http_base_delay: float = 2.0

    # ChEMBL processing settings
    chembl_fp_length: int = 2048
    chembl_max_atoms: int = 50
    chembl_report_interval: int = 50000
    chembl_tautomer_fingerprints: bool = True

    # MCS settings
    mcs_timeout: int = 60
    mcs_threshold: float = 0.7

    # Logging
    log_level: str = "INFO"

    # User agent for HTTP requests
    user_agent: str = "SmartClassBot/1.0 (https://github.com/zamboni-lab/smartclass)"

    # Environment variable prefix
    _ENV_PREFIX: ClassVar[str] = "SMARTCLASS_"

    def __post_init__(self) -> None:
        """Load settings from environment variables after initialization."""
        self._load_from_env()
        self._ensure_directories()

    def _load_from_env(self) -> None:
        """Override settings from environment variables."""
        env_mappings = {
            "OUTPUT_DIR": ("output_dir", Path),
            "CACHE_DIR": ("cache_dir", Path),
            "WIKIDATA_ENDPOINT": ("wikidata_endpoint", str),
            "QLEVER_ENDPOINT": ("qlever_endpoint", str),
            "HTTP_TIMEOUT": ("http_timeout", int),
            "HTTP_MAX_RETRIES": ("http_max_retries", int),
            "HTTP_BASE_DELAY": ("http_base_delay", float),
            "CHEMBL_FP_LENGTH": ("chembl_fp_length", int),
            "CHEMBL_MAX_ATOMS": ("chembl_max_atoms", int),
            "CHEMBL_REPORT_INTERVAL": ("chembl_report_interval", int),
            "MCS_TIMEOUT": ("mcs_timeout", int),
            "MCS_THRESHOLD": ("mcs_threshold", float),
            "LOG_LEVEL": ("log_level", str),
            "USER_AGENT": ("user_agent", str),
        }

        for env_suffix, (attr_name, type_converter) in env_mappings.items():
            env_var = f"{self._ENV_PREFIX}{env_suffix}"
            env_value = os.environ.get(env_var)
            if env_value is not None:
                try:
                    setattr(self, attr_name, type_converter(env_value))
                except (ValueError, TypeError):
                    pass  # Keep default if conversion fails

    def _ensure_directories(self) -> None:
        """Create output and cache directories if they don't exist."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get_output_path(self, filename: str) -> Path:
        """Get full path for an output file.

        :param filename: Name of the output file.
        :returns: Full path to the output file.
        """
        return self.output_dir / filename

    def get_cache_path(self, filename: str) -> Path:
        """Get full path for a cache file.

        :param filename: Name of the cache file.
        :returns: Full path to the cache file.
        """
        return self.cache_dir / filename


# Global configuration instance (singleton pattern)
_config: Config | None = None


def get_config() -> Config:
    """Get the global configuration instance.

    Creates a new Config instance on first call, then returns the same
    instance on subsequent calls.

    :returns: The global Config instance.
    """
    global _config
    if _config is None:
        _config = Config()
    return _config


def reset_config() -> None:
    """Reset the global configuration to None.

    Useful for testing or when you need to reload configuration.
    """
    global _config
    _config = None
