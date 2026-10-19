"""Version information for :mod:`smartclass`.

Run with ``python -m smartclass.version``
"""

from __future__ import annotations

import os
from subprocess import CalledProcessError, check_output

from smartclass.logging import get_logger

__all__ = [
    "VERSION",
    "get_git_hash",
    "get_version",
]

logger = get_logger(__name__)

VERSION = "0.0.1-dev"


def get_git_hash() -> str:
    """Get the :mod:`smartclass` git hash.

    Returns
    -------
    str
        Hash.
    """
    with open(os.devnull, "w") as devnull:
        try:
            ret = check_output(
                ["git", "rev-parse", "HEAD"],
                cwd=os.path.dirname(__file__),
                stderr=devnull,
            )
        except CalledProcessError:
            return "UNHASHED"
        else:
            return ret.strip().decode("utf-8")[:8]


def get_version(with_git_hash: bool = False) -> str:
    """Get the :mod:`smartclass` version string, including a git hash.

    Parameters
    ----------
    with_git_hash : bool
        Whether to append the current git hash. Default is False.

    Returns
    -------
    str
        Version string.
    """
    return f"{VERSION}-{get_git_hash()}" if with_git_hash else VERSION


if __name__ == "__main__":
    logger.debug(get_version(with_git_hash=True))
