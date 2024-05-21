"""Sample list."""

from __future__ import annotations

from random import sample

__all__ = ["sample_list"]


def sample_list(list: list, samples_max: int = 1000) -> list:
    """
    Sample from a list of things.

    :param list: List of things.
    :type list: list

    :param samples_max: Maximum number of samples to return.
    :type samples_max: int

    :returns: Sampled things.
    :rtype: list
    """
    return sample(list, min(len(list), samples_max))
