"""Trivial version test."""

from __future__ import annotations

import unittest

from smartclass.version import get_version


class TestVersion(unittest.TestCase):
    """Trivially test a version."""

    def test_version_type(self):
        """Test the version is a string.

        This is only meant to be an example test.
        """
        version = get_version()
        self.assertIsInstance(version, str)
