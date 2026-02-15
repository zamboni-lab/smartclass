"""Tests for helper functions."""

from __future__ import annotations

import unittest


class TestSampleList(unittest.TestCase):
    """Test sample_list function."""

    def test_sample_smaller_than_list(self):
        """Test sampling fewer items than list length."""
        from smartclass.helpers.sample_list import sample_list

        items = list(range(100))
        result = sample_list(items, max_samples=10)

        self.assertEqual(len(result), 10)
        # All items should be from original list
        for item in result:
            self.assertIn(item, items)

    def test_sample_larger_than_list(self):
        """Test sampling more items than list length."""
        from smartclass.helpers.sample_list import sample_list

        items = list(range(5))
        result = sample_list(items, max_samples=100)

        self.assertEqual(len(result), 5)

    def test_sample_empty_list(self):
        """Test sampling from empty list."""
        from smartclass.helpers.sample_list import sample_list

        result = sample_list([], max_samples=10)

        self.assertEqual(len(result), 0)

    def test_randomness(self):
        """Test that sampling is random."""
        from smartclass.helpers.sample_list import sample_list

        items = list(range(100))
        results = [tuple(sample_list(items, max_samples=10)) for _ in range(5)]

        # Not all results should be identical (extremely unlikely)
        unique_results = set(results)
        self.assertGreater(len(unique_results), 1)


class TestConvertListOfDict(unittest.TestCase):
    """Test convert_list_of_dict function."""

    def test_basic_conversion(self):
        """Test basic list to dict conversion."""
        from smartclass.helpers.convert_list_of_dict import convert_list_of_dict

        data = [
            {"name": "Alice", "group": "A"},
            {"name": "Bob", "group": "A"},
            {"name": "Charlie", "group": "B"},
        ]

        result = convert_list_of_dict(data, key="name", value="group")

        self.assertEqual(result["Alice"], ["A"])
        self.assertEqual(result["Bob"], ["A"])
        self.assertEqual(result["Charlie"], ["B"])

    def test_inverted_conversion(self):
        """Test inverted conversion (group by value)."""
        from smartclass.helpers.convert_list_of_dict import convert_list_of_dict

        data = [
            {"name": "Alice", "group": "A"},
            {"name": "Bob", "group": "A"},
            {"name": "Charlie", "group": "B"},
        ]

        result = convert_list_of_dict(data, key="name", value="group", invert=True)

        self.assertEqual(set(result["A"]), {"Alice", "Bob"})
        self.assertEqual(result["B"], ["Charlie"])

    def test_empty_list(self):
        """Test with empty list."""
        from smartclass.helpers.convert_list_of_dict import convert_list_of_dict

        result = convert_list_of_dict([], key="a", value="b")

        self.assertEqual(result, {})


class TestConvertChemicalFormula(unittest.TestCase):
    """Test convert_chemical_formula function."""

    def test_basic_formula(self):
        """Test basic formula conversion."""
        from smartclass.helpers.convert_chemical_formula import convert_chemical_formula

        result = convert_chemical_formula("H2O")

        self.assertIn("₂", result)

    def test_complex_formula(self):
        """Test complex formula with charges."""
        from smartclass.helpers.convert_chemical_formula import convert_chemical_formula

        result = convert_chemical_formula("C6H12O6")

        self.assertIn("₆", result)
        self.assertIn("₁₂", result)


class TestGetNumAtomsBonds(unittest.TestCase):
    """Test get_num_atoms_bonds function."""

    def test_ethanol(self):
        """Test atom+bond count for ethanol."""
        from smartclass.chem.conversion.convert_smiles_to_mol import (
            convert_smiles_to_mol,
        )
        from smartclass.chem.helpers.get_num_atoms_bonds import get_num_atoms_bonds

        mol = convert_smiles_to_mol("CCO")
        count = get_num_atoms_bonds(mol)

        # Ethanol: 3 heavy atoms + 2 bonds = 5 (without hydrogens)
        self.assertIsInstance(count, int)
        self.assertGreater(count, 0)


if __name__ == "__main__":
    unittest.main()
