"""Tests for helper functions."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path


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


class TestSplitCsv(unittest.TestCase):
    """Test split_csv function."""

    def test_splits_into_chunks_with_header(self):
        """Each output chunk keeps the header and respects lines_per_file."""
        from smartclass.helpers.split_csv import split_csv

        with tempfile.TemporaryDirectory() as tmp:
            input_file = Path(tmp) / "input.csv"
            input_file.write_text("a,b\n1,2\n3,4\n5,6\n", encoding="utf-8")

            output_dir = Path(tmp) / "out"
            split_csv(str(input_file), str(output_dir), lines_per_file=2)

            chunks = sorted(output_dir.glob("split_*.csv"))
            self.assertEqual(len(chunks), 2)
            for chunk in chunks:
                first_line = chunk.read_text(encoding="utf-8").splitlines()[0]
                self.assertEqual(first_line, "a,b")

    def test_rejects_non_positive_lines_per_file(self):
        """lines_per_file must be positive."""
        from smartclass.helpers.split_csv import split_csv

        with tempfile.TemporaryDirectory() as tmp:
            input_file = Path(tmp) / "input.csv"
            input_file.write_text("a\n1\n", encoding="utf-8")

            with self.assertRaises(ValueError):
                split_csv(str(input_file), str(Path(tmp) / "out"), lines_per_file=0)


class TestCombineCsvFiles(unittest.TestCase):
    """Test combine_csv_files function."""

    def test_combines_files(self):
        """Rows from all input files are concatenated."""
        import polars

        from smartclass.helpers.combine_csv_files import combine_csv_files

        with tempfile.TemporaryDirectory() as tmp:
            file_a = Path(tmp) / "a.tsv"
            file_b = Path(tmp) / "b.tsv"
            file_a.write_text("x\n1\n", encoding="utf-8")
            file_b.write_text("x\n2\n", encoding="utf-8")
            output = Path(tmp) / "combined.tsv"

            combine_csv_files([str(file_a), str(file_b)], str(output))

            df = polars.read_csv(output, separator="\t")
            self.assertEqual(df["x"].to_list(), [1, 2])

    def test_rejects_empty_input(self):
        """An empty input list is rejected with a clear error."""
        from smartclass.helpers.combine_csv_files import combine_csv_files

        with self.assertRaises(ValueError):
            combine_csv_files([], "out.tsv")


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
        assert mol is not None
        count = get_num_atoms_bonds(mol)

        # Ethanol: 3 heavy atoms + 2 bonds = 5 (without hydrogens)
        self.assertIsInstance(count, int)
        self.assertGreater(count, 0)


if __name__ == "__main__":
    unittest.main()
