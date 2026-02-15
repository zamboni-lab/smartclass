"""Tests for chemical conversion functions."""

from __future__ import annotations

import unittest

from rdkit.Chem import Mol


class TestConvertSmilesToMol(unittest.TestCase):
    """Test SMILES to Mol conversion."""

    def test_valid_smiles(self):
        """Test conversion of valid SMILES."""
        from smartclass.chem.conversion.convert_smiles_to_mol import (
            convert_smiles_to_mol,
        )

        mol = convert_smiles_to_mol("CCO")

        self.assertIsInstance(mol, Mol)

    def test_invalid_smiles_returns_none(self):
        """Test invalid SMILES returns None."""
        from smartclass.chem.conversion.convert_smiles_to_mol import (
            convert_smiles_to_mol,
        )

        mol = convert_smiles_to_mol("invalid_smiles_xyz")

        self.assertIsNone(mol)

    def test_empty_string_returns_none(self):
        """Test empty string returns None."""
        from smartclass.chem.conversion.convert_smiles_to_mol import (
            convert_smiles_to_mol,
        )

        mol = convert_smiles_to_mol("")

        self.assertIsNone(mol)

    def test_none_input_returns_none(self):
        """Test None input returns None."""
        from smartclass.chem.conversion.convert_smiles_to_mol import (
            convert_smiles_to_mol,
        )

        mol = convert_smiles_to_mol(None)

        self.assertIsNone(mol)

    def test_complex_smiles(self):
        """Test conversion of complex organic molecule."""
        from smartclass.chem.conversion.convert_smiles_to_mol import (
            convert_smiles_to_mol,
        )

        # Caffeine SMILES
        caffeine = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        mol = convert_smiles_to_mol(caffeine)

        self.assertIsInstance(mol, Mol)


class TestConvertMolToSmiles(unittest.TestCase):
    """Test Mol to SMILES conversion."""

    def test_mol_to_smiles(self):
        """Test conversion from Mol to SMILES."""
        from smartclass.chem.conversion.convert_smiles_to_mol import (
            convert_smiles_to_mol,
        )
        from smartclass.chem.conversion.convert_mol_to_smiles import (
            convert_mol_to_smiles,
        )

        mol = convert_smiles_to_mol("CCO")
        smiles = convert_mol_to_smiles(mol)

        self.assertIsInstance(smiles, str)
        self.assertTrue(len(smiles) > 0)


class TestConvertMolToInchikey(unittest.TestCase):
    """Test Mol to InChIKey conversion."""

    def test_mol_to_inchikey(self):
        """Test conversion from Mol to InChIKey."""
        from smartclass.chem.conversion.convert_smiles_to_mol import (
            convert_smiles_to_mol,
        )
        from smartclass.chem.conversion.convert_mol_to_inchikey import (
            convert_mol_to_inchikey,
        )

        mol = convert_smiles_to_mol("CCO")
        inchikey = convert_mol_to_inchikey(mol)

        self.assertIsInstance(inchikey, str)
        # InChIKey format: 14 chars + hyphen + 10 chars + hyphen + 1 char
        self.assertEqual(len(inchikey), 27)
        self.assertIn("-", inchikey)


class TestConvertInchiToSmiles(unittest.TestCase):
    """Test InChI to SMILES conversion."""

    def test_valid_inchi(self):
        """Test conversion of valid InChI."""
        from smartclass.chem.conversion.convert_inchi_to_smiles import (
            convert_inchi_to_smiles,
        )

        # Ethanol InChI
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        smiles = convert_inchi_to_smiles(inchi)

        self.assertIsInstance(smiles, str)
        self.assertTrue(len(smiles) > 0)

    def test_invalid_inchi_returns_none(self):
        """Test invalid InChI returns None."""
        from smartclass.chem.conversion.convert_inchi_to_smiles import (
            convert_inchi_to_smiles,
        )

        smiles = convert_inchi_to_smiles("invalid_inchi")

        self.assertIsNone(smiles)


class TestRoundTripConversions(unittest.TestCase):
    """Test round-trip conversions preserve structure identity."""

    def test_smiles_to_mol_to_smiles(self):
        """Test SMILES -> Mol -> SMILES preserves structure."""
        from smartclass.chem.conversion.convert_smiles_to_mol import (
            convert_smiles_to_mol,
        )
        from smartclass.chem.conversion.convert_mol_to_smiles import (
            convert_mol_to_smiles,
        )
        from smartclass.chem.conversion.convert_mol_to_inchikey import (
            convert_mol_to_inchikey,
        )

        original_smiles = "c1ccccc1"  # Benzene
        mol = convert_smiles_to_mol(original_smiles)
        result_smiles = convert_mol_to_smiles(mol)

        # Convert both to InChIKey for comparison (canonical)
        mol1 = convert_smiles_to_mol(original_smiles)
        mol2 = convert_smiles_to_mol(result_smiles)

        self.assertEqual(
            convert_mol_to_inchikey(mol1),
            convert_mol_to_inchikey(mol2),
        )


if __name__ == "__main__":
    unittest.main()
