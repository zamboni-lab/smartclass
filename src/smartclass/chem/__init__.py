"""Smartclass classifies structures using SMARTS.chem."""

from __future__ import annotations

# from smartclass.chem.classification import check_classification
from smartclass.chem.classification import (
    bfs_search_classes_generator,  # noqa:F401
    get_class_structures,  # noqa:F401
    search_class,  # noqa:F401
    search_classes,  # noqa:F401
)
from smartclass.chem.conversion.convert_inchi_to_inchikey import (  # noqa:F401
    convert_inchi_to_inchikey,
)
from smartclass.chem.conversion.convert_inchi_to_mass import convert_inchi_to_mass  # noqa:F401
from smartclass.chem.conversion.convert_inchi_to_mol import convert_inchi_to_mol  # noqa:F401
from smartclass.chem.conversion.convert_inchi_to_smiles import convert_inchi_to_smiles  # noqa:F401
from smartclass.chem.conversion.convert_mol_to_cxsmiles import convert_mol_to_cxsmiles  # noqa:F401
from smartclass.chem.conversion.convert_mol_to_inchi import convert_mol_to_inchi  # noqa:F401
from smartclass.chem.conversion.convert_mol_to_inchikey import convert_mol_to_inchikey  # noqa:F401
from smartclass.chem.conversion.convert_mol_to_smarts import convert_mol_to_smarts  # noqa:F401
from smartclass.chem.conversion.convert_mol_to_smiles import convert_mol_to_smiles  # noqa:F401
from smartclass.chem.conversion.convert_molblock_to_mol import convert_molblock_to_mol  # noqa:F401

# from smartclass.chem.conversion.convert_name_and_molblock_to_cxsmiles import convert_mol_to_cxsmiles
from smartclass.chem.conversion.convert_smarts_to_mol import convert_smarts_to_mol  # noqa:F401
from smartclass.chem.conversion.convert_smiles_to_canonical_smiles import (  # noqa:F401
    convert_smiles_to_canonical_smiles,
)
from smartclass.chem.conversion.convert_smiles_to_formula import (  # noqa:F401
    convert_smiles_to_formula,
)
from smartclass.chem.conversion.convert_smiles_to_inchi import convert_smiles_to_inchi  # noqa:F401
from smartclass.chem.conversion.convert_smiles_to_mass import convert_smiles_to_mass  # noqa:F401
from smartclass.chem.conversion.convert_smiles_to_mol import convert_smiles_to_mol  # noqa:F401
from smartclass.chem.helpers.check_missing_stereochemistry import (  # noqa:F401
    check_missing_stereochemistry,
)
from smartclass.chem.helpers.enumerate_structures import enumerate_structures  # noqa:F401
from smartclass.chem.helpers.fix_inchi_tautomerization import fix_inchi_tautomerization  # noqa:F401
from smartclass.chem.helpers.get_num_atoms_bonds import get_num_atoms_bonds  # noqa:F401
from smartclass.chem.helpers.get_num_matched_atoms_bonds import (  # noqa:F401
    get_num_matched_atoms_bonds,
)
from smartclass.chem.helpers.standardize import standardize  # noqa:F401
from smartclass.chem.similarity.calculate_mcs import calculate_mcs  # noqa:F401

# from smartclass.chem.complexity import *


# from smartclass.chem.similarity import calculate_myopic_mces
# from smartclass.chem.similarity import calculate_structural_similarity
# from smartclass.chem.similarity import measure_mhfp
