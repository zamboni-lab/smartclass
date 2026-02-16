"""Wikidata transform functions for chemical data."""
# ruff: noqa: E501
# Long module paths are unavoidable due to package structure

from __future__ import annotations

from smartclass.resources.wikidata.transforms.check_smiles import check_smiles
from smartclass.resources.wikidata.transforms.transform_entities_to_stereoisomers import (
    transform_entities_to_stereoisomers,
)
from smartclass.resources.wikidata.transforms.transform_formula_to_formula import (
    transform_formula_to_formula,
)
from smartclass.resources.wikidata.transforms.transform_inchi_mass_to_inchi_mass import (
    transform_inchi_mass_to_inchi_mass,
)
from smartclass.resources.wikidata.transforms.transform_inchi_to_inchikey import (
    transform_inchi_to_inchikey,
)
from smartclass.resources.wikidata.transforms.transform_inchi_to_mass import (
    transform_inchi_to_mass,
)
from smartclass.resources.wikidata.transforms.transform_inchi_to_smiles_canonical import (
    transform_inchi_to_smiles_canonical,
)
from smartclass.resources.wikidata.transforms.transform_inchi_to_smiles_isomeric import (
    transform_inchi_to_smiles_isomeric,
)
from smartclass.resources.wikidata.transforms.transform_smiles_c_to_smiles_c_tauto import (
    transform_smiles_c_to_smiles_c_tauto,
)
from smartclass.resources.wikidata.transforms.transform_smiles_i_to_smiles_c import (
    transform_smiles_i_to_smiles_c,
)
from smartclass.resources.wikidata.transforms.transform_smiles_i_to_smiles_i import (
    transform_smiles_i_to_smiles_i,
)
from smartclass.resources.wikidata.transforms.transform_smiles_i_to_smiles_i_tauto import (
    transform_smiles_i_to_smiles_i_tauto,
)
from smartclass.resources.wikidata.transforms.transform_smiles_mass_to_smiles_mass import (
    transform_smiles_mass_to_smiles_mass,
)
from smartclass.resources.wikidata.transforms.transform_smiles_to_formula import (
    transform_smiles_to_formula,
)
from smartclass.resources.wikidata.transforms.transform_smiles_to_inchi import (
    transform_smiles_to_inchi,
)
from smartclass.resources.wikidata.transforms.transform_smiles_to_mass import (
    transform_smiles_to_mass,
)
from smartclass.resources.wikidata.transforms.transform_stereoisomers_to_entities import (
    transform_stereoisomers_to_entities,
)


__all__ = [
    "check_smiles",
    "transform_entities_to_stereoisomers",
    "transform_formula_to_formula",
    "transform_inchi_mass_to_inchi_mass",
    "transform_inchi_to_inchikey",
    "transform_inchi_to_mass",
    "transform_inchi_to_smiles_canonical",
    "transform_inchi_to_smiles_isomeric",
    "transform_smiles_c_to_smiles_c_tauto",
    "transform_smiles_i_to_smiles_c",
    "transform_smiles_i_to_smiles_i",
    "transform_smiles_i_to_smiles_i_tauto",
    "transform_smiles_mass_to_smiles_mass",
    "transform_smiles_to_formula",
    "transform_smiles_to_inchi",
    "transform_smiles_to_mass",
    "transform_stereoisomers_to_entities",
]
