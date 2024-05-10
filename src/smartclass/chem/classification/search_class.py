"""Search for a single class and return the results as a list of dictionaries."""

from __future__ import annotations

from rdkit.Chem import FilterCatalog, SubstructMatchParameters

from smartclass.chem.conversion.convert_mol_to_inchikey import convert_mol_to_inchikey
from smartclass.chem.conversion.convert_smarts_to_mol import convert_smarts_to_mol
from smartclass.chem.helpers.enumerate_structures import enumerate_structures
from smartclass.chem.helpers.get_num_atoms_bonds import get_num_atoms_bonds
from smartclass.chem.helpers.get_num_matched_atoms_bonds import get_num_matched_atoms_bonds


def search_class(
    class_dict: dict[str, list[str]],
    structures: list,
    params: SubstructMatchParameters,
) -> list[dict[str, str]]:
    """
    Search for a single class and return the results as a list of dictionaries.

    :param class_dict: Dictionary of the chemical classes.
    :type class_dict: dict[str, list[str]]

    :param structures: List of structures.
    :type structures: list

    :param params: Parameters for matching.
    :type params: SubstructMatchParameters


    :returns: A list of dictionaries containing matches.
    :rtype: list[dict[str, str]]
    """
    results: list = []

    # COMMENT: This is the actual version. TODO params could be removed
    # TODO See https://greglandrum.github.io/rdkit-blog/posts/2021-05-13-intro-to-the-molecule-enumerator.html#using-molbundles-for-substructure-search  # noqa: E501
    for class_id, class_structures in class_dict.items():
        for class_structure in class_structures:
            catalog = FilterCatalog.FilterCatalog()
            # patterns = [convert_smarts_to_mol(class_structure)]
            # TODO looks important
            # Chem.Kekulize(mol)
            # Chem.SetGenericQueriesFromProperties(mol)
            # mol = Chem.AdjustQueryProperties(mol)
            patterns = enumerate_structures(mol=convert_smarts_to_mol(class_structure))
            # TODO index as we don't know if the catalog can have multiple matchers
            for index, pattern in enumerate(patterns):
                catalog.AddEntry(
                    FilterCatalog.FilterCatalogEntry(
                        name=f"({class_id}_{index})", matcher=FilterCatalog.SmartsMatcher(pattern)
                    )
                )
            for structure in structures:
                if catalog.HasMatch(structure):
                    result = {
                        "class_id": class_id,
                        "class_structure": class_structure,
                        "structure_inchikey": convert_mol_to_inchikey(structure),
                        "structure_ab": get_num_atoms_bonds(structure),
                        "matched_ab": get_num_matched_atoms_bonds(mol_1=structure, mol_2=pattern),
                    }
                    results.append(result)

    return results
