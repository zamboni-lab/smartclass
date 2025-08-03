"""Query Wikidata using SPARQL and export results."""

from __future__ import annotations

from collections.abc import Callable

from tqdm import tqdm

from smartclass.helpers import get_request, read_query
from smartclass.io import export_results
from smartclass.resources.wikidata.transforms import (
    check_smiles,
    transform_entities_to_stereoisomers,
    transform_formula_to_formula,
    transform_inchi_mass_to_inchi_mass,
    transform_inchi_to_inchikey,
    transform_inchi_to_mass,
    transform_inchi_to_smiles_canonical,
    transform_inchi_to_smiles_isomeric,
    transform_smiles_c_to_smiles_c_tauto,
    transform_smiles_i_to_smiles_c,
    transform_smiles_i_to_smiles_i,
    transform_smiles_i_to_smiles_i_tauto,
    transform_smiles_mass_to_smiles_mass,
    transform_smiles_to_formula,
    transform_smiles_to_inchi,
    transform_smiles_to_mass,
    transform_stereoisomers_to_entities,
)

__all__ = ["query_wikidata"]

TRANSFORM_FUNCTIONS: dict[str, Callable] = {
    "check_smiles": check_smiles,
    "transform_inchi_to_inchikey": transform_inchi_to_inchikey,
    "transform_inchi_to_mass": transform_inchi_to_mass,
    "transform_inchi_to_smiles_canonical": transform_inchi_to_smiles_canonical,
    "transform_inchi_to_smiles_isomeric": transform_inchi_to_smiles_isomeric,
    "transform_smiles_to_formula": transform_smiles_to_formula,
    "transform_smiles_to_inchi": transform_smiles_to_inchi,
    "transform_smiles_to_mass": transform_smiles_to_mass,
    "transform_smiles_i_to_smiles_c": transform_smiles_i_to_smiles_c,
    "transform_smiles_i_to_smiles_i": transform_smiles_i_to_smiles_i,
    "transform_formula_to_formula": transform_formula_to_formula,
    "transform_entities_to_stereoisomers": transform_entities_to_stereoisomers,
    "transform_stereoisomers_to_entities": transform_stereoisomers_to_entities,
    "transform_smiles_c_to_smiles_c_tauto": transform_smiles_c_to_smiles_c_tauto,
    "transform_smiles_i_to_smiles_i_tauto": transform_smiles_i_to_smiles_i_tauto,
    "transform_smiles_mass_to_smiles_mass": transform_smiles_mass_to_smiles_mass,
    "transform_inchi_mass_to_inchi_mass": transform_inchi_mass_to_inchi_mass,
}


def query_wikidata(
    query: str,
    output: str,
    remove_prefix: bool = True,
    transform: str | None = None,
    url: str = "https://query.wikidata.org/sparql",
) -> None:
    """
    Execute a SPARQL query on Wikidata and export the results.

    :param query: Path to the SPARQL query file.
    :type query: str

    :param output: Path to the output file where results will be exported.
    :type output: str

    :param remove_prefix: Flag to remove wikidata prefix.
    :type remove_prefix: bool

    :param transform: The transform to apply
    :type transform: Union[None,str]

    :param url: URL of the Wikidata SPARQL endpoint.
    :type url: str
    """
    query_str = read_query(query)

    results = get_request(url=url, query=query_str)
    results = [x for x in results if x is not None]

    if remove_prefix:
        prefix = "http://www.wikidata.org/entity/"
        prefix_len = len(prefix)
        for record in results:
            for key, value in record.items():
                if isinstance(value, str) and value.startswith(prefix):
                    record[key] = value[prefix_len:]

    if transform is not None:
        if transform in TRANSFORM_FUNCTIONS:
            transformed_results = []
            for item in tqdm(results, desc=f"Transforming with {transform}"):
                transformed = TRANSFORM_FUNCTIONS[transform](item)
                if transformed is not None:
                    transformed_results.append(transformed)
            results = transformed_results

    export_results(output=output, results=results)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)

    queries = [
        (
            "smartclass/data/queries/classes_cxsmiles.rq",
            "scratch/wikidata_classes_cxsmiles.tsv",
        ),
        (
            "smartclass/data/queries/classes_smarts.rq",
            "scratch/wikidata_classes_smarts.tsv",
        ),
        (
            "smartclass/data/queries/classes_smiles.rq",
            "scratch/wikidata_classes_smiles.tsv",
        ),
        (
            "smartclass/data/queries/classes_smiles_isomeric.rq",
            "scratch/wikidata_classes_smiles_isomeric.tsv",
        ),
        (
            "smartclass/data/queries/classes_taxonomy.rq",
            "scratch/wikidata_classes_taxonomy.tsv",
        ),
    ]

    for query, output in queries:
        query_wikidata(query, output)
    logging.info(output)
