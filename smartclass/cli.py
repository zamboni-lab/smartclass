"""Command line interface for :mod:`smartclass`.

Why does this file exist, and why not put this in ``__main__``?
You might be tempted to import things from ``__main__`` later,
but that will cause problems--the code will get executed twice:

- When you run ``python3 -m smartclass`` python will execute``__main__.py`` as a script.
  That means there won't be any ``smartclass.__main__`` in ``sys.modules``.
- When you import __main__ it will get executed again (as a module) because
  there's no ``smartclass.__main__`` in ``sys.modules``.

.. seealso:: https://click.palletsprojects.com/en/8.1.x/setuptools/#setuptools-integration
"""

from __future__ import annotations

from pathlib import Path

import click

from smartclass.logging import configure_logging, get_logger


__all__ = [
    "main",
]

logger = get_logger(__name__)

# Default values as constants for maintainability
DEFAULT_ID_COLUMN = "class"
DEFAULT_SMARTS_COLUMN = "structure"
DEFAULT_WIKIDATA_ENDPOINT = "https://query.wikidata.org/sparql"


@click.group()
@click.version_option()
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity (-v for INFO, -vv for DEBUG).",
)
def main(verbose: int) -> None:
    """Smartclass - Classify chemical structures using SMARTS patterns.

    A tool for classifying chemical structures against SMARTS-based
    chemical class definitions from Wikidata or custom sources.

    \b
    Examples:
        # Classify a single molecule
        smartclass searchclasses -s "CCO" -c classes.tsv

        # Query Wikidata for chemical classes
        smartclass querywikidata -q query.rq -o output.tsv
    """
    # Configure logging based on verbosity
    log_levels = {0: "WARNING", 1: "INFO", 2: "DEBUG"}
    level = log_levels.get(min(verbose, 2), "DEBUG")
    configure_logging(level=level, force=True)


@main.command()
@click.option(
    "-i",
    "--input-file",
    type=click.Path(exists=True, path_type=Path),
    help="Input CSV/TSV file(s) to combine.",
    multiple=True,
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Output file path.",
    required=True,
)
def combinecsvfiles(input_file: tuple[Path, ...], output: Path) -> None:
    """Combine multiple CSV/TSV files into one.

    Merges rows from multiple input files with the same schema.
    """
    from smartclass.helpers import combine_csv_files

    combine_csv_files(input_file, output)
    click.echo(f"Combined {len(input_file)} files -> {output}")


@main.command()
@click.option(
    "-f",
    "--fp-len",
    type=int,
    default=2048,
    show_default=True,
    help="Fingerprint length for molecular fingerprints.",
)
@click.option(
    "-m",
    "--max-atoms",
    type=int,
    default=50,
    show_default=True,
    help="Maximum number of atoms in molecules to process.",
)
@click.option(
    "-r",
    "--report-interval",
    type=int,
    default=50000,
    show_default=True,
    help="Progress reporting interval (number of molecules).",
)
@click.option(
    "-t",
    "--tautomer-fingerprints/--no-tautomer-fingerprints",
    default=True,
    show_default=True,
    help="Include tautomer fingerprints in output.",
)
def getlatestchembl(
    fp_len: int,
    max_atoms: int,
    report_interval: int,
    tautomer_fingerprints: bool,
) -> None:
    """Download and process the latest ChEMBL database.

    Generates fingerprints for molecules in ChEMBL for use in
    classification and similarity searches.
    """
    from smartclass.resources.chembl import get_latest_chembl

    click.echo("Downloading and processing ChEMBL...")
    get_latest_chembl(fp_len, max_atoms, report_interval, tautomer_fingerprints)
    click.echo("ChEMBL processing complete.")


@main.command()
def loadpkgdata() -> None:
    """Load and display bundled package data.

    Shows the chemical classes, mappings, and MIA data included
    with the smartclass package.
    """
    from smartclass.io import load_pkg_data

    classes, mappings, mia = load_pkg_data()
    click.echo("=== Chemical Classes ===")
    click.echo(classes)
    click.echo("\n=== Mappings ===")
    click.echo(mappings)
    click.echo("\n=== MIA Data ===")
    click.echo(mia)


@main.command()
@click.option(
    "-q",
    "--query",
    type=click.Path(exists=True, path_type=Path),
    help="Path to SPARQL query file (.rq).",
    required=True,
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Output file path (CSV or TSV based on extension).",
    required=True,
)
@click.option(
    "-r",
    "--remove-prefix/--keep-prefix",
    "remove_prefix",
    default=True,
    show_default=True,
    help="Remove Wikidata entity prefix from results.",
)
@click.option(
    "-t",
    "--transform",
    type=click.Choice(
        [
            "check_smiles",
            "transform_inchi_to_inchikey",
            "transform_inchi_to_mass",
            "transform_inchi_to_smiles_canonical",
            "transform_inchi_to_smiles_isomeric",
            "transform_smiles_to_formula",
            "transform_smiles_to_inchi",
            "transform_smiles_to_mass",
            "transform_smiles_i_to_smiles_c",
            "transform_formula_to_formula",
        ],
        case_sensitive=False,
    ),
    help="Apply a transformation to query results.",
)
@click.option(
    "-u",
    "--url",
    type=str,
    default=DEFAULT_WIKIDATA_ENDPOINT,
    show_default=True,
    help="SPARQL endpoint URL.",
)
def querywikidata(
    query: Path,
    output: Path,
    remove_prefix: bool,
    transform: str | None,
    url: str,
) -> None:
    """Query Wikidata using SPARQL and export results.

    Execute a SPARQL query against Wikidata and optionally apply
    chemical transformations to the results.

    Examples:

    \b
        # Get chemical classes with SMARTS
        smartclass querywikidata -q classes_smarts.rq -o classes.tsv

    \b
        # Generate InChIKeys from InChI
        smartclass querywikidata -q inchi.rq -t transform_inchi_to_inchikey -o keys.csv
    """
    from smartclass.resources.wikidata import query_wikidata

    query_wikidata(str(query), str(output), remove_prefix, transform, url)
    click.echo(f"Results saved to: {output}")


@main.command()
@click.option(
    "-c",
    "--classes-file",
    type=click.Path(exists=True, path_type=Path),
    help="TSV file with chemical class definitions (class, structure columns).",
)
@click.option(
    "-d",
    "--classes-name-id",
    type=str,
    default=DEFAULT_ID_COLUMN,
    show_default=True,
    help="Column name for class identifiers.",
)
@click.option(
    "-e",
    "--classes-name-smarts",
    type=str,
    default=DEFAULT_SMARTS_COLUMN,
    show_default=True,
    help="Column name for SMARTS patterns.",
)
@click.option(
    "-f",
    "--include-hierarchy/--no-hierarchy",
    default=False,
    show_default=True,
    help="Use chemical hierarchy for faster BFS-based searching.",
)
@click.option(
    "-i",
    "--input-smiles",
    type=click.Path(exists=True, path_type=Path),
    help="Input file containing SMILES (CSV/TSV with 'smiles' column).",
)
@click.option(
    "-s",
    "--smiles",
    type=str,
    multiple=True,
    help="SMILES string(s) to classify. Can be specified multiple times.",
)
@click.option(
    "-z",
    "--closest-only/--all-matches",
    default=True,
    show_default=True,
    help="Return only closest matching class per structure.",
)
def searchclasses(
    classes_file: Path | None,
    classes_name_id: str,
    classes_name_smarts: str,
    closest_only: bool,
    include_hierarchy: bool,
    input_smiles: Path | None,
    smiles: tuple[str, ...],
) -> None:
    """Classify chemical structures against SMARTS-based classes.

    Match input SMILES strings against chemical class definitions
    using substructure searching. Results include the matched class
    and structural similarity metrics.

    Examples:

    \b
        # Classify a single molecule
        smartclass searchclasses -s "CCO" -c classes.tsv -v

    \b
        # Classify molecules from a file
        smartclass searchclasses -i molecules.tsv -c classes.tsv

    \b
        # Get all matches, not just closest
        smartclass searchclasses -s "CCO" -c classes.tsv --all-matches
    """
    from smartclass.chem.classification import search_classes

    if not smiles and not input_smiles:
        raise click.UsageError("Must provide either --smiles or --input-smiles")

    results = search_classes(
        classes_file=str(classes_file) if classes_file else None,
        classes_name_id=classes_name_id,
        classes_name_smarts=classes_name_smarts,
        closest_only=closest_only,
        include_hierarchy=include_hierarchy,
        input_smiles=str(input_smiles) if input_smiles else None,
        smiles=list(smiles) if smiles else None,
    )

    click.echo(f"Classification complete. Found {len(results)} matches.")
    logger.info(f"Results: {results}")


if __name__ == "__main__":
    main()
