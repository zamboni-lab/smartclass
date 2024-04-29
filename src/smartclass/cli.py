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

import logging

import click

__all__ = [
    "main",
]

logger = logging.getLogger(__name__)


@click.group()
@click.version_option()
def main():
    """CLI for smartclass."""


@main.command()
@click.option(
    "-i",
    "--input-file",
    type=click.Path(exists=True),
    help="Input file.",
    multiple=True,
    required=True,
)
@click.option("-o", "--output", type=click.Path(), help="Output file.", required=True)
def combinecsvfiles(input_file, output):
    """CLI command that calls combine_csv_files function."""
    from smartclass.helpers import combine_csv_files

    combine_csv_files(input_file, output)
    click.echo(output)


# TODO add parameters
@main.command()
def getlatestchembl():
    """CLI command that calls get_latest_chembl function."""
    from smartclass.resources.chembl import get_latest_chembl

    get_latest_chembl()


@main.command()
def loadpkgdata():
    """CLI command that calls load_pkg_data function."""
    from smartclass.io import load_pkg_data

    classes, mappings, mia = load_pkg_data()
    click.echo(classes)
    click.echo(mappings)
    click.echo(mia)


@main.command()
@click.option("-q", "--query", type=click.Path(exists=True), help="Query file.", required=True)
@click.option("-o", "--output", type=click.Path(), help="Output file.", required=True)
@click.option(
    "-r", "--remove_prefix", type=bool, help="Remove prefix.", required=False, default=True
)
@click.option("-t", "--transform", type=str, help="Transform to apply.", required=False)
def querywikidata(query, output, remove_prefix, transform):
    """CLI command that calls query_wikidata function."""
    from smartclass.resources.wikidata import query_wikidata

    query_wikidata(query, output, remove_prefix, transform)
    click.echo(output)


@main.command()
@click.option("-c", "--classes-file", type=click.Path(exists=True), help="Input file (classes).")
@click.option(
    "-d", "--classes-name-id", type=str, help="Name of the ID column in the classes file."
)
@click.option(
    "-e", "--classes-name-smarts", type=str, help="Name of the SMARTS column in the classes file."
)
@click.option("-f", "--include-hierarchy", type=bool, help="Use a chemical hierarchy to go faster.")
@click.option("-i", "--input-smiles", type=click.Path(exists=True), help="Input file (SMILES).")
@click.option("-s", "--smiles", type=str, multiple=True, help="(List of) SMILES string(s)")
@click.option("-z", "--closest-only", type=bool, default=True, help="Return closest only.")
def searchclasses(
    classes_file,
    classes_name_id,
    classes_name_smarts,
    closest_only,
    include_hierarchy,
    input_smiles,
    smiles,
):
    """CLI command that calls search_classes function."""
    from smartclass.chem.classification import search_classes

    if not classes_name_id:
        classes_name_id = "class"
    if not classes_name_smarts:
        classes_name_smarts = "structure"

    search_classes(
        classes_file,
        classes_name_id,
        classes_name_smarts,
        closest_only,
        include_hierarchy,
        input_smiles,
        smiles,
    )
    click.echo("Done")


if __name__ == "__main__":
    main()
