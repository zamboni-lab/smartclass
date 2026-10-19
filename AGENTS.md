# AGENTS.md — smartclass development guide for AI coding agents

smartclass classifies chemical structures using SMARTS patterns against Wikidata-sourced chemical chemical classes.

## Stack

- Python ≥3.12 · uv · Click · rdkit · polars · pytest · ruff · ty · great-docs

## Quick start

```sh
uv sync                        # install all dependencies
uv run smartclass --help       # list CLI subcommands
```

## Verify changes (run both)

```sh
uv run prek run --all-files    # lint + format + type-check + security
uv run pytest tests/ -x        # 56 tests, ~3s, 47% coverage
```

Single-file test: `uv run pytest tests/test_conversions.py -v`

## Project layout

```
smartclass/                  # source package
├── chem/                    # chemical operations
│   ├── classification/      # core: search_classes.py, bfs_search_classes_generator.py
│   ├── conversion/          # SMILES↔InChI↔MOL↔SMARTS↔formula/mass (16 modules)
│   ├── helpers/             # stereochemistry, standardization, structure enumeration
│   ├── similarity/          # MCS, MCES, structural similarity
│   └── complexity/          # (stub)
├── io/                      # file I/O: TSV/CSV/JSON, package data loading, results export
├── helpers/                 # utilities: CSV, formula, HTTP requests, logging
├── resources/               # external data: chebi, chembl, chemont, wikidata, lotus
├── data/                    # bundled TSV lookup tables and SPARQL queries
├── cli.py                   # Click CLI entry point (6 subcommands)
├── config.py                # configuration (env vars: SMARTCLASS_*)
├── version.py               # VERSION constant
└── __init__.py              # package exports
tests/                       # pytest tests (8 files, 56 tests)
notebooks/                   # standalone scripts + Quarto notebook
scratch/                     # generated data (gitignored)
great-docs/                  # generated docs site (gitignored)
```

## CLI commands

| Command | Purpose |
|---------|---------|
| `querywikidata` | Query Wikidata SPARQL endpoint, optionally transform results |
| `searchclasses` | Classify SMILES against SMARTS-based classes |
| `searchclasses-all-sources` | Classify against SMARTS + SMILES + CXSMILES sources |
| `combinecsvfiles` | Combine multiple CSV/TSV files |
| `getlatestchembl` | Download and process latest ChEMBL database |
| `loadpkgdata` | Load and display bundled package data |

CLI convention: `-v`/`--verbose` is a **group-level** option — place it BEFORE the subcommand:

```sh
smartclass -v searchclasses -s "CCO" -c classes.tsv
```

## Conventions

- Docstrings: numpy-style (enforced by ruff, documented in great-docs.yml)
- Logging: use `smartclass.logging.get_logger(__name__)` — never `print()` except in CLI `click.echo()`
- Testing: pytest (tests use `unittest.TestCase` classes, run via pytest runner)
- Type checking: `ty` (not mypy, not pyright — pyright config in pyproject.toml is legacy)
- Version: managed in `pyproject.toml` (`version`) and `smartclass/version.py` (`VERSION`)
- Dependencies: `uv` (lockfile: `uv.lock`)

## Environment variables

All configurable via `SMARTCLASS_` prefix (see `smartclass/config.py`):

| Variable | Default | Description |
|----------|---------|-------------|
| `SMARTCLASS_OUTPUT_DIR` | `output` | Output directory |
| `SMARTCLASS_CACHE_DIR` | `.smartclass_cache` | Cache directory |
| `SMARTCLASS_WIKIDATA_ENDPOINT` | Wikidata SPARQL | SPARQL endpoint URL |
| `SMARTCLASS_HTTP_TIMEOUT` | `60` | HTTP timeout (seconds) |
| `SMARTCLASS_HTTP_MAX_RETRIES` | `3` | Max retries |
| `SMARTCLASS_LOG_LEVEL` | `INFO` | Log level |

## Important paths

- Package data: `smartclass/data/` (TSV tables, SPARQL queries)
- Test data: `smartclass/data/bitter_smiles.tsv`
- Generated scratch data: `scratch/` (gitignored — must be generated before use)

## Gotchas

- The bundled `smartclass/data/classes_smarts.tsv` uses columns `class_id` and `class_smarts`, not the CLI defaults `class` and `structure`. Use `--classes-name-id class_id --classes-name-smarts class_smarts`.
- `smartclass/chem/similarity/measure_mhfp.py` is entirely stub code — do not reference it in new code.
- `uv run ty check` is not available directly — use `uv run prek run ty-check --all-files`.
