# Architecture

## Overview

smartclass is a Python package for classifying chemical structures using SMARTS
(SMiles ARbitrary Transformation and Search) patterns against chemical class
definitions sourced from Wikidata.

## Data Flow

```
Input SMILES/InChI
  → convert_smiles_to_mol()   (rdkit)
  → substructure search       (rdkit SubstructMatch)
  → BFS hierarchy traversal   (bfs_search_classes_generator.py)
  → MatchResult list          (class_id, matched atoms/bonds, etc.)
  → export_results()          (TSV/CSV)
```

## Package Structure

### `smartclass/chem/` --- Chemical operations

- **`classification/`** --- Core substructure search and BFS hierarchy
  traversal:
  - `search_classes.py` (306 lines) --- main entry point for `search_classes()`
  - `bfs_search_classes_generator.py` (258 lines) --- BFS traversal of class
    hierarchies
  - `search_class.py` --- per-class substructure matching
  - `check_classification.py` --- validation utilities (0% coverage)
  - `get_class_structures.py` --- class-to-SMARTS mapping builder

- **`conversion/`** --- 16 modules for format conversion:
  - SMILES ↔ InChI ↔ MOL ↔ SMARTS ↔ formula/mass
  - Each module is a single function (5--18 lines)

- **`helpers/`** --- Stereochemistry, standardization, enumeration:
  - `canonicalize_smarts.py`, `check_missing_stereochemistry.py`,
    `enumerate_structures.py`, `fix_inchi_tauterization.py`, `standardize.py`

- **`similarity/`** --- MCS, MCES, structural similarity:
  - `calculate_mcs.py` --- Maximum Common Substructure
  - `calculate_myopic_mces.py` (0% coverage, stub)
  - `calculate_structural_similarity.py` (0% coverage, stub)
  - `measure_mhfp.py` --- entirely commented-out stub code

- **`complexity/`** --- (stub, empty `__init__.py`)

### `smartclass/io/` --- File I/O

- `load_pkg_*` modules for bundled package data (TSV, JSON, mappings)
- `load_external_classes_file.py` --- load TSV classes from external files
- `export_results.py` --- write results to TSV/CSV
- `export_dict_to_json.py` --- write dict to JSON

### `smartclass/resources/` --- External data sources

- **wikidata/** --- SPARQL queries and result transforms
  - `query_wikidata.py` --- SPARQL query execution
  - `transforms/` --- 16 transform modules (check_smiles,
    transform_inchi_to_inchikey, etc.)
- **chebi/** --- Chemical Entities of Biological Interest extraction
- **chembl/** --- ChEMBL database download and fingerprint generation
- **chemont/** --- Chemical Ontology hierarchy
- **lotus/** --- Natural products from LOTUS

### `smartclass/helpers/` --- Utility functions

- `combine_csv_files.py`, `split_csv.py` --- CSV manipulation
- `convert_chemical_formula.py` --- superscript/subscript conversion
- `get_request.py` --- HTTP requests with retry/backoff
- `read_query.py` --- SPARQL query file loading
- `sample_list.py`, `download_file_if_not_exists.py`

### `smartclass/data/` --- Bundled data

- `classes_smarts.tsv` --- 15 chemical classes with SMARTS patterns
- `mia_smarts.tsv` --- Mono Indole Alkaloids substructures
- `chemontids__wdids.tsv` --- ChEMONT↔Wikidata ID mappings
- `existing_classes.json`, `existing_classified_datasets.json`
- `bitter_smiles.tsv` --- test data
- `queries/` --- 26+ SPARQL query templates (.rq files)

## CLI Architecture

The CLI (`smartclass/cli.py`, 526 lines) is built with Click and organized as a
group with 6 subcommands:

  | Command                     | Purpose                                      |
  | ---------                   | ---------                                    |
  | `querywikidata`             | Run SPARQL queries against Wikidata          |
  | `searchclasses`             | Classify SMILES against SMARTS-based classes |
  | `searchclasses-all-sources` | Classify against SMARTS + SMILES + CXSMILES  |
  | `combinecsvfiles`           | Merge CSV/TSV files                          |
  | `getlatestchembl`           | Download/process ChEMBL database             |
  | `loadpkgdata`               | Display bundled package data                 |

**Important:** `-v`/`--verbose` is a group-level option. Place it BEFORE the
subcommand: `smartclass -v searchclasses ...`, not
`smartclass searchclasses ... -v`.

## Testing

- Framework: pytest (test classes use `unittest.TestCase`, run via pytest)
- Test files: 8 in `tests/` (test_config, test_conversions, test_exceptions,
  test_get_request, test_helpers, test_logging, test_version)
- 56 tests, ~3s runtime, 47% coverage
- **Gaps:** CLI (0% coverage), classification (0--17%), resources (0--32%)

## Tooling

  | Tool | Purpose | Command |
  | ------ | --------- | --------- |
  | uv | Dependency management + build | `uv sync`, `uv run ...` |
  | ruff | Lint + format | `uv run ruff check .`, `uv run ruff format .` |
  | ty | Type checker | `uv run prek run ty-check` |
  | prek | Pre-commit hooks | `uv run prek run --all-files` |
  | pytest | Tests | `uv run pytest tests/` |
  | great-docs | Documentation | `uv run great-docs build` |
  | zizmor | GitHub workflow security scanner | (via prek) |
