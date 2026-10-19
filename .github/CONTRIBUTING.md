# Contributing

Contributions to this repository are welcomed and encouraged.

## Code Contribution

This project uses the [GitHub Flow](https://guides.github.com/introduction/flow)
model for code contributions. Follow these steps:

1. [Create a fork](https://help.github.com/articles/fork-a-repo) of the upstream
   repository at
   [`zambonilab/smartclass`](https://github.com/zamboni-lab/smartclass) on your
   GitHub account (or in one of your organizations)
2. [Clone your
   fork](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
   with `git clone https://github.com/<your namespace here>/smartclass.git`
3. Make and commit changes to your fork with `git commit`
4. Push changes to your fork with `git push`
5. Repeat steps 3 and 4 as needed
6. Submit a pull request back to the upstream repository

### Merge Model

This repository uses [squash
merges](https://docs.github.com/en/github/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/about-pull-request-merges#squash-and-merge-your-pull-request-commits)
to group all related commits in a given pull request into a single commit upon
acceptance and merge into the main branch. This has several benefits:

1. Keeps the commit history on the main branch focused on high-level narrative
2. Enables people to make lots of small commits without worrying about muddying
   up the commit history
3. Commits correspond 1-to-1 with pull requests

### Code Style

This project uses [static
typing](https://docs.python.org/3/library/typing.html). It uses
[`ty`](https://github.com/astral-sh/ty) as a type checker (configured in
`pyproject.toml` under `[tool.ty.rules]`). You can check your code with:

```sh
uv run prek run ty-check --all-files
```

This project uses [`ruff`](https://github.com/astral-sh/ruff) for linting and
formatting (configured in `pyproject.toml` under `[tool.ruff]`). Run the
following commands to check and fix code style issues:

```sh
# Check for lint issues
uv run ruff check .

# Auto-fix lint issues
uv run ruff check --fix .

# Format code
uv run ruff format .
```

### Logging

Python's builtin `print()` should not be used (except when writing to files). In
a command line setting or `main()` function for a module, you can use
`click.echo()`. Otherwise, you should use the centralized logging module:

```python
from smartclass.logging import get_logger

logger = get_logger(__name__)
logger.info("Processing started")
```

### Pre-commit Hooks

All checks (ruff, ty, typos, tombi, deptry, zizmor, and more) are run via
[`prek`](https://prek.j178.dev) before each commit. To install the hooks:

```sh
uv run prek install
```

To run all checks manually:

```sh
uv run prek run --all-files
```

### Documentation

All public functions (i.e., not starting with an underscore `_`) must be
documented using numpy-style docstrings. Example:

```python
def classify_smiles(smiles: str, closest_only: bool = True) -> list[dict]:
    """Classify a SMILES string against chemical classes.

    Parameters
    ----------
    smiles : str
        The SMILES string to classify.
    closest_only : bool
        If True, return only the closest match. Default is True.

    Returns
    -------
    list[dict]
        List of classification results.

    Raises
    ------
    SMILESError
        If the SMILES cannot be parsed.
    """
```

Documentation is built with
[`great-docs`](https://github.com/posit-dev/great-docs). Build locally with:

```sh
uv sync --group docs
uv run great-docs build
```

### Testing

Tests are written using the `unittest` framework (executed via pytest) in the
`tests/` directory. Run tests with:

```sh
# Run all tests
uv run pytest tests/

# Run with coverage
uv run pytest tests/ --cov=smartclass --cov-report=term-missing

# Run specific test file
uv run pytest tests/test_conversions.py -v
```

All tests must pass before a contribution can be accepted.

### Syncing your fork

If other code is updated before your contribution gets merged, you might need to
resolve conflicts against the main branch. After cloning, you should add the
upstream repository with:

```sh
git remote add upstream https://github.com/zamboni-lab/smartclass.git
git fetch upstream
git merge upstream/main
```

### Python Version Compatibility

This project supports Python 3.12 and later.

See https://endoflife.date/python for a timeline of Python release and
end-of-life dates.

## Acknowledgements

These code contribution guidelines are derived from the
[cthoyt/cookiecutter-snekpack](https://github.com/cthoyt/cookiecutter-snekpack)
Python package template. They're free to reuse and modify as long as they're
properly acknowledged.
