# Contributing

Contributions to this repository are welcomed and encouraged.

## Code Contribution

This project uses the [GitHub Flow](https://guides.github.com/introduction/flow)
model for code contributions. Follow these steps:

1. [Create a fork](https://help.github.com/articles/fork-a-repo) of the upstream
   repository at [`zambonilab/smartclass`](https://gitlab.ethz.ch/zambonilab/smartclass)
   on your GitHub account (or in one of your organizations)
2. [Clone your fork](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
   with `git clone https://gitlab.ethz.ch/<your namespace here>/smartclass.git`
3. Make and commit changes to your fork with `git commit`
4. Push changes to your fork with `git push`
5. Repeat steps 3 and 4 as needed
6. Submit a pull request back to the upstream repository

### Merge Model

This repository
uses [squash merges](https://docs.github.com/en/github/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/about-pull-request-merges#squash-and-merge-your-pull-request-commits)
to group all related commits in a given pull request into a single commit upon
acceptance and merge into the main branch. This has several benefits:

1. Keeps the commit history on the main branch focused on high-level narrative
2. Enables people to make lots of small commits without worrying about muddying
   up the commit history
3. Commits correspond 1-to-1 with pull requests

### Code Style

This project encourages the use of optional static typing. It
uses [`mypy`](http://mypy-lang.org/) as a type checker
and [`sphinx_autodoc_typehints`](https://github.com/agronholm/sphinx-autodoc-typehints)
to automatically generate documentation based on type hints. You can check if
your code passes `mypy` with `tox -e mypy`.

This project uses [`black`](https://github.com/psf/black) to automatically
enforce a consistent code style. You can apply `black` and other pre-configured
linters with `tox -e lint`.

This project uses [`flake8`](https://flake8.pycqa.org) and several plugins for
additional checks of documentation style, security issues, good variable
nomenclature, and more (
see [`tox.ini`](tox.ini) for a list of flake8 plugins). You can check if your
code passes `flake8` with `tox -e flake8`.

Each of these checks are run on each commit using GitHub Actions as a continuous
integration service. Passing all of them is required for accepting a
contribution. If you're unsure how to address the feedback from one of these
tools, please say so either in the description of your pull request or in a
comment, and we will help you.

### Logging

Python's builtin `print()` should not be used (except when writing to files),
it's checked by the
[`flake8-print`](https://github.com/jbkahn/flake8-print) plugin to `flake8`. If
you're in a command line setting or `main()` function for a module, you can use
`click.echo()`. Otherwise, you should use the centralized logging module:

```python
from smartclass.logging import get_logger

logger = get_logger(__name__)
logger.info("Processing started")
```

### Code Quality

This project uses [`ruff`](https://github.com/astral-sh/ruff) for linting and formatting.
Run the following commands to check and fix code style issues:

```shell
# Check for issues
uv run ruff check .

# Auto-fix issues
uv run ruff check --fix .

# Format code
uv run ruff format .
```

### Documentation

All public functions (i.e., not starting with an underscore `_`) must be
documented using Google-style docstrings. Example:

```python
def classify_smiles(smiles: str, closest_only: bool = True) -> list[dict]:
    """Classify a SMILES string against chemical classes.

    Args:
        smiles: The SMILES string to classify.
        closest_only: If True, return only the closest match.

    Returns:
        List of classification results.

    Raises:
        SMILESError: If the SMILES cannot be parsed.
    """
```

This project uses [`sphinx`](https://www.sphinx-doc.org) to automatically build
documentation. You can build it locally with:

```shell
uv run sphinx-build docs docs/_build
```

### Testing

Functions in this repository should be unit tested. Tests are written using
the `unittest` framework in the `tests/` directory. Run tests with:

```shell
# Run all tests
uv run pytest

# Run with coverage
uv run pytest --cov=smartclass --cov-report=html

# Run specific test file
uv run pytest tests/test_api.py -v
```

All tests must pass before a contribution can be accepted.

### Syncing your fork

If other code is updated before your contribution gets merged, you might need to
resolve conflicts against the main branch. After cloning, you should add the
upstream repository with:

```shell
git remote add upstream https://gitlab.ethz.ch/zambonilab/smartclass.git
git fetch upstream
git merge upstream/main
```

### Python Version Compatibility

This project supports Python 3.11 and later. We aim to support all versions
of Python that have not passed their end-of-life dates.

See https://endoflife.date/python for a timeline of Python release and
end-of-life dates.

## Acknowledgements

These code contribution guidelines are derived from
the [cthoyt/cookiecutter-snekpack](https://github.com/cthoyt/cookiecutter-snekpack)
Python package template. They're free to reuse and modify as long as they're properly acknowledged.
