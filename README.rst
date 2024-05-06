Smartclass |zenodo| |build| |coverage| |documentation| |black|
=======================================================================
Short Description TODO

To get a (very) small overview

 .. code-block:: sh

    poetry run quarto render notebooks/smartclass.qmd

🚀 Installation
---------------
..
 Uncomment this section after your first ``tox -e finish``
 The most recent release can be installed from
 `PyPI <https://pypi.org/project/smartclass/>`_ with:

 .. code-block:: sh

    pip install smartclass

The most recent code and data can be installed directly using poetry:

.. code-block:: sh

    poetry install

💪 Getting Started
------------------
For now:

Get defined chemical classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    poetry run smartclass querywikidata -q src/smartclass/data/queries/classes_cxsmiles.rq -o scratch/wikidata_classes_cxsmiles.tsv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/classes_smarts.rq -o scratch/wikidata_classes_smarts.tsv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/classes_smiles_canonical.rq -o scratch/wikidata_classes_smiles_canonical.tsv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/classes_smiles_isomeric.rq -o scratch/wikidata_classes_smiles_isomeric.tsv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/classes_taxonomy.rq -o scratch/wikidata_classes_taxonomy.tsv

Combine them
~~~~~~~~~~~~
.. code-block:: sh

    poetry run smartclass combinecsvfiles -i scratch/wikidata_classes_smiles_canonical.tsv -i scratch/wikidata_classes_smiles_isomeric.tsv -i scratch/wikidata_classes_smarts.tsv -i scratch/wikidata_classes_cxsmiles.tsv -o scratch/wikidata_classes_full.tsv

Get some structures
~~~~~~~~~~~~~~~~~~~
.. code-block:: sh
    
    poetry run smartclass getlatestchembl
    # poetry run python3 smartclass/extract_smiles_lotus.py

Load package data
~~~~~~~~~~~~~~~~~~~
.. code-block:: sh
    
    poetry run smartclass loadpkgdata

Classify compounds
~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    poetry run smartclass searchclasses --help
    # poetry run smartclass searchclasses
    poetry run smartclass searchclasses -s "O=C1OC2CC3C(C=C(OC)C(=O)C3(C)C4C(=O)C(OC)=C(C)C(C1)C24C)C" -c scratch/wikidata_classes_full.tsv --verbose
    poetry run smartclass searchclasses -s "O=C1OC2CC3C(C=C(OC)C(=O)C3(C)C4C(=O)C(OC)=C(C)C(C1)C24C)C" -c scratch/wikidata_classes_full.tsv --closest-only False --verbose
    poetry run smartclass searchclasses -i src/smartclass/data/bitter_smiles.tsv -c scratch/wikidata_classes_full.tsv

Measure substructures distances
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    poetry run python3 smartclass/measure_mhfp.py -i scratch/wikidata_classes_smarts.tsv -o scratch/distances_wikidata.tsv
    poetry run python3 smartclass/measure_mhfp.py -i data/mia_smarts.tsv -o scratch/distances_mia.tsv

Wikidata maintenance
~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_canonical_no_formula.rq -t transform_smiles_to_formula -o scratch/formulas_canonical.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_no_smiles_isomeric_no_mass.rq -t transform_inchi_to_mass -o scratch/masses_inchi.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_isomeric_no_mass.rq -t transform_smiles_to_mass -o scratch/masses_isomeric.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_isomeric_no_inchi.rq -t transform_smiles_to_inchi -o scratch/inchis_isomeric.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_isomeric_smiles_inchi_no_stereo.rq -t transform_smiles_to_inchi -o scratch/incorrect_inchis.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_isomeric_no_formula.rq -t transform_smiles_to_formula -o scratch/formulas_isomeric.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_isomeric_no_canonical.rq -t transform_smiles_i_to_smiles_c -o scratch/smiles_isomeric.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_no_smiles_isomeric_inchi.rq -t transform_inchi_to_smiles_isomeric -o scratch/smiles_i_inchi.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_inchi_no_isomeric_smiles.rq -t transform_inchi_to_smiles_isomeric -o scratch/smiles_i_inchi_2.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_inchi_no_canonical_smiles.rq -t transform_inchi_to_smiles_canonical -o scratch/smiles_c_inchi.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_inchi_no_inchikey.rq -t transform_inchi_to_inchikey -o scratch/inchikeys.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_formula_malformed.rq -t transform_formula_to_formula -o scratch/formulas_malformed.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_isomeric_incorrect_mass.rq -t transform_smiles_mass_to_smiles_mass -o scratch/masses_incorrect_isomeric.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/maintenance_no_smiles_isomeric_incorrect_mass.rq -t transform_inchi_mass_to_inchi_mass -o scratch/masses_incorrect_inchi.csv

    poetry run smartclass querywikidata -q src/smartclass/data/queries/chemicals_smiles_canonical_no_ref.rq -t transform_smiles_c_to_smiles_c_tauto -o scratch/smiles_c_tauto.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/chemicals_smiles_isomeric_no_ref.rq -t transform_smiles_i_to_smiles_i -o scratch/smiles_i.csv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/chemicals_smiles_isomeric_no_ref.rq -t transform_smiles_i_to_smiles_i_tauto -o scratch/smiles_i_tauto.csv
    # not working for now
    # poetry run smartclass querywikidata -q src/smartclass/data/queries/stereoisomers_smiles_isomeric_inchi.rq -t transform_stereoisomers_to_entities -o scratch/stereo_i_to_entities.csv
    # poetry run smartclass querywikidata -q src/smartclass/data/queries/stereoisomers_smiles_canonical_no_isomeric.rq -t transform_stereoisomers_to_entities -o scratch/stereo_c_to_entities.csv
    # poetry run smartclass querywikidata -q src/smartclass/data/queries/chemical_entities_smiles_isomeric_inchi.rq -t transform_entities_to_stereoisomers -o scratch/entities_i_to_stereo.csv
    # poetry run smartclass querywikidata -q src/smartclass/data/queries/chemical_entities_smiles_canonical_no_isomeric.rq -t transform_entities_to_stereoisomers -o scratch/entities_c_to_stereo.csv

Improve current classes
~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    poetry run smartclass querywikidata -q src/smartclass/data/queries/chemicals_classes.rq -o scratch/wikidata_chemicals_classes.tsv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/chemical_entities_smiles_isomeric_inchi.rq -o scratch/wikidata_chemical_entities_smiles_inchi.tsv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/stereoisomers_smiles_isomeric_inchi.rq -o scratch/wikidata_stereoisomers_smiles_isomeric_inchi.tsv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/stereoisomers_smiles_canonical_no_isomeric_inchi.rq -o scratch/wikidata_stereoisomers_smiles_canonical_inchi.tsv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/chemicals_tautomer_of.rq -o scratch/wikidata_chemicals_tautomer_of.tsv
    poetry run python3 notebooks/improve_classes.py
    poetry run python3 src/smartclass/helpers/split_csv.py

Add stereoisomers pairs
~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    poetry run smartclass querywikidata -q src/smartclass/data/queries/chemical_entities_inchi.rq -o scratch/wikidata_chemical_entities_inchis.tsv
    poetry run smartclass querywikidata -q src/smartclass/data/queries/chemicals_stereosiomer_of.rq -o scratch/wikidata_chemicals_stereoisomer_of.tsv
    poetry run python3 notebooks/pair_stereoisomers.py

Command Line Interface
~~~~~~~~~~~~~~~~~~~~~~
The smartclass command line tool is automatically installed. It can
be used from the shell with the ``--help`` flag to show all subcommands:

.. code-block:: sh

    poetry run smartclass --help


TODO show the most useful thing the CLI does! The CLI will have documentation auto-generated by ``sphinx``.


👐 Contributing
---------------
Contributions, whether filing an issue, making a pull request, or forking, are appreciated.
See `CONTRIBUTING.md <https://github.com/zamboni-lab/smartclass/blob/main/.github/CONTRIBUTING.md>`_ for more information on getting involved.

👋 Attribution
--------------

⚖️ License
~~~~~~~~~~
The code in this package is licensed under the MIT License.

..
 📖 Citation
 ~~~~~~~~~~~
 Citation goes here!

..
 🎁 Support
 ~~~~~~~~~~
 This project has been supported by the following organizations (in alphabetical order):
 - [Harvard Program in Therapeutic Science - Laboratory of Systems Pharmacology](https://hits.harvard.edu/the-program/laboratory-of-systems-pharmacology/)

..
 💰 Funding
 ~~~~~~~~~~
 This project has been supported by the following grants:
 - [Harvard Program in Therapeutic Science - Laboratory of Systems Pharmacology](https://hits.harvard.edu/the-program/laboratory-of-systems-pharmacology/)

🍪 Cookiecutter
~~~~~~~~~~~~~~~
This package was created with `@audreyfeldroy <https://github.com/audreyfeldroy>`_'s
`cookiecutter <https://github.com/cookiecutter/cookiecutter>`_ package using `@cthoyt <https://github.com/cthoyt>`_'s
`cookiecutter-snekpack <https://github.com/cthoyt/cookiecutter-snekpack>`_ template.

🛠️ For Developers
-----------------
Development Installation
~~~~~~~~~~~~~~~~~~~~~~~~
To install in development mode, use the following:

.. code-block:: sh

    git clone git+https://github.com/zamboni-lab/smartclass.git
    cd smartclass
    pip install -e .

🥼 Testing
~~~~~~~~~~
After cloning the repository and installing ``tox`` with ``pip install tox``, the unit tests in the ``tests/`` folder can be
run reproducibly with:

.. code-block:: sh

    tox

Additionally, these tests are automatically re-run with each commit in a `GitHub Action <https://github.com/zamboni-lab/smartclass/actions?query=workflow%3ACI>`_.

📖 Building the Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The documentation can be built locally using the following:

.. code-block:: sh

    tox -e docs
    open docs/build/html/index.html

The documentation automatically installs the package as well as the ``docs``
extra specified in the `setup.cfg <setup.cfg>`_. ``sphinx`` plugins
like ``texext`` can be added there. Additionally, they need to be added to the
``extensions`` list in `docs/source/conf.py <docs/source/conf.py>`_.

📦 Making a Release
~~~~~~~~~~~~~~~~~~~
After installing the package in development mode and installing
``tox`` with ``pip install tox``, the commands for making a new release are contained within the ``finish`` environment
in ``tox.ini``. Run the following from the shell:

.. code-block:: sh

    tox -e finish

This script does the following:

1. Uses `Bump2Version <https://github.com/c4urself/bump2version>`_ to switch the version number in the ``setup.cfg``,
   ``src/smartclass/version.py``, and `docs/source/conf.py <docs/source/conf.py>`_ to not have the ``-dev`` suffix
2. Packages the code in both a tar archive and a wheel using `build <https://github.com/pypa/build>`_
3. Uploads to PyPI using `twine <https://github.com/pypa/twine>`_. Be sure to have a ``.pypirc`` file configured to avoid the need for manual input at this
   step
4. Push to GitHub. You'll need to make a release going with the commit where the version was bumped.
5. Bump the version to the next patch. If you made big changes and want to bump the version by minor, you can
   use ``tox -e bumpversion -- minor`` after.

Logo
~~~~
The smartclass `logo <https://github.com/smartclass/smartclass-art>`_ was designed by `TODO <https://github.com/TODO>`_.

.. |build| image:: https://github.com/smartclass/smartclass/workflows/Tests/badge.svg
    :target: https://github.com/smartclass/smartclass/actions
    :alt: Build Status

.. |coverage| image:: https://codecov.io/gh/smartclass/smartclass/coverage.svg?branch=develop
    :target: https://codecov.io/gh/smartclass/smartclass/branch/develop
    :alt: Development Coverage Status

.. |documentation| image:: https://readthedocs.org/projects/smartclass/badge/?version=latest
    :target: http://smartclass.readthedocs.io/en/latest/
    :alt: Development Documentation Status

.. |climate| image:: https://codeclimate.com/github/smartclass/smartclass/badges/gpa.svg
    :target: https://codeclimate.com/github/smartclass/smartclass
    :alt: Code Climate

.. |python_versions| image:: https://img.shields.io/pypi/pyversions/smartclass.svg
    :target: https://pypi.python.org/pypi/smartclass
    :alt: Stable Supported Python Versions

.. |pypi_version| image:: https://img.shields.io/pypi/v/smartclass.svg
    :target: https://pypi.python.org/pypi/smartclass
    :alt: Current version on PyPI

.. |pypi_license| image:: https://img.shields.io/pypi/l/smartclass.svg
    :target: https://github.com/smartclass/smartclass/blob/main/LICENSE
    :alt: MIT License

.. |zenodo| image:: https://zenodo.org/badge/TODO.svg
    :target: https://zenodo.org/badge/latestdoi/TODO

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: Code style: black
