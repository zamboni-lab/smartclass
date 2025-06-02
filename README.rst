Smartclass
==========
Smartclass classifies structures using SMARTS.

🚀 Installation
---------------
.. code-block:: sh

    uv sync


💪 Getting Started
------------------
A mini notebook is available.

.. code-block:: sh

    uv run quarto render notebooks/smartclass.qmd

🌟 Main steps
-------------

The main steps are briefly shown below.

Get defined chemical classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    uv run smartclass querywikidata -q src/smartclass/data/queries/classes_smarts.rq -o scratch/wikidata_classes_smarts.tsv


Canonicalize them
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    uv run notebooks/canonicalize_smarts.py

Get some SMILES to classify
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    uv run smartclass querywikidata -q src/smartclass/data/queries/chemicals_smiles_canonical.rq -o scratch/wikidata_chemicals_smiles.tsv
    uv run smartclass querywikidata -q src/smartclass/data/queries/classes_fake_smiles.rq -o scratch/wikidata_classes_smiles.tsv

Classify compounds
~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    uv run smartclass searchclasses --help
    # uv run smartclass searchclasses
    uv run smartclass searchclasses -s "O=C1OC2CC3C(C=C(OC)C(=O)C3(C)C4C(=O)C(OC)=C(C)C(C1)C24C)C" -c scratch/wikidata_classes_smarts.tsv --verbose
    uv run smartclass searchclasses -s "O=C1OC2CC3C(C=C(OC)C(=O)C3(C)C4C(=O)C(OC)=C(C)C(C1)C24C)C" -c scratch/wikidata_classes_smarts.tsv --closest-only False --verbose
    uv run smartclass searchclasses -i src/smartclass/data/bitter_smiles.tsv -c scratch/wikidata_classes_smarts.tsv
    # uv run smartclass searchclasses -i scratch/wikidata_chemicals_smiles.tsv -c scratch/wikidata_classes_smarts.tsv --closest-only False
    # TODO Improve classes taxonomy
    # uv run smartclass searchclasses -i scratch/wikidata_classes_smiles.tsv -c scratch/wikidata_classes_smarts.tsv --closest-only False

🤯 Future steps (not fully available now)
-----------------------------------------

Measure substructures distances
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    uv run python3 smartclass/measure_mhfp.py -i scratch/wikidata_classes_smarts.tsv -o scratch/distances_wikidata.tsv
    uv run python3 smartclass/measure_mhfp.py -i data/mia_smarts.tsv -o scratch/distances_mia.tsv

🛠️ Wikidata maintenance
-----------------------

Generic maintenance
~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_canonical_no_formula.rq -t transform_smiles_to_formula -o scratch/formulas_canonical.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_no_smiles_isomeric_no_mass.rq -t transform_inchi_to_mass -o scratch/masses_inchi.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_isomeric_no_mass.rq -t transform_smiles_to_mass -o scratch/masses_isomeric.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_isomeric_no_inchi.rq -t transform_smiles_to_inchi -o scratch/inchis_isomeric.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_isomeric_smiles_inchi_no_stereo.rq -t transform_smiles_to_inchi -o scratch/incorrect_inchis.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_isomeric_no_formula.rq -t transform_smiles_to_formula -o scratch/formulas_isomeric.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_isomeric_no_canonical.rq -t transform_smiles_i_to_smiles_c -o scratch/smiles_isomeric.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_no_smiles_isomeric_inchi.rq -t transform_inchi_to_smiles_isomeric -o scratch/smiles_i_inchi.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_inchi_no_isomeric_smiles.rq -t transform_inchi_to_smiles_isomeric -o scratch/smiles_i_inchi_2.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_inchi_no_canonical_smiles.rq -t transform_inchi_to_smiles_canonical -o scratch/smiles_c_inchi.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_inchi_no_inchikey.rq -t transform_inchi_to_inchikey -o scratch/inchikeys.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_formula_malformed.rq -t transform_formula_to_formula -o scratch/formulas_malformed.csv
    # WDQS is not parsing the REGEX well
    # uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_formula_malformed.rq -t transform_formula_to_formula -o scratch/formulas_malformed.csv -u https://qlever.cs.uni-freiburg.de/api/wikidata
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_smiles_isomeric_incorrect_mass.rq -t transform_smiles_mass_to_smiles_mass -o scratch/masses_incorrect_isomeric.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/maintenance_no_smiles_isomeric_incorrect_mass.rq -t transform_inchi_mass_to_inchi_mass -o scratch/masses_incorrect_inchi.csv

    uv run smartclass querywikidata -q src/smartclass/data/queries/chemicals_smiles_canonical.rq -t check_smiles -o scratch/smiles_c_invalid.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/chemicals_smiles_isomeric.rq -t check_smiles -o scratch/smiles_i_invalid.csv

    uv run smartclass querywikidata -q src/smartclass/data/queries/chemicals_smiles_canonical_no_ref.rq -t transform_smiles_c_to_smiles_c_tauto -o scratch/smiles_c_tauto.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/chemicals_smiles_isomeric_no_ref.rq -t transform_smiles_i_to_smiles_i -o scratch/smiles_i.csv
    uv run smartclass querywikidata -q src/smartclass/data/queries/chemicals_smiles_isomeric_no_ref.rq -t transform_smiles_i_to_smiles_i_tauto -o scratch/smiles_i_tauto.csv
    # not working for now
    # uv run smartclass querywikidata -q src/smartclass/data/queries/stereoisomers_smiles_isomeric_inchi.rq -t transform_stereoisomers_to_entities -o scratch/stereo_i_to_entities.csv
    # uv run smartclass querywikidata -q src/smartclass/data/queries/stereoisomers_smiles_canonical_no_isomeric.rq -t transform_stereoisomers_to_entities -o scratch/stereo_c_to_entities.csv
    # uv run smartclass querywikidata -q src/smartclass/data/queries/chemical_entities_smiles_isomeric_inchi.rq -t transform_entities_to_stereoisomers -o scratch/entities_i_to_stereo.csv
    # uv run smartclass querywikidata -q src/smartclass/data/queries/chemical_entities_smiles_canonical_no_isomeric.rq -t transform_entities_to_stereoisomers -o scratch/entities_c_to_stereo.csv

Improve current classes
~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    uv run smartclass querywikidata -q src/smartclass/data/queries/chemicals_inchikey.rq -o scratch/wikidata_chemicals_inchikeys.tsv
    uv run smartclass querywikidata -q src/smartclass/data/queries/chemicals_classes.rq -o scratch/wikidata_chemicals_classes.tsv
    uv run smartclass querywikidata -q src/smartclass/data/queries/chemical_entities_smiles_isomeric_inchi.rq -o scratch/wikidata_chemical_entities_smiles_inchi.tsv
    uv run smartclass querywikidata -q src/smartclass/data/queries/stereoisomers_smiles_isomeric_inchi.rq -o scratch/wikidata_stereoisomers_smiles_isomeric_inchi.tsv
    uv run smartclass querywikidata -q src/smartclass/data/queries/stereoisomers_smiles_canonical_no_isomeric_inchi.rq -o scratch/wikidata_stereoisomers_smiles_canonical_inchi.tsv
    uv run smartclass querywikidata -q src/smartclass/data/queries/chemicals_tautomer_of.rq -o scratch/wikidata_chemicals_tautomer_of.tsv
    uv run python3 notebooks/improve_classes.py
    uv run python3 notebooks/improve_subclasses_inchikeys.py
    uv run python3 src/smartclass/helpers/split_csv.py

Add stereoisomers pairs
~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: sh

    uv run smartclass querywikidata -q src/smartclass/data/queries/chemical_entities_inchi.rq -o scratch/wikidata_chemical_entities_inchis.tsv
    uv run smartclass querywikidata -q src/smartclass/data/queries/chemicals_stereosiomer_of.rq -o scratch/wikidata_chemicals_stereoisomer_of.tsv
    uv run python3 notebooks/pair_stereoisomers.py

🖥 Command Line Interface
-------------------------
The smartclass command line tool is automatically installed. It can
be used from the shell with the ``--help`` flag to show all subcommands:

.. code-block:: sh

    uv run smartclass --help

👐 Contributing
---------------
Contributions, whether filing an issue, making a pull request, or forking, are appreciated.
See `CONTRIBUTING.md <https://github.com/zamboni-lab/smartclass/blob/main/.github/CONTRIBUTING.md>`_ for more information on getting involved.

👋 Attribution
--------------
A lot of inspiration and initial data has been taken by the huge work done by `@rwst <https://github.com/rwst>`_ with `YACCL <https://github.com/rwst/yaccl>`_.

⚖️ License
~~~~~~~~~~
See `LICENSE <https://github.com/zamboni-lab/smartclass/blob/main/LICENSE>`_

..
 📖 Citation
 ~~~~~~~~~~~
 Citation goes here!

..
 🎁 Support
 ~~~~~~~~~~
 This project has been supported by the following organizations (in alphabetical order):
 - [TODO](TODO)

..
 💰 Funding
 ~~~~~~~~~~
 This project has been supported by the following grants:
 - [TODO](TODO)

🍪 Cookiecutter
~~~~~~~~~~~~~~~
This package was created with `@audreyfeldroy <https://github.com/audreyfeldroy>`_'s
`cookiecutter <https://github.com/cookiecutter/cookiecutter>`_ package using `@cthoyt <https://github.com/cthoyt>`_'s
`cookiecutter-snekpack <https://github.com/cthoyt/cookiecutter-snekpack>`_ template.
