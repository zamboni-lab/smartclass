"""Check stereochemistry."""

from __future__ import annotations

import logging

from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import GetStereoisomerCount, StereoEnumerationOptions

__all__ = [
    "check_missing_stereochemistry",
]


def check_missing_stereochemistry(smiles: str, use_legacy: bool = False) -> bool | None:
    """
    Check stereochemistry.

    :param smiles: SMILES.
    :type smiles: str

    :param use_legacy: Flag indicating to use legacy or not.
    :type use_legacy: bool

    :returns: flag.
    :rtype: Union[bool, None]
    """
    # See https://github.com/rdkit/rdkit/issues/6217
    Chem.SetAllowNontetrahedralChirality(False)
    # Chem.SetUseLegacyStereoPerception(use_legacy)
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        Chem.FindPotentialStereoBonds(mol, cleanIt=False)
        si = Chem.FindPotentialStereo(mol, cleanIt=False, flagPossible=True)
        specified_count = 0
        unspecified_count = 0
        unknown_count = 0
        for element in si:
            # logging.warning(
            #     f"Type: {element.type}, Which: {element.centeredOn},
            # Specified: {element.specified}, Descriptor: {element.descriptor}"
            #     )
            if element.specified:
                if str(element.specified) == "Unknown":
                    unknown_count += 1
                else:
                    specified_count += 1
            else:
                unspecified_count += 1

        # TODO this is almost working for now
        # logging.warning(specified_count)
        # logging.warning(unspecified_count)
        # logging.warning(unknown_count)
        if unspecified_count + unknown_count == 0:
            # Fully defined
            return False
        elif specified_count > 0:
            # TODO additional check leading to almost everything OK
            return count_multiple_stereoisomers(mol)
        else:
            # TODO additional check leading to almost everything OK
            return count_multiple_stereoisomers(mol)
    else:
        return None


def count_multiple_stereoisomers(
    mol,
):
    """
    Count multiple stereoisomers.

    :param mol: Mol.
    :type mol: TODO

    :returns: flag.
    :rtype: bool
    """
    if (
        GetStereoisomerCount(
            mol,
            options=StereoEnumerationOptions(
                # onlyUnassigned=False,
                # onlyStereoGroups=True,
                tryEmbedding=True,
            ),
        )
        == 1
    ):
        # logging.warning(GetStereoisomerCount(mol, options=options))
        return False
    else:
        # logging.warning(GetStereoisomerCount(mol, options=options))
        return True


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Check stereochemistry")
    parser.add_argument("--smiles", type=str, help="SMILES string (optional)")
    parser.add_argument(
        "--use-legacy",
        action="store_true",
        help="Use legacy stereochemistry perception (default: False)",
    )
    args = parser.parse_args()

    use_legacy = args.use_legacy

    def check_expectations(smiles, expectation, use_legacy=use_legacy):
        """Check expectations."""
        if expectation is None:
            logging.info("\033[33m" + smiles)
        elif check_missing_stereochemistry(smiles, use_legacy) == expectation:
            logging.info("\033[32m" + smiles)
        else:
            logging.info("\033[31m" + smiles)
        logging.info("\033[0m")

    smiles_undefined_to_test = [
        "CCOC(=O)C=C(CC(=O)O)C(=O)O",
        "C1CCCCCC=CCCCC1",
        "CC(CCCF)(F)Cl",
        "CC(F)C(F)(Cl)C(F)C",
        "C[C@H]1C=C[C@@H](O)C=C(Cl)CCOC(=O)C2=CC[C@H]3CCC[C@]4(CCC[C@H]14)N3C2",
        "COc1cc(C(=O)O[C@@H]2CC(O)(C(=O)O)C[C@@H](OC(=O)c3ccc(O)c(OC)c3)C2O)ccc1O",
        "C=CCSS/C=C/CS(=O)CC=C",
    ]

    smiles_defined_to_test = [
        "C[C@H](CN1C=NC2=C1N=CN=C2N)OCP(=O)(O)O",
        "C1CCCCC/C=C/CCCC1",
        "Cc1cc(C(=C2C=CC(=Nc3ccc(S(=O)(=O)[O-])cc3)C=C2)c2ccc(Nc3ccc(S(=O)(=O)[O-])cc3)cc2)cc(S(=O)(=O)O)c1N",
        "CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)c3ccc(cc3)N=C=S)C(=O)OC",
        "O=[As]O[As]=O",
        "CN[C@H](C)[C@H](OB1C2CCCC(CCC2)[C@H]1[Si](C)(C)C)c1ccccc1",
        "C1CN2C(=N)N([C@H]([C@H]3[C@]2(C1(O)O)N=C(N3)N)COC(=O)N)O",
        "CC12CCC(CC1)C(C)(C)O2",
        "COC(=O)C1CCC(C)CC1",
        "COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O",
    ]

    smiles_unclear_to_test = [
        "C1CN2CCC1[C@@]3(C2)CC4=C(O3)N=CC=C4",
        "C1c2cccnc2O[C@]11CN2CCC1CC2",
        "C1CN2CCC1[C@H](C2)OC(=O)N3CCC4=CC=CC=C4[C@@H]3C5=CC=CC=C5",
    ]

    logging.info("Undefined:")
    for smiles in smiles_undefined_to_test:
        check_expectations(smiles, True)
    logging.info("")

    logging.info("Defined:")
    for smiles in smiles_defined_to_test:
        check_expectations(smiles, False)
    logging.info("")

    logging.info("Unclear:")
    for smiles in smiles_unclear_to_test:
        check_expectations(smiles, None)
