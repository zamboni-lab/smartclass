"""Fix InChI tautomerization."""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import AllChem

__all__ = [
    "fix_inchi_tautomerization",
]


def fix_inchi_tautomerization(smiles: str) -> str:
    """
    Fix InChI tautomerization.

    :param smiles: SMILES.
    :type smiles: str

    :returns: SMILES.
    :rtype: str
    """
    #  Courtesy of Richy Leroy (via Jean-Marc Nuzillard)
    ####################################################
    transformations = {
        "Secondary Iminol": {
            "SMARTS": "[NH0]=C([OH])([!O])",
            "Reaction": "[C:1]([OH:2])=[NH0:3]>>[C:1](=[OH0:2])[NH:3]",
        },
        "Primary Iminol": {
            "SMARTS": "[NH1]=C([OH])([!O])",
            "Reaction": "[CD3:1]([OH:2])=[NH:3]>>[CD3:1](=[OH0D1:2])[NH2:3]",
        },
        "Secondary Carbamate": {
            "SMARTS": "[A&!H][NH0]=C([OH])O",
            "Reaction": "[O:4][CD3:1]([OH:2])=[NH0:3]>>[O:4][CD3:1](=[OH0D1:2])[NH:3]",
        },
        "Primary Carbamate": {
            "SMARTS": "[NH]=C([OH])O",
            "Reaction": "[O:2][CD3:3]([OH1:4])=[NH:5]>>[O:2][CD3:3](=[OH0D1:4])[NH2:5]",
        },
        "N-Formyl": {
            "SMARTS": "[CH]([OH])=N",
            "Reaction": "[CH1:1]([OH:2])=[N:3]>>[CH1:1](=[OH0D1:2])[NH:3]",
        },
        "Enol": {
            "SMARTS": "[C&!c]([!CX3&!O])([!CX3&!O])=[C&!c]([OH])([!N])",
            "Reaction": "[!c&C:1]([OH:2])=[!c&C:3]>>[C:1](=[OH0:2])[CH:3]",
        },
        "Enethiol": {
            "SMARTS": "[C&!c]([!CX3&!SHX2])([!CX3&!SHX2])=[C&!c]([SHX2])",
            "Reaction": "[!c&C:1]([SH:2])=[!c&C:3]>>[C:1](=[SH0:2])[CH:3]",
        },
        # "Enethiol_n": {
        #     "SMARTS": "[N&!n]([!CX3&!SHX2])([!CX3&!SHX2])=[C&!c]([SHX2])",
        #     "Reaction": "[N&!n:1]([!CX3&!SHX2:2])([!CX3&!SHX2:3])=[C&!c:4]([SHX2])>>[C:4](=[SH0:2])[CH:3]",
        # },
    }
    ####################################################

    # TODO fix this later
    if "." in smiles:
        return smiles
    # TODO deuterated compounds
    if "2H" in smiles:
        return smiles

    else:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            changes_count = 0
            changed = True
            while changed:
                changed = False
                for _, data in transformations.items():
                    target = Chem.MolFromSmarts(data["SMARTS"])
                    if mol.HasSubstructMatch(target):
                        rxn = AllChem.ReactionFromSmarts(data["Reaction"])
                        reactants = (mol,)
                        mol = rxn.RunReactants(reactants)
                        mol = mol[0][0]
                        Chem.SanitizeMol(mol)
                        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
                        changed = True
                        changes_count += 1
                        break
            if changes_count > 0:
                return Chem.MolToSmiles(mol)
            else:
                return smiles
        else:
            return smiles


if __name__ == "__main__":
    smiles_to_test = [
        "CC(C)C1=CC[C@H]2C(=C1)CC[C@@H]3[C@@]2(CCC[C@@]3(C)C(=O)O)C",
        "O=C(O)CCC(N=C(O)C(N=C(O)C(N=C(O)C(N)CC=1C=CC=CC1)CC(C)C)CC2=CN=CN2)C(O)=NC(C(O)=NC3C(O)=NCCCC3)C(C)CC",
        "CC[C@H](C)C(N=C(O)[C@@H](CCC(=O)O)N=C(O)[C@H](Cc1cnc[nH]1)N=C(O)[C@H](CC(C)C)N=C(O)[C@H](N)Cc1ccccc1)C(O)=N[C@H]1CCCCN=C1O",  # noqa:E501
    ]

    for smiles in smiles_to_test:
        fix_inchi_tautomerization(smiles)
