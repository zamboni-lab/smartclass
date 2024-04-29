from rdkit import Chem

m = Chem.MolFromSmiles(
    "O=C1C2=C(O)C=C(OC3OC(COC4OC(C)C(O)C(O)C4O)C(O)C(O)C3O)C=C2OC(C5=CC=C(OC)C(O)=C5)C1"
)

s1 = Chem.MolFromSmiles(
    "COc1ccc([C@@H]2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@H](C)[C@@H](O)[C@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O"
)

s2 = Chem.MolFromSmiles(
    "COc1ccc([C@@H]2CC(=O)c3c(O)cc(O[C@@H]4O[C@H](CO[C@@H]5O[C@@H](C)[C@H](O)[C@H](O)[C@H]5O)[C@@H](O)[C@H](O)[C@H]4O)cc3O2)cc1O"
)

s3 = Chem.MolFromSmiles(
    "COc1ccc(cc1O)[C@@H]1CC(=O)c2c(O)cc(O[C@@H]3OC(CO[C@@H]4OC(C)[C@H](O)[C@H](O)C4O)[C@@H](O)[C@H](O)C3O)cc2O1"
)


print(s1.HasSubstructMatch(m, useChirality=True))
print(s2.HasSubstructMatch(m, useChirality=True))
print(s1.HasSubstructMatch(s2, useChirality=True))
print(s1.HasSubstructMatch(s3, useChirality=True))
print(s2.HasSubstructMatch(s3, useChirality=True))
