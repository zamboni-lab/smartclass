"""TODO.

TODO.
"""

# MISSING DEF
# __all__ = [
#     "measure_mhfp",
# ]

# from itertools import combinations

# from rdkit import Chem
# from rdkit.Chem import rdMHFPFingerprint

# from smartclass.io import export_results

# import csv

# __all__ = [
#     "load_substructures_smarts",
# ]


# def load_substructures_smarts(substructures: str) -> dict:
#     """
#     Description.
#
#     :param substructures: Description.
#     :type substructures: str
#
#     :return: Description.
#     :rtype: dict
#     """
#     substructures_dict = {}
#     try:
#         with open(substructures, "r") as f:
#             reader = csv.DictReader(f, delimiter="\t")
#             for row in reader:
#                 class_name = row["class_id"]
#                 smarts = row["class_smarts"]
#                 if class_name and smarts:
#                     substructures_dict[class_name] = smarts
#     except FileNotFoundError:
#         print(f"File not found: {substructures}")
#     except Exception as e:
#         print(f"An error occurred while loading substructures: {str(e)}")
#
#     return substructures_dict

# # TODO
# input = "smartclass/data/mia_smarts.tsv"
# output = "distances_mia.tsv"

# substructures_dict = load_substructures_smarts(substructures = input)
# print(substructures_dict)
# substructures = list(substructures_dict.values())

# encoder = rdMHFPFingerprint.MHFPEncoder(2048)

# mols = [Chem.MolFromSmarts(smarts) for smarts in substructures]
# fps = [
#     encoder.EncodeMol(mol, kekulize=False) for mol in mols
# ]  # Some SMARTS in Wikidata lead to errors if `kekulize = True`, investigate

# results = [
#     {
#         "query": list(substructures_dict.keys())[i],
#         "target": list(substructures_dict.keys())[j],
#         "distance": rdMHFPFingerprint.MHFPEncoder.Distance(fps[i], fps[j]),
#     }
#     for i, j in combinations(range(len(substructures)), 2)
# ]
# print(results)

# export_results(output=output, results=results)
