from __future__ import annotations

import polars as pl
from rdkit import Chem

from smartclass.chem.helpers.remove_layers_from_inchi import remove_layers_from_inchi


def match_smiles(smiles_1: str, smiles_2: str) -> bool:
    mol_1 = Chem.MolFromSmiles(smiles_1)
    mol_2 = Chem.MolFromSmiles(smiles_2)
    if mol_1 is not None and mol_2 is not None:
        return mol_1.HasSubstructMatch(mol_2, useChirality=True)
    else:
        return False


LAYERS = ["b", "t", "m", "s"]

stereoisomers_canonical_df = pl.read_csv(
    "scratch/wikidata_stereoisomers_smiles_canonical_inchi.tsv",
    separator="\t",
)
stereoisomers_isomeric_df = pl.read_csv(
    "scratch/wikidata_stereoisomers_smiles_isomeric_inchi.tsv",
    separator="\t",
)
tautomers_df = pl.read_csv(
    "scratch/wikidata_chemicals_tautomer_of.tsv",
    separator="\t",
)
chemicals_df = pl.read_csv(
    "scratch/wikidata_chemical_entities_smiles_inchi.tsv",
    separator="\t",
)
classes_df = pl.read_csv(
    "scratch/wikidata_chemicals_classes.tsv",
    separator="\t",
).drop("instance")

# Combine both stereoisomers and chemicals dfs
stereoisomers_canonical_df = stereoisomers_canonical_df.select(sorted(stereoisomers_canonical_df.columns))
stereoisomers_isomeric_df = stereoisomers_isomeric_df.select(sorted(stereoisomers_isomeric_df.columns))
chemicals_df = chemicals_df.select(sorted(chemicals_df.columns))
stereoisomers_df = pl.concat([stereoisomers_canonical_df, stereoisomers_isomeric_df], rechunk=True).unique()

stereoisomers_df = pl.concat([stereoisomers_canonical_df, stereoisomers_isomeric_df], rechunk=True).unique()
chemicals_df = pl.concat([chemicals_df, stereoisomers_isomeric_df], rechunk=True).unique()

# Remove InChI stereo layers
stereoisomers_df = stereoisomers_df.with_columns(
    pl.col("inchi")
    .map_elements(lambda x: remove_layers_from_inchi(x, layers=LAYERS), return_dtype=str)
    .alias("inchi_no_stereo"),
).drop("inchi")
chemicals_df = chemicals_df.with_columns(
    pl.col("inchi")
    .map_elements(lambda x: remove_layers_from_inchi(x, layers=LAYERS), return_dtype=str)
    .alias("inchi_no_stereo"),
).drop("inchi")

# Merge stereoisomers_df with chemicals_df based on the "inchi_no_stereo" column
merged_df = chemicals_df.join(
    stereoisomers_df,
    left_on="inchi_no_stereo",
    right_on="inchi_no_stereo",
    how="inner",
).drop("inchi_no_stereo")
print(merged_df)

# Merge merged_df with tautomers_df to avoid tautomers
merged_df = merged_df.join(tautomers_df, left_on="structure", right_on="structure_1", how="anti")
print(merged_df)

merged_df = (
    merged_df.with_columns(
        pl.struct(["smiles", "smiles_right"])
        .map_elements(
            lambda cols: match_smiles(cols["smiles"], cols["smiles_right"]),
            return_dtype=bool,
        )
        .alias("match"),
    )
    .filter(pl.col("smiles") != pl.col("smiles_right"))
    .drop("smiles", "smiles_right")
    .filter(pl.col("match"))
    .drop("match")
)
print(merged_df)

# TODO check if this part could be improved

# Merge with classes_df based on the "structure" columns
final_merged_df = merged_df.join(classes_df, left_on="structure", right_on="structure").join(
    classes_df,
    left_on="structure_right",
    right_on="structure",
    how="inner",
)
print(final_merged_df)

# Filter rows where classes match
matched_classes_df = final_merged_df.filter(pl.col("class_right") == pl.col("class")).filter(
    pl.col("structure") != pl.col("structure_right"),
)
matched_classes_df_2 = final_merged_df.filter(pl.col("class") == "Q11173").filter(
    pl.col("structure") != pl.col("structure_right"),
)

# Renaming
matched_classes_df = matched_classes_df.select([
    pl.col("structure").alias("qid"),
    pl.col("structure_right").alias("P279"),
    pl.lit("Q113907573").alias("S887"),
    pl.col("class").alias("-P279"),
]).unique()
matched_classes_df_2 = matched_classes_df_2.select([
    pl.col("structure").alias("qid"),
    pl.col("structure_right").alias("P279"),
    pl.lit("Q113907573").alias("S887"),
    pl.col("class").alias("-P279"),
]).unique()

# Removing wikidata url
matched_classes_df = matched_classes_df.with_columns(
    qid=pl.col("qid").str.replace("http://www.wikidata.org/entity/(.)", "${1}"),
    P279=pl.col("P279").str.replace("http://www.wikidata.org/entity/(.)", "${1}"),
    P268=pl.col("-P279").str.replace("http://www.wikidata.org/entity/(.)", "${1}"),
)
matched_classes_df_2 = matched_classes_df_2.with_columns(
    qid=pl.col("qid").str.replace("http://www.wikidata.org/entity/(.)", "${1}"),
    P279=pl.col("P279").str.replace("http://www.wikidata.org/entity/(.)", "${1}"),
    P268=pl.col("-P279").str.replace("http://www.wikidata.org/entity/(.)", "${1}"),
)

# Renaming (again)
matched_classes_df = matched_classes_df.select([
    pl.col("qid"),
    pl.col("P279"),
    pl.col("S887"),
    pl.col("P268").alias("-P279"),
]).unique()
matched_classes_df_2 = matched_classes_df_2.select([
    pl.col("qid"),
    pl.col("P279"),
    pl.col("S887"),
    pl.col("P268").alias("-P279"),
]).unique()

# Export
print(matched_classes_df)
print(matched_classes_df_2)
matched_classes_df.write_csv("scratch/maintenance_improve_subclasses.csv")
matched_classes_df_2.write_csv("scratch/maintenance_improve_subclasses_2.csv")
