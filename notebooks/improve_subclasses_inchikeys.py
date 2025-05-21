from __future__ import annotations

import polars as pl

# Paths
file_inchis = "scratch/wikidata_chemical_entities_inchis.tsv"
file_inchikeys = "scratch/wikidata_chemicals_inchikeys.tsv"
file_classes = "scratch/wikidata_chemicals_classes.tsv"

# Import
df_inchis = pl.read_csv(file_inchis, separator="\t")
df_inchikeys = pl.read_csv(file_inchikeys, separator="\t")
df_subclass = pl.read_csv(file_classes, separator="\t").drop("instance")

# Filter out isotopes
df_inchis_filtered = df_inchis.filter(~pl.col("inchi").str.contains("/i")).drop("inchi")

# Create a new column with the first 14 characters of the inchikey
df_inchikeys = df_inchikeys.with_columns(
    [pl.col("inchikey").str.slice(0, 14).alias("short_inchikey")]
)

# Perform a self-join on the short_inchikey column
df_joined = df_inchikeys.join(
    df_inchikeys, on="short_inchikey", how="inner", suffix="_right"
).drop("short_inchikey")

# Filter rows where the first inchikey contains "UHFFFAOYSA-N"
df_filtered = df_joined.filter(pl.col("inchikey").str.contains("UHFFFAOYSA-N"))

# Ensure that inchikey and inchikey_right are different
df_filtered = df_filtered.filter(pl.col("inchikey") != pl.col("inchikey_right")).drop(
    "inchikey", "inchikey_right"
)

# Join the filtered result with the InChI DataFrame on the InChI key
df_filtered = df_filtered.join(
    df_inchis_filtered, left_on="structure_right", right_on="structure", how="inner"
)

# Join the filtered result with the subclass DataFrame on 'structure'
df_filtered = df_filtered.join(
    df_subclass, left_on="structure", right_on="structure", how="inner"
)
df_filtered = df_filtered.join(
    df_subclass,
    left_on="structure_right",
    right_on="structure",
    how="inner",
    suffix="_right",
)

# Ensure the pair is not already modeled
df_final = df_filtered.filter(pl.col("structure") != pl.col("class_right"))

# Create new DataFrames based on specific class conditions
df_1 = df_final.filter(
    (pl.col("class") == "Q11173") & (pl.col("class_right") != "Q11173")
).select(
    [
        pl.col("structure").alias("qid"),
        pl.col("class").alias("-P279"),
        pl.col("class_right").alias("P279"),
    ]
)

df_2 = (
    df_final.filter(pl.col("class_right") == "Q11173")
    .select(
        [
            pl.col("structure_right").alias("qid"),
            pl.col("class_right").alias("-P279"),
            pl.col("structure").alias("P279"),
        ]
    )
    .with_columns([pl.lit("Q113993940").alias("S887")])
)

# Ensure structure is not the same as P279 (might be the case due to multiple InChIKeys)
df_1 = df_1.filter(pl.col("qid") != pl.col("P279"))
df_2 = df_2.filter(pl.col("qid") != pl.col("P279"))

print(df_1)
print(df_2)

# Export
df_1.write_csv("scratch/maintenance_improve_subclasses_inchikeys_1.csv")
df_2.write_csv("scratch/maintenance_improve_subclasses_inchikeys_2.csv")
