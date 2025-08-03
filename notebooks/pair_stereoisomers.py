from __future__ import annotations

import polars as pl

from smartclass.chem.helpers.check_layers_from_inchi import check_layers_from_inchi
from smartclass.chem.helpers.remove_layers_from_inchi import remove_layers_from_inchi

LAYERS = ["b", "t", "m", "s"]

stereo_df = pl.read_csv(
    "scratch/wikidata_chemicals_stereoisomer_of.tsv",
    separator="\t",
).select([
    pl.col("structure_1").alias("structure"),
    pl.col("structure_2").alias("structure_right"),
])

inchis_df = pl.read_csv(
    "scratch/wikidata_chemical_entities_inchis.tsv",
    separator="\t",
)

# Filter out InchIs that contain the specified layers
inchis_df = inchis_df.filter(
    pl.col("inchi").map_elements(
        lambda x: check_layers_from_inchi(x, layers=LAYERS),
        return_dtype=bool,
    ),
)

# Remove InChI stereo layers
inchis_df = inchis_df.with_columns(
    pl.col("inchi")
    .map_elements(
        lambda x: remove_layers_from_inchi(x, layers=LAYERS),
        return_dtype=str,
    )
    .alias("inchi_no_stereo"),
)

# Self-join on the inchikey_prefix column
merged_df = inchis_df.join(
    inchis_df.clone(),
    on=[
        "inchi_no_stereo",
    ],
).filter(pl.col("structure") != pl.col("structure_right"))

# Remove pairs already present
add_df = merged_df.join(
    stereo_df,
    on=[
        "structure",
        "structure_right",
    ],
    how="anti",
).filter(pl.col("structure") != pl.col("structure_right"))

# Remove wrong pairs
remove_df = stereo_df.join(
    merged_df,
    on=[
        "structure",
        "structure_right",
    ],
    how="anti",
).filter(pl.col("structure") != pl.col("structure_right"))

# Renaming
add_df = (
    add_df.select([
        pl.col("structure").alias("qid"),
        pl.col("structure_right").alias("P3364"),
        pl.lit("Q123137214").alias("S887"),
    ])
    .unique()
    .with_columns(
        qid=pl.col("qid").str.replace("http://www.wikidata.org/entity/(.)", "${1}"),
        P3364=pl.col("P3364").str.replace("http://www.wikidata.org/entity/(.)", "${1}"),
        S887=pl.col("S887").str.replace("http://www.wikidata.org/entity/(.)", "${1}"),
    )
    .select([
        pl.col("qid"),
        pl.col("P3364"),
        pl.col("S887"),
    ])
    .unique()
)

remove_df = (
    remove_df.select([
        pl.col("structure").alias("qid"),
        pl.col("structure_right").alias("P3364"),
    ])
    .unique()
    .with_columns(
        qid=pl.col("qid").str.replace("http://www.wikidata.org/entity/(.)", "${1}"),
        P3364=pl.col("P3364").str.replace("http://www.wikidata.org/entity/(.)", "${1}"),
    )
    .select([
        pl.col("qid"),
        pl.col("P3364").alias("-P3364"),
    ])
    .unique()
)

# Export
print(add_df)
print(remove_df)
add_df.write_csv("scratch/maintenance_pair_stereoisomers_add.csv")
remove_df.write_csv("scratch/maintenance_pair_stereoisomers_remove.csv")
