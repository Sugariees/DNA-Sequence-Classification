# ============================================================
# EXACT DUPLICATE SEQUENCE ANALYSIS
#
# PURPOSE:
# Identify identical DNA sequences using SHA256 fingerprints.
# Compare duplicate records within and between species.
# ============================================================

import pandas as pd
from pathlib import Path


# STEP 1: LOAD SEQUENCE METADATA

metadata_file = Path(
    "data/processed/final/sequence_metadata.csv"
)

df = pd.read_csv(metadata_file)


# STEP 2: CHECK REQUIRED COLUMNS

required_columns = [
    "species",
    "accession",
    "sequence_length",
    "sha256"
]

for column in required_columns:

    if column not in df.columns:
        raise ValueError(f"Missing column: {column}")


# STEP 3: FIND IDENTICAL SEQUENCE GROUPS

duplicate_groups = (
    df.groupby("sha256")
    .agg(
        sequence_count=("accession", "count"),
        species_count=("species", "nunique"),
        species=("species", lambda x: ", ".join(sorted(set(x)))),
        accessions=("accession", lambda x: ", ".join(x)),
        sequence_length=("sequence_length", "first")
    )
    .reset_index()
)

duplicate_groups = duplicate_groups[
    duplicate_groups["sequence_count"] > 1
]


# STEP 4: DISPLAY EXACT DUPLICATES

print("\nEXACT DUPLICATE SEQUENCE ANALYSIS")
print("=" * 65)

print("Total records:", len(df))

print("Unique sequence fingerprints:", df["sha256"].nunique())

print("Number of duplicate sequence groups:", len(duplicate_groups))

if len(duplicate_groups) > 0:

    print("\nDuplicate sequence groups:")

    print(
        duplicate_groups[
            [
                "sequence_count",
                "species_count",
                "species",
                "accessions",
                "sequence_length"
            ]
        ].to_string(index=False)
    )

else:

    print("\nNo exact duplicate sequences found.")


# STEP 5: CHECK CROSS-SPECIES IDENTICAL SEQUENCES

cross_species_duplicates = duplicate_groups[
    duplicate_groups["species_count"] > 1
]

print("\nCross-species exact duplicate groups:",
      len(cross_species_duplicates))


# STEP 6: SAVE RESULTS

output_file = Path(
    "results/final/reports/exact_duplicate_analysis.csv"
)

output_file.parent.mkdir(parents=True, exist_ok=True)

duplicate_groups.to_csv(output_file, index=False)

print("\nReport saved to:", output_file)

print("\nANALYSIS COMPLETED!")
