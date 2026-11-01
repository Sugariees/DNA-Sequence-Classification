# ============================================================
# SEQUENCE SIMILARITY GROUPING
# PROJECT: BACTERIAL 16S rRNA CLASSIFICATION
#
# PURPOSE:
# 1. Load within-species pairwise sequence similarities.
# 2. Group highly similar sequences together.
# 3. Count the groups present in each species.
# 4. Save the group assignments for leakage-aware ML splitting.
#
# IMPORTANT:
# Sequences in the same group must stay together when
# creating training and testing datasets.
# ============================================================


# STEP 1: IMPORT LIBRARIES

import pandas as pd
from pathlib import Path


# STEP 2: DEFINE FILE LOCATIONS

input_file = Path(
    "results/final/reports/pairwise_sequence_similarity.csv"
)

output_folder = Path("data/processed/final")

output_file = output_folder / "sequence_groups.csv"


# STEP 3: SET SIMILARITY THRESHOLD

# Sequences with at least 99% identity and sufficient
# alignment coverage will be connected into the same group.

SIMILARITY_THRESHOLD = 99.0
MIN_COVERAGE = 90.0


# STEP 4: LOAD PAIRWISE SIMILARITY DATA

print("\nLoading pairwise sequence similarity data...")

similarity_df = pd.read_csv(input_file)

print("Total pairwise comparisons:", len(similarity_df))

print("\nAvailable columns:")

print(similarity_df.columns.tolist())


# STEP 5: IDENTIFY REQUIRED COLUMNS

required_columns = [
    "species",
    "accession_1",
    "accession_2",
    "identity_percent",
    "coverage_1_percent",
    "coverage_2_percent"
]

missing_columns = [
    column
    for column in required_columns
    if column not in similarity_df.columns
]

if missing_columns:
    raise ValueError(
        f"Missing required columns: {missing_columns}"
    )


# STEP 6: CREATE SEQUENCE GROUPS

all_group_records = []

species_list = similarity_df["species"].unique()

for species in species_list:

    print("\n" + "=" * 60)

    print("Processing:", species)

    species_df = similarity_df[
        similarity_df["species"] == species
    ]

    # Collect every accession present in this species.
    accessions = sorted(
        set(species_df["accession_1"])
        | set(species_df["accession_2"])
    )

    # Each sequence initially belongs to its own group.
    parent = {
        accession: accession
        for accession in accessions
    }

    # Find the root parent of a sequence.
    def find_root(sequence):

        if parent[sequence] != sequence:

            parent[sequence] = find_root(parent[sequence])

        return parent[sequence]

    # Connect two sequences into the same group.
    def connect_sequences(sequence_1, sequence_2):

        root_1 = find_root(sequence_1)

        root_2 = find_root(sequence_2)

        if root_1 != root_2:

            parent[root_2] = root_1

    # STEP 6A: CONNECT HIGHLY SIMILAR SEQUENCES

    qualifying_pairs = species_df[
        (species_df["identity_percent"] >= SIMILARITY_THRESHOLD)
        & (species_df["coverage_1_percent"] >= MIN_COVERAGE)
        & (species_df["coverage_2_percent"] >= MIN_COVERAGE)
    ]

    for _, row in qualifying_pairs.iterrows():

        connect_sequences(
            row["accession_1"],
            row["accession_2"]
        )

    # STEP 6B: ASSIGN FINAL GROUP NUMBERS

    root_to_group = {}

    for accession in accessions:

        root = find_root(accession)

        if root not in root_to_group:

            root_to_group[root] = (
                f"{species.replace(' ', '_')}_G"
                f"{len(root_to_group) + 1:03d}"
            )

        group_id = root_to_group[root]

        all_group_records.append({
            "accession": accession,
            "species": species,
            "group_id": group_id
        })

    # STEP 6C: DISPLAY GROUP INFORMATION

    group_sizes = pd.DataFrame([
        record
        for record in all_group_records
        if record["species"] == species
    ])["group_id"].value_counts()

    print("Sequences:", len(accessions))

    print("Qualifying similar pairs:", len(qualifying_pairs))

    print("Distinct groups:", len(group_sizes))

    print("Largest group size:", group_sizes.max())

    print("Smallest group size:", group_sizes.min())


# STEP 7: CREATE AND SAVE GROUP ASSIGNMENTS

groups_df = pd.DataFrame(all_group_records)

output_folder.mkdir(parents=True, exist_ok=True)

groups_df.to_csv(output_file, index=False)


# STEP 8: DISPLAY FINAL SUMMARY

print("\n" + "=" * 60)

print("SEQUENCE GROUPING SUMMARY")

print("=" * 60)

summary = (
    groups_df.groupby("species")
    .agg(
        total_sequences=("accession", "count"),
        distinct_groups=("group_id", "nunique")
    )
)

print(summary)

print("\nTotal sequences:", len(groups_df))

print("Total distinct groups:", groups_df["group_id"].nunique())

print("\nSaved group assignments to:", output_file)

print("\nSEQUENCE GROUPING COMPLETED!")
