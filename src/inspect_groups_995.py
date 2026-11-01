# ============================================================
# PURPOSE: Inspect sequence similarity groups at 99.5%
# METHOD: Connected components using Union-Find
# ============================================================

import pandas as pd

from pathlib import Path
from collections import defaultdict


# SECTION 1: File paths and thresholds

input_file = Path(
    "results/final/reports/pairwise_sequence_similarity.csv"
)

output_file = Path(
    "results/final/reports/groups_995_summary.csv"
)

IDENTITY_THRESHOLD = 99.5
COVERAGE_THRESHOLD = 90.0


# SECTION 2: Load similarity data

data = pd.read_csv(input_file)


# SECTION 3: Union-Find functions

parent = {}


def find(sequence):
    """Find the representative of a sequence group."""

    if parent[sequence] != sequence:
        parent[sequence] = find(parent[sequence])

    return parent[sequence]


def union(sequence_1, sequence_2):
    """Join two sequences into the same connected group."""

    root_1 = find(sequence_1)
    root_2 = find(sequence_2)

    if root_1 != root_2:
        parent[root_2] = root_1


# SECTION 4: Process each species independently

all_groups = []

for species in sorted(data["species"].unique()):

    species_data = data[data["species"] == species]

    # Collect every accession belonging to this species.
    accessions = set(species_data["accession_1"])
    accessions.update(species_data["accession_2"])

    parent = {accession: accession for accession in accessions}

    # Join sequences that meet both similarity and coverage thresholds.
    for _, row in species_data.iterrows():

        if (
            row["identity_percent"] >= IDENTITY_THRESHOLD
            and row["coverage_1_percent"] >= COVERAGE_THRESHOLD
            and row["coverage_2_percent"] >= COVERAGE_THRESHOLD
        ):
            union(row["accession_1"], row["accession_2"])

    # Collect accessions into connected components.
    groups = defaultdict(list)

    for accession in sorted(accessions):
        groups[find(accession)].append(accession)

    # Sort groups from largest to smallest.
    sorted_groups = sorted(
        groups.values(),
        key=len,
        reverse=True
    )

    print("\nSPECIES:", species)
    print("=" * 60)
    print("Total sequences:", len(accessions))
    print("Number of groups:", len(sorted_groups))

    for number, members in enumerate(sorted_groups, start=1):

        group_id = f"G{number:03d}"

        print(
            f"{group_id}: {len(members)} sequences"
        )

        all_groups.append({
            "species": species,
            "group_id": group_id,
            "group_size": len(members),
            "accessions": ";".join(members)
        })


# SECTION 5: Save group information

group_table = pd.DataFrame(all_groups)

group_table.to_csv(output_file, index=False)

print("\nGROUP SIZE SUMMARY")
print("=" * 60)

print(
    group_table.groupby("species")["group_size"]
    .agg(["count", "min", "max", "mean"])
)

print("\nReport saved:", output_file)
print("\nANALYSIS COMPLETED!")
