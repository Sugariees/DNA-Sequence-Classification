# ============================================================
# COMPARE SEQUENCE GROUPING THRESHOLDS
#
# PURPOSE:
# Examine how the number and size of sequence groups change
# at different pairwise identity thresholds.
# ============================================================

import pandas as pd
from pathlib import Path


# STEP 1: LOAD PAIRWISE SIMILARITY RESULTS

input_file = Path(
    "results/final/reports/pairwise_sequence_similarity.csv"
)

similarity_df = pd.read_csv(input_file)


# STEP 2: DEFINE THRESHOLDS

thresholds = [99.5, 99.0, 98.0, 97.0, 95.0]

MIN_COVERAGE = 90.0


# STEP 3: GROUPING FUNCTION

def count_sequence_groups(species_df, threshold):

    accessions = sorted(
        set(species_df["accession_1"])
        | set(species_df["accession_2"])
    )

    parent = {accession: accession for accession in accessions}

    def find_root(sequence):

        if parent[sequence] != sequence:
            parent[sequence] = find_root(parent[sequence])

        return parent[sequence]

    def connect(sequence_1, sequence_2):

        root_1 = find_root(sequence_1)
        root_2 = find_root(sequence_2)

        if root_1 != root_2:
            parent[root_2] = root_1

    qualifying_pairs = species_df[
        (species_df["identity_percent"] >= threshold)
        & (species_df["coverage_1_percent"] >= MIN_COVERAGE)
        & (species_df["coverage_2_percent"] >= MIN_COVERAGE)
    ]

    for _, row in qualifying_pairs.iterrows():

        connect(row["accession_1"], row["accession_2"])

    group_sizes = {}

    for accession in accessions:

        root = find_root(accession)

        group_sizes[root] = group_sizes.get(root, 0) + 1

    return {
        "total_sequences": len(accessions),
        "number_of_groups": len(group_sizes),
        "largest_group": max(group_sizes.values()),
        "smallest_group": min(group_sizes.values()),
        "qualifying_pairs": len(qualifying_pairs)
    }


# STEP 4: COMPARE THRESHOLDS FOR EACH SPECIES

results = []

for species in sorted(similarity_df["species"].unique()):

    species_df = similarity_df[
        similarity_df["species"] == species
    ]

    for threshold in thresholds:

        summary = count_sequence_groups(
            species_df,
            threshold
        )

        results.append({
            "species": species,
            "identity_threshold": threshold,
            **summary
        })


# STEP 5: DISPLAY RESULTS

results_df = pd.DataFrame(results)

print("\nSEQUENCE GROUPING THRESHOLD COMPARISON")
print("=" * 80)

print(results_df.to_string(index=False))


# STEP 6: SAVE RESULTS

output_file = Path(
    "results/final/reports/group_threshold_comparison.csv"
)

output_file.parent.mkdir(parents=True, exist_ok=True)

results_df.to_csv(output_file, index=False)

print("\nResults saved to:", output_file)
