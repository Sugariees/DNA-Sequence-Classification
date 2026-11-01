# ============================================================
# PROJECT: DNA Sequence Classification
# PURPOSE: Investigate similarity within large sequence groups
# METHOD: Medoid representative + pairwise identity analysis
# ============================================================

import pandas as pd

from pathlib import Path


# SECTION 1: Define file paths

similarity_file = Path(
    "results/final/reports/pairwise_sequence_similarity.csv"
)

groups_file = Path(
    "results/final/reports/groups_995_summary.csv"
)

output_directory = Path("results/final/reports")

output_directory.mkdir(parents=True, exist_ok=True)


# SECTION 2: Load the datasets

similarity_data = pd.read_csv(similarity_file)

groups_data = pd.read_csv(groups_file)


# SECTION 3: Store results

representative_results = []
group_results = []


# SECTION 4: Analyse each sequence group

for _, group_row in groups_data.iterrows():

    species = group_row["species"]
    group_id = group_row["group_id"]

    accessions = group_row["accessions"].split(";")

    # Ignore singleton groups because they have no pairwise comparisons.
    if len(accessions) < 2:
        continue

    # Select pairwise comparisons belonging to this species.
    species_pairs = similarity_data[
        similarity_data["species"] == species
    ]

    # Keep only comparisons where both sequences belong to this group.
    group_pairs = species_pairs[
        species_pairs["accession_1"].isin(accessions)
        & species_pairs["accession_2"].isin(accessions)
    ].copy()

    if group_pairs.empty:
        continue

    # SECTION 5: Calculate average identity for each sequence.

    identity_scores = {}

    for accession in accessions:

        relevant_pairs = group_pairs[
            (group_pairs["accession_1"] == accession)
            | (group_pairs["accession_2"] == accession)
        ]

        identity_scores[accession] = (
            relevant_pairs["identity_percent"].mean()
        )

    # The medoid is the sequence with the highest average identity.
    representative = max(
        identity_scores,
        key=identity_scores.get
    )

    representative_average = identity_scores[representative]

    # SECTION 6: Compare every sequence with the representative.

    representative_pairs = []

    for accession in accessions:

        if accession == representative:
            continue

        matching_pair = group_pairs[
            (
                (group_pairs["accession_1"] == representative)
                & (group_pairs["accession_2"] == accession)
            )
            |
            (
                (group_pairs["accession_1"] == accession)
                & (group_pairs["accession_2"] == representative)
            )
        ]

        if matching_pair.empty:
            continue

        row = matching_pair.iloc[0]

        representative_pairs.append(
            row["identity_percent"]
        )

        representative_results.append({
            "species": species,
            "group_id": group_id,
            "representative": representative,
            "accession": accession,
            "identity_percent": row["identity_percent"]
        })

    # SECTION 7: Summarise all pairwise comparisons within the group.

    group_results.append({
        "species": species,
        "group_id": group_id,
        "group_size": len(accessions),
        "representative": representative,
        "representative_average_identity": round(
            representative_average, 3
        ),
        "representative_min_identity": round(
            min(representative_pairs), 3
        ),
        "representative_max_identity": round(
            max(representative_pairs), 3
        ),
        "within_group_pair_count": len(group_pairs),
        "within_group_min_identity": round(
            group_pairs["identity_percent"].min(), 3
        ),
        "within_group_mean_identity": round(
            group_pairs["identity_percent"].mean(), 3
        ),
        "within_group_max_identity": round(
            group_pairs["identity_percent"].max(), 3
        )
    })


# SECTION 8: Create result tables

representative_table = pd.DataFrame(representative_results)

group_summary = pd.DataFrame(group_results)


# SECTION 9: Display large-group results

print("\nLARGE SEQUENCE GROUP ANALYSIS")
print("=" * 75)

if not group_summary.empty:

    for _, row in group_summary.iterrows():

        print(f"\nSpecies: {row['species']}")
        print(f"Group: {row['group_id']}")
        print(f"Number of sequences: {row['group_size']}")

        print(
            "Representative:",
            row["representative"]
        )

        print(
            "Representative average identity:",
            row["representative_average_identity"], "%"
        )

        print(
            "Representative minimum identity:",
            row["representative_min_identity"], "%"
        )

        print(
            "Within-group mean identity:",
            row["within_group_mean_identity"], "%"
        )

        print(
            "Within-group minimum identity:",
            row["within_group_min_identity"], "%"
        )

        print(
            "Within-group maximum identity:",
            row["within_group_max_identity"], "%"
        )


# SECTION 10: Save reports

group_summary.to_csv(
    output_directory / "large_group_similarity_summary.csv",
    index=False
)

representative_table.to_csv(
    output_directory / "representative_similarity_details.csv",
    index=False
)

print("\nREPORTS SAVED")
print("=" * 75)

print("1. large_group_similarity_summary.csv")
print("2. representative_similarity_details.csv")

print("\nANALYSIS COMPLETED!")
