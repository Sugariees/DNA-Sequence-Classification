# ============================================================
# INSPECT SEQUENCE GROUPS
# ============================================================

import pandas as pd
from pathlib import Path

# Load the sequence group assignments.
input_file = Path("data/processed/final/sequence_groups.csv")

groups_df = pd.read_csv(input_file)

# Display the number of sequences in each group.
print("\nGROUP SIZES")
print("=" * 60)

group_sizes = (
    groups_df.groupby(["species", "group_id"])
    .size()
    .reset_index(name="sequence_count")
    .sort_values(["species", "sequence_count"], ascending=[True, False])
)

for species in groups_df["species"].unique():

    print(f"\n{species}")

    species_groups = group_sizes[
        group_sizes["species"] == species
    ]

    print(species_groups.to_string(index=False))


# Display all accessions belonging to each group.
print("\n\nACCESSIONS IN EACH GROUP")
print("=" * 60)

for _, row in group_sizes.iterrows():

    group_id = row["group_id"]

    members = groups_df[
        groups_df["group_id"] == group_id
    ]["accession"].tolist()

    print(f"\n{group_id} ({len(members)} sequences)")

    for accession in members:
        print("  ", accession)


# Save detailed group sizes for later analysis.
output_file = Path(
    "results/final/reports/sequence_group_sizes.csv"
)

output_file.parent.mkdir(parents=True, exist_ok=True)

group_sizes.to_csv(output_file, index=False)

print("\nGroup size report saved to:", output_file)
