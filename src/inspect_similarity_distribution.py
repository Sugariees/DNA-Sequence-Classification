# ============================================================
# PAIRWISE SIMILARITY DISTRIBUTION
# ============================================================

import pandas as pd
from pathlib import Path

# Load within-species pairwise similarity results.
input_file = Path(
    "results/final/reports/pairwise_sequence_similarity.csv"
)

df = pd.read_csv(input_file)

# Define identity ranges for comparison.
bins = [0, 90, 95, 97, 98, 99, 99.5, 100.0001]

labels = [
    "<=90%",
    "90-95%",
    "95-97%",
    "97-98%",
    "98-99%",
    "99-99.5%",
    "99.5-100%"
]

# Assign each pair to an identity range.
df["identity_range"] = pd.cut(
    df["identity_percent"],
    bins=bins,
    labels=labels,
    include_lowest=True
)

# Display distribution for each species.
print("\nPAIRWISE IDENTITY DISTRIBUTION")
print("=" * 65)

for species in sorted(df["species"].unique()):

    species_df = df[df["species"] == species]

    print(f"\n{species}")

    print(
        species_df["identity_range"]
        .value_counts(sort=False)
        .to_string()
    )

    print("\nIdentity statistics:")

    print(
        species_df["identity_percent"]
        .describe()
        .round(2)
        .to_string()
    )

# Save the distribution.
summary = (
    df.groupby(
        ["species", "identity_range"],
        observed=False
    )
    .size()
    .reset_index(name="pair_count")
)

output_file = Path(
    "results/final/reports/similarity_distribution.csv"
)

summary.to_csv(output_file, index=False)

print("\nSaved report to:", output_file)
