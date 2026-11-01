# ============================================================
# INSPECT HIGHLY SIMILAR SEQUENCE PAIRS
# ============================================================

import pandas as pd
from pathlib import Path

# Load pairwise similarity results.
input_file = Path(
    "results/final/reports/pairwise_sequence_similarity.csv"
)

df = pd.read_csv(input_file)

# Define identity thresholds to inspect.
thresholds = [100.0, 99.9, 99.5, 99.0]

print("\nHIGH-SIMILARITY PAIR ANALYSIS")
print("=" * 70)

for species in sorted(df["species"].unique()):

    species_df = df[df["species"] == species]

    print(f"\n{species}")

    for threshold in thresholds:

        count = (
            species_df["identity_percent"] >= threshold
        ).sum()

        print(
            f"Pairs with identity >= {threshold}%: {count}"
        )

# Extract pairs with identity >= 99.5%.
high_similarity = df[
    df["identity_percent"] >= 99.5
].copy()

# Sort from highest to lowest identity.
high_similarity = high_similarity.sort_values(
    ["species", "identity_percent"],
    ascending=[True, False]
)

# Save the detailed results.
output_file = Path(
    "results/final/reports/high_similarity_pairs.csv"
)

high_similarity.to_csv(output_file, index=False)

print("\nTotal pairs with identity >= 99.5%:",
      len(high_similarity))

print("Detailed results saved to:", output_file)

print("\nTop 20 highly similar pairs:")

print(
    high_similarity.head(20).to_string(index=False)
)
