# ============================================================
# FINAL SEQUENCE SIMILARITY SUMMARY
# PROJECT: BACTERIAL 16S rRNA CLASSIFICATION
#
# PURPOSE:
# Summarize pairwise sequence identity values for each species.
# Calculate minimum, maximum, mean and median identity.
# Count pairs above selected identity thresholds.
# Save the summary as a CSV report.
# ============================================================


# STEP 1: IMPORT REQUIRED LIBRARIES

# csv reads and writes CSV files.
import csv

# statistics calculates mean and median.
import statistics

# defaultdict organizes results separately by species.
from collections import defaultdict

# Path manages file and folder locations.
from pathlib import Path


# STEP 2: DEFINE INPUT AND OUTPUT FILES

# Existing pairwise similarity report.
input_file = Path(
    "results/final/reports/pairwise_sequence_similarity.csv"
)

# Output summary report.
output_file = Path(
    "results/final/reports/sequence_similarity_summary.csv"
)


# STEP 3: CHECK WHETHER THE INPUT FILE EXISTS

if not input_file.exists():
    print("ERROR: Pairwise similarity report not found.")
    print("Expected location:", input_file)
    raise SystemExit


# STEP 4: CREATE A DICTIONARY TO STORE IDENTITY VALUES

# Each species will have its own list of identity percentages.
species_identities = defaultdict(list)


# STEP 5: READ THE PAIRWISE SIMILARITY REPORT

with open(input_file, "r", newline="", encoding="utf-8") as csv_file:

    reader = csv.DictReader(csv_file)

    for row in reader:

        # Extract the species name.
        species_name = row["species"]

        # Convert identity percentage from text to a number.
        identity = float(row["identity_percent"])

        # Store the identity value under the corresponding species.
        species_identities[species_name].append(identity)


# STEP 6: DEFINE IDENTITY THRESHOLDS

# We will count how many sequence pairs exceed each threshold.
thresholds = [97, 98, 99, 99.5]


# STEP 7: CALCULATE SUMMARY STATISTICS

summary_results = []

for species_name, identities in species_identities.items():

    print("\n" + "=" * 60)
    print("SPECIES:", species_name)
    print("=" * 60)

    # Calculate basic statistics.
    minimum_identity = min(identities)
    maximum_identity = max(identities)
    mean_identity = statistics.mean(identities)
    median_identity = statistics.median(identities)

    print("Total pairwise comparisons:", len(identities))
    print(f"Minimum identity: {minimum_identity:.2f}%")
    print(f"Maximum identity: {maximum_identity:.2f}%")
    print(f"Mean identity: {mean_identity:.2f}%")
    print(f"Median identity: {median_identity:.2f}%")

    # Store the results in a dictionary.
    summary = {
        "species": species_name,
        "total_pairs": len(identities),
        "minimum_identity": round(minimum_identity, 2),
        "maximum_identity": round(maximum_identity, 2),
        "mean_identity": round(mean_identity, 2),
        "median_identity": round(median_identity, 2)
    }

    # Count pairs at or above each identity threshold.
    for threshold in thresholds:

        count = sum(
            identity >= threshold for identity in identities
        )

        percentage = (count / len(identities)) * 100

        print(
            f"Pairs with at least {threshold}% identity:",
            count,
            f"({percentage:.2f}%)"
        )

        # Add threshold statistics to the summary.
        summary[f"pairs_at_least_{threshold}_percent"] = count
        summary[f"percentage_at_least_{threshold}_percent"] = round(
            percentage, 2
        )

    summary_results.append(summary)


# STEP 8: SAVE SUMMARY RESULTS TO CSV

fieldnames = list(summary_results[0].keys())

with open(output_file, "w", newline="", encoding="utf-8") as csv_file:

    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

    writer.writeheader()
    writer.writerows(summary_results)


# STEP 9: DISPLAY COMPLETION MESSAGE

print("\n" + "=" * 60)
print("SIMILARITY SUMMARY COMPLETED")
print("=" * 60)

print("Summary saved to:", output_file)
