# ============================================================
# FINAL SEQUENCE SIMILARITY ANALYSIS
# PROJECT: BACTERIAL 16S rRNA CLASSIFICATION
#
# PURPOSE:
# Compare every possible pair of sequences within each species.
# Calculate sequence identity and coverage.
# Save all comparison results in a CSV file.
# ============================================================

# STEP 1: IMPORT REQUIRED LIBRARIES

# AlignIO reads the aligned FASTA files produced by MAFFT.
from Bio import AlignIO

# combinations generates every possible pair of sequences.
from itertools import combinations

# csv allows us to save results in a spreadsheet-friendly format.
import csv

# Path helps us manage file and folder locations.
from pathlib import Path


# STEP 2: DEFINE INPUT AND OUTPUT LOCATIONS

# Folder containing our three MAFFT alignment files.
alignment_folder = Path("results/final/alignment")

# Folder where we will save the similarity report.
report_folder = Path("results/final/reports")

# Create the report folder if it does not already exist.
report_folder.mkdir(parents=True, exist_ok=True)

# Final CSV report location.
output_file = report_folder / "pairwise_sequence_similarity.csv"


# STEP 3: DEFINE THE THREE SPECIES FILES

species_files = {
    "Escherichia coli": "Escherichia_coli_aligned.fasta",
    "Bacillus subtilis": "Bacillus_subtilis_aligned.fasta",
    "Staphylococcus aureus": "Staphylococcus_aureus_aligned.fasta"
}


# STEP 4: CREATE A LIST TO STORE ALL COMPARISON RESULTS

all_results = []


# STEP 5: ANALYSE EACH SPECIES

for species_name, filename in species_files.items():

    print("\n" + "=" * 65)
    print("SPECIES:", species_name)
    print("=" * 65)

    # Construct the complete path to the alignment file.
    alignment_file = alignment_folder / filename

    # Read the aligned FASTA file.
    alignment = AlignIO.read(alignment_file, "fasta")

    # Display the number of sequences.
    print("Number of sequences:", len(alignment))

    # Calculate the total number of alignment positions.
    print("Alignment length:", alignment.get_alignment_length())

    # STEP 6: COMPARE EVERY POSSIBLE SEQUENCE PAIR

    for sequence_1, sequence_2 in combinations(alignment, 2):

        # Count identical nucleotide positions.
        matches = 0

        # Count positions where both sequences have a nucleotide.
        compared_positions = 0

        # Original sequence lengths, excluding alignment gaps.
        length_1 = len(sequence_1.seq) - sequence_1.seq.count("-")
        length_2 = len(sequence_2.seq) - sequence_2.seq.count("-")

        # Compare the two aligned sequences position by position.
        for base_1, base_2 in zip(sequence_1.seq, sequence_2.seq):

            # Ignore positions containing a gap in either sequence.
            if base_1 == "-" or base_2 == "-":
                continue

            compared_positions += 1

            # Check whether both nucleotides are identical.
            if base_1 == base_2:
                matches += 1

        # Calculate percentage sequence identity.
        if compared_positions > 0:
            identity = (matches / compared_positions) * 100
        else:
            identity = 0

        # Calculate how much of each sequence was compared.
        coverage_1 = (compared_positions / length_1) * 100
        coverage_2 = (compared_positions / length_2) * 100

        # Store the comparison results.
        all_results.append({
            "species": species_name,
            "accession_1": sequence_1.id,
            "accession_2": sequence_2.id,
            "identity_percent": round(identity, 2),
            "compared_positions": compared_positions,
            "coverage_1_percent": round(coverage_1, 2),
            "coverage_2_percent": round(coverage_2, 2)
        })


# STEP 7: SAVE ALL RESULTS INTO A CSV FILE

fieldnames = [
    "species",
    "accession_1",
    "accession_2",
    "identity_percent",
    "compared_positions",
    "coverage_1_percent",
    "coverage_2_percent"
]

with open(output_file, "w", newline="", encoding="utf-8") as csv_file:

    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

    # Write column headings.
    writer.writeheader()

    # Write all sequence comparison results.
    writer.writerows(all_results)


# STEP 8: DISPLAY A SUMMARY

print("\n" + "=" * 65)
print("SIMILARITY ANALYSIS COMPLETED")
print("=" * 65)

print("Total sequence pairs compared:", len(all_results))

print("Results saved to:", output_file)

# Display the number of pairs with identity of at least 99%.
high_similarity_pairs = [
    result for result in all_results
    if result["identity_percent"] >= 99
]

print(
    "Pairs with at least 99% identity:",
    len(high_similarity_pairs)
)

print("\nFirst 10 high-similarity pairs:")

for result in high_similarity_pairs[:10]:

    print(
        result["species"],
        "|",
        result["accession_1"],
        "vs",
        result["accession_2"],
        "| Identity:",
        result["identity_percent"],
        "%"
    )
