# ============================================================
# K-MER FEATURE EXTRACTION
# PROJECT: BACTERIAL 16S rRNA CLASSIFICATION
#
# PURPOSE:
# Convert DNA sequences into numerical features using
# k-mer frequencies.
#
# INPUT:
# 150 cleaned bacterial 16S rRNA sequences.
#
# OUTPUT:
# 1. Numerical feature matrix (CSV).
# 2. Species labels (CSV).
# ============================================================


# STEP 1: IMPORT LIBRARIES

# SeqIO reads biological sequence files in FASTA format.
from Bio import SeqIO

# Counter counts how frequently each k-mer occurs.
from collections import Counter

# Path manages file and folder locations.
from pathlib import Path

# csv saves numerical features into CSV files.
import csv


# STEP 2: DEFINE PROJECT LOCATIONS

# Folder containing our cleaned FASTA files.
input_folder = Path("data/processed/final")

# Folder where numerical features will be stored.
output_folder = Path("data/processed/final")

output_folder.mkdir(parents=True, exist_ok=True)


# STEP 3: DEFINE SPECIES AND THEIR FASTA FILES

species_files = {
    "Escherichia coli": "Escherichia_coli_16S_clean.fasta",
    "Bacillus subtilis": "Bacillus_subtilis_16S_clean.fasta",
    "Staphylococcus aureus": "Staphylococcus_aureus_16S_clean.fasta"
}


# STEP 4: DEFINE K-MER SIZE

# k = 4 means we examine groups of four consecutive nucleotides.
K = 4


# STEP 5: GENERATE ALL POSSIBLE K-MERS

# DNA contains four possible nucleotides:
# A = Adenine
# T = Thymine
# G = Guanine
# C = Cytosine

# For k = 4, there are 4^4 = 256 possible k-mers.

nucleotides = ["A", "T", "G", "C"]

all_kmers = []

for base_1 in nucleotides:
    for base_2 in nucleotides:
        for base_3 in nucleotides:
            for base_4 in nucleotides:

                kmer = base_1 + base_2 + base_3 + base_4

                all_kmers.append(kmer)


print("K-mer size:", K)
print("Total possible k-mers:", len(all_kmers))


# STEP 6: CREATE STORAGE FOR FEATURES AND LABELS

feature_rows = []

label_rows = []


# STEP 7: PROCESS EACH SPECIES

for species_name, filename in species_files.items():

    print("\nProcessing:", species_name)

    input_file = input_folder / filename

    if not input_file.exists():
        raise FileNotFoundError(
            f"FASTA file not found: {input_file}"
        )

    # Read all sequences from the FASTA file.
    sequences = SeqIO.parse(input_file, "fasta")

    sequence_count = 0

    for record in sequences:

        sequence_count += 1

        # Convert sequence to uppercase text.
        dna_sequence = str(record.seq).upper()

        # Generate overlapping k-mers.
        sequence_kmers = [
            dna_sequence[i:i + K]
            for i in range(len(dna_sequence) - K + 1)
        ]

        # Count how frequently each k-mer occurs.
        kmer_counts = Counter(sequence_kmers)

        # Calculate total number of k-mers.
        total_kmers = len(sequence_kmers)

        # Calculate normalized frequency for every possible k-mer.
        features = {}

        for kmer in all_kmers:

            frequency = kmer_counts[kmer] / total_kmers

            features[kmer] = frequency

        # Store accession and numerical features.
        feature_row = {
            "accession": record.id,
            **features
        }

        feature_rows.append(feature_row)

        # Store the correct species label.
        label_rows.append({
            "accession": record.id,
            "species": species_name
        })

    print("Sequences processed:", sequence_count)


# STEP 8: SAVE FEATURE MATRIX

feature_file = output_folder / "kmer_features.csv"

feature_fieldnames = ["accession"] + all_kmers

with open(
    feature_file,
    "w",
    newline="",
    encoding="utf-8"
) as csv_file:

    writer = csv.DictWriter(
        csv_file,
        fieldnames=feature_fieldnames
    )

    writer.writeheader()

    writer.writerows(feature_rows)


# STEP 9: SAVE SPECIES LABELS

label_file = output_folder / "sequence_labels.csv"

with open(
    label_file,
    "w",
    newline="",
    encoding="utf-8"
) as csv_file:

    writer = csv.DictWriter(
        csv_file,
        fieldnames=["accession", "species"]
    )

    writer.writeheader()

    writer.writerows(label_rows)


# STEP 10: DISPLAY FINAL SUMMARY

print("\n" + "=" * 60)
print("K-MER FEATURE EXTRACTION COMPLETED")
print("=" * 60)

print("Total sequences processed:", len(feature_rows))

print("K-mer size:", K)

print("Features per sequence:", len(all_kmers))

print("Feature matrix:", feature_file)

print("Species labels:", label_file)

print("Expected feature matrix dimensions: 150 x 256")

print("\nAll features successfully extracted!")
