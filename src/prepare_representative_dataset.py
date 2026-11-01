# ============================================================
# PROJECT: DNA Sequence Classification
# SCRIPT: Prepare Representative ML Dataset
#
# PURPOSE:
# Extract normalized 4-mer frequencies from the 37
# CD-HIT representative DNA sequences.
#
# OUTPUT:
# data/processed/diversity_controlled/
#     representative_ml_dataset.csv
# ============================================================


# -------------------- IMPORT LIBRARIES ----------------------

from Bio import SeqIO
from collections import Counter
from itertools import product
from pathlib import Path

import pandas as pd


# -------------------- PROJECT PATHS --------------------------

PROJECT_DIR = Path(__file__).resolve().parent.parent

INPUT_DIR = PROJECT_DIR / "data" / "processed" / "diversity_controlled"

OUTPUT_FILE = INPUT_DIR / "representative_ml_dataset.csv"


# -------------------- CONFIGURATION --------------------------

K = 4

SPECIES_FILES = {
    "Escherichia coli": "Escherichia_coli_clustered.fasta",
    "Bacillus subtilis": "Bacillus_subtilis_clustered.fasta",
    "Staphylococcus aureus": "Staphylococcus_aureus_clustered.fasta"
}


# -------------------- GENERATE ALL K-MERS --------------------

# Generate every possible combination of A, T, G and C
# for a k-mer length of 4.

BASES = "ACGT"

ALL_KMERS = [
    "".join(kmer)
    for kmer in product(BASES, repeat=K)
]


# -------------------- EXTRACT K-MER FEATURES -----------------

def extract_kmer_features(sequence):

    # Convert sequence into uppercase format
    sequence = str(sequence).upper()

    # Ensure sequence contains only valid DNA bases
    if not sequence or any(base not in BASES for base in sequence):
        raise ValueError("Invalid DNA sequence encountered.")

    # Count overlapping k-mers
    kmer_counts = Counter(
        sequence[i:i + K]
        for i in range(len(sequence) - K + 1)
    )

    # Calculate total number of k-mer windows
    total_windows = len(sequence) - K + 1

    # Normalize frequencies for all 256 possible k-mers
    features = {
        kmer: kmer_counts.get(kmer, 0) / total_windows
        for kmer in ALL_KMERS
    }

    return features


# -------------------- PROCESS REPRESENTATIVE FASTA FILES ----

all_records = []

for species, filename in SPECIES_FILES.items():

    fasta_path = INPUT_DIR / filename

    # Check whether the representative FASTA file exists
    if not fasta_path.exists():
        raise FileNotFoundError(
            f"Missing input file: {fasta_path}"
        )

    records = list(SeqIO.parse(fasta_path, "fasta"))

    print(f"\nProcessing: {species}")
    print(f"Representative sequences: {len(records)}")

    # Process every representative sequence
    for record in records:

        sequence = str(record.seq).upper()

        # Extract normalized k-mer frequencies
        features = extract_kmer_features(sequence)

        # Store accession, species label and features
        row = {
            "accession": record.id,
            "species": species,
            "sequence_length": len(sequence)
        }

        row.update(features)

        all_records.append(row)


# -------------------- CREATE DATAFRAME -----------------------

dataset = pd.DataFrame(all_records)


# -------------------- VALIDATE DATASET ------------------------

expected_features = 4 ** K

assert len(dataset) == 37, (
    f"Expected 37 sequences, found {len(dataset)}"
)

assert dataset["accession"].is_unique, (
    "Duplicate accessions detected."
)

assert dataset["species"].notna().all()

feature_columns = ALL_KMERS

assert len(feature_columns) == expected_features

# Verify that normalized frequencies sum approximately to 1
frequency_sums = dataset[feature_columns].sum(axis=1)

assert ((frequency_sums - 1.0).abs() < 1e-9).all(), (
    "Some sequences have incorrect normalized frequencies."
)


# -------------------- SAVE DATASET ---------------------------

dataset.to_csv(OUTPUT_FILE, index=False)


# -------------------- DISPLAY RESULTS ------------------------

print("\n==========================================")
print("REPRESENTATIVE DATASET CREATED")
print("==========================================")

print(f"\nTotal sequences: {len(dataset)}")

print(f"Number of k-mer features: {len(feature_columns)}")

print("\nSpecies distribution:")
print(dataset["species"].value_counts())

print("\nDataset dimensions:")
print(dataset.shape)

print(f"\nSaved to:\n{OUTPUT_FILE}")

print("\nAll validation checks passed!")
