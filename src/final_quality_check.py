# ============================================================
# FINAL DATASET QUALITY CONTROL
# PROJECT: BACTERIAL 16S rRNA CLASSIFICATION
#
# PURPOSE:
# Check sequence lengths, nucleotide composition,
# exact duplicates, and metadata consistency.
# ============================================================


# ------------------------------------------------------------
# STEP 1: IMPORT LIBRARIES
# ------------------------------------------------------------

from Bio import SeqIO
from pathlib import Path
import csv


# ------------------------------------------------------------
# STEP 2: DEFINE INPUT FILES AND QUALITY LIMITS
# ------------------------------------------------------------

DATA_FOLDER = Path("data/processed/final")

SPECIES_FILES = {
    "Escherichia coli":
        DATA_FOLDER / "Escherichia_coli_16S_clean.fasta",

    "Bacillus subtilis":
        DATA_FOLDER / "Bacillus_subtilis_16S_clean.fasta",

    "Staphylococcus aureus":
        DATA_FOLDER / "Staphylococcus_aureus_16S_clean.fasta"
}

METADATA_FILE = DATA_FOLDER / "sequence_metadata.csv"

MIN_LENGTH = 1300
MAX_LENGTH = 1600

VALID_BASES = set("ATGC")


# ------------------------------------------------------------
# STEP 3: STORAGE FOR QUALITY CONTROL RESULTS
# ------------------------------------------------------------

total_sequences = 0

# Store sequence fingerprints separately for each species
species_sequences = {}

# Store all fingerprints to identify cross-species duplicates
all_sequence_owners = {}

# Store accession numbers found in FASTA files
fasta_accessions = set()

# Store sequence IDs and their species
fasta_records = []


# ------------------------------------------------------------
# STEP 4: ANALYSE EACH SPECIES
# ------------------------------------------------------------

print("\nFINAL DATASET QUALITY REPORT")
print("=" * 60)

for species, fasta_file in SPECIES_FILES.items():

    print("\nSpecies:", species)
    print("-" * 50)

    # Check whether the FASTA file exists
    if not fasta_file.exists():

        print("ERROR: FASTA file not found:", fasta_file)
        continue

    # Read all sequences from the FASTA file
    records = list(
        SeqIO.parse(fasta_file, "fasta")
    )

    total_sequences += len(records)

    # Store fingerprints for this species
    species_sequences[species] = {}

    # Storage for length statistics
    lengths = []

    # Count sequences containing ambiguous bases
    ambiguous_count = 0

    # Count sequences outside the expected length range
    length_failure_count = 0

    # Count exact duplicates within this species
    duplicate_count = 0

    # Check every sequence
    for record in records:

        sequence = str(record.seq).upper()

        lengths.append(len(sequence))

        # Record accession number
        fasta_accessions.add(record.id)

        fasta_records.append({
            "accession": record.id,
            "species": species
        })

        # Check sequence length
        if not MIN_LENGTH <= len(sequence) <= MAX_LENGTH:

            length_failure_count += 1

        # Check nucleotide composition
        if not set(sequence).issubset(VALID_BASES):

            ambiguous_count += 1

        # Check for exact duplicates within this species
        if sequence in species_sequences[species]:

            duplicate_count += 1

        else:

            species_sequences[species][sequence] = record.id

        # Track which species contain this exact sequence
        all_sequence_owners.setdefault(
            sequence, set()
        ).add(species)

    # Display statistics
    if lengths:

        print("Number of sequences:", len(records))

        print("Shortest sequence:", min(lengths), "bp")

        print("Longest sequence:", max(lengths), "bp")

        print(
            "Average length:",
            round(sum(lengths) / len(lengths), 2),
            "bp"
        )

        print(
            "Sequences outside length range:",
            length_failure_count
        )

        print(
            "Sequences with ambiguous bases:",
            ambiguous_count
        )

        print(
            "Exact duplicates within species:",
            duplicate_count
        )


# ------------------------------------------------------------
# STEP 5: CHECK CROSS-SPECIES EXACT DUPLICATES
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("CROSS-SPECIES DUPLICATE CHECK")
print("=" * 60)

cross_species_duplicates = []

for sequence, owners in all_sequence_owners.items():

    if len(owners) > 1:

        cross_species_duplicates.append({
            "sequence": sequence,
            "species": sorted(owners)
        })

print(
    "Exact sequences shared across species:",
    len(cross_species_duplicates)
)

if cross_species_duplicates:

    for duplicate in cross_species_duplicates:

        print(
            "Shared by:",
            ", ".join(duplicate["species"]),
            "| Length:",
            len(duplicate["sequence"])
        )


# ------------------------------------------------------------
# STEP 6: CHECK METADATA CONSISTENCY
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("METADATA CONSISTENCY CHECK")
print("=" * 60)

if not METADATA_FILE.exists():

    print("ERROR: Metadata file not found.")

else:

    with open(
        METADATA_FILE,
        "r",
        newline="",
        encoding="utf-8"
    ) as file:

        metadata = list(
            csv.DictReader(file)
        )

    metadata_accessions = {
        row["accession"] for row in metadata
    }

    print("FASTA sequence count:", total_sequences)

    print("Metadata record count:", len(metadata))

    # Check whether FASTA accessions exist in metadata
    missing_metadata = (
        fasta_accessions - metadata_accessions
    )

    # Check whether metadata contains unexpected accessions
    extra_metadata = (
        metadata_accessions - fasta_accessions
    )

    print(
        "FASTA accessions missing from metadata:",
        len(missing_metadata)
    )

    print(
        "Metadata accessions absent from FASTA:",
        len(extra_metadata)
    )

    if missing_metadata:

        print("Missing accessions:", missing_metadata)

    if extra_metadata:

        print("Unexpected accessions:", extra_metadata)


# ------------------------------------------------------------
# STEP 7: FINAL SUMMARY
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("QUALITY CONTROL SUMMARY")
print("=" * 60)

print("Total FASTA sequences:", total_sequences)

print(
    "Total unique sequences across all species:",
    len(all_sequence_owners)
)

print(
    "Cross-species exact duplicate sequences:",
    len(cross_species_duplicates)
)

print("\nQUALITY CHECK COMPLETED")
