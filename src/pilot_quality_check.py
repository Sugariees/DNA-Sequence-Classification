# Import SeqIO to read biological sequences from FASTA files
from Bio import SeqIO

# Import os to check whether the input file exists
import os


# --------------------------------------------------
# STEP 1: DEFINE INPUT FILE
# --------------------------------------------------

# Location of the E. coli pilot FASTA file
input_file = "data/raw/pilot/E_coli_pilot.fasta"


# --------------------------------------------------
# STEP 2: CHECK FILE EXISTENCE
# --------------------------------------------------

# Stop the program if the input file does not exist
if not os.path.exists(input_file):

    print("ERROR: FASTA file not found.")

    exit()


# --------------------------------------------------
# STEP 3: DEFINE VALID DNA BASES
# --------------------------------------------------

# Only these four nucleotides are accepted
valid_bases = {"A", "T", "G", "C"}


# --------------------------------------------------
# STEP 4: READ AND CHECK SEQUENCES
# --------------------------------------------------

# Read all sequences from the FASTA file
records = list(SeqIO.parse(input_file, "fasta"))

# Store sequences that pass quality control
accepted_records = []

# Store sequences that fail quality control
rejected_records = []


# Examine each sequence individually
for record in records:

    # Convert sequence into uppercase text
    sequence = str(record.seq).upper()

    # Check whether every nucleotide belongs to A, T, G or C
    if set(sequence).issubset(valid_bases):

        print(record.id, "PASS")

        accepted_records.append(record)

    else:

        print(record.id, "FAIL: Contains ambiguous bases")

        rejected_records.append(record)


# --------------------------------------------------
# STEP 5: DISPLAY QUALITY CONTROL SUMMARY
# --------------------------------------------------

print("\nQUALITY CONTROL SUMMARY")
print("-" * 35)

print("Total sequences:", len(records))
print("Accepted sequences:", len(accepted_records))
print("Rejected sequences:", len(rejected_records))


# --------------------------------------------------
# STEP 6: SAVE ACCEPTED SEQUENCES
# --------------------------------------------------

# Save accepted sequences separately
output_file = "data/raw/pilot/E_coli_pilot_clean.fasta"

SeqIO.write(accepted_records, output_file, "fasta")

print("\nClean sequences saved to:", output_file)
