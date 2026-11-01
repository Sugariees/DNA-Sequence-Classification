# Import SeqIO to read FASTA sequences
from Bio import SeqIO

# Import combinations to compare sequences in pairs
from itertools import combinations


# --------------------------------------------------
# STEP 1: LOAD CLEAN PILOT DATASET
# --------------------------------------------------

input_file = "data/raw/pilot/E_coli_pilot_clean.fasta"

records = list(SeqIO.parse(input_file, "fasta"))


# --------------------------------------------------
# STEP 2: CHECK EXACT DUPLICATES
# --------------------------------------------------

print("EXACT DUPLICATE CHECK")
print("-" * 40)

duplicate_found = False

# Compare every possible pair of sequences
for seq1, seq2 in combinations(records, 2):

    # Convert sequences into uppercase strings
    sequence1 = str(seq1.seq).upper()
    sequence2 = str(seq2.seq).upper()

    # Check whether both sequences are identical
    if sequence1 == sequence2:

        print("Duplicate found:", seq1.id, "and", seq2.id)

        duplicate_found = True


# --------------------------------------------------
# STEP 3: DISPLAY RESULTS
# --------------------------------------------------

if not duplicate_found:

    print("No exact duplicates found.")

print("\nTotal sequences checked:", len(records))
