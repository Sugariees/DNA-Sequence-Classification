# Import AlignIO to read multiple sequence alignments
from Bio import AlignIO

# Import combinations to compare every possible sequence pair


from itertools import combinations


# --------------------------------------------------
# STEP 1: LOAD THE ALIGNMENT
# --------------------------------------------------

# Read the alignment file provided as a command-line argument
import sys

if len(sys.argv) != 2:
    print("Usage: python src/pilot_similarity.py <alignment_file>")
    sys.exit(1)

alignment_file = sys.argv[1]

alignment = AlignIO.read(alignment_file, "fasta")


# --------------------------------------------------
# STEP 2: DISPLAY ALIGNMENT INFORMATION
# --------------------------------------------------

print("PILOT SEQUENCE SIMILARITY ANALYSIS")
print("-" * 50)

print("Number of sequences:", len(alignment))
print("Alignment length:", alignment.get_alignment_length())


# --------------------------------------------------
# STEP 3: COMPARE EVERY SEQUENCE PAIR
# --------------------------------------------------

for seq1, seq2 in combinations(alignment, 2):

    matches = 0
    compared = 0

    # Calculate original sequence lengths without alignment gaps
    length1 = len(seq1.seq) - seq1.seq.count("-")
    length2 = len(seq2.seq) - seq2.seq.count("-")

    # Compare positions where both sequences contain nucleotides
    for base1, base2 in zip(seq1.seq, seq2.seq):

        if base1 == "-" or base2 == "-":
            continue

        compared += 1

        if base1 == base2:
            matches += 1

    # Calculate percentage identity
    identity = (matches / compared) * 100

    # Calculate approximate coverage
    coverage1 = (compared / length1) * 100
    coverage2 = (compared / length2) * 100

    # Display pairwise results
    print(f"\n{seq1.id} vs {seq2.id}")
    print(f"Identity: {identity:.2f}%")
    print("Compared positions:", compared)

    print(f"Coverage of {seq1.id}: {coverage1:.2f}%")
    print(f"Coverage of {seq2.id}: {coverage2:.2f}%")
