from Bio import AlignIO
from itertools import combinations

alignment = AlignIO.read(
    "results/alignment/S_aureus_aligned.fasta",
    "fasta"
)

print("PAIRWISE SEQUENCE IDENTITY AND COVERAGE")
print("-" * 55)

for seq1, seq2 in combinations(alignment, 2):

    matches = 0
    compared = 0

    length1 = len(seq1.seq) - seq1.seq.count("-")
    length2 = len(seq2.seq) - seq2.seq.count("-")

    for base1, base2 in zip(seq1.seq, seq2.seq):

        if base1 == "-" or base2 == "-":
            continue

        compared += 1

        if base1 == base2:
            matches += 1

    identity = (matches / compared) * 100

    coverage1 = (compared / length1) * 100
    coverage2 = (compared / length2) * 100

    print(f"\n{seq1.id} vs {seq2.id}")
    print(f"Identity: {identity:.2f}%")
    print(f"Compared positions: {compared}")
    print(f"Coverage of {seq1.id}: {coverage1:.2f}%")
    print(f"Coverage of {seq2.id}: {coverage2:.2f}%")
