from Bio import AlignIO

alignment = AlignIO.read(
    "results/alignment/S_aureus_aligned.fasta",
    "fasta"
)

print("Number of sequences:", len(alignment))
print("Alignment length:", alignment.get_alignment_length())

print("\nGAP ANALYSIS")

for record in alignment:

    sequence = str(record.seq)

    gaps = sequence.count("-")

    print(record.id)
    print("Gaps:", gaps)
    print("Nucleotides:", len(sequence) - gaps)
