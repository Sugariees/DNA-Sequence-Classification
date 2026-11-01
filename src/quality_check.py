from Bio import SeqIO
import glob

files = glob.glob("data/raw/*_16S.fasta")

print("=" * 55)
print("       16S rRNA DATASET QUALITY CONTROL")
print("=" * 55)

total_sequences = 0

for file in files:

    print("\nSpecies file:", file)

    records = list(SeqIO.parse(file, "fasta"))

    print("Number of sequences:", len(records))

    for record in records:

        sequence = str(record.seq).upper()

        total_sequences += 1

        print("\nAccession:", record.id)
        print("Length:", len(sequence), "bp")

        if set(sequence).issubset({"A", "T", "G", "C"}):
            print("Base quality: PASS")
        else:
            print("Base quality: FAIL")

        if 1300 <= len(sequence) <= 1600:
            print("Length check: PASS")
        else:
            print("Length check: REVIEW")

print("\n" + "=" * 55)
print("TOTAL SEQUENCES:", total_sequences)
print("=" * 55)
