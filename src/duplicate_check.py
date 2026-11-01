from Bio import SeqIO
import glob

sequence_database = {}

files = glob.glob("data/raw/*_16S.fasta")

for file in files:

    records = SeqIO.parse(file, "fasta")

    for record in records:

        sequence = str(record.seq).upper()

        if sequence in sequence_database:

            print("DUPLICATE FOUND!")
            print("Current:", record.id)
            print("Previously found:", sequence_database[sequence])

        else:

            sequence_database[sequence] = record.id

print("\nTotal unique sequences:", len(sequence_database))
