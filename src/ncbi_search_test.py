from Bio import Entrez, SeqIO

Entrez.email = os.environ.get("NCBI_EMAIL", "your_email@example.com")

handle = Entrez.esearch(
    db="nucleotide",
    term='"Escherichia coli"[Organism] AND 16S ribosomal RNA[Title] AND 1300:1600[SLEN]',
    retmax=5
)

results = Entrez.read(handle)
handle.close()

sequence_ids = results["IdList"]

handle = Entrez.efetch(
    db="nucleotide",
    id=",".join(sequence_ids),
    rettype="gb",
    retmode="text"
)

records = list(SeqIO.parse(handle, "genbank"))
handle.close()

for record in records:
    print("\nAccession:", record.id)
    print("Organism:", record.annotations.get("organism"))
    print("Record length:", len(record.seq))

    for feature in record.features:

        if feature.type == "rRNA":

            product = feature.qualifiers.get("product", [])

            if any("16S" in item for item in product):

                print("16S gene length:", len(feature.extract(record.seq)))
                print("Feature location:", feature.location)
import os

os.makedirs("data/raw/pilot", exist_ok=True)

SeqIO.write(
    records,
    "data/raw/pilot/E_coli_pilot.fasta",
    "fasta"
)

print("\nSaved five candidate sequences successfully.")
