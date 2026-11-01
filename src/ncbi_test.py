from Bio import Entrez, SeqIO

Entrez.email = os.environ.get("NCBI_EMAIL", "your_email@example.com")

accessions = [
    "QB220964.1",
    "QB226258.1",
    "QB200050.1"
]

for accession in accessions:

    print("\n" + "=" * 50)
    print("ACCESSION:", accession)

    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="gb",
        retmode="text"
    )

    record = SeqIO.read(handle, "genbank")
    handle.close()

    print("Organism:", record.annotations.get("organism"))
    print("Record length:", len(record.seq))

    for feature in record.features:

        if feature.type == "rRNA":

            product = feature.qualifiers.get("product", [])

            if any("16S" in item for item in product):

                print("Feature location:", feature.location)
                print("Product:", product)
                print("Gene length:", len(feature.extract(record.seq)))
