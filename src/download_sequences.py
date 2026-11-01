from Bio import Entrez, SeqIO
import os
import time

# Identify ourselves to NCBI
Entrez.email = os.environ.get("NCBI_EMAIL", "your_email@example.com")

# Define the three bacterial species
species_list = {
    "Escherichia_coli": '"Escherichia coli"[Organism]',
    "Bacillus_subtilis": '"Bacillus subtilis"[Organism]',
    "Staphylococcus_aureus": '"Staphylococcus aureus"[Organism]'
}


def download_sequences(species_name, search_term, number=50):

    print(f"\nDownloading sequences for {species_name}...")

    # Search NCBI
    handle = Entrez.esearch(
        db="nucleotide",
        term=f'{search_term} AND 16S ribosomal RNA[Title] AND 1300:1600[SLEN]',
        retmax=number
    )

    results = Entrez.read(handle)
    handle.close()

    sequence_ids = results["IdList"]

    print("Records found:", results["Count"])
    print("Records selected:", len(sequence_ids))

    # Download sequences from NCBI
    handle = Entrez.efetch(
        db="nucleotide",
        id=",".join(sequence_ids),
        rettype="gb",
        retmode="text"
    )

    records = list(SeqIO.parse(handle, "genbank"))
    handle.close()

    print("Sequences downloaded:", len(records))

    # Filter and validate sequences
    valid_records = []
    seen_sequences = set()
    for record in records:

        # Verify organism
        organism = record.annotations.get("organism", "").strip()

        if not organism:
            continue

        if species_name.replace("_", " ").lower() != organism.lower():
            continue

        # Find the annotated 16S rRNA gene
        for feature in record.features:

            if feature.type != "rRNA":
                continue

            product = " ".join(feature.qualifiers.get("product", [])).lower()

            if "16s ribosomal rna" not in product:
                continue

            # Extract only the 16S rRNA gene region
            gene_sequence = feature.extract(record.seq)
            sequence = str(gene_sequence).upper()

            # Remove ambiguous bases
            if not set(sequence).issubset({"A", "T", "G", "C"}):
                continue

            # Remove sequences that are too short
            if len(sequence) < 500:
                continue

            # Remove exact duplicate gene sequences
            if sequence in seen_sequences:
                continue

            seen_sequences.add(sequence)

            from Bio.SeqRecord import SeqRecord

            extracted_record = SeqRecord(
                gene_sequence,
                id=record.id,
                description=f"{organism} 16S ribosomal RNA"
            )

            valid_records.append(extracted_record)
    print("Valid unique sequences:", len(valid_records))

    # Create output folder
    os.makedirs("data/raw", exist_ok=True)

    # Save sequences in FASTA format
    output_file = f"data/raw/{species_name}_16S.fasta"

    SeqIO.write(valid_records, output_file, "fasta")

    print("Saved to:", output_file)


# Download sequences for all three species
for species_name, search_term in species_list.items():

    download_sequences(species_name, search_term, number=5)

    time.sleep(1)
