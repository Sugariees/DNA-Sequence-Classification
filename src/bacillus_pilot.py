# Import Entrez to search and download records from NCBI
from Bio import Entrez

# Import SeqIO to read and write biological sequences
from Bio import SeqIO

# Import os to create directories
import os


# --------------------------------------------------
# STEP 1: CONNECT TO NCBI
# --------------------------------------------------

# Identify ourselves to NCBI using our email address
Entrez.email = os.environ.get("NCBI_EMAIL", "your_email@example.com")


# --------------------------------------------------
# STEP 2: SEARCH FOR BACILLUS SUBTILIS 16S rRNA
# --------------------------------------------------

# Search for Bacillus subtilis 16S rRNA records
# Select records with lengths between 1300 and 1600 bp
# Retrieve five candidate records for this pilot

search_term = (
    '"Bacillus subtilis"[Organism] '
    'AND 16S ribosomal RNA[Title] '
    'AND 1300:1600[SLEN]'
)

handle = Entrez.esearch(
    db="nucleotide",
    term=search_term,
    retmax=5
)

# Read the search results
results = Entrez.read(handle)

# Close the NCBI connection
handle.close()

# Extract the NCBI record IDs
sequence_ids = results["IdList"]

print("Total matching NCBI records:", results["Count"])
print("Candidate records selected:", len(sequence_ids))


# --------------------------------------------------
# STEP 3: DOWNLOAD GENBANK RECORDS
# --------------------------------------------------

if not sequence_ids:
    print("No candidate records found.")
    exit()

handle = Entrez.efetch(
    db="nucleotide",
    id=",".join(sequence_ids),
    rettype="gb",
    retmode="text"
)

# Parse the downloaded records
records = list(SeqIO.parse(handle, "genbank"))

# Close the connection
handle.close()


# --------------------------------------------------
# STEP 4: EXTRACT AND VALIDATE 16S rRNA SEQUENCES
# --------------------------------------------------

# Store accepted sequences here
accepted_records = []

# Store exact sequences to detect duplicates
seen_sequences = set()

for record in records:

    # Retrieve the organism annotation
    organism = record.annotations.get("organism", "").strip()

    # Confirm that the record belongs to Bacillus subtilis
    if organism.lower() != "bacillus subtilis":
        continue

    # Examine the annotations to locate the 16S rRNA gene
    for feature in record.features:

        if feature.type != "rRNA":
            continue

        product = " ".join(
            feature.qualifiers.get("product", [])
        ).lower()

        if "16s ribosomal rna" not in product:
            continue

        # Extract only the annotated 16S rRNA region
        gene_sequence = feature.extract(record.seq)

        sequence = str(gene_sequence).upper()

        # Check the extracted gene length
        if not 1300 <= len(sequence) <= 1600:
            continue

        # Accept only sequences containing A, T, G and C
        if not set(sequence).issubset({"A", "T", "G", "C"}):
            continue

        # Remove exact duplicate sequences
        if sequence in seen_sequences:
            continue

        seen_sequences.add(sequence)

        # Create a FASTA record containing the extracted gene
        from Bio.SeqRecord import SeqRecord

        extracted_record = SeqRecord(
            gene_sequence,
            id=record.id,
            description=f"{organism} 16S ribosomal RNA"
        )

        accepted_records.append(extracted_record)

        # Stop checking features after finding the 16S gene
        break


# --------------------------------------------------
# STEP 5: DISPLAY QUALITY CONTROL SUMMARY
# --------------------------------------------------

print("\nBACILLUS SUBTILIS PILOT SUMMARY")
print("-" * 40)

print("Downloaded GenBank records:", len(records))
print("Accepted sequences:", len(accepted_records))


# --------------------------------------------------
# STEP 6: SAVE THE PILOT DATASET
# --------------------------------------------------

# Create the pilot directory if it does not exist
os.makedirs("data/raw/pilot", exist_ok=True)

# Define the output FASTA file
output_file = "data/raw/pilot/B_subtilis_pilot_clean.fasta"

# Save the accepted sequences
SeqIO.write(accepted_records, output_file, "fasta")

print("\nSaved to:", output_file)
