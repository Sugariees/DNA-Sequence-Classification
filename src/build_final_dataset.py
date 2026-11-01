# ============================================================
# PROJECT: BACTERIAL 16S rRNA SEQUENCE CLASSIFICATION
# FILE: build_final_dataset.py
#
# PURPOSE:
# Collect and validate 16S rRNA gene sequences from NCBI.
# Reuse previously collected pilot sequences.
# Save accepted sequences and metadata for further analysis.
# ============================================================


# ------------------------------------------------------------
# STEP 1: IMPORT REQUIRED LIBRARIES
# ------------------------------------------------------------

from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import os
import csv
import time
import random
import hashlib


# ------------------------------------------------------------
# STEP 2: CONFIGURATION
# ------------------------------------------------------------

# Identify ourselves when accessing NCBI
Entrez.email = os.environ.get("NCBI_EMAIL", "your_email@example.com")

# Number of sequences we ultimately want per species
TARGET_PER_SPECIES = 50

# Maximum number of new NCBI candidates to examine per species
MAX_CANDIDATES = 200

# Number of GenBank records downloaded in one request
BATCH_SIZE = 20

# Preliminary 16S sequence length filter
MIN_LENGTH = 1300
MAX_LENGTH = 1600

# Fixed random seed makes candidate selection reproducible
RANDOM_SEED = 42

# Define the bacterial species and their NCBI search terms
SPECIES = {
    "Escherichia_coli": '"Escherichia coli"[Organism]',
    "Bacillus_subtilis": '"Bacillus subtilis"[Organism]',
    "Staphylococcus_aureus": '"Staphylococcus aureus"[Organism]'
}

# Pilot FASTA files that we want to reuse
PILOT_FILES = {
    "Escherichia_coli":
        "data/raw/pilot/E_coli_pilot_clean.fasta",

    "Bacillus_subtilis":
        "data/raw/pilot/B_subtilis_pilot_clean.fasta",

    "Staphylococcus_aureus":
        "data/raw/pilot/S_aureus_pilot_clean.fasta"
}

# Create output folders
RAW_FOLDER = "data/raw/final"
PROCESSED_FOLDER = "data/processed/final"

os.makedirs(RAW_FOLDER, exist_ok=True)
os.makedirs(PROCESSED_FOLDER, exist_ok=True)


# ------------------------------------------------------------
# STEP 3: STORAGE FOR RESULTS
# ------------------------------------------------------------

# Store accepted sequence metadata
all_metadata = []

# Store rejected candidate information
all_rejections = []


# ------------------------------------------------------------
# STEP 4: CREATE A SEQUENCE FINGERPRINT
# ------------------------------------------------------------

def sequence_fingerprint(sequence):

    # SHA-256 generates a unique fingerprint for a sequence
    return hashlib.sha256(
        sequence.encode("utf-8")
    ).hexdigest()


# ------------------------------------------------------------
# STEP 5: EXTRACT STRAIN AND ISOLATE INFORMATION
# ------------------------------------------------------------

def get_source_information(record):

    strain = ""
    isolate = ""

    # Source features contain organism and isolate information
    for feature in record.features:

        if feature.type == "source":

            qualifiers = feature.qualifiers

            strain = "; ".join(
                qualifiers.get("strain", [])
            )

            isolate = "; ".join(
                qualifiers.get("isolate", [])
            )

            break

    return strain, isolate


# ------------------------------------------------------------
# STEP 6: IDENTIFY THE 16S rRNA FEATURE
# ------------------------------------------------------------

def find_16s_feature(record):

    # Search the GenBank annotations for 16S rRNA
    matching_features = []

    for feature in record.features:

        if feature.type != "rRNA":
            continue

        product = " ".join(
            feature.qualifiers.get("product", [])
        ).lower()

        if "16s ribosomal rna" in product:
            matching_features.append(feature)

    # Return a feature only when exactly one is found
    if len(matching_features) == 1:
        return matching_features[0]

    return None


# ------------------------------------------------------------
# STEP 7: PROCESS ONE GENBANK RECORD
# ------------------------------------------------------------

def process_record(
    record,
    species_name,
    accepted_records,
    seen_sequences,
    source_type
):

    # Obtain the organism annotation
    organism = record.annotations.get("organism", "")

    expected_organism = species_name.replace(
        "_", " "
    ).lower()

    # Verify that the record belongs to the expected species
    if organism.lower() != expected_organism:

        all_rejections.append({
            "species": species_name,
            "accession": record.id,
            "reason": "Organism mismatch"
        })

        return

    # Find the annotated 16S rRNA feature
    feature = find_16s_feature(record)

    if feature is None:

        all_rejections.append({
            "species": species_name,
            "accession": record.id,
            "reason": "Missing or ambiguous 16S annotation"
        })

        return

    # Extract only the annotated gene sequence
    extracted_sequence = feature.extract(record.seq)

    sequence = str(extracted_sequence).upper()

    # Check the extracted sequence length
    if not MIN_LENGTH <= len(sequence) <= MAX_LENGTH:

        all_rejections.append({
            "species": species_name,
            "accession": record.id,
            "reason": "Length outside accepted range"
        })

        return

    # Reject sequences containing ambiguous nucleotides
    if not set(sequence).issubset({"A", "T", "G", "C"}):

        all_rejections.append({
            "species": species_name,
            "accession": record.id,
            "reason": "Ambiguous nucleotide present"
        })

        return

    # Calculate a fingerprint for duplicate detection
    fingerprint = sequence_fingerprint(sequence)

    # Remove exact duplicates within this species
    if fingerprint in seen_sequences:

        all_rejections.append({
            "species": species_name,
            "accession": record.id,
            "reason": "Exact duplicate sequence"
        })

        return

    # Save the fingerprint
    seen_sequences.add(fingerprint)

    # Extract strain and isolate annotations
    strain, isolate = get_source_information(record)

    # Store the extracted sequence as a FASTA record
    fasta_record = SeqRecord(
        extracted_sequence,
        id=record.id,
        description=f"{organism} 16S rRNA"
    )

    accepted_records.append(fasta_record)

    # Record whether the feature has partial boundaries
    feature_location = str(feature.location)

    partial_boundary = (
        "<" in feature_location or ">" in feature_location
    )

    # Store information about the accepted sequence
    all_metadata.append({

        "species": organism,
        "accession": record.id,
        "strain": strain,
        "isolate": isolate,
        "sequence_length": len(sequence),
        "source_type": source_type,
        "feature_location": feature_location,
        "partial_boundary": partial_boundary,
        "record_description": record.description,
        "sha256": fingerprint

    })


# ------------------------------------------------------------
# STEP 8: LOAD PILOT ACCESSION NUMBERS
# ------------------------------------------------------------

def get_pilot_accessions(pilot_file):

    # Return an empty list if a pilot file is missing
    if not os.path.exists(pilot_file):
        print("Pilot file not found:", pilot_file)
        return []

    # Read previously collected FASTA sequences
    pilot_records = list(
        SeqIO.parse(pilot_file, "fasta")
    )

    # Extract their accession numbers
    accessions = [
        record.id for record in pilot_records
    ]

    return accessions


# ------------------------------------------------------------
# STEP 9: SEARCH NCBI AND COLLECT CANDIDATES
# ------------------------------------------------------------

def collect_species(species_name, search_term):

    print("\n" + "=" * 60)
    print("PROCESSING:", species_name)
    print("=" * 60)

    # Store accepted sequences for this species
    accepted_records = []

    # Store exact sequence fingerprints for this species
    seen_sequences = set()

    # Track accession numbers already processed
    processed_accessions = set()

    # --------------------------------------------------------
    # PART A: REUSE PILOT SEQUENCES
    # --------------------------------------------------------

    print("\nChecking existing pilot sequences...")

    pilot_file = PILOT_FILES[species_name]

    pilot_accessions = get_pilot_accessions(
        pilot_file
    )

    print("Pilot accessions found:", len(pilot_accessions))

    # Download the corresponding GenBank records
    for start in range(0, len(pilot_accessions), BATCH_SIZE):

        batch = pilot_accessions[
            start:start + BATCH_SIZE
        ]

        handle = Entrez.efetch(
            db="nucleotide",
            id=",".join(batch),
            rettype="gb",
            retmode="text"
        )

        pilot_records = list(
            SeqIO.parse(handle, "genbank")
        )

        handle.close()

        # Process every pilot record
        for record in pilot_records:

            if record.id in processed_accessions:
                continue

            processed_accessions.add(record.id)

            process_record(
                record,
                species_name,
                accepted_records,
                seen_sequences,
                "Reused pilot"
            )

        time.sleep(1)

    print(
        "Accepted after pilot reuse:",
        len(accepted_records)
    )

    # --------------------------------------------------------
    # PART B: SEARCH NCBI FOR NEW CANDIDATES
    # --------------------------------------------------------

    # Stop if the target has already been reached
    if len(accepted_records) < TARGET_PER_SPECIES:

        print("\nSearching NCBI for additional candidates...")

        query = (
            search_term
            + " AND 16S ribosomal RNA[Title]"
            + f" AND {MIN_LENGTH}:{MAX_LENGTH}[SLEN]"
        )

        # Obtain the total number of search results
        handle = Entrez.esearch(
            db="nucleotide",
            term=query,
            retmax=0
        )

        search_results = Entrez.read(handle)
        handle.close()

        total_available = int(
            search_results["Count"]
        )

        print(
            "Total matching NCBI records:",
            total_available
        )

        # Retrieve all candidate IDs in pages
        all_ids = []

        PAGE_SIZE = 10000

        for start in range(
            0,
            total_available,
            PAGE_SIZE
        ):

            page_size = min(
                PAGE_SIZE,
                total_available - start
            )

            handle = Entrez.esearch(
                db="nucleotide",
                term=query,
                retstart=start,
                retmax=page_size
            )

            page_results = Entrez.read(handle)
            handle.close()

            all_ids.extend(
                page_results["IdList"]
            )

            time.sleep(1)

        print("Candidate IDs collected:", len(all_ids))

        # Select a reproducible random sample from the results
        random_generator = random.Random(
            RANDOM_SEED + len(species_name)
        )

        candidate_limit = min(
            MAX_CANDIDATES,
            len(all_ids)
        )

        selected_ids = random_generator.sample(
            all_ids,
            candidate_limit
        )

        print(
            "New candidate records selected:",
            len(selected_ids)
        )

        # ----------------------------------------------------
        # PART C: DOWNLOAD AND PROCESS NEW RECORDS
        # ----------------------------------------------------

        for start in range(
            0,
            len(selected_ids),
            BATCH_SIZE
        ):

            # Stop when the target is reached
            if len(accepted_records) >= TARGET_PER_SPECIES:
                break

            batch = selected_ids[
                start:start + BATCH_SIZE
            ]

            handle = Entrez.efetch(
                db="nucleotide",
                id=",".join(batch),
                rettype="gb",
                retmode="text"
            )

            downloaded_records = list(
                SeqIO.parse(handle, "genbank")
            )

            handle.close()

            # Process each downloaded record
            for record in downloaded_records:

                # Avoid processing the same accession twice
                if record.id in processed_accessions:
                    continue

                processed_accessions.add(record.id)

                process_record(
                    record,
                    species_name,
                    accepted_records,
                    seen_sequences,
                    "New NCBI candidate"
                )

                if len(accepted_records) >= TARGET_PER_SPECIES:
                    break

            print(
                "Accepted sequences:",
                len(accepted_records)
            )

            # Pause between NCBI requests
            time.sleep(1)

    # --------------------------------------------------------
    # PART D: SAVE ACCEPTED FASTA FILE
    # --------------------------------------------------------

    output_file = os.path.join(
        PROCESSED_FOLDER,
        f"{species_name}_16S_clean.fasta"
    )

    SeqIO.write(
        accepted_records,
        output_file,
        "fasta"
    )

    # --------------------------------------------------------
    # PART E: DISPLAY SUMMARY
    # --------------------------------------------------------

    print("\nSPECIES COLLECTION SUMMARY")
    print("-" * 40)

    print("Accepted sequences:", len(accepted_records))
    print("Target sequences:", TARGET_PER_SPECIES)

    print("Saved FASTA file:", output_file)

    if len(accepted_records) < TARGET_PER_SPECIES:

        print(
            "NOTE: Fewer than 50 sequences passed the filters."
        )


# ------------------------------------------------------------
# STEP 10: RUN COLLECTION FOR ALL THREE SPECIES
# ------------------------------------------------------------

for species_name, search_term in SPECIES.items():

    collect_species(
        species_name,
        search_term
    )

    time.sleep(1)


# ------------------------------------------------------------
# STEP 11: SAVE ACCEPTED SEQUENCE METADATA
# ------------------------------------------------------------

metadata_file = os.path.join(
    PROCESSED_FOLDER,
    "sequence_metadata.csv"
)

metadata_columns = [
    "species",
    "accession",
    "strain",
    "isolate",
    "sequence_length",
    "source_type",
    "feature_location",
    "partial_boundary",
    "record_description",
    "sha256"
]

with open(
    metadata_file,
    "w",
    newline="",
    encoding="utf-8"
) as csv_file:

    writer = csv.DictWriter(
        csv_file,
        fieldnames=metadata_columns
    )

    writer.writeheader()
    writer.writerows(all_metadata)


# ------------------------------------------------------------
# STEP 12: SAVE REJECTION REPORT
# ------------------------------------------------------------

rejection_file = os.path.join(
    PROCESSED_FOLDER,
    "rejected_candidates.csv"
)

with open(
    rejection_file,
    "w",
    newline="",
    encoding="utf-8"
) as csv_file:

    writer = csv.DictWriter(
        csv_file,
        fieldnames=[
            "species",
            "accession",
            "reason"
        ]
    )

    writer.writeheader()
    writer.writerows(all_rejections)


# ------------------------------------------------------------
# STEP 13: FINAL DATASET SUMMARY
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("FINAL DATASET SUMMARY")
print("=" * 60)

for species_name in SPECIES:

    expected_name = species_name.replace(
        "_", " "
    )

    count = sum(
        1 for item in all_metadata
        if item["species"].lower() == expected_name.lower()
    )

    print(
        species_name,
        ":",
        count,
        "accepted sequences"
    )

print("\nTotal accepted sequences:", len(all_metadata))

print("\nMetadata file:", metadata_file)
print("Rejection report:", rejection_file)

print("\nDATASET COLLECTION FINISHED")
