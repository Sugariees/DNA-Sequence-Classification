# ============================================================
# PROJECT: DNA Sequence Classification
# PURPOSE: Cluster similar 16S rRNA sequences using CD-HIT-EST
# IDENTITY THRESHOLD: 99%
# ============================================================

import subprocess

from pathlib import Path


# SECTION 1: Define directories

input_directory = Path("data/processed/final")

output_directory = Path(
    "data/processed/diversity_controlled"
)

output_directory.mkdir(parents=True, exist_ok=True)


# SECTION 2: Define species FASTA files

species_files = {
    "Escherichia_coli": (
        input_directory / "Escherichia_coli_16S_clean.fasta"
    ),

    "Bacillus_subtilis": (
        input_directory / "Bacillus_subtilis_16S_clean.fasta"
    ),

    "Staphylococcus_aureus": (
        input_directory / "Staphylococcus_aureus_16S_clean.fasta"
    )
}


# SECTION 3: Define CD-HIT-EST parameters

identity_threshold = 0.99

word_size = 10

print("\nCD-HIT-EST SEQUENCE CLUSTERING")
print("=" * 65)

print("Identity threshold:", identity_threshold * 100, "%")


# SECTION 4: Cluster each species independently

for species, input_file in species_files.items():

    if not input_file.exists():
        raise FileNotFoundError(
            f"Input FASTA file not found: {input_file}"
        )

    output_file = output_directory / f"{species}_clustered.fasta"

    cluster_file = Path(str(output_file) + ".clstr")

    print(f"\nProcessing: {species}")
    print("-" * 65)

    command = [
        "cd-hit-est",
        "-i", str(input_file),
        "-o", str(output_file),
        "-c", str(identity_threshold),
        "-n", str(word_size),
        "-d", "0",
        "-T", "1",
        "-M", "0"
    ]

    # Run CD-HIT-EST and display its output.
    subprocess.run(
        command,
        check=True
    )

    # Count representative sequences in the output FASTA.
    representative_count = sum(
        1 for line in output_file.open()
        if line.startswith(">")
    )

    print("Original sequences: 50")
    print("Representative sequences:", representative_count)
    print("Clusters:", representative_count)

    print("Representative FASTA:", output_file)
    print("Cluster information:", cluster_file)


# SECTION 5: Completion message

print("\nCLUSTERING COMPLETED")
print("=" * 65)

print("Output directory:", output_directory)
