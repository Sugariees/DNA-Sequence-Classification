# ============================================================
# FINAL MULTIPLE SEQUENCE ALIGNMENT
# PROJECT: BACTERIAL 16S rRNA CLASSIFICATION
#
# PURPOSE:
# Align the 16S rRNA sequences of each bacterial species
# using the MAFFT alignment program.
# ============================================================


# ------------------------------------------------------------
# STEP 1: IMPORT LIBRARIES
# ------------------------------------------------------------

import subprocess
from pathlib import Path


# ------------------------------------------------------------
# STEP 2: DEFINE INPUT AND OUTPUT FOLDERS
# ------------------------------------------------------------

input_folder = Path("data/processed/final")

output_folder = Path("results/final/alignment")

# Create the output folder if it does not exist
output_folder.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# STEP 3: DEFINE SPECIES FILES
# ------------------------------------------------------------

species_files = {
    "Escherichia_coli":
        "Escherichia_coli_16S_clean.fasta",

    "Bacillus_subtilis":
        "Bacillus_subtilis_16S_clean.fasta",

    "Staphylococcus_aureus":
        "Staphylococcus_aureus_16S_clean.fasta"
}


# ------------------------------------------------------------
# STEP 4: RUN MAFFT FOR EACH SPECIES
# ------------------------------------------------------------

for species, filename in species_files.items():

    print("\n" + "=" * 60)
    print("ALIGNING:", species)
    print("=" * 60)

    # Define the input FASTA file
    input_file = input_folder / filename

    # Define the output alignment file
    output_file = (
        output_folder /
        f"{species}_aligned.fasta"
    )

    # Stop if the input file does not exist
    if not input_file.exists():

        print("ERROR: Input file not found:", input_file)
        continue

    # Run MAFFT using its automatic algorithm selection
    command = [
        "mafft",
        "--auto",
        "--thread",
        "1",
        str(input_file)
    ]

    try:

        # Execute MAFFT and capture its alignment output
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True
        )

        # Save the alignment as a FASTA file
        with open(
            output_file,
            "w",
            encoding="utf-8"
        ) as output_handle:

            output_handle.write(result.stdout)

        print("Alignment completed successfully.")

        print("Saved to:", output_file)

        # Display MAFFT diagnostic messages if available
        if result.stderr.strip():

            print("\nMAFFT information:")
            print(result.stderr.strip())

    except subprocess.CalledProcessError as error:

        print("MAFFT returned an error:")
        print(error.stderr)

    except FileNotFoundError:

        print("ERROR: MAFFT is not installed or not in PATH.")


# ------------------------------------------------------------
# STEP 5: COMPLETION MESSAGE
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("ALL ALIGNMENT TASKS FINISHED")
print("=" * 60)
