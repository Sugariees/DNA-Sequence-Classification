# ============================================================
# PURPOSE: Inspect CD-HIT-EST cluster sizes
# ============================================================

from pathlib import Path
import re
import pandas as pd


# SECTION 1: Define directories

input_directory = Path("data/processed/diversity_controlled")

output_directory = Path("results/final/reports")

output_directory.mkdir(parents=True, exist_ok=True)


# SECTION 2: Define cluster files

species_files = [
    "Escherichia_coli_clustered.fasta.clstr",
    "Bacillus_subtilis_clustered.fasta.clstr",
    "Staphylococcus_aureus_clustered.fasta.clstr"
]


# SECTION 3: Read cluster information

all_results = []

for filename in species_files:

    cluster_file = input_directory / filename

    species = filename.replace("_clustered.fasta.clstr", "")

    print("\nSPECIES:", species)
    print("=" * 60)

    cluster_sizes = []
    current_cluster_size = 0

    with cluster_file.open() as file:

        for line in file:

            if line.startswith(">Cluster"):

                if current_cluster_size > 0:
                    cluster_sizes.append(current_cluster_size)

                current_cluster_size = 0

            else:
                current_cluster_size += 1

        if current_cluster_size > 0:
            cluster_sizes.append(current_cluster_size)

    for number, size in enumerate(cluster_sizes, start=1):

        print(f"Cluster {number:02d}: {size} sequences")

        all_results.append({
            "species": species,
            "cluster_id": f"Cluster_{number:02d}",
            "cluster_size": size
        })


# SECTION 4: Save results

results = pd.DataFrame(all_results)

results.to_csv(
    output_directory / "cdhit_cluster_sizes.csv",
    index=False
)


# SECTION 5: Display summary

print("\nCLUSTER SUMMARY")
print("=" * 60)

print(
    results.groupby("species")["cluster_size"]
    .agg(["count", "min", "max", "mean", "sum"])
)

print("\nReport saved: cdhit_cluster_sizes.csv")
print("\nANALYSIS COMPLETED!")
