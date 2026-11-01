# ============================================================
# MACHINE LEARNING DATASET PREPARATION
# PROJECT: BACTERIAL 16S rRNA CLASSIFICATION
#
# PURPOSE:
# 1. Load k-mer numerical features.
# 2. Load bacterial species labels.
# 3. Match features and labels using accession numbers.
# 4. Separate input features (X) and target labels (y).
# 5. Validate the prepared dataset.
# ============================================================


# STEP 1: IMPORT LIBRARIES

import pandas as pd
from pathlib import Path


# STEP 2: DEFINE FILE LOCATIONS

data_folder = Path("data/processed/final")

feature_file = data_folder / "kmer_features.csv"

label_file = data_folder / "sequence_labels.csv"

output_file = data_folder / "ml_dataset.csv"


# STEP 3: LOAD THE FEATURE MATRIX

print("\nLoading k-mer features...")

features_df = pd.read_csv(feature_file)

print("Feature data shape:", features_df.shape)


# STEP 4: LOAD SPECIES LABELS

print("\nLoading species labels...")

labels_df = pd.read_csv(label_file)

print("Label data shape:", labels_df.shape)


# STEP 5: CHECK FOR MISSING VALUES

print("\nChecking missing values...")

print("Missing feature values:", features_df.isnull().sum().sum())

print("Missing label values:", labels_df.isnull().sum().sum())


# STEP 6: CHECK ACCESSION NUMBERS

print("\nChecking accession numbers...")

# Each accession should appear only once in each file.
if features_df["accession"].duplicated().any():
    raise ValueError("Duplicate accession numbers in feature data.")

if labels_df["accession"].duplicated().any():
    raise ValueError("Duplicate accession numbers in label data.")

# Check whether both files contain exactly the same accessions.
feature_accessions = set(features_df["accession"])

label_accessions = set(labels_df["accession"])

if feature_accessions != label_accessions:
    raise ValueError("Feature and label accessions do not match.")

print("All accession numbers match.")


# STEP 7: COMBINE FEATURES AND LABELS

print("\nCombining features and labels...")

# Merge using accession numbers, not row positions.
dataset = pd.merge(
    features_df,
    labels_df,
    on="accession",
    how="inner",
    validate="one_to_one"
)

print("Combined dataset shape:", dataset.shape)


# STEP 8: CHECK SPECIES DISTRIBUTION

print("\nSpecies distribution:")

print(dataset["species"].value_counts())


# STEP 9: SEPARATE INPUT FEATURES (X) AND TARGET (y)

# X contains the 256 numerical k-mer frequency features.
X = dataset.drop(columns=["accession", "species"])

# y contains the correct bacterial species for each sequence.
y = dataset["species"]


# STEP 10: DISPLAY DATASET INFORMATION

print("\n" + "=" * 60)
print("MACHINE LEARNING DATASET PREPARATION")
print("=" * 60)

print("Total samples:", len(dataset))

print("Number of input features:", X.shape[1])

print("Number of target classes:", y.nunique())

print("\nInput feature matrix shape:", X.shape)

print("Target label shape:", y.shape)

print("\nTarget classes:")

for species in sorted(y.unique()):
    print("-", species)


# STEP 11: VALIDATE THE DATASET

assert len(X) == len(y), "Feature and label counts differ."

assert X.shape[1] == 256, "Unexpected number of k-mer features."

assert len(dataset) == 150, "Unexpected number of sequences."

assert not X.isnull().any().any(), "Missing values detected."

assert X.dtypes.apply(
    lambda dtype: pd.api.types.is_numeric_dtype(dtype)
).all(), "Non-numeric features detected."

print("\nAll dataset validation checks passed!")


# STEP 12: SAVE THE COMBINED DATASET

dataset.to_csv(output_file, index=False)

print("\nCombined dataset saved to:", output_file)

print("\nDATASET PREPARATION COMPLETED!")
