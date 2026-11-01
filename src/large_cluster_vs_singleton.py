# ============================================================
# PROJECT: DNA Sequence Classification
# SCRIPT: Large Cluster vs Singleton Evaluation
#
# PURPOSE:
# Train a Random Forest using sequences from the largest
# similarity clusters and test it on singleton sequences.
#
# TRAINING:
# Sequences belonging to groups with more than one member.
#
# TESTING:
# Sequences belonging to singleton groups.
#
# IMPORTANT:
# This is an exploratory evaluation, not independent
# biological validation.
# ============================================================


# -------------------- IMPORT LIBRARIES ----------------------

from pathlib import Path

import pandas as pd
import joblib

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix
)


# -------------------- PROJECT PATHS --------------------------

PROJECT_DIR = Path(__file__).resolve().parent.parent

DATA_FILE = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "final"
    / "ml_dataset.csv"
)

GROUP_FILE = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "final"
    / "sequence_groups.csv"
)

OUTPUT_DIR = (
    PROJECT_DIR
    / "results"
    / "final"
    / "singleton_evaluation"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# -------------------- LOAD DATA ------------------------------

dataset = pd.read_csv(DATA_FILE)
groups_data = pd.read_csv(GROUP_FILE)

print("ML dataset:", dataset.shape)
print("Group dataset:", groups_data.shape)


# -------------------- MATCH GROUP INFORMATION ----------------

# Match each sequence with its similarity group.

merged = dataset.merge(
    groups_data,
    on=["accession", "species"],
    how="left",
    validate="one_to_one"
)

if merged["group_id"].isna().any():
    raise ValueError("Missing group assignments detected.")


# -------------------- CALCULATE GROUP SIZES ------------------

# Count how many sequences belong to each similarity group.

group_sizes = (
    merged.groupby("group_id")["accession"]
    .transform("count")
)

merged["group_size"] = group_sizes


# -------------------- SEPARATE TRAINING AND TESTING DATA ------

# Training: sequences belonging to groups with multiple members.
# Testing: sequences belonging to singleton groups.

train_data = merged[merged["group_size"] > 1].copy()

test_data = merged[merged["group_size"] == 1].copy()


# -------------------- PREPARE FEATURES AND LABELS ------------

excluded_columns = [
    "accession",
    "species",
    "group_id",
    "group_size"
]

X_train = train_data.drop(columns=excluded_columns)
y_train = train_data["species"]

X_test = test_data.drop(columns=excluded_columns)
y_test = test_data["species"]


# -------------------- DISPLAY DATASET INFORMATION ------------

print("\nTraining sequences:", len(X_train))
print("Testing sequences:", len(X_test))

print("\nTraining class distribution:")
print(y_train.value_counts())

print("\nTesting class distribution:")
print(y_test.value_counts())


# -------------------- INITIALIZE RANDOM FOREST ---------------

model = RandomForestClassifier(
    n_estimators=300,
    random_state=42,
    class_weight="balanced"
)


# -------------------- TRAIN MODEL ----------------------------

print("\nTraining Random Forest...")

model.fit(X_train, y_train)

print("Training completed!")


# -------------------- PREDICT TEST SEQUENCES -----------------

predictions = model.predict(X_test)


# -------------------- EVALUATE MODEL -------------------------

accuracy = accuracy_score(y_test, predictions)

print("\n==========================================")
print("SINGLETON EVALUATION RESULTS")
print("==========================================")

print(f"\nAccuracy: {accuracy:.4f}")

print("\nClassification Report:")

report = classification_report(
    y_test,
    predictions,
    zero_division=0
)

print(report)


# -------------------- CONFUSION MATRIX -----------------------

labels = sorted(dataset["species"].unique())

matrix = confusion_matrix(
    y_test,
    predictions,
    labels=labels
)

print("\nConfusion Matrix:")

print(pd.DataFrame(
    matrix,
    index=labels,
    columns=labels
))


# -------------------- SAVE RESULTS ---------------------------

pd.DataFrame({
    "accession": test_data["accession"].values,
    "actual_species": y_test.values,
    "predicted_species": predictions,
    "correct": y_test.values == predictions
}).to_csv(
    OUTPUT_DIR / "singleton_predictions.csv",
    index=False
)

pd.DataFrame(
    matrix,
    index=labels,
    columns=labels
).to_csv(
    OUTPUT_DIR / "confusion_matrix.csv"
)

pd.DataFrame(
    classification_report(
        y_test,
        predictions,
        output_dict=True,
        zero_division=0
    )
).transpose().to_csv(
    OUTPUT_DIR / "classification_report.csv"
)

# Save the model separately
joblib.dump(
    model,
    OUTPUT_DIR / "singleton_evaluation_model.pkl"
)


# -------------------- FINAL MESSAGE -------------------------

print("\nResults saved in:")
print(OUTPUT_DIR)

print("\nSingleton evaluation completed!")
