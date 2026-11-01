# ============================================================
# PROJECT: DNA Sequence Classification
# SCRIPT: Train Random Forest on Representative Dataset
#
# PURPOSE:
# Train and evaluate a Random Forest classifier using
# the 37 CD-HIT representative sequences.
# ============================================================


# -------------------- IMPORT LIBRARIES ----------------------

from pathlib import Path

import pandas as pd
import joblib

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
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
    / "diversity_controlled"
    / "representative_ml_dataset.csv"
)

OUTPUT_DIR = (
    PROJECT_DIR
    / "results"
    / "final"
    / "representative_model"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# -------------------- LOAD DATASET ---------------------------

dataset = pd.read_csv(DATA_FILE)

print("Dataset loaded successfully.")
print("Dataset shape:", dataset.shape)


# -------------------- PREPARE FEATURES AND LABELS ------------

# Remove accession, species and sequence length.
# Only the 256 numerical k-mer features are used for training.

X = dataset.drop(
    columns=["accession", "species", "sequence_length"]
)

# Species names are the classification labels.
y = dataset["species"]

print("\nFeature matrix:", X.shape)
print("Label vector:", y.shape)

print("\nClass distribution:")
print(y.value_counts())


# -------------------- SPLIT DATASET --------------------------

# Reserve 20% for testing.
# Stratification preserves approximate class proportions.

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    random_state=42,
    stratify=y
)

print("\nTraining sequences:", len(X_train))
print("Testing sequences:", len(X_test))


# -------------------- INITIALIZE RANDOM FOREST ---------------

model = RandomForestClassifier(
    n_estimators=300,
    random_state=42,
    class_weight="balanced"
)


# -------------------- TRAIN MODEL ----------------------------

print("\nTraining Random Forest model...")

model.fit(X_train, y_train)

print("Training completed!")


# -------------------- MAKE PREDICTIONS -----------------------

predictions = model.predict(X_test)


# -------------------- EVALUATE MODEL -------------------------

accuracy = accuracy_score(y_test, predictions)

print("\n==========================================")
print("REPRESENTATIVE MODEL RESULTS")
print("==========================================")

print(f"\nAccuracy: {accuracy:.4f}")

print("\nClassification Report:")

report = classification_report(
    y_test,
    predictions,
    zero_division=0
)

print(report)

print("\nConfusion Matrix:")

labels = sorted(y.unique())

matrix = confusion_matrix(
    y_test,
    predictions,
    labels=labels
)

print(pd.DataFrame(
    matrix,
    index=labels,
    columns=labels
))


# -------------------- SAVE RESULTS ---------------------------

# Save classification report
report_dict = classification_report(
    y_test,
    predictions,
    output_dict=True,
    zero_division=0
)

pd.DataFrame(report_dict).transpose().to_csv(
    OUTPUT_DIR / "classification_report.csv"
)

# Save confusion matrix
pd.DataFrame(
    matrix,
    index=labels,
    columns=labels
).to_csv(OUTPUT_DIR / "confusion_matrix.csv")

# Save individual predictions
pd.DataFrame({
    "actual_species": y_test.values,
    "predicted_species": predictions
}).to_csv(
    OUTPUT_DIR / "predictions.csv",
    index=False
)

# Save trained model
joblib.dump(
    model,
    OUTPUT_DIR / "representative_random_forest.pkl"
)


# -------------------- FINAL MESSAGE -------------------------

print("\nResults saved in:")
print(OUTPUT_DIR)

print("\nRepresentative model training completed!")
