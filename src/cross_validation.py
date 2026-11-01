# ============================================================
# PROJECT: DNA Sequence Classification
# SCRIPT: Random Forest Cross-Validation
#
# PURPOSE:
# Evaluate the Random Forest model using 5-fold
# stratified cross-validation.
# ============================================================


# -------------------- IMPORT LIBRARIES ----------------------

from pathlib import Path

import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate


# -------------------- PROJECT PATHS --------------------------

PROJECT_DIR = Path(__file__).resolve().parent.parent

DATA_FILE = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "final"
    / "ml_dataset.csv"
)

OUTPUT_DIR = (
    PROJECT_DIR
    / "results"
    / "final"
    / "cross_validation"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# -------------------- LOAD DATASET ---------------------------

dataset = pd.read_csv(DATA_FILE)

print("Dataset loaded:", dataset.shape)


# -------------------- PREPARE FEATURES AND LABELS ------------

X = dataset.drop(
    columns=["accession", "species"]
)

y = dataset["species"]

print("Feature matrix:", X.shape)
print("Labels:", y.shape)


# -------------------- INITIALIZE RANDOM FOREST ---------------

model = RandomForestClassifier(
    n_estimators=300,
    random_state=42,
    class_weight="balanced"
)


# -------------------- DEFINE CROSS-VALIDATION ----------------

# Split the dataset into five stratified folds.
# Each fold maintains approximately the same class proportions.

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)


# -------------------- DEFINE EVALUATION METRICS --------------

scoring_metrics = {
    "accuracy": "accuracy",
    "precision": "precision_macro",
    "recall": "recall_macro",
    "f1": "f1_macro"
}


# -------------------- PERFORM CROSS-VALIDATION ---------------

print("\nRunning 5-fold cross-validation...")

results = cross_validate(
    model,
    X,
    y,
    cv=cv,
    scoring=scoring_metrics,
    return_train_score=True,
    n_jobs=-1
)


# -------------------- DISPLAY RESULTS ------------------------

print("\n==========================================")
print("5-FOLD CROSS-VALIDATION RESULTS")
print("==========================================")

for metric in scoring_metrics:

    scores = results[f"test_{metric}"]

    print(f"\n{metric.upper()}")
    print("Individual fold scores:", scores)
    print(f"Mean: {scores.mean():.4f}")
    print(f"Standard deviation: {scores.std():.4f}")


# -------------------- SAVE FOLD RESULTS ----------------------

fold_results = pd.DataFrame({
    "fold": range(1, 6),
    "accuracy": results["test_accuracy"],
    "precision": results["test_precision"],
    "recall": results["test_recall"],
    "f1_score": results["test_f1"]
})

fold_results.to_csv(
    OUTPUT_DIR / "fold_results.csv",
    index=False
)


# -------------------- SAVE SUMMARY ---------------------------

summary = fold_results.drop(columns="fold").agg(
    ["mean", "std"]
)

summary.to_csv(
    OUTPUT_DIR / "cross_validation_summary.csv"
)


# -------------------- FINAL MESSAGE --------------------------

print("\nResults saved in:")
print(OUTPUT_DIR)

print("\nCross-validation completed!")
