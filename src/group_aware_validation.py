# ============================================================
# PROJECT: DNA Sequence Classification
# SCRIPT: Group-Aware Cross-Validation
#
# PURPOSE:
# Evaluate the Random Forest model while preventing
# sequences from the same similarity group from appearing
# in both training and testing folds.
# ============================================================


# -------------------- IMPORT LIBRARIES ----------------------

from pathlib import Path

import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score
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
    / "group_aware_validation"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# -------------------- LOAD DATA ------------------------------

dataset = pd.read_csv(DATA_FILE)

groups_data = pd.read_csv(GROUP_FILE)

print("ML dataset:", dataset.shape)
print("Group dataset:", groups_data.shape)


# -------------------- VERIFY ACCESSIONS ----------------------

# Match the group information to the ML dataset.
# Validate that every accession has exactly one group assignment.

if dataset["accession"].duplicated().any():
    raise ValueError("Duplicate accessions in ML dataset.")

if groups_data["accession"].duplicated().any():
    raise ValueError("Duplicate accessions in group dataset.")

if len(dataset) != len(groups_data):
    raise ValueError("Dataset and group file have different record counts.")

merged = dataset.merge(
    groups_data,
    on=["accession", "species"],
    how="left",
    validate="one_to_one"
)

if merged["group_id"].isna().any():
    raise ValueError("Some sequences do not have group assignments.")

print("\nAll accessions successfully matched.")


# -------------------- PREPARE FEATURES AND LABELS ------------

# Remove identifiers and labels from the feature matrix.

X = merged.drop(
    columns=["accession", "species", "group_id"]
)

y = merged["species"]

# Group IDs must be globally unique across species.
# The existing group IDs already include species names.

groups = merged["group_id"]

print("\nFeature matrix:", X.shape)
print("Number of sequences:", len(y))
print("Number of similarity groups:", groups.nunique())

print("\nSpecies distribution:")
print(y.value_counts())

print("\nGroups per species:")
print(
    merged.groupby("species")["group_id"].nunique()
)


# -------------------- INITIALIZE MODEL -----------------------

model = RandomForestClassifier(
    n_estimators=300,
    random_state=42,
    class_weight="balanced"
)


# -------------------- DEFINE GROUP-AWARE VALIDATION ----------

# Try five folds while keeping each similarity group intact.
# The number of groups per species must be sufficient.

cv = StratifiedGroupKFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)


# -------------------- RUN VALIDATION -------------------------

fold_results = []
prediction_results = []

print("\nStarting group-aware cross-validation...")

for fold_number, (train_indices, test_indices) in enumerate(
    cv.split(X, y, groups),
    start=1
):

    # Select training and testing records
    X_train = X.iloc[train_indices]
    X_test = X.iloc[test_indices]

    y_train = y.iloc[train_indices]
    y_test = y.iloc[test_indices]

    train_groups = set(groups.iloc[train_indices])
    test_groups = set(groups.iloc[test_indices])

    # Confirm that no similarity group is shared
    overlap = train_groups.intersection(test_groups)

    if overlap:
        raise ValueError(
            f"Group leakage detected in fold {fold_number}: {overlap}"
        )

    # Train a fresh model for this fold
    model.fit(X_train, y_train)

    # Predict species for the held-out sequences
    predictions = model.predict(X_test)

    # Calculate evaluation metrics
    fold_accuracy = accuracy_score(y_test, predictions)

    fold_precision = precision_score(
        y_test,
        predictions,
        average="macro",
        zero_division=0
    )

    fold_recall = recall_score(
        y_test,
        predictions,
        average="macro",
        zero_division=0
    )

    fold_f1 = f1_score(
        y_test,
        predictions,
        average="macro",
        zero_division=0
    )

    fold_results.append({
        "fold": fold_number,
        "training_sequences": len(train_indices),
        "testing_sequences": len(test_indices),
        "training_groups": len(train_groups),
        "testing_groups": len(test_groups),
        "accuracy": fold_accuracy,
        "precision": fold_precision,
        "recall": fold_recall,
        "f1_score": fold_f1
    })

    # Save individual predictions for later inspection
    for accession, actual, predicted in zip(
        merged.iloc[test_indices]["accession"],
        y_test,
        predictions
    ):

        prediction_results.append({
            "fold": fold_number,
            "accession": accession,
            "actual_species": actual,
            "predicted_species": predicted,
            "correct": actual == predicted
        })

    print(
        f"Fold {fold_number}: "
        f"Train={len(train_indices)}, "
        f"Test={len(test_indices)}, "
        f"Accuracy={fold_accuracy:.4f}, "
        f"Groups separated successfully"
    )


# -------------------- CREATE RESULTS TABLE -------------------

results_df = pd.DataFrame(fold_results)

predictions_df = pd.DataFrame(prediction_results)


# -------------------- CALCULATE SUMMARY ----------------------

summary_df = results_df[
    ["accuracy", "precision", "recall", "f1_score"]
].agg(["mean", "std"])


# -------------------- SAVE RESULTS ---------------------------

results_df.to_csv(
    OUTPUT_DIR / "group_aware_fold_results.csv",
    index=False
)

predictions_df.to_csv(
    OUTPUT_DIR / "group_aware_predictions.csv",
    index=False
)

summary_df.to_csv(
    OUTPUT_DIR / "group_aware_summary.csv"
)


# -------------------- DISPLAY RESULTS ------------------------

print("\n==========================================")
print("GROUP-AWARE CROSS-VALIDATION RESULTS")
print("==========================================")

print("\nFold results:")
print(results_df.to_string(index=False))

print("\nMean and standard deviation:")
print(summary_df)

print("\nAll folds completed without group overlap.")

print("\nResults saved in:")
print(OUTPUT_DIR)
