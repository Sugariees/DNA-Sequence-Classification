# ============================================================
# PROJECT: Machine Learning-Based DNA Sequence Classification
# SCRIPT: Final Experiment Comparison
#
# PURPOSE:
# Consolidate results from all four completed experiments.
# No model training is performed here.
# ============================================================

from pathlib import Path
import pandas as pd


# -------------------- PROJECT PATHS --------------------------

PROJECT_DIR = Path(__file__).resolve().parent.parent

RESULTS_DIR = PROJECT_DIR / "results" / "final"

OUTPUT_DIR = RESULTS_DIR / "final_comparison"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# -------------------- LOAD CLASSIFICATION REPORTS ------------

def load_classification_metrics(report_path):
    """Extract accuracy and macro-average metrics."""

    report = pd.read_csv(report_path, index_col=0)

    macro = report.loc["macro avg"]

    accuracy = float(report.loc["accuracy"].iloc[0])

    return {
        "accuracy": accuracy,
        "precision": float(macro["precision"]),
        "recall": float(macro["recall"]),
        "f1_score": float(macro["f1-score"])
    }


# -------------------- EXPERIMENT 1 ---------------------------

original = load_classification_metrics(
    RESULTS_DIR / "reports" /
    "random_forest_classification_report.csv"
)


# -------------------- EXPERIMENT 2 ---------------------------

representative = load_classification_metrics(
    RESULTS_DIR / "representative_model" /
    "classification_report.csv"
)


# -------------------- EXPERIMENT 3 ---------------------------

# Cross-validation summary stores mean and standard deviation
# in rows, with evaluation metrics in columns.

cv_file = (
    RESULTS_DIR
    / "cross_validation"
    / "cross_validation_summary.csv"
)

cv_data = pd.read_csv(cv_file, index_col=0)

cv_mean = cv_data.loc["mean"]


# -------------------- EXPERIMENT 4 ---------------------------

singleton = load_classification_metrics(
    RESULTS_DIR / "singleton_evaluation" /
    "classification_report.csv"
)


# -------------------- BUILD COMPARISON TABLE -----------------

experiments = [
    {
        "experiment": "Original Random Holdout",
        "dataset_size": 150,
        "training_size": 120,
        "testing_size": 30,
        "evaluation_method": "Stratified 80/20 random holdout",
        **original
    },
    {
        "experiment": "Representative Dataset",
        "dataset_size": 37,
        "training_size": 29,
        "testing_size": 8,
        "evaluation_method": "Stratified 80/20 random holdout",
        **representative
    },
    {
        "experiment": "Five-Fold Cross Validation",
        "dataset_size": 150,
        "training_size": None,
        "testing_size": None,
        "evaluation_method": "Stratified 5-fold cross-validation",
        "accuracy": float(cv_mean["accuracy"]),
        "precision": float(cv_mean["precision"]),
        "recall": float(cv_mean["recall"]),
        "f1_score": float(cv_mean["f1_score"])
    },
    {
        "experiment": "Singleton Evaluation",
        "dataset_size": 150,
        "training_size": 132,
        "testing_size": 18,
        "evaluation_method": "Large clusters vs singleton groups",
        **singleton
    }
]


comparison_df = pd.DataFrame(experiments)

comparison_df = comparison_df.rename(columns={
    "precision": "macro_precision",
    "recall": "macro_recall",
    "f1_score": "macro_f1_score"
})


# -------------------- SAVE COMPARISON ------------------------

comparison_file = (
    OUTPUT_DIR / "final_experiment_comparison.csv"
)

comparison_df.to_csv(comparison_file, index=False)


# -------------------- SAVE CROSS-VALIDATION DETAILS ----------

cv_data.to_csv(
    OUTPUT_DIR / "cross_validation_details.csv"
)


# -------------------- DISPLAY RESULTS ------------------------

print("\nFINAL EXPERIMENT COMPARISON")
print("=" * 90)

print(comparison_df.to_string(index=False))

print("\nCROSS-VALIDATION MEAN AND STANDARD DEVIATION")
print("=" * 60)

print(cv_data.to_string())

print("\nReports saved successfully:")
print(comparison_file)
print(OUTPUT_DIR / "cross_validation_details.csv")

print("\nExperiment comparison completed!")


# ============================================================
# END OF SCRIPT
# ============================================================
