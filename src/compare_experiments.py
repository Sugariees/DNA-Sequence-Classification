# ============================================================
# PROJECT: DNA Sequence Classification
# SCRIPT: Compare Experiments
#
# PURPOSE:
# Compare the baseline Random Forest model (150 sequences)
# with the representative Random Forest model (37 sequences).
#
# INPUT:
# Experiment 1 and Experiment 2 classification reports.
#
# OUTPUT:
# Overall comparison and per-species comparison CSV files.
# ============================================================


# -------------------- IMPORT LIBRARIES ----------------------

from pathlib import Path
import pandas as pd


# -------------------- PROJECT PATHS --------------------------

PROJECT_DIR = Path(__file__).resolve().parent.parent

# Experiment 1: Original 150-sequence model
EXPERIMENT_1_REPORT = (
    PROJECT_DIR
    / "results"
    / "final"
    / "reports"
    / "random_forest_classification_report.csv"
)

# Experiment 2: Representative 37-sequence model
EXPERIMENT_2_REPORT = (
    PROJECT_DIR
    / "results"
    / "final"
    / "representative_model"
    / "classification_report.csv"
)

OUTPUT_DIR = (
    PROJECT_DIR
    / "results"
    / "final"
    / "experiment_comparison"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# -------------------- EXPERIMENT DETAILS --------------------

EXPERIMENT_DETAILS = {
    "Experiment 1 - Original": {
        "total_sequences": 150,
        "training_sequences": 120,
        "testing_sequences": 30
    },

    "Experiment 2 - Representatives": {
        "total_sequences": 37,
        "training_sequences": 29,
        "testing_sequences": 8
    }
}


# -------------------- LOAD CLASSIFICATION REPORTS ------------

def load_report(report_path):

    if not report_path.exists():
        raise FileNotFoundError(
            f"Classification report not found: {report_path}"
        )

    # Load the saved classification report
    report = pd.read_csv(report_path, index_col=0)

    return report


report_1 = load_report(EXPERIMENT_1_REPORT)
report_2 = load_report(EXPERIMENT_2_REPORT)


# -------------------- EXTRACT OVERALL METRICS ----------------

def extract_overall_metrics(report):

    # Extract overall accuracy
    accuracy = report.loc["accuracy", "f1-score"]

    # Extract macro-average metrics
    macro_precision = report.loc["macro avg", "precision"]
    macro_recall = report.loc["macro avg", "recall"]
    macro_f1 = report.loc["macro avg", "f1-score"]

    # Extract weighted-average metrics
    weighted_precision = report.loc["weighted avg", "precision"]
    weighted_recall = report.loc["weighted avg", "recall"]
    weighted_f1 = report.loc["weighted avg", "f1-score"]

    return {
        "accuracy": accuracy,
        "macro_precision": macro_precision,
        "macro_recall": macro_recall,
        "macro_f1": macro_f1,
        "weighted_precision": weighted_precision,
        "weighted_recall": weighted_recall,
        "weighted_f1": weighted_f1
    }


metrics_1 = extract_overall_metrics(report_1)
metrics_2 = extract_overall_metrics(report_2)


# -------------------- CREATE OVERALL COMPARISON --------------

comparison_rows = []

for experiment_name, metrics in [
    ("Experiment 1 - Original", metrics_1),
    ("Experiment 2 - Representatives", metrics_2)
]:

    details = EXPERIMENT_DETAILS[experiment_name]

    row = {
        "experiment": experiment_name,
        **details,
        **metrics
    }

    comparison_rows.append(row)


comparison_df = pd.DataFrame(comparison_rows)


# -------------------- CREATE PER-SPECIES COMPARISON ----------

species_names = [
    "Bacillus subtilis",
    "Escherichia coli",
    "Staphylococcus aureus"
]

species_rows = []

for species in species_names:

    # Extract metrics for each species from both experiments
    for experiment_name, report in [
        ("Experiment 1 - Original", report_1),
        ("Experiment 2 - Representatives", report_2)
    ]:

        species_rows.append({
            "experiment": experiment_name,
            "species": species,
            "precision": report.loc[species, "precision"],
            "recall": report.loc[species, "recall"],
            "f1_score": report.loc[species, "f1-score"],
            "support": report.loc[species, "support"]
        })


species_comparison_df = pd.DataFrame(species_rows)


# -------------------- SAVE COMPARISON REPORTS ----------------

comparison_df.to_csv(
    OUTPUT_DIR / "overall_experiment_comparison.csv",
    index=False
)

species_comparison_df.to_csv(
    OUTPUT_DIR / "species_wise_comparison.csv",
    index=False
)


# -------------------- DISPLAY RESULTS ------------------------

print("\n==========================================")
print("OVERALL EXPERIMENT COMPARISON")
print("==========================================")

print(comparison_df.to_string(index=False))

print("\n==========================================")
print("SPECIES-WISE COMPARISON")
print("==========================================")

print(species_comparison_df.to_string(index=False))

print("\n==========================================")
print("COMPARISON COMPLETED")
print("==========================================")

print(f"\nReports saved in:\n{OUTPUT_DIR}")
