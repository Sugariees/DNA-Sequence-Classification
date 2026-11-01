# ============================================================
# PROJECT: DNA Sequence Classification Using Machine Learning
# PURPOSE: Visualise Random Forest model performance
# ============================================================

# SECTION 1: Import libraries

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
from sklearn.metrics import confusion_matrix


# SECTION 2: Define file paths

results_directory = Path("results/final")

reports_directory = results_directory / "reports"

confusion_file = reports_directory / "random_forest_confusion_matrix.csv"

features_file = Path("data/processed/final/kmer_features.csv")

model_output = results_directory / "random_forest_model.pkl"


# SECTION 3: Load the confusion matrix

print("\nLOADING CONFUSION MATRIX")

confusion_data = pd.read_csv(confusion_file, index_col=0)

plt.figure(figsize=(9, 7))

sns.heatmap(
    confusion_data,
    annot=True,
    fmt="d",
    cmap="Blues",
    linewidths=0.5
)

plt.title("Random Forest - Confusion Matrix")
plt.xlabel("Predicted Species")
plt.ylabel("Actual Species")
plt.tight_layout()

confusion_plot = results_directory / "confusion_matrix.png"

plt.savefig(confusion_plot, dpi=300)
plt.close()

print("Confusion matrix saved:", confusion_plot)


# SECTION 4: Load the trained Random Forest model

import joblib

if not model_output.exists():
    raise FileNotFoundError(
        "Trained model not found. Update the training script "
        "to save the model before running this visualisation."
    )

model = joblib.load(model_output)


# SECTION 5: Load k-mer feature names

features_data = pd.read_csv(features_file)

# Exclude accession because it is an identifier, not a feature.

feature_names = [
    column for column in features_data.columns
    if column != "accession"
]


# SECTION 6: Extract feature importance

importance_data = pd.DataFrame({
    "kmer": feature_names,
    "importance": model.feature_importances_
})

importance_data = importance_data.sort_values(
    by="importance",
    ascending=False
)

importance_data.to_csv(
    reports_directory / "kmer_feature_importance.csv",
    index=False
)


# SECTION 7: Plot the top 20 important k-mers

top_features = importance_data.head(20)

plt.figure(figsize=(12, 8))

sns.barplot(
    data=top_features,
    x="importance",
    y="kmer"
)

plt.title("Top 20 Important K-mers - Random Forest")
plt.xlabel("Feature Importance")
plt.ylabel("DNA 4-mer")

plt.tight_layout()

importance_plot = results_directory / "top_20_kmer_importance.png"

plt.savefig(importance_plot, dpi=300)
plt.close()


# SECTION 8: Completion message

print("\nVISUALISATION COMPLETED")
print("=" * 60)

print("1. Confusion matrix:", confusion_plot)
print("2. Feature importance:", importance_plot)
print("3. Feature importance data:",
      reports_directory / "kmer_feature_importance.csv")
