# ============================================================
# PROJECT: DNA Sequence Classification Using Machine Learning
# MODEL: Random Forest Classifier
# PURPOSE: Classify bacterial 16S rRNA sequences into
#          Escherichia coli, Bacillus subtilis, and
#          Staphylococcus aureus.
# ============================================================

# SECTION 1: Import required libraries

import pandas as pd

from pathlib import Path

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier

from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix
)


# SECTION 2: Define input and output file paths

dataset_file = Path("data/processed/final/ml_dataset.csv")

results_directory = Path("results/final/reports")

results_directory.mkdir(parents=True, exist_ok=True)


# SECTION 3: Load the prepared machine-learning dataset

print("\nLOADING DATASET")
print("=" * 60)

dataset = pd.read_csv(dataset_file)

print("Dataset loaded successfully.")
print("Dataset shape:", dataset.shape)


# SECTION 4: Separate features (X) and labels (y)

# The accession column identifies each sequence.
# The species column contains the class we want to predict.
# All remaining numerical columns are k-mer features.

X = dataset.drop(columns=["accession", "species"])

y = dataset["species"]

print("\nFEATURE AND LABEL INFORMATION")
print("=" * 60)

print("Number of sequences:", len(X))
print("Number of k-mer features:", X.shape[1])

print("\nClass distribution:")
print(y.value_counts())


# SECTION 5: Divide data into training and testing sets

# 80% of the data is used for training.
# 20% is reserved for testing.
# Stratification preserves the class proportions.

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    random_state=42,
    stratify=y
)

print("\nDATA SPLIT")
print("=" * 60)

print("Training sequences:", len(X_train))
print("Testing sequences:", len(X_test))


# SECTION 6: Create the Random Forest model

# n_estimators = 300 decision trees.
# random_state = 42 makes the experiment reproducible.
# class_weight = balanced adjusts for class imbalance if present.

model = RandomForestClassifier(
    n_estimators=300,
    random_state=42,
    class_weight="balanced"
)


# SECTION 7: Train the model

print("\nTRAINING RANDOM FOREST MODEL")
print("=" * 60)

model.fit(X_train, y_train)

print("Model training completed.")


# SECTION 8: Predict the species of test sequences

predictions = model.predict(X_test)


# SECTION 9: Evaluate model performance

accuracy = accuracy_score(y_test, predictions)

print("\nMODEL EVALUATION")
print("=" * 60)

print(f"Accuracy: {accuracy * 100:.2f}%")

print("\nClassification Report:")
print(classification_report(y_test, predictions, zero_division=0))

print("\nConfusion Matrix:")

labels = sorted(y.unique())

matrix = confusion_matrix(
    y_test,
    predictions,
    labels=labels
)

confusion_table = pd.DataFrame(
    matrix,
    index=labels,
    columns=labels
)

print(confusion_table)


# SECTION 10: Save evaluation results

confusion_table.to_csv(
    results_directory / "random_forest_confusion_matrix.csv"
)

report = classification_report(
    y_test,
    predictions,
    zero_division=0,
    output_dict=True
)

pd.DataFrame(report).transpose().to_csv(
    results_directory / "random_forest_classification_report.csv"
)


# SECTION 11: Save predictions

prediction_results = pd.DataFrame({
    "actual_species": y_test.values,
    "predicted_species": predictions
})

prediction_results.to_csv(
    results_directory / "random_forest_predictions.csv",
    index=False
)
# SECTION 11A: Save the trained Random Forest model

import joblib

joblib.dump(
    model,
    "results/final/random_forest_model.pkl"
)

print("Trained model saved successfully.")

# SECTION 12: Display completion message

print("\nRESULTS SAVED")
print("=" * 60)

print("1. random_forest_confusion_matrix.csv")
print("2. random_forest_classification_report.csv")
print("3. random_forest_predictions.csv")

print("\nBASELINE MODEL TRAINING COMPLETED!")
