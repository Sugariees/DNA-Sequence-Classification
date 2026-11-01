# Machine Learning-Based Classification of Bacterial 16S rRNA Gene Sequences

## 1. Project Overview

This project investigates the use of machine learning to classify bacterial 16S rRNA gene sequences into three bacterial species:

- *Escherichia coli*
- *Bacillus subtilis*
- *Staphylococcus aureus*

DNA sequences were retrieved programmatically from the NCBI Nucleotide database using Biopython. Sequence processing, k-mer feature extraction, similarity analysis and Random Forest classification were performed using Python.

## 2. Research Question

Can machine learning distinguish between three bacterial species using their 16S rRNA gene sequence patterns?

## 3. Dataset

The dataset contains 150 curated 16S rRNA sequences.

| Species | Number of sequences |
|---|---:|
| *Escherichia coli* | 50 |
| *Bacillus subtilis* | 50 |
| *Staphylococcus aureus* | 50 |
| **Total** | **150** |

Sequences were retrieved from the NCBI Nucleotide database and filtered using the following criteria:

- Annotated 16S rRNA gene feature.
- Sequence length between 1,300 and 1,600 bp.
- Only A, T, G and C nucleotides.
- Removal of exact duplicate sequences within each species.

**Dataset considerations:**

- Sequence length filtering does not guarantee that every sequence represents a complete, full-length 16S rRNA gene.
- The dataset is a curated research dataset, not an exhaustive collection of these bacterial species.
- The dataset represents only three selected bacterial species.

## 4. Methodology

The project follows these steps:

1. Retrieve nucleotide sequences from NCBI using Biopython Entrez.
2. Extract annotated 16S rRNA gene sequences.
3. Perform sequence quality control.
4. Identify exact duplicates and examine sequence similarity.
5. Perform multiple sequence alignment using MAFFT.
6. Cluster similar sequences using CD-HIT-EST.
7. Extract normalized tetranucleotide (4-mer) frequencies.
8. Train a Random Forest classifier.
9. Evaluate model performance using classification metrics and confusion matrices.
10. Investigate sequence similarity and potential data leakage using different evaluation strategies.

## 5. Feature Extraction

Each DNA sequence was converted into a numerical feature vector using tetranucleotide frequencies.

A 4-mer consists of four consecutive nucleotides.

Since DNA has four possible nucleotides:

**4⁴ = 256 possible tetranucleotides**

Each sequence is represented by 256 normalized frequency features.

The frequencies are calculated using overlapping four-nucleotide windows and normalized by the total number of valid 4-mer windows.

## 6. Machine Learning Model

A Random Forest classifier was implemented using scikit-learn.

Main model settings:

- Number of trees: 300
- Random state: 42
- Class weighting: Balanced
- Input features: 256 normalized 4-mer frequencies
- Classification task: Three-class bacterial species classification

## 7. Experimental Results

Four main evaluation experiments were performed.

| Experiment | Dataset size | Test sequences | Accuracy |
|---|---:|---:|---:|
| Original random holdout | 150 | 30 | 100% |
| Representative dataset | 37 | 8 | 87.50% |
| Five-fold cross-validation | 150 | 5 folds | 99.33% mean |
| Large-cluster vs singleton evaluation | 150 | 18 | 94.44% |

### 7.1 Original Random Holdout

The original dataset contained 150 sequences. An 80:20 stratified random split was used for training and testing.

The model achieved 100% accuracy on the 30-sequence test set.

### 7.2 Representative Dataset Evaluation

CD-HIT-EST clustering at 99% sequence identity reduced the original dataset from 150 sequences to 37 representative sequences:

- *Escherichia coli*: 13 representatives
- *Bacillus subtilis*: 12 representatives
- *Staphylococcus aureus*: 12 representatives

The representative dataset model achieved 87.50% accuracy on an eight-sequence test set.

This test set is small, so the result should be interpreted cautiously.

### 7.3 Five-Fold Cross-Validation

Stratified five-fold cross-validation was performed on the original 150-sequence dataset.

The mean accuracy was 99.33%.

However, ordinary stratified cross-validation does not prevent sequence similarity leakage because highly similar sequences may occur in different folds.

### 7.4 Large-Cluster vs Singleton Evaluation

To investigate the effect of sequence similarity, an additional experiment trained the model on 132 sequences belonging to multi-sequence similarity groups and evaluated it on 18 singleton sequences.

The model achieved:

- Accuracy: 94.44%
- Correct predictions: 17 out of 18
- Incorrect predictions: 1

Singleton status at a 99% sequence identity threshold does not establish biological independence. Therefore, this experiment is exploratory and should not be interpreted as independent external validation.

### Important Interpretation

The high accuracy obtained from random holdout and ordinary cross-validation may be optimistic because similar sequences can occur in both training and testing sets.

The representative dataset and singleton evaluation provide additional exploratory assessments, but neither establishes generalisation to independently collected bacterial strains.

## 8. Tools and Technologies

- Python
- Biopython
- NCBI Entrez
- pandas
- NumPy
- scikit-learn
- Matplotlib
- Seaborn
- MAFFT
- CD-HIT-EST
- Git and GitHub

## 9. Project Structure

```text
Project-01-DNA-Sequence-Classification/
│
├── data/
│   ├── raw/
│   └── processed/
│       ├── final/
│       └── diversity_controlled/
│
├── notebooks/
│
├── results/
│   └── final/
│
├── src/
│
├── requirements.txt
├── .gitignore
└── README.md
```

## 10. Limitations

- The dataset contains only three bacterial species.
- The dataset size is relatively small.
- Sequence length filtering does not guarantee complete, full-length 16S rRNA sequences.
- 16S rRNA is a conserved gene, and closely related sequences can share highly similar patterns.
- Random splitting can cause sequence similarity leakage.
- Ordinary stratified cross-validation does not eliminate this leakage.
- The representative dataset contains only 37 sequences.
- The singleton evaluation contains only 18 test sequences.
- Singleton sequences are not necessarily biologically independent.
- The experiments do not establish generalisation to unseen bacterial species or independently collected clinical isolates.

## 11. Future Scope

- Expand the dataset to additional bacterial species.
- Include independently collected strain-level sequences.
- Compare different k-mer sizes.
- Evaluate additional machine learning algorithms.
- Investigate more rigorous similarity-aware validation strategies.
- Explore deep learning approaches for DNA sequence classification.

## 12. Conclusion

This project demonstrates an end-to-end bioinformatics and machine learning workflow for classifying three bacterial species using 16S rRNA gene sequences.

It combines automated sequence retrieval from NCBI, biological sequence processing, k-mer feature extraction, sequence similarity analysis and supervised machine learning.

The experiments achieved promising results on the selected dataset. However, sequence redundancy, limited dataset size and validation constraints highlight the need for more rigorous evaluation using independently collected bacterial sequences.

## 13. Installation

### 13.1 Python Dependencies

Create a virtual environment:

```bash
python3 -m venv .venv
```

Activate the virtual environment:

```bash
source .venv/bin/activate
```

Install the required Python packages:

```bash
pip install -r requirements.txt
```

### 13.2 External Bioinformatics Tools

This project uses:

- **MAFFT** for multiple sequence alignment.
- **CD-HIT-EST** for nucleotide sequence clustering.

Install these tools on Ubuntu:

```bash
sudo apt update
sudo apt install mafft cd-hit
```

Verify the installations:

```bash
mafft --version
cd-hit-est -h
```
