# Sample Analysis

Analyzes clustering output from map_clustering pipeline using trained models to predict optimal clades and recall thresholds for new metagenomic samples.

## Overview

The `analyze_samples.py` script processes sample directories containing clustering results and uses pre-trained:

- **CompositionModeller**: Predicts optimal clade selection
- **RecallModeller**: Predicts recall cutoff thresholds

## Installation

```bash
# Activate the virtual environment
source .venv/bin/activate

# Set Python path
export PYTHONPATH=$(pwd)
```

## Usage

```bash
python analyze_samples.py \
    --samples-dir /path/to/samples \
    --training-dir /path/to/training \
    --output-dir analysis_output \
    --tax-level order \
    --generate-plots
```

### Arguments

| Argument           | Required | Default                            | Description                                |
| ------------------ | -------- | ---------------------------------- | ------------------------------------------ |
| `--samples-dir`    | Yes      | -                                  | Directory containing sample subdirectories |
| `--training-dir`   | Yes      | -                                  | Directory containing trained models        |
| `--output-dir`     | No       | `analysis_output`                  | Output directory                           |
| `--tax-level`      | No       | `order`                            | Taxonomic level for analysis               |
| `--generate-plots` | No       | False                              | Generate visualization plots               |
| `--taxids-file`    | No       | `{training-dir}/taxids_to_use.tsv` | Path to taxids file                        |
| `--taxonomy-db`    | No       | `taxa.db`                          | Path to taxonomy database                  |

### Input Structure

Expected directory structure:

```
samples_dir/
├── ERR13488565/
│   └── output_clustering/
│       └── {run_id}/
│           ├── classification/
│           │   └── {sample}_merged_classification.tsv
│           ├── output/
│           │   └── matched_assemblies.tsv
│           └── clustering/
│               └── ...
└── ERR13488567/
    └── ...
```

### Training Directory Structure

```
training_dir/
├── models/
│   ├── recall_modeller.joblib
│   ├── composition_modeller.joblib
│   └── crosshit_modeller.joblib
└── taxids_to_use.tsv
```

## Output Structure

```
output_dir/
├── summary.tsv              # Per-sample summary statistics
├── all_predictions.csv     # All predictions concatenated
├── metadata.json           # Run metadata
├── plots/                  # Visualizations (if --generate-plots)
│   ├── precision_distribution.png
│   ├── sample_comparison.png
│   ├── taxa_frequency.png
│   └── clade_size_distribution.png
└── samples/
    ├── ERR13488565/
    │   └── predictions.tsv
    └── ERR13488567/
        └── predictions.tsv
```

### Summary Columns

| Column                  | Description                                          |
| ----------------------- | ---------------------------------------------------- |
| `sample`                | Sample identifier                                    |
| `n_detections`          | Number of detections                                 |
| `mean_precision`        | Mean precision across nodes                          |
| `median_precision`      | Median precision                                     |
| `unique_taxa`           | Number of unique taxids detected                     |
| `n_clades`              | Total number of clades                               |
| `high_confidence_count` | Number of high-confidence detections (precision=1.0) |
| `low_confidence_count`  | Number of low-confidence detections (precision<0.5)  |
| `mean_n_leaves`         | Mean number of leaves per node                       |

### Predictions Columns

| Column             | Description              |
| ------------------ | ------------------------ |
| `data_set`         | Sample identifier        |
| `node`             | Node identifier          |
| `n_leaves`         | Number of leaves in node |
| `leaves`           | List of leaf identifiers |
| `best_taxid_match` | Best matching taxid      |
| `node_precision`   | Predicted precision      |
| `node_taxids`      | Taxids in node           |
| `leaf`             | Representative leaf      |
| `best_match_taxid` | Best match taxid         |
| `description`      | Taxonomic description    |

## Example

```bash
# Basic analysis
python deployment/analysis/analyze_samples.py \
    --samples-dir /data/samples \
    --training-dir /models/virus_training \
    --output-dir results

# With plots
python deployment/analysis/analyze_samples.py \
    --samples-dir /data/samples \
    --training-dir /models/virus_training \
    --output-dir results \
    --tax-level order \
    --generate-plots
```

## Dependencies

- Python 3.10+
- metagenomics_utils (in project root)
- See `.venv` for Python package dependencies
