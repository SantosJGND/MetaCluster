# Sample Analysis

Analyzes clustering output from map_clustering pipeline using trained models to predict optimal clades and recall thresholds for new metagenomic samples.

## Overview

The `analyze_samples.py` script processes sample directories containing clustering results and uses pre-trained:

- **CompositionModeller** (XGB, RF, GB, LR, or Optuna-XGB): Predicts optimal clade selection
- **RecallModeller** (GP+CLF, Direct XGBoost, XGBoost Multi-output, or Cutoff RF): Predicts recall cutoff thresholds

For each sample, the script:
1. Builds an overlap graph from clustering output
2. Uses the recall model to determine a percentile cutoff for filtering
3. Runs the composition model on both full and pruned (filtered) data
4. Outputs per-sample predictions and aggregate statistics

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

| Argument            | Required | Default                            | Description                                         |
| ------------------- | -------- | ---------------------------------- | --------------------------------------------------- |
| `--samples-dir`     | Yes      | -                                  | Directory containing sample subdirectories           |
| `--training-dir`    | Yes      | -                                  | Directory containing trained models (`models/` subdir) |
| `--output-dir`      | No       | `analysis_output`                  | Output directory                                     |
| `--tax-level`       | No       | `order`                            | Taxonomic level for analysis                         |
| `--generate-plots`  | No       | False                              | Generate visualization plots                         |
| `--taxids-file`     | No       | `{training-dir}/models/taxids_to_use.parquet` | Path to taxids file                 |
| `--taxonomy-db`     | No       | `taxa.db`                          | Path to taxonomy database (SQLite)                   |
| `--target-recall`   | No       | `1.0`                              | Target recall threshold for filtering                |
| `--data-set-divide` | No       | `20`                               | Number of recall divisions used in model training    |
| `--recall-model`    | No       | `direct_xgb`                      | Recall model variant: `gp_clf`, `direct`, `xgb`, `direct_xgb` |
| `--composition-model` | No     | `rf`                               | Composition model variant: `xgb`, `rf`, `gb`, `lr`, `xgb_optimized` |

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
│               └── (distance matrices, node data)
└── ERR13488567/
    └── ...
```

### Training Directory Structure

```
training_dir/
├── models/
│   ├── recall_gp_clf_pipeline.pkl    (GP+CLF recall model)
│   ├── recall_xgb_bundle.pkl         (XGBoost multi-output recall model)
│   ├── direct_xgb_bundle.pkl         (Direct XGBoost recall model)
│   ├── cutoff_recall_bundle.pkl      (Cutoff RF recall model)
│   ├── composition_xgb_bundle.pkl    (XGBoost composition model)
│   ├── composition_rf_bundle.pkl     (Random Forest composition model)
│   ├── composition_gb_bundle.pkl     (Gradient Boosting composition model)
│   ├── composition_lr_bundle.pkl     (Logistic Regression composition model)
│   ├── composition_optuna_bundle.pkl (Optuna-tuned XGBoost composition model)
│   └── taxids_to_use.parquet
```

One recall and one composition model are loaded; the script picks the first matching bundle found in the directory.

## Output Structure

```
output_dir/
├── summary.tsv                   # Per-sample summary statistics
├── all_predictions.csv           # All full-data predictions concatenated
├── all_predictions_pruned.csv    # All pruned predictions concatenated
├── cluster_statistics.tsv        # Per-sample cluster/filtering metrics
├── metadata.json                 # Run metadata
├── plots/                        # Visualizations (if --generate-plots)
│   ├── precision_distribution.png
│   ├── sample_comparison.png
│   ├── taxa_frequency.png
│   ├── clade_size_distribution.png
│   └── filtering_benefit.png
└── samples/
    ├── ERR13488565/
    │   ├── predictions.tsv           # Full-data clade predictions
    │   └── predictions_pruned.tsv    # Pruned clade predictions
    └── ERR13488567/
        ├── predictions.tsv
        └── predictions_pruned.tsv
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

### Cluster Statistics Columns

| Column                              | Description                                           |
| ----------------------------------- | ----------------------------------------------------- |
| `sample`                            | Sample identifier                                     |
| `n_total_intermediate_classifications` | Total intermediate classifications before pruning  |
| `n_classifications_with_coverage`   | Classifications with coverage > 0                     |
| `n_final_clusters`                  | Clades from full-data composition prediction          |
| `n_final_clusters_pruned`           | Clades from pruned composition prediction             |
| `n_indexes_kept_from_recall`        | Number of classifications kept after recall filtering |
| `target_percentile`                 | Predicted percentile cutoff                           |
| `target_recall`                     | Target recall threshold used                          |
| `prop_coverage_above_cutoff`        | Proportion of references with coverage retained       |
| `prop_coverage_below_cutoff`        | Proportion of references with coverage lost           |

## Example

```bash
# Basic analysis
python deployment/analysis/analyze_samples.py \
    --samples-dir /data/samples \
    --training-dir /models/virus_training \
    --output-dir results

# With plots and custom recall target
python deployment/analysis/analyze_samples.py \
    --samples-dir /data/samples \
    --training-dir /models/virus_training \
    --output-dir results \
    --tax-level order \
    --target-recall 0.95 \
    --generate-plots
```

## Dependencies

- Python 3.10+
- metagenomics_utils (in project root)
- See `.venv` for Python package dependencies
