# Composition Model Comparison

Experimental side-by-side comparison of classifiers for the `stop_traversal` (clade composition) prediction task. Loads cached training data, replicates production preprocessing, trains and evaluates multiple classifiers.

## Task

Predict whether traversal should stop at a given clade node based on statistical and taxonomic features. The production `CompositionModeller` uses one model — this script compares 6+ alternatives with uniform preprocessing.

## Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--output_dir` | No | `composition_comparison_outputs` | Output directory for results |
| `--input_cache` | No | `deployment/ml_api/training_results_cache.parquet` | Path to training data cache |

## Example

```bash
python deployment/model_evaluation/analysis_scripts/composition_model_comparison.py \
    --output_dir /path/to/composition_comparison_outputs \
    --input_cache deployment/ml_api/training_results_cache.parquet
```

## Methods

### Dataset

Training data comes from `training_results_cache.parquet`, produced by the evaluation pipeline's clade-traversal step (`data_set_traversal_with_precision()`). The dataset comprises 4,920 rows and 23 columns: 11 family-level viral taxonomy proportions, 4 statistical features (`n_leaves`, `tax_diversity`, `Min_Dist`, `Min_Shared`), 7 metadata columns, and the binary target `stop_traversal`. The target distribution is 65.3% positive (stop) and 34.7% negative (continue).

### Preprocessing

Metadata columns are dropped. The remaining 15 features are split into two groups: **stat features** (4 numeric columns) and **taxonomy features** (11 composition columns). A stratified 80/20 train-test split is applied (`random_state=42`). Stat features are standardised via `StandardScaler`; taxonomy features (already [0,1] proportions) are left unaltered.

### Models Compared

| Model | Features | Key Parameters |
|---|---|---|
| FixedThreshold(min_shared >= 0.6) | All | Rule-based baseline, no training |
| XGBoost+Optuna | All | 50-trial Optuna search over n_estimators[200-800], max_depth[3-8], lr[0.01-0.2], subsample[0.6-1.0], colsample[0.6-1.0], reg_alpha/lambda[1e-3-10], gamma[0-5], min_child_weight[1-10], pos_weight |
| XGBoost | All | 300 trees, depth=6, lr=0.1, subsample=0.8, colsample=0.8 |
| RandomForest | All | 300 trees, depth=12, min_samples_leaf=3, balanced class_weight |
| GradientBoosting | All | 300 trees, depth=5, lr=0.1, subsample=0.8 |
| LogisticRegression | Stats only (4) | L2-regularised, balanced class_weight |
| LightGBM (optional) | All | 300 trees, depth=6, lr=0.1, balanced class_weight |

### Evaluation Metrics

Accuracy, precision, recall, F1, ROC-AUC, PR-AUC (positive class = stop_traversal). ROC and PR curves are generated for each model offering predicted probabilities.

## Outputs

### Tables

#### `comparison_results.csv`

One row per model, sorted by F1 descending.

| Column | Description |
|---|---|
| `Model` | Model name |
| `Accuracy` | `(TP+TN) / (TP+TN+FP+FN)` |
| `Precision` | `TP / (TP+FP)`, positive class = stop_traversal |
| `Recall` | `TP / (TP+FN)`, positive class = stop_traversal |
| `F1` | `2 × Precision × Recall / (Precision + Recall)` |
| `ROC-AUC` | Area under the ROC curve computed from predicted probabilities |
| `PR-AUC` | Average precision score (area under Precision-Recall curve) |
| `Train_Time_s` | Wall-clock training time in seconds |

#### `classification_report.txt`

Full sklearn `classification_report` for each model, including per-class precision/recall/F1, support counts, ROC-AUC, PR-AUC, and training time.

### Plots

#### `roc_curves.png`

All models overlaid on a single axis. X-axis: False Positive Rate, Y-axis: True Positive Rate. Each curve labelled with model name and AUC. Black dashed diagonal = random classifier. Legend in lower-right.

#### `pr_curves.png`

All models overlaid on a single axis. X-axis: Recall, Y-axis: Precision. Each curve labelled with model name and PR-AUC. Grey dashed horizontal line at the class prevalence ratio (~0.35) = no-skill baseline.

#### `f1_comparison.png`

Horizontal bar chart. Models on the Y-axis, F1 score on X-axis [0, 1]. The best-performing model is highlighted in green; others in grey. Value label at the end of each bar.

#### `precision_recall_bars.png`

Grouped bar chart with models on the X-axis. Three bars per model: Precision (blue), Recall (red), F1 (green). X-axis labels rotated 30°.

#### `confusion_matrices.png`

Grid of normalised confusion matrices (one per model, max 3 columns). Rows: true labels (Continue, Stop). Columns: predicted labels. Cell values are row-normalised fractions. Colormap: Blues.

#### `feature_importance.png`

Top-15 horizontal bar charts for each tree-based model (XGBoost, RandomForest, GradientBoosting, LightGBM). Models without feature importance (FixedThreshold, LogisticRegression) are omitted.

#### `training_time.png`

Horizontal bar chart. X-axis: training time in seconds. The slowest model is highlighted in orange; others in grey. Time value label at the end of each bar.
