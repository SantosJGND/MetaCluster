# Composition Prediction

Predicts the number of unique input taxids per dataset from m-stats features.
Employs both continuous regression (raw count) and binned ordinal classification
(very_low / low / medium / high / very_high). The predicted count is also used
as an additional feature for the last-TP division prediction pipeline.

## Overview

Two parallel tasks:

- **Regression (continuous)** — predicts the exact number of unique taxids
- **Classification (binned)** — predicts which ordinal bin the count falls into

### Bin Edges

| Bin | Label | Range (n_taxids) |
|---|---|---|
| 0 | `very_low` | 0–3 |
| 1 | `low` | 3–6 |
| 2 | `medium` | 6–9 |
| 3 | `high` | 9–12 |
| 4 | `very_high` | 12+ |

## Data Loading

Two modes, controlled by `--recall_cache`:

**1. From cache** (`--recall_cache` provided): reads pre-computed features from
`recall_results_cache.parquet` (6 stat columns + taxonomy proportions from the
intersection with `ref_taxa`). Targets are computed from the input summary TSVs
on disk.

**2. From study** (no cache): loads raw m_stats matrices per dataset via
`OverlapManager`, transforms with `RecallFeatureTransformer`, computes targets
from input summaries.

## Features

### Statistical Features (6)

| Feature | Description |
|---|---|
| `number_hits` | Total number of hits in the clustering output |
| `counts_kurtosis` | Kurtosis of read-count distribution across leaves |
| `counts_skewness` | Skewness of read-count distribution across leaves |
| `tax_diversity_shannon` | Shannon diversity index of taxonomic assignments |
| `max_uniq_reads` | Maximum unique reads assigned to a single leaf |
| `total_uniq_reads` | Total unique reads across all leaves |

### Taxonomy Features

Proportion columns at the configured `--tax_level` (default: family). When using
cache mode, these are the family-level proportions present in the cache
(e.g. Baculoviridae, Coronaviridae, Flaviviridae, etc.). In study mode, they
are derived from the `taxids_to_use` plan via `get_subset_composition_counts()`.

## Models

### Regression (continuous count)

| Model | Grid Search | Key Parameters |
|---|---|---|
| MeanBaseline | — | Predicts training-set mean |
| **XGBoost+Grid** | 243 combinations | n_est=[100,300,500], depth=[4,6,8], lr=[0.05,0.1,0.2], subsample=[0.7,0.85,1.0], colsample=[0.7,0.85,1.0] |
| XGBoost | — | 300 trees, depth=6, lr=0.1, subsample=0.8 |
| RandomForest | — | 300 trees, depth=12, min_samples_leaf=3 |
| GradientBoosting | — | 300 trees, depth=5, lr=0.1 |
| LinearRegression | — | Stats-only features |

### Classification (binned count)

| Model | Grid Search | Key Parameters |
|---|---|---|
| MajorityClassBaseline | — | Predicts majority bin |
| **XGBClassifier+Grid** | 81 combinations | n_est=[100,300], depth=[4,6,8], lr=[0.05,0.1,0.2], subsample=[0.7,0.85,1.0] |
| XGBClassifier | — | 300 trees, depth=6, lr=0.1 |
| RandomForestClassifier | — | 300 trees, depth=12 |
| GradientBoostingClassifier | — | 300 trees, depth=5, lr=0.1 |

## Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--study_output_filepath` | Yes | — | Path to study output directory |
| `--taxid_plan_filepath` | Yes | — | Path to taxid plan TSV |
| `--analysis_output_filepath` | Yes | — | Root path for outputs |
| `--output_dir` | No | `composition_prediction` | Subdirectory name |
| `--tax_level` | No | `family` | Taxonomic level for composition features |
| `--data_set_divide` | No | `16` | Number of recall divisions |
| `--max_training` | No | — | Max training datasets to use |
| `--recall_cache` | No | — | Path to recall_results_cache.parquet |
| `--verbose` | No | false | Verbose logging |

## Example

```bash
python deployment/model_evaluation/analysis_scripts/input_composition_prediction.py \
    --study_output_filepath /path/to/study \
    --taxid_plan_filepath /path/to/taxid_plan.tsv \
    --analysis_output_filepath /path/to/output \
    --output_dir composition_prediction \
    --tax_level family
```

## Outputs

Saved to `{analysis_output_filepath}/{output_dir}/`:

### Tables

#### `metrics_regression.csv`

One row per regression model, sorted by R² descending.

| Column | Description |
|---|---|
| `Model` | Model name (MeanBaseline, XGBoost+Grid, XGBoost, RandomForest, GradientBoosting, LinearRegression) |
| `R2` | Coefficient of determination |
| `MAE` | Mean absolute error (predicted − true count) |
| `RMSE` | Root mean squared error |
| `MAPE` | Mean absolute percentage error |
| `Acc@1` | Percentage of predictions within ±1 of the true count |
| `Acc@2` | Percentage of predictions within ±2 of the true count |
| `Acc@5` | Percentage of predictions within ±5 of the true count |
| `Train_Time_s` | Wall-clock training time in seconds |

#### `metrics_classification.csv`

One row per classification model, sorted by Accuracy descending.

| Column | Description |
|---|---|
| `Model` | Model name (MajorityClassBaseline, XGBClassifier+Grid, XGBClassifier, RandomForestClassifier, GradientBoostingClassifier) |
| `Accuracy` | Exact bin-match rate |
| `Balanced_Accuracy` | Average recall per bin (accounts for class imbalance) |
| `Macro_F1` | Unweighted mean F1 across bins |
| `Macro_Precision` | Unweighted mean precision across bins |
| `Macro_Recall` | Unweighted mean recall across bins |
| `Ordinal_MAE` | Mean absolute error in bin index (0–4) |
| `Ordinal_Acc@0` | Percentage of exact bin matches (same as Accuracy) |
| `Ordinal_Acc@1` | Percentage of predictions within ±1 bin |
| `Train_Time_s` | Wall-clock training time in seconds |

#### `detailed_predictions.tsv`

Per-dataset predictions for both regression and classification tasks.

| Column | Description |
|---|---|
| `dataset_name` | Dataset identifier |
| `true_count` | Actual number of unique input taxids |
| `predicted_count` | Regression prediction of unique taxid count |
| `true_bin` | True ordinal bin (0 = very_low … 4 = very_high) |
| `predicted_bin` | Predicted ordinal bin (0–4) |

#### `best_params.txt`

Best GridSearchCV hyperparameters for XGBoost+Grid (regression, 243 combinations) and XGBClassifier+Grid (classification, 81 combinations). Includes `n_estimators`, `max_depth`, `learning_rate`, `subsample`, `colsample_bytree`, and the corresponding CV score.

### Plots

#### `predictions_vs_actual.png`

Facet grid with one panel per regression model. X-axis: true count of unique input taxids. Y-axis: predicted count. Diagonal dashed line = perfect prediction. Points are individual datasets. Shows systematic bias (points consistently above/below the diagonal).

#### `residuals.png`

One histogram per regression model. X-axis: residual (true − predicted count). Y-axis: frequency (number of datasets with that residual). Vertical dashed line at zero. Skew or spread indicates bias or variance issues.

#### `confusion_matrix.png`

Grid of normalised confusion matrices — one panel per classification model. Rows: true bin label (very_low → very_high). Columns: predicted bin label. Cell values are row-normalised fractions. Colormap: Blues. Ideal model has a dark diagonal.

#### `classification_metrics.png`

Grouped bar chart per classification model. Each model has three bars: Accuracy, MacroF1, OrdinalMAE (rescaled to [0,1] for display). X-axis: model name. Y-axis: score.

#### `feature_importance_regression.png`

Top-15 horizontal bar charts for each tree-based regression model (XGBoost+Grid, XGBoost, RandomForest, GradientBoosting). X-axis: importance score. Y-axis: feature name. LinearRegression (coefficient-based) and MeanBaseline are excluded.

#### `feature_importance_classification.png`

Same layout as regression, for tree-based classification models (XGBClassifier+Grid, XGBClassifier, RandomForestClassifier, GradientBoostingClassifier).

#### `training_time.png`

Bar chart comparing training times across all models. X-axis: model name (grouped by regression vs classification). Y-axis: training time in seconds. Bars coloured by model family for visual grouping.

## Integration with Last-TP Division Prediction

The predicted N-taxid count is used as an additional feature in
[`last_tp_division_prediction_second.py`](last_tp_division_prediction.md).

### Data Flow

1. `last_tp_division_prediction_second.py` loads the recall cache and the
   study output directory
2. For each dataset in the cache, loads `{dataset}/input/{dataset}.tsv` and
   computes `taxid.nunique()`
3. Merges this true N-taxid count into the cache DataFrame (column
   `n_taxids_true`)
4. After the standard train/test split and feature scaling, trains an
   **XGBoost regressor** (300 trees, max_depth=6, log1p-transformed) to
   predict N-taxids from the same stat + taxonomy features
5. Appends the predicted N-taxid count as an additional feature column to
   both train and test feature matrices
6. All 13 model variants (GP kernels, Hurdle, RF, XGBoost, NGBoost,
   DirectXGB, DirectRF, GP-DirectXGB, GPCLFThreshold) are trained and
   evaluated on the augmented feature set

### Usage with Integration

```bash
python last_tp_division_prediction_second.py \
    --study_output_filepath /path/to/study \
    --input_cache /path/to/recall_results_cache.parquet \
    --output_dir /path/to/outputs
```

Without `--study_output_filepath`, the behaviour is unchanged (no N-taxid
augmentation).

### Diagnostics

At runtime, the N-taxid predictor prints:

- Test R² and MAE of the N-taxid predictor on the held-out set
- Range of predicted N-taxids (train and test)
- Range of true N-taxids (train and test)

---

## Evaluation Metrics

### Regression

| Metric | Description |
|---|---|
| R² | Coefficient of determination |
| MAE | Mean absolute error |
| RMSE | Root mean squared error |
| MAPE | Mean absolute percentage error |
| Acc@1 | % of predictions within ±1 of true count |
| Acc@2 | % of predictions within ±2 of true count |
| Acc@5 | % of predictions within ±5 of true count |

### Classification (binned)

| Metric | Description |
|---|---|
| Accuracy | Exact bin-match rate |
| Balanced Accuracy | Average recall per bin |
| Macro F1 | Unweighted F1 across bins |
| Macro Precision | Unweighted precision across bins |
| Macro Recall | Unweighted recall across bins |
| Ordinal MAE | Mean absolute error in bin index |
| Ordinal Acc@0 | % exact bin match |
| Ordinal Acc@1 | % within ±1 bin |
