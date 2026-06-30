# Analysis Scripts

Experimental comparison scripts for evaluating model variants on viral
metagenomic data. Each script tests alternative approaches against the
production pipeline and reports side-by-side metrics.

## Setup

```bash
source .venv/bin/activate
export PYTHONPATH=$(pwd)
```

All scripts require prior evaluation pipeline output (study output directory
or cached training data).

---

## 1. `compare_sort_strategies.py`

Sweeps `sort_strategy` × `recall_model_interface` over cached or freshly
retrieved training data, then evaluates each combination on held-out test
datasets. Reports `target_percentile`, `keep_index`, `actual_recall`,
`max_possible_recall`, and `standardized_recall`.

**Sort strategies:** `reads`, `taxid_roundrobin`, `rarity_boost`, `tax_level_stratified`

**Recall models:** `xgb`, `morf`, `direct_xgb`, `gp_clf`, `fixed_12_taxids`
(plus any `fixed_N` count)

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--study_output_filepath` | Yes | - | Path to study output directory |
| `--taxid_plan_filepath` | Yes | - | Path to taxid plan TSV |
| `--analysis_output_filepath` | Yes | - | Root path for analysis outputs (plots, models, etc.) |
| `--output_dir` | No | `sort_comparison` | Subdirectory name within analysis_output_filepath for comparison results |
| `--sort_strategies` | No | `reads taxid_roundrobin rarity_boost tax_level_stratified` | Sort strategies to compare |
| `--recall_models` | No | `xgb morf direct_xgb gp_clf fixed_12_taxids` | Recall model types to compare |
| `--target_recall` | No | `1.0` | Target recall threshold |
| `--data_set_divide` | No | `16` | Number of recall divisions |
| `--tax_level` | No | `genus` | Taxonomic level |
| `--max_training` | No | - | Max training datasets to use |
| `--verbose` | No | false | Verbose logging |

### Example

```bash
python deployment/model_evaluation/analysis_scripts/compare_sort_strategies.py \
    --study_output_filepath /path/to/study \
    --taxid_plan_filepath /path/to/taxid_plan.tsv \
    --analysis_output_filepath /path/to/analysis/output \
    --sort_strategies reads taxid_roundrobin \
    --recall_models xgb direct_xgb gp_clf \
    --target_recall 0.95 \
    --output_dir sort_comparison_rerun
```

### Outputs

Saved to `{analysis_output_filepath}/{output_dir}/`:

| File | Description |
|------|-------------|
| `sort_strategy_comparison_results.tsv` | Per-dataset results for every (strategy, model) combination |
| `sort_strategy_comparison_summary.tsv` | Aggregated statistics by (strategy, model) |
| `comparison_target_percentile.png` | Boxplot of target percentile by strategy, faceted by model |
| `comparison_actual_recall.png` | Boxplot of actual recall at predicted cutoff |
| `comparison_standardized_recall.png` | Boxplot of standardized recall (actual / max possible) |
| `comparison_heatmap.png` | Mean target percentile heatmap |
| `comparison_percentile_vs_recall.png` | Density grid of target percentile vs standardized recall |
| `comparison_true_vs_predicted_cutoff.png` | Density grid of true vs predicted cutoff |

---

## 2. `composition_model_comparison.py`

Experimental side-by-side comparison of classifiers for the `stop_traversal`
(clade composition) prediction task. Loads cached training data, replicates
production preprocessing, trains and evaluates multiple classifiers.

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--output_dir` | No | `composition_comparison_outputs` | Output directory for results |
| `--input_cache` | No | `deployment/ml_api/training_results_cache.parquet` | Path to training data cache |

### Example

```bash
python deployment/model_evaluation/analysis_scripts/composition_model_comparison.py \
    --output_dir /path/to/composition_comparison_outputs \
    --input_cache deployment/ml_api/training_results_cache.parquet
```

### Methods

**Dataset.** Training data comes from `deployment/ml_api/training_results_cache.parquet`,
produced by the evaluation pipeline's clade-traversal step
(`data_set_traversal_with_precision()`). The dataset comprises 4,920 rows and
23 columns: 11 family-level viral taxonomy proportions (Baculoviridae,
Coronaviridae, Flaviviridae, Herelleviridae, Orthoherpesviridae,
Papillomaviridae, Peduoviridae, Picornaviridae, Rhabdoviridae, Schitoviridae,
Straboviridae), 4 statistical features (`n_leaves`, `tax_diversity`, `Min_Dist`,
`Min_Shared`), 7 metadata columns (`data_set`, `node`, `n_true_leaves`,
`precision`, `new_precision`, `precision_increased`, `unclassified`), and the
binary target `stop_traversal`. The target distribution is 65.3% positive
(stop) and 34.7% negative (continue).

**Preprocessing.** Metadata columns are dropped. The remaining 15 features
are split into two groups: **stat features** (4 numeric columns) and
**taxonomy features** (11 composition columns). A stratified 80/20
train-test split is applied (random_state=42). Stat features are standardised
via `StandardScaler`, taxonomy features (already [0,1] proportions) are left
unaltered.

**Models compared:**
- **FixedThreshold(min_shared >= 0.6)** — rule-based baseline, no training
- **Random Forest** — 300 trees, max_depth=12, balanced class_weight
- **XGBoost** — 300 trees, max_depth=6, lr=0.1
- **Gradient Boosting** — 300 trees, max_depth=5, lr=0.1
- **Logistic Regression** — L2-regularised, stats-only features
- **XGBoost + Optuna** (if optuna available) — 50-trial hyperparameter search
- **LightGBM** (if lightgbm available) — 300 trees, max_depth=6

**Evaluation metrics:** Accuracy, precision, recall, F1, ROC-AUC, PR-AUC
(positive class = stop_traversal). ROC and PR curves generated for each model
offering predicted probabilities.

### Outputs

| File | Description |
|------|-------------|
| `comparison_results.csv` | Metrics table sorted by F1 |
| `classification_report.txt` | Full sklearn report per model |
| `roc_curves.png` | ROC curves with AUC labels |
| `pr_curves.png` | Precision-Recall curves |
| `f1_comparison.png` | F1-score bar chart |
| `precision_recall_bars.png` | Grouped precision/recall/F1 bars |
| `confusion_matrices.png` | Normalised confusion matrix grid |
| `feature_importance.png` | Top-15 feature importances (tree models) |
| `training_time.png` | Training-time comparison |

---

## 3. `last_tp_division_prediction_second.py`

Predicts the fraction (division) along a sorted ranked list where recall
first crosses threshold τ. Compares 4 decision approaches (CLF, Likelihood,
Direct, Direct+ISO) across 13 model variants including 6 Baseline GP kernels,
Hurdle (logistic + GP), multi-output regressors (RF, XGBoost, NGBoost), and
direct fraction regressors (DirectXGB, DirectRF, GP-DirectXGB).

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--output_dir` | No | `last_tp_division_outputs` | Output directory for results |
| `--input_cache` | No | `/path/to/recall_results_cache.parquet` | Path to recall results cache parquet |

### Example

```bash
python deployment/model_evaluation/analysis_scripts/last_tp_division_prediction_second.py \
    --output_dir /path/to/last_tp_division_outputs \
    --input_cache /path/to/recall_results_cache.parquet
```

### Methods

**Dataset.** Training data is sourced from the GPCLF model cache
(`recall_results_cache.parquet`), produced by the evaluation pipeline's
recall-modelling step. The raw dataset comprises 933 samples with 20 recall
divisions (`index_recall_1` … `index_recall_20`), 5 statistical features
(`counts_kurtosis`, `counts_skewness`, `tax_diversity_shannon`,
`max_uniq_reads`, `total_uniq_reads`), 16 taxonomy (family-level) proportion
columns, `last_best_match_relindex`, and metadata.

Samples where `last_best_match_relindex ∈ {0, 1}` (degenerate edge cases) are
removed, along with samples whose recall never exceeds 0 across all divisions.
This yields 787 samples (630 train, 157 test).

**Preprocessing.** A stratified 80/20 train-test split is applied
(random_state=42). The 5 statistical features are standardised via
`StandardScaler`; taxonomy proportions are left as raw [0, 1] values. The
division grid uses 20 bins at fractions {0.05, 0.10, …, 1.00}. The last 4
divisions (fractions > 0.80) are effectively constant at recall = 1.0 and are
treated as inactive — models predict on the first 16 active divisions and the
remaining 4 are padded.

**Model families:**
- **Baseline GP (6 kernels)** — Independent Gaussian process regressor per
  active division. Kernels: RBF, Matern-0.5, Matern-2.5, RationalQuadratic,
  RBF+RQ, Matern×Matern.
- **Hurdle model** — Logistic regression (spike classifier) + Matern GP
  (conditional continuous model), combined as a mixture.
- **Deterministic multi-output regressors** — Random Forest, XGBoost, NGBoost
  (all with per-division NGBoost via CRPS grid search).
- **Direct fraction regressors** — DirectXGB, DirectRF, GP-DirectXGB (stacked
  meta-model with GP features + XGBoost).
- **GPCLFThreshold** — Production baseline with RationalQuadratic GP and CLF
  threshold search.

**Decision approaches:** Each probabilistic model is evaluated with two rules:
- **CLF** — first division where P(recall ≥ τ) ≥ X (sweep X ∈ [0.025, 0.975])
- **Likelihood** — first division where pdf(τ; μ, σ) ≥ τ_lik
- **Direct** — deterministic models compute τ-crossing fraction from mean
- **Direct+ISO** — Direct + isotonic regression calibration

**Evaluation protocol.** An asymmetric squared-error loss penalises over- and
under-estimation of the τ-crossing fraction with configurable ratio (default
over:under = 1:3). Models are evaluated at τ ∈ {0.85, 0.90, 0.95, 0.98, 1.00}.
For each (model, approach, ratio, τ) combination, mean and median actual recall,
percentage below τ, and asymmetric loss are reported.

### Outputs

**Summary Tables:**
- `summary_results_tau{τ}.csv` — Best loss per (model, approach, asymmetry ratio)
- `full_grid_results_tau{τ}.parquet` — Full grid search results
- `actual_recall_at_index_summary.csv` — All (model, approach) pairs at τ = 1.00
- `tau_sweep_summary.csv` — Best model per τ

**Loss Analysis Plots:**
- `loss_curves_{model}_tau{τ}.png` — Loss vs. decision parameter
- `best_loss_comparison_tau{τ}.png` — Heatmap of best loss
- `sensitivity_ratio.png` — Loss sensitivity to asymmetry ratio

**Calibration & Readout:**
- `calibration_{model}_{approach}_tau{τ}.png` — Predicted vs. true
- `actual_recall_at_index_{model}_{approach}_tau{τ}.png` — Recall histogram
- `recall_landscape.png` — Recall distribution across divisions
- `list_savings_GPCLFThreshold.png` — Saved list fraction distribution

**Cross-τ Comparison:**
- `roc_curves.png` — TPR vs. FPR sweep
- `tau_sweep_actual_recall.png` — Mean actual recall vs. target τ
- `tau_sweep_tradeoff.png` — Loss vs. recall tradeoff
