# Last-TP Division Prediction

Predicts the fraction (division) along a sorted ranked list where recall first crosses threshold τ. Compares 4 decision approaches (CLF, Likelihood, Direct, Direct+ISO) across 13 model variants including 6 Baseline GP kernels, Hurdle (logistic + GP), multi-output regressors (RF, XGBoost, NGBoost), and direct fraction regressors (DirectXGB, DirectRF, GP-DirectXGB). Features an integrated N-taxid predictor for feature augmentation.

## Task

Given statistical and taxonomic features of a dataset, predict at which percentile of sorted results the recall reaches a target τ. This determines how much of the ranked list needs to be inspected to achieve a given recall level.

## Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--output_dir` | No | `last_tp_division_outputs` | Output directory for results |
| `--input_cache` | No | `recall_results_cache.parquet` | Path to recall results cache parquet |
| `--study_output_filepath` | No | — | Study output directory (enables N-taxid augmentation) |

When `--study_output_filepath` is provided, the script additionally trains an XGBoost regressor to predict the number of unique input taxids (`n_taxids_true`) and appends this prediction as an extra feature for all 13 models.

## Example

```bash
# Without N-taxid augmentation
python deployment/model_evaluation/analysis_scripts/last_tp_division_prediction_second.py \
    --output_dir /path/to/last_tp_division_outputs \
    --input_cache /path/to/recall_results_cache.parquet

# With N-taxid augmentation
python deployment/model_evaluation/analysis_scripts/last_tp_division_prediction_second.py \
    --study_output_filepath /path/to/study \
    --input_cache /path/to/recall_results_cache.parquet \
    --output_dir /path/to/outputs
```

## Methods

### Dataset

Training data is sourced from the GPCLF model cache (`recall_results_cache.parquet`), produced by the evaluation pipeline's recall-modelling step. The raw dataset comprises 933 samples with 20 recall divisions (`index_recall_1`…`index_recall_20`), 5 statistical features (`counts_kurtosis`, `counts_skewness`, `tax_diversity_shannon`, `max_uniq_reads`, `total_uniq_reads`), 16 taxonomy (family-level) proportion columns, `last_best_match_relindex`, and metadata.

Samples where `last_best_match_relindex ∈ {0, 1}` (degenerate edge cases) are removed, along with samples whose recall never exceeds 0 across all divisions. This yields 787 samples (630 train, 157 test).

### Preprocessing

A stratified 80/20 train-test split is applied (`random_state=42`). The 5 statistical features are standardised via `StandardScaler`; taxonomy proportions are left as raw [0, 1] values. The division grid uses 20 bins at fractions {0.05, 0.10, …, 1.00}. The last 4 divisions (fractions > 0.80) are effectively constant at recall = 1.0 and are treated as inactive — models predict on the first 16 active divisions and the remaining 4 are padded.

### Model Families

#### Baseline GP (6 kernels)

Independent Gaussian process regressor per active division. Each division gets its own GP with one of:

| Kernel | Description |
|---|---|
| RBF | Radial basis function |
| Matern-0.5 | Matern(ν=0.5) — equivalent to exponential |
| Matern-2.5 | Matern(ν=2.5) |
| RationalQuadratic | Rational quadratic |
| RBF+RQ | Sum of RBF and Rational Quadratic |
| Matern×Matern | Product of two Matern kernels (different length scales) |

#### Hurdle Model

Two-stage: logistic regression (spike classifier for P(recall=0)) + Matern GP (conditional continuous model for P(recall>0)), combined as a mixture.

#### Deterministic Multi-Output Regressors

| Model | Description |
|---|---|
| RandomForest | Multi-output RF, 300 trees |
| XGBoost | Multi-output XGBoost, 300 trees |
| NGBoost | Natural Gradient Boosting, per-division via CRPS grid search |

#### Direct Fraction Regressors

| Model | Description |
|---|---|
| DirectXGB | XGBoost predicting τ-crossing fraction directly |
| DirectRF | Random Forest predicting τ-crossing fraction directly |
| GP-DirectXGB | Stacked meta-model: GP features → XGBoost |

#### Production Baseline

GPCLFThreshold — RationalQuadratic GP with CLF threshold search.

### Decision Approaches

Each probabilistic model is evaluated with two rules:

| Approach | Description |
|---|---|
| CLF | First division where P(recall ≥ τ) ≥ X (sweep X ∈ [0.025, 0.975]) |
| Likelihood | First division where pdf(τ; μ, σ) ≥ τ_lik |
| Direct | Deterministic models compute τ-crossing fraction from mean prediction |
| Direct+ISO | Direct + isotonic regression calibration |

### N-Taxid Feature Augmentation

When `--study_output_filepath` is provided, the script:

1. Loads `{dataset}/input/{dataset}.tsv` for each dataset in the cache
2. Computes `taxid.nunique()` as `n_taxids_true`
3. Merges this into the cache DataFrame
4. Trains an XGBoost regressor (300 trees, max_depth=6, log1p-transformed) to predict N-taxids from stat + taxonomy features
5. Appends the predicted N-taxid count as an additional feature column to both train and test matrices

At runtime, the N-taxid predictor prints test R², MAE, and value ranges for diagnostics.

### Evaluation Protocol

An asymmetric squared-error loss penalises over- and under-estimation of the τ-crossing fraction with configurable ratio (default over:under = 1:3). Models are evaluated at τ ∈ {0.85, 0.90, 0.95, 0.98, 1.00}. For each (model, approach, ratio, τ) combination, mean and median actual recall, percentage below τ, and asymmetric loss are reported.

## Outputs

### Summary Tables

#### `summary_results_tau{τ}.csv` (generated for each τ ∈ {0.85, 0.90, 0.95, 0.98, 1.00})

Best-performing combination per (model, approach, asymmetry_ratio).

| Column | Description |
|---|---|
| `model` | Model variant name |
| `approach` | Decision approach (CLF, Likelihood, Direct, Direct+ISO) |
| `asymmetry_ratio` | Over:under penalty ratio tested |
| `opt_decision_param` | Optimal decision parameter value (CLF threshold X or likelihood threshold) |
| `best_loss` | Achieved asymmetric loss at optimal decision parameter |
| `mean_actual_recall` | Mean actual recall at predicted cutoff across test datasets |
| `median_actual_recall` | Median actual recall at predicted cutoff |
| `pct_below_tau` | Percentage of test samples where actual recall < target τ |

#### `full_grid_results_tau{τ}.parquet`

Full grid search results: every (model, approach, decision_parameter_value, asymmetry_ratio) combination with the corresponding loss, mean actual recall, median actual recall, and below_tau percentage. Used for detailed analysis and custom filtering.

#### `actual_recall_at_index_summary.csv`

All (model, approach) combinations evaluated at τ = 1.00. Columns include model, approach, mean_actual_recall, median_actual_recall, pct_achieved_target, pct_below_tau. Provides a comprehensive view of recall attainment at the strictest threshold.

#### `tau_sweep_summary.csv`

Best-performing model per τ value.

| Column | Description |
|---|---|
| `tau` | Target recall threshold |
| `model` | Best model variant |
| `approach` | Best decision approach |
| `asymmetry_ratio` | Best asymmetry ratio |
| `best_loss` | Best loss achieved |
| `mean_actual_recall` | Mean actual recall at best configuration |

### Plots

#### Loss Analysis

##### `loss_curves_{model}_tau{τ}.png`

X-axis: decision parameter (CLF threshold X or likelihood threshold across the sweep range). Y-axis: asymmetric loss. One curve per asymmetry ratio (e.g., 1:1, 1:3, 1:5). Shows how loss varies and where the minimum is located.

##### `best_loss_comparison_tau{τ}.png`

Heatmap with model variants on rows and decision approaches on columns. Cell colour = minimum asymmetric loss achieved. Lower (darker) is better. Colour bar with loss scale.

##### `sensitivity_ratio.png`

X-axis: over:under asymmetry ratio. Y-axis: best loss. One line per (model, approach) or aggregated summary. Shows sensitivity of the optimal loss to the choice of asymmetry penalty.

#### Calibration & Readout

##### `calibration_{model}_{approach}_tau{τ}.png`

Scatter plot. X-axis: predicted recall at the model's predicted cutoff, Y-axis: actual recall achieved at that cutoff. Points represent individual test datasets. A diagonal reference line indicates perfect calibration.

##### `actual_recall_at_index_{model}_{approach}_tau{τ}.png`

Histogram of actual recall values at the predicted cutoff across test datasets. X-axis: actual recall [0, 1], Y-axis: count of datasets. Vertical line at target τ.

##### `recall_landscape.png`

X-axis: division index (1–20, corresponding to fractions 0.05–1.00). Y-axis: recall value. Shows mean recall ± standard deviation across datasets at each division, illustrating how recall accumulates through the sorted list.

##### `list_savings_GPCLFThreshold.png`

Histogram of `keep_index / total_leaves` (the fraction of the ranked list retained). X-axis: saved fraction [0, 1], Y-axis: count of datasets. Shows how much of the list the production baseline retains.

#### Cross-τ Comparison

##### `roc_curves.png`

X-axis: False Positive Rate (sweeping the decision threshold), Y-axis: True Positive Rate (budget = τ). One curve per τ value. Diagonal = random. Used to evaluate the binary recall-attainment decision across thresholds.

##### `tau_sweep_actual_recall.png`

X-axis: target τ (0.85, 0.90, 0.95, 0.98, 1.00). Y-axis: mean actual recall. One line per (model, approach) combination. Shows how well each combination tracks the target across difficulty levels.

##### `tau_sweep_tradeoff.png`

X-axis: mean actual recall at predicted cutoff, Y-axis: asymmetric loss. One point per (model, approach, τ) combination. Points from the same (model, approach) are connected across τ values. Ideal models cluster in the upper-left (high recall, low loss).
