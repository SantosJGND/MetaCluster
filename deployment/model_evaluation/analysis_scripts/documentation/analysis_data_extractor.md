# Analysis Data Extractor

Extracts structured analysis data from simulation outputs for manuscript result generation. Processes study output directories, computes per-dataset metrics (recall, precision, cross-hit/spurious counts, relindex), breaks down recall by classifier, and optionally runs explanatory analyses using mixed-effects models, quantile regression, and GEE-based recall probability models.

**Note:** This script lives at `deployment/model_evaluation/analysis_data_extractor.py` (the parent of `analysis_scripts/`), not in `analysis_scripts/` itself.

## Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--study-output` | Yes | — | Path to the study output directory (contains per-dataset folders) |
| `--ncbi-db` | Yes | — | Path to the NCBI taxonomy database (`taxa.db`) |
| `--output-dir` | Yes | — | Directory where all analysis outputs will be saved |
| `--cross-hit-threshold` | No | `0.3` | Cross-hit threshold for m-stats matrix |
| `--min-taxonomic-score` | No | `0.7` | Minimum taxonomic score for m-stats matrix |
| `--explanatory` | No | false | If set, runs explanatory analysis (mixed-effects, quantile regression, GEE) |

## Example

```bash
PYTHONPATH=/path/to/project/root python deployment/model_evaluation/analysis_data_extractor.py \
    --study-output /path/to/study_output \
    --ncbi-db /path/to/taxa.db \
    --output-dir /path/to/output
```

## Methods

### Data Extraction

Per dataset, the script loads:
- Input summary from `{dataset}/{dataset}.tsv` (taxid, reads, mutation_rate, ranks)
- Clustering data via `OverlapManager`
- Matched assemblies from `output/matched_assemblies.tsv`
- Classification results from `classification/{dataset}_merged_classification.tsv`

The m-stats matrix is computed via `get_m_stats_matrix`, producing per-leaf statistics including `best_match_is_best`, `is_crosshit`, `is_trash`, `best_match_taxid`, `numreads`, `total_uniq_reads`, `coverage`, and taxonomic ranks.

### Metrics Computed

- **Recall**: Fraction of input taxids recovered in output (via `compute_recall`)
- **Recall (cov-filtered)**: Recall after coverage filtering
- **Relindex**: Relative index of the last "best match is best" entry, normalised by total m-stats length (`last_best_match_relindex`)
- **Cross-hit / Spurious counts**: Per-dataset counts of cross-hits and spurious assignments
- **Precision**: Strict (exact) and approximate precision
- **Per-classifier recall**: Recall broken down by individual classifier

### Explanatory Analysis (when `--explanatory` is set)

#### Mixed-Effects Models

Fitted via `statsmodels.MixedLM` for TP read count targets (`numreads`, `total_uniq_reads`) and `taxid_relindex`. Formula:
```
log_target ~ log_reads + mutation_rate + I(mutation_rate**2) + coverage + taxid_relindex
```
Random intercepts for taxonomic order and optional variance components for dataset.

#### Quantile Regression

Upper 95th percentile quantile models for the same targets, via `statsmodels.QuantReg`.

#### GEE Recall Models

Three formula variants (A, B, C) predicting binary recall (`recalled`) using `log_reads`, `mutation_rate`, `order` as predictors, with `Binomial` family and `Independence` covariance structure, grouped by `data_set`.

> **Target definition.** `recalled` is **classification-aware**: 1 if the input taxid is best-matched in the m-stats matrix *or* the classifier called it with `uniq_reads > 0`. The assembly-only indicator is kept alongside as `recalled_assembly` (and `recalled_classification` for the classifier-only view) so an assembly-based GEE can be refit from the same `recall_data.tsv` by swapping the response in `RECALL_FORMULAS`. Coefficient tables and probability surfaces generated before this change used the assembly-only target (which was additionally coverage-filtered) and are not comparable. See `output_description.md` → "Recall definitions".

## Outputs

### Core Tables

#### `per_dataset_metrics.tsv`

One row per dataset.

| Column | Description |
|---|---|
| `data_set` | Dataset name |
| `n_input_taxids` | Number of unique taxids in the ground-truth input |
| `n_output_leaves` | Number of leaves in the raw m-stats matrix (output size) |
| `n_tp` | True-positive taxids (input taxids found in output) |
| `n_cross_hit` | Cross-hit count |
| `n_spurious` | Spurious assignment count |
| `overall_recall` | Classification-aware recall: input taxids that are best-matched in the m-stats matrix **or** called by the classifier with `uniq_reads > 0`, divided by `n_input_taxids` |
| `recall_cov_filtered` | Assembly-based recall restricted to best-matched leaves with coverage > 0 |
| `recall_assembly` | Assembly-based recall (best-matched leaves, zero-coverage included) |
| `recall_classification` | Classifier-only recall: input taxids present in `<ds>_merged_classification.tsv` with `uniq_reads > 0` |
| `last_best_match_relindex` | Relative index [0, 1] of the last correct `best_match_is_best` entry, normalised by total m-stats length |
| `precision_strict` | Strict precision (output matches must be exact) |
| `precision_approx` | Approximate precision (allowing taxonomic relaxation) |

#### `aggregate_statistics.tsv`

Summary statistics computed across all datasets for each numeric column in `per_dataset_metrics.tsv`.

| Column | Description |
|---|---|
| `metric` | Name of the metric column |
| `count` | Number of non-null values |
| `mean` | Arithmetic mean |
| `std` | Standard deviation |
| `min` | Minimum value |
| `25%` | First quartile |
| `50%` | Median |
| `75%` | Third quartile |
| `max` | Maximum value |

#### `recall_per_classifier.tsv`

One row per (dataset, classifier) pair.

| Column | Description |
|---|---|
| `data_set` | Dataset name |
| `classifier` | Classifier name |
| `tp_taxids` | Number of input taxids found by this classifier |
| `total_input` | Total number of input taxids in ground truth |
| `recall` | Recall fraction for this classifier on this dataset |

#### `recall_per_classifier_aggregate.tsv`

One row per classifier, aggregated across datasets.

| Column | Description |
|---|---|
| `classifier` | Classifier name |
| `n_datasets` | Number of datasets evaluated |
| `mean_recall` | Mean recall across datasets |
| `median_recall` | Median recall |
| `std_recall` | Standard deviation of recall |
| `min_recall` | Minimum recall |
| `max_recall` | Maximum recall |

#### `classifier_hit_counts.tsv`

One row per (dataset, classifier) pair.

| Column | Description |
|---|---|
| `data_set` | Dataset name |
| `classifier` | Classifier name |
| `total_hits` | Total assignments by this classifier |
| `tp_hits` | True-positive assignments |
| `cross_hits` | Cross-hit assignments |
| `trash` | Trash (discarded) assignments |

### Core Plots

#### `relindex_distribution.png`

Histogram of `last_best_match_relindex` across all datasets. X-axis: relindex [0, 1]. Y-axis: number of datasets. Shows how far into each sorted output list the last correct assignment appears — values near 0 indicate the last correct taxid is found early; values near 1 indicate nearly the entire list must be inspected.

#### `dataset_size_distribution.png`

Histogram of raw output leaf counts (`n_output_leaves`). X-axis: number of leaves. Y-axis: number of datasets. Illustrates the size distribution of clustering outputs across datasets.

#### `tp_cross_spurious_distribution.png`

Boxplot with three boxes. Y-axis: count. The three boxes show the distribution of TP counts, cross-hit counts, and spurious counts across datasets. Outliers plotted as individual points.

#### `recall_per_classifier.png`

Boxplot comparing recall across classifiers. X-axis: each classifier name, plus an "overall" box (using all classifiers jointly). Y-axis: recall fraction [0, 1]. Green box = overall recall. Each box shows the distribution of per-dataset recall values.

#### `precision_global_distribution.png`

Overlaid histograms of strict precision and approximate precision across all datasets. X-axis: precision [0, 1]. Y-axis: number of datasets. Two fills distinguish strict vs approximate. Shows how much precision changes when taxonomic relaxation is allowed.

#### `precision_per_classifier.png`

Boxplot comparing precision across classifiers. X-axis: each classifier name, plus overall. Y-axis: precision [0, 1]. Two panels or grouped boxes per classifier for strict and approximate precision.

### Explanatory Outputs (when `--explanatory` is set)

#### `explanatory/tp_hit_data.tsv`

TP-level predictor data for mixed-effects and quantile modelling. One row per (dataset, taxid) TP hit.

| Column | Description |
|---|---|
| `data_set` | Dataset name |
| `taxid` | Taxonomic ID |
| `numreads` | Read count for this taxid |
| `total_uniq_reads` | Total unique reads assigned to this leaf |
| `coverage` | Coverage value |
| `log_reads` | Log-transformed reads |
| `mutation_rate` | Mutation rate |
| `taxid_relindex` | Relative index of this taxid in the sorted list |
| `order` | Taxonomic order (random intercept group) |
| `recalled` | Binary recall indicator (1 = taxid best-matched in output **or** detected by the classifier with `uniq_reads > 0`) |
| `recalled_assembly` | 1 = taxid best-matched in the m-stats matrix (assembly-based) |
| `recalled_classification` | 1 = taxid called by the classifier with `uniq_reads > 0` |

#### `explanatory/model_{target}_summary.tsv`

Mixed-effects model coefficients for the given target (`numreads`, `total_uniq_reads`, or `taxid_relindex`).

| Column | Description |
|---|---|
| `coef` | Fixed-effect coefficient estimate |
| `std_err` | Standard error of the coefficient |
| `z` | Z-statistic |
| `P>|z|` | Two-sided p-value |
| `[0.025` | Lower bound of 95% confidence interval |
| `0.975]` | Upper bound of 95% confidence interval |

#### `explanatory/model_{target}_random_effects.tsv`

Random intercepts per taxonomic order.

| Column | Description |
|---|---|
| `group` | Taxonomic order name |
| `intercept` | Random intercept deviation from the fixed-effect mean |

#### `explanatory/model_{target}_variance_components.tsv`

Variance component estimates for each random effect term (e.g., order intercept, dataset).

#### `explanatory/diagnostic_{target}.png`

2 × 2 diagnostic grid for the mixed-effects model:
- **Top-left**: Residuals vs fitted values (check for homoscedasticity and patterns)
- **Top-right**: Q-Q plot of residuals (check for normality)
- **Bottom-left**: Observed vs predicted values (calibration)
- **Bottom-right**: Histogram of random intercepts (check for normality of random effects)

#### `explanatory/quantile_{target}_summary.tsv`

Quantile regression (upper 95th percentile) coefficients.

| Column | Description |
|---|---|
| `coef` | Coefficient estimate |
| `std_err` | Standard error |
| `t` | T-statistic |
| `P>|t|` | Two-sided p-value |
| `[0.025` | Lower bound of 95% confidence interval |
| `0.975]` | Upper bound of 95% confidence interval |

#### `explanatory/diagnostic_quantile_{target}.png`

Diagnostic plot for the quantile regression model. Includes residuals vs fitted and Q-Q plot.

#### `explanatory/scatter_predictors_vs_target.png`

PairGrid scatter matrix of predictors (log_reads, mutation_rate, coverage, taxid_relindex) with the target on the diagonal histograms. Points are coloured by taxonomic order. Shows pairwise relationships and potential interactions.

#### `explanatory/recall_model/recall_data.tsv`

Binary recall data for GEE modelling. One row per (dataset, taxid).

| Column | Description |
|---|---|
| `data_set` | Dataset name |
| `taxid` | Taxonomic ID |
| `reads_simulated` | Reads simulated for this taxid |
| `mutation_rate` | Mutation rate |
| `order` | Taxonomic order |
| `recalled` | Binary outcome (1 = taxid best-matched in the m-stats matrix **or** detected by the classifier with `uniq_reads > 0`) — GEE target |
| `recalled_assembly` | 1 = taxid best-matched in the m-stats matrix only (assembly-based, secondary) |
| `recalled_classification` | 1 = taxid called by the classifier with `uniq_reads > 0` |

> `log_reads = log1p(reads_simulated)` is derived at fit time (`fit_recall_gee`); it is **not** stored in this file. This file is sufficient to replay the recall surfaces (see below).

#### `explanatory/recall_model/recall_variant_{A,B,C}_summary.tsv`

GEE coefficient summaries for each formula variant. Terms are keyed by their
patsy formula name (first column, e.g. `Intercept`, `log_reads`, `mutation_rate`,
`I(mutation_rate ** 2)`, `C(order)[T.X]`).

| Column | Description |
|---|---|
| (first, term names) | Patsy term name |
| `coef` | Coefficient estimate |
| `se` | Robust standard error (sandwich estimator) |
| `z` | Z-statistic |
| `p_value` | Two-sided p-value |
| `ci_lower` | Lower bound of the 95% confidence interval |
| `ci_upper` | Upper bound of the 95% confidence interval |

#### `explanatory/recall_model/recall_calibration.png`

Calibration curve (X-axis: predicted probability, Y-axis: observed frequency) showing how well predicted recall probabilities match empirical rates, with ROC curve inset. Perfect calibration follows the diagonal.

#### `explanatory/recall_model/recall_probability_surface.png`

Contour plot panel per fitted GEE formula variant (A, B, C). X-axis: `log_reads`,
Y-axis: `mutation_rate`. Colour fill: predicted recall probability. Shows the
interaction surface between read depth and mutation rate on recall probability.
Lighter regions indicate higher recall probability. Each panel is evaluated at the
most frequent taxonomic order by default (`order_mode="reference"`).

**Replaying from a saved directory (no refit):** the surface can be reconstructed
purely from `recall_data.tsv` + the three `recall_variant_{A,B,C}_summary.tsv` files
(the design matrix is rebuilt via `patsy` from the training data's `design_info`,
so predictions match the fitted models):

```bash
python deployment/model_evaluation/analysis_data_extractor.py \
    --replay-recall <output>/explanatory/recall_model \
    --output-dir <target> \
    [--replay-order-mode reference|average]
```

`--replay-order-mode average` instead evaluates the surface at **every** taxonomic
order in `recall_data.tsv` and plots their mean (average across classes).
