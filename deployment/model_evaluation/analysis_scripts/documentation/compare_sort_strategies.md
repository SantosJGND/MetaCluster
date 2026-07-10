# Sort Strategy Comparison

Sweeps `sort_strategy` × `recall_model_interface` over cached or freshly retrieved training data, then evaluates each combination on held-out test datasets. Includes a rule-based baseline (`fixed_12_taxids`) that keeps leaves until N unique `best_match_taxids` are covered after sorting.

## Task

For every (strategy, model) pair, the script:
1. Loads raw m-stats matrices for training and test datasets
2. Applies the sort strategy and extracts features via `RecallFeatureTransformer`
3. Trains the recall model on training data (or skips for rule-based baselines)
4. Predicts a cutoff percentile on each test dataset
5. Evaluates actual recall and standardized recall (actual / max_possible) at that cutoff

## Arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--study_output_filepath` | Yes | — | Path to study output directory |
| `--taxid_plan_filepath` | Yes | — | Path to taxid plan TSV |
| `--analysis_output_filepath` | Yes | — | Root path for analysis outputs |
| `--output_dir` | No | `sort_comparison` | Subdirectory name within analysis_output_filepath |
| `--sort_strategies` | No | `reads taxid_roundrobin rarity_boost tax_level_stratified` | Sort strategies to compare |
| `--recall_models` | No | `xgb morf direct_xgb gp_clf fixed_12_taxids` | Recall model types to compare |
| `--target_recall` | No | `1.0` | Target recall threshold |
| `--data_set_divide` | No | `16` | Number of recall divisions |
| `--tax_level` | No | `genus` | Taxonomic level |
| `--max_training` | No | — | Max training datasets to use |
| `--verbose` | No | false | Verbose logging |

## Example

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

## Sort Strategies

| Strategy | Description |
|---|---|
| `reads` | Sort by read count descending |
| `taxid_roundrobin` | Round-robin across unique taxids |
| `rarity_boost` | Boost rare assignments |
| `tax_level_stratified` | Stratified by taxonomic level |

## Recall Models

| Model | Type | Description |
|---|---|---|
| `xgb` | Trained | Standard XGBoost recall modeller |
| `morf` | Trained | Morf (Random Forest) recall modeller |
| `direct_xgb` | Trained | Direct XGBoost cutoff predictor |
| `gp_clf` | Trained | GP classifier with CLF threshold |
| `fixed_12_taxids` | Rule-based | Keep leaves until 12 unique best_match_taxids covered |

Any `fixed_N` pattern is accepted (e.g. `fixed_20`), creating a rule-based baseline that keeps N unique taxids.

## Methods

### Data Loading

Loads datasets from the study output directory. For each dataset, raw m-stats matrices are loaded via `_get_raw_m_stats_matrix` (using `OverlapManager` and `get_m_stats_matrix`). Input summaries are loaded from `{dataset}/input/{dataset}.tsv`.

### Feature Extraction

`RecallFeatureTransformer` applies the chosen sort strategy to the m-stats matrix, then computes statistical features (kurtosis, skewness, diversity, read counts, etc.) and taxonomy features at the configured `--tax_level`.

### Training & Evaluation Loop

1. Training matrices are collected and transformed
2. For each model type, a `CutoffRecallModeller`, `DirectXGBRecallModeller`, `GPCLFRecallModeller`, `RecallModeller`, or `FixedCountModeller` is instantiated
3. If trainable, an 80/20 train-test split is applied internally
4. The modeller is fitted on training data
5. For each test dataset, the modeller predicts a cutoff percentile, which is converted to a `keep_index`
6. Actual recall, max possible recall, and standardized recall are computed at that keep index

## Outputs

Saved to `{analysis_output_filepath}/{output_dir}/`:

### Tables

#### `sort_strategy_comparison_results.tsv`

One row per (sort_strategy, recall_model, test_dataset) combination.

| Column | Description |
|---|---|
| `sort_strategy` | Sort strategy name (reads, taxid_roundrobin, rarity_boost, tax_level_stratified) |
| `recall_model` | Recall model name (xgb, morf, direct_xgb, gp_clf, fixed_N) |
| `data_set` | Test dataset identifier |
| `total_leaves` | Total number of leaves in the sorted m-stats matrix for this dataset |
| `n_input_taxids` | Number of unique input taxids in the ground-truth input summary |
| `target_percentile` | Predicted cutoff percentile by the model |
| `keep_index` | Integer keep index = `round(target_percentile × total_leaves)`, clamped to ≥ 1 |
| `actual_recall` | Intersection-based `compute_recall` at the selected `keep_index` |
| `max_possible_recall` | Recall using ALL leaves (no cutoff) — theoretical ceiling for standardisation |
| `standardized_recall` | `actual_recall / max_possible_recall`, clipped to [0, 1] |
| `recall_ceiling_gap` | `max_possible_recall − actual_recall` |
| `true_cutoff` | Ground-truth fractional position of the last valid unique `best_match_taxid` |
| `n_taxids_found` | Number of input taxids recovered in the first `keep_index` leaves |
| `recall_gap` | `actual_recall − target_recall` |

#### `sort_strategy_comparison_summary.tsv`

Aggregated per (sort_strategy, recall_model). Each grouped row contains:

| Column | Description |
|---|---|
| `sort_strategy` | Sort strategy name |
| `recall_model` | Recall model name |
| `n_datasets` | Number of test datasets |
| `mean_percentile` | Mean target percentile |
| `std_percentile` | Standard deviation of target percentile |
| `median_percentile` | Median target percentile |
| `mean_keep_index` | Mean keep index |
| `std_keep_index` | Standard deviation of keep index |
| `median_keep_index` | Median keep index |
| `mean_actual_recall` | Mean actual recall |
| `std_actual_recall` | Standard deviation of actual recall |
| `mean_recall_gap` | Mean recall gap (actual − target) |
| `mean_taxids_found` | Mean number of input taxids found |
| `mean_input_taxids` | Mean number of input taxids per dataset |
| `mean_standardized_recall` | Mean standardized recall (available when all models produce max_possible) |
| `std_standardized_recall` | Standard deviation of standardized recall |
| `mean_max_possible_recall` | Mean max possible recall |

Sorted by (sort_strategy, mean_percentile).

### Plots

#### `comparison_target_percentile.png`

Boxplot faceted by recall model (columns). Each panel shows `target_percentile` on the Y-axis with one box per sort strategy. Palette: Set2. X-axis labels rotated 30°. Title: "Target percentile by strategy and model".

#### `comparison_actual_recall.png`

Same faceted boxplot layout. Y-axis: `actual_recall`. Red dashed horizontal line at the median of each panel.

#### `comparison_standardized_recall.png`

Same faceted boxplot layout. Y-axis: `standardized_recall`. Red dashed line at median. Only shown when standardized_recall data is available.

#### `comparison_heatmap.png`

Heatmap with sort_strategy on rows, recall_model on columns. Cell colour = `mean_percentile` (summary table). Colormap: YlOrRd. Annotated with value (.3f). Colour bar labelled "Mean target percentile".

#### `comparison_percentile_vs_recall.png`

2D density grid with rows = recall_model, cols = sort_strategy. Each panel: X-axis = `target_percentile`, Y-axis = `standardized_recall`. KDE fill (alpha 0.35) + black scatter points (alpha 0.15). Grey dashed diagonal. Axis limits [0, 1].

#### `comparison_true_vs_predicted_cutoff.png`

Same density grid layout. X-axis = `true_cutoff` (ground-truth fractional position of the last valid unique best_match_taxid), Y-axis = `target_percentile` (predicted cutoff). Diagonal dashed reference line.
