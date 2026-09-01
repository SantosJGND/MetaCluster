# Model Evaluation Output Description

All outputs are written under `{analysis_output_filepath}` (set via `--analysis_output_filepath`).

---

## Directory Structure

```
{analysis_output_filepath}/
├── models/
│   ├── recall_xgb_bundle.pkl          (or variant: cutoff_recall_bundle.pkl, direct_xgb_bundle.pkl, recall_gp_clf_pipeline.pkl)
│   ├── composition_xgb_bundle.pkl     (or variant: composition_{rf,gb,lr,optuna}_bundle.pkl)
│   ├── cross_hit_xgb_bundle.pkl
│   ├── taxids_to_use.parquet
│   ├── recall_model_summary.png
│   ├── recall_landscape.png           (gp_clf only)
│   ├── recall_calibration.png         (gp_clf only)
│   ├── recall_actual_at_index.png     (gp_clf only)
│   └── cache/
│       ├── training_results_cache.parquet
│       ├── prediction_results_cache.parquet
│       ├── recall_results_cache.parquet
│       ├── taxids_to_use_cache.parquet
│       └── recall_matrices_cache.joblib
│
├── cross_hit_metrics_per_dataset.tsv
├── cross_hit_metrics_stats.tsv
├── spurious_hit_metrics_per_dataset.tsv
├── spurious_hit_metrics_stats.tsv
├── training_stats_matrices.tsv          (--no-cache only)
├── evaluation_results.json
├── pipeline_metadata.tsv
├── test_datasets_input_df.tsv
├── test_datasets_overall_precision.tsv
├── test_datasets_summary_results.tsv
├── test_datasets_spurious_composition.tsv
├── test_datasets_cross_hit_composition.tsv
│
├── precision_clade_post_histogram.png
├── precision_metrics_boxplot.png
├── precision_metrics_histogram.png
├── recall_metrics_boxplot.png
├── recall_improvement_histogram.png
├── probability_metrics_boxplot.png
├── cross_hit_composition_heatmap.png
├── spurious_composition_heatmap.png
├── filtering_benefit.png
├── recall_error_boxplot.png
├── recall_rmse_distribution.png
├── last_best_match_vs_rmse.png
├── cutoff_error_histogram.png
├── cutoff_confusion_matrix.png
├── cross_hit_metrics_boxplot.png
├── cross_hit_distribution_histogram.png
├── cross_hit_improvement_scatter.png
├── mutation_rate_vs_crosshits.png
├── cluster_size_distribution.png
├── cross_hit_{metric}_distribution.png       (4 boxplots)
├── cross_hit_dotplot_reads_vs_count.png
├── cross_hit_dotplot_reads_vs_mapped.png
├── spurious_hit_{metric}_distribution.png    (3 boxplots)
├── spurious_hit_dotplot_reads_vs_count.png
├── spurious_hit_dotplot_reads_vs_mapped.png
│
└── evaluation_report.html
```

---

## Models Directory

### `models/*.bundle.pkl` — Trained Model Bundles (Format A)

Serialized model bundles saved via each modeller's `save_model()` method. Files are `joblib.dump()` of a Python dict (not a raw model object). See [Serialization Formats](#serialization-formats) for details.

| `--recall_model_interface` | Saved file | Class |
|---|---|---|
| `xgb`, `morf`, `moxgb_optimized`, `morf_optimized`, `monn_optimized` | `recall_xgb_bundle.pkl` | `RecallModeller` — multi-output regressor predicting full recall curve |
| `direct` | `cutoff_recall_bundle.pkl` | `CutoffRecallModeller` — RF classifier predicting k_min directly |
| `direct_xgb` | `direct_xgb_bundle.pkl` | `DirectXGBRecallModeller` — XGBoost regressor predicting tau-crossing fraction directly |
| `gp_clf` | `recall_gp_clf_pipeline.pkl` | `GPCLFRecallModeller` — per-division GP regressors + CLF threshold |

| `--composition_model_interface` | Saved file | Class |
|---|---|---|
| `xgb` | `composition_xgb_bundle.pkl` | `XGBCompositionModeller` |
| `xgb_optimized` | `composition_optuna_bundle.pkl` | `OptunaXGBCompositionModeller` |
| `rf` | `composition_rf_bundle.pkl` | `RFCompositionModeller` |
| `gb` | `composition_gb_bundle.pkl` | `GBCompositionModeller` |
| `lr` | `composition_lr_bundle.pkl` | `LRCompositionModeller` |

| Model | Saved file | Class |
|---|---|---|
| Cross-hit (always trained) | `cross_hit_xgb_bundle.pkl` | `CrossHitModeller` |

### `models/taxids_to_use.parquet`

The `taxids_to_use` DataFrame at training time — columns: `taxid`, `order`, `family`, `genus`.

### Serialization Formats

Two persistence formats exist, serving different use cases:

| Aspect | Format A — Dict Bundle | Format B — Full Object Dump |
|--------|------------------------|-----------------------------|
| **Extension** | `.pkl` | `.joblib` |
| **Saved via** | Each modeller's `save_model()` | `ModelRegistry.save_models()` |
| **Content** | `joblib.dump(dict)` with specific keys | `joblib.dump(modeller_object)` — entire instance |
| **Files** | `recall_xgb_bundle.pkl`, `composition_xgb_bundle.pkl`, `cross_hit_xgb_bundle.pkl`, etc. | `recall_model.joblib`, `composition_model.joblib`, `crosshit_model.joblib` |
| **Used by** | `ModelTrainer.save_models()` → evaluation pipeline; `ml_api` local file loading | `ModelRegistry.load_models()` → custom deployments |

**Format A dict keys (varies by modeller):**

| Modeller | Bundle keys |
|---|---|
| `RecallModeller` / `CutoffRecallModeller` / `DirectXGBRecallModeller` | `model_type`, `model`, `feature_names`, `target_names`, `data_set_divide`, `transformer`, (optional) `target_recall` |
| `GPCLFRecallModeller` | Above + `pipeline`, `optimal_tau`, `optimal_X`, `optimal_loss` |
| `BaseCompositionModeller` subclasses | `model_type`, `pipeline`, `feature_names`, `X_train`, `X_test`, `y_train`, `y_test` |
| `CrossHitModeller` | `model`, `scaler`, `pca` |

**Important:** The two formats are **not interchangeable**. The `ml_api` (FastAPI inference server) expects Format A dict bundles with a `"transformer"` key at the top level. Format B files cannot be loaded by the API.

### `models/cache/` — Cached Training Data

Written on first run (or when `--no-cache` is absent). Read on subsequent runs to skip re-processing.

| File | Format | Contents |
|------|--------|----------|
| `training_results_cache.parquet` | Parquet | Merged output of `data_set_traversal_with_precision()` — per-dataset traversal results with columns: `data_set`, `node`, `best_match_taxid`, `precision`, `stop_traversal`, `n_true_leaves`, plus one column per tax-level class |
| `prediction_results_cache.parquet` | Parquet | Merged output of `cross_hit_prediction_matrix()` — training features for cross-hit model |
| `recall_results_cache.parquet` | Parquet | Merged output of `predict_recall_cutoff_vars()` — recall feature vectors with bin counts per division |
| `taxids_to_use_cache.parquet` | Parquet | Snapshot of `taxids_to_use` at training time |
| `recall_matrices_cache.joblib` | joblib | `List[pd.DataFrame]` — one raw `m_stats_matrix` per training dataset |

---

## Training Data Spurious / Cross-Hit Analysis

Generated by `evaluate.py:285-286` → `analyze_cross_hit_distribution()` (lines 18-82) and `analyze_spurious_hit_distribution()` (lines 85-146). These process `DatasetResult` objects from the training data to produce per-class metrics and plots.

### `cross_hit_metrics_per_dataset.tsv`

One row per (class, dataset) combination. Columns:

| Column | Type | Description |
|--------|------|-------------|
| `class` | str | Taxonomic class at the configured `--tax_level` |
| `data_set` | str | Dataset name |
| `reads_simulated` | int | Reads simulated for this class |
| `cross_hit_count` | int | Number of cross-hit leaves matching this class |
| `cross_hit_reads_mapped` | int | Total reads mapped to cross-hit leaves in this class |
| `ratio` | float | `cross_hit_reads_mapped / reads_simulated` |

### `cross_hit_metrics_stats.tsv`

Aggregated by class across all datasets. Columns: `class`, then for each metric `(reads_simulated, cross_hit_count, cross_hit_reads_mapped, ratio)`: `_{mean, median, std, min, max}`, plus `n_simulations` (count).

### `spurious_hit_metrics_per_dataset.tsv`

Same structure as cross-hit but with `spurious_hit_count`, `spurious_hit_reads_mapped`, `ratio`.

### `spurious_hit_metrics_stats.tsv`

Same aggregation as cross-hit stats for spurious columns.

---

## Evaluation Results (TSV)

Generated by `evaluate.py:302` → `results.save_tsv()` (`result_models.py:211-217`) and `results.save_metadata()` (`result_models.py:340-356`).

### `pipeline_metadata.tsv`

Pipeline run summary. Two columns: `metric`, `value`.

| Metric | Description |
|--------|-------------|
| `total_datasets` | Total datasets attempted (successful + failed + skipped) |
| `successful` | Datasets that produced a `DatasetResult` |
| `failed` | Datasets that raised an exception |
| `skipped_count` | Datasets silently skipped (no mapped reads) |
| `skipped_datasets` | Semicolon-separated list of skipped dataset names (omitted if none) |
| `failed_datasets` | Semicolon-separated list of error messages (omitted if none) |

The same two-column `metric`/`value` schema is emitted to
`evaluation_results.json`→`metadata` (`skipped_count`, `failed_datasets` included),
so the TSV and JSON serializations always agree.

### EDA Cohort Metadata — `analysis_data_extractor.py`

The EDA producer (`deployment/model_evaluation/analysis_data_extractor.py`) writes a
cohort-scoped `pipeline_metadata.tsv` in the same output directory as
`per_dataset_metrics.tsv` (e.g. `<domain>_eda/pipeline_metadata.tsv`). This logs the
full study cohort, so downstream statistics can state their denominator honestly.

| Metric | Description |
|--------|-------------|
| `total_attempted` | All dataset folders scanned from `--study-output` |
| `extracted` | Rows written to `per_dataset_metrics.tsv` (= successful) |
| `dropped` | `total_attempted − extracted` |
| `failed` | Datasets that raised an exception during extraction |
| `skipped` | Datasets with no usable m-stats / no mapped reads (recall 0) |
| `skipped_datasets` | Semicolon-separated names (omitted if none) |
| `failed_datasets` | Semicolon-separated error messages (omitted if none) |
| `study_gaps` | Folders missing required files (input/clustering/matched); not attempted |
| `study_gap_datasets` | Semicolon-separated names (omitted if none) |

Invariant: `dropped == skipped + failed + study_gaps`.

### `test_datasets_input_df.tsv`

Merged input summary for all test datasets. Columns: `sample`, `taxid`, `reads`, `mutation_rate`, `order`, `family`, `genus`, `data_set`, `found_in_recall_filter`, `found_in_fixed_filter`.

### `test_datasets_overall_precision.tsv`

One row per test dataset.

| Column | Type | Description |
|--------|------|-------------|
| `recall_baseline` | float | Raw recall |
| `precision_fixed` | float | Clade precision after fixed filter (max 12 taxids, distance threshold 0.6) |
| `precision_clade_post` | float | Clade precision after all cleanup (recall filter + cross-hit cleanup + composition prediction) |
| `data_set` | str | Dataset name |

### `test_datasets_summary_results.tsv`

One row per test dataset with ~39+ columns:

**Precision metrics:**
- `purity` — fraction of non-trash leaves
- `purity_cov_filtered` — non-trash leaves with coverage > 0
- `precision_best_match` — best-match taxids that are in input / total best-match taxids
- `precision_clade_full` — clade precision before cross-hit cleanup
- `precision_clade_post_cleanup` — clade precision after cross-hit cleanup
- `precision_clade_fixed` — clade precision after fixed filter (min_dist=0.6, max 12 taxids)

**Recall metrics:**
- `recall_baseline` — |output ∩ input| / |input|
- `recall_baseline_cov_filtered` — after coverage > 0 filter
- `recall_clade_pre_cleanup` — clade recall before cross-hit cleanup
- `recall_clade_post_cleanup` — clade recall after cross-hit cleanup
- `recall_after_recall_filter` — recall after applying the predicted recall cutoff
- `recall_fixed_max_12` — recall after fixed filter

**Cross-hit metrics:**
- `predicted_cross_hits` — number of cross-hits detected by model
- `cross_hit_specificity` — TP / (TP + FP)
- `cross_hit_precision` — same as specificity
- `cross_hit_recall` — TP / total_true_cross_hits
- `cross_hit_f1` — harmonic mean of cross_hit_precision and cross_hit_recall
- `total_true_cross_hits` — leaves where best_match_is_best is False
- `total_cross_hit_reads_mapped` — reads mapped to cross-hit leaves
- `cross_hit_counts_per_class` — dict of {class: count}
- `cross_hit_reads_per_class` — dict of {class: reads}

**Recall model diagnostics** (present when recall_metrics dict is populated):
- `recall_metric_target_recall` — requested recall target
- `recall_metric_target_percentile` — predicted percentile
- `recall_metric_keep_index` — predicted keep index
- `recall_metric_predicted_k_min` — predicted minimum bin count
- `recall_metric_actual_k_min` — actual minimum bin count
- `recall_metric_cutoff_error` — predicted - actual k_min
- `recall_metric_last_best_match_relindex` — relative index of last best match
- `recall_metric_per_division_recall_rmse` — RMSE across divisions
- `recall_metric_error_index_recall_{i}` — per-division recall error

### `test_datasets_spurious_composition.tsv`

Per-dataset spurious (trash) composition — columns per taxonomic class at the configured tax level, with `data_set` and `tax_level`.

### `test_datasets_cross_hit_composition.tsv`

Same structure for cross-hit composition.

---

## JSON Output

### `evaluation_results.json`

Generated by `evaluate.py:303` → `results.to_json()` (`result_models.py:184-195`).

```json
{
  "generated_at": "2026-06-25T12:00:00",
  "metadata": {
    "total_datasets": 50,
    "successful": 48,
    "failed": 2,
    "errors": ["..."]

  },
  "test_results": [
    {"precision_fixed": 0.85, "precision_clade_post": 0.92, "data_set": "dataset_001"}
  ],
  "summary_results": [
    {/* all columns from the summary TSV */}
  ],
  "spurious_composition": [{/* per-dataset spurious class fractions */}],
  "cross_hit_composition": [{/* per-dataset cross-hit class fractions */}]
}
```

---

## PNG Images

All generated by `evaluate.py:307` → `ResultVisualizer.plot_all()` (`visualization.py:32-55`), plus training analysis plots from `analyze_cross_hit_distribution` and `analyze_spurious_hit_distribution`. Dimensions ~10×6 to 12×8 inches, matplotlib tight_layout, default 100dpi (unless noted).

### Precision Distribution

**`precision_clade_post_histogram.png`** — Histogram with KDE overlay of `precision_clade_post` across all test datasets. xlim [0, 2], 20 bins.

### Precision Comparison

| Image | Description |
|-------|-------------|
| `precision_metrics_boxplot.png` | Side-by-side boxplot of 6 precision metrics: `precision_best_match`, `purity`, `purity_cov_filtered`, `precision_clade_full`, `precision_clade_post_cleanup`, `precision_clade_fixed`. ylim [0, 3]. |
| `precision_metrics_histogram.png` | Overlaid histograms (hue=Metric) of all precision metrics, 30 bins, xlim [0, 3]. |

### Recall Comparison

| Image | Description |
|-------|-------------|
| `recall_metrics_boxplot.png` | Boxplot of 6 recall metrics: `baseline`, `baseline_cov_filtered`, `clade_pre_cleanup`, `clade_post_cleanup`, `after_recall_filter`, `fixed_max_12`. |

### Recall Improvement

**`recall_improvement_histogram.png`** — Boxplot of 4 difference metrics (despite the name):
- `recall_clade_diff` = clade_pre_cleanup - baseline
- `recall_clade_diff_clean` = clade_post_cleanup - baseline
- `recall_clade_diff_predicted_leaves` = after_recall_filter - baseline
- `recall_method_comparison` = fixed_max_12 - after_recall_filter

### Probability Metrics

**`probability_metrics_boxplot.png`** — Violin plot (despite name) of 4 probability products:
- `Prob_Find_any` = recall_raw × purity
- `Prob_Find_true` = recall_raw × precision_best_match
- `Prob_Find_true_clade_full` = clade_recall_pre × clade_precision_full
- `Prob_Find_true_clade_clean` = clade_recall_post × clade_precision_post

### Composition Heatmaps

| Image | Description |
|-------|-------------|
| `cross_hit_composition_heatmap.png` | Heatmap of mean cross-hit composition fraction per `tax_level` row. Rows sum to 1. viridis colormap. |
| `spurious_composition_heatmap.png` | Same for spurious/trash composition. |

### Filtering Benefit

**`filtering_benefit.png`** — Horizontal stacked bar chart. One bar per sample. Blue = coverage retained above cutoff, coral = coverage lost. xlim [0, 1].

### Recall Model Diagnostics

These plots are only generated when the recall model produces `recall_metric_*` columns:

| Image | Description |
|-------|-------------|
| `recall_error_boxplot.png` | Per-division recall prediction error (predicted - true), red dashed line at 0. |
| `recall_rmse_distribution.png` | Histogram with KDE of per-division recall RMSE across datasets. |
| `last_best_match_vs_rmse.png` | Scatter of last_best_match_relindex vs per-division RMSE, colored by `recall_after_recall_filter`. |
| `cutoff_error_histogram.png` | Histogram with KDE of `cutoff_error` (predicted_k_min - actual_k_min). |
| `cutoff_confusion_matrix.png` | Heatmap crosstab of predicted vs actual k_min bins. Blues colormap. |

### Cross-Hit Plots

Only generated when `cross_hit_precision` column exists in summary results:

| Image | Description |
|-------|-------------|
| `cross_hit_metrics_boxplot.png` | Boxplot of precision, recall, F1, specificity. ylim [0, 1.1]. 300dpi. |
| `cross_hit_distribution_histogram.png` | Two side-by-side histograms: predicted cross-hits (steelblue) and true cross-hits (coral). |
| `cross_hit_improvement_scatter.png` | Scatter of predicted vs true cross-hits with red dashed perfect-prediction line. 300dpi. |

### Additional Plots

| Image | Description |
|-------|-------------|
| `mutation_rate_vs_crosshits.png` | Scatter of mutation_rate vs cross_hit_count from merged input+cross-hit data, point size proportional to count. |
| `cluster_size_distribution.png` | Two panels: histogram of cluster size distribution (coral) + barplot of top 20 largest clusters (viridis). |

### Training Cross-Hit Analysis Plots

Generated by `analyze_cross_hit_distribution()`:

| Image | Description |
|-------|-------------|
| `cross_hit_reads_simulated_distribution.png` | Boxplot per class of reads simulated. |
| `cross_hit_cross_hit_count_distribution.png` | Boxplot per class of cross-hit count. |
| `cross_hit_cross_hit_reads_mapped_distribution.png` | Boxplot per class of cross-hit reads mapped. |
| `cross_hit_ratio_distribution.png` | Boxplot per class of cross-hit reads / simulated reads ratio. |
| `cross_hit_dotplot_reads_vs_count.png` | Log-log scatter of reads_simulated vs cross_hit_count. |
| `cross_hit_dotplot_reads_vs_mapped.png` | Log-log scatter of reads_simulated vs cross_hit_reads_mapped. |

### Training Spurious-Hit Analysis Plots

Generated by `analyze_spurious_hit_distribution()`:

| Image | Description |
|-------|-------------|
| `spurious_hit_spurious_hit_count_distribution.png` | Boxplot per class of spurious-hit count (salmon). |
| `spurious_hit_spurious_hit_reads_mapped_distribution.png` | Boxplot per class of spurious-hit reads mapped. |
| `spurious_hit_ratio_distribution.png` | Boxplot per class of spurious-hit reads / simulated reads ratio. |
| `spurious_hit_dotplot_reads_vs_count.png` | Log-log scatter of reads_simulated vs spurious_hit_count. |
| `spurious_hit_dotplot_reads_vs_mapped.png` | Log-log scatter of reads_simulated vs spurious_hit_reads_mapped. |

---

## HTML Report

### `evaluation_report.html`

Generated by `evaluate.py:311` → `generate_report()` (`visualization.py:995-1036`).

Structure:
- **Summary Metrics** section — metric cards showing mean/median/std/min/max for each precision column
- **Visualizations** section — embedded `<img>` tags for each plot PNG that exists on disk. Referenced plots: `precision_clade_post_histogram.png`, `precision_metrics_boxplot.png`, `recall_metrics_boxplot.png`, `recall_improvement_histogram.png`, `probability_metrics_boxplot.png`, `cross_hit_composition_heatmap.png`, `spurious_composition_heatmap.png`, `recall_error_boxplot.png`, `recall_error_trend.png`, `mutation_rate_vs_crosshits.png`, `cluster_size_distribution.png`
- **Footer** with pipeline name

---

## Pipeline Data Flow

```
study_output/                      (input)
├── dataset_001/
│   ├── input/dataset_001.tsv
│   ├── output/clade_report_with_references.tsv
│   └── clustering/               (OverlapManager data)
├── dataset_002/
│   └── ...
└── training_stats_matrices.tsv   (only with --no-cache)

  ↓ DataLoader
  ↓ ModelTrainer (cached to models/cache/*.parquet)
  ↓ evaluate.py::main()

analysis_output/                   (output)
├── models/*.joblib               trained model artifacts
├── models/cache/*                cached training data
├── cross_hit_metrics_*.tsv       training analysis
├── spurious_hit_metrics_*.tsv    training analysis
├── test_datasets_*.tsv           evaluation results
├── evaluation_results.json       aggregated results
├── *.png                         all plots
└── evaluation_report.html        HTML report
```
