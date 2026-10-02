# Model Evaluation Output Description

All outputs are written under `{analysis_output_filepath}` (set via `--analysis_output_filepath`).

---

## Directory Structure

```
{analysis_output_filepath}/
├── models/
│   ├── recall_xgb_bundle.pkl          (or variant: cutoff_recall_bundle.pkl, direct_xgb_bundle.pkl, recall_gp_clf_pipeline.pkl)
│   ├── composition_xgb_bundle.pkl     (or variant: composition_{rf,gb,lr,optuna}_bundle.pkl)
│   ├── composition_feature_importances.tsv
│   ├── recall_feature_importances.tsv (not for gp_clf / monn_optimized — no native importances)
│   ├── cross_hit_xgb_bundle.pkl       (--enable-cross-hit only)
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
├── cross_hit_metrics_per_dataset.tsv   (--enable-cross-hit only)
├── cross_hit_metrics_stats.tsv         (--enable-cross-hit only)
├── spurious_hit_metrics_per_dataset.tsv
├── spurious_hit_metrics_stats.tsv
├── training_stats_matrices.tsv          (--no-cache only)
├── evaluation_results.json
├── evaluation_results_agent.json
├── pipeline_metadata.tsv
├── test_datasets_input_df.tsv
├── test_datasets_reference_leaves.tsv
├── test_datasets_overall_precision.tsv
├── test_datasets_summary_results.tsv
├── test_datasets_spurious_composition.tsv
├── test_datasets_cross_hit_composition.tsv   (--enable-cross-hit only)
├── precision_summary_statistics.tsv
├── recall_summary_statistics.tsv
├── recall_complete_summary_statistics.tsv    (only when some datasets fail the gate)
├── cross_hit_summary_statistics.tsv          (--enable-cross-hit only)
│
├── precision_clade_post_histogram.png
├── precision_metrics_boxplot.png
├── precision_metrics_histograms.png
├── recall_metrics_boxplot.png
├── recall_improvement_histogram.png
├── probability_metrics_boxplot.png
├── cross_hit_composition_heatmap.png   (--enable-cross-hit only)
├── spurious_composition_heatmap.png
├── filtering_benefit.png
├── recall_error_boxplot.png
├── recall_error_trend.png
├── recall_rmse_distribution.png
├── last_best_match_vs_rmse.png
├── cutoff_error_histogram.png
├── cutoff_confusion_matrix.png
├── cross_hit_metrics_boxplot.png            (--enable-cross-hit only)
├── cross_hit_distribution_histogram.png     (--enable-cross-hit only)
├── cross_hit_improvement_scatter.png        (--enable-cross-hit only)
├── mutation_rate_vs_crosshits.png
├── cluster_size_distribution.png
├── cross_hit_{metric}_distribution.png       (4 boxplots; --enable-cross-hit only)
├── cross_hit_dotplot_reads_vs_count.png      (--enable-cross-hit only)
├── cross_hit_dotplot_reads_vs_mapped.png     (--enable-cross-hit only)
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
| Cross-hit (`--enable-cross-hit` only) | `cross_hit_xgb_bundle.pkl` | `CrossHitModeller` |

> Cross-hit modelling is optional: with the default `--no-enable-cross-hit` the
> model is **not trained**, the bundle is not written, and every cross-hit-specific
> output below is skipped. The skip is logged in `evaluate.log` and recorded as a
> `cross_hit_enabled` row in `pipeline_metadata.tsv` / `evaluation_results.json`
> metadata.

### Composition Model Feature Importances

The composition model is trained during evaluation (`evaluate.py:299-301` → `trainer.train_models()` builds `composition_modeller`, `models.py:385,421`), then `trainer.evaluate_models(config.models_dir)` (`models.py:482-489`) calls `composition_modeller.eval_and_plot(X_test, y_test, output_dir, X_train=...)` (`om_models.py:2396-2406`). The feature-importance plots and TSV are written into `config.models_dir` (`{analysis_output_filepath}/models/`).

| Composition modeller | Output file | Source attribute |
|---|---|---|
| `XGBCompositionModeller` (`xgb`) | `composition_feature_importance.png` (top-20 horizontal bar, `om_models.py:2300-2314`) | `classifier.feature_importances_` |
| `RFCompositionModeller` (`rf`) | `composition_feature_importance.png` | `classifier.feature_importances_` |
| `GBCompositionModeller` (`gb`) | `composition_feature_importance.png` | `classifier.feature_importances_` |
| `OptunaXGBCompositionModeller` (`xgb_optimized`) | `composition_feature_importance.png` | `classifier.feature_importances_` |
| `LRCompositionModeller` (`lr`) | `lr_coefficients.png` (absolute coefficient bar, `om_models.py:2651-2667`) | `classifier.coef_` |

All composition backends additionally write **`composition_feature_importances.tsv`** via
`save_feature_importances()` (`om_models.py:2465-2475`), called from `eval_and_plot()`
(`om_models.py:2403`). It is a two-column TSV (feature name, importance magnitude),
sorted descending — tree backends use `feature_importances_`, the linear backend the
absolute `coef_`, resolved through `feature_importance_dataframe()` (`om_models.py:2437-2461`)
with feature names aligned via `feature_names_in_` / preprocessor output names
(`om_models.py:2418-2435`).

**Notes:**
- `eval_and_plot` also attempts SHAP plots (`shap_summary_plot.png`, `shap_bar_plot.png`, `shap_dependence_plot.png`, plus interaction plots) via `shap_eval_plot()` / `shap_interaction_plot()` (`om_models.py:2316-2394`). These are wrapped in a silent `try/except pass`, so they will **not** appear if `shap` is not installed, lacks the required dependence, or errors — they are best-effort only. The `LRCompositionModeller` overrides both SHAP methods to no-ops.
- Tree-based backends (`feature_importances_`) and the linear backend (`coef_`) are mutually exclusive: only the plot matching the chosen `--composition_model_interface` is emitted.

### Recall Model Feature Importances

After recall training, `trainer.evaluate_models(config.models_dir)` calls
`self.recall_modeller.save_feature_importances(output_dir)` (`models.py:478-479`).
This writes **`models/recall_feature_importances.tsv`** — a sorted, two-column TSV
(feature name, importance magnitude) of the **aggregated per-feature importances**
(feature means across all targets), via `RecallModeller.save_feature_importances()`
(`om_models.py:718-728`).

| `--recall_model_interface` | Importance source |
|---|---|
| `xgb`, `morf`, `moxgb_optimized`, `morf_optimized` | Mean of `est.feature_importances_` over the `MultiOutputRegressor` estimators (`om_models.py:705-717`) |
| `direct` | `model.feature_importances_` from the RF classifier (`om_models.py:981-988`) |
| `direct_xgb` | `model.feature_importances_` from the single XGBoost regressor (`om_models.py:1144-1151`) |
| `gp_clf`, `monn_optimized` | None — per-division GP regressors and the MLP backend (`MultiOutputRegressor(MLPRegressor(...))`, `om_models.py:430-431`) have no native importances; the TSV is **skipped** (logged only) |

The base `RecallModeller.model_summary` no longer emits the legacy
`recall_model_analysis_results_feature_importances.tsv`.

### `models/taxids_to_use.parquet`

The `taxids_to_use` DataFrame at training time — columns: `taxid`, `order`, `family`, `genus`.

> **Important — not to be confused with the input panel (taxid plan).**
> This table is *not* the `taxid_plan` / `bacterial_assess.tsv` input panel. The input panel lists the reference taxa configured for the simulation (loaded from `--taxid-plan-filepath`). `taxids_to_use.parquet` is a **derived artifact** computed from the study's *output* (the clade reports of the simulated datasets) — it records which taxids actually appear in the analysis outputs, standardised at a chosen taxonomic level. They are distinct objects with different sources.

**Construction** (`data_loader.py`):

1. **Entry point** — `evaluate.py:274` calls `DataLoader(config).initialize()`, which runs `_establish_taxids()` last in its loading sequence (`data_loader.py:269-281`).
2. **Gather output taxids** — `output_parse()` (`data_loader.py:75-110`) scans **every** dataset folder's `output/clade_report_with_references.tsv`, unions all output taxids found, and adds the input taxids from `all_input_data["taxid"].dropna().unique()`.
3. **Resolve lineages** — `ncbi_wrapper.resolve_lineages()` is called, then `get_level()` fills the `order`, `family`, `genus` columns from the NCBI taxonomy DB for each taxid.
4. **Filter by frequency** — `establish_taxids_to_use()` (`data_loader.py:113-159`) drops rows whose tax-level value-count frequency `<= min_tax_count` (i.e. `config.taxa_threshold`, default `0.02`).
5. **Clean + sentinel** — drops NaNs on the chosen tax level, dedupes on it, then appends an `unclassified` sentinel row with `taxid=0`.
6. **Persist** — `evaluate.py:300` → `trainer.save_models(config.models_dir)` → `models.py:458-461` writes this DataFrame to `{analysis_output_filepath}/models/taxids_to_use.parquet`. A copy is also cached as `models/cache/taxids_to_use_cache.parquet` (see below).

**Consumers**: `loader.get_taxids_to_use()` is passed to `ModelTrainer`, `TrainingCrossHitAnalyzer`, and `BatchEvaluator` (`evaluate.py:295, 308, 321`); it is reloaded by `models.py:641-644`, `ml_api/train_recall.py:41`, and `deployment/analysis/analyze_samples.py:521,534`.

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

Generated by `evaluate.py:285-286` → `analyze_cross_hit_distribution()` (lines 18-82) and `analyze_spurious_hit_distribution()` (lines 85-146). These process `DatasetResult` objects from the training data to produce per-class metrics and plots. The cross-hit portion (TSVs + plots) runs only with `--enable-cross-hit`; the spurious-hit portion always runs.

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
| `cross_hit_enabled` | `true`/`false` — whether cross-hit modelling+filtering ran (`--enable-cross-hit`). Present so a skipped run is unambiguous when reading the results |

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
| `recall_analysis_datasets` | Datasets passing the assembly-completeness gate (used for recall-facing analyses) |
| `assembly_incomplete_excluded` | Datasets skipped from recall-facing analyses by the gate (kept in precision/composition outputs) |

Each `per_dataset_metrics.tsv` row also carries an `assembly_complete` flag marking gate status.

Invariant: `dropped == skipped + failed + study_gaps`.

### Recall definitions (shared by both pipelines)

Recall is always `|detected ∩ input_taxids| / |input_taxids|`; the variants differ only in
what counts as *detected*. The same three definitions are emitted by the evaluation
pipeline (`test_datasets_summary_results.tsv`) and by the EDA extractor
(`per_dataset_metrics.tsv`, `recall_data.tsv`):

| Definition | Detected = | Evaluation column | Extractor per-dataset | Extractor per-taxid |
|---|---|---|---|---|
| **Classification-aware** (headline) | best-matched, non-trash leaf in the m-stats matrix **or** classifier row with `uniq_reads > 0` in `<ds>_merged_classification.tsv` | `recall_baseline` | `overall_recall` | `recalled` (GEE target) |
| **Assembly** | best-matched, non-trash leaf only (zero-coverage leaves included) | `recall_baseline_assembly` | `recall_assembly` | `recalled_assembly` |
| **Assembly, coverage-filtered** | best-matched leaf with `coverage > 0` | `recall_baseline_cov_filtered` | `recall_cov_filtered` | — |
| **Classification only** | classifier row with `uniq_reads > 0` | `recall_baseline_classification` | `recall_classification` | `recalled_classification` |

Notes:
- Post-filter metrics (`recall_after_recall_filter`, `recall_fixed_max_12`, `recall_clade_*`)
  measure what survives the filtering/clustering stages and therefore stay **assembly-based**.
- The classification-aware definition credits taxa the classifier found but for which no
  reference genome was matched (`classified_no_assembly` in the diagnostic below); the
  assembly definitions measure the full classify → map → cluster chain.
- Datasets without a classification file fall back to assembly-based values for the
  headline column (the classifier set is empty).
- Historical note: before this change `recall_baseline` only counted best-matched leaves and
  the best-marking step skipped zero-coverage leaves, so it was identical to
  `recall_baseline_cov_filtered`. Values regenerated with the current code are therefore not
  comparable to older summary files without recomputation.

#### Known issue — `max_taxids_fixed_filter` is not wired up

`recall_fixed_max_12` and `found_in_fixed_filter` come from a **hard-coded** limit of 12.
The `--max_taxids_fixed_filter` CLI flag and `EvaluatorConfig.max_taxids_fixed_filter`
(default 15) are accepted but never read: `DatasetProcessor._apply_fixed_filter` is called
with a literal `12`.

Two further caveats about what the limit actually means:

- The value truncates **assembly/leaf rows**, not unique taxids, and it is applied twice
  (`OverlapManager(max_taxids=...)` and an explicit `head(max_taxids)`), so the effective
  row count depends on how many leaves the top assemblies contribute.
- Recall after it is therefore comparable across datasets only in the sense that every
  dataset used the same constant.

Behaviour is left unchanged pending a decision on how the flag should be defined.

### Recall-gap diagnostic — `analysis_scripts/diagnose_recall_gap.py`

File-based diagnostic (no lineage lookups, no network) that assigns every
`(dataset, input taxid)` to exactly one loss bucket:

| Bucket | Meaning |
|---|---|
| `absent_classification` | taxid never appears in `<ds>_merged_classification.tsv` (classifier DB gap / unmappable k-mers) |
| `classified_no_assembly` | classified, but no row in `output/matched_assemblies.tsv` (taxonomy → genome lookup failed; reference-store gap or accession/filename mismatch) |
| `assembly_no_coverage` | matched to an assembly that has no row in `merged_coverage_statistics.tsv` (never mapped) |
| `assembly_zero_reads` | coverage row exists but `numreads == 0` (mapping failed, e.g. high mutation rate) |
| `assembly_with_reads` | reads mapped to the matched assembly — recallable by the m-stats pipeline; any further loss is score/level gating |

Outputs (in `--output-dir`):

| File | Description |
|---|---|
| `recall_gap_per_taxid.tsv` | One row per `(data_set, taxid)`: `reads_simulated`, `mutation_rate`, `bucket`, `in_classification`, `classifier_uniq_reads`, `matched_accession`, `mapped_reads`, plus `has_*_file` flags for missing inputs |
| `recall_gap_per_dataset.tsv` | Bucket `count` / `share` per dataset |
| `recall_gap_summary.tsv` | Bucket `count` / `share` over all input taxa |
| `recall_gap_buckets.png` | Overall bucket shares and stacked shares by mutation-rate bin (skipped with `--no-plot`) |

Usage: `python -m deployment.model_evaluation.analysis_scripts.diagnose_recall_gap --study-output <dir> --output-dir <dir> [--limit N] [--no-plot]`.

### `test_datasets_input_df.tsv`

One row per `(data_set, input taxid)` — the wide per-reference table. It answers
*where did this reference go?* and is the input to `test_datasets_reference_leaves.tsv`'s
funnel.

The first 12 columns are the historical schema, unchanged and in their original
order, so existing readers that select by name or by leading position keep working:

| # | Column | Meaning |
|---|---|---|
| 1 | `sample` | Simulated sample name |
| 2 | `taxid` | Input reference taxid (join key) |
| 3 | `reads` | Reads simulated for this reference |
| 4 | `mutation_rate` | Simulated mutation rate |
| 5–7 | `order`, `family`, `genus` | Resolved taxonomy |
| 8 | `match_index` | Best-matching assembly accession |
| 9 | `match_coverage` | Coverage of that assembly |
| 10 | `found_in_recall_filter` | Survived the recall filter |
| 11 | `found_in_fixed_filter` | Survived the fixed-`max_taxids` filter |
| 12 | `data_set` | Dataset name |

Everything after column 12 is new tracking detail, in the order given by
`reference_details.REFERENCE_DETAIL_COLUMNS`:

| Group | Columns |
|---|---|
| Classifier evidence | `in_classification`, `classifier_reads_available`, `classifier_uniq_reads`, `classifier_names`, `has_coverage_row` |
| Best-match leaf | `best_accession`, `best_leaf`, `leaf_rank`, `numreads`, `total_uniq_reads`, `covbases`, `meanmapq`, `error_rate`, `best_match_score`, `best_match_level`, `max_coverage`, `max_shared` |
| Leaf aggregates | `n_leaves`, `n_clean_matches`, `n_covered_leaves`, `is_trash_best` |
| Clade membership | `found_in_clade_{pre,post,fixed}`, `clade_{pre,post,fixed}_node`, `clade_{pre,post,fixed}_precision`, `clade_{pre,post,fixed}_n_nodes` |
| Attribution | `loss_stage`, then the independent `recalled_*` / `lost_*` / `survived_post_cleanup` flags |

Notes:
- Everything keys on `m_stats["best_match_taxid"]` — the reference a leaf was matched to —
  which is what `metrics.compute_recall` and `compute_clade_recall` use.
  `output/matched_assemblies.tsv` is **not** the join key: its `taxid` column holds the
  *matched assembly's* taxid, not the input reference's.
- `n_clean_matches` mirrors the recall criterion exactly (≥1 leaf with
  `is_trash == False and best_match_is_best == True`), so
  `n_clean_matches > 0` reproduces `compute_recall`'s recalled-taxid set.
- `classifier_uniq_reads` is **empty**, not `0`, when the classification table has no
  read-count column (older raw-virus studies), so "no count available" is distinguishable
  from "zero reads". `classifier_reads_available` is `False` in exactly that case, because
  `analysis_data_extractor.load_classifier_detected_taxids()` drops its `uniq_reads > 0`
  filter when the column is absent and therefore credits every listed taxid.
  `loss_stage` follows the helper rather than assuming a zero count.
- `clade_*_n_nodes` can exceed 1: a reference whose leaves land in more than one predicted
  clade is reported in each, with `clade_*_node` / `clade_*_precision` describing the
  best-matched leaf's clade.
- Building these tables is best-effort. If it fails, a warning is logged and the row is
  emitted without the new columns; metrics are unaffected.

#### `loss_stage`

An ordered, **first-match** cascade in pipeline order. Each reference lands in exactly
one bucket, so the buckets sum to the input count and form a proper funnel. The names
`absent_classification`, `classified_no_assembly`, `assembly_no_coverage` and
`assembly_zero_reads` are kept identical to `analysis_scripts/diagnose_recall_gap.py`.

`loss_stage` uses the **same** detection criterion as the headline `recall_baseline`
(classification-aware recall): a reference is recalled if it has a clean assembly match
(`n_clean_matches > 0`) **or** the classifier called it. Consequently only
`absent_classification` and `classified_no_assembly` are actual recall losses, and the
table satisfies

```python
recalled_classification_aware == loss_stage not in ("absent_classification", "classified_no_assembly")
```

The remaining buckets describe how far an already-recalled reference got.

| `loss_stage` | Recall loss? | Meaning |
|---|---|---|
| `absent_classification` | yes | taxid is in neither source: absent from `<ds>_merged_classification.tsv` and without a clean leaf |
| `survived_classifier_only` | no | classifier called it (`classifier_uniq_reads > 0`, or listed when no read count exists) but no clean assembly leaf matched |
| `classified_no_assembly` | yes | listed by the classifier with no usable read count, and no clean assembly leaf matched |
| `assembly_no_coverage` | no | clean leaf matched but has no row in `merged_coverage_statistics.tsv` |
| `assembly_zero_reads` | no | coverage row exists but `numreads == 0` |
| `lost_recall_truncation` | no | best leaf's rank fell outside the recall filter's `keep_index` |
| `lost_clade_prediction` | no | not present in any post-cleanup predicted clade |
| `survived_post_cleanup` | no | present in a post-cleanup predicted clade |

`no_best_match` and `best_match_trash` were removed from the cascade because neither can
ever win it: a reference with leaves always has a best match, and once a clean match
exists the best leaf is by definition not trash. Both are still reported, via the
`lost_no_best_match` / `lost_best_match_trash` flags and the `n_clean_matches` /
`is_trash_best` columns.

#### Recall and independent loss flags

`recalled_classification_aware`, `recalled_assembly` and `recalled_classifier_only`
decompose the recall set; `recalled_classification_aware` is
`recalled_assembly or recalled_classifier_only`. `recalled_classifier_only` is also the
`survived_classifier_only` stage.

`lost_before_classification`, `lost_assembly_lookup`, `lost_no_coverage`,
`lost_zero_reads`, `lost_no_best_match`, `lost_best_match_trash`,
`lost_recall_truncation`, `lost_clade_prediction`, `survived_post_cleanup` are **not** a
second copy of `loss_stage`. A reference that fails several conditions at once reports
every applicable flag, so the flags can be used to verify the cascade rather than being a
copy of it. Note `lost_assembly_lookup` means "no *clean* assembly match", which is
independent of recall — a classifier-credited reference with only trash leaves is
`lost_assembly_lookup` yet still recalled.

### `test_datasets_reference_leaves.tsv`

One row per `(data_set, input taxid, matched assembly)` — the long per-leaf table that
explains *why* a reference aggregated to the `loss_stage` above.

Key columns:

| Column | Meaning |
|---|---|
| `data_set` | Dataset name |
| `taxid` | The **input reference** this leaf was matched to (join key to `test_datasets_input_df.tsv`) |
| `assembly_taxid` | The matched assembly's own taxid (m_stats' own `taxid` column) |
| `assembly_accession` | Matched assembly accession |
| `leaf` | Alignment leaf file |
| `leaf_rank` | Position in the `total_uniq_reads`-descending order the filters truncate on |
| `is_trash`, `max_shared` | Per-leaf match flags |
| `coverage`, `covbases`, `numreads`, `meanmapq`, `error_rate` | Mapping statistics |
| `is_input_taxid` | `False` for leaves whose best match is outside the input set |
| `survives_recall_filter`, `survives_fixed_filter` | Rank **and** coverage survival, matching `_apply_recall_filter` / `_apply_fixed_filter` |
| `in_clade_pre`, `in_clade_post`, `in_clade_fixed` | Leaf present in a predicted clade at that stage |

Notes:
- `taxid` is the input taxid, not the assembly's; the assembly's taxid is kept separately
  as `assembly_taxid` because m_stats already carries a `taxid` column of its own.
- Leaves whose best match is outside the input set are dropped, except for leaves with no
  best match at all, which are retained with an empty `taxid` and `is_input_taxid = False`.
- `survives_*_filter` is `False` only when the corresponding filter ran and dropped the
  leaf. With `--no-apply-recall-filter` the rank cutoff is inactive and only the coverage
  condition applies — a disabled filter never reports a leaf as lost.
- The file always carries its header, even when a run produced no leaves.

### `test_datasets_overall_precision.tsv`

One row per test dataset.

| Column | Type | Description |
|--------|------|-------------|
| `recall_baseline` | float | Classification-aware baseline recall (see [Recall definitions](#recall-definitions-shared-by-both-pipelines)) |
| `precision_fixed` | float | Clade precision after fixed filter (max 12 taxids, distance threshold 0.6) |
| `precision_clade_post` | float | Clade precision after all cleanup (post-cleanup tree: recall-filtered and, when `--enable-cross-hit`, cross-hit-cleaned; see the [pre/post caveat](#test_datasets_summary_resultstsv)) |
| `data_set` | str | Dataset name |

### `test_datasets_summary_results.tsv`

One row per test dataset with ~39+ columns. When cross-hit modelling is disabled
(`--no-enable-cross-hit`, the default), the cross-hit metrics below are **absent**
(not zero-filled), which skips the cross-hit plots and statistics.

**Assembly-completeness gate:**
- `assembly_complete` — whether the dataset passed the assembly-completeness gate (at most `--max-missing-refs-pct` of detected references with `uniq_reads >= --min-uniq-reads` are missing from `output/matched_assemblies.tsv`). Present on every row; recall-facing plots and recall/denominator analyses filter on it (`visualization.py:50-51`), while precision/composition plots keep all datasets.

**Precision metrics:**
- `purity` — fraction of non-trash leaves
- `purity_cov_filtered` — non-trash leaves with coverage > 0
- `precision_best_match` — best-match taxids that are in input / total best-match taxids
- `precision_clade_full` — composition-model clade precision on the **full tree** (before recall filter or cross-hit cleanup)
- `precision_clade_post_cleanup` — composition-model clade precision on the **post-cleanup tree**: recall-filtered (top `keep_index` leaves after the predicted recall cutoff) and, when `--enable-cross-hit`, additionally cross-hit-cleaned
- `precision_clade_fixed` — clade precision after the fixed filter (`min_dist=0.6`, max 12 taxids) on the post-cleanup tree

**Recall metrics** (definitions in [Recall definitions](#recall-definitions-shared-by-both-pipelines)):
- `recall_baseline` — classification-aware: best-matched leaves ∪ classifier-detected taxids (`uniq_reads > 0`)
- `recall_baseline_assembly` — assembly-based, best-matched leaves (zero-coverage included)
- `recall_baseline_cov_filtered` — assembly-based, best-matched leaves with coverage > 0
- `recall_baseline_classification` — classifier-only
- `recall_clade_pre_cleanup` — composition-model clade recall on the **full tree** (assembly-based)
- `recall_clade_post_cleanup` — composition-model clade recall on the **post-cleanup tree** (assembly-based): recall-filtered and, when `--enable-cross-hit`, additionally cross-hit-cleaned
- `recall_after_recall_filter` — recall after applying the predicted recall cutoff (assembly-based). Falls back to `0.0` when the recall filter is skipped via `--no-apply-recall-filter`
- `recall_fixed_max_12` — recall after fixed filter (assembly-based)

**Baseline recall-gap decomposition** (present on every row, `batch_evaluator.py`):
- `recall_baseline_cov_gap` — `recall_baseline − recall_baseline_cov_filtered`; the total drop between the headline and coverage-filtered baselines
- `recall_classification_credit` — `recall_baseline − recall_baseline_assembly`; input taxids credited solely by classifier evidence (`uniq_reads > 0`) with no best-matched assembly
- `recall_zero_coverage_loss` — `recall_baseline_assembly − recall_baseline_cov_filtered`; best-matched assemblies whose leaf has `coverage == 0` (mapping-stage loss)

> **Best-match coverage preference** — since the covered-assembly change,
> `_compute_best_matches` (`metagenomics_utils/overlap_manager/node_stats.py`) sorts each
> `best_match_taxid` group by `coverage>0` first, then `best_match_score`/`coverage`/`error_rate`,
> so a **covered** assembly is marked `best_match_is_best` whenever one exists — a taxid whose
> top-scored assembly has `coverage == 0` no longer drops out of `recall_baseline_cov_filtered`
> if a covered assembly is available. `recall_baseline` is unchanged (every group still emits one
> best row). Values computed before this change are not directly comparable because best rows also
> feed precision, `is_trash`/`is_crosshit`, purity, and composition.

> **Pre/post naming caveat** — `precision_clade_post_cleanup` / `recall_clade_post_cleanup`
> are computed on the **post-cleanup tree** (`dataset_processor.py:543`). By default that
> is the recall-filtered tree (truncated `OverlapManager(max_proportion=target_percentile)`),
> so post-cleanup differs from pre-cleanup even when `enable_cross_hit=False`. When
> `--enable-cross-hit`, `_apply_crosshit_cleanup` (`dataset_processor.py:151-152`) mutates
> that same tree in place, so the cleaned tree is what `_predict_clades_postcleanup`
> (`dataset_processor.py:156`) sees — the post-cleanup columns then reflect cross-hit
> removal. With `--no-apply-recall-filter` the post-cleanup tree is the original full OM
> (optionally cross-hit-cleaned), isolating the cross-hit cleanup effect; in that mode the
> recall model is not invoked and `recall_metric_*` diagnostics are absent. Cross-hit
> metrics below are scoped to the kept leaves of the post-cleanup tree.

**Cross-hit metrics** (present only with `--enable-cross-hit`):
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

Same structure for cross-hit composition. Written only with `--enable-cross-hit`.

### Input set vs. `taxids_to_use` — which is used for TP / spurious / cross-hit estimates

TP, cross-hit, spurious, recall and precision estimates are computed **against the full input taxid set of each dataset** (`result.input_df`), **not** against the filtered `taxids_to_use` table:

- **`input_df`** (`dataset_processor.py:155-157`) is built from each dataset's own `input/{dataset}.tsv`, then merged with `input_tax_df` (the expanded lineage table of **all** input taxids) for its `order`/`family`/`genus` columns.
- **Recall** (`metrics.py:270-296`) and **precision** (`metrics.py:240-267`) both take `input_summary` and use `set(input_summary["taxid"].unique())` as the denominator / reference: e.g. `recall = |output ∩ input| / |input|`.
- **Cross-hit TP/FP** (`dataset_processor.py:460-468`) derive solely from `m_stats` flags (`best_match_is_best`, `is_trash`) and the predicted `filtered` rows — no `taxids_to_use` involved.

`taxids_to_use` is used **only as the composition class schema**, not as the TP reference:
- `get_spurious_composition()` / `get_cross_hit_composition()` (`om_models.py:3055,3074`) receive `tax_df=taxids_to_use`, but use it purely to supply the taxonomic columns in `get_subset_composition()`. The rows are selected by `is_trash` / `is_crosshit` flags and matched to the dataset's own input taxids (`cross_hit_match == taxid`), **not** restricted to `taxids_to_use`.

> Since `taxids_to_use` is the frequency-filtered subset (rows with tax-level frequency `> taxa_threshold`), the TP/recall/precision estimators intentionally bypass it to avoid under-counting. Only the composition breakdown reuses it as the class schema.

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

### `evaluation_results_agent.json`

Agent-parseable aggregate JSON written by `evaluate.py` via
`BatchEvaluationResult.write_agent_output()` (`result_models.py:253-336`). Same schema as
`evaluation_results.json` plus an `aggregate_summary` block — per-column `mean`/`median`/`std`/
`q1`/`q3`/`min`/`max` for `precision_cols`, `recall_cols` (including the baseline-gap
decomposition), and `cross_hit_cols`, plus the `prob_find_true` / `prob_find_any` /
`prob_find_true_clade_clean` products.

### `recall_summary_statistics.tsv`

`pd.DataFrame.describe().T` over the 8 `recall_*` columns + the baseline-gap decomposition
(see `test_datasets_summary_results.tsv`), written by `evaluate.py:423` →
`BatchEvaluator.save_summary_statistics()` (`batch_evaluator.py`). This is the aggregate
counterpart to `precision_summary_statistics.tsv`.

### `recall_complete_summary_statistics.tsv`

Same as `recall_summary_statistics.tsv` but restricted to rows with
`assembly_complete == True`, mirroring the gate applied to recall-facing plots. Written only
when at least one dataset fails the gate (when zero datasets pass, recall plots also fall back
to the full set with a warning).

---

## PNG Images

All generated by `evaluate.py:307` → `ResultVisualizer.plot_all()` (`visualization.py:32-55`), plus training analysis plots from `analyze_cross_hit_distribution` and `analyze_spurious_hit_distribution`. Dimensions ~10×6 to 12×8 inches, matplotlib tight_layout, default 100dpi (unless noted).

### Precision Distribution

**`precision_clade_post_histogram.png`** — Histogram with KDE overlay of `precision_clade_post` across all test datasets. xlim [0, 2], 20 bins.

### Precision Comparison

| Image | Description |
|-------|-------------|
| `precision_metrics_boxplot.png` | Side-by-side boxplot of 6 precision metrics: `precision_best_match`, `purity`, `purity_cov_filtered`, `precision_clade_full`, `precision_clade_post_cleanup`, `precision_clade_fixed`. ylim [0, 3]. |
| `precision_metrics_histograms.png` | Overlaid histograms (hue=Metric) of all precision metrics, 30 bins, xlim [0, 3]. |

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

> The clade products inherit the [pre/post semantics](#test_datasets_summary_resultstsv):
> `*_clade_full` uses the full-tree prediction, `*_clade_clean` uses the post-cleanup tree
> prediction (recall-filtered and, when `--enable-cross-hit`, cross-hit-cleaned).

### Composition Heatmaps

| Image | Description |
|-------|-------------|
| `cross_hit_composition_heatmap.png` | Heatmap of mean cross-hit composition fraction per `tax_level` row. Rows sum to 1. viridis colormap. `--enable-cross-hit` only |
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

Only generated when cross-hit modelling is enabled (`--enable-cross-hit`) and the
`cross_hit_precision` column exists in summary results:

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

Generated by `analyze_cross_hit_distribution()`, only with `--enable-cross-hit`:

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
│   ├── test_datasets_input_df.tsv         per-reference tracking + loss funnel
│   └── test_datasets_reference_leaves.tsv per-leaf tracking
├── evaluation_results.json       aggregated results
├── *.png                         all plots
└── evaluation_report.html        HTML report
```

---

## Regenerating Outputs

Almost everything documented above is produced by the main evaluation run
(`evaluate.py::main()`: cross-hit/spurious training analysis at
`evaluate.py:285-286`, summary statistics at `evaluate.py:411`, results TSV/JSON
at `evaluate.py:302-303`, plots at `evaluate.py:307`, HTML report at
`evaluate.py:311`). Rerun with the pipeline defaults:

```bash
export PYTHONPATH=$(pwd)
source .venv/bin/activate

python deployment/model_evaluation/evaluate.py \
    --study_output_filepath /path/to/study_output \
    --taxid_plan_filepath /path/to/taxid_plan.tsv \
    --analysis_output_filepath /path/to/output
```

- Post-cleanup clade columns are recomputed from the post-cleanup tree by default
  (recall-filtered; see the [pre/post caveat](#test_datasets_summary_resultstsv)).
  Add `--enable-cross-hit` to also apply cross-hit cleanup, or
  `--no-apply-recall-filter` to isolate its effect on the full tree.
- Cross-hit artefacts (model bundle, `cross_hit_*` metrics/composition/plots,
  `cross_hit_summary_statistics.tsv`) are produced **only** with `--enable-cross-hit`.
  Default runs skip them and record `cross_hit_enabled=false` in
  `pipeline_metadata.tsv` / `evaluation_results.json` metadata and in `evaluate.log`.
- Subsequent runs reuse `models/cache/`; pass `--no-cache` to recompute training
  data (also re-emits `training_stats_matrices.tsv`).
- `_predict_clades_postcleanup`, cross-hit cleanup, and the recall filter all run
  per dataset inside `dataset_processor.py::process()`.

Two secondary, documented producers are rerun independently against the same
`--study-output`:

| Producer | Command |
|---|---|
| EDA / cohort metadata | `python deployment/model_evaluation/analysis_data_extractor.py --study-output <dir> --ncbi-db <db> --output-dir <dir> [--explanatory]` — writes `per_dataset_metrics.tsv`, `recall_data.tsv`, and a cohort-scoped `pipeline_metadata.tsv` |
| Recall-gap diagnostic | `python -m deployment.model_evaluation.analysis_scripts.diagnose_recall_gap --study-output <dir> --output-dir <dir> [--limit N] [--no-plot]` — writes `recall_gap_*.tsv` and `recall_gap_buckets.png` |

The experimental scripts under `analysis_scripts/` (`compare_sort_strategies.py`,
`composition_model_comparison.py`, `last_tp_division_prediction_second.py`,
`input_composition_prediction.py`, `debug_metrics_test2.py`) implement their own
pipelines and are **not** documented here — see
[`analysis_scripts/README.md`](analysis_scripts/README.md) and
[`analysis_scripts/documentation/`](analysis_scripts/documentation/).
