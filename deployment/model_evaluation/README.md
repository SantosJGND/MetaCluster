# Model Evaluation

Evaluates clustering and classification results using trained models to assess precision, recall, and composition metrics.

## Overview

The model evaluation pipeline processes study output directories containing:
- Input simulation tables
- Clustering results
- Classification outputs

It then generates comprehensive metrics and visualizations comparing predicted vs. actual results.

## Modules

### Core Modules

| Module | Description |
|--------|-------------|
| `evaluate.py` | Main evaluation script |
| `models.py` | Model training (RecallModeller, BaseCompositionModeller + 5 variants, CrossHitModeller) |
| `batch_evaluator.py` | Batch evaluation across multiple datasets |
| `dataset_processor.py` | Individual dataset processing |
| `data_loader.py` | Input data loading utilities |
| `visualization.py` | Plot generation |
| `metrics.py` | Metric calculations |
| `result_models.py` | Data structures for results |
| `reference_details.py` | Per-reference (`loss_stage`) and per-leaf tracking tables |
| `config.py` | Configuration classes |
| `analysis_scripts/` | Experimental analysis scripts (sort-strategy comparison, composition model comparison, last-TP-division prediction) |

## Available Models

Three model types are trained per pipeline run: recall, composition, and cross-hit. The recall model has multiple variants selectable via `--recall_model_interface`.

### Recall Models

The recall model predicts the fraction of reads/leaves to keep to achieve a target recall. All variants share the same API (`train_model`, `predict_cutoff`, `save_model`, `load_model`).

| `--recall_model_interface` | Class | Description |
|---|---|---|
| `xgb` (default) | `RecallModeller` | Multi-output XGBoost regressor (100 trees) predicting the full recall curve, then finding the percentile for the target recall |
| `direct_xgb` | `DirectXGBRecallModeller` | Direct-fraction XGBoost regressor (200 trees) with asymmetric sample weights (underestimation penalised 3:1). Predicts τ-crossing fraction directly without bin-index search. |
| `morf` | `RecallModeller` | Multi-output Random Forest regressor (100 trees, no hyperparameter tuning) |
| `moxgb_optimized` | `RecallModeller` | XGBoost with Optuna optimisation (50 trials, 10-fold CV) |
| `morf_optimized` | `RecallModeller` | Random Forest with Optuna optimisation (50 trials, 10-fold CV) |
| `monn_optimized` | `RecallModeller` | MLP neural network with Optuna optimisation (50 trials, 10-fold CV) |
| `gp_clf` | `GPCLFRecallModeller` | Per-division Gaussian Process regressors with automatic grid-search for optimal τ (recall target) and X (confidence threshold). Enables diagnostic plots (`recall_landscape.png`, `recall_calibration.png`, `recall_actual_at_index.png`). |
| `direct` | `CutoffRecallModeller` | Random forest classifier predicting minimum bin count `k_min` to reach the target recall directly. Supports probability-guided cutoff via `--cutoff_confidence`. Legacy; superseded by `direct_xgb`. |

### Composition Models

Five composition model variants are selectable via `--composition_model_interface`. All share the same `BaseCompositionModeller` ABC with `fit`, `predict_proba`, `save_model`, `load_model`, and `eval_and_plot`.

| `--composition_model_interface` | Class | Description |
|---|---|---|
| `xgb` (default) | `XGBCompositionModeller` | XGBoost classifier (300 trees, max_depth=6) via sklearn Pipeline + ColumnTransformer |
| `xgb_optimized` | `OptunaXGBCompositionModeller` | XGBoost + Optuna hyperparameter search (50 trials, 10-fold CV) via existing `ClusteringPipeline` |
| `rf` | `RFCompositionModeller` | Random Forest classifier (300 trees, max_depth=12, balanced class_weight) |
| `gb` | `GBCompositionModeller` | Gradient Boosting classifier (300 trees, max_depth=5, lr=0.1, subsample=0.8) |
| `lr` | `LRCompositionModeller` | Logistic Regression (C=1.0, balanced, stats-only features via `remainder='drop'`) |

### Cross-Hit Model

Predicts cross-hit probability and is used to clean predicted cross-hit leaves
from the post-cleanup tree. The model is **trained and applied only when
`--enable-cross-hit`** (default: disabled). With the default `--no-enable-cross-hit`
the model is not trained and all cross-hit-specific outputs (metrics TSVs, plots,
composition, summary statistics) are skipped; this is recorded as a
`cross_hit_enabled` row in `pipeline_metadata.tsv` / `evaluation_results.json`
metadata and in `evaluate.log`.

## Usage

### Command Line

```bash
# Activate environment
source .venv/bin/activate
export PYTHONPATH=$(pwd)

# Run evaluation
python deployment/model_evaluation/evaluate.py \
    --study_output_filepath /path/to/study_output \
    --taxid_plan_filepath /path/to/taxid_plan.tsv \
    --analysis_output_filepath /path/to/output
```

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--study_output_filepath` | Yes | - | Path to study output directory |
| `--taxid_plan_filepath` | Yes | - | Path to taxid plan TSV |
| `--analysis_output_filepath` | Yes | - | Path for analysis output |
| `--recall_model_interface` | No | `xgb` | Recall model variant (see Available Models) |
| `--composition_model_interface` | No | `xgb` | Composition model variant: `xgb`, `xgb_optimized`, `rf`, `gb`, `lr` (see Composition Models) |
| `--target_recall` | No | 1.0 | Target recall threshold for cutoff prediction |
| `--cutoff_confidence` | No | - | Confidence level for prob-guided cutoff (`direct` only) |
| `--threshold` | No | 0.3 | Threshold for cross-hit filtering |
| `--taxa_threshold` | No | 0.02 | Minimum taxa proportion |
| `--tax_level_to_use` | No | `order` | Taxonomic level |
| `--data_set_divide` | No | 16 | Dataset division for training |
| `--holdout_proportion` | No | 0.3 | Test set proportion |
| `--enable-cross-hit` | No | false | Enable cross-hit modelling and cleanup during evaluation (default: disabled). When disabled the cross-hit model is not trained and cross-hit-specific outputs are skipped |
| `--apply-recall-filter` | No | true | Apply recall-model leaf truncation before post-cleanup clade prediction. Use `--no-apply-recall-filter` to isolate the cross-hit cleanup effect on the full tree |

### Caching

The pipeline caches parsed training data to avoid re-processing datasets on subsequent runs.

```bash
# First run - computes and caches training data
python deployment/model_evaluation/evaluate.py \
    --study_output_filepath /path/to/study_output \
    --taxid_plan_filepath /path/to/taxid_plan.tsv \
    --analysis_output_filepath /path/to/output

# Subsequent runs - uses cached data (fast)
python deployment/model_evaluation/evaluate.py \
    --study_output_filepath /path/to/study_output \
    --taxid_plan_filepath /path/to/taxid_plan.tsv \
    --analysis_output_filepath /path/to/output

# Force recompute (ignore cache)
python deployment/model_evaluation/evaluate.py \
    --study_output_filepath /path/to/study_output \
    --taxid_plan_filepath /path/to/taxid_plan.tsv \
    --analysis_output_filepath /path/to/output \
    --no-cache
```

Cache is stored in: `{analysis_output_filepath}/models/cache/`

| Argument | Default | Description |
|----------|---------|-------------|
| `--no-cache` | false | Force recompute cached training data |

### Python API

```python
from deployment.model_evaluation import BatchEvaluator, EvaluatorConfig
from deployment.model_evaluation.visualization import ResultVisualizer
from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper

# Configure
config = EvaluatorConfig(
    study_output_filepath="/path/to/output",
    tax_level="order"
)

# Run evaluation
evaluator = BatchEvaluator(config, models, ncbi_wrapper, input_tax_df, taxids_to_use)
results = evaluator.evaluate(test_datasets)

# Generate visualizations
visualizer = ResultVisualizer("/path/to/output/plots")
visualizer.plot_all(results)
```

## Assembly-Completeness Gate

Recall/denominator-facing analyses only use datasets where **at most
`max_missing_pct` (default 5%) of detected references
(`uniq_reads >= min_uniq_reads`) lack a matched assembly**.
A dataset passes when at most `--max-missing-refs-pct` percent of the taxids in
`classification/<ds>_merged_classification.tsv` with `uniq_reads >= --min-uniq-reads`
(default 1) are missing from `output/matched_assemblies.tsv` — i.e. the fraction in
`unmatched_classified_taxids.tsv` is `<= max_missing_pct`.
This mirrors the `retrieve` gate (`--max_missing_pct` in `reference_management/main.py`):
a dataset where >5% of references could not be mapped to assemblies is **skipped**
(excluded), so 0-recall outputs caused by a taxonomy→genome lookup failure are never
treated as genuine, while a small tolerated fraction does not discard the dataset.

- **On by default** (`--require-complete-assemblies`); add `--allow-incomplete`
  to include failing datasets. Applies to: `analysis_data_extractor.py`
  (`recall_per_classifier.tsv`, `recall_data.tsv` for the GEE recall model,
  `aggregate_statistics.tsv`), `evaluate.py` (recall model training matrices +
  recall-facing plots), `diagnose_recall_gap.py`, and `compare_sort_strategies.py`
  (recall-model training). `--min-uniq-reads` sets the detection threshold;
  `--max-missing-refs-pct` (default 5) sets the tolerated unmatched fraction.
- **At the classify stage** (`deployment/classify/classify.nf` +
  `deploy_classify_map.sh`): `params.max_missing_references_pct` (default 5) is
  forwarded as `retrieve --max_missing_pct`. Above-threshold datasets exit with
  code 3, are recorded in `incomplete_datasets.tsv` with reason
  `too_many_missing_references`, and are excluded from analysis without failing
  the deployment (`exit 0`); datasets at or below the threshold proceed with the
  matched subset.
- **Not filtered** (precision/composition): TP/cross-hit/spurious counts,
  precision and purity, classifier hit counts, cross-hit/spurious composition.
- Every `per_dataset_metrics.tsv` / `test_datasets_summary_results.tsv` row carries
  an `assembly_complete` flag; `pipeline_metadata.tsv` records
  `recall_analysis_datasets` and `assembly_incomplete_excluded`.

As of the virus study (1690 datasets, `--min-uniq-reads 1`, strict 0%): 128 (7.8%)
datasets pass study-wide; 35 of the 431 final-analysis test datasets present
locally pass. All 18 viral orders remain represented. A 5% tolerance raises these
counts. Caveat: the pre-fix runs mostly fail the gate because of the environmental
store/NCBI mismatch (see CHANGELOG), not bad data — re-running classify with the
corrected store flips most datasets to pass.

## Output

### Files Generated

| File | Description |
|------|-------------|
| `test_datasets_input_df.tsv` | Per-reference tracking table; one row per `(data_set, input taxid)`. Historical 12-column schema unchanged, followed by classifier/mapping/clade evidence and a `loss_stage` funnel |
| `test_datasets_reference_leaves.tsv` | Per-leaf tracking table; one row per `(data_set, input taxid, matched assembly)`. Explains each `loss_stage` |
| `test_datasets_overall_precision.tsv` | Per-dataset precision scores |
| `test_datasets_summary_results.tsv` | Detailed metrics per dataset |
| `pipeline_metadata.tsv` | Run summary (dataset counts, skipped/failed names) |
| `test_datasets_spurious_composition.tsv` | Spurious (unclassified) composition |
| `test_datasets_cross_hit_composition.tsv` | Cross-hit composition (`--enable-cross-hit` only) |
| `precision_summary_statistics.tsv` | Summary statistics for precision metrics |
| `recall_summary_statistics.tsv` | Summary statistics for recall metrics (incl. baseline gap decomposition) |
| `recall_complete_summary_statistics.tsv` | Same, restricted to datasets passing the assembly-completeness gate (only when some datasets fail it) |
| `cross_hit_summary_statistics.tsv` | Summary statistics for cross-hit metrics (`--enable-cross-hit` only) |
| `evaluation_results_agent.json` | Agent-parseable JSON with per-column aggregate statistics (mean/median/std/quartiles) |
| `models/` | Trained model files (`cross_hit_xgb_bundle.pkl` only with `--enable-cross-hit`) |
| `models/` | Trained model files |
| `models/cache/` | Cached training data (parquet files) |

> The EDA producer (`analysis_data_extractor.py`) emits the same cohort metadata
> schema (`pipeline_metadata.tsv`: `total_attempted`, `extracted`, `dropped`,
> `failed`, `skipped`, plus name/message lists and `study_gaps`) scoped to the full
> study cohort in its output directory. `evaluate.py`'s TSV + `evaluation_results.json`
> report identical `successful`/`failed`/`skipped_count`, plus a
> `cross_hit_enabled` flag documenting whether cross-hit modelling/filtering ran
> (referenced in `evaluate.log`).

### Visualizations (PNG)

Cross-hit-specific plots (`cross_hit_composition_heatmap.png`, cross-hit metrics
boxplot/distribution/improvement, and the training cross-hit distribution plots)
are generated **only with `--enable-cross-hit`**.

| Plot | Description |
|------|-------------|
| `precision_clade_post_histogram.png` | Distribution of precision scores |
| `precision_metrics_boxplot.png` | Comparison of precision metrics |
| `precision_metrics_histogram.png` | Precision metrics distribution |
| `recall_metrics_boxplot.png` | Comparison of recall metrics |
| `recall_improvement_histogram.png` | Recall improvement distribution |
| `probability_metrics_boxplot.png` | Probability metrics comparison |
| `cross_hit_composition_heatmap.png` | Cross-hit composition heatmap |
| `trash_composition_heatmap.png` | Trash composition heatmap |
| `clade_precision.png` | Clade precision by taxonomic level |

### HTML Report

Generate HTML report embedding all plots:

```python
from deployment.model_evaluation.visualization import ResultVisualizer

visualizer = ResultVisualizer(output_dir)
visualizer.plot_all(results)
html_path = visualizer.generate_html_report(results)
```

## Metrics

### Precision Metrics

| Metric | Description |
|--------|-------------|
| `overall_precision_raw` | Unique correct predictions / total predictions |
| `fuzzy_precision_raw` | Predictions with >0 coverage / total |
| `fuzzy_precision_cov_filtered` | After coverage filtering |
| `clade_precision_full` | Composition-model clade prediction on the full tree (before recall filter or cross-hit cleanup) |
| `clade_precision_post` | Composition-model clade prediction on the post-cleanup tree: recall-filtered (top `keep_index` leaves) and, when `--enable-cross-hit`, additionally cross-hit-cleaned. With `--no-apply-recall-filter`, the full tree (optionally cross-hit-cleaned) |

### Recall Metrics

| Metric | Description |
|--------|-------------|
| `recall_raw` | Correct predictions / total expected |
| `recall_cov_filtered` | After coverage filter |
| `clade_recall` | Composition-model clade prediction on the full tree (before recall filter or cross-hit cleanup) |
| `recall_filtered_leaves` | After leaf filtering |

#### Baseline recall decomposition

The gap between `recall_baseline` and `recall_baseline_cov_filtered` is decomposed
into two attributable columns on every `test_datasets_summary_results.tsv` row:

| Column | Definition | Meaning |
|--------|-----------|---------|
| `recall_classification_credit` | `recall_baseline − recall_baseline_assembly` | Input taxids credited only by classifier evidence (`uniq_reads > 0`) with no best-matched assembly |
| `recall_zero_coverage_loss` | `recall_baseline_assembly − recall_baseline_cov_filtered` | Best-matched assemblies whose leaf has `coverage == 0` (mapping-stage loss) |
| `recall_baseline_cov_gap` | `recall_baseline − recall_baseline_cov_filtered` | The two components summed |

> Best-match marking (`_compute_best_matches`, `metagenomics_utils/overlap_manager/node_stats.py`)
> prefers a **covered** assembly within a `best_match_taxid` group: coverage>0 rows sort ahead of
> zero-coverage rows (score still decides within each coverage class), so a taxid whose top-scoring
> assembly has `coverage == 0` is still counted by `recall_baseline_cov_filtered` when a covered
> assembly exists. Single-candidate groups (the legacy F4 case) are unchanged.
> `recall_baseline` itself is unaffected — every group still emits one best row.

### Per-Reference Loss Attribution

Aggregate recall says *how much* was recovered but not *which references were lost or
why*. `reference_details.py` adds that attribution to two tables (full column reference
in `output_description.md`):

- `test_datasets_input_df.tsv` — the existing per-`(data_set, taxid)` table, widened in
  place. The historical 12 columns keep their names and order; classifier evidence,
  best-match leaf statistics, clade membership and a `loss_stage` bucket are appended.
- `test_datasets_reference_leaves.tsv` — one row per `(data_set, taxid, matched
  assembly)`, so each `loss_stage` can be traced to the leaves that caused it.

`loss_stage` is an ordered **first-match** cascade in pipeline order
(`absent_classification` → `survived_classifier_only` → `classified_no_assembly` →
`assembly_no_coverage` → `assembly_zero_reads` → `lost_recall_truncation` →
`lost_clade_prediction` → `survived_post_cleanup`). Each reference lands in exactly one
bucket, so the counts form a proper funnel. The names `absent_classification`,
`classified_no_assembly`, `assembly_no_coverage` and `assembly_zero_reads` are kept
identical to `analysis_scripts/diagnose_recall_gap.py` so the two diagnostics can be
compared directly.

**`loss_stage` uses the same definition as the headline recall.** A reference is recalled
(classification-aware, matching `recall_baseline`) when it has a clean assembly match
(`n_clean_matches > 0`) **or** the classifier called it. Only two buckets are therefore
actual recall losses — `absent_classification` and `classified_no_assembly` — which gives
the invariant:

```python
recalled_classification_aware == loss_stage not in ("absent_classification", "classified_no_assembly")
```

Every other bucket describes how far an already-recalled reference got. In particular a
reference the classifier found but that has no clean assembly match is
`survived_classifier_only`, **not** `classified_no_assembly`: labelling it a loss would
contradict the reported recall. The `recalled_classification_aware`,
`recalled_assembly` and `recalled_classifier_only` flags make the decomposition explicit
so the invariant can be asserted directly on the table.

Two conditions were dropped from the cascade because they can never win it, and are
reported instead through the `lost_no_best_match` / `lost_best_match_trash` flags and the
`n_clean_matches` / `is_trash_best` columns:

- `no_best_match` — a reference with leaves always has a best match.
- `best_match_trash` — a clean match means the best leaf is not trash, and a trash-only
  reference is recalled by the classifier.

The boolean flags alongside `loss_stage` are deliberately **independent** rather than a
copy of it: a reference that fails several conditions reports every applicable flag, so
the flags can be used to verify the cascade.

| Flag group | Flags |
|---|---|
| Recall decomposition | `recalled_classification_aware`, `recalled_assembly`, `recalled_classifier_only` |
| Upstream evidence | `lost_before_classification`, `lost_assembly_lookup` (no *clean* assembly match, independent of recall), `lost_no_coverage`, `lost_zero_reads` |
| Diagnostics (unreachable as stages) | `lost_no_best_match`, `lost_best_match_trash` |
| Downstream | `lost_recall_truncation`, `lost_clade_prediction`, `survived_post_cleanup` |

Note `lost_assembly_lookup` is independent of recall: a classifier-credited reference whose
only leaves are trash is `lost_assembly_lookup` yet still recalled.

Key design points:

- Everything keys on `m_stats["best_match_taxid"]`, which is what `metrics.compute_recall`
  and `compute_clade_recall` use. `output/matched_assemblies.tsv` is **not** the join key —
  its `taxid` column holds the matched *assembly's* taxid, so joining input taxids on it
  reports that nothing was ever matched.
- `n_clean_matches` reproduces the recall criterion exactly (≥1 leaf with
  `is_trash == False and best_match_is_best == True`), so `n_clean_matches > 0` yields
  the same taxid set as `compute_recall`.
- `classifier_uniq_reads` is empty rather than `0` when the classification table carries no
  read-count column, so "unavailable" stays distinguishable from "zero reads".
  `classifier_reads_available` makes the same distinction explicit, because
  `analysis_data_extractor.load_classifier_detected_taxids()` drops its `uniq_reads > 0`
  filter entirely when the column is absent — on such (legacy) studies *every* taxid in
  the classification table counts as detected. `loss_stage` reproduces that behaviour
  rather than assuming a zero count.
- Building these tables is best-effort: a failure logs a warning and emits the row without
  the new columns rather than aborting a run whose metrics are already computed.

> Known issue: `--max_taxids_fixed_filter` / `EvaluatorConfig.max_taxids_fixed_filter`
> are accepted but never read — `_apply_fixed_filter` uses a hard-coded `12`. The limit
> also truncates leaf rows, not unique taxids. See `output_description.md`
> § "Known issue — `max_taxids_fixed_filter` is not wired up".

### Pre- vs Post-Cleanup Semantics

The "pre-cleanup" / "post-cleanup" clade metrics measure the composition model's
clade prediction on two different trees. The actual pipeline stages
(see `dataset_processor.py`) are:

1. **Pre-cleanup** — the composition model predicts clades on the **full tree**
   (all leaves of the original OverlapManager). Populates `clade_precision_full`
   and `clade_recall` / `recall_clade_pre_cleanup`.
2. **Recall filter** (on by default) — `cut_off_recall_prediction` builds a
   **new, truncated** OverlapManager with `max_proportion=target_percentile`,
   keeping only the top `keep_index` leaves sorted by `total_uniq_reads`. The
   tree topology itself changes (fewer leaves, different root nodes). With
   `--no-apply-recall-filter` this step is skipped and the post-cleanup tree is
   the original full OM.
3. **Cross-hit cleanup** (when `--enable-cross-hit`) — `_apply_crosshit_cleanup`
   mutates the post-cleanup tree in place (zeroes `numreads`/`coverage` of
   predicted cross-hit leaves, prunes, rebuilds). Its metrics are scoped to the
   kept leaves of that tree.
4. **Post-cleanup** — the composition model predicts clades on the resulting
   tree. Populates `clade_precision_post` / `recall_clade_post_cleanup`, plus
   the fixed (`min_dist=0.6`) variant `clade_precision_fixed` /
   `recall_fixed_max_12`.

Consequences:

- `clade_precision_post` / `recall_clade_post_cleanup` run on the
  **recall-filtered tree** by default, so they differ from their pre-cleanup
  counterparts **even when `enable_cross_hit=False`** — the difference reflects
  the recall-filter truncation.
- With `--enable-cross-hit`, the cross-hit cleanup now **feeds the post-cleanup
  columns**: the cleaned tree is the one used for `clade_precision_post` /
  `recall_clade_post_cleanup`. Toggling it therefore changes those values.
- Use `--no-apply-recall-filter` together with `--enable-cross-hit` to isolate
  the cross-hit cleanup effect on the full tree. Note that with
  `--no-apply-recall-filter` the recall model is not invoked, so the
  `recall_metric_*` diagnostics and `recall_filtered_leaves` are absent.
- The cross-hit metrics themselves (`predicted_cross_hits`,
  `cross_hit_specificity`, `cross_hit_precision`, `cross_hit_recall`,
  `cross_hit_f1`) describe the cleanup of the post-cleanup tree and are
  therefore scoped to its kept leaves.

### Probability Metrics

| Metric | Description |
|--------|-------------|
| `Prob_Find_any` | recall_raw × fuzzy_precision_raw |
| `Prob_Find_true` | recall_raw × overall_precision_raw |
| `Prob_Find_true_clade_full` | clade_recall × clade_precision_full |

## Data Structures

### Input Table Format

TSV with columns:
- `sample`: Sample identifier
- `taxid`: NCBI taxid
- `reads`: Number of reads
- `mutation_rate`: Mutation rate (0.0-1.0)
- `accid`: Assembly accession

### Study Output Structure

```
study_output/
├── dataset_001/
│   ├── input/
│   │   └── dataset_001.tsv
│   ├── output/
│   │   └── clade_report_with_references.tsv
│   └── clustering/
│       └── (clustering output files)
├── dataset_002/
│   └── ...
└── ...
```

## Changelog

### 2026-10-02

- **Per-reference loss attribution.** `test_datasets_input_df.tsv` is widened in place
  (no new flag): the historical 12 columns keep their names and order, and classifier
  evidence, best-match leaf statistics, clade membership (pre/post/fixed), an ordered
  `loss_stage` funnel and independent `lost_*` flags are appended. New
  `test_datasets_reference_leaves.tsv` gives one row per
  `(data_set, taxid, matched assembly)` so each funnel bucket can be traced to its
  leaves. See § "Per-Reference Loss Attribution". Implementation in the new
  `reference_details.py`.
- **`dataset_processor` internal signature changes.** `_predict_clades_precleanup` now
  returns `(result, clades_df)`, `_predict_clades_postcleanup` returns
  `(result, clades_df, fixed_clades_df)`, and `_apply_recall_filter` returns a third
  element, `keep_index`. These are private methods; `process()` collects the frames
  explicitly rather than through a mutable out-parameter.
- **`input_df["genus"]` is now populated.** `data_loader.expand_input_data` resolves the
  genus when the input table does not carry one, matching how it already handled
  `order`/`family`.
- Documented that `--max_taxids_fixed_filter` / `EvaluatorConfig.max_taxids_fixed_filter`
  are accepted but never read (`_apply_fixed_filter` hard-codes `12`), and that the limit
  truncates leaf rows rather than unique taxids. Behaviour unchanged.

### 2026-09-21

- **Cross-hit analysis is now optional (default off).** Opt in with
  `--enable-cross-hit`. When disabled: the cross-hit model is not trained, the 9
  `cross_hit_*` columns are absent from `test_datasets_summary_results.tsv` (not
  zero-filled), and `cross_hit_xgb_bundle.pkl`, `crosshit_model.joblib`,
  `test_datasets_cross_hit_composition.tsv`, `cross_hit_metrics_*.tsv`,
  `cross_hit_summary_statistics.tsv`, and the cross-hit PNG plots are not
  generated. The skip is recorded as a `cross_hit_enabled` row in
  `pipeline_metadata.tsv` / `evaluation_results.json` metadata and as a startup
  marker in `evaluate.log`.
- **Coherence fix: post-cleanup clade columns now reflect the actual cleanups.**
  `precision_clade_post_cleanup`, `recall_clade_post_cleanup`,
  `precision_clade_fixed`, and `recall_fixed_max_12` are computed on the
  recall-filtered tree by default (optionally cross-hit-cleaned with
  `--enable-cross-hit`) instead of the full overlap matrix. Added
  `recall_after_recall_filter` / `recall_fixed_max_12` summary columns.
  `--no-apply-recall-filter` restores full-tree post-cleanup behavior, in which
  case the recall model is not invoked and the `recall_metric_*` diagnostics and
  `recall_filtered_leaves` are absent.

## Dependencies

- Python 3.10+
- See `.venv` for packages:
  - pandas, numpy
  - scikit-learn
  - xgboost
  - matplotlib, seaborn
  - biopython
