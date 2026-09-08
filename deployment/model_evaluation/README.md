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
| `test_datasets_overall_precision.tsv` | Per-dataset precision scores |
| `test_datasets_summary_results.tsv` | Detailed metrics per dataset |
| `pipeline_metadata.tsv` | Run summary (dataset counts, skipped/failed names) |
| `test_datasets_spurious_composition.tsv` | Spurious (unclassified) composition |
| `test_datasets_cross_hit_composition.tsv` | Cross-hit composition |
| `precision_summary_statistics.tsv` | Summary statistics for precision metrics |
| `cross_hit_summary_statistics.tsv` | Summary statistics for cross-hit metrics |
| `models/` | Trained model files |
| `models/cache/` | Cached training data (parquet files) |

> The EDA producer (`analysis_data_extractor.py`) emits the same cohort metadata
> schema (`pipeline_metadata.tsv`: `total_attempted`, `extracted`, `dropped`,
> `failed`, `skipped`, plus name/message lists and `study_gaps`) scoped to the full
> study cohort in its output directory. `evaluate.py`'s TSV + `evaluation_results.json`
> report identical `successful`/`failed`/`skipped_count`.

### Visualizations (PNG)

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
| `clade_precision_full` | Precision with predicted clades |
| `clade_precision_post` | After cross-hit cleanup |

### Recall Metrics

| Metric | Description |
|--------|-------------|
| `recall_raw` | Correct predictions / total expected |
| `recall_cov_filtered` | After coverage filter |
| `clade_recall` | With predicted clades |
| `recall_filtered_leaves` | After leaf filtering |

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

## Dependencies

- Python 3.10+
- See `.venv` for packages:
  - pandas, numpy
  - scikit-learn
  - xgboost
  - matplotlib, seaborn
  - biopython
