# Debug Metrics Test

Standalone debugging script for comparing old (buggy) vs new (fixed) metrics (recall, precision, purity) computed from simulation outputs. Tests m-stats matrix computation, performs set intersection analysis (TP/FP/FN), evaluates filter effects from recall prediction models, compares clade precision between a "selected" (composition-model) method and a "fixed" (distance-threshold) method, and runs a batch clade precision diagnosis across many datasets.

**Important:** This script uses **hardcoded paths** — edit the constants at the top of the file (`STUDY_OUTPUT`, `TAXID_PLAN`, `ANALYSIS_OUT`) before running. There is no CLI argument parser.

## Hardcoded Configuration

```python
STUDY_OUTPUT = Path("/home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/virus/output_study")
TAXID_PLAN = Path("/home/bioinf/Desktop/INSA/Manuscripts/Clustering_Clinical_Metagenomics/data/Panel/viral_assess.tsv")
ANALYSIS_OUT = Path("/home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/virus/filter_xgb_family/")
```

## Methods

### Data Loading

Uses `DataLoader` / `DatasetProcessor` to load study data, NCBI taxonomy, and input taxids. Models are loaded from `ANALYSIS_OUT / "models"` directory:
- Recall models: `recall_gp_clf_pipeline.pkl`, `recall_xgb_bundle.pkl`, `direct_xgb_bundle.pkl`
- Composition models: `composition_xgb_bundle.pkl`, `composition_rf_bundle.pkl`, `composition_gb_bundle.pkl`, `composition_lr_bundle.pkl`, `composition_optuna_bundle.pkl`
- Cross-hit model: loaded via `CrossHitModeller`

### Comparison Methods

#### Old (Buggy) vs New (Fixed) Metrics

Custom functions compute purity, precision, and recall using incorrect denominator logic (counting rows instead of set intersections). These are compared against the correct implementations from `deployment.model_evaluation.metrics` (`compute_purity`, `compute_mstats_precision`, `compute_recall`).

#### Set Intersection Analysis

Direct `set()` intersection of input vs output taxids to compute true TP, FP, FN.

#### Filter Effect Evaluation

Uses `cut_off_recall_prediction` to apply a trained recall model (XGBoost, GP CLF, or Direct XGBoost) as a filter on the overlap manager, then recomputes m-stats to measure the effect on metrics.

#### Clade Prediction Methods

| Method | Description |
|---|---|
| `predict_data_set_clades_composition` | Uses a trained composition model (XGBoost, RF, GB, LR, or OptunaXGB) |
| `predict_data_set_clades_fixed` | Uses a simple distance threshold (`min_dist_threshold=0.6`) |

Clade precision and recall are computed via `calculate_clade_precision` / `compute_clade_recall`.

#### Batch Diagnosis

Iterates over datasets (up to 200), computing clade precision for three pipelines:
- **full**: Pre-crosshit-cleanup, composition model
- **post**: Post-crosshit-cleanup, composition model
- **fixed**: Post-crosshit-cleanup, distance threshold 0.6

Failure modes are recorded:
- **A**: Empty predicted set
- **B**: No intersection (precision=0)
- **C**: Exception

## Outputs

### Tables

#### `clade_precision_diagnostics.csv`

One row per dataset in the batch diagnosis loop. Records clade precision for each pipeline and any failure mode.

| Column | Description |
|---|---|
| `data_set` | Dataset name |
| `precision_full` | Clade precision pre-crosshit-cleanup, composition model |
| `precision_post` | Clade precision post-crosshit-cleanup, composition model |
| `precision_fixed` | Clade precision post-crosshit-cleanup, distance threshold 0.6 |
| `failure_full` | Failure mode for full pipeline: `A` (empty predicted set), `B` (no intersection ⇒ precision=0), `C` (exception), or blank (success) |
| `failure_post` | Same failure mode for post-crosshit-cleanup composition pipeline |
| `failure_fixed` | Same failure mode for post-crosshit-cleanup fixed-threshold pipeline |

### Plots

#### `clade_precision_comparison.png`

Bar chart comparing clade precision between "selected" (composition model, post-crosshit-cleanup) and "fixed" (distance threshold 0.6, post-crosshit-cleanup). X-axis: dataset index or name. Y-axis: clade precision [0, 1]. Two bars per dataset: selected (blue) and fixed (orange). Shows which method achieves higher precision on each dataset.

#### `old_vs_new_metrics.png`

Side-by-side bar chart comparing old (buggy) vs new (fixed) metrics. Three metric groups on the X-axis: purity, precision, recall. Two bars per group: "old" (buggy, counting rows instead of set intersections) and "new" (corrected, intersection-based). Y-axis: metric value [0, 1]. Demonstrates the magnitude of the bug in each metric.

#### `batch_clade_precision_diagnosis.png`

2 × 3 grid of diagnostic plots:
- **Top row**: Histograms of clade precision for each pipeline (full, post, fixed). X-axis: precision [0, 1]. Y-axis: number of datasets. Shows the precision distribution and prevalence of low-precision failures.
- **Bottom row**: Bar charts of failure mode counts per pipeline. X-axis: failure mode (A = empty predicted set, B = no intersection, C = exception). Y-axis: count of datasets. Helps diagnose which failure type dominates each pipeline.
