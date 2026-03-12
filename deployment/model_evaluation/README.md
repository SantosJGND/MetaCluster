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
| `models.py` | Model training (RecallModeller, CompositionModeller, CrossHitModeller) |
| `batch_evaluator.py` | Batch evaluation across multiple datasets |
| `dataset_processor.py` | Individual dataset processing |
| `data_loader.py` | Input data loading utilities |
| `visualization.py` | Plot generation |
| `metrics.py` | Metric calculations |
| `result_models.py` | Data structures for results |
| `config.py` | Configuration classes |

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
| `--threshold` | No | 0.3 | Threshold for cross-hit filtering |
| `--taxa_threshold` | No | 0.02 | Minimum taxa proportion |
| `--tax_level_to_use` | No | `order` | Taxonomic level |
| `--data_set_divide` | No | 5 | Dataset division for training |
| `--holdout_proportion` | No | 0.3 | Test set proportion |

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

## Output

### Files Generated

| File | Description |
|------|-------------|
| `test_datasets_overall_precision.tsv` | Per-dataset precision scores |
| `test_datasets_summary_results.tsv` | Detailed metrics per dataset |
| `test_datasets_trash_composition.tsv` | Trash (unclassified) composition |
| `test_datasets_cross_hit_composition.tsv` | Cross-hit composition |
| `precision_summary_statistics.tsv` | Summary statistics |
| `models/` | Trained model files |

### Visualizations (PNG)

| Plot | Description |
|------|-------------|
| `overall_precision_histogram.png` | Distribution of precision scores |
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
