# Mapping Cluster

[![CI](https://github.com/insa/metagenomics-evaluation-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/insa/metagenomics-evaluation-pipeline/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

Metagenomic analysis creates a wealth of data that can be challenging to interpret. One of the main challenges in this field is to discriminate between Cross-Hits — metagenomic classifiers often return multiple related hits, making it difficult to determine which hit is the most relevant for a given sample. This project benchmarks the effectiveness of clustering methods in resolving these Cross-Hits.

## Quick Start

```bash
# Install Python dependencies
pip install -e .

# Run the Nextflow pipeline
nextflow run main.nf -profile conda --params_file deployment/params.json
```

See [run_tests/](run_tests/) for working configuration examples.

## Project Structure

The project is structured to facilitate the evaluation of clustering methods on metagenomic data:

- **`simulation/`** — Simulating reads from reference assemblies to create controlled datasets
- **`classification/`** — Classifying metagenomic reads using Centrifuge and Kraken2 classifiers
- **`reference_management/`** — Managing reference assemblies (fetching from NCBI, local DB management)
- **`clustering/`** — Clustering classified reads and evaluating the clusters
- **`deployment/`** — Nextflow pipeline workflows, model evaluation, and ML API
  - **`workflows/`** — Nextflow workflow definitions (`benchmark.nf`, `classify_only.nf`)
  - **`modules/`** — Reusable Nextflow process modules
  - **`model_evaluation/`** — Evaluation pipeline with metrics, visualization, and trained models
  - **`ml_api/`** — FastAPI-based ML inference server
  - **`ml_api_client/`** — CLI and Python client for the ML API
- **`metagenomics_utils/`** — Shared utilities for NCBI, data processing, and overlap analysis
- **`tests/`** — pytest test suite (137+ tests)

### Shared Utilities

The `metagenomics_utils/` package contains shared utilities:

- **`ncbi_tools.py`** — NCBI tools for taxonomy and sequence data retrieval
- **`reference_utils.py`** — Reference management utilities
- **`dataframe_utils.py`** — DataFrame utilities for tabular data
- **`overlap_manager/`** — Overlap analysis for clustering evaluation
  - `core.py` — Core overlap analysis
  - `manager.py` — OverlapManager class
  - `stats.py` — Statistics functions
  - `diversity.py` — Diversity metrics
  - `node_stats.py` — Node-level statistics

### Deployments

- [**deployment/simulation/**](deployment/simulation) — Simulate reads from reference assemblies
- [**deployment/classify/**](deployment/classify) — Classify metagenomic reads
- [**deployment/model_evaluation/**](deployment/model_evaluation) — Evaluate classification/clustering models
- [**deployment/ml_api/**](deployment/ml_api) — ML model serving API
- [**deployment/ml_api_client/**](deployment/ml_api_client) — ML API client SDK

## Requirements

- Nextflow
- Python 3.10+
- Conda (for classifier and mapping tool environments)

## Installation

### Python Environment

```bash
pip install -e .
```

### Conda Environments

Each pipeline stage runs in its own conda environment. The `-profile conda` Nextflow flag activates the correct environment per process.

| Params key | YAML file | Conda env name | Tools |
|---|---|---|---|
| `prinseq_conda_env` | `envs/preproc_env.yaml` | `preproc` | prinseq-plus-plus, minimap2, seqtk |
| `minimap2_conda_env` | `envs/mapping_env.yaml` | `mapping` | minimap2, samtools |
| `msamtools_conda_env` | `envs/msamtools_env.yaml` | `msamtools` | msamtools |
| `wgsim_conda_env` | `envs/simulate_env.yaml` | `simulate` | wgsim |
| `kraken2_conda_env`, `centrifuge_conda_env`, `diamond_conda_env`, `krakenunique_conda_env` | `envs/classification_env.yaml` | `classification` | kraken2, centrifuge, diamond, krakenuniq |

```bash
conda env create -f envs/preproc_env.yaml
conda env create -f envs/mapping_env.yaml
conda env create -f envs/msamtools_env.yaml
conda env create -f envs/simulate_env.yaml
conda env create -f envs/classification_env.yaml
```

After creating the environments, set each `*_conda_env` parameter in your params JSON file.

### Running the Pipeline

```bash
nextflow run main.nf -profile conda --params_file deployment/params.json
```

## Running Benchmark Analysis

```bash
python deployment/model_evaluation/evaluate.py \
  --study_output_filepath /path/to/nextflow/output \
  --taxid_plan_filepath /path/to/taxid_plan.tsv \
  --analysis_output_filepath /path/to/analysis/output
```

Requirements:
- `ncbi-taxonomist` CLI tool installed and in PATH
- `NCBI_EMAIL` environment variable set (e.g., in `.env` file)

## Reference Management

References used for mapping are stored locally in the `assembly_store` directory (specified in `params.json`). Reference identification uses the `Entrez` tool suite via `Bio.Entrez`.

The software first checks for a specified taxid in the `nucleotide` database, falling back to the `assembly` database.

To check reference availability before running the full pipeline:

```bash
python reference_management/main.py check \
  --input_table path/to/taxid_table.tsv \
  --assessment /path/to/assessment.tsv
```

## Model Evaluation Output

See [deployment/model_evaluation/output_description.md](deployment/model_evaluation/output_description.md) for a complete description of evaluation outputs, including metrics, plots, and model artifacts.

## ML API

A FastAPI-based inference server serves trained models for recall prediction, composition stop-traversal prediction, and clustering threshold prediction.

```bash
docker build -t ml_api deployment/ml_api/
docker run -p 8000:8000 ml_api
```

See [deployment/ml_api/README.md](deployment/ml_api/README.md) for endpoint documentation.

## Tests

```bash
# Run all tests
pytest

# Run unit tests only (fast, no external dependencies)
pytest -m "not slow and not integration and not requires_ncbi"
```

## Citation

If you use this software in your research, please cite:

```bibtex
@software{dourado2026mapping,
  author = {Dourado, Joao},
  title = {{Mapping Cluster: Metagenomics Cross-Hit Resolution Benchmark}},
  year = {2026},
  doi = {10.5281/zenodo.XXXXXXX},
  url = {https://github.com/insa/metagenomics-evaluation-pipeline}
}
```

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
