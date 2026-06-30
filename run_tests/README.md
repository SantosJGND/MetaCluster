# Test Run Configuration

This directory contains configuration files and test scripts for running the metagenomics evaluation pipeline.

## Quick Start

1. Copy `.env.example` to `.env` and set required environment variables:
   ```bash
   cp .env.example .env
   # Edit .env with your settings (e.g., NCBI_EMAIL)
   ```

2. Validate parameters before running tests:
   ```bash
   ./run_tests/check_params.sh simulation
   ```

3. Run tests:
   ```bash
   ./run_tests/test_simulation.sh
   ./run_tests/test_classify.sh
   ./run_tests/test_full_pipeline.sh
   ```

## Params Files

### `params_test.json`

Test configuration with relative paths and placeholder database paths. Used for automated testing.

**Key features:**
- Uses relative paths (e.g., `run_tests/assembly_store` instead of absolute paths)
- Uses conda environment names instead of full paths
- Database indices use placeholder paths (`path/to/...`)

## Input Tables

Sample input tables are provided in `tables/minimal/`:
- `minimal_single.tsv` - Single sample dataset
- `minimal_multi.tsv` - Multiple sample dataset

**Note:** These tables reference `run_tests/assembly_store/` which is not included in the repository. You will need to either:
1. Generate test data using the simulation workflow
2. Provide your own assembly store directory
3. Update the paths in the table to point to your data

## Required Parameters

### Simulation Workflow

| Parameter | Description |
|-----------|-------------|
| `input_table_dir` | Path to input TSV table directory |
| `output_dir` | Output directory |
| `assembly_store` | Path to assembly store |
| `python_bin` | Python binary path |
| `references_extract_script` | Path to reference extraction script |
| `wgsim_python_path` | Path to wgsim simulation script |
| `wgsim_args` | WGSim arguments |
| `mess_conda_env` | Conda environment for mess |
| `minimap2_conda_env` | Conda environment for minimap2 |

### Classification Workflow

| Parameter | Description |
|-----------|-------------|
| `reads` | Path to input reads directory |
| `output_dir` | Output directory |
| `python_bin` | Python binary path |
| `classifier_process_script` | Path to classifier processing script |
| `references_extract_script` | Path to reference extraction script |
| `prinseq_params` | Prinseq QC parameters |
| `prinseq_conda_env` | Conda environment for prinseq |
| `kraken2_index` | Kraken2 database index path |
| `kraken2_conda_env` | Conda environment for kraken2 |
| `centrifuge_index` | Centrifuge database index path |
| `centrifuge_conda_env` | Conda environment for centrifuge |
| `diamond_index` | Diamond database index path |
| `diamond_conda_env` | Conda environment for diamond |
| `krakenunique_index` | KrakenUnique database index path |
| `krakenunique_conda_env` | Conda environment for krakenunique |
| `minimap2_illumina_params` | Minimap2 parameters for Illumina |
| `minimap2_conda_env` | Conda environment for minimap2 |
| `msamtools_params` | Msamtools filter parameters |
| `msamtools_conda_env` | Conda environment for msamtools |
| `mosdepth_conda_env` | Conda environment for mosdepth |

## Notes

- All test scripts accept `PARAMS_FILE` environment variable to specify a custom params file
- The `-params-file` argument is used to pass parameters to Nextflow workflows
- Parameter validation runs automatically before each test script executes
