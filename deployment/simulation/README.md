## Simulation Workflow for Metagenomics Benchmarking

This directory contains a simulation-only workflow for generating synthetic metagenomic reads from reference sequences.

# Overview

The simulation workflow:
1. Extracts reference sequences from taxonomic IDs in the input table
2. Simulates paired-end reads from the references

The output can then be used by:
- `deployment/classify/` - for classification
- `deployment/model_evaluation/` - for analysis

# Running the Simulation

## Basic Usage

```bash
nextflow run simulate.nf -profile conda \
  --input_table /path/to/input_table.tsv \
  --output_dir /path/to/output
```

# Output Structure

```
{output_dir}/
└── {dataset_name}/
    ├── input/
    │   └── {dataset_name}.tsv          # Input table
    └── fastq/
        ├── {dataset_name}_R1.fq.gz    # Simulated reads (forward)
        └── {dataset_name}_R2.fq.gz    # Simulated reads (reverse)
```

# Complete Workflow

To run the full benchmark pipeline:

1. **Simulate reads:**
   ```bash
   nextflow run deployment/simulation/simulate.nf -profile conda \
     --input_table /path/to/input_table.tsv \
     --output_dir /path/to/simulation_output
   ```

2. **Run classification and clustering:**
   ```bash
   nextflow run deployment/classify/classify.nf \
     -profile conda \
     --reads /path/to/simulation_output/{dataset_name}/fastq \
     --output_dir /path/to/classification_output
   ```

3. **Analyze results:**
   ```bash
   python deployment/model_evaluation/evaluate.py \
     --study_output_filepath /path/to/classification_output \
     --taxid_plan_filepath /path/to/taxid_plan.tsv \
     --analysis_output_filepath /path/to/analysis_output
   ```

# Requirements

- Python packages (see root README)
- Conda environments configured in `params.json`
