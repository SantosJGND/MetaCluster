## Classification and Clustering Workflow

This directory contains a unified workflow for classifying metagenomic reads and clustering the results.

# Overview

The classify workflow:
1. Reads input FASTQ files (compatible with output from `deployment/simulation/simulate.nf`)
2. Optionally performs QC filtering
3. Runs classifiers (configurable: Centrifuge, Kraken2, Diamond, KrakenUnique)
4. Each classifier output is annotated with a `software` column to track source
5. Merges all classifier results
6. Maps reads to references
7. Clusters mapped reads

# Running the Classification

## Basic Usage

```bash
nextflow run classify.nf -profile conda \
  --reads /path/to/fastq \
  --output_dir /path/to/output
```

This runs with default parameters:
- `qc = true` (enable read filtering)
- `centrifuge = true` (enable Centrifuge)
- `kraken2 = true` (enable Kraken2)
- `diamond = false` (disable Diamond)
- `krakenunique = false` (disable KrakenUnique)

## Options

### QC Filtering

Enable or disable read quality control:

```bash
# With QC (default)
nextflow run classify.nf -profile conda --qc true

# Without QC
nextflow run classify.nf -profile conda --qc false
```

### Classifiers

Enable or disable specific classifiers:

```bash
# Only Kraken2
nextflow run classify.nf -profile conda --centrifuge false --kraken2 true

# All classifiers
nextflow run classify.nf -profile conda \
  --centrifuge true \
  --kraken2 true \
  --diamond true \
  --krakenunique true

# Default (Centrifuge + Kraken2)
nextflow run classify.nf -profile conda
```

### Combined Options

```bash
# Full configuration
nextflow run classify.nf -profile conda \
  --reads /path/to/fastq \
  --output_dir /path/to/output \
  --qc true \
  --centrifuge true \
  --kraken2 true \
  --diamond true \
  --krakenunique true
```

# Output Structure

```
{output_dir}/
└── {analysis_id}/
    ├── classification/
    │   ├── kraken2/
    │   ├── centrifuge/
    │   ├── diamond/          (if --diamond true)
    │   └── krakenunique/     (if --krakenunique true)
    ├── {analysis_id}_merged_classification.tsv   # Contains 'software' column
    └── ...
```

# Strict Assembly Validation Gate

Before reads are mapped, `ExtractReferenceSequences` must find an assembly for
every classified taxid with `uniq_reads >= min_uniq_reads`. Lookup order:
local assembly store, then one round of NCBI fallback. Taxids that cannot be
matched after both attempts fail the dataset:

- `--fail_on_missing_assemblies true` (default): the workflow errors, so the
  dataset produces no output and is treated as **incomplete** (its name is
  appended to `{output_dir}/incomplete_datasets.tsv` by
  `deployment/workflows/deploy_classify_map.sh`).
- `--fail_on_missing_assemblies false`: missing assemblies are reported in
  `unmatched_classified_taxids.tsv` but the dataset still maps what it can
  (diagnostic / escaping hatch).

Incomplete datasets are **excluded from the final analysis**: the analysis
extractor reports them as an explicit `incomplete` bucket in
`pipeline_metadata.tsv` and never counts them in recall/precision denominators.

Config (`params.json`):
```json
"min_uniq_reads": 1,
"fail_on_missing_assemblies": true
```

# Software Column

Each classifier adds a `software` column to identify the source:

| Classifier | software value |
|------------|----------------|
| Centrifuge | `centrifuge` |
| Kraken2 | `kraken2` |
| Diamond | `diamond` |
| KrakenUnique | `krakenunique` |

This allows tracking which classifier(s) identified each taxid in the merged results.

# Complete Workflow

Deployment is split into two independent, resumable steps so a validation gate
runs before mapping (see `deployment/workflows/deploy_simulate.sh` and
`deployment/workflows/deploy_classify_map.sh`).

1. **Simulate reads (`deploy_simulate.sh`):**

   ```bash
   nextflow run deployment/simulation/simulate.nf -profile conda \
     --input_table /path/to/input_table.tsv \
     --output_dir /path/to/simulation_output
   ```

2. **Classify, validate assemblies, and cluster (`deploy_classify_map.sh`):**

   ```bash
   nextflow run deployment/classify/classify.nf -profile conda \
     --reads /path/to/simulation_output/{dataset_name}/fastq \
     --output_dir /path/to/classification_output
   ```

3. **Evaluate results:**
   ```bash
   python deployment/model_evaluation/evaluate.py \
     --study_output_filepath /path/to/classification_output \
     --taxid_plan_filepath /path/to/taxid_plan.tsv \
     --analysis_output_filepath /path/to/analysis_output
   ```

# Requirements

- Python packages (see root README)
- Conda environments configured in `params.json`
- Reference databases:
  - Kraken2 index
  - Centrifuge index
  - (Optional) Diamond index
  - (Optional) KrakenUnique index
