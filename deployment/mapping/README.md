# Reference-truth mapping workflow

Maps a dataset's simulated reads directly against the references that were used
to simulate it — no classification, no clustering, no clade report. It exists to
study **mapping recall** (can the reads reach their own reference at all?) and
**recall recovery** (how much of the production pipeline's missed recall is
attributable to mapping vs. classification/assembly selection).

## Why this is faithful to the production pipeline

* The reference set is reconstructed with the **same retrieval** used by
  `deployment/simulation/simulate.nf` (`retrieve` with no `--include_term` /
  `--exclude_term`), so the references are exactly the ones wgsim simulated
  from.
* Reads are mapped with the **same mapping / filtering / coverage parameters**
  as `deployment/classify/classify.nf`
  (`minimap2_illumina_params`, `msamtools_params`, `samtools coverage`).
* Only the classification, clustering and clade-matching steps are removed.

## Inputs

| Parameter | Description |
|-----------|-------------|
| `input_table` | Study table (`taxid`, `reads`, ...). Same table fed to simulation. |
| `reads` | Directory holding the dataset's simulated `*_R1.fq.gz` / `*_R2.fq.gz`. |
| `output_dir` | Study output root; each dataset writes to `<output_dir>/<analysis_id>/mapping/`. |
| `analysis_id` | Dataset name. |
| `assembly_store` | Assembly store used by `retrieve`. |
| `references_extract_script` | `reference_management/main.py`. |
| `python_bin` | Python interpreter for `retrieve` and the statistics script. |
| `minimap2_illumina_params`, `msamtools_params` | Mapping parameters (must match classification). |
| `*_conda_env` | Conda environments for the mapping processes. |

## Outputs

Written inside each dataset (non-destructive; never touches `output/`):

```
<output_dir>/<dataset>/mapping/
├── matched_assemblies.tsv
├── merged_coverage_statistics.tsv
├── mapping_statistics.tsv                 # one row per dataset
└── mapping_statistics_per_reference.tsv   # one row per reference
```

`mapping_statistics.tsv` reports two recall definitions, so both are available
for downstream analysis:

* `mapping_recall_numreads` — fraction of input taxids with `numreads > 0`.
* `mapping_recall_coverage` — fraction of input taxids with `coverage > 0`.

plus reference counts, reads simulated/mapped, and coverage summaries.
`n_taxids_no_assembly` captures input taxids for which `retrieve` found no
assembly (they can never be recalled — a hard ceiling on any mapping recall).

The study-level roll-up is produced by `deploy_map_references.sh` into
`<output_dir>/mapping_statistics_all.tsv`.

## Usage

```bash
./deployment/workflows/deploy_map_references.sh
```

Or run a single dataset directly:

```bash
nextflow run deployment/mapping/map_references.nf \
    -profile conda \
    -params-file deployment/mapping/params.json \
    --input_table "$OUTPUT_DIR/$DS/input/$DS.tsv" \
    --reads "$OUTPUT_DIR/$DS/fastq" \
    --output_dir "$OUTPUT_DIR" \
    --analysis_id "$DS"
```
