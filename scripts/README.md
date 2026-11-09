# EsViritu Simulation Scripts

Two scripts run EsViritu on simulated datasets and evaluate its detections against
the simulation ground truth:

| Script | Role |
|---|---|
| `run_esviritu_sim.sh` | Discovers simulated datasets, runs EsViritu on each, verifies outputs, then invokes the evaluator. |
| `evaluate_esviritu_sim.py` | Resolves EsViritu identifiers to NCBI taxids, classifies each prediction as TP/FP/FN/CROSSHIT, and writes metrics, tables, plots, and a log. |

---

## Expected Input Layout

```
<sim_root>/
└── <dataset>/
    ├── fastq/
    │   ├── <dataset>_R1.fq.gz        # or <dataset>_1.fq.gz, or *R1*.fq.gz
    │   └── <dataset>_R2.fq.gz
    └── input/
        └── <dataset>.tsv             # ground truth
```

The ground-truth table is read column-name-wise. The taxid is taken from the first
available of `taxid`, `TaxID`, `taxID`, `taxon`; if none is present, a numeric prefix
of the form `<digits>_` is extracted from the `path` / `fasta` / `assembly_file`
basename. The simulated assembly accession is taken from the first available of
`assembly_accession`, `accid`, `accession_id`, `accession`, and `Accession`.

When an input taxid is above species rank, it is rewritten to the species-rank taxid
of the simulated accession. The rewrite is applied only when the input taxid lies in
the resolved species lineage; otherwise the input taxid is kept. Inputs that resolve
to the same species taxid are counted as one ground-truth taxon.

EsViritu outputs are expected in the run output directory as
`<dataset>.detected_virus.info.tsv` (plus `<dataset>_esviritu.readstats.yaml`).

---

## `run_esviritu_sim.sh`

### Output Directory Layout

Defaults: `--output-dir <sim_root>/esviritu_sim`, `--eval-dir <output-dir>/evaluation`.

```
<output_dir>/
├── esviritu_sim_manifest.tsv
├── <dataset>.detected_virus.info.tsv
├── <dataset>_esviritu.readstats.yaml
├── <dataset>_esviritu_sim_driver.out
└── evaluation/
    ├── esviritu_sim_evaluation.log
    ├── esviritu_identifier_taxid_cache.json
    ├── esviritu_gt_accession_species_taxid_cache.json
    ├── esviritu_gt_taxid_lineage_cache.json
    ├── esviritu_sim_evaluation_detailed.tsv
    ├── esviritu_sim_evaluation_per_dataset.tsv
    ├── esviritu_sim_evaluation_overall.tsv
    ├── esviritu_sim_evaluation_aggregate.tsv
    ├── esviritu_sim_evaluation_by_order.tsv
    ├── esviritu_sim_evaluation_metadata.json
    └── *.png
```

### Files Written by the Driver

| File | Description |
|---|---|
| `<dataset>.detected_virus.info.tsv` | EsViritu per-virus detection table. The evaluator consumes `Assembly`, `Accession`, `read_count`, `covered_bases`, `mean_coverage`, and uses `species`/`order` as taxonomy fallbacks. |
| `<dataset>_esviritu.readstats.yaml` | EsViritu per-sample read statistics. Used only to confirm a run completed. |
| `<dataset>_esviritu_sim_driver.out` | Combined stdout/stderr of the EsViritu process for that dataset. |
| `esviritu_sim_manifest.tsv` | Append-only run ledger with columns `dataset`, `status` (`ok`/`fail`), `note` (`ok`, `exists`, `preflight`, `run`, `verify`, `unknown`), `started`, `finished`. |
| `evaluation/*` | Produced by `evaluate_esviritu_sim.py` (documented below). |

If EsViritu exits non-zero but both output tables exist, the dataset is accepted with a
warning. If the log contains `No reads aligned`, header-only EsViritu outputs are
created so the dataset still evaluates:

- `detected_virus.info.tsv` with the standard EsViritu header and no data rows.
- `_esviritu.readstats.yaml` containing `reads for EsViritu denominator: 0`.

With `--keep-temp`, EsViritu's `--keep` flag is passed through so temporary files are
retained. With `--dry-run`, no commands are executed and no outputs are written.

### Dataset Discovery and Skipping

Datasets come from `--datasets` (comma separated), `--datasets-file` (one per line),
or automatic discovery under `--sim-root` (directories named `fastq`, falling back to
`*_R1.fq.gz` / `*_1.fq.gz` files), filtered by `--dataset-glob`.

A dataset is considered done, and skipped unless `--force` is given, when both
`<dataset>.detected_virus.info.tsv` and `<dataset>_esviritu.readstats.yaml` are
non-empty. `--eval-only` instead discovers datasets from existing EsViritu tables in
the output directory and validates that both the EsViritu table and the ground-truth
table exist.

The script always runs a final full evaluation over all discovered datasets after the
per-dataset loop. That final run overwrites the evaluation tables/plots with the full
dataset set. Evaluation failures are suppressed by the driver (`|| true`); consult
`evaluation/esviritu_sim_evaluation.log` for errors.

---

## `evaluate_esviritu_sim.py`

### Identifier and Taxonomy Resolution

Ground-truth taxids are resolved from the simulated assembly accession. For each
dataset, unique accessions are looked up with
`metagenomics_utils.ncbi_tools.resolve_accession_species_taxids`, which returns the
nearest species-rank taxid. An above-species input taxid is rewritten to that species
taxid only when the input taxid is an ancestor of the resolved species taxid. Input
taxids that are already at or below species rank, unresolved accessions, and
non-descendant taxids are kept unchanged. Because the ground-truth set is keyed by the
resolved taxid, multiple inputs resolving to the same species are counted once.

Accession lookups are cached in
`evaluation/esviritu_gt_accession_species_taxid_cache.json`; lineage checks are cached
in `evaluation/esviritu_gt_taxid_lineage_cache.json`. Rewrite outcomes are counted by
reason in the evaluation metadata.

Each detection's taxid is resolved from its identifiers through NCBI. The EsViritu
`species` and `order` columns are retained as fallbacks for rows whose identifier
cannot be resolved:

1. `Accession` is tried first, then `Assembly`.
2. Each identifier is normalised with `metagenomics_utils.ncbi_tools.normalize_accession`
   and looked up with `retrieve_taxid_from_accession`.
3. Successful lookups are cached in
   `evaluation/esviritu_identifier_taxid_cache.json` under keys of the form
   `<database>:<accession>`, so repeated runs avoid network lookups.
4. A row whose identifier cannot be resolved gets `pred_taxid = None`; its EsViritu
   `species` and `order` values are used as fallbacks. If neither is usable, the row
   is counted as a false positive (FP).
5. Both ground-truth and predicted taxids are mapped to species-level match keys and
   orders through the NCBI taxonomy wrapper, using the repository-level `taxa.db` when
   it exists. The species-rank name is preferred. When a taxid has no species-rank
   lineage entry, the deepest available taxonomic name is used as the match key.
   Rank prefixes (for example `s__` or `o__`) are stripped and species names are
   lower-cased for comparison.

Ground-truth taxids without a species-rank entry can therefore still match through
the deepest-name fallback; only taxids with no usable match key are counted as false
negatives (FN).

### Classification Rule

For each dataset, with `--min-reads` (default `1.0`):

1. Predictions with `read_count < min_reads` are ignored entirely.
2. Ground-truth taxids are grouped by resolved species-level match key (species-rank
   name, or the deepest available name when no species rank exists).
3. A predicted row is a candidate match when its resolved species-level key equals a
   ground-truth key. All other rows are **FP**.
4. Within each matched species group, rows are ranked by
   `(covered_bases > 0, read_count, mean_coverage)` descending.
   - The best row is the match for the species.
   - Every additional row for that species is a **CROSSHIT**.
5. For every ground-truth taxid:
   - **TP** if its species-level key has a best matching row.
   - **FN** otherwise.

A crosshit is never counted as an FP and is excluded from the primary precision and
recall. When several ground-truth taxids resolve to the same species, the species is
counted once. FP records retain the predicted species and order rather than
discarding them.

### Status Semantics

| Status | Meaning |
|---|---|
| `TP` | Ground-truth taxid with a matching predicted species; the best row for that species is recorded. |
| `FN` | Ground-truth taxid with no matching predicted species-level key (or no usable taxonomy). |
| `FP` | Predicted row whose species-level key does not match any ground-truth taxon, including unresolved identifiers without an EsViritu fallback. |
| `CROSSHIT` | Additional predicted row for a species that already has a best (TP) row. Reported but excluded from precision/recall. |

### Metrics

Per dataset:

```
precision              = TP / (TP + FP)
precision_strict       = TP / (TP + FP + CROSSHIT)
recall                 = TP / (TP + FN)
f1                     = 2 * precision * recall / (precision + recall)
crosshit_rate_matched  = CROSSHIT / (TP + CROSSHIT)
crosshit_rate_all      = CROSSHIT / (TP + FP + CROSSHIT)
crosshits_per_tp       = CROSSHIT / TP
```

The primary `precision` treats crosshits as excluded. `precision_strict` also counts
crosshits in the denominator. `crosshit_rate_matched` uses matched predictions
(`TP + CROSSHIT`) as the denominator; `crosshit_rate_all` uses all predictions
(`TP + FP + CROSSHIT`).

- **Overall / micro metrics** pool the `detailed` records across all datasets before
  computing precision, precision_strict, recall, F1, and the crosshit rates.
- **Aggregate / macro metrics** are the distribution of the per-dataset precision,
  recall, and F1: count, mean, median, population standard deviation, minimum, and
  maximum.
- **Order-level metrics** recompute precision, recall, and F1 from the TP/FP/FN records
  grouped by taxonomic order. Crosshits are not included. FP records retain their
  predicted order; only records without a resolved order fall into the `unknown` group.

---

## Evaluation Outputs (`<eval_dir>/`)

### `esviritu_sim_evaluation.log`

Run log (also echoed to stderr). Records dataset discovery, per-dataset TP/FP/FN/
crosshit counts, unresolved identifiers, taxonomy failures, macro statistics, and the
final overall metrics.

### `esviritu_identifier_taxid_cache.json`

Persistent `database:accession -> taxid` cache. Written atomically; a malformed cache
is ignored with a warning and rebuilt.

### `esviritu_sim_evaluation_detailed.tsv`

One row per classified record.

| Column | Description |
|---|---|
| `dataset` | Dataset name. |
| `taxid` | Ground-truth taxid (semicolon-joined when several taxids share a species); empty for FP. |
| `order` | Taxonomic order (predicted order for FP records), or empty when unresolved. |
| `species` | Resolved, normalised species-level match key. |
| `assembly` | EsViritu `Assembly` identifier (empty for FN). |
| `read_count` | EsViritu read count (0 for FN). |
| `status` | `TP`, `FP`, `FN`, or `CROSSHIT`. |

### `esviritu_sim_evaluation_per_dataset.tsv`

One row per dataset.

| Column | Description |
|---|---|
| `dataset` | Dataset name. |
| `n_gt` | Number of ground-truth taxids. |
| `n_pred_matched` | Number of TP taxids (equal to `tp`). |
| `tp`, `fp`, `fn`, `crosshit` | Status counts. |
| `precision`, `recall`, `f1` | Per-dataset primary metrics. |
| `precision_strict` | Per-dataset precision with crosshits in the denominator. |
| `crosshit_rate_matched`, `crosshit_rate_all` | Per-dataset crosshit rates using matched and all predictions, respectively. |
| `crosshits_per_tp` | Per-dataset crosshits divided by TP (`0` when TP is zero). |

### `esviritu_sim_evaluation_overall.tsv`

Two columns, `metric` and `value`, with rows for `precision`, `precision_strict`,
`recall`, `f1`, `tp`, `fp`, `fn`, `crosshit`, `crosshit_rate_matched`,
`crosshit_rate_all`, and `crosshits_per_tp`. Values are pooled (micro) over all
datasets.

### `esviritu_sim_evaluation_aggregate.tsv`

One row per metric (`precision`, `recall`, `f1`) with columns `n_datasets`, `mean`,
`median`, `std` (population), `min`, and `max`, computed across per-dataset values.

### `esviritu_sim_evaluation_by_order.tsv`

One row per taxonomic order with columns `order`, `tp`, `fp`, `fn`, `precision`,
`recall`, and `f1`, sorted by order name.

### `esviritu_sim_evaluation_metadata.json`

```json
{
  "datasets": ["<dataset>", "..."],
  "min_reads": 1.0,
  "n_datasets": 1,
  "taxonomy_species_fallback": "deepest_available_name",
  "prediction_taxonomy_fallback": "esviritu_species_and_order",
  "gt_taxid_resolution": "accession_species_taxid",
  "gt_accession_columns": ["assembly_accession", "accid", "accession_id", "accession", "Accession"],
  "gt_accession_cache": "esviritu_gt_accession_species_taxid_cache.json",
  "gt_taxid_lineage_cache": "esviritu_gt_taxid_lineage_cache.json",
  "gt_rewrite_counts": {"rewritten": 0, "already_species": 0},
  "fp_orders_preserved": true
}
```

### Plots

Written only when `--plots` is given (the shell driver defaults to plots), `--no-plots`
is absent, matplotlib is importable, and at least one dataset was evaluated. Plot
failures are logged and do not abort the run.

| File | Description |
|---|---|
| `esviritu_sim_precision_recall_per_dataset.png` | Paired bars of per-dataset recall and precision. |
| `esviritu_sim_recall_by_order.png` | Recall per taxonomic order. |
| `esviritu_sim_aggregate_metrics.png` | Mean precision/recall/F1 across datasets with standard-deviation error bars. |
| `esviritu_sim_metric_distribution.png` | Boxplot of per-dataset precision (TP/(TP+FP)), strict precision (TP/(TP+CH+FP)), recall, and F1. Datasets with no matched predictions contribute 0 to precision and strict precision. |
| `esviritu_sim_tp_fp_fn_totals.png` | Stacked bar of pooled TP/FP/FN totals. |

### Conditional Outputs

Table files with no rows are not written (the writer skips empty tables). In practice:

- `overall` and `metadata` are always written.
- `detailed`, `per_dataset`, `aggregate`, and `by_order` are omitted when the run
  produced no corresponding records.

---

## Command-Line Reference

### `run_esviritu_sim.sh`

| Option | Default | Description |
|---|---|---|
| `--sim-root PATH` | required | Root containing one directory per simulated dataset. |
| `--output-dir PATH` | `<sim_root>/esviritu_sim` | EsViritu run outputs. |
| `--eval-dir PATH` | `<output_dir>/evaluation` | Evaluation outputs. |
| `--datasets LIST` | auto-discover | Comma-separated dataset names. |
| `--datasets-file PATH` | none | File with one dataset name per line. |
| `--dataset-glob PATTERN` | `*` | Filter for auto-discovered datasets. |
| `--esviritu-bin PATH` | `EsViritu` | EsViritu executable. |
| `--esviritu-db PATH` | `$ESVIRITU_DB` | EsViritu database directory. |
| `--conda-env NAME` | none | Conda environment to activate. |
| `--conda-bin PATH` | `conda` from PATH | Conda executable. |
| `--cpus N` | `8` | EsViritu threads. |
| `--mmp PRESET` | `sr` | EsViritu `-mmP` preset. |
| `--min-align-length N` | `100` | EsViritu minimum alignment length. |
| `--quality True\|False` | `False` | EsViritu `-q`. |
| `--filter-seqs True\|False` | `False` | EsViritu `-f`. |
| `--memory-max SIZE` | `24G` | Memory cap applied when `systemd-run` is available. |
| `--nice N` | none | Run EsViritu under `nice`. |
| `--eval-min-reads N` | `1` | Passed to the evaluator as `--min-reads`. |
| `--eval-only` | off | Evaluate existing outputs without running EsViritu. |
| `--plots` / `--no-plots` | plots on | Enable/disable evaluation plots. |
| `--keep-temp` | off | Pass `--keep` to EsViritu. |
| `--force` | off | Re-run datasets that already have outputs. |
| `--stop-on-error` | off | Stop the dataset loop on the first failure. |
| `--dry-run` | off | Print intended actions only. |
| `-v\|--verbose` | off | Pass verbose logging to the evaluator. |

### `evaluate_esviritu_sim.py`

| Option | Default | Description |
|---|---|---|
| `--esviritu-dir PATH` | required | Directory containing `<dataset>.detected_virus.info.tsv`. |
| `--sim-root PATH` | required | Root containing ground-truth inputs. |
| `--output-dir PATH` | required | Evaluation output directory. |
| `--datasets LIST` | auto-discover | Comma-separated dataset names. |
| `--datasets-file PATH` | none | File with one dataset name per line. |
| `--min-reads N` | `1.0` | Ignore prediction rows below this read count. |
| `--plots` | off | Write PNG plots when matplotlib is available. |
| `--no-plots` | off | Force plots off. |
| `--verbose` | off | Enable debug-level logging. |

### Example

```bash
# Run EsViritu and evaluate everything under a simulation root
./scripts/run_esviritu_sim.sh \
  --sim-root /path/to/simulations \
  --esviritu-db /path/to/esviritu_db \
  --conda-env esviritu

# Re-evaluate existing EsViritu outputs only
./scripts/run_esviritu_sim.sh \
  --sim-root /path/to/simulations \
  --eval-only \
  --plots
```
