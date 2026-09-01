# Upstream Change Requests — Cohort / Dataset-Count Metadata

**Status:** Requested (not yet implemented upstream)
**Owner:** `metagenomics-evaluation-pipeline` maintainers
**Context for the upstream agent:** These changes are prerequisites for downstream
manuscript scripts (`writing/results/scripts/`) to compute summary statistics with
correct denominators. The downstream scripts were written before dataset counts
(analysed / empty / skipped / failed) were logged; they now need that metadata to
report honest statistics. This file documents exactly what is needed and **why**.

---

## Background

Raw simulation outputs contain datasets that are sometimes **empty** (no mapped
reads), **skipped**, or **fail to process**. Every summary statistic depends on
which datasets are in the denominator:

- a **failed** dataset produces no result → it must be _discarded_ (never counted);
- a **skipped / empty** dataset genuinely has recall 0 → it should be _excluded from
  per-dataset distribution means_ but _counted in an attempted-cohort lower bound_;
- the number of datasets that actually yielded per-dataset metrics is the
  denominator for all means.

Two producers write dataset-level outputs, and **only one of them logs cohort
counts**:

| Producer                                                 | Writes                                      | Logs cohort counts?       |
| -------------------------------------------------------- | ------------------------------------------- | ------------------------- |
| `deployment/model_evaluation/evaluate.py`                | `full_analysis/<run>/test_datasets_*.tsv`   | ✔ `pipeline_metadata.tsv` |
| `deployment/model_evaluation/analysis_data_extractor.py` | `<domain>_eda/per_dataset_metrics.tsv` etc. | ✘ **none**                |

The EDA producer covers the **full 2000-dataset study cohort**, while `evaluate.py`
runs the **600-dataset train/test split**. Their denominators are different cohorts
and must not be conflated. The missing piece is the EDA producer's cohort log.

---

## Request 1 — `analysis_data_extractor.py`: emit a cohort metadata file

### What is requested

Alongside `per_dataset_metrics.tsv`, write a `pipeline_metadata.tsv` in the same
output directory, two columns `metric` / `value`, with rows:

| metric             | value (example)   | meaning                                                  |
| ------------------ | ----------------- | -------------------------------------------------------- |
| `total_attempted`  | 2000              | datasets scanned from `--study-output`                   |
| `extracted`        | 1650              | rows written to `per_dataset_metrics.tsv` (= successful) |
| `dropped`          | 350               | `total_attempted − extracted`                            |
| `failed`           | 0                 | datasets that raised an exception during extraction      |
| `skipped`          | 350               | datasets with no usable m-stats / no reads (recall 0)    |
| `skipped_datasets` | `dataset_...;...` | semicolon-separated names (omit if none)                 |
| `failed_datasets`  | (messages)        | error messages (omit if none)                            |

### Why it is requested

The downstream EDA scripts compute statistics whose denominator is the number of
extracted datasets, but **nothing records how many datasets were dropped or why**.
Without `total_attempted` / `extracted` / `dropped`, it is impossible to:

- state the denominator unambiguously in the manuscript Methods;
- distinguish _empty/skipped_ datasets (genuine recall 0) from _failed_ datasets
  (must be discarded) at extraction time;
- reproduce which datasets contributed to the means.

The `evaluate.py` producer already emits `pipeline_metadata.tsv`
(`total_datasets`, `successful`, `failed`, `skipped_count`, `skipped_datasets`);
the EDA producer should mirror the same schema **scoped to the full study cohort**
so both cohort levels are documented.

### Exact contract

1. `extracted == len(per_dataset_metrics.tsv) rows`.
2. `dropped == total_attempted − extracted`.
3. `dropped == skipped + failed + (study gaps, if any)` — document the breakdown.
4. `skipped` = datasets with zero usable m-stats / no mapped reads (recall 0).
5. `failed` = datasets that raised an exception; **excluded** from all statistics.
6. File written _atomically_ (same writer as `evaluate.py`'s metadata writer).

---

## Request 2 — `evaluate.py`: already correct; confirm + extend JSON

`pipeline_metadata.tsv` already records `total_datasets`, `successful`, `failed`,
`skipped_count`, `skipped_datasets`, `failed_datasets` — this matches the contract
and is already consumed downstream.

### What is requested

Confirm that the JSON metadata block
(`evaluation_results.json` → `metadata`) and the TSV report the **same** counts,
and that `successful == len(test_datasets_summary_results.tsv) − 1`.
If `skipped_count` / `failed` are not present in the JSON, add them for parity.

### Why

Downstream pipeline scripts parse `pipeline_metadata.tsv` (TSV is authoritative).
This is a consistency check so the two serializations never diverge.

---

## Dependency map — which local scripts consume which metadata

| Local script file                                                                                                                         | Metadata consumed                                    | Purpose                                                                                          |
| ----------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
| `bacteria_pipeline_recall.py`, `bacteria_pipeline_summary.py`, `viral_pipeline_recall.py`, `viral_pipeline_summary.py`                    | `full_analysis/<run>/pipeline_metadata.tsv`          | use `successful` as the `N` denominator; assert row count == `successful`; report failed/skipped |
| `verify_bacteria_aggregates.py`                                                                                                           | `full_analysis/<run>/pipeline_metadata.tsv`          | derive expected test-set `N` and cohort-sum invariant instead of hard-coding                     |
| `bacteria_dataset_characteristics.py`, `dataset_characteristics.py`, `classifier_summary.py`, `combine_tables.py`, `hit_type_analysis.py` | **`<domain>_eda/pipeline_metadata.tsv` (Request 1)** | report `n_extracted`, `total_attempted`, `dropped`; state denominator as extracted-N             |
| `verify_viral_aggregates.py`                                                                                                              | `full_analysis` metadata (optional, parity)          | cohort-sum invariant                                                                             |

Any `pipeline_metadata.tsv` read must be **defensive**: if the file is absent
(e.g. extraction not yet re-run), local scripts fall back to row-count `N` and
emit a note. This lets downstream deploy before regenerating.

---

## Acceptance criteria (for the upstream agent)

1. Running `analysis_data_extractor.py` after this change produces
   `<domain>_eda/pipeline_metadata.tsv` with all rows from Request 1, where
   `extracted` matches the row count of the corresponding `per_dataset_metrics.tsv`.
2. `evaluate.py`'s `pipeline_metadata.tsv` and `evaluation_results.json` report
   identical `successful` / `failed` / `skipped_count`.
3. Neither change alters any existing output besides adding metadata / JSON fields.
4. Re-running both producers regenerates metadata that is internally consistent.
