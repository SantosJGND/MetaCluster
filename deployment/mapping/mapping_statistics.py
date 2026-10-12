#!/usr/bin/env python3
"""Reference-truth mapping statistics for the stripped-down mapping workflow.

For a single dataset this joins three sources:

  * the input table (``taxid``, ``reads`` -- reads simulated per reference),
  * ``matched_assemblies.tsv`` (``taxid``, ``assembly_accession``, ``assembly_file``),
  * ``merged_coverage_statistics.tsv`` (``file``, ``numreads``, ``coverage``, ...).

It writes two tables into ``--output-dir``:

  * ``mapping_statistics.tsv`` -- one row per dataset, with two mapping-recall
    definitions (recovered when ``numreads > 0`` and recovered when
    ``coverage > 0``).
  * ``mapping_statistics_per_reference.tsv`` -- one row per reference, with the
    per-reference mapping / coverage values used to build the summary.

Coverage rows are matched to references by accession encoded in the coverage
``file`` name (``<sample>_<taxid>_<accession>.fna.sorted.coverage.txt``), the
same convention used by the production pipeline.
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

TAXID_CANDIDATES = ("taxid", "TaxID", "taxID", "taxon")
READS_CANDIDATES = ("reads", "num_reads", "read_count", "n_reads", "numreads")


def get_args():
    parser = argparse.ArgumentParser(description="Compute reference-truth mapping statistics for one dataset.")
    parser.add_argument("--input-table", required=True, help="Study/dataset input table (taxid + reads).")
    parser.add_argument("--matched-assemblies", required=True, help="matched_assemblies.tsv from `retrieve`.")
    parser.add_argument("--coverage-statistics", required=True, help="merged_coverage_statistics.tsv.")
    parser.add_argument("--dataset", required=True, help="Dataset / analysis id.")
    parser.add_argument("--output-dir", default=".", help="Directory where the statistics tables are written.")
    return parser.parse_args()


def _detect_column(df, candidates, context):
    for name in candidates:
        if name in df.columns:
            return name
    raise ValueError(f"{context}: none of {list(candidates)} found in columns {list(df.columns)}")


def _normalize_taxid(value):
    """Return a stable string key for a taxonomic id (handles ints/floats/strings)."""
    if value is None or (isinstance(value, float) and np.isnan(value)) or pd.isna(value):
        return None
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value).strip()


def _normalize_coverage(df):
    """Normalize coverage columns to a known set; aggregate per reference file."""
    df = df.rename(columns={"#rname": "rname", "reference": "rname"})

    for col in ("numreads", "coverage"):
        if col not in df.columns:
            raise ValueError(f"coverage statistics missing required column '{col}' (columns: {list(df.columns)})")

    # Aggregate contigs belonging to the same reference file (mirrors node_stats).
    agg = {"numreads": "sum", "coverage": "max"}
    for col, how in (("covbases", "sum"), ("meanmapq", "max"), ("error_rate", "max")):
        if col in df.columns:
            agg[col] = how

    if "file" not in df.columns:
        raise ValueError("coverage statistics missing required 'file' column")

    df = df.groupby("file", as_index=False).agg(agg)
    for col in ("covbases", "meanmapq", "error_rate"):
        if col not in df.columns:
            df[col] = np.nan
    return df


def _match_coverage(coverage_df, taxid, accession):
    """Return the coverage slice for a reference, keyed by taxid+accession in the file name."""
    files = coverage_df["file"].astype(str)
    if taxid is not None and accession:
        pattern = f"_{taxid}_{accession}.fna."
        subset = coverage_df[files.str.contains(re.escape(pattern), regex=True, na=False)]
        if not subset.empty:
            return subset
    if accession:
        subset = coverage_df[files.str.contains(re.escape(str(accession)), regex=True, na=False)]
        if not subset.empty:
            return subset
    return coverage_df.iloc[0:0]


def main():
    args = get_args()
    os.makedirs(args.output_dir, exist_ok=True)

    input_df = pd.read_csv(args.input_table, sep="\t")
    taxid_col = _detect_column(input_df, TAXID_CANDIDATES, "input table")
    reads_col = _detect_column(input_df, READS_CANDIDATES, "input table")

    input_df = input_df.copy()
    input_df["_taxid"] = input_df[taxid_col].apply(_normalize_taxid)
    input_df["_reads"] = pd.to_numeric(input_df[reads_col], errors="coerce").fillna(0.0)

    reads_by_taxid = input_df.dropna(subset=["_taxid"]).groupby("_taxid")["_reads"].sum()
    input_taxids = list(reads_by_taxid.index)
    reads_simulated_total = float(reads_by_taxid.sum())

    matched = pd.read_csv(args.matched_assemblies, sep="\t")
    coverage = _normalize_coverage(pd.read_csv(args.coverage_statistics, sep="\t"))

    records = []
    if {"taxid", "assembly_accession"}.issubset(matched.columns):
        for _, row in matched.iterrows():
            taxid = _normalize_taxid(row.get("taxid"))
            accession = row.get("assembly_accession")
            if pd.isna(accession):
                continue
            accession = str(accession)
            assembly_file = row.get("assembly_file")
            assembly_file = None if pd.isna(assembly_file) else str(assembly_file)

            subset = _match_coverage(coverage, taxid, accession)
            numreads = float(subset["numreads"].sum()) if not subset.empty else 0.0
            cov = float(subset["coverage"].max()) if not subset.empty else 0.0
            covbases = float(subset["covbases"].sum()) if "covbases" in subset and not subset.empty else 0.0
            meanmapq = float(subset["meanmapq"].max()) if "meanmapq" in subset and not subset.empty else np.nan
            error_rate = float(subset["error_rate"].max()) if "error_rate" in subset and not subset.empty else np.nan

            records.append(
                {
                    "data_set": args.dataset,
                    "taxid": taxid,
                    "assembly_accession": accession,
                    "assembly_file": assembly_file,
                    "reads_simulated": float(reads_by_taxid.get(taxid, 0.0)),
                    "numreads_mapped": numreads,
                    "coverage": cov,
                    "covbases": covbases,
                    "meanmapq": meanmapq,
                    "error_rate": error_rate,
                    "recovered_numreads": numreads > 0,
                    "recovered_coverage": cov > 0,
                }
            )

    per_reference = pd.DataFrame(
        records,
        columns=[
            "data_set",
            "taxid",
            "assembly_accession",
            "assembly_file",
            "reads_simulated",
            "numreads_mapped",
            "coverage",
            "covbases",
            "meanmapq",
            "error_rate",
            "recovered_numreads",
            "recovered_coverage",
        ],
    )

    n_input_refs = len(input_taxids)
    if not per_reference.empty:
        per_taxid = per_reference.groupby("taxid").agg(
            recovered_numreads=("recovered_numreads", "any"),
            recovered_coverage=("recovered_coverage", "any"),
        )
    else:
        per_taxid = pd.DataFrame(columns=["recovered_numreads", "recovered_coverage"])

    n_taxids_with_assembly = int(per_reference["taxid"].nunique()) if not per_reference.empty else 0
    n_recovered_numreads = int(per_taxid["recovered_numreads"].sum()) if not per_taxid.empty else 0
    n_recovered_coverage = int(per_taxid["recovered_coverage"].sum()) if not per_taxid.empty else 0

    if not per_reference.empty:
        positive = per_reference[per_reference["coverage"] > 0]
        mean_coverage_all = float(per_reference["coverage"].mean())
        median_coverage_all = float(per_reference["coverage"].median())
        mean_coverage_mapped = float(positive["coverage"].mean()) if not positive.empty else 0.0
    else:
        mean_coverage_all = median_coverage_all = mean_coverage_mapped = 0.0

    def _ratio(num, den):
        return float(num) / float(den) if den else 0.0

    summary = pd.DataFrame(
        [
            {
                "data_set": args.dataset,
                "n_input_refs": n_input_refs,
                "n_reference_files": int(len(per_reference)),
                "n_taxids_with_assembly": n_taxids_with_assembly,
                "n_taxids_no_assembly": max(n_input_refs - n_taxids_with_assembly, 0),
                "n_recovered_numreads": n_recovered_numreads,
                "n_recovered_coverage": n_recovered_coverage,
                "mapping_recall_numreads": _ratio(n_recovered_numreads, n_input_refs),
                "mapping_recall_coverage": _ratio(n_recovered_coverage, n_input_refs),
                "reads_simulated_total": reads_simulated_total,
                "reads_mapped_total": float(per_reference["numreads_mapped"].sum()) if not per_reference.empty else 0.0,
                "mean_coverage": mean_coverage_all,
                "median_coverage": median_coverage_all,
                "mean_coverage_mapped": mean_coverage_mapped,
            }
        ]
    )

    per_reference.to_csv(os.path.join(args.output_dir, "mapping_statistics_per_reference.tsv"), sep="\t", index=False)
    summary.to_csv(os.path.join(args.output_dir, "mapping_statistics.tsv"), sep="\t", index=False)

    print(
        f"[{args.dataset}] input_refs={n_input_refs} "
        f"recall(numreads)={summary['mapping_recall_numreads'].iloc[0]:.3f} "
        f"recall(coverage)={summary['mapping_recall_coverage'].iloc[0]:.3f}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
