"""
Diagnose where input taxa are lost before they can be recalled.

For every ``(dataset, input taxid)`` pair the script assigns exactly one bucket,
using only the flat files produced by the study pipeline (no lineage lookups,
no network access):

======================  ============================================================
bucket                  meaning
======================  ============================================================
absent_classification   taxid never appears in ``<ds>_merged_classification.tsv``
                        (classifier database gap / unmappable k-mers)
classified_no_assembly  classified, but no row in ``output/matched_assemblies.tsv``
                        (taxonomy -> genome lookup failed; reference-store gap or
                        accession/filename mismatch)
assembly_no_coverage    matched to an assembly, but that assembly has no row in
                        ``merged_coverage_statistics.tsv`` (never mapped)
assembly_zero_reads     coverage row exists but ``numreads == 0`` (mapping failed,
                        e.g. high mutation rate)
assembly_with_reads     reads mapped to the matched assembly -> recallable by the
                        m-stats pipeline (further losses are score/level gating)
======================  ============================================================

Usage::

    python -m deployment.model_evaluation.analysis_scripts.diagnose_recall_gap \
        --study-output /path/to/output_study --output-dir /path/to/diag [--limit N]
"""

from __future__ import annotations

import argparse
import logging
import os
import re
from collections.abc import Iterable

import pandas as pd

from deployment.model_evaluation.data_loader import passes_assembly_completeness

logger = logging.getLogger(__name__)

BUCKETS = [
    "absent_classification",
    "classified_no_assembly",
    "assembly_no_coverage",
    "assembly_zero_reads",
    "assembly_with_reads",
]

_ACCESSION_RE = re.compile(r"(GC[AF]_\d+\.\d+|[A-Z]{1,2}_?\d{5,}\.\d+)")


# ---------------------------------------------------------------------------
# Core (pure) logic
# ---------------------------------------------------------------------------


def _to_int_set(values: Iterable) -> set[int]:
    out: set[int] = set()
    for v in values:
        try:
            if pd.notna(v):
                out.add(int(v))
        except (TypeError, ValueError):
            continue
    return out


def _accession_key(value: str | None) -> str | None:
    """Normalise an accession or filename to a bare accession string."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    m = _ACCESSION_RE.search(str(value))
    return m.group(0) if m else str(value)


def classify_recall_gap(
    input_df: pd.DataFrame,
    classification_df: pd.DataFrame | None,
    matched_df: pd.DataFrame | None,
    coverage_df: pd.DataFrame | None,
) -> pd.DataFrame:
    """
    Bucket each input taxid according to how far it survives through the pipeline.

    Returns one row per unique input taxid with columns
    ``taxid, reads_simulated, mutation_rate, bucket, in_classification,
    classifier_uniq_reads, matched_accession, mapped_reads``.
    """
    if input_df is None or input_df.empty or "taxid" not in input_df.columns:
        return pd.DataFrame(columns=["taxid", "reads_simulated", "mutation_rate", "bucket"])

    inputs = input_df.dropna(subset=["taxid"]).copy()
    inputs["taxid"] = inputs["taxid"].astype(int)
    agg = {}
    if "reads" in inputs.columns:
        agg["reads_simulated"] = ("reads", "sum")
    if "mutation_rate" in inputs.columns:
        agg["mutation_rate"] = ("mutation_rate", "mean")
    per_taxid = inputs.groupby("taxid").agg(**agg).reset_index() if agg else inputs[["taxid"]].drop_duplicates()
    for col in ("reads_simulated", "mutation_rate"):
        if col not in per_taxid.columns:
            per_taxid[col] = float("nan")

    # classification evidence
    cls_reads: dict[int, float] = {}
    if classification_df is not None and not classification_df.empty and "taxid" in classification_df.columns:
        cls = classification_df.dropna(subset=["taxid"]).copy()
        cls["taxid"] = cls["taxid"].astype(int)
        if "uniq_reads" in cls.columns:
            cls["uniq_reads"] = pd.to_numeric(cls["uniq_reads"], errors="coerce").fillna(0)
            cls_reads = cls.groupby("taxid")["uniq_reads"].sum().to_dict()
        else:
            cls_reads = {t: float("nan") for t in cls["taxid"].unique()}

    # matched assemblies (taxid -> accession)
    matched_acc: dict[int, str] = {}
    if matched_df is not None and not matched_df.empty and "taxid" in matched_df.columns:
        m = matched_df.dropna(subset=["taxid"]).copy()
        m["taxid"] = m["taxid"].astype(int)
        acc_col = "assembly_accession" if "assembly_accession" in m.columns else None
        for _, row in m.iterrows():
            acc = _accession_key(row[acc_col]) if acc_col else None
            if acc is None and "assembly_file" in m.columns:
                acc = _accession_key(os.path.basename(str(row["assembly_file"])))
            matched_acc.setdefault(int(row["taxid"]), acc)

    # coverage stats (accession -> mapped reads)
    cov_reads: dict[str, float] = {}
    if coverage_df is not None and not coverage_df.empty:
        cov = coverage_df.copy()
        name_col = "#rname" if "#rname" in cov.columns else ("rname" if "rname" in cov.columns else None)
        reads_col = "numreads" if "numreads" in cov.columns else None
        if name_col is not None:
            for _, row in cov.iterrows():
                key = _accession_key(row[name_col])
                if key is None and "file" in cov.columns:
                    key = _accession_key(os.path.basename(str(row["file"])))
                if key is None:
                    continue
                n = float(row[reads_col]) if reads_col and pd.notna(row[reads_col]) else 0.0
                cov_reads[key] = cov_reads.get(key, 0.0) + n

    rows = []
    for _, r in per_taxid.iterrows():
        tid = int(r["taxid"])
        in_cls = tid in cls_reads
        uniq = cls_reads.get(tid, float("nan"))
        acc = matched_acc.get(tid)
        mapped = cov_reads.get(acc) if acc is not None else None

        if not in_cls:
            bucket = "absent_classification"
        elif acc is None:
            bucket = "classified_no_assembly"
        elif mapped is None:
            bucket = "assembly_no_coverage"
        elif mapped <= 0:
            bucket = "assembly_zero_reads"
        else:
            bucket = "assembly_with_reads"

        rows.append(
            {
                "taxid": tid,
                "reads_simulated": r["reads_simulated"],
                "mutation_rate": r["mutation_rate"],
                "bucket": bucket,
                "in_classification": int(in_cls),
                "classifier_uniq_reads": uniq,
                "matched_accession": acc,
                "mapped_reads": mapped if mapped is not None else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def summarise_buckets(per_taxid: pd.DataFrame) -> pd.DataFrame:
    """Aggregate bucket counts and shares (overall and per dataset if present)."""
    if per_taxid.empty:
        return pd.DataFrame(columns=["data_set", "bucket", "count", "share"])
    group_cols = ["data_set"] if "data_set" in per_taxid.columns else []
    counts = per_taxid.groupby(group_cols + ["bucket"]).size().rename("count").reset_index()
    totals = per_taxid.groupby(group_cols).size().rename("total").reset_index() if group_cols else None
    if totals is not None:
        counts = counts.merge(totals, on=group_cols)
    else:
        counts["total"] = len(per_taxid)
    counts["share"] = counts["count"] / counts["total"]
    return counts.drop(columns=["total"])


# ---------------------------------------------------------------------------
# File-system driver
# ---------------------------------------------------------------------------


def _read_tsv(path: str) -> pd.DataFrame | None:
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception as e:  # pragma: no cover - defensive
        logger.warning(f"Could not read {path}: {e}")
        return None


def diagnose_dataset(study_output: str, dataset: str) -> pd.DataFrame:
    base = os.path.join(study_output, dataset)
    input_df = _read_tsv(os.path.join(base, "input", f"{dataset}.tsv"))
    if input_df is None:
        return pd.DataFrame()
    classification_df = _read_tsv(os.path.join(base, "classification", f"{dataset}_merged_classification.tsv"))
    matched_df = _read_tsv(os.path.join(base, "output", "matched_assemblies.tsv"))
    coverage_df = _read_tsv(os.path.join(base, "output", "merged_coverage_statistics.tsv"))

    per_taxid = classify_recall_gap(input_df, classification_df, matched_df, coverage_df)
    per_taxid.insert(0, "data_set", dataset)
    per_taxid["has_classification_file"] = classification_df is not None
    per_taxid["has_matched_file"] = matched_df is not None
    per_taxid["has_coverage_file"] = coverage_df is not None
    return per_taxid


def diagnose_study(
    study_output: str,
    limit: int | None = None,
    require_complete_assemblies: bool = True,
    min_uniq_reads: int = 1,
    max_missing_pct: float = 5.0,
) -> pd.DataFrame:
    datasets = sorted(d for d in os.listdir(study_output) if os.path.isdir(os.path.join(study_output, d)))
    if require_complete_assemblies:
        datasets = [
            d
            for d in datasets
            if passes_assembly_completeness(study_output, d, min_uniq_reads=min_uniq_reads, max_missing_pct=max_missing_pct)
        ]
        logger.info(
            f"Assembly-completeness gate ON: restricting recall-gap diagnosis to {len(datasets)} complete datasets"
        )
    if limit:
        datasets = datasets[:limit]
    frames = []
    for i, ds in enumerate(datasets, 1):
        df = diagnose_dataset(study_output, ds)
        if not df.empty:
            frames.append(df)
        if i % 100 == 0:
            logger.info(f"Processed {i}/{len(datasets)} datasets")
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _plot_buckets(summary: pd.DataFrame, per_taxid: pd.DataFrame, output_dir: str) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover
        logger.warning("matplotlib not available; skipping plots")
        return

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    overall = summary.set_index("bucket").reindex(BUCKETS).fillna(0)
    axes[0].bar(overall.index, overall["share"], color="steelblue")
    axes[0].set_ylabel("share of input taxa")
    axes[0].set_title("Where input taxa are lost (all datasets)")
    axes[0].tick_params(axis="x", rotation=30)

    if "mutation_rate" in per_taxid.columns and per_taxid["mutation_rate"].notna().any():
        bins = [-0.001, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 1.0]
        tmp = per_taxid.copy()
        tmp["mut_bin"] = pd.cut(tmp["mutation_rate"], bins=bins)
        shares = tmp.groupby("mut_bin", observed=True)["bucket"].value_counts(normalize=True).unstack().fillna(0)
        shares = shares.reindex(columns=BUCKETS, fill_value=0)
        shares.plot(kind="bar", stacked=True, ax=axes[1], colormap="viridis")
        axes[1].set_xlabel("mutation rate")
        axes[1].set_ylabel("share")
        axes[1].set_title("Loss bucket by mutation rate")
        axes[1].legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1.0))
        axes[1].tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "recall_gap_buckets.png"), dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--study-output", required=True, help="Study output directory (one sub-dir per dataset)")
    parser.add_argument("--output-dir", required=True, help="Where to write the diagnostic TSVs/plot")
    parser.add_argument("--limit", type=int, default=None, help="Only process the first N datasets")
    parser.add_argument("--no-plot", action="store_true", help="Skip PNG generation")
    parser.add_argument(
        "--require-complete-assemblies",
        action="store_true",
        default=True,
        help="Restrict diagnosis to datasets passing the assembly-completeness gate (at most "
        "--max-missing-refs-pct% of detected references with uniq_reads>=--min-uniq-reads are "
        "unmatched). Default: on.",
    )
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help="Include datasets failing the assembly-completeness gate (equivalent to "
        "--no-require-complete-assemblies).",
    )
    parser.add_argument(
        "--min-uniq-reads", type=int, default=1, help="Minimum uniq_reads for a detected taxid. Default: 1."
    )
    parser.add_argument(
        "--max-missing-refs-pct",
        type=float,
        default=5.0,
        help="Max tolerated percent of qualified taxids lacking a matched assembly before a dataset "
        "fails the completeness gate (default: 5.0).",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    os.makedirs(args.output_dir, exist_ok=True)

    per_taxid = diagnose_study(
        args.study_output,
        limit=args.limit,
        require_complete_assemblies=args.require_complete_assemblies and not args.allow_incomplete,
        min_uniq_reads=args.min_uniq_reads,
        max_missing_pct=args.max_missing_refs_pct,
    )
    if per_taxid.empty:
        logger.error("No datasets with input files found")
        return 1

    per_taxid_path = os.path.join(args.output_dir, "recall_gap_per_taxid.tsv")
    per_taxid.to_csv(per_taxid_path, sep="\t", index=False)

    overall = summarise_buckets(per_taxid.drop(columns=["data_set"]))
    overall.to_csv(os.path.join(args.output_dir, "recall_gap_summary.tsv"), sep="\t", index=False)
    per_dataset = summarise_buckets(per_taxid)
    per_dataset.to_csv(os.path.join(args.output_dir, "recall_gap_per_dataset.tsv"), sep="\t", index=False)

    logger.info("Recall-gap buckets (share of all input taxa):")
    for _, row in overall.set_index("bucket").reindex(BUCKETS).dropna().reset_index().iterrows():
        logger.info(f"  {row['bucket']:<24} {int(row['count']):>7}  {row['share']:6.1%}")
    n_missing_cls = int((~per_taxid["has_classification_file"]).sum())
    n_missing_matched = int((~per_taxid["has_matched_file"]).sum())
    if n_missing_cls or n_missing_matched:
        logger.warning(
            f"Rows from datasets lacking files: classification={n_missing_cls}, matched_assemblies={n_missing_matched}"
        )

    if not args.no_plot:
        _plot_buckets(overall, per_taxid, args.output_dir)
    logger.info(f"Saved diagnostics to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
