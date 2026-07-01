#!/usr/bin/env python3
"""
Sample analysis script using trained models to analyze new metagenomic data.

Analyzes clustering output from map_clustering pipeline and generates predictions
using trained RecallModeller and CompositionModeller.

Usage:
    python analyze_samples.py --samples-dir /path/to/samples --training-dir /path/to/training
"""

import argparse
import json
import logging
import os
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper
from metagenomics_utils.overlap_manager import OverlapManager
from metagenomics_utils.overlap_manager.node_stats import (
    dataframe_update_with_lineage,
)
from metagenomics_utils.overlap_manager.om_models import (
    CutoffRecallModeller,
    GPCLFRecallModeller,
    RecallModeller,
    XGBCompositionModeller,
    predict_data_set_clades_composition,
)

logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def _load_recall_modeller(model_dir: str, data_set_divide: int) -> RecallModeller:
    """
    Load the correct recall modeller based on which bundle file exists.

    Inspects the model directory for known bundle filenames and instantiates
    the appropriate subclass of RecallModeller, keeping the call site unified
    via predict_cutoff().

    Returns:
        An instance of a RecallModeller subclass (GPCLFRecallModeller,
        CutoffRecallModeller, or base RecallModeller).

    Raises:
        FileNotFoundError: If no recognised bundle file is found.
    """
    bundle_map = {
        "recall_gp_clf_pipeline.pkl": (GPCLFRecallModeller, {"data_set_divide": data_set_divide}),
        "cutoff_recall_bundle.pkl": (CutoffRecallModeller, {"data_set_divide": data_set_divide, "target_recall": 1.0}),
        "recall_xgb_bundle.pkl": (RecallModeller, {"data_set_divide": data_set_divide}),
    }

    for bundle, (cls, kwargs) in bundle_map.items():
        bundle_path = os.path.join(model_dir, bundle)
        if os.path.exists(bundle_path):
            modeller = cls(**kwargs)
            modeller.load_model(model_dir)
            logger.info(f"Loaded {cls.__name__} from {bundle}")
            return modeller

    available = [f for f in os.listdir(model_dir) if f.endswith("_bundle.pkl")]
    raise FileNotFoundError(
        f"No recognised recall bundle in {model_dir}. Expected one of: {list(bundle_map)}. Found: {available}"
    )


def _add_best_match_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add best match columns to a dataframe."""
    df = df.copy()
    df["best_match_score"] = 0
    df["best_match_taxid"] = df["taxid"]
    df["best_match_is_best"] = False
    return df


def process_sample(
    sample: str,
    results_dir: Path,
    model_dir: Path,
    ncbi_wrapper: NCBITaxonomistWrapper,
    taxids_to_use: pd.DataFrame,
    tax_level: str,
    output_db: str,
    target_recall: float = 0.95,
    data_set_divide: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Process a single sample and return predictions and cluster metrics.

    Args:
        sample: Sample identifier.
        results_dir: Path to sample results directory.
        model_dir: Path to trained models.
        ncbi_wrapper: NCBI taxonomy wrapper.
        taxids_to_use: DataFrame of taxids to use for analysis.
        tax_level: Taxonomic level for analysis.
        output_db: Path to taxonomy database.
        target_recall: Target recall threshold (default 0.95).
        data_set_divide: Number of recall divisions used in model (default 5).

    Returns:
        Tuple of (predictions DataFrame, pruned predictions DataFrame, cluster metrics dict)
    """
    merged_classification_dir = results_dir / "classification"
    merged_classification_files = list(merged_classification_dir.glob("*merged_classification.tsv"))

    matched_assemblies = results_dir / "output" / "matched_assemblies.tsv"

    if not merged_classification_files or not matched_assemblies.exists():
        logger.error(f"Missing files for sample {sample}")
        return pd.DataFrame(), pd.DataFrame(), {}

    merged_classification = merged_classification_files[0]
    merged_classification_df = pd.read_csv(merged_classification, sep="\t")
    matched_assemblies_df = pd.read_csv(matched_assemblies, sep="\t")

    if not matched_assemblies_df["taxid"].isin(merged_classification_df["taxid"]).all():
        logger.error(f"Not all taxids in matched assemblies are present in merged classification for sample {sample}")
        return pd.DataFrame(), pd.DataFrame(), {}

    ncbi_wrapper.resolve_lineages(taxids_to_use["taxid"].tolist())

    composition_modeller = XGBCompositionModeller()
    composition_modeller.load_model(os.path.join(model_dir))

    recall_modeller = _load_recall_modeller(os.path.join(model_dir), data_set_divide)

    overlap_manager = OverlapManager(os.path.join(results_dir, "clustering"), max_proportion=1.0)
    if not overlap_manager.check_data_available():
        logger.error(f"No overlapping data available for sample {sample}")
        return pd.DataFrame(), pd.DataFrame(), {}

    result_taxids = overlap_manager.m_stats_matrix["taxid"].unique().tolist()

    ncbi_wrapper.resolve_lineages(result_taxids)

    m_stats_stats_matrix = dataframe_update_with_lineage(overlap_manager.m_stats_matrix, ncbi_wrapper=ncbi_wrapper)
    m_stats_stats_matrix = _add_best_match_columns(m_stats_stats_matrix)

    raw_pred = recall_modeller.predict(m_stats_stats_matrix)
    target_percentile = recall_modeller.predict_cutoff(m_stats_stats_matrix, target_recall=target_recall)

    if hasattr(raw_pred, "shape") and raw_pred.ndim == 2:
        raw_pred_df = pd.DataFrame(raw_pred, columns=recall_modeller.RecP_target_cols)
        last_best_match_relindex = raw_pred_df.iloc[0].get("last_best_match_relindex", 0.0)
    else:
        raw_pred_df = pd.DataFrame(raw_pred, columns=["k_min"])
        last_best_match_relindex = 0.0

    keep_index = round(target_percentile * m_stats_stats_matrix.shape[0])
    if keep_index == 0:
        keep_index = 1
    logger.info(
        f"{sample}: keeping top {keep_index} predictions (target_recall={target_recall}, percentile={target_percentile:.2f})"
    )

    # Calculate coverage proportions for filtering analysis
    original_matrix = overlap_manager.original_m_stats_matrix
    if original_matrix is not None and not original_matrix.empty:
        total_with_coverage = (original_matrix["coverage"] > 0).sum()
        kept_matrix = original_matrix.head(keep_index)
        kept_with_coverage = (kept_matrix["coverage"] > 0).sum() if not kept_matrix.empty else 0
        filtered_with_coverage = total_with_coverage - kept_with_coverage

        prop_coverage_above = kept_with_coverage / total_with_coverage if total_with_coverage > 0 else 0
        prop_coverage_below = filtered_with_coverage / total_with_coverage if total_with_coverage > 0 else 0
        resulting_percentile = keep_index / original_matrix.shape[0] if original_matrix.shape[0] > 0 else 0
    else:
        prop_coverage_above = 0
        prop_coverage_below = 0
        resulting_percentile = 0

    pruned_overlap_manager = OverlapManager(os.path.join(results_dir, "clustering"), max_proportion=target_percentile)

    logger.info(
        f"[{sample}] Pruned: {len(pruned_overlap_manager.leaves)} leaves, {len(pruned_overlap_manager.m_stats_matrix)} matrix rows"
    )

    results_pred = predict_data_set_clades_composition(
        sample,
        m_stats_stats_matrix=m_stats_stats_matrix,
        overlap_manager=overlap_manager,
        modeller=composition_modeller,
        input_taxa=taxids_to_use,
        tax_level=tax_level,
    )

    results_pred_sample_df = results_pred.merge(
        m_stats_stats_matrix[["best_match_taxid", "description"]].reset_index(),
        left_on="best_taxid_match",
        right_on="best_match_taxid",
    )

    m_stats_stats_matrix_pruned = dataframe_update_with_lineage(
        pruned_overlap_manager.m_stats_matrix, ncbi_wrapper=ncbi_wrapper
    )
    m_stats_stats_matrix_pruned = _add_best_match_columns(m_stats_stats_matrix_pruned)

    pruned_results_pred_first = predict_data_set_clades_composition(
        sample,
        m_stats_stats_matrix=m_stats_stats_matrix_pruned,
        overlap_manager=pruned_overlap_manager,
        modeller=composition_modeller,
        input_taxa=taxids_to_use,
        tax_level=tax_level,
    )
    logger.info(
        f"[{sample}] Predictions: {len(pruned_results_pred_first)} clades, {len(m_stats_stats_matrix_pruned)} matrix rows"
    )

    pruned_results_pred = pruned_results_pred_first.merge(
        m_stats_stats_matrix_pruned[["best_match_taxid", "description"]].reset_index(),
        left_on="best_taxid_match",
        right_on="best_match_taxid",
    )

    if len(pruned_results_pred) == 0:
        logger.warning(f"No predictions left after pruning for sample {sample}")

    cluster_metrics = {
        "sample": sample,
        "n_total_intermediate_classifications": overlap_manager.original_m_stats_matrix.shape[0],
        "n_classifications_with_coverage": overlap_manager.m_stats_matrix[
            overlap_manager.m_stats_matrix["coverage"] > 0
        ].shape[0],
        "n_final_clusters": len(results_pred),
        "n_final_clusters_pruned": len(pruned_results_pred),
        "n_indexes_kept_from_recall": keep_index,
        "target_percentile": target_percentile,
        "target_recall": target_recall,
        "last_best_match_relindex": last_best_match_relindex,
        "resulting_percentile": resulting_percentile,
        "prop_coverage_above_cutoff": prop_coverage_above,
        "prop_coverage_below_cutoff": prop_coverage_below,
    }

    return results_pred_sample_df, pruned_results_pred, cluster_metrics


def generate_summary(sample_results: pd.DataFrame) -> pd.DataFrame:
    """
    Generate summary statistics for each sample.
    """
    summaries = []

    for sample in sample_results["data_set"].unique():
        sample_data = sample_results[sample_results["data_set"] == sample]

        if sample_data.empty:
            continue

        representative_row = sample_data.loc[sample_data["node_precision"].idxmax()]

        summaries.append(
            {
                "sample": sample,
                "n_detections": len(sample_data),
                "mean_precision": sample_data["node_precision"].mean(),
                "median_precision": sample_data["node_precision"].median(),
                "unique_taxa": sample_data["best_match_taxid"].nunique(),
                "n_clades": len(sample_data),
                "high_confidence_count": len(sample_data[sample_data["node_precision"] == 1.0]),
                "low_confidence_count": len(sample_data[sample_data["node_precision"] < 0.5]),
                "mean_n_leaves": sample_data["n_leaves"].mean(),
                "best_match_taxid": int(representative_row["best_match_taxid"]),
                "description": representative_row["description"],
                "n_leaves": int(representative_row["n_leaves"]),
            }
        )

    return pd.DataFrame(summaries)


def generate_cluster_stats(cluster_metrics: list[dict]) -> pd.DataFrame:
    """
    Generate cluster statistics for each sample.

    Args:
        cluster_metrics: List of cluster metric dictionaries from process_sample()

    Returns:
        DataFrame with cluster statistics per sample
    """
    if not cluster_metrics:
        return pd.DataFrame()

    return pd.DataFrame(cluster_metrics)


def generate_plots(sample_results: pd.DataFrame, output_dir: Path) -> None:
    """
    Generate visualization plots.
    """
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    sns.set_style("whitegrid")

    plt.figure(figsize=(10, 6))
    sns.histplot(
        sample_results["node_precision"],
        bins=20,
        kde=True,
        color="steelblue",
    )
    plt.xlabel("Precision")
    plt.ylabel("Count")
    plt.title("Distribution of Node Precision Across All Samples")
    plt.tight_layout()
    plt.savefig(plots_dir / "precision_distribution.png", dpi=150)
    plt.close()
    logger.info("Generated: precision_distribution.png")

    summary = generate_summary(sample_results)
    plt.figure(figsize=(12, 6))
    sns.barplot(data=summary, x="sample", y="n_detections", palette="viridis")
    plt.xlabel("Sample")
    plt.ylabel("Number of Detections")
    plt.title("Detection Counts by Sample")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(plots_dir / "sample_comparison.png", dpi=150)
    plt.close()
    logger.info("Generated: sample_comparison.png")

    top_taxa = sample_results.groupby("description")["data_set"].count().sort_values(ascending=False).head(20)
    plt.figure(figsize=(12, 8))
    sns.barplot(x=top_taxa.values, y=top_taxa.index, palette="magma")
    plt.xlabel("Frequency")
    plt.ylabel("Taxonomy")
    plt.title("Top 20 Most Frequently Detected Taxa")
    plt.tight_layout()
    plt.savefig(plots_dir / "taxa_frequency.png", dpi=150)
    plt.close()
    logger.info("Generated: taxa_frequency.png")

    clade_sizes = sample_results["n_leaves"].value_counts().sort_index()
    single_leaf = clade_sizes.get(1, 0)
    multi_leaf = clade_sizes[clade_sizes.index > 1].sum()

    plt.figure(figsize=(8, 8))
    plt.pie(
        [single_leaf, multi_leaf],
        labels=["Single-leaf clades", "Multi-leaf clades"],
        autopct="%1.1f%%",
        colors=["#2ecc71", "#e74c3c"],
        explode=(0.05, 0),
    )
    plt.title("Clade Size Distribution")
    plt.tight_layout()
    plt.savefig(plots_dir / "clade_size_distribution.png", dpi=150)
    plt.close()
    logger.info("Generated: clade_size_distribution.png")


def plot_filtering_benefit(cluster_stats_df: pd.DataFrame, output_dir: Path) -> None:
    """Plot coverage retention vs filtering for each sample."""
    if cluster_stats_df.empty:
        return

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))

    samples = cluster_stats_df["sample"].tolist()
    kept = cluster_stats_df["prop_coverage_above_cutoff"].tolist()
    lost = cluster_stats_df["prop_coverage_below_cutoff"].tolist()

    y_pos = range(len(samples))

    ax.barh(y_pos, kept, label="Coverage retained (kept)", color="steelblue")
    ax.barh(y_pos, lost, left=kept, label="Coverage lost (filtered)", color="coral")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(samples)
    ax.set_xlabel("Proportion of references with coverage")
    ax.set_title(
        "Filtering Benefit: Coverage Retained vs Lost\n(Higher kept = more time saved with minimal sensitivity loss)"
    )
    ax.legend(loc="lower right")
    ax.set_xlim(0, 1)

    plt.tight_layout()
    plt.savefig(plots_dir / "filtering_benefit.png", dpi=150)
    plt.close()
    logger.info("Generated: filtering_benefit.png")


def main():
    parser = argparse.ArgumentParser(description="Analyze samples using trained composition and recall models.")
    parser.add_argument(
        "--samples-dir",
        required=True,
        help="Directory containing sample subdirectories",
    )
    parser.add_argument(
        "--training-dir",
        required=True,
        help="Directory containing trained models",
    )
    parser.add_argument(
        "--output-dir",
        default="analysis_output",
        help="Output directory (default: analysis_output)",
    )
    parser.add_argument(
        "--tax-level",
        default="order",
        help="Taxonomic level for analysis (default: order)",
    )
    parser.add_argument(
        "--generate-plots",
        action="store_true",
        help="Generate visualization plots",
    )
    parser.add_argument(
        "--taxids-file",
        default=None,
        help="Path to taxids_to_use.tsv file (default: {training-dir}/taxids_to_use.tsv)",
    )
    parser.add_argument(
        "--taxonomy-db",
        default="taxa.db",
        help="Path to taxonomy database (default: taxa.db)",
    )
    parser.add_argument(
        "--target-recall",
        type=float,
        default=1.0,
        help="Target recall threshold (default: 1.0)",
    )
    parser.add_argument(
        "--data-set-divide",
        type=int,
        default=20,
        help="Number of recall divisions used in model training (default: 20)",
    )

    args = parser.parse_args()

    samples_dir = Path(args.samples_dir)
    training_dir = Path(args.training_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    model_dir = training_dir / "models"
    taxids_file = model_dir / "taxids_to_use.parquet"

    if not taxids_file.exists():
        raise FileNotFoundError(f"Taxids file not found: {taxids_file}")

    samples = [d for d in os.listdir(samples_dir) if os.path.isdir(samples_dir / d) and d.startswith("ERR")]

    if not samples:
        logger.error(f"No sample directories found in {samples_dir}")
        return

    logger.info(f"Found {len(samples)} samples to process")

    taxids_to_use = pd.read_parquet(taxids_file)

    ncbi_wrapper = NCBITaxonomistWrapper(db=args.taxonomy_db)

    sample_results = []
    sample_results_pruned = []
    cluster_metrics_list = []
    samples_dir.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        results_dir = samples_dir / sample / "output_clustering"

        if not results_dir.exists():
            logger.error(f"No output_clustering directory for sample {sample}")
            continue

        sdir_results = os.listdir(results_dir)
        if len(sdir_results) == 0:
            logger.error(f"No results for sample {sample}")
            continue

        results_dir = results_dir / sdir_results[0]

        results_pred, pruned_results_pred, cluster_metrics = process_sample(
            sample=sample,
            results_dir=results_dir,
            model_dir=model_dir,
            ncbi_wrapper=ncbi_wrapper,
            taxids_to_use=taxids_to_use,
            tax_level=args.tax_level,
            output_db=args.taxonomy_db,
            target_recall=args.target_recall,
            data_set_divide=args.data_set_divide,
        )

        if not results_pred.empty:
            sample_results.append(results_pred)
            cluster_metrics_list.append(cluster_metrics)

            sample_output_dir = output_dir / "samples" / sample
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            results_pred.to_csv(sample_output_dir / "predictions.tsv", sep="\t", index=False)

        if not pruned_results_pred.empty:
            sample_results_pruned.append(pruned_results_pred)
            sample_output_dir = output_dir / "samples" / sample
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            pruned_results_pred.to_csv(sample_output_dir / "predictions_pruned.tsv", sep="\t", index=False)

    if not sample_results:
        logger.error("No results to save")
        return

    final_results_df = pd.concat(sample_results, ignore_index=True)
    final_pruned_results_df = pd.concat(sample_results_pruned, ignore_index=True)

    final_results_df.to_csv(output_dir / "all_predictions.csv", index=False)
    logger.info(f"Saved: {output_dir / 'all_predictions.csv'}")

    final_pruned_results_df.to_csv(output_dir / "all_predictions_pruned.csv", index=False)
    logger.info(f"Saved: {output_dir / 'all_predictions_pruned.csv'}")

    summary_df = generate_summary(final_results_df)
    summary_df.to_csv(output_dir / "summary.tsv", sep="\t", index=False)
    logger.info(f"Saved: {output_dir / 'summary.tsv'}")

    cluster_stats_df = generate_cluster_stats(cluster_metrics_list)
    if not cluster_stats_df.empty:
        cluster_stats_df.to_csv(output_dir / "cluster_statistics.tsv", sep="\t", index=False)
        logger.info(f"Saved: {output_dir / 'cluster_statistics.tsv'}")

    metadata = {
        "generated_at": datetime.now().isoformat(),
        "samples_dir": str(samples_dir),
        "training_dir": str(training_dir),
        "tax_level": args.tax_level,
        "n_samples": len(samples),
        "n_successful": len(sample_results),
        "model_path": str(model_dir),
    }

    with open(output_dir / "metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    logger.info(f"Saved: {output_dir / 'metadata.json'}")

    if args.generate_plots:
        logger.info("Generating plots...")
        generate_plots(final_results_df, output_dir)
        plot_filtering_benefit(cluster_stats_df, output_dir)

    logger.info(f"Analysis complete! Results saved to: {output_dir}")
    logger.info(f"  - {output_dir / 'summary.tsv'}")
    logger.info(f"  - {output_dir / 'all_predictions.csv'}")
    logger.info(f"  - {output_dir / 'cluster_statistics.tsv'}")
    logger.info(f"  - {output_dir / 'metadata.json'}")
    logger.info(f"  - {output_dir / 'samples/'}")


if __name__ == "__main__":
    main()
