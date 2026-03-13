#!/usr/bin/env python3
"""
Sample analysis script using trained models to analyze new metagenomic data.

Analyzes clustering output from map_clustering pipeline and generates predictions
using trained RecallModeller and CompositionModeller.
"""

import argparse
import json
import os
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper
from metagenomics_utils.overlap_manager import OverlapManager
from metagenomics_utils.overlap_manager.diversity import (
    kurtosis,
    shannon_diversity_from_list,
    skewness,
)
from metagenomics_utils.overlap_manager.node_stats import (
    dataframe_update_with_lineage,
    get_subset_composition_counts,
)
from metagenomics_utils.overlap_manager.om_models import (
    CompositionModeller,
    RecallModeller,
    predict_data_set_clades,
)


def predict_recall_cutoff_vars(
    data_set_name: str,
    m_stats_stats_matrix: pd.DataFrame,
    input_tax_df: pd.DataFrame,
    tax_level: str = "order",
) -> pd.DataFrame:
    """
    Predict recall at various cutoffs and other composition statistics.
    """
    m_stats_stats_matrix = m_stats_stats_matrix.sort_values(
        by="total_uniq_reads", ascending=False
    ).reset_index(drop=True)

    composition = get_subset_composition_counts(
        m_stats_stats_matrix, input_tax_df, tax_level=tax_level
    )[["tax_level", "proportion"]].set_index("tax_level").T

    tax_diversity_shannon = shannon_diversity_from_list(
        m_stats_stats_matrix[tax_level].dropna().tolist()
    )

    counts_skewness = skewness(m_stats_stats_matrix["total_uniq_reads"].tolist())
    counts_kurtosis = kurtosis(m_stats_stats_matrix["total_uniq_reads"].tolist())

    composition.insert(0, "tax_diversity_shannon", tax_diversity_shannon)
    composition.insert(0, "counts_skewness", counts_skewness)
    composition.insert(0, "counts_kurtosis", counts_kurtosis)

    composition.insert(0, "data_set", data_set_name)

    composition.reset_index(drop=True, inplace=True)

    return composition


def process_sample(
    sample: str,
    results_dir: Path,
    model_dir: Path,
    ncbi_wrapper: NCBITaxonomistWrapper,
    taxids_to_use: pd.DataFrame,
    tax_level: str,
    output_db: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Process a single sample and return predictions and cluster metrics.
    
    Returns:
        Tuple of (predictions DataFrame, pruned predictions DataFrame, cluster metrics dict)
    """
    merged_classification = results_dir / "classification" / f"{sample}_merged_classification.tsv"
    matched_assemblies = results_dir / "output" / "matched_assemblies.tsv"

    if not merged_classification.exists() or not matched_assemblies.exists():
        print(f"Missing files for sample {sample}")
        return pd.DataFrame(), pd.DataFrame(), {}

    merged_classification_df = pd.read_csv(merged_classification, sep="\t")
    matched_assemblies_df = pd.read_csv(matched_assemblies, sep="\t")

    if not matched_assemblies_df["taxid"].isin(merged_classification_df["taxid"]).all():
        print(
            f"Not all taxids in matched assemblies are present in merged classification for sample {sample}"
        )
        return pd.DataFrame(), pd.DataFrame(), {}

    ncbi_wrapper.resolve_lineages(taxids_to_use["taxid"].tolist())

    mock_composition_stats = pd.DataFrame(
        columns=[
            "data_set",
            "node",
            "n_true_leaves",
            "precision_increased",
            "new_precision",
            "precision",
            "stop_traversal",
            "unclassified",
        ]
    )
    composition_modeller = CompositionModeller(mock_composition_stats)
    composition_modeller.load_model(os.path.join(model_dir))

    mock_recall_stats = pd.DataFrame(columns=["data_set"])
    recall_modeller = RecallModeller(mock_recall_stats, 5)
    recall_modeller.load_model(os.path.join(model_dir))

    overlap_manager = OverlapManager(
        os.path.join(results_dir, "clustering"), max_taxids=1000
    )
    if not overlap_manager.check_data_available():
        print(f"No overlapping data available for sample {sample}")
        return pd.DataFrame(), {}

    result_taxids = overlap_manager.m_stats_matrix["taxid"].unique().tolist()
    ncbi_wrapper.resolve_lineages(result_taxids)

    m_stats_stats_matrix = dataframe_update_with_lineage(
        overlap_manager.m_stats_matrix, ncbi_wrapper=ncbi_wrapper
    )
    m_stats_stats_matrix["best_match_score"] = 1
    m_stats_stats_matrix["best_match_taxid"] = m_stats_stats_matrix["taxid"]
    m_stats_stats_matrix["best_match_is_best"] = True

    recall_stats = predict_recall_cutoff_vars(
        sample, m_stats_stats_matrix, taxids_to_use, tax_level=tax_level
    )

    recall_pred = recall_modeller.model.predict(
        recall_stats[recall_modeller.RecP_feature_cols]
    )
    recall_pred_df = pd.DataFrame(recall_pred, columns=recall_modeller.RecP_target_cols)
    keep_index = (
        recall_pred_df.iloc[0]["last_best_match_relindex"]
        * overlap_manager.m_stats_matrix.shape[0]
    )
    keep_index = round(keep_index + 1)
    print(f"{sample}: keeping top {keep_index} predictions")

    pruned_overlap_manager = OverlapManager(
        os.path.join(results_dir, "clustering"), max_taxids=keep_index
    )

    results_pred = predict_data_set_clades(
        sample,
        m_stats_stats_matrix=m_stats_stats_matrix,
        overlap_manager=overlap_manager,
        modeller=composition_modeller,
        input_taxa=taxids_to_use,
        tax_level=tax_level,
    )

    pruned_results_pred = predict_data_set_clades(
        sample,
        m_stats_stats_matrix=m_stats_stats_matrix,
        overlap_manager=pruned_overlap_manager,
        modeller=composition_modeller,
        input_taxa=taxids_to_use,
        tax_level=tax_level,
    )

    results_pred_sample_df = results_pred.merge(
        m_stats_stats_matrix[["best_match_taxid", "description"]].reset_index(),
        left_on="best_taxid_match",
        right_on="best_match_taxid",
    )

    pruned_results_pred = pruned_results_pred.merge(
        m_stats_stats_matrix[["best_match_taxid", "description"]].reset_index(),
        left_on="best_taxid_match",
        right_on="best_match_taxid",
    )

    cluster_metrics = {
        "sample": sample,
        "n_total_intermediate_classifications": overlap_manager.m_stats_matrix.shape[0],
        "n_classifications_with_coverage": overlap_manager.m_stats_matrix[
            overlap_manager.m_stats_matrix["coverage"] > 0
        ].shape[0],
        "n_final_clusters": len(results_pred),
        "n_final_clusters_pruned": len(pruned_results_pred),
        "n_indexes_kept_from_recall": keep_index,
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
                "high_confidence_count": len(
                    sample_data[sample_data["node_precision"] == 1.0]
                ),
                "low_confidence_count": len(
                    sample_data[sample_data["node_precision"] < 0.5]
                ),
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
    print("Generated: precision_distribution.png")

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
    print("Generated: sample_comparison.png")

    top_taxa = (
        sample_results.groupby("description")["data_set"]
        .count()
        .sort_values(ascending=False)
        .head(20)
    )
    plt.figure(figsize=(12, 8))
    sns.barplot(x=top_taxa.values, y=top_taxa.index, palette="magma")
    plt.xlabel("Frequency")
    plt.ylabel("Taxonomy")
    plt.title("Top 20 Most Frequently Detected Taxa")
    plt.tight_layout()
    plt.savefig(plots_dir / "taxa_frequency.png", dpi=150)
    plt.close()
    print("Generated: taxa_frequency.png")

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
    print("Generated: clade_size_distribution.png")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze samples using trained composition and recall models."
    )
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

    args = parser.parse_args()

    samples_dir = Path(args.samples_dir)
    training_dir = Path(args.training_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    taxids_file = args.taxids_file or training_dir / "taxids_to_use.tsv"
    model_dir = training_dir / "models"

    if not taxids_file.exists():
        raise FileNotFoundError(f"Taxids file not found: {taxids_file}")

    samples = [
        d
        for d in os.listdir(samples_dir)
        if os.path.isdir(samples_dir / d) and d.startswith("ERR")
    ]

    if not samples:
        print(f"No sample directories found in {samples_dir}")
        return

    print(f"Found {len(samples)} samples to process")

    taxids_to_use = pd.read_csv(taxids_file, sep="\t")

    ncbi_wrapper = NCBITaxonomistWrapper(db=args.taxonomy_db)

    sample_results = []
    sample_results_pruned = []
    cluster_metrics_list = []
    samples_dir.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        results_dir = samples_dir / sample / "output_clustering"

        if not results_dir.exists():
            print(f"No output_clustering directory for sample {sample}")
            continue

        sdir_results = os.listdir(results_dir)
        if len(sdir_results) == 0:
            print(f"No results for sample {sample}")
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
        )

        if not results_pred.empty:
            sample_results.append(results_pred)
            cluster_metrics_list.append(cluster_metrics)

            sample_output_dir = output_dir / "samples" / sample
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            results_pred.to_csv(
                sample_output_dir / "predictions.tsv", sep="\t", index=False
            )

        if not pruned_results_pred.empty:
            sample_results_pruned.append(pruned_results_pred)
            sample_output_dir = output_dir / "samples" / sample
            sample_output_dir.mkdir(parents=True, exist_ok=True)
            pruned_results_pred.to_csv(
                sample_output_dir / "predictions_pruned.tsv", sep="\t", index=False
            )

    if not sample_results:
        print("No results to save")
        return

    final_results_df = pd.concat(sample_results, ignore_index=True)

    final_results_df.to_csv(output_dir / "all_predictions.csv", index=False)
    print(f"Saved: {output_dir / 'all_predictions.csv'}")

    summary_df = generate_summary(final_results_df)
    summary_df.to_csv(output_dir / "summary.tsv", sep="\t", index=False)
    print(f"Saved: {output_dir / 'summary.tsv'}")

    cluster_stats_df = generate_cluster_stats(cluster_metrics_list)
    if not cluster_stats_df.empty:
        cluster_stats_df.to_csv(output_dir / "cluster_statistics.tsv", sep="\t", index=False)
        print(f"Saved: {output_dir / 'cluster_statistics.tsv'}")

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
    print(f"Saved: {output_dir / 'metadata.json'}")

    if args.generate_plots:
        print("Generating plots...")
        generate_plots(final_results_df, output_dir)

    print(f"\nAnalysis complete! Results saved to: {output_dir}")
    print(f"  - {output_dir / 'summary.tsv'}")
    print(f"  - {output_dir / 'all_predictions.csv'}")
    print(f"  - {output_dir / 'cluster_statistics.tsv'}")
    print(f"  - {output_dir / 'metadata.json'}")
    print(f"  - {output_dir / 'samples/'}")


if __name__ == "__main__":
    main()
