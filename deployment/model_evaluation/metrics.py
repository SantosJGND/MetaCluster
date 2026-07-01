"""
Pure metric calculation functions for evaluator module.

These functions have no side effects and operate only on their inputs.
"""

import logging

import pandas as pd

logger = logging.getLogger(__name__)


def compute_cross_hit_stats(results: list["DatasetResult"]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute cross-hit statistics from DatasetResult objects.

    For each class, computes:
    - Reads simulated for class
    - Cross-hits found (count)
    - Total cross-hit reads mapped
    - Ratio of cross-hit reads to simulated reads

    Args:
        results: List of DatasetResult objects with cross-hit metrics populated

    Returns:
        Tuple of (per_dataset_df, class_stats):
        - per_dataset_df: DataFrame with per-dataset metrics (class, data_set, reads_simulated,
                        cross_hit_count, cross_hit_reads_mapped, ratio)
        - class_stats: Aggregated statistics by class (mean, median, std, min, max, count)
    """
    records = []
    for result in results:
        if result.reads_simulated_per_class is None:
            continue

        cross_hit_counts = {}
        if result.cross_hit.cross_hit_counts_per_class:
            for item in result.cross_hit.cross_hit_counts_per_class:
                cross_hit_counts[item["class"]] = item["count"]

        cross_hit_reads = {}
        if result.cross_hit.cross_hit_reads_per_class:
            for item in result.cross_hit.cross_hit_reads_per_class:
                cross_hit_reads[item["class"]] = item["reads"]

        reads_simulated_lookup = {}
        for item in result.reads_simulated_per_class:
            reads_simulated_lookup[item["class"]] = item["reads"]

        for item in result.reads_simulated_per_class:
            cls = item["class"]
            reads_sim = reads_simulated_lookup.get(cls, 0)
            cross_hit_mapped = cross_hit_reads.get(cls, 0)
            ratio = cross_hit_mapped / reads_sim if reads_sim > 0 else 0.0
            records.append(
                {
                    "class": cls,
                    "data_set": result.data_set,
                    "reads_simulated": reads_sim,
                    "cross_hit_count": cross_hit_counts.get(cls, 0),
                    "cross_hit_reads_mapped": cross_hit_mapped,
                    "ratio": ratio,
                }
            )

    per_dataset_df = pd.DataFrame(records)
    if per_dataset_df.empty:
        return per_dataset_df, pd.DataFrame()

    class_stats = (
        per_dataset_df.groupby("class")
        .agg(
            {
                "reads_simulated": ["mean", "median", "std", "min", "max", "count"],
                "cross_hit_count": ["mean", "median", "std", "min", "max"],
                "cross_hit_reads_mapped": ["mean", "median", "std", "min", "max"],
                "ratio": ["mean", "median", "std", "min", "max"],
            }
        )
        .reset_index()
    )

    class_stats.columns = [
        "class",
        "reads_simulated_mean",
        "reads_simulated_median",
        "reads_simulated_std",
        "reads_simulated_min",
        "reads_simulated_max",
        "n_simulations",
        "cross_hit_count_mean",
        "cross_hit_count_median",
        "cross_hit_count_std",
        "cross_hit_count_min",
        "cross_hit_count_max",
        "cross_hit_reads_mapped_mean",
        "cross_hit_reads_mapped_median",
        "cross_hit_reads_mapped_std",
        "cross_hit_reads_mapped_min",
        "cross_hit_reads_mapped_max",
        "ratio_mean",
        "ratio_median",
        "ratio_std",
        "ratio_min",
        "ratio_max",
    ]

    return per_dataset_df, class_stats


def compute_spurious_hit_stats(results: list["DatasetResult"]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute spurious-hit statistics from DatasetResult objects.

    For each class, computes:
    - Reads simulated for class
    - Spurious-hits found (count)
    - Total spurious-hit reads mapped
    - Ratio of spurious-hit reads to simulated reads

    Args:
        results: List of DatasetResult objects with spurious-hit metrics populated
    Returns:
        Tuple of (per_dataset_df, class_stats):
        - per_dataset_df: DataFrame with per-dataset metrics (class, data_set, reads_simulated, spurious_hit_count, spurious_hit_reads_mapped, ratio)
        - class_stats: Aggregated statistics by class (mean, median, std, min, max, count)
    """
    records = []
    for result in results:
        if result.reads_simulated_per_class is None:
            continue

        spurious_hit_counts = {}
        if result.cross_hit.spurious_hit_counts_per_class:
            for item in result.cross_hit.spurious_hit_counts_per_class:
                spurious_hit_counts[item["class"]] = item["count"]

        spurious_hit_reads = {}
        if result.cross_hit.spurious_hit_reads_per_class:
            for item in result.cross_hit.spurious_hit_reads_per_class:
                spurious_hit_reads[item["class"]] = item["reads"]

        reads_simulated_lookup = {}
        for item in result.reads_simulated_per_class:
            reads_simulated_lookup[item["class"]] = item["reads"]

        for item in result.reads_simulated_per_class:
            cls = item["class"]
            reads_sim = reads_simulated_lookup.get(cls, 0)
            spurious_hit_mapped = spurious_hit_reads.get(cls, 0)
            ratio = spurious_hit_mapped / reads_sim if reads_sim > 0 else 0.0
            records.append(
                {
                    "class": cls,
                    "data_set": result.data_set,
                    "reads_simulated": reads_sim,
                    "spurious_hit_count": spurious_hit_counts.get(cls, 0),
                    "spurious_hit_reads_mapped": spurious_hit_mapped,
                    "ratio": ratio,
                }
            )

    per_dataset_df = pd.DataFrame(records)
    if per_dataset_df.empty:
        return per_dataset_df, pd.DataFrame()

    class_stats = (
        per_dataset_df.groupby("class")
        .agg(
            {
                "reads_simulated": ["mean", "median", "std", "min", "max", "count"],
                "spurious_hit_count": ["mean", "median", "std", "min", "max"],
                "spurious_hit_reads_mapped": ["mean", "median", "std", "min", "max"],
                "ratio": ["mean", "median", "std", "min", "max"],
            }
        )
        .reset_index()
    )

    class_stats.columns = [
        "class",
        "reads_simulated_mean",
        "reads_simulated_median",
        "reads_simulated_std",
        "reads_simulated_min",
        "reads_simulated_max",
        "n_simulations",
        "spurious_hit_count_mean",
        "spurious_hit_count_median",
        "spurious_hit_count_std",
        "spurious_hit_count_min",
        "spurious_hit_count_max",
        "spurious_hit_reads_mapped_mean",
        "spurious_hit_reads_mapped_median",
        "spurious_hit_reads_mapped_std",
        "spurious_hit_reads_mapped_min",
        "spurious_hit_reads_mapped_max",
        "ratio_mean",
        "ratio_median",
        "ratio_std",
        "ratio_min",
        "ratio_max",
    ]

    return per_dataset_df, class_stats


def compute_purity(m_stats_matrix: pd.DataFrame, is_trash_column: str = "is_trash") -> tuple[float, float, float]:
    """
    Compute matching purity at raw and coverage-filtered levels.

    Measures the proportion of non-trash leaves (leaves that are either
    a best-match or a cross-hit) relative to total leaves.

    Args:
        m_stats_matrix: DataFrame with leaf statistics
        is_trash_column: Name of column indicating if leaf is spurious

    Returns:
        Tuple of (purity_raw, purity_cov_filtered)
    """
    total = len(m_stats_matrix)
    if total == 0:
        return 0.0, 0.0, 0.0

    raw = len(m_stats_matrix[m_stats_matrix[is_trash_column] == False]) / total
    raw_best_match = (
        len(m_stats_matrix[(m_stats_matrix[is_trash_column] == False) & (m_stats_matrix["best_match_is_best"] == True)])
        / total
    )
    cov_filtered = (
        len(m_stats_matrix[(m_stats_matrix[is_trash_column] == False) & (m_stats_matrix["coverage"] > 0)]) / total
    )

    return raw, raw_best_match, cov_filtered


def compute_mstats_precision(m_stats_matrix: pd.DataFrame, input_summary: pd.DataFrame) -> float:
    """
    Compute overall precision as |predicted_taxids ∩ input_taxids| / |predicted_taxids|.

    Measures what fraction of predicted (best-match) taxids actually
    appear in the simulated input.

    Args:
        m_stats_matrix: DataFrame with leaf statistics
        input_summary: Input data summary with taxid information

    Returns:
        Overall precision value
    """
    if len(m_stats_matrix) == 0:
        return 0.0
    unique_matches = m_stats_matrix.dropna(subset=["best_match_taxid"])
    unique_matches = unique_matches[
        (unique_matches["is_trash"] == False) & (unique_matches["best_match_is_best"] == True)
    ]
    if len(unique_matches) == 0:
        return 0.0
    input_taxids = set(input_summary["taxid"].unique())
    output_taxids = set(unique_matches["best_match_taxid"].unique())
    if len(output_taxids) == 0:
        return 0.0
    correct = len(output_taxids & input_taxids)
    return correct / len(output_taxids)


def compute_recall(clean_m_stats: pd.DataFrame, input_summary: pd.DataFrame) -> tuple[float, float, list, list]:
    """
    Compute all recall metrics.

    Args:
        clean_m_stats: Filtered m_stats matrix (non-spurious)
        input_summary: Input data summary with taxid information

    Returns:
        Tuple of (recall_raw, recall_cov_filtered, clade_recall, recall_filtered_leaves)
    """
    input_taxids = set(input_summary["taxid"].unique())
    unique_taxids = len(input_taxids)

    if unique_taxids == 0:
        return 0.0, 0.0, [], []
    clean_unique = clean_m_stats.dropna(subset=["best_match_taxid"])
    clean_unique = clean_unique[(clean_unique["is_trash"] == False) & (clean_unique["best_match_is_best"] == True)]
    output_taxids = set(clean_unique["best_match_taxid"].dropna().unique())
    recall_raw = len(output_taxids & input_taxids) / unique_taxids
    clean_cov_filtered = clean_unique[(clean_unique["coverage"] > 0)]
    output_taxids_cov = set(clean_cov_filtered["best_match_taxid"].dropna().astype(int).unique())
    recall_cov = len(output_taxids_cov & input_taxids) / unique_taxids
    logger.debug(
        f"Unique taxids in input: {unique_taxids}, Unique matches in clean m_stats: {len(clean_unique)}, Recall raw: {recall_raw}, Recall coverage filtered: {recall_cov}"
    )
    return recall_raw, recall_cov, sorted(output_taxids & input_taxids), sorted(output_taxids_cov & input_taxids)


def compute_clade_recall(results_df: pd.DataFrame, input_summary: pd.DataFrame) -> float:
    """
    Compute clade recall as |predicted_taxids ∩ input_taxids| / |input_taxids|.

    Measures what fraction of input taxids are recovered in the
    predicted clades.

    Args:
        results_df: Predicted clades DataFrame
        input_summary: Input data summary with taxid information

    Returns:
        Recall value
    """
    input_taxids = set(input_summary["taxid"].unique())
    unique_taxids = len(input_taxids)
    if unique_taxids == 0:
        return 0.0
    predicted_taxids = set(results_df["best_taxid_match"].dropna().astype(int).unique())
    return len(predicted_taxids & input_taxids) / unique_taxids


def safe_divide(numerator: float, denominator: float) -> float:
    """
    Safe division that handles zero denominator.

    Args:
        numerator: Numerator value
        denominator: Denominator value

    Returns:
        Result of division or 0.0 if denominator is 0
    """
    if denominator == 0:
        return 0.0
    return numerator / denominator
