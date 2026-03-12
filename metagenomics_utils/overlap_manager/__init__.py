"""
Overlap manager module for clustering analysis.
"""

from metagenomics_utils.overlap_manager.manager import OverlapManager, merge_by_assembly_ID, merge_to_matched
from metagenomics_utils.overlap_manager.node_stats import (
    get_m_stats_matrix,
    node_composition_level,
    node_leaf_shannon_tax_diversity,
    node_total_true_leaves,
    node_leaves_best_taxids,
    get_subset_composition,
    get_composition_by_leaf,
    get_subset_composition_counts,
    normalize_by_taxlevel,
    compute_node_stats,
    compute_node_purity,
    update_df_best_match,
    match_leaf,
)
from metagenomics_utils.overlap_manager.diversity import (
    shannon_diversity,
    shannon_diversity_from_counts,
    shannon_diversity_from_list,
    skewness,
    kurtosis,
)
