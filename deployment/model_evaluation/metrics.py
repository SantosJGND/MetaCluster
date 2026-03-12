"""
Pure metric calculation functions for evaluator module.

These functions have no side effects and operate only on their inputs.
"""

from typing import Tuple, Optional
import pandas as pd
import numpy as np


def compute_fuzzy_precision(
    m_stats_matrix: pd.DataFrame,
    is_trash_column: str = 'is_trash'
) -> Tuple[float, float]:
    """
    Compute fuzzy precision at raw and coverage-filtered levels.
    
    Args:
        m_stats_matrix: DataFrame with leaf statistics
        is_trash_column: Name of column indicating if leaf is trash
        
    Returns:
        Tuple of (fuzzy_precision_raw, fuzzy_precision_cov_filtered)
    """
    total = len(m_stats_matrix)
    if total == 0:
        return 0.0, 0.0
    
    clean = m_stats_matrix[m_stats_matrix[is_trash_column] == False]
    raw = len(clean) / total
    
    cov_filtered = len(clean[clean['coverage'] > 0]) / total
    return raw, cov_filtered


def compute_overall_precision(m_stats_matrix: pd.DataFrame) -> float:
    """
    Compute overall precision as unique taxids / total leaves.
    
    Args:
        m_stats_matrix: DataFrame with leaf statistics
        
    Returns:
        Overall precision value
    """
    if len(m_stats_matrix) == 0:
        return 0.0
    unique_matches = m_stats_matrix.dropna(subset=['best_match_taxid'])
    return unique_matches['best_match_taxid'].nunique() / len(m_stats_matrix)


def compute_recall(
    clean_m_stats: pd.DataFrame,
    input_summary: pd.DataFrame
) -> Tuple[float, float, float, float]:
    """
    Compute all recall metrics.
    
    Args:
        clean_m_stats: Filtered m_stats matrix (non-trash)
        input_summary: Input data summary with taxid information
        
    Returns:
        Tuple of (recall_raw, recall_cov_filtered, clade_recall, recall_filtered_leaves)
    """
    unique_taxids = input_summary['taxid'].nunique()
    if unique_taxids == 0:
        return 0.0, 0.0, 0.0, 0.0
    
    clean_unique = clean_m_stats.drop_duplicates(subset=['best_match_taxid']
    ).dropna(subset=['best_match_taxid'])
    
    recall_raw = len(clean_unique) / unique_taxids
    
    cov_filtered = clean_m_stats[clean_m_stats['coverage'] > 0]
    recall_cov = len(cov_filtered.drop_duplicates(subset=['best_match_taxid'])
    ) / unique_taxids
    
    return recall_raw, recall_cov, 0.0, 0.0


def compute_trash_flags(m_stats_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Add is_trash column to m_stats_matrix.
    
    A leaf is considered trash if:
    - best_match_is_best is False AND
    - is_crosshit is False
    
    Args:
        m_stats_matrix: DataFrame with leaf statistics
        
    Returns:
        Copy of m_stats_matrix with is_trash column added
    """
    m_stats = m_stats_matrix.copy()
    m_stats['is_trash'] = (m_stats['best_match_is_best'] == False) & (m_stats['is_crosshit'] == False)
    return m_stats


def compute_raw_pred_accuracy(
    m_stats_matrix: pd.DataFrame,
    input_summary: pd.DataFrame
) -> pd.Series:
    """
    Compute raw prediction accuracy for each taxid.
    
    Args:
        m_stats_matrix: DataFrame with leaf statistics
        input_summary: Input data summary with taxid information
        
    Returns:
        Series with accuracy values indexed by taxid
    """
    accuracy = input_summary['taxid'].apply(
        lambda x: (m_stats_matrix['best_match_taxid'] == x).sum()
    )
    return accuracy


def compute_clade_accuracy(
    results_df: pd.DataFrame,
    input_summary: pd.DataFrame
) -> int:
    """
    Compute clade prediction accuracy.
    
    Args:
        results_df: Predicted clades DataFrame
        input_summary: Input data summary with taxid information
        
    Returns:
        Total number of correct predictions
    """
    return input_summary['taxid'].apply(
        lambda x: (results_df['best_taxid_match'] == x).sum()
    ).sum()


def compute_clade_recall(
    results_df: pd.DataFrame,
    input_summary: pd.DataFrame
) -> float:
    """
    Compute clade recall.
    
    Args:
        results_df: Predicted clades DataFrame
        input_summary: Input data summary with taxid information
        
    Returns:
        Recall value
    """
    unique_taxids = input_summary['taxid'].nunique()
    if unique_taxids == 0:
        return 0.0
    return results_df['best_taxid_match'].dropna().nunique() / unique_taxids


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


def compute_precision_stats(m_stats_matrix: pd.DataFrame) -> dict:
    """
    Compute all precision-related statistics.
    
    Args:
        m_stats_matrix: DataFrame with leaf statistics
        
    Returns:
        Dictionary with precision statistics
    """
    total = len(m_stats_matrix)
    if total == 0:
        return {
            'fuzzy_precision_raw': 0.0,
            'fuzzy_precision_cov_filtered': 0.0,
            'overall_precision_raw': 0.0,
        }
    
    m_stats = compute_trash_flags(m_stats_matrix)
    clean = m_stats[m_stats['is_trash'] == False]
    
    fuzzy_raw = safe_divide(len(clean), total)
    fuzzy_cov = safe_divide(len(clean[clean['coverage'] > 0]), total)
    overall = compute_overall_precision(m_stats)
    
    return {
        'fuzzy_precision_raw': fuzzy_raw,
        'fuzzy_precision_cov_filtered': fuzzy_cov,
        'overall_precision_raw': overall,
    }


def compute_recall_stats(
    clean_m_stats: pd.DataFrame,
    input_summary: pd.DataFrame
) -> dict:
    """
    Compute all recall-related statistics.
    
    Args:
        clean_m_stats: Filtered m_stats matrix (non-trash)
        input_summary: Input data summary with taxid information
        
    Returns:
        Dictionary with recall statistics
    """
    unique_taxids = input_summary['taxid'].nunique()
    if unique_taxids == 0:
        return {
            'recall_raw': 0.0,
            'recall_cov_filtered': 0.0,
        }
    
    clean_unique = clean_m_stats.drop_duplicates(subset=['best_match_taxid']
    ).dropna(subset=['best_match_taxid'])
    
    recall_raw = safe_divide(len(clean_unique), unique_taxids)
    
    cov_filtered = clean_m_stats[clean_m_stats['coverage'] > 0]
    recall_cov = safe_divide(
        len(cov_filtered.drop_duplicates(subset=['best_match_taxid'])),
        unique_taxids
    )
    
    return {
        'recall_raw': recall_raw,
        'recall_cov_filtered': recall_cov,
    }


def fill_missing_tax_levels(
    composition_df: pd.DataFrame,
    missing_levels: set,
    fill_value: float = np.nan
) -> pd.DataFrame:
    """
    Fill missing tax levels in composition DataFrame.
    
    Args:
        composition_df: Composition DataFrame
        missing_levels: Set of missing tax levels
        fill_value: Value to use for missing levels
        
    Returns:
        Updated composition DataFrame
    """
    if not missing_levels:
        return composition_df
    
    result = composition_df.copy()
    for level in missing_levels:
        result[level] = fill_value
    return result
