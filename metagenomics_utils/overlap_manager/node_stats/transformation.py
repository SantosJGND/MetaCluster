"""
Data transformation functions for node statistics.
"""
import pandas as pd


def normalize_by_taxlevel(prediction_matrix: pd.DataFrame, tax_level: str = 'genus'):
    """Normalize statistics by taxonomic level."""
    normalized_matrix = []
    
    for _, group in prediction_matrix.groupby(tax_level):
        stats_cols = ['coverage', 'covbases', 'meanmapq', 'error_rate']
        group_stats = group[stats_cols]
        group_scaled = group.copy()
        group_nocols = group.drop(columns=stats_cols)
        group_scaled = pd.concat([group_nocols, group_stats], axis=1)

        normalized_matrix.append(group_scaled)
    if len(normalized_matrix) == 0:
        return pd.DataFrame()
    normalized_matrix = pd.concat(normalized_matrix, ignore_index=True)

    return normalized_matrix
