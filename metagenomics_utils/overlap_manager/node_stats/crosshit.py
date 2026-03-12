"""
Cross-hit detection functions.
"""
import pandas as pd

from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper


def id_crosshit(row, overlap_manager, m_stats_stats_matrix: pd.DataFrame, ncbi_wrapper: NCBITaxonomistWrapper, best_hits: pd.DataFrame, cross_hit_threshold: float = 0.3):
    """Identify cross-hit matches."""
    row['is_crosshit'] = False
    row['cross_hit_match'] = None
    row['crosshit_match_score'] = 0.0
    row['crosshit_match_level'] = None
    row['crossh_hit_dist'] = 1.0

    if row['best_match_is_best'] is True:
        return row

    if row['leaf'] not in overlap_manager.distance_mat.index:
        return row

    best_match = overlap_manager.distance_mat.loc[row['leaf']]

    if best_match.empty or len(best_hits) == 0:
        return row

    if best_match[best_match.index.isin(best_hits['leaf'])].empty:
        return row

    best_match_dist = best_match[best_match.index.isin(best_hits['leaf'])].min()
    best_match = best_match[best_match.index.isin(best_hits['leaf'])].idxmin()

    best_match_taxid = m_stats_stats_matrix[m_stats_stats_matrix['leaf'] == best_match]['best_match_taxid'].values
    if len(best_match_taxid) > 0:
        best_match_taxid = best_match_taxid[0]
        score, level = ncbi_wrapper.compare_lineages_relative(best_match_taxid, row['taxid'])
        if best_match_dist < cross_hit_threshold:
            row['is_crosshit'] = True
        row['cross_hit_match'] = best_match_taxid
        row['crosshit_match_score'] = score
        row['crosshit_match_level'] = level
        row['crossh_hit_dist'] = best_match_dist

    return row


def get_max_shared(row, overlap_manager, best_hits: pd.DataFrame):
    """Calculate maximum shared reads."""
    if row['leaf'] not in overlap_manager.distance_mat.index:
        return row

    best_match = overlap_manager.distance_mat.loc[row['leaf']]
    if best_match.empty or len(best_hits) == 0:
        return row
    
    if best_match[best_match.index != row['leaf']].empty:
        return row

    best_match_dist = best_match[best_match.index != row['leaf']].min()

    row['max_shared'] = 1.0 - best_match_dist
    return row
