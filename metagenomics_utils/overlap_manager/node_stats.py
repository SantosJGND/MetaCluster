import pandas as pd
import networkx as nx
import numpy as np
from collections import Counter
import os

from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper
from metagenomics_utils.overlap_manager.diversity import (
    shannon_diversity,
    shannon_diversity_from_counts,
    shannon_diversity_from_list,
    skewness,
    kurtosis,
)

from typing import List

################################################# node stats and composition functions #################################################



def normalize_by_taxlevel(prediction_matrix:pd.DataFrame, tax_level: str = 'genus'):
    from sklearn.preprocessing import StandardScaler
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


def id_crosshit(row, overlap_manager, m_stats_stats_matrix: pd.DataFrame, ncbi_wrapper: NCBITaxonomistWrapper, best_hits: pd.DataFrame, cross_hit_threshold: float = 0.3):
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


def dataframe_update_with_lineage(df: pd.DataFrame, ncbi_wrapper: NCBITaxonomistWrapper) -> pd.DataFrame:
    df['order'] = df.apply(lambda row: ncbi_wrapper.get_level(row['taxid'], 'order'), axis=1)
    df['order'] = df['order'].fillna('unclassified')
    df['family'] = df.apply(lambda row: ncbi_wrapper.get_level(row['taxid'], 'family'), axis=1)
    df['family'] = df['family'].fillna('unclassified')
    df['genus'] = df.apply(lambda row: ncbi_wrapper.get_level(row['taxid'], 'genus'), axis=1)
    df['genus'] = df['genus'].fillna('unclassified')
    return df


def get_m_stats_matrix(data_set_name, study_output_filepath, ncbi_wrapper: NCBITaxonomistWrapper, overlap_manager, cross_hit_threshold: float = 0.3):
    from metagenomics_utils.overlap_manager.manager import merge_to_matched, merge_by_assembly_ID
    
    matched_assemblies_file = os.path.join(study_output_filepath, f"{data_set_name}", "output", "matched_assemblies.tsv")
    merged_stats_file = os.path.join(study_output_filepath, f"{data_set_name}", "output", "merged_coverage_statistics.tsv")
    output_filepath = os.path.join(study_output_filepath, f"{data_set_name}", "input", f"{data_set_name}.tsv")

    if not os.path.exists(matched_assemblies_file) or not os.path.exists(merged_stats_file) or not os.path.exists(output_filepath):
        return pd.DataFrame()

    input_df = pd.read_csv(output_filepath, sep="\t")
    input_taxids = input_df['taxid'].dropna().unique().tolist()
    matched = pd.read_csv(matched_assemblies_file, sep="\t")
    m_stats_stats_matrix = pd.read_csv(merged_stats_file, sep="\t").rename(columns={"#rname": "assembly_accession"})

    if "error_rate" in m_stats_stats_matrix.columns:
        m_stats_stats_matrix = m_stats_stats_matrix[['assembly_accession', "coverage", "covbases", "meanmapq", "error_rate", "file"]]
    else:
        m_stats_stats_matrix = m_stats_stats_matrix[['assembly_accession', "coverage", "covbases", "meanmapq", "file"]]
        m_stats_stats_matrix['error_rate'] = 0.0

    matched = matched[matched['assembly_file'].notna()]
    matched['assembly_file'] = matched['assembly_file'].apply(lambda x: os.path.basename(x))
    matched = matched.sort_values(by=['total_uniq_reads'], ascending=False).drop_duplicates(subset=['assembly_accession'], keep= 'first')

    m_stats_stats_matrix = m_stats_stats_matrix.apply(merge_to_matched, axis=1, matched=matched)
    m_stats_stats_matrix = m_stats_stats_matrix[m_stats_stats_matrix['taxid'].notna()]
    m_stats_stats_matrix.drop(columns= ['file'], inplace=True)
    
    m_stats_stats_matrix = m_stats_stats_matrix.apply(lambda row: update_df_best_match(row, input_taxids, ncbi_wrapper), axis=1)
    m_stats_stats_matrix['leaf'] = m_stats_stats_matrix['assid'].apply(lambda x: match_leaf(x, overlap_manager.leaves))

    m_stats_stats_matrix = merge_by_assembly_ID(m_stats_stats_matrix)
    m_stats_stats_matrix = m_stats_stats_matrix.drop_duplicates(subset=['assid'])

    m_stats_stats_matrix = dataframe_update_with_lineage(m_stats_stats_matrix, ncbi_wrapper)

    mstats = []
    m_stats_stats_matrix['best_match_is_best'] = False
    for _, data in m_stats_stats_matrix.groupby('best_match_taxid'):

        data = data.sort_values(by=['best_match_score', 'coverage', 'error_rate'], ascending=[False, False, True])
        found = False
        data.loc[:, 'best_match_is_best'] = False
        for ix, row in data.iterrows():
            if row['coverage'] == 0.0:
                continue
            if not found:
                data.at[ix, 'best_match_is_best'] = True
                found = True

        mstats.append(data)

    mstats.append(m_stats_stats_matrix[m_stats_stats_matrix['best_match_taxid'].isna()])
    m_stats_stats_matrix = pd.concat(mstats, ignore_index=True)

    if m_stats_stats_matrix.empty:
        return pd.DataFrame()

    best_hits = m_stats_stats_matrix[m_stats_stats_matrix['best_match_is_best'] == True]

    m_stats_stats_matrix = m_stats_stats_matrix.apply(id_crosshit, axis=1, overlap_manager=overlap_manager, m_stats_stats_matrix=m_stats_stats_matrix, ncbi_wrapper=ncbi_wrapper,
                                                      best_hits=best_hits, cross_hit_threshold=cross_hit_threshold)
    m_stats_stats_matrix['max_shared'] = 1.0 
    m_stats_stats_matrix = m_stats_stats_matrix.apply(get_max_shared, overlap_manager=overlap_manager, best_hits=best_hits, axis=1)

    return m_stats_stats_matrix 



def node_total_true_leaves(overlap_manager, node, m_stats_stats_matrix) -> list:
    tree = overlap_manager.tree

    descendants = nx.descendants(tree, node)
    leaves = [n for n in descendants if tree.out_degree(n) == 0]
    if tree.out_degree(node) == 0:
        leaves.append(node)

    leaf_taxids = m_stats_stats_matrix[(m_stats_stats_matrix['leaf'].isin(leaves)) & (m_stats_stats_matrix['best_match_is_best'] == True)]['best_match_taxid'].dropna().tolist()
    
    return leaf_taxids


def node_leaves_best_taxids(overlap_manager, node, m_stats_stats_matrix) -> list:
    tree = overlap_manager.tree

    descendants = nx.descendants(tree, node)
    leaves = [n for n in descendants if tree.out_degree(n) == 0]
    if tree.out_degree(node) == 0:
        leaves.append(node)

    leaf_taxids = m_stats_stats_matrix[(m_stats_stats_matrix['leaf'].isin(leaves))]['best_match_taxid'].dropna().tolist()
    
    return leaf_taxids


def node_leaf_shannon_tax_diversity(overlap_manager, node, m_stats_stats_matrix, tax_level: str = "order") -> float:
    
    tree = overlap_manager.tree

    descendants = nx.descendants(tree, node)
    leaves = [n for n in descendants if tree.out_degree(n) == 0]
    if tree.out_degree(node) == 0:
        leaves.append(node)

    leaf_taxa = m_stats_stats_matrix[(m_stats_stats_matrix['leaf'].isin(leaves))][tax_level].dropna().tolist()

    return shannon_diversity_from_list(leaf_taxa)


def get_subset_composition(node_data, tax_data: pd.DataFrame, tax_level: str = "order"):
    valid_values = set(tax_data[tax_level].dropna().astype(str).tolist())
    node_data[tax_level] = (
        node_data[tax_level]
        .fillna('unclassified')
        .astype(str)
        .apply(lambda x: x if x in valid_values else 'unclassified')
    )

    composition = node_data[tax_level].value_counts(normalize=True).reset_index()

    composition.columns = ['tax_level', 'proportion']

    input_taxa = tax_data.drop_duplicates(subset=[tax_level])
    missing_orders = [tax for tax in input_taxa[tax_level] if tax not in composition['tax_level'].tolist()]

    composition = pd.concat([composition, pd.DataFrame({'tax_level': missing_orders, 'proportion': [0.0]*len(missing_orders)})], ignore_index=True)
    composition.loc[:,'tax_level'] = pd.Categorical(composition['tax_level'], categories=input_taxa[tax_level], ordered=True)
    composition = composition.sort_values('tax_level').reset_index(drop=True)
    return composition


def get_subset_composition_counts(node_data, tax_data: pd.DataFrame, tax_level: str = "order", count_column: str = "total_uniq_reads"):
    valid_values = set(tax_data[tax_level].dropna().astype(str).tolist())
    node_data[tax_level] = (
        node_data[tax_level]
        .fillna('unclassified')
        .astype(str)
        .apply(lambda x: x if x in valid_values else 'unclassified')
    )

    counts_by_taxa = node_data.groupby(tax_level)[count_column].sum().reset_index()
    total_counts = counts_by_taxa[count_column].sum()
    counts_by_taxa['proportion'] = counts_by_taxa[count_column] / total_counts
    input_taxa = tax_data.drop_duplicates(subset=[tax_level])
    input_taxa = input_taxa[input_taxa[tax_level].notna()]
    missing_orders = [tax for tax in input_taxa[tax_level] if tax not in counts_by_taxa[tax_level].tolist()]
    counts_by_taxa = pd.concat([counts_by_taxa, pd.DataFrame({'order': missing_orders, 'total_uniq_reads': [0]*len(missing_orders), 'proportion': [0.0]*len(missing_orders)})], ignore_index=True)

    counts_by_taxa.loc[:,'tax_level'] = pd.Categorical(counts_by_taxa[tax_level], categories=input_taxa[tax_level], ordered=True)
    counts_by_taxa = counts_by_taxa.sort_values('tax_level').reset_index(drop=True)

    return counts_by_taxa


def node_composition_level(overlap_manager, node, m_stats_stats_matrix, tax_data: pd.DataFrame, tax_level: str = "order"):
    node_leaves = overlap_manager.get_node_leaves(node)
    node_data = m_stats_stats_matrix[m_stats_stats_matrix['leaf'].isin(node_leaves)].copy()
    composition = get_subset_composition(node_data, tax_data, tax_level=tax_level)
    return composition



def get_composition_by_leaf(overlap_manager, m_stats_stats_matrix, tax_data: pd.DataFrame, tax_level: str = "order"):
    leaves = overlap_manager.leaves
    compositions = []
    for leaf in leaves:
        composition = node_composition_level(overlap_manager, leaf, m_stats_stats_matrix, tax_data, tax_level=tax_level).set_index('tax_level').T
        composition.insert(0, 'leaf', leaf)
        compositions.append(composition)

    return pd.concat(compositions, axis=0)



def compute_node_stats(overlap_manager) -> pd.DataFrame:
    all_node_stats = overlap_manager.all_node_stats.copy()
    dist_cache = {n: nx.single_source_shortest_path_length(overlap_manager.tree.to_undirected(), n) for n in overlap_manager.all_nodes}
    summary_node_stats = pd.concat([clade_graph_metrics(overlap_manager, node, dist_cache) for node in overlap_manager.all_nodes], ignore_index=True)
    summary_node_stats = all_node_stats.merge(summary_node_stats, on='Node', how='left')

    return summary_node_stats

def compute_node_purity(overlap_manager, m_stats_stats_matrix) -> pd.DataFrame:
    tree = overlap_manager.tree
    all_node_stats = overlap_manager.all_node_stats.copy()
    all_node_stats['leaf'] = all_node_stats['Node'].apply(lambda x: overlap_manager.get_node_leaves(x))
    all_node_stats = all_node_stats.explode('leaf')

    all_node_stats = all_node_stats.merge(m_stats_stats_matrix[['leaf', 'best_match_taxid']], on='leaf', how='left')

    node_purity = {}

    for node in tree.nodes:
        descendants = nx.descendants(tree, node)
        leaves = [n for n in descendants if tree.out_degree(n) == 0]
        if tree.out_degree(node) == 0:
            leaves.append(node)

        leaf_taxids = all_node_stats[all_node_stats['leaf'].isin(leaves)]['best_match_taxid'].dropna().tolist()
        total_leaves = len(leaf_taxids)
        if total_leaves == 0:
            node_purity[node] = (0, 0.0)
            continue

        taxid_counts = Counter(leaf_taxids)
        most_common_taxid, most_common_count = taxid_counts.most_common(1)[0]
        purity = most_common_count / total_leaves

        node_purity[node] = (most_common_taxid, purity)

    purity_df = pd.DataFrame.from_dict(node_purity, orient='index', columns=['Most_Common_TaxID', 'Purity'])
    purity_df.reset_index(inplace=True)
    purity_df.rename(columns={'index': 'Node'}, inplace=True)

    return purity_df


def all_antichains_covering_leaves(overlap_manager, purity_df: pd.DataFrame) -> list[list]:
    antichains_keep = []
    tree = overlap_manager.tree
    for antichain in nx.antichains(tree):
        anti_chain_df = pd.DataFrame(antichain, columns=['Node'])
        anti_chain_df = anti_chain_df.merge(purity_df, on='Node', how='left')
        if sum(anti_chain_df['Num_Leaves'].fillna(0)) < len(overlap_manager.leaves):
            continue
        node_left = purity_df[~purity_df['Node'].isin(anti_chain_df['Node'])]
        precision = anti_chain_df.drop_duplicates(subset = ['Node', 'Most_Common_TaxID'])['Node'].nunique() / purity_df['Most_Common_TaxID'].nunique()
        precision_balanced = -1 * abs(precision - 1)
        nodes = anti_chain_df['Node'].tolist()
        internal_nodes = anti_chain_df[anti_chain_df['nleaves'] > 1]
        min_pair_dist_internal_nodes = internal_nodes['Min_Pairwise_Dist'].min() if not internal_nodes.empty else 1 
        max_pair_dist_internal_nodes = internal_nodes['Min_Pairwise_Dist'].max() if not internal_nodes.empty else 0
        min_shared_max = internal_nodes['Min_Shared'].max() if 'Min_Shared' in internal_nodes.columns else 0
        parents = [tree.predecessors(n) for n in nodes if tree.in_degree(n) > 0]
        parents = [item for sublist in parents for item in sublist]
        parents_min_pair_dist = purity_df[purity_df['Node'].isin(parents)]['Min_Pairwise_Dist'].max()
        
        anti_chain_df = anti_chain_df.aggregate(
            {
                'Num_Leaves': 'sum',
                'Private_Reads': 'sum',
                'Min_Pairwise_Dist': 'max',
                'Private_Proportion': 'mean',
                'Purity': 'mean',
                'Min_Pairwise_Dist_tree': 'min',
                'Avg_Pairwise_Dist': 'mean',
                'Med_Pairwise_Dist': 'mean',
                'Dist_to_Selected': 'min',
                'Dist_to_NonSelected': 'min',
            }, axis=0
        )

        anti_chain_df['nnodes'] = len(nodes)
        anti_chain_df['Node'] = list(nodes)
        anti_chain_df['Precision'] = precision
        anti_chain_df['MinPDist_internal'] = min_pair_dist_internal_nodes
        anti_chain_df['MaxPDist_internal'] = max_pair_dist_internal_nodes
        anti_chain_df['MaxShared_internal'] = min_shared_max
        anti_chain_df['Parents_MaxPDist'] = parents_min_pair_dist
        anti_chain_df['Precision_Balanced'] = precision_balanced

        antichains_keep.append(anti_chain_df)
        
    return antichains_keep



def clade_graph_metrics(overlap_manager, node, dist_cache: dict) -> pd.DataFrame:
    selected_nodes = overlap_manager.get_node_leaves(node)
    sel_list = list(selected_nodes)
    pairwise = []
    for i in range(len(sel_list)):
        for j in range(i + 1, len(sel_list)):
            n1 = sel_list[i]
            n2 = sel_list[j]
            dist = dist_cache[n1].get(n2, np.inf)
            pairwise.append(dist)

    if len(pairwise) == 0:
        avg_pairwise = 0
        min_pairwise = 0
        med_pairwise = 0
    else:
        arr = np.array(pairwise)
        avg_pairwise = arr.mean()
        min_pairwise = arr.min()
        med_pairwise = np.median(arr)
    
    dist_to_selected = min((dist_cache[node][s] for s in selected_nodes), default=np.inf)
    dist_to_nonselected = min((dist_cache[node][m] for m in overlap_manager.all_nodes if m not in selected_nodes), default=np.inf)
    average_coverage = overlap_manager.m_stats_matrix[overlap_manager.m_stats_matrix.index.isin(selected_nodes)]['coverage'].mean()
    private_reads = overlap_manager.all_node_stats[overlap_manager.all_node_stats['Node'].isin(selected_nodes)]['Private_Reads'].sum()

    df = pd.DataFrame({
        'Node': [node],
        'nleaves': [len(selected_nodes)],
        'nuniq': [len(set(selected_nodes))],
        'NuniqReads': [private_reads],
        'Avg_Coverage': [average_coverage],
        'Min_Pairwise_Dist_tree': [min_pairwise],
        'Avg_Pairwise_Dist': [avg_pairwise],
        'Med_Pairwise_Dist': [med_pairwise],
        'Dist_to_Selected': [dist_to_selected],
        'Dist_to_NonSelected': [dist_to_nonselected]
    })
    return df


def antichain_classification(overlap_manager, m_stats_stats_matrix: pd.DataFrame) -> pd.DataFrame:
    summary_node_stats = compute_node_stats(overlap_manager)
    purity_df = compute_node_purity(overlap_manager, m_stats_stats_matrix)
    purity_df = purity_df.merge(summary_node_stats, on='Node', how='left')

    antichains = all_antichains_covering_leaves(overlap_manager, purity_df)
    antichains_df = pd.DataFrame(antichains)
    antichains_df.sort_values(by=['Precision_Balanced', 'MinPDist_internal', 'MaxShared_internal', 'Purity', 'Private_Reads'], ascending=[False, False, False, False], inplace=True)
    antichains_df = antichains_df.reset_index(drop=True)

    antichain_selected = antichains_df[(antichains_df['Precision_Balanced'] == 0.0) & (antichains_df['MaxPDist_internal'] > 0.05) & (antichains_df['MaxShared_internal'] <= 0.8)]
    if antichain_selected.empty is True:
        antichain_selected = pd.DataFrame(antichains_df.iloc[[0]])

    antichain_selected = antichain_selected.explode('Node').reset_index(drop=True)

    classified_purity= purity_df.copy()

    def selected_condition(row):
        if row['nleaves'] > 1 and row['Min_Pairwise_Dist'] <= 0.05:
            return False
        return row['Node'] in antichain_selected['Node'].values

    classified_purity['selected'] = classified_purity.apply(selected_condition, axis=1)
    classified_purity.sort_values(by=['selected', 'Purity', 'Min_Pairwise_Dist'], ascending=[False, False, False], inplace=True)
    internal_nodes = classified_purity[classified_nodes > 1]
    classified_purity['Z_min_dist'] = (internal_nodes['Min_Pairwise_Dist'] - internal_nodes['Min_Pairwise_Dist'].mean()) / internal_nodes['Min_Pairwise_Dist'].std() if internal_nodes['Min_Pairwise_Dist'].std() > 0 else 0

    classified_purity = classified_purity.reset_index(drop=True)

    return classified_purity


def find_assembly_mapping(row, stats_matrix):
    accession = row['assembly_accession']
    if accession is None or pd.isna(accession):
        row['clade'] = 'unmapped'
        row['nuniq'] = 0
        row['freq'] = 0
        row['Min_Pairwise_Dist'] = 0
        row['nleaves'] = 0
        return row

    match = stats_matrix[stats_matrix['leaves'].str.contains(accession, na=False) | (stats_matrix['clade'].str.contains(accession, na=False))]

    if match.empty:
        row['clade'] = None
        row['nuniq'] = 0
        row['freq'] = 0
        row['Min_Pairwise_Dist'] = 0
        row['nleaves'] = 0
    else:
        row['clade'] = match['clade'].values[0]
        row['nuniq'] = match['nuniq'].values[0]
        row['freq'] = match['freq'].values[0]
        row['Min_Pairwise_Dist'] = match['Min_Pairwise_Dist'].values[0]
        row['nleaves'] = match['nleaves'].values[0]

    return row


def find_best_match(taxid1, taxid_list, ncbi_wrapper: NCBITaxonomistWrapper):
    best_taxid = None
    best_level = None
    best_score = 0.0
    for taxid2 in taxid_list:
        score, level = ncbi_wrapper.compare_lineages_relative(taxid2, taxid1)
        if level is not None and score > best_score:
            best_level = level
            best_taxid = taxid2
            best_score = score
    return best_taxid, best_level, best_score


def update_df_best_match(row, taxid_list, ncbi_wrapper: NCBITaxonomistWrapper, min_score = 0.7):
    from typing import Optional
    taxid = row['taxid']
    if pd.isna(taxid):
        row['best_match_taxid'] = None
        row['best_match_level'] = -1
        row['best_match_score'] = 0.0
        row['name'] = None

    best_taxid, best_level, best_score = find_best_match(taxid, taxid_list, ncbi_wrapper)

    if best_score < min_score:
        best_taxid = None
        best_level = -1
        best_score = 0.0

    row['best_match_taxid'] = best_taxid
    row['best_match_level'] = best_level
    row['best_match_score'] = best_score
    
    if best_taxid is not None:
        row['name'] = ncbi_wrapper.get_name(best_taxid)
    else:
        row['name'] = None

    return row


def match_leaf(accid, list):
    for item in list:
        if accid in item:
            return item
    return None
