from __future__ import annotations

import os
from collections import Counter

import numpy as np
import pandas as pd

from metagenomics_utils.overlap_manager.diversity import (
    kurtosis,
    shannon_diversity_from_list,
    skewness,
)

################################################# node stats and composition functions #################################################


def normalize_by_taxlevel(prediction_matrix: pd.DataFrame, tax_level: str = "genus"):
    normalized = prediction_matrix.copy()
    stats_cols = ["coverage", "covbases", "meanmapq", "error_rate"]

    for col in stats_cols:
        normalized[col] = normalized[col].astype(float)
        for _, group in prediction_matrix.groupby(tax_level):
            idx = group.index

            col_mean = group[col].mean()
            col_std = group[col].std()
            if col_std > 0:
                normalized.loc[idx, col] = (group[col] - col_mean) / col_std
            else:
                normalized.loc[idx, col] = 0.0

    return normalized


def id_crosshit(
    row,
    overlap_manager: OverlapManager,
    m_stats_stats_matrix: pd.DataFrame,
    ncbi_wrapper: NCBITaxonomistWrapper,
    best_hits: pd.DataFrame,
    cross_hit_threshold: float = 0.3,
):
    row["is_crosshit"] = False
    row["cross_hit_match"] = None
    row["crosshit_match_score"] = 0.0
    row["crosshit_match_level"] = None
    row["crossh_hit_dist"] = 1.0

    if row["best_match_is_best"] is True:
        return row

    if row["leaf"] not in overlap_manager.distance_mat.index:
        return row

    best_match = overlap_manager.distance_mat.loc[row["leaf"]]

    if best_match.empty or len(best_hits) == 0:
        return row

    if best_match[best_match.index.isin(best_hits["leaf"])].empty:
        return row

    best_match_dist = best_match[best_match.index.isin(best_hits["leaf"])].min()
    # index of all best matches
    # best_match = best_match[best_match.index.isin(best_hits['leaf'])].idxmin()
    best_match = best_match[best_match.index.isin(best_hits["leaf"])].idxmin()

    best_match_taxid = m_stats_stats_matrix[m_stats_stats_matrix["leaf"] == best_match]["best_match_taxid"].values
    if len(best_match_taxid) > 0:
        best_match_taxid = best_match_taxid[0]
        score, level = ncbi_wrapper.compare_lineages_relative(best_match_taxid, row["taxid"])
        if best_match_dist < cross_hit_threshold:
            row["is_crosshit"] = True
        row["cross_hit_match"] = best_match_taxid
        row["crosshit_match_score"] = score
        row["crosshit_match_level"] = level
        row["crossh_hit_dist"] = best_match_dist

    return row


def get_max_shared(row, overlap_manager, best_hits: pd.DataFrame):
    if row["leaf"] not in overlap_manager.distance_mat.index:
        return row

    best_match = overlap_manager.distance_mat.loc[row["leaf"]]
    if best_match.empty or len(best_hits) == 0:
        return row

    if best_match[best_match.index != row["leaf"]].empty:
        return row

    best_match_dist = best_match[best_match.index != row["leaf"]].min()

    row["max_shared"] = 1.0 - best_match_dist
    return row


def dataframe_update_with_lineage(df: pd.DataFrame, ncbi_wrapper: NCBITaxonomistWrapper) -> pd.DataFrame:
    df["order"] = df.apply(lambda row: ncbi_wrapper.get_level(row["taxid"], "order"), axis=1)
    df["order"] = df["order"].fillna("unclassified")
    df["family"] = df.apply(lambda row: ncbi_wrapper.get_level(row["taxid"], "family"), axis=1)
    df["family"] = df["family"].fillna("unclassified")
    df["genus"] = df.apply(lambda row: ncbi_wrapper.get_level(row["taxid"], "genus"), axis=1)
    df["genus"] = df["genus"].fillna("unclassified")
    return df


def get_m_stats_matrix(
    data_set_name,
    study_output_filepath,
    ncbi_wrapper: NCBITaxonomistWrapper,
    overlap_manager: OverlapManager,
    cross_hit_threshold: float = 0.3,
    min_taxonomic_score: float = 0.7,
    best_match_level: str = "species",
    filter_no_leaf=True,
    log=False,
):
    """
    Build the m-stats matrix by assembling coverage statistics, classification
    results, best-match annotation, cross-hit identification, and trash labelling
    into a single DataFrame.

    min_taxonomic_score = minimum taxonomic match.
    """
    paths = _get_m_stats_paths(data_set_name, study_output_filepath)
    if not _validate_m_stats_paths(paths):
        return pd.DataFrame()

    input_df, matched, m_stats, classification = _load_m_stats_inputs(paths)
    matched = _resolve_assembly_classifications(matched, classification)
    m_stats = _merge_matched_stats(m_stats, matched)
    if log:
        print(f"  Dataset {data_set_name}: {len(m_stats)} matched assemblies with taxid")
    if m_stats.empty:
        return pd.DataFrame()

    input_taxids = input_df["taxid"].dropna().unique().tolist()
    m_stats = _compute_best_matches(m_stats, input_taxids, ncbi_wrapper, overlap_manager, min_taxonomic_score)
    if log:
        print(f"  Dataset {data_set_name}: {len(m_stats)} assemblies with best-match taxid")
    if m_stats.empty:
        return pd.DataFrame()

    m_stats = _flag_cross_hits(m_stats, overlap_manager, ncbi_wrapper, cross_hit_threshold)
    if log:
        print(m_stats.head())
        print(m_stats["leaf"])
        print(f"  Dataset {data_set_name}: {len(m_stats)} assemblies after cross-hit and trash filtering")
    m_stats = _finalize_m_stats(m_stats, filter_no_leaf)
    if log:
        print(m_stats.head())
        print(m_stats["leaf"])
        print(f"  Dataset {data_set_name}: {len(m_stats)} assemblies after cross-hit and trash filtering")
    return m_stats


# ---------------------------------------------------------------------------
# Pipeline stage helpers
# ---------------------------------------------------------------------------


def _get_m_stats_paths(data_set_name, study_output_filepath):
    """Return dict of file paths for a single dataset's m-stats pipeline."""
    return {
        "matched_assemblies": os.path.join(study_output_filepath, data_set_name, "output", "matched_assemblies.tsv"),
        "merged_stats": os.path.join(study_output_filepath, data_set_name, "output", "merged_coverage_statistics.tsv"),
        "input": os.path.join(study_output_filepath, data_set_name, "input", f"{data_set_name}.tsv"),
        "classification": os.path.join(
            study_output_filepath, data_set_name, "classification", f"{data_set_name}_merged_classification.tsv"
        ),
    }


def _validate_m_stats_paths(paths):
    """Return True when all required input files exist."""
    return all(os.path.exists(paths[k]) for k in ("matched_assemblies", "merged_stats", "input"))


def _load_m_stats_inputs(paths):
    """Read and return (input_df, matched, m_stats, classification) DataFrames."""
    input_df = pd.read_csv(paths["input"], sep="\t")
    matched = pd.read_csv(paths["matched_assemblies"], sep="\t")
    classification = pd.read_csv(paths["classification"], sep="\t")

    m_stats = pd.read_csv(paths["merged_stats"], sep="\t").rename(columns={"#rname": "assembly_accession"})
    if "error_rate" in m_stats.columns:
        m_stats = m_stats[["assembly_accession", "coverage", "covbases", "meanmapq", "numreads", "error_rate", "file"]]
    else:
        m_stats = m_stats[["assembly_accession", "coverage", "covbases", "meanmapq", "numreads", "file"]]
        m_stats["error_rate"] = 0.0

    return input_df, matched, m_stats, classification


def _resolve_assembly_classifications(matched, classification):
    """Annotate matched assemblies with classification info — taxid, reads, classifiers."""

    def _classify_row(row):
        taxid = row["taxid"]
        match = classification[classification["taxid"] == taxid]
        if match.empty:
            row["classifiers"] = "unclassified"
            row["total_uniq_reads"] = 0
        else:
            if "classifiers" in match.columns:
                row["classifiers"] = match["classifiers"].values[0]
            elif "classification" in match.columns:
                row["classifiers"] = match["classification"].values[0]
            else:
                raise ValueError("No classification found")
            if "total_uniq_reads" in match.columns:
                row["total_uniq_reads"] = match["total_uniq_reads"].values[0]
            elif "uniq_reads" in match.columns:
                row["total_uniq_reads"] = match["uniq_reads"].values[0]
            else:
                row["total_uniq_reads"] = 0
        return row

    matched = matched[matched["assembly_file"].notna()]
    matched = matched.apply(_classify_row, axis=1)
    matched["assembly_file"] = matched["assembly_file"].apply(lambda x: os.path.basename(x))
    matched = matched.sort_values(by=["total_uniq_reads"], ascending=False).drop_duplicates(
        subset=["assembly_accession"], keep="first"
    )
    return matched


def _merge_matched_stats(m_stats, matched):
    """Merge coverage stats with matched assemblies; keep only rows with a taxid."""
    from metagenomics_utils.overlap_manager.manager import _merge_matched_vectorized

    m_stats = _merge_matched_vectorized(m_stats, matched)
    m_stats = m_stats[m_stats["taxid"].notna()]
    m_stats.drop(columns=["file"], inplace=True)
    return m_stats


def _compute_best_matches(m_stats, input_taxids, ncbi_wrapper, overlap_manager, min_taxonomic_score):
    """Annotate best-match taxid, leaf, lineage, and mark best-of-group."""
    from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper
    from metagenomics_utils.overlap_manager.manager import merge_by_assembly_ID

    m_stats = m_stats.apply(
        lambda row: update_df_best_match(
            row, input_taxids, ncbi_wrapper, min_score=min_taxonomic_score, best_level=NCBITaxonomistWrapper.TAX_SPECIES
        ),
        axis=1,
    )
    print("updated")

    m_stats["leaf"] = m_stats["assid"].apply(lambda x: match_leaf(x, overlap_manager.leaves))

    m_stats = merge_by_assembly_ID(m_stats)
    m_stats = m_stats.drop_duplicates(subset=["assid"])
    m_stats = dataframe_update_with_lineage(m_stats, ncbi_wrapper)
    print("updated lineage")

    # mark one best match per group
    ncbi_tools = NCBITaxonomistWrapper()
    groups = []
    m_stats["best_match_is_best"] = False
    for _, group in m_stats.groupby("best_match_taxid"):
        group = group.sort_values(
            by=["best_match_score", "coverage", "error_rate"], ascending=[False, False, True]
        ).reset_index(drop=True)
        found = False
        for ix, row in group.iterrows():
            if row["coverage"] == 0.0:
                continue
            if ncbi_tools.level_is_below(row["best_match_level"], NCBITaxonomistWrapper.TAX_SPECIES):
                continue
            if not found:
                group.at[ix, "best_match_is_best"] = True
                found = True
        groups.append(group)

    groups.append(m_stats[m_stats["best_match_taxid"].isna()])
    m_stats = pd.concat(groups, ignore_index=True)

    return m_stats


def _flag_cross_hits(m_stats, overlap_manager, ncbi_wrapper, cross_hit_threshold):
    """Identify cross-hits and compute max shared proportion for each leaf."""
    best_hits = m_stats[m_stats["best_match_is_best"] == True]

    m_stats = m_stats.apply(
        id_crosshit,
        axis=1,
        overlap_manager=overlap_manager,
        m_stats_stats_matrix=m_stats,
        ncbi_wrapper=ncbi_wrapper,
        best_hits=best_hits,
        cross_hit_threshold=cross_hit_threshold,
    )
    m_stats["max_shared"] = 1.0
    m_stats = m_stats.apply(get_max_shared, axis=1, overlap_manager=overlap_manager, best_hits=best_hits)
    return m_stats


def _finalize_m_stats(m_stats, filter_no_leaf):
    """Apply leaf/assid indexing, sort, and label trash rows."""
    if filter_no_leaf:
        m_stats = m_stats[m_stats["leaf"].notna()]
        m_stats.set_index("leaf", inplace=True)
    else:
        m_stats.set_index("assid", inplace=True)
    m_stats = m_stats.sort_values("total_uniq_reads", ascending=False)

    def is_trash(row):
        if pd.isna(row["best_match_taxid"]):
            return True
        return not row["best_match_is_best"] and not row["is_crosshit"]

    m_stats["is_trash"] = m_stats.apply(is_trash, axis=1)

    # move total_uniq_reads to last column
    cols = m_stats.columns.tolist()
    cols.remove("total_uniq_reads")
    cols.append("total_uniq_reads")
    m_stats = m_stats[cols]

    return m_stats


def node_total_true_leaves(overlap_manager, node, m_stats_stats_matrix) -> list:
    import networkx as nx
    tree = overlap_manager.tree

    descendants = nx.descendants(tree, node)
    leaves = [n for n in descendants if tree.out_degree(n) == 0]

    if tree.out_degree(node) == 0:
        leaves.append(node)

    leaf_taxids = (
        m_stats_stats_matrix[
            (m_stats_stats_matrix.index.isin(leaves)) & (m_stats_stats_matrix["best_match_is_best"] == True)
        ]["best_match_taxid"]
        .dropna()
        .tolist()
    )

    return leaf_taxids


def node_leaves_best_taxids(overlap_manager, node, m_stats_stats_matrix) -> list:
    import networkx as nx
    tree = overlap_manager.tree

    descendants = nx.descendants(tree, node)
    leaves = [n for n in descendants if tree.out_degree(n) == 0]
    if tree.out_degree(node) == 0:
        leaves.append(node)

    leaf_taxids = m_stats_stats_matrix[(m_stats_stats_matrix.index.isin(leaves))]["best_match_taxid"].dropna().tolist()

    return leaf_taxids


def node_leaf_shannon_tax_diversity(overlap_manager, node, m_stats_stats_matrix, tax_level: str = "order") -> float:
    import networkx as nx

    tree = overlap_manager.tree

    descendants = nx.descendants(tree, node)
    leaves = [n for n in descendants if tree.out_degree(n) == 0]
    if tree.out_degree(node) == 0:
        leaves.append(node)

    leaf_taxa = m_stats_stats_matrix[(m_stats_stats_matrix.index.isin(leaves))][tax_level].dropna().tolist()

    return shannon_diversity_from_list(leaf_taxa)


def get_subset_composition(node_data, tax_data: pd.DataFrame, tax_level: str = "order"):
    valid_values = set(tax_data[tax_level].dropna().astype(str).tolist())
    node_data[tax_level] = (
        node_data[tax_level]
        .fillna("unclassified")
        .astype(str)
        .apply(lambda x: x if x in valid_values else "unclassified")
    )

    composition = node_data[tax_level].value_counts(normalize=True).reset_index()

    composition.columns = ["tax_level", "proportion"]

    input_taxa = tax_data.drop_duplicates(subset=[tax_level])
    missing_orders = [tax for tax in input_taxa[tax_level] if tax not in composition["tax_level"].tolist()]

    composition = pd.concat(
        [composition, pd.DataFrame({"tax_level": missing_orders, "proportion": [0.0] * len(missing_orders)})],
        ignore_index=True,
    )
    composition.loc[:, "tax_level"] = pd.Categorical(
        composition["tax_level"], categories=input_taxa[tax_level], ordered=True
    )
    composition = composition.sort_values("tax_level").reset_index(drop=True)
    return composition


def get_subset_composition_counts(
    node_data, tax_data: pd.DataFrame, tax_level: str = "order", count_column: str = "total_uniq_reads"
):
    valid_values = set(tax_data[tax_level].dropna().astype(str).tolist())
    node_data[tax_level] = (
        node_data[tax_level]
        .fillna("unclassified")
        .astype(str)
        .apply(lambda x: x if x in valid_values else "unclassified")
    )

    counts_by_taxa = node_data.groupby(tax_level)[count_column].sum().reset_index()
    total_counts = counts_by_taxa[count_column].sum()
    counts_by_taxa["proportion"] = counts_by_taxa[count_column] / total_counts
    input_taxa = tax_data.drop_duplicates(subset=[tax_level])
    input_taxa = input_taxa[input_taxa[tax_level].notna()]
    missing_orders = [tax for tax in input_taxa[tax_level] if tax not in counts_by_taxa[tax_level].tolist()]
    counts_by_taxa = pd.concat(
        [
            counts_by_taxa,
            pd.DataFrame(
                {
                    tax_level: missing_orders,
                    "total_uniq_reads": [0] * len(missing_orders),
                    "proportion": [0.0] * len(missing_orders),
                }
            ),
        ],
        ignore_index=True,
    )

    counts_by_taxa.loc[:, "tax_level"] = pd.Categorical(
        counts_by_taxa[tax_level], categories=input_taxa[tax_level], ordered=True
    )
    counts_by_taxa = counts_by_taxa.sort_values("tax_level").reset_index(drop=True)

    return counts_by_taxa


def node_composition_with_stats(
    node_data, tax_data: pd.DataFrame, tax_level: str = "order", count_column: str = "total_uniq_reads"
):
    composition = get_subset_composition_counts(node_data, tax_data, tax_level=tax_level, count_column=count_column)
    composition = composition[["tax_level", "proportion"]].set_index("tax_level").T
    tax_diversity_shannon = shannon_diversity_from_list(node_data[tax_level].dropna().tolist())

    counts_skewness = skewness(node_data["total_uniq_reads"].tolist())
    counts_kurtosis = kurtosis(node_data["total_uniq_reads"].tolist())
    number_hits = len(node_data)

    composition.insert(0, "total_uniq_reads", node_data["total_uniq_reads"].sum())
    composition.insert(0, "max_uniq_reads", node_data["total_uniq_reads"].max())
    composition.insert(0, "tax_diversity_shannon", tax_diversity_shannon)
    composition.insert(0, "counts_skewness", counts_skewness)
    composition.insert(0, "counts_kurtosis", counts_kurtosis)
    composition.insert(0, "number_hits", number_hits)

    return composition


def node_composition_level(
    overlap_manager, node, m_stats_stats_matrix, tax_data: pd.DataFrame, tax_level: str = "order"
):
    node_leaves = overlap_manager.get_node_leaves(node)
    node_data = m_stats_stats_matrix[m_stats_stats_matrix.index.isin(node_leaves)].copy()
    composition = get_subset_composition(node_data, tax_data, tax_level=tax_level)
    return composition


def get_composition_by_leaf(overlap_manager, m_stats_stats_matrix, tax_data: pd.DataFrame, tax_level: str = "order"):
    leaves = overlap_manager.leaves
    compositions = []
    for leaf in leaves:
        composition = (
            node_composition_level(overlap_manager, leaf, m_stats_stats_matrix, tax_data, tax_level=tax_level)
            .set_index("tax_level")
            .T
        )
        composition.insert(0, "leaf", leaf)
        compositions.append(composition)

    return pd.concat(compositions, axis=0)


def compute_node_stats(overlap_manager) -> pd.DataFrame:
    import networkx as nx
    all_node_stats = overlap_manager.all_node_stats.copy()
    dist_cache = {
        n: nx.single_source_shortest_path_length(overlap_manager.tree.to_undirected(), n)
        for n in overlap_manager.all_nodes
    }
    summary_node_stats = pd.concat(
        [clade_graph_metrics(overlap_manager, node, dist_cache) for node in overlap_manager.all_nodes],
        ignore_index=True,
    )
    summary_node_stats = all_node_stats.merge(summary_node_stats, on="Node", how="left")

    return summary_node_stats


def compute_node_purity(overlap_manager, m_stats_stats_matrix) -> pd.DataFrame:
    tree = overlap_manager.tree
    all_node_stats = overlap_manager.all_node_stats.copy()
    all_node_stats["leaf"] = all_node_stats["Node"].apply(lambda x: overlap_manager.get_node_leaves(x))
    all_node_stats = all_node_stats.explode("leaf")
    m_stats_stats_matrix["leaf"] = m_stats_stats_matrix.apply(
        lambda row: overlap_manager.get_leaf(row), axis=1
    )
    all_node_stats = all_node_stats.merge(m_stats_stats_matrix[["leaf", "best_match_taxid"]], on="leaf", how="left")

    node_purity = {}

    for node in tree.nodes:
        import networkx as nx
        descendants = nx.descendants(tree, node)
        leaves = [n for n in descendants if tree.out_degree(n) == 0]
        if tree.out_degree(node) == 0:
            leaves.append(node)

        leaf_taxids = all_node_stats[all_node_stats["leaf"].isin(leaves)]["best_match_taxid"].dropna().tolist()
        total_leaves = len(leaf_taxids)
        if total_leaves == 0:
            node_purity[node] = (0, 0.0)
            continue

        taxid_counts = Counter(leaf_taxids)
        most_common_taxid, most_common_count = taxid_counts.most_common(1)[0]
        purity = most_common_count / total_leaves

        node_purity[node] = (most_common_taxid, purity)

    purity_df = pd.DataFrame.from_dict(node_purity, orient="index", columns=["Most_Common_TaxID", "Purity"])
    purity_df.reset_index(inplace=True)
    purity_df.rename(columns={"index": "Node"}, inplace=True)

    return purity_df


def all_antichains_covering_leaves(overlap_manager, purity_df: pd.DataFrame) -> list[list]:
    antichains_keep = []
    import networkx as nx
    tree = overlap_manager.tree
    for antichain in nx.antichains(tree):
        anti_chain_df = pd.DataFrame(antichain, columns=["Node"])
        anti_chain_df = anti_chain_df.merge(purity_df, on="Node", how="left")
        if sum(anti_chain_df["Num_Leaves"].fillna(0)) < len(overlap_manager.leaves):
            continue
        node_left = purity_df[~purity_df["Node"].isin(anti_chain_df["Node"])]
        precision = (
            anti_chain_df.drop_duplicates(subset=["Node", "Most_Common_TaxID"])["Node"].nunique()
            / purity_df["Most_Common_TaxID"].nunique()
        )
        precision_balanced = -1 * abs(precision - 1)
        nodes = anti_chain_df["Node"].tolist()
        internal_nodes = anti_chain_df[anti_chain_df["nleaves"] > 1]
        min_pair_dist_internal_nodes = internal_nodes["Min_Pairwise_Dist"].min() if not internal_nodes.empty else 1
        max_pair_dist_internal_nodes = internal_nodes["Min_Pairwise_Dist"].max() if not internal_nodes.empty else 0
        min_shared_max = internal_nodes["Min_Shared"].max() if "Min_Shared" in internal_nodes.columns else 0
        parents = [tree.predecessors(n) for n in nodes if tree.in_degree(n) > 0]
        parents = [item for sublist in parents for item in sublist]
        parents_min_pair_dist = purity_df[purity_df["Node"].isin(parents)]["Min_Pairwise_Dist"].max()

        anti_chain_df = anti_chain_df.aggregate(
            {
                "Num_Leaves": "sum",
                "Private_Reads": "sum",
                "Min_Pairwise_Dist": "max",
                "Private_Proportion": "mean",
                "Purity": "mean",
                "Min_Pairwise_Dist_tree": "min",
                "Avg_Pairwise_Dist": "mean",
                "Med_Pairwise_Dist": "mean",
                "Dist_to_Selected": "min",
                "Dist_to_NonSelected": "min",
            },
            axis=0,
        )

        anti_chain_df["nnodes"] = len(nodes)
        anti_chain_df["Node"] = list(nodes)
        anti_chain_df["Precision"] = precision
        anti_chain_df["MinPDist_internal"] = min_pair_dist_internal_nodes
        anti_chain_df["MaxPDist_internal"] = max_pair_dist_internal_nodes
        anti_chain_df["MaxShared_internal"] = min_shared_max
        anti_chain_df["Parents_MaxPDist"] = parents_min_pair_dist
        anti_chain_df["Precision_Balanced"] = precision_balanced

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
    dist_to_nonselected = min(
        (dist_cache[node][m] for m in overlap_manager.all_nodes if m not in selected_nodes), default=np.inf
    )
    average_coverage = overlap_manager.m_stats_matrix[overlap_manager.m_stats_matrix.index.isin(selected_nodes)][
        "coverage"
    ].mean()
    private_reads = overlap_manager.all_node_stats[overlap_manager.all_node_stats["Node"].isin(selected_nodes)][
        "Private_Reads"
    ].sum()

    df = pd.DataFrame(
        {
            "Node": [node],
            "nleaves": [len(selected_nodes)],
            "nuniq": [len(set(selected_nodes))],
            "NuniqReads": [private_reads],
            "Avg_Coverage": [average_coverage],
            "Min_Pairwise_Dist_tree": [min_pairwise],
            "Avg_Pairwise_Dist": [avg_pairwise],
            "Med_Pairwise_Dist": [med_pairwise],
            "Dist_to_Selected": [dist_to_selected],
            "Dist_to_NonSelected": [dist_to_nonselected],
        }
    )
    return df


def antichain_classification(overlap_manager, m_stats_stats_matrix: pd.DataFrame) -> pd.DataFrame:
    summary_node_stats = compute_node_stats(overlap_manager)
    purity_df = compute_node_purity(overlap_manager, m_stats_stats_matrix)
    purity_df = purity_df.merge(summary_node_stats, on="Node", how="left")

    antichains = all_antichains_covering_leaves(overlap_manager, purity_df)
    antichains_df = pd.DataFrame(antichains)
    antichains_df.sort_values(
        by=["Precision_Balanced", "MinPDist_internal", "MaxShared_internal", "Purity", "Private_Reads"],
        ascending=[False, False, False, False],
        inplace=True,
    )
    antichains_df = antichains_df.reset_index(drop=True)

    antichain_selected = antichains_df[
        (antichains_df["Precision_Balanced"] == 0.0)
        & (antichains_df["MaxPDist_internal"] > 0.05)
        & (antichains_df["MaxShared_internal"] <= 0.8)
    ]
    if antichain_selected.empty is True:
        antichain_selected = pd.DataFrame(antichains_df.iloc[[0]])

    antichain_selected = antichain_selected.explode("Node").reset_index(drop=True)

    classified_purity = purity_df.copy()

    def selected_condition(row):
        if row["nleaves"] > 1 and row["Min_Pairwise_Dist"] <= 0.05:
            return False
        return row["Node"] in antichain_selected["Node"].values

    classified_purity["selected"] = classified_purity.apply(selected_condition, axis=1)
    classified_purity.sort_values(
        by=["selected", "Purity", "Min_Pairwise_Dist"], ascending=[False, False, False], inplace=True
    )
    internal_nodes = classified_purity[classified_nodes > 1]
    classified_purity["Z_min_dist"] = (
        (internal_nodes["Min_Pairwise_Dist"] - internal_nodes["Min_Pairwise_Dist"].mean())
        / internal_nodes["Min_Pairwise_Dist"].std()
        if internal_nodes["Min_Pairwise_Dist"].std() > 0
        else 0
    )

    classified_purity = classified_purity.reset_index(drop=True)

    return classified_purity


def find_assembly_mapping(row, stats_matrix):
    accession = row["assembly_accession"]
    if accession is None or pd.isna(accession):
        row["clade"] = "unmapped"
        row["nuniq"] = 0
        row["freq"] = 0
        row["Min_Pairwise_Dist"] = 0
        row["nleaves"] = 0
        return row

    match = stats_matrix[
        stats_matrix["leaves"].str.contains(accession, na=False)
        | (stats_matrix["clade"].str.contains(accession, na=False))
    ]

    if match.empty:
        row["clade"] = None
        row["nuniq"] = 0
        row["freq"] = 0
        row["Min_Pairwise_Dist"] = 0
        row["nleaves"] = 0
    else:
        row["clade"] = match["clade"].values[0]
        row["nuniq"] = match["nuniq"].values[0]
        row["freq"] = match["freq"].values[0]
        row["Min_Pairwise_Dist"] = match["Min_Pairwise_Dist"].values[0]
        row["nleaves"] = match["nleaves"].values[0]

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


def update_df_best_match(row, taxid_list, ncbi_wrapper: NCBITaxonomistWrapper, min_score=0.7, best_level=None):
    taxid = row["taxid"]
    if pd.isna(taxid):
        row["best_match_taxid"] = None
        row["best_match_level"] = -1
        row["best_match_score"] = 0.0
        row["name"] = None

    best_taxid, best_level, best_score = find_best_match(taxid, taxid_list, ncbi_wrapper)

    if best_score < min_score:
        best_taxid = None
        best_level = -1
        best_score = 0.0

    row["best_match_taxid"] = best_taxid
    row["best_match_level"] = best_level
    row["best_match_score"] = best_score

    if best_taxid is not None:
        row["name"] = ncbi_wrapper.get_name(best_taxid)
    else:
        row["name"] = None

    return row


def match_leaf(accid, list):
    for item in list:
        if accid in item:
            return item
    return None
