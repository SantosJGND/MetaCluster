import logging
import os
import re

import networkx as nx
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def merge_by_assembly_ID(m_stats_matrix: pd.DataFrame) -> pd.DataFrame:
    merged_rows = []
    for assembly_id, group in m_stats_matrix.groupby("assid"):
        group["total_uniq_reads"] = group["total_uniq_reads"].sum()
        group["coverage"] = group["coverage"].mean()
        group["meanmapq"] = group["meanmapq"].mean()
        group["covbases"] = group["covbases"].sum()
        group["error_rate"] = group["error_rate"].mean()
        merged_rows.append(group.sort_values(by="coverage", ascending=False).iloc[0])
    return pd.DataFrame(merged_rows)


def merge_to_matched(row, matched):
    if matched.empty:
        row["taxid"] = None
        row["assid"] = None
        row["description"] = None
        row["total_uniq_reads"] = 0
        return row

    file = row["file"]
    file_basename = os.path.basename(file)
    matched_assids = matched["assembly_accession"].tolist()
    match_assid = [assid for assid in matched_assids if assid in file_basename]
    if match_assid:
        row["taxid"] = matched.loc[matched["assembly_accession"] == match_assid[0], "taxid"].values[0]
        row["assid"] = match_assid[0]
        row["description"] = matched.loc[matched["assembly_accession"] == match_assid[0], "description"].values[0]
        row["total_uniq_reads"] = matched.loc[
            matched["assembly_accession"] == match_assid[0], "total_uniq_reads"
        ].values[0]
    else:
        row["taxid"] = None
        row["assid"] = None
        row["description"] = None
        row["total_uniq_reads"] = 0
    return row


NCBI_ACCESSION_RE = re.compile(r'(?:GC[AF]|NC|NZ|NT|NW|AC|AE|CP)_\d+\.\d+')

def _extract_assid(index_value: str, known_assids: list[str] | None = None) -> str | None:
    """Extract the assembly accession from a distance-matrix index.

    Tries in order:
    1. Match a standard NCBI assembly accession pattern (GCF_, GCA_, NC_, …).
    2. If a list of known assids is given, build a regex union and match
       any of them (longest match wins, same strategy as _merge_matched_vectorized).

    Returns None when nothing matches.
    """
    m = NCBI_ACCESSION_RE.search(index_value)
    if m:
        return m.group(0)
    if known_assids:
        known_sorted = sorted(set(known_assids), key=len, reverse=True)
        known_pat = re.compile("(" + "|".join(re.escape(a) for a in known_sorted) + ")")
        m = known_pat.search(index_value)
        if m:
            return m.group(0)
    return None


def _merge_matched_vectorized(m_stats: pd.DataFrame, matched: pd.DataFrame) -> pd.DataFrame:
 
    if matched.empty:
        result = m_stats.copy()
        result["taxid"] = None
        result["assid"] = None
        result["description"] = None
        result["total_uniq_reads"] = 0
        return result

    known_assids = matched["assembly_accession"].unique().tolist()
    file_basenames = m_stats["file"].apply(os.path.basename)
    matched_acc = file_basenames.apply(lambda x: _extract_assid(x, known_assids))
    matched_acc = matched_acc.where(matched_acc.isin(known_assids), None)

    result = m_stats.copy()
    result["_matched_acc"] = matched_acc

    result = result.merge(
        matched[["assembly_accession", "taxid", "description", "total_uniq_reads"]],
        left_on="_matched_acc",
        right_on="assembly_accession",
        how="left",
        suffixes=("", "_r"),
    )

    result["assid"] = result["_matched_acc"].where(result["_matched_acc"].notna(), None)
    result["taxid"] = result["taxid"].astype(object).where(result["taxid"].notna(), None)
    result["description"] = result["description"].astype(object).where(result["description"].notna(), None)
    result["total_uniq_reads"] = result["total_uniq_reads"].fillna(0).astype(int)

    result.drop(columns=["_matched_acc"], inplace=True, errors="ignore")

    return result


class OverlapManager:
    """
    Manages overlap-based clustering analysis for metagenomic data.

    Handles tree construction, node pruning, and statistical analysis of
    overlapping genomic regions across multiple samples.
    """

    def __init__(
        self, output_dir: str, max_proportion: float = 1.0, max_taxids: int | None = None, skip_build: bool = False
    ):
        self.output_dir = output_dir
        self.data_set_name = os.path.basename(os.path.dirname(output_dir))
        self.distance_matrix_filepath = os.path.join(output_dir, "distance_matrix.tsv")
        self.rundir = os.path.dirname(output_dir)
        self.max_proportion = max_proportion
        self.max_taxids = max_taxids
        self.leaves = []
        self.all_leaves_global = []
        self.taxids_to_keep = []
        self.cluster_map = {}
        self.all_node_stats = pd.DataFrame()
        self.all_edges = pd.DataFrame()
        self.all_nodes = []
        self.tree = nx.DiGraph()
        self.root_nodes = []
        self.distance_mat = pd.DataFrame()
        self.original_m_stats_matrix = pd.DataFrame()
        self.m_stats_matrix = pd.DataFrame()
        self.node_leaves_cache: dict[str, list[str]] = {}

        if not skip_build and self.check_data_available():
            self.build()

    def build(self):
        """Build tree from distance matrix and compute node statistics."""
        self.m_stats_matrix = pd.DataFrame()
        try:
            self.new_tree_from_distance_matrix()
            self.prune_empty_nodes()
            self._rebuild_node_leaves_cache()
            self.recalculate_all_min_pairwise_dist()
            print("recalc")
        except Exception as e:
            import traceback

            logger.error(f"Error creating new tree from distance matrix: {e}")
            traceback.print_exc()
            raise e

        self.root_nodes = [n for n, d in self.tree.in_degree() if d == 0]
        self.leaves = [n for n, d in self.tree.out_degree() if d == 0]

    def get_leaf(self, row):
        accession = row["assid"]
        if accession is None or pd.isna(accession):
            return None

        for leaf in self.all_leaves_global:
            if leaf == accession or leaf.endswith(accession) or accession in leaf:
                return leaf
        return None

    @staticmethod
    def matrix_to_phylotriangle(distance_matrix: pd.DataFrame):
        """
        Convert a square distance DataFrame to a Bio.Phylo DistanceMatrix

        (lower-triangular format required by Bio.Phylo.TreeConstruction).
        """
        distmat = distance_matrix.values.tolist()
        distmat = [x[: i + 1] for i, x in enumerate(distmat)]
        distmat = DistanceMatrix(list(distance_matrix.index), distmat)

        return distmat

    def tree_from_distance_matrix(self, distance_matrix: pd.DataFrame):
        """
        Build a Bio.Phylo Tree from a distance matrix.

        Uses NJ for >=3 taxa, UPGMA for 2 taxa, and returns a singleton
        tree for a single taxon.
        """
        distmat = self.matrix_to_phylotriangle(distance_matrix)
        constructor = DistanceTreeConstructor()
        n = distance_matrix.shape[0] if not distance_matrix.empty else 0

        if n == 0:
            return Phylo.BaseTree.Tree(rooted=True)

        if n == 1:
            tree = Phylo.BaseTree.Tree(rooted=False)
            tree.clade = Phylo.BaseTree.Clade(name=distance_matrix.index[0])
            tree.ladderize()
            return tree

        if n == 2:
            tree = constructor.upgma(distmat)
        else:
            tree = constructor.nj(distmat)

        tree.rooted = False
        tree.ladderize()
        return tree

    def bio_tree_edges_to_dataframe(self, tree) -> pd.DataFrame:
        edges = []
        for clade in tree.find_clades(order="level"):
            for child in clade.clades:
                edges.append((clade.name, child.name))
        return pd.DataFrame(edges, columns=["node1", "node2"])

    def njbio_from_distance_matrix(self, distance_matrix: pd.DataFrame):
        """Build NetworkX DiGraph from a Bio.Phylo distance-based tree."""
        tree = self.tree_from_distance_matrix(distance_matrix)

        self.all_edges = self.bio_tree_edges_to_dataframe(tree)

        self.all_nodes = list(set(self.all_edges["node1"]).union(set(self.all_edges["node2"])))
        self.tree = nx.from_pandas_edgelist(self.all_edges, source="node1", target="node2", create_using=nx.DiGraph())
        self.root_nodes = [n for n, d in self.tree.in_degree() if d == 0]
        self.leaves = [n.name for n in tree.get_terminals()]

    def read_distance_matrix(self) -> pd.DataFrame:
        if os.path.exists(self.distance_matrix_filepath):
            distance_matrix = pd.read_csv(self.distance_matrix_filepath, sep="\t", index_col=0)
            self.all_leaves_global = list(distance_matrix.index)

            distance_matrix = distance_matrix.loc[~distance_matrix.index.duplicated(keep="first")]
            distance_matrix = distance_matrix.loc[:, ~distance_matrix.columns.duplicated(keep="first")]
            nodes_to_filter = self.m_stats_matrix[self.m_stats_matrix["coverage"] > 0]["assid"].tolist()
            all_nodes_in_matrix = distance_matrix.index.to_list()
            nodes_to_keep = [node for node in all_nodes_in_matrix if any(n in str(node) for n in nodes_to_filter)]

            distance_matrix = distance_matrix.loc[
                distance_matrix.index.isin(nodes_to_keep), distance_matrix.columns.isin(nodes_to_keep)
            ]

            if distance_matrix.shape[0] < 1:
                return pd.DataFrame()

            distance_matrix.index = distance_matrix.index.map(str)
            distance_matrix.columns = distance_matrix.columns.map(str)
            return distance_matrix
        
        return pd.DataFrame()

    def new_tree_from_distance_matrix(self):
        """
        Construct a new tree from the distance matrix.

        The distance matrix is asymmetric where dist[i][j] = 1 - (proportion of i reads also present in j)
        """
        if os.path.exists(self.distance_matrix_filepath):
            merged_stats_file = os.path.join(
                os.path.dirname(self.output_dir), "output", "merged_coverage_statistics.tsv"
            )
            matched_assemblies_file = os.path.join(self.rundir, "output", "matched_assemblies.tsv")
            merged_classification_results_filepath = os.path.join(
                self.rundir, "classification", f"{os.path.basename(self.rundir)}_merged_classification.tsv"
            )
            from metagenomics_utils.overlap_manager.node_stats import _resolve_assembly_classifications

            merged_classification_results = pd.read_csv(merged_classification_results_filepath, sep="\t")

            matched = pd.read_csv(matched_assemblies_file, sep="\t")

            matched = _resolve_assembly_classifications(matched, merged_classification_results)
            if os.path.exists(merged_stats_file):
                m = pd.read_csv(merged_stats_file, sep="\t").rename(columns={"#rname": "assembly_accession"})
                if "error_rate" in m.columns:
                    m = m[["assembly_accession", "coverage", "covbases", "meanmapq", "numreads", "error_rate", "file"]]
                else:
                    m = m[["assembly_accession", "coverage", "covbases", "meanmapq", "numreads", "file"]]
                    m["error_rate"] = 0.0
                m = m.groupby("file").agg(
                    {
                        "coverage": "max",
                        "covbases": "sum",
                        "meanmapq": "max",
                        "numreads": "sum",
                        "error_rate": "max",
                        "assembly_accession": "first",
                    }
                ).reset_index()
                self.original_m_stats_matrix = m
            else:
                self.original_m_stats_matrix = pd.DataFrame()
            self.m_stats_matrix = self.original_m_stats_matrix.copy()

            if not self.m_stats_matrix.empty and "total_uniq_reads" in matched.columns:
                self.m_stats_matrix = _merge_matched_vectorized(self.m_stats_matrix, matched)

            self.m_stats_matrix = merge_by_assembly_ID(self.m_stats_matrix)
            if not self.m_stats_matrix.empty and "total_uniq_reads" in matched.columns:
                self.m_stats_matrix = self.m_stats_matrix.sort_values(by=["total_uniq_reads"], ascending=False)
                index_keep = int(self.m_stats_matrix.shape[0] * self.max_proportion)
                if self.max_taxids is not None:
                    index_keep = min(self.m_stats_matrix.shape[0], self.max_taxids)
                self.m_stats_matrix = self.m_stats_matrix[:index_keep].reset_index(drop=True)

            distance_matrix = self.read_distance_matrix()

            self.m_stats_matrix["leaf"] = self.m_stats_matrix.apply(self.get_leaf, axis=1)

            self.m_stats_matrix = self.m_stats_matrix.dropna(subset=["leaf"])
            self.m_stats_matrix.set_index("leaf", inplace=True)

            if distance_matrix.empty is True:
                self.tree = nx.DiGraph()
                self.all_nodes = list(self.tree.nodes())
                self.all_node_stats = pd.DataFrame(
                    columns=["Node", "Num_Leaves", "Min_Dist", "Private_Reads", "Private_Proportion"]
                )
                return
            if distance_matrix.shape[0] == 1:
                self.distance_mat = distance_matrix
                self.leaves = list(distance_matrix.index)
                self.all_nodes = list(distance_matrix.index)
                self.root_nodes = list(distance_matrix.index)
                self.tree = nx.DiGraph()
                self.tree.add_node(self.leaves[0])
                self.recreate_all_node_stats_from_tree()
                return

            proximity_matrix = 1 - distance_matrix
            weighted_proximity_matrix = self.weighted_matrix(proximity_matrix)
            weighted_distance_matrix = 1 - weighted_proximity_matrix

            weighted_distance_matrix.to_csv(
                os.path.join(self.output_dir, "weighted_distance_matrix.tsv"), sep="\t", index=True
            )
            # standardize weights by row
            symweighted_distance_matrix = self.matrix_symmetrize_mean(weighted_distance_matrix)
            self.distance_mat = symweighted_distance_matrix

            self.njbio_from_distance_matrix(symweighted_distance_matrix)
            self.all_nodes = list(self.tree.nodes())
            self.leaves = [n for n, d in self.tree.out_degree() if d == 0]
            self.root_nodes = [n for n, d in self.tree.in_degree() if d == 0]

            self.recreate_all_node_stats_from_tree()

            logger.info(f"Tree created with {len(self.tree.nodes())} nodes and {len(self.leaves)} leaves.")

            self.all_edges = pd.DataFrame(self.tree.edges(), columns=["node1", "node2"])

    def recreate_all_node_stats_from_tree(self):
        self._rebuild_node_leaves_cache()

        new_stats = []
        for node in self.all_nodes:
            leaves = self.node_leaves_cache.get(node, [])

            if not leaves:
                continue

            distance_matrix_subset = self.distance_mat.loc[leaves, leaves]
            min_dist = distance_matrix_subset.min().min()
            new_stats.append(
                {
                    "Node": node,
                    "Num_Leaves": len(leaves),
                    "Min_Dist": min_dist,
                    "Private_Reads": 0,
                    "Private_Proportion": 1,
                }
            )

        self.all_node_stats = pd.DataFrame(new_stats)
        if self.all_node_stats.empty is True:
            self.all_node_stats = pd.DataFrame(
                columns=["Node", "Num_Leaves", "Min_Dist", "Private_Reads", "Private_Proportion"]
            )

    def weighted_matrix(self, distance_matrix: pd.DataFrame) -> pd.DataFrame:
        """
        Weight the distance matrix by the relative error rates of mappings against respective assemblies. 

        """
        known_assids = self.m_stats_matrix["assid"].dropna().unique().tolist() if not self.m_stats_matrix.empty else None
        valid_idx = []
        err_vals = []
        for orig_idx in distance_matrix.index:
            # Legacy parsing — kept for reference:
            # assid = orig_idx.strip("./").strip(self.data_set_name).strip(".fna.bam")
            # assid = "_".join(assid.split("_")[1:])
            assid = _extract_assid(orig_idx, known_assids)
            if assid is None:
                continue
            match = self.m_stats_matrix[self.m_stats_matrix["assid"] == assid]
            if not match.empty:
                valid_idx.append(orig_idx)
                err_vals.append(match["error_rate"].values[0])

        distance_matrix = distance_matrix.loc[valid_idx, valid_idx]
        values = distance_matrix.values
        index = distance_matrix.index
        cols = distance_matrix.columns
        n = distance_matrix.shape[0]

        if n == 0:
            return pd.DataFrame()

        err = np.array(err_vals)
        ratio = err[np.newaxis, :] / np.maximum(err[:, np.newaxis], 1e-10)
        ratio = np.where(ratio > 1.0, 1.0 / ratio, ratio)

        weighted = values * ratio

        mask = weighted >= 1.0
        weighted[mask] = 1.0 / weighted[mask]

        np.fill_diagonal(weighted, 1.0)

        return pd.DataFrame(weighted, index=index, columns=cols)

    def matrix_symmetrize_mean(self, distance_matrix: pd.DataFrame) -> pd.DataFrame:
        vals = distance_matrix.values
        sym = np.maximum(vals, vals.T)
        return pd.DataFrame(sym, index=distance_matrix.index, columns=distance_matrix.columns)

    def recalculate_all_min_pairwise_dist(self):
        if not os.path.exists(self.distance_matrix_filepath):
            return

        distance_matrix = self.read_distance_matrix()
        proximity_matrix = 1 - distance_matrix
        weighted_proximity_matrix = self.weighted_matrix(proximity_matrix)
        wp_values = weighted_proximity_matrix.values
        index_arr = list(weighted_proximity_matrix.index)
        index_to_pos = {name: i for i, name in enumerate(index_arr)}

        def calc_node_stats(node):
            leaves_parted = self.get_leaves_parted(node)
            if len(leaves_parted) < 2:
                return {"Min_Pairwise_Dist": 0, "Min_Shared": 0, "Min_Dist": 0}

            left = leaves_parted[0]
            right = leaves_parted[1]
            left_pos = [index_to_pos[l] for l in left]
            right_pos = [index_to_pos[l] for l in right]
            all_pos = left_pos + right_pos

            sub = wp_values[np.ix_(left_pos, right_pos)]
            pair_dists = sub #np.maximum(sub, sub.T)
            min_dist = float(pair_dists.min()) if pair_dists.size > 0 else 0.0

            all_internal = set(left) | set(right)
            other_keys = [l for l in index_arr if l not in all_internal]
            if other_keys:
                other_pos = [index_to_pos[l] for l in other_keys]
                sub_ab = wp_values[np.ix_(all_pos, other_pos)]
                ext_kept = float(sub_ab.min())
            else:
                ext_kept = 0.0

            return {
                "Min_Dist": min_dist,
                "Min_Shared": ext_kept,
                "Min_Pairwise_Dist": min_dist if min_dist < 1 else 1.0 / max(min_dist, 1e-10),
            }

        self.all_node_stats = self.all_node_stats[self.all_node_stats["Node"].isin(self.all_nodes)]
        node_stats_list = []
        for node in self.all_node_stats["Node"]:
            stats = calc_node_stats(node)
            stats["Node"] = node
            node_stats_list.append(stats)

        if not node_stats_list:
            self.all_node_stats["Min_Pairwise_Dist"] = 0
            self.all_node_stats["Min_Shared"] = 0
            self.all_node_stats["Min_Dist"] = 0
        else:
            stats_df = pd.DataFrame(node_stats_list)
            self.all_node_stats = self.all_node_stats.merge(stats_df, on="Node", how="left")

        internal_nodes = self.all_node_stats[~self.all_node_stats["Node"].isin(self.leaves)]
        if internal_nodes.empty is True:
            self.all_node_stats["Z_min_dist"] = 0
            self.all_node_stats["Z_Min_Shared"] = 0
        else:
            self.all_node_stats["Z_min_dist"] = (
                self.all_node_stats["Min_Pairwise_Dist"] - internal_nodes["Min_Pairwise_Dist"].mean()
            ) / (internal_nodes["Min_Pairwise_Dist"].std() if internal_nodes["Min_Pairwise_Dist"].std() > 0 else 1)
            self.all_node_stats["Z_Min_Shared"] = (
                self.all_node_stats["Min_Shared"] - internal_nodes["Min_Shared"].mean()
            ) / (internal_nodes["Min_Shared"].std() if internal_nodes["Min_Shared"].std() > 0 else 1)

    def determine_min_pairwise_dist(self, node: str) -> float:
        if node not in self.all_nodes:
            return float("nan")
        if self.distance_mat.empty:
            return float("nan")
        leaves = self.get_node_leaves(node)
        if len(leaves) < 2:
            return 0.0
        sub_matrix = self.distance_mat.loc[leaves, leaves]
        tril_indices = np.tril_indices_from(sub_matrix, k=-1)
        pairwise_distances = sub_matrix.values[tril_indices]
        return float(pairwise_distances.min())

    def remove_leaf_and_ascendants(self, node: str):
        """
        Remove node and all its ascendants that become childless.
        State (all_nodes, root_nodes, leaves) must be recomputed by the caller after all removals.
        """
        if node not in self.tree:
            return

        parents = list(self.tree.predecessors(node))
        self.tree.remove_node(node)
        for parent in parents:
            if parent in self.tree and self.tree.out_degree(parent) == 0:
                self.remove_leaf_and_ascendants(parent)

    def prune_empty_nodes(self):
        """Remove nodes from the tree that have no coverage."""
        self.all_node_stats = self.all_node_stats.dropna()
        nodes_to_remove = []

        existing_leaves = self.m_stats_matrix[self.m_stats_matrix["coverage"] > 0.0].index.tolist()
        nodes_to_remove = nodes_to_remove + [node for node in self.leaves if node not in existing_leaves]
        for node in nodes_to_remove:
            self.remove_leaf_and_ascendants(node)

        self.all_nodes = list(self.tree.nodes())
        self.all_node_stats = self.all_node_stats[self.all_node_stats["Node"].isin(self.all_nodes)]
        self.all_edges = pd.DataFrame(self.tree.edges(), columns=["node1", "node2"])
        if self.all_edges.empty:
            self.tree = nx.DiGraph()
            self.tree.add_nodes_from(self.all_nodes)
        else:
            self.tree = nx.from_pandas_edgelist(
                self.all_edges, source="node1", target="node2", create_using=nx.DiGraph()
            )

        self.leaves = [n for n, d in self.tree.out_degree() if d == 0]
        self.root_nodes = [n for n, d in self.tree.in_degree() if d == 0]

    def print_tree(self, threshold=0.5):
        if self.tree.number_of_nodes() == 0:
            return None
        try:
            import matplotlib.pyplot as plt

            pos = nx.nx_agraph.graphviz_layout(self.tree, prog="dot")
            fig, ax = plt.subplots(figsize=(12, 8))
            nx.draw(
                self.tree,
                pos,
                with_labels=True,
                node_size=500,
                node_color="lightblue",
                font_size=10,
                font_weight="bold",
                ax=ax,
            )
            plt.close(fig)
            return fig
        except ImportError:
            print("matplotlib is not installed. Please install it to visualize the tree.")
            return None
        except Exception as e:
            print(f"An error occurred while generating the tree visualization: {e}")
            return None

    def print_tree_given_colors(self, nodes: list):
        if self.tree.number_of_nodes() == 0:
            print("The tree is empty.")
            return None
        try:
            import matplotlib.pyplot as plt

            color_list = plt.get_cmap("tab20").colors

            color_map = {node: color_list[i % len(color_list)] for i, node in enumerate(nodes)}
            for node in nodes:
                descendants = nx.descendants(self.tree, node)
                for desc in descendants:
                    color_map[desc] = color_map[node]

            node_colors = [color_map.get(node, "lightgrey") for node in self.tree.nodes()]

            pos = nx.nx_agraph.graphviz_layout(self.tree, prog="dot")
            fig, ax = plt.subplots(figsize=(12, 8))
            nx.draw(
                self.tree,
                pos,
                with_labels=True,
                node_size=500,
                node_color=node_colors,
                font_size=7,
                font_weight="bold",
                ax=ax,
            )
            plt.close(fig)
            return fig
        except ImportError:
            print("matplotlib is not installed. Please install it to visualize the tree.")
            return None
        except Exception as e:
            print(f"An error occurred while generating the tree visualization: {e}")
            return None

    def print_tree_given_colors_plotly(self, nodes: list):
        if self.tree.number_of_nodes() == 0:
            print("The tree is empty.")
            return None
        try:
            import matplotlib.pyplot as plt
            import plotly.graph_objects as go

            color_list = plt.get_cmap("tab20").colors
            all_nodes_stats = self.all_node_stats.set_index("Node").to_dict("index")
            color_map = {
                node: f"rgb({int(color_list[i][0] * 255)},{int(color_list[i][1] * 255)},{int(color_list[i][2] * 255)})"
                for i, node in enumerate(nodes)
            }
            for node in nodes:
                descendants = nx.descendants(self.tree, node)
                for desc in descendants:
                    color_map[desc] = color_map[node]

            node_colors = [color_map.get(node, "lightgrey") for node in self.tree.nodes()]

            pos = nx.nx_agraph.graphviz_layout(self.tree, prog="dot")

            edge_x = []
            edge_y = []
            for edge in self.tree.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.append(x0)
                edge_x.append(x1)
                edge_x.append(None)
                edge_y.append(y0)
                edge_y.append(y1)
                edge_y.append(None)

            edge_trace = go.Scatter(
                x=edge_x, y=edge_y, line=dict(width=1, color="#888"), hoverinfo="none", mode="lines"
            )

            node_x = []
            node_y = []
            for node in self.tree.nodes():
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)

            node_trace = go.Scatter(
                x=node_x,
                y=node_y,
                mode="markers+text",
                text=[
                    f"{node}\n{all_nodes_stats[node]['Min_Pairwise_Dist']:.2f}" if node in all_nodes_stats else node
                    for node in self.tree.nodes()
                ],
                textposition="bottom center",
                hoverinfo="text",
                marker=dict(showscale=False, colorscale="YlGnBu", color=node_colors, size=20, line_width=2),
            )
            fig = go.Figure(
                data=[edge_trace, node_trace],
                layout=go.Layout(
                    title="Tree Visualization",
                    showlegend=False,
                    hovermode="closest",
                    margin=dict(l=0, r=0, t=40, b=0),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                ),
            )

            return fig
        except ImportError:
            print("Plotly is not installed. Please install it to visualize the tree.")
            return None
        except Exception as e:
            print(f"An error occurred while generating the tree visualization: {e}")
            return None

    def print_tree_with_clades(self, threshold=0.5):
        """Draw tree with clades colored by cluster membership.

        Side effect: populates ``self.cluster_map`` via ``selected_nodes_matrix()``.
        """
        if self.tree.number_of_nodes() == 0:
            print("The tree is empty.")
            return None
        try:
            import matplotlib.pyplot as plt

            color_list = plt.get_cmap("tab20").colors

            node_colors = []
            all_nodes_stats = self.all_node_stats.set_index("Node")
            _ = self.selected_nodes_matrix(threshold)

            color_map = {}

            color_index = 0

            for node in self.tree.nodes():
                cluster = self.cluster_map.get(node, None)
                if cluster:
                    if cluster not in color_map:
                        color_map[cluster] = color_list[color_index % len(color_list)]
                        color_index += 1
                    node_colors.append(color_map[cluster])
                else:
                    node_colors.append("lightgrey")

            pos = nx.nx_agraph.graphviz_layout(self.tree, prog="dot")
            fig, ax = plt.subplots(figsize=(12, 8))
            labels_with_min_dist = {
                node: f"{node}\n{all_nodes_stats.loc[node]['Min_Pairwise_Dist']:.2f}" for node in self.tree.nodes()
            }
            nx.draw(
                self.tree,
                pos,
                with_labels=True,
                node_size=500,
                node_color=node_colors,
                font_size=7,
                font_weight="bold",
                ax=ax,
                labels=labels_with_min_dist,
            )
            plt.close(fig)
            return fig
        except ImportError:
            print("matplotlib is not installed. Please install it to visualize the tree.")
            return None
        except Exception as e:
            print(f"An error occurred while generating the tree visualization: {e}")
            return None

    def print_tree_with_clades_plotly(self, threshold=0.5):
        if self.tree.number_of_nodes() == 0:
            print("The tree is empty.")
            return None
        try:
            import matplotlib.pyplot as plt
            import plotly.graph_objects as go

            color_list = plt.get_cmap("tab20").colors

            node_colors = []
            all_nodes_stats = self.all_node_stats.set_index("Node").to_dict("index")
            selected_nodes_matrix = self.selected_nodes_matrix(threshold)

            color_map = {}

            color_index = 0

            for node in self.tree.nodes():
                cluster = self.cluster_map.get(node, None)
                if cluster:
                    if cluster not in color_map:
                        color_map[cluster] = (
                            f"rgb({int(color_list[color_index][0] * 255)},{int(color_list[color_index][1] * 255)},{int(color_list[color_index][2] * 255)})"
                        )
                        color_index += 1
                    node_colors.append(color_map[cluster])
                else:
                    node_colors.append("lightgrey")

            pos = nx.nx_agraph.graphviz_layout(self.tree, prog="dot")

            edge_x = []
            edge_y = []
            for edge in self.tree.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.append(x0)
                edge_x.append(x1)
                edge_x.append(None)
                edge_y.append(y0)
                edge_y.append(y1)
                edge_y.append(None)

            edge_trace = go.Scatter(
                x=edge_x, y=edge_y, line=dict(width=1, color="#888"), hoverinfo="none", mode="lines"
            )

            node_x = []
            node_y = []
            for node in self.tree.nodes():
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)

            node_trace = go.Scatter(
                x=node_x,
                y=node_y,
                mode="markers+text",
                text=[
                    f"{node}\n{all_nodes_stats[node]['Min_Pairwise_Dist']:.2f}" if node in all_nodes_stats else node
                    for node in self.tree.nodes()
                ],
                textposition="bottom center",
                hoverinfo="text",
                marker=dict(showscale=False, colorscale="YlGnBu", color=node_colors, size=20, line_width=2),
            )
            fig = go.Figure(
                data=[edge_trace, node_trace],
                layout=go.Layout(
                    showlegend=False,
                    hovermode="closest",
                    margin=dict(b=20, l=5, r=5, t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                ),
            )
            return fig
        except ImportError:
            print("plotly is not installed. Please install it to visualize the tree.")
            return None
        except Exception as e:
            print(f"An error occurred while generating the tree visualization: {e}")
            return None

    def check_data_available(self):
        if os.path.exists(self.distance_matrix_filepath) is False:
            return False
        return True

    def node_selector(self, node: str, threshold=0.5):
        if node in self.leaves:
            return True
        node_data = self.all_node_stats[self.all_node_stats["Node"] == node]
        if node_data.empty:
            return False
        return node_data["Min_Pairwise_Dist"].values[0] >= threshold and node_data["Min_Shared"].values[0] <= threshold

    def traverse_graph_recursive(self, node, threshold=0.5, nodes_return=None, _seen=None) -> list[str]:
        if nodes_return is None:
            nodes_return = []
            _seen = set()
        if self.node_selector(node, threshold):
            if node not in _seen:
                nodes_return.append(node)
                _seen.add(node)
                descendents = nx.descendants(self.tree, node)
                for desc in descendents:
                    self.cluster_map[desc] = node
                self.cluster_map[node] = node
        else:
            for neighbor in self.tree.neighbors(node):
                if neighbor not in _seen:
                    self.traverse_graph_recursive(neighbor, threshold, nodes_return, _seen)

        return nodes_return

    def get_node_leaves(self, node):
        return self.node_leaves_cache.get(node, [])

    def get_leaves_parted(self, node):
        if self.tree.out_degree(node) == 0:
            return [self.node_leaves_cache.get(node, [node])]
        return [self.node_leaves_cache.get(child, []) for child in self.tree.successors(node)]

    def _rebuild_node_leaves_cache(self):
        self.node_leaves_cache = {}
        def _dfs(node):
            if self.tree.out_degree(node) == 0:
                leaves = [node]
            else:
                leaves = []
                for child in self.tree.successors(node):
                    leaves.extend(_dfs(child))
            self.node_leaves_cache[node] = leaves
            return leaves
        for root in self.root_nodes:
            _dfs(root)

    def selected_nodes_matrix(self, threshold=0.5):
        self.cluster_map = {}

        selected_nodes = []
        for root in self.root_nodes:
            selected_nodes.extend(self.traverse_graph_recursive(root, threshold, []))
        selected_nodes = list(set(selected_nodes))

        stats_matrix = self.all_node_stats[self.all_node_stats["Node"].isin(selected_nodes)]

        def get_leaves_str(row):
            node = row["Node"]
            leaves = self.get_node_leaves(node)
            row["leaves"] = ",".join(map(str, leaves))
            return row

        stats_matrix = stats_matrix.apply(get_leaves_str, axis=1)
        return stats_matrix
