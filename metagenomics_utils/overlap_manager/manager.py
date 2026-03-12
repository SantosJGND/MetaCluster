import pandas as pd
import networkx as nx
import os
from typing import List
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, to_tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo


def merge_by_assembly_ID(m_stats_matrix: pd.DataFrame) -> pd.DataFrame:
    merged_rows = []
    for assembly_id, group in m_stats_matrix.groupby('assid'):
        group['total_uniq_reads'] = group['total_uniq_reads'].sum()
        group['coverage'] = group['coverage'].mean()
        group['meanmapq'] = group['meanmapq'].mean()
        group['covbases'] = group['covbases'].sum()
        group['error_rate'] = group['error_rate'].mean()
        merged_rows.append(group.sort_values(by='coverage', ascending=False).iloc[0])
    return pd.DataFrame(merged_rows)


def merge_to_matched(row, matched):
    if matched.empty:
        row['taxid'] = None
        row['assid'] = None
        row['description'] = None
        row['total_uniq_reads'] = 0
        return row
    
    file = row['file']
    matched_assids = matched['assembly_accession'].tolist()
    match_assid = [assid for assid in matched_assids if assid in file]
    if match_assid:
        row['taxid'] = matched.loc[matched['assembly_accession'] == match_assid[0], 'taxid'].values[0]
        row['assid'] = match_assid[0]
        row['description'] = matched.loc[matched['assembly_accession'] == match_assid[0], 'description'].values[0]
        row['total_uniq_reads'] = matched.loc[matched['assembly_accession'] == match_assid[0], 'total_uniq_reads'].values[0]
    else:
        row['taxid'] = None
        row['assid'] = None
        row['description'] = None
        row['total_uniq_reads'] = 0
    return row


class OverlapManager:
    def __init__(self, output_dir: str, max_taxids = 20):
        self.output_dir = output_dir
        self.rundir = os.path.dirname(output_dir)
        self.max_taxids = max_taxids
        self.leaves = []
        if not self.check_data_available():
            self.all_node_stats = pd.DataFrame()
            self.all_edges = pd.DataFrame()
            self.all_nodes = []
            self.tree = nx.DiGraph()
            self.root_nodes = []
            self.distance_mat = pd.DataFrame()
            self.distance_matrix_filepath = ""
            self.m_stats_matrix = pd.DataFrame()
            return
        
        self.cluster_map = {}
        self.all_node_stats = pd.read_csv(os.path.join(output_dir, "all_node_statistics.tsv"), sep = "\t")
        matched_assemblies_file = os.path.join(self.rundir, "output", "matched_assemblies.tsv")
        self.distance_matrix_filepath = os.path.join(output_dir, "distance_matrix.tsv")
        self.m_stats_matrix = pd.DataFrame()

        try:
            try:
                self.all_edges = pd.read_csv(os.path.join(output_dir, "nj_tree_edges.txt"), sep = "\t", header=None)
            except pd.errors.EmptyDataError:
                self.all_edges = pd.DataFrame(columns=['node1', 'node2'])

                
            
            self.all_edges.columns = ['node1', 'node2']
            self.all_nodes = list(self.all_node_stats['Node'].unique())
            if self.all_edges.empty:
                self.tree= nx.DiGraph() 
                self.tree.add_nodes_from(self.all_nodes)
            else:
                self.tree = nx.from_pandas_edgelist(self.all_edges, source='node1', target='node2',create_using=nx.DiGraph())

            distance_matrix = pd.read_csv(self.distance_matrix_filepath, sep="\t", index_col=0)
            self.distance_mat = pd.read_csv(self.distance_matrix_filepath, sep="\t", index_col=0)
            merged_stats_file = os.path.join(os.path.dirname(self.output_dir), "output", "merged_coverage_statistics.tsv")
            matched = pd.read_csv(matched_assemblies_file, sep="\t")
            matched['assembly_file'] = matched['assembly_file'].apply(lambda x: os.path.basename(x))
            matched = matched.sort_values(by=['total_uniq_reads'], ascending=False).drop_duplicates(subset=['assembly_accession'], keep= 'first')
            self.m_stats_matrix = pd.read_csv(merged_stats_file, sep="\t").rename(columns={"#rname": "assembly_accession"}) if os.path.exists(merged_stats_file) else pd.DataFrame()


            if not self.m_stats_matrix.empty and 'total_uniq_reads' in matched.columns:
                self.m_stats_matrix = self.m_stats_matrix.apply(merge_to_matched, axis=1, matched=matched)

            self.m_stats_matrix = merge_by_assembly_ID(self.m_stats_matrix)

            if not self.m_stats_matrix.empty and 'total_uniq_reads' in matched.columns:
                self.m_stats_matrix = self.m_stats_matrix.sort_values(by = ['total_uniq_reads'], ascending=False)
                self.m_stats_matrix = self.m_stats_matrix[:self.max_taxids].reset_index(drop=True)

            for n in self.all_nodes:
                if n not in self.tree:
                    self.tree.add_node(n)
        except Exception as e:
            
            import traceback
            traceback.print_exc()
            self.all_edges = pd.DataFrame(columns=['node1', 'node2'])
            distance_matrix = pd.read_csv(self.distance_matrix_filepath, sep="\t", index_col=0)
            self.distance_mat = pd.read_csv(self.distance_matrix_filepath, sep="\t", index_col=0)
            self.all_nodes = [str(n) for n in distance_matrix.index]
            self.tree = nx.DiGraph()
            self.tree.add_nodes_from(self.all_nodes)
            self.m_stats_matrix = pd.DataFrame(columns = ['description','taxid', 'file', 'assid', 'coverage'])
        
        self.leaves = [n for n,d in self.tree.out_degree() if d==0]
        self.m_stats_matrix['leaf'] = self.m_stats_matrix.apply(self.get_leaf, axis=1)
        self.m_stats_matrix = self.m_stats_matrix.dropna(subset=['leaf'])
        self.m_stats_matrix.set_index('leaf', inplace=True)

        for ix in enumerate(self.leaves):
            if ix[1] not in self.m_stats_matrix.index:
                self.leaves.remove(ix[1])
                self.all_nodes.remove(ix[1])
        
        try:
            self.prune_empty_nodes()
            self.recalculate_all_min_pairwise_dist()

        except Exception as e:
            import traceback
            print(f"Error pruning empty nodes: {e}")
            traceback.print_exc()

        self.root_nodes = [n for n,d in self.tree.in_degree() if d==0]
        self.leaves = [n for n,d in self.tree.out_degree() if d==0]

        try:

            if len(self.m_stats_matrix) > 1:
                self.new_tree_from_distance_matrix()
                self.recalculate_all_min_pairwise_dist()
        except Exception as e:
            import traceback
            print(f"Error creating new tree from distance matrix: {e}")
            traceback.print_exc()
            raise e
    
    def get_leaf(self,row):
        accession = row['assid']
        if accession is None or pd.isna(accession):
            return None

        for leaf in self.leaves:
            if accession in leaf:
                return leaf
        return None
    
    @staticmethod
    def matrix_to_phylotriangle(distance_matrix: pd.DataFrame):
        """
        Return tree from matrix
        """
        distmat = distance_matrix.values.tolist()
        distmat = [x[: i + 1] for i, x in enumerate(distmat)]
        distmat = DistanceMatrix(list(distance_matrix.index), distmat)

        return distmat

    def tree_from_distance_matrix(self, distance_matrix: pd.DataFrame):
        """
        Return tree from distance matrix
        """
        distmat = self.matrix_to_phylotriangle(distance_matrix)
        constructor = DistanceTreeConstructor()

        if distance_matrix.empty is True:
            return Phylo.BaseTree.Tree(rooted="True")

        if distance_matrix.shape[0] < 2:
            tree = constructor.nj(distmat)
        else:
            tree = constructor.upgma(distmat)

        tree.rooted = False
        tree.ladderize()
        return tree
    
    def bio_tree_edges_to_dataframe(self, tree) -> pd.DataFrame:
        edges = []
        for clade in tree.find_clades(order='level'):
            for child in clade.clades:
                edges.append((clade.name, child.name))
        return pd.DataFrame(edges, columns=['node1', 'node2'])
    

    def njbio_from_distance_matrix(self, distance_matrix: pd.DataFrame):
        """
        """
        tree = self.tree_from_distance_matrix(distance_matrix)
        nodelist = tree.get_nonterminals() + tree.get_terminals()

        self.all_edges = self.bio_tree_edges_to_dataframe(tree)
        
        self.all_nodes = list(set(self.all_edges['node1']).union(set(self.all_edges['node2'])))
        self.tree = nx.from_pandas_edgelist(self.all_edges, source='node1', target='node2', create_using=nx.DiGraph())
        self.root_nodes = [tree.root.name]
        self.leaves = [n.name for n in tree.get_terminals()]

    def read_distance_matrix(self) -> pd.DataFrame:
        if os.path.exists(self.distance_matrix_filepath):

            distance_matrix = pd.read_csv(self.distance_matrix_filepath, sep="\t", index_col=0)

            distance_matrix = distance_matrix.loc[~distance_matrix.index.duplicated(keep='first')]
            distance_matrix = distance_matrix.loc[:, ~distance_matrix.columns.duplicated(keep='first')]
            nodes_to_filter = self.m_stats_matrix[self.m_stats_matrix['coverage'] > 0].index
            
            distance_matrix = distance_matrix.loc[distance_matrix.index.isin(nodes_to_filter), distance_matrix.columns.isin(nodes_to_filter)]
            
            if distance_matrix.shape[0] < 1:
                return pd.DataFrame()

            distance_matrix.index = distance_matrix.index.map(str)
            distance_matrix.columns = distance_matrix.columns.map(str)
            return distance_matrix

    def new_tree_from_distance_matrix(self):
        """
        original distance matrix :
        - asymetric
        - dist[i][j] = 1 - (proportion of i reads also present in j)
        """
        if os.path.exists(self.distance_matrix_filepath):
            distance_matrix = self.read_distance_matrix()
            
            if distance_matrix.empty is True:
                return
            if distance_matrix.shape[0] == 1:
                self.leaves = list(distance_matrix.index)
                self.all_nodes = list(distance_matrix.index)
                self.root_nodes = list(distance_matrix.index)
                self.tree = nx.DiGraph()

            proximity_matrix = 1 - distance_matrix
            weighted_proximity_matrix = self.weighted_matrix(proximity_matrix)
            weighted_distance_matrix = 1 - weighted_proximity_matrix

            weighted_distance_matrix.to_csv(os.path.join(self.output_dir, "weighted_distance_matrix.tsv"), sep="\t", index=True)
            # standardize weights by row
            symweighted_distance_matrix = self.matrix_symmetrize_mean(weighted_distance_matrix)
            self.distance_mat = weighted_distance_matrix
            self.njbio_from_distance_matrix(symweighted_distance_matrix)
            self.root_nodes = [n for n,d in self.tree.in_degree() if d==0]
            self.all_nodes = list(self.tree.nodes())
            self.leaves = [n for n,d in self.tree.out_degree() if d==0]

            self.recreate_all_node_stats()

            self.all_edges = pd.DataFrame(self.tree.edges(), columns=['node1', 'node2'])

    def recreate_all_node_stats(self):
        
        new_stats = []
        for node in self.all_nodes:
            leaves = self.get_node_leaves(node)
            original_stats = self.all_node_stats[self.all_node_stats['Node'] == node]
            private_reads = sum(original_stats['Private_Reads']) if not original_stats.empty else 0
            new_stats.append({
                'Node': node,
                'Num_Leaves': len(leaves),
                'Private_Reads': private_reads,
                'Private_Proportion': 1,
            })
        
        self.all_node_stats = pd.DataFrame(new_stats)
        if self.all_node_stats.empty is True:
            self.all_node_stats = pd.DataFrame(columns=['Node', 'Num_Leaves', 'Private_Reads', 'Private_Proportion'])


    def weighted_proximity(self, leaf1: str, leaf2: str, distance_matrix: pd.DataFrame) -> float:
        if leaf1 not in distance_matrix.index or leaf2 not in distance_matrix.columns:
            return float('nan')
                
        dist_ij = distance_matrix.loc[leaf1, leaf2]
        if leaf1 == leaf2:
            return 1.0

         ### Weight by coverage and mapping quality
        i_stats = self.m_stats_matrix.loc[leaf1] if leaf1 in self.m_stats_matrix.index else None
        j_stats = self.m_stats_matrix.loc[leaf2] if leaf2 in self.m_stats_matrix.index else None
        
        if i_stats is None or j_stats is None:
            return 0.0

        weight_ij_er = (j_stats['error_rate'] / i_stats['error_rate']) if i_stats['error_rate'] > 0 else 1
        weight_ij_er = 1 / weight_ij_er if weight_ij_er > 1 else weight_ij_er
        weight_ij =  weight_ij_er
        wdist_ij = dist_ij * weight_ij
        #wdist_ij = wdist_ij * 
        #####
        mapped_reads_i = i_stats['numreads'] if i_stats is not None else 0
        shared_reads_ij = mapped_reads_i * dist_ij
        nom = shared_reads_ij * j_stats['error_rate']
        denom = nom + mapped_reads_i * i_stats['error_rate']

        return wdist_ij if wdist_ij < 1 else 1 / wdist_ij
    
    def weighted_matrix(self, distance_matrix: pd.DataFrame) -> pd.DataFrame:
        weighted_mat = distance_matrix.copy()
        # remove duplicated indices and columns
        
        for i in distance_matrix.index:
            for j in distance_matrix.columns:
                weighted_mat.loc[i,j] = self.weighted_proximity(i, j, distance_matrix)

        return weighted_mat
    
    def matrix_symmetrize_mean(self, distance_matrix: pd.DataFrame) -> pd.DataFrame:
        sym_matrix = distance_matrix.copy()
        for i in distance_matrix.index:
            for j in distance_matrix.columns:
                if i != j:
                    sym_matrix.loc[i,j] = np.min([distance_matrix.loc[i,j], distance_matrix.loc[j,i]])
                    sym_matrix.loc[j,i] = sym_matrix.loc[i,j]

        return sym_matrix


    def recalculate_all_min_pairwise_dist(self):
        if os.path.exists(self.distance_matrix_filepath):
            #
            distance_matrix = self.read_distance_matrix()
            proximity_matrix = 1 - distance_matrix
            weighted_proximity_matrix = self.weighted_matrix(proximity_matrix)
            #

            def new_min_distance(row):   
                node = row['Node']
                leaves = self.get_leaves_parted(node)
                all_leaves = [x for sublist in leaves for x in sublist]

                other_leaves = [n for n in self.leaves if n not in all_leaves]
                if len(leaves) < 2:
                    row['Min_Pairwise_Dist'] = 0
                    row['Min_Shared'] = 0
                    row['Min_Dist'] = 0
                    return row

                distances = []
                leaves_left = leaves[0]
                leaves_right = leaves[1]
                distance_external = []



                for i in range(len(leaves_left)):
                    for j in range(len(leaves_right)):

                        wij_dist = weighted_proximity_matrix.loc[leaves_left[i], leaves_right[j]]
                        wji_dist = weighted_proximity_matrix.loc[leaves_right[j], leaves_left[i]]

                        #distances.append(max(ij_dist * weight_ij, j_dist * weight_ji))
                        distances.append(max(wij_dist, wji_dist))
                    
                    for ol in range(len(other_leaves)):
                        wij_dist = weighted_proximity_matrix.loc[leaves_left[i], other_leaves[ol]]
                        wji_dist = weighted_proximity_matrix.loc[other_leaves[ol], leaves_left[i]]

                        #distances.append(max(ij_dist * weight_ij, j_dist * weight_ji))
                        distance_external.append(max(wij_dist, wji_dist))


                min_distance = max(distances) if distances else 0.0

                row['Min_Dist'] = min_distance
                row['Min_Shared'] = max(distance_external) if distance_external else 0.0
                 ### Normalize min distance to be between 0 and 1
                row['Min_Pairwise_Dist'] = min_distance if min_distance < 1 else 1 / min_distance

                return row
            
            self.all_node_stats = self.all_node_stats[self.all_node_stats['Node'].isin(self.all_nodes)]
            self.all_node_stats= self.all_node_stats.apply(lambda row: new_min_distance(row), axis = 1)
            internal_nodes = self.all_node_stats[~self.all_node_stats['Node'].isin(self.leaves)]

            if internal_nodes.empty is True:
                self.all_node_stats['Z_min_dist'] = 0
                self.all_node_stats['Z_Min_Shared'] = 0
            else:
                self.all_node_stats['Z_min_dist'] = (self.all_node_stats['Min_Pairwise_Dist'] - internal_nodes['Min_Pairwise_Dist'].mean()) / (internal_nodes['Min_Pairwise_Dist'].std() if internal_nodes['Min_Pairwise_Dist'].std() > 0 else 1)
                self.all_node_stats['Z_Min_Shared'] = (self.all_node_stats['Min_Shared'] - internal_nodes['Min_Shared'].mean()) / (internal_nodes['Min_Shared'].std() if internal_nodes['Min_Shared'].std() > 0 else 1)

    def determine_min_pairwise_dist(self, node: str) -> float:
        if node not in self.all_nodes:
            return float('nan')
        if os.path.exists(self.distance_matrix_filepath):
            distance_matrix = pd.read_csv(self.distance_matrix_filepath, sep="\t", index_col=0)
            distance_matrix.index = distance_matrix.index.map(str)
            distance_matrix.columns = distance_matrix.columns.map(str)
            leaves = self.get_node_leaves(node)
            if len(leaves) < 2:
                return 0.0
            sub_matrix = distance_matrix.loc[leaves, leaves]
            tril_indices = np.tril_indices_from(sub_matrix, k=-1)
            pairwise_distances = sub_matrix.values[tril_indices]
            min_distance = pairwise_distances.min()
            return min_distance
        else:
            return float('nan')

    def remove_leaf_and_ascendants(self, node: str):
        """
        prevent dangling nodes.
        remove node in the tree and all its ascendants if they become leaves.
        """
        if node not in self.tree:
            return
        
        parents = list(self.tree.predecessors(node))
        self.tree.remove_node(node)
        for parent in parents:
            if self.tree.out_degree(parent) == 0:  # If parent has no other children
                self.remove_leaf_and_ascendants(parent)
        self.all_nodes = list(self.tree.nodes())
        self.root_nodes = [n for n,d in self.tree.in_degree() if d==0]
        self.leaves = [n for n,d in self.tree.out_degree() if d==0]

    def prune_empty_nodes(self):
        self.all_node_stats = self.all_node_stats.dropna()
        nodes_to_remove = [node for node in self.tree.nodes() if node not in self.all_node_stats['Node'].values]    


        self.m_stats_matrix['leaf'] = self.m_stats_matrix.apply(self.get_leaf, axis=1)
        self.m_stats_matrix = self.m_stats_matrix.dropna(subset=['leaf'])

        existing_leaves = self.m_stats_matrix[self.m_stats_matrix['coverage'] > 0.0]['leaf'].unique().tolist()
        nodes_to_remove = nodes_to_remove + [node for node in self.all_nodes if node not in existing_leaves]

        for node in nodes_to_remove:
            self.remove_leaf_and_ascendants(node)
        
        self.all_nodes = list(self.tree.nodes())
        self.all_edges = pd.DataFrame(self.tree.edges(), columns=['node1', 'node2'])
        if self.all_edges.empty:
            self.tree= nx.DiGraph() 
            self.tree.add_nodes_from(self.all_nodes)
        else:
            self.tree = nx.from_pandas_edgelist(self.all_edges, source='node1', target='node2', create_using=nx.DiGraph())

    def print_tree(self, threshold = 0.5):
        if self.tree.number_of_nodes() == 0:
            return None
        try:
            import matplotlib.pyplot as plt
            pos = nx.nx_agraph.graphviz_layout(self.tree, prog="dot")
            fig, ax = plt.subplots(figsize=(12, 8))
            nx.draw(self.tree, pos, with_labels=True, node_size=500, node_color="lightblue", font_size=10, font_weight="bold", ax=ax)
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
            color_list = plt.get_cmap('tab20').colors

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
                ax=ax
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
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
            import matplotlib.pyplot as plt

            color_list = plt.get_cmap('tab20').colors
            all_nodes_stats = self.all_node_stats.set_index('Node').to_dict('index')
            color_map = {node: f'rgb({int(color_list[i][0]*255)},{int(color_list[i][1]*255)},{int(color_list[i][2]*255)})' for i, node in enumerate(nodes)}
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
                x=edge_x, y=edge_y,
                line=dict(width=1, color='#888'),
                hoverinfo='none',
                mode='lines')

            node_x = []
            node_y = []
            for node in self.tree.nodes():
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)

            node_trace = go.Scatter(
                x=node_x, y=node_y,
                mode='markers+text',
                text=[f"{node}\n{all_nodes_stats[node]['Min_Pairwise_Dist']:.2f}" if node in all_nodes_stats else node for node in self.tree.nodes()],
                textposition="bottom center",
                hoverinfo='text',
                marker=dict(
                    showscale=False,
                    colorscale='YlGnBu',
                    color=node_colors,
                    size=20,
                    line_width=2))
            fig = go.Figure(data=[edge_trace, node_trace],
                            layout=go.Layout(
                                title="Tree Visualization",
                                showlegend=False,
                                hovermode='closest',
                                margin=dict(l=0, r=0, t=40, b=0),
                                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                            ))

            
            return fig
        except ImportError:
            print("Plotly is not installed. Please install it to visualize the tree.")
            return None
        except Exception as e:
            print(f"An error occurred while generating the tree visualization: {e}")
            return None

    def print_tree_with_clades(self, threshold=0.5):
        if self.tree.number_of_nodes() == 0:
            print("The tree is empty.")
            return None
        try:
            import matplotlib.pyplot as plt
            color_list = plt.get_cmap('tab20').colors

            import matplotlib.pyplot as plt
            node_colors = []
            all_nodes_stats = self.all_node_stats.set_index('Node')
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
            labels_with_min_dist = {node: f"{node}\n{all_nodes_stats.loc[node]['Min_Pairwise_Dist']:.2f}" for node in self.tree.nodes()}
            nx.draw(
                self.tree,
                pos,
                with_labels=True,
                node_size=500,
                node_color=node_colors,
                font_size=7,
                font_weight="bold",
                ax=ax,
                labels=labels_with_min_dist
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
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
            import matplotlib.pyplot as plt

            color_list = plt.get_cmap('tab20').colors

            node_colors = []
            all_nodes_stats = self.all_node_stats.set_index('Node').to_dict('index')
            selected_nodes_matrix = self.selected_nodes_matrix(threshold)
            
            color_map = {}

            color_index = 0

            for node in self.tree.nodes():
                cluster = self.cluster_map.get(node, None)
                if cluster:
                    if cluster not in color_map:
                        color_map[cluster] = f'rgb({int(color_list[color_index][0]*255)},{int(color_list[color_index][1]*255)},{int(color_list[color_index][2]*255)})'
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
                x=edge_x, y=edge_y,
                line=dict(width=1, color='#888'),
                hoverinfo='none',
                mode='lines')

            node_x = []
            node_y = []
            for node in self.tree.nodes():
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)

            node_trace = go.Scatter(
                x=node_x, y=node_y,
                mode='markers+text',
                text=[f"{node}\n{all_nodes_stats[node]['Min_Pairwise_Dist']:.2f}" if node in all_nodes_stats else node for node in self.tree.nodes()],
                textposition="bottom center",
                hoverinfo='text',
                marker=dict(
                    showscale=False,
                    colorscale='YlGnBu',
                    color=node_colors,
                    size=20,
                    line_width=2))
            fig = go.Figure(data=[edge_trace, node_trace],
                            layout=go.Layout(
                                showlegend=False,
                                hovermode='closest',
                                margin=dict(b=20,l=5,r=5,t=40),
                                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
            return fig
        except ImportError:
            print("plotly is not installed. Please install it to visualize the tree.")
            return None
        except Exception as e:
            print(f"An error occurred while generating the tree visualization: {e}")
            return None     


    def check_data_available(self):
        if os.path.exists(os.path.join(self.output_dir, "all_node_statistics.tsv")) and os.path.exists(os.path.join(self.output_dir, "nj_tree_edges.txt")):
            return True
        return False
    
    def node_selector(self, node: str, threshold = 0.5):
        if node in self.leaves:
            return True
        node_data = self.all_node_stats[self.all_node_stats['Node'] == node]
        if node_data.empty:
            return False
        return node_data["Min_Pairwise_Dist"].values[0] >= threshold and node_data["Min_Shared"].values[0] <= threshold
    
    def traverse_graph_recursive(self, node, threshold = 0.5, nodes_return = []) -> List[str]:
        if self.node_selector(node, threshold) == True:
            if node not in nodes_return:
                nodes_return.append(node)
                descendents = nx.descendants(self.tree, node)
                for desc in descendents:
                    self.cluster_map[desc] = node
                self.cluster_map[node] = node
        else:

            for neighbor in self.tree.neighbors(node):
                if neighbor not in nodes_return:
                    self.traverse_graph_recursive(neighbor, threshold, nodes_return)

        return nodes_return
    
    def get_node_leaves(self, node):
        
        descendants = nx.descendants(self.tree, node)
        leaves = [n for n in descendants if self.tree.out_degree(n) == 0]
        if self.tree.out_degree(node) == 0:
            leaves.append(node)
        return leaves
    
    def get_leaves_parted(self, node):
        if self.tree.out_degree(node) == 0:
            return [node]
        leaves = []
        for child in self.tree.successors(node):
            leaves.append(self.get_node_leaves(child))
        return leaves

    def selected_nodes_matrix(self, threshold = 0.5):
        self.cluster_map = {}

        selected_nodes = []
        for root in self.root_nodes:
            selected_nodes.extend(self.traverse_graph_recursive(root, threshold, []))
        selected_nodes = list(set(selected_nodes))

        stats_matrix = self.all_node_stats[self.all_node_stats['Node'].isin(selected_nodes)]
        def get_leaves_str(row):
            node = row['Node']
            leaves = self.get_node_leaves(node)
            row['leaves'] = ','.join(map(str, leaves))
            return row
        stats_matrix = stats_matrix.apply(get_leaves_str, axis=1)
        return stats_matrix
