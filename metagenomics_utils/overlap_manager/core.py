"""
Core functions for read overlap analysis.
"""
from typing import List

import numpy as np
import pandas as pd


def pairwise_shared_count(read_profile_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Return dataframe of pairwise shared read counts,
    use matrix multiplication to sum shared reads from binary matrix for each pair.
    """
    binary_matrix = np.array(read_profile_matrix)

    shared_reads = []
    for i in range(binary_matrix.shape[0]):
        row = binary_matrix[i].reshape(1, -1)
        prod0 = row @ binary_matrix.T
        shared_reads.append(prod0)

    if len(shared_reads) == 0:
        return pd.DataFrame(
            index=read_profile_matrix.index, columns=read_profile_matrix.index
        )
    shared_reads = np.concatenate(shared_reads, axis=0)

    shared_reads = pd.DataFrame(
        shared_reads,
        index=read_profile_matrix.index,
        columns=read_profile_matrix.index,
    )

    return shared_reads


def square_and_fill_diagonal(clade_read_matrix: pd.DataFrame) -> pd.DataFrame:
    """Calculate pairwise shared proportions with diagonal set to 0."""
    shared_clade_matrix = pairwise_shared_count(clade_read_matrix)

    shared_clade_matrix = shared_clade_matrix.div(
        clade_read_matrix.sum(axis=1),
        axis=0,
    )

    np.fill_diagonal(shared_clade_matrix.values, 0)
    shared_clade_matrix = shared_clade_matrix.fillna(0)

    return shared_clade_matrix


def very_similar_groups_from_dataframe(
    read_profile_matrix_filtered: pd.DataFrame, threshold: float = 0.95
) -> List[tuple]:
    """
    Find very similar entries in pairwise shared clade.
    """
    shared_read_matrix = square_and_fill_diagonal(read_profile_matrix_filtered)

    clusters_assigment_dict = {}
    clusternum = 0

    for i in range(shared_read_matrix.shape[0]):
        for j in range(shared_read_matrix.shape[1]):
            if i == j:
                continue
            if (
                shared_read_matrix.iloc[i, j] >= threshold
                or shared_read_matrix.iloc[j, i] >= threshold
            ):
                assignments = [
                    clusters_assigment_dict.get(i, None),
                    clusters_assigment_dict.get(j, None),
                ]
                assignments = [x for x in assignments if x is not None]

                if not assignments:
                    clusters_assigment_dict[i] = clusternum
                    clusters_assigment_dict[j] = clusternum
                    clusternum += 1
                elif len(assignments) == 1:
                    clusters_assigment_dict[j] = assignments[0]
                else:
                    for a in assignments[1:]:
                        for k, v in clusters_assigment_dict.items():
                            if v == a:
                                clusters_assigment_dict[k] = assignments[0]

    groups = {}
    for idx, cluster in clusters_assigment_dict.items():
        if cluster not in groups:
            groups[cluster] = []
        groups[cluster].append(idx)

    return [tuple(v) for v in groups.values()]
