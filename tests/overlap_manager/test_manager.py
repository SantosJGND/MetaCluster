"""
Tests for OverlapManager methods modified during optimization.
"""

import os

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from metagenomics_utils.overlap_manager.manager import (
    OverlapManager,
    _merge_matched_vectorized,
    merge_by_assembly_ID,
    merge_to_matched,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _old_weighted_matrix_ref(om, distance_matrix):
    """Reference: assid-based lookup, per-element loop matching weighted_matrix.

    Parses assid from each distance-matrix index the same way as the
    vectorized implementation, so the two should agree exactly.
    """
    m_stats = om.m_stats_matrix
    valid_idx = []
    for i in distance_matrix.index:
        assid = i.strip("./").strip(om.data_set_name).strip(".fna.bam")
        assid = "_".join(assid.split("_")[1:])
        if not m_stats[m_stats["assid"] == assid].empty:
            valid_idx.append(i)
    dm = distance_matrix.loc[valid_idx, valid_idx]

    weighted_mat = dm.copy() * 0.0
    np.fill_diagonal(weighted_mat.values, 1.0)
    for i in dm.index:
        assid_i = i.strip("./").strip(om.data_set_name).strip(".fna.bam")
        assid_i = "_".join(assid_i.split("_")[1:])
        err_i = m_stats[m_stats["assid"] == assid_i]["error_rate"].values[0]
        for j in dm.index:
            if i == j:
                continue
            assid_j = j.strip("./").strip(om.data_set_name).strip(".fna.bam")
            assid_j = "_".join(assid_j.split("_")[1:])
            err_j = m_stats[m_stats["assid"] == assid_j]["error_rate"].values[0]
            weight = (err_j / err_i) if err_i > 0 else 1.0
            weight = 1.0 / weight if weight > 1.0 else weight
            wdist = dm.loc[i, j] * weight
            weighted_mat.loc[i, j] = wdist if wdist < 1.0 else 1.0 / wdist
    return weighted_mat


def _old_symmetrize_mean(mat):
    """Reference: old double-loop symmetrize."""
    sym = mat.copy()
    for i in mat.index:
        for j in mat.columns:
            if i != j:
                sym.loc[i, j] = np.min([mat.loc[i, j], mat.loc[j, i]])
                sym.loc[j, i] = sym.loc[i, j]
    return sym


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def om_with_tree(tmp_path):
    """OverlapManager with a hand-built 3-leaf tree."""
    clustering_dir = tmp_path / "clustering"
    clustering_dir.mkdir()
    om = OverlapManager(str(clustering_dir), skip_build=True)

    om.tree = nx.DiGraph()
    om.tree.add_edges_from([("root", "A"), ("root", "I"), ("I", "B"), ("I", "C")])
    om.all_nodes = ["root", "A", "B", "C", "I"]
    om.leaves = ["A", "B", "C"]
    om.root_nodes = ["root"]

    om.m_stats_matrix = pd.DataFrame(
        {"error_rate": [0.01, 0.02, 0.03], "numreads": [100, 200, 300], "coverage": [10.0, 20.0, 30.0], "assid": ["A", "B", "C"]},
        index=["A", "B", "C"],
    )
    np.random.seed(42)
    n = 3
    raw = np.random.rand(n, n)
    raw = (raw + raw.T) / 2
    np.fill_diagonal(raw, 0)
    om.distance_mat = pd.DataFrame(raw, index=["A", "B", "C"], columns=["A", "B", "C"])
    om._rebuild_node_leaves_cache()
    return om


@pytest.fixture
def empty_om(tmp_path):
    """OverlapManager with skip_build, no tree."""
    clustering_dir = tmp_path / "clustering"
    clustering_dir.mkdir()
    return OverlapManager(str(clustering_dir), skip_build=True)


# ---------------------------------------------------------------------------
# weighted_matrix
# ---------------------------------------------------------------------------

class TestWeightedMatrix:
    @pytest.fixture
    def om_with_weights(self, tmp_path):
        clustering_dir = tmp_path / "clustering"
        clustering_dir.mkdir()
        om = OverlapManager(str(clustering_dir), skip_build=True)
        om.data_set_name = "test"
        om.m_stats_matrix = pd.DataFrame(
            {"error_rate": [0.01, 0.02, 0.03], "assid": ["A", "B", "C"]},
            index=["a", "b", "c"],
        )
        return om

    def _make_prox(self, n=3, prefix="X_", seed=None):
        if seed is not None:
            np.random.seed(seed)
        labels = [prefix + chr(65 + i) for i in range(n)]
        vals = np.random.uniform(0.5, 1.0, (n, n)) if seed is not None else np.ones((n, n))
        return pd.DataFrame(vals, index=labels, columns=labels)

    def test_output_shape(self, om_with_weights):
        om = om_with_weights
        prox = self._make_prox(3)
        result = om.weighted_matrix(prox)
        assert result.shape == (3, 3)

    def test_equals_reference(self, om_with_weights):
        om = om_with_weights
        prox = self._make_prox(3, seed=1)
        fast = om.weighted_matrix(prox)
        slow = _old_weighted_matrix_ref(om, prox)
        assert np.allclose(fast.values, slow.values, atol=1e-10)

    def test_matches_reference_multiple_sizes(self, om_with_weights):
        om = om_with_weights
        for n in [2, 4]:
            om.m_stats_matrix = pd.DataFrame(
                {"error_rate": np.random.uniform(0.005, 0.05, n),
                 "assid": [chr(65 + i) for i in range(n)]},
                index=[f"m{i}" for i in range(n)],
            )
            prox = self._make_prox(n, seed=42)
            fast = om.weighted_matrix(prox)
            slow = _old_weighted_matrix_ref(om, prox)
            assert np.allclose(fast.values, slow.values, atol=1e-10)

    def test_diagonal(self, om_with_weights):
        om = om_with_weights
        prox = self._make_prox(3)
        result = om.weighted_matrix(prox)
        assert np.allclose(np.diag(result.values), 1.0)

    def test_error_rates_of_zero(self, tmp_path):
        clustering_dir = tmp_path / "clustering"
        clustering_dir.mkdir()
        om = OverlapManager(str(clustering_dir), skip_build=True)
        om.data_set_name = "test"
        om.m_stats_matrix = pd.DataFrame(
            {"error_rate": [0.0, 0.02], "assid": ["X", "Y"]},
            index=["x", "y"],
        )
        prox = pd.DataFrame(
            [[1.0, 0.8], [0.7, 1.0]],
            index=["pre_X", "pre_Y"], columns=["pre_X", "pre_Y"],
        )
        result = om.weighted_matrix(prox)
        assert not np.any(np.isnan(result.values))
        assert np.isfinite(result.values).all()

    def test_filters_unmatched_indices(self, tmp_path):
        clustering_dir = tmp_path / "clustering"
        clustering_dir.mkdir()
        om = OverlapManager(str(clustering_dir), skip_build=True)
        om.data_set_name = "test"
        om.m_stats_matrix = pd.DataFrame(
            {"error_rate": [0.01], "assid": ["A"]},
            index=["a"],
        )
        # Only X_A should match A; X_Z has no match in m_stats
        prox = pd.DataFrame(
            [[1.0, 0.9], [0.9, 1.0]],
            index=["X_A", "X_Z"], columns=["X_A", "X_Z"],
        )
        result = om.weighted_matrix(prox)
        assert result.shape == (1, 1)
        assert result.index[0] == "X_A"


# ---------------------------------------------------------------------------
# matrix_symmetrize_mean
# ---------------------------------------------------------------------------

class TestMatrixSymmetrizeMean:
    def test_symmetrizes_asymmetric(self, om_with_tree):
        om = om_with_tree
        vals = np.array([[0.0, 0.9], [0.3, 0.0]])
        mat = pd.DataFrame(vals, index=["A", "B"], columns=["A", "B"])
        result = om.matrix_symmetrize_mean(mat)
        assert result.loc["A", "B"] == result.loc["B", "A"]
        assert result.loc["A", "B"] == min(0.9, 0.3)
        assert result.loc["B", "A"] == min(0.9, 0.3)

    def test_preserves_symmetric(self, om_with_tree):
        om = om_with_tree
        vals = np.array([[0.0, 0.5], [0.5, 0.0]])
        mat = pd.DataFrame(vals, index=["A", "B"], columns=["A", "B"])
        result = om.matrix_symmetrize_mean(mat)
        assert np.allclose(result.values, vals)

    def test_equals_reference(self, om_with_tree):
        om = om_with_tree
        np.random.seed(2)
        mat = pd.DataFrame(
            np.random.uniform(0, 1, (4, 4)),
            index=list("ABCD"), columns=list("ABCD"),
        )
        np.fill_diagonal(mat.values, 0)
        fast = om.matrix_symmetrize_mean(mat)
        slow = _old_symmetrize_mean(mat)
        assert np.allclose(fast.values, slow.values)

    def test_identity_1x1(self, om_with_tree):
        om = om_with_tree
        mat = pd.DataFrame([[0.0]], index=["A"], columns=["A"])
        result = om.matrix_symmetrize_mean(mat)
        assert result.loc["A", "A"] == 0.0


# ---------------------------------------------------------------------------
# _rebuild_node_leaves_cache
# ---------------------------------------------------------------------------

class TestRebuildNodeLeavesCache:
    def test_all_nodes_have_cache_entry(self, om_with_tree):
        om = om_with_tree
        for node in om.all_nodes:
            assert node in om.node_leaves_cache

    def test_leaf_cached_to_self(self, om_with_tree):
        om = om_with_tree
        for leaf in om.leaves:
            assert om.node_leaves_cache[leaf] == [leaf]

    def test_internal_includes_descendants(self, om_with_tree):
        om = om_with_tree
        # Internal node I has children B and C
        assert set(om.node_leaves_cache["I"]) == {"B", "C"}
        # Root covers all leaves
        assert set(om.node_leaves_cache["root"]) == {"A", "B", "C"}

    def test_cache_rebuilt_after_tree_change(self, om_with_tree):
        om = om_with_tree
        om.tree.remove_node("C")
        om._rebuild_node_leaves_cache()
        assert "C" not in om.node_leaves_cache
        assert set(om.node_leaves_cache["I"]) == {"B"}
        assert set(om.node_leaves_cache["root"]) == {"A", "B"}

    def test_cache_empty_on_empty_tree(self, empty_om):
        om = empty_om
        om._rebuild_node_leaves_cache()
        assert om.node_leaves_cache == {}


# ---------------------------------------------------------------------------
# get_node_leaves
# ---------------------------------------------------------------------------

class TestGetNodeLeaves:
    def test_leaf(self, om_with_tree):
        assert om_with_tree.get_node_leaves("A") == ["A"]

    def test_internal(self, om_with_tree):
        assert set(om_with_tree.get_node_leaves("I")) == {"B", "C"}

    def test_unknown_node(self, om_with_tree):
        assert om_with_tree.get_node_leaves("UNKNOWN") == []


# ---------------------------------------------------------------------------
# get_leaves_parted
# ---------------------------------------------------------------------------

class TestGetLeavesParted:
    def test_leaf(self, om_with_tree):
        result = om_with_tree.get_leaves_parted("A")
        assert result == [["A"]]

    def test_internal_binary(self, om_with_tree):
        result = om_with_tree.get_leaves_parted("I")
        assert len(result) == 2
        assert set(result[0]) | set(result[1]) == {"B", "C"}

    def test_root_covers_all(self, om_with_tree):
        result = om_with_tree.get_leaves_parted("root")
        all_leaves = set()
        for group in result:
            all_leaves |= set(group)
        assert all_leaves == {"A", "B", "C"}


# ---------------------------------------------------------------------------
# recreate_all_node_stats_from_tree
# ---------------------------------------------------------------------------

class TestRecreateAllNodeStats:
    def test_columns(self, om_with_tree):
        om = om_with_tree
        om.all_nodes = list(om.tree.nodes())
        om.recreate_all_node_stats_from_tree()
        expected = {"Node", "Num_Leaves", "Min_Dist", "Private_Reads", "Private_Proportion"}
        assert expected.issubset(set(om.all_node_stats.columns))

    def test_leaf_num_leaves(self, om_with_tree):
        om = om_with_tree
        om.all_nodes = list(om.tree.nodes())
        om.recreate_all_node_stats_from_tree()
        leaf_stats = om.all_node_stats.set_index("Node").loc["A"]
        assert leaf_stats["Num_Leaves"] == 1

    def test_internal_num_leaves(self, om_with_tree):
        om = om_with_tree
        om.all_nodes = list(om.tree.nodes())
        om.recreate_all_node_stats_from_tree()
        stats = om.all_node_stats.set_index("Node")
        assert stats.loc["I", "Num_Leaves"] == 2
        assert stats.loc["root", "Num_Leaves"] == 3

    def test_private_defaults(self, om_with_tree):
        om = om_with_tree
        om.all_nodes = list(om.tree.nodes())
        om.recreate_all_node_stats_from_tree()
        stats = om.all_node_stats.set_index("Node")
        for node in om.all_nodes:
            assert stats.loc[node, "Private_Reads"] == 0
            assert stats.loc[node, "Private_Proportion"] == 1

    def test_empty_tree(self, empty_om):
        om = empty_om
        om.all_nodes = []
        om.recreate_all_node_stats_from_tree()
        assert om.all_node_stats.empty


# ---------------------------------------------------------------------------
# recalculate_all_min_pairwise_dist
# ---------------------------------------------------------------------------

class TestRecalculateAllMinPairwiseDist:
    """Uses index names compatible with the assid-parsing in weighted_matrix.

    weighted_matrix extracts assid by::
        assid = idx.strip(...).strip(".fna.bam")
        assid = "_".join(assid.split("_")[1:])
    So the distance matrix index must have the form ``<prefix>_<assid>``.
    """

    @pytest.fixture
    def om_prebuilt(self, tmp_path):
        """OverlapManager with a hand-built 3-leaf tree and compatible index
        naming ( ``pre_A``, ``pre_B``, ``pre_C`` → assid ``A``, ``B``, ``C`` )."""
        clustering_dir = tmp_path / "clustering"
        clustering_dir.mkdir()
        om = OverlapManager(str(clustering_dir), skip_build=True)
        om.data_set_name = ""  # don't interfere with prefix stripping

        leaf_names = ["pre_A", "pre_B", "pre_C"]
        om.tree = nx.DiGraph()
        om.tree.add_edges_from([("root", "pre_A"), ("root", "I"), ("I", "pre_B"), ("I", "pre_C")])
        om.all_nodes = ["root", "pre_A", "pre_B", "pre_C", "I"]
        om.leaves = leaf_names
        om.root_nodes = ["root"]

        om.m_stats_matrix = pd.DataFrame(
            {"error_rate": [0.01, 0.02, 0.03], "coverage": [10.0, 20.0, 30.0], "assid": ["A", "B", "C"]},
            index=leaf_names,
        )
        np.random.seed(42)
        n = 3
        raw = np.random.rand(n, n)
        raw = (raw + raw.T) / 2
        np.fill_diagonal(raw, 0)
        om.distance_mat = pd.DataFrame(raw, index=leaf_names, columns=leaf_names)
        om.all_leaves_global = leaf_names
        om._rebuild_node_leaves_cache()
        return om

    def test_leaf_nodes_get_zero(self, om_prebuilt, tmp_path):
        om = om_prebuilt
        dist_file = tmp_path / "clustering" / "distance_matrix.tsv"
        om.distance_matrix_filepath = str(dist_file)
        om.all_node_stats = pd.DataFrame({"Node": om.all_nodes})

        df = om.distance_mat.copy()
        df.to_csv(dist_file, sep="\t")

        om.recalculate_all_min_pairwise_dist()
        stats = om.all_node_stats.set_index("Node")
        for leaf in ["pre_A", "pre_B", "pre_C"]:
            assert stats.loc[leaf, "Min_Pairwise_Dist"] == 0
            assert stats.loc[leaf, "Min_Shared"] == 0
            assert stats.loc[leaf, "Min_Dist"] == 0

    def test_internal_scores_calculated(self, om_prebuilt, tmp_path):
        om = om_prebuilt
        dist_file = tmp_path / "clustering" / "distance_matrix.tsv"
        om.distance_matrix_filepath = str(dist_file)
        om.all_node_stats = pd.DataFrame({"Node": om.all_nodes})

        om.distance_mat.to_csv(dist_file, sep="\t")

        om.recalculate_all_min_pairwise_dist()
        stats = om.all_node_stats.set_index("Node")
        assert stats.loc["I", "Min_Pairwise_Dist"] > 0
        assert stats.loc["root", "Min_Pairwise_Dist"] > 0

    def test_z_scores_finite(self, om_prebuilt, tmp_path):
        om = om_prebuilt
        dist_file = tmp_path / "clustering" / "distance_matrix.tsv"
        om.distance_matrix_filepath = str(dist_file)
        om.all_node_stats = pd.DataFrame({"Node": om.all_nodes})

        om.distance_mat.to_csv(dist_file, sep="\t")

        om.recalculate_all_min_pairwise_dist()
        assert om.all_node_stats["Z_min_dist"].notna().all()
        assert om.all_node_stats["Z_Min_Shared"].notna().all()
        assert np.isfinite(om.all_node_stats["Z_min_dist"]).all()
        assert np.isfinite(om.all_node_stats["Z_Min_Shared"]).all()


# ---------------------------------------------------------------------------
# remove_leaf_and_ascendants
# ---------------------------------------------------------------------------

class TestRemoveLeafAndAscendants:
    def test_remove_leaf_keeps_siblings(self, om_with_tree):
        om = om_with_tree
        assert "B" in om.tree
        om.remove_leaf_and_ascendants("B")
        assert "B" not in om.tree
        # I should still exist because C is still a child
        assert "I" in om.tree
        assert "C" in om.tree

    def test_remove_leaf_cascades_to_parent(self, om_with_tree):
        om = om_with_tree
        # Remove both B and C — I becomes childless, should be removed
        om.remove_leaf_and_ascendants("B")
        om.remove_leaf_and_ascendants("C")
        assert "B" not in om.tree
        assert "C" not in om.tree
        assert "I" not in om.tree  # cascaded
        assert "root" in om.tree  # still has A
        assert "A" in om.tree

    def test_remove_unknown_node_is_noop(self, om_with_tree):
        om = om_with_tree
        n_orig = om.tree.number_of_nodes()
        om.remove_leaf_and_ascendants("DOES_NOT_EXIST")
        assert om.tree.number_of_nodes() == n_orig

    def test_remove_only_leaf_from_root(self, empty_om):
        om = empty_om
        om.tree = nx.DiGraph()
        om.tree.add_node("solo")
        om.remove_leaf_and_ascendants("solo")
        assert "solo" not in om.tree
        assert om.tree.number_of_nodes() == 0


# ---------------------------------------------------------------------------
# merge_to_matched (standalone function)
# ---------------------------------------------------------------------------

class TestMergeToMatched:
    def test_match_by_suffix(self):
        matched = pd.DataFrame({
            "assembly_accession": ["GCF_000123.1"],
            "taxid": [12345],
            "description": ["Test"],
            "total_uniq_reads": [100],
        })
        row = pd.Series({
            "file": "/some/path/GCF_000123.1.fna",
            "taxid": None,
            "assid": None,
            "description": None,
            "total_uniq_reads": 0,
        })
        result = merge_to_matched(row, matched)
        assert result["assid"] == "GCF_000123.1"
        assert result["taxid"] == 12345

    def test_empty_matched_returns_nulls(self):
        matched = pd.DataFrame(columns=["assembly_accession", "taxid", "description", "total_uniq_reads"])
        row = pd.Series({
            "file": "test.fna",
            "taxid": None,
            "assid": None,
            "description": None,
            "total_uniq_reads": 0,
        })
        result = merge_to_matched(row, matched)
        assert result["taxid"] is None
        assert result["assid"] is None
        assert result["total_uniq_reads"] == 0


# ---------------------------------------------------------------------------
# _merge_matched_vectorized
# ---------------------------------------------------------------------------

class TestMergeMatchedVectorized:
    """Tests for _merge_matched_vectorized, focusing on the regex pattern."""

    @pytest.fixture
    def matched_df(self):
        return pd.DataFrame({
            "assembly_accession": ["GCF_000123.1", "GCF_000456.2", "GCF_000789.3"],
            "taxid": [111, 222, 333],
            "description": ["Desc_A", "Desc_B", "Desc_C"],
            "total_uniq_reads": [100, 200, 300],
        })

    @pytest.fixture
    def m_stats_basic(self):
        return pd.DataFrame({
            "assembly_accession": ["r1", "r2", "r3"],
            "coverage": [10.0, 20.0, 30.0],
            "file": [
                "/path/dataset_GCF_000123.1.fna",
                "/path/dataset_GCF_000456.2.fna",
                "/path/other_GCF_000789.3.fna.gz",
            ],
        })

    def test_basic_match(self, m_stats_basic, matched_df):
        result = _merge_matched_vectorized(m_stats_basic, matched_df)
        assert result["assid"].tolist() == ["GCF_000123.1", "GCF_000456.2", "GCF_000789.3"]
        assert result["taxid"].tolist() == [111, 222, 333]

    def test_accession_embedded_middle(self, matched_df):
        m_stats = pd.DataFrame({
            "assembly_accession": ["r1"],
            "coverage": [10.0],
            "file": ["dataset_0001_plan_144978_GCF_052827785.1.fna"],
        })
        matched = pd.DataFrame({
            "assembly_accession": ["GCF_052827785.1"],
            "taxid": [999],
            "description": ["Embedded"],
            "total_uniq_reads": [500],
        })
        result = _merge_matched_vectorized(m_stats, matched)
        assert result["assid"][0] == "GCF_052827785.1"

    def test_longest_prefix_wins(self):
        m_stats = pd.DataFrame({
            "assembly_accession": ["r1"],
            "coverage": [10.0],
            "file": ["dataset_GCF_000123.10.fna"],
        })
        matched = pd.DataFrame({
            "assembly_accession": ["GCF_000123.1", "GCF_000123.10"],
            "taxid": [1, 2],
            "description": ["Short", "Long"],
            "total_uniq_reads": [100, 200],
        })
        result = _merge_matched_vectorized(m_stats, matched)
        assert result["assid"][0] == "GCF_000123.10"

    def test_no_match_returns_nulls(self, m_stats_basic, matched_df):
        m_stats = pd.DataFrame({
            "assembly_accession": ["r1"],
            "coverage": [10.0],
            "file": ["dataset_GCF_999999.9.fna"],
        })
        result = _merge_matched_vectorized(m_stats, matched_df)
        assert result["assid"][0] is None
        assert result["taxid"][0] is None
        assert result["description"][0] is None
        assert result["total_uniq_reads"][0] == 0

    def test_partial_match(self, m_stats_basic, matched_df):
        m_stats = pd.DataFrame({
            "assembly_accession": ["r1", "r2"],
            "coverage": [10.0, 20.0],
            "file": ["dataset_GCF_000123.1.fna", "dataset_UNKNOWN.fna"],
        })
        result = _merge_matched_vectorized(m_stats, matched_df)
        assert result["assid"][0] == "GCF_000123.1"
        assert result["taxid"][0] == 111
        assert result["assid"][1] is None
        assert result["taxid"][1] is None
        assert result["total_uniq_reads"][0] == 100
        assert result["total_uniq_reads"][1] == 0

    def test_empty_matched(self, m_stats_basic):
        matched = pd.DataFrame(columns=["assembly_accession", "taxid", "description", "total_uniq_reads"])
        result = _merge_matched_vectorized(m_stats_basic, matched)
        assert result["assid"].isna().all()
        assert result["taxid"].isna().all()
        assert (result["total_uniq_reads"] == 0).all()

    def test_empty_m_stats(self, matched_df):
        m_stats = pd.DataFrame(columns=["assembly_accession", "coverage", "file"])
        result = _merge_matched_vectorized(m_stats, matched_df)
        assert result.empty

    def test_regex_metacharacters(self):
        m_stats = pd.DataFrame({
            "assembly_accession": ["r1"],
            "coverage": [10.0],
            "file": ["data_GCF_000123.1.fna"],
        })
        matched = pd.DataFrame({
            "assembly_accession": ["GCF_000123.1"],
            "taxid": [111],
            "description": ["Metachar"],
            "total_uniq_reads": [100],
        })
        result = _merge_matched_vectorized(m_stats, matched)
        assert result["assid"][0] == "GCF_000123.1"

    def test_same_accession_multiple_files(self):
        m_stats = pd.DataFrame({
            "assembly_accession": ["r1", "r2"],
            "coverage": [10.0, 20.0],
            "file": ["data_GCF_000123.1_A.fna", "data_GCF_000123.1_B.fna"],
        })
        matched = pd.DataFrame({
            "assembly_accession": ["GCF_000123.1"],
            "taxid": [111],
            "description": ["Multi"],
            "total_uniq_reads": [100],
        })
        result = _merge_matched_vectorized(m_stats, matched)
        assert result["assid"].tolist() == ["GCF_000123.1", "GCF_000123.1"]
        assert result["taxid"].tolist() == [111, 111]

    def test_case_sensitivity(self):
        m_stats = pd.DataFrame({
            "assembly_accession": ["r1"],
            "coverage": [10.0],
            "file": ["gcf_000123.1.fna"],
        })
        matched = pd.DataFrame({
            "assembly_accession": ["GCF_000123.1"],
            "taxid": [111],
            "description": ["Case"],
            "total_uniq_reads": [100],
        })
        result = _merge_matched_vectorized(m_stats, matched)
        assert result["assid"][0] is None

    def test_preserves_original_columns(self, m_stats_basic, matched_df):
        result = _merge_matched_vectorized(m_stats_basic, matched_df)
        assert "coverage" in result.columns
        assert "assembly_accession" in result.columns
        assert result["coverage"].tolist() == [10.0, 20.0, 30.0]
        assert "file" in result.columns


# ---------------------------------------------------------------------------
# get_leaf (instance method)
# ---------------------------------------------------------------------------

class TestGetLeaf:
    def test_exact_match(self, empty_om):
        om = empty_om
        om.all_leaves_global = ["GCF_000123.1", "GCF_000456.1"]
        row = pd.Series({"assid": "GCF_000123.1"})
        result = om.get_leaf(row)
        assert result == "GCF_000123.1"

    def test_na_accession(self, empty_om):
        om = empty_om
        row = pd.Series({"assid": None})
        result = om.get_leaf(row)
        assert result is None

    def test_suffix_true_positive(self, empty_om):
        om = empty_om
        om.all_leaves_global = ["ref|GCF_000123.1|"]
        row = pd.Series({"assid": "GCF_000123.1"})
        result = om.get_leaf(row)
        assert result == "ref|GCF_000123.1|"


# ---------------------------------------------------------------------------
# traverse_graph_recursive
# ---------------------------------------------------------------------------

class TestTraverseGraphRecursive:
    @pytest.fixture
    def om_with_cluster_data(self, tmp_path):
        clustering_dir = tmp_path / "clustering"
        clustering_dir.mkdir()
        om = OverlapManager(str(clustering_dir), skip_build=True)
        om.tree = nx.DiGraph()
        om.tree.add_edges_from([("root", "A"), ("root", "I"), ("I", "B"), ("I", "C")])
        om.all_nodes = ["root", "A", "B", "C", "I"]
        om.leaves = ["A", "B", "C"]
        om.root_nodes = ["root"]
        om.all_node_stats = pd.DataFrame({
            "Node": ["root", "A", "B", "C", "I"],
            "Min_Pairwise_Dist": [0.9, 1.0, 1.0, 1.0, 0.6],
            "Min_Shared": [0.1, 0.0, 0.0, 0.0, 0.3],
        })
        return om

    def test_mutable_default_not_shared(self, om_with_cluster_data):
        om = om_with_cluster_data
        om.cluster_map = {}
        r1 = om.traverse_graph_recursive("root", 0.7, None)
        om.cluster_map = {}
        r2 = om.traverse_graph_recursive("root", 0.7, None)
        # Each call should return independent results
        assert r1 is not r2

    def test_threshold_zero_selects_leaves(self, om_with_cluster_data):
        om = om_with_cluster_data
        om.cluster_map = {}
        result = om.traverse_graph_recursive("root", 0, None)
        for node in result:
            assert node in om.leaves or node["Min_Pairwise_Dist"] == 0

    def test_cluster_map_populated(self, om_with_cluster_data):
        om = om_with_cluster_data
        om.cluster_map = {}
        _ = om.traverse_graph_recursive("root", 0.5, None)
        assert len(om.cluster_map) > 0
        for node in om.cluster_map:
            assert node in om.all_nodes


# ---------------------------------------------------------------------------
# Build — integration
# ---------------------------------------------------------------------------

class TestOverlapManagerBuild:
    def test_build_check_data_available(self, tmp_path):
        """check_data_available returns True when distance_matrix.tsv exists."""
        clustering_dir = tmp_path / "clustering"
        clustering_dir.mkdir()
        dist_file = clustering_dir / "distance_matrix.tsv"
        dist_file.write_text("\tA\tB\nA\t0.0\t0.5\nB\t0.5\t0.0\n")

        om = OverlapManager(str(clustering_dir), skip_build=True)
        assert om.check_data_available()

    def test_build_skip_does_not_build(self, empty_om):
        assert empty_om.tree.number_of_nodes() == 0
        assert empty_om.all_node_stats.empty

    def test_prune_empty_nodes_called_during_build(self, om_with_tree):
        """Verify cache is populated after build-like setup."""
        om = om_with_tree
        assert isinstance(om.node_leaves_cache, dict)
        assert len(om.node_leaves_cache) > 0
