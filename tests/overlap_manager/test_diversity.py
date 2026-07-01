"""
Tests for diversity functions in metagenomics_utils.overlap_manager.diversity
"""

import pytest

from metagenomics_utils.overlap_manager.diversity import (
    kurtosis,
    shannon_diversity,
    shannon_diversity_from_counts,
    shannon_diversity_from_list,
    skewness,
)


class TestShannonDiversity:
    """Tests for shannon_diversity function."""

    def test_uniform_distribution(self):
        """Test with uniform proportions."""
        proportions = [0.25, 0.25, 0.25, 0.25]
        result = shannon_diversity(proportions)
        assert result > 0
        assert result == pytest.approx(1.386, rel=0.01)

    def test_single_element(self):
        """Test with single element (maximum diversity)."""
        proportions = [1.0]
        result = shannon_diversity(proportions)
        assert result == 0.0

    def test_zero_proportion(self):
        """Test with zero proportion."""
        proportions = [0.0, 0.5, 0.5]
        result = shannon_diversity(proportions)
        assert result > 0

    def test_empty_list(self):
        """Test with empty list."""
        proportions = []
        result = shannon_diversity(proportions)
        assert result == 0.0


class TestShannonDiversityFromCounts:
    """Tests for shannon_diversity_from_counts function."""

    def test_equal_counts(self):
        """Test with equal counts."""
        counts = [10, 10, 10, 10]
        result = shannon_diversity_from_counts(counts)
        assert result > 0
        assert result == pytest.approx(1.386, rel=0.01)

    def test_single_count(self):
        """Test with single count."""
        counts = [100]
        result = shannon_diversity_from_counts(counts)
        assert result == 0.0

    def test_zero_total(self):
        """Test with all zeros."""
        counts = [0, 0, 0]
        result = shannon_diversity_from_counts(counts)
        assert result == 0.0


class TestShannonDiversityFromList:
    """Tests for shannon_diversity_from_list function."""

    def test_taxa_list(self):
        """Test with list of taxa."""
        taxa = ["A", "A", "B", "B", "C"]
        result = shannon_diversity_from_list(taxa)
        assert result > 0

    def test_single_taxa(self):
        """Test with single taxa repeated."""
        taxa = ["A", "A", "A"]
        result = shannon_diversity_from_list(taxa)
        assert result == 0.0

    def test_empty_list(self):
        """Test with empty list."""
        taxa = []
        result = shannon_diversity_from_list(taxa)
        assert result == 0.0

    def test_all_unique(self):
        """Test with all unique taxa."""
        taxa = ["A", "B", "C", "D"]
        result = shannon_diversity_from_list(taxa)
        assert result > 0


class TestSkewness:
    """Tests for skewness function."""

    def test_symmetric_distribution(self):
        """Test with symmetric distribution."""
        proportions = [0.25, 0.25, 0.25, 0.25]
        result = skewness(proportions)
        assert result == pytest.approx(0.0, abs=0.01)

    def test_right_skewed(self):
        """Test with right-skewed distribution."""
        proportions = [0.1, 0.2, 0.3, 0.4]
        result = skewness(proportions)
        assert result > 0

    def test_empty_list(self):
        """Test with empty list."""
        proportions = []
        result = skewness(proportions)
        assert result == 0.0


class TestKurtosis:
    """Tests for kurtosis function."""

    def test_normal_like(self):
        """Test with normal-like distribution."""
        proportions = [0.25, 0.25, 0.25, 0.25]
        result = kurtosis(proportions)
        assert result <= 0

    def test_empty_list(self):
        """Test with empty list."""
        proportions = []
        result = kurtosis(proportions)
        assert result == 0.0


class TestOverlapManagerMissingFiles:
    """Tests for OverlapManager handling missing all_node_statistics.tsv."""

    def test_missing_all_node_stats_creates_empty_df(self, tmp_path):
        """Test that missing all_node_statistics.tsv creates empty DataFrame."""
        clustering_dir = tmp_path / "clustering"
        clustering_dir.mkdir()

        distance_matrix = clustering_dir / "distance_matrix.tsv"
        distance_matrix.write_text("""	asm1	asm2
asm1	0.0	0.5
asm2	0.5	0.0
""")

        from metagenomics_utils.overlap_manager.manager import OverlapManager

        om = OverlapManager(str(clustering_dir), skip_build=True)

        assert om.all_node_stats.empty
        assert om.all_edges.empty

    def test_check_data_available_returns_true_when_distance_matrix_exists(self, tmp_path):
        """Test that check_data_available returns True when distance_matrix.tsv exists."""
        clustering_dir = tmp_path / "clustering"
        clustering_dir.mkdir()

        distance_matrix = clustering_dir / "distance_matrix.tsv"
        distance_matrix.write_text("""	asm1	asm2
asm1	0.0	0.5
asm2	0.5	0.0
""")

        from metagenomics_utils.overlap_manager.manager import OverlapManager

        om = OverlapManager(str(clustering_dir), skip_build=True)

        result = om.check_data_available()
        assert result is True

    def test_tree_rebuilt_from_distance_matrix(self, tmp_path):
        """Test that tree can be rebuilt from distance_matrix.tsv."""
        clustering_dir = tmp_path / "clustering"
        clustering_dir.mkdir()

        distance_matrix = clustering_dir / "distance_matrix.tsv"
        distance_matrix.write_text("""	asm1	asm2
asm1	0.0	0.5
asm2	0.5	0.0
""")

        from metagenomics_utils.overlap_manager.manager import OverlapManager

        om = OverlapManager(str(clustering_dir), skip_build=True)
        assert hasattr(om, "tree")
        assert hasattr(om, "recreate_all_node_stats_from_tree")
        assert callable(om.recreate_all_node_stats_from_tree)

    def test_recreate_all_node_stats_method_exists(self, tmp_path):
        """Test that OverlapManager has recreate_all_node_stats method."""
        clustering_dir = tmp_path / "clustering"
        clustering_dir.mkdir()

        distance_matrix = clustering_dir / "distance_matrix.tsv"
        distance_matrix.write_text("""	asm1	asm2
asm1	0.0	0.3
asm2	0.3	0.0
""")

        from metagenomics_utils.overlap_manager.manager import OverlapManager

        om = OverlapManager(str(clustering_dir), skip_build=True)

        assert hasattr(om, "recreate_all_node_stats_from_tree")
        assert callable(om.recreate_all_node_stats_from_tree)
