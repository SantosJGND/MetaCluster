"""
Tests for diversity functions in metagenomics_utils.overlap_manager.diversity
"""
import pytest
import numpy as np
from metagenomics_utils.overlap_manager.diversity import (
    shannon_diversity,
    shannon_diversity_from_counts,
    shannon_diversity_from_list,
    skewness,
    kurtosis,
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
        taxa = ['A', 'A', 'B', 'B', 'C']
        result = shannon_diversity_from_list(taxa)
        assert result > 0

    def test_single_taxa(self):
        """Test with single taxa repeated."""
        taxa = ['A', 'A', 'A']
        result = shannon_diversity_from_list(taxa)
        assert result == 0.0

    def test_empty_list(self):
        """Test with empty list."""
        taxa = []
        result = shannon_diversity_from_list(taxa)
        assert result == 0.0

    def test_all_unique(self):
        """Test with all unique taxa."""
        taxa = ['A', 'B', 'C', 'D']
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
