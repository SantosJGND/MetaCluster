"""
General diversity functions for ecological and taxonomic analysis.
"""
from typing import List
import numpy as np
from collections import Counter


def shannon_diversity(proportions: List[float]) -> float:
    """
    Calculate Shannon diversity index given a list of proportions.
    """
    return -sum(p * (p if p == 0 else np.log(p)) for p in proportions)


def shannon_diversity_from_counts(counts: List[int]) -> float:
    """
    Calculate Shannon diversity index given a list of counts.
    """
    total = sum(counts)
    if total == 0:
        return 0.0
    proportions = [count / total for count in counts]
    return shannon_diversity(proportions)


def shannon_diversity_from_list(taxa: List[str]) -> float:
    """
    Calculate Shannon diversity index given a list of taxa.
    """
    if not taxa:
        return 0.0
    counts = Counter(taxa)
    return shannon_diversity_from_counts(list(counts.values()))


def skewness(proportions: List[float]) -> float:
    """
    Calculate skewness of a distribution given a list of proportions.
    """
    mean = np.mean(proportions)
    std_dev = np.std(proportions)
    if std_dev == 0 or len(proportions) == 0:
        return 0.0
    skewness = sum((p - mean) ** 3 for p in proportions) / (len(proportions) * (std_dev ** 3))
    return skewness


def kurtosis(proportions: List[float]) -> float:
    """
    Calculate kurtosis of a distribution given a list of proportions.
    """
    mean = np.mean(proportions)
    std_dev = np.std(proportions)
    if std_dev == 0 or len(proportions) == 0:
        return 0.0
    kurtosis = sum((p - mean) ** 4 for p in proportions) / (len(proportions) * (std_dev ** 4)) - 3
    return kurtosis
