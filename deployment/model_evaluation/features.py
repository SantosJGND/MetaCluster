"""
Feature extraction transformer for recall model.

Transforms raw m_stats statistics matrices into feature vectors
suitable for training RecallModeller and its subclasses.
"""

import logging
from typing import List, Optional, Literal

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin

logger = logging.getLogger(__name__)

SORT_STRATEGIES = Literal['reads', 'taxid_roundrobin', 'rarity_boost', 'tax_level_stratified']


class RecallFeatureTransformer(TransformerMixin, BaseEstimator):
    """
    Collapses a raw m_stats DataFrame into a single-row feature vector.

    The transformer learns the reference taxon schema from taxids_to_use
    during ``fit()``, then during ``transform()`` computes:

    - **Statistical features**: ``number_hits``, ``counts_kurtosis``,
      ``counts_skewness``, ``tax_diversity_shannon``, ``max_uniq_reads``,
      ``total_uniq_reads``
    - **Taxonomic composition features**: proportion of total_uniq_reads
      per taxon at the configured ``tax_level``. Missing taxa get 0.0.
    - **Target columns**: ``last_best_match_relindex`` and
      ``index_recall_1`` … ``index_recall_{N}`` (recall at each division).

    Parameters
    ----------
    tax_level : str
        Taxonomic level for composition features (default ``'order'``).
    data_set_divide : int
        Number of recall divisions for target columns (default 5).
    sort_strategy : str
        How to sort leaves before computing recall divisions.
        Options:
        - ``'reads'``: sort by ``total_uniq_reads`` descending (current default).
        - ``'taxid_roundrobin'``: interleave by ``taxid`` group.
        - ``'rarity_boost'``: boost under-represented taxa.
        - ``'tax_level_stratified'``: interleave by taxonomic level.
    """

    def __init__(self, tax_level: str = 'order', data_set_divide: int = 5,
                 sort_strategy: SORT_STRATEGIES = 'reads'):
        self.tax_level = tax_level
        self.data_set_divide = data_set_divide
        self.sort_strategy = sort_strategy
        self.reference_taxa_: Optional[List[str]] = None
        self._feature_names: Optional[List[str]] = None

    def fit(
        self,
        X: Optional[List[pd.DataFrame]] = None,
        y: Optional[pd.DataFrame] = None,
        taxids_to_use: Optional[pd.DataFrame] = None,
    ) -> 'RecallFeatureTransformer':
        """
        Learn the reference taxon schema from ``taxids_to_use``.

        Parameters
        ----------
        X : list of DataFrames or None
            Ignored (accepted for sklearn compatibility).
        y : DataFrame or None
            Ignored.
        taxids_to_use : DataFrame or None
            DataFrame containing a column named ``self.tax_level`` with
            the reference set of taxa. Taxa not present in this set are
            mapped to ``'unclassified'``.

        Returns
        -------
        self
        """
        taxa = ['unclassified']
        if taxids_to_use is not None and self.tax_level in taxids_to_use.columns:
            ref = taxids_to_use[self.tax_level].dropna().unique().tolist()
            taxa = sorted(set(taxa + ref))
        self.reference_taxa_ = taxa

        stat_cols = [
            'number_hits', 'counts_kurtosis', 'counts_skewness',
            'tax_diversity_shannon', 'max_uniq_reads', 'total_uniq_reads',
        ]
        comp_cols = [str(t) for t in self.reference_taxa_]
        self._feature_names = stat_cols + comp_cols

        logger.debug(
            "RecallFeatureTransformer fitted: %d stat cols + %d composition cols "
            "(%d total features)",
            len(stat_cols), len(comp_cols), len(self._feature_names),
        )
        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        """
        Transform a single m_stats matrix into a single-row feature + target vector.

        Parameters
        ----------
        X : pd.DataFrame
            A single ``m_stats_stats_matrix`` DataFrame. Must contain at
            least the columns ``total_uniq_reads``, ``best_match_is_best``,
            and ``self.tax_level``.

        Returns
        -------
        pd.DataFrame
            Single-row DataFrame with stat features, composition features,
            and target columns.
        """
        m_stats = X.copy()
        m_stats = self._apply_sort(m_stats).reset_index(drop=True)

        from metagenomics_utils.overlap_manager.node_stats import (
            node_composition_with_stats,
        )

        tax_ref = pd.DataFrame({self.tax_level: self.reference_taxa_})
        composition = node_composition_with_stats(
            m_stats, tax_ref, tax_level=self.tax_level,
        )

        best_match_indices = m_stats.index[
            m_stats['best_match_is_best'] == True
        ].tolist()
        last_best_match_index = (
            best_match_indices[-1] + 1 if best_match_indices else 0
        )
        last_best_match_relindex = (
            last_best_match_index / len(m_stats) if len(m_stats) > 0 else 0.0
        )

        composition.insert(0, 'last_best_match_relindex', last_best_match_relindex)

        for d in reversed(range(1, self.data_set_divide + 1)):
            threshold = int(m_stats.shape[0] * d / self.data_set_divide)
            short = [i for i in best_match_indices if i < threshold]
            recall = len(short) / len(best_match_indices) if best_match_indices else 0.0
            composition.insert(0, f'index_recall_{d}', recall)

        composition.reset_index(drop=True, inplace=True)

        logger.debug(
            "RecallFeatureTransformer: %d features + %d targets -> %d columns",
            len(self._feature_names), len(self.target_columns_),
            composition.shape[1],
        )
        return composition

    # ── Sort strategies ──────────────────────────────────────────────

    def _apply_sort(self, m_stats: pd.DataFrame) -> pd.DataFrame:
        """Dispatch to the appropriate sort method based on sort_strategy."""
        if self.sort_strategy == 'reads':
            return self._sort_reads(m_stats)
        elif self.sort_strategy == 'taxid_roundrobin':
            return self._sort_taxid_roundrobin(m_stats)
        elif self.sort_strategy == 'rarity_boost':
            return self._sort_rarity_boost(m_stats)
        elif self.sort_strategy == 'tax_level_stratified':
            return self._sort_tax_level_stratified(m_stats)
        else:
            raise ValueError(
                f"Unknown sort strategy: {self.sort_strategy}. "
                f"Valid options: reads, taxid_roundrobin, rarity_boost, tax_level_stratified"
            )

    def _sort_reads(self, m_stats: pd.DataFrame) -> pd.DataFrame:
        """Baseline: sort by total_uniq_reads descending."""
        return m_stats.sort_values('total_uniq_reads', ascending=False)

    def _sort_taxid_roundrobin(self, m_stats: pd.DataFrame) -> pd.DataFrame:
        """Interleave leaves by taxid group.

        Within each taxid, leaves are sorted by total_uniq_reads descending.
        Then groups are interleaved: first-best leaf from each taxid,
        second-best from each taxid, etc. Maximises distinct taxids early.
        """
        df = m_stats.copy()
        df['_recall_sort_rank'] = df.groupby('taxid', sort=False).cumcount()
        result = df.sort_values(
            ['_recall_sort_rank', 'total_uniq_reads'],
            ascending=[True, False],
        )
        return result.drop(columns=['_recall_sort_rank'])

    def _sort_rarity_boost(self, m_stats: pd.DataFrame) -> pd.DataFrame:
        """Boost leaves from taxonomically under-represented taxa.

        Score = total_uniq_reads * (1 - taxon_read_fraction).
        Taxa with few total reads get a boost, pushing their leaves up.
        """
        df = m_stats.copy()
        taxon_totals = df.groupby(self.tax_level)['total_uniq_reads'].transform('sum')
        total_reads = df['total_uniq_reads'].sum()
        fraction = taxon_totals / total_reads if total_reads > 0 else 0
        df['_recall_sort_score'] = df['total_uniq_reads'] * (1.0 - fraction)
        result = df.sort_values('_recall_sort_score', ascending=False)
        return result.drop(columns=['_recall_sort_score'])

    def _sort_tax_level_stratified(self, m_stats: pd.DataFrame) -> pd.DataFrame:
        """Interleave leaves by taxonomic level group.

        Within each group (e.g. order), leaves are sorted by
        total_uniq_reads descending. Groups are then interleaved,
        ensuring each taxon is represented before taking second leaves.
        """
        df = m_stats.copy()
        df['_recall_sort_rank'] = df.groupby(self.tax_level, sort=False).cumcount()
        result = df.sort_values(
            ['_recall_sort_rank', 'total_uniq_reads'],
            ascending=[True, False],
        )
        return result.drop(columns=['_recall_sort_rank'])

    @property
    def target_columns_(self) -> List[str]:
        """Names of target columns produced by ``transform()``."""
        return ['last_best_match_relindex'] + [
            f'index_recall_{i}' for i in range(1, self.data_set_divide + 1)
        ]

    def get_feature_names_out(self) -> List[str]:
        """
        Names of feature columns (stat + composition) produced by ``transform()``.

        Only available after ``fit()`` has been called.
        """
        if self._feature_names is None:
            raise RuntimeError(
                "get_feature_names_out() called before fit(). "
                "Call fit() first."
            )
        return list(self._feature_names)
