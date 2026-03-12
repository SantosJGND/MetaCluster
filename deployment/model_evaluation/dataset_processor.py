"""
Dataset processor class for evaluator module.

Processes a single dataset through the full evaluation pipeline.
"""

import os
import logging
from typing import Optional, Any

import pandas as pd
import numpy as np

from metagenomics_utils.overlap_manager import OverlapManager
from metagenomics_utils.overlap_manager.node_stats import get_m_stats_matrix
from metagenomics_utils.overlap_manager.om_models import (
    get_trash_composition,
    get_cross_hit_composition,
    calculate_overall_precision,
    predict_data_set_clades,
    cross_hit_prediction,
    cut_off_recall_prediction,
)

from .config import EvaluatorConfig, TrainedModels
from .result_models import (
    DatasetResult,
    PrecisionMetrics,
    RecallMetrics,
    CrossHitMetrics,
)
from .exceptions import DataLoadError, PredictionError, OverlapManagerError
from .metrics import (
    compute_fuzzy_precision,
    compute_overall_precision,
    compute_recall,
    compute_trash_flags,
    compute_clade_accuracy,
    compute_clade_recall,
    safe_divide,
)


logger = logging.getLogger(__name__)


class DatasetProcessor:
    """
    Processes a single dataset through the full evaluation pipeline.
    
    The pipeline consists of the following phases:
    1. Load data (OverlapManager, input TSV)
    2. Compute baseline metrics
    3. Apply recall filter
    4. Pre-cleanup clade prediction
    5. Cross-hit prediction & cleanup
    6. Post-cleanup prediction
    """
    
    def __init__(
        self,
        config: EvaluatorConfig,
        models: TrainedModels,
        ncbi_wrapper: Any,
        input_tax_df: pd.DataFrame,
        taxids_to_use: pd.DataFrame,
    ):
        """
        Initialize dataset processor.
        
        Args:
            config: EvaluatorConfig instance
            models: TrainedModels instance with trained models
            ncbi_wrapper: NCBI TaxonomistWrapper instance
            input_tax_df: Input taxid DataFrame with lineage info
            taxids_to_use: Taxids to use for evaluation
        """
        self.config = config
        self.models = models
        self.ncbi = ncbi_wrapper
        self.input_tax_df = input_tax_df
        self.taxids_to_use = taxids_to_use
    
    def process(self, data_set_name: str) -> Optional[DatasetResult]:
        """
        Process one dataset end-to-end.
        
        Args:
            data_set_name: Name of the dataset to process
            
        Returns:
            DatasetResult or None if dataset should be skipped
        """
        logger.info(f"Processing dataset: {data_set_name}")
        
        try:
            overlap_manager, input_summary = self._load_data(data_set_name)
            if overlap_manager is None:
                logger.warning(f"Skipping {data_set_name}: failed to load data")
                return None
            
            result = DatasetResult(
                data_set=data_set_name,
                sample=input_summary['sample'].iloc[0] if len(input_summary) > 0 else '',
                input_taxid_count=int(input_summary['taxid'].nunique()),
            )
            
            result = self._compute_baseline_metrics(overlap_manager, input_summary, result)
            result = self._apply_recall_filter(data_set_name, overlap_manager, result)
            result = self._predict_clades_precleanup(data_set_name, overlap_manager, result)
            result = self._apply_crosshit_cleanup(data_set_name, overlap_manager, result)
            result = self._predict_clades_postcleanup(data_set_name, overlap_manager, result)
            
            logger.info(f"Completed processing {data_set_name}")
            return result
            
        except FileNotFoundError as e:
            raise DataLoadError(data_set_name, str(e))
        except Exception as e:
            raise PredictionError(data_set_name, "processing", str(e)) from e
    
    def _load_data(self, data_set_name: str) -> tuple[Optional[OverlapManager], pd.DataFrame]:
        """
        Load OverlapManager and input data.
        
        Args:
            data_set_name: Name of dataset
            
        Returns:
            Tuple of (OverlapManager, input_summary) or (None, empty DataFrame)
        """
        om_path = os.path.join(
            self.config.study_output_filepath, 
            data_set_name, 
            "clustering"
        )
        input_path = os.path.join(
            self.config.study_output_filepath, 
            data_set_name, 
            "input", 
            f"{data_set_name}.tsv"
        )
        
        if not os.path.exists(input_path):
            logger.warning(f"Input file not found: {input_path}")
            raise FileNotFoundError(f"Input file not found: {input_path}")
        
        try:
            overlap_manager = OverlapManager(om_path)
        except Exception as e:
            raise OverlapManagerError(data_set_name, "load", str(e))
        
        if overlap_manager.m_stats_matrix.empty:
            logger.warning(f"No mapped reads for {data_set_name}")
            return None, pd.DataFrame()
        
        input_df = pd.read_csv(input_path, sep="\t")
        input_summary = input_df[['sample', 'taxid', 'reads', 'mutation_rate']].drop_duplicates()
        input_summary = input_summary.merge(self.input_tax_df, on='taxid', how='left')
        
        return overlap_manager, input_summary
    
    def _compute_baseline_metrics(
        self,
        overlap_manager: OverlapManager,
        input_summary: pd.DataFrame,
        result: DatasetResult
    ) -> DatasetResult:
        """
        Compute pre-prediction baseline metrics.
        
        Args:
            overlap_manager: Loaded OverlapManager
            input_summary: Input summary DataFrame
            result: Result object to populate
            
        Returns:
            Updated DatasetResult
        """
        m_stats = get_m_stats_matrix(
            result.data_set,
            self.config.study_output_filepath,
            self.ncbi,
            overlap_manager
        )
        m_stats = compute_trash_flags(m_stats)
        
        result.output_raw = len(m_stats)
        result.output_cov_filtered = int(len(m_stats[m_stats['coverage'] > 0]))
        
        clean = m_stats[m_stats['is_trash'] == False]
        
        fuzzy_raw, fuzzy_cov = compute_fuzzy_precision(m_stats)
        overall_raw = compute_overall_precision(m_stats)
        recall_raw, recall_cov, _, _ = compute_recall(clean, input_summary)
        
        result.precision = PrecisionMetrics(
            fuzzy_precision_raw=fuzzy_raw,
            fuzzy_precision_cov_filtered=fuzzy_cov,
            overall_precision_raw=overall_raw,
        )
        
        result.recall = RecallMetrics(
            recall_raw=recall_raw,
            recall_cov_filtered=recall_cov,
        )
        
        try:
            trash_comp = get_trash_composition(
                input_summary, m_stats, self.input_tax_df,
                tax_level=self.config.tax_level
            )
            result.trash_composition = trash_comp.to_dict(orient='records') if not trash_comp.empty else None
        except Exception as e:
            logger.warning(f"Failed to compute trash composition: {e}")
            result.trash_composition = None
        
        try:
            cross_hit_comp = get_cross_hit_composition(
                input_summary, m_stats, self.input_tax_df,
                tax_level=self.config.tax_level
            )
            result.cross_hit_composition = cross_hit_comp.to_dict(orient='records') if not cross_hit_comp.empty else None
        except Exception as e:
            logger.warning(f"Failed to compute cross-hit composition: {e}")
            result.cross_hit_composition = None
        
        return result
    
    def _apply_recall_filter(
        self,
        data_set_name: str,
        overlap_manager: OverlapManager,
        result: DatasetResult
    ) -> DatasetResult:
        """
        Apply recall prediction model to filter leaves.
        
        Args:
            data_set_name: Name of dataset
            overlap_manager: OverlapManager instance
            result: Result object to populate
            
        Returns:
            Updated DatasetResult
        """
        m_stats = get_m_stats_matrix(
            data_set_name,
            self.config.study_output_filepath,
            self.ncbi,
            overlap_manager
        )
        
        try:
            filtered_om = cut_off_recall_prediction(
                self.config.study_output_filepath,
                data_set_name,
                self.models.recall_modeller,
                self.config.data_set_divide,
                m_stats,
                self.taxids_to_use
            )
        except Exception as e:
            logger.warning(f"Recall filter failed: {e}")

            return result
        
        new_m_stats = get_m_stats_matrix(
            data_set_name,
            self.config.study_output_filepath,
            self.ncbi,
            filtered_om
        )
        new_m_stats = compute_trash_flags(new_m_stats)
        clean = new_m_stats[new_m_stats['is_trash'] == False]
        
        input_summary_temp = pd.DataFrame({'taxid': self.input_tax_df['taxid'].unique()})
        recall_raw, recall_cov, _, _ = compute_recall(clean, input_summary_temp)
        
        result.recall.recall_filtered_leaves = recall_raw
        
        return result
    
    def _predict_clades_precleanup(
        self,
        data_set_name: str,
        overlap_manager: OverlapManager,
        result: DatasetResult
    ) -> DatasetResult:
        """
        Run clade prediction before cross-hit cleanup.
        
        Args:
            data_set_name: Name of dataset
            overlap_manager: OverlapManager instance
            result: Result object to populate
            
        Returns:
            Updated DatasetResult
        """
        m_stats = get_m_stats_matrix(
            data_set_name,
            self.config.study_output_filepath,
            self.ncbi,
            overlap_manager
        )
        
        try:
            results_df = predict_data_set_clades(
                data_set_name,
                m_stats,
                overlap_manager,
                self.models.composition_modeller,
                self.taxids_to_use,
                tax_level=self.config.tax_level
            )
        except Exception as e:
            logger.warning(f"Clade prediction pre-cleanup failed: {e}")
            import traceback
            traceback.print_exc()
            return result
        
        precision = calculate_overall_precision(results_df, m_stats)
        
        result.predicted_clades_pre = len(results_df)
        result.precision.clade_precision_full = precision
        
        input_summary_temp = pd.DataFrame({'taxid': self.input_tax_df['taxid'].unique()})
        result.recall.clade_recall = compute_clade_recall(results_df, input_summary_temp)
        
        return result
    
    def _apply_crosshit_cleanup(
        self,
        data_set_name: str,
        overlap_manager: OverlapManager,
        result: DatasetResult
    ) -> DatasetResult:
        """
        Predict cross-hits and remove them from the tree.
        
        Args:
            data_set_name: Name of dataset
            overlap_manager: OverlapManager instance
            result: Result object to populate
            
        Returns:
            Updated DatasetResult
        """
        try:
            predictions = cross_hit_prediction(
                data_set_name,
                self.config.study_output_filepath,
                self.ncbi,
                self.models.crosshit_modeller,
                overlap_manager,
                self.taxids_to_use,
                tax_level=self.config.tax_level
            )
        except Exception as e:
            logger.warning(f"Cross-hit prediction failed: {e}")
            return result
        
        filtered = predictions[predictions['prob_best_match'] > self.config.cross_hit_threshold]
        
        result.cross_hit.predicted_cross_hits = len(filtered)
        
        if len(filtered) > 0:
            result.cross_hit.cleanup_accuracy = safe_divide(
                sum(filtered['is_trash']),
                len(filtered)
            )
        
        if filtered.empty or len(overlap_manager.leaves) == 0:
            return result
        
        try:
            cross_hits = [
                leaf for leaf in filtered['leaf'].tolist() 
                if leaf in overlap_manager.leaves
            ]
            overlap_manager.m_stats_matrix.loc[
                overlap_manager.m_stats_matrix.index.isin(cross_hits),
                'numreads'
            ] = 0
            overlap_manager.m_stats_matrix.loc[
                overlap_manager.m_stats_matrix.index.isin(cross_hits),
                'coverage'
            ] = 0
            
            if overlap_manager.m_stats_matrix.shape[0] > 0:
                overlap_manager.prune_empty_nodes()
                overlap_manager.new_tree_from_distance_matrix()
                overlap_manager.recalculate_all_min_pairwise_dist()
        except Exception as e:
            logger.warning(f"Cross-hit cleanup failed: {e}")
        
        return result
    
    def _predict_clades_postcleanup(
        self,
        data_set_name: str,
        overlap_manager: OverlapManager,
        result: DatasetResult
    ) -> DatasetResult:
        """
        Final clade prediction after all cleanup.
        
        Args:
            data_set_name: Name of dataset
            overlap_manager: OverlapManager instance
            result: Result object to populate
            
        Returns:
            Updated DatasetResult
        """
        if overlap_manager.m_stats_matrix.shape[0] < 2:
            logger.warning(f"Not enough leaves after cleanup for {data_set_name}")
            result.predicted_clades_post = 0
            result.precision.clade_precision_post = 0.0
            return result
        
        m_stats = get_m_stats_matrix(
            data_set_name,
            self.config.study_output_filepath,
            self.ncbi,
            overlap_manager
        )
        
        try:
            result_df = predict_data_set_clades(
                data_set_name,
                m_stats,
                overlap_manager,
                self.models.composition_modeller,
                self.taxids_to_use,
                tax_level=self.config.tax_level
            )
        except Exception as e:
            logger.warning(f"Clade prediction post-cleanup failed: {e}")
            result.predicted_clades_post = 0
            result.precision.clade_precision_post = 0.0
            return result
        
        if result_df.empty:
            result.predicted_clades_post = 0
            result.precision.clade_precision_post = 0.0
            return result
        
        precision = calculate_overall_precision(result_df, m_stats)
        
        result.predicted_clades_post = len(result_df)
        result.precision.clade_precision_post = precision
        
        return result
