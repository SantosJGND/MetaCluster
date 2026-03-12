"""
Batch evaluator class for evaluator module.

Manages evaluation across multiple test datasets.
"""

import logging
from typing import List, Optional

import pandas as pd

from tqdm import tqdm

from .config import EvaluatorConfig, TrainedModels
from .result_models import BatchEvaluationResult, DatasetResult, create_empty_result
from .exceptions import EvaluatorError, ResultsAggregationError
from .dataset_processor import DatasetProcessor


logger = logging.getLogger(__name__)


class BatchEvaluator:
    """
    Manages evaluation across multiple test datasets.
    
    This class orchestrates the evaluation of multiple datasets,
    handles error collection, and aggregates results.
    """
    
    def __init__(
        self,
        config: EvaluatorConfig,
        models: TrainedModels,
        ncbi_wrapper: any,
        input_tax_df: pd.DataFrame,
        taxids_to_use: pd.DataFrame,
    ):
        """
        Initialize batch evaluator.
        
        Args:
            config: EvaluatorConfig instance
            models: TrainedModels instance
            ncbi_wrapper: NCBI TaxonomistWrapper instance
            input_tax_df: Input taxid DataFrame
            taxids_to_use: Taxids to use for evaluation
        """
        self.config = config
        self.models = models
        self.ncbi = ncbi_wrapper
        self.input_tax_df = input_tax_df
        self.taxids_to_use = taxids_to_use
        
        self.processor = DatasetProcessor(
            config=config,
            models=models,
            ncbi_wrapper=ncbi_wrapper,
            input_tax_df=input_tax_df,
            taxids_to_use=taxids_to_use,
        )
    
    def evaluate(
        self, 
        dataset_names: List[str], 
        progress: bool = True
    ) -> BatchEvaluationResult:
        """
        Evaluate all datasets and aggregate results.
        
        Args:
            dataset_names: List of dataset names to evaluate
            progress: Whether to show progress bar
            
        Returns:
            BatchEvaluationResult with aggregated results
        """
        logger.info(f"Starting evaluation of {len(dataset_names)} datasets")
        
        results: List[DatasetResult] = []
        errors: List[EvaluatorError] = []
        
        iterator = tqdm(dataset_names, desc="Evaluating datasets") if progress else dataset_names
        
        for name in iterator:
            try:
                result = self.processor.process(name)
                if result is not None:
                    results.append(result)
            except EvaluatorError as e:
                errors.append(e)
                logger.warning(f"Error processing {name}: {e}")
            except Exception as e:
                logger.error(f"Unexpected error processing {name}: {e}")
                errors.append(EvaluatorError(f"Unexpected error: {e}"))
        
        if errors:
            logger.warning(f"Failed to process {len(errors)} out of {len(dataset_names)} datasets")
        
        return self._aggregate_results(results, errors)
    
    def _aggregate_results(
        self, 
        results: List[DatasetResult],
        errors: List[EvaluatorError]
    ) -> BatchEvaluationResult:
        """
        Convert list of DatasetResult to BatchEvaluationResult.
        
        Args:
            results: List of successful results
            errors: List of errors encountered
            
        Returns:
            Aggregated BatchEvaluationResult
        """
        if not results:
            logger.warning("No successful results to aggregate")
            result = create_empty_result()
            result.metadata = {
                'total_datasets': 0,
                'successful': 0,
                'failed': len(errors),
                'errors': [str(e) for e in errors]
            }
            return result
        
        try:
            summary_dfs = []
            test_records = []
            trash_dfs = []
            cross_hit_dfs = []
            
            for r in results:
                summary_row = {
                    'data_set': r.data_set,
                    'sample': r.sample,
                    'input_taxid_count': r.input_taxid_count,
                    'output_raw': r.output_raw,
                    'output_cov_filtered': r.output_cov_filtered,
                    'predicted_clades_pre': r.predicted_clades_pre,
                    'predicted_clades_post': r.predicted_clades_post,
                    'raw_pred_accuracy': 0,
                    'fuzzy_precision_raw': r.precision.fuzzy_precision_raw,
                    'fuzzy_precision_cov_filtered': r.precision.fuzzy_precision_cov_filtered,
                    'overall_precision_raw': r.precision.overall_precision_raw,
                    'clade_precision_full': r.precision.clade_precision_full,
                    'clade_precision_post': r.precision.clade_precision_post,
                    'recall_raw': r.recall.recall_raw,
                    'recall_cov_filtered': r.recall.recall_cov_filtered,
                    'clade_recall': r.recall.clade_recall,
                    'recall_filtered_leaves': r.recall.recall_filtered_leaves,
                    'predicted_cross_hits': r.cross_hit.predicted_cross_hits,
                    'cleanup_accuracy': r.cross_hit.cleanup_accuracy,
                }
                summary_dfs.append(pd.DataFrame([summary_row]))
                
                test_records.append({
                    'overall_precision': r.precision.clade_precision_post,
                    'data_set': r.data_set
                })
                
                if r.trash_composition is not None:
                    trash_df = pd.DataFrame(r.trash_composition)
                    trash_df['data_set'] = r.data_set
                    trash_dfs.append(trash_df)
                
                if r.cross_hit_composition is not None:
                    cross_hit_df = pd.DataFrame(r.cross_hit_composition)
                    cross_hit_df['data_set'] = r.data_set
                    cross_hit_dfs.append(cross_hit_df)
            
            summary_results = pd.concat(summary_dfs, ignore_index=True)
            test_results = pd.DataFrame(test_records)
            trash_composition = pd.concat(trash_dfs, ignore_index=True) if trash_dfs else pd.DataFrame()
            cross_hit_composition = pd.concat(cross_hit_dfs, ignore_index=True) if cross_hit_dfs else pd.DataFrame()
            
            result = BatchEvaluationResult(
                test_results=test_results,
                summary_results=summary_results,
                trash_composition=trash_composition,
                cross_hit_composition=cross_hit_composition,
                metadata={
                    'total_datasets': len(results) + len(errors),
                    'successful': len(results),
                    'failed': len(errors),
                    'errors': [str(e) for e in errors]
                }
            )
            
            logger.info(f"Aggregated results from {len(results)} datasets")
            return result
            
        except Exception as e:
            raise ResultsAggregationError(str(e), len(results)) from e
    
    def save_results(
        self, 
        result: BatchEvaluationResult, 
        output_dir: str,
        save_json: bool = True,
        save_tsv: bool = True
    ) -> None:
        """
        Save results to files.
        
        Args:
            result: BatchEvaluationResult to save
            output_dir: Output directory path
            save_json: Whether to save JSON
            save_tsv: Whether to save TSV
        """
        import os
        
        os.makedirs(output_dir, exist_ok=True)
        
        if save_json:
            json_path = os.path.join(output_dir, "evaluation_results.json")
            result.to_json(json_path)
            logger.info(f"Saved JSON results to {json_path}")
        
        if save_tsv:
            result.save_tsv(output_dir)
            logger.info(f"Saved TSV results to {output_dir}")
    
    def save_summary_statistics(
        self,
        result: BatchEvaluationResult,
        output_dir: str
    ) -> None:
        """
        Save summary statistics to file.
        
        Args:
            result: BatchEvaluationResult
            output_dir: Output directory path
        """
        import os
        
        precision_cols = [
            'overall_precision_raw', 'fuzzy_precision_raw', 'fuzzy_precision_cov_filtered',
            'clade_precision_full', 'clade_precision_post'
        ]
        
        available_cols = [c for c in precision_cols if c in result.summary_results.columns]
        
        if not available_cols:
            logger.warning("No precision columns found for summary statistics")
            return
        
        stats = result.summary_results[available_cols].describe().T
        stats_path = os.path.join(output_dir, "precision_summary_statistics.tsv")
        stats.to_csv(stats_path, sep="\t")
        logger.info(f"Saved summary statistics to {stats_path}")
