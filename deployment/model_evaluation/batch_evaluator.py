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
        taxids_to_use: Optional[pd.DataFrame] = None,
    ):
        """
        Initialize batch evaluator.
        
        Args:
            config: EvaluatorConfig instance
            models: TrainedModels instance
            ncbi_wrapper: NCBI TaxonomistWrapper instance
            input_tax_df: Input taxid DataFrame
            taxids_to_use: Taxids to use for evaluation (optional, will use models.taxids_to_use if available)
        """
        self.config = config
        self.models = models
        self.ncbi = ncbi_wrapper
        self.input_tax_df = input_tax_df
        
        self.taxids_to_use = self._resolve_taxids_to_use(taxids_to_use, models)
        
        self.processor = DatasetProcessor(
            config=config,
            models=models,
            ncbi_wrapper=ncbi_wrapper,
            input_tax_df=input_tax_df,
            taxids_to_use=self.taxids_to_use,
        )
    
    def _resolve_taxids_to_use(
        self, 
        provided_taxids: Optional[pd.DataFrame], 
        models: TrainedModels
    ) -> pd.DataFrame:
        """
        Resolve taxids_to_use with validation.
        
        Ensures that the taxids_to_use used for prediction matches the one used for training.
        
        Args:
            provided_taxids: Taxids provided at prediction time
            models: TrainedModels instance (may contain taxids_to_use from training)
        
        Returns:
            Validated taxids_to_use DataFrame
            
        Raises:
            ValueError: If provided taxids don't match model taxids
        """
        model_taxids = models.taxids_to_use
        
        if model_taxids is None and provided_taxids is None:
            raise ValueError(
                "taxids_to_use is required. Either pass it explicitly or ensure models.taxids_to_use is set."
            )
        
        if model_taxids is not None and provided_taxids is None:
            logger.info("Using taxids_to_use from trained models (none provided)")
            return model_taxids
        
        if model_taxids is None and provided_taxids is not None:
            logger.warning(
                "taxids_to_use provided but models.taxids_to_use is not set. "
                "Cannot validate consistency with training data."
            )
            return provided_taxids
        
        if not self._taxids_equal(model_taxids, provided_taxids):
            provided_count = len(provided_taxids) if provided_taxids is not None else 0
            model_count = len(model_taxids) if model_taxids is not None else 0
            raise ValueError(
                f"taxids_to_use mismatch! "
                f"Provided {provided_count} taxids but model was trained with {model_count} taxids. "
                f"Ensure the same taxids_to_use is used for both training and prediction. "
                f"This may indicate inconsistent study_output_filepath, tax_level, or taxa_threshold parameters."
            )
        
        logger.info("taxids_to_use validated against trained model")
        assert provided_taxids is not None
        return provided_taxids
    
    def _taxids_equal(self, df1: Optional[pd.DataFrame], df2: Optional[pd.DataFrame]) -> bool:
        """
        Check if two taxids_to_use DataFrames are equivalent.
        
        Compares by taxid column values only, ignoring order.
        
        Args:
            df1: First DataFrame
            df2: Second DataFrame
            
        Returns:
            True if DataFrames contain the same taxids
        """
        if df1 is None or df2 is None:
            return False
        
        if set(df1.columns) != set(df2.columns):
            return False
        
        taxid_col = 'taxid'
        if taxid_col not in df1.columns or taxid_col not in df2.columns:
            return False
        
        return set(df1[taxid_col]) == set(df2[taxid_col])
    
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
                import traceback
                traceback.print_exc()
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
            input_dfs = []
            test_records = []
            spurious_dfs = []
            cross_hit_dfs = []
            for r in results:
                summary_row = {
                    'data_set': r.data_set,
                    'sample': r.sample,
                    'input_taxid_count': r.input_taxid_count,
                    'input_read_counts': r.input_read_counts,
                    'reads_simulated_per_class': r.reads_simulated_per_class,
                    'output_raw': r.output_raw,
                    'output_cov_filtered': r.output_cov_filtered,
                    'predicted_clades_pre': r.predicted_clades_pre,
                    'predicted_clades_post': r.predicted_clades_post,
                    'predicted_clades_fixed': r.predicted_clades_fixed,
                    'raw_pred_accuracy': 0,
                    'purity': r.precision.purity_raw,
                    'purity_cov_filtered': r.precision.purity_cov_filtered,
                    'precision_best_match': r.precision.precision_best_match,
                    'precision_clade_full': r.precision.clade_precision_full,
                    'precision_clade_post_cleanup': r.precision.clade_precision_post,
                    'precision_clade_fixed': r.precision.clade_precision_fixed,
                    'recall_baseline': r.recall.recall_raw,
                    'recall_baseline_cov_filtered': r.recall.recall_cov_filtered,
                    'recall_clade_pre_cleanup': r.recall.clade_recall_pre_cleanup,
                    'recall_clade_post_cleanup': r.recall.clade_recall_post_cleanup,
                    'recall_after_recall_filter': r.recall.recall_filtered_leaves,
                    'recall_fixed_max_12': r.recall.recall_fixed_filter,
                    'predicted_cross_hits': r.cross_hit.predicted_cross_hits,
                    'cross_hit_specificity': r.cross_hit.cross_hit_specificity,
                    'cross_hit_precision': r.cross_hit.cross_hit_precision,
                    'cross_hit_recall': r.cross_hit.cross_hit_recall,
                    'cross_hit_f1': r.cross_hit.cross_hit_f1,
                    'total_true_cross_hits': r.cross_hit.total_true_cross_hits,
                    'total_cross_hit_reads_mapped': r.cross_hit.total_cross_hit_reads_mapped,
                    'cross_hit_counts_per_class': r.cross_hit.cross_hit_counts_per_class,
                    'cross_hit_reads_per_class': r.cross_hit.cross_hit_reads_per_class,
                }

                if r.recall.recall_metrics:
                    for k, v in r.recall.recall_metrics.items():
                        if isinstance(v, (int, float)):
                            summary_row[f'recall_metric_{k}'] = v
                        elif isinstance(v, dict) and k == 'per_division_recall_errors':
                            for div_k, div_v in v.items():
                                summary_row[f'recall_metric_error_{div_k}'] = div_v

                if r.input_df is not None:
                    input_df = pd.DataFrame(r.input_df)
                    input_df['data_set'] = r.data_set
                    input_dfs.append(input_df)

                summary_dfs.append(pd.DataFrame([summary_row]))
                
                test_records.append({
                    'recall_baseline': r.recall.recall_raw,
                    'precision_fixed': r.precision.clade_precision_fixed,
                    'precision_clade_post': r.precision.clade_precision_post,
                    'data_set': r.data_set
                })


                if r.spurious_composition is not None:
                    spurious_df = pd.DataFrame(r.spurious_composition)
                    spurious_df['data_set'] = r.data_set
                    spurious_dfs.append(spurious_df)
                
                if r.cross_hit_composition is not None:
                    cross_hit_df = pd.DataFrame(r.cross_hit_composition)
                    cross_hit_df['data_set'] = r.data_set
                    cross_hit_dfs.append(cross_hit_df)

            input_df = pd.concat(input_dfs, ignore_index=True) if input_dfs else pd.DataFrame()
            summary_results = pd.concat(summary_dfs, ignore_index=True)
            test_results = pd.DataFrame(test_records)
            spurious_composition = pd.concat(spurious_dfs, ignore_index=True) if spurious_dfs else pd.DataFrame()
            cross_hit_composition = pd.concat(cross_hit_dfs, ignore_index=True) if cross_hit_dfs else pd.DataFrame()
            
            result = BatchEvaluationResult(
                input_df=input_df,
                test_results=test_results,
                summary_results=summary_results,
                spurious_composition=spurious_composition,
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
        
        agent_path = os.path.join(output_dir, "evaluation_results_agent.json")
        result.write_agent_output(agent_path)
        logger.info(f"Saved agent-parseable JSON to {agent_path}")
        
        if save_tsv:
            result.save_tsv(output_dir)
            self.taxids_to_use.to_csv(os.path.join(output_dir, "taxids_to_use.tsv"), sep="\t", index=False)
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
            'precision_best_match', 'purity', 'purity_cov_filtered',
            'precision_clade_full', 'precision_clade_post_cleanup'
        ]
        
        available_cols = [c for c in precision_cols if c in result.summary_results.columns]
        
        if not available_cols:
            logger.warning("No precision columns found for summary statistics")
        else:
            stats = result.summary_results[available_cols].describe().T
            stats_path = os.path.join(output_dir, "precision_summary_statistics.tsv")
            stats.to_csv(stats_path, sep="\t")
            logger.info(f"Saved summary statistics to {stats_path}")
        
        cross_hit_cols = [
            'cross_hit_precision', 'cross_hit_recall', 'cross_hit_f1',
            'cross_hit_specificity', 'predicted_cross_hits', 'total_true_cross_hits'
        ]
        
        available_cross_hit_cols = [c for c in cross_hit_cols if c in result.summary_results.columns]
        
        if not available_cross_hit_cols:
            logger.info("No cross-hit columns found for summary statistics (cross-hit may be disabled)")
        else:
            cross_hit_stats = result.summary_results[available_cross_hit_cols].describe().T
            cross_hit_stats_path = os.path.join(output_dir, "cross_hit_summary_statistics.tsv")
            cross_hit_stats.to_csv(cross_hit_stats_path, sep="\t")
            logger.info(f"Saved cross-hit summary statistics to {cross_hit_stats_path}")
