"""
Dataset processor class for evaluator module.

Processes a single dataset through the full evaluation pipeline.
"""

import logging
import os
from typing import Any

import numpy as np
import pandas as pd

from metagenomics_utils.overlap_manager import OverlapManager
from metagenomics_utils.overlap_manager.node_stats import get_m_stats_matrix
from metagenomics_utils.overlap_manager.om_models import (
    calculate_clade_precision,
    cross_hit_prediction,
    cut_off_recall_prediction,
    get_cross_hit_composition,
    get_spurious_composition,
    predict_data_set_clades_composition,
    predict_data_set_clades_fixed,
)

from .config import EvaluatorConfig, TrainedModels
from .exceptions import DataLoadError, OverlapManagerError, PredictionError
from .metrics import (
    compute_clade_recall,
    compute_mstats_precision,
    compute_purity,
    compute_recall,
    safe_divide,
)
from .result_models import (
    DatasetResult,
    PrecisionMetrics,
    RecallMetrics,
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

    def process(self, data_set_name: str) -> DatasetResult | None:
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
                input_df=input_summary,
                sample=input_summary["sample"].iloc[0] if len(input_summary) > 0 else "",
                input_taxid_count=int(input_summary["taxid"].nunique()),
            )

            result = self._compute_baseline_metrics(overlap_manager, result)
            result = self._predict_clades_precleanup(data_set_name, overlap_manager, result)
            result, filtered_om = self._apply_recall_filter(data_set_name, overlap_manager, result=result)
            result, _ = self._apply_fixed_filter(data_set_name, result=result, max_taxids=12)

            if self.config.enable_cross_hit:
                result = self._apply_crosshit_cleanup(data_set_name, overlap_manager, result)
            else:
                logger.info("Cross-hit cleanup disabled, skipping")

            result = self._predict_clades_postcleanup(data_set_name, filtered_om, result)

            logger.info(f"Completed processing {data_set_name}")
            return result

        except FileNotFoundError as e:
            raise DataLoadError(data_set_name, str(e))
        except Exception as e:
            raise PredictionError(data_set_name, "processing", str(e)) from e

    def _load_data(self, data_set_name: str) -> tuple[OverlapManager | None, pd.DataFrame]:
        """
        Load OverlapManager and input data.

        Args:
            data_set_name: Name of dataset

        Returns:
            Tuple of (OverlapManager, input_summary) or (None, empty DataFrame)
        """
        om_path = os.path.join(self.config.study_output_filepath, data_set_name, "clustering")
        input_path = os.path.join(self.config.study_output_filepath, data_set_name, "input", f"{data_set_name}.tsv")

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
        input_summary = input_df[["sample", "taxid", "reads", "mutation_rate"]].drop_duplicates()
        input_summary = input_summary.merge(self.input_tax_df, on="taxid", how="left")

        return overlap_manager, input_summary

    def _compute_baseline_metrics(self, overlap_manager: OverlapManager, result: DatasetResult) -> DatasetResult:
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
            result.data_set, self.config.study_output_filepath, self.ncbi, overlap_manager, filter_no_leaf=False
        )

        taxid_match = m_stats[m_stats["best_match_is_best"] == True]

        def match_taxid(row: pd.Series) -> pd.Series:
            taxid = row["taxid"]
            if taxid in taxid_match["best_match_taxid"].values:
                match_index = taxid_match[taxid_match["best_match_taxid"] == taxid].index[0]
                match_cov = taxid_match.loc[match_index, "coverage"]
                row["match_index"] = match_index
                row["match_coverage"] = match_cov
            else:
                row["match_index"] = np.nan
                row["match_coverage"] = np.nan
            return row

        result.input_df = result.input_df.apply(match_taxid, axis=1)

        result.output_raw = len(overlap_manager.original_m_stats_matrix)
        result.output_taxid_count = int(len(m_stats["taxid"].dropna().unique()))
        result.output_cov_filtered = int(len(m_stats[m_stats["coverage"] > 0]))

        cross_hit_mask = m_stats["is_crosshit"] == True
        cross_hit_subset = m_stats[cross_hit_mask]
        cross_hit_subset = cross_hit_subset[["best_match_taxid", "numreads"]].merge(
            result.input_df, left_on="best_match_taxid", right_on="taxid", how="left"
        )
        result.cross_hit.total_true_cross_hits = int(cross_hit_subset.shape[0])
        result.cross_hit.total_cross_hit_reads_mapped = int(cross_hit_subset["numreads"].sum())

        tax_level = self.config.tax_level
        if not cross_hit_subset.empty and tax_level in cross_hit_subset.columns:
            cross_hit_by_class = cross_hit_subset.groupby(tax_level)["best_match_taxid"].count().reset_index()
            cross_hit_by_class.columns = ["class", "count"]
            result.cross_hit.cross_hit_counts_per_class = cross_hit_by_class.to_dict(orient="records")

            cross_hit_reads_by_class = cross_hit_subset.groupby(tax_level)["numreads"].sum().reset_index()
            cross_hit_reads_by_class.columns = ["class", "reads"]
            result.cross_hit.cross_hit_reads_per_class = cross_hit_reads_by_class.to_dict(orient="records")

        spurious_mask = m_stats["is_trash"] == True
        spurious_subset = m_stats[spurious_mask]
        spurious_subset = spurious_subset[["best_match_taxid", "numreads"]].merge(
            result.input_df, left_on="best_match_taxid", right_on="taxid", how="left"
        )
        if not spurious_subset.empty and tax_level in spurious_subset.columns:
            spurious_by_class = spurious_subset.groupby(tax_level)["best_match_taxid"].count().reset_index()
            spurious_by_class.columns = ["class", "count"]
            result.cross_hit.spurious_hit_counts_per_class = spurious_by_class.to_dict(orient="records")

            spurious_reads_by_class = spurious_subset.groupby(tax_level)["numreads"].sum().reset_index()
            spurious_reads_by_class.columns = ["class", "reads"]
            result.cross_hit.spurious_hit_reads_per_class = spurious_reads_by_class.to_dict(orient="records")

        fuzzy_raw, raw_best_match, fuzzy_cov = compute_purity(m_stats)
        overall_raw = compute_mstats_precision(m_stats, result.input_df)
        recall_raw, recall_cov, _, _ = compute_recall(m_stats, result.input_df)

        result.recall = RecallMetrics(
            recall_raw=recall_raw,
            recall_cov_filtered=recall_cov,
        )

        result.precision = PrecisionMetrics(
            purity_raw=fuzzy_raw,
            purity_cov_filtered=fuzzy_cov,
            precision_best_match=raw_best_match,
        )

        try:
            spurious_comp = get_spurious_composition(
                result.input_df, m_stats, self.taxids_to_use, tax_level=self.config.tax_level
            )
            result.spurious_composition = spurious_comp.to_dict(orient="records") if not spurious_comp.empty else None
        except Exception as e:
            logger.warning(f"Failed to compute spurious composition: {e}")
            result.spurious_composition = None

        try:
            cross_hit_comp = get_cross_hit_composition(
                result.input_df, m_stats, self.taxids_to_use, tax_level=self.config.tax_level
            )
            result.cross_hit_composition = (
                cross_hit_comp.to_dict(orient="records") if not cross_hit_comp.empty else None
            )
        except Exception as e:
            logger.warning(f"Failed to compute cross-hit composition: {e}")
            result.cross_hit_composition = None

        try:
            from metagenomics_utils.overlap_manager.om_models import get_subset_composition

            input_composition = (
                get_subset_composition(result.input_df, self.taxids_to_use, tax_level=self.config.tax_level)
                .set_index("tax_level")
                .T
            )
            input_composition = input_composition.reset_index(drop=True)
            input_composition.loc[:, "tax_level"] = self.config.tax_level
            input_composition.insert(0, "data_set", result.data_set)
            result.input_composition = input_composition.to_dict(orient="records")

        except Exception as e:
            logger.warning(f"Failed to compute input composition: {e}")
            result.input_composition = None

        try:
            from metagenomics_utils.overlap_manager.node_stats import get_subset_composition_counts

            input_read_counts = get_subset_composition_counts(
                result.input_df, self.taxids_to_use, tax_level=self.config.tax_level, count_column="reads"
            )
            input_read_counts = input_read_counts[["tax_level", "total_uniq_reads"]].rename(
                columns={"tax_level": "class", "total_uniq_reads": "reads"}
            )
            result.input_read_counts = input_read_counts.to_dict(orient="records")

        except Exception as e:
            logger.warning(f"Failed to compute input read counts: {e}")
            result.input_read_counts = None

        return result

    def _apply_fixed_filter(
        self, data_set_name: str, result: DatasetResult, max_taxids: int
    ) -> tuple[DatasetResult, OverlapManager]:

        overlap_manager = OverlapManager(
            os.path.join(self.config.study_output_filepath, f"{data_set_name}", "clustering"), max_taxids=max_taxids
        )

        new_m_stats = get_m_stats_matrix(
            data_set_name, self.config.study_output_filepath, self.ncbi, overlap_manager, filter_no_leaf=False
        )
        new_m_stats = new_m_stats.head(max_taxids)  # Ensure we only take the top max_taxids after filtering
        new_m_stats = new_m_stats[new_m_stats["coverage"] > 0]  # Filter out zero coverage leaves for recall calculation

        logger.debug(f"######### FIXED FILTER M-STATS {data_set_name} #########")
        logger.debug(f"{new_m_stats.shape}, {max_taxids}")
        recall_raw, recall_cov, kept_taxids, kept_taxids_cov = compute_recall(new_m_stats, result.input_df)
        print("#########################################")

        result.input_df["found_in_fixed_filter"] = result.input_df["taxid"].isin(kept_taxids)

        result.recall.recall_fixed_filter = recall_raw

        return result, overlap_manager

    def _apply_recall_filter(
        self, data_set_name: str, overlap_manager: OverlapManager, result: DatasetResult
    ) -> tuple[DatasetResult, OverlapManager]:
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
            data_set_name, self.config.study_output_filepath, self.ncbi, overlap_manager, filter_no_leaf=False
        )

        try:
            filtered_om, metrics_dict = cut_off_recall_prediction(
                self.config.study_output_filepath,
                data_set_name,
                self.models.recall_modeller,
                self.config.data_set_divide,
                m_stats,
                self.taxids_to_use,
                tax_level=self.config.tax_level,
                target_recall=self.config.target_recall,
                confidence=self.config.cutoff_confidence,
            )
        except Exception as e:
            logger.warning(f"Recall filter failed: {e}")
            import traceback

            print(traceback.format_exc())
            raise PredictionError(data_set_name, "recall_filter", str(e)) from e

        new_m_stats = get_m_stats_matrix(
            data_set_name, self.config.study_output_filepath, self.ncbi, filtered_om, filter_no_leaf=False
        )

        new_m_stats = new_m_stats.head(
            metrics_dict.get("keep_index", len(new_m_stats))
        )  # Ensure we only take the top keep_index after filtering
        new_m_stats = new_m_stats[new_m_stats["coverage"] > 0]  # Filter out zero coverage leaves for recall calculation

        logger.debug(f"######### RECALL FILTER M-STATS {data_set_name} #########")
        logger.debug(
            f"{new_m_stats.shape}, {filtered_om.original_m_stats_matrix.shape}, {self.config.target_recall}, {metrics_dict.get('target_percentile', None)}, {metrics_dict.get('keep_index', None)}"
        )
        recall_raw, recall_cov, taxid_found, _taxid_found_cov = compute_recall(new_m_stats, result.input_df)
        print("#########################################")

        result.input_df["found_in_recall_filter"] = result.input_df["taxid"].isin(taxid_found)

        result.recall.recall_filtered_leaves = recall_raw

        result.recall.recall_metrics = metrics_dict

        overlap_manager = filtered_om

        return result, overlap_manager

    def _predict_clades_precleanup(
        self, data_set_name: str, overlap_manager: OverlapManager, result: DatasetResult
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
            data_set_name, self.config.study_output_filepath, self.ncbi, overlap_manager, filter_no_leaf=True
        )

        try:
            results_df = predict_data_set_clades_composition(
                data_set_name,
                m_stats,
                overlap_manager,
                self.models.composition_modeller,
                self.taxids_to_use,
                tax_level=self.config.tax_level,
            )
        except Exception as e:
            logger.warning(f"Clade prediction pre-cleanup failed: {e}")
            import traceback

            traceback.print_exc()
            return result

        precision = calculate_clade_precision(results_df, result.input_df)

        result.predicted_clades_pre = len(results_df)
        result.precision.clade_precision_full = precision

        result.recall.clade_recall_pre_cleanup = compute_clade_recall(results_df, result.input_df)

        return result

    def _apply_crosshit_cleanup(
        self, data_set_name: str, overlap_manager: OverlapManager, result: DatasetResult
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
                tax_level=self.config.tax_level,
            )
        except Exception as e:
            logger.warning(f"Cross-hit prediction failed: {e}")
            return result

        filtered = predictions[predictions["prob_best_match"] > self.config.cross_hit_threshold]

        m_stats = get_m_stats_matrix(data_set_name, self.config.study_output_filepath, self.ncbi, overlap_manager)

        true_cross_mask = (m_stats["best_match_is_best"] == False) | (m_stats["best_match_taxid"].isna())
        total_true_cross_hits = int(true_cross_mask.sum())

        TP = len(filtered[filtered["is_trash"] == True]) if len(filtered) > 0 else 0
        FP = len(filtered[filtered["is_trash"] == False]) if len(filtered) > 0 else 0

        cross_hit_specificity = safe_divide(TP, len(filtered)) if len(filtered) > 0 else 0.0
        cross_hit_precision = cross_hit_specificity
        cross_hit_recall = safe_divide(TP, total_true_cross_hits) if total_true_cross_hits > 0 else 0.0
        cross_hit_f1 = safe_divide(2 * cross_hit_precision * cross_hit_recall, cross_hit_precision + cross_hit_recall)

        result.cross_hit.predicted_cross_hits = len(filtered)
        result.cross_hit.cross_hit_specificity = cross_hit_specificity
        result.cross_hit.cross_hit_precision = cross_hit_precision
        result.cross_hit.cross_hit_recall = cross_hit_recall
        result.cross_hit.cross_hit_f1 = cross_hit_f1
        result.cross_hit.total_true_cross_hits = total_true_cross_hits

        if filtered.empty or len(overlap_manager.leaves) == 0:
            return result

        try:
            cross_hits = [leaf for leaf in filtered["leaf"].tolist() if leaf in overlap_manager.leaves]
            overlap_manager.m_stats_matrix.loc[overlap_manager.m_stats_matrix.index.isin(cross_hits), "numreads"] = 0
            overlap_manager.m_stats_matrix.loc[overlap_manager.m_stats_matrix.index.isin(cross_hits), "coverage"] = 0

            if overlap_manager.m_stats_matrix.shape[0] > 0:
                overlap_manager.prune_empty_nodes()
                overlap_manager.new_tree_from_distance_matrix()
                overlap_manager.recalculate_all_min_pairwise_dist()
        except Exception as e:
            logger.warning(f"Cross-hit cleanup failed: {e}")

        return result

    def _predict_clades_postcleanup(
        self, data_set_name: str, overlap_manager: OverlapManager, result: DatasetResult
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
            data_set_name, self.config.study_output_filepath, self.ncbi, overlap_manager, filter_no_leaf=False
        )

        try:
            fixed_result_df = predict_data_set_clades_fixed(
                data_set_name, m_stats, overlap_manager, min_dist_threshold=0.6
            )
            result.predicted_clades_fixed = len(fixed_result_df)
            result.precision.clade_precision_fixed = calculate_clade_precision(fixed_result_df, result.input_df)

        except Exception as e:
            logger.warning(f"Fixed Clade prediction post-cleanup failed: {e}")
            result.predicted_clades_fixed = 0
            result.precision.clade_precision_fixed = 0.0
            return result

        print(f"######### FIXED CLADES {data_set_name} #########")
        print(f"Predicted clades: {len(fixed_result_df)}, Precision: {result.precision.clade_precision_fixed:.4f}")

        try:
            result_df = predict_data_set_clades_composition(
                data_set_name,
                m_stats,
                overlap_manager,
                self.models.composition_modeller,
                self.taxids_to_use,
                tax_level=self.config.tax_level,
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

        precision = calculate_clade_precision(result_df, result.input_df)
        print(f"######### POST-CLEANUP CLADES {data_set_name} #########")
        print(f"Predicted clades: {len(result_df)}, Precision: {precision:.4f}")

        result.predicted_clades_fixed = len(fixed_result_df)
        result.predicted_clades_post = len(result_df)
        result.precision.clade_precision_post = precision

        result.recall.clade_recall_post_cleanup = compute_clade_recall(result_df, result.input_df)

        return result
