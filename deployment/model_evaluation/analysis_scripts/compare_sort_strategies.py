"""
Compare recall sort strategies across multiple recall model types.

Sweeps sort_strategy x recall_model_interface over cached or freshly
retrieved training data, then evaluates each combination on held-out
test datasets. Reports target_percentile, keep_index, actual recall,
max_possible_recall, and standardized_recall for every combination.

Includes a rule-based baseline (fixed_12_taxids) that keeps leaves until
N unique best_match_taxids are covered after sorting.

Usage:
    python deployment/model_evaluation/analysis_scripts/compare_sort_strategies.py \\
        --study_output_filepath /path/to/study \\
        --taxid_plan_filepath /path/to/taxid_plan.tsv \\
        --analysis_output_filepath /path/to/output/
"""

import argparse
import logging
import os
import re
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from deployment.model_evaluation.config import EvaluatorConfig
from deployment.model_evaluation.data_loader import DataLoader
from deployment.model_evaluation.features import RecallFeatureTransformer
from deployment.model_evaluation.metrics import compute_recall

logger = logging.getLogger(__name__)

# -- Strategies and models to sweep --

DEFAULT_SORT_STRATEGIES = ['reads', 'taxid_roundrobin', 'rarity_boost', 'tax_level_stratified']

DEFAULT_RECALL_MODELS = ['xgb', 'morf', 'direct_xgb', 'gp_clf', 'fixed_12_taxids']


# -- Helpers --

def _get_raw_m_stats_matrix(
    data_set_name: str,
    study_output_filepath: Path,
    ncbi_wrapper,
) -> Optional[pd.DataFrame]:
    """Load the raw m_stats matrix (pre-sorting) for a single dataset."""
    from metagenomics_utils.overlap_manager import OverlapManager
    from metagenomics_utils.overlap_manager.node_stats import get_m_stats_matrix

    om_path = os.path.join(str(study_output_filepath), data_set_name, "clustering")
    try:
        overlap_manager = OverlapManager(om_path)
    except Exception:
        return None

    if overlap_manager.m_stats_matrix.empty:
        return None

    m_stats = get_m_stats_matrix(
        data_set_name,
        str(study_output_filepath),
        ncbi_wrapper,
        overlap_manager,
        filter_no_leaf=False,
    )
    return m_stats if not m_stats.empty else None


def _load_input_summary(
    data_set_name: str,
    study_output_filepath: Path,
) -> Optional[pd.DataFrame]:
    """Load the input summary (taxid column) for a single dataset."""
    input_path = os.path.join(
        str(study_output_filepath), data_set_name, "input", f"{data_set_name}.tsv"
    )
    if not os.path.exists(input_path):
        return None
    input_df = pd.read_csv(input_path, sep="\t")
    return input_df[['sample', 'taxid', 'reads', 'mutation_rate']].drop_duplicates()


def _collect_matrices(
    folder_names: List[str],
    config: EvaluatorConfig,
    ncbi_wrapper,
) -> List[tuple]:
    """Collect raw m_stats matrices for a list of dataset folder names.

    Returns list of (dataset_name, matrix) tuples.
    """
    matrices = []
    for name in folder_names:
        m = _get_raw_m_stats_matrix(name, config.study_output_filepath, ncbi_wrapper)
        if m is not None:
            matrices.append((name, m))
    logger.info("Collected %d raw m_stats matrices from %d folders", len(matrices), len(folder_names))
    return matrices


def _compute_actual_recall(
    sorted_m_stats: pd.DataFrame,
    keep_index: int,
    input_summary: pd.DataFrame,
) -> float:
    """Compute intersection-based recall at a given keep_index.

    Uses the fixed compute_recall from metrics.py: measures what fraction
    of input taxids are recovered in the top keep_index leaves.
    """
    top = sorted_m_stats.iloc[:keep_index].copy()
    top.reset_index(drop=False, inplace=True)
    rec_raw, rec_cov, _, _ = compute_recall(top, input_summary)
    return rec_raw


def _compute_max_possible_recall(
    sorted_m_stats: pd.DataFrame,
    input_summary: pd.DataFrame,
) -> float:
    """Compute recall using ALL leaves (no cutoff) -- the ceiling for standardization."""
    rec_raw, rec_cov, _, _ = compute_recall(sorted_m_stats, input_summary)
    return rec_raw


def _compute_n_taxids_found(
    sorted_m_stats: pd.DataFrame,
    keep_index: int,
    input_summary: pd.DataFrame,
) -> int:
    """Count input taxids recovered in first keep_index leaves (intersection)."""
    input_taxids = set(input_summary['taxid'].unique())
    top = sorted_m_stats.iloc[:keep_index]
    best = top[(top['is_trash'] == False) & (top['best_match_is_best'] == True)]
    output_taxids = set(best['best_match_taxid'].dropna().unique())
    return len(output_taxids & input_taxids)


def _compute_true_cutoff(sorted_m_stats: pd.DataFrame) -> float:
    """Fractional position of the last unique valid best_match_taxid.

    This is the ground-truth cutoff: the minimum percentile of sorted
    data needed to capture every valid (non-trash, best) unique taxid.
    Returns 0.0 if no valid taxids are found.
    """
    seen: set = set()
    last_pos = -1
    for pos, (_, row) in enumerate(sorted_m_stats.iterrows()):
        is_trash = row.get('is_trash', True)
        is_best = row.get('best_match_is_best', False)
        if not is_trash and is_best:
            taxid = row.get('best_match_taxid')
            if pd.notna(taxid) and taxid not in seen:
                seen.add(taxid)
                last_pos = pos
    if last_pos < 0:
        return 0.0
    return (last_pos + 1) / len(sorted_m_stats)


class FixedCountModeller:
    """Baseline: keep leaves until N unique best_match_taxids are covered.

    Predicts the percentile cutoff that captures ``n_taxids`` unique
    best-match taxids after sorting with the configured strategy.
    Requires no training -- purely rule-based.
    """

    def __init__(self, n_taxids: int = 12, feature_transformer=None):
        self.n_taxids = n_taxids
        self.feature_transformer = feature_transformer

    def fit(self, *args, **kwargs):
        pass

    def predict_cutoff(
        self, m_stats: pd.DataFrame, target_recall: float = None,
    ) -> float:
        if self.feature_transformer is not None:
            sorted_stats = self.feature_transformer._apply_sort(m_stats.copy())
        else:
            sorted_stats = m_stats.sort_values('total_uniq_reads', ascending=False)

        seen: set = set()
        if sorted_stats.shape[0] <= self.n_taxids:
            return 1.0
        else:
            return 12 / sorted_stats.shape[0]


def _create_modeller(model_type: str, data_set_divide: int, transformer, target_recall: float):
    """Factory: create a recall modeller with the given transformer injected."""
    from metagenomics_utils.overlap_manager.om_models import (
        RecallModeller,
        CutoffRecallModeller,
        DirectXGBRecallModeller,
        GPCLFRecallModeller,
        InjectModellerInterface,
    )

    if model_type.startswith('fixed_'):
        match = re.match(r'fixed_(\d+)', model_type)
        if not match:
            raise ValueError(f"Invalid fixed-count model type: {model_type}")
        n_taxids = int(match.group(1))
        modeller = FixedCountModeller(
            n_taxids=n_taxids,
            feature_transformer=transformer,
        )
    elif model_type == 'gp_clf':
        model_interface = InjectModellerInterface(model_type='gp_clf')
        modeller = GPCLFRecallModeller(
            data_set_divide=data_set_divide,
            model_interface=model_interface,
            feature_transformer=transformer,
            sort_strategy=transformer.sort_strategy,
        )
    elif model_type == 'direct':
        modeller = CutoffRecallModeller(
            data_set_divide=data_set_divide,
            target_recall=target_recall,
            feature_transformer=transformer,
            sort_strategy=transformer.sort_strategy,
        )
    elif model_type == 'direct_xgb':
        modeller = DirectXGBRecallModeller(
            data_set_divide=data_set_divide,
            target_recall=target_recall,
            feature_transformer=transformer,
            sort_strategy=transformer.sort_strategy,
        )
    else:
        model_interface = InjectModellerInterface(model_type=model_type)
        modeller = RecallModeller(
            data_set_divide=data_set_divide,
            model_interface=model_interface,
            feature_transformer=transformer,
            sort_strategy=transformer.sort_strategy,
        )
    return modeller


# -- Main sweep --

def run_comparison(
    config: EvaluatorConfig,
    sort_strategies: List[str],
    recall_models: List[str],
    target_recall: float,
) -> pd.DataFrame:
    """
    Run sort_strategy x recall_model comparison.

    Returns a DataFrame with one row per (strategy, model, test_dataset).
    """
    from metagenomics_utils.overlap_manager.om_models import InjectModellerInterface
    from sklearn.model_selection import train_test_split

    logger.info("Loading data...")
    loader = DataLoader(config).initialize()
    ncbi = loader.get_ncbi_wrapper()
    input_tax_df = loader.get_input_tax_df()
    taxids_to_use = loader.get_taxids_to_use()
    training_folders = loader.get_training_folders()
    test_folders = loader.get_test_folders()

    logger.info("Collecting training matrices (%d folders) ...", len(training_folders))
    train_matrices = _collect_matrices(training_folders, config, ncbi)

    logger.info("Collecting test matrices (%d folders) ...", len(test_folders))
    test_matrices = _collect_matrices(test_folders, config, ncbi)

    if not train_matrices:
        logger.error("No training matrices found -- cannot train models.")
        return pd.DataFrame()

    if not test_matrices:
        logger.error("No test matrices found -- nothing to evaluate.")
        return pd.DataFrame()

    records = []

    for strategy in sort_strategies:
        logger.info("=" * 60)
        logger.info("Strategy: %s", strategy)
        logger.info("=" * 60)

        transformer = RecallFeatureTransformer(
            tax_level=config.tax_level,
            data_set_divide=config.data_set_divide,
            sort_strategy=strategy,
        )
        transformer.fit(taxids_to_use=taxids_to_use)

        train_rows = []
        for _name, m in train_matrices:
            try:
                row = transformer.transform(m)
                train_rows.append(row)
            except Exception as e:
                logger.warning("Transform failed for training matrix (%s): %s", _name, e)

        if not train_rows:
            logger.warning("No training rows for strategy '%s' -- skipping", strategy)
            continue

        train_df = pd.concat(train_rows, ignore_index=True)
        feat_cols = transformer.get_feature_names_out()
        target_cols = transformer.target_columns_
        X = train_df[feat_cols]
        Y = train_df[target_cols]

        raw_train_matrices = [m for _name, m in train_matrices]

        for model_type in recall_models:
            logger.info("  Model: %s", model_type)

            is_fixed = model_type.startswith('fixed_')

            modeller = _create_modeller(
                model_type, config.data_set_divide, transformer, target_recall,
            )

            if is_fixed:
                pass
            else:
                X_tr, X_te, Y_tr, Y_te = train_test_split(
                    X, Y, test_size=0.2, random_state=42,
                )
                modeller.X_test = X_te
                modeller.y_test = Y_te

                try:
                    if model_type == 'gp_clf':
                        modeller.fit(
                            raw_train_matrices,
                            taxids_to_use,
                            test_size=0.2,
                            random_state=42,
                            optimize=True,
                        )
                    elif model_type == 'direct_xgb':
                        modeller.fit(raw_train_matrices, taxids_to_use, test_size=0.2)
                    elif model_type == 'direct':
                        modeller.fit(raw_train_matrices, taxids_to_use, test_size=0.2)
                    else:
                        modeller.model = modeller.model_interface.train_model(X_tr, Y_tr)
                        modeller.RecP_feature_cols = feat_cols
                        modeller.RecP_target_cols = target_cols

                except Exception as e:
                    logger.warning(
                        "Training failed for %s / %s: %s", strategy, model_type, e,
                    )
                    import traceback
                    traceback.print_exc()
                    continue

            for ds_name, test_m in test_matrices:
                try:
                    input_summary = _load_input_summary(ds_name, config.study_output_filepath)
                    if input_summary is None:
                        logger.warning("No input summary for %s -- skipping", ds_name)
                        continue

                    pred_percentile = modeller.predict_cutoff(
                        test_m, target_recall=target_recall,
                    )
                    keep_idx = max(1, round(pred_percentile * len(test_m)))

                    sorted_test = transformer._apply_sort(test_m.copy()).reset_index(drop=True)

                    actual_recall = _compute_actual_recall(sorted_test, keep_idx, input_summary)
                    max_possible_recall = _compute_max_possible_recall(sorted_test, input_summary)
                    standardized_recall = (
                        actual_recall / max_possible_recall
                        if max_possible_recall > 0 else 0.0
                    )
                    standardized_recall = min(standardized_recall, 1.0)
                    true_cutoff = _compute_true_cutoff(sorted_test)
                    n_taxids_found = _compute_n_taxids_found(sorted_test, keep_idx, input_summary)
                    n_input_taxids = int(input_summary['taxid'].nunique())

                    records.append({
                        'sort_strategy': strategy,
                        'recall_model': model_type,
                        'data_set': ds_name,
                        'total_leaves': len(sorted_test),
                        'n_input_taxids': n_input_taxids,
                        'target_percentile': pred_percentile,
                        'keep_index': keep_idx,
                        'actual_recall': actual_recall,
                        'max_possible_recall': max_possible_recall,
                        'standardized_recall': standardized_recall,
                        'recall_ceiling_gap': max_possible_recall - actual_recall,
                        'true_cutoff': true_cutoff,
                        'n_taxids_found': n_taxids_found,
                        'recall_gap': actual_recall - target_recall,
                    })
                except Exception as e:
                    logger.warning(
                        "Prediction failed for %s / %s on %s: %s",
                        strategy, model_type, ds_name, e,
                    )
                    continue

    return pd.DataFrame(records)


# -- Report helpers --

def generate_summary(results_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate results by (sort_strategy, recall_model)."""
    if results_df.empty:
        return pd.DataFrame()

    group_cols = ['sort_strategy', 'recall_model']
    has_standardized = 'standardized_recall' in results_df.columns
    has_max_recall = 'max_possible_recall' in results_df.columns

    agg_dict = {
        'n_datasets': ('target_percentile', 'count'),
        'mean_percentile': ('target_percentile', 'mean'),
        'std_percentile': ('target_percentile', 'std'),
        'median_percentile': ('target_percentile', 'median'),
        'mean_keep_index': ('keep_index', 'mean'),
        'std_keep_index': ('keep_index', 'std'),
        'median_keep_index': ('keep_index', 'median'),
        'mean_actual_recall': ('actual_recall', 'mean'),
        'std_actual_recall': ('actual_recall', 'std'),
        'mean_recall_gap': ('recall_gap', 'mean'),
        'mean_taxids_found': ('n_taxids_found', 'mean'),
        'mean_input_taxids': ('n_input_taxids', 'mean'),
    }
    if has_standardized:
        agg_dict['mean_standardized_recall'] = ('standardized_recall', 'mean')
        agg_dict['std_standardized_recall'] = ('standardized_recall', 'std')
    if has_max_recall:
        agg_dict['mean_max_possible_recall'] = ('max_possible_recall', 'mean')

    summary = results_df.groupby(group_cols).agg(**agg_dict).reset_index()

    summary = summary.sort_values(['sort_strategy', 'mean_percentile'])
    return summary


def _density_grid(results_df, x_col, y_col, x_label, y_label, suptitle, filename, output_dir,
                  models=None, strategies=None):
    """2D density + scatter grid, rows=recall_model, cols=sort_strategy."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    if models is None:
        models = sorted(results_df['recall_model'].unique())
    if strategies is None:
        strategies = sorted(results_df['sort_strategy'].unique())

    n_models, n_strats = len(models), len(strategies)
    fig, axes = plt.subplots(n_models, n_strats,
                             figsize=(4.5 * n_strats, 4 * n_models),
                             squeeze=False)

    for mi, model in enumerate(models):
        for si, strat in enumerate(strategies):
            ax = axes[mi][si]
            sub = results_df[
                (results_df['recall_model'] == model) &
                (results_df['sort_strategy'] == strat)
            ]

            if len(sub) >= 3:
                sns.kdeplot(data=sub, x=x_col, y=y_col, fill=True,
                            thresh=0.05, alpha=0.35, ax=ax)
            sns.scatterplot(data=sub, x=x_col, y=y_col, color='black',
                            alpha=0.15, s=12, linewidth=0, ax=ax)

            ax.plot([0, 1], [0, 1], '--', color='grey', alpha=0.4, linewidth=0.8)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)

            if mi == 0:
                ax.set_title(strat, fontsize=10, fontweight='bold')
            if si == 0:
                ax.set_ylabel(f'{model}\n{y_label}', fontsize=9)
            if mi == n_models - 1:
                ax.set_xlabel(x_label, fontsize=9)
            ax.tick_params(labelsize=7)
            ax.grid(True, alpha=0.15)

    fig.suptitle(suptitle, fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(str(output_dir / filename), bbox_inches='tight')
    plt.close(fig)


def save_plots(results_df: pd.DataFrame, summary_df: pd.DataFrame, output_dir: Path):
    """Generate comparison boxplots and heatmaps."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib/seaborn not available -- skipping plots")
        return

    sns.set_style('whitegrid')
    plt.rcParams['figure.dpi'] = 120

    if not results_df.empty:
        models = results_df['recall_model'].unique()
        n_models = len(models)
        fig, axes = plt.subplots(1, n_models, figsize=(5 * n_models, 5),
                                 sharey=True, squeeze=False)
        for ax, model in zip(axes[0], models):
            subset = results_df[results_df['recall_model'] == model]
            sns.boxplot(data=subset, x='sort_strategy', y='target_percentile',
                        ax=ax, palette='Set2', order=sorted(subset['sort_strategy'].unique()))
            ax.set_title(f'Model: {model}')
            ax.set_xlabel('Sort strategy')
            ax.tick_params(axis='x', rotation=30)
        fig.suptitle('Target percentile by strategy and model', fontsize=14)
        plt.tight_layout()
        plt.savefig(str(output_dir / 'comparison_target_percentile.png'),
                    bbox_inches='tight')
        plt.close(fig)

        pivot = summary_df.pivot_table(
            index='sort_strategy', columns='recall_model',
            values='mean_percentile',
        )
        if not pivot.empty:
            fig, ax = plt.subplots(figsize=(8, 5))
            sns.heatmap(pivot, annot=True, fmt='.3f', cmap='YlOrRd',
                        ax=ax, cbar_kws={'label': 'Mean target percentile'})
            ax.set_title('Mean target percentile')
            plt.tight_layout()
            plt.savefig(str(output_dir / 'comparison_heatmap.png'),
                        bbox_inches='tight')
            plt.close(fig)

        fig, axes = plt.subplots(1, n_models, figsize=(5 * n_models, 5),
                                 sharey=True, squeeze=False)
        for ax, model in zip(axes[0], models):
            subset = results_df[results_df['recall_model'] == model]
            sns.boxplot(data=subset, x='sort_strategy', y='actual_recall',
                        ax=ax, palette='Set2', order=sorted(subset['sort_strategy'].unique()))
            ax.axhline(subset['actual_recall'].median(), color='red', ls='--', lw=1, alpha=0.5)
            ax.set_title(f'Model: {model}')
            ax.set_xlabel('Sort strategy')
            ax.tick_params(axis='x', rotation=30)
        fig.suptitle('Actual recall at predicted cutoff', fontsize=14)
        plt.tight_layout()
        plt.savefig(str(output_dir / 'comparison_actual_recall.png'),
                    bbox_inches='tight')
        plt.close(fig)

        if 'standardized_recall' in results_df.columns:
            fig, axes = plt.subplots(1, n_models, figsize=(5 * n_models, 5),
                                     sharey=True, squeeze=False)
            for ax, model in zip(axes[0], models):
                subset = results_df[results_df['recall_model'] == model]
                sns.boxplot(data=subset, x='sort_strategy', y='standardized_recall',
                            ax=ax, palette='Set2',
                            order=sorted(subset['sort_strategy'].unique()))
                ax.axhline(subset['standardized_recall'].median(),
                           color='red', ls='--', lw=1, alpha=0.5)
                ax.set_title(f'Model: {model}')
                ax.set_xlabel('Sort strategy')
                ax.tick_params(axis='x', rotation=30)
            fig.suptitle('Standardized recall (actual / max possible)', fontsize=14)
            plt.tight_layout()
            plt.savefig(str(output_dir / 'comparison_standardized_recall.png'),
                        bbox_inches='tight')
            plt.close(fig)

            _density_grid(
                results_df=results_df,
                x_col='target_percentile',
                y_col='standardized_recall',
                x_label='Target percentile',
                y_label='Standardized recall',
                suptitle='Target percentile vs standardized recall\n(model x strategy density)',
                filename='comparison_percentile_vs_recall.png',
                output_dir=output_dir,
            )

        if 'true_cutoff' in results_df.columns:
            _density_grid(
                results_df=results_df,
                x_col='true_cutoff',
                y_col='target_percentile',
                x_label='True cutoff (last valid taxid)',
                y_label='Predicted cutoff',
                suptitle='True vs predicted cutoff\n(model x strategy density)',
                filename='comparison_true_vs_predicted_cutoff.png',
                output_dir=output_dir,
            )

    logger.info("Plots saved to %s", output_dir)


# -- CLI --

def get_args():
    parser = argparse.ArgumentParser(
        description="Compare recall sort strategies across recall model types.",
    )
    parser.add_argument("--study_output_filepath", type=str, required=True)
    parser.add_argument("--taxid_plan_filepath", type=str, required=True)
    parser.add_argument("--analysis_output_filepath", type=str, required=True)
    parser.add_argument("--output_dir", type=str, default="sort_comparison",
                        help="Subdirectory name within analysis_output_filepath for results "
                             "(default: sort_comparison)")
    parser.add_argument("--sort_strategies", type=str, nargs='+',
                        default=DEFAULT_SORT_STRATEGIES,
                        help="Sort strategies to compare")
    parser.add_argument("--recall_models", type=str, nargs='+',
                        default=DEFAULT_RECALL_MODELS,
                        help="Recall model types to compare")
    parser.add_argument("--target_recall", type=float, default=1.0,
                        help="Target recall threshold")
    parser.add_argument("--data_set_divide", type=int, default=16,
                        help="Number of recall divisions")
    parser.add_argument("--tax_level", type=str, default='genus',
                        help="Taxonomic level")
    parser.add_argument("--max_training", type=int, default=None,
                        help="Max training datasets to use")
    parser.add_argument("--verbose", action="store_true",
                        help="Verbose logging")
    return parser.parse_args()


def main():
    args = get_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        stream=sys.stdout,
    )

    analysis_root = Path(args.analysis_output_filepath)
    output_path = analysis_root / args.output_dir
    output_path.mkdir(parents=True, exist_ok=True)

    config = EvaluatorConfig(
        study_output_filepath=Path(args.study_output_filepath),
        taxid_plan_filepath=Path(args.taxid_plan_filepath),
        analysis_output_filepath=analysis_root,
        target_recall=args.target_recall,
        data_set_divide=args.data_set_divide,
        tax_level=args.tax_level,
        max_training=args.max_training,
        holdout_proportion=0.3,
        recall_sort_strategy='reads',
    )

    logger.info("=" * 60)
    logger.info("Sort-strategy comparison sweep")
    logger.info("=" * 60)
    logger.info("Sort strategies: %s", args.sort_strategies)
    logger.info("Recall models:  %s", args.recall_models)
    logger.info("Target recall:  %.2f", args.target_recall)
    logger.info("Output dir:     %s", output_path)
    logger.info("=" * 60)

    results_df = run_comparison(
        config=config,
        sort_strategies=args.sort_strategies,
        recall_models=args.recall_models,
        target_recall=args.target_recall,
    )

    if results_df.empty:
        logger.error("No results produced. Check data availability.")
        sys.exit(1)

    results_path = output_path / "sort_strategy_comparison_results.tsv"
    results_df.to_csv(results_path, sep="\t", index=False)
    logger.info("Detailed results saved to %s", results_path)

    summary_df = generate_summary(results_df)
    summary_path = output_path / "sort_strategy_comparison_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    logger.info("Summary saved to %s", summary_path)

    print("\n" + "=" * 80)
    print("SUMMARY: Mean target_percentile by (sort_strategy, recall_model)")
    print("=" * 80)
    pivot = summary_df.pivot_table(
        index='sort_strategy', columns='recall_model',
        values='mean_percentile', aggfunc='first',
    )
    print(pivot.to_string(float_format="%.4f"))
    print("\nMean actual recall:")
    pivot_recall = summary_df.pivot_table(
        index='sort_strategy', columns='recall_model',
        values='mean_actual_recall', aggfunc='first',
    )
    print(pivot_recall.to_string(float_format="%.4f"))

    if 'mean_standardized_recall' in summary_df.columns:
        print("\nMean standardized recall (actual / max_possible):")
        pivot_std = summary_df.pivot_table(
            index='sort_strategy', columns='recall_model',
            values='mean_standardized_recall', aggfunc='first',
        )
        print(pivot_std.to_string(float_format="%.4f"))

    print("=" * 80)

    save_plots(results_df, summary_df, output_path)

    logger.info("Done.")


if __name__ == "__main__":
    main()
