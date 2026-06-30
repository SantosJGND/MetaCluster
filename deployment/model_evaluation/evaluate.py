"""
Evaluator module entry point.

This module provides the main entry point and orchestration functions
for the model evaluation pipeline.
"""

import logging
import os

import matplotlib.pyplot as plt
import seaborn as sns

from .metrics import compute_cross_hit_stats, compute_spurious_hit_stats
from .models import TrainingCrossHitAnalyzer


def analyze_cross_hit_distribution(
    results: list,
    output_dir: str,
    top_n: int = 15
):
    """
    Analyze cross-hit distribution from DatasetResult objects.
    
    Saves CSV files and generates violin/scatter plots.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    per_dataset_df, class_stats = compute_cross_hit_stats(results)
    
    if per_dataset_df.empty:
        return None, None
    
    per_dataset_df.to_csv(os.path.join(output_dir, 'cross_hit_metrics_per_dataset.tsv'), sep='\t', index=False)
    class_stats.to_csv(os.path.join(output_dir, 'cross_hit_metrics_stats.tsv'), sep='\t', index=False)
    
    top_classes = class_stats[class_stats['n_simulations'] >= 3].nlargest(top_n, 'n_simulations')['class'].tolist()
    
    if top_classes:
        plot_df = per_dataset_df[per_dataset_df['class'].isin(top_classes)]
        
        for metric, ylabel in [
            ('reads_simulated', 'Reads Simulated'),
            ('cross_hit_count', 'Cross-Hit Count'),
            ('cross_hit_reads_mapped', 'Cross-Hit Reads Mapped'),
            ('ratio', 'Cross-Hit Reads / Reads Simulated'),
        ]:
            fig, ax = plt.subplots(figsize=(12, 6))
            sns.boxplot(data=plot_df, x='class', y=metric, order=top_classes,
                          ax=ax, color='steelblue')
            ax.set_xlabel('Class', fontsize=12)
            ax.set_ylabel(ylabel, fontsize=12)
            ax.tick_params(axis='x', rotation=45)
            ax.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'cross_hit_{metric}_distribution.png'), dpi=150)
            plt.close()
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.scatterplot(data=plot_df, x='reads_simulated', y='cross_hit_count',
                       ax=ax, alpha=0.6)
        ax.set_xlabel('Reads Simulated', fontsize=12)
        ax.set_ylabel('Cross-Hit Count', fontsize=12)
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'cross_hit_dotplot_reads_vs_count.png'), dpi=150)
        plt.close()
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.scatterplot(data=plot_df, x='reads_simulated', y='cross_hit_reads_mapped',
                       ax=ax, alpha=0.6)
        ax.set_xlabel('Reads Simulated', fontsize=12)
        ax.set_ylabel('Cross-Hit Reads Mapped', fontsize=12)
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'cross_hit_dotplot_reads_vs_mapped.png'), dpi=150)
        plt.close()
    
    return per_dataset_df, class_stats


def analyze_spurious_hit_distribution(
    results: list,
    output_dir: str,
    top_n: int = 15
):
    """
    Analyze spurious-hit distribution from DatasetResult objects.
    
    Saves CSV files and generates violin/scatter plots.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    per_dataset_df, class_stats = compute_spurious_hit_stats(results)
    
    if per_dataset_df.empty:
        return None, None
    
    per_dataset_df.to_csv(os.path.join(output_dir, 'spurious_hit_metrics_per_dataset.tsv'), sep='\t', index=False)
    class_stats.to_csv(os.path.join(output_dir, 'spurious_hit_metrics_stats.tsv'), sep='\t', index=False)
    
    top_classes = class_stats[class_stats['n_simulations'] >= 3].nlargest(top_n, 'n_simulations')['class'].tolist()
    
    if top_classes:
        plot_df = per_dataset_df[per_dataset_df['class'].isin(top_classes)]
        
        for metric, ylabel in [
            ('spurious_hit_count', 'Spurious-Hit Count'),
            ('spurious_hit_reads_mapped', 'Spurious-Hit Reads Mapped'),
            ('ratio', 'Spurious-Hit Reads / Reads Simulated'),
        ]:
            fig, ax = plt.subplots(figsize=(12, 6))
            sns.boxplot(data=plot_df, x='class', y=metric, order=top_classes,
                          ax=ax, color='salmon')
            ax.set_xlabel('Class', fontsize=12)
            ax.set_ylabel(ylabel, fontsize=12)
            ax.tick_params(axis='x', rotation=45)
            ax.grid(axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'spurious_hit_{metric}_distribution.png'), dpi=150)
            plt.close()
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.scatterplot(data=plot_df, x='reads_simulated', y='spurious_hit_count',
                       ax=ax, alpha=0.6)
        ax.set_xlabel('Reads Simulated', fontsize=12)
        ax.set_ylabel('Spurious-Hit Count', fontsize=12)
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'spurious_hit_dotplot_reads_vs_count.png'), dpi=150)
        plt.close()
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.scatterplot(data=plot_df, x='reads_simulated', y='spurious_hit_reads_mapped',
                       ax=ax, alpha=0.6)
        ax.set_xlabel('Reads Simulated', fontsize=12)
        ax.set_ylabel('Spurious-Hit Reads Mapped', fontsize=12)
        ax.set_xscale('log')
        ax.set_yscale('log')    
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'spurious_hit_dotplot_reads_vs_mapped.png'), dpi=150)
        plt.close() 

def get_args():
    import argparse

    parser = argparse.ArgumentParser(description="Model Deployment for Overlap Manager")
    parser.add_argument("--study_output_filepath", type=str, required=True, help="Path to the study output directory")
    parser.add_argument("--taxid_plan_filepath", type=str, required=True, help="Path to the taxid plan file")
    parser.add_argument("--analysis_output_filepath", type=str, required=True, help="Path to save analysis outputs")
    parser.add_argument("--threshold", type=float, default=0.3, help="Threshold value for model")
    parser.add_argument("--taxa_threshold", type=float, default=0.02, help="Taxa threshold for filtering")
    parser.add_argument("--tax_level_to_use", type=str, default='order', help="Taxonomic level to use")
    parser.add_argument("--data_set_divide", type=int, default=16, help="Data set divide for training/testing")
    parser.add_argument("--target_recall", type=float, default=1.0, help="Target recall threshold for filtering (default: 1.0)")
    parser.add_argument("--holdout_proportion", type=float, default=0.3, help="Proportion of data to hold out for testing")
    parser.add_argument("--output_db_dir", type=str, required=False, help="Path to the output database directory")
    parser.add_argument("--max_training", type=str, default=None, help="Maximum number of training datasets to use")
    parser.add_argument("--max_taxids_fixed_filter", type=int, default=15, help="Maximum number of taxids for fixed filter")
    parser.add_argument("--recall_model_interface", type=str, default='xgb', help="Model interface: 'xgb' (XGBoost multi-output), 'morf' (RF multi-output), 'moxgb_optimized', 'morf_optimized', 'monn_optimized', 'direct' (RF cutoff classifier), 'direct_xgb' (DirectXGB regressor), 'gp_clf' (GP+CLF)")
    parser.add_argument("--composition_model_interface", type=str, default='xgb', help="Composition model: 'xgb' (default, XGBoost), 'xgb_optimized' (XGBoost+Optuna), 'rf' (RandomForest), 'gb' (GradientBoosting), 'lr' (LogisticRegression, stats-only)")
    parser.add_argument("--cutoff_confidence", type=float, default=None, help="Confidence for probability-guided cutoff (0-1, default: point prediction)")
    parser.add_argument("--recall_sort_strategy", type=str, default='reads', help="Sort strategy for recall: 'reads', 'taxid_roundrobin', 'rarity_boost', 'tax_level_stratified'")
    parser.add_argument("--cross_hit_threshold", type=float, default=0.95, help="Cross-hit probability threshold for filtering")
    parser.add_argument("--enable-cross-hit", action=argparse.BooleanOptionalAction, default=False, 
                        help="Enable cross-hit cleanup during evaluation (default: enabled)")
    
    parser.add_argument("--log-format", type=str, default="text", choices=["json", "text"], help="Log output format")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose console logging")
    parser.add_argument("--generate-report", action="store_true", default=True, help="Generate HTML report (default: enabled)")
    parser.add_argument("--no-report", action="store_true", help="Disable HTML report generation")
    parser.add_argument("--use-mlflow", action=argparse.BooleanOptionalAction, default=True, help="Enable MLflow tracking (default: enabled)")
    parser.add_argument("--mlflow-uri", type=str, default=None, help="MLflow tracking URI")
    parser.add_argument("--no-cache", action="store_true", help="Force recompute cached training data")
    
    return parser.parse_args()


def main(args):
    """
    Main function using the refactored architecture.
    
    Uses:
    - Pydantic v2 config
    - Structured logging
    - DataLoader class
    - ModelTrainer class
    - BatchEvaluator class
    - ResultVisualizer class
    - MLflow tracking (optional)
    - HTML report generation (optional)
    """
    from datetime import datetime
    
    from deployment.model_evaluation.config import EvaluatorConfig
    from deployment.model_evaluation.logging_config import setup_logging
    from deployment.model_evaluation.data_loader import DataLoader
    from deployment.model_evaluation.models import ModelTrainer
    from deployment.model_evaluation.batch_evaluator import BatchEvaluator
    from deployment.model_evaluation.visualization import ResultVisualizer, generate_report
    from deployment.model_evaluation.models import MLflowTracker
    
    config = EvaluatorConfig.from_args(args)
    
    if args.no_report:
        config.generate_report = False
    
    logger = setup_logging(
        log_dir=config.analysis_output_filepath,
        log_format=args.log_format,
        console_level=logging.DEBUG if args.verbose else logging.INFO
    )
    
    logger.info("=" * 60)
    logger.info("Starting Model Deployment Analysis")
    logger.info("=" * 60)
    logger.info(f"Study output: {config.study_output_filepath}")
    logger.info(f"Analysis output: {config.analysis_output_filepath}")
    logger.info(f"Tax level: {config.tax_level}")
    logger.info(f"MLflow enabled: {config.use_mlflow}")
    logger.info(f"Report enabled: {config.generate_report}")
    logger.info("=" * 60)
    
    mlflow_tracker = None
    if config.use_mlflow:
        try:
            mlflow_tracker = MLflowTracker(
                experiment_name="metagenomics-evaluation",
                tracking_uri=config.mlflow_uri
            )
            mlflow_tracker.start_run(run_name=f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
            mlflow_tracker.log_params({
                'tax_level': config.tax_level,
                'data_set_divide': config.data_set_divide,
                'cross_hit_threshold': config.cross_hit_threshold,
                'taxa_threshold': config.taxa_threshold,
                'holdout_proportion': config.holdout_proportion,
            })
            logger.info(f"MLflow tracking started: {config.mlflow_uri or 'default'}")
        except ImportError:
            logger.warning("MLflow not installed, disabling tracking")
            config.use_mlflow = False
    
    logger.info("Loading data...")
    loader = DataLoader(config).initialize()
    
    training_folders = loader.get_training_folders()
    test_folders = loader.get_test_folders()
    
    logger.info(f"Training datasets: {len(training_folders)}")
    logger.info(f"Test datasets: {len(test_folders)}")
    
    if mlflow_tracker:
        mlflow_tracker.log_metrics({
            'n_training_datasets': len(training_folders),
            'n_test_datasets': len(test_folders),
        })
    
    logger.info("Training models...")
    trainer = ModelTrainer(
        config=config,
        ncbi_wrapper=loader.get_ncbi_wrapper(),
        input_tax_df=loader.get_input_tax_df(),
        taxids_to_use=loader.get_taxids_to_use(),
        use_cache=config.use_cache,
    )
    
    trainer.train_models(training_folders)
    trainer.save_models(str(config.models_dir))
    trainer.evaluate_models(str(config.models_dir))
    
    logger.info("Analyzing cross-hit distribution from training data...")
    training_analyzer = TrainingCrossHitAnalyzer(
        config=config,
        ncbi_wrapper=loader.get_ncbi_wrapper(),
        input_tax_df=loader.get_input_tax_df(),
        taxids_to_use=loader.get_taxids_to_use(),
    )
    training_results = training_analyzer.analyze_training_data(training_folders)
    if training_results:
        analyze_cross_hit_distribution(training_results, str(config.analysis_output_filepath), top_n=15)
        analyze_spurious_hit_distribution(training_results, str(config.analysis_output_filepath), top_n=15)
        
    logger.info("Evaluating on test datasets...")
    evaluator = BatchEvaluator(
        config=config,
        models=trainer.models,
        ncbi_wrapper=loader.get_ncbi_wrapper(),
        input_tax_df=loader.get_input_tax_df(),
        taxids_to_use=loader.get_taxids_to_use(),
    )
    
    results = evaluator.evaluate(test_folders)
    
    logger.info("Saving results...")
    if config.use_cache is False:
        trainer._stats_matrices.to_csv(str(config.analysis_output_filepath / "training_stats_matrices.tsv"), sep="\t", index=False)
    results.save_tsv(str(config.analysis_output_filepath))
    results.to_json(str(config.analysis_output_filepath / "evaluation_results.json"))
    
    logger.info("Generating visualizations...")
    visualizer = ResultVisualizer(str(config.analysis_output_filepath))
    visualizer.plot_all(results)
    
    if config.generate_report:
        logger.info("Generating HTML report...")
        report_path = generate_report(results, str(config.analysis_output_filepath))
        logger.info(f"Report generated: {report_path}")
    
    if mlflow_tracker:
        final_metrics = {
            'test_precision_mean': results.test_results['precision_clade_post'].mean() if not results.test_results.empty else 0,
            'test_precision_std': results.test_results['precision_clade_post'].std() if not results.test_results.empty else 0,
            'n_datasets_evaluated': results.get_dataset_count(),
        }
        mlflow_tracker.log_metrics(final_metrics)
        mlflow_tracker.end_run()
    
    logger.info("=" * 60)
    logger.info("Pipeline completed successfully!")
    logger.info("=" * 60)
    
    return results


if __name__ == "__main__":
    args = get_args()
    main(args)
