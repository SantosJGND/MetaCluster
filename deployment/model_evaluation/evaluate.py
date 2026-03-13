import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import logging
import os

from metagenomics_utils.overlap_manager import OverlapManager
from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper
from metagenomics_utils.overlap_manager.node_stats import get_m_stats_matrix 
from metagenomics_utils.overlap_manager.om_models import data_set_traversal_with_precision, cross_hit_prediction_matrix, predict_recall_cutoff_vars
from metagenomics_utils.overlap_manager.om_models import get_trash_composition, get_cross_hit_composition, calculate_overall_precision, predict_data_set_clades, cross_hit_prediction, cut_off_recall_prediction
from metagenomics_utils.overlap_manager.om_models import CompositionModeller, RecallModeller, CrossHitModeller


def retrieve_simulation_input(study_output_filepath: str) -> pd.DataFrame:
    """Retrieve simulation input data from study output directories.
    
    Args:
        study_output_filepath: Path to the study output directory containing dataset subfolders.
    
    Returns:
        DataFrame with all input data merged from all datasets.
    """
    input_files = []
    folders = [f for f in os.listdir(study_output_filepath) if os.path.isdir(os.path.join(study_output_filepath, f))]

    for data_set_name in folders:

        output_filepath = os.path.join(study_output_filepath, f"{data_set_name}", "input", f"{data_set_name}.tsv")
        df = pd.read_csv(output_filepath, sep="\t")

        df['data_set'] = data_set_name.replace("_plan", "")
        if df.empty:
            continue
        input_files.append(df)

    if not input_files:
        return pd.DataFrame()
    
    all_input_data = pd.concat(input_files, ignore_index=True)

    return all_input_data



def process_clade_report(data_set_name: str, study_output_filepath: str) -> list:
    """Process clade report for a dataset and extract unique taxids.
    
    Args:
        data_set_name: Name of the dataset.
        study_output_filepath: Path to the study output directory.
    
    Returns:
        List of unique taxids found in the clade report.
    """
    clade_report = os.path.join(study_output_filepath, f"{data_set_name}", "output", "clade_report_with_references.tsv")
    if not os.path.exists(clade_report):
        return []

    clade_report = pd.read_csv(clade_report, sep="\t").rename(columns={"#rname": "assembly_accession"})

    taxids = clade_report['taxid'].dropna().unique().tolist() 
    return taxids


def output_parse(study_output_filepath: str, ncbi_wrapper: NCBITaxonomistWrapper, input_taxids: list) -> pd.DataFrame:
    """Parse output from study and create taxid table with lineage information.
    
    Args:
        study_output_filepath: Path to the study output directory.
        ncbi_wrapper: NCBI TaxonomistWrapper for lineage resolution.
        input_taxids: List of input taxids to include.
    
    Returns:
        DataFrame with taxid lineage information (order, family, genus).
    """
    folders = [f for f in os.listdir(study_output_filepath) if os.path.isdir(os.path.join(study_output_filepath, f))]
    all_output_taxids = set()
    for data_set_name in folders:
        output_taxids = process_clade_report(data_set_name, study_output_filepath)
        all_output_taxids.update(output_taxids)
    
    taxids = list(all_output_taxids) + input_taxids

    ncbi_wrapper.resolve_lineages(taxids)

    output_taxids_table = pd.DataFrame(list(all_output_taxids), columns=['taxid'])
    output_taxids_table['order'] = output_taxids_table.apply(lambda row: ncbi_wrapper.get_level(row['taxid'], 'order'), axis=1)
    output_taxids_table['family'] = output_taxids_table.apply(lambda row: ncbi_wrapper.get_level(row['taxid'], 'family'), axis=1)
    output_taxids_table['genus'] = output_taxids_table.apply(lambda row: ncbi_wrapper.get_level(row['taxid'], 'genus'), axis=1)

    return output_taxids_table

def establish_taxids_to_use(study_output_filepath: str, ncbi_wrapper: NCBITaxonomistWrapper, input_taxids: list, tax_level_to_use: str = 'order', min_tax_count: float = 0.02, normalize: bool = True) -> pd.DataFrame:
    """Establish taxids to use based on frequency threshold at a specific taxonomic level.
    
    Args:
        study_output_filepath: Path to the study output directory.
        ncbi_wrapper: NCBI TaxonomistWrapper for lineage resolution.
        input_taxids: List of input taxids.
        tax_level_to_use: Taxonomic level to use for filtering (default: 'order').
        min_tax_count: Minimum frequency threshold (default: 0.02).
        normalize: Whether to normalize counts (default: True).
    
    Returns:
        DataFrame of taxids to use with lineage information.
    """
    output_taxids_table = output_parse(study_output_filepath, ncbi_wrapper, input_taxids)
    taxids_to_use = output_taxids_table[output_taxids_table[tax_level_to_use].map(output_taxids_table[tax_level_to_use].value_counts(normalize=normalize)) > min_tax_count] 

    taxids_to_use = taxids_to_use.dropna(subset=[tax_level_to_use]).drop_duplicates(subset=[tax_level_to_use]).reset_index(drop=True)
    taxids_to_use = pd.concat([taxids_to_use, pd.DataFrame({'taxid': [0], 'order': ['unclassified'], 'family': ['unclassified'], 'genus': ['unclassified']})], ignore_index=True)
    return taxids_to_use
        

def expand_input_data(all_input_data: pd.DataFrame, ncbi_wrapper: NCBITaxonomistWrapper, taxid_plan:pd.DataFrame) -> pd.DataFrame:
    """
    Expand input data with lineage information from NCBI Taxonomist
    """
    all_input_data = all_input_data.merge(taxid_plan, on='taxid', how='left')
    input_tax_df = pd.DataFrame(all_input_data['taxid'].dropna().unique(), columns=['taxid'])
    input_tax_df['order'] = input_tax_df.apply(lambda row: ncbi_wrapper.get_level(row['taxid'], 'order'), axis=1)
    input_tax_df['family'] = input_tax_df.apply(lambda row: ncbi_wrapper.get_level(row['taxid'], 'family'), axis=1)
    input_tax_df = input_tax_df.drop_duplicates(subset=['taxid'])

    input_tax_df = pd.concat([input_tax_df, pd.DataFrame({'taxid': [0], 'order': ['unclassified'], 'family': ['unclassified'], 'genus': ['unclassified']})], ignore_index=True)
    input_tax_df = input_tax_df.replace({np.nan: 'unclassified'})
    
    return input_tax_df

def run_data_retrieval(trainning_folders:list, study_output_filepath: str, data_set_divide: int, ncbi_wrapper: NCBITaxonomistWrapper, input_tax_df: pd.DataFrame, tax_level_to_use: str, taxids_to_use: pd.DataFrame):
    
    trainning_results = []
    prediction_trainning_results = []
    recall_trainning_results = []
        
    for data_set_name in trainning_folders:
        try:

            overlap_manager = OverlapManager(os.path.join(study_output_filepath, f"{data_set_name}", "clustering"))

            if overlap_manager.m_stats_matrix.empty:
                continue

            result_df = data_set_traversal_with_precision(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager, taxids_to_use, tax_level=tax_level_to_use)
            prediction_matrix = cross_hit_prediction_matrix(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager, taxids_to_use, tax_level=tax_level_to_use)
            m_stats_stats_matrix = get_m_stats_matrix(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager)
            if m_stats_stats_matrix.empty:
                continue
            recall_stats = predict_recall_cutoff_vars(data_set_divide, data_set_name, m_stats_stats_matrix, taxids_to_use, tax_level= tax_level_to_use)
            if not result_df.empty:
                trainning_results.append(result_df)
            if not prediction_matrix.empty:
                prediction_trainning_results.append(prediction_matrix)
            if not recall_stats.empty:
                recall_trainning_results.append(recall_stats)
        except Exception as e:
            print(f"Error: {e}")
            print(f"Processing {data_set_name}")
            import traceback
            traceback.print_exc()


    trainning_results_df = pd.concat(trainning_results, ignore_index=True)
    prediction_trainning_results_df = pd.concat(prediction_trainning_results, ignore_index=True)
    recall_trainning_results = pd.concat(recall_trainning_results, ignore_index=True)

    return trainning_results_df, prediction_trainning_results_df, recall_trainning_results



def compound_eda_function(remaining_folders,
                            study_output_filepath: str, 
                          data_set_divide: int,
                            ncbi_wrapper: NCBITaxonomistWrapper, 
                            input_tax_df: pd.DataFrame, 
                            tax_level_to_use: str, 
                            taxids_to_use: pd.DataFrame, 
                            recall_predict_modeller: RecallModeller, 
                            composition_modeller: CompositionModeller,
                            cross_hit_modeller: CrossHitModeller,
                            cross_hit_threshold: float = 0.9):

    test_results = []
    summary_results = []
    cross_hit_results = []
    trash_results = []

    for data_set_name in remaining_folders:
        try:
            overlap_manager = OverlapManager(os.path.join(study_output_filepath, f"{data_set_name}", "clustering"))
            input_df_filepath = os.path.join(study_output_filepath, f"{data_set_name}", "input", f"{data_set_name}.tsv")
            if not os.path.exists(input_df_filepath):
                print(f"Input file {input_df_filepath} does not exist. Skipping {data_set_name}.")
                continue
            if overlap_manager.m_stats_matrix.empty:
                #print(f"No reads mapped for {data_set_name}. Skipping.")
                continue
            input_df = pd.read_csv(input_df_filepath, sep="\t")
            input_df_summary = input_df[['sample', 'taxid', 'reads', 'mutation_rate']].drop_duplicates()
            input_df_summary = input_df_summary.merge(input_tax_df, on='taxid', how='left')
            input_df_summary.loc[:,'output_raw'] = overlap_manager.m_stats_matrix.shape[0]
            input_df_summary.loc[:,'output_cov_filtered'] = overlap_manager.m_stats_matrix[overlap_manager.m_stats_matrix['coverage'] > 0].shape[0]
            # pre-filter leaves with zero coverage
            m_stats_stats_matrix = get_m_stats_matrix(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager)
            m_stats_stats_matrix['is_trash'] = (m_stats_stats_matrix['best_match_is_best'] == False) & (m_stats_stats_matrix['is_crosshit'] == False)
            trash_composition = get_trash_composition(input_df_summary, m_stats_stats_matrix, input_tax_df, tax_level=tax_level_to_use)
            cross_hit_composition = get_cross_hit_composition(input_df_summary, m_stats_stats_matrix, input_tax_df, tax_level=tax_level_to_use)
            
            clean_m_stats = m_stats_stats_matrix[m_stats_stats_matrix['is_trash'] == False]
            input_df_summary.loc[:, 'raw_pred_accuracy'] = input_df_summary['taxid'].apply(lambda x: (m_stats_stats_matrix['best_match_taxid'] == x).sum())
            input_df_summary.loc[:, 'fuzzy_precision_raw'] = clean_m_stats.shape[0] / m_stats_stats_matrix.shape[0] if m_stats_stats_matrix.shape[0] > 0 else 0.0
            input_df_summary.loc[:, 'fuzzy_precision_cov_filtered'] = clean_m_stats[(clean_m_stats['coverage'] > 0)].shape[0] / m_stats_stats_matrix.shape[0] if m_stats_stats_matrix.shape[0] > 0 else 0.0
            input_df_summary.loc[:, 'overall_precision_raw'] = clean_m_stats.dropna(subset=['best_match_taxid'])["best_match_taxid"].nunique() / m_stats_stats_matrix.shape[0] if m_stats_stats_matrix.shape[0] > 0 else 0.0
            input_df_summary.loc[:,'recall_raw'] = clean_m_stats.drop_duplicates(subset=['best_match_taxid']).dropna(subset=['best_match_taxid']).shape[0] / input_df_summary['taxid'].nunique() if input_df_summary['taxid'].nunique() > 0 else 0.0
            input_df_summary.loc[:,'recall_cov_filtered'] = clean_m_stats[(clean_m_stats['coverage'] > 0)].drop_duplicates(subset=['best_match_taxid']).dropna(subset=['best_match_taxid']).shape[0] / input_df_summary['taxid'].nunique() if input_df_summary['taxid'].nunique() > 0 else 0.0

            ### recall prediction to filter leaves
            overlap_manager = cut_off_recall_prediction(study_output_filepath, data_set_name, recall_predict_modeller, data_set_divide, 
                                                        m_stats_stats_matrix, taxids_to_use)
            m_stats_stats_matrix = get_m_stats_matrix(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager)
            m_stats_stats_matrix['is_trash'] = (m_stats_stats_matrix['best_match_is_best'] == False) & (m_stats_stats_matrix['is_crosshit'] == False)
            clean_m_stats = m_stats_stats_matrix[m_stats_stats_matrix['is_trash'] == False]
            input_df_summary.loc[:,'recall_filtered_leaves'] = clean_m_stats.drop_duplicates(subset=['best_match_taxid']).dropna(subset=['best_match_taxid']).shape[0] / input_df_summary['taxid'].nunique() if input_df_summary['taxid'].nunique() > 0 else 0.0

            missing_tax_levels = set(input_tax_df[tax_level_to_use].unique()) - set(input_df_summary[tax_level_to_use].unique())
            trash_composition.loc[:, list(missing_tax_levels)] = np.nan
            cross_hit_composition.loc[:, list(missing_tax_levels)] = np.nan


            results_df_pre = predict_data_set_clades(data_set_name, m_stats_stats_matrix, overlap_manager, composition_modeller, taxids_to_use, tax_level=tax_level_to_use)
            precision_pre = calculate_overall_precision(results_df_pre, m_stats_stats_matrix)
            
            input_df_summary.loc[:,'predicted_clades'] = results_df_pre.shape[0]
            input_df_summary.loc[:,'clade_pred_accuracy_pre'] = input_df_summary['taxid'].apply(lambda x: (results_df_pre['best_taxid_match'] == x).sum())
            input_df_summary.loc[:,'clade_precision_full'] = precision_pre
            input_df_summary.loc[:,'clade_recall'] = results_df_pre['best_taxid_match'].dropna().nunique() / input_df_summary['taxid'].nunique() if input_df_summary['taxid'].nunique() > 0 else 0.0

            cross_hit_predictions = cross_hit_prediction(data_set_name, 
                                                         study_output_filepath, 
                                                         ncbi_wrapper, 
                                                         cross_hit_modeller,
                                                         overlap_manager, 
                                                         taxids_to_use, tax_level=tax_level_to_use
                                                        )
            cross_hit_predictions = cross_hit_predictions[cross_hit_predictions['prob_best_match'] > cross_hit_threshold]
            input_df_summary.loc[:,'predicted_cross_hits'] = cross_hit_predictions.shape[0]
            input_df_summary.loc[:, 'cleanup_accuracy'] = sum(cross_hit_predictions['is_trash']) / cross_hit_predictions.shape[0] if cross_hit_predictions.shape[0] > 0 else 0.0

            if not cross_hit_predictions.empty:
                if len(overlap_manager.leaves) == 0:
                    continue

                cross_hits = cross_hit_predictions[cross_hit_predictions['prob_best_match'] > cross_hit_threshold]['leaf'].tolist()
                cross_hits = [leaf for leaf in cross_hits if leaf in overlap_manager.leaves]
                overlap_manager.m_stats_matrix.loc[overlap_manager.m_stats_matrix.index.isin(cross_hits), 'numreads'] = 0
                overlap_manager.m_stats_matrix.loc[overlap_manager.m_stats_matrix.index.isin(cross_hits), 'coverage'] = 0

                if overlap_manager.m_stats_matrix.shape[0] > 0:
                    overlap_manager.prune_empty_nodes()
                    overlap_manager.new_tree_from_distance_matrix()
                    overlap_manager.recalculate_all_min_pairwise_dist()

            prediction_matrix = get_m_stats_matrix(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager)
            if prediction_matrix.shape[0] < 2:
                continue
            result_df = predict_data_set_clades(data_set_name, prediction_matrix, overlap_manager, composition_modeller, taxids_to_use, tax_level=tax_level_to_use)
            if result_df.empty:
                input_df_summary.loc[:,'clade_pred_accuracy_post'] = 0
                input_df_summary.loc[:,'predicted_clades_post'] = 0
                input_df_summary.loc[:,'clade_precision_post'] = 0.0
                summary_results.append(input_df_summary)
                test_results.append((0.0, data_set_name))
                continue
            precision = calculate_overall_precision(result_df, prediction_matrix)
            input_df_summary.loc[:,'clade_pred_accuracy_post'] = input_df_summary['taxid'].apply(lambda x: (result_df['best_taxid_match'] == x).sum())
            input_df_summary.loc[:,'predicted_clades_post'] = result_df.shape[0]
            input_df_summary.loc[:,'clade_precision_post'] = precision
            summary_results.append(input_df_summary)
            test_results.append((precision, data_set_name))
            trash_results.append(trash_composition)
            cross_hit_results.append(cross_hit_composition)
        except Exception as e:
            print(f"Error: {e}")
            print(f"Processing {data_set_name}")
            import traceback
            traceback.print_exc()

    if not summary_results:
        summary_results_df = pd.DataFrame()
    else:
        summary_results_df = pd.concat(summary_results, ignore_index=True)
        summary_results_df.to_csv(os.path.join(study_output_filepath, "test_datasets_summary_results.tsv"), sep="\t", index=False)
    
    test_results_df = pd.DataFrame(test_results, columns=['overall_precision', 'data_set'])
    
    if not trash_results:
        trash_results_df = pd.DataFrame()
    else:
        trash_results_df = pd.concat(trash_results, ignore_index=True)
    
    if not cross_hit_results:
        cross_hit_results_df = pd.DataFrame()
    else:
        cross_hit_results_df = pd.concat(cross_hit_results, ignore_index=True)

    return test_results_df, summary_results_df, trash_results_df, cross_hit_results_df


def compound_eda_function_new(
    remaining_folders,
    study_output_filepath: str, 
    data_set_divide: int,
    ncbi_wrapper: NCBITaxonomistWrapper, 
    input_tax_df: pd.DataFrame, 
    tax_level_to_use: str, 
    taxids_to_use: pd.DataFrame, 
    recall_predict_modeller, 
    composition_modeller,
    cross_hit_modeller,
    cross_hit_threshold: float = 0.9):
    """
    Backward-compatible wrapper around new BatchEvaluator.
    
    This function provides the same interface as compound_eda_function
    but uses the new modular architecture.
    
    Args:
        remaining_folders: List of test dataset folders
        study_output_filepath: Path to study output directory
        data_set_divide: Dataset divide for training/testing
        ncbi_wrapper: NCBI TaxonomistWrapper instance
        input_tax_df: Input tax DataFrame
        tax_level_to_use: Taxonomic level to use
        taxids_to_use: Taxids to use for evaluation
        recall_predict_modeller: Trained recall modeller
        composition_modeller: Trained composition modeller
        cross_hit_modeller: Trained cross-hit modeller
        cross_hit_threshold: Cross-hit probability threshold
        
    Returns:
        Tuple of (test_results_df, summary_results_df, trash_results_df, cross_hit_results_df)
    """
    from .config import EvaluatorConfig, TrainedModels
    from .batch_evaluator import BatchEvaluator
    
    config = EvaluatorConfig(
        study_output_filepath=study_output_filepath,
        data_set_divide=data_set_divide,
        tax_level=tax_level_to_use,
        cross_hit_threshold=cross_hit_threshold,
    )
    
    models = TrainedModels(
        recall_modeller=recall_predict_modeller,
        composition_modeller=composition_modeller,
        crosshit_modeller=cross_hit_modeller,
    )
    
    evaluator = BatchEvaluator(
        config=config,
        models=models,
        ncbi_wrapper=ncbi_wrapper,
        input_tax_df=input_tax_df,
        taxids_to_use=taxids_to_use,
    )
    
    result = evaluator.evaluate(remaining_folders, progress=False)
    
    return (
        result.test_results,
        result.summary_results,
        result.trash_composition,
        result.cross_hit_composition,
    )


def get_args():
    import argparse

    parser = argparse.ArgumentParser(description="Model Deployment for Overlap Manager")
    parser.add_argument("--study_output_filepath", type=str, required=True, help="Path to the study output directory")
    parser.add_argument("--taxid_plan_filepath", type=str, required=True, help="Path to the taxid plan file")
    parser.add_argument("--analysis_output_filepath", type=str, required=True, help="Path to save analysis outputs")
    parser.add_argument("--threshold", type=float, default=0.3, help="Threshold value for model")
    parser.add_argument("--taxa_threshold", type=float, default=0.02, help="Taxa threshold for filtering")
    parser.add_argument("--tax_level_to_use", type=str, default='order', help="Taxonomic level to use")
    parser.add_argument("--data_set_divide", type=int, default=5, help="Data set divide for training/testing")
    parser.add_argument("--holdout_proportion", type=float, default=0.3, help="Proportion of data to hold out for testing")
    parser.add_argument("--output_db_dir", type=str, required=False, help="Path to the output database directory")
    parser.add_argument("--max_training", type=str, default=None, help="Maximum number of training datasets to use")
    parser.add_argument("--cross_hit_threshold", type=float, default=0.9, help="Cross-hit probability threshold for filtering")
    
    # New arguments for enhanced features
    parser.add_argument("--log-format", type=str, default="text", choices=["json", "text"], help="Log output format")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose console logging")
    parser.add_argument("--generate-report", action="store_true", default=True, help="Generate HTML report (default: enabled)")
    parser.add_argument("--no-report", action="store_true", help="Disable HTML report generation")
    parser.add_argument("--use-mlflow", action="store_true", default=True, help="Enable MLflow tracking (default: enabled)")
    parser.add_argument("--mlflow-uri", type=str, default=None, help="MLflow tracking URI")
    parser.add_argument("--use-legacy", action="store_true", help="Use legacy implementation (deprecated)")
    
    return parser.parse_args()


def main_new(args):
    """
    New main function using the refactored architecture.
    
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
    import logging
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
    logger.info("Starting Model Deployment Analysis (NEW ARCHITECTURE)")
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
    )
    
    trainer.train_models(training_folders)
    trainer.save_models(str(config.models_dir))
    trainer.evaluate_models(str(config.models_dir))
    
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
            'test_precision_mean': results.test_results['overall_precision'].mean() if not results.test_results.empty else 0,
            'test_precision_std': results.test_results['overall_precision'].std() if not results.test_results.empty else 0,
            'n_datasets_evaluated': results.get_dataset_count(),
        }
        mlflow_tracker.log_metrics(final_metrics)
        mlflow_tracker.end_run()
    
    logger.info("=" * 60)
    logger.info("Pipeline completed successfully!")
    logger.info("=" * 60)
    
    return results


if __name__ == "__main__":
    import warnings
    import sys
    
    args = get_args()
    
    if args.use_legacy:
        warnings.warn(
            "Using legacy implementation. This is deprecated and will be removed "
            "after tests for the new implementation are written.",
            DeprecationWarning,
            stacklevel=2
        )
    else:
        # Use new implementation
        main_new(args)
        sys.exit(0)
    
    # Legacy code below
    study_output_filepath = args.study_output_filepath
    taxid_plan_filepath = args.taxid_plan_filepath
    analysis_output_filepath = args.analysis_output_filepath
    threshold = args.threshold
    taxa_threshold = args.taxa_threshold
    tax_level_to_use = args.tax_level_to_use
    data_set_divide = args.data_set_divide
    holdout_proportion = args.holdout_proportion
    proportion_train = 1 - holdout_proportion
    max_training = args.max_training
    cross_hit_threshold = args.cross_hit_threshold

    # Example usage:
    # study_output_filepath = "/path/to/study_output"
    # taxid_plan_filepath = "/path/to/taxid_plan.tsv"
    # analysis_output_filepath = "/path/to/output"
    #threshold = 0.3
    #taxa_threshold = 0.02
    #tax_level_to_use = 'order'
    #data_set_divide = 5

    #### output logger to file
    models_subdir = os.path.join(analysis_output_filepath, "models")
    os.makedirs(analysis_output_filepath, exist_ok=True)
    os.makedirs(models_subdir, exist_ok=True)

    logging.basicConfig(level=logging.INFO, filename=os.path.join(analysis_output_filepath, "evaluate.log"), filemode='w',
                        format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    logger.info("Starting model deployment analysis")
    logger.info(f"Study output filepath: {study_output_filepath}")
    logger.info(f"Taxid plan filepath: {taxid_plan_filepath}")
    logger.info(f"Analysis output filepath: {analysis_output_filepath}")
    #### Define output files
    output_lineages = os.path.join(study_output_filepath, "lineages.tsv")
    output_db = os.path.join(study_output_filepath, "taxa.db")

    #### Load taxid plan
    taxid_plan = pd.read_csv(taxid_plan_filepath, sep="\t")
    taxid_plan = taxid_plan[['taxid', 'description','lineage']].drop_duplicates(subset=['taxid'])
    folders = [f for f in os.listdir(study_output_filepath) if os.path.isdir(os.path.join(study_output_filepath, f))]
    from random import sample
    if max_training is not None:
        if len(folders) > int(max_training):
            folders = sample(folders, int(max_training))

    all_input_data = retrieve_simulation_input(study_output_filepath)
    ncbi_wrapper = NCBITaxonomistWrapper(db=output_db)
    input_taxids = all_input_data['taxid'].dropna().unique().tolist()
    ncbi_wrapper.resolve_lineages(input_taxids)

    #### load output taxids
    input_tax_df = expand_input_data(all_input_data, ncbi_wrapper, taxid_plan)
    logger.info(f"Number of unique input taxids: {input_tax_df.shape[0]}")
    logger.info(f"ncbi_wrapper loaded with {len(ncbi_wrapper.lineages)} taxid lineages.")
    taxids_to_use = establish_taxids_to_use(study_output_filepath, ncbi_wrapper, input_taxids, tax_level_to_use= tax_level_to_use, min_tax_count= taxa_threshold, normalize = True)
    logger.info(f"Number of taxids to use at level {tax_level_to_use}: {taxids_to_use.shape[0]}")
    logger.info(f"ncbi_wrapper loaded with {len(ncbi_wrapper.lineages)} taxid lineages.")
    ### Run data retrieval
    trainning_folders = folders[:int(len(folders) * proportion_train)]
    logger.info(f"Number of training datasets: {len(trainning_folders)}")
    trainning_results_df, prediction_trainning_results_df, recall_trainning_results = run_data_retrieval(trainning_folders, study_output_filepath, data_set_divide, ncbi_wrapper, input_tax_df, tax_level_to_use, taxids_to_use)
    ndata_sets = trainning_results_df['data_set'].nunique()
    logger.info(f"Number of training datasets retrieved: {ndata_sets}")
    ### Models
    from metagenomics_utils.overlap_manager.om_models import RecallModeller, CompositionModeller, CrossHitModeller
    recall_modeller= RecallModeller(recall_trainning_results= recall_trainning_results, data_set_divide= data_set_divide)
    composition_modeller= CompositionModeller(trainning_results_df= trainning_results_df)
    crosshit_modeller= CrossHitModeller(prediction_trainning_results_df= prediction_trainning_results_df)

    model_recall, X_test_recall, Y_test_recall =recall_modeller.train_model()
    model_composition, X_train_composition, X_test_composition, y_test_composition = composition_modeller.train_model()
    crosshit_modeller.train_model()

    recall_modeller.save_model(models_subdir)
    composition_modeller.save_model(models_subdir)
    crosshit_modeller.save_model(models_subdir)

    recall_modeller.model_summary(recall_modeller.model, X_test_recall, Y_test_recall, analysis_output_filepath)

    composition_modeller.eval_and_plot(X_test_composition, y_test_composition, analysis_output_filepath, X_train=X_train_composition)

    #### Cross-Hit Analysis
    remaining_folders = folders[int(len(folders) * proportion_train):]
    test_results_df, summary_results_df, trash_results_df, cross_hit_results_df = compound_eda_function(
        remaining_folders,
        study_output_filepath, 
        data_set_divide,
        ncbi_wrapper, 
        input_tax_df, 
        tax_level_to_use, 
        taxids_to_use, 
        recall_modeller, 
        composition_modeller,
        crosshit_modeller,
        cross_hit_threshold
    )

    taxids_to_use.to_csv(os.path.join(analysis_output_filepath, "taxids_to_use.tsv"), sep="\t", index=False)
    test_results_df.to_csv(os.path.join(analysis_output_filepath, "test_datasets_overall_precision.tsv"), sep="\t", index=False)
    trash_results_df.to_csv(os.path.join(analysis_output_filepath, "test_datasets_trash_composition.tsv"), sep="\t", index=False)
    cross_hit_results_df.to_csv(os.path.join(analysis_output_filepath, "test_datasets_cross_hit_composition.tsv"), sep="\t", index=False)

    data_set_summary_results = summary_results_df.copy()
    data_set_summary_results = data_set_summary_results.drop_duplicates(subset=['sample']).copy()
    nsets_analysed = summary_results_df['sample'].nunique()
    logger.info(f"Number of test datasets analysed: {nsets_analysed}")

    summary_stats_precision = ['overall_precision_raw','fuzzy_precision_raw', 'fuzzy_precision_cov_filtered', 'clade_precision_full', 'clade_precision_post']
    df = data_set_summary_results[summary_stats_precision]
    stats = df.describe().T
    stats.to_csv(os.path.join(analysis_output_filepath, "precision_summary_statistics.tsv"), sep="\t")

    plt.figure(figsize=(10, 6))
    sns.histplot(test_results_df['overall_precision'], bins=20, kde=True)
    plt.xlabel('Overall Precision')
    plt.xlim(0, 2)
    plt.ylabel('Frequency')
    plt.title('Distribution of Overall Precision Across Test Datasets')
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_output_filepath, "overall_precision_histogram.png"))
    plt.close()

    plt.figure(figsize=(12, 8))
    melted_df = data_set_summary_results.melt(id_vars=['sample'], value_vars=['overall_precision_raw','fuzzy_precision_raw', 'fuzzy_precision_cov_filtered', 'clade_precision_full', 'clade_precision_post'], var_name='Metric', value_name='Value')
    sns.boxplot(x='Metric', y='Value', data=melted_df)
    plt.title('Comparison of Precision Metrics Across Datasets')
    plt.ylabel('Precision')
    plt.xlabel('Metric')
    plt.ylim(0, 3)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_output_filepath, "precision_metrics_boxplot.png"))
    plt.close()

    plt.figure(figsize=(12, 8))
    melted_df = data_set_summary_results.melt(id_vars=['sample'], value_vars=['overall_precision_raw','fuzzy_precision_raw', 'fuzzy_precision_cov_filtered', 'clade_precision_full', 'clade_precision_post'], var_name='Metric', value_name='Value')
    sns.histplot(data=melted_df, x='Value', hue='Metric', element='step', stat='density', common_norm=False, bins=20)
    plt.title('Comparison of Precision Metrics Across Datasets')
    plt.ylabel('Precision')
    plt.xlabel('Metric')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_output_filepath, "precision_metrics_histogram.png"))
    plt.close()

    plt.figure(figsize=(12, 6))
    recall_summary = data_set_summary_results[['sample', 'recall_raw', 'recall_cov_filtered', 'clade_recall', 'recall_filtered_leaves']].copy()
    recall_summary['recall_clade_diff'] = recall_summary['clade_recall'] - recall_summary['recall_raw']
    recall_summary['recall_clade_diff_predicted_leaves'] = recall_summary['recall_filtered_leaves'] - recall_summary['recall_raw']
    sns.histplot(recall_summary['recall_clade_diff'], bins=20, kde=True)
    sns.histplot(recall_summary['recall_clade_diff_predicted_leaves'], bins=20, kde=True, color='orange')
    plt.title('Distribution of Recall Improvement (Clade Recall - Raw Recall)')
    plt.xlabel('Recall Improvement')
    plt.ylabel('Frequency')
    plt.axvline(0, color='red', linestyle='--', label='No Improvement')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_output_filepath, "recall_improvement_histogram.png"))
    plt.close()

    plt.figure(figsize=(12, 6))
    melted_df = data_set_summary_results.melt(id_vars=['sample'], value_vars=['recall_raw', 'recall_cov_filtered' ,'clade_recall', 'recall_filtered_leaves'], var_name='Metric', value_name='Value')
    sns.boxplot(x='Metric', y='Value', data=melted_df)
    plt.title('Comparison of Recall Metrics Across Datasets')
    plt.ylabel('Recall')
    plt.xlabel('Metric')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_output_filepath, "recall_metrics_boxplot.png"))
    plt.close()

    precisions_df = summary_results_df[['sample', 'recall_raw', 'recall_cov_filtered', 'clade_recall','fuzzy_precision_raw', 'fuzzy_precision_cov_filtered','overall_precision_raw', 'clade_precision_full', 'clade_precision_post']].drop_duplicates()

    precisions_df.loc[:, 'Prob_Find_any'] = precisions_df['recall_raw'] * precisions_df['fuzzy_precision_raw']
    precisions_df.loc[:, 'Prob_Find_any_cov'] = precisions_df['recall_cov_filtered'] * precisions_df['fuzzy_precision_cov_filtered']
    precisions_df.loc[:, 'Prob_Find_true'] = precisions_df['recall_raw'] * precisions_df['overall_precision_raw']
    precisions_df.loc[:, 'Prob_Find_true_clade_full'] = precisions_df['clade_recall'] * precisions_df['clade_precision_full']

    precisions_df_melt = precisions_df.melt(id_vars=['sample'], value_vars=['Prob_Find_any', 'Prob_Find_any_cov', 'Prob_Find_true', 'Prob_Find_true_clade_full'], var_name='Metric', value_name='Value')
    plt.figure(figsize=(12, 8))
    sns.boxplot(x='Metric', y='Value', data=precisions_df_melt)
    plt.title('Comparison of Probability Metrics Across Datasets')
    plt.ylabel('Probability')
    plt.xlabel('Metric')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_output_filepath, "probability_metrics_boxplot.png"))
    plt.close()

    plt.figure(figsize=(10, 8))
    plt.ylabel('Raw Precision (Post)')
    sns.boxplot(x=tax_level_to_use, y='raw_pred_accuracy', data=summary_results_df)
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(analysis_output_filepath, "clade_precision.png"))
    plt.close()

    def composition_summary(composition_df: pd.DataFrame):
        summary_list = []
        for tax_level in composition_df['tax_level'].unique():
            subset = composition_df[composition_df['tax_level'] == tax_level]
            mean_values = subset.drop(columns=['taxid', 'tax_level', 'data_set']).mean()
            mean_values = pd.DataFrame(mean_values).T
            
            std_values = subset.drop(columns=['taxid', 'tax_level']).std()
            mean_values.insert(0, 'tax_level', tax_level)
            summary_list.append(mean_values)
        summary_list = pd.concat(summary_list, ignore_index=True)
        summary_list.set_index('tax_level', inplace=True)
        return summary_list

    cross_hit_composition_summary_df = composition_summary(cross_hit_results_df)
    # sort by tax_level - columns and index
    cross_hit_composition_summary_df = cross_hit_composition_summary_df.sort_index()
    cross_hit_composition_summary_df = cross_hit_composition_summary_df.reindex(sorted(cross_hit_composition_summary_df.columns), axis=1)
    # plot heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(cross_hit_composition_summary_df, cmap='viridis', annot=True, fmt=".2f")
    plt.title('Average Cross-Hit Composition by Tax Level')
    plt.xlabel('Taxa')
    plt.ylabel('Tax Level')
    plt.savefig(os.path.join(analysis_output_filepath, "cross_hit_composition_heatmap.png"))
    plt.close()


    def composition_summary(composition_df: pd.DataFrame):
        summary_list = []
        for tax_level in composition_df['tax_level'].unique():
            subset = composition_df[composition_df['tax_level'] == tax_level]
            mean_values = subset.drop(columns=['taxid', 'tax_level', 'data_set']).mean()
            mean_values = pd.DataFrame(mean_values).T
            
            mean_values.insert(0, 'tax_level', tax_level)
            summary_list.append(mean_values)
        summary_list = pd.concat(summary_list, ignore_index=True)
        summary_list.set_index('tax_level', inplace=True)
        return summary_list

    trash_composition_summary_df = composition_summary(trash_results_df)

    # sort by tax_level - columns and index
    trash_composition_summary_df = trash_composition_summary_df.sort_index()
    trash_composition_summary_df = trash_composition_summary_df.reindex(sorted(trash_composition_summary_df.columns), axis=1)
    # plot heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(trash_composition_summary_df, cmap='viridis', annot=True, fmt=".4f")
    plt.title('Average Trash Composition by Tax Level')
    plt.xlabel('Taxa')
    plt.ylabel('Tax Level')
    plt.savefig(os.path.join(analysis_output_filepath, "trash_composition_heatmap.png"))
    plt.close()