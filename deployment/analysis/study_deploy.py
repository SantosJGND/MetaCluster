import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from metagenomics_utils.overlap_manager import OverlapManager
from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper
import os
from metagenomics_utils.overlap_manager.om_models import (
    data_set_traversal_with_precision,
    cross_hit_prediction_matrix,
    get_m_stats_matrix,
    predict_recall_cutoff_vars,
    get_trash_composition,
    get_cross_hit_composition,
    calculate_overall_precision,
    predict_data_set_clades,
    cross_hit_prediction,
    cut_off_recall_prediction,
    CompositionModeller,
    RecallModeller,
    CrossHitModeller,
)
import numpy as np
import logging

def retrieve_simulation_input(study_output_filepath: str) -> pd.DataFrame:
    input_files = []
    folders = [f for f in os.listdir(study_output_filepath) if os.path.isdir(os.path.join(study_output_filepath, f))]

    for data_set_name in folders:

        output_filepath = os.path.join(study_output_filepath, f"{data_set_name}", "input", f"{data_set_name}.tsv")
        df = pd.read_csv(output_filepath, sep="\t")

        df['data_set'] = data_set_name.replace("_plan", "")
        if df.empty:
            continue
        input_files.append(df)

    
    all_input_data = pd.concat(input_files, ignore_index=True)

    return all_input_data



def process_clade_report(data_set_name, study_output_filepath):

    clade_report = os.path.join(study_output_filepath, f"{data_set_name}", "output", "clade_report_with_references.tsv")
    if not os.path.exists(clade_report):
        return []

    clade_report = pd.read_csv(clade_report, sep="\t").rename(columns={"#rname": "assembly_accession"})

    taxids = clade_report['taxid'].dropna().unique().tolist() 
    return taxids


def output_parse(study_output_filepath: str, ncbi_wrapper: NCBITaxonomistWrapper):
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

def establish_taxids_to_use(study_output_filepath: str, ncbi_wrapper: NCBITaxonomistWrapper, tax_level_to_use= 'order', min_tax_count= 0.02, normalize = True):
    
    output_taxids_table = output_parse(study_output_filepath, ncbi_wrapper)
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
            recall_stats = predict_recall_cutoff_vars(data_set_divide, data_set_name, m_stats_stats_matrix, input_tax_df)
            if not result_df.empty:
                trainning_results.append(result_df)
            if not prediction_matrix.empty:
                prediction_trainning_results.append(prediction_matrix)
            if not recall_stats.empty:
                recall_trainning_results.append(recall_stats)
        except Exception as e:
            print(f"s: {e}")
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
                            cross_hit_modeller: CrossHitModeller):

    test_results = []
    # cross-hit threshold
    cross_hit_threshold = 0.9
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
                                                        m_stats_stats_matrix, input_tax_df)
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

            if cross_hit_predictions.empty is False:

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
            print(f"s: {e}")
            print(f"Processing {data_set_name}")
            import traceback
            traceback.print_exc()

    summary_results_df = pd.concat(summary_results, ignore_index=True)
    summary_results_df.to_csv(os.path.join(study_output_filepath, "test_datasets_summary_results.tsv"), sep="\t", index=False)
    test_results_df = pd.DataFrame(test_results, columns=['overall_precision', 'data_set'])
    trash_results_df = pd.concat(trash_results, ignore_index=True)
    cross_hit_results_df = pd.concat(cross_hit_results, ignore_index=True)

    return test_results_df, summary_results_df, trash_results_df, cross_hit_results_df

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

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    study_output_filepath = args.study_output_filepath
    taxid_plan_filepath = args.taxid_plan_filepath
    analysis_output_filepath = args.analysis_output_filepath
    threshold = args.threshold
    taxa_threshold = args.taxa_threshold
    tax_level_to_use = args.tax_level_to_use
    data_set_divide = args.data_set_divide
    holdout_proportion = args.holdout_proportion
    proportion_train = 1 - holdout_proportion

    # Example usage:
    # study_output_filepath = "/path/to/study_output"
    # taxid_plan_filepath = "/path/to/taxid_plan.tsv"
    # analysis_output_filepath = "/path/to/output"
    #threshold = 0.3
    #taxa_threshold = 0.02
    #tax_level_to_use = 'order'
    #data_set_divide = 5

    #### output logger to file
    os.makedirs(analysis_output_filepath, exist_ok=True)

    logging.basicConfig(level=logging.INFO, filename=os.path.join(analysis_output_filepath, "study_deploy.log"), filemode='w',
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

    all_input_data = retrieve_simulation_input(study_output_filepath)
    ncbi_wrapper = NCBITaxonomistWrapper(db=output_db)
    input_taxids = all_input_data['taxid'].dropna().unique().tolist()
    ncbi_wrapper.resolve_lineages(input_taxids)

    #### load output taxids
    input_tax_df = expand_input_data(all_input_data, ncbi_wrapper, taxid_plan)
    logger.info(f"Number of unique input taxids: {input_tax_df.shape[0]}")
    logger.info(f"ncbi_wrapper loaded with {len(ncbi_wrapper.lineages)} taxid lineages.")
    taxids_to_use = establish_taxids_to_use(study_output_filepath, ncbi_wrapper, tax_level_to_use= tax_level_to_use, min_tax_count= taxa_threshold, normalize = True)
    logger.info(f"Number of taxids to use at level {tax_level_to_use}: {taxids_to_use.shape[0]}")
    logger.info(f"ncbi_wrapper loaded with {len(ncbi_wrapper.lineages)} taxid lineages.")
    ### Run data retrieval
    trainning_folders = folders[:int(len(folders) * proportion_train)]
    logger.info(f"Number of training datasets: {len(trainning_folders)}")
    trainning_results_df, prediction_trainning_results_df, recall_trainning_results = run_data_retrieval(trainning_folders, study_output_filepath, data_set_divide, ncbi_wrapper, input_tax_df, tax_level_to_use, taxids_to_use)
    ndata_sets = trainning_results_df['data_set'].nunique()
    logger.info(f"Number of training datasets retrieved: {ndata_sets}")
    ### Models
    recall_modeller= RecallModeller(recall_trainning_results= recall_trainning_results, data_set_divide= data_set_divide)
    composition_modeller= CompositionModeller(trainning_results_df= trainning_results_df)
    crosshit_modeller= CrossHitModeller(prediction_trainning_results_df= prediction_trainning_results_df)

    model_recall, X_test_recall, Y_test_recall =recall_modeller.train_model()
    model_composition, X_train_composition, X_test_composition, y_test_composition = composition_modeller.train_model()
    crosshit_modeller.train_model()

    recall_modeller.model_summary(recall_modeller.model, X_test_recall, Y_test_recall, analysis_output_filepath)
    print("Recall model evaluation completed.")
    print(X_train_composition.shape, X_test_composition.shape, y_test_composition.shape)
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
        crosshit_modeller
    )


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
            mean_values = subset.drop(columns=['taxid', 'tax_level']).mean()
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


    def composition_summarylk(composition_df: pd.DataFrame):
        summary_list = []
        for tax_level in composition_df['tax_level'].unique():
            subset = composition_df[composition_df['tax_level'] == tax_level]
            mean_values = subset.drop(columns=['taxid', 'tax_level']).mean()
            mean_values = pd.DataFrame(mean_values).T
            
            std_values = subset.drop(columns=['taxid', 'tax_level']).std()
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