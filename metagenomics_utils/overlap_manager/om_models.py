

from typing import List, Optional

from sklearn.discriminant_analysis import StandardScaler
from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.multioutput import MultiOutputRegressor
from metagenomics_utils.overlap_manager import OverlapManager
import pandas as pd
from metagenomics_utils.overlap_manager.node_stats import node_composition_level, node_leaf_shannon_tax_diversity, node_total_true_leaves
from metagenomics_utils.overlap_manager.node_stats import node_leaves_best_taxids, node_leaf_shannon_tax_diversity, get_subset_composition
import os
import numpy as np
from metagenomics_utils.overlap_manager.node_stats import get_m_stats_matrix, normalize_by_taxlevel, get_composition_by_leaf
from metagenomics_utils.overlap_manager.node_stats import get_subset_composition_counts
from metagenomics_utils.overlap_manager.diversity import shannon_diversity_from_list, skewness, kurtosis
import joblib

########################################################################################################
########################################################################################################
######################################################################################
### DATA EXTRACTION ######################################################  ##########
###################### PRECISION-BASED TRAVERSAL #####################################


def traversal_with_precision(
        overlap_manager: OverlapManager, 
        node: str, 
        m_stats_stats_matrix, 
        input_taxa: pd.DataFrame, 
        tax_level: str = "order", 
        results = [],
        force_stop = False) -> List[pd.DataFrame]:
    """
    Recursive function.
    Traverse the tree, internal nodes only. At each node:
    - compute node precision.
    - compute node composition at tax_level
    - extract node Min_Dist and Min_Shared
    - compute precision of split children.
    - determine if split increases precision.
    - store results. 
    if precision is increased by splitting, traverse (internal nodes only) children. else stop.
    """

    composition = node_composition_level(overlap_manager, node, m_stats_stats_matrix, input_taxa, tax_level=tax_level).set_index('tax_level').T
    composition.reset_index(drop=True, inplace=True)
    node_true_leaves = node_total_true_leaves(overlap_manager, node, m_stats_stats_matrix)
    node_precision = 1 / len(set(node_true_leaves)) if len(node_true_leaves) > 0 else 0.0
    node_precision = 1 / node_precision if node_precision > 1 else node_precision
    node_row = overlap_manager.all_node_stats[overlap_manager.all_node_stats['Node'] == node]
    min_dist = node_row['Min_Pairwise_Dist'].values[0]
    min_shared = node_row['Min_Shared'].values[0]
    node_total_leaf_taxa_div = node_leaf_shannon_tax_diversity(overlap_manager, node, m_stats_stats_matrix, tax_level=tax_level)

    node_children = list(overlap_manager.tree.successors(node))
    new_precision = len(node_children) / len(set(node_true_leaves)) if len(node_true_leaves) > 0 else 0.0
    new_precision = 1 / new_precision if new_precision > 1 else new_precision

    precision_increased = new_precision > node_precision

    # stop conditions
    stop_traversal = not precision_increased
    ## forcing stop conditions
    if force_stop:
        if min_dist == 1.0:
            stop_traversal = True
        if min_dist == 0.0:
            stop_traversal = False

    local_results = pd.DataFrame({
        'node': [node],
        'n_leaves': [len(overlap_manager.get_node_leaves(node))],
        'tax_diversity': [node_total_leaf_taxa_div],
        'n_true_leaves': [len(set(node_true_leaves))],
        'precision': [node_precision],
        'new_precision': [new_precision],
        'precision_increased': [precision_increased],
        'Min_Dist': [min_dist],
        'Min_Shared': [min_shared],
        'stop_traversal': [stop_traversal],
    })
    local_results = pd.concat([local_results, composition], axis=1)
    results.append(local_results)
    for child in node_children:
        if overlap_manager.tree.out_degree(child) > 0:  # internal node
            traversal_with_precision(overlap_manager, child, m_stats_stats_matrix, input_taxa, tax_level=tax_level, results=results)
    return results


def data_set_traversal_with_precision(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager: OverlapManager, input_taxa: pd.DataFrame, tax_level: str = "order"):
    m_stats_stats_matrix = get_m_stats_matrix(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager)
    results = []
    for root in overlap_manager.root_nodes:

        root_results = traversal_with_precision(overlap_manager, root, m_stats_stats_matrix, input_taxa, tax_level=tax_level, results=[])
        root_results = [r for r in root_results if r.empty == False]
        results.extend(root_results)

    if len(results) == 0:
        return pd.DataFrame()

    results_df = pd.concat(results, ignore_index=True)

    results_df.insert(0, 'data_set', data_set_name)
    return results_df

#####################################################################################
################ CROSS-HIT ANALYSIS ######################################################
######################################################################################


def cross_hit_prediction_matrix(data_set_name, 
                                study_output_filepath, 
                                ncbi_wrapper, 
                                overlap_manager: OverlapManager, 
                                tax_df, tax_level: str = 'order'):
    """
    For each leaf, compute the number of different taxids predicted at tax_level.
    """
    prediction_matrix = get_m_stats_matrix(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager)

    if prediction_matrix.empty:
        return pd.DataFrame()
    
    if len(overlap_manager.leaves) == 0:
        return pd.DataFrame()
    
    prediction_matrix = normalize_by_taxlevel(prediction_matrix, tax_level=tax_level)

    if prediction_matrix.empty:
        return pd.DataFrame()

    composition = get_composition_by_leaf(overlap_manager, prediction_matrix, tax_df, tax_level=tax_level)
    prediction_matrix.loc[:,'is_trash'] = (prediction_matrix['best_match_is_best'] == False) & (prediction_matrix['is_crosshit'] == False) 
    prediction_matrix_stats = prediction_matrix[['leaf', 'is_trash','coverage', 'covbases', 'meanmapq', 'error_rate', 'max_shared', 'total_uniq_reads']].copy()
    
    prediction_matrix_stats = prediction_matrix_stats.merge(composition, on='leaf', how='left')
    prediction_matrix_stats = prediction_matrix_stats[prediction_matrix_stats['leaf'].notna()]

    return prediction_matrix_stats


######################################################################################
################ RECALL CUTOFF #######################################################################
######################################################################################


def predict_recall_cutoff_vars(data_set_divide:int, data_set_name: str, m_stats_stats_matrix: pd.DataFrame, input_tax_df: pd.DataFrame, tax_level: str = "order") -> pd.DataFrame:
    """
    Predict recall at various cutoffs and other composition statistics.
    """
    # index of last True valuie in best_match_is_best
    m_stats_stats_matrix = m_stats_stats_matrix.sort_values(by="total_uniq_reads", ascending=False).reset_index(drop=True)
    best_match_indices = m_stats_stats_matrix.index[m_stats_stats_matrix['best_match_is_best'] == True].tolist()
    last_best_match_index = best_match_indices[-1] + 1 if best_match_indices else -1 
    last_best_match_relindex = (last_best_match_index) / len(m_stats_stats_matrix) if len(m_stats_stats_matrix) > 0 else 0.0

    composition = get_subset_composition_counts(m_stats_stats_matrix, input_tax_df, tax_level=tax_level)[['tax_level', 'proportion']].set_index('tax_level').T


    tax_diversity_shannon = shannon_diversity_from_list(m_stats_stats_matrix['order'].dropna().tolist())

    counts_skewness = skewness(m_stats_stats_matrix['total_uniq_reads'].tolist())
    counts_kurtosis = kurtosis(m_stats_stats_matrix['total_uniq_reads'].tolist())

    composition.insert(0, 'tax_diversity_shannon', tax_diversity_shannon)
    composition.insert(0, 'counts_skewness', counts_skewness)
    composition.insert(0, 'counts_kurtosis', counts_kurtosis)

    for set_divide in reversed(range(1, data_set_divide + 1)):
        threshold_index = int(len(m_stats_stats_matrix) * set_divide / data_set_divide)
        best_match_indices_short = [idx for idx in best_match_indices if idx < threshold_index]
        index_recall = len(best_match_indices_short) / len(best_match_indices) if len(best_match_indices) > 0 else 0.0
        composition.insert(0, f'index_recall_{set_divide}', index_recall)

    composition.insert(0, 'last_best_match_relindex', last_best_match_relindex)
    composition.insert(0, 'data_set', data_set_name)

    composition.reset_index(drop=True, inplace=True)

    #composition.drop(columns=['tax_level'], inplace=True)
    return composition




########################################################################################################
########################################################################################################
######################################################################################
### PREDICTION ######################################################  ##########
###################### RECALL #####################################


class RecallModeller:

    model_save_filename = "recall_xgb_bundle.pkl"

    def __init__(self, recall_trainning_results, data_set_divide: int):
        self.recall_trainning_results = recall_trainning_results.drop(columns=['data_set'])
        self.data_set_divide = data_set_divide

        self.RecP_target_cols = [
            'last_best_match_relindex',
        ] + [f'index_recall_{i}' for i in range(1, data_set_divide + 1)]

        self.RecP_feature_cols = self.recall_trainning_results.columns.difference(self.RecP_target_cols).tolist()
        self.model: Optional[MultiOutputRegressor] = None
        self.X_test = None
        self.y_test = None


    def prep_data(self):
        X = self.recall_trainning_results[self.RecP_feature_cols]
        Y = self.recall_trainning_results[self.RecP_target_cols]
        return X, Y

    def split_data(self, test_size=0.2, random_state=42):
        from sklearn.model_selection import train_test_split
        X, Y = self.prep_data()
        X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, random_state=random_state)
        return X_train, X_test, Y_train, Y_test


    def multioutput_regressor(self, X_train, Y_train):

        rf = RandomForestRegressor(n_estimators=100, random_state=42)
        multi_rf = MultiOutputRegressor(rf)
        multi_rf.fit(X_train, Y_train)
        return multi_rf
    
    def train_model(self):
        X_train, X_test, Y_train, Y_test = self.split_data()
        model = self.multioutput_regressor(X_train, Y_train)
        self.model = model
        self.X_test = X_test
        self.y_test = Y_test
        return model, X_test, Y_test

    def save_model(self, output_directory: str):
        if self.model is not None:
            joblib.dump({
                'model': self.model,
                'feature_names': self.RecP_feature_cols,
                'target_names': self.RecP_target_cols
            }, os.path.join(output_directory, self.model_save_filename))
        else:
            print("No model to save.")

    def load_model(self, input_directory: str):
        try:
            bundle = joblib.load(os.path.join(input_directory, self.model_save_filename))
            self.model = bundle['model']
            self.RecP_feature_cols = bundle['feature_names']
            self.RecP_target_cols = bundle['target_names']
        except Exception as e:
            print(f"Error loading model: {e}")

    def evaluate_model(self, model, X_test, Y_test, ouptput_filepath):
        from sklearn.metrics import r2_score, mean_squared_error

        Y_pred = model.predict(X_test)
        r2_scores = {}
        mse_scores = {}
        with open(ouptput_filepath, 'w') as f:
            f.write("\t".join(["Target_Column", "R2_Score", "MSE"]) + "\n")
            for i, col in enumerate(self.RecP_target_cols):
                r2 = r2_score(Y_test.iloc[:, i], Y_pred[:, i])
                mse = mean_squared_error(Y_test.iloc[:, i], Y_pred[:, i])
                r2_scores[col] = r2
                mse_scores[col] = mse
                f.write(f"{col}\t{r2}\t{mse}\n")
        return r2_scores, mse_scores
    
    def feature_importances(self, model, output_filepath):
        importances = np.mean([est.feature_importances_ for est in model.estimators_], axis=0)
        feat_importance = pd.Series(importances, index=self.RecP_feature_cols).sort_values(ascending=False)
        feat_importance.to_csv(output_filepath)

    def plot_eval(self, X_test, Y_test, analysis_output_filepath):
        import matplotlib.pyplot as plt
        Y_pred = self.model.predict(X_test)
        differences = Y_test.values - Y_pred
        avg_differences = differences.mean(axis=0)

        # plot distribution of differences
        plt.boxplot(differences, labels=self.RecP_target_cols)
        plt.title("Distribution of Differences between True and Predicted Recall")
        plt.ylabel("Difference")
        plt.savefig(analysis_output_filepath)
        plt.close()


    def model_summary(self, model, X_test, Y_test, analysis_output_filedir):
        Y_pred = model.predict(X_test)
        from sklearn.metrics import r2_score, mean_squared_error
        r2 = r2_score(Y_test, Y_pred, multioutput='uniform_average')
        mse = mean_squared_error(Y_test, Y_pred, multioutput='uniform_average')

        print(f"model_summary R² = {r2:.3f}, MSE = {mse:.3f}")

        analysis_output_filepath = os.path.join(analysis_output_filedir, "recall_model_analysis_results.txt")
        r2_scores, mse_scores = self.evaluate_model(model, X_test, Y_test, analysis_output_filepath)
        feat_importance_filepath = analysis_output_filepath.replace(".txt", "_feature_importances.tsv")
        self.feature_importances(model, feat_importance_filepath)
        self.plot_eval(X_test, Y_test, analysis_output_filepath.replace(".txt", "_recall_prediction_differences.png"))
        import matplotlib.pyplot as plt

        for i in range(2,4):  # first 3 test samples
            plt.plot(range(1, 6), Y_test.iloc[i, 1:], 'o-', label='True')
            plt.plot(range(1, 6), Y_pred[i, 1:], 's--', label='Predicted')
            plt.title(f"Sample {i}: recall curve")
            plt.xlabel("Recall index (1–5)")
            plt.ylabel("Recall value")
            plt.legend()
            plt.savefig(analysis_output_filepath.replace(".txt", f"_recall_curve_sample_{i}.png"))
            plt.close()

        return r2_scores, mse_scores



def cut_off_recall_prediction(study_output_filepath: str, data_set_name: str, modeller: RecallModeller, data_set_divide:int, m_stats_stats_matrix: pd.DataFrame, taxids_to_use: pd.DataFrame, tax_level = 'order') -> OverlapManager:
    """
    Predict recall at various cutoffs and filter leaves based on threshold.
    """
    recall_stats = predict_recall_cutoff_vars(data_set_divide, data_set_name, m_stats_stats_matrix, taxids_to_use, tax_level=tax_level)
    recall_pred = modeller.model.predict(recall_stats[modeller.RecP_feature_cols])
    recall_pred_df = pd.DataFrame(recall_pred, columns=modeller.RecP_target_cols)
    keep_index = recall_pred_df.iloc[0]['last_best_match_relindex'] * m_stats_stats_matrix.shape[0] 
    keep_index = round(keep_index + 1)
    #print(f"Keeping top {keep_index} leaves based on predicted last_best_match_relindex.")
    overlap_manager = OverlapManager(os.path.join(study_output_filepath, f"{data_set_name}", "clustering"), max_taxids=keep_index)

    return overlap_manager


######################################################################################
################ NCBI TAXONOMIST UTILITIES ###########################################
######################################################################################

class ScalerProxy(StandardScaler):

    def __init__(self):
        super().__init__()

    def fit_transform(self, x_data: pd.DataFrame):
        return x_data

    def inverse_transform(self, x_data: pd.DataFrame):
        return x_data

    def fit(self, x_data: pd.DataFrame):
        pass


class CompositionModeller:

    model_save_filename = "composition_xgb_bundle.pkl"

    def __init__(self, trainning_results_df: pd.DataFrame):
        self.trainning_results_df = trainning_results_df
        self.X = trainning_results_df.drop(columns=['data_set', 'node', 'n_true_leaves', 'precision_increased', 'new_precision', 'precision', 'stop_traversal', 'unclassified'])
        self.X_stats_cols= ['n_leaves', 'tax_diversity', 'Min_Dist', 'Min_Shared']
        self.X_tax_cols = self.X.columns.difference(self.X_stats_cols).tolist()

        self.y = trainning_results_df['stop_traversal'].astype(int)
        self.X_test = None
        self.X_train = None
        self.y_test = None
        self.scaler = ScalerProxy()
        self.pca = None
        self.model = None
    
    def split_data(self, test_size=0.2, random_state=42):
        from sklearn.model_selection import train_test_split
        X_train, X_test, y_train, y_test = train_test_split(self.X, self.y, test_size=test_size, random_state=random_state)
        return X_train, X_test, y_train, y_test
    
    def prep_data(self, scale: bool= True, pca_transform: bool = False):

        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA
        X_train, X_test, y_train, y_test = self.split_data()

        if scale:
            scaler = StandardScaler()
            X_train[self.X_stats_cols] = scaler.fit_transform(X_train[self.X_stats_cols])
            X_test[self.X_stats_cols] = scaler.transform(X_test[self.X_stats_cols])
            self.scaler = scaler

        if pca_transform:
            pca = PCA(n_components=0.95)
            X_train_pca = pca.fit_transform(X_train)
            X_test_pca = pca.transform(X_test)
            self.pca = pca
            return X_train_pca, X_test_pca, y_train, y_test

        return X_train, X_test, y_train, y_test

    def xgbc_model(self, X_train, y_train, **kwargs):
        from xgboost import XGBClassifier

        model = XGBClassifier(use_label_encoder=False, eval_metric='logloss', **kwargs)
        model.fit(X_train, y_train)
        return model
    
    def xgbc_model_bayes_optimized(self, X_train, y_train):
        import optuna
        from sklearn.model_selection import cross_val_score, StratifiedKFold
        from xgboost import XGBClassifier

        def objective(trial):
            params = {
                "objective": "binary:logistic",
                "eval_metric": "auc",
                "n_estimators": trial.suggest_int("n_estimators", 200, 800),
                "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.2, log=True),
                "max_depth": trial.suggest_int("max_depth", 3, 8),
                "subsample": trial.suggest_float("subsample", 0.6, 1.0),
                "colsample_bytree": trial.suggest_float("colsample_bytree", 0.6, 1.0),
                "reg_alpha": trial.suggest_float("reg_alpha", 1e-3, 10.0, log=True),
                "reg_lambda": trial.suggest_float("reg_lambda", 1e-3, 10.0, log=True),
                "gamma": trial.suggest_float("gamma", 0, 5),
                "min_child_weight": trial.suggest_int("min_child_weight", 1, 10),
                "random_state": 42,
                "n_jobs": -1
            }

            model = XGBClassifier(**params)
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            scores = cross_val_score(model, X_train, y_train, cv=cv, scoring="roc_auc", n_jobs=-1)
            return scores.mean()

        study = optuna.create_study(direction="maximize")
        study.optimize(objective, n_trials=50, show_progress_bar=True)

        # Refit model on full training data
        best_params = study.best_trial.params
        best_model = XGBClassifier(**best_params)
    
        best_model.fit(X_train, y_train)
        return best_model, study

    def train_model(self, bayes_optimized: bool = True, **kwargs) -> tuple:
        X_train, X_test, y_train, y_test = self.prep_data()
        self.X_train = X_train
        self.X_test = X_test
        self.y_test = y_test

        if bayes_optimized:
            model, _study = self.xgbc_model_bayes_optimized(X_train, y_train)
            self.model = model
            return model, X_train, X_test, y_test

        model = self.xgbc_model(X_train, y_train, **kwargs)
        self.model = model
        return model, X_train, X_test, y_test
    
    def save_model(self, output_directory: str):
        if self.model is not None:
            joblib.dump({
                'model': self.model, 
                'scaler': self.scaler,
                'taxa': self.X_tax_cols,
            }, os.path.join(output_directory, self.model_save_filename))
        else:
            print("No model to save.")
    
    def load_model(self, output_directory: str):
        try:
            bundle = joblib.load(os.path.join(output_directory, self.model_save_filename))
            self.model = bundle['model']
            self.scaler = bundle['scaler']
            self.X_tax_cols = bundle['taxa']
        except Exception as e:
            print(f"Error loading model: {e}")

    def evaluate_model(self, model, X_test, y_test):
        from sklearn.metrics import classification_report, confusion_matrix

        y_pred = model.predict(X_test)
        print(y_pred)
        print(y_test)
        report = classification_report(y_test, y_pred, output_dict=True)
        cm = confusion_matrix(y_test, y_pred)
        return report, cm
    
    def plot_eval(self, model, X_test, y_test, output_directory):
        import matplotlib.pyplot as plt
        import seaborn as sns

        # ========================
        # 5. Global feature importance (gain-based)
        # ========================

        xgb_importance = pd.Series(model.feature_importances_, index=self.X.columns)
        xgb_importance.sort_values(ascending=True).tail(20).plot.barh(figsize=(8, 8))
        plt.title("Top 20 XGBoost Feature Importances")
        plt.tight_layout()
        plt.savefig(f"{output_directory}/xgb_feature_importance.png")
        plt.close()

    def shap_eval_plot(self, model: XGBClassifier, X_train, X_val, output_directory):
        import shap
        import matplotlib.pyplot as plt

        # ========================
        # 6. SHAP Interpretation
        # ========================

        # Create SHAP explainer (tree-based method)
        explainer = shap.Explainer(model.predict, X_train)
        shap_values = explainer(X_val)

        # --- a) Global summary plot
        plt.figure(figsize=(10, 7))
        shap.summary_plot(shap_values, X_val, plot_type="dot", show=False)
        plt.savefig(f"{output_directory}/shap_summary_plot.png")
        plt.close()

        # --- b) Global bar plot (mean absolute SHAP values)
        plt.figure(figsize=(10, 7))
        shap.summary_plot(shap_values, X_val, plot_type="bar", show=False)
        plt.savefig(f"{output_directory}/shap_bar_plot.png")
        plt.close()

        # --- c) Individual feature dependence
        # Example: examine how Min_Shared influences precision_increased
        plt.figure(figsize=(10, 7))
        shap.dependence_plot("Min_Shared", shap_values.values, X_val, show=False)
        plt.savefig(f"{output_directory}/shap_dependence_plot.png")
        plt.close()

        # --- d) Inspect a single prediction
        # Pick an example (e.g., first sample)
        idx = 0
        plt.figure(figsize=(10, 7))
        shap.plots.waterfall(shap_values[idx], show=False)
        plt.savefig(f"{output_directory}/shap_waterfall_plot.png")
        plt.close()

    def shap_interaction_plot(self, model, X_train, X_val, output_directory: str):
        import shap
        import matplotlib.pyplot as plt
        # Explicitly create a TreeExplainer (important!)
        explainer = shap.TreeExplainer(model)

        # Compute SHAP interaction values (can be memory heavy for large data)
        interaction_values = explainer.shap_interaction_values(X_val)

        import numpy as np
        interaction_strength = np.abs(interaction_values).mean(axis=0)

        # Convert to DataFrame for readability
        interaction_df = pd.DataFrame(interaction_strength, index=X_train.columns, columns=X_train.columns)
        # lower triangle only
        interaction_df = interaction_df.where(np.tril(np.ones(interaction_df.shape), k=-1).astype(bool))
        import seaborn as sns
        plt.figure(figsize=(10, 8))
        sns.heatmap(interaction_df, cmap="viridis")
        plt.title("Mean absolute SHAP interaction values")
        plt.savefig(f"{output_directory}/shap_interaction_heatmap.png")
        interaction_df = interaction_df.drop(columns=self.X_stats_cols, index=self.X_stats_cols)
        distance_matrix = 1 - interaction_df
        # fill diagonal with 0
        np.fill_diagonal(distance_matrix.values, 0)
        distance_matrix = distance_matrix.fillna(0)
        # create tree from distance matrix
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import squareform
        import matplotlib.pyplot as plt

        #condensed_dist = squareform(distance_matrix.fillna(0).values)
        Z = linkage(distance_matrix, method='average')
        plt.figure(figsize=(6, 7))
        dendrogram(Z, labels=distance_matrix.index, orientation='left', leaf_rotation=0)
        plt.title("Feature Clustering based on SHAP Interaction Values")
        plt.xlabel("Distance")
        plt.ylabel("Features")
        plt.tight_layout()
        plt.savefig(f"{output_directory}/shap_interaction_dendrogram.png")
    
    def eval_and_plot(self, X_test, y_test, output_directory, X_train=None):
        report, cm = self.evaluate_model(self.model, X_test, y_test)
        self.plot_eval(self.model, X_test, y_test, output_directory)
        if X_train is not None:
            self.shap_eval_plot(self.model, X_train, X_test, output_directory)
            self.shap_interaction_plot(self.model, X_train, X_test, output_directory)
        return report, cm


#######################################################################################
################ CROSS-HIT MODELLING ######################################################
#######################################################################################

class CrossHitModeller:

    model_save_filename = "cross_hit_xgb_bundle.pkl"

    def __init__(self, prediction_trainning_results_df):
        self.prediction_trainning_results_df = prediction_trainning_results_df
        self.X = self.prediction_trainning_results_df.drop(columns=['leaf', 'is_trash'])
        self.pred_stats_cols = ['coverage', 'covbases', 'meanmapq', 'error_rate', 'max_shared', 'total_uniq_reads']
        self.y = self.prediction_trainning_results_df['is_trash'].astype(int)
        self.scaler = None
        self.pca = None
        self.model = None

    def split_data(self, test_size=0.2, random_state=42):
        from sklearn.model_selection import train_test_split
        X_train, X_test, y_train, y_test = train_test_split(self.X, self.y, test_size=test_size, random_state=random_state)
        return X_train, X_test, y_train, y_test

    def prep_data(self, scale: bool= True, transform: bool = False):

        from sklearn.preprocessing import StandardScaler
        X_train, X_test, y_train, y_test = self.split_data()

        if scale:
            scaler = StandardScaler()
            X_train[self.pred_stats_cols] = scaler.fit_transform(X_train[self.pred_stats_cols])
            X_test[self.pred_stats_cols] = scaler.transform(X_test[self.pred_stats_cols])
            self.scaler = scaler
        
        
        if transform:
            from sklearn.decomposition import PCA
            pca = PCA(n_components=0.95)
            X_train_pca = pca.fit_transform(X_train)
            X_test_pca = pca.transform(X_test)
            self.pca = pca
            return X_train_pca, X_test_pca, y_train, y_test

        return X_train, X_test, y_train, y_test

    def xgbc_model(self, X_train, y_train, **kwargs):
        from xgboost import XGBClassifier

        model = XGBClassifier(use_label_encoder=False, eval_metric='logloss', **kwargs)
        model.fit(X_train, y_train)
        return model
    
    def train_model(self, optimized: bool = True, **kwargs):
        if optimized:
            return self.train_model_bayes_optimized()
        X_train, X_test, y_train, y_test = self.prep_data()
        model = self.xgbc_model(X_train, y_train, **kwargs)
        self.model = model
        return model, X_test, y_test
    
    def train_model_bayes_optimized(self):
        import optuna
        from sklearn.model_selection import cross_val_score, StratifiedKFold
        from xgboost import XGBClassifier

        X_train, X_test, y_train, y_test = self.prep_data()

        def objective(trial):
            params = {
                "objective": "binary:logistic",
                "eval_metric": "auc",
                "n_estimators": trial.suggest_int("n_estimators", 200, 800),
                "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.2, log=True),
                "max_depth": trial.suggest_int("max_depth", 3, 8),
                "subsample": trial.suggest_float("subsample", 0.6, 1.0),
                "colsample_bytree": trial.suggest_float("colsample_bytree", 0.6, 1.0),
                "reg_alpha": trial.suggest_float("reg_alpha", 1e-3, 10.0, log=True),
                "reg_lambda": trial.suggest_float("reg_lambda", 1e-3, 10.0, log=True),
                "gamma": trial.suggest_float("gamma", 0, 5),
                "min_child_weight": trial.suggest_int("min_child_weight", 1, 10),
                "random_state": 42,
                "n_jobs": -1
            }

            model = XGBClassifier(**params)
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            scores = cross_val_score(model, X_train, y_train, cv=cv, scoring="roc_auc", n_jobs=-1)
            return scores.mean()

        study = optuna.create_study(direction="maximize")
        study.optimize(objective, n_trials=50, show_progress_bar=True)

        # Refit model on full training data
        best_params = study.best_trial.params
        best_model = XGBClassifier(**best_params)
    
        best_model.fit(X_train, y_train)
        self.model = best_model
        return best_model, X_test, y_test, study

    def save_model(self, output_directory: str):

        if self.model is not None:
            joblib.dump({
                'model': self.model,
                'scaler': self.scaler,
            }, os.path.join(output_directory, self.model_save_filename))
        else:
            print("No model to save.")

    def load_model(self, input_directory: str):

        try:
            bundle = joblib.load(os.path.join(input_directory, self.model_save_filename))
            self.model = bundle['model']
            self.scaler = bundle['scaler']
        except Exception as e:
            print(f"Error loading model: {e}")  


########################################################################################
################ TRAVERSAL ######################################################
########################################################################################

def cross_hit_prediction(data_set_name, 
                         study_output_filepath, 
                         ncbi_wrapper, 
                         modeller: CrossHitModeller,
                         overlap_manager: OverlapManager, 
                         tax_df, 
                         tax_level: str = 'order'):

    prediction_matrix = cross_hit_prediction_matrix(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager, tax_df, tax_level=tax_level)
    

    if prediction_matrix.empty or len(overlap_manager.leaves) == 0:
        return pd.DataFrame(columns=['leaf', 'is_trash', 'prob_best_match', 'pred_best_match'])

    
    X_pred = prediction_matrix.drop(columns=['leaf', 'is_trash'])
    pred_stats_cols = modeller.pred_stats_cols
    
    if modeller.scaler is not None:
        X_pred_scaled = X_pred
        X_pred_scaled[modeller.pred_stats_cols] = modeller.scaler.transform(X_pred_scaled[modeller.pred_stats_cols])
    else:
        X_pred_scaled = X_pred
    
    if modeller.pca is not None:
        X_pred_stats = X_pred_scaled[pred_stats_cols]
        X_pred_tax = X_pred_scaled.drop(columns=pred_stats_cols)
        X_pred_tax_pca = pd.DataFrame(modeller.pca.transform(X_pred_tax))
        X_pred_tax_pca.columns = [f'pca_{i+1}' for i in range(X_pred_tax_pca.shape[1])]
        X_pred_scaled = pd.concat([X_pred_stats.reset_index(drop=True), X_pred_tax_pca.reset_index(drop=True)], axis=1)
    
    if modeller.model is None:
        y_prob = np.zeros(X_pred_scaled.shape[0])
    else:
        y_prob = modeller.model.predict_proba(X_pred_scaled)[:, 1]

    output = prediction_matrix[['leaf', 'is_trash']].copy()
    output.loc[:, 'prob_best_match'] = y_prob
    output.loc[:, 'pred_best_match'] = (y_prob > 0.5).astype(int)

    return output

def traversal_with_prediction(overlap_manager: OverlapManager, node: str, modeller: CompositionModeller, m_stats_stats_matrix, tax_df: pd.DataFrame, tax_level: str = "order", results = []) -> List[pd.DataFrame]:
    """
    Recursive function.
    Traverse the tree, internal nodes only. At each node:
    - compute node precision.
    - compute node composition at tax_level
    - extract node Min_Dist and Min_Shared
    - compute precision of split children.
    - use model to predict if split increases precision.
    - store results. 
    if model predicts precision is increased by splitting, traverse (internal nodes only) children. else stop.
    """

    composition = node_composition_level(overlap_manager, node, m_stats_stats_matrix, tax_df, tax_level=tax_level).set_index('tax_level').T
    composition= composition.reset_index(drop=True)

    node_true_leaves = node_total_true_leaves(overlap_manager, node, m_stats_stats_matrix)
    node_precision = 1 / len(set(node_true_leaves)) if len(node_true_leaves) > 0 else 0.0
    node_precision = 1 / node_precision if node_precision > 1 else node_precision
    node_leaf_taxids = node_leaves_best_taxids(overlap_manager, node, m_stats_stats_matrix)
    node_row = overlap_manager.all_node_stats[overlap_manager.all_node_stats['Node'] == node]
    min_dist = node_row['Min_Pairwise_Dist'].values[0]
    min_shared = node_row['Min_Shared'].values[0]

    node_total_leaf_taxa_div = node_leaf_shannon_tax_diversity(overlap_manager, node, m_stats_stats_matrix, tax_level=tax_level)

    input_features = pd.DataFrame({
        'n_leaves': [len(overlap_manager.get_node_leaves(node))],
        'tax_diversity': [node_total_leaf_taxa_div],
        'Min_Dist': [min_dist],
        'Min_Shared': [min_shared],

    })

    input_features = pd.concat([input_features, composition], axis=1).drop(columns=['unclassified'], errors='ignore')
    if modeller.scaler is not None:
        input_stats = input_features[modeller.X_stats_cols]
        input_tax = input_features.drop(columns=modeller.X_stats_cols)
        input_stats = pd.DataFrame(modeller.scaler.transform(input_stats), columns=input_stats.columns)
        input_features = pd.concat([input_stats.reset_index(drop=True),input_tax.reset_index(drop=True)], axis=1)

    if modeller.pca is not None:
        input_tax = input_features.drop(columns=modeller.X_stats_cols)
        input_tax = modeller.pca.transform(input_tax)
        input_tax = pd.DataFrame(input_tax, columns=[f'pca_{i+1}' for i in range(input_tax.shape[1])])
        input_features = pd.concat([input_tax.reset_index(drop=True), input_features[modeller.X_stats_cols].reset_index(drop=True)], axis=1)

    stop_traversal_pred = modeller.model.predict(input_features)[0]
    if stop_traversal_pred is False:
        print(f"Stopping at node {node} with {len(overlap_manager.get_node_leaves(node))} leaves.")

    best_taxid_match = m_stats_stats_matrix[m_stats_stats_matrix['leaf'].isin(overlap_manager.get_node_leaves(node))].copy()
    best_taxid_match = best_taxid_match.sort_values(by=['best_match_score'], ascending=False)['best_match_taxid'].tolist()
    best_taxid_match = best_taxid_match[0] if len(best_taxid_match) > 0 else None

    if stop_traversal_pred == True or overlap_manager.tree.out_degree(node) == 0:  # stop condition or leaf
        node_leaves = overlap_manager.get_node_leaves(node)
        results.append({
            'node': node,
            'n_leaves': len(node_leaves),
            'leaves': node_leaves,
            'best_taxid_match': best_taxid_match,
            'node_precision': node_precision,
            'node_taxids': node_leaf_taxids,
        })

    # Traverse children if prediction is positive
    else:
        for child in overlap_manager.tree.successors(node):
            if overlap_manager.tree.out_degree(child) > 0:  # internal node
                traversal_with_prediction(overlap_manager, child, modeller, m_stats_stats_matrix, tax_df, tax_level=tax_level, results=results)
            else:
                best_taxid_match = m_stats_stats_matrix[m_stats_stats_matrix['leaf'] == child]['best_match_taxid'].tolist()
                best_taxid_match = best_taxid_match[0] if len(best_taxid_match) > 0 else None
                results.append(
                    {
                        'node': child,
                        'n_leaves': 1,
                        'leaves': [child],
                        'best_taxid_match': best_taxid_match,
                        'node_precision': 1.0,
                        'node_taxids': [best_taxid_match] if best_taxid_match is not None else [],
                    }
                )
        

    return results

def predict_data_set_clades(data_set_name, m_stats_stats_matrix, overlap_manager: OverlapManager, modeller: CompositionModeller, input_taxa: pd.DataFrame, tax_level: str = "order"):

    results = []
    for root in overlap_manager.root_nodes:

        root_results = traversal_with_prediction(overlap_manager, root, modeller, m_stats_stats_matrix, input_taxa, tax_level=tax_level, results=[])
        results.extend(root_results)

    if len(results) == 0:
        return pd.DataFrame(columns = ['data_set', 'node', 'n_leaves', 'leaves', 'best_taxid_match', 'node_precision', 'node_taxids'])

    results_df = pd.DataFrame(results)

    results_df.insert(0, 'data_set', data_set_name)
    return results_df


def calculate_overall_precision(result_df, m_stats_stats_matrix):
    true_taxids = m_stats_stats_matrix[(m_stats_stats_matrix['best_match_is_best']) | (m_stats_stats_matrix['is_crosshit'])].copy()
    true_taxids = true_taxids[(true_taxids['leaf'].isna() == False)]['best_match_taxid'].tolist()
    n_true = len(set(true_taxids))
    n_predicted = result_df.drop_duplicates(subset=['node', 'best_taxid_match']).shape[0]
    overall_precision = n_true / n_predicted if n_predicted > 0 else 0.0
    #overall_precision = 1 / overall_precision if overall_precision > 1 else overall_precision
    return overall_precision


def data_set_merged_output(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager: OverlapManager,  modeller: CompositionModeller, input_taxa: pd.DataFrame, tax_level: str = "order"):
    m_stats_stats_matrix = get_m_stats_matrix(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager)
    if m_stats_stats_matrix.empty:
        return pd.DataFrame(), 0.0
    results_df = predict_data_set_clades(data_set_name, m_stats_stats_matrix, overlap_manager, modeller, input_taxa, tax_level=tax_level)
    if results_df.empty:
        return pd.DataFrame(), 0.0
    overall_precision = calculate_overall_precision(results_df, m_stats_stats_matrix)

    return results_df, overall_precision

def get_trash_composition(input_df, m_stats_matrix: pd.DataFrame, tax_df: pd.DataFrame, tax_level: str = "order"):
    compositions = []
    
    for _, row in input_df.iterrows():
        taxid = row['taxid']
        tax = row[tax_level]
        trash_subset = m_stats_matrix[(m_stats_matrix['best_match_is_best'] == False) & (m_stats_matrix['is_crosshit'] == False)]
        trash_subset = trash_subset[trash_subset['cross_hit_match'] == taxid]
        subset_composition = get_subset_composition(trash_subset, tax_df, tax_level=tax_level).set_index('tax_level').T
        
        subset_composition = subset_composition.reset_index(drop=True)
        subset_composition.loc[:, 'tax_level'] = tax
        
        subset_composition.insert(0, 'taxid', taxid)

        compositions.append(subset_composition)
    return pd.concat(compositions, ignore_index=True, axis=0)


def get_cross_hit_composition(input_df, m_stats_matrix: pd.DataFrame, tax_df: pd.DataFrame, tax_level: str = "order"):
    compositions = []

    for _, row in input_df.iterrows():
        taxid = row['taxid']
        tax = row[tax_level]
        crosshit_subset = m_stats_matrix[m_stats_matrix['is_crosshit'] == True]
        crosshit_subset = crosshit_subset[crosshit_subset['cross_hit_match'] == taxid]
        subset_composition = get_subset_composition(crosshit_subset, tax_df, tax_level=tax_level).set_index('tax_level').T
        
        subset_composition = subset_composition.reset_index(drop=True)
        subset_composition.loc[:, 'tax_level'] = tax
        
        subset_composition.insert(0, 'taxid', taxid)

        compositions.append(subset_composition)
    return pd.concat(compositions, ignore_index=True, axis=0)
