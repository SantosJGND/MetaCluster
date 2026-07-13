from __future__ import annotations

import logging
import os
import warnings
from copy import deepcopy
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.base import BaseEstimator, ClassifierMixin, RegressorMixin
from sklearn.discriminant_analysis import StandardScaler
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

from datetime import datetime, timezone

from metagenomics_utils.overlap_manager.node_stats import (
    get_composition_by_leaf,
    get_m_stats_matrix,
    get_subset_composition,
    node_composition_level,
    node_composition_with_stats,
    node_leaf_shannon_tax_diversity,
    node_leaves_best_taxids,
    node_total_true_leaves,
    normalize_by_taxlevel,
)

logger = logging.getLogger(__name__)

########################################################################################################
########################################################################################################
########################################################################################################
### DATA EXTRACTION ######################################################  ############################
###################### PRECISION-BASED TRAVERSAL #######################################################


def traversal_with_precision(
    overlap_manager: OverlapManager,
    node: str,
    m_stats_stats_matrix,
    input_taxa: pd.DataFrame,
    tax_level: str = "order",
    results=None,
    force_stop=False,
) -> list[pd.DataFrame]:
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
    if results is None:
        results = []

    composition = (
        node_composition_level(overlap_manager, node, m_stats_stats_matrix, input_taxa, tax_level=tax_level)
        .set_index("tax_level")
        .T
    )
    composition.reset_index(drop=True, inplace=True)
    node_true_leaves = node_total_true_leaves(overlap_manager, node, m_stats_stats_matrix)
    node_precision = 1 / len(set(node_true_leaves)) if len(node_true_leaves) > 0 else 0.0
    node_precision = 1 / node_precision if node_precision > 1 else node_precision
    node_row = overlap_manager.all_node_stats[overlap_manager.all_node_stats["Node"] == node]
    min_dist = node_row["Min_Pairwise_Dist"].values[0]
    min_shared = node_row["Min_Shared"].values[0]
    node_total_leaf_taxa_div = node_leaf_shannon_tax_diversity(
        overlap_manager, node, m_stats_stats_matrix, tax_level=tax_level
    )

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

    local_results = pd.DataFrame(
        {
            "node": [node],
            "n_leaves": [len(overlap_manager.get_node_leaves(node))],
            "tax_diversity": [node_total_leaf_taxa_div],
            "n_true_leaves": [len(set(node_true_leaves))],
            "precision": [node_precision],
            "new_precision": [new_precision],
            "precision_increased": [precision_increased],
            "Min_Dist": [min_dist],
            "Min_Shared": [min_shared],
            "stop_traversal": [stop_traversal],
        }
    )
    local_results = pd.concat([local_results, composition], axis=1)
    results.append(local_results)
    for child in node_children:
        if overlap_manager.tree.out_degree(child) > 0:  # internal node
            traversal_with_precision(
                overlap_manager, child, m_stats_stats_matrix, input_taxa, tax_level=tax_level, results=results
            )
    return results


def data_set_traversal_with_precision(
    data_set_name,
    study_output_filepath,
    ncbi_wrapper,
    overlap_manager: OverlapManager,
    input_taxa: pd.DataFrame,
    tax_level: str = "order",
):
    m_stats_stats_matrix = get_m_stats_matrix(
        data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager, filter_no_leaf=True
    )
    results = []

    for root in overlap_manager.root_nodes:
        root_results = traversal_with_precision(
            overlap_manager, root, m_stats_stats_matrix, input_taxa, tax_level=tax_level, results=[]
        )
        root_results = [r for r in root_results if r.empty == False]
        results.extend(root_results)

    if len(results) == 0:
        return pd.DataFrame()

    results_df = pd.concat(results, ignore_index=True)

    results_df.insert(0, "data_set", data_set_name)
    return results_df


#####################################################################################
################ CROSS-HIT ANALYSIS ######################################################
######################################################################################


def cross_hit_prediction_matrix(
    data_set_name,
    study_output_filepath,
    ncbi_wrapper,
    overlap_manager: OverlapManager,
    tax_df,
    tax_level: str = "order",
):
    """
    For each leaf, compute the number of different taxids predicted at tax_level.
    """
    prediction_matrix = get_m_stats_matrix(data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager)
    prediction_matrix.reset_index(inplace=True)

    if prediction_matrix.empty:
        return pd.DataFrame()

    if len(overlap_manager.leaves) == 0:
        return pd.DataFrame()

    prediction_matrix = normalize_by_taxlevel(prediction_matrix, tax_level=tax_level)

    if prediction_matrix.empty:
        return pd.DataFrame()

    composition = get_composition_by_leaf(overlap_manager, prediction_matrix, tax_df, tax_level=tax_level)
    # is_trash is already added by get_m_stats_matrix() and preserved by normalize_by_taxlevel()
    prediction_matrix_stats = prediction_matrix[
        ["leaf", "is_trash", "coverage", "covbases", "meanmapq", "error_rate", "max_shared", "total_uniq_reads"]
    ].copy()

    prediction_matrix_stats = prediction_matrix_stats.merge(composition, on="leaf", how="left")
    prediction_matrix_stats = prediction_matrix_stats[prediction_matrix_stats["leaf"].notna()]

    return prediction_matrix_stats


######################################################################################
################ RECALL CUTOFF #######################################################################
######################################################################################


def predict_recall_cutoff_vars(
    data_set_divide: int,
    data_set_name: str,
    m_stats_stats_matrix: pd.DataFrame,
    taxa_df: pd.DataFrame,
    tax_level: str = "order",
) -> pd.DataFrame:
    """
    Predict recall at various cutoffs and other composition statistics.

    .. deprecated::
       Use ``RecallFeatureTransformer`` instead. This function is kept for
       backward compatibility with ``deployment/analysis/`` scripts.
       The new ``RecallModeller.fit()`` / ``predict_cutoff()`` methods
       handle feature extraction internally.
    """
    import warnings

    warnings.warn(
        "predict_recall_cutoff_vars is deprecated. "
         "Use metagenomics_utils.overlap_manager.feature_transformer.RecallFeatureTransformer instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    # index of last True valuie in best_match_is_best
    m_stats_stats_matrix = m_stats_stats_matrix.sort_values(by="total_uniq_reads", ascending=False).reset_index(
        drop=True
    )
    best_match_indices = m_stats_stats_matrix.index[m_stats_stats_matrix["best_match_is_best"] == True].tolist()
    last_best_match_index = best_match_indices[-1] + 1 if best_match_indices else 0
    last_best_match_relindex = (
        (last_best_match_index) / len(m_stats_stats_matrix) if len(m_stats_stats_matrix) > 0 else 0.0
    )

    composition = node_composition_with_stats(m_stats_stats_matrix, taxa_df, tax_level=tax_level)

    for set_divide in reversed(range(1, data_set_divide + 1)):
        threshold_index = int(m_stats_stats_matrix.shape[0] * set_divide / data_set_divide)
        best_match_indices_short = [idx for idx in best_match_indices if idx < threshold_index]
        index_recall = len(best_match_indices_short) / len(best_match_indices) if len(best_match_indices) > 0 else 0.0
        composition.insert(0, f"index_recall_{set_divide}", index_recall)

    composition.insert(0, "last_best_match_relindex", last_best_match_relindex)
    composition.insert(0, "data_set", data_set_name)

    composition.reset_index(drop=True, inplace=True)

    # composition.drop(columns=['tax_level'], inplace=True)
    return composition


EPS = 1e-6


def compute_ground_truth_division(recall_at_divs, tau, n_divisions, fractions):
    """Find first fraction where recall >= tau."""
    for i in range(n_divisions):
        if recall_at_divs[i] >= tau:
            if i == 0:
                return float(tau * fractions[0] / max(recall_at_divs[0], EPS))
            r_low, r_high = recall_at_divs[i - 1], recall_at_divs[i]
            f_low, f_high = fractions[i - 1], fractions[i]
            if r_high == r_low:
                return f_high
            frac = f_low + (tau - r_low) * (f_high - f_low) / (r_high - r_low)
            return float(np.clip(frac, f_low, f_high))
    return 1.0


def asymmetric_loss(y_true, y_pred, over_cost=1.0, under_cost=3.0):
    """Asymmetric squared loss: overestimation penalised more."""
    diff = y_pred - y_true
    loss = np.where(diff > 0, over_cost * diff**2, under_cost * diff**2)
    return float(np.mean(loss))


def predict_division_clf(probs, X_cutoff, fractions):
    """First division where prob >= X_cutoff; interpolate."""
    if probs.ndim == 1:
        probs = probs.reshape(-1, 1)
    n = probs.shape[1]
    predictions = np.full(n, np.nan)
    for s in range(n):
        idxs = np.where(probs[:, s] >= X_cutoff)[0]
        if len(idxs) == 0:
            predictions[s] = 1.0
            continue
        i = idxs[0]
        if i == 0:
            if probs[0, s] == X_cutoff:
                predictions[s] = fractions[0]
            else:
                predictions[s] = float(X_cutoff * fractions[0] / max(probs[0, s], EPS))
        else:
            p_low, p_high = probs[i - 1, s], probs[i, s]
            f_low, f_high = fractions[i - 1], fractions[i]
            if p_high == p_low:
                predictions[s] = f_high
            else:
                frac = f_low + (X_cutoff - p_low) * (f_high - f_low) / (p_high - p_low)
                predictions[s] = float(np.clip(frac, f_low, f_high))
    return np.clip(predictions, 0.0, 1.0)


def actual_recall_at_fraction(recall_curve, fraction, n_divisions, fractions):
    """Interpolate recall at an arbitrary fraction from a recall curve."""
    if fraction <= 0.0:
        return 0.0
    if fraction >= 1.0:
        return float(recall_curve[-1])
    div_idx = int(np.floor(fraction * n_divisions))
    div_idx = min(div_idx, n_divisions - 2)
    f_low, f_high = fractions[div_idx], fractions[div_idx + 1]
    r_low, r_high = recall_curve[div_idx], recall_curve[div_idx + 1]
    if f_high == f_low:
        return float(r_low)
    frac = (fraction - f_low) / (f_high - f_low)
    return float(r_low + frac * (r_high - r_low))


########################################################################################################
########################################################################################################
######################################################################################
### PREDICTION ######################################################  ##########
###################### RECALL #####################################


def multioutput_regressor(X_train, Y_train):

    rf = RandomForestRegressor(n_estimators=100, random_state=42)
    multi_rf = MultiOutputRegressor(rf)
    multi_rf.fit(X_train, Y_train)
    return multi_rf, None


def xgb_multioutput(X_train, Y_train):
    from xgboost import XGBRegressor

    xgb = XGBRegressor(n_estimators=100, random_state=42)
    multi = MultiOutputRegressor(xgb)
    multi.fit(X_train, Y_train)
    return multi, None


def moxgb_bayes_optimized(X_train, Y_train):
    import optuna
    from sklearn.model_selection import KFold, cross_val_score
    from sklearn.multioutput import MultiOutputRegressor
    from xgboost import XGBRegressor

    cv = KFold(n_splits=10, shuffle=True, random_state=42)

    def objective_rf(trial):

        n_estimators = trial.suggest_int("n_estimators", 50, 300)
        max_depth = trial.suggest_int("max_depth", 3, 20)
        learning_rate = trial.suggest_float("learning_rate", 0.01, 0.3)
        reg_alpha = trial.suggest_float("reg_alpha", 0.0, 1.0)
        reg_lambda = trial.suggest_float("reg_lambda", 0.0, 1.0)

        model = MultiOutputRegressor(
            XGBRegressor(
                learning_rate=learning_rate,
                reg_alpha=reg_alpha,
                reg_lambda=reg_lambda,
                n_estimators=n_estimators,
                max_depth=max_depth,
                random_state=42,
            )
        )

        scores = cross_val_score(model, X_train, Y_train, cv=cv, scoring="neg_mean_absolute_error").mean()
        return scores

    study_rf = optuna.create_study(direction="minimize")
    study_rf.optimize(objective_rf, n_trials=50)
    rf_optuna = MultiOutputRegressor(XGBRegressor(**study_rf.best_params))
    rf_optuna.fit(X_train, Y_train)

    return rf_optuna, study_rf


def random_forest_reg_multioutput_bayes_optimized(X_train, Y_train):
    import optuna
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import KFold, cross_val_score
    from sklearn.multioutput import MultiOutputRegressor

    cv = KFold(n_splits=10, shuffle=True, random_state=42)

    def objective_rf(trial):
        n_estimators = trial.suggest_int("n_estimators", 50, 300)
        max_depth = trial.suggest_int("max_depth", 3, 20)
        min_samples_split = trial.suggest_int("min_samples_split", 2, 20)
        min_samples_leaf = trial.suggest_int("min_samples_leaf", 1, 20)

        model = MultiOutputRegressor(
            RandomForestRegressor(
                n_estimators=n_estimators,
                max_depth=max_depth,
                min_samples_split=min_samples_split,
                min_samples_leaf=min_samples_leaf,
                random_state=42,
            )
        )

        scores = cross_val_score(model, X_train, Y_train, cv=cv, scoring="neg_mean_absolute_error").mean()
        return scores

    study_rf = optuna.create_study(direction="minimize")
    study_rf.optimize(objective_rf, n_trials=50)
    rf_optuna = MultiOutputRegressor(RandomForestRegressor(**study_rf.best_params))
    rf_optuna.fit(X_train, Y_train)

    best_params = study_rf.best_params
    best_model = MultiOutputRegressor(RandomForestRegressor(**best_params))
    best_model.fit(X_train, Y_train)
    return best_model, study_rf


def neural_network_multioutput_bayes_optimized(X_train, Y_train):
    import optuna
    from sklearn.model_selection import KFold, cross_val_score
    from sklearn.multioutput import MultiOutputRegressor
    from sklearn.neural_network import MLPRegressor

    cv = KFold(n_splits=10, shuffle=True, random_state=42)

    def objective_nn(trial):
        hidden_layer_sizes = trial.suggest_categorical(
            "hidden_layer_sizes", [(20, 50, 30), (20, 50, 50), (100, 50, 50), (40, 60, 20)]
        )
        activation = trial.suggest_categorical("activation", ["relu", "tanh"])
        alpha = trial.suggest_float("alpha", 1e-5, 1e-2, log=True)

        model = MultiOutputRegressor(
            MLPRegressor(
                hidden_layer_sizes=hidden_layer_sizes, activation=activation, alpha=alpha, max_iter=500, random_state=42
            )
        )

        scores = cross_val_score(model, X_train, Y_train, cv=cv, scoring="neg_mean_absolute_error").mean()
        return scores

    study_nn = optuna.create_study(direction="minimize")
    study_nn.optimize(objective_nn, n_trials=50)
    nn_optuna = MultiOutputRegressor(MLPRegressor(**study_nn.best_params))
    nn_optuna.fit(X_train, Y_train)

    best_params = study_nn.best_params
    best_model = MultiOutputRegressor(MLPRegressor(**best_params))
    best_model.fit(X_train, Y_train)
    return best_model, study_nn


class GPModelWrapper:
    """Container for per-division GP regressors with feature scaling."""

    def __init__(self):
        self.scaler = StandardScaler()
        self.models = []  # list of dicts with 'type': 'gp'|'constant', 'model', etc.
        self.feature_names_ = []
        self.target_names_ = []
        self.recall_col_indices_ = []
        self.meta_col_indices_ = []

    def fit(self, X, Y):
        """Train per-column GP regressors. X, Y are DataFrames."""
        self.feature_names_ = list(X.columns)
        self.target_names_ = list(Y.columns)

        self.recall_col_indices_ = [i for i, n in enumerate(self.target_names_) if n.startswith("index_recall_")]
        self.meta_col_indices_ = [i for i, n in enumerate(self.target_names_) if not n.startswith("index_recall_")]

        X_scaled = self.scaler.fit_transform(X.values)
        Y_array = Y.values
        n_targets = Y_array.shape[1]
        kernel_template = ConstantKernel(1.0) * RBF(length_scale=1.0) + WhiteKernel(noise_level=0.1)

        for i in range(n_targets):
            y_col = Y_array[:, i]
            if i in self.meta_col_indices_:
                self.models.append(
                    {
                        "type": "constant",
                        "mean": float(np.mean(y_col)),
                        "std": float(max(np.std(y_col), 0.01)),
                    }
                )
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    gp = GaussianProcessRegressor(
                        kernel=deepcopy(kernel_template), normalize_y=True, n_restarts_optimizer=3, random_state=42
                    )
                    gp.fit(X_scaled, y_col)
                    self.models.append({"type": "gp", "model": gp})

    def predict(self, X):
        """Return posterior means — same shape as MultiOutputRegressor."""
        means, _ = self._predict_all(X)
        return means

    def predict_with_std(self, X):
        """Return (means, stds) arrays."""
        return self._predict_all(X)

    def _predict_all(self, X):
        X_scaled = self.scaler.transform(X.values if hasattr(X, "values") else np.asarray(X))
        n = X_scaled.shape[0]
        n_targets = len(self.models)
        means = np.zeros((n, n_targets))
        stds = np.ones((n, n_targets)) * 0.01

        for i, info in enumerate(self.models):
            if info["type"] == "constant":
                means[:, i] = info["mean"]
                stds[:, i] = info["std"]
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    mu, sigma = info["model"].predict(X_scaled, return_std=True)
                    means[:, i] = mu.flatten()
                    stds[:, i] = np.maximum(sigma.flatten(), 1e-6)
        return means, stds


def gp_clf_training_function(X_train, Y_train):
    """Train per-division GPs; returns (GPModelWrapper, None) for InjectModellerInterface."""
    wrapper = GPModelWrapper()
    wrapper.fit(X_train, Y_train)
    return wrapper, None


class InjectModellerInterface:
    model_save_filename = "inject_xgb_bundle.pkl"

    model_map = {
        "xgb": xgb_multioutput,
        "morf": multioutput_regressor,
        "moxgb_optimized": moxgb_bayes_optimized,
        "morf_optimized": random_forest_reg_multioutput_bayes_optimized,
        "monn_optimized": neural_network_multioutput_bayes_optimized,
        "gp_clf": gp_clf_training_function,
    }

    def __init__(self, model_type: str = "xgb"):

        self.training_func = self.model_map.get(model_type, multioutput_regressor)
        self.model_type = model_type

        if self.training_func is None:
            print(f"Model {model_type} not found. Defaulting to multioutput_regressor.")
            self.training_func = multioutput_regressor

    def train_model(self, X_train, Y_train, bayes_optimize: bool = True) -> MultiOutputRegressor:

        model, _study = self.training_func(X_train, Y_train)
        return model


class RecallModeller:
    model_save_filename = "recall_xgb_bundle.pkl"
    model_category = "recall"

    def __getstate__(self):
        state = self.__dict__.copy()
        state.pop("X_test", None)
        state.pop("y_test", None)
        state.pop("model_interface", None)
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.X_test = None
        self.y_test = None
        self.model_interface = None

    def __init__(
        self,
        data_set_divide: int = 5,
        model_interface: InjectModellerInterface | None = None,
        tax_level: str = "order",
        feature_transformer=None,
        sort_strategy: str = "reads",
        description: str | None = None,
    ):
        self.data_set_divide = data_set_divide
        self.tax_level = tax_level
        self.model: MultiOutputRegressor | None = None
        self.X_test = None
        self.y_test = None
        self.date_trained: str | None = None
        self.model_interface = model_interface or InjectModellerInterface(model_type="xgb")
        model_type = self.model_interface.model_type if self.model_interface else "xgb"
        self.description = description or f"Multi-output recall model ({model_type})"
        from metagenomics_utils.overlap_manager.feature_transformer import RecallFeatureTransformer

        self.transformer = feature_transformer or RecallFeatureTransformer(
            tax_level=tax_level,
            data_set_divide=data_set_divide,
            sort_strategy=sort_strategy,
        )

    def fit(
        self,
        m_stats_matrices,
        taxids_to_use,
        test_size=0.2,
        random_state=42,
    ):
        self.transformer.fit(m_stats_matrices, taxids_to_use=taxids_to_use)
        rows = [self.transformer.transform(m) for m in m_stats_matrices]
        full = pd.concat(rows, ignore_index=True)

        feat_cols = self.transformer.get_feature_names_out()
        target_cols = self.transformer.target_columns_

        X = full[feat_cols]
        Y = full[target_cols]

        logger = logging.getLogger(__name__)
        logger.info(
            "RecallModeller.fit: %d samples, %d features, %d targets",
            len(full),
            len(feat_cols),
            len(target_cols),
        )

        from sklearn.model_selection import train_test_split

        X_train, X_test, Y_train, Y_test = train_test_split(
            X,
            Y,
            test_size=test_size,
            random_state=random_state,
        )
        self.X_test = X_test
        self.y_test = Y_test

        self.model = self.model_interface.train_model(X_train, Y_train)

        self.RecP_feature_cols = feat_cols
        self.RecP_target_cols = target_cols
        self.date_trained = datetime.now(timezone.utc).isoformat()

        return self

    # ── Persistence ──

    def save_model(self, output_directory: str):
        if self.model is None:
            print("No model to save.")
            return
        bundle = {
            "model_type": "xgb_multi",
            "model_category": self.model_category,
            "description": self.description,
            "date_trained": self.date_trained,
            "model": self.model,
            "feature_names": self.RecP_feature_cols,
            "target_names": self.RecP_target_cols,
            "data_set_divide": self.data_set_divide,
            "tax_level": self.tax_level,
            "transformer": self.transformer,
        }
        if hasattr(self, "target_recall"):
            bundle["target_recall"] = self.target_recall
        joblib.dump(bundle, os.path.join(output_directory, self.model_save_filename))

    def load_model(self, input_directory: str):
        try:
            bundle = joblib.load(os.path.join(input_directory, self.model_save_filename))
            if isinstance(bundle, dict):
                self.model = bundle["model"]
                self.RecP_feature_cols = bundle["feature_names"]
                self.RecP_target_cols = bundle["target_names"]
                self.data_set_divide = bundle.get("data_set_divide", self.data_set_divide)
                self.tax_level = bundle.get("tax_level", self.tax_level)
                self.description = bundle.get("description", self.description)
                self.date_trained = bundle.get("date_trained", self.date_trained)
                self.transformer = bundle["transformer"]
                if "target_recall" in bundle:
                    self.target_recall = bundle["target_recall"]
            else:
                self.model = bundle
        except Exception as e:
            print(f"Error loading model: {e}")

    # ── Evaluation ──

    def evaluate_model(self, model, X_test, Y_test, ouptput_filepath):
        from sklearn.metrics import mean_squared_error, r2_score

        Y_pred = model.predict(X_test)
        r2_scores = {}
        mse_scores = {}
        with open(ouptput_filepath, "w") as f:
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

        plt.boxplot(differences, labels=self.RecP_target_cols)
        plt.title("Distribution of Differences between True and Predicted Recall")
        plt.ylabel("Difference")
        plt.savefig(analysis_output_filepath)
        plt.close()

    def model_summary(self, model, X_test, Y_test, analysis_output_filedir):
        Y_pred = model.predict(X_test)
        from sklearn.metrics import mean_squared_error, r2_score

        r2 = r2_score(Y_test, Y_pred, multioutput="uniform_average")
        mse = mean_squared_error(Y_test, Y_pred, multioutput="uniform_average")

        print(f"model_summary R² = {r2:.3f}, MSE = {mse:.3f}")

        analysis_output_filepath = os.path.join(analysis_output_filedir, "recall_model_analysis_results.txt")
        r2_scores, mse_scores = self.evaluate_model(model, X_test, Y_test, analysis_output_filepath)
        feat_importance_filepath = analysis_output_filepath.replace(".txt", "_feature_importances.tsv")
        self.feature_importances(model, feat_importance_filepath)
        self.plot_eval(X_test, Y_test, analysis_output_filepath.replace(".txt", "_recall_prediction_differences.png"))
        import matplotlib.pyplot as plt

        print(Y_test.shape)

        for i in range(3):
            if i > Y_test.shape[0] - 1:
                continue
            plt.plot(range(1, Y_test.shape[1]), Y_test.iloc[i, 1:], "o-", label="True")
            plt.plot(range(1, Y_test.shape[1]), Y_pred[i, 1:], "s--", label="Predicted")
            plt.title(f"Sample {i}: recall curve")
            plt.xlabel("Recall index (1–5)")
            plt.ylabel("Recall value")
            plt.legend()
            plt.savefig(analysis_output_filepath.replace(".txt", f"_recall_curve_sample_{i}.png"))
            plt.close()

        return r2_scores, mse_scores

    def predict_cutoff(self, X, target_recall, confidence=None):
        X_feat = self.transformer.transform(X)
        X_feat = X_feat[self.transformer.get_feature_names_out()]
        recall_pred = self.model.predict(X_feat)
        recall_pred_df = pd.DataFrame(recall_pred, columns=self.RecP_target_cols)
        return find_percentile_for_recall(recall_pred_df.iloc[0], self.data_set_divide, target_recall)

    def predict(self, X):
        X_feat = self.transformer.transform(X)
        X_feat = X_feat[self.transformer.get_feature_names_out()]
        return self.model.predict(X_feat)

    def evaluate_prediction(self, X_raw, target_recall=1.0):
        X_full = self.transformer.transform(X_raw)
        feat_cols = self.transformer.get_feature_names_out()
        target_cols = self.transformer.target_columns_

        truth = X_full[target_cols].iloc[0]
        X_feat = X_full[feat_cols]
        pred_raw = self.model.predict(X_feat)

        metrics = {
            "last_best_match_relindex": float(truth["last_best_match_relindex"]),
        }

        recall_cols = [c for c in target_cols if c.startswith("index_recall_")]

        if pred_raw.ndim == 2 and pred_raw.shape[1] == len(target_cols):
            pred = pd.Series(pred_raw[0], index=target_cols)
            errors = {c: float(pred[c] - truth[c]) for c in recall_cols}
            metrics["per_division_recall_errors"] = errors
            metrics["per_division_recall_rmse"] = float(np.sqrt(np.mean([e**2 for e in errors.values()])))
        else:
            pred_k_min = int(pred_raw[0])
            actual_k_min = next(
                (i + 1 for i, c in enumerate(recall_cols) if truth[c] >= target_recall), self.data_set_divide
            )
            metrics["predicted_k_min"] = pred_k_min
            metrics["actual_k_min"] = actual_k_min
            metrics["cutoff_error"] = pred_k_min - actual_k_min

        return metrics


class CutoffRecallModeller(RecallModeller):
    """
    Predict minimum bins k_min needed to achieve target_recall.

    Trains a RandomForest classifier: features -> k_min (1..data_set_divide).
    Supports predict_proba for probability-guided cutoff decisions.
    """

    model_save_filename = "cutoff_recall_bundle.pkl"

    def __init__(
        self,
        data_set_divide: int = 5,
        target_recall: float = 1.0,
        tax_level: str = "order",
        feature_transformer=None,
        sort_strategy: str = "reads",
        description: str | None = None,
    ):
        super().__init__(
            data_set_divide=data_set_divide,
            tax_level=tax_level,
            feature_transformer=feature_transformer,
            sort_strategy=sort_strategy,
            description=description or "Cutoff recall model (RandomForest)",
        )
        self.target_recall = target_recall
        self.model: RandomForestClassifier | None = None

    def _compute_k_min(self, row):
        for k in range(1, self.data_set_divide + 1):
            if row.get(f"index_recall_{k}", 0.0) >= self.target_recall:
                return k
        return self.data_set_divide

    def fit(
        self,
        m_stats_matrices,
        taxids_to_use,
        test_size=0.2,
        random_state=42,
    ):
        self.transformer.fit(m_stats_matrices, taxids_to_use=taxids_to_use)
        rows = [self.transformer.transform(m) for m in m_stats_matrices]
        full = pd.concat(rows, ignore_index=True)

        feat_cols = self.transformer.get_feature_names_out()
        target_cols = self.transformer.target_columns_

        X = full[feat_cols]
        recall_cols = [c for c in target_cols if c.startswith("index_recall_")]
        y = full[recall_cols].apply(self._compute_k_min, axis=1)

        logger = logging.getLogger(__name__)
        logger.info(
            "CutoffRecallModeller.fit: %d samples, %d features",
            len(full),
            len(feat_cols),
        )

        X_train, X_test, y_train, y_test = train_test_split(
            X,
            y,
            test_size=test_size,
            random_state=random_state,
        )
        self.X_test = X_test
        self.y_test = y_test

        clf = RandomForestClassifier(
            n_estimators=300,
            max_depth=12,
            min_samples_leaf=3,
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        )
        clf.fit(X_train, y_train)
        self.model = clf
        self.RecP_feature_cols = feat_cols
        self.RecP_target_cols = target_cols
        self.date_trained = datetime.now(timezone.utc).isoformat()
        return self

    def predict_cutoff(self, X, target_recall=None, confidence=None):
        X_feat = self.transformer.transform(X)
        X_feat = X_feat[self.transformer.get_feature_names_out()]
        X_arr = X_feat.values

        if target_recall is not None and abs(target_recall - self.target_recall) > 1e-6:
            raise ValueError(
                f"CutoffRecallModeller trained for target_recall={self.target_recall}, "
                f"but predict_cutoff called with target_recall={target_recall}. "
                "Train a separate model per target recall or use --target_recall "
                f"matching the training value ({self.target_recall})."
            )
        if confidence is not None:
            proba = self.model.predict_proba(X_arr)[0]
            classes = list(self.model.classes_)
            full = np.zeros(self.data_set_divide)
            for c, p in zip(classes, proba):
                full[int(c) - 1] = p
            cum = np.cumsum(full)
            k = np.argmax(cum >= confidence) + 1 if np.any(cum >= confidence) else self.data_set_divide
        else:
            k = int(self.model.predict(X_arr)[0])
        return k / self.data_set_divide

    def save_model(self, output_directory: str):
        if self.model is None:
            print("No model to save.")
            return
        bundle = {
            "model_type": "xgb_direct",
            "model_category": self.model_category,
            "description": self.description,
            "date_trained": self.date_trained,
            "model": self.model,
            "feature_names": self.RecP_feature_cols,
            "target_recall": self.target_recall,
            "data_set_divide": self.data_set_divide,
            "tax_level": self.tax_level,
            "transformer": self.transformer,
        }
        joblib.dump(bundle, os.path.join(output_directory, self.model_save_filename))

    def load_model(self, input_directory: str):
        try:
            bundle = joblib.load(os.path.join(input_directory, self.model_save_filename))
            if isinstance(bundle, dict):
                self.model = bundle["model"]
                self.RecP_feature_cols = bundle["feature_names"]
                self.target_recall = bundle["target_recall"]
                self.data_set_divide = bundle["data_set_divide"]
                self.tax_level = bundle.get("tax_level", self.tax_level)
                self.description = bundle.get("description", self.description)
                self.date_trained = bundle.get("date_trained", self.date_trained)
                self.transformer = bundle["transformer"]
            else:
                self.model = bundle
        except Exception as e:
            print(f"Error loading model: {e}")

    def evaluate_model(self, model, X_test, Y_test, output_filepath):
        from sklearn.metrics import accuracy_score

        Y_pred = model.predict(X_test)
        acc = accuracy_score(Y_test, Y_pred)
        with open(output_filepath, "w") as f:
            f.write(f"accuracy\t{acc}\n")
        return {"accuracy": acc}, {}

    def model_summary(self, model, X_test, Y_test, analysis_output_filedir):
        from sklearn.metrics import accuracy_score

        Y_pred = model.predict(X_test)
        acc = accuracy_score(Y_test, Y_pred)
        print(f"CutoffRecallModeller test accuracy = {acc:.3f}")
        output_filepath = os.path.join(analysis_output_filedir, "recall_model_analysis_results.txt")
        self.evaluate_model(model, X_test, Y_test, output_filepath)
        return {"accuracy": acc}, {}


class DirectXGBRecallModeller(RecallModeller):
    """
    Direct XGBoost regression for recall cutoff fraction prediction.

    Trains an XGBRegressor to directly predict the tau-crossing fraction
    (not bin indices). Uses asymmetric sample weights to penalise
    underestimation more heavily (3:1 under:over).

    Inherits prep_data, split_data, predict from RecallModeller.
    """

    model_save_filename = "direct_xgb_bundle.pkl"

    def __init__(
        self,
        data_set_divide: int = 5,
        target_recall: float = 1.0,
        tax_level: str = "order",
        feature_transformer=None,
        sort_strategy: str = "reads",
        description: str | None = None,
    ):
        super().__init__(
            data_set_divide=data_set_divide,
            tax_level=tax_level,
            feature_transformer=feature_transformer,
            sort_strategy=sort_strategy,
            description=description or "Direct XGBoost recall model",
        )
        self.target_recall = target_recall
        self.model: RegressorMixin | None = None
        self.fractions = np.arange(1, data_set_divide + 1) / data_set_divide

    def _compute_ground_truth(self, row):
        recall_vals = np.array([row.get(f"index_recall_{d}", 0.0) for d in range(1, self.data_set_divide + 1)])
        return compute_ground_truth_division(
            recall_vals,
            self.target_recall,
            self.data_set_divide,
            self.fractions,
        )

    def _asymmetric_weights(self, gt, over_cost=1.0, under_cost=3.0, floor=0.05):
        base = 1.0 / np.clip(gt, floor, 1.0)
        return base * under_cost + (1.0 - base / base.max()) * over_cost

    def fit(
        self,
        m_stats_matrices,
        taxids_to_use,
        test_size=0.2,
        random_state=42,
    ):
        from xgboost import XGBRegressor

        self.transformer.fit(m_stats_matrices, taxids_to_use=taxids_to_use)
        rows = [self.transformer.transform(m) for m in m_stats_matrices]
        full = pd.concat(rows, ignore_index=True)

        feat_cols = self.transformer.get_feature_names_out()
        target_cols = self.transformer.target_columns_

        X = full[feat_cols]
        recall_cols = [c for c in target_cols if c.startswith("index_recall_")]
        y = full[recall_cols].apply(self._compute_ground_truth, axis=1)

        logger = logging.getLogger(__name__)
        logger.info(
            "DirectXGBRecallModeller.fit: %d samples, %d features",
            len(full),
            len(feat_cols),
        )

        X_train, X_test, y_train, y_test = train_test_split(
            X,
            y,
            test_size=test_size,
            random_state=random_state,
        )
        self.X_test = X_test
        self.y_test = y_test

        sw = self._asymmetric_weights(y_train.values, over_cost=1.0, under_cost=3.0)
        reg = XGBRegressor(n_estimators=200, random_state=42, verbosity=0)
        reg.fit(X_train, y_train, sample_weight=sw)
        self.model = reg
        self.RecP_feature_cols = feat_cols
        self.RecP_target_cols = target_cols
        self.date_trained = datetime.now(timezone.utc).isoformat()
        return self

    def predict_cutoff(self, X, target_recall=None, confidence=None):
        X_feat = self.transformer.transform(X)
        X_feat = X_feat[self.transformer.get_feature_names_out()]
        return float(np.clip(self.model.predict(X_feat)[0], 0.0, 1.0))

    def save_model(self, output_directory: str):
        if self.model is None:
            print("No model to save.")
            return
        bundle = {
            "model_type": "direct_xgb",
            "model_category": self.model_category,
            "description": self.description,
            "date_trained": self.date_trained,
            "model": self.model,
            "feature_names": self.RecP_feature_cols,
            "target_recall": self.target_recall,
            "data_set_divide": self.data_set_divide,
            "tax_level": self.tax_level,
            "transformer": self.transformer,
        }
        joblib.dump(bundle, os.path.join(output_directory, self.model_save_filename))

    def load_model(self, input_directory: str):
        try:
            bundle = joblib.load(os.path.join(input_directory, self.model_save_filename))
            if isinstance(bundle, dict):
                self.model = bundle["model"]
                self.RecP_feature_cols = bundle["feature_names"]
                self.target_recall = bundle["target_recall"]
                self.data_set_divide = bundle["data_set_divide"]
                self.tax_level = bundle.get("tax_level", self.tax_level)
                self.description = bundle.get("description", self.description)
                self.date_trained = bundle.get("date_trained", self.date_trained)
                self.fractions = np.arange(1, self.data_set_divide + 1) / self.data_set_divide
                self.transformer = bundle["transformer"]
            else:
                self.model = bundle
        except Exception as e:
            print(f"Error loading model: {e}")

    def evaluate_model(self, model, X_test, Y_test, output_filepath):
        from sklearn.metrics import mean_squared_error, r2_score

        Y_pred = model.predict(X_test)
        r2 = r2_score(Y_test, Y_pred)
        mse = mean_squared_error(Y_test, Y_pred)
        with open(output_filepath, "w") as f:
            f.write(f"r2\t{r2}\nmse\t{mse}\n")
        return {"r2": r2, "mse": mse}, {}

    def model_summary(self, model, X_test, Y_test, analysis_output_filedir):
        from sklearn.metrics import mean_squared_error, r2_score

        Y_pred = model.predict(X_test)
        r2 = r2_score(Y_test, Y_pred)
        mse = mean_squared_error(Y_test, Y_pred)
        print(f"DirectXGB model_summary R² = {r2:.3f}, MSE = {mse:.3f}")
        output_filepath = os.path.join(analysis_output_filedir, "recall_model_analysis_results.txt")
        self.evaluate_model(model, X_test, Y_test, output_filepath)
        return {"r2": r2, "mse": mse}, {}


class GPCLFThreshold(RegressorMixin, BaseEstimator):
    """
    Gaussian Process ensemble + CLF threshold for recall cutoff prediction.

    Trains per-division GPs to predict recall at each division,
    then applies a CLF threshold to determine the target_percentile.

    Improves GP fit quality by:
    - Normalising recall curves per-sample to [0, 1] (each sample's max recall = 1)
    - Filtering degenerate samples (max recall = 0 or already >= τ at division 1)
    - Optionally binarising taxonomy features to 0/1
    - Optimising τ/X on a held-out validation split to avoid overfitting

    Parameters
    ----------
    tau : float, default=0.9
        Target recall threshold (τ), in normalised space (fraction of per-sample max).
    x_thresh : float, default=0.5
        Confidence threshold (X) for the CLF decision rule.
    n_divisions : int, default=5
        Number of recall divisions.
    optimize : bool, default=True
        Whether to grid-search τ and X during fit to minimise asymmetric loss.
    filter_degenerate : bool, default=True
        Whether to remove samples with recall max = 0 or recall already >= τ
        at the first division. These are edge cases the GP cannot model well.
    binarize_taxonomy : bool, default=True
        Whether to binarise taxonomy features (all features not in the known
        stats set) to 0/1 before scaling.
    val_split : float, default=0.2
        Fraction of training data held out for threshold optimisation.
        If 0.0, optimisation uses all training data (in-sample).
    taxonomy_cols : list of str or None, default=None
        Explicit list of taxonomy column names for binarisation.
        If None and binarize_taxonomy=True, auto-detect as all feature names
        not in the known stats feature set.
    random_state : int, default=42
        Random state for validation split shuffling.
    """

    def __init__(
        self,
        tau=0.9,
        x_thresh=0.5,
        n_divisions=5,
        optimize=True,
        filter_degenerate=True,
        binarize_taxonomy=True,
        val_split=0.2,
        taxonomy_cols=None,
        random_state=42,
    ):
        self.tau = tau
        self.x_thresh = x_thresh
        self.n_divisions = n_divisions
        self.optimize = optimize
        self.filter_degenerate = filter_degenerate
        self.binarize_taxonomy = binarize_taxonomy
        self.val_split = val_split
        self.taxonomy_cols = taxonomy_cols
        self.random_state = random_state
        self.fractions_ = np.arange(1, n_divisions + 1) / n_divisions
        self._stats_features = {
            "counts_kurtosis",
            "counts_skewness",
            "tax_diversity_shannon",
            "max_uniq_reads",
            "total_uniq_reads",
        }

    def __getstate__(self):
        state = super().__getstate__()
        for attr in ("_stats_features", "taxonomy_cols", "optimize",
                     "filter_degenerate", "binarize_taxonomy", "val_split",
                     "random_state", "feature_names_", "target_names_"):
            state.pop(attr, None)
        return state

    def __setstate__(self, state):
        super().__setstate__(state)
        self._stats_features = {
            "counts_kurtosis", "counts_skewness", "tax_diversity_shannon",
            "max_uniq_reads", "total_uniq_reads",
        }
        self.taxonomy_cols = None
        self.optimize = True
        self.filter_degenerate = True
        self.binarize_taxonomy = True
        self.val_split = 0.2
        self.random_state = 42
        self.feature_names_ = None
        self.target_names_ = None

    def fit(self, X, y, target_names=None):
        """
        Fit per-division GPs and optionally optimize thresholds.

        Applies recall curve normalisation, degenerate sample filtering,
        taxonomy binarisation, and validation-based threshold optimisation.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        y : array-like, shape (n_samples, n_targets)
            Columns: index_recall_1..N, plus optional meta columns.
        target_names : list of str, optional
            Names of target columns. If provided, stored as target_names_.
            Columns NOT starting with 'index_recall_' are treated as
            constant-output (no GP fitted).
        """
        self.feature_names_ = list(X.columns) if hasattr(X, "columns") else None

        # ── Binarise taxonomy features (before numpy conversion) ──
        if self.binarize_taxonomy and self.feature_names_ is not None and hasattr(X, "columns"):
            if self.taxonomy_cols is not None:
                tax_cols = [c for c in self.taxonomy_cols if c in X.columns]
            else:
                tax_cols = [c for c in self.feature_names_ if c not in self._stats_features]
            for col in tax_cols:
                X[col] = (X[col] > 0).astype(float)

        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        n_targets = y.shape[1]

        if target_names is not None:
            self.target_names_ = list(target_names)
            self.recall_indices_ = [i for i, n in enumerate(self.target_names_) if str(n).startswith("index_recall_")]
            self.meta_indices_ = [i for i, n in enumerate(self.target_names_) if not str(n).startswith("index_recall_")]
        else:
            self.recall_indices_ = list(range(n_targets))
            self.meta_indices_ = []

        # ── Filter degenerate samples ──
        if self.filter_degenerate:
            if "last_best_match_relindex" in self.target_names_:
                lbmr_idx = next(i for i in self.meta_indices_ if self.target_names_[i] == "last_best_match_relindex")
                degenerate = (y[:, lbmr_idx] <= 0.0) | (y[:, lbmr_idx] >= 1.0)
            else:
                recall_y = y[:, self.recall_indices_]
                recall_max = recall_y.max(axis=1)
                degenerate = recall_max <= 0.0
            valid = ~degenerate
            if valid.sum() < X.shape[0]:
                X = X[valid]
                y = y[valid]

        # ── Split train/val before GP training ──
        n_samples = y.shape[0]
        if self.optimize and self.val_split > 0:
            rng = np.random.RandomState(self.random_state)
            perm = rng.permutation(n_samples)
            split_idx = int(n_samples * (1.0 - self.val_split))
            train_idx = perm[:split_idx]
            val_idx = perm[split_idx:]
        else:
            train_idx = np.arange(n_samples)
            val_idx = np.arange(n_samples)

        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]

        # ── Scale ──
        self.scaler_ = StandardScaler()
        X_train_scaled = self.scaler_.fit_transform(X_train)
        X_val_scaled = self.scaler_.transform(X_val) if self.val_split > 0 else X_train_scaled
        X_val_scaled[np.isnan(X_val_scaled)] = 0.0

        # ── Normalise recall curves per-sample ──
        recall_y_train = y_train[:, self.recall_indices_]
        sample_max = recall_y_train.max(axis=1, keepdims=True)
        sample_max = np.maximum(sample_max, EPS)
        recall_y_train_norm = recall_y_train / sample_max
        y_train_norm = y_train.copy()
        y_train_norm[:, self.recall_indices_] = recall_y_train_norm

        if self.val_split > 0:
            recall_y_val = y_val[:, self.recall_indices_]
            val_sample_max = recall_y_val.max(axis=1, keepdims=True)
            val_sample_max = np.maximum(val_sample_max, EPS)
            recall_y_val_norm = recall_y_val / val_sample_max
            y_val_norm = y_val.copy()
            y_val_norm[:, self.recall_indices_] = recall_y_val_norm
        else:
            y_val_norm = y_train_norm

        # ── Fit GPs on training split (normalised recall) ──
        n_targets = y_train_norm.shape[1]
        kernel_template = ConstantKernel(1.0) * RBF(length_scale=1.0) + WhiteKernel(noise_level=0.1)
        self.models_ = []

        for i in range(n_targets):
            y_col = y_train_norm[:, i]
            if i in self.meta_indices_:
                self.models_.append(
                    {
                        "type": "constant",
                        "mean": float(np.mean(y_col)),
                        "std": float(max(np.std(y_col), 0.01)),
                    }
                )
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    gp = GaussianProcessRegressor(
                        kernel=deepcopy(kernel_template),
                        normalize_y=True,
                        n_restarts_optimizer=3,
                        random_state=self.random_state,
                    )
                    gp.fit(X_train_scaled, y_col)
                    self.models_.append({"type": "gp", "model": gp})

        # ── Optimise thresholds on validation split ──
        if self.optimize:
            self._optimize_thresholds(X_val_scaled, y_val_norm)

        return self

    def _predict_all(self, X):
        """Return (means, stds) arrays for all targets."""
        means = np.zeros((X.shape[0], len(self.models_)))
        stds = np.ones((X.shape[0], len(self.models_))) * 0.01

        for i, info in enumerate(self.models_):
            if info["type"] == "constant":
                means[:, i] = info["mean"]
                stds[:, i] = info["std"]
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    mu, sigma = info["model"].predict(X, return_std=True)
                    means[:, i] = mu.flatten()
                    stds[:, i] = np.maximum(sigma.flatten(), 1e-6)
        return means, stds

    def predict_raw(self, X):
        """
        Return per-target posterior means for recall columns (NORMALISED scale,
        i.e. each sample's recall curve is a fraction of its per-sample max).
        Shape matches MultiOutputRegressor: (n_samples, n_targets).
        """
        X = np.asarray(X, dtype=float)
        X_scaled = self.scaler_.transform(X)
        
        X_scaled[np.isnan(X_scaled)] = 0.0

        means, _ = self._predict_all(X_scaled)
        return means

    def _optimize_thresholds(self, X_scaled, y):
        """
        Grid search over τ and X to minimise cutoff fraction error.

        For each (τ, X) pair:
          1. Compute P(recall > τ) per division from the GP posterior.
          2. Predict the cutoff fraction where P(recall > τ) >= X.
          3. Compute the true cutoff fraction where recall >= τ.
          4. Loss = asymmetric squared gap between true and predicted
             fractions. Over-prediction (too late) penalised 3× heavier
             than under-prediction (too early), incentivising slightly
             conservative cutoffs that preserve recall.
        """
        tau_values = np.linspace(0.10, 0.975, 36)
        X_values = np.linspace(0.025, 0.975, 39)
        n_val = y.shape[0]
        y_recall = y[:, self.recall_indices_]

        means, stds = self._predict_all(X_scaled)
        recall_means = means[:, self.recall_indices_]
        recall_stds = np.maximum(stds[:, self.recall_indices_], 1e-6)

        best_loss = float("inf")
        best_tau = None
        best_X = None
        
        for tau in tau_values:
            probs = np.zeros((self.n_divisions, n_val))
            for d in range(self.n_divisions):
                mu = recall_means[:, d]
                sigma = recall_stds[:, d]
                probs[d, :] = 1.0 - norm.cdf(tau, loc=mu, scale=sigma)

            y_true_frac = np.array(
                [
                    compute_ground_truth_division(y_recall[s], tau, self.n_divisions, self.fractions_)
                    for s in range(n_val)
                ]
            )

            for xc in X_values:
                y_pred = predict_division_clf(probs, xc, self.fractions_)
                loss = asymmetric_loss(y_true_frac, y_pred, over_cost=3, under_cost=1)
                if loss < best_loss:
                    best_loss = loss
                    best_tau = tau
                    best_X = xc

        self.optimal_tau_ = best_tau
        self.optimal_X_ = best_X
        self.optimal_loss_ = best_loss

    def predict(self, X, tau=None, x_thresh=None):
        """
        Predict target_percentile using CLF threshold on GP posterior.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        tau : float or None
            Override τ. Defaults to optimal_tau_ if optimised, else self.tau.
        x_thresh : float or None
            Override X. Defaults to optimal_X_ if optimised, else self.x_thresh.

        Returns
        -------
        ndarray of shape (n_samples,)
        """
        tau = tau if tau is not None else (getattr(self, "optimal_tau_", None) or self.tau)
        xc = x_thresh if x_thresh is not None else (getattr(self, "optimal_X_", None) or self.x_thresh)

        X = np.asarray(X, dtype=float)
        X_scaled = self.scaler_.transform(X)
        means, stds = self._predict_all(X_scaled)

        recall_means = means[:, self.recall_indices_]
        recall_stds = np.maximum(stds[:, self.recall_indices_], 1e-6)

        n = X.shape[0]
        probs = np.zeros((self.n_divisions, n))
        for d in range(self.n_divisions):
            mu = recall_means[:, d]
            sigma = recall_stds[:, d]
            probs[d, :] = 1.0 - norm.cdf(tau, loc=mu, scale=sigma)

        preds = predict_division_clf(probs, xc, self.fractions_)
        return preds


class GPCLFRecallModeller(RecallModeller):
    """
    Per-division Gaussian Process with CLF threshold for recall cutoff prediction.

    Trains per-division GP regressors, then uses posterior probabilities
    P(recall > τ) and a CLF threshold X to predict the cutoff fraction.
    """

    model_save_filename = "recall_gp_clf_pipeline.pkl"

    def __init__(
        self,
        data_set_divide=5,
        model_interface=None,
        tax_level="order",
        feature_transformer=None,
        sort_strategy="reads",
        description: str | None = None,
    ):
        super().__init__(
            data_set_divide=data_set_divide,
            model_interface=model_interface,
            tax_level=tax_level,
            feature_transformer=feature_transformer,
            sort_strategy=sort_strategy,
            description=description or "GP+CLF recall model",
        )
        self._pipeline: GPCLFThreshold | None = None
        self.optimal_tau = None
        self.optimal_X = None
        self.optimal_loss = None
        self.fractions = np.arange(1, data_set_divide + 1) / data_set_divide

    @property
    def pipeline(self):
        return self._pipeline

    @pipeline.setter
    def pipeline(self, value):
        self._pipeline = value

    @property
    def model(self):
        return self._pipeline

    @model.setter
    def model(self, value):
        self._pipeline = value

    def fit(
        self,
        m_stats_matrices,
        taxids_to_use,
        test_size=0.2,
        random_state=42,
        optimize=True,
    ):
        self.transformer.fit(m_stats_matrices, taxids_to_use=taxids_to_use)
        rows = [self.transformer.transform(m) for m in m_stats_matrices]
        full = pd.concat(rows, ignore_index=True)

        feat_cols = self.transformer.get_feature_names_out()
        target_cols = self.transformer.target_columns_

        X = full[feat_cols]
        Y = full[target_cols]

        logger = logging.getLogger(__name__)
        logger.info(
            "GPCLFRecallModeller.fit: %d samples, %d features, %d targets",
            len(full),
            len(feat_cols),
            len(target_cols),
        )

        from sklearn.model_selection import train_test_split

        X_train, X_test, Y_train, Y_test = train_test_split(
            X,
            Y,
            test_size=test_size,
            random_state=random_state,
        )
        self.X_test = X_test
        self.y_test = Y_test

        self.pipeline = GPCLFThreshold(
            tau=0.9,
            x_thresh=0.5,
            n_divisions=self.data_set_divide,
            optimize=optimize,
        )
        self.pipeline.fit(
            X_train,
            Y_train,
            target_names=target_cols,
        )

        self.optimal_tau = self.pipeline.optimal_tau_
        self.optimal_X = self.pipeline.optimal_X_
        self.optimal_loss = self.pipeline.optimal_loss_

        self.RecP_feature_cols = feat_cols
        self.RecP_target_cols = target_cols
        self.date_trained = datetime.now(timezone.utc).isoformat()
        return self

    def predict_cutoff(self, X, target_recall=None, confidence=None):
        X_feat = self.transformer.transform(X)
        X_feat = X_feat[self.transformer.get_feature_names_out()]
        tau = target_recall if target_recall is not None else (self.optimal_tau or 0.9)
        x_thresh = confidence if confidence is not None else (self.optimal_X or 0.5)
        preds = self.pipeline.predict(
            X_feat,
            tau=tau,
            x_thresh=x_thresh,
        )
        return float(preds[0])

    def predict(self, X):
        X_feat = self.transformer.transform(X)
        X_feat = X_feat[self.transformer.get_feature_names_out()]
        return self.pipeline.predict_raw(X_feat)

    def save_model(self, output_directory):
        if self.pipeline is not None:
            bundle = {
                "model_type": "gp_clf",
                "model_category": self.model_category,
                "description": self.description,
                "date_trained": self.date_trained,
                "pipeline": self.pipeline,
                "feature_names": list(self.RecP_feature_cols),
                "target_names": list(self.RecP_target_cols),
                "data_set_divide": self.data_set_divide,
                "tax_level": self.tax_level,
                "optimal_tau": self.optimal_tau,
                "optimal_X": self.optimal_X,
                "optimal_loss": self.optimal_loss,
                "transformer": self.transformer,
            }
            joblib.dump(bundle, os.path.join(output_directory, self.model_save_filename))

    def load_model(self, input_directory):
        try:
            bundle = joblib.load(os.path.join(input_directory, self.model_save_filename))
            if isinstance(bundle, dict):
                self.pipeline = bundle["pipeline"]
                self.RecP_feature_cols = bundle["feature_names"]
                self.RecP_target_cols = bundle["target_names"]
                self.data_set_divide = bundle["data_set_divide"]
                self.tax_level = bundle.get("tax_level", self.tax_level)
                self.description = bundle.get("description", self.description)
                self.date_trained = bundle.get("date_trained", self.date_trained)
                self.optimal_tau = bundle["optimal_tau"]
                self.optimal_X = bundle["optimal_X"]
                self.optimal_loss = bundle["optimal_loss"]
                self.fractions = self.pipeline.fractions_
                self.transformer = bundle["transformer"]
            else:
                self.pipeline = bundle
        except Exception as e:
            print(f"Error loading model: {e}")

    def model_summary(self, model, X_test, Y_test, analysis_output_filedir):
        from sklearn.metrics import mean_squared_error, r2_score

        Y_pred = model.predict_raw(X_test)
        recall_cols = [c for c in self.RecP_target_cols if c.startswith("index_recall_")]
        if recall_cols:
            Y_test_recall = Y_test[recall_cols].values
            sample_max = Y_test_recall.max(axis=1, keepdims=True)
            sample_max = np.maximum(sample_max, 1e-6)
            Y_test_norm = Y_test_recall / sample_max
            Y_pred_recall = Y_pred[:, [i for i, n in enumerate(model.target_names_) if n.startswith("index_recall_")]]
            r2 = r2_score(Y_test_norm, Y_pred_recall, multioutput="uniform_average")
            mse = mean_squared_error(Y_test_norm, Y_pred_recall, multioutput="uniform_average")
        else:
            r2 = r2_score(Y_test, Y_pred, multioutput="uniform_average")
            mse = mean_squared_error(Y_test, Y_pred, multioutput="uniform_average")

        print(f"GPCLF model_summary R² = {r2:.3f}, MSE = {mse:.3f} (normalised recall)")
        if self.optimal_tau is not None:
            print(f"  Optimal τ={self.optimal_tau:.3f}, X={self.optimal_X:.3f}, loss={self.optimal_loss:.4f}")

        # ── True recall at predicted best cutoff ──
        if recall_cols:
            y_pred_frac = model.predict(X_test)
            y_recall = Y_test_recall.T
            n_test = y_recall.shape[1]
            n_div = self.data_set_divide
            tau = self.optimal_tau or 0.9

            actual_recalls = np.array(
                [
                    actual_recall_at_fraction(y_recall[:, s], y_pred_frac[s], n_div, self.fractions)
                    for s in range(n_test)
                ]
            )

            frac_below_tau = (actual_recalls < tau).mean()

            y_true_frac = np.array(
                [compute_ground_truth_division(y_recall[:, s], tau, n_div, self.fractions) for s in range(n_test)]
            )
            test_loss = asymmetric_loss(y_true_frac, y_pred_frac, over_cost=3, under_cost=1)

            print(
                f"  True recall at predicted cutoff: mean={actual_recalls.mean():.3f}, "
                f"median={np.median(actual_recalls):.3f}, std={actual_recalls.std():.3f}"
            )
            print(f"  % below τ={tau:.3f}: {frac_below_tau:.1%}")
            print(f"  Test loss (3:1) = {test_loss:.4f}")

        output_path = os.path.join(analysis_output_filedir, "recall_model_analysis_results.txt")
        with open(output_path, "w") as f:
            f.write(f"GP-CLF Recall Model\nR² = {r2:.4f}\nMSE = {mse:.4f}\n")
            if self.optimal_tau is not None:
                f.write(
                    f"Optimal τ = {self.optimal_tau:.4f}\n"
                    f"Optimal X = {self.optimal_X:.4f}\n"
                    f"Optimal loss (3:1) = {self.optimal_loss:.4f}\n"
                )
            if recall_cols:
                f.write(
                    f"True recall at predicted cutoff (test set):\n"
                    f"  mean = {actual_recalls.mean():.4f}\n"
                    f"  median = {np.median(actual_recalls):.4f}\n"
                    f"  std = {actual_recalls.std():.4f}\n"
                    f"  min = {actual_recalls.min():.4f}\n"
                    f"  % below τ={tau:.3f} = {frac_below_tau:.1%}\n"
                    f"  Test asymmetric loss (3:1) = {test_loss:.4f}\n"
                )
        return {
            "r2": r2,
            "mse": mse,
            "test_loss": test_loss,
            "recall_at_cutoff_mean": float(actual_recalls.mean()),
            "recall_at_cutoff_below_tau": float(frac_below_tau),
        }, {}

    def plot_diagnostics(self, output_dir):
        """Generate diagnostic plots: recall landscape, calibration, actual-recall-at-index."""
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns

        if self.X_test is None or self.y_test is None:
            print("No test data available for diagnostics plots")
            return

        sns.set_style("whitegrid")
        plt.rcParams["figure.dpi"] = 120

        recall_cols = [c for c in self.RecP_target_cols if c.startswith("index_recall_")]
        y_recall = self.y_test[recall_cols].values.T
        n_test = y_recall.shape[1]
        n_div = self.data_set_divide

        # Posterior predictions on test set
        means, stds = self.pipeline._predict_all(
            self.pipeline.scaler_.transform(
                self.X_test.values if hasattr(self.X_test, "values") else np.asarray(self.X_test)
            )
        )
        tau = self.optimal_tau or 0.9
        x_thresh = self.optimal_X or 0.5

        recall_indices = self.pipeline.recall_indices_
        probs = np.zeros((n_div, n_test))
        for d in range(n_div):
            idx = recall_indices[d]
            mu = means[:, idx]
            sigma = np.maximum(stds[:, idx], 1e-6)
            probs[d, :] = 1.0 - norm.cdf(tau, loc=mu, scale=sigma)

        y_pred = predict_division_clf(probs, x_thresh, self.fractions)
        y_true = np.array(
            [compute_ground_truth_division(y_recall[:, s], tau, n_div, self.fractions) for s in range(n_test)]
        )

        # —— 1. Recall landscape ——
        fig, ax = plt.subplots(figsize=(12, 5))
        bp_data = [y_recall[d, :] for d in range(n_div)]
        ax.boxplot(
            bp_data,
            positions=range(1, n_div + 1),
            widths=0.6,
            patch_artist=True,
            boxprops=dict(facecolor="lightblue", alpha=0.7),
            medianprops=dict(color="red", lw=2),
        )
        means_line = [np.mean(y_recall[d, :]) for d in range(n_div)]
        ax.plot(range(1, n_div + 1), means_line, "ko-", lw=2, markersize=5, label="Mean recall")
        ax.axhline(tau, color="green", ls="--", lw=1.5, label=f"τ = {tau}")
        ax.set_xlabel("Division")
        ax.set_ylabel("Recall")
        ax.set_title("Recall landscape (test set)")
        ax.set_xticks(range(1, n_div + 1))
        ax.set_xticklabels([f"{i}\n({self.fractions[i - 1]:.0%})" for i in range(1, n_div + 1)], fontsize=7)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, axis="y")
        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, "recall_landscape.png"), dpi=120, bbox_inches="tight")
        plt.close(fig)
        print("  Saved recall_landscape.png")

        # —— 2. Calibration scatter ——
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        ax = axes[0]
        over, under = y_pred > y_true, y_pred <= y_true
        ax.scatter(y_true[over], y_pred[over], c="red", alpha=0.5, s=20, label=f"Over ({over.sum()})")
        ax.scatter(y_true[under], y_pred[under], c="blue", alpha=0.5, s=20, label=f"Under ({under.sum()})")
        ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
        ax.set_xlabel("True fraction")
        ax.set_ylabel("Predicted fraction")
        ax.set_title(f"GP-CLF (τ={tau:.2f}, X={x_thresh:.2f})")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal")

        ax = axes[1]
        errors = y_pred - y_true
        ax.hist(errors, bins=30, alpha=0.7, color="gray", edgecolor="black")
        ax.axvline(0, color="black", ls="--", lw=1)
        ax.axvline(np.median(errors), color="red", ls="-", lw=1.5, label=f"Median: {np.median(errors):.3f}")
        ax.set_xlabel("Predicted − True")
        ax.set_ylabel("Count")
        loss_val = self.optimal_loss if self.optimal_loss is not None else asymmetric_loss(y_true, y_pred)
        ax.set_title(f"Error distribution (loss={loss_val:.4f})")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, "recall_calibration.png"), dpi=120, bbox_inches="tight")
        plt.close(fig)
        print("  Saved recall_calibration.png")

        # —— 3. Actual recall at predicted index ——
        actual_recalls = np.array(
            [actual_recall_at_fraction(y_recall[:, s], y_pred[s], n_div, self.fractions) for s in range(n_test)]
        )

        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        ax = axes[0]
        ax.hist(actual_recalls, bins=25, alpha=0.7, color="steelblue", edgecolor="black", density=True)
        ax.axvline(tau, color="red", ls="--", lw=2, label=f"τ = {tau}")
        ax.axvline(
            np.mean(actual_recalls), color="darkorange", ls="-", lw=2, label=f"Mean = {np.mean(actual_recalls):.3f}"
        )
        ax.axvline(
            np.median(actual_recalls), color="green", ls=":", lw=2, label=f"Median = {np.median(actual_recalls):.3f}"
        )
        ax.set_xlabel("Actual recall at predicted last-index")
        ax.set_ylabel("Density")
        ax.set_title(f"GP-CLF (τ={tau:.2f}, X={x_thresh:.2f})")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, axis="y")

        ax = axes[1]
        n_show = min(50, n_test)
        rng = np.random.RandomState(42)
        sample_idxs = rng.choice(n_test, n_show, replace=False)
        for s in sample_idxs:
            ax.plot(self.fractions, y_recall[:, s], color="grey", alpha=0.15, lw=0.5)
        ax.plot(self.fractions, means_line, "k-", lw=2, label="Mean recall")
        jitter = rng.uniform(-0.01, 0.01, n_test)
        for s in range(n_test):
            color = "red" if actual_recalls[s] < tau else "blue"
            ax.scatter(y_pred[s] + jitter[s], actual_recalls[s], c=color, alpha=0.4, s=8, edgecolors="none")
        ax.axhline(tau, color="red", ls="--", lw=1.5)
        ax.set_xlabel("Predicted last-index (fraction)")
        ax.set_ylabel("Actual recall at that index")
        ax.set_title("Readout: predicted vs actual recall")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        fig.savefig(os.path.join(output_dir, "recall_actual_at_index.png"), dpi=120, bbox_inches="tight")
        plt.close(fig)
        print("  Saved recall_actual_at_index.png")


def find_recall_for_target(
    recall_pred_row: pd.Series,
    data_set_divide: int,
    target_recall: float = 0.95,
) -> float:
    """
    Find the index percentile at which recall reaches target_recall.
    Uses round-up to next percentile.
    """
    for i in range(1, data_set_divide + 1):
        recall_at_percentile = recall_pred_row.get(f"index_recall_{i}", 0)
        if recall_at_percentile >= target_recall:
            return recall_at_percentile
    return 1.0


def find_percentile_for_recall(
    recall_pred_row: pd.Series,
    data_set_divide: int,
    target_recall: float = 0.95,
) -> float:
    """
    Find the index percentile at which recall reaches target_recall.
    Uses round-up to next percentile.
    """
    for i in range(1, data_set_divide + 1):
        recall_at_percentile = recall_pred_row.get(f"index_recall_{i}", 0)
        if recall_at_percentile >= target_recall:
            return i / data_set_divide
    return 1.0


def cut_off_recall_prediction(
    study_output_filepath: Path,
    data_set_name: str,
    modeller: RecallModeller | CutoffRecallModeller,
    data_set_divide: int,
    m_stats_stats_matrix: pd.DataFrame,
    taxids_to_use: pd.DataFrame,
    tax_level: str = "order",
    target_recall: float = 1.0,
    confidence: float | None = None,
) -> tuple[OverlapManager, dict]:
    """
    Predict recall at various cutoffs and filter leaves based on threshold.

    Works with both RecallModeller (per-bin) and CutoffRecallModeller (direct cutoff).

    Args:
        study_output_filepath: Path to study output directory.
        data_set_name: Name of the dataset.
        modeller: Trained RecallModeller or CutoffRecallModeller.
        data_set_divide: Number of recall divisions used in model.
        m_stats_stats_matrix: Stats matrix for the dataset.
        taxids_to_use: DataFrame of taxids to use.
        tax_level: Taxonomic level for analysis.
        target_recall: Target recall threshold.
        confidence: Confidence level for prob-guided cutoff (CutoffRecallModeller only).

    Returns:
        Tuple of (OverlapManager, metrics dict).
    """
    from metagenomics_utils.overlap_manager import OverlapManager
    target_percentile = modeller.predict_cutoff(m_stats_stats_matrix, target_recall, confidence)
    keep_index = round(target_percentile * m_stats_stats_matrix.shape[0])
    if keep_index == 0:
        keep_index = 1

    original_matrix = m_stats_stats_matrix.copy()
    total_with_coverage = (original_matrix["coverage"] > 0).sum()
    kept_matrix = original_matrix.head(keep_index)
    kept_with_coverage = (kept_matrix["coverage"] > 0).sum() if not kept_matrix.empty else 0
    filtered_with_coverage = total_with_coverage - kept_with_coverage

    prop_coverage_above = kept_with_coverage / total_with_coverage if total_with_coverage > 0 else 0
    prop_coverage_below = filtered_with_coverage / total_with_coverage if total_with_coverage > 0 else 0
    resulting_percentile = keep_index / original_matrix.shape[0] if original_matrix.shape[0] > 0 else 0

    overlap_manager = OverlapManager(
        os.path.join(study_output_filepath, f"{data_set_name}", "clustering"), max_proportion=target_percentile
    )

    metrics = {
        "target_percentile": target_percentile,
        "resulting_percentile": resulting_percentile,
        "keep_index": keep_index,
        "prop_coverage_above_cutoff": prop_coverage_above,
        "prop_coverage_below_cutoff": prop_coverage_below,
    }

    eval_metrics = modeller.evaluate_prediction(m_stats_stats_matrix, target_recall)
    metrics.update(eval_metrics)

    return overlap_manager, metrics


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


class ClusteringPipeline(BaseEstimator, ClassifierMixin):
    """Sklearn-compatible XGBoost classifier with StandardScaler and Optuna support.

    Scales numeric columns, then trains/predicts with XGBoost for clade traversal decisions.
    Supports Optuna-based hyperparameter optimization via the same search space as
    CompositionModeller.xgbc_model_bayes_optimized().
    """

    def __init__(
        self,
        numeric_cols=None,
        taxon_cols=None,
        feature_names=None,
        scale_pos_weight="auto",
        optuna_trials=50,
        use_optuna=True,
        random_state=42,
        xgb_params=None,
    ):
        self.numeric_cols = numeric_cols or []
        self.taxon_cols = taxon_cols or []
        self.feature_names = feature_names
        self.scale_pos_weight = scale_pos_weight
        self.optuna_trials = optuna_trials
        self.use_optuna = use_optuna
        self.random_state = random_state
        self.xgb_params = xgb_params or {}

    def fit(self, X: pd.DataFrame | np.ndarray, y: pd.Series | np.ndarray):
        if isinstance(X, pd.DataFrame):
            self.feature_names_ = X.columns.tolist()
            self.n_features_in_ = X.shape[1]
        else:
            self.feature_names_ = self.feature_names
            self.n_features_in_ = X.shape[1] if X.ndim == 2 else 1

        if self.scale_pos_weight == "auto":
            n_neg = (y == 0).sum()
            n_pos = (y == 1).sum()
            pos_weight = n_neg / max(n_pos, 1)
        else:
            pos_weight = self.scale_pos_weight

        self.scaler_ = StandardScaler()
        if isinstance(X, pd.DataFrame) and self.numeric_cols:
            X_scaled = X.values.copy()
            num_idx = [self.feature_names_.index(c) for c in self.numeric_cols]
            X_scaled[:, num_idx] = self.scaler_.fit_transform(X_scaled[:, num_idx])
        else:
            n_num = len(self.numeric_cols)
            X_arr = np.asarray(X, dtype=float).copy()
            if n_num > 0:
                X_arr[:, :n_num] = self.scaler_.fit_transform(X_arr[:, :n_num])
            else:
                X_arr = self.scaler_.fit_transform(X_arr)
            X_scaled = X_arr

        if self.use_optuna:
            self.classifier_, self.study_ = self._train_with_optuna(X_scaled, y, pos_weight)
        else:
            self.classifier_ = XGBClassifier(
                scale_pos_weight=pos_weight,
                random_state=self.random_state,
                use_label_encoder=False,
                eval_metric="logloss",
                **self.xgb_params,
            )
            self.classifier_.fit(X_scaled, y)

        return self

    def __getstate__(self):
        state = super().__getstate__()
        state.pop("study_", None)
        state.pop("taxon_cols", None)
        state.pop("optuna_trials", None)
        state.pop("use_optuna", None)
        state.pop("scale_pos_weight", None)
        state.pop("xgb_params", None)
        state.pop("feature_names", None)
        state.pop("feature_names_", None)
        state.pop("random_state", None)
        return state

    def __setstate__(self, state):
        super().__setstate__(state)
        self.study_ = None
        self.taxon_cols = []
        self.optuna_trials = 0
        self.use_optuna = False
        self.scale_pos_weight = "auto"
        self.xgb_params = {}
        self.feature_names = None
        self.feature_names_ = None
        self.random_state = 42

    def _train_with_optuna(self, X, y, pos_weight):
        import optuna
        from sklearn.model_selection import StratifiedKFold, cross_val_score

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
                "scale_pos_weight": pos_weight,
                "random_state": self.random_state,
                "n_jobs": -1,
            }
            model = XGBClassifier(**params)
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=self.random_state)
            scores = cross_val_score(model, X, y, cv=cv, scoring="roc_auc", n_jobs=-1)
            return scores.mean()

        study = optuna.create_study(direction="maximize")
        study.optimize(objective, n_trials=self.optuna_trials, show_progress_bar=True)

        best_params = study.best_trial.params
        best_model = XGBClassifier(**best_params)
        best_model.fit(X, y)
        return best_model, study

    def predict(self, X):
        X_scaled = self._scale(X)
        clf = getattr(self, "classifier_", None) or getattr(self, "model", None)
        return clf.predict(X_scaled)

    def predict_proba(self, X):
        X_scaled = self._scale(X)
        clf = getattr(self, "classifier_", None) or getattr(self, "model", None)
        return clf.predict_proba(X_scaled)

    def _scale(self, X):
        X = np.asarray(X, dtype=float).copy()
        if X.ndim == 1:
            X = X.reshape(1, -1)
        n_num = len(self.numeric_cols)
        scaler = getattr(self, "scaler_", None) or getattr(self, "scaler", None)
        if scaler is None:
            return X
        if n_num > 0:
            X[:, :n_num] = scaler.transform(X[:, :n_num])
        return X

    def get_params(self, deep=True):
        return {
            "numeric_cols": self.numeric_cols,
            "taxon_cols": self.taxon_cols,
            "feature_names": self.feature_names,
            "scale_pos_weight": self.scale_pos_weight,
            "optuna_trials": self.optuna_trials,
            "use_optuna": self.use_optuna,
            "random_state": self.random_state,
            "xgb_params": self.xgb_params,
        }

    def set_params(self, **params):
        for k, v in params.items():
            setattr(self, k, v)
        return self


from abc import ABC, abstractmethod

STATS_COLS = ["n_leaves", "tax_diversity", "Min_Dist", "Min_Shared"]
DROPPED_COLS = [
    "data_set",
    "node",
    "n_true_leaves",
    "precision_increased",
    "new_precision",
    "precision",
    "stop_traversal",
    "unclassified",
]


class BaseCompositionModeller(ABC):
    """Base class for composition (stop_traversal) classifiers.

    Subclasses define _build_pipeline() to return a sklearn-compatible
    estimator with .fit(), .predict(), .predict_proba().
    """

    model_save_filename: str = "composition_bundle.pkl"
    model_category: str = "composition"

    def __init__(self, tax_level: str = "order", description: str | None = None):
        self.tax_level = tax_level
        self.description = description or f"{type(self).__name__} composition model"
        self.date_trained: str | None = None
        self.pipeline: BaseEstimator | Pipeline | None = None
        self.X_train = None
        self.X_test = None
        self.y_train = None
        self.y_test = None
        self._feature_names = None

    def __getstate__(self):
        state = self.__dict__.copy()
        state.pop("X_train", None)
        state.pop("X_test", None)
        state.pop("y_train", None)
        state.pop("y_test", None)
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.X_train = None
        self.X_test = None
        self.y_train = None
        self.y_test = None

    @abstractmethod
    def _build_pipeline(self, X_train, y_train):
        """Return a fitted sklearn Pipeline or compatible estimator."""

    def fit(self, X_train, y_train, X_test=None, y_test=None):
        self.X_train = X_train
        self.y_train = y_train
        self.X_test = X_test
        self.y_test = y_test
        self._feature_names = list(X_train.columns) if hasattr(X_train, "columns") else None
        self.pipeline = self._build_pipeline(X_train, y_train)
        self.date_trained = datetime.now(timezone.utc).isoformat()
        return self

    @property
    def model(self):
        """Return the inner estimator (for backward compat with traversal_with_prediction)."""
        return self.pipeline

    def predict_proba(self, X):
        return self.pipeline.predict_proba(X)

    def predict(self, X):
        return self.pipeline.predict(X)

    def save_model(self, output_directory: str):
        if self.pipeline is not None:
            bundle = {
                "model_type": type(self).__name__,
                "model_category": self.model_category,
                "description": self.description,
                "date_trained": self.date_trained,
                "pipeline": self.pipeline,
                "feature_names": self._feature_names,
                "tax_level": self.tax_level,
                "X_train": self.X_train,
                "X_test": self.X_test,
                "y_train": self.y_train,
                "y_test": self.y_test,
            }
            joblib.dump(bundle, os.path.join(output_directory, self.model_save_filename))
        else:
            print("No model to save.")

    def load_model(self, output_directory: str):
        try:
            data = joblib.load(os.path.join(output_directory, self.model_save_filename))
            if isinstance(data, dict):
                self.pipeline = data["pipeline"]
                self._feature_names = data.get("feature_names")
                self.tax_level = data.get("tax_level", self.tax_level)
                self.description = data.get("description", self.description)
                self.date_trained = data.get("date_trained", self.date_trained)
                self.X_train = data.get("X_train")
                self.X_test = data.get("X_test")
                self.y_train = data.get("y_train")
                self.y_test = data.get("y_test")
            else:
                self.pipeline = data
                self._feature_names = None
        except Exception as e:
            print(f"Error loading model: {e}")

    def evaluate_model(self, X_test=None, y_test=None):
        from sklearn.metrics import classification_report, confusion_matrix

        X_test = X_test if X_test is not None else self.X_test
        y_test = y_test if y_test is not None else self.y_test
        if X_test is None or y_test is None:
            raise ValueError("X_test and y_test must be provided for evaluation.")

        if self.pipeline is None:
            raise ValueError("Model pipeline is not fitted. Call fit() before evaluation.")
        y_pred = self.pipeline.predict(X_test)
        report = classification_report(y_test, y_pred, output_dict=True)
        cm = confusion_matrix(y_test, y_pred)
        return report, cm

    def _pos_weight(self, y):
        n_neg = (y == 0).sum()
        n_pos = (y == 1).sum()
        return n_neg / max(n_pos, 1)

    def _make_column_transformer(self, X_train, stats_cols=None, scaler=True):
        """Build ColumnTransformer that scales stats cols, leaves tax cols raw."""
        from sklearn.compose import ColumnTransformer
        from sklearn.preprocessing import StandardScaler

        stats_cols = stats_cols or [c for c in STATS_COLS if c in X_train.columns]
        tax_cols = [c for c in X_train.columns if c not in stats_cols]
        transformers = []
        if stats_cols:
            transformers.append(("scaler", StandardScaler() if scaler else "passthrough", stats_cols))
        if tax_cols:
            transformers.append(("passthrough", "passthrough", tax_cols))
        return ColumnTransformer(transformers, remainder="drop")

    def plot_eval(self, output_directory):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        classifier = self._get_classifier()
        if classifier is None or not hasattr(classifier, "feature_importances_"):
            return

        if self._feature_names is None:
            return

        importances = pd.Series(classifier.feature_importances_, index=self._feature_names)
        importances.sort_values(ascending=True).tail(20).plot.barh(figsize=(8, 8))
        plt.title("Top 20 Feature Importances")
        plt.tight_layout()
        plt.savefig(os.path.join(output_directory, "composition_feature_importance.png"))
        plt.close()

    def shap_eval_plot(self, output_directory):
        try:
            import matplotlib
            import shap

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            classifier = self._get_classifier()
            if classifier is None or not hasattr(classifier, "feature_importances_"):
                return
            if self.X_train is None or self.X_test is None:
                return

            explainer = shap.Explainer(classifier.predict, self.X_train)
            shap_values = explainer(self.X_test)

            plt.figure(figsize=(10, 7))
            shap.summary_plot(shap_values, self.X_test, plot_type="dot", show=False)
            plt.savefig(os.path.join(output_directory, "shap_summary_plot.png"))
            plt.close()

            plt.figure(figsize=(10, 7))
            shap.summary_plot(shap_values, self.X_test, plot_type="bar", show=False)
            plt.savefig(os.path.join(output_directory, "shap_bar_plot.png"))
            plt.close()

            plt.figure(figsize=(10, 7))
            shap.dependence_plot("Min_Shared", shap_values.values, self.X_test, show=False)
            plt.savefig(os.path.join(output_directory, "shap_dependence_plot.png"))
            plt.close()
        except Exception:
            pass

    def shap_interaction_plot(self, output_directory):
        try:
            import matplotlib
            import shap

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import numpy as np
            import seaborn as sns

            classifier = self._get_classifier()
            if classifier is None or not hasattr(classifier, "feature_importances_"):
                return
            if self.X_train is None or self.X_test is None:
                return

            explainer = shap.TreeExplainer(classifier)
            interaction_values = explainer.shap_interaction_values(self.X_test)
            interaction_strength = np.abs(interaction_values).mean(axis=0)
            interaction_df = pd.DataFrame(
                interaction_strength,
                index=self.X_train.columns,
                columns=self.X_train.columns,
            )
            interaction_df = interaction_df.where(np.tril(np.ones(interaction_df.shape), k=-1).astype(bool))
            plt.figure(figsize=(10, 8))
            sns.heatmap(interaction_df, cmap="viridis")
            plt.title("Mean absolute SHAP interaction values")
            plt.savefig(os.path.join(output_directory, "shap_interaction_heatmap.png"))
            plt.close()

            from scipy.cluster.hierarchy import dendrogram, linkage

            dist_df = interaction_df.drop(columns=STATS_COLS, index=STATS_COLS, errors="ignore")
            np.fill_diagonal(dist_df.values, 0)
            dist_df = dist_df.fillna(0)
            Z = linkage(dist_df, method="average")
            plt.figure(figsize=(6, 7))
            dendrogram(Z, labels=dist_df.index, orientation="left", leaf_rotation=0)
            plt.title("Feature Clustering based on SHAP Interaction Values")
            plt.tight_layout()
            plt.savefig(os.path.join(output_directory, "shap_interaction_dendrogram.png"))
            plt.close()
        except Exception:
            pass

    def eval_and_plot(self, X_test, y_test, output_directory, X_train=None):
        self.X_test = X_test
        self.y_test = y_test
        if X_train is not None:
            self.X_train = X_train
        report, cm = self.evaluate_model(X_test, y_test)
        self.plot_eval(output_directory)
        self.shap_eval_plot(output_directory)
        self.shap_interaction_plot(output_directory)
        return report, cm

    def _get_classifier(self) -> BaseEstimator | Pipeline | None:
        """Return the inner classifier from the pipeline."""
        if self.pipeline is None:
            return None
        if hasattr(self.pipeline, "classifier_"):
            return self.pipeline.classifier_
        if hasattr(self.pipeline, "named_steps"):
            return self.pipeline[-1] if hasattr(self.pipeline[-1], "predict") else self.pipeline
        return self.pipeline


class XGBCompositionModeller(BaseCompositionModeller):
    """XGBoost classifier with fixed hyperparams (no Optuna)."""

    model_save_filename = "composition_xgb_bundle.pkl"

    def __init__(
        self, n_estimators=300, max_depth=6, learning_rate=0.1, subsample=0.8, colsample_bytree=0.8, random_state=42,
        tax_level="order",
    ):
        super().__init__(tax_level=tax_level)
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.subsample = subsample
        self.colsample_bytree = colsample_bytree
        self.random_state = random_state

    def _build_pipeline(self, X_train, y_train) -> Pipeline:
        pos_weight = self._pos_weight(y_train)
        ct = self._make_column_transformer(X_train)
        return Pipeline(
            [
                ("preprocessor", ct),
                (
                    "classifier",
                    XGBClassifier(
                        n_estimators=self.n_estimators,
                        max_depth=self.max_depth,
                        learning_rate=self.learning_rate,
                        subsample=self.subsample,
                        colsample_bytree=self.colsample_bytree,
                        scale_pos_weight=pos_weight,
                        random_state=self.random_state,
                        eval_metric="logloss",
                    ),
                ),
            ]
        ).fit(X_train, y_train)


class OptunaXGBCompositionModeller(BaseCompositionModeller):
    """XGBoost classifier with Optuna hyperparameter optimisation
    (replicates the original CompositionModeller behaviour).
    """

    model_save_filename = "composition_optuna_bundle.pkl"

    def __init__(self, optuna_trials=50, random_state=42, tax_level="order"):
        super().__init__(tax_level=tax_level)
        self.optuna_trials = optuna_trials
        self.random_state = random_state

    def _build_pipeline(self, X_train, y_train):
        stats_cols = [c for c in STATS_COLS if c in X_train.columns]
        tax_cols = [c for c in X_train.columns if c not in stats_cols]
        pipeline = ClusteringPipeline(
            numeric_cols=stats_cols,
            taxon_cols=tax_cols,
            feature_names=list(X_train.columns),
            use_optuna=True,
            optuna_trials=self.optuna_trials,
            random_state=self.random_state,
        )
        pipeline.fit(X_train, y_train)
        return pipeline


class RFCompositionModeller(BaseCompositionModeller):
    """Random Forest classifier."""

    model_save_filename = "composition_rf_bundle.pkl"

    def __init__(self, n_estimators=300, max_depth=12, min_samples_leaf=3, random_state=42, tax_level="order"):
        super().__init__(tax_level=tax_level)
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.min_samples_leaf = min_samples_leaf
        self.random_state = random_state

    def _build_pipeline(self, X_train, y_train):
        from sklearn.pipeline import Pipeline

        ct = self._make_column_transformer(X_train)
        return Pipeline(
            [
                ("preprocessor", ct),
                (
                    "classifier",
                    RandomForestClassifier(
                        n_estimators=self.n_estimators,
                        max_depth=self.max_depth,
                        min_samples_leaf=self.min_samples_leaf,
                        class_weight="balanced",
                        random_state=self.random_state,
                        n_jobs=-1,
                    ),
                ),
            ]
        ).fit(X_train, y_train)


class GBCompositionModeller(BaseCompositionModeller):
    """Gradient Boosting classifier."""

    model_save_filename = "composition_gb_bundle.pkl"

    def __init__(self, n_estimators=300, max_depth=5, learning_rate=0.1, subsample=0.8, random_state=42, tax_level="order"):
        super().__init__(tax_level=tax_level)
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.subsample = subsample
        self.random_state = random_state

    def _build_pipeline(self, X_train, y_train):
        from sklearn.ensemble import GradientBoostingClassifier
        from sklearn.pipeline import Pipeline

        ct = self._make_column_transformer(X_train)
        return Pipeline(
            [
                ("preprocessor", ct),
                (
                    "classifier",
                    GradientBoostingClassifier(
                        n_estimators=self.n_estimators,
                        max_depth=self.max_depth,
                        learning_rate=self.learning_rate,
                        subsample=self.subsample,
                        random_state=self.random_state,
                    ),
                ),
            ]
        ).fit(X_train, y_train)


class LRCompositionModeller(BaseCompositionModeller):
    """Logistic Regression classifier (stats-only features)."""

    model_save_filename = "composition_lr_bundle.pkl"

    def __init__(self, C=1.0, max_iter=1000, random_state=42, tax_level="order"):
        super().__init__(tax_level=tax_level)
        self.C = C
        self.max_iter = max_iter
        self.random_state = random_state

    def _build_pipeline(self, X_train, y_train):
        from sklearn.compose import ColumnTransformer
        from sklearn.linear_model import LogisticRegression
        from sklearn.pipeline import Pipeline

        stats_cols = [c for c in STATS_COLS if c in X_train.columns]
        return Pipeline(
            [
                (
                    "preprocessor",
                    ColumnTransformer(
                        [
                            ("scaler", StandardScaler(), stats_cols),
                        ],
                        remainder="drop",
                    ),
                ),
                (
                    "classifier",
                    LogisticRegression(
                        C=self.C,
                        class_weight="balanced",
                        max_iter=self.max_iter,
                        random_state=self.random_state,
                    ),
                ),
            ]
        ).fit(X_train, y_train)

    def plot_eval(self, output_directory):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        classifier = self._get_classifier()
        if classifier is None or not hasattr(classifier, "coef_"):
            return

        coef = classifier.coef_.flatten()
        names = ["n_leaves", "tax_diversity", "Min_Dist", "Min_Shared"]
        pd.Series(np.abs(coef), index=names).sort_values().plot.barh(figsize=(7, 4))
        plt.title("Logistic Regression Coefficients (absolute)")
        plt.tight_layout()
        plt.savefig(os.path.join(output_directory, "lr_coefficients.png"))
        plt.close()

    def shap_eval_plot(self, output_directory):
        pass

    def shap_interaction_plot(self, output_directory):
        pass


#######################################################################################
################ CROSS-HIT MODELLING ######################################################
#######################################################################################


class CrossHitModeller:
    model_save_filename = "cross_hit_xgb_bundle.pkl"
    model_category: str = "crosshit"

    def __init__(self, prediction_trainning_results_df, description: str | None = None):
        self.description = description or "Cross-hit XGBoost classifier"
        self.date_trained: str | None = None
        self.prediction_trainning_results_df = prediction_trainning_results_df
        self.X = self.prediction_trainning_results_df.drop(columns=["leaf", "is_trash"])
        self.pred_stats_cols = ["coverage", "covbases", "meanmapq", "error_rate", "max_shared", "total_uniq_reads"]
        self.y = self.prediction_trainning_results_df["is_trash"].astype(int)
        self.scaler = None
        self.pca = None
        self.model = None

    def split_data(self, test_size=0.2, random_state=42):
        from sklearn.model_selection import train_test_split

        X_train, X_test, y_train, y_test = train_test_split(
            self.X, self.y, test_size=test_size, random_state=random_state
        )
        return X_train, X_test, y_train, y_test

    def prep_data(self, scale: bool = True, transform: bool = False):

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

        model = XGBClassifier(use_label_encoder=False, eval_metric="logloss", **kwargs)
        model.fit(X_train, y_train)
        return model

    def train_model(self, optimized: bool = True, **kwargs):
        if optimized:
            return self.train_model_bayes_optimized()
        X_train, X_test, y_train, y_test = self.prep_data()
        model = self.xgbc_model(X_train, y_train, **kwargs)
        self.model = model
        self.date_trained = datetime.now(timezone.utc).isoformat()
        return model, X_test, y_test

    def train_model_bayes_optimized(self):
        import optuna
        from sklearn.model_selection import StratifiedKFold, cross_val_score
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
                "n_jobs": -1,
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
        self.date_trained = datetime.now(timezone.utc).isoformat()
        return best_model, X_test, y_test, study

    def save_model(self, output_directory: str):

        if self.model is not None:
            joblib.dump(
                {
                    "model": self.model,
                    "model_category": self.model_category,
                    "description": self.description,
                    "date_trained": self.date_trained,
                    "scaler": self.scaler,
                    "pca": self.pca,
                },
                os.path.join(output_directory, self.model_save_filename),
            )
        else:
            print("No model to save.")

    def load_model(self, input_directory: str):

        try:
            bundle = joblib.load(os.path.join(input_directory, self.model_save_filename))
            if isinstance(bundle, dict):
                self.model = bundle["model"]
                self.description = bundle.get("description", self.description)
                self.date_trained = bundle.get("date_trained", self.date_trained)
                self.scaler = bundle.get("scaler")
                self.pca = bundle.get("pca")
            else:
                self.model = bundle
        except Exception as e:
            print(f"Error loading model: {e}")


########################################################################################
################ TRAVERSAL ######################################################
########################################################################################


def cross_hit_prediction(
    data_set_name,
    study_output_filepath,
    ncbi_wrapper,
    modeller: CrossHitModeller,
    overlap_manager: OverlapManager,
    tax_df,
    tax_level: str = "order",
):

    prediction_matrix = cross_hit_prediction_matrix(
        data_set_name, study_output_filepath, ncbi_wrapper, overlap_manager, tax_df, tax_level=tax_level
    )

    if prediction_matrix.empty or len(overlap_manager.leaves) == 0:
        return pd.DataFrame(columns=["leaf", "is_trash", "prob_best_match", "pred_best_match"])

    X_pred = prediction_matrix.drop(columns=["leaf", "is_trash"])
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
        X_pred_tax_pca.columns = [f"pca_{i + 1}" for i in range(X_pred_tax_pca.shape[1])]
        X_pred_scaled = pd.concat([X_pred_stats.reset_index(drop=True), X_pred_tax_pca.reset_index(drop=True)], axis=1)

    if modeller.model is None:
        y_prob = np.zeros(X_pred_scaled.shape[0])
    else:
        y_prob = modeller.model.predict_proba(X_pred_scaled)[:, 1]

    output = prediction_matrix[["leaf", "is_trash"]].copy()
    output.loc[:, "prob_best_match"] = y_prob
    output.loc[:, "pred_best_match"] = (y_prob > 0.5).astype(int)

    return output


def traversal_with_clustering_fixed(
    overlap_manager: OverlapManager,
    node: str,
    stats_matrix: pd.DataFrame,
    min_dist_threshold: float = 0.6,
    results=None,
) -> list[pd.DataFrame]:
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
    if results is None:
        results = []

    node_true_leaves = node_total_true_leaves(overlap_manager, node, stats_matrix)
    node_precision = 1 / len(set(node_true_leaves)) if len(node_true_leaves) > 0 else 0.0
    node_precision = 1 / node_precision if node_precision > 1 else node_precision
    node_leaf_taxids = node_leaves_best_taxids(overlap_manager, node, stats_matrix)
    node_row = overlap_manager.all_node_stats[overlap_manager.all_node_stats["Node"] == node]
    min_dist = node_row["Min_Pairwise_Dist"].values[0]
    min_shared = node_row["Min_Shared"].values[0]

    best_taxid_match = stats_matrix[stats_matrix.index.isin(overlap_manager.get_node_leaves(node))].copy()
    best_taxid_match = best_taxid_match.sort_values(by=["coverage"], ascending=False)["best_match_taxid"].tolist()
    best_taxid_match = best_taxid_match[0] if len(best_taxid_match) > 0 else None

    if min_shared >= min_dist_threshold:  # stop condition or leaf
        node_leaves = overlap_manager.get_node_leaves(node)
        results.append(
            {
                "node": node,
                "n_leaves": len(node_leaves),
                "leaves": node_leaves,
                "best_taxid_match": best_taxid_match,
                "node_precision": node_precision,
                "node_taxids": node_leaf_taxids,
            }
        )

    # Traverse children if prediction is positive
    else:
        for child in overlap_manager.tree.successors(node):
            if overlap_manager.tree.out_degree(child) > 0:  # internal node
                traversal_with_clustering_fixed(
                    overlap_manager, child, stats_matrix, min_dist_threshold=min_dist_threshold, results=results
                )
            else:
                best_taxid_match = stats_matrix[stats_matrix.index == child]["best_match_taxid"].tolist()
                best_taxid_match = best_taxid_match[0] if len(best_taxid_match) > 0 else None
                results.append(
                    {
                        "node": child,
                        "n_leaves": 1,
                        "leaves": [child],
                        "best_taxid_match": best_taxid_match,
                        "node_precision": 1.0,
                        "node_taxids": [best_taxid_match] if best_taxid_match is not None else [],
                    }
                )

    return results


def traversal_with_prediction(
    overlap_manager: OverlapManager,
    node: str,
    modeller: BaseCompositionModeller,
    stats_matrix,
    tax_df: pd.DataFrame,
    tax_level: str = "order",
    results=None,
) -> list[pd.DataFrame]:
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
    if results is None:
        results = []
    if modeller.model is None:
        raise ValueError("Modeller has no trained model. Please train the model before calling this function.")

    composition = (
        node_composition_level(overlap_manager, node, stats_matrix, tax_df, tax_level=tax_level)
        .set_index("tax_level")
        .T
    )
    composition = composition.reset_index(drop=True)

    node_true_leaves = node_total_true_leaves(overlap_manager, node, stats_matrix)
    node_precision = 1 / len(set(node_true_leaves)) if len(node_true_leaves) > 0 else 0.0
    node_precision = 1 / node_precision if node_precision > 1 else node_precision
    node_leaf_taxids = node_leaves_best_taxids(overlap_manager, node, stats_matrix)
    node_row = overlap_manager.all_node_stats[overlap_manager.all_node_stats["Node"] == node]
    min_dist = node_row["Min_Pairwise_Dist"].values[0]
    min_shared = node_row["Min_Shared"].values[0]

    node_total_leaf_taxa_div = node_leaf_shannon_tax_diversity(overlap_manager, node, stats_matrix, tax_level=tax_level)

    input_features = pd.DataFrame(
        {
            "n_leaves": [len(overlap_manager.get_node_leaves(node))],
            "tax_diversity": [node_total_leaf_taxa_div],
            "Min_Dist": [min_dist],
            "Min_Shared": [min_shared],
        }
    )

    input_features = pd.concat([input_features, composition], axis=1).drop(
        columns=["unclassified"], errors="ignore", axis=1
    )

    expected_columns = getattr(modeller.model, "feature_names_in_", None)
    if expected_columns is not None:
        input_features = input_features.reindex(columns=expected_columns, fill_value=0)

    stop_traversal_pred = modeller.model.predict(input_features)[0]
    if stop_traversal_pred is False:
        print(f"Stopping at node {node} with {len(overlap_manager.get_node_leaves(node))} leaves.")

    best_taxid_match = stats_matrix[stats_matrix.index.isin(overlap_manager.get_node_leaves(node))].copy()
    # best_taxid_match = best_taxid_match[best_taxid_match['is_trash'] == False]
    best_taxid_match = best_taxid_match.sort_values(by=["coverage"], ascending=False)["best_match_taxid"].tolist()

    best_taxid_match = best_taxid_match[0] if len(best_taxid_match) > 0 else None

    if stop_traversal_pred == True or overlap_manager.tree.out_degree(node) == 0:  # stop condition or leaf
        node_leaves = overlap_manager.get_node_leaves(node)
        results.append(
            {
                "node": node,
                "n_leaves": len(node_leaves),
                "leaves": node_leaves,
                "best_taxid_match": best_taxid_match,
                "node_precision": node_precision,
                "node_taxids": node_leaf_taxids,
            }
        )

    # Traverse children if prediction is positive
    else:
        for child in overlap_manager.tree.successors(node):
            if overlap_manager.tree.out_degree(child) > 0:  # internal node
                traversal_with_prediction(
                    overlap_manager, child, modeller, stats_matrix, tax_df, tax_level=tax_level, results=results
                )
            else:
                best_taxid_match = stats_matrix[stats_matrix.index == child]["best_match_taxid"].tolist()
                best_taxid_match = best_taxid_match[0] if len(best_taxid_match) > 0 else None
                results.append(
                    {
                        "node": child,
                        "n_leaves": 1,
                        "leaves": [child],
                        "best_taxid_match": best_taxid_match,
                        "node_precision": 1.0,
                        "node_taxids": [best_taxid_match] if best_taxid_match is not None else [],
                    }
                )

    return results


def predict_data_set_clades_composition(
    data_set_name,
    m_stats_stats_matrix,
    overlap_manager: OverlapManager,
    modeller: BaseCompositionModeller,
    input_taxa: pd.DataFrame,
    tax_level: str = "order",
):

    results = []
    for root in overlap_manager.root_nodes:
        root_results = traversal_with_prediction(
            overlap_manager, root, modeller, m_stats_stats_matrix, input_taxa, tax_level=tax_level, results=[]
        )
        results.extend(root_results)

    if len(results) == 0:
        return pd.DataFrame(
            columns=["data_set", "node", "n_leaves", "leaves", "best_taxid_match", "node_precision", "node_taxids"]
        )

    results_df = pd.DataFrame(results)

    results_df.insert(0, "data_set", data_set_name)
    return results_df


def predict_data_set_clades_fixed(
    data_set_name, m_stats_stats_matrix, overlap_manager: OverlapManager, min_dist_threshold: float = 0.6
):

    results = []
    for root in overlap_manager.root_nodes:
        root_results = traversal_with_clustering_fixed(
            overlap_manager, root, m_stats_stats_matrix, min_dist_threshold=min_dist_threshold, results=[]
        )
        results.extend(root_results)

    if len(results) == 0:
        return pd.DataFrame(
            columns=["data_set", "node", "n_leaves", "leaves", "best_taxid_match", "node_precision", "node_taxids"]
        )

    results_df = pd.DataFrame(results)

    results_df.insert(0, "data_set", data_set_name)
    return results_df


def calculate_clade_precision(result_df, input_summary=None):
    """
    Calculate the overall precision of the predicted clades.

    When input_summary is provided, computes intersection-based precision:
    |predicted_taxids ∩ input_taxids| / |predicted_taxids|.
    Otherwise falls back to legacy behavior (fraction of predictions with a match).
    """

    n_predicted = result_df.drop_duplicates(subset=["node", "best_taxid_match"])

    if input_summary is not None:
        input_taxids = set(input_summary["taxid"].unique())
        predicted_taxids = set(n_predicted["best_taxid_match"].dropna().astype(int).unique())
        if len(predicted_taxids) == 0:
            n_nodes = len(result_df.drop_duplicates(subset=["node"]))
            logger.debug(
                f"calculate_clade_precision: precision=0 (empty predicted set), "
                f"n_nodes={n_nodes}, n_input_taxids={len(input_taxids)}"
            )
            return 0.0
        intersection = predicted_taxids & input_taxids
        prec = len(intersection) / len(predicted_taxids)
        if prec == 0.0:
            sample = list(predicted_taxids)[:5]
            logger.debug(
                f"calculate_clade_precision: precision=0 (no intersection), "
                f"n_predicted_taxids={len(predicted_taxids)}, "
                f"n_input_taxids={len(input_taxids)}, "
                f"sample_predicted={sample}"
            )
        return prec

    overall_precision = (
        n_predicted.dropna(subset=["best_taxid_match"]).shape[0] / n_predicted.shape[0]
        if n_predicted.shape[0] > 0
        else 0.0
    )
    return overall_precision


def get_spurious_composition(input_df, m_stats_matrix: pd.DataFrame, tax_df: pd.DataFrame, tax_level: str = "order"):
    compositions = []

    for _, row in input_df.iterrows():
        taxid = row["taxid"]
        tax = row[tax_level]
        trash_subset = m_stats_matrix[m_stats_matrix["is_trash"] == True]
        trash_subset = trash_subset[trash_subset["cross_hit_match"] == taxid]
        subset_composition = get_subset_composition(trash_subset, tax_df, tax_level=tax_level).set_index("tax_level").T

        subset_composition = subset_composition.reset_index(drop=True)
        subset_composition.loc[:, "tax_level"] = tax

        subset_composition.insert(0, "taxid", taxid)

        compositions.append(subset_composition)
    return pd.concat(compositions, ignore_index=True, axis=0)


def get_cross_hit_composition(input_df, m_stats_matrix: pd.DataFrame, tax_df: pd.DataFrame, tax_level: str = "order"):
    compositions = []

    for _, row in input_df.iterrows():
        taxid = row["taxid"]
        tax = row[tax_level]
        crosshit_subset = m_stats_matrix[m_stats_matrix["is_crosshit"] == True]
        crosshit_subset = crosshit_subset[crosshit_subset["cross_hit_match"] == taxid]
        subset_composition = (
            get_subset_composition(crosshit_subset, tax_df, tax_level=tax_level).set_index("tax_level").T
        )

        subset_composition = subset_composition.reset_index(drop=True)
        subset_composition.loc[:, "tax_level"] = tax

        subset_composition.insert(0, "taxid", taxid)

        compositions.append(subset_composition)
    return pd.concat(compositions, ignore_index=True, axis=0)
