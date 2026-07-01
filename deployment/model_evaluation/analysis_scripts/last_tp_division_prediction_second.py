"""
Predict the last true-positive division: where recall crosses threshold τ.

Compares 4 decision approaches × 13 model variants:
  Approaches:
    - CLF:         first division where P(recall ≥ τ) ≥ X (probabilistic)
    - Likelihood:  first division where pdf(τ; μ, σ) ≥ threshold (probabilistic)
    - Direct:      predict τ-crossing fraction directly from predicted mean (deterministic)
    - Direct+ISO:  Direct + isotonic regression calibration
  Models:
    - Baseline GP (6 kernels: RBF, Matern-0.5, Matern-2.5, RQ, RBF+RQ, MaternxMatern)
    - Hurdle (CLF + Likelihood)
    - RF Direct
    - XGBoost Direct
    - NGBoost Direct
    - DirectXGB
    - DirectRF
    - GP-DirectXGB (stacked meta-model)
    - GPCLFThreshold (production baseline)

With configurable asymmetric loss (over- vs. under-prediction penalty).
"""

import argparse
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

import sys
import time
from copy import deepcopy
from itertools import product
from pathlib import Path

import matplotlib
from ngboost import NGBRegressor
from ngboost.distns import Normal
from scipy.stats import norm
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, Matern, RationalQuadratic, WhiteKernel
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputRegressor
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")
plt.rcParams["figure.dpi"] = 120


# ── CLI args ──────────────────────────────────────────────────
def get_args():
    parser = argparse.ArgumentParser(
        description="Predict the last true-positive division: where recall crosses threshold tau.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="last_tp_division_outputs",
        help="Output directory for results (default: last_tp_division_outputs)",
    )
    parser.add_argument(
        "--input_cache",
        type=str,
        default="/home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/"
        "study/analysis/virus/filter_gp_clf_family/models/"
        "cache/recall_results_cache.parquet",
        help="Path to recall results cache parquet file",
    )
    return parser.parse_args()


_args = get_args()

# ── globals ──────────────────────────────────────────────────
random_state = 42
EPS = 1e-6
SPIKE_THRESHOLD = 0.999
MIN_NONSPIKE_SAMPLES = 10
N_DIVISIONS = 20
DIVISION_FRACTIONS = np.arange(1, N_DIVISIONS + 1) / N_DIVISIONS  # [0.05, 0.10, ..., 1.0]
BASELINE_TAU = 1.0  # reference for comparison plots
TARGET_TAUS = [0.85, 0.90, 0.95, 0.98, 1.0]  # sweep these target recall values

# ── Stacked GP→XGBoost meta-model config ──
GP_SUBSET_DIVISIONS = [0, 4, 9, 14]  # fractions [0.05, 0.25, 0.50, 0.75]
GP_KERNEL = "Matern-2.5"

# ── Active divisions (last 4 are recall=1.0 constant) ──
N_ACTIVE_DIVISIONS = 16

# ── Kernel variants for the Baseline GP ──
KERNELS = {
    "RBF": ConstantKernel(1.0) * RBF(length_scale=1.0),
    "Matern-0.5": ConstantKernel(1.0) * Matern(length_scale=1.0, nu=0.5),
    "Matern-2.5": ConstantKernel(1.0) * Matern(length_scale=1.0, nu=2.5),
    "RationalQuadratic": ConstantKernel(1.0) * RationalQuadratic(length_scale=1.0, alpha=1.0),
    "RBF+RQ": ConstantKernel(1.0) * (RBF(length_scale=1.0) + RationalQuadratic(length_scale=1.0, alpha=1.0)),
    "MaternxMatern": ConstantKernel(1.0) * Matern(length_scale=1.0, nu=1.5) * Matern(length_scale=1.0, nu=1.5),
}

SAVE_DIR = Path(_args.output_dir)


def model_file_safe(name):
    return name.replace(" ", "_").replace("(", "").replace(")", "")


SAVE_DIR.mkdir(parents=True, exist_ok=True)

# ═══════════════════════════════════════════════════════════════
# 1. DATA LOADING
# ═══════════════════════════════════════════════════════════════
print("=" * 60)
print("1. DATA LOADING")
print("=" * 60)

trainning_data_cache = Path(_args.input_cache)

training_data = pd.read_parquet(trainning_data_cache)
print(f"Training data loaded with shape: {training_data.shape}")

y_columns = [x for x in training_data.columns if x.startswith("index_recall_")]
y_columns_active = y_columns[:N_ACTIVE_DIVISIONS]
y_second = "last_best_match_relindex"
features = [col for col in training_data.columns if col not in y_columns + [y_second]]
features = [x for x in features if x not in ["data_set", "unclassified"]]
stats_features = ["counts_kurtosis", "counts_skewness", "tax_diversity_shannon", "max_uniq_reads", "total_uniq_reads"]
taxonomy_features = [col for col in features if col not in stats_features]

n_features = len(features)
print(f"Features ({n_features}): {len(stats_features)} stats + {len(taxonomy_features)} taxonomy")

# Filter out samples where last_best_match_relindex == 0 or 1
# (these are degenerate / edge cases)
training_data = training_data[
    (training_data["last_best_match_relindex"] > 0.0) & (training_data["last_best_match_relindex"] < 1.0)
]

X = training_data[features].copy()
# for ftax in taxonomy_features:
#    X[ftax] = (X[ftax] > 0).astype(float)

Y = training_data[y_columns]
# filter samples where max y_columns = 0.0
filter_nullmask = Y.max(axis=1) > 0.0
X = X.loc[filter_nullmask]
Y = Y.loc[filter_nullmask]
# recall does not always reach 1.0, so we divide each row by its max to get a fraction in [0, 1]
# Y = Y.div(Y.max(axis=1), axis=0).fillna(0.0)


print(f"X shape: {X.shape}, Y shape: {Y.shape}")

X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=random_state)
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

print(f"Train: {X_train_scaled.shape[0]} samples, Test: {X_test_scaled.shape[0]} samples")


# ═══════════════════════════════════════════════════════════════
# 2. HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════
print("=" * 60)
print("2. HELPER FUNCTIONS")
print("=" * 60)


# ── Shared helpers from om_models.py (with standalone fallback) ──
try:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
    from metagenomics_utils.overlap_manager.om_models import (
        actual_recall_at_fraction as _actual_recall_at_fraction,
    )
    from metagenomics_utils.overlap_manager.om_models import (
        asymmetric_loss,
        predict_division_clf,
    )
    from metagenomics_utils.overlap_manager.om_models import (
        compute_ground_truth_division as _compute_ground_truth_division,
    )

    def compute_ground_truth_division(recall_at_divs, tau):
        return _compute_ground_truth_division(recall_at_divs, tau, N_DIVISIONS, DIVISION_FRACTIONS)

    def actual_recall_at_fraction(recall_curve, fraction):
        return _actual_recall_at_fraction(recall_curve, fraction, N_DIVISIONS, DIVISION_FRACTIONS)
except ImportError:
    # ── Fallback local definitions for standalone mode ──
    def compute_ground_truth_division(recall_at_divs, tau):
        for i in range(N_DIVISIONS):
            if recall_at_divs[i] >= tau:
                if i == 0:
                    return float(tau * DIVISION_FRACTIONS[0] / max(recall_at_divs[0], EPS))
                r_low, r_high = recall_at_divs[i - 1], recall_at_divs[i]
                f_low, f_high = DIVISION_FRACTIONS[i - 1], DIVISION_FRACTIONS[i]
                if r_high == r_low:
                    return f_high
                frac = f_low + (tau - r_low) * (f_high - f_low) / (r_high - r_low)
                return float(np.clip(frac, f_low, f_high))
        return 1.0

    def asymmetric_loss(y_true, y_pred, over_cost=1.0, under_cost=3.0):
        diff = y_pred - y_true
        loss = np.where(diff > 0, over_cost * diff**2, under_cost * diff**2)
        return float(np.mean(loss))

    def actual_recall_at_fraction(recall_curve, fraction):
        if fraction <= 0.0:
            return 0.0
        if fraction >= 1.0:
            return float(recall_curve[-1])
        div_idx = int(np.floor(fraction * N_DIVISIONS))
        div_idx = min(div_idx, N_DIVISIONS - 2)
        f_low, f_high = DIVISION_FRACTIONS[div_idx], DIVISION_FRACTIONS[div_idx + 1]
        r_low, r_high = recall_curve[div_idx], recall_curve[div_idx + 1]
        if f_high == f_low:
            return float(r_low)
        frac = (fraction - f_low) / (f_high - f_low)
        return float(r_low + frac * (r_high - r_low))

    def predict_division_clf(probs, X_cutoff, fractions):
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


def compute_ece(pred_probs, binary_outcomes, n_bins=10):
    """Expected Calibration Error."""
    bins = np.linspace(0, 1, n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    ece = 0.0
    bin_counts = []
    bin_accs = []
    for i in range(n_bins):
        mask = (pred_probs >= bins[i]) & (pred_probs < bins[i + 1])
        if i == n_bins - 1:
            mask = (pred_probs >= bins[i]) & (pred_probs <= bins[i + 1])
        if mask.sum() == 0:
            bin_counts.append(0)
            bin_accs.append(0.0)
            continue
        bin_prob = pred_probs[mask].mean()
        bin_acc = binary_outcomes[mask].mean()
        weight = mask.sum() / len(pred_probs)
        ece += weight * abs(bin_prob - bin_acc)
        bin_counts.append(mask.sum())
        bin_accs.append(bin_acc)
    return ece, bin_centers, np.array(bin_counts), np.array(bin_accs)


# ── Tau-sweep helpers ──────────────────────────────────────────
def compute_ground_truth_for_tau(y_test_array, tau, n_test):
    """Ground-truth fraction where recall >= tau for each test sample."""
    return np.array([compute_ground_truth_division(y_test_array[:, s], tau) for s in range(n_test)])


def precompute_prob_lik(posterior, tau, model_keys):
    """Precompute P(recall > tau) and scaled likelihood for each model at given tau."""
    pre = {}
    for mk in model_keys:
        mu = posterior[mk]["mu"]
        std = posterior[mk]["std"]
        ps = posterior[mk]["p_spike"]
        if mk == "Hurdle":
            p_cond = 1.0 - norm.cdf(tau, loc=mu, scale=std)
            prob = ps + (1.0 - ps) * p_cond
        else:
            prob = 1.0 - norm.cdf(tau, loc=mu, scale=std)
        z = (tau - mu) / std
        lik = np.exp(-0.5 * z**2)
        pre[mk] = {"prob": prob, "lik": lik}
    return pre


def run_grid_search_for_tau(
    precomputed, deterministic_preds, ground_truth, model_keys, n_test, tau, X_vals, lik_tau_vals, asym_ratios
):
    """Grid search over CLF + Likelihood approaches for a given tau."""
    rows = []
    for mk in model_keys:
        prob = precomputed[mk]["prob"]
        lik = precomputed[mk]["lik"]
        for over_cost, under_cost in asym_ratios:
            for X_val in X_vals:
                yp = predict_division_clf(prob, X_val, DIVISION_FRACTIONS)
                loss = asymmetric_loss(ground_truth, yp, over_cost, under_cost)
                rows.append(
                    {
                        "model": mk,
                        "approach": "CLF",
                        "param": X_val,
                        "over_cost": over_cost,
                        "under_cost": under_cost,
                        "loss": loss,
                        "tau": tau,
                    }
                )
            for lik_tau in lik_tau_vals:
                yp = predict_division_lik(lik, lik_tau, DIVISION_FRACTIONS)
                loss = asymmetric_loss(ground_truth, yp, over_cost, under_cost)
                rows.append(
                    {
                        "model": mk,
                        "approach": "Likelihood",
                        "param": lik_tau,
                        "over_cost": over_cost,
                        "under_cost": under_cost,
                        "loss": loss,
                        "tau": tau,
                    }
                )
    for mk, mu in deterministic_preds.items():
        yp = np.array([compute_ground_truth_division(mu[:, s], tau) for s in range(n_test)])
        for over_cost, under_cost in asym_ratios:
            loss = asymmetric_loss(ground_truth, yp, over_cost, under_cost)
            rows.append(
                {
                    "model": mk,
                    "approach": "Direct",
                    "param": 0.0,
                    "over_cost": over_cost,
                    "under_cost": under_cost,
                    "loss": loss,
                    "tau": tau,
                }
            )
    return pd.DataFrame(rows)


def find_optimal_params(results_df):
    """Best param per (model, approach, over_cost, under_cost, tau)."""
    summary = []
    for keys, group in results_df.groupby(["model", "approach", "over_cost", "under_cost", "tau"]):
        best = group.loc[group["loss"].idxmin()]
        summary.append(best)
    return pd.DataFrame(summary)


def evaluate_actual_recall_at_prediction(y_pred, y_test_array, n_test):
    """Actual recall achieved at each predicted fraction."""
    return np.array([actual_recall_at_fraction(y_test_array[:, s], y_pred[s]) for s in range(n_test)])


# ── Filter degenerate samples (never reach τ=1.0 before last division) ──
print("\nFiltering samples that never reach τ=1.0 before the last division...")
for split_name, y_df, X_df, X_scaled in [
    ("Train", y_train, X_train, X_train_scaled),
    ("Test", y_test, X_test, X_test_scaled),
]:
    y_arr = y_df.values if hasattr(y_df, "values") else np.asarray(y_df)
    gt_arr = np.array([compute_ground_truth_division(y_arr[s], BASELINE_TAU) for s in range(len(y_arr))])
    valid = gt_arr < 1.0
    n_before = len(y_df)
    n_after = int(valid.sum())
    n_removed = n_before - n_after
    if n_removed > 0:
        print(f"  {split_name}: removed {n_removed}/{n_before} samples (keep {n_after})")
        idx = np.where(valid)[0]
        if split_name == "Train":
            y_train = y_df.iloc[idx]
            X_train = X_df.iloc[idx]
            X_train_scaled = X_scaled[idx]
        else:
            y_test = y_df.iloc[idx]
            X_test = X_df.iloc[idx]
            X_test_scaled = X_scaled[idx]
    else:
        print(f"  {split_name}: all {n_before} samples valid")


# ═══════════════════════════════════════════════════════════════
# 3. HURDLE MODEL CLASS
# ═══════════════════════════════════════════════════════════════
print("=" * 60)
print("3. HURDLE MODEL CLASS")
print("=" * 60)


class HurdleRecallModel:
    """Two-stage hurdle: spike classifier (LR) + conditional GP (Matern)."""

    def __init__(self, percentile_label=""):
        self.clf = LogisticRegression(C=1.0, max_iter=1000, class_weight="balanced", random_state=42)
        self.clf_fitted = True
        self.fixed_p_spike = None
        self.gp = None
        self.cond_mean = 0.0
        self.cond_std = 0.01
        self.has_gp = False
        self.label = percentile_label

    def fit(self, X, y):
        y_spike = (y >= SPIKE_THRESHOLD).astype(int)
        n_classes = len(np.unique(y_spike))
        if n_classes >= 2:
            self.clf.fit(X, y_spike)
            self.clf_fitted = True
            self.fixed_p_spike = None
        else:
            self.clf_fitted = False
            self.fixed_p_spike = float(y_spike.mean())

        mask = y < SPIKE_THRESHOLD
        n_nonspike = mask.sum()
        if n_nonspike >= MIN_NONSPIKE_SAMPLES:
            kernel = ConstantKernel(1.0) * Matern(length_scale=1.0, nu=1.5) + WhiteKernel(noise_level=0.1)
            self.gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=3, random_state=42)
            self.gp.fit(X[mask], y[mask])
            self.has_gp = True
        elif n_nonspike > 0:
            self.cond_mean = np.mean(y[mask])
            self.cond_std = max(np.std(y[mask]), 0.01)

    def _get_p_spike(self, X):
        if self.clf_fitted:
            return self.clf.predict_proba(X)[:, 1]
        else:
            return np.full(len(X), self.fixed_p_spike)

    def predict_mean(self, X):
        p_spike = self._get_p_spike(X)
        if self.has_gp:
            mu_c = self.gp.predict(X)
        else:
            mu_c = np.full(len(X), self.cond_mean)
        return p_spike * 1.0 + (1.0 - p_spike) * mu_c

    def predict_attainment(self, X, thr):
        """Return P(recall >= thr) for each test point."""
        p_spike = self._get_p_spike(X)
        if self.has_gp:
            mu_c, std_c = self.gp.predict(X, return_std=True)
        else:
            mu_c = np.full(len(X), self.cond_mean)
            std_c = np.full(len(X), self.cond_std)
        std_c = np.maximum(std_c, 1e-6)
        p_cond = 1.0 - norm.cdf(thr, loc=mu_c, scale=std_c)
        return p_spike + (1.0 - p_spike) * p_cond

    def predict_pdf_scaled(self, X, thr):
        """
        Return standardized likelihood pdf_scaled_at_thr
        = φ(thr; μ_c, σ_c) / max(φ) = exp(-0.5 * ((thr - μ_c) / σ_c)²)
        Uses conditional GP (non-spike part) only.
        """
        if self.has_gp:
            mu_c, std_c = self.gp.predict(X, return_std=True)
        else:
            mu_c = np.full(len(X), self.cond_mean)
            std_c = np.full(len(X), self.cond_std)
        std_c = np.maximum(std_c, 1e-6)
        z = (thr - mu_c) / std_c
        return np.exp(-0.5 * z**2)


# ═══════════════════════════════════════════════════════════════
# 4. TRAIN MODELS ON ALL 20 DIVISIONS
# ═══════════════════════════════════════════════════════════════
print("=" * 60)
print("4. TRAINING MODELS ON ALL 20 DIVISIONS")
print("=" * 60)

# 4a. Baseline GP variants (6 kernels)
base_gp_models = {}
for kernel_name, base_kernel in KERNELS.items():
    kernel = base_kernel + WhiteKernel(noise_level=0.1)
    print(f"\n--- Baseline GP ({kernel_name}) ---")
    models = {}
    t0 = time.time()
    for col in y_columns_active:
        y_tr = y_train[col].values.copy()
        gp = GaussianProcessRegressor(
            kernel=deepcopy(kernel), normalize_y=True, n_restarts_optimizer=3, random_state=42
        )
        gp.fit(X_train_scaled, y_tr)
        models[col] = gp
    t = time.time() - t0
    print(f"  Trained in {t:.1f}s")
    base_gp_models[kernel_name] = models

# ── Extract GP-Matern features for stacking ──
print("\n--- Extracting GP-Matern features for stacking ---")
gp_feat_train_list, gp_feat_test_list = [], []
for d in GP_SUBSET_DIVISIONS:
    col = y_columns_active[d]
    mu_tr, std_tr = base_gp_models[GP_KERNEL][col].predict(X_train_scaled, return_std=True)
    mu_te, std_te = base_gp_models[GP_KERNEL][col].predict(X_test_scaled, return_std=True)
    gp_feat_train_list.extend([mu_tr.flatten(), np.maximum(std_tr.flatten(), 1e-6)])
    gp_feat_test_list.extend([mu_te.flatten(), np.maximum(std_te.flatten(), 1e-6)])
gp_feat_train = np.column_stack(gp_feat_train_list)
gp_feat_test = np.column_stack(gp_feat_test_list)
print(f"  GP features shape: train {gp_feat_train.shape}, test {gp_feat_test.shape}")

# 4b. Hurdle model (LogisticReg + Matern GP)
print("\n--- Hurdle Model ---")
hurdle_models = {}
t0 = time.time()
for col in y_columns_active:
    y_tr = y_train[col].values.copy()
    hrm = HurdleRecallModel(col)
    hrm.fit(X_train_scaled, y_tr)
    hurdle_models[col] = hrm
hurdle_time = time.time() - t0
print(f"  Hurdle models trained in {hurdle_time:.1f}s")

# 4c. Random Forest (multi-output, direct recall prediction)
print("\n--- Random Forest ---")
rf_model = MultiOutputRegressor(RandomForestRegressor(n_estimators=100, random_state=42))
t0 = time.time()
rf_model.fit(X_train_scaled, y_train[y_columns_active].values)
print(f"  RF trained in {time.time() - t0:.1f}s")

# 4d. XGBoost (multi-output, direct recall prediction)
print("\n--- XGBoost ---")
xgb_model = MultiOutputRegressor(XGBRegressor(n_estimators=100, random_state=42))
t0 = time.time()
xgb_model.fit(X_train_scaled, y_train[y_columns_active].values)
print(f"  XGBoost trained in {time.time() - t0:.1f}s")

# 4e. NGBoost (probabilistic gradient boosting, one model per division)
#     Grid search over hyperparams via CRPS, then train on non-constant divisions
print("\n--- NGBoost: grid search + training ---")

y_std_cols = y_train[y_columns_active].std(axis=0).values
trainable = np.where(y_std_cols > 1e-4)[0]
constant = np.where(y_std_cols <= 1e-4)[0]
print(f"  Trainable: {len(trainable)} divisions, Constant (skip): {len(constant)} divisions")

X_gs_tr, X_gs_val, y_gs_tr, y_gs_val = train_test_split(
    X_train_scaled,
    y_train[y_columns_active],
    test_size=0.2,
    random_state=42,
)

rep_divs = [d for d in [0, 4, 8, 12, 15] if d in trainable]

best_crps = np.inf
best_params_ngb = None

for n_est, lr, mb in product([200, 500], [0.01, 0.05], [0.5, 1.0]):
    scores = []
    t0 = time.time()
    for d in rep_divs:
        col = y_columns_active[d]
        model = NGBRegressor(
            Dist=Normal,
            n_estimators=n_est,
            learning_rate=lr,
            minibatch_frac=mb,
            verbose=False,
            random_state=42,
        )
        y_vals = y_gs_tr[col].values
        finite = np.isfinite(y_vals)
        model.fit(X_gs_tr[finite], y_vals[finite])
        dist = model.pred_dist(X_gs_val)
        mu = dist.params["loc"]
        sigma = np.maximum(dist.params["scale"], 1e-8)
        z = (y_gs_val[col].values - mu) / sigma
        crps = sigma * (z * (2 * norm.cdf(z) - 1) + 2 * norm.pdf(z) - 1.0 / np.sqrt(np.pi))
        scores.append(np.mean(crps))
    avg_crps = np.mean(scores)
    elapsed = time.time() - t0
    print(f"  n_est={n_est:3d}  lr={lr:.3f}  mb={mb:.1f}  CRPS={avg_crps:.5f}  [{elapsed:.0f}s]")
    if avg_crps < best_crps:
        best_crps = avg_crps
        best_params_ngb = {"n_estimators": n_est, "learning_rate": lr, "minibatch_frac": mb}

print(f"Best NGBoost params: {best_params_ngb} (CRPS={best_crps:.5f})")

print("--- Training NGBoost models with best params ---")
ngb_models = {}
t0 = time.time()
for d in trainable:
    col = y_columns_active[d]
    model = NGBRegressor(
        Dist=Normal,
        **best_params_ngb,
        verbose=False,
        random_state=42,
    )
    y_vals = y_train[col].values
    finite = np.isfinite(y_vals)
    model.fit(X_train_scaled[finite], y_vals[finite])
    ngb_models[col] = model
for d in constant:
    col = y_columns_active[d]
    ngb_models[col] = None
print(f"  NGBoost trained in {time.time() - t0:.1f}s")

# ── Extract posterior predictions on test set ──────────────
print("\nExtracting posterior predictions on test set...")

# Structure:  {model_key: {mu: (n_div, n_test), std: (n_div, n_test), p_spike: (n_div, n_test)}}
# 6 kernel variants + Hurdle = 7 model keys
posterior = {}
for kernel_name in KERNELS:
    posterior[f"Baseline GP {kernel_name}"] = {
        "mu": [],
        "std": [],
        "p_spike": [np.zeros(len(X_test_scaled))] * N_DIVISIONS,
    }
posterior["Hurdle"] = {"mu": [], "std": [], "p_spike": []}

for d in range(N_ACTIVE_DIVISIONS):
    col = y_columns_active[d]
    y_te = y_test[col].values.copy()

    # Baseline GP variants
    for kernel_name in KERNELS:
        key = f"Baseline GP {kernel_name}"
        mu_gp, std_gp = base_gp_models[kernel_name][col].predict(X_test_scaled, return_std=True)
        posterior[key]["mu"].append(mu_gp.flatten())
        posterior[key]["std"].append(np.maximum(std_gp.flatten(), 1e-6))

    # Hurdle
    hrm = hurdle_models[col]
    p_spike = hrm._get_p_spike(X_test_scaled).flatten()
    posterior["Hurdle"]["p_spike"].append(p_spike)
    if hrm.has_gp:
        mu_c, std_c = hrm.gp.predict(X_test_scaled, return_std=True)
        posterior["Hurdle"]["mu"].append(mu_c.flatten())
        posterior["Hurdle"]["std"].append(np.maximum(std_c.flatten(), 1e-6))
    else:
        posterior["Hurdle"]["mu"].append(np.full(len(X_test_scaled), hrm.cond_mean))
        posterior["Hurdle"]["std"].append(np.full(len(X_test_scaled), max(hrm.cond_std, 1e-6)))

# Pad inactive divisions (recall=1.0 constant, no uncertainty)
n_remain = N_DIVISIONS - N_ACTIVE_DIVISIONS
if n_remain > 0:
    for model_key in posterior:
        posterior[model_key]["mu"].extend([np.ones(len(X_test_scaled))] * n_remain)
        posterior[model_key]["std"].extend([np.full(len(X_test_scaled), 1e-6)] * n_remain)
    posterior["Hurdle"]["p_spike"].extend([np.zeros(len(X_test_scaled))] * n_remain)

# Convert lists to arrays: (n_divisions, n_test)
for model_key in posterior:
    for key in ["mu", "std"]:
        posterior[model_key][key] = np.array(posterior[model_key][key])
    posterior[model_key]["p_spike"] = np.array(posterior[model_key]["p_spike"])

n_test = len(X_test_scaled)
y_test_array = y_test.values.T  # (n_divisions, n_test)

print(f"Posterior arrays: ({N_DIVISIONS}, {n_test}) per model ({len(posterior)} models)")

# ── Extract deterministic predictions (RF, XGBoost, NGBoost) ──
print("\nExtracting deterministic predictions...")


def _pad_recall(preds):
    n_remain = N_DIVISIONS - N_ACTIVE_DIVISIONS
    return np.vstack([preds, np.ones((n_remain, preds.shape[1]))])


rf_preds = _pad_recall(rf_model.predict(X_test_scaled).T)
xgb_preds = _pad_recall(xgb_model.predict(X_test_scaled).T)
ngb_preds = _pad_recall(
    np.array(
        [
            ngb_models[col].predict(X_test_scaled) if ngb_models[col] is not None else np.full(len(X_test_scaled), 1.0)
            for col in y_columns_active
        ]
    )
)
deterministic_preds = {"RF": rf_preds, "XGBoost": xgb_preds, "NGBoost": ngb_preds}

# ═══════════════════════════════════════════════════════════════
# 4b. GPCLFTHRESHOLD FIT (fit once, evaluate per tau later)
# ═══════════════════════════════════════════════════════════════
print("=" * 60)
print("4b. GPCLFTHRESHOLD FIT (production implementation)")
print("=" * 60)

from metagenomics_utils.overlap_manager.om_models import GPCLFThreshold

gpclf_model = GPCLFThreshold(
    n_divisions=N_ACTIVE_DIVISIONS,
    optimize=False,
    tau=0.95,
    filter_degenerate=False,
    binarize_taxonomy=False,
    val_split=0.0,
    random_state=random_state,
)
gpclf_model.fit(X_train, y_train[y_columns_active], target_names=list(y_columns_active))
gpclf_y_pred = gpclf_model.predict(X_test)
gpclf_actual_recalls = np.array([actual_recall_at_fraction(y_test_array[:, s], gpclf_y_pred[s]) for s in range(n_test)])

# Diagnostics
print(
    f"  y_pred: min={gpclf_y_pred.min():.3f}, max={gpclf_y_pred.max():.3f}, "
    f"mean={gpclf_y_pred.mean():.3f}, std={gpclf_y_pred.std():.3f}"
)
raw = gpclf_model.predict_raw(X_test)
print("  Raw recall preds (mean \u00b1 std across test samples per division):")
for d in range(N_ACTIVE_DIVISIONS):
    print(f"    div {d:2d} ({DIVISION_FRACTIONS[d]:.2f}): {raw[:, d].mean():.4f} \u00b1 {raw[:, d].std():.4f}")
tau_gpclf = gpclf_model.tau
for div_idx, label in [(0, "first"), (-1, "last")]:
    scale = max(raw[:, div_idx].std(), 0.01)
    p = 1.0 - norm.cdf(tau_gpclf, loc=raw[:, div_idx], scale=scale)
    print(f"  P(recall>{tau_gpclf:.2f}) @ {label} div: mean={p.mean():.3f} \u00b1 {p.std():.3f}")
print()

# ── Grid search constants (shared across all tau) ──
ASYMMETRIC_RATIOS = [(1, 1), (1.5, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1), (8, 1), (10, 1)]
ASYMMETRIC_RATIOS = [(1, r) for r in [1.5, 2, 3, 4, 5, 6, 8, 10]]

X_VALUES = np.linspace(0.025, 0.975, 39)
LIK_TAU_VALUES = np.linspace(0.025, 0.975, 39)


def predict_division_lik(lik_values, lik_cutoff, fractions):
    """Find first division where lik >= lik_cutoff, else argmax."""
    n = lik_values.shape[1]
    predictions = np.full(n, np.nan)
    for s in range(n):
        idxs = np.where(lik_values[:, s] >= lik_cutoff)[0]
        if len(idxs) == 0:
            i = np.argmax(lik_values[:, s])
            predictions[s] = fractions[i]
            continue
        i = idxs[0]
        if i == 0:
            predictions[s] = fractions[0] * lik_cutoff / max(lik_values[0, s], EPS)
        else:
            lik_low, lik_high = lik_values[i - 1, s], lik_values[i, s]
            f_low, f_high = fractions[i - 1], fractions[i]
            if lik_high == lik_low:
                predictions[s] = f_high
            else:
                frac = f_low + (lik_cutoff - lik_low) * (f_high - f_low) / (lik_high - lik_low)
                predictions[s] = float(np.clip(frac, f_low, f_high))
    return np.clip(predictions, 0.0, 1.0)


# ═══════════════════════════════════════════════════════════════
# 5–11. TAU LOOP
# ═══════════════════════════════════════════════════════════════
model_keys = list(posterior.keys())
TARGET_OVER = 1
TARGET_UNDER = 3
all_results_dfs = []
all_summary_dfs = []
perf_rows = []

for tau in TARGET_TAUS:
    print("\n" + "=" * 60)
    print(f"  TARGET τ = {tau:.2f}")
    print("=" * 60)

    # ══════════════════════════════════════════════════════════
    # 5. GROUND TRUTH
    # ══════════════════════════════════════════════════════════
    print(f"\n5. COMPUTING GROUND TRUTH (τ = {tau})")
    ground_truth = compute_ground_truth_for_tau(y_test_array, tau, n_test)

    print(f"\nGround truth (fraction where recall ≥ {tau}):")
    gt_series = pd.Series(ground_truth)
    print(gt_series.describe())
    never_reach = (ground_truth >= 1.0).mean()
    print(f"  Samples reaching τ only at last division (gt >= 1.0): {never_reach:.1%}")

    # ══════════════════════════════════════════════════════════
    # 6. GRID SEARCH
    # ══════════════════════════════════════════════════════════
    print(f"\n6. GRID SEARCH (τ = {tau})")
    precomputed = precompute_prob_lik(posterior, tau, model_keys)

    results_df = run_grid_search_for_tau(
        precomputed,
        deterministic_preds,
        ground_truth,
        model_keys,
        n_test,
        tau,
        X_VALUES,
        LIK_TAU_VALUES,
        ASYMMETRIC_RATIOS,
    )
    all_results_dfs.append(results_df)
    print(f"  Grid search complete. Results shape: {results_df.shape}")

    # ── 6b. GPCLFThreshold evaluation ──
    gpclf_below = (gpclf_actual_recalls < tau).mean()
    gpclf_loss_tau = asymmetric_loss(ground_truth, gpclf_y_pred, over_cost=1, under_cost=3)
    print(f"\nGPCLFThreshold (\u03c4={gpclf_model.tau:.3f}, X={gpclf_model.x_thresh:.3f})")
    print(
        f"  Actual recall at cutoff: mean={gpclf_actual_recalls.mean():.3f}, "
        f"median={np.median(gpclf_actual_recalls):.3f}"
    )
    print(f"  % below \u03c4={tau}: {gpclf_below:.1%}")

    # ══════════════════════════════════════════════════════════
    # 7. OPTIMAL PARAMETERS
    # ══════════════════════════════════════════════════════════
    print(f"\n7. OPTIMAL PARAMETERS (τ = {tau})")
    summary_df = find_optimal_params(results_df)
    all_summary_dfs.append(summary_df)

    print("\n--- Best loss per (model, approach, ratio) ---")
    for (m, a, oc, uc), g in summary_df.groupby(["model", "approach", "over_cost", "under_cost"]):
        best = g.loc[g["loss"].idxmin()]
        print(f"  {m:12s} | {a:11s} | over:{oc} under:{uc} | param={best['param']:.2f} | loss={best['loss']:.4f}")

    # ══════════════════════════════════════════════════════════
    # 8a. LOSS CURVES (only for BASELINE_TAU)
    # ══════════════════════════════════════════════════════════
    if tau == BASELINE_TAU:
        print("\n--- 8a. Loss vs param curves ---")
        for mk in model_keys:
            fig, ax = plt.subplots(figsize=(10, 6))
            for approach in ["CLF", "Likelihood"]:
                for oc, uc in [(1, 1), (3, 1), (10, 1)]:
                    sub = summary_df[
                        (summary_df["model"] == mk)
                        & (summary_df["approach"] == approach)
                        & (summary_df["over_cost"] == oc)
                        & (summary_df["under_cost"] == uc)
                    ]
                    if sub.empty:
                        continue
                    best = sub.loc[sub["loss"].idxmin()]
                    ax.plot(best["param"], best["loss"], "o", label=f"{approach} {oc}:{uc}")
            ax.set_xlabel("Parameter value (X or \u03c4_lik)")
            ax.set_ylabel("Best loss")
            ax.set_title(f"{mk} \u2014 optimal loss by approach and ratio (\u03c4={tau})")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            fname = f"loss_curves_{model_file_safe(mk)}_tau{tau:.2f}.png"
            fig.savefig(SAVE_DIR / fname, dpi=120, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved {fname}")

    # ══════════════════════════════════════════════════════════
    # 8b. HEATMAP (per tau)
    # ══════════════════════════════════════════════════════════
    print("--- 8b. Best loss heatmap ---")
    heat_data = summary_df.groupby(["model", "approach", "over_cost", "under_cost"])["loss"].min().reset_index()
    heat_pivot = heat_data.pivot_table(
        index="model", columns=["approach", "over_cost", "under_cost"], values="loss", aggfunc="min"
    )
    n_models = len(heat_pivot)
    fig_h = max(6, n_models * 0.45)
    fig, ax = plt.subplots(figsize=(16, fig_h))
    sns.heatmap(
        heat_pivot,
        annot=True,
        fmt=".4f",
        cmap="YlOrRd_r",
        ax=ax,
        cbar_kws={"label": "Best loss"},
        linewidths=0.5,
        linecolor="gray",
    )
    ax.set_title(f"Best asymmetric loss by model, approach, and ratio (\u03c4={tau})", fontsize=12)
    ax.set_ylabel("")
    plt.tight_layout()
    fname = f"best_loss_comparison_tau{tau:.2f}.png"
    fig.savefig(SAVE_DIR / fname, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {fname}")

    # ══════════════════════════════════════════════════════════
    # 9. SENSITIVITY (only for BASELINE_TAU)
    # ══════════════════════════════════════════════════════════
    if tau == BASELINE_TAU:
        print("--- 9. Sensitivity analysis ---")
        fig, ax = plt.subplots(figsize=(14, 7))
        for (m, a), g in summary_df.groupby(["model", "approach"]):
            ratio_data = g.groupby(["over_cost", "under_cost"])["loss"].min().reset_index()
            ratio_data["ratio"] = ratio_data.apply(lambda r: f"{r['over_cost']}:{r['under_cost']}", axis=1)
            ax.plot(ratio_data["ratio"], ratio_data["loss"], "o-", lw=1.5, label=f"{m} \u2014 {a}")
        ax.set_xlabel("Asymmetry ratio (over:under)")
        ax.set_ylabel("Best loss (min over param)")
        ax.set_title("Sensitivity to asymmetry ratio")
        ax.legend(fontsize=7, ncol=2)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        fig.savefig(SAVE_DIR / "sensitivity_ratio.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        print("  Saved sensitivity_ratio.png")

    # ══════════════════════════════════════════════════════════
    # 10. SUMMARY TABLE (per tau)
    # ══════════════════════════════════════════════════════════
    print("--- 10. Summary table ---")
    tau_summary = []
    for (m, a, oc, uc), g in summary_df.groupby(["model", "approach", "over_cost", "under_cost"]):
        best = g.loc[g["loss"].idxmin()]
        tau_summary.append(
            {
                "model": m,
                "approach": a,
                "ratio": f"{oc}:{uc}",
                "optimal_param": round(best["param"], 2),
                "best_loss": round(best["loss"], 4),
                "tau": round(tau, 2),
            }
        )
    tau_final = pd.DataFrame(tau_summary)
    print(tau_final.to_string(index=False))
    tau_final.to_csv(SAVE_DIR / f"summary_results_tau{tau:.2f}.csv", index=False)
    results_df.to_parquet(SAVE_DIR / f"full_grid_results_tau{tau:.2f}.parquet")
    print(f"  Results saved for \u03c4={tau:.2f}")

    # ══════════════════════════════════════════════════════════
    # COLLECT PERFORMANCE (for cross-tau comparison)
    # ══════════════════════════════════════════════════════════
    for mk in model_keys + list(deterministic_preds.keys()):
        for approach in ["CLF", "Likelihood", "Direct"]:
            sub = summary_df[
                (summary_df["model"] == mk)
                & (summary_df["approach"] == approach)
                & (summary_df["over_cost"] == TARGET_OVER)
                & (summary_df["under_cost"] == TARGET_UNDER)
            ]
            if sub.empty:
                continue
            best = sub.loc[sub["loss"].idxmin()]
            param_opt = best["param"]
            best_loss = best["loss"]

            if approach == "Direct":
                mu = deterministic_preds[mk]
                y_pred = np.array([compute_ground_truth_division(mu[:, s], tau) for s in range(n_test)])
            else:
                prob = precomputed[mk]["prob"]
                lik = precomputed[mk]["lik"]
                if approach == "CLF":
                    y_pred = predict_division_clf(prob, param_opt, DIVISION_FRACTIONS)
                else:
                    y_pred = predict_division_lik(lik, param_opt, DIVISION_FRACTIONS)

            actual_recalls = evaluate_actual_recall_at_prediction(y_pred, y_test_array, n_test)

            perf_rows.append(
                {
                    "tau": round(tau, 2),
                    "model": mk,
                    "approach": approach,
                    "param_opt": param_opt,
                    "loss": best_loss,
                    "mean_actual_recall": np.mean(actual_recalls),
                    "median_actual_recall": np.median(actual_recalls),
                    "pct_below_tau": (actual_recalls < tau).mean(),
                }
            )

    # GPCLFThreshold per-tau performance
    perf_rows.append(
        {
            "tau": round(tau, 2),
            "model": "GPCLFThreshold",
            "approach": "CLF",
            "param_opt": gpclf_model.x_thresh,
            "loss": gpclf_loss_tau,
            "mean_actual_recall": np.mean(gpclf_actual_recalls),
            "median_actual_recall": np.median(gpclf_actual_recalls),
            "pct_below_tau": gpclf_below,
        }
    )

    # ══════════════════════════════════════════════════════════
    # OPTION A: ISOTONIC CALIBRATION of Direct model predictions
    # ══════════════════════════════════════════════════════════
    from sklearn.isotonic import IsotonicRegression

    n_train = len(y_train)
    y_train_array = y_train.values.T
    gt_train = np.array([compute_ground_truth_division(y_train_array[:, s], tau) for s in range(n_train)])

    for mk_direct in ["RF", "XGBoost"]:
        mu_train = _pad_recall((rf_model if mk_direct == "RF" else xgb_model).predict(X_train_scaled).T)
        fracs_train = np.array([compute_ground_truth_division(mu_train[:, s], tau) for s in range(n_train)])
        mu_test = deterministic_preds[mk_direct]
        fracs_test = np.array([compute_ground_truth_division(mu_test[:, s], tau) for s in range(n_test)])

        if np.unique(fracs_train).size >= 5:
            iso = IsotonicRegression(out_of_bounds="clip")
            iso.fit(fracs_train, gt_train)
            fracs_cal = iso.transform(fracs_test)
        else:
            fracs_cal = fracs_test.copy()

        cal_loss = asymmetric_loss(ground_truth, fracs_cal, TARGET_OVER, TARGET_UNDER)
        cal_recalls = evaluate_actual_recall_at_prediction(fracs_cal, y_test_array, n_test)

        perf_rows.append(
            {
                "tau": round(tau, 2),
                "model": mk_direct,
                "approach": "Direct+ISO",
                "param_opt": round(float(fracs_cal.mean() - fracs_test.mean()), 4),
                "loss": cal_loss,
                "mean_actual_recall": np.mean(cal_recalls),
                "median_actual_recall": np.median(cal_recalls),
                "pct_below_tau": (cal_recalls < tau).mean(),
            }
        )
        print(
            f"  {mk_direct} + isotonic: loss={cal_loss:.4f}, "
            f"mean_recall={np.mean(cal_recalls):.3f} "
            f"(shift={fracs_cal.mean() - fracs_test.mean():+.4f})"
        )

    # ══════════════════════════════════════════════════════════
    # OPTION B: DIRECT FRACTION PREDICTION with asymmetric weights
    # ══════════════════════════════════════════════════════════
    # Train single-output models to predict the tau-crossing fraction directly.
    # Asymmetric sample weights penalise overestimation more.

    def _asymmetric_weights(gt, over_cost=3.0, under_cost=1.0, floor=0.05):
        """Samples with small true fractions get higher weight to fight overestimation."""
        base = 1.0 / np.clip(gt, floor, 1.0)
        return base * under_cost + (1.0 - base / base.max()) * over_cost

    for frac_name, model_cls, frac_params in [
        ("DirectXGB", XGBRegressor, {"n_estimators": 200, "random_state": 42, "verbosity": 0}),
        ("DirectRF", RandomForestRegressor, {"n_estimators": 200, "random_state": 42}),
    ]:
        sw = _asymmetric_weights(gt_train, over_cost=TARGET_OVER, under_cost=TARGET_UNDER)
        frac_model = model_cls(**frac_params)
        frac_model.fit(X_train_scaled, gt_train, sample_weight=sw)
        frac_preds = np.clip(frac_model.predict(X_test_scaled), 0.0, 1.0)

        frac_loss = asymmetric_loss(ground_truth, frac_preds, TARGET_OVER, TARGET_UNDER)
        frac_recalls = evaluate_actual_recall_at_prediction(frac_preds, y_test_array, n_test)

        perf_rows.append(
            {
                "tau": round(tau, 2),
                "model": frac_name,
                "approach": "Direct",
                "param_opt": 0.0,
                "loss": frac_loss,
                "mean_actual_recall": np.mean(frac_recalls),
                "median_actual_recall": np.median(frac_recalls),
                "pct_below_tau": (frac_recalls < tau).mean(),
            }
        )
        print(f"  {frac_name}: loss={frac_loss:.4f}, mean_recall={np.mean(frac_recalls):.3f}")

    # GP-DirectXGB: stacked on GP-Matern-2.5 predictions at subset divisions
    sw = _asymmetric_weights(gt_train, over_cost=TARGET_OVER, under_cost=TARGET_UNDER)
    gp_meta = XGBRegressor(n_estimators=200, random_state=42, verbosity=0)
    gp_meta.fit(gp_feat_train, gt_train, sample_weight=sw)
    gp_meta_preds = np.clip(gp_meta.predict(gp_feat_test), 0.0, 1.0)

    gp_meta_loss = asymmetric_loss(ground_truth, gp_meta_preds, TARGET_OVER, TARGET_UNDER)
    gp_meta_recalls = evaluate_actual_recall_at_prediction(gp_meta_preds, y_test_array, n_test)

    perf_rows.append(
        {
            "tau": round(tau, 2),
            "model": "GP-DirectXGB",
            "approach": "Direct",
            "param_opt": 0.0,
            "loss": gp_meta_loss,
            "mean_actual_recall": np.mean(gp_meta_recalls),
            "median_actual_recall": np.median(gp_meta_recalls),
            "pct_below_tau": (gp_meta_recalls < tau).mean(),
        }
    )
    print(f"  GP-DirectXGB: loss={gp_meta_loss:.4f}, mean_recall={np.mean(gp_meta_recalls):.3f}")

    # ══════════════════════════════════════════════════════════
    # 8c. CALIBRATION SCATTER (only for BASELINE_TAU)
    # ══════════════════════════════════════════════════════════
    if tau == BASELINE_TAU:
        print("\n--- 8c. Calibration scatter ---")
        # Standard models: GP-based + deterministic
        for mk in list(posterior.keys()) + list(deterministic_preds.keys()):
            for approach in ["CLF", "Likelihood", "Direct"]:
                sub = summary_df[
                    (summary_df["model"] == mk)
                    & (summary_df["approach"] == approach)
                    & (summary_df["over_cost"] == TARGET_OVER)
                    & (summary_df["under_cost"] == TARGET_UNDER)
                ]
                if sub.empty:
                    continue
                best = sub.loc[sub["loss"].idxmin()]
                param_opt = best["param"]

                if approach == "Direct":
                    mu = deterministic_preds[mk]
                    y_pred = np.array([compute_ground_truth_division(mu[:, s], tau) for s in range(n_test)])
                else:
                    prob = precomputed[mk]["prob"]
                    lik = precomputed[mk]["lik"]
                    if approach == "CLF":
                        y_pred = predict_division_clf(prob, param_opt, DIVISION_FRACTIONS)
                    else:
                        y_pred = predict_division_lik(lik, param_opt, DIVISION_FRACTIONS)

                fig, axes = plt.subplots(1, 2, figsize=(12, 5))
                ax = axes[0]
                over = y_pred > ground_truth
                under = y_pred <= ground_truth
                ax.scatter(ground_truth[over], y_pred[over], c="red", alpha=0.5, s=20, label=f"Over ({over.sum()})")
                ax.scatter(
                    ground_truth[under], y_pred[under], c="blue", alpha=0.5, s=20, label=f"Under ({under.sum()})"
                )
                ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
                ax.set_xlabel("True fraction")
                ax.set_ylabel("Predicted fraction")
                ax.set_title(f"{mk} \u2014 {approach} (\u03c4={tau}, param={param_opt:.2f})")
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)
                ax.set_aspect("equal")

                ax = axes[1]
                errors = y_pred - ground_truth
                ax.hist(errors, bins=30, alpha=0.7, color="gray", edgecolor="black")
                ax.axvline(0, color="black", ls="--", lw=1)
                ax.axvline(np.median(errors), color="red", ls="-", lw=1.5, label=f"Median: {np.median(errors):.3f}")
                ax.set_xlabel("Predicted \u2212 True")
                ax.set_ylabel("Count")
                ax.set_title(f"Error distribution (loss={best['loss']:.4f})")
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)

                plt.tight_layout()
                fname = f"calibration_{model_file_safe(mk)}_{approach}_tau{tau:.2f}.png"
                fig.savefig(SAVE_DIR / fname, dpi=120, bbox_inches="tight")
                plt.close(fig)
                print(f"  Saved {fname}")

        # GPCLFThreshold
        y_true_gpclf_sc = compute_ground_truth_for_tau(y_test_array, tau, n_test)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        ax = axes[0]
        over = gpclf_y_pred > y_true_gpclf_sc
        under = gpclf_y_pred <= y_true_gpclf_sc
        ax.scatter(y_true_gpclf_sc[over], gpclf_y_pred[over], c="red", alpha=0.5, s=20, label=f"Over ({over.sum()})")
        ax.scatter(
            y_true_gpclf_sc[under], gpclf_y_pred[under], c="blue", alpha=0.5, s=20, label=f"Under ({under.sum()})"
        )
        ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
        ax.set_xlabel("True fraction")
        ax.set_ylabel("Predicted fraction")
        ax.set_title(f"GPCLFThreshold (\u03c4={gpclf_model.tau:.2f}, X={gpclf_model.x_thresh:.2f})")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal")
        ax = axes[1]
        errors = gpclf_y_pred - y_true_gpclf_sc
        ax.hist(errors, bins=30, alpha=0.7, color="gray", edgecolor="black")
        ax.axvline(0, color="black", ls="--", lw=1)
        ax.axvline(np.median(errors), color="red", ls="-", lw=1.5, label=f"Median: {np.median(errors):.3f}")
        ax.set_xlabel("Predicted \u2212 True")
        ax.set_ylabel("Count")
        ax.set_title("Error distribution")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        fname = f"calibration_GPCLFThreshold_CLF_tau{tau:.2f}.png"
        fig.savefig(SAVE_DIR / fname, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {fname}")

        # Direct+ISO (Option A)
        for mk_direct in ["RF", "XGBoost"]:
            mu_train = _pad_recall((rf_model if mk_direct == "RF" else xgb_model).predict(X_train_scaled).T)
            fracs_train = np.array([compute_ground_truth_division(mu_train[:, s], tau) for s in range(n_train)])
            mu_test = deterministic_preds[mk_direct]
            fracs_test = np.array([compute_ground_truth_division(mu_test[:, s], tau) for s in range(n_test)])

            if np.unique(fracs_train).size >= 5:
                iso = IsotonicRegression(out_of_bounds="clip")
                iso.fit(fracs_train, gt_train)
                fracs_cal = iso.transform(fracs_test)
            else:
                fracs_cal = fracs_test.copy()

            cal_loss = asymmetric_loss(ground_truth, fracs_cal, TARGET_OVER, TARGET_UNDER)

            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            ax = axes[0]
            over = fracs_cal > ground_truth
            under = fracs_cal <= ground_truth
            ax.scatter(ground_truth[over], fracs_cal[over], c="red", alpha=0.5, s=20, label=f"Over ({over.sum()})")
            ax.scatter(ground_truth[under], fracs_cal[under], c="blue", alpha=0.5, s=20, label=f"Under ({under.sum()})")
            ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
            ax.set_xlabel("True fraction")
            ax.set_ylabel("Predicted fraction")
            ax.set_title(f"{mk_direct} \u2014 Direct+ISO (\u03c4={tau})")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            ax.set_aspect("equal")
            ax = axes[1]
            errors = fracs_cal - ground_truth
            ax.hist(errors, bins=30, alpha=0.7, color="gray", edgecolor="black")
            ax.axvline(0, color="black", ls="--", lw=1)
            ax.axvline(np.median(errors), color="red", ls="-", lw=1.5, label=f"Median: {np.median(errors):.3f}")
            ax.set_xlabel("Predicted \u2212 True")
            ax.set_ylabel("Count")
            ax.set_title(f"Error distribution (loss={cal_loss:.4f})")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            fname = f"calibration_{mk_direct}_Direct+ISO_tau{tau:.2f}.png"
            fig.savefig(SAVE_DIR / fname, dpi=120, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved {fname}")

        # DirectXGB, DirectRF (Option B)
        for frac_name, model_cls, frac_params in [
            ("DirectXGB", XGBRegressor, {"n_estimators": 200, "random_state": 42, "verbosity": 0}),
            ("DirectRF", RandomForestRegressor, {"n_estimators": 200, "random_state": 42}),
        ]:
            sw = _asymmetric_weights(gt_train, over_cost=TARGET_OVER, under_cost=TARGET_UNDER)
            frac_model = model_cls(**frac_params)
            frac_model.fit(X_train_scaled, gt_train, sample_weight=sw)
            frac_preds = np.clip(frac_model.predict(X_test_scaled), 0.0, 1.0)
            frac_loss = asymmetric_loss(ground_truth, frac_preds, TARGET_OVER, TARGET_UNDER)

            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            ax = axes[0]
            over = frac_preds > ground_truth
            under = frac_preds <= ground_truth
            ax.scatter(ground_truth[over], frac_preds[over], c="red", alpha=0.5, s=20, label=f"Over ({over.sum()})")
            ax.scatter(
                ground_truth[under], frac_preds[under], c="blue", alpha=0.5, s=20, label=f"Under ({under.sum()})"
            )
            ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
            ax.set_xlabel("True fraction")
            ax.set_ylabel("Predicted fraction")
            ax.set_title(f"{frac_name} \u2014 Direct (\u03c4={tau})")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            ax.set_aspect("equal")
            ax = axes[1]
            errors = frac_preds - ground_truth
            ax.hist(errors, bins=30, alpha=0.7, color="gray", edgecolor="black")
            ax.axvline(0, color="black", ls="--", lw=1)
            ax.axvline(np.median(errors), color="red", ls="-", lw=1.5, label=f"Median: {np.median(errors):.3f}")
            ax.set_xlabel("Predicted \u2212 True")
            ax.set_ylabel("Count")
            ax.set_title(f"Error distribution (loss={frac_loss:.4f})")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            fname = f"calibration_{frac_name}_Direct_tau{tau:.2f}.png"
            fig.savefig(SAVE_DIR / fname, dpi=120, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved {fname}")

        # GP-DirectXGB calibration scatter
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        ax = axes[0]
        over = gp_meta_preds > ground_truth
        under = gp_meta_preds <= ground_truth
        ax.scatter(ground_truth[over], gp_meta_preds[over], c="red", alpha=0.5, s=20, label=f"Over ({over.sum()})")
        ax.scatter(ground_truth[under], gp_meta_preds[under], c="blue", alpha=0.5, s=20, label=f"Under ({under.sum()})")
        ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
        ax.set_xlabel("True fraction")
        ax.set_ylabel("Predicted fraction")
        ax.set_title(f"GP-DirectXGB \u2014 Direct (\u03c4={tau})")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal")
        ax = axes[1]
        errors = gp_meta_preds - ground_truth
        ax.hist(errors, bins=30, alpha=0.7, color="gray", edgecolor="black")
        ax.axvline(0, color="black", ls="--", lw=1)
        ax.axvline(np.median(errors), color="red", ls="-", lw=1.5, label=f"Median: {np.median(errors):.3f}")
        ax.set_xlabel("Predicted \u2212 True")
        ax.set_ylabel("Count")
        ax.set_title(f"Error distribution (loss={gp_meta_loss:.4f})")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        fname = f"calibration_GP-DirectXGB_Direct_tau{tau:.2f}.png"
        fig.savefig(SAVE_DIR / fname, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {fname}")

    # ══════════════════════════════════════════════════════════
    # 11b. ACTUAL RECALL AT INDEX (only for BASELINE_TAU)
    # ══════════════════════════════════════════════════════════
    if tau == BASELINE_TAU:
        # 11a. Recall landscape
        print("\n--- 11a. Recall landscape ---")
        fig, ax = plt.subplots(figsize=(12, 5))
        bp_data = [y_test_array[d, :] for d in range(N_DIVISIONS)]
        ax.boxplot(
            bp_data,
            positions=range(1, N_DIVISIONS + 1),
            widths=0.6,
            patch_artist=True,
            boxprops=dict(facecolor="lightblue", alpha=0.7),
            medianprops=dict(color="red", lw=2),
            whiskerprops=dict(lw=1.5),
            capprops=dict(lw=1.5),
        )
        means = [np.mean(y_test_array[d, :]) for d in range(N_DIVISIONS)]
        ax.plot(range(1, N_DIVISIONS + 1), means, "ko-", lw=2, markersize=5, label="Mean recall")
        ax.axhline(tau, color="green", ls="--", lw=1.5, label=f"\u03c4 = {tau}")
        ax.set_xlabel("Division p", fontsize=12)
        ax.set_ylabel("Recall", fontsize=12)
        ax.set_title("Recall landscape across test samples", fontsize=13, fontweight="bold")
        ax.set_xticks(range(1, N_DIVISIONS + 1))
        ax.set_xticklabels([f"{p}\n({DIVISION_FRACTIONS[p - 1]:.0%})" for p in range(1, N_DIVISIONS + 1)], fontsize=7)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, axis="y")
        plt.tight_layout()
        fig.savefig(SAVE_DIR / "recall_landscape.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        print("  Saved recall_landscape.png")

        # 11b. Actual recall at predicted last-index
        print("\n--- 11b. Actual recall at predicted last-index ---")
        dist_summary_rows = []
        for mk in model_keys + list(deterministic_preds.keys()):
            for approach in ["CLF", "Likelihood", "Direct"]:
                sub = summary_df[
                    (summary_df["model"] == mk)
                    & (summary_df["approach"] == approach)
                    & (summary_df["over_cost"] == TARGET_OVER)
                    & (summary_df["under_cost"] == TARGET_UNDER)
                ]
                if sub.empty:
                    continue
                best = sub.loc[sub["loss"].idxmin()]
                param_opt = best["param"]
                best_loss = best["loss"]

                if approach == "Direct":
                    mu = deterministic_preds[mk]
                    y_pred = np.array([compute_ground_truth_division(mu[:, s], tau) for s in range(n_test)])
                else:
                    prob = precomputed[mk]["prob"]
                    lik = precomputed[mk]["lik"]
                    if approach == "CLF":
                        y_pred = predict_division_clf(prob, param_opt, DIVISION_FRACTIONS)
                    else:
                        y_pred = predict_division_lik(lik, param_opt, DIVISION_FRACTIONS)

                actual_recalls = evaluate_actual_recall_at_prediction(y_pred, y_test_array, n_test)
                below_tau = (actual_recalls < tau).mean()

                dist_summary_rows.append(
                    {
                        "model": mk,
                        "approach": approach,
                        "param_opt": param_opt,
                        "loss": best_loss,
                        "mean_actual_recall": np.mean(actual_recalls),
                        "median_actual_recall": np.median(actual_recalls),
                        "pct_below_tau": below_tau,
                    }
                )

                fig, axes = plt.subplots(1, 2, figsize=(13, 5))
                ax = axes[0]
                ax.hist(actual_recalls, bins=25, alpha=0.7, color="steelblue", edgecolor="black", density=True)
                ax.axvline(tau, color="red", ls="--", lw=2, label=f"\u03c4 = {tau}")
                ax.axvline(
                    np.mean(actual_recalls),
                    color="darkorange",
                    ls="-",
                    lw=2,
                    label=f"Mean = {np.mean(actual_recalls):.3f}",
                )
                ax.axvline(
                    np.median(actual_recalls),
                    color="green",
                    ls=":",
                    lw=2,
                    label=f"Median = {np.median(actual_recalls):.3f}",
                )
                ax.set_xlabel("Actual recall at predicted last-index", fontsize=11)
                ax.set_ylabel("Density", fontsize=11)
                ax.set_title(
                    f"{mk} \u2014 {approach}\n\u03c4={tau}, param={param_opt:.2f}, loss={best_loss:.4f}", fontsize=10
                )
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3, axis="y")

                ax = axes[1]
                n_show = min(50, n_test)
                rng = np.random.RandomState(42)
                sample_idxs = rng.choice(n_test, n_show, replace=False)
                for s in sample_idxs:
                    ax.plot(DIVISION_FRACTIONS, y_test_array[:, s], color="grey", alpha=0.15, lw=0.5)
                means = [np.mean(y_test_array[d, :]) for d in range(N_DIVISIONS)]
                ax.plot(DIVISION_FRACTIONS, means, "k-", lw=2, label="Mean recall")
                jitter = rng.uniform(-0.01, 0.01, n_test)
                for s in range(n_test):
                    color = "red" if actual_recalls[s] < tau else "blue"
                    ax.scatter(y_pred[s] + jitter[s], actual_recalls[s], c=color, alpha=0.4, s=8, edgecolors="none")
                ax.axhline(tau, color="red", ls="--", lw=1.5)
                ax.set_xlabel("Predicted last-index (fraction)", fontsize=11)
                ax.set_ylabel("Actual recall at that index", fontsize=11)
                ax.set_title("Readout: actual recall at each predicted index", fontsize=10)
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)
                plt.tight_layout()
                fname = f"actual_recall_at_index_{model_file_safe(mk)}_{approach}_tau{tau:.2f}.png"
                fig.savefig(SAVE_DIR / fname, dpi=120, bbox_inches="tight")
                plt.close(fig)
                print(f"  Saved {fname}")

        # Extra models in dist_summary_rows: GPCLFThreshold, Direct+ISO, DirectXGB, DirectRF
        dist_summary_rows.append(
            {
                "model": "GPCLFThreshold",
                "approach": "CLF",
                "param_opt": gpclf_model.x_thresh,
                "loss": gpclf_loss_tau,
                "mean_actual_recall": np.mean(gpclf_actual_recalls),
                "median_actual_recall": np.median(gpclf_actual_recalls),
                "pct_below_tau": gpclf_below,
            }
        )
        for mk_direct in ["RF", "XGBoost"]:
            mu_train = _pad_recall((rf_model if mk_direct == "RF" else xgb_model).predict(X_train_scaled).T)
            fracs_train = np.array([compute_ground_truth_division(mu_train[:, s], tau) for s in range(n_train)])
            mu_test = deterministic_preds[mk_direct]
            fracs_test = np.array([compute_ground_truth_division(mu_test[:, s], tau) for s in range(n_test)])
            if np.unique(fracs_train).size >= 5:
                iso = IsotonicRegression(out_of_bounds="clip")
                iso.fit(fracs_train, gt_train)
                fracs_cal = iso.transform(fracs_test)
            else:
                fracs_cal = fracs_test.copy()
            cal_loss = asymmetric_loss(ground_truth, fracs_cal, TARGET_OVER, TARGET_UNDER)
            cal_recalls = evaluate_actual_recall_at_prediction(fracs_cal, y_test_array, n_test)
            dist_summary_rows.append(
                {
                    "model": mk_direct,
                    "approach": "Direct+ISO",
                    "param_opt": round(float(fracs_cal.mean() - fracs_test.mean()), 4),
                    "loss": cal_loss,
                    "mean_actual_recall": np.mean(cal_recalls),
                    "median_actual_recall": np.median(cal_recalls),
                    "pct_below_tau": (cal_recalls < tau).mean(),
                }
            )
        for frac_name, model_cls, frac_params in [
            ("DirectXGB", XGBRegressor, {"n_estimators": 200, "random_state": 42, "verbosity": 0}),
            ("DirectRF", RandomForestRegressor, {"n_estimators": 200, "random_state": 42}),
        ]:
            sw = _asymmetric_weights(gt_train, over_cost=TARGET_OVER, under_cost=TARGET_UNDER)
            frac_model = model_cls(**frac_params)
            frac_model.fit(X_train_scaled, gt_train, sample_weight=sw)
            frac_preds = np.clip(frac_model.predict(X_test_scaled), 0.0, 1.0)
            frac_loss = asymmetric_loss(ground_truth, frac_preds, TARGET_OVER, TARGET_UNDER)
            frac_recalls = evaluate_actual_recall_at_prediction(frac_preds, y_test_array, n_test)
            dist_summary_rows.append(
                {
                    "model": frac_name,
                    "approach": "Direct",
                    "param_opt": 0.0,
                    "loss": frac_loss,
                    "mean_actual_recall": np.mean(frac_recalls),
                    "median_actual_recall": np.median(frac_recalls),
                    "pct_below_tau": (frac_recalls < tau).mean(),
                }
            )
        dist_summary_rows.append(
            {
                "model": "GP-DirectXGB",
                "approach": "Direct",
                "param_opt": 0.0,
                "loss": gp_meta_loss,
                "mean_actual_recall": np.mean(gp_meta_recalls),
                "median_actual_recall": np.median(gp_meta_recalls),
                "pct_below_tau": (gp_meta_recalls < tau).mean(),
            }
        )

        # 11c. GPCLFThreshold calibration
        print("\n--- 11c. GPCLFThreshold calibration ---")
        y_true_gpclf = compute_ground_truth_for_tau(y_test_array, tau, n_test)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        ax = axes[0]
        over = gpclf_y_pred > y_true_gpclf
        under = gpclf_y_pred <= y_true_gpclf
        ax.scatter(y_true_gpclf[over], gpclf_y_pred[over], c="red", alpha=0.5, s=20, label=f"Over ({over.sum()})")
        ax.scatter(y_true_gpclf[under], gpclf_y_pred[under], c="blue", alpha=0.5, s=20, label=f"Under ({under.sum()})")
        ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
        ax.set_xlabel("True fraction")
        ax.set_ylabel("Predicted fraction")
        ax.set_title(
            f"GPCLFThreshold (\u03c4={gpclf_model.tau:.2f}, X={gpclf_model.x_thresh:.2f})\nvs target \u03c4={tau}"
        )
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal")

        ax = axes[1]
        errors = gpclf_y_pred - y_true_gpclf
        ax.hist(errors, bins=30, alpha=0.7, color="gray", edgecolor="black")
        ax.axvline(0, color="black", ls="--", lw=1)
        ax.axvline(np.median(errors), color="red", ls="-", lw=1.5, label=f"Median: {np.median(errors):.3f}")
        ax.set_xlabel("Predicted \u2212 True")
        ax.set_ylabel("Count")
        ax.set_title("Error distribution")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        fname = f"calibration_GPCLFThreshold_CLF_tau{tau:.2f}.png"
        fig.savefig(SAVE_DIR / fname, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {fname}")

        # List savings for GPCLFThreshold
        print("\n--- GPCLFThreshold fraction of list saved ---")
        savings = 1.0 - gpclf_y_pred
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(savings, bins=30, alpha=0.7, color="seagreen", edgecolor="black", density=True)
        ax.axvline(np.median(savings), color="darkgreen", ls="-", lw=2, label=f"Median: {np.median(savings):.3f}")
        ax.axvline(np.mean(savings), color="lime", ls="--", lw=2, label=f"Mean: {np.mean(savings):.3f}")
        ax.set_xlabel("Fraction of ranked list saved (1 \u2212 predicted index)")
        ax.set_ylabel("Density")
        ax.set_title("GPCLFThreshold \u2014 fraction of list saved")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, axis="y")
        plt.tight_layout()
        fig.savefig(SAVE_DIR / "list_savings_GPCLFThreshold.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        print("  Saved list_savings_GPCLFThreshold.png")

        # 11d. ROC-style analysis
        print("\n--- 11d. ROC-style analysis ---")

        def _tpr_fpr(y_pred):
            tpr = float(np.mean([actual_recall_at_fraction(y_test_array[:, s], y_pred[s]) for s in range(n_test)]))
            fpr = float(np.mean(y_pred))
            return tpr, fpr

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, label="Random (AUC=0.500)")
        colors = plt.cm.tab10(np.linspace(0, 1, 10))
        color_idx = 0
        roc_summary = {}

        for mk in model_keys:
            prob = precomputed[mk]["prob"]
            lik = precomputed[mk]["lik"]
            base_color = colors[color_idx % len(colors)]
            color_idx += 1
            for approach, param_values, is_clf in [
                ("CLF", X_VALUES, True),
                ("Likelihood", LIK_TAU_VALUES, False),
            ]:
                fpr_list, tpr_list = [], []
                for p in param_values:
                    if is_clf:
                        y_pred = predict_division_clf(prob, p, DIVISION_FRACTIONS)
                    else:
                        y_pred = predict_division_lik(lik, p, DIVISION_FRACTIONS)
                    tpr, fpr = _tpr_fpr(y_pred)
                    fpr_list.append(fpr)
                    tpr_list.append(tpr)
                auc_val = np.trapz(tpr_list, fpr_list)
                ls = "-" if is_clf else "--"
                ax.plot(
                    fpr_list,
                    tpr_list,
                    color=base_color,
                    ls=ls,
                    lw=1.5,
                    label=f"{mk} \u2014 {approach} (AUC={auc_val:.3f})",
                )
                best_sub = summary_df[
                    (summary_df["model"] == mk)
                    & (summary_df["approach"] == approach)
                    & (summary_df["over_cost"] == TARGET_OVER)
                    & (summary_df["under_cost"] == TARGET_UNDER)
                ]
                if not best_sub.empty:
                    best_row = best_sub.loc[best_sub["loss"].idxmin()]
                    best_param = best_row["param"]
                    if is_clf:
                        y_pred_opt = predict_division_clf(prob, best_param, DIVISION_FRACTIONS)
                    else:
                        y_pred_opt = predict_division_lik(lik, best_param, DIVISION_FRACTIONS)
                    tpr_opt, fpr_opt = _tpr_fpr(y_pred_opt)
                else:
                    tpr_opt, fpr_opt = float("nan"), float("nan")
                roc_summary[(mk, approach)] = {
                    "auc": auc_val,
                    "mean_fpr": fpr_opt,
                    "mean_tpr": tpr_opt,
                }

        for mk, mu in deterministic_preds.items():
            y_pred = np.array([compute_ground_truth_division(mu[:, s], tau) for s in range(n_test)])
            tpr, fpr = _tpr_fpr(y_pred)
            ax.scatter(fpr, tpr, s=80, marker="D", label=f"{mk} \u2014 Direct (TPR={tpr:.3f}, FPR={fpr:.3f})")
            roc_summary[(mk, "Direct")] = {
                "auc": float("nan"),
                "mean_fpr": fpr,
                "mean_tpr": tpr,
            }

        tpr_gpclf = float(np.mean(gpclf_actual_recalls))
        fpr_gpclf = float(np.mean(gpclf_y_pred))
        ax.scatter(
            fpr_gpclf,
            tpr_gpclf,
            s=80,
            marker="s",
            label=f"GPCLFThreshold \u2014 CLF (TPR={tpr_gpclf:.3f}, FPR={fpr_gpclf:.3f})",
        )
        roc_summary[("GPCLFThreshold", "CLF")] = {
            "auc": float("nan"),
            "mean_fpr": fpr_gpclf,
            "mean_tpr": tpr_gpclf,
        }

        # Direct+ISO (Option A)
        for mk_direct in ["RF", "XGBoost"]:
            mu_train = _pad_recall((rf_model if mk_direct == "RF" else xgb_model).predict(X_train_scaled).T)
            fracs_train = np.array([compute_ground_truth_division(mu_train[:, s], tau) for s in range(n_train)])
            mu_test = deterministic_preds[mk_direct]
            fracs_test = np.array([compute_ground_truth_division(mu_test[:, s], tau) for s in range(n_test)])
            if np.unique(fracs_train).size >= 5:
                iso = IsotonicRegression(out_of_bounds="clip")
                iso.fit(fracs_train, gt_train)
                fracs_cal = iso.transform(fracs_test)
            else:
                fracs_cal = fracs_test.copy()
            tpr_iso, fpr_iso = _tpr_fpr(fracs_cal)
            roc_summary[(mk_direct, "Direct+ISO")] = {
                "auc": float("nan"),
                "mean_fpr": fpr_iso,
                "mean_tpr": tpr_iso,
            }
            ax.scatter(
                fpr_iso,
                tpr_iso,
                s=80,
                marker="^",
                label=f"{mk_direct} \u2014 Direct+ISO (TPR={tpr_iso:.3f}, FPR={fpr_iso:.3f})",
            )

        # DirectXGB, DirectRF (Option B)
        for frac_name, model_cls, frac_params in [
            ("DirectXGB", XGBRegressor, {"n_estimators": 200, "random_state": 42, "verbosity": 0}),
            ("DirectRF", RandomForestRegressor, {"n_estimators": 200, "random_state": 42}),
        ]:
            sw = _asymmetric_weights(gt_train, over_cost=TARGET_OVER, under_cost=TARGET_UNDER)
            frac_model = model_cls(**frac_params)
            frac_model.fit(X_train_scaled, gt_train, sample_weight=sw)
            frac_preds = np.clip(frac_model.predict(X_test_scaled), 0.0, 1.0)
            tpr_frac, fpr_frac = _tpr_fpr(frac_preds)
            roc_summary[(frac_name, "Direct")] = {
                "auc": float("nan"),
                "mean_fpr": fpr_frac,
                "mean_tpr": tpr_frac,
            }
            ax.scatter(
                fpr_frac,
                tpr_frac,
                s=80,
                marker="v",
                label=f"{frac_name} \u2014 Direct (TPR={tpr_frac:.3f}, FPR={fpr_frac:.3f})",
            )

        # GP-DirectXGB (stacked meta-model)
        tpr_gpmeta, fpr_gpmeta = _tpr_fpr(gp_meta_preds)
        roc_summary[("GP-DirectXGB", "Direct")] = {
            "auc": float("nan"),
            "mean_fpr": fpr_gpmeta,
            "mean_tpr": tpr_gpmeta,
        }
        ax.scatter(
            fpr_gpmeta,
            tpr_gpmeta,
            s=80,
            marker="*",
            label=f"GP-DirectXGB \u2014 Direct (TPR={tpr_gpmeta:.3f}, FPR={fpr_gpmeta:.3f})",
        )

        ax.set_xlabel("Mean FPR (predicted fraction)", fontsize=12)
        ax.set_ylabel("Mean TPR (actual recall)", fontsize=12)
        ax.set_title("ROC-style: TPR vs FPR across parameter sweep", fontsize=13, fontweight="bold")
        ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        plt.tight_layout()
        fig.savefig(SAVE_DIR / "roc_curves.png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        print("  Saved roc_curves.png")
        print(f"  ROC summary computed for {len(roc_summary)} model/approach pairs")

        # Back-fill AUC, FPR, TPR into dist_summary_rows
        for row in dist_summary_rows:
            key = (row["model"], row["approach"])
            if key in roc_summary:
                s = roc_summary[key]
                row["auc"] = round(s["auc"], 4) if not np.isnan(s["auc"]) else float("nan")
                row["mean_fpr"] = round(s["mean_fpr"], 4)
                row["mean_tpr"] = round(s["mean_tpr"], 4)
            else:
                row["auc"] = float("nan")
                row["mean_fpr"] = float("nan")
                row["mean_tpr"] = float("nan")

        dist_summary_df = pd.DataFrame(dist_summary_rows)
        print("\nActual recall at predicted last-index \u2014 summary:")
        print(dist_summary_df.to_string(index=False))
        dist_summary_df.to_csv(SAVE_DIR / "actual_recall_at_index_summary.csv", index=False)


# ═══════════════════════════════════════════════════════════════
# 12. CROSS-TAU COMPARISON
# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("12. CROSS-TAU COMPARISON")
print("=" * 60)

perf_df = pd.DataFrame(perf_rows)

# 12a. Best performance per tau
print("\n--- 12a. Best mean_actual_recall per tau ---")
best_per_tau = perf_df.loc[perf_df.groupby("tau")["mean_actual_recall"].idxmax()].reset_index(drop=True)
print(best_per_tau.to_string(index=False))
best_per_tau.to_csv(SAVE_DIR / "tau_sweep_summary.csv", index=False)
print("  Saved tau_sweep_summary.csv")

# Identify optimal tau
optimal_row = best_per_tau.loc[best_per_tau["mean_actual_recall"].idxmax()]
opt_tau = optimal_row["tau"]
print(
    f"\nOptimal \u03c4 = {opt_tau:.2f} "
    f"(mean actual recall = {optimal_row['mean_actual_recall']:.4f}, "
    f"model = {optimal_row['model']}, approach = {optimal_row['approach']})"
)

# 12b. Mean actual recall vs tau
print("\n--- 12b. Actual recall vs \u03c4 ---")
fig, ax = plt.subplots(figsize=(10, 6))
for (mk, approach), g in perf_df.groupby(["model", "approach"]):
    g = g.sort_values("tau")
    ax.plot(g["tau"], g["mean_actual_recall"], "o-", lw=1.5, label=f"{mk} \u2014 {approach}")
ax.plot([0.82, 1.02], [0.82, 1.02], "k--", alpha=0.3, label="y = x (perfect)")
ax.axvline(opt_tau, color="red", ls=":", lw=1, alpha=0.7, label=f"Optimal \u03c4 = {opt_tau:.2f}")
ax.set_xlabel("Target recall \u03c4", fontsize=12)
ax.set_ylabel("Mean actual recall at prediction", fontsize=12)
ax.set_title("Actual recall vs target \u03c4", fontsize=13, fontweight="bold")
ax.legend(fontsize=7, ncol=2)
ax.grid(True, alpha=0.3)
ax.set_xlim(0.83, 1.02)
ax.set_ylim(0.83, 1.02)
plt.tight_layout()
fig.savefig(SAVE_DIR / "tau_sweep_actual_recall.png", dpi=120, bbox_inches="tight")
plt.close(fig)
print("  Saved tau_sweep_actual_recall.png")

# 12c. Tradeoff: loss vs actual recall
print("\n--- 12c. Loss vs actual recall tradeoff ---")
fig, ax = plt.subplots(figsize=(10, 7))
tau_colors = {t: plt.cm.viridis(i / max(len(TARGET_TAUS) - 1, 1)) for i, t in enumerate(sorted(TARGET_TAUS))}
markers = {"CLF": "o", "Likelihood": "s", "Direct": "D"}
seen_labels = set()
for _, row in perf_df.iterrows():
    c = tau_colors[row["tau"]]
    m = markers.get(row["approach"], "o")
    label = f"\u03c4={row['tau']:.2f}"
    if label not in seen_labels:
        seen_labels.add(label)
        ax.scatter(row["loss"], row["mean_actual_recall"], c=[c], marker=m, s=60, alpha=0.8, label=label)
    else:
        ax.scatter(row["loss"], row["mean_actual_recall"], c=[c], marker=m, s=60, alpha=0.8)
    ax.annotate(str(row["model"])[:10], (row["loss"], row["mean_actual_recall"]), fontsize=5, alpha=0.7)
ax.set_xlabel("Best loss (asymmetric, 3:1)", fontsize=12)
ax.set_ylabel("Mean actual recall", fontsize=12)
ax.set_title("Tradeoff: loss vs actual recall across \u03c4 values", fontsize=13, fontweight="bold")
ax.legend(fontsize=8, ncol=2)
ax.grid(True, alpha=0.3)
plt.tight_layout()
fig.savefig(SAVE_DIR / "tau_sweep_tradeoff.png", dpi=120, bbox_inches="tight")
plt.close(fig)
print("  Saved tau_sweep_tradeoff.png")

# 12d. Improvement over baseline
print("\n--- 12d. Improvement over baseline ---")
baseline_perf = perf_df[perf_df["tau"] == BASELINE_TAU].copy()
best_baseline = baseline_perf.loc[baseline_perf["mean_actual_recall"].idxmax()]
improvement = optimal_row["mean_actual_recall"] - best_baseline["mean_actual_recall"]
print(
    f"  Baseline (\u03c4={BASELINE_TAU:.2f}): "
    f"best mean actual recall = {best_baseline['mean_actual_recall']:.4f} "
    f"({best_baseline['model']} \u2014 {best_baseline['approach']})"
)
print(
    f"  Optimal (\u03c4={opt_tau:.2f}): "
    f"best mean actual recall = {optimal_row['mean_actual_recall']:.4f} "
    f"({optimal_row['model']} \u2014 {optimal_row['approach']})"
)
print(f"  Improvement: {improvement:+.4f} ({improvement * 100:+.2f}%)")


print("\n" + "=" * 60)
print("DONE")
print("=" * 60)
