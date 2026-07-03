"""
Predict input taxonomic composition from m_stats features.

This script merges taxon-count regression and binned-oridnal classification.
It predicts both the continuous number of unique input taxids per dataset
and a binned window (very_low / low / medium / high / very_high).

Regression outputs (continuous):
  - MeanBaseline, XGBoost+GridSearch, XGBoost, RandomForest,
    GradientBoosting, LinearRegression (stats-only)

Classification outputs (binned):
  - MajorityClassBaseline, XGBClassifier+GridSearch, XGBClassifier,
    RandomForestClassifier, GradientBoostingClassifier

Outputs to ``--output_dir``:
  - metrics_regression.csv          — continuous model comparison
  - metrics_classification.csv      — binned model comparison
  - predictions_vs_actual.png       — scatter per model (regression)
  - residuals.png                   — residual histograms (regression)
  - confusion_matrix.png            — heatmap grid per classification model
  - classification_metrics.png      — bar chart comparing models
  - feature_importance.png          — top-15 (tree models)
  - training_time.png               — all models
  - detailed_predictions.tsv        — per-dataset predictions
  - best_params.txt                 — best GridSearchCV params
"""

import argparse
import json
import time
import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import (
    GradientBoostingClassifier,
    GradientBoostingRegressor,
    RandomForestClassifier,
    RandomForestRegressor,
)
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    precision_score,
    r2_score,
    recall_score,
)
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.preprocessing import StandardScaler

from deployment.model_evaluation.config import EvaluatorConfig

warnings.filterwarnings("ignore")

matplotlib.use("Agg")
sns.set_style("whitegrid")

# ── Constants ──────────────────────────────────────────────────────────────

RANDOM_STATE = 42
TEST_SIZE = 0.2

STATS_FEATURES = [
    "number_hits",
    "counts_kurtosis",
    "counts_skewness",
    "tax_diversity_shannon",
    "max_uniq_reads",
    "total_uniq_reads",
]

BIN_EDGES = [0, 3, 6, 9, 12, 100]
BIN_LABELS = ["very_low", "low", "medium", "high", "very_high"]

REGRESSION_MODELS_CFG = [
    {
        "name": "MeanBaseline",
        "module": None,
        "class_name": None,
        "use_grid": False,
        "stats_only": False,
    },
    {
        "name": "XGBoost+Grid",
        "module": "xgboost",
        "class_name": "XGBRegressor",
        "use_grid": True,
        "stats_only": False,
        "param_grid": {
            "n_estimators": [100, 300, 500],
            "max_depth": [4, 6, 8],
            "learning_rate": [0.05, 0.1, 0.2],
            "subsample": [0.7, 0.85, 1.0],
            "colsample_bytree": [0.7, 0.85, 1.0],
        },
    },
    {
        "name": "XGBoost",
        "module": "xgboost",
        "class_name": "XGBRegressor",
        "use_grid": False,
        "stats_only": False,
        "params": {
            "n_estimators": 300,
            "max_depth": 6,
            "learning_rate": 0.1,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
            "random_state": RANDOM_STATE,
        },
    },
    {
        "name": "RandomForest",
        "module": "sklearn.ensemble",
        "class_name": "RandomForestRegressor",
        "use_grid": False,
        "stats_only": False,
        "params": {
            "n_estimators": 300,
            "max_depth": 12,
            "min_samples_leaf": 3,
            "random_state": RANDOM_STATE,
            "n_jobs": -1,
        },
    },
    {
        "name": "GradientBoosting",
        "module": "sklearn.ensemble",
        "class_name": "GradientBoostingRegressor",
        "use_grid": False,
        "stats_only": False,
        "params": {
            "n_estimators": 300,
            "max_depth": 5,
            "learning_rate": 0.1,
            "subsample": 0.8,
            "random_state": RANDOM_STATE,
        },
    },
    {
        "name": "LinearRegression",
        "module": "sklearn.linear_model",
        "class_name": "LinearRegression",
        "use_grid": False,
        "stats_only": True,
        "params": {},
    },
]

CLASSIFICATION_MODELS_CFG = [
    {
        "name": "MajorityClassBaseline",
        "module": None,
        "class_name": None,
        "use_grid": False,
        "stats_only": False,
    },
    {
        "name": "XGBClassifier+Grid",
        "module": "xgboost",
        "class_name": "XGBClassifier",
        "use_grid": True,
        "stats_only": False,
        "param_grid": {
            "n_estimators": [100, 300],
            "max_depth": [4, 6, 8],
            "learning_rate": [0.05, 0.1, 0.2],
            "subsample": [0.7, 0.85, 1.0],
        },
    },
    {
        "name": "XGBClassifier",
        "module": "xgboost",
        "class_name": "XGBClassifier",
        "use_grid": False,
        "stats_only": False,
        "params": {
            "n_estimators": 300,
            "max_depth": 6,
            "learning_rate": 0.1,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
            "random_state": RANDOM_STATE,
        },
    },
    {
        "name": "RandomForest",
        "module": "sklearn.ensemble",
        "class_name": "RandomForestClassifier",
        "use_grid": False,
        "stats_only": False,
        "params": {
            "n_estimators": 300,
            "max_depth": 12,
            "min_samples_leaf": 3,
            "random_state": RANDOM_STATE,
            "n_jobs": -1,
        },
    },
    {
        "name": "GradientBoosting",
        "module": "sklearn.ensemble",
        "class_name": "GradientBoostingClassifier",
        "use_grid": False,
        "stats_only": False,
        "params": {
            "n_estimators": 300,
            "max_depth": 5,
            "learning_rate": 0.1,
            "subsample": 0.8,
            "random_state": RANDOM_STATE,
        },
    },
]


# ── Helpers ────────────────────────────────────────────────────────────────


def _import_model(cfg):
    mod = __import__(cfg["module"], fromlist=[cfg["class_name"]])
    return getattr(mod, cfg["class_name"])


def _safe_metric(y_true, y_pred, fn, **kw):
    try:
        return fn(y_true, y_pred, **kw)
    except Exception:
        return float("nan")


def evaluate_regression(y_true, y_pred):
    residual = y_true - y_pred
    return {
        "R2": _safe_metric(y_true, y_pred, r2_score),
        "MAE": _safe_metric(y_true, y_pred, mean_absolute_error),
        "RMSE": _safe_metric(y_true, y_pred, lambda a, b: mean_squared_error(a, b) ** 0.5),
        "MAPE": float(np.mean(np.abs(residual / (y_true + 1e-8))) * 100),
        "acc_1": float(np.mean(np.abs(residual) <= 1) * 100),
        "acc_2": float(np.mean(np.abs(residual) <= 2) * 100),
        "acc_5": float(np.mean(np.abs(residual) <= 5) * 100),
    }


def evaluate_classification(y_true, y_pred):
    labels = sorted(set(y_true) | set(y_pred))
    ordinal_diff = np.abs(np.asarray(y_true) - np.asarray(y_pred))
    return {
        "Accuracy": _safe_metric(y_true, y_pred, accuracy_score),
        "BalancedAccuracy": _safe_metric(y_true, y_pred, balanced_accuracy_score),
        "MacroF1": _safe_metric(y_true, y_pred, f1_score, average="macro"),
        "MacroPrecision": _safe_metric(y_true, y_pred, precision_score, average="macro", zero_division=0),
        "MacroRecall": _safe_metric(y_true, y_pred, recall_score, average="macro", zero_division=0),
        "OrdinalMAE": float(np.mean(ordinal_diff)),
        "OrdinalAcc0": float(np.mean(ordinal_diff == 0) * 100),
        "OrdinalAcc1": float(np.mean(ordinal_diff <= 1) * 100),
    }


def bin_count(y):
    return pd.cut(y, bins=BIN_EDGES, labels=BIN_LABELS, right=True, include_lowest=True).astype(str)


class MeanBaseline:
    def fit(self, X, y):
        self.mean_ = y.mean()
        return self

    def predict(self, X):
        return np.full(X.shape[0], self.mean_)


class MajorityClassBaseline:
    def fit(self, X, y):
        self.majority_ = pd.Series(y).value_counts().index[0]
        return self

    def predict(self, X):
        n = X.shape[0] if hasattr(X, "shape") else len(X)
        return np.full(n, self.majority_)


# ── Plotting ──────────────────────────────────────────────────────────────


def plot_predictions_vs_actual(results, save_path):
    n = len(results)
    cols = min(3, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4.5))
    axes = axes.flatten() if n > 1 else [axes]

    all_true = np.concatenate([r["y_true"] for r in results])
    global_min, global_max = all_true.min(), all_true.max()
    lims = (global_min - 1, global_max + 1)

    for i, r in enumerate(results):
        ax = axes[i]
        ax.scatter(r["y_true"], r["y_pred"], alpha=0.5, s=15, c="steelblue", edgecolors="none")
        ax.plot(lims, lims, "r--", lw=1, alpha=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("True n_taxids")
        ax.set_ylabel("Predicted n_taxids")
        ax.set_title(f"{r['name']}  (R²={r['R2']:.3f})")
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Predicted vs Actual — input taxon count", fontsize=13)
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


def plot_residuals(results, save_path):
    n = len(results)
    cols = min(3, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4.5, rows * 3.5))
    axes = axes.flatten() if n > 1 else [axes]

    for i, r in enumerate(results):
        ax = axes[i]
        res = r["y_true"] - r["y_pred"]
        ax.hist(res, bins=25, color="steelblue", edgecolor="black", alpha=0.7)
        ax.axvline(0, color="red", ls="--", lw=1)
        ax.set_xlabel("Residual (true − pred)")
        ax.set_ylabel("Count")
        ax.set_title(f"{r['name']}  (MAE={r['MAE']:.2f})")
        ax.grid(True, alpha=0.3)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Residuals — input taxon count", fontsize=13)
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


def plot_confusion_matrix(results, save_path, label_names):
    n = len(results)
    cols = min(3, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4.5, rows * 4))
    axes = axes.flatten() if n > 1 else [axes]
    int_labels = list(range(len(label_names)))

    for i, r in enumerate(results):
        ax = axes[i]
        cm = confusion_matrix(r["y_true_bin"], r["y_pred_bin"], labels=int_labels)
        sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", ax=ax,
                    xticklabels=label_names, yticklabels=label_names)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        ax.set_title(f"{r['name']}  (Acc={r['Accuracy']:.3f})")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Confusion Matrix — binned taxon count", fontsize=13)
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


def plot_classification_metrics(results, save_path):
    metrics_df = pd.DataFrame([
        {"Model": r["name"], **{k: r[k] for k in ["Accuracy", "MacroF1", "OrdinalMAE"]}}
        for r in results
    ])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, metric in zip(axes, ["Accuracy", "MacroF1", "OrdinalMAE"]):
        vals = metrics_df.set_index("Model")[metric]
        order = vals.sort_values(ascending=False if metric != "OrdinalMAE" else True).index
        sns.barplot(x=vals.values, y=order, ax=ax, palette="viridis")
        ax.set_title(f"{metric}")
        ax.set_xlabel("")
    fig.suptitle("Classification metrics — binned taxon count", fontsize=13)
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


def plot_feature_importance(results, save_path, top_n=15, prefix="regression"):
    tree_results = [r for r in results if r["importance"] is not None]
    if not tree_results:
        return

    n = len(tree_results)
    fig, axes = plt.subplots(1, n, figsize=(n * 5, 5))
    if n == 1:
        axes = [axes]

    for ax, r in zip(axes, tree_results):
        feat_names = r.get("feature_names", [])
        imp = pd.Series(r["importance"], index=feat_names)
        imp = imp.sort_values(ascending=True).tail(min(top_n, len(imp)))
        ax.barh(range(len(imp)), imp.values, color="steelblue", edgecolor="black")
        ax.set_yticks(range(len(imp)))
        ax.set_yticklabels(imp.index, fontsize=8)
        ax.set_xlabel("Importance")
        ax.set_title(r["name"], fontsize=10)
        ax.grid(True, alpha=0.3, axis="x")

    fig.suptitle(f"Top-{top_n} Feature Importances — {prefix}", fontsize=13)
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


def plot_training_time(all_results, save_path):
    names = [r["name"] for r in all_results]
    times = [r["train_time"] for r in all_results]
    colors = ["#e67e22" if t == max(times) else "#bdc3c7" for t in times]

    fig, ax = plt.subplots(figsize=(9, 4))
    bars = ax.barh(names, times, color=colors, edgecolor="black", height=0.6)
    for bar, val in zip(bars, times):
        ax.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height() / 2, f"{val:.2f}s", va="center", fontsize=9)
    ax.set_xlabel("Training Time (seconds)")
    ax.set_title("Training Time — all models")
    ax.grid(True, alpha=0.3, axis="x")
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


# ── Data loading ──────────────────────────────────────────────────────────


def load_data(config, recall_cache=None):
    from deployment.model_evaluation.analysis_scripts.predictor_inputs import (
        collect_training_data,
    )
    data = collect_training_data(config, sort_strategy="reads", recall_cache_path=recall_cache)
    return data


# ── Main ──────────────────────────────────────────────────────────────────


def get_args():
    parser = argparse.ArgumentParser(
        description="Predict input taxon count (continuous + binned) from m_stats features.",
    )
    parser.add_argument("--study_output_filepath", type=str, required=True)
    parser.add_argument("--taxid_plan_filepath", type=str, required=True)
    parser.add_argument("--analysis_output_filepath", type=str, required=True)
    parser.add_argument(
        "--output_dir",
        type=str,
        default="composition_prediction",
        help="Subdirectory within analysis_output_filepath (default: composition_prediction)",
    )
    parser.add_argument("--tax_level", type=str, default="family")
    parser.add_argument("--data_set_divide", type=int, default=16)
    parser.add_argument("--max_training", type=int, default=None)
    parser.add_argument("--recall_cache", type=str, default=None,
                        help="Path to recall_results_cache.parquet (bypasses study output loading)")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main():
    args = get_args()

    import logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    analysis_root = Path(args.analysis_output_filepath)
    output_path = analysis_root / args.output_dir
    output_path.mkdir(parents=True, exist_ok=True)

    config = EvaluatorConfig(
        study_output_filepath=Path(args.study_output_filepath),
        taxid_plan_filepath=Path(args.taxid_plan_filepath),
        analysis_output_filepath=analysis_root,
        data_set_divide=args.data_set_divide,
        tax_level=args.tax_level,
        max_training=args.max_training,
        holdout_proportion=0.3,
    )

    print("=" * 60)
    print("Input Taxonomix Prediction (continuous + binned)")
    print(f"  Output dir: {output_path}")
    print(f"  tax_level:  {config.tax_level}")
    print(f"  Bins:       {list(zip(BIN_LABELS, [(BIN_EDGES[i]+1, BIN_EDGES[i+1]) for i in range(len(BIN_LABELS))]))}")
    print("=" * 60)

    data = load_data(config, recall_cache=args.recall_cache)
    X, y = data["X_train"], data["y_count"]
    X_test, y_test = data["X_test"], data["y_count_test"]

    if X.empty or len(X) < 5:
        print("ERROR: Too few training samples. Check study output data.")
        return

    # Bin target
    label_to_int = {lbl: i for i, lbl in enumerate(BIN_LABELS)}
    y_binned = bin_count(y).map(label_to_int)
    y_test_binned = (bin_count(y_test).map(label_to_int)
                     if not y_test.empty else pd.Series(dtype=float))

    # Split train/val
    X_tr, X_va, y_tr, y_va = train_test_split(
        X, y, test_size=TEST_SIZE, random_state=RANDOM_STATE
    )
    y_tr_binned = bin_count(y_tr).map(label_to_int)
    y_va_binned = bin_count(y_va).map(label_to_int)

    # Scale stats features
    scaler = StandardScaler()
    X_tr_scaled = X_tr.copy()
    X_va_scaled = X_va.copy()
    avail_stats = [c for c in STATS_FEATURES if c in X_tr.columns]
    X_tr_scaled[avail_stats] = scaler.fit_transform(X_tr[avail_stats])
    X_va_scaled[avail_stats] = scaler.transform(X_va[avail_stats])

    if not X_test.empty:
        X_test_scaled = X_test.copy()
        X_test_scaled[avail_stats] = scaler.transform(X_test[avail_stats])

    print(f"Train/Val split: {len(X_tr)} / {len(X_va)}")
    print(f"y range: {y.min():.0f} – {y.max():.0f} (mean={y.mean():.1f}, std={y.std():.1f})")
    print(f"Bin distribution (train):")
    for i, lbl in enumerate(BIN_LABELS):
        cnt = (y_tr_binned == i).sum()
        print(f"  {lbl}: {cnt} ({cnt/len(y_tr_binned)*100:.1f}%)")

    # ── Filter available models ─────────────────────────────────────────

    def _filter_models(cfgs):
        available = []
        for cfg in cfgs:
            if cfg["module"] is None:
                available.append(cfg)
                continue
            try:
                _import_model(cfg)
            except (ImportError, ModuleNotFoundError):
                print(f"  Skipping '{cfg['name']}' (not installed)")
                continue
            available.append(cfg)
        return available

    reg_models = _filter_models(REGRESSION_MODELS_CFG)
    clf_models = _filter_models(CLASSIFICATION_MODELS_CFG)

    bin_labels = BIN_LABELS
    bin_label_set = set(bin_labels)

    # ══════════════════════════════════════════════════════════════════════
    # PART 1: Regression (continuous count)
    # ══════════════════════════════════════════════════════════════════════

    print("\n" + "=" * 60)
    print("PART 1: Continuous Count Regression")
    print("=" * 60)

    reg_results = []

    for cfg in reg_models:
        print(f"\n--- {cfg['name']} ---")
        stats_only = cfg.get("stats_only", False)

        if stats_only:
            X_tr_use = X_tr_scaled[avail_stats]
            X_va_use = X_va_scaled[avail_stats]
        else:
            X_tr_use = X_tr_scaled
            X_va_use = X_va_scaled

        if cfg["module"] is None:
            model = MeanBaseline()
            t0 = time.time()
            model.fit(X_tr_use, y_tr)
            train_time = time.time() - t0
            best_params = None
        elif cfg.get("use_grid"):
            cls = _import_model(cfg)
            gs = GridSearchCV(
                cls(),
                cfg["param_grid"],
                cv=5,
                scoring="neg_mean_absolute_error",
                n_jobs=-1,
                verbose=0,
            )
            t0 = time.time()
            gs.fit(X_tr_use, y_tr)
            train_time = time.time() - t0
            model = gs.best_estimator_
            best_params = gs.best_params_
            print(f"  Best params: {best_params}")
            print(f"  Best CV MAE: {-gs.best_score_:.3f}")
        else:
            cls = _import_model(cfg)
            params = cfg["params"].copy()
            model = cls(**params)
            t0 = time.time()
            model.fit(X_tr_use, y_tr)
            train_time = time.time() - t0
            best_params = None

        y_pred = model.predict(X_va_use)
        metrics = evaluate_regression(y_va.values, y_pred)

        # Test set
        if not X_test.empty and not y_test.empty:
            if stats_only:
                X_test_use = X_test_scaled[avail_stats]
            else:
                X_test_use = X_test_scaled
            y_pred_test = model.predict(X_test_use)
            test_metrics = evaluate_regression(y_test.values, y_pred_test)
            for k, v in test_metrics.items():
                metrics[f"test_{k}"] = v
            print(f"  Test R²={test_metrics['R2']:.4f}, MAE={test_metrics['MAE']:.3f}")
        else:
            test_metrics = None

        if hasattr(model, "feature_importances_"):
            importance = model.feature_importances_
        elif hasattr(model, "coef_"):
            importance = np.abs(model.coef_).flatten() if model.coef_ is not None else None
        else:
            importance = None

        result = {
            "name": cfg["name"],
            **metrics,
            "train_time": train_time,
            "y_true": y_va.values,
            "y_pred": y_pred,
            "importance": importance,
            "feature_names": list(X_tr_use.columns),
            "best_params": best_params,
            "model": model,
            "test_metrics": test_metrics,
            "task": "regression",
        }
        reg_results.append(result)

        print(f"  R²:   {metrics['R2']:.4f}")
        print(f"  MAE:  {metrics['MAE']:.3f}")
        print(f"  RMSE: {metrics['RMSE']:.3f}")
        print(f"  MAPE: {metrics['MAPE']:.2f}%")
        print(f"  Acc±1: {metrics['acc_1']:.1f}%")
        print(f"  Time: {train_time:.2f}s")

    # ══════════════════════════════════════════════════════════════════════
    # PART 2: Classification (binned count)
    # ══════════════════════════════════════════════════════════════════════

    print("\n" + "=" * 60)
    print("PART 2: Binned Count Classification")
    print("=" * 60)

    clf_results = []

    for cfg in clf_models:
        print(f"\n--- {cfg['name']} ---")
        stats_only = cfg.get("stats_only", False)

        if stats_only:
            X_tr_use = X_tr_scaled[avail_stats]
            X_va_use = X_va_scaled[avail_stats]
        else:
            X_tr_use = X_tr_scaled
            X_va_use = X_va_scaled

        if cfg["module"] is None:
            model = MajorityClassBaseline()
            t0 = time.time()
            model.fit(X_tr_use, y_tr_binned)
            train_time = time.time() - t0
            best_params = None
        elif cfg.get("use_grid"):
            cls = _import_model(cfg)
            gs = GridSearchCV(
                cls(),
                cfg["param_grid"],
                cv=5,
                scoring="accuracy",
                n_jobs=-1,
                verbose=0,
            )
            t0 = time.time()
            gs.fit(X_tr_use, y_tr_binned)
            train_time = time.time() - t0
            model = gs.best_estimator_
            best_params = gs.best_params_
            print(f"  Best params: {best_params}")
            print(f"  Best CV Acc: {gs.best_score_:.3f}")
        else:
            cls = _import_model(cfg)
            params = cfg["params"].copy()
            model = cls(**params)
            t0 = time.time()
            model.fit(X_tr_use, y_tr_binned)
            train_time = time.time() - t0
            best_params = None

        y_pred_bin = model.predict(X_va_use)
        y_pred_bin_int = np.asarray(y_pred_bin, dtype=int)
        y_va_bin_int = np.asarray(y_va_binned, dtype=int)

        metrics = evaluate_classification(y_va_bin_int, y_pred_bin_int)

        # Test set evaluation
        test_clf_metrics = None
        if not X_test.empty and not y_test.empty:
            if stats_only:
                X_test_use = X_test_scaled[avail_stats]
            else:
                X_test_use = X_test_scaled
            y_pred_test_bin = model.predict(X_test_use)
            y_pred_test_bin_int = np.asarray(y_pred_test_bin, dtype=int)
            y_test_bin_int = np.asarray(y_test_binned, dtype=int)
            test_clf_metrics = evaluate_classification(y_test_bin_int, y_pred_test_bin_int)
            for k, v in test_clf_metrics.items():
                metrics[f"test_{k}"] = v
            print(f"  Test Acc={test_clf_metrics['Accuracy']:.4f}, OrdinalMAE={test_clf_metrics['OrdinalMAE']:.3f}")

        # Importance
        if hasattr(model, "feature_importances_"):
            importance = model.feature_importances_
        elif hasattr(model, "coef_"):
            importance = np.abs(model.coef_).mean(axis=0) if model.coef_ is not None else None
        else:
            importance = None

        result = {
            "name": cfg["name"],
            **metrics,
            "train_time": train_time,
            "y_true": y_va.values,
            "y_pred": None,
            "y_true_bin": y_va_bin_int,
            "y_pred_bin": y_pred_bin_int,
            "importance": importance,
            "feature_names": list(X_tr_use.columns),
            "best_params": best_params,
            "model": model,
            "test_metrics": test_clf_metrics,
            "task": "classification",
        }
        clf_results.append(result)

        print(f"  Accuracy:    {metrics['Accuracy']:.4f}")
        print(f"  Macro F1:    {metrics['MacroF1']:.4f}")
        print(f"  Ordinal MAE: {metrics['OrdinalMAE']:.3f} (bins off)")
        print(f"  Time:        {train_time:.2f}s")

    # ══════════════════════════════════════════════════════════════════════
    # Save metrics
    # ══════════════════════════════════════════════════════════════════════

    # Regression summary
    reg_summary_rows = []
    for r in reg_results:
        row = {
            "Model": r["name"],
            "R2": round(r["R2"], 4),
            "MAE": round(r["MAE"], 4),
            "RMSE": round(r["RMSE"], 4),
            "MAPE": round(r["MAPE"], 2),
            "Acc_1": round(r["acc_1"], 1),
            "Acc_2": round(r["acc_2"], 1),
            "Acc_5": round(r["acc_5"], 1),
            "Train_Time_s": round(r["train_time"], 2),
        }
        if r.get("test_metrics"):
            row["Test_R2"] = round(r["test_metrics"]["R2"], 4)
            row["Test_MAE"] = round(r["test_metrics"]["MAE"], 4)
        reg_summary_rows.append(row)

    reg_summary_df = pd.DataFrame(reg_summary_rows).sort_values("R2", ascending=False)
    csv_path = output_path / "metrics_regression.csv"
    reg_summary_df.to_csv(csv_path, index=False)
    print(f"\n{'=' * 60}")
    print("Regression Metrics (sorted by R²):")
    print(reg_summary_df.to_string(index=False))

    # Classification summary
    clf_summary_rows = []
    for r in clf_results:
        row = {
            "Model": r["name"],
            "Accuracy": round(r["Accuracy"], 4),
            "BalancedAcc": round(r["BalancedAccuracy"], 4),
            "MacroF1": round(r["MacroF1"], 4),
            "OrdinalMAE": round(r["OrdinalMAE"], 3),
            "OrdinalAcc0": round(r["OrdinalAcc0"], 1),
            "OrdinalAcc1": round(r["OrdinalAcc1"], 1),
            "Train_Time_s": round(r["train_time"], 2),
        }
        if r.get("test_metrics"):
            row["Test_Acc"] = round(r["test_metrics"]["Accuracy"], 4)
            row["Test_OrdinalMAE"] = round(r["test_metrics"]["OrdinalMAE"], 3)
        clf_summary_rows.append(row)

    clf_summary_df = pd.DataFrame(clf_summary_rows).sort_values("Accuracy", ascending=False)
    csv_path = output_path / "metrics_classification.csv"
    clf_summary_df.to_csv(csv_path, index=False)
    print(f"\nClassification Metrics (sorted by Accuracy):")
    print(clf_summary_df.to_string(index=False))

    # Best params
    all_results = reg_results + clf_results
    for r in all_results:
        if r["best_params"]:
            params_path = output_path / "best_params.txt"
            with open(params_path, "w") as f:
                f.write(f"Model: {r['name']}\n")
                f.write(json.dumps(r["best_params"], indent=2))
            print(f"\n  Best params saved to {params_path}")
            break

    # ── Detailed predictions ────────────────────────────────────────────

    detail_rows = []
    for r in reg_results:
        for i in range(len(r["y_true"])):
            detail_rows.append({
                "Model": r["name"],
                "Task": "regression",
                "sample": i,
                "true_n_taxids": int(r["y_true"][i]),
                "pred_n_taxids": round(r["y_pred"][i], 2),
                "residual": r["y_true"][i] - r["y_pred"][i],
                "true_bin": "",
                "pred_bin": "",
            })

    int_to_label = {i: lbl for i, lbl in enumerate(BIN_LABELS)}
    for r in clf_results:
        for i in range(len(r["y_true_bin"])):
            ti = int(r["y_true_bin"][i])
            pi = int(r["y_pred_bin"][i])
            detail_rows.append({
                "Model": r["name"],
                "Task": "classification",
                "sample": i,
                "true_n_taxids": int(r["y_true"][i]) if i < len(r["y_true"]) else "",
                "pred_n_taxids": "",
                "residual": "",
                "true_bin": int_to_label.get(ti, str(ti)),
                "pred_bin": int_to_label.get(pi, str(pi)),
            })

    pd.DataFrame(detail_rows).to_csv(
        output_path / "detailed_predictions.tsv", sep="\t", index=False
    )
    print(f"\n  Detailed predictions saved")

    # ── Plots ───────────────────────────────────────────────────────────

    print("\nGenerating plots...")
    plot_predictions_vs_actual(reg_results, output_path / "predictions_vs_actual.png")
    plot_residuals(reg_results, output_path / "residuals.png")
    plot_confusion_matrix(clf_results, output_path / "confusion_matrix.png", bin_labels)
    plot_classification_metrics(clf_results, output_path / "classification_metrics.png")
    plot_feature_importance(reg_results, output_path / "feature_importance_regression.png", prefix="regression")
    plot_feature_importance(clf_results, output_path / "feature_importance_classification.png", prefix="classification")
    plot_training_time(all_results, output_path / "training_time.png")

    print(f"\n{'=' * 60}")
    print(f"All outputs saved to {output_path}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
