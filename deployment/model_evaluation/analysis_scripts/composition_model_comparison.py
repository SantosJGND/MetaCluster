"""
Experimental comparison of multiple models for the clade-composition
(stop_traversal) prediction task.

Loads the cached training data produced by
``data_set_traversal_with_precision()``, applies the same preprocessing as
``CompositionModeller``, then trains and evaluates several classifiers
side-by-side.

Model families compared:
  - FixedThreshold(min_shared>=0.6)  -- rule-based baseline (no training)
  - XGBoost (fixed hparams)
  - Random Forest
  - Gradient Boosting (sklearn)
  - Logistic Regression (stats-only features)
Optional (skipped if unavailable):
  - XGBoost + Optuna  (replicates the production ``ClusteringPipeline``)
  - LightGBM

Outputs (saved to ``--output_dir``):
  - ``comparison_results.csv``     -- side-by-side metrics table
  - ``roc_curves.png``             -- ROC curves, all models
  - ``pr_curves.png``              -- Precision-Recall curves
  - ``f1_comparison.png``          -- F1-score bar chart
  - ``precision_recall_bars.png``  -- Precision / Recall / F1 grouped bars
  - ``confusion_matrices.png``     -- Grid of normalised confusion matrices
  - ``feature_importance.png``     -- Top-15 importance per tree model
  - ``training_time.png``          -- Training-time comparison
  - ``classification_report.txt``  -- Full sklearn report per model
"""

import argparse
import os
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, average_precision_score, roc_curve, precision_recall_curve,
    classification_report, confusion_matrix, ConfusionMatrixDisplay,
)

warnings.filterwarnings("ignore")


class FixedThresholdClustering:
    """Baseline: stop_traversal = (Min_Shared >= 0.6).
    Mirrors the decision rule in traversal_with_clustering_fixed.
    """
    def fit(self, X, y):
        return self
    def predict(self, X):
        return (X['Min_Shared'] >= 0.6).astype(int)
    def predict_proba(self, X):
        preds = self.predict(X)
        proba = np.zeros((len(preds), 2))
        proba[np.arange(len(preds)), preds] = 1.0
        return proba


###############################################################################
# CONSTANTS
###############################################################################

RANDOM_STATE = 42
TEST_SIZE = 0.2

STATS_FEATURES = ["n_leaves", "tax_diversity", "Min_Dist", "Min_Shared"]
DROPPED_COLS = [
    "data_set", "node", "n_true_leaves", "precision",
    "new_precision", "precision_increased", "stop_traversal", "unclassified",
]

MODELS_CFG = [
    {
        "name": "FixedThreshold(min_shared>=0.6)",
        "module": "__main__",
        "class_name": "FixedThresholdClustering",
        "use_optuna": False,
        "params": {},
    },
    {
        "name": "XGBoost+Optuna",
        "module": "xgboost",
        "class_name": "XGBClassifier",
        "use_optuna": True,
        "params": {},
    },
    {
        "name": "XGBoost",
        "module": "xgboost",
        "class_name": "XGBClassifier",
        "use_optuna": False,
        "params": {
            "n_estimators": 300,
            "max_depth": 6,
            "learning_rate": 0.1,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
            "random_state": RANDOM_STATE,
            "eval_metric": "logloss",
        },
    },
    {
        "name": "RandomForest",
        "module": "sklearn.ensemble",
        "class_name": "RandomForestClassifier",
        "use_optuna": False,
        "params": {
            "n_estimators": 300,
            "max_depth": 12,
            "min_samples_leaf": 3,
            "class_weight": "balanced",
            "random_state": RANDOM_STATE,
            "n_jobs": -1,
        },
    },
    {
        "name": "GradientBoosting",
        "module": "sklearn.ensemble",
        "class_name": "GradientBoostingClassifier",
        "use_optuna": False,
        "params": {
            "n_estimators": 300,
            "max_depth": 5,
            "learning_rate": 0.1,
            "subsample": 0.8,
            "random_state": RANDOM_STATE,
        },
    },
    {
        "name": "LogisticRegression",
        "module": "sklearn.linear_model",
        "class_name": "LogisticRegression",
        "use_optuna": False,
        "stats_only": True,
        "params": {
            "C": 1.0,
            "class_weight": "balanced",
            "random_state": RANDOM_STATE,
            "max_iter": 1000,
        },
    },
    {
        "name": "LightGBM",
        "module": "lightgbm",
        "class_name": "LGBMClassifier",
        "use_optuna": False,
        "optional": True,
        "params": {
            "n_estimators": 300,
            "max_depth": 6,
            "learning_rate": 0.1,
            "class_weight": "balanced",
            "random_state": RANDOM_STATE,
            "verbose": -1,
        },
    },
]


###############################################################################
# HELPERS
###############################################################################

def _import_model(cfg):
    """Dynamically import a model class from a (module, class_name) spec."""
    mod = __import__(cfg["module"], fromlist=[cfg["class_name"]])
    return getattr(mod, cfg["class_name"])


def _optuna_trial(trial, X, y, pos_weight):
    """Optuna objective for XGBoost hyperparameter search."""
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
        "random_state": RANDOM_STATE,
        "n_jobs": -1,
    }
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    from xgboost import XGBClassifier
    model = XGBClassifier(**params)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    scores = cross_val_score(model, X, y, cv=cv, scoring="roc_auc", n_jobs=-1)
    return scores.mean()


def _train_with_optuna(X_train, y_train, pos_weight, n_trials=50):
    """Run Optuna optimisation and return best model."""
    import optuna
    study = optuna.create_study(direction="maximize")
    study.optimize(
        lambda trial: _optuna_trial(trial, X_train, y_train, pos_weight),
        n_trials=n_trials,
        show_progress_bar=True,
    )
    from xgboost import XGBClassifier
    best = XGBClassifier(**study.best_trial.params)
    best.fit(X_train, y_train)
    return best, study


def _safe_metric(y_true, y_pred, y_prob, fn, name="", **kw):
    """Compute metric; return NaN on failure (e.g., only one class in batch)."""
    try:
        if y_prob is not None and name in ("roc_auc", "average_precision"):
            return fn(y_true, y_prob, **kw)
        return fn(y_true, y_pred, **kw)
    except Exception:
        return float("nan")


def evaluate_model(y_true, y_pred, y_prob):
    """Return a dict of classification metrics."""
    return {
        "accuracy": _safe_metric(y_true, y_pred, None, accuracy_score),
        "precision": _safe_metric(y_true, y_pred, None, precision_score, zero_division=0),
        "recall": _safe_metric(y_true, y_pred, None, recall_score, zero_division=0),
        "f1": _safe_metric(y_true, y_pred, None, f1_score, zero_division=0),
        "roc_auc": _safe_metric(y_true, y_pred, y_prob, roc_auc_score, name="roc_auc"),
        "pr_auc": _safe_metric(y_true, y_pred, y_prob, average_precision_score, name="average_precision"),
    }


###############################################################################
# PLOTTING
###############################################################################

def plot_roc_curves(results, save_path):
    """Plot ROC curves for all models on a single axis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 7))
    for r in results:
        if r["fpr"] is None or r["tpr"] is None:
            continue
        ax.plot(r["fpr"], r["tpr"], lw=1.5, label=f"{r['name']} (AUC={r['roc_auc']:.3f})")
    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.4)
    ax.set_ylabel("True Positive Rate")
    ax.set_xlabel("False Positive Rate")
    ax.set_title("ROC Curves -- stop_traversal prediction")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


def plot_pr_curves(results, save_path):
    """Plot Precision-Recall curves for all models."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 7))
    for r in results:
        if r["precision_curve"] is None or r["recall_curve"] is None:
            continue
        ax.plot(r["recall_curve"], r["precision_curve"], lw=1.5,
                label=f"{r['name']} (AUC={r['pr_auc']:.3f})")
    ax.axhline(y=results[0].get("pos_ratio", 0.35), color="gray", ls="--", lw=1, alpha=0.5,
               label=f"Prevalence ({results[0].get('pos_ratio', 0.35):.2f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision-Recall Curves -- stop_traversal prediction")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


def plot_f1_comparison(results, save_path):
    """Bar chart of F1 scores."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = [r["name"] for r in results]
    f1_vals = [r["f1"] for r in results]
    colors = ["#2ecc71" if v >= max(f1_vals) else "#95a5a6" for v in f1_vals]

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.barh(names, f1_vals, color=colors, edgecolor="black", height=0.6)
    for bar, val in zip(bars, f1_vals):
        ax.text(bar.get_width() + 0.005, bar.get_y() + bar.get_height() / 2,
                f"{val:.3f}", va="center", fontsize=9)
    ax.set_xlabel("F1 Score")
    ax.set_title("F1 Score (positive class = stop_traversal)")
    ax.set_xlim(0, 1.0)
    ax.grid(True, alpha=0.3, axis="x")
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


def plot_precision_recall_bars(results, save_path):
    """Grouped bar chart of Precision / Recall / F1 per model."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = [r["name"] for r in results]
    x = np.arange(len(names))
    width = 0.25

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x - width, [r["precision"] for r in results], width, label="Precision", color="#3498db")
    ax.bar(x, [r["recall"] for r in results], width, label="Recall", color="#e74c3c")
    ax.bar(x + width, [r["f1"] for r in results], width, label="F1", color="#2ecc71")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Score")
    ax.set_ylim(0, 1.0)
    ax.set_title("Precision / Recall / F1 by Model")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


def plot_confusion_matrices(results, save_path):
    """Grid of normalised confusion matrices."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(results)
    cols = min(3, n)
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3.5, rows * 3.5))
    axes = axes.flatten() if n > 1 else [axes]

    for i, r in enumerate(results):
        cm = confusion_matrix(r["y_true"], r["y_pred"])
        cm_norm = cm.astype("float") / cm.sum(axis=1, keepdims=True)
        disp = ConfusionMatrixDisplay(cm_norm, display_labels=["Continue", "Stop"])
        disp.plot(ax=axes[i], cmap="Blues", values_format=".2f", colorbar=False)
        axes[i].set_title(r["name"], fontsize=10)
        axes[i].set_xlabel("Predicted")
        axes[i].set_ylabel("True")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Normalised Confusion Matrices", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


def plot_feature_importance(results, save_path, top_n=15):
    """Top-N feature importance for each model with feature importances."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

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

    fig.suptitle(f"Top-{top_n} Feature Importances", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


def plot_training_time(results, save_path):
    """Bar chart of training time per model."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = [r["name"] for r in results]
    times = [r["train_time"] for r in results]
    colors = ["#e67e22" if t == max(times) else "#bdc3c7" for t in times]

    fig, ax = plt.subplots(figsize=(9, 4))
    bars = ax.barh(names, times, color=colors, edgecolor="black", height=0.6)
    for bar, val in zip(bars, times):
        ax.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height() / 2,
                f"{val:.2f}s", va="center", fontsize=9)
    ax.set_xlabel("Training Time (seconds)")
    ax.set_title("Training Time by Model")
    ax.grid(True, alpha=0.3, axis="x")
    plt.tight_layout()
    fig.savefig(save_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {save_path}")


###############################################################################
# MAIN
###############################################################################

def get_args():
    parser = argparse.ArgumentParser(
        description="Compare composition (stop_traversal) models side-by-side.",
    )
    parser.add_argument("--output_dir", type=str,
                        default="composition_comparison_outputs",
                        help="Output directory for results "
                             "(default: composition_comparison_outputs)")
    parser.add_argument("--input_cache", type=str,
                        default="deployment/ml_api/training_results_cache.parquet",
                        help="Path to training data cache parquet file "
                             "(default: deployment/ml_api/training_results_cache.parquet)")
    return parser.parse_args()


def main():
    args = get_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cache_path = Path(args.input_cache)
    if not cache_path.exists():
        print(f"ERROR: cache not found at {cache_path}")
        print("Run the evaluation pipeline first to generate training data.")
        return

    print("=" * 60)
    print("Composition Model Comparison")
    print(f"  Output dir: {output_dir}")
    print(f"  Input cache: {cache_path}")
    print("=" * 60)

    df = pd.read_parquet(cache_path)
    print(f"  Shape: {df.shape}")
    print(f"  Columns: {list(df.columns)}")

    missing = [c for c in DROPPED_COLS if c not in df.columns]
    if missing:
        print(f"  WARNING: expected columns not found: {missing}")

    available_drop = [c for c in DROPPED_COLS if c in df.columns]
    X_raw = df.drop(columns=available_drop)
    y = df["stop_traversal"].astype(int)

    tax_cols = [c for c in X_raw.columns if c not in STATS_FEATURES]

    print(f"\nFeature breakdown:")
    print(f"  Stats cols ({len(STATS_FEATURES)}): {STATS_FEATURES}")
    print(f"  Taxonomy cols ({len(tax_cols)}): {tax_cols}")
    print(f"  Target distribution:\n    {y.value_counts().to_dict()}")

    X_train, X_test, y_train, y_test = train_test_split(
        X_raw, y, test_size=TEST_SIZE, random_state=RANDOM_STATE,
        stratify=y,
    )
    print(f"\nTrain/Test split: {len(X_train)} / {len(X_test)}")

    scaler = StandardScaler()
    X_train_scaled = X_train.copy()
    X_test_scaled = X_test.copy()
    X_train_scaled[STATS_FEATURES] = scaler.fit_transform(X_train[STATS_FEATURES])
    X_test_scaled[STATS_FEATURES] = scaler.transform(X_test[STATS_FEATURES])

    n_neg = (y_train == 0).sum()
    n_pos = (y_train == 1).sum()
    pos_weight = n_neg / max(n_pos, 1)
    pos_ratio_test = y_test.mean()
    print(f"  Positive weight: {pos_weight:.2f}  (train prevalence: {y_train.mean():.3f})")

    available_models = []
    for cfg in MODELS_CFG:
        try:
            _import_model(cfg)
        except (ImportError, ModuleNotFoundError):
            msg = f"  Skipping '{cfg['name']}' (not installed)"
            if cfg.get("optional"):
                print(msg)
            else:
                print(msg + " -- required models may be missing")
            continue
        if cfg.get("use_optuna"):
            try:
                import optuna  # noqa: F401
            except ImportError:
                print(f"  Skipping '{cfg['name']}' (optuna not available)")
                continue
        available_models.append(cfg)

    all_results = []

    for cfg in available_models:
        print(f"\n--- {cfg['name']} ---")
        cls = _import_model(cfg)
        stats_only = cfg.get("stats_only", False)

        if stats_only:
            X_tr = X_train_scaled[STATS_FEATURES]
            X_te = X_test_scaled[STATS_FEATURES]
        else:
            X_tr = X_train_scaled
            X_te = X_test_scaled

        if cfg["use_optuna"]:
            t0 = time.time()
            model, study = _train_with_optuna(
                X_tr.values, y_train.values, pos_weight, n_trials=50,
            )
            train_time = time.time() - t0
            print(f"  Best params: {study.best_trial.params}")
            print(f"  Best CV AUC: {study.best_value:.4f}")
        else:
            params = cfg["params"].copy()
            if "scale_pos_weight" in cfg["params"]:
                params.setdefault("scale_pos_weight", pos_weight)
            elif "class_weight" in params and params["class_weight"] == "balanced":
                pass
            model = cls(**params)
            t0 = time.time()
            model.fit(X_tr, y_train)
            train_time = time.time() - t0

        y_pred = model.predict(X_te)
        y_prob = model.predict_proba(X_te)[:, 1] if hasattr(model, "predict_proba") else None

        metrics = evaluate_model(y_test, y_pred, y_prob)

        fpr, tpr, _ = roc_curve(y_test, y_prob) if y_prob is not None else (None, None, None)
        prec_curve, rec_curve, _ = precision_recall_curve(y_test, y_prob) if y_prob is not None else (None, None)

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
            "y_true": y_test,
            "y_pred": y_pred,
            "y_prob": y_prob,
            "fpr": fpr,
            "tpr": tpr,
            "precision_curve": prec_curve,
            "recall_curve": rec_curve,
            "importance": importance,
            "feature_names": list(X_tr.columns),
            "model": model,
            "pos_ratio": float(pos_ratio_test),
        }
        all_results.append(result)

        print(f"  Accuracy:  {metrics['accuracy']:.4f}")
        print(f"  Precision: {metrics['precision']:.4f}")
        print(f"  Recall:    {metrics['recall']:.4f}")
        print(f"  F1:        {metrics['f1']:.4f}")
        print(f"  ROC-AUC:   {metrics['roc_auc']:.4f}")
        print(f"  PR-AUC:    {metrics['pr_auc']:.4f}")
        print(f"  Time:      {train_time:.2f}s")

    summary_rows = []
    for r in all_results:
        summary_rows.append({
            "Model": r["name"],
            "Accuracy": round(r["accuracy"], 4),
            "Precision": round(r["precision"], 4),
            "Recall": round(r["recall"], 4),
            "F1": round(r["f1"], 4),
            "ROC-AUC": round(r["roc_auc"], 4),
            "PR-AUC": round(r["pr_auc"], 4),
            "Train_Time_s": round(r["train_time"], 2),
        })
    summary_df = pd.DataFrame(summary_rows).sort_values("F1", ascending=False)
    csv_path = output_dir / "comparison_results.csv"
    summary_df.to_csv(csv_path, index=False)
    print(f"\n{'=' * 60}")
    print("Comparison Results (sorted by F1):")
    print(f"{'=' * 60}")
    print(summary_df.to_string(index=False))

    report_path = output_dir / "classification_report.txt"
    with open(report_path, "w") as f:
        for r in all_results:
            f.write(f"\n{'=' * 60}\n")
            f.write(f"Model: {r['name']}\n")
            f.write(f"{'=' * 60}\n")
            f.write(classification_report(r["y_true"], r["y_pred"], digits=4))
            f.write(f"ROC-AUC:  {r['roc_auc']:.4f}\n")
            f.write(f"PR-AUC:   {r['pr_auc']:.4f}\n")
            f.write(f"Train time: {r['train_time']:.2f}s\n")
    print(f"\nDetailed report saved to {report_path}")

    print(f"\nGenerating plots...")
    plot_roc_curves(all_results, output_dir / "roc_curves.png")
    plot_pr_curves(all_results, output_dir / "pr_curves.png")
    plot_f1_comparison(all_results, output_dir / "f1_comparison.png")
    plot_precision_recall_bars(all_results, output_dir / "precision_recall_bars.png")
    plot_confusion_matrices(all_results, output_dir / "confusion_matrices.png")
    plot_feature_importance(all_results, output_dir / "feature_importance.png")
    plot_training_time(all_results, output_dir / "training_time.png")

    print(f"\n{'=' * 60}")
    print(f"All outputs saved to {output_dir}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
