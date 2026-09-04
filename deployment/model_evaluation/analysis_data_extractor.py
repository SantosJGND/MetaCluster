#!/usr/bin/env python3
"""
Extract structured analysis data from simulation outputs for Results section.

Produces per-dataset metrics, recall-by-classifier breakdown, relindex
distribution, and aggregate statistics.

Usage:
    PYTHONPATH=/path/to/project/root python analysis_data_extractor.py \
        --study-output /path/to/study_output \
        --ncbi-db /path/to/taxa.db \
        --output-dir /path/to/output

    Or via main pipeline: PYTHONPATH is set automatically by the runner.
"""

import argparse
import gc
import logging
import os
import resource
import sys
from pathlib import Path

_script_dir = Path(__file__).resolve().parent
_project_root = _script_dir.parent.parent
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

from deployment.model_evaluation.metrics import compute_recall
from deployment.model_evaluation.result_models import write_pipeline_metadata
from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper
from metagenomics_utils.overlap_manager import OverlapManager
from metagenomics_utils.overlap_manager.node_stats import get_m_stats_matrix

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(description="Extract structured analysis data from simulation outputs")
    parser.add_argument("--study-output", type=str, required=True, help="Path to study output directory")
    parser.add_argument("--ncbi-db", type=str, required=True, help="Path to NCBI taxonomy database (taxa.db)")
    parser.add_argument("--output-dir", type=str, required=True, help="Directory to save analysis outputs")
    parser.add_argument("--cross-hit-threshold", type=float, default=0.3, help="Cross-hit threshold (default: 0.3)")
    parser.add_argument("--min-taxonomic-score", type=float, default=0.7, help="Minimum taxonomic score (default: 0.7)")
    parser.add_argument(
        "--explanatory", action="store_true", help="Run explanatory analysis: mixed-effects model on TP read counts"
    )
    return parser.parse_args()


def discover_datasets(study_output_filepath: str) -> list[str]:
    folders = sorted(
        f for f in os.listdir(study_output_filepath) if os.path.isdir(os.path.join(study_output_filepath, f))
    )
    valid = []
    for ds in folders:
        input_path = os.path.join(study_output_filepath, ds, "input", f"{ds}.tsv")
        clustering_path = os.path.join(study_output_filepath, ds, "clustering")
        matched_path = os.path.join(study_output_filepath, ds, "output", "matched_assemblies.tsv")
        if os.path.exists(input_path) and os.path.exists(clustering_path) and os.path.exists(matched_path):
            valid.append(ds)
        else:
            logger.debug(f"Skipping {ds}: missing required files")
    logger.info(f"Discovered {len(valid)} valid datasets out of {len(folders)} total folders")
    return valid


def collect_all_input_taxids(study_output_filepath: str, datasets: list[str]) -> list[int]:
    all_taxids: set[int] = set()
    for ds in datasets:
        input_path = os.path.join(study_output_filepath, ds, "input", f"{ds}.tsv")
        try:
            df = pd.read_csv(input_path, sep="\t")
            all_taxids.update(df["taxid"].dropna().unique().tolist())
        except Exception as e:
            logger.warning(f"Could not read input for {ds}: {e}")
    logger.info(f"Collected {len(all_taxids)} unique input taxids across all datasets")
    return sorted(all_taxids)


def compute_relindex(m_stats: pd.DataFrame) -> float:
    if m_stats.empty:
        return float("nan")
    best_indices = m_stats.index[m_stats["best_match_is_best"] == True].tolist()
    if not best_indices:
        return 0.0
    last_idx = best_indices[-1] + 1
    return last_idx / len(m_stats)


def load_classifier_map(study_output_filepath: str, dataset: str) -> pd.DataFrame:
    cls_path = os.path.join(study_output_filepath, dataset, "classification", f"{dataset}_merged_classification.tsv")
    if not os.path.exists(cls_path):
        logger.warning(f"Classification file not found for {dataset}: {cls_path}")
        return pd.DataFrame({"taxid": pd.Series(dtype=int), "classifiers": pd.Series(dtype=str)})
    df = pd.read_csv(cls_path, sep="\t")
    if "classifiers" in df.columns:
        pass
    elif "classification" in df.columns:
        df = df.rename(columns={"classification": "classifiers"})
    else:
        df["classifiers"] = "unclassified"
    return df[["taxid", "classifiers"]].dropna(subset=["taxid"])


def compute_per_classifier_recall(
    m_stats: pd.DataFrame, classifier_map: pd.DataFrame, input_taxids: set[int]
) -> list[dict]:
    tp_hits = m_stats[m_stats["best_match_is_best"] == True].copy()
    if tp_hits.empty:
        return []

    tp_hits = tp_hits.merge(classifier_map, on="taxid", how="left")
    tp_hits["classifiers"] = tp_hits["classifiers"].fillna("unclassified")

    classifier_taxids: dict[str, set] = {}
    for _, row in tp_hits.iterrows():
        bmt = row.get("best_match_taxid")
        if pd.isna(bmt):
            continue
        tools = str(row["classifiers"]).split("/")
        for tool in tools:
            tool = tool.strip()
            if not tool:
                continue
            if tool not in classifier_taxids:
                classifier_taxids[tool] = set()
            classifier_taxids[tool].add(bmt)

    total_input = len(input_taxids)
    records = []
    for classifier, found_taxids in sorted(classifier_taxids.items()):
        intersection = found_taxids & input_taxids
        recall = len(intersection) / total_input if total_input > 0 else 0.0
        records.append(
            {
                "classifier": classifier,
                "tp_taxids_found": len(intersection),
                "total_input_taxids": total_input,
                "recall": recall,
            }
        )
    return records


def aggregate_input_predictors(input_df: pd.DataFrame) -> pd.DataFrame:
    required = [c for c in ["taxid", "reads"] if c in input_df.columns]
    if len(required) < 2:
        return pd.DataFrame()
    has_mutation = "mutation_rate" in input_df.columns
    agg_dict = {"reads_simulated": ("reads", "sum")}
    if has_mutation:
        agg_dict["mutation_rate"] = ("mutation_rate", "mean")
    result = input_df.groupby("taxid", as_index=False).agg(**agg_dict)
    if not has_mutation:
        result["mutation_rate"] = float("nan")
    return result


def compute_per_taxid_relindex(m_stats: pd.DataFrame) -> dict:
    best = m_stats[m_stats["best_match_is_best"] == True]
    if best.empty:
        return {}
    total = len(m_stats)
    result = {}
    for taxid, grp in best.groupby("best_match_taxid"):
        last_idx = grp.index[-1]
        result[taxid] = {
            "taxid_relindex": (last_idx + 1) / total,
            "taxid_best_count": len(grp),
        }
    return result


def collect_tp_hit_data(
    m_stats: pd.DataFrame,
    input_predictors: pd.DataFrame,
    dataset_name: str,
    per_taxid_relindex: dict = None,
) -> list[dict]:
    tp = m_stats[m_stats["best_match_is_best"] == True].copy()
    if tp.empty or input_predictors.empty:
        return []
    merged = tp.merge(input_predictors, left_on="best_match_taxid", right_on="taxid", how="left")
    if per_taxid_relindex is None:
        per_taxid_relindex = {}
    records = []
    for _, row in merged.iterrows():
        tid = row.get("best_match_taxid")
        rel = per_taxid_relindex.get(tid, {})
        records.append(
            {
                "data_set": dataset_name,
                "taxid": tid,
                "numreads": row.get("numreads", 0),
                "total_uniq_reads": row.get("total_uniq_reads", 0),
                "coverage": row.get("coverage", 0),
                "order": str(row.get("order", "unclassified")),
                "family": str(row.get("family", "unclassified")),
                "genus": str(row.get("genus", "unclassified")),
                "mutation_rate": row.get("mutation_rate", float("nan")),
                "reads_simulated": row.get("reads_simulated", 0),
                "taxid_relindex": rel.get("taxid_relindex", float("nan")),
                "taxid_best_count": rel.get("taxid_best_count", 0),
            }
        )
    return records


def collect_precision_hits(
    m_stats: pd.DataFrame,
    classifier_map: pd.DataFrame,
    dataset_name: str,
) -> list[dict]:
    merged = m_stats.merge(classifier_map, on="taxid", how="left")
    merged["classifiers"] = merged["classifiers"].fillna("unclassified")
    records = []
    for _, row in merged.iterrows():
        tools = str(row["classifiers"]).split("/")
        for tool in tools:
            tool = tool.strip()
            if not tool:
                continue
            records.append(
                {
                    "data_set": dataset_name,
                    "classifier": tool,
                    "best_match_is_best": bool(row["best_match_is_best"]),
                    "is_crosshit": bool(row["is_crosshit"]),
                    "is_trash": bool(row["is_trash"]),
                }
            )
    return records


def process_dataset(
    dataset: str,
    study_output_filepath: str,
    ncbi_wrapper: NCBITaxonomistWrapper,
    cross_hit_threshold: float,
    min_taxonomic_score: float,
    explanatory: bool = False,
) -> tuple[dict, list[dict], list[dict], list[dict], list[dict]]:
    logger.info(f"Processing: {dataset}")

    om_path = os.path.join(study_output_filepath, dataset, "clustering")
    input_path = os.path.join(study_output_filepath, dataset, "input", f"{dataset}.tsv")

    input_df = pd.read_csv(input_path, sep="\t")
    input_taxids = set(input_df["taxid"].dropna().unique())

    overlap_manager = OverlapManager(om_path)
    if not overlap_manager.leaves:
        logger.warning(f"{dataset}: OverlapManager has no leaves, skipping")
        return None, None, None, None, None
    m_stats = get_m_stats_matrix(
        dataset,
        study_output_filepath,
        ncbi_wrapper,
        overlap_manager,
        cross_hit_threshold=cross_hit_threshold,
        min_taxonomic_score=min_taxonomic_score,
        filter_no_leaf=False,
    )

    if m_stats.empty:
        logger.warning(f"{dataset}: m-stats matrix is empty")
        return None, None, None, None, None

    m_stats = m_stats.reset_index()

    metric_relindex = compute_relindex(m_stats)

    tp_count = int((m_stats["best_match_is_best"] == True).sum())
    cross_hit_count = int((m_stats["is_crosshit"] == True).sum())
    spurious_count = int((m_stats["is_trash"] == True).sum())

    recall_raw, recall_cov, _, _ = compute_recall(m_stats, input_df)

    classifier_map = load_classifier_map(study_output_filepath, dataset)
    per_classifier_records = compute_per_classifier_recall(m_stats, classifier_map, input_taxids)

    raw_size = len(m_stats)
    raw_taxids = len(m_stats["taxid"].dropna().unique())

    per_dataset = {
        "data_set": dataset,
        "input_taxid_count": len(input_taxids),
        "output_raw": raw_size,
        "output_taxid_count": raw_taxids,
        "tp_count": tp_count,
        "cross_hit_count": cross_hit_count,
        "spurious_count": spurious_count,
        "overall_recall": recall_raw,
        "recall_cov_filtered": recall_cov,
        "last_best_match_relindex": metric_relindex,
    }

    tp_hit_records = []
    if explanatory:
        input_predictors = aggregate_input_predictors(input_df)
        per_taxid_relindex = compute_per_taxid_relindex(m_stats)
        tp_hit_records = collect_tp_hit_data(m_stats, input_predictors, dataset, per_taxid_relindex)

    classifier_hit_records = collect_precision_hits(m_stats, classifier_map, dataset)

    recall_records = collect_recall_data(m_stats, input_df, dataset, ncbi_wrapper)

    return per_dataset, per_classifier_records, tp_hit_records, classifier_hit_records, recall_records


def collect_recall_data(
    m_stats: pd.DataFrame,
    input_df: pd.DataFrame,
    dataset_name: str,
    ncbi_wrapper: NCBITaxonomistWrapper,
) -> list[dict]:
    tp_taxids: set[int] = set()
    if not m_stats.empty and "best_match_is_best" in m_stats.columns:
        tp_taxids = set(m_stats.loc[m_stats["best_match_is_best"] == True, "best_match_taxid"].dropna().unique())

    records = []
    for _, row in input_df.iterrows():
        tid = row.get("taxid")
        if pd.isna(tid):
            continue
        tid = int(tid)
        order = str(ncbi_wrapper.get_level(tid, "order") or "unclassified")
        records.append(
            {
                "data_set": dataset_name,
                "taxid": tid,
                "reads_simulated": int(row.get("reads", 0)),
                "mutation_rate": float(row.get("mutation_rate", float("nan"))),
                "order": order,
                "recalled": int(tid in tp_taxids),
            }
        )
    return records


# ---------------------------------------------------------------------------
# Explanatory analysis: mixed-effects models for TP read counts
# ---------------------------------------------------------------------------


def fit_tp_mixed_model(
    data: pd.DataFrame,
    target_col: str,
) -> object | None:
    from statsmodels.regression.mixed_linear_model import MixedLM

    df = data.dropna(subset=[target_col, "mutation_rate", "reads_simulated", "coverage", "order", "taxid_relindex"])
    df = df[df["reads_simulated"] > 0].copy()
    if df.empty:
        logger.warning(f"No valid rows for {target_col} model")
        return None

    df["log_target"] = np.log1p(df[target_col])
    df["log_reads"] = np.log1p(df["reads_simulated"])

    n_orders = df["order"].nunique()
    n_datasets = df["data_set"].nunique()

    if n_orders < 2:
        logger.warning(f"Only {n_orders} unique orders for {target_col}, need >= 2")
        return None

    formula = "log_target ~ log_reads + mutation_rate + I(mutation_rate**2) + coverage + taxid_relindex"
    vc_formula = {}
    if n_datasets >= 2:
        vc_formula["dataset"] = "0 + C(data_set)"

    try:
        model = MixedLM.from_formula(formula, groups="order", vc_formula=vc_formula if vc_formula else None, data=df)
        result = model.fit(reml=True, maxiter=200)
        result._fit_data = df
        return result
    except Exception as e:
        print(f"Error fitting mixed model for {target_col}: {e}")
        logger.warning(f"Mixed model failed for {target_col}: {e}")
        return None


def fit_tp_relindex_mixed_model(data: pd.DataFrame) -> object | None:
    from statsmodels.regression.mixed_linear_model import MixedLM

    eps = 1e-6
    target_col = "taxid_relindex"
    df = data.dropna(subset=[target_col, "mutation_rate", "reads_simulated", "coverage", "order"])
    df = df[df["reads_simulated"] > 0].copy()
    if df.empty:
        logger.warning("No valid rows for taxid_relindex model")
        return None

    p = np.clip(df[target_col], eps, 1 - eps)
    df["log_target"] = np.log(p / (1 - p))
    df["log_reads"] = np.log1p(df["reads_simulated"])

    n_orders = df["order"].nunique()
    n_datasets = df["data_set"].nunique()

    if n_orders < 2:
        logger.warning(f"Only {n_orders} unique orders for taxid_relindex, need >= 2")
        return None

    formula = "log_target ~ log_reads + mutation_rate + I(mutation_rate**2) + coverage"
    vc_formula = {}
    if n_datasets >= 2:
        vc_formula["dataset"] = "0 + C(data_set)"

    try:
        model = MixedLM.from_formula(formula, groups="order", vc_formula=vc_formula if vc_formula else None, data=df)
        result = model.fit(reml=True, maxiter=200)
        result._fit_data = df
        from scipy.special import expit

        result._back_transform = expit
        return result
    except Exception as e:
        print(f"Error fitting mixed model for taxid_relindex: {e}")
        logger.warning(f"Mixed model failed for taxid_relindex: {e}")
        return None


def fit_tp_upper_quantile(
    data: pd.DataFrame,
    target_col: str,
    logit: bool = False,
    quantile: float = 0.95,
) -> object | None:
    from statsmodels.regression.quantile_regression import QuantReg

    df = data.dropna(subset=[target_col, "mutation_rate", "reads_simulated", "coverage"])
    df = df[df["reads_simulated"] > 0].copy()
    if df.empty or len(df) < 50:
        logger.warning(f"Too few rows for quantile model ({target_col})")
        return None

    df["log_reads"] = np.log1p(df["reads_simulated"])
    if logit:
        eps = 1e-6
        p = np.clip(df[target_col], eps, 1 - eps)
        df["log_target"] = np.log(p / (1 - p))
    else:
        df["log_target"] = np.log1p(df[target_col])

    formula = "log_target ~ log_reads + mutation_rate + I(mutation_rate**2) + coverage"

    try:
        model = QuantReg.from_formula(formula, data=df)
        result = model.fit(q=quantile, max_iter=2000)
        result._fit_data = df
        result._is_logit = logit
        result._quantile = quantile
        return result
    except Exception as e:
        print(f"Error fitting quantile model for {target_col}: {e}")
        logger.warning(f"Quantile model failed for {target_col}: {e}")
        return None


def save_quantile_summary(result, target_name: str, output_dir: str):
    import scipy.stats as scipy_stats

    tvals = result.params / result.bse
    pvals = 2 * (1 - scipy_stats.norm.cdf(np.abs(tvals)))
    fe = pd.DataFrame(
        {
            "coef": result.params,
            "se": result.bse,
            "tvalue": tvals,
            "p_value": pvals,
        }
    )
    fe["ci_lower"] = fe["coef"] - 1.96 * fe["se"]
    fe["ci_upper"] = fe["coef"] + 1.96 * fe["se"]
    fe.to_csv(os.path.join(output_dir, f"quantile_{target_name}_summary.tsv"), sep="\t")


def plot_quantile_diagnostics(result, target_name: str, output_dir: str):
    df = result._fit_data
    fitted = result.fittedvalues
    observed = df["log_target"]
    residuals = observed - fitted

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].scatter(fitted, residuals, alpha=0.3, s=8)
    axes[0].axhline(0, color="red", linestyle="--")
    axes[0].set_xlabel("Fitted (transformed)")
    axes[0].set_ylabel("Residuals")
    axes[0].set_title(f"Quantile {result._quantile:.0%} — Residuals vs Fitted")

    axes[1].scatter(fitted, observed, alpha=0.3, s=8)
    lo = min(fitted.min(), observed.min())
    hi = max(fitted.max(), observed.max())
    axes[1].plot([lo, hi], [lo, hi], "r--")
    axes[1].set_xlabel("Predicted")
    axes[1].set_ylabel("Observed")
    axes[1].set_title(f"Quantile {result._quantile:.0%} — Observed vs Predicted")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"diagnostic_quantile_{target_name}.png"), dpi=150)
    plt.close()


RECALL_FORMULAS = {
    "A": "recalled ~ log_reads + mutation_rate + I(mutation_rate**2) + C(order)",
    "B": "recalled ~ log_reads * mutation_rate + I(mutation_rate**2) + C(order)",
    "C": "recalled ~ log_reads + mutation_rate + C(order)",
}


def fit_recall_gee(data: pd.DataFrame, formula: str) -> object | None:
    from statsmodels.genmod.cov_struct import Independence
    from statsmodels.genmod.families import Binomial
    from statsmodels.genmod.generalized_estimating_equations import GEE

    df = data.dropna(subset=["reads_simulated", "mutation_rate", "order"]).copy()
    if df.empty or df["recalled"].nunique() < 2:
        return None

    df["log_reads"] = np.log1p(df["reads_simulated"])
    n_datasets = df["data_set"].nunique()
    if n_datasets < 2:
        return None

    try:
        model = GEE.from_formula(formula, groups="data_set", data=df, family=Binomial(), cov_struct=Independence())
        result = model.fit()
        result._fit_data = df
        return result
    except Exception as e:
        logger.warning(f"Recall GEE failed: {e}")
        return None


def save_recall_summary(result, label: str, output_dir: str):
    coef = result.params
    se = result.bse
    tvals = coef / se
    import scipy.stats as scipy_stats

    pvals = 2 * (1 - scipy_stats.norm.cdf(np.abs(tvals)))
    fe = pd.DataFrame(
        {
            "coef": coef,
            "se": se,
            "z": tvals,
            "p_value": pvals,
        }
    )
    fe["ci_lower"] = fe["coef"] - 1.96 * fe["se"]
    fe["ci_upper"] = fe["coef"] + 1.96 * fe["se"]
    fe.to_csv(os.path.join(output_dir, f"recall_variant_{label}_summary.tsv"), sep="\t")


def plot_recall_calibration(results: dict, output_dir: str):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for label, result in results.items():
        df = result._fit_data.copy()
        df["predicted"] = result.fittedvalues
        df["pred_bin"] = pd.cut(df["predicted"], bins=20, include_lowest=True)
        cal = df.groupby("pred_bin", observed=True).agg(
            n=("recalled", "count"),
            obs_rate=("recalled", "mean"),
            pred_rate=("predicted", "mean"),
        )
        axes[0].plot(cal["pred_rate"], cal["obs_rate"], "o-", label=f"Variant {label}", ms=4)
        axes[0].plot([0, 1], [0, 1], "k--", lw=1)
        axes[0].set_xlabel("Predicted probability")
        axes[0].set_ylabel("Observed rate")
        axes[0].set_title("Calibration")
        axes[0].legend()

        from sklearn.metrics import auc, roc_curve

        if df["recalled"].nunique() >= 2:
            fpr, tpr, _ = roc_curve(df["recalled"], df["predicted"])
            roc_auc = auc(fpr, tpr)
            axes[1].plot(fpr, tpr, lw=1.5, label=f"Variant {label} (AUC={roc_auc:.3f})")
    axes[1].plot([0, 1], [0, 1], "k--", lw=1)
    axes[1].set_xlabel("False Positive Rate")
    axes[1].set_ylabel("True Positive Rate")
    axes[1].set_title("ROC Curve")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "recall_calibration.png"), dpi=150)
    plt.close()


def plot_recall_surface(results: dict, output_dir: str):
    result = list(results.values())[0]
    df = result._fit_data
    mr_grid = np.linspace(df["mutation_rate"].min(), df["mutation_rate"].max(), 50)
    lr_grid = np.linspace(df["log_reads"].min(), df["log_reads"].max(), 50)
    mr_mesh, lr_mesh = np.meshgrid(mr_grid, lr_grid)

    order_levels = df["order"].value_counts().index.tolist()
    ref_order = order_levels[0] if order_levels else "unclassified"

    grid_df = pd.DataFrame(
        {
            "log_reads": lr_mesh.ravel(),
            "mutation_rate": mr_mesh.ravel(),
            "order": ref_order,
        }
    )

    for col in df.columns:
        if col not in grid_df.columns and col != "log_reads":
            if pd.api.types.is_numeric_dtype(df[col]) and col not in ["recalled", "data_set", "taxid"]:
                grid_df[col] = df[col].median()

    for label, result in results.items():
        try:
            grid_df["pred"] = result.predict(grid_df)
            prob = grid_df["pred"].values.reshape(50, 50)
        except Exception:
            prob = np.zeros((50, 50))

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.contourf(lr_grid, mr_grid, prob.T, levels=20, cmap="viridis")
    plt.colorbar(im, ax=ax, label="Predicted recall probability")
    ax.set_xlabel("log(reads_simulated + 1)")
    ax.set_ylabel("Mutation rate")
    ax.set_title(f"Recall probability (order={ref_order})")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "recall_probability_surface.png"), dpi=150)
    plt.close()


def save_model_summary(result, target_name: str, output_dir: str):
    fe = pd.DataFrame(
        {
            "coef": result.fe_params,
            "se": result.bse_fe,
            "z": result.tvalues,
            "p_value": result.pvalues,
        }
    )
    fe["ci_lower"] = fe["coef"] - 1.96 * fe["se"]
    fe["ci_upper"] = fe["coef"] + 1.96 * fe["se"]
    fe.to_csv(os.path.join(output_dir, f"model_{target_name}_summary.tsv"), sep="\t")

    try:
        re = result.random_effects
        rows = []
        for k, v in re.items():
            val = v[0] if isinstance(v, (np.ndarray, list)) else v
            rows.append({"group": str(k), "intercept": float(val)})
        pd.DataFrame(rows).to_csv(
            os.path.join(output_dir, f"model_{target_name}_random_effects.tsv"), sep="\t", index=False
        )
    except Exception:
        pass

    try:
        vc = result.variance_components
        if vc:
            pd.DataFrame([{"component": k, "variance": v} for k, v in vc.items()]).to_csv(
                os.path.join(output_dir, f"model_{target_name}_variance_components.tsv"), sep="\t", index=False
            )
    except Exception:
        pass


def plot_tp_diagnostics(result, target_name: str, output_dir: str):
    df = result._fit_data
    fitted = result.fittedvalues
    residuals = result.resid
    back_transform = getattr(result, "_back_transform", np.expm1)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    axes[0, 0].scatter(fitted, residuals, alpha=0.3, s=8)
    axes[0, 0].axhline(0, color="red", linestyle="--", linewidth=1)
    axes[0, 0].set_xlabel("Fitted")
    axes[0, 0].set_ylabel("Residuals")
    axes[0, 0].set_title("Residuals vs Fitted")

    from scipy import stats as scipy_stats

    scipy_stats.probplot(residuals, dist="norm", plot=axes[0, 1])
    axes[0, 1].set_title("Q-Q Plot")

    axes[1, 0].scatter(fitted, back_transform(df["log_target"]), alpha=0.3, s=8)
    lo = min(fitted.min(), df["log_target"].min())
    hi = max(fitted.max(), df["log_target"].max())
    axes[1, 0].plot([lo, hi], back_transform([lo, hi]), color="red", linestyle="--", linewidth=1)
    axes[1, 0].set_xlabel("Predicted (transformed scale)")
    axes[1, 0].set_ylabel("Observed")
    axes[1, 0].set_title("Observed vs Predicted")

    try:
        re = result.random_effects
        order_names = sorted(re.keys(), key=lambda k: re[k][0] if isinstance(re[k], (np.ndarray, list)) else re[k])
        order_vals = [re[o][0] if isinstance(re[o], (np.ndarray, list)) else re[o] for o in order_names]
        axes[1, 1].scatter(order_vals, range(len(order_vals)), alpha=0.7)
        axes[1, 1].axvline(0, color="red", linestyle="--", linewidth=1)
        axes[1, 1].set_xlabel("Random intercept")
        axes[1, 1].set_ylabel("Order (sorted by effect)")
        axes[1, 1].set_title(f"Random Effects — Order (n={len(order_names)})")
    except Exception:
        pass

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"diagnostic_{target_name}.png"), dpi=150)
    plt.close()


def plot_tp_scatters(data: pd.DataFrame, output_dir: str):
    cols = [
        c
        for c in ["numreads", "total_uniq_reads", "reads_simulated", "mutation_rate", "coverage", "taxid_relindex"]
        if c in data.columns
    ]
    if len(cols) < 3:
        return
    plot_df = data[cols + ["order"]].dropna().copy()
    plot_df = plot_df[plot_df["order"] != "unclassified"]
    top = plot_df["order"].value_counts().nlargest(8).index
    plot_df = plot_df[plot_df["order"].isin(top)]
    if plot_df.empty:
        return

    g = sns.PairGrid(plot_df, vars=cols, hue="order", palette="Set2")
    g.map_diag(sns.histplot)
    g.map_offdiag(sns.scatterplot, alpha=0.4, s=6)
    g.add_legend()
    g.savefig(os.path.join(output_dir, "scatter_predictors_vs_target.png"), dpi=150)
    plt.close("all")


# ---------------------------------------------------------------------------


def plot_relindex_distribution(df: pd.DataFrame, output_dir: str):
    data = df["last_best_match_relindex"].dropna()
    if data.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(data, kde=True, bins=30, color="steelblue", edgecolor="white", ax=ax)
    ax.axvline(data.median(), color="red", linestyle="--", label=f"Median: {data.median():.3f}")
    ax.axvline(data.mean(), color="orange", linestyle="--", label=f"Mean: {data.mean():.3f}")
    ax.set_xlabel("Last Best Match Relindex", fontsize=12)
    ax.set_ylabel("Number of Datasets", fontsize=12)
    ax.set_title("Distribution of Last Best Match Relative Index", fontsize=14)
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "relindex_distribution.png"), dpi=150)
    plt.close()


def plot_dataset_size_distribution(df: pd.DataFrame, output_dir: str):
    data = df["output_taxid_count"].dropna()
    if data.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(data, kde=True, bins=30, color="seagreen", edgecolor="white", ax=ax)
    ax.axvline(data.median(), color="red", linestyle="--", label=f"Median: {data.median():.0f}")
    ax.set_xlabel("Raw Output Size (number of leaves)", fontsize=12)
    ax.set_ylabel("Number of Datasets", fontsize=12)
    ax.set_title("Distribution of Raw Dataset Size", fontsize=14)
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "dataset_size_distribution.png"), dpi=150)
    plt.close()


def plot_tp_cross_spurious(df: pd.DataFrame, output_dir: str):
    fig, ax = plt.subplots(figsize=(12, 6))
    plot_df = df.melt(
        id_vars=["data_set"],
        value_vars=["tp_count", "cross_hit_count", "spurious_count"],
        var_name="hit_type",
        value_name="count",
    )
    sns.boxplot(data=plot_df, x="hit_type", y="count", ax=ax, palette="Set2")
    ax.set_xlabel("Hit Type", fontsize=12)
    ax.set_ylabel("Count per Dataset", fontsize=12)
    ax.set_title("Distribution of TP / Cross-Hit / Spurious Counts", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "tp_cross_spurious_distribution.png"), dpi=150)
    plt.close()


def plot_recall_comparison(classifier_df: pd.DataFrame, per_dataset_df: pd.DataFrame, output_dir: str):
    if classifier_df.empty:
        return
    overall = per_dataset_df[["data_set", "overall_recall"]].copy()
    overall = overall.rename(columns={"overall_recall": "recall"})
    overall["classifier"] = "overall"
    combined = pd.concat([classifier_df, overall], ignore_index=True)

    fig, ax = plt.subplots(figsize=(14, 6))
    top_classifiers = classifier_df.groupby("classifier")["recall"].median().index.tolist()
    order = ["overall"] + top_classifiers
    plot_df = combined[combined["classifier"].isin(order)]
    palette = ["#e74c3c" if c == "overall" else "lightblue" for c in order]
    sns.boxplot(data=plot_df, x="classifier", y="recall", ax=ax, order=order, palette=palette)
    ax.set_xlabel("Classifier", fontsize=12)
    ax.set_ylabel("Recall", fontsize=12)
    ax.set_title("Recall: Overall vs Per-Classifier", fontsize=14)
    ax.tick_params(axis="x", rotation=45)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "recall_per_classifier.png"), dpi=150)
    plt.close()


def plot_global_precision_histogram(per_dataset_df: pd.DataFrame, output_dir: str):
    if "tp_count" not in per_dataset_df.columns or "output_taxid_count" not in per_dataset_df.columns:
        return
    df = per_dataset_df.copy()
    df["precision_strict"] = df["tp_count"] / df["output_taxid_count"].clip(lower=1)
    df["precision_approximate"] = (df["tp_count"] + df["cross_hit_count"]) / df["output_taxid_count"].clip(lower=1)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.histplot(
        df["precision_approximate"],
        kde=True,
        bins=30,
        color="seagreen",
        alpha=0.5,
        label="Approximate (incl. CrossHits)",
        ax=ax,
    )
    sns.histplot(
        df["precision_strict"], kde=True, bins=30, color="steelblue", alpha=0.5, label="Strict (Best Match only)", ax=ax
    )
    ax.axvline(df["precision_approximate"].median(), color="seagreen", linestyle="--", linewidth=1.5)
    ax.axvline(df["precision_strict"].median(), color="steelblue", linestyle="--", linewidth=1.5)
    ax.set_xlabel("Precision", fontsize=12)
    ax.set_ylabel("Number of Datasets", fontsize=12)
    ax.set_title("Global Precision Distribution: Strict vs Approximate", fontsize=14)
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "precision_global_distribution.png"), dpi=150)
    plt.close()


def plot_precision_per_classifier(classifier_hit_df: pd.DataFrame, per_dataset_df: pd.DataFrame, output_dir: str):
    if classifier_hit_df.empty:
        return
    per_ds_clf = (
        classifier_hit_df.groupby(["data_set", "classifier"])
        .agg(
            total_hits=("classifier", "count"),
            tp_hits=("best_match_is_best", "sum"),
        )
        .reset_index()
    )
    per_ds_clf["precision"] = per_ds_clf["tp_hits"] / per_ds_clf["total_hits"].clip(lower=1)

    global_df = per_dataset_df.copy()
    global_df["classifier"] = "overall"
    global_df["precision"] = global_df["tp_count"] / global_df["output_taxid_count"].clip(lower=1)
    global_plot = global_df[["data_set", "classifier", "precision"]]

    combined = pd.concat([per_ds_clf[["data_set", "classifier", "precision"]], global_plot], ignore_index=True)

    fig, ax = plt.subplots(figsize=(14, 6))
    top_classifiers = per_ds_clf.groupby("classifier")["precision"].median().nlargest(19).index.tolist()
    order = ["overall"] + top_classifiers
    plot_df = combined[combined["classifier"].isin(order)]
    palette = ["#e74c3c" if c == "overall" else "lightblue" for c in order]
    sns.boxplot(data=plot_df, x="classifier", y="precision", ax=ax, order=order, palette=palette)
    ax.set_xlabel("Classifier", fontsize=12)
    ax.set_ylabel("Precision (Strict)", fontsize=12)
    ax.set_title("Precision: Overall vs Per-Classifier (Top 19)", fontsize=14)
    ax.tick_params(axis="x", rotation=45)
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "precision_per_classifier.png"), dpi=150)
    plt.close()


def main():
    args = parse_args()

    study_output = args.study_output
    ncbi_db = args.ncbi_db
    output_dir = args.output_dir

    os.makedirs(output_dir, exist_ok=True)

    logger.info(f"Study output: {study_output}")
    logger.info(f"NCBI DB: {ncbi_db}")
    logger.info(f"Output dir: {output_dir}")

    all_folders = sorted(
        f for f in os.listdir(study_output) if os.path.isdir(os.path.join(study_output, f))
    )
    total_attempted = len(all_folders)

    datasets = discover_datasets(study_output)
    if not datasets:
        logger.error("No valid datasets found")
        sys.exit(1)

    study_gaps = [f for f in all_folders if f not in datasets]
    if study_gaps:
        logger.warning(f"{len(study_gaps)} folders missing required files (study gaps): {', '.join(study_gaps[:10])}")

    logger.info("Initializing NCBI TaxonomistWrapper and resolving lineages...")
    ncbi_wrapper = NCBITaxonomistWrapper(db=ncbi_db)
    all_taxids = collect_all_input_taxids(study_output, datasets)
    ncbi_wrapper.resolve_lineages(all_taxids)

    per_dataset_records = []
    classifier_records = []
    tp_data_records = []
    classifier_hit_records = []
    recall_records = []
    skipped_names = []
    failed_names = []
    failed_messages = []
    explanatory = args.explanatory
    total_datasets = len(datasets)

    def _write_partial_metadata():
        """Write whatever metadata we have so far — called on interrupt/signal."""
        partial_extracted = len(per_dataset_records)
        partial_skipped = len(skipped_names)
        partial_failed = len(failed_names)
        rows = [
            ("total_attempted", total_attempted),
            ("extracted", partial_extracted),
            ("dropped", total_attempted - partial_extracted),
            ("failed", partial_failed),
            ("skipped", partial_skipped),
        ]
        if skipped_names:
            rows.append(("skipped_datasets", ";".join(skipped_names)))
        if failed_messages:
            rows.append(("failed_datasets", ";".join(failed_messages)))
        if study_gaps:
            rows.append(("study_gaps", len(study_gaps)))
            rows.append(("study_gap_datasets", ";".join(study_gaps)))
        rows.append(("note", "partial — process interrupted before completion"))
        write_pipeline_metadata(rows, output_dir)

    import atexit
    import signal

    _original_sigint = signal.getsignal(signal.SIGINT)
    _original_sigterm = signal.getsignal(signal.SIGTERM)

    def _signal_handler(signum, frame):
        logger.warning(f"Received signal {signum}, writing partial results before exit...")
        _write_partial_metadata()
        signal.signal(signal.SIGINT, _original_sigint)
        signal.signal(signal.SIGTERM, _original_sigterm)
        raise SystemExit(1)

    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)
    atexit.register(_write_partial_metadata)

    for i, ds in enumerate(datasets, 1):
        try:
            result = process_dataset(
                ds,
                study_output,
                ncbi_wrapper,
                args.cross_hit_threshold,
                args.min_taxonomic_score,
                explanatory=explanatory,
            )
            if result[0] is not None:
                per_dataset_records.append(result[0])
                classifier_records.extend(result[1])
                tp_data_records.extend(result[2])
                classifier_hit_records.extend(result[3])
                recall_records.extend(result[4])
            else:
                skipped_names.append(ds)
        except Exception as e:
            logger.warning(f"Error processing {ds}: {e}")
            failed_names.append(ds)
            failed_messages.append(str(e))
        finally:
            gc.collect()

        if i % 10 == 0 or i == total_datasets:
            rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
            logger.info(
                f"[{i}/{total_datasets}] Progress: extracted={len(per_dataset_records)}, "
                f"skipped={len(skipped_names)}, failed={len(failed_names)}, RSS={rss_mb:.0f}MB"
            )

    atexit.unregister(_write_partial_metadata)
    signal.signal(signal.SIGINT, _original_sigint)
    signal.signal(signal.SIGTERM, _original_sigterm)

    extracted = len(per_dataset_records)
    skipped = len(skipped_names)
    failed = len(failed_names)

    logger.info(f"Total attempted: {total_attempted}, Extracted: {extracted}, Skipped: {skipped}, Failed: {failed}, Study gaps: {len(study_gaps)}")

    metadata_rows = [
        ("total_attempted", total_attempted),
        ("extracted", extracted),
        ("dropped", total_attempted - extracted),
        ("failed", failed),
        ("skipped", skipped),
    ]
    if skipped_names:
        metadata_rows.append(("skipped_datasets", ";".join(skipped_names)))
    if failed_messages:
        metadata_rows.append(("failed_datasets", ";".join(failed_messages)))
    if study_gaps:
        metadata_rows.append(("study_gaps", len(study_gaps)))
        metadata_rows.append(("study_gap_datasets", ";".join(study_gaps)))
    
    write_pipeline_metadata(metadata_rows, output_dir)

    if not per_dataset_records:
        logger.error("No datasets were successfully processed")
        sys.exit(1)

    per_dataset_df = pd.DataFrame(per_dataset_records)
    classifier_df = pd.DataFrame(classifier_records)

    per_dataset_path = os.path.join(output_dir, "per_dataset_metrics.tsv")
    per_dataset_df.to_csv(per_dataset_path, sep="\t", index=False)
    logger.info(f"Saved: {per_dataset_path}")

    numeric_cols = per_dataset_df.select_dtypes(include=[np.number]).columns.tolist()
    stats_df = per_dataset_df[numeric_cols].describe().T.reset_index()
    stats_df.columns = ["metric", "count", "mean", "std", "min", "q25", "q50", "q75", "max"]
    stats_path = os.path.join(output_dir, "aggregate_statistics.tsv")
    stats_df.to_csv(stats_path, sep="\t", index=False)
    logger.info(f"Saved: {stats_path}")

    if not classifier_df.empty:
        classifier_path = os.path.join(output_dir, "recall_per_classifier.tsv")
        classifier_df.to_csv(classifier_path, sep="\t", index=False)
        logger.info(f"Saved: {classifier_path}")

        agg_classifier = (
            classifier_df.groupby("classifier")
            .agg(
                dataset_count=("recall", "count"),
                mean_recall=("recall", "mean"),
                median_recall=("recall", "median"),
                std_recall=("recall", "std"),
                min_recall=("recall", "min"),
                max_recall=("recall", "max"),
            )
            .reset_index()
        )
        agg_path = os.path.join(output_dir, "recall_per_classifier_aggregate.tsv")
        agg_classifier.to_csv(agg_path, sep="\t", index=False)
        logger.info(f"Saved: {agg_path}")

        plot_recall_comparison(classifier_df, per_dataset_df, output_dir)

    plot_relindex_distribution(per_dataset_df, output_dir)
    plot_dataset_size_distribution(per_dataset_df, output_dir)
    plot_tp_cross_spurious(per_dataset_df, output_dir)

    plot_global_precision_histogram(per_dataset_df, output_dir)

    classifier_hit_df = pd.DataFrame(classifier_hit_records)
    if not classifier_hit_df.empty:
        plot_precision_per_classifier(classifier_hit_df, per_dataset_df, output_dir)
        clf_hit_path = os.path.join(output_dir, "classifier_hit_counts.tsv")
        per_ds_clf = (
            classifier_hit_df.groupby(["data_set", "classifier"])
            .agg(
                total_hits=("classifier", "count"),
                tp_hits=("best_match_is_best", "sum"),
                cross_hits=("is_crosshit", "sum"),
                trash=("is_trash", "sum"),
            )
            .reset_index()
        )
        per_ds_clf.to_csv(clf_hit_path, sep="\t", index=False)
        logger.info(f"Saved: {clf_hit_path}")

    logger.info(f"All outputs saved to {output_dir}")

    # ---- explanatory analysis ----
    if explanatory:
        try:
            import statsmodels.api as sm
            from statsmodels.regression.mixed_linear_model import MixedLM
        except ImportError:
            logger.error("statsmodels is required for explanatory analysis. Install: pip install statsmodels")
            sys.exit(1)

        exp_dir = os.path.join(output_dir, "explanatory")
        os.makedirs(exp_dir, exist_ok=True)

        if not tp_data_records:
            logger.warning("No TP hit data collected for explanatory analysis")
        else:
            tp_data = pd.DataFrame(tp_data_records)
            tp_path = os.path.join(exp_dir, "tp_hit_data.tsv")
            tp_data.to_csv(tp_path, sep="\t", index=False)
            logger.info(f"Saved: {tp_path} ({len(tp_data)} TP rows)")

            for target in ["numreads", "total_uniq_reads"]:
                logger.info(f"Fitting mixed model for {target}...")
                result = fit_tp_mixed_model(tp_data, target)
                if result is not None:
                    save_model_summary(result, target, exp_dir)
                    plot_tp_diagnostics(result, target, exp_dir)
                    logger.info(f"  Model converged for {target}")

            logger.info("Fitting mixed model for taxid_relindex...")
            relindex_result = fit_tp_relindex_mixed_model(tp_data)
            if relindex_result is not None:
                save_model_summary(relindex_result, "taxid_relindex", exp_dir)
                plot_tp_diagnostics(relindex_result, "taxid_relindex", exp_dir)
                logger.info("  Model converged for taxid_relindex")

            for target in ["numreads", "total_uniq_reads", "taxid_relindex"]:
                logit_flag = target == "taxid_relindex"
                logger.info(f"Fitting upper quantile model for {target}...")
                q_result = fit_tp_upper_quantile(tp_data, target, logit=logit_flag)
                if q_result is not None:
                    save_quantile_summary(q_result, target, exp_dir)
                    plot_quantile_diagnostics(q_result, target, exp_dir)
                    logger.info(f"  Quantile model converged for {target}")

            plot_tp_scatters(tp_data, exp_dir)

            logger.info(f"Explanatory outputs saved to {exp_dir}")

        # ---- recall probability models (GEE) ----
        if not recall_records:
            logger.warning("No recall data collected")
        else:
            recall_dir = os.path.join(exp_dir, "recall_model")
            os.makedirs(recall_dir, exist_ok=True)
            recall_df = pd.DataFrame(recall_records)
            recall_path = os.path.join(recall_dir, "recall_data.tsv")
            recall_df.to_csv(recall_path, sep="\t", index=False)
            logger.info(f"Saved: {recall_path} ({len(recall_df)} rows)")

            recall_results = {}
            for label, formula in RECALL_FORMULAS.items():
                logger.info(f"Fitting recall GEE variant {label}...")
                result = fit_recall_gee(recall_df, formula)
                if result is not None:
                    recall_results[label] = result
                    save_recall_summary(result, label, recall_dir)
                    logger.info(f"  Recall GEE variant {label} converged")

            if recall_results:
                plot_recall_calibration(recall_results, recall_dir)
                plot_recall_surface(recall_results, recall_dir)
                logger.info(f"Recall model outputs saved to {recall_dir}")

        logger.info(f"Explanatory outputs saved to {exp_dir}")


if __name__ == "__main__":
    main()
