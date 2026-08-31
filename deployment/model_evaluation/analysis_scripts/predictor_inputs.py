"""
Shared data-loading helpers for input-target prediction scripts.

Provides:
  - _get_raw_m_stats_matrix   — load raw leaf stats for one dataset
  - _load_input_summary       — load ground-truth input for one dataset
  - reconstruct_recall_feature_frame — rebuild recall feature matrix from a model cache dir
  - collect_training_data     — assemble features, count-target, comp-target
"""

import logging
import os
import warnings
from pathlib import Path

import joblib
import pandas as pd

from deployment.model_evaluation.config import EvaluatorConfig
from deployment.model_evaluation.data_loader import DataLoader
from metagenomics_utils.overlap_manager.feature_transformer import RecallFeatureTransformer
from metagenomics_utils.overlap_manager.node_stats import (
    get_m_stats_matrix,
    get_subset_composition_counts,
)

logger = logging.getLogger(__name__)

STAT_FEATURES = [
    "number_hits",
    "counts_kurtosis",
    "counts_skewness",
    "tax_diversity_shannon",
    "max_uniq_reads",
    "total_uniq_reads",
]

RECALL_DIVISIONS = 20

RECALL_CACHE_FILES = (
    "recall_matrices_cache.joblib",
    "taxids_to_use_cache.parquet",
    "training_results_cache.parquet",
)


def reconstruct_recall_feature_frame(
    cache_dir,
    data_set_divide: int = 20,
    tax_level: str = "order",
    sort_strategy: str = "reads",
) -> pd.DataFrame:
    """
    Rebuild the recall feature matrix from an evaluation model cache directory.

    The current pipeline stores raw per-dataset m_stats matrices in
    ``models/cache/recall_matrices_cache.joblib`` (not the older
    ``recall_results_cache.parquet``). This helper reconstructs the feature
    vectors exactly as ``RecallModeller.fit()`` does via
    ``RecallFeatureTransformer``, producing columns ``index_recall_1..N``,
    ``last_best_match_relindex``, the 6 stat features, and taxonomy
    proportions, plus a ``data_set`` column for N-taxid augmentation.

    Parameters
    ----------
    cache_dir : str | Path
        The model cache directory (e.g. ``<analysis>/models/cache``).
    data_set_divide : int
        Number of recall divisions (produces ``index_recall_1..N``; must be
        at least the active division count used downstream).
    tax_level : str
        Taxonomic level for composition features.
    sort_strategy : str
        Sort strategy for ``RecallFeatureTransformer``.

    Returns
    -------
    pd.DataFrame
        Reconstructed feature + target matrix, one row per training dataset,
        with a ``data_set`` column.
    """
    cache_dir = Path(cache_dir)
    for fname in RECALL_CACHE_FILES:
        if not (cache_dir / fname).exists():
            raise FileNotFoundError(
                f"Cache file not found: {cache_dir / fname}. Expected an evaluation model "
                f"cache directory containing {', '.join(RECALL_CACHE_FILES)}."
            )

    matrices = joblib.load(cache_dir / "recall_matrices_cache.joblib")
    taxids_df = pd.read_parquet(cache_dir / "taxids_to_use_cache.parquet")

    if not matrices:
        raise ValueError("recall_matrices_cache.joblib contains no training matrices.")

    transformer = RecallFeatureTransformer(
        tax_level=tax_level,
        data_set_divide=data_set_divide,
        sort_strategy=sort_strategy,
    )
    transformer.fit(taxids_to_use=taxids_df)
    rows = [transformer.transform(m) for m in matrices]
    frame = pd.concat(rows, ignore_index=True)

    train_names_df = pd.read_parquet(cache_dir / "training_results_cache.parquet")
    ds_names = (
        train_names_df["data_set"].unique()
        if "data_set" in train_names_df.columns
        else []
    )
    if len(ds_names) == len(matrices):
        frame["data_set"] = list(ds_names)
    else:
        frame["data_set"] = [f"dataset_{i:04d}" for i in range(len(matrices))]
        warnings.warn(
            f"Cache dataset-name count ({len(ds_names)}) != matrices ({len(matrices)}); "
            "using positional dataset names. N-taxid augmentation may misalign.",
            stacklevel=2,
        )

    return frame


def _get_raw_m_stats_matrix(
    data_set_name: str,
    study_output_filepath: Path,
    ncbi_wrapper,
) -> pd.DataFrame | None:
    """Load the raw m_stats matrix (pre-sorting) for a single dataset."""
    from metagenomics_utils.overlap_manager import OverlapManager

    om_path = os.path.join(str(study_output_filepath), data_set_name, "clustering")
    try:
        overlap_manager = OverlapManager(om_path)
    except Exception:
        return None

    if overlap_manager.m_stats_matrix.empty:
        return None

    try:
        m_stats = get_m_stats_matrix(
            data_set_name,
            str(study_output_filepath),
            ncbi_wrapper,
            overlap_manager,
            filter_no_leaf=False,
        )
    except Exception as e:
        logger.debug("get_m_stats_matrix failed for %s: %s", data_set_name, e)
        return None

    return m_stats if not m_stats.empty else None


def _load_input_summary(
    data_set_name: str,
    study_output_filepath: Path,
) -> pd.DataFrame | None:
    """Load the input summary (taxid column) for a single dataset."""
    input_path = os.path.join(
        str(study_output_filepath), data_set_name, "input", f"{data_set_name}.tsv"
    )
    if not os.path.exists(input_path):
        return None
    input_df = pd.read_csv(input_path, sep="\t")
    return input_df[["sample", "taxid", "reads", "mutation_rate"]].drop_duplicates()


def _compute_targets(
    dataset_names: list[str],
    config: EvaluatorConfig,
    input_tax_df: pd.DataFrame,
    taxids_to_use: pd.DataFrame,
    ref_taxa: list[str],
):
    """Load input summaries for a list of datasets and compute count + comp targets."""
    count_vals, comp_rows = [], []
    for name in dataset_names:
        inp = _load_input_summary(name, config.study_output_filepath)
        if inp is None or inp.empty:
            logger.debug("Skipping %s (no input summary)", name)
            count_vals.append(None)
            comp_rows.append(None)
            continue
        inp = inp.merge(input_tax_df, on="taxid", how="left")
        count_vals.append(inp["taxid"].nunique())
        comp = get_subset_composition_counts(
            inp, taxids_to_use,
            tax_level=config.tax_level, count_column="reads",
        )
        comp_row = comp.set_index("tax_level")["proportion"]
        comp_rows.append(comp_row)

    return count_vals, comp_rows


def collect_training_data(
    config: EvaluatorConfig,
    sort_strategy: str = "reads",
    data_set_divide: int = 16,
    recall_cache_path: str | Path | None = None,
):
    """
    Assemble feature/target matrices for both prediction tasks.

    Two modes:
      1. *From cache* — if ``recall_cache_path`` points to a model cache
         directory (``<analysis>/models/cache``) or an existing
         ``recall_results_cache.parquet``, pre-computed features are reused
         (a cache dir is reconstructed via ``reconstruct_recall_feature_frame``).
      2. *From study* — loads raw m_stats matrices from study output and
         runs ``RecallFeatureTransformer`` per dataset.

    Returns
    -------
    dict with keys:
        X_train, y_count, y_comp,
        X_test, y_count_test, y_comp_test,
        feat_names, ref_taxa,
        transformer, loader, input_tax_df
    """
    loader = DataLoader(config).initialize()
    ncbi = loader.get_ncbi_wrapper()
    input_tax_df = loader.get_input_tax_df()
    taxids_to_use = loader.get_taxids_to_use()
    train_folders = loader.get_training_folders()
    test_folders = loader.get_test_folders()

    all_train_set = set(train_folders)
    all_test_set = set(test_folders)

    ref_taxa = sorted(taxids_to_use[config.tax_level].dropna().unique())
    logger.info(
        "Reference taxa at '%s': %d (used as composition columns)",
        config.tax_level,
        len(ref_taxa),
    )

    use_cache = recall_cache_path and Path(recall_cache_path).exists()

    if use_cache:
        cache_path = Path(recall_cache_path)
        if cache_path.is_dir():
            logger.info("Loading pre-computed features from cache dir: %s", cache_path)
            cache = reconstruct_recall_feature_frame(
                cache_path,
                data_set_divide=RECALL_DIVISIONS,
                tax_level=config.tax_level,
                sort_strategy=sort_strategy,
            )
        else:
            logger.info("Loading pre-computed features from cache: %s", cache_path)
            cache = pd.read_parquet(cache_path)
        cache_datasets = set(cache["data_set"].unique())

        train_in_cache = all_train_set & cache_datasets
        test_in_cache = all_test_set & cache_datasets

        train_cache = cache[cache["data_set"].isin(train_in_cache)].copy()
        test_cache = cache[cache["data_set"].isin(test_in_cache)].copy() if test_in_cache else pd.DataFrame()

        # Identify taxonomy columns in cache: columns not in stat/target/data_set
        target_cols = {"data_set", "last_best_match_relindex"} | {
            f"index_recall_{i}" for i in range(1, 21)
        }
        cache_tax_cols = [c for c in cache.columns if c not in set(STAT_FEATURES) | target_cols]

        # Use intersection of ref_taxa with cache taxonomy columns
        common_taxa = [t for t in ref_taxa if t in cache_tax_cols]
        missing_taxa = [t for t in ref_taxa if t not in cache_tax_cols]
        if missing_taxa:
            logger.warning(
                "Cache missing %d reference taxa (will be zero-filled): %s",
                len(missing_taxa), missing_taxa[:5],
            )

        feat_names = STAT_FEATURES + common_taxa

        X_train = train_cache[feat_names].reset_index(drop=True) if not train_cache.empty else pd.DataFrame()
        X_test = test_cache[feat_names].reset_index(drop=True) if not test_cache.empty else pd.DataFrame()

        # Compute targets from input summaries
        train_names = sorted(train_in_cache, key=lambda n: list(train_folders).index(n) if n in train_folders else 0)
        test_names = sorted(test_in_cache, key=lambda n: list(test_folders).index(n) if n in test_folders else 0)

        train_counts, train_comps = _compute_targets(train_names, config, input_tax_df, taxids_to_use, ref_taxa)
        test_counts, test_comps = _compute_targets(test_names, config, input_tax_df, taxids_to_use, ref_taxa)

        # Filter out rows where input loading failed
        valid_train = [(c, p) for c, p in zip(train_counts, train_comps) if c is not None]
        valid_test = [(c, p) for c, p in zip(test_counts, test_comps) if c is not None]

        if valid_train:
            X_train = X_train.iloc[:len(valid_train)].reset_index(drop=True)
            yc_train = pd.Series([v[0] for v in valid_train], name="n_taxids")
            yp_train = pd.DataFrame([v[1] for v in valid_train], columns=ref_taxa).fillna(0.0)
        else:
            X_train, yc_train, yp_train = pd.DataFrame(), pd.Series(name="n_taxids", dtype=float), pd.DataFrame()

        if valid_test:
            X_test = X_test.iloc[:len(valid_test)].reset_index(drop=True)
            yc_test = pd.Series([v[0] for v in valid_test], name="n_taxids")
            yp_test = pd.DataFrame([v[1] for v in valid_test], columns=ref_taxa).fillna(0.0)
        else:
            X_test, yc_test, yp_test = pd.DataFrame(), pd.Series(name="n_taxids", dtype=float), pd.DataFrame()

        transformer = None

        logger.info(
            "Cache: %d train / %d test samples, %d features (stat=%d, tax=%d)",
            len(X_train), len(X_test), len(feat_names), len(STAT_FEATURES), len(common_taxa),
        )

    else:
        # Original approach: load m_stats from study output
        logger.info("Loading raw m_stats matrices from study output")
        transformer = RecallFeatureTransformer(
            tax_level=config.tax_level,
            data_set_divide=data_set_divide,
            sort_strategy=sort_strategy,
        )
        transformer.fit(taxids_to_use=taxids_to_use)
        feat_names = transformer.get_feature_names_out()
        logger.info("Feature dimension: %d (%d stat + %d composition)",
                    len(feat_names), 6, len(feat_names) - 6)

        def _collect(folders):
            X_rows, count_vals, comp_rows = [], [], []
            for name in folders:
                m_stats = _get_raw_m_stats_matrix(name, config.study_output_filepath, ncbi)
                if m_stats is None or m_stats.empty:
                    logger.debug("Skipping %s (no m_stats)", name)
                    continue
                inp = _load_input_summary(name, config.study_output_filepath)
                if inp is None or inp.empty:
                    logger.debug("Skipping %s (no input summary)", name)
                    continue
                inp = inp.merge(input_tax_df, on="taxid", how="left")

                try:
                    feats = transformer.transform(m_stats)
                except Exception as e:
                    logger.warning("Transform failed for %s: %s", name, e)
                    continue
                X_rows.append(feats)
                count_vals.append(inp["taxid"].nunique())

                comp = get_subset_composition_counts(
                    inp, taxids_to_use,
                    tax_level=config.tax_level, count_column="reads",
                )
                comp_row = comp.set_index("tax_level")["proportion"]
                comp_rows.append(comp_row)

            if not X_rows:
                return pd.DataFrame(), pd.Series(name="n_taxids", dtype=float), pd.DataFrame()

            X = pd.concat(X_rows, ignore_index=True)
            yc = pd.Series(count_vals, name="n_taxids")
            yp = pd.DataFrame(comp_rows, columns=ref_taxa).fillna(0.0)
            return X, yc, yp

        X_train, yc_train, yp_train = _collect(train_folders)
        X_test, yc_test, yp_test = _collect(test_folders)

        logger.info(
            "Study: %d train / %d test samples, %d features",
            len(X_train), len(X_test), len(feat_names),
        )

    return dict(
        X_train=X_train,
        y_count=yc_train,
        y_comp=yp_train,
        X_test=X_test,
        y_count_test=yc_test,
        y_comp_test=yp_test,
        feat_names=feat_names,
        ref_taxa=ref_taxa,
        transformer=transformer,
        loader=loader,
        input_tax_df=input_tax_df,
    )
