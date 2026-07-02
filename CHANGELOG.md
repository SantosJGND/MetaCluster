# Changelog

## [1.0.0] - 2026-07-01

### Added

- **`main.nf` entry point** — root-level workflow dispatcher supporting `benchmark` (default) and `classify_only` workflows via the `--workflow` parameter.
- **Root `nextflow.config`** — consolidated DSL2 configuration with conda profile, resource defaults, and process-specific environment mappings, replacing per-module config scatter.
- **Known limitations section** in README — documents the skipped Kraken tree processing test as a known issue.
- **Data availability section** in README — placeholder for publication results archive.
- **Reproducibility section** in README — documents tested software versions and pinned dependencies.

### Changed

- **Version bump** — `0.1.0` → `1.0.0` (first stable release matching publication).
- **CITATION.cff** — updated `date-released` to `2026-07-01`, version to `1.0.0`.
- **`.gitignore`** — added `nextflow_dag.html`, `clade_precision_diagnostics.csv`, `local_evaluation/` to ignore list.

### Removed

- **Tracked artifacts** — `nextflow_dag.html` removed from git tracking.
- **`local_evaluation/`** — empty directory removed.

## [Unreleased]

### Changed

- **`__getstate__`/`__setstate__` on 4 custom model classes** (`RecallModeller`, `GPCLFThreshold`, `ClusteringPipeline`, `BaseCompositionModeller`) — strips training-only attributes (`X_test`, `y_test`, `study_`, training DataFrames, etc.) from serialized pickles; `__setstate__` re-inits defaults. Existing models must be retrained.
- **Lazy imports in `metagenomics_utils`** — `networkx`, `NCBITaxonomistWrapper`, `OverlapManager` moved from module-level to inside training-only functions; `overlap_manager/__init__.py` uses `__getattr__` for `OverlapManager`/`merge_*`; `metagenomics_utils/__init__.py` no longer re-exports `ncbi_tools` or `OverlapManager`. `from __future__ import annotations` added to `node_stats.py` and `om_models.py`.

### Removed

- **`biopython`, `matplotlib`, `optuna` from `ml_api/requirements.txt`** — no longer needed at inference time.
- **Heavy files from Docker image** — `ncbi_tools.py`, `reference_utils.py`, `core.py`, `dataframe_utils.py` deleted from `COPY` in `deployment/ml_api/Dockerfile`.

### Added

- **Per-taxid recall GEE models** (`deployment/model_evaluation/analysis_data_extractor.py`) — `collect_recall_data()` builds per-taxid binary recall records inside `process_dataset()`; `fit_recall_gee()` fits Binomial GEE (`statsmodels.GEE`, `family=Binomial`, `cov_struct=Independence`, `groups=data_set`); `save_recall_summary()`, `plot_recall_calibration()`, `plot_recall_surface()` produce coefficient tables, calibration/ROC curves, and probability surface plots. Three formula variants in `RECALL_FORMULAS`: A (`~ log_reads + mutation_rate + I(mutation_rate²) + C(order)`), B (`~ log_reads * mutation_rate + I(mutation_rate²) + C(order)`), C (`~ log_reads + mutation_rate + C(order)`).
- **`I(mutation_rate²)` quadratic term** — added to all three mixed model formulas (`fit_tp_mixed_model`, `fit_tp_relindex_mixed_model`) and all three quantile formulas (`fit_tp_upper_quantile`).
- **TP relindex mixed model** (`fit_tp_relindex_mixed_model()`) — MixedLM with logit-transformed `taxid_relindex` as target, formula `log_target ~ log_reads + mutation_rate + I(mutation_rate²) + coverage`.
- **Upper quantile regression** (`fit_tp_upper_quantile()`, `save_quantile_summary()`, `plot_quantile_diagnostics()`) — 95th-percentile QuantReg for all three targets with logit (relindex) or log1p (count) transforms.
- **`taxid_relindex` to scatter pairplot** — added to `plot_tp_scatters()` column list.
- **Recall model README** (`explanatory/recall_model/README.md`) — documents data construction, GEE method, coefficient tables for all three variants with order effects, calibration/ROC, probability surface, and file manifest.

### Fixed

- **`collect_recall_data()` silent order-resolution bug** (`analysis_data_extractor.py:337`) — replaced `ncbi_wrapper.get_lineage(tid).get('order', 'unclassified')` (returns a list → `AttributeError` silently caught by bare `except Exception: pass`, defaulting all taxids to `unclassified`) with `ncbi_wrapper.get_level(tid, 'order')`; removed bare `except`. All 40 viral taxids now resolve to 15 named orders via local `taxa.db`.
- **`taxid_best_count` removed from mixed model formulas** — column was constant (=1) across all TP rows, causing `LinAlgError: Singular matrix` in `MixedLM`. Dropped from count and relindex formulas.
- **`predict_composition_stop_traversal`: missing taxonomy features default to 0.0** — was HTTP 422 when request omitted taxonomy columns
- **`predict_televir_clustering_threshold`: `taxon_hits` defaults to 12 zeros** — was `[]` → `KeyError: 'Baculoviridae'`
- **`predict_recall_cutoff_from_table`: `tax_level` column forwarded to transformer** — was `KeyError: 'order'` in `node_composition_with_stats`
- **`ClusteringPipeline._scale`/`predict`/`predict_proba`: backward compat for old pickles** — pickles with `model`/`scaler` attrs (instead of `classifier_`/`scaler_`) now load and predict correctly

### Changed

- **`ml_api` build instructions** — replaced `docker compose up ml_app` (no compose file exists) with `docker build -t ml_api . && docker run -p 8000:8000`
- **`plot_tp_diagnostics()` back-transform** — uses `result._back_transform` (expit for logit, expm1 for log1p) instead of hard-coded `np.expm1`, enabling correct diagnostic plots for the logit-transformed `taxid_relindex` target.
- **`process_dataset()` return type** — 5-tuple `(per_dataset, per_classifier, tp_hit_data, classifier_hits, recall_records)`; skip returns emit 5 `None`s.
- **Recall models run outside `tp_data_records` guard** — GEE models process recall data unconditionally (recall records are collected from input TSVs, independent of TP hit data).

- **`FixedCountModeller` baseline** (`deployment/model_evaluation/compare_sort_strategies.py:144`) — rule-based recall model that keeps leaves until N unique `best_match_taxid` values are covered after strategy-based sorting. No training required. Factory accepts any `fixed_N` model type (e.g. `fixed_12taxids`, `fixed_5`).
- **Standardized recall metric** in sort-strategy comparison — `max_possible_recall` (recall on all leaves, no cutoff) and `standardized_recall` (actual / max, clipped to [0,1]) added to per-dataset records, summary aggregation, printed pivot tables, and a new faceted boxplot `comparison_standardized_recall.png`.
- **Target-percentile vs actual-recall scatter plot** (`comparison_percentile_vs_recall.png`) — 2×2 subplots faceted by sort strategy, raw scatter points per (model, dataset) with mean ± std errorbar overlay per model type, plus y=x reference line.
- **Dedicated output subdirectory** — comparison results now saved to `{analysis_output}/sort_comparison/` (separate from `models/`).

### Fixed

- **Model name parsing in `_create_modeller`** — replaced naive `split('_')[1]` with `re.match(r'fixed_(\d+)')` to correctly extract the taxid count from names like `fixed_12taxids` (was raising `ValueError` on `int('12taxids')`).
- **`FixedCountModeller.predict_cutoff()` positional counter** — replaced `iterrows()` index label (DataFrame index, not a row position) with `enumerate()` for correct percentile calculation; added row filtering (`is_trash == False` and `best_match_is_best == True`) so only valid best-match taxids count toward the N-taxid target, matching the recall computation in `_compute_actual_recall`.

- **`--composition_model_interface` CLI argument** (`deployment/model_evaluation/evaluate.py:83`) — selects composition model variant: `xgb` (default), `xgb_optimized`, `rf`, `gb`, `lr`
- **`composition_model_interface` config field** (`deployment/model_evaluation/config.py:42`) — wired into `EvaluatorConfig`, `from_args()`, `from_dict()`
- **`BaseCompositionModeller` ABC** (`metagenomics_utils/overlap_manager/om_models.py:2140`) — defines shared API (`fit`, `predict_proba`, `predict`, `save_model`, `load_model`, `evaluate_model`, `eval_and_plot`); provides `_pos_weight()`, `_make_column_transformer()`, `_get_classifier()` helpers; `plot_eval()`, `shap_eval_plot()`, `shap_interaction_plot()` moved from old class
- **5 composition model implementations** (`metagenomics_utils/overlap_manager/om_models.py:2340`):
  - `XGBCompositionModeller` (`xgb`) — 300 trees, max_depth=6, sklearn Pipeline + ColumnTransformer
  - `OptunaXGBCompositionModeller` (`xgb_optimized`) — wraps existing `ClusteringPipeline` with Optuna 50-trial CV
  - `RFCompositionModeller` (`rf`) — 300 trees, max_depth=12, balanced class_weight
  - `GBCompositionModeller` (`gb`) — 300 trees, max_depth=5, lr=0.1, subsample=0.8
  - `LRCompositionModeller` (`lr`) — LogisticRegression C=1.0, balanced, stats-only via `remainder='drop'`
- **`_init_composition_modeller()` factory** in `ModelTrainer` (`deployment/model_evaluation/models.py:283`) — dispatches to the correct variant based on config

- **`POST /predict_composition_stop_traversal` endpoint** (`deployment/ml_api/app.py`) — new ML API endpoint that loads the composition dict-bundle, validates input features against stored `feature_names`, runs `pipeline.predict()` / `pipeline.predict_proba()`, and returns `stop_traversal` boolean + probability
- **`CompositionStopTraversalRequest` / `CompositionStopTraversalResult` schemas** (`deployment/ml_api/validation/schemas.py`) — request accepts `features: Dict[str, float]`, response returns `stop_traversal`, `probability`, `confidence_score`
- **`"composition"` to `PROJECT_MODELS`** and `ModelFile.project_files` (`deployment/ml_api/config.py`) — composition model is now cacheable and warm-loadable alongside recall variants and televir
- **`predict_composition_stop_traversal()` client method** (`deployment/ml_api_client/client.py`) and CLI subcommand `predict-composition-stop` (`deployment/ml_api_client/cli.py`) — wraps the new endpoint
- **`ml_api/README.md` documentation** — updated endpoints table (removed stale entries, added `predict_recall_cutoff_from_table` and `predict_composition_stop_traversal`), added model types table, corrected reload and curl examples
- **`ml_api_client/README.md` documentation** — added composition-stop Python example, JSON input format, and workflow step

### Removed

- **`CompositionModeller` class** (`metagenomics_utils/overlap_manager/om_models.py`) — monolithic class removed; replaced by `BaseCompositionModeller` ABC + 5 implementations

### Changed

- **`ModelTrainer.train_models()`** (`deployment/model_evaluation/models.py`) — now calls `_init_composition_modeller()` + `train_test_split` + `fit()` instead of `CompositionModeller(training_df).train_model()`
- **`BaseCompositionModeller.save_model()` / `load_model()`** (`metagenomics_utils/overlap_manager/om_models.py:2189`) — now writes/reads a dict bundle (`model_type`, `pipeline`, `feature_names`) instead of a flat pipeline dump; `load_model` still accepts legacy flat pipelines for backward compatibility
- **`traversal_with_prediction()` and `predict_data_set_clades_composition()` type hints** (`metagenomics_utils/overlap_manager/om_models.py`) — parameter type changed from `CompositionModeller` to `BaseCompositionModeller`
- **`TrainedModels.composition_modeller` type hint** (`deployment/model_evaluation/models.py:77`) — `Optional[BaseCompositionModeller]`
- **`analyze_samples.py` and `study_deploy.py`** — updated to `XGBCompositionModeller()` with no mock DataFrame; `study_deploy.py` uses `train_test_split` + `fit()` + `eval_and_plot()` pattern
- **`models.py` docstrings** — updated from `CompositionModeller` to `BaseCompositionModeller + 5 variants`
- **`composition_comparison_outputs/README.md`** — `CompositionModeller` reference → `BaseCompositionModeller`
- **README.md** — modules table, composition models section (5-variant table), arguments table (`--composition_model_interface`) updated

### Fixed

- **`predict_data_set_clades_composition()`** — no longer accesses `.pipeline.predict_proba()` directly; calls `.predict_proba()` on the ABC, making it compatible with all variants

### Removed

- **Legacy `GET /predict_recall_cutoff` endpoint** — removed from `app.py`; the only recall endpoint is now `POST /predict_recall_cutoff_from_table`
- **Legacy `POST /predict/{model_type}` generic endpoint** — removed from `app.py`
- **`PredictRequest`/`PredictResponse` schemas** — removed from `app.py` (no longer needed after generic predict deletion)
- **`TAXON_ORDER` and `RecallCutoffRequest` schemas** — removed from `validation/schemas.py` (only used by the deleted legacy endpoint)
- **`"recall"` from `PROJECT_MODELS` and `ModelFile.project_files`** — removed from `config.py` (all recall variants are now addressed by `recall_{variant}` keys via `RECALL_MODEL_VARIANTS`)
- **Old data files** — deleted `recall_results_cache.parquet`, `taxids_to_use_cache.parquet`, and `models/recall_gp_clf_pipeline.pkl` (flat-format bundle, replaced by dict-bundle format)
- **Dead `"recall"` branch in `registry.py::all_model_keys()`** — simplified to iterate `PROJECT_MODELS` directly

### Changed

- **`/reload` and `/reload/{model_type}`** — updated to use `all_model_keys()` instead of `PROJECT_MODELS` for variant-aware cache management

- **`RecallFeatureTransformer`** — new sklearn `TransformerMixin`/`BaseEstimator` in `deployment/model_evaluation/features.py` that learns reference taxa from `taxids_to_use` and collapses raw m_stats matrices into single-row feature+target vectors, replacing the standalone `predict_recall_cutoff_vars` feature extraction

- **`xgb_multioutput()` function** in `metagenomics_utils/overlap_manager/om_models.py` — plain XGBoost multi-output regressor (`MultiOutputRegressor(XGBRegressor(n_estimators=100))`), registered as `'xgb'` in `InjectModellerInterface.model_map` (was previously falling through to RandomForest)

- **`DirectXGBRecallModeller`** class in `metagenomics_utils/overlap_manager/om_models.py` — direct fraction regression with `XGBRegressor(n_estimators=200)` and asymmetric sample weights (3:1 under:over penalty)

- **`direct_xgb` model interface** in `deployment/model_evaluation/models.py` — `ModelTrainer` routes `'direct_xgb'` to `DirectXGBRecallModeller`

- **`composition_model_comparison.py`** — experimental comparison script for the clade-composition (`stop_traversal`) prediction task, analogous to `last_tp_division_prediction_second.py`; loads cached training data, replicates production preprocessing, trains and evaluates multiple classifiers (RandomForest, GradientBoosting, LogisticRegression, XGBoost, LightGBM) side-by-side; outputs metrics table, ROC/PR curves, confusion matrices, feature importance, and training-time comparison to `composition_comparison_outputs/`

- **`composition_comparison_outputs/README.md`** — documents pipeline, input data, algorithms compared, results summary, and output file descriptions for the composition model comparison

- **`last_tp_division_outputs/README.md`** — documents all models (6 Baseline GP kernels, Hurdle, RF, XGBoost, NGBoost, GPCLFThreshold, DirectXGB, DirectRF, GP-DirectXGB), decision approaches (CLF/Likelihood/Direct), pipeline steps, input data, results summary with τ=1.00 metrics and cross-τ tradeoff, and output file descriptions

### Changed

- **`data_set_divide` default 5 → 16** across config and CLI — `EvaluatorConfig.data_set_divide` default in `config.py` and `--data_set_divide` argparse default in `evaluate.py` changed from 5 to 16; `gp_clf` mode no longer overrides `data_set_divide` to 20, now respects the config default

- **`composition_model_comparison.py` — ROC axes swapped** — TPR now on x-axis, FPR on y-axis; axis labels updated accordingly

- **`composition_model_comparison.py` — optuna import crash fixed** — moved `import optuna` to the model-availability check so that XGBoost+Optuna is gracefully skipped when optuna is unimportable (e.g. missing `_sqlite3`), rather than crashing mid-run

- **Recall model training pipeline — feature preprocessing moved into the model class**:
  - `metagenomics_utils/overlap_manager/om_models.py`: `RecallModeller`, `CutoffRecallModeller`, and `GPCLFRecallModeller` now accept either `recall_training_results` (legacy, for backward compat) or omit it for the new transformer-based path; added `fit(m_stats_matrices, taxids_to_use)` method; `save_model`/`load_model` include the transformer when available; `predict_cutoff` and `predict` handle both raw m_stats and processed inputs
  - `cut_off_recall_prediction` simplified — passes raw m_stats directly to `modeller.predict_cutoff()`; removed internal `predict_recall_cutoff_vars` call
  - `predict_recall_cutoff_vars` now emits a `DeprecationWarning` pointing to `RecallFeatureTransformer`
  - `deployment/model_evaluation/models.py`: `run_data_retrieval` collects raw recall matrices alongside processed features; `train_models` calls new `fit()` API when computing fresh data, falls back to `train_model()` for cached data
  - `deployment/model_evaluation/__init__.py`: exports `RecallFeatureTransformer`
- **Analysis scripts updated for transformer-based API**:
  - `deployment/analysis/analyze_samples.py`: removed local `predict_recall_cutoff_vars`; `_load_recall_modeller` constructs with `recall_training_results=None`; inference passes raw m_stats to `predict()`/`predict_cutoff()` with fallback for legacy models via `_extract_legacy_features()`
  - `deployment/ml_api/train_recall.py`: rewritten to accept CLI args and collect raw m_stats matrices via `OverlapManager`; uses `GPCLFRecallModeller.fit()` instead of `train_model()`

- **Terminology: `trash` → `spurious`** across output-facing code and documentation:
  - `metagenomics_utils/overlap_manager/om_models.py`: renamed `get_trash_composition()` → `get_spurious_composition()`
  - `deployment/model_evaluation/metrics.py`: renamed `compute_trash_hit_stats()` → `compute_spurious_hit_stats()`; output columns `trash_hit_count`/`trash_hit_reads_mapped` → `spurious_hit_count`/`spurious_hit_reads_mapped`
  - `deployment/model_evaluation/evaluate.py`: renamed `analyze_trash_hit_distribution()` → `analyze_spurious_hit_distribution()`; output filenames and plot labels updated
  - `deployment/model_evaluation/visualization.py`: renamed `plot_trash_heatmap()` → `plot_spurious_heatmap()`; plot title and filename updated
  - `deployment/model_evaluation/result_models.py`: `DatasetResult.trash_composition` field → `spurious_composition`; `CrossHitMetrics.trash_hit_counts_per_class` → `spurious_hit_counts_per_class`; JSON keys and TSV filenames updated
  - `deployment/model_evaluation/batch_evaluator.py`, `dataset_processor.py`, `models.py`: all variable/field references updated
  - `deployment/analysis/study_deploy.py`: updated import and call site
  - Internal `is_trash` DataFrame column name retained for backward compatibility

- **Output format: added agent-parseable JSON**:
  - Added `BatchEvaluationResult.write_agent_output(filepath)` in `result_models.py`
  - Produces `evaluation_results_agent.json` with aggregate statistics (mean/median/std/quartiles for precision, recall, cross-hit metrics, probability metrics) alongside per-dataset records
  - Integrated into `BatchEvaluator.save_results()`

- **Consistent column naming in TSV output**:
  | Old name | New name |
  |---|---|
  | `overall_precision_raw` | `precision_best_match` |
  | `fuzzy_precision_raw` | `precision_fuzzy` |
  | `fuzzy_precision_cov_filtered` | `precision_fuzzy_cov_filtered` |
  | `clade_precision_full` | `precision_clade_full` |
  | `clade_precision_post` | `precision_clade_post_cleanup` |
  | `clade_precision_fixed` | `precision_clade_fixed` |
  | `recall_raw` | `recall_baseline` |
  | `recall_cov_filtered` | `recall_baseline_cov_filtered` |
  | `clade_recall` | `recall_clade_pre_cleanup` |
  | `clade_recall_clean` | `recall_clade_post_cleanup` |
  | `recall_filtered_leaves` | `recall_after_recall_filter` |
  | `recall_fixed_filter` | `recall_fixed_max_12` |
  - All downstream plot code in `visualization.py`, `batch_evaluator.save_summary_statistics()`, and `result_models.get_summary_stats()` updated accordingly

### Added

- **`RecallModeller.evaluate_prediction(X_raw, target_recall)`** (`metagenomics_utils/overlap_manager/om_models.py:685`) — new method that compares model predictions against ground truth targets after `RecallFeatureTransformer.transform()`; returns `last_best_match_relindex` for all types; per-division modellers additionally get `per_division_recall_errors` (dict) and `per_division_recall_rmse`; cutoff modellers get `predicted_k_min`, `actual_k_min`, and `cutoff_error`
- **`cut_off_recall_prediction()` evaluation wiring** (`metagenomics_utils/overlap_manager/om_models.py:1712`) — now calls `modeller.evaluate_prediction()` and merges results into the returned metrics dict
- **`RecallMetrics.recall_metrics` field** (`deployment/model_evaluation/result_models.py:69`) — new `dict` field storing per-dataset recall evaluation metrics alongside existing scalar recall values
- **metrics storage in `dataset_processor.py:410`** — stores the full `metrics_dict` as `result.recall.recall_metrics` after the recall filter step
- **summary_results flattening** (`deployment/model_evaluation/batch_evaluator.py:250-256`) — scalar metrics flattened into `recall_metric_*` columns; `per_division_recall_errors` dict flattened into `recall_metric_error_index_recall_*` columns for direct plotting
- **5 recall modeller evaluation plots** (`deployment/model_evaluation/visualization.py:469-541`):
  - `plot_recall_error_boxplot` — boxplot of pred − true error per division index (with zero line)
  - `plot_recall_rmse_distribution` — histogram + KDE of per-division RMSE
  - `plot_last_best_match_vs_rmse` — scatter of relindex vs RMSE colored by `recall_after_recall_filter`
  - `plot_cutoff_error_histogram` — histogram of `predicted_k_min − actual_k_min`
  - `plot_cutoff_confusion_matrix` — heatmap of predicted vs actual k_min
  - All wired into `ResultVisualizer.plot_all()`, self-guarding on column existence
