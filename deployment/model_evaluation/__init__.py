"""
Evaluator module for deployment benchmark.

This module provides classes and functions for evaluating taxonomic
classification models on metagenomics datasets.

All sub-module imports are lazy to avoid pulling in heavy dependencies
(e.g. BioPython) when only a subset of the module is needed — such as
when unpickling model files that reference ``RecallFeatureTransformer``.
"""


_LazyImports = {
    "BaseEvaluator": ".base",
    "DataLoaderBase": ".base",
    "DatasetProcessorBase": ".base",
    "EvaluationMetadata": ".base",
    "ModelTrainerBase": ".base",
    "ProgressTracker": ".base",
    "TqdmProgressTracker": ".base",
    "VisualizationBase": ".base",
    "BatchEvaluator": ".batch_evaluator",
    "EvaluatorConfig": ".config",
    "LoggingConfig": ".config",
    "MLflowConfig": ".config",
    "ModelConfig": ".config",
    "TrainedModels": ".config",
    "VisualizationConfig": ".config",
    "DataLoader": ".data_loader",
    "establish_taxids_to_use": ".data_loader",
    "expand_input_data": ".data_loader",
    "get_dataset_folders": ".data_loader",
    "load_taxid_plan": ".data_loader",
    "output_parse": ".data_loader",
    "process_clade_report": ".data_loader",
    "retrieve_simulation_input": ".data_loader",
    "split_train_test_folders": ".data_loader",
    "DatasetProcessor": ".dataset_processor",
    "ConfigurationError": ".exceptions",
    "DataLoadError": ".exceptions",
    "EvaluatorError": ".exceptions",
    "ModelError": ".exceptions",
    "OverlapManagerError": ".exceptions",
    "PredictionError": ".exceptions",
    "ResultsAggregationError": ".exceptions",
    "RecallFeatureTransformer": ".features",
    "JSONFormatter": ".logging_config",
    "LogContext": ".logging_config",
    "StructuredLogger": ".logging_config",
    "get_logger": ".logging_config",
    "setup_logging": ".logging_config",
    "compute_clade_recall": ".metrics",
    "compute_cross_hit_stats": ".metrics",
    "compute_mstats_precision": ".metrics",
    "compute_purity": ".metrics",
    "compute_recall": ".metrics",
    "safe_divide": ".metrics",
    "CrossValidator": ".models",
    "MLflowTracker": ".models",
    "ModelEvaluator": ".models",
    "ModelRegistry": ".models",
    "ModelTrainer": ".models",
    "TrainingCrossHitAnalyzer": ".models",
    "train_and_evaluate": ".models",
    "BatchEvaluationResult": ".result_models",
    "CrossHitMetrics": ".result_models",
    "DatasetResult": ".result_models",
    "PrecisionMetrics": ".result_models",
    "RecallMetrics": ".result_models",
    "create_empty_result": ".result_models",
    "PlotlyVisualizer": ".visualization",
    "ReportGenerator": ".visualization",
    "ResultVisualizer": ".visualization",
    "generate_report": ".visualization",
}


def __getattr__(name):
    import importlib

    mod_path = _LazyImports.get(name)
    if mod_path is not None:
        mod = importlib.import_module(mod_path, __name__)
        return getattr(mod, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    # Config
    "EvaluatorConfig",
    "TrainedModels",
    "ModelConfig",
    "VisualizationConfig",
    "MLflowConfig",
    "LoggingConfig",
    # Exceptions
    "EvaluatorError",
    "DataLoadError",
    "PredictionError",
    "OverlapManagerError",
    "ConfigurationError",
    "ModelError",
    "ResultsAggregationError",
    # Result Models
    "PrecisionMetrics",
    "RecallMetrics",
    "CrossHitMetrics",
    "DatasetResult",
    "BatchEvaluationResult",
    "create_empty_result",
    # Metrics
    "compute_purity",
    "compute_mstats_precision",
    "compute_recall",
    "compute_clade_recall",
    "compute_cross_hit_stats",
    "safe_divide",
    # Data Loader
    "retrieve_simulation_input",
    "process_clade_report",
    "output_parse",
    "establish_taxids_to_use",
    "expand_input_data",
    "load_taxid_plan",
    "get_dataset_folders",
    "split_train_test_folders",
    "DataLoader",
    # Processing
    "DatasetProcessor",
    "BatchEvaluator",
    # Visualization
    "ResultVisualizer",
    "PlotlyVisualizer",
    "ReportGenerator",
    "generate_report",
    # Features
    "RecallFeatureTransformer",
    # Models
    "ModelTrainer",
    "train_and_evaluate",
    "ModelRegistry",
    "MLflowTracker",
    "CrossValidator",
    "ModelEvaluator",
    "TrainingCrossHitAnalyzer",
    # Logging
    "setup_logging",
    "get_logger",
    "StructuredLogger",
    "JSONFormatter",
    "LogContext",
    # Base Classes
    "BaseEvaluator",
    "DatasetProcessorBase",
    "ModelTrainerBase",
    "DataLoaderBase",
    "VisualizationBase",
    "ProgressTracker",
    "TqdmProgressTracker",
    "EvaluationMetadata",
]
