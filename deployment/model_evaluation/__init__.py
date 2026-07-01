"""
Evaluator module for deployment benchmark.

This module provides classes and functions for evaluating taxonomic
classification models on metagenomics datasets.
"""

from .base import (
    BaseEvaluator,
    DataLoaderBase,
    DatasetProcessorBase,
    EvaluationMetadata,
    ModelTrainerBase,
    ProgressTracker,
    TqdmProgressTracker,
    VisualizationBase,
)
from .batch_evaluator import BatchEvaluator
from .config import (
    EvaluatorConfig,
    LoggingConfig,
    MLflowConfig,
    ModelConfig,
    TrainedModels,
    VisualizationConfig,
)
from .data_loader import (
    DataLoader,
    establish_taxids_to_use,
    expand_input_data,
    get_dataset_folders,
    load_taxid_plan,
    output_parse,
    process_clade_report,
    retrieve_simulation_input,
    split_train_test_folders,
)
from .dataset_processor import DatasetProcessor
from .exceptions import (
    ConfigurationError,
    DataLoadError,
    EvaluatorError,
    ModelError,
    OverlapManagerError,
    PredictionError,
    ResultsAggregationError,
)
from .features import (
    RecallFeatureTransformer,
)
from .logging_config import (
    JSONFormatter,
    LogContext,
    StructuredLogger,
    get_logger,
    setup_logging,
)
from .metrics import (
    compute_clade_recall,
    compute_cross_hit_stats,
    compute_mstats_precision,
    compute_purity,
    compute_recall,
    safe_divide,
)
from .models import (
    CrossValidator,
    MLflowTracker,
    ModelEvaluator,
    ModelRegistry,
    ModelTrainer,
    TrainingCrossHitAnalyzer,
    train_and_evaluate,
)
from .result_models import (
    BatchEvaluationResult,
    CrossHitMetrics,
    DatasetResult,
    PrecisionMetrics,
    RecallMetrics,
    create_empty_result,
)
from .visualization import (
    PlotlyVisualizer,
    ReportGenerator,
    ResultVisualizer,
    generate_report,
)

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
