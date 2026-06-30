"""
Evaluator module for deployment benchmark.

This module provides classes and functions for evaluating taxonomic
classification models on metagenomics datasets.
"""

from .config import (
    EvaluatorConfig,
    TrainedModels,
    ModelConfig,
    VisualizationConfig,
    MLflowConfig,
    LoggingConfig,
)
from .exceptions import (
    EvaluatorError,
    DataLoadError,
    PredictionError,
    OverlapManagerError,
    ConfigurationError,
    ModelError,
    ResultsAggregationError,
)
from .result_models import (
    PrecisionMetrics,
    RecallMetrics,
    CrossHitMetrics,
    DatasetResult,
    BatchEvaluationResult,
    create_empty_result,
)
from .metrics import (
    compute_purity,
    compute_mstats_precision,
    compute_recall,
    compute_clade_recall,
    compute_cross_hit_stats,
    safe_divide,
)
from .data_loader import (
    retrieve_simulation_input,
    process_clade_report,
    output_parse,
    establish_taxids_to_use,
    expand_input_data,
    load_taxid_plan,
    get_dataset_folders,
    split_train_test_folders,
    DataLoader,
)
from .dataset_processor import DatasetProcessor
from .batch_evaluator import BatchEvaluator
from .visualization import (
    ResultVisualizer,
    PlotlyVisualizer,
    ReportGenerator,
    generate_report,
)
from .features import (
    RecallFeatureTransformer,
)
from .models import (
    ModelTrainer,
    train_and_evaluate,
    ModelRegistry,
    MLflowTracker,
    CrossValidator,
    ModelEvaluator,
    TrainingCrossHitAnalyzer,
)
from .logging_config import (
    setup_logging,
    get_logger,
    StructuredLogger,
    JSONFormatter,
    LogContext,
)
from .base import (
    BaseEvaluator,
    DatasetProcessorBase,
    ModelTrainerBase,
    DataLoaderBase,
    VisualizationBase,
    ProgressTracker,
    TqdmProgressTracker,
    EvaluationMetadata,
)


__all__ = [
    # Config
    'EvaluatorConfig',
    'TrainedModels',
    'ModelConfig',
    'VisualizationConfig',
    'MLflowConfig',
    'LoggingConfig',
    # Exceptions
    'EvaluatorError',
    'DataLoadError',
    'PredictionError',
    'OverlapManagerError',
    'ConfigurationError',
    'ModelError',
    'ResultsAggregationError',
    # Result Models
    'PrecisionMetrics',
    'RecallMetrics',
    'CrossHitMetrics',
    'DatasetResult',
    'BatchEvaluationResult',
    'create_empty_result',
    # Metrics
    'compute_purity',
    'compute_mstats_precision',
    'compute_recall',
    'compute_clade_recall',
    'compute_cross_hit_stats',
    'safe_divide',
    # Data Loader
    'retrieve_simulation_input',
    'process_clade_report',
    'output_parse',
    'establish_taxids_to_use',
    'expand_input_data',
    'load_taxid_plan',
    'get_dataset_folders',
    'split_train_test_folders',
    'DataLoader',
    # Processing
    'DatasetProcessor',
    'BatchEvaluator',
    # Visualization
    'ResultVisualizer',
    'PlotlyVisualizer',
    'ReportGenerator',
    'generate_report',
    # Features
    'RecallFeatureTransformer',
    # Models
    'ModelTrainer',
    'train_and_evaluate',
    'ModelRegistry',
    'MLflowTracker',
    'CrossValidator',
    'ModelEvaluator',
    'TrainingCrossHitAnalyzer',
    # Logging
    'setup_logging',
    'get_logger',
    'StructuredLogger',
    'JSONFormatter',
    'LogContext',
    # Base Classes
    'BaseEvaluator',
    'DatasetProcessorBase',
    'ModelTrainerBase',
    'DataLoaderBase',
    'VisualizationBase',
    'ProgressTracker',
    'TqdmProgressTracker',
    'EvaluationMetadata',
]
