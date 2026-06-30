"""
Custom exceptions for the evaluator module.
"""


class EvaluatorError(Exception):
    """Base exception for evaluator module."""
    
    def __init__(self, message: str, details: dict = None):
        self.message = message
        self.details = details or {}
        super().__init__(self.message)
    
    def __str__(self):
        if self.details:
            details_str = ", ".join(f"{k}={v}" for k, v in self.details.items())
            return f"{self.message} ({details_str})"
        return self.message


class DataLoadError(EvaluatorError):
    """Failed to load dataset files."""
    
    def __init__(self, dataset: str, reason: str, filepath: str = None):
        details = {"dataset": dataset, "reason": reason}
        if filepath:
            details["filepath"] = filepath
        super().__init__(f"Failed to load dataset '{dataset}': {reason}", details)
        self.dataset = dataset
        self.reason = reason
        self.filepath = filepath


class PredictionError(EvaluatorError):
    """Failed during model prediction."""
    
    def __init__(self, dataset: str, phase: str, reason: str):
        details = {"dataset": dataset, "phase": phase, "reason": reason}
        super().__init__(f"Prediction failed for '{dataset}' at {phase}: {reason}", details)
        self.dataset = dataset
        self.phase = phase
        self.reason = reason


class OverlapManagerError(EvaluatorError):
    """OverlapManager operation failed."""
    
    def __init__(self, dataset: str, operation: str, reason: str):
        details = {"dataset": dataset, "operation": operation, "reason": reason}
        super().__init__(f"OverlapManager error for '{dataset}' during {operation}: {reason}", details)
        self.dataset = dataset
        self.operation = operation
        self.reason = reason


class ConfigurationError(EvaluatorError):
    """Invalid configuration parameters."""
    
    def __init__(self, parameter: str, value, reason: str):
        details = {"parameter": parameter, "value": str(value), "reason": reason}
        super().__init__(f"Invalid configuration for '{parameter}': {reason}", details)
        self.parameter = parameter
        self.value = value
        self.reason = reason


class ModelError(EvaluatorError):
    """Model-related error."""
    
    def __init__(self, model_type: str, operation: str, reason: str):
        details = {"model_type": model_type, "operation": operation, "reason": reason}
        super().__init__(f"Model error ({model_type}) during {operation}: {reason}", details)
        self.model_type = model_type
        self.operation = operation
        self.reason = reason


class ResultsAggregationError(EvaluatorError):
    """Failed to aggregate results."""
    
    def __init__(self, reason: str, dataset_count: int = None):
        details = {"reason": reason}
        if dataset_count is not None:
            details["dataset_count"] = dataset_count
        super().__init__(f"Failed to aggregate results: {reason}", details)
        self.dataset_count = dataset_count
