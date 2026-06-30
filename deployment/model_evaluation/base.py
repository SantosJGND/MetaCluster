"""
Abstract base classes for evaluator module.

Provides interfaces for extensible evaluator, processor, and trainer implementations.
"""

from abc import ABC, abstractmethod
from typing import Any, Generic, TypeVar, Optional, Protocol
from dataclasses import dataclass

import pandas as pd


T = TypeVar('T')
TResult = TypeVar('TResult')


class BaseEvaluator(ABC):
    """
    Abstract base class for evaluators.
    
    Subclasses must implement the evaluate method to process
    datasets and return results.
    """
    
    @abstractmethod
    def evaluate(self, datasets: list[str]) -> 'BatchEvaluationResult':
        """
        Evaluate all datasets.
        
        Args:
            datasets: List of dataset identifiers to evaluate
            
        Returns:
            BatchEvaluationResult containing aggregated results
        """
        pass
    
    @abstractmethod
    def save_results(self, results: 'BatchEvaluationResult', output_dir: str) -> None:
        """
        Save evaluation results to disk.
        
        Args:
            results: Evaluation results to save
            output_dir: Output directory path
        """
        pass


class DatasetProcessorBase(ABC, Generic[TResult]):
    """
    Abstract base class for dataset processors.
    
    Defines the interface for processing individual datasets.
    
    Type Parameters:
        TResult: Type of result returned by processing
    """
    
    @abstractmethod
    def process(self, dataset_name: str) -> Optional[TResult]:
        """
        Process a single dataset.
        
        Args:
            dataset_name: Identifier for the dataset to process
            
        Returns:
            Processing result, or None if dataset should be skipped
        """
        pass
    
    @abstractmethod
    def validate_dataset(self, dataset_name: str) -> bool:
        """
        Validate that a dataset can be processed.
        
        Args:
            dataset_name: Identifier for the dataset
            
        Returns:
            True if dataset is valid for processing
        """
        pass


class ModelTrainerBase(ABC):
    """
    Abstract base class for model trainers.
    
    Defines the interface for training ML models.
    """
    
    @abstractmethod
    def train(self, datasets: list[str]) -> 'TrainedModels':
        """
        Train models on the provided datasets.
        
        Args:
            datasets: List of dataset identifiers for training
            
        Returns:
            TrainedModels containing all trained model instances
        """
        pass
    
    @abstractmethod
    def save_models(self, output_path: str) -> None:
        """
        Save trained models to disk.
        
        Args:
            output_path: Directory to save models
        """
        pass
    
    @abstractmethod
    def load_models(self, input_path: str) -> 'TrainedModels':
        """
        Load pre-trained models from disk.
        
        Args:
            input_path: Directory containing saved models
            
        Returns:
            TrainedModels loaded from disk
        """
        pass


class DataLoaderBase(ABC):
    """
    Abstract base class for data loaders.
    
    Defines the interface for loading and preparing data.
    """
    
    @abstractmethod
    def load(self) -> 'DataLoader':
        """
        Load all data.
        
        Returns:
            Self for method chaining
        """
        pass
    
    @abstractmethod
    def get_training_data(self) -> list[str]:
        """
        Get list of training dataset identifiers.
        
        Returns:
            List of training dataset names
        """
        pass
    
    @abstractmethod
    def get_test_data(self) -> list[str]:
        """
        Get list of test dataset identifiers.
        
        Returns:
            List of test dataset names
        """
        pass


class VisualizationBase(ABC):
    """
    Abstract base class for result visualization.
    
    Defines the interface for generating visualizations.
    """
    
    @abstractmethod
    def plot_all(self, results: 'BatchEvaluationResult') -> None:
        """
        Generate all visualizations.
        
        Args:
            results: Evaluation results to visualize
        """
        pass
    
    @abstractmethod
    def save(self, output_dir: str) -> None:
        """
        Save all visualizations to disk.
        
        Args:
            output_dir: Directory to save visualizations
        """
        pass


class ResultAggregator(ABC, Generic[T]):
    """
    Abstract base class for result aggregation.
    
    Defines how individual results are combined into
    batch results.
    """
    
    @abstractmethod
    def aggregate(self, results: list[T]) -> 'BatchEvaluationResult':
        """
        Aggregate individual results into batch result.
        
        Args:
            results: List of individual results
            
        Returns:
            Aggregated batch result
        """
        pass


class Configurable(Protocol):
    """
    Protocol for configuration objects.
    
    Objects implementing this protocol can be configured
    via YAML, JSON, or environment variables.
    """
    
    @classmethod
    def from_dict(cls, config: dict) -> Any:
        """Create instance from dictionary."""
        ...
    
    @classmethod
    def from_file(cls, path: str) -> Any:
        """Load instance from file."""
        ...
    
    def to_dict(self) -> dict:
        """Convert configuration to dictionary."""
        ...


@dataclass
class EvaluationMetadata:
    """Metadata about an evaluation run."""
    start_time: Optional[str] = None
    end_time: Optional[str] = None
    dataset_count: int = 0
    successful_count: int = 0
    failed_count: int = 0
    config: Optional[dict] = None
    errors: list[str] = None
    
    def __post_init__(self):
        if self.errors is None:
            self.errors = []


class ProgressTracker(ABC):
    """
    Abstract base class for progress tracking.
    
    Implementations can track progress via tqdm, rich,
    or other progress display libraries.
    """
    
    @abstractmethod
    def start(self, total: int, description: str = "") -> None:
        """Start progress tracking."""
        pass
    
    @abstractmethod
    def update(self, n: int = 1) -> None:
        """Update progress by n steps."""
        pass
    
    @abstractmethod
    def finish(self) -> None:
        """Finish and clean up progress tracking."""
        pass
    
    @abstractmethod
    def set_description(self, description: str) -> None:
        """Update the progress description."""
        pass


class TqdmProgressTracker(ProgressTracker):
    """Progress tracker implementation using tqdm."""
    
    def __init__(self):
        self._progress = None
    
    def start(self, total: int, description: str = "") -> None:
        from tqdm import tqdm
        self._progress = tqdm(total=total, desc=description, unit="item")
    
    def update(self, n: int = 1) -> None:
        if self._progress:
            self._progress.update(n)
    
    def finish(self) -> None:
        if self._progress:
            self._progress.close()
    
    def set_description(self, description: str) -> None:
        if self._progress:
            self._progress.set_description(description)
