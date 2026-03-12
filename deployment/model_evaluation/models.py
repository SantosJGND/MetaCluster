"""
Model training module for evaluator.

Handles training of RecallModeller, CompositionModeller, and CrossHitModeller.
"""

import os
import logging
import tempfile
from pathlib import Path
from typing import Tuple, Any, Optional

import pandas as pd
import joblib

from .config import TrainedModels, ModelConfig

from metagenomics_utils.overlap_manager import OverlapManager
from metagenomics_utils.overlap_manager.om_models import (
    RecallModeller,
    CompositionModeller, 
    CrossHitModeller,
    data_set_traversal_with_precision,
    cross_hit_prediction_matrix,
    predict_recall_cutoff_vars,
)

from .config import EvaluatorConfig
from .result_models import BatchEvaluationResult
from .exceptions import ModelError


logger = logging.getLogger(__name__)


class ModelTrainer:
    """
    Handles model training for the evaluator.
    
    Trains three models:
    - RecallModeller: Predicts recall cutoff
    - CompositionModeller: Predicts clade composition  
    - CrossHitModeller: Predicts cross-hit probability
    """
    
    def __init__(
        self,
        config: EvaluatorConfig,
        ncbi_wrapper: Any,
        input_tax_df: pd.DataFrame,
        taxids_to_use: pd.DataFrame,
    ):
        """
        Initialize model trainer.
        
        Args:
            config: EvaluatorConfig instance
            ncbi_wrapper: NCBI TaxonomistWrapper instance
            input_tax_df: Input tax DataFrame
            taxids_to_use: Taxids to use for evaluation
        """
        self.config = config
        self.ncbi = ncbi_wrapper
        self.input_tax_df = input_tax_df
        self.taxids_to_use = taxids_to_use
        
        self.recall_modeller: Optional[RecallModeller] = None
        self.composition_modeller: Optional[CompositionModeller] = None
        self.crosshit_modeller: Optional[CrossHitModeller] = None
        
        self._training_results: Optional[pd.DataFrame] = None
        self._prediction_results: Optional[pd.DataFrame] = None
        self._recall_results: Optional[pd.DataFrame] = None

        self.models: Optional[TrainedModels] = None
    
    def run_data_retrieval(self, training_folders: list) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Run data retrieval for training datasets.
        
        Args:
            training_folders: List of training dataset folder names
            
        Returns:
            Tuple of (training_results_df, prediction_results_df, recall_results_df)
        """
        training_results = []
        prediction_results = []
        recall_results = []
        
        for data_set_name in training_folders:
            try:
                overlap_manager = OverlapManager(
                    os.path.join(self.config.study_output_filepath, data_set_name, "clustering")
                )
                
                if overlap_manager.m_stats_matrix.empty:
                    logger.warning(f"Skipping {data_set_name}: no mapped reads")
                    continue
                
                result_df = data_set_traversal_with_precision(
                    data_set_name, 
                    self.config.study_output_filepath, 
                    self.ncbi, 
                    overlap_manager, 
                    self.taxids_to_use, 
                    tax_level=self.config.tax_level
                )
                
                prediction_matrix = cross_hit_prediction_matrix(
                    data_set_name, 
                    self.config.study_output_filepath, 
                    self.ncbi, 
                    overlap_manager, 
                    self.taxids_to_use, 
                    tax_level=self.config.tax_level
                )
                
                m_stats_stats_matrix = self._get_m_stats_matrix(data_set_name, overlap_manager)
                
                if m_stats_stats_matrix.empty:
                    continue
                
                recall_stats = predict_recall_cutoff_vars(
                    self.config.data_set_divide, 
                    data_set_name, 
                    m_stats_stats_matrix, 
                    self.taxids_to_use, 
                    tax_level=self.config.tax_level
                )
                
                if result_df is not None and not result_df.empty:
                    training_results.append(result_df)
                if prediction_matrix is not None and not prediction_matrix.empty:
                    prediction_results.append(prediction_matrix)
                if recall_stats is not None and not recall_stats.empty:
                    recall_results.append(recall_stats)
                    
            except Exception as e:
                logger.error(f"Error processing {data_set_name}: {e}")
                continue
        
        training_results_df = pd.concat(training_results, ignore_index=True) if training_results else pd.DataFrame()
        prediction_results_df = pd.concat(prediction_results, ignore_index=True) if prediction_results else pd.DataFrame()
        recall_results_df = pd.concat(recall_results, ignore_index=True) if recall_results else pd.DataFrame()
        
        self._training_results = training_results_df
        self._prediction_results = prediction_results_df
        self._recall_results = recall_results_df
        
        logger.info(f"Retrieved training data: {len(training_results_df)} results, {len(prediction_results_df)} predictions")
        
        return training_results_df, prediction_results_df, recall_results_df
    
    def _get_m_stats_matrix(self, data_set_name: str, overlap_manager: OverlapManager) -> pd.DataFrame:
        """
        Get m_stats matrix for a dataset.
        
        Args:
            data_set_name: Name of dataset
            overlap_manager: OverlapManager instance
            
        Returns:
            m_stats matrix DataFrame
        """
        from metagenomics_utils.overlap_manager.node_stats import get_m_stats_matrix
        return get_m_stats_matrix(
            data_set_name,
            self.config.study_output_filepath,
            self.ncbi,
            overlap_manager
        )
    
    def train_models(self, training_folders: list) -> Tuple[Any, Any, Any]:
        """
        Train all models.
        
        Args:
            training_folders: List of training dataset folder names
            
        Returns:
            Tuple of (recall_model, composition_model, crosshit_model)
        """
        logger.info(f"Starting model training with {len(training_folders)} datasets")
        
        training_df, prediction_df, recall_df = self.run_data_retrieval(training_folders)
        
        if training_df.empty:
            raise ModelError("composition", "train", "No training data available")
        
        self.recall_modeller = RecallModeller(
            recall_trainning_results=recall_df,
            data_set_divide=self.config.data_set_divide
        )
        
        self.composition_modeller = CompositionModeller(
            trainning_results_df=training_df
        )
        
        self.crosshit_modeller = CrossHitModeller(
            prediction_trainning_results_df=prediction_df
        )
        
        logger.info("Training recall model...")
        model_recall, X_test_recall, Y_test_recall = self.recall_modeller.train_model()
        
        logger.info("Training composition model...")
        model_composition, X_train_composition, X_test_composition, y_test_composition = self.composition_modeller.train_model()
        
        logger.info("Training crosshit model...")
        self.crosshit_modeller.train_model()

        self.models = TrainedModels(
            recall_modeller=self.recall_modeller,
            composition_modeller=self.composition_modeller,
            crosshit_modeller=self.crosshit_modeller
        )

        return model_recall, model_composition, self.crosshit_modeller
    
    def save_models(self, output_dir: str) -> None:
        """
        Save trained models to disk.
        
        Args:
            output_dir: Directory to save models
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        if self.recall_modeller is not None:
            self.recall_modeller.save_model(output_dir)
            logger.info(f"Saved recall model to {output_dir}")
        
        if self.composition_modeller is not None:
            self.composition_modeller.save_model(output_dir)
            logger.info(f"Saved composition model to {output_dir}")
        
        if self.crosshit_modeller is not None:
            self.crosshit_modeller.save_model(output_dir)
            logger.info(f"Saved crosshit model to {output_dir}")
    
    def evaluate_models(self, output_dir: str) -> None:
        """
        Evaluate models and save evaluation results.
        
        Args:
            output_dir: Directory to save evaluation results
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        if self.recall_modeller is not None and hasattr(self.recall_modeller, 'model'):
            logger.info("Evaluating recall model...")
            self.recall_modeller.model_summary(
                self.recall_modeller.model,
                self.recall_modeller.X_test,
                self.recall_modeller.y_test,
                output_dir
            )
        
        if self.composition_modeller is not None:
            logger.info("Evaluating composition model...")
            self.composition_modeller.eval_and_plot(
                self.composition_modeller.X_test,
                self.composition_modeller.y_test,
                output_dir,
                X_train=self.composition_modeller.X_train
            )
    
    def get_trained_models(self) -> Tuple[RecallModeller, CompositionModeller, CrossHitModeller]:
        """
        Get trained model instances.
        
        Returns:
            Tuple of (recall_modeller, composition_modeller, crosshit_modeller)
        """
        if self.recall_modeller is None or self.composition_modeller is None or self.crosshit_modeller is None:
            raise ModelError("all", "get", "Models not trained yet")
        
        return self.recall_modeller, self.composition_modeller, self.crosshit_modeller


def train_and_evaluate(
    config: EvaluatorConfig,
    training_folders: list,
    ncbi_wrapper: Any,
    input_tax_df: pd.DataFrame,
    taxids_to_use: pd.DataFrame,
    models_output_dir: str,
) -> Tuple[RecallModeller, CompositionModeller, CrossHitModeller]:
    """
    Convenience function to train and evaluate all models.
    
    Args:
        config: EvaluatorConfig instance
        training_folders: List of training dataset folder names
        ncbi_wrapper: NCBI TaxonomistWrapper instance
        input_tax_df: Input tax DataFrame
        taxids_to_use: Taxids to use for evaluation
        models_output_dir: Directory to save models and evaluations
        
    Returns:
        Tuple of trained (recall_modeller, composition_modeller, crosshit_modeller)
    """
    trainer = ModelTrainer(
        config=config,
        ncbi_wrapper=ncbi_wrapper,
        input_tax_df=input_tax_df,
        taxids_to_use=taxids_to_use,
    )
    
    trainer.train_models(training_folders)
    trainer.save_models(models_output_dir)
    trainer.evaluate_models(models_output_dir)
    
    return trainer.get_trained_models()


try:
    import mlflow
    from mlflow.tracking import MlflowClient
    MLFLOW_AVAILABLE = True
except ImportError:
    MLFLOW_AVAILABLE = False


class ModelRegistry:
    """
    Registry for loading and saving trained models.
    
    Provides organized model persistence with metadata tracking.
    """
    
    def __init__(self, base_path: str = "models"):
        self.base_path = Path(base_path)
        self.models_dir = self.base_path / "trained"
        self.metadata_file = self.base_path / "metadata.json"
    
    def save_models(
        self,
        models: 'TrainedModels',
        name: str,
        metadata: Optional[dict] = None,
    ) -> Path:
        """
        Save trained models to registry.
        
        Args:
            models: TrainedModels instance
            name: Name for the model version
            metadata: Optional metadata dict
        
        Returns:
            Path to saved models
        """
        import json
        import joblib
        
        model_path = self.models_dir / name
        model_path.mkdir(parents=True, exist_ok=True)
        
        if models.recall_modeller is not None:
            joblib.dump(models.recall_modeller, model_path / "recall_model.joblib")
        
        if models.composition_modeller is not None:
            joblib.dump(models.composition_modeller, model_path / "composition_model.joblib")
        
        if models.crosshit_modeller is not None:
            joblib.dump(models.crosshit_modeller, model_path / "crosshit_model.joblib")
        
        metadata = metadata or {}
        metadata['name'] = name
        metadata['path'] = str(model_path)
        
        existing = []
        if self.metadata_file.exists():
            with open(self.metadata_file) as f:
                existing = json.load(f)
        
        existing.append(metadata)
        
        with open(self.metadata_file, 'w') as f:
            json.dump(existing, f, indent=2)
        
        return model_path
    
    def load_models(self, name: str) -> 'TrainedModels':
        """
        Load models from registry by name.
        
        Args:
            name: Model name to load
        
        Returns:
            TrainedModels instance
        """
        import joblib
        
        model_path = self.models_dir / name
        if not model_path.exists():
            raise FileNotFoundError(f"Model not found: {name}")
        
        models = TrainedModels()
        
        recall_path = model_path / "recall_model.joblib"
        if recall_path.exists():
            models.recall_modeller = joblib.load(recall_path)
        
        comp_path = model_path / "composition_model.joblib"
        if comp_path.exists():
            models.composition_modeller = joblib.load(comp_path)
        
        cross_path = model_path / "crosshit_model.joblib"
        if cross_path.exists():
            models.crosshit_modeller = joblib.load(cross_path)
        
        return models
    
    def list_models(self) -> list[dict]:
        """List all registered models."""
        import json
        
        if not self.metadata_file.exists():
            return []
        
        with open(self.metadata_file) as f:
            return json.load(f)
    
    def delete_model(self, name: str) -> None:
        """Delete a model from the registry."""
        import shutil
        
        model_path = self.models_dir / name
        if model_path.exists():
            shutil.rmtree(model_path)
        
        import json
        if self.metadata_file.exists():
            with open(self.metadata_file) as f:
                models = json.load(f)
            models = [m for m in models if m.get('name') != name]
            with open(self.metadata_file, 'w') as f:
                json.dump(models, f, indent=2)


class MLflowTracker:
    """
    MLflow integration for experiment tracking.
    
    Provides automatic logging of parameters, metrics, and models
    to MLflow.
    """
    
    def __init__(
        self,
        experiment_name: str = "metagenomics-evaluation",
        tracking_uri: Optional[str] = None,
    ):
        if not MLFLOW_AVAILABLE:
            raise ImportError("MLflow is not installed. Install with: pip install mlflow")
        
        self.experiment_name = experiment_name
        self.tracking_uri = tracking_uri
        self._run = None
        
        if tracking_uri:
            mlflow.set_tracking_uri(tracking_uri)
        
        mlflow.set_experiment(experiment_name)
    
    def __enter__(self):
        self._run = mlflow.start_run()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._run:
            mlflow.end_run()
    
    def start_run(self, run_name: Optional[str] = None) -> Any:
        """Start a new MLflow run."""
        self._run = mlflow.start_run(run_name=run_name)
        return self._run
    
    def end_run(self) -> None:
        """End the current MLflow run."""
        mlflow.end_run()
    
    def log_params(self, params: dict) -> None:
        """Log parameters to MLflow."""
        mlflow.log_params(params)
    
    def log_metrics(self, metrics: dict, step: Optional[int] = None) -> None:
        """Log metrics to MLflow."""
        mlflow.log_metrics(metrics, step=step)
    
    def log_model(self, model: Any, name: str, **kwargs) -> None:
        """Log a model to MLflow."""
        mlflow.sklearn.log_model(model, name, **kwargs)
    
    def log_artifact(self, local_path: str, artifact_path: Optional[str] = None) -> None:
        """Log an artifact to MLflow."""
        mlflow.log_artifact(local_path, artifact_path)
    
    def log_dataframe(self, df: pd.DataFrame, name: str) -> None:
        """Log a DataFrame as an artifact."""
        import tempfile
        import os
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, f"{name}.parquet")
            df.to_parquet(path)
            mlflow.log_artifact(path)


class CrossValidator:
    """
    K-fold cross-validation for model evaluation.
    """
    
    def __init__(self, n_splits: int = 5, random_state: int = 42):
        self.n_splits = n_splits
        self.random_state = random_state
    
    def cross_validate(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        model_class,
        model_params: dict,
    ) -> dict:
        """
        Perform k-fold cross-validation.
        
        Args:
            X: Feature DataFrame
            y: Target Series
            model_class: Model class to instantiate
            model_params: Parameters for model
        
        Returns:
            Dictionary with cross-validation results
        """
        from sklearn.model_selection import cross_val_score
        
        model = model_class(**model_params)
        
        scores = cross_val_score(
            model, X, y, 
            cv=self.n_splits,
            scoring='accuracy',
        )
        
        return {
            'scores': scores.tolist(),
            'mean': scores.mean(),
            'std': scores.std(),
            'n_splits': self.n_splits,
        }


class ModelEvaluator:
    """
    Advanced model evaluation with cross-validation.
    """
    
    def __init__(self, config: Optional['ModelConfig'] = None):
        self.config = config or ModelConfig()
        self.cross_validator = CrossValidator(
            n_splits=self.config.k_folds,
            random_state=self.config.random_state,
        )
    
    def evaluate_with_cv(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        model_class,
    ) -> dict:
        """Evaluate model with cross-validation."""
        model_params = self._get_model_params(model_class)
        cv_results = self.cross_validator.cross_validate(X, y, model_class, model_params)
        return cv_results
    
    def _get_model_params(self, model_class) -> dict:
        """Get appropriate params based on model class."""
        from xgboost import XGBClassifier
        from sklearn.ensemble import RandomForestRegressor
        
        if issubclass(model_class, XGBClassifier):
            return self.config.recall_hyperparams
        elif issubclass(model_class, RandomForestRegressor):
            return self.config.composition_hyperparams
        return {}
