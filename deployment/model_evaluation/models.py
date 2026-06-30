"""
Model training module for evaluator.

Handles training of RecallModeller, BaseCompositionModeller (5 variants), and CrossHitModeller.
"""

import os
import logging
from pathlib import Path
from typing import Tuple, Any, Optional, List
import joblib
from metagenomics_utils.overlap_manager.node_stats import get_m_stats_matrix
from .result_models import DatasetResult

import pandas as pd

from .config import TrainedModels, ModelConfig

from metagenomics_utils.overlap_manager import OverlapManager
from metagenomics_utils.overlap_manager.om_models import (
    RecallModeller,
    InjectModellerInterface,
    BaseCompositionModeller,
    XGBCompositionModeller,
    OptunaXGBCompositionModeller,
    RFCompositionModeller,
    GBCompositionModeller,
    LRCompositionModeller,
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
    - BaseCompositionModeller (5 variants via --composition_model_interface): Predicts clade composition
    - CrossHitModeller: Predicts cross-hit probability
    """
    
    def __init__(
        self,
        config: EvaluatorConfig,
        ncbi_wrapper: Any,
        input_tax_df: pd.DataFrame,
        taxids_to_use: pd.DataFrame,
        use_cache: bool = True,
    ):
        """
        Initialize model trainer.
        
        Args:
            config: EvaluatorConfig instance
            ncbi_wrapper: NCBI TaxonomistWrapper instance
            input_tax_df: Input tax DataFrame
            taxids_to_use: Taxids to use for evaluation
            use_cache: Whether to use cached training data (default: True)
        """
        self.config = config
        self.ncbi = ncbi_wrapper
        self.input_tax_df = input_tax_df
        self.taxids_to_use = taxids_to_use
        self.use_cache = use_cache
        
        self.recall_modeller: Optional[RecallModeller] = None
        self.composition_modeller: Optional[BaseCompositionModeller] = None
        self.crosshit_modeller: Optional[CrossHitModeller] = None
        
        self._training_results: Optional[pd.DataFrame] = None
        self._prediction_results: Optional[pd.DataFrame] = None
        self._recall_results: Optional[pd.DataFrame] = None
        self._recall_matrices: Optional[List[pd.DataFrame]] = None

        self.models: Optional[TrainedModels] = None
    
    @property
    def cache_dir(self) -> Path:
        """Return the cache directory path."""
        return self.config.models_dir / "cache"
    
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
        stats_matrices = []
        recall_matrices = []
        
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

                m_stats_stats_matrix = self._get_m_stats_matrix(data_set_name, overlap_manager, filter_no_leaf=False)

                if m_stats_stats_matrix.empty:
                    continue
                
                # Processed features (for cache / composition models)
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
                if m_stats_stats_matrix is not None and not m_stats_stats_matrix.empty:
                    stats_matrices.append(m_stats_stats_matrix)
                    recall_matrices.append(m_stats_stats_matrix)

            except Exception as e:
                logger.error(f"Error processing {data_set_name}: {e}")
                import traceback
                traceback.print_exc()
                continue


        stats_matrices = pd.concat(stats_matrices, ignore_index=True) if stats_matrices else pd.DataFrame()
        training_results_df = pd.concat(training_results, ignore_index=True) if training_results else pd.DataFrame()
        prediction_results_df = pd.concat(prediction_results, ignore_index=True) if prediction_results else pd.DataFrame()
        recall_results_df = pd.concat(recall_results, ignore_index=True) if recall_results else pd.DataFrame()
        
        self._training_results = training_results_df
        self._prediction_results = prediction_results_df
        self._recall_results = recall_results_df
        self._stats_matrices = stats_matrices
        self._recall_matrices = recall_matrices

        logger.info(f"Retrieved training data: {len(training_results_df)} results, {len(prediction_results_df)} predictions")
        logger.info(f"Raw recall matrices: {len(recall_matrices)} datasets")
        
        return training_results_df, prediction_results_df, recall_results_df

    def _get_m_stats_matrix(self, data_set_name: str, overlap_manager: OverlapManager, filter_no_leaf: bool= False) -> pd.DataFrame:
        """
        Get m_stats matrix for a dataset.
        
        Args:
            data_set_name: Name of dataset
            overlap_manager: OverlapManager instance
            
        Returns:
            m_stats matrix DataFrame
        """
        return get_m_stats_matrix(
            data_set_name,
            self.config.study_output_filepath,
            self.ncbi,
            overlap_manager,
            filter_no_leaf=filter_no_leaf
        )
    
    def _ensure_cache_dir(self) -> Path:
        """Ensure the cache directory exists."""
        cache_dir = self.cache_dir
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir
    
    def save_cached_data(
        self,
        training_df: pd.DataFrame,
        prediction_df: pd.DataFrame,
        recall_df: pd.DataFrame,
    ) -> None:
        """
        Save parsed training data to cache for future use.
        
        Args:
            training_df: Training results DataFrame
            prediction_df: Prediction results DataFrame
            recall_df: Recall results DataFrame
        """
        cache_dir = self._ensure_cache_dir()
        
        training_path = cache_dir / "training_results_cache.parquet"
        prediction_path = cache_dir / "prediction_results_cache.parquet"
        recall_path = cache_dir / "recall_results_cache.parquet"
        taxid_path = cache_dir / "taxids_to_use_cache.parquet"

        training_df.to_parquet(training_path, index=False)
        prediction_df.to_parquet(prediction_path, index=False)
        recall_df.to_parquet(recall_path, index=False)
        self.taxids_to_use.to_parquet(taxid_path, index=False)

        if self._recall_matrices:
            matrices_path = cache_dir / "recall_matrices_cache.joblib"
            joblib.dump(self._recall_matrices, matrices_path)

        logger.info(f"Saved cached training data to {cache_dir}")

    def load_cached_data(self) -> Optional[Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, List[pd.DataFrame]]]:
        """
        Load cached training data if available.
        
        Returns:
            Tuple of (training_df, prediction_df, recall_df, taxids_df, recall_matrices) if cache exists, None otherwise
        """
        if not self.use_cache:
            
            logger.info("Cache disabled, skipping load")
            return None
        
        cache_dir = self.cache_dir
        
        training_path = cache_dir / "training_results_cache.parquet"
        prediction_path = cache_dir / "prediction_results_cache.parquet"
        recall_path = cache_dir / "recall_results_cache.parquet"
        taxids_path = cache_dir / "taxids_to_use_cache.parquet"
        matrices_path = cache_dir / "recall_matrices_cache.joblib"

        if not all(p.exists() for p in [training_path, prediction_path, recall_path, taxids_path, matrices_path]):
            logger.info("Cache not found, will compute fresh data")
            return None
        
        try:
            training_df = pd.read_parquet(training_path)
            prediction_df = pd.read_parquet(prediction_path)
            recall_df = pd.read_parquet(recall_path)
            taxids_df = pd.read_parquet(taxids_path)
            recall_matrices = joblib.load(matrices_path)

            logger.info(f"Loaded cached training data: {len(training_df)} training, {len(prediction_df)} prediction, {len(recall_df)} recall, {len(taxids_df)} taxids records, {len(recall_matrices)} matrices")
            return training_df, prediction_df, recall_df, taxids_df, recall_matrices 
        except Exception as e:
            logger.warning(f"Failed to load cache: {e}, will compute fresh data")
            return None
    
    def clear_cache(self) -> None:
        """Delete cached training data."""
        cache_dir = self.cache_dir
        if cache_dir.exists():
            for cache_file in cache_dir.glob("*_cache*"):
                cache_file.unlink()
            logger.info(f"Cleared cache at {cache_dir}")
        else:
            logger.info("No cache to clear")
    
    def _init_composition_modeller(self, training_df: pd.DataFrame) -> BaseCompositionModeller:
        """Factory: create composition modeller based on config."""
        interface = self.config.composition_model_interface
        logger.info(f"Composition model interface: {interface}")
        if interface == 'xgb':
            return XGBCompositionModeller()
        elif interface == 'xgb_optimized':
            return OptunaXGBCompositionModeller()
        elif interface == 'rf':
            return RFCompositionModeller()
        elif interface == 'gb':
            return GBCompositionModeller()
        elif interface == 'lr':
            return LRCompositionModeller()
        else:
            logger.warning(f"Unknown composition model interface '{interface}', falling back to xgb")
            return XGBCompositionModeller()
    
    def train_models(self, training_folders: list, force_refresh: bool = False) -> Tuple[Any, Any]:
        logger.info(f"Starting model training with {len(training_folders)} datasets")

        if force_refresh:
            logger.info("Force refresh enabled, clearing cache")
            self.clear_cache()

        cached_data = self.load_cached_data()

        if cached_data is not None:
            training_df, prediction_df, recall_df, taxids_df, recall_matrices = cached_data
            logger.info("Using cached training data")
            self._recall_matrices = recall_matrices
        else:
            logger.info("Computing fresh training data")
            training_df, prediction_df, recall_df = self.run_data_retrieval(training_folders)
            if not training_df.empty:
                self.save_cached_data(training_df, prediction_df, recall_df)
            taxids_df = self.taxids_to_use

        self.taxids_to_use = taxids_df

        if training_df.empty:
            raise ModelError("composition", "train", "No training data available")

        if self.config.recall_model_interface == 'gp_clf':
            from metagenomics_utils.overlap_manager.om_models import GPCLFRecallModeller
            model_interface = InjectModellerInterface(model_type='gp_clf')
            self.recall_modeller = GPCLFRecallModeller(
                data_set_divide=self.config.data_set_divide,
                model_interface=model_interface,
                sort_strategy=self.config.recall_sort_strategy,
            )
            logger.info(f"GP-CLF mode: using per-division GPs with CLF threshold, data_set_divide={self.config.data_set_divide}")
        elif self.config.recall_model_interface == 'direct':
            from metagenomics_utils.overlap_manager.om_models import CutoffRecallModeller
            self.recall_modeller = CutoffRecallModeller(
                data_set_divide=self.config.data_set_divide,
                target_recall=self.config.target_recall,
                sort_strategy=self.config.recall_sort_strategy,
            )
        elif self.config.recall_model_interface == 'direct_xgb':
            from metagenomics_utils.overlap_manager.om_models import DirectXGBRecallModeller
            self.recall_modeller = DirectXGBRecallModeller(
                data_set_divide=self.config.data_set_divide,
                target_recall=self.config.target_recall,
                sort_strategy=self.config.recall_sort_strategy,
            )
            logger.info(f"DirectXGB mode: direct fraction regression with XGBoost, data_set_divide={self.config.data_set_divide}")
        else:
            model_interface = InjectModellerInterface(model_type=self.config.recall_model_interface)
            self.recall_modeller = RecallModeller(
                data_set_divide=self.config.data_set_divide,
                model_interface=model_interface,
                sort_strategy=self.config.recall_sort_strategy,
            )

        self.composition_modeller = self._init_composition_modeller(training_df)

        self.crosshit_modeller = CrossHitModeller(
            prediction_trainning_results_df=prediction_df,
        )

        logger.info("Training recall model...")
        model_recall = self.recall_modeller.fit(
            self._recall_matrices, taxids_df,
        )

        logger.info("Training composition model...")
        DROPPED_COLS = ['data_set', 'node', 'n_true_leaves', 'precision_increased',
                        'new_precision', 'precision', 'stop_traversal', 'unclassified']
        available_drop = [c for c in DROPPED_COLS if c in training_df.columns]
        X = training_df.drop(columns=available_drop)
        y = training_df['stop_traversal'].astype(int)
        from sklearn.model_selection import train_test_split
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42, stratify=y,
        )
        self.composition_modeller.fit(X_train, y_train, X_test, y_test)

        logger.info("Training crosshit model...")
        self.crosshit_modeller.train_model()

        self.models = TrainedModels(
            recall_modeller=self.recall_modeller,
            composition_modeller=self.composition_modeller,
            crosshit_modeller=self.crosshit_modeller,
            taxids_to_use=self.taxids_to_use,
        )

        return model_recall, self.crosshit_modeller
    
    def save_models(self, output_dir: str) -> None:
        """
        Save trained models to disk.
        
        Args:
            output_dir: Directory to save models
        """
        import pandas as pd
        
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
        
        if self.taxids_to_use is not None:
            taxids_path = os.path.join(output_dir, "taxids_to_use.parquet")
            self.taxids_to_use.to_parquet(taxids_path, index=False)
            logger.info(f"Saved taxids_to_use ({len(self.taxids_to_use)} taxids) to {taxids_path}")
    
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
            if hasattr(self.recall_modeller, 'plot_diagnostics'):
                logger.info("Generating recall diagnostic plots...")
                self.recall_modeller.plot_diagnostics(output_dir)
        
        if self.composition_modeller is not None:
            logger.info("Evaluating composition model...")
            self.composition_modeller.eval_and_plot(
                self.composition_modeller.X_test,
                self.composition_modeller.y_test,
                output_dir,
                X_train=self.composition_modeller.X_train
            )
    
    def get_trained_models(self) -> Tuple[RecallModeller, BaseCompositionModeller, CrossHitModeller]:
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
) -> Tuple[RecallModeller, BaseCompositionModeller, CrossHitModeller]:
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
        import pandas as pd
        
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
        
        taxids_path = model_path / "taxids_to_use.parquet"
        if taxids_path.exists():
            models.taxids_to_use = pd.read_parquet(taxids_path)
            logger.info(f"Loaded taxids_to_use ({len(models.taxids_to_use)} taxids) from {taxids_path}")
        
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


class TrainingCrossHitAnalyzer:
    """
    Processes training datasets to compute cross-hit metrics.
    """
    
    def __init__(
        self,
        config: 'EvaluatorConfig',
        ncbi_wrapper,
        input_tax_df,
        taxids_to_use,
    ):
        self.config = config
        self.ncbi = ncbi_wrapper
        self.input_tax_df = input_tax_df
        self.taxids_to_use = taxids_to_use

        self.logger = logging.getLogger(__name__)

    def analyze_training_data(self, training_folders: list) -> list:
        """
        Process all training folders and return DatasetResult objects.
        
        Returns:
            List of DatasetResult with cross-hit metrics populated
        """
        results = []
        for folder in training_folders:
            try:
                result = self._process_dataset(folder)
                if result:
                    results.append(result)
            except Exception as e:
                self.logger.warning(f"Error processing {folder}: {e}")
                import traceback
                traceback.print_exc()
        return results
    
    def _process_dataset(self, data_set_name: str):
        """Process single training dataset for cross-hit metrics only."""
        from metagenomics_utils.overlap_manager import OverlapManager
        
        om_path = os.path.join(self.config.study_output_filepath, data_set_name, "clustering")
        input_path = os.path.join(self.config.study_output_filepath, data_set_name, "input", f"{data_set_name}.tsv")
        
        if not os.path.exists(input_path):
            return None
        
        try:
            overlap_manager = OverlapManager(om_path)
        except Exception:
            return None
        
        if overlap_manager.m_stats_matrix.empty:
            return None
        
        input_df = pd.read_csv(input_path, sep="\t")
        input_summary = input_df[['sample', 'taxid', 'reads', 'mutation_rate']].drop_duplicates()
        input_summary = input_summary.merge(self.input_tax_df, on='taxid', how='left')
        
        result = DatasetResult(data_set=data_set_name, input_df=input_summary)
        result = self._compute_cross_hit_metrics(overlap_manager, result)
        return result

    def _compute_cross_hit_metrics(self, overlap_manager, result: DatasetResult):
        """Compute cross-hit related metrics."""

        m_stats = get_m_stats_matrix(
            result.data_set,
            self.config.study_output_filepath,
            self.ncbi,
            overlap_manager,
            filter_no_leaf=False
        )
        
        cross_hit_mask = m_stats['is_crosshit'] == True
        cross_hit_subset = m_stats[cross_hit_mask]
        tax_level = self.config.tax_level
        
        if not cross_hit_subset.empty:
            merged = cross_hit_subset[['best_match_taxid', 'numreads']].merge(
                result.input_df,
                left_on='best_match_taxid',
                right_on='taxid',
                how='left'
            )
            result.cross_hit.total_true_cross_hits = int(merged.shape[0])
            result.cross_hit.total_cross_hit_reads_mapped = int(merged['numreads'].sum())
            
            if tax_level in merged.columns:
                cross_hit_by_class = merged.groupby(tax_level)['best_match_taxid'].count().reset_index()
                cross_hit_by_class.columns = ['class', 'count']
                result.cross_hit.cross_hit_counts_per_class = cross_hit_by_class.to_dict(orient='records')
                
                cross_hit_reads_by_class = merged.groupby(tax_level)['numreads'].sum().reset_index()
                cross_hit_reads_by_class.columns = ['class', 'reads']
                result.cross_hit.cross_hit_reads_per_class = cross_hit_reads_by_class.to_dict(orient='records')
        
        spurious_mask = m_stats['is_trash'] == True
        spurious_subset = m_stats[spurious_mask]
        spurious_subset = spurious_subset[['best_match_taxid', 'numreads']].merge(
            result.input_df,
            left_on='best_match_taxid',
            right_on='taxid',
            how='left'
        )
        if not spurious_subset.empty and tax_level in spurious_subset.columns:
            spurious_by_class = spurious_subset.groupby(tax_level)['best_match_taxid'].count().reset_index()
            spurious_by_class.columns = ['class', 'count']
            result.cross_hit.spurious_hit_counts_per_class = spurious_by_class.to_dict(orient='records')

            spurious_reads_by_class = spurious_subset.groupby(tax_level)['numreads'].sum().reset_index()
            spurious_reads_by_class.columns = ['class', 'reads']
            result.cross_hit.spurious_hit_reads_per_class = spurious_reads_by_class.to_dict(orient='records')
        
        
        tax_level = self.config.tax_level
        reads_simulated_by_class = result.input_df.groupby(tax_level)['reads'].sum().reset_index()
        reads_simulated_by_class.columns = ['class', 'reads']
        result.reads_simulated_per_class = reads_simulated_by_class.to_dict(orient='records')
        
        return result
