"""
Tests for deployment/model_evaluation/models.py
"""
import pytest
import pandas as pd
import numpy as np
import os
import tempfile
import shutil
import json
import logging
from unittest.mock import Mock, patch, mock_open
from pathlib import Path


class TestModelTrainer:
    """Tests for ModelTrainer class."""

    def test_initialization(self, tmp_path):
        """Verify trainer initializes with correct config, ncbi_wrapper, dataframes."""
        from deployment.model_evaluation.models import ModelTrainer
        from deployment.model_evaluation.config import EvaluatorConfig

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n")
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        config = EvaluatorConfig(
            study_output_filepath=str(study_dir),
            taxid_plan_filepath=str(plan_file),
            analysis_output_filepath=str(output_dir)
        )

        ncbi_wrapper = Mock()
        input_tax_df = pd.DataFrame({'taxid': [9606], 'order': ['Primates']})
        taxids_to_use = pd.DataFrame({'taxid': [9606], 'order': ['Primates']})

        trainer = ModelTrainer(
            config=config,
            ncbi_wrapper=ncbi_wrapper,
            input_tax_df=input_tax_df,
            taxids_to_use=taxids_to_use
        )

        assert trainer.config == config
        assert trainer.ncbi == ncbi_wrapper
        assert trainer.input_tax_df.equals(input_tax_df)
        assert trainer.taxids_to_use.equals(taxids_to_use)
        assert trainer.recall_modeller is None
        assert trainer.composition_modeller is None
        assert trainer.crosshit_modeller is None

    def test_run_data_retrieval_empty_folders(self, tmp_path):
        """Handle empty training folder list."""
        from deployment.model_evaluation.models import ModelTrainer
        from deployment.model_evaluation.config import EvaluatorConfig

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n")
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        config = EvaluatorConfig(
            study_output_filepath=str(study_dir),
            taxid_plan_filepath=str(plan_file),
            analysis_output_filepath=str(output_dir)
        )

        trainer = ModelTrainer(
            config=config,
            ncbi_wrapper=Mock(),
            input_tax_df=pd.DataFrame(),
            taxids_to_use=pd.DataFrame()
        )

        training_df, prediction_df, recall_df = trainer.run_data_retrieval([])

        assert training_df.empty
        assert prediction_df.empty
        assert recall_df.empty

    def test_run_data_retrieval_skips_invalid_datasets(self, tmp_path):
        """Skip datasets without mapped reads."""
        from deployment.model_evaluation.models import ModelTrainer
        from deployment.model_evaluation.config import EvaluatorConfig

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n")
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        config = EvaluatorConfig(
            study_output_filepath=str(study_dir),
            taxid_plan_filepath=str(plan_file),
            analysis_output_filepath=str(output_dir)
        )

        with patch('deployment.model_evaluation.models.OverlapManager') as mock_om:
            mock_om_instance = Mock()
            mock_om_instance.m_stats_matrix = pd.DataFrame()
            mock_om.return_value = mock_om_instance

            trainer = ModelTrainer(
                config=config,
                ncbi_wrapper=Mock(),
                input_tax_df=pd.DataFrame(),
                taxids_to_use=pd.DataFrame()
            )

            training_df, prediction_df, recall_df = trainer.run_data_retrieval(['dataset1'])

            assert training_df.empty

    def test_train_models_empty_data_raises(self, tmp_path):
        """Raise ModelError when no training data."""
        from deployment.model_evaluation.models import ModelTrainer
        from deployment.model_evaluation.config import EvaluatorConfig
        from deployment.model_evaluation.exceptions import ModelError

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n")
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        config = EvaluatorConfig(
            study_output_filepath=str(study_dir),
            taxid_plan_filepath=str(plan_file),
            analysis_output_filepath=str(output_dir)
        )

        with patch('deployment.model_evaluation.models.OverlapManager') as mock_om, \
             patch('deployment.model_evaluation.models.data_set_traversal_with_precision') as mock_traversal, \
             patch('deployment.model_evaluation.models.cross_hit_prediction_matrix') as mock_cross, \
             patch('deployment.model_evaluation.models.predict_recall_cutoff_vars') as mock_recall, \
             patch('metagenomics_utils.overlap_manager.node_stats.get_m_stats_matrix') as mock_get_m:

            mock_om_instance = Mock()
            mock_om_instance.m_stats_matrix = pd.DataFrame({'leaf': ['l1']})
            mock_om.return_value = mock_om_instance
            mock_traversal.return_value = pd.DataFrame()
            mock_cross.return_value = pd.DataFrame()
            mock_get_m.return_value = pd.DataFrame()
            mock_recall.return_value = pd.DataFrame()

            trainer = ModelTrainer(
                config=config,
                ncbi_wrapper=Mock(),
                input_tax_df=pd.DataFrame(),
                taxids_to_use=pd.DataFrame()
            )

            with pytest.raises(ModelError) as exc_info:
                trainer.train_models(['dataset1'])

            assert 'composition' in str(exc_info.value).lower()

    def test_save_models_creates_directory(self, tmp_path):
        """Create output directory if not exists."""
        from deployment.model_evaluation.models import ModelTrainer
        from deployment.model_evaluation.config import EvaluatorConfig

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n")
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        config = EvaluatorConfig(
            study_output_filepath=str(study_dir),
            taxid_plan_filepath=str(plan_file),
            analysis_output_filepath=str(output_dir)
        )

        trainer = ModelTrainer(
            config=config,
            ncbi_wrapper=Mock(),
            input_tax_df=pd.DataFrame({'taxid': [9606]}),
            taxids_to_use=pd.DataFrame({'taxid': [9606]})
        )

        with tempfile.TemporaryDirectory() as tmptest_dir:
            trainer.save_models(tmptest_dir)
            assert os.path.exists(tmptest_dir)

    def test_get_trained_models_requires_training(self, tmp_path):
        """Raise error if models not trained."""
        from deployment.model_evaluation.models import ModelTrainer
        from deployment.model_evaluation.config import EvaluatorConfig
        from deployment.model_evaluation.exceptions import ModelError

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n")
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        config = EvaluatorConfig(
            study_output_filepath=str(study_dir),
            taxid_plan_filepath=str(plan_file),
            analysis_output_filepath=str(output_dir)
        )

        trainer = ModelTrainer(
            config=config,
            ncbi_wrapper=Mock(),
            input_tax_df=pd.DataFrame(),
            taxids_to_use=pd.DataFrame()
        )

        with pytest.raises(ModelError) as exc_info:
            trainer.get_trained_models()

        assert 'not trained' in str(exc_info.value).lower()


class TestModelRegistry:
    """Tests for ModelRegistry class."""

    def test_load_models_not_found_raises(self, tmp_path):
        """Raise FileNotFoundError for unknown model."""
        from deployment.model_evaluation.models import ModelRegistry

        registry = ModelRegistry(base_path=str(tmp_path / "models"))

        with pytest.raises(FileNotFoundError):
            registry.load_models("nonexistent_model")

    def test_list_models_empty(self, tmp_path):
        """Return empty list when no metadata file."""
        from deployment.model_evaluation.models import ModelRegistry

        registry = ModelRegistry(base_path=str(tmp_path / "models"))

        assert registry.list_models() == []

    def test_registry_initialization(self, tmp_path):
        """Test registry initialization."""
        from deployment.model_evaluation.models import ModelRegistry

        registry = ModelRegistry(base_path=str(tmp_path / "models"))

        assert "models" in str(registry.base_path)


class TestMLflowTracker:
    """Tests for MLflowTracker class."""

    def test_init_requires_mlflow(self):
        """Raise ImportError when mlflow not installed."""
        from deployment.model_evaluation.models import MLflowTracker

        with patch('deployment.model_evaluation.models.MLFLOW_AVAILABLE', False):
            with pytest.raises(ImportError):
                MLflowTracker()

    @patch('deployment.model_evaluation.models.mlflow')
    def test_context_manager(self, mock_mlflow):
        """Test __enter__/__exit__ methods."""
        from deployment.model_evaluation.models import MLflowTracker

        mock_run = Mock()
        mock_mlflow.start_run.return_value = mock_run

        tracker = MLflowTracker(experiment_name="test")

        with tracker as t:
            assert t == tracker

        mock_mlflow.start_run.assert_called_once()
        mock_mlflow.end_run.assert_called_once()

    @patch('deployment.model_evaluation.models.mlflow')
    def test_start_run(self, mock_mlflow):
        """Start new MLflow run."""
        from deployment.model_evaluation.models import MLflowTracker

        mock_run = Mock()
        mock_mlflow.start_run.return_value = mock_run

        tracker = MLflowTracker()
        result = tracker.start_run("test_run")

        mock_mlflow.start_run.assert_called_once_with(run_name="test_run")
        assert result == mock_run

    @patch('deployment.model_evaluation.models.mlflow')
    def test_log_params(self, mock_mlflow):
        """Log parameters to MLflow."""
        from deployment.model_evaluation.models import MLflowTracker

        tracker = MLflowTracker()
        tracker.log_params({'learning_rate': 0.01, 'n_estimators': 100})

        mock_mlflow.log_params.assert_called_once_with({'learning_rate': 0.01, 'n_estimators': 100})

    @patch('deployment.model_evaluation.models.mlflow')
    def test_log_metrics(self, mock_mlflow):
        """Log metrics with optional step."""
        from deployment.model_evaluation.models import MLflowTracker

        tracker = MLflowTracker()
        tracker.log_metrics({'accuracy': 0.95}, step=1)

        mock_mlflow.log_metrics.assert_called_once_with({'accuracy': 0.95}, step=1)

    @patch('deployment.model_evaluation.models.mlflow')
    def test_log_model(self, mock_mlflow):
        """Log sklearn model."""
        from deployment.model_evaluation.models import MLflowTracker

        tracker = MLflowTracker()
        mock_model = Mock()
        tracker.log_model(mock_model, "my_model")

        mock_mlflow.sklearn.log_model.assert_called_once_with(mock_model, "my_model")

    @patch('deployment.model_evaluation.models.mlflow')
    def test_log_dataframe(self, mock_mlflow, tmp_path):
        """Log DataFrame as parquet artifact."""
        from deployment.model_evaluation.models import MLflowTracker

        tracker = MLflowTracker()
        df = pd.DataFrame({'a': [1, 2], 'b': [3, 4]})

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test_df.parquet")
            df.to_parquet(path)
            mock_mlflow.log_artifact.reset_mock()
            mock_mlflow.log_artifact.return_value = None
            tracker.log_dataframe(df, "test_df")

        mock_mlflow.log_artifact.assert_called_once()


class TestCrossValidator:
    """Tests for CrossValidator class."""

    def test_cross_validate_returns_scores(self):
        """Return dict with scores, mean, std."""
        from sklearn.model_selection import cross_val_score
        from deployment.model_evaluation.models import CrossValidator
        from sklearn.ensemble import RandomForestClassifier

        validator = CrossValidator(n_splits=3, random_state=42)

        X = pd.DataFrame({'feature': [1, 2, 3, 4, 5, 6]})
        y = pd.Series([0, 1, 0, 1, 0, 1])

        with patch.object(__import__('sklearn.model_selection', fromlist=['cross_val_score']), 'cross_val_score') as mock_cv:
            mock_cv.return_value = np.array([0.8, 0.9, 0.85])

            result = validator.cross_validate(X, y, RandomForestClassifier, {'n_estimators': 10})

            assert 'scores' in result
            assert 'mean' in result
            assert 'std' in result
            assert 'n_splits' in result
            assert result['mean'] == pytest.approx(0.85)

    def test_cross_validate_with_xgboost(self):
        """Test with XGBClassifier."""
        from deployment.model_evaluation.models import CrossValidator

        validator = CrossValidator(n_splits=3)

        X = pd.DataFrame({'feature': [1, 2, 3, 4, 5, 6]})
        y = pd.Series([0, 1, 0, 1, 0, 1])

        with patch('sklearn.model_selection.cross_val_score') as mock_cv:
            mock_cv.return_value = np.array([0.8, 0.9, 0.85])

            from xgboost import XGBClassifier
            result = validator.cross_validate(X, y, XGBClassifier, {'n_estimators': 10})

            assert result['mean'] == pytest.approx(0.85)


class TestModelEvaluator:
    """Tests for ModelEvaluator class."""

    def test_evaluate_with_cv(self):
        """Test cross-validation evaluation."""
        from deployment.model_evaluation.models import ModelEvaluator
        from sklearn.ensemble import RandomForestClassifier

        evaluator = ModelEvaluator()

        X = pd.DataFrame({'feature': [1, 2, 3, 4, 5, 6]})
        y = pd.Series([0, 1, 0, 1, 0, 1])

        with patch.object(evaluator.cross_validator, 'cross_validate') as mock_cv:
            mock_cv.return_value = {'scores': [0.8], 'mean': 0.8, 'std': 0.1, 'n_splits': 3}

            result = evaluator.evaluate_with_cv(X, y, RandomForestClassifier)

            assert result['mean'] == 0.8

    def test_get_model_params_xgboost(self):
        """Return correct params for XGBoost."""
        from deployment.model_evaluation.models import ModelEvaluator
        from xgboost import XGBClassifier

        evaluator = ModelEvaluator()

        params = evaluator._get_model_params(XGBClassifier)

        assert 'n_estimators' in params

    def test_get_model_params_rf(self):
        """Return correct params for RandomForest."""
        from deployment.model_evaluation.models import ModelEvaluator
        from sklearn.ensemble import RandomForestRegressor

        evaluator = ModelEvaluator()

        params = evaluator._get_model_params(RandomForestRegressor)

        assert 'n_estimators' in params


class TestTrainedModels:
    """Tests for TrainedModels class."""

    def test_trained_models_initialization(self):
        """Test TrainedModels initialization."""
        from deployment.model_evaluation.config import TrainedModels

        recall = Mock()
        composition = Mock()
        crosshit = Mock()

        models = TrainedModels(
            recall_modeller=recall,
            composition_modeller=composition,
            crosshit_modeller=crosshit
        )

        assert models.recall_modeller == recall
        assert models.composition_modeller == composition
        assert models.crosshit_modeller == crosshit

    def test_is_complete(self):
        """Test is_complete property."""
        from deployment.model_evaluation.config import TrainedModels

        models = TrainedModels(
            recall_modeller=Mock(),
            composition_modeller=Mock(),
            crosshit_modeller=Mock()
        )

        assert models.is_complete is True

        models_incomplete = TrainedModels(
            recall_modeller=Mock(),
            composition_modeller=None,
            crosshit_modeller=Mock()
        )

        assert models_incomplete.is_complete is False

    def test_repr(self):
        """Test __repr__ method."""
        from deployment.model_evaluation.config import TrainedModels

        models = TrainedModels(
            recall_modeller=Mock(),
            composition_modeller=None,
            crosshit_modeller=Mock()
        )

        repr_str = repr(models)
        assert 'recall=True' in repr_str
        assert 'composition=False' in repr_str
        assert 'crosshit=True' in repr_str


class TestTrainAndEvaluate:
    """Tests for train_and_evaluate function."""

    @patch('deployment.model_evaluation.models.ModelTrainer')
    def test_train_and_evaluate_full_pipeline(self, mock_trainer_class, tmp_path):
        """Test complete training pipeline with mocks."""
        from deployment.model_evaluation.models import train_and_evaluate
        from deployment.model_evaluation.config import EvaluatorConfig

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n")
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        config = EvaluatorConfig(
            study_output_filepath=str(study_dir),
            taxid_plan_filepath=str(plan_file),
            analysis_output_filepath=str(output_dir)
        )

        mock_trainer_instance = Mock()
        mock_trainer_class.return_value = mock_trainer_instance

        mock_recall = Mock()
        mock_composition = Mock()
        mock_crosshit = Mock()
        mock_trainer_instance.get_trained_models.return_value = (mock_recall, mock_composition, mock_crosshit)

        with tempfile.TemporaryDirectory() as tmptest_dir:
            result = train_and_evaluate(
                config=config,
                training_folders=['dataset1'],
                ncbi_wrapper=Mock(),
                input_tax_df=pd.DataFrame(),
                taxids_to_use=pd.DataFrame(),
                models_output_dir=tmptest_dir
            )

            mock_trainer_instance.train_models.assert_called_once_with(['dataset1'])
            mock_trainer_instance.save_models.assert_called_once_with(tmptest_dir)
            mock_trainer_instance.evaluate_models.assert_called_once_with(tmptest_dir)
            mock_trainer_instance.get_trained_models.assert_called_once()

            assert result == (mock_recall, mock_composition, mock_crosshit)
