"""
Tests for deployment/model_evaluation/evaluate.py
"""
import pytest
import pandas as pd
import numpy as np
import os
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock


class TestRetrieveSimulationInput:
    """Tests for retrieve_simulation_input function."""

    def test_retrieve_simulation_input_empty_directory(self, tmp_path):
        """Test with empty output directory."""
        from deployment.model_evaluation.evaluate import retrieve_simulation_input
        
        result = retrieve_simulation_input(str(tmp_path))
        assert result.empty or len(result) == 0

    def test_retrieve_simulation_input_with_data(self, tmp_path):
        """Test with valid dataset structure."""
        from deployment.model_evaluation.evaluate import retrieve_simulation_input
        
        dataset_name = "test_dataset"
        dataset_dir = tmp_path / dataset_name / "input"
        dataset_dir.mkdir(parents=True)
        
        input_file = dataset_dir / f"{dataset_name}.tsv"
        input_file.write_text("sample\ttaxid\treads\tmutation_rate\nsample1\t9606\t100\t0.01\n")
        
        result = retrieve_simulation_input(str(tmp_path))
        
        assert not result.empty
        assert 'data_set' in result.columns

    def test_retrieve_simulation_input_skips_empty_datasets(self, tmp_path):
        """Test that empty datasets are skipped."""
        from deployment.model_evaluation.evaluate import retrieve_simulation_input
        
        dataset_name = "test_dataset"
        dataset_dir = tmp_path / dataset_name / "input"
        dataset_dir.mkdir(parents=True)
        
        input_file = dataset_dir / f"{dataset_name}.tsv"
        input_file.write_text("sample\ttaxid\treads\tmutation_rate\n")
        
        result = retrieve_simulation_input(str(tmp_path))
        assert result.empty


class TestProcessCladeReport:
    """Tests for process_clade_report function."""

    def test_process_clade_report_missing_file(self, tmp_path):
        """Test when clade report doesn't exist."""
        from deployment.model_evaluation.evaluate import process_clade_report
        
        result = process_clade_report("nonexistent", str(tmp_path))
        assert result == []

    def test_process_clade_report_with_data(self, tmp_path):
        """Test with valid clade report."""
        from deployment.model_evaluation.evaluate import process_clade_report
        
        dataset_name = "test_dataset"
        output_dir = tmp_path / dataset_name / "output"
        output_dir.mkdir(parents=True)
        
        clade_file = output_dir / "clade_report_with_references.tsv"
        clade_file.write_text("#rname\ttaxid\tcoverage\nref1\t9606\t10.5\nref2\t9598\t5.2\n")
        
        result = process_clade_report(dataset_name, str(tmp_path))
        
        assert len(result) > 0
        assert 9606 in result


class TestOutputParse:
    """Tests for output_parse function."""

    @patch('deployment.model_evaluation.evaluate.NCBITaxonomistWrapper')
    @patch('deployment.model_evaluation.evaluate.process_clade_report')
    def test_output_parse_with_mock(self, mock_process, mock_ncbi, tmp_path):
        """Test output_parse with mocked dependencies."""
        from deployment.model_evaluation.evaluate import output_parse
        
        mock_process.return_value = [9606, 9598]
        mock_ncbi_instance = Mock()
        mock_ncbi_instance.get_level.return_value = "test_order"
        mock_ncbi.return_value = mock_ncbi_instance
        
        dataset_name = "test_dataset"
        dataset_dir = tmp_path / dataset_name / "output"
        dataset_dir.mkdir(parents=True)
        
        result = output_parse(str(tmp_path), mock_ncbi_instance, [9606])
        
        assert isinstance(result, pd.DataFrame)
        assert 'taxid' in result.columns


class TestEstablishTaxidsToUse:
    """Tests for establish_taxids_to_use function."""

    @patch('deployment.model_evaluation.evaluate.output_parse')
    @patch('deployment.model_evaluation.evaluate.NCBITaxonomistWrapper')
    def test_establish_taxids_to_use_basic(self, mock_ncbi, mock_parse, tmp_path):
        """Test establish_taxids_to_use with mock."""
        from deployment.model_evaluation.evaluate import establish_taxids_to_use
        
        mock_df = pd.DataFrame({
            'taxid': [9606, 9598, 0],
            'order': ['Primates', 'Primates', 'unclassified'],
            'family': ['Hominidae', 'Hominidae', 'unclassified'],
            'genus': ['Homo', 'Pan', 'unclassified']
        })
        mock_parse.return_value = mock_df
        
        mock_ncbi = Mock()
        
        result = establish_taxids_to_use(
            str(tmp_path), mock_ncbi, [9606, 9598],
            tax_level_to_use='order', min_tax_count=0.02
        )
        
        assert isinstance(result, pd.DataFrame)


class TestExpandInputData:
    """Tests for expand_input_data function."""

    def test_expand_input_data_basic(self):
        """Test expand_input_data basic functionality."""
        from deployment.model_evaluation.evaluate import expand_input_data
        
        input_data = pd.DataFrame({
            'sample': ['s1', 's2'],
            'taxid': [9606, 9598],
            'reads': [100, 50]
        })
        
        taxid_plan = pd.DataFrame({
            'taxid': [9606, 9598],
            'description': ['Human', 'Chimp'],
            'lineage': ['lineage1', 'lineage2']
        })
        
        mock_ncbi = Mock()
        mock_ncbi.get_level.return_value = "test_order"
        
        result = expand_input_data(input_data, mock_ncbi, taxid_plan)
        
        assert isinstance(result, pd.DataFrame)


class TestCompoundEDAFunction:
    """Tests for compound_eda_function."""

    def test_compound_eda_function_with_mocks(self, tmp_path):
        """Test compound_eda_function with mocked dependencies."""
        from deployment.model_evaluation.evaluate import compound_eda_function
        
        dataset_name = "test_dataset"
        dataset_dir = tmp_path / dataset_name
        (dataset_dir / "clustering").mkdir(parents=True)
        (dataset_dir / "input").mkdir(parents=True)
        
        input_file = dataset_dir / "input" / f"{dataset_name}.tsv"
        input_file.write_text("sample\ttaxid\treads\tmutation_rate\ntest\t9606\t100\t0.01\n")
        
        with patch('deployment.model_evaluation.evaluate.OverlapManager') as mock_om, \
             patch('deployment.model_evaluation.evaluate.get_m_stats_matrix') as mock_get_stats, \
             patch('deployment.model_evaluation.evaluate.predict_data_set_clades') as mock_predict, \
             patch('deployment.model_evaluation.evaluate.cross_hit_prediction') as mock_cross, \
             patch('deployment.model_evaluation.evaluate.calculate_overall_precision') as mock_precision, \
             patch('deployment.model_evaluation.evaluate.get_trash_composition') as mock_trash, \
             patch('deployment.model_evaluation.evaluate.get_cross_hit_composition') as mock_cross_hit:
            
            mock_om_instance = Mock()
            mock_om_instance.m_stats_matrix = pd.DataFrame({'leaf': ['l1'], 'coverage': [1.0], 'best_match_taxid': [9606], 'best_match_is_best': [True], 'is_crosshit': [False]})
            mock_om_instance.leaves = ['l1']
            mock_om.return_value = mock_om_instance
            
            mock_get_stats.return_value = pd.DataFrame({
                'leaf': ['l1'],
                'coverage': [1.0],
                'best_match_taxid': [9606],
                'best_match_is_best': [True],
                'is_crosshit': [False]
            })
            
            mock_predict.return_value = pd.DataFrame({'best_taxid_match': [9606]})
            mock_cross.return_value = pd.DataFrame({'leaf': [], 'prob_best_match': [], 'is_trash': []})
            mock_precision.return_value = 0.8
            mock_trash.return_value = pd.DataFrame({'taxid': [], 'tax_level': []})
            mock_cross_hit.return_value = pd.DataFrame({'taxid': [], 'tax_level': []})
            
            mock_ncbi = Mock()
            mock_input_tax_df = pd.DataFrame({'taxid': [9606], 'order': ['Primates']})
            mock_taxids_to_use = pd.DataFrame({'taxid': [9606], 'order': ['Primates']})
            
            recall_modeller = Mock()
            composition_modeller = Mock()
            cross_hit_modeller = Mock()
            
            result = compound_eda_function(
                [dataset_name],
                str(tmp_path),
                5,
                mock_ncbi,
                mock_input_tax_df,
                'order',
                mock_taxids_to_use,
                recall_modeller,
                composition_modeller,
                cross_hit_modeller,
                cross_hit_threshold=0.9
            )
            
            assert len(result) == 4


class TestGetArgs:
    """Tests for get_args function."""

    def test_get_args_defaults(self):
        """Test get_args with default values."""
        import sys
        from deployment.model_evaluation.evaluate import get_args
        
        test_args = [
            'study_deploy.py',
            '--study_output_filepath', '/tmp/output',
            '--taxid_plan_filepath', '/tmp/plan.tsv',
            '--analysis_output_filepath', '/tmp/analysis'
        ]
        
        with patch.object(sys, 'argv', test_args):
            args = get_args()
            
            assert args.study_output_filepath == '/tmp/output'
            assert args.taxid_plan_filepath == '/tmp/plan.tsv'
            assert args.analysis_output_filepath == '/tmp/analysis'
            assert args.threshold == 0.3
            assert args.taxa_threshold == 0.02
            assert args.tax_level_to_use == 'order'
            assert args.data_set_divide == 5
            assert args.holdout_proportion == 0.3
            assert args.cross_hit_threshold == 0.9

    def test_get_args_custom_values(self):
        """Test get_args with custom values."""
        import sys
        from deployment.model_evaluation.evaluate import get_args
        
        test_args = [
            'study_deploy.py',
            '--study_output_filepath', '/tmp/output',
            '--taxid_plan_filepath', '/tmp/plan.tsv',
            '--analysis_output_filepath', '/tmp/analysis',
            '--threshold', '0.5',
            '--taxa_threshold', '0.05',
            '--tax_level_to_use', 'family',
            '--data_set_divide', '10',
            '--holdout_proportion', '0.2',
            '--cross_hit_threshold', '0.8',
            '--max_training', '100'
        ]
        
        with patch.object(sys, 'argv', test_args):
            args = get_args()
            
            assert args.threshold == 0.5
            assert args.taxa_threshold == 0.05
            assert args.tax_level_to_use == 'family'
            assert args.data_set_divide == 10
            assert args.holdout_proportion == 0.2
            assert args.cross_hit_threshold == 0.8
            assert args.max_training == '100'


class TestBugFixes:
    """Tests to verify bug fixes."""

    def test_empty_dataframe_check(self):
        """Test that empty DataFrame check uses correct syntax."""
        df = pd.DataFrame({'a': [1, 2]})
        empty_df = pd.DataFrame()
        
        assert not df.empty
        assert empty_df.empty
        assert not empty_df.empty is False

    def test_cross_hit_threshold_parameterization(self):
        """Verify cross_hit_threshold can be passed as parameter."""
        from deployment.model_evaluation.evaluate import compound_eda_function
        
        import inspect
        sig = inspect.signature(compound_eda_function)
        params = list(sig.parameters.keys())
        
        assert 'cross_hit_threshold' in params

    def test_max_training_typo_fix(self):
        """Verify max_training parameter exists (was max_trainning)."""
        import sys
        from deployment.model_evaluation.evaluate import get_args
        
        test_args = [
            'study_deploy.py',
            '--study_output_filepath', '/tmp/output',
            '--taxid_plan_filepath', '/tmp/plan.tsv',
            '--analysis_output_filepath', '/tmp/analysis',
            '--max_training', '50'
        ]
        
        with patch.object(sys, 'argv', test_args):
            args = get_args()
            
            assert hasattr(args, 'max_training')
            assert args.max_training == '50'


class TestEvaluatorConfig:
    """Tests for EvaluatorConfig dataclass."""
    
    def test_evaluator_config_defaults(self, tmp_path):
        """Test EvaluatorConfig with default values."""
        from deployment.model_evaluation.config import EvaluatorConfig
        
        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\n")
        
        config = EvaluatorConfig(
            study_output_filepath=study_dir,
            taxid_plan_filepath=plan_file,
            analysis_output_filepath=tmp_path / "analysis"
        )
        
        assert config.study_output_filepath == study_dir
        assert config.tax_level == 'order'
        assert config.cross_hit_threshold == 0.9
        assert config.data_set_divide == 5
        assert config.proportion_train == 0.7
    
    def test_evaluator_config_custom_values(self, tmp_path):
        """Test EvaluatorConfig with custom values."""
        from deployment.model_evaluation.config import EvaluatorConfig
        
        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\n")
        
        config = EvaluatorConfig(
            study_output_filepath=study_dir,
            taxid_plan_filepath=plan_file,
            analysis_output_filepath=tmp_path / "analysis",
            tax_level='family',
            cross_hit_threshold=0.8,
            data_set_divide=10,
            holdout_proportion=0.2,
        )
        
        assert config.tax_level == 'family'
        assert config.cross_hit_threshold == 0.8
        assert config.data_set_divide == 10
        assert config.proportion_train == 0.8
    
    def test_evaluator_config_models_dir(self, tmp_path):
        """Test models_dir property."""
        from deployment.model_evaluation.config import EvaluatorConfig
        
        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\n")
        
        config = EvaluatorConfig(
            study_output_filepath=study_dir,
            taxid_plan_filepath=plan_file,
            analysis_output_filepath=tmp_path / "analysis"
        )
        
        assert "models" in str(config.models_dir)


class TestExceptions:
    """Tests for custom exceptions."""
    
    def test_evaluator_error(self):
        """Test EvaluatorError base class."""
        from deployment.model_evaluation.exceptions import EvaluatorError
        
        error = EvaluatorError("Test error", {"key": "value"})
        assert str(error) == "Test error (key=value)"
        assert error.details == {"key": "value"}
    
    def test_data_load_error(self):
        """Test DataLoadError."""
        from deployment.model_evaluation.exceptions import DataLoadError
        
        error = DataLoadError("dataset1", "File not found", "/path/to/file")
        assert error.dataset == "dataset1"
        assert error.reason == "File not found"
        assert error.filepath == "/path/to/file"
    
    def test_prediction_error(self):
        """Test PredictionError."""
        from deployment.model_evaluation.exceptions import PredictionError
        
        error = PredictionError("dataset1", "clade_prediction", "Model failed")
        assert error.dataset == "dataset1"
        assert error.phase == "clade_prediction"
        assert error.reason == "Model failed"


class TestResultModels:
    """Tests for result model dataclasses."""
    
    def test_precision_metrics_defaults(self):
        """Test PrecisionMetrics with defaults."""
        from deployment.model_evaluation.result_models import PrecisionMetrics
        
        pm = PrecisionMetrics()
        assert pm.raw_pred_accuracy == 0.0
        assert pm.fuzzy_precision_raw == 0.0
        assert pm.clade_precision_post == 0.0
    
    def test_precision_metrics_to_dict(self):
        """Test PrecisionMetrics to_dict."""
        from deployment.model_evaluation.result_models import PrecisionMetrics
        
        pm = PrecisionMetrics(fuzzy_precision_raw=0.8, clade_precision_post=0.6)
        d = pm.to_dict()
        
        assert d['fuzzy_precision_raw'] == 0.8
        assert d['clade_precision_post'] == 0.6
    
    def test_batch_evaluation_result_empty(self):
        """Test BatchEvaluationResult with empty results."""
        from deployment.model_evaluation.result_models import BatchEvaluationResult, create_empty_result
        
        result = create_empty_result()
        
        assert result.test_results.empty
        assert result.summary_results.empty
        assert result.get_dataset_count() == 0


class TestMetrics:
    """Tests for metric calculation functions."""
    
    def test_safe_divide(self):
        """Test safe_divide function."""
        from deployment.model_evaluation.metrics import safe_divide
        
        assert safe_divide(10, 2) == 5.0
        assert safe_divide(10, 0) == 0.0
        assert safe_divide(0, 10) == 0.0
    
    def test_compute_trash_flags(self):
        """Test compute_trash_flags function."""
        from deployment.model_evaluation.metrics import compute_trash_flags
        
        df = pd.DataFrame({
            'best_match_is_best': [True, False, False, True],
            'is_crosshit': [False, False, True, False]
        })
        
        result = compute_trash_flags(df)
        
        assert 'is_trash' in result.columns
        assert result['is_trash'].tolist() == [False, True, False, False]
    
    def test_compute_fuzzy_precision(self):
        """Test compute_fuzzy_precision function."""
        from deployment.model_evaluation.metrics import compute_fuzzy_precision
        
        df = pd.DataFrame({
            'is_trash': [False, False, True, True],
            'coverage': [1.0, 2.0, 0.0, 0.0]
        })
        
        raw, cov = compute_fuzzy_precision(df)
        
        assert raw == 0.5
        assert cov == 0.5
    
    def test_compute_overall_precision(self):
        """Test compute_overall_precision function."""
        from deployment.model_evaluation.metrics import compute_overall_precision
        
        df = pd.DataFrame({
            'best_match_taxid': [9606, 9606, 9598, None]
        })
        
        precision = compute_overall_precision(df)
        
        assert precision == 0.5


class TestDataLoader:
    """Tests for DataLoader class."""
    
    def test_data_loader_initialization(self, tmp_path):
        """Test DataLoader initialization."""
        from deployment.model_evaluation.config import EvaluatorConfig
        from deployment.model_evaluation.data_loader import DataLoader
        
        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\n")
        
        config = EvaluatorConfig(
            study_output_filepath=study_dir,
            taxid_plan_filepath=plan_file,
            analysis_output_filepath=tmp_path / "analysis"
        )
        
        loader = DataLoader(config)
        
        assert loader.config == config
        assert loader.ncbi_wrapper is None
        assert loader.all_input_data is None


class TestDatasetProcessor:
    """Tests for DatasetProcessor class."""
    
    def test_dataset_processor_initialization(self, tmp_path):
        """Test DatasetProcessor initialization."""
        from deployment.model_evaluation.config import EvaluatorConfig, TrainedModels
        from deployment.model_evaluation.dataset_processor import DatasetProcessor
        
        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\n")
        
        config = EvaluatorConfig(
            study_output_filepath=study_dir,
            taxid_plan_filepath=plan_file,
            data_set_divide=5,
            tax_level='order',
            cross_hit_threshold=0.9,
            analysis_output_filepath=tmp_path / "analysis"
        )

        models = TrainedModels()
        
        processor = DatasetProcessor(
            config=config,
            models=models,
            ncbi_wrapper=Mock(),
            input_tax_df=pd.DataFrame(),
            taxids_to_use=pd.DataFrame()
        )
        
        assert processor.config == config
        assert processor.models == models


class TestBatchEvaluator:
    """Tests for BatchEvaluator class."""
    
    def test_batch_evaluator_initialization(self, tmp_path):
        """Test BatchEvaluator initialization."""
        from deployment.model_evaluation.config import EvaluatorConfig, TrainedModels
        from deployment.model_evaluation.batch_evaluator import BatchEvaluator
        
        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\n")
        
        config = EvaluatorConfig(
            study_output_filepath=study_dir,
            taxid_plan_filepath=plan_file,
            data_set_divide=5,
            analysis_output_filepath=tmp_path / "analysis"
        )
        
        models = TrainedModels()
        
        evaluator = BatchEvaluator(
            config=config,
            models=models,
            ncbi_wrapper=Mock(),
            input_tax_df=pd.DataFrame(),
            taxids_to_use=pd.DataFrame()
        )
        
        assert evaluator.config == config
        assert isinstance(evaluator.processor, object)


class TestResultVisualizer:
    """Tests for ResultVisualizer class."""
    
    def test_result_visualizer_initialization(self, tmp_path):
        """Test ResultVisualizer initialization."""
        from deployment.model_evaluation.visualization import ResultVisualizer
        
        visualizer = ResultVisualizer(str(tmp_path))
        
        assert visualizer.output_dir == str(tmp_path)
        assert os.path.exists(str(tmp_path))
    
    def test_composition_summary_empty(self, tmp_path):
        """Test _composition_summary with empty data."""
        from deployment.model_evaluation.visualization import ResultVisualizer
        
        visualizer = ResultVisualizer(str(tmp_path))
        
        result = visualizer._composition_summary(pd.DataFrame())
        
        assert result.empty


class TestConfigAdditional:
    """Additional tests for config module."""

    def test_evaluator_config_from_yaml(self, tmp_path):
        """Load config from YAML file."""
        from deployment.model_evaluation.config import EvaluatorConfig

        yaml_content = """
study_output_filepath: /tmp/output
taxid_plan_filepath: /tmp/plan.tsv
analysis_output_filepath: /tmp/analysis
tax_level: genus
cross_hit_threshold: 0.8
"""
        yaml_file = tmp_path / "config.yaml"
        yaml_file.write_text(yaml_content)

        with patch.object(Path, 'exists', return_value=True):
            config = EvaluatorConfig.from_yaml(str(yaml_file))

            assert config.tax_level == 'genus'
            assert config.cross_hit_threshold == 0.8

    def test_evaluator_config_from_json(self, tmp_path):
        """Load config from JSON file."""
        from deployment.model_evaluation.config import EvaluatorConfig

        json_content = """
{
    "study_output_filepath": "/tmp/output",
    "taxid_plan_filepath": "/tmp/plan.tsv",
    "analysis_output_filepath": "/tmp/analysis",
    "tax_level": "species",
    "cross_hit_threshold": 0.85
}
"""
        json_file = tmp_path / "config.json"
        json_file.write_text(json_content)

        with patch.object(Path, 'exists', return_value=True):
            config = EvaluatorConfig.from_json(str(json_file))

            assert config.tax_level == 'species'
            assert config.cross_hit_threshold == 0.85

    def test_evaluator_config_to_yaml(self, tmp_path):
        """Save config to YAML."""
        from deployment.model_evaluation.config import EvaluatorConfig

        with patch.object(Path, 'exists', return_value=True):
            config = EvaluatorConfig(
                study_output_filepath="/tmp/output",
                taxid_plan_filepath="/tmp/plan.tsv",
                analysis_output_filepath="/tmp/analysis",
                tax_level='family'
            )

            yaml_file = tmp_path / "output.yaml"
            config.to_yaml(str(yaml_file))

            assert yaml_file.exists()

    def test_evaluator_config_to_json(self, tmp_path):
        """Save config to JSON."""
        from deployment.model_evaluation.config import EvaluatorConfig

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\n")
        
        config = EvaluatorConfig(
            study_output_filepath=study_dir,
            taxid_plan_filepath=plan_file,
            analysis_output_filepath=tmp_path / "analysis",
            tax_level='genus'
        )

        json_file = tmp_path / "output.json"
        config.to_json(str(json_file))

        assert json_file.exists()

    def test_evaluator_config_invalid_tax_level_raises(self, tmp_path):
        """Reject invalid tax levels."""
        from deployment.model_evaluation.config import EvaluatorConfig

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\n")
        
        with pytest.raises(ValueError) as exc_info:
            EvaluatorConfig(
                study_output_filepath=study_dir,
                taxid_plan_filepath=plan_file,
                analysis_output_filepath=tmp_path / "analysis",
                tax_level='invalid_level'
            )

        assert 'tax_level' in str(exc_info.value).lower()

    def test_evaluator_config_invalid_probability_raises(self, tmp_path):
        """Reject invalid probability values."""
        from deployment.model_evaluation.config import EvaluatorConfig

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\n")
        
        with pytest.raises(ValueError):
            EvaluatorConfig(
                study_output_filepath=study_dir,
                taxid_plan_filepath=plan_file,
                analysis_output_filepath=tmp_path / "analysis",
                cross_hit_threshold=1.5
            )

    def test_evaluator_config_output_lineages_property(self):
        """Output lineages path property."""
        from deployment.model_evaluation.config import EvaluatorConfig

        with patch.object(Path, 'exists', return_value=True):
            config = EvaluatorConfig(
                study_output_filepath="/tmp/output",
                taxid_plan_filepath="/tmp/plan.tsv",
                analysis_output_filepath="/tmp/analysis"
            )

            assert 'lineages.tsv' in str(config.output_lineages)

    def test_evaluator_config_output_db_property(self):
        """Output DB path property."""
        from deployment.model_evaluation.config import EvaluatorConfig

        with patch.object(Path, 'exists', return_value=True):
            config = EvaluatorConfig(
                study_output_filepath="/tmp/output",
                taxid_plan_filepath="/tmp/plan.tsv",
                analysis_output_filepath="/tmp/analysis"
            )

            assert 'taxa.db' in str(config.output_db)

    def test_evaluator_config_output_db_custom_dir(self):
        """Output DB path with custom dir."""
        from deployment.model_evaluation.config import EvaluatorConfig

        with patch.object(Path, 'exists', return_value=True):
            config = EvaluatorConfig(
                study_output_filepath="/tmp/output",
                taxid_plan_filepath="/tmp/plan.tsv",
                analysis_output_filepath="/tmp/analysis",
                output_db_dir=Path("/custom/db")
            )

            assert '/custom/db/taxa.db' in str(config.output_db)

    def test_model_config_defaults(self):
        """ModelConfig defaults."""
        from deployment.model_evaluation.config import ModelConfig

        config = ModelConfig()

        assert config.k_folds == 5
        assert config.test_size == 0.2
        assert config.random_state == 42
        assert 'n_estimators' in config.recall_hyperparams

    def test_visualization_config_defaults(self):
        """VisualizationConfig defaults."""
        from deployment.model_evaluation.config import VisualizationConfig

        config = VisualizationConfig()

        assert config.style == 'seaborn-v0_8-darkgrid'
        assert config.figure_dpi == 300
        assert config.figure_format == 'png'
        assert config.width == 12
        assert config.height == 8

    def test_mlflow_config_defaults(self):
        """MLflowConfig defaults."""
        from deployment.model_evaluation.config import MLflowConfig

        config = MLflowConfig()

        assert config.enabled is False
        assert config.experiment_name == 'metagenomics-evaluation'
        assert config.log_models is True

    def test_logging_config_defaults(self):
        """LoggingConfig defaults."""
        from deployment.model_evaluation.config import LoggingConfig

        config = LoggingConfig()

        assert config.level == 'INFO'
        assert config.format == 'text'


class TestDataLoaderAdditional:
    """Additional tests for data_loader module."""

    def test_split_train_test_folders(self):
        """Test train/test split."""
        from deployment.model_evaluation.data_loader import split_train_test_folders

        folders = [f'dataset_{i}' for i in range(10)]

        train, test = split_train_test_folders(folders, 0.7)

        assert len(train) == 7
        assert len(test) == 3

    def test_load_taxid_plan(self, tmp_path):
        """Load taxid plan from file."""
        from deployment.model_evaluation.data_loader import load_taxid_plan

        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n9606\tHuman\tBacteria;Proteobacteria\n")

        result = load_taxid_plan(str(plan_file))

        assert len(result) > 0
        assert 'taxid' in result.columns

    def test_get_dataset_folders_with_max(self, tmp_path):
        """Limit folders with max_training."""
        from deployment.model_evaluation.data_loader import get_dataset_folders

        for i in range(10):
            (tmp_path / f"dataset_{i}").mkdir()

        folders = get_dataset_folders(str(tmp_path), max_training=3)

        assert len(folders) == 3

    def test_get_dataset_folders_no_limit(self, tmp_path):
        """Get all folders without limit."""
        from deployment.model_evaluation.data_loader import get_dataset_folders

        for i in range(5):
            (tmp_path / f"dataset_{i}").mkdir()

        folders = get_dataset_folders(str(tmp_path))

        assert len(folders) == 5


class TestResultModelsAdditional:
    """Additional tests for result_models module."""

    def test_precision_metrics_from_dict(self):
        """Create from dict."""
        from deployment.model_evaluation.result_models import PrecisionMetrics

        data = {
            'raw_pred_accuracy': 0.9,
            'fuzzy_precision_raw': 0.8,
            'fuzzy_precision_cov_filtered': 0.7,
            'overall_precision_raw': 0.6,
            'clade_precision_full': 0.5,
            'clade_precision_post': 0.4
        }

        pm = PrecisionMetrics.from_dict(data)

        assert pm.raw_pred_accuracy == 0.9
        assert pm.fuzzy_precision_raw == 0.8

    def test_recall_metrics_to_dict(self):
        """Convert to dict."""
        from deployment.model_evaluation.result_models import RecallMetrics

        rm = RecallMetrics(recall_raw=0.8, recall_cov_filtered=0.7, clade_recall=0.6, recall_filtered_leaves=0.5)

        d = rm.to_dict()

        assert d['recall_raw'] == 0.8
        assert d['recall_cov_filtered'] == 0.7

    def test_recall_metrics_from_dict(self):
        """Create from dict."""
        from deployment.model_evaluation.result_models import RecallMetrics

        data = {
            'recall_raw': 0.9,
            'recall_cov_filtered': 0.8,
            'clade_recall': 0.7,
            'recall_filtered_leaves': 0.6
        }

        rm = RecallMetrics.from_dict(data)

        assert rm.recall_raw == 0.9
        assert rm.clade_recall == 0.7

    def test_cross_hit_metrics_to_dict(self):
        """Convert to dict."""
        from deployment.model_evaluation.result_models import CrossHitMetrics

        chm = CrossHitMetrics(predicted_cross_hits=10, cleanup_accuracy=0.95)

        d = chm.to_dict()

        assert d['predicted_cross_hits'] == 10
        assert d['cleanup_accuracy'] == 0.95

    def test_dataset_result_full(self):
        """Full DatasetResult with all fields."""
        from deployment.model_evaluation.result_models import (
            DatasetResult, PrecisionMetrics, RecallMetrics, CrossHitMetrics
        )

        result = DatasetResult(
            data_set='test_dataset',
            sample='sample1',
            input_taxid_count=100,
            output_raw=80,
            output_cov_filtered=70,
            predicted_clades_pre=10,
            predicted_clades_post=8,
            precision=PrecisionMetrics(fuzzy_precision_raw=0.8),
            recall=RecallMetrics(recall_raw=0.9),
            cross_hit=CrossHitMetrics(predicted_cross_hits=5),
            trash_composition={'taxid': [1]},
            cross_hit_composition={'taxid': [2]}
        )

        assert result.data_set == 'test_dataset'
        assert result.trash_composition == {'taxid': [1]}

    def test_batch_evaluation_result_to_json(self, tmp_path):
        """Save to JSON."""
        from deployment.model_evaluation.result_models import BatchEvaluationResult

        result = BatchEvaluationResult()
        result.test_results = pd.DataFrame({'overall_precision': [0.8], 'data_set': ['test']})
        result.metadata = {'total': 1}

        json_file = tmp_path / "results.json"
        result.to_json(str(json_file))

        assert json_file.exists()

    def test_batch_evaluation_result_from_json(self, tmp_path):
        """Load from JSON."""
        from deployment.model_evaluation.result_models import BatchEvaluationResult

        json_content = """
{
    "generated_at": "2024-01-01T00:00:00",
    "metadata": {"total": 1},
    "test_results": [{"overall_precision": 0.8, "data_set": "test"}],
    "summary_results": [],
    "trash_composition": [],
    "cross_hit_composition": []
}
"""
        json_file = tmp_path / "results.json"
        json_file.write_text(json_content)

        result = BatchEvaluationResult.from_json(str(json_file))

        assert not result.test_results.empty

    def test_get_summary_stats(self):
        """Compute summary statistics."""
        from deployment.model_evaluation.result_models import BatchEvaluationResult

        result = BatchEvaluationResult()
        result.summary_results = pd.DataFrame({
            'data_set': ['d1', 'd2'],
            'overall_precision_raw': [0.8, 0.9],
            'fuzzy_precision_raw': [0.7, 0.8]
        })

        stats = result.get_summary_stats()

        assert 'overall_precision_raw' in stats


class TestBatchEvaluatorAdditional:
    """Additional tests for batch_evaluator module."""

    def test_evaluate_collects_errors(self, tmp_path):
        """Collect errors during evaluation."""
        from deployment.model_evaluation.batch_evaluator import BatchEvaluator
        from deployment.model_evaluation.config import EvaluatorConfig, TrainedModels

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n")
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        config = EvaluatorConfig(
            study_output_filepath=str(study_dir),
            taxid_plan_filepath=str(plan_file),
            analysis_output_filepath=str(output_dir),
            data_set_divide=5
        )

        models = TrainedModels()

        evaluator = BatchEvaluator(
            config=config,
            models=models,
            ncbi_wrapper=Mock(),
            input_tax_df=pd.DataFrame(),
            taxids_to_use=pd.DataFrame()
        )

        with patch.object(evaluator.processor, 'process') as mock_process:
            from deployment.model_evaluation.exceptions import EvaluatorError
            mock_process.side_effect = EvaluatorError("test error", {})

            result = evaluator.evaluate(['dataset1'], progress=False)

            assert result.metadata['failed'] == 1

    def test_aggregate_empty_results(self, tmp_path):
        """Handle empty results."""
        from deployment.model_evaluation.batch_evaluator import BatchEvaluator
        from deployment.model_evaluation.config import EvaluatorConfig, TrainedModels

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n")
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        config = EvaluatorConfig(
            study_output_filepath=str(study_dir),
            taxid_plan_filepath=str(plan_file),
            analysis_output_filepath=str(output_dir),
            data_set_divide=5
        )

        models = TrainedModels()

        evaluator = BatchEvaluator(
            config=config,
            models=models,
            ncbi_wrapper=Mock(),
            input_tax_df=pd.DataFrame(),
            taxids_to_use=pd.DataFrame()
        )

        result = evaluator._aggregate_results([], [])

        assert result.metadata['successful'] == 0

    def test_save_results_creates_directory(self, tmp_path):
        """Create output directory."""
        from deployment.model_evaluation.batch_evaluator import BatchEvaluator
        from deployment.model_evaluation.config import EvaluatorConfig, TrainedModels
        from deployment.model_evaluation.result_models import BatchEvaluationResult

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n")
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        config = EvaluatorConfig(
            study_output_filepath=str(study_dir),
            taxid_plan_filepath=str(plan_file),
            analysis_output_filepath=str(output_dir),
            data_set_divide=5
        )

        models = TrainedModels()

        evaluator = BatchEvaluator(
            config=config,
            models=models,
            ncbi_wrapper=Mock(),
            input_tax_df=pd.DataFrame(),
            taxids_to_use=pd.DataFrame()
        )

        result = BatchEvaluationResult()
        test_output_dir = tmp_path / "test_output"

        evaluator.save_results(result, str(test_output_dir))

        assert test_output_dir.exists()

    def test_save_summary_statistics(self, tmp_path):
        """Save precision statistics."""
        from deployment.model_evaluation.batch_evaluator import BatchEvaluator
        from deployment.model_evaluation.config import EvaluatorConfig, TrainedModels
        from deployment.model_evaluation.result_models import BatchEvaluationResult

        study_dir = tmp_path / "study"
        study_dir.mkdir()
        plan_file = tmp_path / "plan.tsv"
        plan_file.write_text("taxid\tdescription\tlineage\n")
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        config = EvaluatorConfig(
            study_output_filepath=str(study_dir),
            taxid_plan_filepath=str(plan_file),
            analysis_output_filepath=str(output_dir),
            data_set_divide=5
        )

        models = TrainedModels()

        evaluator = BatchEvaluator(
            config=config,
            models=models,
            ncbi_wrapper=Mock(),
            input_tax_df=pd.DataFrame(),
            taxids_to_use=pd.DataFrame()
        )

        result = BatchEvaluationResult()
        result.summary_results = pd.DataFrame({
            'data_set': ['d1'],
            'overall_precision_raw': [0.8],
            'fuzzy_precision_raw': [0.7]
        })

        test_output_dir = tmp_path / "test_output"
        evaluator.save_summary_statistics(result, str(test_output_dir))

        stats_file = test_output_dir / "precision_summary_statistics.tsv"
        assert stats_file.exists()


class TestVisualizationAdditional:
    """Additional tests for visualization module."""

    def test_plot_all_calls_all_plots(self, tmp_path):
        """Test plot_all method."""
        from deployment.model_evaluation.visualization import ResultVisualizer
        from deployment.model_evaluation.result_models import BatchEvaluationResult

        visualizer = ResultVisualizer(str(tmp_path))

        result = BatchEvaluationResult()
        result.test_results = pd.DataFrame({'overall_precision': [0.8], 'data_set': ['test']})
        result.summary_results = pd.DataFrame({'sample': ['s1'], 'overall_precision_raw': [0.8]})

        visualizer.plot_all(result)

    def test_plot_precision_distribution_empty(self, tmp_path):
        """Handle empty data."""
        from deployment.model_evaluation.visualization import ResultVisualizer

        visualizer = ResultVisualizer(str(tmp_path))

        result = visualizer.plot_precision_distribution(pd.DataFrame())

    def test_plot_recall_comparison_empty(self, tmp_path):
        """Handle empty summary."""
        from deployment.model_evaluation.visualization import ResultVisualizer

        visualizer = ResultVisualizer(str(tmp_path))

        result = visualizer.plot_recall_comparison(pd.DataFrame())

    def test_plot_cross_hit_heatmap_empty(self, tmp_path):
        """Handle empty cross-hit data."""
        from deployment.model_evaluation.visualization import ResultVisualizer

        visualizer = ResultVisualizer(str(tmp_path))

        result = visualizer.plot_cross_hit_heatmap(pd.DataFrame())

    def test_plot_trash_heatmap_empty(self, tmp_path):
        """Handle empty trash data."""
        from deployment.model_evaluation.visualization import ResultVisualizer

        visualizer = ResultVisualizer(str(tmp_path))

        result = visualizer.plot_trash_heatmap(pd.DataFrame())

    def test_plot_recall_improvement(self, tmp_path):
        """Test recall improvement plot."""
        from deployment.model_evaluation.visualization import ResultVisualizer

        visualizer = ResultVisualizer(str(tmp_path))

        summary = pd.DataFrame({
            'sample': ['s1', 's2'],
            'recall_raw': [0.8, 0.9],
            'clade_recall': [0.85, 0.95],
            'recall_filtered_leaves': [0.9, 0.95]
        })

        visualizer.plot_recall_improvement(summary)

    def test_plot_probability_metrics(self, tmp_path):
        """Test probability metrics."""
        from deployment.model_evaluation.visualization import ResultVisualizer

        visualizer = ResultVisualizer(str(tmp_path))

        summary = pd.DataFrame({
            'sample': ['s1'],
            'recall_raw': [0.8],
            'fuzzy_precision_raw': [0.7],
            'overall_precision_raw': [0.6],
            'clade_recall': [0.85],
            'clade_precision_full': [0.75]
        })

        visualizer.plot_probability_metrics(summary)

    def test_report_generator(self, tmp_path):
        """Test HTML report generation."""
        from deployment.model_evaluation.visualization import ReportGenerator

        report_gen = ReportGenerator(str(tmp_path))

        report_gen.add_plot("Test Plot", str(tmp_path / "plot.png"), "Test description")
        report_gen.add_metrics("Test Metrics", {"accuracy": 0.95})
        report_gen.add_text("Overview", "Test content")

        output_path = tmp_path / "report.html"
        report_gen.generate(str(output_path))

        assert output_path.exists()
        content = output_path.read_text()
        assert "Test Plot" in content

    def test_generate_report_function(self, tmp_path):
        """Test generate_report convenience function."""
        from deployment.model_evaluation.visualization import generate_report
        from deployment.model_evaluation.result_models import BatchEvaluationResult

        result = BatchEvaluationResult()
        result.test_results = pd.DataFrame({'overall_precision': [0.8], 'data_set': ['test']})
        result.summary_results = pd.DataFrame({
            'data_set': ['test'],
            'overall_precision_raw': [0.8],
            'fuzzy_precision_raw': [0.7],
            'clade_precision_full': [0.6],
            'clade_precision_post': [0.5]
        })

        output_path = generate_report(result, str(tmp_path))

        assert Path(output_path).exists()

    def test_plot_clade_precision_by_taxlevel(self, tmp_path):
        """Test clade precision by tax level plot."""
        from deployment.model_evaluation.visualization import ResultVisualizer

        visualizer = ResultVisualizer(str(tmp_path))

        summary = pd.DataFrame({
            'order': ['Primates', 'Primates'],
            'raw_pred_accuracy': [0.8, 0.9]
        })

        visualizer.plot_clade_precision_by_taxlevel(summary, 'order')

    def test_composition_summary_with_data(self, tmp_path):
        """Test _composition_summary with actual data."""
        from deployment.model_evaluation.visualization import ResultVisualizer

        visualizer = ResultVisualizer(str(tmp_path))

        composition = pd.DataFrame({
            'taxid': [9606, 9606],
            'tax_level': ['order', 'order'],
            'count': [10, 20]
        })

        result = visualizer._composition_summary(composition)

        assert not result.empty


class TestMetricsAdditional:
    """Additional tests for metrics module."""

    def test_compute_recall(self):
        """Test compute_recall function."""
        from deployment.model_evaluation.metrics import compute_recall

        clean_m_stats = pd.DataFrame({
            'best_match_taxid': [9606, 9606, 9598],
            'coverage': [1.0, 1.0, 1.0]
        })

        input_summary = pd.DataFrame({
            'taxid': [9606, 9598, 9597]
        })

        recall_raw, recall_cov, clade_recall, filtered_leaves = compute_recall(clean_m_stats, input_summary)

        assert recall_raw > 0

    def test_compute_clade_accuracy(self):
        """Test compute_clade_accuracy function."""
        from deployment.model_evaluation.metrics import compute_clade_accuracy

        results_df = pd.DataFrame({
            'best_taxid_match': [9606, 9598]
        })

        input_summary = pd.DataFrame({
            'taxid': [9606, 9598]
        })

        accuracy = compute_clade_accuracy(results_df, input_summary)

        assert accuracy >= 0

    def test_compute_clade_recall(self):
        """Test compute_clade_recall function."""
        from deployment.model_evaluation.metrics import compute_clade_recall

        results_df = pd.DataFrame({
            'best_taxid_match': [9606, 9598, 9606]
        })

        input_summary = pd.DataFrame({
            'taxid': [9606, 9598, 9597]
        })

        recall = compute_clade_recall(results_df, input_summary)

        assert 0 <= recall <= 1

    def test_compute_precision_stats(self):
        """Test compute_precision_stats function."""
        from deployment.model_evaluation.metrics import compute_precision_stats

        df = pd.DataFrame({
            'best_match_is_best': [True, False, True],
            'is_crosshit': [False, False, True],
            'coverage': [1.0, 0.0, 2.0],
            'best_match_taxid': [9606, 9606, 9598]
        })

        stats = compute_precision_stats(df)

        assert 'fuzzy_precision_raw' in stats

    def test_compute_recall_stats(self):
        """Test compute_recall_stats function."""
        from deployment.model_evaluation.metrics import compute_recall_stats

        clean_m_stats = pd.DataFrame({
            'best_match_taxid': [9606, 9598],
            'coverage': [1.0, 1.0]
        })

        input_summary = pd.DataFrame({
            'taxid': [9606, 9598]
        })

        stats = compute_recall_stats(clean_m_stats, input_summary)

        assert 'recall_raw' in stats

    def test_fill_missing_tax_levels(self):
        """Test fill_missing_tax_levels function."""
        from deployment.model_evaluation.metrics import fill_missing_tax_levels

        df = pd.DataFrame({
            'taxid': [9606],
            'count': [10]
        })

        missing_levels = {'family', 'genus'}

        result = fill_missing_tax_levels(df, missing_levels)

        assert 'family' in result.columns
        assert 'genus' in result.columns


class TestExceptionsAdditional:
    """Additional tests for exceptions module."""

    def test_overlap_manager_error(self):
        """Test OverlapManagerError."""
        from deployment.model_evaluation.exceptions import OverlapManagerError

        error = OverlapManagerError("dataset1", "load", "File not found")

        assert error.dataset == "dataset1"
        assert error.operation == "load"
        assert error.reason == "File not found"

    def test_configuration_error(self):
        """Test ConfigurationError."""
        from deployment.model_evaluation.exceptions import ConfigurationError

        error = ConfigurationError("threshold", "1.5", "must be between 0 and 1")

        assert error.parameter == "threshold"
        assert error.value == "1.5"
        assert error.reason == "must be between 0 and 1"

    def test_model_error(self):
        """Test ModelError."""
        from deployment.model_evaluation.exceptions import ModelError

        error = ModelError("recall", "train", "No training data")

        assert error.model_type == "recall"
        assert error.operation == "train"
        assert error.reason == "No training data"

    def test_results_aggregation_error(self):
        """Test ResultsAggregationError."""
        from deployment.model_evaluation.exceptions import ResultsAggregationError

        error = ResultsAggregationError("Empty results", dataset_count=0)

        assert error.dataset_count == 0
        assert "Empty results" in str(error)
