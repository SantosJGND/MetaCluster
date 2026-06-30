"""
Test fixtures for model_evaluation tests.
"""

import pytest
import sys
from pathlib import Path

FIXTURES_DIR = Path(__file__).parent / "fixtures"
PROJECT_ROOT = Path(__file__).parent.parent.parent.parent


@pytest.fixture
def fixtures_dir():
    """Return the fixtures directory path."""
    return FIXTURES_DIR


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary output directory."""
    out = tmp_path / "output"
    out.mkdir()
    return out


@pytest.fixture
def minimal_study_dir(fixtures_dir):
    """Return the minimal study directory path."""
    return fixtures_dir / "minimal_study"


@pytest.fixture
def taxid_plan_file(fixtures_dir):
    """Return the taxid plan file path."""
    return fixtures_dir / "minimal_taxid_plan.tsv"


@pytest.fixture
def sample_m_stats():
    """Create a sample m_stats DataFrame."""
    import pandas as pd
    return pd.DataFrame({
        'best_match_taxid': [9606, 9606, 9598, None],
        'best_match_is_best': [True, True, False, False],
        'is_crosshit': [False, False, True, False],
        'coverage': [1.0, 2.0, 0.0, 0.0],
    })


@pytest.fixture
def sample_input_summary():
    """Create a sample input summary DataFrame."""
    import pandas as pd
    return pd.DataFrame({
        'sample': ['sample1', 'sample1', 'sample2'],
        'taxid': [9606, 9598, 9606],
        'reads': [100, 50, 120],
        'mutation_rate': [0.01, 0.02, 0.01],
        'order': ['Primates', 'Primates', 'Primates'],
        'family': ['Hominidae', 'Hominidae', 'Hominidae'],
    })


@pytest.fixture
def mock_ncbi_wrapper():
    """Create a mock NCBI wrapper."""
    from unittest.mock import Mock
    wrapper = Mock()
    wrapper.resolve_lineages = Mock()
    wrapper.get_level = Mock(return_value="Primates")
    wrapper.lineages = {}
    return wrapper


@pytest.fixture
def mock_models():
    """Create mock trained models."""
    from unittest.mock import Mock
    models = Mock()
    models.recall_modeller = Mock()
    models.composition_modeller = Mock()
    models.crosshit_modeller = Mock()
    models.is_complete = True
    return models


@pytest.fixture
def project_root():
    """Return the project root directory."""
    return PROJECT_ROOT
