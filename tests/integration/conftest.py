"""
Pytest configuration for integration tests.
"""
import os
import shutil
import tempfile
from pathlib import Path

import pytest

# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent.parent.absolute()
TEST_RUN_DIR = PROJECT_ROOT / "test_run"
TABLES_DIR = TEST_RUN_DIR / "tables"
MINIMAL_TABLES_DIR = TABLES_DIR / "minimal"
READS_DIR = TEST_RUN_DIR / "reads"


@pytest.fixture
def project_root():
    """Return project root directory."""
    return PROJECT_ROOT


@pytest.fixture
def test_run_dir():
    """Return test_run directory."""
    return TEST_RUN_DIR


@pytest.fixture
def minimal_tables_dir():
    """Return minimal test tables directory."""
    return MINIMAL_TABLES_DIR


@pytest.fixture
def reads_dir():
    """Return test reads directory."""
    return READS_DIR


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary output directory for tests."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    yield output_dir
    # Cleanup is automatic with tmp_path


@pytest.fixture
def cleanup_nextflow():
    """Cleanup Nextflow work directories after tests."""
    yield
    # Cleanup nextflow cache/work if needed
    work_dir = PROJECT_ROOT / "work"
    if work_dir.exists():
        shutil.rmtree(work_dir)
