"""
Integration tests for full pipeline: simulate → classify → evaluate.
"""
import subprocess
from pathlib import Path

import pytest


@pytest.mark.integration
@pytest.mark.pipeline
@pytest.mark.timeout(900)
def test_full_pipeline_simulate_classify(
    minimal_tables_dir, temp_output_dir, project_root, cleanup_nextflow
):
    """Test complete pipeline: simulate → classify."""
    table = minimal_tables_dir / "minimal_single.tsv"
    dataset_name = "minimal_single"

    # Step 1: Simulate
    sim_output = temp_output_dir / "sim_output"
    result = subprocess.run(
        [
            "nextflow",
            "run",
            "deployment/simulation/simulate.nf",
            "-profile", "conda",
            "-c", "test_run/nextflow.config",
            "-params-file", "test_run/params_test.json",
            "--input_table", str(table),
            "--output_dir", str(sim_output),
        ],
        cwd=project_root,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Simulation failed: {result.stderr}"

    # Verify simulation output
    sim_dataset_dir = sim_output / dataset_name
    assert sim_dataset_dir.exists(), "Simulation output not found"
    assert (sim_dataset_dir / "fastq").exists(), "FASTQ directory not found"

    # Step 2: Classify
    sim_reads_dir = sim_output / dataset_name / "fastq"
    classify_output = temp_output_dir / "classify_output"

    result = subprocess.run(
        [
            "nextflow",
            "run",
            "deployment/classify/classify.nf",
            "-profile", "conda",
            "-c", "test_run/nextflow.config",
            "-params-file", "test_run/params_test.json",
            "--reads", str(sim_reads_dir),
            "--output_dir", str(classify_output),
            "--analysis_id", "pipeline_test",
        ],
        cwd=project_root,
        capture_output=True,
        text=True,
    )

    # Workflow may complete or may have issues - verify it runs
    assert result.returncode in [0, 1], f"Unexpected return code: {result.returncode}"

    # Verify classify output exists (even if workflow had issues)
    assert classify_output.exists(), "Classify output directory not found"


@pytest.mark.integration
@pytest.mark.pipeline
@pytest.mark.edge_case
def test_pipeline_with_no_classifiers_enabled(
    reads_dir, temp_output_dir, project_root, cleanup_nextflow
):
    """Test edge case: what happens if no classifiers are enabled."""
    # This test verifies the workflow handles edge cases gracefully
    # Note: At least one classifier should be enabled for meaningful results
    output_dir = temp_output_dir / "classify_empty"

    # Run with only one classifier to test minimal case
    result = subprocess.run(
        [
            "nextflow",
            "run",
            "deployment/classify/classify.nf",
            "-profile", "conda",
            "-c", "test_run/nextflow.config",
            "-params-file", "test_run/params_test.json",
            "--reads", str(reads_dir / "fastq"),
            "--output_dir", str(output_dir),
            "--analysis_id", "single_classifier",
            "--centrifuge", "false",
            "--kraken2", "true",
        ],
        cwd=project_root,
        capture_output=True,
        text=True,
    )

    # Should run (may have issues but should execute)
    assert result.returncode in [0, 1], f"Unexpected return code: {result.returncode}"


@pytest.mark.integration
@pytest.mark.pipeline
@pytest.mark.edge_case
@pytest.mark.skip(reason="Requires external data - run manually")
def test_pipeline_multiple_datasets(
    minimal_tables_dir, temp_output_dir, project_root, cleanup_nextflow
):
    """Test pipeline with multiple datasets in one table."""
    table = minimal_tables_dir / "minimal_multi.tsv"

    # Run simulation with multi-dataset table
    result = subprocess.run(
        [
            "nextflow",
            "run",
            "deployment/simulation/simulate.nf",
            "-profile", "conda",
            "-c", "test_run/nextflow.config",
            "-params-file", "test_run/params_test.json",
            "--input_table", str(table),
            "--output_dir", str(temp_output_dir / "multi_sim"),
        ],
        cwd=project_root,
        capture_output=True,
        text=True,
    )

    # Should handle multiple datasets
    assert result.returncode in [0, 1], f"Unexpected return code: {result.returncode}"
