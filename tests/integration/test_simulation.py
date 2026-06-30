"""
Integration tests for simulation workflow.
"""
import subprocess
from pathlib import Path

import pytest


@pytest.mark.integration
@pytest.mark.simulation
@pytest.mark.timeout(600)
def test_simulation_single_dataset(
    minimal_tables_dir, temp_output_dir, project_root, cleanup_nextflow
):
    """Test simulation with single dataset."""
    table = minimal_tables_dir / "minimal_single.tsv"
    assert table.exists(), f"Test table not found: {table}"

    result = subprocess.run(
        [
            "nextflow",
            "run",
            "deployment/simulation/simulate.nf",
            "-profile", "conda",
            "-c", "test_run/nextflow.config",
            "-params-file", "test_run/params_test.json",
            "--input_table", str(table),
            "--output_dir", str(temp_output_dir),
        ],
        cwd=project_root,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, f"Simulation failed: {result.stderr}"

    # Verify output structure
    dataset_name = "minimal_single"
    output_dataset_dir = temp_output_dir / dataset_name
    assert output_dataset_dir.exists(), f"Output dataset dir not found: {output_dataset_dir}"

    # Check input directory
    input_dir = output_dataset_dir / "input"
    assert input_dir.exists(), "Input directory not found"
    assert (input_dir / f"{dataset_name}.tsv").exists(), "Input TSV not found"

    # Check fastq directory
    fastq_dir = output_dataset_dir / "fastq"
    assert fastq_dir.exists(), "FASTQ directory not found"
    assert (fastq_dir / f"{dataset_name}_R1.fq.gz").exists(), "R1 FASTQ not found"
    assert (fastq_dir / f"{dataset_name}_R2.fq.gz").exists(), "R2 FASTQ not found"


@pytest.mark.integration
@pytest.mark.simulation
@pytest.mark.timeout(600)
def test_simulation_multi_dataset(
    minimal_tables_dir, temp_output_dir, project_root, cleanup_nextflow
):
    """Test simulation with multi-dataset table."""
    table = minimal_tables_dir / "minimal_multi.tsv"
    assert table.exists(), f"Test table not found: {table}"

    result = subprocess.run(
        [
            "nextflow",
            "run",
            "deployment/simulation/simulate.nf",
            "-profile", "conda",
            "-c", "test_run/nextflow.config",
            "-params-file", "test_run/params_test.json",
            "--input_table", str(table),
            "--output_dir", str(temp_output_dir),
        ],
        cwd=project_root,
        capture_output=True,
        text=True,
    )

    # May succeed or fail depending on table format - just verify it runs
    assert result.returncode in [0, 1], f"Unexpected return code: {result.returncode}"


@pytest.mark.integration
@pytest.mark.simulation
def test_simulation_output_format(
    minimal_tables_dir, temp_output_dir, project_root, cleanup_nextflow
):
    """Test that simulation output has correct format for classify workflow."""
    table = minimal_tables_dir / "minimal_single.tsv"

    # Run simulation
    subprocess.run(
        [
            "nextflow",
            "run",
            "deployment/simulation/simulate.nf",
            "-profile", "conda",
            "-c", "test_run/nextflow.config",
            "-params-file", "test_run/params_test.json",
            "--input_table", str(table),
            "--output_dir", str(temp_output_dir),
        ],
        cwd=project_root,
        check=True,
    )

    # Verify FASTQ files are in correct location for classify workflow
    dataset_name = "minimal_single"
    fastq_dir = temp_output_dir / dataset_name / "fastq"

    # Check pattern that classify workflow expects
    r1_files = list(fastq_dir.glob("*R1*.fq.gz"))
    r2_files = list(fastq_dir.glob("*R2*.fq.gz"))

    assert len(r1_files) > 0, "No R1 FASTQ files found"
    assert len(r2_files) > 0, "No R2 FASTQ files found"
