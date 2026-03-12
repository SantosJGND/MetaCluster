#!/bin/bash

# Conda and environment configuration
CONDA_PATH="${CONDA_PATH:-/home/bioinf/miniforge3/bin/conda}"
eval "$($CONDA_PATH shell.bash hook)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VENV_PATH="${VENV_PATH:-$PROJECT_ROOT/.venv}"
source "$VENV_PATH/bin/activate"

cd "$PROJECT_ROOT"

PARAMS_FILE="${PARAMS_FILE:-test_run/params_local.json}"
NEXTFLOW_CONFIG="test_run/nextflow.config"
READS_DIR="/home/bioinf/Desktop/INSA/Projectos/metagenomics-evaluation-pipeline/test_run/reads/fastq"
NEXTFLOW_SCRIPT="/home/bioinf/Desktop/INSA/Projectos/metagenomics-evaluation-pipeline/deployment/deployment_map_cluster/classify_map_cluster.nf"
ANALYSIS_ID="test_analysis"

echo "=== Deploy Cluster ==="
echo "Params file: $PARAMS_FILE"

# Validate required parameters
if ! ./test_run/check_params.sh classification "$PARAMS_FILE"; then
    echo "ERROR: Parameter validation failed"
    exit 1
fi

nextflow run $NEXTFLOW_SCRIPT  \
    -params-file "$PARAMS_FILE" \
    --analysis_id $ANALYSIS_ID \
    -profile conda \
    -c "$NEXTFLOW_CONFIG" \
    --reads "$READS_DIR"
if [ $? -ne 0 ]; then
    echo "Nextflow run failed for $ANALYSIS_ID"
fi
