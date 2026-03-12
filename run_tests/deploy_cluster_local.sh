#!/bin/bash
# Deploy cluster - local version with hardcoded paths

# Conda and environment configuration (hardcoded for local testing)
eval "$(/home/bioinf/miniforge3/bin/conda shell.bash hook)"
source /home/bioinf/Desktop/CODE/INSA/TOOLS/metagenomics-evaluation-pipeline/.venv/bin/activate

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

PARAMS_FILE="test_run/params_local.json"
NEXTFLOW_CONFIG="test_run/nextflow.config"
READS_DIR="/home/bioinf/Desktop/INSA/Projectos/metagenomics-evaluation-pipeline/test_run/reads/fastq"
NEXTFLOW_SCRIPT="/home/bioinf/Desktop/INSA/Projectos/metagenomics-evaluation-pipeline/deployment/deployment_map_cluster/classify_map_cluster.nf"
ANALYSIS_ID="test_analysis"

echo "=== Deploy Cluster (Local) ==="
echo "Params file: $PARAMS_FILE"

nextflow run $NEXTFLOW_SCRIPT  \
    -params-file "$PARAMS_FILE" \
    --analysis_id $ANALYSIS_ID \
    -profile conda \
    -c "$NEXTFLOW_CONFIG" \
    --reads "$READS_DIR"
if [ $? -ne 0 ]; then
    echo "Nextflow run failed for $ANALYSIS_ID"
fi
