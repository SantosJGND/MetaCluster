#!/bin/bash

# Step 1 of 2: simulate reads for each input table (independent and resumable).
# Run this BEFORE deploy_classify_map.sh. Outputs land in $OUTPUT_DIR/<table>/fastq.
#
# All remote paths are configurable via environment variables so the same
# script can target the INSAnode defaults or a local checkout.

set -u

eval "$(${CONDA_PATH:-/home/insaflu/miniforge3}/bin/conda shell.bash hook)"
source "${VENV_PATH:-/home/insaflu/work/DDI/Projects/MetaCluster/.venv}/bin/activate"
export PYTHONPATH="${PROJECT_PYTHONPATH:-/home/insaflu/work/DDI/Projects/MetaCluster}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT=$(pwd)

TABLES_DIR="${TABLES_DIR:-/home/insaflu/work/DDI/Projects/manuscript_cluster_analysis/viruses/viruses_ebinger}"
PARAMS_FILE="${PARAMS_FILE:-params_simulate.json}"
NEXTFLOW_CONFIG="${NEXTFLOW_CONFIG:-nextflow.config}"
OUTPUT_DIR="${OUTPUT_DIR:-output_study}"
DEPLOYMENT_HOME="${DEPLOYMENT_HOME:-/home/insaflu/work/DDI/Projects/MetaCluster/deployment}"

if [ -f ./check_params.sh ]; then
    ./check_params.sh full_pipeline "$PARAMS_FILE"
fi

echo "=== Deploy Study: SIMULATION STEP ==="
echo "Tables dir : $TABLES_DIR"
echo "Params file: $PARAMS_FILE"
echo "Output dir : $OUTPUT_DIR"

FAILED=0

for TABLE in $(ls "$TABLES_DIR"); do
    TABLE_FILE="$TABLES_DIR/$TABLE"
    [ -f "$TABLE_FILE" ] || continue

    ANALYSIS_ID=$(basename "$TABLE" .tsv)
    SUBSIM="$OUTPUT_DIR/${TABLE%.*}"

    if [ -d "$SUBSIM/fastq" ] && ls "$SUBSIM/fastq"/*.gz > /dev/null 2>&1; then
        echo "Simulation outputs already exist for $TABLE. Skipping."
        continue
    fi

    if pgrep -f "$ANALYSIS_ID" > /dev/null; then
        echo "Process is currently running for $ANALYSIS_ID. Skipping."
        continue
    fi

    echo "=== [simulate] $ANALYSIS_ID ==="
    nextflow run "$DEPLOYMENT_HOME/simulation/simulate.nf" \
        -profile conda \
        -params-file "$PARAMS_FILE" \
        --input_table "$TABLE_FILE" \
        --output_dir "$OUTPUT_DIR" \
        --analysis_id "$ANALYSIS_ID" \
        -ansi-log false

    if [ $? -ne 0 ]; then
        echo "WARNING: simulation failed for $ANALYSIS_ID"
        FAILED=1
    else
        echo "OK: simulation finished for $ANALYSIS_ID"
    fi

    rm -rf work
    rm -rf .nextflow*
done

if [ $FAILED -ne 0 ]; then
    echo "=== Simulation step finished with failures (see WARNING lines above) ==="
fi
exit $FAILED