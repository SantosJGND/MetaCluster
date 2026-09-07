#!/bin/bash

# Step 2 of 2: classify reads, merge results, map to references, cluster.
# Run AFTER deploy_simulate.sh. Reads are consumed from $OUTPUT_DIR/<table>/fastq.
#
# Assembly validation is fail-hard: classify.nf fails the dataset when any
# classified taxid (uniq_reads >= min_uniq_reads) has no matched assembly.
# Failed datasets are recorded in $OUTPUT_DIR/incomplete_datasets.tsv and are
# EXCLUDED from the final analysis (they are not counted as completed runs).
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

INCOMPLETE_FILE="$OUTPUT_DIR/incomplete_datasets.tsv"

if [ ! -f "$INCOMPLETE_FILE" ]; then
    mkdir -p "$OUTPUT_DIR"
    printf 'dataset\tstage\treason\ttimestamp\n' > "$INCOMPLETE_FILE"
fi

if [ -f ./check_params.sh ]; then
    ./check_params.sh full_pipeline "$PARAMS_FILE"
fi

echo "=== Deploy Study: CLASSIFY + MERGE + MAP STEP ==="
echo "Tables dir : $TABLES_DIR"
echo "Params file: $PARAMS_FILE"
echo "Output dir : $OUTPUT_DIR"
echo "Incomplete : $INCOMPLETE_FILE"

FAILED=0

for TABLE in $(ls "$TABLES_DIR"); do
    TABLE_FILE="$TABLES_DIR/$TABLE"
    [ -f "$TABLE_FILE" ] || continue
    
    ANALYSIS_ID=$(basename "$TABLE" .tsv)
    SUBSIM="$OUTPUT_DIR/${TABLE%.*}"
    READS_DIR="$SUBSIM/fastq"
    
    if [ -d "$SUBSIM/output" ]; then
        echo "Output files already exist for $ANALYSIS_ID. Skipping."
        continue
    fi
    if [ ! -d "$READS_DIR" ] || ! ls "$READS_DIR"/*.gz > /dev/null 2>&1; then
        echo "Reads not found for $ANALYSIS_ID ($READS_DIR). Run deploy_simulate.sh first."
        continue
    fi
    
    if pgrep -f "$ANALYSIS_ID" > /dev/null; then
        echo "Process is currently running for $ANALYSIS_ID. Skipping."
        continue
    fi
    
    echo "=== [classify+map] $ANALYSIS_ID ==="
    nextflow run "$DEPLOYMENT_HOME/classify/classify.nf" \
    -profile conda \
    -params-file "$PARAMS_FILE" \
    --reads "$READS_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --analysis_id "$ANALYSIS_ID" \
    -ansi-log false
    
    if [ $? -ne 0 ]; then
        echo "INCOMPLETE: classify/map failed for $ANALYSIS_ID (likely assembly validation or runtime error)"
        printf '%s\tclassify\tassembly_validation_or_runtime_error\t%s\n' \
        "$ANALYSIS_ID" "$(date -Iseconds)" >> "$INCOMPLETE_FILE"
        FAILED=1
    else
        echo "OK: classify/map finished for $ANALYSIS_ID"
    fi
    
    rm -rf work
    rm -rf .nextflow*
done

echo "=== Classify/map step finished ==="
if [ $FAILED -ne 0 ]; then
    echo "Incomplete datasets were recorded in $INCOMPLETE_FILE (excluded from final analysis)."
fi
exit $FAILED