#!/bin/bash


eval "$(/home/bioinf/miniforge3/bin/conda shell.bash hook)"
source /home/bioinf/Desktop/CODE/INSA/TOOLS/metagenomics-evaluation-pipeline/.venv/bin/activate

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

TABLES_DIR="tables/study_datasets"
PARAMS_FILE="${PARAMS_FILE:-test_run/params_local.json}"
NEXTFLOW_CONFIG="test_run/nextflow.config"
OUTPUT_DIR="output_study"

echo "=== Deploy Study ==="
echo "Params file: $PARAMS_FILE"

# Validate required parameters
if ! ./test_run/check_params.sh full_pipeline "$PARAMS_FILE"; then
    echo "ERROR: Parameter validation failed"
    exit 1
fi

for table in $(ls $TABLES_DIR); do
    OUTPUT_SUBDIR="${table%.*}"
    OUTPUT_SUBDIR="${OUTPUT_DIR}/${OUTPUT_SUBDIR}"
    if [ -d "$OUTPUT_SUBDIR" ]; then
        echo "Output directory $OUTPUT_SUBDIR already exists. Skipping $table"
        continue 2  # Skip to next iteration of the loo

    fi
    echo "Deploying $table"
    nextflow run ../deployment/deployment_benchmark/simulate_map_cluster_nofilter_plus.nf \
        -params-file "$PARAMS_FILE" \
        -profile conda \
        -c "$NEXTFLOW_CONFIG" \
        --input_table $TABLES_DIR/$table 
    if [ $? -ne 0 ]; then
        echo "Nextflow run failed for $table"
    fi
    #rm -rf work
    #rm -rf .nextflow*
done
