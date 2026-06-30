#!/bin/bash


eval "$(/home/insaflu/miniforge3/bin/conda shell.bash hook)"
source /home/insaflu/work/DDI/Projects/MetaCluster/.venv/bin/activate
PYTHONPATH=/home/insaflu/work/DDI/Projects/MetaCluster
export PYTHONPATH

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT=`pwd`
#cd "$PROJECT_ROOT"

TABLES_DIR="/home/insaflu/work/DDI/Projects/manuscript_cluster_analysis/viruses/viruses_ebinger"
PARAMS_FILE="params_simulate.json"
NEXTFLOW_CONFIG="nextflow.config"
OUTPUT_DIR="output_study"
DEPLOYMENT_HOME=/home/insaflu/work/DDI/Projects/MetaCluster/deployment

echo "=== Deploy Study ==="
echo "Params file: $PARAMS_FILE"
echo "Project Dir: $PROJECT_ROOT"
# Validate required parameters
#if ! ./check_params.sh full_pipeline "$PARAMS_FILE"; then
#    echo "ERROR: Parameter validation failed"
#    exit 1
#fi

./check_params.sh full_pipeline "$PARAMS_FILE"


TABLE=$1

OUTPUT_SUBDIR="${TABLE%.*}"
OUTPUT_SUBDIR="${OUTPUT_DIR}/${OUTPUT_SUBDIR}"
READS_DIR="$OUTPUT_SUBDIR/fastq"
if [ -d "$OUTPUT_SUBDIR" ]; then
    echo "Output directory $OUTPUT_SUBDIR already exists."
    
    if [ -d "$OUTPUT_SUBDIR/output" ]; then
        echo "Output files already exist for $TABLE. Skipping."
        continue 2  # Skip to next iteration of the loop
    fi
    
    # check if currently any process with table name
    if pgrep -f "$TABLE" > /dev/null; then
        echo "Process is currently running for $TABLE. Skipping."
        continue 2  # Skip to next iteration of the loop
    fi
    
fi

echo $TABLE
ANALYSIS_ID=$(basename "$TABLE" .tsv)
TABLE_FILE="$TABLES_DIR/$TABLE"

echo "Deploying $table"
nextflow run $DEPLOYMENT_HOME/simulation/simulate.nf \
-profile conda \
-params-file "$PARAMS_FILE" \
--input_table "$TABLE_FILE" \
--output_dir "$OUTPUT_DIR" \
--analysis_id "$ANALYSIS_ID" \
-w single-work \
-ansi-log false

nextflow run $DEPLOYMENT_HOME/classify/classify.nf \
-profile conda \
-params-file "$PARAMS_FILE" \
--reads "$READS_DIR" \
--output_dir "$OUTPUT_DIR" \
--analysis_id "$ANALYSIS_ID" \
-w single-work \
-ansi-log false
