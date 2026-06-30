#!/bin/bash


eval "$(/home/bioinf/miniforge3/bin/conda shell.bash hook)"
source /home/bioinf/Desktop/CODE/INSA/TOOLS/metagenomics-evaluation-pipeline/.venv/bin/activate

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

TABLES_DIR="test_run/tables/study_datasets"
PARAMS_FILE="${PARAMS_FILE:-test_run/params_local.json}"
NEXTFLOW_CONFIG="test_run/nextflow.config"
OUTPUT_DIR="output_study"
DEPLOYMENT_HOME="$PROJECT_ROOT/deployment"
echo "=== Deploy Study ==="
echo "Params file: $PARAMS_FILE"

# Validate required parameters
./test_run/check_params.sh full_pipeline "$PARAMS_FILE"

for TABLE in $(ls $TABLES_DIR); do
    OUTPUT_SUBDIR="${TABLE%.*}"
    OUTPUT_SUBDIR="${OUTPUT_DIR}/${OUTPUT_SUBDIR}"
    READS_DIR="$OUTPUT_SUBDIR/fastq"
    if [ -d "$OUTPUT_SUBDIR" ]; then
        echo "Output directory $OUTPUT_SUBDIR already exists."
        
        if [ -d "$OUTPUT_SUBDIR/output" ]; then
            echo "Output files already exist for $TABLE. Skipping."
            continue 2  # Skip to next iteration of the loop
        else
            echo "Output directory $OUTPUT_SUBDIR is empty."
        fi
        
    fi
    
    echo $TABLE
    ANALYSIS_ID=$(basename "$TABLE" .tsv)
    TABLE_FILE="$TABLES_DIR/$TABLE"
    
    # check if currently any process with table name
    if pgrep -f "$TABLE" > /dev/null; then
        echo "Process is currently running for $TABLE. Skipping."
        continue 2  # Skip to next iteration of the loop
    else
        echo "No process is currently running for $TABLE. Proceeding with deployment."
    fi
    
    if pgrep -f "$ANALYSIS_ID" > /dev/null; then
        echo "Process is currently running for $ANALYSIS_ID. Skipping."
        continue 2  # Skip to next iteration of the loop
    else
        echo "No process is currently running for $ANALYSIS_ID. Proceeding with deployment."
    fi
    
    echo "Deploying $table"
    nextflow run $DEPLOYMENT_HOME/simulation/simulate.nf \
    -profile conda \
    -params-file "$PARAMS_FILE" \
    --input_table "$TABLE_FILE" \
    --output_dir "$OUTPUT_DIR" \
    --analysis_id "$ANALYSIS_ID" \
    -ansi-log false
    
    nextflow run $DEPLOYMENT_HOME/classify/classify.nf \
    -profile conda \
    -params-file "$PARAMS_FILE" \
    --reads "$READS_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --analysis_id "$ANALYSIS_ID" \
    -ansi-log false
    
    if [ $? -ne 0 ]; then
        echo "Nextflow run failed for $table"
    fi
    rm -rf work
    rm -rf .nextflow*
done
