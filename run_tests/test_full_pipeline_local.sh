#!/bin/bash
# Test complete pipeline: simulate → classify → evaluate - local version with hardcoded paths

set -e

# Conda and environment configuration (hardcoded for local testing)
eval "$(/home/bioinf/miniforge3/bin/conda shell.bash hook)"
source /home/bioinf/Desktop/CODE/INSA/TOOLS/metagenomics-evaluation-pipeline/.venv/bin/activate

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

TABLE="test_run/tables/minimal/minimal_single.tsv"
DATASET_NAME="minimal_single"
SIM_OUTPUT="test_run/output_pipeline_sim"
CLASSIFY_OUTPUT="test_run/output_pipeline_classify"
EVAL_OUTPUT="test_run/output_pipeline_eval"
PARAMS_FILE="test_run/params_local.json"

echo "=== Testing Full Pipeline (Local) ==="
echo "Input table: $TABLE"
echo "Params file: $PARAMS_FILE"

# Step 1: Simulate
echo ""
echo "=== Step 1: Simulation ==="
if [ -d "$SIM_OUTPUT" ]; then
    rm -rf "$SIM_OUTPUT"
fi

nextflow run deployment/simulation/simulate.nf \
    -profile conda \
    -params-file "$PARAMS_FILE" \
    --input_table "$TABLE" \
    --output_dir "$SIM_OUTPUT"

# Verify simulation output
if [ ! -d "$SIM_OUTPUT/$DATASET_NAME" ]; then
    echo "FAIL: Simulation output not found"
    exit 1
fi
echo "OK: Simulation complete"

# Step 2: Classify
echo ""
echo "=== Step 2: Classification ==="
SIM_READS="$SIM_OUTPUT/$DATASET_NAME/fastq"

if [ -d "$CLASSIFY_OUTPUT" ]; then
    rm -rf "$CLASSIFY_OUTPUT"
fi

nextflow run deployment/classify/classify.nf \
    -profile conda \
    -params-file "$PARAMS_FILE" \
    --reads "$SIM_READS" \
    --output_dir "$CLASSIFY_OUTPUT" \
    --analysis_id "pipeline_test"

# Verify classify output
if [ ! -d "$CLASSIFY_OUTPUT" ]; then
    echo "FAIL: Classification output not found"
    exit 1
fi
echo "OK: Classification complete"

# Step 3: Evaluate
echo ""
echo "=== Step 3: Evaluation ==="

python deployment/model_evaluation/evaluate.py \
    --study_output_filepath "$CLASSIFY_OUTPUT" \
    --taxid_plan_filepath "$TABLE" \
    --analysis_output_filepath "$EVAL_OUTPUT"

# Verify evaluate output
if [ ! -d "$EVAL_OUTPUT" ]; then
    echo "FAIL: Evaluation output not found"
    exit 1
fi
echo "OK: Evaluation complete"

echo ""
echo "=== PASS: Full pipeline test ==="
