#!/bin/bash
# Test simulation workflow

set -e

# Conda and environment configuration
CONDA_PATH="${CONDA_PATH:-/home/bioinf/miniforge3/bin/conda}"
eval "$($CONDA_PATH shell.bash hook)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VENV_PATH="${VENV_PATH:-$PROJECT_ROOT/.venv}"
source "$VENV_PATH/bin/activate"

cd "$PROJECT_ROOT"

TABLE="test_run/tables/minimal/minimal_single.tsv"
OUTPUT_DIR="test_run/output_test_simulation"
PARAMS_FILE="${PARAMS_FILE:-test_run/params_test.json}"

echo "=== Testing Simulation Workflow ==="
echo "Input table: $TABLE"
echo "Output directory: $OUTPUT_DIR"
echo "Params file: $PARAMS_FILE"

# Validate required parameters
if ! ./test_run/check_params.sh simulation "$PARAMS_FILE"; then
    echo "ERROR: Parameter validation failed"
    exit 1
fi

# Clean output directory if exists
if [ -d "$OUTPUT_DIR" ]; then
    echo "Cleaning existing output directory..."
    rm -rf "$OUTPUT_DIR"
fi

# Run simulation
nextflow run deployment/simulation/simulate.nf \
    -profile conda \
    -params-file "$PARAMS_FILE" \
    --input_table "$TABLE" \
    --output_dir "$OUTPUT_DIR"

# Verify output
echo "=== Verifying Output ==="
DATASET_NAME="minimal_single"
EXPECTED_FILES=(
    "$OUTPUT_DIR/$DATASET_NAME/input/$DATASET_NAME.tsv"
    "$OUTPUT_DIR/$DATASET_NAME/fastq/${DATASET_NAME}_R1.fq.gz"
    "$OUTPUT_DIR/$DATASET_NAME/fastq/${DATASET_NAME}_R2.fq.gz"
)

FAILED=0
for f in "${EXPECTED_FILES[@]}"; do
    if [ -f "$f" ]; then
        echo "OK: $f exists"
    else
        echo "FAIL: Missing $f"
        FAILED=1
    fi
done

if [ $FAILED -eq 0 ]; then
    echo "=== PASS: Simulation test ==="
else
    echo "=== FAIL: Simulation test ==="
    exit 1
fi
