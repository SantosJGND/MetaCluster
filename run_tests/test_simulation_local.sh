#!/bin/bash
# Test simulation workflow - local version with hardcoded paths

set -e

# Conda and environment configuration (hardcoded for local testing)
eval "$(/home/bioinf/miniforge3/bin/conda shell.bash hook)"
source /home/bioinf/Desktop/CODE/INSA/TOOLS/metagenomics-evaluation-pipeline/.venv/bin/activate

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

TABLE="test_run/tables/minimal/minimal_single.tsv"
OUTPUT_DIR="test_run/output_test_simulation"
PARAMS_FILE="test_run/params_local.json"

echo "=== Testing Simulation Workflow (Local) ==="
echo "Input table: $TABLE"
echo "Output directory: $OUTPUT_DIR"
echo "Params file: $PARAMS_FILE"

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
