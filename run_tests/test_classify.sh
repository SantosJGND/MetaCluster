#!/bin/bash
# Test classify workflow with different configurations

set -e

# Conda and environment configuration
CONDA_PATH="${CONDA_PATH:-/home/bioinf/miniforge3/bin/conda}"
eval "$($CONDA_PATH shell.bash hook)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VENV_PATH="${VENV_PATH:-$PROJECT_ROOT/.venv}"
source "$VENV_PATH/bin/activate"

cd "$PROJECT_ROOT"

READS_DIR="test_run/reads/fastq"
ANALYSIS_ID="test_classify"
PARAMS_FILE="${PARAMS_FILE:-test_run/params_test.json}"

echo "=== Testing Classify Workflow Variants ==="
echo "Reads directory: $READS_DIR"
echo "Params file: $PARAMS_FILE"

# Validate required parameters
if ! ./test_run/check_params.sh classification "$PARAMS_FILE"; then
    echo "ERROR: Parameter validation failed"
    exit 1
fi

# Test configurations: (name, qc, centrifuge, kraken2, diamond, krakenunique)
declare -a TESTS=(
    "default:true:true:true:false:false"
    "no_qc:false:true:true:false:false"
    "kraken2_only:false:false:true:false:false"
)

run_test() {
    local name=$1
    local qc=$2
    local centrifuge=$3
    local kraken2=$4
    local diamond=$5
    local krakenunique=$6
    
    local output_dir="test_run/output_classify_${name}"
    
    echo ""
    echo "=== Test: $name ==="
    echo "QC: $qc, Centrifuge: $centrifuge, Kraken2: $kraken2, Diamond: $diamond, KrakenUnique: $krakenunique"
    
    # Clean output directory if exists
    if [ -d "$output_dir" ]; then
        rm -rf "$output_dir"
    fi
    
    # Run classify
    nextflow run deployment/classify/classify.nf \
        -profile conda \
        -params-file "$PARAMS_FILE" \
        --reads "$READS_DIR" \
        --output_dir "$output_dir" \
        --analysis_id "$ANALYSIS_ID" \
        --qc "$qc" \
        --centrifuge "$centrifuge" \
        --kraken2 "$kraken2" \
        --diamond "$diamond" \
        --krakenunique "$krakenunique"
    
    # Verify output exists
    if [ -d "$output_dir" ]; then
        echo "OK: Output directory created: $output_dir"
    else
        echo "FAIL: Output directory not created: $output_dir"
        return 1
    fi
    
    return 0
}

FAILED=0
for test in "${TESTS[@]}"; do
    IFS=':' read -r name qc centrifuge kraken2 diamond krakenunique <<< "$test"
    if ! run_test "$name" "$qc" "$centrifuge" "$kraken2" "$diamond" "$krakenunique"; then
        FAILED=1
        echo "FAIL: Test '$name' failed"
    fi
done

if [ $FAILED -eq 0 ]; then
    echo ""
    echo "=== PASS: All classify variant tests ==="
else
    echo ""
    echo "=== FAIL: Some classify variant tests failed ==="
    exit 1
fi
