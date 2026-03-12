#!/bin/bash
# Parameter validation script for Nextflow workflows
# Validates required parameters exist in the params file before running tests

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
    echo "Usage: $0 <simulation|classification|full_pipeline> [params_file]"
    echo ""
    echo "Arguments:"
    echo "  simulation     - Validate only simulation parameters"
    echo "  classification - Validate only classification parameters"
    echo "  full_pipeline  - Validate all parameters"
    echo "  params_file    - Optional: path to params file (default: test_run/params_test.json)"
    echo ""
    echo "Examples:"
    echo "  $0 simulation"
    echo "  $0 classification test_run/params_local.json"
    echo "  $0 full_pipeline"
}

if [ $# -lt 1 ]; then
    usage
    exit 1
fi

WORKFLOW_TYPE=$1
PARAMS_FILE="${2:-test_run/params_test.json}"

echo "=== Parameter Validation ==="
echo "Workflow: $WORKFLOW_TYPE"
echo "Params file: $PARAMS_FILE"

if [ ! -f "$PARAMS_FILE" ]; then
    echo -e "${RED}ERROR: Params file not found: $PARAMS_FILE${NC}"
    exit 1
fi

check_param() {
    local param_name=$1
    local param_value=$(python3 -c "import json; params=json.load(open('$PARAMS_FILE')); print(params.get('$param_name', ''))")
    
    if [ -z "$param_value" ]; then
        echo -e "${RED}  MISSING: $param_name${NC}"
        return 1
    else
        echo -e "${GREEN}  OK: $param_name = $param_value${NC}"
        return 0
    fi
}

check_conda_env() {
    local env_name=$1
    local env_path=$(python3 -c "import json; params=json.load(open('$PARAMS_FILE')); print(params.get('$env_name', ''))")
    
    if [ -z "$env_path" ]; then
        echo -e "${RED}  MISSING: $env_name${NC}"
        return 1
    fi
    
    if [[ "$env_path" == "path/to/"* ]]; then
        echo -e "${YELLOW}  PLACEHOLDER: $env_name (not a real path)${NC}"
        return 0
    fi
    
    if [ -d "$env_path" ]; then
        echo -e "${GREEN}  OK: $env_name = $env_path${NC}"
        return 0
    else
        echo -e "${RED}  NOT FOUND: $env_name = $env_path (directory does not exist)${NC}"
        return 1
    fi
}

validate_simulation() {
    echo ""
    echo "=== Validating Simulation Parameters ==="
    local failed=0
    
    check_param "input_table" || failed=1
    check_param "output_dir" || failed=1
    check_param "assembly_store" || failed=1
    check_param "python_bin" || failed=1
    check_param "references_extract_script" || failed=1
    check_param "wgsim_python_path" || failed=1
    check_param "wgsim_args" || failed=1
    check_conda_env "mess_conda_env" || failed=1
    check_conda_env "minimap2_conda_env" || failed=1
    
    return $failed
}

validate_classification() {
    echo ""
    echo "=== Validating Classification Parameters ==="
    local failed=0
    
    check_param "reads" || failed=1
    check_param "output_dir" || failed=1
    check_param "python_bin" || failed=1
    check_param "classifier_process_script" || failed=1
    check_param "references_extract_script" || failed=1
    check_param "prinseq_params" || failed=1
    check_conda_env "prinseq_conda_env" || failed=1
    check_conda_env "kraken2_conda_env" || failed=1
    check_conda_env "centrifuge_conda_env" || failed=1
    check_conda_env "diamond_conda_env" || failed=1
    check_conda_env "krakenunique_conda_env" || failed=1
    check_conda_env "minimap2_conda_env" || failed=1
    check_conda_env "msamtools_conda_env" || failed=1
    check_conda_env "mosdepth_conda_env" || failed=1
    check_param "minimap2_illumina_params" || failed=1
    check_param "msamtools_params" || failed=1
    check_param "kraken2_index" || failed=1
    check_param "centrifuge_index" || failed=1
    
    return $failed
}

validate_full_pipeline() {
    echo ""
    echo "=== Validating Full Pipeline Parameters ==="
    local failed=0
    
    validate_simulation || failed=1
    validate_classification || failed=1
    
    return $failed
}

case $WORKFLOW_TYPE in
    simulation)
        validate_simulation
        ;;
    classification)
        validate_classification
        ;;
    full_pipeline)
        validate_full_pipeline
        ;;
    *)
        echo -e "${RED}Unknown workflow type: $WORKFLOW_TYPE${NC}"
        usage
        exit 1
        ;;
esac

exit $?
