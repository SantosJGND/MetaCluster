#!/bin/bash

# Reference-truth mapping step: map each dataset's simulated reads against the
# references used to simulate it (no classification / clustering), using the
# same mapping parameters as the production pipeline.
#
# References are recovered from the published study tables
# ($OUTPUT_DIR/<dataset>/input/<dataset>.tsv, produced by simulate.nf), so the
# mapped references match the ones wgsim simulated from. Statistics are written
# inside each dataset at $OUTPUT_DIR/<dataset>/mapping/ and rolled up into
# $OUTPUT_DIR/mapping_statistics_all.tsv.
#
# All remote paths are configurable via environment variables so the same
# script can target the INSAnode defaults or a local checkout.

set -u

eval "$(${CONDA_PATH:-/home/insaflu/miniforge3}/bin/conda shell.bash hook)"
source "${VENV_PATH:-/home/insaflu/work/DDI/Projects/MetaCluster/.venv}/bin/activate"
export PYTHONPATH="${PROJECT_PYTHONPATH:-/home/insaflu/work/DDI/Projects/MetaCluster}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

TABLES_DIR="${TABLES_DIR:-/home/insaflu/work/DDI/Projects/manuscript_cluster_analysis/viruses/viruses_ebinger}"
PARAMS_FILE="${PARAMS_FILE:-deployment/mapping/params.json}"
WORKFLOW="${WORKFLOW:-deployment/mapping/map_references.nf}"
OUTPUT_DIR="${OUTPUT_DIR:-output_study}"
AGGREGATE_SCRIPT="${AGGREGATE_SCRIPT:-deployment/mapping/aggregate_mapping_statistics.py}"

INCOMPLETE_FILE="$OUTPUT_DIR/incomplete_mapping.tsv"

if [ ! -f "$INCOMPLETE_FILE" ]; then
    mkdir -p "$OUTPUT_DIR"
    printf 'dataset\tstage\treason\ttimestamp\n' > "$INCOMPLETE_FILE"
fi

if [ -f ./check_params.sh ]; then
    ./check_params.sh mapping "$PARAMS_FILE"
fi

echo "=== Deploy Study: REFERENCE-TRUTH MAPPING STEP ==="
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
    INPUT_TABLE="$SUBSIM/input/$ANALYSIS_ID.tsv"
    MAPPING_DIR="$SUBSIM/mapping"

    # Prefer the published study table (references match the simulation); fall
    # back to the raw study table if simulation did not publish one.
    if [ ! -f "$INPUT_TABLE" ]; then
        if [ -f "$TABLE_FILE" ]; then
            INPUT_TABLE="$TABLE_FILE"
        else
            echo "Input table not found for $ANALYSIS_ID. Run deploy_simulate.sh first."
            continue
        fi
    fi

    if [ -f "$MAPPING_DIR/mapping_statistics.tsv" ]; then
        echo "Mapping outputs already exist for $ANALYSIS_ID. Skipping."
        continue
    fi
    if [ ! -d "$READS_DIR" ] || ! ls "$READS_DIR"/*.gz > /dev/null 2>&1; then
        echo "Reads not found for $ANALYSIS_ID ($READS_DIR). Run deploy_simulate.sh first."
        printf '%s\tmapping\tmissing_reads\t%s\n' "$ANALYSIS_ID" "$(date -Iseconds)" >> "$INCOMPLETE_FILE"
        continue
    fi
    if pgrep -f "$ANALYSIS_ID" > /dev/null; then
        echo "Process is currently running for $ANALYSIS_ID. Skipping."
        continue
    fi

    echo "=== [mapping] $ANALYSIS_ID ==="
    nextflow run "$WORKFLOW" \
        -profile conda \
        -params-file "$PARAMS_FILE" \
        --input_table "$INPUT_TABLE" \
        --reads "$READS_DIR" \
        --output_dir "$OUTPUT_DIR" \
        --analysis_id "$ANALYSIS_ID" \
        -ansi-log false

    RET=$?
    if [ $RET -ne 0 ]; then
        echo "INCOMPLETE: reference-truth mapping failed for $ANALYSIS_ID"
        printf '%s\tmapping\truntime_error\t%s\n' "$ANALYSIS_ID" "$(date -Iseconds)" >> "$INCOMPLETE_FILE"
        FAILED=1
    else
        echo "OK: reference-truth mapping finished for $ANALYSIS_ID"
    fi

    rm -rf work
    rm -rf .nextflow*
done

if [ -f "$AGGREGATE_SCRIPT" ]; then
    echo "=== Aggregating mapping statistics ==="
    python3 "$AGGREGATE_SCRIPT" --study-output "$OUTPUT_DIR" || true
fi

echo "=== Reference-truth mapping step finished ==="
if [ $FAILED -ne 0 ]; then
    echo "Incomplete datasets were recorded in $INCOMPLETE_FILE."
fi
exit $FAILED
