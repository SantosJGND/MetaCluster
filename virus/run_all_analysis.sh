#!/usr/bin/env bash
set -euo pipefail

# -------------------------------------------------------
# Run all virus analysis scripts
# -------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

export PYTHONPATH="${PYTHONPATH:+$PYTHONPATH:}$PROJECT_DIR"
export PATH="$PROJECT_DIR/.venv/bin:$PATH"
PYTHON="$PROJECT_DIR/.venv/bin/python"

STUDY_OUTPUT="/home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/virus/output_study"
TAXID_PLAN="/home/bioinf/Desktop/INSA/Manuscripts/Clustering_Clinical_Metagenomics/data/Panel/viral_assess.tsv"
ANALYSIS_OUTPUT="/home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/virus"
RECALL_CACHE="$ANALYSIS_OUTPUT/filter_direct_xgb_family/models/cache/recall_results_cache.parquet"

ANALYSIS_SCRIPTS="$PROJECT_DIR/deployment/model_evaluation/analysis_scripts"

echo "========================================================================"
echo " Run all virus analysis scripts"
echo "  Python:   $PYTHON"
echo "  Study:    $STUDY_OUTPUT"
echo "  Taxid:    $TAXID_PLAN"
echo "  Output:   $ANALYSIS_OUTPUT"
echo "  Cache:    $RECALL_CACHE"
echo "========================================================================"

# Input taxon composition prediction (continuous + binned)
echo ""
echo ">>> [1/1] Input taxon composition prediction (continuous + binned)"
$PYTHON "$ANALYSIS_SCRIPTS/input_composition_prediction.py" \
    --study_output_filepath "$STUDY_OUTPUT" \
    --taxid_plan_filepath "$TAXID_PLAN" \
    --analysis_output_filepath "$ANALYSIS_OUTPUT" \
    --tax_level family \
    --recall_cache "$RECALL_CACHE" \
    --verbose

echo ""
echo "Done."
