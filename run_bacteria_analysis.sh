#!/usr/bin/env bash
set -e

ROOT=/home/bioinf/Desktop/CODE/INSA/TOOLS/metagenomics-evaluation-pipeline
STUDY=/home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/bacteria/output_study
TAXID_PLAN=/home/bioinf/Desktop/INSA/Manuscripts/Clustering_Clinical_Metagenomics/data/Panel/bacterial_assess.tsv
NCBI_DB=$STUDY/taxa.db
CACHE_DIR=/home/bioinf/Desktop/INSA/Manuscripts/Clustering_Clinical_Metagenomics/manuscript_package/data/results/bacteria/full_analysis/filter_direct_xgb_family/models/cache
RESULTS_DIR=/home/bioinf/Desktop/INSA/Manuscripts/Clustering_Clinical_Metagenomics/manuscript_package/data/results/bacteria

cd "$ROOT"
source .venv/bin/activate
export PYTHONPATH=$(pwd)

# Verify taxa.db exists in output_study
if [ ! -f "$NCBI_DB" ]; then
    echo "ERROR: taxa.db not found at $NCBI_DB"
    exit 1
fi
echo "taxa.db found at $NCBI_DB ($(du -h "$NCBI_DB" | cut -f1))"

# Clean previous outputs for fresh run
rm -rf "$RESULTS_DIR"/compare_sort_strategies/sort_comparison
rm -rf "$RESULTS_DIR"/bacteria_eda/*
rm -rf "$RESULTS_DIR"/composition_comparison_outputs/*
rm -rf "$RESULTS_DIR"/last_tp_division_outputs/*

mkdir -p "$RESULTS_DIR"/{compare_sort_strategies,bacteria_eda,composition_comparison_outputs,last_tp_division_outputs}

echo "[1/6] Sort comparison (tax_level=family)..."
python3 deployment/model_evaluation/analysis_scripts/compare_sort_strategies.py \
  --study_output_filepath "$STUDY" \
  --taxid_plan_filepath "$TAXID_PLAN" \
  --analysis_output_filepath "$RESULTS_DIR/compare_sort_strategies" \
  --tax_level family &
PID1=$!

echo "[2/6] EDA + explanatory..."
python3 deployment/model_evaluation/analysis_data_extractor.py \
  --study-output "$STUDY" \
  --ncbi-db "$NCBI_DB" \
  --output-dir "$RESULTS_DIR/bacteria_eda" \
  --explanatory &
PID2=$!

echo "[3/6] Composition model comparison..."
python3 deployment/model_evaluation/analysis_scripts/composition_model_comparison.py \
  --output_dir "$RESULTS_DIR/composition_comparison_outputs" \
  --input_cache "$CACHE_DIR/training_results_cache.parquet" &
PID3=$!

echo "[4/6] Last TP division prediction..."
python3 deployment/model_evaluation/analysis_scripts/last_tp_division_prediction_second.py \
  --output_dir "$RESULTS_DIR/last_tp_division_outputs" \
  --input_cache "$CACHE_DIR/recall_results_cache.parquet" &
PID4=$!

echo "PIDs: $PID1 $PID2 $PID3 $PID4"
echo "Waiting for steps 1-4 to finish..."
wait
echo "Steps 1-4 complete."

echo "[5/6] Input composition prediction (run 1)..."
python3 -m deployment.model_evaluation.analysis_scripts.input_composition_prediction \
  --study_output_filepath "$STUDY" \
  --taxid_plan_filepath "$TAXID_PLAN" \
  --analysis_output_filepath /home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/bacteria/composition_prediction \
  --tax_level family \
  --recall_cache /home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/bacteria/filter_direct_xgb_family/models/cache/recall_results_cache.parquet \
  --verbose

echo "[6/6] Input composition prediction (run 2)..."
python3 -m deployment.model_evaluation.analysis_scripts.input_composition_prediction \
  --study_output_filepath "$STUDY" \
  --taxid_plan_filepath "$TAXID_PLAN" \
  --analysis_output_filepath /home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/bacteria/composition_prediction \
  --tax_level family \
  --recall_cache /home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/bacteria/filter_direct_xgb_family/models/cache/recall_results_cache.parquet \
  --verbose

echo "All done."
