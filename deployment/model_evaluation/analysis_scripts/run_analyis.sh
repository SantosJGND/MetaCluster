#!/usr/bin/env bash
set -e

ROOT=/home/ddi/TOOLS/MetaCluster
STUDY=/home/ddi/STUDIES/MetaCluster/Sims/bacteria/filter_01/output_study
TAXID_PLAN=/home/ddi/STUDIES/MetaCluster/Sims/panels/bacterial_assess.tsv
NCBI_DB=$STUDY/taxa.db
CACHE_DIR=/home/ddi/STUDIES/MetaCluster/evaluation/filter_direct_xgb_order/models/cache
RESULTS_DIR=/home/ddi/STUDIES/MetaCluster/analysis

cd "$ROOT"
source .venv/bin/activate
export PYTHONPATH=$(pwd)

mkdir -p "$RESULTS_DIR"/{compare_sort_strategies,bacteria_eda,composition_comparison_outputs,last_tp_division_outputs}

echo "[1/5] Sort comparison..."
python3 deployment/model_evaluation/analysis_scripts/compare_sort_strategies.py \
--study_output_filepath "$STUDY" \
--taxid_plan_filepath "$TAXID_PLAN" \
--analysis_output_filepath "$RESULTS_DIR/compare_sort_strategies" \
--tax_level family &
PID1=$!

echo "[2/5] EDA + explanatory..."
python3 deployment/model_evaluation/analysis_data_extractor.py \
--study-output "$STUDY" \
--ncbi-db "$NCBI_DB" \
--output-dir "$RESULTS_DIR/bacteria_eda" \
--explanatory &
PID2=$!

echo "[3/5] Composition model comparison..."
python3 deployment/model_evaluation/analysis_scripts/composition_model_comparison.py \
--output_dir "$RESULTS_DIR/composition_comparison_outputs" \
--input_cache "$CACHE_DIR/training_results_cache.parquet" &
PID3=$!

echo "[4/5] Last TP division prediction..."
python3 deployment/model_evaluation/analysis_scripts/last_tp_division_prediction_second.py \
--output_dir "$RESULTS_DIR/last_tp_division_outputs" \
--analysis_output_filepath "$(dirname "$(dirname "$CACHE_DIR")")" \
--study_output_filepath "$STUDY" &
PID4=$!

echo "[5/5] input_taxon_count_prediction.py"
#/home/bioinf/Desktop/CODE/INSA/TOOLS/metagenomics-evaluation-pipeline/.venv/bin/
python3 -m deployment.model_evaluation.analysis_scripts.input_composition_prediction \
--study_output_filepath "$STUDY" \
--taxid_plan_filepath "$TAXID_PLAN" \
--analysis_output_filepath "$RESULTS_DIR/compare_sort_strategies" \
--tax_level order \
--recall_cache "$CACHE_DIR" \
--verbose &
PID5=$!


echo "PIDs: $PID1 $PID2 $PID3 $PID4 $PID5"
echo "Waiting for all to finish..."
wait
echo "All done."
