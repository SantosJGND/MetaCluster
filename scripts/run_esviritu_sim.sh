#!/usr/bin/env bash
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SIM_ROOT=""
OUTPUT_DIR=""
EVAL_DIR=""
DATASETS_ARG=""
DATASETS_FILE=""
DATASET_GLOB="*"
ESVIRITU_BIN="EsViritu"
ESVIRITU_DB="${ESVIRITU_DB:-}"
CONDA_ENV=""
CONDA_BIN=""
CPUS=8
MMP="sr"
MIN_ALIGN_LENGTH=100
QUAL="False"
FILTER_SEQS="False"
MEMORY_MAX="24G"
NICE_VAL=""
KEEP_TEMP=0
FORCE=0
STOP_ON_ERROR=0
DRY_RUN=0
VERBOSE=0
MANIFEST=""
EVAL_SCRIPT="$SCRIPT_DIR/evaluate_esviritu_sim.py"
EVAL_MIN_READS=1
PLOTS=1
DATASETS=()
FAIL_REASON=""
log(){ printf '%s [%s] %s\n' "$(date '+%F %T')" "$1" "$2"; }
info(){ log INFO "$*"; }
warn(){ log WARN "$*" >&2; }
err(){ log ERROR "$*" >&2; }
die(){ err "$*"; exit 2; }
usage(){
cat <<EOF
Usage: $(basename "$0") [options]
  --sim-root PATH
  --output-dir PATH
  --eval-dir PATH
  --datasets LIST
  --datasets-file PATH
  --dataset-glob PATTERN
  --esviritu-bin PATH
  --esviritu-db PATH
  --conda-env NAME
  --conda-bin PATH
  --cpus N
  --mmp PRESET
  --min-align-length N
  --quality True|False
  --filter-seqs True|False
  --memory-max SIZE
  --nice N
  --eval-min-reads N
  --no-plots
  --keep-temp
  --force
  --stop-on-error
  --dry-run
  -v|--verbose
  -h|--help
EOF
}
parse_args(){
while (($#)); do
case "$1" in
--sim-root) SIM_ROOT="${2:?}"; shift 2;;
--output-dir) OUTPUT_DIR="${2:?}"; shift 2;;
--eval-dir) EVAL_DIR="${2:?}"; shift 2;;
--datasets) DATASETS_ARG="${2:?}"; shift 2;;
--datasets-file) DATASETS_FILE="${2:?}"; shift 2;;
--dataset-glob) DATASET_GLOB="${2:?}"; shift 2;;
--esviritu-bin) ESVIRITU_BIN="${2:?}"; shift 2;;
--esviritu-db) ESVIRITU_DB="${2:?}"; shift 2;;
--conda-env) CONDA_ENV="${2:?}"; shift 2;;
--conda-bin) CONDA_BIN="${2:?}"; shift 2;;
--cpus) CPUS="${2:?}"; shift 2;;
--mmp) MMP="${2:?}"; shift 2;;
--min-align-length) MIN_ALIGN_LENGTH="${2:?}"; shift 2;;
--quality) QUAL="${2:?}"; shift 2;;
--filter-seqs) FILTER_SEQS="${2:?}"; shift 2;;
--memory-max) MEMORY_MAX="${2:?}"; shift 2;;
--nice) NICE_VAL="${2:?}"; shift 2;;
--eval-min-reads) EVAL_MIN_READS="${2:?}"; shift 2;;
--no-plots) PLOTS=0; shift;;
--keep-temp) KEEP_TEMP=1; shift;;
--force) FORCE=1; shift;;
--stop-on-error) STOP_ON_ERROR=1; shift;;
--dry-run) DRY_RUN=1; shift;;
-v|--verbose) VERBOSE=1; shift;;
-h|--help) usage; exit 0;;
*) usage >&2; die "unknown $1";;
esac
done
[[ -n $SIM_ROOT ]] || die "--sim-root required"
SIM_ROOT="${SIM_ROOT%/}"
[[ -n $OUTPUT_DIR ]] || OUTPUT_DIR="$SIM_ROOT/esviritu_sim"
OUTPUT_DIR="${OUTPUT_DIR%/}"
[[ -n $EVAL_DIR ]] || EVAL_DIR="$OUTPUT_DIR/evaluation"
EVAL_DIR="${EVAL_DIR%/}"
MANIFEST="$OUTPUT_DIR/esviritu_sim_manifest.tsv"
[[ -n $CONDA_BIN ]] || CONDA_BIN="$(command -v conda 2>/dev/null || echo /home/insaflu/miniforge3/bin/conda)"
}
discover_datasets(){
DATASETS=()
local -a found=()
if [[ -n $DATASETS_ARG ]]; then
  IFS=',' read -r -a found <<<"$DATASETS_ARG"
elif [[ -n $DATASETS_FILE ]]; then
  [[ -f $DATASETS_FILE ]] || die "datasets file not found"
  mapfile -t found < <(grep -v '^[[:space:]]*$' "$DATASETS_FILE" | tr -d ' \r')
else
  # try fastq dirs first
  while IFS= read -r fdir; do
    [[ -z $fdir ]] && continue
    local ddir; ddir="$(dirname "$fdir")"
    local dname; dname="$(basename "$ddir")"
    [[ "$dname" == $DATASET_GLOB || $DATASET_GLOB == "*" ]] && found+=("$dname")
  done < <(find "$SIM_ROOT" -type d -name fastq 2>/dev/null | head -50000)
  if ((${#found[@]} == 0)); then
    while IFS= read -r f; do
      [[ -z $f ]] && continue
      local fdir2; fdir2="$(dirname "$f")"
      local ddir2; ddir2="$(dirname "$fdir2")"
      local dname2; dname2="$(basename "$ddir2")"
      [[ "$dname2" == $DATASET_GLOB || $DATASET_GLOB == "*" ]] && found+=("$dname2")
    done < <(find "$SIM_ROOT" -type f \( -name "*_R1.fq.gz" -o -name "*_1.fq.gz" \) 2>/dev/null | head -50000)
  fi
fi
local -A seen=()
local d
for d in "${found[@]}"; do
  [[ -n $d ]] || continue
  [[ -n ${seen[$d]:-} ]] && continue
  seen[$d]=1
  DATASETS+=("$d")
done
((${#DATASETS[@]})) || die "no datasets discovered under $SIM_ROOT"
}
READ_FILES=(); READ_FMT=""; READ_NOTE=""
resolve_sim_reads(){
local id=$1 dir="$SIM_ROOT/$id/fastq"; READ_FILES=(); READ_FMT=""; READ_NOTE=""
[[ -d $dir ]] || { READ_NOTE="no $dir"; return 1; }
local r1="" r2=""
if [[ -f "$dir/${id}_R1.fq.gz" && -f "$dir/${id}_R2.fq.gz" ]]; then r1="$dir/${id}_R1.fq.gz"; r2="$dir/${id}_R2.fq.gz"
elif [[ -f "$dir/${id}_1.fq.gz" && -f "$dir/${id}_2.fq.gz" ]]; then r1="$dir/${id}_1.fq.gz"; r2="$dir/${id}_2.fq.gz"
else
  local -a files=()
  while IFS= read -r f; do [[ -n $f ]] && files+=("$f"); done < <(find "$dir" -maxdepth 1 -type f \( -name "*R1*.fq.gz" -o -name "*_1.fq.gz" \) | sort)
  ((${#files[@]}>=1)) && r1="${files[0]}"
  if [[ -n $r1 ]]; then
    local base; base="$(basename "$r1")"
    if [[ $base == *R1* ]]; then r2="${r1/R1/R2}"; fi
    if [[ -z $r2 && $base == *_1* ]]; then r2="${r1/_1./_2.}"; fi
  fi
fi
if [[ -n $r1 && -n $r2 && -f $r1 && -f $r2 ]]; then READ_FILES=("$r1" "$r2"); READ_FMT="paired"; return 0; fi
READ_NOTE="no paired"; return 1
}
sample_done(){ local id=$1; [[ -s "$OUTPUT_DIR/$id.detected_virus.info.tsv" ]] && [[ -s "$OUTPUT_DIR/${id}_esviritu.readstats.yaml" ]]; }
create_empty_outputs(){ local id=$1; mkdir -p "$OUTPUT_DIR"; printf 'sample_ID\tName\tdescription\tLength\tSegment\tAccession\tAssembly\tAsm_length\tkingdom\tphylum\ttclass\torder\tfamily\tgenus\tspecies\tsubspecies\tadj_taxonomy\tRPKMF\tread_count\tcovered_bases\tmean_coverage\tavg_read_identity\tconsensus_ref_identity\tPi\tfiltered_reads_in_sample\n' > "$OUTPUT_DIR/$id.detected_virus.info.tsv"; [[ -s "$OUTPUT_DIR/${id}_esviritu.readstats.yaml" ]]||printf 'reads for EsViritu denominator: 0\n'> "$OUTPUT_DIR/${id}_esviritu.readstats.yaml"; }
CAP_MODE=""; CAP_DETECTED=0
detect_cap_mode(){ (( CAP_DETECTED ))&&return 0; CAP_DETECTED=1; [[ -z $MEMORY_MAX || $MEMORY_MAX == off || $MEMORY_MAX == 0 ]]&&{ CAP_MODE=disabled; return 0; }; if command -v systemd-run >/dev/null 2>&1; then systemd-run --user --scope --quiet --collect -p "MemoryMax=$MEMORY_MAX" true >/dev/null 2>&1&&{ CAP_MODE=user; return 0; }; systemd-run --scope --quiet --collect -p "MemoryMax=$MEMORY_MAX" true >/dev/null 2>&1&&{ CAP_MODE=system; return 0; }; fi; CAP_MODE=unavailable; }
run_capped(){ detect_cap_mode; local -a cmd; if [[ -n $NICE_VAL ]]&&command -v nice >/dev/null 2>&1; then cmd=(nice -n "$NICE_VAL" "$@"); else cmd=("$@"); fi; case $CAP_MODE in user) systemd-run --user --scope --quiet --collect -p "MemoryMax=$MEMORY_MAX" -p MemorySwapMax=0 -- "${cmd[@]}";; system) systemd-run --scope --quiet --collect -p "MemoryMax=$MEMORY_MAX" -p MemorySwapMax=0 -- "${cmd[@]}";; disabled) "${cmd[@]}";; *) "${cmd[@]}";; esac; }
run_in_env(){ if [[ -z $CONDA_ENV ]]; then run_capped "$@"; return $?; fi; run_capped bash -c 'eval "$("$1" shell.bash hook)"||exit 1; "$1" activate "$2"||exit 1; shift 2; exec "$@"' _ "$CONDA_BIN" "$CONDA_ENV" "$@"; }
stage_preflight(){ local id=$1; [[ -z $ESVIRITU_DB ]]&&{ err "$id no db"; return 1; }; [[ -d $ESVIRITU_DB ]]||{ err "$id db missing"; return 1; }; local m=(); [[ -s "$ESVIRITU_DB/virus_pathogen_database.fna" ]]||m+=("fna"); [[ -s "$ESVIRITU_DB/virus_pathogen_database.all_metadata.tsv" ]]||m+=("tsv"); ((${#m[@]}))&&{ err "$id db bad"; return 1; }; resolve_sim_reads "$id"||{ err "$id reads bad"; return 1; }; return 0; }
stage_run(){ local id=$1; local -a cmd=("$ESVIRITU_BIN" -r "${READ_FILES[@]}" -s "$id" -o "$OUTPUT_DIR" -p "$READ_FMT" -t "$CPUS" -mmP "$MMP" --min-align-length "$MIN_ALIGN_LENGTH" -q "$QUAL" -f "$FILTER_SEQS" --db "$ESVIRITU_DB"); (( KEEP_TEMP ))&&cmd+=("--keep"); (( DRY_RUN ))&&{ info "$id would run"; return 0; }; local logf="$OUTPUT_DIR/${id}_esviritu_sim_driver.out"; run_in_env "${cmd[@]}" > "$logf" 2>&1; local rc=$?; if (( rc!=0 )); then if sample_done "$id"; then warn "$id tables ok despite rc"; return 0; fi; if grep -q "No reads aligned" "$logf" 2>/dev/null; then create_empty_outputs "$id"; sample_done "$id"&&return 0; fi; err "$id rc $rc"; return 1; fi; return 0; }
stage_verify(){ local id=$1; (( DRY_RUN ))&&return 0; sample_done "$id"||{ err "$id missing"; return 1; }; return 0; }
stage_evaluate(){ local id=$1; (( DRY_RUN ))&&return 0; [[ -f $EVAL_SCRIPT ]]||return 0; local args=(--esviritu-dir "$OUTPUT_DIR" --sim-root "$SIM_ROOT" --output-dir "$EVAL_DIR" --datasets "$id" --min-reads "$EVAL_MIN_READS"); (( PLOTS ))&&args+=(--plots)||args+=(--no-plots); (( VERBOSE ))&&args+=(--verbose); python3 "$EVAL_SCRIPT" "${args[@]}" >/dev/null 2>&1||true; return 0; }
process_sample(){ local id=$1; FAIL_REASON=""; stage_preflight "$id"||{ FAIL_REASON="preflight"; return 1; }; stage_run "$id"||{ FAIL_REASON="run"; return 1; }; stage_verify "$id"||{ FAIL_REASON="verify"; return 1; }; stage_evaluate "$id"||true; return 0; }
manifest_init(){ [[ -f $MANIFEST ]]||{ mkdir -p "$(dirname "$MANIFEST")"; printf 'dataset\tstatus\tnote\tstarted\tfinished\n'> "$MANIFEST"; }; }
manifest_add(){ printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5">> "$MANIFEST"; }
main(){
parse_args "$@"; discover_datasets; info "sim:$SIM_ROOT out:$OUTPUT_DIR eval:$EVAL_DIR n:${#DATASETS[@]}"; (( DRY_RUN ))&&{ for d in "${DATASETS[@]}"; do info "would $d"; done; exit 0; }; mkdir -p "$OUTPUT_DIR" "$EVAL_DIR"; manifest_init; local n_ok=0 n_skip=0; local -a failures=(); for id in "${DATASETS[@]}"; do local started finished status="ok" note="-"; started=$(date '+%F %T'); if sample_done "$id" && (( ! FORCE )); then info "$id skip"; note="exists"; ((n_skip++)); else info "$id start"; process_sample "$id"; local rc=$?; if (( rc==0 )); then note="ok"; ((n_ok++)); else status="fail"; note="${FAIL_REASON:-unknown}"; failures+=("$id"); err "$id fail"; (( STOP_ON_ERROR ))&&break; fi; fi; finished=$(date '+%F %T'); manifest_add "$id" "$status" "$note" "$started" "$finished"; done; info "ok:$n_ok skip:$n_skip fail:${#failures[@]}"; [[ -f $EVAL_SCRIPT ]]&&python3 "$EVAL_SCRIPT" --esviritu-dir "$OUTPUT_DIR" --sim-root "$SIM_ROOT" --output-dir "$EVAL_DIR" --min-reads "$EVAL_MIN_READS" $([ $PLOTS -eq 1 ]&&echo --plots||echo --no-plots) >/dev/null 2>&1||true; ((${#failures[@]}))&&exit 1; exit 0; }
main "$@"
