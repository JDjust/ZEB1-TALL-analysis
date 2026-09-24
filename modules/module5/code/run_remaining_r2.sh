#!/usr/bin/env bash
set -uo pipefail
PY=/opt/bioinfo/envs/omicverse-core/bin/python
RS=/opt/bioinfo/envs/r-bio-2026.08/bin/Rscript
ROOT=/data-b/liangfuhua/projects/ZEB1_TALL_analysis
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

run_one() {
  local mod="$1" py="$2" rs="$3"
  mkdir -p "$ROOT/$mod/logs" "$ROOT/$mod/figures" "$ROOT/$mod/tables" "$ROOT/$mod/processed"
  cd "$ROOT/$mod"
  echo "[$(date)] START $mod"
  if [ -n "$py" ] && [ -f "code/$py" ]; then
    "$PY" "code/$py" > "logs/${py%.py}.log" 2>&1
    echo "$mod PREP_EXIT:$?"
    tail -n 40 "logs/${py%.py}.log"
  fi
  if [ -n "$rs" ] && [ -f "code/$rs" ]; then
    "$RS" "code/$rs" > "logs/${rs%.R}.log" 2>&1
    echo "$mod FIG_EXIT:$?"
    tail -n 30 "logs/${rs%.R}.log"
  fi
  echo "[$(date)] FIG_N_$mod:$(ls figures/M*_r2_* 2>/dev/null | wc -l)"
}

run_one module4 "" "06_leftover_r2.R"
run_one module5 "01c_round2.py" "02c_analyze_r2.R"
run_one module7 "" "02c_analyze_r2.R"
run_one module6 "01c_chip_scatac.py" "02c_analyze_r2.R"
run_one module10 "01c_scrna.py" "02c_analyze_r2.R"
echo "[$(date)] ALL_REMAINING_DONE"
