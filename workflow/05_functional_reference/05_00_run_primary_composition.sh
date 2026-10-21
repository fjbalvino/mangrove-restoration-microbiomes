#!/usr/bin/env bash
# Reproduce read-based composition, reference contrasts, cores and figures.
# Outputs remain local. Progress is reported every 30 seconds.
set +e
run_functional_main() {
  local task_dir task_repo task_out task_r task_py task_script
  task_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd) || return 1
  task_repo=$(cd -- "$task_dir/../.." && pwd) || return 1
  task_out="${FUNCTIONAL_OUT_ROOT:-$task_repo/local_runs/functional_reference}"
  task_r="${FUNCTIONAL_RSCRIPT:-Rscript}"
  task_py="${FUNCTIONAL_PYTHON:-python3}"
  if [ "${1:-}" != "--pipeline" ]; then
    command -v "$task_r" >/dev/null 2>&1 || { printf 'Set FUNCTIONAL_RSCRIPT to an installed Rscript.\n'; return 1; }
    "$task_py" "$task_repo/checks/check_minimal_repository.py" || return $?
    "$task_py" "$task_repo/checks/check_functional_inputs.py" || return $?
    "$task_r" -e 'for (x in c("vegan","permute","ggplot2","dplyr")) if (!requireNamespace(x,quietly=TRUE)) stop("Missing R package: ",x)' || return $?
    "$task_py" -c 'import numpy, pandas, scipy, matplotlib' || return $?
    AUDIT_LOG_DIR="$task_repo/local_runs/logs" bash "$task_repo/tools/run_monitored.sh" \
      bash "$task_dir/05_00_run_primary_composition.sh" --pipeline
    return $?
  fi
  "$task_r" "$task_dir/05_01_composition_reference.R" \
    --base "$task_repo" --out_root "$task_out" || return $?
  for task_script in 05_02_partial_dbRDA.R 05_03_depth_reference_contrasts.R; do
    "$task_r" "$task_dir/$task_script" --out_root "$task_out" || return $?
  done
  "$task_r" "$task_dir/05_04_reference_core.R" \
    --base "$task_repo" --out_root "$task_out" || return $?
  "$task_r" "$task_dir/05_05_validate_results.R" \
    --base "$task_repo" --out_root "$task_out" || return $?
  "$task_py" "$task_dir/05_06_assemble_figures.py" --out_root "$task_out" || return $?
  printf 'FUNCTIONAL_REFERENCE_COMPLETE | Results: %s\n' "$task_out"
}
run_functional_main "$@"
