#!/usr/bin/env bash
# Corrected input preparation, network reconstruction, sensitivity, audit and figures.
set -euo pipefail
if [[ ${1:-} == --help || $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: bash $0 MATRIX_DIRECTORY METADATA_CSV [NEW_OUTPUT_DIRECTORY]"
  echo "The two corrected KO matrices are external; an existing output directory is refused."
  [[ ${1:-} == --help ]] && exit 0 || exit 2
fi
network_scripts="$(cd "$(dirname "${BASH_SOURCE[0]}")/read_based" && pwd)"
network_matrices="$(realpath "$1")"
network_metadata="$(realpath "$2")"
network_out="$(realpath -m "${3:-local_runs/read_based_networks}")"
network_python="${PYTHON:-python3}"
if [[ -e "$network_out" ]]; then
  echo "Output already exists: $network_out. Choose a fresh directory." >&2
  exit 2
fi
mkdir -p "$network_out/logs"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
network_pid=""
trap 'if [[ -n "$network_pid" ]]; then kill "$network_pid" 2>/dev/null || true; fi' INT TERM EXIT
run_step() {
  local step="$1"; shift
  local log="$network_out/logs/$step.log"
  echo "$(date -u +%FT%TZ) | $step started"
  "$@" > "$log" 2>&1 &
  network_pid=$!
  while kill -0 "$network_pid" 2>/dev/null; do
    sleep 30
    if kill -0 "$network_pid" 2>/dev/null; then
      echo "$(date -u +%FT%TZ) | $step running"
      tail -n 2 "$log"
    fi
  done
  if wait "$network_pid"; then
    network_pid=""
    tail -n 3 "$log"
  else
    network_pid=""
    tail -n 30 "$log" >&2
    return 1
  fi
}
run_step prepare "$network_python" -u "$network_scripts/prepare_inputs.py" \
  --matrices "$network_matrices" --metadata "$network_metadata" --bundle "$network_out/runtime"
run_step primary "$network_python" -u "$network_scripts/network_analysis.py" \
  --bundle "$network_out/runtime" --out-root "$network_out/results"
network_run="$(cat "$network_out/results/LATEST_809.txt")"
run_step sensitivity "$network_python" -u "$network_scripts/parameter_sensitivity.py" \
  --bundle "$network_out/runtime" --run "$network_run"
run_step deletion "$network_python" -u "$network_scripts/profile_deletion.py" \
  --bundle "$network_out/runtime" --run "$network_run"
run_step audit "$network_python" -u "$network_scripts/validate_networks.py" \
  --bundle "$network_out/runtime" --run "$network_run"
run_step algorithms "$network_python" -u "$network_scripts/algorithm_stability.py" \
  --bundle "$network_out/runtime" --run "$network_run"
run_step algorithm_audit "$network_python" -u "$network_scripts/validate_algorithms.py" \
  --bundle "$network_out/runtime" --run "$network_run"
run_step graphical_lasso "$network_python" -u "$network_scripts/graphical_lasso_screen.py" \
  --run "$network_run"
run_step figures "$network_python" -u "$network_scripts/render_figures.py" \
  --run "$network_run" --output "$network_out/figures"
run_step algorithm_figure "$network_python" -u "$network_scripts/render_algorithm_stability.py" \
  --run "$network_run" --output "$network_out/figures"
echo "PASS | Network outputs: $network_run"
echo "PNG, PDF and SVG: $network_out/figures"
