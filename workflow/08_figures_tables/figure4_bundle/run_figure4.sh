#!/usr/bin/env bash
# Run as: bash run_figure4.sh local | landalab
# No background R jobs, no changes to original scripts, no interactive-shell exit.
set +e

run_figure4_main() {
  local fig4_mode="${1:-local}"
  local fig4_root
  local fig4_python="${FIG4_PYTHON:-python3}"
  local fig4_output
  fig4_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)" || return 1
  case "$fig4_mode" in
    local) fig4_output="$fig4_root/outputs" ;;
    landalab) fig4_output="/home/fjbalvino/Tipping_points/resultados_finales" ;;
    *) printf 'Uso: bash run_figure4.sh local|landalab\n'; return 2 ;;
  esac
  fig4_output="${FIG4_OUT_ROOT:-$fig4_output}"
  if ! command -v "$fig4_python" >/dev/null 2>&1; then
    printf 'Python no encontrado: %s. Define FIG4_PYTHON con la ruta del entorno.\n' "$fig4_python"
    return 1
  fi
  "$fig4_python" "$fig4_root/code/verify_package.py" --root "$fig4_root" || return 1
  printf 'Modo: %s\nPython: %s\nSalidas: %s\nMonitoreo interno cada 30 s.\n' "$fig4_mode" "$fig4_python" "$fig4_output"
  # Bound numerical-library threads even when the user's server is busy.
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
    "$fig4_python" "$fig4_root/code/reproduce_figure4.py" \
      --source-dir "$fig4_root/inputs" \
      --config "$fig4_root/config/plot.json" \
      --out-root "$fig4_output"
  local fig4_status=$?
  printf 'Codigo de retorno: %s\n' "$fig4_status"
  return "$fig4_status"
}
run_figure4_main "$@"
