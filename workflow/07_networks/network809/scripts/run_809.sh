#!/usr/bin/env bash
# Inputs: network809/inputs/{clr_200_KO.tsv.gz,metadata.csv}.
# Diagnostics additionally consume the completed 809 run (objects, draws and tables).
# Outputs: configured results root/809_residualized_functional_networks_*/; see PROTOCOL.md.
# Provenance: fixed-lambda ridge + LIONESS (Kuijjer et al., 2019, doi:10.1016/j.isci.2019.03.021).
# Curation: OpenAI Codex (OpenAI, 2026); documentation only; algorithm body preserved.
# Module filenames remain unchanged to preserve Python imports and historical source-hash checks.
set +e
run_809_main() {
BUNDLE809=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
PY809=${PY809:-python3}
OUT809=${OUT809:-"$BUNDLE809/results"}
mkdir -p "$OUT809"
LOG809="$OUT809/809_launch_$(date -u +%Y%m%d_%H%M%S).log"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$PY809" "$BUNDLE809/scripts/809_residual_networks.py" --bundle "$BUNDLE809" --out-root "$OUT809" > "$LOG809" 2>&1 &
PID809=$!
(
 while kill -0 "$PID809" 2>/dev/null; do
  date -u '+%F %T UTC'
  ps -p "$PID809" -o pid,stat,etime,time,pcpu,pmem,comm 2>/dev/null || printf 'Analysis PID %s active\n' "$PID809"
  tail -n 8 "$LOG809"
  sleep 30
 done
) &
MON809=$!
wait "$PID809"
STATUS809=$?
kill "$MON809" 2>/dev/null
wait "$MON809" 2>/dev/null
printf '\nAnalysis status: %s\nLog: %s\n' "$STATUS809" "$LOG809"
tail -n 25 "$LOG809"

# Curation correction: propagate the analysis status; no interactive-shell exit.
return "$STATUS809"
}
run_809_main "$@"
