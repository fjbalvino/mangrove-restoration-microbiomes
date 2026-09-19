#!/usr/bin/env bash
# Inputs: network809/inputs/{clr_200_KO.tsv.gz,metadata.csv}.
# Diagnostics additionally consume the completed 809 run (objects, draws and tables).
# Outputs: configured results root/809_residualized_functional_networks_*/; see PROTOCOL.md.
# Provenance: fixed-lambda ridge + LIONESS (Kuijjer et al., 2019, doi:10.1016/j.isci.2019.03.021).
# Curation: OpenAI Codex (OpenAI, 2026); launcher propagates failures; analysis unchanged.
# Module filenames remain unchanged to preserve Python imports and historical source-hash checks.
set +e
run_diagnostic_main() {
BUNDLE809=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
PY809=${PY809:-python3}
TASK809=${1:-audit}
OUT809=${OUT809:-"$BUNDLE809/results"}
RUN809=${RUN809:-$(head -n 1 "$OUT809/LATEST_809.txt")}
case "$TASK809" in
 audit) SCRIPT809=809_validate.py ;;
 deletion) SCRIPT809=809_profile_deletion.py ;;
 *) printf 'Use audit or deletion.\n'; return 2 ;;
esac
if [ -n "$SCRIPT809" ] && [ -d "$RUN809" ]; then
 mkdir -p "$RUN809/logs"
 LOG809="$RUN809/logs/${TASK809}_$(date -u +%Y%m%d_%H%M%S).log"
 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$PY809" "$BUNDLE809/scripts/$SCRIPT809" --bundle "$BUNDLE809" --run "$RUN809" > "$LOG809" 2>&1 &
 PID809=$!
 (
  while kill -0 "$PID809" 2>/dev/null; do
   date -u '+%F %T UTC'
   ps -p "$PID809" -o pid,stat,etime,time,pcpu,pmem,comm 2>/dev/null || printf 'Diagnostic PID %s active\n' "$PID809"
   tail -n 6 "$LOG809"
   sleep 30
  done
 ) &
 MON809=$!
 wait "$PID809"
 STATUS809=$?
 kill "$MON809" 2>/dev/null
 wait "$MON809" 2>/dev/null
 printf '\nDiagnostic status: %s\nLog: %s\n' "$STATUS809" "$LOG809"
 tail -n 20 "$LOG809"
else
 printf 'No valid run directory: %s\n' "$RUN809"
fi

return "${STATUS809:-1}"
}
run_diagnostic_main "$@"
