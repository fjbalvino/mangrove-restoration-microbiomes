#!/usr/bin/env bash
# Inputs: convergence808/inputs/; validator also needs a completed run.
# Outputs: results root/808_functional_convergence_<timestamp>/ or reference-run audit.
# Algorithmic provenance: CLR/Aitchison; Freedman-Lane profile permutations;
# locality-stratified profile bootstrap; Holm/BH families in PROTOCOL.md.
# Source SHA-256: ad2e3f1aa548a4f204823c059e3ac8452c61c75fb77fc237f3a0c7f424f51dc2
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Launcher correction: return the analysis status after displaying the log.
set +e
run_808_main() {
BUNDLE808=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
PY808=${PY808:-python3}
OUT808=${OUT808:-"$BUNDLE808/results"}
mkdir -p "$OUT808"
LOG808="$OUT808/808_launch_$(date -u +%Y%m%d_%H%M%S).log"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 "$PY808" "$BUNDLE808/scripts/808_functional_convergence.py" --bundle "$BUNDLE808" --out-root "$OUT808" > "$LOG808" 2>&1 &
PID808=$!
(
 while kill -0 "$PID808" 2>/dev/null; do
  date -u '+%F %T UTC'
  ps -p "$PID808" -o pid,stat,etime,time,pcpu,pmem,comm 2>/dev/null
  tail -n 8 "$LOG808"
  sleep 30
 done
) &
MON808=$!
wait "$PID808"
STATUS808=$?
kill "$MON808" 2>/dev/null
wait "$MON808" 2>/dev/null
printf '\nAnalysis status: %s\nLog: %s\n' "$STATUS808" "$LOG808"
tail -n 25 "$LOG808"

return "$STATUS808"
}
run_808_main "$@"
