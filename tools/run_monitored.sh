#!/usr/bin/env bash
# Inputs: command and arguments; optional AUDIT_LOG_DIR.
# Outputs: timestamped log; actual child return status.
# Method: process monitoring every 30 seconds, one numerical thread.
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
set +e
run_monitored_main() {
  if [ "$#" -eq 0 ]; then printf 'Usage: bash run_monitored.sh command [arguments]\n'; return 2; fi
  local task_logs="${AUDIT_LOG_DIR:-local_runs/logs}"
  mkdir -p -- "$task_logs" || return 1
  local task_log="$task_logs/run_$(date -u +%Y%m%d_%H%M%S)_$$.log"
  env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
    "$@" >"$task_log" 2>&1 &
  local task_pid=$!
  printf 'PID: %s | Log: %s | Monitor: 30 s\n' "$task_pid" "$task_log"
  (
    while kill -0 "$task_pid" 2>/dev/null; do
      date -u '+%F %T UTC'
      ps -p "$task_pid" -o pid,stat,etime,time,pcpu,pmem,comm 2>/dev/null
      tail -n 8 "$task_log"
      sleep 30
    done
  ) &
  local task_monitor=$!
  wait "$task_pid"
  local task_status=$?
  kill "$task_monitor" 2>/dev/null
  wait "$task_monitor" 2>/dev/null
  tail -n 25 "$task_log"
  printf '\nReturn status: %s | Log: %s\n' "$task_status" "$task_log"
  return "$task_status"
}
run_monitored_main "$@"
