#!/usr/bin/env bash
# Inputs: convergence808/inputs/ (KO/PFAM counts and aligned metadata).
# Outputs: OUT808/808_functional_convergence_<timestamp>/ (default local_runs/808).
# Algorithmic provenance: CLR/Aitchison; Freedman-Lane profile permutations;
# within-locality profile bootstrap; declared Holm and BH families.
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
set +e
run_convergence_main() {
  local task_dir task_repo task_out
  task_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd) || return 1
  task_repo=$(cd -- "$task_dir/../.." && pwd) || return 1
  task_out="${OUT808:-$task_repo/local_runs/808}"
  bash "$task_repo/tools/run_monitored.sh" "${PY808:-python3}" \
    "$task_dir/convergence808/scripts/808_functional_convergence.py" \
    --bundle "$task_dir/convergence808" --out-root "$task_out"
}
run_convergence_main "$@"
