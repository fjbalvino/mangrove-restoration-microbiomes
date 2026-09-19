#!/usr/bin/env bash
# Inputs: source_data/functional/tables/082A_*tsv.gz; indices/metadata_003b.csv.
# Outputs: FUNCTIONAL_OUT_ROOT/602_validar_composicion_funcional_splitplot_*/.
# Algorithmic provenance: technical filters, CLR/Aitchison, profile-restricted
# permutations, dispersion and matched-reference distances.
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
# Calls the unchanged primary 602 producer with explicit included inputs.
set +e
run_primary_functional_main() {
  local task_dir task_repo task_r
  task_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd) || return 1
  task_repo=$(cd -- "$task_dir/../.." && pwd) || return 1
  task_r="${FUNCTIONAL_RSCRIPT:-Rscript}"
  if ! command -v "$task_r" >/dev/null 2>&1; then
    printf 'Rscript is unavailable. Set FUNCTIONAL_RSCRIPT to the R interpreter.\n'
    return 1
  fi
  bash "$task_repo/tools/run_monitored.sh" "$task_r" \
    "$task_dir/05_01_composition_reference.R" \
    --run082 "$task_repo/source_data/functional" \
    --meta "$task_repo/source_data/indices/metadata_003b.csv" \
    --out_root "${FUNCTIONAL_OUT_ROOT:-$task_repo/local_runs/functional_reference}"
}
run_primary_functional_main "$@"
