#!/usr/bin/env bash
# Reassembled read-derived export -> full-catalogue CLR -> gene models.
# No reads, databases, fitted models or figures are uploaded by this script.
set -euo pipefail
if [[ $# -ne 2 ]]; then
  echo 'Usage: bash workflow/06_gene_associations/06_05_run_read_based.sh EXTERNAL_INPUT_DIR NEW_ANALYSIS_DIR' >&2
  exit 2
fi
task_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ -e "$2" ]]; then
  echo "Analysis directory already exists; choose a fresh directory: $2" >&2
  exit 2
fi
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
python3 "$task_dir/read_based/prepare_gene_matrix.py" --input "$1" --analysis-root "$2"
python3 "$task_dir/read_based/annotate_selected.py" --input "$1" --analysis-root "$2"
python3 "$task_dir/read_based/fit_random_intercept.py" --analysis-root "$2"
python3 "$task_dir/read_based/summarize_models.py" --analysis-root "$2"
echo 'Gene-model tables written. Independent numerical validation remains required for a new full fit.'
