#!/usr/bin/env bash
# External inputs only. All results and figures remain local.
set -euo pipefail
if [[ $# -ne 3 ]]; then
  echo 'Usage: bash workflow/09_environment_enrichment/09_00_run_enrichment.sh ANALYSIS_ROOT REFERENCE_DIR OUTPUT_DIR' >&2
  exit 2
fi
task_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$task_dir/../.." && pwd)
analysis_root=$1
reference_dir=$2
output_dir=$3
if [[ -e "$output_dir" ]]; then
  echo "Output already exists; choose a fresh directory: $output_dir" >&2
  exit 2
fi
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
python3 "$repo_root/checks/check_environment_enrichment_inputs.py" \
  --analysis-root "$analysis_root" --reference "$reference_dir"
python3 "$task_dir/explore_enrichment.py" --source "$analysis_root" --output "$output_dir"
python3 "$task_dir/explore_modules.py" --source "$analysis_root" --reference "$reference_dir" --output "$output_dir"
python3 "$task_dir/explore_enrichment.py" --source "$analysis_root" --output "$output_dir/single_KO_functions" --single-ko
python3 "$task_dir/explore_modules.py" --source "$analysis_root" --reference "$reference_dir" --output "$output_dir/single_KO_modules" --single-ko
python3 "$task_dir/prepare_tables.py" --results "$output_dir" --reference "$reference_dir"
python3 "$repo_root/checks/verify_environment_enrichment.py" --results "$output_dir"
python3 "$repo_root/workflow/08_figures_tables/08_04_environment_enrichment_figure5.py" \
  --results "$output_dir" --output "$output_dir/figures"
