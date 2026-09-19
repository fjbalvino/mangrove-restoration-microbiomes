#!/usr/bin/env bash
# Inputs: figure4_bundle/inputs and config/plot.json (derived tables must be generated first).
# Outputs: FIG4_OUT_ROOT/604d_reproducir_figura4_A_D_*/.
# Method: plots frozen mixed-model results; no model refitting.
# AI-assisted curation: OpenAI Codex (OpenAI, 2026).
set +e
FIG4_TASK_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
bash "$FIG4_TASK_DIR/figure4_bundle/run_figure4.sh" "${1:-local}"
