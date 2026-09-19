#!/usr/bin/env bash
# Inputs: network809/inputs/; Outputs: OUT809 or network809/results/.
# Method: residualisation, fixed ridge, LIONESS and profile resampling; see network809/PROTOCOL.md.
# AI-assisted curation: OpenAI Codex (OpenAI, 2026). Canonical entry point.
set +e
TASK_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
bash "$TASK_DIR/network809/scripts/run_809.sh"
