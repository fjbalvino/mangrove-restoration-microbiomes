#!/usr/bin/env bash
# Inputs: convergence808/inputs/; validator also needs a completed run.
# Outputs: results root/808_functional_convergence_<timestamp>/ or reference-run audit.
# Algorithmic provenance: CLR/Aitchison; Freedman-Lane profile permutations;
# locality-stratified profile bootstrap; Holm/BH families in PROTOCOL.md.
# Source SHA-256: 59421fe696e5b503b350113a451515b4144e46bf41f3b90432f7c3464dfbf026
# AI-assisted code curation: OpenAI Codex (OpenAI, 2026).
set +e
BUNDLE808=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export OUT808=${OUT808:-/home/fjbalvino/Tipping_points/resultados_finales}
bash "$BUNDLE808/scripts/run_808.sh"
