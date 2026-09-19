#!/usr/bin/env bash
# Inputs: network809/inputs/{clr_200_KO.tsv.gz,metadata.csv}.
# Diagnostics additionally consume the completed 809 run (objects, draws and tables).
# Outputs: configured results root/809_residualized_functional_networks_*/; see PROTOCOL.md.
# Provenance: fixed-lambda ridge + LIONESS (Kuijjer et al., 2019, doi:10.1016/j.isci.2019.03.021).
# Curation: OpenAI Codex (OpenAI, 2026); documentation only; algorithm body preserved.
# Module filenames remain unchanged to preserve Python imports and historical source-hash checks.
set +e
BUNDLE809=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export OUT809=${OUT809:-/home/fjbalvino/Tipping_points/resultados_finales}
bash "$BUNDLE809/scripts/run_809.sh"
