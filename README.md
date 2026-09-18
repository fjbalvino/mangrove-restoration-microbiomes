# Mangrove restoration microbiomes

Code and provenance for **Locality and depth structure sediment microbiomes across
mangrove restoration stages**: 51 metagenomes, 17 profiles, three localities and
three sediment depths (5, 20 and 40 cm).

## Analysis workflow

| Directory | Analysis | Manuscript output |
|---|---|---|
| `00_preflight` | Sample identity and complete-profile checks | Cohort |
| `01_environment_indices` | Environmental axes, calibrated MHI, archived HI | Methods; index supplements |
| `02_metagenomics` | Gene-count filtering and KO/PFAM aggregation | Functional inputs |
| `03_taxonomic_structure` | Screening, diversity, composition, variance partition | Fig. 1 |
| `04_taxonomic_change_points` | Global and depth-resolved TITAN2 | Fig. 2 |
| `05_functional_reference` | Functional composition and local reference proximity | Fig. 3 |
| `06_gene_associations` | Protein CLR, mixed models and annotation | Fig. 4 |
| `07_networks` | Final residualized ridge/LIONESS analysis, 809 | Network supplementary analyses |
| `08_figures_tables` | Frozen-table figure reproduction | Fig. 4 |
| `supplementary/hfr` | Historical HFR audit and held-out evaluation | Excluded recovery claim; provenance only |

Directories are under `workflow/`. New filenames make manuscript order explicit;
historical output IDs remain unchanged to preserve lineage. The index scripts
are not a shell pipeline: HI uses an archived calibration, not the current
taxon-selection output. Consult the contracts before executing individual files.

## Reproduce the available components

These commands check or reproduce three components of the analysis. Run them
from the root of the complete analysis package, with its scripts and input
files present. The GitHub `main` branch currently contains only this README;
the directories referenced below must also be available locally.

### 1. Verify the historical HI

This command reconstructs the archived HI from its original CLR matrix and
metadata. It removes the fitted effects of locality and categorical depth,
obtains the first principal component of the residuals, and orients its sign
using the preserved and degraded reference groups. It then compares the
reconstructed values with the archived HI and checks that the archived values
match the HI column in the integrated metadata.

Requires Python, NumPy and pandas, plus the four CSV files in
`source_data/indices/`.

```bash
set +e
bash tools/run_monitored.sh python3 checks/verify_historical_hi.py
```

**Output:** `validation/HI_reconstruction_audit.json`, with numerical differences
and a PASS/FAIL result at a maximum absolute-error tolerance of `1e-9`.
This checks consistency with the historical calibration; it does not provide
independent biological validation of the index.

### 2. Regenerate Figure 4

This command verifies the bundled input files and redraws Figure 4, including
panels A–D, from the saved statistical tables and plotting configuration.
It uses the existing model results without refitting the upstream gene models.

Requires Python, NumPy, pandas and Matplotlib, plus the contents of
`workflow/08_figures_tables/figure4_bundle/`.

```bash
set +e
bash workflow/08_figures_tables/08_02_run_figure4.sh local
```

**Output:** a new timestamped directory under
`workflow/08_figures_tables/figure4_bundle/outputs/`, containing the figure,
individual panels, input copies and run metadata. The `local` argument selects
this output location; `FIG4_OUT_ROOT` can override it.

### 3. Recalculate the functional networks

This command runs the 809 analysis using the included 200-KO CLR matrix and
metadata for 51 samples. It compares three network constructions: unadjusted;
adjusted for locality and depth; and additionally adjusted for the four
environmental axes. It estimates ridge partial correlations and sample-specific
LIONESS networks, then evaluates network metrics, their associations with
HI/MHI, and their sensitivity using profile-level resampling and permutations.

Requires the input files in `workflow/07_networks/network809/inputs/` and the
package versions specified in
`workflow/07_networks/network809/requirements.txt`.

```bash
set +e
bash workflow/07_networks/07_01_run_residual_networks.sh
```

**Output:** a timestamped analysis directory and execution log under
`workflow/07_networks/network809/results/`. Set `OUT809` to change the output
root or `PY809` to select a Python interpreter. The resampling makes this a
substantial calculation; repository checks do not launch it automatically.

The launchers limit numerical libraries to one thread and report progress
every 30 seconds. `set +e` disables automatic shell termination on command
failure; it does not suppress errors or make a failed analysis successful.
In the complete package, [reproduction instructions](docs/REPRODUCIBILITY.md)
describe landalab paths and additional data requirements.

## Inputs, outputs and transparency

Every selected scientific script has a documentation header, its source hash and
method attribution. Dynamic expressions are explicitly labelled; they are not
invented resolved file paths. Complete static contracts are in `docs/contracts/`.
`docs/INPUTS.tsv`, `OUTPUTS.tsv`, `FROZEN_FILES.tsv` and `MISSING_INPUTS.tsv`
distinguish input expressions, saved evidence and unavailable bytes.

The original capsule SHA-256 is
`9c3707137004532fa2bc5fe582b4a1ca9e46c8a714987ee6e7a6db5cd905e59b`.
All 3,291 entries in its checksum manifest were verified on receipt.

This review archive retains the authors' source code and does not grant a new
licence over pre-existing work. See `docs/AI_ASSISTANCE.md` for curation provenance.
Raw sequence data are reported under BioProject PRJNA1502944 in the manuscript;
public availability of those reads has not been checked in this code audit.
