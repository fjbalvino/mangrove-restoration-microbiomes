# Mangrove restoration microbiomes

Code and provenance for **Locality and depth structure sediment microbiomes across
mangrove restoration stages**: 51 metagenomes, 17 profiles, three localities and
three sediment depths (5, 20 and 40 cm).

**Status: audited review candidate, not a complete end-to-end reproduction release.**
The current manuscript and supplementary PDFs were compared with the collected
code and saved results. Two predictor-lineage discrepancies require resolution:
the TITAN2 MHI branch used pre-003b values, and the depth-wise variance partition
used an ECI vector different from the current integrated metadata. See
[audit findings](docs/AUDIT_FINDINGS.md) before using archived statistics.

## Scientific order

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

## Which code is current?

`scripts_finales/` and code snapshots inside `resultados_finales/` take precedence.
The exact code saved with a result establishes its producer; a later working copy
is not retrospectively attributed to that result. Older `Scripts/` and `Results/`
files are retained only for verified upstream provenance (notably archived HI and
the 044/082A matrices). The separately delivered **809 residualization package**
is authoritative for networks. It ran in GPT, not landalab.

The audit captured 523 script copies representing 410 unique SHA-256 versions.
All are catalogued in `docs/SCRIPT_INDEX.tsv`; 33 analysis/figure source versions
are explicitly mapped to manuscript sections, in addition to the final 809 bundle
and new audit utilities. Static cataloguing is not line-by-line validation of all
410 programs. The mapping and exact-copy identities are available in
`docs/METHODS_CODE_CROSSWALK.tsv` and `docs/CURRENT_SOURCE_LINEAGE.tsv`.

## Reproduce the available components

With Python, NumPy and pandas available:

```bash
set +e
bash tools/run_monitored.sh python3 checks/verify_historical_hi.py
```

Figure 4 from its included audited tables:

```bash
set +e
bash workflow/08_figures_tables/08_02_run_figure4.sh local
```

Final networks (requires the versions in the bundled `requirements.txt`):

```bash
set +e
bash workflow/07_networks/07_01_run_residual_networks.sh
```

The network calculation is substantial; it is not launched automatically by
repository checks. Launchers use one numerical thread and monitor every 30 s.
For landalab paths and the next targeted data collection, see
[reproduction instructions](docs/REPRODUCIBILITY.md).

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
