# Inputs and execution

## Delivery scope

The current delivery comprises code, configuration, metadata and the shipped
inputs listed in `INPUTS_MINIMAL.tsv`. The corrected read-based KO/PFAM matrices
required by section 05 are supplied separately; their expected names, sizes and
hashes are recorded in `FUNCTIONAL_INPUTS_REQUIRED.tsv`. The current gene/environment
update adds code and documentation only. Gene matrices, annotation snapshots,
model-result tables and figures remain external. Their transfer is not required
to finish this code delivery, but execution still requires those files.
`GENE_ENVIRONMENT.md` describes the current sections 06, 08 and 09.

The canonical phyloseq and upstream metadata/QC files listed in
`MISSING_INPUTS.tsv` also remain to be supplied. The collection command in the
README creates a local archive on the server; it does not publish that archive.
After the data upload, record the actual filenames, download locations, byte
sizes and SHA-256 checksums. Archived inventory sizes should be checked against
the files actually transferred.

## Available entry points

| Component | Required source input | Current availability |
|---|---|---|
| Historical HI | Historical CLR, metadata and fixed HI values | Included; numerical reconstruction verified |
| Environmental MHI | Original 57-row calibration and 51-row metadata | Included; numerical reconstruction verified |
| Taxonomic analyses | Canonical phyloseq plus metadata | Phyloseq still required from the analysis server |
| Corrected functional composition and core | Read-based KO/PFAM fragment counts, length-adjusted values and canonical metadata | Code included and validated with local inputs; four read-based matrices must be supplied separately |
| Historical supplementary 808 | Assembled-ORF matrices, current and historical labels | Historical reproduction only; read-based update pending |
| Historical networks 809 | Fixed historical 200-KO CLR and metadata | Historical reproduction only |
| Current networks | Corrected KO RPK/count matrices and calibrated 51-row metadata | Current code, algorithm benchmark and independent validation included; matrices external |
| Current read-based gene models | Per-sample gene export, QC controls and historical annotation dictionary | Code included; full source inputs supplied externally; legacy 06_01–06_04 retained separately |
| Environmental enrichment | Complete verified common-slope tables, tested-gene annotations and KEGG/PFAM snapshots | Current code and input hashes included; data and databases external |
| Current Figures 4 and 5 | Verified gene-model/figure tables and enrichment tables | Current renderers included; PNG/PDF/SVG generated locally |

`tools/prepare_runtime_inputs.py` makes byte-identical local copies for the two
Python analysis bundles; it refuses to replace conflicting runtime inputs.
Primary functional composition reads the canonical matrix paths directly.

The corrected section 05 is a runnable pipeline after supplying the four verified matrices;
see `FUNCTIONAL_REFERENCE.md`. Its launcher does not read historical 082A matrices.
Other R analyses retain historical CLI defaults and output names. Consult their
command-line arguments and `docs/contracts/` before executing those stages.

## Current network execution

Use `07_02_run_read_based.sh` with the two corrected KO matrices and canonical
metadata. It selects current features, rebuilds networks and resampling, compares
parameter settings and algorithms, independently audits the primary results and
renders PNG/PDF/SVG locally. See `NETWORKS.md` for commands and limitations.
The old `07_00`, `07_01` and `network809/` remain historical.

## Current gene/environment execution

Use `06_05_run_read_based.sh` for corrected gene preparation and models and
`09_00_run_enrichment.sh` for the verified Figure 4A enrichment. The current
renderers are `08_03_read_based_figure4.py` and
`08_04_environment_enrichment_figure5.py`. Their commands, external input layout,
validation scope and interpretation limits are in `GENE_ENVIRONMENT.md`.
The earlier gene scripts and Figure 4 bundle described below are historical.

## Historical protein CLR preparation

`06_01_full_catalogue_CLR.py` requires `--run044` containing the full filtered
count matrix and four metadata/QC tables under `tables/`. Its `--run003b` must
contain `003_metadata_integrada_canon_51.csv` and `003b_run_summary.tsv` under
`tables/`. The former is byte-identical to the included `metadata_003b.csv`.
Pass the included `202_metadata_raw1031_alineada.csv` to `--canon`.
These are checked against the original hashes and QC totals; missing controls
must not be fabricated. It generates the selected 100,000-gene CLR and count matrices.

## Historical annotation and Figure 4 bundle

The annotation script exposes `--annotation`, `--run-603c` and `--audit-dir`.
Its archived catalogue is
`FINAL_ANNOTATED_ABUNDANCE_noMZ_bacteria.tsv.gz`; provide its local location through
`--annotation`.
The full catalogue is not bundled. A reduced annotation table must retain all
required fields, all tested gene IDs and validated mapping/duplicate semantics.
The script's verified 603c audit inputs must be regenerated before annotation.

Figure rendering requires `Figure4_counts_q010_q005.tsv`,
`Figure4_interactions_original_statistics.tsv`,
`Figure4_two_candidates_sample_data.tsv` and their `SHA256SUMS.json` in the
figure bundle's runtime `inputs/` directory. These are derived analysis products,
not independent source data, and are intentionally excluded from this repository.
This describes only the historical bundle. The corrected renderer can read
current model tables directly; see the current gene/environment execution above.

## Historical predictor lineage

The archived TITAN MHI branch used pre-003b MHI values, and the archived depth-wise
variance partition used an ECI vector different from current integrated metadata.
Choose and document the intended calibration before rerunning those branches.
Publishing source code and inputs does not resolve these scientific lineage differences.
