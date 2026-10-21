# Inputs and execution

## Delivery scope

The current delivery comprises code, configuration, metadata and the shipped
inputs listed in `INPUTS_MINIMAL.tsv`. The corrected read-based KO/PFAM matrices
required by section 05 are supplied separately; their expected names, sizes and
hashes are recorded in `FUNCTIONAL_INPUTS_REQUIRED.tsv`. Large gene-level source matrices remain
on the analysis server for a separate later upload, as requested by the project
owner. Their transfer is not required to finish this code delivery, but the
analyses that depend on them still require access to those files.

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
| Historical networks 809 | Fixed 200-KO CLR matrix and integrated metadata | Read-based update pending; prior fitted objects not required |
| Protein models | Full filtered protein counts and four 044 metadata/QC controls | Still required from the analysis server |
| Gene annotation | Matching annotation catalogue and regenerated model-audit inputs | Server input/preparation still required |
| Figure 4 | Regenerated model-statistic tables and plotting configuration | Plotting code included; fitted result tables deliberately absent |

`tools/prepare_runtime_inputs.py` makes byte-identical local copies for the two
Python analysis bundles; it refuses to replace conflicting runtime inputs.
Primary functional composition reads the canonical matrix paths directly.

The corrected section 05 is a runnable pipeline after supplying the four verified matrices;
see `FUNCTIONAL_REFERENCE.md`. Its launcher does not read historical 082A matrices.
Other R analyses retain historical CLI defaults and output names. Consult their
command-line arguments and `docs/contracts/` before executing those stages.

## Protein CLR preparation

`06_01_full_catalogue_CLR.py` requires `--run044` containing the full filtered
count matrix and four metadata/QC tables under `tables/`. Its `--run003b` must
contain `003_metadata_integrada_canon_51.csv` and `003b_run_summary.tsv` under
`tables/`. The former is byte-identical to the included `metadata_003b.csv`.
Pass the included `202_metadata_raw1031_alineada.csv` to `--canon`.
These are checked against the original hashes and QC totals; missing controls
must not be fabricated. It generates the selected 100,000-gene CLR and count matrices.

## Annotation and figure preparation

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
No automated conversion from a fresh complete gene-model run to every figure
input has been demonstrated in this minimal publication.

## Historical predictor lineage

The archived TITAN MHI branch used pre-003b MHI values, and the archived depth-wise
variance partition used an ECI vector different from current integrated metadata.
Choose and document the intended calibration before rerunning those branches.
Publishing source code and inputs does not resolve these scientific lineage differences.
