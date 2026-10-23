# Mangrove restoration microbiomes

Analysis code and essential available inputs for **Locality and depth structure
sediment microbiomes across mangrove restoration stages**: 51 metagenomes,
17 sediment profiles, three localities and depths of 5, 20 and 40 cm.

This repository contains scripts, configuration, metadata, fixed index-calibration
inputs and functional matrices. Figures, PDFs, fitted models, permutation/bootstrap
draws and full statistical result tables are generated locally. A compact numerical
fixture is included to verify the corrected functional section.

Section 05 now uses the corrected read-based KO/PFAM matrices, with separate
fragment-count and length-adjusted inputs supplied separately. With those inputs,
it reproduces Figure 3, functional cores,
Supplementary Tables 6 and 10, and Supplementary Figure 9. See the
[section-specific methods and execution guide](docs/FUNCTIONAL_REFERENCE.md).
The corrected read-based gene models and environmental functional enrichment
are now in sections 06 and 09. They reproduce the current Figures 4 and 5 using
external inputs. See [GENE_ENVIRONMENT.md](docs/GENE_ENVIRONMENT.md) for commands,
input requirements, validation and interpretation limits. This update adds only
code, configuration and documentation: no data matrices, annotation databases,
figures or binary archives. Section 10 adds the targeted carbon, nitrogen and
sulfur marker analysis that follows environmental enrichment in the manuscript,
with Figure 6, Supplementary Figure 10 and Supplementary Tables 16–17. See
[FUNCTIONAL_CYCLES.md](docs/FUNCTIONAL_CYCLES.md). This exploratory section
replaces networks in the manuscript; all section 07 code remains historical.
Ordered restoration trends remain a separate branch.

## Analysis workflow

| Directory under `workflow/` | Purpose |
|---|---|
| `00_preflight` | Sample identity and complete-profile checks |
| `01_environment_indices` | Environmental axes and MHI; historical HI calibration retained |
| `02_metagenomics` | Protein filtering and KO/PFAM aggregation |
| `03_taxonomic_structure` | Taxonomic screening, diversity and composition |
| `04_taxonomic_change_points` | Global and depth-wise TITAN2 |
| `05_functional_reference` | Functional composition and reference proximity |
| `06_gene_associations` | Corrected read-based gene CLR, mixed models and annotation; older entry points retained as historical |
| `07_networks` | Historical networks and algorithm audits, withdrawn from the current manuscript |
| `08_figures_tables` | Current Figures 4 and 5 from verified tables; historical Figure 4 bundle retained separately |
| `09_environment_enrichment` | KO/PFAM overrepresentation and potential KEGG modules by environmental axis and sign |
| `10_functional_cycles` | Candidate C/N/S markers, profile-level inference, quantification sensitivities and local Figure 6 / Supplementary Figure 10 |

## Included inputs

- Integrated 51-sample metadata and the fixed historical HI/MHI calibration inputs.
- Quantification provenance for corrected KO/PFAM matrices supplied separately.
- Historical assembled-ORF matrices used only by the unrecalculated 808 branch.
- The historical 51-sample × 200-KO CLR matrix, retained for archived networks only.
- Small sample-order and provenance controls required by the gene-CLR producer.

[INPUTS_MINIMAL.tsv](docs/INPUTS_MINIMAL.tsv) records each input's role, bytes and
SHA-256. Historical calibration matrices are retained because fitting new indices
would change the study predictors.

## Run the available components

Install the versions in `requirements-analysis.txt`, then run from this root.
Before section 05, supply the four matrices listed in
[FUNCTIONAL_INPUTS_REQUIRED.tsv](docs/FUNCTIONAL_INPUTS_REQUIRED.tsv) as described
in the section guide:

```bash
set +e
python3 checks/check_minimal_repository.py
python3 checks/verify_historical_hi.py
python3 checks/verify_environmental_mhi.py
python3 tools/prepare_runtime_inputs.py
```

The staging command copies verified inputs to the legacy bundle paths expected
by the Python modules; those duplicate runtime copies are ignored by Git.

```bash
set +e
# Corrected composition, reference proximity, cores, validation and local figures.
bash workflow/05_functional_reference/05_00_run_primary_composition.sh
# Historical 808 reproduction; still uses assembled-ORF matrices.
bash workflow/05_functional_reference/05_04_run_convergence_sensitivity.sh
# Targeted C/N/S markers; supply the compact verified inputs externally.
python3 workflow/10_functional_cycles/10_00_run_cycles.py \
  --input-dir /path/to/CNS/inputs --output-dir local_runs/functional_cycles
```

The corrected section writes to `local_runs/functional_reference/` and refuses
to overwrite an existing primary run. Set `FUNCTIONAL_OUT_ROOT` to a fresh local
directory for a repeat. Progress is reported every 30 seconds; failures stop the
pipeline. Other launchers retain their historical output conventions.
The marker launcher verifies input hashes and the 17 complete profiles, then
creates a new local run with its exact code and configuration snapshot. It does
not require network outputs. Optional historical network commands are in
[NETWORKS.md](docs/NETWORKS.md).

## External inputs

The four corrected KO/PFAM input matrices are external inputs, not included in
this code update. Their expected filenames, sizes and SHA-256 values are in
`docs/FUNCTIONAL_INPUTS_REQUIRED.tsv`; quantification provenance is included.
Place the verified TSV files in `source_data/functional/read_based/` before
running section 05. The launcher reports missing or altered matrices and stops.
Once those matrices are supplied, reproducing section 05 does not require
the large gene catalogue or another read-mapping run. Rebuilding quantification
from FASTQ still requires the reads, reference sequences, mapping outputs and
annotation links held on the analysis server.

The canonical `001_phyloseq_canon_51.rds` and full
`044_protein_filtered_counts_51samples.tsv.gz` were not available for inclusion.
The unchanged gene-CLR producer also checks four upstream 044 metadata/QC files.
[MISSING_INPUTS.tsv](docs/MISSING_INPUTS.tsv) lists the six exact paths.

```bash
set +e
# Later, on landalab: collect the six listed inputs and compute SHA-256.
# This command creates a local archive; it does not upload data.
bash tools/run_monitored.sh python3 tools/collect_required_inputs.py --group core
```

The `06_01`–`06_04` entry points describe the historical assembled-ORF analysis.
For the current gene results, use `06_05_run_read_based.sh` and its `read_based/`
scripts. They centre against 12,540,513 retained read-derived genes before selecting
100,000 by CLR variance. Rebuilding these inputs needs the external per-sample
gene export and historical annotation dictionary; enrichment can instead start
from verified complete model tables and tested-gene annotations.

```bash
# Paths below refer to separately supplied files; the repository does not contain them.
bash workflow/09_environment_enrichment/09_00_run_enrichment.sh \
  /path/to/Figura4 /path/to/Supplementary_Data1/reference \
  local_runs/environment_enrichment

# Render current Figure 4 from its external, verified source tables.
python3 workflow/08_figures_tables/08_03_read_based_figure4.py \
  --figure-data /path/to/Supplementary_Data1/figure_data \
  --output local_runs/current_figures
```

The enrichment launcher also renders Figure 5 locally in PNG, PDF and SVG.
The older `08_01`, `08_02` and `figure4_bundle/` are historical and do not produce
the current manuscript figure. No script uploads its local outputs.

See [execution requirements](docs/REPRODUCIBILITY.md),
[methods](docs/METHODS.md), [code changes](docs/MINIMAL_CHANGES.md), and
[Español](README_ES.md). A complete raw-read-to-paper rerun is not claimed.
