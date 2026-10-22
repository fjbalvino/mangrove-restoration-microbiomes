# Read-based genes and environmental functional enrichment

This is the code release for the current Figure 4 gene associations and Figure 5
functional interpretation in *Locality and depth structure sediment microbiomes
across mangrove restoration stages*. It contains no new matrices, databases,
figures, fitted results or binary archives. All analyses write local outputs.

## Current and historical entry points

| Component | Current code | Scope |
|---|---|---|
| Gene preparation and models | `workflow/06_gene_associations/06_05_run_read_based.sh` | Read-derived RPK, full retained-catalogue CLR, selection, REML and BH |
| Functional enrichment | `workflow/09_environment_enrichment/09_00_run_enrichment.sh` | Figure 4A common slopes; KO/PFAM terms and potential modules |
| Figure 4 | `workflow/08_figures_tables/08_03_read_based_figure4.py` | Corrected association counts, depth slopes and two gene examples |
| Figure 5 | `workflow/08_figures_tables/08_04_environment_enrichment_figure5.py` | KO enrichment and exploratory module candidates, blue gradient |

The earlier `06_01`–`06_04` scripts and `08_01`, `08_02` and `figure4_bundle/`
remain historical. They do not produce the current manuscript statistics.
Ordered restoration trends remain separate. Current network execution is
documented in `NETWORKS.md`; it does not change the gene-model/enrichment branch.

## External inputs

Install `requirements-analysis.txt` with Python 3.12. The verified environment
uses NumPy 2.3.5, pandas 2.2.3, SciPy 1.17.0 and Matplotlib 3.10.8. R is needed
only to regenerate the independent numerical benchmark (`lmerTest`, `lme4`,
`numDeriv`); comparison with its saved output uses Python.

Two separately supplied packages contain the runtime inputs:

- `Analisis_Figura4.zip`: complete gene-model tables, metadata, numerical
  validation inputs/results and the original producer scripts.
- `Supplementary_Data1.zip`: tested-gene annotations, frozen Figure 4 source
  tables, enrichment results and reference annotation snapshots.

These packages are external. They are not downloaded by the scripts or included
in this commit. `ENVIRONMENT_ENRICHMENT_INPUTS.json` records expected sizes and
SHA-256 hashes for the inputs used by the enrichment launcher. The reference
snapshots were retrieved on 22 October 2026; do not silently substitute a newer
KEGG or PFAM database. No database content is committed here.

The analysis directory must contain:

```text
prepared/selected_genes.tsv
prepared/metadata.csv
models/{axis}_all_results.tsv.gz
```

The five axes are `vegetation_landscape_PC1`, `water_inundation_PC1`,
`moisture_stress_PC1`, `physicochemical_PC1` and `nutrients_redox_PC1`.
Each full result file contains 100,000 common slopes, 100,000 interaction tests
and 300,000 depth slopes. A list of significant genes alone is insufficient:
the nonsingular tested background is required.

If restoring from the two packages, use `gene_associations/tested_genes_and_annotations.tsv`
from Supplementary Data 1 as `prepared/selected_genes.tsv`. Recreate the prepared
metadata from the unchanged 51-sample export as below. The supplied hashes check
this restored layout, including metadata serialization in the pinned environment.

```bash
python3 - /path/to/Figura4 /path/to/Supplementary_Data1 <<'PY'
from pathlib import Path
import shutil, sys
import pandas as pd
analysis, supplement = map(Path, sys.argv[1:])
prepared = analysis / 'prepared'
prepared.mkdir(exist_ok=True)
shutil.copy2(supplement / 'gene_associations/tested_genes_and_annotations.tsv',
             prepared / 'selected_genes.tsv')
metadata = pd.read_csv(analysis / 'input/analysis_metadata/metadata_003b.csv')
metadata['profile_id'] = metadata.sample_id.str.replace(r'_[123]$', '', regex=True)
metadata.to_csv(prepared / 'metadata.csv', index=False)
PY
```

## Reproduce enrichment and figures

Run from the repository root and choose a new output directory:

```bash
bash workflow/09_environment_enrichment/09_00_run_enrichment.sh \
  /path/to/Figura4 /path/to/Supplementary_Data1/reference \
  local_runs/environment_enrichment

python3 workflow/08_figures_tables/08_03_read_based_figure4.py \
  --figure-data /path/to/Supplementary_Data1/figure_data \
  --output local_runs/current_figures

# Optional comparison of complete output tables with the verified external package.
python3 checks/verify_environment_enrichment.py \
  --results local_runs/environment_enrichment \
  --reference-results /path/to/Supplementary_Data1/results
```

The launcher checks external inputs, runs primary and single-KO sensitivity
analyses, joins labels and gene membership, verifies the results and renders
Figure 5. Its PNG, PDF and SVG files stay in the local output directory. Figure 4
can also read a complete new gene-model run with `--analysis-root` instead of
`--figure-data`; that route additionally needs `selected_unique_counts.npy`.
Both renderers use the approved blue gradient. Neither refits models.

Expected checks are 42,770 KO/PFAM tests and 2,606 module tests per gene cutoff;
102 term–axis–sign contrasts (15 KO, 87 PFAM; 96 distinct terms); and 11 module
contrasts (10 distinct modules). Of the primary contrasts, 45 term contrasts and
one module contrast persist when genes are selected at q <=0.05. Excluding
multi-KO clusters retains 11 of 15 primary KO contrasts and no primary module
contrast. The sensitivity has its own complete, smaller eligible background.

## Reconstruct the read-based gene models

Full preparation requires the reassembled read-derived export, supplied outside
GitHub. The input directory must contain:

```text
analysis_metadata/metadata_003b.csv
samples/{sample_id}/historical_gene_counts.tsv.gz
samples/{sample_id}/COMPLETE.json
02_annotations/historical_gene_dictionary.tsv.gz
```

Each gene table has `gene_numeric_id`, `unique_total`, `RPK` and `both_RPK`.
`COMPLETE.json` supplies the independently exported sample totals. The dictionary
maps numeric IDs to historical cluster IDs and KO/PFAM sets. This release starts
from those exports; it does not rerun FASTQ processing, mapping or annotation.
Allow several GB of memory and sufficient disk space for sample caches and
the selected matrices. The full-catalogue pass is much more expensive than ORA.

```bash
bash workflow/06_gene_associations/06_05_run_read_based.sh \
  /path/to/reassembled_gene_export local_runs/read_based_genes
```

The script refuses an existing analysis directory. It retains genes detected in
at least five samples with at least ten unique fragments in total. `ln(RPK+1)`
is centred over all 12,540,513 retained genes before the top 100,000 CLR variances
are selected. All 51 samples and 17 profiles are retained. Common and
depth-interaction models use locality, categorical depth and a profile random
intercept. The predictor and response are standardized; REML uses analytic
Satterthwaite derivatives. See `METHODS.md` for adjustment families.

The model status deliberately remains validation-pending after a new full fit.
Do not relabel it as validated merely because files were produced. The enrichment
launcher pins the archived, verified result tables byte for byte; newly refitted
or recompressed tables require an independent comparison before changing their
manifest. The individual enrichment programs expose explicit `--source`,
`--output` and, where needed, `--reference` arguments for reviewed new runs.

## Numerical validation and release checks

The external analysis archive contains 200 standardized gene responses and a
prespecified selection of 357 gene–axis pairs. Both model types give 714 fits and
1,785 contrasts. The following refits this subset in Python and compares it with
the archived independent R results:

```bash
OPENBLAS_NUM_THREADS=1 python3 checks/verify_read_based_reml.py \
  --analysis-root /path/to/Figura4 --output local_runs/reml_validation
```

The comparison preserves the original handling of 20 fits close to the
zero-variance boundary. R derivatives evaluated in variance-ratio coordinates
replace the corresponding degrees of freedom and P values; other estimates are
unchanged. The two R benchmark producers are included in
`workflow/06_gene_associations/read_based/r_validation/`. They take the external
analysis directory as their first argument and require its saved selection and
response tables. They do not create or upload those inputs.

The code release was checked by rerunning the complete enrichment and both
single-KO sensitivities, comparing all six statistical output tables with the
external reference, refitting the 714-model validation subset, rendering both
figures locally, and checking Python/Bash syntax and the text-only commit list.
The full catalogue was not reconstructed and all 1.2 million gene models were
not refitted during this repository update.

## Interpretation and manuscript mapping

Enrichment uses Figure 4A common slopes only, separately for each environmental
axis and positive/negative sign. The primary gene cutoff is q <=0.10; q <=0.05
is a sensitivity. The universe is restricted to tested, usable, nonsingular
genes annotated in the relevant ontology. Zero-hit eligible terms stay in the
BH family. Minimum selected-gene requirements apply after BH adjustment.

The module analysis has a separate multiple-testing family and requires at least
two distinct KOs as well as minimum gene counts. Modules remain **potential
functional candidates**, including candidates sensitive to multiple-KO
annotation. The results do not establish complete pathways, expression, enzyme
activity or ecosystem process rates. Gene dependence can inflate nominal
enrichment evidence; correlated axes and shared genes are not independent
ecological mechanisms. K05889 annotation does not establish plastic pollution.

Figure 4 remains the gene-association figure. Figure 5 is a separate functional
interpretation figure, preserving readable panel sizes. Supplementary Table 19
summarizes the potential modules; full PFAM, gene membership, background and
sensitivity tables belong in external Supplementary Data 1. The manuscript and
supplement remain their existing online documents; no copies are added to GitHub.
