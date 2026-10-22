# Functional composition and proximity to conserved references

This is the corrected read-based analysis for Figure 3, Supplementary Tables 6
and 14, and Supplementary Figure 9 in the current online supplement. It replaces
the earlier section 05 implementation at the same code paths. Generated images,
PDFs, fitted objects and full result archives are not committed.

## Inputs and measurement

The four external text matrices required in `source_data/functional/read_based/`
contain 51 sample columns, with functions in rows: 8,416 KOs and 8,034 PFAM features.
This code update does not publish these matrices. Their filenames, byte sizes and
SHA-256 checksums are in `FUNCTIONAL_INPUTS_REQUIRED.tsv`. Sample identifiers, site, locality, depth
and stage come from `source_data/indices/metadata_003b.csv`. `MHI_local` is retained
in that file but is not a predictor or selection variable in this section.

Counts represent mapped paired-end fragments. Each fragment is counted once
when its reported accepted matches identify a single reference unit; either or
both reads can provide support. Ambiguous fragments are excluded. KO and PFAM
are aggregated separately using complete annotation agreement among linked
protein clusters. When an agreed annotation contains multiple terms, its count
is divided equally among them.

Before functional aggregation, each reference-unit count is divided by sequence
length in kilobases: `assigned fragments * 1000 / nucleotide length`. Output files
call this RPK, although the numerator here is fragments, not individual reads.
No genomic-equivalent correction is applied. `QUANTIFICATION.json` records the
upstream assignment rules and aggregate audit provenance. The separate source
export contains `annotation_coverage_by_sample.tsv` for sample-level coverage;
that audit table is not bundled or required by this analysis launcher. These inputs describe DNA-based functional
potential, not expression or process rates. They do not report assembly quality.

The old `082A_*` matrices count assembled ORF cluster members. They remain only
because the unrecalculated 808 branch consumes them. The fixed network input
also remains historical. Neither is used by the corrected primary launcher.
The corrected gene associations and environmental enrichment are documented in
`GENE_ENVIRONMENT.md`; corrected network execution is documented in `NETWORKS.md`.
Ordered restoration trends remain separate.
The section 05 launcher does not rerun read mapping or those other blocks.

## Statistical design

- Retain functions occurring in at least 26 samples, with summed fractional
  fragment count at least 10 and positive variance after log(1 + length-adjusted
  count). Counts determine the total-count filter; length-adjusted values supply
  composition. This retains 7,802 KOs and 7,190 PFAM features.
- Add one to length-adjusted values, take natural logarithms, and subtract the
  sample mean within each retained annotation layer. Euclidean distances between
  these CLR profiles are Aitchison distances.
- Test stage after locality and depth with partial dbRDA/marginal PERMANOVA.
  Use 9,999 distinct nonidentity permutations of complete three-depth profiles
  within locality, preserving depth order (seed 123). There are 17 independent
  sites. The fitted stage and residual degrees of freedom are 4 and 42; the
  51 rows are not treated as 51 independently permutable sites. Dispersion uses
  stage centroids, bias adjustment and the same permutation matrix.
- Compare each non-conserved sample with Conserved from the same locality and
  depth. Average across depths within each of 14 non-conserved sites, then test
  stage after locality (F degrees of freedom 3, 8; 9,999 within-locality
  permutations; seeds 125 and 126). References are fixed, one profile per locality.
- Figure 3c,d uses `distance(sample, Conserved) - distance(Degraded, Conserved)`.
  Negative values mean closer reference proximity. These depth-specific contrasts
  are descriptive; most locality-stage cells contain one site.
- Report all four preprocessing checks: added value 0.5, occurrence in at least
  five samples, Hellinger/Bray-Curtis, and rescaling each complete annotation layer
  to one million before adding one and calculating CLR. Rescaling changes the
  relative contribution of a fixed pseudocount. KO and PFAM are related annotation
  layers from the same samples, not independent biological replications.

Core detection uses at least 0.01% relative abundance in every one of the nine
Conserved samples. Relative abundance uses the complete annotation layer before
composition filtering. The core contains 1,647 KOs and 1,620 PFAM features;
Degraded median coverage is 99.4536% and 99.4444%. Coverage counts detected core
functions; it is not the fraction of all reads belonging to the core. Stage cores
use the same threshold in every sample of that stage, so unequal sample numbers
limit direct comparisons of core size.

Core CLR is recalculated within the core. Core proximity uses
`distance(Degraded, Conserved) - distance(sample, Conserved)`, averaged across
depths within site: positive values mean closer proximity. This is an absolute
distance difference, not a percentage, and its sign is opposite to Figure 3.
Stage tests use the 14 site means and locality-only residual permutations.

## Reproduce

Validated versions: R 4.3.3, vegan 2.6-4, permute 0.9-7, ggplot2 3.4.4 and dplyr
1.1.4. Python versions are in `requirements-analysis.txt`. Install these packages
in your chosen environment. Supply these four files from the verified functional
matrix export in `source_data/functional/read_based/`:

- `KO_complete_counts.tsv`
- `KO_complete_RPK.tsv`
- `PFAM_complete_counts.tsv`
- `PFAM_complete_RPK.tsv`

These files are ignored by Git. The launcher verifies their recorded sizes and
SHA-256 checksums before fitting models; missing or changed inputs stop execution.
Then run from the repository root:

```bash
set +e
bash workflow/05_functional_reference/05_00_run_primary_composition.sh
```

The launcher validates shipped and externally supplied input hashes, reports
progress every 30 seconds, and runs
602 composition/dispersion/reference tests, 603 ordination, 604 descriptive
contrasts, 605 core analyses, numerical validation and figure assembly. It stops
on an error. It refuses to replace an existing 602 output directory. For another
execution, set `FUNCTIONAL_OUT_ROOT` to a fresh local directory. Optional
`FUNCTIONAL_RSCRIPT` and `FUNCTIONAL_PYTHON` select installed interpreters.

Default outputs are under `local_runs/functional_reference/`, ignored by Git.
The `figures/` directory contains Figure 3 (185 × 230 mm, 8–11 pt) and the two-page
Supplementary Figure 9. To redraw from existing fitted tables without refitting:

```bash
python3 workflow/05_functional_reference/05_06_assemble_figures.py \
  --out_root local_runs/functional_reference
```

The validator compares primary and sensitivity tests, core statistics and all
84 Figure 3 contrasts with the reviewed numerical expectations. It also checks
permutation restrictions and reconstructs the 84 contrasts directly from the two
source matrices, independently of the plotting tables. Full tables, sample
scores, permutations and R session information are generated locally.

## Interpretation retained in the current text

Stage explains 12.7% of total KO variation and 13.4% of PFAM variation
(P = 0.0303 and 0.0255); dispersion P = 0.6428 and 0.5801. With rescaling before
the fixed pseudocount, stage P = 0.0558 and 0.0626. The stage association therefore
depends on preprocessing. Mean reference distance has no supported common stage
association (P = 0.9310 and 0.7320). Locality- and depth-specific contrasts remain
descriptive. These comparisons do not estimate temporal recovery or ecosystem
process rates.

Validation completed from the repository paths on 21 October 2026: all 200
numerical expectations and 84 independently reconstructed contrasts passed.
See `FUNCTIONAL_REFERENCE_VALIDATION.json` for the recorded scope.
