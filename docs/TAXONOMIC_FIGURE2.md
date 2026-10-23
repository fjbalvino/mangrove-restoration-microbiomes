# Taxonomic change points and Figure 2

The current Figure 2 uses the fixed candidate set of 1,031 taxa, the canonical
51-sample phyloseq object and the environmental values in
`source_data/indices/metadata_003b.csv`. Functional matrices are not inputs.
All large inputs, fitted objects, result tables and figure files remain external.
The completed rerun and its numerical checks are summarised in
[TAXONOMIC_FIGURE2_VALIDATION.json](TAXONOMIC_FIGURE2_VALIDATION.json).

## Source scripts and corrections

The original Figure 2 combined the global forest plot from
`workflow/04_taxonomic_change_points/04_04_plot_global_TITAN.R` (historical 404)
and the depth forest plot from `04_05_depth_TITAN.R` (historical 405).
The global models are fitted by `04_01_global_TITAN.R` (historical 401).
Steps 402 and 403 subsequently extracted and annotated the global indicators.
The current assembled renderer is `04_07_assemble_figure2.py`.

Both fitting scripts now default to the canonical 003b metadata and check all
six predictors, sample IDs, locality, stage, depth and profile identifiers by
sample ID. They reject historical MHI values rather than silently accepting a
column merely named `MHI_local`. The five other environmental axes have the
same values in the historical and canonical tables.

The global extractor now reads `sppmax` and `sumz.cp` directly. The former call
to the unavailable `TITAN2::titan.out()` could turn an error into a one-row
indicator table. The global summary also reversed the numeric direction codes;
TITAN2 uses 1 for decreasing and 2 for increasing responses. Later 402/403
tables had already repaired these aspects of the archived results. Current
counts, thresholds and plot symbols all derive from the same fitted objects.
The depth summary distinguishes all tested taxa from retained indicators.
No viral taxa are excluded; the misleading historical exclusion counter and
the global taxon-list filename were corrected.

## Analysis design

- Global: sample Hellinger transformation within the complete 1,031-taxon
  candidate set, then the mean of the three depths for each of 17 profiles.
  Each environmental predictor is averaged over the same three samples.
- Separate depths: 17 observations per horizon, one per profile. Retained taxa
  are 1,015 at 5 cm, 978 at 20 cm and 997 at 40 cm, using occurrence in at least
  three samples, total count at least one and positive transformed variance.
  The Hellinger denominator is computed before these depth filters.
- TITAN2 2.4.4: `minSplt=5`, `numPerm=999`, `nBoot=500`, `imax=FALSE`,
  `ivTot=FALSE`, `pur.cut=rel.cut=0.95`, `ncpus=1`, `memory=FALSE`.
- Global RNG: seed 123, continuous stream in the historical order MHI,
  nutrients/sulfide, physicochemistry, NDMI, surface water, vegetation/landscape.
- Depth seeds: base 123 plus the fixed position in the historical 3-depth by
  6-gradient loop. Within each depth the order is MHI, NDMI, nutrients/sulfide,
  physicochemistry, vegetation/landscape and surface water. Seeds are 124–129,
  130–135 and 136–141. Requesting a subset does not change those assignments.

These are exploratory univariate associations. TITAN purity and reliability
criteria are not an across-taxon/gradient false-discovery-rate correction.
Locality and other gradients are not jointly adjusted. Reusing the same
profiles across horizons does not provide a formal test of depth differences.
Overlapping bootstrap ranges do not establish equality. The analysis does not
identify causal environmental effects, a measured hydroperiod, moisture stress
or irreversible ecosystem tipping points.

## Repeat the analysis

Install R 4.3.3 and TITAN2 2.4.4 with the dependencies used by the fitting scripts
(including phyloseq, dplyr, tidyr, readr, tibble, stringr, purrr, ggplot2 and scales).
Python plotting requires pandas, NumPy and Matplotlib. Keep the following data
outside this repository:

1. `001_phyloseq_canon_51.rds`.
2. A 202 output directory containing
   `tables/202_CLR_raw1031_muestras_por_taxa.csv`; only its fixed taxon IDs are
   used, while abundance values come from phyloseq.

Run from the repository root, substituting the external input paths:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla \
  workflow/04_taxonomic_change_points/04_01_global_TITAN.R \
  --ps /path/to/001_phyloseq_canon_51.rds --out202 /path/to/202_run \
  --out_root /path/to/results --skip_plots true

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla \
  workflow/04_taxonomic_change_points/04_05_depth_TITAN.R \
  --ps /path/to/001_phyloseq_canon_51.rds --out202 /path/to/202_run \
  --results_root /path/to/results --skip_plots true
```

The default uses the unmodified package calculations. Optional
`--stream_permutations true` bounds memory by processing one candidate
partition's permutation matrix at a time. It retains the native IndVal,
permutation, mean, standard-deviation and probability calculations. The
execution-equivalence check compares all 13 fitted-object components and the
final RNG state with the unmodified package, including actual 17-profile,
1,031-taxon inputs when supplied.

Optional `--checkpoint_root /path/to/checkpoints` saves completed bootstrap
replicates and RNG state without changing resampling. Interrupted runs can
resume only with identical ordered inputs, parameters, R/TITAN2 versions and
initial RNG state. A checkpoint is not a completed analysis.

For independent scheduling of global gradients, `prepare_global_rng.R` provides
`advance_titan_rng()` for the pinned TITAN2 version. The `--rng_state` option
accepts a record containing the gradient, initial and expected final RNG states,
RNG kind, R/TITAN2 versions and resampling parameters. Each completed fit must
match the expected final state. Standard sequential execution needs no such
record. Do not replace the continuous stream with a reset to seed 123 for every
global gradient.

## Extract and assemble

Create an external CSV manifest with columns `stratum`, `gradient`, `rds` and
one row for every fit. `stratum` is `Global`, `5`, `20` or `40`. The six gradient
names are the canonical column names. Paths in `rds` identify completed 401 or
405 TITAN objects. Relative paths are resolved from the manifest directory. Then run:

```bash
Rscript --vanilla workflow/04_taxonomic_change_points/04_06_extract_figure2.R \
  /path/to/run_manifest.csv /path/to/figure2_tables

python workflow/04_taxonomic_change_points/04_07_assemble_figure2.py \
  --data /path/to/figure2_tables/Figure2_source_data.csv \
  --output /path/to/figures/Figure2
```

The extractor verifies all 24 fits, the 500-replicate arrays, model arguments,
taxon denominators and agreement between stored filters and recomputed
purity/reliability. It suppresses filtered thresholds for directions with zero
retained indicators. Thin intervals are the 5th–95th percentiles (central 90%);
thick intervals are the 10th–90th percentiles (central 80%).

The assembled figure has six panels, one per gradient. Each uses a shared
horizontal scale for its global and three depth estimates. Its endpoint labels
follow the verified environmental loadings; scales are not comparable across
different gradients. PNG (600 dpi), vector PDF and editable SVG are generated
locally. Decorative artwork is not required by or included in this workflow.

Run execution checks with:

```bash
OPENBLAS_NUM_THREADS=1 Rscript --vanilla checks/check_titan_execution.R
```

An optional argument accepts an external RDS list with `env` and `taxa` to
extend the equivalence check to real analysis inputs. The checks also reject
obsolete MHI values and duplicate sample IDs while accepting reordered rows
that match by ID.


## Numerical precision across historical runtimes

All retained taxon identities, response directions, counts and community
thresholds/ranges reproduce the reference analysis (maximum threshold/range
difference <5e-15). Raw reliability values are not bitwise identical:
1436 of 24,126 taxon-gradient rows differ, by at most
0.020. Among retained indicators,
2 reliability values differ, by at most
0.002; no selection changes.

A controlled replay of the first five MHI bootstraps at 5 cm used the same
ordered inputs and seed 124. With reference BLAS it reproduced every archived
bootstrap metric exactly. With OpenBLAS, small changes in matrix products
altered permutation tie comparisons and some P values. This confirms numerical
backend sensitivity in the native package, independently of the optional
execution aids. The current rerun used OpenBLAS 0.3.26; the corrected archived
MHI analysis used reference BLAS. Historical five-axis runs used a different
OpenBLAS build. Full historical reliability differences remain in the external
audit; identical figure values do not imply bitwise identity of every metric.

To repeat that diagnostic using the archived MHI 5-cm fit supplied with the
external analysis package (the five diagnostic replicates are not a replacement
for the 500-replicate inferential runs):

```bash
OPENBLAS_NUM_THREADS=1 Rscript --vanilla checks/check_titan_blas_precision.R \
  /path/to/reference_MHI5_titan.rds 124 openblas /path/to/audit/openblas

OPENBLAS_NUM_THREADS=1 \
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 \
Rscript --vanilla checks/check_titan_blas_precision.R \
  /path/to/reference_MHI5_titan.rds 124 reference /path/to/audit/reference
```

The second command requires the reference BLAS library at the shown Linux
path; it selects that library for that process only.
