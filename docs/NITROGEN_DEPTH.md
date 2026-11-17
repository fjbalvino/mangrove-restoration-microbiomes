# Nitrogen paired-depth pilot

This exploratory extension adds overall paired depth comparisons to section 10.
It answers whether nitrogen components differ between 20 or 40 cm and 5 cm
across the 17 profiles. It is distinct from the existing question of whether
those vertical changes differ among restoration stages. Restoration remains
the manuscript's primary comparison. This update changes Results and Methods;
Discussion is unchanged.

## Run

Extract the external `Piloto_nitrogeno_profundidad_CF_51.zip` and provide its
root directory. The archive includes the five 49-KO normalised matrices,
metadata, frozen scope and eligibility, all preceding tests, and cached KEGG
nitrogen map assets. No FASTQ processing, remapping or annotation download is
needed. The input manifest records exact filenames, sizes and hashes.

```bash
python3 workflow/10_functional_cycles/10_02_run_nitrogen_depth.py \
  --pilot-dir /path/to/extracted/pilot --check-only
python3 workflow/10_functional_cycles/10_02_run_nitrogen_depth.py \
  --pilot-dir /path/to/extracted/pilot \
  --out-root local_runs/nitrogen_depth
```

Dependencies: Python, numpy, pandas, scipy, matplotlib and Pillow. The launcher
creates a fresh `10_02_run_nitrogen_depth_YYYYMMDD_HHMMSS` directory, verifies
inputs before execution and records code hashes. The preserved producer writes
`19_nitrogen_depth_YYYYMMDD_HHMMSS`, matching its script number. Failure stops
the run. Figures, complete tables, input matrices and KEGG assets stay external
to the repository.

## Design and quantification

Use 51 samples in 17 profiles with depths 5, 20 and 40 cm. The 49 eligible nitrogen
KOs were fixed from the expanded panel before their depth tests; they are the
union of the nitrogen map and operational nitrogen rules after primary
positivity filtering. Preserve the primary single-KO/fractional numerator and
median of 76 reference KOs. Reuse the four existing quantification alternatives.
No pseudocount is introduced. Unestimable sensitivity tests receive P=1 and are
marked explicitly, retaining the same family size.

Each profile contributes one log2 difference per depth contrast. Exact two-sided
sign tests assess directional consistency. Differences within 1e-12 of zero are
treated as ties and omitted from the sign test. Median log2 differences are
exponentiated to paired abundance ratios. Pointwise order-statistic intervals
have 95.0958% nominal coverage for 17 independent profiles under a common-median
model; they are not simultaneous intervals. Locality/stage heterogeneity limits
that common-median interpretation. Pairing accounts for repeated horizons, but
the test is not a separate adjusted effect for each locality or stage.

Holm98 includes 49 KOs × two depth contrasts. Holm1409 combines these 98 tests
with the preceding 1,311 stage and stage-by-depth-change tests, separately for
each quantification variant. Neither is a study-wide correction. Extending
paired depth tests to other cycles requires recalculating the expanded family.
Within-stage summaries remain descriptive (3–5 profiles per stage). Omission
ranges describe influence, not confidence intervals.

## Findings and limits

The primary Holm98 analysis supports declines in one KO at 20 cm and 12 at 40 cm.
Ten 40 cm signals retain support in all five variants. Under Holm1409, three
40 cm contrasts remain supported in the primary analysis and two remain
supported in all variants:

| KO | Component | Median 40/5 | Pointwise 95.1% interval | Primary Holm1409 P | Quantification |
|---|---|---:|---|---:|---|
| K15864 | nirS; nitrite/hydroxylamine transformations | 0.320 | 0.095–0.468 | 0.0215 | Recoverable single-KO |
| K00284 | gltS; ferredoxin-dependent glutamate synthase | 0.559 | 0.439–0.759 | 0.0215 | Provisional fractional |

Both decline in all 17 profiles in the primary analysis. K00262 is the third
primary Holm1409 signal, but does not retain support in all alternatives.
These are changes in annotated genetic representation, not evidence of changed
rates, complete pathways, or restoration recovery. In particular, nirS is not
exclusive to denitrification, and gltS retains allocation uncertainty.
No supported stage-by-depth-change interaction was established in the primary
expanded analysis; this does not negate a pooled paired depth effect.

## Outputs and verification

`depth_tests.tsv` contains all 490 variant–KO–depth tests, both multiplicity
corrections, paired median ratios and intervals. `paired_profiles.tsv`,
`within_stage.tsv`, `leave_out.tsv`, `summary.tsv` and `robustness.tsv` retain
individual profiles and sensitivity information. Figure and map producers use
the same numerical tables. The full plot shows all eligible KOs rather than
selecting by P. The HTML map shows pooled paired changes, not stage-specific
tests, and retains KEGG attribution. Browser interaction was not verified.

The delivered standalone producer was rerun and reproduced all statistical
columns exactly. Sign-test probabilities were checked independently with the
binomial distribution, family sizes were checked, and plots inspected.
See [validation](NITROGEN_DEPTH_VALIDATION.json) and the exact
[manuscript additions](NITROGEN_DEPTH_MANUSCRIPT.md).

External archive SHA-256: `2315de7d25fe7b69b64e3913b018c4eaf081e366e56ac1c19147dfdee7e1f225`.
