# Expanded biogeochemical analysis

This is the current section 10 workflow and Figure 6 producer. It covers carbon,
nitrogen, sulfur, methane/C1 and phosphorus. The preceding 11-marker panel and
its 33-test correction remain historical; do not combine their adjusted P values
with the expanded analysis. The expanded secondary question is whether paired
depth changes differ among stages, not whether an overall depth effect exists.

## Inputs and execution

Use the `reproducibilidad` directory from the external archive
`Panel_biogeoquimico_CF_51_completo.zip` (SHA-256 is recorded below). It contains
five RPK matrices, the annotation-pattern dictionary, canonical metadata,
reference list, frozen KO/map links, curated rules, previously audited per-sample
rule evidence, nitrogen comparison controls, and seven cached KEGG maps.
`workflow/10_functional_cycles/expanded/inputs.json` records all required hashes.
No source matrices, annotation databases, map images or complete result tables
are added to GitHub by this update.

```bash
python3 workflow/10_functional_cycles/10_01_run_expanded.py \
  --bundle-root /path/to/extracted/reproducibilidad --check-only
python3 workflow/10_functional_cycles/10_01_run_expanded.py \
  --bundle-root /path/to/extracted/reproducibilidad \
  --out-root local_runs/expanded_cycles
```

Python dependencies: numpy, pandas, scipy, matplotlib and Pillow. The launcher
checks hashes and the complete 51-sample/17-profile design, creates a fresh
`10_01_run_expanded_YYYYMMDD_HHMMSS` directory, links the immutable external
inputs, copies the preserved producer scripts, and runs them sequentially.
A failed producer stops the run. The inner statistical folder is named
`15_compare_expanded_YYYYMMDD_HHMMSS`; its prefix matches its producer.
`--prepare-only` verifies and stages the run without fitting models.

The producers preserve the executed analysis rather than silently recalculating
annotations against a current database. Cached maps are verified before launch;
the normal run needs no new KEGG download. KEGG assets remain attributed to
Kanehisa Laboratories. Consult KEGG terms for reuse.

## Methods and estimands

The manuscript wording is mirrored in
[EXPANDED_BIOGEOCHEMICAL_METHODS.md](EXPANDED_BIOGEOCHEMICAL_METHODS.md).
The full 8,416-KO matrix was queried against map01200, map00710, map00720,
map00680, map00910, map00920 and map00440, together with all operational rules
for the five cycles. This defines 748 unique KOs, of which 497 have source rows.
The primary positivity criterion retains 437 KOs: 180 recoverable single-KO
contributions and 257 provisional fractional estimates. Shared KOs enter the
statistical family only once. Phosphate transport and polyphosphate rules
supplement the phosphonate map; no map coverage fraction is called pathway
completeness.

Fractional abundance R and full-term abundance F constrain single-KO abundance.
For compatible multi-label annotation sizes s, the lower bound is
`max(0, F - (F-R)/(1-1/min(s)))` and the upper bound substitutes `max(s)`.
Without multi-label patterns, the original contribution is used. Coincident
bounds identify a recoverable single-KO contribution. Bounds are algebraic
assignment constraints, not confidence intervals or proof of biochemical
specificity. Annotation patterns and RPK are restricted to the matching
consensus status. Complete-consensus RPK is primary; supported partial consensus
is one sensitivity. Thresholds and numerical tolerances are preserved in code.

The primary numerator is divided by the original median RPK of 76 reference KOs.
No pseudocount is used. An unavailable KO remains unavailable, not zero.
Rules retain mandatory complexes (AND), alternative enzymes (OR), optional
components, and met/not-met/unknown states per sample. Community-level rules do
not establish genomic linkage. The producer reads the previously audited rule
layer, retaining its hashes; it does not reannotate sequences or rebuild that
layer from protein records.

For each KO, the primary response is the mean log2 ratio over the three depths
within each profile. The design contains an intercept, two locality indicators
and four categorical stage indicators relative to conserved (n=17; rank=7;
residual df=10). Secondary responses are paired log2 20/5 and 40/5 ratios, using
the same design. Null-restricted wild bootstrap uses 99,999 Rademacher vectors
(seed 2026111701), HC2-scaled locality-only residuals and a four-df HC3 Wald
statistic. Each vector weights whole profiles. P=(exceedances+1)/100000.
Holm correction jointly covers 437 × 3 = 1,311 tests in each variant.

Four alternatives change fractional allocation, the reference denominator to
56 recoverable reference KOs, both-mate support, or supported partial consensus.
They keep the primary eligible KO universe. An unestimable sensitivity test
receives P=1 and is marked explicitly. Pairwise coefficients are descriptive;
no pairwise tests or uncalibrated pairwise confidence intervals are supplied.
Leave-profile/locality-out coefficients are reported only when estimable.

## Outputs and validation

The run produces the updated Figure 6 (PNG/PDF), seven-map HTML atlas, within-stage
depth plots, complete tests and estimates for Supplementary Table 16, per-sample
rule states, diagnostics, omission sensitivities, source hashes and an explanatory
report. Figures use a common 1/4–4 stage/conserved ratio scale; saturation is
explicit. Missing and ineligible entries are distinct from observed zero.

All 1,311 primary adjusted P values equal 1; the minimum unadjusted stage P is
0.00818. No stage omnibus test passes Holm correction in any of the five variants.
K01659 has a depth40-by-stage signal only in the fractional variant
(Holm P=0.01311 versus primary P=1); this is not a robust restoration finding.
No support does not demonstrate equivalence, especially with 3–5 profiles per
stage. DNA establishes genetic potential, not activity or recovery of fluxes.

Independent least-squares calculations reproduce all primary stage coefficients.
The earlier nitrogen raw P values are unchanged; expanded multiplicity is
recomputed. The 84 atlas map/stage/contrast combinations generated valid SVG
under a JavaScript DOM stub. Interactive browser behavior was not verified
because Chromium could not launch in the execution environment. Manuscript
Figure 6 and surrounding text were checked in the exported PDF.
There was no new finite-sample calibration of the expanded bootstrap; inference
remains exploratory. See [validation](EXPANDED_BIOGEOCHEMICAL_VALIDATION.json).

External archive SHA-256: `2c3715162de20de9ec49eafe8cc45af21b208858322b2c422eb061f0de4f54f8`.
