# Carbon, nitrogen and sulfur markers

Section `10_functional_cycles` produces the current manuscript section
**Carbon, nitrogen and sulfur markers vary with sediment depth**. It follows
environmental gene associations (Figure 4) and functional enrichment (Figure 5).
It uses the same corrected read-based KO export and the same 51 samples from
17 complete sediment profiles. It does not require network results or rerun
gene selection. All section 07 network code is now historical.

This is a targeted exploratory follow-up, added after the broad functional
analyses. The panel was fixed before its own empirical tests; it was not
preregistered before the study and is not validation with independent samples.
The code release adds no sample data, reference lists, annotation databases,
images, result tables, fitted objects or manuscript copies. Those files remain
external or are generated locally. No script accesses Landa Lab or uploads files.

## Connection to the current manuscript

| Manuscript item | Local producer/output |
|---|---|
| Marker Methods and Results | `analysis_specification.json`; all numerical tables under `results/` |
| Figure 6 | `figuras_cns.py` → `figures/Figura_CNS_profundidad.{png,pdf,svg}` |
| Supplementary Figure 10 | `figuras_cns.py` → `figures/Figura_CNS_sensibilidad.{png,pdf,svg}` |
| Supplementary Table 16 | `report_tables.py` → `results/Supplementary_Table16_CNS.tsv`; all 11 markers, restoration tests and both depth comparisons |
| Supplementary Table 17 | `report_tables.py` → `results/Supplementary_Table17_CNS.tsv`; all 19 environmental models |

The figures label both the cycle and candidate function. Five unresolved
fractional assignments carry a dagger. Intervals are pointwise; asterisks denote
the separate multiple-testing correction, not simultaneous interval coverage.
The top three Figure 6 examples were selected after reviewing the complete
panel; its lower panel retains all 11 markers. No new ecological Discussion
text is supplied in this release.

Current supplementary numbering is: functional tests, Table 6; core, Table 10;
potential enrichment modules, Table 15; C/N/S markers, Tables 16–17. The old
network Figure 10 and Tables 10–13 no longer belong to the current manuscript.
Historical network file names are preserved only for traceability.

## Marker definitions

| Cycle | Marker | KO | Candidate function | Allocation in the primary analysis |
|---|---|---|---|---|
| Carbon | rbcL | K01601 | CO2 incorporation by RuBisCO | Recoverable single KO |
| Carbon | acsB | K14138 | Acetyl-CoA synthesis | Provisional fractional |
| Nitrogen | nifH | K02588 | Nitrogen fixation | Provisional fractional |
| Nitrogen | nirK | K00368 | Nitrite to nitric oxide | Provisional fractional |
| Nitrogen | nirS | K15864 | Nitrite / hydroxylamine reduction | Recoverable single KO |
| Nitrogen | norB | K04561 | Nitric oxide to N2O | Recoverable single KO |
| Nitrogen | nosZ | K00376 | N2O to nitrogen gas | Provisional fractional |
| Nitrogen | nrfA | K03385 | Nitrite to ammonium | Recoverable single KO |
| Sulfur | dsrA | K11180 | Sulfite transformations | Recoverable single KO |
| Sulfur | aprA | K00394 | Activated sulfate transformations | Recoverable single KO |
| Sulfur | soxB | K17224 | Thiosulfate oxidation | Provisional fractional |

`code/marker_annotations.py` records KEGG entry links and the specific limits
of each assignment. A unique KO allocation is not sequence-level functional
validation. K15864 includes hydroxylamine reduction and is not exclusive to
denitrification; K01601 alone does not establish a complete Calvin cycle.
dsrA/aprA do not resolve net sulfur reduction versus oxidation. Supporting
subunits are descriptive concordance checks, not independent confirmations.
Sparse or absent targets, including methane markers, remain in the detection
report and were not converted into zero-filled hypothesis tests.

## Inputs and reproducible execution

Use Python 3.12 and `requirements-analysis.txt`. The manuscript run used Python
3.12.14, NumPy 2.3.5, pandas 2.2.3, SciPy 1.17.0 and Matplotlib 3.10.8. Neither R,
python-docx nor scikit-learn is required for this section. The source numerical
calculations and figure layout are preserved. Eight source scripts are copied
byte for byte; two change only binary-output writing to verify complete files
before accepting a run. Source and release hashes are recorded in
`FUNCTIONAL_CYCLES_PROVENANCE.json`. Execution-status wording in the configuration
was updated without changing markers, hypotheses, seeds or numerical settings.

`FUNCTIONAL_CYCLES_INPUTS.json` pins 15 required compact external files: 12 KO
quantification matrices, canonical metadata, an annotation-pattern subset and
the official 76-KO reference list. These compact matrices retain the panel,
supporting/descriptive targets and reference KOs. They must not be substituted
with the prevalence-filtered composition or 200-KO network matrices. The two
primary full source KO hashes match the corrected export used by sections 05
and 07; extra quantification layers support annotation and read-pair sensitivity.

Choose an unused output directory and run from the repository root:

```bash
python3 checks/check_functional_cycles_inputs.py --input-dir /path/to/CNS/inputs
python3 workflow/10_functional_cycles/10_00_run_cycles.py \
  --input-dir /path/to/CNS/inputs \
  --output-dir local_runs/functional_cycles
```

The launcher checks every required hash and the sample/profile design before
creating outputs. It stages only the declared inputs and exact code/configuration
in a new run directory, runs all five quantification variants, checks the
results, writes both essential supplementary tables and renders figures locally.
Binary resampling files and figure exports are written atomically and checked
against their complete in-memory contents; PNGs are decoded before acceptance.
The launcher stops on failure, rejects existing output directories and reports progress
every 30 seconds during long steps. Logs and an execution manifest are retained.
Python assertions must remain enabled; do not use `python -O`.

The 15 compact inputs can also be reconstructed from the corrected export:

```bash
python3 workflow/10_functional_cycles/10_01_prepare_inputs.py \
  --source /path/to/corrected/input \
  --reference-list /path/to/uscg_76_kegg_min_2313.lst \
  --output-dir external_inputs/CNS
```

The source directory must contain `matrices/`,
`03_unit_annotations/annotation_patterns.json` and `metadata_003b.csv`. The
preparer verifies all full-source hashes, extracts the compact rows and patterns,
then verifies the resulting package against the compact manifest. It does not
download the reference list or regenerate read mapping.

Optional flags on `10_00_run_cycles.py`:

- `--with-previous-comparison`: also check three external inputs and run the
  descriptive comparison with earlier analyses and total-KO normalization.
- `--reference-results /path/to/reviewed/results`: compare all generated
  scientific TSV tables with the external reviewed result directory. The two
  new supplementary table assemblies are checked against their source values.
- `--no-figures`: calculate and verify all scientific tables without rendering.
- `--calibrate`: repeat design simulations and interval diagnostics before
  fitting empirical outcomes. This preserves the rejected unrestricted
  bootstrap-t interval diagnostic; that interval is not used for reporting.

Existing results can be checked or redrawn without refitting:

```bash
python3 checks/verify_functional_cycles.py --self-test
python3 checks/verify_functional_cycles.py \
  --run-root local_runs/functional_cycles \
  --reference-results /path/to/reviewed/results
python3 local_runs/functional_cycles/code/figuras_cns.py
```

## Statistical specification

The primary quantity is marker RPK divided by the median RPK of 76 universal
single-copy reference KOs. It approximates genomic representation, not absolute
copies per sediment mass or complete genome-resolved pathways. It is not the
full MUSiCC inter- and intra-sample correction. Six markers have recoverable
single-KO numerators; five use fractional numerators and remain provisional.

Restoration models use each profile's mean log2 ratio over its three depths,
with equal horizon weights, fixed locality and five categorical conditions.
They have 17 observations, rank 7 and 10 residual degrees of freedom. A
null-restricted wild bootstrap uses 99,999 Rademacher draws, one independent
weight per profile, HC2-scaled reduced-model residuals, an HC3-studentized
four-coefficient Wald statistic and the plus-one P-value correction. It tests
an overall condition contrast; individual condition coefficients are descriptive.

Depth tests compare 20 and 40 cm with 5 cm within each of the 17 profiles.
Effects are median paired log2 differences, reported as ratios. Two-sided exact
sign tests discard exact zero differences and report the effective sample size.
The fifth and thirteenth ordered differences provide a 95.0958% pointwise
interval under independent observations with a continuous common median.
The 11 restoration tests and 22 depth tests share one 33-test Holm family.

The 19 declared environmental models compare profile means, with fixed locality
and one standardized predictor (13 residual degrees of freedom). Concentrations
are log2-transformed after averaging depths; redox is not log-transformed.
Restricted wild-bootstrap P values share one 19-test BH family; intervals use
HC3-t13. These models do not adjust for restoration condition and do not test
within-profile environmental responses. The ambiguous redox field name does
not establish its measurement depth. Pointwise intervals are not adjusted
confidence intervals for the Holm or BH families.

The five variants are primary, all-fractional markers, a 56-reference subset,
fragments supported by both reads, and broader annotation with fractional
allocation. Each uses the same 17 profiles and its own complete correction
families. They are sensitivity analyses, not five independent studies. The
broader-annotation variant should also be compared with the fractional column
to hold allocation treatment constant. Algebraic allocation bounds describe
label uncertainty, not confidence intervals or all annotation errors.

Leave-one-profile/locality checks remain descriptive. A random-intercept REML
crosscheck estimates mean depth changes, whereas the primary depth estimand is
a median. The Python REML implementation comes from the previously checked
gene-analysis code; a new independent R comparison was not made for this panel.
All inference assumes profiles are independent after the stated adjustment and
is conditional on the three sampled localities. Bootstrap count never enters
biological sample size or model degrees of freedom.

## Results continuity and validation scope

nirS, dsrA and rbcL have lower normalized representation at 40 than 5 cm across
all five specifications. Their primary median paired ratios are approximately
0.320, 0.500 and 0.366. nirK and soxB retain their provisional status and have
processing-dependent support. No overall restoration test or environmental
association passes its declared multiple-testing correction. Null results and
all five variants are retained.

This does not contradict earlier restoration or network analyses: the main
signals concern paired depth changes, the response is normalized to reference
genes, and the correction family is smaller. Earlier broad gene models did
detect environmental associations. The marker panel does not establish recovery,
equivalence between conditions, measured activity, or causal effects. Holm
controls the declared marker family, not all searches in the project.

Release validation reconstructs the compact inputs, runs the complete panel
and optional comparison, verifies 110 depth tests independently from sample IDs,
checks Holm/BH families and 17-profile resampling, compares complete generated
scientific tables with the reviewed outputs and renders the two figures locally.
Invalid input hashes, duplicate profile-depth rows and inconsistent profile IDs
must fail. `FUNCTIONAL_CYCLES_VALIDATION.json` records the executed checks.
Raw-read reconstruction, new statistical calibration, and independent sequence
validation are outside this repository update.
