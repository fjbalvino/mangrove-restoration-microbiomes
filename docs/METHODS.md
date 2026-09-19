# Methods represented by the repository

This map follows the final manuscript and supplementary PDFs supplied on
18 September 2026 (`SOURCE_DOCUMENTS.tsv`). `METHODS_CODE_CROSSWALK.tsv`
identifies renamed scripts, original hashes and the scope of review.

## Study design and indices

The analytical cohort comprises 51 samples in 17 complete profiles from Carmen,
Cozumel and Tuxpan. Nominal depths 5, 20 and 40 cm represent sediment intervals
5–15, 20–30 and 40–50 cm. Profiles are the independent units; layers remain
together in the specified profile resampling schemes.

Environmental axes use median imputation, mean/SD standardisation and PCA.
The global ECI uses 19 variables, excluding depth and landscape Shannon
diversity. Five blocks describe vegetation/landscape, surface-water signal,
NDMI, physicochemistry and nutrients/sulfide. Stored PC signs define direction.

MHI uses 20 environmental variables, including landscape Shannon diversity,
excluding depth, microbial data and read counts. Within each locality, variables
are centred by their medians and divided by their MADs. Euclidean distance is
measured to the mean of local Conserved observations, pooling depths. The fixed
57-observation calibration gives
`MHI = 1 - 2*(distance - locality_min)/(locality_max - locality_min)`.
Values are imported by sample ID into the 51-sample cohort. The independent
implementation in `checks/verify_environmental_mhi.py` reproduces archived raw
and final scores. +1 denotes the smallest observed distance in that locality;
zero is the midpoint of its distance range. They are not ecological thresholds
or a common absolute scale across localities.

HI/RBMA1 is the fixed first principal component of taxonomic CLR residuals after
categorical locality/depth adjustment, using 1,467 historical taxa and centred,
unscaled PCA. Its sign places the Conserved median above Degraded. Zero is an
ordination origin. The current 1,031-taxon screen does not redefine HI.
`checks/verify_historical_hi.py` reconstructs HI using QR residualisation and SVD.

## Metagenomic inputs

The manuscript describes read QC, taxonomic classification, assembly, ORF
prediction, clustering and annotation. The received capsule begins downstream
of those processes; it does not establish all original raw-read commands or
database snapshots. Section 02 inventories, filters and aggregates the existing
annotated gene-abundance data. Current 502/504 scripts are maintained successors;
saved 044/082A matrices retain their original provenance; the minimal publication
keeps the current source implementations and documents the source hashes.
Recovered 082A counts are included in `source_data/functional/tables`.

## Taxonomic structure and change points

Section 03 contains Hill diversity, CLR/Aitchison composition, the Sarle/dip
screen and depth-specific variance partitioning. Diversity tests use profile
means. Composition tests preserve profile layers under locality-restricted
permutations. The joint raw-CLR screen selects 1,031 taxa using Sarle >0.555
and BH-adjusted dip P <=0.05. Other selections are historical/sensitivity
branches. Variation-partition plots normalise positive unique fractions,
excluding shared and unexplained variation.

Section 04 implements global TITAN2 on 17 profile means of Hellinger-transformed
abundances and depth-specific TITAN2 on 17 samples per horizon. Settings are
999 permutations, 500 bootstraps, minimum split 5 and purity/reliability >=0.95.
Means are taken after transformation; the 1,031-taxon denominator is preserved.
Gradients are tested separately. The predictor discrepancies in
the historical predictor-lineage section of `REPRODUCIBILITY.md` affect archived MHI thresholds and ECI variance fractions.

## Functional composition and reference proximity

Primary 602/603/604 analyses retain 7,126 KOs and 6,681 PFAMs with prevalence
>=26 samples, total abundance >=10 and positive log(1+abundance) variance.
CLR is calculated within each retained layer with pseudocount 1. Models condition
on locality and depth, preserve profiles in resampling and test dispersion
separately. Reference distances use the complete feature space and a conserved
sample matched by locality and depth. Conserved self-distances are excluded from
inference. Negative restored-minus-Degraded contrasts indicate greater proximity
to that reference, not measured ecosystem-process recovery.

The recovered 808 sensitivity uses 7,129 prevalence-filtered KOs. It fits stage,
depth and locality to 42 non-conserved samples, with 9,999 Freedman–Lane profile
permutations and 4,999 profile bootstraps. Its five-test Holm family and nine
exploratory depth-contrast BH family are specified in the bundled protocol.
It is not the producer of Figure 3. One reference profile per locality leaves
reference uncertainty unestimated.

## Gene-level associations

Section 06 centres CLR against the full filtered catalogue of 1,475,487 genes
before selecting 100,000 variable genes. REML mixed models use locality,
categorical depth and a profile random intercept, comparing common slopes and
predictor-by-depth interactions. Satterthwaite tests and BH adjustment use
declared families: five environmental axes pooled separately for common slopes,
interactions and depth slopes; MHI remains separate. Primary discoveries exclude
unusable/singular focal fits. Annotation supplies putative functions, not evidence
of expression or activity. Figure 4 reproduces saved statistics without refitting.

## Functional networks

Final analysis 809 holds 200 KOs fixed. Three constructions compare original CLR,
locality/depth residuals, and residuals additionally adjusted for four environmental
axes. NDMI is not added to that four-axis design. Nuisance fits and standardisation
are repeated within leave-one-out and full bootstrap reconstructions. Precision
matrices invert `0.9*S + (0.1 + 1e-8)*I`: regularisation is fixed, not estimated.

LIONESS sample weights are `n*w_full - (n-1)*w_leave_one_out`. Each sample retains
995 of 19,900 edges by absolute weight. Mean weight and weighted natural
connectivity describe weights, not changing edge counts or measured resilience.
Primary tests use 99,999 within-locality profile permutations; fixed-network
intervals use 4,999 profile draws, with 499 additional full reconstructions.
Joint correction families, consensus masks and global tests are specified in
`workflow/07_networks/network809/PROTOCOL.md`.

All three constructions and corrected results are retained. Network associations
do not identify confirmed individual KO drivers. HFR has insufficient held-out
coverage and remains an excluded ancillary branch in `workflow/supplementary/hfr`.

## Method references

- Aitchison (1982), compositional analysis: https://doi.org/10.1111/j.2517-6161.1982.tb01195.x
- McMurdie and Holmes (2013), phyloseq: https://doi.org/10.1371/journal.pone.0061217
- Baker and King (2010), TITAN: https://doi.org/10.1111/j.2041-210X.2009.00007.x
- Bates et al. (2015), lme4: https://doi.org/10.18637/jss.v067.i01
- Kuznetsova et al. (2017), lmerTest: https://doi.org/10.18637/jss.v082.i13
- Benjamini and Hochberg (1995), FDR: https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Kuijjer et al. (2019), LIONESS: https://doi.org/10.1016/j.isci.2019.03.021

References identify methods, not validation of study-specific indices.
