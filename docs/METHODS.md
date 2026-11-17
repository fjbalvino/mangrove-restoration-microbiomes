# Methods represented by the repository

The map was established from the September 2026 manuscript. Sections 05, 06 and
09 now document the corrected October read-based analyses and environmental
enrichment. `CODE_PROVENANCE.tsv` records current code hashes and changes.

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

The historical HI/RBMA1 is the fixed first principal component of taxonomic CLR residuals after
categorical locality/depth adjustment, using 1,467 historical taxa and centred,
unscaled PCA. Its sign places the Conserved median above Degraded. Zero is an
ordination origin. The current 1,031-taxon screen does not redefine HI.
`checks/verify_historical_hi.py` reconstructs HI using QR residualisation and SVD.
HI is retained for archived network reproduction and is not used by the current
expanded biogeochemical section or described as a current manuscript index.

## Metagenomic inputs

The manuscript describes read QC, taxonomic classification, assembly, ORF
prediction, clustering and annotation. The received capsule begins downstream
of those processes; it does not establish all original raw-read commands or
database snapshots. Section 02 inventories, filters and aggregates the existing
annotated gene-abundance data. Current 502/504 scripts are maintained successors;
saved 044/082A matrices retain their original provenance; the minimal publication
keeps the current source implementations and documents the source hashes.
Recovered 082A values count assembled ORF cluster members and remain in
`source_data/functional/tables` only for historical reproduction. Corrected
read-based counts and length-adjusted values must be supplied separately in
`source_data/functional/read_based`; section 05 uses only these corrected inputs.
Its quantification rules are documented in `FUNCTIONAL_REFERENCE.md` and
`source_data/functional/read_based/QUANTIFICATION.json`.

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

The corrected 602/603/604 analyses retain 7,802 KOs and 7,190 PFAM features.
Filtering requires occurrence in at least 26 samples, summed fractional fragment
count at least 10 and positive log(1 + length-adjusted value) variance. CLR uses
length-adjusted values plus one separately within each layer. The models account
for locality and depth and use 9,999 distinct whole-profile permutations within
locality. Dispersion is tested separately with the same permutation matrix.

Reference tests use three-depth means from 14 non-conserved sites, adjusting for
locality (stage/residual degrees of freedom 3, 8). Each reference is Conserved
from the same locality and depth. Figure 3 uses sample distance minus the local
Degraded distance, so negative values mean closer proximity. Contrasts are
descriptive and do not measure temporal or ecosystem-process recovery.

All four preprocessing checks are reported. Stage significance depends on
preprocessing: after rescaling the complete layer to one million before adding
one and taking CLR, P values are 0.0558 and 0.0626. A normalization-independent
stage effect is not established.

The reference core uses relative length-adjusted values at least 0.01% in all
nine Conserved samples, with the complete layer as denominator. There are 1,647
KO and 1,620 PFAM core features; Degraded median coverage is 99.4536% and 99.4444%.
Core distance reductions use Degraded distance minus sample distance, the opposite
sign to Figure 3. These are distance differences, not percentages. The new
`05_04_reference_core.R` reproduces core coverage, overlap and site-level tests.
Full methods and commands are in `FUNCTIONAL_REFERENCE.md`.

The historical, unrecalculated 808 sensitivity uses 7,129 prevalence-filtered KOs from assembled-ORF counts. It fits stage,
depth and locality to 42 non-conserved samples, with 9,999 Freedman–Lane profile
permutations and 4,999 profile bootstraps. Its five-test Holm family and nine
exploratory depth-contrast BH family are specified in the bundled protocol.
It is not the producer of Figure 3. One reference profile per locality leaves
reference uncertainty unestimated.

## Gene-level associations

The current section 06 uses uniquely assigned fragments and length-adjusted
abundance (RPK). Retention requires detection in at least five samples and at
least ten fragments across the cohort. Of 22,456,116 detected clusters, 12,540,513
are retained. Natural-log CLR with pseudocount one is centred over this entire
retained catalogue before selecting the 100,000 genes with greatest CLR variance.
Selection does not use stage or any predictor.

Gaussian REML models adjust for locality and categorical depth, with a profile
random intercept. Each predictor is tested separately, with a common slope and
predictor-by-depth interaction. Slopes are expressed per predictor SD and response
SD. Satterthwaite t tests assess slopes; a two-degree-of-freedom F test assesses
interaction. The five environmental axes share separate BH families of 500,000
common slopes, 500,000 interactions and 1,500,000 depth slopes. MHI uses separate
families of 100,000, 100,000 and 300,000 tests. Singular/unusable fits receive
P=1 within adjustment and are excluded from primary discoveries.

At q <=0.10, the current results contain 42,222 common associations, 63,942
depth-specific associations and 1,234 interactions, involving 48,580 distinct
genes. MHI has no primary discoveries. Figure 4 uses these statistics and the
K00619 and K03621 examples. The older top-level `06_01`–`06_04` programs and
figure bundle remain historical; `06_05_run_read_based.sh` is the current entry.

## Functional enrichment by environmental axis

Section 09 tests common-slope genes from Figure 4A separately by axis and slope
sign. The background contains the tested, usable, nonsingular genes with an
annotation in the relevant ontology. Terms with at least ten background genes
enter one-sided Fisher tests, including terms with zero selected genes. BH is
applied to all 42,770 KO/PFAM tests across five axes and two signs, separately
for gene q cutoffs of 0.10 and 0.05. Terms are reported at enrichment q <=0.10
with at least five selected genes, after adjustment. Fold enrichment compares
the fraction of selected genes assigned to a term with its background fraction;
it is not an abundance fold change.

KEGG modules form a separate exploratory family of 2,606 tests per gene cutoff.
Eligibility requires ten background genes and two distinct KOs; reporting also
requires five selected genes and two selected KOs. Genes assigned multiple KOs
are excluded from both the selected set and background in an annotation sensitivity.
There are 102 KO/PFAM contrasts (96 terms) and 11 module contrasts (10 modules).
None of the module contrasts survives the single-KO sensitivity; all remain
potential functional candidates. Figure 5 presents KO terms and these candidates
separately from Figure 4. Full PFAM results remain external supplementary data.

Stored axis signs define direction. Descriptive correlations help interpret these
signs; they are not refitted PCA loadings. Shared genes and correlated axes do not
represent independent processes. Gene dependence can make nominal enrichment
FDR optimistic. DNA abundance supports functional potential, not pathway
completeness, expression, activity or process rates. See `GENE_ENVIRONMENT.md`.


## Environmental module component evidence

The original gene models and ORA remain distinct from the new descriptive component assessment. The current Figure 5a retains all 15 primary KO contrasts. Figure 5b retains all 11 candidate module contrasts and displays their primary and single-KO enrichment q values beside block coverage and per-sample representation. The original module bubble plot remains in the external archive.

The primary ORA was rerun from the pinned inputs. All six complete statistical tables matched the archived reference at the declared numerical tolerances. The result remains 102 KO/PFAM contrasts (96 terms) and 11 module contrasts (10 modules). Ten of those module contrasts remain eligible but nonsignificant in the single-KO sensitivity; M00417 becomes ineligible. The original thresholds and multiplicity families were not relaxed.

The new analysis evaluates the full required/alternative components of the 10 modules, selected after the original enrichment. It uses the complete CF KO matrix for representation and the tested 100,000 genes for associations. Its definitions were retrieved on 16 November 2026 (17 November UTC), separately from the October enrichment reference. All ten definition KO sets match the frozen October KO–module links. Existing annotations were not reassigned.


The evaluation preserves gene associations, module enrichment and community representation as separate evidence levels. It adds no block-level hypothesis tests. Full details and commands: [ENVIRONMENT_COMPONENTS.md](ENVIRONMENT_COMPONENTS.md).

## Expanded biogeochemical pathway and component analysis

We examined carbon, nitrogen, sulfur, methane and phosphorus components using the full matrix of 8,416 KOs. Seven KEGG pathway maps and curated rules defined a panel of 748 unique KOs. The rules specified required subunits and alternative components. We evaluated them separately in each sample as supported, unsupported by the available evidence, or unresolved. Map identifiers and rule definitions are provided with Supplementary Table 16. The panel was fixed before its own tests, after the nitrogen pilot, and remains exploratory. We used the existing annotations without reannotating sequences.

We retained the original fragment-based RPK values and normalised them by the sample median across 76 reference KOs 42. We used recoverable single-KO contributions where these could be determined, and otherwise retained fractional contributions as provisional. The 437 KOs with positive primary values in all 51 samples were log₂-transformed without adding a pseudocount. Other components remained in the descriptive inventory. These ratios approximate genomic representation; they do not measure absolute abundance or activity. Components detected in a community need not occur in the same organism.

Restoration stage was the primary comparison. We averaged log₂-normalised values across the three depths of each profile and fitted models with locality and categorical stage, using conserved profiles as the reference. The 17 profiles were the independent units. We then analysed paired 20/5 and 40/5 cm changes to ask whether depth patterns differed among stages. Overall stage tests used a null-restricted wild bootstrap with 99,999 draws and heteroscedasticity-robust statistics 43. Holm correction 44 covered 1,311 tests: 437 stage comparisons and 874 comparisons of depth changes among stages. Shared KOs were counted once. Stage-to-conserved ratios describe estimated differences; they were not tested separately.

We checked sensitivity to annotation allocation, a denominator of 56 reference KOs with recoverable single-KO contributions, inclusion of only fragments supported by both mates, and inclusion of supported partial-consensus annotations. We also examined residuals and influential profiles, and repeated stage estimates after omitting each profile or locality where the comparison remained estimable. These omission ranges describe sensitivity, not confidence intervals. Bootstrap construction, diagnostics and all quantification alternatives accompany Supplementary Table 16. With only three to five profiles per stage, a lack of statistical support does not establish equivalence.


Technical details, input hashes and execution: [EXPANDED_BIOGEOCHEMICAL.md](EXPANDED_BIOGEOCHEMICAL.md). The former 11-marker analysis remains historical.

## Nitrogen paired-depth follow-up

For the exploratory nitrogen depth analysis, we retained the 49 eligible nitrogen KOs from the expanded panel and compared 20 and 40 cm with 5 cm within each profile. Exact two-sided sign tests assessed the direction of paired log₂ changes across the 17 profiles. We report median ratios and pointwise 95.1% order-statistic intervals under a common-median model; these intervals are not simultaneous, and heterogeneity among profiles limits their interpretation. Holm correction covered the 98 nitrogen depth tests. We also corrected them jointly with the preceding 1,311 tests, giving a combined family of 1,409. This does not cover every analysis explored in the study. We repeated the comparisons with all five quantification variants and retained the same eligible KO set. Within-stage depth summaries were descriptive because each stage contained only three to five profiles. Full results and sensitivity checks accompany Supplementary Table 16.

See [NITROGEN_DEPTH.md](NITROGEN_DEPTH.md) for execution, multiplicity and interpretation.

## Historical network methods

All network analyses have been withdrawn from the current manuscript. Their
source code and detailed methods remain in `workflow/07_networks/` and
`NETWORKS.md` for historical reproduction. They do not feed section 10.

## Method references

- Aitchison (1982), compositional analysis: https://doi.org/10.1111/j.2517-6161.1982.tb01195.x
- McMurdie and Holmes (2013), phyloseq: https://doi.org/10.1371/journal.pone.0061217
- Baker and King (2010), TITAN: https://doi.org/10.1111/j.2041-210X.2009.00007.x
- Bates et al. (2015), lme4: https://doi.org/10.18637/jss.v067.i01
- Kuznetsova et al. (2017), lmerTest: https://doi.org/10.18637/jss.v082.i13
- Benjamini and Hochberg (1995), FDR: https://doi.org/10.1111/j.2517-6161.1995.tb02031.x
- Manor and Borenstein (2015), MUSiCC reference normalization: https://doi.org/10.1186/s13059-015-0610-8
- Davidson and Flachaire (2008), wild bootstrap: https://doi.org/10.1016/j.jeconom.2008.08.003

References identify methods, not validation of study-specific indices.


