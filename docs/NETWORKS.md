# Historical corrected functional networks

This analysis has been withdrawn from the current manuscript and replaced by
the C/N/S marker section in [FUNCTIONAL_CYCLES.md](FUNCTIONAL_CYCLES.md).
Code, methods and audit records are retained for traceability. Figure numbers
below refer to the former network draft; current Supplementary Figure 10 is
the marker quantification sensitivity. These scripts are not prerequisites
for any current marker result.

`workflow/07_networks/07_02_run_read_based.sh` is the read-based historical entry point. It
prepares networks from corrected read-derived KO abundances. The older `07_00`,
`07_01` and `network809/` inputs and scripts remain historical reproductions.
No new input matrix, fitted object, statistical table or figure is published in
this code update. Outputs are generated locally; no script uploads them.

## Inputs and execution

Supply `KO_complete_RPK.tsv` and `KO_complete_counts.tsv` from the corrected
51-sample export. They contain 8,416 KOs; counts are fractional assigned fragments
and RPK is length-adjusted abundance. The metadata and HI/MHI calibrations are
unchanged. `prepare_inputs.py` verifies all three input SHA-256 values before
selection. Expected matrices are the same as section 05; see
`FUNCTIONAL_INPUTS_REQUIRED.tsv` and `FUNCTIONAL_REFERENCE.md`.

From the repository root, with `requirements-analysis.txt` installed:

```bash
bash workflow/07_networks/07_02_run_read_based.sh \
  /path/to/corrected/matrices \
  source_data/indices/metadata_003b.csv \
  local_runs/read_based_networks
```

Choose a new output directory for each execution. The launcher refuses an existing
directory, reports progress every 30 seconds, stops on errors and runs preparation,
primary inference, parameter sensitivity, profile deletion, algorithm comparison,
numerical feasibility, independent validation and plotting. `PYTHON` can select
the Python interpreter. The run takes minutes
rather than requiring a new mapping job. Inspect logs for the actual duration.

The `runtime/selection/` directory records feature eligibility, ranking and
selection stability. `results/LATEST_809.txt` identifies the complete current
run. The `809_` output prefix preserves the audited table schema; it does not
mean that historical input objects were reused. PNG, PDF and SVG are produced
in `figures/`. Rerender a completed run with:

```bash
python3 workflow/07_networks/read_based/render_figures.py \
  --run /path/to/completed/network/run --output local_runs/network_figures
```

## Fixed design and interpretation

Eligibility requires detection in at least 26 of 51 samples, total fractional
counts >=10 and positive log1p(RPK) variance, retaining 7,802 KOs. Rank variance
after centering over eligible KOs, select 200, then recenter log1p(RPK) over those
200. Selection uses no index, environmental axis or restoration label. Only 27
of these 200 overlap the historical selected set. Deleting one profile before
reselection retains 87–94.5% (median 91%) of the primary set.

Three matched constructions use original CLR, locality/depth residuals, and
residuals additionally adjusted for all five environmental axes, including NDMI.
Design ranks are 1, 5 and 10. Standardization and nuisance fits are repeated
inside each leave-one-out or bootstrap network. Primary settings remain 200
KOs, ridge lambda 0.10, epsilon 1e-8, and 5% density (995/19,900 edges per sample).
LIONESS pseudovalues can exceed [-1, 1]. The mean unweighted degree is fixed at
9.95 and cannot be interpreted as a changing network property.

The protocol was fixed before inspecting the corrected network–index results.
It retains 99,999 primary profile permutations, 4,999 fixed-network bootstrap
draws, 499 complete reconstructions and 200 consensus bootstrap draws. Profiles
are resampled within locality and retain their three depths. Features, HI/MHI
and environmental-axis definitions remain fixed during reconstruction. The
exploratory reconstruction intervals have uncalibrated coverage; duplicated
profiles can shift the LIONESS reference distribution.

Eight one-factor checks use ridge 0.05/0.20, density 2.5%/10%, 100/300 KOs,
the historical four-axis adjustment, and per-million RPK before adding one.
The latter keeps the primary 200 features. Its different result reflects the
interaction of rescaling with a fixed pseudocount; CLR without an additive
pseudocount would be scale invariant. Each family uses synchronized permutations
and jointly corrects all its variants: 60 index, 45 configuration and 30
magnitude/pattern tests. These are sensitivity checks, not parameter selection.
Density variants leave full edge vectors unchanged, so identical configuration
tests are expected and are not independent confirmations.

## Verified outcome, 22 October 2026

All 17 primary acceptance checks and 122 independent numerical checks passed.
The independent QR/Woodbury implementation reconstructed all 3,044,700 LIONESS
weights, all 200 consensus draws per construction and the 99,999 index null
statistics. It also checked selected full reconstructions, metrics, projection
statistics, correction families, resampling restrictions and source hashes.
All 499 reconstruction replicates per construction were valid. Parameter
variants had full-rank designs and no zero-variance KO columns. Execution used
Python with the package versions in `requirements-analysis.txt`; R was not run.

The earlier unconditional headline is not supported by the corrected primary
comparison. Absolute-weight configuration explained 2.52%/2.44% of residual
variation with HI after L+D/L+D+E network construction (maxT P=0.2000/0.3522;
Holm P=0.5112/0.5384). Signed weights passed Holm (0.0030/0.0012), but not maxT
(0.2120/0.1248). No primary HI/MHI summary-metric association passed either
joint eight-test correction; all reconstruction intervals included zero.

Two exploratory signals must remain visible. Stable-consensus mean-weight HI
slopes were -0.3405/-0.3165, with maxT P=0.0292/0.0374 but Holm P=0.0752 for
both. Per-million abundance gave absolute-configuration partial R2=0.03239/0.03153,
maxT P=0.0198/0.0284, but Holm P=0.5642/0.4992. Its L+D log-L2 magnitude test
also passed maxT (0.0424), but not Holm (0.5046). The two corrections target the
same family through different procedures; both are retained. Neither is chosen
after seeing the results. Summary-metric tests across all parameter variants
passed neither joint correction. These outcomes support dependence on
preprocessing, edge definition and multiplicity procedure, not a robust general
relationship between HI and network strength or demonstrated ecological recovery.

## Figures and reporting

`Functional_network_configuration` shows unsupervised edge ordination, full-vector
configuration tests with both corrections, and absolute-configuration parameter
sensitivity. HI colors points but does not fit or orient the PCA. The first two
PCs explain only about 9% of variation; inference uses the complete edge vectors.
`Supplementary_Figure10` separates conditional and reconstructed intervals.
`Functional_network_consensus` uses shared node coordinates and retains all 200
KOs; its geometry is a display aid, not ecological distance. Consensus counts
(388/323/257) differ from the fixed 995 selected LIONESS edges per sample.

Keep associations descriptive: coabundance edges are not confirmed biochemical
interactions, natural connectivity is not measured resilience, and HI is derived
from the same microbiomes. No individual KO driver is confirmed here. Full methods,
seeds, correction families and limits are in `read_based/PROTOCOL.md`.

## Alternative algorithm stability

At the user's subsequent request, the separate `ALGORITHM_BENCHMARK.md` fixed
comparisons before inspecting alternative-method stability. Ledoit–Wolf and
OAS were evaluated on the same 200 KOs, 995 edges, 200 paired profile bootstraps,
17 complete-profile deletions and abundance-scaling check. Adaptive shrinkage
was refitted inside every network and LIONESS deletion. HI/MHI did not select
or evaluate methods. Original CLR controls are retained in the full outputs.

| Method | Bootstrap edge Jaccard L+D / L+D+E | Deletion sample-edge Jaccard L+D / L+D+E | Deletion mean-weight rank correlation L+D / L+D+E |
|---|---|---|---|
| Ridge 0.10 | 0.257 / 0.222 | 0.725 / 0.708 | 0.921 / 0.888 |
| Ledoit–Wolf | 0.271 / 0.234 | 0.746 / 0.723 | 0.958 / 0.960 |
| OAS | 0.271 / 0.234 | 0.745 / 0.721 | 0.956 / 0.962 |

These are medians of paired diagnostics, not independent replicate experiments.
Both adaptive methods modestly improve repeatability, particularly the ranking
of sample mean weights, but neither reaches the prespecified >=0.02 improvement
in bootstrap Jaccard in both adjusted constructions. Normalization sample-edge
Jaccard remains only 0.208–0.222 across these methods. The primary estimator is
therefore retained; adaptive shrinkage is a useful sensitivity comparison, not
a demonstrated solution to preprocessing dependence. No alternative HI test was
run or used to choose a method. All 6,089,400 alternative LIONESS weights were
independently checked through QR residuals and direct sklearn covariance output.

Graphical lasso converged in only 4/24 screened combinations (three constructions,
four alpha values, two convergence settings); no one alpha converged across all
constructions. This is a solver feasibility result, not an empirical stability
ranking of sparse networks. SPIEC-EASI/StARS was not executed, and the screen
does not establish that those methods are intrinsically unstable. Each covariance
has 200 features but residual rank only 50/46/41; there are 17 complete profiles.
Sampling stability cannot establish the biological truth of inferred edges.
