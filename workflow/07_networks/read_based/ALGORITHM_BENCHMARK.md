# Algorithm stability comparison

Added at the user's request after the corrected ridge primary analysis was
completed. The criteria below are fixed before comparing alternative-method
stability or testing their HI associations. This is an exploratory benchmark,
not an independent confirmation of the original scientific hypothesis.

## Comparable estimators

Compare fixed ridge (lambda 0.10), Ledoit–Wolf adaptive covariance shrinkage and
Oracle Approximating Shrinkage (OAS). Each estimates partial correlations of the
same 200 selected KO CLR columns with the same three nuisance designs. Refit
residualization, scaling and any adaptive shrinkage on every subset, including
each LIONESS deletion. Do not alter KO selection or HI/MHI calibrations.
The closed-form shrinkage formulas use sample rows; their iid optimality is not
assumed for these clustered profiles. All empirical stability assessment keeps
profiles intact. Shrinkage strength alone is not evidence of ecological accuracy.

## Paired, profile-preserving diagnostics

Use the same 200 within-locality bootstrap draws as the completed primary
analysis. For each full aggregate graph, retain exactly 995 strongest absolute
partial correlations. Report bootstrap Jaccard similarity to the full graph,
chance-adjusted edge overlap, signed-weight correlation, sign agreement among
shared edges, edge frequencies, and the number passing frequency >=0.50 and
sign consistency >=0.80. Shared density prevents the empty-graph stability trap.
Selection frequencies are computed from resampled graphs, not treated as P values.

Delete each of the 17 complete profiles in turn. Reconstruct aggregate graphs
and all 48 LIONESS sample contributions on the remaining profiles. Compare them
with the corresponding original 51-sample estimates: aggregate edge Jaccard,
median sample-edge Jaccard, weight correlation and rank correlation of mean
absolute edge weight/natural connectivity. No HI labels enter these diagnostics.

Assess normalization sensitivity by comparing the original RPK+1 and per-million
RPK+1 inputs on exactly the same 200 KOs. Report aggregate and sample-edge overlap,
weight correlations and sample-metric rank correlations. Do not mistake repeated
calculations on the same profiles for independent biological replication.

An alternative will be described as consistently more stable only if its median
bootstrap Jaccard improves by >=0.02 in BOTH adjusted constructions, its profile
deletion and normalization median sample-edge Jaccard do not decline by more
than 0.02 in either, and its profile-deletion mean-weight rank correlation does
not decline by more than 0.02. Otherwise report the tradeoff rather than naming
a universal winner. This practical threshold is descriptive, not a formal test.
Even a successful candidate requires separate inferential evaluation before
replacing the primary scientific conclusion. No HI P value selects a method.

## Sparse precision feasibility

Additionally screen graphical lasso on all three full-data covariance matrices
at alpha 0.50, 0.20, 0.10 and 0.05. Assess convergence, dual gap, numerical
conditioning and edge count, separately from sampling stability. A default-solver
pilot showed convergence/positive-definiteness failures; retain these failures
and attempt stricter convergence controls, without using HI. A failed screen
does not prove that the mathematical algorithm or SPIEC-EASI is intrinsically
unstable. No result will be labelled SPIEC-EASI: that R package is not executed.

References: Ledoit & Wolf (2004), doi:10.1016/S0047-259X(03)00096-4;
Chen et al. (2010), doi:10.1109/TSP.2010.2053029;
Kurtz et al. (2015), doi:10.1371/journal.pcbi.1004226;
Liu et al. (2010), https://arxiv.org/abs/1006.3316.
Implementations: scikit-learn 1.8.0 covariance functions.
