# 809 Residualized functional networks

This is an exploratory sensitivity analysis specified before inspecting the new network–index results. The hypotheses and resampling budgets below are fixed for this execution; no scenario will be selected on significance. Networks are not assumed to measure recovery or activity.

## Three matched constructions

Use the exact 200 KO CLR matrix selected by 801. No reselection, second logarithm, replacement of indices or change to sample membership. Use all 51 samples from 17 profiles. The fixed feature universe isolates the influence of upstream covariate adjustment and does not test whether a residual-variance-based feature selection would behave differently.

1. `unadjusted`: historical CLR networks, reproduced as a numerical control.
2. `space_depth`: residuals of each KO CLR regressed on intercept, categorical locality and categorical depth (rank 5).
3. `space_depth_environment`: additionally include the four existing vegetation/landscape, water/inundation, physicochemical and nutrients/redox PC1 axes (rank 9).

Continuous environmental predictors are standardized to aid conditioning; with an intercept this does not change the fitted linear subspace. No HI, MHI or restoration stage enters network residualization. All 200 KO regressions use the same sample-level design. Adjustment removes fitted linear mean effects, not every possible nonlinear, interaction, spatial or taxonomic effect.

## Network inference

For every full network, leave-one-sample-out network and reconstructed bootstrap network, refit nuisance coefficients on that subset before standardizing residual KO columns and estimating ridge partial correlations. Use lambda 0.10, epsilon 1e-8 and the original covariance convention. LIONESS pseudovalue = n * full edge − (n−1) * leave-one-out edge, where both fits include their own residualization. Standardizing full-data residuals once and deleting their rows is not the implemented method.

Retain the 995 strongest absolute pseudovalues in each sample, with deterministic KO-ID tie breaking. Calculate mean absolute retained weight and weighted natural connectivity of the absolute adjacency matrix. LIONESS pseudovalues are not restricted to [−1,1].

## Association comparison

For direct comparison, all three constructions use the same downstream standardized models: metric ~ HI or MHI, separately, + depth + locality. This isolates upstream network construction. The environment-residualized construction is not labelled an estimated causal direct effect.

Use the archived 99,999 within-locality profile permutation plan from 803, preserving depth alignment. Networks and nuisance fits are fixed during these index-label permutations; the network algorithm does not use the permuted index. Conditional null inference retains the original profile-exchangeability assumptions and does not establish independent validation of microbiome-derived HI.

Report four-test Holm/maxT families per construction for historical comparability. The principal new-comparison family jointly covers eight tests (two adjusted constructions × two indices × two metrics), with Holm and maxT. Baseline results are a control, not part of that eight-test family. Do not declare an association robust solely because one adjustment passes a within-scenario threshold.

Use 4,999 profile bootstrap draws within locality for intervals conditional on the sample networks. Additionally perform 499 paired bootstrap reconstructions of the entire residualization–ridge–LIONESS–edge-selection procedure for ALL THREE constructions on the first 499 of those same draws. Derive exploratory reconstruction-aware coefficient intervals and paired changes relative to baseline. Features, indices and the environmental axes themselves remain fixed; the original HI/MHI construction is not retrained. The 499-reconstruction intervals have limited tail precision and are reported alongside, not substituted for, 4,999 conditional intervals. Count rank-deficient or nonfinite replicates explicitly.

## Stable edges and depth

Construct stable consensus masks for each construction using the same 200 profile bootstrap draws stratified within locality, refitting residualization each time. Selection frequency >=0.50 and sign consistency >=0.80. These common draws differ from the historical unstratified 802 draws, so consensus counts must be compared with the new baseline consensus, not attributed solely to residualization relative to historical 356 edges. Retain all nodes. Test HI against both stable-consensus metrics with 4,999 permutations and conditional bootstrap intervals; four new tests across the two adjusted constructions form a separate Holm/maxT family. Empty or numerically degenerate consensus sets produce unavailable results, never an alternate cutoff.

Report depth-specific stable-consensus mean-weight HI slopes and HI×depth interaction as exploratory. Interaction tests use Freedman–Lane profile residual permutations within locality (4,999), retaining main effects. Holm across the two new interaction tests. Within-depth HI slopes adjust for locality; BH spans the six slopes of the two adjusted constructions. Baseline strata are descriptive controls. Compare interactions, not significance differences among strata.

## Global configuration and decomposition

Repeat absolute-weight and signed-weight Euclidean configuration tests and retained-edge square-root-Jaccard tests versus HI, adjusted for depth/locality, with 4,999 synchronized profile permutations. Report per-construction three-test corrections and a six-test family across the two new constructions. Repeat log-L2 magnitude and L2-normalized absolute pattern tests with 99,999 profile permutations, with two-test within-construction and four-test across-new-construction Holm/maxT families. Distinct response R² values are not additive. All use the full 19,900-edge vectors except retained-edge presence.

Environmental correlations are descriptive Spearman diagnostics only in this stage. Since environmental axes directly enter one network construction, do not present their marginal P values as independent validation or reuse an unqualified environmental-permutation test. Individual-KO functional attribution is not added in this sensitivity stage; it would require a separately specified interpretation and multiplicity family.

## Limits and reporting

Keep baseline, both residualizations, all correction families, predictor variance diagnostics, global/leave-one-out/bootstrapped design ranks, zero-variance counts and feature identities. A large reduction in MHI residual variance under environmental adjustment reflects overlap with its environmental construction, not proof that MHI is irrelevant. The two residualizations ask different biological questions; neither is inherently the correct recovery network. Report the complete comparison even if no association remains significant. No change to the 808 matched-reference analysis is made.
