# Taxonomic Figure 1

`workflow/03_taxonomic_structure/03_05_figure1.py` regenerates Figure 1 B-D
from the unchanged taxonomic object and verified archived numerical tables.
It checks current metadata against the original 301 sample labels and uses
the current `metadata_003b.csv` environmental predictor for panel D.
This section precedes taxonomic change points and the functional analyses.

## Inputs and execution

All analysis inputs were recovered and verified on 2026-10-23. They are
external files: no abundance matrices, R objects, databases, plots or result
tables are added to GitHub by this update. Exact sizes and SHA-256 values are
in `workflow/03_taxonomic_structure/figure1/inputs_required.tsv`. Keep that
manifest's relative layout in an external input directory.

The canonical phyloseq object matches the original 301 MD5
`4480069cd8bb7914b2fd2e81bcb70b89`. The original 201 candidate-selection table
and 202 membership table also match the archived 301 input hashes. Candidate
membership is checked, not reselected. No taxonomic labels are changed.

```bash
python3 -m pip install -r workflow/03_taxonomic_structure/figure1/requirements.txt
python3 workflow/03_taxonomic_structure/03_05_figure1.py \
  --input-dir /path/to/Figure1_inputs \
  --output-dir local_runs/taxonomic_figure1
```

Python 3.12 and the listed package versions were used for verification.
CairoSVG also needs the Cairo shared library. The output directory must be new.
For the original supplied 1605 x 2048 PNG, optionally append
`--source-figure /path/to/Figure_1.png --panel-a-bottom 460` to preserve panel A
and create the complete assembly. The crop retains all non-white pixels of A
and stops before the old B label. B-D are vector artwork in SVG and PDF; A
remains an embedded raster from the supplied figure, not an editable reconstruction
of the Inkscape source. PNG outputs are 4380 pixels wide at 600 dpi.

## What is recomputed and what is retained

The audit checks 51 unique samples, 17 complete profiles, three depths per
profile, and agreement of locality, stage, depth and profile identifiers with
the archived 301 ordinations. Stage profile counts are 3, 3, 5, 3 and 3.
It independently recomputes every profile-mean Hill value, locality-adjusted
mean and confidence interval, stage F statistic, marginal compositional R2/F,
two-dimensional PCoA coordinates and axis fractions, and historical depth-wise
unique adjusted fractions. It verifies the spatial-median dispersion F to
1e-4 absolute tolerance because the independent median solver and vegan use
different convergence criteria. The 39 numerical comparisons passed.

The original profile-level permutation P values and Hill BH q values are
retained from the archived R 301 run. This plotting update does not claim a
new execution of those permutations, bootstrap analyses, or the screening step.
The archived R code remains available in 03_01 through 03_04. Its paths should
be supplied explicitly if running the entire historical workflow on Landa Lab.

For D, the calculation is multivariate least squares on CLR counts with
pseudocount 1. A predictor's unique adjusted fraction is the full model's
adjusted R2 minus that of the model containing the other two predictors.
The implementation first reproduces the archived `vegan::varpart` fractions
using the historical ECI, then uses the current ECI. No new permutation P
values are attached to these descriptive fractions.

## Corrections that matter for interpretation

- The old 302 ECI differs from the current integrated predictor (correlation
  -0.937708, not a simple sign reversal). Using the current predictor changes
  D: locality's unique fractions are 6.4%, 11.9% and 27.4% at 5, 20 and 40 cm.
  The corresponding stage fractions are 8.1%, 3.3% and 3.8%; environment is
  8.5%, -1.1% and 3.4%. These percentages refer to total variation.
- D no longer rescales the three fractions to 100%. The negative estimate is
  retained as an adjusted estimate, not interpreted as a negative ecological
  effect. Shared and residual fractions are not shown.
- The former 40-cm asterisk came from an unrestricted, sequential stage test
  with stage entered before ECI. It did not test the unique adjusted fraction
  displayed in D. Its archived across-depth BH q was 0.060. The new figure
  does not attach this test to the unique fraction.
- 301 uses spatial medians with bias adjustment for dispersion, not centroids.
  Correct the figure legend accordingly; the archived P values remain 0.558
  and 0.244.
- The 12,761-taxon object contains 421 entries assigned to viral kingdoms and
  43 with no kingdom assignment; the candidate set contains 25 viral entries
  and four unassigned. The input therefore must not be described as exclusively
  bacterial and archaeal. This audit inventories existing labels only; it does
  not reclassify or remove taxa.
- The candidate P = 0.030 remains nominal and exploratory. A supported result
  in one matrix and an unsupported result in the other is not a formal test
  of a difference between their stage effects.

Local outputs include the version comparison, numerical validation, input hashes,
existing taxonomy inventory, source tables and PNG/PDF/SVG figures. Online
Methods and Discussion were not modified as part of this figure revision.

The figure caption, Results and relevant supplementary labels should be
reconciled together when integrating the revision into the manuscript. The
old 03_04 script and its contract are retained for historical reproduction;
use 03_05 for the current figure.
