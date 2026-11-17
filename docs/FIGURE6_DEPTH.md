# Figure 6 and paired nitrogen trajectories

Run `python3 21_figure6_depth.py --expanded-dir EXPANDED_RUN --depth-dir NITROGEN_DEPTH_RUN --out-root OUTPUT_ROOT` with Python, numpy, pandas and matplotlib. Inputs are external outputs of scripts 15 and 19; no statistics are refitted. The output directory starts with script number 21 and a UTC timestamp.

The 12 displayed KOs were selected after testing because the primary 40/5 cm contrast passed Holm98. This is a selected display, not an independent validation. All 49 eligible KOs remain in the full nitrogen test table. Holm1409 combines 98 depth tests with the 1,311 expanded-panel tests; it is not a correction over every analysis in the manuscript. Pointwise intervals are not simultaneous intervals.

## Main figure caption

Figure 6 | Biogeochemical gene representation across restoration stages and depths. a, Twelve examples selected by process identity. b, All 437 eligible KOs, counted once and ordered by identifier. Colours in a–b show locality-adjusted stage/conserved ratios from mean log₂ normalised abundance across three horizons per profile; conserved equals 1. Both panels share a scale saturated at ¼ and 4. D, E, I and A denote degraded, early, intermediate and advanced stages. Grey cells lack eligible estimates. No primary stage test retained significance after Holm correction across 1,311 tests; displayed pairwise estimates were not separately tested. c, Twelve nitrogen KOs selected after testing because their 40/5 cm contrasts passed Holm correction across 98 depth tests. Points show median paired ratios for 17 profiles; bars show pointwise 95.1% order-statistic intervals under a common-median model. Filled symbols indicate Holm-adjusted P < 0.05 across 98 tests; diamonds also pass the combined 1,409-test correction. Open symbols do not pass the 98-test correction. Bold labels retain combined-family support across all five quantification variants. Daggers mark provisional fractional contributions. Ratios below 1 indicate lower representation at depth. Quantities describe genetic potential, not process rates.

## Supplementary figure caption (number to be assigned in the supplement)

Nitrogen depth trajectories in 17 sediment profiles. Panels show the same 12 KOs selected for Figure 6c. Each line connects the three observed horizons of one profile; colours identify locality. Normalised representation (KO RPK divided by the median RPK of 76 reference genes) is expressed relative to the same profile at 5 cm. All profiles therefore start at 1. Black diamonds show the median ratio across profiles. Axes share the same logarithmic scale, with actual depth spacing. Lines connect observations and are not fitted linear trends. Daggers identify provisional fractional contributions. These plots describe vertical variation; they do not test stage-specific differences or process rates.

## Validation

The plot checks that selection matches all 12 supported primary 40/5 cm contrasts, preserves 17 profiles at each depth and exports the exact 24 contrasts and 612 plotted observations. Figure 6a–b preserve the expanded-panel estimates. Discussion text was not edited. Figures and source data remain external to the GitHub repository.
