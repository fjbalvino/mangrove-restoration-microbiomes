# Carbon, nitrogen and sulfur markers

This is the current Figure 6 / Supplementary Figure 10 analysis, following
environmental enrichment. It replaces networks in the manuscript and retains
all 11 markers, provisional candidates, null findings and sensitivity analyses.

See [methods, external inputs and commands](../../docs/FUNCTIONAL_CYCLES.md),
the [analysis specification](analysis_specification.json), and the
[validation record](../../docs/FUNCTIONAL_CYCLES_VALIDATION.json).

Run `10_00_run_cycles.py --help` for the launcher. `10_01_prepare_inputs.py`
can reconstruct the compact input package from the corrected external KO export.
All data and generated figures stay local. The `code/` directory preserves the
reviewed calculations and figure layout, with verified binary-output writing;
`report_tables.py` adds the two essential supplementary source-table assemblies.
