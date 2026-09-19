# Code curation provenance

Repository preparation and input/dependency refactoring were assisted by OpenAI
Codex (OpenAI, 2026). Existing scientific source authorship is retained. No
unverified backend model identifier is claimed.

This minimal publication preserves the numerical algorithms in the canonical
analysis scripts. The documented 809 input refactor recreates edge ordering and
permutations rather than reading fitted historical objects. Checks against the
recovered reference run found identical weights and per-sample metrics in all
three constructions and exactly identical 99,999 permutations. Full reconstruction
bootstraps and R analyses were not rerun for this publication.

The input staging utility copies explicit files and rejects conflicting existing
copies. The collector reads only listed server paths, without analyses or source
modification. File hashes and source correspondence are in CODE_PROVENANCE.tsv
and INPUTS_MINIMAL.tsv. Further changes are documented in MINIMAL_CHANGES.md.

The independent HI/MHI checks reconstruct their fixed calibrations. Historical
calibrations and currently selected taxonomic subcompositions are not interchangeable.
