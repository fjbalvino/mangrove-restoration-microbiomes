# Environmental associations: enrichment and component evidence

## What changed

The original gene models and ORA remain distinct from the new descriptive component assessment. The current Figure 5a retains all 15 primary KO contrasts. Figure 5b retains all 11 candidate module contrasts and displays their primary and single-KO enrichment q values beside block coverage and per-sample representation. The original module bubble plot remains in the external archive.

The primary ORA was rerun from the pinned inputs. All six complete statistical tables matched the archived reference at the declared numerical tolerances. The result remains 102 KO/PFAM contrasts (96 terms) and 11 module contrasts (10 modules). Ten of those module contrasts remain eligible but nonsignificant in the single-KO sensitivity; M00417 becomes ineligible. The original thresholds and multiplicity families were not relaxed.

The new analysis evaluates the full required/alternative components of the 10 modules, selected after the original enrichment. It uses the complete CF KO matrix for representation and the tested 100,000 genes for associations. Its definitions were retrieved on 16 November 2026 (17 November UTC), separately from the October enrichment reference. All ten definition KO sets match the frozen October KO–module links. Existing annotations were not reassigned.

## Evidence levels

- Gene association: standardized common-slope estimate, confidence interval and original BH q; interaction and depth-slope families remain separate.
- ORA: overrepresentation relative to tested, usable, ontology-annotated genes; enrichment is not abundance change.
- Representation: three-state required/alternative rule in one community sample, using fractional or conservative single-KO bounds. Co-detection does not establish a pathway in one organism.
- Associated block coverage: requirements covered by same-sign qualifying gene assignments. This is descriptive, with no new block-level P value, and is not a test of coordinated activity.

Six module rules are met in all 51 samples under primary single-KO evidence and retain that representation across both-mate and supported-partial-consensus variants. Glycolysis and ectoine have negative single-KO gene associations covering 9/9 and 5/5 blocks along NDMI and water-signal axes. This does not restore their nonsignificant single-KO ORA. The complete inventory retains incomplete and unresolved blocks, opposite signs, nonsignificant gene estimates and untested components. KO abundance summaries preserve restoration order and depth without inventing a pathway abundance or new stage test. No inference of equivalence or measured activity follows.

## Reproduce

Use the original external analysis archive and the updated Supplementary_Data1.zip, plus the original CF input export already required by the expanded panel. The CF data root contains metadata_003b.csv, matrices/KO_*.tsv and 03_unit_annotations/annotation_patterns.json. Definitions and reference genes are external files; no database or figure is added to GitHub.

```bash
python3 workflow/09_environment_enrichment/09_01_run_component_evidence.py \
  --analysis /path/to/Figura4 \
  --reference /path/to/Supplementary_Data1/reference \
  --reference-results /path/to/Supplementary_Data1/results \
  --cf-data /path/to/CF_export \
  --definitions /path/to/Supplementary_Data1/component_evidence/reference/definitions.json \
  --reference-genes /path/to/uscg_76_kegg_min_2313.lst \
  --out-root local_runs
```

09_01 orchestrates the independently executed producer steps. 09_00 alone reproduces original ORA and its historical module bubble plot; it does not create the revised Figure 5. 22 writes a timestamped evidence directory and LATEST.txt; 24 writes a separate timestamped figure/HTML directory and LATEST_24.txt. The new launcher was syntax-checked; its component commands were executed directly in this update.

## Verification and limits

Eleven original input hashes passed. All archived BH families were recalculated from original gene-model P values: 42,222 common slopes, 63,942 depth slopes, 1,234 interactions, 48,580 distinct environmental genes and no MHI discoveries. A fresh 714-fit comparison (1,785 contrasts) with the independent R benchmark passed; the full catalogue and all 1.2 million models were not refitted. Forty-four logic checks cover AND/OR, optional terms and unknown values. The component output contains 3,060 sample/rule/variant/evidence rows.

The ordered-trend producer underlying Supplementary Table 11 was recovered from Analisis_Figura4.zip and is preserved in workflow/06_gene_associations/ordered_trends/. Its KO input equals the CF complete-consensus RPK matrix after transposition. Existing results have minimum q=0.408183 for KO and 0.219292 for PFAM. Independent exact-test unit checks were rerun; the full ordered-trend analysis was not rerun in this update.

The manuscript, supplementary methods, Table 15 explanation and Figure 5 were updated together. Discussion was preserved byte-for-byte in the paragraph text. Sparse replication motivates explicit uncertainty, not weaker thresholds or an assumption that every unsupported pattern is caused by sample size.

## Method reference

KEGG module syntax: https://www.kegg.jp/kegg/module.html . Required pathway blocks and subunits use AND, alternatives use OR, and optional components do not determine rule satisfaction. Definitions, raw snapshots, hashes and the crosswalk against frozen links are in the external component-evidence archive.
