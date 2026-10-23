# programDE source audit for Figures 2–4, S4, and S6

## Scope

This audit records the source notebooks inspected in `edtnguyen/programDE` on branch `sherlock-sync-20261013-070509` and the implementation state in `edtnguyen/cc-perturbseq-manuscript-analysis` on branch `curation/figure-reproduction-v1`.

The port is limited to analyses supported by identifiable source code. It does not infer missing source tables, reverse-engineer values from rendered figures, or mark numerical reproduction as complete before Sherlock inputs are resolved and outputs are compared with the manuscript panels.

## Source mapping and implementation state

| Requested panel | programDE source located | Implemented in curation branch | Remaining dependency |
|---|---|---|---|
| Figure 2a | `notebooks/program_level_temporal_map_notebooks/02_gene_program_effect_temporal_heatmap.ipynb` | Pairwise target correlation from gene-by-program/time effect vectors; correlation-distance/average-linkage ordering; pathway zooms; source-data export | Run against the final all-day programDE export and compare target ordering and group membership with the manuscript |
| Figure 2b | Same temporal-map notebook and curated program sets | Regulator-by-program heatmap, q-value stars, pathway blocks, and Day 3 EC-output column | Confirm the canonical EC-output table, displayed effect column, and sign convention |
| Figure 3d | `notebooks/program_annotation_D2_P15_P17_D3_P2_P34.ipynb` and `notebooks/program_annotation_panels_D0_D1_D2_D3_rotated_logos_D2_extra_D3_trimmed_markers.ipynb` | Aligned P2/P34 programDE effects, normalized cNMF loadings, motif enrichment, significance stars, and source-data export | Confirm exact motif columns/order and signed motif-score input; match motif-logo placement |
| Figure 4e | Same program-annotation notebooks | Aligned P15/P17 programDE effects, normalized cNMF loadings, motif enrichment, significance stars, and source-data export | Confirm exact motif columns/order and signed motif-score input; match motif-logo placement |
| Figure S4a | Evaluation logic was documented, but the exact final K-selection table was not located | Generic four-metric K-selection plot | Supply the final K, stability, reconstruction-error, explained-variance, and GO/gene-set-enrichment table |
| Figure S4b | `notebooks/iPSC_EC_QQ_full_covariates_with_replicate_correlation_bundle/` | Day-specific real-target and matched-NTC Q–Q plots | Resolve final p-value files/directories in `paths.local.yaml` |
| Figure S4c | Same full-covariate bundle | D0–D3 replicate scatterplots, Pearson correlations, counts, and source-data export | Run against final pair table and compare inclusion filters |
| Figure S4d | K=50 cNMF spectra | Pairwise top-300-gene overlap | Run against final spectra and verify program ordering |
| Figure S4e | K=50 cNMF usage matrix | Pairwise Pearson correlation of program usage | Resolve final usage table or AnnData and verify all retained cells |
| Figure S6 | Four-day integrated annotation notebook | Day-grouped programDE effects, cNMF loadings, motif enrichment, significance stars, and per-layer source data | Freeze regulator, marker-gene, motif, and motif-logo selections/order per day |

## Preserved source conventions

- Figure 2a uses correlation distance and average linkage on target-level perturbation-effect vectors.
- ProgramDE effects and q-values are consumed as outputs; figure notebooks do not rerun the CRT.
- cNMF loadings are normalized within each displayed program.
- Significance thresholds are `* q < 0.05`, `** q < 0.01`, and `*** q < 0.001`.
- Integrated annotation aligns perturbation effects, cNMF gene loadings, and TF-motif enrichment.
- Implemented panels write PDF, SVG, PNG, source-data TSVs, and a run receipt.
- Receipts remain `IMPLEMENTED_NOT_YET_NUMERICALLY_VALIDATED` until final Sherlock execution and panel comparison.

## Inputs requiring confirmation

1. Figure S4a final cNMF K-selection metrics table.
2. Figures 3d, 4e, and S6 signed motif-score matrix, final motif identities/order, and logo selection.
3. Figure 2b canonical Day 3 EC-output table, effect column, and sign convention.
4. Figure S4b/e final real-target/NTC p-value exports and final K=50 usage matrix.

## Validation gate

Promote a panel from `IMPLEMENTED_NOT_YET_NUMERICALLY_VALIDATED` only after all configured inputs exist, the notebook executes without hidden state, source-data dimensions and identifiers are checked, statistics and ordering are compared with the manuscript, and vector output is inspected for editable text and intended rasterization.
