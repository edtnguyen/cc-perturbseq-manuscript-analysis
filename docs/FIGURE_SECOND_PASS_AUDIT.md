# Figure reproduction second-pass audit

## Scope

This audit covers the second implementation pass on branch `curation/figure-reproduction-v1`. It completes the remaining figure panels for which the source analysis, calculation, or required data schema could be identified without choosing among conflicting candidate inputs.

The implementation remains separate from numerical validation. Panels are recorded as `IMPLEMENTED_NOT_YET_NUMERICALLY_VALIDATED` until they execute on Sherlock against the final `/oak` inputs and are compared with the manuscript panels.

## Implemented in this pass

| Figure | Panels | Implementation |
|---|---|---|
| Figure 1 | b | Genome-wide screen volcano using U-test log2 effect and direction-specific MaGECK FDR; pathway highlights; rasterized background points |
| Figure 1 | c | Guide-level Day 3 EC-percentage distributions for the selected regulator groups, including control-reference line and FDR annotations when present |
| Figure 1 | d | Integrated program atlas using the highest activity among the curated D0–D3 cNMF programs |
| Figure 2 | d | Guide-level Day 3 EC fractions for ETV2, KDR, PLCG1 and PRKCD |
| Figure 2 | e | Shared-coordinate UMAPs for P17, P27, SLC8A1 and GATA4 |
| Figure 2 | f | Carrier-seeded same-day neighborhood enrichment for ETV2, KDR, PRKCD and PLCG1 using 300 nearest cells in z-scored 50-program usage space |
| Figure 3 | b | P2/P34 and CDH5/ACTC1 UMAPs plus Hippo-target carrier-neighborhood maps |
| Figure 3 | c | Guide-level Day 3 EC fractions for NF2, LATS2, AMOTL2 and TAOK1 |
| Figure 4 | d | P15/P17 and FLI1/SLC8A1 UMAPs |
| Figure S1 | c | Validation-screen biological-replicate comparison |
| Figure S1 | d | Genome-wide versus validation-screen comparison |
| Figure S2 | b–f | hiPSC circular-capture versus CROP-seq guide recovery, detection, transcriptome recovery, saturation and on-target knockdown |
| Figure S2 | g–h | TeloHAEC circular-capture versus CROP-seq guide and transcriptome recovery |
| Figure S2 | i–n | Parse Evercode v3 primer-condition, guide-pool specificity, transcriptome-recovery, saturation and guide-pool-mixing panels |
| Figure S3 | a–e | Barcode-rank, UMI/gene, QC-distribution, retained-cell and guide-capture panels from standardized preprocessing-QC inputs |
| Figure S3 | f | All-target SCEPTRE knockdown-efficiency distributions by day |
| Figure S5 | a–c | Multiome ATAC-fragment, retained-cell and aggregate TSS-enrichment summaries |
| Figure S5 | d–e | ChromBPNet observed-versus-predicted counts and fold-level count/profile performance |
| Figure S5 | f | TF-MoDISco seqlet and QC-passing motif summaries |
| Figure S5 | h | Ward-ordered motif-enrichment heatmaps for D2 P15/P17 and D3 P2/P24/P27/P34 |
| Figure S7 | a | Genome-wide screen effects for canonical PKC isoforms |
| Figure S7 | b | Day-level pseudobulk TPM from raw counts and GENCODE v43 union-exon lengths |
| Figure S8 | b, e | Guide-level Day 3 EC fractions for BMP-pathway and mesoderm-regulator groups |
| Figure S8 | c, f | Program and marker UMAPs for the two regulator groups |
| Figure S8 | d, g | Carrier-seeded neighborhood maps for the two regulator groups |
| Figure S9 | all | NF2, LATS2, AMOTL2 and TAOK1 relative-expression estimates across D0–D3 with model-based intervals |

Figure 2a and the non-deferred portions of Figure S4 were retained from the first pass.

## Deliberately deferred because inputs remain ambiguous

| Figure | Reason |
|---|---|
| Figure 2b | The canonical Day 3 EC-output field and displayed sign convention have not been uniquely identified. |
| Figure 3d | The final signed motif-score matrix, motif order and motif-logo selection have not been uniquely identified. |
| Figure 4e | Same unresolved motif inputs as Figure 3d. |
| Figure S4a | The final combined K-selection metric table has not been uniquely identified. |
| Figure S4b | The final real-target and frequency-matched NTC p-value exports have not been uniquely identified. |
| Figure S4e | The final K=50 cell-by-program usage input for this panel has not been confirmed. |
| Figure S6 | The final signed motif matrix and frozen regulator/gene/motif/logo selections remain unresolved. |

## Panels retained as source assets

Figure 1a, Figure 2c, Figure S1a–b, Figure S2a, Figure S5g and Figure S8a require editable original artwork, FCS/FlowJo outputs, or deliberate redraws rather than inferred computational reconstruction.

## Shared implementation rules

- Upstream `/oak` inputs are read-only.
- Generated outputs are written under `outputs.root` with a panel-specific directory.
- Implemented panels export PDF, SVG, PNG, source-data TSVs and a JSON run receipt.
- PDF text uses TrueType font embedding and SVG text remains editable.
- Dense background point layers are rasterized where appropriate; labels and foreground layers remain vector/editable.
- Carrier-neighborhood analysis uses same-day neighborhoods, Euclidean distance in z-scored K=50 cNMF usage, `k=300`, a `0.5` pseudocount, two-sided Fisher tests and target-wise BH correction.
- The integrated AnnData, K=50 usage and guide-metadata paths recovered from source notebooks are included in `config/paths.example.yaml`; unresolved inputs remain explicit placeholders.

## Static validation

GitHub Actions workflow `Validate figure notebooks` completed successfully for commit `396ba4afb03491b1cda7cd1ffa38e754fc1094fe` (run `38432876612`). The workflow parsed and compiled every Python module under `src/cc_perturbseq` and every Python code cell in the figure notebooks.

This validation establishes valid Python/notebook syntax. It does not establish that every Sherlock path currently resolves, that every input table matches the inferred adapter schema, or that plotted numerical values match the manuscript. Those checks require fresh-kernel execution on Sherlock followed by source-data and rendered-panel comparison.
