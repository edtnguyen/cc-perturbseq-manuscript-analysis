# Portable processed-data release: Hippo and mixed hiPSC–K562

## Run order

1. On Sherlock, run `notebooks/submission/Prepare_GEO_Processed_Data.ipynb` in the environment used for the accepted figure notebooks. The existing Hippo reconstruction must already be complete.
2. Run `notebooks/submission/Test_GEO_Processed_Data.ipynb`. It reads only the two H5MU files, uses the existing panel producers, and writes numerical comparison reports and figure previews.
3. Review those previews against the accepted panels before submission. Deposit the validated H5MU files with the raw FASTQs and the required GEO study/sample metadata.

Default output is `<repository>/geo_release/`, which must remain outside Git. Set `CCPERTURBSEQ_GEO_DIR` to use another durable directory. On another computer, clone the same repository revision, install the analysis dependencies, and put the downloaded objects in `geo_release/`. The test notebook does not load the Sherlock path configuration.

No actual study H5MU is supplied by the code commit: the count matrices are on Oak, and the export notebook creates the study objects there. Synthetic CI fixtures are tests only and must never be submitted as study data.

## Hippo_timecourse_integrated.h5mu

`mod['rna']` retains the original combined broad-QC AnnData, including X, gene annotations, all available layers, embeddings, guide-related arrays and unstructured metadata. X remains the original unnormalized expression matrix; values are not coerced to integers. Metadata are checked against the frozen projection and augmented with `used_in_Figure_3f` and `used_in_Figure_4c`.

`mod['cnmf']` retains the original cell-by-50-program projected AnnData. X is normalized top_specific_200 usage. The original obsm raw usage, normalized usage, and CLR matrices and their explicit P0–P49 order are retained. No cNMF fitting occurs during export. The complete TPM reference is embedded in `varm['reference_gene_spectra_tpm']`, with gene names in `uns['reference_gene_names']`; the selected normalized fixed-H matrix is in `varm['fixed_H_selected_normalized']`, with ordered genes in `uns['fixed_H_gene_names']`.

All broad-QC cells remain available, including cells outside the displayed time points and sublibrary s8. The original condition-specific QC and s8 exclusion are applied only by the original panel producers. Do not treat all deposited cells as panel-included cells.

This file reproduces **Figure 3f and Figure 4c**. It does not claim to reproduce the main perturbation-screen UMAPs, regulator/motif panels, or flow-cytometry panels in Figures 3–4, which depend on other experiments.

The original marker-expression normalization is performed separately for each panel's selected cells. A generic global normalized-expression layer would not be an equivalent substitute. The exporter instead embeds the original marker summaries and the actual Scanpy dotplot color/size matrices as reference tables.

## iPSC_K562_mixed_integrated.h5mu

`mod['rna']` contains the source-QC-passing raw RNA matrix and gene/cell annotations. `obs` includes primer condition, frozen transcriptomic lineage, and every per-cell guide metric used for S2i/j/l/m/n.

`mod['crispr']` contains a **sparse** cell-by-guide matrix of original guide UMIs, summed over duplicate cell-guide rows but not thresholded. Its `var` retains the source guide reference and the guide-pool identity. All QC-passing RNA cells are retained, including zero-guide cells. Cells and order match RNA exactly. No binary guide-assignment rule is invented.

`uns['saturation']` contains all four parsed Split-pipe saturation reports. `uns['raw_saturation_reports_json']` also preserves their original text. These files are required for S2m: sequencing-depth saturation cannot be inferred from only the final expression matrix.

`uns['cluster_annotation']` and `uns['manuscript_median_diagnostics']` preserve the available classification summary and non-blocking historical diagnostics. Supply `S2_FROZEN_METRICS` in the export notebook to use labels from an accepted S2 run; it must include every QC-passing cell, not one filtered panel's table. From that accepted session:

```python
result.cell_metrics.to_csv('Figure_S02_all_cell_metrics.csv.gz')
```

Without this optional table, the existing S2 classifier is run once during export. Public reproduction uses the stored lineage labels; it never reclusters.

The independent public-file check recalculates total and top-guide counts, matched-pool totals, specificity, and the >=3-UMI per-cell-guide threshold from the sparse CRISPR matrix, then calls the existing panel producers. It does not reconstruct an enormous dense guide matrix.

The raw RNA count matrix from this experiment is deposited in full for the QC-passing cells. The source notebook removed features without gene names before calculating QC; this source filtering is retained rather than silently changing the feature universe.

## Shared contents and validation

Both objects contain sample summaries, source-file checksums, software versions, implementation checksums, a schema identifier, full numeric pre-export panel references, and checksums for the stored matrices. Inherited paths are historical provenance, not active input dependencies.

Files are written to a temporary H5MU, reopened, validated, and renamed to the requested final filename only on success. Existing outputs are not overwritten unless explicitly requested. The exporter emits `.validation.json`, `.h5mu.md5`, and `.h5mu.sha256` sidecars. The test notebook adds the original-style PDF/SVG/PNG panels and regenerated CSV source tables.

Tests distinguish numerical reproduction from scientific acceptance. Roundtrip agreement confirms the deposit reconstructs the current accepted producers; it does not independently establish that their results match every manuscript claim. Minor historical cell-count differences do not cause failure. Broken cell correspondence, a mismatched cNMF basis, invalid counts, missing required tables, or changed reconstructed panel values do cause failure.

## GEO submission boundary

The raw FASTQs remain a separate data set for each experiment. Preserve the reads that carry cell barcodes/UMIs and guide identities, not just the transcript read. Submit complete study/sample metadata and explain how pooled sequencing libraries, sublibraries, primer conditions and biological samples map to the processed cell annotations. The software does not invent GSM/SRR accessions, FASTQ filenames, treatment doses, instrument details, or sample-to-FASTQ mappings.

GEO's submission portal describes raw data, processed data and metadata as submission components. H5AD and H5MU are present in public GEO processed-data records, but a portable object alone is not a complete submission.

Reference documentation:
- https://www.ncbi.nlm.nih.gov/geo/info/seq.html
- https://submit.ncbi.nlm.nih.gov/
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305370 (example of H5MU and H5AD processed files, not evidence about this study)

## Dependencies and tests

Use the analysis environment that already runs Figures 3f/4c and S2, plus `mudata`. Synthetic CI additionally tests file roundtrips and public regeneration after the original inputs have been removed.

```bash
PYTHONPATH=src python -m pytest -q tests/test_geo_submission.py
```

The public test needs the repository's existing plotting dependencies (`scanpy`, `plotnine`, `seaborn`, etc.), but no external scientific data files. Hippo regeneration materializes temporary local H5ADs to call the unchanged path-based producers; allow disk space for that extraction. Run the large study export on a compute node.

## Figure-port revision

S2 plotting now uses the re-scripted plotnine paper-panel source rather than the
previous seaborn approximation. Run S2, review its outputs, and pass the saved
`Figure_S02_all_cell_metrics.csv.gz` as `S2_FROZEN_METRICS` when preparing the
release. Existing mixed-cell H5MU exports retain the labels/reference tables of
the producer version that created them; regenerate that release after accepting
the corrected S2 run. Hippo calculations are unchanged.
