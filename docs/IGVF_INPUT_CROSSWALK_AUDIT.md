# IGVF manuscript-input audit

Public API checked at 2026-11-01 04:49 UTC (October 31 evening in California). Figure code audited at main commit `ebdc6a74fc7ed5296c671e249aa78eda5772f4d4`. Main and its notebooks were not modified.

Scope excludes the Hippo time-course and mixed hiPSC-K562 datasets on GEO, and the motif-enrichment/map/logo files that will be bundled separately. Non-motif dependencies of Figures 3d, 4e and S6 remain in scope.

## Verification

Executed preparation/upload notebooks were compared with current public IGVF metadata. All 16 original CC-Perturb-seq staged-file size/MD5 receipts agree with current portal metadata. All 35 mapped core files below are currently released with upload_status validated. Separately, 29 public tabular payloads were downloaded and their MD5s checked. D0 H5AD and H5MU annotations and HDF5 schemas were inspected using HTTP Range requests, not full matrix downloads. Exact equivalence of the independent usage/consensus TSV representations to H5MU matrices is not established.

Primary API records are `https://api.data.igvf.org/<accession>/?format=json`. Read-only workflow runs 41274245832 and 41274454185 retain the metadata, table headers, checksum checks and HDF5 schema evidence as artifacts. No portal credentials were used.

## Day-specific analysis and preprocessing sets

|Day|Supplied analysis set|Linked preprocessing analysis set|Cells|
|---|---|---|---:|
|D0|IGVFDS0463XTUH|IGVFDS5257JGQO|38609|
|D1|IGVFDS0991DBQZ|IGVFDS3530HHEU|16925|
|D2|IGVFDS4402ZGWJ|IGVFDS7380DUFB|21980|
|D3|IGVFDS2307REVJ|IGVFDS6826YKCZ|44102|

Current portal ownership places the H5MU in the supplied analysis set; annotated H5AD, guide quantifications and per-cell QC are in its linked preprocessing set.

|Day|Annotated H5AD|cNMF H5MU|Guide quantifications|Per-cell QC|
|---|---|---|---|---|
|D0|IGVFFI7756RKHT|IGVFFI8722NAWN|IGVFFI1608GIRE|IGVFFI3905TBLN|
|D1|IGVFFI4783YGLZ|IGVFFI6740UUOG|IGVFFI6431OXJN|IGVFFI0627NEZT|
|D2|IGVFFI5165BSGC|IGVFFI7155BENM|IGVFFI1924OSEM|IGVFFI1451CYGW|
|D3|IGVFFI9086SWZC|IGVFFI2345GXDF|IGVFFI7077LERW|IGVFFI8282ATZQ|

The executed `cc_perturbseq_upload_ready_builder/notebooks/01_make_igvf_upload_ready_files.ipynb` reads the integrated `3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad` and March `030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/adata/cNMF_50_0_2.h5mu`. It partitions by day while preserving cell IDs and requiring exact H5AD/H5MU cell-set agreement. It does not independently recompute UMAP or refit cNMF. The subsequent `02_upload_igvf_analysis_files_reuse_existing_software.ipynb` records the submitted-file accessions and receipts.

Staging root is `/oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/<DAY>/`. File suffixes are `cc_perturbseq/CCPerturbSeq_<DAY>_annotated_cell_gene_guide.h5ad`, `perturbnmf/PerturbNMF_<DAY>_cell_gene_program.h5mu`, `cc_perturbseq/guide_quantifications.tsv.gz`, and `cc_perturbseq/per_cell_quality_report.tsv.gz`.

The public D0 H5AD has `X`, `layers/raw_counts`, `obsm/X_umap`, `obsm/X_pca`, `obsm/guide_assignment` (38609 x 1922), `uns/guide_names`, biological-sample/day/QC fields and Leiden labels. The public D0 H5MU has `rna` and `cNMF` modalities: cNMF.X is 38609 x 50, var names are 0..49, and cNMF.varm/loadings is 50 x 17538 with genes in cNMF.uns/var_names. Usage row sums are not uniformly one: do not impose a new normalization. Concatenation must preserve IDs, columns and embeddings.

## Program-analysis tables

|Day|Program-regulator associations|Regulator-regulator correlations|Program-program correlations|
|---|---|---|---|
|D0|IGVFFI8214DMKC|IGVFFI1127PXEP|IGVFFI5671ZTPF|
|D1|IGVFFI4053SPMD|IGVFFI3541QUUN|IGVFFI0574VVJA|
|D2|IGVFFI9940ZZPJ|IGVFFI6275JIJY|IGVFFI3823LYPH|
|D3|IGVFFI5389NMWB|IGVFFI5745PRRI|IGVFFI3188WQFI|

These files are in the supplied day sets. Program-regulator tables have 14900 rows per day, with log2FC, p_nominal_nlog10, fdr_nlog10 and fdr_method='Storey Q-value'. They are pooled results, not biological-replicate beta matrices or NTC pseudo-target p-values. Within-day regulator correlations are not the Figure 2a four-day-profile correlation or S4 replicate concordance.

Shared set IGVFDS9440URJB contains IGVFFI3590IELV (17546-gene universe: 17538 fitted genes plus eight additional perturbed regulators) and IGVFFI6391GCJR (15000 rows, top300 genes per program). In the downloaded gene-program table all program_annotation, program_annotation_term and program_annotation_id values are missing. It does not replace curated Programs_annotations.xlsx or complete consensus spectra. Table labels GP_1..GP_50 require explicit correspondence to H5MU labels 0..49; no silent offset is assumed.

## Genome-wide and validation screens

|Set|File role|Accession|
|---|---|---|
|IGVFDS9629VXQK|Genome-wide guide counts|IGVFFI8424MVZE|
|IGVFDS9629VXQK|Genome-wide MAGeCK differential statistics|IGVFFI9001HIYI|
|IGVFDS4708PCVX|Validation guide counts|IGVFFI5157STFU|
|IGVFDS4708PCVX|Validation MAGeCK differential statistics|IGVFFI3975KRIX|

The executed `prepare_and_submit_both_crispri_screens_igvf_analysis_files.ipynb` establishes their preparation lineage. Genome-wide statistics come from `Pos_vs_Neg.gene_summary.txt`; they retain native directional MAGeCK fields but do not include the separate Figure 1b U-test result table `combined_utest_P1P2_adj_test.txt`.

Validation IGVFFI5157STFU contains 1926 guide rows. The builder reconciled UTestAnalysis/guideCounts.csv against MageckAnalysis/mageck_sgRNA_readcount.txt: all nonmissing values agreed, with two missing positive-sort counts filled from the corresponding MAGeCK columns. It recorded 91 guide-ID rewrites and one NEDD8/NEDD8-MDP1 alias correction. Preserve the gene analysis grouping. Deposited-to-local column mapping is D01->D0R1, D02->D0R2, D3neg1->D3NR1, D3neg2->D3NR2, D3pos1->D3PR1, D3pos2->D3PR2. Validation differential statistics are not the guide-count input used by Figure 1c's Welch calculations.

## Guide reference

IGVFFI0180ZBTX in CuratedSet IGVFDS1427MQWJ is the CC-Perturb-seq 1926-guide reference: 1806 targeting, 60 non-targeting and 60 safe-targeting. The deposited intended_target_name uses Ensembl identifiers. It is a counterpart to 300genes_guide_metadata_v43.tsv, not an asserted byte-identical replacement. Preserve exact guide identities, gene-symbol analysis groups and pooled control definitions when adapting it.

## Still unresolved

1. Figure 1b: no corresponding accession identified for `combined_utest_P1P2_adj_test.txt`. Need its accession if elsewhere, or publish the exact small table. Do not replace U-test log2FC with MAGeCK LFC.
2. S4: no corresponding deposit identified for the complete numerical-result bundle containing real/NTC empirical Q-Q points, replicate-beta matrices and manifest, pooled full-day p-values, and target-cell-with-cNMF-usage counts. The public pooled program/correlation tables cannot replace it. Existing numerical results suffice for plotting; another CRT run is not needed simply to publish them.
3. S6: curated `Programs_annotations.xlsx` remains required; the public program table's annotation columns are empty. This workbook or equivalent TSV can be bundled in repository data.
4. March separate usage/consensus TSVs: deposited H5MU contains candidate full matrices, but labels/scaling/value equivalence to the exact TSV inputs still needs checking. The top300 table is not a full-spectrum replacement.
5. Optional S4 raw recomputation: current backend uses January `012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all`, whereas deposited H5MU provenance is March 030526. Original January-model and per-day raw inputs need separate resolution for exact raw recomputation. Do not silently swap model versions.

## No new accession inherently required

`data/figure_02/cNMF_50_0_2_032526_6co.xlsx` is already bundled on main for Figures 2a/b, 3d, 4e and S6. Its SHA256 is d2fe4834b1efd3c00e01062b56b907a6178bfbc48a40beb5ea7b347538bfc824. Exact equivalence to portal program tables was not established; keep the bundled workbook for now.

Figure 1d's current renderer uses repository-fixed manuscript labels/order. The external program_activity_long_v4.csv and program_order.txt requirements are redundant for that rendering and can be removed or bundled; they do not inherently require new IGVF accessions.
