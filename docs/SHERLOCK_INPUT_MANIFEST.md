# Sherlock input manifest for figure reproduction

## Scope

This manifest records the current Sherlock inputs for remaking:

- Main Figures 1–2, Figure 3b–d and Figure 4d–e.
- Supplementary Figures S1–S9.

The genome-wide and focused validation screens are treated as fixed upstream analyses. Their processed result tables will be used directly unless an audit reveals a discrepancy that cannot be resolved from the original analysis documentation.

## A. Genome-wide CRISPRi screen

Base directory:

```text
/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis
```

| Logical input | Sherlock path | Panels | Status / use |
|---|---|---|---|
| MAGeCK, CDH5+ versus CDH5− | `/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Pos_vs_Neg/Pos_vs_Neg.gene_summary.txt` | Fig. 1b; Fig. S7a; Fig. S1d | Primary processed genome-wide result table. Audit promoter identifiers and directional positive/negative fields before plotting. |
| U-test, CDH5+ versus CDH5− | `/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Pos_vs_Neg/combined_utest_P1P2_adj_test.txt` | Fig. 1b candidate x-axis source; Fig. S1d | Upstream U-test result table. Do not substitute for MAGeCK significance. Determine whether the final Fig. 1b x-axis uses this `log2FC` field or the MAGeCK directional LFC. |
| MAGeCK, D3 total versus D0 | `/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/D0_vs_Total/D0_vs_EC_total.gene_summary.txt` | Supporting screen QC / survival audit | Not currently a plotted panel in the requested scope unless needed to reproduce a supplemental QC annotation. |
| U-test, D3 total versus D0 | `/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/D0_vs_Total/combined_utest_P1P2_adj_test.txt` | Supporting screen QC / survival audit | Same restriction as above. |
| MAGeCK output-column README | `/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Mageck_output_README.xlsx` | Audit | Use to verify directional field definitions. |
| Combined guide-count input | `/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Files_used_for_mageck/combined_mageck_input.txt` | Provenance only | Retain for traceability; no need to rerun MAGeCK for initial figure reproduction. |
| Original plotting script | `/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/generatePlots.R` | Fig. 1b provenance | Historical plotting code. It is not automatically canonical because the final manuscript styling and axis definitions changed. |
| Historical volcano outputs | `/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/figures_250115` | Visual comparison | Use only to compare data and labels against the final manuscript. |
| Additional screen plots | `/oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/mageckflute/230920_MAGECK_Flute.R` | Provenance | Rank and related legacy plots. |

### Genome-wide screen caveats

1. The offboarding document reports a duplicated `NEDD8` entry caused by inconsistent underscore usage. The audit must identify and reconcile this before counting genes or joining tables.
2. Some identifiers include promoter suffixes such as `A1BG_-P1`. They must be checked against the remapped hg38 guide metadata before collapsing to gene symbols.
3. Fig. 1b uses MAGeCK FDR on the y-axis. The precise x-axis field must be confirmed by matching the final panel and legend; recovered candidate scripts disagree between U-test `log2FC` and MAGeCK directional LFC.
4. Significance must use the direction-appropriate MAGeCK field (`neg|fdr` for negative effects, `pos|fdr` for positive effects) unless the original result documentation establishes a different canonical summary.

## B. Focused validation screen

Base directory:

```text
/oak/stanford/groups/engreitz/Users/dulguun/231023_Validation_Screen_Files
```

| Logical input | Sherlock path | Panels | Status / use |
|---|---|---|---|
| Main validation result table | `/oak/stanford/groups/engreitz/Users/dulguun/231023_Validation_Screen_Files/d3pos_d3neg_mageckLFC.csv` | Fig. S1c–d | Primary table. It includes validation-screen effects and original genome-wide effect columns from MAGeCK and U-test analyses. |
| Gene-type annotation | `/oak/stanford/groups/engreitz/Users/dulguun/231023_Validation_Screen_Files/validation_genes_type.xlsx` | Fig. S1d | Defines `genome_wide_hit`, `genome_wide_non_hit` and `tri_perturb_seq` categories. |
| Screen-hit summary | `/oak/stanford/groups/engreitz/Users/dulguun/231023_Validation_Screen_Files/ec_screen_hits.csv` | Fig. S1d; regulator-set audit | Lists all 297 targets, original and validation effects, and screen-hit indicators. |
| Original plotting script | `/oak/stanford/groups/engreitz/Users/dulguun/231023_Validation_Screen_Files/generatePlots.R` | Fig. S1c–d provenance | Historical script. It computes linear-model R-squared, whereas the current results text specifies Pearson correlation for Fig. S1d. |
| MAGeCK analysis directory | `/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/MageckAnalysis/` | Provenance only | No rerun needed initially. |
| MAGeCK guide-count table | `/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/MageckAnalysis/mageck_sgRNA_readcount.txt` | Provenance only | Retain for traceability. |
| U-test analysis directory | `/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/UTestAnalysis/` | Provenance only | No rerun needed initially. |
| U-test guide-count table | `/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/UTestAnalysis/guideCounts.csv` | Provenance only | Retain for traceability. |

### Validation-screen caveats

1. The current manuscript states R² = 0.85 for the two validation replicates in Fig. S1c, but describes Fig. S1d as Pearson's r = 0.72 with P = 1.85 × 10⁻40. The historical R script calculates R² for both comparisons. These statistics must not be conflated.
2. Confirm which original-screen effect column is used in Fig. S1d: `Original.FC.Mageck` or `Original.FC.U.Test`.
3. Confirm which validation effect column is used: the historical script uses `Average.U.Test.FC` even in a plot labeled as an original-MAGeCK comparison.
4. Colors in Fig. S1d must match the current legend: significant in the genome-wide screen, validation screen, or both.

## C. Flow-cytometry assets for Fig. S1a–b

The offboarding documentation points to the genome-wide-screen FACS record `230410_EC_Screen_FACS` and to validation-screen project documentation, but it does not provide an exact Sherlock path to the canonical FCS files, FlowJo workspace or final vector export.

Required assets:

- Fig. S1a: unstained, isotype-control and antibody-stained genome-wide-screen FCS files plus final gates.
- Fig. S1b: focused validation-screen biological replicate 1 and 2 FCS files plus final gates.
- Prefer the original FlowJo workspace and editable PDF/SVG exports.

Status: `BLOCKED_INPUT_PATH` until an exact source location is supplied or found.

## D. CC-Perturb-seq core inputs

These paths remain to be consolidated from the recovered code and `programDE` branch:

- final day-specific AnnData/H5MU objects for D0–D3;
- final SCEPTRE guide-assignment matrices;
- final integrated D0–D3 UMAP coordinates;
- K = 50, density-threshold 2.0 cNMF usages and spectra;
- six-covariate programDE CRT result tables;
- guide-level D3 EC-percentage table and 100,000-permutation results;
- carrier-seeded neighborhood tables or inputs required to regenerate them.

These support Fig. 1c–d, Fig. 2a–b and 2d–f, Fig. 3b–c, Fig. 4d, Fig. S3–S4 and Fig. S7b–S9.

## E. Multiome and motif inputs

Recovered code points to the following candidate analysis root:

```text
/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec
```

Candidate final products include:

```text
analysis/2.tf_enrichment_in_programs/enhancer_enrichment_results.tsv
results/compendium/modisco_compiled/modisco_compiled.h5
analysis/0.0.motif_id_annotation/tables/cluster_to_tf_motif_id_mapping.tsv
```

These support Fig. 3d, Fig. 4e, Fig. S5 and Fig. S6. They must be checked against the exact multiome/scE2G/ChromBPNet analysis version described in the current manuscript before being treated as canonical.

The older offboarding paths for Seurat and ArchR objects document earlier 5-time-point multiome processing but are not assumed to be the final motif-analysis inputs.

## F. Immediate audit command

Run the repository script on Sherlock:

```bash
python3 scripts/audit_screen_inputs.py \
  --output reproduction_audit/screen_inputs
```

It performs read-only existence, size, checksum, delimiter/header and identifier checks and writes a machine-readable receipt. It does not modify any upstream result file.
