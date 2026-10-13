# Collaborator source requests for manuscript figures

This file tracks panels that should be recovered from the original analyst or experimentalist rather than reconstructed from the manuscript raster image.

## Minimum package to request for each panel

Ask for the smallest complete package that permits verification and later revision:

1. original plotting script or notebook;
2. source table used by the plot, preferably with one row per biological replicate or sample;
3. raw source when practical (`.fcs`, microscopy image, count matrix, etc.);
4. editable vector output or editable panel source;
5. statistical-test output, including test name, sidedness, contrasts, multiple-testing correction, and exact values;
6. sample labels, biological-replicate definitions, and plotted ordering;
7. relevant package/software versions;
8. a short note identifying the source data location and the panel version used in the current manuscript.

Do not request large raw sequencing files when a count matrix and analysis receipt are sufficient.

## Genome-wide and validation screen team

Panels:

- Fig. 1b
- Fig. S1a–d
- genome-wide-screen contribution to Fig. S7a

Request:

- canonical gene-level MAGeCK result table used in the current manuscript;
- guide-level screen counts or normalized guide ratios underlying the gene-level analysis;
- focused-validation-screen replicate tables;
- scripts for the volcano and cross-screen/replicate correlation panels;
- original FCS/FlowJo workspaces or final editable flow panels for Fig. S1a,b;
- exact fields used for directional FDR and gene-effect values.

## Hippo experimental team

Panels:

- Fig. 3e
- Fig. 4a,b
- Fig. 5a–c
- Fig. S10
- Fig. S11
- Fig. S13a
- editable schematics in Fig. 3a and Fig. 5a where applicable

Request:

- replicate-level flow and qPCR tables;
- FCS files and FlowJo workspaces for representative histograms/gating panels;
- exact gates used for CD144 and ETV2-mCherry positivity;
- qPCR normalization tables before and after YWHAZ normalization;
- substrate stiffness, density, treatment-window, and replicate metadata;
- exact Welch-test outputs and the mapping from comparisons to significance bars;
- full-resolution immunostaining images and editable assembly for Fig. S10;
- editable vector outputs for all plots and schematics.

## Multiome / scE2G / ChromBPNet team

Panels:

- motif portions of Fig. 3d and Fig. 4e
- Fig. S5
- motif portions of Fig. S6

Request:

- final gene–motif score and enrichment tables used for plotting;
- motif-name and motif-family mapping table;
- exact program/day selections and row/column ordering;
- scripts for ChromBPNet performance, seqlet/hit QC, and enrichment heatmaps;
- motif-logo assets and the code mapping logos to columns;
- the commit or software versions for scE2G, ChromBPNet, TF-MoDISco, Fi-NeMo, and motif annotation;
- editable output for the analysis schematic in Fig. S5g.

## Zebrafish team

Panels:

- Fig. 5d
- Fig. S13b–d

Request:

- bulk RNA-seq count matrix, sample metadata, DESeq2 result table, and plotting script;
- exact sample names and treatment labels used in the correlation heatmap;
- original ISH images at full resolution;
- blinded embryo-scoring table for posterior and anterior views;
- script or statistical receipt for Fisher exact tests;
- editable plot outputs.

## Tracking fields

Record each response in `PANEL_REPRODUCTION_MATRIX.csv` by updating:

- `candidate_source`
- `status`
- `blocking_item`

A collaborator panel is complete only when the plotted values and statistical annotations can be traced to a source table or analysis output, even when the final panel itself is assembled manually.
