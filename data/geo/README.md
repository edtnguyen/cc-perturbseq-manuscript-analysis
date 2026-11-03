# Optional processed GEO objects

Default figure execution uses the compact tables in
`data/reproduction/geo_plots/`; no GEO download is required.

| Dataset | GEO accession | Notebook / panels |
|---|---|---|
| Parse v3 CC-Perturb-seq, iPSC-K562 mix | [GSE352893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352893) | Figure S2i, j, l, m, n |
| Small molecule Hippo scRNA-Seq | [GSE352939](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352939) | Figures 3f and 4c |

These accessions are prefilled in the corresponding notebooks. Set
`DATA_SOURCE="geo"` to resolve and download the processed supplementary file and
regenerate plotting tables from the submitted integrated H5MU. `GEO_FILENAME`
may be set to an exact filename if discovery reports ambiguous candidates.

Downloads are cached in `data/geo/downloads/` and remain ignored by Git.
Regenerated tables go to `data/derived/geo_plots/`, never over published tables.
Released GEO files need no login. Before release, select `local_file` and provide
`Hippo_timecourse_integrated.h5mu` or `iPSC_K562_mixed_integrated.h5mu` from the
submission workflow. Assigned accessions do not by themselves guarantee public
file availability. Regeneration does not download FASTQs, refit cNMF, or
recluster cells.
