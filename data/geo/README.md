# Optional processed GEO objects

The default figure notebooks use the precomputed plotting tables in
`data/reproduction/geo_plots/`, not files in this directory.

To regenerate those tables, select `DATA_SOURCE="geo"` and enter the experiment's
public GSE/GSM accession, or select `local_file` and provide the integrated
`Hippo_timecourse_integrated.h5mu` or `iPSC_K562_mixed_integrated.h5mu` produced by
the submission workflow. Automatic downloads are cached in `data/geo/downloads/`.
These large objects remain ignored by Git. Regeneration outputs are stored in
`data/derived/geo_plots/`, never over the published reference tables.
