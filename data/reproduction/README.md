# Figure reproduction companion data

Run `scripts/collect_publication_data.py --apply` in the maintainer's environment
to collect and verify the exact non-IGVF inputs, including motif files, S4 numerical
results, and plotting tables for the two GEO experiments. `--geo-plots-only`
updates just the Hippo and mixed-cell plotting tables.

Publish the data files and manifests together. The main `bundle_manifest.json`
checksums the screen/program/motif/S4 companions. Each GEO plotting bundle under
`geo_plots/hippo_timecourse/` and `geo_plots/ipsc_k562_mix/` has its own manifest.
Their default use requires no original matrix, experiment folder, or GEO login.
No placeholder scientific tables are distributed. Installing the code alone
does not populate this directory. See `docs/PUBLIC_DATA.md`.
