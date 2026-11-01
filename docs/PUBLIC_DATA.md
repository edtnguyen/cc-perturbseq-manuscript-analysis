# Public inputs for manuscript figure reproduction

The seven figure notebooks use IGVF deposits, versioned companion data under
`data/reproduction/`, including precomputed plotting data for the two GEO experiments.
The integrated GEO H5MU files are optional inputs for regenerating those tables.
They do not read original institutional filesystem locations or old notebooks.
Previously saved notebook outputs are retained; those outputs may contain paths
from the recorded execution. They are not inputs to a new run.

## IGVF downloads and authentication

Each notebook's setup downloads only its required deposits to `data/igvf/` and
verifies the accession, file size and MD5 against
`config/public_data_manifest.json`. Downloads are atomic; partial or
checksum-mismatched files are not consumed. Verified cached downloads are reused.

The selected files are released and normally download anonymously. Authentication
is optional. Supply **portal access keys**, not a portal login password, through
`IGVF_ACCESS_KEY` and `IGVF_SECRET_ACCESS_KEY`; or allow the hidden `getpass`
prompts when a request returns HTTP 401/403. Keys remain in memory, are not put
in command-line arguments, notebook outputs, receipts or committed files, and
are not forwarded to another origin on storage redirects. Authentication cannot
grant access beyond the permissions of the supplied key.

Official API documentation: https://github.com/IGVF-DACC/igvf-python-client

## Repository paths

All paths in `config/paths.example.yaml` are resolved from the repository root,
not from the directory containing the notebook. Legacy `paths.local.yaml` and
`CCPERTURBSEQ_PATHS` are deliberately ignored, so an old Sherlock configuration
cannot silently override public inputs. Use `CCPERTURBSEQ_PUBLIC_PATHS` for an
explicit portable override file. Existing local files are not overwritten.

Day-specific annotated H5AD and cNMF H5MU partitions are assembled in
`data/derived/` in the original, published cell order. Existing UMAP coordinates,
program values and guide-column identities are preserved. cNMF and UMAP are not
refitted. The Figure 1d derived H5MU contains only its required cNMF modality.

The exact separate March usage and consensus TSV representations are retained in
the companion bundle as losslessly compressed files. Their equivalence to H5MU
matrices has not been assumed. The original guide-metadata table is retained as
an analysis-label companion rather than replacing gene symbols with portal
Ensembl labels. Figure 1c uses the deposited six count columns through a verified
one-to-one bridge that preserves the original guide names, order and QC totals.

## One-time maintainer collection (run where the originals exist)

Run from an up-to-date **main** checkout using the existing analysis environment:

```bash
python scripts/collect_publication_data.py
python scripts/collect_publication_data.py --apply
```

The first command is a dry run. The second copies and verifies the exact data.
Source locations are isolated in `config/publication_source_files.json`; this
maintainer-only inventory is never read by a figure notebook. Override that
inventory with `--source-config`. Select a completed S4 result directory using
`--s4-results` and, when needed, `--target-counts` for its matching support table.

The collector includes the genome-wide U-test table; the validation guide-name
and QC bridge; cNMF usage/spectra; program annotations and display metadata;
original expression/cNMF cell order; exact guide metadata; the three selected
motif files; the complete S4 plotting inputs; and precomputed Hippo/mixed-cell plotting tables
generated from the two submitted integrated H5MUs. It excludes FASTQs, expression
matrices already deposited to IGVF, old notebooks, generated figure files and
unrelated analyses. It does not rerun CRT. If S4 target-cell counts need rebuilding,
it uses the original recorded inputs on the source machine once, then bundles
the resulting numerical counts.

The motif files are the active `enhancer_enrichment_results_min10.tsv`, current
`cluster_to_tf_motif_id_mapping.tsv`, and **standard** `modisco_compiled.h5` (not
the profile compendium). Runtime validation pins their bytes, not an institution
path.

Copying is verified by SHA-256; gzip copies are also checked after decompression.
No files are staged to Git automatically. A 95 MiB per-file guard prevents an
accidental oversized Git publication. `--overwrite` backs up an existing bundle.
The script emits `data/reproduction/bundle_manifest.json` only after all required
inputs are collected successfully. Publish the entire `data/reproduction/`
directory, including that manifest and the independently checksummed
`geo_plots/hippo_timecourse/` and `geo_plots/ipsc_k562_mix/` bundles. Their
matrices/metrics are regenerated with the existing numerical producers, not
from figure images. Its contents have **not** been published
merely by installing this code.

## S4

S4 reads the completed model's Q–Q data, replicate beta matrices and relocated
manifest, pooled p-values, and target-cell support counts from
`data/reproduction/s4/k_50.dt_0_2/`. Beta file references are relative; absolute
references are rejected. The >=20 support rule and pooled full-day q<=0.10 rule
are unchanged, with empirical Q–Q values for S4b and skew-preferred pooled
selection for S4c. The support table is read directly, without probing raw source
paths.

Exact raw-CRT recomputation is a separate, optional CLI workflow requiring the
original January-model inputs and an explicit programDE installation. These
inputs are not deposited in the supplied March IGVF sets. The figure notebook
therefore does not silently swap models or attempt raw recomputation.

## GEO experiments: precomputed by default

Figures 3f, 4c and S2i/j/l/m/n default to `DATA_SOURCE="precomputed"`.
They read the following versioned directories, with no GEO download or processed
single-cell object required:

```text
data/reproduction/geo_plots/hippo_timecourse/
data/reproduction/geo_plots/ipsc_k562_mix/
```

The Hippo bundle includes 12 compressed tables: gene/program summaries and the
actual dot-color and dot-size matrices for each of Figure 3f and Figure 4c. The
Scanpy renderer receives these matrices directly; it does not re-normalize the
already scaled gene colors. The mixed-cell bundle includes all QC-passing cells'
metrics and frozen lineage labels, all four saturation curves, cluster annotations,
and median diagnostics. Zero-guide cells are retained. No old notebook or rendered
figure is used as an input.

### Prepare the precomputed tables for publication

The normal collector includes both bundles. It reads the exact integrated H5MU
files created for GEO, with defaults recorded only in the maintainer source inventory.
The processed objects themselves are not copied into Git. To specify their locations:

```bash
python scripts/collect_publication_data.py --apply --hippo-h5mu /path/to/Hippo_timecourse_integrated.h5mu --mixed-h5mu /path/to/iPSC_K562_mixed_integrated.h5mu
```

To add only these two bundles without recopying the other companion data:

```bash
python scripts/collect_publication_data.py --geo-plots-only
python scripts/collect_publication_data.py --geo-plots-only --apply
```

The first command is a dry run. Use the two path options above when the submitted
files are in another location. Use `--overwrite` to back up existing plot bundles
before replacing them. Copying/regeneration is not a Git push; after validation,
commit and push `data/reproduction/geo_plots/` with its two `manifest.json` files.
No integrated H5MU or raw sequencing file is uploaded by this command.

### Optional regeneration from GEO

In the relevant notebook settings:

```python
DATA_SOURCE = "geo"
GEO_ACCESSION = "GSE123456"  # Replace with this experiment's public GSE or GSM
GEO_FILENAME = None         # Or the exact processed supplementary filename
GEO_EXPECTED_SHA256 = None  # Optional SHA-256 of the decompressed H5MU
GEO_LOCAL_FILE = None
```

The example accession above is illustrative, not an accession assigned to this
study. Public study accessions have not yet been supplied in this repository.
The resolver reads NCBI SOFT metadata, selects the integrated processed file,
and downloads it over HTTPS. Series- and sample-level attachments are supported.
Ambiguous matches require an explicit filename. A `.h5mu.gz` is decompressed
before use. Raw TAR/FASTQ downloads are not substituted for processed data.
Released GEO records require no authentication. Private/unreleased records give
an actionable error; no NCBI password or reviewer token is requested or stored.

Before release or when an integrated file was obtained separately:

```python
DATA_SOURCE = "local_file"
GEO_LOCAL_FILE = "data/geo/Hippo_timecourse_integrated.h5mu"
GEO_ACCESSION = ""
```

For S2 use `data/geo/iPSC_K562_mixed_integrated.h5mu` instead. Relative paths resolve
from the repository root. This study's `ccperturbseq.geo.v1` integrated H5MU schema
is required. An arbitrary RNA-only H5AD cannot replace the cNMF/guide modalities,
frozen labels, or sequencing-saturation tables; unsupported objects are rejected.

Regeneration recalculates panel inputs from the submitted modalities and verifies
embedded matrix checksums and reference tables when present. Hippo uses the
unchanged condition-specific QC, s8 exclusion, panel-specific RNA normalization,
and fixed-reference CLR values. S2 recomputes sparse guide metrics without
reclustering. Regenerated tables go to a content-addressed directory under
`data/derived/geo_plots/`; published tables are never overwritten by notebook
execution. `comparison_to_published.json` records whether all tables agree when
a reference bundle is present. Differences are reported, not silently ignored.
The plots in that run use the regenerated tables.

A command-line equivalent is available:

```bash
python scripts/prepare_geo_plot_data.py --dataset hippo_timecourse --geo-accession GSE123456
python scripts/prepare_geo_plot_data.py --dataset ipsc_k562_mix --local-file data/geo/iPSC_K562_mixed_integrated.h5mu
```

NCBI interface documentation:
https://www.ncbi.nlm.nih.gov/geo/info/download.html
https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html

## Checks and execution

```bash
python scripts/audit_public_inputs.py
python scripts/audit_public_inputs.py --check-data
python scripts/prepare_public_data.py --figure Figure_01 --no-prompt
```

The static audit checks all seven notebook sources and default-path providers.
`--check-data` verifies companion and precomputed plotting bundles; it does not
require the optional GEO H5MUs and reports missing/invalid plotting tables honestly; it does not generate
placeholder data. The preparation command downloads and assembles actual IGVF
inputs. Every notebook performs the same preparation in its setup cell.

The tests use checksum validation, synthetic day partitions, the real table
adapter and original plotting functions. They also execute the supplementary
notebooks in a fresh process with access to institutional paths blocked. They
are portability tests, not a run on the unpublished study companion files.
Additional tests compare original Hippo producers with table-based rendering,
verify optional GEO/local regeneration, and execute the affected plotting cells
without institutional paths. Actual accession-specific downloads require the
public study accession; they cannot be certified before it is available.

## Small intermediate files, without creating an integrated H5MU

Only plotting tables are published for the GEO experiments. The integrated H5MU
files belong on GEO, not in Git. The default figure notebooks already read
`data/reproduction/geo_plots/` and accept the same compact tables described here.

Run on the source machine, using the environment that runs the accepted figures:

```bash
python scripts/export_compact_plot_data.py
python scripts/export_compact_plot_data.py --apply
```

The first command is a dry run. By default the exporter runs the unchanged
Hippo numerical producers on the existing query/projection H5AD pair, and reads
S2's saved `Figure_S02_all_cell_metrics.csv.gz`, transcript saturation table,
cluster annotations and median diagnostics. It does not build an H5MU, rerun
cNMF, recluster S2 cells, or read a rendered figure. Hippo's original numerical
producers still need the same memory as when generating the accepted panels.

Only the maintainer inventory `config/publication_source_files.json` contains
source-machine paths. They can be overridden without editing code:

```bash
python scripts/export_compact_plot_data.py --apply \
  --hippo-query /path/to/Hippo_timecourse_combined.h5ad \
  --hippo-projection /path/to/projected.cell_x_program.h5ad \
  --mixed-figure-dir /path/to/Figure_S02
```

`--dataset hippo_timecourse` or `--dataset ipsc_k562_mix` exports just one dataset.
Missing complete cell metrics are not substituted with a filtered panel table.
Use `--overwrite` only to back up and replace existing plot bundles; unrelated
data are preserved.

For already-created GEO objects, an alternative reads only small saved HDF5
annotation/reference-table groups. It does not load count matrices or write any
intermediate H5AD/H5MU:

```bash
python scripts/export_compact_plot_data.py --apply \
  --hippo-h5mu /path/to/Hippo_timecourse_integrated.h5mu \
  --mixed-h5mu /path/to/iPSC_K562_mixed_integrated.h5mu
```

This is **saved-table extraction**, not independent verification of the count
matrices. Its source digest covers the selected plotting tables, not the entire
H5MU. The existing optional GEO regeneration mode still recomputes the metrics
and checks embedded matrix/reference data. It remains available in the notebooks.

The outputs are 12 gzip-compressed Hippo tables, four gzip-compressed S2 tables,
and a `manifest.json` per dataset. Hippo stores summaries and exact dot-color and
dot-size matrices. S2 retains per-cell values (including zero-guide cells) for
its violin, box and mixing plots; it does not replace distributions with averages.
Every table is read back and checked. The exporter prints total bytes and the
largest table for each dataset. Actual study sizes are measured at export time,
not estimated from synthetic tests.

The existing full companion collector accepts the same route:

```bash
python scripts/collect_publication_data.py --compact-geo-inputs --apply
```

To collect only the two plot bundles with that interface, add `--geo-plots-only`.
The compact exporter does not stage or push files. Publish only:

```bash
git add data/reproduction/geo_plots
git diff --cached --stat
git commit -m "Add compact Hippo and mixed-cell plotting tables"
git push origin HEAD:main
```

Run `python scripts/audit_public_inputs.py --check-data` after the remaining
companion data are present. No full GEO object is required for default plotting.
