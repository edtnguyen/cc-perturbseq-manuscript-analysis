# Public inputs for manuscript figure reproduction

The seven figure notebooks use IGVF deposits, versioned companion data under
`data/reproduction/`, and the two integrated GEO H5MU files under `data/geo/`.
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
motif files; and the complete S4 plotting inputs. It excludes FASTQs, expression
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
directory, including that manifest. Its contents have **not** been published
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

## GEO datasets

The processed-object reader is portable, but the public GEO accessions/download
URLs have not been supplied. For now place the submitted files at:

```text
data/geo/Hippo_timecourse_integrated.h5mu
data/geo/iPSC_K562_mixed_integrated.h5mu
```

They are **not copied or uploaded** by the IGVF companion collector. The Hippo
RNA and cNMF modalities are extracted into `data/derived/hippo/` without refitting.
S2 uses frozen lineage labels and recomputes guide metrics from its sparse guide
modality; it does not recluster or read Split-pipe source files. Both require the
`ccperturbseq.geo.v1` schema produced by the submission workflow.

Missing GEO files produce an explicit missing-input error, never a fallback to
Sherlock. Having portable code is not the same as having all public data
published: the companion bundle and these two H5MU files must be materialized
before all seven notebooks can run.

## Checks and execution

```bash
python scripts/audit_public_inputs.py
python scripts/audit_public_inputs.py --check-data
python scripts/prepare_public_data.py --figure Figure_01 --no-prompt
```

The static audit checks all seven notebook sources and default-path providers.
`--check-data` reports missing companion/GEO files honestly; it does not generate
placeholder data. The preparation command downloads and assembles actual IGVF
inputs. Every notebook performs the same preparation in its setup cell.

The tests use checksum validation, synthetic day partitions, the real table
adapter and original plotting functions. They also execute the supplementary
notebooks in a fresh process with access to institutional paths blocked. They
are portability tests, not a numerical validation of the unpublished study
companion files.
