# Hippo time-course reconstruction and figure execution

## Run on Sherlock

Fetch the updated curation branch in the existing checkout:

```bash
cd /oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction
git fetch origin
git switch curation/figure-reproduction-v1
git pull --ff-only origin curation/figure-reproduction-v1
```

Open [`notebooks/preprocessing/Rebuild_Hippo_Timecourse.ipynb`](../notebooks/preprocessing/Rebuild_Hippo_Timecourse.ipynb)
in the original Sherlock `sc-dl` analysis environment, review its settings cell,
and choose **Run All**. Run within a high-memory compute allocation: the original
concatenation materializes 2,118,711 cells before broad QC. The earlier suggested
256 GB allocation is a starting estimate, not a measured memory requirement.

The notebook performs these stages in order:

1. Import the required packages, check all eight processed sublibrary inputs, and
   check the gene-spectra TPM reference.
2. Rebuild and reopen the combined gene-expression H5AD using the production
   concatenation, broad QC, and sample mapping.
3. Rebuild and reopen the selected `top_specific_200` fixed-reference projection.
4. Validate complete cell/metadata alignment, historical cohort counts, program
   labels, matrix semantics, and the inputs needed by the two Hippo panels.
5. Update only the two `hippo_timecourse` entries in the active local YAML;
   preserve other settings and save the previous exact YAML under the run folder.
6. Run the existing main-figure data/runtime preflight, then execute
   `Figure_03.ipynb` and `Figure_04.ipynb` sequentially in separate fresh kernels
   using the same Python executable as the rebuild notebook.

`RUN_FIGURE_NOTEBOOKS = True` is the default. Set it to `False` to stop after the
validated inputs have been configured. Downstream execution covers Figure 3
panels b–d and f and Figure 4 panels c–e. Figures 1–2 are included in the existing
preflight but are not automatically executed. Supplementary notebooks retain
their prior scope and execution status.

## Paths and settings

The settings cell reads `config/paths.example.yaml`, merged with the optional
`config/paths.local.yaml`, or an explicit override selected by
`CCPERTURBSEQ_PATHS`. Existing local overrides take precedence. The rebuild writes
its validated output paths back into the active override file, so an old local
Hippo scratch path is replaced after successful validation. It never overwrites
the tracked example YAML during execution.

| Setting | Default / role |
|---|---|
| `INPUT_ROOT` | `/scratch/users/edtng/Hippo_timecourse`; parent of the eight processed sublibrary folders |
| `REBUILD_ROOT` | `/oak/stanford/groups/engreitz/Users/tri/Hippo_timecourse_rebuild` |
| `SPECTRA_TPM` | The fixed 030526 reference's `gene_spectra_tpm.k_50.dt_0_2.txt` |
| `OVERWRITE` | `False`; reuse a validated matching output, reject unknown or changed outputs |
| `STRICT_HISTORICAL` | `True`; enforce recorded cohort, sample, reference, and selected-gene counts |
| `COMPRESSION` | `None`, matching the source write; `gzip` or `lzf` is optional |
| `RUN_FIGURE_NOTEBOOKS` | `True`; execute the affected main notebooks after reconstruction |

If the eight processed files have moved, edit `INPUT_ROOT` in the notebook or
`hippo_rebuild.input_root` in the local YAML. The exact relative template is:

```text
Hippo:eTN8_GEX_sub{i}/kb_all_main_raw/counts_filtered/adata.h5ad
```

The eight values of `{i}` are 1 through 8, in that order. Missing inputs are
reported together before cohort loading. The producer starts from these
processed objects; it cannot recover missing quantification/cell-calling outputs
from a notebook alone. Their current existence must be checked on Sherlock.

The reference path is:

```text
/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50.gene_spectra_tpm.k_50.dt_0_2.txt
```

This differs from `ccperturbseq.cnmf_spectra`, which names the consensus spectra
used by program-annotation panels. The overdispersed-gene list and annotation
workbook support other comparison modes or diagnostics in the historical
reprojection notebook; they are not inputs to this selected projection. The
tracked six-covariate workbook is still required for the main figure notebooks.

## Outputs

The canonical H5ADs are written to:

```text
/oak/stanford/groups/engreitz/Users/tri/Hippo_timecourse_rebuild/Hippo_timecourse_combined.h5ad
/oak/stanford/groups/engreitz/Users/tri/Hippo_timecourse_rebuild/perturbseq_cNMF_projection/reprojection_gene_set_comparison/Hippo_timecourse_reprojected_cNMF_compare_gene_sets.cell_x_program.h5ad
```

Each has an adjacent rebuild receipt. The combined receipt ends in
`.h5ad.rebuild_receipt.json`; the projection receipt replaces `.h5ad` with
`.rebuild_receipt.json`. A timestamped `runs/<UTC timestamp>/` directory retains:

- Settings, package versions, input records, stage reports, and validation.
- The prior local path configuration when an update was necessary.
- The full main-figure preflight log.
- `Figure_03.executed.ipynb` and `Figure_04.executed.ipynb`.
- `figure_execution_receipt.json`, including the current Git commit and Python
  executable, execution status, timestamps, and any failure.
- Copies of successful figure output directories under `figures/`.

The ordinary figure-output location still follows the existing `outputs.root`
setting. Successful outputs are additionally copied to the Oak run directory.
Tracked notebooks remain free of execution outputs. These files are not uploaded
to GitHub by the notebook.

## Scientific provenance

### Combined expression query

The uploaded production notebook has SHA-256
`df24b42905e800ffe501ae8e1d2c334c40c684c2306eab36a44dc9f869565a3c`.
Its relevant source cells exactly match the same zero-based cells in the existing
[`recovery/raw_downloads/Hippo_timecourse.ipynb`](../recovery/raw_downloads/Hippo_timecourse.ipynb).
The archived notebook and uploaded notebook have different full-file hashes;
the equality claim concerns the producer source cells listed below.

| Source cell | Transformation |
|---|---|
| 5 | Read eight sublibraries into keys `s1`…`s8`; concatenate with `label="sublibID"`, `index_unique="_"`, `join="outer"`, `merge="same"`, `uns_merge="same"` |
| 14 | Retain `1200 < n_genes_by_counts < 12000`, `1200 < total_counts < 50000`, `pct_counts_mt < 10`, and `biological_sample != "empty"` |
| 20 | Remove unused biological-sample categories after empty-sample exclusion |
| 23 | Strip biological-sample labels, apply the original 54-entry mapping, and make `sample_name` categorical |
| 27 | Write the combined object before later condition-specific filtering |

Historical dimensions are **2,118,711 × 62,757** before broad QC and
**212,891 × 62,757** afterward. The validator also checks all 54 historical
per-sample counts from the production notebook, rather than only an aggregate
cell count. It preserves raw `X`, `ambiguous`/`mature`/`nascent` layers, guide
matrices, gene annotations, and sample/sublibrary metadata. Quantification values
may be fractional: no rounding or integer-count assertion is introduced.

The additional condition-specific mitochondrial/ribosomal filter occurs later
in the source. It must not be applied before writing the combined query.

### Final selected projection

The implementation follows zero-based cells 3, 5, 7, 9, 11, and 13 of
[`recovery/raw_downloads/reproject_cnmf_compare_gene_sets.ipynb`](../recovery/raw_downloads/reproject_cnmf_compare_gene_sets.ipynb).
Those cells also match the archived Sherlock source-text mirror. The mirror
identifies its source notebook SHA-256 as
`7a57ce043047b360c5a0384b13d3cbcd4113674f89690bfbb4621720960ce694`.

The preserved steps are:

- Read the **50 × 17,538** program-by-gene TPM reference and map its programs
  explicitly to `P0`…`P49`.
- For each program, rank genes by its fraction of each gene's total reference
  weight, then by program weight, using the original two-column pandas sort.
- Take each program's first 200 genes, preserve first-occurrence union order,
  then intersect with query symbols. Historically this selects **7,268 genes**.
- Sum duplicate query symbols and normalize using the original full-library
  `total_counts` to 1,000,000 per cell. Renormalize each selected reference
  spectrum exactly as in the source.
- Fit cell usage with fixed spectra, float32 multiplicative updates, Frobenius
  loss, 1,000 maximum iterations, tolerance `1e-5`, and random seed `0`.
- Normalize usage; floor at `1e-6`; renormalize; apply centered log2 CLR.
- Preserve all **212,891** query cells in exactly the same order and add the
  original condition, replicate, treatment, and condition-specific QC flag.
  Historically **99,416** cells pass that flag.

Only the final selected mode is rebuilt. The H5AD contains normalized selected
usage in `X` and these explicitly named matrices:

```python
obsm["X_cnmf_usage_raw_top_specific_200"]
obsm["X_cnmf_usage_top_specific_200"]
obsm["X_cnmf_clr_top_specific_200"]
```

The full historical comparison wrote multiple modes and used its first mode in
`X`; that convenience field therefore differs in this selected-only artifact.
The downstream figures consume the named selected CLR matrix. The older
coordinate-descent all-shared-gene projection and the late Hunter delta-to-control
export are not substitutes for that matrix. No cNMF basis is refitted.

The source selected solve recorded **140 iterations**. Runtime iterations,
convergence warnings, and package versions are saved; the solver settings are
never changed automatically to force a historical iteration count. Narrow tests
compare the new numerical helpers directly with source-notebook functions in the
same environment. Matching those tests does not imply bitwise identity across
unrecorded historical software versions.

## Downstream compatibility and validation

Figure 3f and Figure 4c read gene expression from the raw combined object and
program scores from the selected projection. They require complete cell coverage,
align explicitly by cell ID, verify relevant metadata, and verify program-column
labels. Missing metadata cannot silently disable QC. Every requested condition
must retain cells after the source condition-specific filter and `s8` exclusion.
Marker lists, expression normalization, condition order, and plotting scales are
unchanged.

The new workflow also checks the current figure marker/condition requirements
before replacing local path entries. Projection validation recomputes CLR from
the stored usage and checks full cell order and metadata agreement. H5AD writes
are staged, reopened, and validated before replacing final files.

Completed outputs are reused only when their receipts, input fingerprints,
output fingerprints, and scientific checks agree. Large-file fingerprints record
path, size, and nanosecond modification time; they are not described as content
hashes. The small reference and implementation are hashed where recorded. An
unknown or changed output requires deliberately setting `OVERWRITE = True`.

If execution stops after reconstruction, the validated H5ADs remain available.
Correct the reported downstream problem and rerun the final cell, or open the
affected main notebook directly. A failed figure notebook is saved with the
failed cell, and later figure notebooks do not execute after that failure.

## Dependencies and verification

Use the original analysis environment. The reconstruction uses `anndata`,
`h5py`, `numpy`, `pandas`, `scipy`, `scikit-learn`, and `PyYAML`. Figure validation
and execution also use the repository's plotting dependencies, `nbformat`,
`nbclient`, and `ipykernel`. The notebook imports its required packages before
loading the cohort and reports missing or broken imports together.

The repository preflight supports the new entry point:

```bash
python scripts/preflight_main_figure_notebooks.py --include-hippo-rebuild --check-runtime
```

After reconstruction, the usual Sherlock check remains:

```bash
python scripts/preflight_main_figure_notebooks.py --check-data --check-runtime
```

The focused local suite is:

```bash
PYTHONPATH=src python -m pytest tests/test_hippo_rebuild.py tests/test_hippo_projection.py tests/test_hippo_consumers.py tests/test_hippo_workflow.py
```

It covers production QC boundaries, sample mappings, preserved fractional counts,
layer/guide/barcode alignment, source numerical parity, output round trips,
consumer alignment and expression normalization, safe configuration updates,
fresh-kernel execution, and failure records. Production source-file existence,
the full reconstruction, and manuscript figure comparison must be verified on
Sherlock; local fixture tests do not establish those outcomes.

### Implementation checks completed on 2026-10-29

- The 33 focused tests passed, including the fresh-kernel success/failure tests.
- The main-plus-rebuild preflight passed for five notebooks, including runtime dependency checks and the tracked workbook SHA-256 check. This local run did not use `--check-data` for Sherlock paths.
- All 14 curated notebooks passed notebook-schema validation and code-cell compilation.
- The actual rebuild notebook ran successfully in a fresh kernel on a bounded synthetic eight-sublibrary fixture: 432 cells × 1,305 genes, all 54 sample labels, and a 432 × 50 projected object. Only `STRICT_HISTORICAL` and downstream figure execution were disabled in the temporary test copy. Marker/condition checks and local-config preservation passed. Separate tests exercised downstream notebook execution and both Hippo panel builders.

These checks validate the implementation and its integration. The historical 212,891-cell dataset was not available in the development environment, and full figure execution/reference comparison remains a Sherlock task.
