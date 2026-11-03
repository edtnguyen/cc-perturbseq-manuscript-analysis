# CC-Perturb-seq manuscript figure reproduction

This public `main` branch contains the computational notebooks and supporting code used to reproduce the selected manuscript figures.

## Figure notebooks

Main figures:

- `notebooks/main/Figure_01.ipynb`
- `notebooks/main/Figure_02.ipynb`
- `notebooks/main/Figure_03.ipynb`
- `notebooks/main/Figure_04.ipynb`

Supplementary figures:

- `notebooks/supplementary/Figure_S02.ipynb`
- `notebooks/supplementary/Figure_S04.ipynb`
- `notebooks/supplementary/Figure_S06.ipynb`

## Supporting code

- `src/cc_perturbseq/` — panel calculations, plotting functions, data readers, validation, CRT and post-CRT utilities
- `scripts/` — repository and figure preflight utilities
- `config/` — figure manifests, default input paths, and public-environment specifications
- `data/figure_02/` — bundled Figure 2 workbook and documentation
- `docs/` — input provenance, source audits, motif-input records, and figure-reproduction notes

## Public data and execution

Run notebooks from the repository checkout using the Python 3.11 environment
specified in `config/geo-reproduction-requirements.txt`. Their setup cells download
the required checksum-pinned IGVF files automatically, anonymously by default.
Optional portal access-key authentication is supported without storing secrets.

The notebooks use repository-relative paths from `config/paths.example.yaml`.
Legacy institutional `paths.local.yaml` overrides are not used. Large downloaded
objects are cached under `data/igvf/`; derived representations are written under
`data/derived/`; figures and source tables are written under `results/`.

Additional manuscript inputs are provided through `data/reproduction/`, including
the U-test table, exact program representations and annotations, motif files, and
completed S4 numerical results. A one-time maintainer collector prepares and
checks these files. Hippo and mixed-cell panels use compact precomputed tables
by default. Optional regeneration downloads the submitted integrated H5MUs from
GEO: **GSE352939** (Hippo) and **GSE352893** (iPSC-K562 mix).

**Publication readiness:** default plotting requires the companion data and
precomputed tables, not the large GEO objects. Optional GEO regeneration requires
publicly accessible processed supplementary files. An assigned accession alone
does not guarantee that downloads are released. Missing inputs produce explicit
errors; no workflow falls back to institutional filesystem paths.

```bash
python scripts/audit_public_inputs.py
python scripts/audit_public_inputs.py --check-data
python scripts/prepare_public_data.py --figure Figure_01
```

The first command checks all seven notebook sources. The second checks actual
companion and precomputed plotting-table checksums; it reports incomplete data
when files are absent. The third prepares inputs for one figure (repeat
`--figure` for additional figures).

See [Public data and collection instructions](docs/PUBLIC_DATA.md) for the
maintainer collection command, authentication, dataset accession manifest,
S4 model-version constraints, and validation scope. Existing scientific
calculations, plotting styles, and saved notebook figures are retained.

## Hippo and mixed-cell plotting data

Figures 3f, 4c and S2 default to the checksummed tables under
`data/reproduction/geo_plots/`. No GEO download or original single-cell file is
needed for this route. The optional regeneration settings are already configured:

| Experiment | GEO accession | Notebook / panels |
|---|---|---|
| Parse v3 CC-Perturb-seq, iPSC-K562 mix | [GSE352893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352893) | `Figure_S02.ipynb`, panels i, j, l, m, n |
| Small molecule Hippo scRNA-Seq | [GSE352939](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352939) | `Figure_03.ipynb` panel f and `Figure_04.ipynb` panel c |

To download the submitted object and regenerate plotting inputs, change only:

```python
DATA_SOURCE = "geo"
```

The notebook retains the appropriate `GEO_ACCESSION` above. `GEO_FILENAME=None`
automatically discovers the processed supplementary file; set its exact filename
only if multiple matching files are reported. Released records need no GEO login.
An unreleased record requires `DATA_SOURCE="local_file"` and `GEO_LOCAL_FILE`
pointing to the submitted integrated H5MU until the files become public.

Downloads are cached in `data/geo/downloads/`. Regenerated tables go to
`data/derived/geo_plots/` and are compared with published tables when available;
the published precomputed bundle is never overwritten. Regeneration starts from
the submitted processed RNA plus guide/program modalities, not FASTQ files; it
does not refit cNMF, recluster cells, or reprocess sequencing reads.

Command-line equivalents:

```bash
python scripts/prepare_geo_plot_data.py --dataset ipsc_k562_mix --geo-accession GSE352893
python scripts/prepare_geo_plot_data.py --dataset hippo_timecourse --geo-accession GSE352939
```

Maintainers can run `python scripts/collect_publication_data.py --geo-plots-only --apply`
to generate just these two plotting bundles, then commit and push them. This does
not upload the H5MU files. See `docs/PUBLIC_DATA.md` for input-path options and
publication instructions. The actual data are not present until that collection
step has been run and its outputs published.

Small plotting intermediates for the two GEO experiments can be prepared with
`python scripts/export_compact_plot_data.py --apply`; no new integrated H5MU is
needed. Only the compressed tables in `data/reproduction/geo_plots/` are committed.
The optional GEO regeneration route remains available in the notebooks.
