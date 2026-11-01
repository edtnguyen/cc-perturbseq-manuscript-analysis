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
checks these files. The Hippo and mixed-cell integrated GEO H5MUs belong in
`data/geo/`; their download URLs are not configured yet.

**Publication readiness:** adding the download code does not populate the
companion bundle or the GEO files. Until the actual companion data is committed
and the GEO files are available, affected notebooks stop with explicit missing-
input errors rather than using an inaccessible institutional fallback.

```bash
python scripts/audit_public_inputs.py
python scripts/audit_public_inputs.py --check-data
python scripts/prepare_public_data.py --figure Figure_01
```

The first command checks all seven notebook sources. The second checks actual
companion-file checksums and GEO-file availability; it reports incomplete data
when files are absent. The third prepares inputs for one figure (repeat
`--figure` for additional figures).

See [Public data and collection instructions](docs/PUBLIC_DATA.md) for the
maintainer collection command, authentication, dataset accession manifest,
S4 model-version constraints, and validation scope. Existing scientific
calculations, plotting styles, and saved notebook figures are retained.
