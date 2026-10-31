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

No other analysis notebooks are included on `main`.

## Supporting code

- `src/cc_perturbseq/` — panel calculations, plotting functions, data readers, validation, CRT and post-CRT utilities
- `scripts/` — repository and figure preflight utilities
- `config/` — figure manifests, default input paths, and public-environment specifications
- `data/figure_02/` — bundled Figure 2 workbook and documentation
- `docs/` — input provenance, source audits, motif-input records, and figure-reproduction notes

## Running the notebooks

Run notebooks from the repository checkout so they can locate `src/cc_perturbseq`.

Machine-specific paths should be placed in an untracked `config/paths.local.yaml`, using `config/paths.example.yaml` as the starting point. Local values are merged over the tracked defaults.

Large single-cell objects and other primary analysis inputs are not stored in Git. Public processed-data objects and raw sequencing files are intended to be distributed through GEO; institutional paths in the example configuration document the original analysis inputs.

Before running the figure notebooks, the main static preflight can be executed with:

```bash
python scripts/preflight_main_figure_notebooks.py --check-data --check-runtime
```
