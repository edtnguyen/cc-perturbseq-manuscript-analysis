# Figure reproduction notebooks

This directory contains one notebook per requested manuscript figure.

## Main figures

- `main/Figure_01.ipynb`
- `main/Figure_02.ipynb`
- `main/Figure_03_panels_b_d.ipynb`
- `main/Figure_04_panels_d_e.ipynb`

## Supplementary figures

- `supplementary/Figure_S01.ipynb`
- `supplementary/Figure_S02.ipynb`
- `supplementary/Figure_S03.ipynb`
- `supplementary/Figure_S04.ipynb`
- `supplementary/Figure_S05.ipynb`
- `supplementary/Figure_S06.ipynb`
- `supplementary/Figure_S07.ipynb`
- `supplementary/Figure_S08.ipynb`
- `supplementary/Figure_S09.ipynb`

The notebooks load their panel registries from `config/figure_notebooks.yaml` and machine-specific input paths from `config/paths.local.yaml`.

## Setup on Sherlock

```bash
cp config/paths.example.yaml config/paths.local.yaml
$EDITOR config/paths.local.yaml
```

Keep `config/paths.local.yaml` untracked. Large H5AD, H5, RDS and intermediate result objects remain outside Git.

## Current status

These are structured figure drivers. They contain input audits, panel definitions, output directories and source-provenance leads. The canonical historical analysis code still needs to be selected, ported and validated panel by panel. A panel is complete only when it runs from a fresh kernel and exports vector/raster outputs, plotted source data and a run receipt.
