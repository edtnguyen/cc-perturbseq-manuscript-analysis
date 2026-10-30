# Figure reproduction notebooks

## Hippo reconstruction and one-go execution

Open [`preprocessing/Rebuild_Hippo_Timecourse.ipynb`](preprocessing/Rebuild_Hippo_Timecourse.ipynb) on Sherlock and **Run All**. Its settings cell controls the eight-source input directory, Oak output directory, and downstream execution. It rebuilds both Hippo inputs, validates them, updates the active local configuration while preserving other settings, and executes Figure 3 and Figure 4 in fresh kernels. See [the detailed workflow](../docs/HIPPO_REBUILD.md).

## Main figures

- `main/Figure_01.ipynb`
- `main/Figure_02.ipynb`
- `main/Figure_03.ipynb` — panels b–d and f
- `main/Figure_04.ipynb` — panels c–e

## Supplementary figures retained

Only the supplementary notebooks currently needed for manuscript reproduction are retained:

- `supplementary/Figure_S02.ipynb` — mixed hiPSC–K562 Parse v3 panels i, j, l, m and n
- `supplementary/Figure_S04.ipynb` — program CRT calibration and biological-replicate effect concordance, panels b and c
- `supplementary/Figure_S06.ipynb` — integrated D0–D3 program annotation

The notebooks for Figures S1, S3, S5, S7, S8 and S9 are intentionally not part of this curated branch.

## Path configuration

The notebooks use the tracked Sherlock defaults in `config/paths.example.yaml`. An optional, untracked `config/paths.local.yaml` supplies only overrides; existing entries are merged over the defaults. `CCPERTURBSEQ_PATHS`, if set, selects an explicit override file. Do not recopy the example YAML over an existing local configuration.

The rebuild notebook writes the two validated Hippo paths into the active override file and backs up the previous file under its run directory. Large H5ADs and generated figures stay outside Git.

## Validation

```bash
python scripts/preflight_main_figure_notebooks.py --check-data --check-runtime
```

To include the reconstruction notebook in the static/runtime audit, add `--include-hippo-rebuild`. Source notebook archives under `recovery/` are provenance records, not Run All entry points. Execution receipts and comparison with manuscript references determine the actual reproduction status.
