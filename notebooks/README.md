# Figure reproduction notebooks

## Hippo reconstruction and one-go execution

Open [`preprocessing/Rebuild_Hippo_Timecourse.ipynb`](preprocessing/Rebuild_Hippo_Timecourse.ipynb) on Sherlock and **Run All**. Its settings cell controls the eight-source input directory, Oak output directory, and downstream execution. It rebuilds both Hippo inputs, validates historical counts and program identities, updates the active local configuration while preserving other settings, and executes Figure 3 and Figure 4 in fresh kernels. See [the detailed workflow](../docs/HIPPO_REBUILD.md).

## Main figures

- `main/Figure_01.ipynb`
- `main/Figure_02.ipynb`
- `main/Figure_03.ipynb` — panels b–d and f
- `main/Figure_04.ipynb` — panels c–e

## Supplementary figures

`notebooks/supplementary/Figure_S01.ipynb` through `Figure_S09.ipynb` retain their existing scope and execution status. The Hippo reconstruction does not establish that supplementary panels have been executed or matched to the manuscript.

## Path configuration

The notebooks use the tracked Sherlock defaults in `config/paths.example.yaml`. An optional, untracked `config/paths.local.yaml` supplies only overrides; existing entries are merged over the defaults. `CCPERTURBSEQ_PATHS`, if set, selects an explicit override file. Do not recopy the example YAML over an existing local configuration.

The rebuild notebook writes the two validated Hippo paths into the active override file and backs up the previous file under its run directory. Large H5ADs and generated figures stay outside Git.

## Validation

```bash
python scripts/preflight_main_figure_notebooks.py --check-data --check-runtime
```

To include the reconstruction notebook in the static/runtime audit, add `--include-hippo-rebuild`. Source notebook archives under `recovery/` are provenance records, not Run All entry points. Execution receipts and comparison with manuscript references determine the actual reproduction status.
