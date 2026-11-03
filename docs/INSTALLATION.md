# Installation and notebook execution

## Scope

`environment.yml` installs a new **Python 3.11**, CPU-only environment for all
seven published notebooks: Figures 1–4, S2, S4 and S6. It includes numerical and
statistical libraries; AnnData/MuData/HDF5 readers; Scanpy, plotnine, Matplotlib
and motif-logo rendering; Excel/YAML readers; JupyterLab and its Python kernel;
and the repository's test dependencies. Scientific package pins are retained
from `config/geo-reproduction-requirements.txt`; tests keep the two lists in sync.

No R, rpy2, CUDA, GPU, separate SCEPTRE environment, or external programDE checkout
is needed to render these notebooks from the published inputs. S4 plots completed
CRT results. Optional **raw CRT estimation**, FASTQ preprocessing and model fitting
are outside this environment's scope and require their own software/data setup.
The optional GEO route recomputes plotting tables from the submitted processed
H5MU, not from FASTQs.

The primary installation target is Linux x86-64 with Conda/Miniforge. Other
platforms may have different wheel availability; do not silently replace pinned
scientific versions when an installation cannot be solved.

## Create a separate Conda environment

From the repository root (the directory containing `environment.yml`):

```bash
conda env create --file environment.yml
conda activate ccperturbseq-figures
python -m ipykernel install --user --name ccperturbseq-figures --display-name "Python (ccperturbseq figures)"
python scripts/check_environment.py
```

The final check must report `ENVIRONMENT PASS`. It checks Python, requirement
versions, imports of installed libraries and repository modules, and `pip check`.
It does not download study datasets. Do not run these installs over the existing
`sc-dl` or `sceptre_env` environments. If `ccperturbseq-figures` already exists,
inspect it first, or create a distinct name with
`conda env create --name ccperturbseq-figures-new --file environment.yml` and
activate/register that name instead.

Start local JupyterLab with `jupyter lab`. In the notebook interface choose
**Python (ccperturbseq figures)**, restart the kernel, and run from the first cell.
Installing a package in a terminal's environment does not switch an already
running notebook kernel.

## Sherlock

Use the updated main-based checkout that contains your published companion data.
For example, in your existing test checkout:

```bash
cd /oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-test-20261103_050510
conda env create --prefix /oak/stanford/groups/engreitz/Users/tri/envs/ccperturbseq-figures --file environment.yml
conda activate /oak/stanford/groups/engreitz/Users/tri/envs/ccperturbseq-figures
python -m ipykernel install --user --name ccperturbseq-figures --display-name "Python (ccperturbseq figures)"
python scripts/check_environment.py
```

Change only the checkout path if you are using a different checkout. The prefix
above keeps the environment separate from `sc-dl` and outside purgeable scratch.
It is an installation example, not a hard-coded runtime dependency. Use a
compute-node allocation or Sherlock's allocated Jupyter session for scientific
execution, rather than running the notebooks on a login node. The code does not
launch a Jupyter server on the login node or modify your existing kernels.

When the terminal and notebook disagree, run this in the notebook:

```python
import sys
print(sys.executable)
```

It should point into `ccperturbseq-figures`. After pulling module updates, restart
the notebook kernel so Python does not retain previously imported functions.

## Without Conda

Use an existing Python 3.11 installation to create an isolated virtual environment:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r config/geo-reproduction-requirements.txt
python -m ipykernel install --user --name ccperturbseq-figures --display-name "Python (ccperturbseq figures)"
python scripts/check_environment.py
```

This uses the same Python requirements as the YAML. `.venv/` is ignored by Git.
The exact transitive environment can be recorded after a successful install:

```bash
mkdir -p output/environment
python -m pip freeze > output/environment/pip-freeze.txt
conda list --explicit > output/environment/conda-explicit.txt
```

Run the last command only for a Conda environment. It is platform-specific and
is an audit record, not a cross-platform substitute for `environment.yml`.

## Check data separately

```bash
python scripts/audit_public_inputs.py
python scripts/audit_public_inputs.py --check-data
```

The first command checks notebook sources and paths. The second verifies the
actual companion/plotting bundles. A package install cannot fix missing data or
an unreleased GEO accession. The maintainer collector is not needed by readers
when the published `data/reproduction/` bundle is already in their clone.
See [Public inputs](PUBLIC_DATA.md) for downloads and checksums.

## Output locations

Every figure notebook writes to the same checkout's `output/<figure_id>/`, with
panel subdirectories such as `output/Figure_01/panel_b/`. Numerical tables,
PDF/SVG/PNG exports and run receipts use this root. Scanpy's default figure
and generated cache directories follow the same convention. New processed-data
intermediates and optional GEO-regenerated tables use `output/cache/derived/`.

This is independent of the launch directory, `$SCRATCH` and legacy
`paths.local.yaml`. Public path overrides cannot move generated files outside
`output/`, and a symlink redirecting `output/` outside the checkout is rejected.
The setup cell prints the exact output directory. Downloads remain cached in
`data/igvf/` and `data/geo/downloads/`; these are input files, not plotted outputs.
The original deposited or bundled files are not modified.

Old `results/`, `data/derived/`, and scratch outputs are not moved, deleted or
reused as figure outputs. Existing cache files can remain; the new intermediate
paths are rebuilt as needed. Saved notebook outputs can display historical paths
until the notebook is rerun and saved. `output/` is ignored by Git; notebooks
still retain their saved embedded figures. Matplotlib's existing editable-text
PDF/SVG settings are unchanged. Arial is not distributed with the repository;
Matplotlib can use its bundled DejaVu Sans fallback when Arial is unavailable.

## Reference documentation

- Conda environment creation: https://docs.conda.io/projects/conda/en/stable/commands/env/create.html
- Environment YAML specification: https://conda.org/learn/specifications/exchange/environment-yml/
- IPython kernel registration: https://ipython.readthedocs.io/en/stable/install/kernel_install.html
- Sherlock job guidance: https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/
