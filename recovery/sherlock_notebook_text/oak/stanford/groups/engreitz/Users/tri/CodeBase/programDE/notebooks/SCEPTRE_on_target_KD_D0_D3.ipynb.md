# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/SCEPTRE_on_target_KD_D0_D3.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/SCEPTRE_on_target_KD_D0_D3.ipynb`
- SHA-256: `3132ce1a8c4ca7c11cf9236c7b8474a29436820edec5038f5822361974aa44a3`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# SCEPTRE on-target CRISPRi knockdown across D0–D3 — corrected SCEPTRE environment

This notebook reproduces the on-target SCEPTRE analysis design in `ipsc_ec_sceptre.ipynb` and makes one overlaid plotnine figure for **NF2, LATS2, AMOTL2 and TAOK1** across **D0, D1, D2 and D3**.

## Environment design

Use the following two Sherlock environments:

- **Notebook kernel and Python/plotnine:** `/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl`
- **R and SCEPTRE:** `/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/bin/Rscript`

The notebook runs with the `sc-dl` Python kernel and launches the SCEPTRE stage through the **direct `Rscript` executable** above. It no longer invokes `conda run`, and it does not use the `rpy2` environment.

The installed SCEPTRE in `sceptre_env` is version 0.10.2, whereas official `fold_change` and `se_fold_change` output begins in version 0.10.3. To avoid modifying or destabilizing the existing environment, the notebook installs the exact **SCEPTRE v0.10.3 tag** into this private R library when needed:

```text
/oak/stanford/groups/engreitz/Users/tri/Rlibs/sceptre_0.10.3
```

The subprocess receives that directory through `R_LIBS_USER`, so SCEPTRE v0.10.3 is loaded ahead of the environment's v0.10.2 copy. The install is skipped on subsequent runs after the exact version is verified. Set `INSTALL_PRIVATE_SCEPTRE_IF_NEEDED=False` only after v0.10.3 is already accessible to this Rscript.

## Figure encoding

- All four time points are displayed on **one axis**, horizontally offset within each target gene.
- **Shape** identifies time point: D0 circle, D1 square, D2 triangle and D3 diamond.
- **Open symbols at 1** show the control reference.
- **Filled symbols** show the SCEPTRE fold-change estimate for cells assigned a targeting guide.
- **Orange** indicates within-day BH-adjusted SCEPTRE `q ≤ 0.05`; black indicates `q > 0.05`.
- Grey bars are approximate 95% model-based intervals derived from SCEPTRE v0.10.3 `fold_change` and `se_fold_change` outputs.

## Analysis retained from the original notebook

For each day separately, the notebook uses the original raw response-count and gRNA-count transfer files, imports the data as high MOI, tests the on-target target–response pairs, uses two-sided tests and complement controls, assigns gRNAs with the mixture method, applies the original mitochondrial QC threshold, runs the calibration check, and then runs the discovery analysis. It writes new result files and never overwrites the historical `TargetKD.csv` files.

The SCEPTRE p-values remain conditional-resampling-based. The displayed intervals are model-based standard-error intervals and are not obtained by inverting the resampling test; the two uncertainty summaries need not agree exactly.

## Code cell 1

```python
import os
from pathlib import Path

class cfg:
    TARGETS = ["NF2", "LATS2", "AMOTL2", "TAOK1"]
    DAYS = ["d0", "d1", "d2", "d3"]

    # Correct Sherlock environments.
    SCEPTRE_ENV = Path("/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env")
    SCEPTRE_RSCRIPT = SCEPTRE_ENV / "bin/Rscript"
    SC_DL_ENV = Path("/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl")
    REQUIRE_SC_DL_KERNEL = True

    # Use an exact private SCEPTRE 0.10.3 installation without replacing 0.10.2
    # inside sceptre_env. The first run installs it when absent; later runs verify/reuse it.
    SCEPTRE_PRIVATE_R_LIB = Path(
        "/oak/stanford/groups/engreitz/Users/tri/Rlibs/sceptre_0.10.3"
    )
    SCEPTRE_REQUIRED_VERSION = "0.10.3"
    SCEPTRE_GITHUB_REPO = "Katsevich-Lab/sceptre"
    SCEPTRE_GITHUB_REF = "v0.10.3"
    INSTALL_PRIVATE_SCEPTRE_IF_NEEDED = True
    CRAN_REPOSITORY = "https://cloud.r-project.org"

    # Paths copied from the original SCEPTRE notebook.
    NORMALIZED_DATA_DIR = Path(
        "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
        "2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data"
    )
    BASEDIR = Path(
        "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks"
    )
    GUIDE_METADATA_TSV = BASEDIR / "300genes_guide_metadata_v43.tsv"
    OLD_COMBINED_CSV = BASEDIR / "combined_TargetKD.csv"
    OUTDIR = BASEDIR / "sceptre_on_target_KD_D0_D3_overlay_plotnine_sceptre_env_v2"

    # Input preparation and SCEPTRE execution.
    FORCE_REBUILD_TRANSFER = False
    FORCE_RERUN_SCEPTRE = False
    RUN_CALIBRATION_CHECK = True
    SAVE_SCEPTRE_RDS = False
    R_PARALLEL = False  # Matches the original notebook and avoids cluster oversubscription.
    RANDOM_SEED = 20260922
    RESULT_FILENAME = "TargetKD_v0.10.3_with_SE.csv"

    # Statistical display.
    Q_THRESHOLD = 0.05
    CI_LEVEL = 0.95
    DAY_OFFSETS = {"d0": -0.27, "d1": -0.09, "d2": 0.09, "d3": 0.27}
    DAY_LABELS = {"d0": "D0", "d1": "D1", "d2": "D2", "d3": "D3"}
    DAY_SHAPES = {"D0": "o", "D1": "s", "D2": "^", "D3": "D"}

    # Plot settings.
    FIGSIZE = (8.7, 3.6)
    EXPORT_DPI = 600
    FONT_SIZE = 9
    POINT_SIZE = 2.8
    POINT_STROKE = 0.8
    CI_LINEWIDTH = 0.65
    CI_CAP_WIDTH = 0.07
    SIGNIFICANT_COLOR = "#D95F02"
    NS_COLOR = "#222222"
    CONTROL_EDGE = "#777777"
    CI_COLOR = "#777777"
    REFERENCE_LINE_COLOR = "#B0B0B0"

    TITLE = "On-target CRISPRi knockdown across D0–D3"
    Y_LABEL = "Target-gene expression relative to control"
    X_LABEL = "Target gene"
    PANEL_LETTER = ""

for day in cfg.DAYS:
    if day not in cfg.DAY_OFFSETS or day not in cfg.DAY_LABELS:
        raise ValueError(f"Missing plotting configuration for {day}.")
if len(set(cfg.TARGETS)) != len(cfg.TARGETS):
    raise ValueError("TARGETS must be unique.")
```

## Markdown cell 2

## 1. Validate the Python kernel and direct SCEPTRE Rscript

Run this notebook with the **`sc-dl` Python kernel**. The plotting and transfer stages use Python packages from:

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl
```

The SCEPTRE stage is a separate subprocess launched directly through:

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/bin/Rscript
```

The notebook does not use `conda run`, does not use `$OAK/Users/tri/envs/rpy2`, and does not require R cell magic.

If `sc-dl` is not registered as a Jupyter kernel, run this once in a Sherlock terminal:

```bash
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python -m ipykernel install --user \
  --name sc-dl --display-name "Python (sc-dl)"
```

On first execution, the notebook checks the SCEPTRE version visible to the direct Rscript. If exact v0.10.3 is not available, it installs the official `v0.10.3` GitHub tag into the private library configured above. This does **not** replace the v0.10.2 package inside `sceptre_env`. A new R process is used for the analysis, so no Jupyter-kernel restart is required after installation.

## Code cell 3

```python
import csv
import gc
import hashlib
import importlib.metadata
import json
import os
import re
import subprocess
import sys
import textwrap
import warnings
import zipfile
from datetime import datetime, timezone

import matplotlib as mpl
from matplotlib import font_manager
import muon as mu
import numpy as np
import pandas as pd
import scipy
import scipy.io as sio
from scipy import sparse
from scipy.stats import norm
from IPython.display import display, Markdown

try:
    import plotnine as p9
    from plotnine import (
        ggplot, aes, geom_point, geom_errorbar, geom_hline,
        scale_fill_manual, scale_color_manual, scale_shape_manual,
        scale_x_continuous, scale_y_continuous, coord_cartesian,
        labs, theme_classic, theme, element_text, element_blank,
        element_line, guides, guide_legend,
    )
except ImportError as exc:
    raise ImportError(
        "plotnine is required in the sc-dl kernel. Select the kernel backed by "
        f"{cfg.SC_DL_ENV / 'bin/python'} and restart the notebook."
    ) from exc

OUTDIR = Path(cfg.OUTDIR)
OUTDIR.mkdir(parents=True, exist_ok=True)


def _resolved(pathlike):
    return Path(pathlike).expanduser().resolve()


current_prefix = _resolved(sys.prefix)
expected_prefix = _resolved(cfg.SC_DL_ENV)
if cfg.REQUIRE_SC_DL_KERNEL and current_prefix != expected_prefix:
    raise RuntimeError(
        "This notebook must run with the sc-dl Python kernel.\n"
        f"Current sys.prefix: {current_prefix}\n"
        f"Required environment: {expected_prefix}\n"
        "Select/register the kernel using the command in the preceding markdown cell."
    )
if not cfg.SCEPTRE_ENV.is_dir():
    raise FileNotFoundError(f"SCEPTRE environment not found: {cfg.SCEPTRE_ENV}")
if not cfg.SCEPTRE_RSCRIPT.is_file():
    raise FileNotFoundError(f"Direct SCEPTRE Rscript not found: {cfg.SCEPTRE_RSCRIPT}")
if not cfg.SC_DL_ENV.is_dir():
    raise FileNotFoundError(f"sc-dl environment not found: {cfg.SC_DL_ENV}")

available_fonts = {font.name for font in font_manager.fontManager.ttflist}
AI_FONT = next(
    (name for name in ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]
     if name in available_fonts),
    "DejaVu Sans",
)
mpl.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": [AI_FONT],
})

print("Output directory:", OUTDIR)
print("Notebook Python:", sys.executable)
print("Notebook environment:", current_prefix)
print("SCEPTRE environment:", cfg.SCEPTRE_ENV)
print("Direct Rscript:", cfg.SCEPTRE_RSCRIPT)
print("Private R library:", cfg.SCEPTRE_PRIVATE_R_LIB)
print("Required SCEPTRE version:", cfg.SCEPTRE_REQUIRED_VERSION)
print("Font:", AI_FONT)
print("plotnine:", p9.__version__)
print("Python:", sys.version.split()[0])
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/muon/_core/preproc.py:31: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
  if Version(scanpy.__version__) < Version("1.10"):

```

### Text output (stream)

```text
Output directory: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/sceptre_on_target_KD_D0_D3_overlay_plotnine_sceptre_env_v2
Notebook Python: /oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python
Notebook environment: /oak/stanford/groups/engreitz/Users/tri/envs/sc-dl
SCEPTRE environment: /oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env
Direct Rscript: /oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/bin/Rscript
Private R library: /oak/stanford/groups/engreitz/Users/tri/Rlibs/sceptre_0.10.3
Required SCEPTRE version: 0.10.3
Font: DejaVu Sans
plotnine: 0.15.4
Python: 3.11.14

```

## Markdown cell 4

## 2. Validate or recreate the four SCEPTRE transfer folders

The original notebook transferred each day's MuData into Matrix Market files plus names, covariates and mappings. This cell reuses those files when complete. If a folder is missing or incomplete, it recreates it from the original normalized H5MU file while preserving the original definitions:

- response matrix: `GEX.layers['raw_counts']`, transposed to genes × cells;
- gRNA matrix: `guide.X`, transposed to guides × cells;
- batch: `GEX.obs['sample_id']`;
- response/gRNA depth and complexity fields copied exactly as in the original transfer table;
- mitochondrial percentage: `GEX.obs['pct_counts_mt']`;
- gRNA target: `intended_target_name` from `300genes_guide_metadata_v43.tsv`;
- response mapping: `GEX.var['symbol']` paired with the Ensembl response ID.

No cells are silently removed during transfer, matching the original notebook.


## Code cell 5

```python
TRANSFER_FILENAMES = [
    "response_matrix.mtx",
    "grna_matrix.mtx",
    "gene_names.csv",
    "guide_names.csv",
    "cell_names.csv",
    "covariates.csv",
    "guide_mapping.csv",
    "gene_symbol_map.csv",
]


def transfer_dir(day: str) -> Path:
    return cfg.BASEDIR / f"{day}_sceptre_transfer"


def transfer_paths(day: str) -> dict[str, Path]:
    directory = transfer_dir(day)
    return {name: directory / name for name in TRANSFER_FILENAMES}


def complete_transfer(day: str) -> bool:
    return all(path.is_file() and path.stat().st_size > 0 for path in transfer_paths(day).values())


def prepare_transfer(day: str, guide_metadata: pd.DataFrame | None = None) -> dict[str, Path]:
    paths = transfer_paths(day)
    directory = transfer_dir(day)
    directory.mkdir(parents=True, exist_ok=True)

    if complete_transfer(day) and not cfg.FORCE_REBUILD_TRANSFER:
        print(f"{day.upper()}: reusing complete transfer folder: {directory}")
        return paths

    if guide_metadata is None:
        if not cfg.GUIDE_METADATA_TSV.is_file():
            raise FileNotFoundError(cfg.GUIDE_METADATA_TSV)
        guide_metadata = pd.read_csv(cfg.GUIDE_METADATA_TSV, sep="\t")

    required_meta = {"guide_id", "intended_target_name"}
    missing_meta = required_meta - set(guide_metadata.columns)
    if missing_meta:
        raise ValueError(f"Guide metadata lacks columns: {sorted(missing_meta)}")
    if guide_metadata["guide_id"].duplicated().any():
        duplicated = guide_metadata.loc[guide_metadata["guide_id"].duplicated(), "guide_id"].head().tolist()
        raise ValueError(f"Duplicated guide_id values in metadata, e.g. {duplicated}")

    h5mu_path = cfg.NORMALIZED_DATA_DIR / f"{day}_normalized_gex_and_guide.h5mu"
    if not h5mu_path.is_file():
        raise FileNotFoundError(h5mu_path)

    print(f"{day.upper()}: rebuilding transfer files from {h5mu_path}")
    mdata = mu.read(h5mu_path)
    if "GEX" not in mdata.mod or "guide" not in mdata.mod:
        raise KeyError(f"{h5mu_path} must contain GEX and guide modalities.")
    gex = mdata["GEX"]
    guide = mdata["guide"]

    if not gex.obs_names.equals(guide.obs_names):
        raise ValueError(f"GEX and guide cell orders differ for {day}; do not guess an alignment.")
    if "raw_counts" not in gex.layers:
        raise KeyError(f"GEX.layers['raw_counts'] is missing for {day}.")
    required_obs = {
        "sample_id", "guides_per_cell", "guide_umi_counts",
        "n_genes_by_counts", "total_counts", "pct_counts_mt",
    }
    missing_obs = required_obs - set(gex.obs.columns)
    if missing_obs:
        raise KeyError(f"Missing GEX.obs fields for {day}: {sorted(missing_obs)}")
    if "symbol" not in gex.var.columns:
        raise KeyError(f"GEX.var['symbol'] is missing for {day}.")

    sio.mmwrite(paths["response_matrix.mtx"], sparse.csr_matrix(gex.layers["raw_counts"]).T)
    sio.mmwrite(paths["grna_matrix.mtx"], sparse.csr_matrix(guide.X).T)
    pd.Series(gex.var_names.astype(str)).to_csv(paths["gene_names.csv"], index=False, header=False)
    pd.Series(guide.var_names.astype(str)).to_csv(paths["guide_names.csv"], index=False, header=False)
    pd.Series(gex.obs_names.astype(str)).to_csv(paths["cell_names.csv"], index=False, header=False)

    covariates = pd.DataFrame({
        "batch": gex.obs["sample_id"].to_numpy(),
        "grna_n_nonzero": gex.obs["guides_per_cell"].to_numpy(),
        "grna_n_umis": gex.obs["guide_umi_counts"].to_numpy(),
        "response_n_nonzero": gex.obs["n_genes_by_counts"].to_numpy(),
        "response_n_umis": gex.obs["total_counts"].to_numpy(),
        "response_p_mito": gex.obs["pct_counts_mt"].to_numpy(),
    }, index=gex.obs_names.astype(str))
    covariates.to_csv(paths["covariates.csv"])

    lookup = guide_metadata.set_index("guide_id")
    guide_ids = pd.Index(guide.var_names.astype(str), name="grna_id")
    missing_guides = guide_ids.difference(lookup.index.astype(str))
    if len(missing_guides):
        raise ValueError(f"{len(missing_guides)} guides are absent from metadata; examples: {missing_guides[:5].tolist()}")
    # Reindex after converting the lookup index to string to avoid mixed-type surprises.
    lookup = lookup.copy()
    lookup.index = lookup.index.astype(str)
    mapped_targets = lookup.loc[guide_ids, "intended_target_name"]
    if mapped_targets.isna().any():
        bad = guide_ids[mapped_targets.isna().to_numpy()][:5].tolist()
        raise ValueError(f"Unmapped guide targets remain for {day}; examples: {bad}")
    guide_mapping = pd.DataFrame({
        "grna_id": guide_ids,
        "grna_target": mapped_targets.astype(str).to_numpy(),
    })
    guide_mapping.to_csv(paths["guide_mapping.csv"], index=False)

    if gex.var["symbol"].isna().any():
        raise ValueError(f"Missing gene symbols remain in {day}; do not create ambiguous on-target pairs.")
    gene_map = pd.DataFrame({
        "response_id": gex.var_names.astype(str),
        "grna_target": gex.var["symbol"].astype(str).to_numpy(),
    })
    gene_map.to_csv(paths["gene_symbol_map.csv"], index=False)

    missing_targets = sorted(set(cfg.TARGETS) - set(gene_map["grna_target"]))
    if missing_targets:
        raise ValueError(f"Target response genes absent from {day}: {missing_targets}")

    del mdata, gex, guide
    gc.collect()
    if not complete_transfer(day):
        missing = [str(path) for path in paths.values() if not path.is_file() or path.stat().st_size == 0]
        raise IOError(f"Transfer creation failed for {day}: {missing}")
    print(f"{day.upper()}: transfer files written and verified.")
    return paths


guide_metadata = None
if cfg.FORCE_REBUILD_TRANSFER or not all(complete_transfer(day) for day in cfg.DAYS):
    if not cfg.GUIDE_METADATA_TSV.is_file():
        raise FileNotFoundError(cfg.GUIDE_METADATA_TSV)
    guide_metadata = pd.read_csv(cfg.GUIDE_METADATA_TSV, sep="\t")

TRANSFER_PATHS = {day: prepare_transfer(day, guide_metadata) for day in cfg.DAYS}
print("All transfer folders are complete.")
```

### Text output (stream)

```text
D0: reusing complete transfer folder: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d0_sceptre_transfer
D1: reusing complete transfer folder: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d1_sceptre_transfer
D2: reusing complete transfer folder: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer
D3: reusing complete transfer folder: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d3_sceptre_transfer
All transfer folders are complete.

```

## Markdown cell 6

## 3. Install/verify SCEPTRE v0.10.3, then run or reuse each day's discovery analysis

The direct R executable is:

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/bin/Rscript
```

The next cell first runs a small installer/preflight script. It checks `.libPaths()`, verifies the package version and package location, and—when needed—installs the exact official `v0.10.3` tag into:

```text
/oak/stanford/groups/engreitz/Users/tri/Rlibs/sceptre_0.10.3
```

The existing v0.10.2 installation in `sceptre_env` is not removed or overwritten. The analysis subprocess receives the private directory through `R_LIBS_USER` and verifies both `packageVersion("sceptre") == "0.10.3"` and the resolved package path before loading the analysis inputs.

The SCEPTRE analysis retains the original choices:

- high-MOI import;
- on-target response–target pairs only;
- two-sided testing;
- union integration of guides targeting the same gene;
- complement controls;
- conditional randomization testing with skew-normal approximation;
- the original/default high-MOI adjustment formula;
- mixture-model gRNA assignment;
- `p_mito_threshold = 0.25` in QC;
- calibration check before discovery analysis;
- BH correction within each day's complete on-target discovery family.

The official fold-change standard errors are available from v0.10.3. The notebook deliberately does not manufacture a standard error from the SCEPTRE p-value or from `z_orig`. Existing complete v0.10.3 result files are reused unless `FORCE_RERUN_SCEPTRE=True`.

## Code cell 7

```python
def run_streamed(command: list[str], *, env: dict[str, str], log_path: Path) -> None:
    """Run a subprocess, stream combined stdout/stderr, and save the complete log."""
    print("Running external command:")
    print(" ".join(command))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log_handle:
        process = subprocess.Popen(
            command,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="")
            log_handle.write(line)
            log_handle.flush()
        return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)


def direct_rscript_command(script_path: Path) -> list[str]:
    """Use the user's exact sceptre_env Rscript; never route through conda run."""
    script_path = Path(script_path).expanduser().resolve()
    if not cfg.SCEPTRE_RSCRIPT.is_file():
        raise FileNotFoundError(cfg.SCEPTRE_RSCRIPT)
    return [str(cfg.SCEPTRE_RSCRIPT), "--vanilla", str(script_path)]


# Keep numerical libraries single-threaded inside the R subprocess. SCEPTRE itself
# is also run with parallel=False, matching the original notebook.
R_SUBPROCESS_ENV = os.environ.copy()
R_SUBPROCESS_ENV.update({
    "R_LIBS_USER": str(cfg.SCEPTRE_PRIVATE_R_LIB),
    "SCEPTRE_PRIVATE_R_LIB": str(cfg.SCEPTRE_PRIVATE_R_LIB),
    "SCEPTRE_REQUIRED_VERSION": cfg.SCEPTRE_REQUIRED_VERSION,
    "SCEPTRE_GITHUB_REPO": cfg.SCEPTRE_GITHUB_REPO,
    "SCEPTRE_GITHUB_REF": cfg.SCEPTRE_GITHUB_REF,
    "SCEPTRE_INSTALL_IF_NEEDED": str(cfg.INSTALL_PRIVATE_SCEPTRE_IF_NEEDED).lower(),
    "SCEPTRE_CRAN_REPOSITORY": cfg.CRAN_REPOSITORY,
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
})

INSTALL_SCRIPT_PATH = OUTDIR / "ensure_sceptre_0.10.3_private.R"
INSTALL_LOG_PATH = OUTDIR / "ensure_sceptre_0.10.3_private.log"

INSTALL_SCRIPT_TEXT = r"""
parse_bool <- function(x) tolower(trimws(x)) %in% c("true", "t", "1", "yes", "y")
env_or_default <- function(name, default) {
  value <- Sys.getenv(name)
  if (nzchar(value)) value else default
}
private_lib <- env_or_default(
  "SCEPTRE_PRIVATE_R_LIB",
  "/oak/stanford/groups/engreitz/Users/tri/Rlibs/sceptre_0.10.3"
)
required_version_text <- env_or_default("SCEPTRE_REQUIRED_VERSION", "0.10.3")
required_version <- package_version(required_version_text)
repo <- env_or_default("SCEPTRE_GITHUB_REPO", "Katsevich-Lab/sceptre")
ref <- env_or_default("SCEPTRE_GITHUB_REF", "v0.10.3")
install_if_needed <- parse_bool(env_or_default("SCEPTRE_INSTALL_IF_NEEDED", "true"))
cran_repository <- env_or_default("SCEPTRE_CRAN_REPOSITORY", "https://cloud.r-project.org")
output_dir <- env_or_default("SCEPTRE_OUTPUT_DIR", private_lib)

if (!nzchar(private_lib)) stop("SCEPTRE_PRIVATE_R_LIB is empty.")
dir.create(private_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(private_lib, .libPaths())))

package_location <- function(package_name) {
  tryCatch(find.package(package_name, quiet = TRUE), error = function(e) "")
}
package_version_from_path <- function(path) {
  if (!nzchar(path)) return(package_version("0.0.0"))
  description_path <- file.path(path, "DESCRIPTION")
  if (!file.exists(description_path)) return(package_version("0.0.0"))
  package_version(read.dcf(description_path, fields = "Version")[1, 1])
}

current_path <- package_location("sceptre")
current_version <- package_version_from_path(current_path)
cat("R executable:", file.path(R.home("bin"), "R"), "\n")
cat("Initial .libPaths():", paste(.libPaths(), collapse = " | "), "\n")
cat("Initially resolved SCEPTRE path:", ifelse(nzchar(current_path), current_path, "<not found>"), "\n")
cat("Initially resolved SCEPTRE version:", as.character(current_version), "\n")

need_install <- current_version != required_version
if (need_install) {
  if (!install_if_needed) {
    stop(
      paste0(
        "Exact SCEPTRE ", required_version_text, " is not available. ",
        "Set INSTALL_PRIVATE_SCEPTRE_IF_NEEDED=True or install ", repo, "@", ref,
        " into ", private_lib, "."
      )
    )
  }
  cat("Installing exact SCEPTRE tag", ref, "into", private_lib, "\n")
  remotes_path <- package_location("remotes")
  if (!nzchar(remotes_path)) {
    install.packages(
      "remotes", repos = cran_repository, lib = private_lib,
      dependencies = c("Depends", "Imports", "LinkingTo")
    )
  }
  remotes::install_github(
    repo = repo,
    ref = ref,
    lib = private_lib,
    dependencies = FALSE,
    upgrade = "never",
    force = TRUE,
    build_vignettes = FALSE,
    quiet = FALSE
  )
}

# Re-resolve after installation. This script never attaches the old namespace.
.libPaths(unique(c(private_lib, .libPaths())))
final_path <- package_location("sceptre")
final_version <- package_version_from_path(final_path)
if (!nzchar(final_path)) stop("SCEPTRE is still unavailable after preflight/install.")
if (final_version != required_version) {
  stop(
    paste0(
      "Expected exact SCEPTRE ", required_version_text,
      " but resolved ", final_version, " at ", final_path, "."
    )
  )
}
private_norm <- normalizePath(private_lib, winslash = "/", mustWork = TRUE)
final_norm <- normalizePath(final_path, winslash = "/", mustWork = TRUE)
if (!startsWith(final_norm, paste0(private_norm, "/"))) {
  warning(
    paste0(
      "Exact SCEPTRE version was found outside the private library at ", final_norm,
      ". The version is acceptable, but the package location differs from the preferred private library."
    )
  )
}
cat("Verified SCEPTRE version:", as.character(final_version), "\n")
cat("Verified SCEPTRE package path:", final_path, "\n")
writeLines(
  c(
    paste0("sceptre_version=", as.character(final_version)),
    paste0("sceptre_path=", final_path),
    paste0("private_library=", private_lib),
    paste0("r_version=", R.version.string)
  ),
  file.path(output_dir, "sceptre_private_library_preflight.txt")
)
"""

R_SUBPROCESS_ENV["SCEPTRE_OUTPUT_DIR"] = str(OUTDIR)
INSTALL_SCRIPT_PATH.write_text(INSTALL_SCRIPT_TEXT, encoding="utf-8")
INSTALL_COMMAND_USED = direct_rscript_command(INSTALL_SCRIPT_PATH)
run_streamed(INSTALL_COMMAND_USED, env=R_SUBPROCESS_ENV, log_path=INSTALL_LOG_PATH)

preflight_record = OUTDIR / "sceptre_private_library_preflight.txt"
if not preflight_record.is_file() or preflight_record.stat().st_size == 0:
    raise IOError("SCEPTRE private-library preflight did not produce its verification record.")
print(preflight_record.read_text(encoding="utf-8"))

# Variables consumed by the main analysis R script.
R_SUBPROCESS_ENV.update({
    "SCEPTRE_BASEDIR": str(cfg.BASEDIR),
    "SCEPTRE_DAYS": ",".join(cfg.DAYS),
    "SCEPTRE_RESULT_FILENAME": cfg.RESULT_FILENAME,
    "SCEPTRE_FORCE_RERUN": str(cfg.FORCE_RERUN_SCEPTRE).lower(),
    "SCEPTRE_RUN_CALIBRATION": str(cfg.RUN_CALIBRATION_CHECK).lower(),
    "SCEPTRE_SAVE_RDS": str(cfg.SAVE_SCEPTRE_RDS).lower(),
    "SCEPTRE_PARALLEL": str(cfg.R_PARALLEL).lower(),
    "SCEPTRE_Q_THRESHOLD": repr(cfg.Q_THRESHOLD),
    "SCEPTRE_RANDOM_SEED": str(cfg.RANDOM_SEED),
})

R_SCRIPT_PATH = OUTDIR / "run_sceptre_with_official_se.R"
R_LOG_PATH = OUTDIR / "run_sceptre_with_official_se.log"
```

### Text output (stream)

```text
Running external command:
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/bin/Rscript --vanilla /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/sceptre_on_target_KD_D0_D3_overlay_plotnine_sceptre_env_v2/ensure_sceptre_0.10.3_private.R
R executable: /oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/R/bin/R 
Initial .libPaths(): /oak/stanford/groups/engreitz/Users/tri/Rlibs/sceptre_0.10.3 | /oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/R/library 
Initially resolved SCEPTRE path: /oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/R/library/sceptre 
Initially resolved SCEPTRE version: 0.10.3 
Warning message:
Exact SCEPTRE version was found outside the private library at /oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/R/library/sceptre. The version is acceptable, but the package location differs from the preferred private library. 
Verified SCEPTRE version: 0.10.3 
Verified SCEPTRE package path: /oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/R/library/sceptre 
sceptre_version=0.10.3
sceptre_path=/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/R/library/sceptre
private_library=/oak/stanford/groups/engreitz/Users/tri/Rlibs/sceptre_0.10.3
r_version=R version 4.5.1 (2025-06-13)


```

## Code cell 8

```python
R_SCRIPT_TEXT = r'''
private_lib <- Sys.getenv("SCEPTRE_PRIVATE_R_LIB")
if (!nzchar(private_lib)) stop("SCEPTRE_PRIVATE_R_LIB is empty.")
.libPaths(unique(c(private_lib, .libPaths())))

suppressPackageStartupMessages({
  library(sceptre)
  library(Matrix)
})

parse_bool <- function(x) tolower(trimws(x)) %in% c("true", "t", "1", "yes", "y")
base_dir <- Sys.getenv("SCEPTRE_BASEDIR")
days <- strsplit(Sys.getenv("SCEPTRE_DAYS"), ",", fixed = TRUE)[[1]]
result_filename <- Sys.getenv("SCEPTRE_RESULT_FILENAME")
force_rerun <- parse_bool(Sys.getenv("SCEPTRE_FORCE_RERUN"))
run_calibration <- parse_bool(Sys.getenv("SCEPTRE_RUN_CALIBRATION"))
save_rds <- parse_bool(Sys.getenv("SCEPTRE_SAVE_RDS"))
run_parallel <- parse_bool(Sys.getenv("SCEPTRE_PARALLEL"))
q_threshold <- as.numeric(Sys.getenv("SCEPTRE_Q_THRESHOLD"))
base_seed <- as.integer(Sys.getenv("SCEPTRE_RANDOM_SEED"))
required_version <- package_version(Sys.getenv("SCEPTRE_REQUIRED_VERSION"))
installed_version <- packageVersion("sceptre")
installed_path <- find.package("sceptre")

cat("SCEPTRE version:", as.character(installed_version), "\n")
cat("SCEPTRE package path:", installed_path, "\n")
cat("R library paths:", paste(.libPaths(), collapse = " | "), "\n")
if (installed_version != required_version) {
  stop(paste0("Expected exact SCEPTRE ", required_version, " but loaded ", installed_version, " from ", installed_path, "."))
}
if (!is.finite(q_threshold) || q_threshold <= 0 || q_threshold >= 1) {
  stop("SCEPTRE_Q_THRESHOLD must lie strictly between 0 and 1.")
}

required_result_columns <- c(
  "response_id", "grna_target", "n_nonzero_trt", "n_nonzero_cntrl",
  "pass_qc", "p_value", "fold_change", "se_fold_change",
  "log_2_fold_change", "significant"
)

to_logical_strict <- function(x, column_name) {
  if (is.logical(x)) return(x)
  values <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(values))
  out[values %in% c("true", "t", "1")] <- TRUE
  out[values %in% c("false", "f", "0")] <- FALSE
  unknown <- is.na(out) & !is.na(x)
  if (any(unknown)) {
    stop(paste("Cannot parse", column_name, "as logical; examples:", paste(head(unique(values[unknown])), collapse = ", ")))
  }
  out
}

normalise_result <- function(result, day_label) {
  result <- as.data.frame(result, stringsAsFactors = FALSE)
  missing <- setdiff(required_result_columns, names(result))
  if (length(missing) > 0) {
    stop(paste(day_label, "result lacks required columns:", paste(missing, collapse = ", ")))
  }
  numeric_columns <- c(
    "n_nonzero_trt", "n_nonzero_cntrl", "p_value",
    "fold_change", "se_fold_change", "log_2_fold_change"
  )
  for (column in numeric_columns) result[[column]] <- as.numeric(result[[column]])
  result$pass_qc <- to_logical_strict(result$pass_qc, "pass_qc")
  result$significant <- to_logical_strict(result$significant, "significant")
  result$q_bh <- NA_real_
  eligible <- which(result$pass_qc & is.finite(result$p_value))
  if (length(eligible) == 0) stop(paste(day_label, "has no QC-passing finite p-values."))
  result$q_bh[eligible] <- p.adjust(result$p_value[eligible], method = "BH")
  result$q_significant <- !is.na(result$q_bh) & result$q_bh <= q_threshold
  result$day <- day_label

  finite_fc <- !is.na(result$pass_qc) & result$pass_qc &
    is.finite(result$fold_change) & result$fold_change > 0 &
    is.finite(result$log_2_fold_change)
  if (any(finite_fc, na.rm = TRUE)) {
    discrepancy <- max(abs(log2(result$fold_change[finite_fc]) - result$log_2_fold_change[finite_fc]))
    if (!is.finite(discrepancy) || discrepancy > 1e-6) {
      stop(paste(day_label, "fold_change and log_2_fold_change are inconsistent; max discrepancy =", discrepancy))
    }
  }
  result
}

read_names <- function(path) {
  as.character(read.csv(path, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)[, 1])
}

for (day_index in seq_along(days)) {
  day <- days[[day_index]]
  day_dir <- file.path(base_dir, paste0(day, "_sceptre_transfer"))
  updated_path <- file.path(day_dir, result_filename)
  historical_path <- file.path(day_dir, "TargetKD.csv")
  result <- NULL
  result_source <- NULL

  if (!force_rerun) {
    for (candidate in c(updated_path, historical_path)) {
      if (!file.exists(candidate)) next
      candidate_result <- tryCatch(
        read.csv(candidate, check.names = FALSE, stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (!is.null(candidate_result) && all(required_result_columns %in% names(candidate_result))) {
        result <- normalise_result(candidate_result, day)
        result_source <- candidate
        break
      }
    }
  }

  if (is.null(result)) {
    cat("\n", day, ": running SCEPTRE from transfer inputs.\n", sep = "")
    required_inputs <- c(
      "response_matrix.mtx", "grna_matrix.mtx", "gene_names.csv",
      "guide_names.csv", "cell_names.csv", "covariates.csv",
      "guide_mapping.csv", "gene_symbol_map.csv"
    )
    missing_inputs <- required_inputs[!file.exists(file.path(day_dir, required_inputs))]
    if (length(missing_inputs) > 0) {
      stop(paste(day, "is missing transfer files:", paste(missing_inputs, collapse = ", ")))
    }

    response_mat <- as(readMM(file.path(day_dir, "response_matrix.mtx")), "dgCMatrix")
    grna_mat <- as(readMM(file.path(day_dir, "grna_matrix.mtx")), "dgCMatrix")
    gene_names <- read_names(file.path(day_dir, "gene_names.csv"))
    guide_names <- read_names(file.path(day_dir, "guide_names.csv"))
    cell_names <- read_names(file.path(day_dir, "cell_names.csv"))
    if (length(gene_names) != nrow(response_mat)) stop(paste(day, "gene-name count mismatch."))
    if (length(guide_names) != nrow(grna_mat)) stop(paste(day, "guide-name count mismatch."))
    if (length(cell_names) != ncol(response_mat) || length(cell_names) != ncol(grna_mat)) {
      stop(paste(day, "cell-name count mismatch."))
    }
    rownames(response_mat) <- gene_names
    colnames(response_mat) <- cell_names
    rownames(grna_mat) <- guide_names
    colnames(grna_mat) <- cell_names

    covariates <- read.csv(
      file.path(day_dir, "covariates.csv"), row.names = 1,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    if (!all(cell_names %in% rownames(covariates))) stop(paste(day, "covariates omit cells."))
    covariates <- covariates[cell_names, , drop = FALSE]
    if ("response_p_mito" %in% names(covariates)) {
      names(covariates)[names(covariates) == "response_p_mito"] <- "mt_percent"
    }
    reserved_columns <- c(
      "response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis"
    )
    covariates <- covariates[, !(names(covariates) %in% reserved_columns), drop = FALSE]

    guide_mapping <- read.csv(
      file.path(day_dir, "guide_mapping.csv"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    gene_map <- read.csv(
      file.path(day_dir, "gene_symbol_map.csv"),
      check.names = FALSE, stringsAsFactors = FALSE
    )

    set.seed(base_seed + day_index - 1L)
    sceptre_object <- import_data(
      response_matrix = response_mat,
      grna_matrix = grna_mat,
      grna_target_data_frame = guide_mapping,
      extra_covariates = covariates,
      moi = "high"
    )

    discovery_pairs <- merge(guide_mapping, gene_map, by = "grna_target")
    discovery_pairs <- unique(discovery_pairs[, c("grna_target", "response_id")])
    if (nrow(discovery_pairs) == 0) stop(paste(day, "has no on-target discovery pairs."))

    analysis_formula <- formula(
      ~ log(response_n_nonzero) + log(response_n_umis) +
        log(grna_n_nonzero) + log(grna_n_umis) + batch
    )
    sceptre_object <- set_analysis_parameters(
      sceptre_object = sceptre_object,
      discovery_pairs = discovery_pairs,
      side = "both",
      grna_integration_strategy = "union",
      formula_object = analysis_formula,
      resampling_approximation = "skew_normal",
      control_group = "complement",
      resampling_mechanism = "crt",
      multiple_testing_method = "BH",
      multiple_testing_alpha = q_threshold
    )
    sceptre_object <- assign_grnas(
      sceptre_object = sceptre_object,
      method = "mixture",
      parallel = run_parallel
    )
    sceptre_object <- run_qc(
      sceptre_object = sceptre_object,
      p_mito_threshold = 0.25
    )

    if (run_calibration) {
      sceptre_object <- run_calibration_check(
        sceptre_object = sceptre_object,
        parallel = run_parallel
      )
      calibration_result <- get_result(
        sceptre_object = sceptre_object,
        analysis = "run_calibration_check"
      )
      write.csv(
        as.data.frame(calibration_result),
        file.path(day_dir, "calibration_result_updated_sceptre.csv"),
        row.names = FALSE
      )
    }

    sceptre_object <- run_discovery_analysis(
      sceptre_object = sceptre_object,
      output_amount = 1,
      parallel = run_parallel
    )
    result <- get_result(
      sceptre_object = sceptre_object,
      analysis = "run_discovery_analysis"
    )
    result <- normalise_result(result, day)
    result_source <- "fresh_sceptre_rerun"

    if (save_rds) {
      saveRDS(sceptre_object, file.path(day_dir, "sceptre_object_v0.10.3_with_SE.rds"))
    }
    writeLines(
      capture.output(sessionInfo()),
      file.path(day_dir, "sceptre_sessionInfo_v0.10.3.txt")
    )
    writeLines(
      capture.output(print(sceptre_object)),
      file.path(day_dir, "sceptre_analysis_summary_v0.10.3.txt")
    )
    rm(
      response_mat, grna_mat, covariates, guide_mapping, gene_map,
      discovery_pairs, sceptre_object
    )
    gc()
  } else {
    cat(day, ": reusing result with standard errors from ", result_source, "\n", sep = "")
  }

  write.csv(result, updated_path, row.names = FALSE)
  writeLines(
    c(
      paste0("sceptre_version=", as.character(installed_version)),
      paste0("sceptre_path=", installed_path),
      paste0("result_source=", result_source),
      paste0("q_threshold=", q_threshold),
      paste0("n_pairs=", nrow(result)),
      paste0("n_qc_passing=", sum(result$pass_qc, na.rm = TRUE))
    ),
    file.path(day_dir, "TargetKD_with_SE_v0.10.3_metadata.txt")
  )
  cat(day, ": wrote ", updated_path, "\n", sep = "")
}
'''

R_SCRIPT_PATH.write_text(R_SCRIPT_TEXT, encoding="utf-8")
if not R_SCRIPT_PATH.is_file() or R_SCRIPT_PATH.stat().st_size == 0:
    raise IOError(f"Failed to write R script: {R_SCRIPT_PATH}")

R_COMMAND_USED = direct_rscript_command(R_SCRIPT_PATH)
run_streamed(R_COMMAND_USED, env=R_SUBPROCESS_ENV, log_path=R_LOG_PATH)

if not R_LOG_PATH.is_file() or R_LOG_PATH.stat().st_size == 0:
    raise IOError(f"SCEPTRE subprocess did not produce a log: {R_LOG_PATH}")
print("SCEPTRE subprocess completed successfully.")
print("R script:", R_SCRIPT_PATH)
print("R log:", R_LOG_PATH)

```

### Text output (stream)

```text
Running external command:
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/bin/Rscript --vanilla /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/sceptre_on_target_KD_D0_D3_overlay_plotnine_sceptre_env_v2/run_sceptre_with_official_se.R
SCEPTRE version: 0.10.3 
SCEPTRE package path: /oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/R/library/sceptre 
R library paths: /oak/stanford/groups/engreitz/Users/tri/Rlibs/sceptre_0.10.3 | /oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/R/library 

d0: running SCEPTRE from transfer inputs.
'as(<dgTMatrix>, "dgCMatrix")' is deprecated.
Use 'as(., "CsparseMatrix")' instead.
See help("Deprecated") and help("Matrix-deprecated").
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178488.23-P1P2 (1 of 1926)
Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178620.23-P1P2 (5 of 1926)
Performing gRNA-to-cell assignments for gRNA ACTL7B_-_111617995.23-P1P2 (10 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR1B_+_52345762.23-P1P2 (15 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2A_+_148602343.23-P1P2 (20 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_+_38495339.23-P1P2 (25 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_-_38495819.23-P1P2 (30 of 1926)
Performing gRNA-to-cell assignments for gRNA ADCY2_-_7396429.23-P1P2 (35 of 1926)
Performing gRNA-to-cell assignments for gRNA AMOTL2_+_134093279.23-P1P2 (40 of 1926)
Performing gRNA-to-cell assignments for gRNA ANP32E_+_150208342.23-P1P2 (45 of 1926)
Performing gRNA-to-cell assignments for gRNA APAF1_-_99039071.23-P1P2 (50 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_+_72766725.23-P1P2 (55 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_-_72766755.23-P1P2 (60 of 1926)
Performing gRNA-to-cell assignments for gRNA ARTN_-_44399429.23-P1 (65 of 1926)
Performing gRNA-to-cell assignments for gRNA ASCC3_-_101329193.23-P1P2 (70 of 1926)
Performing gRNA-to-cell assignments for gRNA ASH2L_-_37963059.23-P1P2 (75 of 1926)
Performing gRNA-to-cell assignments for gRNA ATMIN_+_81069661.23-P1P2 (80 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_+_111992112.23-P1P2 (85 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_-_111992129.23-P1P2 (90 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP6AP1_-_153657075.23-P1P2 (95 of 1926)
Performing gRNA-to-cell assignments for gRNA ATXN2L_-_28834668.23-P1P2 (100 of 1926)
Performing gRNA-to-cell assignments for gRNA BAD_+_64052102.23-P1P2 (105 of 1926)
Performing gRNA-to-cell assignments for gRNA BAG6_+_31620398.23-P1P2 (110 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_+_49458268.23-P1P2 (115 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_-_49458323.23-P1P2 (120 of 1926)
Performing gRNA-to-cell assignments for gRNA BCAS2_-_115124045.23-P1P2 (125 of 1926)
Performing gRNA-to-cell assignments for gRNA BMPR1A_-_88516386.23-P1P2 (130 of 1926)
Performing gRNA-to-cell assignments for gRNA BRD3_-_136932911.23-P1P2 (135 of 1926)
Performing gRNA-to-cell assignments for gRNA C17orf49_+_6918114.23-P2 (140 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_+_31532436.23-P1P2 (145 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_-_31532650.23-P1P2 (150 of 1926)
Performing gRNA-to-cell assignments for gRNA CACFD1_-_136325160.23-P1P2 (155 of 1926)
Performing gRNA-to-cell assignments for gRNA CACNG1_+_65040756.23-P1P2 (160 of 1926)
Performing gRNA-to-cell assignments for gRNA CAPN5_+_76778225.23-P1P2 (165 of 1926)
Performing gRNA-to-cell assignments for gRNA CASP3_+_185570560.23-P1P2 (170 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_+_107384793.23-P1P2 (175 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_-_107384635.23-P1P2 (180 of 1926)
Performing gRNA-to-cell assignments for gRNA CCM2_-_45039941.23-P1P2 (185 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNC_-_100016267.23-P1P2 (190 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNT2_-_135676416.23-P1P2 (195 of 1926)
Performing gRNA-to-cell assignments for gRNA CCT6B_+_33288381.23-P1P2 (200 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_+_30366413.23-P1P2 (205 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_-_30366630.23-P1P2 (210 of 1926)
Performing gRNA-to-cell assignments for gRNA CD99L2_-_150067120.23-P1P2 (215 of 1926)
Performing gRNA-to-cell assignments for gRNA CDH5_-_66400630.23-P1P2 (220 of 1926)
Performing gRNA-to-cell assignments for gRNA CDKN1A_-_36646538.23-P1P2 (225 of 1926)
Performing gRNA-to-cell assignments for gRNA CENPC_+_68411190.23-P1P2 (230 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_+_30773828.23-P1P2 (235 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_-_30773391.23-P1P2 (240 of 1926)
Performing gRNA-to-cell assignments for gRNA CHMP3_-_86790451.23-P1P2 (245 of 1926)
Performing gRNA-to-cell assignments for gRNA CHORDC1_-_89955894.23-P1P2 (250 of 1926)
Performing gRNA-to-cell assignments for gRNA CNEP1R1_+_50059284.23-P1P2 (255 of 1926)
Performing gRNA-to-cell assignments for gRNA CNOT2_+_70637651.23-P1P2 (260 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_+_83956339.23-P1P2 (265 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_-_83956572.23-P1P2 (270 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS6_-_99686676.23-P1P2 (275 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_-_31634088.23-P1P2 (280 of 1926)
Performing gRNA-to-cell assignments for gRNA CTR9_+_10772881.23-P1P2 (285 of 1926)
Performing gRNA-to-cell assignments for gRNA CUL3_-_225449667.23-P1P2 (290 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_+_145150267.23-P1P2 (295 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_-_145149996.23-P1P2 (300 of 1926)
Performing gRNA-to-cell assignments for gRNA CYCS_+_25164869.23-P1P2 (305 of 1926)
Performing gRNA-to-cell assignments for gRNA DACH1_-_72441101.23-P1P2 (310 of 1926)
Performing gRNA-to-cell assignments for gRNA DBR1_+_137893673.23-P1P2 (315 of 1926)
Performing gRNA-to-cell assignments for gRNA DCTN2_-_57940910.23-P1P2 (320 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_+_61851532.23-P1P2 (325 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_-_61851812.23-P1P2 (330 of 1926)
Performing gRNA-to-cell assignments for gRNA DGCR8_-_20068158.23-P1P2 (335 of 1926)
Performing gRNA-to-cell assignments for gRNA DNAJC17_-_41099215.23-P1P2 (340 of 1926)
Performing gRNA-to-cell assignments for gRNA DNM2_+_10829157.23-P1P2 (345 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM1_+_49575018.23-P1P2 (350 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_+_130700056.23-P1P2 (355 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_-_130700076.23-P1P2 (360 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM3_+_155113020.23-P1P2 (365 of 1926)
Performing gRNA-to-cell assignments for gRNA DROSHA_-_31531982.23-P1P2 (370 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_+_3985081.23-P1 (375 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_-_3985443.23-P1 (380 of 1926)
Performing gRNA-to-cell assignments for gRNA EHMT2_-_31865199.23-P1P2 (385 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3F_+_8009347.23-P1P2 (390 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3M_-_32605405.23-P1P2 (395 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4E_+_99849831.23-P1P2 (400 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184032364.23-P1P2 (405 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184033042.23-P1P2 (410 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G2_+_10830447.23-P1 (415 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF5A_-_7211217.23-P1P2 (420 of 1926)
Performing gRNA-to-cell assignments for gRNA TCEB3_-_24069975.23-P1P2 (425 of 1926)
Performing gRNA-to-cell assignments for gRNA EMC2_+_109456080.23-P1P2 (430 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823220.23-P1P2 (435 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823252.23-P1P2 (440 of 1926)
Performing gRNA-to-cell assignments for gRNA EOMES_-_27763342.23-P1P2 (445 of 1926)
Performing gRNA-to-cell assignments for gRNA ERCC2_+_45873802.23-P1P2 (450 of 1926)
Performing gRNA-to-cell assignments for gRNA ERO1L_+_53162425.23-P1P2 (455 of 1926)
Performing gRNA-to-cell assignments for gRNA ETS1_+_128392044.23-P1P2 (460 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_+_36134592.23-P1P2 (465 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_-_36134789.23-P1P2 (470 of 1926)
Performing gRNA-to-cell assignments for gRNA EWSR1_-_29664347.23-P1P2 (475 of 1926)
Performing gRNA-to-cell assignments for gRNA EXOSC8_-_37574979.23-P1P2 (480 of 1926)
Performing gRNA-to-cell assignments for gRNA FBXW11_-_171433300.23-P1P2 (485 of 1926)
Performing gRNA-to-cell assignments for gRNA FEN1_+_61560782.23-P1P2 (490 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_+_75532080.23-P1P2 (495 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_-_75532461.23-P1P2 (500 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT7_-_139926495.23-P1P2 (505 of 1926)
Performing gRNA-to-cell assignments for gRNA G3BP1_-_151151499.23-P1P2 (510 of 1926)
Performing gRNA-to-cell assignments for gRNA GCK_+_44229021.23-P1P2 (515 of 1926)
Performing gRNA-to-cell assignments for gRNA GOLT1B_+_21654777.23-P1P2 (520 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_+_45000529.23-P1P2 (525 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_-_45000948.23-P1P2 (530 of 1926)
Performing gRNA-to-cell assignments for gRNA GPI_+_34856284.23-P1P2 (535 of 1926)
Performing gRNA-to-cell assignments for gRNA GPR107_-_132816044.23-P1P2 (540 of 1926)
Performing gRNA-to-cell assignments for gRNA GSPT1_-_12009714.23-P1P2 (545 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A1_+_81687215.23-P1P2 (550 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_+_59949355.23-P1P2 (555 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_-_59949657.23-P1P2 (560 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2H4_-_30876185.23-P1P2 (565 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF3C2_+_27579604.23-P1P2 (570 of 1926)
Performing gRNA-to-cell assignments for gRNA HAND1_+_153857732.23-P1P2 (575 of 1926)
Performing gRNA-to-cell assignments for gRNA HDAC2_-_114292253.23-P1P2 (580 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_+_81399312.23-P1P2 (585 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_-_81399351.23-P1P2 (590 of 1926)
Performing gRNA-to-cell assignments for gRNA HINT2_+_35814995.23-P1P2 (595 of 1926)
Performing gRNA-to-cell assignments for gRNA HNF4A_+_43029971.23-P1P2 (600 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPM_-_8509885.23-P1 (605 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPU_+_245027753.23-P1P2 (610 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_+_43702292.23-P1P2 (615 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_-_43702305.23-P1P2 (620 of 1926)
Performing gRNA-to-cell assignments for gRNA HSPA5_-_128003534.23-P1P2 (625 of 1926)
Performing gRNA-to-cell assignments for gRNA HUWE1_-_53710540.23-P1P2 (630 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_-_30007588.23-P1P2 (635 of 1926)
Performing gRNA-to-cell assignments for gRNA INTS10_+_19675000.23-P1P2 (640 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_+_65903046.23-P1P2 (645 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_-_65903523.23-P1P2 (650 of 1926)
Performing gRNA-to-cell assignments for gRNA C7orf26_+_6630358.23-P1P2 (655 of 1926)
Performing gRNA-to-cell assignments for gRNA IPO8_+_30848771.23-P1P2 (660 of 1926)
Performing gRNA-to-cell assignments for gRNA ITGB1BP1_+_9563543.23-P1P2 (665 of 1926)
Performing gRNA-to-cell assignments for gRNA SSFA2_+_182756743.23-P1P2 (670 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_+_109914946.23-P1P2 (675 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_-_109915112.23-P1P2 (680 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM1A_-_23346409.23-P1P2 (685 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM3B_-_137688431.23-P1P2 (690 of 1926)
Performing gRNA-to-cell assignments for gRNA KDR_-_55991659.23-P1P2 (695 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA0232_+_6784651.23-P1P2 (700 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153706791.23-P1P2 (705 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153707248.23-P1P2 (710 of 1926)
Performing gRNA-to-cell assignments for gRNA LAMTOR4_-_99746921.23-P1P2 (715 of 1926)
Performing gRNA-to-cell assignments for gRNA LARP4_-_50794864.23-P1P2 (720 of 1926)
Performing gRNA-to-cell assignments for gRNA LATS2_-_21635336.23-P1P2 (725 of 1926)
Performing gRNA-to-cell assignments for gRNA LEMD3_+_65563441.23-P1P2 (730 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_+_54960428.23-P1P2 (735 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_-_54960616.23-P1P2 (740 of 1926)
Performing gRNA-to-cell assignments for gRNA LHFPL2_-_77944519.23-P1P2 (745 of 1926)
Performing gRNA-to-cell assignments for gRNA LSM12_-_42143603.23-P1P2 (750 of 1926)
Performing gRNA-to-cell assignments for gRNA MANBAL_+_35918136.23-P1P2 (755 of 1926)
Performing gRNA-to-cell assignments for gRNA MAU2_+_19431732.23-P1P2 (760 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_+_47807678.23-P1P2 (765 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_-_47807877.23-P1P2 (770 of 1926)
Performing gRNA-to-cell assignments for gRNA MEAF6_-_37980319.23-P1P2 (775 of 1926)
Performing gRNA-to-cell assignments for gRNA MED1_-_37607140.23-P1P2 (780 of 1926)
Performing gRNA-to-cell assignments for gRNA MED12_+_70338541.23-P1P2 (785 of 1926)
Performing gRNA-to-cell assignments for gRNA MED15_+_20862189.23-P1 (790 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_+_57479335.23-P1P2 (795 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_-_57479437.23-P1P2 (800 of 1926)
Performing gRNA-to-cell assignments for gRNA MEMO1_-_32264882.23-P1 (805 of 1926)
Performing gRNA-to-cell assignments for gRNA MIXL1_+_226411425.23-P1P2 (810 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL28_+_420487.23-P1P2 (815 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL32_+_42972018.23-P1P2 (820 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_+_179322181.23-P1P2 (825 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_-_179322213.23-P1P2 (830 of 1926)
Performing gRNA-to-cell assignments for gRNA LUZP6_-_135612750.23-P1P2 (835 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA20_-_19998023.23-P1P2 (840 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA25_+_112546538.23-P1P2 (845 of 1926)
Performing gRNA-to-cell assignments for gRNA NCBP1_+_100396082.23-P1P2 (850 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_+_172428630.23-P1P2 (855 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_-_172428841.23-P1P2 (860 of 1926)
Performing gRNA-to-cell assignments for gRNA NEDD8_-_24701363.23-P1P2 (865 of 1926)
Performing gRNA-to-cell assignments for gRNA NF2_+_29999635.23-P1P2 (870 of 1926)
Performing gRNA-to-cell assignments for gRNA NFIB_+_14314559.23-P1P2 (875 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_+_139418329.23-P2 (880 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139418349.23-P2 (885 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139440503.23-P1 (890 of 1926)
Performing gRNA-to-cell assignments for gRNA NUFIP2_+_27621073.23-P1P2 (895 of 1926)
Performing gRNA-to-cell assignments for gRNA NUP93_-_56764090.23-P1P2 (900 of 1926)
Performing gRNA-to-cell assignments for gRNA OGT_+_70752948.23-P1P2 (905 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_+_59383170.23-P1P2 (910 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_-_59383348.23-P1P2 (915 of 1926)
Performing gRNA-to-cell assignments for gRNA PAM16_-_4401236.23-P1P2 (920 of 1926)
Performing gRNA-to-cell assignments for gRNA PAXBP1_-_34143614.23-P1P2 (925 of 1926)
Performing gRNA-to-cell assignments for gRNA PBLD_-_70092436.23-P1 (930 of 1926)
Performing gRNA-to-cell assignments for gRNA PCBP2_+_53835614.23-P2 (935 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_+_94292752.23-P2 (940 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_-_94293027.23-P2 (945 of 1926)
Performing gRNA-to-cell assignments for gRNA PGRMC2_+_129208913.23-P1P2 (950 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_-_178866211.23-P1P2 (955 of 1926)
Performing gRNA-to-cell assignments for gRNA PITPNC1_-_65373433.23-P1P2 (960 of 1926)
Performing gRNA-to-cell assignments for gRNA PLCG1_+_39766176.23-P1P2 (965 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_+_151042273.23-P2 (970 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_-_151042422.23-P2 (975 of 1926)
Performing gRNA-to-cell assignments for gRNA PLRG1_-_155471154.23-P1P2 (980 of 1926)
Performing gRNA-to-cell assignments for gRNA PLXNA2_-_208417275.23-P1P2 (985 of 1926)
Performing gRNA-to-cell assignments for gRNA POLE4_-_75185725.23-P1P2 (990 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR2M_+_57998962.23-P1P2 (995 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940045.23-P1P2 (1000 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940414.23-P1P2 (1005 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3K_-_103144.23-P1P2 (1010 of 1926)
Performing gRNA-to-cell assignments for gRNA PPHLN1_-_42720007.23-P1P2 (1015 of 1926)
Performing gRNA-to-cell assignments for gRNA PPIL4_+_149867146.23-P1P2 (1020 of 1926)
Performing gRNA-to-cell assignments for gRNA PPP2R1A_+_52693317.23-P1 (1025 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_+_53195217.23-P1P2 (1030 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_-_53195332.23-P1P2 (1035 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKD1_-_30396996.23-P1P2 (1040 of 1926)
Performing gRNA-to-cell assignments for gRNA PRMT5_+_23398554.23-P1 (1045 of 1926)
Performing gRNA-to-cell assignments for gRNA C2orf53_-_27362230.23-P1P2 (1050 of 1926)
Performing gRNA-to-cell assignments for gRNA PRSS8_+_31147010.23-P1P2 (1055 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_+_36909062.23-P1P2 (1060 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_-_36909230.23-P1P2 (1065 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMD4_-_151227289.23-P1P2 (1070 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMG4_-_3259227.23-P1P2 (1075 of 1926)
Performing gRNA-to-cell assignments for gRNA PTAR1_+_72374835.23-P1P2 (1080 of 1926)
Performing gRNA-to-cell assignments for gRNA PTBP1_+_797463.23-P1P2 (1085 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_+_89623388.23-P1P2 (1090 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_-_89623526.23-P1P2 (1095 of 1926)
Performing gRNA-to-cell assignments for gRNA PTMS_-_6875591.23-P1P2 (1100 of 1926)
Performing gRNA-to-cell assignments for gRNA PYGO2_-_154934135.23-P1P2 (1105 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB18_-_27793292.23-P1P2 (1110 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB6A_+_73471983.23-P1P2 (1115 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_+_76251918.23-P1P2 (1120 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_-_76252008.23-P1P2 (1125 of 1926)
Performing gRNA-to-cell assignments for gRNA GNB2L1_-_180670840.23-P1P2 (1130 of 1926)
Performing gRNA-to-cell assignments for gRNA RAD21_-_117886635.23-P1P2 (1135 of 1926)
Performing gRNA-to-cell assignments for gRNA RANGRF_+_8192257.23-P1P2 (1140 of 1926)
Performing gRNA-to-cell assignments for gRNA RBBP5_+_205091050.23-P1P2 (1145 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM14-RBM4_+_66406280.23-P2 (1150 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM33_+_155437802.23-P1P2 (1155 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM4_+_66406319.23-P2 (1160 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_+_114271352.23-P1P2 (1165 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_-_114271700.23-P1P2 (1170 of 1926)
Performing gRNA-to-cell assignments for gRNA RBX1_-_41347852.23-P1P2 (1175 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1551_+_32112399.23-P1P2 (1180 of 1926)
Performing gRNA-to-cell assignments for gRNA RFX3_-_3525750.23-P1P2 (1185 of 1926)
Performing gRNA-to-cell assignments for gRNA RGP1_+_35749328.23-P1P2 (1190 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_+_95067102.23-P1P2 (1195 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_-_95067164.23-P1P2 (1200 of 1926)
Performing gRNA-to-cell assignments for gRNA RIBC1_-_53449606.23-P1P2 (1205 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1432_+_5629173.23-P1P2 (1210 of 1926)
Performing gRNA-to-cell assignments for gRNA RLF_+_40627522.23-P1P2 (1215 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF115_-_145611148.23-P1P2 (1220 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_+_104296142.23-P1P2 (1225 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_-_104296250.23-P1P2 (1230 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF40_-_30773681.23-P1P2 (1235 of 1926)
Performing gRNA-to-cell assignments for gRNA RPL37A_-_217363600.23-P1P2 (1240 of 1926)
Performing gRNA-to-cell assignments for gRNA RPS15A_-_18801590.23-P1P2 (1245 of 1926)
Performing gRNA-to-cell assignments for gRNA SAV1_+_51134992.23-P1P2 (1250 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_+_33104038.23-P1P2 (1255 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_-_33104343.23-P1P2 (1260 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF8_+_155054976.23-P1P2 (1265 of 1926)
Performing gRNA-to-cell assignments for gRNA EPT1_+_26569050.23-P1P2 (1270 of 1926)
Performing gRNA-to-cell assignments for gRNA SERF2_+_44069659.23-P2 (1275 of 1926)
Performing gRNA-to-cell assignments for gRNA SERPING1_+_57365137.23-P1P2 (1280 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_+_122242041.23-P1P2 (1285 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_-_122242294.23-P1P2 (1290 of 1926)
Performing gRNA-to-cell assignments for gRNA SETDB1_-_150898815.23-P1P2 (1295 of 1926)
Performing gRNA-to-cell assignments for gRNA SF1_-_64545937.23-P1P2 (1300 of 1926)
Performing gRNA-to-cell assignments for gRNA CCDC101_-_28565268.23-P1P2 (1305 of 1926)
Performing gRNA-to-cell assignments for gRNA SKA3_+_21750519.23-P1P2 (1310 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_+_78591497.23-P1P2 (1315 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_-_78591737.23-P1P2 (1320 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC29A1_-_44187438.23-P2 (1325 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC35C1_-_45826711.23-P1P2 (1330 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD1_+_146403428.23-P1P2 (1335 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD3_+_67358229.23-P1P2 (1340 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_+_48556623.23-P1P2 (1345 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_-_48556657.23-P1P2 (1350 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD5_-_135468694.23-P1P2 (1355 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA4_-_11071741.23-P1P2 (1360 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA5_+_144435243.23-P1P2 (1365 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCB1_+_24129228.23-P1 (1370 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_+_47822975.23-P1P2 (1375 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_-_47822893.23-P1P2 (1380 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCD3_-_150974212.23-P1 (1385 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCE1_-_38803669.23-P1P2 (1390 of 1926)
Performing gRNA-to-cell assignments for gRNA SMC2_-_106856631.23-P1P2 (1395 of 1926)
Performing gRNA-to-cell assignments for gRNA SMG5_+_156252537.23-P1P2 (1400 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_+_103852.23-P1P2 (1405 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_-_104242.23-P1P2 (1410 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP27_-_70121181.23-P1P2 (1415 of 1926)
Performing gRNA-to-cell assignments for gRNA SNUPN_-_75917960.23-P1P2 (1420 of 1926)
Performing gRNA-to-cell assignments for gRNA SOX4_+_21594213.23-P1 (1425 of 1926)
Performing gRNA-to-cell assignments for gRNA SPIN1_+_91003737.23-P1P2 (1430 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_+_2328680.23-P1P2 (1435 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_-_2328912.23-P1P2 (1440 of 1926)
Performing gRNA-to-cell assignments for gRNA SPTY2D1_-_18655846.23-P1P2 (1445 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP19_+_112197188.23-P1P2 (1450 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP54_-_35452244.23-P1P2 (1455 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP68_+_74068540.23-P1P2 (1460 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_+_126138730.23-P1P2 (1465 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_-_126138742.23-P1P2 (1470 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPRB_-_133524719.23-P1P2 (1475 of 1926)
Performing gRNA-to-cell assignments for gRNA SRRM2_-_2802714.23-P1P2 (1480 of 1926)
Performing gRNA-to-cell assignments for gRNA SRSF6_+_42086993.23-P1 (1485 of 1926)
Performing gRNA-to-cell assignments for gRNA SSB_+_170655450.23-P1P2 (1490 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_+_57102902.23-P1P2 (1495 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_-_57103355.23-P1P2 (1500 of 1926)
Performing gRNA-to-cell assignments for gRNA STAG2_-_123095658.23-P1 (1505 of 1926)
Performing gRNA-to-cell assignments for gRNA STK40_+_36851507.23-P1P2 (1510 of 1926)
Performing gRNA-to-cell assignments for gRNA SUPT20H_+_37633798.23-P1P2 (1515 of 1926)
Performing gRNA-to-cell assignments for gRNA SUZ12_+_30264157.23-P1P2 (1520 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_+_64901973.23-P1P2 (1525 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_-_64901981.23-P1P2 (1530 of 1926)
Performing gRNA-to-cell assignments for gRNA TADA2B_-_7045121.23-P1P2 (1535 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF10_+_6633406.23-P1P2 (1540 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF2_+_120845032.23-P1P2 (1545 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF3_-_7860540.23-P1P2 (1550 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_+_60640814.23-P1P2 (1555 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_-_60640498.23-P1P2 (1560 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF6L_-_62538954.23-P1P2 (1565 of 1926)
Performing gRNA-to-cell assignments for gRNA TAOK1_-_27717907.23-P1P2 (1570 of 1926)
Performing gRNA-to-cell assignments for gRNA FAM208A_+_56717273.23-P1P2 (1575 of 1926)
Performing gRNA-to-cell assignments for gRNA TBC1D9_+_141677481.23-P1P2 (1580 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_+_170863524.23-P1P2 (1585 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_-_170863481.23-P1P2 (1590 of 1926)
Performing gRNA-to-cell assignments for gRNA TBX3_-_115121998.23-P1P2 (1595 of 1926)
Performing gRNA-to-cell assignments for gRNA T_-_166581333.23-P1P2 (1600 of 1926)
Performing gRNA-to-cell assignments for gRNA TGFBR2_+_30648526.23-P1P2 (1605 of 1926)
Performing gRNA-to-cell assignments for gRNA THBD_+_23029196.23-P2 (1610 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_+_2977592.23-P1P2 (1615 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_-_2977643.23-P1P2 (1620 of 1926)
Performing gRNA-to-cell assignments for gRNA MESDC1_-_81293340.23-P1P2 (1625 of 1926)
Performing gRNA-to-cell assignments for gRNA TMEM41B_+_9336156.23-P1P2 (1630 of 1926)
Performing gRNA-to-cell assignments for gRNA TMTC3_-_88536127.23-P1P2 (1635 of 1926)
Performing gRNA-to-cell assignments for gRNA TOP1_+_39657538.23-P1P2 (1640 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_+_7590531.23-P1P2 (1645 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_-_7590704.23-P1P2 (1650 of 1926)
Performing gRNA-to-cell assignments for gRNA TRIM24_-_138145308.23-P1P2 (1655 of 1926)
Performing gRNA-to-cell assignments for gRNA TRRAP_-_98476465.23-P1P2 (1660 of 1926)
Performing gRNA-to-cell assignments for gRNA TSTD2_+_100395594.23-P1P2 (1665 of 1926)
Performing gRNA-to-cell assignments for gRNA UBE2E1_+_23847541.23-P1P2 (1670 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_+_38924011.23-P1P2 (1675 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_-_38924421.23-P1P2 (1680 of 1926)
Performing gRNA-to-cell assignments for gRNA UPF2_+_12084799.23-P1P2 (1685 of 1926)
Performing gRNA-to-cell assignments for gRNA UXT_-_47518505.23-P1P2 (1690 of 1926)
Performing gRNA-to-cell assignments for gRNA VHL_-_10183525.23-P1P2 (1695 of 1926)
Performing gRNA-to-cell assignments for gRNA VPS54_+_64246104.23-P1P2 (1700 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_+_28821698.23-P1P2 (1705 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_-_28822074.23-P1P2 (1710 of 1926)
Performing gRNA-to-cell assignments for gRNA WAPAL_+_88281533.23-P1P2 (1715 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR12_-_203776383.23-P1P2 (1720 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR77_+_111991860.23-P1P2 (1725 of 1926)
Performing gRNA-to-cell assignments for gRNA WIZ_-_15543016.23-P1P2 (1730 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_+_9983328.23-P1P2 (1735 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_-_9983510.23-P1P2 (1740 of 1926)
Performing gRNA-to-cell assignments for gRNA XPO5_-_43543670.23-P1P2 (1745 of 1926)
Performing gRNA-to-cell assignments for gRNA XRN1_-_142166648.23-P1P2 (1750 of 1926)
Performing gRNA-to-cell assignments for gRNA YARS2_+_32908787.23-P1P2 (1755 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB22_+_33285454.23-P1P2 (1760 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_+_6640093.23-P1P2 (1765 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_-_6640138.23-P1P2 (1770 of 1926)
Performing gRNA-to-cell assignments for gRNA ZC3H8_-_113012435.23-P1P2 (1775 of 1926)
Performing gRNA-to-cell assignments for gRNA ZNF281_+_200379151.23-P1P2 (1780 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00117 (1785 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00233 (1790 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00518 (1795 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00656 (1800 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01226 (1805 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01730 (1810 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02009 (1815 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02288 (1820 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02777 (1825 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03026 (1830 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03252 (1835 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03620 (1840 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_SGTN_04 (1845 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_SGTN_25 (1850 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_SGTN_26 (1855 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1044 (1860 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1414 (1865 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1554 (1870 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1785 (1875 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2009 (1880 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2271 (1885 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2425 (1890 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2614 (1895 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_293 (1900 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_541 (1905 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_725 (1910 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_860 (1915 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_SGTN_14 (1920 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_SGTN_01 (1925 of 1926)
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Constructing negative control pairs. ✓
Running precomputation on response ENSG00000160741.17 (1 of 259)
Running precomputation on response ENSG00000286221.1 (5 of 259)
Running precomputation on response ENSG00000162971.11 (10 of 259)
Running precomputation on response ENSG00000183087.15 (15 of 259)
Running precomputation on response ENSG00000228566.2 (20 of 259)
Running precomputation on response ENSG00000213235.3 (25 of 259)
Running precomputation on response ENSG00000120645.12 (30 of 259)
Running precomputation on response ENSG00000149100.13 (35 of 259)
Running precomputation on response ENSG00000142046.15 (40 of 259)
Running precomputation on response ENSG00000143748.18 (45 of 259)
Running precomputation on response ENSG00000149311.20 (50 of 259)
Running precomputation on response ENSG00000133980.5 (55 of 259)
Running precomputation on response ENSG00000196312.16 (60 of 259)
Running precomputation on response ENSG00000113522.14 (65 of 259)
Running precomputation on response ENSG00000259673.6 (70 of 259)
Running precomputation on response ENSG00000100580.8 (75 of 259)
Running precomputation on response ENSG00000235655.3 (80 of 259)
Running precomputation on response ENSG00000196323.14 (85 of 259)
Running precomputation on response ENSG00000204520.14 (90 of 259)
Running precomputation on response ENSG00000146834.15 (95 of 259)
Running precomputation on response ENSG00000102978.13 (100 of 259)
Running precomputation on response ENSG00000164758.7 (105 of 259)
Running precomputation on response ENSG00000171097.14 (110 of 259)
Running precomputation on response ENSG00000111704.11 (115 of 259)
Running precomputation on response ENSG00000008086.13 (120 of 259)
Running precomputation on response ENSG00000175390.15 (125 of 259)
Running precomputation on response ENSG00000100300.18 (130 of 259)
Running precomputation on response ENSG00000196678.14 (135 of 259)
Running precomputation on response ENSG00000119638.14 (140 of 259)
Running precomputation on response ENSG00000116120.11 (145 of 259)
Running precomputation on response ENSG00000112667.13 (150 of 259)
Running precomputation on response ENSG00000111203.13 (155 of 259)
Running precomputation on response ENSG00000160606.11 (160 of 259)
Running precomputation on response ENSG00000242612.7 (165 of 259)
Running precomputation on response ENSG00000103404.15 (170 of 259)
Running precomputation on response ENSG00000097046.13 (175 of 259)
Running precomputation on response ENSG00000213722.10 (180 of 259)
Running precomputation on response ENSG00000133104.14 (185 of 259)
Running precomputation on response ENSG00000107771.17 (190 of 259)
Running precomputation on response ENSG00000187720.15 (195 of 259)
Running precomputation on response ENSG00000108344.15 (200 of 259)
Running precomputation on response ENSG00000135424.19 (205 of 259)
Running precomputation on response ENSG00000198408.14 (210 of 259)
Running precomputation on response ENSG00000162374.18 (215 of 259)
Running precomputation on response ENSG00000187987.10 (220 of 259)
Running precomputation on response ENSG00000126464.15 (225 of 259)
Running precomputation on response ENSG00000125149.12 (230 of 259)
Running precomputation on response ENSG00000154079.6 (235 of 259)
Running precomputation on response ENSG00000162769.13 (240 of 259)
Running precomputation on response ENSG00000135763.10 (245 of 259)
Running precomputation on response ENSG00000290018.1 (250 of 259)
Running precomputation on response ENSG00000102241.12 (255 of 259)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00259&non-targeting_02159&non-targeting_02571&non-targeting_03029&non-targeting_03081 (1 of 95)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00355&non-targeting_02025&non-targeting_02089&non-targeting_03081&non-targeting_03384 (5 of 95)
Analyzing pairs containing gRNA group non-targeting_00518&non-targeting_03081&non-targeting_03342&non-targeting_03427&non-targeting_03635&non-targeting_03664 (10 of 95)
Analyzing pairs containing gRNA group non-targeting_01509&non-targeting_01812&non-targeting_02025&non-targeting_02677&non-targeting_02777&non-targeting_03384 (15 of 95)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_00231&non-targeting_01185&non-targeting_02025&non-targeting_03384&non-targeting_03427 (20 of 95)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_02193&non-targeting_02543&non-targeting_02837&non-targeting_03026&non-targeting_03342 (25 of 95)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_00117&non-targeting_00196&non-targeting_00555&non-targeting_02025&non-targeting_02193 (30 of 95)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00355&non-targeting_00518&non-targeting_01601&non-targeting_02677&non-targeting_03427 (35 of 95)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00390&non-targeting_01601&non-targeting_02009&non-targeting_02288&non-targeting_03029 (40 of 95)
Analyzing pairs containing gRNA group non-targeting_00231&non-targeting_00614&non-targeting_00629&non-targeting_00918&non-targeting_02193&non-targeting_02543 (45 of 95)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00390&non-targeting_01732&non-targeting_02193&non-targeting_03081&non-targeting_03384 (50 of 95)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_00196&non-targeting_00390&non-targeting_01730&non-targeting_03081&non-targeting_03252 (55 of 95)
Analyzing pairs containing gRNA group non-targeting_00043&non-targeting_00355&non-targeting_01601&non-targeting_02814&non-targeting_03081&non-targeting_03620 (60 of 95)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_01185&non-targeting_03029&non-targeting_03081&non-targeting_03180&non-targeting_03384 (65 of 95)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00259&non-targeting_00972&non-targeting_01601&non-targeting_02089&non-targeting_03620 (70 of 95)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_01858&non-targeting_02777&non-targeting_03158&non-targeting_03252&non-targeting_03483 (75 of 95)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00614&non-targeting_00625&non-targeting_01601&non-targeting_02837&non-targeting_03180 (80 of 95)
Analyzing pairs containing gRNA group non-targeting_00629&non-targeting_00656&non-targeting_02193&non-targeting_02543&non-targeting_03026&non-targeting_03483 (85 of 95)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_00213&non-targeting_00656&non-targeting_02911&non-targeting_03158&non-targeting_03483 (90 of 95)
Analyzing pairs containing gRNA group non-targeting_00629&non-targeting_00794&non-targeting_02962&non-targeting_03029&non-targeting_03620&non-targeting_03635 (95 of 95)
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Running precomputation on response ENSG00000004487.18 (1 of 261)
Running precomputation on response ENSG00000117000.9 (5 of 261)
Running precomputation on response ENSG00000265491.5 (10 of 261)
Running precomputation on response ENSG00000179085.8 (15 of 261)
Running precomputation on response ENSG00000119185.13 (20 of 261)
Running precomputation on response ENSG00000143952.20 (25 of 261)
Running precomputation on response ENSG00000138385.16 (30 of 261)
Running precomputation on response ENSG00000134086.9 (35 of 261)
Running precomputation on response ENSG00000173473.11 (40 of 261)
Running precomputation on response ENSG00000144959.11 (45 of 261)
Running precomputation on response ENSG00000128052.10 (50 of 261)
Running precomputation on response ENSG00000153147.6 (55 of 261)
Running precomputation on response ENSG00000113360.17 (60 of 261)
Running precomputation on response ENSG00000113658.18 (65 of 261)
Running precomputation on response ENSG00000180822.12 (70 of 261)
Running precomputation on response ENSG00000204371.12 (75 of 261)
Running precomputation on response ENSG00000112249.14 (80 of 261)
Running precomputation on response ENSG00000112592.15 (85 of 261)
Running precomputation on response ENSG00000242265.6 (90 of 261)
Running precomputation on response ENSG00000105887.11 (95 of 261)
Running precomputation on response ENSG00000104412.8 (100 of 261)
Running precomputation on response ENSG00000107036.12 (105 of 261)
Running precomputation on response ENSG00000136925.15 (110 of 261)
Running precomputation on response ENSG00000136908.18 (115 of 261)
Running precomputation on response ENSG00000151461.20 (120 of 261)
Running precomputation on response ENSG00000107779.14 (125 of 261)
Running precomputation on response ENSG00000198730.9 (130 of 261)
Running precomputation on response ENSG00000149136.9 (135 of 261)
Running precomputation on response ENSG00000162227.8 (140 of 261)
Running precomputation on response ENSG00000173933.21 (145 of 261)
Running precomputation on response ENSG00000182934.12 (150 of 261)
Running precomputation on response ENSG00000174718.12 (155 of 261)
Running precomputation on response ENSG00000197111.16 (160 of 261)
Running precomputation on response ENSG00000120868.14 (165 of 261)
Running precomputation on response ENSG00000165480.16 (170 of 261)
Running precomputation on response ENSG00000100462.16 (175 of 261)
Running precomputation on response ENSG00000151748.15 (180 of 261)
Running precomputation on response ENSG00000255529.9 (185 of 261)
Running precomputation on response ENSG00000169371.14 (190 of 261)
Running precomputation on response ENSG00000167978.18 (195 of 261)
Running precomputation on response ENSG00000168488.19 (200 of 261)
Running precomputation on response ENSG00000102900.13 (205 of 261)
Running precomputation on response ENSG00000108961.14 (210 of 261)
Running precomputation on response ENSG00000125686.12 (215 of 261)
Running precomputation on response ENSG00000154217.16 (220 of 261)
Running precomputation on response ENSG00000005206.17 (225 of 261)
Running precomputation on response ENSG00000011451.21 (230 of 261)
Running precomputation on response ENSG00000105568.19 (235 of 261)
Running precomputation on response ENSG00000124181.15 (240 of 261)
Running precomputation on response ENSG00000159086.15 (245 of 261)
Running precomputation on response ENSG00000186575.19 (250 of 261)
Running precomputation on response ENSG00000086758.17 (255 of 261)
Running precomputation on response ENSG00000071553.18 (260 of 261)
Analyzing pairs containing gRNA group KDM1A (1 of 261)
Analyzing pairs containing gRNA group RLF (5 of 261)
Analyzing pairs containing gRNA group RNF115 (10 of 261)
Analyzing pairs containing gRNA group DPM3 (15 of 261)
Analyzing pairs containing gRNA group ITGB1BP1 (20 of 261)
Analyzing pairs containing gRNA group VPS54 (25 of 261)
Analyzing pairs containing gRNA group SSB (30 of 261)
Analyzing pairs containing gRNA group VHL (35 of 261)
Analyzing pairs containing gRNA group SMARCC1 (40 of 261)
Analyzing pairs containing gRNA group NCEH1 (45 of 261)
Analyzing pairs containing gRNA group KDR (50 of 261)
Analyzing pairs containing gRNA group SMARCA5 (55 of 261)
Analyzing pairs containing gRNA group DROSHA (60 of 261)
Analyzing pairs containing gRNA group SMAD5 (65 of 261)
Analyzing pairs containing gRNA group PSMG4 (70 of 261)
Analyzing pairs containing gRNA group EHMT2 (75 of 261)
Analyzing pairs containing gRNA group ASCC3 (80 of 261)
Analyzing pairs containing gRNA group TBP (85 of 261)
Analyzing pairs containing gRNA group PEG10 (90 of 261)
Analyzing pairs containing gRNA group MTPN (95 of 261)
Analyzing pairs containing gRNA group EMC2 (100 of 261)
Analyzing pairs containing gRNA group RIC1 (105 of 261)
Analyzing pairs containing gRNA group TSTD2 (110 of 261)
Analyzing pairs containing gRNA group DPM2 (115 of 261)
Analyzing pairs containing gRNA group UPF2 (120 of 261)
Analyzing pairs containing gRNA group BMPR1A (125 of 261)
Analyzing pairs containing gRNA group CTR9 (130 of 261)
Analyzing pairs containing gRNA group SSRP1 (135 of 261)
Analyzing pairs containing gRNA group TAF6L (140 of 261)
Analyzing pairs containing gRNA group RBM4 (145 of 261)
Analyzing pairs containing gRNA group SRPRA (150 of 261)
Analyzing pairs containing gRNA group RESF1 (155 of 261)
Analyzing pairs containing gRNA group PCBP2 (160 of 261)
Analyzing pairs containing gRNA group APAF1 (165 of 261)
Analyzing pairs containing gRNA group SKA3 (170 of 261)
Analyzing pairs containing gRNA group PRMT5 (175 of 261)
Analyzing pairs containing gRNA group SAV1 (180 of 261)
Analyzing pairs containing gRNA group POLR2M (185 of 261)
Analyzing pairs containing gRNA group SNUPN (190 of 261)
Analyzing pairs containing gRNA group SRRM2 (195 of 261)
Analyzing pairs containing gRNA group ATXN2L (200 of 261)
Analyzing pairs containing gRNA group NUP93 (205 of 261)
Analyzing pairs containing gRNA group RANGRF (210 of 261)
Analyzing pairs containing gRNA group MED1 (215 of 261)
Analyzing pairs containing gRNA group PITPNC1 (220 of 261)
Analyzing pairs containing gRNA group SPPL2B (225 of 261)
Analyzing pairs containing gRNA group WIZ (230 of 261)
Analyzing pairs containing gRNA group PPP2R1A (235 of 261)
Analyzing pairs containing gRNA group PLCG1 (240 of 261)
Analyzing pairs containing gRNA group PAXBP1 (245 of 261)
Analyzing pairs containing gRNA group NF2 (250 of 261)
Analyzing pairs containing gRNA group HUWE1 (255 of 261)
Analyzing pairs containing gRNA group ATP6AP1 (260 of 261)
d0: wrote /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d0_sceptre_transfer/TargetKD_v0.10.3_with_SE.csv

d1: running SCEPTRE from transfer inputs.
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178488.23-P1P2 (1 of 1926)
Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178620.23-P1P2 (5 of 1926)
Performing gRNA-to-cell assignments for gRNA ACTL7B_-_111617995.23-P1P2 (10 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR1B_+_52345762.23-P1P2 (15 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2A_+_148602343.23-P1P2 (20 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_+_38495339.23-P1P2 (25 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_-_38495819.23-P1P2 (30 of 1926)
Performing gRNA-to-cell assignments for gRNA ADCY2_-_7396429.23-P1P2 (35 of 1926)
Performing gRNA-to-cell assignments for gRNA AMOTL2_+_134093279.23-P1P2 (40 of 1926)
Performing gRNA-to-cell assignments for gRNA ANP32E_+_150208342.23-P1P2 (45 of 1926)
Performing gRNA-to-cell assignments for gRNA APAF1_-_99039071.23-P1P2 (50 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_+_72766725.23-P1P2 (55 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_-_72766755.23-P1P2 (60 of 1926)
Performing gRNA-to-cell assignments for gRNA ARTN_-_44399429.23-P1 (65 of 1926)
Performing gRNA-to-cell assignments for gRNA ASCC3_-_101329193.23-P1P2 (70 of 1926)
Performing gRNA-to-cell assignments for gRNA ASH2L_-_37963059.23-P1P2 (75 of 1926)
Performing gRNA-to-cell assignments for gRNA ATMIN_+_81069661.23-P1P2 (80 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_+_111992112.23-P1P2 (85 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_-_111992129.23-P1P2 (90 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP6AP1_-_153657075.23-P1P2 (95 of 1926)
Performing gRNA-to-cell assignments for gRNA ATXN2L_-_28834668.23-P1P2 (100 of 1926)
Performing gRNA-to-cell assignments for gRNA BAD_+_64052102.23-P1P2 (105 of 1926)
Performing gRNA-to-cell assignments for gRNA BAG6_+_31620398.23-P1P2 (110 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_+_49458268.23-P1P2 (115 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_-_49458323.23-P1P2 (120 of 1926)
Performing gRNA-to-cell assignments for gRNA BCAS2_-_115124045.23-P1P2 (125 of 1926)
Performing gRNA-to-cell assignments for gRNA BMPR1A_-_88516386.23-P1P2 (130 of 1926)
Performing gRNA-to-cell assignments for gRNA BRD3_-_136932911.23-P1P2 (135 of 1926)
Performing gRNA-to-cell assignments for gRNA C17orf49_+_6918114.23-P2 (140 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_+_31532436.23-P1P2 (145 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_-_31532650.23-P1P2 (150 of 1926)
Performing gRNA-to-cell assignments for gRNA CACFD1_-_136325160.23-P1P2 (155 of 1926)
Performing gRNA-to-cell assignments for gRNA CACNG1_+_65040756.23-P1P2 (160 of 1926)
Performing gRNA-to-cell assignments for gRNA CAPN5_+_76778225.23-P1P2 (165 of 1926)
Performing gRNA-to-cell assignments for gRNA CASP3_+_185570560.23-P1P2 (170 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_+_107384793.23-P1P2 (175 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_-_107384635.23-P1P2 (180 of 1926)
Performing gRNA-to-cell assignments for gRNA CCM2_-_45039941.23-P1P2 (185 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNC_-_100016267.23-P1P2 (190 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNT2_-_135676416.23-P1P2 (195 of 1926)
Performing gRNA-to-cell assignments for gRNA CCT6B_+_33288381.23-P1P2 (200 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_+_30366413.23-P1P2 (205 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_-_30366630.23-P1P2 (210 of 1926)
Performing gRNA-to-cell assignments for gRNA CD99L2_-_150067120.23-P1P2 (215 of 1926)
Performing gRNA-to-cell assignments for gRNA CDH5_-_66400630.23-P1P2 (220 of 1926)
Performing gRNA-to-cell assignments for gRNA CDKN1A_-_36646538.23-P1P2 (225 of 1926)
Performing gRNA-to-cell assignments for gRNA CENPC_+_68411190.23-P1P2 (230 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_+_30773828.23-P1P2 (235 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_-_30773391.23-P1P2 (240 of 1926)
Performing gRNA-to-cell assignments for gRNA CHMP3_-_86790451.23-P1P2 (245 of 1926)
Performing gRNA-to-cell assignments for gRNA CHORDC1_-_89955894.23-P1P2 (250 of 1926)
Performing gRNA-to-cell assignments for gRNA CNEP1R1_+_50059284.23-P1P2 (255 of 1926)
Performing gRNA-to-cell assignments for gRNA CNOT2_+_70637651.23-P1P2 (260 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_+_83956339.23-P1P2 (265 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_-_83956572.23-P1P2 (270 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS6_-_99686676.23-P1P2 (275 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_-_31634088.23-P1P2 (280 of 1926)
Performing gRNA-to-cell assignments for gRNA CTR9_+_10772881.23-P1P2 (285 of 1926)
Performing gRNA-to-cell assignments for gRNA CUL3_-_225449667.23-P1P2 (290 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_+_145150267.23-P1P2 (295 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_-_145149996.23-P1P2 (300 of 1926)
Performing gRNA-to-cell assignments for gRNA CYCS_+_25164869.23-P1P2 (305 of 1926)
Performing gRNA-to-cell assignments for gRNA DACH1_-_72441101.23-P1P2 (310 of 1926)
Performing gRNA-to-cell assignments for gRNA DBR1_+_137893673.23-P1P2 (315 of 1926)
Performing gRNA-to-cell assignments for gRNA DCTN2_-_57940910.23-P1P2 (320 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_+_61851532.23-P1P2 (325 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_-_61851812.23-P1P2 (330 of 1926)
Performing gRNA-to-cell assignments for gRNA DGCR8_-_20068158.23-P1P2 (335 of 1926)
Performing gRNA-to-cell assignments for gRNA DNAJC17_-_41099215.23-P1P2 (340 of 1926)
Performing gRNA-to-cell assignments for gRNA DNM2_+_10829157.23-P1P2 (345 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM1_+_49575018.23-P1P2 (350 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_+_130700056.23-P1P2 (355 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_-_130700076.23-P1P2 (360 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM3_+_155113020.23-P1P2 (365 of 1926)
Performing gRNA-to-cell assignments for gRNA DROSHA_-_31531982.23-P1P2 (370 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_+_3985081.23-P1 (375 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_-_3985443.23-P1 (380 of 1926)
Performing gRNA-to-cell assignments for gRNA EHMT2_-_31865199.23-P1P2 (385 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3F_+_8009347.23-P1P2 (390 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3M_-_32605405.23-P1P2 (395 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4E_+_99849831.23-P1P2 (400 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184032364.23-P1P2 (405 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184033042.23-P1P2 (410 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G2_+_10830447.23-P1 (415 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF5A_-_7211217.23-P1P2 (420 of 1926)
Performing gRNA-to-cell assignments for gRNA TCEB3_-_24069975.23-P1P2 (425 of 1926)
Performing gRNA-to-cell assignments for gRNA EMC2_+_109456080.23-P1P2 (430 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823220.23-P1P2 (435 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823252.23-P1P2 (440 of 1926)
Performing gRNA-to-cell assignments for gRNA EOMES_-_27763342.23-P1P2 (445 of 1926)
Performing gRNA-to-cell assignments for gRNA ERCC2_+_45873802.23-P1P2 (450 of 1926)
Performing gRNA-to-cell assignments for gRNA ERO1L_+_53162425.23-P1P2 (455 of 1926)
Performing gRNA-to-cell assignments for gRNA ETS1_+_128392044.23-P1P2 (460 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_+_36134592.23-P1P2 (465 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_-_36134789.23-P1P2 (470 of 1926)
Performing gRNA-to-cell assignments for gRNA EWSR1_-_29664347.23-P1P2 (475 of 1926)
Performing gRNA-to-cell assignments for gRNA EXOSC8_-_37574979.23-P1P2 (480 of 1926)
Performing gRNA-to-cell assignments for gRNA FBXW11_-_171433300.23-P1P2 (485 of 1926)
Performing gRNA-to-cell assignments for gRNA FEN1_+_61560782.23-P1P2 (490 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_+_75532080.23-P1P2 (495 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_-_75532461.23-P1P2 (500 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT7_-_139926495.23-P1P2 (505 of 1926)
Performing gRNA-to-cell assignments for gRNA G3BP1_-_151151499.23-P1P2 (510 of 1926)
Performing gRNA-to-cell assignments for gRNA GCK_+_44229021.23-P1P2 (515 of 1926)
Performing gRNA-to-cell assignments for gRNA GOLT1B_+_21654777.23-P1P2 (520 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_+_45000529.23-P1P2 (525 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_-_45000948.23-P1P2 (530 of 1926)
Performing gRNA-to-cell assignments for gRNA GPI_+_34856284.23-P1P2 (535 of 1926)
Performing gRNA-to-cell assignments for gRNA GPR107_-_132816044.23-P1P2 (540 of 1926)
Performing gRNA-to-cell assignments for gRNA GSPT1_-_12009714.23-P1P2 (545 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A1_+_81687215.23-P1P2 (550 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_+_59949355.23-P1P2 (555 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_-_59949657.23-P1P2 (560 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2H4_-_30876185.23-P1P2 (565 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF3C2_+_27579604.23-P1P2 (570 of 1926)
Performing gRNA-to-cell assignments for gRNA HAND1_+_153857732.23-P1P2 (575 of 1926)
Performing gRNA-to-cell assignments for gRNA HDAC2_-_114292253.23-P1P2 (580 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_+_81399312.23-P1P2 (585 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_-_81399351.23-P1P2 (590 of 1926)
Performing gRNA-to-cell assignments for gRNA HINT2_+_35814995.23-P1P2 (595 of 1926)
Performing gRNA-to-cell assignments for gRNA HNF4A_+_43029971.23-P1P2 (600 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPM_-_8509885.23-P1 (605 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPU_+_245027753.23-P1P2 (610 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_+_43702292.23-P1P2 (615 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_-_43702305.23-P1P2 (620 of 1926)
Performing gRNA-to-cell assignments for gRNA HSPA5_-_128003534.23-P1P2 (625 of 1926)
Performing gRNA-to-cell assignments for gRNA HUWE1_-_53710540.23-P1P2 (630 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_-_30007588.23-P1P2 (635 of 1926)
Performing gRNA-to-cell assignments for gRNA INTS10_+_19675000.23-P1P2 (640 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_+_65903046.23-P1P2 (645 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_-_65903523.23-P1P2 (650 of 1926)
Performing gRNA-to-cell assignments for gRNA C7orf26_+_6630358.23-P1P2 (655 of 1926)
Performing gRNA-to-cell assignments for gRNA IPO8_+_30848771.23-P1P2 (660 of 1926)
Performing gRNA-to-cell assignments for gRNA ITGB1BP1_+_9563543.23-P1P2 (665 of 1926)
Performing gRNA-to-cell assignments for gRNA SSFA2_+_182756743.23-P1P2 (670 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_+_109914946.23-P1P2 (675 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_-_109915112.23-P1P2 (680 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM1A_-_23346409.23-P1P2 (685 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM3B_-_137688431.23-P1P2 (690 of 1926)
Performing gRNA-to-cell assignments for gRNA KDR_-_55991659.23-P1P2 (695 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA0232_+_6784651.23-P1P2 (700 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153706791.23-P1P2 (705 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153707248.23-P1P2 (710 of 1926)
Performing gRNA-to-cell assignments for gRNA LAMTOR4_-_99746921.23-P1P2 (715 of 1926)
Performing gRNA-to-cell assignments for gRNA LARP4_-_50794864.23-P1P2 (720 of 1926)
Performing gRNA-to-cell assignments for gRNA LATS2_-_21635336.23-P1P2 (725 of 1926)
Performing gRNA-to-cell assignments for gRNA LEMD3_+_65563441.23-P1P2 (730 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_+_54960428.23-P1P2 (735 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_-_54960616.23-P1P2 (740 of 1926)
Performing gRNA-to-cell assignments for gRNA LHFPL2_-_77944519.23-P1P2 (745 of 1926)
Performing gRNA-to-cell assignments for gRNA LSM12_-_42143603.23-P1P2 (750 of 1926)
Performing gRNA-to-cell assignments for gRNA MANBAL_+_35918136.23-P1P2 (755 of 1926)
Performing gRNA-to-cell assignments for gRNA MAU2_+_19431732.23-P1P2 (760 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_+_47807678.23-P1P2 (765 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_-_47807877.23-P1P2 (770 of 1926)
Performing gRNA-to-cell assignments for gRNA MEAF6_-_37980319.23-P1P2 (775 of 1926)
Performing gRNA-to-cell assignments for gRNA MED1_-_37607140.23-P1P2 (780 of 1926)
Performing gRNA-to-cell assignments for gRNA MED12_+_70338541.23-P1P2 (785 of 1926)
Performing gRNA-to-cell assignments for gRNA MED15_+_20862189.23-P1 (790 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_+_57479335.23-P1P2 (795 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_-_57479437.23-P1P2 (800 of 1926)
Performing gRNA-to-cell assignments for gRNA MEMO1_-_32264882.23-P1 (805 of 1926)
Performing gRNA-to-cell assignments for gRNA MIXL1_+_226411425.23-P1P2 (810 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL28_+_420487.23-P1P2 (815 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL32_+_42972018.23-P1P2 (820 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_+_179322181.23-P1P2 (825 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_-_179322213.23-P1P2 (830 of 1926)
Performing gRNA-to-cell assignments for gRNA LUZP6_-_135612750.23-P1P2 (835 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA20_-_19998023.23-P1P2 (840 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA25_+_112546538.23-P1P2 (845 of 1926)
Performing gRNA-to-cell assignments for gRNA NCBP1_+_100396082.23-P1P2 (850 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_+_172428630.23-P1P2 (855 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_-_172428841.23-P1P2 (860 of 1926)
Performing gRNA-to-cell assignments for gRNA NEDD8_-_24701363.23-P1P2 (865 of 1926)
Performing gRNA-to-cell assignments for gRNA NF2_+_29999635.23-P1P2 (870 of 1926)
Performing gRNA-to-cell assignments for gRNA NFIB_+_14314559.23-P1P2 (875 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_+_139418329.23-P2 (880 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139418349.23-P2 (885 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139440503.23-P1 (890 of 1926)
Performing gRNA-to-cell assignments for gRNA NUFIP2_+_27621073.23-P1P2 (895 of 1926)
Performing gRNA-to-cell assignments for gRNA NUP93_-_56764090.23-P1P2 (900 of 1926)
Performing gRNA-to-cell assignments for gRNA OGT_+_70752948.23-P1P2 (905 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_+_59383170.23-P1P2 (910 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_-_59383348.23-P1P2 (915 of 1926)
Performing gRNA-to-cell assignments for gRNA PAM16_-_4401236.23-P1P2 (920 of 1926)
Performing gRNA-to-cell assignments for gRNA PAXBP1_-_34143614.23-P1P2 (925 of 1926)
Performing gRNA-to-cell assignments for gRNA PBLD_-_70092436.23-P1 (930 of 1926)
Performing gRNA-to-cell assignments for gRNA PCBP2_+_53835614.23-P2 (935 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_+_94292752.23-P2 (940 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_-_94293027.23-P2 (945 of 1926)
Performing gRNA-to-cell assignments for gRNA PGRMC2_+_129208913.23-P1P2 (950 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_-_178866211.23-P1P2 (955 of 1926)
Performing gRNA-to-cell assignments for gRNA PITPNC1_-_65373433.23-P1P2 (960 of 1926)
Performing gRNA-to-cell assignments for gRNA PLCG1_+_39766176.23-P1P2 (965 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_+_151042273.23-P2 (970 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_-_151042422.23-P2 (975 of 1926)
Performing gRNA-to-cell assignments for gRNA PLRG1_-_155471154.23-P1P2 (980 of 1926)
Performing gRNA-to-cell assignments for gRNA PLXNA2_-_208417275.23-P1P2 (985 of 1926)
Performing gRNA-to-cell assignments for gRNA POLE4_-_75185725.23-P1P2 (990 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR2M_+_57998962.23-P1P2 (995 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940045.23-P1P2 (1000 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940414.23-P1P2 (1005 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3K_-_103144.23-P1P2 (1010 of 1926)
Performing gRNA-to-cell assignments for gRNA PPHLN1_-_42720007.23-P1P2 (1015 of 1926)
Performing gRNA-to-cell assignments for gRNA PPIL4_+_149867146.23-P1P2 (1020 of 1926)
Performing gRNA-to-cell assignments for gRNA PPP2R1A_+_52693317.23-P1 (1025 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_+_53195217.23-P1P2 (1030 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_-_53195332.23-P1P2 (1035 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKD1_-_30396996.23-P1P2 (1040 of 1926)
Performing gRNA-to-cell assignments for gRNA PRMT5_+_23398554.23-P1 (1045 of 1926)
Performing gRNA-to-cell assignments for gRNA C2orf53_-_27362230.23-P1P2 (1050 of 1926)
Performing gRNA-to-cell assignments for gRNA PRSS8_+_31147010.23-P1P2 (1055 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_+_36909062.23-P1P2 (1060 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_-_36909230.23-P1P2 (1065 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMD4_-_151227289.23-P1P2 (1070 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMG4_-_3259227.23-P1P2 (1075 of 1926)
Performing gRNA-to-cell assignments for gRNA PTAR1_+_72374835.23-P1P2 (1080 of 1926)
Performing gRNA-to-cell assignments for gRNA PTBP1_+_797463.23-P1P2 (1085 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_+_89623388.23-P1P2 (1090 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_-_89623526.23-P1P2 (1095 of 1926)
Performing gRNA-to-cell assignments for gRNA PTMS_-_6875591.23-P1P2 (1100 of 1926)
Performing gRNA-to-cell assignments for gRNA PYGO2_-_154934135.23-P1P2 (1105 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB18_-_27793292.23-P1P2 (1110 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB6A_+_73471983.23-P1P2 (1115 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_+_76251918.23-P1P2 (1120 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_-_76252008.23-P1P2 (1125 of 1926)
Performing gRNA-to-cell assignments for gRNA GNB2L1_-_180670840.23-P1P2 (1130 of 1926)
Performing gRNA-to-cell assignments for gRNA RAD21_-_117886635.23-P1P2 (1135 of 1926)
Performing gRNA-to-cell assignments for gRNA RANGRF_+_8192257.23-P1P2 (1140 of 1926)
Performing gRNA-to-cell assignments for gRNA RBBP5_+_205091050.23-P1P2 (1145 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM14-RBM4_+_66406280.23-P2 (1150 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM33_+_155437802.23-P1P2 (1155 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM4_+_66406319.23-P2 (1160 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_+_114271352.23-P1P2 (1165 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_-_114271700.23-P1P2 (1170 of 1926)
Performing gRNA-to-cell assignments for gRNA RBX1_-_41347852.23-P1P2 (1175 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1551_+_32112399.23-P1P2 (1180 of 1926)
Performing gRNA-to-cell assignments for gRNA RFX3_-_3525750.23-P1P2 (1185 of 1926)
Performing gRNA-to-cell assignments for gRNA RGP1_+_35749328.23-P1P2 (1190 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_+_95067102.23-P1P2 (1195 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_-_95067164.23-P1P2 (1200 of 1926)
Performing gRNA-to-cell assignments for gRNA RIBC1_-_53449606.23-P1P2 (1205 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1432_+_5629173.23-P1P2 (1210 of 1926)
Performing gRNA-to-cell assignments for gRNA RLF_+_40627522.23-P1P2 (1215 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF115_-_145611148.23-P1P2 (1220 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_+_104296142.23-P1P2 (1225 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_-_104296250.23-P1P2 (1230 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF40_-_30773681.23-P1P2 (1235 of 1926)
Performing gRNA-to-cell assignments for gRNA RPL37A_-_217363600.23-P1P2 (1240 of 1926)
Performing gRNA-to-cell assignments for gRNA RPS15A_-_18801590.23-P1P2 (1245 of 1926)
Performing gRNA-to-cell assignments for gRNA SAV1_+_51134992.23-P1P2 (1250 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_+_33104038.23-P1P2 (1255 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_-_33104343.23-P1P2 (1260 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF8_+_155054976.23-P1P2 (1265 of 1926)
Performing gRNA-to-cell assignments for gRNA EPT1_+_26569050.23-P1P2 (1270 of 1926)
Performing gRNA-to-cell assignments for gRNA SERF2_+_44069659.23-P2 (1275 of 1926)
Performing gRNA-to-cell assignments for gRNA SERPING1_+_57365137.23-P1P2 (1280 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_+_122242041.23-P1P2 (1285 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_-_122242294.23-P1P2 (1290 of 1926)
Performing gRNA-to-cell assignments for gRNA SETDB1_-_150898815.23-P1P2 (1295 of 1926)
Performing gRNA-to-cell assignments for gRNA SF1_-_64545937.23-P1P2 (1300 of 1926)
Performing gRNA-to-cell assignments for gRNA CCDC101_-_28565268.23-P1P2 (1305 of 1926)
Performing gRNA-to-cell assignments for gRNA SKA3_+_21750519.23-P1P2 (1310 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_+_78591497.23-P1P2 (1315 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_-_78591737.23-P1P2 (1320 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC29A1_-_44187438.23-P2 (1325 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC35C1_-_45826711.23-P1P2 (1330 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD1_+_146403428.23-P1P2 (1335 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD3_+_67358229.23-P1P2 (1340 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_+_48556623.23-P1P2 (1345 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_-_48556657.23-P1P2 (1350 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD5_-_135468694.23-P1P2 (1355 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA4_-_11071741.23-P1P2 (1360 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA5_+_144435243.23-P1P2 (1365 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCB1_+_24129228.23-P1 (1370 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_+_47822975.23-P1P2 (1375 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_-_47822893.23-P1P2 (1380 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCD3_-_150974212.23-P1 (1385 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCE1_-_38803669.23-P1P2 (1390 of 1926)
Performing gRNA-to-cell assignments for gRNA SMC2_-_106856631.23-P1P2 (1395 of 1926)
Performing gRNA-to-cell assignments for gRNA SMG5_+_156252537.23-P1P2 (1400 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_+_103852.23-P1P2 (1405 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_-_104242.23-P1P2 (1410 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP27_-_70121181.23-P1P2 (1415 of 1926)
Performing gRNA-to-cell assignments for gRNA SNUPN_-_75917960.23-P1P2 (1420 of 1926)
Performing gRNA-to-cell assignments for gRNA SOX4_+_21594213.23-P1 (1425 of 1926)
Performing gRNA-to-cell assignments for gRNA SPIN1_+_91003737.23-P1P2 (1430 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_+_2328680.23-P1P2 (1435 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_-_2328912.23-P1P2 (1440 of 1926)
Performing gRNA-to-cell assignments for gRNA SPTY2D1_-_18655846.23-P1P2 (1445 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP19_+_112197188.23-P1P2 (1450 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP54_-_35452244.23-P1P2 (1455 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP68_+_74068540.23-P1P2 (1460 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_+_126138730.23-P1P2 (1465 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_-_126138742.23-P1P2 (1470 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPRB_-_133524719.23-P1P2 (1475 of 1926)
Performing gRNA-to-cell assignments for gRNA SRRM2_-_2802714.23-P1P2 (1480 of 1926)
Performing gRNA-to-cell assignments for gRNA SRSF6_+_42086993.23-P1 (1485 of 1926)
Performing gRNA-to-cell assignments for gRNA SSB_+_170655450.23-P1P2 (1490 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_+_57102902.23-P1P2 (1495 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_-_57103355.23-P1P2 (1500 of 1926)
Performing gRNA-to-cell assignments for gRNA STAG2_-_123095658.23-P1 (1505 of 1926)
Performing gRNA-to-cell assignments for gRNA STK40_+_36851507.23-P1P2 (1510 of 1926)
Performing gRNA-to-cell assignments for gRNA SUPT20H_+_37633798.23-P1P2 (1515 of 1926)
Performing gRNA-to-cell assignments for gRNA SUZ12_+_30264157.23-P1P2 (1520 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_+_64901973.23-P1P2 (1525 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_-_64901981.23-P1P2 (1530 of 1926)
Performing gRNA-to-cell assignments for gRNA TADA2B_-_7045121.23-P1P2 (1535 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF10_+_6633406.23-P1P2 (1540 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF2_+_120845032.23-P1P2 (1545 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF3_-_7860540.23-P1P2 (1550 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_+_60640814.23-P1P2 (1555 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_-_60640498.23-P1P2 (1560 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF6L_-_62538954.23-P1P2 (1565 of 1926)
Performing gRNA-to-cell assignments for gRNA TAOK1_-_27717907.23-P1P2 (1570 of 1926)
Performing gRNA-to-cell assignments for gRNA FAM208A_+_56717273.23-P1P2 (1575 of 1926)
Performing gRNA-to-cell assignments for gRNA TBC1D9_+_141677481.23-P1P2 (1580 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_+_170863524.23-P1P2 (1585 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_-_170863481.23-P1P2 (1590 of 1926)
Performing gRNA-to-cell assignments for gRNA TBX3_-_115121998.23-P1P2 (1595 of 1926)
Performing gRNA-to-cell assignments for gRNA T_-_166581333.23-P1P2 (1600 of 1926)
Performing gRNA-to-cell assignments for gRNA TGFBR2_+_30648526.23-P1P2 (1605 of 1926)
Performing gRNA-to-cell assignments for gRNA THBD_+_23029196.23-P2 (1610 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_+_2977592.23-P1P2 (1615 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_-_2977643.23-P1P2 (1620 of 1926)
Performing gRNA-to-cell assignments for gRNA MESDC1_-_81293340.23-P1P2 (1625 of 1926)
Performing gRNA-to-cell assignments for gRNA TMEM41B_+_9336156.23-P1P2 (1630 of 1926)
Performing gRNA-to-cell assignments for gRNA TMTC3_-_88536127.23-P1P2 (1635 of 1926)
Performing gRNA-to-cell assignments for gRNA TOP1_+_39657538.23-P1P2 (1640 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_+_7590531.23-P1P2 (1645 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_-_7590704.23-P1P2 (1650 of 1926)
Performing gRNA-to-cell assignments for gRNA TRIM24_-_138145308.23-P1P2 (1655 of 1926)
Performing gRNA-to-cell assignments for gRNA TRRAP_-_98476465.23-P1P2 (1660 of 1926)
Performing gRNA-to-cell assignments for gRNA TSTD2_+_100395594.23-P1P2 (1665 of 1926)
Performing gRNA-to-cell assignments for gRNA UBE2E1_+_23847541.23-P1P2 (1670 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_+_38924011.23-P1P2 (1675 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_-_38924421.23-P1P2 (1680 of 1926)
Performing gRNA-to-cell assignments for gRNA UPF2_+_12084799.23-P1P2 (1685 of 1926)
Performing gRNA-to-cell assignments for gRNA UXT_-_47518505.23-P1P2 (1690 of 1926)
Performing gRNA-to-cell assignments for gRNA VHL_-_10183525.23-P1P2 (1695 of 1926)
Performing gRNA-to-cell assignments for gRNA VPS54_+_64246104.23-P1P2 (1700 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_+_28821698.23-P1P2 (1705 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_-_28822074.23-P1P2 (1710 of 1926)
Performing gRNA-to-cell assignments for gRNA WAPAL_+_88281533.23-P1P2 (1715 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR12_-_203776383.23-P1P2 (1720 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR77_+_111991860.23-P1P2 (1725 of 1926)
Performing gRNA-to-cell assignments for gRNA WIZ_-_15543016.23-P1P2 (1730 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_+_9983328.23-P1P2 (1735 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_-_9983510.23-P1P2 (1740 of 1926)
Performing gRNA-to-cell assignments for gRNA XPO5_-_43543670.23-P1P2 (1745 of 1926)
Performing gRNA-to-cell assignments for gRNA XRN1_-_142166648.23-P1P2 (1750 of 1926)
Performing gRNA-to-cell assignments for gRNA YARS2_+_32908787.23-P1P2 (1755 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB22_+_33285454.23-P1P2 (1760 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_+_6640093.23-P1P2 (1765 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_-_6640138.23-P1P2 (1770 of 1926)
Performing gRNA-to-cell assignments for gRNA ZC3H8_-_113012435.23-P1P2 (1775 of 1926)
Performing gRNA-to-cell assignments for gRNA ZNF281_+_200379151.23-P1P2 (1780 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00117 (1785 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00233 (1790 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00518 (1795 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00656 (1800 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01226 (1805 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01730 (1810 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02009 (1815 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02288 (1820 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02777 (1825 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03026 (1830 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03252 (1835 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03620 (1840 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_SGTN_04 (1845 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_SGTN_25 (1850 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_SGTN_26 (1855 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1044 (1860 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1414 (1865 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1554 (1870 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1785 (1875 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2009 (1880 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2271 (1885 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2425 (1890 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2614 (1895 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_293 (1900 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_541 (1905 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_725 (1910 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_860 (1915 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_SGTN_14 (1920 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_SGTN_01 (1925 of 1926)
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Constructing negative control pairs. ✓
Running precomputation on response ENSG00000197747.9 (1 of 242)
Running precomputation on response ENSG00000103035.11 (5 of 242)
Running precomputation on response ENSG00000077150.21 (10 of 242)
Running precomputation on response ENSG00000175198.17 (15 of 242)
Running precomputation on response ENSG00000236859.8 (20 of 242)
Running precomputation on response ENSG00000176171.12 (25 of 242)
Running precomputation on response ENSG00000162733.19 (30 of 242)
Running precomputation on response ENSG00000269547.1 (35 of 242)
Running precomputation on response ENSG00000145996.11 (40 of 242)
Running precomputation on response ENSG00000116141.17 (45 of 242)
Running precomputation on response ENSG00000166025.18 (50 of 242)
Running precomputation on response ENSG00000089775.12 (55 of 242)
Running precomputation on response ENSG00000099204.22 (60 of 242)
Running precomputation on response ENSG00000157224.16 (65 of 242)
Running precomputation on response ENSG00000173207.13 (70 of 242)
Running precomputation on response ENSG00000008869.12 (75 of 242)
Running precomputation on response ENSG00000111145.8 (80 of 242)
Running precomputation on response ENSG00000137815.15 (85 of 242)
Running precomputation on response ENSG00000213614.11 (90 of 242)
Running precomputation on response ENSG00000048342.18 (95 of 242)
Running precomputation on response ENSG00000197302.11 (100 of 242)
Running precomputation on response ENSG00000065243.20 (105 of 242)
Running precomputation on response ENSG00000080546.14 (110 of 242)
Running precomputation on response ENSG00000138835.24 (115 of 242)
Running precomputation on response ENSG00000118971.9 (120 of 242)
Running precomputation on response ENSG00000182287.15 (125 of 242)
Running precomputation on response ENSG00000143786.8 (130 of 242)
Running precomputation on response ENSG00000034713.8 (135 of 242)
Running precomputation on response ENSG00000133874.2 (140 of 242)
Running precomputation on response ENSG00000101413.12 (145 of 242)
Running precomputation on response ENSG00000269609.6 (150 of 242)
Running precomputation on response ENSG00000067829.20 (155 of 242)
Running precomputation on response ENSG00000269955.2 (160 of 242)
Running precomputation on response ENSG00000121897.15 (165 of 242)
Running precomputation on response ENSG00000110651.13 (170 of 242)
Running precomputation on response ENSG00000030419.17 (175 of 242)
Running precomputation on response ENSG00000136866.14 (180 of 242)
Running precomputation on response ENSG00000181830.9 (185 of 242)
Running precomputation on response ENSG00000176681.15 (190 of 242)
Running precomputation on response ENSG00000148730.7 (195 of 242)
Running precomputation on response ENSG00000069869.17 (200 of 242)
Running precomputation on response ENSG00000109670.17 (205 of 242)
Running precomputation on response ENSG00000107521.21 (210 of 242)
Running precomputation on response ENSG00000086015.23 (215 of 242)
Running precomputation on response ENSG00000112343.11 (220 of 242)
Running precomputation on response ENSG00000181027.11 (225 of 242)
Running precomputation on response ENSG00000171314.9 (230 of 242)
Running precomputation on response ENSG00000044090.13 (235 of 242)
Running precomputation on response ENSG00000143727.16 (240 of 242)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00259&non-targeting_02159&non-targeting_02571&non-targeting_03029&non-targeting_03081 (1 of 92)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00355&non-targeting_02025&non-targeting_02089&non-targeting_03081&non-targeting_03384 (5 of 92)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_01730&non-targeting_02193&non-targeting_03180&non-targeting_03620&non-targeting_03664 (10 of 92)
Analyzing pairs containing gRNA group non-targeting_00043&non-targeting_00184&non-targeting_02677&non-targeting_02777&non-targeting_02814&non-targeting_03427 (15 of 92)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00355&non-targeting_00390&non-targeting_00918&non-targeting_01509&non-targeting_01599 (20 of 92)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00972&non-targeting_02193&non-targeting_02574&non-targeting_02777&non-targeting_02962 (25 of 92)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00355&non-targeting_01601&non-targeting_03026&non-targeting_03081&non-targeting_03158 (30 of 92)
Analyzing pairs containing gRNA group non-targeting_00231&non-targeting_00518&non-targeting_00629&non-targeting_00794&non-targeting_02777&non-targeting_02837 (35 of 92)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00390&non-targeting_01601&non-targeting_02009&non-targeting_02288&non-targeting_03029 (40 of 92)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00231&non-targeting_00324&non-targeting_02025&non-targeting_03252&non-targeting_03384 (45 of 92)
Analyzing pairs containing gRNA group non-targeting_00043&non-targeting_00213&non-targeting_00259&non-targeting_00355&non-targeting_00614&non-targeting_01934 (50 of 92)
Analyzing pairs containing gRNA group non-targeting_00629&non-targeting_01185&non-targeting_01732&non-targeting_02089&non-targeting_02677&non-targeting_03483 (55 of 92)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00518&non-targeting_00614&non-targeting_02025&non-targeting_03342&non-targeting_03427 (60 of 92)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00231&non-targeting_00656&non-targeting_02025&non-targeting_02288&non-targeting_02777 (65 of 92)
Analyzing pairs containing gRNA group non-targeting_00614&non-targeting_00972&non-targeting_02911&non-targeting_03081&non-targeting_03180&non-targeting_03384 (70 of 92)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_00555&non-targeting_01934&non-targeting_02193&non-targeting_02677&non-targeting_03252 (75 of 92)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00656&non-targeting_00918&non-targeting_01934&non-targeting_02574&non-targeting_02814 (80 of 92)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_00355&non-targeting_02837&non-targeting_02911&non-targeting_03342&non-targeting_03620 (85 of 92)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00518&non-targeting_00614&non-targeting_01226&non-targeting_02025&non-targeting_02837 (90 of 92)
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Running precomputation on response ENSG00000004487.18 (1 of 243)
Running precomputation on response ENSG00000117000.9 (5 of 243)
Running precomputation on response ENSG00000265491.5 (10 of 243)
Running precomputation on response ENSG00000198952.8 (15 of 243)
Running precomputation on response ENSG00000153187.21 (20 of 243)
Running precomputation on response ENSG00000162961.14 (25 of 243)
Running precomputation on response ENSG00000121989.15 (30 of 243)
Running precomputation on response ENSG00000036257.14 (35 of 243)
Running precomputation on response ENSG00000114739.14 (40 of 243)
Running precomputation on response ENSG00000114127.11 (45 of 243)
Running precomputation on response ENSG00000128052.10 (50 of 243)
Running precomputation on response ENSG00000153147.6 (55 of 243)
Running precomputation on response ENSG00000113360.17 (60 of 243)
Running precomputation on response ENSG00000113658.18 (65 of 243)
Running precomputation on response ENSG00000180822.12 (70 of 243)
Running precomputation on response ENSG00000124571.18 (75 of 243)
Running precomputation on response ENSG00000131013.4 (80 of 243)
Running precomputation on response ENSG00000106591.4 (85 of 243)
Running precomputation on response ENSG00000188186.11 (90 of 243)
Running precomputation on response ENSG00000104613.12 (95 of 243)
Running precomputation on response ENSG00000179091.5 (100 of 243)
Running precomputation on response ENSG00000106723.17 (105 of 243)
Running precomputation on response ENSG00000136908.18 (110 of 243)
Running precomputation on response ENSG00000151461.20 (115 of 243)
Running precomputation on response ENSG00000107779.14 (120 of 243)
Running precomputation on response ENSG00000198730.9 (125 of 243)
Running precomputation on response ENSG00000149136.9 (130 of 243)
Running precomputation on response ENSG00000162298.19 (135 of 243)
Running precomputation on response ENSG00000076053.12 (140 of 243)
Running precomputation on response ENSG00000133704.10 (145 of 243)
Running precomputation on response ENSG00000197111.16 (150 of 243)
Running precomputation on response ENSG00000110906.13 (155 of 243)
Running precomputation on response ENSG00000165480.16 (160 of 243)
Running precomputation on response ENSG00000100462.16 (165 of 243)
Running precomputation on response ENSG00000165417.12 (170 of 243)
Running precomputation on response ENSG00000166949.17 (175 of 243)
Running precomputation on response ENSG00000161981.11 (180 of 243)
Running precomputation on response ENSG00000134419.15 (185 of 243)
Running precomputation on response ENSG00000052344.16 (190 of 243)
Running precomputation on response ENSG00000141510.18 (195 of 243)
Running precomputation on response ENSG00000125686.12 (200 of 243)
Running precomputation on response ENSG00000154217.16 (205 of 243)
Running precomputation on response ENSG00000005206.17 (210 of 243)
Running precomputation on response ENSG00000011451.21 (215 of 243)
Running precomputation on response ENSG00000105568.19 (220 of 243)
Running precomputation on response ENSG00000124193.16 (225 of 243)
Running precomputation on response ENSG00000128191.17 (230 of 243)
Running precomputation on response ENSG00000100387.9 (235 of 243)
Running precomputation on response ENSG00000147162.15 (240 of 243)
Analyzing pairs containing gRNA group KDM1A (1 of 243)
Analyzing pairs containing gRNA group RLF (5 of 243)
Analyzing pairs containing gRNA group RNF115 (10 of 243)
Analyzing pairs containing gRNA group SMG5 (15 of 243)
Analyzing pairs containing gRNA group HNRNPU (20 of 243)
Analyzing pairs containing gRNA group DPY30 (25 of 243)
Analyzing pairs containing gRNA group ACVR2A (30 of 243)
Analyzing pairs containing gRNA group CUL3 (35 of 243)
Analyzing pairs containing gRNA group ACVR2B (40 of 243)
Analyzing pairs containing gRNA group XRN1 (45 of 243)
Analyzing pairs containing gRNA group KDR (50 of 243)
Analyzing pairs containing gRNA group SMARCA5 (55 of 243)
Analyzing pairs containing gRNA group DROSHA (60 of 243)
Analyzing pairs containing gRNA group SMAD5 (65 of 243)
Analyzing pairs containing gRNA group PSMG4 (70 of 243)
Analyzing pairs containing gRNA group XPO5 (75 of 243)
Analyzing pairs containing gRNA group PPIL4 (80 of 243)
Analyzing pairs containing gRNA group MRPL32 (85 of 243)
Analyzing pairs containing gRNA group LAMTOR4 (90 of 243)
Analyzing pairs containing gRNA group INTS10 (95 of 243)
Analyzing pairs containing gRNA group CYC1 (100 of 243)
Analyzing pairs containing gRNA group SPIN1 (105 of 243)
Analyzing pairs containing gRNA group DPM2 (110 of 243)
Analyzing pairs containing gRNA group UPF2 (115 of 243)
Analyzing pairs containing gRNA group BMPR1A (120 of 243)
Analyzing pairs containing gRNA group CTR9 (125 of 243)
Analyzing pairs containing gRNA group SSRP1 (130 of 243)
Analyzing pairs containing gRNA group SYVN1 (135 of 243)
Analyzing pairs containing gRNA group RBM7 (140 of 243)
Analyzing pairs containing gRNA group IPO8 (145 of 243)
Analyzing pairs containing gRNA group PCBP2 (150 of 243)
Analyzing pairs containing gRNA group KCTD10 (155 of 243)
Analyzing pairs containing gRNA group SKA3 (160 of 243)
Analyzing pairs containing gRNA group PRMT5 (165 of 243)
Analyzing pairs containing gRNA group GTF2A1 (170 of 243)
Analyzing pairs containing gRNA group SMAD3 (175 of 243)
Analyzing pairs containing gRNA group SNRNP25 (180 of 243)
Analyzing pairs containing gRNA group RPS15A (185 of 243)
Analyzing pairs containing gRNA group PRSS8 (190 of 243)
Analyzing pairs containing gRNA group TP53 (195 of 243)
Analyzing pairs containing gRNA group MED1 (200 of 243)
Analyzing pairs containing gRNA group PITPNC1 (205 of 243)
Analyzing pairs containing gRNA group SPPL2B (210 of 243)
Analyzing pairs containing gRNA group WIZ (215 of 243)
Analyzing pairs containing gRNA group PPP2R1A (220 of 243)
Analyzing pairs containing gRNA group SRSF6 (225 of 243)
Analyzing pairs containing gRNA group DGCR8 (230 of 243)
Analyzing pairs containing gRNA group RBX1 (235 of 243)
Analyzing pairs containing gRNA group OGT (240 of 243)
d1: wrote /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d1_sceptre_transfer/TargetKD_v0.10.3_with_SE.csv

d2: running SCEPTRE from transfer inputs.
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178488.23-P1P2 (1 of 1926)
Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178620.23-P1P2 (5 of 1926)
Performing gRNA-to-cell assignments for gRNA ACTL7B_-_111617995.23-P1P2 (10 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR1B_+_52345762.23-P1P2 (15 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2A_+_148602343.23-P1P2 (20 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_+_38495339.23-P1P2 (25 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_-_38495819.23-P1P2 (30 of 1926)
Performing gRNA-to-cell assignments for gRNA ADCY2_-_7396429.23-P1P2 (35 of 1926)
Performing gRNA-to-cell assignments for gRNA AMOTL2_+_134093279.23-P1P2 (40 of 1926)
Performing gRNA-to-cell assignments for gRNA ANP32E_+_150208342.23-P1P2 (45 of 1926)
Performing gRNA-to-cell assignments for gRNA APAF1_-_99039071.23-P1P2 (50 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_+_72766725.23-P1P2 (55 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_-_72766755.23-P1P2 (60 of 1926)
Performing gRNA-to-cell assignments for gRNA ARTN_-_44399429.23-P1 (65 of 1926)
Performing gRNA-to-cell assignments for gRNA ASCC3_-_101329193.23-P1P2 (70 of 1926)
Performing gRNA-to-cell assignments for gRNA ASH2L_-_37963059.23-P1P2 (75 of 1926)
Performing gRNA-to-cell assignments for gRNA ATMIN_+_81069661.23-P1P2 (80 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_+_111992112.23-P1P2 (85 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_-_111992129.23-P1P2 (90 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP6AP1_-_153657075.23-P1P2 (95 of 1926)
Performing gRNA-to-cell assignments for gRNA ATXN2L_-_28834668.23-P1P2 (100 of 1926)
Performing gRNA-to-cell assignments for gRNA BAD_+_64052102.23-P1P2 (105 of 1926)
Performing gRNA-to-cell assignments for gRNA BAG6_+_31620398.23-P1P2 (110 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_+_49458268.23-P1P2 (115 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_-_49458323.23-P1P2 (120 of 1926)
Performing gRNA-to-cell assignments for gRNA BCAS2_-_115124045.23-P1P2 (125 of 1926)
Performing gRNA-to-cell assignments for gRNA BMPR1A_-_88516386.23-P1P2 (130 of 1926)
Performing gRNA-to-cell assignments for gRNA BRD3_-_136932911.23-P1P2 (135 of 1926)
Performing gRNA-to-cell assignments for gRNA C17orf49_+_6918114.23-P2 (140 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_+_31532436.23-P1P2 (145 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_-_31532650.23-P1P2 (150 of 1926)
Performing gRNA-to-cell assignments for gRNA CACFD1_-_136325160.23-P1P2 (155 of 1926)
Performing gRNA-to-cell assignments for gRNA CACNG1_+_65040756.23-P1P2 (160 of 1926)
Performing gRNA-to-cell assignments for gRNA CAPN5_+_76778225.23-P1P2 (165 of 1926)
Performing gRNA-to-cell assignments for gRNA CASP3_+_185570560.23-P1P2 (170 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_+_107384793.23-P1P2 (175 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_-_107384635.23-P1P2 (180 of 1926)
Performing gRNA-to-cell assignments for gRNA CCM2_-_45039941.23-P1P2 (185 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNC_-_100016267.23-P1P2 (190 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNT2_-_135676416.23-P1P2 (195 of 1926)
Performing gRNA-to-cell assignments for gRNA CCT6B_+_33288381.23-P1P2 (200 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_+_30366413.23-P1P2 (205 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_-_30366630.23-P1P2 (210 of 1926)
Performing gRNA-to-cell assignments for gRNA CD99L2_-_150067120.23-P1P2 (215 of 1926)
Performing gRNA-to-cell assignments for gRNA CDH5_-_66400630.23-P1P2 (220 of 1926)
Performing gRNA-to-cell assignments for gRNA CDKN1A_-_36646538.23-P1P2 (225 of 1926)
Performing gRNA-to-cell assignments for gRNA CENPC_+_68411190.23-P1P2 (230 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_+_30773828.23-P1P2 (235 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_-_30773391.23-P1P2 (240 of 1926)
Performing gRNA-to-cell assignments for gRNA CHMP3_-_86790451.23-P1P2 (245 of 1926)
Performing gRNA-to-cell assignments for gRNA CHORDC1_-_89955894.23-P1P2 (250 of 1926)
Performing gRNA-to-cell assignments for gRNA CNEP1R1_+_50059284.23-P1P2 (255 of 1926)
Performing gRNA-to-cell assignments for gRNA CNOT2_+_70637651.23-P1P2 (260 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_+_83956339.23-P1P2 (265 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_-_83956572.23-P1P2 (270 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS6_-_99686676.23-P1P2 (275 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_-_31634088.23-P1P2 (280 of 1926)
Performing gRNA-to-cell assignments for gRNA CTR9_+_10772881.23-P1P2 (285 of 1926)
Performing gRNA-to-cell assignments for gRNA CUL3_-_225449667.23-P1P2 (290 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_+_145150267.23-P1P2 (295 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_-_145149996.23-P1P2 (300 of 1926)
Performing gRNA-to-cell assignments for gRNA CYCS_+_25164869.23-P1P2 (305 of 1926)
Performing gRNA-to-cell assignments for gRNA DACH1_-_72441101.23-P1P2 (310 of 1926)
Performing gRNA-to-cell assignments for gRNA DBR1_+_137893673.23-P1P2 (315 of 1926)
Performing gRNA-to-cell assignments for gRNA DCTN2_-_57940910.23-P1P2 (320 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_+_61851532.23-P1P2 (325 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_-_61851812.23-P1P2 (330 of 1926)
Performing gRNA-to-cell assignments for gRNA DGCR8_-_20068158.23-P1P2 (335 of 1926)
Performing gRNA-to-cell assignments for gRNA DNAJC17_-_41099215.23-P1P2 (340 of 1926)
Performing gRNA-to-cell assignments for gRNA DNM2_+_10829157.23-P1P2 (345 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM1_+_49575018.23-P1P2 (350 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_+_130700056.23-P1P2 (355 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_-_130700076.23-P1P2 (360 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM3_+_155113020.23-P1P2 (365 of 1926)
Performing gRNA-to-cell assignments for gRNA DROSHA_-_31531982.23-P1P2 (370 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_+_3985081.23-P1 (375 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_-_3985443.23-P1 (380 of 1926)
Performing gRNA-to-cell assignments for gRNA EHMT2_-_31865199.23-P1P2 (385 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3F_+_8009347.23-P1P2 (390 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3M_-_32605405.23-P1P2 (395 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4E_+_99849831.23-P1P2 (400 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184032364.23-P1P2 (405 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184033042.23-P1P2 (410 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G2_+_10830447.23-P1 (415 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF5A_-_7211217.23-P1P2 (420 of 1926)
Performing gRNA-to-cell assignments for gRNA TCEB3_-_24069975.23-P1P2 (425 of 1926)
Performing gRNA-to-cell assignments for gRNA EMC2_+_109456080.23-P1P2 (430 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823220.23-P1P2 (435 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823252.23-P1P2 (440 of 1926)
Performing gRNA-to-cell assignments for gRNA EOMES_-_27763342.23-P1P2 (445 of 1926)
Performing gRNA-to-cell assignments for gRNA ERCC2_+_45873802.23-P1P2 (450 of 1926)
Performing gRNA-to-cell assignments for gRNA ERO1L_+_53162425.23-P1P2 (455 of 1926)
Performing gRNA-to-cell assignments for gRNA ETS1_+_128392044.23-P1P2 (460 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_+_36134592.23-P1P2 (465 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_-_36134789.23-P1P2 (470 of 1926)
Performing gRNA-to-cell assignments for gRNA EWSR1_-_29664347.23-P1P2 (475 of 1926)
Performing gRNA-to-cell assignments for gRNA EXOSC8_-_37574979.23-P1P2 (480 of 1926)
Performing gRNA-to-cell assignments for gRNA FBXW11_-_171433300.23-P1P2 (485 of 1926)
Performing gRNA-to-cell assignments for gRNA FEN1_+_61560782.23-P1P2 (490 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_+_75532080.23-P1P2 (495 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_-_75532461.23-P1P2 (500 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT7_-_139926495.23-P1P2 (505 of 1926)
Performing gRNA-to-cell assignments for gRNA G3BP1_-_151151499.23-P1P2 (510 of 1926)
Performing gRNA-to-cell assignments for gRNA GCK_+_44229021.23-P1P2 (515 of 1926)
Performing gRNA-to-cell assignments for gRNA GOLT1B_+_21654777.23-P1P2 (520 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_+_45000529.23-P1P2 (525 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_-_45000948.23-P1P2 (530 of 1926)
Performing gRNA-to-cell assignments for gRNA GPI_+_34856284.23-P1P2 (535 of 1926)
Performing gRNA-to-cell assignments for gRNA GPR107_-_132816044.23-P1P2 (540 of 1926)
Performing gRNA-to-cell assignments for gRNA GSPT1_-_12009714.23-P1P2 (545 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A1_+_81687215.23-P1P2 (550 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_+_59949355.23-P1P2 (555 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_-_59949657.23-P1P2 (560 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2H4_-_30876185.23-P1P2 (565 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF3C2_+_27579604.23-P1P2 (570 of 1926)
Performing gRNA-to-cell assignments for gRNA HAND1_+_153857732.23-P1P2 (575 of 1926)
Performing gRNA-to-cell assignments for gRNA HDAC2_-_114292253.23-P1P2 (580 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_+_81399312.23-P1P2 (585 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_-_81399351.23-P1P2 (590 of 1926)
Performing gRNA-to-cell assignments for gRNA HINT2_+_35814995.23-P1P2 (595 of 1926)
Performing gRNA-to-cell assignments for gRNA HNF4A_+_43029971.23-P1P2 (600 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPM_-_8509885.23-P1 (605 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPU_+_245027753.23-P1P2 (610 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_+_43702292.23-P1P2 (615 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_-_43702305.23-P1P2 (620 of 1926)
Performing gRNA-to-cell assignments for gRNA HSPA5_-_128003534.23-P1P2 (625 of 1926)
Performing gRNA-to-cell assignments for gRNA HUWE1_-_53710540.23-P1P2 (630 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_-_30007588.23-P1P2 (635 of 1926)
Performing gRNA-to-cell assignments for gRNA INTS10_+_19675000.23-P1P2 (640 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_+_65903046.23-P1P2 (645 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_-_65903523.23-P1P2 (650 of 1926)
Performing gRNA-to-cell assignments for gRNA C7orf26_+_6630358.23-P1P2 (655 of 1926)
Performing gRNA-to-cell assignments for gRNA IPO8_+_30848771.23-P1P2 (660 of 1926)
Performing gRNA-to-cell assignments for gRNA ITGB1BP1_+_9563543.23-P1P2 (665 of 1926)
Performing gRNA-to-cell assignments for gRNA SSFA2_+_182756743.23-P1P2 (670 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_+_109914946.23-P1P2 (675 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_-_109915112.23-P1P2 (680 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM1A_-_23346409.23-P1P2 (685 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM3B_-_137688431.23-P1P2 (690 of 1926)
Performing gRNA-to-cell assignments for gRNA KDR_-_55991659.23-P1P2 (695 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA0232_+_6784651.23-P1P2 (700 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153706791.23-P1P2 (705 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153707248.23-P1P2 (710 of 1926)
Performing gRNA-to-cell assignments for gRNA LAMTOR4_-_99746921.23-P1P2 (715 of 1926)
Performing gRNA-to-cell assignments for gRNA LARP4_-_50794864.23-P1P2 (720 of 1926)
Performing gRNA-to-cell assignments for gRNA LATS2_-_21635336.23-P1P2 (725 of 1926)
Performing gRNA-to-cell assignments for gRNA LEMD3_+_65563441.23-P1P2 (730 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_+_54960428.23-P1P2 (735 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_-_54960616.23-P1P2 (740 of 1926)
Performing gRNA-to-cell assignments for gRNA LHFPL2_-_77944519.23-P1P2 (745 of 1926)
Performing gRNA-to-cell assignments for gRNA LSM12_-_42143603.23-P1P2 (750 of 1926)
Performing gRNA-to-cell assignments for gRNA MANBAL_+_35918136.23-P1P2 (755 of 1926)
Performing gRNA-to-cell assignments for gRNA MAU2_+_19431732.23-P1P2 (760 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_+_47807678.23-P1P2 (765 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_-_47807877.23-P1P2 (770 of 1926)
Performing gRNA-to-cell assignments for gRNA MEAF6_-_37980319.23-P1P2 (775 of 1926)
Performing gRNA-to-cell assignments for gRNA MED1_-_37607140.23-P1P2 (780 of 1926)
Performing gRNA-to-cell assignments for gRNA MED12_+_70338541.23-P1P2 (785 of 1926)
Performing gRNA-to-cell assignments for gRNA MED15_+_20862189.23-P1 (790 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_+_57479335.23-P1P2 (795 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_-_57479437.23-P1P2 (800 of 1926)
Performing gRNA-to-cell assignments for gRNA MEMO1_-_32264882.23-P1 (805 of 1926)
Performing gRNA-to-cell assignments for gRNA MIXL1_+_226411425.23-P1P2 (810 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL28_+_420487.23-P1P2 (815 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL32_+_42972018.23-P1P2 (820 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_+_179322181.23-P1P2 (825 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_-_179322213.23-P1P2 (830 of 1926)
Performing gRNA-to-cell assignments for gRNA LUZP6_-_135612750.23-P1P2 (835 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA20_-_19998023.23-P1P2 (840 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA25_+_112546538.23-P1P2 (845 of 1926)
Performing gRNA-to-cell assignments for gRNA NCBP1_+_100396082.23-P1P2 (850 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_+_172428630.23-P1P2 (855 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_-_172428841.23-P1P2 (860 of 1926)
Performing gRNA-to-cell assignments for gRNA NEDD8_-_24701363.23-P1P2 (865 of 1926)
Performing gRNA-to-cell assignments for gRNA NF2_+_29999635.23-P1P2 (870 of 1926)
Performing gRNA-to-cell assignments for gRNA NFIB_+_14314559.23-P1P2 (875 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_+_139418329.23-P2 (880 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139418349.23-P2 (885 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139440503.23-P1 (890 of 1926)
Performing gRNA-to-cell assignments for gRNA NUFIP2_+_27621073.23-P1P2 (895 of 1926)
Performing gRNA-to-cell assignments for gRNA NUP93_-_56764090.23-P1P2 (900 of 1926)
Performing gRNA-to-cell assignments for gRNA OGT_+_70752948.23-P1P2 (905 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_+_59383170.23-P1P2 (910 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_-_59383348.23-P1P2 (915 of 1926)
Performing gRNA-to-cell assignments for gRNA PAM16_-_4401236.23-P1P2 (920 of 1926)
Performing gRNA-to-cell assignments for gRNA PAXBP1_-_34143614.23-P1P2 (925 of 1926)
Performing gRNA-to-cell assignments for gRNA PBLD_-_70092436.23-P1 (930 of 1926)
Performing gRNA-to-cell assignments for gRNA PCBP2_+_53835614.23-P2 (935 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_+_94292752.23-P2 (940 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_-_94293027.23-P2 (945 of 1926)
Performing gRNA-to-cell assignments for gRNA PGRMC2_+_129208913.23-P1P2 (950 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_-_178866211.23-P1P2 (955 of 1926)
Performing gRNA-to-cell assignments for gRNA PITPNC1_-_65373433.23-P1P2 (960 of 1926)
Performing gRNA-to-cell assignments for gRNA PLCG1_+_39766176.23-P1P2 (965 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_+_151042273.23-P2 (970 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_-_151042422.23-P2 (975 of 1926)
Performing gRNA-to-cell assignments for gRNA PLRG1_-_155471154.23-P1P2 (980 of 1926)
Performing gRNA-to-cell assignments for gRNA PLXNA2_-_208417275.23-P1P2 (985 of 1926)
Performing gRNA-to-cell assignments for gRNA POLE4_-_75185725.23-P1P2 (990 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR2M_+_57998962.23-P1P2 (995 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940045.23-P1P2 (1000 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940414.23-P1P2 (1005 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3K_-_103144.23-P1P2 (1010 of 1926)
Performing gRNA-to-cell assignments for gRNA PPHLN1_-_42720007.23-P1P2 (1015 of 1926)
Performing gRNA-to-cell assignments for gRNA PPIL4_+_149867146.23-P1P2 (1020 of 1926)
Performing gRNA-to-cell assignments for gRNA PPP2R1A_+_52693317.23-P1 (1025 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_+_53195217.23-P1P2 (1030 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_-_53195332.23-P1P2 (1035 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKD1_-_30396996.23-P1P2 (1040 of 1926)
Performing gRNA-to-cell assignments for gRNA PRMT5_+_23398554.23-P1 (1045 of 1926)
Performing gRNA-to-cell assignments for gRNA C2orf53_-_27362230.23-P1P2 (1050 of 1926)
Performing gRNA-to-cell assignments for gRNA PRSS8_+_31147010.23-P1P2 (1055 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_+_36909062.23-P1P2 (1060 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_-_36909230.23-P1P2 (1065 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMD4_-_151227289.23-P1P2 (1070 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMG4_-_3259227.23-P1P2 (1075 of 1926)
Performing gRNA-to-cell assignments for gRNA PTAR1_+_72374835.23-P1P2 (1080 of 1926)
Performing gRNA-to-cell assignments for gRNA PTBP1_+_797463.23-P1P2 (1085 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_+_89623388.23-P1P2 (1090 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_-_89623526.23-P1P2 (1095 of 1926)
Performing gRNA-to-cell assignments for gRNA PTMS_-_6875591.23-P1P2 (1100 of 1926)
Performing gRNA-to-cell assignments for gRNA PYGO2_-_154934135.23-P1P2 (1105 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB18_-_27793292.23-P1P2 (1110 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB6A_+_73471983.23-P1P2 (1115 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_+_76251918.23-P1P2 (1120 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_-_76252008.23-P1P2 (1125 of 1926)
Performing gRNA-to-cell assignments for gRNA GNB2L1_-_180670840.23-P1P2 (1130 of 1926)
Performing gRNA-to-cell assignments for gRNA RAD21_-_117886635.23-P1P2 (1135 of 1926)
Performing gRNA-to-cell assignments for gRNA RANGRF_+_8192257.23-P1P2 (1140 of 1926)
Performing gRNA-to-cell assignments for gRNA RBBP5_+_205091050.23-P1P2 (1145 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM14-RBM4_+_66406280.23-P2 (1150 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM33_+_155437802.23-P1P2 (1155 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM4_+_66406319.23-P2 (1160 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_+_114271352.23-P1P2 (1165 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_-_114271700.23-P1P2 (1170 of 1926)
Performing gRNA-to-cell assignments for gRNA RBX1_-_41347852.23-P1P2 (1175 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1551_+_32112399.23-P1P2 (1180 of 1926)
Performing gRNA-to-cell assignments for gRNA RFX3_-_3525750.23-P1P2 (1185 of 1926)
Performing gRNA-to-cell assignments for gRNA RGP1_+_35749328.23-P1P2 (1190 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_+_95067102.23-P1P2 (1195 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_-_95067164.23-P1P2 (1200 of 1926)
Performing gRNA-to-cell assignments for gRNA RIBC1_-_53449606.23-P1P2 (1205 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1432_+_5629173.23-P1P2 (1210 of 1926)
Performing gRNA-to-cell assignments for gRNA RLF_+_40627522.23-P1P2 (1215 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF115_-_145611148.23-P1P2 (1220 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_+_104296142.23-P1P2 (1225 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_-_104296250.23-P1P2 (1230 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF40_-_30773681.23-P1P2 (1235 of 1926)
Performing gRNA-to-cell assignments for gRNA RPL37A_-_217363600.23-P1P2 (1240 of 1926)
Performing gRNA-to-cell assignments for gRNA RPS15A_-_18801590.23-P1P2 (1245 of 1926)
Performing gRNA-to-cell assignments for gRNA SAV1_+_51134992.23-P1P2 (1250 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_+_33104038.23-P1P2 (1255 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_-_33104343.23-P1P2 (1260 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF8_+_155054976.23-P1P2 (1265 of 1926)
Performing gRNA-to-cell assignments for gRNA EPT1_+_26569050.23-P1P2 (1270 of 1926)
Performing gRNA-to-cell assignments for gRNA SERF2_+_44069659.23-P2 (1275 of 1926)
Performing gRNA-to-cell assignments for gRNA SERPING1_+_57365137.23-P1P2 (1280 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_+_122242041.23-P1P2 (1285 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_-_122242294.23-P1P2 (1290 of 1926)
Performing gRNA-to-cell assignments for gRNA SETDB1_-_150898815.23-P1P2 (1295 of 1926)
Performing gRNA-to-cell assignments for gRNA SF1_-_64545937.23-P1P2 (1300 of 1926)
Performing gRNA-to-cell assignments for gRNA CCDC101_-_28565268.23-P1P2 (1305 of 1926)
Performing gRNA-to-cell assignments for gRNA SKA3_+_21750519.23-P1P2 (1310 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_+_78591497.23-P1P2 (1315 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_-_78591737.23-P1P2 (1320 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC29A1_-_44187438.23-P2 (1325 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC35C1_-_45826711.23-P1P2 (1330 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD1_+_146403428.23-P1P2 (1335 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD3_+_67358229.23-P1P2 (1340 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_+_48556623.23-P1P2 (1345 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_-_48556657.23-P1P2 (1350 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD5_-_135468694.23-P1P2 (1355 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA4_-_11071741.23-P1P2 (1360 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA5_+_144435243.23-P1P2 (1365 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCB1_+_24129228.23-P1 (1370 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_+_47822975.23-P1P2 (1375 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_-_47822893.23-P1P2 (1380 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCD3_-_150974212.23-P1 (1385 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCE1_-_38803669.23-P1P2 (1390 of 1926)
Performing gRNA-to-cell assignments for gRNA SMC2_-_106856631.23-P1P2 (1395 of 1926)
Performing gRNA-to-cell assignments for gRNA SMG5_+_156252537.23-P1P2 (1400 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_+_103852.23-P1P2 (1405 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_-_104242.23-P1P2 (1410 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP27_-_70121181.23-P1P2 (1415 of 1926)
Performing gRNA-to-cell assignments for gRNA SNUPN_-_75917960.23-P1P2 (1420 of 1926)
Performing gRNA-to-cell assignments for gRNA SOX4_+_21594213.23-P1 (1425 of 1926)
Performing gRNA-to-cell assignments for gRNA SPIN1_+_91003737.23-P1P2 (1430 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_+_2328680.23-P1P2 (1435 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_-_2328912.23-P1P2 (1440 of 1926)
Performing gRNA-to-cell assignments for gRNA SPTY2D1_-_18655846.23-P1P2 (1445 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP19_+_112197188.23-P1P2 (1450 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP54_-_35452244.23-P1P2 (1455 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP68_+_74068540.23-P1P2 (1460 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_+_126138730.23-P1P2 (1465 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_-_126138742.23-P1P2 (1470 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPRB_-_133524719.23-P1P2 (1475 of 1926)
Performing gRNA-to-cell assignments for gRNA SRRM2_-_2802714.23-P1P2 (1480 of 1926)
Performing gRNA-to-cell assignments for gRNA SRSF6_+_42086993.23-P1 (1485 of 1926)
Performing gRNA-to-cell assignments for gRNA SSB_+_170655450.23-P1P2 (1490 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_+_57102902.23-P1P2 (1495 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_-_57103355.23-P1P2 (1500 of 1926)
Performing gRNA-to-cell assignments for gRNA STAG2_-_123095658.23-P1 (1505 of 1926)
Performing gRNA-to-cell assignments for gRNA STK40_+_36851507.23-P1P2 (1510 of 1926)
Performing gRNA-to-cell assignments for gRNA SUPT20H_+_37633798.23-P1P2 (1515 of 1926)
Performing gRNA-to-cell assignments for gRNA SUZ12_+_30264157.23-P1P2 (1520 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_+_64901973.23-P1P2 (1525 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_-_64901981.23-P1P2 (1530 of 1926)
Performing gRNA-to-cell assignments for gRNA TADA2B_-_7045121.23-P1P2 (1535 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF10_+_6633406.23-P1P2 (1540 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF2_+_120845032.23-P1P2 (1545 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF3_-_7860540.23-P1P2 (1550 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_+_60640814.23-P1P2 (1555 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_-_60640498.23-P1P2 (1560 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF6L_-_62538954.23-P1P2 (1565 of 1926)
Performing gRNA-to-cell assignments for gRNA TAOK1_-_27717907.23-P1P2 (1570 of 1926)
Performing gRNA-to-cell assignments for gRNA FAM208A_+_56717273.23-P1P2 (1575 of 1926)
Performing gRNA-to-cell assignments for gRNA TBC1D9_+_141677481.23-P1P2 (1580 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_+_170863524.23-P1P2 (1585 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_-_170863481.23-P1P2 (1590 of 1926)
Performing gRNA-to-cell assignments for gRNA TBX3_-_115121998.23-P1P2 (1595 of 1926)
Performing gRNA-to-cell assignments for gRNA T_-_166581333.23-P1P2 (1600 of 1926)
Performing gRNA-to-cell assignments for gRNA TGFBR2_+_30648526.23-P1P2 (1605 of 1926)
Performing gRNA-to-cell assignments for gRNA THBD_+_23029196.23-P2 (1610 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_+_2977592.23-P1P2 (1615 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_-_2977643.23-P1P2 (1620 of 1926)
Performing gRNA-to-cell assignments for gRNA MESDC1_-_81293340.23-P1P2 (1625 of 1926)
Performing gRNA-to-cell assignments for gRNA TMEM41B_+_9336156.23-P1P2 (1630 of 1926)
Performing gRNA-to-cell assignments for gRNA TMTC3_-_88536127.23-P1P2 (1635 of 1926)
Performing gRNA-to-cell assignments for gRNA TOP1_+_39657538.23-P1P2 (1640 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_+_7590531.23-P1P2 (1645 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_-_7590704.23-P1P2 (1650 of 1926)
Performing gRNA-to-cell assignments for gRNA TRIM24_-_138145308.23-P1P2 (1655 of 1926)
Performing gRNA-to-cell assignments for gRNA TRRAP_-_98476465.23-P1P2 (1660 of 1926)
Performing gRNA-to-cell assignments for gRNA TSTD2_+_100395594.23-P1P2 (1665 of 1926)
Performing gRNA-to-cell assignments for gRNA UBE2E1_+_23847541.23-P1P2 (1670 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_+_38924011.23-P1P2 (1675 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_-_38924421.23-P1P2 (1680 of 1926)
Performing gRNA-to-cell assignments for gRNA UPF2_+_12084799.23-P1P2 (1685 of 1926)
Performing gRNA-to-cell assignments for gRNA UXT_-_47518505.23-P1P2 (1690 of 1926)
Performing gRNA-to-cell assignments for gRNA VHL_-_10183525.23-P1P2 (1695 of 1926)
Performing gRNA-to-cell assignments for gRNA VPS54_+_64246104.23-P1P2 (1700 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_+_28821698.23-P1P2 (1705 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_-_28822074.23-P1P2 (1710 of 1926)
Performing gRNA-to-cell assignments for gRNA WAPAL_+_88281533.23-P1P2 (1715 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR12_-_203776383.23-P1P2 (1720 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR77_+_111991860.23-P1P2 (1725 of 1926)
Performing gRNA-to-cell assignments for gRNA WIZ_-_15543016.23-P1P2 (1730 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_+_9983328.23-P1P2 (1735 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_-_9983510.23-P1P2 (1740 of 1926)
Performing gRNA-to-cell assignments for gRNA XPO5_-_43543670.23-P1P2 (1745 of 1926)
Performing gRNA-to-cell assignments for gRNA XRN1_-_142166648.23-P1P2 (1750 of 1926)
Performing gRNA-to-cell assignments for gRNA YARS2_+_32908787.23-P1P2 (1755 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB22_+_33285454.23-P1P2 (1760 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_+_6640093.23-P1P2 (1765 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_-_6640138.23-P1P2 (1770 of 1926)
Performing gRNA-to-cell assignments for gRNA ZC3H8_-_113012435.23-P1P2 (1775 of 1926)
Performing gRNA-to-cell assignments for gRNA ZNF281_+_200379151.23-P1P2 (1780 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00117 (1785 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00233 (1790 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00518 (1795 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00656 (1800 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01226 (1805 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01730 (1810 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02009 (1815 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02288 (1820 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02777 (1825 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03026 (1830 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03252 (1835 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03620 (1840 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_SGTN_04 (1845 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_SGTN_25 (1850 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_SGTN_26 (1855 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1044 (1860 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1414 (1865 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1554 (1870 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1785 (1875 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2009 (1880 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2271 (1885 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2425 (1890 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2614 (1895 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_293 (1900 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_541 (1905 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_725 (1910 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_860 (1915 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_SGTN_14 (1920 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_SGTN_01 (1925 of 1926)
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Constructing negative control pairs. ✓
Running precomputation on response ENSG00000143545.10 (1 of 226)
Running precomputation on response ENSG00000132603.15 (5 of 226)
Running precomputation on response ENSG00000155252.14 (10 of 226)
Running precomputation on response ENSG00000154175.18 (15 of 226)
Running precomputation on response ENSG00000230590.13 (20 of 226)
Running precomputation on response ENSG00000290146.1 (25 of 226)
Running precomputation on response ENSG00000035115.22 (30 of 226)
Running precomputation on response ENSG00000115207.15 (35 of 226)
Running precomputation on response ENSG00000140259.7 (40 of 226)
Running precomputation on response ENSG00000169180.12 (45 of 226)
Running precomputation on response ENSG00000110330.10 (50 of 226)
Running precomputation on response ENSG00000226361.5 (55 of 226)
Running precomputation on response ENSG00000100150.20 (60 of 226)
Running precomputation on response ENSG00000129696.13 (65 of 226)
Running precomputation on response ENSG00000133983.16 (70 of 226)
Running precomputation on response ENSG00000170949.18 (75 of 226)
Running precomputation on response ENSG00000167930.17 (80 of 226)
Running precomputation on response ENSG00000145901.16 (85 of 226)
Running precomputation on response ENSG00000163214.21 (90 of 226)
Running precomputation on response ENSG00000249803.6 (95 of 226)
Running precomputation on response ENSG00000237350.1 (100 of 226)
Running precomputation on response ENSG00000185627.19 (105 of 226)
Running precomputation on response ENSG00000290263.1 (110 of 226)
Running precomputation on response ENSG00000100441.10 (115 of 226)
Running precomputation on response ENSG00000147799.12 (120 of 226)
Running precomputation on response ENSG00000182378.15 (125 of 226)
Running precomputation on response ENSG00000108424.11 (130 of 226)
Running precomputation on response ENSG00000196850.6 (135 of 226)
Running precomputation on response ENSG00000127948.17 (140 of 226)
Running precomputation on response ENSG00000159899.16 (145 of 226)
Running precomputation on response ENSG00000172530.22 (150 of 226)
Running precomputation on response ENSG00000065029.15 (155 of 226)
Running precomputation on response ENSG00000158623.14 (160 of 226)
Running precomputation on response ENSG00000177731.16 (165 of 226)
Running precomputation on response ENSG00000127947.16 (170 of 226)
Running precomputation on response ENSG00000166037.11 (175 of 226)
Running precomputation on response ENSG00000250802.7 (180 of 226)
Running precomputation on response ENSG00000038358.15 (185 of 226)
Running precomputation on response ENSG00000214425.11 (190 of 226)
Running precomputation on response ENSG00000150275.20 (195 of 226)
Running precomputation on response ENSG00000030582.18 (200 of 226)
Running precomputation on response ENSG00000109118.14 (205 of 226)
Running precomputation on response ENSG00000187778.14 (210 of 226)
Running precomputation on response ENSG00000182919.15 (215 of 226)
Running precomputation on response ENSG00000254726.3 (220 of 226)
Running precomputation on response ENSG00000285816.1 (225 of 226)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00259&non-targeting_02159&non-targeting_02571&non-targeting_03029&non-targeting_03081 (1 of 90)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00355&non-targeting_02025&non-targeting_02089&non-targeting_03081&non-targeting_03384 (5 of 90)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_01730&non-targeting_02193&non-targeting_03180&non-targeting_03620&non-targeting_03664 (10 of 90)
Analyzing pairs containing gRNA group non-targeting_00043&non-targeting_00117&non-targeting_01509&non-targeting_02543&non-targeting_03026&non-targeting_03252 (15 of 90)
Analyzing pairs containing gRNA group non-targeting_00259&non-targeting_01732&non-targeting_01934&non-targeting_02193&non-targeting_02837&non-targeting_03081 (20 of 90)
Analyzing pairs containing gRNA group non-targeting_00259&non-targeting_01504&non-targeting_02009&non-targeting_02574&non-targeting_03081&non-targeting_03252 (25 of 90)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_02543&non-targeting_02777&non-targeting_02814&non-targeting_03252&non-targeting_03483 (30 of 90)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00918&non-targeting_01730&non-targeting_02025&non-targeting_02911&non-targeting_03384 (35 of 90)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00390&non-targeting_00518&non-targeting_00625&non-targeting_03081&non-targeting_03252 (40 of 90)
Analyzing pairs containing gRNA group non-targeting_00231&non-targeting_00614&non-targeting_00629&non-targeting_00918&non-targeting_02193&non-targeting_02543 (45 of 90)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00390&non-targeting_01732&non-targeting_02193&non-targeting_03081&non-targeting_03384 (50 of 90)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00518&non-targeting_00614&non-targeting_01226&non-targeting_02025&non-targeting_02837 (55 of 90)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_02193&non-targeting_02571&non-targeting_02911&non-targeting_03081&non-targeting_03180 (60 of 90)
Analyzing pairs containing gRNA group non-targeting_00324&non-targeting_00629&non-targeting_00972&non-targeting_02677&non-targeting_03252&non-targeting_03664 (65 of 90)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00656&non-targeting_02543&non-targeting_02574&non-targeting_03081&non-targeting_03384 (70 of 90)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_00259&non-targeting_00629&non-targeting_01601&non-targeting_01732&non-targeting_02837 (75 of 90)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_00355&non-targeting_02837&non-targeting_02911&non-targeting_03342&non-targeting_03620 (80 of 90)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00614&non-targeting_02009&non-targeting_02574&non-targeting_03026&non-targeting_03180 (85 of 90)
Analyzing pairs containing gRNA group non-targeting_00625&non-targeting_01226&non-targeting_01858&non-targeting_02009&non-targeting_02962&non-targeting_03158 (90 of 90)
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Running precomputation on response ENSG00000004487.18 (1 of 227)
Running precomputation on response ENSG00000116459.11 (5 of 227)
Running precomputation on response ENSG00000163348.4 (10 of 227)
Running precomputation on response ENSG00000138018.19 (15 of 227)
Running precomputation on response ENSG00000124380.11 (20 of 227)
Running precomputation on response ENSG00000138385.16 (25 of 227)
Running precomputation on response ENSG00000170142.12 (30 of 227)
Running precomputation on response ENSG00000114019.15 (35 of 227)
Running precomputation on response ENSG00000128052.10 (40 of 227)
Running precomputation on response ENSG00000109436.8 (45 of 227)
Running precomputation on response ENSG00000113360.17 (50 of 227)
Running precomputation on response ENSG00000120733.14 (55 of 227)
Running precomputation on response ENSG00000124766.7 (60 of 227)
Running precomputation on response ENSG00000112237.13 (65 of 227)
Running precomputation on response ENSG00000213079.10 (70 of 227)
Running precomputation on response ENSG00000136280.17 (75 of 227)
Running precomputation on response ENSG00000188186.11 (80 of 227)
Running precomputation on response ENSG00000104613.12 (85 of 227)
Running precomputation on response ENSG00000179091.5 (90 of 227)
Running precomputation on response ENSG00000106723.17 (95 of 227)
Running precomputation on response ENSG00000136908.18 (100 of 227)
Running precomputation on response ENSG00000151461.20 (105 of 227)
Running precomputation on response ENSG00000171862.14 (110 of 227)
Running precomputation on response ENSG00000179119.15 (115 of 227)
Running precomputation on response ENSG00000156603.20 (120 of 227)
Running precomputation on response ENSG00000248643.5 (125 of 227)
Running precomputation on response ENSG00000182934.12 (130 of 227)
Running precomputation on response ENSG00000134283.18 (135 of 227)
Running precomputation on response ENSG00000174106.3 (140 of 227)
Running precomputation on response ENSG00000111300.10 (145 of 227)
Running precomputation on response ENSG00000120699.14 (150 of 227)
Running precomputation on response ENSG00000129559.13 (155 of 227)
Running precomputation on response ENSG00000165417.12 (160 of 227)
Running precomputation on response ENSG00000166949.17 (165 of 227)
Running precomputation on response ENSG00000167978.18 (170 of 227)
Running precomputation on response ENSG00000103549.22 (175 of 227)
Running precomputation on response ENSG00000141510.18 (180 of 227)
Running precomputation on response ENSG00000125686.12 (185 of 227)
Running precomputation on response ENSG00000154217.16 (190 of 227)
Running precomputation on response ENSG00000005206.17 (195 of 227)
Running precomputation on response ENSG00000011451.21 (200 of 227)
Running precomputation on response ENSG00000105568.19 (205 of 227)
Running precomputation on response ENSG00000124193.16 (210 of 227)
Running precomputation on response ENSG00000128191.17 (215 of 227)
Running precomputation on response ENSG00000100387.9 (220 of 227)
Running precomputation on response ENSG00000147162.15 (225 of 227)
Analyzing pairs containing gRNA group KDM1A (1 of 227)
Analyzing pairs containing gRNA group ATP5PB (5 of 227)
Analyzing pairs containing gRNA group PYGO2 (10 of 227)
Analyzing pairs containing gRNA group SELENOI (15 of 227)
Analyzing pairs containing gRNA group SNRNP27 (20 of 227)
Analyzing pairs containing gRNA group SSB (25 of 227)
Analyzing pairs containing gRNA group UBE2E1 (30 of 227)
Analyzing pairs containing gRNA group AMOTL2 (35 of 227)
Analyzing pairs containing gRNA group KDR (40 of 227)
Analyzing pairs containing gRNA group TBC1D9 (45 of 227)
Analyzing pairs containing gRNA group DROSHA (50 of 227)
Analyzing pairs containing gRNA group KDM3B (55 of 227)
Analyzing pairs containing gRNA group SOX4 (60 of 227)
Analyzing pairs containing gRNA group CCNC (65 of 227)
Analyzing pairs containing gRNA group SCAF8 (70 of 227)
Analyzing pairs containing gRNA group CCM2 (75 of 227)
Analyzing pairs containing gRNA group LAMTOR4 (80 of 227)
Analyzing pairs containing gRNA group INTS10 (85 of 227)
Analyzing pairs containing gRNA group CYC1 (90 of 227)
Analyzing pairs containing gRNA group SPIN1 (95 of 227)
Analyzing pairs containing gRNA group DPM2 (100 of 227)
Analyzing pairs containing gRNA group UPF2 (105 of 227)
Analyzing pairs containing gRNA group PTEN (110 of 227)
Analyzing pairs containing gRNA group SPTY2D1 (115 of 227)
Analyzing pairs containing gRNA group MED19 (120 of 227)
Analyzing pairs containing gRNA group RBM14-RBM4 (125 of 227)
Analyzing pairs containing gRNA group SRPRA (130 of 227)
Analyzing pairs containing gRNA group PPHLN1 (135 of 227)
Analyzing pairs containing gRNA group LEMD3 (140 of 227)
Analyzing pairs containing gRNA group NAA25 (145 of 227)
Analyzing pairs containing gRNA group EXOSC8 (150 of 227)
Analyzing pairs containing gRNA group NEDD8 (155 of 227)
Analyzing pairs containing gRNA group GTF2A1 (160 of 227)
Analyzing pairs containing gRNA group SMAD3 (165 of 227)
Analyzing pairs containing gRNA group SRRM2 (170 of 227)
Analyzing pairs containing gRNA group RNF40 (175 of 227)
Analyzing pairs containing gRNA group TP53 (180 of 227)
Analyzing pairs containing gRNA group MED1 (185 of 227)
Analyzing pairs containing gRNA group PITPNC1 (190 of 227)
Analyzing pairs containing gRNA group SPPL2B (195 of 227)
Analyzing pairs containing gRNA group WIZ (200 of 227)
Analyzing pairs containing gRNA group PPP2R1A (205 of 227)
Analyzing pairs containing gRNA group SRSF6 (210 of 227)
Analyzing pairs containing gRNA group DGCR8 (215 of 227)
Analyzing pairs containing gRNA group RBX1 (220 of 227)
Analyzing pairs containing gRNA group OGT (225 of 227)
d2: wrote /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/TargetKD_v0.10.3_with_SE.csv

d3: running SCEPTRE from transfer inputs.
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178488.23-P1P2 (1 of 1926)
Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178620.23-P1P2 (5 of 1926)
Performing gRNA-to-cell assignments for gRNA ACTL7B_-_111617995.23-P1P2 (10 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR1B_+_52345762.23-P1P2 (15 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2A_+_148602343.23-P1P2 (20 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_+_38495339.23-P1P2 (25 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_-_38495819.23-P1P2 (30 of 1926)
Performing gRNA-to-cell assignments for gRNA ADCY2_-_7396429.23-P1P2 (35 of 1926)
Performing gRNA-to-cell assignments for gRNA AMOTL2_+_134093279.23-P1P2 (40 of 1926)
Performing gRNA-to-cell assignments for gRNA ANP32E_+_150208342.23-P1P2 (45 of 1926)
Performing gRNA-to-cell assignments for gRNA APAF1_-_99039071.23-P1P2 (50 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_+_72766725.23-P1P2 (55 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_-_72766755.23-P1P2 (60 of 1926)
Performing gRNA-to-cell assignments for gRNA ARTN_-_44399429.23-P1 (65 of 1926)
Performing gRNA-to-cell assignments for gRNA ASCC3_-_101329193.23-P1P2 (70 of 1926)
Performing gRNA-to-cell assignments for gRNA ASH2L_-_37963059.23-P1P2 (75 of 1926)
Performing gRNA-to-cell assignments for gRNA ATMIN_+_81069661.23-P1P2 (80 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_+_111992112.23-P1P2 (85 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_-_111992129.23-P1P2 (90 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP6AP1_-_153657075.23-P1P2 (95 of 1926)
Performing gRNA-to-cell assignments for gRNA ATXN2L_-_28834668.23-P1P2 (100 of 1926)
Performing gRNA-to-cell assignments for gRNA BAD_+_64052102.23-P1P2 (105 of 1926)
Performing gRNA-to-cell assignments for gRNA BAG6_+_31620398.23-P1P2 (110 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_+_49458268.23-P1P2 (115 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_-_49458323.23-P1P2 (120 of 1926)
Performing gRNA-to-cell assignments for gRNA BCAS2_-_115124045.23-P1P2 (125 of 1926)
Performing gRNA-to-cell assignments for gRNA BMPR1A_-_88516386.23-P1P2 (130 of 1926)
Performing gRNA-to-cell assignments for gRNA BRD3_-_136932911.23-P1P2 (135 of 1926)
Performing gRNA-to-cell assignments for gRNA C17orf49_+_6918114.23-P2 (140 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_+_31532436.23-P1P2 (145 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_-_31532650.23-P1P2 (150 of 1926)
Performing gRNA-to-cell assignments for gRNA CACFD1_-_136325160.23-P1P2 (155 of 1926)
Performing gRNA-to-cell assignments for gRNA CACNG1_+_65040756.23-P1P2 (160 of 1926)
Performing gRNA-to-cell assignments for gRNA CAPN5_+_76778225.23-P1P2 (165 of 1926)
Performing gRNA-to-cell assignments for gRNA CASP3_+_185570560.23-P1P2 (170 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_+_107384793.23-P1P2 (175 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_-_107384635.23-P1P2 (180 of 1926)
Performing gRNA-to-cell assignments for gRNA CCM2_-_45039941.23-P1P2 (185 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNC_-_100016267.23-P1P2 (190 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNT2_-_135676416.23-P1P2 (195 of 1926)
Performing gRNA-to-cell assignments for gRNA CCT6B_+_33288381.23-P1P2 (200 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_+_30366413.23-P1P2 (205 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_-_30366630.23-P1P2 (210 of 1926)
Performing gRNA-to-cell assignments for gRNA CD99L2_-_150067120.23-P1P2 (215 of 1926)
Performing gRNA-to-cell assignments for gRNA CDH5_-_66400630.23-P1P2 (220 of 1926)
Performing gRNA-to-cell assignments for gRNA CDKN1A_-_36646538.23-P1P2 (225 of 1926)
Performing gRNA-to-cell assignments for gRNA CENPC_+_68411190.23-P1P2 (230 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_+_30773828.23-P1P2 (235 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_-_30773391.23-P1P2 (240 of 1926)
Performing gRNA-to-cell assignments for gRNA CHMP3_-_86790451.23-P1P2 (245 of 1926)
Performing gRNA-to-cell assignments for gRNA CHORDC1_-_89955894.23-P1P2 (250 of 1926)
Performing gRNA-to-cell assignments for gRNA CNEP1R1_+_50059284.23-P1P2 (255 of 1926)
Performing gRNA-to-cell assignments for gRNA CNOT2_+_70637651.23-P1P2 (260 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_+_83956339.23-P1P2 (265 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_-_83956572.23-P1P2 (270 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS6_-_99686676.23-P1P2 (275 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_-_31634088.23-P1P2 (280 of 1926)
Performing gRNA-to-cell assignments for gRNA CTR9_+_10772881.23-P1P2 (285 of 1926)
Performing gRNA-to-cell assignments for gRNA CUL3_-_225449667.23-P1P2 (290 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_+_145150267.23-P1P2 (295 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_-_145149996.23-P1P2 (300 of 1926)
Performing gRNA-to-cell assignments for gRNA CYCS_+_25164869.23-P1P2 (305 of 1926)
Performing gRNA-to-cell assignments for gRNA DACH1_-_72441101.23-P1P2 (310 of 1926)
Performing gRNA-to-cell assignments for gRNA DBR1_+_137893673.23-P1P2 (315 of 1926)
Performing gRNA-to-cell assignments for gRNA DCTN2_-_57940910.23-P1P2 (320 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_+_61851532.23-P1P2 (325 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_-_61851812.23-P1P2 (330 of 1926)
Performing gRNA-to-cell assignments for gRNA DGCR8_-_20068158.23-P1P2 (335 of 1926)
Performing gRNA-to-cell assignments for gRNA DNAJC17_-_41099215.23-P1P2 (340 of 1926)
Performing gRNA-to-cell assignments for gRNA DNM2_+_10829157.23-P1P2 (345 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM1_+_49575018.23-P1P2 (350 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_+_130700056.23-P1P2 (355 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_-_130700076.23-P1P2 (360 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM3_+_155113020.23-P1P2 (365 of 1926)
Performing gRNA-to-cell assignments for gRNA DROSHA_-_31531982.23-P1P2 (370 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_+_3985081.23-P1 (375 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_-_3985443.23-P1 (380 of 1926)
Performing gRNA-to-cell assignments for gRNA EHMT2_-_31865199.23-P1P2 (385 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3F_+_8009347.23-P1P2 (390 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3M_-_32605405.23-P1P2 (395 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4E_+_99849831.23-P1P2 (400 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184032364.23-P1P2 (405 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184033042.23-P1P2 (410 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G2_+_10830447.23-P1 (415 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF5A_-_7211217.23-P1P2 (420 of 1926)
Performing gRNA-to-cell assignments for gRNA TCEB3_-_24069975.23-P1P2 (425 of 1926)
Performing gRNA-to-cell assignments for gRNA EMC2_+_109456080.23-P1P2 (430 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823220.23-P1P2 (435 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823252.23-P1P2 (440 of 1926)
Performing gRNA-to-cell assignments for gRNA EOMES_-_27763342.23-P1P2 (445 of 1926)
Performing gRNA-to-cell assignments for gRNA ERCC2_+_45873802.23-P1P2 (450 of 1926)
Performing gRNA-to-cell assignments for gRNA ERO1L_+_53162425.23-P1P2 (455 of 1926)
Performing gRNA-to-cell assignments for gRNA ETS1_+_128392044.23-P1P2 (460 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_+_36134592.23-P1P2 (465 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_-_36134789.23-P1P2 (470 of 1926)
Performing gRNA-to-cell assignments for gRNA EWSR1_-_29664347.23-P1P2 (475 of 1926)
Performing gRNA-to-cell assignments for gRNA EXOSC8_-_37574979.23-P1P2 (480 of 1926)
Performing gRNA-to-cell assignments for gRNA FBXW11_-_171433300.23-P1P2 (485 of 1926)
Performing gRNA-to-cell assignments for gRNA FEN1_+_61560782.23-P1P2 (490 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_+_75532080.23-P1P2 (495 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_-_75532461.23-P1P2 (500 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT7_-_139926495.23-P1P2 (505 of 1926)
Performing gRNA-to-cell assignments for gRNA G3BP1_-_151151499.23-P1P2 (510 of 1926)
Performing gRNA-to-cell assignments for gRNA GCK_+_44229021.23-P1P2 (515 of 1926)
Performing gRNA-to-cell assignments for gRNA GOLT1B_+_21654777.23-P1P2 (520 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_+_45000529.23-P1P2 (525 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_-_45000948.23-P1P2 (530 of 1926)
Performing gRNA-to-cell assignments for gRNA GPI_+_34856284.23-P1P2 (535 of 1926)
Performing gRNA-to-cell assignments for gRNA GPR107_-_132816044.23-P1P2 (540 of 1926)
Performing gRNA-to-cell assignments for gRNA GSPT1_-_12009714.23-P1P2 (545 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A1_+_81687215.23-P1P2 (550 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_+_59949355.23-P1P2 (555 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_-_59949657.23-P1P2 (560 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2H4_-_30876185.23-P1P2 (565 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF3C2_+_27579604.23-P1P2 (570 of 1926)
Performing gRNA-to-cell assignments for gRNA HAND1_+_153857732.23-P1P2 (575 of 1926)
Performing gRNA-to-cell assignments for gRNA HDAC2_-_114292253.23-P1P2 (580 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_+_81399312.23-P1P2 (585 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_-_81399351.23-P1P2 (590 of 1926)
Performing gRNA-to-cell assignments for gRNA HINT2_+_35814995.23-P1P2 (595 of 1926)
Performing gRNA-to-cell assignments for gRNA HNF4A_+_43029971.23-P1P2 (600 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPM_-_8509885.23-P1 (605 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPU_+_245027753.23-P1P2 (610 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_+_43702292.23-P1P2 (615 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_-_43702305.23-P1P2 (620 of 1926)
Performing gRNA-to-cell assignments for gRNA HSPA5_-_128003534.23-P1P2 (625 of 1926)
Performing gRNA-to-cell assignments for gRNA HUWE1_-_53710540.23-P1P2 (630 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_-_30007588.23-P1P2 (635 of 1926)
Performing gRNA-to-cell assignments for gRNA INTS10_+_19675000.23-P1P2 (640 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_+_65903046.23-P1P2 (645 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_-_65903523.23-P1P2 (650 of 1926)
Performing gRNA-to-cell assignments for gRNA C7orf26_+_6630358.23-P1P2 (655 of 1926)
Performing gRNA-to-cell assignments for gRNA IPO8_+_30848771.23-P1P2 (660 of 1926)
Performing gRNA-to-cell assignments for gRNA ITGB1BP1_+_9563543.23-P1P2 (665 of 1926)
Performing gRNA-to-cell assignments for gRNA SSFA2_+_182756743.23-P1P2 (670 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_+_109914946.23-P1P2 (675 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_-_109915112.23-P1P2 (680 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM1A_-_23346409.23-P1P2 (685 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM3B_-_137688431.23-P1P2 (690 of 1926)
Performing gRNA-to-cell assignments for gRNA KDR_-_55991659.23-P1P2 (695 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA0232_+_6784651.23-P1P2 (700 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153706791.23-P1P2 (705 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153707248.23-P1P2 (710 of 1926)
Performing gRNA-to-cell assignments for gRNA LAMTOR4_-_99746921.23-P1P2 (715 of 1926)
Performing gRNA-to-cell assignments for gRNA LARP4_-_50794864.23-P1P2 (720 of 1926)
Performing gRNA-to-cell assignments for gRNA LATS2_-_21635336.23-P1P2 (725 of 1926)
Performing gRNA-to-cell assignments for gRNA LEMD3_+_65563441.23-P1P2 (730 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_+_54960428.23-P1P2 (735 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_-_54960616.23-P1P2 (740 of 1926)
Performing gRNA-to-cell assignments for gRNA LHFPL2_-_77944519.23-P1P2 (745 of 1926)
Performing gRNA-to-cell assignments for gRNA LSM12_-_42143603.23-P1P2 (750 of 1926)
Performing gRNA-to-cell assignments for gRNA MANBAL_+_35918136.23-P1P2 (755 of 1926)
Performing gRNA-to-cell assignments for gRNA MAU2_+_19431732.23-P1P2 (760 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_+_47807678.23-P1P2 (765 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_-_47807877.23-P1P2 (770 of 1926)
Performing gRNA-to-cell assignments for gRNA MEAF6_-_37980319.23-P1P2 (775 of 1926)
Performing gRNA-to-cell assignments for gRNA MED1_-_37607140.23-P1P2 (780 of 1926)
Performing gRNA-to-cell assignments for gRNA MED12_+_70338541.23-P1P2 (785 of 1926)
Performing gRNA-to-cell assignments for gRNA MED15_+_20862189.23-P1 (790 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_+_57479335.23-P1P2 (795 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_-_57479437.23-P1P2 (800 of 1926)
Performing gRNA-to-cell assignments for gRNA MEMO1_-_32264882.23-P1 (805 of 1926)
Performing gRNA-to-cell assignments for gRNA MIXL1_+_226411425.23-P1P2 (810 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL28_+_420487.23-P1P2 (815 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL32_+_42972018.23-P1P2 (820 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_+_179322181.23-P1P2 (825 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_-_179322213.23-P1P2 (830 of 1926)
Performing gRNA-to-cell assignments for gRNA LUZP6_-_135612750.23-P1P2 (835 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA20_-_19998023.23-P1P2 (840 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA25_+_112546538.23-P1P2 (845 of 1926)
Performing gRNA-to-cell assignments for gRNA NCBP1_+_100396082.23-P1P2 (850 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_+_172428630.23-P1P2 (855 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_-_172428841.23-P1P2 (860 of 1926)
Performing gRNA-to-cell assignments for gRNA NEDD8_-_24701363.23-P1P2 (865 of 1926)
Performing gRNA-to-cell assignments for gRNA NF2_+_29999635.23-P1P2 (870 of 1926)
Performing gRNA-to-cell assignments for gRNA NFIB_+_14314559.23-P1P2 (875 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_+_139418329.23-P2 (880 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139418349.23-P2 (885 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139440503.23-P1 (890 of 1926)
Performing gRNA-to-cell assignments for gRNA NUFIP2_+_27621073.23-P1P2 (895 of 1926)
Performing gRNA-to-cell assignments for gRNA NUP93_-_56764090.23-P1P2 (900 of 1926)
Performing gRNA-to-cell assignments for gRNA OGT_+_70752948.23-P1P2 (905 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_+_59383170.23-P1P2 (910 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_-_59383348.23-P1P2 (915 of 1926)
Performing gRNA-to-cell assignments for gRNA PAM16_-_4401236.23-P1P2 (920 of 1926)
Performing gRNA-to-cell assignments for gRNA PAXBP1_-_34143614.23-P1P2 (925 of 1926)
Performing gRNA-to-cell assignments for gRNA PBLD_-_70092436.23-P1 (930 of 1926)
Performing gRNA-to-cell assignments for gRNA PCBP2_+_53835614.23-P2 (935 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_+_94292752.23-P2 (940 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_-_94293027.23-P2 (945 of 1926)
Performing gRNA-to-cell assignments for gRNA PGRMC2_+_129208913.23-P1P2 (950 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_-_178866211.23-P1P2 (955 of 1926)
Performing gRNA-to-cell assignments for gRNA PITPNC1_-_65373433.23-P1P2 (960 of 1926)
Performing gRNA-to-cell assignments for gRNA PLCG1_+_39766176.23-P1P2 (965 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_+_151042273.23-P2 (970 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_-_151042422.23-P2 (975 of 1926)
Performing gRNA-to-cell assignments for gRNA PLRG1_-_155471154.23-P1P2 (980 of 1926)
Performing gRNA-to-cell assignments for gRNA PLXNA2_-_208417275.23-P1P2 (985 of 1926)
Performing gRNA-to-cell assignments for gRNA POLE4_-_75185725.23-P1P2 (990 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR2M_+_57998962.23-P1P2 (995 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940045.23-P1P2 (1000 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940414.23-P1P2 (1005 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3K_-_103144.23-P1P2 (1010 of 1926)
Performing gRNA-to-cell assignments for gRNA PPHLN1_-_42720007.23-P1P2 (1015 of 1926)
Performing gRNA-to-cell assignments for gRNA PPIL4_+_149867146.23-P1P2 (1020 of 1926)
Performing gRNA-to-cell assignments for gRNA PPP2R1A_+_52693317.23-P1 (1025 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_+_53195217.23-P1P2 (1030 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_-_53195332.23-P1P2 (1035 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKD1_-_30396996.23-P1P2 (1040 of 1926)
Performing gRNA-to-cell assignments for gRNA PRMT5_+_23398554.23-P1 (1045 of 1926)
Performing gRNA-to-cell assignments for gRNA C2orf53_-_27362230.23-P1P2 (1050 of 1926)
Performing gRNA-to-cell assignments for gRNA PRSS8_+_31147010.23-P1P2 (1055 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_+_36909062.23-P1P2 (1060 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_-_36909230.23-P1P2 (1065 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMD4_-_151227289.23-P1P2 (1070 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMG4_-_3259227.23-P1P2 (1075 of 1926)
Performing gRNA-to-cell assignments for gRNA PTAR1_+_72374835.23-P1P2 (1080 of 1926)
Performing gRNA-to-cell assignments for gRNA PTBP1_+_797463.23-P1P2 (1085 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_+_89623388.23-P1P2 (1090 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_-_89623526.23-P1P2 (1095 of 1926)
Performing gRNA-to-cell assignments for gRNA PTMS_-_6875591.23-P1P2 (1100 of 1926)
Performing gRNA-to-cell assignments for gRNA PYGO2_-_154934135.23-P1P2 (1105 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB18_-_27793292.23-P1P2 (1110 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB6A_+_73471983.23-P1P2 (1115 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_+_76251918.23-P1P2 (1120 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_-_76252008.23-P1P2 (1125 of 1926)
Performing gRNA-to-cell assignments for gRNA GNB2L1_-_180670840.23-P1P2 (1130 of 1926)
Performing gRNA-to-cell assignments for gRNA RAD21_-_117886635.23-P1P2 (1135 of 1926)
Performing gRNA-to-cell assignments for gRNA RANGRF_+_8192257.23-P1P2 (1140 of 1926)
Performing gRNA-to-cell assignments for gRNA RBBP5_+_205091050.23-P1P2 (1145 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM14-RBM4_+_66406280.23-P2 (1150 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM33_+_155437802.23-P1P2 (1155 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM4_+_66406319.23-P2 (1160 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_+_114271352.23-P1P2 (1165 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_-_114271700.23-P1P2 (1170 of 1926)
Performing gRNA-to-cell assignments for gRNA RBX1_-_41347852.23-P1P2 (1175 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1551_+_32112399.23-P1P2 (1180 of 1926)
Performing gRNA-to-cell assignments for gRNA RFX3_-_3525750.23-P1P2 (1185 of 1926)
Performing gRNA-to-cell assignments for gRNA RGP1_+_35749328.23-P1P2 (1190 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_+_95067102.23-P1P2 (1195 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_-_95067164.23-P1P2 (1200 of 1926)
Performing gRNA-to-cell assignments for gRNA RIBC1_-_53449606.23-P1P2 (1205 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1432_+_5629173.23-P1P2 (1210 of 1926)
Performing gRNA-to-cell assignments for gRNA RLF_+_40627522.23-P1P2 (1215 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF115_-_145611148.23-P1P2 (1220 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_+_104296142.23-P1P2 (1225 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_-_104296250.23-P1P2 (1230 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF40_-_30773681.23-P1P2 (1235 of 1926)
Performing gRNA-to-cell assignments for gRNA RPL37A_-_217363600.23-P1P2 (1240 of 1926)
Performing gRNA-to-cell assignments for gRNA RPS15A_-_18801590.23-P1P2 (1245 of 1926)
Performing gRNA-to-cell assignments for gRNA SAV1_+_51134992.23-P1P2 (1250 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_+_33104038.23-P1P2 (1255 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_-_33104343.23-P1P2 (1260 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF8_+_155054976.23-P1P2 (1265 of 1926)
Performing gRNA-to-cell assignments for gRNA EPT1_+_26569050.23-P1P2 (1270 of 1926)
Performing gRNA-to-cell assignments for gRNA SERF2_+_44069659.23-P2 (1275 of 1926)
Performing gRNA-to-cell assignments for gRNA SERPING1_+_57365137.23-P1P2 (1280 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_+_122242041.23-P1P2 (1285 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_-_122242294.23-P1P2 (1290 of 1926)
Performing gRNA-to-cell assignments for gRNA SETDB1_-_150898815.23-P1P2 (1295 of 1926)
Performing gRNA-to-cell assignments for gRNA SF1_-_64545937.23-P1P2 (1300 of 1926)
Performing gRNA-to-cell assignments for gRNA CCDC101_-_28565268.23-P1P2 (1305 of 1926)
Performing gRNA-to-cell assignments for gRNA SKA3_+_21750519.23-P1P2 (1310 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_+_78591497.23-P1P2 (1315 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_-_78591737.23-P1P2 (1320 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC29A1_-_44187438.23-P2 (1325 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC35C1_-_45826711.23-P1P2 (1330 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD1_+_146403428.23-P1P2 (1335 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD3_+_67358229.23-P1P2 (1340 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_+_48556623.23-P1P2 (1345 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_-_48556657.23-P1P2 (1350 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD5_-_135468694.23-P1P2 (1355 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA4_-_11071741.23-P1P2 (1360 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA5_+_144435243.23-P1P2 (1365 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCB1_+_24129228.23-P1 (1370 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_+_47822975.23-P1P2 (1375 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_-_47822893.23-P1P2 (1380 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCD3_-_150974212.23-P1 (1385 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCE1_-_38803669.23-P1P2 (1390 of 1926)
Performing gRNA-to-cell assignments for gRNA SMC2_-_106856631.23-P1P2 (1395 of 1926)
Performing gRNA-to-cell assignments for gRNA SMG5_+_156252537.23-P1P2 (1400 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_+_103852.23-P1P2 (1405 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_-_104242.23-P1P2 (1410 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP27_-_70121181.23-P1P2 (1415 of 1926)
Performing gRNA-to-cell assignments for gRNA SNUPN_-_75917960.23-P1P2 (1420 of 1926)
Performing gRNA-to-cell assignments for gRNA SOX4_+_21594213.23-P1 (1425 of 1926)
Performing gRNA-to-cell assignments for gRNA SPIN1_+_91003737.23-P1P2 (1430 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_+_2328680.23-P1P2 (1435 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_-_2328912.23-P1P2 (1440 of 1926)
Performing gRNA-to-cell assignments for gRNA SPTY2D1_-_18655846.23-P1P2 (1445 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP19_+_112197188.23-P1P2 (1450 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP54_-_35452244.23-P1P2 (1455 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP68_+_74068540.23-P1P2 (1460 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_+_126138730.23-P1P2 (1465 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_-_126138742.23-P1P2 (1470 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPRB_-_133524719.23-P1P2 (1475 of 1926)
Performing gRNA-to-cell assignments for gRNA SRRM2_-_2802714.23-P1P2 (1480 of 1926)
Performing gRNA-to-cell assignments for gRNA SRSF6_+_42086993.23-P1 (1485 of 1926)
Performing gRNA-to-cell assignments for gRNA SSB_+_170655450.23-P1P2 (1490 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_+_57102902.23-P1P2 (1495 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_-_57103355.23-P1P2 (1500 of 1926)
Performing gRNA-to-cell assignments for gRNA STAG2_-_123095658.23-P1 (1505 of 1926)
Performing gRNA-to-cell assignments for gRNA STK40_+_36851507.23-P1P2 (1510 of 1926)
Performing gRNA-to-cell assignments for gRNA SUPT20H_+_37633798.23-P1P2 (1515 of 1926)
Performing gRNA-to-cell assignments for gRNA SUZ12_+_30264157.23-P1P2 (1520 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_+_64901973.23-P1P2 (1525 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_-_64901981.23-P1P2 (1530 of 1926)
Performing gRNA-to-cell assignments for gRNA TADA2B_-_7045121.23-P1P2 (1535 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF10_+_6633406.23-P1P2 (1540 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF2_+_120845032.23-P1P2 (1545 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF3_-_7860540.23-P1P2 (1550 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_+_60640814.23-P1P2 (1555 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_-_60640498.23-P1P2 (1560 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF6L_-_62538954.23-P1P2 (1565 of 1926)
Performing gRNA-to-cell assignments for gRNA TAOK1_-_27717907.23-P1P2 (1570 of 1926)
Performing gRNA-to-cell assignments for gRNA FAM208A_+_56717273.23-P1P2 (1575 of 1926)
Performing gRNA-to-cell assignments for gRNA TBC1D9_+_141677481.23-P1P2 (1580 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_+_170863524.23-P1P2 (1585 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_-_170863481.23-P1P2 (1590 of 1926)
Performing gRNA-to-cell assignments for gRNA TBX3_-_115121998.23-P1P2 (1595 of 1926)
Performing gRNA-to-cell assignments for gRNA T_-_166581333.23-P1P2 (1600 of 1926)
Performing gRNA-to-cell assignments for gRNA TGFBR2_+_30648526.23-P1P2 (1605 of 1926)
Performing gRNA-to-cell assignments for gRNA THBD_+_23029196.23-P2 (1610 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_+_2977592.23-P1P2 (1615 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_-_2977643.23-P1P2 (1620 of 1926)
Performing gRNA-to-cell assignments for gRNA MESDC1_-_81293340.23-P1P2 (1625 of 1926)
Performing gRNA-to-cell assignments for gRNA TMEM41B_+_9336156.23-P1P2 (1630 of 1926)
Performing gRNA-to-cell assignments for gRNA TMTC3_-_88536127.23-P1P2 (1635 of 1926)
Performing gRNA-to-cell assignments for gRNA TOP1_+_39657538.23-P1P2 (1640 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_+_7590531.23-P1P2 (1645 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_-_7590704.23-P1P2 (1650 of 1926)
Performing gRNA-to-cell assignments for gRNA TRIM24_-_138145308.23-P1P2 (1655 of 1926)
Performing gRNA-to-cell assignments for gRNA TRRAP_-_98476465.23-P1P2 (1660 of 1926)
Performing gRNA-to-cell assignments for gRNA TSTD2_+_100395594.23-P1P2 (1665 of 1926)
Performing gRNA-to-cell assignments for gRNA UBE2E1_+_23847541.23-P1P2 (1670 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_+_38924011.23-P1P2 (1675 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_-_38924421.23-P1P2 (1680 of 1926)
Performing gRNA-to-cell assignments for gRNA UPF2_+_12084799.23-P1P2 (1685 of 1926)
Performing gRNA-to-cell assignments for gRNA UXT_-_47518505.23-P1P2 (1690 of 1926)
Performing gRNA-to-cell assignments for gRNA VHL_-_10183525.23-P1P2 (1695 of 1926)
Performing gRNA-to-cell assignments for gRNA VPS54_+_64246104.23-P1P2 (1700 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_+_28821698.23-P1P2 (1705 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_-_28822074.23-P1P2 (1710 of 1926)
Performing gRNA-to-cell assignments for gRNA WAPAL_+_88281533.23-P1P2 (1715 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR12_-_203776383.23-P1P2 (1720 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR77_+_111991860.23-P1P2 (1725 of 1926)
Performing gRNA-to-cell assignments for gRNA WIZ_-_15543016.23-P1P2 (1730 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_+_9983328.23-P1P2 (1735 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_-_9983510.23-P1P2 (1740 of 1926)
Performing gRNA-to-cell assignments for gRNA XPO5_-_43543670.23-P1P2 (1745 of 1926)
Performing gRNA-to-cell assignments for gRNA XRN1_-_142166648.23-P1P2 (1750 of 1926)
Performing gRNA-to-cell assignments for gRNA YARS2_+_32908787.23-P1P2 (1755 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB22_+_33285454.23-P1P2 (1760 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_+_6640093.23-P1P2 (1765 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_-_6640138.23-P1P2 (1770 of 1926)
Performing gRNA-to-cell assignments for gRNA ZC3H8_-_113012435.23-P1P2 (1775 of 1926)
Performing gRNA-to-cell assignments for gRNA ZNF281_+_200379151.23-P1P2 (1780 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00117 (1785 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00233 (1790 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00518 (1795 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00656 (1800 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01226 (1805 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01730 (1810 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02009 (1815 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02288 (1820 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02777 (1825 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03026 (1830 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03252 (1835 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03620 (1840 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_SGTN_04 (1845 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_SGTN_25 (1850 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_SGTN_26 (1855 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1044 (1860 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1414 (1865 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1554 (1870 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1785 (1875 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2009 (1880 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2271 (1885 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2425 (1890 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2614 (1895 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_293 (1900 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_541 (1905 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_725 (1910 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_860 (1915 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_SGTN_14 (1920 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_SGTN_01 (1925 of 1926)
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Constructing negative control pairs. ✓
Running precomputation on response ENSG00000163346.17 (1 of 233)
Running precomputation on response ENSG00000253767.3 (5 of 233)
Running precomputation on response ENSG00000234287.1 (10 of 233)
Running precomputation on response ENSG00000115808.12 (15 of 233)
Running precomputation on response ENSG00000136699.20 (20 of 233)
Running precomputation on response ENSG00000117448.14 (25 of 233)
Running precomputation on response ENSG00000264324.1 (30 of 233)
Running precomputation on response ENSG00000113456.20 (35 of 233)
Running precomputation on response ENSG00000178971.17 (40 of 233)
Running precomputation on response ENSG00000177663.15 (45 of 233)
Running precomputation on response ENSG00000171204.13 (50 of 233)
Running precomputation on response ENSG00000131002.15 (55 of 233)
Running precomputation on response ENSG00000165061.15 (60 of 233)
Running precomputation on response ENSG00000181827.16 (65 of 233)
Running precomputation on response ENSG00000291064.1 (70 of 233)
Running precomputation on response ENSG00000290192.1 (75 of 233)
Running precomputation on response ENSG00000171132.14 (80 of 233)
Running precomputation on response ENSG00000067225.21 (85 of 233)
Running precomputation on response ENSG00000007129.18 (90 of 233)
Running precomputation on response ENSG00000177156.11 (95 of 233)
Running precomputation on response ENSG00000156482.11 (100 of 233)
Running precomputation on response ENSG00000258526.8 (105 of 233)
Running precomputation on response ENSG00000230409.3 (110 of 233)
Running precomputation on response ENSG00000047410.14 (115 of 233)
Running precomputation on response ENSG00000196839.14 (120 of 233)
Running precomputation on response ENSG00000106624.11 (125 of 233)
Running precomputation on response ENSG00000163382.13 (130 of 233)
Running precomputation on response ENSG00000126821.8 (135 of 233)
Running precomputation on response ENSG00000164823.11 (140 of 233)
Running precomputation on response ENSG00000273748.1 (145 of 233)
Running precomputation on response ENSG00000181264.9 (150 of 233)
Running precomputation on response ENSG00000128487.19 (155 of 233)
Running precomputation on response ENSG00000172183.16 (160 of 233)
Running precomputation on response ENSG00000186260.17 (165 of 233)
Running precomputation on response ENSG00000137942.18 (170 of 233)
Running precomputation on response ENSG00000175155.11 (175 of 233)
Running precomputation on response ENSG00000136997.22 (180 of 233)
Running precomputation on response ENSG00000137992.15 (185 of 233)
Running precomputation on response ENSG00000189120.5 (190 of 233)
Running precomputation on response ENSG00000255245.4 (195 of 233)
Running precomputation on response ENSG00000086289.13 (200 of 233)
Running precomputation on response ENSG00000182359.15 (205 of 233)
Running precomputation on response ENSG00000189058.9 (210 of 233)
Running precomputation on response ENSG00000186834.4 (215 of 233)
Running precomputation on response ENSG00000116138.13 (220 of 233)
Running precomputation on response ENSG00000101230.6 (225 of 233)
Running precomputation on response ENSG00000171724.3 (230 of 233)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00259&non-targeting_02159&non-targeting_02571&non-targeting_03029&non-targeting_03081 (1 of 91)
Analyzing pairs containing gRNA group non-targeting_01226&non-targeting_01509&non-targeting_01732&non-targeting_02543&non-targeting_02837&non-targeting_03026 (5 of 91)
Analyzing pairs containing gRNA group non-targeting_00518&non-targeting_03081&non-targeting_03342&non-targeting_03427&non-targeting_03635&non-targeting_03664 (10 of 91)
Analyzing pairs containing gRNA group non-targeting_01509&non-targeting_01812&non-targeting_02025&non-targeting_02677&non-targeting_02777&non-targeting_03384 (15 of 91)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_00324&non-targeting_00390&non-targeting_00629&non-targeting_01934&non-targeting_03180 (20 of 91)
Analyzing pairs containing gRNA group non-targeting_00518&non-targeting_00656&non-targeting_01730&non-targeting_02574&non-targeting_03029&non-targeting_03158 (25 of 91)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_02543&non-targeting_02777&non-targeting_02814&non-targeting_03252&non-targeting_03483 (30 of 91)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00918&non-targeting_01730&non-targeting_02025&non-targeting_02911&non-targeting_03384 (35 of 91)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00390&non-targeting_01601&non-targeting_02009&non-targeting_02288&non-targeting_03029 (40 of 91)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_02159&non-targeting_02193&non-targeting_02677&non-targeting_02962&non-targeting_03029 (45 of 91)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00390&non-targeting_01732&non-targeting_02193&non-targeting_03081&non-targeting_03384 (50 of 91)
Analyzing pairs containing gRNA group non-targeting_00629&non-targeting_01185&non-targeting_01732&non-targeting_02089&non-targeting_02677&non-targeting_03483 (55 of 91)
Analyzing pairs containing gRNA group non-targeting_01185&non-targeting_02543&non-targeting_02677&non-targeting_03029&non-targeting_03342&non-targeting_03664 (60 of 91)
Analyzing pairs containing gRNA group non-targeting_00614&non-targeting_01185&non-targeting_01226&non-targeting_02089&non-targeting_02159&non-targeting_02962 (65 of 91)
Analyzing pairs containing gRNA group non-targeting_00355&non-targeting_00656&non-targeting_01730&non-targeting_02574&non-targeting_02911&non-targeting_03664 (70 of 91)
Analyzing pairs containing gRNA group non-targeting_00355&non-targeting_00518&non-targeting_00555&non-targeting_00794&non-targeting_01185&non-targeting_02159 (75 of 91)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00196&non-targeting_00213&non-targeting_00518&non-targeting_02574&non-targeting_02837 (80 of 91)
Analyzing pairs containing gRNA group non-targeting_00614&non-targeting_01599&non-targeting_02009&non-targeting_02543&non-targeting_02837&non-targeting_03081 (85 of 91)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00555&non-targeting_01226&non-targeting_01812&non-targeting_02911&non-targeting_03180 (90 of 91)
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Running precomputation on response ENSG00000004487.18 (1 of 238)
Running precomputation on response ENSG00000117000.9 (5 of 238)
Running precomputation on response ENSG00000265491.5 (10 of 238)
Running precomputation on response ENSG00000198952.8 (15 of 238)
Running precomputation on response ENSG00000119185.13 (20 of 238)
Running precomputation on response ENSG00000124380.11 (25 of 238)
Running precomputation on response ENSG00000138434.17 (30 of 238)
Running precomputation on response ENSG00000163513.19 (35 of 238)
Running precomputation on response ENSG00000138231.14 (40 of 238)
Running precomputation on response ENSG00000170871.12 (45 of 238)
Running precomputation on response ENSG00000164040.18 (50 of 238)
Running precomputation on response ENSG00000164305.20 (55 of 238)
Running precomputation on response ENSG00000164292.13 (60 of 238)
Running precomputation on response ENSG00000113196.3 (65 of 238)
Running precomputation on response ENSG00000204435.15 (70 of 238)
Running precomputation on response ENSG00000112237.13 (75 of 238)
Running precomputation on response ENSG00000213079.10 (80 of 238)
Running precomputation on response ENSG00000196367.15 (85 of 238)
Running precomputation on response ENSG00000082014.17 (90 of 238)
Running precomputation on response ENSG00000064313.13 (95 of 238)
Running precomputation on response ENSG00000107185.10 (100 of 238)
Running precomputation on response ENSG00000155827.12 (105 of 238)
Running precomputation on response ENSG00000169925.17 (110 of 238)
Running precomputation on response ENSG00000095787.26 (115 of 238)
Running precomputation on response ENSG00000175390.15 (120 of 238)
Running precomputation on response ENSG00000149100.13 (125 of 238)
Running precomputation on response ENSG00000110048.12 (130 of 238)
Running precomputation on response ENSG00000173933.21 (135 of 238)
Running precomputation on response ENSG00000134954.14 (140 of 238)
Running precomputation on response ENSG00000134283.18 (145 of 238)
Running precomputation on response ENSG00000174106.3 (150 of 238)
Running precomputation on response ENSG00000111300.10 (155 of 238)
Running precomputation on response ENSG00000102710.21 (160 of 238)
Running precomputation on response ENSG00000184304.17 (165 of 238)
Running precomputation on response ENSG00000140264.21 (170 of 238)
Running precomputation on response ENSG00000140395.9 (175 of 238)
Running precomputation on response ENSG00000134419.15 (180 of 238)
Running precomputation on response ENSG00000102900.13 (185 of 238)
Running precomputation on response ENSG00000141510.18 (190 of 238)
Running precomputation on response ENSG00000125686.12 (195 of 238)
Running precomputation on response ENSG00000154217.16 (200 of 238)
Running precomputation on response ENSG00000005206.17 (205 of 238)
Running precomputation on response ENSG00000011451.21 (210 of 238)
Running precomputation on response ENSG00000167615.17 (215 of 238)
Running precomputation on response ENSG00000124193.16 (220 of 238)
Running precomputation on response ENSG00000128191.17 (225 of 238)
Running precomputation on response ENSG00000100387.9 (230 of 238)
Running precomputation on response ENSG00000147162.15 (235 of 238)
Analyzing pairs containing gRNA group KDM1A (1 of 238)
Analyzing pairs containing gRNA group RLF (5 of 238)
Analyzing pairs containing gRNA group RNF115 (10 of 238)
Analyzing pairs containing gRNA group SMG5 (15 of 238)
Analyzing pairs containing gRNA group ITGB1BP1 (20 of 238)
Analyzing pairs containing gRNA group SNRNP27 (25 of 238)
Analyzing pairs containing gRNA group ITPRID2 (30 of 238)
Analyzing pairs containing gRNA group TGFBR2 (35 of 238)
Analyzing pairs containing gRNA group DBR1 (40 of 238)
Analyzing pairs containing gRNA group KIAA0232 (45 of 238)
Analyzing pairs containing gRNA group PGRMC2 (50 of 238)
Analyzing pairs containing gRNA group CASP3 (55 of 238)
Analyzing pairs containing gRNA group RHOBTB3 (60 of 238)
Analyzing pairs containing gRNA group HAND1 (65 of 238)
Analyzing pairs containing gRNA group CSNK2B (70 of 238)
Analyzing pairs containing gRNA group CCNC (75 of 238)
Analyzing pairs containing gRNA group SCAF8 (80 of 238)
Analyzing pairs containing gRNA group TRRAP (85 of 238)
Analyzing pairs containing gRNA group SMARCD3 (90 of 238)
Analyzing pairs containing gRNA group TAF2 (95 of 238)
Analyzing pairs containing gRNA group RGP1 (100 of 238)
Analyzing pairs containing gRNA group RNF20 (105 of 238)
Analyzing pairs containing gRNA group BRD3 (110 of 238)
Analyzing pairs containing gRNA group WAC (115 of 238)
Analyzing pairs containing gRNA group EIF3F (120 of 238)
Analyzing pairs containing gRNA group EIF3M (125 of 238)
Analyzing pairs containing gRNA group OSBP (130 of 238)
Analyzing pairs containing gRNA group RBM4 (135 of 238)
Analyzing pairs containing gRNA group ETS1 (140 of 238)
Analyzing pairs containing gRNA group PPHLN1 (145 of 238)
Analyzing pairs containing gRNA group LEMD3 (150 of 238)
Analyzing pairs containing gRNA group NAA25 (155 of 238)
Analyzing pairs containing gRNA group SUPT20H (160 of 238)
Analyzing pairs containing gRNA group PRKD1 (165 of 238)
Analyzing pairs containing gRNA group SERF2 (170 of 238)
Analyzing pairs containing gRNA group SKIC8 (175 of 238)
Analyzing pairs containing gRNA group RPS15A (180 of 238)
Analyzing pairs containing gRNA group NUP93 (185 of 238)
Analyzing pairs containing gRNA group TP53 (190 of 238)
Analyzing pairs containing gRNA group MED1 (195 of 238)
Analyzing pairs containing gRNA group PITPNC1 (200 of 238)
Analyzing pairs containing gRNA group SPPL2B (205 of 238)
Analyzing pairs containing gRNA group WIZ (210 of 238)
Analyzing pairs containing gRNA group LENG8 (215 of 238)
Analyzing pairs containing gRNA group SRSF6 (220 of 238)
Analyzing pairs containing gRNA group DGCR8 (225 of 238)
Analyzing pairs containing gRNA group RBX1 (230 of 238)
Analyzing pairs containing gRNA group OGT (235 of 238)
d3: wrote /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d3_sceptre_transfer/TargetKD_v0.10.3_with_SE.csv
There were 50 or more warnings (use warnings() to see the first 50)
SCEPTRE subprocess completed successfully.
R script: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/sceptre_on_target_KD_D0_D3_overlay_plotnine_sceptre_env_v2/run_sceptre_with_official_se.R
R log: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/sceptre_on_target_KD_D0_D3_overlay_plotnine_sceptre_env_v2/run_sceptre_with_official_se.log

```

## Markdown cell 9

## 4. Construct SCEPTRE-aligned estimates and official-SE confidence intervals

For each day and target gene, let

\[
\widehat F
\]

be SCEPTRE v0.10.3's `fold_change` estimate, and let

\[
s_F
\]

be its official `se_fold_change` output. The displayed control reference is fixed at 1, and the filled point is \(\widehat F\).

To keep the interval positive, the notebook forms an approximate Wald interval on the log2 scale. By the delta method,

\[
\widehat\theta = \log_2(\widehat F),
\qquad
s_\theta \approx \frac{s_F}{\widehat F\ln 2}.
\]

For confidence level \(1-\alpha\), the plotted relative-expression interval is

\[
\left[
2^{\widehat\theta-z_{1-\alpha/2}s_\theta},
\;2^{\widehat\theta+z_{1-\alpha/2}s_\theta}
\right].
\]

This is a model-based interval for the fold-change estimator. The orange/black significance call is based separately on the conditional-resampling SCEPTRE p-value followed by BH correction across **all QC-passing on-target pairs within the same day**, before selecting the four displayed genes. The notebook does not infer an SE from the p-value or from `z_orig`.

## Code cell 10

```python


def parse_bool_series(series: pd.Series, name: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    text = series.astype(str).str.strip().str.lower()
    mapping = {"true": True, "t": True, "1": True, "false": False, "f": False, "0": False}
    parsed = text.map(mapping)
    if parsed.isna().any():
        examples = series[parsed.isna()].astype(str).drop_duplicates().head().tolist()
        raise ValueError(f"Cannot parse {name} as Boolean; examples: {examples}")
    return parsed.astype(bool)

def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p values, preserving original order."""
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or len(p) == 0 or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError("BH input must be a nonempty finite vector in [0, 1].")
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    n = len(p)
    adjusted_ranked = ranked * n / np.arange(1, n + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.clip(adjusted_ranked, 0, 1)
    adjusted = np.empty(n, dtype=float)
    adjusted[order] = adjusted_ranked
    return adjusted


required_columns = {
    "response_id", "grna_target", "n_nonzero_trt", "n_nonzero_cntrl",
    "pass_qc", "p_value", "fold_change", "se_fold_change",
    "log_2_fold_change", "q_bh", "q_significant", "day",
}
all_results = []
q_audits = []
for day in cfg.DAYS:
    path = transfer_dir(day) / cfg.RESULT_FILENAME
    if not path.is_file():
        raise FileNotFoundError(path)
    result = pd.read_csv(path)
    missing = required_columns - set(result.columns)
    if missing:
        raise ValueError(f"{path} lacks columns: {sorted(missing)}")
    if set(result["day"].astype(str)) != {day}:
        raise ValueError(f"{path} has an unexpected day label.")
    for column in ["p_value", "fold_change", "se_fold_change", "log_2_fold_change", "q_bh"]:
        result[column] = pd.to_numeric(result[column], errors="coerce")
    result["pass_qc"] = parse_bool_series(result["pass_qc"], "pass_qc")
    result["q_significant"] = parse_bool_series(result["q_significant"], "q_significant")

    eligible = result["pass_qc"] & np.isfinite(result["p_value"])
    q_check = bh_adjust(result.loc[eligible, "p_value"].to_numpy(float))
    max_q_difference = float(np.max(np.abs(q_check - result.loc[eligible, "q_bh"].to_numpy(float))))
    if max_q_difference > 1e-12:
        raise ValueError(f"BH audit failed for {day}: max |q_python - q_R| = {max_q_difference}")
    q_audits.append({"day": day, "n_family": int(eligible.sum()), "max_q_difference": max_q_difference})
    all_results.append(result)

all_results = pd.concat(all_results, ignore_index=True)
pd.DataFrame(q_audits).to_csv(OUTDIR / "q_value_audit.csv", index=False)
all_results.to_csv(OUTDIR / "all_on_target_sceptre_results_v0.10.3_with_SE.csv", index=False)

selected = all_results[all_results["grna_target"].isin(cfg.TARGETS)].copy()
expected_keys = {(day, gene) for day in cfg.DAYS for gene in cfg.TARGETS}
actual_keys = set(zip(selected["day"], selected["grna_target"]))
if actual_keys != expected_keys or len(selected) != len(expected_keys):
    missing = sorted(expected_keys - actual_keys)
    duplicates = selected.duplicated(["day", "grna_target"], keep=False)
    raise ValueError(
        f"Expected exactly one row per day/gene. Missing={missing}; "
        f"duplicated rows={selected.loc[duplicates, ['day', 'grna_target']].to_dict('records')}"
    )
if not selected["pass_qc"].all():
    failed = selected.loc[~selected["pass_qc"], ["day", "grna_target"]].to_dict("records")
    raise ValueError(f"Displayed pairs failed SCEPTRE pairwise QC: {failed}")

numeric = selected[["fold_change", "se_fold_change", "log_2_fold_change", "q_bh"]].to_numpy(float)
if not np.isfinite(numeric).all():
    raise ValueError("Selected pairs have missing or nonfinite effect/SE/q values.")
if np.any(selected["fold_change"].to_numpy(float) <= 0):
    raise ValueError("SCEPTRE fold changes must be positive.")
if np.any(selected["se_fold_change"].to_numpy(float) < 0):
    raise ValueError("SCEPTRE fold-change standard errors must be nonnegative.")

selected["log2_fc_recomputed"] = np.log2(selected["fold_change"].to_numpy(float))
selected["log2_fc_discrepancy"] = (
    selected["log2_fc_recomputed"] - selected["log_2_fold_change"].to_numpy(float)
).abs()
if selected["log2_fc_discrepancy"].max() > 1e-6:
    raise ValueError("fold_change and log_2_fold_change are inconsistent for displayed pairs.")

z_value = float(norm.ppf(1 - (1 - cfg.CI_LEVEL) / 2))
selected["se_log2_fold_change_delta"] = (
    selected["se_fold_change"] /
    (selected["fold_change"] * np.log(2.0))
)
selected["ci_low"] = 2.0 ** (
    selected["log_2_fold_change"] - z_value * selected["se_log2_fold_change_delta"]
)
selected["ci_high"] = 2.0 ** (
    selected["log_2_fold_change"] + z_value * selected["se_log2_fold_change_delta"]
)
selected["ci_method"] = "SCEPTRE_v0.10.3_se_fold_change_delta_method_log2_Wald"
selected["q_significant_figure"] = selected["q_bh"] <= cfg.Q_THRESHOLD
selected["knockdown_pct"] = 100.0 * (1.0 - selected["fold_change"])
selected["day_label"] = selected["day"].map(cfg.DAY_LABELS)
selected["gene_order"] = selected["grna_target"].map({gene: i for i, gene in enumerate(cfg.TARGETS)})
selected["day_order"] = selected["day"].map({day: i for i, day in enumerate(cfg.DAYS)})
selected = selected.sort_values(["gene_order", "day_order"], kind="stable").reset_index(drop=True)
selected.to_csv(OUTDIR / "sceptre_4genes_effects_and_CIs.csv", index=False)

display_columns = [
    "day", "grna_target", "fold_change", "se_fold_change",
    "ci_low", "ci_high", "p_value", "q_bh", "knockdown_pct", "ci_method",
]
display(selected[display_columns])
```

### Text output (text/plain)

```text
   day grna_target  fold_change  se_fold_change    ci_low   ci_high  \
0   d0         NF2     0.417021        0.076402  0.291212  0.597181   
1   d1         NF2     0.445419        0.119781  0.262946  0.754520   
2   d2         NF2     0.371508        0.097529  0.222080  0.621479   
3   d3         NF2     0.237905        0.067893  0.135985  0.416214   
4   d0       LATS2     0.282267        0.058146  0.188501  0.422674   
5   d1       LATS2     0.484775        0.092219  0.333900  0.703824   
6   d2       LATS2     0.253111        0.068305  0.149144  0.429554   
7   d3       LATS2     0.273726        0.045401  0.197758  0.378877   
8   d0      AMOTL2     0.499626        0.082285  0.361791  0.689972   
9   d1      AMOTL2     0.552278        0.048238  0.465382  0.655399   
10  d2      AMOTL2     0.446415        0.081121  0.312652  0.637405   
11  d3      AMOTL2     0.372636        0.055916  0.277689  0.500047   
12  d0       TAOK1     0.572794        0.042520  0.495236  0.662499   
13  d1       TAOK1     0.467619        0.050442  0.378507  0.577711   
14  d2       TAOK1     0.460382        0.049497  0.372909  0.568373   
15  d3       TAOK1     0.404188        0.038502  0.335352  0.487155   

         p_value          q_bh  knockdown_pct  \
0   1.185497e-08  2.417303e-08      58.297948   
1   1.356444e-04  2.242285e-04      55.458123   
2   2.168602e-05  4.922727e-05      62.849239   
3   3.867726e-09  1.394726e-08      76.209515   
4   4.658434e-13  1.336100e-12      71.773336   
5   1.891607e-05  3.619373e-05      51.522471   
6   2.315815e-09  1.095188e-08      74.688873   
7   1.906900e-20  1.891009e-19      72.627408   
8   2.917218e-06  4.426709e-06      50.037440   
9   7.786174e-14  6.306801e-13      44.772202   
10  1.104052e-06  3.056340e-06      55.358550   
11  2.996841e-16  2.377494e-15      62.736419   
12  3.156089e-16  1.328611e-15      42.720563   
13  3.300347e-17  5.346562e-16      53.238065   
14  1.386274e-15  1.748246e-14      53.961786   
15  3.475918e-30  6.363604e-29      59.581180   

                                            ci_method  
0   SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
1   SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
2   SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
3   SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
4   SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
5   SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
6   SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
7   SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
8   SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
9   SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
10  SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
11  SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
12  SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
13  SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
14  SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
15  SCEPTRE_v0.10.3_se_fold_change_delta_method_lo...  
```

## Markdown cell 11

## 5. Make one overlaid plotnine panel

Each gene occupies one x-axis group. The four days are horizontally offset within the group and distinguished by shape. For each day, the open symbol marks the control reference at 1 and the filled symbol marks the SCEPTRE fold-change estimate. Confidence bars are drawn only for the estimated fold change because the control value is the reference denominator, not a separately estimated SCEPTRE point in this display.


## Code cell 12

```python
plot_rows = []
gene_positions = {gene: index for index, gene in enumerate(cfg.TARGETS)}
for row in selected.itertuples(index=False):
    x_position = gene_positions[row.grna_target] + cfg.DAY_OFFSETS[row.day]
    common = {
        "target": row.grna_target,
        "day": row.day,
        "day_label": row.day_label,
        "x_plot": x_position,
        "p_value": row.p_value,
        "q_bh": row.q_bh,
        "q_significant": bool(row.q_significant_figure),
        "sceptre_fold_change": row.fold_change,
        "sceptre_se_fold_change": row.se_fold_change,
        "sceptre_log_2_fold_change": row.log_2_fold_change,
        "knockdown_pct": row.knockdown_pct,
        "ci_method": row.ci_method,
    }
    plot_rows.append({
        **common,
        "group": "control",
        "estimate": 1.0,
        "ci_low": np.nan,
        "ci_high": np.nan,
    })
    plot_rows.append({
        **common,
        "group": "targeting",
        "estimate": row.fold_change,
        "ci_low": row.ci_low,
        "ci_high": row.ci_high,
    })

plot_data = pd.DataFrame(plot_rows)
significant_label = f"q ≤ {cfg.Q_THRESHOLD:g}"
legend_order = ["control", "NS", significant_label]
plot_data["point_class"] = np.where(
    plot_data["group"].eq("control"),
    "control",
    np.where(plot_data["q_significant"], significant_label, "NS"),
)
plot_data["point_class"] = pd.Categorical(
    plot_data["point_class"], categories=legend_order, ordered=True
)
day_labels = [cfg.DAY_LABELS[day] for day in cfg.DAYS]
plot_data["day_label"] = pd.Categorical(
    plot_data["day_label"], categories=day_labels, ordered=True
)
plot_data = plot_data.sort_values(
    ["target", "day_label", "group"], kind="stable"
).reset_index(drop=True)
plot_data.to_csv(OUTDIR / "plot_source_data_control_and_targeting.csv", index=False)

targeting_data = plot_data[plot_data["group"].eq("targeting")].copy()
fill_values = {
    "control": "white",
    "NS": cfg.NS_COLOR,
    significant_label: cfg.SIGNIFICANT_COLOR,
}
edge_values = {
    "control": cfg.CONTROL_EDGE,
    "NS": cfg.NS_COLOR,
    significant_label: cfg.SIGNIFICANT_COLOR,
}

upper = float(max(1.0, targeting_data["ci_high"].max()))
y_max = max(1.08, upper * 1.08)
title = cfg.TITLE if not cfg.PANEL_LETTER else f"{cfg.PANEL_LETTER}  {cfg.TITLE}"

panel = (
    ggplot(plot_data, aes(x="x_plot", y="estimate"))
    + geom_hline(
        yintercept=1.0, linetype="dashed",
        color=cfg.REFERENCE_LINE_COLOR, size=0.5,
    )
    + geom_point(
        aes(shape="day_label", fill="point_class", color="point_class"),
        size=cfg.POINT_SIZE,
        stroke=cfg.POINT_STROKE,
        show_legend=True,
    )
    # Draw literal CI endpoints after the points so narrow intervals remain visible.
    + geom_errorbar(
        data=targeting_data,
        mapping=aes(x="x_plot", ymin="ci_low", ymax="ci_high"),
        inherit_aes=False,
        width=cfg.CI_CAP_WIDTH,
        color=cfg.CI_COLOR,
        size=cfg.CI_LINEWIDTH,
        show_legend=False,
    )
    + scale_shape_manual(values=cfg.DAY_SHAPES, name="Time point")
    + scale_fill_manual(
        values=fill_values,
        limits=legend_order,
        breaks=legend_order,
        drop=False,
        name="",
    )
    + scale_color_manual(
        values=edge_values,
        limits=legend_order,
        breaks=legend_order,
        drop=False,
        name="",
        guide=None,
    )
    + guides(
        fill=guide_legend(
            title="", order=1, nrow=1, byrow=True,
            override_aes={
                "shape": "o", "size": cfg.POINT_SIZE,
                "stroke": cfg.POINT_STROKE,
            },
        ),
        shape=guide_legend(
            title="Time point", order=2, nrow=1, byrow=True,
            override_aes={
                "fill": "white", "color": "#444444",
                "size": cfg.POINT_SIZE, "stroke": cfg.POINT_STROKE,
            },
        ),
    )
    + scale_x_continuous(
        breaks=list(range(len(cfg.TARGETS))),
        labels=cfg.TARGETS,
        expand=(0.08, 0.08),
    )
    + scale_y_continuous(expand=(0, 0))
    + coord_cartesian(
        xlim=(-0.55, len(cfg.TARGETS) - 0.45),
        ylim=(0, y_max),
        expand=False,
    )
    + labs(title=title, x=cfg.X_LABEL, y=cfg.Y_LABEL)
    + theme_classic(base_size=cfg.FONT_SIZE, base_family=AI_FONT)
    + theme(
        figure_size=cfg.FIGSIZE,
        dpi=120,
        svg_usefonts=True,
        panel_grid_major=element_blank(),
        panel_grid_minor=element_blank(),
        axis_text_x=element_text(rotation=35, ha="right", style="italic"),
        plot_title=element_text(weight="bold", ha="left"),
        legend_position="top",
        legend_direction="horizontal",
        legend_box="horizontal",
        legend_title=element_text(size=cfg.FONT_SIZE - 1),
        legend_text=element_text(size=cfg.FONT_SIZE - 1),
        axis_line=element_line(color="black", size=0.7),
        axis_ticks=element_line(color="black", size=0.5),
    )
)

figure_paths = []
for extension in ["pdf", "svg", "png"]:
    path = OUTDIR / f"sceptre_on_target_KD_D0_D3_overlay.{extension}"
    panel.save(
        filename=str(path),
        width=cfg.FIGSIZE[0],
        height=cfg.FIGSIZE[1],
        units="in",
        dpi=cfg.EXPORT_DPI,
        verbose=False,
    )
    if not path.is_file() or path.stat().st_size == 0:
        raise IOError(f"Figure export failed: {path}")
    figure_paths.append(path)

svg_text = (OUTDIR / "sceptre_on_target_KD_D0_D3_overlay.svg").read_text(
    encoding="utf-8", errors="ignore"
)
if "<text" not in svg_text:
    raise IOError("SVG export does not contain editable text elements.")

display(panel)
print("Saved figure files:")
for path in figure_paths:
    print(path)
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f8ef3b80310>
```

### Text output (stream)

```text
Saved figure files:
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/sceptre_on_target_KD_D0_D3_overlay_plotnine_sceptre_env_v2/sceptre_on_target_KD_D0_D3_overlay.pdf
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/sceptre_on_target_KD_D0_D3_overlay_plotnine_sceptre_env_v2/sceptre_on_target_KD_D0_D3_overlay.svg
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/sceptre_on_target_KD_D0_D3_overlay_plotnine_sceptre_env_v2/sceptre_on_target_KD_D0_D3_overlay.png

```

## Markdown cell 13

## 6. Save figure legend, provenance and a recoverable output bundle


## Code cell 14

```python
control_description = "the control reference, fixed at a relative-expression value of 1"
legend_text = (
    "On-target CRISPRi knockdown of NF2, LATS2, AMOTL2 and TAOK1 across D0–D3. "
    "Time points are overlaid within each target-gene group and distinguished by symbol shape. "
    f"Open symbols show {control_description}; filled symbols show SCEPTRE v0.10.3 fold-change estimates for cells assigned guides targeting the indicated gene. "
    f"Grey bars show approximate {cfg.CI_LEVEL * 100:g}% Wald confidence intervals formed on the log2 fold-change scale by applying the delta method to SCEPTRE's official fold_change and se_fold_change outputs. "
    "These model-based intervals are not obtained by inverting SCEPTRE's conditional-resampling test. "
    f"Orange symbols indicate BH-adjusted SCEPTRE q ≤ {cfg.Q_THRESHOLD:g}; black symbols indicate q > {cfg.Q_THRESHOLD:g}. "
    "BH adjustment was performed separately within each time point across all QC-passing on-target discovery pairs before selecting the four displayed genes. "
    "The dashed horizontal line marks the control reference of 1."
)
(OUTDIR / "figure_legend.txt").write_text(legend_text + "\n", encoding="utf-8")


def file_record(path: Path) -> dict:
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }

result_files = {day: transfer_dir(day) / cfg.RESULT_FILENAME for day in cfg.DAYS}
provenance = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "analysis_source": "ipsc_ec_sceptre.ipynb",
    "environment_strategy": "sc-dl Jupyter kernel plus direct sceptre_env/bin/Rscript subprocess",
    "sc_dl_env": str(cfg.SC_DL_ENV),
    "sc_dl_python": sys.executable,
    "sceptre_env": str(cfg.SCEPTRE_ENV),
    "sceptre_rscript": str(cfg.SCEPTRE_RSCRIPT),
    "sceptre_private_r_library": str(cfg.SCEPTRE_PRIVATE_R_LIB),
    "sceptre_required_version": cfg.SCEPTRE_REQUIRED_VERSION,
    "sceptre_github_ref": cfg.SCEPTRE_GITHUB_REF,
    "install_command": INSTALL_COMMAND_USED,
    "analysis_r_command": R_COMMAND_USED,
    "install_script": file_record(INSTALL_SCRIPT_PATH),
    "install_log": file_record(INSTALL_LOG_PATH),
    "sceptre_preflight": file_record(OUTDIR / "sceptre_private_library_preflight.txt"),
    "analysis_r_script": file_record(R_SCRIPT_PATH),
    "analysis_r_log": file_record(R_LOG_PATH),
    "sceptre_result_files": {day: file_record(path) for day, path in result_files.items()},
    "targets": cfg.TARGETS,
    "days": cfg.DAYS,
    "q_threshold": cfg.Q_THRESHOLD,
    "q_family": "all QC-passing on-target discovery pairs separately within each day",
    "ci_level": cfg.CI_LEVEL,
    "ci_method": (
        "delta-method Wald interval on log2(fold_change), using official SCEPTRE "
        "v0.10.3 se_log2 = se_fold_change / (fold_change * ln(2))"
    ),
    "control_reference": 1.0,
    "plotting_library": "plotnine",
    "plotnine_version": p9.__version__,
    "python_packages": {
        name: importlib.metadata.version(name)
        for name in ["numpy", "pandas", "scipy", "muon", "matplotlib", "plotnine", "mizani"]
    },
    "font": AI_FONT,
    "figure_files": [file_record(path) for path in figure_paths],
    "original_combined_result_present": cfg.OLD_COMBINED_CSV.is_file(),
}
(OUTDIR / "analysis_provenance.json").write_text(
    json.dumps(provenance, indent=2, default=str) + "\n",
    encoding="utf-8",
)

# Audit against the historical 0.10.2-derived combined table. Historical values
# are never substituted into the new plot.
if cfg.OLD_COMBINED_CSV.is_file():
    historical = pd.read_csv(cfg.OLD_COMBINED_CSV)
    if {"day", "grna_target", "log_2_fold_change"}.issubset(historical.columns):
        historical = historical[historical["grna_target"].isin(cfg.TARGETS)].copy()
        historical["day"] = historical["day"].map(
            lambda x: f"d{int(x)}" if str(x).replace(".0", "").isdigit() else str(x).lower()
        )
        historical = historical[["day", "grna_target", "log_2_fold_change"]].rename(
            columns={"log_2_fold_change": "historical_v0.10.2_log_2_fold_change"}
        )
        comparison = selected.merge(historical, on=["day", "grna_target"], how="left")
        comparison["v0.10.3_minus_historical_log2fc"] = (
            comparison["log_2_fold_change"]
            - comparison["historical_v0.10.2_log_2_fold_change"]
        )
        comparison.to_csv(OUTDIR / "historical_v0.10.2_vs_v0.10.3_effect_audit.csv", index=False)

manifest_rows = []
for path in sorted(OUTDIR.iterdir()):
    if not path.is_file() or path.suffix == ".zip" or path.name == "output_manifest.csv":
        continue
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": digest})
pd.DataFrame(manifest_rows).to_csv(OUTDIR / "output_manifest.csv", index=False)

archive = OUTDIR / "sceptre_on_target_KD_D0_D3_overlay_outputs.zip"
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as package:
    for path in sorted(OUTDIR.iterdir()):
        if path.is_file() and path.suffix != ".zip":
            package.write(path, arcname=path.name)
with zipfile.ZipFile(archive) as package:
    bad_member = package.testzip()
    if bad_member is not None:
        raise IOError(f"ZIP integrity failure at {bad_member}")

print("Figure legend:\n", legend_text)
print("\nVerified output archive:", archive)
```

### Text output (stream)

```text
Figure legend:
 On-target CRISPRi knockdown of NF2, LATS2, AMOTL2 and TAOK1 across D0–D3. Time points are overlaid within each target-gene group and distinguished by symbol shape. Open symbols show the control reference, fixed at a relative-expression value of 1; filled symbols show SCEPTRE v0.10.3 fold-change estimates for cells assigned guides targeting the indicated gene. Grey bars show approximate 95% Wald confidence intervals formed on the log2 fold-change scale by applying the delta method to SCEPTRE's official fold_change and se_fold_change outputs. These model-based intervals are not obtained by inverting SCEPTRE's conditional-resampling test. Orange symbols indicate BH-adjusted SCEPTRE q ≤ 0.05; black symbols indicate q > 0.05. BH adjustment was performed separately within each time point across all QC-passing on-target discovery pairs before selecting the four displayed genes. The dashed horizontal line marks the control reference of 1.

Verified output archive: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/sceptre_on_target_KD_D0_D3_overlay_plotnine_sceptre_env_v2/sceptre_on_target_KD_D0_D3_overlay_outputs.zip

```

## Markdown cell 15

## Run instructions

1. Register/select the `sc-dl` Python kernel if needed:

```bash
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python -m ipykernel install --user \
  --name sc-dl --display-name "Python (sc-dl)"
```

2. Open this notebook, select **Python (sc-dl)** and run from the top.
3. The notebook calls the exact executable:

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/bin/Rscript
```

4. On the first run, it verifies or installs exact SCEPTRE v0.10.3 in:

```text
/oak/stanford/groups/engreitz/Users/tri/Rlibs/sceptre_0.10.3
```

This leaves the environment's v0.10.2 installation unchanged. If GitHub access is blocked on a compute node, run the preflight/install cell once on a Sherlock login node, then rerun the notebook. No Python-kernel restart is needed because SCEPTRE is launched in a new external R process.

5. Any day without a complete v0.10.3 `TargetKD_v0.10.3_with_SE.csv` is reanalyzed. Later runs reuse these files unless `FORCE_RERUN_SCEPTRE=True`.
6. The exact install script, install log, resolved package path/version, analysis R script, subprocess command and full R log are saved with the outputs.
7. The editable combined figure is written to:

```text
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/
sceptre_on_target_KD_D0_D3_overlay_plotnine_sceptre_env_v2/
sceptre_on_target_KD_D0_D3_overlay.svg
```

The corresponding PDF, 600-dpi PNG, source-data tables, legend, provenance record and verified ZIP archive are saved in the same directory.

