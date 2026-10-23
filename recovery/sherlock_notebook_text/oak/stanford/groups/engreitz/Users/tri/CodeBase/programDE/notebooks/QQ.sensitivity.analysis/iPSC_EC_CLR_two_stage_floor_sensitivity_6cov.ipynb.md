# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/iPSC_EC_CLR_two_stage_floor_sensitivity_6cov.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/iPSC_EC_CLR_two_stage_floor_sensitivity_6cov.ipynb`
- SHA-256: `80e024d58bb8b0c599901f41952495f4cc4923b17795e8c7039e90f88775b942`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# iPSC-EC two-stage CLR flooring sensitivity — six-covariate empirical Bernoulli CRT

This notebook independently evaluates both flooring operations used before the centered log-ratio (CLR) transformation of cNMF program usages:

1. **Stage 1:** day-specific replacement of exact zero usages, followed by row normalization.
2. **Stage 2:** the lower-tail `eps_quantile` floor inside `programDE.prepare_crt_inputs`, followed by a second row normalization and CLR.

The full 3 × 4 flooring grid is audited at the exact transformed-response level. The default inferential analysis runs six axis-isolating conditions: all internal quantiles at the current initial floor, and all initial-floor factors at the current internal quantile. Set `CRT_GRID_MODE = "full"` to run empirical CRT inference for all 12 combinations.

Only the **empirical Bernoulli CRT with the Bernoulli index sampler** is used. Skew-normal calibration and stratified resampling are excluded. The full six-source-covariate model is used throughout.

## Markdown cell 1

## Analysis design

The current reference condition is expected to be:

\[
\varepsilon_{0,d}=\max\left(10^{-8},\;0.1\,Q_{0.01}(U_d^{\mathrm{raw}}\mid U_d^{\mathrm{raw}}>0)\right),
\]

with a Stage-1 factor of 1×, followed by the local `programDE` internal floor at its default quantile, expected to be \(q=10^{-4}\).

The transformation audit includes:

- initial-floor factors: 0.1×, 1×, 10×;
- internal quantiles: 0, \(10^{-5}\), \(10^{-4}\), \(10^{-3}\);
- `q=0` as the Stage-1-only condition, because the Stage-1 matrix is already strictly positive and flooring at its minimum changes no entries.

For every condition, the notebook records how many entries Stage 2 modifies, whether those entries were originally positive, pseudomass added at each stage, exact CLR changes, beta and p-value stability, discovery overlap, and NTC pseudo-target calibration.

## Code cell 2

```python
from pathlib import Path
import json
import os
import shlex
import subprocess
import sys

# Limit nested BLAS/OpenMP threading before numerical libraries are imported.
for variable in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
]:
    os.environ.setdefault(variable, "1")

# Sherlock source paths used by the prior iPSC-EC notebooks.
PROGRAMDE_ROOT = "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE"
DATA_DIR = "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data"
CNMF_DIR = "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all"
GUIDE_LONG_TEMPLATE = "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/{day}_sceptre_transfer/{day}_grna_assignments_long.csv"
GUIDE_METADATA_TSV = "./300genes_guide_metadata_v43.tsv"

OUT_DIR = "results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov"
DAYS = ["d0", "d1", "d2", "d3"]
K = 50
DENSITY_THRESHOLD = "0_2"
NTC_LABELS = ["SAFE", "non-targeting", "NTC"]

# Stage 1: scale the current day-specific exact-zero replacement floor.
INITIAL_FLOOR_FACTORS = [0.1, 1.0, 10.0]

# Stage 2: programDE lower-tail floor before CLR.
# q=0 is the important "stage-1 floor only" condition because the stage-1
# matrix is already strictly positive; max(U, min(U)) changes no entries.
INTERNAL_EPS_QUANTILES = [0.0, 1e-5, 1e-4, 1e-3]

REFERENCE_INITIAL_FLOOR_FACTOR = 1.0

# None means: inspect the local prepare_crt_inputs default. The current repo
# implementation is expected to use 1e-4, and the backend records/validates it.
REFERENCE_INTERNAL_EPS_QUANTILE = None

# "primary" runs six axis-isolating CRT conditions:
#   - all q values at the current 1x initial floor
#   - all initial-floor factors at the current q
# "full" runs all 3 x 4 = 12 combinations for every day.
CRT_GRID_MODE = "primary"

# Final empirical Bernoulli CRT setting. A quick code check can temporarily use
# B=127 or B=1023, but final sensitivity results should use B=10000.
B = 10000
N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))

REQUIRE_COVARIATES = True
SKIP_CRT = False

# Stale-output protection. OVERWRITE=True is the safest default after code edits.
OVERWRITE = True
REUSE_MATCHING = False

REAL_CRT_SEED = 123
NTC_GROUP_SEED = 7
NTC_CRT_SEED = 23

MODEL_ROOT = Path(OUT_DIR) / f"k_{K}.dt_{DENSITY_THRESHOLD}"
print("Output root:", MODEL_ROOT)
print("Days:", DAYS)
print("B:", B)
print("N_JOBS:", N_JOBS)
print("Initial-floor factors:", INITIAL_FLOOR_FACTORS)
print("Internal eps quantiles:", INTERNAL_EPS_QUANTILES)
print("CRT grid mode:", CRT_GRID_MODE)

```

### Text output (stream)

```text
Output root: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2
Days: ['d0', 'd1', 'd2', 'd3']
B: 10000
N_JOBS: 16
Initial-floor factors: [0.1, 1.0, 10.0]
Internal eps quantiles: [0.0, 1e-05, 0.0001, 0.001]
CRT grid mode: primary

```

## Markdown cell 3

## Write the shared iPSC-EC data-loading and six-covariate helper backend

## Code cell 4

```python
%%writefile iPSC_EC_QQ_export_all_days.py
#!/usr/bin/env python
"""
Export iPSC-EC CRT QQ-plot p-value points for d0, d1, d2, and d3.

This version uses the full requested covariate adjustment set:
pct_counts_mt, biological_sample, total_counts, guides_per_cell,
guide_umi_counts, and n_genes_by_counts.

This script adapts the simple QQ export pattern from the IP/GV K30 dt2
bundle to the iPSC-EC notebooks.  It uses the source paths and setup used in
qq_sweep_test.ipynb and ipsc_ec_qq.ipynb:

  * H5MU source files: {day}_normalized_gex_and_guide.h5mu
  * cNMF usage files: *.usages.k_{K}.dt_{density_threshold}.consensus.txt
  * guide assignment long tables: {day}_grna_assignments_long.csv

The main exported files are per-day CSVs containing the exact p-values and the
ranked expected/observed -log10(p) coordinates needed to replot the QQ plots.
"""

from __future__ import annotations

# Set thread-count defaults before importing NumPy/SciPy.
import os

for _variable in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
]:
    os.environ.setdefault(_variable, "1")

import argparse
import inspect
import json
import re
import sys
import time
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Iterable

import numpy as np
import pandas as pd


DEFAULT_PROGRAMDE_ROOT = "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE"
DEFAULT_DATA_DIR = (
    "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
    "2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data"
)
DEFAULT_CNMF_DIR = (
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/"
    "012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all"
)
DEFAULT_GUIDE_LONG_TEMPLATE = (
    "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/"
    "{day}_sceptre_transfer/{day}_grna_assignments_long.csv"
)
DEFAULT_DAYS = ["d0", "d1", "d2", "d3"]
DEFAULT_NTC_LABELS = ["SAFE", "non-targeting", "NTC"]


@dataclass
class DayRunSummary:
    day: str
    status: str
    output_dir: str
    h5mu_path: str | None = None
    usage_path: str | None = None
    guide_long_csv: str | None = None
    n_cells: int | None = None
    n_programs: int | None = None
    n_guides: int | None = None
    n_covariates: int | None = None
    matched_usage_cells: int | None = None
    n_observed_target_empirical_points: int | None = None
    n_ntc_pseudotarget_empirical_points: int | None = None
    n_observed_target_skew_points: int | None = None
    n_ntc_pseudotarget_skew_points: int | None = None
    message: str | None = None


def safe_name(value: object) -> str:
    """Return a filesystem-safe label."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "value"


def normalize_density_token(value: object) -> str:
    """Normalize density threshold labels to the dt token used by cNMF filenames."""
    text = str(value).strip()
    text = re.sub(r"^dt[_-]?", "", text, flags=re.IGNORECASE)
    text = text.replace(".", "_").replace("-", "_")
    text = re.sub(r"_+", "_", text).strip("_")
    if not text:
        raise ValueError("Density threshold token is empty.")
    return text


def model_tag(k: int, density_token: str) -> str:
    return f"k_{k}.dt_{density_token}"


def read_table_auto(path: Path) -> pd.DataFrame:
    """Read TSV/CSV files with a conservative delimiter guess."""
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine="python")


def load_programde(programde_root: str | None) -> SimpleNamespace:
    """Import programDE functions after adding the requested repository path."""
    if programde_root:
        root = str(Path(programde_root).expanduser())
        if root not in sys.path:
            sys.path.insert(0, root)

    from src.sceptre import (  # type: ignore
        build_ntc_group_inputs,
        crt_pvals_for_ntc_groups_ensemble,
        crt_pvals_for_ntc_groups_ensemble_skew,
        make_ntc_groups_ensemble,
        prepare_crt_inputs,
        run_all_genes_union_crt,
    )

    try:
        from src.sceptre import compute_ntc_group_null_pvals_parallel  # type: ignore
    except (ImportError, AttributeError):
        compute_ntc_group_null_pvals_parallel = None

    try:
        from src.sceptre import limit_threading  # type: ignore

        limit_threading()
    except (ImportError, AttributeError):
        pass

    try:
        from src.visualization import qq_plot_ntc_pvals  # type: ignore
    except (ImportError, AttributeError):
        qq_plot_ntc_pvals = None

    return SimpleNamespace(
        build_ntc_group_inputs=build_ntc_group_inputs,
        crt_pvals_for_ntc_groups_ensemble=crt_pvals_for_ntc_groups_ensemble,
        crt_pvals_for_ntc_groups_ensemble_skew=crt_pvals_for_ntc_groups_ensemble_skew,
        make_ntc_groups_ensemble=make_ntc_groups_ensemble,
        prepare_crt_inputs=prepare_crt_inputs,
        run_all_genes_union_crt=run_all_genes_union_crt,
        compute_ntc_group_null_pvals_parallel=compute_ntc_group_null_pvals_parallel,
        qq_plot_ntc_pvals=qq_plot_ntc_pvals,
    )


def resolve_modality(mdata: Any, preferred: str | None, fallbacks: list[str], label: str) -> str:
    """Resolve a modality name in a MuData object."""
    modalities = list(mdata.mod.keys())
    candidates = []
    if preferred:
        candidates.append(preferred)
    candidates.extend(fallbacks)

    for candidate in candidates:
        if candidate in mdata.mod:
            return candidate

    lowered = {name.lower(): name for name in modalities}
    for candidate in candidates:
        if candidate and candidate.lower() in lowered:
            return lowered[candidate.lower()]

    raise KeyError(
        f"Could not find {label} modality. Tried {candidates}; available modalities: {modalities}"
    )


def infer_guide2gene_from_names(
    guide_names: Iterable[str],
    ntc_labels: list[str],
) -> dict[str, str]:
    """Infer guide-to-target names using the ipsc_ec_qq.ipynb split rule."""
    labels_lower = {label.lower(): label for label in ntc_labels}
    mapping: dict[str, str] = {}
    for guide in map(str, guide_names):
        cleaned = guide.strip()
        lowered = cleaned.lower()
        matched_control = None
        for label_lower, label in labels_lower.items():
            if (
                lowered == label_lower
                or lowered.startswith(label_lower + "_")
                or lowered.startswith(label_lower + "-")
                or lowered.startswith(label_lower + ".")
            ):
                matched_control = label
                break
        mapping[guide] = matched_control if matched_control is not None else cleaned.split("_", 1)[0]
    return mapping


def load_guide2gene_metadata(
    path: Path,
    guide_names: list[str],
    guide_col: str | None,
    target_col: str | None,
) -> dict[str, str]:
    """Optionally load guide-to-target mapping from a metadata table."""
    table = read_table_auto(path)
    normalized_columns = {column.lower().replace(" ", "_"): column for column in table.columns}

    guide_candidates = [
        guide_col,
        "guide_name",
        "guide_id",
        "guide",
        "grna",
        "grna_name",
        "grna_id",
        "sgrna",
        "protospacer",
    ]
    target_candidates = [
        target_col,
        "target_name",
        "target",
        "gene",
        "gene_symbol",
        "symbol",
        "grna_target",
    ]

    resolved_guide_col = None
    for candidate in guide_candidates:
        if not candidate:
            continue
        key = candidate.lower().replace(" ", "_")
        if key in normalized_columns:
            resolved_guide_col = normalized_columns[key]
            break

    resolved_target_col = None
    for candidate in target_candidates:
        if not candidate:
            continue
        key = candidate.lower().replace(" ", "_")
        if key in normalized_columns and normalized_columns[key] != resolved_guide_col:
            resolved_target_col = normalized_columns[key]
            break

    if resolved_guide_col is None or resolved_target_col is None:
        raise KeyError(
            "Could not infer guide/target columns in guide metadata. "
            f"Columns were: {list(table.columns)}. "
            "Use --guide-metadata-guide-col and --guide-metadata-target-col."
        )

    raw_mapping = dict(
        zip(
            table[resolved_guide_col].astype(str),
            table[resolved_target_col].astype(str),
        )
    )
    return {guide: raw_mapping[guide] for guide in guide_names if guide in raw_mapping}


def build_guide2gene(
    guide_names: list[str],
    ntc_labels: list[str],
    metadata_tsv: str | None,
    guide_col: str | None,
    target_col: str | None,
) -> dict[str, str]:
    """Build guide-to-target mapping, defaulting to guide_name.split('_', 1)[0]."""
    mapping = infer_guide2gene_from_names(guide_names, ntc_labels)

    if metadata_tsv:
        metadata_path = Path(metadata_tsv).expanduser()
        if metadata_path.exists():
            try:
                metadata_mapping = load_guide2gene_metadata(
                    metadata_path,
                    guide_names,
                    guide_col,
                    target_col,
                )
                mapping.update(metadata_mapping)
                print(
                    f"guide metadata: mapped {len(metadata_mapping):,}/{len(guide_names):,} guides "
                    f"from {metadata_path}"
                )
            except Exception as exc:  # keep fallback mapping usable
                warnings.warn(
                    f"Guide metadata {metadata_path} could not be used ({exc}). "
                    "Falling back to guide_name.split('_', 1)[0]."
                )
        else:
            warnings.warn(
                f"Guide metadata file {metadata_path} was not found. "
                "Falling back to guide_name.split('_', 1)[0]."
            )
    return mapping


def h5mu_path_for_day(data_dir: Path, day: str) -> Path:
    return data_dir / f"{day}_normalized_gex_and_guide.h5mu"


def guide_long_path_for_day(template: str, day: str) -> Path:
    return Path(template.format(day=day)).expanduser()


def usage_candidates_from_args(args: argparse.Namespace, day: str, density_token: str) -> list[Path]:
    """Resolve or discover cNMF usage-file candidates."""
    if args.cnmf_file:
        return [Path(args.cnmf_file).expanduser()]

    if args.usage_file_template:
        return [
            Path(
                args.usage_file_template.format(
                    day=day,
                    k=args.k,
                    K=args.k,
                    dt=density_token,
                    density_threshold=density_token,
                    density=density_token,
                )
            ).expanduser()
        ]

    cnmf_dir = Path(args.cnmf_dir).expanduser()
    if cnmf_dir.is_file():
        return [cnmf_dir]

    exact_patterns = [
        f"**/*.usages.k_{args.k}.dt_{density_token}.consensus.txt",
        f"**/*usages*k_{args.k}*dt_{density_token}*consensus*.txt",
        f"**/*k_{args.k}*dt_{density_token}*usage*.txt",
        f"**/*k_{args.k}*dt_{density_token}*.txt",
    ]
    candidates: list[Path] = []
    for pattern in exact_patterns:
        candidates.extend(cnmf_dir.glob(pattern))

    # De-duplicate while preserving order.
    seen = set()
    unique = []
    for path in candidates:
        resolved = path.resolve()
        if resolved not in seen:
            unique.append(path)
            seen.add(resolved)
    return unique


def score_usage_file(path: Path, obs_names: pd.Index) -> tuple[int, int]:
    """Score a usage file by overlap with the day-specific cells."""
    first_column = pd.read_csv(path, sep="\t", usecols=[0], dtype=str).iloc[:, 0]
    usage_index = pd.Index(first_column.astype(str))
    overlap = int(usage_index.intersection(obs_names.astype(str)).size)
    size_penalty = -abs(len(usage_index) - len(obs_names))
    return overlap, size_penalty


def choose_usage_file(candidates: list[Path], obs_names: pd.Index, day: str) -> tuple[Path, int]:
    """Choose the usage file with the best cell-overlap score."""
    existing = [path for path in candidates if path.exists()]
    if not existing:
        raise FileNotFoundError(
            f"No cNMF usage file exists for {day}. Candidate paths were: {candidates}"
        )

    scores = []
    for path in existing:
        try:
            overlap, size_penalty = score_usage_file(path, obs_names)
        except Exception as exc:
            warnings.warn(f"Could not score usage file {path}: {exc}")
            overlap, size_penalty = -1, -10**12
        scores.append((overlap, size_penalty, path))

    scores.sort(key=lambda item: (item[0], item[1], str(item[2])), reverse=True)
    best_overlap, _, best_path = scores[0]
    if best_overlap <= 0:
        raise ValueError(
            f"Best usage file for {day} had no overlap with the H5MU cells: {best_path}"
        )
    if len(scores) > 1:
        print("usage candidates scored for", day)
        for overlap, _, path in scores[:5]:
            print(f"  overlap={overlap:,}  {path}")
        print(f"selected usage file for {day}: {best_path}")
    return best_path, best_overlap


def floor_and_normalize_usage(usage: pd.DataFrame) -> pd.DataFrame:
    """Apply the cNMF zero flooring and row normalization used in ipsc_ec_qq.ipynb."""
    values = usage.to_numpy(dtype=float, copy=True)
    if values.ndim != 2 or values.shape[0] == 0 or values.shape[1] == 0:
        raise ValueError(f"Invalid usage matrix shape: {values.shape}")

    positive = values[values > 0]
    eps = max(1e-8, float(np.quantile(positive, 0.01)) * 0.1) if positive.size else 1e-8
    values = np.where(values > 0, values, eps)
    row_sums = values.sum(axis=1, keepdims=True)
    if np.any(row_sums <= 0) or not np.isfinite(row_sums).all():
        raise ValueError("Usage matrix has invalid row sums after flooring.")
    values = values / row_sums
    return pd.DataFrame(values, index=usage.index, columns=usage.columns)


def load_and_attach_usage(
    adata: Any,
    usage_path: Path,
) -> tuple[pd.DataFrame, int]:
    """Load, align, floor, normalize, and attach cNMF usage to adata.obsm."""
    usage_raw = pd.read_csv(usage_path, sep="\t", index_col=0)
    usage_raw.index = usage_raw.index.astype(str)
    usage_raw.columns = usage_raw.columns.astype(str)

    obs_names = pd.Index(adata.obs_names.astype(str))
    matched = int(usage_raw.index.intersection(obs_names).size)
    usage_aligned = usage_raw.reindex(obs_names).fillna(0.0)
    usage_aligned.index = adata.obs_names
    usage_aligned = floor_and_normalize_usage(usage_aligned)

    adata.obsm["cnmf_usage"] = usage_aligned
    adata.uns["program_names"] = list(usage_aligned.columns)
    return usage_aligned, matched


def load_and_attach_guides(
    adata: Any,
    mdata: Any,
    guide_modality: str,
    guide_long_csv: Path,
    ntc_labels: list[str],
    metadata_tsv: str | None,
    metadata_guide_col: str | None,
    metadata_target_col: str | None,
) -> tuple[list[str], dict[str, str]]:
    """Load guide assignments from the day-specific long CSV and attach them."""
    guide_names = list(map(str, mdata[guide_modality].var_names))
    adata.uns["guide_names"] = guide_names

    if guide_long_csv.exists():
        guide_long = pd.read_csv(guide_long_csv)
        required = {"cell_id", "grna_target", "value"}
        missing = required.difference(guide_long.columns)
        if missing:
            raise KeyError(
                f"Guide long CSV {guide_long_csv} is missing columns {sorted(missing)}; "
                f"observed columns were {list(guide_long.columns)}"
            )

        guide_dense = guide_long.pivot_table(
            index="cell_id",
            columns="grna_target",
            values="value",
            fill_value=0,
            aggfunc="max",
        )
        guide_dense.index = guide_dense.index.astype(str)
        guide_dense.columns = guide_dense.columns.astype(str)

        extra_columns = sorted(set(guide_dense.columns) - set(guide_names))
        if extra_columns:
            warnings.warn(
                f"{len(extra_columns):,} guide columns in {guide_long_csv} were not in "
                "mdata['guide'].var_names and will be ignored. First few: "
                f"{extra_columns[:5]}"
            )

        guide_aligned = guide_dense.reindex(
            index=adata.obs_names.astype(str),
            columns=guide_names,
            fill_value=0,
        )
        guide_aligned.index = adata.obs_names
        guide_aligned = (guide_aligned > 0).astype(np.int8)
        adata.obsm["guide_assignment"] = guide_aligned
    elif "guide_assignment" in adata.obsm:
        warnings.warn(
            f"Guide long CSV {guide_long_csv} was not found; using existing adata.obsm['guide_assignment']."
        )
    elif "guide_assignments" in adata.obsm:
        warnings.warn(
            f"Guide long CSV {guide_long_csv} was not found; using existing adata.obsm['guide_assignments']."
        )
        adata.obsm["guide_assignment"] = adata.obsm["guide_assignments"]
    else:
        raise FileNotFoundError(
            f"Guide long CSV {guide_long_csv} was not found and no guide assignment matrix exists in GEX.obsm."
        )

    guide2gene = build_guide2gene(
        guide_names,
        ntc_labels,
        metadata_tsv,
        metadata_guide_col,
        metadata_target_col,
    )
    adata.uns["guide2gene"] = guide2gene
    return guide_names, guide2gene


def attach_covariates(adata: Any, require_covariates: bool = False) -> pd.DataFrame:
    """Attach the full requested iPSC-EC CRT covariate set.

    Source columns in adata.obs:
      * pct_counts_mt: percent mitochondrial counts, kept on its original scale
      * biological_sample: biological sample, one-hot encoded with one reference level dropped
      * total_counts: total RNA UMI counts, log1p transformed
      * guides_per_cell: guide count burden, log1p transformed
      * guide_umi_counts: guide UMI burden, log1p transformed
      * n_genes_by_counts: detected genes per cell, log1p transformed

    The CRT backend receives only numeric covariate columns.  This avoids relying
    on downstream categorical handling and makes the exact model matrix explicit
    in each day's run_config.json.
    """
    numeric_definitions = [
        ("pct_counts_mt", "pct_counts_mt", "identity"),
        ("log1p_total_counts", "total_counts", "log1p"),
        ("log1p_guides_per_cell", "guides_per_cell", "log1p"),
        ("log1p_guide_umi_counts", "guide_umi_counts", "log1p"),
        ("log1p_n_genes_by_counts", "n_genes_by_counts", "log1p"),
    ]
    categorical_definitions = [
        ("biological_sample", "biological_sample"),
    ]

    covar = pd.DataFrame(index=adata.obs_names)
    missing_keys: list[str] = []
    covariate_sources: dict[str, dict[str, str]] = {}

    for output_name, obs_key, transform in numeric_definitions:
        if obs_key not in adata.obs:
            missing_keys.append(obs_key)
            continue

        values = pd.to_numeric(adata.obs[obs_key], errors="coerce")
        if values.isna().any():
            bad_count = int(values.isna().sum())
            raise ValueError(
                f"Covariate source {obs_key!r} has {bad_count:,} missing/non-numeric values."
            )

        numeric = values.to_numpy(dtype=float)
        if transform == "log1p":
            if np.any(numeric < 0):
                raise ValueError(f"Covariate source {obs_key!r} contains negative values.")
            numeric = np.log1p(numeric)
        elif transform != "identity":
            raise ValueError(f"Unknown transform {transform!r} for covariate {obs_key!r}.")

        if not np.isfinite(numeric).all():
            raise ValueError(f"Covariate {output_name!r} contains non-finite values.")

        covar[output_name] = numeric
        covariate_sources[output_name] = {"obs_key": obs_key, "transform": transform}

    for prefix, obs_key in categorical_definitions:
        if obs_key not in adata.obs:
            missing_keys.append(obs_key)
            continue

        values = pd.Series(adata.obs[obs_key], index=adata.obs_names)
        if values.isna().any():
            bad_count = int(values.isna().sum())
            raise ValueError(
                f"Categorical covariate source {obs_key!r} has {bad_count:,} missing values."
            )

        values = values.astype(str)
        n_levels = int(values.nunique(dropna=False))
        if n_levels <= 1:
            warnings.warn(
                f"Categorical covariate {obs_key!r} has only one level in this day; "
                "no dummy columns will be added."
            )
            continue

        dummies = pd.get_dummies(values, prefix=prefix, drop_first=True, dtype=float)
        dummies.index = adata.obs_names

        # Make dummy names file/model safe and unique after sanitization.
        sanitized = [safe_name(column) for column in dummies.columns]
        seen: dict[str, int] = {}
        unique_columns: list[str] = []
        for column in sanitized:
            count = seen.get(column, 0)
            seen[column] = count + 1
            unique_columns.append(column if count == 0 else f"{column}_{count + 1}")
        dummies.columns = unique_columns

        for column in dummies.columns:
            covariate_sources[column] = {
                "obs_key": obs_key,
                "transform": "one_hot_drop_first",
            }
        covar = pd.concat([covar, dummies], axis=1)

    if missing_keys and require_covariates:
        raise KeyError(f"Missing required covariate columns in adata.obs: {missing_keys}")
    if missing_keys:
        warnings.warn(
            f"Skipping missing covariate columns: {missing_keys}. "
            "Use --require-covariates to make this an error."
        )

    keep = [column for column in covar.columns if covar[column].nunique(dropna=False) > 1]
    dropped = [column for column in covar.columns if column not in keep]
    if dropped:
        warnings.warn(f"Dropping constant covariate columns for this day: {dropped}")

    covar = covar[keep]
    if covar.empty:
        raise ValueError("No non-constant covariate columns remain after filtering.")

    adata.obsm["covar"] = covar
    adata.uns["covariate_sources"] = {
        column: covariate_sources.get(column, {})
        for column in covar.columns
    }
    return covar
def prepare_inputs(programde: SimpleNamespace, adata: Any) -> Any:
    """Call prepare_crt_inputs with the numeric covariate handling used in the IP/GV notebook."""
    kwargs = dict(
        adata=adata,
        usage_key="cnmf_usage",
        covar_key="covar",
        guide_assignment_key="guide_assignment",
        guide2gene_key="guide2gene",
    )
    signature = inspect.signature(programde.prepare_crt_inputs)
    if "numeric_as_category_threshold" in signature.parameters:
        kwargs["numeric_as_category_threshold"] = None
    return programde.prepare_crt_inputs(**kwargs)


def pvalue_frame_to_long(
    frame: pd.DataFrame,
    *,
    day: str,
    pvalue_type: str,
    source: str,
) -> pd.DataFrame:
    """Convert a target x program p-value matrix to long format."""
    current = frame.copy()
    current.index = current.index.astype(str)
    current.columns = current.columns.astype(str)
    current.index.name = "target_name"
    long = current.reset_index().melt(
        id_vars="target_name",
        var_name="program_name",
        value_name="pvalue",
    )
    long["pvalue"] = pd.to_numeric(long["pvalue"], errors="coerce")
    long = long.loc[np.isfinite(long["pvalue"].to_numpy(dtype=float))].copy()
    long["pvalue"] = long["pvalue"].clip(1e-300, 1.0)
    long.insert(0, "source", source)
    long.insert(0, "pvalue_type", pvalue_type)
    long.insert(0, "day", day)
    return long


def add_qq_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Add QQ rank and expected/observed -log10 p-value columns."""
    if frame.empty:
        for column in [
            "rank",
            "n_points",
            "expected_pvalue",
            "expected_neglog10_p",
            "observed_neglog10_p",
        ]:
            frame[column] = []
        return frame

    current = frame.sort_values("pvalue", kind="mergesort").reset_index(drop=True)
    n = len(current)
    current["rank"] = np.arange(1, n + 1)
    current["n_points"] = n
    current["expected_pvalue"] = (current["rank"] - 0.5) / n
    current["expected_neglog10_p"] = -np.log10(current["expected_pvalue"])
    current["observed_neglog10_p"] = -np.log10(current["pvalue"])
    return current


def observed_target_qq_points(
    pvalue_frame: pd.DataFrame,
    *,
    day: str,
    pvalue_type: str,
    ntc_labels: list[str],
    source: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return real-target QQ points and observed NTC-label row points from a p-value matrix."""
    long = pvalue_frame_to_long(
        pvalue_frame,
        day=day,
        pvalue_type=pvalue_type,
        source=source,
    )
    control_keys = {label.strip().lower() for label in ntc_labels}
    is_control = long["target_name"].astype(str).str.strip().str.lower().isin(control_keys)

    real = long.loc[~is_control].copy()
    real.insert(2, "series_type", "observed_target")
    real.insert(2, "series", "Observed target p-values")
    real = add_qq_columns(real)

    control = long.loc[is_control].copy()
    control.insert(2, "series_type", "observed_control_label")
    control.insert(2, "series", "Observed NTC-label row p-values")
    control = add_qq_columns(control)
    return real, control


def ntc_frames_to_long(
    frames: dict[Any, pd.DataFrame],
    *,
    day: str,
    pvalue_type: str,
    source: str,
) -> pd.DataFrame:
    """Convert NTC pseudo-target ensemble frames to long format."""
    rows = []
    for ensemble, frame in frames.items():
        current = frame.copy()
        current.index = current.index.astype(str)
        current.columns = current.columns.astype(str)
        current.index.name = "group_id"
        long = current.reset_index().melt(
            id_vars="group_id",
            var_name="program_name",
            value_name="pvalue",
        )
        long.insert(0, "ensemble", ensemble)
        rows.append(long)

    if not rows:
        return pd.DataFrame(
            columns=[
                "day",
                "pvalue_type",
                "series",
                "series_type",
                "source",
                "ensemble",
                "group_id",
                "program_name",
                "pvalue",
                "rank",
                "n_points",
                "expected_pvalue",
                "expected_neglog10_p",
                "observed_neglog10_p",
            ]
        )

    long = pd.concat(rows, ignore_index=True, sort=False)
    long["pvalue"] = pd.to_numeric(long["pvalue"], errors="coerce")
    long = long.loc[np.isfinite(long["pvalue"].to_numpy(dtype=float))].copy()
    long["pvalue"] = long["pvalue"].clip(1e-300, 1.0)
    long.insert(0, "source", source)
    long.insert(0, "series_type", "ntc_pseudotarget")
    long.insert(0, "series", "Non-targeting pseudo-target p-values")
    long.insert(0, "pvalue_type", pvalue_type)
    long.insert(0, "day", day)
    return add_qq_columns(long)


def array_to_qq_points(
    values: np.ndarray,
    *,
    day: str,
    pvalue_type: str,
    series: str,
    series_type: str,
    source: str,
) -> pd.DataFrame:
    """Convert a flat p-value array to QQ points."""
    vector = np.asarray(values, dtype=float).ravel()
    finite = np.isfinite(vector)
    frame = pd.DataFrame({"source_index": np.flatnonzero(finite), "pvalue": vector[finite]})
    frame["pvalue"] = frame["pvalue"].clip(1e-300, 1.0)
    frame.insert(0, "source", source)
    frame.insert(0, "series_type", series_type)
    frame.insert(0, "series", series)
    frame.insert(0, "pvalue_type", pvalue_type)
    frame.insert(0, "day", day)
    return add_qq_columns(frame)


def save_ntc_frames_wide(frames: dict[Any, pd.DataFrame], path: Path) -> None:
    """Save NTC pseudo-target p-values with one row per ensemble/group."""
    rows = []
    for ensemble, frame in frames.items():
        current = frame.copy()
        current.index = current.index.astype(str)
        current.columns = current.columns.astype(str)
        current.index.name = "group_id"
        current = current.reset_index()
        current.insert(0, "ensemble", ensemble)
        rows.append(current)
    if rows:
        pd.concat(rows, ignore_index=True, sort=False).to_csv(path, sep="\t", index=False)


def save_simple_qq_plot(
    points: pd.DataFrame,
    output_stem: Path,
    *,
    day: str,
    k: int,
    density_token: str,
    pvalue_type: str = "empirical",
) -> None:
    """Save a point-only QQ plot for the empirical p-value points."""
    if points.empty:
        return

    import matplotlib.pyplot as plt

    current = points.loc[points["pvalue_type"] == pvalue_type].copy()
    current = current.loc[current["series_type"].isin(["observed_target", "ntc_pseudotarget"])]
    if current.empty:
        return

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(6.4, 6.0))
    maxima = [1.0]
    styles = [
        ("ntc_pseudotarget", "Non-targeting pseudo-targets", "#4D4D4D", 0.55, 15),
        ("observed_target", "Observed targets", "#17365D", 0.78, 18),
    ]
    for series_type, label, color, alpha, size in styles:
        subset = current.loc[current["series_type"] == series_type]
        if subset.empty:
            continue
        axis.scatter(
            subset["expected_neglog10_p"],
            subset["observed_neglog10_p"],
            s=size,
            color=color,
            alpha=alpha,
            linewidths=0,
            label=label,
            rasterized=False,
        )
        maxima.extend(
            [
                float(subset["expected_neglog10_p"].max()),
                float(subset["observed_neglog10_p"].max()),
            ]
        )

    limit = max(maxima) * 1.03
    axis.plot(
        [0, limit],
        [0, limit],
        linestyle="--",
        linewidth=1.0,
        color="#A6A6A6",
        label="Expected under null",
        zorder=0,
    )
    axis.set_xlim(0, limit)
    axis.set_ylim(0, limit)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("Expected -log10(p)")
    axis.set_ylabel("Observed -log10(p)")
    axis.set_title(f"iPSC-EC {day}: K={k}, dt={density_token}, {pvalue_type} CRT QQ")
    axis.legend(frameon=False, loc="upper left")
    figure.tight_layout()
    figure.savefig(Path(f"{output_stem}.pdf"), bbox_inches="tight")
    figure.savefig(Path(f"{output_stem}.png"), dpi=220, bbox_inches="tight")
    plt.close(figure)


def save_original_programde_qq_plot(
    programde: SimpleNamespace,
    *,
    out: dict[str, Any],
    adata: Any,
    ntc_labels: list[str],
    ntc_empirical: dict[Any, pd.DataFrame],
    ntc_skew: dict[Any, pd.DataFrame] | None,
    null_pvals: np.ndarray | None,
    output_path: Path,
    title: str,
) -> None:
    """Optionally save a qq_plot_ntc_pvals version similar to ipsc_ec_qq.ipynb."""
    if programde.qq_plot_ntc_pvals is None:
        warnings.warn("src.visualization.qq_plot_ntc_pvals is unavailable; skipping original QQ plot.")
        return

    import matplotlib.pyplot as plt

    skew_frame = out.get("pvals_df")
    if skew_frame is None:
        skew_frame = out.get("pvals_skew_df")

    ax = programde.qq_plot_ntc_pvals(
        pvals_raw_df=out["pvals_raw_df"],
        guide2gene=adata.uns["guide2gene"],
        ntc_genes=ntc_labels,
        pvals_skew_df=skew_frame,
        null_pvals=null_pvals,
        ntc_group_pvals_ens=ntc_empirical,
        ntc_group_pvals_skew_ens=ntc_skew,
        show_ntc_ensemble_band=True,
        show_all_pvals=True,
        title=title,
    )
    ax.figure.tight_layout()
    ax.figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(ax.figure)


def save_wide_frame(frame: pd.DataFrame, path: Path) -> None:
    current = frame.copy()
    current.index = current.index.astype(str)
    current.columns = current.columns.astype(str)
    current.index.name = "target_name"
    current.to_csv(path, sep="\t")


def save_long_frame(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False)


def run_day(
    args: argparse.Namespace,
    programde: SimpleNamespace,
    day: str,
    density_token: str,
    model_output_root: Path,
) -> DayRunSummary:
    """Run or skip one day and export all QQ point files."""
    day_output_dir = model_output_root / day
    day_output_dir.mkdir(parents=True, exist_ok=True)
    combined_points_path = day_output_dir / f"{day}_QQ_points.csv"

    if args.skip_if_exists and not args.overwrite and combined_points_path.exists():
        return DayRunSummary(
            day=day,
            status="skipped_existing",
            output_dir=str(day_output_dir),
            message=f"{combined_points_path} already exists",
        )

    h5mu_path = h5mu_path_for_day(Path(args.data_dir).expanduser(), day)
    guide_csv = guide_long_path_for_day(args.guide_long_csv_template, day)

    if args.dry_run:
        # We cannot select a usage file without loading H5MU obs_names, so report candidates only.
        candidates = usage_candidates_from_args(args, day, density_token)
        return DayRunSummary(
            day=day,
            status="dry_run",
            output_dir=str(day_output_dir),
            h5mu_path=str(h5mu_path),
            usage_path=";".join(str(path) for path in candidates[:10]),
            guide_long_csv=str(guide_csv),
            message="dry run only; no files written",
        )

    import muon as mu

    print(f"\n=== {day}: K={args.k}, dt={density_token} ===")
    if not h5mu_path.exists():
        raise FileNotFoundError(f"Missing H5MU for {day}: {h5mu_path}")
    mdata = mu.read(h5mu_path)
    gex_modality = resolve_modality(mdata, args.gex_modality, ["GEX", "gex", "rna", "RNA"], "GEX/RNA")
    guide_modality = resolve_modality(mdata, args.guide_modality, ["guide", "Guide", "gRNA", "grna"], "guide")
    adata = mdata[gex_modality].copy()

    usage_candidates = usage_candidates_from_args(args, day, density_token)
    usage_path, matched_usage_cells = choose_usage_file(usage_candidates, pd.Index(adata.obs_names), day)
    usage, matched_usage_cells = load_and_attach_usage(adata, usage_path)

    guide_names, _ = load_and_attach_guides(
        adata,
        mdata,
        guide_modality,
        guide_csv,
        args.ntc_label,
        args.guide_metadata_tsv,
        args.guide_metadata_guide_col,
        args.guide_metadata_target_col,
    )
    covariates = attach_covariates(adata, require_covariates=args.require_covariates)
    inputs = prepare_inputs(programde, adata)

    print(
        f"cells={adata.n_obs:,}; programs={usage.shape[1]:,}; guides={len(guide_names):,}; "
        f"usage_matched={matched_usage_cells:,}; covariates={list(covariates.columns)}"
    )

    run_started = time.time()
    out = programde.run_all_genes_union_crt(
        inputs=inputs,
        B=args.B,
        n_jobs=args.n_jobs,
        calibrate_skew_normal=not args.no_skew,
        return_raw_pvals=True,
        return_skew_normal=not args.no_skew,
    )

    pvals_empirical = out["pvals_raw_df"]
    pvals_skew = out.get("pvals_df")
    if pvals_skew is None:
        pvals_skew = out.get("pvals_skew_df")
    betas = out.get("betas_df")

    save_wide_frame(pvals_empirical, day_output_dir / f"{day}_target_pvalues_empirical_wide.tsv")
    save_long_frame(
        pvalue_frame_to_long(
            pvals_empirical,
            day=day,
            pvalue_type="empirical",
            source="run_all_genes_union_crt pvals_raw_df",
        ),
        day_output_dir / f"{day}_target_pvalues_empirical_long.tsv",
    )
    if pvals_skew is not None:
        save_wide_frame(pvals_skew, day_output_dir / f"{day}_target_pvalues_skew_wide.tsv")
        save_long_frame(
            pvalue_frame_to_long(
                pvals_skew,
                day=day,
                pvalue_type="skew",
                source="run_all_genes_union_crt pvals_df/pvals_skew_df",
            ),
            day_output_dir / f"{day}_target_pvalues_skew_long.tsv",
        )
    if betas is not None:
        save_wide_frame(betas, day_output_dir / f"{day}_target_betas_wide.tsv")

    ntc_guides, guide_freq, guide_to_bin, real_sigs = programde.build_ntc_group_inputs(
        inputs=inputs,
        ntc_label=args.ntc_label,
        group_size=args.ntc_group_size,
        n_bins=args.ntc_bins,
    )
    groups = programde.make_ntc_groups_ensemble(
        ntc_guides=ntc_guides,
        ntc_freq=guide_freq,
        real_gene_bin_sigs=real_sigs,
        guide_to_bin=guide_to_bin,
        n_ensemble=args.ntc_ensembles,
        seed0=args.ntc_group_seed,
        group_size=args.ntc_group_size,
        max_groups=args.ntc_max_groups,
    )
    ntc_empirical = programde.crt_pvals_for_ntc_groups_ensemble(
        inputs=inputs,
        ntc_groups_ens=groups,
        B=args.B,
        seed0=args.ntc_crt_seed,
    )
    save_ntc_frames_wide(ntc_empirical, day_output_dir / f"{day}_ntc_pseudotarget_pvalues_empirical_wide.tsv")
    np.save(
        day_output_dir / f"{day}_ntc_pseudotarget_pvalues_empirical.npy",
        np.concatenate([
            frame.to_numpy(dtype=float, copy=False).ravel()
            for frame in ntc_empirical.values()
        ]) if ntc_empirical else np.empty(0, dtype=float),
    )

    ntc_skew = None
    if not args.no_skew:
        ntc_skew = programde.crt_pvals_for_ntc_groups_ensemble_skew(
            inputs=inputs,
            ntc_groups_ens=groups,
            B=args.B,
            seed0=args.ntc_crt_seed,
        )
        save_ntc_frames_wide(ntc_skew, day_output_dir / f"{day}_ntc_pseudotarget_pvalues_skew_wide.tsv")
        np.save(
            day_output_dir / f"{day}_ntc_pseudotarget_pvalues_skew.npy",
            np.concatenate([
                frame.to_numpy(dtype=float, copy=False).ravel()
                for frame in ntc_skew.values()
            ]) if ntc_skew else np.empty(0, dtype=float),
        )

    null_pvals = None
    if args.compute_null:
        if programde.compute_ntc_group_null_pvals_parallel is None:
            warnings.warn("compute_ntc_group_null_pvals_parallel unavailable; skipping CRT-null export.")
        else:
            null_pvals = programde.compute_ntc_group_null_pvals_parallel(
                inputs=inputs,
                ntc_groups_ens=groups,
                B=args.B,
                n_jobs=args.n_jobs,
                backend=args.null_backend,
            )
            null_pvals = np.asarray(null_pvals, dtype=float).ravel()
            np.save(day_output_dir / f"{day}_crt_null_pvalues.npy", null_pvals)

    real_empirical, control_empirical = observed_target_qq_points(
        pvals_empirical,
        day=day,
        pvalue_type="empirical",
        ntc_labels=args.ntc_label,
        source="observed target rows from pvals_raw_df",
    )
    ntc_empirical_points = ntc_frames_to_long(
        ntc_empirical,
        day=day,
        pvalue_type="empirical",
        source="grouped non-targeting pseudo-target CRT p-values",
    )

    point_frames = [real_empirical, ntc_empirical_points]
    if not control_empirical.empty:
        point_frames.append(control_empirical)

    real_skew = pd.DataFrame()
    ntc_skew_points = pd.DataFrame()
    control_skew = pd.DataFrame()
    if pvals_skew is not None:
        real_skew, control_skew = observed_target_qq_points(
            pvals_skew,
            day=day,
            pvalue_type="skew",
            ntc_labels=args.ntc_label,
            source="observed target rows from pvals_df/pvals_skew_df",
        )
        point_frames.append(real_skew)
        if not control_skew.empty:
            point_frames.append(control_skew)
    if ntc_skew is not None:
        ntc_skew_points = ntc_frames_to_long(
            ntc_skew,
            day=day,
            pvalue_type="skew",
            source="grouped non-targeting pseudo-target skew-calibrated CRT p-values",
        )
        point_frames.append(ntc_skew_points)

    if null_pvals is not None:
        point_frames.append(
            array_to_qq_points(
                null_pvals,
                day=day,
                pvalue_type="crt_null",
                series="CRT null p-values",
                series_type="crt_null",
                source="compute_ntc_group_null_pvals_parallel",
            )
        )

    real_empirical.to_csv(day_output_dir / f"{day}_observed_target_QQ_points.csv", index=False)
    ntc_empirical_points.to_csv(day_output_dir / f"{day}_nontargeting_pseudotarget_QQ_points.csv", index=False)
    if not control_empirical.empty:
        control_empirical.to_csv(day_output_dir / f"{day}_observed_ntc_label_row_QQ_points.csv", index=False)
    if not real_skew.empty:
        real_skew.to_csv(day_output_dir / f"{day}_observed_target_skew_QQ_points.csv", index=False)
    if not ntc_skew_points.empty:
        ntc_skew_points.to_csv(day_output_dir / f"{day}_nontargeting_pseudotarget_skew_QQ_points.csv", index=False)

    combined_points = pd.concat(point_frames, ignore_index=True, sort=False)
    combined_points.to_csv(combined_points_path, index=False)

    if not args.no_plot:
        save_simple_qq_plot(
            combined_points,
            day_output_dir / f"{day}_QQ_empirical_simple",
            day=day,
            k=args.k,
            density_token=density_token,
            pvalue_type="empirical",
        )
        if not args.no_skew:
            save_simple_qq_plot(
                combined_points,
                day_output_dir / f"{day}_QQ_skew_simple",
                day=day,
                k=args.k,
                density_token=density_token,
                pvalue_type="skew",
            )

    if args.save_original_qq:
        if null_pvals is None and args.compute_null is False:
            warnings.warn(
                "--save-original-qq works best with --compute-null. Saving without CRT-null p-values."
            )
        save_original_programde_qq_plot(
            programde,
            out=out,
            adata=adata,
            ntc_labels=args.ntc_label,
            ntc_empirical=ntc_empirical,
            ntc_skew=ntc_skew,
            null_pvals=null_pvals,
            output_path=day_output_dir / f"{day}_programDE_original_QQ.png",
            title=f"iPSC-EC {day}: K={args.k}, dt={density_token} CRT QQ",
        )

    config = {
        "day": day,
        "k": args.k,
        "density_threshold_token": density_token,
        "B": args.B,
        "n_jobs": args.n_jobs,
        "h5mu_path": str(h5mu_path),
        "gex_modality": gex_modality,
        "guide_modality": guide_modality,
        "usage_path": str(usage_path),
        "guide_long_csv": str(guide_csv),
        "ntc_labels": args.ntc_label,
        "ntc_group_size": args.ntc_group_size,
        "ntc_ensembles": args.ntc_ensembles,
        "ntc_bins": args.ntc_bins,
        "ntc_max_groups": args.ntc_max_groups,
        "ntc_group_seed": args.ntc_group_seed,
        "ntc_crt_seed": args.ntc_crt_seed,
        "covariate_columns": list(covariates.columns),
        "covariate_sources": adata.uns.get("covariate_sources", {}),
        "pvalue_outputs": {
            "observed_target_QQ_points": f"{day}_observed_target_QQ_points.csv",
            "nontargeting_pseudotarget_QQ_points": f"{day}_nontargeting_pseudotarget_QQ_points.csv",
            "combined_QQ_points": f"{day}_QQ_points.csv",
            "target_empirical_wide": f"{day}_target_pvalues_empirical_wide.tsv",
            "target_skew_wide": f"{day}_target_pvalues_skew_wide.tsv" if pvals_skew is not None else None,
            "ntc_empirical_wide": f"{day}_ntc_pseudotarget_pvalues_empirical_wide.tsv",
            "ntc_skew_wide": f"{day}_ntc_pseudotarget_pvalues_skew_wide.tsv" if ntc_skew is not None else None,
        },
        "elapsed_minutes": round((time.time() - run_started) / 60, 3),
    }
    (day_output_dir / f"{day}_run_config.json").write_text(json.dumps(config, indent=2) + "\n")

    print(
        f"saved {day}: observed target empirical points={len(real_empirical):,}; "
        f"NTC empirical points={len(ntc_empirical_points):,}; output={combined_points_path}"
    )

    return DayRunSummary(
        day=day,
        status="completed",
        output_dir=str(day_output_dir),
        h5mu_path=str(h5mu_path),
        usage_path=str(usage_path),
        guide_long_csv=str(guide_csv),
        n_cells=int(adata.n_obs),
        n_programs=int(usage.shape[1]),
        n_guides=int(len(guide_names)),
        n_covariates=int(covariates.shape[1]),
        matched_usage_cells=int(matched_usage_cells),
        n_observed_target_empirical_points=int(len(real_empirical)),
        n_ntc_pseudotarget_empirical_points=int(len(ntc_empirical_points)),
        n_observed_target_skew_points=int(len(real_skew)) if not real_skew.empty else 0,
        n_ntc_pseudotarget_skew_points=int(len(ntc_skew_points)) if not ntc_skew_points.empty else 0,
    )


def consolidate_outputs(model_output_root: Path, summaries: list[DayRunSummary]) -> None:
    """Create combined all-day QQ CSV and manifest."""
    point_tables = []
    manifest_rows = []
    for summary in summaries:
        row = asdict(summary)
        day_dir = Path(summary.output_dir)
        points_path = day_dir / f"{summary.day}_QQ_points.csv"
        row["qq_points_csv"] = str(points_path) if points_path.exists() else None
        if points_path.exists():
            points = pd.read_csv(points_path)
            points.insert(0, "source_csv", str(points_path.relative_to(model_output_root)))
            point_tables.append(points)
            row["qq_points_rows"] = int(len(points))
        else:
            row["qq_points_rows"] = 0
        manifest_rows.append(row)

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(model_output_root / "QQ_export_manifest.csv", index=False)

    if point_tables:
        combined = pd.concat(point_tables, ignore_index=True, sort=False)
        combined.to_csv(model_output_root / "all_days_QQ_points.csv", index=False)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run iPSC-EC CRT for d0-d3 and export observed-target and non-targeting QQ points."
    )
    parser.add_argument("--programde-root", default=DEFAULT_PROGRAMDE_ROOT)
    parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR)
    parser.add_argument("--cnmf-dir", default=DEFAULT_CNMF_DIR)
    parser.add_argument("--cnmf-file", default=None, help="Use one exact cNMF usage file for all days.")
    parser.add_argument(
        "--usage-file-template",
        default=None,
        help=(
            "Optional template for usage files. Available fields: {day}, {k}, {K}, {dt}, "
            "{density_threshold}. If unset, files are discovered under --cnmf-dir."
        ),
    )
    parser.add_argument("--guide-long-csv-template", default=DEFAULT_GUIDE_LONG_TEMPLATE)
    parser.add_argument(
        "--guide-metadata-tsv",
        default="./300genes_guide_metadata_v43.tsv",
        help="Optional guide metadata table. If absent, guide names are mapped with split('_', 1)[0].",
    )
    parser.add_argument("--guide-metadata-guide-col", default=None)
    parser.add_argument("--guide-metadata-target-col", default=None)
    parser.add_argument("--out-dir", default="results/ipsc_ec_qq_exports")
    parser.add_argument("--days", nargs="+", default=DEFAULT_DAYS)
    parser.add_argument("--k", type=int, default=50)
    parser.add_argument(
        "--density-threshold",
        default="0_2",
        help="cNMF density threshold token or number, e.g. 0_2, 0.2, 2_0, or 2.0.",
    )
    parser.add_argument("--gex-modality", default=None)
    parser.add_argument("--guide-modality", default=None)
    parser.add_argument("--B", type=int, default=10000)
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")),
    )
    parser.add_argument("--ntc-label", action="append", default=None)
    parser.add_argument("--ntc-group-size", type=int, default=6)
    parser.add_argument("--ntc-ensembles", type=int, default=10)
    parser.add_argument("--ntc-bins", type=int, default=10)
    parser.add_argument("--ntc-max-groups", type=int, default=None)
    parser.add_argument("--ntc-group-seed", type=int, default=7)
    parser.add_argument("--ntc-crt-seed", type=int, default=23)
    parser.add_argument("--no-skew", action="store_true", help="Do not compute skew-calibrated p-values.")
    parser.add_argument("--compute-null", action="store_true", help="Also compute CRT-null p-values for NTC groups.")
    parser.add_argument("--null-backend", default="threading")
    parser.add_argument("--require-covariates", action="store_true")
    parser.add_argument("--skip-if-exists", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--save-original-qq", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    if args.ntc_label is None:
        args.ntc_label = DEFAULT_NTC_LABELS.copy()
    if args.B < 1:
        raise ValueError("--B must be >= 1")
    if args.n_jobs < 1:
        raise ValueError("--n-jobs must be >= 1")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    density_token = normalize_density_token(args.density_threshold)
    output_root = Path(args.out_dir).expanduser()
    model_output_root = output_root / model_tag(args.k, density_token)
    model_output_root.mkdir(parents=True, exist_ok=True)

    print("iPSC-EC QQ export")
    print("days:", args.days)
    print("model:", model_tag(args.k, density_token))
    print("data_dir:", args.data_dir)
    print("cnmf_dir:", args.cnmf_dir)
    print("out_dir:", model_output_root)
    print("B:", args.B, "n_jobs:", args.n_jobs)

    programde = None if args.dry_run else load_programde(args.programde_root)

    summaries: list[DayRunSummary] = []
    for day in args.days:
        try:
            summary = run_day(args, programde, day, density_token, model_output_root)  # type: ignore[arg-type]
        except Exception as exc:
            if args.dry_run:
                raise
            summary = DayRunSummary(
                day=day,
                status="failed",
                output_dir=str(model_output_root / day),
                message=repr(exc),
            )
            summaries.append(summary)
            consolidate_outputs(model_output_root, summaries)
            raise
        summaries.append(summary)

    consolidate_outputs(model_output_root, summaries)
    print("\nCombined outputs:")
    print("  manifest:", model_output_root / "QQ_export_manifest.csv")
    if (model_output_root / "all_days_QQ_points.csv").exists():
        print("  all days QQ points:", model_output_root / "all_days_QQ_points.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

```

### Text output (stream)

```text
Overwriting iPSC_EC_QQ_export_all_days.py

```

## Markdown cell 5

## Write the two-stage flooring sensitivity backend

## Code cell 6

```python
%%writefile iPSC_EC_CLR_two_stage_floor_sensitivity_6cov.py
#!/usr/bin/env python
"""Two-stage CLR flooring sensitivity analysis for iPSC-EC programDE CRT.

This backend uses the six-source-covariate full-day model:
  * pct_counts_mt (identity)
  * biological_sample (one-hot encoded, one reference level dropped)
  * total_counts (log1p)
  * guides_per_cell (log1p)
  * guide_umi_counts (log1p)
  * n_genes_by_counts (log1p)

It independently varies:
  1. The initial exact-zero replacement floor applied to raw cNMF usages.
  2. The internal programDE lower-tail quantile floor used before CLR.

The transformation audit always evaluates the full Cartesian grid. The CRT can
run either an axis-isolating primary grid (all internal quantiles at the current
initial floor plus all initial-floor factors at the current internal quantile)
or the complete factorial grid. Only the empirical Bernoulli CRT is used; no
skew-normal or stratified-resampling analysis is performed.
"""

from __future__ import annotations

import os

for _variable in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
]:
    os.environ.setdefault(_variable, "1")

import argparse
import gc
import hashlib
import inspect
import json
import subprocess
import time
import warnings
from dataclasses import asdict, dataclass
from itertools import product
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

import iPSC_EC_QQ_export_all_days as qqb


DEFAULT_DAYS = ["d0", "d1", "d2", "d3"]
DEFAULT_NTC_LABELS = ["SAFE", "non-targeting", "NTC"]


@dataclass
class RunRecord:
    day: str
    initial_floor_factor: float
    internal_eps_quantile: float
    initial_floor_label: str
    internal_quantile_label: str
    combination_label: str
    is_reference: bool
    status: str
    output_dir: str
    n_cells: int | None = None
    n_programs: int | None = None
    n_guides: int | None = None
    n_covariate_columns: int | None = None
    n_real_target_program_pairs: int | None = None
    n_ntc_pvalues_empirical: int | None = None
    raw_usage_zero_fraction: float | None = None
    base_raw_zero_floor: float | None = None
    applied_raw_zero_floor: float | None = None
    effective_internal_floor: float | None = None
    stage2_fraction_entries_modified: float | None = None
    stage2_fraction_positive_entries_modified: float | None = None
    stage2_fraction_zero_derived_entries_modified: float | None = None
    median_stage2_pseudomass_fraction: float | None = None
    runtime_seconds: float | None = None
    config_fingerprint: str | None = None
    message: str | None = None


def safe_float_token(value: float) -> str:
    text = f"{float(value):.12g}"
    return (
        text.replace("-", "m")
        .replace("+", "p")
        .replace(".", "_")
    )


def initial_floor_label(factor: float) -> str:
    return f"{safe_float_token(factor)}x"


def internal_quantile_label(quantile: float) -> str:
    value = float(quantile)
    if float_equal(value, 0.0):
        return "q0"
    if abs(value) < 0.01:
        text = f"{value:.0e}"
        text = text.replace("-", "m").replace("+", "p")
        return f"q{text}"
    return f"q{safe_float_token(value)}"


def combination_label(initial_factor: float, internal_quantile: float) -> str:
    return (
        f"initial_{initial_floor_label(initial_factor)}"
        f"__internal_{internal_quantile_label(internal_quantile)}"
    )


def float_equal(a: float, b: float, atol: float = 1e-15) -> bool:
    return bool(np.isclose(float(a), float(b), rtol=0.0, atol=atol))


def unique_floats(values: Iterable[float]) -> list[float]:
    result: list[float] = []
    for value in values:
        current = float(value)
        if not any(float_equal(current, existing) for existing in result):
            result.append(current)
    return result


def json_normalize(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, set):
        return sorted(json_normalize(v) for v in value)
    if isinstance(value, Mapping):
        return {str(k): json_normalize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_normalize(v) for v in value]
    return value


def config_fingerprint(config: Mapping[str, Any]) -> str:
    payload = json.dumps(
        json_normalize(config),
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(json_normalize(payload), indent=2, sort_keys=True) + "\n"
    )
    temporary.replace(path)


def git_commit(path: str | Path) -> str | None:
    try:
        completed = subprocess.run(
            ["git", "-C", str(Path(path).expanduser()), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return None
    commit = completed.stdout.strip()
    return commit or None


def bh_qvalues(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float).ravel()
    q = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return q
    pv = np.clip(p[finite], 0.0, 1.0)
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q_finite = np.empty_like(pv)
    q_finite[order] = np.clip(ranked, 0.0, 1.0)
    q[finite] = q_finite
    return q


def get_prepare_eps_quantile(programde: SimpleNamespace) -> tuple[bool, float | None]:
    signature = inspect.signature(programde.prepare_crt_inputs)
    parameter = signature.parameters.get("eps_quantile")
    if parameter is None:
        return False, None
    if parameter.default is inspect.Parameter.empty:
        return True, None
    try:
        return True, float(parameter.default)
    except (TypeError, ValueError):
        return True, None


def resolve_reference_quantile(
    requested: float | None,
    local_default: float | None,
) -> float:
    if requested is not None:
        reference = float(requested)
    elif local_default is not None:
        reference = float(local_default)
    else:
        raise RuntimeError(
            "Could not infer the current programDE eps_quantile default. "
            "Supply --reference-internal-eps-quantile explicitly."
        )
    if not np.isfinite(reference) or not 0.0 <= reference <= 1.0:
        raise ValueError(
            f"Reference internal eps quantile must be in [0, 1]; got {reference}."
        )
    return reference


def build_audit_combinations(
    initial_factors: Sequence[float],
    internal_quantiles: Sequence[float],
    reference_initial: float,
    reference_quantile: float,
) -> list[tuple[float, float]]:
    combinations = [
        (float(initial), float(quantile))
        for initial, quantile in product(initial_factors, internal_quantiles)
    ]
    combinations.sort(
        key=lambda pair: (
            0
            if float_equal(pair[0], reference_initial)
            and float_equal(pair[1], reference_quantile)
            else 1,
            pair[0],
            pair[1],
        )
    )
    return combinations


def build_crt_combinations(
    initial_factors: Sequence[float],
    internal_quantiles: Sequence[float],
    reference_initial: float,
    reference_quantile: float,
    mode: str,
) -> list[tuple[float, float]]:
    if mode == "full":
        return build_audit_combinations(
            initial_factors,
            internal_quantiles,
            reference_initial,
            reference_quantile,
        )
    if mode != "primary":
        raise ValueError(f"Unknown CRT grid mode: {mode!r}")

    combinations: list[tuple[float, float]] = []

    # Hold the initial floor fixed and vary the internal quantile.
    for quantile in internal_quantiles:
        combinations.append((float(reference_initial), float(quantile)))

    # Hold the internal quantile fixed and vary the initial floor.
    for initial in initial_factors:
        pair = (float(initial), float(reference_quantile))
        if not any(
            float_equal(pair[0], existing[0])
            and float_equal(pair[1], existing[1])
            for existing in combinations
        ):
            combinations.append(pair)

    combinations.sort(
        key=lambda pair: (
            0
            if float_equal(pair[0], reference_initial)
            and float_equal(pair[1], reference_quantile)
            else 1,
            0 if float_equal(pair[0], reference_initial) else 1,
            pair[0],
            pair[1],
        )
    )
    return combinations


def floor_and_normalize_raw_usage(
    usage_raw: pd.DataFrame,
    factor: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Replace exact zeros with a scaled day-specific floor, then row-normalize."""
    if not np.isfinite(factor) or factor <= 0:
        raise ValueError(f"Initial floor factor must be positive and finite; got {factor!r}")

    values = usage_raw.to_numpy(dtype=float, copy=True)
    if values.ndim != 2 or values.size == 0:
        raise ValueError(f"Invalid usage matrix shape: {values.shape}")
    if not np.isfinite(values).all():
        raise ValueError("Raw usage matrix contains non-finite values.")
    if np.any(values < 0):
        raise ValueError("Raw usage matrix contains negative values.")

    positive = values[values > 0]
    positive_q01 = float(np.quantile(positive, 0.01)) if positive.size else np.nan
    base_floor = max(1e-8, positive_q01 * 0.1) if positive.size else 1e-8
    applied_floor = base_floor * float(factor)

    zero_mask = values == 0
    zeros_per_cell = zero_mask.sum(axis=1)
    raw_row_sums = values.sum(axis=1)
    if np.any(raw_row_sums <= 0):
        bad = int((raw_row_sums <= 0).sum())
        raise ValueError(f"Raw usage contains {bad:,} all-zero rows.")

    floored_values = np.where(zero_mask, applied_floor, values)
    row_sums = floored_values.sum(axis=1, keepdims=True)
    if np.any(row_sums <= 0) or not np.isfinite(row_sums).all():
        raise ValueError("Invalid row sums after initial exact-zero replacement.")
    normalized = floored_values / row_sums

    initial_added_mass = zeros_per_cell * applied_floor
    initial_pseudomass_fraction = initial_added_mass / (
        raw_row_sums + initial_added_mass
    )

    zero_derived_values = normalized[zero_mask]

    metadata = {
        "initial_floor_factor": float(factor),
        "initial_floor_label": initial_floor_label(factor),
        "n_cells": int(normalized.shape[0]),
        "n_programs": int(normalized.shape[1]),
        "n_positive_raw_entries": int(positive.size),
        "n_exact_zero_raw_entries": int(zero_mask.sum()),
        "fraction_exact_zero_raw_entries": float(zero_mask.mean()),
        "pct_cells_with_any_raw_zero": float(100.0 * (zeros_per_cell > 0).mean()),
        "median_raw_zeros_per_cell": float(np.median(zeros_per_cell)),
        "p90_raw_zeros_per_cell": float(np.quantile(zeros_per_cell, 0.90)),
        "p95_raw_zeros_per_cell": float(np.quantile(zeros_per_cell, 0.95)),
        "max_raw_zeros_per_cell": int(zeros_per_cell.max()),
        "positive_raw_usage_q01": positive_q01,
        "base_raw_zero_floor": float(base_floor),
        "applied_raw_zero_floor": float(applied_floor),
        "median_initial_pseudomass_fraction": float(
            np.median(initial_pseudomass_fraction)
        ),
        "p95_initial_pseudomass_fraction": float(
            np.quantile(initial_pseudomass_fraction, 0.95)
        ),
        "max_initial_pseudomass_fraction": float(
            np.max(initial_pseudomass_fraction)
        ),
        "median_zero_derived_value_after_stage1": float(
            np.median(zero_derived_values)
        )
        if zero_derived_values.size
        else np.nan,
        "p05_zero_derived_value_after_stage1": float(
            np.quantile(zero_derived_values, 0.05)
        )
        if zero_derived_values.size
        else np.nan,
        "p95_zero_derived_value_after_stage1": float(
            np.quantile(zero_derived_values, 0.95)
        )
        if zero_derived_values.size
        else np.nan,
    }

    return (
        pd.DataFrame(
            normalized,
            index=usage_raw.index,
            columns=usage_raw.columns,
        ),
        metadata,
    )


def stage2_floor_metrics_and_manual_clr(
    stage1_usage: pd.DataFrame,
    usage_raw: pd.DataFrame,
    internal_quantile: float,
) -> tuple[np.ndarray, dict[str, Any], pd.DataFrame, dict[str, np.ndarray]]:
    """Reproduce programDE's internal floor and return stage-specific diagnostics."""
    quantile = float(internal_quantile)
    if not np.isfinite(quantile) or not 0.0 <= quantile <= 1.0:
        raise ValueError(f"Internal eps quantile must be in [0, 1]; got {quantile}.")

    values = stage1_usage.to_numpy(dtype=float, copy=False)
    raw_values = usage_raw.to_numpy(dtype=float, copy=False)
    if values.shape != raw_values.shape:
        raise ValueError("Stage-1 and raw usage matrices have different shapes.")

    internal_floor = float(np.quantile(values, quantile))
    if not np.isfinite(internal_floor) or internal_floor <= 0.0:
        internal_floor = float(np.finfo(np.float64).tiny)

    raw_zero_mask = raw_values == 0
    raw_positive_mask = raw_values > 0
    modified_mask = values < internal_floor
    modified_positive_mask = modified_mask & raw_positive_mask
    modified_zero_derived_mask = modified_mask & raw_zero_mask

    added = np.maximum(internal_floor - values, 0.0)
    added_per_cell = added.sum(axis=1)
    stage2_pseudomass_fraction = added_per_cell / (1.0 + added_per_cell)

    floored = np.maximum(values, internal_floor)
    floored = floored / floored.sum(axis=1, keepdims=True)
    log_values = np.log(floored)
    manual_clr = log_values - log_values.mean(axis=1, keepdims=True)

    modified_per_cell = modified_mask.sum(axis=1)
    modified_positive_per_cell = modified_positive_mask.sum(axis=1)
    modified_zero_per_cell = modified_zero_derived_mask.sum(axis=1)

    zero_derived_values = values[raw_zero_mask]
    modified_positive_values = values[modified_positive_mask]
    positive_uplift = (
        internal_floor / modified_positive_values
        if modified_positive_values.size
        else np.empty(0, dtype=float)
    )

    median_zero_value = (
        float(np.median(zero_derived_values))
        if zero_derived_values.size
        else np.nan
    )

    summary = {
        "internal_eps_quantile": quantile,
        "internal_quantile_label": internal_quantile_label(quantile),
        "effective_internal_floor": internal_floor,
        "stage2_n_entries_modified": int(modified_mask.sum()),
        "stage2_fraction_entries_modified": float(modified_mask.mean()),
        "stage2_n_originally_positive_entries_modified": int(
            modified_positive_mask.sum()
        ),
        "stage2_fraction_originally_positive_entries_modified": float(
            modified_positive_mask.sum() / raw_positive_mask.sum()
        )
        if raw_positive_mask.any()
        else np.nan,
        "stage2_n_zero_derived_entries_modified": int(
            modified_zero_derived_mask.sum()
        ),
        "stage2_fraction_zero_derived_entries_modified": float(
            modified_zero_derived_mask.sum() / raw_zero_mask.sum()
        )
        if raw_zero_mask.any()
        else np.nan,
        "stage2_pct_cells_with_any_entry_modified": float(
            100.0 * (modified_per_cell > 0).mean()
        ),
        "stage2_pct_cells_with_any_originally_positive_entry_modified": float(
            100.0 * (modified_positive_per_cell > 0).mean()
        ),
        "stage2_median_entries_modified_per_cell": float(
            np.median(modified_per_cell)
        ),
        "stage2_p95_entries_modified_per_cell": float(
            np.quantile(modified_per_cell, 0.95)
        ),
        "stage2_median_positive_entries_modified_per_cell": float(
            np.median(modified_positive_per_cell)
        ),
        "stage2_p95_positive_entries_modified_per_cell": float(
            np.quantile(modified_positive_per_cell, 0.95)
        ),
        "median_stage2_pseudomass_fraction": float(
            np.median(stage2_pseudomass_fraction)
        ),
        "p95_stage2_pseudomass_fraction": float(
            np.quantile(stage2_pseudomass_fraction, 0.95)
        ),
        "max_stage2_pseudomass_fraction": float(
            np.max(stage2_pseudomass_fraction)
        ),
        "median_zero_derived_value_before_stage2": median_zero_value,
        "internal_floor_over_median_zero_derived_value": float(
            internal_floor / median_zero_value
        )
        if np.isfinite(median_zero_value) and median_zero_value > 0
        else np.nan,
        "median_uplift_factor_for_modified_positive_entries": float(
            np.median(positive_uplift)
        )
        if positive_uplift.size
        else np.nan,
        "p95_uplift_factor_for_modified_positive_entries": float(
            np.quantile(positive_uplift, 0.95)
        )
        if positive_uplift.size
        else np.nan,
    }

    programs = list(map(str, stage1_usage.columns))
    program_rows: list[dict[str, Any]] = []
    for j, program in enumerate(programs):
        raw_pos_j = raw_positive_mask[:, j]
        raw_zero_j = raw_zero_mask[:, j]
        modified_j = modified_mask[:, j]
        modified_pos_j = modified_positive_mask[:, j]
        modified_zero_j = modified_zero_derived_mask[:, j]
        added_j = added[:, j]
        program_rows.append(
            {
                "program_name": program,
                "stage2_fraction_entries_modified": float(modified_j.mean()),
                "stage2_fraction_originally_positive_entries_modified": float(
                    modified_pos_j.sum() / raw_pos_j.sum()
                )
                if raw_pos_j.any()
                else np.nan,
                "stage2_fraction_zero_derived_entries_modified": float(
                    modified_zero_j.sum() / raw_zero_j.sum()
                )
                if raw_zero_j.any()
                else np.nan,
                "stage2_total_added_mass_before_renormalization": float(
                    added_j.sum()
                ),
                "stage2_median_added_mass_per_cell": float(np.median(added_j)),
                "stage2_p95_added_mass_per_cell": float(
                    np.quantile(added_j, 0.95)
                ),
            }
        )

    cell_diagnostics = {
        "raw_zeros_per_cell": raw_zero_mask.sum(axis=1).astype(float),
        "stage2_modified_entries_per_cell": modified_per_cell.astype(float),
        "stage2_modified_positive_entries_per_cell": modified_positive_per_cell.astype(
            float
        ),
        "stage2_modified_zero_derived_entries_per_cell": modified_zero_per_cell.astype(
            float
        ),
        "stage2_pseudomass_fraction_per_cell": stage2_pseudomass_fraction.astype(
            float
        ),
    }

    return manual_clr, summary, pd.DataFrame(program_rows), cell_diagnostics


def align_raw_usage(
    usage_path: Path,
    obs_names: pd.Index,
) -> tuple[pd.DataFrame, int]:
    usage_raw = pd.read_csv(usage_path, sep="\t", index_col=0)
    usage_raw.index = usage_raw.index.astype(str)
    usage_raw.columns = usage_raw.columns.astype(str)
    obs_as_str = pd.Index(obs_names.astype(str))
    matched = int(usage_raw.index.intersection(obs_as_str).size)

    aligned = usage_raw.reindex(obs_as_str)
    fully_missing = aligned.isna().all(axis=1)
    if fully_missing.any():
        examples = list(aligned.index[fully_missing][:5])
        raise ValueError(
            f"{int(fully_missing.sum()):,} H5MU cells are missing entirely from "
            f"{usage_path}; first few: {examples}."
        )
    if aligned.isna().any().any():
        n_missing = int(aligned.isna().sum().sum())
        raise ValueError(f"Usage alignment produced {n_missing:,} partial missing entries.")

    aligned.index = obs_names
    return aligned, matched


def prepare_inputs(
    programde: SimpleNamespace,
    adata: Any,
    eps_quantile: float,
) -> Any:
    signature = inspect.signature(programde.prepare_crt_inputs)
    if "eps_quantile" not in signature.parameters:
        raise TypeError(
            "The local programDE prepare_crt_inputs() has no eps_quantile parameter."
        )
    kwargs: dict[str, Any] = {
        "adata": adata,
        "usage_key": "cnmf_usage",
        "covar_key": "covar",
        "guide_assignment_key": "guide_assignment",
        "guide2gene_key": "guide2gene",
        "eps_quantile": float(eps_quantile),
    }
    if "numeric_as_category_threshold" in signature.parameters:
        kwargs["numeric_as_category_threshold"] = None
    return programde.prepare_crt_inputs(**kwargs)


def iter_object_fields(obj: Any) -> list[tuple[str, Any]]:
    fields: list[tuple[str, Any]] = []
    if isinstance(obj, Mapping):
        fields.extend((str(k), v) for k, v in obj.items())
    if hasattr(obj, "_asdict"):
        try:
            fields.extend((str(k), v) for k, v in obj._asdict().items())
        except Exception:
            pass
    if hasattr(obj, "__dataclass_fields__"):
        for name in getattr(obj, "__dataclass_fields__", {}):
            try:
                fields.append((str(name), getattr(obj, name)))
            except Exception:
                pass
    try:
        fields.extend((str(k), v) for k, v in vars(obj).items())
    except TypeError:
        pass

    seen: set[str] = set()
    unique: list[tuple[str, Any]] = []
    for name, value in fields:
        if name in seen:
            continue
        seen.add(name)
        unique.append((name, value))
    return unique


def extract_response_matrix(
    inputs: Any,
    n_cells: int,
    n_programs: int,
) -> tuple[np.ndarray, str]:
    candidates: list[tuple[int, str, np.ndarray]] = []
    preferred_tokens = ["y", "response", "outcome", "clr", "usage"]

    for name, value in iter_object_fields(inputs):
        try:
            if hasattr(value, "to_numpy"):
                array = value.to_numpy(dtype=float)
            else:
                array = np.asarray(value, dtype=float)
        except Exception:
            continue
        if array.ndim != 2:
            continue
        if array.shape == (n_cells, n_programs):
            normalized = array
        elif array.shape == (n_programs, n_cells):
            normalized = array.T
        else:
            continue
        lower_name = name.lower()
        score = 0
        for index, token in enumerate(preferred_tokens):
            if lower_name == token:
                score += 100 - index
            elif token in lower_name:
                score += 50 - index
        candidates.append((score, name, normalized))

    if not candidates:
        observed_shapes = []
        for name, value in iter_object_fields(inputs):
            try:
                observed_shapes.append((name, np.asarray(value).shape))
            except Exception:
                pass
        raise RuntimeError(
            "Could not locate the response matrix in prepare_crt_inputs output. "
            f"Expected {(n_cells, n_programs)}; observed {observed_shapes}."
        )

    candidates.sort(key=lambda item: item[0], reverse=True)
    _, field_name, matrix = candidates[0]
    matrix = np.asarray(matrix, dtype=float)
    if not np.isfinite(matrix).all():
        raise ValueError(
            f"Extracted response matrix {field_name!r} contains non-finite values."
        )
    return matrix, field_name


def rowwise_cosine(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    numerator = np.sum(x * y, axis=1)
    denominator = np.sqrt(np.sum(x * x, axis=1) * np.sum(y * y, axis=1))
    result = np.full(x.shape[0], np.nan, dtype=float)
    valid = denominator > 0
    result[valid] = numerator[valid] / denominator[valid]
    return result


def safe_flat_pearson(x: np.ndarray, y: np.ndarray) -> float:
    a = np.asarray(x, dtype=float).ravel()
    b = np.asarray(y, dtype=float).ravel()
    keep = np.isfinite(a) & np.isfinite(b)
    if keep.sum() < 3 or np.std(a[keep]) == 0 or np.std(b[keep]) == 0:
        return np.nan
    return float(np.corrcoef(a[keep], b[keep])[0, 1])


def safe_flat_spearman(
    x: np.ndarray,
    y: np.ndarray,
    max_entries: int = 300_000,
) -> float:
    a = np.asarray(x, dtype=float).ravel()
    b = np.asarray(y, dtype=float).ravel()
    keep = np.isfinite(a) & np.isfinite(b)
    a = a[keep]
    b = b[keep]
    if len(a) > max_entries:
        rng = np.random.default_rng(1)
        indices = rng.choice(len(a), size=max_entries, replace=False)
        a = a[indices]
        b = b[indices]
    if len(a) < 3:
        return np.nan
    a_rank = pd.Series(a).rank(method="average").to_numpy(dtype=float)
    b_rank = pd.Series(b).rank(method="average").to_numpy(dtype=float)
    if np.std(a_rank) == 0 or np.std(b_rank) == 0:
        return np.nan
    return float(np.corrcoef(a_rank, b_rank)[0, 1])


def column_response_metrics(
    current: np.ndarray,
    reference: np.ndarray,
    program_names: Sequence[str],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for j, program in enumerate(program_names):
        current_j = current[:, j]
        reference_j = reference[:, j]
        delta = current_j - reference_j
        rows.append(
            {
                "program_name": str(program),
                "response_pearson_vs_reference": safe_flat_pearson(
                    current_j,
                    reference_j,
                ),
                "response_spearman_vs_reference": safe_flat_spearman(
                    current_j,
                    reference_j,
                    max_entries=len(current_j),
                ),
                "response_rmse_vs_reference": float(
                    np.sqrt(np.mean(delta**2))
                ),
                "response_median_abs_change_vs_reference": float(
                    np.median(np.abs(delta))
                ),
                "response_p95_abs_change_vs_reference": float(
                    np.quantile(np.abs(delta), 0.95)
                ),
            }
        )
    return pd.DataFrame(rows)


def make_target_program_zero_burden(
    adata: Any,
    raw_usage: pd.DataFrame,
    ntc_labels: list[str],
) -> pd.DataFrame:
    guide_assignment = adata.obsm["guide_assignment"]
    guide_names = list(map(str, adata.uns["guide_names"]))
    guide2gene = {str(k): str(v) for k, v in adata.uns["guide2gene"].items()}

    if isinstance(guide_assignment, pd.DataFrame):
        guide_df = guide_assignment.copy()
        guide_df.index = adata.obs_names
        guide_df.columns = guide_df.columns.astype(str)
        guide_df = guide_df.reindex(columns=guide_names, fill_value=0)
        guide_values = guide_df.to_numpy(dtype=bool, copy=False)
    else:
        guide_values = np.asarray(guide_assignment, dtype=bool)
        if guide_values.shape[1] != len(guide_names):
            raise ValueError(
                "Guide assignment matrix column count does not match guide_names."
            )

    target_to_indices: dict[str, list[int]] = {}
    ntc_keys = {str(value).strip().lower() for value in ntc_labels}
    for index, guide in enumerate(guide_names):
        target = guide2gene.get(guide, guide.split("_", 1)[0])
        if str(target).strip().lower() in ntc_keys:
            continue
        target_to_indices.setdefault(str(target), []).append(index)

    zero_mask = raw_usage.to_numpy(dtype=float, copy=False) == 0
    programs = list(map(str, raw_usage.columns))
    rows: list[pd.DataFrame] = []

    for target, indices in target_to_indices.items():
        treated = guide_values[:, indices].any(axis=1)
        n_treated = int(treated.sum())
        n_control = int((~treated).sum())
        if n_treated == 0 or n_control == 0:
            continue
        treated_fraction = zero_mask[treated].mean(axis=0)
        control_fraction = zero_mask[~treated].mean(axis=0)
        rows.append(
            pd.DataFrame(
                {
                    "target_name": target,
                    "program_name": programs,
                    "n_treated_cells": n_treated,
                    "n_control_cells": n_control,
                    "zero_fraction_treated": treated_fraction,
                    "zero_fraction_control": control_fraction,
                    "zero_fraction_treated_minus_control": (
                        treated_fraction - control_fraction
                    ),
                }
            )
        )

    if not rows:
        return pd.DataFrame(
            columns=[
                "target_name",
                "program_name",
                "n_treated_cells",
                "n_control_cells",
                "zero_fraction_treated",
                "zero_fraction_control",
                "zero_fraction_treated_minus_control",
            ]
        )
    return pd.concat(rows, ignore_index=True)


def save_wide(frame: pd.DataFrame, path: Path) -> None:
    current = frame.copy()
    current.index = current.index.astype(str)
    current.columns = current.columns.astype(str)
    current.index.name = "target_name"
    path.parent.mkdir(parents=True, exist_ok=True)
    current.to_csv(path, sep="\t")


def flatten_ntc_frames(
    frames: Mapping[Any, pd.DataFrame] | None,
) -> np.ndarray:
    if not frames:
        return np.empty(0, dtype=float)
    arrays = [
        frame.to_numpy(dtype=float, copy=False).ravel()
        for frame in frames.values()
    ]
    return np.concatenate(arrays) if arrays else np.empty(0, dtype=float)


def call_empirical_real_crt(
    programde: SimpleNamespace,
    inputs: Any,
    B: int,
    n_jobs: int,
    real_crt_seed: int,
) -> tuple[dict[str, Any], str | None]:
    """Run the uncalibrated empirical CRT with the Bernoulli index sampler.

    Current programDE returns the empirical p-value matrix as ``pvals_df`` when
    skew-normal calibration is disabled.  Some older checkouts expose extra
    return flags, so they are set only when present and are kept false to avoid
    requesting skew-normal/raw-calibration outputs.
    """
    kwargs: dict[str, Any] = {
        "inputs": inputs,
        "B": B,
        "n_jobs": n_jobs,
        "calibrate_skew_normal": False,
    }
    signature = inspect.signature(programde.run_all_genes_union_crt)
    if "return_raw_pvals" in signature.parameters:
        kwargs["return_raw_pvals"] = False
    if "return_skew_normal" in signature.parameters:
        kwargs["return_skew_normal"] = False
    if "resampling_method" in signature.parameters:
        kwargs["resampling_method"] = "bernoulli_index"
    seed_parameter = None
    for candidate in ["base_seed", "seed0", "seed", "random_state"]:
        if candidate in signature.parameters:
            kwargs[candidate] = real_crt_seed
            seed_parameter = candidate
            break
    return programde.run_all_genes_union_crt(**kwargs), seed_parameter


def expected_run_files(run_dir: Path, day: str) -> list[Path]:
    return [
        run_dir / f"{day}_target_betas_wide.tsv",
        run_dir / f"{day}_target_pvalues_empirical_wide.tsv",
        run_dir / f"{day}_target_treated_counts.tsv",
        run_dir / f"{day}_ntc_pseudotarget_pvalues_empirical.npy",
        run_dir / f"{day}_run_config.json",
    ]


def matching_complete_run(
    run_dir: Path,
    day: str,
    fingerprint: str,
) -> bool:
    config_path = run_dir / f"{day}_run_config.json"
    if not all(path.exists() for path in expected_run_files(run_dir, day)):
        return False
    try:
        config = json.loads(config_path.read_text())
    except Exception:
        return False
    return config.get("config_fingerprint") == fingerprint


def count_ntc_groups(ntc_groups: Any) -> dict[str, int | None]:
    n_ensembles: int | None = None
    total_groups: int | None = None
    try:
        n_ensembles = len(ntc_groups)
    except Exception:
        pass
    try:
        if isinstance(ntc_groups, Mapping):
            total_groups = sum(len(value) for value in ntc_groups.values())
        else:
            total_groups = sum(len(value) for value in ntc_groups)
    except Exception:
        pass
    return {
        "n_ensembles": n_ensembles,
        "total_groups_across_ensembles": total_groups,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run two-stage cNMF CLR flooring sensitivity with the six-covariate "
            "iPSC-EC Bernoulli CRT."
        )
    )
    parser.add_argument("--programde-root", default=qqb.DEFAULT_PROGRAMDE_ROOT)
    parser.add_argument("--data-dir", default=qqb.DEFAULT_DATA_DIR)
    parser.add_argument("--cnmf-dir", default=qqb.DEFAULT_CNMF_DIR)
    parser.add_argument("--cnmf-file", default=None)
    parser.add_argument("--usage-file-template", default=None)
    parser.add_argument(
        "--guide-long-csv-template",
        default=qqb.DEFAULT_GUIDE_LONG_TEMPLATE,
    )
    parser.add_argument(
        "--guide-metadata-tsv",
        default="./300genes_guide_metadata_v43.tsv",
    )
    parser.add_argument("--guide-metadata-guide-col", default=None)
    parser.add_argument("--guide-metadata-target-col", default=None)
    parser.add_argument(
        "--out-dir",
        default="results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov",
    )
    parser.add_argument("--days", nargs="+", default=DEFAULT_DAYS)
    parser.add_argument("--k", type=int, default=50)
    parser.add_argument("--density-threshold", default="0_2")
    parser.add_argument("--gex-modality", default=None)
    parser.add_argument("--guide-modality", default=None)
    parser.add_argument(
        "--initial-floor-factors",
        nargs="+",
        type=float,
        default=[0.1, 1.0, 10.0],
    )
    parser.add_argument(
        "--internal-eps-quantiles",
        nargs="+",
        type=float,
        default=[0.0, 1e-5, 1e-4, 1e-3],
    )
    parser.add_argument(
        "--reference-initial-floor-factor",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "--reference-internal-eps-quantile",
        type=float,
        default=None,
        help=(
            "Current internal programDE flooring quantile. Omit to use the local "
            "prepare_crt_inputs default."
        ),
    )
    parser.add_argument(
        "--crt-grid-mode",
        choices=["primary", "full"],
        default="primary",
        help=(
            "primary: vary each flooring stage separately around the current "
            "reference; full: run the entire Cartesian grid."
        ),
    )
    parser.add_argument("--B", type=int, default=10000)
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")),
    )
    parser.add_argument("--real-crt-seed", type=int, default=123)
    parser.add_argument("--ntc-label", action="append", default=None)
    parser.add_argument("--ntc-group-size", type=int, default=6)
    parser.add_argument("--ntc-ensembles", type=int, default=10)
    parser.add_argument("--ntc-bins", type=int, default=10)
    parser.add_argument("--ntc-max-groups", type=int, default=None)
    parser.add_argument("--ntc-group-seed", type=int, default=7)
    parser.add_argument("--ntc-crt-seed", type=int, default=23)
    parser.add_argument("--require-covariates", action="store_true")
    parser.add_argument("--reuse-matching", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-crt", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    if args.ntc_label is None:
        args.ntc_label = DEFAULT_NTC_LABELS.copy()
    if args.B < 1:
        raise ValueError("--B must be >= 1")
    if args.n_jobs < 1:
        raise ValueError("--n-jobs must be >= 1")

    args.initial_floor_factors = unique_floats(args.initial_floor_factors)
    args.internal_eps_quantiles = unique_floats(args.internal_eps_quantiles)

    if not all(
        np.isfinite(value) and value > 0
        for value in args.initial_floor_factors
    ):
        raise ValueError("All initial-floor factors must be positive and finite.")
    if not all(
        np.isfinite(value) and 0.0 <= value <= 1.0
        for value in args.internal_eps_quantiles
    ):
        raise ValueError("All internal eps quantiles must lie in [0, 1].")
    if not np.isfinite(args.reference_initial_floor_factor) or (
        args.reference_initial_floor_factor <= 0
    ):
        raise ValueError("Reference initial-floor factor must be positive.")
    if args.reuse_matching and args.overwrite:
        raise ValueError("Choose either --reuse-matching or --overwrite, not both.")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    programde = qqb.load_programde(args.programde_root)
    has_eps_quantile, local_default_quantile = get_prepare_eps_quantile(programde)
    if not has_eps_quantile:
        raise RuntimeError(
            "The local programDE prepare_crt_inputs() does not expose eps_quantile."
        )

    reference_quantile = resolve_reference_quantile(
        args.reference_internal_eps_quantile,
        local_default_quantile,
    )
    reference_initial = float(args.reference_initial_floor_factor)

    initial_factors = unique_floats(
        [*args.initial_floor_factors, reference_initial]
    )
    internal_quantiles = unique_floats(
        [*args.internal_eps_quantiles, reference_quantile]
    )
    initial_factors.sort()
    internal_quantiles.sort()

    audit_combinations = build_audit_combinations(
        initial_factors,
        internal_quantiles,
        reference_initial,
        reference_quantile,
    )
    crt_combinations = build_crt_combinations(
        initial_factors,
        internal_quantiles,
        reference_initial,
        reference_quantile,
        args.crt_grid_mode,
    )

    density_token = qqb.normalize_density_token(args.density_threshold)
    model_root = Path(args.out_dir).expanduser() / qqb.model_tag(
        args.k,
        density_token,
    )
    model_root.mkdir(parents=True, exist_ok=True)

    programde_commit = git_commit(args.programde_root)
    programde_prepare_source = inspect.getsourcefile(programde.prepare_crt_inputs)
    programde_prepare_signature = str(
        inspect.signature(programde.prepare_crt_inputs)
    )

    plan_rows = []
    crt_set = {
        (float(initial), float(quantile))
        for initial, quantile in crt_combinations
    }
    for initial, quantile in audit_combinations:
        plan_rows.append(
            {
                "initial_floor_factor": initial,
                "initial_floor_label": initial_floor_label(initial),
                "internal_eps_quantile": quantile,
                "internal_quantile_label": internal_quantile_label(quantile),
                "combination_label": combination_label(initial, quantile),
                "is_reference": bool(
                    float_equal(initial, reference_initial)
                    and float_equal(quantile, reference_quantile)
                ),
                "run_transformation_audit": True,
                "run_empirical_crt": (initial, quantile) in crt_set
                and not args.skip_crt,
            }
        )
    plan = pd.DataFrame(plan_rows)
    plan_path = model_root / "two_stage_flooring_planned_combinations.csv"
    plan.to_csv(plan_path, index=False)

    print("iPSC-EC two-stage CLR floor sensitivity: six-covariate model")
    print("days:", args.days)
    print("initial-floor factors:", initial_factors)
    print("internal eps quantiles:", internal_quantiles)
    print("reference initial factor:", reference_initial)
    print("reference internal quantile:", reference_quantile)
    print("local programDE default quantile:", local_default_quantile)
    print("CRT grid mode:", args.crt_grid_mode)
    print("audit combinations:", len(audit_combinations))
    print("CRT combinations:", 0 if args.skip_crt else len(crt_combinations))
    print("B:", args.B, "n_jobs:", args.n_jobs)
    print("model root:", model_root)
    print("programDE commit:", programde_commit)
    print("planned combinations:", plan_path)

    helper_args = SimpleNamespace(
        cnmf_file=args.cnmf_file,
        usage_file_template=args.usage_file_template,
        cnmf_dir=args.cnmf_dir,
        k=args.k,
    )

    if args.dry_run:
        path_rows = []
        for day in args.days:
            path_rows.append(
                {
                    "day": day,
                    "h5mu": str(
                        qqb.h5mu_path_for_day(Path(args.data_dir), day)
                    ),
                    "guide_csv": str(
                        qqb.guide_long_path_for_day(
                            args.guide_long_csv_template,
                            day,
                        )
                    ),
                    "usage_candidates": ";".join(
                        str(path)
                        for path in qqb.usage_candidates_from_args(
                            helper_args,
                            day,
                            density_token,
                        )[:10]
                    ),
                }
            )
        path_table = pd.DataFrame(path_rows)
        path_table.to_csv(model_root / "dry_run_paths.csv", index=False)
        print(path_table.to_string(index=False))
        print("\nCombination plan:")
        print(plan.to_string(index=False))
        return 0

    import muon as mu

    run_records: list[RunRecord] = []
    audit_rows: list[dict[str, Any]] = []
    program_audit_rows: list[pd.DataFrame] = []
    zero_burden_paths: list[dict[str, str]] = []
    ntc_group_rows: list[dict[str, Any]] = []

    for day in args.days:
        print(f"\n{'=' * 100}\n{day}\n{'=' * 100}")
        day_started = time.time()

        h5mu_path = qqb.h5mu_path_for_day(
            Path(args.data_dir).expanduser(),
            day,
        )
        guide_csv = qqb.guide_long_path_for_day(
            args.guide_long_csv_template,
            day,
        )
        if not h5mu_path.exists():
            raise FileNotFoundError(f"Missing H5MU for {day}: {h5mu_path}")

        mdata = mu.read(h5mu_path)
        gex_modality = qqb.resolve_modality(
            mdata,
            args.gex_modality,
            ["GEX", "gex", "rna", "RNA"],
            "GEX/RNA",
        )
        guide_modality = qqb.resolve_modality(
            mdata,
            args.guide_modality,
            ["guide", "Guide", "gRNA", "grna"],
            "guide",
        )
        adata = mdata[gex_modality].copy()

        usage_candidates = qqb.usage_candidates_from_args(
            helper_args,
            day,
            density_token,
        )
        usage_path, matched_usage_cells = qqb.choose_usage_file(
            usage_candidates,
            pd.Index(adata.obs_names),
            day,
        )
        usage_raw, matched_usage_cells = align_raw_usage(
            usage_path,
            pd.Index(adata.obs_names),
        )

        guide_names, _ = qqb.load_and_attach_guides(
            adata,
            mdata,
            guide_modality,
            guide_csv,
            args.ntc_label,
            args.guide_metadata_tsv,
            args.guide_metadata_guide_col,
            args.guide_metadata_target_col,
        )
        covariates = qqb.attach_covariates(
            adata,
            require_covariates=args.require_covariates,
        )

        print(
            f"cells={adata.n_obs:,}; programs={usage_raw.shape[1]:,}; "
            f"guides={len(guide_names):,}; matched_usage_cells={matched_usage_cells:,}"
        )
        print("covariate columns:", list(covariates.columns))

        zero_burden = make_target_program_zero_burden(
            adata,
            usage_raw,
            args.ntc_label,
        )
        zero_burden.insert(0, "day", day)
        zero_burden_path = (
            model_root
            / "zero_burden"
            / f"{day}_target_program_zero_burden.tsv"
        )
        zero_burden_path.parent.mkdir(parents=True, exist_ok=True)
        zero_burden.to_csv(zero_burden_path, sep="\t", index=False)
        zero_burden_paths.append({"day": day, "path": str(zero_burden_path)})

        def prepare_variant(
            initial_factor: float,
            internal_quantile: float,
        ) -> dict[str, Any]:
            stage1_usage, stage1_meta = floor_and_normalize_raw_usage(
                usage_raw,
                initial_factor,
            )
            stage1_usage.index = adata.obs_names
            manual_clr, stage2_meta, stage2_program, cell_diag = (
                stage2_floor_metrics_and_manual_clr(
                    stage1_usage,
                    usage_raw,
                    internal_quantile,
                )
            )

            adata.obsm["cnmf_usage"] = stage1_usage
            adata.uns["program_names"] = list(map(str, stage1_usage.columns))
            inputs = prepare_inputs(
                programde,
                adata,
                eps_quantile=internal_quantile,
            )
            response, response_field = extract_response_matrix(
                inputs,
                n_cells=adata.n_obs,
                n_programs=stage1_usage.shape[1],
            )
            response = np.asarray(response, dtype=float).copy()

            manual_difference = np.abs(response - manual_clr)
            max_manual_difference = float(np.max(manual_difference))
            rmse_manual_difference = float(
                np.sqrt(np.mean(manual_difference**2))
            )
            if not np.allclose(
                response,
                manual_clr,
                rtol=1e-10,
                atol=1e-10,
                equal_nan=True,
            ):
                raise RuntimeError(
                    "The extracted programDE response does not match the manually "
                    "reconstructed two-stage CLR transform. Maximum absolute "
                    f"difference={max_manual_difference:.6g}. Inspect the local "
                    "clr_from_usage implementation before interpreting this audit."
                )

            combo = combination_label(initial_factor, internal_quantile)
            stage2_program = stage2_program.copy()
            stage2_program.insert(0, "combination_label", combo)
            stage2_program.insert(
                0,
                "internal_quantile_label",
                internal_quantile_label(internal_quantile),
            )
            stage2_program.insert(0, "internal_eps_quantile", internal_quantile)
            stage2_program.insert(
                0,
                "initial_floor_label",
                initial_floor_label(initial_factor),
            )
            stage2_program.insert(0, "initial_floor_factor", initial_factor)
            stage2_program.insert(0, "day", day)

            return {
                "inputs": inputs,
                "response": response,
                "response_field": response_field,
                "stage1_usage": stage1_usage,
                "stage1_meta": stage1_meta,
                "stage2_meta": stage2_meta,
                "stage2_program": stage2_program,
                "cell_diag": cell_diag,
                "manual_vs_programde_max_abs_difference": max_manual_difference,
                "manual_vs_programde_rmse": rmse_manual_difference,
            }

        reference = prepare_variant(reference_initial, reference_quantile)
        reference_response = reference["response"].copy()
        program_names = list(map(str, usage_raw.columns))

        ntc_guides, guide_freq, guide_to_bin, real_signatures = (
            programde.build_ntc_group_inputs(
                inputs=reference["inputs"],
                ntc_label=args.ntc_label,
                group_size=args.ntc_group_size,
                n_bins=args.ntc_bins,
            )
        )
        ntc_groups = programde.make_ntc_groups_ensemble(
            ntc_guides=ntc_guides,
            ntc_freq=guide_freq,
            real_gene_bin_sigs=real_signatures,
            guide_to_bin=guide_to_bin,
            n_ensemble=args.ntc_ensembles,
            seed0=args.ntc_group_seed,
            group_size=args.ntc_group_size,
            max_groups=args.ntc_max_groups,
        )
        ntc_group_summary = count_ntc_groups(ntc_groups)
        ntc_group_rows.append(
            {
                "day": day,
                "n_ntc_guides": int(len(ntc_guides)),
                "ntc_group_size": args.ntc_group_size,
                "ntc_bins": args.ntc_bins,
                "ntc_ensembles_requested": args.ntc_ensembles,
                "ntc_group_seed": args.ntc_group_seed,
                **ntc_group_summary,
            }
        )

        for initial_factor, internal_quantile in audit_combinations:
            is_reference = bool(
                float_equal(initial_factor, reference_initial)
                and float_equal(internal_quantile, reference_quantile)
            )
            current = (
                reference
                if is_reference
                else prepare_variant(initial_factor, internal_quantile)
            )

            response = current["response"]
            delta = response - reference_response
            absolute_delta = np.abs(delta)
            cell_rmse = np.sqrt(np.mean(delta**2, axis=1))
            cell_similarity = rowwise_cosine(response, reference_response)

            cell_diag = current["cell_diag"]
            audit_rows.append(
                {
                    "day": day,
                    "initial_floor_factor": float(initial_factor),
                    "initial_floor_label": initial_floor_label(initial_factor),
                    "internal_eps_quantile": float(internal_quantile),
                    "internal_quantile_label": internal_quantile_label(
                        internal_quantile
                    ),
                    "combination_label": combination_label(
                        initial_factor,
                        internal_quantile,
                    ),
                    "is_reference": is_reference,
                    "reference_initial_floor_factor": reference_initial,
                    "reference_internal_eps_quantile": reference_quantile,
                    **current["stage1_meta"],
                    **current["stage2_meta"],
                    "response_matrix_field": current["response_field"],
                    "manual_vs_programde_max_abs_difference": current[
                        "manual_vs_programde_max_abs_difference"
                    ],
                    "manual_vs_programde_rmse": current[
                        "manual_vs_programde_rmse"
                    ],
                    "response_entrywise_pearson_vs_reference": safe_flat_pearson(
                        response,
                        reference_response,
                    ),
                    "response_entrywise_spearman_vs_reference": safe_flat_spearman(
                        response,
                        reference_response,
                    ),
                    "response_rmse_vs_reference": float(
                        np.sqrt(np.mean(delta**2))
                    ),
                    "response_median_abs_change_vs_reference": float(
                        np.median(absolute_delta)
                    ),
                    "response_p95_abs_change_vs_reference": float(
                        np.quantile(absolute_delta, 0.95)
                    ),
                    "response_p99_abs_change_vs_reference": float(
                        np.quantile(absolute_delta, 0.99)
                    ),
                    "median_cell_response_cosine_vs_reference": float(
                        np.nanmedian(cell_similarity)
                    ),
                    "p05_cell_response_cosine_vs_reference": float(
                        np.nanquantile(cell_similarity, 0.05)
                    ),
                    "spearman_raw_zeros_per_cell_vs_cell_response_rmse": safe_flat_spearman(
                        cell_diag["raw_zeros_per_cell"],
                        cell_rmse,
                        max_entries=len(cell_rmse),
                    ),
                    "spearman_stage2_modified_entries_per_cell_vs_cell_response_rmse": safe_flat_spearman(
                        cell_diag["stage2_modified_entries_per_cell"],
                        cell_rmse,
                        max_entries=len(cell_rmse),
                    ),
                    "spearman_stage2_modified_positive_entries_per_cell_vs_cell_response_rmse": safe_flat_spearman(
                        cell_diag["stage2_modified_positive_entries_per_cell"],
                        cell_rmse,
                        max_entries=len(cell_rmse),
                    ),
                    "spearman_stage2_pseudomass_vs_cell_response_rmse": safe_flat_spearman(
                        cell_diag["stage2_pseudomass_fraction_per_cell"],
                        cell_rmse,
                        max_entries=len(cell_rmse),
                    ),
                    "usage_path": str(usage_path),
                    "h5mu_path": str(h5mu_path),
                }
            )

            program_metrics = column_response_metrics(
                response,
                reference_response,
                program_names,
            )
            program_metrics = current["stage2_program"].merge(
                program_metrics,
                on="program_name",
                how="left",
                validate="one_to_one",
            )
            program_metrics["is_reference"] = is_reference
            program_audit_rows.append(program_metrics)

            if not is_reference:
                del current
                gc.collect()

        if not args.skip_crt:
            for initial_factor, internal_quantile in crt_combinations:
                run_started = time.time()
                is_reference = bool(
                    float_equal(initial_factor, reference_initial)
                    and float_equal(internal_quantile, reference_quantile)
                )
                label = combination_label(initial_factor, internal_quantile)
                run_dir = model_root / "crt_runs" / label / day
                run_dir.mkdir(parents=True, exist_ok=True)

                current = (
                    reference
                    if is_reference
                    else prepare_variant(initial_factor, internal_quantile)
                )

                config_core = {
                    "analysis": "iPSC-EC two-stage CLR floor sensitivity",
                    "covariate_model": "six-source-covariate full-day model",
                    "resampling_method": "bernoulli_index",
                    "skew_normal_calibration": False,
                    "stratified_resampling": False,
                    "day": day,
                    "initial_floor_factor": initial_factor,
                    "initial_floor_label": initial_floor_label(initial_factor),
                    "internal_eps_quantile": internal_quantile,
                    "internal_quantile_label": internal_quantile_label(
                        internal_quantile
                    ),
                    "combination_label": label,
                    "is_reference": is_reference,
                    "reference_initial_floor_factor": reference_initial,
                    "reference_internal_eps_quantile": reference_quantile,
                    "local_programde_default_eps_quantile": local_default_quantile,
                    "B": args.B,
                    "n_jobs": args.n_jobs,
                    "real_crt_seed_requested": args.real_crt_seed,
                    "ntc_labels": args.ntc_label,
                    "ntc_group_size": args.ntc_group_size,
                    "ntc_ensembles": args.ntc_ensembles,
                    "ntc_bins": args.ntc_bins,
                    "ntc_max_groups": args.ntc_max_groups,
                    "ntc_group_seed": args.ntc_group_seed,
                    "ntc_crt_seed": args.ntc_crt_seed,
                    "stage1_floor_metadata": current["stage1_meta"],
                    "stage2_floor_metadata": current["stage2_meta"],
                    "manual_vs_programde_max_abs_difference": current[
                        "manual_vs_programde_max_abs_difference"
                    ],
                    "covariate_columns": list(covariates.columns),
                    "covariate_sources": adata.uns.get("covariate_sources", {}),
                    "h5mu_path": str(h5mu_path),
                    "usage_path": str(usage_path),
                    "guide_long_csv": str(guide_csv),
                    "guide_metadata_tsv": args.guide_metadata_tsv,
                    "n_cells": int(adata.n_obs),
                    "n_programs": int(usage_raw.shape[1]),
                    "n_guides": int(len(guide_names)),
                    "matched_usage_cells": int(matched_usage_cells),
                    "programde_root": str(args.programde_root),
                    "programde_git_commit": programde_commit,
                    "programde_prepare_crt_inputs_source": programde_prepare_source,
                    "programde_prepare_crt_inputs_signature": programde_prepare_signature,
                }
                fingerprint = config_fingerprint(config_core)

                if args.reuse_matching and matching_complete_run(
                    run_dir,
                    day,
                    fingerprint,
                ):
                    print(f"{day} / {label}: reusing config-matched complete run")
                    run_records.append(
                        RunRecord(
                            day=day,
                            initial_floor_factor=initial_factor,
                            internal_eps_quantile=internal_quantile,
                            initial_floor_label=initial_floor_label(initial_factor),
                            internal_quantile_label=internal_quantile_label(
                                internal_quantile
                            ),
                            combination_label=label,
                            is_reference=is_reference,
                            status="reused_matching",
                            output_dir=str(run_dir),
                            n_cells=int(adata.n_obs),
                            n_programs=int(usage_raw.shape[1]),
                            n_guides=int(len(guide_names)),
                            n_covariate_columns=int(covariates.shape[1]),
                            raw_usage_zero_fraction=current["stage1_meta"][
                                "fraction_exact_zero_raw_entries"
                            ],
                            base_raw_zero_floor=current["stage1_meta"][
                                "base_raw_zero_floor"
                            ],
                            applied_raw_zero_floor=current["stage1_meta"][
                                "applied_raw_zero_floor"
                            ],
                            effective_internal_floor=current["stage2_meta"][
                                "effective_internal_floor"
                            ],
                            stage2_fraction_entries_modified=current[
                                "stage2_meta"
                            ]["stage2_fraction_entries_modified"],
                            stage2_fraction_positive_entries_modified=current[
                                "stage2_meta"
                            ][
                                "stage2_fraction_originally_positive_entries_modified"
                            ],
                            stage2_fraction_zero_derived_entries_modified=current[
                                "stage2_meta"
                            ]["stage2_fraction_zero_derived_entries_modified"],
                            median_stage2_pseudomass_fraction=current[
                                "stage2_meta"
                            ]["median_stage2_pseudomass_fraction"],
                            config_fingerprint=fingerprint,
                        )
                    )
                    if not is_reference:
                        del current
                        gc.collect()
                    continue

                if not args.overwrite and any(
                    path.exists() for path in expected_run_files(run_dir, day)
                ):
                    raise FileExistsError(
                        f"Existing outputs found in {run_dir}. Use --overwrite or "
                        "--reuse-matching after checking the saved configuration."
                    )

                print(
                    f"{day} / {label}: initial floor="
                    f"{current['stage1_meta']['applied_raw_zero_floor']:.6g}; "
                    f"internal q={internal_quantile:.6g}; "
                    f"effective internal floor="
                    f"{current['stage2_meta']['effective_internal_floor']:.6g}; "
                    "running empirical Bernoulli CRT"
                )

                out, real_seed_parameter = call_empirical_real_crt(
                    programde,
                    current["inputs"],
                    B=args.B,
                    n_jobs=args.n_jobs,
                    real_crt_seed=args.real_crt_seed,
                )
                p_empirical = out.get("pvals_raw_df")
                if p_empirical is None:
                    p_empirical = out.get("pvals_df")
                if p_empirical is None:
                    raise RuntimeError(
                        f"{day} / {label}: CRT output lacks empirical p-values."
                    )
                betas = out.get("betas_df")
                treated = out.get("treated_df")
                if betas is None or treated is None:
                    raise RuntimeError(
                        f"{day} / {label}: CRT output lacks betas_df or treated_df."
                    )

                ntc_frames = programde.crt_pvals_for_ntc_groups_ensemble(
                    inputs=current["inputs"],
                    ntc_groups_ens=ntc_groups,
                    B=args.B,
                    seed0=args.ntc_crt_seed,
                )
                ntc_empirical = flatten_ntc_frames(ntc_frames)

                save_wide(
                    betas,
                    run_dir / f"{day}_target_betas_wide.tsv",
                )
                save_wide(
                    p_empirical,
                    run_dir / f"{day}_target_pvalues_empirical_wide.tsv",
                )
                treated_frame = pd.DataFrame(
                    {
                        "target_name": treated.index.astype(str),
                        "n_union_positive_cells": treated.to_numpy(),
                    }
                )
                treated_frame.to_csv(
                    run_dir / f"{day}_target_treated_counts.tsv",
                    sep="\t",
                    index=False,
                )
                np.save(
                    run_dir / f"{day}_ntc_pseudotarget_pvalues_empirical.npy",
                    ntc_empirical,
                )

                config = {
                    **config_core,
                    "config_fingerprint": fingerprint,
                    "real_crt_seed_parameter_used": real_seed_parameter,
                    "response_matrix_field": current["response_field"],
                    "runtime_seconds": round(time.time() - run_started, 3),
                    "output_dir": str(run_dir),
                }
                atomic_write_json(
                    run_dir / f"{day}_run_config.json",
                    config,
                )

                run_records.append(
                    RunRecord(
                        day=day,
                        initial_floor_factor=initial_factor,
                        internal_eps_quantile=internal_quantile,
                        initial_floor_label=initial_floor_label(initial_factor),
                        internal_quantile_label=internal_quantile_label(
                            internal_quantile
                        ),
                        combination_label=label,
                        is_reference=is_reference,
                        status="written",
                        output_dir=str(run_dir),
                        n_cells=int(adata.n_obs),
                        n_programs=int(usage_raw.shape[1]),
                        n_guides=int(len(guide_names)),
                        n_covariate_columns=int(covariates.shape[1]),
                        n_real_target_program_pairs=int(np.prod(betas.shape)),
                        n_ntc_pvalues_empirical=int(ntc_empirical.size),
                        raw_usage_zero_fraction=current["stage1_meta"][
                            "fraction_exact_zero_raw_entries"
                        ],
                        base_raw_zero_floor=current["stage1_meta"][
                            "base_raw_zero_floor"
                        ],
                        applied_raw_zero_floor=current["stage1_meta"][
                            "applied_raw_zero_floor"
                        ],
                        effective_internal_floor=current["stage2_meta"][
                            "effective_internal_floor"
                        ],
                        stage2_fraction_entries_modified=current["stage2_meta"][
                            "stage2_fraction_entries_modified"
                        ],
                        stage2_fraction_positive_entries_modified=current[
                            "stage2_meta"
                        ]["stage2_fraction_originally_positive_entries_modified"],
                        stage2_fraction_zero_derived_entries_modified=current[
                            "stage2_meta"
                        ]["stage2_fraction_zero_derived_entries_modified"],
                        median_stage2_pseudomass_fraction=current["stage2_meta"][
                            "median_stage2_pseudomass_fraction"
                        ],
                        runtime_seconds=round(time.time() - run_started, 3),
                        config_fingerprint=fingerprint,
                    )
                )

                if not is_reference:
                    del current
                    gc.collect()

        del reference
        del mdata, adata, usage_raw
        gc.collect()
        print(f"{day} complete in {(time.time() - day_started) / 60:.2f} minutes")

    audit = pd.DataFrame(audit_rows).sort_values(
        ["day", "initial_floor_factor", "internal_eps_quantile"]
    )
    audit_path = model_root / "two_stage_CLR_transformation_audit.csv"
    audit.to_csv(audit_path, index=False)

    program_audit = (
        pd.concat(program_audit_rows, ignore_index=True)
        if program_audit_rows
        else pd.DataFrame()
    )
    program_audit_path = model_root / "two_stage_CLR_program_audit.csv"
    program_audit.to_csv(program_audit_path, index=False)

    manifest = pd.DataFrame(asdict(record) for record in run_records)
    if not manifest.empty:
        manifest = manifest.sort_values(
            ["day", "initial_floor_factor", "internal_eps_quantile"]
        )
    manifest_path = model_root / "two_stage_floor_sensitivity_run_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    pd.DataFrame(zero_burden_paths).to_csv(
        model_root / "target_program_zero_burden_manifest.csv",
        index=False,
    )
    pd.DataFrame(ntc_group_rows).to_csv(
        model_root / "ntc_group_summary.csv",
        index=False,
    )

    overall_config = {
        "analysis": "iPSC-EC two-stage CLR floor sensitivity",
        "covariate_model": "six-source-covariate full-day model",
        "initial_floor_factors": initial_factors,
        "internal_eps_quantiles": internal_quantiles,
        "reference_initial_floor_factor": reference_initial,
        "reference_internal_eps_quantile": reference_quantile,
        "local_programde_default_eps_quantile": local_default_quantile,
        "crt_grid_mode": args.crt_grid_mode,
        "audit_combinations": audit_combinations,
        "crt_combinations": [] if args.skip_crt else crt_combinations,
        "B": args.B,
        "n_jobs": args.n_jobs,
        "days": args.days,
        "programde_root": args.programde_root,
        "programde_git_commit": programde_commit,
        "programde_prepare_crt_inputs_source": programde_prepare_source,
        "programde_prepare_crt_inputs_signature": programde_prepare_signature,
        "skew_normal_calibration": False,
        "stratified_resampling": False,
        "resampling_method": "bernoulli_index",
    }
    atomic_write_json(
        model_root / "two_stage_floor_sensitivity_overall_config.json",
        overall_config,
    )

    print("\nSaved:")
    print(" ", audit_path)
    print(" ", program_audit_path)
    print(" ", manifest_path)
    print(" ", model_root / "target_program_zero_burden_manifest.csv")
    print(" ", model_root / "ntc_group_summary.csv")
    print(" ", model_root / "two_stage_floor_sensitivity_overall_config.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

```

### Text output (stream)

```text
Overwriting iPSC_EC_CLR_two_stage_floor_sensitivity_6cov.py

```

## Markdown cell 7

## Syntax-check both backend scripts

## Code cell 8

```python
for script in [
    "iPSC_EC_QQ_export_all_days.py",
    "iPSC_EC_CLR_two_stage_floor_sensitivity_6cov.py",
]:
    subprocess.run([sys.executable, "-m", "py_compile", script], check=True)
    print("Syntax OK:", script)

```

### Text output (stream)

```text
Syntax OK: iPSC_EC_QQ_export_all_days.py
Syntax OK: iPSC_EC_CLR_two_stage_floor_sensitivity_6cov.py

```

## Markdown cell 9

## Build the dry-run or full-run command

## Code cell 10

```python
def build_two_stage_command(*, dry_run=False):
    cmd = [
        sys.executable,
        "iPSC_EC_CLR_two_stage_floor_sensitivity_6cov.py",
        "--programde-root", PROGRAMDE_ROOT,
        "--data-dir", DATA_DIR,
        "--cnmf-dir", CNMF_DIR,
        "--guide-long-csv-template", GUIDE_LONG_TEMPLATE,
        "--guide-metadata-tsv", GUIDE_METADATA_TSV,
        "--out-dir", OUT_DIR,
        "--days", *DAYS,
        "--k", str(K),
        "--density-threshold", DENSITY_THRESHOLD,
        "--initial-floor-factors", *[str(x) for x in INITIAL_FLOOR_FACTORS],
        "--internal-eps-quantiles", *[str(x) for x in INTERNAL_EPS_QUANTILES],
        "--reference-initial-floor-factor", str(REFERENCE_INITIAL_FLOOR_FACTOR),
        "--crt-grid-mode", CRT_GRID_MODE,
        "--B", str(B),
        "--n-jobs", str(N_JOBS),
        "--real-crt-seed", str(REAL_CRT_SEED),
        "--ntc-group-seed", str(NTC_GROUP_SEED),
        "--ntc-crt-seed", str(NTC_CRT_SEED),
    ]

    if REFERENCE_INTERNAL_EPS_QUANTILE is not None:
        cmd += [
            "--reference-internal-eps-quantile",
            str(REFERENCE_INTERNAL_EPS_QUANTILE),
        ]

    for label in NTC_LABELS:
        cmd += ["--ntc-label", label]

    if REQUIRE_COVARIATES:
        cmd.append("--require-covariates")
    if SKIP_CRT:
        cmd.append("--skip-crt")
    if OVERWRITE:
        cmd.append("--overwrite")
    elif REUSE_MATCHING:
        cmd.append("--reuse-matching")
    if dry_run:
        cmd.append("--dry-run")

    return cmd


def print_command(cmd):
    print(" \\\n  ".join(shlex.quote(str(part)) for part in cmd))

```

## Markdown cell 11

## Dry run: inspect paths and the planned flooring combinations

## Code cell 12

```python
dry_command = build_two_stage_command(dry_run=True)
print_command(dry_command)
subprocess.run(dry_command, check=True)

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python \
  iPSC_EC_CLR_two_stage_floor_sensitivity_6cov.py \
  --programde-root \
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE \
  --data-dir \
  /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data \
  --cnmf-dir \
  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all \
  --guide-long-csv-template \
  '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/{day}_sceptre_transfer/{day}_grna_assignments_long.csv' \
  --guide-metadata-tsv \
  ./300genes_guide_metadata_v43.tsv \
  --out-dir \
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov \
  --days \
  d0 \
  d1 \
  d2 \
  d3 \
  --k \
  50 \
  --density-threshold \
  0_2 \
  --initial-floor-factors \
  0.1 \
  1.0 \
  10.0 \
  --internal-eps-quantiles \
  0.0 \
  1e-05 \
  0.0001 \
  0.001 \
  --reference-initial-floor-factor \
  1.0 \
  --crt-grid-mode \
  primary \
  --B \
  10000 \
  --n-jobs \
  16 \
  --real-crt-seed \
  123 \
  --ntc-group-seed \
  7 \
  --ntc-crt-seed \
  23 \
  --ntc-label \
  SAFE \
  --ntc-label \
  non-targeting \
  --ntc-label \
  NTC \
  --require-covariates \
  --overwrite \
  --dry-run

```

### Text output (stream)

```text
iPSC-EC two-stage CLR floor sensitivity: six-covariate model
days: ['d0', 'd1', 'd2', 'd3']
initial-floor factors: [0.1, 1.0, 10.0]
internal eps quantiles: [0.0, 1e-05, 0.0001, 0.001]
reference initial factor: 1.0
reference internal quantile: 0.0001
local programDE default quantile: 0.0001
CRT grid mode: primary
audit combinations: 12
CRT combinations: 6
B: 10000 n_jobs: 16
model root: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2
programDE commit: None
planned combinations: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_flooring_planned_combinations.csv
day                                                                                                                                                      h5mu                                                                                                             guide_csv                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             usage_candidates
 d0 /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data/d0_normalized_gex_and_guide.h5mu /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d0_sceptre_transfer/d0_grna_assignments_long.csv /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt
 d1 /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data/d1_normalized_gex_and_guide.h5mu /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d1_sceptre_transfer/d1_grna_assignments_long.csv /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt
 d2 /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data/d2_normalized_gex_and_guide.h5mu /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/d2_grna_assignments_long.csv /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt
 d3 /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data/d3_normalized_gex_and_guide.h5mu /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d3_sceptre_transfer/d3_grna_assignments_long.csv /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt;/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt

Combination plan:
 initial_floor_factor initial_floor_label  internal_eps_quantile internal_quantile_label             combination_label  is_reference  run_transformation_audit  run_empirical_crt
                  1.0                  1x                0.00010                  q1em04   initial_1x__internal_q1em04          True                      True               True
                  0.1                0_1x                0.00000                      q0     initial_0_1x__internal_q0         False                      True              False
                  0.1                0_1x                0.00001                  q1em05 initial_0_1x__internal_q1em05         False                      True              False
                  0.1                0_1x                0.00010                  q1em04 initial_0_1x__internal_q1em04         False                      True               True
                  0.1                0_1x                0.00100                  q1em03 initial_0_1x__internal_q1em03         False                      True              False
                  1.0                  1x                0.00000                      q0       initial_1x__internal_q0         False                      True               True
                  1.0                  1x                0.00001                  q1em05   initial_1x__internal_q1em05         False                      True               True
                  1.0                  1x                0.00100                  q1em03   initial_1x__internal_q1em03         False                      True               True
                 10.0                 10x                0.00000                      q0      initial_10x__internal_q0         False                      True              False
                 10.0                 10x                0.00001                  q1em05  initial_10x__internal_q1em05         False                      True              False
                 10.0                 10x                0.00010                  q1em04  initial_10x__internal_q1em04         False                      True               True
                 10.0                 10x                0.00100                  q1em03  initial_10x__internal_q1em03         False                      True              False

```

### Text output (text/plain)

```text
CompletedProcess(args=['/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python', 'iPSC_EC_CLR_two_stage_floor_sensitivity_6cov.py', '--programde-root', '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE', '--data-dir', '/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data', '--cnmf-dir', '/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all', '--guide-long-csv-template', '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/{day}_sceptre_transfer/{day}_grna_assignments_long.csv', '--guide-metadata-tsv', './300genes_guide_metadata_v43.tsv', '--out-dir', 'results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov', '--days', 'd0', 'd1', 'd2', 'd3', '--k', '50', '--density-threshold', '0_2', '--initial-floor-factors', '0.1', '1.0', '10.0', '--internal-eps-quantiles', '0.0', '1e-05', '0.0001', '0.001', '--reference-initial-floor-factor', '1.0', '--crt-grid-mode', 'primary', '--B', '10000', '--n-jobs', '16', '--real-crt-seed', '123', '--ntc-group-seed', '7', '--ntc-crt-seed', '23', '--ntc-label', 'SAFE', '--ntc-label', 'non-targeting', '--ntc-label', 'NTC', '--require-covariates', '--overwrite', '--dry-run'], returncode=0)
```

## Markdown cell 13

## Run the full analysis

With the default `CRT_GRID_MODE = "primary"`, this performs 24 full empirical CRT runs: six flooring combinations for each of four days. Each run uses `B = 10,000`, so this cell can be computationally expensive. The full 12-condition factorial CRT requires 48 runs.

## Code cell 14

```python
run_command = build_two_stage_command(dry_run=False)
print_command(run_command)
subprocess.run(run_command, check=True)

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python \
  iPSC_EC_CLR_two_stage_floor_sensitivity_6cov.py \
  --programde-root \
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE \
  --data-dir \
  /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data \
  --cnmf-dir \
  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all \
  --guide-long-csv-template \
  '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/{day}_sceptre_transfer/{day}_grna_assignments_long.csv' \
  --guide-metadata-tsv \
  ./300genes_guide_metadata_v43.tsv \
  --out-dir \
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov \
  --days \
  d0 \
  d1 \
  d2 \
  d3 \
  --k \
  50 \
  --density-threshold \
  0_2 \
  --initial-floor-factors \
  0.1 \
  1.0 \
  10.0 \
  --internal-eps-quantiles \
  0.0 \
  1e-05 \
  0.0001 \
  0.001 \
  --reference-initial-floor-factor \
  1.0 \
  --crt-grid-mode \
  primary \
  --B \
  10000 \
  --n-jobs \
  16 \
  --real-crt-seed \
  123 \
  --ntc-group-seed \
  7 \
  --ntc-crt-seed \
  23 \
  --ntc-label \
  SAFE \
  --ntc-label \
  non-targeting \
  --ntc-label \
  NTC \
  --require-covariates \
  --overwrite
iPSC-EC two-stage CLR floor sensitivity: six-covariate model
days: ['d0', 'd1', 'd2', 'd3']
initial-floor factors: [0.1, 1.0, 10.0]
internal eps quantiles: [0.0, 1e-05, 0.0001, 0.001]
reference initial factor: 1.0
reference internal quantile: 0.0001
local programDE default quantile: 0.0001
CRT grid mode: primary
audit combinations: 12
CRT combinations: 6
B: 10000 n_jobs: 16
model root: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2
programDE commit: None
planned combinations: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_flooring_planned_combinations.csv

```

### Text output (stream)

```text
Matplotlib is building the font cache; this may take a moment.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/muon/_core/preproc.py:31: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
  if Version(scanpy.__version__) < Version("1.10"):

```

### Text output (stream)

```text

====================================================================================================
d0
====================================================================================================

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

### Text output (stream)

```text
usage candidates scored for d0
  overlap=38,609  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt
selected usage file for d0: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/iPSC_EC_QQ_export_all_days.py:306: UserWarning: Guide metadata file 300genes_guide_metadata_v43.tsv was not found. Falling back to guide_name.split('_', 1)[0].
  warnings.warn(

```

### Text output (stream)

```text
cells=38,609; programs=50; guides=1,926; matched_usage_cells=38,609
covariate columns: ['pct_counts_mt', 'log1p_total_counts', 'log1p_guides_per_cell', 'log1p_guide_umi_counts', 'log1p_n_genes_by_counts', 'biological_sample_d0_rep2']
d0 / initial_1x__internal_q1em04: initial floor=8.00193e-06; internal q=0.0001; effective internal floor=1.56359e-06; running empirical Bernoulli CRT
d0 / initial_1x__internal_q0: initial floor=8.00193e-06; internal q=0; effective internal floor=1.86384e-09; running empirical Bernoulli CRT
d0 / initial_1x__internal_q1em05: initial floor=8.00193e-06; internal q=1e-05; effective internal floor=1.04167e-07; running empirical Bernoulli CRT
d0 / initial_1x__internal_q1em03: initial floor=8.00193e-06; internal q=0.001; effective internal floor=6.00119e-06; running empirical Bernoulli CRT
d0 / initial_0_1x__internal_q1em04: initial floor=8.00193e-07; internal q=0.0001; effective internal floor=5.38116e-07; running empirical Bernoulli CRT
d0 / initial_10x__internal_q1em04: initial floor=8.00193e-05; internal q=0.0001; effective internal floor=1.56171e-06; running empirical Bernoulli CRT
d0 complete in 163.40 minutes

====================================================================================================
d1
====================================================================================================

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

### Text output (stream)

```text
usage candidates scored for d1
  overlap=16,925  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt
selected usage file for d1: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/iPSC_EC_QQ_export_all_days.py:306: UserWarning: Guide metadata file 300genes_guide_metadata_v43.tsv was not found. Falling back to guide_name.split('_', 1)[0].
  warnings.warn(

```

### Text output (stream)

```text
cells=16,925; programs=50; guides=1,926; matched_usage_cells=16,925
covariate columns: ['pct_counts_mt', 'log1p_total_counts', 'log1p_guides_per_cell', 'log1p_guide_umi_counts', 'log1p_n_genes_by_counts', 'biological_sample_d1_rep2']
d1 / initial_1x__internal_q1em04: initial floor=8.32582e-06; internal q=0.0001; effective internal floor=1.56642e-06; running empirical Bernoulli CRT
d1 / initial_1x__internal_q0: initial floor=8.32582e-06; internal q=0; effective internal floor=4.33192e-08; running empirical Bernoulli CRT
d1 / initial_1x__internal_q1em05: initial floor=8.32582e-06; internal q=1e-05; effective internal floor=2.23303e-07; running empirical Bernoulli CRT
d1 / initial_1x__internal_q1em03: initial floor=8.32582e-06; internal q=0.001; effective internal floor=6.54578e-06; running empirical Bernoulli CRT
d1 / initial_0_1x__internal_q1em04: initial floor=8.32582e-07; internal q=0.0001; effective internal floor=6.35829e-07; running empirical Bernoulli CRT
d1 / initial_10x__internal_q1em04: initial floor=8.32582e-05; internal q=0.0001; effective internal floor=1.56359e-06; running empirical Bernoulli CRT
d1 complete in 54.91 minutes

====================================================================================================
d2
====================================================================================================

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

### Text output (stream)

```text
usage candidates scored for d2
  overlap=21,980  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt
selected usage file for d2: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/iPSC_EC_QQ_export_all_days.py:306: UserWarning: Guide metadata file 300genes_guide_metadata_v43.tsv was not found. Falling back to guide_name.split('_', 1)[0].
  warnings.warn(

```

### Text output (stream)

```text
cells=21,980; programs=50; guides=1,926; matched_usage_cells=21,980
covariate columns: ['pct_counts_mt', 'log1p_total_counts', 'log1p_guides_per_cell', 'log1p_guide_umi_counts', 'log1p_n_genes_by_counts', 'biological_sample_d2_rep2']
d2 / initial_1x__internal_q1em04: initial floor=6.92259e-06; internal q=0.0001; effective internal floor=1.68952e-06; running empirical Bernoulli CRT
d2 / initial_1x__internal_q0: initial floor=6.92259e-06; internal q=0; effective internal floor=3.14284e-08; running empirical Bernoulli CRT
d2 / initial_1x__internal_q1em05: initial floor=6.92259e-06; internal q=1e-05; effective internal floor=1.62524e-07; running empirical Bernoulli CRT
d2 / initial_1x__internal_q1em03: initial floor=6.92259e-06; internal q=0.001; effective internal floor=4.96387e-06; running empirical Bernoulli CRT
d2 / initial_0_1x__internal_q1em04: initial floor=6.92259e-07; internal q=0.0001; effective internal floor=4.77165e-07; running empirical Bernoulli CRT
d2 / initial_10x__internal_q1em04: initial floor=6.92259e-05; internal q=0.0001; effective internal floor=1.68632e-06; running empirical Bernoulli CRT
d2 complete in 92.07 minutes

====================================================================================================
d3
====================================================================================================

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

### Text output (stream)

```text
usage candidates scored for d3
  overlap=44,102  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt
selected usage file for d3: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/iPSC_EC_QQ_export_all_days.py:306: UserWarning: Guide metadata file 300genes_guide_metadata_v43.tsv was not found. Falling back to guide_name.split('_', 1)[0].
  warnings.warn(

```

### Text output (stream)

```text
cells=44,102; programs=50; guides=1,926; matched_usage_cells=44,102
covariate columns: ['pct_counts_mt', 'log1p_total_counts', 'log1p_guides_per_cell', 'log1p_guide_umi_counts', 'log1p_n_genes_by_counts', 'biological_sample_d3_rep2']
d3 / initial_1x__internal_q1em04: initial floor=6.97082e-06; internal q=0.0001; effective internal floor=1.39913e-06; running empirical Bernoulli CRT
d3 / initial_1x__internal_q0: initial floor=6.97082e-06; internal q=0; effective internal floor=5.88065e-13; running empirical Bernoulli CRT
d3 / initial_1x__internal_q1em05: initial floor=6.97082e-06; internal q=1e-05; effective internal floor=1.88459e-07; running empirical Bernoulli CRT
d3 / initial_1x__internal_q1em03: initial floor=6.97082e-06; internal q=0.001; effective internal floor=4.98179e-06; running empirical Bernoulli CRT
d3 / initial_0_1x__internal_q1em04: initial floor=6.97082e-07; internal q=0.0001; effective internal floor=4.51189e-07; running empirical Bernoulli CRT
d3 / initial_10x__internal_q1em04: initial floor=6.97082e-05; internal q=0.0001; effective internal floor=1.39701e-06; running empirical Bernoulli CRT
d3 complete in 222.77 minutes

Saved:
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_CLR_transformation_audit.csv
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_CLR_program_audit.csv
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_floor_sensitivity_run_manifest.csv
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/target_program_zero_burden_manifest.csv
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/ntc_group_summary.csv
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_floor_sensitivity_overall_config.json

```

### Text output (text/plain)

```text
CompletedProcess(args=['/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python', 'iPSC_EC_CLR_two_stage_floor_sensitivity_6cov.py', '--programde-root', '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE', '--data-dir', '/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data', '--cnmf-dir', '/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all', '--guide-long-csv-template', '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/{day}_sceptre_transfer/{day}_grna_assignments_long.csv', '--guide-metadata-tsv', './300genes_guide_metadata_v43.tsv', '--out-dir', 'results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov', '--days', 'd0', 'd1', 'd2', 'd3', '--k', '50', '--density-threshold', '0_2', '--initial-floor-factors', '0.1', '1.0', '10.0', '--internal-eps-quantiles', '0.0', '1e-05', '0.0001', '0.001', '--reference-initial-floor-factor', '1.0', '--crt-grid-mode', 'primary', '--B', '10000', '--n-jobs', '16', '--real-crt-seed', '123', '--ntc-group-seed', '7', '--ntc-crt-seed', '23', '--ntc-label', 'SAFE', '--ntc-label', 'non-targeting', '--ntc-label', 'NTC', '--require-covariates', '--overwrite'], returncode=0)
```

## Markdown cell 15

## Validate exact transformations, the six-covariate design, and run configurations

## Code cell 16

```python
from pathlib import Path
import json
import numpy as np
import pandas as pd
from IPython.display import display

PLAN_PATH = MODEL_ROOT / "two_stage_flooring_planned_combinations.csv"
AUDIT_PATH = MODEL_ROOT / "two_stage_CLR_transformation_audit.csv"
PROGRAM_AUDIT_PATH = MODEL_ROOT / "two_stage_CLR_program_audit.csv"
MANIFEST_PATH = MODEL_ROOT / "two_stage_floor_sensitivity_run_manifest.csv"
OVERALL_CONFIG_PATH = MODEL_ROOT / "two_stage_floor_sensitivity_overall_config.json"

for required_path in [PLAN_PATH, AUDIT_PATH, PROGRAM_AUDIT_PATH, OVERALL_CONFIG_PATH]:
    if not required_path.exists():
        raise FileNotFoundError(f"Missing expected output: {required_path}")

plan = pd.read_csv(PLAN_PATH)
audit = pd.read_csv(AUDIT_PATH)
program_audit = pd.read_csv(PROGRAM_AUDIT_PATH)
overall_config = json.loads(OVERALL_CONFIG_PATH.read_text())
manifest = pd.read_csv(MANIFEST_PATH) if MANIFEST_PATH.exists() else pd.DataFrame()

REFERENCE_Q = float(overall_config["reference_internal_eps_quantile"])
REFERENCE_INITIAL = float(overall_config["reference_initial_floor_factor"])
LOCAL_PROGRAMDE_DEFAULT_Q = overall_config.get("local_programde_default_eps_quantile")

print("Reference initial-floor factor:", REFERENCE_INITIAL)
print("Reference internal eps quantile:", REFERENCE_Q)
print("Local programDE default eps quantile:", LOCAL_PROGRAMDE_DEFAULT_Q)
print("Empirical CRT B:", overall_config["B"])
print("CRT grid mode:", overall_config["crt_grid_mode"])
print("programDE commit:", overall_config.get("programde_git_commit"))

print("\nPlanned combinations:")
display(plan)

key_audit_columns = [
    "day",
    "initial_floor_factor",
    "internal_eps_quantile",
    "is_reference",
    "applied_raw_zero_floor",
    "effective_internal_floor",
    "fraction_exact_zero_raw_entries",
    "stage2_fraction_entries_modified",
    "stage2_fraction_originally_positive_entries_modified",
    "stage2_fraction_zero_derived_entries_modified",
    "median_initial_pseudomass_fraction",
    "median_stage2_pseudomass_fraction",
    "response_entrywise_spearman_vs_reference",
    "response_p95_abs_change_vs_reference",
    "manual_vs_programde_max_abs_difference",
]
print("\nTransformation audit:")
display(audit[key_audit_columns])

max_manual_difference = float(audit["manual_vs_programde_max_abs_difference"].max())
print("Maximum manual-vs-programDE CLR difference:", max_manual_difference)
if max_manual_difference > 1e-8:
    raise RuntimeError(
        "The stage-specific audit does not reproduce the exact programDE response matrix."
    )

q0 = audit[np.isclose(audit["internal_eps_quantile"], 0.0)].copy()
if not q0.empty:
    max_q0_modified = float(q0["stage2_fraction_entries_modified"].max())
    print("Maximum fraction modified by stage 2 when q=0:", max_q0_modified)
    if max_q0_modified > 1e-12:
        raise RuntimeError(
            "q=0 unexpectedly changed entries. Inspect the local CLR implementation."
        )

if not manifest.empty:
    print("\nCRT run manifest:")
    display(manifest)

    covariate_rows = []
    for _, row in manifest[manifest["is_reference"].astype(bool)].iterrows():
        config_path = Path(row["output_dir"]) / f"{row['day']}_run_config.json"
        config = json.loads(config_path.read_text())
        columns = config.get("covariate_columns", [])
        sources = config.get("covariate_sources", {})
        covariate_rows.append({
            "day": row["day"],
            "n_covariate_columns": len(columns),
            "covariate_columns": ", ".join(columns),
            "source_obs_keys": ", ".join(sorted({
                str(value.get("obs_key"))
                for value in sources.values()
                if isinstance(value, dict) and value.get("obs_key") is not None
            })),
            "B": config.get("B"),
            "skew_normal_calibration": config.get("skew_normal_calibration"),
            "resampling_method": config.get("resampling_method"),
        })

    covariate_check = pd.DataFrame(covariate_rows)
    print("\nReference-run covariate audit:")
    display(covariate_check)

    expected_source_keys = {
        "pct_counts_mt",
        "biological_sample",
        "total_counts",
        "guides_per_cell",
        "guide_umi_counts",
        "n_genes_by_counts",
    }
    for _, row in covariate_check.iterrows():
        observed = {x.strip() for x in row["source_obs_keys"].split(",") if x.strip()}
        missing = expected_source_keys.difference(observed)
        if missing:
            raise RuntimeError(f"{row['day']} is missing six-covariate sources: {sorted(missing)}")
else:
    print("\nNo CRT manifest was produced because SKIP_CRT=True.")

```

### Text output (stream)

```text
Reference initial-floor factor: 1.0
Reference internal eps quantile: 0.0001
Local programDE default eps quantile: 0.0001
Empirical CRT B: 10000
CRT grid mode: primary
programDE commit: None

Planned combinations:

```

### Text output (text/plain)

```text
    initial_floor_factor initial_floor_label  internal_eps_quantile  \
0                    1.0                  1x                0.00010   
1                    0.1                0_1x                0.00000   
2                    0.1                0_1x                0.00001   
3                    0.1                0_1x                0.00010   
4                    0.1                0_1x                0.00100   
5                    1.0                  1x                0.00000   
6                    1.0                  1x                0.00001   
7                    1.0                  1x                0.00100   
8                   10.0                 10x                0.00000   
9                   10.0                 10x                0.00001   
10                  10.0                 10x                0.00010   
11                  10.0                 10x                0.00100   

   internal_quantile_label              combination_label  is_reference  \
0                   q1em04    initial_1x__internal_q1em04          True   
1                       q0      initial_0_1x__internal_q0         False   
2                   q1em05  initial_0_1x__internal_q1em05         False   
3                   q1em04  initial_0_1x__internal_q1em04         False   
4                   q1em03  initial_0_1x__internal_q1em03         False   
5                       q0        initial_1x__internal_q0         False   
6                   q1em05    initial_1x__internal_q1em05         False   
7                   q1em03    initial_1x__internal_q1em03         False   
8                       q0       initial_10x__internal_q0         False   
9                   q1em05   initial_10x__internal_q1em05         False   
10                  q1em04   initial_10x__internal_q1em04         False   
11                  q1em03   initial_10x__internal_q1em03         False   

    run_transformation_audit  run_empirical_crt  
0                       True               True  
1                       True              False  
2                       True              False  
3                       True               True  
4                       True              False  
5                       True               True  
6                       True               True  
7                       True               True  
8                       True              False  
9                       True              False  
10                      True               True  
11                      True              False  
```

### Text output (stream)

```text

Transformation audit:

```

### Text output (text/plain)

```text
   day  initial_floor_factor  internal_eps_quantile  is_reference  \
0   d0                   0.1                0.00000         False   
1   d0                   0.1                0.00001         False   
2   d0                   0.1                0.00010         False   
3   d0                   0.1                0.00100         False   
4   d0                   1.0                0.00000         False   
5   d0                   1.0                0.00001         False   
6   d0                   1.0                0.00010          True   
7   d0                   1.0                0.00100         False   
8   d0                  10.0                0.00000         False   
9   d0                  10.0                0.00001         False   
10  d0                  10.0                0.00010         False   
11  d0                  10.0                0.00100         False   
12  d1                   0.1                0.00000         False   
13  d1                   0.1                0.00001         False   
14  d1                   0.1                0.00010         False   
15  d1                   0.1                0.00100         False   
16  d1                   1.0                0.00000         False   
17  d1                   1.0                0.00001         False   
18  d1                   1.0                0.00010          True   
19  d1                   1.0                0.00100         False   
20  d1                  10.0                0.00000         False   
21  d1                  10.0                0.00001         False   
22  d1                  10.0                0.00010         False   
23  d1                  10.0                0.00100         False   
24  d2                   0.1                0.00000         False   
25  d2                   0.1                0.00001         False   
26  d2                   0.1                0.00010         False   
27  d2                   0.1                0.00100         False   
28  d2                   1.0                0.00000         False   
29  d2                   1.0                0.00001         False   
30  d2                   1.0                0.00010          True   
31  d2                   1.0                0.00100         False   
32  d2                  10.0                0.00000         False   
33  d2                  10.0                0.00001         False   
34  d2                  10.0                0.00010         False   
35  d2                  10.0                0.00100         False   
36  d3                   0.1                0.00000         False   
37  d3                   0.1                0.00001         False   
38  d3                   0.1                0.00010         False   
39  d3                   0.1                0.00100         False   
40  d3                   1.0                0.00000         False   
41  d3                   1.0                0.00001         False   
42  d3                   1.0                0.00010          True   
43  d3                   1.0                0.00100         False   
44  d3                  10.0                0.00000         False   
45  d3                  10.0                0.00001         False   
46  d3                  10.0                0.00010         False   
47  d3                  10.0                0.00100         False   

    applied_raw_zero_floor  effective_internal_floor  \
0             8.001929e-07              1.864177e-09   
1             8.001929e-07              1.041811e-07   
2             8.001929e-07              5.381165e-07   
3             8.001929e-07              6.149873e-07   
4             8.001929e-06              1.863838e-09   
5             8.001929e-06              1.041671e-07   
6             8.001929e-06              1.563586e-06   
7             8.001929e-06              6.001187e-06   
8             8.001929e-05              1.860452e-09   
9             8.001929e-05              1.040267e-07   
10            8.001929e-05              1.561712e-06   
11            8.001929e-05              1.396906e-05   
12            8.325819e-07              4.332716e-08   
13            8.325819e-07              2.233497e-07   
14            8.325819e-07              6.358292e-07   
15            8.325819e-07              6.659473e-07   
16            8.325819e-06              4.331922e-08   
17            8.325819e-06              2.233028e-07   
18            8.325819e-06              1.566422e-06   
19            8.325819e-06              6.545783e-06   
20            8.325819e-05              4.324007e-08   
21            8.325819e-05              2.228349e-07   
22            8.325819e-05              1.563594e-06   
23            8.325819e-05              1.590146e-05   
24            6.922590e-07              3.143510e-08   
25            6.922590e-07              1.625514e-07   
26            6.922590e-07              4.771652e-07   
27            6.922590e-07              5.020718e-07   
28            6.922590e-06              3.142843e-08   
29            6.922590e-06              1.625242e-07   
30            6.922590e-06              1.689523e-06   
31            6.922590e-06              4.963874e-06   
32            6.922590e-05              3.136189e-08   
33            6.922590e-05              1.622528e-07   
34            6.922590e-05              1.686318e-06   
35            6.922590e-05              1.605626e-05   
36            6.970817e-07              5.881640e-13   
37            6.970817e-07              1.884983e-07   
38            6.970817e-07              4.511892e-07   
39            6.970817e-07              5.085351e-07   
40            6.970817e-06              5.880648e-13   
41            6.970817e-06              1.884589e-07   
42            6.970817e-06              1.399135e-06   
43            6.970817e-06              4.981785e-06   
44            6.970817e-05              5.870746e-13   
45            6.970817e-05              1.880778e-07   
46            6.970817e-05              1.397006e-06   
47            6.970817e-05              1.552487e-05   

    fraction_exact_zero_raw_entries  stage2_fraction_entries_modified  \
0                          0.446590                          0.000000   
1                          0.446590                          0.000010   
2                          0.446590                          0.000080   
3                          0.446590                          0.000999   
4                          0.446590                          0.000000   
5                          0.446590                          0.000010   
6                          0.446590                          0.000100   
7                          0.446590                          0.000995   
8                          0.446590                          0.000000   
9                          0.446590                          0.000010   
10                         0.446590                          0.000100   
11                         0.446590                          0.001000   
12                         0.507464                          0.000000   
13                         0.507464                          0.000011   
14                         0.507464                          0.000083   
15                         0.507464                          0.000980   
16                         0.507464                          0.000000   
17                         0.507464                          0.000011   
18                         0.507464                          0.000100   
19                         0.507464                          0.000981   
20                         0.507464                          0.000000   
21                         0.507464                          0.000011   
22                         0.507464                          0.000100   
23                         0.507464                          0.001001   
24                         0.565127                          0.000000   
25                         0.565127                          0.000010   
26                         0.565127                          0.000074   
27                         0.565127                          0.000971   
28                         0.565127                          0.000000   
29                         0.565127                          0.000010   
30                         0.565127                          0.000100   
31                         0.565127                          0.000979   
32                         0.565127                          0.000000   
33                         0.565127                          0.000010   
34                         0.565127                          0.000100   
35                         0.565127                          0.001000   
36                         0.548154                          0.000000   
37                         0.548154                          0.000010   
38                         0.548154                          0.000088   
39                         0.548154                          0.001000   
40                         0.548154                          0.000000   
41                         0.548154                          0.000010   
42                         0.548154                          0.000100   
43                         0.548154                          0.000989   
44                         0.548154                          0.000000   
45                         0.548154                          0.000010   
46                         0.548154                          0.000100   
47                         0.548154                          0.001000   

    stage2_fraction_originally_positive_entries_modified  \
0                                            0.000000      
1                                            0.000019      
2                                            0.000063      
3                                            0.000073      
4                                            0.000000      
5                                            0.000019      
6                                            0.000182      
7                                            0.000795      
8                                            0.000000      
9                                            0.000019      
10                                           0.000182      
11                                           0.001807      
12                                           0.000000      
13                                           0.000022      
14                                           0.000082      
15                                           0.000082      
16                                           0.000000      
17                                           0.000022      
18                                           0.000204      
19                                           0.000821      
20                                           0.000000      
21                                           0.000022      
22                                           0.000204      
23                                           0.002032      
24                                           0.000000      
25                                           0.000023      
26                                           0.000075      
27                                           0.000080      
28                                           0.000000      
29                                           0.000023      
30                                           0.000230      
31                                           0.000707      
32                                           0.000000      
33                                           0.000023      
34                                           0.000230      
35                                           0.002300      
36                                           0.000000      
37                                           0.000023      
38                                           0.000067      
39                                           0.000083      
40                                           0.000000      
41                                           0.000023      
42                                           0.000222      
43                                           0.000728      
44                                           0.000000      
45                                           0.000023      
46                                           0.000222      
47                                           0.002214      

    stage2_fraction_zero_derived_entries_modified  \
0                                        0.000000   
1                                        0.000000   
2                                        0.000101   
3                                        0.002147   
4                                        0.000000   
5                                        0.000000   
6                                        0.000000   
7                                        0.001243   
8                                        0.000000   
9                                        0.000000   
10                                       0.000000   
11                                       0.000000   
12                                       0.000000   
13                                       0.000000   
14                                       0.000084   
15                                       0.001851   
16                                       0.000000   
17                                       0.000000   
18                                       0.000000   
19                                       0.001136   
20                                       0.000000   
21                                       0.000000   
22                                       0.000000   
23                                       0.000000   
24                                       0.000000   
25                                       0.000000   
26                                       0.000072   
27                                       0.001657   
28                                       0.000000   
29                                       0.000000   
30                                       0.000000   
31                                       0.001188   
32                                       0.000000   
33                                       0.000000   
34                                       0.000000   
35                                       0.000000   
36                                       0.000000   
37                                       0.000000   
38                                       0.000104   
39                                       0.001756   
40                                       0.000000   
41                                       0.000000   
42                                       0.000000   
43                                       0.001205   
44                                       0.000000   
45                                       0.000000   
46                                       0.000000   
47                                       0.000000   

    median_initial_pseudomass_fraction  median_stage2_pseudomass_fraction  \
0                             0.000018                                0.0   
1                             0.000018                                0.0   
2                             0.000018                                0.0   
3                             0.000018                                0.0   
4                             0.000176                                0.0   
5                             0.000176                                0.0   
6                             0.000176                                0.0   
7                             0.000176                                0.0   
8                             0.001756                                0.0   
9                             0.001756                                0.0   
10                            0.001756                                0.0   
11                            0.001756                                0.0   
12                            0.000021                                0.0   
13                            0.000021                                0.0   
14                            0.000021                                0.0   
15                            0.000021                                0.0   
16                            0.000209                                0.0   
17                            0.000209                                0.0   
18                            0.000209                                0.0   
19                            0.000209                                0.0   
20                            0.002091                                0.0   
21                            0.002091                                0.0   
22                            0.002091                                0.0   
23                            0.002091                                0.0   
24                            0.000019                                0.0   
25                            0.000019                                0.0   
26                            0.000019                                0.0   
27                            0.000019                                0.0   
28                            0.000194                                0.0   
29                            0.000194                                0.0   
30                            0.000194                                0.0   
31                            0.000194                                0.0   
32                            0.001939                                0.0   
33                            0.001939                                0.0   
34                            0.001939                                0.0   
35                            0.001939                                0.0   
36                            0.000019                                0.0   
37                            0.000019                                0.0   
38                            0.000019                                0.0   
39                            0.000019                                0.0   
40                            0.000190                                0.0   
41                            0.000190                                0.0   
42                            0.000190                                0.0   
43                            0.000190                                0.0   
44                            0.001899                                0.0   
45                            0.001899                                0.0   
46                            0.001899                                0.0   
47                            0.001899                                0.0   

    response_entrywise_spearman_vs_reference  \
0                                   0.998457   
1                                   0.998457   
2                                   0.998455   
3                                   0.998454   
4                                   0.999999   
5                                   0.999999   
6                                   1.000000   
7                                   0.999964   
8                                   0.993063   
9                                   0.993064   
10                                  0.993068   
11                                  0.993089   
12                                  0.998622   
13                                  0.998621   
14                                  0.998616   
15                                  0.998615   
16                                  0.999998   
17                                  0.999998   
18                                  1.000000   
19                                  0.999955   
20                                  0.992403   
21                                  0.992403   
22                                  0.992410   
23                                  0.992445   
24                                  0.998672   
25                                  0.998672   
26                                  0.998667   
27                                  0.998667   
28                                  0.999996   
29                                  0.999996   
30                                  1.000000   
31                                  0.999965   
32                                  0.992077   
33                                  0.992078   
34                                  0.992088   
35                                  0.992145   
36                                  0.998670   
37                                  0.998669   
38                                  0.998667   
39                                  0.998666   
40                                  0.999990   
41                                  0.999991   
42                                  1.000000   
43                                  0.999960   
44                                  0.991792   
45                                  0.991796   
46                                  0.991804   
47                                  0.991859   

    response_p95_abs_change_vs_reference  \
0                               1.519706   
1                               1.519706   
2                               1.519706   
3                               1.519706   
4                               0.000000   
5                               0.000000   
6                               0.000000   
7                               0.000000   
8                               1.519706   
9                               1.519706   
10                              1.519706   
11                              1.519706   
12                              1.473654   
13                              1.473654   
14                              1.473654   
15                              1.473654   
16                              0.000000   
17                              0.000000   
18                              0.000000   
19                              0.000000   
20                              1.473654   
21                              1.473654   
22                              1.473654   
23                              1.473654   
24                              1.519706   
25                              1.519706   
26                              1.519706   
27                              1.519706   
28                              0.000000   
29                              0.000000   
30                              0.000000   
31                              0.000000   
32                              1.519706   
33                              1.519706   
34                              1.519706   
35                              1.519706   
36                              1.473654   
37                              1.473654   
38                              1.473654   
39                              1.473654   
40                              0.000000   
41                              0.000000   
42                              0.000000   
43                              0.000000   
44                              1.473654   
45                              1.473654   
46                              1.473654   
47                              1.473654   

    manual_vs_programde_max_abs_difference  
0                                      0.0  
1                                      0.0  
2                                      0.0  
3                                      0.0  
4                                      0.0  
5                                      0.0  
6                                      0.0  
7                                      0.0  
8                                      0.0  
9                                      0.0  
10                                     0.0  
11                                     0.0  
12                                     0.0  
13                                     0.0  
14                                     0.0  
15                                     0.0  
16                                     0.0  
17                                     0.0  
18                                     0.0  
19                                     0.0  
20                                     0.0  
21                                     0.0  
22                                     0.0  
23                                     0.0  
24                                     0.0  
25                                     0.0  
26                                     0.0  
27                                     0.0  
28                                     0.0  
29                                     0.0  
30                                     0.0  
31                                     0.0  
32                                     0.0  
33                                     0.0  
34                                     0.0  
35                                     0.0  
36                                     0.0  
37                                     0.0  
38                                     0.0  
39                                     0.0  
40                                     0.0  
41                                     0.0  
42                                     0.0  
43                                     0.0  
44                                     0.0  
45                                     0.0  
46                                     0.0  
47                                     0.0  
```

### Text output (stream)

```text
Maximum manual-vs-programDE CLR difference: 0.0
Maximum fraction modified by stage 2 when q=0: 0.0

CRT run manifest:

```

### Text output (text/plain)

```text
   day  initial_floor_factor  internal_eps_quantile initial_floor_label  \
0   d0                   0.1                0.00010                0_1x   
1   d0                   1.0                0.00000                  1x   
2   d0                   1.0                0.00001                  1x   
3   d0                   1.0                0.00010                  1x   
4   d0                   1.0                0.00100                  1x   
5   d0                  10.0                0.00010                 10x   
6   d1                   0.1                0.00010                0_1x   
7   d1                   1.0                0.00000                  1x   
8   d1                   1.0                0.00001                  1x   
9   d1                   1.0                0.00010                  1x   
10  d1                   1.0                0.00100                  1x   
11  d1                  10.0                0.00010                 10x   
12  d2                   0.1                0.00010                0_1x   
13  d2                   1.0                0.00000                  1x   
14  d2                   1.0                0.00001                  1x   
15  d2                   1.0                0.00010                  1x   
16  d2                   1.0                0.00100                  1x   
17  d2                  10.0                0.00010                 10x   
18  d3                   0.1                0.00010                0_1x   
19  d3                   1.0                0.00000                  1x   
20  d3                   1.0                0.00001                  1x   
21  d3                   1.0                0.00010                  1x   
22  d3                   1.0                0.00100                  1x   
23  d3                  10.0                0.00010                 10x   

   internal_quantile_label              combination_label  is_reference  \
0                   q1em04  initial_0_1x__internal_q1em04         False   
1                       q0        initial_1x__internal_q0         False   
2                   q1em05    initial_1x__internal_q1em05         False   
3                   q1em04    initial_1x__internal_q1em04          True   
4                   q1em03    initial_1x__internal_q1em03         False   
5                   q1em04   initial_10x__internal_q1em04         False   
6                   q1em04  initial_0_1x__internal_q1em04         False   
7                       q0        initial_1x__internal_q0         False   
8                   q1em05    initial_1x__internal_q1em05         False   
9                   q1em04    initial_1x__internal_q1em04          True   
10                  q1em03    initial_1x__internal_q1em03         False   
11                  q1em04   initial_10x__internal_q1em04         False   
12                  q1em04  initial_0_1x__internal_q1em04         False   
13                      q0        initial_1x__internal_q0         False   
14                  q1em05    initial_1x__internal_q1em05         False   
15                  q1em04    initial_1x__internal_q1em04          True   
16                  q1em03    initial_1x__internal_q1em03         False   
17                  q1em04   initial_10x__internal_q1em04         False   
18                  q1em04  initial_0_1x__internal_q1em04         False   
19                      q0        initial_1x__internal_q0         False   
20                  q1em05    initial_1x__internal_q1em05         False   
21                  q1em04    initial_1x__internal_q1em04          True   
22                  q1em03    initial_1x__internal_q1em03         False   
23                  q1em04   initial_10x__internal_q1em04         False   

     status                                         output_dir  n_cells  ...  \
0   written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    38609  ...   
1   written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    38609  ...   
2   written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    38609  ...   
3   written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    38609  ...   
4   written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    38609  ...   
5   written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    38609  ...   
6   written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    16925  ...   
7   written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    16925  ...   
8   written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    16925  ...   
9   written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    16925  ...   
10  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    16925  ...   
11  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    16925  ...   
12  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    21980  ...   
13  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    21980  ...   
14  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    21980  ...   
15  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    21980  ...   
16  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    21980  ...   
17  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    21980  ...   
18  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    44102  ...   
19  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    44102  ...   
20  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    44102  ...   
21  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    44102  ...   
22  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    44102  ...   
23  written  results/ipsc_ec_clr_two_stage_floor_sensitivit...    44102  ...   

    base_raw_zero_floor  applied_raw_zero_floor  effective_internal_floor  \
0              0.000008            8.001929e-07              5.381165e-07   
1              0.000008            8.001929e-06              1.863838e-09   
2              0.000008            8.001929e-06              1.041671e-07   
3              0.000008            8.001929e-06              1.563586e-06   
4              0.000008            8.001929e-06              6.001187e-06   
5              0.000008            8.001929e-05              1.561712e-06   
6              0.000008            8.325819e-07              6.358292e-07   
7              0.000008            8.325819e-06              4.331922e-08   
8              0.000008            8.325819e-06              2.233028e-07   
9              0.000008            8.325819e-06              1.566422e-06   
10             0.000008            8.325819e-06              6.545783e-06   
11             0.000008            8.325819e-05              1.563594e-06   
12             0.000007            6.922590e-07              4.771652e-07   
13             0.000007            6.922590e-06              3.142843e-08   
14             0.000007            6.922590e-06              1.625242e-07   
15             0.000007            6.922590e-06              1.689523e-06   
16             0.000007            6.922590e-06              4.963874e-06   
17             0.000007            6.922590e-05              1.686318e-06   
18             0.000007            6.970817e-07              4.511892e-07   
19             0.000007            6.970817e-06              5.880648e-13   
20             0.000007            6.970817e-06              1.884589e-07   
21             0.000007            6.970817e-06              1.399135e-06   
22             0.000007            6.970817e-06              4.981785e-06   
23             0.000007            6.970817e-05              1.397006e-06   

    stage2_fraction_entries_modified  \
0                           0.000080   
1                           0.000000   
2                           0.000010   
3                           0.000100   
4                           0.000995   
5                           0.000100   
6                           0.000083   
7                           0.000000   
8                           0.000011   
9                           0.000100   
10                          0.000981   
11                          0.000100   
12                          0.000074   
13                          0.000000   
14                          0.000010   
15                          0.000100   
16                          0.000979   
17                          0.000100   
18                          0.000088   
19                          0.000000   
20                          0.000010   
21                          0.000100   
22                          0.000989   
23                          0.000100   

    stage2_fraction_positive_entries_modified  \
0                                    0.000063   
1                                    0.000000   
2                                    0.000019   
3                                    0.000182   
4                                    0.000795   
5                                    0.000182   
6                                    0.000082   
7                                    0.000000   
8                                    0.000022   
9                                    0.000204   
10                                   0.000821   
11                                   0.000204   
12                                   0.000075   
13                                   0.000000   
14                                   0.000023   
15                                   0.000230   
16                                   0.000707   
17                                   0.000230   
18                                   0.000067   
19                                   0.000000   
20                                   0.000023   
21                                   0.000222   
22                                   0.000728   
23                                   0.000222   

    stage2_fraction_zero_derived_entries_modified  \
0                                        0.000101   
1                                        0.000000   
2                                        0.000000   
3                                        0.000000   
4                                        0.001243   
5                                        0.000000   
6                                        0.000084   
7                                        0.000000   
8                                        0.000000   
9                                        0.000000   
10                                       0.001136   
11                                       0.000000   
12                                       0.000072   
13                                       0.000000   
14                                       0.000000   
15                                       0.000000   
16                                       0.001188   
17                                       0.000000   
18                                       0.000104   
19                                       0.000000   
20                                       0.000000   
21                                       0.000000   
22                                       0.001205   
23                                       0.000000   

    median_stage2_pseudomass_fraction  runtime_seconds  \
0                                 0.0         1599.020   
1                                 0.0         1552.113   
2                                 0.0         1544.002   
3                                 0.0         1533.698   
4                                 0.0         1750.737   
5                                 0.0         1736.398   
6                                 0.0          537.502   
7                                 0.0          538.319   
8                                 0.0          549.701   
9                                 0.0          546.957   
10                                0.0          537.073   
11                                0.0          544.949   
12                                0.0          898.681   
13                                0.0          884.337   
14                                0.0          885.654   
15                                0.0          892.934   
16                                0.0          889.202   
17                                0.0         1032.403   
18                                0.0         2313.004   
19                                0.0         2378.145   
20                                0.0         2127.831   
21                                0.0         2126.671   
22                                0.0         2173.527   
23                                0.0         2157.918   

                                   config_fingerprint  message  
0   38cdc11ddbdd8ec65db3032cc7866ff272e3da44410e0a...      NaN  
1   3a9584f94bfbcd982ae94bba8b34cdd7fd270a8bea8f45...      NaN  
2   c97637ea034b2dd325dd47e3820a6dbbb9ff011670fb37...      NaN  
3   eaa0889e56e22e259ed301a2f6a20dbb1674551536884a...      NaN  
4   e1fe344404cfaefab932419a5083dff0d26e233713b2aa...      NaN  
5   5321c0193d0f98c47f6370cc41f6e738a9ec2d4cd8c92c...      NaN  
6   e864d4fffb246a775da19e5e34985ea221927ee04b397c...      NaN  
7   1fa5df00a23fdf50ca4845479e0880627f65e192182068...      NaN  
8   8f7e6bb82213157bbc5da7d11ed615c33d82a1bf88661b...      NaN  
9   3cfffd4c766985ee5106aed5dc0106af7315e1a4863d36...      NaN  
10  f10701d2910ba3348a32f9438b06ff09ce564cc7318b75...      NaN  
11  2fc0bfe85ee6f07106be3d733a5a180586584c5af94324...      NaN  
12  cdea97525d2cb7383a32bb0f33e4d64360ca179e22caca...      NaN  
13  f767d45c1afcaed50e53dc200ba2ff1cc2e28efb4dc0ee...      NaN  
14  1d5fcc0b33e5c10203b9b71caeab63507fadb1d638cfee...      NaN  
15  344a0e46cd70eb9d7baf15f40b13ccc5604006d33646d4...      NaN  
16  0aba296cf84b49e1afe5b2e8fba1de687db78d6e948da5...      NaN  
17  318bbe18be9929385d2e4840c9f8b21a090d430d8640e0...      NaN  
18  e4730d349613abece855fcc353185ce8dd15dc4b2bd766...      NaN  
19  03ebe4000c9b2e61b58f2a80e2d43b9e9adcb8811bd389...      NaN  
20  0b93d9607cff70222f912e1b311885a4cc84e64fe03386...      NaN  
21  077911e1d3183c7665d2e99a16351b6393cccad7063f78...      NaN  
22  3f5e54b2f14a457b02da4d11e58fe4a53edfd6d9dcf1b5...      NaN  
23  ebdaf7d9038361a0f6c35c9eeaafadad1ea4ccaf1fbb0e...      NaN  

[24 rows x 26 columns]
```

### Text output (stream)

```text

Reference-run covariate audit:

```

### Text output (text/plain)

```text
  day  n_covariate_columns                                  covariate_columns  \
0  d0                    6  pct_counts_mt, log1p_total_counts, log1p_guide...   
1  d1                    6  pct_counts_mt, log1p_total_counts, log1p_guide...   
2  d2                    6  pct_counts_mt, log1p_total_counts, log1p_guide...   
3  d3                    6  pct_counts_mt, log1p_total_counts, log1p_guide...   

                                     source_obs_keys      B  \
0  biological_sample, guide_umi_counts, guides_pe...  10000   
1  biological_sample, guide_umi_counts, guides_pe...  10000   
2  biological_sample, guide_umi_counts, guides_pe...  10000   
3  biological_sample, guide_umi_counts, guides_pe...  10000   

   skew_normal_calibration resampling_method  
0                    False   bernoulli_index  
1                    False   bernoulli_index  
2                    False   bernoulli_index  
3                    False   bernoulli_index  
```

## Markdown cell 17

## Plot Stage-1 and Stage-2 transformation behavior

These plots answer three separate questions:

- Does Stage 2 truncate genuinely positive usages, or only values created from raw zeros?
- How much compositional mass does Stage 2 add before renormalization?
- How much does the exact CLR response matrix change relative to the current implementation?

## Code cell 18

```python
from pathlib import Path
import matplotlib as mpl
import numpy as np
import pandas as pd
from IPython.display import display
from plotnine import (
    ggplot, aes,
    geom_tile, geom_text, geom_line, geom_point, geom_hline,
    facet_wrap,
    scale_fill_gradient, scale_color_manual,
    labs, theme_bw, theme,
    element_text, element_line, element_rect, element_blank,
)

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"

PLOT_DIR = MODEL_ROOT / "illustrator_editable_two_stage_floor_sensitivity"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

DAY_COLORS = {
    "d0": "#17365D",
    "d1": "#7A5195",
    "d2": "#D17A22",
    "d3": "#2A7F62",
}
GRID_COLOR = "#E8EDF3"
BORDER_COLOR = "#182635"


def save_vector_plot(plot, stem, width, height):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    plot.save(stem.with_suffix(".pdf"), width=width, height=height, units="in", verbose=False)
    plot.save(stem.with_suffix(".svg"), width=width, height=height, units="in", verbose=False)
    plot.save(stem.with_suffix(".png"), width=width, height=height, units="in", dpi=300, verbose=False)
    return {
        "pdf": stem.with_suffix(".pdf"),
        "svg": stem.with_suffix(".svg"),
        "png": stem.with_suffix(".png"),
    }


def initial_display(value):
    return f"{float(value):g}×"


def q_display(value, reference_q=REFERENCE_Q):
    value = float(value)
    if np.isclose(value, 0.0):
        label = "0\n(stage-1 only)"
    else:
        label = f"{value:.0e}"
    if np.isclose(value, reference_q):
        label += "\n(current)"
    return label


initial_levels = [initial_display(value) for value in sorted(audit["initial_floor_factor"].unique())]
q_levels = [q_display(value) for value in sorted(audit["internal_eps_quantile"].unique())]

heatmap_base = audit.copy()
heatmap_base["initial_display"] = heatmap_base["initial_floor_factor"].map(initial_display)
heatmap_base["q_display"] = heatmap_base["internal_eps_quantile"].map(q_display)
heatmap_base["initial_display"] = pd.Categorical(
    heatmap_base["initial_display"], categories=initial_levels, ordered=True
)
heatmap_base["q_display"] = pd.Categorical(
    heatmap_base["q_display"], categories=q_levels, ordered=True
)


def make_heatmap(value_column, title, subtitle, fill_label, label_formatter, stem_name):
    current = heatmap_base.copy()
    current["tile_label"] = current[value_column].map(label_formatter)
    plot = (
        ggplot(current, aes("initial_display", "q_display", fill=value_column))
        + geom_tile(color="white", size=0.6)
        + geom_text(aes(label="tile_label"), size=7.2)
        + facet_wrap("~day", ncol=2)
        + scale_fill_gradient(low="#F7FBFF", high="#17365D")
        + labs(
            title=title,
            subtitle=subtitle,
            x="Initial exact-zero floor factor",
            y="Internal programDE eps quantile",
            fill=fill_label,
        )
        + theme_bw(base_size=10)
        + theme(
            figure_size=(8.2, 6.6),
            plot_title=element_text(size=13, weight="bold", ha="center"),
            plot_subtitle=element_text(size=8.5, ha="center"),
            axis_text_x=element_text(size=9),
            axis_text_y=element_text(size=8),
            strip_text=element_text(size=9),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(color=BORDER_COLOR, size=0.65),
            legend_position="right",
        )
    )
    display(plot)
    saved = save_vector_plot(plot, PLOT_DIR / stem_name, width=8.2, height=6.6)
    return plot, saved


p_positive, positive_saved = make_heatmap(
    "stage2_fraction_originally_positive_entries_modified",
    "How often the second floor truncates genuinely positive cNMF usages",
    "Full 3×4 transformation grid; values are percentages of originally positive entries",
    "% positive entries modified",
    lambda value: f"{100 * value:.3f}%",
    "stage2_positive_entries_modified_heatmap",
)

p_pseudomass, pseudomass_saved = make_heatmap(
    "median_stage2_pseudomass_fraction",
    "Median compositional mass introduced by the second flooring step",
    "Mass is calculated before the post-floor row renormalization",
    "Median stage-2 mass",
    lambda value: f"{100 * value:.3f}%",
    "stage2_pseudomass_heatmap",
)

p_response, response_saved = make_heatmap(
    "response_entrywise_spearman_vs_reference",
    "Sensitivity of the exact CLR response matrix to both flooring choices",
    f"Reference: initial floor {REFERENCE_INITIAL:g}× and internal q={REFERENCE_Q:.0e}",
    "CLR rank correlation",
    lambda value: f"{value:.5f}",
    "CLR_response_spearman_heatmap",
)

# Axis-isolating plots: one panel set varies q at the current initial floor;
# the other varies the initial floor at the current q.
axis_rows = []
for _, row in audit.iterrows():
    if np.isclose(row["initial_floor_factor"], REFERENCE_INITIAL):
        current = row.to_dict()
        current["axis"] = "Vary internal floor; initial floor fixed at 1×"
        current["setting"] = q_display(row["internal_eps_quantile"])
        axis_rows.append(current)
    if np.isclose(row["internal_eps_quantile"], REFERENCE_Q):
        current = row.to_dict()
        current["axis"] = f"Vary initial floor; internal q fixed at {REFERENCE_Q:.0e}"
        current["setting"] = initial_display(row["initial_floor_factor"])
        axis_rows.append(current)

axis_data = pd.DataFrame(axis_rows)
metrics = {
    "response_entrywise_spearman_vs_reference": "CLR rank correlation vs current",
    "response_p95_abs_change_vs_reference": "95th percentile |CLR change|",
    "stage2_fraction_originally_positive_entries_modified": "Fraction of positive usages truncated",
    "median_stage2_pseudomass_fraction": "Median stage-2 pseudomass fraction",
}
axis_long = axis_data.melt(
    id_vars=["day", "axis", "setting"],
    value_vars=list(metrics),
    var_name="metric",
    value_name="value",
)
axis_long["metric_label"] = axis_long["metric"].map(metrics)
axis_long["panel"] = axis_long["axis"] + "\n" + axis_long["metric_label"]

p_axis = (
    ggplot(axis_long, aes("setting", "value", color="day", group="day"))
    + geom_hline(yintercept=0, color="#D5DCE3", size=0.3)
    + geom_line(size=0.75)
    + geom_point(size=2.0)
    + facet_wrap("~panel", scales="free", ncol=2)
    + scale_color_manual(values=DAY_COLORS)
    + labs(
        title="Axis-isolating transformation sensitivity",
        subtitle="One flooring stage is varied while the other is held at the current value",
        x="Flooring setting",
        y="Metric value",
        color="Day",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(9.2, 10.2),
        plot_title=element_text(size=13, weight="bold", ha="center"),
        plot_subtitle=element_text(size=8.5, ha="center"),
        axis_text_x=element_text(size=8, rotation=20, ha="right"),
        strip_text=element_text(size=8),
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
        legend_position="top",
    )
)
display(p_axis)
axis_saved = save_vector_plot(
    p_axis,
    PLOT_DIR / "axis_isolating_transformation_sensitivity",
    width=9.2,
    height=10.2,
)

print("Saved transformation plots under:", PLOT_DIR)
for name, paths in {
    "positive-entry modification": positive_saved,
    "stage-2 pseudomass": pseudomass_saved,
    "CLR response correlation": response_saved,
    "axis-isolating summary": axis_saved,
}.items():
    print(f"\n{name}")
    for fmt, path in paths.items():
        print(f"  {fmt.upper()}: {path}")

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c5475b350>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c463fcc10>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c4509d410>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c445f1a50>
```

### Text output (stream)

```text
Saved transformation plots under: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity

positive-entry modification
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/stage2_positive_entries_modified_heatmap.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/stage2_positive_entries_modified_heatmap.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/stage2_positive_entries_modified_heatmap.png

stage-2 pseudomass
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/stage2_pseudomass_heatmap.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/stage2_pseudomass_heatmap.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/stage2_pseudomass_heatmap.png

CLR response correlation
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/CLR_response_spearman_heatmap.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/CLR_response_spearman_heatmap.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/CLR_response_spearman_heatmap.png

axis-isolating summary
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/axis_isolating_transformation_sensitivity.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/axis_isolating_transformation_sensitivity.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/axis_isolating_transformation_sensitivity.png

```

## Markdown cell 19

## Compare empirical CRT betas, p-values, discoveries, and NTC calibration

## Code cell 20

```python
from pathlib import Path
import numpy as np
import pandas as pd
from IPython.display import display

if "manifest" not in globals():
    manifest = pd.read_csv(MODEL_ROOT / "two_stage_floor_sensitivity_run_manifest.csv")
if manifest.empty:
    raise RuntimeError("No CRT outputs are available. Run the full backend with SKIP_CRT=False.")
if "audit" not in globals():
    audit = pd.read_csv(MODEL_ROOT / "two_stage_CLR_transformation_audit.csv")
if "program_audit" not in globals():
    program_audit = pd.read_csv(MODEL_ROOT / "two_stage_CLR_program_audit.csv")
program_audit = program_audit.copy()
program_audit["day"] = program_audit["day"].astype(str)
program_audit["combination_label"] = program_audit["combination_label"].astype(str)
program_audit["program_name"] = program_audit["program_name"].astype(str)

NTC_KEYS = {str(value).strip().lower() for value in NTC_LABELS}
QVALUE_CUTOFFS = [0.05, 0.10]
TOP_N_REFERENCE_EFFECTS = 500


def as_bool(series):
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def bh_qvalues(values):
    p = np.asarray(values, dtype=float).ravel()
    q = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return q
    pv = np.clip(p[finite], 0.0, 1.0)
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q_finite = np.empty_like(pv)
    q_finite[order] = np.clip(ranked, 0.0, 1.0)
    q[finite] = q_finite
    return q


def safe_pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 3 or np.std(x[keep]) == 0 or np.std(y[keep]) == 0:
        return np.nan
    return float(np.corrcoef(x[keep], y[keep])[0, 1])


def safe_spearman(x, y):
    x = pd.Series(np.asarray(x, dtype=float))
    y = pd.Series(np.asarray(y, dtype=float))
    keep = np.isfinite(x) & np.isfinite(y)
    if int(keep.sum()) < 3:
        return np.nan
    x = x[keep]
    y = y[keep]
    if x.nunique() < 2 or y.nunique() < 2:
        return np.nan
    return float(x.corr(y, method="spearman"))


def concordance_correlation(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    if len(x) < 3:
        return np.nan
    covariance = np.mean((x - x.mean()) * (y - y.mean()))
    denominator = x.var() + y.var() + (x.mean() - y.mean()) ** 2
    if denominator == 0:
        return np.nan
    return float(2.0 * covariance / denominator)


def regression_slope_intercept(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    if len(x) < 3 or np.std(x) == 0:
        return np.nan, np.nan
    slope, intercept = np.polyfit(x, y, deg=1)
    return float(slope), float(intercept)


def jaccard(a, b):
    a = np.asarray(a, dtype=bool)
    b = np.asarray(b, dtype=bool)
    union = int((a | b).sum())
    return 1.0 if union == 0 else float((a & b).sum() / union)


def sign_concordance(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y) & (x != 0) & (y != 0)
    if not keep.any():
        return np.nan
    return float((np.sign(x[keep]) == np.sign(y[keep])).mean())


def ntc_calibration_metrics(pvalues):
    p = np.asarray(pvalues, dtype=float).ravel()
    p = p[np.isfinite(p)]
    p = np.clip(p, 1e-300, 1.0)
    if len(p) == 0:
        return {
            "n_ntc_pvalues": 0,
            "ntc_fraction_p_lt_0_01": np.nan,
            "ntc_fraction_p_lt_0_05": np.nan,
            "ntc_fraction_p_lt_0_10": np.nan,
            "ntc_median_p": np.nan,
            "ntc_qq_rmse": np.nan,
            "ntc_qq_tail_rmse": np.nan,
            "ntc_ks_distance_uniform": np.nan,
            "ntc_nominal_absolute_error_sum": np.nan,
        }

    p_sorted = np.sort(p)
    n = len(p_sorted)
    expected_p = (np.arange(1, n + 1) - 0.5) / n
    observed_log = -np.log10(p_sorted)
    expected_log = -np.log10(expected_p)
    error = observed_log - expected_log
    tail_n = max(10, int(np.ceil(0.01 * n)))

    i = np.arange(1, n + 1)
    d_plus = np.max(i / n - p_sorted)
    d_minus = np.max(p_sorted - (i - 1) / n)

    f01 = float((p < 0.01).mean())
    f05 = float((p < 0.05).mean())
    f10 = float((p < 0.10).mean())

    return {
        "n_ntc_pvalues": int(n),
        "ntc_fraction_p_lt_0_01": f01,
        "ntc_fraction_p_lt_0_05": f05,
        "ntc_fraction_p_lt_0_10": f10,
        "ntc_median_p": float(np.median(p)),
        "ntc_qq_rmse": float(np.sqrt(np.mean(error ** 2))),
        "ntc_qq_tail_rmse": float(np.sqrt(np.mean(error[:tail_n] ** 2))),
        "ntc_ks_distance_uniform": float(max(d_plus, d_minus)),
        "ntc_nominal_absolute_error_sum": float(
            abs(f01 - 0.01) + abs(f05 - 0.05) + abs(f10 - 0.10)
        ),
    }


def read_wide_as_series(path, value_name):
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    keep = ~frame.index.to_series().str.strip().str.lower().isin(NTC_KEYS).to_numpy()
    frame = frame.loc[keep]
    try:
        series = frame.stack(future_stack=True)
    except TypeError:
        series = frame.stack(dropna=False)
    series.index = series.index.set_names(["target_name", "program_name"])
    series = pd.to_numeric(series, errors="coerce")
    series.name = value_name
    return series


manifest = manifest.copy()
manifest["is_reference"] = as_bool(manifest["is_reference"])
manifest = manifest[manifest["status"].isin(["written", "reused_matching"])].copy()

summary_rows = []
effect_rows = []
scaling_rows = []

for day in DAYS:
    day_manifest = manifest[manifest["day"].astype(str).eq(day)].copy()
    if day_manifest.empty:
        print(f"{day}: no CRT runs found; skipping")
        continue

    reference_rows = day_manifest[day_manifest["is_reference"]]
    if len(reference_rows) != 1:
        raise RuntimeError(f"{day}: expected exactly one reference run, found {len(reference_rows)}")

    reference_row = reference_rows.iloc[0]
    reference_dir = Path(reference_row["output_dir"])
    beta_reference = read_wide_as_series(
        reference_dir / f"{day}_target_betas_wide.tsv",
        "beta_reference",
    )
    p_reference = read_wide_as_series(
        reference_dir / f"{day}_target_pvalues_empirical_wide.tsv",
        "pvalue_reference",
    )
    reference_table = pd.concat([beta_reference, p_reference], axis=1, join="inner").dropna()
    reference_table["qvalue_reference"] = bh_qvalues(reference_table["pvalue_reference"])

    ntc_reference = np.load(
        reference_dir / f"{day}_ntc_pseudotarget_pvalues_empirical.npy"
    )

    top_n = min(TOP_N_REFERENCE_EFFECTS, len(reference_table))
    top_reference_index = reference_table["beta_reference"].abs().nlargest(top_n).index

    zero_burden_path = MODEL_ROOT / "zero_burden" / f"{day}_target_program_zero_burden.tsv"
    zero_burden = pd.read_csv(zero_burden_path, sep="\t")
    zero_burden["target_name"] = zero_burden["target_name"].astype(str)
    zero_burden["program_name"] = zero_burden["program_name"].astype(str)

    for _, current_row in day_manifest.iterrows():
        current_dir = Path(current_row["output_dir"])
        combination = str(current_row["combination_label"])

        beta_current = read_wide_as_series(
            current_dir / f"{day}_target_betas_wide.tsv",
            "beta_current",
        )
        p_current = read_wide_as_series(
            current_dir / f"{day}_target_pvalues_empirical_wide.tsv",
            "pvalue_current",
        )

        comparison = pd.concat(
            [
                reference_table,
                beta_current,
                p_current,
            ],
            axis=1,
            join="inner",
        ).dropna()
        comparison["qvalue_current"] = bh_qvalues(comparison["pvalue_current"])
        comparison["beta_delta"] = comparison["beta_current"] - comparison["beta_reference"]
        comparison["abs_beta_delta"] = comparison["beta_delta"].abs()
        comparison["abs_beta_reference"] = comparison["beta_reference"].abs()
        comparison["day"] = day
        comparison["initial_floor_factor"] = float(current_row["initial_floor_factor"])
        comparison["internal_eps_quantile"] = float(current_row["internal_eps_quantile"])
        comparison["combination_label"] = combination
        comparison["is_reference"] = bool(current_row["is_reference"])
        comparison = comparison.reset_index()
        effect_rows.append(comparison)

        beta_ref_values = comparison["beta_reference"].to_numpy(dtype=float)
        beta_cur_values = comparison["beta_current"].to_numpy(dtype=float)
        beta_delta = beta_cur_values - beta_ref_values
        slope, intercept = regression_slope_intercept(beta_ref_values, beta_cur_values)
        reference_sd = float(np.std(beta_ref_values))

        top = comparison.set_index(["target_name", "program_name"]).reindex(top_reference_index)
        top = top.dropna(subset=["beta_reference", "beta_current"])

        p_ref_log = -np.log10(np.clip(comparison["pvalue_reference"], 1e-300, 1.0))
        p_cur_log = -np.log10(np.clip(comparison["pvalue_current"], 1e-300, 1.0))
        p_log_delta = p_cur_log.to_numpy() - p_ref_log.to_numpy()

        ntc_current = np.load(
            current_dir / f"{day}_ntc_pseudotarget_pvalues_empirical.npy"
        )
        ntc_metrics = ntc_calibration_metrics(ntc_current)
        ntc_pair_spearman = np.nan
        if np.asarray(ntc_current).size == np.asarray(ntc_reference).size:
            ntc_pair_spearman = safe_spearman(
                np.asarray(ntc_reference).ravel(),
                np.asarray(ntc_current).ravel(),
            )

        row = {
            "day": day,
            "initial_floor_factor": float(current_row["initial_floor_factor"]),
            "initial_floor_label": current_row["initial_floor_label"],
            "internal_eps_quantile": float(current_row["internal_eps_quantile"]),
            "internal_quantile_label": current_row["internal_quantile_label"],
            "combination_label": combination,
            "is_reference": bool(current_row["is_reference"]),
            "n_target_program_pairs": int(len(comparison)),
            "beta_pearson_vs_reference": safe_pearson(beta_ref_values, beta_cur_values),
            "beta_spearman_vs_reference": safe_spearman(beta_ref_values, beta_cur_values),
            "beta_concordance_correlation_vs_reference": concordance_correlation(
                beta_ref_values, beta_cur_values
            ),
            "beta_regression_slope_vs_reference": slope,
            "beta_regression_intercept_vs_reference": intercept,
            "beta_rmse_vs_reference": float(np.sqrt(np.mean(beta_delta ** 2))),
            "beta_normalized_rmse_vs_reference_sd": float(
                np.sqrt(np.mean(beta_delta ** 2)) / reference_sd
            ) if reference_sd > 0 else np.nan,
            "beta_median_abs_change_vs_reference": float(np.median(np.abs(beta_delta))),
            "beta_p95_abs_change_vs_reference": float(np.quantile(np.abs(beta_delta), 0.95)),
            "beta_max_abs_change_vs_reference": float(np.max(np.abs(beta_delta))),
            "beta_sign_concordance_all_nonzero": sign_concordance(
                beta_ref_values, beta_cur_values
            ),
            "top_reference_effect_sign_concordance": sign_concordance(
                top["beta_reference"], top["beta_current"]
            ),
            "pvalue_neglog10_spearman_vs_reference": safe_spearman(
                p_ref_log, p_cur_log
            ),
            "pvalue_neglog10_median_abs_change_vs_reference": float(
                np.median(np.abs(p_log_delta))
            ),
            "pvalue_neglog10_p95_abs_change_vs_reference": float(
                np.quantile(np.abs(p_log_delta), 0.95)
            ),
            "fraction_tests_more_significant_than_reference": float(
                (p_log_delta > 0).mean()
            ),
            "ntc_pvalue_pair_spearman_vs_reference": ntc_pair_spearman,
            **ntc_metrics,
            "beta_path": str(current_dir / f"{day}_target_betas_wide.tsv"),
            "pvalue_path": str(current_dir / f"{day}_target_pvalues_empirical_wide.tsv"),
            "ntc_path": str(current_dir / f"{day}_ntc_pseudotarget_pvalues_empirical.npy"),
        }

        for cutoff in QVALUE_CUTOFFS:
            suffix = str(cutoff).replace(".", "_")
            ref_discovery = comparison["qvalue_reference"].to_numpy() <= cutoff
            current_discovery = comparison["qvalue_current"].to_numpy() <= cutoff
            shared = ref_discovery & current_discovery
            row[f"n_reference_discoveries_q_le_{suffix}"] = int(ref_discovery.sum())
            row[f"n_current_discoveries_q_le_{suffix}"] = int(current_discovery.sum())
            row[f"n_shared_discoveries_q_le_{suffix}"] = int(shared.sum())
            row[f"n_reference_only_discoveries_q_le_{suffix}"] = int(
                (ref_discovery & ~current_discovery).sum()
            )
            row[f"n_current_only_discoveries_q_le_{suffix}"] = int(
                (~ref_discovery & current_discovery).sum()
            )
            row[f"discovery_jaccard_q_le_{suffix}"] = jaccard(
                ref_discovery, current_discovery
            )
            row[f"shared_discovery_sign_concordance_q_le_{suffix}"] = (
                sign_concordance(
                    beta_ref_values[shared],
                    beta_cur_values[shared],
                ) if shared.any() else np.nan
            )

        summary_rows.append(row)

        merged_zero = comparison.merge(
            zero_burden,
            on=["target_name", "program_name"],
            how="inner",
            validate="one_to_one",
        )

        current_program_audit = program_audit[
            (program_audit["day"].astype(str) == day)
            & (program_audit["combination_label"].astype(str) == combination)
        ].copy()
        program_delta = (
            comparison.groupby("program_name", as_index=False)
            .agg(
                mean_abs_beta_change=("abs_beta_delta", "mean"),
                median_abs_beta_change=("abs_beta_delta", "median"),
                p95_abs_beta_change=("abs_beta_delta", lambda x: x.quantile(0.95)),
            )
            .merge(
                current_program_audit,
                on="program_name",
                how="left",
                validate="one_to_one",
            )
        )

        scaling_rows.append({
            "day": day,
            "initial_floor_factor": float(current_row["initial_floor_factor"]),
            "internal_eps_quantile": float(current_row["internal_eps_quantile"]),
            "combination_label": combination,
            "is_reference": bool(current_row["is_reference"]),
            "n_pairs_with_zero_burden": int(len(merged_zero)),
            "spearman_abs_beta_change_vs_zero_fraction_treated": safe_spearman(
                merged_zero["abs_beta_delta"],
                merged_zero["zero_fraction_treated"],
            ),
            "spearman_abs_beta_change_vs_abs_zero_fraction_difference": safe_spearman(
                merged_zero["abs_beta_delta"],
                merged_zero["zero_fraction_treated_minus_control"].abs(),
            ),
            "spearman_signed_beta_change_vs_zero_fraction_difference": safe_spearman(
                merged_zero["beta_delta"],
                merged_zero["zero_fraction_treated_minus_control"],
            ),
            "spearman_program_mean_abs_beta_change_vs_stage2_positive_fraction_modified": safe_spearman(
                program_delta["mean_abs_beta_change"],
                program_delta["stage2_fraction_originally_positive_entries_modified"],
            ),
            "spearman_program_mean_abs_beta_change_vs_program_response_rmse": safe_spearman(
                program_delta["mean_abs_beta_change"],
                program_delta["response_rmse_vs_reference"],
            ),
        })

inference_summary = pd.DataFrame(summary_rows).sort_values(
    ["day", "initial_floor_factor", "internal_eps_quantile"]
)
effect_comparisons = pd.concat(effect_rows, ignore_index=True)
scaling_summary = pd.DataFrame(scaling_rows).sort_values(
    ["day", "initial_floor_factor", "internal_eps_quantile"]
)

INFERENCE_SUMMARY_PATH = MODEL_ROOT / "two_stage_floor_inference_summary.csv"
EFFECT_COMPARISON_PATH = MODEL_ROOT / "two_stage_floor_effect_level_comparison.csv"
SCALING_SUMMARY_PATH = MODEL_ROOT / "two_stage_floor_zero_and_program_scaling_summary.csv"

inference_summary.to_csv(INFERENCE_SUMMARY_PATH, index=False)
effect_comparisons.to_csv(EFFECT_COMPARISON_PATH, index=False)
scaling_summary.to_csv(SCALING_SUMMARY_PATH, index=False)

print("Inferential sensitivity summary:")
display(inference_summary)

print("\nDoes sensitivity scale with raw zeros or stage-2 modification?")
display(scaling_summary)

print("\nSaved:")
print(" ", INFERENCE_SUMMARY_PATH)
print(" ", EFFECT_COMPARISON_PATH)
print(" ", SCALING_SUMMARY_PATH)

```

### Text output (stream)

```text
Inferential sensitivity summary:

```

### Text output (text/plain)

```text
   day  initial_floor_factor initial_floor_label  internal_eps_quantile  \
0   d0                   0.1                0_1x                0.00010   
1   d0                   1.0                  1x                0.00000   
2   d0                   1.0                  1x                0.00001   
3   d0                   1.0                  1x                0.00010   
4   d0                   1.0                  1x                0.00100   
5   d0                  10.0                 10x                0.00010   
6   d1                   0.1                0_1x                0.00010   
7   d1                   1.0                  1x                0.00000   
8   d1                   1.0                  1x                0.00001   
9   d1                   1.0                  1x                0.00010   
10  d1                   1.0                  1x                0.00100   
11  d1                  10.0                 10x                0.00010   
12  d2                   0.1                0_1x                0.00010   
13  d2                   1.0                  1x                0.00000   
14  d2                   1.0                  1x                0.00001   
15  d2                   1.0                  1x                0.00010   
16  d2                   1.0                  1x                0.00100   
17  d2                  10.0                 10x                0.00010   
18  d3                   0.1                0_1x                0.00010   
19  d3                   1.0                  1x                0.00000   
20  d3                   1.0                  1x                0.00001   
21  d3                   1.0                  1x                0.00010   
22  d3                   1.0                  1x                0.00100   
23  d3                  10.0                 10x                0.00010   

   internal_quantile_label              combination_label  is_reference  \
0                   q1em04  initial_0_1x__internal_q1em04         False   
1                       q0        initial_1x__internal_q0         False   
2                   q1em05    initial_1x__internal_q1em05         False   
3                   q1em04    initial_1x__internal_q1em04          True   
4                   q1em03    initial_1x__internal_q1em03         False   
5                   q1em04   initial_10x__internal_q1em04         False   
6                   q1em04  initial_0_1x__internal_q1em04         False   
7                       q0        initial_1x__internal_q0         False   
8                   q1em05    initial_1x__internal_q1em05         False   
9                   q1em04    initial_1x__internal_q1em04          True   
10                  q1em03    initial_1x__internal_q1em03         False   
11                  q1em04   initial_10x__internal_q1em04         False   
12                  q1em04  initial_0_1x__internal_q1em04         False   
13                      q0        initial_1x__internal_q0         False   
14                  q1em05    initial_1x__internal_q1em05         False   
15                  q1em04    initial_1x__internal_q1em04          True   
16                  q1em03    initial_1x__internal_q1em03         False   
17                  q1em04   initial_10x__internal_q1em04         False   
18                  q1em04  initial_0_1x__internal_q1em04         False   
19                      q0        initial_1x__internal_q0         False   
20                  q1em05    initial_1x__internal_q1em05         False   
21                  q1em04    initial_1x__internal_q1em04          True   
22                  q1em03    initial_1x__internal_q1em03         False   
23                  q1em04   initial_10x__internal_q1em04         False   

    n_target_program_pairs  beta_pearson_vs_reference  \
0                    14900                   0.997879   
1                    14900                   0.999994   
2                    14900                   0.999996   
3                    14900                   1.000000   
4                    14900                   0.999993   
5                    14900                   0.993319   
6                    14900                   0.997401   
7                    14900                   0.999996   
8                    14900                   0.999997   
9                    14900                   1.000000   
10                   14900                   0.999988   
11                   14900                   0.991322   
12                   14900                   0.997913   
13                   14900                   0.999993   
14                   14900                   0.999995   
15                   14900                   1.000000   
16                   14900                   0.999996   
17                   14900                   0.993742   
18                   14900                   0.998501   
19                   14900                   0.999997   
20                   14900                   0.999999   
21                   14900                   1.000000   
22                   14900                   0.999997   
23                   14900                   0.995988   

    beta_spearman_vs_reference  ...  n_current_only_discoveries_q_le_0_05  \
0                     0.997481  ...                                    20   
1                     0.999990  ...                                    22   
2                     0.999992  ...                                    27   
3                     1.000000  ...                                     0   
4                     0.999987  ...                                    21   
5                     0.991274  ...                                    73   
6                     0.997185  ...                                    18   
7                     0.999993  ...                                    20   
8                     0.999994  ...                                    27   
9                     1.000000  ...                                     0   
10                    0.999981  ...                                    27   
11                    0.990165  ...                                    25   
12                    0.997651  ...                                    11   
13                    0.999987  ...                                    13   
14                    0.999990  ...                                    17   
15                    1.000000  ...                                     0   
16                    0.999992  ...                                    16   
17                    0.991493  ...                                    50   
18                    0.997919  ...                                    25   
19                    0.999990  ...                                    18   
20                    0.999996  ...                                    19   
21                    1.000000  ...                                     0   
22                    0.999991  ...                                    24   
23                    0.993104  ...                                    82   

    discovery_jaccard_q_le_0_05  shared_discovery_sign_concordance_q_le_0_05  \
0                      0.895137                                          1.0   
1                      0.933333                                          1.0   
2                      0.929323                                          1.0   
3                      1.000000                                          1.0   
4                      0.937785                                          1.0   
5                      0.867792                                          1.0   
6                      0.771429                                          1.0   
7                      0.811321                                          1.0   
8                      0.840183                                          1.0   
9                      1.000000                                          1.0   
10                     0.835616                                          1.0   
11                     0.797235                                          1.0   
12                     0.873656                                          1.0   
13                     0.925134                                          1.0   
14                     0.907407                                          1.0   
15                     1.000000                                          1.0   
16                     0.917772                                          1.0   
17                     0.846715                                          1.0   
18                     0.919802                                          1.0   
19                     0.952144                                          1.0   
20                     0.954183                                          1.0   
21                     1.000000                                          1.0   
22                     0.951437                                          1.0   
23                     0.886598                                          1.0   

    n_reference_discoveries_q_le_0_1  n_current_discoveries_q_le_0_1  \
0                                890                             852   
1                                890                             898   
2                                890                             909   
3                                890                             890   
4                                890                             897   
5                                890                             961   
6                                276                             262   
7                                276                             278   
8                                276                             277   
9                                276                             276   
10                               276                             279   
11                               276                             298   
12                               509                             464   
13                               509                             522   
14                               509                             512   
15                               509                             509   
16                               509                             493   
17                               509                             529   
18                              1312                            1259   
19                              1312                            1303   
20                              1312                            1301   
21                              1312                            1312   
22                              1312                            1308   
23                              1312                            1361   

    n_shared_discoveries_q_le_0_1  n_reference_only_discoveries_q_le_0_1  \
0                             830                                     60   
1                             869                                     21   
2                             872                                     18   
3                             890                                      0   
4                             866                                     24   
5                             863                                     27   
6                             254                                     22   
7                             264                                     12   
8                             267                                      9   
9                             276                                      0   
10                            268                                      8   
11                            260                                     16   
12                            452                                     57   
13                            496                                     13   
14                            488                                     21   
15                            509                                      0   
16                            479                                     30   
17                            477                                     32   
18                           1238                                     74   
19                           1286                                     26   
20                           1280                                     32   
21                           1312                                      0   
22                           1285                                     27   
23                           1260                                     52   

    n_current_only_discoveries_q_le_0_1  discovery_jaccard_q_le_0_1  \
0                                    22                    0.910088   
1                                    29                    0.945593   
2                                    37                    0.940669   
3                                     0                    1.000000   
4                                    31                    0.940282   
5                                    98                    0.873482   
6                                     8                    0.894366   
7                                    14                    0.910345   
8                                    10                    0.933566   
9                                     0                    1.000000   
10                                   11                    0.933798   
11                                   38                    0.828025   
12                                   12                    0.867562   
13                                   26                    0.927103   
14                                   24                    0.915572   
15                                    0                    1.000000   
16                                   14                    0.915870   
17                                   52                    0.850267   
18                                   21                    0.928732   
19                                   17                    0.967645   
20                                   21                    0.960240   
21                                    0                    1.000000   
22                                   23                    0.962547   
23                                  101                    0.891720   

    shared_discovery_sign_concordance_q_le_0_1  
0                                          1.0  
1                                          1.0  
2                                          1.0  
3                                          1.0  
4                                          1.0  
5                                          1.0  
6                                          1.0  
7                                          1.0  
8                                          1.0  
9                                          1.0  
10                                         1.0  
11                                         1.0  
12                                         1.0  
13                                         1.0  
14                                         1.0  
15                                         1.0  
16                                         1.0  
17                                         1.0  
18                                         1.0  
19                                         1.0  
20                                         1.0  
21                                         1.0  
22                                         1.0  
23                                         1.0  

[24 rows x 51 columns]
```

### Text output (stream)

```text

Does sensitivity scale with raw zeros or stage-2 modification?

```

### Text output (text/plain)

```text
   day  initial_floor_factor  internal_eps_quantile  \
0   d0                   0.1                0.00010   
1   d0                   1.0                0.00000   
2   d0                   1.0                0.00001   
3   d0                   1.0                0.00010   
4   d0                   1.0                0.00100   
5   d0                  10.0                0.00010   
6   d1                   0.1                0.00010   
7   d1                   1.0                0.00000   
8   d1                   1.0                0.00001   
9   d1                   1.0                0.00010   
10  d1                   1.0                0.00100   
11  d1                  10.0                0.00010   
12  d2                   0.1                0.00010   
13  d2                   1.0                0.00000   
14  d2                   1.0                0.00001   
15  d2                   1.0                0.00010   
16  d2                   1.0                0.00100   
17  d2                  10.0                0.00010   
18  d3                   0.1                0.00010   
19  d3                   1.0                0.00000   
20  d3                   1.0                0.00001   
21  d3                   1.0                0.00010   
22  d3                   1.0                0.00100   
23  d3                  10.0                0.00010   

                combination_label  is_reference  n_pairs_with_zero_burden  \
0   initial_0_1x__internal_q1em04         False                     14900   
1         initial_1x__internal_q0         False                     14900   
2     initial_1x__internal_q1em05         False                     14900   
3     initial_1x__internal_q1em04          True                     14900   
4     initial_1x__internal_q1em03         False                     14900   
5    initial_10x__internal_q1em04         False                     14900   
6   initial_0_1x__internal_q1em04         False                     14900   
7         initial_1x__internal_q0         False                     14900   
8     initial_1x__internal_q1em05         False                     14900   
9     initial_1x__internal_q1em04          True                     14900   
10    initial_1x__internal_q1em03         False                     14900   
11   initial_10x__internal_q1em04         False                     14900   
12  initial_0_1x__internal_q1em04         False                     14900   
13        initial_1x__internal_q0         False                     14900   
14    initial_1x__internal_q1em05         False                     14900   
15    initial_1x__internal_q1em04          True                     14900   
16    initial_1x__internal_q1em03         False                     14900   
17   initial_10x__internal_q1em04         False                     14900   
18  initial_0_1x__internal_q1em04         False                     14900   
19        initial_1x__internal_q0         False                     14900   
20    initial_1x__internal_q1em05         False                     14900   
21    initial_1x__internal_q1em04          True                     14900   
22    initial_1x__internal_q1em03         False                     14900   
23   initial_10x__internal_q1em04         False                     14900   

    spearman_abs_beta_change_vs_zero_fraction_treated  \
0                                            0.079068   
1                                            0.164578   
2                                            0.173769   
3                                                 NaN   
4                                            0.373217   
5                                            0.079078   
6                                            0.092044   
7                                            0.090382   
8                                            0.107120   
9                                                 NaN   
10                                           0.253547   
11                                           0.092187   
12                                          -0.046619   
13                                           0.100313   
14                                           0.102883   
15                                                NaN   
16                                           0.139269   
17                                          -0.046523   
18                                          -0.116332   
19                                           0.085315   
20                                           0.077589   
21                                                NaN   
22                                           0.151240   
23                                          -0.116294   

    spearman_abs_beta_change_vs_abs_zero_fraction_difference  \
0                                            0.329620          
1                                            0.034036          
2                                            0.029353          
3                                                 NaN          
4                                            0.012053          
5                                            0.329632          
6                                            0.357314          
7                                            0.014885          
8                                            0.012804          
9                                                 NaN          
10                                          -0.032230          
11                                           0.357332          
12                                           0.260747          
13                                          -0.018931          
14                                          -0.014740          
15                                                NaN          
16                                          -0.029062          
17                                           0.260703          
18                                           0.355522          
19                                          -0.056668          
20                                          -0.047734          
21                                                NaN          
22                                          -0.095808          
23                                           0.355524          

    spearman_signed_beta_change_vs_zero_fraction_difference  \
0                                           -0.729867         
1                                           -0.008715         
2                                           -0.003018         
3                                                 NaN         
4                                           -0.028787         
5                                            0.729922         
6                                           -0.726739         
7                                           -0.016057         
8                                           -0.013972         
9                                                 NaN         
10                                          -0.039974         
11                                           0.726740         
12                                          -0.639520         
13                                          -0.009284         
14                                          -0.012423         
15                                                NaN         
16                                           0.009212         
17                                           0.639579         
18                                          -0.669016         
19                                          -0.039407         
20                                          -0.043118         
21                                                NaN         
22                                           0.022998         
23                                           0.669013         

    spearman_program_mean_abs_beta_change_vs_stage2_positive_fraction_modified  \
0                                           -0.095062                            
1                                                 NaN                            
2                                            0.605692                            
3                                                 NaN                            
4                                            0.903088                            
5                                           -0.001064                            
6                                           -0.018672                            
7                                                 NaN                            
8                                            0.435606                            
9                                                 NaN                            
10                                           0.933781                            
11                                           0.129326                            
12                                           0.042506                            
13                                                NaN                            
14                                           0.584729                            
15                                                NaN                            
16                                           0.892900                            
17                                          -0.048062                            
18                                          -0.213925                            
19                                                NaN                            
20                                           0.674229                            
21                                                NaN                            
22                                           0.891066                            
23                                          -0.245482                            

    spearman_program_mean_abs_beta_change_vs_program_response_rmse  
0                                            0.352413               
1                                            0.934506               
2                                            0.938375               
3                                                 NaN               
4                                            0.971765               
5                                            0.352413               
6                                           -0.116062               
7                                            0.905386               
8                                            0.903365               
9                                                 NaN               
10                                           0.951784               
11                                          -0.115102               
12                                           0.453061               
13                                           0.932369               
14                                           0.930231               
15                                                NaN               
16                                           0.934246               
17                                           0.455366               
18                                           0.653782               
19                                           0.960911               
20                                           0.961205               
21                                                NaN               
22                                           0.919808               
23                                           0.652245               
```

### Text output (stream)

```text

Saved:
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_floor_inference_summary.csv
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_floor_effect_level_comparison.csv
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_floor_zero_and_program_scaling_summary.csv

```

## Markdown cell 21

## Plot inferential stability across both flooring dimensions

## Code cell 22

```python
from pathlib import Path
import numpy as np
import pandas as pd
from IPython.display import display
from plotnine import (
    ggplot, aes,
    geom_tile, geom_text, geom_line, geom_point, geom_hline,
    facet_wrap,
    scale_fill_gradient, scale_color_manual,
    labs, theme_bw, theme,
    element_text, element_line, element_rect, element_blank,
)

if "inference_summary" not in globals():
    inference_summary = pd.read_csv(MODEL_ROOT / "two_stage_floor_inference_summary.csv")

infer = inference_summary.copy()
infer["initial_display"] = infer["initial_floor_factor"].map(initial_display)
infer["q_display"] = infer["internal_eps_quantile"].map(q_display)
infer["initial_display"] = pd.Categorical(
    infer["initial_display"], categories=initial_levels, ordered=True
)
infer["q_display"] = pd.Categorical(
    infer["q_display"], categories=q_levels, ordered=True
)


def make_inference_heatmap(value_column, title, subtitle, fill_label, label_formatter, stem_name):
    current = infer.copy()
    current["tile_label"] = current[value_column].map(label_formatter)
    plot = (
        ggplot(current, aes("initial_display", "q_display", fill=value_column))
        + geom_tile(color="white", size=0.6)
        + geom_text(aes(label="tile_label"), size=7.2)
        + facet_wrap("~day", ncol=2)
        + scale_fill_gradient(low="#F7FBFF", high="#17365D")
        + labs(
            title=title,
            subtitle=subtitle,
            x="Initial exact-zero floor factor",
            y="Internal programDE eps quantile",
            fill=fill_label,
        )
        + theme_bw(base_size=10)
        + theme(
            figure_size=(8.2, 6.6),
            plot_title=element_text(size=13, weight="bold", ha="center"),
            plot_subtitle=element_text(size=8.5, ha="center"),
            axis_text_x=element_text(size=9),
            axis_text_y=element_text(size=8),
            strip_text=element_text(size=9),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(color=BORDER_COLOR, size=0.65),
            legend_position="right",
        )
    )
    display(plot)
    return save_vector_plot(plot, PLOT_DIR / stem_name, width=8.2, height=6.6)


beta_heatmap_saved = make_inference_heatmap(
    "beta_spearman_vs_reference",
    "DE program effect stability across both flooring choices",
    f"Empirical Bernoulli CRT, B={B:,}; reference is 1× and q={REFERENCE_Q:.0e}",
    "Beta rank correlation",
    lambda value: f"{value:.3f}",
    "beta_spearman_two_stage_heatmap",
)

jaccard_heatmap_saved = make_inference_heatmap(
    "discovery_jaccard_q_le_0_1",
    "Stability of q≤0.10 discoveries across flooring choices",
    "Jaccard overlap with the current two-stage flooring implementation",
    "Discovery Jaccard",
    lambda value: f"{value:.3f}",
    "discovery_jaccard_two_stage_heatmap",
)

ntc_rate_heatmap_saved = make_inference_heatmap(
    "ntc_fraction_p_lt_0_05",
    "NTC pseudo-target type-I error across flooring choices",
    "Ideal calibration is 0.05",
    "NTC fraction p<0.05",
    lambda value: f"{value:.3f}",
    "NTC_fraction_p_lt_0_05_two_stage_heatmap",
)

ntc_qq_heatmap_saved = make_inference_heatmap(
    "ntc_qq_rmse",
    "NTC pseudo-target QQ deviation across flooring choices",
    "Lower QQ RMSE indicates closer agreement with the uniform null",
    "NTC QQ RMSE",
    lambda value: f"{value:.3f}",
    "NTC_QQ_RMSE_two_stage_heatmap",
)

# Axis-isolating inferential curves.
axis_rows = []
for _, row in inference_summary.iterrows():
    if np.isclose(row["initial_floor_factor"], REFERENCE_INITIAL):
        current = row.to_dict()
        current["axis"] = "Vary internal floor; initial floor fixed at 1×"
        current["setting"] = q_display(row["internal_eps_quantile"])
        axis_rows.append(current)
    if np.isclose(row["internal_eps_quantile"], REFERENCE_Q):
        current = row.to_dict()
        current["axis"] = f"Vary initial floor; internal q fixed at {REFERENCE_Q:.0e}"
        current["setting"] = initial_display(row["initial_floor_factor"])
        axis_rows.append(current)

axis_infer = pd.DataFrame(axis_rows)
infer_metrics = {
    "beta_spearman_vs_reference": "Beta rank correlation vs current",
    "top_reference_effect_sign_concordance": "Top-500 effect sign concordance",
    "pvalue_neglog10_spearman_vs_reference": "-log10(p) rank correlation",
    "discovery_jaccard_q_le_0_1": "q≤0.10 discovery Jaccard",
    "ntc_fraction_p_lt_0_05": "NTC fraction p<0.05",
    "ntc_qq_rmse": "NTC QQ RMSE",
}
axis_infer_long = axis_infer.melt(
    id_vars=["day", "axis", "setting"],
    value_vars=list(infer_metrics),
    var_name="metric",
    value_name="value",
)
axis_infer_long["metric_label"] = axis_infer_long["metric"].map(infer_metrics)
axis_infer_long["panel"] = axis_infer_long["axis"] + "\n" + axis_infer_long["metric_label"]

p_axis_infer = (
    ggplot(axis_infer_long, aes("setting", "value", color="day", group="day"))
    + geom_hline(yintercept=0, color="#D5DCE3", size=0.3)
    + geom_line(size=0.75)
    + geom_point(size=2.0)
    + facet_wrap("~panel", scales="free", ncol=2)
    + scale_color_manual(values=DAY_COLORS)
    + labs(
        title="Axis-isolating inferential sensitivity",
        subtitle="Only empirical Bernoulli CRT p-values are compared",
        x="Flooring setting",
        y="Metric value",
        color="Day",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(9.2, 13.0),
        plot_title=element_text(size=13, weight="bold", ha="center"),
        plot_subtitle=element_text(size=8.5, ha="center"),
        axis_text_x=element_text(size=8, rotation=20, ha="right"),
        strip_text=element_text(size=8),
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
        legend_position="top",
    )
)
display(p_axis_infer)
axis_infer_saved = save_vector_plot(
    p_axis_infer,
    PLOT_DIR / "axis_isolating_inferential_sensitivity",
    width=9.2,
    height=13.0,
)

print("Saved inferential summary plots under:", PLOT_DIR)
for name, paths in {
    "beta rank correlation": beta_heatmap_saved,
    "discovery Jaccard": jaccard_heatmap_saved,
    "NTC p<0.05": ntc_rate_heatmap_saved,
    "NTC QQ RMSE": ntc_qq_heatmap_saved,
    "axis-isolating inference": axis_infer_saved,
}.items():
    print(f"\n{name}")
    for fmt, path in paths.items():
        print(f"  {fmt.upper()}: {path}")

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5d1032c090>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c3c4bce50>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c3c3628d0>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c3a870090>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c3c2bd290>
```

### Text output (stream)

```text
Saved inferential summary plots under: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity

beta rank correlation
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/beta_spearman_two_stage_heatmap.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/beta_spearman_two_stage_heatmap.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/beta_spearman_two_stage_heatmap.png

discovery Jaccard
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/discovery_jaccard_two_stage_heatmap.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/discovery_jaccard_two_stage_heatmap.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/discovery_jaccard_two_stage_heatmap.png

NTC p<0.05
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/NTC_fraction_p_lt_0_05_two_stage_heatmap.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/NTC_fraction_p_lt_0_05_two_stage_heatmap.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/NTC_fraction_p_lt_0_05_two_stage_heatmap.png

NTC QQ RMSE
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/NTC_QQ_RMSE_two_stage_heatmap.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/NTC_QQ_RMSE_two_stage_heatmap.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/NTC_QQ_RMSE_two_stage_heatmap.png

axis-isolating inference
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/axis_isolating_inferential_sensitivity.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/axis_isolating_inferential_sensitivity.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/axis_isolating_inferential_sensitivity.png

```

## Markdown cell 23

## NTC QQ curves and per-day beta scatterplots

## Code cell 24

```python
from pathlib import Path
import numpy as np
import pandas as pd
from IPython.display import display
from plotnine import (
    ggplot, aes,
    geom_abline, geom_line, geom_point, geom_hline, geom_vline,
    facet_wrap,
    scale_color_manual, scale_x_continuous, scale_y_continuous,
    coord_equal, labs, guides, guide_legend,
    theme_bw, theme,
    element_text, element_line, element_rect, element_blank,
)


def combination_display(initial_factor, internal_q, is_reference=False):
    initial = f"{float(initial_factor):g}×"
    if np.isclose(float(internal_q), 0.0):
        q_text = "q=0 (stage-1 only)"
    else:
        q_text = f"q={float(internal_q):.0e}"
    label = f"{initial}, {q_text}"
    if bool(is_reference):
        label += " (current)"
    return label


manifest_plot = manifest.copy()
manifest_plot["is_reference"] = as_bool(manifest_plot["is_reference"])
manifest_plot = manifest_plot[manifest_plot["status"].isin(["written", "reused_matching"])].copy()
manifest_plot["combination_display"] = manifest_plot.apply(
    lambda row: combination_display(
        row["initial_floor_factor"],
        row["internal_eps_quantile"],
        row["is_reference"],
    ),
    axis=1,
)

combination_order = (
    manifest_plot[[
        "initial_floor_factor",
        "internal_eps_quantile",
        "is_reference",
        "combination_display",
    ]]
    .drop_duplicates()
    .sort_values(
        ["is_reference", "initial_floor_factor", "internal_eps_quantile"],
        ascending=[False, True, True],
    )["combination_display"]
    .tolist()
)

palette_values = [
    "#1B3A57",
    "#3F88C5",
    "#74B3CE",
    "#8C6BB1",
    "#D17A22",
    "#B2182B",
    "#2A7F62",
    "#6B6B6B",
    "#C05A9D",
    "#8A9A5B",
    "#E3A018",
    "#4D4D4D",
]
COMBINATION_COLORS = {
    label: palette_values[index % len(palette_values)]
    for index, label in enumerate(combination_order)
}

# NTC QQ curves.
ntc_frames = []
for _, row in manifest_plot.iterrows():
    run_dir = Path(row["output_dir"])
    day = str(row["day"])
    pvalues = np.asarray(
        np.load(run_dir / f"{day}_ntc_pseudotarget_pvalues_empirical.npy"),
        dtype=float,
    ).ravel()
    pvalues = pvalues[np.isfinite(pvalues)]
    pvalues = np.clip(pvalues, 1e-300, 1.0)
    pvalues.sort()
    n = len(pvalues)
    expected_p = (np.arange(1, n + 1) - 0.5) / n
    ntc_frames.append(pd.DataFrame({
        "day": day,
        "combination_display": row["combination_display"],
        "expected_neglog10_p": -np.log10(expected_p),
        "observed_neglog10_p": -np.log10(pvalues),
    }))

ntc_qq = pd.concat(ntc_frames, ignore_index=True)
ntc_qq["combination_display"] = pd.Categorical(
    ntc_qq["combination_display"],
    categories=combination_order,
    ordered=True,
)
qq_limit = float(np.ceil(max(
    ntc_qq["expected_neglog10_p"].max(),
    ntc_qq["observed_neglog10_p"].max(),
) * 10) / 10)

p_ntc = (
    ggplot(
        ntc_qq,
        aes(
            "expected_neglog10_p",
            "observed_neglog10_p",
            color="combination_display",
            group="combination_display",
        ),
    )
    + geom_abline(
        intercept=0,
        slope=1,
        linetype="dashed",
        color="#AEB7C2",
        size=0.8,
    )
    + geom_line(size=0.75, alpha=0.88)
    + facet_wrap("~day", ncol=2)
    + scale_color_manual(values=COMBINATION_COLORS, drop=False)
    + scale_x_continuous(limits=(0, qq_limit))
    + scale_y_continuous(limits=(0, qq_limit))
    + coord_equal()
    + labs(
        title="NTC pseudo-target calibration across two-stage flooring choices",
        subtitle=f"Empirical Bernoulli CRT; B={B:,}; same cells, covariates, NTC groups, and seeds",
        x="Expected -log10(p)",
        y="Observed -log10(p)",
        color="Flooring combination",
    )
    + guides(color=guide_legend(ncol=2, override_aes={"size": 1.6, "alpha": 1}))
    + theme_bw(base_size=10)
    + theme(
        figure_size=(9.0, 8.0),
        plot_title=element_text(size=13, weight="bold", ha="center"),
        plot_subtitle=element_text(size=8.5, ha="center"),
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
        legend_position="top",
        legend_background=element_blank(),
        legend_key=element_blank(),
        legend_text=element_text(size=7.8),
    )
)
display(p_ntc)
ntc_saved = save_vector_plot(
    p_ntc,
    PLOT_DIR / "NTC_QQ_two_stage_floor_sensitivity",
    width=9.0,
    height=8.0,
)

# Per-day beta scatterplots against the current implementation.
effect_plot = effect_comparisons.copy()
effect_plot = effect_plot[~effect_plot["is_reference"].astype(bool)].copy()
effect_plot = effect_plot.merge(
    manifest_plot[["day", "combination_label", "combination_display"]].drop_duplicates(),
    on=["day", "combination_label"],
    how="left",
    validate="many_to_one",
)

summary_for_panels = inference_summary.copy()
summary_for_panels["is_reference"] = as_bool(summary_for_panels["is_reference"])
summary_for_panels = summary_for_panels[~summary_for_panels["is_reference"]].copy()
summary_for_panels = summary_for_panels.merge(
    manifest_plot[["day", "combination_label", "combination_display"]].drop_duplicates(),
    on=["day", "combination_label"],
    how="left",
    validate="one_to_one",
)
summary_for_panels["panel_label"] = summary_for_panels.apply(
    lambda row: (
        f"{row['combination_display']}\n"
        f"Spearman={row['beta_spearman_vs_reference']:.3f}; "
        f"CCC={row['beta_concordance_correlation_vs_reference']:.3f}"
    ),
    axis=1,
)
effect_plot = effect_plot.merge(
    summary_for_panels[["day", "combination_label", "panel_label"]],
    on=["day", "combination_label"],
    how="left",
    validate="many_to_one",
)

beta_saved_by_day = {}
for day in DAYS:
    day_data = effect_plot[effect_plot["day"].astype(str).eq(day)].copy()
    if day_data.empty:
        continue
    axis_limit = float(np.nanmax(np.abs(
        day_data[["beta_reference", "beta_current"]].to_numpy(dtype=float)
    )))
    axis_limit = axis_limit * 1.05 if np.isfinite(axis_limit) and axis_limit > 0 else 1.0
    n_panels = int(day_data["panel_label"].nunique())
    ncol = 2 if n_panels > 1 else 1
    height = max(4.5, 3.5 * np.ceil(n_panels / ncol))

    p_beta = (
        ggplot(day_data, aes("beta_reference", "beta_current"))
        + geom_hline(yintercept=0, color="#D5DCE3", size=0.3)
        + geom_vline(xintercept=0, color="#D5DCE3", size=0.3)
        + geom_abline(
            intercept=0,
            slope=1,
            linetype="dashed",
            color="#AEB7C2",
            size=0.75,
        )
        + geom_point(size=0.65, alpha=0.35, color="#17365D", stroke=0)
        + facet_wrap("~panel_label", ncol=ncol)
        + scale_x_continuous(limits=(-axis_limit, axis_limit))
        + scale_y_continuous(limits=(-axis_limit, axis_limit))
        + coord_equal()
        + labs(
            title=f"{day}: DE program beta sensitivity to two-stage flooring",
            subtitle="Each point is one real target × cNMF program effect",
            x="Beta under current two-stage flooring",
            y="Beta under alternative flooring",
        )
        + theme_bw(base_size=10)
        + theme(
            figure_size=(8.5, float(height)),
            plot_title=element_text(size=13, weight="bold", ha="center"),
            plot_subtitle=element_text(size=8.5, ha="center"),
            strip_text=element_text(size=7.8),
            panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(color=BORDER_COLOR, size=0.65),
        )
    )
    display(p_beta)
    beta_saved_by_day[day] = save_vector_plot(
        p_beta,
        PLOT_DIR / f"{day}_beta_scatter_two_stage_floor_sensitivity",
        width=8.5,
        height=float(height),
    )

print("Saved NTC QQ files:")
for fmt, path in ntc_saved.items():
    print(f"  {fmt.upper()}: {path}")

print("\nSaved beta scatter files:")
for day, paths in beta_saved_by_day.items():
    print(day)
    for fmt, path in paths.items():
        print(f"  {fmt.upper()}: {path}")

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c3e05a6d0>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c3ac54910>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c3ad34310>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c3c21ebd0>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c3c2a29d0>
```

### Text output (stream)

```text
Saved NTC QQ files:
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/NTC_QQ_two_stage_floor_sensitivity.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/NTC_QQ_two_stage_floor_sensitivity.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/NTC_QQ_two_stage_floor_sensitivity.png

Saved beta scatter files:
d0
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d0_beta_scatter_two_stage_floor_sensitivity.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d0_beta_scatter_two_stage_floor_sensitivity.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d0_beta_scatter_two_stage_floor_sensitivity.png
d1
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d1_beta_scatter_two_stage_floor_sensitivity.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d1_beta_scatter_two_stage_floor_sensitivity.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d1_beta_scatter_two_stage_floor_sensitivity.png
d2
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d2_beta_scatter_two_stage_floor_sensitivity.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d2_beta_scatter_two_stage_floor_sensitivity.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d2_beta_scatter_two_stage_floor_sensitivity.png
d3
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d3_beta_scatter_two_stage_floor_sensitivity.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d3_beta_scatter_two_stage_floor_sensitivity.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/d3_beta_scatter_two_stage_floor_sensitivity.png

```

## Markdown cell 25

## Test whether sensitivity scales with the amount of flooring

This diagnostic asks whether programs whose genuinely positive usages are more frequently truncated by Stage 2 also show larger changes in estimated target-program effects. The accompanying table also tests whether beta changes scale with raw-zero burden in treated versus control cells.

## Code cell 26

```python
from pathlib import Path
import numpy as np
import pandas as pd
from IPython.display import display
from plotnine import (
    ggplot, aes,
    geom_point, geom_smooth, geom_hline,
    facet_wrap,
    scale_color_manual,
    labs, theme_bw, theme,
    element_text, element_line, element_rect, element_blank,
)

# Aggregate effect sensitivity by program and relate it to how often stage 2
# truncates genuinely positive usages in that program.
program_beta_change = (
    effect_comparisons
    .groupby(["day", "combination_label", "program_name"], as_index=False)
    .agg(
        mean_abs_beta_change=("abs_beta_delta", "mean"),
        median_abs_beta_change=("abs_beta_delta", "median"),
        p95_abs_beta_change=("abs_beta_delta", lambda x: x.quantile(0.95)),
    )
)
program_sensitivity = program_beta_change.merge(
    program_audit,
    on=["day", "combination_label", "program_name"],
    how="left",
    validate="one_to_one",
)
program_sensitivity = program_sensitivity.merge(
    manifest_plot[["day", "combination_label", "combination_display"]].drop_duplicates(),
    on=["day", "combination_label"],
    how="left",
    validate="many_to_one",
)
program_sensitivity["is_reference"] = as_bool(program_sensitivity["is_reference"])
program_sensitivity = program_sensitivity[~program_sensitivity["is_reference"]].copy()
program_sensitivity["panel"] = program_sensitivity["day"].astype(str) + ": " + program_sensitivity["combination_display"].astype(str)

p_program = (
    ggplot(
        program_sensitivity,
        aes(
            "stage2_fraction_originally_positive_entries_modified",
            "mean_abs_beta_change",
        ),
    )
    + geom_point(size=1.35, alpha=0.70, color="#17365D", stroke=0)
    + geom_smooth(method="lm", se=False, linetype="dashed", size=0.65, color="#B2182B")
    + facet_wrap("~panel", scales="free", ncol=2)
    + labs(
        title="Does second-floor truncation predict beta sensitivity?",
        subtitle="Each point is one cNMF program; y is averaged across real targets",
        x="Fraction of originally positive usages modified by stage 2",
        y="Mean absolute beta change versus current flooring",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(9.0, 15.0),
        plot_title=element_text(size=13, weight="bold", ha="center"),
        plot_subtitle=element_text(size=8.5, ha="center"),
        strip_text=element_text(size=7.4),
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
    )
)
display(p_program)
program_saved = save_vector_plot(
    p_program,
    PLOT_DIR / "program_beta_change_vs_stage2_positive_truncation",
    width=9.0,
    height=15.0,
)

PROGRAM_SENSITIVITY_PATH = MODEL_ROOT / "two_stage_program_level_beta_sensitivity.csv"
program_sensitivity.to_csv(PROGRAM_SENSITIVITY_PATH, index=False)

print("Zero/stage-2 scaling summary:")
display(scaling_summary)

print("\nSaved:")
print(" ", PROGRAM_SENSITIVITY_PATH)
for fmt, path in program_saved.items():
    print(f"  {fmt.upper()}: {path}")

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f5c3c354d50>
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/stats/stat_smooth.py:213: PlotnineWarning: Smoothing requires 2 or more points. Got 1. Not enough points for smoothing. If this message a surprise, make sure the column mapped to the x aesthetic has the right dtype.

```

### Text output (stream)

```text
Zero/stage-2 scaling summary:

```

### Text output (text/plain)

```text
   day  initial_floor_factor  internal_eps_quantile  \
0   d0                   0.1                0.00010   
1   d0                   1.0                0.00000   
2   d0                   1.0                0.00001   
3   d0                   1.0                0.00010   
4   d0                   1.0                0.00100   
5   d0                  10.0                0.00010   
6   d1                   0.1                0.00010   
7   d1                   1.0                0.00000   
8   d1                   1.0                0.00001   
9   d1                   1.0                0.00010   
10  d1                   1.0                0.00100   
11  d1                  10.0                0.00010   
12  d2                   0.1                0.00010   
13  d2                   1.0                0.00000   
14  d2                   1.0                0.00001   
15  d2                   1.0                0.00010   
16  d2                   1.0                0.00100   
17  d2                  10.0                0.00010   
18  d3                   0.1                0.00010   
19  d3                   1.0                0.00000   
20  d3                   1.0                0.00001   
21  d3                   1.0                0.00010   
22  d3                   1.0                0.00100   
23  d3                  10.0                0.00010   

                combination_label  is_reference  n_pairs_with_zero_burden  \
0   initial_0_1x__internal_q1em04         False                     14900   
1         initial_1x__internal_q0         False                     14900   
2     initial_1x__internal_q1em05         False                     14900   
3     initial_1x__internal_q1em04          True                     14900   
4     initial_1x__internal_q1em03         False                     14900   
5    initial_10x__internal_q1em04         False                     14900   
6   initial_0_1x__internal_q1em04         False                     14900   
7         initial_1x__internal_q0         False                     14900   
8     initial_1x__internal_q1em05         False                     14900   
9     initial_1x__internal_q1em04          True                     14900   
10    initial_1x__internal_q1em03         False                     14900   
11   initial_10x__internal_q1em04         False                     14900   
12  initial_0_1x__internal_q1em04         False                     14900   
13        initial_1x__internal_q0         False                     14900   
14    initial_1x__internal_q1em05         False                     14900   
15    initial_1x__internal_q1em04          True                     14900   
16    initial_1x__internal_q1em03         False                     14900   
17   initial_10x__internal_q1em04         False                     14900   
18  initial_0_1x__internal_q1em04         False                     14900   
19        initial_1x__internal_q0         False                     14900   
20    initial_1x__internal_q1em05         False                     14900   
21    initial_1x__internal_q1em04          True                     14900   
22    initial_1x__internal_q1em03         False                     14900   
23   initial_10x__internal_q1em04         False                     14900   

    spearman_abs_beta_change_vs_zero_fraction_treated  \
0                                            0.079068   
1                                            0.164578   
2                                            0.173769   
3                                                 NaN   
4                                            0.373217   
5                                            0.079078   
6                                            0.092044   
7                                            0.090382   
8                                            0.107120   
9                                                 NaN   
10                                           0.253547   
11                                           0.092187   
12                                          -0.046619   
13                                           0.100313   
14                                           0.102883   
15                                                NaN   
16                                           0.139269   
17                                          -0.046523   
18                                          -0.116332   
19                                           0.085315   
20                                           0.077589   
21                                                NaN   
22                                           0.151240   
23                                          -0.116294   

    spearman_abs_beta_change_vs_abs_zero_fraction_difference  \
0                                            0.329620          
1                                            0.034036          
2                                            0.029353          
3                                                 NaN          
4                                            0.012053          
5                                            0.329632          
6                                            0.357314          
7                                            0.014885          
8                                            0.012804          
9                                                 NaN          
10                                          -0.032230          
11                                           0.357332          
12                                           0.260747          
13                                          -0.018931          
14                                          -0.014740          
15                                                NaN          
16                                          -0.029062          
17                                           0.260703          
18                                           0.355522          
19                                          -0.056668          
20                                          -0.047734          
21                                                NaN          
22                                          -0.095808          
23                                           0.355524          

    spearman_signed_beta_change_vs_zero_fraction_difference  \
0                                           -0.729867         
1                                           -0.008715         
2                                           -0.003018         
3                                                 NaN         
4                                           -0.028787         
5                                            0.729922         
6                                           -0.726739         
7                                           -0.016057         
8                                           -0.013972         
9                                                 NaN         
10                                          -0.039974         
11                                           0.726740         
12                                          -0.639520         
13                                          -0.009284         
14                                          -0.012423         
15                                                NaN         
16                                           0.009212         
17                                           0.639579         
18                                          -0.669016         
19                                          -0.039407         
20                                          -0.043118         
21                                                NaN         
22                                           0.022998         
23                                           0.669013         

    spearman_program_mean_abs_beta_change_vs_stage2_positive_fraction_modified  \
0                                           -0.095062                            
1                                                 NaN                            
2                                            0.605692                            
3                                                 NaN                            
4                                            0.903088                            
5                                           -0.001064                            
6                                           -0.018672                            
7                                                 NaN                            
8                                            0.435606                            
9                                                 NaN                            
10                                           0.933781                            
11                                           0.129326                            
12                                           0.042506                            
13                                                NaN                            
14                                           0.584729                            
15                                                NaN                            
16                                           0.892900                            
17                                          -0.048062                            
18                                          -0.213925                            
19                                                NaN                            
20                                           0.674229                            
21                                                NaN                            
22                                           0.891066                            
23                                          -0.245482                            

    spearman_program_mean_abs_beta_change_vs_program_response_rmse  
0                                            0.352413               
1                                            0.934506               
2                                            0.938375               
3                                                 NaN               
4                                            0.971765               
5                                            0.352413               
6                                           -0.116062               
7                                            0.905386               
8                                            0.903365               
9                                                 NaN               
10                                           0.951784               
11                                          -0.115102               
12                                           0.453061               
13                                           0.932369               
14                                           0.930231               
15                                                NaN               
16                                           0.934246               
17                                           0.455366               
18                                           0.653782               
19                                           0.960911               
20                                           0.961205               
21                                                NaN               
22                                           0.919808               
23                                           0.652245               
```

### Text output (stream)

```text

Saved:
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_program_level_beta_sensitivity.csv
  PDF: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/program_beta_change_vs_stage2_positive_truncation.pdf
  SVG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/program_beta_change_vs_stage2_positive_truncation.svg
  PNG: results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/illustrator_editable_two_stage_floor_sensitivity/program_beta_change_vs_stage2_positive_truncation.png

```

## Markdown cell 27

## Manuscript-facing summary table and reporting text

## Code cell 28

```python
# Compact manuscript-facing summary and interpretation flags.
summary_for_report = inference_summary.merge(
    audit[[
        "day",
        "combination_label",
        "applied_raw_zero_floor",
        "effective_internal_floor",
        "fraction_exact_zero_raw_entries",
        "stage2_fraction_entries_modified",
        "stage2_fraction_originally_positive_entries_modified",
        "stage2_fraction_zero_derived_entries_modified",
        "median_initial_pseudomass_fraction",
        "median_stage2_pseudomass_fraction",
        "response_entrywise_spearman_vs_reference",
        "response_p95_abs_change_vs_reference",
    ]],
    on=["day", "combination_label"],
    how="left",
    validate="one_to_one",
)
summary_for_report = summary_for_report.merge(
    manifest_plot[["day", "combination_label", "combination_display"]].drop_duplicates(),
    on=["day", "combination_label"],
    how="left",
    validate="one_to_one",
)

summary_for_report["axis_tested"] = np.select(
    [
        np.isclose(summary_for_report["initial_floor_factor"], REFERENCE_INITIAL)
        & ~np.isclose(summary_for_report["internal_eps_quantile"], REFERENCE_Q),
        ~np.isclose(summary_for_report["initial_floor_factor"], REFERENCE_INITIAL)
        & np.isclose(summary_for_report["internal_eps_quantile"], REFERENCE_Q),
        np.isclose(summary_for_report["initial_floor_factor"], REFERENCE_INITIAL)
        & np.isclose(summary_for_report["internal_eps_quantile"], REFERENCE_Q),
    ],
    [
        "internal floor varied",
        "initial floor varied",
        "current reference",
    ],
    default="both floors varied",
)

summary_for_report["response_rank_stable_ge_0_99"] = (
    summary_for_report["response_entrywise_spearman_vs_reference"] >= 0.99
)
summary_for_report["beta_rank_stable_ge_0_95"] = (
    summary_for_report["beta_spearman_vs_reference"] >= 0.95
)
summary_for_report["discovery_jaccard_ge_0_80"] = (
    summary_for_report["discovery_jaccard_q_le_0_1"] >= 0.80
)
summary_for_report["ntc_p05_within_0_02_of_nominal"] = (
    (summary_for_report["ntc_fraction_p_lt_0_05"] - 0.05).abs() <= 0.02
)

report_columns = [
    "day",
    "combination_display",
    "axis_tested",
    "is_reference",
    "applied_raw_zero_floor",
    "effective_internal_floor",
    "fraction_exact_zero_raw_entries",
    "stage2_fraction_originally_positive_entries_modified",
    "median_initial_pseudomass_fraction",
    "median_stage2_pseudomass_fraction",
    "response_entrywise_spearman_vs_reference",
    "response_p95_abs_change_vs_reference",
    "beta_spearman_vs_reference",
    "beta_concordance_correlation_vs_reference",
    "top_reference_effect_sign_concordance",
    "pvalue_neglog10_spearman_vs_reference",
    "discovery_jaccard_q_le_0_05",
    "discovery_jaccard_q_le_0_1",
    "ntc_fraction_p_lt_0_05",
    "ntc_qq_rmse",
    "ntc_qq_tail_rmse",
    "ntc_ks_distance_uniform",
    "response_rank_stable_ge_0_99",
    "beta_rank_stable_ge_0_95",
    "discovery_jaccard_ge_0_80",
    "ntc_p05_within_0_02_of_nominal",
]

manuscript_summary = summary_for_report[report_columns].sort_values(
    ["day", "axis_tested", "combination_display"]
)
MANUSCRIPT_SUMMARY_PATH = MODEL_ROOT / "two_stage_floor_manuscript_summary.csv"
manuscript_summary.to_csv(MANUSCRIPT_SUMMARY_PATH, index=False)

display(manuscript_summary)

methods_text = f"""
Exact-zero cNMF usages were first replaced using a day-specific floor,
ε0,d = max[1e-8, 0.1 × Q0.01(Uraw,d | Uraw,d > 0)], and each cell's usage
vector was renormalized. The current analysis then applied the programDE
lower-tail floor at quantile q={REFERENCE_Q:.0e} before a second row
renormalization and CLR transformation. Sensitivity analyses independently
varied the initial floor over {INITIAL_FLOOR_FACTORS}× and the internal
quantile over {INTERNAL_EPS_QUANTILES}, including q=0 as a stage-1-only
condition. The full 3×4 grid was evaluated at the transformed-response level.
Empirical Bernoulli CRT inference used B={B:,} resamples for the axis-isolating
{CRT_GRID_MODE!r} grid, while cells, six-source covariates, guide assignments,
NTC pseudo-target groups, and random seeds were held fixed.
""".strip()

interpretation_text = """
Interpretation guide:
- Stage-2 positive-entry modification quantifies whether the internal floor
  truncates genuinely positive cNMF usages rather than only zero-derived values.
- CLR and beta correlations assess numerical and inferential stability.
- Discovery Jaccard assesses stability of q-value conclusions.
- NTC nominal rates, QQ RMSE, tail RMSE, and KS distance assess calibration.
- Correlation of beta changes with raw-zero burden or stage-2 truncation helps
  identify whether sensitivity scales with the amount of flooring applied.
""".strip()

REPORT_TEXT_PATH = MODEL_ROOT / "two_stage_floor_reporting_text.txt"
REPORT_TEXT_PATH.write_text(methods_text + "\n\n" + interpretation_text + "\n")

print("\nSuggested Methods text:\n")
print(methods_text)
print("\nSaved:")
print(" ", MANUSCRIPT_SUMMARY_PATH)
print(" ", REPORT_TEXT_PATH)

```

### Text output (text/plain)

```text
   day     combination_display            axis_tested  is_reference  \
3   d0   1×, q=1e-04 (current)      current reference          True   
0   d0           0.1×, q=1e-04   initial floor varied         False   
5   d0            10×, q=1e-04   initial floor varied         False   
1   d0  1×, q=0 (stage-1 only)  internal floor varied         False   
4   d0             1×, q=1e-03  internal floor varied         False   
2   d0             1×, q=1e-05  internal floor varied         False   
9   d1   1×, q=1e-04 (current)      current reference          True   
6   d1           0.1×, q=1e-04   initial floor varied         False   
11  d1            10×, q=1e-04   initial floor varied         False   
7   d1  1×, q=0 (stage-1 only)  internal floor varied         False   
10  d1             1×, q=1e-03  internal floor varied         False   
8   d1             1×, q=1e-05  internal floor varied         False   
15  d2   1×, q=1e-04 (current)      current reference          True   
12  d2           0.1×, q=1e-04   initial floor varied         False   
17  d2            10×, q=1e-04   initial floor varied         False   
13  d2  1×, q=0 (stage-1 only)  internal floor varied         False   
16  d2             1×, q=1e-03  internal floor varied         False   
14  d2             1×, q=1e-05  internal floor varied         False   
21  d3   1×, q=1e-04 (current)      current reference          True   
18  d3           0.1×, q=1e-04   initial floor varied         False   
23  d3            10×, q=1e-04   initial floor varied         False   
19  d3  1×, q=0 (stage-1 only)  internal floor varied         False   
22  d3             1×, q=1e-03  internal floor varied         False   
20  d3             1×, q=1e-05  internal floor varied         False   

    applied_raw_zero_floor  effective_internal_floor  \
3             8.001929e-06              1.563586e-06   
0             8.001929e-07              5.381165e-07   
5             8.001929e-05              1.561712e-06   
1             8.001929e-06              1.863838e-09   
4             8.001929e-06              6.001187e-06   
2             8.001929e-06              1.041671e-07   
9             8.325819e-06              1.566422e-06   
6             8.325819e-07              6.358292e-07   
11            8.325819e-05              1.563594e-06   
7             8.325819e-06              4.331922e-08   
10            8.325819e-06              6.545783e-06   
8             8.325819e-06              2.233028e-07   
15            6.922590e-06              1.689523e-06   
12            6.922590e-07              4.771652e-07   
17            6.922590e-05              1.686318e-06   
13            6.922590e-06              3.142843e-08   
16            6.922590e-06              4.963874e-06   
14            6.922590e-06              1.625242e-07   
21            6.970817e-06              1.399135e-06   
18            6.970817e-07              4.511892e-07   
23            6.970817e-05              1.397006e-06   
19            6.970817e-06              5.880648e-13   
22            6.970817e-06              4.981785e-06   
20            6.970817e-06              1.884589e-07   

    fraction_exact_zero_raw_entries  \
3                          0.446590   
0                          0.446590   
5                          0.446590   
1                          0.446590   
4                          0.446590   
2                          0.446590   
9                          0.507464   
6                          0.507464   
11                         0.507464   
7                          0.507464   
10                         0.507464   
8                          0.507464   
15                         0.565127   
12                         0.565127   
17                         0.565127   
13                         0.565127   
16                         0.565127   
14                         0.565127   
21                         0.548154   
18                         0.548154   
23                         0.548154   
19                         0.548154   
22                         0.548154   
20                         0.548154   

    stage2_fraction_originally_positive_entries_modified  \
3                                            0.000182      
0                                            0.000063      
5                                            0.000182      
1                                            0.000000      
4                                            0.000795      
2                                            0.000019      
9                                            0.000204      
6                                            0.000082      
11                                           0.000204      
7                                            0.000000      
10                                           0.000821      
8                                            0.000022      
15                                           0.000230      
12                                           0.000075      
17                                           0.000230      
13                                           0.000000      
16                                           0.000707      
14                                           0.000023      
21                                           0.000222      
18                                           0.000067      
23                                           0.000222      
19                                           0.000000      
22                                           0.000728      
20                                           0.000023      

    median_initial_pseudomass_fraction  median_stage2_pseudomass_fraction  \
3                             0.000176                                0.0   
0                             0.000018                                0.0   
5                             0.001756                                0.0   
1                             0.000176                                0.0   
4                             0.000176                                0.0   
2                             0.000176                                0.0   
9                             0.000209                                0.0   
6                             0.000021                                0.0   
11                            0.002091                                0.0   
7                             0.000209                                0.0   
10                            0.000209                                0.0   
8                             0.000209                                0.0   
15                            0.000194                                0.0   
12                            0.000019                                0.0   
17                            0.001939                                0.0   
13                            0.000194                                0.0   
16                            0.000194                                0.0   
14                            0.000194                                0.0   
21                            0.000190                                0.0   
18                            0.000019                                0.0   
23                            0.001899                                0.0   
19                            0.000190                                0.0   
22                            0.000190                                0.0   
20                            0.000190                                0.0   

    ...  discovery_jaccard_q_le_0_05  discovery_jaccard_q_le_0_1  \
3   ...                     1.000000                    1.000000   
0   ...                     0.895137                    0.910088   
5   ...                     0.867792                    0.873482   
1   ...                     0.933333                    0.945593   
4   ...                     0.937785                    0.940282   
2   ...                     0.929323                    0.940669   
9   ...                     1.000000                    1.000000   
6   ...                     0.771429                    0.894366   
11  ...                     0.797235                    0.828025   
7   ...                     0.811321                    0.910345   
10  ...                     0.835616                    0.933798   
8   ...                     0.840183                    0.933566   
15  ...                     1.000000                    1.000000   
12  ...                     0.873656                    0.867562   
17  ...                     0.846715                    0.850267   
13  ...                     0.925134                    0.927103   
16  ...                     0.917772                    0.915870   
14  ...                     0.907407                    0.915572   
21  ...                     1.000000                    1.000000   
18  ...                     0.919802                    0.928732   
23  ...                     0.886598                    0.891720   
19  ...                     0.952144                    0.967645   
22  ...                     0.951437                    0.962547   
20  ...                     0.954183                    0.960240   

    ntc_fraction_p_lt_0_05  ntc_qq_rmse  ntc_qq_tail_rmse  \
3                 0.060656     0.032617          0.073421   
0                 0.059836     0.031654          0.069202   
5                 0.060328     0.033657          0.089626   
1                 0.060820     0.032563          0.073898   
4                 0.060820     0.032654          0.073980   
2                 0.060820     0.032579          0.073986   
9                 0.057204     0.022570          0.086654   
6                 0.056774     0.021458          0.089591   
11                0.058710     0.025248          0.086614   
7                 0.056989     0.022459          0.085629   
10                0.057204     0.022574          0.087874   
8                 0.056989     0.022472          0.085629   
15                0.053543     0.008694          0.064003   
12                0.053071     0.010517          0.085927   
17                0.050709     0.010799          0.083190   
13                0.053386     0.008777          0.064375   
16                0.053071     0.008640          0.064023   
14                0.053386     0.008767          0.064437   
21                0.055742     0.032500          0.254527   
18                0.055613     0.028006          0.219553   
23                0.057032     0.037830          0.272255   
19                0.055871     0.032574          0.254151   
22                0.056000     0.032351          0.253952   
20                0.055613     0.032520          0.254011   

    ntc_ks_distance_uniform  response_rank_stable_ge_0_99  \
3                  0.019763                          True   
0                  0.018315                          True   
5                  0.023588                          True   
1                  0.019681                          True   
4                  0.020922                          True   
2                  0.019681                          True   
9                  0.016353                          True   
6                  0.015569                          True   
11                 0.020477                          True   
7                  0.015823                          True   
10                 0.016068                          True   
8                  0.015823                          True   
15                 0.006253                          True   
12                 0.006328                          True   
17                 0.011335                          True   
13                 0.006423                          True   
16                 0.006725                          True   
14                 0.006523                          True   
21                 0.016035                          True   
18                 0.012921                          True   
23                 0.017524                          True   
19                 0.015906                          True   
22                 0.016052                          True   
20                 0.015874                          True   

    beta_rank_stable_ge_0_95  discovery_jaccard_ge_0_80  \
3                       True                       True   
0                       True                       True   
5                       True                       True   
1                       True                       True   
4                       True                       True   
2                       True                       True   
9                       True                       True   
6                       True                       True   
11                      True                       True   
7                       True                       True   
10                      True                       True   
8                       True                       True   
15                      True                       True   
12                      True                       True   
17                      True                       True   
13                      True                       True   
16                      True                       True   
14                      True                       True   
21                      True                       True   
18                      True                       True   
23                      True                       True   
19                      True                       True   
22                      True                       True   
20                      True                       True   

    ntc_p05_within_0_02_of_nominal  
3                             True  
0                             True  
5                             True  
1                             True  
4                             True  
2                             True  
9                             True  
6                             True  
11                            True  
7                             True  
10                            True  
8                             True  
15                            True  
12                            True  
17                            True  
13                            True  
16                            True  
14                            True  
21                            True  
18                            True  
23                            True  
19                            True  
22                            True  
20                            True  

[24 rows x 26 columns]
```

### Text output (stream)

```text

Suggested Methods text:

Exact-zero cNMF usages were first replaced using a day-specific floor,
ε0,d = max[1e-8, 0.1 × Q0.01(Uraw,d | Uraw,d > 0)], and each cell's usage
vector was renormalized. The current analysis then applied the programDE
lower-tail floor at quantile q=1e-04 before a second row
renormalization and CLR transformation. Sensitivity analyses independently
varied the initial floor over [0.1, 1.0, 10.0]× and the internal
quantile over [0.0, 1e-05, 0.0001, 0.001], including q=0 as a stage-1-only
condition. The full 3×4 grid was evaluated at the transformed-response level.
Empirical Bernoulli CRT inference used B=10,000 resamples for the axis-isolating
'primary' grid, while cells, six-source covariates, guide assignments,
NTC pseudo-target groups, and random seeds were held fixed.

Saved:
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_floor_manuscript_summary.csv
  results/ipsc_ec_clr_two_stage_floor_sensitivity_6cov/k_50.dt_0_2/two_stage_floor_reporting_text.txt

```

## Markdown cell 29

## Interpreting the result

A strong robustness result would have all of the following:

- very small Stage-2 modification of originally positive usages;
- near-unity CLR and beta rank correlations across the tested settings;
- high sign concordance for the strongest effects;
- stable discovery sets at the chosen q-value threshold;
- similar NTC nominal rates and QQ deviations;
- little relationship between beta changes and raw-zero burden or Stage-2 truncation.

If `q=0` differs materially from the current condition, the second flooring step is consequential and should be reported explicitly. If changing the initial factor at fixed internal `q` has larger effects, sensitivity is mainly driven by exact-zero replacement. If both axes matter, report the two-stage procedure and the full sensitivity grid rather than describing a single pseudocount.

