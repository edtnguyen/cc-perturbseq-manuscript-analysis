# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/fig3_carrierseeded_neighborhood_plotnine_style_notebook_package/fig3_carrierseeded_neighborhood_panels_plotnine_style.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/fig3_carrierseeded_neighborhood_plotnine_style_notebook_package/fig3_carrierseeded_neighborhood_panels_plotnine_style.ipynb`
- SHA-256: `212a7910cc939459f91719d27b8425bd0408a4ca7842c4e61f0a2c03bb8576b6`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0


# Fig. 3 panel f/j carrier-seeded neighborhood maps

This notebook is the interactive version of the Sherlock script package. It matches the UMAP-neighborhood logic from `ipsc-ec.Hippo.figure.ipynb`:

1. Use **target-carrier cells as seed cells**.
2. Find nearest neighbors in **z-scored cNMF program-usage space**.
3. Restrict neighbors to the **same day** as each seed.
4. Compare local target-carrier counts to **same-day non-targeting-control carriers**, excluding cells carrying the tested target from the NT comparison.
5. Compute local enrichment and Fisher/BH FDR.
6. Plot UMAP maps for panels **f** and **j**, displaying the figures inline with the same plotnine style used in the previous Hippo UMAP-neighborhood notebook.

Run this notebook cell-by-cell on Sherlock so you can monitor QC tables and plots as they are generated.


**Plotting style update:** the panel f/j UMAPs now use the earlier plotnine simple-color style: gray all-cell background, colored carrier-seeded seed cells, diverging blue-white-red log2OR scale, shared all-day UMAP coordinates, and Illustrator-editable SVG/PDF text settings.

## Markdown cell 1


## 1. Editable configuration

The paths below were extracted from your uploaded `ipsc-ec.Hippo.figure.ipynb`. Update `GUIDE_METADATA_TSV` if the relative path is not valid from the notebook working directory.


## Code cell 2

```python

from pathlib import Path

class cfg:
    # Extracted from uploaded notebook
    ADATA_H5AD = "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad"
    CNMF_USAGE_TXT = "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50.usages.k_50.dt_0_2.consensus.txt"

    # The notebook used this as a relative path. Change to absolute path if needed.
    GUIDE_METADATA_TSV = "300genes_guide_metadata_v43.tsv"

    OUTDIR = "fig3_carrierseeded_neighborhood_outputs"

    # AnnData keys / columns
    DAY_COL = "day"
    CELLTYPE_COL = "cell_type"
    GUIDE_NAMES_KEY = "guide_names"
    GUIDE_ASSIGNMENT_KEY = "guide_assignment"
    UMAP_KEY = "X_umap"
    CNMF_KEY = "cnmf_usage"

    # Guide metadata columns
    GUIDE_ID_COL = "guide_id"
    TARGET_COL_CANDIDATES = ["putative_target_genes", "intended_target_name", "target", "target_name", "gene"]

    # Panel target lists
    PANEL_F_TARGETS = ["ACVR2A", "BMPR1A", "SMAD1", "SMAD5"]
    PANEL_J_TARGETS = ["ETV2", "HAND1", "PRKCD", "PLCG1"]

    # Controls
    CONTROL_TARGET = "non-targeting_control"
    CONTROL_PATTERNS = [
        r"non[-_ ]?target",
        r"non[-_ ]?targeting",
        r"\bntc\b",
        r"negative",
        r"neg[-_ ]?ctrl",
        r"scram",
        r"scrambled",
        r"safe",
        r"aavs1",
        r"control",
    ]

    # KNN neighborhood settings matching the uploaded notebook section
    PROGRAMS_FOR_KNN = None  # None = all cNMF programs
    CNMF_TRANSFORM = "zscore"
    K = 300
    PSEUDOCOUNT = 0.5
    RESTRICT_NEIGHBORS_TO_SAME_DAY = True
    SAME_DAY_BASELINE = True
    MODE = "all_carriers"
    DAYS_TO_USE = ["d0", "d1", "d2", "d3"]
    LOG2FC_CLIP = 3.0
    RANDOM_SEED = 1
    EXPORT_NEIGHBOR_MEMBERSHIP = False

    # Plot settings copied from previous Hippo UMAP-neighborhood notebook style
    MAX_BACKGROUND_CELLS = 60000
    BACKGROUND_POINT_SIZE = 0.13
    BACKGROUND_ALPHA = 0.12
    FOREGROUND_POINT_SIZE = 1.7
    CMAP_LOW = "#053061"
    CMAP_MID = "#F7F7F7"
    CMAP_HIGH = "#67001F"
    PLOT_CLIP = 2.0
    FACET_NCOL = 2
    FIGSIZE_GRID = (8.4, 7.2)
    RASTERIZE_BACKGROUND_FOR_SPEED = True  # set False for fully vector background points

OUTDIR = Path(cfg.OUTDIR)
OUTDIR.mkdir(parents=True, exist_ok=True)
print("Output directory:", OUTDIR.resolve())
print("adata:", cfg.ADATA_H5AD)
print("cNMF usage:", cfg.CNMF_USAGE_TXT)
print("guide metadata:", cfg.GUIDE_METADATA_TSV)

```

## Markdown cell 3

## 2. Imports

## Code cell 4

```python

import os
import re
import zipfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from matplotlib.colors import LinearSegmentedColormap
from IPython.display import display, Image, Markdown

warnings.filterwarnings("ignore")

mpl.rcParams.update({
    "svg.fonttype": "none",   # preserve text as editable text in SVG
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
})

CONTROL_TARGET = cfg.CONTROL_TARGET


# plotnine style matching previous Hippo UMAP-neighborhood notebook
from plotnine import (
    ggplot, aes, geom_point, facet_wrap,
    scale_color_gradient2, coord_equal, labs,
    theme_bw, theme,
    element_blank, element_rect, element_text
)
from matplotlib import font_manager

available_fonts = {f.name for f in font_manager.fontManager.ttflist}
FONT_CANDIDATES = ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]
AI_FONT = next((f for f in FONT_CANDIDATES if f in available_fonts), "DejaVu Sans")
print("Using plot font:", AI_FONT)

```

## Markdown cell 5

## 3. Helper functions: guide/target parsing and assignment

## Code cell 6

```python

def looks_control(x):
    s = str(x).lower()
    return any(re.search(p, s) for p in cfg.CONTROL_PATTERNS)


def clean_target_string(x):
    if pd.isna(x):
        return []
    s = str(x).strip()
    if s.lower() in ["", "nan", "none", "na"]:
        return []
    s = s.replace("[", "").replace("]", "").replace("'", "").replace('"', "")
    parts = re.split(r"[;,|/]+", s)
    return [p.strip() for p in parts if p.strip()]


def canonical_guide_id(x):
    return str(x).strip().replace(" ", "")


def clean_day(x):
    s = str(x).strip().lower()
    if s in ["0", "d0", "day0", "day 0"]:
        return "d0"
    if s in ["1", "d1", "day1", "day 1"]:
        return "d1"
    if s in ["2", "d2", "day2", "day 2"]:
        return "d2"
    if s in ["3", "d3", "day3", "day 3"]:
        return "d3"
    return s


def get_guide_assignment_matrix(adata):
    ga = adata.obsm[cfg.GUIDE_ASSIGNMENT_KEY]
    if sparse.issparse(ga):
        return ga.tocsr().astype(bool)
    return np.asarray(ga).astype(bool)


def any_guides_per_cell(guide_assignment, guide_indices):
    if len(guide_indices) == 0:
        return np.zeros(guide_assignment.shape[0], dtype=bool)
    if sparse.issparse(guide_assignment):
        return np.asarray(guide_assignment[:, guide_indices].sum(axis=1)).ravel() > 0
    return guide_assignment[:, guide_indices].any(axis=1)


def choose_target_col(meta):
    for col in cfg.TARGET_COL_CANDIDATES:
        if col in meta.columns:
            return col
    raise ValueError(f"None of TARGET_COL_CANDIDATES found in metadata. Columns={meta.columns.tolist()}")


def load_guide_target_from_metadata(adata):
    guide_names = list(map(str, adata.uns[cfg.GUIDE_NAMES_KEY]))
    guide_name_to_idx = {g: i for i, g in enumerate(guide_names)}

    meta = pd.read_csv(cfg.GUIDE_METADATA_TSV, sep="\t")
    if cfg.GUIDE_ID_COL not in meta.columns:
        raise ValueError(f"{cfg.GUIDE_ID_COL!r} missing from {cfg.GUIDE_METADATA_TSV}")
    target_col = choose_target_col(meta)
    print(f"Using metadata target column: {target_col}")

    meta = meta[[cfg.GUIDE_ID_COL, target_col]].copy()
    meta.columns = ["guide_id", "target_entry"]
    meta["guide_id"] = meta["guide_id"].astype(str)

    meta["guide_index"] = meta["guide_id"].map(guide_name_to_idx)
    meta["match_method"] = np.where(meta["guide_index"].notna(), "exact_guide_id_match", None)

    # fallback: match after removing spaces
    if meta["guide_index"].isna().any():
        canonical_to_guide = {}
        duplicated = set()
        for g in guide_names:
            cg = canonical_guide_id(g)
            if cg in canonical_to_guide:
                duplicated.add(cg)
            else:
                canonical_to_guide[cg] = g
        for cg in duplicated:
            canonical_to_guide.pop(cg, None)
        for idx in meta.index[meta["guide_index"].isna()]:
            cg = canonical_guide_id(meta.at[idx, "guide_id"])
            if cg in canonical_to_guide:
                matched = canonical_to_guide[cg]
                meta.at[idx, "guide_index"] = guide_name_to_idx[matched]
                meta.at[idx, "match_method"] = "canonical_guide_id_match"

    mapped = meta[meta["guide_index"].notna()].copy()
    mapped["guide_index"] = mapped["guide_index"].astype(int)
    mapped["guide_name"] = [guide_names[i] for i in mapped["guide_index"]]

    records = []
    for _, row in mapped.iterrows():
        targets = clean_target_string(row["target_entry"])
        if looks_control(row["guide_name"]) or looks_control(row["target_entry"]):
            targets = [CONTROL_TARGET]
        for target in targets:
            records.append({
                "guide_index": int(row["guide_index"]),
                "guide_name": str(row["guide_name"]),
                "target_name": str(target),
                "raw_target_entry": row["target_entry"],
                "match_method": row["match_method"],
            })
    guide_target = pd.DataFrame(records).drop_duplicates()
    guide_target.to_csv(OUTDIR / "guide_to_target_map_from_300genes_metadata.csv", index=False)

    print("Guide metadata mapping QC")
    print(f"adata guide columns: {len(guide_names):,}")
    print(f"metadata guide rows: {len(meta):,}")
    print(f"metadata rows matched to adata: {mapped['guide_id'].nunique():,}")
    print(f"expanded guide-target rows: {len(guide_target):,}")
    print(f"unique mapped targets: {guide_target['target_name'].nunique():,}")
    return guide_target


def build_target_assignment(adata, guide_target):
    guide_assignment = get_guide_assignment_matrix(adata)
    if guide_assignment.shape[1] != len(adata.uns[cfg.GUIDE_NAMES_KEY]):
        raise ValueError("guide_assignment column count does not match guide_names length")

    series_by_target = {}
    target_to_guide_indices = {}
    for target, sub in guide_target.groupby("target_name"):
        idx = sorted(set(sub["guide_index"].astype(int).tolist()))
        target_to_guide_indices[target] = idx
        series_by_target[target] = any_guides_per_cell(guide_assignment, idx)

    target_assignment = pd.DataFrame(series_by_target, index=adata.obs_names).astype(bool)

    if CONTROL_TARGET not in target_assignment.columns:
        control_like = [c for c in target_assignment.columns if looks_control(c)]
        if len(control_like) == 0:
            raise ValueError("No non-targeting/control guides found")
        target_assignment[CONTROL_TARGET] = target_assignment[control_like].any(axis=1)

    carrier_counts = pd.DataFrame({
        "target": target_assignment.columns,
        "n_carrier_cells": target_assignment.sum(axis=0).astype(int).values,
        "n_guides": [len(target_to_guide_indices.get(t, [])) for t in target_assignment.columns],
    }).sort_values("n_carrier_cells", ascending=False)
    carrier_counts.to_csv(OUTDIR / "target_carrier_counts.csv", index=False)
    print("Control carrier cells:", int(target_assignment[CONTROL_TARGET].sum()))
    return target_assignment, target_to_guide_indices

```

## Markdown cell 7

## 4. Helper functions: cNMF KNN representation and base inputs

## Code cell 8

```python

def load_or_make_cnmf_usage(adata):
    if cfg.CNMF_KEY in adata.obsm:
        usage = adata.obsm[cfg.CNMF_KEY]
        if isinstance(usage, pd.DataFrame):
            df = usage.copy()
        else:
            df = pd.DataFrame(np.asarray(usage), index=adata.obs_names)
            df.columns = [f"Program {i}" for i in range(df.shape[1])]
        return df

    print("Loading cNMF usage:", cfg.CNMF_USAGE_TXT)
    usage_df = pd.read_csv(cfg.CNMF_USAGE_TXT, sep="\t", index_col=0)
    usage_aligned = usage_df.reindex(adata.obs_names).fillna(0)
    usage_aligned.columns = [f"Program {i}" for i in range(usage_aligned.shape[1])]
    adata.obsm[cfg.CNMF_KEY] = usage_aligned
    return usage_aligned


def make_cnmf_knn_matrix(program_df):
    if cfg.PROGRAMS_FOR_KNN is None:
        cols = list(program_df.columns)
    else:
        cols = []
        for p in cfg.PROGRAMS_FOR_KNN:
            col = f"Program {int(p)}"
            if col not in program_df.columns:
                raise KeyError(f"{col} not found in cNMF usage columns")
            cols.append(col)
    X = program_df[cols].values.astype(float)
    if cfg.CNMF_TRANSFORM == "zscore":
        X = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
    elif cfg.CNMF_TRANSFORM == "clr_zscore":
        X = np.log1p(X / np.maximum(X.sum(axis=1, keepdims=True), 1e-12))
        X = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
    elif cfg.CNMF_TRANSFORM in [None, "none"]:
        pass
    else:
        raise ValueError(f"Unknown CNMF_TRANSFORM={cfg.CNMF_TRANSFORM}")
    pd.DataFrame({"program_column_used": cols}).to_csv(OUTDIR / "cnmf_programs_used_for_knn.csv", index=False)
    print(f"Using {len(cols)} cNMF programs for KNN; transform={cfg.CNMF_TRANSFORM}")
    return X, cols


def prepare_base_inputs(adata):
    adata_use = adata.copy()
    if not adata_use.obs_names.is_unique:
        adata_use.obs_names_make_unique()

    adata_use.obs["day_clean"] = adata_use.obs[cfg.DAY_COL].map(clean_day)
    if cfg.CELLTYPE_COL in adata_use.obs.columns:
        adata_use.obs["cell_type_clean"] = adata_use.obs[cfg.CELLTYPE_COL].astype(str)
    else:
        adata_use.obs["cell_type_clean"] = "unknown"

    guide_target = load_guide_target_from_metadata(adata_use)
    target_assignment, _ = build_target_assignment(adata_use, guide_target)

    guide_assignment = get_guide_assignment_matrix(adata_use)
    adata_use.obs["n_total_guides"] = (
        np.asarray(guide_assignment.sum(axis=1)).ravel()
        if sparse.issparse(guide_assignment)
        else guide_assignment.sum(axis=1)
    ).astype(int)
    targeting_cols = [c for c in target_assignment.columns if c != CONTROL_TARGET]
    adata_use.obs["n_targeting_targets"] = target_assignment[targeting_cols].sum(axis=1).astype(int).values

    program_df = load_or_make_cnmf_usage(adata_use)
    X_knn, cols = make_cnmf_knn_matrix(program_df)

    if cfg.UMAP_KEY not in adata_use.obsm:
        raise KeyError(f"adata.obsm[{cfg.UMAP_KEY!r}] not found")
    umap = np.asarray(adata_use.obsm[cfg.UMAP_KEY])
    adata_use.obs["UMAP1"] = umap[:, 0]
    adata_use.obs["UMAP2"] = umap[:, 1]

    if cfg.DAYS_TO_USE is not None:
        days_clean = [clean_day(d) for d in cfg.DAYS_TO_USE]
        day_mask = adata_use.obs["day_clean"].isin(days_clean).values
    else:
        day_mask = np.ones(adata_use.n_obs, dtype=bool)

    background = adata_use.obs[["UMAP1", "UMAP2", "day_clean", "cell_type_clean", "n_total_guides", "n_targeting_targets"]].copy()
    background.index.name = "cell_id"
    background.reset_index().to_csv(OUTDIR / "all_days_background_cells_for_plotting.csv", index=False)
    return adata_use, target_assignment, X_knn, day_mask, program_df

```

## Markdown cell 9

## 5. Helper functions: carrier-seeded neighborhood scoring

## Code cell 10

```python

def get_masks_for_target(adata, target_assignment, target, day_mask):
    if target not in target_assignment.columns:
        matches = [c for c in target_assignment.columns if c.upper() == target.upper()]
        if not matches:
            raise KeyError(f"Target {target} not found in target_assignment")
        target = matches[0]
    target_raw = target_assignment[target].values.astype(bool)
    nt_raw = target_assignment[CONTROL_TARGET].values.astype(bool)

    if cfg.MODE == "all_carriers":
        target_mask = target_raw & day_mask
        nt_mask = nt_raw & (~target_raw) & day_mask
    elif cfg.MODE == "low_guide_load_leq2":
        low = adata.obs["n_total_guides"].values <= 2
        target_mask = target_raw & low & day_mask
        nt_mask = nt_raw & (~target_raw) & low & day_mask
    elif cfg.MODE == "single_selected_target":
        n_targeting = adata.obs["n_targeting_targets"].values
        target_mask = target_raw & (n_targeting == 1) & day_mask
        nt_mask = nt_raw & (~target_raw) & (n_targeting == 0) & day_mask
    else:
        raise ValueError(f"Unknown MODE={cfg.MODE}")
    return target, target_mask, nt_mask


def fdr_bh(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() > 0:
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def score_one_seed_neighborhood(adata, target, seed_cell_idx, neigh, target_mask, nt_mask, day_vec, cell_type_vec):
    seed_day = day_vec[seed_cell_idx]
    bg_mask = (day_vec == seed_day) if cfg.SAME_DAY_BASELINE else np.ones(len(day_vec), dtype=bool)
    bg_target_total = max(int(np.sum(target_mask & bg_mask)) - 1, 0)
    bg_nt_total = int(np.sum(nt_mask & bg_mask))

    local_target_n = int(np.sum(target_mask[neigh]))
    local_nt_n = int(np.sum(nt_mask[neigh]))

    local_ratio = (local_target_n + cfg.PSEUDOCOUNT) / (local_nt_n + cfg.PSEUDOCOUNT)
    bg_ratio = (bg_target_total + cfg.PSEUDOCOUNT) / (bg_nt_total + cfg.PSEUDOCOUNT)
    log2_local_enrichment = np.log2(local_ratio / bg_ratio)

    bg_target_outside = max(bg_target_total - local_target_n, 0)
    bg_nt_outside = max(bg_nt_total - local_nt_n, 0)
    try:
        fisher_p = fisher_exact([[local_target_n, local_nt_n], [bg_target_outside, bg_nt_outside]], alternative="two-sided")[1]
    except Exception:
        fisher_p = np.nan

    n_target_nt = local_target_n + local_nt_n
    return {
        "target": target,
        "seed_cell_idx": int(seed_cell_idx),
        "seed_cell_id": str(adata.obs_names[seed_cell_idx]),
        "day": seed_day,
        "cell_type": cell_type_vec[seed_cell_idx],
        "mode": cfg.MODE,
        "k_requested": int(cfg.K),
        "k_used": int(len(neigh)),
        "restrict_neighbors_to_same_day": bool(cfg.RESTRICT_NEIGHBORS_TO_SAME_DAY),
        "same_day_baseline": bool(cfg.SAME_DAY_BASELINE),
        "local_target_n": local_target_n,
        "local_nt_n": local_nt_n,
        "local_target_plus_nt_n": n_target_nt,
        "background_target_total": int(bg_target_total),
        "background_nt_total": int(bg_nt_total),
        "local_target_fraction_among_target_nt": local_target_n / n_target_nt if n_target_nt > 0 else np.nan,
        "background_target_fraction_among_target_nt": bg_target_total / (bg_target_total + bg_nt_total) if (bg_target_total + bg_nt_total) > 0 else np.nan,
        "log2_local_enrichment": float(log2_local_enrichment),
        "fisher_p": float(fisher_p) if np.isfinite(fisher_p) else np.nan,
        "UMAP1": float(adata.obs["UMAP1"].iloc[seed_cell_idx]),
        "UMAP2": float(adata.obs["UMAP2"].iloc[seed_cell_idx]),
        "n_total_guides": int(adata.obs["n_total_guides"].iloc[seed_cell_idx]),
        "n_targeting_targets": int(adata.obs["n_targeting_targets"].iloc[seed_cell_idx]),
    }


def compute_carrier_seeded_enrichment(adata, target_assignment, X_knn, target, day_mask):
    resolved_target, target_mask, nt_mask = get_masks_for_target(adata, target_assignment, target, day_mask)
    seed_idx_all = np.where(target_mask)[0]
    if len(seed_idx_all) == 0:
        print(f"[skip] {target}: no carrier seeds")
        return pd.DataFrame()

    day_vec = adata.obs["day_clean"].astype(str).values
    cell_type_vec = adata.obs["cell_type_clean"].astype(str).values
    rows = []

    if cfg.RESTRICT_NEIGHBORS_TO_SAME_DAY:
        for day in sorted(pd.Series(day_vec[seed_idx_all]).unique()):
            day_cell_idx = np.where((day_vec == day) & day_mask)[0]
            seed_idx = seed_idx_all[day_vec[seed_idx_all] == day]
            if len(day_cell_idx) <= 1 or len(seed_idx) == 0:
                continue
            local_pos = {cell_idx: pos for pos, cell_idx in enumerate(day_cell_idx)}
            seed_local = np.array([local_pos[i] for i in seed_idx])
            X_day = X_knn[day_cell_idx]
            k_query = min(cfg.K + 1, len(day_cell_idx))
            nbrs = NearestNeighbors(n_neighbors=k_query, metric="euclidean")
            nbrs.fit(X_day)
            neigh_local = nbrs.kneighbors(X_day[seed_local], return_distance=False)
            for seed_cell_idx, neigh_l in zip(seed_idx, neigh_local):
                neigh = day_cell_idx[neigh_l]
                neigh = neigh[neigh != seed_cell_idx][:cfg.K]
                rows.append(score_one_seed_neighborhood(adata, resolved_target, seed_cell_idx, neigh, target_mask, nt_mask, day_vec, cell_type_vec))
    else:
        selected_idx = np.where(day_mask)[0]
        pos = {cell_idx: i for i, cell_idx in enumerate(selected_idx)}
        seed_idx = np.array([i for i in seed_idx_all if i in pos])
        seed_local = np.array([pos[i] for i in seed_idx])
        X_sel = X_knn[selected_idx]
        k_query = min(cfg.K + 1, len(selected_idx))
        nbrs = NearestNeighbors(n_neighbors=k_query, metric="euclidean")
        nbrs.fit(X_sel)
        neigh_local = nbrs.kneighbors(X_sel[seed_local], return_distance=False)
        for seed_cell_idx, neigh_l in zip(seed_idx, neigh_local):
            neigh = selected_idx[neigh_l]
            neigh = neigh[neigh != seed_cell_idx][:cfg.K]
            rows.append(score_one_seed_neighborhood(adata, resolved_target, seed_cell_idx, neigh, target_mask, nt_mask, day_vec, cell_type_vec))

    out = pd.DataFrame(rows)
    if len(out) == 0:
        return out
    out["fisher_fdr"] = fdr_bh(out["fisher_p"].values)
    out["fisher_fdr_by_day"] = np.nan
    for day, sub in out.groupby("day"):
        out.loc[sub.index, "fisher_fdr_by_day"] = fdr_bh(sub["fisher_p"].values)
    out["minus_log10_fdr"] = -np.log10(out["fisher_fdr"].clip(lower=1e-300))
    out["minus_log10_fdr_by_day"] = -np.log10(out["fisher_fdr_by_day"].clip(lower=1e-300))
    out["significant_fdr_0_05"] = out["fisher_fdr"] < 0.05
    out["significant_fdr_by_day_0_05"] = out["fisher_fdr_by_day"] < 0.05
    return out


def load_or_compute_scores(adata_use, target_assignment, X_knn, targets, day_mask):
    frames = []
    for target in targets:
        out_file = OUTDIR / f"carrier_seeded_scores_{target}_{cfg.MODE}.csv"
        if out_file.exists():
            print("Loading existing:", out_file)
            scores = pd.read_csv(out_file)
        else:
            print("Computing carrier-seeded enrichment:", target)
            scores = compute_carrier_seeded_enrichment(adata_use, target_assignment, X_knn, target, day_mask)
            if len(scores):
                scores.to_csv(out_file, index=False)
        if len(scores):
            frames.append(scores)
    if not frames:
        raise RuntimeError("No carrier-seeded scores produced")
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(OUTDIR / f"combined_carrier_seeded_scores_{cfg.MODE}.csv", index=False)

    summary_day = combined.groupby(["target", "day"]).agg(
        n_seed_cells=("seed_cell_id", "count"),
        median_log2_enrichment=("log2_local_enrichment", "median"),
        mean_log2_enrichment=("log2_local_enrichment", "mean"),
        frac_positive=("log2_local_enrichment", lambda x: np.mean(x > 0)),
        n_fdr05=("fisher_fdr", lambda x: int(np.sum(x < 0.05))),
        n_fdr_by_day05=("fisher_fdr_by_day", lambda x: int(np.sum(x < 0.05))),
        median_k_used=("k_used", "median"),
        median_local_target_n=("local_target_n", "median"),
        median_local_nt_n=("local_nt_n", "median"),
    ).reset_index()
    summary_day.to_csv(OUTDIR / f"summary_by_target_day_{cfg.MODE}.csv", index=False)
    return combined, summary_day

```

## Markdown cell 11

## 6. Helper functions: plotting and inline display

## Code cell 12

```python

def _geom_point_optional_raster(data, mapping, raster=False, **kwargs):
    """
    Helper for plotnine versions that may or may not support raster=True.
    If raster=True is unsupported, it falls back to vector points.
    """
    if raster:
        try:
            return geom_point(data, mapping, raster=True, **kwargs)
        except TypeError:
            return geom_point(data, mapping, **kwargs)
    return geom_point(data, mapping, **kwargs)


def prepare_background_for_plotting(bg):
    """Stratified background subsampling by day, matching previous notebook style."""
    bg = bg.copy()
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    if len(bg) > cfg.MAX_BACKGROUND_CELLS:
        bg_parts = []
        if "day_clean" in bg.columns:
            n_days = max(1, bg["day_clean"].nunique())
            n_per_day = max(1, cfg.MAX_BACKGROUND_CELLS // n_days)
            for day, sub in bg.groupby("day_clean"):
                n = min(len(sub), n_per_day)
                idx = rng.choice(sub.index.values, size=n, replace=False) if len(sub) > n else sub.index.values
                bg_parts.append(sub.loc[idx])
            bg_plot = pd.concat(bg_parts, ignore_index=True)
        else:
            idx = rng.choice(bg.index.values, size=cfg.MAX_BACKGROUND_CELLS, replace=False)
            bg_plot = bg.loc[idx].copy()
    else:
        bg_plot = bg.copy()
    return bg_plot


def get_shared_limits(bg, scores):
    xmin = min(bg["UMAP1"].min(), scores["UMAP1"].min())
    xmax = max(bg["UMAP1"].max(), scores["UMAP1"].max())
    ymin = min(bg["UMAP2"].min(), scores["UMAP2"].min())
    ymax = max(bg["UMAP2"].max(), scores["UMAP2"].max())
    pad_x = (xmax - xmin) * 0.03
    pad_y = (ymax - ymin) * 0.03
    return (xmin - pad_x, xmax + pad_x), (ymin - pad_y, ymax + pad_y)


def make_plotnine_neighborhood_facet(scores, bg, targets, title):
    """
    Plotnine faceted UMAP style copied from the previous Hippo notebook.

    - gray all-day background cells
    - carrier seed cells colored by local log2OR
    - shared UMAP coordinates and shared color scale
    - no dot size/outline encoding, to keep small panels clean
    """
    bg_plot = prepare_background_for_plotting(bg)
    scores_plot = scores[scores["target"].astype(str).str.upper().isin([t.upper() for t in targets])].copy()
    scores_plot["target"] = pd.Categorical(scores_plot["target"], categories=targets, ordered=True)
    scores_plot["log2OR_plot"] = scores_plot["log2_local_enrichment"].clip(-cfg.PLOT_CLIP, cfg.PLOT_CLIP)

    XLIM, YLIM = get_shared_limits(bg, scores_plot)

    p = (
        ggplot()
        + _geom_point_optional_raster(
            bg_plot,
            aes("UMAP1", "UMAP2"),
            raster=cfg.RASTERIZE_BACKGROUND_FOR_SPEED,
            color="#D0D0D0",
            size=cfg.BACKGROUND_POINT_SIZE,
            alpha=cfg.BACKGROUND_ALPHA,
            show_legend=False,
        )
        + geom_point(
            scores_plot,
            aes("UMAP1", "UMAP2", color="log2OR_plot"),
            size=cfg.FOREGROUND_POINT_SIZE,
            alpha=0.95,
        )
        + facet_wrap("~target", ncol=cfg.FACET_NCOL)
        + scale_color_gradient2(
            low=cfg.CMAP_LOW,
            mid=cfg.CMAP_MID,
            high=cfg.CMAP_HIGH,
            midpoint=0,
            limits=(-cfg.PLOT_CLIP, cfg.PLOT_CLIP),
            breaks=[-2, -1, 0, 1, 2],
            name="Carrier-seeded\nenrichment\nlog2OR",
        )
        + coord_equal(xlim=XLIM, ylim=YLIM, expand=False)
        + labs(
            title=title,
            subtitle="All D0-D3 cells in one shared UMAP; color = local enrichment vs non-targeting controls",
            x="single-cell state embedding 1",
            y="single-cell state embedding 2",
        )
        + theme_bw(base_size=12, base_family=AI_FONT)
        + theme(
            figure_size=cfg.FIGSIZE_GRID,
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(color="#7f7f7f", fill=None, size=0.8),
            axis_title=element_text(size=13, family=AI_FONT),
            axis_text=element_text(size=10, family=AI_FONT),
            plot_title=element_text(size=17, weight="bold", family=AI_FONT),
            plot_subtitle=element_text(size=10.5, family=AI_FONT),
            strip_background=element_rect(fill="#D9D9D9", color="#7f7f7f", size=0.8),
            strip_text=element_text(size=14, weight="bold", family=AI_FONT),
            legend_title=element_text(size=11, weight="bold", family=AI_FONT),
            legend_text=element_text(size=10, family=AI_FONT),
            legend_key=element_blank(),
        )
    )
    return p


def save_plotnine_panel(p, out_prefix):
    saved = {}
    width, height = cfg.FIGSIZE_GRID
    for ext in ["png", "pdf", "svg"]:
        path = OUTDIR / f"{out_prefix}.{ext}"
        if ext == "png":
            p.save(path, width=width, height=height, dpi=450, verbose=False)
        else:
            p.save(path, width=width, height=height, verbose=False)
        saved[ext] = path
    return saved


def show_png(path, width=900):
    path = Path(path)
    if path.exists():
        display(Image(filename=str(path), width=width))
    else:
        print("Missing:", path)

```

## Markdown cell 13

## 7. Load AnnData and prepare guide/cNMF inputs

This is the most time-consuming setup cell. It writes QC files into `OUTDIR`.

## Code cell 14

```python

print("Reading AnnData:", cfg.ADATA_H5AD)
adata = sc.read_h5ad(cfg.ADATA_H5AD)
print(adata)

adata_use, target_assignment, X_knn, day_mask, program_df = prepare_base_inputs(adata)
print("Prepared adata_use:", adata_use)
print("target_assignment:", target_assignment.shape)
print("X_knn:", X_knn.shape)

# quick QC display
carrier_counts = pd.read_csv(OUTDIR / "target_carrier_counts.csv")
display(carrier_counts.head(20))

```

## Markdown cell 15

## 8. Compute or load carrier-seeded enrichment scores

This cell computes scores for panel f and panel j targets. If the per-target CSVs already exist in `OUTDIR`, they are loaded instead of recomputed.

## Code cell 16

```python

all_targets = list(dict.fromkeys(cfg.PANEL_F_TARGETS + cfg.PANEL_J_TARGETS))
combined_scores, summary_day = load_or_compute_scores(adata_use, target_assignment, X_knn, all_targets, day_mask)

print("Combined scores:", combined_scores.shape)
display(summary_day)

```

## Markdown cell 17

## 9. Plot panel f and display inline

## Code cell 18

```python

bg = pd.read_csv(OUTDIR / "all_days_background_cells_for_plotting.csv")

p_panel_f = make_plotnine_neighborhood_facet(
    combined_scores,
    bg,
    cfg.PANEL_F_TARGETS,
    "Panel f: carrier-seeded CRISPRi neighborhood enrichment"
)

panel_f_paths = save_plotnine_panel(
    p_panel_f,
    "panel_f_carrierseeded_neighborhood_grid_plotnine_previous_style"
)
print(panel_f_paths)
show_png(panel_f_paths["png"], width=900)

# Also display as a plotnine object inline in notebooks that support it.
p_panel_f

```

## Markdown cell 19

## 10. Plot panel j and display inline

## Code cell 20

```python

p_panel_j = make_plotnine_neighborhood_facet(
    combined_scores,
    bg,
    cfg.PANEL_J_TARGETS,
    "Panel j: carrier-seeded CRISPRi neighborhood enrichment"
)

panel_j_paths = save_plotnine_panel(
    p_panel_j,
    "panel_j_carrierseeded_neighborhood_grid_plotnine_previous_style"
)
print(panel_j_paths)
show_png(panel_j_paths["png"], width=900)

# Also display as a plotnine object inline in notebooks that support it.
p_panel_j

```

## Markdown cell 21

## 11. Optional: inspect one target/day table

Use this to inspect whether enrichment is concentrated at D2/D3 or another day.

## Code cell 22

```python

# Edit target below interactively
TARGET_TO_INSPECT = "NF2"  # or ACVR2A, BMPR1A, SMAD1, SMAD5, KDR, HAND1, PRKCD, PLCG1

if TARGET_TO_INSPECT in combined_scores["target"].unique():
    display(
        combined_scores.query("target == @TARGET_TO_INSPECT")
        .groupby(["day", "cell_type"])
        .agg(
            n_seed_cells=("seed_cell_id", "count"),
            median_log2_enrichment=("log2_local_enrichment", "median"),
            mean_log2_enrichment=("log2_local_enrichment", "mean"),
            n_fdr_by_day05=("fisher_fdr_by_day", lambda x: int(np.sum(x < 0.05)))
        )
        .reset_index()
        .sort_values(["day", "median_log2_enrichment"], ascending=[True, False])
    )
else:
    print(f"{TARGET_TO_INSPECT} not in combined_scores. Available targets:")
    print(sorted(combined_scores["target"].unique()))

```

## Markdown cell 23

## 12. Package outputs

## Code cell 24

```python

zip_path = Path(f"{cfg.OUTDIR}.zip")
with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
    for fn in OUTDIR.iterdir():
        if fn.is_file():
            z.write(fn, arcname=fn.name)
print("Saved outputs to", OUTDIR.resolve())
print("ZIP:", zip_path.resolve())

```

