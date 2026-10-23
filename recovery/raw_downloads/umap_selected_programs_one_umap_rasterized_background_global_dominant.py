# %% [markdown]
# # One UMAP for selected programs with rasterized background
#
# This script makes a single UMAP colored only for cells whose **globally dominant**
# cNMF program is one of:
# - P13
# - P30
# - P15
# - P17
#
# All other cells are shown as background.
#
# Background points are rasterized to keep the SVG file size small.
# Foreground highlighted points remain vector, and text stays editable in Illustrator.

# %%
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

# %%
def find_repo_root(start: Path = Path.cwd()) -> Path:
    start = start.resolve()
    required_paths = [Path("data/program_activity_long_v4.csv")]
    for candidate in [start, *start.parents]:
        if all((candidate / rel_path).exists() for rel_path in required_paths):
            return candidate
    raise FileNotFoundError("Could not auto-detect REPO_ROOT. Set REPO_ROOT manually.")

REPO_ROOT = find_repo_root()
OUT_DIR = REPO_ROOT / "analysis" / "regulator_dotplots" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FONT_DIR = REPO_ROOT / ".pixi" / "envs" / "default" / "fonts"
if FONT_DIR.exists():
    for font_file in fm.findSystemFonts(fontpaths=[str(FONT_DIR)]):
        fm.fontManager.addfont(font_file)

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

print(f"REPO_ROOT = {REPO_ROOT}")
print(f"OUT_DIR   = {OUT_DIR}")

# %%
def align_umap_to_cnmf_cells(*, cnmf, adata_int, obsm_key="X_umap", require_exact_match=True):
    cnmf_obs = pd.Index(cnmf.obs_names)
    adata_obs = pd.Index(adata_int.obs_names)

    if obsm_key not in adata_int.obsm_keys():
        raise KeyError(f"{obsm_key!r} not found in adata_int.obsm")

    missing_in_adata = cnmf_obs[~cnmf_obs.isin(adata_obs)]
    extra_in_adata = adata_obs[~adata_obs.isin(cnmf_obs)]
    if len(missing_in_adata) > 0:
        raise ValueError(f"Some cNMF cells are missing from adata_int: {missing_in_adata[:10].tolist()}")
    if require_exact_match and len(extra_in_adata) > 0:
        raise ValueError(f"adata_int contains extra cells not in cNMF: {extra_in_adata[:10].tolist()}")

    adata_umap = np.asarray(adata_int.obsm[obsm_key])
    if cnmf_obs.equals(adata_obs):
        aligned_umap = adata_umap[:, :2]
    else:
        positions = adata_obs.get_indexer(cnmf_obs)
        aligned_umap = adata_umap[positions, :2]

    print(f"Transferred {obsm_key}: {aligned_umap.shape[0]} cells x {aligned_umap.shape[1]} dims")
    return aligned_umap


def parse_program_id(x):
    s = str(x)
    return int(s[1:]) if s.startswith("P") else int(s)

# %%
H5MU_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/adata/cNMF_50_0_2.h5mu"
)
ADATA_INT_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
    "2.preprocessing_tonys_pipeline/3.integrated_days/"
    "3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad"
)
PROGRAM_META_PATH = REPO_ROOT / "data" / "program_activity_long_v4.csv"
# PROGRAM_META_PATH = Path("/mnt/data/program_activity_long_v4_updated_colors_labels.csv")

SELECTED_PROGRAMS = [13, 30, 15, 17]
THRESHOLD = 0.02

BACKGROUND_COLOR = "#C8B6FF"  # lavender-purple, distinct from revised program palette
BACKGROUND_ALPHA = 0.23
BACKGROUND_SIZE = 0.22
FOREGROUND_ALPHA = 0.95
FOREGROUND_SIZE = 0.7
FOREGROUND_EDGEWIDTH = 0.0
FIG_W = 4.2
FIG_H = 4.2
OUTPUT_STEM = OUT_DIR / "umap_P13_P30_P15_P17_one_umap_global_dominant_rasterized_background"

# %%
print("Loading h5mu...")
mdata = mu.read_h5mu(H5MU_PATH)
cnmf = mdata.mod["cNMF"]
print(f"Loaded cNMF modality: {cnmf.n_obs} cells x {cnmf.n_vars} programs")

print("Loading integrated AnnData for UMAP...")
adata_int = sc.read_h5ad(ADATA_INT_PATH, backed="r")
print(f"Loaded integrated AnnData: {adata_int.n_obs} cells x {adata_int.n_vars} genes")

umap_coords = align_umap_to_cnmf_cells(cnmf=cnmf, adata_int=adata_int, obsm_key="X_umap", require_exact_match=True)
scores = cnmf.X.toarray() if hasattr(cnmf.X, "toarray") else np.asarray(cnmf.X)
print(f"scores shape = {scores.shape}")
print(f"umap_coords shape = {umap_coords.shape}")

# %%
prog_meta = pd.read_csv(PROGRAM_META_PATH).drop_duplicates("program_name")
prog_meta["program_name"] = prog_meta["program_name"].astype(int)
label_map = dict(zip(prog_meta["program_name"], prog_meta.get("program_label", prog_meta["program_name"].astype(str)).astype(str)))

preferred_color_columns = ["distinct_color", "program_color", "color", "class_color"]
color_column = next((c for c in preferred_color_columns if c in prog_meta.columns), None)
if color_column is None:
    raise ValueError(f"No usable color column found. Available columns: {list(prog_meta.columns)}")
prog_color = {}
for program, color in zip(prog_meta["program_name"], prog_meta[color_column]):
    if pd.isna(color):
        continue
    color = str(color).strip()
    try:
        mpl.colors.to_rgba(color)
        prog_color[int(program)] = color
    except Exception:
        pass

cnmf_program_ids = [parse_program_id(v) for v in cnmf.var_names] if len(cnmf.var_names) == scores.shape[1] else list(range(scores.shape[1]))
program_to_col = {pid: idx for idx, pid in enumerate(cnmf_program_ids)}

for p in SELECTED_PROGRAMS:
    if p not in program_to_col:
        raise ValueError(f"Selected program P{p} not found among cNMF columns")
    if p not in prog_color:
        raise ValueError(f"Selected program P{p} missing a valid color in metadata")
    label_map.setdefault(p, f"P{p}")

print("Selected programs:")
for p in SELECTED_PROGRAMS:
    print(f"  P{p}: column {program_to_col[p]}, label={label_map[p]}, color={prog_color[p]}")

# %%
# Global dominant program across all cNMF programs.
dominant_col_idx = np.argmax(scores, axis=1)
dominant_score = scores[np.arange(scores.shape[0]), dominant_col_idx]
dominant_program_all = np.array([cnmf_program_ids[idx] for idx in dominant_col_idx], dtype=int)

# Only color cells whose globally dominant program is one of the selected programs.
selected_mask = np.isin(dominant_program_all, SELECTED_PROGRAMS) & (dominant_score > THRESHOLD)

summary_records = []
for p in SELECTED_PROGRAMS:
    mask = selected_mask & (dominant_program_all == p)
    summary_records.append({
        "program": p,
        "program_label": label_map[p],
        "hex_color": prog_color[p],
        "n_cells": int(mask.sum()),
    })
summary_df = pd.DataFrame(summary_records)
summary_df.to_csv(OUTPUT_STEM.with_suffix('.summary.csv'), index=False)
print(summary_df)
print("selected cells:", int(selected_mask.sum()))
print("background cells:", int((~selected_mask).sum()))

# %%
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), constrained_layout=True)
x = umap_coords[:, 0]
y = umap_coords[:, 1]

# Background first, rasterized.
bg_mask = ~selected_mask
ax.scatter(
    x[bg_mask], y[bg_mask],
    s=BACKGROUND_SIZE, c=BACKGROUND_COLOR, alpha=BACKGROUND_ALPHA,
    linewidths=0, edgecolors="none", rasterized=True, zorder=1,
)

# Foreground selected programs, vector.
for p in SELECTED_PROGRAMS:
    mask = selected_mask & (dominant_program_all == p)
    ax.scatter(
        x[mask], y[mask],
        s=FOREGROUND_SIZE, c=prog_color[p], alpha=FOREGROUND_ALPHA,
        linewidths=FOREGROUND_EDGEWIDTH, edgecolors="none", rasterized=False,
        zorder=2, label=f"P{p}: {label_map[p]} (n={int(mask.sum()):,})",
    )

ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect("equal")
for spine in ax.spines.values():
    spine.set_visible(False)

legend_elements = [
    Line2D([0], [0], marker='o', linestyle='', markersize=4.5,
           markerfacecolor=BACKGROUND_COLOR, markeredgewidth=0, alpha=0.75,
           label=f"Background (n={int(bg_mask.sum()):,})")
]
for p in SELECTED_PROGRAMS:
    mask = selected_mask & (dominant_program_all == p)
    legend_elements.append(
        Line2D([0], [0], marker='o', linestyle='', markersize=4.5,
               markerfacecolor=prog_color[p], markeredgewidth=0,
               label=f"P{p}: {label_map[p]} (n={int(mask.sum()):,})")
    )

ax.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7, handletextpad=0.5, borderaxespad=0.0)
ax.set_title(
    f"Integrated UMAP colored by globally dominant selected program\n"
    f"Only cells whose global dominant program is in {', '.join([f'P{p}' for p in SELECTED_PROGRAMS])} are colored; threshold > {THRESHOLD}",
    fontsize=8, pad=6,
)

for ext in ["svg", "pdf", "png"]:
    fig.savefig(OUTPUT_STEM.with_suffix(f".{ext}"), dpi=300 if ext == "png" else None)

plt.show()
plt.close(fig)

print("Saved outputs:")
print(OUTPUT_STEM.with_suffix('.svg'))
print(OUTPUT_STEM.with_suffix('.pdf'))
print(OUTPUT_STEM.with_suffix('.png'))
print(OUTPUT_STEM.with_suffix('.summary.csv'))
