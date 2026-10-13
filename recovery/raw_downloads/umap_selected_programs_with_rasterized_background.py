# %% [markdown]
# # UMAP panels for selected programs with rasterized background
#
# This notebook/script makes a 2x2 UMAP panel figure for P13, P30, P15, and P17.
#
# - Background cells are shown in a single distinct background color that is intentionally
#   different from all program colors in the revised `program_activity_long_v4.csv`.
# - Background points are rasterized to keep the SVG file size small.
# - Highlighted cells remain vector objects, and text stays editable in Illustrator.
#
# Highlighting rule (default):
# - A cell is highlighted in a panel if the selected program is that cell's dominant cNMF
#   program and its score is above `THRESHOLD`.
#
# You can switch to score-based highlighting by setting:
# - `HIGHLIGHT_MODE = "score"`
#   which highlights cells with selected-program score > `THRESHOLD` regardless of dominance.

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

warnings.filterwarnings("ignore")

# %%
def find_repo_root(start: Path = Path.cwd()) -> Path:
    """Find repository root by looking upward for expected project files."""
    start = start.resolve()
    required_paths = [Path("data/program_activity_long_v4.csv")]
    for candidate in [start, *start.parents]:
        if all((candidate / rel_path).exists() for rel_path in required_paths):
            return candidate
    raise FileNotFoundError(
        "Could not auto-detect REPO_ROOT. Set REPO_ROOT manually in this cell."
    )


# In a notebook, __file__ is usually not defined, so use an auto-detected repo root.
REPO_ROOT = find_repo_root()
OUT_DIR = REPO_ROOT / "analysis" / "regulator_dotplots" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FONT_DIR = REPO_ROOT / ".pixi" / "envs" / "default" / "fonts"
if FONT_DIR.exists():
    for font_file in fm.findSystemFonts(fontpaths=[str(FONT_DIR)]):
        fm.fontManager.addfont(font_file)

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        # Keep text editable in Illustrator when saving SVG.
        "svg.fonttype": "none",
    }
)

print(f"REPO_ROOT = {REPO_ROOT}")
print(f"OUT_DIR   = {OUT_DIR}")

# %%
# -----------------------------
# Configuration
# -----------------------------

H5MU_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/adata/cNMF_50_0_2.h5mu"
)

ADATA_INT_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
    "2.preprocessing_tonys_pipeline/3.integrated_days/"
    "3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad"
)

# Use the revised metadata file if you have replaced the repo copy.
# Otherwise point this to the updated CSV manually.
PROGRAM_META_PATH = REPO_ROOT / "data" / "program_activity_long_v4.csv"
# Example manual override:
# PROGRAM_META_PATH = Path("/mnt/data/program_activity_long_v4_updated_colors_labels.csv")

SELECTED_PROGRAMS = [13, 30, 15, 17]
THRESHOLD = 0.02
HIGHLIGHT_MODE = "dominant"  # "dominant" or "score"

# Chosen to be clearly different from the revised program palette (grey/yellow/orange/red/blue/green).
# This soft purple/lavender is visually distinct and still recedes behind highlighted cells.
BACKGROUND_COLOR = "#C8B6FF"
BACKGROUND_ALPHA = 0.22
BACKGROUND_SIZE = 0.18

HIGHLIGHT_ALPHA = 0.95
HIGHLIGHT_SIZE = 1.25
HIGHLIGHT_EDGEWIDTH = 0.0

N_COLS = 2
FIG_W = 9.4
FIG_H = 8.0

OUTPUT_STEM = OUT_DIR / "umap_P13_P30_P15_P17_selected_with_rasterized_background"

# %%
def align_umap_to_cnmf_cells(
    *,
    cnmf,
    adata_int,
    obsm_key: str = "X_umap",
    require_exact_match: bool = True,
) -> np.ndarray:
    """Return adata_int.obsm[obsm_key] reordered to cnmf.obs_names."""
    cnmf_obs = pd.Index(cnmf.obs_names, name="cnmf_obs_names")
    adata_obs = pd.Index(adata_int.obs_names, name="adata_int_obs_names")

    if not cnmf_obs.is_unique:
        duplicated = cnmf_obs[cnmf_obs.duplicated()].unique()[:10].tolist()
        raise ValueError(f"cnmf.obs_names contains duplicate cell IDs, e.g. {duplicated}")

    if not adata_obs.is_unique:
        duplicated = adata_obs[adata_obs.duplicated()].unique()[:10].tolist()
        raise ValueError(f"adata_int.obs_names contains duplicate cell IDs, e.g. {duplicated}")

    if obsm_key not in adata_int.obsm_keys():
        raise KeyError(
            f"{obsm_key!r} was not found in adata_int.obsm. "
            f"Available obsm keys: {list(adata_int.obsm_keys())}"
        )

    missing_in_adata = cnmf_obs[~cnmf_obs.isin(adata_obs)]
    extra_in_adata = adata_obs[~adata_obs.isin(cnmf_obs)]

    if len(missing_in_adata) > 0:
        raise ValueError(
            "Some cNMF cells are missing from adata_int. "
            f"n_missing={len(missing_in_adata)}; first few: {missing_in_adata[:10].tolist()}"
        )

    if require_exact_match and len(extra_in_adata) > 0:
        raise ValueError(
            "adata_int contains cells that are not present in cNMF. "
            f"n_extra={len(extra_in_adata)}; first few: {extra_in_adata[:10].tolist()}"
        )

    adata_umap = np.asarray(adata_int.obsm[obsm_key])
    if adata_umap.ndim != 2 or adata_umap.shape[1] < 2:
        raise ValueError(
            f"adata_int.obsm[{obsm_key!r}] should be a 2D array with at least 2 columns; "
            f"got shape {adata_umap.shape}"
        )

    if adata_umap.shape[0] != adata_int.n_obs:
        raise ValueError(
            f"adata_int.obsm[{obsm_key!r}] has {adata_umap.shape[0]} rows, "
            f"but adata_int has {adata_int.n_obs} cells"
        )

    if cnmf_obs.equals(adata_obs):
        print("Cell indexes match in the same order; using adata_int X_umap directly.")
        aligned_umap = adata_umap[:, :2]
    else:
        print("Cell indexes match as IDs but not in the same order; reordering adata_int X_umap to cnmf.obs_names.")
        adata_positions = adata_obs.get_indexer(cnmf_obs)
        if np.any(adata_positions < 0):
            bad = cnmf_obs[np.where(adata_positions < 0)[0]][:10].tolist()
            raise ValueError(f"Unexpected unmatched cells after index check: {bad}")
        aligned_umap = adata_umap[adata_positions, :2]

    if aligned_umap.shape[0] != cnmf.n_obs:
        raise ValueError(
            f"Aligned UMAP has {aligned_umap.shape[0]} rows, but cNMF has {cnmf.n_obs} cells"
        )

    print(f"Transferred {obsm_key} from adata_int: {aligned_umap.shape[0]} cells x {aligned_umap.shape[1]} dimensions")
    return aligned_umap

# %%
def parse_program_id(x):
    s = str(x)
    if s.startswith("P"):
        s = s[1:]
    return int(s)


# %%
if not H5MU_PATH.exists():
    raise FileNotFoundError(f"Could not find h5mu file: {H5MU_PATH}")
if not ADATA_INT_PATH.exists():
    raise FileNotFoundError(f"Could not find integrated h5ad file: {ADATA_INT_PATH}")
if not PROGRAM_META_PATH.exists():
    raise FileNotFoundError(f"Could not find program metadata CSV: {PROGRAM_META_PATH}")

print("Loading h5mu...")
mdata = mu.read_h5mu(H5MU_PATH)
cnmf = mdata.mod["cNMF"]
print(f"Loaded cNMF modality: {cnmf.n_obs} cells x {cnmf.n_vars} programs")

print("Loading integrated AnnData for UMAP...")
adata_int = sc.read_h5ad(ADATA_INT_PATH, backed="r")
print(f"Loaded integrated AnnData: {adata_int.n_obs} cells x {adata_int.n_vars} genes")

umap_coords = align_umap_to_cnmf_cells(
    cnmf=cnmf,
    adata_int=adata_int,
    obsm_key="X_umap",
    require_exact_match=True,
)

scores = cnmf.X.toarray() if hasattr(cnmf.X, "toarray") else np.asarray(cnmf.X)
print(f"scores shape = {scores.shape}")
print(f"umap_coords shape = {umap_coords.shape}")

# %%
# Load revised program metadata
prog_meta = pd.read_csv(PROGRAM_META_PATH).drop_duplicates("program_name")
prog_meta["program_name"] = prog_meta["program_name"].astype(int)

label_map = {}
if "program_label" in prog_meta.columns:
    label_map = dict(zip(prog_meta["program_name"], prog_meta["program_label"].astype(str)))
else:
    label_map = {int(program): f"P{int(program)}" for program in prog_meta["program_name"]}

preferred_color_columns = [
    "distinct_color",
    "program_color",
    "color",
    "class_color",
]
color_column = next((col for col in preferred_color_columns if col in prog_meta.columns), None)
if color_column is None:
    raise ValueError(
        "No usable color column found in program metadata. "
        f"Available columns: {list(prog_meta.columns)}"
    )

prog_color = {}
for program, color in zip(prog_meta["program_name"], prog_meta[color_column]):
    if pd.isna(color):
        continue
    color = str(color).strip()
    if not color:
        continue
    try:
        mpl.colors.to_rgba(color)
    except ValueError:
        continue
    prog_color[int(program)] = color

# Validate selected programs
missing = [p for p in SELECTED_PROGRAMS if p not in label_map]
if missing:
    raise ValueError(f"Selected programs missing from metadata: {missing}")
missing_color = [p for p in SELECTED_PROGRAMS if p not in prog_color]
if missing_color:
    raise ValueError(f"Selected programs missing valid colors: {missing_color}")

# %%
# Determine how cNMF columns map to program ids.
# Preferred: use cnmf.var_names as program ids. Fallback: assume columns are already in 0..n-1 order.
var_name_parse_ok = True
try:
    cnmf_program_ids = [parse_program_id(v) for v in cnmf.var_names]
except Exception:
    var_name_parse_ok = False
    cnmf_program_ids = list(range(scores.shape[1]))

if len(cnmf_program_ids) != scores.shape[1]:
    raise ValueError("Could not determine cNMF program ids for score matrix columns.")

program_to_col = {pid: idx for idx, pid in enumerate(cnmf_program_ids)}
for p in SELECTED_PROGRAMS:
    if p not in program_to_col:
        raise ValueError(
            f"Selected program P{p} was not found among cNMF score columns. "
            f"First few available ids: {cnmf_program_ids[:10]}"
        )

print("Selected programs:")
for p in SELECTED_PROGRAMS:
    print(f"  P{p}: column {program_to_col[p]}, label={label_map[p]}, color={prog_color[p]}")

# %%
# Determine dominant program over the full cNMF score matrix.
dominant_col_idx = np.argmax(scores, axis=1)
dominant_score = scores[np.arange(scores.shape[0]), dominant_col_idx]
dominant_program = np.array([cnmf_program_ids[idx] for idx in dominant_col_idx], dtype=int)

summary_records = []
panel_data = {}
for p in SELECTED_PROGRAMS:
    col_idx = program_to_col[p]
    target_score = scores[:, col_idx]

    if HIGHLIGHT_MODE == "dominant":
        highlight_mask = (dominant_program == p) & (dominant_score > THRESHOLD)
    elif HIGHLIGHT_MODE == "score":
        highlight_mask = target_score > THRESHOLD
    else:
        raise ValueError("HIGHLIGHT_MODE must be 'dominant' or 'score'")

    background_mask = ~highlight_mask

    panel_data[p] = {
        "score": target_score,
        "highlight_mask": highlight_mask,
        "background_mask": background_mask,
        "n_highlight": int(highlight_mask.sum()),
    }

    summary_records.append(
        {
            "program": p,
            "program_label": label_map[p],
            "hex_color": prog_color[p],
            "highlight_mode": HIGHLIGHT_MODE,
            "threshold": THRESHOLD,
            "n_highlight": int(highlight_mask.sum()),
            "n_background": int(background_mask.sum()),
        }
    )

summary_df = pd.DataFrame(summary_records)
print(summary_df)
summary_df.to_csv(OUTPUT_STEM.with_suffix(".summary.csv"), index=False)

# %%
# Build multi-panel figure.
n_panels = len(SELECTED_PROGRAMS)
n_cols = N_COLS
n_rows = int(np.ceil(n_panels / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(FIG_W, FIG_H), constrained_layout=True)
axes = np.array(axes).reshape(-1)

# UMAP limits shared across panels
x = umap_coords[:, 0]
y = umap_coords[:, 1]
x_pad = 0.02 * (x.max() - x.min())
y_pad = 0.02 * (y.max() - y.min())

for ax_idx, program in enumerate(SELECTED_PROGRAMS):
    ax = axes[ax_idx]
    dat = panel_data[program]
    bg = dat["background_mask"]
    fg = dat["highlight_mask"]
    color = prog_color[program]
    label = label_map[program]

    # Background first; rasterized in SVG/PDF to keep file size small.
    ax.scatter(
        x[bg],
        y[bg],
        s=BACKGROUND_SIZE,
        c=BACKGROUND_COLOR,
        alpha=BACKGROUND_ALPHA,
        linewidths=0,
        edgecolors="none",
        rasterized=True,
        zorder=1,
    )

    # Foreground selected-program cells as vector points.
    ax.scatter(
        x[fg],
        y[fg],
        s=HIGHLIGHT_SIZE,
        c=color,
        alpha=HIGHLIGHT_ALPHA,
        linewidths=HIGHLIGHT_EDGEWIDTH,
        edgecolors="none",
        rasterized=False,
        zorder=2,
    )

    ax.set_title(f"{label}\n(n={dat['n_highlight']:,})", color=color, fontsize=10, fontweight="bold", pad=4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    ax.set_xlim(x.min() - x_pad, x.max() + x_pad)
    ax.set_ylim(y.min() - y_pad, y.max() + y_pad)
    for spine in ax.spines.values():
        spine.set_visible(False)

# Hide any extra axes if panel count does not fill the grid.
for extra_ax in axes[n_panels:]:
    extra_ax.axis("off")

fig.suptitle(
    f"Selected cNMF programs on integrated UMAP\nBackground color = {BACKGROUND_COLOR}; mode = {HIGHLIGHT_MODE}; threshold = {THRESHOLD}",
    fontsize=11,
    y=1.01,
)

for ext in ["svg", "pdf", "png"]:
    fig.savefig(OUTPUT_STEM.with_suffix(f".{ext}"), dpi=300 if ext == "png" else None)

plt.show()
plt.close(fig)

print("Saved outputs:")
for ext in ["svg", "pdf", "png", "summary.csv"]:
    if ext == "summary.csv":
        print(OUTPUT_STEM.with_suffix(".summary.csv"))
    else:
        print(OUTPUT_STEM.with_suffix(f".{ext}"))
