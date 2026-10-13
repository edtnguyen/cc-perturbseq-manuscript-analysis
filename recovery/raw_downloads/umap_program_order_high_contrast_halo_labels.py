# %% [markdown]
# # All-cells UMAP colored by dominant selected cNMF program
#
# Improvements relative to the previous program-order UMAP:
#
# 1. high-contrast day-constrained palette
# 2. lighter, less conspicuous background
# 3. white halo under colored foreground cells
# 4. abundant programs plotted first, rare programs plotted last
# 5. larger foreground points and smaller background points
# 6. direct P## labels on the UMAP at each program's median location
#
# Background cells are rasterized to keep SVG size smaller. Foreground cells are
# vector by default so the colored program points remain Illustrator-editable.
# Set RASTERIZE_FOREGROUND = True if the SVG becomes too large.

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

try:
    from IPython.display import display
except Exception:
    def display(x):
        print(x.to_string() if hasattr(x, "to_string") else x)

warnings.filterwarnings("ignore")

# %%
def find_repo_root(start: Path = Path.cwd()) -> Path:
    """Find the repository root by looking upward for expected project files."""
    start = start.resolve()
    required_paths = [
        Path("data/program_activity_long_v4.csv"),
        Path("analysis/regulator_dotplots/program_order.txt"),
    ]

    for candidate in [start, *start.parents]:
        if all((candidate / rel_path).exists() for rel_path in required_paths):
            return candidate

    raise FileNotFoundError(
        "Could not auto-detect REPO_ROOT. Set REPO_ROOT manually in this cell."
    )


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
        "font.size": 7,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",  # keep text editable in Illustrator
    }
)

print(f"REPO_ROOT = {REPO_ROOT}")
print(f"OUT_DIR   = {OUT_DIR}")

# %%
# -----------------------------
# Input paths
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

PROGRAM_META_PATH = REPO_ROOT / "data" / "program_activity_long_v4.csv"
PROGRAM_ORDER_PATH = REPO_ROOT / "analysis" / "regulator_dotplots" / "program_order.txt"

# %%
# -----------------------------
# Plot settings
# -----------------------------
THRESHOLD = 0.02

# "global": color a cell only if its dominant program across ALL cNMF programs
#           is listed in program_order.txt.
# "program_order": old behavior; color by highest-scoring program among only
#                  the programs listed in program_order.txt.
DOMINANCE_MODE = "global"

# Lighter and less conspicuous than the previous gray/lavender background.
# This also avoids competing with D0 gray programs.
BACKGROUND_COLOR = "#F2F5FA"
BACKGROUND_ALPHA = 0.42
BACKGROUND_SIZE = 0.10

# Larger foreground points. Foreground is vector by default for Illustrator editing.
FOREGROUND_SIZE = 0.85
FOREGROUND_ALPHA = 0.95
RASTERIZE_FOREGROUND = False

# White halo under colored cells to separate similar within-day shades.
USE_WHITE_HALO = True
HALO_COLOR = "white"
HALO_ALPHA = 0.80
HALO_SIZE_MULTIPLIER = 2.6
RASTERIZE_HALO = RASTERIZE_FOREGROUND

# Direct labels near each program region.
ANNOTATE_PROGRAMS = True
LABEL_FONT_SIZE = 6.0
LABEL_ALPHA = 0.92
LABEL_BOX_ALPHA = 0.78

# Optional manual label offsets in UMAP units, if two labels overlap.
# Example: LABEL_OFFSETS = {23: (0.25, -0.1), 15: (-0.15, 0.2)}
LABEL_OFFSETS = {}

FIGSIZE = (5.2, 5.2)
OUTPUT_STEM = OUT_DIR / "umap_program_order_high_contrast_halo_labels"

# %%
# High-contrast palette constrained by day/color-family rules.
# D0: gray; D1: orange; D2: darker orange/brown;
# D3 arterial EC: red; D3 non-EC: blue.
HIGH_CONTRAST_PROGRAM_COLORS = {
    # Day 0 programs: grays
    22: "#1F1F1F",
    29: "#5C5C5C",
    7:  "#9E9E9E",
    49: "#C7C7C7",

    # Day 1 programs: oranges
    30: "#FFD166",
    28: "#FFA600",
    13: "#F97316",
    48: "#C2410C",

    # Day 2 programs: darker oranges/browns
    15: "#B45309",
    23: "#92400E",
    19: "#7C2D12",
    17: "#5F2608",
    14: "#3F1D05",

    # Day 3 arterial EC programs: reds
    11: "#FB6A4A",
    2:  "#A50F15",

    # Day 3 non-EC programs: blues
    34: "#B3E5FC",
    27: "#56B4E9",
    24: "#1F78B4",
    16: "#005AB5",
    44: "#003F7F",
    33: "#001D4A",
}

# Optional family labels for summary output.
PROGRAM_COLOR_FAMILY = {
    22: "D0_gray", 29: "D0_gray", 7: "D0_gray", 49: "D0_gray",
    30: "D1_orange", 28: "D1_orange", 13: "D1_orange", 48: "D1_orange",
    15: "D2_dark_orange", 23: "D2_dark_orange", 19: "D2_dark_orange", 17: "D2_dark_orange", 14: "D2_dark_orange",
    11: "D3_arterial_red", 2: "D3_arterial_red",
    34: "D3_nonEC_blue", 27: "D3_nonEC_blue", 24: "D3_nonEC_blue", 16: "D3_nonEC_blue", 44: "D3_nonEC_blue", 33: "D3_nonEC_blue",
}

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

    print(f"Transferred {obsm_key} from adata_int: {aligned_umap.shape[0]} cells x {aligned_umap.shape[1]} dimensions")
    return aligned_umap


def parse_program_id(x):
    """Parse cNMF program IDs from var_names like '13' or 'P13'."""
    s = str(x)
    if s.startswith("P"):
        s = s[1:]
    return int(s)

# %%
for path in [H5MU_PATH, ADATA_INT_PATH, PROGRAM_META_PATH, PROGRAM_ORDER_PATH]:
    if not path.exists():
        raise FileNotFoundError(path)

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
# Load metadata labels and default colors, then override selected program colors
# with the high-contrast day-constrained palette above.
prog_meta = pd.read_csv(PROGRAM_META_PATH).drop_duplicates("program_name")
prog_meta["program_name"] = prog_meta["program_name"].astype(int)
program_ids = prog_meta["program_name"].astype(int)

if "program_label" in prog_meta.columns:
    label_map = dict(zip(program_ids, prog_meta["program_label"].astype(str)))
else:
    label_map = {int(program): f"P{int(program)}" for program in program_ids}

preferred_color_columns = ["distinct_color", "program_color", "color", "class_color"]
color_column = next((col for col in preferred_color_columns if col in prog_meta.columns), None)

prog_color = {}
if color_column is not None:
    print(f"Using `{color_column}` from program metadata as baseline colors")
    for program, color in zip(program_ids, prog_meta[color_column]):
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

# Override with the high-contrast colors for the selected program-order programs.
prog_color.update(HIGH_CONTRAST_PROGRAM_COLORS)

# Parse program_order.txt.
programs = []
for line in PROGRAM_ORDER_PATH.read_text().splitlines():
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    programs.append(int(line))

unique_programs = list(dict.fromkeys(programs))
print(f"Loaded {len(unique_programs)} unique programs from program_order.txt:")
print(unique_programs)

# Map cNMF score columns to zero-based program IDs.
try:
    cnmf_program_ids = [parse_program_id(v) for v in cnmf.var_names]
except Exception:
    cnmf_program_ids = list(range(scores.shape[1]))

if len(cnmf_program_ids) != scores.shape[1]:
    raise ValueError("Could not map cNMF score columns to program IDs.")

program_to_col = {pid: idx for idx, pid in enumerate(cnmf_program_ids)}
missing_from_scores = [p for p in unique_programs if p not in program_to_col]
if missing_from_scores:
    raise ValueError(f"These program_order programs are missing from cNMF score matrix: {missing_from_scores}")

# Fallback colors for any program in program_order not explicitly assigned.
fallback_cmap = plt.get_cmap("tab20", max(len(unique_programs), 1))
for idx, program in enumerate(unique_programs):
    if program not in prog_color:
        prog_color[program] = mpl.colors.to_hex(fallback_cmap(idx))

print(f"Color entries available for {len(prog_color)} programs")

# Save the exact color table used in this plot.
color_table = pd.DataFrame(
    {
        "program": unique_programs,
        "program_label": [label_map.get(p, f"P{p}") for p in unique_programs],
        "hex_color": [prog_color[p] for p in unique_programs],
        "color_family": [PROGRAM_COLOR_FAMILY.get(p, "fallback_or_metadata") for p in unique_programs],
    }
)
color_table.to_csv(OUTPUT_STEM.with_suffix(".colors_used.csv"), index=False)
display(color_table)

# %%
# Assign dominant program per cell.
if DOMINANCE_MODE == "global":
    dominant_col_idx = np.argmax(scores, axis=1)
    max_score = scores[np.arange(scores.shape[0]), dominant_col_idx]
    dominant_program = np.array([cnmf_program_ids[idx] for idx in dominant_col_idx], dtype=int)
    selected_mask = np.isin(dominant_program, unique_programs) & (max_score > THRESHOLD)
elif DOMINANCE_MODE == "program_order":
    selected_cols = [program_to_col[p] for p in unique_programs]
    sub_scores = scores[:, selected_cols]
    max_prog_idx = np.argmax(sub_scores, axis=1)
    max_score = np.max(sub_scores, axis=1)
    dominant_program = np.array([unique_programs[idx] for idx in max_prog_idx], dtype=int)
    selected_mask = max_score > THRESHOLD
else:
    raise ValueError("DOMINANCE_MODE must be 'global' or 'program_order'")

summary_records = []
for program in unique_programs:
    mask = selected_mask & (dominant_program == program)
    summary_records.append(
        {
            "program": program,
            "program_label": label_map.get(program, f"P{program}"),
            "hex_color": prog_color.get(program, "#999999"),
            "color_family": PROGRAM_COLOR_FAMILY.get(program, "fallback_or_metadata"),
            "n_cells": int(mask.sum()),
        }
    )
summary_df = pd.DataFrame(summary_records)
summary_df.to_csv(OUTPUT_STEM.with_suffix(".summary.csv"), index=False)
display(summary_df)

print(f"DOMINANCE_MODE = {DOMINANCE_MODE!r}")
print(f"Assigned colors to {int(selected_mask.sum()):,} cells above threshold {THRESHOLD}")
print(f"Background cells: {int((~selected_mask).sum()):,}")

# %%
# Plot abundant programs first and rare programs last so rare programs stay visible.
program_counts = {
    program: int((selected_mask & (dominant_program == program)).sum())
    for program in unique_programs
}
plot_order = sorted(unique_programs, key=lambda p: program_counts[p], reverse=True)
legend_order = unique_programs  # keep legend in program_order.txt order

fig, ax = plt.subplots(figsize=FIGSIZE)

x = umap_coords[:, 0]
y = umap_coords[:, 1]

# Background cells: rasterized, small, light.
bg_mask = ~selected_mask
ax.scatter(
    x[bg_mask],
    y[bg_mask],
    c=BACKGROUND_COLOR,
    s=BACKGROUND_SIZE,
    alpha=BACKGROUND_ALPHA,
    rasterized=True,
    linewidths=0,
    edgecolors="none",
    zorder=1,
)

# Foreground cells: halo + color.
for program in plot_order:
    mask = selected_mask & (dominant_program == program)
    if mask.sum() == 0:
        continue

    if USE_WHITE_HALO:
        ax.scatter(
            x[mask],
            y[mask],
            c=HALO_COLOR,
            s=FOREGROUND_SIZE * HALO_SIZE_MULTIPLIER,
            alpha=HALO_ALPHA,
            rasterized=RASTERIZE_HALO,
            linewidths=0,
            edgecolors="none",
            zorder=2,
        )

    ax.scatter(
        x[mask],
        y[mask],
        c=prog_color.get(program, "#999999"),
        s=FOREGROUND_SIZE,
        alpha=FOREGROUND_ALPHA,
        rasterized=RASTERIZE_FOREGROUND,
        linewidths=0,
        edgecolors="none",
        zorder=3,
    )

# Direct P labels near robust median position of each program.
if ANNOTATE_PROGRAMS:
    for program in legend_order:
        mask = selected_mask & (dominant_program == program)
        if mask.sum() == 0:
            continue
        x_med = float(np.median(x[mask]))
        y_med = float(np.median(y[mask]))
        dx, dy = LABEL_OFFSETS.get(program, (0.0, 0.0))
        ax.text(
            x_med + dx,
            y_med + dy,
            f"P{program}",
            color=prog_color.get(program, "#999999"),
            fontsize=LABEL_FONT_SIZE,
            fontweight="bold",
            ha="center",
            va="center",
            alpha=LABEL_ALPHA,
            bbox=dict(
                boxstyle="round,pad=0.18",
                facecolor="white",
                edgecolor="none",
                alpha=LABEL_BOX_ALPHA,
            ),
            zorder=10,
        )

ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect("equal")
for spine in ax.spines.values():
    spine.set_visible(False)

handles = []
for program in legend_order:
    label = label_map.get(program, f"P{program}")
    short_label = label.split(" — ")[-1] if " — " in label else label
    handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=prog_color.get(program, "#999999"),
            markeredgecolor="white",
            markeredgewidth=0.8,
            markersize=4.2,
            label=f"P{program} {short_label}",
        )
    )

ax.legend(
    handles=handles,
    fontsize=4.2,
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.0, 0.5),
    ncol=1,
    handletextpad=0.35,
    labelspacing=0.35,
)

plt.tight_layout()

pdf_path = OUTPUT_STEM.with_suffix(".pdf")
png_path = OUTPUT_STEM.with_suffix(".png")
svg_path = OUTPUT_STEM.with_suffix(".svg")

fig.savefig(pdf_path, bbox_inches="tight")
fig.savefig(png_path, bbox_inches="tight")
fig.savefig(svg_path, format="svg", bbox_inches="tight")

print(f"Saved: {pdf_path}")
print(f"Saved: {png_path}")
print(f"Saved editable SVG: {svg_path}")
print(f"Saved color table: {OUTPUT_STEM.with_suffix('.colors_used.csv')}")
print(f"Saved summary table: {OUTPUT_STEM.with_suffix('.summary.csv')}")
plt.show()
