# %% [markdown]
# # All-cells UMAP colored by dominant selected cNMF program
#
# Improvements relative to the previous program-order UMAP:
#
# 1. high-contrast day-constrained palette
# 2. lighter, less conspicuous background
# 3. one colored marker per selected cell; no halo layer
# 4. foreground cells plotted as one collection, sorted dark-to-light
# 5. smaller opaque foreground points to reduce apparent stacking
# 6. direct P## labels on the UMAP at each program's median location
# 7. adaptive legend marker outlines
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
# Keep the background very faint so P22/P29 and P8/P11 remain easier to find.
# Background cells are all cells not colored as selected programs. Use a neutral,
# low-alpha gray so they provide UMAP context without reading as another program.
BACKGROUND_COLOR = "#666A70"
BACKGROUND_ALPHA = 0.55
BACKGROUND_SIZE = 0.060

# Foreground cells. This version uses one colored marker per selected cell.
# Slightly smaller opaque points reduce overplotting in dense regions.
FOREGROUND_SIZE = 0.50
TOP_OVERLAY_PROGRAMS = [22]
TOP_OVERLAY_SIZE_MULTIPLIER = 1.25
FOREGROUND_ALPHA = 1.00
RASTERIZE_FOREGROUND = False

# Plot all selected cells as one scatter collection rather than one collection
# per program. This avoids program-by-program layer artifacts in SVG/PDF.
PLOT_FOREGROUND_AS_SINGLE_COLLECTION = True

# The difficult programs are visually buried by nearby high-density programs.
# Draw all non-target programs first, then P22/P29 and P8/P11 last. This improves
# pair visibility without assigning any cell more than one program color.
DRAW_TARGET_PROGRAMS_LAST = True
TARGET_PROGRAMS_TO_DRAW_LAST = {22, 29, 8, 11}

# Halo under colored cells.
# Dark foreground colors keep a white halo. Light foreground colors get a
# low-alpha dark halo so P7/P48-like colors do not disappear on a white/light UMAP.
USE_WHITE_HALO = False
HALO_COLOR = "white"
HALO_ALPHA = 0.76
HALO_SIZE_MULTIPLIER = 2.6
ADAPTIVE_HALO_FOR_LIGHT_COLORS = True
LIGHT_COLOR_LUMINANCE_CUTOFF = 0.45
VERY_LIGHT_COLOR_LUMINANCE_CUTOFF = 0.75
DARK_HALO_COLOR = "#111111"
DARK_HALO_ALPHA = 0.30
DARK_HALO_SIZE_MULTIPLIER = 2.65
VERY_LIGHT_DARK_HALO_ALPHA = 0.42
VERY_LIGHT_DARK_HALO_SIZE_MULTIPLIER = 3.2
RASTERIZE_HALO = RASTERIZE_FOREGROUND

# Program-specific halo overrides for light colors that otherwise blend into
# the white/light UMAP background or into neighboring programs.
# tuple = (halo_color, halo_alpha, halo_size_multiplier)
PROGRAM_HALO_OVERRIDES = {
    22: ("#303030", 0.36, 2.85),  # very light D0 P22 if halo mode is enabled
    48: ("#174A33", 0.32, 2.75),  # pale green P48 on light background
    11: ("#5A2D00", 0.22, 2.35),  # darker orange P11
    34: ("#004E64", 0.22, 2.35),  # cyan P34
}

# Legend marker edges. Light markers need a dark edge on the white legend area;
# dark markers keep a white edge.
LEGEND_MID_LIGHT_EDGE_LUMINANCE_CUTOFF = 0.35
LEGEND_LIGHT_EDGE_LUMINANCE_CUTOFF = 0.45
LEGEND_VERY_LIGHT_EDGE_LUMINANCE_CUTOFF = 0.75
LEGEND_LIGHT_EDGE_COLOR = "#222222"
LEGEND_DARK_EDGE_COLOR = "white"
LEGEND_MID_LIGHT_EDGE_WIDTH = 0.85
LEGEND_LIGHT_EDGE_WIDTH = 1.05
LEGEND_VERY_LIGHT_EDGE_WIDTH = 1.20
LEGEND_DARK_EDGE_WIDTH = 0.65

# Direct labels near each program region.
# Keep labels off the UMAP so text does not obscure points.
ANNOTATE_PROGRAMS = False
LABEL_PROGRAMS = None
LABEL_FONT_SIZE = 6.2
LABEL_ALPHA = 0.96
LABEL_BOX_ALPHA = 0.84

# Optional manual label offsets in UMAP units. These prevent the pair labels
# from sitting directly on top of each other when two programs overlap.
LABEL_OFFSETS = {}

# Extra pair-only outputs. These do not change the main UMAP; they create
# two focused helper panels where all other cells are greyed out.
MAKE_PAIR_SPOTLIGHTS = True
PAIR_SPOTLIGHT_PROGRAMS = [(22, 29), (7, 29)]
SPOTLIGHT_BACKGROUND_COLOR = "#666A70"
SPOTLIGHT_OTHER_SELECTED_COLOR = "#D8DDE3"
SPOTLIGHT_OTHER_SELECTED_ALPHA = 0.14
SPOTLIGHT_PAIR_SIZE = 0.82
SPOTLIGHT_PAIR_ALPHA = 1.00
SPOTLIGHT_LABEL_FONT_SIZE = 7.0
ANNOTATE_SPOTLIGHTS = False

FIGSIZE = (5.2, 5.2)
OUTPUT_STEM = OUT_DIR / "umap_program_order_pair_visibility_optimized_v11_p22_blue_darkest_bg"

# %%
# Distinct palette constrained by the requested day/color-family rules.
#
# I treated the first gray/blue palette as Day 0 because Day 1 is specified
# again below with the green palette. Pure white was not assigned because the
# UMAP background is white/light.
HIGH_CONTRAST_PROGRAM_COLORS = {
    # Day 0 programs: keep P29/P7 in the warm brown family, but move P22 to a
    # saturated royal blue because the tan P22 was getting buried in the D0 island.
    22: "#0057ff",
    29: "#6f1d1b",
    7:  "#99582a",

    # Day 1 programs: keep the previously working green palette.
    28: "#081c15",
    30: "#74c69d",
    13: "#2d6a4f",
    48: "#d8f3dc",

    # Day 2 programs: P14 remains the perceptual midpoint between P15 and P17.
    15: "#bb9457",
    14: "#ae8698",
    17: "#907ad6",

    # Day 3 P2/P3/P8/P11: keep P8 bright yellow and P11 darker orange.
    # This increases separation between the nearby yellow/orange arterial programs.
    2:  "#372b2b",
    3:  "#850a0a",
    8:  "#f6de26",
    11: "#d97706",

    # Day 3 P33/P34: blue palette, dark + cyan.
    33: "#03045e",
    34: "#48cae4",
}

# Optional family labels for summary output.
PROGRAM_COLOR_FAMILY = {
    22: "D0_high_contrast", 29: "D0_high_contrast", 7: "D0_high_contrast",
    28: "D1_green", 30: "D1_green", 13: "D1_green", 48: "D1_green",
    15: "D2_brown_purple", 14: "D2_brown_purple", 17: "D2_brown_purple",
    2: "D3_red_orange", 3: "D3_red_orange", 8: "D3_red_orange", 11: "D3_red_orange",
    33: "D3_blue", 34: "D3_blue",
}


def relative_luminance(color: str) -> float:
    """Return WCAG-style relative luminance for a Matplotlib color."""
    rgb = np.array(mpl.colors.to_rgb(color), dtype=float)
    linear_rgb = np.where(
        rgb <= 0.04045,
        rgb / 12.92,
        ((rgb + 0.055) / 1.055) ** 2.4,
    )
    return float(np.dot(linear_rgb, [0.2126, 0.7152, 0.0722]))


def halo_style_for_program(program: int, color: str) -> tuple:
    """Return halo color, alpha, and size multiplier for one program."""
    if program in PROGRAM_HALO_OVERRIDES:
        return PROGRAM_HALO_OVERRIDES[program]

    luminance = relative_luminance(color)
    if ADAPTIVE_HALO_FOR_LIGHT_COLORS and luminance >= VERY_LIGHT_COLOR_LUMINANCE_CUTOFF:
        return DARK_HALO_COLOR, VERY_LIGHT_DARK_HALO_ALPHA, VERY_LIGHT_DARK_HALO_SIZE_MULTIPLIER
    if ADAPTIVE_HALO_FOR_LIGHT_COLORS and luminance >= LIGHT_COLOR_LUMINANCE_CUTOFF:
        return DARK_HALO_COLOR, DARK_HALO_ALPHA, DARK_HALO_SIZE_MULTIPLIER
    return HALO_COLOR, HALO_ALPHA, HALO_SIZE_MULTIPLIER


def legend_style_for_color(color: str) -> tuple:
    """Return legend marker edge color and edge width for one fill color."""
    luminance = relative_luminance(color)
    if luminance >= LEGEND_VERY_LIGHT_EDGE_LUMINANCE_CUTOFF:
        return LEGEND_LIGHT_EDGE_COLOR, LEGEND_VERY_LIGHT_EDGE_WIDTH
    if luminance >= LEGEND_LIGHT_EDGE_LUMINANCE_CUTOFF:
        return LEGEND_LIGHT_EDGE_COLOR, LEGEND_LIGHT_EDGE_WIDTH
    if luminance >= LEGEND_MID_LIGHT_EDGE_LUMINANCE_CUTOFF:
        return LEGEND_LIGHT_EDGE_COLOR, LEGEND_MID_LIGHT_EDGE_WIDTH
    return LEGEND_DARK_EDGE_COLOR, LEGEND_DARK_EDGE_WIDTH

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
        "relative_luminance": [relative_luminance(prog_color[p]) for p in unique_programs],
        "halo_color": [halo_style_for_program(p, prog_color[p])[0] for p in unique_programs],
        "halo_alpha": [halo_style_for_program(p, prog_color[p])[1] for p in unique_programs],
        "halo_size_multiplier": [halo_style_for_program(p, prog_color[p])[2] for p in unique_programs],
        "legend_edge_color": [legend_style_for_color(prog_color[p])[0] for p in unique_programs],
        "legend_edge_width": [legend_style_for_color(prog_color[p])[1] for p in unique_programs],
        "drawn_last_target": [p in TARGET_PROGRAMS_TO_DRAW_LAST for p in unique_programs],
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
            "relative_luminance": relative_luminance(prog_color.get(program, "#999999")),
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
# Draw darker non-target colors first and lighter non-target colors later.
# The hard-to-find target programs are handled specially in the foreground sort
# below so they are not buried by abundant light programs such as P7.
program_counts = {
    program: int((selected_mask & (dominant_program == program)).sum())
    for program in unique_programs
}
program_luminance = {
    program: relative_luminance(prog_color.get(program, "#999999"))
    for program in unique_programs
}
plot_order = sorted(
    unique_programs,
    key=lambda p: (p in TARGET_PROGRAMS_TO_DRAW_LAST, program_luminance[p], -program_counts[p]),
)
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

# Foreground cells.
# Default mode here: one selected cell = one colored marker.
# The points are sorted dark-to-light so light colors are still drawn later,
# but there is no extra halo layer and no separate scatter layer per program.
if PLOT_FOREGROUND_AS_SINGLE_COLLECTION:
    fg_idx = np.flatnonzero(selected_mask)
    if len(fg_idx) > 0:
        fg_luminance = np.array(
            [program_luminance.get(int(dominant_program[i]), 0.5) for i in fg_idx],
            dtype=float,
        )

        if DRAW_TARGET_PROGRAMS_LAST:
            target_flags = np.array(
                [int(dominant_program[i]) in TARGET_PROGRAMS_TO_DRAW_LAST for i in fg_idx],
                dtype=bool,
            )
            non_target_idx = fg_idx[~target_flags]
            non_target_luminance = fg_luminance[~target_flags]
            target_idx = fg_idx[target_flags]

            # Non-target cells keep the old dark-to-light ordering. Target cells
            # are appended in original observation order, which avoids forcing
            # one member of each target pair to fully cover the other.
            non_target_idx = non_target_idx[np.argsort(non_target_luminance, kind="stable")]
            fg_idx = np.concatenate([non_target_idx, target_idx])
        else:
            fg_idx = fg_idx[np.argsort(fg_luminance, kind="stable")]

        fg_colors = [prog_color.get(int(dominant_program[i]), "#999999") for i in fg_idx]

        ax.scatter(
            x[fg_idx],
            y[fg_idx],
            c=fg_colors,
            s=FOREGROUND_SIZE,
            alpha=FOREGROUND_ALPHA,
            rasterized=RASTERIZE_FOREGROUND,
            linewidths=0,
            edgecolors="none",
            zorder=3,
        )

        # Re-draw the most buried programs after everything else so they sit fully
        # in front. This does not change the one-cell/one-program assignment; it
        # only affects visual stacking order.
        for top_program in TOP_OVERLAY_PROGRAMS:
            top_mask = selected_mask & (dominant_program == top_program)
            if top_mask.sum() == 0:
                continue
            ax.scatter(
                x[top_mask],
                y[top_mask],
                c=prog_color.get(top_program, "#999999"),
                s=FOREGROUND_SIZE * TOP_OVERLAY_SIZE_MULTIPLIER,
                alpha=FOREGROUND_ALPHA,
                rasterized=RASTERIZE_FOREGROUND,
                linewidths=0,
                edgecolors="none",
                zorder=4,
            )
else:
    # Optional legacy mode: one scatter collection per program, with optional halo.
    # This can make overplotting appear more layered, but can be useful if you
    # want program groups to be easier to select separately in Illustrator.
    for program in plot_order:
        mask = selected_mask & (dominant_program == program)
        if mask.sum() == 0:
            continue

        this_color = prog_color.get(program, "#999999")
        halo_color, halo_alpha, halo_size_multiplier = halo_style_for_program(program, this_color)

        if USE_WHITE_HALO:
            ax.scatter(
                x[mask],
                y[mask],
                c=halo_color,
                s=FOREGROUND_SIZE * halo_size_multiplier,
                alpha=halo_alpha,
                rasterized=RASTERIZE_HALO,
                linewidths=0,
                edgecolors="none",
                zorder=2,
            )

        ax.scatter(
            x[mask],
            y[mask],
            c=this_color,
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
        if LABEL_PROGRAMS is not None and program not in LABEL_PROGRAMS:
            continue
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
    legend_color = prog_color.get(program, "#999999")
    legend_edge_color, legend_edge_width = legend_style_for_color(legend_color)
    handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=legend_color,
            markeredgecolor=legend_edge_color,
            markeredgewidth=legend_edge_width,
            markersize=4.4,
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


# %%
def save_pair_spotlight(pair_programs: tuple[int, int]) -> None:
    """Save a focused UMAP for one hard-to-see program pair."""
    pair_programs = tuple(int(p) for p in pair_programs)
    pair_label = "_".join(f"P{p}" for p in pair_programs)

    fig_pair, ax_pair = plt.subplots(figsize=FIGSIZE)

    # Very faint all-cell background.
    ax_pair.scatter(
        x,
        y,
        c=SPOTLIGHT_BACKGROUND_COLOR,
        s=BACKGROUND_SIZE,
        alpha=BACKGROUND_ALPHA,
        rasterized=True,
        linewidths=0,
        edgecolors="none",
        zorder=1,
    )

    # Grey out selected cells from other programs so the pair locations are visible.
    other_selected_mask = selected_mask & ~np.isin(dominant_program, pair_programs)
    ax_pair.scatter(
        x[other_selected_mask],
        y[other_selected_mask],
        c=SPOTLIGHT_OTHER_SELECTED_COLOR,
        s=FOREGROUND_SIZE,
        alpha=SPOTLIGHT_OTHER_SELECTED_ALPHA,
        rasterized=True,
        linewidths=0,
        edgecolors="none",
        zorder=2,
    )

    # Plot the pair. Use a thin adaptive edge only in the spotlight panels so
    # very light colors remain visible; the main UMAP remains one face color per cell.
    for program in pair_programs:
        mask = selected_mask & (dominant_program == program)
        if mask.sum() == 0:
            continue
        this_color = prog_color.get(program, "#999999")
        edge_color, edge_width = legend_style_for_color(this_color)
        ax_pair.scatter(
            x[mask],
            y[mask],
            c=this_color,
            s=SPOTLIGHT_PAIR_SIZE,
            alpha=SPOTLIGHT_PAIR_ALPHA,
            rasterized=RASTERIZE_FOREGROUND,
            linewidths=edge_width * 0.20,
            edgecolors=edge_color,
            zorder=4,
        )

        if ANNOTATE_SPOTLIGHTS:
            x_med = float(np.median(x[mask]))
            y_med = float(np.median(y[mask]))
            dx, dy = LABEL_OFFSETS.get(program, (0.0, 0.0))
            ax_pair.text(
                x_med + dx,
                y_med + dy,
                f"P{program}",
                color=this_color,
                fontsize=SPOTLIGHT_LABEL_FONT_SIZE,
                fontweight="bold",
                ha="center",
                va="center",
                alpha=LABEL_ALPHA,
                bbox=dict(
                    boxstyle="round,pad=0.20",
                    facecolor="white",
                    edgecolor="none",
                    alpha=LABEL_BOX_ALPHA,
                ),
                zorder=10,
            )

    ax_pair.set_title(pair_label.replace("_", " vs "), fontsize=7, pad=2)
    ax_pair.set_xticks([])
    ax_pair.set_yticks([])
    ax_pair.set_aspect("equal")
    for spine in ax_pair.spines.values():
        spine.set_visible(False)

    pair_handles = []
    for program in pair_programs:
        label = label_map.get(program, f"P{program}")
        short_label = label.split(" — ")[-1] if " — " in label else label
        legend_color = prog_color.get(program, "#999999")
        legend_edge_color, legend_edge_width = legend_style_for_color(legend_color)
        pair_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=legend_color,
                markeredgecolor=legend_edge_color,
                markeredgewidth=legend_edge_width,
                markersize=4.8,
                label=f"P{program} {short_label}",
            )
        )

    ax_pair.legend(
        handles=pair_handles,
        fontsize=5.2,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        ncol=1,
        handletextpad=0.35,
        labelspacing=0.35,
    )

    plt.tight_layout()
    pair_stem = OUTPUT_STEM.with_name(f"{OUTPUT_STEM.name}_spotlight_{pair_label}")
    pair_pdf_path = pair_stem.with_suffix(".pdf")
    pair_png_path = pair_stem.with_suffix(".png")
    pair_svg_path = pair_stem.with_suffix(".svg")
    fig_pair.savefig(pair_pdf_path, bbox_inches="tight")
    fig_pair.savefig(pair_png_path, bbox_inches="tight")
    fig_pair.savefig(pair_svg_path, format="svg", bbox_inches="tight")
    print(f"Saved pair spotlight: {pair_pdf_path}")
    print(f"Saved pair spotlight: {pair_png_path}")
    print(f"Saved pair spotlight editable SVG: {pair_svg_path}")
    plt.show()


if MAKE_PAIR_SPOTLIGHTS:
    for pair in PAIR_SPOTLIGHT_PROGRAMS:
        save_pair_spotlight(pair)
