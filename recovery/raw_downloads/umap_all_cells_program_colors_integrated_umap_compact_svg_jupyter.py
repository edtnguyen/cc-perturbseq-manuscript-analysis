# %% [markdown]
# # All-cells UMAP colored by dominant gene program
#
# Each cell is colored by its highest-scoring program from the `program_order.txt` list.
#
# This version loads cNMF scores from `mdata`, but uses `X_umap` from the integrated AnnData object (`adata_int`). The notebook checks that `cnmf.obs_names` and `adata_int.obs_names` match, then reorders the integrated UMAP coordinates to the cNMF cell order before plotting.
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


# In a notebook, __file__ is usually not defined, so use an auto-detected repo root.
# If auto-detection fails, replace this with: REPO_ROOT = Path("/path/to/your/repo")
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
        # Keep text editable in Illustrator when saving SVG.
        "svg.fonttype": "none",
    }
)

print(f"REPO_ROOT = {REPO_ROOT}")
print(f"OUT_DIR   = {OUT_DIR}")

# %%
def align_umap_to_cnmf_cells(
    *,
    cnmf,
    adata_int,
    obsm_key: str = "X_umap",
    require_exact_match: bool = True,
) -> np.ndarray:
    """Return adata_int.obsm[obsm_key] reordered to cnmf.obs_names.

    The plotting code uses cNMF scores, so the UMAP coordinate rows must be in
    the same order as cnmf.obs_names. This function verifies that mapping.
    """
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
H5MU_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/adata/cNMF_50_0_2.h5mu"
)

ADATA_INT_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
    "2.preprocessing_tonys_pipeline/3.integrated_days/"
    "3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad"
)

if not H5MU_PATH.exists():
    raise FileNotFoundError(f"Could not find h5mu file: {H5MU_PATH}")

if not ADATA_INT_PATH.exists():
    raise FileNotFoundError(f"Could not find integrated h5ad file: {ADATA_INT_PATH}")

print("Loading h5mu...")
mdata = mu.read_h5mu(H5MU_PATH)
cnmf = mdata.mod["cNMF"]
print(f"Loaded cNMF modality: {cnmf.n_obs} cells x {cnmf.n_vars} programs")

print("Loading integrated AnnData for UMAP...")
# backed="r" avoids loading the full expression matrix; obs/obsm are enough here.
adata_int = sc.read_h5ad(ADATA_INT_PATH, backed="r")
print(f"Loaded integrated AnnData: {adata_int.n_obs} cells x {adata_int.n_vars} genes")

umap_coords = align_umap_to_cnmf_cells(
    cnmf=cnmf,
    adata_int=adata_int,
    obsm_key="X_umap",
    require_exact_match=True,
)

# Store the transferred coordinates in memory too, so downstream cells can inspect them.
cnmf.obsm["X_umap_from_adata_int"] = umap_coords

scores = cnmf.X.toarray() if hasattr(cnmf.X, "toarray") else np.asarray(cnmf.X)
print(f"scores shape = {scores.shape}")
print(f"umap_coords shape = {umap_coords.shape}")

# %%
prog_meta = pd.read_csv(REPO_ROOT / "data" / "program_activity_long_v4.csv").drop_duplicates(
    "program_name"
)
program_ids = prog_meta["program_name"].astype(int)

if "program_label" in prog_meta.columns:
    label_map = dict(zip(program_ids, prog_meta["program_label"]))
else:
    label_map = {int(program): f"P{int(program)}" for program in program_ids}

preferred_color_columns = [
    "distinct_color",
    "program_color",
    "color",
    "class_color",
]
color_column = next((col for col in preferred_color_columns if col in prog_meta.columns), None)

prog_color = {}
if color_column is not None:
    print(f"Using `{color_column}` from program metadata for colors")
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
else:
    print(
        "No color column found in program metadata; generating fallback colors. "
        f"Available columns: {list(prog_meta.columns)}"
    )

program_order_path = REPO_ROOT / "analysis" / "regulator_dotplots" / "program_order.txt"
programs = []
for line in program_order_path.read_text().splitlines():
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    programs.append(int(line))

unique_programs = list(dict.fromkeys(programs))

fallback_cmap = plt.get_cmap("tab20", max(len(unique_programs), 1))
for idx, program in enumerate(unique_programs):
    if program not in prog_color:
        prog_color[program] = mpl.colors.to_hex(fallback_cmap(idx))

print(f"Loaded {cnmf.n_obs} cells, {len(unique_programs)} programs")
print(f"Color entries available for {len(prog_color)} programs")

# %%
# Assign dominant program per cell.
threshold = 0.02

sub_scores = scores[:, unique_programs]
max_prog_idx = np.argmax(sub_scores, axis=1)
max_score = np.max(sub_scores, axis=1)
dominant_program = np.array([unique_programs[idx] for idx in max_prog_idx])

n_cells = len(max_prog_idx)
cell_colors = np.array(["#e0e0e0"] * n_cells, dtype=object)

for program in unique_programs:
    mask = (max_score > threshold) & (dominant_program == program)
    cell_colors[mask] = prog_color.get(program, "#999999")

print(f"Assigned colors to {(max_score > threshold).sum()} cells above threshold {threshold}")

# %%
fig, ax = plt.subplots(figsize=(5, 5))

# Keep the dense UMAP point cloud rasterized in SVG as embedded image layers.
# This keeps the SVG much smaller while preserving editable text/legend elements.

grey_mask = cell_colors == "#e0e0e0"
ax.scatter(
    umap_coords[grey_mask, 0],
    umap_coords[grey_mask, 1],
    c="#e8e8e8",
    s=0.1,
    rasterized=True,
    linewidths=0,
    zorder=1,
)
for program in unique_programs:
    color = prog_color.get(program, "#999999")
    mask = (max_score > threshold) & (dominant_program == program)

    if mask.sum() > 0:
        ax.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            c=color,
            s=0.3,
            rasterized=True,
            linewidths=0,
            zorder=2,
            alpha=0.7,
        )
ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect("equal")
for spine in ax.spines.values():
    spine.set_visible(False)

handles = []
for program in unique_programs:
    label = label_map.get(program, f"P{program}")
    short_label = label.split(" — ")[-1] if " — " in label else label
    handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=prog_color.get(program, "#999999"),
            markersize=3,
            label=f"P{program} {short_label}",
        )
    )

ax.legend(
    handles=handles,
    fontsize=3.5,
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.0, 0.5),
    ncol=1,
    handletextpad=0.3,
    labelspacing=0.3,
)

plt.tight_layout()

pdf_path = OUT_DIR / "umap_all_cells_program_colors_integrated_umap.pdf"
png_path = OUT_DIR / "umap_all_cells_program_colors_integrated_umap.png"
svg_path = OUT_DIR / "umap_all_cells_program_colors_integrated_umap_compact.svg"

# PDF/PNG/SVG: keep the dense point cloud rasterized so files stay manageable.
# Because `svg.fonttype` is set to `none`, text and legend labels stay editable in Illustrator.
fig.savefig(pdf_path, bbox_inches="tight")
fig.savefig(png_path, bbox_inches="tight")
fig.savefig(svg_path, format="svg", bbox_inches="tight")

print(f"Saved: {pdf_path}")
print(f"Saved: {png_path}")
print(f"Saved compact SVG: {svg_path}")
plt.show()

