# %% [markdown]
# # All-cells UMAP colored by dominant gene program
#
# Each cell is colored by its highest-scoring program from the `program_order.txt` list. Colors are hand-picked to be distinguishable across all programs, with related programs using similar hue families.

# %%
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
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
    }
)

print(f"REPO_ROOT = {REPO_ROOT}")
print(f"OUT_DIR   = {OUT_DIR}")

# %%
H5MU_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/adata/cNMF_50_0_2.h5mu"
)

if not H5MU_PATH.exists():
    raise FileNotFoundError(f"Could not find h5mu file: {H5MU_PATH}")

print("Loading h5mu...")
mdata = mu.read_h5mu(H5MU_PATH)
cnmf = mdata.mod["cNMF"]

umap_coords = np.asarray(cnmf.obsm["X_umap"])
scores = cnmf.X.toarray() if hasattr(cnmf.X, "toarray") else np.asarray(cnmf.X)

print(f"Loaded cNMF modality: {cnmf.n_obs} cells x {cnmf.n_vars} programs")

# %%
prog_meta = pd.read_csv(REPO_ROOT / "data" / "program_activity_long_v4.csv").drop_duplicates(
    "program_name"
)

label_map = dict(
    zip(prog_meta["program_name"].astype(int), prog_meta["program_label"])
)
prog_color = dict(
    zip(prog_meta["program_name"].astype(int), prog_meta["distinct_color"])
)

program_order_path = REPO_ROOT / "analysis" / "regulator_dotplots" / "program_order.txt"
programs = []
for line in program_order_path.read_text().splitlines():
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    programs.append(int(line))

unique_programs = list(dict.fromkeys(programs))

print(f"Loaded {cnmf.n_obs} cells, {len(unique_programs)} programs")

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

pdf_path = OUT_DIR / "umap_all_cells_program_colors.pdf"
png_path = OUT_DIR / "umap_all_cells_program_colors.png"
fig.savefig(pdf_path, bbox_inches="tight")
fig.savefig(png_path, bbox_inches="tight")

print(f"Saved: {pdf_path}")
print(f"Saved: {png_path}")
plt.show()