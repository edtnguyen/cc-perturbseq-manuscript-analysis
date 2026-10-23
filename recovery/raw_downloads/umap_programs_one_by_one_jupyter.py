# %% [markdown]
# # Display programs from the metadata CSV one by one on UMAP
#
# This notebook plots each program listed in `data/program_activity_long_v4.csv` on the same UMAP.
#
# Default behavior:
# - uses cNMF scores from the `.h5mu` file;
# - uses `X_umap` transferred from the integrated `.h5ad` file, aligned by cell index;
# - displays every program one by one;
# - saves compact SVG/PNG files for each program.

# %%
import re
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# %%
def find_repo_root(start: Path = Path.cwd()) -> Path:
    """Find the repository root by looking upward for the program metadata CSV."""
    start = start.resolve()
    required_paths = [Path("data/program_activity_long_v4.csv")]

    for candidate in [start, *start.parents]:
        if all((candidate / rel_path).exists() for rel_path in required_paths):
            return candidate

    raise FileNotFoundError(
        "Could not auto-detect REPO_ROOT. Set REPO_ROOT manually in this cell."
    )


# In a notebook, __file__ is usually not defined, so use an auto-detected repo root.
# If auto-detection fails, replace this with: REPO_ROOT = Path("/path/to/your/repo")
REPO_ROOT = find_repo_root()

OUT_DIR = REPO_ROOT / "analysis" / "regulator_dotplots" / "figures" / "umap_programs_one_by_one"
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
        # Keep labels/text editable in Illustrator while rasterizing dense point layers.
        "svg.fonttype": "none",
    }
)

print(f"REPO_ROOT = {REPO_ROOT}")
print(f"OUT_DIR   = {OUT_DIR}")

# %%
H5MU_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/adata/cNMF_50_0_2.h5mu"
)

# Set to False to use cnmf.obsm["X_umap"] instead of the integrated h5ad UMAP.
USE_ADATA_INT_UMAP = True
ADATA_INT_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
    "2.preprocessing_tonys_pipeline/3.integrated_days/"
    "3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad"
)

# If True, error when adata_int has cells not present in cnmf or vice versa.
# If False, adata_int may contain extra cells, but every cnmf cell must still be present.
REQUIRE_EXACT_OBS_MATCH = True

PROGRAM_CSV_PATH = REPO_ROOT / "data" / "program_activity_long_v4.csv"

# %%
def dense_array(x):
    """Convert sparse/dask-like arrays to a NumPy array when needed."""
    return x.toarray() if hasattr(x, "toarray") else np.asarray(x)


def transfer_umap_from_adata_int(cnmf, adata_int, *, require_exact_match: bool = True) -> np.ndarray:
    """Return adata_int X_umap reordered to match cnmf.obs_names."""
    if "X_umap" not in adata_int.obsm:
        raise KeyError("adata_int.obsm does not contain 'X_umap'.")

    cnmf_index = pd.Index(cnmf.obs_names, name="cnmf_obs_names")
    adata_index = pd.Index(adata_int.obs_names, name="adata_int_obs_names")

    if not cnmf_index.is_unique:
        duplicated = cnmf_index[cnmf_index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"cnmf.obs_names is not unique. Examples: {duplicated}")
    if not adata_index.is_unique:
        duplicated = adata_index[adata_index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"adata_int.obs_names is not unique. Examples: {duplicated}")

    missing_in_adata = cnmf_index.difference(adata_index)
    extra_in_adata = adata_index.difference(cnmf_index)

    if len(missing_in_adata) > 0:
        examples = missing_in_adata[:10].tolist()
        raise ValueError(
            f"{len(missing_in_adata):,} cnmf cells are missing from adata_int. "
            f"Examples: {examples}"
        )

    if require_exact_match and len(extra_in_adata) > 0:
        examples = extra_in_adata[:10].tolist()
        raise ValueError(
            f"adata_int has {len(extra_in_adata):,} extra cells not present in cnmf. "
            f"Examples: {examples}. If this is expected, set REQUIRE_EXACT_OBS_MATCH = False."
        )

    positions = adata_index.get_indexer(cnmf_index)
    if np.any(positions < 0):
        raise RuntimeError("Unexpected failed indexer values after missing-cell check.")

    adata_umap = np.asarray(adata_int.obsm["X_umap"])
    umap_coords = adata_umap[positions, :]

    if umap_coords.shape[0] != cnmf.n_obs or umap_coords.shape[1] != 2:
        raise ValueError(
            "Transferred UMAP has unexpected shape: "
            f"{umap_coords.shape}; expected ({cnmf.n_obs}, 2)."
        )

    print("Transferred X_umap from adata_int to cnmf order")
    print(f"  cnmf cells:      {len(cnmf_index):,}")
    print(f"  adata_int cells: {len(adata_index):,}")
    print(f"  extra adata_int cells ignored: {len(extra_in_adata):,}")
    return umap_coords


if not H5MU_PATH.exists():
    raise FileNotFoundError(f"Could not find h5mu file: {H5MU_PATH}")

print("Loading h5mu...")
mdata = mu.read_h5mu(H5MU_PATH)
cnmf = mdata.mod["cNMF"]
scores = dense_array(cnmf.X)

if USE_ADATA_INT_UMAP:
    if not ADATA_INT_PATH.exists():
        raise FileNotFoundError(f"Could not find integrated h5ad file: {ADATA_INT_PATH}")
    import scanpy as sc

    print("Loading integrated h5ad for UMAP transfer...")
    adata_int = sc.read_h5ad(ADATA_INT_PATH, backed="r")
    umap_coords = transfer_umap_from_adata_int(
        cnmf,
        adata_int,
        require_exact_match=REQUIRE_EXACT_OBS_MATCH,
    )
    cnmf.obsm["X_umap_from_adata_int"] = umap_coords
    try:
        adata_int.file.close()
    except Exception:
        pass
else:
    if "X_umap" not in cnmf.obsm:
        raise KeyError("cnmf.obsm does not contain 'X_umap'.")
    umap_coords = np.asarray(cnmf.obsm["X_umap"])

print(f"Loaded cNMF modality: {cnmf.n_obs:,} cells x {cnmf.n_vars:,} programs")
print(f"UMAP coordinates shape: {umap_coords.shape}")
print(f"Scores shape: {scores.shape}")

# %%
def coerce_program_id(value):
    """Convert program_name values such as 3, 3.0, or 'P3' to integer IDs."""
    if pd.isna(value):
        return None
    try:
        return int(value)
    except Exception:
        match = re.search(r"\d+", str(value))
        return int(match.group(0)) if match else None


def is_valid_color(value) -> bool:
    if pd.isna(value):
        return False
    try:
        mpl.colors.to_hex(value)
        return True
    except ValueError:
        return False


def safe_filename(value: str, max_len: int = 90) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return cleaned[:max_len] if cleaned else "program"


def short_label(label: str) -> str:
    label = str(label)
    return label.split(" — ")[-1] if " — " in label else label


prog_meta_raw = pd.read_csv(PROGRAM_CSV_PATH)
if "program_name" not in prog_meta_raw.columns:
    raise KeyError(
        f"Expected a 'program_name' column in {PROGRAM_CSV_PATH}. "
        f"Available columns: {list(prog_meta_raw.columns)}"
    )

prog_meta = prog_meta_raw.copy()
prog_meta["program_id"] = prog_meta["program_name"].map(coerce_program_id)
prog_meta = prog_meta.dropna(subset=["program_id"]).copy()
prog_meta["program_id"] = prog_meta["program_id"].astype(int)
prog_meta = prog_meta.drop_duplicates("program_id", keep="first").reset_index(drop=True)

programs_from_csv = prog_meta["program_id"].tolist()

if "program_label" in prog_meta.columns:
    label_map = dict(zip(prog_meta["program_id"], prog_meta["program_label"].astype(str)))
else:
    label_map = {program: f"P{program}" for program in programs_from_csv}

preferred_color_columns = ["distinct_color", "program_color", "color", "class_color"]
color_column = next((col for col in preferred_color_columns if col in prog_meta.columns), None)

prog_color = {}
if color_column is not None:
    for program, color in zip(prog_meta["program_id"], prog_meta[color_column]):
        if is_valid_color(color):
            prog_color[int(program)] = mpl.colors.to_hex(color)
    print(f"Using '{color_column}' from program metadata for program colors")
else:
    print("No color column found in program metadata; using fallback tab20 colors")

fallback_cmap = plt.get_cmap("tab20")
for idx, program in enumerate(programs_from_csv):
    if program not in prog_color:
        prog_color[program] = mpl.colors.to_hex(fallback_cmap(idx % fallback_cmap.N))

print(f"Loaded {len(programs_from_csv):,} unique programs from {PROGRAM_CSV_PATH}")
print(f"Program metadata columns: {list(prog_meta_raw.columns)}")

# %%
def resolve_program_columns(program_ids, cnmf, scores) -> tuple[dict[int, int], list[int]]:
    """Map program IDs from the CSV to columns in cnmf.X."""
    n_program_cols = scores.shape[1]
    var_names = [str(name) for name in getattr(cnmf, "var_names", [])]
    var_lookup = {name: idx for idx, name in enumerate(var_names)}

    program_to_col = {}
    skipped = []

    for program in program_ids:
        candidates = [
            str(program),
            f"P{program}",
            f"p{program}",
            f"program_{program}",
            f"Program_{program}",
            f"cNMF_{program}",
        ]

        matched_col = None
        for candidate in candidates:
            if candidate in var_lookup:
                matched_col = var_lookup[candidate]
                break

        # Fallback: previous scripts used program IDs directly as integer column indices.
        if matched_col is None and 0 <= int(program) < n_program_cols:
            matched_col = int(program)

        if matched_col is None:
            skipped.append(int(program))
        else:
            program_to_col[int(program)] = int(matched_col)

    return program_to_col, skipped


program_to_col, skipped_programs = resolve_program_columns(programs_from_csv, cnmf, scores)
programs = [program for program in programs_from_csv if program in program_to_col]

if len(programs) == 0:
    raise ValueError(
        "None of the programs from the CSV could be mapped to columns in cnmf.X. "
        "Check whether program_name is 0-based, 1-based, or stored in cnmf.var_names."
    )

if skipped_programs:
    print(
        f"Skipped {len(skipped_programs):,} programs that could not be mapped to cnmf.X columns: "
        f"{skipped_programs[:20]}"
    )

print(f"Will plot {len(programs):,} programs from the CSV")
print("First mapped programs:")
for program in programs[:10]:
    print(f"  P{program} -> scores[:, {program_to_col[program]}] | {label_map.get(program, f'P{program}')}")

# %%
# Plot controls.
PLOT_MODE = "score"  # "score" or "dominant"

# Score mode: plot the activity score of one program across cells.
SCORE_MIN = 0.0
SCORE_VMAX_QUANTILE = 0.99

# Dominant mode: highlight cells for which the current program is the highest-scoring program.
DOMINANT_THRESHOLD = 0.02

FIGSIZE = (4.2, 4.2)
BACKGROUND_SIZE = 0.12
POINT_SIZE = 0.35
BACKGROUND_COLOR = "#e8e8e8"
HIGHLIGHT_ALPHA = 0.85

SAVE_FORMATS = ("svg", "png")

# Precompute dominant program among the CSV-listed programs for optional dominant mode.
program_cols = [program_to_col[program] for program in programs]
sub_scores = scores[:, program_cols]
max_program_idx = np.argmax(sub_scores, axis=1)
max_score = np.max(sub_scores, axis=1)
dominant_program = np.array([programs[idx] for idx in max_program_idx])

print("Plot controls loaded")

# %%
def make_program_cmap(program: int):
    """Make a light-to-program-color colormap for score plots."""
    color = prog_color.get(program, "#377eb8")
    if not is_valid_color(color):
        color = "#377eb8"
    return mpl.colors.LinearSegmentedColormap.from_list(
        f"program_{program}_cmap",
        ["#f0f0f0", mpl.colors.to_hex(color)],
    )


def style_umap_axis(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_program(
    program: int,
    *,
    mode: str = PLOT_MODE,
    save: bool = False,
    show: bool = True,
    close: bool = False,
):
    """Plot one program on the UMAP.

    Parameters
    ----------
    program
        Integer program ID from `program_activity_long_v4.csv`.
    mode
        "score" plots continuous cNMF score for the program.
        "dominant" highlights cells where this program is the dominant program.
    save
        Save files into OUT_DIR using SAVE_FORMATS.
    show
        Display the figure in the notebook.
    close
        Close the figure after display/save to reduce memory use in loops.
    """
    program = int(program)
    if program not in program_to_col:
        raise KeyError(f"Program P{program} is not mapped to a cnmf.X column.")

    label = label_map.get(program, f"P{program}")
    label_short = short_label(label)
    color = prog_color.get(program, "#377eb8")

    fig, ax = plt.subplots(figsize=FIGSIZE)

    ax.scatter(
        umap_coords[:, 0],
        umap_coords[:, 1],
        c=BACKGROUND_COLOR,
        s=BACKGROUND_SIZE,
        rasterized=True,
        linewidths=0,
        zorder=1,
    )

    if mode == "score":
        score_values = np.asarray(scores[:, program_to_col[program]], dtype=float)
        mask = np.isfinite(score_values) & (score_values > SCORE_MIN)
        n_cells = int(mask.sum())

        if n_cells > 0:
            mask_idx = np.flatnonzero(mask)
            order = mask_idx[np.argsort(score_values[mask_idx])]
            positive_scores = score_values[mask]
            vmax = float(np.quantile(positive_scores, SCORE_VMAX_QUANTILE))
            if not np.isfinite(vmax) or vmax <= SCORE_MIN:
                vmax = float(np.nanmax(positive_scores))
            if not np.isfinite(vmax) or vmax <= SCORE_MIN:
                vmax = SCORE_MIN + 1e-12

            scatter = ax.scatter(
                umap_coords[order, 0],
                umap_coords[order, 1],
                c=score_values[order],
                cmap=make_program_cmap(program),
                vmin=SCORE_MIN,
                vmax=vmax,
                s=POINT_SIZE,
                rasterized=True,
                linewidths=0,
                zorder=2,
                alpha=HIGHLIGHT_ALPHA,
            )
            cbar = fig.colorbar(scatter, ax=ax, fraction=0.045, pad=0.015)
            cbar.set_label("cNMF score", fontsize=6)
            cbar.ax.tick_params(labelsize=5, length=2)

        title_suffix = f"score > {SCORE_MIN:g}: {n_cells:,} cells"

    elif mode == "dominant":
        mask = (max_score > DOMINANT_THRESHOLD) & (dominant_program == program)
        n_cells = int(mask.sum())

        if n_cells > 0:
            ax.scatter(
                umap_coords[mask, 0],
                umap_coords[mask, 1],
                c=color,
                s=POINT_SIZE,
                rasterized=True,
                linewidths=0,
                zorder=2,
                alpha=HIGHLIGHT_ALPHA,
            )

        title_suffix = f"dominant cells > {DOMINANT_THRESHOLD:g}: {n_cells:,}"

    else:
        raise ValueError("mode must be either 'score' or 'dominant'.")

    style_umap_axis(ax)
    ax.set_title(f"P{program} {label_short}\n{title_suffix}", fontsize=7)
    plt.tight_layout()

    saved_paths = []
    if save:
        file_stem = f"umap_P{program:02d}_{safe_filename(label_short)}_{mode}"
        for file_format in SAVE_FORMATS:
            path = OUT_DIR / f"{file_stem}.{file_format}"
            fig.savefig(path, format=file_format, bbox_inches="tight")
            saved_paths.append(path)

    if show:
        plt.show()

    if close:
        plt.close(fig)

    return fig, ax, saved_paths


# Preview the first program.
plot_program(programs[0], mode=PLOT_MODE, save=False, show=True, close=False)

# %% [markdown]
# ## Optional: interactive single-program viewer
#
# This cell shows one program at a time using a dropdown. If `ipywidgets` is unavailable, use `plot_program(program_id)` directly.

# %%
try:
    import ipywidgets as widgets
    from IPython.display import clear_output, display

    program_options = [
        (f"P{program} {short_label(label_map.get(program, f'P{program}'))}", program)
        for program in programs
    ]

    program_dropdown = widgets.Dropdown(
        options=program_options,
        value=programs[0],
        description="Program",
        layout=widgets.Layout(width="420px"),
    )
    mode_dropdown = widgets.ToggleButtons(
        options=["score", "dominant"],
        value=PLOT_MODE,
        description="Mode",
    )
    output = widgets.Output()

    def update_program_plot(*_):
        with output:
            clear_output(wait=True)
            plot_program(
                program_dropdown.value,
                mode=mode_dropdown.value,
                save=False,
                show=True,
                close=True,
            )

    program_dropdown.observe(update_program_plot, names="value")
    mode_dropdown.observe(update_program_plot, names="value")
    display(program_dropdown, mode_dropdown, output)
    update_program_plot()

except Exception as exc:
    print(f"ipywidgets viewer not available: {exc}")
    print("Use plot_program(program_id, mode='score') instead.")

# %% [markdown]
# ## Display and save every CSV-listed program one by one
#
# Running this cell will show every program in the notebook and save per-program SVG/PNG files.
# Dense point layers are rasterized inside SVG files to keep them compact; text remains editable.

# %%
RUN_ALL_PROGRAMS = True
SAVE_EACH_PROGRAM = True
DISPLAY_EACH_PROGRAM = True

all_saved_paths = []
if RUN_ALL_PROGRAMS:
    for program in programs:
        _, _, saved_paths = plot_program(
            program,
            mode=PLOT_MODE,
            save=SAVE_EACH_PROGRAM,
            show=DISPLAY_EACH_PROGRAM,
            close=True,
        )
        all_saved_paths.extend(saved_paths)

    print(f"Finished plotting {len(programs):,} programs")
    if SAVE_EACH_PROGRAM:
        print(f"Saved {len(all_saved_paths):,} files to: {OUT_DIR}")
else:
    print("RUN_ALL_PROGRAMS is False. Set it to True to display/save every program.")
