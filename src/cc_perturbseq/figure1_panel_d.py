from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# -----------------------------------------------------------------------------
# Figure 1d: integrated UMAP with dominant selected cNMF program
# Source: programDE/notebooks/per_day_umaps/jesse_umap.ipynb
# -----------------------------------------------------------------------------


def align_umap_to_cnmf_cells(cnmf, adata_int, obsm_key: str = "X_umap") -> np.ndarray:
    cnmf_obs = pd.Index(cnmf.obs_names)
    adata_obs = pd.Index(adata_int.obs_names)
    if not cnmf_obs.is_unique or not adata_obs.is_unique:
        raise ValueError("Cell identifiers must be unique in both objects.")
    if obsm_key not in adata_int.obsm_keys():
        raise KeyError(f"{obsm_key!r} is absent from the integrated AnnData object.")
    missing = cnmf_obs[~cnmf_obs.isin(adata_obs)]
    extra = adata_obs[~adata_obs.isin(cnmf_obs)]
    if len(missing) or len(extra):
        raise ValueError(
            f"Cell-set mismatch: {len(missing)} missing from AnnData; "
            f"{len(extra)} extra in AnnData."
        )
    coordinates = np.asarray(adata_int.obsm[obsm_key])[:, :2]
    if cnmf_obs.equals(adata_obs):
        return coordinates
    positions = adata_obs.get_indexer(cnmf_obs)
    if np.any(positions < 0):
        raise ValueError("Unexpected unmatched cells after the cell-set check.")
    return coordinates[positions]


def make_figure_01d(
    h5mu_path: str | Path,
    adata_path: str | Path,
    metadata_path: str | Path,
    order_path: str | Path,
):
    import muon as mu
    import scanpy as sc

    h5mu_path = Path(h5mu_path)
    adata_path = Path(adata_path)
    metadata_path = Path(metadata_path)
    order_path = Path(order_path)

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
            "svg.fonttype": "none",
        }
    )

    mdata = mu.read_h5mu(h5mu_path)
    cnmf = mdata.mod["cNMF"]
    adata_int = sc.read_h5ad(adata_path, backed="r")
    coordinates = align_umap_to_cnmf_cells(cnmf, adata_int)
    scores = cnmf.X.toarray() if hasattr(cnmf.X, "toarray") else np.asarray(cnmf.X)
    if scores.shape[0] != coordinates.shape[0]:
        raise ValueError(
            f"cNMF/UMAP row mismatch: {scores.shape[0]} versus {coordinates.shape[0]}."
        )

    metadata = pd.read_csv(metadata_path).drop_duplicates("program_name")
    program_ids = metadata["program_name"].astype(int)
    label_map = (
        dict(zip(program_ids, metadata["program_label"]))
        if "program_label" in metadata.columns
        else {int(program): f"P{int(program)}" for program in program_ids}
    )
    color_column = next(
        (
            column
            for column in (
                "distinct_color",
                "program_color",
                "color",
                "class_color",
            )
            if column in metadata.columns
        ),
        None,
    )
    program_colors: dict[int, str] = {}
    if color_column is not None:
        for program, color in zip(program_ids, metadata[color_column]):
            if pd.isna(color):
                continue
            color = str(color).strip()
            try:
                mpl.colors.to_rgba(color)
            except ValueError:
                continue
            program_colors[int(program)] = color

    programs = []
    for line in order_path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            programs.append(int(line))
    unique_programs = list(dict.fromkeys(programs))
    if not unique_programs:
        raise ValueError("program_order.txt did not contain any program IDs.")
    if max(unique_programs) >= scores.shape[1] or min(unique_programs) < 0:
        raise IndexError(
            f"Program IDs {unique_programs} are incompatible with "
            f"{scores.shape[1]} cNMF columns."
        )

    fallback = plt.get_cmap("tab20", max(len(unique_programs), 1))
    for index, program in enumerate(unique_programs):
        program_colors.setdefault(program, mpl.colors.to_hex(fallback(index)))

    threshold = 0.02
    selected_scores = scores[:, unique_programs]
    maximum_index = np.argmax(selected_scores, axis=1)
    maximum_score = np.max(selected_scores, axis=1)
    dominant_program = np.asarray(
        [unique_programs[index] for index in maximum_index], dtype=int
    )
    above_threshold = maximum_score > threshold
    colors = np.asarray(["#e0e0e0"] * len(maximum_score), dtype=object)
    for program in unique_programs:
        mask = above_threshold & (dominant_program == program)
        colors[mask] = program_colors[program]

    fig, ax = plt.subplots(figsize=(5, 5))
    gray = ~above_threshold
    ax.scatter(
        coordinates[gray, 0],
        coordinates[gray, 1],
        c="#e8e8e8",
        s=0.1,
        rasterized=True,
        linewidths=0,
        zorder=1,
    )
    for program in unique_programs:
        mask = above_threshold & (dominant_program == program)
        if mask.any():
            ax.scatter(
                coordinates[mask, 0],
                coordinates[mask, 1],
                c=program_colors[program],
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
        short_label = (
            str(label).split(" — ")[-1] if " — " in str(label) else str(label)
        )
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=program_colors[program],
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

    source = pd.DataFrame(
        {
            "cell_id": cnmf.obs_names.astype(str),
            "umap_1": coordinates[:, 0],
            "umap_2": coordinates[:, 1],
            "dominant_selected_program": dominant_program,
            "maximum_selected_program_score": maximum_score,
            "above_0.02_threshold": above_threshold,
            "plot_color": colors,
        }
    )
    if getattr(adata_int, "file", None) is not None:
        adata_int.file.close()
    return fig, source, unique_programs, threshold
