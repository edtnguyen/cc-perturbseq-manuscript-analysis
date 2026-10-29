from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# -----------------------------------------------------------------------------
# Figure 1d: integrated UMAP with dominant selected cNMF program
# Source: programDE/notebooks/per_day_umaps/jesse_umap.ipynb
#
# The archived notebook contains an initial renderer that reads `distinct_color`
# from program metadata and a later manuscript renderer under "Try different
# pallette". The latter is the source of the final Figure 1 palette and program
# set. The previous repository port inadvertently retained the initial renderer.
# -----------------------------------------------------------------------------

FIGURE1D_THRESHOLD = 0.02

# Order follows the grouped labels in the manuscript Figure 1 program atlas.
FIGURE1D_MANUSCRIPT_PROGRAMS: tuple[int, ...] = (
    # Day 0
    22,
    29,
    7,
    # Day 1
    28,
    30,
    13,
    48,
    # Day 2
    15,
    17,
    14,
    19,
    23,
    # Day 3 arterial / endothelial cells
    2,
    3,
    8,
    11,
    42,
    # Day 3 alternative fates / non-EC
    24,
    33,
    34,
    27,
)

# Exact high-contrast, stage-constrained palette from the later manuscript
# renderer in jesse_umap.ipynb.
FIGURE1D_MANUSCRIPT_COLORS: dict[int, str] = {
    # Day 0: muted gray / pink
    29: "#4A4A4A",
    22: "#C8B4BC",
    7: "#8A8F96",
    # Day 1: green
    28: "#3E9A5F",
    13: "#86C86A",
    30: "#1F6F4A",
    48: "#4D9C72",
    # Day 2: orange / yellow / purple transition
    15: "#F39C12",
    23: "#F1C40F",
    14: "#B794C8",
    17: "#8E44AD",
    19: "#6F78D7",
    # Day 3 EC: red family
    42: "#8C2D04",
    8: "#C51B8A",
    11: "#F04E37",
    2: "#D73027",
    3: "#FB6A4A",
    # Day 3 non-EC: blue family
    27: "#2171B5",
    34: "#084594",
    33: "#9ECAE1",
    24: "#6BAED6",
}

# Labels used in the manuscript program list. Metadata remains the fallback for
# any future program not explicitly defined here.
FIGURE1D_MANUSCRIPT_LABELS: dict[int, str] = {
    22: "Pluripotency",
    29: "Pluripotency",
    7: "Pluripotency",
    28: "Primitive streak",
    30: "Primitive streak, BMP-repressed 1",
    13: "Posterior primitive streak",
    48: "Primitive streak, BMP-repressed 2",
    15: "Vascular progenitor",
    17: "Lat mesoderm, YAP/TAZ-induced",
    14: "Lateral mesoderm",
    19: "Posterior mesodermal",
    23: "Endothelial induction",
    2: "Artery endothelium 1",
    3: "Endothelial function 1",
    8: "Artery endothelium 2",
    11: "Endothelial function 2",
    42: "Basement membrane, adhesion",
    24: "Mesenchymal, HOX",
    33: "Anterior-like mesoderm",
    34: "Mesenchyme, YAP/TAZ-induced",
    27: "Mesenchyme, ECM",
}

# Programs explicitly re-drawn last in the recovered manuscript renderer so
# they remain visible in dense regions. This changes only draw order, never cell
# assignment.
FIGURE1D_TOP_OVERLAY_PROGRAMS: tuple[int, ...] = (7, 2, 34, 30)


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


def _parse_program_id(value: object) -> int:
    text = str(value).strip()
    if text.upper().startswith("P"):
        text = text[1:]
    return int(text)


def _program_ids_from_cnmf(cnmf, n_columns: int) -> list[int]:
    """Return explicit zero-based program IDs for cNMF score columns.

    The recovered source uses cNMF `var_names` when they form a complete P0..P49
    mapping and otherwise falls back to column position. This prevents a silent
    mismatch between program ID and matrix column.
    """
    try:
        parsed = [_parse_program_id(value) for value in cnmf.var_names]
    except (TypeError, ValueError):
        return list(range(n_columns))
    if len(parsed) != n_columns or len(set(parsed)) != n_columns:
        return list(range(n_columns))
    if set(parsed) != set(range(n_columns)):
        return list(range(n_columns))
    return parsed


def assign_figure1d_programs(
    scores: np.ndarray,
    *,
    program_ids: Sequence[int] | None = None,
    selected_programs: Sequence[int] = FIGURE1D_MANUSCRIPT_PROGRAMS,
    threshold: float = FIGURE1D_THRESHOLD,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Assign cells using the manuscript renderer's global-dominance rule.

    A cell is colored only when its highest-scoring program across the complete
    cNMF basis is one of the manuscript programs and the maximum score exceeds
    0.02. This is intentionally different from forcing every cell to the best
    program within a reduced selected-program subset.
    """
    matrix = np.asarray(scores)
    if matrix.ndim != 2 or matrix.shape[1] == 0:
        raise ValueError("scores must be a nonempty cell-by-program matrix.")
    if not np.isfinite(matrix).all():
        raise ValueError("cNMF scores contain non-finite values.")

    ids = list(range(matrix.shape[1])) if program_ids is None else [int(x) for x in program_ids]
    if len(ids) != matrix.shape[1] or len(set(ids)) != len(ids):
        raise ValueError("program_ids must uniquely identify every score column.")

    dominant_column = np.argmax(matrix, axis=1)
    maximum_score = matrix[np.arange(matrix.shape[0]), dominant_column]
    dominant_program = np.asarray([ids[index] for index in dominant_column], dtype=int)
    displayed = np.isin(dominant_program, np.asarray(selected_programs, dtype=int)) & (
        maximum_score > float(threshold)
    )
    return dominant_program, maximum_score, displayed


def _relative_luminance(color: str) -> float:
    rgb = np.asarray(mpl.colors.to_rgb(color), dtype=float)
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    return float(np.dot(linear, [0.2126, 0.7152, 0.0722]))


def _legend_edge(color: str) -> tuple[str, float]:
    luminance = _relative_luminance(color)
    if luminance >= 0.75:
        return "#222222", 1.20
    if luminance >= 0.45:
        return "#222222", 1.05
    if luminance >= 0.35:
        return "#222222", 0.85
    return "white", 0.65


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

    # Retain the path as a provenance input for the notebook, but do not use its
    # obsolete reduced program list. The final manuscript renderer defines the
    # complete 21-program set explicitly above.
    if not order_path.is_file():
        raise FileNotFoundError(order_path)

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
    try:
        coordinates = align_umap_to_cnmf_cells(cnmf, adata_int)
        scores = cnmf.X.toarray() if hasattr(cnmf.X, "toarray") else np.asarray(cnmf.X)
        if scores.shape[0] != coordinates.shape[0]:
            raise ValueError(
                f"cNMF/UMAP row mismatch: {scores.shape[0]} versus {coordinates.shape[0]}."
            )

        metadata = pd.read_csv(metadata_path).drop_duplicates("program_name")
        metadata["program_name"] = pd.to_numeric(
            metadata["program_name"], errors="raise"
        ).astype(int)
        metadata_labels = (
            dict(zip(metadata["program_name"], metadata["program_label"].astype(str)))
            if "program_label" in metadata.columns
            else {}
        )
        label_map = dict(metadata_labels)
        label_map.update(FIGURE1D_MANUSCRIPT_LABELS)

        programs = list(FIGURE1D_MANUSCRIPT_PROGRAMS)
        if set(programs) != set(FIGURE1D_MANUSCRIPT_COLORS):
            raise RuntimeError("The manuscript Figure 1d program and color sets disagree.")

        program_ids = _program_ids_from_cnmf(cnmf, scores.shape[1])
        missing = sorted(set(programs) - set(program_ids))
        if missing:
            raise ValueError(f"Manuscript programs are missing from cNMF scores: {missing}")

        dominant_program, maximum_score, displayed = assign_figure1d_programs(
            scores,
            program_ids=program_ids,
            selected_programs=programs,
            threshold=FIGURE1D_THRESHOLD,
        )

        background_color = "#666A70"
        cell_colors = np.asarray([background_color] * len(maximum_score), dtype=object)
        for program in programs:
            mask = displayed & (dominant_program == program)
            cell_colors[mask] = FIGURE1D_MANUSCRIPT_COLORS[program]

        fig, ax = plt.subplots(figsize=(5.2, 5.2))
        background = ~displayed
        ax.scatter(
            coordinates[background, 0],
            coordinates[background, 1],
            c=background_color,
            s=0.060,
            alpha=0.55,
            rasterized=True,
            linewidths=0,
            edgecolors="none",
            zorder=1,
        )

        foreground_index = np.flatnonzero(displayed)
        if len(foreground_index):
            overlay_set = set(FIGURE1D_TOP_OVERLAY_PROGRAMS)
            foreground_luminance = np.asarray(
                [
                    _relative_luminance(
                        FIGURE1D_MANUSCRIPT_COLORS[int(dominant_program[index])]
                    )
                    for index in foreground_index
                ],
                dtype=float,
            )
            overlay = np.asarray(
                [int(dominant_program[index]) in overlay_set for index in foreground_index],
                dtype=bool,
            )
            ordinary_index = foreground_index[~overlay]
            ordinary_luminance = foreground_luminance[~overlay]
            ordinary_index = ordinary_index[
                np.argsort(ordinary_luminance, kind="stable")
            ]
            overlay_index = foreground_index[overlay]
            foreground_index = np.concatenate([ordinary_index, overlay_index])
            foreground_colors = [
                FIGURE1D_MANUSCRIPT_COLORS[int(dominant_program[index])]
                for index in foreground_index
            ]
            ax.scatter(
                coordinates[foreground_index, 0],
                coordinates[foreground_index, 1],
                c=foreground_colors,
                s=0.62,
                alpha=1.0,
                rasterized=True,
                linewidths=0,
                edgecolors="none",
                zorder=3,
            )

            for program in FIGURE1D_TOP_OVERLAY_PROGRAMS:
                mask = displayed & (dominant_program == program)
                if not mask.any():
                    continue
                ax.scatter(
                    coordinates[mask, 0],
                    coordinates[mask, 1],
                    c=FIGURE1D_MANUSCRIPT_COLORS[program],
                    s=0.62 * 1.25,
                    alpha=1.0,
                    rasterized=True,
                    linewidths=0,
                    edgecolors="none",
                    zorder=4,
                )

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        for spine in ax.spines.values():
            spine.set_visible(False)

        handles = []
        for program in programs:
            color = FIGURE1D_MANUSCRIPT_COLORS[program]
            edge_color, edge_width = _legend_edge(color)
            handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="none",
                    markerfacecolor=color,
                    markeredgecolor=edge_color,
                    markeredgewidth=edge_width,
                    markersize=3.2,
                    label=f"P{program} {label_map.get(program, f'P{program}')}",
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

        selected_program = pd.array(
            [int(value) if keep else pd.NA for value, keep in zip(dominant_program, displayed)],
            dtype="Int64",
        )
        selected_score = np.where(displayed, maximum_score, np.nan)
        source = pd.DataFrame(
            {
                "cell_id": cnmf.obs_names.astype(str),
                "umap_1": coordinates[:, 0],
                "umap_2": coordinates[:, 1],
                # Compatibility fields used by the current Figure_01 notebook.
                "dominant_selected_program": selected_program,
                "maximum_selected_program_score": selected_score,
                "above_0.02_threshold": displayed,
                # Explicit manuscript semantics.
                "dominant_global_program": dominant_program,
                "maximum_global_program_score": maximum_score,
                "dominant_program_is_manuscript_program": np.isin(
                    dominant_program, programs
                ),
                "plot_color": cell_colors,
            }
        )
        return fig, source, programs, FIGURE1D_THRESHOLD
    finally:
        if getattr(adata_int, "file", None) is not None:
            adata_int.file.close()
