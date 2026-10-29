from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from typing import Any, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


# -----------------------------------------------------------------------------
# Source-faithful feature UMAPs used by Figures 2e, 3b and 4d.
#
# Recovered source:
#   programDE/notebooks/per_day_umaps/jesse_umap.ipynb
#   later block: "Seurat-style individual UMAPs for selected genes and cNMF
#   programs"
#
# The source uses:
#   * integrated AnnData X_umap after exact cell-ID alignment;
#   * cNMF program columns by zero-based positional ID, i.e. cnmf.X[:, P];
#   * expression from integrated AnnData X, summing duplicate gene symbols;
#   * a light-gray -> dark-gray color scale;
#   * no lower clipping and upper-tail clipping only at p99;
#   * ascending draw order so the highest-value cells are plotted last;
#   * point size 0.30 and a separate colorbar for every feature.
# -----------------------------------------------------------------------------

SOURCE_UPPER_PERCENTILE = 99.0
SOURCE_POINT_SIZE = 0.30
SOURCE_INDIVIDUAL_FIGSIZE = (3.25, 3.0)
SOURCE_CMAP_COLORS = ("#D3D3D3", "#4D4D4D")
SOURCE_CMAP = mpl.colors.LinearSegmentedColormap.from_list(
    "seurat_grayblue",
    list(SOURCE_CMAP_COLORS),
)
SOURCE_GENE_ALIASES = {
    # The recovered notebook explicitly accepts this historical typo.
    "ACT1C1": "ACTC1",
}

# Executed values printed by the recovered source notebook. They are retained as
# provenance checks for a full Sherlock execution, not as hard-coded scaling.
SOURCE_EXECUTED_P99 = {
    ("gene", "ACTC1"): 3.287982,
    ("gene", "GATA4"): 2.463099,
    ("gene", "SLC8A1"): 2.587773,
    ("gene", "FLI1"): 3.283608,
    ("gene", "CDH5"): 3.359469,
    ("program", 2): 0.402983,
    ("program", 34): 0.184445,
    ("program", 15): 0.347766,
    ("program", 17): 0.331715,
    ("program", 27): 0.359415,
}


@dataclass
class SourceMatplotlibPlot:
    """Matplotlib figure with the small plotnine-like save surface notebooks use."""

    figure: Any

    def save(
        self,
        filename: str | Path,
        *,
        dpi: int = 300,
        verbose: bool = False,
        limitsize: bool = False,
        **kwargs: Any,
    ) -> None:
        del verbose, limitsize
        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.figure.savefig(
            path,
            dpi=dpi,
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
            **kwargs,
        )

    def draw(self, show: bool = False):
        if show:
            plt.show()
        return self.figure

    def _repr_png_(self) -> bytes:
        buffer = BytesIO()
        self.figure.savefig(
            buffer,
            format="png",
            dpi=150,
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
        )
        return buffer.getvalue()


def upper_clip_only(
    values: Sequence[float] | np.ndarray,
    percentile: float = SOURCE_UPPER_PERCENTILE,
) -> tuple[np.ndarray, float, float]:
    """Match the recovered source's upper-tail-only clipping exactly."""
    vector = np.asarray(values, dtype=float).reshape(-1)
    finite = vector[np.isfinite(vector)]
    if finite.size == 0:
        return np.zeros_like(vector), 0.0, 1.0

    upper = float(np.percentile(finite, percentile))
    if not np.isfinite(upper):
        upper = float(np.nanmax(finite))
    lower = float(np.nanmin(finite))
    if upper <= lower:
        upper = float(np.nanmax(finite))
        if upper <= lower:
            upper = lower + 1e-12

    plotted = np.minimum(vector, upper)
    plotted[~np.isfinite(plotted)] = lower
    return plotted, lower, upper


def source_program_vector(
    usage: pd.DataFrame | np.ndarray,
    program_id: int,
) -> tuple[np.ndarray, str]:
    """Read a program by zero-based column position, as cnmf.X[:, program_id].

    The original source does not resolve a label such as ``P17``. It indexes the
    cNMF matrix positionally. This distinction matters when a text export uses
    one-based labels while retaining zero-based matrix order.
    """
    position = int(program_id)
    matrix = usage.to_numpy() if isinstance(usage, pd.DataFrame) else np.asarray(usage)
    if matrix.ndim != 2:
        raise ValueError("cNMF usage must be a two-dimensional cell-by-program matrix.")
    if position < 0 or position >= matrix.shape[1]:
        raise IndexError(
            f"Program P{position} is outside a cNMF matrix with {matrix.shape[1]} columns."
        )
    vector = np.asarray(matrix[:, position], dtype=float).reshape(-1)
    if not np.isfinite(vector).all():
        raise ValueError(f"Program P{position} contains missing or non-finite values.")
    source_column = (
        str(usage.columns[position])
        if isinstance(usage, pd.DataFrame)
        else f"column_{position}"
    )
    return vector, source_column


def _candidate_gene_columns(adata) -> list[pd.Index]:
    candidates = [pd.Index(adata.var_names.astype(str))]
    for column in ("gene", "gene_symbol", "gene_name", "symbol", "feature_name"):
        if column in adata.var.columns:
            candidates.append(pd.Index(adata.var[column].astype(str)))
    return candidates


def _find_gene_positions(adata, gene: str) -> tuple[np.ndarray, str]:
    requested = str(gene)
    aliases = [requested]
    canonical = SOURCE_GENE_ALIASES.get(requested.upper())
    if canonical and canonical not in aliases:
        aliases.append(canonical)

    candidates = _candidate_gene_columns(adata)
    for alias in aliases:
        exact: list[int] = []
        for candidate in candidates:
            exact.extend(np.flatnonzero(candidate.to_numpy() == alias).tolist())
        if exact:
            return np.asarray(sorted(set(exact)), dtype=int), alias

    for alias in aliases:
        alias_upper = alias.upper()
        matches: list[int] = []
        for candidate in candidates:
            matches.extend(
                np.flatnonzero(candidate.str.upper().to_numpy() == alias_upper).tolist()
            )
        if matches:
            return np.asarray(sorted(set(matches)), dtype=int), alias

    raise KeyError(
        f"Gene {requested!r} was not found in var_names or supported gene-symbol columns."
    )


def source_gene_vector(adata, gene: str) -> tuple[np.ndarray, str, np.ndarray]:
    """Match the recovered source: expression from X, summing duplicate symbols."""
    positions, matched_name = _find_gene_positions(adata, gene)
    matrix = adata[:, positions].X
    if sparse.issparse(matrix):
        vector = np.asarray(matrix.sum(axis=1)).ravel()
    else:
        vector = np.asarray(matrix).sum(axis=1)
    vector = np.asarray(vector, dtype=float).reshape(-1)
    if not np.isfinite(vector).all():
        raise ValueError(f"Expression values for {gene!r} contain non-finite values.")
    return vector, matched_name, positions


def _validate_alignment(adata, usage: pd.DataFrame | np.ndarray) -> None:
    n_usage = usage.shape[0]
    if n_usage != adata.n_obs:
        raise ValueError(
            f"cNMF/AnnData row mismatch: {n_usage} usage rows versus {adata.n_obs} cells."
        )
    if isinstance(usage, pd.DataFrame):
        usage_index = pd.Index(usage.index.astype(str))
        adata_index = pd.Index(adata.obs_names.astype(str))
        if not usage_index.equals(adata_index):
            missing = adata_index[~adata_index.isin(usage_index)]
            extra = usage_index[~usage_index.isin(adata_index)]
            raise ValueError(
                "cNMF usage and AnnData cells are not in the exact same order "
                f"({len(missing)} missing, {len(extra)} extra)."
            )


def _feature_vector(
    adata,
    usage: pd.DataFrame | np.ndarray,
    kind: str,
    feature: int | str,
) -> tuple[np.ndarray, dict[str, Any]]:
    kind_clean = str(kind).lower()
    if kind_clean == "program":
        program_id = int(feature)
        values, source_column = source_program_vector(usage, program_id)
        return values, {
            "feature_kind": "program",
            "feature": f"P{program_id}",
            "program_id_zero_based": program_id,
            "source_usage_column": source_column,
            "source_semantics": f"positional cNMF column {program_id}",
        }
    if kind_clean == "gene":
        values, matched_name, positions = source_gene_vector(adata, str(feature))
        return values, {
            "feature_kind": "gene",
            "feature": str(feature),
            "matched_gene_name": matched_name,
            "matched_feature_positions": ",".join(map(str, positions.tolist())),
            "source_semantics": "integrated AnnData X; duplicate gene symbols summed",
        }
    raise ValueError(f"Unknown feature kind: {kind!r}")


def render_source_feature_umap_grid(
    adata,
    usage: pd.DataFrame | np.ndarray,
    specs: Sequence[tuple[str, int | str, str]],
    *,
    upper_percentile: float = SOURCE_UPPER_PERCENTILE,
    point_size: float = SOURCE_POINT_SIZE,
    ncols: int = 2,
    figure_size: tuple[float, float] | None = None,
) -> tuple[SourceMatplotlibPlot, dict[str, pd.DataFrame]]:
    """Render a composite from the recovered individual-UMAP source grammar."""
    if "X_umap" not in adata.obsm:
        raise KeyError("AnnData lacks obsm['X_umap'].")
    _validate_alignment(adata, usage)

    coordinates = np.asarray(adata.obsm["X_umap"], dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError(f"X_umap must be n_cells x >=2; found {coordinates.shape}.")
    coordinates = coordinates[:, :2]
    if not np.isfinite(coordinates).all():
        raise ValueError("X_umap contains missing or non-finite coordinates.")

    n_features = len(specs)
    if n_features == 0:
        raise ValueError("At least one feature specification is required.")
    ncols = max(1, min(int(ncols), n_features))
    nrows = int(np.ceil(n_features / ncols))
    if figure_size is None:
        figure_size = (
            SOURCE_INDIVIDUAL_FIGSIZE[0] * ncols,
            SOURCE_INDIVIDUAL_FIGSIZE[1] * nrows,
        )

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    figure, axes = plt.subplots(
        nrows,
        ncols,
        figsize=figure_size,
        squeeze=False,
        constrained_layout=True,
    )
    source: dict[str, pd.DataFrame] = {}

    for axis, (kind, feature, title) in zip(axes.ravel(), specs):
        values, provenance = _feature_vector(adata, usage, kind, feature)
        plotted, lower, upper = upper_clip_only(values, upper_percentile)
        order = np.argsort(plotted, kind="stable")
        scatter = axis.scatter(
            coordinates[order, 0],
            coordinates[order, 1],
            c=plotted[order],
            cmap=SOURCE_CMAP,
            vmin=lower,
            vmax=upper,
            s=point_size,
            linewidths=0,
            edgecolors="none",
            rasterized=True,
        )
        axis.set_aspect("equal")
        axis.axis("off")
        axis.set_title(str(title), fontsize=11, fontweight="bold", pad=3)

        colorbar = figure.colorbar(
            scatter,
            ax=axis,
            fraction=0.035,
            pad=0.015,
            aspect=25,
        )
        legend_title = (
            "cNMF usage" if str(kind).lower() == "program" else "X expression"
        )
        colorbar.set_label(legend_title, fontsize=7)
        colorbar.ax.tick_params(labelsize=6, length=2)
        colorbar.outline.set_linewidth(0.5)
        axis.text(
            0.02,
            0.02,
            f"upper clip: p{upper_percentile:g} = {upper:.3g}",
            transform=axis.transAxes,
            fontsize=5.5,
            color="#4D4D4D",
            ha="left",
            va="bottom",
        )

        key = (
            f"program_{int(feature)}"
            if str(kind).lower() == "program"
            else f"gene_{feature}"
        )
        display_fraction = np.divide(
            plotted - lower,
            upper - lower,
            out=np.zeros_like(plotted),
            where=(upper - lower) > 0,
        )
        table = pd.DataFrame(
            {
                "umap_1": coordinates[:, 0],
                "umap_2": coordinates[:, 1],
                "value": values,
                "display_value": plotted,
                "display_vmin": lower,
                "display_vmax_p99": upper,
                "display_fraction": display_fraction,
                "draw_rank": np.argsort(order, kind="stable"),
            },
            index=adata.obs_names,
        )
        for name, value in provenance.items():
            table[name] = value
        table["upper_clip_percentile"] = float(upper_percentile)
        table["colormap_low"] = SOURCE_CMAP_COLORS[0]
        table["colormap_high"] = SOURCE_CMAP_COLORS[1]
        source[key] = table

    for axis in axes.ravel()[n_features:]:
        axis.set_visible(False)

    return SourceMatplotlibPlot(figure), source


FIGURE2E_SPECS = (
    ("program", 17, "Programs\nP17\nLat. mesoderm, YAP/TAZ-induced"),
    ("gene", "SLC8A1", "Markers\nSLC8A1"),
    ("program", 27, "Programs\nP27\nMesenchyme, ECM"),
    ("gene", "GATA4", "Markers\nGATA4"),
)

FIGURE3B_SPECS = (
    ("program", 2, "Programs\nP2\nArtery endothelium 1"),
    ("gene", "CDH5", "Markers\nCDH5"),
    ("program", 34, "Programs\nP34\nMesenchyme, YAP/TAZ-induced"),
    ("gene", "ACTC1", "Markers\nACTC1"),
)

FIGURE4D_SPECS = (
    ("program", 15, "Programs\nP15\nVascular progenitor"),
    ("gene", "FLI1", "Markers\nFLI1"),
    ("program", 17, "Programs\nP17\nLat. mesoderm, YAP/TAZ-induced"),
    ("gene", "SLC8A1", "Markers\nSLC8A1"),
)


def plot_figure2e_umaps(adata, usage: pd.DataFrame, **_: Any):
    return render_source_feature_umap_grid(adata, usage, FIGURE2E_SPECS)


def plot_figure3b_feature_umaps(adata, usage: pd.DataFrame, **_: Any):
    return render_source_feature_umap_grid(adata, usage, FIGURE3B_SPECS)


def plot_figure4d_feature_umaps(adata, usage: pd.DataFrame, **_: Any):
    return render_source_feature_umap_grid(adata, usage, FIGURE4D_SPECS)
