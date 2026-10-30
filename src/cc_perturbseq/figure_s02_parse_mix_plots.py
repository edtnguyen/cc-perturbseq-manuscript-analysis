"""Plot producers for Figure S2 mixed hiPSC–K562 panels."""
from __future__ import annotations

import inspect
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import UnivariateSpline

from .figure_s02_parse_mix_data import (
    CAPTURE_PRIMER_ORDER,
    CELL_TYPE_COLORS,
    CELL_TYPE_ORDER,
    MANUSCRIPT_TRANSCRIPT_ORDER,
    PRIMER_COLORS,
    PRIMER_ORDER,
)


def _violin_kwargs() -> dict[str, Any]:
    """Return width/bandwidth arguments supported by the installed seaborn.

    The source notebook used ``scale='width'`` and ``bw_adjust=1``. Seaborn
    renamed ``scale`` to ``density_norm`` in 0.13, so inspect the installed
    function instead of assuming one API and failing at runtime.
    """
    parameters = inspect.signature(sns.violinplot).parameters
    kwargs: dict[str, Any] = {
        "cut": 0,
        "inner": "box",
        "linewidth": 0.8,
    }
    if "density_norm" in parameters:
        kwargs["density_norm"] = "width"
    elif "scale" in parameters:
        kwargs["scale"] = "width"
    if "bw_adjust" in parameters:
        kwargs["bw_adjust"] = 1.0
    return kwargs


def _draw_violin(
    *,
    data: pd.DataFrame,
    category_axis: str,
    value_column: str,
    order,
    palette,
    ax,
):
    """Draw a categorical violin across seaborn 0.12 and 0.13 APIs."""
    if category_axis not in {"x", "y"}:
        raise ValueError(
            f"category_axis must be 'x' or 'y'; observed {category_axis!r}."
        )
    parameters = inspect.signature(sns.violinplot).parameters
    kwargs: dict[str, Any] = {
        "data": data,
        "order": list(order),
        "palette": list(palette),
        "ax": ax,
        **_violin_kwargs(),
    }
    if category_axis == "x":
        kwargs["x"] = "sample"
        kwargs["y"] = value_column
    else:
        kwargs["y"] = "sample"
        kwargs["x"] = value_column

    # Seaborn 0.13 deprecates palette-without-hue. Older releases do not
    # support the legend argument, so use redundant hue only when supported.
    if "legend" in parameters:
        kwargs.update(
            {
                "hue": "sample",
                "hue_order": list(order),
                "legend": False,
            }
        )
    return sns.violinplot(**kwargs)


def _format_median(value: float, *, decimals: int = 0) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.{decimals}f}"


def _annotate_medians_above(
    ax,
    data: pd.DataFrame,
    category: str,
    value: str,
    order,
    *,
    decimals=0,
):
    medians = data.groupby(category, observed=False)[value].median()
    for position, group in enumerate(order):
        median = float(medians.get(group, np.nan))
        ax.text(
            position,
            1.02,
            _format_median(median, decimals=decimals),
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=8,
        )


def _cell_type_violin(
    cell_metrics: pd.DataFrame,
    *,
    value: str,
    primers,
    ylabel: str,
    log_plus_one: bool,
    decimals: int,
):
    data = cell_metrics.loc[
        cell_metrics["cell_type"].isin(CELL_TYPE_ORDER)
        & cell_metrics["sample"].astype(str).isin(primers)
    ].copy()
    data["sample"] = pd.Categorical(
        data["sample"].astype(str), primers, ordered=True
    )
    plot_value = value
    if log_plus_one:
        plot_value = f"{value}_plot"
        data[plot_value] = pd.to_numeric(data[value], errors="coerce") + 1.0
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4), sharey=True)
    for ax, cell_type in zip(axes, CELL_TYPE_ORDER):
        subset = data.loc[data["cell_type"].eq(cell_type)].dropna(subset=[value])
        if subset.empty:
            raise ValueError(
                f"No {cell_type} cells are available for violin metric {value!r}."
            )
        _draw_violin(
            data=subset,
            category_axis="x",
            value_column=plot_value,
            order=primers,
            palette=[PRIMER_COLORS[primer] for primer in primers],
            ax=ax,
        )
        if log_plus_one:
            ax.set_yscale("log")
            ax.set_ylim(bottom=1)
        ax.set_title(cell_type, fontsize=10, fontweight="bold")
        ax.set_xlabel("Round 1 guide primer")
        ax.set_ylabel(ylabel if ax is axes[0] else "")
        ax.tick_params(axis="x", rotation=35)
        _annotate_medians_above(
            ax,
            subset,
            "sample",
            value,
            primers,
            decimals=decimals,
        )
        sns.despine(ax=ax)
    fig.tight_layout()
    return fig, data


def plot_panel_i(cell_metrics: pd.DataFrame):
    return _cell_type_violin(
        cell_metrics,
        value="matched_pool_guide_umis",
        primers=CAPTURE_PRIMER_ORDER,
        ylabel="Matched-pool guide UMIs per cell (+1, log scale)",
        log_plus_one=True,
        decimals=0,
    )


def plot_panel_j(cell_metrics: pd.DataFrame):
    fig, data = _cell_type_violin(
        cell_metrics,
        value="guide_pool_specificity",
        primers=CAPTURE_PRIMER_ORDER,
        ylabel="Guide-pool specificity",
        log_plus_one=False,
        decimals=2,
    )
    for ax in fig.axes[:2]:
        ax.set_ylim(0, 1.05)
    return fig, data


def plot_panel_l(cell_metrics: pd.DataFrame):
    return _cell_type_violin(
        cell_metrics,
        value="top_guide_umis",
        primers=PRIMER_ORDER,
        ylabel="Top-guide UMIs per cell (+1, log scale)",
        log_plus_one=True,
        decimals=0,
    )


def plot_panel_m(cell_metrics: pd.DataFrame, saturation: pd.DataFrame):
    transcript = cell_metrics.loc[
        cell_metrics["sample"].astype(str).isin(MANUSCRIPT_TRANSCRIPT_ORDER),
        ["sample", "total_counts"],
    ].copy()
    transcript["sample"] = pd.Categorical(
        transcript["sample"].astype(str),
        MANUSCRIPT_TRANSCRIPT_ORDER,
        ordered=True,
    )
    if transcript.empty:
        raise ValueError("No transcript-count cells are available for panel S2m.")
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(9.4, 3.8),
        gridspec_kw={"width_ratios": [1.0, 1.15]},
    )
    counts = transcript.groupby("sample", observed=False).size()
    _draw_violin(
        data=transcript,
        category_axis="y",
        value_column="total_counts",
        order=MANUSCRIPT_TRANSCRIPT_ORDER,
        palette=[PRIMER_COLORS[p] for p in MANUSCRIPT_TRANSCRIPT_ORDER],
        ax=axes[0],
    )
    axes[0].set_yticks(
        np.arange(len(MANUSCRIPT_TRANSCRIPT_ORDER)),
        [
            f"{primer} (n={int(counts.get(primer, 0)):,})"
            for primer in MANUSCRIPT_TRANSCRIPT_ORDER
        ],
    )
    axes[0].set_xlabel("Total transcriptome UMIs per cell")
    axes[0].set_ylabel("")
    sns.despine(ax=axes[0])

    marker_map = {"CS1": "o", "DC1": "s", "Dual": "^", "no_guide": "D"}
    for primer in MANUSCRIPT_TRANSCRIPT_ORDER:
        frame = saturation.loc[
            saturation["primer"].astype(str).eq(primer)
        ].sort_values("mean_assigned_reads_per_cell")
        if frame.empty:
            raise ValueError(f"No saturation rows were found for primer {primer!r}.")
        x = frame["mean_assigned_reads_per_cell"].to_numpy(float)
        y = frame["median_transcript_umis_per_cell"].to_numpy(float)
        x_unique, first = np.unique(x, return_index=True)
        y_unique = y[first]
        if len(x_unique) >= 4:
            smoothing = len(x_unique) * (np.std(y_unique) ** 2) * 0.001
            spline = UnivariateSpline(x_unique, y_unique, s=smoothing)
            x_plot = np.linspace(x_unique.min(), x_unique.max(), 400)
            y_plot = spline(x_plot)
            markevery = 30
        else:
            x_plot, y_plot, markevery = x_unique, y_unique, 1
        axes[1].plot(
            x_plot,
            y_plot,
            marker=marker_map[primer],
            markevery=markevery,
            linewidth=1.6,
            markersize=4,
            label=primer,
            color=PRIMER_COLORS[primer],
        )
    axes[1].set_xlabel("Mean assigned reads per cell")
    axes[1].set_ylabel("Median transcriptome UMIs per cell")
    axes[1].legend(frameon=False, fontsize=8)
    sns.despine(ax=axes[1])
    fig.tight_layout()
    return fig, transcript, saturation.copy()


def plot_panel_n(cell_metrics: pd.DataFrame):
    source = cell_metrics[
        ["sample", "cell_type", "hiPSC_guide_umis_ge3", "K562_guide_umis_ge3"]
    ].copy()
    keep = (
        source["hiPSC_guide_umis_ge3"].gt(0)
        | source["K562_guide_umis_ge3"].gt(0)
    )
    source = source.loc[keep].copy()
    if source.empty:
        raise ValueError(
            "No cells retain a >=3-UMI guide observation for panel S2n."
        )
    x = source["hiPSC_guide_umis_ge3"].to_numpy(float) + 1.0
    y = source["K562_guide_umis_ge3"].to_numpy(float) + 1.0
    high = float(max(np.max(x), np.max(y)))
    high_decade = max(10.0, 10 ** np.ceil(np.log10(high)))
    bins = np.logspace(0, np.log10(high_decade), 60)

    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    divider = make_axes_locatable(ax)
    ax_top = divider.append_axes("top", size="20%", pad=0.22, sharex=ax)
    ax_right = divider.append_axes("right", size="20%", pad=0.22, sharey=ax)
    histogram = ax.hist2d(
        x,
        y,
        bins=[bins, bins],
        cmap="viridis",
        norm=LogNorm(vmin=1),
        rasterized=True,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, high_decade)
    ax.set_ylim(1, high_decade)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(
        True,
        which="both",
        linestyle="--",
        linewidth=0.45,
        alpha=0.35,
    )
    ax.set_xlabel("Total hiPSC-pool guide UMIs (+1, log)")
    ax.set_ylabel("Total K562-pool guide UMIs (+1, log)")

    ax_top.set_xscale("log")
    ax_top.hist(
        x,
        bins=bins,
        histtype="stepfilled",
        alpha=0.6,
        color=CELL_TYPE_COLORS["hiPSC"],
    )
    ax_top.set_ylabel("Cells")
    ax_top.tick_params(axis="x", labelbottom=False)
    ax_top.grid(False)
    ax_right.set_yscale("log")
    ax_right.hist(
        y,
        bins=bins,
        orientation="horizontal",
        histtype="stepfilled",
        alpha=0.6,
        color=CELL_TYPE_COLORS["K562"],
    )
    ax_right.set_xlabel("Cells")
    ax_right.tick_params(axis="y", labelleft=False)
    ax_right.grid(False)

    color_axis = inset_axes(
        ax,
        width="4%",
        height="34%",
        loc="lower left",
        bbox_to_anchor=(1.30, 0.64, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    colorbar = fig.colorbar(histogram[3], cax=color_axis)
    colorbar.set_label("Cells per bin")
    fig.suptitle(
        f"Guide-pool mixing (n={len(source):,} cells)",
        y=0.96,
        fontsize=11,
    )
    return fig, source


__all__ = [
    "plot_panel_i",
    "plot_panel_j",
    "plot_panel_l",
    "plot_panel_m",
    "plot_panel_n",
]
