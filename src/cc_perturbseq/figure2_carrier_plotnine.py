from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matplotlib as mpl
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    coord_equal,
    element_blank,
    element_rect,
    element_text,
    facet_wrap,
    geom_point,
    ggplot,
    labs,
    scale_color_gradient2,
    theme,
    theme_bw,
)
from scipy import sparse
from scipy.stats import fisher_exact
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests

from .figure2_ec import _assignment_matrix, _typed_guide_map
from .program_io import clean_day
from .single_cell_panels import _umap


FIGURE2F_TARGETS = ["ETV2", "KDR", "PRKCD", "PLCG1"]

# Recovered settings from:
# programDE/notebooks/fig3_carrierseeded_neighborhood_plotnine_style_notebook_package/
# fig3_carrierseeded_neighborhood_panels_interactive_plotnine_style.ipynb
FIGURE2F_K = 300
FIGURE2F_PSEUDOCOUNT = 0.5
FIGURE2F_MAX_BACKGROUND_CELLS = 60_000
FIGURE2F_PLOT_SEED = 17
FIGURE2F_LOG2OR_CLIP = 2.0
FIGURE2F_CMAP_LOW = "#053061"
FIGURE2F_CMAP_MID = "#F7F7F7"
FIGURE2F_CMAP_HIGH = "#67001F"
FIGURE2F_BACKGROUND_SIZE = 0.005
FIGURE2F_FOREGROUND_SIZE = 0.02


def _zscore_usage(usage: pd.DataFrame) -> np.ndarray:
    values = usage.to_numpy(float)
    mean = np.nanmean(values, axis=0)
    std = np.nanstd(values, axis=0)
    std[~np.isfinite(std) | (std == 0)] = 1.0
    return np.nan_to_num((values - mean) / std)


def _any_guides(assignment, columns: np.ndarray) -> np.ndarray:
    if not len(columns):
        return np.zeros(assignment.shape[0], dtype=bool)
    if sparse.issparse(assignment):
        return np.asarray(assignment[:, columns].sum(axis=1)).ravel() > 0
    return np.asarray(assignment[:, columns]).any(axis=1)


def _bh(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    result = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p)
    if finite.any():
        result[finite] = multipletests(p[finite], method="fdr_bh")[1]
    return result


def compute_figure2f_carrier_seeded(
    adata,
    usage: pd.DataFrame,
    guide_metadata_path: str | Path,
    *,
    targets: Sequence[str] = FIGURE2F_TARGETS,
    k: int = FIGURE2F_K,
    pseudocount: float = FIGURE2F_PSEUDOCOUNT,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Repeat the recovered carrier-seeded neighborhood analysis.

    This function ports the analysis choices from
    `fig3_carrierseeded_neighborhood_panels_interactive_plotnine_style.ipynb`:

    * all cells carrying a target guide are seed cells;
    * pooled control carriers comprise all 120 non-targeting/safe-targeting
      guides and exclude cells carrying the target under evaluation;
    * neighbors are restricted to the seed cell's differentiation day;
    * Euclidean KNN uses z-scored activity across all 50 cNMF programs;
    * local enrichment is the target/control odds in the KNN divided by the
      same-day background target/control odds, with a 0.5 pseudocount;
    * Fisher exact P values are BH-adjusted across seed cells for each target,
      with an additional within-day FDR column.
    """
    if pseudocount <= 0:
        raise ValueError("pseudocount must be positive")

    assignment = _assignment_matrix(adata)
    guide_map = _typed_guide_map(adata, guide_metadata_path)
    control_columns = guide_map.index[
        guide_map["is_non_targeting"] | guide_map["is_safe_targeting"]
    ].to_numpy(int)
    if len(control_columns) != 120:
        raise ValueError(
            "Expected 120 pooled negative-control guides "
            "(60 non-targeting plus 60 safe-targeting); "
            f"found {len(control_columns)}."
        )

    control_raw = _any_guides(assignment, control_columns)
    days = (
        pd.Series(adata.obs["day"], index=adata.obs_names)
        .map(clean_day)
        .astype(str)
        .to_numpy()
    )
    if "cell_type" in adata.obs.columns:
        cell_types = (
            pd.Series(adata.obs["cell_type"], index=adata.obs_names)
            .astype(str)
            .to_numpy()
        )
    else:
        cell_types = np.repeat("unknown", adata.n_obs)

    coords = _umap(adata)
    embedding = _zscore_usage(usage)

    n_total_guides = np.asarray(assignment.sum(axis=1)).ravel().astype(int)
    targeting_columns = guide_map.index[guide_map["is_targeting"]].to_numpy(int)
    targeting_matrix = assignment[:, targeting_columns]
    n_targeting_guides = (
        np.asarray(targeting_matrix.sum(axis=1)).ravel().astype(int)
    )

    all_results: list[pd.DataFrame] = []
    for target in targets:
        target_columns = guide_map.index[
            guide_map["is_targeting"] & guide_map["target"].eq(target)
        ].to_numpy(int)
        if not len(target_columns):
            raise KeyError(f"No guide columns found for panel-f target {target}")

        target_raw = _any_guides(assignment, target_columns)
        target_mask = target_raw
        control_mask = control_raw & (~target_raw)
        seed_indices = np.flatnonzero(target_mask)
        rows: list[dict[str, Any]] = []

        for day in sorted(pd.unique(days[seed_indices])):
            day_cells = np.flatnonzero(days == day)
            day_seeds = seed_indices[days[seed_indices] == day]
            if len(day_cells) <= 1 or not len(day_seeds):
                continue

            k_query = min(k + 1, len(day_cells))
            model = NearestNeighbors(
                n_neighbors=k_query,
                metric="euclidean",
            )
            model.fit(embedding[day_cells])
            neighbor_local = model.kneighbors(
                embedding[day_seeds],
                return_distance=False,
            )

            day_mask = days == day
            background_target_total_all = int(np.sum(target_mask & day_mask))
            background_control_total = int(np.sum(control_mask & day_mask))

            for seed_index, local_neighbors in zip(day_seeds, neighbor_local):
                neighbors = day_cells[local_neighbors]
                neighbors = neighbors[neighbors != seed_index][:k]

                local_target_n = int(np.sum(target_mask[neighbors]))
                local_control_n = int(np.sum(control_mask[neighbors]))
                background_target_total = max(
                    background_target_total_all - 1,
                    0,
                )

                local_ratio = (
                    local_target_n + pseudocount
                ) / (
                    local_control_n + pseudocount
                )
                background_ratio = (
                    background_target_total + pseudocount
                ) / (
                    background_control_total + pseudocount
                )
                log2_enrichment = float(
                    np.log2(local_ratio / background_ratio)
                )

                target_outside = max(
                    background_target_total - local_target_n,
                    0,
                )
                control_outside = max(
                    background_control_total - local_control_n,
                    0,
                )
                fisher_p = float(
                    fisher_exact(
                        [
                            [local_target_n, local_control_n],
                            [target_outside, control_outside],
                        ],
                        alternative="two-sided",
                    ).pvalue
                )
                local_total = local_target_n + local_control_n
                background_total = (
                    background_target_total + background_control_total
                )
                rows.append(
                    {
                        "target": target,
                        "seed_cell_idx": int(seed_index),
                        "seed_cell_id": str(adata.obs_names[seed_index]),
                        "day": day,
                        "cell_type": cell_types[seed_index],
                        "mode": "all_carriers",
                        "k_requested": int(k),
                        "k_used": int(len(neighbors)),
                        "restrict_neighbors_to_same_day": True,
                        "same_day_baseline": True,
                        "local_target_n": local_target_n,
                        "local_nt_n": local_control_n,
                        "local_target_plus_nt_n": local_total,
                        "background_target_total": background_target_total,
                        "background_nt_total": background_control_total,
                        "local_target_fraction_among_target_nt": (
                            local_target_n / local_total
                            if local_total > 0
                            else np.nan
                        ),
                        "background_target_fraction_among_target_nt": (
                            background_target_total / background_total
                            if background_total > 0
                            else np.nan
                        ),
                        "log2_local_enrichment": log2_enrichment,
                        "fisher_p": fisher_p,
                        "UMAP1": float(coords[seed_index, 0]),
                        "UMAP2": float(coords[seed_index, 1]),
                        "n_total_guides": int(n_total_guides[seed_index]),
                        "n_targeting_targets": int(
                            n_targeting_guides[seed_index]
                        ),
                    }
                )

        target_scores = pd.DataFrame(rows)
        if target_scores.empty:
            continue

        target_scores["fisher_fdr"] = _bh(target_scores["fisher_p"])
        target_scores["fisher_fdr_by_day"] = np.nan
        for _, subset in target_scores.groupby("day"):
            target_scores.loc[
                subset.index,
                "fisher_fdr_by_day",
            ] = _bh(subset["fisher_p"])
        target_scores["minus_log10_fdr"] = -np.log10(
            target_scores["fisher_fdr"].clip(lower=1e-300)
        )
        target_scores["minus_log10_fdr_by_day"] = -np.log10(
            target_scores["fisher_fdr_by_day"].clip(lower=1e-300)
        )
        target_scores["significant_fdr_0_05"] = (
            target_scores["fisher_fdr"] < 0.05
        )
        target_scores["significant_fdr_by_day_0_05"] = (
            target_scores["fisher_fdr_by_day"] < 0.05
        )
        all_results.append(target_scores)

    if not all_results:
        raise RuntimeError("No carrier-seeded neighborhood scores were produced.")

    scores = pd.concat(all_results, ignore_index=True)
    background = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "UMAP1": coords[:, 0],
            "UMAP2": coords[:, 1],
            "day_clean": days,
            "cell_type_clean": cell_types,
            "n_total_guides": n_total_guides,
            "n_targeting_targets": n_targeting_guides,
        }
    )
    return scores, background


def _geom_point_optional_raster(
    data,
    mapping,
    *,
    raster: bool = False,
    **kwargs,
):
    """Use raster=True for background points when supported by plotnine."""
    if raster:
        try:
            return geom_point(data, mapping, raster=True, **kwargs)
        except TypeError:
            return geom_point(data, mapping, **kwargs)
    return geom_point(data, mapping, **kwargs)


def plot_figure2f(
    scores: pd.DataFrame,
    background: pd.DataFrame,
    *,
    targets: Sequence[str] = FIGURE2F_TARGETS,
    max_background_cells: int = FIGURE2F_MAX_BACKGROUND_CELLS,
    seed: int = FIGURE2F_PLOT_SEED,
    log2or_clip: float = FIGURE2F_LOG2OR_CLIP,
):
    """Repeat the recovered faceted plotnine UMAP renderer."""
    mpl.rcParams.update(
        {
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "text.usetex": False,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
            "savefig.transparent": False,
        }
    )

    rng = np.random.default_rng(seed)
    background = background.copy()
    if len(background) > max_background_cells:
        parts: list[pd.DataFrame] = []
        n_days = max(1, background["day_clean"].nunique())
        n_per_day = max_background_cells // n_days
        for _, subset in background.groupby("day_clean", sort=True):
            if len(subset) > n_per_day:
                indices = rng.choice(
                    subset.index.to_numpy(),
                    size=n_per_day,
                    replace=False,
                )
                parts.append(background.loc[indices])
            else:
                parts.append(subset)
        background_plot = pd.concat(parts, ignore_index=True)
    else:
        background_plot = background.copy()

    scores_plot = scores.loc[
        scores["target"].astype(str).isin(
            [str(target) for target in targets]
        )
    ].copy()
    scores_plot["target"] = pd.Categorical(
        scores_plot["target"].astype(str),
        categories=list(targets),
        ordered=True,
    )
    scores_plot["log2OR_plot"] = scores_plot[
        "log2_local_enrichment"
    ].clip(-log2or_clip, log2or_clip)
    scores_plot["carrier_fraction_pct"] = (
        100.0
        * scores_plot["local_target_n"]
        / scores_plot["k_used"].replace(0, np.nan)
    )

    xmin = min(background["UMAP1"].min(), scores_plot["UMAP1"].min())
    xmax = max(background["UMAP1"].max(), scores_plot["UMAP1"].max())
    ymin = min(background["UMAP2"].min(), scores_plot["UMAP2"].min())
    ymax = max(background["UMAP2"].max(), scores_plot["UMAP2"].max())
    pad_x = (xmax - xmin) * 0.03
    pad_y = (ymax - ymin) * 0.03

    plot = (
        ggplot()
        + _geom_point_optional_raster(
            background_plot,
            aes("UMAP1", "UMAP2"),
            raster=True,
            color="#D0D0D0",
            size=FIGURE2F_BACKGROUND_SIZE,
            alpha=0.12,
            show_legend=False,
        )
        + geom_point(
            scores_plot,
            aes("UMAP1", "UMAP2", color="log2OR_plot"),
            size=FIGURE2F_FOREGROUND_SIZE,
            alpha=0.95,
        )
        + facet_wrap("~target", ncol=2)
        + scale_color_gradient2(
            low=FIGURE2F_CMAP_LOW,
            mid=FIGURE2F_CMAP_MID,
            high=FIGURE2F_CMAP_HIGH,
            midpoint=0,
            limits=(-log2or_clip, log2or_clip),
            breaks=[-2, -1, 0, 1, 2],
            name="Carrier-seeded\nenrichment\nlog2OR",
        )
        + coord_equal(
            xlim=(xmin - pad_x, xmax + pad_x),
            ylim=(ymin - pad_y, ymax + pad_y),
            expand=False,
        )
        + labs(
            title="CRISPRi effects",
            subtitle=(
                "Carrier-seeded neighborhoods; all D0–D3 cells shown "
                "in one shared embedding"
            ),
            x="single-cell state embedding 1",
            y="single-cell state embedding 2",
        )
        + theme_bw(base_size=12, base_family="DejaVu Sans")
        + theme(
            figure_size=(8.2, 7.0),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(
                color="#7f7f7f",
                fill=None,
                size=0.8,
            ),
            axis_title=element_text(size=13, family="DejaVu Sans"),
            axis_text=element_text(size=10, family="DejaVu Sans"),
            plot_title=element_text(
                size=17,
                weight="bold",
                family="DejaVu Sans",
            ),
            plot_subtitle=element_text(
                size=10.5,
                family="DejaVu Sans",
            ),
            strip_background=element_rect(
                fill="#D9D9D9",
                color="#7f7f7f",
                size=0.8,
            ),
            strip_text=element_text(
                size=14,
                weight="bold",
                family="DejaVu Sans",
            ),
            legend_title=element_text(
                size=11,
                weight="bold",
                family="DejaVu Sans",
            ),
            legend_text=element_text(
                size=10,
                family="DejaVu Sans",
            ),
        )
    )
    return plot, scores_plot, background_plot
