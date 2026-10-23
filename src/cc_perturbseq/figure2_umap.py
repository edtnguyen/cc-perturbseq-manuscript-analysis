from __future__ import annotations

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    coord_fixed,
    element_rect,
    element_text,
    facet_wrap,
    geom_point,
    ggplot,
    scale_color_cmap,
    theme,
    theme_void,
)

from .single_cell_panels import expression_vector, program_vector


def plot_figure2e_umaps(
    adata,
    usage: pd.DataFrame,
    *,
    point_size: float = 0.42,
    background_size: float = 0.18,
    vmax_quantile: float = 0.99,
):
    """Render P17/P27 and SLC8A1/GATA4 UMAPs with plotnine and YlGnBu."""
    if "X_umap" not in adata.obsm:
        raise KeyError("AnnData lacks obsm['X_umap']")
    coords = np.asarray(adata.obsm["X_umap"])
    specs = [
        ("program", 17, "Programs\nP17\nLat. mesoderm, YAP/TAZ-induced"),
        ("gene", "SLC8A1", "Markers\nSLC8A1"),
        ("program", 27, "Programs\nP27\nMesenchyme, ECM"),
        ("gene", "GATA4", "Markers\nGATA4"),
    ]
    feature_order = [label for _, _, label in specs]
    background_frames: list[pd.DataFrame] = []
    foreground_frames: list[pd.DataFrame] = []
    source: dict[str, pd.DataFrame] = {}

    for kind, feature, label in specs:
        values = (
            program_vector(usage, feature)
            if kind == "program"
            else expression_vector(adata, str(feature))
        )
        values = np.asarray(values, dtype=float)
        finite = values[np.isfinite(values)]
        vmax = float(np.quantile(finite, vmax_quantile)) if finite.size else 1.0
        vmax = max(vmax, np.nextafter(0, 1))
        display = np.clip(np.nan_to_num(values, nan=0.0) / vmax, 0, 1)
        key = f"{kind}_{feature}"
        source[key] = pd.DataFrame(
            {
                "umap_1": coords[:, 0],
                "umap_2": coords[:, 1],
                "value": values,
                "display_vmax_p99": vmax,
                "display_fraction": display,
            },
            index=adata.obs_names,
        )
        background_frames.append(
            pd.DataFrame({"umap_1": coords[:, 0], "umap_2": coords[:, 1], "feature": label})
        )
        positive = np.flatnonzero(np.isfinite(values) & (values > 0))
        positive = positive[np.argsort(values[positive])]
        foreground_frames.append(
            pd.DataFrame(
                {
                    "umap_1": coords[positive, 0],
                    "umap_2": coords[positive, 1],
                    "display_fraction": display[positive],
                    "feature": label,
                }
            )
        )

    background = pd.concat(background_frames, ignore_index=True)
    foreground = pd.concat(foreground_frames, ignore_index=True)
    background["feature"] = pd.Categorical(background["feature"], feature_order, ordered=True)
    foreground["feature"] = pd.Categorical(foreground["feature"], feature_order, ordered=True)

    plot = (
        ggplot()
        + geom_point(
            background, aes("umap_1", "umap_2"), inherit_aes=False,
            color="#D9D9D9", size=background_size, stroke=0, raster=True,
        )
        + geom_point(
            foreground, aes("umap_1", "umap_2", color="display_fraction"),
            inherit_aes=False, size=point_size, stroke=0, raster=True,
        )
        + scale_color_cmap(cmap_name="YlGnBu", limits=(0, 1), guide=None)
        + facet_wrap("~feature", ncol=2)
        + coord_fixed()
        + theme_void()
        + theme(
            figure_size=(5.7, 5.2),
            strip_text=element_text(family="Liberation Sans", size=8.2, color="black"),
            strip_background=element_rect(fill="white", color="white"),
            panel_spacing=0.10,
            plot_background=element_rect(fill="white", color="white"),
            panel_background=element_rect(fill="white", color="white"),
            plot_margin=0.02,
        )
    )
    return plot, source
