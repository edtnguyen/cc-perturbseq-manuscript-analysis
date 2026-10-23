from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

from .single_cell_panels import expression_vector, program_vector

def plot_figure2e_umaps(
    adata,
    usage: pd.DataFrame,
    *,
    point_size: float = 1.0,
    background_size: float = 0.25,
    vmax_quantile: float = 0.99,
):
    if "X_umap" not in adata.obsm:
        raise KeyError("AnnData lacks obsm['X_umap']")
    coords = np.asarray(adata.obsm["X_umap"])
    features = [
        ("program", 17, "P17", "Lat. mesoderm, YAP/TAZ-induced"),
        ("gene", "SLC8A1", "SLC8A1", ""),
        ("program", 27, "P27", "Mesenchyme, ECM"),
        ("gene", "GATA4", "GATA4", ""),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(5.7, 5.2))
    fig.subplots_adjust(left=0.03, right=0.94, bottom=0.03, top=0.90, wspace=0.15, hspace=0.18)
    source: dict[str, pd.DataFrame] = {}
    for axis, (kind, feature, title, subtitle) in zip(axes.ravel(), features):
        values = program_vector(usage, feature) if kind == "program" else expression_vector(adata, str(feature))
        values = np.asarray(values, dtype=float)
        finite = values[np.isfinite(values)]
        vmax = float(np.quantile(finite, vmax_quantile)) if finite.size else 1.0
        vmax = max(vmax, np.nextafter(0, 1))
        axis.scatter(coords[:, 0], coords[:, 1], s=background_size, c="#D9D9D9", linewidths=0, rasterized=True)
        positive = np.flatnonzero(np.isfinite(values) & (values > 0))
        positive = positive[np.argsort(values[positive])]
        scatter = axis.scatter(
            coords[positive, 0], coords[positive, 1], c=np.clip(values[positive], 0, vmax),
            cmap="YlGnBu", norm=Normalize(0, vmax), s=point_size, linewidths=0, rasterized=True,
        )
        panel_title = title if not subtitle else f"{title}\n{subtitle}"
        axis.set_title(
            panel_title, fontsize=8.2,
            fontstyle="italic" if kind == "gene" else "normal",
            pad=3, linespacing=1.05,
        )
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_visible(False)
        cb = fig.colorbar(scatter, ax=axis, fraction=0.040, pad=0.01)
        cb.ax.tick_params(labelsize=5, length=1.5)
        source[f"{kind}_{feature}"] = pd.DataFrame({
            "umap_1": coords[:, 0], "umap_2": coords[:, 1], "value": values,
        }, index=adata.obs_names)
    fig.text(0.25, 0.985, "Programs", ha="center", va="top", fontsize=10, fontweight="bold")
    fig.text(0.74, 0.985, "Markers", ha="center", va="top", fontsize=10, fontweight="bold")
    fig.text(0.005, 0.985, "e", ha="left", va="top", fontsize=15, fontweight="bold")
    return fig, source
