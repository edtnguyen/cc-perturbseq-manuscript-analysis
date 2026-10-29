"""Figure 4 panel producers with source-faithful feature UMAPs.

The retained legacy module contains all non-UMAP calculations. The feature
UMAP producer is overridden by the exact later `jesse_umap.ipynb`
Seurat-style implementation shared across Figures 2–4.
"""
from . import figure4_panels_legacy as _legacy

for _name in dir(_legacy):
    if not _name.startswith("__"):
        globals()[_name] = getattr(_legacy, _name)

from .manuscript_feature_umaps import plot_figure4d_feature_umaps

del _name
