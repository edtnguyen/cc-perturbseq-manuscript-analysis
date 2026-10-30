"""Figure 4 panel producers with source-faithful feature UMAPs.

The retained legacy module contains all non-UMAP calculations. The feature
UMAP producer is overridden by the exact later ``jesse_umap.ipynb``
Seurat-style implementation shared across Figures 2–4. Motif enrichment and
logo inputs are locked to the author-selected active files.
"""
from . import figure4_panels_legacy as _legacy
from .motif_sources import (
    canonical_enhancer_enrichment_path,
    canonical_motif_logo_inputs,
)

for _name in dir(_legacy):
    if not _name.startswith("__"):
        globals()[_name] = getattr(_legacy, _name)

from .manuscript_feature_umaps import plot_figure4d_feature_umaps


def make_figure4e_annotation(
    *,
    workbook_path,
    spectra_path,
    motif_path,
    cluster_map_path,
    modisco_h5_path,
    top_motifs_per_program=5,
):
    """Run the recovered Figure 4e workflow with canonical motif inputs."""
    canonical_map, canonical_h5 = canonical_motif_logo_inputs(
        configured_cluster_map=cluster_map_path,
        configured_modisco_h5=modisco_h5_path,
    )
    return _legacy.make_figure4e_annotation(
        workbook_path=workbook_path,
        spectra_path=spectra_path,
        motif_path=canonical_enhancer_enrichment_path(motif_path),
        cluster_map_path=canonical_map,
        modisco_h5_path=canonical_h5,
        top_motifs_per_program=top_motifs_per_program,
    )


del _name
