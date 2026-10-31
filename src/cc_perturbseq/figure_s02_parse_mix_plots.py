"""S2 panel adapters around the ported original plotnine functions.
Only column naming, figure dimensions, and return values are adapted here.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from plotnine import theme
from . import figure_s02_source_plotnine as original


def _source_columns(metrics):
    out = metrics.copy()
    out["cell_id"] = out.index.astype(str)
    names = {"total_guide_umis":"total_gRNA", "K562_guide_umis":"K562_total",
        "hiPSC_guide_umis":"iPSC_total", "K562_guide_umis_ge3":"K562_total_ge3",
        "hiPSC_guide_umis_ge3":"iPSC_total_ge3", "matched_pool_guide_umis":"expected_species_total_gRNA",
        "guide_pool_specificity":"correct_species_guide_fraction"}
    for current, old in names.items():
        if current in out:
            out[old] = out[current]
    return out


def _draw(plot,width,height):
    # Equivalent to the source's explicit width/height export request.
    return (plot + theme(figure_size=(width,height))).draw()


def plot_panel_i(cell_metrics):
    expected, _, source = original.source_specificity_panels(_source_columns(cell_metrics))
    return _draw(expected,7.4,4.2), cell_metrics.loc[source.index].copy()


def plot_panel_j(cell_metrics):
    _, specificity, source = original.source_specificity_panels(_source_columns(cell_metrics))
    return _draw(specificity,7.4,4.2), cell_metrics.loc[source.index].copy()


def plot_panel_l(cell_metrics):
    data = _source_columns(cell_metrics)
    # Source function explicitly supports the no-guide control; retain S2l's scope.
    a = original.make_top_guide_panel(data.loc[data.cell_type.eq("K562")],"Guide capture in K562 cells",include_no_guide=True)
    b = original.make_top_guide_panel(data.loc[data.cell_type.eq("hiPSC")],"Guide capture in hiPSC cells",include_no_guide=True)
    # The original theme hides minor grids. Explicitly suppress their automatic
    # calculation: mizani otherwise fails when the observed log range has only
    # one major tick. Major ticks, limits, values and visible source style stay intact.
    for plot in (a, b):
        plot.scales.get_scales("y").minor_breaks = []
    keep = cell_metrics.cell_type.isin(["K562","hiPSC"]) & cell_metrics["sample"].astype(str).isin(original.SAMPLE_ORDER)
    return _draw(a | b,9.4,4.6), cell_metrics.loc[keep].copy()


def plot_panel_m(cell_metrics,saturation):
    sat = saturation.rename(columns={"primer":"sample", "mean_assigned_reads_per_cell":"read_count", "median_transcript_umis_per_cell":"transcript_count"})
    plot, _ = original.source_transcriptome_panel(_source_columns(cell_metrics),sat)
    data = cell_metrics.loc[cell_metrics["sample"].astype(str).isin(original.SAMPLE_ORDER),["sample","total_counts"]].copy()
    return _draw(plot,11,4.2), data, saturation.copy()


def plot_panel_n(cell_metrics):
    plot, _, _, selected = original.source_mixing_panels(_source_columns(cell_metrics))
    cols = ["sample","cell_type","hiPSC_guide_umis_ge3","K562_guide_umis_ge3"]
    return _draw(plot,6.2,5.4), cell_metrics.loc[selected.index,cols].copy()


def plot_panel_n_marginals(cell_metrics):
    _, x, y, _ = original.source_mixing_panels(_source_columns(cell_metrics))
    return _draw(x,5.2,3.2), _draw(y,3.8,4.8)


__all__ = ["plot_panel_i","plot_panel_j","plot_panel_l","plot_panel_m","plot_panel_n","plot_panel_n_marginals"]
