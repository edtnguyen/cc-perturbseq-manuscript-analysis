"""S2 plotnine source port. No source notebooks are read at runtime.
Paper panels: original cells 39/41; yield/saturation: 21; mixing: 27.
Original source: K562_hiPSC_ParseV3_plotnine_figure_reproduction.ipynb.
"""
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Sequence
import json
import os
import re
import sys
import warnings
import zipfile
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import mannwhitneyu
import matplotlib as mpl
from plotnine import (
    ggplot, aes, annotate,
    geom_abline, geom_boxplot, geom_histogram, geom_jitter,
    geom_line, geom_point, geom_text, geom_violin,
    stat_bin_2d,
    facet_grid, facet_wrap,
    coord_equal, coord_flip,
    labs,
    scale_color_cmap, scale_color_manual,
    scale_fill_cmap, scale_fill_manual,
    scale_x_continuous, scale_x_log10,
    scale_y_continuous, scale_y_log10,
    position_dodge, position_jitterdodge,
    theme, theme_bw,
    element_blank, element_line, element_rect, element_text,
)
from mizani.formatters import comma_format, percent_format
from IPython.display import display

SAMPLE_ORDER = ["CS1", "DC1", "Dual", "no_guide"]
PRIMER_SAMPLES = ["CS1", "DC1", "Dual"]
SAMPLE_COLORS = {"CS1":"#4E79A7", "DC1":"#F28E2B", "Dual":"#E15759", "no_guide":"#76B7B2"}
CELLTYPE_COLORS = {"K562":"#4E79A7", "hiPSC":"#F28E2B"}
GUIDE_UMI_THRESHOLD = 3
mpl.rcParams.update({"svg.fonttype":"none", "pdf.fonttype":42, "ps.fonttype":42})
def normalize_cell_id(value: object) -> str:
    """Remove a split-pipe sublibrary suffix such as '__s1'."""
    return re.sub(r"__s\d+$", "", str(value))

def ordered_sample(series: pd.Series) -> pd.Categorical:
    return pd.Categorical(series.astype(str), categories=SAMPLE_ORDER, ordered=True)

def first_existing(paths: Sequence[Path]) -> Path | None:
    return next((p for p in paths if p.exists()), None)

def theme_publication(base_size: int = 10):
    return (
        theme_bw(base_size=base_size)
        + theme(
            panel_grid_minor=element_blank(),
            panel_grid_major=element_line(color="#D8D8D8", size=0.45),
            panel_border=element_rect(color="black", fill=None, size=0.7),
            axis_title=element_text(size=base_size + 1),
            axis_text=element_text(size=base_size),
            legend_title=element_text(size=base_size),
            legend_text=element_text(size=max(base_size - 1, 7)),
            strip_background=element_rect(fill="white", color="black", size=0.7),
            strip_text=element_text(size=base_size, weight="bold"),
            plot_title=element_text(size=base_size + 3, weight="bold"),
            plot_subtitle=element_text(size=base_size + 1),
        )
    )

def make_top_guide_panel(data: pd.DataFrame, title: str, include_no_guide: bool = False):
    samples = SAMPLE_ORDER if include_no_guide else PRIMER_SAMPLES
    x = data.loc[data["sample"].astype(str).isin(samples)].copy()
    x["sample"] = pd.Categorical(x["sample"].astype(str), categories=samples, ordered=True)
    x["top_plot"] = x["top_guide_umis"] + 1
    counts = x.groupby("sample", observed=True).size()
    x["sample_label"] = x["sample"].astype(str).map(
        {s: f"{s}\n(n={int(counts.get(s, 0)):,})" for s in samples}
    )
    x["sample_label"] = pd.Categorical(
        x["sample_label"],
        categories=[f"{s}\n(n={int(counts.get(s, 0)):,})" for s in samples],
        ordered=True,
    )
    ann = (
        x.groupby(["sample", "sample_label"], observed=True)
        .agg(median=("top_guide_umis", "median"), y=("top_plot", lambda z: z.quantile(0.995) * 1.25))
        .reset_index()
    )
    ann["label"] = ann["median"].round().astype(int).astype(str)
    return (
        ggplot(x, aes(x="sample_label", y="top_plot", fill="sample"))
        + geom_violin(trim=True, width=0.90, alpha=0.88)
        + geom_boxplot(width=0.16, fill="white", color="black", outlier_alpha=0)
        + geom_text(data=ann, mapping=aes(x="sample_label", y="y", label="label"),
                    inherit_aes=False, size=8)
        + scale_y_log10(labels=comma_format())
        + scale_fill_manual(values=SAMPLE_COLORS, guide=None)
        + labs(x=None, y="Top guide UMIs per cell", title=title, subtitle="Labels show medians")
        + theme_publication(10)
    )
def source_specificity_panels(cell_metrics):
    paper_type_col = "cell_type"
    # Expected-species guide specificity and expected-species UMIs.
    specificity = cell_metrics.loc[
        cell_metrics[paper_type_col].isin(["K562", "hiPSC"])
        & cell_metrics["sample"].astype(str).isin(PRIMER_SAMPLES)
    ].copy()
    specificity["plot_cell_type"] = specificity[paper_type_col]
    specificity["sample"] = pd.Categorical(specificity["sample"].astype(str), categories=PRIMER_SAMPLES, ordered=True)

    if "correct_species_guide_fraction" not in specificity.columns:
        expected = np.where(
            specificity["plot_cell_type"].eq("K562"), specificity["K562_total"], specificity["iPSC_total"]
        )
        off = np.where(
            specificity["plot_cell_type"].eq("K562"), specificity["iPSC_total"], specificity["K562_total"]
        )
        denom = expected + off
        specificity["correct_species_guide_fraction"] = np.divide(
            expected, denom, out=np.full(len(specificity), np.nan), where=denom > 0
        )
        specificity["expected_species_total_gRNA"] = expected

    p_specificity = (
        ggplot(specificity, aes(x="sample", y="correct_species_guide_fraction", fill="sample"))
        + geom_violin(trim=True, alpha=0.85)
        + geom_boxplot(width=0.15, fill="white", color="black", outlier_alpha=0)
        + facet_wrap("~plot_cell_type", nrow=1)
        + scale_fill_manual(values=SAMPLE_COLORS, guide=None)
        + scale_y_continuous(limits=(0, 1.02), labels=percent_format())
        + labs(x=None, y="Expected-species guide fraction", title="Expected-species guide specificity")
        + theme_publication(10)
    )

    expected_df = specificity.copy()
    expected_df["expected_plot"] = expected_df["expected_species_total_gRNA"] + 1
    p_expected = (
        ggplot(expected_df, aes(x="sample", y="expected_plot", fill="sample"))
        + geom_violin(trim=True, alpha=0.85)
        + geom_boxplot(width=0.15, fill="white", color="black", outlier_alpha=0)
        + facet_wrap("~plot_cell_type", nrow=1)
        + scale_fill_manual(values=SAMPLE_COLORS, guide=None)
        + scale_y_log10(labels=comma_format())
        + labs(x=None, y="Expected-species guide UMIs per cell", title="Expected-species guide UMIs")
        + theme_publication(10)
    )
    return p_expected, p_specificity, specificity

def source_transcriptome_panel(cell_metrics, saturation):
    yield_df = cell_metrics[["sample", "total_counts"]].dropna().copy()
    yield_df["sample"] = ordered_sample(yield_df["sample"])
    counts = yield_df.groupby("sample", observed=True).size()
    yield_df["sample_label"] = yield_df["sample"].astype(str).map(
        {s: f"{s} (n={int(counts.get(s, 0)):,})" for s in SAMPLE_ORDER}
    )
    yield_order = [f"{s} (n={int(counts.get(s, 0)):,})" for s in SAMPLE_ORDER]
    yield_df["sample_label"] = pd.Categorical(yield_df["sample_label"], categories=yield_order, ordered=True)

    p_yield = (
        ggplot(yield_df, aes(x="sample_label", y="total_counts", fill="sample"))
        + geom_violin(trim=True, width=0.9, alpha=0.75)
        + geom_boxplot(width=0.12, outlier_alpha=0, fill="white", color="black")
        + coord_flip()
        + scale_fill_manual(values=SAMPLE_COLORS, guide=None)
        + labs(x=None, y="Total UMI counts per cell", title="Transcriptome yield")
        + theme_publication(10)
        + theme(axis_title_y=element_blank())
    )

    if not saturation.empty:
        p_saturation = (
            ggplot(saturation, aes(x="read_count", y="transcript_count", color="sample", group="sample"))
            + geom_line(size=0.8)
            + geom_point(size=1.5)
            + scale_color_manual(values=SAMPLE_COLORS)
            + labs(
                x="Mean reads per cell",
                y="Median transcript count",
                color="Primer type",
                title="Sequencing saturation",
            )
            + theme_publication(10)
            + theme(legend_position=(0.18, 0.76))
        )
        p_transcriptome = p_yield | p_saturation
    else:
        p_transcriptome = p_yield
    return p_transcriptome, yield_df

def source_mixing_panels(cell_metrics):
    i_col = f"iPSC_total_ge{GUIDE_UMI_THRESHOLD}"
    k_col = f"K562_total_ge{GUIDE_UMI_THRESHOLD}"

    mix_df = cell_metrics[["cell_id", "sample", i_col, k_col]].copy()
    mix_df = mix_df.loc[(mix_df[i_col] >= GUIDE_UMI_THRESHOLD) | (mix_df[k_col] >= GUIDE_UMI_THRESHOLD)].copy()
    mix_df["iPSC_plot"] = mix_df[i_col] + 1
    mix_df["K562_plot"] = mix_df[k_col] + 1

    p_mix = (
        ggplot(mix_df, aes(x="iPSC_plot", y="K562_plot"))
        + stat_bin_2d(bins=60)
        + geom_abline(slope=1, intercept=0, linetype="dashed", color="white", alpha=0.7)
        + scale_x_log10(labels=comma_format())
        + scale_y_log10(labels=comma_format())
        + scale_fill_cmap(cmap_name="viridis", trans="log10", name="Cells per bin")
        + coord_equal()
        + labs(
            x="Total iPSC guide UMIs + 1 (log scale)",
            y="Total K562 guide UMIs + 1 (log scale)",
            title=f"iPSC versus K562 guide UMIs per cell (n={len(mix_df):,})",
            subtitle=f"Guide counts were retained when the per-guide UMI count was ≥ {GUIDE_UMI_THRESHOLD}",
        )
        + theme_publication(10)
    )

    p_mix_x = (
        ggplot(mix_df, aes(x="iPSC_plot"))
        + geom_histogram(bins=60, fill="#4E79A7", alpha=0.8)
        + scale_x_log10(labels=comma_format())
        + labs(x="iPSC guide UMIs + 1", y="Cells", title="iPSC-guide marginal distribution")
        + theme_publication(9)
    )

    p_mix_y = (
        ggplot(mix_df, aes(x="K562_plot"))
        + geom_histogram(bins=60, fill="#76B7B2", alpha=0.8)
        + scale_x_log10(labels=comma_format())
        + coord_flip()
        + labs(x="K562 guide UMIs + 1", y="Cells", title="K562-guide marginal distribution")
        + theme_publication(9)
    )
    return p_mix, p_mix_x, p_mix_y, mix_df
