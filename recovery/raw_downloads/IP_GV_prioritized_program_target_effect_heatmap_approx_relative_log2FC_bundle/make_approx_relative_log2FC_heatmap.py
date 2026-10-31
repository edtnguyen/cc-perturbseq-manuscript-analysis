#!/usr/bin/env python3
# Recreate the approximate-relative-log2FC heatmap from:
# /mnt/data/IP_GV_prioritized_program_target_effect_heatmap_approx_log2FC/IP_GV_prioritized_program_nominal_target_effect_heatmap_approx_relative_log2FC_data.csv

from pathlib import Path
import pandas as pd
import matplotlib as mpl
from plotnine import (
    aes, coord_fixed, element_blank, element_text, geom_text,
    geom_tile, geom_vline, ggplot, labs, scale_fill_gradient2,
    theme, theme_classic,
)

DATA = Path('/mnt/data/IP_GV_prioritized_program_target_effect_heatmap_approx_log2FC/IP_GV_prioritized_program_nominal_target_effect_heatmap_approx_relative_log2FC_data.csv')
OUT = DATA.parent
df = pd.read_csv(DATA)

program_order = ["P15\nHomeostatic EC", "P9\nAngiogenic migration", "P6\nInflammatory adhesion", "P12\nsGC-cGMP remodeling"]
target_order = ["ATP13A3", "ACVRL1", "EIF2AK4", "CAV1", "KDR", "SMAD4\u2020", "ENG", "SMAD9", "TRP53", "NOTCH1", "CDKN1A"]

df["target_label"] = pd.Categorical(
    df["target_label"], categories=target_order, ordered=True
)
df["program_label"] = pd.Categorical(
    df["program_label"],
    categories=list(reversed(program_order)),
    ordered=True,
)
df["star"] = df["star"].fillna("")

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"

separator_x = target_order.index("TRP53") + 0.5

p = (
    ggplot(
        df,
        aes(
            x="target_label",
            y="program_label",
            fill="approx_relative_log2FC_equal_mouse_mean",
        ),
    )
    + geom_tile(color="white", size=0.8)
    + geom_text(aes(label="star"), size=13, color="#111111")
    + geom_vline(
        xintercept=separator_x,
        color="#666666",
        linetype="dashed",
        size=0.55,
    )
    + scale_fill_gradient2(
        low="#2166AC",
        mid="#F7F7F7",
        high="#B2182B",
        midpoint=0,
        limits=(-2.1, 2.1),
        breaks=[-2, -1, 0, 1, 2],
        na_value="#E6E6E6",
        name="Approx. relative\nprogram log2FC",
    )
    + coord_fixed(ratio=1)
    + labs(
        title="Nominal cross-mouse effects on prioritized cNMF programs",
        subtitle=(
            "Color = approximate relative program log2FC, averaged equally "
            "across IP_GFP and GV_GFP"
        ),
        x="Perturbed target gene",
        y="",
        caption=(
            "Approx. log2FC = (30/29) x beta_CLR / ln(2); "
            "* = same-sign signed-Stouffer p < 0.05. "
            "Gray = absent from the filtered nominal-pair table."
        ),
    )
    + theme_classic()
    + theme(
        figure_size=(9.8, 5.9),
        plot_title=element_text(size=14, weight="bold", ha="center"),
        plot_subtitle=element_text(size=8.8, ha="center", color="#4D4D4D"),
        plot_caption=element_text(size=7.0, ha="left", color="#555555"),
        axis_title_x=element_text(size=9),
        axis_title_y=element_blank(),
        axis_text_x=element_text(
            rotation=45, ha="right", size=9.3, style="italic"
        ),
        axis_text_y=element_text(size=9.3, ha="right"),
        axis_ticks=element_blank(),
        axis_line=element_blank(),
        panel_grid=element_blank(),
        legend_position="bottom",
        legend_direction="horizontal",
        legend_title=element_text(size=8.2, ha="center"),
        legend_text=element_text(size=7.7),
        legend_key_width=140,
        legend_key_height=10,
    )
)

p.save(
    OUT / "IP_GV_prioritized_program_nominal_target_effect_heatmap_"
          "approx_relative_log2FC_square_cells_editable.pdf",
    verbose=False,
)
p.save(
    OUT / "IP_GV_prioritized_program_nominal_target_effect_heatmap_"
          "approx_relative_log2FC_square_cells_editable.svg",
    verbose=False,
)
p.save(
    OUT / "IP_GV_prioritized_program_nominal_target_effect_heatmap_"
          "approx_relative_log2FC_preview.png",
    dpi=300,
    verbose=False,
)
