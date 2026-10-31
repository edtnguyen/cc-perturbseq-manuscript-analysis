from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from plotnine import (
    aes, arrow, coord_cartesian, element_blank, element_line, element_rect,
    element_text, geom_hline, geom_point, geom_segment, geom_text, geom_vline,
    ggplot, labs, scale_color_identity, scale_x_continuous, scale_y_continuous,
    theme, theme_classic,
)

# -----------------------------------------------------------------------------
# Figure 1b: genome-wide screen volcano
# Source: Pos_vs_Neg_volcano_plotnine_MAGECK_FDR_pathway_colored_v8.py
# -----------------------------------------------------------------------------

PANEL_B_PATHWAY_COLORS = {
    "SMAD1": "#7A43A2",
    "SMAD5": "#7A43A2",
    "BMPR1A": "#7A43A2",
    "MIXL1": "#D97A34",
    "EOMES": "#D97A34",
    "KDR": "#059C48",
    "ETV2": "#059C48",
    "PLCG1": "#059C48",
    "KCTD10": "#059C48",
    "PRKCD": "#059C48",
    "NF2": "#ED4A97",
    "TAOK1": "#ED4A97",
    "LATS2": "#ED4A97",
    "AMOTL2": "#ED4A97",
}

PANEL_B_LABEL_POSITIONS = {
    "KDR": (-3.28, 4.15),
    "NF2": (-2.22, 4.56),
    "ETV2": (-2.22, 4.38),
    "PLCG1": (-2.22, 4.20),
    "SMAD1": (-2.82, 3.75),
    "SMAD5": (-2.76, 3.48),
    "BMPR1A": (-2.84, 3.21),
    "KCTD10": (-1.10, 4.60),
    "MIXL1": (-1.10, 4.38),
    "PRKCD": (-1.08, 4.18),
    "TAOK1": (0.10, 4.40),
    "LATS2": (0.10, 4.14),
    "AMOTL2": (0.10, 3.88),
    "EOMES": (-0.34, 2.18),
}


def make_figure_01b(
    mageck_path: str | Path,
    utest_path: str | Path,
):
    """Reproduce Figure 1b from the verified Plotnine script."""
    mageck_path = Path(mageck_path)
    utest_path = Path(utest_path)
    font_family = "Liberation Sans"
    fdr_cutoff = 0.05
    y_cutoff = -np.log10(fdr_cutoff)
    x_limits = (-3.45, 2.55)
    y_limits = (-0.05, 4.85)

    mpl.rcParams.update(
        {
            "font.family": font_family,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": True,
        }
    )

    mageck = pd.read_csv(mageck_path, sep="\t")
    utest = pd.read_csv(utest_path, sep="\t").rename(
        columns={
            "p_value": "u_test_p_value",
            "log2FC": "u_test_log2FC",
            "p_value_corrected": "u_test_fdr",
        }
    )
    required_mageck = {"id", "neg|fdr", "pos|fdr"}
    required_utest = {"Target", "u_test_log2FC"}
    if missing := required_mageck.difference(mageck.columns):
        raise ValueError(f"MAGeCK table is missing columns: {sorted(missing)}")
    if missing := required_utest.difference(utest.columns):
        raise ValueError(f"U-test table is missing columns: {sorted(missing)}")

    plot_df = mageck.merge(utest, left_on="id", right_on="Target", how="outer")
    plot_df["id"] = plot_df["id"].fillna(plot_df["Target"])
    plot_df = plot_df.drop(columns=["Target"])
    plot_df["gene"] = (
        plot_df["id"]
        .astype(str)
        .str.replace(r"_-.*$", "", regex=True)
        .str.replace(r"-P1P2$", "", regex=True)
    )
    plot_df["x"] = pd.to_numeric(plot_df["u_test_log2FC"], errors="coerce")
    neg_fdr = pd.to_numeric(plot_df["neg|fdr"], errors="coerce")
    pos_fdr = pd.to_numeric(plot_df["pos|fdr"], errors="coerce")
    plot_df["mageck_fdr"] = np.where(plot_df["x"] < 0, neg_fdr, pos_fdr)
    raw_fdr = pd.to_numeric(plot_df["mageck_fdr"], errors="coerce").clip(
        lower=np.finfo(float).tiny
    )
    plot_df["y"] = -np.log10(raw_fdr)
    plot_df = plot_df.dropna(subset=["x", "y"]).copy()

    control_names = {"non-targeting", "safe-targeting"}
    plot_df["category"] = "Gene"
    plot_df.loc[plot_df["gene"].isin(control_names), "category"] = "Control"
    mageck_hit = (neg_fdr < fdr_cutoff) | (pos_fdr < fdr_cutoff)
    plot_df.loc[
        mageck_hit & ~plot_df["gene"].isin(control_names), "category"
    ] = "MAGeCK hit"
    plot_df.loc[
        plot_df["gene"].isin(PANEL_B_PATHWAY_COLORS), "category"
    ] = "Gene highlight"

    background = plot_df[plot_df["category"] == "Gene"].copy()
    controls = plot_df[plot_df["category"] == "Control"].copy()
    hits = plot_df[plot_df["category"] == "MAGeCK hit"].copy()
    negative_hits = hits[hits["x"] < 0].copy()
    positive_hits = hits[hits["x"] > 0].copy()
    highlighted = plot_df[plot_df["category"] == "Gene highlight"].copy()
    for gene, color in PANEL_B_PATHWAY_COLORS.items():
        highlighted.loc[highlighted["gene"] == gene, "point_color"] = color

    best = (
        highlighted.sort_values("mageck_fdr")
        .drop_duplicates("gene")
        .set_index("gene")
    )
    label_rows = []
    for gene, (label_x, label_y) in PANEL_B_LABEL_POSITIONS.items():
        if gene not in best.index:
            continue
        row = best.loc[gene]
        point_x, point_y = float(row["x"]), float(row["y"])
        dx, dy = point_x - label_x, point_y - label_y
        distance = float(np.hypot(dx, dy))
        if distance > 0:
            gap = min(0.18, distance * 0.45)
            line_xend = label_x + gap * dx / distance
            line_yend = label_y + gap * dy / distance
        else:
            line_xend, line_yend = label_x, label_y
        label_rows.append(
            {
                "gene": gene,
                "point_x": point_x,
                "point_y": point_y,
                "label_x": label_x,
                "label_y": label_y,
                "line_xend": line_xend,
                "line_yend": line_yend,
                "label_color": PANEL_B_PATHWAY_COLORS[gene],
            }
        )
    labels_df = pd.DataFrame(label_rows)

    arrow_df = pd.DataFrame(
        [
            {"x": -1.35, "y": 0.45, "xend": x_limits[0] + 0.03, "yend": 0.45},
            {"x": 0.75, "y": 0.45, "xend": x_limits[1] - 0.03, "yend": 0.45},
        ]
    )
    arrow_text_df = pd.DataFrame(
        [
            {"x": -2.40, "y": 0.17, "label": "Enriched in\nnegative fraction"},
            {"x": 1.68, "y": 0.17, "label": "Enriched in\npositive fraction"},
        ]
    )

    volcano = (
        ggplot()
        + geom_vline(
            xintercept=0, color="#B3BAC3", linetype="dashed", size=1.1
        )
        + geom_hline(
            yintercept=y_cutoff, color="#B3BAC3", linetype="dashed", size=1.1
        )
        + geom_point(
            background,
            aes("x", "y"),
            color="#D5D5D5",
            size=2.2,
            alpha=0.52,
            stroke=0,
            raster=True,
        )
        + geom_point(
            controls,
            aes("x", "y"),
            color="#AEB3D2",
            size=2.6,
            alpha=0.72,
            stroke=0,
        )
        + geom_point(
            negative_hits,
            aes("x", "y"),
            color="#111111",
            size=2.8,
            alpha=0.95,
            stroke=0,
        )
        + geom_point(
            positive_hits,
            aes("x", "y"),
            color="#8F0000",
            size=2.8,
            alpha=0.95,
            stroke=0,
        )
        + geom_segment(
            labels_df,
            aes(
                x="point_x",
                y="point_y",
                xend="line_xend",
                yend="line_yend",
                color="label_color",
            ),
            inherit_aes=False,
            size=0.48,
            alpha=0.95,
        )
        + geom_point(
            highlighted,
            aes("x", "y", color="point_color"),
            size=6.6,
            alpha=1.0,
            stroke=0,
        )
        + geom_text(
            labels_df,
            aes(
                x="label_x",
                y="label_y",
                label="gene",
                color="label_color",
            ),
            inherit_aes=False,
            family=font_family,
            fontstyle="italic",
            size=11.2,
        )
        + geom_segment(
            arrow_df,
            aes(x="x", y="y", xend="xend", yend="yend"),
            inherit_aes=False,
            color="black",
            size=0.55,
            arrow=arrow(length=0.10, ends="last", type="closed"),
        )
        + geom_text(
            arrow_text_df,
            aes(x="x", y="y", label="label"),
            inherit_aes=False,
            family=font_family,
            size=10.2,
            lineheight=0.95,
            color="#171717",
        )
        + scale_color_identity()
        + scale_x_continuous(breaks=np.arange(-3, 3, 1), expand=(0, 0))
        + scale_y_continuous(breaks=np.arange(0, 5, 1), expand=(0, 0))
        + coord_cartesian(xlim=x_limits, ylim=y_limits)
        + labs(x="Log₂(positive / negative)", y="−Log₁₀(MaGECK FDR)")
        + theme_classic(base_family=font_family, base_size=14)
        + theme(
            figure_size=(12.0, 6.2),
            axis_title_x=element_text(size=17, margin={"t": 10}),
            axis_title_y=element_text(size=17, margin={"r": 10}),
            axis_text=element_text(size=13, color="#111111"),
            axis_ticks=element_line(color="#111111", size=0.8),
            axis_line=element_line(color="#111111", size=1.0),
            legend_position="none",
            panel_grid=element_blank(),
            plot_background=element_rect(fill="white", color="white"),
            panel_background=element_rect(fill="white", color="white"),
        )
    )

    figure = volcano.draw(show=False)
    figure.subplots_adjust(top=0.82, left=0.12, right=0.985, bottom=0.17)
    axis = figure.axes[0]
    axis_box = axis.get_position()
    header_y = axis_box.y1 + 0.02
    header_height = 0.072
    figure.add_artist(
        Rectangle(
            (axis_box.x0, header_y),
            axis_box.width,
            header_height,
            transform=figure.transFigure,
            facecolor="#C8D7DF",
            edgecolor="none",
            clip_on=False,
            zorder=2,
        )
    )
    figure.text(
        axis_box.x0 + axis_box.width / 2,
        header_y + header_height / 2,
        "Positive versus negative fraction genome-wide screen",
        ha="center",
        va="center",
        family=font_family,
        fontsize=20,
        color="#171717",
        zorder=3,
    )
    figure.text(
        axis_box.x0 - 0.065,
        header_y + header_height / 2,
        "b",
        ha="center",
        va="center",
        family=font_family,
        fontsize=22,
        fontweight="bold",
        color="#111111",
    )
    return figure, plot_df, labels_df
