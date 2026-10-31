from __future__ import annotations

from pathlib import Path
import hashlib

import matplotlib as mpl
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    arrow,
    coord_cartesian,
    element_blank,
    element_line,
    element_rect,
    element_text,
    geom_boxplot,
    geom_hline,
    geom_point,
    geom_rect,
    geom_segment,
    geom_text,
    geom_vline,
    ggplot,
    labs,
    scale_color_identity,
    scale_fill_manual,
    scale_x_continuous,
    scale_y_continuous,
    theme,
    theme_classic,
)

from .figure1_panel_b import PANEL_B_LABEL_POSITIONS, PANEL_B_PATHWAY_COLORS
from .figure1_panel_c import (
    PANEL_C_CONTROL_TYPES,
    PANEL_C_GROUP_COLORS,
    PANEL_C_GROUP_LABELS,
    PANEL_C_GROUPS,
    PANEL_C_VALUE_COL,
    compute_panel_c_welch_statistics,
    load_validation_screen_counts,
    reconstruct_validation_guide_ec_percentages,
    sha256_file,
    split_panel_c_guides,
)


# These functions retain the previously recovered Figure 1 calculations and
# plotting choices, but return plotnine ggplot objects directly. Matplotlib is
# used only by plotnine as its rendering backend and for editable font settings;
# no pyplot figure construction or post-draw artist editing is used here.


def _set_editable_vector_fonts(font_family: str = "Liberation Sans") -> None:
    mpl.rcParams.update(
        {
            "font.family": font_family,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": True,
            "savefig.transparent": False,
        }
    )


def save_plotnine_panel(
    plot,
    outdir: str | Path,
    stem: str,
    *,
    png_dpi: int = 600,
) -> dict[str, str]:
    """Save one plotnine object without converting it to a pyplot Figure."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, str] = {}
    for suffix in ("pdf", "svg", "png"):
        path = outdir / f"{stem}.{suffix}"
        plot.save(
            filename=str(path),
            dpi=png_dpi if suffix == "png" else 300,
            verbose=False,
            limitsize=False,
        )
        outputs[suffix] = str(path)
    return outputs


def make_figure_01b_plotnine(
    mageck_path: str | Path,
    utest_path: str | Path,
):
    """Return Figure 1b as the source-recovered plotnine ggplot object."""
    mageck_path = Path(mageck_path)
    utest_path = Path(utest_path)
    font_family = "Liberation Sans"
    fdr_cutoff = 0.05
    y_cutoff = -np.log10(fdr_cutoff)
    x_limits = (-3.45, 2.55)
    data_y_limits = (-0.05, 4.85)
    plot_y_limits = (-0.05, 5.20)
    _set_editable_vector_fonts(font_family)

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
    mageck_hit = (
        pd.to_numeric(plot_df["neg|fdr"], errors="coerce").lt(fdr_cutoff)
        | pd.to_numeric(plot_df["pos|fdr"], errors="coerce").lt(fdr_cutoff)
    )
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
    highlighted["point_color"] = highlighted["gene"].map(PANEL_B_PATHWAY_COLORS)

    best = (
        highlighted.sort_values("mageck_fdr")
        .drop_duplicates("gene")
        .set_index("gene")
    )
    label_rows: list[dict[str, object]] = []
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
    header_df = pd.DataFrame(
        [
            {
                "xmin": x_limits[0],
                "xmax": x_limits[1],
                "ymin": 4.92,
                "ymax": plot_y_limits[1],
            }
        ]
    )
    header_text_df = pd.DataFrame(
        [
            {
                "x": float(np.mean(x_limits)),
                "y": 5.06,
                "label": "Positive versus negative fraction genome-wide screen",
            },
            {"x": x_limits[0] + 0.12, "y": 5.06, "label": "b"},
        ]
    )

    plot = (
        ggplot()
        + geom_rect(
            header_df,
            aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"),
            inherit_aes=False,
            fill="#C8D7DF",
            color=None,
        )
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
        + geom_text(
            header_text_df.iloc[[0]],
            aes(x="x", y="y", label="label"),
            inherit_aes=False,
            family=font_family,
            size=14.0,
            color="#171717",
            fontweight="normal",
        )
        + geom_text(
            header_text_df.iloc[[1]],
            aes(x="x", y="y", label="label"),
            inherit_aes=False,
            family=font_family,
            size=15.0,
            color="#111111",
            fontweight="bold",
        )
        + scale_color_identity()
        + scale_x_continuous(breaks=np.arange(-3, 3, 1), expand=(0, 0))
        + scale_y_continuous(breaks=np.arange(0, 5, 1), expand=(0, 0))
        + coord_cartesian(xlim=x_limits, ylim=plot_y_limits)
        + labs(x="Log₂(positive / negative)", y="−Log₁₀(MaGECK FDR)")
        + theme_classic(base_family=font_family, base_size=14)
        + theme(
            figure_size=(12.0, 6.6),
            axis_title_x=element_text(size=17, margin={"t": 10}),
            axis_title_y=element_text(size=17, margin={"r": 10}),
            axis_text=element_text(size=13, color="#111111"),
            axis_ticks=element_line(color="#111111", size=0.8),
            axis_line=element_line(color="#111111", size=1.0),
            legend_position="none",
            panel_grid=element_blank(),
            plot_background=element_rect(fill="white", color="white"),
            panel_background=element_rect(fill="white", color="white"),
            plot_margin=0.02,
        )
    )
    plot_df.attrs["source_data_y_limits"] = data_y_limits
    return plot, plot_df, labels_df


def plot_figure_01c_plotnine(
    controls: pd.DataFrame,
    targeting: dict[str, pd.DataFrame],
    stats: pd.DataFrame,
    *,
    control_baseline_pct: float = 70.0,
):
    """Render the validation-screen guide distributions directly with plotnine."""
    font_family = "Liberation Sans"
    _set_editable_vector_fonts(font_family)
    stats_index = stats.set_index("gene")

    guide_frames: list[pd.DataFrame] = []
    jitter_rows: list[dict[str, object]] = []
    position_rows: list[dict[str, object]] = []
    span_rows: list[dict[str, object]] = []
    star_rows: list[dict[str, object]] = []
    x = 0.0

    for group_name, genes in PANEL_C_GROUPS:
        start = x
        for gene in genes:
            source = controls if gene == "NTC" else targeting[gene]
            frame = source.copy()
            digest = int.from_bytes(
                hashlib.sha256(gene.encode()).digest()[:4], "little"
            )
            rng = np.random.default_rng(np.random.SeedSequence([7, digest]))
            spread = 0.09 if gene == "NTC" else 0.055
            jitter = np.clip(rng.normal(0, spread, len(frame)), -0.16, 0.16)
            frame["plot_gene"] = gene
            frame["group"] = group_name
            frame["x_center"] = x
            frame["x_jitter"] = x + jitter
            guide_frames.append(frame)
            position_rows.append(
                {
                    "plot_gene": gene,
                    "group": group_name,
                    "x": x,
                    "label": f"{gene}\nn = {len(frame)}",
                }
            )
            for (_, row), offset in zip(frame.iterrows(), jitter):
                jitter_rows.append(
                    {
                        "plot_gene": gene,
                        "guide_id": row["guide_id"],
                        PANEL_C_VALUE_COL: float(row[PANEL_C_VALUE_COL]),
                        "x_center": x,
                        "x_jitter": float(offset),
                        "x": float(x + offset),
                        "group": group_name,
                    }
                )
            if gene != "NTC":
                annotation = str(stats_index.loc[gene, "significance"])
                if annotation:
                    star_rows.append(
                        {"x": x, "y": 96.5, "label": annotation}
                    )
            x += 1.0
        end = x - 1.0
        span_rows.append(
            {
                "group": group_name,
                "xmin": start - 0.34,
                "xmax": end + 0.34,
                "x": (start + end) / 2,
                "ymin": 102.0,
                "ymax": 106.5,
                "label_y": 109.0,
                "label": PANEL_C_GROUP_LABELS[group_name],
            }
        )
        x += 0.90

    guide_plot = pd.concat(guide_frames, ignore_index=True)
    positions = pd.DataFrame(position_rows)
    spans = pd.DataFrame(span_rows)
    stars = pd.DataFrame(star_rows, columns=["x", "y", "label"])
    jitter_table = pd.DataFrame(jitter_rows)
    x_min = float(positions["x"].min() - 0.6)
    x_max = float(positions["x"].max() + 0.6)

    plot = (
        ggplot()
        + geom_hline(
            yintercept=control_baseline_pct,
            color="#9EA3A8",
            linetype="dashed",
            size=0.8,
        )
        + geom_rect(
            spans,
            aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="group"),
            inherit_aes=False,
            color="black",
            size=0.35,
            alpha=0.22,
        )
        + geom_boxplot(
            guide_plot,
            aes(
                x="x_center",
                y=PANEL_C_VALUE_COL,
                group="plot_gene",
                fill="group",
            ),
            inherit_aes=False,
            width=0.42,
            outlier_alpha=0,
            color="black",
            size=0.75,
            alpha=0.86,
        )
        + geom_point(
            guide_plot,
            aes(x="x_jitter", y=PANEL_C_VALUE_COL),
            inherit_aes=False,
            color="#707070",
            size=1.55,
            alpha=0.50,
            stroke=0,
        )
        + geom_text(
            stars,
            aes(x="x", y="y", label="label"),
            inherit_aes=False,
            family=font_family,
            size=10.0,
            fontweight="bold",
        )
        + geom_text(
            spans,
            aes(x="x", y="label_y", label="label"),
            inherit_aes=False,
            family=font_family,
            size=8.5,
            fontweight="bold",
            lineheight=0.95,
        )
        + scale_fill_manual(values=PANEL_C_GROUP_COLORS, guide=None)
        + scale_x_continuous(
            breaks=positions["x"].tolist(),
            labels=positions["label"].tolist(),
            limits=(x_min, x_max),
            expand=(0, 0),
        )
        + scale_y_continuous(
            breaks=np.arange(0, 101, 20),
            limits=(0, 112),
            expand=(0, 0),
        )
        + labs(
            title="c    Validation screen",
            subtitle=(
                "Two-sided Welch t-test; Benjamini–Hochberg correction across "
                "all 297 target genes"
            ),
            x="",
            y="Estimated D3 EC cells (%)",
            caption="* FDR < 0.05    ** FDR < 0.01    *** FDR < 0.001",
        )
        + theme_classic(base_family=font_family, base_size=12)
        + theme(
            figure_size=(13.4, 6.5),
            axis_text_x=element_text(rotation=58, ha="right", size=9.5),
            axis_text_y=element_text(size=10.5),
            axis_title_x=element_blank(),
            axis_title_y=element_text(size=15, margin={"r": 8}),
            axis_ticks=element_line(color="#111111", size=0.8),
            axis_line=element_line(color="#111111", size=1.0),
            plot_title=element_text(size=17, weight="bold", color="#111111"),
            plot_subtitle=element_text(size=9.5, color="#222222"),
            plot_caption=element_text(size=9.0, ha="right", color="#222222"),
            legend_position="none",
            panel_grid=element_blank(),
            plot_background=element_rect(fill="white", color="white"),
            panel_background=element_rect(fill="white", color="white"),
            plot_margin=0.02,
        )
    )
    return plot, jitter_table


def run_figure_01c_plotnine(
    validation_count_path: str | Path,
    outdir: str | Path,
    *,
    control_baseline_pct: float = 70.0,
    pseudocount: float = 0.5,
):
    """Recalculate Figure 1c and export its plotnine rendering and source tables."""
    validation_count_path = Path(validation_count_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    counts = load_validation_screen_counts(validation_count_path)
    guide_table, offsets = reconstruct_validation_guide_ec_percentages(
        counts,
        control_baseline_pct=control_baseline_pct,
        pseudocount=pseudocount,
    )
    controls, targeting = split_panel_c_guides(guide_table)
    stats = compute_panel_c_welch_statistics(controls, targeting)
    selected_genes = [gene for _, genes in PANEL_C_GROUPS[1:] for gene in genes]
    selected_stats = stats.set_index("gene").loc[selected_genes].reset_index()

    plot, jitter = plot_figure_01c_plotnine(
        controls,
        targeting,
        stats,
        control_baseline_pct=control_baseline_pct,
    )
    outputs = save_plotnine_panel(
        plot,
        outdir,
        "Figure_01c_validation_screen_EC_fraction",
        png_dpi=600,
    )

    guide_path = outdir / "Figure_01c_reconstructed_guide_EC_percentages.csv"
    stats_path = outdir / "Figure_01c_Welch_statistics_all_297_genes.csv"
    selected_path = outdir / "Figure_01c_Welch_statistics_plotted_genes.csv"
    jitter_path = outdir / "Figure_01c_plot_jitter_coordinates.csv"
    guide_table.to_csv(guide_path, index=False)
    stats.to_csv(stats_path, index=False)
    selected_stats.to_csv(selected_path, index=False)
    jitter.to_csv(jitter_path, index=False)
    outputs.update(
        {
            "guide_level_reconstruction": str(guide_path),
            "all_gene_statistics": str(stats_path),
            "plotted_gene_statistics": str(selected_path),
            "jitter_coordinates": str(jitter_path),
        }
    )
    return (
        plot,
        stats,
        selected_stats,
        guide_table,
        outputs,
        sha256_file(validation_count_path),
        offsets,
    )
