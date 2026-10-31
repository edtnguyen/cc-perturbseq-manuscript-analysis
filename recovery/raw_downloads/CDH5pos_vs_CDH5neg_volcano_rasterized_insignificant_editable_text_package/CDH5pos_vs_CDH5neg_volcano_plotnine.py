#!/usr/bin/env python3
"""Create a hybrid vector/raster CDH5-positive versus CDH5-negative volcano plot.

All 20,655 targets are shown. The 20,416 nonsignificant gray background points
are rasterized to reduce PDF and SVG size. Significant black/red points,
selected-gene points, leader lines, axes, legend, and all text remain vector.
SVG text is retained as text, and PDF fonts are embedded as CID TrueType fonts
for editing in Adobe Illustrator. Labels use Plotnine's
``geom_text(adjust_text=...)`` workflow.
"""

from __future__ import annotations

import argparse
import math
import zipfile
from pathlib import Path

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    coord_cartesian,
    element_blank,
    element_line,
    element_rect,
    element_text,
    geom_hline,
    geom_point,
    geom_text,
    geom_vline,
    ggplot,
    labs,
    guide_legend,
    guides,
    scale_color_manual,
    scale_x_continuous,
    scale_y_continuous,
    theme,
)

BASE_NAME = "CDH5pos_vs_CDH5neg_volcano_selected_genes_black_gray_red_adjusttext_rasterized_insignificant"
FONT_FAMILY = "Liberation Sans"

SELECTED_GENES = [
    "BMPR1A", "SMAD1", "SMAD5",
    "MIXL1", "EOMES",
    "KDR", "KCTD10", "PRKCD", "ETV2", "PLCG1",
    "TAOK1", "AMOTL2", "LATS2", "NF2",
]

# Labels begin at one common vertical offset above their points. Their final
# positions are determined automatically by adjustText rather than by a
# gene-specific coordinate table.
LABEL_INITIAL_Y_OFFSET = 0.80


def configure_matplotlib() -> None:
    """Keep SVG text editable and embed TrueType fonts in PDF."""
    mpl.rcParams["font.family"] = FONT_FAMILY
    mpl.rcParams["figure.dpi"] = 300
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["pdf.compression"] = 9
    mpl.rcParams["text.usetex"] = False
    # Keep any Matplotlib fallback artists (including legend marker stubs)
    # within the requested black/gray/red palette.
    mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
        color=["#000000", "#B22222", "#BDBDBD"]
    )


def load_input(path: str | Path) -> pd.DataFrame:
    """Read a MAGeCK gene summary from a TSV or from a ZIP archive."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as zf:
            candidates = [
                n for n in zf.namelist()
                if n.endswith("Pos_vs_Neg.gene_summary.txt")
                and "__MACOSX" not in n
                and not Path(n).name.startswith("._")
            ]
            if len(candidates) != 1:
                raise ValueError(
                    f"Expected exactly one Pos_vs_Neg.gene_summary.txt in {path}; "
                    f"found {candidates}"
                )
            with zf.open(candidates[0]) as handle:
                return pd.read_csv(handle, sep="\t")

    return pd.read_csv(path, sep="\t")


def prepare_data(raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Construct plotting fields and the selected-gene annotation table."""
    required = {"id", "neg|lfc", "neg|fdr", "pos|fdr"}
    missing = required.difference(raw.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df = raw.copy()
    df["gene"] = df["id"].str.replace(r"_-.*$", "", regex=True)
    df["log2FC"] = pd.to_numeric(df["neg|lfc"], errors="coerce")
    neg_fdr = pd.to_numeric(df["neg|fdr"], errors="coerce")
    pos_fdr = pd.to_numeric(df["pos|fdr"], errors="coerce")
    df["FDR"] = np.where(df["log2FC"] < 0, neg_fdr, pos_fdr)

    if df[["log2FC", "FDR"]].isna().any().any():
        bad = int(df[["log2FC", "FDR"]].isna().any(axis=1).sum())
        raise ValueError(f"Found {bad} rows with non-numeric log2FC or FDR values")

    positive_fdr = df.loc[df["FDR"] > 0, "FDR"]
    floor = float(positive_fdr.min()) if not positive_fdr.empty else np.finfo(float).tiny
    df["FDR_plot"] = df["FDR"].clip(lower=floor, upper=1.0)
    df["neg_log10_FDR"] = -np.log10(df["FDR_plot"])
    significance_order = [
        "Significant decrease",
        "Significant increase",
        "Not significant",
    ]
    df["significance"] = np.select(
        [
            (df["FDR"] < 0.05) & (df["log2FC"] < 0),
            (df["FDR"] < 0.05) & (df["log2FC"] > 0),
        ],
        ["Significant decrease", "Significant increase"],
        default="Not significant",
    )
    df["significance"] = pd.Categorical(
        df["significance"], categories=significance_order, ordered=True
    )

    selected = df.loc[df["gene"].isin(SELECTED_GENES)].copy()
    duplicate_counts = selected["gene"].value_counts()
    bad_duplicates = duplicate_counts[duplicate_counts != 1]
    missing_selected = sorted(set(SELECTED_GENES).difference(selected["gene"]))
    if not bad_duplicates.empty or missing_selected:
        raise ValueError(
            "Selected-gene mapping is not one-to-one. "
            f"Missing={missing_selected}; duplicate counts={bad_duplicates.to_dict()}"
        )

    selected["label_start_x"] = selected["log2FC"]
    selected["label_start_y"] = selected["neg_log10_FDR"] + LABEL_INITIAL_Y_OFFSET
    selected["requested_order"] = pd.Categorical(
        selected["gene"], categories=SELECTED_GENES, ordered=True
    )
    selected = selected.sort_values("requested_order")
    return df, selected


def build_volcano(df: pd.DataFrame, selected: pd.DataFrame):
    """Build the publication-style Plotnine volcano plot."""
    configure_matplotlib()
    fdr_threshold = -math.log10(0.05)
    significance_order = [
        "Significant decrease",
        "Significant increase",
        "Not significant",
    ]
    counts = df["significance"].value_counts().reindex(significance_order, fill_value=0)
    legend_labels = {
        category: f"{category} (n={int(counts[category]):,})"
        for category in significance_order
    }
    significance_colors = {
        "Significant decrease": "#000000",
        "Significant increase": "#B22222",
        "Not significant": "#BDBDBD",
    }

    insignificant = df.loc[df["significance"] == "Not significant"].copy()
    significant = df.loc[df["significance"] != "Not significant"].copy()
    return (
        ggplot(df, aes("log2FC", "neg_log10_FDR"))
        + geom_hline(
            yintercept=fdr_threshold,
            color="#555555",
            linetype="dashed",
            size=0.65,
        )
        + geom_vline(
            xintercept=0,
            color="#555555",
            linetype="dashed",
            size=0.65,
        )
        # Rasterize only the dense gray background cloud.
        + geom_point(
            data=insignificant,
            mapping=aes("log2FC", "neg_log10_FDR"),
            inherit_aes=False,
            color="#BDBDBD",
            alpha=0.30,
            size=1.10,
            stroke=0,
            raster=True,
            show_legend=False,
        )
        # Keep significant points as editable vector objects.
        + geom_point(
            data=significant,
            mapping=aes("log2FC", "neg_log10_FDR", color="significance"),
            inherit_aes=False,
            alpha=1.0,
            size=1.70,
            stroke=0,
            show_legend=True,
        )
        # Thin black underlay plus colored overlay makes requested points visible.
        + geom_point(
            data=selected,
            mapping=aes("log2FC", "neg_log10_FDR"),
            inherit_aes=False,
            color="black",
            size=3.85,
            stroke=0,
        )
        + geom_point(
            data=selected,
            mapping=aes(
                "log2FC", "neg_log10_FDR", color="significance"
            ),
            inherit_aes=False,
            size=2.80,
            alpha=1.0,
            stroke=0,
            show_legend=False,
        )
        + geom_text(
            data=selected,
            mapping=aes(x="label_start_x", y="label_start_y", label="gene"),
            inherit_aes=False,
            family=FONT_FAMILY,
            size=7.6,
            color="black",
            ha="center",
            va="center",
            adjust_text={
                "x": df.loc[df["significance"] != "Not significant", "log2FC"].to_numpy(),
                "y": df.loc[df["significance"] != "Not significant", "neg_log10_FDR"].to_numpy(),
                "target_x": selected["log2FC"].to_numpy(),
                "target_y": selected["neg_log10_FDR"].to_numpy(),
                "avoid_self": True,
                "prevent_crossings": True,
                "force_text": (0.80, 1.05),
                "force_static": (0.20, 0.45),
                "force_pull": (0.0, 0.0),
                "force_explode": (0.55, 0.90),
                "pull_threshold": 12,
                "expand": (1.18, 1.45),
                "max_move": (15, 15),
                "ensure_inside_axes": True,
                "only_move": {
                    "text": "xy",
                    "static": "xy",
                    "explode": "xy",
                    "pull": "xy",
                },
                "min_arrow_len": 0,
                "iter_lim": 800,
                "arrowprops": {
                    "arrowstyle": "-",
                    "linewidth": 0.45,
                    "color": "#333333",
                    "shrinkA": 2.0,
                    "shrinkB": 4.0,
                },
            },
        )
        + scale_color_manual(
            values=significance_colors,
            breaks=significance_order,
            labels=[legend_labels[x] for x in significance_order],
            name=None,
            drop=False,
        )
        + guides(
            color=guide_legend(
                override_aes={"size": 3.6, "alpha": 1.0, "stroke": 0}
            )
        )
        + scale_x_continuous(
            breaks=[-3, -2, -1, 0, 1, 2],
            expand=(0, 0),
        )
        + scale_y_continuous(
            breaks=[0, 1, 2, 3, 4, 5, 6],
            expand=(0, 0),
        )
        + coord_cartesian(xlim=(-3.85, 2.75), ylim=(-0.05, 6.34))
        + labs(
            x=r"$\log_{2}\mathrm{FC}$ (gRNA enrichment in CDH5$^{+}$ vs CDH5$^{-}$)",
            y=r"$-\log_{10}(\mathrm{FDR})$",
        )
        + theme(
            figure_size=(6.5, 4.5),
            panel_background=element_rect(fill="white", color=None),
            plot_background=element_rect(fill="white", color=None),
            panel_border=element_rect(fill=None, color="black", size=0.75),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            axis_line=element_blank(),
            axis_ticks=element_line(color="black", size=0.75),
            axis_ticks_length=4,
            axis_text=element_text(
                family=FONT_FAMILY, size=9.8, color="black"
            ),
            axis_title_x=element_text(
                family=FONT_FAMILY,
                size=11.8,
                color="black",
                margin={"t": 8},
            ),
            axis_title_y=element_text(
                family=FONT_FAMILY,
                size=11.8,
                color="black",
                margin={"r": 8},
            ),
            plot_margin=0.025,
            legend_position=(0.825, 0.865),
            legend_justification=(0.5, 0.5),
            legend_direction="vertical",
            legend_background=element_rect(
                fill="white", color="#666666", size=0.55
            ),
            legend_key=element_rect(fill="white", color=None),
            legend_title=element_blank(),
            legend_text=element_text(
                family=FONT_FAMILY, size=8.6, color="black"
            ),
            legend_margin=5,
            legend_box_margin=0,
        )
    )


def save_outputs(
    plot,
    df: pd.DataFrame,
    selected: pd.DataFrame,
    outdir: str | Path,
) -> dict[str, Path]:
    """Save vector, high-resolution raster, and source-data outputs."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    paths = {
        "pdf": outdir / f"{BASE_NAME}.pdf",
        "svg": outdir / f"{BASE_NAME}.svg",
        "png": outdir / f"{BASE_NAME}.png",
        "plot_data": outdir / f"{BASE_NAME}_plot_data.csv.gz",
        "selected_stats": outdir / f"{BASE_NAME}_selected_gene_statistics.csv",
    }

    # Draw once so adjustText computes one deterministic arrangement, then
    # reuse that exact arrangement for PDF, SVG, and the 600-dpi PNG.
    fig = plot.draw(show=False)
    fig.set_size_inches(6.5, 4.5, forward=True)
    fig.savefig(paths["pdf"], format="pdf", dpi=300)
    fig.savefig(paths["svg"], format="svg", dpi=300)
    fig.savefig(paths["png"], format="png", dpi=600)
    plt.close(fig)

    plot_cols = [
        "id", "gene", "log2FC", "FDR", "neg_log10_FDR", "significance"
    ]
    df.loc[:, plot_cols].to_csv(paths["plot_data"], index=False, compression="gzip")

    stats_cols = [
        "gene", "id", "log2FC", "FDR", "neg_log10_FDR", "significance"
    ]
    selected.loc[:, stats_cols].to_csv(paths["selected_stats"], index=False)
    return paths


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=script_dir / "source" / "Pos_vs_Neg.gene_summary.txt",
        help="MAGeCK gene-summary TSV or ZIP containing it.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=script_dir,
        help="Output directory.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raw = load_input(args.input)
    df, selected = prepare_data(raw)
    plot = build_volcano(df, selected)
    paths = save_outputs(plot, df, selected, args.outdir)

    counts = df["significance"].value_counts().to_dict()
    print(f"Rows plotted: {len(df):,}")
    print(f"Selected genes labeled: {len(selected)}")
    print(f"Significance counts: {counts}")
    for key, path in paths.items():
        print(f"{key}: {path}")


if __name__ == "__main__":
    main()
