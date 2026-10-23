from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    annotate,
    coord_fixed,
    element_rect,
    element_text,
    geom_raster,
    geom_rect,
    geom_segment,
    geom_text,
    geom_tile,
    ggplot,
    scale_fill_gradientn,
    theme,
    theme_void,
)

from .figure2_common import FIGURE2_PROGRAMS, FIGURE2_TARGET_GROUPS, significance_stars
from .figure2_data import cluster_effect_matrix

DIVERGING_COLORS = [
    "#213460", "#2A62A5", "#71A6CA", "#CDE2ED", "#FAFBFA",
    "#F2DAD2", "#E3AA97", "#CA6E5D", "#9A2A30", "#621422",
]


def _gene_y_positions():
    y_by_gene: dict[str, float] = {}
    center_by_group: dict[str, float] = {}
    title_y_by_group: dict[str, float] = {}
    y = 48.8
    for group_index, (group_name, genes) in enumerate(FIGURE2_TARGET_GROUPS.items()):
        ys = []
        for gene in genes:
            y_by_gene[gene] = y
            ys.append(y)
            y -= 1.9
        center_by_group[group_name] = float(np.mean(ys))
        title_y_by_group[group_name] = max(ys) + 1.7
        if group_index < len(FIGURE2_TARGET_GROUPS) - 1:
            y -= 2.0
    return y_by_gene, center_by_group, title_y_by_group


def plot_figure2a(matrix: pd.DataFrame):
    """Plot Figure 2a with plotnine."""
    correlation, order, _ = cluster_effect_matrix(matrix)
    rank = {target: index for index, target in enumerate(order)}
    y_by_gene, center_by_group, title_y_by_group = _gene_y_positions()

    panel_width, panel_height = 68.0, 60.0
    main_x0, main_y0, main_size = 2.2, 17.3, 31.5
    inset_x0, inset_cell, gene_label_x = 43.5, 1.78, 54.7
    n = len(order)
    x = main_x0 + (np.arange(n) + 0.5) * main_size / n
    y = main_y0 + main_size - (np.arange(n) + 0.5) * main_size / n
    xx, yy = np.meshgrid(x, y)
    main_df = pd.DataFrame(
        {"x": xx.ravel(), "y": yy.ravel(), "fill_value": correlation.to_numpy().ravel()}
    )

    inset_rows: list[dict[str, Any]] = []
    box_rows: list[dict[str, Any]] = []
    segment_rows: list[dict[str, Any]] = []
    inset_records: list[dict[str, Any]] = []
    elbow_1 = [35.3, 35.8, 36.3, 36.8]
    elbow_2 = [40.8, 41.0, 41.2, 41.4]

    for group_index, (group_name, genes) in enumerate(FIGURE2_TARGET_GROUPS.items()):
        group_ranks = sorted(rank[gene] for gene in genes)
        if group_ranks != list(range(group_ranks[0], group_ranks[-1] + 1)):
            raise ValueError(f"Highlighted group is not contiguous after clustering: {group_name}")
        subset = correlation.loc[genes, genes]
        for row_index, row_gene in enumerate(genes):
            for col_index, col_gene in enumerate(genes):
                value = float(subset.loc[row_gene, col_gene])
                inset_rows.append(
                    {
                        "group": group_name,
                        "row_gene": row_gene,
                        "column_gene": col_gene,
                        "x": inset_x0 + (col_index + 0.5) * inset_cell,
                        "y": y_by_gene[row_gene],
                        "fill_value": value,
                    }
                )
                inset_records.append(
                    {
                        "group": group_name,
                        "row_gene": row_gene,
                        "column_gene": col_gene,
                        "pearson_r": value,
                    }
                )

        start, end = group_ranks[0], group_ranks[-1]
        cell = main_size / n
        xmin = main_x0 + start * cell
        xmax = main_x0 + (end + 1) * cell
        ymax = main_y0 + main_size - start * cell
        ymin = main_y0 + main_size - (end + 1) * cell
        box_rows.append(
            {"group": group_name, "xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax,
             "start_rank": start, "end_rank": end}
        )
        points = [
            (xmax, (ymin + ymax) / 2),
            (elbow_1[group_index], (ymin + ymax) / 2),
            (elbow_2[group_index], center_by_group[group_name]),
            (inset_x0 - 0.55, center_by_group[group_name]),
        ]
        for (x1, y1), (x2, y2) in zip(points[:-1], points[1:]):
            segment_rows.append(
                {"group": group_name, "x": x1, "y": y1, "xend": x2, "yend": y2}
            )

    inset_df = pd.DataFrame(inset_rows)
    box_df = pd.DataFrame(box_rows)
    segment_df = pd.DataFrame(segment_rows)
    cb_values = np.linspace(-1.0, 1.0, 181)
    cb_x0, cb_x1 = 43.5, 55.8
    cb_step = (cb_x1 - cb_x0) / len(cb_values)
    cb_df = pd.DataFrame(
        {
            "x": cb_x0 + (np.arange(len(cb_values)) + 0.5) * cb_step,
            "y": np.full(len(cb_values), 11.1),
            "fill_value": cb_values,
        }
    )

    plot = (
        ggplot()
        + geom_raster(main_df, aes("x", "y", fill="fill_value"), inherit_aes=False)
        + geom_rect(
            pd.DataFrame([{"xmin": main_x0, "xmax": main_x0 + main_size,
                           "ymin": main_y0, "ymax": main_y0 + main_size}]),
            aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"),
            inherit_aes=False, fill=None, color="#B8B8B8", size=0.30,
        )
        + geom_tile(
            inset_df, aes("x", "y", fill="fill_value"), inherit_aes=False,
            width=inset_cell, height=inset_cell, color="white", size=0.35,
        )
        + geom_segment(
            segment_df, aes(x="x", y="y", xend="xend", yend="yend"),
            inherit_aes=False, color="black", size=0.48,
        )
        + geom_rect(
            box_df, aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"),
            inherit_aes=False, fill=None, color="#555555", size=0.56,
        )
        + geom_tile(
            cb_df, aes("x", "y", fill="fill_value"), inherit_aes=False,
            width=cb_step + 0.006, height=1.2, color=None,
        )
        + scale_fill_gradientn(colors=DIVERGING_COLORS, limits=(-1.0, 1.0), guide=None)
        + coord_fixed(ratio=1.0, xlim=(0, panel_width), ylim=(0, panel_height), expand=False)
        + theme_void()
        + theme(
            figure_size=(6.8, 6.0),
            plot_background=element_rect(fill="white", color="white"),
            panel_background=element_rect(fill="white", color="white"),
            text=element_text(family="Liberation Sans", color="black"),
            plot_margin=0,
        )
        + annotate("text", x=1.2, y=57.1, label="a", family="Liberation Sans", fontweight="bold", size=25)
        + annotate("text", x=17.95, y=53.65,
                   label="297 target genes plus\npooled non-targeting control",
                   family="Liberation Sans", size=13.0)
        + annotate("text", x=17.95, y=14.25,
                   label="Pearson r — effect on program expression",
                   family="Liberation Sans", size=8.9)
        + annotate("text", x=cb_x0, y=9.95, label="-1", family="Liberation Sans", size=7.8)
        + annotate("text", x=(cb_x0 + cb_x1) / 2, y=9.95, label="0", family="Liberation Sans", size=7.8)
        + annotate("text", x=cb_x1, y=9.95, label="1", family="Liberation Sans", size=7.8)
        + annotate("text", x=(cb_x0 + cb_x1) / 2, y=8.15,
                   label="Pearson r", family="Liberation Sans", size=11.8)
    )
    for group_name, genes in FIGURE2_TARGET_GROUPS.items():
        plot += annotate(
            "text", x=inset_x0 + len(genes) * inset_cell / 2,
            y=title_y_by_group[group_name], label=group_name,
            family="Liberation Sans", size=9.7,
        )
        for gene in genes:
            plot += annotate(
                "text", x=gene_label_x, y=y_by_gene[gene], label=gene,
                family="Liberation Sans", size=9.2, ha="left",
            )

    return (
        plot,
        correlation,
        pd.DataFrame({"cluster_rank": np.arange(1, n + 1), "target_name": order}),
        pd.DataFrame(inset_records),
        box_df[["group", "start_rank", "end_rank"]],
    )


def plot_figure2b(selected_associations: pd.DataFrame, ec_stats: pd.DataFrame):
    """Plot Figure 2b with plotnine, including the 120-control EC-output column."""
    targets = [gene for genes in FIGURE2_TARGET_GROUPS.values() for gene in genes]
    y_by_gene, center_by_group, _ = _gene_y_positions()
    x_by_key = {
        ("D1", 13): 13.0, ("D1", 30): 16.2, ("D1", 48): 19.4,
        ("D2", 15): 27.0, ("D2", 14): 30.2, ("D2", 17): 33.4,
        ("D3", 2): 41.0, ("D3", 3): 44.2, ("D3", 8): 47.4,
        ("D3", 11): 50.6, ("D3", 33): 53.8, ("D3", 34): 57.0,
        ("D3", 27): 60.2,
    }
    ec_x = 70.0
    lookup = selected_associations.set_index(["target_name", "Sample", "program_name"])
    ec_lookup = ec_stats.set_index("gene")
    tile_rows: list[dict[str, Any]] = []
    star_rows: list[dict[str, Any]] = []
    program_source_rows: list[dict[str, Any]] = []
    effect_table = pd.DataFrame(index=targets)
    q_table = pd.DataFrame(index=targets)

    for target in targets:
        for item in FIGURE2_PROGRAMS:
            day, program = item["day"], int(item["program"])
            row = lookup.loc[(target, day, program)]
            effect = float(row["approx_log2FC"])
            qvalue = float(row["adj_pval"])
            x = x_by_key[(day, program)]
            column = f"{day}_P{program}"
            effect_table.loc[target, column] = effect
            q_table.loc[target, column] = qvalue
            record = {
                "target": target, "day": day, "program": program,
                "program_label": item["label"], "effect": effect,
                "FDR": qvalue, "significance": significance_stars(qvalue),
            }
            program_source_rows.append(record)
            tile_rows.append(
                {**record, "x": x, "y": y_by_gene[target],
                 "display_value": float(np.clip(effect, -3, 3))}
            )
            if significance_stars(qvalue):
                star_rows.append(
                    {"x": x, "y": y_by_gene[target],
                     "label": significance_stars(qvalue),
                     "text_color": "white" if abs(effect) >= 1.45 else "black"}
                )

        ec_row = ec_lookup.loc[target]
        effect = float(ec_row["effect_log2_target_over_control"])
        qvalue = float(ec_row["FDR_BH_297_genes"])
        effect_table.loc[target, "D3_EC_log2_target_over_120_controls"] = effect
        q_table.loc[target, "D3_EC_log2_target_over_120_controls"] = qvalue
        tile_rows.append(
            {"target": target, "day": "D3", "program": np.nan,
             "program_label": "D3 EC log2(target / 120 controls)",
             "effect": effect, "FDR": qvalue,
             "significance": significance_stars(qvalue),
             "x": ec_x, "y": y_by_gene[target],
             "display_value": float(np.clip(effect, -3, 3))}
        )
        if significance_stars(qvalue):
            star_rows.append(
                {"x": ec_x, "y": y_by_gene[target],
                 "label": significance_stars(qvalue),
                 "text_color": "white" if abs(effect) >= 1.45 else "black"}
            )

    tile_df = pd.DataFrame(tile_rows)
    stars = pd.DataFrame(star_rows, columns=["x", "y", "label", "text_color"])
    white_stars = stars.loc[stars["text_color"].eq("white")]
    black_stars = stars.loc[stars["text_color"].eq("black")]

    cb_values = np.linspace(-3.0, 3.0, 241)
    cb_x0, cb_x1 = 77.0, 98.0
    cb_step = (cb_x1 - cb_x0) / len(cb_values)
    cb_df = pd.DataFrame(
        {"x": cb_x0 + (np.arange(len(cb_values)) + 0.5) * cb_step,
         "y": np.full(len(cb_values), 41.1), "display_value": cb_values}
    )
    header_lines = pd.DataFrame(
        [
            {"x": 11.3, "xend": 21.1, "y": 52.4, "yend": 52.4},
            {"x": 25.3, "xend": 35.1, "y": 52.4, "yend": 52.4},
            {"x": 39.3, "xend": 51.9, "y": 52.4, "yend": 52.4},
            {"x": 52.3, "xend": 61.9, "y": 52.4, "yend": 52.4},
            {"x": 68.3, "xend": 71.7, "y": 52.4, "yend": 52.4},
        ]
    )

    plot = (
        ggplot()
        + geom_tile(
            tile_df, aes("x", "y", fill="display_value"), inherit_aes=False,
            width=3.05, height=1.78, color="white", size=0.35,
        )
        + geom_tile(
            cb_df, aes("x", "y", fill="display_value"), inherit_aes=False,
            width=cb_step + 0.004, height=1.3, color=None,
        )
        + geom_text(
            white_stars, aes("x", "y", label="label"), inherit_aes=False,
            color="white", family="Liberation Sans", size=5.6, fontweight="bold",
        )
        + geom_text(
            black_stars, aes("x", "y", label="label"), inherit_aes=False,
            color="black", family="Liberation Sans", size=5.6, fontweight="bold",
        )
        + geom_segment(
            header_lines, aes(x="x", y="y", xend="xend", yend="yend"),
            inherit_aes=False, color="#5A5A5A", size=0.48,
        )
        + scale_fill_gradientn(colors=DIVERGING_COLORS, limits=(-3.0, 3.0), guide=None)
        + coord_fixed(ratio=1.0, xlim=(0, 106), ylim=(0, 60), expand=False)
        + theme_void()
        + theme(
            figure_size=(10.6, 6.0),
            plot_background=element_rect(fill="white", color="white"),
            panel_background=element_rect(fill="white", color="white"),
            text=element_text(family="Liberation Sans", color="black"),
            plot_margin=0,
        )
        + annotate("text", x=1.35, y=56.55, label="b",
                   family="Liberation Sans", fontweight="bold", size=25)
        + annotate("text", x=36.5, y=57.2, label="Effects on gene programs",
                   family="Liberation Sans", fontweight="bold", size=15.5)
        + annotate("text", x=70.0, y=57.25, label="Effects on EC\ndifferentiation",
                   family="Liberation Sans", fontweight="bold", size=15.0)
        + annotate("text", x=16.2, y=54.0, label="Day 1",
                   family="Liberation Sans", fontweight="bold", size=10.5)
        + annotate("text", x=30.2, y=54.0, label="Day 2",
                   family="Liberation Sans", fontweight="bold", size=10.5)
        + annotate("text", x=46.6, y=54.0, label="Day 3 — EC",
                   family="Liberation Sans", fontweight="bold", size=10.5)
        + annotate("text", x=57.0, y=54.0, label="Non-EC",
                   family="Liberation Sans", fontweight="bold", size=10.0)
        + annotate("text", x=77.0, y=45.0, label="Effect / EC ratio",
                   family="Liberation Sans", fontweight="bold", size=10.8, ha="left")
        + annotate("text", x=77.0, y=43.75, label="log2 scale",
                   family="Liberation Sans", size=9.7, ha="left")
        + annotate("text", x=cb_x0, y=39.85, label="-3", family="Liberation Sans", size=7.2)
        + annotate("text", x=(cb_x0 + cb_x1) / 2, y=39.85, label="0", family="Liberation Sans", size=7.2)
        + annotate("text", x=cb_x1, y=39.85, label="+3", family="Liberation Sans", size=7.2)
        + annotate("text", x=77.0, y=36.8, label="FDR",
                   family="Liberation Sans", fontweight="bold", size=10.8, ha="left")
        + annotate("text", x=77.0, y=35.45,
                   label="* < 0.05\n** < 0.01\n*** < 0.001",
                   family="Liberation Sans", size=10.2, ha="left", va="top")
        + annotate("text", x=77.0, y=29.1, label="D3 EC ratio",
                   family="Liberation Sans", fontweight="bold", size=10.5, ha="left")
        + annotate("text", x=77.0, y=27.8,
                   label="log2(target mean / 120-control mean)",
                   family="Liberation Sans", size=9.2, ha="left")
        + annotate("text", x=77.0, y=26.5,
                   label="Welch t-test; BH across 297 targets",
                   family="Liberation Sans", size=8.7, ha="left")
    )

    for group_name, genes in FIGURE2_TARGET_GROUPS.items():
        plot += annotate(
            "text", x=0.35, y=center_by_group[group_name], label=group_name,
            family="Liberation Sans", fontweight="bold", size=9.2, ha="left",
        )
        for gene in genes:
            plot += annotate(
                "text", x=9.7, y=y_by_gene[gene], label=gene,
                family="Liberation Sans", size=10.4, ha="right",
            )
    for item in FIGURE2_PROGRAMS:
        plot += annotate(
            "text", x=x_by_key[(item["day"], int(item["program"]))] + 0.2,
            y=13.25, label=item["label"], family="Liberation Sans",
            size=8.0, angle=55, ha="right", va="top",
        )
    plot += annotate(
        "text", x=ec_x + 0.2, y=13.25,
        label="D3 EC\nlog2 ratio vs 120 controls",
        family="Liberation Sans", size=8.5, angle=55, ha="right", va="top",
    )

    return (
        plot,
        effect_table,
        q_table,
        pd.DataFrame(program_source_rows),
        ec_stats.copy(),
    )
