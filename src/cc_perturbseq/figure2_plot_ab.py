from __future__ import annotations

from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import ConnectionPatch, Rectangle
import numpy as np
import pandas as pd

from .figure2_common import (
    FIGURE2_DIVERGING, FIGURE2_PROGRAMS, FIGURE2_TARGET_GROUPS, significance_stars,
)
from .figure2_data import cluster_effect_matrix

def plot_figure2a(matrix: pd.DataFrame):
    correlation, order, _ = cluster_effect_matrix(matrix)
    ranks = {gene: rank for rank, gene in enumerate(order)}
    mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})

    fig = plt.figure(figsize=(8.4, 5.8), facecolor="white")
    main_ax = fig.add_axes([0.06, 0.18, 0.48, 0.70])
    image = main_ax.imshow(
        correlation.to_numpy(float), cmap=FIGURE2_DIVERGING, vmin=-1, vmax=1,
        interpolation="nearest", rasterized=True,
    )
    main_ax.set_xticks([])
    main_ax.set_yticks([])
    main_ax.set_title("297 target genes plus\npooled non-targeting control", fontsize=10, pad=8)
    for spine in main_ax.spines.values():
        spine.set_edgecolor("#A0A0A0")
        spine.set_linewidth(0.6)

    inset_positions = [
        [0.61, 0.62, 0.16, 0.22],
        [0.82, 0.65, 0.11, 0.16],
        [0.60, 0.23, 0.21, 0.25],
        [0.83, 0.27, 0.16, 0.21],
    ]
    inset_records: list[dict[str, Any]] = []
    box_records: list[dict[str, Any]] = []

    n = len(order)
    for (group_name, genes), position in zip(FIGURE2_TARGET_GROUPS.items(), inset_positions):
        group_ranks = sorted(ranks[gene] for gene in genes)
        if group_ranks != list(range(group_ranks[0], group_ranks[-1] + 1)):
            raise ValueError(f"Highlighted group is not contiguous after clustering: {group_name}")
        start, end = group_ranks[0], group_ranks[-1]
        rect = Rectangle((start - 0.5, start - 0.5), end - start + 1, end - start + 1,
                         fill=False, edgecolor="#303030", linewidth=0.8)
        main_ax.add_patch(rect)
        box_records.append({"group": group_name, "start_rank": start, "end_rank": end})

        ax = fig.add_axes(position)
        subset = correlation.loc[genes, genes]
        ax.imshow(subset.to_numpy(float), cmap=FIGURE2_DIVERGING, vmin=-1, vmax=1, interpolation="nearest")
        ax.set_xticks(np.arange(len(genes)), genes, rotation=90, fontsize=5.3)
        ax.set_yticks(np.arange(len(genes)), genes, fontsize=5.3)
        ax.set_title(group_name, fontsize=7.2, pad=3)
        for spine in ax.spines.values():
            spine.set_edgecolor("#707070")
            spine.set_linewidth(0.5)
        for row_gene in genes:
            for col_gene in genes:
                inset_records.append({
                    "group": group_name,
                    "row_gene": row_gene,
                    "column_gene": col_gene,
                    "pearson_r": float(subset.loc[row_gene, col_gene]),
                })

        y_center = (start + end) / 2
        connection = ConnectionPatch(
            xyA=(end + 0.5, y_center), coordsA=main_ax.transData,
            xyB=(-0.08, 0.5), coordsB=ax.transAxes,
            color="#505050", linewidth=0.55,
        )
        fig.add_artist(connection)

    cax = fig.add_axes([0.12, 0.09, 0.30, 0.025])
    cb = fig.colorbar(image, cax=cax, orientation="horizontal")
    cb.set_ticks([-1, 0, 1])
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("Pearson r — effect on program expression", fontsize=8)
    fig.text(0.015, 0.94, "a", fontsize=15, fontweight="bold", ha="left", va="top")
    return fig, correlation, pd.DataFrame({"cluster_rank": np.arange(1, n + 1), "target_name": order}), pd.DataFrame(inset_records), pd.DataFrame(box_records)

def plot_figure2b(selected_associations: pd.DataFrame, ec_stats: pd.DataFrame):
    targets = [gene for genes in FIGURE2_TARGET_GROUPS.values() for gene in genes]
    program_keys = [(item["day"], item["program"]) for item in FIGURE2_PROGRAMS]
    program_labels = [item["label"] for item in FIGURE2_PROGRAMS]
    lookup = selected_associations.set_index(["target_name", "Sample", "program_name"])
    values = np.full((len(targets), len(program_keys) + 1), np.nan)
    qvalues = np.full_like(values, np.nan)
    source_rows: list[dict[str, Any]] = []
    for row_index, target in enumerate(targets):
        for col_index, (day, program) in enumerate(program_keys):
            row = lookup.loc[(target, day, program)]
            effect = float(row["approx_log2FC"])
            qvalue = float(row["adj_pval"])
            values[row_index, col_index] = effect
            qvalues[row_index, col_index] = qvalue
            source_rows.append({
                "target": target, "day": day, "program": program,
                "program_label": program_labels[col_index],
                "effect": effect, "FDR": qvalue,
                "significance": significance_stars(qvalue),
            })

    ec_lookup = ec_stats.set_index("gene")
    for row_index, target in enumerate(targets):
        if target not in ec_lookup.index:
            raise KeyError(f"Missing EC-output statistics for {target}.")
        values[row_index, -1] = float(ec_lookup.loc[target, "effect_log2_target_over_control"])
        qvalues[row_index, -1] = float(ec_lookup.loc[target, "FDR_BH"])

    fig, ax = plt.subplots(figsize=(10.6, 6.0))
    norm = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
    image = ax.imshow(values, cmap=FIGURE2_DIVERGING, norm=norm, aspect="auto", interpolation="nearest")
    xlabels = [*program_labels, "D3 EC\nlog2(target / NTC)"]
    ax.set_xticks(np.arange(len(xlabels)), xlabels, rotation=55, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(targets)), targets, fontsize=8)
    ax.tick_params(length=0)
    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            star = significance_stars(qvalues[row, col])
            if star:
                color = "white" if abs(values[row, col]) >= 1.45 else "black"
                ax.text(col, row, star, ha="center", va="center", fontsize=7, fontweight="bold", color=color)

    row_start = 0
    for group_name, genes in FIGURE2_TARGET_GROUPS.items():
        row_end = row_start + len(genes)
        center = (row_start + row_end - 1) / 2
        ax.text(-2.55, center, group_name, ha="right", va="center", fontsize=8, fontweight="bold")
        if row_end < len(targets):
            ax.axhline(row_end - 0.5, color="white", linewidth=2.2)
        row_start = row_end

    spans = [(0, 2, "Day 1"), (3, 5, "Day 2"), (6, 9, "Day 3 — EC"), (10, 12, "Non-EC"), (13, 13, "EC\ndifferentiation")]
    for left, right, label in spans:
        center = (left + right) / 2
        ax.text(center, -1.35, label, ha="center", va="bottom", fontsize=8.5, fontweight="bold")
        if right < values.shape[1] - 1:
            ax.axvline(right + 0.5, color="white", linewidth=2.2)
    ax.set_xlim(-0.5, values.shape[1] - 0.5)
    ax.set_ylim(values.shape[0] - 0.5, -1.8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(image, ax=ax, fraction=0.028, pad=0.035)
    cb.set_label("Effect / EC ratio (log2 scale)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.text(0.012, 0.975, "b", fontsize=15, fontweight="bold", ha="left", va="top")
    fig.tight_layout(rect=[0.25, 0.03, 0.95, 0.95])

    effect_table = pd.DataFrame(values, index=targets, columns=[f"{d}_P{p}" for d, p in program_keys] + ["D3_EC_log2_target_over_NTC"])
    q_table = pd.DataFrame(qvalues, index=targets, columns=effect_table.columns)
    ec_source = ec_stats.copy()
    return fig, effect_table, q_table, pd.DataFrame(source_rows), ec_source
