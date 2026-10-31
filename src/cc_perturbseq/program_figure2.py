from __future__ import annotations

from collections.abc import Mapping, Sequence

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from .program_io import (
    effect_matrix,
    infer_column,
    read_table,
    significance_stars,
)


def _cluster_order(matrix: pd.DataFrame) -> list[str]:
    matrix = matrix.replace([np.inf, -np.inf], np.nan).fillna(0)
    if len(matrix) < 3:
        return matrix.index.astype(str).tolist()
    distance = pdist(matrix.to_numpy(float), metric="correlation")
    if not np.isfinite(distance).all() or np.allclose(distance, 0):
        return matrix.index.astype(str).tolist()
    order = leaves_list(linkage(distance, method="average"))
    return matrix.index[order].astype(str).tolist()


def plot_target_effect_correlation(
    effects: pd.DataFrame,
    *,
    target_groups: Mapping[str, Sequence[str]] | None = None,
    min_features: int = 10,
):
    matrix = effect_matrix(effects)
    matrix = matrix.loc[matrix.notna().sum(axis=1) >= min_features]
    order = _cluster_order(matrix)
    correlation = matrix.T.corr(min_periods=min_features).reindex(
        index=order, columns=order
    )
    target_groups = target_groups or {
        "BMP/TGFβ": ["BMPR1A", "SMAD1", "SMAD5", "SMAD4"],
        "Mesoderm": ["MIXL1", "EOMES"],
        "LM + EC": ["KDR", "KCTD10", "PRKCD", "ETV2", "PLCG1"],
        "Hippo": ["TAOK1", "AMOTL2", "LATS2", "NF2"],
    }

    figure = plt.figure(figsize=(8.8, 6), constrained_layout=True)
    grid = figure.add_gridspec(2, 3, width_ratios=[2.7, 1, 1])
    axis = figure.add_subplot(grid[:, 0])
    image = axis.imshow(
        correlation,
        cmap="RdBu_r",
        vmin=-1,
        vmax=1,
        aspect="auto",
        interpolation="nearest",
    )
    axis.set_title("Target perturbation-effect correlation", fontsize=9)
    axis.set_xticks([])
    axis.set_yticks([])
    figure.colorbar(image, ax=axis, fraction=0.03, pad=0.02, label="Pearson r")

    zoom_axes = [
        figure.add_subplot(grid[0, 1]),
        figure.add_subplot(grid[0, 2]),
        figure.add_subplot(grid[1, 1]),
        figure.add_subplot(grid[1, 2]),
    ]
    for zoom_axis, (label, targets) in zip(zoom_axes, target_groups.items()):
        present = [target for target in targets if target in correlation.index]
        zoom_axis.set_title(label, fontsize=8)
        if not present:
            zoom_axis.text(
                0.5,
                0.5,
                "No configured targets found",
                ha="center",
                va="center",
                transform=zoom_axis.transAxes,
                fontsize=6,
            )
            zoom_axis.set_axis_off()
            continue
        subset = correlation.loc[present, present]
        zoom_axis.imshow(subset, cmap="RdBu_r", vmin=-1, vmax=1)
        zoom_axis.set_xticks(
            np.arange(len(present)), present, rotation=90, fontsize=6
        )
        zoom_axis.set_yticks(np.arange(len(present)), present, fontsize=6)
        for row in range(len(present)):
            for column in range(len(present)):
                zoom_axis.text(
                    column,
                    row,
                    f"{subset.iloc[row, column]:.2f}",
                    ha="center",
                    va="center",
                    fontsize=5,
                )
    return figure, correlation, order


def load_ec_output(path: str) -> pd.Series:
    table = read_table(path)
    target_column = infer_column(table, ["target", "gene", "gene_symbol"])
    effect_column = infer_column(
        table,
        [
            "log2_ec_ratio",
            "ec% vs negative control",
            "ec_effect",
            "ec differentiation effect",
            "log2 effect",
            "effect",
        ],
    )
    return pd.to_numeric(table[effect_column], errors="coerce").groupby(
        table[target_column].astype(str)
    ).mean()


def plot_selected_program_effects(
    effects: pd.DataFrame,
    ec_output: pd.Series,
    *,
    target_groups: Mapping[str, Sequence[str]],
    programs_by_day: Mapping[str, Sequence[int]],
):
    targets = [target for group in target_groups.values() for target in group]
    available = set(zip(effects.day, effects.program))
    features = [
        (day, int(program))
        for day, programs in programs_by_day.items()
        for program in programs
        if (day, int(program)) in available
    ]
    subset = effects[effects.target.isin(targets)].copy()
    subset["feature"] = subset.apply(
        lambda row: f"{row.day}\nP{int(row.program)}", axis=1
    )
    columns = [f"{day}\nP{program}" for day, program in features]
    values = subset.pivot_table(
        index="target", columns="feature", values="effect", aggfunc="mean"
    ).reindex(index=targets, columns=columns)
    qvalues = subset.pivot_table(
        index="target", columns="feature", values="q_value", aggfunc="min"
    ).reindex(index=targets, columns=columns)
    values["D3 EC\noutput"] = ec_output.reindex(targets)

    figure, axis = plt.subplots(
        figsize=(max(7.5, 0.42 * values.shape[1]), max(4.5, 0.32 * values.shape[0]))
    )
    finite = np.abs(values.to_numpy(float))
    finite = finite[np.isfinite(finite)]
    vmax = max(0.5, float(np.percentile(finite, 98))) if finite.size else 1.0
    image = axis.imshow(
        values,
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax),
        aspect="auto",
    )
    axis.set_xticks(
        np.arange(values.shape[1]), values.columns, rotation=90, fontsize=7
    )
    axis.set_yticks(np.arange(values.shape[0]), values.index, fontsize=7)
    figure.colorbar(
        image, ax=axis, fraction=0.025, pad=0.02, label="Effect / EC ratio (log scale)"
    )
    for row, target in enumerate(values.index):
        for column, feature in enumerate(columns):
            star = significance_stars(qvalues.loc[target, feature])
            if star:
                axis.text(
                    column, row, star, ha="center", va="center", fontsize=6
                )

    start = 0
    for label, group in target_groups.items():
        end = start + len(group)
        axis.text(
            -0.75,
            (start + end - 1) / 2,
            label,
            ha="right",
            va="center",
            fontsize=7,
        )
        if end < len(targets):
            axis.axhline(end - 0.5, color="white", linewidth=2)
        start = end
    return figure, values, qvalues
