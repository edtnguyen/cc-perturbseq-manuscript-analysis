from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage

from .program_io import clean_day, infer_column, read_table, significance_stars


def load_multiome_qc(path: str | Path) -> dict[str, pd.DataFrame]:
    path = Path(path)
    if path.is_file():
        return {"qc": read_table(path)}
    output: dict[str, pd.DataFrame] = {}
    for file in path.rglob("*"):
        if file.suffix.lower() not in {".csv", ".tsv", ".txt", ".parquet", ".pq"}:
            continue
        name = file.stem.lower()
        if "fragment" in name:
            output["fragments"] = read_table(file)
        elif "cell" in name:
            output["cells"] = read_table(file)
        elif "tss" in name:
            output["tss"] = read_table(file)
    return output


def plot_multiome_qc(qc: dict[str, pd.DataFrame], tss_path: str | Path | None = None):
    fragments = qc.get("fragments", qc.get("qc"))
    cells = qc.get("cells", qc.get("qc"))
    tss = read_table(tss_path) if tss_path is not None and Path(tss_path).exists() else qc.get("tss")
    if fragments is None or cells is None or tss is None:
        raise KeyError("Multiome QC requires fragment, cell-count and TSS-profile tables")

    day_f = infer_column(fragments, ["day", "timepoint"])
    fragment_col = infer_column(fragments, ["atac_fragments", "fragments", "n_fragments", "total_fragments"])
    day_c = infer_column(cells, ["day", "timepoint"])
    cell_col = infer_column(cells, ["cells", "n_cells", "cell_count", "retained_cells"])
    day_t = infer_column(tss, ["day", "timepoint"])
    distance_col = infer_column(tss, ["distance", "distance_from_tss", "position"])
    enrichment_col = infer_column(tss, ["tss_enrichment", "enrichment", "normalized_insertions"])

    fragments = fragments.assign(day_clean=fragments[day_f].map(clean_day))
    cells = cells.assign(day_clean=cells[day_c].map(clean_day))
    tss = tss.assign(day_clean=tss[day_t].map(clean_day))
    order = ["D0", "D1", "D2", "D3", "D4"]
    fragment_summary = fragments.groupby("day_clean")[fragment_col].sum().reindex(order).dropna()
    cell_summary = cells.groupby("day_clean")[cell_col].sum().reindex(order).dropna()

    fig, axes = plt.subplots(1, 3, figsize=(9.4, 2.9), constrained_layout=True)
    axes[0].bar(fragment_summary.index, fragment_summary.values / 1e6)
    axes[0].set_ylabel("ATAC fragments (millions)")
    axes[1].bar(cell_summary.index, cell_summary.values)
    axes[1].set_ylabel("RNA-QC cells")
    for day, subset in tss.groupby("day_clean", sort=False):
        subset = subset.sort_values(distance_col)
        axes[2].plot(subset[distance_col], subset[enrichment_col], linewidth=1, label=day)
    axes[2].set_xlabel("Distance from TSS (bp)")
    axes[2].set_ylabel("TSS enrichment")
    axes[2].legend(frameon=False, fontsize=6)
    for axis in axes:
        axis.tick_params(labelsize=6)
    return fig, {"fragments": fragment_summary.to_frame("fragments"), "cells": cell_summary.to_frame("cells"), "tss": tss}


def standardize_chrombpnet_performance(path: str | Path) -> pd.DataFrame:
    table = read_table(path)
    day_col = infer_column(table, ["day", "timepoint"])
    fold_col = infer_column(table, ["fold", "cv_fold"], required=False)
    observed_col = infer_column(table, ["observed", "observed_logcounts", "observed_counts"], required=False)
    predicted_col = infer_column(table, ["predicted", "predicted_logcounts", "predicted_counts"], required=False)
    count_col = infer_column(table, ["count_pearson", "pearson_r", "count_correlation"], required=False)
    jsd_col = infer_column(table, ["profile_jsd", "jsd", "median_jsd"], required=False)
    output = pd.DataFrame(
        {
            "day": table[day_col].map(clean_day),
            "fold": table[fold_col].astype(str) if fold_col else "pooled",
            "observed": pd.to_numeric(table[observed_col], errors="coerce") if observed_col else np.nan,
            "predicted": pd.to_numeric(table[predicted_col], errors="coerce") if predicted_col else np.nan,
            "count_pearson": pd.to_numeric(table[count_col], errors="coerce") if count_col else np.nan,
            "profile_jsd": pd.to_numeric(table[jsd_col], errors="coerce") if jsd_col else np.nan,
        }
    )
    return output


def plot_chrombpnet_performance(table: pd.DataFrame):
    days = [day for day in ("D0", "D1", "D2", "D3", "D4") if day in set(table.day)]
    fig_scatter, axes = plt.subplots(1, len(days), figsize=(2.2 * len(days), 2.4), squeeze=False, constrained_layout=True)
    scatter_source: list[pd.DataFrame] = []
    for axis, day in zip(axes.ravel(), days):
        subset = table[(table.day == day) & table.observed.notna() & table.predicted.notna()]
        if len(subset) > 100000:
            subset = subset.sample(100000, random_state=7)
        axis.hexbin(subset.observed, subset.predicted, gridsize=45, mincnt=1, cmap="viridis", rasterized=True)
        r = subset[["observed", "predicted"]].corr().iloc[0, 1] if len(subset) > 1 else np.nan
        axis.set_title(f"{day}\nr = {r:.3f}", fontsize=7)
        axis.set_xlabel("Observed log(counts+1)", fontsize=6)
        axis.set_ylabel("Predicted log(counts+1)", fontsize=6)
        axis.tick_params(labelsize=5.5)
        scatter_source.append(subset.assign(day=day))

    summary = table.groupby(["day", "fold"], as_index=False)[["count_pearson", "profile_jsd"]].mean()
    fig_summary, axes2 = plt.subplots(1, 2, figsize=(5.2, 2.7), constrained_layout=True)
    for axis, metric, label in zip(axes2, ["count_pearson", "profile_jsd"], ["Count prediction Pearson r", "Profile prediction median JSD"]):
        for index, day in enumerate(days):
            values = summary.loc[summary.day == day, metric].dropna().to_numpy()
            if len(values):
                axis.scatter(np.full(len(values), index), values, s=12)
                axis.bar(index, values.mean(), alpha=0.25)
        axis.set_xticks(range(len(days)), days)
        axis.set_ylabel(label, fontsize=7)
        axis.tick_params(labelsize=6)
    return fig_scatter, fig_summary, pd.concat(scatter_source, ignore_index=True), summary


def standardize_motif_qc(path: str | Path) -> pd.DataFrame:
    table = read_table(path)
    day_col = infer_column(table, ["day", "timepoint"])
    seqlet_col = infer_column(table, ["seqlets", "n_seqlets", "modisco_seqlets"])
    instance_col = infer_column(table, ["motif_instances", "finemo_instances", "instances"])
    motif_col = infer_column(table, ["motifs", "n_motifs", "motif_count"])
    sign_col = infer_column(table, ["sign", "pattern_sign", "positive_negative"], required=False)
    return pd.DataFrame(
        {
            "day": table[day_col].map(clean_day),
            "seqlets": pd.to_numeric(table[seqlet_col], errors="coerce"),
            "instances": pd.to_numeric(table[instance_col], errors="coerce"),
            "motifs": pd.to_numeric(table[motif_col], errors="coerce"),
            "sign": table[sign_col].astype(str) if sign_col else "all",
        }
    )


def plot_motif_qc(table: pd.DataFrame):
    order = [day for day in ("D0", "D1", "D2", "D3", "D4") if day in set(table.day)]
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.8), constrained_layout=True)
    seqlets = table.groupby("day").seqlets.sum().reindex(order)
    axes[0].bar(order, seqlets)
    axes[0].set_ylabel("TF-MoDISco seqlets")
    pivot = table.pivot_table(index="day", columns="sign", values="motifs", aggfunc="sum").reindex(order)
    bottom = np.zeros(len(order))
    for sign in pivot.columns:
        values = pivot[sign].fillna(0).to_numpy()
        axes[1].bar(order, values, bottom=bottom, label=sign)
        bottom += values
    axes[1].set_ylabel("QC-passing motifs")
    axes[1].legend(frameon=False, fontsize=6)
    for axis in axes:
        axis.tick_params(labelsize=6)
    return fig


def standardize_enrichment(path: str | Path) -> pd.DataFrame:
    table = read_table(path)
    day_col = infer_column(table, ["day", "timepoint"])
    program_col = infer_column(table, ["program", "program_id"])
    tf_col = infer_column(table, ["tf", "motif", "tf_motif_id"])
    q_col = infer_column(table, ["qval", "q_value", "qvalue", "FDR"])
    value_col = infer_column(table, ["neg_log10_qval", "minus_log10_q", "score"], required=False)
    auc_col = infer_column(table, ["AUC", "auc", "rank_biserial"], required=False)
    output = pd.DataFrame(
        {
            "day": table[day_col].map(clean_day),
            "program": pd.to_numeric(table[program_col], errors="coerce"),
            "tf": table[tf_col].astype(str),
            "q_value": pd.to_numeric(table[q_col], errors="coerce"),
            "value": pd.to_numeric(table[value_col], errors="coerce") if value_col else -np.log10(pd.to_numeric(table[q_col], errors="coerce").clip(lower=np.nextafter(0, 1))),
            "auc": pd.to_numeric(table[auc_col], errors="coerce") if auc_col else np.nan,
        }
    )
    return output.dropna(subset=["day", "program", "tf", "q_value"])


def _motif_heatmap(table: pd.DataFrame, day: str, programs: list[int]):
    subset = table[(table.day == day) & table.program.isin(programs)].copy()
    significant = subset.loc[subset.q_value < 0.05, "tf"].unique().tolist()
    values = subset.pivot_table(index="tf", columns="program", values="value", aggfunc="max").reindex(index=significant, columns=programs).fillna(0)
    qvalues = subset.pivot_table(index="tf", columns="program", values="q_value", aggfunc="min").reindex(index=significant, columns=programs).fillna(1)
    if len(values) > 1:
        order = leaves_list(linkage(values.to_numpy(), method="ward", metric="euclidean"))
        values = values.iloc[order]
        qvalues = qvalues.iloc[order]
    return values, qvalues


def plot_selected_motif_enrichment(table: pd.DataFrame):
    panels = [("D2", [15, 17]), ("D3", [2, 24, 27, 34])]
    matrices = [_motif_heatmap(table, day, programs) for day, programs in panels]
    vmax_values = np.concatenate([matrix.to_numpy().ravel() for matrix, _ in matrices if matrix.size])
    vmax = max(1.0, float(np.nanpercentile(vmax_values[vmax_values > 0], 99))) if np.any(vmax_values > 0) else 1.0
    fig, axes = plt.subplots(1, 2, figsize=(8.8, max(4.2, 0.22 * max(len(matrix) for matrix, _ in matrices))), constrained_layout=True)
    source: dict[str, pd.DataFrame] = {}
    for axis, (day, programs), (values, qvalues) in zip(axes, panels, matrices):
        image = axis.imshow(values, cmap="Blues", norm=Normalize(0, vmax), aspect="auto")
        axis.set_xticks(np.arange(len(programs)), [f"P{program}" for program in programs], fontsize=7)
        axis.set_yticks(np.arange(len(values)), [re.sub(r"_.*$", "", value) for value in values.index], fontsize=5.5)
        axis.set_title(day, fontsize=8)
        for row, tf in enumerate(values.index):
            for column, program in enumerate(programs):
                star = significance_stars(qvalues.loc[tf, program])
                if star:
                    axis.text(column, row, star, ha="center", va="center", fontsize=5)
        source[f"{day}_values"] = values
        source[f"{day}_qvalues"] = qvalues
    fig.colorbar(image, ax=axes, fraction=0.025, pad=0.02, label="−log₁₀(q-value)")
    return fig, source
