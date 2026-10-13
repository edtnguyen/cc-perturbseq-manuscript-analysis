from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .program_io import clean_day, infer_column, read_table


def load_qc_bundle(path: str | Path) -> dict[str, Any]:
    """Load a standardized preprocessing-QC bundle.

    A directory may contain ``barcode_rank``, ``cell_metrics``, ``retention``
    and ``guide_metrics`` CSV/TSV/Parquet files. A h5ad input supplies
    post-QC cell and guide metrics but cannot reconstruct discarded barcodes.
    """
    path = Path(path)
    if path.suffix.lower() == ".h5ad":
        import anndata as ad

        adata = ad.read_h5ad(path, backed="r")
        return {"adata": adata, "cell_metrics": adata.obs.copy()}
    if path.is_file():
        return {"cell_metrics": read_table(path)}
    aliases = {
        "barcode_rank": ["barcode_rank", "barcode-rank", "rank"],
        "cell_metrics": ["cell_metrics", "qc_metrics", "per_cell"],
        "retention": ["retention", "cells_retained", "cell_counts"],
        "guide_metrics": ["guide_metrics", "guide_qc", "per_cell_guide"],
    }
    output: dict[str, Any] = {}
    files = [file for file in path.rglob("*") if file.suffix.lower() in {".csv", ".tsv", ".txt", ".parquet", ".pq"}]
    for key, patterns in aliases.items():
        match = next((file for file in files if any(pattern in file.stem.lower() for pattern in patterns)), None)
        if match:
            output[key] = read_table(match)
    if not output:
        raise FileNotFoundError(f"No recognized QC tables below {path}")
    return output


def plot_barcode_rank(table: pd.DataFrame):
    rank_col = infer_column(table, ["rank", "barcode_rank"])
    umi_col = infer_column(table, ["total_umis", "total_counts", "umi_count"])
    day_col = infer_column(table, ["day", "timepoint"])
    library_col = infer_column(table, ["sublibrary", "pool", "library", "sample"], required=False)
    call_col = infer_column(table, ["is_cell", "called_cell", "cell_call"], required=False)
    threshold_col = infer_column(table, ["threshold", "inflection", "cell_call_threshold"], required=False)
    data = table.copy()
    data["day_clean"] = data[day_col].map(clean_day)
    fig, axes = plt.subplots(2, 2, figsize=(6.8, 5.8), constrained_layout=True)
    for axis, day in zip(axes.ravel(), ["D0", "D1", "D2", "D3"]):
        subset = data[data.day_clean == day]
        groups = subset.groupby(library_col) if library_col else [("all", subset)]
        for label, group in groups:
            axis.plot(
                pd.to_numeric(group[rank_col], errors="coerce"),
                pd.to_numeric(group[umi_col], errors="coerce"),
                linewidth=0.8,
                label=str(label),
            )
        if call_col:
            called = subset[subset[call_col].astype(bool)]
            if len(called):
                axis.axvline(pd.to_numeric(called[rank_col], errors="coerce").max(), color="0.4", linestyle="--", linewidth=0.7)
        elif threshold_col:
            values = pd.to_numeric(subset[threshold_col], errors="coerce").dropna()
            if len(values):
                axis.axhline(values.iloc[0], color="0.4", linestyle="--", linewidth=0.7)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(day, fontsize=8)
        axis.set_xlabel("Barcode rank", fontsize=7)
        axis.set_ylabel("UMI count", fontsize=7)
        axis.tick_params(labelsize=6)
    return fig, data


def _cell_metrics(bundle: dict[str, Any]) -> pd.DataFrame:
    if "cell_metrics" not in bundle:
        raise KeyError("QC bundle lacks cell_metrics")
    data = bundle["cell_metrics"].copy()
    day_col = infer_column(data, ["day", "timepoint", "diff_day"])
    data["day_clean"] = data[day_col].map(clean_day)
    return data


def plot_umi_gene_scatter(bundle: dict[str, Any]):
    data = _cell_metrics(bundle)
    umi_col = infer_column(data, ["total_counts", "total_umis", "nCount_RNA"])
    gene_col = infer_column(data, ["n_genes_by_counts", "detected_genes", "nFeature_RNA"])
    barcode_col = infer_column(data, ["round1", "round_1", "well", "bci"], required=False)
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.8), constrained_layout=True)
    for axis, day in zip(axes.ravel(), ["D0", "D1", "D2", "D3"]):
        subset = data[data.day_clean == day]
        color = pd.Categorical(subset[barcode_col]).codes if barcode_col else "0.4"
        axis.scatter(
            pd.to_numeric(subset[umi_col], errors="coerce"),
            pd.to_numeric(subset[gene_col], errors="coerce"),
            c=color,
            cmap="viridis" if barcode_col else None,
            color=None if barcode_col else "0.4",
            s=2,
            alpha=0.35,
            linewidths=0,
            rasterized=True,
        )
        axis.set_title(day, fontsize=8)
        axis.set_xlabel("UMIs per cell", fontsize=7)
        axis.set_ylabel("Genes per cell", fontsize=7)
        axis.tick_params(labelsize=6)
    return fig, data


def plot_metric_distributions(bundle: dict[str, Any]):
    data = _cell_metrics(bundle)
    stage_col = infer_column(data, ["qc_stage", "stage", "before_after", "retained"], required=False)
    metrics = [
        ("UMIs", infer_column(data, ["total_counts", "total_umis"])),
        ("Genes", infer_column(data, ["n_genes_by_counts", "detected_genes"])),
        ("Mitochondrial (%)", infer_column(data, ["pct_counts_mt", "mitochondrial_fraction"])),
        ("Ribosomal (%)", infer_column(data, ["pct_counts_ribo", "ribosomal_fraction"])),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(11.2, 3.0), constrained_layout=True)
    day_order = ["D0", "D1", "D2", "D3"]
    for axis, (label, column) in zip(axes, metrics):
        positions = []
        values = []
        labels = []
        position = 0
        for day in day_order:
            day_data = data[data.day_clean == day]
            stages = day_data[stage_col].dropna().unique() if stage_col else ["retained"]
            for stage in stages:
                subset = day_data if not stage_col else day_data[day_data[stage_col] == stage]
                vector = pd.to_numeric(subset[column], errors="coerce").dropna().to_numpy()
                if len(vector):
                    positions.append(position)
                    values.append(vector)
                    labels.append(f"{day}\n{stage}" if stage_col else day)
                    position += 1
            position += 0.5
        axis.violinplot(values, positions=positions, showextrema=False, widths=0.8)
        axis.set_xticks(positions, labels, rotation=90, fontsize=5.5)
        axis.set_title(label, fontsize=8)
        axis.tick_params(axis="y", labelsize=6)
    return fig, data


def plot_cells_retained(bundle: dict[str, Any]):
    if "retention" in bundle:
        table = bundle["retention"].copy()
        day_col = infer_column(table, ["day", "timepoint"])
        count_col = infer_column(table, ["cells", "n_cells", "retained_cells", "count"])
        table["day_clean"] = table[day_col].map(clean_day)
        counts = table.groupby("day_clean")[count_col].sum().reindex(["D0", "D1", "D2", "D3"])
    else:
        data = _cell_metrics(bundle)
        counts = data.groupby("day_clean").size().reindex(["D0", "D1", "D2", "D3"])
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    ax.bar(counts.index, counts.values)
    for index, value in enumerate(counts.values):
        ax.text(index, value, f"{int(value):,}", ha="center", va="bottom", fontsize=6)
    ax.set_ylabel("Cells after gene-expression QC")
    return fig, counts.rename("cells").to_frame()


def plot_guide_metrics(bundle: dict[str, Any]):
    data = bundle.get("guide_metrics", _cell_metrics(bundle)).copy()
    day_col = infer_column(data, ["day", "timepoint", "diff_day"])
    data["day_clean"] = data[day_col].map(clean_day)
    guide_umi_col = infer_column(data, ["guide_umi_counts", "total_guide_umis", "guide_umis"])
    multiplicity_col = infer_column(data, ["guides_per_cell", "n_assigned_guides", "guide_count"])
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.0), constrained_layout=True)
    for axis, column, title in zip(
        axes,
        [guide_umi_col, multiplicity_col],
        ["Guide UMIs per cell", "Provisional guides per cell"],
    ):
        values = [pd.to_numeric(data.loc[data.day_clean == day, column], errors="coerce").dropna() for day in ["D0", "D1", "D2", "D3"]]
        axis.violinplot(values, positions=np.arange(4), showextrema=False)
        axis.set_xticks(np.arange(4), ["D0", "D1", "D2", "D3"])
        axis.set_title(title, fontsize=8)
        axis.tick_params(labelsize=6)
    return fig, data
