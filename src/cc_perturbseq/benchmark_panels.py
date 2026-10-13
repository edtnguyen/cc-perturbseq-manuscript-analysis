from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from .program_io import infer_column, read_table


def _load_adata(path: str | Path):
    import anndata as ad

    return ad.read_h5ad(path)


def _numeric_obs(adata, candidates: list[str]) -> pd.Series | None:
    for column in candidates:
        if column in adata.obs.columns:
            values = pd.to_numeric(adata.obs[column], errors="coerce")
            if values.notna().any():
                return values
    return None


def _condition(adata) -> pd.Series:
    for column in ("capture_method", "method", "condition", "protocol", "capture", "primer_condition"):
        if column in adata.obs.columns:
            return adata.obs[column].astype(str)
    raise KeyError("Could not infer benchmark condition column from adata.obs")


def _guide_matrix(adata):
    for key in ("guide_umi_counts", "guide_counts", "guide_umis", "guide_capture_counts"):
        if key in adata.obsm:
            return adata.obsm[key]
    guide_positions = [index for index, name in enumerate(adata.var_names.astype(str)) if re_search_guide(name)]
    if guide_positions:
        return adata[:, guide_positions].X
    return None


def re_search_guide(name: str) -> bool:
    text = str(name).lower()
    return text.startswith("sg_") or "guide" in text or "grna" in text


def benchmark_cell_metrics(path: str | Path) -> pd.DataFrame:
    adata = _load_adata(path)
    condition = _condition(adata)
    transcriptome = _numeric_obs(adata, ["total_counts", "transcriptome_umis", "gex_umis", "nCount_RNA"])
    if transcriptome is None:
        transcriptome = pd.Series(np.asarray(adata.X.sum(axis=1)).ravel(), index=adata.obs_names)
    matrix = _guide_matrix(adata)
    if matrix is None:
        top_guide = _numeric_obs(adata, ["top_guide_umis", "max_guide_umis", "most_abundant_guide_umis"])
        total_guide = _numeric_obs(adata, ["guide_umi_counts", "total_guide_umis"])
    else:
        if sparse.issparse(matrix):
            top_guide = pd.Series(matrix.max(axis=1).toarray().ravel(), index=adata.obs_names)
            total_guide = pd.Series(np.asarray(matrix.sum(axis=1)).ravel(), index=adata.obs_names)
        else:
            array = np.asarray(matrix)
            top_guide = pd.Series(array.max(axis=1), index=adata.obs_names)
            total_guide = pd.Series(array.sum(axis=1), index=adata.obs_names)
    output = pd.DataFrame(
        {
            "condition": condition.to_numpy(),
            "transcriptome_umis": transcriptome.to_numpy(),
            "top_guide_umis": top_guide.to_numpy() if top_guide is not None else np.nan,
            "total_guide_umis": total_guide.to_numpy() if total_guide is not None else np.nan,
        },
        index=adata.obs_names,
    )
    for column in ("cell_type", "assigned_cell_type", "lineage", "primer_condition"):
        if column in adata.obs.columns:
            output[column] = adata.obs[column].astype(str).to_numpy()
    return output


def _violin(axis, table: pd.DataFrame, value: str, order: list[str], ylabel: str):
    vectors = [pd.to_numeric(table.loc[table.condition == condition, value], errors="coerce").dropna() for condition in order]
    axis.violinplot(vectors, positions=np.arange(len(order)), showextrema=False, showmedians=True)
    axis.set_xticks(np.arange(len(order)), order, rotation=30, ha="right", fontsize=6)
    axis.set_ylabel(ylabel, fontsize=7)
    axis.tick_params(axis="y", labelsize=6)
    for index, vector in enumerate(vectors):
        if len(vector):
            axis.text(index, np.nanmax(vector), f"{np.nanmedian(vector):.1f}", ha="center", va="bottom", fontsize=5.5)


def plot_capture_benchmark(cell_metrics: pd.DataFrame, *, order: list[str] | None = None):
    order = order or list(pd.unique(cell_metrics.condition))
    fig, axes = plt.subplots(1, 3, figsize=(8.2, 2.8), constrained_layout=True)
    _violin(axes[0], cell_metrics, "top_guide_umis", order, "Guide UMIs per cell")
    detection = cell_metrics.assign(detected=cell_metrics.top_guide_umis > 0).groupby("condition").detected.mean().reindex(order) * 100
    axes[1].bar(np.arange(len(order)), detection.values)
    axes[1].set_xticks(np.arange(len(order)), order, rotation=30, ha="right", fontsize=6)
    axes[1].set_ylabel("Cells with ≥1 guide (%)", fontsize=7)
    for index, value in enumerate(detection.values):
        axes[1].text(index, value, f"{value:.1f}", ha="center", va="bottom", fontsize=5.5)
    _violin(axes[2], cell_metrics, "transcriptome_umis", order, "Transcriptome UMIs per cell")
    return fig, detection.rename("detection_percent").to_frame()


def standardize_downsampling(path: str | Path) -> pd.DataFrame:
    table = read_table(path)
    condition_col = infer_column(table, ["condition", "method", "capture_method", "primer"])
    reads_col = infer_column(table, ["assigned_reads_per_cell", "mean_reads_per_cell", "reads_per_cell", "reads"])
    transcript_col = infer_column(table, ["median_transcript_umis", "median_umis", "transcript_umis"], required=False)
    gene_col = infer_column(table, ["median_genes", "detected_genes", "genes_per_cell"], required=False)
    return pd.DataFrame(
        {
            "condition": table[condition_col].astype(str),
            "reads_per_cell": pd.to_numeric(table[reads_col], errors="coerce"),
            "transcript_umis": pd.to_numeric(table[transcript_col], errors="coerce") if transcript_col else np.nan,
            "detected_genes": pd.to_numeric(table[gene_col], errors="coerce") if gene_col else np.nan,
        }
    )


def plot_downsampling(table: pd.DataFrame, *, outcome: str, ylabel: str):
    fig, ax = plt.subplots(figsize=(4.1, 3.2))
    for condition, subset in table.groupby("condition", sort=False):
        subset = subset.sort_values("reads_per_cell")
        ax.plot(subset.reads_per_cell, subset[outcome], marker="o", markersize=2.5, linewidth=1, label=condition)
    ax.set_xlabel("Assigned reads per cell")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False, fontsize=6)
    return fig


def standardize_pilot_knockdown(path: str | Path) -> pd.DataFrame:
    table = read_table(path)
    target_col = infer_column(table, ["target", "gene", "target_gene"])
    condition_col = infer_column(table, ["condition", "method", "capture_method"])
    expression_col = infer_column(table, ["normalized_expression", "mean_expression", "expression", "relative_expression"])
    q_col = infer_column(table, ["q_value", "qvalue", "FDR", "padj"], required=False)
    lower_col = infer_column(table, ["lower", "ci_lower", "lower_ci"], required=False)
    upper_col = infer_column(table, ["upper", "ci_upper", "upper_ci"], required=False)
    return pd.DataFrame(
        {
            "target": table[target_col].astype(str),
            "condition": table[condition_col].astype(str),
            "expression": pd.to_numeric(table[expression_col], errors="coerce"),
            "q_value": pd.to_numeric(table[q_col], errors="coerce") if q_col else np.nan,
            "lower": pd.to_numeric(table[lower_col], errors="coerce") if lower_col else np.nan,
            "upper": pd.to_numeric(table[upper_col], errors="coerce") if upper_col else np.nan,
        }
    )


def plot_pilot_knockdown(table: pd.DataFrame):
    targets = table.groupby("target").expression.max().sort_values().index.tolist()
    conditions = list(pd.unique(table.condition))
    fig, ax = plt.subplots(figsize=(max(7, 0.35 * len(targets)), 3.3))
    offsets = np.linspace(-0.18, 0.18, len(conditions))
    for offset, condition in zip(offsets, conditions):
        subset = table[table.condition == condition].set_index("target").reindex(targets)
        x = np.arange(len(targets)) + offset
        ax.errorbar(
            x,
            subset.expression,
            yerr=[subset.expression - subset.lower, subset.upper - subset.expression] if subset.lower.notna().any() else None,
            fmt="o",
            markersize=3,
            linewidth=0.8,
            label=condition,
        )
    ax.set_xticks(np.arange(len(targets)), targets, rotation=60, ha="right", fontsize=6)
    ax.set_ylabel("Normalized target-gene expression")
    ax.legend(frameon=False, fontsize=6)
    return fig


def standardize_parse_metrics(
    adata_path: str | Path,
    guide_counts_path: str | Path | None = None,
) -> pd.DataFrame:
    metrics = benchmark_cell_metrics(adata_path)
    adata = _load_adata(adata_path)
    primer = next((column for column in ("primer_condition", "condition", "primer", "rt_primer") if column in adata.obs.columns), None)
    cell_type = next((column for column in ("cell_type", "assigned_cell_type", "lineage") if column in adata.obs.columns), None)
    if primer:
        metrics["primer"] = adata.obs[primer].astype(str).to_numpy()
    else:
        metrics["primer"] = metrics.condition
    if cell_type:
        metrics["cell_type"] = adata.obs[cell_type].astype(str).to_numpy()
    if guide_counts_path is not None and Path(guide_counts_path).exists():
        counts = read_table(guide_counts_path)
        cell_col = infer_column(counts, ["cell", "barcode", "cell_barcode"])
        hi_col = infer_column(counts, ["hipsc_pool_umis", "hipsc_umis", "wTC11_umis"])
        k_col = infer_column(counts, ["k562_pool_umis", "k562_umis"])
        pool = pd.DataFrame(
            {
                "hipsc_pool_umis": pd.to_numeric(counts[hi_col], errors="coerce"),
                "k562_pool_umis": pd.to_numeric(counts[k_col], errors="coerce"),
            },
            index=counts[cell_col].astype(str),
        )
        metrics = metrics.join(pool, how="left")
        matched = np.where(metrics.cell_type.str.contains("K562", case=False, na=False), metrics.k562_pool_umis, metrics.hipsc_pool_umis)
        cross = np.where(metrics.cell_type.str.contains("K562", case=False, na=False), metrics.hipsc_pool_umis, metrics.k562_pool_umis)
        metrics["matched_pool_umis"] = matched
        metrics["cross_pool_umis"] = cross
        metrics["guide_pool_specificity"] = matched / (matched + cross).replace(0, np.nan)
    return metrics


def plot_parse_benchmark(metrics: pd.DataFrame):
    primers = list(pd.unique(metrics.primer))
    cell_types = list(pd.unique(metrics.cell_type.dropna())) if "cell_type" in metrics else ["All"]
    fig, axes = plt.subplots(2, 3, figsize=(9.4, 5.5), constrained_layout=True)
    for cell_type in cell_types[:2]:
        subset = metrics[metrics.cell_type == cell_type]
        for primer, part in subset.groupby("primer", sort=False):
            x = primers.index(primer)
            values = part.get("matched_pool_umis", part.top_guide_umis).dropna()
            axes[0, 0].scatter(np.full(min(len(values), 3000), x), values.iloc[:3000], s=1, alpha=0.15, rasterized=True)
        axes[0, 0].set_title("Matched-pool guide UMIs", fontsize=8)
    if "guide_pool_specificity" in metrics:
        data = [metrics.loc[metrics.primer == primer, "guide_pool_specificity"].dropna() for primer in primers]
        axes[0, 1].violinplot(data, positions=np.arange(len(primers)), showextrema=False, showmedians=True)
        axes[0, 1].set_title("Guide-pool specificity", fontsize=8)
    data = [metrics.loc[metrics.primer == primer, "top_guide_umis"].dropna() for primer in primers]
    axes[0, 2].violinplot(data, positions=np.arange(len(primers)), showextrema=False, showmedians=True)
    axes[0, 2].set_title("Top-guide UMIs", fontsize=8)
    data = [metrics.loc[metrics.primer == primer, "transcriptome_umis"].dropna() for primer in primers]
    axes[1, 0].violinplot(data, positions=np.arange(len(primers)), showextrema=False, showmedians=True)
    axes[1, 0].set_title("Transcriptome UMIs", fontsize=8)
    if {"hipsc_pool_umis", "k562_pool_umis"}.issubset(metrics.columns):
        axes[1, 1].hist2d(
            np.log10(metrics.hipsc_pool_umis.fillna(0) + 1),
            np.log10(metrics.k562_pool_umis.fillna(0) + 1),
            bins=60,
            cmap="viridis",
        )
        axes[1, 1].set_xlabel("log₁₀(hiPSC-pool UMIs + 1)", fontsize=6)
        axes[1, 1].set_ylabel("log₁₀(K562-pool UMIs + 1)", fontsize=6)
        axes[1, 1].set_title("Guide-pool mixing", fontsize=8)
    axes[1, 2].axis("off")
    for axis in axes.ravel()[:4]:
        axis.set_xticks(np.arange(len(primers)), primers, rotation=45, ha="right", fontsize=5.5)
        axis.tick_params(axis="y", labelsize=5.5)
    return fig
