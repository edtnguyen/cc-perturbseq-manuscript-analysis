from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .benchmark_panels import _load_adata, benchmark_cell_metrics
from .program_io import infer_column, read_table


def standardize_parse_metrics(
    adata_path: str | Path,
    guide_counts_path: str | Path | None = None,
) -> pd.DataFrame:
    """Standardize Parse benchmark metrics with safe vectorized specificity."""
    metrics = benchmark_cell_metrics(adata_path)
    adata = _load_adata(adata_path)
    primer = next(
        (
            column
            for column in ("primer_condition", "condition", "primer", "rt_primer")
            if column in adata.obs.columns
        ),
        None,
    )
    cell_type = next(
        (
            column
            for column in ("cell_type", "assigned_cell_type", "lineage")
            if column in adata.obs.columns
        ),
        None,
    )
    metrics["primer"] = (
        adata.obs[primer].astype(str).to_numpy() if primer else metrics.condition
    )
    if cell_type:
        metrics["cell_type"] = adata.obs[cell_type].astype(str).to_numpy()

    if guide_counts_path is not None and Path(guide_counts_path).exists():
        counts = read_table(guide_counts_path)
        cell_col = infer_column(counts, ["cell", "barcode", "cell_barcode"])
        hi_col = infer_column(counts, ["hipsc_pool_umis", "hipsc_umis", "wtc11_umis"])
        k_col = infer_column(counts, ["k562_pool_umis", "k562_umis"])
        pool = pd.DataFrame(
            {
                "hipsc_pool_umis": pd.to_numeric(counts[hi_col], errors="coerce"),
                "k562_pool_umis": pd.to_numeric(counts[k_col], errors="coerce"),
            },
            index=counts[cell_col].astype(str),
        )
        metrics = metrics.join(pool, how="left")
        is_k562 = metrics.cell_type.str.contains("K562", case=False, na=False)
        matched = np.where(is_k562, metrics.k562_pool_umis, metrics.hipsc_pool_umis).astype(float)
        cross = np.where(is_k562, metrics.hipsc_pool_umis, metrics.k562_pool_umis).astype(float)
        denominator = matched + cross
        metrics["matched_pool_umis"] = matched
        metrics["cross_pool_umis"] = cross
        metrics["guide_pool_specificity"] = np.divide(
            matched,
            denominator,
            out=np.full_like(matched, np.nan, dtype=float),
            where=np.isfinite(denominator) & (denominator > 0),
        )
    return metrics


def plot_knockdown_distribution_by_day(table: pd.DataFrame):
    """Plot the Figure S3f all-target knockdown distribution by day."""
    data = table.copy()
    data["knockdown_percent"] = 100.0 * (1.0 - data.fold_change)
    days = ["D0", "D1", "D2", "D3"]
    rng = np.random.default_rng(13)
    fig, ax = plt.subplots(figsize=(5.8, 3.3))
    vectors = []
    for index, day in enumerate(days):
        subset = data[data.day == day].dropna(subset=["knockdown_percent"])
        vectors.append(subset.knockdown_percent.to_numpy())
        jitter = rng.uniform(-0.22, 0.22, len(subset))
        significant = subset.q_value <= 0.05
        ax.scatter(
            index + jitter[~significant.to_numpy()],
            subset.loc[~significant, "knockdown_percent"],
            s=7,
            color="0.65",
            alpha=0.5,
            linewidths=0,
            rasterized=True,
        )
        ax.scatter(
            index + jitter[significant.to_numpy()],
            subset.loc[significant, "knockdown_percent"],
            s=8,
            color="0.15",
            alpha=0.7,
            linewidths=0,
            rasterized=True,
        )
    ax.boxplot(
        vectors,
        positions=np.arange(len(days)),
        widths=0.42,
        showfliers=False,
        patch_artist=False,
        medianprops={"linewidth": 1.0, "color": "black"},
        boxprops={"linewidth": 0.8},
        whiskerprops={"linewidth": 0.8},
        capprops={"linewidth": 0.8},
    )
    ax.axhline(0, color="0.55", linestyle="--", linewidth=0.8)
    ax.set_xticks(np.arange(len(days)), days)
    ax.set_ylabel("Target knockdown efficiency (%)")
    source = data[["target", "day", "fold_change", "knockdown_percent", "q_value"]].copy()
    return fig, source
