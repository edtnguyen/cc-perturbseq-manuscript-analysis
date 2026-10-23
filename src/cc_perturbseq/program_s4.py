from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .program_io import clean_day, infer_column, infer_day_from_path, read_table


def plot_cnmf_diagnostics(metrics: pd.DataFrame, *, selected_k: int = 50):
    k_column = infer_column(metrics, ["k", "n_components", "components"])
    specifications = [
        ("stability", ["stability", "silhouette"]),
        ("reconstruction error", ["reconstruction_error", "reconstruction", "loss"]),
        ("explained variance", ["explained_variance", "variance_explained"]),
        ("GO term enrichment", ["go_enrichment", "enriched_gene_sets", "gene_set_enrichment", "unique_go"]),
    ]
    figure, axes = plt.subplots(1, 4, figsize=(10.5, 2.6), constrained_layout=True)
    for axis, (title, candidates) in zip(axes, specifications):
        column = infer_column(metrics, candidates)
        x = pd.to_numeric(metrics[k_column], errors="coerce")
        y = pd.to_numeric(metrics[column], errors="coerce")
        order = np.argsort(x)
        axis.plot(x.iloc[order], y.iloc[order], marker="o", markersize=3, linewidth=1)
        axis.axvline(selected_k, color="black", linestyle="--", linewidth=0.8)
        axis.set_title(title, fontsize=8)
        axis.set_xlabel("Number of components (K)", fontsize=7)
        axis.tick_params(labelsize=6)
    return figure


def _pvalues(path: str | Path, kind: str) -> pd.DataFrame:
    path = Path(path)
    files = (
        sorted(path.rglob("*.csv")) + sorted(path.rglob("*.tsv"))
        if path.is_dir()
        else [path]
    )
    pieces: list[pd.DataFrame] = []
    for file in files:
        table = read_table(file)
        p_column = infer_column(table, ["p_value", "pvalue", "pval", "p"], required=False)
        if p_column is None:
            numeric = [column for column in table.columns if pd.api.types.is_numeric_dtype(table[column])]
            if len(numeric) != 1:
                continue
            p_column = numeric[0]
        day_column = infer_column(table, ["day", "timepoint"], required=False)
        frame = pd.DataFrame({"p_value": pd.to_numeric(table[p_column], errors="coerce")})
        frame["day"] = table[day_column].map(clean_day) if day_column else infer_day_from_path(file)
        frame["kind"] = kind
        pieces.append(frame)
    if not pieces:
        raise ValueError(f"No p-value column found in {path}")
    return pd.concat(pieces, ignore_index=True).dropna(subset=["p_value", "day"])


def plot_crt_qq(real_path: str | Path, ntc_path: str | Path):
    data = pd.concat([_pvalues(real_path, "Real target"), _pvalues(ntc_path, "Matched NTC")], ignore_index=True)
    figure, axes = plt.subplots(1, 4, figsize=(10, 2.6), constrained_layout=True)
    rows: list[pd.DataFrame] = []
    for axis, day in zip(axes, ["D0", "D1", "D2", "D3"]):
        maxima = [1.0]
        for kind in ["Real target", "Matched NTC"]:
            values = data.loc[(data.day == day) & (data.kind == kind), "p_value"].clip(lower=np.nextafter(0, 1), upper=1).sort_values().to_numpy()
            if not len(values):
                continue
            expected = (np.arange(1, len(values) + 1) - 0.5) / len(values)
            rows.append(pd.DataFrame({"day": day, "kind": kind, "expected": expected, "observed": values}))
            x = -np.log10(expected)
            y = -np.log10(values)
            maxima.extend([x.max(), y.max()])
            axis.plot(x, y, linewidth=1, label=kind)
        limit = max(maxima)
        axis.plot([0, limit], [0, limit], color="0.5", linestyle="--", linewidth=0.8)
        axis.set_title(day, fontsize=8)
        axis.set_xlabel("Expected −log10(p)", fontsize=7)
        axis.tick_params(labelsize=6)
        if day == "D0":
            axis.set_ylabel("Observed −log10(p)", fontsize=7)
            axis.legend(frameon=False, fontsize=6)
    return figure, pd.concat(rows, ignore_index=True)


def plot_replicate_effects(path: str | Path):
    table = read_table(path)
    day_column = infer_column(table, ["day", "timepoint"], required=False)
    x_column = infer_column(table, ["replicate_1", "rep1", "beta_1", "effect_1", "x"], required=False)
    y_column = infer_column(table, ["replicate_2", "rep2", "beta_2", "effect_2", "y"], required=False)
    if not x_column or not y_column:
        replicate_column = infer_column(table, ["replicate", "biological_sample", "sample"])
        effect_column = infer_column(table, ["effect", "beta", "obs_coef", "coefficient"])
        target_column = infer_column(table, ["target", "gene"])
        program_column = infer_column(table, ["program", "program_id"])
        table["day"] = table[day_column].map(clean_day) if day_column else None
        wide = table.pivot_table(index=["day", target_column, program_column], columns=replicate_column, values=effect_column, aggfunc="mean")
        if wide.shape[1] < 2:
            raise ValueError("Fewer than two replicate columns.")
        pairs = wide.iloc[:, :2].dropna().reset_index().rename(columns={wide.columns[0]: "replicate_1", wide.columns[1]: "replicate_2"})
    else:
        pairs = pd.DataFrame({
            "day": table[day_column].map(clean_day) if day_column else None,
            "replicate_1": pd.to_numeric(table[x_column], errors="coerce"),
            "replicate_2": pd.to_numeric(table[y_column], errors="coerce"),
        }).dropna()

    figure, axes = plt.subplots(2, 2, figsize=(5.3, 5), constrained_layout=True)
    summary: list[dict[str, float | int | str]] = []
    for axis, day in zip(axes.ravel(), ["D0", "D1", "D2", "D3"]):
        subset = pairs[pairs.day == day]
        axis.scatter(subset.replicate_1, subset.replicate_2, s=5, alpha=0.45, rasterized=True)
        if len(subset) >= 2:
            correlation = subset[["replicate_1", "replicate_2"]].corr().iloc[0, 1]
            summary.append({"day": day, "pearson_r": correlation, "n": len(subset)})
            axis.text(0.03, 0.95, f"Pearson r={correlation:.2f}\nn={len(subset):,}", transform=axis.transAxes, va="top", fontsize=7)
        limits = (
            np.nanmin(subset[["replicate_1", "replicate_2"]]) if len(subset) else -1,
            np.nanmax(subset[["replicate_1", "replicate_2"]]) if len(subset) else 1,
        )
        axis.plot(limits, limits, color="0.6", linestyle="--", linewidth=0.8)
        axis.set_title(day, fontsize=8)
        axis.set_xlabel("Biological replicate 1", fontsize=7)
        axis.set_ylabel("Biological replicate 2", fontsize=7)
        axis.tick_params(labelsize=6)
    return figure, pd.DataFrame(summary)


def plot_program_gene_overlap(spectra: pd.DataFrame, *, top_n: int = 300):
    programs = spectra.index.astype(int).tolist()
    sets = {program: set(spectra.loc[program].nlargest(top_n).index.astype(str)) for program in programs}
    overlap = pd.DataFrame(index=programs, columns=programs, dtype=float)
    for first in programs:
        for second in programs:
            overlap.loc[first, second] = len(sets[first] & sets[second])
    figure, axis = plt.subplots(figsize=(6.2, 5.6))
    image = axis.imshow(overlap, cmap="Reds", aspect="auto")
    axis.set_title(f"Shared top-{top_n} genes between programs", fontsize=9)
    axis.set_xticks(np.arange(len(programs)), [f"P{x}" for x in programs], rotation=90, fontsize=5)
    axis.set_yticks(np.arange(len(programs)), [f"P{x}" for x in programs], fontsize=5)
    figure.colorbar(image, ax=axis, fraction=0.03, pad=0.02, label="Shared genes")
    return figure, overlap


def plot_usage_correlation(usage: pd.DataFrame):
    correlation = usage.corr(method="pearson")
    figure, axis = plt.subplots(figsize=(6.2, 5.6))
    image = axis.imshow(correlation, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    axis.set_title(f"Program–program usage correlation across all cells (n={len(usage):,})", fontsize=9)
    axis.set_xticks(np.arange(len(correlation)), correlation.columns, rotation=90, fontsize=5)
    axis.set_yticks(np.arange(len(correlation)), correlation.index, fontsize=5)
    figure.colorbar(image, ax=axis, fraction=0.03, pad=0.02, label="Pearson r")
    return figure, correlation
