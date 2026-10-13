from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
import re
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from .program_io import clean_day, infer_column, parse_program, read_table


def load_integrated_anndata(
    adata_path: str | Path,
    usage_path: str | Path | None = None,
):
    import anndata as ad

    adata = ad.read_h5ad(adata_path)
    if usage_path is not None and Path(usage_path).exists():
        usage = pd.read_csv(usage_path, sep=r"\s+|\t|,", engine="python", index_col=0)
        usage = usage.reindex(adata.obs_names)
        usage.columns = [
            f"P{parse_program(column)}" if parse_program(column) is not None else str(column)
            for column in usage.columns
        ]
        adata.obsm["cnmf_usage"] = usage
    elif "cnmf_usage" in adata.obsm:
        matrix = adata.obsm["cnmf_usage"]
        if isinstance(matrix, pd.DataFrame):
            usage = matrix.copy()
        else:
            names = adata.uns.get("cnmf_program_names")
            if names is None or len(names) != matrix.shape[1]:
                names = [f"P{i + 1}" for i in range(matrix.shape[1])]
            usage = pd.DataFrame(matrix, index=adata.obs_names, columns=names)
        usage.columns = [
            f"P{parse_program(column)}" if parse_program(column) is not None else str(column)
            for column in usage.columns
        ]
        adata.obsm["cnmf_usage"] = usage
    else:
        usage = None
    return adata, usage


def _umap(adata) -> np.ndarray:
    if "X_umap" not in adata.obsm:
        raise KeyError("AnnData lacks obsm['X_umap']")
    return np.asarray(adata.obsm["X_umap"])


def _gene_symbols(adata) -> pd.Index:
    for column in ("gene_symbol", "gene_name", "symbol", "feature_name"):
        if column in adata.var.columns:
            return pd.Index(adata.var[column].astype(str))
    return pd.Index(adata.var_names.astype(str))


def expression_vector(adata, gene: str) -> np.ndarray:
    symbols = _gene_symbols(adata)
    positions = np.flatnonzero(symbols == gene)
    if not len(positions):
        raise KeyError(f"Gene {gene!r} is absent from the AnnData feature metadata")
    matrix = adata[:, positions].X
    if sparse.issparse(matrix):
        values = np.asarray(matrix.sum(axis=1)).ravel()
    else:
        values = np.asarray(matrix).sum(axis=1)
    return values.astype(float)


def program_vector(usage: pd.DataFrame, program: int | str) -> np.ndarray:
    number = parse_program(program)
    candidates = [str(program)]
    if number is not None:
        candidates.extend([f"P{number}", f"Program {number}", f"program_{number}"])
    for candidate in candidates:
        if candidate in usage.columns:
            return pd.to_numeric(usage[candidate], errors="coerce").fillna(0).to_numpy(float)
    for column in usage.columns:
        if parse_program(column) == number:
            return pd.to_numeric(usage[column], errors="coerce").fillna(0).to_numpy(float)
    raise KeyError(f"Program {program!r} is absent from usage columns")


def plot_umap_feature_grid(
    adata,
    usage: pd.DataFrame,
    features: Sequence[tuple[str, str | int]],
    *,
    point_size: float = 1.2,
    background_size: float = 0.08,
    vmax_quantile: float = 0.99,
    max_background_cells: int = 70000,
    seed: int = 17,
):
    """Plot a grid of program-usage and expression layers on shared UMAP limits.

    Each feature tuple is (kind, name), where kind is ``program`` or ``gene``.
    """
    coords = _umap(adata)
    rng = np.random.default_rng(seed)
    background_index = np.arange(len(coords))
    if len(background_index) > max_background_cells:
        background_index = rng.choice(background_index, max_background_cells, replace=False)

    ncols = min(4, len(features))
    nrows = int(np.ceil(len(features) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(2.5 * ncols, 2.25 * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    source: dict[str, pd.DataFrame] = {}
    for axis, (kind, name) in zip(axes.ravel(), features):
        if kind.lower() == "program":
            values = program_vector(usage, name)
            title = f"P{parse_program(name)}"
        elif kind.lower() == "gene":
            values = expression_vector(adata, str(name))
            title = str(name)
        else:
            raise ValueError(f"Unknown feature kind: {kind}")
        finite = values[np.isfinite(values)]
        vmax = float(np.quantile(finite, vmax_quantile)) if finite.size else 1.0
        vmax = max(vmax, np.nextafter(0, 1))
        order = np.argsort(values)
        axis.scatter(
            coords[background_index, 0],
            coords[background_index, 1],
            s=background_size,
            color="0.88",
            linewidths=0,
            rasterized=True,
        )
        scatter = axis.scatter(
            coords[order, 0],
            coords[order, 1],
            c=np.clip(values[order], 0, vmax),
            cmap="YlGnBu",
            norm=Normalize(0, vmax),
            s=point_size,
            linewidths=0,
            rasterized=True,
        )
        axis.set_title(title, fontsize=8, fontstyle="italic" if kind.lower() == "gene" else "normal")
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        fig.colorbar(scatter, ax=axis, fraction=0.045, pad=0.01)
        source[f"{kind}_{name}"] = pd.DataFrame(
            {"umap_1": coords[:, 0], "umap_2": coords[:, 1], "value": values},
            index=adata.obs_names,
        )
    for axis in axes.ravel()[len(features) :]:
        axis.set_axis_off()
    return fig, source


def plot_program_atlas(
    adata,
    usage: pd.DataFrame,
    *,
    programs_by_day: Mapping[str, Sequence[int]],
    program_labels: Mapping[int, str] | None = None,
    max_cells: int = 100000,
    seed: int = 7,
):
    programs = [int(program) for values in programs_by_day.values() for program in values]
    matrix = np.column_stack([program_vector(usage, program) for program in programs])
    winner_index = matrix.argmax(axis=1)
    winner = np.asarray(programs)[winner_index]
    maximum = matrix.max(axis=1)
    winner = winner.astype(object)
    winner[maximum <= 0] = "Other"
    coords = _umap(adata)
    rng = np.random.default_rng(seed)
    draw = np.arange(len(coords))
    if len(draw) > max_cells:
        draw = rng.choice(draw, max_cells, replace=False)
    draw = draw[np.argsort(maximum[draw])]

    palette = plt.get_cmap("tab20", len(programs))
    colors = {program: palette(index) for index, program in enumerate(programs)}
    colors["Other"] = (0.82, 0.82, 0.82, 1)
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    for category in ["Other", *programs]:
        mask = winner[draw] == category
        if not mask.any():
            continue
        label = "Other" if category == "Other" else (program_labels or {}).get(category, f"P{category}")
        ax.scatter(
            coords[draw[mask], 0],
            coords[draw[mask], 1],
            s=1.0 if category != "Other" else 0.35,
            color=colors[category],
            linewidths=0,
            rasterized=True,
            label=label,
        )
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(
        frameon=False,
        fontsize=5.5,
        markerscale=4,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        ncol=1,
    )
    source = pd.DataFrame(
        {
            "umap_1": coords[:, 0],
            "umap_2": coords[:, 1],
            "program": winner,
            "max_selected_usage": maximum,
        },
        index=adata.obs_names,
    )
    return fig, source


def load_guide_target_map(
    adata,
    metadata_path: str | Path | None = None,
) -> pd.DataFrame:
    names = pd.Index([str(value) for value in adata.uns["guide_names"]], name="guide")
    output = pd.DataFrame({"guide": names})
    if metadata_path is not None and Path(metadata_path).exists():
        metadata = read_table(metadata_path)
        guide_col = infer_column(metadata, ["guide_id", "guide", "sgRNA", "grna"])
        target_col = infer_column(
            metadata,
            ["putative_target_genes", "intended_target_name", "target", "target_name", "gene"],
        )
        mapping = metadata[[guide_col, target_col]].rename(
            columns={guide_col: "guide", target_col: "target"}
        )
        mapping["guide"] = mapping.guide.astype(str)
        mapping["target"] = mapping.target.astype(str).str.split(r"[;,|/]").str[0].str.strip()
        output = output.merge(mapping.drop_duplicates("guide"), on="guide", how="left")
    else:
        output["target"] = output.guide.str.replace(r"(?:_guide|_g|[-_]sg)\d+.*$", "", regex=True)
        output["target"] = output.target.str.replace(r"[-_]\d+$", "", regex=True)
    control_pattern = r"non[-_ ]?target|safe[-_ ]?target|\bntc\b|negative|scram|aavs1|control"
    output["is_control"] = output.guide.str.contains(control_pattern, case=False, regex=True)
    output["is_control"] |= output.target.str.contains(control_pattern, case=False, regex=True, na=False)
    return output


def _assignment(adata):
    if "guide_assignment" not in adata.obsm:
        raise KeyError("AnnData lacks obsm['guide_assignment']")
    assignment = adata.obsm["guide_assignment"]
    return assignment.tocsr().astype(bool) if sparse.issparse(assignment) else np.asarray(assignment, dtype=bool)


def _any_assignment(assignment, columns: np.ndarray) -> np.ndarray:
    if not len(columns):
        return np.zeros(assignment.shape[0], dtype=bool)
    if sparse.issparse(assignment):
        return np.asarray(assignment[:, columns].sum(axis=1)).ravel() > 0
    return assignment[:, columns].any(axis=1)


def _zscore_usage(usage: pd.DataFrame) -> np.ndarray:
    values = usage.to_numpy(float)
    mean = np.nanmean(values, axis=0)
    std = np.nanstd(values, axis=0)
    std[~np.isfinite(std) | (std == 0)] = 1.0
    return np.nan_to_num((values - mean) / std)


def carrier_seeded_neighborhood(
    adata,
    usage: pd.DataFrame,
    guide_map: pd.DataFrame,
    targets: Sequence[str],
    *,
    k: int = 300,
    pseudocount: float = 0.5,
) -> pd.DataFrame:
    """Compute the manuscript carrier-seeded same-day neighborhood statistic."""
    from sklearn.neighbors import NearestNeighbors

    assignment = _assignment(adata)
    days = pd.Series(adata.obs["day"], index=adata.obs_names).map(clean_day)
    embedding = _zscore_usage(usage)
    control_columns = guide_map.index[guide_map.is_control].to_numpy(int)
    has_control = _any_assignment(assignment, control_columns)
    results: list[pd.DataFrame] = []

    for target in targets:
        target_columns = guide_map.index[guide_map.target.eq(target)].to_numpy(int)
        if not len(target_columns):
            continue
        has_target = _any_assignment(assignment, target_columns)
        target_seed = has_target & ~has_control
        target_control = has_control & ~has_target
        for day in sorted(days.dropna().unique()):
            day_mask = days.eq(day).to_numpy()
            day_indices = np.flatnonzero(day_mask)
            seeds = np.flatnonzero(target_seed & day_mask)
            if not len(seeds):
                continue
            T = int((target_seed & day_mask).sum())
            C = int((target_control & day_mask).sum())
            if T == 0 or C == 0:
                continue
            n_neighbors = min(k + 1, len(day_indices))
            model = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
            model.fit(embedding[day_indices])
            _, local_indices = model.kneighbors(embedding[seeds])
            neighbors = day_indices[local_indices]
            rows: list[dict[str, Any]] = []
            for seed, neighbor_row in zip(seeds, neighbors):
                neighbor_row = neighbor_row[neighbor_row != seed][:k]
                t = int(target_seed[neighbor_row].sum())
                c = int(target_control[neighbor_row].sum())
                enrichment = np.log2(
                    ((t + pseudocount) / (c + pseudocount))
                    / ((T + pseudocount) / (C + pseudocount))
                )
                odds_table = [[t, c], [max(T - t, 0), max(C - c, 0)]]
                p_value = fisher_exact(odds_table, alternative="two-sided").pvalue
                rows.append(
                    {
                        "cell": adata.obs_names[seed],
                        "cell_index": seed,
                        "target": target,
                        "day": day,
                        "target_neighbors": t,
                        "control_neighbors": c,
                        "neighbor_target_fraction": t / max(len(neighbor_row), 1),
                        "log2_enrichment": enrichment,
                        "p_value": p_value,
                    }
                )
            frame = pd.DataFrame(rows)
            if not frame.empty:
                frame["fdr"] = multipletests(frame.p_value, method="fdr_bh")[1]
                results.append(frame)
    if not results:
        return pd.DataFrame(
            columns=[
                "cell", "cell_index", "target", "day", "target_neighbors",
                "control_neighbors", "neighbor_target_fraction", "log2_enrichment",
                "p_value", "fdr",
            ]
        )
    return pd.concat(results, ignore_index=True)


def plot_carrier_neighborhood_grid(
    adata,
    neighborhoods: pd.DataFrame,
    targets: Sequence[str],
    *,
    clip: float = 2.0,
    max_background_cells: int = 70000,
    seed: int = 19,
):
    coords = _umap(adata)
    rng = np.random.default_rng(seed)
    background = np.arange(len(coords))
    if len(background) > max_background_cells:
        background = rng.choice(background, max_background_cells, replace=False)
    ncols = min(4, len(targets))
    nrows = int(np.ceil(len(targets) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(2.5 * ncols, 2.3 * nrows),
        squeeze=False,
        constrained_layout=True,
    )
    for axis, target in zip(axes.ravel(), targets):
        axis.scatter(
            coords[background, 0],
            coords[background, 1],
            s=0.08,
            color="0.86",
            linewidths=0,
            rasterized=True,
        )
        data = neighborhoods[neighborhoods.target == target].copy()
        if not data.empty:
            values = data.log2_enrichment.clip(-clip, clip).to_numpy()
            indices = data.cell_index.to_numpy(int)
            sizes = 4 + 35 * np.sqrt(data.neighbor_target_fraction.clip(lower=0).to_numpy())
            significant = data.fdr.to_numpy() < 0.05
            scatter = axis.scatter(
                coords[indices, 0],
                coords[indices, 1],
                c=values,
                cmap="RdBu_r",
                norm=TwoSlopeNorm(vmin=-clip, vcenter=0, vmax=clip),
                s=sizes,
                linewidths=np.where(significant, 0.45, 0),
                edgecolors=np.where(significant, "black", "none"),
                rasterized=True,
            )
            fig.colorbar(scatter, ax=axis, fraction=0.045, pad=0.01)
        axis.set_title(f"{target} CRISPRi", fontsize=8, fontstyle="italic")
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
    for axis in axes.ravel()[len(targets) :]:
        axis.set_axis_off()
    return fig
