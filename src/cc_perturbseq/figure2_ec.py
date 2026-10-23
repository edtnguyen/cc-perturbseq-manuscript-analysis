from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
import hashlib
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import ttest_ind

from .figure2_common import FIGURE2D_TARGETS, bh_adjust, significance_stars
from .program_io import clean_day
from .single_cell_panels import load_guide_target_map

def _assignment_matrix(adata):
    if "guide_assignment" not in adata.obsm:
        raise KeyError("AnnData lacks obsm['guide_assignment']")
    matrix = adata.obsm["guide_assignment"]
    return matrix.tocsr().astype(bool) if sparse.issparse(matrix) else np.asarray(matrix, dtype=bool)

def _column_carriers(assignment, column: int) -> np.ndarray:
    if sparse.issparse(assignment):
        return np.asarray(assignment[:, column].toarray()).ravel().astype(bool)
    return np.asarray(assignment[:, column], dtype=bool).ravel()

def _any_carriers(assignment, columns: np.ndarray) -> np.ndarray:
    if len(columns) == 0:
        return np.zeros(assignment.shape[0], dtype=bool)
    if sparse.issparse(assignment):
        return np.asarray(assignment[:, columns].sum(axis=1)).ravel() > 0
    return np.asarray(assignment[:, columns], dtype=bool).any(axis=1)

def _d3_mask(adata) -> np.ndarray:
    if "day" not in adata.obs:
        raise KeyError("AnnData lacks obs['day']")
    day = pd.Series(adata.obs["day"], index=adata.obs_names).map(clean_day)
    return day.eq("D3").to_numpy()

def _ec_mask(adata) -> tuple[np.ndarray, str]:
    if "cell_type" in adata.obs:
        labels = adata.obs["cell_type"].astype(str)
        mask = labels.str.contains(r"(?:^|[^A-Za-z])art[-_ ]?EC(?:$|[^A-Za-z])", case=False, regex=True)
        if mask.any():
            return mask.to_numpy(), "obs['cell_type'] contains artEC"
    if "leiden_res_0.5" in adata.obs:
        return adata.obs["leiden_res_0.5"].astype(str).eq("0").to_numpy(), "obs['leiden_res_0.5'] == '0' (historical artEC mapping)"
    raise KeyError("Could not identify the arterial-EC state from cell_type or leiden_res_0.5.")

def _non_targeting_columns(guide_map: pd.DataFrame) -> np.ndarray:
    text = guide_map["guide"].astype(str) + " " + guide_map["target"].fillna("").astype(str)
    non_targeting = text.str.contains(
        r"non[-_ ]?target(?:ing)?|\bntc\b|negative|neg[-_ ]?ctrl",
        case=False, regex=True,
    )
    safe = text.str.contains(r"safe[-_ ]?target|\bsafe\b", case=False, regex=True)
    return guide_map.index[non_targeting & ~safe].to_numpy(int)

def compute_d3_guide_ec_statistics(
    adata,
    guide_metadata_path: str | Path,
    targets: Sequence[str],
    *,
    nt_exclusion_targets: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Calculate one D3 EC percentage per guide and compare targets with NTC guides.

    Target-guide percentages include all D3 carriers of that guide. For each NTC
    guide, D3 cells carrying any guide targeting ``nt_exclusion_targets`` are
    excluded before calculating the guide-specific EC percentage.
    """
    targets = [str(target).upper() for target in targets]
    exclusion_targets = [str(target).upper() for target in nt_exclusion_targets]
    assignment = _assignment_matrix(adata)
    guide_map = load_guide_target_map(adata, guide_metadata_path).copy()
    guide_map["target_upper"] = guide_map["target"].fillna("").astype(str).str.upper()
    if assignment.shape[1] != len(guide_map):
        raise ValueError("Guide-assignment columns do not match adata.uns['guide_names'].")

    d3 = _d3_mask(adata)
    ec, ec_rule = _ec_mask(adata)
    nt_columns = _non_targeting_columns(guide_map)
    if len(nt_columns) == 0:
        raise ValueError("No non-targeting guide columns were identified.")

    exclusion_columns = guide_map.index[guide_map["target_upper"].isin(exclusion_targets)].to_numpy(int)
    has_excluded_target = _any_carriers(assignment, exclusion_columns)
    rows: list[dict[str, Any]] = []

    for column in nt_columns:
        carriers = _column_carriers(assignment, int(column)) & d3 & ~has_excluded_target
        n = int(carriers.sum())
        n_ec = int((carriers & ec).sum())
        rows.append({
            "guide_id": str(guide_map.loc[column, "guide"]),
            "target": "Non-targeting",
            "guide_type": "non-targeting",
            "n_D3_cells": n,
            "n_D3_EC_cells": n_ec,
            "n_D3_nonEC_cells": n - n_ec,
            "pct_EC_D3": 100.0 * n_ec / n if n else np.nan,
            "exclusion_rule": "exclude NTC-carrier cells with any guide targeting " + ", ".join(exclusion_targets),
        })

    for target in targets:
        columns = guide_map.index[guide_map["target_upper"].eq(target)].to_numpy(int)
        if len(columns) == 0:
            raise ValueError(f"No guide columns were mapped to {target}.")
        for column in columns:
            carriers = _column_carriers(assignment, int(column)) & d3
            n = int(carriers.sum())
            n_ec = int((carriers & ec).sum())
            rows.append({
                "guide_id": str(guide_map.loc[column, "guide"]),
                "target": target,
                "guide_type": "targeting",
                "n_D3_cells": n,
                "n_D3_EC_cells": n_ec,
                "n_D3_nonEC_cells": n - n_ec,
                "pct_EC_D3": 100.0 * n_ec / n if n else np.nan,
                "exclusion_rule": "none",
            })

    guide_table = pd.DataFrame(rows)
    control_values = guide_table.loc[
        guide_table["guide_type"].eq("non-targeting"), "pct_EC_D3"
    ].dropna().to_numpy(float)
    if len(control_values) < 2:
        raise ValueError("Fewer than two evaluable non-targeting guides remain.")

    stats_rows: list[dict[str, Any]] = []
    for target in targets:
        target_values = guide_table.loc[
            guide_table["target"].eq(target), "pct_EC_D3"
        ].dropna().to_numpy(float)
        if len(target_values) < 2:
            raise ValueError(f"Fewer than two evaluable guides for {target}.")
        result = ttest_ind(target_values, control_values, equal_var=False, nan_policy="raise")
        mean_target = float(target_values.mean())
        mean_control = float(control_values.mean())
        stats_rows.append({
            "gene": target,
            "n_targeting_guides": len(target_values),
            "n_non_targeting_guides": len(control_values),
            "mean_target_EC_pct": mean_target,
            "median_target_EC_pct": float(np.median(target_values)),
            "sd_target_EC_pct": float(np.std(target_values, ddof=1)),
            "mean_non_targeting_EC_pct": mean_control,
            "median_non_targeting_EC_pct": float(np.median(control_values)),
            "sd_non_targeting_EC_pct": float(np.std(control_values, ddof=1)),
            "difference_target_minus_control_pp": mean_target - mean_control,
            "effect_log2_target_over_control": float(np.log2(mean_target / mean_control)),
            "welch_t_statistic": float(result.statistic),
            "welch_df": float(result.df),
            "welch_pvalue_two_sided": float(result.pvalue),
        })
    stats = pd.DataFrame(stats_rows)
    stats["FDR_BH"] = bh_adjust(stats["welch_pvalue_two_sided"].to_numpy(float))
    stats["significance"] = stats["FDR_BH"].map(significance_stars)
    stats["BH_number_tests"] = len(stats)
    metadata = {
        "ec_definition": ec_rule,
        "n_D3_cells": int(d3.sum()),
        "n_non_targeting_guides_designed": int(len(nt_columns)),
        "n_non_targeting_guides_evaluable": int(len(control_values)),
        "targets": targets,
        "nt_exclusion_targets": exclusion_targets,
    }
    return guide_table, stats, metadata

def plot_figure2d(guide_table: pd.DataFrame, stats: pd.DataFrame):
    order = ["Non-targeting", *FIGURE2D_TARGETS]
    colors = {"Non-targeting": "white", **{target: "#059C48" for target in FIGURE2D_TARGETS}}
    datasets = [guide_table.loc[guide_table["target"].eq(target), "pct_EC_D3"].dropna().to_numpy(float) for target in order]
    stats_index = stats.set_index("gene")
    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    positions = np.arange(len(order))
    box = ax.boxplot(
        datasets, positions=positions, widths=0.48, whis=1.5, patch_artist=True,
        showfliers=False,
        boxprops={"edgecolor": "black", "linewidth": 1.2},
        whiskerprops={"color": "black", "linewidth": 1.1},
        capprops={"color": "black", "linewidth": 1.1},
        medianprops={"color": "black", "linewidth": 1.5},
    )
    for patch, target in zip(box["boxes"], order):
        patch.set_facecolor(colors[target])
        patch.set_alpha(0.9)

    jitter_rows: list[dict[str, Any]] = []
    for x, target, values in zip(positions, order, datasets):
        seed = int.from_bytes(hashlib.sha256(target.encode()).digest()[:4], "little")
        rng = np.random.default_rng(seed)
        jitter = np.clip(rng.normal(0, 0.06 if target != "Non-targeting" else 0.10, len(values)), -0.18, 0.18)
        point_color = "#111111" if target == "Non-targeting" else "#059C48"
        ax.scatter(np.full(len(values), x) + jitter, values, s=16, c=point_color, alpha=0.55, linewidths=0, zorder=4)
        guide_rows = guide_table.loc[guide_table["target"].eq(target)].dropna(subset=["pct_EC_D3"])
        for (_, row), offset in zip(guide_rows.iterrows(), jitter):
            jitter_rows.append({"target": target, "guide_id": row["guide_id"], "pct_EC_D3": row["pct_EC_D3"], "x": x + float(offset)})
        if target != "Non-targeting":
            star = str(stats_index.loc[target, "significance"])
            if star:
                ax.text(x, 96, star, ha="center", va="bottom", fontsize=11, fontweight="bold")

    control_mean = float(np.mean(datasets[0]))
    ax.axhline(control_mean, color="#999999", linestyle="--", linewidth=1.2)
    ax.set_xticks(positions, [f"{target}\nn={len(values)}" for target, values in zip(order, datasets)], rotation=55, ha="right", fontsize=8)
    ax.set_ylabel("EC differentiation efficiency (%)", fontsize=10)
    ax.set_xlabel("gRNAs", fontsize=9)
    ax.set_ylim(0, 105)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(0.01, 0.98, "d", fontsize=15, fontweight="bold", ha="left", va="top")
    fig.tight_layout(rect=[0.06, 0.02, 1, 0.96])
    return fig, pd.DataFrame(jitter_rows)
