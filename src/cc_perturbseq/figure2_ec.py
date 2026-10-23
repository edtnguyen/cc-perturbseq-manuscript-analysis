from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
import hashlib
import re
from typing import Any

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    annotate,
    element_blank,
    element_rect,
    element_text,
    geom_boxplot,
    geom_hline,
    geom_point,
    geom_text,
    ggplot,
    labs,
    scale_fill_manual,
    scale_x_continuous,
    scale_y_continuous,
    theme,
    theme_classic,
)
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


def _column_counts(
    matrix,
    eligible: np.ndarray,
    ec_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if sparse.issparse(matrix):
        n_total = np.asarray(matrix[eligible].sum(axis=0)).ravel().astype(int)
        n_ec = np.asarray(matrix[eligible & ec_mask].sum(axis=0)).ravel().astype(int)
    else:
        n_total = np.asarray(matrix[eligible].sum(axis=0)).ravel().astype(int)
        n_ec = np.asarray(matrix[eligible & ec_mask].sum(axis=0)).ravel().astype(int)
    pct = np.divide(
        100.0 * n_ec,
        n_total,
        out=np.full(len(n_total), np.nan, dtype=float),
        where=n_total > 0,
    )
    return n_total, n_ec, pct


def _d3_mask(adata) -> np.ndarray:
    if "day" not in adata.obs:
        raise KeyError("AnnData lacks obs['day']")
    day = pd.Series(adata.obs["day"], index=adata.obs_names).map(clean_day)
    return day.eq("D3").to_numpy()


def _ec_mask(adata) -> tuple[np.ndarray, str]:
    if "cell_type" in adata.obs:
        labels = adata.obs["cell_type"].astype(str)
        mask = labels.str.contains(
            r"(?:^|[^A-Za-z])art[-_ ]?EC(?:$|[^A-Za-z])",
            case=False,
            regex=True,
        )
        if mask.any():
            return mask.to_numpy(), "obs['cell_type'] contains artEC"
    if "leiden_res_0.5" in adata.obs:
        return (
            adata.obs["leiden_res_0.5"].astype(str).eq("0").to_numpy(),
            "obs['leiden_res_0.5'] == '0' (historical artEC mapping)",
        )
    raise KeyError("Could not identify the arterial-EC state from cell_type or leiden_res_0.5.")


def _canonical_target(value: Any) -> str:
    target = str(value).strip().upper()
    return "NEDD8" if target == "NEDD8-MDP1" else target


def _typed_guide_map(adata, guide_metadata_path: str | Path) -> pd.DataFrame:
    guide_map = load_guide_target_map(adata, guide_metadata_path).copy()
    guide_map["guide"] = guide_map["guide"].astype(str)
    guide_map["target"] = guide_map["target"].fillna("").map(_canonical_target)
    combined = (
        guide_map["guide"] + " " + guide_map["target"]
    ).str.upper().str.replace("_", "-", regex=False)
    guide_map["is_non_targeting"] = combined.str.contains(
        r"NON[- ]?TARGET(?:ING)?|\bNTC\b|NEG[- ]?CTRL|NEGATIVE",
        case=False,
        regex=True,
    )
    guide_map["is_safe_targeting"] = combined.str.contains(
        r"SAFE[- ]?TARGET(?:ING)?|\bSAFE\b",
        case=False,
        regex=True,
    )
    guide_map["is_control"] = guide_map["is_non_targeting"] | guide_map["is_safe_targeting"]
    guide_map["is_targeting"] = ~guide_map["is_control"] & guide_map["target"].ne("")
    return guide_map


def _d3_inputs(adata, guide_metadata_path: str | Path):
    assignment = _assignment_matrix(adata)
    guide_map = _typed_guide_map(adata, guide_metadata_path)
    if assignment.shape[1] != len(guide_map):
        raise ValueError("Guide-assignment columns do not match adata.uns['guide_names'].")
    d3 = _d3_mask(adata)
    ec_full, ec_rule = _ec_mask(adata)
    assignment_d3 = assignment[d3]
    if sparse.issparse(assignment_d3):
        assignment_d3 = assignment_d3.tocsr().astype(bool)
    else:
        assignment_d3 = np.asarray(assignment_d3, dtype=bool)
    return assignment_d3, ec_full[d3], guide_map, int(d3.sum()), ec_rule


def compute_figure2b_ec_statistics(
    adata,
    guide_metadata_path: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Calculate Figure 2b D3 EC effects for all 297 targets.

    The control distribution contains all 120 controls: 60 non-targeting and
    60 safe-targeting guides. For each target gene, targeting-guide percentages
    exclude cells carrying any control guide, while each control-guide
    percentage excludes cells carrying a guide to that target. Welch tests are
    performed on guide-level percentages and BH-adjusted across 297 targets.
    """
    assignment, ec, guide_map, n_d3, ec_rule = _d3_inputs(adata, guide_metadata_path)
    control_columns = guide_map.index[guide_map["is_control"]].to_numpy(int)
    nt_columns = guide_map.index[guide_map["is_non_targeting"]].to_numpy(int)
    safe_columns = guide_map.index[guide_map["is_safe_targeting"]].to_numpy(int)
    if len(control_columns) != 120 or len(nt_columns) != 60 or len(safe_columns) != 60:
        raise ValueError(
            "Expected 120 controls (60 non-targeting and 60 safe-targeting); "
            f"found total={len(control_columns)}, non-targeting={len(nt_columns)}, safe={len(safe_columns)}."
        )

    targeting_map = guide_map.loc[guide_map["is_targeting"]].copy()
    genes = sorted(targeting_map["target"].unique())
    if len(genes) != 297:
        raise ValueError(f"Expected 297 canonical target genes; found {len(genes)}.")

    control_matrix = assignment[:, control_columns]
    has_any_control = _any_carriers(assignment, control_columns)
    target_eligible = ~has_any_control

    stats_rows: list[dict[str, Any]] = []
    guide_rows: list[dict[str, Any]] = []
    for gene in genes:
        target_columns = targeting_map.index[targeting_map["target"].eq(gene)].to_numpy(int)
        target_matrix = assignment[:, target_columns]
        has_target = _any_carriers(assignment, target_columns)

        target_total, target_ec, target_pct = _column_counts(target_matrix, target_eligible, ec)
        control_total, control_ec, control_pct = _column_counts(control_matrix, ~has_target, ec)
        target_values = target_pct[np.isfinite(target_pct)]
        control_values = control_pct[np.isfinite(control_pct)]
        if len(target_values) < 2 or len(control_values) < 2:
            raise ValueError(
                f"Insufficient evaluable guides for {gene}: target={len(target_values)}, controls={len(control_values)}."
            )

        result = ttest_ind(target_values, control_values, equal_var=False, nan_policy="raise")
        mean_target = float(target_values.mean())
        mean_control = float(control_values.mean())
        stats_rows.append(
            {
                "gene": gene,
                "n_targeting_guides": len(target_values),
                "n_control_guides": len(control_values),
                "mean_target_EC_pct": mean_target,
                "median_target_EC_pct": float(np.median(target_values)),
                "sd_target_EC_pct": float(np.std(target_values, ddof=1)),
                "mean_control_EC_pct": mean_control,
                "median_control_EC_pct": float(np.median(control_values)),
                "sd_control_EC_pct": float(np.std(control_values, ddof=1)),
                "difference_target_minus_control_pp": mean_target - mean_control,
                "effect_log2_target_over_control": float(np.log2(mean_target / mean_control)),
                "welch_t_statistic": float(result.statistic),
                "welch_df": float(result.df),
                "welch_pvalue_two_sided": float(result.pvalue),
            }
        )

        for column, n, n_ec, pct in zip(target_columns, target_total, target_ec, target_pct):
            guide_rows.append(
                {
                    "comparison_gene": gene,
                    "guide_id": str(guide_map.loc[column, "guide"]),
                    "guide_target": gene,
                    "guide_type": "targeting",
                    "n_D3_cells": int(n),
                    "n_D3_EC_cells": int(n_ec),
                    "n_D3_nonEC_cells": int(n - n_ec),
                    "pct_EC_D3": float(pct) if np.isfinite(pct) else np.nan,
                    "exclusion_rule": "exclude target-guide carriers with any of 120 control guides",
                }
            )
        for column, n, n_ec, pct in zip(control_columns, control_total, control_ec, control_pct):
            guide_rows.append(
                {
                    "comparison_gene": gene,
                    "guide_id": str(guide_map.loc[column, "guide"]),
                    "guide_target": "NON-TARGETING" if guide_map.loc[column, "is_non_targeting"] else "SAFE",
                    "guide_type": "control",
                    "n_D3_cells": int(n),
                    "n_D3_EC_cells": int(n_ec),
                    "n_D3_nonEC_cells": int(n - n_ec),
                    "pct_EC_D3": float(pct) if np.isfinite(pct) else np.nan,
                    "exclusion_rule": f"exclude control-guide carriers with any guide targeting {gene}",
                }
            )

    stats = pd.DataFrame(stats_rows).sort_values("gene").reset_index(drop=True)
    stats["FDR_BH_297_genes"] = bh_adjust(stats["welch_pvalue_two_sided"].to_numpy(float))
    stats["significance"] = stats["FDR_BH_297_genes"].map(significance_stars)
    stats["BH_number_tests"] = len(stats)
    metadata = {
        "ec_definition": ec_rule,
        "n_D3_cells": n_d3,
        "n_non_targeting_guides_designed": len(nt_columns),
        "n_safe_targeting_guides_designed": len(safe_columns),
        "n_control_guides_designed": len(control_columns),
        "n_target_genes_tested": len(stats),
        "target_exclusion": "exclude cells carrying any of 120 control guides",
        "control_exclusion": "target-specific exclusion of cells carrying the target under evaluation",
    }
    return pd.DataFrame(guide_rows), stats, metadata


def compute_d3_guide_ec_statistics(
    adata,
    guide_metadata_path: str | Path,
    targets: Sequence[str],
    *,
    nt_exclusion_targets: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Calculate panel-d guide percentages and Welch tests against 60 NTC guides.

    Targeting-guide percentages include all D3 carriers of the individual guide.
    For every non-targeting guide, D3 cells carrying any guide targeting
    ``nt_exclusion_targets`` are excluded before computing the guide percentage.
    """
    targets = [_canonical_target(target) for target in targets]
    exclusion_targets = [_canonical_target(target) for target in nt_exclusion_targets]
    assignment, ec, guide_map, n_d3, ec_rule = _d3_inputs(adata, guide_metadata_path)
    nt_columns = guide_map.index[guide_map["is_non_targeting"]].to_numpy(int)
    if len(nt_columns) != 60:
        raise ValueError(f"Expected 60 non-targeting guides; found {len(nt_columns)}.")

    exclusion_columns = guide_map.index[
        guide_map["is_targeting"] & guide_map["target"].isin(exclusion_targets)
    ].to_numpy(int)
    has_excluded_target = _any_carriers(assignment, exclusion_columns)
    rows: list[dict[str, Any]] = []

    for column in nt_columns:
        carriers = _column_carriers(assignment, int(column)) & ~has_excluded_target
        n = int(carriers.sum())
        n_ec = int((carriers & ec).sum())
        rows.append(
            {
                "guide_id": str(guide_map.loc[column, "guide"]),
                "target": "Non-targeting",
                "guide_type": "non-targeting",
                "n_D3_cells": n,
                "n_D3_EC_cells": n_ec,
                "n_D3_nonEC_cells": n - n_ec,
                "pct_EC_D3": 100.0 * n_ec / n if n else np.nan,
                "exclusion_rule": "exclude NTC-carrier cells with any guide targeting " + ", ".join(exclusion_targets),
            }
        )

    for target in targets:
        columns = guide_map.index[
            guide_map["is_targeting"] & guide_map["target"].eq(target)
        ].to_numpy(int)
        if len(columns) == 0:
            raise ValueError(f"No guide columns were mapped to {target}.")
        for column in columns:
            carriers = _column_carriers(assignment, int(column))
            n = int(carriers.sum())
            n_ec = int((carriers & ec).sum())
            rows.append(
                {
                    "guide_id": str(guide_map.loc[column, "guide"]),
                    "target": target,
                    "guide_type": "targeting",
                    "n_D3_cells": n,
                    "n_D3_EC_cells": n_ec,
                    "n_D3_nonEC_cells": n - n_ec,
                    "pct_EC_D3": 100.0 * n_ec / n if n else np.nan,
                    "exclusion_rule": "none",
                }
            )

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
        stats_rows.append(
            {
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
            }
        )
    stats = pd.DataFrame(stats_rows)
    stats["FDR_BH"] = bh_adjust(stats["welch_pvalue_two_sided"].to_numpy(float))
    stats["significance"] = stats["FDR_BH"].map(significance_stars)
    stats["BH_number_tests"] = len(stats)
    metadata = {
        "ec_definition": ec_rule,
        "n_D3_cells": n_d3,
        "n_non_targeting_guides_designed": len(nt_columns),
        "n_non_targeting_guides_evaluable": len(control_values),
        "targets": targets,
        "nt_exclusion_targets": exclusion_targets,
    }
    return guide_table, stats, metadata


def plot_figure2d(guide_table: pd.DataFrame, stats: pd.DataFrame):
    """Render Figure 2d with plotnine."""
    order = ["Non-targeting", *FIGURE2D_TARGETS]
    positions = {target: index for index, target in enumerate(order)}
    plot_data = guide_table.loc[guide_table["target"].isin(order)].dropna(subset=["pct_EC_D3"]).copy()
    plot_data["x_center"] = plot_data["target"].map(positions).astype(float)

    jitter_parts: list[pd.DataFrame] = []
    for target in order:
        frame = plot_data.loc[plot_data["target"].eq(target)].copy()
        seed = int.from_bytes(hashlib.sha256(target.encode()).digest()[:4], "little")
        rng = np.random.default_rng(seed)
        spread = 0.10 if target == "Non-targeting" else 0.06
        frame["x"] = frame["x_center"] + np.clip(
            rng.normal(0, spread, len(frame)), -0.18, 0.18
        )
        jitter_parts.append(frame)
    plot_data = pd.concat(jitter_parts, ignore_index=True)

    stats_index = stats.set_index("gene")
    star_rows = [
        {"x": positions[target], "y": 96.0, "label": str(stats_index.loc[target, "significance"])}
        for target in FIGURE2D_TARGETS
        if str(stats_index.loc[target, "significance"])
    ]
    star_df = pd.DataFrame(star_rows, columns=["x", "y", "label"])
    control_mean = float(
        plot_data.loc[plot_data["target"].eq("Non-targeting"), "pct_EC_D3"].mean()
    )
    fill_values = {"Non-targeting": "white", **{target: "#059C48" for target in FIGURE2D_TARGETS}}

    plot = (
        ggplot()
        + geom_hline(yintercept=control_mean, color="#999999", linetype="dashed", size=0.8)
        + geom_boxplot(
            plot_data,
            aes(x="x_center", y="pct_EC_D3", group="target", fill="target"),
            inherit_aes=False,
            width=0.48,
            outlier_alpha=0,
            color="black",
            size=0.65,
        )
        + geom_point(
            plot_data,
            aes(x="x", y="pct_EC_D3"),
            inherit_aes=False,
            color="#666666",
            size=1.6,
            alpha=0.58,
            stroke=0,
        )
        + geom_text(
            star_df,
            aes(x="x", y="y", label="label"),
            inherit_aes=False,
            family="Liberation Sans",
            size=10,
            fontweight="bold",
        )
        + scale_fill_manual(values=fill_values, guide=None)
        + scale_x_continuous(
            breaks=list(positions.values()),
            labels=[f"{target}\nn={int((plot_data['target'] == target).sum())}" for target in order],
            limits=(-0.65, len(order) - 0.35),
        )
        + scale_y_continuous(breaks=np.arange(0, 101, 20), limits=(0, 105))
        + labs(x="gRNAs", y="EC differentiation efficiency (%)")
        + theme_classic(base_family="Liberation Sans", base_size=10)
        + theme(
            figure_size=(4.8, 4.2),
            axis_text_x=element_text(rotation=55, ha="right", size=8),
            axis_text_y=element_text(size=8),
            axis_title=element_text(size=10),
            legend_position="none",
            panel_grid=element_blank(),
            plot_background=element_rect(fill="white", color="white"),
            panel_background=element_rect(fill="white", color="white"),
        )
        + annotate(
            "text", x=-0.55, y=102.5, label="d", family="Liberation Sans", fontweight="bold", size=15
        )
    )
    jitter = plot_data[["target", "guide_id", "pct_EC_D3", "x"]].copy()
    return plot, jitter
