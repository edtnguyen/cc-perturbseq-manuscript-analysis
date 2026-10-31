from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
import hashlib
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
from scipy.stats import ttest_ind

from .figure2_common import FIGURE2D_TARGETS, bh_adjust, significance_stars
from .figure2_ec import _any_carriers, _canonical_target, _column_carriers, _d3_inputs


def compute_figure2d_guide_ec_120_controls(
    adata,
    guide_metadata_path: str | Path,
    targets: Sequence[str] = FIGURE2D_TARGETS,
    *,
    control_exclusion_targets: Sequence[str] = FIGURE2D_TARGETS,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Calculate Figure 2d guide-level EC percentages using all 120 controls.

    The control group combines 60 non-targeting and 60 safe-targeting guides.
    For each individual control guide, D3 cells carrying any guide targeting
    KDR, PLCG1, PRKCD, or ETV2 are excluded before calculating its EC
    percentage. Targeting-guide percentages use all D3 carriers of the
    corresponding guide. Each target-guide distribution is compared with the
    120-guide control distribution by a two-sided Welch t-test, followed by BH
    correction across the four displayed genes.
    """
    targets = [_canonical_target(target) for target in targets]
    exclusion_targets = [
        _canonical_target(target) for target in control_exclusion_targets
    ]
    assignment, ec, guide_map, n_d3, ec_rule = _d3_inputs(
        adata, guide_metadata_path
    )

    control_columns = guide_map.index[guide_map["is_control"]].to_numpy(int)
    nt_columns = guide_map.index[guide_map["is_non_targeting"]].to_numpy(int)
    safe_columns = guide_map.index[guide_map["is_safe_targeting"]].to_numpy(int)
    if (
        len(control_columns) != 120
        or len(nt_columns) != 60
        or len(safe_columns) != 60
    ):
        raise ValueError(
            "Expected 120 controls (60 non-targeting and 60 safe-targeting); "
            f"found total={len(control_columns)}, non-targeting={len(nt_columns)}, "
            f"safe-targeting={len(safe_columns)}."
        )

    exclusion_columns = guide_map.index[
        guide_map["is_targeting"]
        & guide_map["target"].isin(exclusion_targets)
    ].to_numpy(int)
    has_excluded_target = _any_carriers(assignment, exclusion_columns)

    rows: list[dict[str, Any]] = []
    for column in control_columns:
        carriers = (
            _column_carriers(assignment, int(column))
            & ~has_excluded_target
        )
        n_cells = int(carriers.sum())
        n_ec = int((carriers & ec).sum())
        control_type = (
            "non-targeting"
            if bool(guide_map.loc[column, "is_non_targeting"])
            else "safe-targeting"
        )
        rows.append(
            {
                "guide_id": str(guide_map.loc[column, "guide"]),
                "target": "Negative controls",
                "guide_type": "control",
                "control_type": control_type,
                "n_D3_cells": n_cells,
                "n_D3_EC_cells": n_ec,
                "n_D3_nonEC_cells": n_cells - n_ec,
                "pct_EC_D3": (
                    100.0 * n_ec / n_cells if n_cells else np.nan
                ),
                "exclusion_rule": (
                    "exclude control-guide carrier cells with any guide "
                    "targeting " + ", ".join(exclusion_targets)
                ),
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
            n_cells = int(carriers.sum())
            n_ec = int((carriers & ec).sum())
            rows.append(
                {
                    "guide_id": str(guide_map.loc[column, "guide"]),
                    "target": target,
                    "guide_type": "targeting",
                    "control_type": np.nan,
                    "n_D3_cells": n_cells,
                    "n_D3_EC_cells": n_ec,
                    "n_D3_nonEC_cells": n_cells - n_ec,
                    "pct_EC_D3": (
                        100.0 * n_ec / n_cells if n_cells else np.nan
                    ),
                    "exclusion_rule": "none beyond D3 restriction",
                }
            )

    guide_table = pd.DataFrame(rows)
    control_values = guide_table.loc[
        guide_table["guide_type"].eq("control"), "pct_EC_D3"
    ].dropna().to_numpy(float)
    if len(control_values) < 2:
        raise ValueError("Fewer than two evaluable control guides remain.")

    stats_rows: list[dict[str, Any]] = []
    for target in targets:
        target_values = guide_table.loc[
            guide_table["target"].eq(target), "pct_EC_D3"
        ].dropna().to_numpy(float)
        if len(target_values) < 2:
            raise ValueError(f"Fewer than two evaluable guides for {target}.")
        result = ttest_ind(
            target_values,
            control_values,
            equal_var=False,
            nan_policy="raise",
        )
        mean_target = float(target_values.mean())
        mean_control = float(control_values.mean())
        stats_rows.append(
            {
                "gene": target,
                "n_targeting_guides": len(target_values),
                "n_control_guides": len(control_values),
                "n_non_targeting_guides_designed": len(nt_columns),
                "n_safe_targeting_guides_designed": len(safe_columns),
                "mean_target_EC_pct": mean_target,
                "median_target_EC_pct": float(np.median(target_values)),
                "sd_target_EC_pct": float(np.std(target_values, ddof=1)),
                "mean_control_EC_pct": mean_control,
                "median_control_EC_pct": float(np.median(control_values)),
                "sd_control_EC_pct": float(np.std(control_values, ddof=1)),
                "difference_target_minus_control_pp": (
                    mean_target - mean_control
                ),
                "effect_log2_target_over_control": float(
                    np.log2(mean_target / mean_control)
                ),
                "welch_t_statistic": float(result.statistic),
                "welch_df": float(result.df),
                "welch_pvalue_two_sided": float(result.pvalue),
            }
        )

    stats = pd.DataFrame(stats_rows)
    stats["FDR_BH"] = bh_adjust(
        stats["welch_pvalue_two_sided"].to_numpy(float)
    )
    stats["significance"] = stats["FDR_BH"].map(significance_stars)
    stats["BH_number_tests"] = len(stats)
    metadata = {
        "ec_definition": ec_rule,
        "n_D3_cells": n_d3,
        "n_control_guides_designed": len(control_columns),
        "n_control_guides_evaluable": len(control_values),
        "n_non_targeting_guides_designed": len(nt_columns),
        "n_safe_targeting_guides_designed": len(safe_columns),
        "targets": targets,
        "control_exclusion_targets": exclusion_targets,
    }
    return guide_table, stats, metadata


def plot_figure2d_120_controls(
    guide_table: pd.DataFrame,
    stats: pd.DataFrame,
):
    """Render Figure 2d with plotnine and a 120-guide control box."""
    order = ["Negative controls", *FIGURE2D_TARGETS]
    positions = {target: index for index, target in enumerate(order)}
    plot_data = guide_table.loc[
        guide_table["target"].isin(order)
    ].dropna(subset=["pct_EC_D3"]).copy()
    plot_data["x_center"] = plot_data["target"].map(positions).astype(float)

    jitter_parts: list[pd.DataFrame] = []
    for target in order:
        frame = plot_data.loc[plot_data["target"].eq(target)].copy()
        seed = int.from_bytes(
            hashlib.sha256(target.encode()).digest()[:4], "little"
        )
        rng = np.random.default_rng(seed)
        spread = 0.10 if target == "Negative controls" else 0.06
        frame["x"] = frame["x_center"] + np.clip(
            rng.normal(0, spread, len(frame)), -0.18, 0.18
        )
        jitter_parts.append(frame)
    plot_data = pd.concat(jitter_parts, ignore_index=True)

    stats_index = stats.set_index("gene")
    star_rows = [
        {
            "x": positions[target],
            "y": 96.0,
            "label": str(stats_index.loc[target, "significance"]),
        }
        for target in FIGURE2D_TARGETS
        if str(stats_index.loc[target, "significance"])
    ]
    star_df = pd.DataFrame(star_rows, columns=["x", "y", "label"])
    control_mean = float(
        plot_data.loc[
            plot_data["target"].eq("Negative controls"), "pct_EC_D3"
        ].mean()
    )
    fill_values = {
        "Negative controls": "white",
        **{target: "#059C48" for target in FIGURE2D_TARGETS},
    }

    plot = (
        ggplot()
        + geom_hline(
            yintercept=control_mean,
            color="#999999",
            linetype="dashed",
            size=0.8,
        )
        + geom_boxplot(
            plot_data,
            aes(
                x="x_center",
                y="pct_EC_D3",
                group="target",
                fill="target",
            ),
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
            labels=[
                f"{target}\nn={int((plot_data['target'] == target).sum())}"
                for target in order
            ],
            limits=(-0.65, len(order) - 0.35),
        )
        + scale_y_continuous(
            breaks=np.arange(0, 101, 20),
            limits=(0, 105),
        )
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
            "text",
            x=-0.55,
            y=102.5,
            label="d",
            family="Liberation Sans",
            fontweight="bold",
            size=15,
        )
    )
    jitter = plot_data[[
        "target", "guide_id", "control_type", "pct_EC_D3", "x"
    ]].copy()
    return plot, jitter
