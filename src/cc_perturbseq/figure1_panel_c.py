from __future__ import annotations

from pathlib import Path
import hashlib
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import ttest_ind


# -----------------------------------------------------------------------------
# Figure 1c: validation-screen guide-level Day 3 EC percentages
#
# The plot is rebuilt from Chloe McCreery's validation-screen guide-count table:
#   /oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/
#   MageckAnalysis/mageck_sgRNA_readcount.txt
#
# For each biological replicate, guide enrichment in the CDH5-positive versus
# CDH5-negative fraction is converted to a relative log2 odds shift. The shift
# is centered so that the arithmetic mean reconstructed EC percentage across
# the 120 control guides is 70%. Replicate-specific EC percentages are averaged
# per guide. Each target gene is compared with the 120 control guides using a
# two-sided Welch t-test; P values are BH-adjusted across all 297 target genes.
# -----------------------------------------------------------------------------

PANEL_C_VALUE_COL = "estimated_EC_pct_control70"
PANEL_C_CONTROL_TYPES = {"NON-TARGETING", "SAFE"}
PANEL_C_GROUPS = [
    ("Control", ["NTC"]),
    ("BMP/TGFβ pathway", ["BMPR1A", "SMAD1", "SMAD5"]),
    ("Mesoderm regulators", ["MIXL1", "EOMES"]),
    ("EC + LM regulators", ["KDR", "KCTD10", "PRKCD", "ETV2", "PLCG1"]),
    ("Hippo pathway", ["TAOK1", "AMOTL2", "LATS2", "NF2"]),
]
PANEL_C_GROUP_LABELS = {
    "Control": "Control",
    "BMP/TGFβ pathway": "BMP/TGFβ\npathway",
    "Mesoderm regulators": "Mesoderm\nregulators",
    "EC + LM regulators": "EC + LM\nregulators",
    "Hippo pathway": "Hippo\npathway",
}
PANEL_C_GROUP_COLORS = {
    "Control": "white",
    "BMP/TGFβ pathway": "#7A43A2",
    "Mesoderm regulators": "#D97A34",
    "EC + LM regulators": "#059C48",
    "Hippo pathway": "#ED4A97",
}


def sha256_file(path: str | Path) -> str:
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def panel_c_bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or not np.all(np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("BH adjustment requires finite p-values in [0, 1].")
    order = np.argsort(p, kind="stable")
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.empty(len(p), dtype=float)
    q[order] = np.clip(np.minimum.accumulate(ranked[::-1])[::-1], 0, 1)
    return q


def panel_c_stars(q: float) -> str:
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""


def _read_count_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".csv":
        table = pd.read_csv(path)
    else:
        try:
            table = pd.read_csv(path, sep="\t")
        except Exception:
            table = pd.read_csv(path, sep=None, engine="python")
    table.columns = [str(column).strip() for column in table.columns]
    return table


def _infer_column(table: pd.DataFrame, candidates: list[str]) -> str:
    normalized = {str(column).strip().lower(): str(column) for column in table.columns}
    for candidate in candidates:
        if candidate.lower() in normalized:
            return normalized[candidate.lower()]
    raise ValueError(
        f"Could not identify any of {candidates}. Available columns: {list(table.columns)}"
    )


def _guide_type(guide_id: str, target: str) -> str:
    text = f"{guide_id} {target}".upper().replace("_", "-")
    if re.search(r"NON[- ]?TARGET", text):
        return "NON-TARGETING"
    if re.search(r"SAFE[- ]?TARGET", text) or target.strip().upper() == "SAFE":
        return "SAFE"
    return "TARGETING"


def load_validation_screen_counts(path: str | Path) -> pd.DataFrame:
    """Load and validate Chloe's validation-screen guide-count table."""
    raw = _read_count_table(path)
    guide_col = _infer_column(raw, ["sgRNA", "guide", "guide_id", "grna"])
    target_col = _infer_column(raw, ["Gene", "target", "target_gene", "gene_symbol"])
    required_counts = [
        "D3PR1", "D3PR2", "D0R1", "D0R2",
        "D3NR1", "D3NR2", "D3TOTAL1", "D3TOTAL2",
    ]
    missing = [column for column in required_counts if column not in raw.columns]
    if missing:
        raise ValueError(
            "Validation-screen count table is missing required columns: "
            + ", ".join(missing)
        )

    table = pd.DataFrame(
        {
            "guide_id": raw[guide_col].astype(str).str.strip(),
            "target_gene": raw[target_col].astype(str).str.strip(),
        }
    )
    for column in required_counts:
        table[column] = pd.to_numeric(raw[column], errors="coerce")
    if table["guide_id"].eq("").any() or table["guide_id"].duplicated().any():
        raise ValueError("Guide identifiers must be nonempty and unique.")
    if table[required_counts].isna().any().any():
        bad = int(table[required_counts].isna().sum().sum())
        raise ValueError(f"Guide-count table contains {bad} nonnumeric count values.")
    if (table[required_counts] < 0).any().any():
        raise ValueError("Guide counts must be nonnegative.")

    table["guide_type"] = [
        _guide_type(guide, target)
        for guide, target in zip(table["guide_id"], table["target_gene"])
    ]
    table.loc[table["guide_type"].eq("NON-TARGETING"), "target_gene"] = "NON-TARGETING"
    table.loc[table["guide_type"].eq("SAFE"), "target_gene"] = "SAFE"
    table.loc[table["guide_type"].eq("TARGETING"), "target_gene"] = (
        table.loc[table["guide_type"].eq("TARGETING"), "target_gene"]
        .astype(str).str.strip().str.upper()
    )

    counts_by_type = table["guide_type"].value_counts().to_dict()
    expected = {"TARGETING": 1806, "NON-TARGETING": 60, "SAFE": 60}
    if counts_by_type != expected:
        raise ValueError(
            f"Unexpected guide-type counts: {counts_by_type}; expected {expected}."
        )
    n_genes = table.loc[table["guide_type"].eq("TARGETING"), "target_gene"].nunique()
    if n_genes != 297:
        raise ValueError(f"Expected 297 targeting genes; found {n_genes}.")
    return table


def _ec_probability_from_log2_odds(log2_odds: np.ndarray) -> np.ndarray:
    odds = np.exp2(np.asarray(log2_odds, dtype=float))
    return odds / (1.0 + odds)


def _calibrate_replicate_to_control_mean(
    raw_log2_positive_over_negative: np.ndarray,
    control_mask: np.ndarray,
    baseline_fraction: float,
) -> tuple[np.ndarray, float]:
    if not 0 < baseline_fraction < 1:
        raise ValueError("The control baseline fraction must lie strictly between 0 and 1.")
    baseline_log2_odds = np.log2(baseline_fraction / (1.0 - baseline_fraction))

    def objective(offset: float) -> float:
        probabilities = _ec_probability_from_log2_odds(
            baseline_log2_odds + raw_log2_positive_over_negative - offset
        )
        return float(probabilities[control_mask].mean() - baseline_fraction)

    offset = float(brentq(objective, -100.0, 100.0))
    centered_log2_odds = raw_log2_positive_over_negative - offset
    reconstructed = _ec_probability_from_log2_odds(
        baseline_log2_odds + centered_log2_odds
    )
    return reconstructed, offset


def reconstruct_validation_guide_ec_percentages(
    count_table: pd.DataFrame,
    *,
    control_baseline_pct: float = 70.0,
    pseudocount: float = 0.5,
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Reconstruct guide-level EC percentages from D3 positive/negative counts."""
    if pseudocount <= 0:
        raise ValueError("pseudocount must be positive.")
    baseline_fraction = float(control_baseline_pct) / 100.0
    control_mask = count_table["guide_type"].isin(PANEL_C_CONTROL_TYPES).to_numpy()
    if int(control_mask.sum()) != 120:
        raise ValueError(f"Expected 120 controls; found {int(control_mask.sum())}.")

    output = count_table.copy()
    offsets: dict[str, float] = {}
    for replicate in (1, 2):
        positive = output[f"D3PR{replicate}"].to_numpy(float)
        negative = output[f"D3NR{replicate}"].to_numpy(float)
        raw_log2_odds = np.log2(
            (positive + pseudocount) / (negative + pseudocount)
        )
        reconstructed, offset = _calibrate_replicate_to_control_mean(
            raw_log2_odds,
            control_mask,
            baseline_fraction,
        )
        output[f"log2_odds_vs_control_R{replicate}"] = raw_log2_odds - offset
        output[f"estimated_EC_pct_R{replicate}_control70"] = 100.0 * reconstructed
        offsets[f"replicate_{replicate}_log2_odds_offset"] = offset

    output[PANEL_C_VALUE_COL] = output[
        ["estimated_EC_pct_R1_control70", "estimated_EC_pct_R2_control70"]
    ].mean(axis=1)
    output["EC_change_vs_control_percentage_points"] = (
        output[PANEL_C_VALUE_COL] - control_baseline_pct
    )
    output["replicate_absolute_difference_percentage_points"] = (
        output["estimated_EC_pct_R1_control70"]
        - output["estimated_EC_pct_R2_control70"]
    ).abs()
    delta_1 = output["estimated_EC_pct_R1_control70"] - control_baseline_pct
    delta_2 = output["estimated_EC_pct_R2_control70"] - control_baseline_pct
    output["replicate_direction_concordant"] = (delta_1 * delta_2) > 0
    output["replicate_difference_gt20pp"] = (
        output["replicate_absolute_difference_percentage_points"] > 20
    )
    output["low_D3_representation_lt100_in_either_rep"] = (
        output[["D3TOTAL1", "D3TOTAL2"]].min(axis=1) < 100
    )
    output["low_D0_representation_lt20_in_either_rep"] = (
        output[["D0R1", "D0R2"]].min(axis=1) < 20
    )

    controls = output.loc[control_mask]
    for replicate in (1, 2):
        value = float(controls[f"estimated_EC_pct_R{replicate}_control70"].mean())
        if not np.isclose(value, control_baseline_pct, atol=1e-10, rtol=0):
            raise AssertionError(
                f"Replicate {replicate} control mean is {value}, not {control_baseline_pct}."
            )
    if not np.isclose(
        controls[PANEL_C_VALUE_COL].mean(), control_baseline_pct, atol=1e-10, rtol=0
    ):
        raise AssertionError("Combined control-guide EC mean does not equal the baseline.")
    return output, offsets


def split_panel_c_guides(table: pd.DataFrame):
    controls = table.loc[table["guide_type"].isin(PANEL_C_CONTROL_TYPES)].copy()
    targeting = {
        gene: group.sort_values("guide_id").copy()
        for gene, group in table.loc[table["guide_type"].eq("TARGETING")].groupby(
            "target_gene", sort=True
        )
    }
    for _, genes in PANEL_C_GROUPS[1:]:
        for gene in genes:
            if gene not in targeting:
                raise ValueError(f"Requested Figure 1c gene is missing: {gene}")
    return controls.sort_values("guide_id"), targeting


def compute_panel_c_welch_statistics(
    controls: pd.DataFrame,
    targeting: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Run two-sided Welch t-tests and BH correction across all target genes."""
    control_values = controls[PANEL_C_VALUE_COL].to_numpy(float)
    control_mean = float(control_values.mean())
    rows: list[dict[str, object]] = []
    for gene, group in targeting.items():
        target_values = group[PANEL_C_VALUE_COL].to_numpy(float)
        result = ttest_ind(
            target_values,
            control_values,
            equal_var=False,
            nan_policy="raise",
        )
        pvalue = float(result.pvalue)
        if not np.isfinite(pvalue):
            raise ValueError(f"Welch t-test returned a nonfinite P value for {gene}.")
        target_mean = float(target_values.mean())
        rows.append(
            {
                "gene": gene,
                "n_targeting_guides": len(target_values),
                "n_control_guides": len(control_values),
                "mean_target_EC_pct": target_mean,
                "median_target_EC_pct": float(np.median(target_values)),
                "sd_target_EC_pct": float(np.std(target_values, ddof=1)),
                "mean_control_EC_pct": control_mean,
                "median_control_EC_pct": float(np.median(control_values)),
                "sd_control_EC_pct": float(np.std(control_values, ddof=1)),
                "difference_target_minus_control_pp": target_mean - control_mean,
                "effect_log2_target_over_control": float(
                    np.log2(target_mean / control_mean)
                ),
                "welch_t_statistic": float(result.statistic),
                "welch_pvalue_two_sided": pvalue,
            }
        )

    stats = pd.DataFrame(rows).sort_values("gene").reset_index(drop=True)
    if len(stats) != 297:
        raise ValueError(f"Expected 297 target-gene tests; found {len(stats)}.")
    stats["FDR_BH"] = panel_c_bh_adjust(
        stats["welch_pvalue_two_sided"].to_numpy()
    )
    stats["significance"] = stats["FDR_BH"].map(panel_c_stars)
    stats["BH_number_genes"] = len(stats)
    stats["test"] = "two-sided Welch t-test"
    return stats


def plot_panel_c(
    controls: pd.DataFrame,
    targeting: dict[str, pd.DataFrame],
    stats: pd.DataFrame,
    *,
    control_baseline_pct: float = 70.0,
):
    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "DejaVu Sans",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.transparent": False,
        }
    )
    stats_index = stats.set_index("gene")
    datasets: list[np.ndarray] = []
    labels: list[str] = []
    positions: list[float] = []
    spans: list[tuple[str, float, float]] = []
    names: list[str] = []
    group_by_gene: dict[str, str] = {"NTC": "Control"}
    for group_name, genes in PANEL_C_GROUPS[1:]:
        for gene in genes:
            group_by_gene[gene] = group_name

    jitter_rows: list[dict[str, object]] = []
    x = 0.0
    for group_name, genes in PANEL_C_GROUPS:
        start = x
        for gene in genes:
            group = controls if gene == "NTC" else targeting[gene]
            values = group[PANEL_C_VALUE_COL].to_numpy(float)
            positions.append(x)
            datasets.append(values)
            names.append(gene)
            labels.append(f"{gene}\nn = {len(values)}")
            x += 1.0
        spans.append((group_name, start, x - 1.0))
        x += 0.90

    x_min, x_max = -0.6, positions[-1] + 0.6
    width_per_position = 0.58
    ax_width = (x_max - x_min) * width_per_position
    ax_height, left, right, bottom, top = 3.85, 1.15, 0.15, 1.80, 0.62
    fig_w, fig_h = left + ax_width + right, bottom + ax_height + top
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes(
        [left / fig_w, bottom / fig_h, ax_width / fig_w, ax_height / fig_h]
    )
    ax.axhline(
        control_baseline_pct,
        color="#9EA3A8",
        linestyle="--",
        linewidth=1.6,
        zorder=1,
    )
    boxplot = ax.boxplot(
        datasets,
        positions=positions,
        widths=0.42,
        whis=1.5,
        patch_artist=True,
        showfliers=False,
        manage_ticks=False,
        medianprops={"color": "black", "linewidth": 2.0, "zorder": 5},
        whiskerprops={"color": "black", "linewidth": 1.8, "zorder": 5},
        capprops={"color": "black", "linewidth": 1.8, "zorder": 5},
        boxprops={"edgecolor": "black", "linewidth": 1.8, "zorder": 3},
    )
    for box, gene in zip(boxplot["boxes"], names):
        box.set_facecolor(PANEL_C_GROUP_COLORS[group_by_gene[gene]])
        box.set_alpha(0.86 if gene != "NTC" else 1.0)

    for gene, position, values in zip(names, positions, datasets):
        digest = int.from_bytes(hashlib.sha256(gene.encode()).digest()[:4], "little")
        rng = np.random.default_rng(np.random.SeedSequence([7, digest]))
        jitter = np.clip(
            rng.normal(0, 0.09 if gene == "NTC" else 0.055, len(values)),
            -0.16,
            0.16,
        )
        ax.scatter(
            np.full(len(values), position) + jitter,
            values,
            s=16 if gene == "NTC" else 18,
            c="#7A7A7A",
            alpha=0.48,
            linewidths=0,
            edgecolors="none",
            zorder=6,
            rasterized=False,
        )
        group = controls if gene == "NTC" else targeting[gene]
        for (_, guide), offset in zip(group.iterrows(), jitter):
            jitter_rows.append(
                {
                    "plot_gene": gene,
                    "guide_id": guide["guide_id"],
                    PANEL_C_VALUE_COL: guide[PANEL_C_VALUE_COL],
                    "x_center": position,
                    "x_jitter": float(offset),
                }
            )
        if gene != "NTC":
            annotation = str(stats_index.loc[gene, "significance"])
            if annotation:
                ax.text(
                    position,
                    95,
                    annotation,
                    ha="center",
                    va="bottom",
                    fontsize=12,
                    fontweight="bold",
                    zorder=8,
                )

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, 105)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_xticks(positions)
    ax.set_xticklabels(
        labels,
        rotation=58,
        ha="right",
        rotation_mode="anchor",
        fontsize=12,
    )
    ax.set_ylabel("Estimated D3 EC cells (%)", fontsize=18, labelpad=8)
    ax.tick_params(axis="y", labelsize=12, width=2.0, length=6)
    ax.tick_params(axis="x", width=2.0, length=6, pad=6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(2.2)

    transform = ax.get_xaxis_transform()
    for group_name, start, end in spans:
        ax.plot(
            [start - 0.32, end + 0.32],
            [-0.290, -0.290],
            transform=transform,
            clip_on=False,
            color="black",
            linewidth=1,
        )
        ax.text(
            (start + end) / 2,
            -0.323,
            PANEL_C_GROUP_LABELS[group_name],
            transform=transform,
            ha="center",
            va="top",
            fontsize=10.5,
            linespacing=1.15,
            fontweight="bold",
            clip_on=False,
        )

    fig.text(
        left / fig_w,
        (bottom + ax_height + 0.24) / fig_h,
        "Two-sided Welch t-test; BH correction across all 297 target genes",
        ha="left",
        va="bottom",
        fontsize=9.5,
    )
    fig.text(
        (left + ax_width) / fig_w,
        (bottom + ax_height + 0.24) / fig_h,
        "* FDR < 0.05   ** < 0.01   *** < 0.001",
        ha="right",
        va="bottom",
        fontsize=9.5,
    )
    return fig, pd.DataFrame(jitter_rows)


def _save_figure(fig, outdir: Path, stem: str) -> dict[str, str]:
    outputs: dict[str, str] = {}
    for suffix in ("pdf", "svg", "png"):
        path = outdir / f"{stem}.{suffix}"
        fig.savefig(
            path,
            dpi=600 if suffix == "png" else None,
            bbox_inches="tight",
            pad_inches=0.03,
            facecolor="white",
        )
        outputs[suffix] = str(path)
    return outputs


def run_figure_01c(
    validation_count_path: str | Path,
    outdir: str | Path,
    *,
    control_baseline_pct: float = 70.0,
    pseudocount: float = 0.5,
):
    """Recalculate Figure 1c from Chloe's validation-screen count table."""
    validation_count_path = Path(validation_count_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    counts = load_validation_screen_counts(validation_count_path)
    guide_table, offsets = reconstruct_validation_guide_ec_percentages(
        counts,
        control_baseline_pct=control_baseline_pct,
        pseudocount=pseudocount,
    )
    controls, targeting = split_panel_c_guides(guide_table)
    stats = compute_panel_c_welch_statistics(controls, targeting)
    selected_genes = [gene for _, genes in PANEL_C_GROUPS[1:] for gene in genes]
    selected_stats = stats.set_index("gene").loc[selected_genes].reset_index()

    fig, jitter = plot_panel_c(
        controls,
        targeting,
        stats,
        control_baseline_pct=control_baseline_pct,
    )
    outputs = _save_figure(fig, outdir, "Figure_01c_validation_screen_EC_fraction")

    guide_path = outdir / "Figure_01c_reconstructed_guide_EC_percentages.csv"
    stats_path = outdir / "Figure_01c_Welch_statistics_all_297_genes.csv"
    selected_path = outdir / "Figure_01c_Welch_statistics_plotted_genes.csv"
    jitter_path = outdir / "Figure_01c_plot_jitter_coordinates.csv"
    guide_table.to_csv(guide_path, index=False)
    stats.to_csv(stats_path, index=False)
    selected_stats.to_csv(selected_path, index=False)
    jitter.to_csv(jitter_path, index=False)
    outputs.update(
        {
            "guide_level_reconstruction": str(guide_path),
            "all_gene_statistics": str(stats_path),
            "plotted_gene_statistics": str(selected_path),
            "jitter_coordinates": str(jitter_path),
        }
    )
    return (
        fig,
        stats,
        selected_stats,
        guide_table,
        outputs,
        sha256_file(validation_count_path),
        offsets,
    )
