from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
import re
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from .program_io import infer_column, read_table, significance_stars


PATHWAY_COLORS = {
    "BMP/TGFβ pathway": "#7A43A2",
    "Mesoderm regulators": "#D97A34",
    "LM + EC regulators": "#059C48",
    "Hippo pathway": "#ED4A97",
}

FIGURE1_TARGET_GROUPS: dict[str, list[str]] = {
    "BMP/TGFβ pathway": ["BMPR1A", "SMAD1", "SMAD5"],
    "Mesoderm regulators": ["MIXL1", "EOMES"],
    "LM + EC regulators": ["KDR", "KCTD10", "PRKCD", "ETV2", "PLCG1"],
    "Hippo pathway": ["TAOK1", "AMOTL2", "LATS2", "NF2"],
}


def _clean_gene(value: Any) -> str:
    return re.sub(r"(?:_-.*|-P1P2)$", "", str(value).strip())


def load_genomewide_screen(
    mageck_path: str | Path,
    utest_path: str | Path,
) -> pd.DataFrame:
    """Load the genome-wide screen using the source-script sign convention.

    The x-axis is the U-test log2 fold change. The MaGECK FDR displayed on the
    y-axis is selected from the negative or positive direction according to
    the sign of that fold change.
    """
    mageck = read_table(mageck_path)
    utest = read_table(utest_path)

    mageck_id = infer_column(mageck, ["id", "target", "gene"])
    utest_id = infer_column(utest, ["target", "id", "gene"])
    utest_effect = infer_column(utest, ["log2FC", "log2_fc", "effect"])
    neg_fdr_col = infer_column(mageck, ["neg|fdr", "negative_fdr", "neg_fdr"])
    pos_fdr_col = infer_column(mageck, ["pos|fdr", "positive_fdr", "pos_fdr"])

    merged = mageck.merge(
        utest,
        left_on=mageck_id,
        right_on=utest_id,
        how="outer",
        suffixes=("_mageck", "_utest"),
    )
    raw_id = merged[mageck_id].where(merged[mageck_id].notna(), merged[utest_id])
    merged["gene"] = raw_id.map(_clean_gene)
    merged["effect"] = pd.to_numeric(merged[utest_effect], errors="coerce")
    neg_fdr = pd.to_numeric(merged[neg_fdr_col], errors="coerce")
    pos_fdr = pd.to_numeric(merged[pos_fdr_col], errors="coerce")
    merged["mageck_fdr"] = np.where(merged["effect"] < 0, neg_fdr, pos_fdr)
    merged["mageck_fdr"] = pd.to_numeric(merged["mageck_fdr"], errors="coerce")
    merged["neg_fdr"] = neg_fdr
    merged["pos_fdr"] = pos_fdr
    merged["is_hit"] = (neg_fdr < 0.05) | (pos_fdr < 0.05)
    return merged.dropna(subset=["gene", "effect", "mageck_fdr"]).copy()


def plot_genomewide_volcano(
    screen: pd.DataFrame,
    *,
    target_groups: Mapping[str, Sequence[str]] = FIGURE1_TARGET_GROUPS,
    fdr_threshold: float = 0.05,
    effect_reference: float = 0.5,
):
    plot = screen.copy()
    plot["minus_log10_fdr"] = -np.log10(
        plot["mageck_fdr"].clip(lower=np.nextafter(0, 1))
    )
    color_by_gene = {
        gene: PATHWAY_COLORS[group]
        for group, genes in target_groups.items()
        for gene in genes
    }
    highlighted = plot[plot.gene.isin(color_by_gene)].copy()
    background = plot[~plot.is_hit & ~plot.gene.isin(color_by_gene)]
    hits = plot[plot.is_hit & ~plot.gene.isin(color_by_gene)]

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.scatter(
        background.effect,
        background.minus_log10_fdr,
        s=5,
        alpha=0.45,
        color="0.82",
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        hits.effect,
        hits.minus_log10_fdr,
        s=7,
        alpha=0.8,
        color="0.08",
        linewidths=0,
        rasterized=True,
    )
    for gene, group in highlighted.groupby("gene", sort=False):
        color = color_by_gene[gene]
        ax.scatter(
            group.effect,
            group.minus_log10_fdr,
            s=28,
            color=color,
            linewidths=0,
            zorder=4,
        )
        row = group.sort_values("mageck_fdr").iloc[0]
        ax.annotate(
            gene,
            (row.effect, row.minus_log10_fdr),
            xytext=(3, 4),
            textcoords="offset points",
            fontsize=6.5,
            fontstyle="italic",
            color=color,
        )

    ax.axhline(-np.log10(fdr_threshold), color="0.6", linestyle="--", linewidth=0.8)
    ax.axvline(-effect_reference, color="#E99A9A", linestyle=":", linewidth=0.8)
    ax.axvline(effect_reference, color="#E99A9A", linestyle=":", linewidth=0.8)
    ax.set_xlabel("log₂ guide enrichment in CD144⁺ vs CD144⁻")
    ax.set_ylabel("−log₁₀ MaGECK FDR")
    ax.text(0.02, 0.04, "Enriched in\nCD144⁻", transform=ax.transAxes, fontsize=7)
    ax.text(0.98, 0.04, "Enriched in\nCD144⁺", transform=ax.transAxes, ha="right", fontsize=7)
    return fig, plot


def standardize_validation_summary(
    summary_path: str | Path,
    gene_types_path: str | Path | None = None,
) -> pd.DataFrame:
    table = read_table(summary_path)
    gene_col = infer_column(table, ["Gene", "Target", "gene_symbol"])
    output = pd.DataFrame({"gene": table[gene_col].astype(str)})
    candidates = {
        "replicate_1": ["U.Test.FC.BR1", "BR1", "replicate_1", "biorep1"],
        "replicate_2": ["U.Test.FC.BR2", "BR2", "replicate_2", "biorep2"],
        "validation_effect": ["Average.U.Test.FC", "validation_effect", "average_fc"],
        "genomewide_effect": ["Original.FC.Mageck", "original_mageck", "genomewide_effect"],
        "validation_fdr": ["FDR.corrected", "FDR", "validation_fdr", "padj"],
    }
    for name, options in candidates.items():
        col = infer_column(table, options, required=False)
        output[name] = pd.to_numeric(table[col], errors="coerce") if col else np.nan

    if gene_types_path is not None and Path(gene_types_path).exists():
        types = read_table(gene_types_path)
        type_gene = infer_column(types, ["Target", "Gene", "gene"])
        type_col = infer_column(types, ["Type", "gene_type", "category"])
        mapping = types[[type_gene, type_col]].rename(
            columns={type_gene: "gene", type_col: "gene_type"}
        )
        mapping["gene"] = mapping.gene.astype(str)
        output = output.merge(mapping, on="gene", how="left")
    else:
        output["gene_type"] = np.nan
    return output


def _identity_limits(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    values = pd.concat([x, y]).replace([np.inf, -np.inf], np.nan).dropna()
    if values.empty:
        return -1.0, 1.0
    pad = max(0.1, 0.05 * (values.max() - values.min()))
    return float(values.min() - pad), float(values.max() + pad)


def plot_validation_replicates(table: pd.DataFrame):
    data = table.dropna(subset=["replicate_1", "replicate_2"]).copy()
    r = data[["replicate_1", "replicate_2"]].corr().iloc[0, 1] if len(data) > 1 else np.nan
    fig, ax = plt.subplots(figsize=(4.1, 4.0))
    ax.scatter(data.replicate_1, data.replicate_2, s=7, color="0.15", alpha=0.55, rasterized=True)
    lower, upper = _identity_limits(data.replicate_1, data.replicate_2)
    ax.plot([lower, upper], [lower, upper], color="0.5", linestyle="--", linewidth=0.8)
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_xlabel("Validation-screen bioreplicate 1\nlog₂ enrichment")
    ax.set_ylabel("Validation-screen bioreplicate 2\nlog₂ enrichment")
    ax.text(0.04, 0.95, f"Pearson r = {r:.2f}\nn = {len(data):,}", transform=ax.transAxes, va="top", fontsize=7)
    return fig, data, r


def plot_genomewide_validation(table: pd.DataFrame):
    data = table.dropna(subset=["genomewide_effect", "validation_effect"]).copy()
    r = data[["genomewide_effect", "validation_effect"]].corr().iloc[0, 1] if len(data) > 1 else np.nan
    is_hit = data.gene_type.astype(str).str.contains("genome_wide_hit", case=False, na=False)
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    ax.scatter(
        data.loc[~is_hit, "genomewide_effect"],
        data.loc[~is_hit, "validation_effect"],
        s=7,
        color="0.65",
        alpha=0.45,
        linewidths=0,
        rasterized=True,
    )
    ax.scatter(
        data.loc[is_hit, "genomewide_effect"],
        data.loc[is_hit, "validation_effect"],
        s=9,
        color="0.05",
        alpha=0.7,
        linewidths=0,
        rasterized=True,
    )
    lower, upper = _identity_limits(data.genomewide_effect, data.validation_effect)
    ax.axhline(0, color="0.75", linewidth=0.6)
    ax.axvline(0, color="0.75", linewidth=0.6)
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_xlabel("Genome-wide screen\nlog₂ enrichment")
    ax.set_ylabel("Validation screen\nlog₂ enrichment")
    ax.text(0.04, 0.95, f"Pearson r = {r:.2f}\nn = {len(data):,}", transform=ax.transAxes, va="top", fontsize=7)
    return fig, data, r


def standardize_ec_fraction_table(path: str | Path) -> pd.DataFrame:
    table = read_table(path)
    guide_col = infer_column(table, ["guide", "guide_id", "sgRNA", "grna"])
    target_col = infer_column(table, ["target", "gene", "gene_symbol", "target_gene"])
    ec_col = infer_column(
        table,
        ["D3 EC percentage", "D3_EC_percentage", "EC_percent", "EC%", "ec_fraction", "ec_percentage"],
    )
    fdr_col = infer_column(
        table,
        ["BH_FDR", "FDR", "q_value", "qvalue", "adjusted_p_value", "padj"],
        required=False,
    )
    control_col = infer_column(table, ["is_control", "control", "guide_type", "type"], required=False)
    output = pd.DataFrame(
        {
            "guide": table[guide_col].astype(str),
            "target": table[target_col].astype(str),
            "ec_percent": pd.to_numeric(table[ec_col], errors="coerce"),
            "fdr": pd.to_numeric(table[fdr_col], errors="coerce") if fdr_col else np.nan,
        }
    )
    if output.ec_percent.dropna().max() <= 1.5:
        output["ec_percent"] *= 100.0
    control_pattern = r"non[-_ ]?target|safe[-_ ]?target|\bntc\b|negative|control"
    output["is_control"] = output.target.str.contains(control_pattern, case=False, regex=True)
    output["is_control"] |= output.guide.str.contains(control_pattern, case=False, regex=True)
    if control_col:
        control_value = table[control_col].astype(str).str.lower()
        output["is_control"] |= control_value.isin(["true", "1", "control", "ntc", "safe-targeting", "non-targeting"])
    output.loc[output.is_control, "target"] = "Control"
    return output.dropna(subset=["ec_percent"]).copy()


def plot_guide_ec_fraction(
    table: pd.DataFrame,
    *,
    targets: Sequence[str],
    target_groups: Mapping[str, Sequence[str]] | None = None,
    control_label: str = "Control",
    seed: int = 11,
):
    order = [control_label, *targets]
    data = table[table.target.isin(order)].copy()
    if data.empty:
        raise ValueError(f"None of the requested targets were present: {order}")
    rng = np.random.default_rng(seed)
    fig, ax = plt.subplots(figsize=(max(4.8, 0.55 * len(order)), 4.1))
    colors: dict[str, str] = {control_label: "0.55"}
    if target_groups:
        for group, genes in target_groups.items():
            for gene in genes:
                colors[gene] = PATHWAY_COLORS.get(group, "0.25")
    for index, target in enumerate(order):
        values = data.loc[data.target == target, "ec_percent"].to_numpy(float)
        if not len(values):
            continue
        ax.boxplot(
            values,
            positions=[index],
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            boxprops={"facecolor": colors.get(target, "0.35"), "alpha": 0.8, "linewidth": 0.7},
            medianprops={"color": "black", "linewidth": 0.9},
            whiskerprops={"linewidth": 0.7},
            capprops={"linewidth": 0.7},
        )
        jitter = rng.uniform(-0.16, 0.16, size=len(values))
        ax.scatter(
            index + jitter,
            values,
            s=12,
            facecolors="white",
            edgecolors=colors.get(target, "0.2"),
            linewidths=0.7,
            zorder=3,
        )
        if target != control_label:
            q = data.loc[data.target == target, "fdr"].dropna()
            if not q.empty:
                star = significance_stars(float(q.min()))
                if star:
                    ax.text(index, min(99, values.max() + 5), star, ha="center", va="bottom", fontsize=7)
    control_median = data.loc[data.target == control_label, "ec_percent"].median()
    if np.isfinite(control_median):
        ax.axhline(control_median, color="0.55", linestyle="--", linewidth=0.8)
    ax.set_xticks(range(len(order)), order, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Day 3 EC cells (%)")
    ax.set_ylim(bottom=0)
    return fig, data


PKC_GROUPS: dict[str, list[str]] = {
    "Conventional PKC": ["PRKCA", "PRKCB", "PRKCG"],
    "Novel PKC": ["PRKCD", "PRKCE", "PRKCH", "PRKCQ"],
    "Atypical PKC": ["PRKCI", "PRKCZ"],
}


def plot_pkc_screen(screen: pd.DataFrame):
    order = [gene for genes in PKC_GROUPS.values() for gene in genes]
    data = (
        screen[screen.gene.isin(order)]
        .sort_values("mageck_fdr")
        .drop_duplicates("gene")
        .set_index("gene")
        .reindex(order)
        .reset_index()
    )
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    y = np.arange(len(order))[::-1]
    significant = data.mageck_fdr < 0.05
    ax.hlines(y, 0, data.effect, color="0.7", linewidth=0.9)
    ax.scatter(data.effect[~significant], y[~significant], s=25, color="0.55", zorder=3)
    ax.scatter(data.effect[significant], y[significant], s=34, color="#B2182B", zorder=3)
    ax.axvline(0, color="0.25", linestyle="--", linewidth=0.8)
    ax.set_yticks(y, order, fontsize=7)
    ax.set_xlabel("log₂ sgRNA abundance in CD144⁺ / CD144⁻ cells")
    start = 0
    for group, genes in PKC_GROUPS.items():
        center = y[start : start + len(genes)].mean()
        ax.text(1.02, center, group, transform=ax.get_yaxis_transform(), va="center", fontsize=7)
        if start + len(genes) < len(order):
            ax.axhline(y[start + len(genes) - 1] - 0.5, color="0.8", linewidth=0.7)
        start += len(genes)
    ax.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="", color="#B2182B", label="FDR < 0.05"),
            Line2D([], [], marker="o", linestyle="", color="0.55", label="Not significant"),
        ],
        frameon=False,
        fontsize=6,
        loc="upper left",
    )
    return fig, data
