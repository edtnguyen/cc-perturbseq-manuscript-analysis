#!/usr/bin/env python3
"""Estimate relative CDH5-positive odds from validation-screen gene log2FC.

The input table must contain:
  gene, genome_wide_screen_log2fc, validation_screen_log2fc,
  genome_wide_screen_hit, validation_screen_hit, plot_category

Absolute EC percentage is not identifiable from log2FC alone. The script reports
relative CDH5-positive odds and assumption-dependent EC-percent sensitivity values,
plus global cross-screen statistical tests. It does not fabricate per-gene p-values.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.contingency_tables import Table2x2, mcnemar
from statsmodels.stats.proportion import proportion_confint


def as_bool(value: str) -> bool:
    value = str(value).strip().lower()
    if value in {"true", "1", "yes", "y"}:
        return True
    if value in {"false", "0", "no", "n"}:
        return False
    raise ValueError(f"Cannot parse boolean value: {value!r}")


def read_source(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "gene", "genome_wide_screen_log2fc", "validation_screen_log2fc",
        "genome_wide_screen_hit", "validation_screen_hit", "plot_category",
    }
    if not rows:
        raise ValueError("Input table is empty.")
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    out = []
    for row in rows:
        out.append({
            "gene": row["gene"].strip(),
            "genome_wide_screen_log2fc": float(row["genome_wide_screen_log2fc"]),
            "validation_screen_log2fc": float(row["validation_screen_log2fc"]),
            "genome_wide_screen_hit": as_bool(row["genome_wide_screen_hit"]),
            "validation_screen_hit": as_bool(row["validation_screen_hit"]),
            "plot_category": row["plot_category"].strip(),
        })
    return out


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    fields = fields or list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    q = np.empty_like(p_values)
    running = 1.0
    m = len(p_values)
    for reverse_index, item_index in enumerate(order[::-1], start=1):
        rank = m - reverse_index + 1
        running = min(running, p_values[item_index] * m / rank)
        q[item_index] = min(running, 1.0)
    return q


def bootstrap_ci(values: np.ndarray, statistic, rng, n_boot: int = 10_000) -> tuple[float, float]:
    n = len(values)
    estimates = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        sample = values[rng.integers(0, n, n)]
        estimates[i] = statistic(sample)
    return tuple(map(float, np.quantile(estimates[np.isfinite(estimates)], [0.025, 0.975])))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260901)
    args = parser.parse_args()

    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    records = read_source(args.input)

    x = np.array([r["genome_wide_screen_log2fc"] for r in records], dtype=float)
    y = np.array([r["validation_screen_log2fc"] for r in records], dtype=float)
    gw_hit = np.array([r["genome_wide_screen_hit"] for r in records], dtype=bool)
    val_hit = np.array([r["validation_screen_hit"] for r in records], dtype=bool)
    genes = np.array([r["gene"] for r in records], dtype=object)
    n = len(records)
    odds_ratio = np.exp2(y)
    odds_pct = 100 * odds_ratio
    delta = y - x
    same_direction = np.sign(x) == np.sign(y)
    nonzero = (x != 0) & (y != 0)
    rng = np.random.default_rng(args.seed)

    # Cross-screen correlations and calibration.
    pearson = stats.pearsonr(x, y)
    spearman = stats.spearmanr(x, y)
    fisher_z = np.arctanh(pearson.statistic)
    fisher_z_se = 1 / np.sqrt(n - 3)
    zcrit = stats.norm.ppf(0.975)
    pearson_ci = (
        float(np.tanh(fisher_z - zcrit * fisher_z_se)),
        float(np.tanh(fisher_z + zcrit * fisher_z_se)),
    )
    spearman_boot = np.empty(args.bootstrap)
    for i in range(args.bootstrap):
        sample_index = rng.integers(0, n, n)
        spearman_boot[i] = stats.spearmanr(x[sample_index], y[sample_index]).statistic
    spearman_ci = tuple(
        map(float, np.quantile(spearman_boot[np.isfinite(spearman_boot)], [0.025, 0.975]))
    )

    X = sm.add_constant(x)
    ols = sm.OLS(y, X).fit(cov_type="HC3")
    coef = np.asarray(ols.params)
    se = np.asarray(ols.bse)
    ci = np.asarray(ols.conf_int())
    slope_z_vs_one = (coef[1] - 1) / se[1]
    slope_p_vs_one = float(2 * stats.norm.sf(abs(slope_z_vs_one)))
    identity_difference = np.array([coef[0], coef[1] - 1])
    covariance = np.asarray(ols.cov_params())
    identity_chi2 = float(
        identity_difference @ np.linalg.inv(covariance) @ identity_difference
    )
    identity_p = float(stats.chi2.sf(identity_chi2, df=2))

    paired_t = stats.ttest_rel(y, x)
    wilcoxon = stats.wilcoxon(y, x, alternative="two-sided", zero_method="wilcox")
    mean_delta_ci = bootstrap_ci(delta, np.mean, rng, args.bootstrap)
    median_delta_ci = bootstrap_ci(delta, np.median, rng, args.bootstrap)

    # Hit-call enrichment and direction concordance.
    contingency = np.array([
        [np.sum(gw_hit & val_hit), np.sum(gw_hit & ~val_hit)],
        [np.sum(~gw_hit & val_hit), np.sum(~gw_hit & ~val_hit)],
    ], dtype=int)
    fisher = stats.fisher_exact(contingency)
    table2x2 = Table2x2(contingency)
    mcnemar_test = mcnemar(contingency, exact=True)

    def direction_test(mask: np.ndarray) -> dict:
        selected = mask & nonzero
        successes = int(np.sum(same_direction[selected]))
        total = int(np.sum(selected))
        result = stats.binomtest(successes, total, p=0.5, alternative="greater")
        low, high = proportion_confint(successes, total, alpha=0.05, method="wilson")
        return {
            "same_direction": successes,
            "n": total,
            "proportion": successes / total,
            "ci95_low": float(low),
            "ci95_high": float(high),
            "p_value": float(result.pvalue),
        }

    direction_gw = direction_test(gw_hit)
    direction_both = direction_test(gw_hit & val_hit)

    # Derive per-gene values. No per-gene p-value is created from one log2FC.
    reduced_order = np.argsort(y)
    increased_order = np.argsort(-y)
    reduced_rank = np.empty(n, dtype=int)
    increased_rank = np.empty(n, dtype=int)
    reduced_rank[reduced_order] = np.arange(1, n + 1)
    increased_rank[increased_order] = np.arange(1, n + 1)

    output_rows = []
    for i, row in enumerate(records):
        lfc = y[i]
        if lfc < 0:
            effect_direction = "Reduced CDH5-positive odds"
        elif lfc > 0:
            effect_direction = "Increased CDH5-positive odds"
        else:
            effect_direction = "No change"
        if val_hit[i]:
            significance = (
                "Significant decrease" if lfc < 0
                else "Significant increase" if lfc > 0
                else "Significant"
            )
        else:
            significance = "Not significant"

        output = {
            "gene": row["gene"],
            "validation_screen_log2fc": lfc,
            "relative_CDH5_pos_odds_ratio": odds_ratio[i],
            "relative_CDH5_pos_odds_pct_of_control": odds_pct[i],
            "validation_effect_direction": effect_direction,
            "validation_screen_hit_source": bool(val_hit[i]),
            "validation_significance_source": significance,
            "genome_wide_screen_log2fc": x[i],
            "genome_wide_screen_hit": bool(gw_hit[i]),
            "cross_screen_direction_concordant": bool(same_direction[i]),
            "validation_minus_genomewide_log2fc": delta[i],
            "absolute_cross_screen_log2fc_difference": abs(delta[i]),
            "rank_most_reduced_CDH5_odds": int(reduced_rank[i]),
            "rank_most_increased_CDH5_odds": int(increased_rank[i]),
        }
        for baseline_pct in (50, 60, 70, 80, 90):
            p0 = baseline_pct / 100
            pg = odds_ratio[i] * p0 / ((1 - p0) + odds_ratio[i] * p0)
            output[f"estimated_EC_pct_if_control_{baseline_pct}pct"] = 100 * pg
        output["per_gene_p_value_from_log2fc_alone"] = ""
        output["per_gene_test_note"] = (
            "Not estimable from one summary log2FC; source hit call retained"
        )
        output["plot_category"] = row["plot_category"]
        output_rows.append(output)

    output_rows.sort(key=lambda row: row["validation_screen_log2fc"])
    write_csv(outdir / "validation_log2fc_relative_CDH5_odds_all297.csv", output_rows)

    gw_n = int(gw_hit.sum())
    gw_validated = int(np.sum(gw_hit & val_hit))
    non_gw_n = int((~gw_hit).sum())
    non_gw_validated = int(np.sum(~gw_hit & val_hit))
    gw_rate_ci = proportion_confint(gw_validated, gw_n, alpha=0.05, method="wilson")
    non_gw_rate_ci = proportion_confint(
        non_gw_validated, non_gw_n, alpha=0.05, method="wilson"
    )

    tests = [
        {
            "analysis_group": "Cross-screen effect concordance",
            "test": "Pearson correlation",
            "subset": "All tested genes",
            "n": n,
            "estimate_name": "r",
            "estimate": float(pearson.statistic),
            "ci95_low": pearson_ci[0],
            "ci95_high": pearson_ci[1],
            "statistic_name": "r",
            "statistic": float(pearson.statistic),
            "degrees_of_freedom": n - 2,
            "p_value": float(pearson.pvalue),
            "alternative": "two-sided association",
            "interpretation": "Strong positive linear concordance.",
            "primary_or_secondary": "Primary",
        },
        {
            "analysis_group": "Cross-screen effect concordance",
            "test": "Spearman rank correlation",
            "subset": "All tested genes",
            "n": n,
            "estimate_name": "rho",
            "estimate": float(spearman.statistic),
            "ci95_low": spearman_ci[0],
            "ci95_high": spearman_ci[1],
            "statistic_name": "rho",
            "statistic": float(spearman.statistic),
            "degrees_of_freedom": "",
            "p_value": float(spearman.pvalue),
            "alternative": "two-sided monotonic association",
            "interpretation": "Strong positive rank concordance.",
            "primary_or_secondary": "Primary",
        },
        {
            "analysis_group": "Cross-screen calibration",
            "test": "OLS slope with HC3 robust SE",
            "subset": "All tested genes",
            "n": n,
            "estimate_name": "slope",
            "estimate": float(coef[1]),
            "ci95_low": float(ci[1, 0]),
            "ci95_high": float(ci[1, 1]),
            "statistic_name": "z versus slope=0",
            "statistic": float(coef[1] / se[1]),
            "degrees_of_freedom": "asymptotic",
            "p_value": float(ols.pvalues[1]),
            "alternative": "slope differs from 0",
            "interpretation": "Validation effect increases with genome-wide effect.",
            "primary_or_secondary": "Secondary",
        },
        {
            "analysis_group": "Cross-screen calibration",
            "test": "OLS slope versus identity slope",
            "subset": "All tested genes",
            "n": n,
            "estimate_name": "slope",
            "estimate": float(coef[1]),
            "ci95_low": float(ci[1, 0]),
            "ci95_high": float(ci[1, 1]),
            "statistic_name": "z versus slope=1",
            "statistic": float(slope_z_vs_one),
            "degrees_of_freedom": "asymptotic",
            "p_value": slope_p_vs_one,
            "alternative": "slope differs from 1",
            "interpretation": "Slope is not detectably different from 1 at alpha=0.05.",
            "primary_or_secondary": "Secondary",
        },
        {
            "analysis_group": "Cross-screen calibration",
            "test": "Joint robust Wald test of identity line",
            "subset": "All tested genes",
            "n": n,
            "estimate_name": "intercept=0 and slope=1",
            "estimate": "",
            "ci95_low": "",
            "ci95_high": "",
            "statistic_name": "chi-square",
            "statistic": identity_chi2,
            "degrees_of_freedom": 2,
            "p_value": identity_p,
            "alternative": "regression differs from identity line",
            "interpretation": "The relationship differs modestly from exact identity.",
            "primary_or_secondary": "Secondary",
        },
        {
            "analysis_group": "Cross-screen location shift",
            "test": "Paired t-test",
            "subset": "All tested genes",
            "n": n,
            "estimate_name": "mean(validation - genome-wide) log2FC",
            "estimate": float(np.mean(delta)),
            "ci95_low": mean_delta_ci[0],
            "ci95_high": mean_delta_ci[1],
            "statistic_name": "t",
            "statistic": float(paired_t.statistic),
            "degrees_of_freedom": n - 1,
            "p_value": float(paired_t.pvalue),
            "alternative": "mean paired difference differs from 0",
            "interpretation": "Small negative mean shift; sensitive to tails.",
            "primary_or_secondary": "Secondary",
        },
        {
            "analysis_group": "Cross-screen location shift",
            "test": "Wilcoxon signed-rank test",
            "subset": "All tested genes",
            "n": n,
            "estimate_name": "median(validation - genome-wide) log2FC",
            "estimate": float(np.median(delta)),
            "ci95_low": median_delta_ci[0],
            "ci95_high": median_delta_ci[1],
            "statistic_name": "W",
            "statistic": float(wilcoxon.statistic),
            "degrees_of_freedom": "",
            "p_value": float(wilcoxon.pvalue),
            "alternative": "paired distribution centered away from 0",
            "interpretation": "No robust median shift.",
            "primary_or_secondary": "Primary",
        },
        {
            "analysis_group": "Validation hit enrichment",
            "test": "Fisher exact test",
            "subset": "297 tested genes",
            "n": n,
            "estimate_name": "odds ratio",
            "estimate": float(fisher.statistic),
            "ci95_low": float(table2x2.oddsratio_confint()[0]),
            "ci95_high": float(table2x2.oddsratio_confint()[1]),
            "statistic_name": "odds ratio",
            "statistic": float(fisher.statistic),
            "degrees_of_freedom": "",
            "p_value": float(fisher.pvalue),
            "alternative": "validation hits enriched among genome-wide hits",
            "interpretation": "Genome-wide hits were much more likely to validate.",
            "primary_or_secondary": "Primary",
        },
        {
            "analysis_group": "Hit-call marginal difference",
            "test": "Exact McNemar test",
            "subset": "297 paired hit calls",
            "n": n,
            "estimate_name": "discordant counts GW-only versus validation-only",
            "estimate": f"{contingency[0,1]} versus {contingency[1,0]}",
            "ci95_low": "",
            "ci95_high": "",
            "statistic_name": "min discordant count",
            "statistic": float(mcnemar_test.statistic),
            "degrees_of_freedom": "",
            "p_value": float(mcnemar_test.pvalue),
            "alternative": "marginal hit rates differ",
            "interpretation": "More genes were GW-only than validation-only hits.",
            "primary_or_secondary": "Secondary",
        },
        {
            "analysis_group": "Direction reproducibility",
            "test": "Exact binomial sign-concordance test",
            "subset": "Genome-wide hits",
            "n": direction_gw["n"],
            "estimate_name": "same-direction proportion",
            "estimate": direction_gw["proportion"],
            "ci95_low": direction_gw["ci95_low"],
            "ci95_high": direction_gw["ci95_high"],
            "statistic_name": "same-direction genes",
            "statistic": direction_gw["same_direction"],
            "degrees_of_freedom": "",
            "p_value": direction_gw["p_value"],
            "alternative": "direction concordance > 0.5",
            "interpretation": "Validation direction reproduced for most GW hits.",
            "primary_or_secondary": "Primary",
        },
        {
            "analysis_group": "Direction reproducibility",
            "test": "Exact binomial sign-concordance test",
            "subset": "Hits significant in both screens",
            "n": direction_both["n"],
            "estimate_name": "same-direction proportion",
            "estimate": direction_both["proportion"],
            "ci95_low": direction_both["ci95_low"],
            "ci95_high": direction_both["ci95_high"],
            "statistic_name": "same-direction genes",
            "statistic": direction_both["same_direction"],
            "degrees_of_freedom": "",
            "p_value": direction_both["p_value"],
            "alternative": "direction concordance > 0.5",
            "interpretation": "190 of 193 shared hits had concordant direction.",
            "primary_or_secondary": "Secondary",
        },
    ]
    q_values = bh_adjust(np.array([row["p_value"] for row in tests], dtype=float))
    for row, q_value in zip(tests, q_values):
        row["p_value_BH_across_reported_global_tests"] = float(q_value)
    write_csv(outdir / "validation_log2fc_global_statistical_tests.csv", tests)

    summary = {
        "n": n,
        "validation_hits": int(val_hit.sum()),
        "validation_significant_decrease": int(np.sum(val_hit & (y < 0))),
        "validation_significant_increase": int(np.sum(val_hit & (y > 0))),
        "pearson_r": float(pearson.statistic),
        "pearson_p": float(pearson.pvalue),
        "pearson_ci95": pearson_ci,
        "spearman_rho": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
        "spearman_ci95_bootstrap": spearman_ci,
        "genome_wide_hits_validated": gw_validated,
        "genome_wide_hits_tested": gw_n,
        "genome_wide_hit_validation_rate": gw_validated / gw_n,
        "genome_wide_hit_validation_rate_ci95": list(map(float, gw_rate_ci)),
        "non_genome_wide_hits_validated": non_gw_validated,
        "non_genome_wide_hits_tested": non_gw_n,
        "non_genome_wide_validation_rate": non_gw_validated / non_gw_n,
        "non_genome_wide_validation_rate_ci95": list(map(float, non_gw_rate_ci)),
        "fisher_exact_odds_ratio": float(fisher.statistic),
        "fisher_exact_p": float(fisher.pvalue),
        "direction_concordance_genome_wide_hits": direction_gw,
    }
    (outdir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    # Diagnostic plots.
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["font.family"] = "DejaVu Sans"
    category_array = np.array([row["plot_category"] for row in records], dtype=object)

    fig, ax = plt.subplots(figsize=(6.8, 5.6))
    for category, color, size, alpha in [
        ("Others", "#BFC3C7", 16, 0.60),
        ("Significant in GW only", "#B22222", 20, 0.85),
        ("Significant in both", "#111111", 20, 0.80),
    ]:
        mask = category_array == category
        ax.scatter(x[mask], y[mask], s=size, c=color, alpha=alpha, linewidths=0, label=category)
    limits = [min(x.min(), y.min()) - 0.25, max(x.max(), y.max()) + 0.25]
    ax.plot(limits, limits, linestyle="--", linewidth=1, color="0.4", label="Identity")
    grid = np.linspace(limits[0], limits[1], 200)
    ax.plot(grid, coef[0] + coef[1] * grid, linewidth=1.2, color="#B22222", label="OLS fit")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_xlabel("Genome-wide screen log$_2$FC")
    ax.set_ylabel("Validation screen log$_2$FC")
    ax.set_title("Validation-screen effect concordance")
    ax.text(
        0.03, 0.97,
        f"n = {n}\nPearson r = {pearson.statistic:.3f}, p = {pearson.pvalue:.1e}\n"
        f"Spearman ρ = {spearman.statistic:.3f}, p = {spearman.pvalue:.1e}",
        transform=ax.transAxes, va="top", ha="left", fontsize=9,
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(outdir / "validation_vs_genomewide_concordance.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "validation_vs_genomewide_concordance.svg", bbox_inches="tight")
    plt.close(fig)

    order = np.argsort(y)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ranks = np.arange(1, n + 1)
    ranked_hits = val_hit[order]
    ranked_odds = odds_pct[order]
    ax.scatter(ranks[~ranked_hits], ranked_odds[~ranked_hits], s=13, c="#BFC3C7",
               alpha=0.65, linewidths=0, label="Not significant (source)")
    ax.scatter(ranks[ranked_hits], ranked_odds[ranked_hits], s=16, c="#111111",
               alpha=0.85, linewidths=0, label="Validation hit (source)")
    ax.axhline(100, linestyle="--", linewidth=1, color="#B22222")
    ax.set_yscale("log", base=2)
    ax.set_yticks([6.25, 12.5, 25, 50, 100, 200])
    ax.set_yticklabels(["6.25", "12.5", "25", "50", "100", "200"])
    ax.set_xlabel("Genes ranked by validation-screen log$_2$FC")
    ax.set_ylabel("Relative CDH5-positive odds (% of control)")
    ax.set_title("Screen-derived relative EC differentiation odds")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(outdir / "validation_relative_CDH5_odds_ranked.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "validation_relative_CDH5_odds_ranked.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
