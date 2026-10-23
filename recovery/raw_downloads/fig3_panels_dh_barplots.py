#!/usr/bin/env python
"""Make Fig. 3 panels d and h: D3 EC fraction bar plots for selected regulators.

This script computes target-specific D3 EC fraction:
  target group  = cells carrying guides targeting the gene
  control group = cells carrying negative-control guides and NOT guides for tested gene
Cells carrying both NT and tested target are excluded from the control group.

Outputs: PNG/PDF/SVG for panel d and h, plus source CSVs.
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.stats import fisher_exact

import fig3_config as cfg
from fig3_utils import ensure_outdir, load_guide_metadata, get_masks_for_target, is_d3, pick_sample_col, savefig_all


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def load_or_compute_fraction_table(adata, guide_meta, targets, panel_name):
    # If user already has a fraction table, try to read compatible rows.
    if cfg.PRECOMPUTED_STATE_FRACTION_CSV:
        df = pd.read_csv(cfg.PRECOMPUTED_STATE_FRACTION_CSV)
        return df

    if cfg.DAY_COL not in adata.obs.columns:
        raise KeyError(f"Missing day column {cfg.DAY_COL} in adata.obs")
    if cfg.CELL_TYPE_COL not in adata.obs.columns:
        raise KeyError(f"Missing cell type column {cfg.CELL_TYPE_COL} in adata.obs")

    d3_mask = is_d3(adata.obs[cfg.DAY_COL].values, cfg.D3_DAY_VALUES)
    ec_re = re.compile(cfg.EC_CELLTYPE_REGEX, re.I)
    ec_mask = adata.obs[cfg.CELL_TYPE_COL].astype(str).map(lambda x: bool(ec_re.search(x))).values
    sample_col = pick_sample_col(adata, cfg.SAMPLE_COL_CANDIDATES)

    rows = []
    all_targets = ["NT"] + list(targets)
    pvals = []
    pval_targets = []

    # compute global NT for first target-specific controls? For NT bar, use NT control from first target if possible,
    # or all NT carriers. We compute NT cells as target-specific control for each target in pairwise tests, and one
    # global NT bar using all NT cells from the first target call.
    global_nt_mask = None

    for target in targets:
        target_mask, control_mask, nt_mask = get_masks_for_target(adata, target, guide_meta, cfg)
        if global_nt_mask is None:
            global_nt_mask = nt_mask

        group_defs = {
            "target": target_mask,
            "control_for_target": control_mask,
        }
        for group, mask in group_defs.items():
            use = d3_mask & mask
            n_total = int(use.sum())
            n_ec = int((use & ec_mask).sum())
            rows.append({
                "panel": panel_name,
                "target": target,
                "group": group,
                "sample": "all_cells",
                "n_total_D3": n_total,
                "n_EC_D3": n_ec,
                "fraction_EC": n_ec / n_total if n_total > 0 else np.nan,
            })
            if sample_col is not None:
                for sample, idx in adata.obs.groupby(sample_col).indices.items():
                    s_mask = np.zeros(adata.n_obs, dtype=bool)
                    s_mask[np.asarray(idx)] = True
                    s_use = use & s_mask
                    s_total = int(s_use.sum())
                    s_ec = int((s_use & ec_mask).sum())
                    rows.append({
                        "panel": panel_name,
                        "target": target,
                        "group": group,
                        "sample": str(sample),
                        "n_total_D3": s_total,
                        "n_EC_D3": s_ec,
                        "fraction_EC": s_ec / s_total if s_total > 0 else np.nan,
                    })

        # Fisher exact target vs target-specific control.
        a = int((d3_mask & target_mask & ec_mask).sum())
        b = int((d3_mask & target_mask & ~ec_mask).sum())
        c = int((d3_mask & control_mask & ec_mask).sum())
        d = int((d3_mask & control_mask & ~ec_mask).sum())
        _, pval = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        pvals.append(pval)
        pval_targets.append(target)

    # Add global NT bar.
    if global_nt_mask is not None:
        use = d3_mask & global_nt_mask
        rows.append({
            "panel": panel_name,
            "target": "NT",
            "group": "global_NT",
            "sample": "all_cells",
            "n_total_D3": int(use.sum()),
            "n_EC_D3": int((use & ec_mask).sum()),
            "fraction_EC": (use & ec_mask).sum() / use.sum() if use.sum() > 0 else np.nan,
        })
        if sample_col is not None:
            for sample, idx in adata.obs.groupby(sample_col).indices.items():
                s_mask = np.zeros(adata.n_obs, dtype=bool)
                s_mask[np.asarray(idx)] = True
                s_use = use & s_mask
                rows.append({
                    "panel": panel_name,
                    "target": "NT",
                    "group": "global_NT",
                    "sample": str(sample),
                    "n_total_D3": int(s_use.sum()),
                    "n_EC_D3": int((s_use & ec_mask).sum()),
                    "fraction_EC": (s_use & ec_mask).sum() / s_use.sum() if s_use.sum() > 0 else np.nan,
                })

    df = pd.DataFrame(rows)
    stats = pd.DataFrame({"target": pval_targets, "p_value": pvals})
    stats["fdr"] = bh_fdr(stats["p_value"])
    return df, stats


def plot_bar_panel(frac_df, stats_df, targets, panel_name, title, out_prefix):
    order = ["NT"] + list(targets)
    # all-cell fraction for bars
    bar = frac_df[(frac_df["sample"] == "all_cells") & (frac_df["target"].isin(order))].copy()
    # For target rows, use group target; for NT row, use global_NT
    bar = bar[((bar["target"] == "NT") & (bar["group"] == "global_NT")) | ((bar["target"] != "NT") & (bar["group"] == "target"))]
    bar["target"] = pd.Categorical(bar["target"], categories=order, ordered=True)
    bar = bar.sort_values("target")

    dots = frac_df[(frac_df["sample"] != "all_cells") & (frac_df["target"].isin(order))].copy()
    dots = dots[((dots["target"] == "NT") & (dots["group"] == "global_NT")) | ((dots["target"] != "NT") & (dots["group"] == "target"))]
    dots["target"] = pd.Categorical(dots["target"], categories=order, ordered=True)

    fig, ax = plt.subplots(figsize=cfg.BAR_FIGSIZE)
    x = np.arange(len(order))
    colors = ["#B0B0B0"] + ["#D8753E"] * len(targets)
    heights = []
    for t in order:
        val = bar.loc[bar["target"].astype(str) == t, "fraction_EC"]
        heights.append(float(val.iloc[0]) * 100 if len(val) else np.nan)
    ax.bar(x, heights, color=colors, edgecolor="none", width=0.72, zorder=2)

    rng = np.random.default_rng(3)
    for i, t in enumerate(order):
        sub = dots[dots["target"].astype(str) == t].dropna(subset=["fraction_EC"])
        if len(sub):
            jitter = rng.normal(0, 0.055, len(sub))
            ax.scatter(np.full(len(sub), i) + jitter, sub["fraction_EC"] * 100, s=6, c="black", alpha=0.85, linewidths=0, zorder=3)
            sem = sub["fraction_EC"].std(ddof=1) / np.sqrt(len(sub)) * 100 if len(sub) > 1 else 0
            mean = sub["fraction_EC"].mean() * 100
            ax.errorbar(i, mean, yerr=sem, color="black", lw=0.8, capsize=2, zorder=4)

    ax.axhline(heights[0], color="#777777", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("Day 3 EC cells (%)", fontsize=8)
    ax.set_title(title, fontsize=9, pad=4)
    ax.tick_params(axis="y", labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, max(100, np.nanmax(heights) * 1.15 if np.isfinite(heights).any() else 100))

    # add FDR labels, if available
    for i, t in enumerate(order[1:], start=1):
        row = stats_df[stats_df["target"] == t]
        if len(row):
            fdr = row["fdr"].iloc[0]
            if fdr < 0.001:
                txt = "***"
            elif fdr < 0.01:
                txt = "**"
            elif fdr < 0.05:
                txt = "*"
            else:
                txt = ""
            if txt:
                ax.text(i, heights[i] + 4, txt, ha="center", va="bottom", fontsize=8)

    fig.tight_layout()
    savefig_all(fig, out_prefix)
    plt.close(fig)


def main():
    ensure_outdir(cfg.OUTDIR)
    guide_meta = load_guide_metadata(cfg.GUIDE_METADATA_TSV, cfg.GUIDE_ID_COL, cfg.TARGET_COL_CANDIDATES)
    adata = sc.read_h5ad(cfg.ADATA_H5AD)

    for panel_name, targets, title in [
        ("d", cfg.PANEL_D_TARGETS, "BMP/TGFb-pathway CRISPRi"),
        ("h", cfg.PANEL_H_TARGETS, "Selected regulator CRISPRi"),
    ]:
        frac, stats = load_or_compute_fraction_table(adata, guide_meta, targets, panel_name)
        frac_path = os.path.join(cfg.OUTDIR, f"panel_{panel_name}_D3_EC_fraction_source.csv")
        stats_path = os.path.join(cfg.OUTDIR, f"panel_{panel_name}_D3_EC_fraction_stats.csv")
        frac.to_csv(frac_path, index=False)
        stats.to_csv(stats_path, index=False)
        plot_bar_panel(frac, stats, targets, panel_name, title, os.path.join(cfg.OUTDIR, f"panel_{panel_name}_D3_EC_fraction_barplot"))
        print("Saved", frac_path, stats_path)


if __name__ == "__main__":
    main()
