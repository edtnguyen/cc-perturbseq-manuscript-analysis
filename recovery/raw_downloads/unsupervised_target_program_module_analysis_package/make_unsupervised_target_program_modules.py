#!/usr/bin/env python
"""
Unsupervised target-module and post-hoc Hippo coherence/enrichment analyses.

Goal: identify coherent target-program effect modules without manually grouping regulators first.
Manual target sets such as Hippo are used only post hoc for labeling/enrichment/coherence tests.
"""
import os
import zipfile
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'  # keep text editable in Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import pdist
from scipy.stats import zscore, fisher_exact, mannwhitneyu
from scipy.stats import pearsonr

# --------------------------
# Inputs
# --------------------------
BASE = "/mnt/data"
EFFECT_CSV = os.path.join(BASE, "outputs_program_first_advanced_full", "target_by_curated_program_effect_matrix.csv")
FDR_CSV = os.path.join(BASE, "outputs_program_first_advanced_full", "target_by_curated_program_fdr_matrix.csv")
OUTDIR = os.path.join(BASE, "unsupervised_target_program_module_analysis")
os.makedirs(OUTDIR, exist_ok=True)

PROGRAM_ORDER = ["P7","P22","P29","P28","P13","P30","P48","P33","P15","P14","P17","P34","P27","P2","P3","P8","P11"]
PROGRAM_LABELS = {
    "P7": "P7\npluri", "P22": "P22\npluri", "P29": "P29\npluri",
    "P28": "P28\nD1", "P13": "P13\nPS", "P30": "P30\nPS", "P48": "P48\nneuro", "P33": "P33\nD1",
    "P15": "P15\nvascular", "P14": "P14\nLM", "P17": "P17\nYAP-LM",
    "P34": "P34\nYAP-mesen", "P27": "P27\nECM-mesen",
    "P2": "P2\nEC", "P3": "P3\nEC", "P8": "P8\nEC", "P11": "P11\nEC",
}
PROGRAM_GROUPS = {
    "D0/pluri": ["P7","P22","P29"],
    "D1/early": ["P28","P13","P30","P48","P33"],
    "D2/transition": ["P15","P14","P17"],
    "D3/alternative": ["P34","P27"],
    "D3/EC": ["P2","P3","P8","P11"],
}
GROUP_COLORS = {
    "D0/pluri": "#BDBDBD",
    "D1/early": "#E69F00",
    "D2/transition": "#56B4E9",
    "D3/alternative": "#D55E00",
    "D3/EC": "#0072B2",
}

# Target sets are used only post hoc.
POSTHOC_TARGET_SETS = {
    "Hippo/YAP": ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"],
    "BMP/TGF-SMAD": ["SMAD1", "SMAD4", "SMAD5", "LEMD3", "BMPR1A", "ACVR2A", "ACVR1B", "ACVR2B"],
    "EC regulators": ["KDR", "ETV2", "ERG", "FLI1", "TAL1", "SOX17", "CDH5", "PECAM1"],
    "Candidate/other": ["PRKCD"],
}
HIGHLIGHT_TARGETS = sorted(set(sum(POSTHOC_TARGET_SETS.values(), [])))
HIPPO = POSTHOC_TARGET_SETS["Hippo/YAP"]

# --------------------------
# Load and preprocess
# --------------------------
effects = pd.read_csv(EFFECT_CSV, index_col=0)
fdr = pd.read_csv(FDR_CSV, index_col=0)
# normalize index name
if effects.index.name is None:
    effects.index.name = "target"
if fdr.index.name is None:
    fdr.index.name = "target"

# Keep program order present in data.
programs = [p for p in PROGRAM_ORDER if p in effects.columns]
effects = effects[programs].copy()
fdr = fdr[programs].copy()
effects = effects.loc[effects.index.notna()]
effects = effects.replace([np.inf, -np.inf], np.nan).fillna(0)
fdr = fdr.reindex(index=effects.index, columns=programs)

# Row-normalized matrix for pattern-based clustering.
X_raw = effects.values.astype(float)
X_z = zscore(X_raw, axis=1, nan_policy="omit")
X_z = np.nan_to_num(X_z, nan=0.0, posinf=0.0, neginf=0.0)

# Correlation distance clustering on row-normalized profiles.
# This is pattern-focused; the heatmaps below display raw effects.
row_link = linkage(pdist(X_z, metric="correlation"), method="average")
row_order = leaves_list(row_link)
ordered_targets = effects.index[row_order].tolist()

N_MODULES = 8
cluster_ids = fcluster(row_link, t=N_MODULES, criterion="maxclust")
cluster_series = pd.Series(cluster_ids, index=effects.index, name="unsup_module_raw_id")
# Re-label modules by dendrogram order for readability.
module_order_by_tree = []
for t in ordered_targets:
    m = cluster_series.loc[t]
    if m not in module_order_by_tree:
        module_order_by_tree.append(m)
module_map = {old: i+1 for i, old in enumerate(module_order_by_tree)}
modules = cluster_series.map(module_map).rename("module")

# Module summary tables.
module_table = pd.DataFrame({"target": effects.index, "module": modules.values}).set_index("target")
module_table["is_hippo"] = module_table.index.isin(HIPPO)
for set_name, genes in POSTHOC_TARGET_SETS.items():
    module_table[f"is_{set_name.replace('/', '_').replace('-', '_').replace(' ', '_')}"] = module_table.index.isin(genes)
module_table.to_csv(os.path.join(OUTDIR, "unsupervised_target_module_assignments.csv"))

# --------------------------
# Helpers
# --------------------------
def savefig(fig, stem):
    paths = []
    for ext in ["png", "pdf", "svg"]:
        path = os.path.join(OUTDIR, f"{stem}.{ext}")
        fig.savefig(path, dpi=350 if ext == "png" else None, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def program_group_for(p):
    for g, plist in PROGRAM_GROUPS.items():
        if p in plist:
            return g
    return "other"


def mean_pairwise_corr(mat):
    # mat shape = genes x programs
    if mat.shape[0] < 2:
        return np.nan
    cc = np.corrcoef(mat)
    iu = np.triu_indices_from(cc, 1)
    vals = cc[iu]
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return np.nan
    return float(np.mean(vals))

# --------------------------
# 1. Unsupervised all-target heatmap with post-hoc Hippo ticks
# --------------------------
ordered_effects = effects.loc[ordered_targets, programs]
ordered_modules = modules.loc[ordered_targets]

fig_w = 9.6
fig_h = 11.0
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
clip_val = 3.0
imshow_data = ordered_effects.clip(-clip_val, clip_val).values
norm = TwoSlopeNorm(vmin=-clip_val, vcenter=0, vmax=clip_val)
im = ax.imshow(imshow_data, aspect="auto", cmap="RdBu_r", norm=norm, interpolation="nearest")

ax.set_xticks(np.arange(len(programs)))
ax.set_xticklabels([PROGRAM_LABELS.get(p, p) for p in programs], rotation=45, ha="right", fontsize=8)
ax.set_yticks([])
ax.set_title("Unsupervised clustering of target-program effect profiles", fontsize=13, weight="bold", pad=22)
ax.set_xlabel("Expert-curated cNMF programs", fontsize=10)
ax.set_ylabel("Perturbed targets clustered without pathway labels", fontsize=10)

# module side bar
module_colors = plt.cm.tab20(np.linspace(0, 1, N_MODULES))
mod_color_map = {i+1: module_colors[i] for i in range(N_MODULES)}
for y, m in enumerate(ordered_modules.values):
    ax.add_patch(Rectangle((-1.05, y-0.5), 0.28, 1.0, facecolor=mod_color_map[m], edgecolor="none", clip_on=False))

# program group top bar
for x, p in enumerate(programs):
    g = program_group_for(p)
    ax.add_patch(Rectangle((x-0.5, -6.0), 1.0, 2.0, facecolor=GROUP_COLORS.get(g, "#AAAAAA"), edgecolor="none", clip_on=False))

# post-hoc Hippo ticks/labels
hippo_y = []
for y, t in enumerate(ordered_targets):
    if t in HIPPO:
        hippo_y.append(y)
        ax.plot([len(programs)-0.25, len(programs)+0.15], [y, y], color="#000000", lw=1.2, clip_on=False)
        ax.text(len(programs)+0.25, y, t, va="center", ha="left", fontsize=8, weight="bold", clip_on=False)

# annotate module boundaries
last_m = ordered_modules.iloc[0]
start = 0
for y, m in enumerate(ordered_modules.values.tolist() + [None]):
    if y == len(ordered_modules) or m != last_m:
        ax.axhline(y - 0.5, color="white", lw=1.2)
        mid = (start + y - 1) / 2
        ax.text(-1.18, mid, f"M{last_m}", va="center", ha="right", fontsize=7, color="black", clip_on=False)
        start = y
        last_m = m

cbar = fig.colorbar(im, ax=ax, shrink=0.55, pad=0.015)
cbar.set_label("program effect size\n(clipped)", fontsize=8)
cbar.ax.tick_params(labelsize=7)

# legends for program groups / module bar
legend_y = -0.13
ax.text(0.0, legend_y, "Top bar: program group", transform=ax.transAxes, fontsize=8, va="center")
x0 = 0.25
for i, (g, col) in enumerate(GROUP_COLORS.items()):
    ax.add_patch(Rectangle((x0 + i*0.13, legend_y-0.015), 0.018, 0.03, transform=ax.transAxes, facecolor=col, edgecolor="none", clip_on=False))
    ax.text(x0 + i*0.13 + 0.022, legend_y, g, transform=ax.transAxes, fontsize=7, va="center")
ax.text(0.0, legend_y-0.055, "Left bar: unsupervised target module; right ticks: Hippo genes highlighted post hoc", transform=ax.transAxes, fontsize=8, va="center")

savefig(fig, "unsupervised_all_target_program_heatmap_hippo_posthoc")

# --------------------------
# 2. Module centroid heatmap
# --------------------------
centroids = effects.join(modules).groupby("module")[programs].mean()
module_sizes = modules.value_counts().sort_index()
hippo_counts = module_table.groupby("module")["is_hippo"].sum().astype(int)
fig, ax = plt.subplots(figsize=(8.6, 4.2))
cent = centroids.loc[range(1, N_MODULES+1), programs]
clip_cent = cent.clip(-2.5, 2.5)
im = ax.imshow(clip_cent.values, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5))
ax.set_xticks(np.arange(len(programs)))
ax.set_xticklabels([PROGRAM_LABELS.get(p,p) for p in programs], rotation=45, ha="right", fontsize=8)
ax.set_yticks(np.arange(N_MODULES))
ax.set_yticklabels([f"M{m}  n={module_sizes.get(m,0)}  Hippo={hippo_counts.get(m,0)}" for m in range(1,N_MODULES+1)], fontsize=8)
ax.set_title("Unsupervised target modules summarized by mean program effects", fontsize=12, weight="bold", pad=12)
for y in np.arange(N_MODULES+1)-0.5:
    ax.axhline(y, color="white", lw=1)
for x in np.arange(len(programs)+1)-0.5:
    ax.axvline(x, color="white", lw=0.6)
cbar = fig.colorbar(im, ax=ax, shrink=0.72, pad=0.02)
cbar.set_label("mean effect", fontsize=8)
cbar.ax.tick_params(labelsize=7)
fig.tight_layout()
savefig(fig, "unsupervised_module_centroid_program_effects")
centroids.to_csv(os.path.join(OUTDIR, "unsupervised_module_centroid_program_effects.csv"))

# --------------------------
# 3. Module enrichment table and heatmap (post-hoc tests)
# --------------------------
rows = []
all_targets = set(effects.index)
for m in range(1, N_MODULES+1):
    in_module = set(module_table.index[module_table["module"] == m])
    for set_name, genes in POSTHOC_TARGET_SETS.items():
        genes_present = set([g for g in genes if g in all_targets])
        a = len(in_module & genes_present)
        b = len(in_module - genes_present)
        c = len(genes_present - in_module)
        d = len(all_targets - in_module - genes_present)
        odds, p = fisher_exact([[a,b],[c,d]], alternative="greater") if len(genes_present) > 0 else (np.nan, np.nan)
        rows.append({
            "module": m,
            "target_set": set_name,
            "overlap": a,
            "module_size": len(in_module),
            "set_size": len(genes_present),
            "odds_ratio": odds,
            "p_value": p,
            "genes_in_module": ";".join(sorted(in_module & genes_present)),
        })
enrich = pd.DataFrame(rows)
# BH FDR within all module x set tests
pvals = enrich["p_value"].values
order = np.argsort(pvals)
rank = np.empty_like(order)
rank[order] = np.arange(1, len(pvals)+1)
enrich["fdr"] = np.minimum.accumulate((pvals[order][::-1] * len(pvals) / np.arange(len(pvals),0,-1)))[::-1][np.argsort(order)] if len(pvals)>0 else []
# simpler correction to avoid weird order issues
sorted_p = pvals[order]
adj = np.empty_like(sorted_p, dtype=float)
prev = 1.0
for i in range(len(sorted_p)-1, -1, -1):
    val = sorted_p[i] * len(sorted_p) / (i+1)
    prev = min(prev, val)
    adj[i] = prev
enrich["fdr"] = pd.Series(adj, index=enrich.index[order]).sort_index().values
enrich.to_csv(os.path.join(OUTDIR, "posthoc_target_set_enrichment_by_unsupervised_module.csv"), index=False)

pivot = enrich.pivot(index="module", columns="target_set", values="overlap").reindex(range(1,N_MODULES+1))
pivot_fdr = enrich.pivot(index="module", columns="target_set", values="fdr").reindex(range(1,N_MODULES+1))
fig, ax = plt.subplots(figsize=(5.8, 4.2))
vals = pivot.fillna(0).values
im = ax.imshow(vals, aspect="auto", cmap="Blues", interpolation="nearest")
ax.set_xticks(np.arange(pivot.shape[1]))
ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
ax.set_yticks(np.arange(pivot.shape[0]))
ax.set_yticklabels([f"M{m}" for m in pivot.index], fontsize=8)
ax.set_title("Post-hoc target-set overlap by unsupervised module", fontsize=11, weight="bold")
for i in range(vals.shape[0]):
    for j in range(vals.shape[1]):
        overlap = int(vals[i,j])
        f = pivot_fdr.iloc[i,j]
        txt = str(overlap)
        if np.isfinite(f) and f < 0.05 and overlap > 0:
            txt += "*"
        ax.text(j, i, txt, ha="center", va="center", fontsize=8, color="black" if overlap < vals.max()/2 else "white")
fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02, label="overlap count")
fig.tight_layout()
savefig(fig, "posthoc_target_set_overlap_by_unsupervised_module")

# --------------------------
# 4. Zoom heatmap of the Hippo-containing module
# --------------------------
hippo_module_counts = module_table[module_table.index.isin(HIPPO)].groupby("module").size()
if len(hippo_module_counts) > 0:
    hippo_module = int(hippo_module_counts.sort_values(ascending=False).index[0])
else:
    hippo_module = int(hippo_counts.sort_values(ascending=False).index[0])
module_targets = module_table.index[module_table["module"] == hippo_module].tolist()
# sort within module by P34 descending, then EC mean ascending if available
mod_df = effects.loc[module_targets, programs].copy()
if all(p in mod_df.columns for p in ["P2","P3","P8","P11"]):
    mod_df["EC_mean"] = mod_df[["P2","P3","P8","P11"]].mean(axis=1)
else:
    mod_df["EC_mean"] = 0
sort_cols = []
if "P34" in mod_df.columns: sort_cols.append("P34")
sort_cols.append("EC_mean")
ascending = [False, True]
mod_order = mod_df.sort_values(sort_cols, ascending=ascending).index.tolist()
zoom = effects.loc[mod_order, programs]
fig_h = max(4.5, 0.22*len(mod_order) + 1.5)
fig, ax = plt.subplots(figsize=(8.2, fig_h))
im = ax.imshow(zoom.clip(-3,3).values, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3))
ax.set_xticks(np.arange(len(programs)))
ax.set_xticklabels([PROGRAM_LABELS.get(p,p) for p in programs], rotation=45, ha="right", fontsize=8)
# label only module members if module not too huge; otherwise label highlighted genes only
if len(mod_order) <= 45:
    ylabels = [t if (t in HIGHLIGHT_TARGETS or t in HIPPO) else t for t in mod_order]
    ax.set_yticks(np.arange(len(mod_order)))
    ax.set_yticklabels(ylabels, fontsize=7)
    for label in ax.get_yticklabels():
        if label.get_text() in HIPPO:
            label.set_fontweight("bold"); label.set_color("#000000")
        elif label.get_text() in HIGHLIGHT_TARGETS:
            label.set_color("#444444")
else:
    ax.set_yticks([])
    for y, t in enumerate(mod_order):
        if t in HIGHLIGHT_TARGETS:
            ax.text(len(programs)+0.2, y, t, va="center", ha="left", fontsize=8, weight="bold" if t in HIPPO else "normal", clip_on=False)
            ax.plot([len(programs)-0.2, len(programs)+0.1], [y,y], color="black", lw=0.8, clip_on=False)
ax.set_title(f"Zoom: unsupervised module M{hippo_module} containing Hippo genes", fontsize=12, weight="bold")
fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02, label="program effect")
fig.tight_layout()
savefig(fig, "zoom_hippo_containing_unsupervised_module")
pd.DataFrame({"target": mod_order, "module": hippo_module}).to_csv(os.path.join(OUTDIR, "targets_in_hippo_containing_unsupervised_module.csv"), index=False)

# --------------------------
# 5. Coherence permutation test for Hippo target profiles
# --------------------------
hippo_present = [g for g in HIPPO if g in effects.index]
obs_coherence = mean_pairwise_corr(effects.loc[hippo_present, programs].values)
rng = np.random.default_rng(12345)
all_targets_list = effects.index.tolist()
N_PERM = 20000
null_vals = np.empty(N_PERM)
for i in range(N_PERM):
    sample = rng.choice(all_targets_list, size=len(hippo_present), replace=False)
    null_vals[i] = mean_pairwise_corr(effects.loc[sample, programs].values)
emp_p = (np.sum(null_vals >= obs_coherence) + 1) / (N_PERM + 1)
coh_stats = pd.DataFrame([{
    "target_set": "Hippo/YAP",
    "n_targets": len(hippo_present),
    "observed_mean_pairwise_corr": obs_coherence,
    "empirical_p_greater": emp_p,
    "null_mean": np.mean(null_vals),
    "null_sd": np.std(null_vals),
    "n_permutations": N_PERM,
    "targets": ";".join(hippo_present),
}])
coh_stats.to_csv(os.path.join(OUTDIR, "hippo_profile_coherence_permutation_test.csv"), index=False)

fig, ax = plt.subplots(figsize=(5.8, 4.2))
ax.hist(null_vals, bins=55, color="#BDBDBD", edgecolor="white", alpha=0.95)
ax.axvline(obs_coherence, color="#D73027", lw=2.2)
ax.text(obs_coherence, ax.get_ylim()[1]*0.88, f" Hippo observed\n mean r = {obs_coherence:.2f}\n empirical P = {emp_p:.4f}", ha="left", va="top", fontsize=9, color="#D73027")
ax.set_title("Post-hoc coherence test for Hippo program-effect profiles", fontsize=11, weight="bold")
ax.set_xlabel("Mean pairwise correlation of target effect profiles")
ax.set_ylabel("Random target sets")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
savefig(fig, "hippo_effect_profile_coherence_permutation")

# --------------------------
# 6. Ranked/rug enrichment plots for individual program metrics
# --------------------------
metrics = pd.DataFrame(index=effects.index)
metrics["P34 increase"] = effects["P34"] if "P34" in effects.columns else np.nan
metrics["P17 increase"] = effects["P17"] if "P17" in effects.columns else np.nan
metrics["P15 decrease"] = -effects["P15"] if "P15" in effects.columns else np.nan
if all(p in effects.columns for p in ["P2","P3","P8","P11"]):
    metrics["EC endpoint decrease"] = -effects[["P2","P3","P8","P11"]].mean(axis=1)
else:
    metrics["EC endpoint decrease"] = np.nan
if "P27" in effects.columns:
    metrics["P27 increase"] = effects["P27"]

rank_rows = []
fig, axes = plt.subplots(len(metrics.columns), 1, figsize=(7.2, 1.65*len(metrics.columns)), sharex=False)
if len(metrics.columns) == 1:
    axes = [axes]
for ax, metric in zip(axes, metrics.columns):
    ser = metrics[metric].dropna().sort_values(ascending=False)
    ranks = np.arange(1, len(ser)+1)
    ax.plot(ranks, ser.values, color="#333333", lw=1.2)
    ax.axhline(0, color="#999999", lw=0.8, ls="--")
    # Hippo rug ticks
    hippo_ranks = []
    for g in hippo_present:
        if g in ser.index:
            r = int(np.where(ser.index.values == g)[0][0]) + 1
            hippo_ranks.append(r)
            y = ser.loc[g]
            ax.plot([r, r], [ax.get_ylim()[0], y], color="#D73027", lw=1.0, alpha=0.8)
            ax.scatter([r], [y], s=28, color="#D73027", zorder=3, edgecolor="white", linewidth=0.4)
            ax.text(r, y, f" {g}", fontsize=7, va="center", ha="left", color="#D73027")
    hippo_vals = ser.loc[[g for g in hippo_present if g in ser.index]].values
    non_vals = ser.loc[[g for g in ser.index if g not in hippo_present]].values
    u_stat, p_mwu = mannwhitneyu(hippo_vals, non_vals, alternative="greater")
    mean_rank = float(np.mean(hippo_ranks)) if hippo_ranks else np.nan
    rank_rows.append({
        "metric": metric,
        "mean_hippo_metric": np.mean(hippo_vals),
        "mean_nonhippo_metric": np.mean(non_vals),
        "mannwhitney_p_greater": p_mwu,
        "mean_hippo_rank": mean_rank,
        "hippo_ranks": ";".join(map(str, hippo_ranks)),
    })
    ax.set_ylabel(metric, fontsize=9)
    ax.text(0.99, 0.88, f"MWU P={p_mwu:.2g}", transform=ax.transAxes, ha="right", va="top", fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
axes[-1].set_xlabel("Targets ranked by metric (high to low)")
fig.suptitle("Post-hoc Hippo placement on ranked target-effect metrics", fontsize=12, weight="bold", y=1.01)
fig.tight_layout()
savefig(fig, "hippo_rank_rug_enrichment_across_program_metrics")
pd.DataFrame(rank_rows).to_csv(os.path.join(OUTDIR, "hippo_rank_rug_enrichment_statistics.csv"), index=False)
metrics.to_csv(os.path.join(OUTDIR, "target_rank_metrics_for_hippo_enrichment.csv"))

# --------------------------
# 7. Scatter of first two module-relevant raw programs (not MDS)
# --------------------------
# These are interpretable raw axes rather than embedding axes.
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
plot_specs = [
    ("P17", "P34", "P17 effect", "P34 effect", "P17 vs P34 raw program effects"),
    ("P15", "P34", "P15 effect", "P34 effect", "P15 vs P34 raw program effects"),
]
for ax, (xcol, ycol, xlabel, ylabel, title) in zip(axes, plot_specs):
    x = effects[xcol]
    y = effects[ycol]
    ax.scatter(x, y, s=16, c="#BDBDBD", alpha=0.55, edgecolors="none")
    # highlighted Hippo posthoc
    for g in hippo_present:
        ax.scatter(effects.loc[g, xcol], effects.loc[g, ycol], s=58, color="#D73027", edgecolor="white", linewidth=0.5, zorder=3)
        ax.text(effects.loc[g, xcol], effects.loc[g, ycol], f" {g}", fontsize=8, va="center", ha="left", color="#D73027")
    rr, pp = pearsonr(x, y)
    ax.axhline(0, color="#999999", ls="--", lw=0.8); ax.axvline(0, color="#999999", ls="--", lw=0.8)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title, fontsize=10, weight="bold")
    ax.text(0.03, 0.97, f"all targets r={rr:.2f}", transform=ax.transAxes, va="top", fontsize=8)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
savefig(fig, "posthoc_hippo_on_interpretable_raw_program_axes")

# --------------------------
# README and package
# --------------------------
readme = f"""
Unsupervised target-module / post-hoc Hippo analysis
====================================================

Inputs:
- {EFFECT_CSV}
- {FDR_CSV}

Key design:
- Targets were clustered using effect profiles across expert-curated cNMF programs.
- No pathway/regulator labels were used during clustering.
- Hippo/YAP, BMP/TGF-SMAD, EC regulator, and PRKCD labels were used only after clustering for highlighting/enrichment.

Main outputs:
1. unsupervised_all_target_program_heatmap_hippo_posthoc.*
   All targets clustered without pathway labels; Hippo genes highlighted post hoc.
2. unsupervised_module_centroid_program_effects.*
   Mean raw program effects per unsupervised module.
3. posthoc_target_set_overlap_by_unsupervised_module.*
   Post-hoc overlap/enrichment of target sets in unsupervised modules.
4. zoom_hippo_containing_unsupervised_module.*
   Raw program effects for the module containing the most Hippo genes.
5. hippo_effect_profile_coherence_permutation.*
   Random-set test asking whether Hippo effect profiles are more coherent than random target sets of same size.
6. hippo_rank_rug_enrichment_across_program_metrics.*
   Rank/rug plots showing where Hippo genes fall among all targets for P34 up, P17 up, P15 down, EC endpoint down, etc.

Observed Hippo coherence:
- mean pairwise r = {obs_coherence:.3f}
- empirical P = {emp_p:.4g}
"""
with open(os.path.join(OUTDIR, "README_unsupervised_module_analysis.txt"), "w") as f:
    f.write(readme)

script_out = os.path.join(OUTDIR, "make_unsupervised_target_program_modules.py")
try:
    import shutil
    shutil.copyfile(__file__, script_out)
except Exception:
    pass

zip_path = os.path.join(BASE, "unsupervised_target_program_module_analysis_package.zip")
with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
    for root, _, files in os.walk(OUTDIR):
        for fn in files:
            fp = os.path.join(root, fn)
            z.write(fp, arcname=os.path.relpath(fp, OUTDIR))

print("Generated outputs in", OUTDIR)
print("ZIP", zip_path)
print("Hippo module:", hippo_module)
print(coh_stats.to_string(index=False))
print(enrich[(enrich['target_set']=='Hippo/YAP') & (enrich['overlap']>0)].sort_values('p_value').to_string(index=False))

# --------------------------
# 8. Robustness of Hippo co-clustering across dendrogram cuts
# --------------------------
cut_rows = []
for K in range(5, 61):
    cl = pd.Series(fcluster(row_link, t=K, criterion='maxclust'), index=effects.index)
    best = None
    for m in sorted(cl.unique()):
        mod = set(cl.index[cl == m])
        a = len(mod & set(hippo_present))
        if a == 0:
            continue
        b = len(mod - set(hippo_present))
        c = len(set(hippo_present) - mod)
        d = len(set(effects.index) - mod - set(hippo_present))
        odds, p = fisher_exact([[a,b],[c,d]], alternative='greater')
        row = {
            'n_modules_cut': K,
            'module_raw_id': int(m),
            'hippo_overlap': a,
            'module_size': len(mod),
            'fisher_p_greater': p,
            'hippo_genes': ';'.join(sorted(mod & set(hippo_present))),
        }
        if best is None or (a > best['hippo_overlap']) or (a == best['hippo_overlap'] and len(mod) < best['module_size']):
            best = row
    cut_rows.append(best)
cut_df = pd.DataFrame(cut_rows)
cut_df.to_csv(os.path.join(OUTDIR, 'hippo_coclustering_across_dendrogram_cuts.csv'), index=False)

fig, ax1 = plt.subplots(figsize=(6.2, 4.2))
ax1.plot(cut_df['n_modules_cut'], cut_df['module_size'], color='#333333', lw=1.8)
ax1.scatter(cut_df['n_modules_cut'], cut_df['module_size'], s=12, color='#333333')
ax1.set_xlabel('Number of modules from hierarchical cut')
ax1.set_ylabel('Size of module containing most Hippo genes')
ax1.spines['top'].set_visible(False)
ax2 = ax1.twinx()
ax2.plot(cut_df['n_modules_cut'], cut_df['hippo_overlap'], color='#D73027', lw=1.8)
ax2.scatter(cut_df['n_modules_cut'], cut_df['hippo_overlap'], s=18, color='#D73027')
ax2.set_ylabel('Hippo genes in same module', color='#D73027')
ax2.tick_params(axis='y', colors='#D73027')
ax2.set_ylim(0, len(hippo_present)+0.5)
ax1.set_title('Post-hoc robustness of Hippo co-clustering', fontsize=11, weight='bold')
ax1.text(0.03, 0.05, 'Clustering used all targets and curated programs\nHippo labels applied only after clustering', transform=ax1.transAxes, fontsize=8, va='bottom')
fig.tight_layout()
savefig(fig, 'hippo_coclustering_robustness_across_dendrogram_cuts')

# --------------------------
# 9. Pairwise Hippo profile correlation heatmap
# --------------------------
hippo_mat = effects.loc[hippo_present, programs]
hippo_corr = hippo_mat.T.corr()
hippo_corr.to_csv(os.path.join(OUTDIR, 'hippo_pairwise_program_effect_profile_correlations.csv'))
fig, ax = plt.subplots(figsize=(4.2, 3.8))
im = ax.imshow(hippo_corr.values, cmap='Blues', vmin=0, vmax=1)
ax.set_xticks(np.arange(len(hippo_present)))
ax.set_xticklabels(hippo_present, rotation=45, ha='right', fontsize=8)
ax.set_yticks(np.arange(len(hippo_present)))
ax.set_yticklabels(hippo_present, fontsize=8)
for i in range(len(hippo_present)):
    for j in range(len(hippo_present)):
        ax.text(j, i, f'{hippo_corr.iloc[i,j]:.2f}', ha='center', va='center', fontsize=8, color='white' if hippo_corr.iloc[i,j] > 0.65 else 'black')
ax.set_title('Pairwise similarity of Hippo program-effect profiles', fontsize=10.5, weight='bold')
fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02, label='Pearson r')
fig.tight_layout()
savefig(fig, 'hippo_pairwise_program_effect_profile_correlations')

# Repack after appended analyses
zip_path = os.path.join(BASE, 'unsupervised_target_program_module_analysis_package.zip')
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
    for root, _, files in os.walk(OUTDIR):
        for fn in files:
            fp = os.path.join(root, fn)
            z.write(fp, arcname=os.path.relpath(fp, OUTDIR))
print('Updated ZIP with robustness/coherence figures:', zip_path)
