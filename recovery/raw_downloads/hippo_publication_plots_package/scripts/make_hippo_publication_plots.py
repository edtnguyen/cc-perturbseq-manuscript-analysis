from pathlib import Path
import math
import re
import zipfile

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats

# Illustrator-editable text and vector output
mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

IN_DIR = Path("/mnt/data")
OUT_DIR = Path("/mnt/data/hippo_publication_editable_plots")
OUT_DIR.mkdir(exist_ok=True)
DATA_DIR = OUT_DIR / "figure_source_tables"
DATA_DIR.mkdir(exist_ok=True)

paired_path = IN_DIR / "projected_cNMF_program_score_paired_replicate_differences.tsv"
summary_path = IN_DIR / "projected_cNMF_program_score_matched_contrasts.tsv"
key_summary_path = IN_DIR / "key_programs_P15_P17_P34_matched_contrasts.tsv"
enrich_path = IN_DIR / "program_enrichment_from_PyDESeq2_log2FC.tsv"
key_enrich_path = IN_DIR / "key_programs_P15_P17_P34_enrichment_from_PyDESeq2.tsv"
de_summary_path = IN_DIR / "PyDESeq2_contrast_summary.tsv"

paired = pd.read_csv(paired_path, sep="\t")
summary = pd.read_csv(summary_path, sep="\t")
key_summary = pd.read_csv(key_summary_path, sep="\t")
enrich = pd.read_csv(enrich_path, sep="\t")
key_enrich = pd.read_csv(key_enrich_path, sep="\t")
de_summary = pd.read_csv(de_summary_path, sep="\t")

KEY_PROGRAMS = ["P15", "P17", "P34"]

CONTRAST_ORDER = [
    "TDI_at_D1_MPS",
    "IAG_at_D1_MPS",
    "TDI_at_D1.33_LM",
    "IAG_at_D1.33_LM",
    "TDI_at_D1.67_LM",
    "IAG_at_D1.67_LM",
    "TDI_at_D2_LM",
    "IAG_at_D2_LM",
    "TDI_at_D3_AE",
    "IAG_at_D3_AE",
    "TDI_D2_to_D3_vs_D3_AE",
    "IAG_D2_to_D3_vs_D3_AE",
]
CONTRAST_LABELS = {
    "TDI_at_D1_MPS": "D1\nTDI",
    "IAG_at_D1_MPS": "D1\nIAG",
    "TDI_at_D1.33_LM": "D1.33\nTDI",
    "IAG_at_D1.33_LM": "D1.33\nIAG",
    "TDI_at_D1.67_LM": "D1.67\nTDI",
    "IAG_at_D1.67_LM": "D1.67\nIAG",
    "TDI_at_D2_LM": "D2\nTDI",
    "IAG_at_D2_LM": "D2\nIAG",
    "TDI_at_D3_AE": "D3\nTDI",
    "IAG_at_D3_AE": "D3\nIAG",
    "TDI_D2_to_D3_vs_D3_AE": "D2→D3\nTDI",
    "IAG_D2_to_D3_vs_D3_AE": "D2→D3\nIAG",
}
CHEM = {c: ("TDI" if c.startswith("TDI") else "IAG") for c in CONTRAST_ORDER}
CHEM_COLOR = {"TDI": "#D55E00", "IAG": "#0072B2", "control": "#666666"}
REPLICATE_MARKERS = {"rep1": "o", "rep2": "s", "rep3": "^"}

# Filter and order helper columns
paired["contrast_name"] = pd.Categorical(paired["contrast_name"], categories=CONTRAST_ORDER, ordered=True)
paired["program"] = pd.Categorical(paired["program"], categories=[f"P{i}" for i in range(50)], ordered=True)
paired["chem"] = paired["contrast_name"].astype(str).map(CHEM)
paired["contrast_label"] = paired["contrast_name"].astype(str).map(CONTRAST_LABELS)
summary["contrast_name"] = pd.Categorical(summary["contrast_name"], categories=CONTRAST_ORDER, ordered=True)
summary["chem"] = summary["contrast_name"].astype(str).map(CHEM)
summary["contrast_label"] = summary["contrast_name"].astype(str).map(CONTRAST_LABELS)
key_summary["contrast_name"] = pd.Categorical(key_summary["contrast_name"], categories=CONTRAST_ORDER, ordered=True)
key_summary["chem"] = key_summary["contrast_name"].astype(str).map(CHEM)
key_summary["contrast_label"] = key_summary["contrast_name"].astype(str).map(CONTRAST_LABELS)
key_enrich["contrast_name"] = pd.Categorical(key_enrich["contrast_name"], categories=CONTRAST_ORDER, ordered=True)
key_enrich["chem"] = key_enrich["contrast_name"].astype(str).map(CHEM)
key_enrich["contrast_label"] = key_enrich["contrast_name"].astype(str).map(CONTRAST_LABELS)
de_summary["contrast_name"] = pd.Categorical(de_summary["contrast_name"], categories=CONTRAST_ORDER, ordered=True)
de_summary["contrast_label"] = de_summary["contrast_name"].astype(str).map(CONTRAST_LABELS)
de_summary["chem"] = de_summary["contrast_name"].astype(str).map(CHEM)


def save_fig(fig, name):
    for ext in ["svg", "pdf", "png"]:
        fig.savefig(OUT_DIR / f"{name}.{ext}")


def clean_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", length=2, width=0.5)


def sem_ci_95(vals):
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    n = len(vals)
    if n == 0:
        return np.nan, np.nan, np.nan, n
    mean = vals.mean()
    if n == 1:
        return mean, np.nan, np.nan, n
    sem = vals.std(ddof=1) / math.sqrt(n)
    ci = stats.t.ppf(0.975, df=n-1) * sem
    return mean, sem, ci, n

# Figure 1: paired treated-control replicate lines for key programs
sub = paired[paired["program"].astype(str).isin(KEY_PROGRAMS)].copy().sort_values(["program", "contrast_name", "replicate_de"])
fig, axes = plt.subplots(len(KEY_PROGRAMS), 1, figsize=(9.6, 6.2), sharex=True)
if len(KEY_PROGRAMS) == 1:
    axes = [axes]

x_base = np.arange(len(CONTRAST_ORDER))
for ax, program in zip(axes, KEY_PROGRAMS):
    dfp = sub[sub["program"].astype(str) == program]
    for i, contrast in enumerate(CONTRAST_ORDER):
        dfc = dfp[dfp["contrast_name"].astype(str) == contrast]
        color = CHEM_COLOR[CHEM[contrast]]
        for _, row in dfc.iterrows():
            marker = REPLICATE_MARKERS.get(str(row["replicate_de"]), "o")
            ax.plot([i - 0.16, i + 0.16], [row["control_value"], row["treated_value"]], color="#888888", lw=0.7, alpha=0.85, zorder=1)
            ax.scatter(i - 0.16, row["control_value"], s=22, marker=marker, facecolor="white", edgecolor="#666666", lw=0.8, zorder=2)
            ax.scatter(i + 0.16, row["treated_value"], s=22, marker=marker, facecolor=color, edgecolor="black", lw=0.4, zorder=3)
    ax.set_ylabel(f"{program}\nusage")
    ax.axhline(0, color="#dddddd", lw=0.5)
    clean_axes(ax)
axes[-1].set_xticks(x_base)
axes[-1].set_xticklabels([CONTRAST_LABELS[c] for c in CONTRAST_ORDER], rotation=0)
axes[0].set_title("Matched replicate projected cNMF usage: control vs treatment")
legend_handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#666666", label="control", markersize=4),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=CHEM_COLOR["TDI"], markeredgecolor="black", label="TDI", markersize=4),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=CHEM_COLOR["IAG"], markeredgecolor="black", label="IAG", markersize=4),
    Line2D([0], [0], marker="o", color="#888888", markerfacecolor="white", markeredgecolor="black", label="matched replicate", markersize=4, lw=0.7),
]
axes[0].legend(handles=legend_handles, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.05))
fig.tight_layout()
save_fig(fig, "fig1_key_programs_paired_control_vs_treatment")
plt.close(fig)
sub.to_csv(DATA_DIR / "fig1_key_programs_paired_control_vs_treatment_source.tsv", sep="\t", index=False)

# Figure 2: paired differences, individual replicate points + mean 95% CI
fig, axes = plt.subplots(len(KEY_PROGRAMS), 1, figsize=(9.6, 6.0), sharex=True)
if len(KEY_PROGRAMS) == 1:
    axes = [axes]
records = []
for ax, program in zip(axes, KEY_PROGRAMS):
    dfp = sub[sub["program"].astype(str) == program]
    for i, contrast in enumerate(CONTRAST_ORDER):
        dfc = dfp[dfp["contrast_name"].astype(str) == contrast]
        color = CHEM_COLOR[CHEM[contrast]]
        diffs = dfc["paired_diff"].to_numpy(dtype=float)
        mean, sem, ci, n = sem_ci_95(diffs)
        records.append({"program": program, "contrast_name": contrast, "contrast_label": CONTRAST_LABELS[contrast], "chem": CHEM[contrast], "mean_paired_diff": mean, "sem": sem, "ci95": ci, "n": n})
        # jitter replicate points deterministically by rep label
        jitters = {"rep1": -0.08, "rep2": 0.0, "rep3": 0.08}
        for _, row in dfc.iterrows():
            marker = REPLICATE_MARKERS.get(str(row["replicate_de"]), "o")
            ax.scatter(i + jitters.get(str(row["replicate_de"]), 0.0), row["paired_diff"], s=26, marker=marker, facecolor=color, edgecolor="black", lw=0.4, zorder=3)
        if np.isfinite(ci):
            ax.errorbar(i, mean, yerr=ci, fmt="none", ecolor="black", elinewidth=0.9, capsize=2, zorder=4)
        ax.scatter(i, mean, s=46, marker="D", facecolor="white", edgecolor="black", lw=0.8, zorder=5)
    ax.axhline(0, color="#444444", lw=0.7)
    ax.set_ylabel(f"{program}\nΔ usage")
    clean_axes(ax)
axes[-1].set_xticks(x_base)
axes[-1].set_xticklabels([CONTRAST_LABELS[c] for c in CONTRAST_ORDER], rotation=0)
axes[0].set_title("Treatment effect on projected cNMF usage (treated − matched control)")
rep_handles = [Line2D([0], [0], marker=m, color="none", markerfacecolor="#bbbbbb", markeredgecolor="black", label=rep, markersize=4) for rep, m in REPLICATE_MARKERS.items()]
chem_handles = [Patch(facecolor=CHEM_COLOR["TDI"], edgecolor="black", label="TDI"), Patch(facecolor=CHEM_COLOR["IAG"], edgecolor="black", label="IAG")]
mean_handle = Line2D([0], [0], marker="D", color="black", markerfacecolor="white", markeredgecolor="black", label="mean ± 95% CI", markersize=4, lw=0.8)
axes[0].legend(handles=chem_handles + rep_handles + [mean_handle], frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.05))
fig.tight_layout()
save_fig(fig, "fig2_key_programs_paired_differences")
plt.close(fig)
pd.DataFrame(records).to_csv(DATA_DIR / "fig2_key_programs_paired_difference_summary_source.tsv", sep="\t", index=False)

# Figure 3: key program effect heatmap / dot plot of mean paired differences
summary_key = summary[summary["program"].astype(str).isin(KEY_PROGRAMS)].copy()
summary_key["contrast_name"] = pd.Categorical(summary_key["contrast_name"].astype(str), categories=CONTRAST_ORDER, ordered=True)
mat = summary_key.pivot(index="program", columns="contrast_name", values="diff_treated_minus_control").reindex(KEY_PROGRAMS)[CONTRAST_ORDER]
fig, ax = plt.subplots(figsize=(9.6, 2.1))
vals = mat.to_numpy(dtype=float)
vmax = np.nanmax(np.abs(vals))
mesh = ax.pcolormesh(np.arange(vals.shape[1]+1), np.arange(vals.shape[0]+1), vals, cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="flat", linewidth=0.4, edgecolors="white")
ax.set_xlim(0, vals.shape[1]); ax.set_ylim(vals.shape[0], 0)
ax.set_xticks(np.arange(vals.shape[1]) + 0.5); ax.set_xticklabels([CONTRAST_LABELS[c] for c in CONTRAST_ORDER], rotation=0)
ax.set_yticks(np.arange(vals.shape[0]) + 0.5); ax.set_yticklabels(KEY_PROGRAMS)
ax.set_title("Mean treatment effect on projected cNMF usage")
cbar = fig.colorbar(mesh, ax=ax, pad=0.02); cbar.set_label("treated − control")
for r in range(vals.shape[0]):
    for c in range(vals.shape[1]):
        ax.text(c+0.5, r+0.5, f"{vals[r,c]:.3f}", ha="center", va="center", fontsize=5.5)
for spine in ax.spines.values(): spine.set_visible(False)
fig.tight_layout()
save_fig(fig, "fig3_key_programs_effect_heatmap")
plt.close(fig)
mat.to_csv(DATA_DIR / "fig3_key_programs_effect_heatmap_source.tsv", sep="\t")

# Figure 4: enrichment dot/heatmap for key programs
ke = key_enrich.copy()
ke["contrast_name"] = pd.Categorical(ke["contrast_name"].astype(str), categories=CONTRAST_ORDER, ordered=True)
ke = ke[ke["program"].isin(KEY_PROGRAMS)].sort_values(["program", "contrast_name"])
mat_enr = ke.pivot(index="program", columns="contrast_name", values="weighted_mean_log2FC_by_spectra").reindex(KEY_PROGRAMS)[CONTRAST_ORDER]
mat_p = ke.pivot(index="program", columns="contrast_name", values="spearman_padj_bh_within_contrast").reindex(KEY_PROGRAMS)[CONTRAST_ORDER]
fig, ax = plt.subplots(figsize=(9.6, 2.2))
evals = mat_enr.to_numpy(dtype=float)
evmax = np.nanmax(np.abs(evals))
mesh = ax.pcolormesh(np.arange(evals.shape[1]+1), np.arange(evals.shape[0]+1), evals, cmap="RdBu_r", vmin=-evmax, vmax=evmax, shading="flat", linewidth=0.4, edgecolors="white")
ax.set_xlim(0, evals.shape[1]); ax.set_ylim(evals.shape[0], 0)
ax.set_xticks(np.arange(evals.shape[1]) + 0.5); ax.set_xticklabels([CONTRAST_LABELS[c] for c in CONTRAST_ORDER], rotation=0)
ax.set_yticks(np.arange(evals.shape[0]) + 0.5); ax.set_yticklabels(KEY_PROGRAMS)
ax.set_title("DE log2FC enrichment in perturb-seq cNMF spectra")
# annotate stars based on adjusted p
for r in range(evals.shape[0]):
    for c in range(evals.shape[1]):
        p = mat_p.to_numpy(dtype=float)[r, c]
        star = ""
        if np.isfinite(p):
            if p < 0.001: star = "***"
            elif p < 0.01: star = "**"
            elif p < 0.05: star = "*"
        ax.text(c+0.5, r+0.5, f"{evals[r,c]:.2f}{star}", ha="center", va="center", fontsize=5.5)
cbar = fig.colorbar(mesh, ax=ax, pad=0.02); cbar.set_label("spectra-weighted mean DE log2FC")
for spine in ax.spines.values(): spine.set_visible(False)
fig.tight_layout()
save_fig(fig, "fig4_key_programs_de_enrichment_heatmap")
plt.close(fig)
mat_enr.to_csv(DATA_DIR / "fig4_key_programs_de_enrichment_weighted_log2fc_source.tsv", sep="\t")
mat_p.to_csv(DATA_DIR / "fig4_key_programs_de_enrichment_padj_source.tsv", sep="\t")

# Figure 5: gene DE summary up/down counts
fig, ax = plt.subplots(figsize=(8.6, 2.9))
ds = de_summary.sort_values("contrast_name").copy()
x = np.arange(len(CONTRAST_ORDER))
ds = ds.set_index("contrast_name").reindex(CONTRAST_ORDER).reset_index()
up = ds["n_sig_up"].to_numpy(dtype=float)
down = -ds["n_sig_down"].to_numpy(dtype=float)
for i, contrast in enumerate(CONTRAST_ORDER):
    color = CHEM_COLOR[CHEM[contrast]]
    ax.bar(i, up[i], color=color, edgecolor="black", lw=0.3, width=0.72)
    ax.bar(i, down[i], color=color, edgecolor="black", lw=0.3, width=0.72, alpha=0.45)
ax.axhline(0, color="black", lw=0.7)
ax.set_xticks(x); ax.set_xticklabels([CONTRAST_LABELS[c] for c in CONTRAST_ORDER], rotation=0)
ax.set_ylabel("DE genes (padj < 0.05)\nup positive / down negative")
ax.set_title("Pseudobulk PyDESeq2 differential expression summary")
clean_axes(ax)
legend = [Patch(facecolor=CHEM_COLOR["TDI"], edgecolor="black", label="TDI up"), Patch(facecolor=CHEM_COLOR["TDI"], alpha=0.45, edgecolor="black", label="TDI down"), Patch(facecolor=CHEM_COLOR["IAG"], edgecolor="black", label="IAG up"), Patch(facecolor=CHEM_COLOR["IAG"], alpha=0.45, edgecolor="black", label="IAG down")]
ax.legend(handles=legend, frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(1.0, 1.0))
fig.tight_layout()
save_fig(fig, "fig5_pydeseq2_de_gene_counts_summary")
plt.close(fig)
ds.to_csv(DATA_DIR / "fig5_pydeseq2_de_gene_counts_summary_source.tsv", sep="\t", index=False)

# Figure 6: volcano panels for key contrasts (vector but simplified by filtering/alpha)
# Create 6 key volcano plots with vector points. To keep files reasonable, draw all nonsig as light grey, sig as colored.
VOLCANO_CONTRASTS = ["TDI_at_D2_LM", "IAG_at_D2_LM", "TDI_at_D3_AE", "IAG_at_D3_AE", "TDI_D2_to_D3_vs_D3_AE", "IAG_D2_to_D3_vs_D3_AE"]
volcano_files = {Path(p).stem.replace("PyDESeq2_", ""): Path(p) for p in map(str, IN_DIR.glob("PyDESeq2_*.tsv")) if "contrast_summary" not in p}
fig, axes = plt.subplots(2, 3, figsize=(8.2, 5.2), sharex=False, sharey=False)
axes = axes.ravel()
volcano_sources=[]
for ax, contrast in zip(axes, VOLCANO_CONTRASTS):
    pth = volcano_files.get(contrast)
    if pth is None:
        ax.set_axis_off(); continue
    df = pd.read_csv(pth, sep="\t")
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FoldChange", "pvalue"])
    df["minus_log10_p"] = -np.log10(df["pvalue"].clip(lower=1e-300))
    df["sig"] = (df["padj"].fillna(1) < 0.05) & (df["log2FoldChange"].abs() >= 1)
    sig_up = df["sig"] & (df["log2FoldChange"] > 0)
    sig_down = df["sig"] & (df["log2FoldChange"] < 0)
    ax.scatter(df.loc[~df["sig"], "log2FoldChange"], df.loc[~df["sig"], "minus_log10_p"], s=2, color="#bdbdbd", alpha=0.45, linewidths=0)
    ax.scatter(df.loc[sig_up, "log2FoldChange"], df.loc[sig_up, "minus_log10_p"], s=3, color=CHEM_COLOR[CHEM[contrast]], alpha=0.75, linewidths=0)
    ax.scatter(df.loc[sig_down, "log2FoldChange"], df.loc[sig_down, "minus_log10_p"], s=3, color="#444444", alpha=0.75, linewidths=0)
    ax.axvline(0, color="#555555", lw=0.5)
    ax.axvline(1, color="#cccccc", lw=0.5, ls="--")
    ax.axvline(-1, color="#cccccc", lw=0.5, ls="--")
    ax.set_title(CONTRAST_LABELS[contrast].replace("\n", " "), fontsize=7)
    ax.set_xlabel("log2FC")
    ax.set_ylabel("−log10(p)")
    clean_axes(ax)
    volcano_sources.append(df.assign(contrast_name=contrast)[["contrast_name","gene_id","gene_symbol","baseMean","log2FoldChange","pvalue","padj","sig"]])
fig.suptitle("PyDESeq2 volcano plots for key contrasts", y=1.02)
fig.tight_layout()
save_fig(fig, "fig6_key_contrast_volcano_plots")
plt.close(fig)
if volcano_sources:
    pd.concat(volcano_sources, ignore_index=True).to_csv(DATA_DIR / "fig6_key_contrast_volcano_plots_source.tsv.gz", sep="\t", index=False)

# Write a small README
readme = OUT_DIR / "README.txt"
readme.write_text("""Publication-quality Illustrator-editable plots generated from uploaded Hippo time-course result tables.\n\nSVG and PDF files use vector elements and svg.fonttype='none' so text is editable in Illustrator.\n\nFigures:\nfig1_key_programs_paired_control_vs_treatment: Paired replicate lines for P15/P17/P34.\nfig2_key_programs_paired_differences: Treated-control replicate differences and mean ± 95% CI.\nfig3_key_programs_effect_heatmap: Heatmap of mean treated-control effects for key programs.\nfig4_key_programs_de_enrichment_heatmap: Key program DE enrichment heatmap; stars mark spearman padj.\nfig5_pydeseq2_de_gene_counts_summary: Number of up/down DE genes per contrast.\nfig6_key_contrast_volcano_plots: Volcano plots for key D2/D3 contrasts.\n\nSource tables for each figure are in figure_source_tables/.\n""")

# Zip outputs
zip_path = IN_DIR / "hippo_publication_editable_plots.zip"
if zip_path.exists():
    zip_path.unlink()
with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
    for p in OUT_DIR.rglob("*"):
        if p.is_file():
            z.write(p, p.relative_to(IN_DIR))
print("Wrote", zip_path)
print("Output files:")
for p in sorted(OUT_DIR.glob("*.svg")):
    print(p.name)
