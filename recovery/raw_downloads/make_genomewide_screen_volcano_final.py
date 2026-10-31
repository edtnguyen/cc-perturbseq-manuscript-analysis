#!/usr/bin/env python3
"""Final Nature-style, Illustrator-editable volcano plot.

- Uses directional MAGeCK FDR: neg|fdr for negative log2FC and pos|fdr for positive log2FC.
- Uses editable TrueType text in PDF and SVG.
- Rasterizes only the dense nonsignificant background point layer for Illustrator performance.
"""
import os, re, zipfile
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "legend.title_fontsize": 6,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.55,
    "ytick.major.width": 0.55,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

INFILE = "/mnt/data/Pos_vs_Neg.gene_summary.txt"
OUTDIR = "/mnt/data/genomewide_screen_volcano_nature_final"
os.makedirs(OUTDIR, exist_ok=True)

FDR_CUTOFF = 0.05
COL_NON = "#C9C9C9"
COL_NEG = "#2166AC"
COL_POS = "#D6604D"

# Load data.
df = pd.read_csv(INFILE, sep="\t")
df["gene"] = df["id"].astype(str).str.replace(r"(_)?-P\d.*$", "", regex=True)
df["log2FC"] = pd.to_numeric(df["neg|lfc"], errors="coerce")
df["directional_p"] = np.where(
    df["log2FC"] < 0,
    pd.to_numeric(df["neg|p-value"], errors="coerce"),
    np.where(df["log2FC"] > 0, pd.to_numeric(df["pos|p-value"], errors="coerce"),
             df[["neg|p-value", "pos|p-value"]].apply(pd.to_numeric, errors="coerce").min(axis=1))
)
df["directional_fdr"] = np.where(
    df["log2FC"] < 0,
    pd.to_numeric(df["neg|fdr"], errors="coerce"),
    np.where(df["log2FC"] > 0, pd.to_numeric(df["pos|fdr"], errors="coerce"),
             df[["neg|fdr", "pos|fdr"]].apply(pd.to_numeric, errors="coerce").min(axis=1))
)
df["neg_log10_fdr"] = -np.log10(np.clip(df["directional_fdr"].astype(float), 1e-300, None))
df["selection_class"] = "Not significant"
df.loc[(df["directional_fdr"] < FDR_CUTOFF) & (df["log2FC"] < 0), "selection_class"] = "Negative-selection hit"
df.loc[(df["directional_fdr"] < FDR_CUTOFF) & (df["log2FC"] > 0), "selection_class"] = "Positive-selection hit"

# Save processed data.
source_cols = ["id", "gene", "log2FC", "directional_p", "directional_fdr", "neg_log10_fdr", "selection_class", "neg|rank", "pos|rank", "neg|goodsgrna", "pos|goodsgrna"]
df[source_cols].to_csv(os.path.join(OUTDIR, "volcano_processed_source.csv"), index=False)
summary = pd.DataFrame([
    {"metric": "n_total", "value": int(len(df))},
    {"metric": "n_negative_selection_FDR_lt_0.05", "value": int((df["selection_class"] == "Negative-selection hit").sum())},
    {"metric": "n_positive_selection_FDR_lt_0.05", "value": int((df["selection_class"] == "Positive-selection hit").sum())},
    {"metric": "FDR_cutoff", "value": FDR_CUTOFF},
])
summary.to_csv(os.path.join(OUTDIR, "volcano_summary_counts.csv"), index=False)

# Label manually selected interpretable / top hits to avoid overcrowding.
label_genes = ["KDR", "ACVR2A", "SMAD4", "BMPR1A", "NF2", "CDH5", "ETV2", "LATS2", "AMOTL2", "SAV1", "SKA3", "CBLL1", "RPS15A", "SMC2", "CNOT2"]
label_df = df[df["gene"].isin(label_genes)].drop_duplicates("gene")
label_df[source_cols].to_csv(os.path.join(OUTDIR, "volcano_labeled_genes.csv"), index=False)

# Print-size figure: Nature double-column width is 183 mm ~= 7.2 in.
fig, ax = plt.subplots(figsize=(7.20, 4.75), dpi=350)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Background and significant points.
plot_order = ["Not significant", "Negative-selection hit", "Positive-selection hit"]
colors = {"Not significant": COL_NON, "Negative-selection hit": COL_NEG, "Positive-selection hit": COL_POS}
legend_labels = {
    "Not significant": "Not significant",
    "Negative-selection hit": "Negative selection, FDR < 0.05",
    "Positive-selection hit": "Positive selection, FDR < 0.05",
}
for cls in plot_order:
    sub = df[df["selection_class"] == cls]
    ax.scatter(
        sub["log2FC"], sub["neg_log10_fdr"],
        s=4.5 if cls == "Not significant" else 9,
        c=colors[cls],
        alpha=0.35 if cls == "Not significant" else 0.9,
        linewidths=0,
        label=f"{legend_labels[cls]} (n={len(sub):,})",
        rasterized=(cls == "Not significant"),
        zorder=1 if cls == "Not significant" else 3,
    )

ax.axhline(-np.log10(FDR_CUTOFF), color="#6E6E6E", lw=0.65, ls=(0, (3, 2)), zorder=0)
ax.axvline(0, color="#6E6E6E", lw=0.65, ls=(0, (3, 2)), zorder=0)
ax.text(0.995, -np.log10(FDR_CUTOFF) + 0.055, "FDR = 0.05", transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=5.4, color="#555555")

# Manual label target positions for a cleaner final figure.
manual_xytext = {
    "KDR":     (-3.56, 4.44),
    "ACVR2A": (-3.22, 4.67),
    "SMAD4":  (-3.03, 4.50),
    "BMPR1A": (-2.79, 4.32),
    "NF2":    (-2.38, 4.63),
    "CDH5":   (-2.08, 4.50),
    "ETV2":   (-1.64, 4.63),
    "LATS2":  (-1.10, 4.72),
    "AMOTL2": (-0.86, 4.47),
    "SAV1":   (-0.57, 2.43),
    "SKA3":   (2.63, 2.74),
    "CBLL1":  (1.47, 2.74),
    "RPS15A": (2.12, 1.98),
    "SMC2":   (1.27, 1.94),
    "CNOT2":  (1.02, 1.75),
}
for _, row in label_df.iterrows():
    gene = row["gene"]
    xy = (row["log2FC"], row["neg_log10_fdr"])
    xytext = manual_xytext.get(gene, (row["log2FC"], row["neg_log10_fdr"] + 0.12))
    ax.annotate(
        gene, xy=xy, xytext=xytext,
        fontsize=5.25, color="#222222",
        ha="left", va="center",
        arrowprops=dict(arrowstyle="-", color="#8A8A8A", lw=0.35, shrinkA=0, shrinkB=2, alpha=0.75),
        zorder=5,
    )

x_abs = max(3.7, np.nanpercentile(np.abs(df["log2FC"]), 99.8) * 1.05)
ax.set_xlim(-x_abs, x_abs)
ax.set_ylim(-0.05, 5.25)
ax.set_xlabel("log$_2$ fold change (positive vs negative)")
ax.set_ylabel("−log$_{10}$(FDR)")
ax.set_title("Genome-wide screen volcano plot", loc="left", fontweight="bold", pad=6)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
ax.spines["left"].set_linewidth(0.6)
ax.spines["bottom"].set_linewidth(0.6)
ax.tick_params(axis="both", width=0.55, length=2.5, pad=2)
ax.legend(loc="upper right", frameon=False, handletextpad=0.4, borderpad=0.2, labelspacing=0.45, markerscale=1.2)
ax.text(0.015, 0.985, f"{len(df):,} perturbation groups", transform=ax.transAxes, ha="left", va="top", fontsize=6, color="#333333")
fig.tight_layout(pad=0.7)

out_pdf = os.path.join(OUTDIR, "genomewide_screen_volcano_nature_style_final_AI_editable.pdf")
out_svg = os.path.join(OUTDIR, "genomewide_screen_volcano_nature_style_final_AI_editable.svg")
out_png = os.path.join(OUTDIR, "genomewide_screen_volcano_nature_style_final.png")
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_svg, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=600)
plt.close(fig)

zip_path = os.path.join("/mnt/data", "genomewide_screen_volcano_nature_final_package.zip")
with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
    for fn in [out_pdf, out_svg, out_png, os.path.join(OUTDIR, "volcano_processed_source.csv"), os.path.join(OUTDIR, "volcano_labeled_genes.csv"), os.path.join(OUTDIR, "volcano_summary_counts.csv"), os.path.abspath(__file__)]:
        z.write(fn, arcname=os.path.basename(fn))

print("Saved:")
for fn in [out_pdf, out_svg, out_png, os.path.join(OUTDIR, "volcano_processed_source.csv"), os.path.join(OUTDIR, "volcano_labeled_genes.csv"), os.path.join(OUTDIR, "volcano_summary_counts.csv"), zip_path]:
    print(fn)
