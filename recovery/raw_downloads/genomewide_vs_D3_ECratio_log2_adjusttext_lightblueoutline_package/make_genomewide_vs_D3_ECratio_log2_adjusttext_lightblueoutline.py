from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import pearsonr
from adjustText import adjust_text

source_csv = Path(r"/mnt/data/genomewide_vs_D3_ECratio_log2_adjusttext_source.csv")
df = pd.read_csv(source_csv).copy()

df["screen_sig"] = df["screen_fdr"] < 0.05
df["d3_sig"] = df["D3_EC_fdr"] < 0.05
df["plot_group"] = "Not significant"
df.loc[df["d3_sig"] & ~df["screen_sig"], "plot_group"] = "Perturb-seq FDR < 0.05"
df.loc[df["screen_sig"] & ~df["d3_sig"], "plot_group"] = "GW screen FDR < 0.05"
df.loc[df["screen_sig"] & df["d3_sig"], "plot_group"] = "Both FDR < 0.05"

pearson_r, pearson_p = pearsonr(df["screen_log2FC"], df["D3_EC_log2_ratio"])

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 1.2,
})

colors = {'Not significant': '#C9C9C9', 'Perturb-seq FDR < 0.05': '#4F80C1', 'GW screen FDR < 0.05': '#E1A04F', 'Both FDR < 0.05': '#7B5DAE'}
order = ['Not significant', 'Perturb-seq FDR < 0.05', 'GW screen FDR < 0.05', 'Both FDR < 0.05']
sizes = {'Not significant': 18, 'Perturb-seq FDR < 0.05': 34, 'GW screen FDR < 0.05': 34, 'Both FDR < 0.05': 34}
alphas = {'Not significant': 0.6, 'Perturb-seq FDR < 0.05': 0.92, 'GW screen FDR < 0.05': 0.92, 'Both FDR < 0.05': 0.95}
outline_color = "#9FCBEA"

label_targets = ['CBLL1', 'LEMD3', 'SMAD4', 'ACVR2A', 'KDR', 'ETV2', 'BMPR1A', 'PLCG1', 'TAOK1', 'AMOTL2', 'LATS2', 'NF2', 'SAV1', 'ACVR2B', 'SMAD1', 'SMAD5']
label_df = df[df["target"].isin(label_targets)].copy()

fig, ax = plt.subplots(figsize=(6.2, 5.6), dpi=350)

for group in order:
    sub = df[df["plot_group"] == group]
    ax.scatter(
        sub["screen_log2FC"],
        sub["D3_EC_log2_ratio"],
        s=sizes[group],
        c=colors[group],
        alpha=alphas[group],
        edgecolors=outline_color,
        linewidths=0.45,
        label=group,
        zorder=2 if group == "Not significant" else 3,
    )

coef = np.polyfit(df["screen_log2FC"], df["D3_EC_log2_ratio"], 1)
xfit = np.linspace(df["screen_log2FC"].min(), df["screen_log2FC"].max(), 300)
yfit = coef[0] * xfit + coef[1]
ax.plot(xfit, yfit, color="black", lw=1.7, zorder=1)

ax.axhline(0, color="#A8A8A8", lw=1.0, ls=(0, (5, 3)), zorder=0)
ax.axvline(0, color="#A8A8A8", lw=1.0, ls=(0, (5, 3)), zorder=0)

texts = []
for _, row in label_df.iterrows():
    texts.append(ax.text(
        row["screen_log2FC"],
        row["D3_EC_log2_ratio"],
        row["target"],
        fontsize=9.5,
        color="#2B2B2B",
        zorder=5,
    ))

adjust_text(
    texts,
    ax=ax,
    arrowprops=dict(
        arrowstyle="-",
        color="#777777",
        lw=0.45,
        alpha=0.75,
        shrinkA=2,
        shrinkB=2,
    ),
    expand=(1.10, 1.18),
    force_text=(0.45, 0.55),
    force_static=(0.12, 0.18),
    force_pull=(0.03, 0.03),
    only_move={"text": "xy", "static": "xy", "pull": "xy", "explode": "xy"},
)

ax.set_title(f"n = {len(df)} targets; Pearson r = {pearson_r:.2f}", pad=10)
ax.set_xlabel(
    "Genome-wide screen:\n"
    "gRNA enrichment in CDH5$^{+}$ vs CDH5$^{-}$\n"
    "(log$_2$ fold-change)"
)
ax.set_ylabel(
    "CC-Perturb-seq: D3 EC ratio\n"
    "log$_2$(target artEC% / NT-control artEC%)"
)

xpad = 0.10 * (df["screen_log2FC"].max() - df["screen_log2FC"].min())
ypad = 0.12 * (df["D3_EC_log2_ratio"].max() - df["D3_EC_log2_ratio"].min())
ax.set_xlim(df["screen_log2FC"].min() - xpad, df["screen_log2FC"].max() + xpad)
ax.set_ylim(df["D3_EC_log2_ratio"].min() - ypad, df["D3_EC_log2_ratio"].max() + ypad)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
ax.spines["left"].set_linewidth(1.2)
ax.spines["bottom"].set_linewidth(1.2)
ax.tick_params(width=1.1, length=5)

ax.legend(frameon=False, loc="lower right", borderaxespad=0.4, handletextpad=0.4, labelspacing=0.3)
fig.tight_layout()

fig.savefig("genomewide_vs_D3_ECratio_log2_adjusttext_lightblueoutline_AI_editable.svg", bbox_inches="tight", facecolor="white")
fig.savefig("genomewide_vs_D3_ECratio_log2_adjusttext_lightblueoutline_AI_editable.pdf", bbox_inches="tight", facecolor="white")
fig.savefig("genomewide_vs_D3_ECratio_log2_adjusttext_lightblueoutline.png", dpi=350, bbox_inches="tight", facecolor="white")
plt.close(fig)
