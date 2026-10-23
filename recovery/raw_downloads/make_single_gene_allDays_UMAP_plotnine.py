
import os, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

from plotnine import (
    ggplot, aes, geom_point, scale_fill_gradient2, scale_size_continuous,
    coord_equal, labs, theme_bw, theme, element_text, element_blank, element_rect
)

BASE = Path("/mnt/data")
OUTDIR = BASE / "carrier_seeded_single_gene_allDays_UMAP_plotnine"
OUTDIR.mkdir(parents=True, exist_ok=True)

BACKGROUND_CSV = BASE / "all_days_background_cells_for_plotting.csv"
GENES = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]

# Plot settings
LOG2FC_CLIP = 3.0
MAX_BACKGROUND_CELLS = 60000  # stratified subsample for clean rendering; all days are represented
RANDOM_SEED = 17

bg = pd.read_csv(BACKGROUND_CSV)
bg = bg[["UMAP1", "UMAP2", "day_clean"]].copy()

frames = []
for gene in GENES:
    f = BASE / f"carrier_seeded_scores_{gene}_all_carriers.csv"
    df = pd.read_csv(f)
    df["target"] = gene
    frames.append(df)
scores = pd.concat(frames, ignore_index=True)

scores["log2OR_plot"] = scores["log2_local_enrichment"].clip(-LOG2FC_CLIP, LOG2FC_CLIP)
scores["is_sig"] = scores["significant_fdr_0_05"].fillna(False).astype(bool)

# TRUE carrier fraction in neighborhood:
# local target carriers divided by all kNN cells in that seed-cell neighborhood.
scores["carrier_fraction_pct"] = 100 * scores["local_target_n"] / scores["k_used"].replace(0, np.nan)
scores["carrier_fraction_pct"] = scores["carrier_fraction_pct"].clip(lower=0)

# Stratified background sample by day for speed and readability.
rng = np.random.default_rng(RANDOM_SEED)
if len(bg) > MAX_BACKGROUND_CELLS:
    bg_parts = []
    n_days = max(1, bg["day_clean"].nunique())
    n_per_day = MAX_BACKGROUND_CELLS // n_days
    for day, sub in bg.groupby("day_clean"):
        n = min(len(sub), n_per_day)
        idx = rng.choice(sub.index.values, size=n, replace=False) if len(sub) > n else sub.index.values
        bg_parts.append(bg.loc[idx])
    bg_plot = pd.concat(bg_parts, ignore_index=True)
else:
    bg_plot = bg.copy()

# Global limits so each single-gene panel is directly comparable.
xmin = min(bg["UMAP1"].min(), scores["UMAP1"].min())
xmax = max(bg["UMAP1"].max(), scores["UMAP1"].max())
ymin = min(bg["UMAP2"].min(), scores["UMAP2"].min())
ymax = max(bg["UMAP2"].max(), scores["UMAP2"].max())

pad_x = (xmax - xmin) * 0.03
pad_y = (ymax - ymin) * 0.03
XLIM = (xmin - pad_x, xmax + pad_x)
YLIM = (ymin - pad_y, ymax + pad_y)

size_upper = max(30, float(scores["carrier_fraction_pct"].quantile(0.995)))

def make_plot(gene):
    df = scores[scores["target"] == gene].copy()
    nonsig = df[~df["is_sig"]]
    sig = df[df["is_sig"]]

    p = (
        ggplot()
        + geom_point(
            bg_plot,
            aes("UMAP1", "UMAP2"),
            color="#D0D0D0",
            size=0.22,
            alpha=0.18,
            raster=True,
            show_legend=False,
        )
        + geom_point(
            nonsig,
            aes("UMAP1", "UMAP2", fill="log2OR_plot", size="carrier_fraction_pct"),
            shape="o",
            color="#666666",
            stroke=0.30,
            alpha=0.88,
        )
        + geom_point(
            sig,
            aes("UMAP1", "UMAP2", fill="log2OR_plot", size="carrier_fraction_pct"),
            shape="o",
            color="black",
            stroke=0.75,
            alpha=0.96,
        )
        + scale_fill_gradient2(
            low="#2166AC",
            mid="#F7F7F7",
            high="#B2182B",
            midpoint=0,
            limits=(-LOG2FC_CLIP, LOG2FC_CLIP),
            breaks=[-3, -2, -1, 0, 1, 2, 3],
            name="Carrier-seeded\nenrichment\nlog2OR",
        )
        + scale_size_continuous(
            range=(1.2, 10.0),
            breaks=[2, 5, 10, 20, 30],
            limits=(0, size_upper),
            name="Carrier fraction\nin neighborhood (%)",
        )
        + coord_equal(xlim=XLIM, ylim=YLIM, expand=False)
        + labs(
            title=gene,
            subtitle="All D0-D3 cells in one shared UMAP; circles are carrier seed cells; black outline = FDR < 0.05",
            x="single-cell state embedding 1",
            y="single-cell state embedding 2",
        )
        + theme_bw(base_size=12)
        + theme(
            figure_size=(7.2, 6.3),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(color="#7f7f7f", fill=None, size=0.8),
            axis_title=element_text(size=13),
            axis_text=element_text(size=10),
            plot_title=element_text(size=22, weight="bold", ha="left"),
            plot_subtitle=element_text(size=10.5),
            legend_title=element_text(size=11, weight="bold"),
            legend_text=element_text(size=10),
            legend_key=element_blank(),
        )
    )
    return p

created = []
for gene in GENES:
    p = make_plot(gene)
    for ext in ["png", "pdf"]:
        fp = OUTDIR / f"carrier_seeded_{gene}_allDays_singleUMAP_plotnine.{ext}"
        if ext == "png":
            p.save(fp, width=7.2, height=6.3, dpi=450, verbose=False)
        else:
            p.save(fp, width=7.2, height=6.3, verbose=False)
        created.append(fp)

# Export plotting sources and summary.
scores.to_csv(OUTDIR / "source_carrier_seeded_scores_for_plotting.csv", index=False)
bg_plot.to_csv(OUTDIR / "source_background_cells_sampled_for_plotting.csv", index=False)

summary = scores.groupby("target", observed=True).agg(
    n_seed_cells=("seed_cell_id", "count"),
    median_log2_enrichment=("log2_local_enrichment", "median"),
    mean_log2_enrichment=("log2_local_enrichment", "mean"),
    median_carrier_fraction_pct=("carrier_fraction_pct", "median"),
    n_fdr05=("is_sig", "sum"),
).reset_index()
summary.to_csv(OUTDIR / "carrier_seeded_plot_summary.csv", index=False)

# Zip package.
zip_path = BASE / "carrier_seeded_single_gene_allDays_UMAP_plotnine_package.zip"
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
    for fp in OUTDIR.iterdir():
        z.write(fp, arcname=fp.name)

print("Created files:")
for fp in created:
    print(fp)
print("ZIP:", zip_path)
print(summary.to_string(index=False))
