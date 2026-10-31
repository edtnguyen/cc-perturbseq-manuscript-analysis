from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl

# Preserve text objects in SVG and editable TrueType text in PDF.
mpl.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.usetex": False,
    "font.family": "Liberation Sans",
})

from plotnine import (
    aes, element_blank, element_rect, element_text, facet_grid,
    geom_point, geom_segment, geom_text, ggplot, guides, labs,
    scale_color_manual, scale_size_manual, scale_x_continuous,
    theme, theme_classic,
)

HERE = Path(__file__).resolve().parent
df = pd.read_csv(HERE / "PKC_isoform_genomewide_screen_log2FC.csv")

gene_order = [
    "PRKCA", "PRKCB", "PRKCG",
    "PRKCD", "PRKCE", "PRKCH", "PRKCQ",
    "PRKCI", "PRKCZ",
]
family_order = ["Conventional PKC", "Novel PKC", "Atypical PKC"]

df["gene"] = pd.Categorical(
    df["gene"], categories=list(reversed(gene_order)), ordered=True,
)
df["PKC_family"] = pd.Categorical(
    df["PKC_family"], categories=family_order, ordered=True,
)
df["significance"] = np.where(
    df["directional_FDR"] < 0.05, "FDR < 0.05", "Not significant",
)
df["point_size"] = np.where(
    df["directional_FDR"] < 0.05, "Significant", "Other",
)
df["x_start"] = 0.0
df["label"] = df["MAGeCK_log2FC_CDH5pos_vs_CDH5neg"].map(
    lambda value: f"{value:.2f}"
)
df["label_x"] = np.where(
    df["MAGeCK_log2FC_CDH5pos_vs_CDH5neg"] < 0,
    df["MAGeCK_log2FC_CDH5pos_vs_CDH5neg"] - 0.06,
    df["MAGeCK_log2FC_CDH5pos_vs_CDH5neg"] + 0.06,
)

plot = (
    ggplot(df, aes(y="gene", x="MAGeCK_log2FC_CDH5pos_vs_CDH5neg"))
    + geom_segment(
        aes(
            x="x_start",
            xend="MAGeCK_log2FC_CDH5pos_vs_CDH5neg",
            y="gene",
            yend="gene",
        ),
        color="#B5B5B5", size=0.8,
    )
    + geom_point(aes(color="significance", size="point_size"), alpha=0.98)
    + geom_text(
        aes(x="label_x", label="label"),
        size=8.2, va="center", color="#303030",
        family="Liberation Sans",
    )
    + facet_grid("PKC_family ~ .", scales="free_y", space="free_y")
    + scale_color_manual(
        values={"FDR < 0.05": "#B2182B", "Not significant": "#4D4D4D"}
    )
    + scale_size_manual(values={"Significant": 4.5, "Other": 3.0})
    + guides(size="none")
    + scale_x_continuous(
        limits=(-1.67, 0.50),
        breaks=[-1.5, -1.0, -0.5, 0.0],
        expand=(0.01, 0.01),
    )
    + labs(
        title="PKC isoforms in the genome-wide CRISPRi screen",
        subtitle="Guide enrichment in CDH5-positive versus CDH5-negative cells",
        # Plain Unicode is used instead of mathtext so Illustrator keeps every label as text.
        x="log2(sgRNA abundance in CDH5+ / CDH5- cells)",
        y=None,
        color=None,
        caption=(
            "Negative values indicate depletion from the CDH5-positive "
            "fraction. Directional MAGeCK FDR was used for significance."
        ),
    )
    + theme_classic(base_size=10.5, base_family="Liberation Sans")
    + theme(
        figure_size=(5.5, 5.4),
        text=element_text(family="Liberation Sans"),
        plot_title=element_text(size=14, weight="bold", ha="left"),
        plot_subtitle=element_text(size=9.5, color="#555555", ha="left"),
        plot_caption=element_text(size=8, color="#555555", ha="left"),
        axis_text_y=element_text(size=10),
        axis_ticks_major_y=element_blank(),
        strip_background=element_rect(fill="#F2F2F2", color="#D9D9D9", size=0.6),
        strip_text_y=element_text(size=9.5, weight="bold", angle=0),
        legend_position="top",
    )
)

for extension in ["png", "pdf", "svg"]:
    plot.save(
        HERE / f"PKC_isoform_genomewide_screen_log2FC_plotnine_editable_text.{extension}",
        format=extension,
        dpi=350 if extension == "png" else 300,
        verbose=False,
    )

# Verify that SVG labels remain native text objects.
svg_path = HERE / "PKC_isoform_genomewide_screen_log2FC_plotnine_editable_text.svg"
svg_text = svg_path.read_text(encoding="utf-8")
if "<text" not in svg_text:
    raise RuntimeError("SVG does not contain native text objects.")
print(f"Editable SVG text objects: {svg_text.count('<text')}")
