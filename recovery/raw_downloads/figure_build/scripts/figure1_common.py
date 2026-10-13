from __future__ import annotations

from pathlib import Path
import math
import textwrap
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, Rectangle, Polygon
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
CLEAN = ROOT / "cleaned_data"
CANDIDATE = ROOT / "candidate_panels"
PANELS = ROOT / "individual_panels"
ASSEMBLED = ROOT / "assembled_figures"
SOURCE_DATA = ROOT / "source_data"
LEGENDS = ROOT / "legends"

for d in [CANDIDATE, PANELS, ASSEMBLED, SOURCE_DATA, LEGENDS]:
    d.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "Liberation Sans",
    "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
})

# Stable semantic palette.
COL_REQUIRED = "#0072B2"      # blue: knockdown reduces EC differentiation
COL_RESTRICT = "#D55E00"      # vermillion: knockdown shifts toward EC
COL_NEUTRAL = "#B7B7B7"
COL_DARK = "#222222"
COL_TEMP = "#CC6677"
COL_PANEL_BG = "#FAFAFA"
COL_DAY = {"D0": "#6E6E6E", "D1": "#E69F00", "D2": "#56B4E9", "D3": "#009E73"}


def load_screen() -> pd.DataFrame:
    return pd.read_csv(CLEAN / "gw_screen_hits_244_gene_level.tsv", sep="\t")


def _bbox(ax, xy, w, h, fc="white", ec="#555555", lw=0.8, radius=0.08, alpha=1.0, z=2):
    x, y = xy
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={radius}",
        facecolor=fc, edgecolor=ec, linewidth=lw, alpha=alpha, zorder=z,
    )
    ax.add_patch(p)
    return p


def _arrow(ax, start, end, color="#333333", lw=1.2, style="-|>", mutation=10, z=3, connectionstyle="arc3,rad=0"):
    a = FancyArrowPatch(
        start, end,
        arrowstyle=style,
        mutation_scale=mutation,
        linewidth=lw,
        color=color,
        connectionstyle=connectionstyle,
        shrinkA=0,
        shrinkB=0,
        zorder=z,
    )
    ax.add_patch(a)
    return a


def _temporary_tag(ax, text="TEMPORARY — incomplete Batch A inputs", x=0.99, y=0.99, compact=False):
    fs = 5.8 if compact else 7.3
    ax.text(
        x, y, text,
        transform=ax.transAxes, ha="right", va="top",
        fontsize=fs, fontweight="bold", color="#8B1A1A",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="#FFF0F0", edgecolor=COL_TEMP, linewidth=0.7),
        zorder=10,
    )


def _title(ax, title, compact=False, y=1.02):
    ax.text(
        0.0, y, title,
        transform=ax.transAxes, ha="left", va="bottom",
        fontsize=7.0 if compact else 10.5,
        fontweight="bold", color=COL_DARK,
        clip_on=False,
    )


def draw_panel_1a(ax, compact=False):
    """Conceptual genome-wide CRISPRi workflow based on manuscript-supported timing."""
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 5.2)
    ax.axis("off")
    _title(ax, "Genome-wide CRISPRi screen", compact)
    ax.text(0.0, 1.005, "Regulators of endothelial differentiation", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=5.3 if compact else 7.2, color="#555555")

    fs = 5.0 if compact else 6.7
    fs_small = 4.2 if compact else 5.7
    lw = 0.85 if compact else 1.05

    # 1. Library and inducible CRISPRi cells.
    _bbox(ax, (0.15, 2.55), 1.45, 1.30, fc="#F5F5F5", ec="#666666", lw=0.75)
    ax.text(0.875, 3.55, "Genome-wide\ngRNA library", ha="center", va="center", fontsize=fs, fontweight="bold")
    for i in range(4):
        ax.plot([0.45, 1.28], [2.82 + i*0.11, 2.82 + i*0.11], color=COL_DAY["D1"], lw=0.9)
        ax.add_patch(Circle((0.40, 2.82 + i*0.11), 0.025, color=COL_DAY["D2"], zorder=3))
    _arrow(ax, (1.70, 3.20), (2.15, 3.20), lw=lw, mutation=8 if compact else 9)
    _bbox(ax, (2.20, 2.55), 1.50, 1.30, fc="#F4EEF9", ec="#7451A6", lw=0.75)
    ax.text(2.95, 3.46, "Inducible CRISPRi\nhiPSCs", ha="center", va="center", fontsize=fs, fontweight="bold")
    ax.text(2.95, 2.78, "dCas9–KRAB", ha="center", va="center", fontsize=fs_small, color="#555555")

    # 2. Compact timeline.
    x0, x1, yline = 4.10, 8.25, 3.20
    ax.plot([x0, x1], [yline, yline], color="#555555", lw=0.9)
    events = [
        (4.25, "−4", "infect"),
        (4.95, "−2", "+dox"),
        (5.75, "D0", "hiPSC"),
        (6.52, "D1", "PS"),
        (7.28, "D2", "LM"),
        (8.05, "D3", "EC"),
    ]
    for x, top, bottom in events:
        is_day = top.startswith("D")
        col = COL_DAY[top] if is_day else "#555555"
        ax.add_patch(Circle((x, yline), 0.070 if is_day else 0.045, facecolor=col, edgecolor="white", lw=0.35, zorder=5))
        ax.plot([x, x], [yline-0.10, yline+0.10], color="#555555", lw=0.65)
        ax.text(x, 3.62, top, ha="center", va="bottom", fontsize=fs, fontweight="bold", color=col if is_day else COL_DARK)
        ax.text(x, 2.84, bottom, ha="center", va="top", fontsize=fs_small, color="#555555")
    ax.text(4.95, 2.28, "CRISPRi maintained through differentiation", ha="center", va="top",
            fontsize=fs_small, color="#555555")

    # 3. D3 sorting.
    _arrow(ax, (8.35, 3.20), (8.75, 3.20), lw=lw, mutation=8 if compact else 9)
    ax.text(9.18, 4.12, "D3 MACS", ha="center", va="center", fontsize=fs, fontweight="bold")
    _arrow(ax, (8.80, 3.20), (9.05, 3.68), color="#444444", lw=lw, mutation=8 if compact else 9)
    _arrow(ax, (8.80, 3.20), (9.05, 2.72), color="#444444", lw=lw, mutation=8 if compact else 9)
    _bbox(ax, (9.12, 3.42), 1.05, 0.62, fc="#E7F5EC", ec=COL_DAY["D3"], lw=0.75)
    _bbox(ax, (9.12, 2.36), 1.05, 0.62, fc="#F0F0F0", ec="#777777", lw=0.75)
    ax.text(9.645, 3.73, "CDH5+  EC", ha="center", va="center", fontsize=fs, fontweight="bold", color="#137A53")
    ax.text(9.645, 2.67, "CDH5−  non-EC", ha="center", va="center", fontsize=fs, fontweight="bold", color="#555555")

    # 4. Sequencing/readout.
    _arrow(ax, (10.25, 3.20), (10.65, 3.20), color="#444444", lw=lw, mutation=8 if compact else 9)
    _bbox(ax, (10.72, 2.53), 1.10, 1.34, fc="#F7F7F7", ec="#555555", lw=0.75)
    ax.text(11.27, 3.48, "gRNA\nsequencing", ha="center", va="center", fontsize=fs, fontweight="bold")
    ax.text(11.27, 2.79, "regulator\nidentification", ha="center", va="center", fontsize=fs_small, color="#444444")

    ax.text(0.18, 1.25,
            "Infect at day −4; induce dCas9–KRAB with doxycycline at day −2; sort CDH5+ and CDH5− populations on day 3.",
            ha="left", va="top", fontsize=fs_small, color="#4A4A4A")

def draw_candidate_a(ax, screen: pd.DataFrame, compact=False):
    x = screen["u_test_effect_log2_CDH5pos_vs_CDH5neg"]
    y = screen["minus_log10_u_test_fdr"]
    direction = screen["screen_direction"]
    colors = np.where(direction.eq("Required for EC differentiation"), COL_REQUIRED, COL_RESTRICT)
    ax.scatter(x, y, s=15 if compact else 28, c=colors, alpha=0.82, edgecolors="white", linewidths=0.25)
    ax.axvline(0, color="#999999", lw=0.7, ls="--")
    ax.axhline(-math.log10(0.05), color="#999999", lw=0.7, ls=":")
    for gene in ["CDH5", "KDR", "ETV2"]:
        row = screen.loc[screen["gene"].eq(gene)]
        if len(row):
            r = row.iloc[0]
            ax.text(r["u_test_effect_log2_CDH5pos_vs_CDH5neg"], r["minus_log10_u_test_fdr"], gene, fontsize=5.8 if compact else 8.0, fontweight="bold")
    ax.set_xlabel("U-test effect: log2 enrichment in CDH5+ vs CDH5−", fontsize=6 if compact else 8)
    ax.set_ylabel("−log10(U-test FDR)", fontsize=6 if compact else 8)
    ax.tick_params(labelsize=5 if compact else 7)
    ax.spines[["top", "right"]].set_visible(False)
    _title(ax, "Candidate A — coherent U-test volcano (available hits only)", compact)
    _temporary_tag(ax, "TEMPORARY — full tested set missing", compact=compact)


def draw_candidate_b(ax, compact=False):
    ax.set_xlim(-3, 3)
    ax.set_ylim(0, 6)
    ax.axvline(0, color="#C0C0C0", ls="--", lw=0.8)
    ax.axhline(-math.log10(0.05), color="#C0C0C0", ls=":", lw=0.8)
    ax.set_xlabel("MaGECK effect / rank statistic", fontsize=6 if compact else 8)
    ax.set_ylabel("−log10(MaGECK FDR)", fontsize=6 if compact else 8)
    ax.tick_params(labelsize=5 if compact else 7)
    ax.spines[["top", "right"]].set_visible(False)
    _title(ax, "Candidate B — coherent MaGECK result plot", compact)
    ax.text(
        0.5, 0.55,
        "TEMPORARY PLACEHOLDER\n\nPer-target MaGECK effect/rank and FDR\nwere not accessible in Batch A.",
        transform=ax.transAxes, ha="center", va="center",
        fontsize=7 if compact else 10, color="#7A1A1A", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFF3F3", edgecolor=COL_TEMP, linewidth=1.0, linestyle="--"),
    )


def draw_candidate_c(ax, screen: pd.DataFrame, compact=False, selected=False):
    data = screen.sort_values("u_test_effect_log2_CDH5pos_vs_CDH5neg").reset_index(drop=True).copy()
    data["rank"] = np.arange(1, len(data)+1)
    neg = data["u_test_effect_log2_CDH5pos_vs_CDH5neg"] < 0
    ax.scatter(data.loc[neg, "rank"], data.loc[neg, "u_test_effect_log2_CDH5pos_vs_CDH5neg"],
               s=11 if compact else 22, c=COL_REQUIRED, alpha=0.78, edgecolors="none")
    ax.scatter(data.loc[~neg, "rank"], data.loc[~neg, "u_test_effect_log2_CDH5pos_vs_CDH5neg"],
               s=18 if compact else 32, c=COL_RESTRICT, alpha=0.90, edgecolors="white", linewidths=0.25)
    ax.axhline(0, color="#888888", lw=0.7)
    texts = []
    for gene in ["KDR", "CDH5", "ETV2"]:
        row = data.loc[data["gene"].eq(gene)]
        if len(row):
            r = row.iloc[0]
            texts.append(ax.text(r["rank"], r["u_test_effect_log2_CDH5pos_vs_CDH5neg"], gene,
                                 fontsize=5.8 if compact else 8.0, fontweight="bold", color=COL_DARK))
    if texts and not compact:
        # Deterministic small offsets; only three predefined labels are shown.
        for i, txt in enumerate(texts):
            x0, y0 = txt.get_position()
            txt.set_position((x0 + (3 if i % 2 == 0 else -3), y0 + 0.10 + 0.06*i))
            txt.set_ha("left" if i % 2 == 0 else "right")
    ax.set_xlabel("244 final screen-derived regulators ranked by U-test effect", fontsize=6 if compact else 8)
    ax.set_ylabel("U-test effect\nlog2(CDH5+ / CDH5−)", fontsize=6 if compact else 8)
    ax.tick_params(labelsize=5 if compact else 7)
    ax.spines[["top", "right"]].set_visible(False)
    _title(ax, "Screen-derived regulators ranked by coherent U-test effect", compact)
    _temporary_tag(ax, "TEMPORARY — hit-only display", x=0.99, y=0.16, compact=compact)
    # Direction labels.
    ax.text(0.01, 0.04, "Required for EC differentiation", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=5.2 if compact else 7.0, color=COL_REQUIRED, fontweight="bold")
    ax.text(0.99, 0.95, "Restricts EC differentiation", transform=ax.transAxes,
            ha="right", va="top", fontsize=5.2 if compact else 7.0, color=COL_RESTRICT, fontweight="bold")
    if selected:
        ax.text(0.01, 0.95, "Selected temporary main-text representation", transform=ax.transAxes,
                ha="left", va="top", fontsize=5.0 if compact else 6.8, color="#555555", style="italic")


def draw_panel_1c(ax, screen: pd.DataFrame, compact=False):
    subset = screen.loc[screen["gene"].isin(["KDR", "CDH5", "ETV2"])].copy()
    subset = subset.sort_values("u_test_effect_log2_CDH5pos_vs_CDH5neg", ascending=True)
    y = np.arange(len(subset))
    ax.hlines(y, subset["u_test_effect_log2_CDH5pos_vs_CDH5neg"], 0, color="#A7CBE1", lw=2.0 if compact else 3.0)
    ax.scatter(subset["u_test_effect_log2_CDH5pos_vs_CDH5neg"], y, s=36 if compact else 70, c=COL_REQUIRED, edgecolors="white", linewidths=0.6, zorder=3)
    ax.axvline(0, color="#777777", lw=0.7)
    ax.set_yticks(y, subset["gene"])
    ax.set_xlim(min(-3.35, subset["u_test_effect_log2_CDH5pos_vs_CDH5neg"].min()-0.2), 0.15)
    ax.set_xlabel("U-test log2(CDH5+ / CDH5−)", fontsize=5.7 if compact else 8)
    ax.tick_params(labelsize=5.5 if compact else 8)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    _title(ax, "Recovery of known EC regulators" if compact else "Known EC regulators rank among strongest required hits", compact)
    # rank labels
    for yi, (_, r) in zip(y, subset.iterrows()):
        ax.text(0.02, yi, f"rank {int(r['rank_by_u_test_effect_ascending'])} / 244",
                ha="left", va="center", fontsize=5.0 if compact else 7.0, color="#555555")
    ax.text(0.0, -0.28, "Limited manuscript-defined set; final author-approved regulator list pending.",
            transform=ax.transAxes, ha="left", va="top", fontsize=4.8 if compact else 6.6,
            color="#666666", style="italic", clip_on=False)


def draw_panel_1d_placeholder(ax, compact=False):
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-0.7, 5.7)
    cats = ["EC / mesoderm", "Signaling", "Chromatin", "RNA / translation", "Cellular processes"]
    ypos = np.arange(5)[::-1]
    # faint diverging-bar scaffold only; not data.
    left = np.array([0.55, 0.80, 0.48, 0.68, 0.90])
    right = np.array([0.18, 0.25, 0.12, 0.20, 0.16])
    ax.barh(ypos, -left, color="#DDECF4", edgecolor="#B6D5E6", height=0.55)
    ax.barh(ypos, right, color="#FBE4D8", edgecolor="#F0BFA6", height=0.55)
    ax.axvline(0, color="#777777", lw=0.7)
    ax.set_yticks(ypos, cats)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=5 if compact else 7, length=0)
    ax.spines[:].set_visible(False)
    _title(ax, "Functional organization of 244 regulators" if compact else "Functional organization of the 244 screen-derived regulators", compact)
    ax.text(1.08, 5.28, "n = 244 unique genes", ha="right", va="center",
            fontsize=5.2 if compact else 7.0, color="#555555")
    # overlay placeholder notice
    ax.text(
        0.0, 2.35,
        "TEMPORARY PLACEHOLDER\nFinal author-reviewed 244-gene categories not attached\nReview template generated:\nscreen_regulators_244_annotation_template.tsv",
        ha="center", va="center",
        fontsize=5.0 if compact else 7.3, color="#7A1A1A", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="white", edgecolor=COL_TEMP, linewidth=1.0, linestyle="--", alpha=0.96),
        zorder=5,
    )


def draw_panel_1e_placeholder(ax, compact=False):
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-3.2, 3.2)
    ax.axhline(0, color="#C4C4C4", lw=0.7, ls="--")
    ax.axvline(0, color="#C4C4C4", lw=0.7, ls="--")
    ax.plot([-3, 3], [-3, 3], color="#D0D0D0", lw=0.9, ls=":")
    ax.set_xlabel("Genome-wide screen effect", fontsize=5.7 if compact else 8)
    ax.set_ylabel("Validation-screen effect", fontsize=5.7 if compact else 8)
    ax.tick_params(labelsize=5 if compact else 7)
    ax.spines[["top", "right"]].set_visible(False)
    _title(ax, "Secondary validation and reproducibility" if compact else "Secondary screening validates genome-wide effects", compact)
    # replicate inset placeholder
    inset = ax.inset_axes([0.62, 0.55, 0.32, 0.32])
    inset.set_xlim(-2, 2); inset.set_ylim(-2, 2)
    inset.axhline(0, color="#D0D0D0", lw=0.5, ls="--")
    inset.axvline(0, color="#D0D0D0", lw=0.5, ls="--")
    inset.set_xticks([]); inset.set_yticks([])
    inset.set_title("replicate 1 vs 2", fontsize=4.6 if compact else 6.2)
    for sp in inset.spines.values():
        sp.set_color("#BBBBBB")
    ax.text(
        0.45, 0.38,
        "TEMPORARY PLACEHOLDER\nValidation and replicate-level tables not attached\nRequired: gene, validation effect, p/q, replicate 1 and 2",
        transform=ax.transAxes, ha="center", va="center",
        fontsize=5.0 if compact else 7.2, color="#7A1A1A", fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.42", facecolor="white", edgecolor=COL_TEMP, linewidth=1.0, linestyle="--", alpha=0.96),
    )


def panel_letter(ax, letter, compact=True):
    ax.text(-0.065, 1.035, letter, transform=ax.transAxes,
            ha="left", va="bottom", fontsize=9.5 if compact else 13,
            fontweight="bold", color="#111111", clip_on=False)


def export_figure(fig, stem: Path, transparent_png=True, dpi=600):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", transparent=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", transparent=True)
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight", dpi=dpi, transparent=transparent_png)

