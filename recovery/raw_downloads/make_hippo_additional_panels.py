# %% [markdown]
# # Additional Hippo figure panels: volcanoes, top-loading genes, marker-gene DE, and optional UMAPs
#
# This script is designed to run on Sherlock. It uses:
#
# - cNMF workbook: `/scratch/users/edtng/20260325_Hippo_timecourse/cNMF_50_0_2_032526_6co.xlsx`
# - PyDESeq2 contrast tables
# - CLR-matched small-molecule program effects
# - projected zero-based h5ad for optional UMAPs
#
# Outputs are SVG/PDF/PNG with editable text in Illustrator.

# %%
from pathlib import Path
import re
import math
import warnings

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr, pearsonr

warnings.filterwarnings("ignore")

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",  # text editable in Illustrator
})

# %% [markdown]
# ## Settings

# %%
BASE = Path("/scratch/users/edtng/20260325_Hippo_timecourse")
OUT_DIR = BASE / "hippo_additional_panels"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CNMF_WORKBOOK_PATH = BASE / "cNMF_50_0_2_032526_6co.xlsx"

SPECTRA_TPM_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50."
    "gene_spectra_tpm.k_50.dt_0_2.txt"
)

DE_DIR = BASE / "perturbseq_cNMF_matched_DE" / "pydeseq2_gene_DE"
CLR_EFFECT_PATH = BASE / "perturbseq_cNMF_matched_DE" / "program_score_DE_clr" / "small_molecule_projected_cNMF_CLR_matched_effects.tsv"
H5AD_PATH = BASE / "perturbseq_cNMF_projection" / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad"

# Optional cell-level UMAP panels. Set True only if H5AD_PATH has X_umap and you want UMAPs.
RUN_UMAP_PANELS = False
UMAP_MAX_CELLS_PER_CONDITION = 8000

# Selected programs from current Hippo model.
D2_PROGRAMS = ["P14", "P15", "P17", "P19", "P23"]
D3_EC_PROGRAMS = ["P2", "P3", "P8", "P11", "P42"]
D3_ALT_PROGRAMS = ["P16", "P24", "P27", "P33", "P34", "P35"]
SELECTED_PROGRAMS = D2_PROGRAMS + D3_EC_PROGRAMS + D3_ALT_PROGRAMS
KEY_PROGRAMS = ["P15", "P17", "P23", "P2", "P34", "P24", "P33"]

# Curated representative marker genes from annotation notes.
PROGRAM_MARKERS = {
    "P15": ["TAL1", "FLI1", "SOX7"],
    "P23": ["KDR", "ETV2"],
    "P17": ["KRT8", "KRT18"],
    "P2": ["MECOM", "NRP1", "EFNB2"],
    "P34": ["ANKRD1", "KRT8", "KRT18"],
    "P24": ["GATA4", "RSPO2"],
    "P33": ["CER1", "LHX1", "OTX2", "CYP26A1"],
}

# Main contrasts for small molecule DE. Keys are display names.
DE_CONTRAST_FILES = {
    "D2_LM+TDI": "PyDESeq2_TDI_at_D2_LM.tsv",
    "D2_LM+IAG": "PyDESeq2_IAG_at_D2_LM.tsv",
    "D3_AE+TDI": "PyDESeq2_TDI_at_D3_AE.tsv",
    "D3_AE+IAG": "PyDESeq2_IAG_at_D3_AE.tsv",
    "D2_LM+TDI→D3_AE+TDI": "PyDESeq2_TDI_D2_to_D3_vs_D3_AE.tsv",
    "D2_LM+IAG→D3_AE+IAG": "PyDESeq2_IAG_D2_to_D3_vs_D3_AE.tsv",
}

# Contrast to use when coloring top-loading genes for a given program.
PROGRAM_TO_MAIN_DE_CONTRAST = {
    "P15": "D2_LM+TDI",
    "P17": "D2_LM+TDI",
    "P23": "D2_LM+TDI",
    "P2": "D3_AE+TDI",
    "P34": "D3_AE+TDI",
    "P24": "D2_LM+TDI→D3_AE+TDI",
    "P33": "D2_LM+TDI→D3_AE+TDI",
    "P14": "D2_LM+TDI",
    "P19": "D2_LM+TDI",
    "P3": "D3_AE+TDI",
    "P8": "D3_AE+TDI",
    "P11": "D3_AE+TDI",
    "P42": "D3_AE+TDI",
    "P16": "D3_AE+TDI",
    "P27": "D3_AE+TDI",
    "P35": "D3_AE+TDI",
}

print("Output:", OUT_DIR)
print("cNMF workbook:", CNMF_WORKBOOK_PATH, CNMF_WORKBOOK_PATH.exists())
print("DE directory:", DE_DIR, DE_DIR.exists())
print("CLR effects:", CLR_EFFECT_PATH, CLR_EFFECT_PATH.exists())
print("Projected h5ad:", H5AD_PATH, H5AD_PATH.exists())

# %% [markdown]
# ## Helper functions

# %%
def save_figure(fig, prefix: Path):
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in ["svg", "pdf", "png"]:
        fig.savefig(prefix.with_suffix(f".{ext}"), bbox_inches="tight")


def program_to_int(p):
    return int(str(p).replace("P", ""))


def clean_gene_symbol(x):
    return str(x).split(".")[0].strip()


def get_de_gene_col(df):
    for col in ["gene", "Gene", "gene_symbol", "symbol", "gene_name", "response_symbol", "index"]:
        if col in df.columns:
            return col
    # PyDESeq2 often writes gene IDs/symbols as unnamed first column.
    first = df.columns[0]
    if str(first).startswith("Unnamed"):
        return first
    return first


def standardize_de_table(path: Path):
    df = pd.read_csv(path, sep="\t")
    gene_col = get_de_gene_col(df)
    df = df.rename(columns={gene_col: "gene"})
    df["gene"] = df["gene"].astype(str).map(clean_gene_symbol)
    if "log2FoldChange" not in df.columns:
        for c in ["log2FC", "log2_fold_change", "lfc"]:
            if c in df.columns:
                df = df.rename(columns={c: "log2FoldChange"})
                break
    if "padj" not in df.columns:
        for c in ["adj_pval", "qvalue", "FDR"]:
            if c in df.columns:
                df = df.rename(columns={c: "padj"})
                break
    if "pvalue" not in df.columns and "pval" in df.columns:
        df = df.rename(columns={"pval": "pvalue"})
    for c in ["log2FoldChange", "padj", "pvalue", "stat"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_all_de_tables():
    out = {}
    for name, fname in DE_CONTRAST_FILES.items():
        p = DE_DIR / fname
        if p.exists():
            out[name] = standardize_de_table(p)
        else:
            print(f"Missing DE file for {name}: {p}")
    return out


def load_program_loadings():
    loadings = pd.read_excel(CNMF_WORKBOOK_PATH, sheet_name="Program Loadings")
    loadings = loadings.rename(columns={c: str(c).strip() for c in loadings.columns})
    loadings["program"] = "P" + loadings["Program"].astype(int).astype(str)
    loadings["Gene"] = loadings["Gene"].astype(str).map(clean_gene_symbol)
    loadings["Rank"] = pd.to_numeric(loadings["Rank"], errors="coerce")
    return loadings


def load_spectra():
    if not SPECTRA_TPM_PATH.exists():
        print("Spectra file not found; top-loading barplots will use ranks instead of weights.")
        return None
    spectra = pd.read_csv(SPECTRA_TPM_PATH, sep="\t", index_col=0)
    # source rows are 1..50; convert to zero-based P0..P49 if applicable
    row_int = [int(float(x)) for x in spectra.index.astype(str)]
    if row_int == list(range(1, len(row_int) + 1)):
        spectra.index = [f"P{i}" for i in range(len(row_int))]
    else:
        spectra.index = [f"P{i}" for i in row_int]
    spectra.columns = spectra.columns.astype(str).map(clean_gene_symbol)
    spectra = spectra.apply(pd.to_numeric, errors="coerce").fillna(0)
    return spectra


def add_de_to_gene_table(gene_df, de_df):
    if de_df is None:
        gene_df["log2FoldChange"] = np.nan
        gene_df["padj"] = np.nan
        return gene_df
    sub = de_df[["gene", "log2FoldChange", "padj"]].copy()
    return gene_df.merge(sub, on="gene", how="left")


def star_from_padj(p):
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def signed_neglog10_padj(df):
    padj = df["padj"].clip(lower=1e-300)
    return np.sign(df["log2FoldChange"]) * -np.log10(padj)


de_tables = load_all_de_tables()
program_loadings = load_program_loadings()
spectra = load_spectra()
print("Loaded DE tables:", list(de_tables))
print("Program loadings:", program_loadings.shape)
print("Spectra:", None if spectra is None else spectra.shape)

# %% [markdown]
# ## Panel 1: volcano plots with cNMF marker genes highlighted

# %%
def plot_volcano_grid(de_tables, marker_genes_by_program, contrasts_to_plot=None):
    if contrasts_to_plot is None:
        contrasts_to_plot = list(de_tables.keys())
    contrasts_to_plot = [c for c in contrasts_to_plot if c in de_tables]
    n = len(contrasts_to_plot)
    if n == 0:
        print("No DE tables available for volcano plot.")
        return

    highlight_genes = sorted(set(g for genes in marker_genes_by_program.values() for g in genes))
    highlight_gene_to_programs = {}
    for p, genes in marker_genes_by_program.items():
        for g in genes:
            highlight_gene_to_programs.setdefault(g, []).append(p)

    ncols = 3
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.5 * nrows), squeeze=False)

    for ax, contrast in zip(axes.ravel(), contrasts_to_plot):
        df = de_tables[contrast].copy()
        df = df.dropna(subset=["log2FoldChange", "padj"])
        df["neglog10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
        sig = (df["padj"] < 0.05) & (df["log2FoldChange"].abs() >= 0.5)

        # background points are rasterized to keep the SVG compact; highlighted genes and text are vector/editable.
        ax.scatter(
            df.loc[~sig, "log2FoldChange"],
            df.loc[~sig, "neglog10_padj"],
            s=2,
            alpha=0.25,
            rasterized=True,
            linewidths=0,
        )
        ax.scatter(
            df.loc[sig, "log2FoldChange"],
            df.loc[sig, "neglog10_padj"],
            s=3,
            alpha=0.45,
            rasterized=True,
            linewidths=0,
        )

        hdf = df[df["gene"].isin(highlight_genes)].copy()
        ax.scatter(
            hdf["log2FoldChange"],
            hdf["neglog10_padj"],
            s=18,
            linewidths=0.5,
            edgecolors="black",
            zorder=5,
        )
        for _, r in hdf.iterrows():
            label = r["gene"]
            programs = "/".join(highlight_gene_to_programs.get(label, []))
            ax.text(
                r["log2FoldChange"],
                r["neglog10_padj"],
                f" {label}",
                fontsize=5,
                va="center",
                zorder=6,
            )

        ax.axvline(0, linewidth=0.5)
        ax.axhline(-np.log10(0.05), linewidth=0.5, linestyle="--")
        ax.set_title(contrast)
        ax.set_xlabel("log2FC treatment/control")
        ax.set_ylabel("-log10 adj. p")

    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "Panel_volcano_marker_genes_key_contrasts")
    plt.show()

plot_volcano_grid(de_tables, PROGRAM_MARKERS)

# %% [markdown]
# ## Panel 2: top-loading gene barplots, colored by DE log2FC in relevant contrast

# %%
def get_top_loading_table(program, top_n=15):
    if spectra is not None and program in spectra.index:
        s = spectra.loc[program].sort_values(ascending=False).head(top_n)
        df = pd.DataFrame({"gene": s.index, "loading": s.values})
        df["rank"] = np.arange(1, df.shape[0] + 1)
    else:
        sub = program_loadings[program_loadings["program"] == program].sort_values("Rank").head(top_n)
        df = pd.DataFrame({
            "gene": sub["Gene"].to_numpy(),
            "loading": top_n + 1 - np.arange(1, sub.shape[0] + 1),
            "rank": sub["Rank"].to_numpy(),
        })
    return df


def plot_top_loading_gene_bars(programs=KEY_PROGRAMS, top_n=12):
    n = len(programs)
    ncols = 2
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.5, max(3, 2.2 * nrows)), squeeze=False)

    source_tables = []
    for ax, program in zip(axes.ravel(), programs):
        tab = get_top_loading_table(program, top_n=top_n)
        contrast = PROGRAM_TO_MAIN_DE_CONTRAST.get(program)
        de = de_tables.get(contrast)
        tab = add_de_to_gene_table(tab, de)
        tab["program"] = program
        tab["de_contrast"] = contrast
        source_tables.append(tab)

        # Sort so highest loading appears at top.
        plot_tab = tab.iloc[::-1].copy()
        values = plot_tab["loading"].to_numpy()
        colors = plot_tab["log2FoldChange"].to_numpy()
        cmap = plt.get_cmap("RdBu_r")
        if np.isfinite(colors).any():
            vmax = np.nanmax(np.abs(colors))
            norm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
            bar_colors = cmap(norm(colors))
        else:
            bar_colors = None

        ax.barh(plot_tab["gene"], values, color=bar_colors)
        ax.set_title(f"{program} top loading genes\nDE: {contrast}", fontsize=7)
        ax.set_xlabel("cNMF loading" if spectra is not None else "rank score")
        ax.tick_params(axis="y", labelsize=6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    for ax in axes.ravel()[n:]:
        ax.axis("off")

    # shared colorbar for DE log2FC
    sm = mpl.cm.ScalarMappable(cmap=plt.get_cmap("RdBu_r"), norm=mpl.colors.Normalize(vmin=-3, vmax=3))
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), fraction=0.015, pad=0.02)
    cbar.set_label("PyDESeq2 log2FC in selected contrast")

    fig.tight_layout()
    save_figure(fig, OUT_DIR / "Panel_top_loading_gene_bars_colored_by_DE")
    plt.show()

    src = pd.concat(source_tables, ignore_index=True)
    src.to_csv(OUT_DIR / "Panel_top_loading_gene_bars_source.tsv", sep="\t", index=False)

plot_top_loading_gene_bars(KEY_PROGRAMS, top_n=12)

# %% [markdown]
# ## Panel 3: marker-gene DE dot plot
#
# Rows are curated representative marker genes, columns are treatment contrasts.
# Color is PyDESeq2 log2FC; dot size is -log10 adjusted p-value.

# %%
def build_marker_gene_de_matrix():
    rows = []
    for program, genes in PROGRAM_MARKERS.items():
        for gene in genes:
            rows.append({"program": program, "gene": gene})
    marker_df = pd.DataFrame(rows).drop_duplicates()

    recs = []
    for contrast, de in de_tables.items():
        sub = de[["gene", "log2FoldChange", "padj"]].copy()
        tmp = marker_df.merge(sub, on="gene", how="left")
        tmp["contrast"] = contrast
        recs.append(tmp)
    out = pd.concat(recs, ignore_index=True)
    out["neglog10_padj"] = -np.log10(out["padj"].clip(lower=1e-300))
    out["star"] = out["padj"].map(star_from_padj)
    return out


def plot_marker_gene_de_dotplot():
    df = build_marker_gene_de_matrix()
    df.to_csv(OUT_DIR / "Panel_marker_gene_DE_dotplot_source.tsv", sep="\t", index=False)

    # Define row and column order.
    row_order = []
    for p, genes in PROGRAM_MARKERS.items():
        for g in genes:
            key = f"{p} | {g}"
            if key not in row_order:
                row_order.append(key)
    col_order = list(DE_CONTRAST_FILES.keys())
    col_order = [c for c in col_order if c in df["contrast"].unique()]

    df["row_label"] = df["program"] + " | " + df["gene"]
    df["x"] = df["contrast"].map({c:i for i,c in enumerate(col_order)})
    df["y"] = df["row_label"].map({r:i for i,r in enumerate(row_order)})
    df = df.dropna(subset=["x", "y"])

    fig, ax = plt.subplots(figsize=(max(7, 0.75 * len(col_order) + 3), max(5, 0.25 * len(row_order) + 2)))
    vmax = np.nanmax(np.abs(df["log2FoldChange"])) if np.isfinite(df["log2FoldChange"]).any() else 1
    sizes = 10 + 35 * np.clip(df["neglog10_padj"].fillna(0) / 10, 0, 1)

    sc = ax.scatter(
        df["x"], df["y"],
        c=df["log2FoldChange"],
        s=sizes,
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        edgecolors="black",
        linewidths=0.25,
    )
    for _, r in df.iterrows():
        if r["star"]:
            ax.text(r["x"], r["y"], r["star"], ha="center", va="center", fontsize=5)

    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_order, rotation=50, ha="right")
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(row_order)
    ax.set_xlim(-0.6, len(col_order) - 0.4)
    ax.set_ylim(len(row_order) - 0.5, -0.5)
    ax.set_title("Representative marker gene DE")
    ax.set_xlabel("small molecule contrast")
    ax.set_ylabel("program | marker gene")
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("PyDESeq2 log2FC")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "Panel_marker_gene_DE_dotplot")
    plt.show()

plot_marker_gene_de_dotplot()

# %% [markdown]
# ## Panel 4: GSEA-like barcode plots for program top genes on ranked DE lists

# %%
def plot_barcode_for_programs(contrast, programs=("P15", "P17", "P23", "P2", "P34", "P24"), top_n=100):
    if contrast not in de_tables:
        print(f"Missing DE table for {contrast}")
        return
    de = de_tables[contrast].dropna(subset=["log2FoldChange"]).copy()
    if "stat" in de.columns and de["stat"].notna().any():
        de["rank_score"] = de["stat"]
    elif "padj" in de.columns:
        de["rank_score"] = signed_neglog10_padj(de)
    else:
        de["rank_score"] = de["log2FoldChange"]
    de = de.sort_values("rank_score", ascending=False).reset_index(drop=True)
    rank_pos = pd.Series(np.arange(de.shape[0]), index=de["gene"].astype(str)).to_dict()

    n = len(programs)
    fig, axes = plt.subplots(n, 1, figsize=(8.0, max(2.5, 0.75 * n)), sharex=True)
    if n == 1:
        axes = [axes]

    for ax, program in zip(axes, programs):
        top = get_top_loading_table(program, top_n=top_n)["gene"].astype(str).tolist()
        xs = [rank_pos[g] for g in top if g in rank_pos]
        ax.vlines(xs, 0, 1, linewidth=0.6)
        ax.set_yticks([])
        ax.set_ylabel(program, rotation=0, ha="right", va="center")
        ax.set_xlim(0, de.shape[0])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.text(1.0, 0.5, f"{len(xs)}/{top_n}", transform=ax.transAxes, va="center", ha="right", fontsize=6)

    axes[0].set_title(f"Top cNMF loading genes along ranked DE list: {contrast}")
    axes[-1].set_xlabel("Genes ranked from up in treatment to down in treatment")
    fig.tight_layout()
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", contrast)
    save_figure(fig, OUT_DIR / f"Panel_barcode_top_loading_genes_{safe}")
    plt.show()

for c in ["D2_LM+TDI", "D3_AE+TDI", "D2_LM+TDI→D3_AE+TDI", "D2_LM+IAG"]:
    plot_barcode_for_programs(c, programs=["P15", "P17", "P23", "P2", "P34", "P24"], top_n=100)

# %% [markdown]
# ## Panel 5: loading-vs-DE scatter for selected program/contrast pairs

# %%
def plot_loading_vs_de(program_contrast_pairs=None, top_label_n=12):
    if program_contrast_pairs is None:
        program_contrast_pairs = [
            ("P15", "D2_LM+TDI"),
            ("P17", "D2_LM+TDI"),
            ("P23", "D2_LM+TDI"),
            ("P2", "D3_AE+TDI"),
            ("P34", "D3_AE+TDI"),
            ("P24", "D2_LM+TDI→D3_AE+TDI"),
        ]
    n = len(program_contrast_pairs)
    ncols = 3
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.4 * nrows), squeeze=False)

    all_records = []
    for ax, (program, contrast) in zip(axes.ravel(), program_contrast_pairs):
        if spectra is None or program not in spectra.index or contrast not in de_tables:
            ax.axis("off")
            continue
        s = spectra.loc[program]
        tab = pd.DataFrame({"gene": s.index.astype(str), "loading": s.values})
        tab = tab[tab["loading"] > 0]
        de = de_tables[contrast][["gene", "log2FoldChange", "padj"]]
        tab = tab.merge(de, on="gene", how="inner").dropna(subset=["log2FoldChange"])
        if tab.empty:
            ax.axis("off")
            continue

        r_spear, p_spear = spearmanr(tab["loading"], tab["log2FoldChange"])
        ax.scatter(tab["loading"], tab["log2FoldChange"], s=4, alpha=0.35, rasterized=True, linewidths=0)
        top = tab.sort_values("loading", ascending=False).head(top_label_n)
        ax.scatter(top["loading"], top["log2FoldChange"], s=18, edgecolors="black", linewidths=0.3)
        for _, row in top.iterrows():
            ax.text(row["loading"], row["log2FoldChange"], " " + row["gene"], fontsize=5, va="center")
        ax.axhline(0, linewidth=0.5)
        ax.set_title(f"{program} vs {contrast}\nSpearman r={r_spear:.2f}")
        ax.set_xlabel("cNMF spectra loading")
        ax.set_ylabel("DE log2FC")
        tab["program"] = program
        tab["contrast"] = contrast
        tab["spearman_r"] = r_spear
        tab["spearman_p"] = p_spear
        all_records.append(tab)

    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "Panel_loading_vs_DE_scatter_selected_programs")
    plt.show()
    if all_records:
        pd.concat(all_records, ignore_index=True).to_csv(OUT_DIR / "Panel_loading_vs_DE_scatter_source.tsv", sep="\t", index=False)

plot_loading_vs_de()

# %% [markdown]
# ## Panel 6: small-molecule vs Hippo KD selected program effect dot panel
#
# This complements the scatter/correlation plots by showing the specific selected programs.

# %%
def load_perturbseq_assoc(effect_col="approx_log2FC"):
    assoc = pd.read_excel(CNMF_WORKBOOK_PATH, sheet_name="Perturbation Association 1")
    assoc = assoc.copy()
    assoc["program"] = "P" + assoc["program_name"].astype(int).astype(str)
    assoc[effect_col] = pd.to_numeric(assoc[effect_col], errors="coerce")
    return assoc


def plot_selected_program_effect_dotpanel():
    if not CLR_EFFECT_PATH.exists():
        print("CLR effect table missing; skipping selected program effect dotpanel.")
        return
    sm = pd.read_csv(CLR_EFFECT_PATH, sep="\t")
    assoc = load_perturbseq_assoc()
    hippo_targets = ["LATS2", "NF2", "SAV1", "TAOK1", "WWC3", "AMOTL2"]
    hippo = (
        assoc[(assoc["target_name"].isin(hippo_targets)) & (assoc["Sample"].isin(["D2", "D3"])) & (assoc["program"].isin(SELECTED_PROGRAMS))]
        .groupby(["Sample", "program"], observed=True)["approx_log2FC"]
        .mean()
        .reset_index()
    )
    hippo["source"] = "Hippo CRISPRi " + hippo["Sample"].astype(str)
    hippo = hippo.rename(columns={"approx_log2FC": "effect"})

    contrasts_keep = ["D2_LM+TDI", "D2_LM+IAG", "D3_AE+TDI", "D3_AE+IAG", "D2_LM+TDI->D3_AE+TDI", "D2_LM+IAG->D3_AE+IAG"]
    sm = sm[sm["program"].isin(SELECTED_PROGRAMS) & sm["contrast"].isin(contrasts_keep)].copy()
    sm["source"] = sm["contrast"].str.replace("->", "→", regex=False)
    sm = sm.rename(columns={"mean_diff_clr": "effect"})
    plot_df = pd.concat([hippo[["source", "program", "effect"]], sm[["source", "program", "effect"]]], ignore_index=True)

    source_order = ["Hippo CRISPRi D2", "Hippo CRISPRi D3", "D2_LM+TDI", "D2_LM+IAG", "D3_AE+TDI", "D3_AE+IAG", "D2_LM+TDI→D3_AE+TDI", "D2_LM+IAG→D3_AE+IAG"]
    prog_order = SELECTED_PROGRAMS
    plot_df["x"] = plot_df["source"].map({s:i for i,s in enumerate(source_order)})
    plot_df["y"] = plot_df["program"].map({p:i for i,p in enumerate(prog_order)})
    plot_df = plot_df.dropna(subset=["x", "y"])
    vmax = np.nanmax(np.abs(plot_df["effect"]))

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    sc = ax.scatter(plot_df["x"], plot_df["y"], c=plot_df["effect"], s=70, cmap="RdBu_r", vmin=-vmax, vmax=vmax, edgecolors="black", linewidths=0.3)
    ax.set_xticks(range(len(source_order)))
    ax.set_xticklabels(source_order, rotation=50, ha="right")
    ax.set_yticks(range(len(prog_order)))
    ax.set_yticklabels(prog_order)
    ax.set_ylim(len(prog_order)-0.5, -0.5)
    ax.set_title("Selected program effects: Hippo CRISPRi vs small molecules")
    ax.set_xlabel("effect source")
    ax.set_ylabel("program")
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("effect size (programDE or CLR delta)")
    fig.tight_layout()
    save_figure(fig, OUT_DIR / "Panel_selected_program_effect_dotpanel_HippoCRISPRi_vs_small_molecule")
    plt.show()
    plot_df.to_csv(OUT_DIR / "Panel_selected_program_effect_dotpanel_source.tsv", sep="\t", index=False)

plot_selected_program_effect_dotpanel()

# %% [markdown]
# ## Optional Panel 7: UMAP feature panels for projected programs
#
# This is skipped by default. Set `RUN_UMAP_PANELS = True` at the top if your projected h5ad has `X_umap`.

# %%
def plot_umap_program_panels():
    if not RUN_UMAP_PANELS:
        print("RUN_UMAP_PANELS is False; skipping UMAP panels.")
        return
    import scanpy as sc
    adata = sc.read_h5ad(H5AD_PATH, backed="r")
    if "X_umap" not in adata.obsm:
        print("X_umap not found; skipping UMAP panels.")
        adata.file.close()
        return
    if "X_perturbseq_cNMF_usage" not in adata.obsm:
        print("X_perturbseq_cNMF_usage not found; skipping UMAP panels.")
        adata.file.close()
        return
    obs = adata.obs.copy()
    coords = np.asarray(adata.obsm["X_umap"])
    usage = np.asarray(adata.obsm["X_perturbseq_cNMF_usage"])
    meta = adata.uns.get("perturbseq_cNMF_projection", {})
    program_names = [str(x) for x in meta.get("program_names", [f"P{i}" for i in range(usage.shape[1])])]
    p_to_idx = {p:i for i,p in enumerate(program_names)}
    if "condition" not in obs.columns and "sample_name" in obs.columns:
        obs["condition"] = obs["sample_name"].astype(str).str.replace(r"_rep\d+$", "", regex=True)
    conditions = ["D2_LM", "D2_LM+TDI", "D2_LM+IAG", "D3_AE", "D3_AE+TDI", "D3_AE+IAG", "D2_LM+TDI->D3_AE+TDI", "D2_LM+IAG->D3_AE+IAG"]
    rng = np.random.default_rng(0)
    for program in ["P15", "P17", "P23", "P2", "P34", "P24"]:
        if program not in p_to_idx:
            continue
        fig, axes = plt.subplots(2, 4, figsize=(9, 4.8), squeeze=False)
        vals = usage[:, p_to_idx[program]]
        vmax = np.quantile(vals[vals > 0], 0.99) if np.any(vals > 0) else 1
        for ax, cond in zip(axes.ravel(), conditions):
            idx = np.where(obs["condition"].astype(str).to_numpy() == cond)[0]
            if idx.size > UMAP_MAX_CELLS_PER_CONDITION:
                idx = rng.choice(idx, size=UMAP_MAX_CELLS_PER_CONDITION, replace=False)
            ax.scatter(coords[idx,0], coords[idx,1], c=vals[idx], s=1, cmap="viridis", vmin=0, vmax=vmax, linewidths=0, rasterized=True)
            ax.set_title(cond, fontsize=7)
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values(): sp.set_visible(False)
        fig.suptitle(f"Projected {program} usage on UMAP", y=1.02)
        fig.tight_layout()
        save_figure(fig, OUT_DIR / f"Panel_UMAP_projected_{program}_by_condition")
        plt.show()
    adata.file.close()

plot_umap_program_panels()

# %%
print("Saved outputs to:", OUT_DIR)
for p in sorted(OUT_DIR.glob("Panel_*")):
    print(p.name)
