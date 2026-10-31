# %% [markdown]
# # cNMF program overlap and Hippo KD concordance analysis
#
# This notebook uses the cNMF Excel workbook:
#
# - `Program Loadings`: top 300 loading genes per cNMF program
# - `Perturbation Association 1`: effect size of each perturbation/KD on each program on each day
#
# It performs two analyses:
#
# 1. **Top program gene overlap with small-molecule DE genes**
# 2. **Program-level concordance between small-molecule effects and perturb-seq Hippo KD effects**
#
# This does not require the large h5ad/h5mu files.

# %%
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr, pearsonr
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})

# ---- Paths ----
# On Sherlock, change BASE_DIR if needed.
BASE_DIR = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE")

# Use uploaded/local paths if running this notebook beside the uploaded files.
LOCAL_DIR = Path("/mnt/data")

CNMF_XLSX = Path("/mnt/data/cNMF_50_0_2_032526_6co.xlsx")
if not CNMF_XLSX.exists():
    CNMF_XLSX = Path("/scratch/users/edtng/20260325_Hippo_timecourse/cNMF_50_0_2_032526_6co.xlsx")

# Directory containing PyDESeq2_*.tsv files.
DE_DIR = LOCAL_DIR
if not list(DE_DIR.glob("PyDESeq2_*.tsv")):
    DE_DIR = BASE_DIR / "pydeseq2_gene_DE"

# Small-molecule projected program score contrasts.
PROGRAM_SCORE_CONTRASTS = LOCAL_DIR / "projected_cNMF_program_score_matched_contrasts.tsv"
if not PROGRAM_SCORE_CONTRASTS.exists():
    PROGRAM_SCORE_CONTRASTS = BASE_DIR / "program_score_DE" / "projected_cNMF_program_score_matched_contrasts.tsv"

OUT_DIR = Path("/mnt/data/hippo_program_overlap_concordance_outputs")
if not Path("/mnt/data").exists():
    OUT_DIR = BASE_DIR / "program_overlap_concordance"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("CNMF_XLSX:", CNMF_XLSX, CNMF_XLSX.exists())
print("DE_DIR:", DE_DIR, len(list(DE_DIR.glob("PyDESeq2_*.tsv"))))
print("PROGRAM_SCORE_CONTRASTS:", PROGRAM_SCORE_CONTRASTS, PROGRAM_SCORE_CONTRASTS.exists())
print("OUT_DIR:", OUT_DIR)

# ---- Analysis parameters ----
DE_PADJ_CUTOFF = 0.05
DE_ABS_LFC_CUTOFF = 0.0
TOP_N_PROGRAM_GENES = 300

# Hippo-ish KD targets found in the workbook. Edit as needed.
# I would start conservatively with core/regulatory Hippo pathway genes that are present.
HIPPO_KD_TARGETS = ["LATS2", "NF2", "SAV1", "TAOK1", "WWC3", "AMOTL2"]

# Effect column in Perturbation Association 1.
# approx_log2FC is easiest to interpret; beta_hat is also available.
PERTURBSEQ_PROGRAM_EFFECT_COL = "approx_log2FC"

KEY_PROGRAMS = ["P15", "P17", "P34"]

# %% [markdown]
# ## Helper functions

# %%
def bh_adjust(pvals):
    """Benjamini-Hochberg FDR correction. Returns q-values in original order."""
    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan, dtype=float)
    valid = np.isfinite(pvals)
    if valid.sum() == 0:
        return out
    p = pvals[valid]
    order = np.argsort(p)
    ranked = p[order]
    n = len(ranked)
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    temp = np.empty_like(q)
    temp[order] = q
    out[valid] = temp
    return out


def normalize_program_name(x):
    """Return zero-based P# program label."""
    if isinstance(x, str) and x.startswith("P"):
        return x
    return f"P{int(float(x))}"


def contrast_to_stage_day(contrast_name):
    """Map small-molecule contrast names to perturb-seq association Sample day D0/D1/D2/D3."""
    cn = str(contrast_name)
    if "D0" in cn:
        return "D0"
    if "D1_MPS" in cn:
        return "D1"
    if "D1.33_LM" in cn or "D1.67_LM" in cn or "D2_LM" in cn:
        # LM-stage comparisons correspond most closely to D2 in perturb-seq association table.
        return "D2"
    if "D3_AE" in cn:
        return "D3"
    return None


def read_deseq_table(path):
    df = pd.read_csv(path, sep="\t")
    if "contrast_name" not in df.columns:
        name = Path(path).stem.replace("PyDESeq2_", "")
        df["contrast_name"] = name
    if "gene_symbol" in df.columns:
        df["gene_for_overlap"] = df["gene_symbol"].astype(str)
    elif "gene" in df.columns:
        df["gene_for_overlap"] = df["gene"].astype(str)
    elif "gene_id" in df.columns:
        df["gene_for_overlap"] = df["gene_id"].astype(str)
    else:
        raise ValueError(f"Could not find gene column in {path}. Columns: {list(df.columns)}")
    return df

# %% [markdown]
# ## Load cNMF workbook tabs

# %%
program_loadings = pd.read_excel(CNMF_XLSX, sheet_name="Program Loadings")
assoc = pd.read_excel(CNMF_XLSX, sheet_name="Perturbation Association 1")

print("Program Loadings shape:", program_loadings.shape)
print(program_loadings.head())
print("\nPerturbation Association 1 shape:", assoc.shape)
print(assoc.head())

# Clean columns
program_loadings["Program"] = program_loadings["Program"].astype(int)
program_loadings["Rank"] = program_loadings["Rank"].astype(int)
program_loadings["Gene"] = program_loadings["Gene"].astype(str)
program_loadings["program"] = program_loadings["Program"].map(lambda x: f"P{x}")

assoc["program"] = assoc["program_name"].map(normalize_program_name)
assoc["target_name"] = assoc["target_name"].astype(str)
assoc["Sample"] = assoc["Sample"].astype(str)

# Show Hippo-like targets found in workbook to help edit HIPPO_KD_TARGETS.
hippo_keywords = ["LATS", "STK3", "STK4", "SAV", "NF2", "MOB", "YAP", "WWTR", "TEAD", "AMOT", "TAOK", "FAT", "DCHS", "FRMD", "WWC"]
all_targets = sorted(assoc["target_name"].dropna().unique().tolist())
hippo_like_targets = [t for t in all_targets if any(k.lower() in t.lower() for k in hippo_keywords)]
print("\nHippo-like targets found in workbook:")
print(hippo_like_targets)

available_hippo_targets = [t for t in HIPPO_KD_TARGETS if t in set(all_targets)]
missing_hippo_targets = [t for t in HIPPO_KD_TARGETS if t not in set(all_targets)]
print("\nUsing HIPPO_KD_TARGETS available in workbook:", available_hippo_targets)
print("Missing from workbook:", missing_hippo_targets)
if not available_hippo_targets:
    raise ValueError("No HIPPO_KD_TARGETS were found in association table. Edit HIPPO_KD_TARGETS.")

# Top N loading gene sets
program_top_genes = {}
for program, sub in program_loadings.sort_values(["Program", "Rank"]).groupby("program"):
    program_top_genes[program] = sub.head(TOP_N_PROGRAM_GENES)["Gene"].astype(str).tolist()

print("n programs with top genes:", len(program_top_genes))
print("P15 top genes:", program_top_genes.get("P15", [])[:20])

# %% [markdown]
# ## Analysis 1: overlap between small-molecule DE genes and top cNMF loading genes
#
# This is descriptive/supportive. It asks whether top loading genes from a program are overrepresented among genes up- or down-regulated in each small molecule contrast.

# %%
de_files = sorted([p for p in DE_DIR.glob("PyDESeq2_*.tsv") if "contrast_summary" not in p.name])
print("DE files:", len(de_files))
for p in de_files[:5]:
    print(" ", p.name)

all_overlap_records = []

for de_path in de_files:
    de = read_deseq_table(de_path)
    contrast_name = str(de["contrast_name"].iloc[0])

    de = de.dropna(subset=["gene_for_overlap", "log2FoldChange", "padj"]).copy()
    de["gene_for_overlap"] = de["gene_for_overlap"].astype(str)

    background = set(de["gene_for_overlap"])
    de_up = set(de.loc[(de["padj"] < DE_PADJ_CUTOFF) & (de["log2FoldChange"] > DE_ABS_LFC_CUTOFF), "gene_for_overlap"])
    de_down = set(de.loc[(de["padj"] < DE_PADJ_CUTOFF) & (de["log2FoldChange"] < -DE_ABS_LFC_CUTOFF), "gene_for_overlap"])

    for program, genes in program_top_genes.items():
        prog_set = set(genes) & background
        for direction, de_set in [("up", de_up), ("down", de_down)]:
            overlap = prog_set & de_set
            a = len(overlap)
            b = len(prog_set) - a
            c = len(de_set) - a
            d = len(background) - a - b - c
            odds, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
            all_overlap_records.append({
                "contrast_name": contrast_name,
                "program": program,
                "direction": direction,
                "n_background": len(background),
                "n_program_top_genes_in_background": len(prog_set),
                "n_de_genes": len(de_set),
                "n_overlap": a,
                "odds_ratio": odds,
                "pvalue": pvalue,
                "overlap_genes": ";".join(sorted(overlap)),
            })

overlap_df = pd.DataFrame(all_overlap_records)
overlap_df["padj_bh_all"] = bh_adjust(overlap_df["pvalue"].to_numpy())
overlap_df["padj_bh_within_contrast"] = np.nan
for contrast, idx in overlap_df.groupby("contrast_name").groups.items():
    overlap_df.loc[idx, "padj_bh_within_contrast"] = bh_adjust(overlap_df.loc[idx, "pvalue"].to_numpy())

overlap_df = overlap_df.sort_values(["contrast_name", "padj_bh_within_contrast", "pvalue"])
out_path = OUT_DIR / "program_top300_gene_overlap_with_small_molecule_DE.tsv"
overlap_df.to_csv(out_path, sep="\t", index=False)
print("Saved:", out_path)

display(overlap_df.query("program in @KEY_PROGRAMS").sort_values(["contrast_name", "program", "direction"]).head(50))

# %% [markdown]
# ## Analysis 2: program-level concordance between small molecule effects and perturb-seq Hippo KD effects
#
# This uses the `Perturbation Association 1` tab. It compares:
#
# - small-molecule projected cNMF program effects, e.g. TDI vs control
# - perturb-seq Hippo KD effect vector across all 50 programs
#
# For TDI, we expect positive concordance with Hippo KD if TDI mimics Hippo inhibition. For IAG, we may expect weaker or opposite concordance.

# %%
sm_prog = pd.read_csv(PROGRAM_SCORE_CONTRASTS, sep="\t")
sm_prog["program"] = sm_prog["program"].astype(str)
sm_prog["stage_day_for_perturbseq"] = sm_prog["contrast_name"].map(contrast_to_stage_day)

print(sm_prog.shape)
print(sm_prog.head())
print("Stage/day mapping:")
display(sm_prog[["contrast_name", "stage_day_for_perturbseq"]].drop_duplicates().sort_values("contrast_name"))

# Average perturb-seq Hippo KD effects across selected Hippo targets for each Sample x program.
hippo_assoc = assoc.loc[assoc["target_name"].isin(available_hippo_targets)].copy()
hippo_assoc[PERTURBSEQ_PROGRAM_EFFECT_COL] = pd.to_numeric(hippo_assoc[PERTURBSEQ_PROGRAM_EFFECT_COL], errors="coerce")

hippo_avg = (
    hippo_assoc
    .groupby(["Sample", "program"], as_index=False)
    .agg(
        hippo_kd_mean_effect=(PERTURBSEQ_PROGRAM_EFFECT_COL, "mean"),
        hippo_kd_median_effect=(PERTURBSEQ_PROGRAM_EFFECT_COL, "median"),
        n_hippo_targets=("target_name", "nunique"),
        hippo_targets=("target_name", lambda x: ";".join(sorted(set(map(str, x))))),
    )
)

print("Hippo KD averaged effect table:")
display(hippo_avg.head())

# %%
concordance_records = []
merged_records = []

for contrast_name, sub in sm_prog.groupby("contrast_name"):
    stage_day = sub["stage_day_for_perturbseq"].iloc[0]
    if stage_day is None or pd.isna(stage_day):
        continue
    kd_sub = hippo_avg.loc[hippo_avg["Sample"] == stage_day].copy()
    merged = sub.merge(kd_sub, on="program", how="inner")
    if merged.empty:
        continue

    x = merged["hippo_kd_mean_effect"].to_numpy(dtype=float)
    y = merged["diff_treated_minus_control"].to_numpy(dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    merged = merged.loc[valid].copy()

    pearson_r, pearson_p = pearsonr(x, y) if len(x) > 2 else (np.nan, np.nan)
    spearman_r, spearman_p = spearmanr(x, y) if len(x) > 2 else (np.nan, np.nan)

    concordance_records.append({
        "contrast_name": contrast_name,
        "treated_condition": sub["treated_condition"].iloc[0],
        "control_condition": sub["control_condition"].iloc[0],
        "stage_day_for_perturbseq": stage_day,
        "small_molecule_effect_col": "diff_treated_minus_control",
        "perturbseq_effect_col": PERTURBSEQ_PROGRAM_EFFECT_COL,
        "hippo_targets_used": ";".join(available_hippo_targets),
        "n_programs": len(x),
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_r": spearman_r,
        "spearman_p": spearman_p,
    })
    merged["contrast_name"] = contrast_name
    merged_records.append(merged)

concordance_df = pd.DataFrame(concordance_records).sort_values("contrast_name")
concordance_df["spearman_padj_bh"] = bh_adjust(concordance_df["spearman_p"].to_numpy())
concordance_df["pearson_padj_bh"] = bh_adjust(concordance_df["pearson_p"].to_numpy())

program_level_merged_df = pd.concat(merged_records, ignore_index=True) if merged_records else pd.DataFrame()

out1 = OUT_DIR / "small_molecule_vs_perturbseq_HippoKD_program_concordance.tsv"
out2 = OUT_DIR / "small_molecule_vs_perturbseq_HippoKD_program_level_values.tsv"
concordance_df.to_csv(out1, sep="\t", index=False)
program_level_merged_df.to_csv(out2, sep="\t", index=False)

print("Saved:", out1)
print("Saved:", out2)
display(concordance_df)

# Key program details
key_program_details = program_level_merged_df.loc[program_level_merged_df["program"].isin(KEY_PROGRAMS)].copy()
key_path = OUT_DIR / "small_molecule_vs_perturbseq_HippoKD_key_programs.tsv"
key_program_details.to_csv(key_path, sep="\t", index=False)
print("Saved:", key_path)
display(key_program_details[["contrast_name", "program", "diff_treated_minus_control", "hippo_kd_mean_effect", "n_hippo_targets", "hippo_targets"]].head(100))

# %% [markdown]
# ## Quick plots: concordance across all programs and key programs

# %%
# Bar plot of Spearman correlations across all programs
fig, ax = plt.subplots(figsize=(8, 4))
plot_df = concordance_df.copy()
plot_df = plot_df.sort_values("spearman_r")
colors = ["#C2410C" if "TDI" in x else "#1D4ED8" for x in plot_df["contrast_name"]]
ax.barh(plot_df["contrast_name"], plot_df["spearman_r"], color=colors)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlabel("Spearman correlation across programs\nsmall molecule effect vs perturb-seq Hippo KD effect")
ax.set_ylabel("Contrast")
ax.set_title("Program-level concordance with perturb-seq Hippo KD")
fig.tight_layout()
fig.savefig(OUT_DIR / "program_level_concordance_HippoKD_spearman.svg")
fig.savefig(OUT_DIR / "program_level_concordance_HippoKD_spearman.pdf")
fig.savefig(OUT_DIR / "program_level_concordance_HippoKD_spearman.png", dpi=300)
plt.show()

# %%
# Scatter plots for selected contrasts
SELECTED_CONTRASTS = [
    "TDI_at_D2_LM",
    "IAG_at_D2_LM",
    "TDI_at_D3_AE",
    "IAG_at_D3_AE",
    "TDI_D2_to_D3_vs_D3_AE",
    "IAG_D2_to_D3_vs_D3_AE",
]

n = len(SELECTED_CONTRASTS)
fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=False, sharey=False)
axes = axes.ravel()

for ax, contrast in zip(axes, SELECTED_CONTRASTS):
    sub = program_level_merged_df.loc[program_level_merged_df["contrast_name"] == contrast].copy()
    if sub.empty:
        ax.axis("off")
        ax.set_title(contrast + "\nmissing")
        continue

    ax.scatter(sub["hippo_kd_mean_effect"], sub["diff_treated_minus_control"], s=18, alpha=0.7)
    for _, row in sub.loc[sub["program"].isin(KEY_PROGRAMS)].iterrows():
        ax.scatter(row["hippo_kd_mean_effect"], row["diff_treated_minus_control"], s=45, edgecolor="black")
        ax.text(row["hippo_kd_mean_effect"], row["diff_treated_minus_control"], row["program"], fontsize=8)
    ax.axhline(0, color="grey", linewidth=0.7)
    ax.axvline(0, color="grey", linewidth=0.7)
    corr = concordance_df.loc[concordance_df["contrast_name"] == contrast, "spearman_r"]
    corr_txt = f"rho={corr.iloc[0]:.2f}" if len(corr) else "rho=NA"
    ax.set_title(f"{contrast}\n{corr_txt}", fontsize=8)
    ax.set_xlabel("Perturb-seq Hippo KD program effect")
    ax.set_ylabel("Small molecule program effect")

fig.tight_layout()
fig.savefig(OUT_DIR / "program_level_concordance_HippoKD_scatter_selected.svg")
fig.savefig(OUT_DIR / "program_level_concordance_HippoKD_scatter_selected.pdf")
fig.savefig(OUT_DIR / "program_level_concordance_HippoKD_scatter_selected.png", dpi=300)
plt.show()

print("All outputs saved to", OUT_DIR)
