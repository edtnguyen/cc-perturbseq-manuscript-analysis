# %% [markdown]
# # Program-level concordance against perturb-seq D2 and D3 Hippo KD effects
#
# This notebook compares every small-molecule program-effect vector against **both**
# perturb-seq D2 and perturb-seq D3 Hippo KD program-effect vectors.
#
# It keeps the original small-molecule sample/treatment names, e.g.
# `D2_LM+TDI`, `D3_AE+IAG`, `D2_LM+TDI->D3_AE+TDI`, rather than renaming them
# to `TDI_at_D2_LM`.
#
# Each small-molecule effect vector is still:
#
# ```text
# treated_condition - matched control_condition
# ```
#
# For example:
#
# ```text
# D2_LM+TDI means D2_LM+TDI - D2_LM
# D3_AE+IAG means D3_AE+IAG - D3_AE
# D2_LM+TDI->D3_AE+TDI means D2_LM+TDI->D3_AE+TDI - D3_AE
# ```

# %%
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib as mpl
import matplotlib.pyplot as plt
try:
    from IPython.display import display
except Exception:
    def display(x):
        print(x)

mpl.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.family": "Arial",
})

# %% [markdown]
# ## Paths and settings

# %%
LOCAL_DIR = Path("/mnt/data")
SHERLOCK_BASE = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE")

CNMF_XLSX = LOCAL_DIR / "cNMF_50_0_2_032526_6co.xlsx"
if not CNMF_XLSX.exists():
    CNMF_XLSX = Path("/scratch/users/edtng/20260325_Hippo_timecourse/cNMF_50_0_2_032526_6co.xlsx")

PROGRAM_SCORE_CONTRASTS = LOCAL_DIR / "projected_cNMF_program_score_matched_contrasts.tsv"
if not PROGRAM_SCORE_CONTRASTS.exists():
    PROGRAM_SCORE_CONTRASTS = Path(
        "/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE/program_score_DE/"
        "projected_cNMF_program_score_matched_contrasts.tsv"
    )

OUT_DIR = LOCAL_DIR / "hippo_concordance_D2_D3_original_names"
if not LOCAL_DIR.exists():
    OUT_DIR = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE/hippo_concordance_D2_D3_original_names")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Hippo KD targets to average in perturb-seq association table. Edit if desired.
HIPPO_KD_TARGETS = ["LATS2", "NF2", "SAV1", "TAOK1", "WWC3", "AMOTL2"]

# Perturb-seq effect column to use.
PERTURBSEQ_EFFECT_COL = "approx_log2FC"  # alternative: "beta_hat"
SMALL_MOLECULE_EFFECT_COL = "diff_treated_minus_control"

REFERENCE_DAYS = ["D2", "D3"]
KEY_PROGRAMS = ["P15", "P17", "P34"]

# Original small-molecule sample/treatment order.
SMALL_MOLECULE_ORDER = [
    "D1_MPS+TDI",
    "D1_MPS+IAG",
    "D1.33_LM+TDI",
    "D1.33_LM+IAG",
    "D1.67_LM+TDI",
    "D1.67_LM+IAG",
    "D2_LM+TDI",
    "D2_LM+IAG",
    "D3_AE+TDI",
    "D3_AE+IAG",
    "D2_LM+TDI->D3_AE+TDI",
    "D2_LM+IAG->D3_AE+IAG",
]

print("CNMF_XLSX:", CNMF_XLSX, CNMF_XLSX.exists())
print("PROGRAM_SCORE_CONTRASTS:", PROGRAM_SCORE_CONTRASTS, PROGRAM_SCORE_CONTRASTS.exists())
print("OUT_DIR:", OUT_DIR)

# %% [markdown]
# ## Helper functions

# %%
def normalize_program_name(x):
    if isinstance(x, str) and x.startswith("P"):
        return x
    return f"P{int(float(x))}"


def bh_adjust(pvals):
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
    tmp = np.empty_like(q)
    tmp[order] = q
    out[valid] = tmp
    return out


def treatment_class(x):
    x = str(x)
    if "+TDI" in x:
        return "TDI"
    if "+IAG" in x:
        return "IAG"
    return "Other"


def plot_heatmap(df, title, out_prefix, cbar_label="Spearman rho"):
    mat = df.to_numpy(dtype=float)
    vmax = np.nanmax(np.abs(mat))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1
    vmin = -vmax

    fig, ax = plt.subplots(figsize=(5.2, max(4.0, 0.42 * df.shape[0])))
    mesh = ax.pcolormesh(
        np.arange(df.shape[1] + 1),
        np.arange(df.shape[0] + 1),
        mat,
        cmap="RdBu_r",
        vmin=vmin,
        vmax=vmax,
        shading="flat",
        linewidth=0.4,
        edgecolors="white",
    )
    ax.set_xlim(0, df.shape[1])
    ax.set_ylim(df.shape[0], 0)
    ax.set_xticks(np.arange(df.shape[1]) + 0.5)
    ax.set_yticks(np.arange(df.shape[0]) + 0.5)
    ax.set_xticklabels(df.columns.tolist(), rotation=0)
    ax.set_yticklabels(df.index.tolist())
    ax.set_title(title)

    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            val = mat[i, j]
            if np.isfinite(val):
                ax.text(j + 0.5, i + 0.5, f"{val:.2f}", ha="center", va="center", fontsize=7)

    cbar = fig.colorbar(mesh, ax=ax, pad=0.03)
    cbar.set_label(cbar_label)
    fig.tight_layout()
    for ext in ["svg", "pdf", "png"]:
        fig.savefig(OUT_DIR / f"{out_prefix}.{ext}", bbox_inches="tight")
    plt.show()
    plt.close(fig)

# %% [markdown]
# ## Load cNMF workbook and small molecule program-effect table

# %%
assoc = pd.read_excel(CNMF_XLSX, sheet_name="Perturbation Association 1")
assoc = assoc.copy()
assoc["program"] = assoc["program_name"].map(normalize_program_name)
assoc["target_name"] = assoc["target_name"].astype(str)
assoc["Sample"] = assoc["Sample"].astype(str)
assoc[PERTURBSEQ_EFFECT_COL] = pd.to_numeric(assoc[PERTURBSEQ_EFFECT_COL], errors="coerce")

print("Association table shape:", assoc.shape)
print("Samples in perturb-seq association table:", sorted(assoc["Sample"].unique()))
print("Effect column:", PERTURBSEQ_EFFECT_COL)

all_targets = sorted(assoc["target_name"].dropna().unique().tolist())
available_hippo_targets = [t for t in HIPPO_KD_TARGETS if t in set(all_targets)]
missing_hippo_targets = [t for t in HIPPO_KD_TARGETS if t not in set(all_targets)]
print("Available Hippo KD targets:", available_hippo_targets)
print("Missing Hippo KD targets:", missing_hippo_targets)
if not available_hippo_targets:
    raise ValueError("No HIPPO_KD_TARGETS were found in association table.")

sm = pd.read_csv(PROGRAM_SCORE_CONTRASTS, sep="\t")
sm["program"] = sm["program"].astype(str)
sm["treated_condition"] = sm["treated_condition"].astype(str)
sm["control_condition"] = sm["control_condition"].astype(str)
sm[SMALL_MOLECULE_EFFECT_COL] = pd.to_numeric(sm[SMALL_MOLECULE_EFFECT_COL], errors="coerce")

print("Small molecule program score contrast table:", sm.shape)
display(sm[["contrast_name", "treated_condition", "control_condition"]].drop_duplicates())

# %% [markdown]
# ## Average perturb-seq Hippo KD effects at D2 and D3

# %%
hippo_assoc = assoc.loc[
    assoc["target_name"].isin(available_hippo_targets)
    & assoc["Sample"].isin(REFERENCE_DAYS)
].copy()

hippo_avg = (
    hippo_assoc
    .groupby(["Sample", "program"], as_index=False)
    .agg(
        hippo_kd_mean_effect=(PERTURBSEQ_EFFECT_COL, "mean"),
        hippo_kd_median_effect=(PERTURBSEQ_EFFECT_COL, "median"),
        n_hippo_targets=("target_name", "nunique"),
        hippo_targets=("target_name", lambda x: ";".join(sorted(set(map(str, x))))),
    )
)

print("hippo_avg shape:", hippo_avg.shape)
display(hippo_avg.head())

hippo_avg.to_csv(OUT_DIR / "perturbseq_HippoKD_average_program_effects_D2_D3.tsv", sep="\t", index=False)

# %% [markdown]
# ## Compare every small-molecule effect vector to both perturb-seq D2 and D3

# %%
concordance_records = []
merged_records = []

for (contrast_name, treated_condition, control_condition), sub in sm.groupby(["contrast_name", "treated_condition", "control_condition"], sort=False):
    for ref_day in REFERENCE_DAYS:
        kd_sub = hippo_avg.loc[hippo_avg["Sample"] == ref_day].copy()
        merged = sub.merge(kd_sub, on="program", how="inner")
        if merged.empty:
            continue

        x = merged["hippo_kd_mean_effect"].to_numpy(dtype=float)
        y = merged[SMALL_MOLECULE_EFFECT_COL].to_numpy(dtype=float)
        valid = np.isfinite(x) & np.isfinite(y)
        merged = merged.loc[valid].copy()
        x = x[valid]
        y = y[valid]

        pearson_r, pearson_p = pearsonr(x, y) if len(x) > 2 else (np.nan, np.nan)
        spearman_r, spearman_p = spearmanr(x, y) if len(x) > 2 else (np.nan, np.nan)

        concordance_records.append({
            "small_molecule_sample": treated_condition,
            "contrast_name_original": contrast_name,
            "control_condition": control_condition,
            "reference_perturbseq_day": ref_day,
            "small_molecule_effect_col": SMALL_MOLECULE_EFFECT_COL,
            "perturbseq_effect_col": PERTURBSEQ_EFFECT_COL,
            "hippo_targets_used": ";".join(available_hippo_targets),
            "n_programs": len(x),
            "pearson_r": pearson_r,
            "pearson_p": pearson_p,
            "spearman_r": spearman_r,
            "spearman_p": spearman_p,
        })
        merged["small_molecule_sample"] = treated_condition
        merged["contrast_name_original"] = contrast_name
        merged["control_condition"] = control_condition
        merged["reference_perturbseq_day"] = ref_day
        merged_records.append(merged)

concordance_df = pd.DataFrame(concordance_records)
concordance_df["spearman_padj_bh"] = bh_adjust(concordance_df["spearman_p"].to_numpy())
concordance_df["pearson_padj_bh"] = bh_adjust(concordance_df["pearson_p"].to_numpy())

# Preserve requested sample order.
concordance_df["small_molecule_sample"] = pd.Categorical(
    concordance_df["small_molecule_sample"],
    categories=SMALL_MOLECULE_ORDER,
    ordered=True,
)
concordance_df["reference_perturbseq_day"] = pd.Categorical(
    concordance_df["reference_perturbseq_day"],
    categories=REFERENCE_DAYS,
    ordered=True,
)
concordance_df = concordance_df.sort_values(["small_molecule_sample", "reference_perturbseq_day"]).reset_index(drop=True)
concordance_df["small_molecule_sample"] = concordance_df["small_molecule_sample"].astype(str)
concordance_df["reference_perturbseq_day"] = concordance_df["reference_perturbseq_day"].astype(str)

program_level_df = pd.concat(merged_records, ignore_index=True) if merged_records else pd.DataFrame()

concordance_path = OUT_DIR / "small_molecule_vs_perturbseq_HippoKD_concordance_all_samples_against_D2_D3.tsv"
program_level_path = OUT_DIR / "small_molecule_vs_perturbseq_HippoKD_program_values_all_samples_against_D2_D3.tsv"
concordance_df.to_csv(concordance_path, sep="\t", index=False)
program_level_df.to_csv(program_level_path, sep="\t", index=False)

print("Saved:", concordance_path)
print("Saved:", program_level_path)
display(concordance_df)

# %% [markdown]
# ## Key programs P15/P17/P34 against both D2 and D3 references

# %%
key_program_df = program_level_df.loc[program_level_df["program"].isin(KEY_PROGRAMS)].copy()
key_program_df = key_program_df.sort_values(["small_molecule_sample", "reference_perturbseq_day", "program"])
key_path = OUT_DIR / "small_molecule_vs_perturbseq_HippoKD_key_programs_all_samples_against_D2_D3.tsv"
key_program_df.to_csv(key_path, sep="\t", index=False)
print("Saved:", key_path)
display(key_program_df[[
    "small_molecule_sample",
    "control_condition",
    "reference_perturbseq_day",
    "program",
    SMALL_MOLECULE_EFFECT_COL,
    "hippo_kd_mean_effect",
    "n_hippo_targets",
    "hippo_targets",
]].head(200))

# %% [markdown]
# ## Plot 1: Spearman heatmap, all small-molecule samples against perturb-seq D2 and D3

# %%
heatmap_spearman = concordance_df.pivot(
    index="small_molecule_sample",
    columns="reference_perturbseq_day",
    values="spearman_r",
)
heatmap_spearman = heatmap_spearman.reindex(index=SMALL_MOLECULE_ORDER, columns=REFERENCE_DAYS)
plot_heatmap(
    heatmap_spearman,
    title="Program-level concordance with perturb-seq Hippo KD\nall small-molecule samples compared to D2 and D3 references",
    out_prefix="program_level_concordance_HippoKD_D2_D3_reference_heatmap_original_names",
)

# %% [markdown]
# ## Plot 2: Grouped bar plot of Spearman correlations

# %%
plot_df = concordance_df.copy()
plot_df["small_molecule_sample"] = pd.Categorical(plot_df["small_molecule_sample"], categories=SMALL_MOLECULE_ORDER, ordered=True)
plot_df = plot_df.sort_values(["small_molecule_sample", "reference_perturbseq_day"])

x = np.arange(len(SMALL_MOLECULE_ORDER))
width = 0.36
fig, ax = plt.subplots(figsize=(12, 4.5))

for i, ref_day in enumerate(REFERENCE_DAYS):
    sub = plot_df.loc[plot_df["reference_perturbseq_day"] == ref_day].set_index("small_molecule_sample")
    vals = sub.reindex(SMALL_MOLECULE_ORDER)["spearman_r"].to_numpy(dtype=float)
    offset = (i - 0.5) * width
    ax.bar(x + offset, vals, width=width, label=f"Perturb-seq {ref_day}")

ax.axhline(0, color="black", linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels(SMALL_MOLECULE_ORDER, rotation=45, ha="right")
ax.set_ylabel("Spearman correlation across 50 programs")
ax.set_title("Small-molecule program effects vs perturb-seq Hippo KD effects\nD2 and D3 references shown separately")
ax.legend(frameon=False)
fig.tight_layout()
for ext in ["svg", "pdf", "png"]:
    fig.savefig(OUT_DIR / f"program_level_concordance_HippoKD_D2_D3_reference_bars_original_names.{ext}", bbox_inches="tight")
plt.show()
plt.close(fig)

# %% [markdown]
# ## Plot 3: Scatter panels for selected small-molecule samples, against both D2 and D3 references

# %%
SELECTED_SAMPLES = [
    "D2_LM+TDI",
    "D2_LM+IAG",
    "D3_AE+TDI",
    "D3_AE+IAG",
    "D2_LM+TDI->D3_AE+TDI",
    "D2_LM+IAG->D3_AE+IAG",
]

fig, axes = plt.subplots(len(SELECTED_SAMPLES), len(REFERENCE_DAYS), figsize=(8.5, 2.3 * len(SELECTED_SAMPLES)), sharex=False, sharey=False)
if len(SELECTED_SAMPLES) == 1:
    axes = np.asarray([axes])

for row_i, sample in enumerate(SELECTED_SAMPLES):
    for col_i, ref_day in enumerate(REFERENCE_DAYS):
        ax = axes[row_i, col_i]
        sub = program_level_df.loc[
            (program_level_df["small_molecule_sample"] == sample)
            & (program_level_df["reference_perturbseq_day"] == ref_day)
        ].copy()
        if sub.empty:
            ax.axis("off")
            ax.set_title(f"{sample}\nvs perturb-seq {ref_day}: missing")
            continue
        ax.scatter(sub["hippo_kd_mean_effect"], sub[SMALL_MOLECULE_EFFECT_COL], s=18, alpha=0.75)
        for _, r in sub.loc[sub["program"].isin(KEY_PROGRAMS)].iterrows():
            ax.scatter(r["hippo_kd_mean_effect"], r[SMALL_MOLECULE_EFFECT_COL], s=52, edgecolor="black")
            ax.text(r["hippo_kd_mean_effect"], r[SMALL_MOLECULE_EFFECT_COL], r["program"], fontsize=8)
        ax.axhline(0, color="grey", linewidth=0.7)
        ax.axvline(0, color="grey", linewidth=0.7)
        rho = concordance_df.loc[
            (concordance_df["small_molecule_sample"] == sample)
            & (concordance_df["reference_perturbseq_day"] == ref_day),
            "spearman_r",
        ]
        rho_txt = f"rho={rho.iloc[0]:.2f}" if len(rho) else "rho=NA"
        ax.set_title(f"{sample}\nvs perturb-seq {ref_day}; {rho_txt}", fontsize=8)
        ax.set_xlabel(f"Perturb-seq {ref_day} Hippo KD program effect")
        ax.set_ylabel("Small-molecule effect")

fig.tight_layout()
for ext in ["svg", "pdf", "png"]:
    fig.savefig(OUT_DIR / f"program_level_concordance_scatter_selected_D2_D3_reference_original_names.{ext}", bbox_inches="tight")
plt.show()
plt.close(fig)

print("All outputs saved to:", OUT_DIR)
