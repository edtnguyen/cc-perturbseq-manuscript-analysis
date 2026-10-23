# %% [markdown]
# # Hippo CRISPRi perturb-seq vs small-molecule program-effect comparison plots
#
# This script compares Hippo perturbation effects from perturb-seq/programDE with
# TDI/IAG effects from the Hippo inhibitor scRNA-seq experiment projected onto the
# same perturb-seq cNMF programs.
#
# Main assumptions:
# - perturb-seq effects are from the cNMF workbook, sheet `Perturbation Association 1`
# - small-molecule effects are CLR-transformed matched-control effects when available
# - program labels are zero-based: P0..P49
# - output SVG/PDF text is Illustrator-editable

# %%
from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

warnings.filterwarnings("ignore")

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",  # editable text in Illustrator
})

# %% [markdown]
# ## Paths and user settings

# %%
BASE = Path("/scratch/users/edtng/20260325_Hippo_timecourse")
OUT_DIR = BASE / "hippo_perturbseq_vs_small_molecule_plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CNMF_WORKBOOK_PATH = BASE / "cNMF_50_0_2_032526_6co.xlsx"
PROGRAM_EFFECT_SHEET = "Perturbation Association 1"

# perturb-seq effect scale. If you want the regression coefficient itself, use beta_hat.
PERTURBSEQ_EFFECT_COL = "approx_log2FC"  # or "beta_hat"

# Small molecule CLR effect table generated from the previous CLR-matched notebook.
CLR_DIR = BASE / "perturbseq_cNMF_matched_DE" / "program_score_DE_clr"
SM_CLR_EFFECTS = CLR_DIR / "small_molecule_projected_cNMF_CLR_matched_effects.tsv"
SM_CLR_PAIRED = CLR_DIR / "small_molecule_projected_cNMF_CLR_paired_replicate_differences.tsv"

# Fallback old non-CLR program-score effects.
SM_NORM_DIR = BASE / "perturbseq_cNMF_matched_DE" / "program_score_DE"
SM_NORM_EFFECTS = SM_NORM_DIR / "projected_cNMF_program_score_matched_contrasts.tsv"
SM_NORM_PAIRED = SM_NORM_DIR / "projected_cNMF_program_score_paired_replicate_differences.tsv"

PROGRAM_META_PATH = Path("/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/data/program_activity_long_v4.csv")

# Hippo CRISPRi targets to average. Edit if needed.
HIPPO_KD_TARGETS = ["LATS2", "NF2", "SAV1", "TAOK1", "WWC3", "AMOTL2"]

# Program sets for focused panels.
D2_PROGRAMS = ["P14", "P15", "P17", "P19", "P23"]
D3_EC_PROGRAMS = ["P2", "P3", "P8", "P11", "P42"]
D3_ALT_PROGRAMS = ["P16", "P24", "P27", "P33", "P34", "P35"]
SELECTED_PROGRAMS = D2_PROGRAMS + D3_EC_PROGRAMS + D3_ALT_PROGRAMS
KEY_PROGRAMS = ["P15", "P17", "P23", "P2", "P34", "P24"]

# Small molecule contrasts to display. The script handles both these names and older names like TDI_at_D2_LM.
SM_CONTRAST_ORDER = [
    "D1_MPS+TDI", "D1_MPS+IAG",
    "D1.33_LM+TDI", "D1.33_LM+IAG",
    "D1.67_LM+TDI", "D1.67_LM+IAG",
    "D2_LM+TDI", "D2_LM+IAG",
    "D3_AE+TDI", "D3_AE+IAG",
    "D2_LM+TDI->D3_AE+TDI", "D2_LM+IAG->D3_AE+IAG",
]

FOCUSED_SM_CONTRAST_ORDER = [
    "D2_LM+TDI", "D2_LM+IAG",
    "D3_AE+TDI", "D3_AE+IAG",
    "D2_LM+TDI->D3_AE+TDI", "D2_LM+IAG->D3_AE+IAG",
]

# %% [markdown]
# ## Helper functions

# %%
def program_sort_key(program: str) -> int:
    m = re.search(r"(\d+)$", str(program))
    if m is None:
        return 10_000
    return int(m.group(1))


def coerce_program_name(x) -> str:
    s = str(x).strip()
    s = s.replace("P", "")
    return f"P{int(float(s))}"


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_")


def contrast_to_original_name(contrast: str) -> str:
    """Convert older contrast labels to treated sample names where possible."""
    mapping = {
        "TDI_at_D1_MPS": "D1_MPS+TDI",
        "IAG_at_D1_MPS": "D1_MPS+IAG",
        "TDI_at_D1.33_LM": "D1.33_LM+TDI",
        "IAG_at_D1.33_LM": "D1.33_LM+IAG",
        "TDI_at_D1.67_LM": "D1.67_LM+TDI",
        "IAG_at_D1.67_LM": "D1.67_LM+IAG",
        "TDI_at_D2_LM": "D2_LM+TDI",
        "IAG_at_D2_LM": "D2_LM+IAG",
        "TDI_at_D3_AE": "D3_AE+TDI",
        "IAG_at_D3_AE": "D3_AE+IAG",
        "TDI_D2_to_D3_vs_D3_AE": "D2_LM+TDI->D3_AE+TDI",
        "IAG_D2_to_D3_vs_D3_AE": "D2_LM+IAG->D3_AE+IAG",
    }
    return mapping.get(str(contrast), str(contrast))


def short_contrast_label(name: str) -> str:
    name = str(name)
    name = name.replace("D2_LM+TDI->D3_AE+TDI", "D2→D3 TDI")
    name = name.replace("D2_LM+IAG->D3_AE+IAG", "D2→D3 IAG")
    name = name.replace("D1.33_LM", "D1.33")
    name = name.replace("D1.67_LM", "D1.67")
    name = name.replace("D1_MPS", "D1")
    name = name.replace("D2_LM", "D2")
    name = name.replace("D3_AE", "D3")
    return name


def load_program_labels(program_meta_path: Path):
    labels = {f"P{i}": f"P{i}" for i in range(50)}
    if not program_meta_path.exists():
        return labels
    meta = pd.read_csv(program_meta_path)
    if "program_name" not in meta.columns:
        return labels
    label_col = "program_label" if "program_label" in meta.columns else None
    if label_col is None:
        return labels
    for _, row in meta.drop_duplicates("program_name").iterrows():
        try:
            p = f"P{int(float(row['program_name']))}"
        except Exception:
            continue
        label = str(row[label_col]) if pd.notna(row[label_col]) else p
        # Keep compact labels for figures.
        label = label.split(" — ")[-1]
        label = label.replace("program", "prog.")
        if len(label) > 34:
            label = label[:33] + "…"
        labels[p] = f"{p} {label}"
    return labels


def savefig_all(fig, out_prefix: Path):
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".svg"))
    fig.savefig(out_prefix.with_suffix(".pdf"))
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300)


def center_limits(mat, q=1.0):
    x = np.asarray(mat, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return -1, 1
    if q < 1:
        vmax = np.quantile(np.abs(x), q)
    else:
        vmax = np.max(np.abs(x))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1e-9
    return -vmax, vmax


def plot_heatmap(df, title, out_prefix, cmap="RdBu_r", center=True, annotate=False, ylabels=None):
    mat = df.to_numpy(dtype=float)
    vmin, vmax = center_limits(mat, q=1.0) if center else (np.nanmin(mat), np.nanmax(mat))
    n_rows, n_cols = df.shape
    fig_w = max(7, 0.48 * n_cols + 4)
    fig_h = max(4, 0.28 * n_rows + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    mesh = ax.pcolormesh(np.arange(n_cols + 1), np.arange(n_rows + 1), mat, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_xticklabels([short_contrast_label(x) for x in df.columns], rotation=55, ha="right")
    ax.set_yticklabels(ylabels if ylabels is not None else df.index.tolist())
    ax.set_title(title)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if annotate:
        for i in range(n_rows):
            for j in range(n_cols):
                val = mat[i, j]
                if np.isfinite(val):
                    ax.text(j + 0.5, i + 0.5, f"{val:.2f}", ha="center", va="center", fontsize=5)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label("effect size")
    fig.tight_layout()
    savefig_all(fig, out_prefix)
    plt.close(fig)

# %% [markdown]
# ## Load perturb-seq and small-molecule effect tables

# %%
program_labels = load_program_labels(PROGRAM_META_PATH)

if not CNMF_WORKBOOK_PATH.exists():
    raise FileNotFoundError(f"Missing cNMF workbook: {CNMF_WORKBOOK_PATH}")

assoc = pd.read_excel(CNMF_WORKBOOK_PATH, sheet_name=PROGRAM_EFFECT_SHEET)
required = {"target_name", "program_name", "Sample", PERTURBSEQ_EFFECT_COL}
missing = required - set(assoc.columns)
if missing:
    raise ValueError(f"cNMF workbook missing columns: {missing}. Columns: {list(assoc.columns)}")

assoc = assoc.copy()
assoc["program"] = assoc["program_name"].map(coerce_program_name)
assoc["Sample"] = assoc["Sample"].astype(str)
assoc["target_name"] = assoc["target_name"].astype(str)
assoc[PERTURBSEQ_EFFECT_COL] = pd.to_numeric(assoc[PERTURBSEQ_EFFECT_COL], errors="coerce")

hippo_indiv = assoc.loc[
    assoc["target_name"].isin(HIPPO_KD_TARGETS) & assoc["Sample"].isin(["D2", "D3"]),
    ["target_name", "Sample", "program", PERTURBSEQ_EFFECT_COL]
].dropna().copy()

hippo_mean = (
    hippo_indiv
    .groupby(["Sample", "program"], observed=True)[PERTURBSEQ_EFFECT_COL]
    .mean()
    .reset_index()
    .rename(columns={PERTURBSEQ_EFFECT_COL: "hippo_kd_mean_effect"})
)

print("Perturb-seq Hippo individual effects:", hippo_indiv.shape)
print("Perturb-seq Hippo mean effects:", hippo_mean.shape)

# Small molecule effect table. Prefer CLR.
if SM_CLR_EFFECTS.exists():
    sm_effect = pd.read_csv(SM_CLR_EFFECTS, sep="\t")
    sm_effect_source = "CLR"
    sm_effect_col = "mean_diff_clr"
    sm_paired_path = SM_CLR_PAIRED
    sm_paired_effect_col = "paired_diff_clr"
elif SM_NORM_EFFECTS.exists():
    sm_effect = pd.read_csv(SM_NORM_EFFECTS, sep="\t")
    sm_effect_source = "normalized_usage"
    sm_effect_col = "mean_diff" if "mean_diff" in sm_effect.columns else "mean_diff_usage"
    sm_paired_path = SM_NORM_PAIRED
    sm_paired_effect_col = "paired_diff"
else:
    raise FileNotFoundError("Could not find small-molecule effect table. Run the program-score effect analysis first.")

sm_effect = sm_effect.copy()
sm_effect["program"] = sm_effect["program"].map(coerce_program_name)
sm_effect["contrast_original"] = sm_effect["contrast"].map(contrast_to_original_name)
sm_effect[sm_effect_col] = pd.to_numeric(sm_effect[sm_effect_col], errors="coerce")
print("Small molecule effect source:", sm_effect_source)
print("Small molecule effects:", sm_effect.shape)

if sm_paired_path.exists():
    sm_paired = pd.read_csv(sm_paired_path, sep="\t")
    sm_paired = sm_paired.copy()
    sm_paired["program"] = sm_paired["program"].map(coerce_program_name)
    sm_paired["contrast_original"] = sm_paired["contrast"].map(contrast_to_original_name)
    if sm_paired_effect_col not in sm_paired.columns:
        # fallback column names
        candidates = [c for c in sm_paired.columns if "paired" in c and "diff" in c]
        if candidates:
            sm_paired_effect_col = candidates[0]
    sm_paired[sm_paired_effect_col] = pd.to_numeric(sm_paired[sm_paired_effect_col], errors="coerce")
    print("Small molecule paired effects:", sm_paired.shape)
else:
    sm_paired = None
    print("No paired replicate table found:", sm_paired_path)

# Save preprocessed tables
hippo_indiv.to_csv(OUT_DIR / "source_perturbseq_HippoKD_individual_program_effects.tsv", sep="\t", index=False)
hippo_mean.to_csv(OUT_DIR / "source_perturbseq_HippoKD_mean_program_effects.tsv", sep="\t", index=False)
sm_effect.to_csv(OUT_DIR / f"source_small_molecule_{sm_effect_source}_program_effects.tsv", sep="\t", index=False)
if sm_paired is not None:
    sm_paired.to_csv(OUT_DIR / f"source_small_molecule_{sm_effect_source}_paired_program_effects.tsv", sep="\t", index=False)

# %% [markdown]
# ## Panel 1: selected program effect heatmap, Hippo CRISPRi vs small molecule

# %%
# Build effect matrix with rows = selected programs and columns = PerturbSeq D2/D3 + small molecule contrasts.
cols = ["Perturb-seq Hippo KD D2", "Perturb-seq Hippo KD D3"] + FOCUSED_SM_CONTRAST_ORDER
mat_records = []
for program in SELECTED_PROGRAMS:
    row = {"program": program}
    for day in ["D2", "D3"]:
        vals = hippo_mean.loc[(hippo_mean["Sample"] == day) & (hippo_mean["program"] == program), "hippo_kd_mean_effect"]
        row[f"Perturb-seq Hippo KD {day}"] = float(vals.iloc[0]) if len(vals) else np.nan
    for contrast in FOCUSED_SM_CONTRAST_ORDER:
        vals = sm_effect.loc[(sm_effect["contrast_original"] == contrast) & (sm_effect["program"] == program), sm_effect_col]
        row[contrast] = float(vals.iloc[0]) if len(vals) else np.nan
    mat_records.append(row)

effect_mat = pd.DataFrame(mat_records).set_index("program")[cols]
effect_mat.to_csv(OUT_DIR / "Panel1_selected_program_effect_matrix_HippoKD_vs_small_molecule.tsv", sep="\t")
plot_heatmap(
    effect_mat,
    title=f"Selected program effects: Hippo CRISPRi vs small molecule ({sm_effect_source})",
    out_prefix=OUT_DIR / "Panel1_selected_program_effect_heatmap_HippoKD_vs_small_molecule",
    cmap="RdBu_r",
    center=True,
    annotate=True,
    ylabels=[program_labels.get(p, p) for p in effect_mat.index],
)

# %% [markdown]
# ## Panel 2: replicate/target distribution plot for key programs

# %%
def plot_key_program_effect_distributions(programs=KEY_PROGRAMS, contrasts=FOCUSED_SM_CONTRAST_ORDER):
    n_rows = len(programs)
    categories = ["HippoKD D2", "HippoKD D3"] + [short_contrast_label(c) for c in contrasts]
    x = np.arange(len(categories))
    fig_w = max(11, 0.8 * len(categories) + 3)
    fig_h = max(2.4 * n_rows, 5)
    fig, axes = plt.subplots(n_rows, 1, figsize=(fig_w, fig_h), sharex=True)
    if n_rows == 1:
        axes = [axes]
    rng = np.random.default_rng(0)

    for ax, program in zip(axes, programs):
        # collect values per category
        all_vals = []
        # Hippo D2/D3 individual target effects
        for day in ["D2", "D3"]:
            vals = hippo_indiv.loc[(hippo_indiv["Sample"] == day) & (hippo_indiv["program"] == program), PERTURBSEQ_EFFECT_COL].to_numpy(float)
            all_vals.append(vals)
        # Small molecule paired replicate effects
        for contrast in contrasts:
            if sm_paired is not None:
                vals = sm_paired.loc[(sm_paired["contrast_original"] == contrast) & (sm_paired["program"] == program), sm_paired_effect_col].to_numpy(float)
            else:
                vals = sm_effect.loc[(sm_effect["contrast_original"] == contrast) & (sm_effect["program"] == program), sm_effect_col].to_numpy(float)
            all_vals.append(vals)

        for i, vals in enumerate(all_vals):
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            jitter = rng.normal(0, 0.055, size=vals.size)
            ax.scatter(np.full(vals.size, x[i]) + jitter, vals, s=16, alpha=0.85, linewidth=0.3, edgecolors="black")
            ax.plot([x[i] - 0.18, x[i] + 0.18], [np.mean(vals), np.mean(vals)], linewidth=2)
        ax.axhline(0, linewidth=0.8)
        ax.set_ylabel(program)
        ax.set_title(program_labels.get(program, program), loc="left", fontsize=8)
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(categories, rotation=55, ha="right")
    axes[-1].set_xlabel("Effect source / contrast")
    fig.suptitle("Program effects: Hippo KD targets and small-molecule matched replicates", y=0.995)
    fig.tight_layout()
    savefig_all(fig, OUT_DIR / "Panel2_key_program_effect_distributions_targets_and_replicates")
    plt.close(fig)

plot_key_program_effect_distributions()

# %% [markdown]
# ## Panel 3: small molecule vs Hippo CRISPRi scatter plots

# %%
def compute_concordance(sm_contrasts=SM_CONTRAST_ORDER, programs=None):
    if programs is None:
        programs = sorted(sm_effect["program"].unique(), key=program_sort_key)
    records = []
    merged_records = []
    for contrast in sm_contrasts:
        sm_sub = sm_effect.loc[sm_effect["contrast_original"] == contrast, ["program", sm_effect_col]].dropna().copy()
        sm_sub = sm_sub[sm_sub["program"].isin(programs)]
        for ref_day in ["D2", "D3"]:
            ref_sub = hippo_mean.loc[hippo_mean["Sample"] == ref_day, ["program", "hippo_kd_mean_effect"]].dropna().copy()
            ref_sub = ref_sub[ref_sub["program"].isin(programs)]
            merged = ref_sub.merge(sm_sub, on="program", how="inner")
            if merged.shape[0] < 3:
                continue
            x = merged["hippo_kd_mean_effect"].to_numpy(float)
            y = merged[sm_effect_col].to_numpy(float)
            pearson_r, pearson_p = pearsonr(x, y)
            spearman_r, spearman_p = spearmanr(x, y)
            records.append({
                "small_molecule_contrast": contrast,
                "perturbseq_reference_day": ref_day,
                "n_programs": merged.shape[0],
                "pearson_r": pearson_r,
                "pearson_pvalue": pearson_p,
                "spearman_r": spearman_r,
                "spearman_pvalue": spearman_p,
                "program_set": "selected" if set(programs) == set(SELECTED_PROGRAMS) else "all_or_custom",
            })
            tmp = merged.copy()
            tmp["small_molecule_contrast"] = contrast
            tmp["perturbseq_reference_day"] = ref_day
            merged_records.append(tmp)
    return pd.DataFrame(records), pd.concat(merged_records, ignore_index=True) if merged_records else pd.DataFrame()

conc_all, conc_values_all = compute_concordance(SM_CONTRAST_ORDER, programs=sorted(sm_effect["program"].unique(), key=program_sort_key))
conc_sel, conc_values_sel = compute_concordance(FOCUSED_SM_CONTRAST_ORDER, programs=SELECTED_PROGRAMS)
conc_all.to_csv(OUT_DIR / "Panel3_concordance_all_programs.tsv", sep="\t", index=False)
conc_values_all.to_csv(OUT_DIR / "Panel3_concordance_all_programs_values.tsv", sep="\t", index=False)
conc_sel.to_csv(OUT_DIR / "Panel3_concordance_selected_programs.tsv", sep="\t", index=False)
conc_values_sel.to_csv(OUT_DIR / "Panel3_concordance_selected_programs_values.tsv", sep="\t", index=False)

# Heatmap of Spearman rho across all programs.
conc_heat = conc_all.pivot(index="small_molecule_contrast", columns="perturbseq_reference_day", values="spearman_r")
ordered_rows = [c for c in SM_CONTRAST_ORDER if c in conc_heat.index]
conc_heat = conc_heat.loc[ordered_rows, [c for c in ["D2", "D3"] if c in conc_heat.columns]]
plot_heatmap(
    conc_heat,
    title="Small molecule vs perturb-seq Hippo KD concordance\nSpearman r across all programs",
    out_prefix=OUT_DIR / "Panel3_concordance_heatmap_all_programs",
    cmap="RdBu_r",
    center=True,
    annotate=True,
)

# %%
def plot_scatter_grid(conc_values, contrasts=FOCUSED_SM_CONTRAST_ORDER, ref_days=("D2", "D3"), label_programs=KEY_PROGRAMS, suffix="selected_programs"):
    n_rows = len(contrasts)
    n_cols = len(ref_days)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.1 * n_cols, 2.9 * n_rows), squeeze=False)
    for i, contrast in enumerate(contrasts):
        for j, ref_day in enumerate(ref_days):
            ax = axes[i, j]
            sub = conc_values.loc[(conc_values["small_molecule_contrast"] == contrast) & (conc_values["perturbseq_reference_day"] == ref_day)].copy()
            if sub.empty:
                ax.axis("off")
                continue
            x = sub["hippo_kd_mean_effect"].to_numpy(float)
            y = sub[sm_effect_col].to_numpy(float)
            ax.scatter(x, y, s=22, alpha=0.75, linewidth=0.3, edgecolors="black")
            ax.axhline(0, linewidth=0.8)
            ax.axvline(0, linewidth=0.8)
            for _, row in sub.iterrows():
                if row["program"] in label_programs:
                    ax.text(row["hippo_kd_mean_effect"], row[sm_effect_col], row["program"], fontsize=7)
            if len(x) >= 3:
                rho, _ = spearmanr(x, y)
                pr, _ = pearsonr(x, y)
                ax.text(0.03, 0.97, f"rho={rho:.2f}\nr={pr:.2f}\nn={len(x)}", transform=ax.transAxes, ha="left", va="top", fontsize=7)
            ax.set_title(f"{short_contrast_label(contrast)} vs HippoKD {ref_day}", fontsize=8)
            if i == n_rows - 1:
                ax.set_xlabel(f"Perturb-seq Hippo KD effect ({ref_day})")
            if j == 0:
                ax.set_ylabel(f"Small molecule effect\n({sm_effect_source})")
    fig.tight_layout()
    savefig_all(fig, OUT_DIR / f"Panel3_scatter_grid_{suffix}")
    plt.close(fig)

plot_scatter_grid(conc_values_sel, suffix="selected_programs")
plot_scatter_grid(conc_values_all, contrasts=FOCUSED_SM_CONTRAST_ORDER, suffix="all_programs_key_labels")

# %% [markdown]
# ## Panel 4: lollipop comparison for selected programs against D2/D3 Hippo KD

# %%
def plot_selected_program_lollipop(programs=KEY_PROGRAMS, contrasts=FOCUSED_SM_CONTRAST_ORDER):
    # Build a long table with rows = effect sources.
    rows = []
    for p in programs:
        for day in ["D2", "D3"]:
            vals = hippo_mean.loc[(hippo_mean["Sample"] == day) & (hippo_mean["program"] == p), "hippo_kd_mean_effect"]
            if len(vals):
                rows.append({"program": p, "source": f"HippoKD {day}", "effect": float(vals.iloc[0]), "group": "Perturb-seq"})
        for c in contrasts:
            vals = sm_effect.loc[(sm_effect["contrast_original"] == c) & (sm_effect["program"] == p), sm_effect_col]
            if len(vals):
                rows.append({"program": p, "source": short_contrast_label(c), "effect": float(vals.iloc[0]), "group": "Small molecule"})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "Panel4_selected_program_lollipop_source.tsv", sep="\t", index=False)
    sources = ["HippoKD D2", "HippoKD D3"] + [short_contrast_label(c) for c in contrasts]
    fig, axes = plt.subplots(1, len(programs), figsize=(3.4 * len(programs), 4.2), sharey=True)
    if len(programs) == 1:
        axes = [axes]
    for ax, p in zip(axes, programs):
        sub = df[df["program"] == p].copy()
        sub["source"] = pd.Categorical(sub["source"], categories=sources, ordered=True)
        sub = sub.sort_values("source")
        y = np.arange(sub.shape[0])
        x = sub["effect"].to_numpy(float)
        ax.axvline(0, linewidth=0.8)
        for yi, xi in zip(y, x):
            ax.plot([0, xi], [yi, yi], linewidth=1.5)
        ax.scatter(x, y, s=28, linewidth=0.4, edgecolors="black", zorder=3)
        ax.set_yticks(y)
        if ax is axes[0]:
            ax.set_yticklabels(sub["source"].tolist())
        else:
            ax.set_yticklabels([])
        ax.set_title(program_labels.get(p, p), fontsize=8)
        ax.set_xlabel("effect")
    fig.suptitle(f"Selected program effects: perturb-seq HippoKD and small molecules ({sm_effect_source})", y=0.99)
    fig.tight_layout()
    savefig_all(fig, OUT_DIR / "Panel4_selected_program_effect_lollipop")
    plt.close(fig)

plot_selected_program_lollipop()

# %% [markdown]
# ## Done

# %%
print("Saved outputs to:", OUT_DIR)
for p in sorted(OUT_DIR.glob("Panel*.svg")):
    print(p.name)
