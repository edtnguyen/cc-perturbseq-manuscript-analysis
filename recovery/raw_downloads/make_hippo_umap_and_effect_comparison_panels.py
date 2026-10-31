# %% [markdown]
# # Hippo UMAP + effect-comparison panels
#
# This notebook/script makes figure candidates that are less heatmap-heavy:
#
# 1. UMAP feature panels for projected perturb-seq cNMF programs and representative marker genes/modules.
# 2. Dot/lollipop panels comparing small-molecule effects to perturb-seq Hippo CRISPRi effects on selected programs.
#
# The cNMF workbook path is set to the Sherlock location provided by Tri:
# `/scratch/users/edtng/20260325_Hippo_timecourse/cNMF_50_0_2_032526_6co.xlsx`.

# %%
from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import ttest_rel
import matplotlib as mpl
import matplotlib.pyplot as plt

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
    "svg.fonttype": "none",  # keep text editable in Illustrator
})

# %% [markdown]
# ## Paths and user settings

# %%
BASE = Path("/scratch/users/edtng/20260325_Hippo_timecourse")
PROJ_DIR = BASE / "perturbseq_cNMF_projection"
DE_DIR = BASE / "perturbseq_cNMF_matched_DE"

H5AD_PATH = PROJ_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad"
CNMF_WORKBOOK_PATH = BASE / "cNMF_50_0_2_032526_6co.xlsx"

# Preferred small-molecule effect file: CLR scale to match programDE.
SMALL_MOL_CLR_EFFECT_PATH = DE_DIR / "program_score_DE_clr" / "small_molecule_projected_cNMF_CLR_matched_effects.tsv"
SMALL_MOL_CLR_PAIRED_PATH = DE_DIR / "program_score_DE_clr" / "small_molecule_projected_cNMF_CLR_paired_replicate_differences.tsv"

# Fallback if CLR table does not exist.
SMALL_MOL_SCORE_EFFECT_PATH = DE_DIR / "program_score_DE" / "projected_cNMF_program_score_matched_contrasts.tsv"
SMALL_MOL_SCORE_PAIRED_PATH = DE_DIR / "program_score_DE" / "projected_cNMF_program_score_paired_replicate_differences.tsv"

OUT_DIR = BASE / "hippo_umap_and_effect_comparison_panels"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Relevant programs. Adjust if Kyle/Hunter annotation changes.
D2_PROGRAMS = ["P14", "P15", "P17", "P19", "P23"]
D3_EC_PROGRAMS = ["P2", "P3", "P8", "P11", "P42"]
D3_ALT_PROGRAMS = ["P16", "P24", "P27", "P33", "P34", "P35"]
SELECTED_PROGRAMS = D2_PROGRAMS + D3_EC_PROGRAMS + D3_ALT_PROGRAMS
KEY_PROGRAMS = ["P15", "P17", "P23", "P2", "P34", "P24"]

HIPPO_KD_TARGETS = ["LATS2", "NF2", "SAV1", "TAOK1", "WWC3", "AMOTL2"]
PROGRAM_ASSOC_SHEET = "Perturbation Association 1"
PERTURBSEQ_EFFECT_COL = "approx_log2FC"  # use beta_hat if preferred

# Conditions for UMAP and treatment-effect displays.
FOCUS_CONDITIONS = [
    "D2_LM", "D2_LM+TDI", "D2_LM+IAG",
    "D3_AE", "D3_AE+TDI", "D3_AE+IAG",
    "D2_LM+TDI->D3_AE+TDI",
    "D2_LM+IAG->D3_AE+IAG",
]

# UMAP options
UMAP_KEY = "X_umap"
COMPUTE_UMAP_IF_MISSING = False  # turn on only if you want to compute a new UMAP
MAX_CELLS_PER_CONDITION_FOR_UMAP = 8000
RANDOM_SEED = 0
POINT_SIZE = 0.6
POINT_ALPHA = 0.85

# Marker genes/modules for UMAP and marker score panels.
MARKER_PROGRAM_GENES = {
    "P15 vascular progenitor": ["TAL1", "FLI1", "SOX7"],
    "P23 vascular LM/EC induction": ["KDR", "ETV2"],
    "P17 epithelial/LM alt": ["KRT8", "KRT18"],
    "P2 arterial EC": ["MECOM", "NRP1", "EFNB2"],
    "P34 YAP/mesenchyme": ["ANKRD1", "KRT8", "KRT18"],
    "P24 GATA4/RSPO2": ["GATA4", "RSPO2"],
}

# %% [markdown]
# ## Helpers

# %%
def infer_condition_from_sample_name(sample_name: str) -> str:
    """Strip trailing _repN from sample_name to get condition."""
    m = re.match(r"(.+)_rep\d+$", str(sample_name))
    return m.group(1) if m else str(sample_name)


def program_id(program: str) -> int:
    m = re.search(r"(\d+)$", str(program))
    if m is None:
        raise ValueError(f"Cannot parse program id from {program!r}")
    return int(m.group(1))


def sanitize_filename(x: str) -> str:
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))
    return re.sub(r"_+", "_", x).strip("_")


def ensure_condition_column(obs: pd.DataFrame) -> pd.DataFrame:
    obs = obs.copy()
    if "condition" not in obs.columns:
        if "sample_name" not in obs.columns:
            raise KeyError("Need either obs['condition'] or obs['sample_name']")
        obs["condition"] = obs["sample_name"].map(infer_condition_from_sample_name)
    obs["condition"] = obs["condition"].astype(str)
    if "sample_name" in obs.columns:
        obs["sample_name"] = obs["sample_name"].astype(str)
    return obs


def get_program_names(adata, usage):
    meta = adata.uns.get("perturbseq_cNMF_projection", {})
    names = []
    if isinstance(meta, dict) and "program_names" in meta:
        names = [str(x) for x in meta["program_names"]]
    if len(names) != usage.shape[1]:
        names = [f"P{i}" for i in range(usage.shape[1])]
    order = np.argsort([program_id(p) for p in names])
    return [names[i] for i in order], order


def get_var_index_for_gene_symbol(adata, symbol: str):
    """Return var index/indices for a gene symbol using adata.var['gene'] if available, else var_names."""
    symbol = str(symbol)
    if "gene" in adata.var.columns:
        matches = np.flatnonzero(adata.var["gene"].astype(str).to_numpy() == symbol)
        if len(matches) > 0:
            return matches.tolist()
    matches = np.flatnonzero(adata.var_names.astype(str).to_numpy() == symbol)
    return matches.tolist()


def get_log2cpm_for_genes(adata, gene_symbols, cell_idx=None):
    """Extract log2(CPM + 1) for a small list of genes without loading all genes."""
    if cell_idx is None:
        cell_idx = np.arange(adata.n_obs)
    total_counts = pd.to_numeric(adata.obs["total_counts"], errors="coerce").to_numpy()[cell_idx]
    total_counts = np.where(total_counts > 0, total_counts, np.nan)

    values = {}
    missing = []
    for gene in gene_symbols:
        var_idx = get_var_index_for_gene_symbol(adata, gene)
        if len(var_idx) == 0:
            missing.append(gene)
            continue
        sub = adata[cell_idx, var_idx].X
        if sparse.issparse(sub):
            arr = np.asarray(sub.sum(axis=1)).ravel()
        else:
            arr = np.asarray(sub).sum(axis=1).ravel()
        cpm = arr / total_counts * 1e6
        values[gene] = np.log2(cpm + 1)
    return pd.DataFrame(values, index=adata.obs_names[cell_idx].astype(str)), missing


def average_marker_module(log2cpm_df: pd.DataFrame, genes):
    present = [g for g in genes if g in log2cpm_df.columns]
    if len(present) == 0:
        return pd.Series(np.nan, index=log2cpm_df.index), []
    return log2cpm_df[present].mean(axis=1), present


def save_fig(fig, out_prefix: Path):
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in ["svg", "pdf", "png"]:
        fig.savefig(out_prefix.with_suffix(f".{ext}"), dpi=300 if ext == "png" else None, bbox_inches="tight")

# %% [markdown]
# ## Load projected time-course AnnData

# %%
adata = sc.read_h5ad(H5AD_PATH)
adata.obs = ensure_condition_column(adata.obs)
usage = np.asarray(adata.obsm["X_perturbseq_cNMF_usage"])
program_names, program_order = get_program_names(adata, usage)
usage = usage[:, program_order]
program_to_col = {p: i for i, p in enumerate(program_names)}

print(adata)
print("condition values:", sorted(adata.obs["condition"].unique().tolist()))
print("programs:", program_names[:5], "...", program_names[-5:])
print("UMAP present:", UMAP_KEY in adata.obsm)

# %% [markdown]
# ## Optional: compute UMAP if missing
#
# This can take time. It writes the UMAP back into memory only; save manually if needed.

# %%
if UMAP_KEY not in adata.obsm:
    if not COMPUTE_UMAP_IF_MISSING:
        print(f"No {UMAP_KEY} found. Set COMPUTE_UMAP_IF_MISSING=True to compute one, or add UMAP coordinates to the h5ad.")
    else:
        print("Computing UMAP. This can be slow for 200k+ cells.")
        if "X_pca" not in adata.obsm:
            tmp = adata.copy()
            sc.pp.normalize_total(tmp, target_sum=1e4)
            sc.pp.log1p(tmp)
            sc.pp.highly_variable_genes(tmp, n_top_genes=3000, subset=True)
            sc.pp.scale(tmp, max_value=10)
            sc.tl.pca(tmp, n_comps=50)
            adata.obsm["X_pca"] = tmp.obsm["X_pca"]
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)
        sc.tl.umap(adata)

# %% [markdown]
# ## UMAP panel: selected programs and marker modules
#
# This makes UMAP feature plots for selected conditions. Points are rasterized to keep SVGs reasonable; text remains editable.

# %%
def downsample_indices_by_condition(obs, conditions, max_per_condition=8000, seed=0):
    rng = np.random.default_rng(seed)
    keep = []
    cond_arr = obs["condition"].astype(str).to_numpy()
    for cond in conditions:
        idx = np.flatnonzero(cond_arr == cond)
        if len(idx) == 0:
            continue
        if len(idx) > max_per_condition:
            idx = rng.choice(idx, size=max_per_condition, replace=False)
        keep.extend(idx.tolist())
    return np.array(sorted(keep), dtype=int)


def plot_umap_feature_facets(feature_values, feature_name, conditions=FOCUS_CONDITIONS, out_prefix=None, cmap="viridis"):
    if UMAP_KEY not in adata.obsm:
        print(f"Skipping {feature_name}: no {UMAP_KEY} in adata.obsm")
        return

    idx = downsample_indices_by_condition(adata.obs, conditions, MAX_CELLS_PER_CONDITION_FOR_UMAP, RANDOM_SEED)
    if len(idx) == 0:
        print("No cells for selected conditions.")
        return

    coords = np.asarray(adata.obsm[UMAP_KEY])[idx]
    values = np.asarray(feature_values)[idx]
    conds = adata.obs["condition"].astype(str).to_numpy()[idx]
    present_conditions = [c for c in conditions if np.any(conds == c)]

    n = len(present_conditions)
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.0*ncols, 3.0*nrows), squeeze=False)

    finite = np.isfinite(values)
    vmax = np.nanpercentile(values[finite], 99) if finite.any() else 1
    vmin = np.nanpercentile(values[finite], 1) if finite.any() else 0
    if np.isclose(vmin, vmax):
        vmax = vmin + 1e-9

    for ax, cond in zip(axes.ravel(), present_conditions):
        m = conds == cond
        sca = ax.scatter(
            coords[m, 0], coords[m, 1],
            c=values[m],
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            linewidths=0,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
        )
        ax.set_title(cond)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_aspect("equal")

    for ax in axes.ravel()[len(present_conditions):]:
        ax.axis("off")

    fig.suptitle(feature_name, y=1.02)
    cbar = fig.colorbar(sca, ax=axes.ravel().tolist(), fraction=0.02, pad=0.01)
    cbar.set_label(feature_name)
    fig.tight_layout()

    if out_prefix is not None:
        save_fig(fig, out_prefix)
    plt.show()
    plt.close(fig)

# Program UMAPs
for program in ["P15", "P17", "P23", "P2", "P34", "P24"]:
    if program in program_to_col:
        plot_umap_feature_facets(
            usage[:, program_to_col[program]],
            feature_name=f"Projected {program} usage",
            out_prefix=OUT_DIR / "umap_programs" / f"UMAP_projected_{program}_selected_conditions",
        )

# Marker module UMAPs
all_marker_genes = sorted(set(g for genes in MARKER_PROGRAM_GENES.values() for g in genes))
if UMAP_KEY in adata.obsm:
    idx_all = np.arange(adata.n_obs)
    marker_expr, missing_genes = get_log2cpm_for_genes(adata, all_marker_genes, idx_all)
    print("Missing marker genes:", missing_genes)
    for module_name, genes in MARKER_PROGRAM_GENES.items():
        vals, present = average_marker_module(marker_expr, genes)
        print(module_name, "present genes:", present)
        plot_umap_feature_facets(
            vals.to_numpy(),
            feature_name=f"{module_name} marker module log2CPM",
            out_prefix=OUT_DIR / "umap_marker_modules" / f"UMAP_{sanitize_filename(module_name)}_marker_module_selected_conditions",
        )

# %% [markdown]
# ## Load perturb-seq Hippo CRISPRi program effects

# %%
assoc = pd.read_excel(CNMF_WORKBOOK_PATH, sheet_name=PROGRAM_ASSOC_SHEET)
assoc = assoc.copy()
assoc["program"] = "P" + assoc["program_name"].astype(int).astype(str)
assoc = assoc[assoc["program"].isin(SELECTED_PROGRAMS)].copy()
assoc = assoc[assoc["target_name"].isin(HIPPO_KD_TARGETS)].copy()
assoc = assoc[assoc["Sample"].isin(["D2", "D3"])].copy()
assoc[PERTURBSEQ_EFFECT_COL] = pd.to_numeric(assoc[PERTURBSEQ_EFFECT_COL], errors="coerce")

print("assoc shape:", assoc.shape)
print("targets:", sorted(assoc["target_name"].unique().tolist()))
print("programs:", sorted(assoc["program"].unique().tolist(), key=program_id))

hippo_effect_individual = (
    assoc.groupby(["Sample", "target_name", "program"], observed=True)[PERTURBSEQ_EFFECT_COL]
    .mean()
    .reset_index()
    .rename(columns={"Sample": "perturbseq_day", PERTURBSEQ_EFFECT_COL: "perturbseq_HippoKD_effect"})
)

hippo_effect_mean = (
    hippo_effect_individual.groupby(["perturbseq_day", "program"], observed=True)["perturbseq_HippoKD_effect"]
    .agg(["mean", "std", "count"])
    .reset_index()
    .rename(columns={"mean": "perturbseq_HippoKD_mean_effect", "std": "perturbseq_HippoKD_sd", "count": "n_targets"})
)

hippo_effect_individual.to_csv(OUT_DIR / "perturbseq_HippoKD_individual_target_effects_selected_programs.tsv", sep="\t", index=False)
hippo_effect_mean.to_csv(OUT_DIR / "perturbseq_HippoKD_mean_effects_selected_programs.tsv", sep="\t", index=False)
display(hippo_effect_mean.head())

# %% [markdown]
# ## Load small-molecule program effects

# %%
if SMALL_MOL_CLR_EFFECT_PATH.exists():
    sm = pd.read_csv(SMALL_MOL_CLR_EFFECT_PATH, sep="\t")
    effect_col = "mean_diff_clr"
    effect_label = "Small-molecule CLR effect"
    paired_path = SMALL_MOL_CLR_PAIRED_PATH
elif SMALL_MOL_SCORE_EFFECT_PATH.exists():
    sm = pd.read_csv(SMALL_MOL_SCORE_EFFECT_PATH, sep="\t")
    effect_col = "mean_diff"
    effect_label = "Small-molecule mean usage difference"
    paired_path = SMALL_MOL_SCORE_PAIRED_PATH
else:
    raise FileNotFoundError("No small-molecule effect table found.")

sm = sm.copy()
sm = sm[sm["program"].isin(SELECTED_PROGRAMS)].copy()
sm[effect_col] = pd.to_numeric(sm[effect_col], errors="coerce")

# Keep contrasts most relevant for Fig. 4a.
contrast_order = [
    "D2_LM+TDI",
    "D2_LM+IAG",
    "D3_AE+TDI",
    "D3_AE+IAG",
    "D2_LM+TDI->D3_AE+TDI",
    "D2_LM+IAG->D3_AE+IAG",
]
sm = sm[sm["contrast"].isin(contrast_order)].copy()
sm["contrast"] = pd.Categorical(sm["contrast"], categories=contrast_order, ordered=True)
sm.to_csv(OUT_DIR / "small_molecule_effects_selected_programs.tsv", sep="\t", index=False)

print("small molecule effect source:", effect_label)
display(sm.head())

# %% [markdown]
# ## Plot: program effects as dot/lollipop comparison, not heatmap

# %%
def plot_effect_dotpanel(hippo_mean, sm_df, out_prefix):
    programs = [p for p in SELECTED_PROGRAMS if p in set(hippo_mean["program"]) or p in set(sm_df["program"])]
    programs = sorted(programs, key=program_id)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, max(5, 0.32*len(programs))), sharey=True)

    y_pos = np.arange(len(programs))

    # Left: perturb-seq Hippo KD mean effects, D2 and D3.
    ax = axes[0]
    for day, offset in [("D2", -0.12), ("D3", 0.12)]:
        sub = hippo_mean[hippo_mean["perturbseq_day"] == day].set_index("program").reindex(programs)
        vals = sub["perturbseq_HippoKD_mean_effect"].to_numpy(dtype=float)
        ax.axvline(0, color="0.7", lw=0.8)
        ax.scatter(vals, y_pos + offset, s=30, label=f"Perturb-seq Hippo KD {day}")
        for y, v in zip(y_pos + offset, vals):
            if np.isfinite(v):
                ax.plot([0, v], [y, y], lw=0.7, alpha=0.7)

    ax.set_title("Perturb-seq Hippo KD")
    ax.set_xlabel(f"Program effect ({PERTURBSEQ_EFFECT_COL})")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(programs)
    ax.legend(frameon=False, fontsize=6)

    # Right: small molecules.
    ax = axes[1]
    offsets = np.linspace(-0.25, 0.25, len(contrast_order))
    for contrast, offset in zip(contrast_order, offsets):
        sub = sm_df[sm_df["contrast"] == contrast].set_index("program").reindex(programs)
        vals = sub[effect_col].to_numpy(dtype=float)
        ax.axvline(0, color="0.7", lw=0.8)
        ax.scatter(vals, y_pos + offset, s=22, label=contrast)
        for y, v in zip(y_pos + offset, vals):
            if np.isfinite(v):
                ax.plot([0, v], [y, y], lw=0.5, alpha=0.5)

    ax.set_title("Small molecule")
    ax.set_xlabel(effect_label)
    ax.legend(frameon=False, fontsize=5, loc="center left", bbox_to_anchor=(1.0, 0.5))

    for ax in axes:
        ax.grid(axis="y", color="0.92", lw=0.5)
        ax.set_ylim(len(programs)-0.5, -0.5)
        for spine in ax.spines.values():
            spine.set_visible(False)

    fig.suptitle("Effects on selected D2/D3 programs")
    fig.tight_layout()
    save_fig(fig, out_prefix)
    plt.show()
    plt.close(fig)

plot_effect_dotpanel(
    hippo_effect_mean,
    sm,
    OUT_DIR / "Fig4a_program_effect_dotpanel_HippoKD_vs_small_molecule_selected_programs",
)

# %% [markdown]
# ## Plot: selected program scatter, small molecule vs perturb-seq Hippo KD

# %%
def plot_smallmol_vs_hippokd_selected_scatter(hippo_mean, sm_df, out_prefix):
    # Compare small molecule against D2 and D3 HippoKD references for selected programs.
    merged_records = []
    for contrast in contrast_order:
        sub_sm = sm_df[sm_df["contrast"] == contrast][["contrast", "program", effect_col]]
        for day in ["D2", "D3"]:
            sub_h = hippo_mean[hippo_mean["perturbseq_day"] == day][["program", "perturbseq_HippoKD_mean_effect"]]
            m = sub_sm.merge(sub_h, on="program", how="inner")
            m["perturbseq_day"] = day
            merged_records.append(m)
    merged = pd.concat(merged_records, ignore_index=True)
    merged.to_csv(OUT_DIR / "small_molecule_vs_HippoKD_selected_program_scatter_values.tsv", sep="\t", index=False)

    nrows = len(contrast_order)
    fig, axes = plt.subplots(nrows, 2, figsize=(7.6, 2.05*nrows), sharex=True, sharey=True)

    for i, contrast in enumerate(contrast_order):
        for j, day in enumerate(["D2", "D3"]):
            ax = axes[i, j]
            sub = merged[(merged["contrast"] == contrast) & (merged["perturbseq_day"] == day)]
            ax.axhline(0, color="0.75", lw=0.8)
            ax.axvline(0, color="0.75", lw=0.8)
            ax.scatter(sub["perturbseq_HippoKD_mean_effect"], sub[effect_col], s=28)
            for _, row in sub.iterrows():
                if row["program"] in KEY_PROGRAMS:
                    ax.text(row["perturbseq_HippoKD_mean_effect"], row[effect_col], row["program"], fontsize=6)
            ax.set_title(f"{contrast} vs HippoKD {day}", fontsize=7)
            if i == nrows - 1:
                ax.set_xlabel(f"HippoKD {day} effect")
            if j == 0:
                ax.set_ylabel(effect_label)
            for spine in ax.spines.values():
                spine.set_visible(False)

    fig.tight_layout()
    save_fig(fig, out_prefix)
    plt.show()
    plt.close(fig)

plot_smallmol_vs_hippokd_selected_scatter(
    hippo_effect_mean,
    sm,
    OUT_DIR / "Fig4a_small_molecule_vs_HippoKD_selected_program_scatter_D2_D3_refs",
)

# %% [markdown]
# ## Optional: UMAP marker/program panels are already saved above if UMAP exists.

# %%
print("Saved outputs to:", OUT_DIR)
for p in sorted(OUT_DIR.glob("*.svg")):
    print(p.name)
