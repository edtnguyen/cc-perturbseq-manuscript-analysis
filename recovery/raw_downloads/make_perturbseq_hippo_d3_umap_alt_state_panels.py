# %% [markdown]
# # Perturb-seq D3 Hippo UMAP panels: alternative states
#
# This notebook/script makes UMAP panels from the perturb-seq dataset itself,
# focused on whether Hippo perturbations occupy D3 alternative/YAP-like states.
#
# Outputs:
# - D3 UMAP colored by selected perturb-seq cNMF programs
# - D3 UMAP colored by EC score, alternative/YAP score, and alternative-minus-EC score
# - D3 UMAP highlighting Hippo KD carrier cells vs controls/other cells
# - D3 UMAP panels for individual Hippo KD targets
# - D3 target-level summary dot/bar plots for selected programs and EC/alternative scores
#
# Edit the USER SETTINGS cell if your day/target columns differ.

# %%
from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import scanpy as sc
import muon as mu
from scipy import sparse

warnings.filterwarnings("ignore")

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",  # keep text editable in Illustrator
})

# %% [markdown]
# ## User settings

# %%
# Perturb-seq cNMF h5mu from the perturb-seq screen.
H5MU_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/adata/cNMF_50_0_2.h5mu"
)

# Integrated perturb-seq h5ad used previously for UMAP and obs metadata.
# If your cNMF h5mu already has X_umap and target metadata, you can set USE_INTEGRATED_H5AD_FOR_UMAP_AND_OBS = False.
USE_INTEGRATED_H5AD_FOR_UMAP_AND_OBS = True
ADATA_INT_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
    "2.preprocessing_tonys_pipeline/3.integrated_days/"
    "3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad"
)

OUT_DIR = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_hippo_d3_umap_alt_state_panels")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# If set to None, the script tries to infer columns and prints candidates.
DAY_COL = None
TARGET_COL = None

# D3 values are detected by string contains "D3" by default.
D3_VALUE_REGEX = r"D3"

# Hippo KD targets to highlight.
HIPPO_TARGETS = ["LATS2", "NF2", "SAV1", "TAOK1", "WWC3", "AMOTL2"]

# Control labels, if TARGET_COL has a target/perturbation label.
CONTROL_REGEX = r"NT|NON[-_ ]?TARGET|CONTROL|SAFE|SCRAMBLE|INTERGENIC|NEG"

# Selected programs. These are zero-based perturb-seq cNMF program IDs.
D2_BRANCH_PROGRAMS = ["P14", "P15", "P17", "P19", "P23"]
D3_EC_PROGRAMS = ["P2", "P3", "P8", "P11", "P42"]
D3_ALT_PROGRAMS = ["P16", "P24", "P27", "P33", "P34", "P35"]
KEY_PROGRAMS_TO_UMAP = ["P15", "P17", "P23", "P2", "P34", "P24"]

# Use normalized program scores for UMAP features. If False, uses raw cnmf.X.
ROW_NORMALIZE_CNMFSCORES = True

# UMAP plotting options.
POINT_SIZE_ALL = 0.15
POINT_SIZE_HIGHLIGHT = 0.8
UMAP_ALPHA_ALL = 0.65
UMAP_ALPHA_HIGHLIGHT = 0.95
RASTERIZE_POINTS = True  # keeps SVG compact; text remains editable
MAX_CELLS_FOR_CONTINUOUS_UMAP = None  # set e.g. 120000 for faster plotting; None = all D3 cells
RANDOM_SEED = 0

# Color limits for continuous UMAPs.
FEATURE_VMIN_QUANTILE = 0.00
FEATURE_VMAX_QUANTILE = 0.995

# %% [markdown]
# ## Helper functions

# %%
def to_dense_array(x):
    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


def sanitize_name(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x)).strip("_")


def parse_program_id(program):
    return int(str(program).replace("P", ""))


def infer_day_column(obs):
    candidates = []
    for col in obs.columns:
        s = obs[col].astype(str)
        vals = s.dropna().unique()
        if len(vals) > 200:
            continue
        d_hits = sum(bool(re.search(r"D[0-9]", v, flags=re.IGNORECASE)) for v in vals)
        day_hits = sum("day" in v.lower() for v in vals)
        if d_hits > 0 or day_hits > 0:
            candidates.append((col, len(vals), d_hits, day_hits, list(vals[:10])))
    return candidates


def infer_target_column(obs, hippo_targets):
    hippo = [t.upper() for t in hippo_targets]
    candidates = []
    for col in obs.columns:
        if obs[col].dtype.name not in ["category", "object", "string"]:
            continue
        s = obs[col].astype(str).str.upper()
        vals = s.dropna().unique()
        if len(vals) > 2000:
            # target columns should not usually have one unique value per cell.
            continue
        n_hippo = int(sum(s.str.contains("|".join(hippo), regex=True, na=False)))
        n_control = int(sum(s.str.contains(CONTROL_REGEX, regex=True, na=False)))
        if n_hippo > 0 or n_control > 0:
            candidates.append((col, len(vals), n_hippo, n_control, list(vals[:10])))
    candidates = sorted(candidates, key=lambda x: (x[2], x[3]), reverse=True)
    return candidates


def get_target_series(obs, target_col, hippo_targets):
    if target_col is None:
        candidates = infer_target_column(obs, hippo_targets)
        print("Candidate target columns:")
        for c in candidates[:20]:
            print(c)
        if not candidates:
            raise ValueError(
                "Could not infer target column. Set TARGET_COL manually after inspecting obs columns."
            )
        target_col = candidates[0][0]
        print(f"Using inferred TARGET_COL={target_col!r}")
    return obs[target_col].astype(str), target_col


def get_day_series(obs, day_col):
    if day_col is None:
        candidates = infer_day_column(obs)
        print("Candidate day columns:")
        for c in candidates[:20]:
            print(c)
        if not candidates:
            raise ValueError("Could not infer day column. Set DAY_COL manually after inspecting obs columns.")
        day_col = candidates[0][0]
        print(f"Using inferred DAY_COL={day_col!r}")
    return obs[day_col].astype(str), day_col


def make_target_group(target_string_series):
    s = target_string_series.astype(str).str.upper()
    hippo_regex = "|".join([re.escape(t.upper()) for t in HIPPO_TARGETS])
    out = pd.Series("Other", index=target_string_series.index, dtype=object)
    out[s.str.contains(CONTROL_REGEX, regex=True, na=False)] = "Control/NT"
    out[s.str.contains(hippo_regex, regex=True, na=False)] = "Hippo KD"
    return out


def get_individual_hippo_target(target_string):
    s = str(target_string).upper()
    hits = [t for t in HIPPO_TARGETS if re.search(rf"\b{re.escape(t.upper())}\b", s)]
    if hits:
        return hits[0]
    # fallback substring match for combined guide strings
    hits = [t for t in HIPPO_TARGETS if t.upper() in s]
    return hits[0] if hits else None


def robust_feature_limits(values, qmin=0.0, qmax=0.995):
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0, 1
    vmin = float(np.quantile(finite, qmin))
    vmax = float(np.quantile(finite, qmax))
    if np.isclose(vmin, vmax):
        vmax = vmin + 1e-9
    return vmin, vmax


def plot_umap_continuous(df, value_col, title, out_prefix, cmap="viridis", overlay_hippo=True):
    plot_df = df.copy()
    if MAX_CELLS_FOR_CONTINUOUS_UMAP is not None and len(plot_df) > MAX_CELLS_FOR_CONTINUOUS_UMAP:
        plot_df = plot_df.sample(MAX_CELLS_FOR_CONTINUOUS_UMAP, random_state=RANDOM_SEED)

    vals = plot_df[value_col].to_numpy(dtype=float)
    vmin, vmax = robust_feature_limits(vals, FEATURE_VMIN_QUANTILE, FEATURE_VMAX_QUANTILE)
    order = np.argsort(vals)
    plot_df = plot_df.iloc[order]

    fig, ax = plt.subplots(figsize=(4.5, 4.2))
    sca = ax.scatter(
        plot_df["UMAP1"], plot_df["UMAP2"],
        c=plot_df[value_col],
        s=POINT_SIZE_ALL,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        alpha=UMAP_ALPHA_ALL,
        linewidths=0,
        rasterized=RASTERIZE_POINTS,
    )

    if overlay_hippo and "target_group" in plot_df.columns:
        hippo = plot_df[plot_df["target_group"] == "Hippo KD"]
        if len(hippo) > 0:
            ax.scatter(
                hippo["UMAP1"], hippo["UMAP2"],
                facecolors="none",
                edgecolors="black",
                s=POINT_SIZE_HIGHLIGHT * 3,
                linewidths=0.25,
                alpha=0.5,
                rasterized=RASTERIZE_POINTS,
            )

    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(value_col)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT_DIR / f"{out_prefix}.{ext}", bbox_inches="tight")
    plt.close(fig)


def plot_umap_target_groups(df, title, out_prefix):
    fig, ax = plt.subplots(figsize=(4.5, 4.2))
    # Use neutral/standard colors. Edit if needed for journal palette.
    group_order = ["Other", "Control/NT", "Hippo KD"]
    colors = {"Other": "#d9d9d9", "Control/NT": "#4d4d4d", "Hippo KD": "#d62728"}
    sizes = {"Other": 0.10, "Control/NT": 0.20, "Hippo KD": 0.80}
    alphas = {"Other": 0.25, "Control/NT": 0.55, "Hippo KD": 0.9}
    zorders = {"Other": 1, "Control/NT": 2, "Hippo KD": 3}

    for group in group_order:
        sub = df[df["target_group"] == group]
        if len(sub) == 0:
            continue
        ax.scatter(
            sub["UMAP1"], sub["UMAP2"],
            c=colors[group],
            s=sizes[group],
            alpha=alphas[group],
            linewidths=0,
            label=f"{group} (n={len(sub):,})",
            rasterized=RASTERIZE_POINTS,
            zorder=zorders[group],
        )

    ax.legend(frameon=False, loc="upper right", fontsize=6, markerscale=4)
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT_DIR / f"{out_prefix}.{ext}", bbox_inches="tight")
    plt.close(fig)

# %% [markdown]
# ## Load cNMF scores, UMAP, obs metadata, and align cells

# %%
print("Loading cNMF h5mu:", H5MU_PATH)
mdata = mu.read_h5mu(H5MU_PATH)
cnmf = mdata.mod["cNMF"]
cnmf_cells = pd.Index(cnmf.obs_names.astype(str))

scores = to_dense_array(cnmf.X).astype(np.float32)
if ROW_NORMALIZE_CNMFSCORES:
    row_sums = scores.sum(axis=1, keepdims=True)
    scores = np.divide(scores, row_sums, out=np.zeros_like(scores), where=row_sums > 0)

program_names = [f"P{i}" for i in range(scores.shape[1])]
score_df = pd.DataFrame(scores, index=cnmf_cells, columns=program_names)
print("cNMF score matrix:", score_df.shape)

if USE_INTEGRATED_H5AD_FOR_UMAP_AND_OBS:
    print("Loading integrated h5ad for UMAP/obs:", ADATA_INT_PATH)
    adata_int = sc.read_h5ad(ADATA_INT_PATH, backed="r")
    source_cells = pd.Index(adata_int.obs_names.astype(str))
    if not cnmf_cells.isin(source_cells).all():
        missing = cnmf_cells[~cnmf_cells.isin(source_cells)][:10].tolist()
        raise ValueError(f"Some cNMF cells are missing from integrated h5ad. Examples: {missing}")
    pos = source_cells.get_indexer(cnmf_cells)
    obs = adata_int.obs.iloc[pos].copy()
    if "X_umap" not in adata_int.obsm:
        raise KeyError(f"X_umap not found in ADATA_INT. obsm keys: {list(adata_int.obsm.keys())}")
    umap = np.asarray(adata_int.obsm["X_umap"])[pos, :2]
    try:
        adata_int.file.close()
    except Exception:
        pass
else:
    obs = cnmf.obs.copy()
    if "X_umap" not in cnmf.obsm:
        raise KeyError("X_umap not found in cnmf.obsm. Set USE_INTEGRATED_H5AD_FOR_UMAP_AND_OBS=True.")
    umap = np.asarray(cnmf.obsm["X_umap"])

obs.index = cnmf_cells
umap_df = pd.DataFrame(umap, index=cnmf_cells, columns=["UMAP1", "UMAP2"])

print("obs shape:", obs.shape)
print("UMAP shape:", umap_df.shape)
print("obs columns:")
print(list(obs.columns))

# %% [markdown]
# ## Identify D3 cells and Hippo target groups

# %%
day_series, DAY_COL_USED = get_day_series(obs, DAY_COL)
d3_mask = day_series.str.contains(D3_VALUE_REGEX, regex=True, case=False, na=False)
print("DAY_COL_USED:", DAY_COL_USED)
print("D3 cells:", int(d3_mask.sum()), "of", len(d3_mask))
print("Day value counts:")
print(day_series.value_counts(dropna=False).head(20))

target_series, TARGET_COL_USED = get_target_series(obs, TARGET_COL, HIPPO_TARGETS)
target_group = make_target_group(target_series)
hippo_target = target_series.map(get_individual_hippo_target)

print("TARGET_COL_USED:", TARGET_COL_USED)
print("target_group counts, all cells:")
print(target_group.value_counts(dropna=False))
print("target_group counts, D3 cells:")
print(target_group[d3_mask].value_counts(dropna=False))
print("individual Hippo target counts, D3 cells:")
print(hippo_target[d3_mask].value_counts(dropna=False))

# %% [markdown]
# ## Build D3 plotting table

# %%
d3_idx = cnmf_cells[d3_mask.to_numpy()]
plot_df = pd.concat([
    umap_df.loc[d3_idx],
    score_df.loc[d3_idx],
], axis=1)
plot_df["day"] = day_series.loc[d3_idx].astype(str).to_numpy()
plot_df["target_raw"] = target_series.loc[d3_idx].astype(str).to_numpy()
plot_df["target_group"] = target_group.loc[d3_idx].astype(str).to_numpy()
plot_df["hippo_target"] = hippo_target.loc[d3_idx].fillna("not_Hippo").astype(str).to_numpy()

# Composite scores.
for p in D3_EC_PROGRAMS + D3_ALT_PROGRAMS + D2_BRANCH_PROGRAMS:
    if p not in plot_df.columns:
        raise ValueError(f"Program {p} missing from score_df columns")

plot_df["D3_EC_score"] = plot_df[D3_EC_PROGRAMS].mean(axis=1)
plot_df["D3_alt_YAP_score"] = plot_df[D3_ALT_PROGRAMS].mean(axis=1)
plot_df["D3_alt_minus_EC_score"] = plot_df["D3_alt_YAP_score"] - plot_df["D3_EC_score"]
plot_df["P34_minus_P2"] = plot_df["P34"] - plot_df["P2"]

print(plot_df.shape)
print(plot_df[["target_group", "hippo_target"]].value_counts().head(20))
plot_df.to_csv(OUT_DIR / "D3_perturbseq_umap_plotting_table_selected_scores.tsv.gz", sep="\t", compression="gzip")

# %% [markdown]
# ## UMAP: target group and individual Hippo targets

# %%
plot_umap_target_groups(
    plot_df,
    title="Perturb-seq D3 UMAP: Hippo KD carriers",
    out_prefix="D3_UMAP_target_group_HippoKD_vs_control_other",
)

# Individual Hippo target panels.
for target in HIPPO_TARGETS:
    df_t = plot_df.copy()
    df_t["target_group"] = np.where(df_t["hippo_target"] == target, target, "Other")
    # Reuse category plotting with a simple custom plot for this target.
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    other = df_t[df_t["target_group"] == "Other"]
    hit = df_t[df_t["target_group"] == target]
    ax.scatter(other["UMAP1"], other["UMAP2"], c="#d9d9d9", s=0.10, alpha=0.25, linewidths=0, rasterized=RASTERIZE_POINTS)
    ax.scatter(hit["UMAP1"], hit["UMAP2"], c="#d62728", s=0.7, alpha=0.9, linewidths=0, rasterized=RASTERIZE_POINTS, label=f"{target} (n={len(hit):,})")
    ax.legend(frameon=False, loc="upper right", fontsize=6, markerscale=4)
    ax.set_title(f"D3 perturb-seq UMAP: {target} KD carriers")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for spine in ax.spines.values(): spine.set_visible(False)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT_DIR / f"D3_UMAP_Hippo_target_{target}.{ext}", bbox_inches="tight")
    plt.close(fig)

# %% [markdown]
# ## UMAP: selected programs and alternative/EC composite scores

# %%
features_to_plot = []
features_to_plot += KEY_PROGRAMS_TO_UMAP
features_to_plot += ["D3_EC_score", "D3_alt_YAP_score", "D3_alt_minus_EC_score", "P34_minus_P2"]

for feat in features_to_plot:
    plot_umap_continuous(
        plot_df,
        value_col=feat,
        title=f"Perturb-seq D3 UMAP: {feat}",
        out_prefix=f"D3_UMAP_feature_{sanitize_name(feat)}",
        cmap="viridis",
        overlay_hippo=True,
    )

# %% [markdown]
# ## Target-level summary: EC vs alternative score by target group

# %%
# Summarize selected scores by target. This is useful beside the UMAP as a quantitative panel.
summary_features = ["D3_EC_score", "D3_alt_YAP_score", "D3_alt_minus_EC_score", "P2", "P34", "P34_minus_P2"]
summary = (
    plot_df
    .groupby("hippo_target", observed=True)[summary_features]
    .agg(["mean", "median", "count"])
)
summary.to_csv(OUT_DIR / "D3_target_level_selected_score_summary.tsv", sep="\t")
display(summary)

# For cleaner plotting, compare only individual Hippo targets and non-Hippo/background groups.
plot_targets = HIPPO_TARGETS + ["not_Hippo"]
summary_mean = plot_df.groupby("hippo_target", observed=True)[summary_features].mean().reindex(plot_targets)
summary_mean.to_csv(OUT_DIR / "D3_target_level_selected_score_means.tsv", sep="\t")

# Dot/lollipop-style summary.
for feat in ["D3_EC_score", "D3_alt_YAP_score", "D3_alt_minus_EC_score", "P2", "P34", "P34_minus_P2"]:
    vals = summary_mean[feat]
    fig, ax = plt.subplots(figsize=(4.8, 2.8))
    y = np.arange(len(vals))
    ax.axvline(float(vals.loc["not_Hippo"]), color="#777777", lw=0.8, ls="--", label="not_Hippo mean")
    ax.scatter(vals, y, s=25)
    for yi, target in enumerate(vals.index):
        ax.text(vals.loc[target], yi, "  " + target, va="center", fontsize=7)
    ax.set_yticks([])
    ax.set_xlabel(feat)
    ax.set_title(f"D3 perturb-seq: target-level mean {feat}")
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT_DIR / f"D3_target_mean_{sanitize_name(feat)}.{ext}", bbox_inches="tight")
    plt.close(fig)

print("Saved outputs to", OUT_DIR)
