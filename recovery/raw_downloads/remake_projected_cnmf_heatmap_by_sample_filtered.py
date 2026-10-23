# %% [markdown]
# # Projected perturb-seq cNMF program heatmap by sample_name
#
# This notebook/script loads the saved zero-based projected AnnData/MuData,
# applies condition-specific mitochondrial/ribosomal QC filters, and remakes a
# projected cNMF program x sample_name heatmap with all program names shown.

# %%
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt

try:
    import muon as mu
except ImportError:
    mu = None

# Edit paths here if needed.
OUT_DIR = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_projection")
H5AD_PATH = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad"
H5MU_PATH = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5mu"
PROGRAM_META_PATH = Path("/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/data/program_activity_long_v4.csv")

PLOT_OUT_DIR = OUT_DIR / "heatmaps_after_qc_filter"
PLOT_OUT_DIR.mkdir(parents=True, exist_ok=True)

# Use normalized projected usage by default. Change to raw if needed.
SCORE_KEY = "X_perturbseq_cNMF_usage"       # row-normalized cell x program usage
# SCORE_KEY = "X_perturbseq_cNMF_usage_raw" # unnormalized projected score

SAMPLE_COL = "sample_name"
CONDITION_COL = "condition"
SHOW_PROGRAM_LABELS = True
MAX_PROGRAM_LABEL_CHARS = 90

# Heatmap scaling:
# - "zscore_by_program": highlights relative changes across samples per program.
# - "mean": plots mean projected usage directly.
HEATMAP_SCALE = "zscore_by_program"
# HEATMAP_SCALE = "mean"

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",  # editable text in Illustrator
})


def show_df(df, n=20):
    """Display a dataframe in notebooks, or print in plain Python."""
    try:
        display(df)
    except NameError:
        print(df.head(n).to_string())

print("H5AD:", H5AD_PATH)
print("H5MU:", H5MU_PATH)
print("Output:", PLOT_OUT_DIR)

# %% [markdown]
# ## Load saved projected AnnData/MuData
#
# Load the AnnData in backed mode so the large gene-expression matrix is not read into memory.
# The projected cNMF matrix is only cells x 50 and will be loaded separately from `obsm`.

# %%
if not H5AD_PATH.exists():
    raise FileNotFoundError(f"Missing projected AnnData: {H5AD_PATH}")

adata = sc.read_h5ad(H5AD_PATH, backed="r")
print(adata)
print("obs columns:", list(adata.obs.columns))
print("obsm keys:", list(adata.obsm.keys()))

mdata = None
if mu is not None and H5MU_PATH.exists():
    try:
        mdata = mu.read_h5mu(H5MU_PATH, backed="r")
        print("\nLoaded MuData:")
        print(mdata)
        print("modalities:", list(mdata.mod.keys()))
    except Exception as e:
        print("Could not load MuData; continuing with AnnData only.")
        print(type(e).__name__, e)
elif mu is None:
    print("muon is not installed; continuing with AnnData only.")
else:
    print("MuData file not found; continuing with AnnData only:", H5MU_PATH)

# %% [markdown]
# ## Define QC cutoffs
#
# Cells are removed if either:
#
# ```text
# pct_counts_mt >= condition-specific mt cutoff
# OR
# pct_counts_ribo <= condition-specific ribo cutoff
# ```

# %%
cond_pct_ribo_cutofs = {
    'D0_hPSC': 2.4,
    'D1_MPS': 2.3,
    'D1_MPS+TDI': 2.5,
    'D1_MPS+IAG': 2.2,
    'D1.33_LM': 1.9,
    'D1.33_LM+TDI': 2.1,
    'D1.33_LM+IAG': 2.2,
    'D1.67_LM': 2.0,
    'D1.67_LM+TDI': 2.2,
    'D1.67_LM+IAG': 1.9,
    'D2_LM': 2.0,
    'D2_LM+TDI': 2.4,
    'D2_LM+IAG': 1.9,
    'D3_AE': 1.0,
    'D3_AE+TDI': 2.0,
    'D3_AE+IAG': 1.0,
    'D2_LM+TDI->D3_AE+TDI': 2.5,
    'D2_LM+IAG->D3_AE+IAG': 1.0,
}

cond_pct_mt_cutofs = {
    'D0_hPSC': 5.5,
    'D1_MPS': 3.0,
    'D1_MPS+TDI': 3.0,
    'D1_MPS+IAG': 4,
    'D1.33_LM': 3.5,
    'D1.33_LM+TDI': 3.75,
    'D1.33_LM+IAG': 5,
    'D1.67_LM': 4.0,
    'D1.67_LM+TDI': 5.5,
    'D1.67_LM+IAG': 3.5,
    'D2_LM': 6.5,
    'D2_LM+TDI': 8,
    'D2_LM+IAG': 3.5,
    'D3_AE': 3.25,
    'D3_AE+TDI': 7.0,
    'D3_AE+IAG': 3.75,
    'D2_LM+TDI->D3_AE+TDI': 7,
    'D2_LM+IAG->D3_AE+IAG': 3.75,
}

condition_order = list(cond_pct_mt_cutofs.keys())
assert set(cond_pct_mt_cutofs) == set(cond_pct_ribo_cutofs), "MT and ribo cutoff dictionaries have different condition keys."
print("Conditions with cutoffs:", len(condition_order))

# %% [markdown]
# ## Ensure/derive condition labels
#
# The filter uses `adata.obs[condition]`. If the projected file does not contain a `condition`
# column, this cell tries to infer it from `sample_name` by exact match or substring match.

# %%
obs = adata.obs.copy()

if SAMPLE_COL not in obs.columns:
    raise KeyError(f"SAMPLE_COL={SAMPLE_COL!r} is not present in adata.obs. Available columns: {list(obs.columns)}")


def infer_condition_from_sample(value, valid_conditions):
    """Infer condition from a sample string, preferring exact then longest substring match."""
    s = str(value)
    if s in valid_conditions:
        return s
    for cond in sorted(valid_conditions, key=len, reverse=True):
        if cond in s:
            return cond
    return np.nan

if CONDITION_COL not in obs.columns:
    inferred = obs[SAMPLE_COL].map(lambda x: infer_condition_from_sample(x, condition_order))
    frac_mapped = inferred.notna().mean()
    print(f"{CONDITION_COL!r} not found in obs. Inferred condition from {SAMPLE_COL!r}: {frac_mapped:.3%} mapped.")
    if frac_mapped < 1.0:
        unmapped = obs.loc[inferred.isna(), SAMPLE_COL].astype(str).value_counts().head(50)
        print("Unmapped sample_name values:")
        print(unmapped.to_string())
        raise ValueError(
            "Could not infer condition for all cells. Add adata.obs['condition'] or adjust infer_condition_from_sample()."
        )
    obs[CONDITION_COL] = inferred.astype(str)
else:
    obs[CONDITION_COL] = obs[CONDITION_COL].astype(str)

print("Condition counts:")
show_df(obs[CONDITION_COL].value_counts(dropna=False).rename_axis("condition").reset_index(name="n_cells"), n=100)

# %% [markdown]
# ## Apply condition-specific ribo/mito filters

# %%
required_cols = [CONDITION_COL, SAMPLE_COL, "pct_counts_mt", "pct_counts_ribo"]
missing_cols = [c for c in required_cols if c not in obs.columns]
if missing_cols:
    raise KeyError(f"Missing required columns in adata.obs: {missing_cols}")

observed_conditions = set(obs[CONDITION_COL].astype(str).unique())
missing_cutoffs = sorted(observed_conditions - set(condition_order))
if missing_cutoffs:
    raise KeyError(
        "These observed conditions do not have cutoffs: " + ", ".join(missing_cutoffs)
    )

mt = pd.to_numeric(obs["pct_counts_mt"], errors="coerce")
ribo = pd.to_numeric(obs["pct_counts_ribo"], errors="coerce")
if mt.isna().any() or ribo.isna().any():
    raise ValueError("pct_counts_mt or pct_counts_ribo contains NA/non-numeric values after conversion.")

mt_cutoff = obs[CONDITION_COL].map(cond_pct_mt_cutofs).astype(float)
ribo_cutoff = obs[CONDITION_COL].map(cond_pct_ribo_cutofs).astype(float)

qc_fail = (mt >= mt_cutoff) | (ribo <= ribo_cutoff)
keep = ~qc_fail
keep_mask = keep.to_numpy()

print(f"Cells before filtering: {adata.n_obs:,}")
print(f"Cells removed:          {int(qc_fail.sum()):,} ({qc_fail.mean():.2%})")
print(f"Cells kept:             {int(keep.sum()):,} ({keep.mean():.2%})")

filter_annot = obs[[CONDITION_COL, SAMPLE_COL, "pct_counts_mt", "pct_counts_ribo"]].copy()
filter_annot["mt_cutoff"] = mt_cutoff.to_numpy()
filter_annot["ribo_cutoff"] = ribo_cutoff.to_numpy()
filter_annot["qc_fail"] = qc_fail.to_numpy()
filter_annot["qc_keep"] = keep.to_numpy()

filter_summary_condition = (
    filter_annot
    .groupby(CONDITION_COL, observed=True)
    .agg(
        n_cells=("qc_fail", "size"),
        n_removed=("qc_fail", "sum"),
        frac_removed=("qc_fail", "mean"),
        mt_cutoff=("mt_cutoff", "first"),
        ribo_cutoff=("ribo_cutoff", "first"),
        median_pct_mt=("pct_counts_mt", "median"),
        median_pct_ribo=("pct_counts_ribo", "median"),
    )
)
filter_summary_condition["n_kept"] = filter_summary_condition["n_cells"] - filter_summary_condition["n_removed"]
filter_summary_condition = filter_summary_condition.loc[[c for c in condition_order if c in filter_summary_condition.index]]
show_df(filter_summary_condition, n=100)

filter_summary_sample = (
    filter_annot
    .groupby([CONDITION_COL, SAMPLE_COL], observed=True)
    .agg(
        n_cells=("qc_fail", "size"),
        n_removed=("qc_fail", "sum"),
        frac_removed=("qc_fail", "mean"),
    )
)
filter_summary_sample["n_kept"] = filter_summary_sample["n_cells"] - filter_summary_sample["n_removed"]
show_df(filter_summary_sample.reset_index(), n=100)

filter_summary_condition.to_csv(PLOT_OUT_DIR / "qc_filter_summary_by_condition.tsv", sep="\t")
filter_summary_sample.to_csv(PLOT_OUT_DIR / "qc_filter_summary_by_condition_and_sample_name.tsv", sep="\t")

# Backed AnnData/MuData filtered views for downstream code, without loading adata.X.
adata_f = adata[keep_mask, :]
print("Filtered AnnData view:", adata_f)

mdata_f = None
if mdata is not None:
    try:
        if len(mdata.obs_names) == len(keep_mask) and np.array_equal(
            pd.Index(mdata.obs_names).astype(str).to_numpy(),
            pd.Index(adata.obs_names).astype(str).to_numpy(),
        ):
            mdata_f = mdata[keep_mask, :]
            print("Filtered MuData view:", mdata_f)
        else:
            print("MuData obs_names do not exactly match AnnData obs_names; not creating mdata_f view.")
    except Exception as e:
        print("Could not create mdata_f view; heatmap does not require it.")
        print(type(e).__name__, e)

# %% [markdown]
# ## Load projected cNMF scores and zero-based program names

# %%
if SCORE_KEY not in adata.obsm:
    raise KeyError(f"{SCORE_KEY!r} is not in adata.obsm. Available obsm keys: {list(adata.obsm.keys())}")

usage = np.asarray(adata.obsm[SCORE_KEY])
if usage.ndim != 2:
    raise ValueError(f"Expected {SCORE_KEY} to be 2D, got shape {usage.shape}")

n_cells, n_programs = usage.shape
print("usage shape:", usage.shape)
assert n_cells == adata.n_obs

# Program names from uns, with robust zero-based fallback.
projection_uns = adata.uns.get("perturbseq_cNMF_projection", {}) if hasattr(adata, "uns") else {}
program_names = None
if isinstance(projection_uns, dict) and "program_names" in projection_uns:
    program_names = [str(x) for x in list(projection_uns["program_names"])]

if program_names is None or len(program_names) != n_programs:
    score_cols = [c for c in obs.columns if c.startswith("perturbseq_cNMF_P") and c.endswith("_usage")]
    if len(score_cols) == n_programs:
        def score_col_to_program(c):
            return c.replace("perturbseq_cNMF_", "").replace("_usage", "")
        score_cols = sorted(score_cols, key=lambda c: int(score_col_to_program(c).replace("P", "")))
        program_names = [score_col_to_program(c) for c in score_cols]
    else:
        program_names = [f"P{i}" for i in range(n_programs)]

expected_zero = [f"P{i}" for i in range(n_programs)]
expected_one = [f"P{i}" for i in range(1, n_programs + 1)]
if program_names == expected_one:
    print("WARNING: program names look one-based. Relabeling display names to zero-based P0...P49.")
    program_names = expected_zero
elif program_names != expected_zero:
    print("WARNING: program names are not exactly P0...P49. Using names from file:")
    print(program_names)

print("Program names:", program_names[:10], "...", program_names[-5:])

usage_f = usage[keep_mask, :]
obs_f = obs.loc[keep_mask].copy()
print("Filtered usage shape:", usage_f.shape)

# Basic sanity checks.
print("NaN count:", np.isnan(usage_f).sum())
print("Inf count:", np.isinf(usage_f).sum())
print("negative count:", (usage_f < -1e-8).sum())
print("row sum summary:")
print(pd.Series(usage_f.sum(axis=1)).describe())

# %% [markdown]
# ## Attach program labels from program metadata

# %%
label_map = {}
if PROGRAM_META_PATH.exists():
    prog_meta = pd.read_csv(PROGRAM_META_PATH)
    if {"program_name", "program_label"}.issubset(prog_meta.columns):
        tmp = prog_meta[["program_name", "program_label"]].drop_duplicates("program_name").copy()
        tmp["program_name"] = tmp["program_name"].astype(int)
        label_map = {f"P{int(r.program_name)}": str(r.program_label) for r in tmp.itertuples(index=False)}
        print(f"Loaded program labels from {PROGRAM_META_PATH}")
    else:
        print("Program metadata missing program_name/program_label columns; using P labels only.")
else:
    print("Program metadata file not found; using P labels only:", PROGRAM_META_PATH)


def shorten_label(x, max_chars=MAX_PROGRAM_LABEL_CHARS):
    x = str(x)
    if len(x) <= max_chars:
        return x
    return x[: max_chars - 1] + "…"

program_display_names = []
for p in program_names:
    label = label_map.get(p, "")
    if SHOW_PROGRAM_LABELS and label and label.lower() != "nan":
        program_display_names.append(f"{p} | {shorten_label(label)}")
    else:
        program_display_names.append(p)

program_label_df = pd.DataFrame({
    "program": program_names,
    "display_name": program_display_names,
    "program_label": [label_map.get(p, "") for p in program_names],
})
show_df(program_label_df, n=100)
program_label_df.to_csv(PLOT_OUT_DIR / "zero_based_program_labels_used_in_heatmap.tsv", sep="\t", index=False)

# %% [markdown]
# ## Compute mean projected cNMF usage by sample_name

# %%
if SAMPLE_COL not in obs_f.columns:
    raise KeyError(f"{SAMPLE_COL!r} not found in filtered obs.")

score_df = pd.DataFrame(usage_f, index=obs_f.index, columns=program_names)
score_df[SAMPLE_COL] = obs_f[SAMPLE_COL].astype(str).to_numpy()
score_df[CONDITION_COL] = obs_f[CONDITION_COL].astype(str).to_numpy()

group_mean_sample_x_program = score_df.groupby(SAMPLE_COL, observed=True)[program_names].mean()
group_median_sample_x_program = score_df.groupby(SAMPLE_COL, observed=True)[program_names].median()
group_n = score_df.groupby(SAMPLE_COL, observed=True).size().rename("n_cells")

# Order samples by condition order, then sample_name.
def mode_or_first(values):
    vc = pd.Series(values).value_counts(dropna=False)
    return str(vc.index[0]) if len(vc) else ""

sample_condition = score_df.groupby(SAMPLE_COL, observed=True)[CONDITION_COL].agg(mode_or_first)
sample_order_df = pd.DataFrame({
    SAMPLE_COL: group_mean_sample_x_program.index.astype(str),
})
sample_order_df[CONDITION_COL] = sample_order_df[SAMPLE_COL].map(sample_condition)
sample_order_df["condition_rank"] = sample_order_df[CONDITION_COL].map({c: i for i, c in enumerate(condition_order)})
sample_order_df["condition_rank"] = sample_order_df["condition_rank"].fillna(len(condition_order)).astype(int)
sample_order_df = sample_order_df.sort_values(["condition_rank", SAMPLE_COL])
sample_order = sample_order_df[SAMPLE_COL].tolist()

heatmap_mean = group_mean_sample_x_program.loc[sample_order, program_names].T
heatmap_median = group_median_sample_x_program.loc[sample_order, program_names].T
sample_n_ordered = group_n.loc[sample_order]

print("Heatmap mean matrix, programs x sample_name:", heatmap_mean.shape)
show_df(sample_order_df, n=100)
show_df(heatmap_mean.iloc[:10, :10], n=20)

heatmap_mean.to_csv(PLOT_OUT_DIR / "projected_cNMF_mean_usage_by_sample_name.filtered.tsv", sep="\t")
heatmap_median.to_csv(PLOT_OUT_DIR / "projected_cNMF_median_usage_by_sample_name.filtered.tsv", sep="\t")
sample_n_ordered.to_csv(PLOT_OUT_DIR / "filtered_cells_per_sample_name.tsv", sep="\t", header=True)

# %% [markdown]
# ## Make heatmap

# %%
if HEATMAP_SCALE == "zscore_by_program":
    row_mean = heatmap_mean.mean(axis=1)
    row_std = heatmap_mean.std(axis=1).replace(0, np.nan)
    heatmap_plot = heatmap_mean.sub(row_mean, axis=0).div(row_std, axis=0).fillna(0)
    colorbar_label = "Mean projected cNMF usage, z-scored within program"
    file_suffix = "mean_usage_zscore_by_program"
elif HEATMAP_SCALE == "mean":
    heatmap_plot = heatmap_mean.copy()
    colorbar_label = "Mean projected cNMF usage"
    file_suffix = "mean_usage"
else:
    raise ValueError(f"Unknown HEATMAP_SCALE: {HEATMAP_SCALE}")

heatmap_plot.to_csv(PLOT_OUT_DIR / f"projected_cNMF_by_sample_name.filtered.{file_suffix}.tsv", sep="\t")

n_programs, n_samples = heatmap_plot.shape
fig_w = max(8, 0.42 * n_samples + 2.5)
fig_h = max(10, 0.24 * n_programs + 1.5)

xlabels = [f"{s}\n(n={int(sample_n_ordered.loc[s]):,})" for s in heatmap_plot.columns]

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
im = ax.imshow(heatmap_plot.to_numpy(), aspect="auto", interpolation="nearest")

ax.set_xticks(np.arange(n_samples))
ax.set_xticklabels(xlabels, rotation=90, ha="center", va="top", fontsize=6)
ax.set_yticks(np.arange(n_programs))
ax.set_yticklabels(program_display_names, fontsize=5)

ax.set_xlabel(SAMPLE_COL)
ax.set_ylabel("Projected perturb-seq cNMF program")
ax.set_title(f"Projected perturb-seq cNMF programs by {SAMPLE_COL}\nQC-filtered cells only")

cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
cbar.set_label(colorbar_label)

fig.tight_layout()

base = PLOT_OUT_DIR / f"heatmap_projected_cNMF_by_{SAMPLE_COL}.filtered.{file_suffix}.all_program_names"
fig.savefig(base.with_suffix(".png"), bbox_inches="tight")
fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
print("Saved:")
print(base.with_suffix(".png"))
print(base.with_suffix(".pdf"))
print(base.with_suffix(".svg"))
plt.show()

# %% [markdown]
# ## Optional: show the unscaled mean usage heatmap too

# %%
if HEATMAP_SCALE != "mean":
    heatmap_plot_unscaled = heatmap_mean.copy()
    heatmap_plot_unscaled.to_csv(PLOT_OUT_DIR / "projected_cNMF_by_sample_name.filtered.mean_usage.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(heatmap_plot_unscaled.to_numpy(), aspect="auto", interpolation="nearest")

    ax.set_xticks(np.arange(n_samples))
    ax.set_xticklabels(xlabels, rotation=90, ha="center", va="top", fontsize=6)
    ax.set_yticks(np.arange(n_programs))
    ax.set_yticklabels(program_display_names, fontsize=5)

    ax.set_xlabel(SAMPLE_COL)
    ax.set_ylabel("Projected perturb-seq cNMF program")
    ax.set_title(f"Projected perturb-seq cNMF programs by {SAMPLE_COL}\nQC-filtered cells only; unscaled mean usage")

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
    cbar.set_label("Mean projected cNMF usage")

    fig.tight_layout()

    base_unscaled = PLOT_OUT_DIR / f"heatmap_projected_cNMF_by_{SAMPLE_COL}.filtered.mean_usage.all_program_names"
    fig.savefig(base_unscaled.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(base_unscaled.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base_unscaled.with_suffix(".svg"), bbox_inches="tight")
    print("Saved unscaled heatmap:")
    print(base_unscaled.with_suffix(".png"))
    print(base_unscaled.with_suffix(".pdf"))
    print(base_unscaled.with_suffix(".svg"))
    plt.show()

# %% [markdown]
# ## Close backed files

# %%
try:
    adata.file.close()
except Exception:
    pass

if mdata is not None:
    try:
        mdata.file.close()
    except Exception:
        pass

print("Done.")
