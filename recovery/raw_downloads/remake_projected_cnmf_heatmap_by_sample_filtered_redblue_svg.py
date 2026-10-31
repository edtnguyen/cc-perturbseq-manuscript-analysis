# %%
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
import matplotlib as mpl
import matplotlib.pyplot as plt

# Keep text editable in Illustrator SVG output
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# %%
OUT_DIR = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_projection")
H5AD_PATH = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad"
H5MU_PATH = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5mu"
PROGRAM_META_PATH = Path("/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/data/program_activity_long_v4.csv")
HEATMAP_OUT_DIR = OUT_DIR / "heatmaps_after_qc_filter_redblue_svg"
HEATMAP_OUT_DIR.mkdir(parents=True, exist_ok=True)

SCORE_KEY = "X_perturbseq_cNMF_usage"  # normalized projected usage
SAMPLE_COL = "sample_name"
CONDITION_COL = "condition"

# Optional explicit sample order; otherwise uses observed order in the filtered data.
SAMPLE_ORDER = None

# Red-blue diverging heatmap palette
ZSCORE_CMAP = "RdBu_r"
RAW_CMAP = "RdBu_r"

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
    'D2_LM+IAG->D3_AE+IAG': 1.0
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
    'D2_LM+IAG->D3_AE+IAG': 3.75
}

# %%
adata = sc.read_h5ad(H5AD_PATH)
print(adata)
print("obs columns:", list(adata.obs.columns))
print("obsm keys:", list(adata.obsm.keys()))

# %%
# Derive or use condition column
if CONDITION_COL not in adata.obs.columns:
    if SAMPLE_COL not in adata.obs.columns:
        raise ValueError(f"Neither `{CONDITION_COL}` nor `{SAMPLE_COL}` found in adata.obs")
    adata.obs[CONDITION_COL] = adata.obs[SAMPLE_COL].astype(str)
    print(f"`{CONDITION_COL}` not found; copied from `{SAMPLE_COL}`.")

# Validate condition coverage for filtering dictionaries
conditions_present = pd.Index(adata.obs[CONDITION_COL].astype(str).unique())
missing_mt = sorted(set(conditions_present) - set(cond_pct_mt_cutofs))
missing_ribo = sorted(set(conditions_present) - set(cond_pct_ribo_cutofs))
if missing_mt or missing_ribo:
    raise ValueError(
        "Some conditions are missing from the QC cutoff dictionaries.\n"
        f"Missing in cond_pct_mt_cutofs: {missing_mt}\n"
        f"Missing in cond_pct_ribo_cutofs: {missing_ribo}"
    )

# %%
# Apply the requested condition-specific QC filters
keep_mask = np.ones(adata.n_obs, dtype=bool)
summary_records = []
summary_sample_records = []

for cond in adata.obs[CONDITION_COL].astype(str).value_counts().index:
    cond_mask = adata.obs[CONDITION_COL].astype(str).to_numpy() == cond
    fail_mask = cond_mask & (
        (adata.obs["pct_counts_mt"].to_numpy() >= cond_pct_mt_cutofs[cond]) |
        (adata.obs["pct_counts_ribo"].to_numpy() <= cond_pct_ribo_cutofs[cond])
    )
    keep_mask[fail_mask] = False

    n_total = int(cond_mask.sum())
    n_remove = int(fail_mask.sum())
    n_keep = n_total - n_remove
    summary_records.append({
        "condition": cond,
        "n_total": n_total,
        "n_remove": n_remove,
        "n_keep": n_keep,
        "frac_remove": n_remove / n_total if n_total > 0 else np.nan,
        "pct_counts_mt_cutoff": cond_pct_mt_cutofs[cond],
        "pct_counts_ribo_cutoff": cond_pct_ribo_cutofs[cond],
    })

    if SAMPLE_COL in adata.obs.columns:
        for sample_name, sample_idx in adata.obs.loc[cond_mask, SAMPLE_COL].astype(str).value_counts().items():
            sub_mask = cond_mask & (adata.obs[SAMPLE_COL].astype(str).to_numpy() == sample_name)
            sub_fail = fail_mask & (adata.obs[SAMPLE_COL].astype(str).to_numpy() == sample_name)
            n_sub_total = int(sub_mask.sum())
            n_sub_remove = int(sub_fail.sum())
            summary_sample_records.append({
                "condition": cond,
                SAMPLE_COL: sample_name,
                "n_total": n_sub_total,
                "n_remove": n_sub_remove,
                "n_keep": n_sub_total - n_sub_remove,
                "frac_remove": n_sub_remove / n_sub_total if n_sub_total > 0 else np.nan,
            })

adata_f = adata[keep_mask].copy()
print(f"Cells before filter: {adata.n_obs:,}")
print(f"Cells after filter:  {adata_f.n_obs:,}")

qc_summary_df = pd.DataFrame(summary_records)
qc_summary_sample_df = pd.DataFrame(summary_sample_records)
qc_summary_df.to_csv(HEATMAP_OUT_DIR / "qc_filter_summary_by_condition.tsv", sep="\t", index=False)
qc_summary_sample_df.to_csv(HEATMAP_OUT_DIR / f"qc_filter_summary_by_condition_and_{SAMPLE_COL}.tsv", sep="\t", index=False)

# %%
# Pull projected usage matrix after filtering
if SCORE_KEY not in adata_f.obsm:
    raise KeyError(f"{SCORE_KEY} not found in adata.obsm")

usage = np.asarray(adata_f.obsm[SCORE_KEY])
print("usage shape:", usage.shape)

# Program names from saved projection metadata if available; otherwise infer P0..P49
program_names = None
if "perturbseq_cNMF_projection" in adata_f.uns:
    meta_uns = adata_f.uns["perturbseq_cNMF_projection"]
    if isinstance(meta_uns, dict) and "program_names" in meta_uns:
        program_names = [str(x) for x in meta_uns["program_names"]]

if program_names is None:
    program_names = [f"P{i}" for i in range(usage.shape[1])]

if len(program_names) != usage.shape[1]:
    raise ValueError(
        f"Number of program names ({len(program_names)}) does not match usage columns ({usage.shape[1]})."
    )

print("First program names:", program_names[:5])
print("Last program names:", program_names[-5:])

# %%
# Load human-readable program labels if available
program_label_map = {p: p for p in program_names}
label_table = None
if PROGRAM_META_PATH.exists():
    prog_meta = pd.read_csv(PROGRAM_META_PATH).drop_duplicates("program_name")
    if "program_name" in prog_meta.columns:
        for _, row in prog_meta.iterrows():
            pid = int(row["program_name"])
            pname = f"P{pid}"
            if pname in program_label_map:
                if "program_label" in prog_meta.columns and pd.notna(row["program_label"]):
                    program_label_map[pname] = f"{pname}: {row['program_label']}"
                else:
                    program_label_map[pname] = pname

label_table = pd.DataFrame({
    "program": program_names,
    "display_label": [program_label_map[p] for p in program_names],
})
label_table.to_csv(HEATMAP_OUT_DIR / "zero_based_program_labels_used_in_heatmap.tsv", sep="\t", index=False)

# %%
# Make per-cell score dataframe and aggregate by sample_name
if SAMPLE_COL not in adata_f.obs.columns:
    raise KeyError(f"{SAMPLE_COL} not found in adata.obs")

score_df = pd.DataFrame(usage, index=adata_f.obs_names, columns=program_names)
score_df[SAMPLE_COL] = adata_f.obs[SAMPLE_COL].astype(str).to_numpy()
score_df[CONDITION_COL] = adata_f.obs[CONDITION_COL].astype(str).to_numpy()

if SAMPLE_ORDER is None:
    sample_order = score_df[SAMPLE_COL].drop_duplicates().tolist()
else:
    sample_order = [x for x in SAMPLE_ORDER if x in set(score_df[SAMPLE_COL])]

# Aggregate
mean_by_sample = score_df.groupby(SAMPLE_COL)[program_names].mean()
median_by_sample = score_df.groupby(SAMPLE_COL)[program_names].median()
n_by_sample = score_df.groupby(SAMPLE_COL).size().rename("n_cells")

# Reorder samples
mean_by_sample = mean_by_sample.reindex(sample_order)
median_by_sample = median_by_sample.reindex(sample_order)
n_by_sample = n_by_sample.reindex(sample_order)

mean_by_sample.to_csv(HEATMAP_OUT_DIR / f"projected_cNMF_mean_usage_by_{SAMPLE_COL}.filtered.tsv", sep="\t")
median_by_sample.to_csv(HEATMAP_OUT_DIR / f"projected_cNMF_median_usage_by_{SAMPLE_COL}.filtered.tsv", sep="\t")
n_by_sample.to_csv(HEATMAP_OUT_DIR / f"projected_cNMF_ncells_by_{SAMPLE_COL}.filtered.tsv", sep="\t")

print(mean_by_sample.shape)
display(n_by_sample.to_frame())
display(mean_by_sample.head())

# %%
# Build heatmap matrices
heatmap_raw = mean_by_sample.T.copy()  # programs x sample_name
heatmap_raw.index = [program_label_map[p] for p in heatmap_raw.index]

heatmap_z = heatmap_raw.sub(heatmap_raw.mean(axis=1), axis=0)
heatmap_z = heatmap_z.div(heatmap_raw.std(axis=1).replace(0, np.nan), axis=0)

heatmap_raw.to_csv(HEATMAP_OUT_DIR / f"projected_cNMF_mean_usage_by_{SAMPLE_COL}.filtered.programs_x_samples.tsv", sep="\t")
heatmap_z.to_csv(HEATMAP_OUT_DIR / f"projected_cNMF_mean_usage_zscore_by_program_by_{SAMPLE_COL}.filtered.programs_x_samples.tsv", sep="\t")

# %%
def plot_heatmap_pcolormesh(df, title, out_prefix, cmap="RdBu_r", center_zero=False):
    mat = df.to_numpy(dtype=float)
    n_rows, n_cols = mat.shape

    # Symmetric color range for z-scored heatmaps; regular min-max for raw
    if center_zero:
        vmax = np.nanmax(np.abs(mat))
        vmin = -vmax
    else:
        vmin = np.nanmin(mat)
        vmax = np.nanmax(mat)
        if np.isclose(vmin, vmax):
            vmax = vmin + 1e-9

    fig_w = max(10, 0.55 * n_cols)
    fig_h = max(12, 0.33 * n_rows)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    mesh = ax.pcolormesh(
        np.arange(n_cols + 1),
        np.arange(n_rows + 1),
        mat,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="flat",
        linewidth=0,
        antialiased=False,
    )

    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)
    ax.set_xticks(np.arange(n_cols) + 0.5)
    ax.set_yticks(np.arange(n_rows) + 0.5)
    ax.set_xticklabels(df.columns.tolist(), rotation=90)
    ax.set_yticklabels(df.index.tolist())
    ax.set_xlabel(SAMPLE_COL)
    ax.set_ylabel("Projected perturb-seq cNMF program")
    ax.set_title(title)

    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label("z-score across samples" if center_zero else "mean projected usage")

    fig.tight_layout()

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(
            HEATMAP_OUT_DIR / f"{out_prefix}.{ext}",
            dpi=300 if ext == "png" else None,
            bbox_inches="tight",
        )

    plt.show()
    plt.close(fig)

# %%
# Main z-scored red-blue heatmap
plot_heatmap_pcolormesh(
    heatmap_z,
    title=f"Projected perturb-seq cNMF by {SAMPLE_COL} (filtered; z-scored by program)",
    out_prefix=f"heatmap_projected_cNMF_by_{SAMPLE_COL}.filtered.mean_usage_zscore_by_program.all_program_names.redblue_editable",
    cmap=ZSCORE_CMAP,
    center_zero=True,
)

# %%
# Raw mean-usage red-blue heatmap
plot_heatmap_pcolormesh(
    heatmap_raw,
    title=f"Projected perturb-seq cNMF by {SAMPLE_COL} (filtered; raw mean usage)",
    out_prefix=f"heatmap_projected_cNMF_by_{SAMPLE_COL}.filtered.mean_usage.all_program_names.redblue_editable",
    cmap=RAW_CMAP,
    center_zero=False,
)

# %%
print("Saved outputs to:")
print(HEATMAP_OUT_DIR)
for p in sorted(HEATMAP_OUT_DIR.glob("*")):
    print(p.name)
