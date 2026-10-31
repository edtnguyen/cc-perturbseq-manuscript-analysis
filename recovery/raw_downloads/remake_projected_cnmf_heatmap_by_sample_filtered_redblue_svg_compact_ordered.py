# %%
from pathlib import Path
import re
import textwrap
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt

# Keep text editable in Illustrator SVG output
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# %%
OUT_DIR = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_projection")
H5AD_PATH = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad"
PROGRAM_META_PATH = Path("/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/data/program_activity_long_v4.csv")
HEATMAP_OUT_DIR = OUT_DIR / "heatmaps_after_qc_filter_redblue_svg_compact_ordered"
HEATMAP_OUT_DIR.mkdir(parents=True, exist_ok=True)

SCORE_KEY = "X_perturbseq_cNMF_usage"  # normalized projected usage
SAMPLE_COL = "sample_name"
CONDITION_COL = "condition"

# More compact label options
PROGRAM_LABEL_MAX_CHARS = 32
PROGRAM_LABEL_WRAP = 24  # set None to disable wrapping
PROGRAM_LABEL_MODE = "compact"  # compact / full / id_only

# Wider canvas options for long sample names
BASE_FIG_W = 16
PER_SAMPLE_W = 0.9
BASE_FIG_H = 14
PER_PROGRAM_H = 0.28
X_LABEL_ROTATION = 50
X_LABEL_HA = "right"

# Red-blue diverging heatmap palette
ZSCORE_CMAP = "RdBu_r"
RAW_CMAP = "RdBu_r"

# Exact time-course order at the condition level
CONDITION_TIMECOURSE_ORDER = [
    'D0_hPSC',
    'D1_MPS',
    'D1_MPS+TDI',
    'D1_MPS+IAG',
    'D1.33_LM',
    'D1.33_LM+TDI',
    'D1.33_LM+IAG',
    'D1.67_LM',
    'D1.67_LM+TDI',
    'D1.67_LM+IAG',
    'D2_LM',
    'D2_LM+TDI',
    'D2_LM+IAG',
    'D2_LM+TDI->D3_AE+TDI',
    'D2_LM+IAG->D3_AE+IAG',
    'D3_AE',
    'D3_AE+TDI',
    'D3_AE+IAG',
]

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
def infer_condition_from_sample_name(sample_name: str):
    sample_name = str(sample_name)
    # Exact match first
    if sample_name in cond_pct_mt_cutofs:
        return sample_name
    # Otherwise find the longest condition key contained in sample_name
    hits = [cond for cond in cond_pct_mt_cutofs if cond in sample_name]
    if not hits:
        return None
    return sorted(hits, key=len, reverse=True)[0]


def compact_program_label(program_id: str, raw_label: str):
    if PROGRAM_LABEL_MODE == "id_only":
        label = str(program_id)
    elif PROGRAM_LABEL_MODE == "full":
        label = f"{program_id}: {raw_label}"
    else:
        # compact
        label = f"{program_id} {raw_label}"
    label = re.sub(r"\s+", " ", str(label)).strip()
    if PROGRAM_LABEL_MAX_CHARS is not None and len(label) > PROGRAM_LABEL_MAX_CHARS:
        label = label[: PROGRAM_LABEL_MAX_CHARS - 1].rstrip() + "…"
    if PROGRAM_LABEL_WRAP:
        label = "\n".join(textwrap.wrap(label, width=PROGRAM_LABEL_WRAP, break_long_words=False))
    return label


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
    inferred = adata.obs[SAMPLE_COL].astype(str).map(infer_condition_from_sample_name)
    if inferred.isna().any():
        bad = adata.obs.loc[inferred.isna(), SAMPLE_COL].astype(str).unique().tolist()
        raise ValueError(
            f"Could not infer condition from some sample names. Example unmatched sample_name values: {bad[:10]}"
        )
    adata.obs[CONDITION_COL] = inferred.to_numpy()
    print(f"`{CONDITION_COL}` not found; inferred from `{SAMPLE_COL}`.")

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
cond_series = adata.obs[CONDITION_COL].astype(str)
sample_series = adata.obs[SAMPLE_COL].astype(str)

for cond in cond_series.value_counts().index:
    cond_mask = cond_series.to_numpy() == cond
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

    for sample_name in sample_series[cond_mask].value_counts().index:
        sub_mask = cond_mask & (sample_series.to_numpy() == sample_name)
        sub_fail = fail_mask & (sample_series.to_numpy() == sample_name)
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

# %%
# Load human-readable program labels if available
program_display_map = {p: p for p in program_names}
program_label_full_map = {p: p for p in program_names}
if PROGRAM_META_PATH.exists():
    prog_meta = pd.read_csv(PROGRAM_META_PATH).drop_duplicates("program_name")
    if "program_name" in prog_meta.columns:
        for _, row in prog_meta.iterrows():
            pid = int(row["program_name"])
            pname = f"P{pid}"
            if pname in program_display_map:
                if "program_label" in prog_meta.columns and pd.notna(row["program_label"]):
                    raw_label = str(row["program_label"])
                else:
                    raw_label = pname
                program_label_full_map[pname] = raw_label
                program_display_map[pname] = compact_program_label(pname, raw_label)

label_table = pd.DataFrame({
    "program": program_names,
    "full_program_label": [program_label_full_map[p] for p in program_names],
    "display_label": [program_display_map[p] for p in program_names],
})
label_table.to_csv(HEATMAP_OUT_DIR / "zero_based_program_labels_used_in_heatmap.tsv", sep="\t", index=False)

# %%
# Make per-cell score dataframe and aggregate by sample_name
if SAMPLE_COL not in adata_f.obs.columns:
    raise KeyError(f"{SAMPLE_COL} not found in adata.obs")

score_df = pd.DataFrame(usage, index=adata_f.obs_names, columns=program_names)
score_df[SAMPLE_COL] = adata_f.obs[SAMPLE_COL].astype(str).to_numpy()
score_df[CONDITION_COL] = adata_f.obs[CONDITION_COL].astype(str).to_numpy()

# Order sample columns exactly by time course using condition order, then preserve within-condition sample order as observed
observed_pairs = score_df[[CONDITION_COL, SAMPLE_COL]].drop_duplicates()
observed_sample_order = []
for cond in CONDITION_TIMECOURSE_ORDER:
    sub = observed_pairs.loc[observed_pairs[CONDITION_COL] == cond, SAMPLE_COL].tolist()
    observed_sample_order.extend(sub)
# Append any unexpected conditions/samples afterward in first-observed order
remaining = [s for s in observed_pairs[SAMPLE_COL].tolist() if s not in observed_sample_order]
sample_order = observed_sample_order + remaining

sample_order_df = pd.DataFrame({
    SAMPLE_COL: sample_order,
    CONDITION_COL: [
        score_df.loc[score_df[SAMPLE_COL] == s, CONDITION_COL].iloc[0] for s in sample_order
    ],
    "sample_order_index": np.arange(len(sample_order)),
})
sample_order_df.to_csv(HEATMAP_OUT_DIR / f"{SAMPLE_COL}_timecourse_order_used.tsv", sep="\t", index=False)

# Aggregate
mean_by_sample = score_df.groupby(SAMPLE_COL)[program_names].mean().reindex(sample_order)
median_by_sample = score_df.groupby(SAMPLE_COL)[program_names].median().reindex(sample_order)
n_by_sample = score_df.groupby(SAMPLE_COL).size().rename("n_cells").reindex(sample_order)

mean_by_sample.to_csv(HEATMAP_OUT_DIR / f"projected_cNMF_mean_usage_by_{SAMPLE_COL}.filtered.tsv", sep="\t")
median_by_sample.to_csv(HEATMAP_OUT_DIR / f"projected_cNMF_median_usage_by_{SAMPLE_COL}.filtered.tsv", sep="\t")
n_by_sample.to_csv(HEATMAP_OUT_DIR / f"projected_cNMF_ncells_by_{SAMPLE_COL}.filtered.tsv", sep="\t")

print(mean_by_sample.shape)
display(sample_order_df)
display(n_by_sample.to_frame())
display(mean_by_sample.head())

# %%
# Build heatmap matrices
heatmap_raw = mean_by_sample.T.copy()  # programs x sample_name
heatmap_raw.index = [program_display_map[p] for p in heatmap_raw.index]

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

    fig_w = max(BASE_FIG_W, PER_SAMPLE_W * n_cols + 8)
    fig_h = max(BASE_FIG_H, PER_PROGRAM_H * n_rows + 4)

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
    ax.set_xticklabels(df.columns.tolist(), rotation=X_LABEL_ROTATION, ha=X_LABEL_HA)
    ax.set_yticklabels(df.index.tolist(), fontsize=9)
    ax.set_xlabel(SAMPLE_COL)
    ax.set_ylabel("Projected perturb-seq cNMF program")
    ax.set_title(title)

    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label("z-score across samples" if center_zero else "mean projected usage")

    fig.subplots_adjust(left=0.28, bottom=0.22, right=0.96, top=0.96)

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
    out_prefix=f"heatmap_projected_cNMF_by_{SAMPLE_COL}.filtered.mean_usage_zscore_by_program.all_program_names.redblue_editable.compact_ordered",
    cmap=ZSCORE_CMAP,
    center_zero=True,
)

# %%
# Raw mean-usage red-blue heatmap
plot_heatmap_pcolormesh(
    heatmap_raw,
    title=f"Projected perturb-seq cNMF by {SAMPLE_COL} (filtered; raw mean usage)",
    out_prefix=f"heatmap_projected_cNMF_by_{SAMPLE_COL}.filtered.mean_usage.all_program_names.redblue_editable.compact_ordered",
    cmap=RAW_CMAP,
    center_zero=False,
)

# %%
print("Saved outputs to:")
print(HEATMAP_OUT_DIR)
for p in sorted(HEATMAP_OUT_DIR.glob("*")):
    print(p.name)
