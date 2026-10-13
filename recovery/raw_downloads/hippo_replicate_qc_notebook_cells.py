# %% [markdown]
# # Replicate QC for Hippo time-course projected cNMF and pseudobulk DE
#
# These cells help identify a biological replicate or sample that behaves differently from the others.
# The first half uses the projected cNMF program-score tables only.
# The optional second half uses pseudobulk counts if those files are available on Sherlock.
#
# Main questions:
#
# 1. For key programs P15/P17/P34, does one matched replicate go opposite to the other two?
# 2. Across all 50 programs, does one replicate have a low correlation with the other matched replicates?
# 3. Does one sample look far from the other samples in the same condition?
# 4. If a replicate looks suspicious, how much do conclusions change when excluding it?

# %%
from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd
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
    "svg.fonttype": "none",
})

# Optional scipy for t-tests and PCA convenience is avoided where possible.
try:
    from scipy import stats
except Exception:
    stats = None

# %% [markdown]
# ## Settings
#
# On Sherlock, use the default `BASE_DIR`. If you are running this notebook in a directory containing the uploaded TSVs, set `BASE_DIR = Path('.')`.

# %%
BASE_DIR = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_matched_DE")

# If the Sherlock path is not available, fall back to the current directory and /mnt/data for local review.
FALLBACK_DIRS = [Path.cwd(), Path("/mnt/data")]

OUT_DIR = BASE_DIR / "replicate_qc"
if not BASE_DIR.exists():
    OUT_DIR = Path.cwd() / "replicate_qc"
OUT_DIR.mkdir(parents=True, exist_ok=True)

KEY_PROGRAMS = ["P15", "P17", "P34"]

CONTRAST_ORDER = [
    "TDI_at_D1_MPS",
    "IAG_at_D1_MPS",
    "TDI_at_D1.33_LM",
    "IAG_at_D1.33_LM",
    "TDI_at_D1.67_LM",
    "IAG_at_D1.67_LM",
    "TDI_at_D2_LM",
    "IAG_at_D2_LM",
    "TDI_at_D3_AE",
    "IAG_at_D3_AE",
    "TDI_D2_to_D3_vs_D3_AE",
    "IAG_D2_to_D3_vs_D3_AE",
]

# These are optional. If present, the pseudobulk QC cells will run.
PSEUDOBULK_COUNTS_PATH = BASE_DIR / "pseudobulk_gene_counts.filtered.tsv.gz"
PSEUDOBULK_META_PATH = BASE_DIR / "pseudobulk_sample_metadata.filtered.tsv"

# Main program-score QC files.
PAIRED_DIFF_PATH = BASE_DIR / "program_score_DE" / "projected_cNMF_program_score_paired_replicate_differences.tsv"
MATCHED_SUMMARY_PATH = BASE_DIR / "program_score_DE" / "projected_cNMF_program_score_matched_contrasts.tsv"
KEY_SUMMARY_PATH = BASE_DIR / "program_score_DE" / "key_programs_P15_P17_P34_matched_contrasts.tsv"


def resolve_path(path, fallback_dirs=FALLBACK_DIRS):
    path = Path(path)
    if path.exists():
        return path
    for d in fallback_dirs:
        alt = d / path.name
        if alt.exists():
            return alt
    return path

PAIRED_DIFF_PATH = resolve_path(PAIRED_DIFF_PATH)
MATCHED_SUMMARY_PATH = resolve_path(MATCHED_SUMMARY_PATH)
KEY_SUMMARY_PATH = resolve_path(KEY_SUMMARY_PATH)
PSEUDOBULK_COUNTS_PATH = resolve_path(PSEUDOBULK_COUNTS_PATH)
PSEUDOBULK_META_PATH = resolve_path(PSEUDOBULK_META_PATH)

print("PAIRED_DIFF_PATH:", PAIRED_DIFF_PATH, PAIRED_DIFF_PATH.exists())
print("MATCHED_SUMMARY_PATH:", MATCHED_SUMMARY_PATH, MATCHED_SUMMARY_PATH.exists())
print("KEY_SUMMARY_PATH:", KEY_SUMMARY_PATH, KEY_SUMMARY_PATH.exists())
print("PSEUDOBULK_COUNTS_PATH:", PSEUDOBULK_COUNTS_PATH, PSEUDOBULK_COUNTS_PATH.exists())
print("PSEUDOBULK_META_PATH:", PSEUDOBULK_META_PATH, PSEUDOBULK_META_PATH.exists())
print("OUT_DIR:", OUT_DIR)

# %% [markdown]
# ## Helper functions

# %%
def natural_program_sort_key(p):
    m = re.search(r"(\d+)$", str(p))
    return int(m.group(1)) if m else 10**9


def parse_replicate(sample):
    m = re.search(r"_rep(\d+)$", str(sample))
    return f"rep{m.group(1)}" if m else "unknown"


def parse_condition(sample):
    return re.sub(r"_rep\d+$", "", str(sample))


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if valid.sum() == 0:
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    restored = np.empty_like(adj)
    restored[order] = adj
    out[valid] = restored
    return out


def simple_pca(X, n_components=2, center=True, scale=False):
    X = np.asarray(X, dtype=float)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    if center:
        X = X - X.mean(axis=0, keepdims=True)
    if scale:
        sd = X.std(axis=0, keepdims=True)
        sd[sd == 0] = 1
        X = X / sd
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    scores = U[:, :n_components] * S[:n_components]
    var = S**2
    explained = var / var.sum() if var.sum() > 0 else np.zeros_like(var)
    return scores, explained[:n_components]


def corrcoef_rows(df):
    mat = df.to_numpy(dtype=float)
    mat = np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0)
    return pd.DataFrame(np.corrcoef(mat), index=df.index, columns=df.index)


def plot_matrix(mat_df, title, out_prefix=None, symmetric=False, cmap="RdBu_r", figsize=None):
    mat = mat_df.to_numpy(dtype=float)
    if symmetric:
        vmax = np.nanmax(np.abs(mat))
        vmin = -vmax
    else:
        vmin = np.nanmin(mat)
        vmax = np.nanmax(mat)
    if figsize is None:
        figsize = (max(5, 0.25 * mat_df.shape[1] + 2), max(4, 0.25 * mat_df.shape[0] + 2))
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(mat_df.shape[1]))
    ax.set_xticklabels(mat_df.columns, rotation=90)
    ax.set_yticks(np.arange(mat_df.shape[0]))
    ax.set_yticklabels(mat_df.index)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.tight_layout()
    if out_prefix is not None:
        for ext in ["png", "pdf", "svg"]:
            fig.savefig(Path(str(out_prefix) + f".{ext}"))
    plt.show()

# %% [markdown]
# ## Load program-score replicate tables

# %%
paired = pd.read_csv(PAIRED_DIFF_PATH, sep="\t")
summary = pd.read_csv(MATCHED_SUMMARY_PATH, sep="\t") if MATCHED_SUMMARY_PATH.exists() else None
key_summary = pd.read_csv(KEY_SUMMARY_PATH, sep="\t") if KEY_SUMMARY_PATH.exists() else None

paired["program"] = paired["program"].astype(str)
paired["replicate_de"] = paired["replicate_de"].astype(str)
paired["control_condition"] = paired["control_sample"].map(parse_condition)
paired["treated_condition"] = paired["treated_sample"].map(parse_condition)
paired["control_rep"] = paired["control_sample"].map(parse_replicate)
paired["treated_rep"] = paired["treated_sample"].map(parse_replicate)
paired["replicate_match_ok"] = paired["control_rep"] == paired["treated_rep"]
paired["paired_diff_check"] = paired["treated_value"] - paired["control_value"]
paired["paired_diff_abs_error"] = (paired["paired_diff_check"] - paired["paired_diff"]).abs()

programs = sorted(paired["program"].unique(), key=natural_program_sort_key)
contrasts = [c for c in CONTRAST_ORDER if c in set(paired["contrast_name"])] + [
    c for c in paired["contrast_name"].unique() if c not in CONTRAST_ORDER
]

print("paired shape:", paired.shape)
print("n contrasts:", len(contrasts), contrasts)
print("n programs:", len(programs), programs[:5], "...", programs[-5:])
print("replicates:", sorted(paired["replicate_de"].unique()))
print("all replicate matches ok:", paired["replicate_match_ok"].all())
print("max paired_diff numerical check error:", paired["paired_diff_abs_error"].max())

display(paired.head())

# %% [markdown]
# ## Key-program replicate plots: paired differences per replicate
#
# This view shows whether one replicate goes against the other two for P15/P17/P34.

# %%
key_df = paired[paired["program"].isin(KEY_PROGRAMS)].copy()
key_df["contrast_name"] = pd.Categorical(key_df["contrast_name"], categories=contrasts, ordered=True)
key_df["program"] = pd.Categorical(key_df["program"], categories=KEY_PROGRAMS, ordered=True)
key_df = key_df.sort_values(["program", "contrast_name", "replicate_de"])

fig, axes = plt.subplots(len(KEY_PROGRAMS), 1, figsize=(11, 2.5 * len(KEY_PROGRAMS)), sharex=True)
if len(KEY_PROGRAMS) == 1:
    axes = [axes]

x_positions = {c: i for i, c in enumerate(contrasts)}
rep_offsets = {"rep1": -0.18, "rep2": 0.0, "rep3": 0.18}

for ax, program in zip(axes, KEY_PROGRAMS):
    sub = key_df[key_df["program"] == program]
    ax.axhline(0, color="black", linewidth=0.8)
    for rep, rep_sub in sub.groupby("replicate_de"):
        xs = [x_positions[c] + rep_offsets.get(rep, 0) for c in rep_sub["contrast_name"].astype(str)]
        ax.scatter(xs, rep_sub["paired_diff"], s=35, label=rep if program == KEY_PROGRAMS[0] else None)
        for x, y in zip(xs, rep_sub["paired_diff"]):
            ax.plot([x, x], [0, y], linewidth=0.6, alpha=0.5)
    ax.set_ylabel(f"{program}\npaired diff")
    ax.set_title(f"{program}: treated - control per matched replicate")

axes[-1].set_xticks(range(len(contrasts)))
axes[-1].set_xticklabels(contrasts, rotation=45, ha="right")
axes[0].legend(frameon=False, ncol=3, loc="upper right")
fig.tight_layout()
for ext in ["png", "pdf", "svg"]:
    fig.savefig(OUT_DIR / f"key_programs_paired_differences_by_replicate.{ext}")
plt.show()

# %% [markdown]
# ## Quantify replicate disagreement per contrast/program
#
# Flags are conservative diagnostics, not automatic exclusion criteria.
#
# A replicate is flagged if:
#
# - its sign disagrees with the median/majority sign for that contrast-program, or
# - it is the strongest outlier from the replicate median, or
# - its leave-one-out influence on the mean effect is large.

# %%
qc_records = []

for (contrast, program), sub in paired.groupby(["contrast_name", "program"], sort=False):
    sub = sub.copy().sort_values("replicate_de")
    vals = sub["paired_diff"].to_numpy(dtype=float)
    reps = sub["replicate_de"].to_numpy()

    mean_val = np.nanmean(vals)
    median_val = np.nanmedian(vals)
    sd_val = np.nanstd(vals, ddof=1) if len(vals) > 1 else np.nan
    signs = np.sign(vals)
    nonzero_signs = signs[signs != 0]
    majority_sign = np.sign(np.nanmedian(vals)) if np.nanmedian(vals) != 0 else 0
    if majority_sign == 0 and len(nonzero_signs) > 0:
        majority_sign = pd.Series(nonzero_signs).mode().iloc[0]

    mad = np.nanmedian(np.abs(vals - median_val))
    denom = 1.4826 * mad if mad > 0 else (sd_val if np.isfinite(sd_val) and sd_val > 0 else np.nan)

    for i, row in enumerate(sub.itertuples(index=False)):
        val = float(row.paired_diff)
        loo_vals = np.delete(vals, i)
        loo_mean = np.nanmean(loo_vals) if len(loo_vals) else np.nan
        robust_z = (val - median_val) / denom if np.isfinite(denom) and denom > 0 else 0.0
        sign_discordant = (np.sign(val) != majority_sign) and (majority_sign != 0) and (val != 0)
        abs_dev = abs(val - median_val)
        influence = abs(mean_val - loo_mean) if np.isfinite(loo_mean) else np.nan

        qc_records.append({
            "contrast_name": contrast,
            "program": program,
            "replicate_de": row.replicate_de,
            "control_sample": row.control_sample,
            "treated_sample": row.treated_sample,
            "paired_diff": val,
            "mean_diff": mean_val,
            "median_diff": median_val,
            "sd_diff": sd_val,
            "majority_sign": majority_sign,
            "sign_discordant": bool(sign_discordant),
            "abs_deviation_from_median": abs_dev,
            "robust_z_like": robust_z,
            "leave_one_out_mean_without_this_rep": loo_mean,
            "absolute_influence_on_mean": influence,
        })

rep_qc = pd.DataFrame(qc_records)

# Rank outlierness within each contrast-program.
rep_qc["rank_abs_deviation_within_contrast_program"] = (
    rep_qc.groupby(["contrast_name", "program"])["abs_deviation_from_median"]
    .rank(method="first", ascending=False)
)
rep_qc["rank_influence_within_contrast_program"] = (
    rep_qc.groupby(["contrast_name", "program"])["absolute_influence_on_mean"]
    .rank(method="first", ascending=False)
)
rep_qc["is_largest_deviation_rep"] = rep_qc["rank_abs_deviation_within_contrast_program"] == 1

# With only 3 reps, robust z is limited. Treat this as a diagnostic.
rep_qc["flag"] = (
    rep_qc["sign_discordant"] |
    ((rep_qc["robust_z_like"].abs() >= 2.0) & rep_qc["is_largest_deviation_rep"])
)

rep_qc.to_csv(OUT_DIR / "replicate_program_score_disagreement_qc.tsv", sep="\t", index=False)

print("Flagged rows:", int(rep_qc["flag"].sum()))
display(rep_qc[rep_qc["flag"]].sort_values(["contrast_name", "program", "replicate_de"]).head(100))

# %%
# Summarize by replicate ID and by sample.
replicate_summary = (
    rep_qc.groupby("replicate_de")
    .agg(
        n_rows=("flag", "size"),
        n_flags=("flag", "sum"),
        n_sign_discordant=("sign_discordant", "sum"),
        median_abs_deviation=("abs_deviation_from_median", "median"),
        median_abs_influence=("absolute_influence_on_mean", "median"),
        max_abs_influence=("absolute_influence_on_mean", "max"),
    )
    .reset_index()
)
replicate_summary["flag_fraction"] = replicate_summary["n_flags"] / replicate_summary["n_rows"]
replicate_summary.to_csv(OUT_DIR / "replicate_id_summary_program_score_qc.tsv", sep="\t", index=False)
display(replicate_summary)

control_sample_summary = (
    rep_qc.groupby("control_sample")
    .agg(n_rows=("flag", "size"), n_flags=("flag", "sum"), median_abs_influence=("absolute_influence_on_mean", "median"))
    .reset_index()
)
treated_sample_summary = (
    rep_qc.groupby("treated_sample")
    .agg(n_rows=("flag", "size"), n_flags=("flag", "sum"), median_abs_influence=("absolute_influence_on_mean", "median"))
    .reset_index()
)
control_sample_summary.to_csv(OUT_DIR / "control_sample_summary_program_score_qc.tsv", sep="\t", index=False)
treated_sample_summary.to_csv(OUT_DIR / "treated_sample_summary_program_score_qc.tsv", sep="\t", index=False)

print("Control samples with most flags:")
display(control_sample_summary.sort_values(["n_flags", "median_abs_influence"], ascending=False).head(20))
print("Treated samples with most flags:")
display(treated_sample_summary.sort_values(["n_flags", "median_abs_influence"], ascending=False).head(20))

# %% [markdown]
# ## Contrast-wide replicate profile QC across all 50 programs
#
# Here each matched replicate is represented by a 50-dimensional vector of paired program-score differences.
# A suspicious replicate often has low correlation with the other two replicate-difference vectors in the same contrast.

# %%
profile = paired.pivot_table(
    index=["contrast_name", "replicate_de"],
    columns="program",
    values="paired_diff",
    aggfunc="mean",
)
profile = profile.reindex(columns=programs)
profile_index_df = profile.index.to_frame(index=False)
profile.index = [f"{c}|{r}" for c, r in profile.index]

profile.to_csv(OUT_DIR / "paired_diff_profiles_contrast_replicate_x_program.tsv", sep="\t")

# Correlation with other two replicates within the same contrast.
profile_qc_records = []
for contrast, sub_idx in profile_index_df.groupby("contrast_name").groups.items():
    labels = profile.index[list(sub_idx)]
    sub = profile.loc[labels]
    for label in labels:
        rep = label.split("|")[1]
        others = [x for x in labels if x != label]
        x = sub.loc[label].to_numpy(dtype=float)
        mean_other = sub.loc[others].mean(axis=0).to_numpy(dtype=float)
        corr = np.corrcoef(x, mean_other)[0, 1] if np.std(x) > 0 and np.std(mean_other) > 0 else np.nan
        dist = np.linalg.norm(x - mean_other)
        profile_qc_records.append({
            "contrast_name": contrast,
            "replicate_de": rep,
            "profile_label": label,
            "corr_to_mean_of_other_reps": corr,
            "euclidean_distance_to_mean_of_other_reps": dist,
        })

profile_qc = pd.DataFrame(profile_qc_records)
profile_qc["flag_low_corr"] = profile_qc["corr_to_mean_of_other_reps"] < 0.25
profile_qc.to_csv(OUT_DIR / "replicate_profile_qc_across_all_programs.tsv", sep="\t", index=False)

display(profile_qc.sort_values("corr_to_mean_of_other_reps").head(30))

# %%
# PCA of replicate-difference profiles.
coords, ev = simple_pca(profile, n_components=2, center=True, scale=True)
profile_plot_df = profile_index_df.copy()
profile_plot_df["profile_label"] = profile.index
profile_plot_df["PC1"] = coords[:, 0]
profile_plot_df["PC2"] = coords[:, 1]
profile_plot_df = profile_plot_df.merge(profile_qc, on=["contrast_name", "replicate_de", "profile_label"], how="left")

fig, ax = plt.subplots(figsize=(7, 5))
for rep, sub in profile_plot_df.groupby("replicate_de"):
    ax.scatter(sub["PC1"], sub["PC2"], s=45, label=rep)
    for _, row in sub.iterrows():
        if row.get("flag_low_corr", False):
            ax.text(row["PC1"], row["PC2"], row["profile_label"], fontsize=5)
ax.axhline(0, color="black", linewidth=0.5)
ax.axvline(0, color="black", linewidth=0.5)
ax.set_xlabel(f"PC1 ({ev[0]*100:.1f}% var)")
ax.set_ylabel(f"PC2 ({ev[1]*100:.1f}% var)")
ax.set_title("PCA of matched replicate-difference profiles across all programs")
ax.legend(frameon=False)
fig.tight_layout()
for ext in ["png", "pdf", "svg"]:
    fig.savefig(OUT_DIR / f"pca_replicate_difference_profiles.{ext}")
plt.show()

# %%
# Heatmap of replicate profile correlations. Useful for spotting profiles that do not cluster with same contrast.
profile_corr = corrcoef_rows(profile)
profile_corr.to_csv(OUT_DIR / "paired_diff_profile_correlation_matrix.tsv", sep="\t")
plot_matrix(
    profile_corr,
    title="Correlation of replicate-difference profiles across all programs",
    out_prefix=OUT_DIR / "paired_diff_profile_correlation_heatmap",
    symmetric=True,
    figsize=(10, 9),
)

# %% [markdown]
# ## Sample-level program profile QC
#
# Reconstruct one projected-program profile per sample from the paired tables, then check if a sample is far from other samples in the same condition.

# %%
# Convert paired contrast file into sample x program matrix.
control_long = paired[["control_sample", "program", "control_value"]].rename(
    columns={"control_sample": "sample", "control_value": "value"}
)
treated_long = paired[["treated_sample", "program", "treated_value"]].rename(
    columns={"treated_sample": "sample", "treated_value": "value"}
)
sample_long = pd.concat([control_long, treated_long], ignore_index=True).drop_duplicates(["sample", "program"])

sample_program = sample_long.pivot(index="sample", columns="program", values="value").reindex(columns=programs)
sample_meta = pd.DataFrame(index=sample_program.index)
sample_meta["condition"] = sample_meta.index.map(parse_condition)
sample_meta["replicate_de"] = sample_meta.index.map(parse_replicate)
sample_meta["treatment_group"] = np.where(sample_meta["condition"].str.contains("TDI"), "TDI", np.where(sample_meta["condition"].str.contains("IAG"), "IAG", "control"))

sample_program.to_csv(OUT_DIR / "sample_projected_program_profiles.tsv", sep="\t")
sample_meta.to_csv(OUT_DIR / "sample_projected_program_profile_metadata.tsv", sep="\t")

print("sample_program shape:", sample_program.shape)
display(sample_meta.head())

# %%
# PCA of sample program profiles.
coords, ev = simple_pca(sample_program, n_components=2, center=True, scale=True)
sample_plot_df = sample_meta.copy()
sample_plot_df["PC1"] = coords[:, 0]
sample_plot_df["PC2"] = coords[:, 1]

fig, ax = plt.subplots(figsize=(7, 5))
for group, sub in sample_plot_df.groupby("treatment_group"):
    ax.scatter(sub["PC1"], sub["PC2"], s=45, label=group)
    for sample, row in sub.iterrows():
        ax.text(row["PC1"], row["PC2"], sample, fontsize=5)
ax.axhline(0, color="black", linewidth=0.5)
ax.axvline(0, color="black", linewidth=0.5)
ax.set_xlabel(f"PC1 ({ev[0]*100:.1f}% var)")
ax.set_ylabel(f"PC2 ({ev[1]*100:.1f}% var)")
ax.set_title("PCA of sample-level projected cNMF profiles")
ax.legend(frameon=False)
fig.tight_layout()
for ext in ["png", "pdf", "svg"]:
    fig.savefig(OUT_DIR / f"pca_sample_projected_program_profiles.{ext}")
plt.show()

# %%
# Distance to condition centroid, leave-one-out within each condition.
centroid_records = []
for condition, samples in sample_meta.groupby("condition").groups.items():
    samples = list(samples)
    sub = sample_program.loc[samples]
    for sample in samples:
        others = [s for s in samples if s != sample]
        x = sub.loc[sample].to_numpy(dtype=float)
        if len(others) >= 1:
            centroid = sub.loc[others].mean(axis=0).to_numpy(dtype=float)
            corr = np.corrcoef(x, centroid)[0, 1] if np.std(x) > 0 and np.std(centroid) > 0 else np.nan
            dist = np.linalg.norm(x - centroid)
        else:
            corr = np.nan
            dist = np.nan
        centroid_records.append({
            "sample": sample,
            "condition": condition,
            "replicate_de": sample_meta.loc[sample, "replicate_de"],
            "corr_to_condition_leave_one_out_centroid": corr,
            "distance_to_condition_leave_one_out_centroid": dist,
        })

sample_centroid_qc = pd.DataFrame(centroid_records)
sample_centroid_qc["distance_rank"] = sample_centroid_qc["distance_to_condition_leave_one_out_centroid"].rank(ascending=False)
sample_centroid_qc["low_corr_flag"] = sample_centroid_qc["corr_to_condition_leave_one_out_centroid"] < 0.5
sample_centroid_qc.to_csv(OUT_DIR / "sample_condition_centroid_qc_projected_program_profiles.tsv", sep="\t", index=False)

display(sample_centroid_qc.sort_values("distance_to_condition_leave_one_out_centroid", ascending=False).head(30))

# %% [markdown]
# ## Optional: pseudobulk count-level replicate QC
#
# This section runs only if `pseudobulk_gene_counts.filtered.tsv.gz` and `pseudobulk_sample_metadata.filtered.tsv` are available.
# It checks whether a sample is unusual at the whole-transcriptome pseudobulk level, independent of cNMF projection.

# %%
if PSEUDOBULK_COUNTS_PATH.exists() and PSEUDOBULK_META_PATH.exists():
    pb_counts = pd.read_csv(PSEUDOBULK_COUNTS_PATH, sep="\t", index_col=0)
    pb_meta = pd.read_csv(PSEUDOBULK_META_PATH, sep="\t")
    print("pseudobulk counts shape before orientation check:", pb_counts.shape)
    print("metadata columns:", pb_meta.columns.tolist())

    # Determine sample id column.
    sample_col_candidates = ["sample", "sample_name", "pseudobulk_sample", "sample_id"]
    sample_col = next((c for c in sample_col_candidates if c in pb_meta.columns), None)
    if sample_col is None:
        raise ValueError("Could not find sample column in pseudobulk metadata.")
    pb_meta[sample_col] = pb_meta[sample_col].astype(str)
    pb_meta = pb_meta.set_index(sample_col, drop=False)

    # Counts are usually genes x samples from the exporter. Transpose if needed.
    if set(pb_meta.index).issubset(set(pb_counts.columns.astype(str))):
        pb_counts_sxg = pb_counts.T
    elif set(pb_meta.index).issubset(set(pb_counts.index.astype(str))):
        pb_counts_sxg = pb_counts
    else:
        raise ValueError("Could not align pseudobulk counts with metadata sample IDs.")

    pb_counts_sxg.index = pb_counts_sxg.index.astype(str)
    pb_counts_sxg = pb_counts_sxg.loc[pb_meta.index.intersection(pb_counts_sxg.index)]
    pb_meta = pb_meta.loc[pb_counts_sxg.index]

    lib_size = pb_counts_sxg.sum(axis=1)
    n_detected = (pb_counts_sxg > 0).sum(axis=1)
    pb_sample_qc = pb_meta.copy()
    pb_sample_qc["library_size"] = lib_size
    pb_sample_qc["n_detected_genes"] = n_detected
    pb_sample_qc.to_csv(OUT_DIR / "pseudobulk_sample_qc_library_size_detected_genes.tsv", sep="\t", index=False)

    display(pb_sample_qc[["library_size", "n_detected_genes"]].describe())
    display(pb_sample_qc.sort_values("library_size").head(20))

    # logCPM PCA.
    cpm = pb_counts_sxg.div(lib_size.replace(0, np.nan), axis=0) * 1e6
    logcpm = np.log1p(cpm.fillna(0.0))

    coords, ev = simple_pca(logcpm, n_components=2, center=True, scale=False)
    pb_plot_df = pb_meta.copy()
    pb_plot_df["PC1"] = coords[:, 0]
    pb_plot_df["PC2"] = coords[:, 1]
    if "condition" not in pb_plot_df.columns:
        pb_plot_df["condition"] = pb_plot_df.index.map(parse_condition)
    if "replicate_de" not in pb_plot_df.columns:
        pb_plot_df["replicate_de"] = pb_plot_df.index.map(parse_replicate)

    fig, ax = plt.subplots(figsize=(7, 5))
    for rep, sub in pb_plot_df.groupby("replicate_de"):
        ax.scatter(sub["PC1"], sub["PC2"], s=45, label=rep)
        for sample, row in sub.iterrows():
            ax.text(row["PC1"], row["PC2"], sample, fontsize=5)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel(f"PC1 ({ev[0]*100:.1f}% var)")
    ax.set_ylabel(f"PC2 ({ev[1]*100:.1f}% var)")
    ax.set_title("PCA of pseudobulk logCPM profiles")
    ax.legend(frameon=False)
    fig.tight_layout()
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT_DIR / f"pca_pseudobulk_logcpm_profiles.{ext}")
    plt.show()

    # Sample correlation heatmap.
    pb_corr = corrcoef_rows(logcpm)
    pb_corr.to_csv(OUT_DIR / "pseudobulk_logcpm_sample_correlation_matrix.tsv", sep="\t")
    plot_matrix(
        pb_corr,
        title="Pseudobulk logCPM sample correlation",
        out_prefix=OUT_DIR / "pseudobulk_logcpm_sample_correlation_heatmap",
        symmetric=True,
        figsize=(10, 9),
    )
else:
    print("Pseudobulk count/meta files not found. Skipping pseudobulk count-level QC.")
    print("Needed:", PSEUDOBULK_COUNTS_PATH)
    print("Needed:", PSEUDOBULK_META_PATH)

# %% [markdown]
# ## Sensitivity analysis: rerun program-score contrasts excluding suspect replicate/sample
#
# Fill in `SUSPECT_REPLICATES` or `SUSPECT_SAMPLES` after inspecting the QC tables above.
# This recomputes only the projected program-score paired tests. It does not rerun gene DE.

# %%
SUSPECT_REPLICATES = []  # example: ["rep2"]
SUSPECT_SAMPLES = []     # example: ["D2_LM+TDI_rep2"]

sensitivity = paired.copy()
if SUSPECT_REPLICATES:
    sensitivity = sensitivity[~sensitivity["replicate_de"].isin(SUSPECT_REPLICATES)].copy()
if SUSPECT_SAMPLES:
    sensitivity = sensitivity[
        ~sensitivity["control_sample"].isin(SUSPECT_SAMPLES)
        & ~sensitivity["treated_sample"].isin(SUSPECT_SAMPLES)
    ].copy()

print("Original rows:", paired.shape[0])
print("Sensitivity rows:", sensitivity.shape[0])
print("Remaining replicates:", sorted(sensitivity["replicate_de"].unique()))

sensitivity_records = []
for (contrast, program), sub in sensitivity.groupby(["contrast_name", "program"], sort=False):
    diffs = sub["paired_diff"].to_numpy(dtype=float)
    n = len(diffs)
    mean_diff = float(np.mean(diffs)) if n else np.nan
    sd = float(np.std(diffs, ddof=1)) if n > 1 else np.nan
    if stats is not None and n > 1:
        t_res = stats.ttest_1samp(diffs, popmean=0.0, nan_policy="omit")
        statistic = float(t_res.statistic)
        pvalue = float(t_res.pvalue)
    else:
        statistic = np.nan
        pvalue = np.nan
    sensitivity_records.append({
        "contrast_name": contrast,
        "program": program,
        "n_matched_pairs": n,
        "mean_paired_diff": mean_diff,
        "sd_paired_diff": sd,
        "statistic": statistic,
        "pvalue": pvalue,
    })

sensitivity_summary = pd.DataFrame(sensitivity_records)
sensitivity_summary["padj_bh_within_contrast"] = sensitivity_summary.groupby("contrast_name")["pvalue"].transform(bh_adjust)
sensitivity_summary.to_csv(OUT_DIR / "program_score_contrast_sensitivity_excluding_suspects.tsv", sep="\t", index=False)

key_sensitivity = sensitivity_summary[sensitivity_summary["program"].isin(KEY_PROGRAMS)].copy()
display(key_sensitivity.sort_values(["contrast_name", "program"]))

# %% [markdown]
# ## Suggested interpretation logic
#
# A replicate is suspicious if it is repeatedly flagged across unrelated contrasts, has low program-profile correlation to other replicates, and/or is an outlier in pseudobulk logCPM PCA.
#
# Do not exclude a replicate only because it disagrees for one key program. Exclude only if there is independent evidence of low quality, mislabeling, failed differentiation, or a global expression/profile outlier.

# %%
print("QC outputs saved in:", OUT_DIR)
for path in sorted(OUT_DIR.glob("*")):
    print(path.name)
