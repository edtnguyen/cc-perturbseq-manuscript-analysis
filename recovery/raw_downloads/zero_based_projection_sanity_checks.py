# %% [markdown]
# # Zero-based sanity checks for Hippo time-course projection
#
# This notebook checks the saved projection after relabeling perturb-seq cNMF
# programs to zero-based names: P0, P1, ..., P49.
#
# It does not recompute the projection.

# %%
from pathlib import Path
import re

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

try:
    import muon as mu
except Exception:
    mu = None

# -------------------------
# Paths
# -------------------------
OUT_DIR = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_projection")

H5AD_PATH = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad"
H5MU_PATH = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5mu"
SCORE_TSV = OUT_DIR / "Hippo_timecourse_combined.perturbseq_cNMF_50_dt_0_2.program_scores.zero_based.tsv.gz"

OLD_H5AD_PATH = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.h5ad"

SPECTRA_TPM_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50."
    "gene_spectra_tpm.k_50.dt_0_2.txt"
)

print("H5AD_PATH:", H5AD_PATH)
print("H5MU_PATH:", H5MU_PATH)
print("SCORE_TSV:", SCORE_TSV)
print("SPECTRA_TPM_PATH:", SPECTRA_TPM_PATH)

# %%
# Helper functions

def p_to_int(program):
    """Convert 'P38', 38, '38' to integer 38."""
    s = str(program)
    s = re.sub(r"^perturbseq_cNMF_", "", s)
    s = re.sub(r"_usage$", "", s)
    s = s.replace("P", "")
    return int(s)


def expected_zero_based_program_names(n_programs):
    return [f"P{i}" for i in range(n_programs)]


def obs_score_cols_zero_based(adata):
    cols = [
        c for c in adata.obs.columns
        if c.startswith("perturbseq_cNMF_P") and c.endswith("_usage")
    ]
    return sorted(cols, key=p_to_int)


def get_program_names_from_saved_adata(adata, n_programs):
    """Prefer saved uns program names. Fallback to obs score columns. Final fallback P0..Pn-1."""
    uns = adata.uns.get("perturbseq_cNMF_projection", {})

    if "program_names" in uns:
        names = [str(x) for x in uns["program_names"]]
        if len(names) == n_programs:
            return names

    cols = obs_score_cols_zero_based(adata)
    if len(cols) == n_programs:
        return [re.sub(r"^perturbseq_cNMF_", "", c).replace("_usage", "") for c in cols]

    return expected_zero_based_program_names(n_programs)


def load_zero_based_spectra(path):
    """Load spectra and relabel rows to P0..P49 if the source rows are 1..50."""
    spectra = pd.read_csv(path, sep="\t", index_col=0)
    spectra.index = spectra.index.astype(str)
    spectra.columns = spectra.columns.astype(str)
    spectra = spectra.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    spectra[spectra < 0] = 0.0

    source = spectra.index.astype(str).tolist()
    try:
        source_int = [int(x) for x in source]
    except Exception:
        source_int = None

    if source_int == list(range(1, len(source_int) + 1)):
        zero_names = [f"P{i - 1}" for i in source_int]
        mapping = pd.DataFrame({
            "spectra_source_program": source,
            "zero_based_program": zero_names,
            "old_one_based_label": [f"P{i}" for i in source_int],
        })
        spectra.index = zero_names
    elif source_int == list(range(0, len(source_int))):
        zero_names = [f"P{i}" for i in source_int]
        mapping = pd.DataFrame({
            "spectra_source_program": source,
            "zero_based_program": zero_names,
            "old_one_based_label": [f"P{i + 1}" for i in source_int],
        })
        spectra.index = zero_names
    else:
        # Fallback: preserve source labels but force P prefix for numeric names.
        zero_names = [x if x.startswith("P") else f"P{x}" for x in source]
        mapping = pd.DataFrame({
            "spectra_source_program": source,
            "zero_based_program": zero_names,
            "old_one_based_label": pd.NA,
        })
        spectra.index = zero_names

    return spectra, mapping


def matrix_qc(mat, name):
    print(f"\n{name}")
    print("shape:", mat.shape)
    print("dtype:", mat.dtype)
    print("NaN count:", int(np.isnan(mat).sum()))
    print("Inf count:", int(np.isinf(mat).sum()))
    print("min:", float(np.nanmin(mat)))
    print("max:", float(np.nanmax(mat)))
    print("negative count:", int((mat < -1e-8).sum()))
    print("exact zero fraction:", float(np.mean(mat == 0)))

# %%
# Load saved zero-based projection.
# Use backed='r' so adata.X is not loaded into memory.

if not H5AD_PATH.exists():
    raise FileNotFoundError(
        f"Missing zero-based h5ad: {H5AD_PATH}\n"
        "Use the relabel-only notebook/script first, or change H5AD_PATH."
    )

adata = sc.read_h5ad(H5AD_PATH, backed="r")

usage_norm = np.asarray(adata.obsm["X_perturbseq_cNMF_usage"])
usage_raw = np.asarray(adata.obsm["X_perturbseq_cNMF_usage_raw"])

n_cells, n_programs = usage_norm.shape
program_names = get_program_names_from_saved_adata(adata, n_programs)
expected_programs = expected_zero_based_program_names(n_programs)

print(adata)
print("usage_norm shape:", usage_norm.shape)
print("usage_raw shape:", usage_raw.shape)
print("program_names first/last:", program_names[:5], program_names[-5:])
print("expected first/last:", expected_programs[:5], expected_programs[-5:])

assert usage_norm.shape == usage_raw.shape
assert n_programs == 50

if program_names != expected_programs:
    raise ValueError(
        "Program names are not zero-based P0..P49.\n"
        f"Observed first/last: {program_names[:5]} ... {program_names[-5:]}"
    )

print("Zero-based program names OK.")

# %%
# Check obs score columns match zero-based names and match obsm values.

score_cols = obs_score_cols_zero_based(adata)
print("Number of obs score columns:", len(score_cols))
print("first obs score columns:", score_cols[:5])
print("last obs score columns:", score_cols[-5:])

expected_cols = [f"perturbseq_cNMF_{p}_usage" for p in expected_programs]
missing = sorted(set(expected_cols) - set(score_cols), key=lambda x: p_to_int(x))
extra = sorted(set(score_cols) - set(expected_cols), key=lambda x: p_to_int(x))

print("missing expected cols:", missing[:10], "n=", len(missing))
print("extra score cols:", extra[:10], "n=", len(extra))

assert not missing
assert not extra

obs_scores = adata.obs[expected_cols].to_numpy(dtype=np.float32)
print("Max abs difference obs score columns vs obsm usage_norm:", float(np.max(np.abs(obs_scores - usage_norm))))

# %%
# Basic numerical QC.

matrix_qc(usage_raw, "usage_raw")
matrix_qc(usage_norm, "usage_norm")

row_sums_raw = usage_raw.sum(axis=1)
row_sums_norm = usage_norm.sum(axis=1)

print("\nRaw row-sum summary:")
display(pd.Series(row_sums_raw).describe())

print("\nNormalized row-sum summary:")
display(pd.Series(row_sums_norm).describe())

print("Max abs deviation from normalized row sum = 1:", float(np.max(np.abs(row_sums_norm - 1))))

# %%
# Check score TSV, if available.

if SCORE_TSV.exists():
    score_df = pd.read_csv(SCORE_TSV, sep="\t", index_col=0)
    print("score_df shape:", score_df.shape)
    print("score_df first/last columns:", score_df.columns[:5].tolist(), score_df.columns[-5:].tolist())
    print("cell order matches AnnData:", np.array_equal(score_df.index.astype(str), adata.obs_names.astype(str)))
    assert score_df.columns.tolist() == expected_programs
    tsv_values = score_df.to_numpy(dtype=np.float32)
    print("Max abs difference TSV vs obsm usage_norm:", float(np.max(np.abs(tsv_values - usage_norm))))
else:
    print("No score TSV found:", SCORE_TSV)

# %%
# Check h5mu program modality, if available.

if H5MU_PATH.exists() and mu is not None:
    mdata = mu.read_h5mu(H5MU_PATH, backed="r")
    print(mdata)
    print("modalities:", list(mdata.mod.keys()))

    prog = mdata.mod["perturbseq_cNMF"]
    print("program modality shape:", prog.shape)
    print("program var_names first/last:", prog.var_names[:5].tolist(), prog.var_names[-5:].tolist())

    assert prog.n_obs == n_cells
    assert prog.n_vars == n_programs
    assert prog.var_names.tolist() == expected_programs

    prog_X = prog.X
    if sparse.issparse(prog_X):
        prog_X = prog_X.toarray()
    else:
        prog_X = np.asarray(prog_X)

    print("Max abs difference h5mu program X vs h5ad usage_norm:", float(np.max(np.abs(prog_X - usage_norm))))
    mdata.file.close()
else:
    print("Skipping h5mu check. Missing file or muon unavailable.")

# %%
# If old one-based file exists, confirm relabeling was only a label shift.
# Expected: old matrix column j equals new matrix column j, but labels changed
# old P1 -> new P0, old P2 -> new P1, ..., old P50 -> new P49.

if OLD_H5AD_PATH.exists():
    old = sc.read_h5ad(OLD_H5AD_PATH, backed="r")
    old_usage = np.asarray(old.obsm["X_perturbseq_cNMF_usage"])

    print("old usage shape:", old_usage.shape)
    print("new usage shape:", usage_norm.shape)
    print("max abs diff old usage vs new usage:", float(np.max(np.abs(old_usage - usage_norm))))

    mapping = pd.DataFrame({
        "matrix_column": np.arange(n_programs),
        "old_one_based_label": [f"P{i}" for i in range(1, n_programs + 1)],
        "new_zero_based_label": [f"P{i}" for i in range(n_programs)],
    })
    display(mapping.head(10))
    display(mapping.tail(10))
    old.file.close()
else:
    print("Old one-based H5AD not found; skipping old/new comparison:", OLD_H5AD_PATH)

# %%
# Program-level distribution summary with zero-based names.

program_summary = pd.DataFrame({
    "program": program_names,
    "mean": usage_norm.mean(axis=0),
    "median": np.median(usage_norm, axis=0),
    "p90": np.percentile(usage_norm, 90, axis=0),
    "p95": np.percentile(usage_norm, 95, axis=0),
    "p99": np.percentile(usage_norm, 99, axis=0),
    "max": usage_norm.max(axis=0),
    "frac_gt_0.01": (usage_norm > 0.01).mean(axis=0),
    "frac_gt_0.05": (usage_norm > 0.05).mean(axis=0),
    "frac_gt_0.10": (usage_norm > 0.10).mean(axis=0),
})
program_summary["old_one_based_label"] = program_summary["program"].map(lambda p: f"P{p_to_int(p) + 1}")
program_summary = program_summary.sort_values("mean", ascending=False)

display(program_summary)

# %%
# Dominant-program distribution, zero-based.

dominant_idx = np.argmax(usage_norm, axis=1)
dominant_program = np.array(program_names, dtype=object)[dominant_idx]
dominant_score = usage_norm.max(axis=1)

entropy = -np.sum(
    np.where(usage_norm > 0, usage_norm * np.log2(usage_norm), 0),
    axis=1,
)

dominant_df = pd.DataFrame({
    "dominant_program": dominant_program,
    "dominant_score": dominant_score,
    "usage_entropy": entropy,
}, index=adata.obs_names.astype(str))

print("Dominant score summary:")
display(dominant_df["dominant_score"].describe())

print("Usage entropy summary:")
display(dominant_df["usage_entropy"].describe())

program_counts = (
    dominant_df["dominant_program"]
    .value_counts()
    .rename_axis("program")
    .reset_index(name="n_cells")
)
program_counts["frac_cells"] = program_counts["n_cells"] / n_cells
program_counts["old_one_based_label"] = program_counts["program"].map(lambda p: f"P{p_to_int(p) + 1}")
display(program_counts)

plt.figure(figsize=(5, 3))
plt.hist(dominant_score, bins=100)
plt.xlabel("Dominant normalized program score")
plt.ylabel("Number of cells")
plt.title("Dominant program score per cell")
plt.tight_layout()
plt.show()

# %%
# QC correlations using BOTH raw and normalized program scores.
# This matters because usage_norm is compositional: every cell sums to 1.

candidate_qc_cols = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "guide_umi_counts",
    "guides_per_cell",
]
available_qc_cols = [c for c in candidate_qc_cols if c in adata.obs.columns]
print("Available QC columns:", available_qc_cols)


def qc_corr_table(score_mat, score_type, use_spearman=False):
    records = []
    for qc_col in available_qc_cols:
        qc = pd.to_numeric(adata.obs[qc_col], errors="coerce").to_numpy()
        valid_qc = np.isfinite(qc)

        for j, program in enumerate(program_names):
            y = score_mat[:, j]
            valid = valid_qc & np.isfinite(y)
            if valid.sum() <= 10:
                pearson = np.nan
                spearman = np.nan
            else:
                pearson = np.corrcoef(qc[valid], y[valid])[0, 1]
                if use_spearman:
                    spearman = spearmanr(qc[valid], y[valid], nan_policy="omit").correlation
                else:
                    spearman = np.nan

            records.append({
                "score_type": score_type,
                "qc_col": qc_col,
                "program": program,
                "old_one_based_label": f"P{p_to_int(program) + 1}",
                "pearson_corr": pearson,
                "spearman_corr": spearman,
            })

    out = pd.DataFrame(records)
    out["abs_pearson"] = out["pearson_corr"].abs()
    out["abs_spearman"] = out["spearman_corr"].abs()
    return out

qc_raw = qc_corr_table(usage_raw, "raw", use_spearman=False)
qc_norm = qc_corr_table(usage_norm, "normalized", use_spearman=False)
qc_both = pd.concat([qc_raw, qc_norm], ignore_index=True)

# Put normalized first because those are likely what you will use downstream.
display(qc_both.sort_values(["abs_pearson", "score_type"], ascending=[False, True]).head(50))

# %%
# Mito/ribo/rRNA loading fractions in the reference spectra, using zero-based program names.
# This helps distinguish true mitochondrial/ribosomal programs from compositional artifacts.

spectra_zb, spectra_mapping = load_zero_based_spectra(SPECTRA_TPM_PATH)
print("spectra_zb shape:", spectra_zb.shape)
display(spectra_mapping.head(10))

assert spectra_zb.index.tolist() == expected_programs

genes = pd.Index(spectra_zb.columns.astype(str))
mt_mask = genes.str.startswith("MT-")
ribo_mask = genes.str.match(r"^(RPL|RPS|MRPL|MRPS)", na=False)
rRNA_mask = genes.str.contains("rRNA", case=False, regex=False)

spectra_qc = []
for program in expected_programs:
    row = spectra_zb.loc[program]
    denom = float(row.sum())
    if denom <= 0:
        denom = np.nan
    top_genes = row.sort_values(ascending=False).head(25).index.tolist()
    spectra_qc.append({
        "program": program,
        "old_one_based_label": f"P{p_to_int(program) + 1}",
        "mt_loading_fraction": float(row.loc[mt_mask].sum() / denom),
        "ribo_loading_fraction": float(row.loc[ribo_mask].sum() / denom),
        "rRNA_loading_fraction": float(row.loc[rRNA_mask].sum() / denom),
        "top_25_genes": ", ".join(top_genes),
    })

spectra_qc = pd.DataFrame(spectra_qc)

display(
    spectra_qc.sort_values(
        ["mt_loading_fraction", "ribo_loading_fraction", "rRNA_loading_fraction"],
        ascending=False,
    ).head(20)
)

# %%
# Merge normalized QC correlations with spectra QC annotations.
# Interpretation:
# - high pct_counts_mt correlation + high mt_loading_fraction => likely mitochondrial/QC-like program.
# - high pct_counts_ribo correlation + high ribo/rRNA loading => likely ribosomal/rRNA-like program.
# - high correlation but low loading => could be compositional, cell-state, or sample confounding.

norm_qc = qc_norm.copy()
merged_qc = norm_qc.merge(spectra_qc, on=["program", "old_one_based_label"], how="left")
merged_qc["abs_pearson"] = merged_qc["pearson_corr"].abs()

interesting_qc = merged_qc[
    merged_qc["qc_col"].isin(["pct_counts_mt", "pct_counts_ribo", "total_counts", "n_genes_by_counts"])
].sort_values("abs_pearson", ascending=False)

display(interesting_qc.head(50))

# %%
# Plot strongest mito/ribo-associated programs using zero-based labels.

plot_rows = interesting_qc[
    interesting_qc["qc_col"].isin(["pct_counts_mt", "pct_counts_ribo"])
].head(8)

rng = np.random.default_rng(0)
sample_idx = rng.choice(n_cells, size=min(30000, n_cells), replace=False)

for _, row in plot_rows.iterrows():
    qc_col = row["qc_col"]
    program = row["program"]
    j = p_to_int(program)

    x = pd.to_numeric(adata.obs[qc_col], errors="coerce").to_numpy()[sample_idx]
    y = usage_norm[sample_idx, j]

    plt.figure(figsize=(4, 3))
    plt.scatter(x, y, s=1, alpha=0.25)
    plt.xlabel(qc_col)
    plt.ylabel(f"{program} normalized usage")
    plt.title(f"{program} vs {qc_col} | old label {row['old_one_based_label']}")
    plt.tight_layout()
    plt.show()

# %%
# Sample/timepoint/group-level summary.
# The projection notebook did not auto-detect a day/time column, so inspect sample fields manually.

candidate_group_cols = [
    "sample_name",
    "biological_sample",
    "sample_id",
    "pool",
    "round1",
    "plate",
    "sublibID",
    "well",
]

for col in candidate_group_cols:
    if col in adata.obs.columns:
        print("\n" + col)
        print(adata.obs[col].value_counts(dropna=False).head(50))

# %%
# Choose the column encoding time/treatment and summarize program usage by group.
# Edit GROUP_COL after inspecting the previous output.

GROUP_COL = "sample_name"  # change if biological_sample or sample_id is better

if GROUP_COL in adata.obs.columns:
    score_df = pd.DataFrame(usage_norm, index=adata.obs_names.astype(str), columns=program_names)
    score_df[GROUP_COL] = adata.obs[GROUP_COL].astype(str).to_numpy()

    group_mean = score_df.groupby(GROUP_COL, observed=True)[program_names].mean()
    group_n = score_df.groupby(GROUP_COL, observed=True).size().rename("n_cells")

    display(group_n.to_frame())
    display(group_mean)

    # Z-score each program across groups for display only.
    mat_z = (group_mean - group_mean.mean(axis=0)) / group_mean.std(axis=0).replace(0, np.nan)

    plt.figure(figsize=(12, max(3, 0.25 * mat_z.shape[0])))
    plt.imshow(mat_z.to_numpy(), aspect="auto")
    plt.colorbar(label="Group mean usage, z-scored per program")
    plt.yticks(range(mat_z.shape[0]), mat_z.index)
    plt.xticks(range(mat_z.shape[1]), mat_z.columns, rotation=90)
    plt.title(f"Projected perturb-seq cNMF usage by {GROUP_COL}")
    plt.tight_layout()
    plt.show()

    dynamic = pd.DataFrame({
        "program": program_names,
        "old_one_based_label": [f"P{p_to_int(p) + 1}" for p in program_names],
        "min_group_mean": group_mean.min(axis=0).values,
        "max_group_mean": group_mean.max(axis=0).values,
    })
    dynamic["range"] = dynamic["max_group_mean"] - dynamic["min_group_mean"]
    dynamic = dynamic.sort_values("range", ascending=False)
    display(dynamic.head(30))
else:
    print(f"GROUP_COL={GROUP_COL!r} not found in adata.obs")

# %%
# Print top genes for QC-associated or dynamic programs.

PROGRAMS_TO_INSPECT = sorted(set(
    interesting_qc.head(12)["program"].tolist()
    + (dynamic.head(12)["program"].tolist() if "dynamic" in globals() else [])
), key=p_to_int)

for program in PROGRAMS_TO_INSPECT:
    row = spectra_zb.loc[program]
    top = row.sort_values(ascending=False).head(30).index.tolist()
    old_label = f"P{p_to_int(program) + 1}"
    print(f"\n{program}  old_one_based_label={old_label}")
    print(", ".join(top))

# %%
# Close backed AnnData handle.
adata.file.close()
print("Done.")
