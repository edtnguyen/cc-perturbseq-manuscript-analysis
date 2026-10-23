# %% [markdown]
# # Audit Hippo time-course projection onto perturb-seq cNMF
#
# These cells check whether the saved projection was done in the intended space:
# raw/non-log counts -> CPM/TPM-like normalization -> fixed perturb-seq cNMF spectra -> zero-based P0..P49 usage scores.
#
# The last optional cell compares the fast sklearn fixed-H projection to exact NNLS on a small cell subset.

# %%
from pathlib import Path
import re
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.optimize import nnls
import matplotlib.pyplot as plt

ADATA_PATH = Path("/scratch/users/edtng/20260325_Hippo_timecourse/Hippo_timecourse_combined.h5ad")
PROJ_DIR = Path("/scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_projection")
H5AD_PROJ = PROJ_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad"
H5MU_PROJ = PROJ_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5mu"
SCORE_TSV = PROJ_DIR / "Hippo_timecourse_combined.perturbseq_cNMF_50_dt_0_2.program_scores.zero_based.tsv.gz"

SPECTRA_TPM_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50."
    "gene_spectra_tpm.k_50.dt_0_2.txt"
)

QUERY_GENE_COLUMN = "gene"
TARGET_SUM = 1_000_000

print("ADATA_PATH:", ADATA_PATH)
print("H5AD_PROJ:", H5AD_PROJ)
print("SCORE_TSV:", SCORE_TSV)
print("SPECTRA_TPM_PATH:", SPECTRA_TPM_PATH)

# %% [markdown]
# ## 1. Check the saved projection object

# %%
adata_proj = sc.read_h5ad(H5AD_PROJ, backed="r")
print(adata_proj)
print("obsm:", list(adata_proj.obsm.keys()))
print("uns:", list(adata_proj.uns.keys()))

usage_norm = np.asarray(adata_proj.obsm["X_perturbseq_cNMF_usage"])
usage_raw = np.asarray(adata_proj.obsm["X_perturbseq_cNMF_usage_raw"])

program_names = [str(x) for x in adata_proj.uns["perturbseq_cNMF_projection"]["program_names"]]
print("usage_norm shape:", usage_norm.shape)
print("usage_raw shape:", usage_raw.shape)
print("program_names first/last:", program_names[:5], program_names[-5:])

assert usage_norm.shape == (adata_proj.n_obs, 50)
assert usage_raw.shape == (adata_proj.n_obs, 50)
assert program_names == [f"P{i}" for i in range(50)]
assert "P50" not in program_names
print("PASS: projection has 212891 x 50 score matrices and zero-based P0..P49 labels.")

# %%
def score_qc(mat, name):
    print(f"\n{name}")
    print("shape:", mat.shape)
    print("dtype:", mat.dtype)
    print("NaN:", int(np.isnan(mat).sum()))
    print("Inf:", int(np.isinf(mat).sum()))
    print("min:", float(np.nanmin(mat)))
    print("max:", float(np.nanmax(mat)))
    print("negative count:", int((mat < -1e-8).sum()))
    print("zero fraction:", float((mat == 0).mean()))

score_qc(usage_raw, "usage_raw")
score_qc(usage_norm, "usage_norm")
row_sum_norm = usage_norm.sum(axis=1)
print("\nusage_norm row sum summary:")
print(pd.Series(row_sum_norm).describe())
print("max abs row-sum deviation from 1:", float(np.max(np.abs(row_sum_norm - 1))))

# %%
# Check obs columns and score TSV agreement.
score_cols = [f"perturbseq_cNMF_P{i}_usage" for i in range(50)]
missing = [c for c in score_cols if c not in adata_proj.obs.columns]
print("missing obs score columns:", missing)
assert not missing

obs_scores = adata_proj.obs[score_cols].to_numpy(dtype=np.float32)
print("max abs diff obs columns vs obsm usage_norm:", float(np.max(np.abs(obs_scores - usage_norm))))

if SCORE_TSV.exists():
    score_df = pd.read_csv(SCORE_TSV, sep="\t", index_col=0)
    print("score_df shape:", score_df.shape)
    print("score_df columns first/last:", score_df.columns[:5].tolist(), score_df.columns[-5:].tolist())
    print("score_df order matches adata:", np.array_equal(score_df.index.astype(str), adata_proj.obs_names.astype(str)))
    print("max abs diff TSV vs obsm:", float(np.max(np.abs(score_df.to_numpy(dtype=np.float32) - usage_norm))))

# %% [markdown]
# ## 2. Check original count matrix and gene matching
#
# Main questions:
# - Is `adata.X` raw/non-log counts?
# - Does `adata.var['gene']` match the spectra gene symbols?
# - Are duplicate gene symbols handled?

# %%
adata_q = sc.read_h5ad(ADATA_PATH, backed="r")
print(adata_q)
print("layers:", list(adata_q.layers.keys()))
print("var columns:", list(adata_q.var.columns))
print("obs columns:", list(adata_q.obs.columns))

spectra = pd.read_csv(SPECTRA_TPM_PATH, sep="\t", index_col=0)
spectra.index = spectra.index.astype(str)
spectra.columns = spectra.columns.astype(str)
spectra = spectra.apply(pd.to_numeric, errors="coerce").fillna(0)
spectra[spectra < 0] = 0

source_rows = spectra.index.tolist()
print("spectra shape:", spectra.shape)
print("spectra first source rows:", source_rows[:5])
assert [int(x) for x in source_rows] == list(range(1, 51))
print("PASS: spectra source rows are 1..50; relabel to P0..P49 is appropriate.")

query_genes = adata_q.var[QUERY_GENE_COLUMN].astype(str)
ref_genes = pd.Index(spectra.columns.astype(str))
query_gene_index = pd.Index(query_genes)
common = ref_genes.intersection(query_gene_index)
print("reference genes:", len(ref_genes))
print("query vars:", len(query_gene_index))
print("shared genes:", len(common))
print("fraction reference genes shared:", len(common) / len(ref_genes))
print("number of duplicated query gene names:", int(query_gene_index.duplicated().sum()))
print("number of unique duplicated gene symbols:", len(query_gene_index[query_gene_index.duplicated()].unique()))

# Top duplicated query genes among reference genes.
dup_counts = query_genes.value_counts()
dup_counts = dup_counts[(dup_counts > 1) & (dup_counts.index.isin(ref_genes))]
print("top duplicated query gene symbols that are in reference:")
display(dup_counts.head(20).to_frame("n_var_rows"))

# Program coverage after matching.
program_mass_coverage = spectra.loc[:, common].sum(axis=1) / spectra.sum(axis=1)
coverage_df = pd.DataFrame({
    "source_program_one_based": spectra.index.astype(str),
    "zero_based_program": [f"P{i}" for i in range(50)],
    "loading_mass_coverage": program_mass_coverage.to_numpy(),
}).sort_values("loading_mass_coverage")
display(coverage_df)
print(coverage_df["loading_mass_coverage"].describe())

# %%
# Check whether adata.X looks like raw counts on a random subset of cells.
rng = np.random.default_rng(0)
n_check = min(2000, adata_q.n_obs)
idx = np.sort(rng.choice(adata_q.n_obs, size=n_check, replace=False))

X_sub = adata_q[idx, :].X
if sparse.issparse(X_sub):
    X_sub = X_sub.tocsr()
    data = X_sub.data
    row_sums = np.asarray(X_sub.sum(axis=1)).ravel()
else:
    X_sub = np.asarray(X_sub)
    data = X_sub.ravel()
    data = data[data != 0]
    row_sums = X_sub.sum(axis=1)

# Sample nonzero values for integer-ish check.
if data.size > 100000:
    data_sample = rng.choice(data, size=100000, replace=False)
else:
    data_sample = data

frac_integerish = np.mean(np.abs(data_sample - np.round(data_sample)) < 1e-6) if data_sample.size else np.nan
print("adata.X subset row-sum summary:")
print(pd.Series(row_sums).describe())
print("nonzero value min/max:", float(data.min()), float(data.max()))
print("fraction sampled nonzero values integer-like:", float(frac_integerish))

if "total_counts" in adata_q.obs.columns:
    obs_total = pd.to_numeric(adata_q.obs.iloc[idx]["total_counts"], errors="coerce").to_numpy()
    print("correlation X row sums vs obs total_counts:", float(np.corrcoef(row_sums, obs_total)[0, 1]))
    print("median abs difference X row sums - obs total_counts:", float(np.median(np.abs(row_sums - obs_total))))

# %% [markdown]
# ## 3. Optional: validate sklearn fixed-H projection against exact NNLS on a subset
#
# This recomputes the projection for a small random subset of cells using both sklearn fixed-H and exact scipy NNLS. It is slower than the other checks. Use 200--1000 cells.

# %%
RUN_NNLS_SUBSET_CHECK = False
N_NNLS_CELLS = 500

if RUN_NNLS_SUBSET_CHECK:
    from sklearn.decomposition import non_negative_factorization
    from scipy.optimize import nnls

    # Reconstruct the same CPM matrix and H matrix for a subset only.
    rng = np.random.default_rng(1)
    cell_idx = np.sort(rng.choice(adata_q.n_obs, size=min(N_NNLS_CELLS, adata_q.n_obs), replace=False))

    X_counts = adata_q[cell_idx, :].X
    X_counts = X_counts.tocsr() if sparse.issparse(X_counts) else sparse.csr_matrix(X_counts)
    cell_sums = np.asarray(X_counts.sum(axis=1)).ravel().astype(float)
    scale = TARGET_SUM / cell_sums

    common_genes = [g for g in spectra.columns if g in set(query_gene_index)]
    common_to_col = {g: j for j, g in enumerate(common_genes)}
    rows, cols = [], []
    for i, g in enumerate(query_gene_index):
        j = common_to_col.get(g)
        if j is not None:
            rows.append(i)
            cols.append(j)
    G = sparse.coo_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(adata_q.n_vars, len(common_genes))).tocsr()
    X_common = (X_counts @ G).multiply(scale[:, None]).tocsr().astype(np.float32)

    spectra_z = spectra.copy()
    spectra_z.index = [f"P{i}" for i in range(50)]
    H = spectra_z.loc[[f"P{i}" for i in range(50)], common_genes].to_numpy(dtype=np.float32)
    keep = (np.asarray(H.sum(axis=0)).ravel() > 0) & (np.asarray(X_common.sum(axis=0)).ravel() > 0)
    Xp = X_common[:, keep].tocsr()
    Hp = H[:, keep]

    # sklearn fixed-H projection.
    gram = Hp @ Hp.T
    W0 = np.asarray((Xp @ Hp.T) @ np.linalg.pinv(gram).astype(np.float32), dtype=np.float32)
    W0[W0 < 0] = 0
    W_skl, _, n_iter = non_negative_factorization(
        Xp, W=W0, H=Hp.copy(), n_components=50, init="custom", update_H=False,
        solver="cd", beta_loss="frobenius", max_iter=1000, tol=1e-5, random_state=0,
        verbose=0, shuffle=False,
    )
    W_skl = np.asarray(W_skl, dtype=np.float32)
    W_skl[W_skl < 0] = 0
    W_skl_norm = W_skl / W_skl.sum(axis=1, keepdims=True)

    # exact NNLS.
    A = Hp.T.astype(float)
    W_nnls = np.zeros_like(W_skl, dtype=np.float32)
    for i in range(Xp.shape[0]):
        xi = Xp[i, :].toarray().ravel().astype(float)
        W_nnls[i, :], _ = nnls(A, xi)
    W_nnls_norm = W_nnls / W_nnls.sum(axis=1, keepdims=True)

    rows = []
    for j in range(50):
        rows.append({
            "program": f"P{j}",
            "corr_sklearn_vs_nnls": np.corrcoef(W_skl_norm[:, j], W_nnls_norm[:, j])[0, 1],
            "mean_abs_diff": float(np.mean(np.abs(W_skl_norm[:, j] - W_nnls_norm[:, j]))),
        })
    nnls_check = pd.DataFrame(rows).sort_values("corr_sklearn_vs_nnls")
    print("sklearn subset n_iter:", n_iter)
    display(nnls_check)
else:
    print("Skipping NNLS subset check. Set RUN_NNLS_SUBSET_CHECK=True to run it.")

# %%
adata_proj.file.close()
adata_q.file.close()
