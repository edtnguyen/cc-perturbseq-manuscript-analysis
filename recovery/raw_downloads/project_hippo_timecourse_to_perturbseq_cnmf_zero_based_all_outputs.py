# %% [markdown]
# 
# # Project Hippo inhibitor time-course scRNA-seq onto perturb-seq cNMF programs
# 
# This notebook projects a new Hippo inhibitor time-course AnnData object onto the 50 cNMF gene programs discovered in the perturb-seq dataset. Program IDs are relabeled to zero-based `P0`–`P49` everywhere.
# 
# Reference matrix used:
# 
# ```text
# 030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50.gene_spectra_tpm.k_50.dt_0_2.txt
# ```
# 
# The notebook keeps the original `adata.X` unchanged. It computes a CPM/TPM-like working matrix from raw counts, projects cells onto fixed perturb-seq program spectra, then saves:
# 
# - projected `.h5ad`
# - projected program-score table
# - `.h5mu` with modalities `rna` and `perturbseq_cNMF`

# %%

from pathlib import Path
import re
import time
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.optimize import nnls

try:
    from IPython.display import display
except Exception:
    display = print

try:
    import muon as mu
except Exception as e:
    mu = None
    print("WARNING: muon could not be imported; .h5mu saving will be skipped unless muon is installed.")
    print(e)

warnings.filterwarnings("ignore")

# -------------------------
# Input/output paths
# -------------------------
ADATA_PATH = Path("/scratch/users/edtng/20260325_Hippo_timecourse/Hippo_timecourse_combined.h5ad")

SPECTRA_TPM_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50."
    "gene_spectra_tpm.k_50.dt_0_2.txt"
)

# Optional; used only to add labels if the file exists.
PROGRAM_META_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/data/program_activity_long_v4.csv"
)

OUT_DIR = ADATA_PATH.parent / "perturbseq_cNMF_projection"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_H5AD = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad"
OUT_H5MU = OUT_DIR / "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5mu"
OUT_SCORES_TSV = OUT_DIR / "Hippo_timecourse_combined.perturbseq_cNMF_50_dt_0_2.program_scores.zero_based.tsv.gz"

# -------------------------
# Projection settings
# -------------------------
TARGET_SUM = 1_000_000  # cNMF spectra_tpm is TPM/CPM-like, so use non-log normalized counts.

# Gene IDs in query AnnData.
# Use "var_names" if adata.var_names are gene symbols matching spectra columns.
# Otherwise set QUERY_GENE_ID_SOURCE = "column" and QUERY_GENE_COLUMN = "gene_name" or similar.
QUERY_GENE_ID_SOURCE = "column"  # "var_names" or "column"
QUERY_GENE_COLUMN = "gene"

# Counts source. If None, the notebook auto-detects a counts layer; if none exists, it uses adata.X.
COUNTS_LAYER = None
PREFERRED_COUNTS_LAYERS = ["counts", "raw_counts", "count", "raw"]

# Projection method:
# - "sklearn_fixed": fast fixed-H NMF-style projection; recommended default.
# - "fast_lstsq_clip": very fast approximate projection; clips negative least-squares coefficients to 0.
# - "nnls": exact scipy NNLS per cell; slow for large datasets.
PROJECTION_METHOD = "sklearn_fixed"
SKLEARN_MAX_ITER = 1000
SKLEARN_TOL = 1e-5

# Store obs columns for every program. This is convenient but adds 50 columns to adata.obs.
ADD_PROGRAM_COLUMNS_TO_OBS = True

# Keep raw usage and normalized usage. Normalized usage makes programs comparable within each cell.
NORMALIZE_USAGE_TO_SUM1 = True

print("ADATA_PATH:", ADATA_PATH)
print("SPECTRA_TPM_PATH:", SPECTRA_TPM_PATH)
print("OUT_DIR:", OUT_DIR)

# %%

# Load and inspect the new time-course AnnData
adata = sc.read_h5ad(ADATA_PATH)

print(adata)
print("\nobs columns:")
print(list(adata.obs.columns))
print("\nvar columns:")
print(list(adata.var.columns))
print("\nlayers:")
print(list(adata.layers.keys()))
print("\nobsm:")
print(list(adata.obsm.keys()))
print("\nuns keys:")
print(list(adata.uns.keys())[:50])

print("\nFirst obs rows:")
display(adata.obs.head())
print("\nFirst var rows:")
display(adata.var.head())

# %%

# Load reference perturb-seq cNMF spectra: programs x genes
spectra_df = pd.read_csv(SPECTRA_TPM_PATH, sep="\t", index_col=0)
spectra_df.index = spectra_df.index.astype(str)
spectra_df.columns = spectra_df.columns.astype(str)

# Make sure values are non-negative and numeric.
spectra_df = spectra_df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
spectra_df[spectra_df < 0] = 0.0

# cNMF output files label programs as 1..K, but the downstream perturb-seq
# program metadata/mdata columns use 0-based IDs. Convert spectra row labels:
# source row "1" -> program_id 0 -> display name "P0"
# source row "50" -> program_id 49 -> display name "P49"
spectra_source_programs = spectra_df.index.astype(str).tolist()
spectra_program_int = [int(x) for x in spectra_source_programs]

if spectra_program_int == list(range(1, len(spectra_program_int) + 1)):
    program_ids = [x - 1 for x in spectra_program_int]
else:
    # Fallback: if a future spectra file is already zero-based, keep it.
    program_ids = spectra_program_int

program_names = [f"P{i}" for i in program_ids]

# Keep the spectra rows in their original order, but relabel the index to P0..P49.
spectra_df.index = program_names
reference_genes = spectra_df.columns.astype(str).tolist()

program_name_to_id = dict(zip(program_names, program_ids))
program_name_to_source_row = dict(zip(program_names, spectra_source_programs))

expected_zero_based_names = [f"P{i}" for i in range(len(program_names))]
if program_names != expected_zero_based_names:
    raise ValueError(
        "Program relabeling did not produce contiguous zero-based names. "
        f"Got first names {program_names[:10]} and last names {program_names[-5:]}"
    )
if "P50" in program_names:
    raise ValueError("P50 found after relabeling; expected P0..P49 only.")

print("spectra_df shape, programs x genes:", spectra_df.shape)
print("source programs from spectra:", spectra_source_programs[:10])
print("zero-based program names:", program_names[:10])
print("first genes:", reference_genes[:10])
print("row sum range:", float(spectra_df.sum(axis=1).min()), float(spectra_df.sum(axis=1).max()))

# Optional program labels/colors from metadata, if available.
program_meta = None
program_label_map = {p: p for p in program_names}

if PROGRAM_META_PATH.exists():
    program_meta = pd.read_csv(PROGRAM_META_PATH).drop_duplicates("program_name")
    if "program_label" in program_meta.columns:
        # program_activity_long_v4.csv uses 0-based program IDs.
        for _, row in program_meta.iterrows():
            if pd.isna(row["program_name"]):
                continue
            pid = int(row["program_name"])
            pname = f"P{pid}"
            label = str(row["program_label"])
            program_label_map[pname] = label
    print("Loaded program metadata:", PROGRAM_META_PATH)
else:
    print("Program metadata not found; proceeding without labels:", PROGRAM_META_PATH)

# %%

def choose_counts_layer(adata, requested_layer=None, preferred_layers=None):
    """Return a count matrix and a human-readable source name."""
    preferred_layers = preferred_layers or []

    if requested_layer is not None:
        if requested_layer not in adata.layers:
            raise KeyError(f"Requested COUNTS_LAYER={requested_layer!r} not found. Available layers: {list(adata.layers.keys())}")
        return adata.layers[requested_layer], f"layers[{requested_layer!r}]"

    for layer in preferred_layers:
        if layer in adata.layers:
            return adata.layers[layer], f"layers[{layer!r}]"

    # Fall back to adata.X. The user said this matrix is currently unnormalized.
    return adata.X, "adata.X"


def get_query_gene_ids(adata, source="var_names", column=None):
    if source == "var_names":
        return np.asarray(adata.var_names.astype(str))
    if source == "column":
        if column is None:
            raise ValueError("QUERY_GENE_COLUMN must be set when QUERY_GENE_ID_SOURCE='column'.")
        if column not in adata.var.columns:
            raise KeyError(f"QUERY_GENE_COLUMN={column!r} not found in adata.var. Available columns: {list(adata.var.columns)}")
        return np.asarray(adata.var[column].astype(str))
    raise ValueError("QUERY_GENE_ID_SOURCE must be 'var_names' or 'column'.")


def as_csr_counts_matrix(X):
    """Convert an AnnData matrix/layer to CSR sparse format for efficient gene aggregation."""
    if sparse.issparse(X):
        return X.tocsr()
    # scRNA-seq counts should usually be sparse. Dense conversion can be memory-heavy.
    print("WARNING: count matrix is dense; converting to sparse CSR. This may use substantial memory.")
    return sparse.csr_matrix(X)


def build_gene_aggregation_matrix(query_gene_ids, common_genes):
    """Build a sparse var x common_genes matrix that sums duplicate query genes if needed."""
    common_to_col = {g: j for j, g in enumerate(common_genes)}
    rows = []
    cols = []
    data = []
    matched_vars = 0

    for i, gene in enumerate(query_gene_ids):
        j = common_to_col.get(gene)
        if j is not None:
            rows.append(i)
            cols.append(j)
            data.append(1.0)
            matched_vars += 1

    G = sparse.coo_matrix(
        (np.asarray(data, dtype=np.float32), (np.asarray(rows), np.asarray(cols))),
        shape=(len(query_gene_ids), len(common_genes)),
    ).tocsr()

    return G, matched_vars


def sparse_row_sums(X):
    return np.asarray(X.sum(axis=1)).ravel()

# %%

# Extract raw counts, compute full-library-size scaling, then subset/aggregate to reference genes.
counts_matrix_raw, counts_source = choose_counts_layer(
    adata,
    requested_layer=COUNTS_LAYER,
    preferred_layers=PREFERRED_COUNTS_LAYERS,
)
counts_csr = as_csr_counts_matrix(counts_matrix_raw)

query_gene_ids = get_query_gene_ids(
    adata,
    source=QUERY_GENE_ID_SOURCE,
    column=QUERY_GENE_COLUMN,
)
query_gene_ids = np.asarray([str(x) for x in query_gene_ids])

query_gene_set = set(query_gene_ids)
common_genes = [g for g in reference_genes if g in query_gene_set]

print("Counts source:", counts_source)
print("Query matrix shape:", counts_csr.shape)
print("Reference genes:", len(reference_genes))
print("Query genes:", len(query_gene_ids))
print("Shared genes:", len(common_genes))

if len(common_genes) < 500:
    raise ValueError(
        f"Only {len(common_genes)} shared genes. Check whether adata.var_names are gene symbols. "
        "If gene symbols are in adata.var, set QUERY_GENE_ID_SOURCE='column' and QUERY_GENE_COLUMN accordingly."
    )

# Full-library cell sums should be computed over all query genes, not only the shared genes.
cell_sums = sparse_row_sums(counts_csr).astype(np.float64)
zero_cells = int(np.sum(cell_sums <= 0))
print("Cells with zero total counts:", zero_cells)
if zero_cells > 0:
    cell_sums[cell_sums <= 0] = np.nan

scale = TARGET_SUM / cell_sums
scale = np.nan_to_num(scale, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)

# Aggregate duplicate gene symbols if present.
G, matched_vars = build_gene_aggregation_matrix(query_gene_ids, common_genes)
print("Matched query var rows, allowing duplicates:", matched_vars)

X_common_counts = counts_csr @ G
X_common_cpm = X_common_counts.multiply(scale[:, None]).tocsr().astype(np.float32)

print("X_common_cpm shape, cells x shared_genes:", X_common_cpm.shape)
print("X_common_cpm nnz:", X_common_cpm.nnz)
print("Mean nonzero value:", float(X_common_cpm.data.mean()) if X_common_cpm.nnz else 0.0)

# %%

# Align the reference spectra to the shared genes.
H = spectra_df.loc[program_names, common_genes].to_numpy(dtype=np.float32)  # programs x genes
H[H < 0] = 0

# Remove genes that have zero loading across all programs or zero expression across all cells.
h_nonzero = np.asarray(H.sum(axis=0)).ravel() > 0
x_nonzero = np.asarray(X_common_cpm.sum(axis=0)).ravel() > 0
keep_gene_mask = h_nonzero & x_nonzero

kept_genes = [g for g, keep in zip(common_genes, keep_gene_mask) if keep]
X_proj = X_common_cpm[:, keep_gene_mask].tocsr()
H_proj = H[:, keep_gene_mask]

print("Projection matrix X_proj, cells x genes:", X_proj.shape)
print("Reference matrix H_proj, programs x genes:", H_proj.shape)
print("Programs:", len(program_names))
print("Genes removed because zero in query or reference:", len(common_genes) - len(kept_genes))

if X_proj.shape[1] < 500:
    raise ValueError(f"Only {X_proj.shape[1]} usable genes after filtering. Check gene IDs/preprocessing.")

# %%

def initial_usage_lstsq_clip(X, H):
    """Fast nonnegative initialization by unconstrained least squares followed by clipping."""
    # X: cells x genes, H: programs x genes
    # Solve W H ~= X => W ~= X H.T (H H.T)^-1
    gram = H @ H.T
    gram_pinv = np.linalg.pinv(gram).astype(np.float32)
    XHt = X @ H.T
    W0 = np.asarray(XHt @ gram_pinv, dtype=np.float32)
    W0[W0 < 0] = 0
    return W0


def project_sklearn_fixed_H(X, H, max_iter=300, tol=1e-4):
    """Project cells onto fixed nonnegative spectra using sklearn's fixed-H NMF update."""
    try:
        from sklearn.decomposition import non_negative_factorization
    except Exception as e:
        raise ImportError("scikit-learn is not available; use PROJECTION_METHOD='fast_lstsq_clip' or 'nnls'.") from e

    n_cells = X.shape[0]
    n_programs = H.shape[0]
    print("Computing least-squares nonnegative initialization...")
    W0 = initial_usage_lstsq_clip(X, H)
    if not np.any(W0):
        W0 = np.full((n_cells, n_programs), 1e-6, dtype=np.float32)

    print("Running sklearn non_negative_factorization with fixed H...")
    W_fit, H_fit, n_iter = non_negative_factorization(
        X,
        W=W0,
        H=H.copy(),
        n_components=n_programs,
        init="custom",
        update_H=False,
        solver="cd",
        beta_loss="frobenius",
        max_iter=max_iter,
        tol=tol,
        random_state=0,
        verbose=1,
        shuffle=False,
    )
    print("sklearn iterations:", n_iter)
    W_fit = np.asarray(W_fit, dtype=np.float32)
    W_fit[W_fit < 0] = 0
    return W_fit


def project_nnls(X, H, report_every=1000):
    """Exact NNLS per cell. This can be slow for large datasets."""
    A = H.T.astype(np.float64)  # genes x programs
    n_cells = X.shape[0]
    n_programs = H.shape[0]
    W = np.zeros((n_cells, n_programs), dtype=np.float32)
    t0 = time.time()

    for i in range(n_cells):
        xi = X[i, :].toarray().ravel().astype(np.float64) if sparse.issparse(X) else np.asarray(X[i, :]).ravel().astype(np.float64)
        W[i, :], _ = nnls(A, xi)
        if (i + 1) % report_every == 0:
            elapsed = time.time() - t0
            print(f"NNLS projected {i + 1:,}/{n_cells:,} cells; elapsed {elapsed / 60:.1f} min")

    return W


t0 = time.time()

if PROJECTION_METHOD == "sklearn_fixed":
    usage_raw = project_sklearn_fixed_H(X_proj, H_proj, max_iter=SKLEARN_MAX_ITER, tol=SKLEARN_TOL)
elif PROJECTION_METHOD == "fast_lstsq_clip":
    usage_raw = initial_usage_lstsq_clip(X_proj, H_proj)
elif PROJECTION_METHOD == "nnls":
    usage_raw = project_nnls(X_proj, H_proj)
else:
    raise ValueError("PROJECTION_METHOD must be 'sklearn_fixed', 'fast_lstsq_clip', or 'nnls'.")

elapsed = time.time() - t0
print(f"Projection complete in {elapsed / 60:.2f} min")
print("usage_raw shape, cells x programs:", usage_raw.shape)
print("usage_raw min/max:", float(np.min(usage_raw)), float(np.max(usage_raw)))

# Normalize each cell's projected usage to sum to 1. This is usually easier to compare across cells.
usage_sum = usage_raw.sum(axis=1, keepdims=True).astype(np.float32)
usage_norm = np.divide(
    usage_raw,
    usage_sum,
    out=np.zeros_like(usage_raw, dtype=np.float32),
    where=usage_sum > 0,
)

print("usage_norm row sum range:", float(usage_norm.sum(axis=1).min()), float(usage_norm.sum(axis=1).max()))

# %%

# Create a zero-based program-score DataFrame for inspection and plotting.
program_score_df = pd.DataFrame(
    usage_norm if NORMALIZE_USAGE_TO_SUM1 else usage_raw,
    index=adata.obs_names.astype(str),
    columns=[str(p) for p in program_names],
)

print(program_score_df.shape)
display(program_score_df.iloc[:5, :10])

# Optional quick summary by time/day if such fields exist.
possible_time_cols = [c for c in ["day", "time", "timepoint", "treatment", "condition"] if c in adata.obs.columns]
print("Possible time/condition columns found:", possible_time_cols)

if "day" in adata.obs.columns:
    summary_by_day = program_score_df.join(adata.obs[["day"]]).groupby("day", observed=True).mean()
    display(summary_by_day.iloc[:, :10])

# %%

def sanitize_program_name(x):
    x = str(x)
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", x)
    if re.match(r"^[0-9]", x):
        x = "P" + x
    return x

# Attach scores to AnnData.
adata.obsm["X_perturbseq_cNMF_usage_raw"] = usage_raw.astype(np.float32)
adata.obsm["X_perturbseq_cNMF_usage"] = usage_norm.astype(np.float32)

if ADD_PROGRAM_COLUMNS_TO_OBS:
    for j, program in enumerate(program_names):
        col = f"perturbseq_cNMF_{sanitize_program_name(program)}_usage"
        adata.obs[col] = usage_norm[:, j] if NORMALIZE_USAGE_TO_SUM1 else usage_raw[:, j]

adata.uns["perturbseq_cNMF_projection"] = {
    "reference": "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50",
    "k": 50,
    "density_threshold": "0.2",
    "spectra_tpm_path": str(SPECTRA_TPM_PATH),
    "query_adata_path": str(ADATA_PATH),
    "counts_source": counts_source,
    "target_sum": TARGET_SUM,
    "projection_method": PROJECTION_METHOD,
    "program_names": [str(p) for p in program_names],
    "program_column_names": [f"perturbseq_cNMF_{sanitize_program_name(p)}_usage" for p in program_names],
    "program_ids_zero_based": [int(program_name_to_id[p]) for p in program_names],
    "spectra_source_programs_one_based": [str(program_name_to_source_row[p]) for p in program_names],
    "program_labels": [str(program_label_map.get(str(p), str(p))) for p in program_names],
    "n_shared_genes_initial": int(len(common_genes)),
    "n_genes_used": int(len(kept_genes)),
    "genes_used": [str(g) for g in kept_genes],
    "usage_obsm_raw": "X_perturbseq_cNMF_usage_raw",
    "usage_obsm_normalized": "X_perturbseq_cNMF_usage",
}

print("Added obsm keys:", [k for k in adata.obsm.keys() if "perturbseq" in k])
print("Added obs score columns:", sum(c.startswith("perturbseq_cNMF_") for c in adata.obs.columns))

# %%

# Save projected score table and projected AnnData.
program_score_df.to_csv(OUT_SCORES_TSV, sep="\t", compression="gzip")
print("Saved score table:", OUT_SCORES_TSV)

adata.write_h5ad(OUT_H5AD, compression="gzip")
print("Saved projected AnnData:", OUT_H5AD)

# %%

# Save a MuData file with two modalities:
#   - rna: original AnnData plus projection in obsm/obs
#   - perturbseq_cNMF: cells x projected perturb-seq programs
if mu is None:
    print("muon is not installed/importable; skipping .h5mu save.")
else:
    program_var = pd.DataFrame(index=pd.Index([str(p) for p in program_names], name="program"))
    program_var["program_name"] = [str(p) for p in program_names]
    program_var["program_id"] = [int(program_name_to_id[p]) for p in program_names]
    program_var["program_id_zero_based"] = [int(program_name_to_id[p]) for p in program_names]
    program_var["spectra_source_program_one_based"] = [str(program_name_to_source_row[p]) for p in program_names]
    program_var["program_label"] = [str(program_label_map.get(str(p), str(p))) for p in program_names]
    program_var["reference"] = "perturbseq_cNMF_50_dt_0_2"

    # Store minimal obs in the program modality; MuData will align cells by obs_names.
    program_obs = pd.DataFrame(index=adata.obs_names.copy())
    for col in ["day", "time", "timepoint", "treatment", "condition"]:
        if col in adata.obs.columns:
            program_obs[col] = adata.obs[col].values

    prog_adata = sc.AnnData(
        X=usage_norm.astype(np.float32),
        obs=program_obs,
        var=program_var,
    )
    prog_adata.layers["usage_raw"] = usage_raw.astype(np.float32)
    prog_adata.layers["usage_sum_to_1"] = usage_norm.astype(np.float32)
    prog_adata.uns["projection"] = adata.uns["perturbseq_cNMF_projection"]

    # Carry over UMAP coordinates if present.
    for key in ["X_umap", "X_pca"]:
        if key in adata.obsm:
            prog_adata.obsm[key] = adata.obsm[key].copy()

    mdata = mu.MuData({
        "rna": adata,
        "perturbseq_cNMF": prog_adata,
    })
    mdata.uns["description"] = (
        "Hippo inhibitor time-course scRNA-seq projected onto perturb-seq cNMF programs "
        "using gene_spectra_tpm.k_50.dt_0_2."
    )

    mdata.write_h5mu(OUT_H5MU)
    print("Saved projected MuData:", OUT_H5MU)

# %%

# Reload sanity checks.
adata_check = sc.read_h5ad(OUT_H5AD, backed="r")
print(adata_check)
print("Projection obsm keys:", [k for k in adata_check.obsm.keys() if "perturbseq" in k])
score_cols_check = [
    c for c in adata_check.obs.columns
    if c.startswith("perturbseq_cNMF_P") and c.endswith("_usage")
]
print("obs score columns:", len(score_cols_check))
print("first obs score columns:", score_cols_check[:5])
print("last obs score columns:", score_cols_check[-5:])
if any("P50" in c for c in score_cols_check):
    raise ValueError("Found P50 in obs score columns; expected zero-based P0..P49 only.")
adata_check.file.close()

if OUT_H5MU.exists() and mu is not None:
    mdata_check = mu.read_h5mu(OUT_H5MU, backed="r")
    print(mdata_check)
    print("modalities:", list(mdata_check.mod.keys()))
    print("program modality shape:", mdata_check.mod["perturbseq_cNMF"].shape)
    print("program modality first vars:", mdata_check.mod["perturbseq_cNMF"].var_names[:5].tolist())
    print("program modality last vars:", mdata_check.mod["perturbseq_cNMF"].var_names[-5:].tolist())
    if "P50" in set(mdata_check.mod["perturbseq_cNMF"].var_names.astype(str)):
        raise ValueError("Found P50 in MuData program modality; expected zero-based P0..P49 only.")
    mdata_check.file.close()

# %%

# Optional quick plot of a projected program over an existing UMAP.
# Set PROGRAM_TO_PLOT to a program name from program_names, e.g. program_names[0].
import matplotlib.pyplot as plt

PROGRAM_TO_PLOT = program_names[0]

if "X_umap" not in adata.obsm:
    print("No adata.obsm['X_umap'] found; skipping UMAP plot.")
else:
    j = program_names.index(PROGRAM_TO_PLOT)
    vals = usage_norm[:, j]
    coords = adata.obsm["X_umap"]

    fig, ax = plt.subplots(figsize=(5, 5), dpi=150)
    sca = ax.scatter(coords[:, 0], coords[:, 1], c=vals, s=1, linewidths=0, rasterized=True)
    ax.set_title(f"Perturb-seq cNMF program {PROGRAM_TO_PLOT}")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.colorbar(sca, ax=ax, label="Projected usage")
    plt.show()

