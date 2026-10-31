"""Rebuild the manuscript's selected Hippo fixed-reference projection.

Scientific source: ``recovery/raw_downloads/reproject_cnmf_compare_gene_sets.ipynb``,
zero-based cells 3, 5, 7, 9, 11, and 13. These cells match the archived Sherlock
text mirror (source notebook SHA-256
``7a57ce043047b360c5a0384b13d3cbcd4113674f89690bfbb4621720960ce694``).

Only ``top_specific_200`` is rebuilt. The original comparison's other modes and
the later Hunter delta-to-control export are not needed by Figures 3f/4c. Here X
is normalized top-specific usage; the explicitly named obsm matrices retain the
original selected-mode definitions. No cNMF basis is refitted.
"""
from __future__ import annotations

import gc
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import re
import tempfile
import warnings

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import non_negative_factorization
from sklearn.exceptions import ConvergenceWarning


MODE = "top_specific_200"
TOP_SPECIFIC_N = 200
TARGET_SUM = 1e6
SKLEARN_MAX_ITER = 1000
SKLEARN_TOL = 1e-5
RANDOM_STATE = 0
CLR_EPS = 1e-6
LOG_BASE = 2
USAGE_KEY = f"X_cnmf_usage_{MODE}"
RAW_KEY = f"X_cnmf_usage_raw_{MODE}"
CLR_KEY = f"X_cnmf_clr_{MODE}"
SOURCE_NOTEBOOK_SHA256 = "7a57ce043047b360c5a0384b13d3cbcd4113674f89690bfbb4621720960ce694"
EXPECTED_QUERY_SHAPE = (212891, 62757)
EXPECTED_REFERENCE_SHAPE = (50, 17538)
EXPECTED_SELECTED_GENES = 7268
EXPECTED_QC_PASS = 99416

COND_PCT_RIBO_CUTOFFS = {
    "D0_hPSC": 2.4, "D1_MPS": 2.3, "D1_MPS+TDI": 2.5,
    "D1_MPS+IAG": 2.2, "D1.33_LM": 1.9, "D1.33_LM+TDI": 2.1,
    "D1.33_LM+IAG": 2.2, "D1.67_LM": 2.0, "D1.67_LM+TDI": 2.2,
    "D1.67_LM+IAG": 1.9, "D2_LM": 2.0, "D2_LM+TDI": 2.4,
    "D2_LM+IAG": 1.9, "D3_AE": 1.0, "D3_AE+TDI": 2.0,
    "D3_AE+IAG": 1.0, "D2_LM+TDI->D3_AE+TDI": 2.5,
    "D2_LM+IAG->D3_AE+IAG": 1.0,
}
COND_PCT_MT_CUTOFFS = {
    "D0_hPSC": 5.5, "D1_MPS": 3.0, "D1_MPS+TDI": 3.0,
    "D1_MPS+IAG": 4.0, "D1.33_LM": 3.5, "D1.33_LM+TDI": 3.75,
    "D1.33_LM+IAG": 5.0, "D1.67_LM": 4.0, "D1.67_LM+TDI": 5.5,
    "D1.67_LM+IAG": 3.5, "D2_LM": 6.5, "D2_LM+TDI": 8.0,
    "D2_LM+IAG": 3.5, "D3_AE": 3.25, "D3_AE+TDI": 7.0,
    "D3_AE+IAG": 3.75, "D2_LM+TDI->D3_AE+TDI": 7.0,
    "D2_LM+IAG->D3_AE+IAG": 3.75,
}


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _signature(path, *, hash_content=False):
    path = Path(path).resolve()
    stat = path.stat()
    if stat.st_size <= 0:
        raise ValueError(f"Input is empty: {path}")
    result = {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if hash_content:
        result["sha256"] = _sha256(path)
    return result


def _versions():
    result = {}
    for name in ("anndata", "h5py", "numpy", "pandas", "scipy", "scikit-learn"):
        try:
            result[name] = version(name)
        except PackageNotFoundError:
            result[name] = "not installed"
    return result


def _finite_nonnegative(values, label):
    values = values.data if sparse.issparse(values) else np.asarray(values).ravel()
    # Do not create a boolean array the size of the full 604-million-entry query.
    for start in range(0, len(values), 1_000_000):
        block = values[start:start + 1_000_000]
        if not np.isfinite(block).all() or (block < 0).any():
            raise ValueError(f"{label} contains nonfinite or negative values")


def read_spectra(path, *, strict_historical=True):
    """Read the original program-by-gene TPM reference with explicit P0 mapping."""
    ref = pd.read_csv(path, sep="\t", index_col=0)
    ref.index = ref.index.astype(str)
    ref.columns = ref.columns.astype(str)
    if not ref.index.is_unique or not ref.columns.is_unique:
        raise ValueError("Reference program IDs and gene symbols must be unique")
    try:
        numeric = [float(value) for value in ref.index]
        ids = [int(value) for value in numeric]
    except ValueError:
        ids = None
    else:
        if any(float(i) != v for i, v in zip(ids, numeric)):
            raise ValueError("Reference program labels must be integer IDs")
        if ids == list(range(1, len(ids) + 1)):
            ref.index = [f"P{i}" for i in range(len(ids))]
        elif ids == list(range(len(ids))):
            ref.index = [f"P{i}" for i in ids]
        else:
            raise ValueError("Numeric reference rows must be ordered 1..K or 0..K-1")
    expected = [f"P{i}" for i in range(ref.shape[0])]
    if set(ref.index) != set(expected):
        raise ValueError(f"Reference must contain exactly {expected}")
    ref = ref.loc[expected].apply(pd.to_numeric, errors="raise")
    _finite_nonnegative(ref.to_numpy(), "Reference spectra")
    if strict_historical and ref.shape != EXPECTED_REFERENCE_SHAPE:
        raise ValueError(f"Reference shape {ref.shape}; expected {EXPECTED_REFERENCE_SHAPE}")
    return ref


def top_specific_genes(ref, query_gene_ids, *, n_top=TOP_SPECIFIC_N):
    """Original specificity/weight ranking and first-occurrence union order.

    Select top genes in the complete reference BEFORE intersecting with query
    genes. Ties follow the original pandas two-column sort; no alphabetical
    tiebreaker or query-detection filter is introduced.
    """
    gene_total = ref.sum(axis=0).replace(0, np.nan)
    specificity = ref.div(gene_total, axis=1).fillna(0)
    selected = []
    for program in ref.index:
        ranked = pd.DataFrame({
            "gene": ref.columns.astype(str),
            "specificity": specificity.loc[program].to_numpy(),
            "weight": ref.loc[program].to_numpy(),
        }).sort_values(["specificity", "weight"], ascending=False)
        selected.extend(ranked.head(n_top)["gene"].tolist())
    query_set = set(pd.Index(query_gene_ids).astype(str))
    return [gene for gene in pd.Index(selected).drop_duplicates() if gene in query_set]


def build_collapsed_query_matrix(X_counts, query_gene_ids, selected_ref_genes, total_counts):
    """Original full-library CPM, duplicate summation, and float32 conversion."""
    X_counts = sparse.csr_matrix(X_counts)
    selected_ref_genes = pd.Index([str(gene) for gene in selected_ref_genes])
    query_gene_ids = pd.Index(query_gene_ids).astype(str)
    gene_to_positions = {}
    for position, gene in enumerate(query_gene_ids):
        if gene in selected_ref_genes:
            gene_to_positions.setdefault(gene, []).append(position)
    present = [gene for gene in selected_ref_genes if gene in gene_to_positions]
    if not present:
        raise ValueError("No selected reference genes are present in the query")
    query_positions, target_positions = [], []
    for target, gene in enumerate(present):
        for position in gene_to_positions[gene]:
            query_positions.append(position)
            target_positions.append(target)
    X_sub = X_counts[:, query_positions].tocsr()
    collapse = sparse.csr_matrix(
        (np.ones(len(query_positions), dtype=np.float32),
         (np.arange(len(query_positions)), target_positions)),
        shape=(len(query_positions), len(present)),
    )
    X_collapsed = X_sub @ collapse
    total_counts = np.asarray(total_counts, dtype=np.float64)
    if total_counts.shape != (X_counts.shape[0],):
        raise ValueError("total_counts does not match query rows")
    if not np.isfinite(total_counts).all() or (total_counts <= 0).any():
        raise ValueError("Full-library total_counts must be finite and positive")
    scale = np.divide(TARGET_SUM, total_counts,
                      out=np.zeros_like(total_counts, dtype=np.float64),
                      where=total_counts > 0)
    X_cpm = X_collapsed.multiply(scale[:, None]).tocsr().astype(np.float32)
    return X_cpm, present


def project_fixed_spectra(X_cpm, H_program_by_gene):
    """Original fixed-H multiplicative-update solve, without altered stopping."""
    X_input = (X_cpm.astype(np.float32) if sparse.issparse(X_cpm)
               else np.asarray(X_cpm, dtype=np.float32, order="C"))
    H = np.asarray(H_program_by_gene, dtype=np.float32, order="C").copy()
    _finite_nonnegative(H, "Selected spectra")
    sums = H.sum(axis=1, keepdims=True)
    if (sums <= 0).any():
        raise ValueError("A program has zero selected reference weight")
    H = np.divide(H, np.maximum(sums, np.asarray(1e-12, dtype=np.float32)),
                  out=np.zeros_like(H, dtype=np.float32), where=sums > 0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        usage, _, n_iter = non_negative_factorization(
            X_input, H=H, n_components=H.shape[0], update_H=False,
            init="custom", solver="mu", beta_loss="frobenius",
            max_iter=SKLEARN_MAX_ITER, tol=SKLEARN_TOL,
            random_state=RANDOM_STATE, verbose=0,
        )
    usage = np.asarray(usage, dtype=np.float32)
    _finite_nonnegative(usage, "Fitted usage")
    row_sums = usage.sum(axis=1, keepdims=True)
    if (row_sums <= 0).any():
        raise ValueError("Some cells have zero fitted program usage; CLR is undefined")
    normalized = np.divide(usage, row_sums, out=np.zeros_like(usage), where=row_sums > 0)
    diagnostics = {
        "n_iter": int(n_iter), "max_iter": SKLEARN_MAX_ITER,
        "hit_iteration_limit": bool(n_iter >= SKLEARN_MAX_ITER),
        "convergence_warning": any(issubclass(item.category, ConvergenceWarning) for item in caught),
        "warnings": [str(item.message) for item in caught],
        "raw_usage_sum_min": float(row_sums.min()),
        "raw_usage_sum_max": float(row_sums.max()),
    }
    for item in caught:
        warnings.warn(str(item.message), item.category, stacklevel=2)
    return usage.astype(np.float32), normalized.astype(np.float32), diagnostics


def clr_transform_usage(usage, eps=CLR_EPS, log_base=LOG_BASE):
    """Original normalized → floor → renormalized → centered log2 transform."""
    usage = np.asarray(usage, dtype=np.float64)
    _finite_nonnegative(usage, "CLR input usage")
    row_sums = usage.sum(axis=1, keepdims=True)
    if np.any(row_sums <= 0):
        raise ValueError("Some cells have zero total usage")
    usage = usage / row_sums
    usage = np.maximum(usage, eps)
    usage = usage / usage.sum(axis=1, keepdims=True)
    if log_base == 2:
        log_usage = np.log2(usage)
    elif log_base == "natural":
        log_usage = np.log(usage)
    else:
        raise ValueError("log_base must be 2 or 'natural'")
    return (log_usage - log_usage.mean(axis=1, keepdims=True)).astype(np.float32)


def add_projection_metadata(obs):
    """Annotate original condition-specific QC without discarding any cells."""
    obs = obs.copy()
    required = {"sample_name", "pct_counts_mt", "pct_counts_ribo"}
    missing = required.difference(obs.columns)
    if missing:
        raise ValueError(f"Query is missing projection metadata: {sorted(missing)}")
    parsed = obs["sample_name"].astype(str).str.extract(r"^(.*)_rep(\d+)$")
    if parsed.isna().any().any():
        raise ValueError("Every sample_name must end in _rep<number>")
    obs["condition"] = parsed[0].to_numpy()
    obs["replicate"] = parsed[1].astype(int).to_numpy()
    obs["treatment"] = obs["condition"].map(
        lambda condition: "TDI" if "TDI" in condition else "IAG" if "IAG" in condition else "control"
    )
    unknown = set(obs["condition"]) - set(COND_PCT_MT_CUTOFFS)
    if unknown:
        raise ValueError(f"Missing original condition-specific QC cutoffs: {sorted(unknown)}")
    for column in ("pct_counts_mt", "pct_counts_ribo"):
        obs[column] = pd.to_numeric(obs[column], errors="raise")
        if not np.isfinite(obs[column].to_numpy()).all():
            raise ValueError(f"Invalid values in {column}")
    mt_cutoff = obs["condition"].map(COND_PCT_MT_CUTOFFS).to_numpy()
    ribo_cutoff = obs["condition"].map(COND_PCT_RIBO_CUTOFFS).to_numpy()
    obs["qc_pass_condition_specific"] = (
        (obs["pct_counts_mt"].to_numpy() < mt_cutoff)
        & (obs["pct_counts_ribo"].to_numpy() > ribo_cutoff)
    )
    return obs


def _read_tables(path):
    """Read tables only: validation must not reload all query count layers."""
    import h5py
    try:
        from anndata.io import read_elem
    except ImportError:  # Older production anndata releases expose the same reader here.
        from anndata.experimental import read_elem
    with h5py.File(path, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
    return obs, var


def validate_projected(query_h5ad, projected_h5ad, *, strict_historical=True):
    """Validate complete cell/metadata correspondence, program order, and CLR."""
    import anndata as ad
    obs, query_var = _read_tables(query_h5ad)
    projected = ad.read_h5ad(projected_h5ad)
    if not obs.index.is_unique or not projected.obs_names.is_unique:
        raise ValueError("Query and projected cell IDs must be unique")
    if not projected.obs_names.equals(obs.index):
        raise ValueError("Projection must preserve every query cell in the exact query order")
    programs = [f"P{i}" for i in range(projected.n_vars)]
    if list(projected.var_names) != programs:
        raise ValueError("Projected program columns are not ordered P0..P(K-1)")
    for column in ("program", "program_id_zero_based"):
        if column not in projected.var:
            raise ValueError(f"Projected var lacks {column}")
    if list(projected.var["program"].astype(str)) != programs:
        raise ValueError("Projected var program labels disagree with var_names")
    if list(projected.var["program_id_zero_based"]) != list(range(projected.n_vars)):
        raise ValueError("Projected zero-based program IDs disagree with column order")
    declared = list(projected.uns.get("cnmf_program_names", []))
    if declared != programs:
        raise ValueError("Missing or inconsistent explicit obsm program-column labels")
    expected_obs = add_projection_metadata(obs)
    for column in expected_obs.columns:
        if column not in projected.obs:
            raise ValueError(f"Projected obs lost metadata column {column}")
        try:
            pd.testing.assert_series_equal(
                expected_obs[column], projected.obs[column], check_dtype=False,
                check_categorical=False, check_names=False,
            )
        except AssertionError as error:
            raise ValueError(f"Query/projected metadata disagree for {column}") from error
    for key in (RAW_KEY, USAGE_KEY, CLR_KEY):
        if key not in projected.obsm:
            raise ValueError(f"Missing projected matrix: {key}")
        matrix = np.asarray(projected.obsm[key])
        if matrix.shape != projected.shape or not np.isfinite(matrix).all():
            raise ValueError(f"Invalid dimensions or nonfinite values in {key}")
    raw = np.asarray(projected.obsm[RAW_KEY])
    normalized = np.asarray(projected.obsm[USAGE_KEY])
    _finite_nonnegative(raw, "Projected raw usage")
    _finite_nonnegative(normalized, "Projected normalized usage")
    if (raw.sum(axis=1) <= 0).any():
        raise ValueError("Projected raw usage has zero rows")
    expected_normalized = raw / raw.sum(axis=1, keepdims=True)
    if not np.allclose(normalized, expected_normalized, rtol=2e-6, atol=1e-7):
        raise ValueError("Normalized usage does not correspond to raw usage")
    if not np.allclose(normalized.sum(axis=1), 1.0, rtol=1e-6, atol=1e-6):
        raise ValueError("Normalized usage row sums differ from one")
    if not np.allclose(projected.obsm[CLR_KEY], clr_transform_usage(normalized), rtol=1e-6, atol=2e-6):
        raise ValueError("Selected CLR matrix does not match the original CLR transform")
    if projected.uns.get("X_projection_mode") != MODE or not np.array_equal(projected.X, normalized):
        raise ValueError("Rebuild X must explicitly contain normalized top_specific_200 usage")
    if list(projected.uns.get("projection_modes", [])) != [MODE]:
        raise ValueError("Rebuild must declare its single selected projection mode")
    n_qc = int(projected.obs["qc_pass_condition_specific"].sum())
    genes = list(projected.uns.get("genes_used_for_projection", []))
    if not genes or len(genes) != len(set(genes)):
        raise ValueError("Missing or duplicate ordered projection gene list")
    if strict_historical:
        if (len(obs), len(query_var)) != EXPECTED_QUERY_SHAPE:
            raise ValueError(f"Query shape differs from historical {EXPECTED_QUERY_SHAPE}")
        if projected.shape != (EXPECTED_QUERY_SHAPE[0], 50):
            raise ValueError("Projected shape differs from historical 212891 x 50")
        if len(genes) != EXPECTED_SELECTED_GENES or n_qc != EXPECTED_QC_PASS:
            raise ValueError(f"Historical gene/QC counts differ: genes={len(genes)}, QC={n_qc}")
    return {
        "query_shape": [len(obs), len(query_var)], "projected_shape": list(projected.shape),
        "exact_cell_order": True, "program_names": programs, "n_genes": len(genes),
        "qc_passing_cells": n_qc, "selected_clr_key": CLR_KEY,
        "metadata_matches": True, "clr_matches_original_formula": True,
    }


def _atomic_json(path, payload):
    path = Path(path)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                     prefix=f".{path.name}.", suffix=".tmp", delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def project_top_specific_200(query_h5ad, spectra_path, output_h5ad, *, overwrite=False,
                             strict_historical=True, compression=None):
    """Build/validate the selected artifact, or resume a matching verified build.

    Large query files are fingerprinted by path, size and nanosecond mtime (not
    misreported as content hashes). The reference and implementation are hashed.
    Unknown existing outputs are never silently reused or overwritten.
    """
    import anndata as ad
    query_h5ad, spectra_path, output_h5ad = map(Path, (query_h5ad, spectra_path, output_h5ad))
    if query_h5ad.resolve() == output_h5ad.resolve():
        raise ValueError("Query and projection output must be different files")
    for path in (query_h5ad, spectra_path):
        if not path.is_file():
            raise FileNotFoundError(f"Required projection input does not exist: {path}")
    receipt_path = output_h5ad.with_suffix(".rebuild_receipt.json")
    inputs = {"query": _signature(query_h5ad), "spectra": _signature(spectra_path, hash_content=True)}
    implementation_sha256 = _sha256(__file__)
    if output_h5ad.exists() and not overwrite:
        if not receipt_path.is_file():
            raise FileExistsError(f"Existing output lacks a rebuild receipt: {output_h5ad}. Set overwrite=True to rebuild it.")
        receipt = json.loads(receipt_path.read_text())
        if (receipt.get("inputs") != inputs
                or receipt.get("implementation_sha256") != implementation_sha256
                or receipt.get("output") != _signature(output_h5ad)):
            raise FileExistsError("Existing projection provenance does not match current inputs/code; set overwrite=True to rebuild")
        receipt["validation"] = validate_projected(query_h5ad, output_h5ad, strict_historical=strict_historical)
        receipt["reused"] = True
        print(f"Reusing validated {MODE}: {output_h5ad}")
        return receipt
    ref = read_spectra(spectra_path, strict_historical=strict_historical)
    print(f"Loading raw query X: {query_h5ad}", flush=True)
    query = ad.read_h5ad(query_h5ad)
    if not query.obs_names.is_unique or not query.var_names.is_unique:
        raise ValueError("Query cell and feature IDs must be unique")
    if "gene" not in query.var:
        raise ValueError("Query lacks var['gene'] for reference gene matching")
    if strict_historical and query.shape != EXPECTED_QUERY_SHAPE:
        raise ValueError(f"Query shape {query.shape}; expected {EXPECTED_QUERY_SHAPE}")
    obs = add_projection_metadata(query.obs)
    gene_ids = pd.Index(query.var["gene"].astype(str).to_numpy())
    X_counts = sparse.csr_matrix(query.X)
    _finite_nonnegative(X_counts, "Query X")
    if "total_counts" in obs:
        total_counts = pd.to_numeric(obs["total_counts"], errors="raise").to_numpy(dtype=float)
    else:
        total_counts = np.asarray(X_counts.sum(axis=1)).ravel().astype(float)
    del query
    gc.collect()
    genes = top_specific_genes(ref, gene_ids)
    if strict_historical and len(genes) != EXPECTED_SELECTED_GENES:
        raise ValueError(f"Selected {len(genes)} query genes; original selected {EXPECTED_SELECTED_GENES}")
    X_cpm, genes = build_collapsed_query_matrix(X_counts, gene_ids, genes, total_counts)
    del X_counts
    gc.collect()
    print(f"Projecting {X_cpm.shape[0]:,} cells using {len(genes):,} genes and {len(ref):,} fixed programs", flush=True)
    raw, normalized, diagnostics = project_fixed_spectra(X_cpm, ref.loc[:, genes].to_numpy())
    del X_cpm
    clr = clr_transform_usage(normalized)
    programs = list(ref.index)
    var = pd.DataFrame(index=pd.Index(programs, name="program"))
    var["program"] = programs
    var["program_id_zero_based"] = list(range(len(programs)))
    projected = ad.AnnData(X=normalized, obs=obs, var=var)
    projected.obsm[USAGE_KEY] = normalized
    projected.obsm[RAW_KEY] = raw
    projected.obsm[CLR_KEY] = clr
    projected.uns.update({
        "projection_modes": [MODE], "X_projection_mode": MODE,
        "X_semantics": "row-normalized top_specific_200 fixed-reference usage",
        "cnmf_program_names": programs, "reference_spectra": str(spectra_path.resolve()),
        "query_h5ad": str(query_h5ad.resolve()), "genes_used_for_projection": genes,
        "projection_gene_set": MODE, "n_top_specific_per_program": TOP_SPECIFIC_N,
        "clr_eps": CLR_EPS, "clr_log_base": LOG_BASE,
        "rebuild_source_notebook_sha256": SOURCE_NOTEBOOK_SHA256,
        "reference_spectra_sha256": inputs["spectra"]["sha256"],
        "fixed_H_n_iter": diagnostics["n_iter"],
        "qc_cutoffs_pct_mt": COND_PCT_MT_CUTOFFS,
        "qc_cutoffs_pct_ribo": COND_PCT_RIBO_CUTOFFS,
    })
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=output_h5ad.parent, prefix=f".{output_h5ad.stem}.",
                                     suffix=".h5ad", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        projected.write_h5ad(temporary, compression=compression)
        del projected, raw, normalized, clr
        gc.collect()
        validation = validate_projected(query_h5ad, temporary, strict_historical=strict_historical)
        temporary.replace(output_h5ad)
    finally:
        temporary.unlink(missing_ok=True)
    receipt = {
        "stage": "top_specific_200_fixed_reference_projection", "reused": False,
        "inputs": inputs, "output": _signature(output_h5ad),
        "source_notebook_sha256": SOURCE_NOTEBOOK_SHA256,
        "source_cells_zero_based": [3, 5, 7, 9, 11, 13],
        "implementation_sha256": implementation_sha256, "versions": _versions(),
        "parameters": {"counts_source": "X", "query_gene_column": "gene",
                       "target_sum": TARGET_SUM, "n_top_specific_per_program": TOP_SPECIFIC_N,
                       "solver": "mu", "update_H": False, "loss": "frobenius",
                       "max_iter": SKLEARN_MAX_ITER, "tol": SKLEARN_TOL,
                       "random_state": RANDOM_STATE, "dtype": "float32",
                       "clr_eps": CLR_EPS, "clr_log_base": LOG_BASE},
        "diagnostics": diagnostics, "validation": validation,
        "historical_n_iter_for_comparison": 140,
        "strict_historical": strict_historical,
        "scope": "Selected mode only; X is selected normalized usage, all cells retained, no control subtraction",
    }
    _atomic_json(receipt_path, receipt)
    print(f"Wrote and validated: {output_h5ad} (n_iter={diagnostics['n_iter']})", flush=True)
    return receipt
