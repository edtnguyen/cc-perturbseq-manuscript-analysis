"""Reconstruct the expression input used by the Hippo manuscript panels.

The production source is ``Hippo_timecourse.ipynb``, uploaded on 2026-10-29
with SHA-256 ``df24b42905e800ffe501ae8e1d2c334c40c684c2306eab36a44dc9f869565a3c``.
Zero-based source cells 5, 14, 20, 23, and 27 establish the concatenation,
broad QC, sample mapping, and write. The later condition-specific QC is NOT
part of this object. This module does not implement the historical all-gene
projection; the manuscript projection is a separate producer.

The same five producer cells are preserved verbatim in the existing repository
archive ``recovery/raw_downloads/Hippo_timecourse.ipynb``. Surrounding analysis
cells differ between that archive and the later uploaded notebook; the source
hash above identifies the uploaded version, not the entire archived notebook.

Input files are read only. An existing result is reused only with a matching
receipt, unchanged input/output file signatures, and successful validation.
"""
from __future__ import annotations

import gc
import hashlib
import importlib.metadata
import json
import os
import platform
import tempfile
from datetime import datetime, timezone
from pathlib import Path


SOURCE_NOTEBOOK_SHA256 = "df24b42905e800ffe501ae8e1d2c334c40c684c2306eab36a44dc9f869565a3c"
SOURCE_CELLS_ZERO_BASED = [5, 14, 20, 23, 27]
REPOSITORY_SOURCE_NOTEBOOK = "recovery/raw_downloads/Hippo_timecourse.ipynb"
REBUILD_VERSION = "hippo-combined-v1"
EXPECTED_PRE_QC_SHAPE = (2_118_711, 62_757)
EXPECTED_POST_QC_SHAPE = (212_891, 62_757)
REQUIRED_LAYERS = ("ambiguous", "mature", "nascent")
REQUIRED_OBSM = ("guide_assignment", "guide_counts")
REQUIRED_INPUT_OBS = (
    "n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo",
    "biological_sample", "plate", "well",
)
INPUT_RELATIVE_TEMPLATE = "Hippo:eTN8_GEX_sub{i}/kb_all_main_raw/counts_filtered/adata.h5ad"

SAMPLE_MAPPING = {
    "1-1": "D0_hPSC_rep1", "1-2": "D0_hPSC_rep2", "1-3": "D0_hPSC_rep3",
    "1-4": "D1_MPS_rep1", "1-5": "D1_MPS_rep2", "1-6": "D1_MPS_rep3",
    "1-7": "D1_MPS+TDI_rep1", "1-8": "D1_MPS+TDI_rep2", "1-9": "D1_MPS+TDI_rep3",
    "1-10": "D1_MPS+IAG_rep1", "1-11": "D1_MPS+IAG_rep2", "1-12": "D1_MPS+IAG_rep3",
    "2-1": "D1.33_LM_rep1", "2-2": "D1.33_LM_rep2", "2-3": "D1.33_LM_rep3",
    "2-4": "D1.33_LM+TDI_rep1", "2-5": "D1.33_LM+TDI_rep2", "2-6": "D1.33_LM+TDI_rep3",
    "2-7": "D1.33_LM+IAG_rep1", "2-8": "D1.33_LM+IAG_rep2", "2-9": "D1.33_LM+IAG_rep3",
    "3-1": "D1.67_LM_rep1", "3-2": "D1.67_LM_rep2", "3-3": "D1.67_LM_rep3",
    "3-4": "D1.67_LM+TDI_rep1", "3-5": "D1.67_LM+TDI_rep2", "3-6": "D1.67_LM+TDI_rep3",
    "3-7": "D1.67_LM+IAG_rep1", "3-8": "D1.67_LM+IAG_rep2", "3-9": "D1.67_LM+IAG_rep3",
    "4-1": "D2_LM_rep1", "4-2": "D2_LM_rep2", "4-3": "D2_LM_rep3",
    "4-4": "D2_LM+TDI_rep1", "4-5": "D2_LM+TDI_rep2", "4-6": "D2_LM+TDI_rep3",
    "4-7": "D2_LM+IAG_rep1", "4-8": "D2_LM+IAG_rep2", "4-9": "D2_LM+IAG_rep3",
    "5-1": "D3_AE_rep1", "5-2": "D3_AE_rep2", "5-3": "D3_AE_rep3",
    "5-4": "D3_AE+TDI_rep1", "5-5": "D3_AE+TDI_rep2", "5-6": "D3_AE+TDI_rep3",
    "5-7": "D3_AE+IAG_rep1", "5-8": "D3_AE+IAG_rep2", "5-9": "D3_AE+IAG_rep3",
    "5-10": "D2_LM+TDI->D3_AE+TDI_rep1", "5-11": "D2_LM+TDI->D3_AE+TDI_rep2",
    "5-12": "D2_LM+TDI->D3_AE+TDI_rep3", "5-13": "D2_LM+IAG->D3_AE+IAG_rep1",
    "5-14": "D2_LM+IAG->D3_AE+IAG_rep2", "5-15": "D2_LM+IAG->D3_AE+IAG_rep3",
}

# Plain-text output of production cell 20, independently of the new pipeline.
HISTORICAL_SAMPLE_COUNTS = {
    "4-1": 8103, "3-8": 7808, "3-7": 7512, "4-3": 6870, "4-2": 6769,
    "4-7": 6517, "4-8": 6376, "1-1": 6239, "1-2": 6142, "4-5": 6095,
    "4-6": 5960, "5-10": 5661, "5-9": 5359, "5-3": 4906, "5-11": 4833,
    "2-7": 4747, "1-11": 4657, "5-12": 4646, "1-12": 4593, "2-1": 4397,
    "1-3": 4389, "2-2": 4149, "3-9": 4128, "2-3": 4002, "1-9": 3984,
    "4-4": 3969, "5-6": 3967, "2-6": 3916, "1-10": 3881, "2-8": 3832,
    "2-5": 3444, "5-4": 3439, "1-8": 3398, "3-6": 3274, "1-5": 2936,
    "2-4": 2885, "1-7": 2698, "5-2": 2626, "4-9": 2597, "1-6": 2573,
    "5-14": 2485, "2-9": 2355, "3-2": 2311, "3-5": 2259, "1-4": 2075,
    "5-5": 1975, "3-3": 1935, "5-7": 1898, "5-13": 1790, "5-8": 1312,
    "5-1": 1172, "3-1": 1089, "5-15": 1032, "3-4": 926,
}


def combined_receipt_path(path: Path) -> Path:
    """Return the adjacent JSON receipt path for a combined H5AD."""
    path = Path(path)
    return path.with_suffix(path.suffix + ".rebuild_receipt.json")


def _file_record(path: Path) -> dict:
    stat = path.stat()
    return {"path": str(path.resolve()), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def preflight_sublibraries(input_root: Path) -> list[dict]:
    """Check all eight exact processed inputs before loading or writing data."""
    root = Path(input_root).expanduser().resolve()
    paths = [root / INPUT_RELATIVE_TEMPLATE.format(i=i) for i in range(1, 9)]
    missing = [str(p) for p in paths if not p.is_file() or p.stat().st_size == 0]
    if missing:
        raise FileNotFoundError(
            "Missing or empty processed Hippo sublibrary inputs:\n  "
            + "\n  ".join(missing)
            + "\nSet INPUT_ROOT to the directory containing the eight Hippo:eTN8_GEX_sub* "
              "folders. Rebuilding requires these processed inputs; it does not rerun quantification."
        )
    if len({p.resolve() for p in paths}) != 8:
        raise ValueError("The eight sublibrary paths must resolve to eight distinct files.")
    return [{"sublibID": f"s{i}", **_file_record(p)} for i, p in enumerate(paths, start=1)]


def _signature(records: list[dict]) -> str:
    canonical = json.dumps(records, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def _versions() -> dict:
    versions = {"python": platform.python_version()}
    for package in ("anndata", "numpy", "pandas", "scipy", "h5py"):
        versions[package] = importlib.metadata.version(package)
    return versions


def _atomic_json(path: Path, content: dict) -> None:
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
                                         prefix=f".{path.name}.", suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(content, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _numeric_qc(obs, *, require_finite=False):
    import numpy as np
    import pandas as pd

    missing = sorted(set(REQUIRED_INPUT_OBS) - set(obs.columns))
    if missing:
        raise ValueError(f"Missing Hippo observation metadata: {missing}")
    numeric = {}
    for name in ("n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"):
        values = pd.to_numeric(obs[name], errors="raise")
        if require_finite and not np.isfinite(values.to_numpy(dtype=float, na_value=np.nan)).all():
            raise ValueError(f"QC metadata {name!r} contains missing or non-finite values.")
        numeric[name] = values
    if require_finite and obs["biological_sample"].isna().any():
        raise ValueError("biological_sample contains missing values.")
    return numeric


def _broad_qc_mask(obs):
    values = _numeric_qc(obs)
    # Preserve the original ordering: compare 'empty' before stripping labels.
    return (
        (values["n_genes_by_counts"] < 12000)
        & (values["n_genes_by_counts"] > 1200)
        & (values["total_counts"] < 50000)
        & (values["total_counts"] > 1200)
        & (values["pct_counts_mt"] < 10)
        & (obs["biological_sample"] != "empty")
    ).fillna(False)


def _check_metadata(obs, var, *, strict_historical: bool) -> dict:
    import numpy as np

    if not obs.index.is_unique or not var.index.is_unique:
        raise ValueError("Cell IDs and feature IDs must each be unique.")
    if len(obs) == 0 or len(var) == 0:
        raise ValueError("The combined expression object must have cells and features.")
    if "gene" not in var:
        raise ValueError("Gene symbols are required in var['gene'].")
    if var["gene"].isna().any() or (var["gene"].astype(str).str.len() == 0).any():
        raise ValueError("var['gene'] contains missing or empty symbols.")
    missing = {"sample_name", "sublibID"} - set(obs.columns)
    if missing:
        raise ValueError(f"Missing combined-object metadata: {sorted(missing)}")
    _numeric_qc(obs, require_finite=True)
    samples = obs["biological_sample"].astype(str)
    if not np.array_equal(samples.to_numpy(), samples.str.strip().to_numpy()):
        raise ValueError("Combined biological_sample values must be stripped, as in source cell 23.")
    expected_names = samples.map(SAMPLE_MAPPING)
    if expected_names.isna().any():
        unknown = sorted(samples[expected_names.isna()].unique())
        raise ValueError(f"Unmapped biological_sample values: {unknown}")
    if not np.array_equal(expected_names.to_numpy(), obs["sample_name"].astype(str).to_numpy()):
        raise ValueError("sample_name does not match the production 54-sample mapping.")
    if not _broad_qc_mask(obs).all():
        raise ValueError("The combined object contains cells failing the production broad QC.")
    sublibraries = obs["sublibID"].astype(str)
    if not set(sublibraries).issubset({f"s{i}" for i in range(1, 9)}):
        raise ValueError("sublibID must use the production s1...s8 labels.")
    suffixes = obs.index.astype(str).str.rsplit("_", n=1).str[-1]
    if not np.array_equal(np.asarray(suffixes), sublibraries.to_numpy()):
        raise ValueError("Cell IDs do not preserve the original _sN sublibrary suffixes.")
    counts = {str(k): int(v) for k, v in samples.value_counts().items()}
    shape = (len(obs), len(var))
    if strict_historical:
        if shape != EXPECTED_POST_QC_SHAPE:
            raise ValueError(f"Combined shape {shape} differs from recorded {EXPECTED_POST_QC_SHAPE}.")
        if counts != HISTORICAL_SAMPLE_COUNTS:
            keys = sorted(set(counts) | set(HISTORICAL_SAMPLE_COUNTS))
            different = {k: {"observed": counts.get(k, 0), "recorded": HISTORICAL_SAMPLE_COUNTS.get(k, 0)}
                         for k in keys if counts.get(k, 0) != HISTORICAL_SAMPLE_COUNTS.get(k, 0)}
            raise ValueError(f"Per-sample cell counts differ from production cell 20: {different}")
    return {
        "shape": list(shape), "n_samples": len(counts), "sample_counts": counts,
        "sublibrary_counts": {str(k): int(v) for k, v in sublibraries.value_counts().items()},
    }


def _h5_shape(element) -> tuple:
    if hasattr(element, "shape"):
        return tuple(element.shape)
    encoding = element.attrs.get("encoding-type", "")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if encoding == "dataframe":
        # obsm legitimately supports DataFrames. Read only index dimensions and
        # column metadata, without materializing a potentially large guide table.
        index_name = element.attrs.get("_index", "_index")
        if isinstance(index_name, bytes):
            index_name = index_name.decode()
        index = element[index_name]
        if hasattr(index, "shape"):
            n_rows = index.shape[0]
        elif "codes" in index:  # Categorical index encoding.
            n_rows = index["codes"].shape[0]
        elif "values" in index:  # Nullable string/integer index encoding.
            n_rows = index["values"].shape[0]
        else:
            raise ValueError(f"Unsupported DataFrame index encoding in {element.name}.")
        return (n_rows, len(element.attrs["column-order"]))
    return tuple(int(n) for n in element.attrs.get("shape", ()))


def _scan_count_matrix(element, expected_shape: tuple, *, row_totals: bool = False):
    """Validate dense/CSR/CSC HDF5 counts without loading an entire layer."""
    import h5py
    import numpy as np
    from scipy import sparse

    if _h5_shape(element) != expected_shape:
        raise ValueError(f"{element.name} shape {_h5_shape(element)} differs from {expected_shape}.")
    totals = np.zeros(expected_shape[0], dtype=np.float64) if row_totals else None

    def check(values):
        if not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError(f"{element.name} contains negative or non-finite raw counts.")

    if isinstance(element, h5py.Dataset):
        for start in range(0, expected_shape[0], 512):
            block = element[start:start + 512, :]
            check(block)
            if totals is not None:
                totals[start:start + len(block)] = block.sum(axis=1, dtype=np.float64)
        return totals

    encoding = element.attrs.get("encoding-type", "")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if encoding not in {"csr_matrix", "csc_matrix"}:
        raise ValueError(f"Unsupported raw-count encoding for {element.name}: {encoding!r}")
    data, indices = element["data"], element["indices"]
    pointer = element["indptr"][:]
    major = expected_shape[0 if encoding == "csr_matrix" else 1]
    minor = expected_shape[1 if encoding == "csr_matrix" else 0]
    if (len(pointer) != major + 1 or pointer[0] != 0 or pointer[-1] != len(data)
            or len(indices) != len(data) or np.any(np.diff(pointer) < 0)):
        raise ValueError(f"Invalid sparse matrix structure in {element.name}.")
    if row_totals and encoding == "csr_matrix":
        for start in range(0, major, 512):
            stop = min(start + 512, major)
            first, last = int(pointer[start]), int(pointer[stop])
            values, columns = data[first:last], indices[first:last]
            check(values)
            if np.any(columns < 0) or np.any(columns >= minor):
                raise ValueError(f"Out-of-range sparse indices in {element.name}.")
            block = sparse.csr_matrix((values, columns, pointer[start:stop + 1] - first),
                                      shape=(stop - start, minor))
            totals[start:stop] = np.asarray(block.sum(axis=1, dtype=np.float64)).ravel()
    else:
        for start in range(0, len(data), 1_000_000):
            values, rows = data[start:start + 1_000_000], indices[start:start + 1_000_000]
            check(values)
            if np.any(rows < 0) or np.any(rows >= minor):
                raise ValueError(f"Out-of-range sparse indices in {element.name}.")
            if totals is not None:  # CSC entries carry row indices.
                totals += np.bincount(rows, weights=values, minlength=expected_shape[0])
    return totals


def validate_combined(path: Path, *, strict_historical: bool = True) -> dict:
    """Validate source schema, QC, mapping, raw counts, and historical counts.

    Counts may be fractional. HDF5 matrices are inspected in bounded blocks;
    ``anndata.read_h5ad(backed='r')`` is avoided because it eagerly loads layers.
    ``strict_historical=False`` is intended for bounded fixtures or an explicitly
    documented revised source cohort, never a silent production fallback.
    """
    import h5py
    import numpy as np
    try:
        from anndata.io import read_elem
    except ImportError:  # Compatibility with the original older AnnData environment.
        from anndata._io.specs import read_elem

    path = Path(path)
    with h5py.File(path, "r") as handle:
        obs, var = read_elem(handle["obs"]), read_elem(handle["var"])
        report = _check_metadata(obs, var, strict_historical=strict_historical)
        shape = tuple(report["shape"])
        for group, required in (("layers", REQUIRED_LAYERS), ("obsm", REQUIRED_OBSM)):
            available = set(handle[group].keys()) if group in handle else set()
            missing = set(required) - available
            if missing:
                raise ValueError(f"Combined object lacks required {group}: {sorted(missing)}")
            report[group] = sorted(available)
        totals = _scan_count_matrix(handle["X"], shape, row_totals=True)
        for name in REQUIRED_LAYERS:
            _scan_count_matrix(handle["layers"][name], shape)
        for name in REQUIRED_OBSM:
            dimensions = _h5_shape(handle["obsm"][name])
            if len(dimensions) != 2 or dimensions[0] != shape[0]:
                raise ValueError(f"obsm[{name!r}] has incompatible cell dimensions: {dimensions}")
        recorded_totals = obs["total_counts"].to_numpy(dtype=float)
        difference = np.abs(totals - recorded_totals)
        report["max_abs_X_vs_total_counts"] = float(difference.max(initial=0))
        # Float32 summation/quantification can introduce small absolute differences.
        if not np.allclose(totals, recorded_totals, rtol=1e-4, atol=1e-3):
            raise ValueError(
                "X row sums do not match obs['total_counts']; confirm this is the original "
                "raw expression matrix, not normalized/logged expression. "
                f"Maximum absolute difference: {report['max_abs_X_vs_total_counts']:.6g}"
            )
        report["strict_historical"] = bool(strict_historical)
        report["raw_count_validation"] = "finite, nonnegative; X totals agree with source QC metadata"
    return report


def _reuse_existing(output: Path, inputs: list[dict], *, strict_historical: bool) -> dict:
    receipt_path = combined_receipt_path(output)
    reason = "receipt is missing or does not match the current inputs/output"
    try:
        receipt = json.loads(receipt_path.read_text())
        trusted = (
            receipt["producer_version"] == REBUILD_VERSION
            and receipt["source_notebook_sha256"] == SOURCE_NOTEBOOK_SHA256
            and receipt["input_signature"] == _signature(inputs)
            and receipt["output_file"] == _file_record(output)
        )
        if trusted:
            receipt["validation"] = validate_combined(output, strict_historical=strict_historical)
            receipt["status"] = "reused"
            print(f"Reusing validated combined expression object: {output}", flush=True)
            return receipt
    except (KeyError, ValueError, OSError) as exc:
        reason = str(exc)
    raise FileExistsError(
        f"Refusing to reuse existing output {output}: {reason}. "
        "Inspect it or explicitly set OVERWRITE=True to rebuild."
    )


def rebuild_combined(input_root: Path, output_h5ad: Path, *, overwrite: bool = False,
                     strict_historical: bool = True, compression=None) -> dict:
    """Rebuild the combined expression object with exact production broad QC.

    Matching completed outputs are validated and reused when ``overwrite=False``.
    Unknown or changed existing outputs are never silently trusted or overwritten.
    """
    import anndata as ad
    import pandas as pd

    inputs = preflight_sublibraries(input_root)
    output = Path(output_h5ad).expanduser().resolve()
    if str(output) in {record["path"] for record in inputs}:
        raise ValueError("The output must be separate from every read-only sublibrary input.")
    if compression not in (None, "gzip", "lzf"):
        raise ValueError("compression must be None, 'gzip', or 'lzf'.")
    if output.exists() and not overwrite:
        return _reuse_existing(output, inputs, strict_historical=strict_historical)

    adatas = {}
    input_shapes = {}
    for record in inputs:
        print(f"Loading {record['sublibID']}: {record['path']}", flush=True)
        current = ad.read_h5ad(record["path"])
        _numeric_qc(current.obs)
        if "gene" not in current.var:
            raise ValueError(f"{record['path']} lacks var['gene'].")
        if not current.obs_names.is_unique or not current.var_names.is_unique:
            raise ValueError(f"{record['path']} has duplicate cell or feature IDs.")
        for attr, required in (("layers", REQUIRED_LAYERS), ("obsm", REQUIRED_OBSM)):
            missing = set(required) - set(getattr(current, attr).keys())
            if missing:
                raise ValueError(f"{record['path']} lacks {attr}: {sorted(missing)}")
        input_shapes[record["sublibID"]] = list(current.shape)
        adatas[record["sublibID"]] = current
    del current
    print("Concatenating s1...s8 with the production outer join and barcode suffixes.", flush=True)
    combined = ad.concat(adatas, label="sublibID", index_unique="_", join="outer",
                         merge="same", uns_merge="same")
    adatas.clear()
    del adatas
    gc.collect()
    pre_qc_shape = tuple(combined.shape)
    if strict_historical and pre_qc_shape != EXPECTED_PRE_QC_SHAPE:
        raise ValueError(f"Pre-QC shape {pre_qc_shape} differs from recorded {EXPECTED_PRE_QC_SHAPE}.")
    keep = _broad_qc_mask(combined.obs)
    filtered = combined[keep.to_numpy(), :].copy()
    del combined, keep
    gc.collect()
    biological = filtered.obs["biological_sample"]
    if isinstance(biological.dtype, pd.CategoricalDtype):
        biological = biological.cat.remove_unused_categories()
    filtered.obs["biological_sample"] = biological.astype(str).str.strip()
    filtered.obs["sample_name"] = filtered.obs["biological_sample"].map(SAMPLE_MAPPING).astype("category")
    metadata = _check_metadata(filtered.obs, filtered.var, strict_historical=strict_historical)
    print(f"Broad QC: {pre_qc_shape[0]:,} -> {filtered.n_obs:,} cells; {filtered.n_vars:,} features.", flush=True)

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=output.parent, prefix=f".{output.name}.",
                                         suffix=".tmp.h5ad", delete=False) as handle:
            temporary = Path(handle.name)
        filtered.write_h5ad(temporary, compression=compression)
        del filtered
        gc.collect()
        validation = validate_combined(temporary, strict_historical=strict_historical)
        # Inputs must not have changed while they were being read.
        if preflight_sublibraries(input_root) != inputs:
            raise RuntimeError("An input file changed during reconstruction; the staged output was discarded.")
        temporary.replace(output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)

    receipt = {
        "status": "rebuilt", "producer_version": REBUILD_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_notebook": "Hippo_timecourse.ipynb",
        "source_notebook_sha256": SOURCE_NOTEBOOK_SHA256,
        "source_cells_zero_based": SOURCE_CELLS_ZERO_BASED,
        "repository_source_notebook": REPOSITORY_SOURCE_NOTEBOOK,
        "source_cells_match_repository_archive": True,
        "input_files": inputs, "input_shapes": input_shapes, "input_signature": _signature(inputs),
        "output_h5ad": str(output), "output_file": _file_record(output),
        "pre_qc_shape": list(pre_qc_shape), "post_qc_shape": metadata["shape"],
        "recorded_pre_qc_shape": list(EXPECTED_PRE_QC_SHAPE),
        "recorded_post_qc_shape": list(EXPECTED_POST_QC_SHAPE),
        "versions": _versions(), "validation": validation,
        "qc": {"n_genes_by_counts": "1200 < x < 12000", "total_counts": "1200 < x < 50000",
               "pct_counts_mt": "x < 10", "biological_sample": "not equal to 'empty' before stripping"},
        "condition_specific_qc_applied": False,
    }
    _atomic_json(combined_receipt_path(output), receipt)
    print(f"Wrote and validated: {output}", flush=True)
    return receipt
