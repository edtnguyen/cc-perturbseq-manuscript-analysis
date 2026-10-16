# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/cc_perturbseq_upload_ready_builder/notebooks/01_make_igvf_upload_ready_files.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/cc_perturbseq_upload_ready_builder/notebooks/01_make_igvf_upload_ready_files.ipynb`
- SHA-256: `dadda47ef40f434b479ad8bac8abf8ddf6c316c262497f41b0f4d6b2434b59ce`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0


# Build the 16 new CC-Perturb-seq analysis files for IGVF

This notebook **writes the final files intended for upload**, rather than only creating a staging plan.
It produces and validates:

- 4 day-specific annotated `.h5ad` files;
- 4 day-specific PerturbNMF `.h5mu` files;
- 4 day-specific `guide_quantifications.tsv.gz` files;
- 4 day-specific `per_cell_quality_report.tsv.gz` files;
- MD5 checksums, a file-level validation report, upload manifests, column dictionaries, and staged IGVF metadata TSVs.

## How day is resolved

The attached `iPSC_EC_QQ_full_covariates_with_replicate_correlation.ipynb` does **not** read a day column from the combined PerturbNMF h5mu. It loops over `d0`, `d1`, `d2`, and `d3` and opens separate files named `{day}_normalized_gex_and_guide.h5mu`.

For the combined `cNMF_50_0_2.h5mu`, this notebook therefore uses a stricter method:

1. read the combined annotated h5ad;
2. require an unambiguous day column, preferring `obs["day"]`;
3. align the h5ad and the PerturbNMF `cNMF` modality by exact cell IDs;
4. compare any day annotation already present in the h5mu against the h5ad;
5. stop on any mismatch, missing cell, unmapped day, or cross-timepoint contamination;
6. create true cell subsets, not four copies of the combined file.

No Portal POST, PATCH, or upload is performed by this notebook.


## Code cell 1

```python

from pathlib import Path

# -----------------------------
# Source files on Sherlock/OAK
# -----------------------------
SOURCE_H5AD = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
    "2.preprocessing_tonys_pipeline/3.integrated_days/"
    "3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad"
)

SOURCE_H5MU = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/"
    "adata/cNMF_50_0_2.h5mu"
)

# Final submission directory from the analysis plan.
OUTPUT_ROOT = Path(
    "/oak/stanford/groups/engreitz/Projects/IGVF/"
    "CC-Perturb-seq/300gene/analysis_submission"
)

# This file is bundled with the notebook package. Adjust PACKAGE_ROOT only if
# the notebook was moved out of the package directory.
_here = Path.cwd().resolve()
PACKAGE_ROOT = _here.parent if _here.name == "notebooks" else _here
SEQUENCE_ACCESSION_CSV = PACKAGE_ROOT / "inputs" / "sequence_file_accessions.csv"

# The actual Portal AnalysisSet accessions supplied by the user.
ANALYSIS_SETS = {
    "D0": "IGVFDS0463XTUH",
    "D1": "IGVFDS0991DBQZ",
    "D2": "IGVFDS4402ZGWJ",
    "D3": "IGVFDS2307REVJ",
}

DAYS = ("D0", "D1", "D2", "D3")
REFERENCE_FILE = "/reference-files/IGVFFI7217ZMJZ/"
AWARD = "/awards/HG011972/"
LAB = "/labs/jesse-engreitz/"

# Safety controls. EXECUTE_WRITES=True creates the final files.
EXECUTE_WRITES = True
OVERWRITE_EXISTING = False
REQUIRE_EXACT_H5AD_H5MU_CELL_SET = True
H5MU_ANCHOR_MODALITY = "cNMF"

WORK_DIR = OUTPUT_ROOT / "_prep_work"
METADATA_DIR = OUTPUT_ROOT / "metadata"
DOCUMENTATION_DIR = OUTPUT_ROOT / "documentation"

print("SOURCE_H5AD:", SOURCE_H5AD)
print("SOURCE_H5MU:", SOURCE_H5MU)
print("OUTPUT_ROOT:", OUTPUT_ROOT)
print("SEQUENCE_ACCESSION_CSV:", SEQUENCE_ACCESSION_CSV)
print("EXECUTE_WRITES:", EXECUTE_WRITES)
print("OVERWRITE_EXISTING:", OVERWRITE_EXISTING)

```

### Text output (stream)

```text
SOURCE_H5AD: /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad
SOURCE_H5MU: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/adata/cNMF_50_0_2.h5mu
OUTPUT_ROOT: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission
SEQUENCE_ACCESSION_CSV: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/cc_perturbseq_upload_ready_builder/inputs/sequence_file_accessions.csv
EXECUTE_WRITES: True
OVERWRITE_EXISTING: False

```

## Code cell 2

```python

# Imports and environment checks. The notebook deliberately uses only plain-text
# output so it does not depend on the VS Code notebook renderer.
import csv
import gc
import gzip
import hashlib
import io
import json
import os
import platform
import re
import shutil
import tempfile
import time
from datetime import datetime, timezone
from typing import Any, Iterable, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

try:
    import anndata as ad
except ImportError as exc:
    raise ImportError(
        "anndata is required. Create/activate the environment described in README.md."
    ) from exc

try:
    import mudata as md
except ImportError:
    try:
        import muon as md
    except ImportError as exc:
        raise ImportError(
            "mudata or muon is required. Create/activate the environment described in README.md."
        ) from exc

print("Python:", platform.python_version())
print("anndata:", getattr(ad, "__version__", "unknown"))
print("mudata/muon:", getattr(md, "__version__", "unknown"))
print("pandas:", pd.__version__)
print("numpy:", np.__version__)
print("h5py:", h5py.__version__)

for required in (SOURCE_H5AD, SOURCE_H5MU, SEQUENCE_ACCESSION_CSV):
    if not required.is_file():
        raise FileNotFoundError(required)

if set(ANALYSIS_SETS) != set(DAYS) or len(set(ANALYSIS_SETS.values())) != 4:
    raise ValueError("ANALYSIS_SETS must contain four unique D0-D3 Portal accessions.")

for directory in (OUTPUT_ROOT, WORK_DIR, METADATA_DIR, DOCUMENTATION_DIR):
    directory.mkdir(parents=True, exist_ok=True)

print("Input paths exist and output directories are writable.")

```

### Text output (stream)

```text
Python: 3.11.14
anndata: 0.12.6
mudata/muon: 0.3.2
pandas: 2.3.3
numpy: 2.3.5
h5py: 3.15.1
Input paths exist and output directories are writable.

```

### Text output (stream)

```text
/tmp/ipykernel_8316/3517122031.py:41: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('anndata')` instead.
  print("anndata:", getattr(ad, "__version__", "unknown"))

```

## Code cell 3

```python

DAY_CANDIDATES = (
    "day",
    "Day",
    "timepoint",
    "time_point",
    "differentiation_day",
    "collection_day",
)

QC_COLUMNS = (
    "biological_sample",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "guides_per_cell",
    "guide_umi_counts",
)

GUIDE_ASSIGNMENT_KEYS = (
    "guide_assignment",
    "guide_assignments",
)

GUIDE_UMI_MATRIX_KEYS = (
    "guide_umi_matrix",
    "guide_counts",
    "guide_umi_counts",
)


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def normalize_day(value: Any) -> str | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    text = str(value).strip().upper()
    text = text.replace("DAY", "D")
    text = re.sub(r"[^A-Z0-9]", "", text)
    return {
        "0": "D0", "D0": "D0",
        "1": "D1", "D1": "D1",
        "2": "D2", "D2": "D2",
        "3": "D3", "D3": "D3",
    }.get(text)


def normalize_day_series(values: pd.Series, *, label: str) -> pd.Series:
    normalized = values.map(normalize_day)
    if normalized.isna().any():
        raw = sorted(map(str, values.loc[normalized.isna()].dropna().unique()))
        raise ValueError(f"{label}: unmapped day values: {raw[:30]}")
    observed = set(normalized.unique())
    missing = set(DAYS) - observed
    extra = observed - set(DAYS)
    if missing or extra:
        raise ValueError(
            f"{label}: expected exactly D0-D3; missing={sorted(missing)}, extra={sorted(extra)}"
        )
    return normalized.astype(str)


def find_day_column(obs: pd.DataFrame, *, label: str) -> str:
    if "day" in obs.columns:
        return "day"
    exact = [column for column in DAY_CANDIDATES if column in obs.columns]
    if len(exact) == 1:
        return exact[0]
    if len(exact) > 1:
        valid = []
        for column in exact:
            mapped = obs[column].map(normalize_day)
            if mapped.notna().all() and set(mapped.unique()) == set(DAYS):
                valid.append(column)
        if len(valid) == 1:
            return valid[0]
        raise ValueError(f"{label}: ambiguous day columns: {exact}; valid={valid}")

    fuzzy = [
        str(column)
        for column in obs.columns
        if re.search(r"(^|_)(day|timepoint|time_point)($|_)", str(column), re.I)
    ]
    valid = []
    for column in fuzzy:
        mapped = obs[column].map(normalize_day)
        if mapped.notna().all() and set(mapped.unique()) == set(DAYS):
            valid.append(column)
    if len(valid) == 1:
        return valid[0]
    raise KeyError(
        f"{label}: no unambiguous D0-D3 day column found. "
        f"Candidates={fuzzy}; obs columns={list(map(str, obs.columns))}"
    )


def md5sum(path: str | Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_replace(temp_path: Path, final_path: Path) -> None:
    final_path.parent.mkdir(parents=True, exist_ok=True)
    os.replace(temp_path, final_path)


def write_deterministic_tsv_gz(df: pd.DataFrame, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(
        prefix=output_path.name + ".",
        suffix=".tmp",
        dir=output_path.parent,
    )
    os.close(fd)
    temp_path = Path(temp_name)
    try:
        with temp_path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
                    df.to_csv(
                        text,
                        sep="\t",
                        index=False,
                        na_rep="NA",
                        lineterminator="\n",
                    )
        atomic_replace(temp_path, output_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()
    return output_path


def find_modality(mdata: Any, requested: str) -> str:
    if requested in mdata.mod:
        return requested
    matches = [name for name in mdata.mod if str(name).lower() == requested.lower()]
    if len(matches) == 1:
        return matches[0]
    raise KeyError(
        f"Could not find modality {requested!r}. Available modalities: {list(mdata.mod.keys())}"
    )


def read_h5mu(path: Path, *, backed: str | bool | None = None):
    if hasattr(md, "read_h5mu"):
        return md.read_h5mu(path, backed=backed)
    return md.read(path, backed=backed)


def write_h5mu(mdata: Any, path: Path) -> None:
    if hasattr(mdata, "write_h5mu"):
        mdata.write_h5mu(path)
    else:
        mdata.write(path)


def close_backed_anndata(adata: Any) -> None:
    try:
        if getattr(adata, "file", None) is not None:
            adata.file.close()
    except Exception:
        pass


def close_backed_mudata(mdata: Any) -> None:
    # MuData versions expose backing handles differently. Best-effort closure only.
    for obj in [mdata, *list(getattr(mdata, "mod", {}).values())]:
        try:
            if getattr(obj, "file", None) is not None:
                obj.file.close()
        except Exception:
            pass


def validate_hdf5_container(path: Path) -> None:
    with h5py.File(path, "r") as handle:
        if not list(handle.keys()):
            raise ValueError(f"{path} is an empty HDF5 container")


def guide_names_from_adata(adata: Any, matrix: Any) -> list[str]:
    if isinstance(matrix, pd.DataFrame):
        names = list(map(str, matrix.columns))
    else:
        names = []
        for key in ("guide_names", "guide_ids", "guides", "gRNA_names"):
            if key in adata.uns:
                values = list(map(str, adata.uns[key]))
                if len(values) == matrix.shape[1]:
                    names = values
                    break
    if len(names) != matrix.shape[1]:
        raise ValueError(
            "A guide-assignment matrix was found, but guide names could not be resolved "
            "from DataFrame columns or adata.uns['guide_names'/'guide_ids'/'guides']."
        )
    if len(set(names)) != len(names):
        raise ValueError("Guide names are not unique.")
    return names


def dataframe_for_serialization(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for column in out.columns:
        if isinstance(out[column].dtype, pd.CategoricalDtype):
            out[column] = out[column].astype(str)
    return out


def extract_guide_quantifications(adata: Any, day: str) -> pd.DataFrame:
    assignment_key = next((key for key in GUIDE_ASSIGNMENT_KEYS if key in adata.obsm), None)
    if assignment_key is None:
        raise KeyError(
            f"No guide assignment matrix found. Tried {GUIDE_ASSIGNMENT_KEYS}; "
            f"available obsm keys={list(adata.obsm.keys())}"
        )

    assignment = adata.obsm[assignment_key]
    guide_names = guide_names_from_adata(adata, assignment)
    cell_ids = np.asarray(adata.obs_names.astype(str))

    if isinstance(assignment, pd.DataFrame):
        assignment_array = assignment.to_numpy()
    else:
        assignment_array = assignment

    umi_key = next(
        (
            key for key in GUIDE_UMI_MATRIX_KEYS
            if key in adata.obsm and getattr(adata.obsm[key], "shape", None) == assignment.shape
        ),
        None,
    )
    umi_matrix = adata.obsm[umi_key] if umi_key else None
    if isinstance(umi_matrix, pd.DataFrame):
        umi_matrix = umi_matrix.to_numpy()

    if sparse.issparse(assignment_array):
        coo = assignment_array.tocoo()
        keep = np.asarray(coo.data) != 0
        rows = np.asarray(coo.row)[keep]
        cols = np.asarray(coo.col)[keep]
        values = np.asarray(coo.data)[keep]
    else:
        dense = np.asarray(assignment_array)
        rows, cols = np.nonzero(dense)
        values = dense[rows, cols]

    if len(rows) == 0:
        raise ValueError(f"{day}: guide-assignment export would be empty")

    order = np.lexsort((cols, rows))
    rows = rows[order]
    cols = cols[order]
    values = np.asarray(values)[order]

    result = pd.DataFrame({
        "cell_id": cell_ids[rows],
        "guide_id": np.asarray(guide_names, dtype=object)[cols],
        "assignment_value": values,
    })

    if umi_matrix is not None:
        if sparse.issparse(umi_matrix):
            umi_values = np.asarray(umi_matrix.tocsr()[rows, cols]).reshape(-1)
        else:
            umi_values = np.asarray(umi_matrix)[rows, cols]
        result["guide_umi_count"] = umi_values

    if np.all(np.isfinite(pd.to_numeric(result["assignment_value"], errors="coerce"))):
        numeric = pd.to_numeric(result["assignment_value"])
        if np.allclose(numeric, np.rint(numeric)):
            result["assignment_value"] = np.rint(numeric).astype("int64")

    if result.duplicated(["cell_id", "guide_id"]).any():
        raise ValueError(f"{day}: duplicate cell_id + guide_id rows detected")
    if not set(result["cell_id"]).issubset(set(cell_ids)):
        raise ValueError(f"{day}: guide table contains cell IDs absent from the h5ad")
    return result


def extract_qc_report(adata: Any, day: str) -> pd.DataFrame:
    missing = [column for column in QC_COLUMNS if column not in adata.obs.columns]
    if missing:
        raise KeyError(
            f"{day}: required QC columns missing from adata.obs: {missing}. "
            f"Available columns={list(map(str, adata.obs.columns))}"
        )

    obs = adata.obs.loc[:, list(QC_COLUMNS)].copy()
    obs = dataframe_for_serialization(obs)
    obs.insert(0, "day", day.lower())
    obs.insert(0, "cell_id", adata.obs_names.astype(str))

    if obs["cell_id"].duplicated().any():
        raise ValueError(f"{day}: duplicate cell IDs in QC report")
    if len(obs) != adata.n_obs:
        raise ValueError(f"{day}: QC report row count does not equal h5ad n_obs")

    for column in (
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "guides_per_cell",
        "guide_umi_counts",
    ):
        numeric = pd.to_numeric(obs[column], errors="coerce")
        if numeric.isna().any():
            raise ValueError(f"{day}: {column} has missing or non-numeric values")
        obs[column] = numeric

    if obs["biological_sample"].isna().any():
        raise ValueError(f"{day}: biological_sample has missing values")
    return obs.reset_index(drop=True)


def validate_tsv_gz(path: Path) -> dict[str, Any]:
    report: dict[str, Any] = {
        "exists": path.is_file(),
        "gzip_ok": False,
        "row_count": 0,
        "column_count": 0,
        "duplicate_columns": "",
        "malformed_rows": 0,
    }
    if not path.is_file():
        report["error"] = "missing file"
        return report
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if not header:
                raise ValueError("empty gzip stream or missing header")
            report["gzip_ok"] = True
            report["column_count"] = len(header)
            dupes = sorted({column for column in header if header.count(column) > 1})
            report["duplicate_columns"] = ",".join(dupes)
            for row in reader:
                report["row_count"] += 1
                if len(row) != len(header):
                    report["malformed_rows"] += 1
        errors = []
        if report["column_count"] < 2:
            errors.append("not tabular")
        if report["row_count"] == 0:
            errors.append("no data rows")
        if report["duplicate_columns"]:
            errors.append("duplicate columns")
        if report["malformed_rows"]:
            errors.append("malformed rows")
        if errors:
            report["error"] = "; ".join(errors)
    except Exception as exc:
        report["error"] = f"{type(exc).__name__}: {exc}"
    return report

```

## Markdown cell 4

## 1. Resolve and validate day assignments

## Code cell 5

```python

# Read only obs metadata from the combined h5ad.
adata_backed = ad.read_h5ad(SOURCE_H5AD, backed="r")
try:
    if not adata_backed.obs_names.is_unique:
        raise ValueError("Combined h5ad obs_names are not unique")
    H5AD_DAY_COLUMN = find_day_column(adata_backed.obs, label="combined h5ad")
    h5ad_day = normalize_day_series(
        adata_backed.obs[H5AD_DAY_COLUMN],
        label=f"combined h5ad obs[{H5AD_DAY_COLUMN!r}]",
    )
    h5ad_day.index = pd.Index(adata_backed.obs_names.astype(str), name="cell_id")
    h5ad_obs_names = h5ad_day.index
    h5ad_n_vars = int(adata_backed.n_vars)
    h5ad_var_names_unique = bool(adata_backed.var_names.is_unique)
    h5ad_obsm_keys = list(map(str, adata_backed.obsm.keys()))
    h5ad_uns_keys = list(map(str, adata_backed.uns.keys()))
finally:
    close_backed_anndata(adata_backed)

print("h5ad day column:", H5AD_DAY_COLUMN)
print("h5ad cells by day:")
print(h5ad_day.value_counts().reindex(DAYS).to_string())
print("h5ad n_vars:", h5ad_n_vars)
print("h5ad var_names unique:", h5ad_var_names_unique)
print("h5ad obsm keys:", h5ad_obsm_keys)

# Load the combined h5mu. It is expected to be much smaller than the RNA h5ad
# because the cNMF modality contains cell x program scores.
mdata = read_h5mu(SOURCE_H5MU)
ANCHOR_MODALITY = find_modality(mdata, H5MU_ANCHOR_MODALITY)
anchor = mdata.mod[ANCHOR_MODALITY]

if not anchor.obs_names.is_unique:
    raise ValueError(f"h5mu modality {ANCHOR_MODALITY!r} obs_names are not unique")
if not anchor.var_names.is_unique:
    raise ValueError(f"h5mu modality {ANCHOR_MODALITY!r} var_names are not unique")

anchor_cells = pd.Index(anchor.obs_names.astype(str), name="cell_id")
missing_from_h5ad = anchor_cells[~anchor_cells.isin(h5ad_obs_names)]
extra_in_h5ad = h5ad_obs_names[~h5ad_obs_names.isin(anchor_cells)]

if len(missing_from_h5ad):
    raise ValueError(
        f"{len(missing_from_h5ad):,} h5mu cNMF cells are missing from the combined h5ad. "
        f"Examples: {missing_from_h5ad[:10].tolist()}"
    )
if REQUIRE_EXACT_H5AD_H5MU_CELL_SET and len(extra_in_h5ad):
    raise ValueError(
        f"Exact cell-set matching was requested, but the h5ad contains "
        f"{len(extra_in_h5ad):,} cells absent from cNMF. "
        f"Examples: {extra_in_h5ad[:10].tolist()}"
    )

# Day mapping aligned to the cNMF cell order.
anchor_day_from_h5ad = h5ad_day.reindex(anchor_cells)
if anchor_day_from_h5ad.isna().any():
    raise ValueError("Day assignment became missing after h5ad-to-h5mu cell alignment")

# Inspect any day columns already present in the h5mu. They are accepted only
# if they agree exactly with the h5ad mapping.
h5mu_day_sources = []
for location, obs in [("mdata.obs", mdata.obs)] + [
    (f"mod[{name}].obs", modality.obs) for name, modality in mdata.mod.items()
]:
    for candidate in DAY_CANDIDATES:
        if candidate not in obs.columns:
            continue
        mapped = obs[candidate].map(normalize_day)
        if mapped.isna().any():
            continue
        mapped.index = pd.Index(obs.index.astype(str))
        common = mapped.index.intersection(anchor_cells)
        if len(common) == 0:
            continue
        mismatches = int((mapped.reindex(common) != anchor_day_from_h5ad.reindex(common)).sum())
        h5mu_day_sources.append({
            "location": location,
            "column": candidate,
            "n_rows": int(len(mapped)),
            "n_common_anchor_cells": int(len(common)),
            "n_mismatches_vs_h5ad": mismatches,
            "values": ",".join(sorted(set(mapped))),
        })
        if mismatches:
            raise ValueError(
                f"Day mismatch between {location}[{candidate!r}] and "
                f"h5ad obs[{H5AD_DAY_COLUMN!r}]: {mismatches:,} cells"
            )

h5mu_day_sources_df = pd.DataFrame(h5mu_day_sources)
if h5mu_day_sources_df.empty:
    print("No usable day column was found in the h5mu; using the aligned h5ad obs day mapping.")
else:
    print("Usable h5mu day annotations (all agree with h5ad):")
    print(h5mu_day_sources_df.to_string(index=False))

print("h5mu modalities:", list(mdata.mod.keys()))
for name, modality in mdata.mod.items():
    print(f"  {name}: {modality.n_obs:,} observations x {modality.n_vars:,} variables")

print("Aligned cNMF cells by day:")
print(anchor_day_from_h5ad.value_counts().reindex(DAYS).to_string())

# Require global and modality-specific observations to be represented in the h5ad
# day map. This makes whole-MuData day slicing unambiguous across modalities.
global_cells = pd.Index(mdata.obs_names.astype(str), name="cell_id")
missing_global = global_cells[~global_cells.isin(h5ad_obs_names)]
if len(missing_global):
    raise ValueError(
        f"MuData global obs has {len(missing_global):,} cells absent from h5ad day metadata. "
        f"Examples: {missing_global[:10].tolist()}"
    )
global_day = h5ad_day.reindex(global_cells)
if global_day.isna().any():
    raise ValueError("MuData global day mapping contains missing values after alignment")

for name, modality in mdata.mod.items():
    modality_cells = pd.Index(modality.obs_names.astype(str))
    missing = modality_cells[~modality_cells.isin(h5ad_obs_names)]
    if len(missing):
        raise ValueError(
            f"Modality {name!r} has {len(missing):,} cells absent from h5ad day metadata. "
            f"Examples: {missing[:10].tolist()}"
        )
    modality_day = h5ad_day.reindex(modality_cells)
    if modality_day.isna().any():
        raise ValueError(f"Modality {name!r} has missing day values after alignment")

# Attach normalized day to the in-memory source object. The source file is never modified.
mdata.obs["day"] = global_day.str.lower().to_numpy()
for name, modality in mdata.mod.items():
    modality.obs["day"] = h5ad_day.reindex(
        pd.Index(modality.obs_names.astype(str))
    ).str.lower().to_numpy()

alignment_rows = []
for day in DAYS:
    h5ad_cells = pd.Index(h5ad_day.index[h5ad_day.eq(day)])
    cnmf_cells = pd.Index(anchor_day_from_h5ad.index[anchor_day_from_h5ad.eq(day)])
    alignment_rows.append({
        "day": day,
        "h5ad_cells": len(h5ad_cells),
        "h5mu_cnmf_cells": len(cnmf_cells),
        "same_cell_ids": set(h5ad_cells) == set(cnmf_cells),
        "analysis_set": ANALYSIS_SETS[day],
    })

alignment_df = pd.DataFrame(alignment_rows)
alignment_df.to_csv(WORK_DIR / "source_day_alignment.tsv", sep="\t", index=False)
h5mu_day_sources_df.to_csv(WORK_DIR / "h5mu_day_sources.tsv", sep="\t", index=False)
print("\nDay alignment summary:")
print(alignment_df.to_string(index=False))

if REQUIRE_EXACT_H5AD_H5MU_CELL_SET and not alignment_df["same_cell_ids"].all():
    raise ValueError("At least one day has different h5ad and cNMF cell IDs")

```

### Text output (stream)

```text
h5ad day column: day
h5ad cells by day:
day
D0    38609
D1    16925
D2    21980
D3    44102
h5ad n_vars: 17538
h5ad var_names unique: True
h5ad obsm keys: ['X_pca', 'X_umap', 'guide_assignment']

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

### Text output (stream)

```text
Usable h5mu day annotations (all agree with h5ad):
     location column  n_rows  n_common_anchor_cells  n_mismatches_vs_h5ad      values
 mod[rna].obs    day  121616                 121616                     0 D0,D1,D2,D3
mod[cNMF].obs    day  121616                 121616                     0 D0,D1,D2,D3
h5mu modalities: ['rna', 'cNMF']
  rna: 121,616 observations x 17,538 variables
  cNMF: 121,616 observations x 50 variables
Aligned cNMF cells by day:
day
D0    38609
D1    16925
D2    21980
D3    44102

Day alignment summary:
day  h5ad_cells  h5mu_cnmf_cells  same_cell_ids   analysis_set
 D0       38609            38609           True IGVFDS0463XTUH
 D1       16925            16925           True IGVFDS0991DBQZ
 D2       21980            21980           True IGVFDS4402ZGWJ
 D3       44102            44102           True IGVFDS2307REVJ

```

## Markdown cell 6

## 2. Build the exact 16-file output plan

## Code cell 7

```python

# Raw SequenceFile accessions are used only as derived_from links for the four
# h5ad MatrixFiles. These records are already uploaded and are not recreated.
seq = pd.read_csv(SEQUENCE_ACCESSION_CSV)
required_seq_columns = {"accession", "day", "modality", "sample_key", "read"}
missing_seq_columns = required_seq_columns - set(seq.columns)
if missing_seq_columns:
    raise KeyError(f"Sequence accession CSV is missing columns: {sorted(missing_seq_columns)}")

for status_column in ("post_success", "aws_upload_success"):
    if status_column in seq.columns:
        seq = seq[seq[status_column].astype(str).str.lower().isin({"true", "1"})].copy()

seq["day"] = seq["day"].map(normalize_day)
if seq["day"].isna().any():
    raise ValueError("Sequence accession CSV contains unmapped day values")
if seq["accession"].duplicated().any():
    raise ValueError("Sequence accession CSV contains duplicate accessions")

raw_accessions_by_day = {}
for day in DAYS:
    day_seq = seq[seq["day"].eq(day)].copy()
    modalities = set(day_seq["modality"].astype(str).str.lower())
    if not any(value in modalities for value in {"gex", "rna", "scrna"}):
        raise ValueError(f"{day}: no GEX/RNA SequenceFile accessions found")
    if not any(value in modalities for value in {"grna", "guide", "crispr"}):
        raise ValueError(f"{day}: no guide/gRNA SequenceFile accessions found")
    day_seq = day_seq.sort_values(["modality", "sample_key", "read", "accession"])
    raw_accessions_by_day[day] = day_seq["accession"].astype(str).tolist()
    print(f"{day}: {len(raw_accessions_by_day[day])} raw SequenceFile accessions")

PREPROCESS_STEP = (
    "jesse-engreitz:analysis_step_version_WTC11-EC_300gene_"
    "CC-Perturbseq_preprocessing_v0.1.0"
)
PERTURBNMF_STEP = (
    "jesse-engreitz:analysis_step_version_WTC11-EC_300gene_"
    "PerturbNMF_gene_programs_v0.1.0"
)

records = []
for day in DAYS:
    h5ad_alias = f"jesse-engreitz:matrix_file_WTC11-EC_300gene_CC-Perturbseq_{day}_h5ad"
    h5mu_alias = f"jesse-engreitz:matrix_file_WTC11-EC_300gene_PerturbNMF_{day}_h5mu"

    h5ad_path = OUTPUT_ROOT / day / "cc_perturbseq" / f"CCPerturbSeq_{day}_annotated_cell_gene_guide.h5ad"
    h5mu_path = OUTPUT_ROOT / day / "perturbnmf" / f"PerturbNMF_{day}_cell_gene_program.h5mu"
    guide_path = OUTPUT_ROOT / day / "cc_perturbseq" / "guide_quantifications.tsv.gz"
    qc_path = OUTPUT_ROOT / day / "cc_perturbseq" / "per_cell_quality_report.tsv.gz"

    records.extend([
        {
            "day": day,
            "profile": "matrix_file",
            "role": "annotated cell by gene and guide matrix",
            "alias": h5ad_alias,
            "file_set": ANALYSIS_SETS[day],
            "file_format": "h5ad",
            "content_type": "annotated cell by gene and guide matrix",
            "submitted_file_name": str(h5ad_path),
            "analysis_step_version": PREPROCESS_STEP,
            "derived_from": ",".join(raw_accessions_by_day[day]),
            "reference_files": REFERENCE_FILE,
            "filtered": "TRUE",
            "controlled_access": "",
            "derived_manually": "FALSE",
            "description": (
                f"Annotated cell-by-gene-and-guide AnnData matrix for WTC11-EC "
                f"300-gene CC-Perturb-seq at {day}."
            ),
        },
        {
            "day": day,
            "profile": "matrix_file",
            "role": "cell by gene program matrix",
            "alias": h5mu_alias,
            "file_set": ANALYSIS_SETS[day],
            "file_format": "h5mu",
            "content_type": "cell by gene program matrix",
            "submitted_file_name": str(h5mu_path),
            "analysis_step_version": PERTURBNMF_STEP,
            "derived_from": h5ad_alias,
            "reference_files": REFERENCE_FILE,
            "filtered": "TRUE",
            "controlled_access": "",
            "derived_manually": "FALSE",
            "description": (
                f"PerturbNMF cell-by-gene-program MuData output for WTC11-EC "
                f"300-gene CC-Perturb-seq at {day}."
            ),
        },
        {
            "day": day,
            "profile": "tabular_file",
            "role": "guide quantifications",
            "alias": f"jesse-engreitz:tabular_file_WTC11-EC_300gene_{day}_guide_quantifications",
            "file_set": ANALYSIS_SETS[day],
            "file_format": "tsv",
            "content_type": "guide quantifications",
            "submitted_file_name": str(guide_path),
            "analysis_step_version": PREPROCESS_STEP,
            "derived_from": h5ad_alias,
            "reference_files": REFERENCE_FILE,
            "filtered": "TRUE",
            "controlled_access": "FALSE",
            "derived_manually": "FALSE",
            "description": f"Per-cell guide assignment quantifications for the {day} analysis.",
        },
        {
            "day": day,
            "profile": "tabular_file",
            "role": "per-cell quality report",
            "alias": f"jesse-engreitz:tabular_file_WTC11-EC_300gene_{day}_per_cell_quality_report",
            "file_set": ANALYSIS_SETS[day],
            "file_format": "tsv",
            "content_type": "per-cell quality report",
            "submitted_file_name": str(qc_path),
            "analysis_step_version": PREPROCESS_STEP,
            "derived_from": h5ad_alias,
            "reference_files": REFERENCE_FILE,
            "filtered": "FALSE",
            "controlled_access": "FALSE",
            "derived_manually": "FALSE",
            "description": f"Per-cell quality-control metrics for the {day} analysis.",
        },
    ])

plan = pd.DataFrame(records)
if len(plan) != 16:
    raise AssertionError(f"Expected 16 planned files, found {len(plan)}")
if plan["alias"].duplicated().any() or plan["submitted_file_name"].duplicated().any():
    raise ValueError("Duplicate aliases or output paths in the plan")
if not all(plan.groupby("day")["file_set"].nunique().eq(1)):
    raise ValueError("A timepoint maps to more than one AnalysisSet")

plan.to_csv(WORK_DIR / "planned_16_files.tsv", sep="\t", index=False)
print(plan[["day", "profile", "role", "file_set", "submitted_file_name"]].to_string(index=False))

```

### Text output (stream)

```text
D0: 24 raw SequenceFile accessions
D1: 24 raw SequenceFile accessions
D2: 24 raw SequenceFile accessions
D3: 24 raw SequenceFile accessions
day      profile                                    role       file_set                                                                                                                                    submitted_file_name
 D0  matrix_file annotated cell by gene and guide matrix IGVFDS0463XTUH /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D0/cc_perturbseq/CCPerturbSeq_D0_annotated_cell_gene_guide.h5ad
 D0  matrix_file             cell by gene program matrix IGVFDS0463XTUH              /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D0/perturbnmf/PerturbNMF_D0_cell_gene_program.h5mu
 D0 tabular_file                   guide quantifications IGVFDS0463XTUH                   /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D0/cc_perturbseq/guide_quantifications.tsv.gz
 D0 tabular_file                 per-cell quality report IGVFDS0463XTUH                 /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D0/cc_perturbseq/per_cell_quality_report.tsv.gz
 D1  matrix_file annotated cell by gene and guide matrix IGVFDS0991DBQZ /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D1/cc_perturbseq/CCPerturbSeq_D1_annotated_cell_gene_guide.h5ad
 D1  matrix_file             cell by gene program matrix IGVFDS0991DBQZ              /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D1/perturbnmf/PerturbNMF_D1_cell_gene_program.h5mu
 D1 tabular_file                   guide quantifications IGVFDS0991DBQZ                   /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D1/cc_perturbseq/guide_quantifications.tsv.gz
 D1 tabular_file                 per-cell quality report IGVFDS0991DBQZ                 /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D1/cc_perturbseq/per_cell_quality_report.tsv.gz
 D2  matrix_file annotated cell by gene and guide matrix IGVFDS4402ZGWJ /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D2/cc_perturbseq/CCPerturbSeq_D2_annotated_cell_gene_guide.h5ad
 D2  matrix_file             cell by gene program matrix IGVFDS4402ZGWJ              /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D2/perturbnmf/PerturbNMF_D2_cell_gene_program.h5mu
 D2 tabular_file                   guide quantifications IGVFDS4402ZGWJ                   /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D2/cc_perturbseq/guide_quantifications.tsv.gz
 D2 tabular_file                 per-cell quality report IGVFDS4402ZGWJ                 /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D2/cc_perturbseq/per_cell_quality_report.tsv.gz
 D3  matrix_file annotated cell by gene and guide matrix IGVFDS2307REVJ /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D3/cc_perturbseq/CCPerturbSeq_D3_annotated_cell_gene_guide.h5ad
 D3  matrix_file             cell by gene program matrix IGVFDS2307REVJ              /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D3/perturbnmf/PerturbNMF_D3_cell_gene_program.h5mu
 D3 tabular_file                   guide quantifications IGVFDS2307REVJ                   /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D3/cc_perturbseq/guide_quantifications.tsv.gz
 D3 tabular_file                 per-cell quality report IGVFDS2307REVJ                 /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D3/cc_perturbseq/per_cell_quality_report.tsv.gz

```

## Markdown cell 8

## 3. Write the four final day-specific h5ad files

## Code cell 9

```python

def validate_existing_h5ad(path: Path, day: str) -> None:
    validate_hdf5_container(path)
    check = ad.read_h5ad(path, backed="r")
    try:
        if check.n_obs == 0 or check.n_vars == 0:
            raise ValueError(f"{path}: empty matrix")
        if not check.obs_names.is_unique:
            raise ValueError(f"{path}: duplicate obs_names")
        if not check.var_names.is_unique:
            raise ValueError(f"{path}: duplicate var_names")
        day_column = find_day_column(check.obs, label=str(path))
        observed = check.obs[day_column].map(normalize_day)
        if observed.isna().any() or set(observed.unique()) != {day}:
            raise ValueError(f"{path}: expected only {day}, observed {sorted(set(observed.dropna()))}")
        assignment_key = next((key for key in GUIDE_ASSIGNMENT_KEYS if key in check.obsm), None)
        if assignment_key is None:
            raise ValueError(f"{path}: no guide assignment matrix")
        matrix = check.obsm[assignment_key]
        if matrix.shape[0] != check.n_obs:
            raise ValueError(f"{path}: guide assignment rows do not equal n_obs")
        guide_names_from_adata(check, matrix)
    finally:
        close_backed_anndata(check)


h5ad_write_rows = []
source = ad.read_h5ad(SOURCE_H5AD, backed="r")
try:
    source_day = normalize_day_series(
        source.obs[H5AD_DAY_COLUMN],
        label=f"source h5ad obs[{H5AD_DAY_COLUMN!r}]",
    )
    source_day.index = pd.Index(source.obs_names.astype(str))

    for day in DAYS:
        output_path = Path(
            plan.loc[
                (plan["day"].eq(day)) & (plan["file_format"].eq("h5ad")),
                "submitted_file_name",
            ].iloc[0]
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        expected_cells = pd.Index(source_day.index[source_day.eq(day)])

        action = "reused_existing"
        if output_path.exists() and not OVERWRITE_EXISTING:
            validate_existing_h5ad(output_path, day)
        else:
            if not EXECUTE_WRITES:
                raise RuntimeError("EXECUTE_WRITES=False; refusing to create final h5ad files")
            action = "written"
            view = source[expected_cells]
            subset = view.to_memory() if hasattr(view, "to_memory") else view.copy()
            subset.obs[H5AD_DAY_COLUMN] = day.lower()
            subset.obs["igvf_submission_day"] = day
            subset.uns["igvf_submission_conversion"] = {
                "source_file": str(SOURCE_H5AD),
                "day_source": f"obs[{H5AD_DAY_COLUMN}]",
                "day": day,
                "generated_utc": now_utc(),
                "method": "cell subset by exact obs_names and normalized day",
            }
            temp_path = output_path.parent / f".{output_path.stem}.tmp.h5ad"
            if temp_path.exists():
                temp_path.unlink()
            subset.write_h5ad(temp_path, compression="gzip")
            validate_existing_h5ad(temp_path, day)
            atomic_replace(temp_path, output_path)
            del subset
            gc.collect()

        validate_existing_h5ad(output_path, day)
        check = ad.read_h5ad(output_path, backed="r")
        try:
            actual_cells = pd.Index(check.obs_names.astype(str))
            if set(actual_cells) != set(expected_cells):
                raise ValueError(f"{day}: final h5ad cell IDs differ from the source day subset")
            h5ad_write_rows.append({
                "day": day,
                "action": action,
                "output_path": str(output_path),
                "n_obs": int(check.n_obs),
                "n_vars": int(check.n_vars),
            })
        finally:
            close_backed_anndata(check)
finally:
    close_backed_anndata(source)

h5ad_write_report = pd.DataFrame(h5ad_write_rows)
h5ad_write_report.to_csv(WORK_DIR / "h5ad_write_report.tsv", sep="\t", index=False)
print(h5ad_write_report.to_string(index=False))

```

### Text output (stream)

```text
day  action                                                                                                                                            output_path  n_obs  n_vars
 D0 written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D0/cc_perturbseq/CCPerturbSeq_D0_annotated_cell_gene_guide.h5ad  38609   17538
 D1 written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D1/cc_perturbseq/CCPerturbSeq_D1_annotated_cell_gene_guide.h5ad  16925   17538
 D2 written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D2/cc_perturbseq/CCPerturbSeq_D2_annotated_cell_gene_guide.h5ad  21980   17538
 D3 written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D3/cc_perturbseq/CCPerturbSeq_D3_annotated_cell_gene_guide.h5ad  44102   17538

```

## Markdown cell 10

## 4. Write the four final day-specific h5mu files

## Code cell 11

```python

# Capture source modality dimensions/variables before splitting.
source_modality_summary = {
    name: {
        "n_obs": int(modality.n_obs),
        "n_vars": int(modality.n_vars),
        "var_names": list(map(str, modality.var_names)),
    }
    for name, modality in mdata.mod.items()
}


def validate_existing_h5mu(path: Path, day: str, expected_cells: pd.Index) -> dict[str, Any]:
    validate_hdf5_container(path)
    check = read_h5mu(path)
    try:
        anchor_name = find_modality(check, H5MU_ANCHOR_MODALITY)
        anchor_check = check.mod[anchor_name]
        actual_cells = pd.Index(anchor_check.obs_names.astype(str))
        if not actual_cells.is_unique:
            raise ValueError(f"{path}: duplicate cNMF obs_names")
        if set(actual_cells) != set(expected_cells):
            missing = expected_cells[~expected_cells.isin(actual_cells)]
            extra = actual_cells[~actual_cells.isin(expected_cells)]
            raise ValueError(
                f"{path}: cNMF cell IDs do not equal the expected {day} subset; "
                f"missing={len(missing)}, extra={len(extra)}"
            )
        if "day" not in anchor_check.obs.columns:
            raise ValueError(f"{path}: cNMF obs does not contain day")
        observed = anchor_check.obs["day"].map(normalize_day)
        if observed.isna().any() or set(observed.unique()) != {day}:
            raise ValueError(f"{path}: expected only {day} in cNMF obs['day']")
        if anchor_check.n_vars == 0:
            raise ValueError(f"{path}: cNMF modality has zero programs")
        for name, modality in check.mod.items():
            if name not in source_modality_summary:
                raise ValueError(f"{path}: unexpected modality {name!r}")
            source_info = source_modality_summary[name]
            if int(modality.n_vars) != source_info["n_vars"]:
                raise ValueError(f"{path}: modality {name!r} n_vars changed")
            if list(map(str, modality.var_names)) != source_info["var_names"]:
                raise ValueError(f"{path}: modality {name!r} var_names changed")
            if "day" not in modality.obs.columns:
                raise ValueError(f"{path}: modality {name!r} lacks obs['day']")
            modality_days = modality.obs["day"].map(normalize_day)
            if modality_days.isna().any() or set(modality_days.unique()) != {day}:
                raise ValueError(f"{path}: modality {name!r} has cross-day observations")
        return {
            "n_obs": int(anchor_check.n_obs),
            "n_programs": int(anchor_check.n_vars),
            "modalities": ",".join(map(str, check.mod.keys())),
        }
    finally:
        close_backed_mudata(check)


h5mu_write_rows = []
for day in DAYS:
    output_path = Path(
        plan.loc[
            (plan["day"].eq(day)) & (plan["file_format"].eq("h5mu")),
            "submitted_file_name",
        ].iloc[0]
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    expected_cells = pd.Index(anchor_day_from_h5ad.index[anchor_day_from_h5ad.eq(day)])

    action = "reused_existing"
    if output_path.exists() and not OVERWRITE_EXISTING:
        details = validate_existing_h5mu(output_path, day, expected_cells)
    else:
        if not EXECUTE_WRITES:
            raise RuntimeError("EXECUTE_WRITES=False; refusing to create final h5mu files")
        action = "written"

        # Whole-MuData slicing subsets observations across modalities. The exact
        # day cell IDs were resolved above from the h5ad and checked against cNMF.
        subset = mdata[expected_cells, :].copy()
        subset.obs["day"] = day.lower()
        subset.obs["igvf_submission_day"] = day
        for name, modality in subset.mod.items():
            modality.obs["day"] = day.lower()
            modality.obs["igvf_submission_day"] = day
        subset.uns["igvf_submission_conversion"] = {
            "source_file": str(SOURCE_H5MU),
            "day_source": f"{SOURCE_H5AD}:obs[{H5AD_DAY_COLUMN}] aligned by obs_names",
            "day": day,
            "generated_utc": now_utc(),
            "method": "MuData cell subset by exact observation IDs",
        }

        temp_path = output_path.parent / f".{output_path.stem}.tmp.h5mu"
        if temp_path.exists():
            temp_path.unlink()
        write_h5mu(subset, temp_path)
        validate_existing_h5mu(temp_path, day, expected_cells)
        atomic_replace(temp_path, output_path)
        del subset
        gc.collect()
        details = validate_existing_h5mu(output_path, day, expected_cells)

    h5mu_write_rows.append({
        "day": day,
        "action": action,
        "output_path": str(output_path),
        **details,
    })

h5mu_write_report = pd.DataFrame(h5mu_write_rows)
h5mu_write_report.to_csv(WORK_DIR / "h5mu_write_report.tsv", sep="\t", index=False)
print(h5mu_write_report.to_string(index=False))

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

### Text output (stream)

```text
day  action                                                                                                                               output_path  n_obs  n_programs modalities
 D0 written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D0/perturbnmf/PerturbNMF_D0_cell_gene_program.h5mu  38609          50   rna,cNMF
 D1 written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D1/perturbnmf/PerturbNMF_D1_cell_gene_program.h5mu  16925          50   rna,cNMF
 D2 written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D2/perturbnmf/PerturbNMF_D2_cell_gene_program.h5mu  21980          50   rna,cNMF
 D3 written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D3/perturbnmf/PerturbNMF_D3_cell_gene_program.h5mu  44102          50   rna,cNMF

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

## Markdown cell 12

## 5. Create the final guide-quantification and per-cell QC tables

## Code cell 13

```python

tabular_write_rows = []
column_dictionary_rows = []

for day in DAYS:
    h5ad_path = Path(
        plan.loc[
            (plan["day"].eq(day)) & (plan["file_format"].eq("h5ad")),
            "submitted_file_name",
        ].iloc[0]
    )
    guide_path = Path(
        plan.loc[
            (plan["day"].eq(day)) & (plan["role"].eq("guide quantifications")),
            "submitted_file_name",
        ].iloc[0]
    )
    qc_path = Path(
        plan.loc[
            (plan["day"].eq(day)) & (plan["role"].eq("per-cell quality report")),
            "submitted_file_name",
        ].iloc[0]
    )

    day_adata = ad.read_h5ad(h5ad_path)
    try:
        guide_df = extract_guide_quantifications(day_adata, day)
        qc_df = extract_qc_report(day_adata, day)

        if guide_path.exists() and not OVERWRITE_EXISTING:
            old = pd.read_csv(guide_path, sep="\t", compression="gzip")
            if not old.equals(guide_df):
                raise FileExistsError(
                    f"{guide_path} exists but differs from the newly derived guide table. "
                    "Set OVERWRITE_EXISTING=True only after reviewing the difference."
                )
            guide_action = "reused_existing"
        else:
            if not EXECUTE_WRITES:
                raise RuntimeError("EXECUTE_WRITES=False; refusing to create guide tables")
            write_deterministic_tsv_gz(guide_df, guide_path)
            guide_action = "written"

        if qc_path.exists() and not OVERWRITE_EXISTING:
            old = pd.read_csv(qc_path, sep="\t", compression="gzip")
            if not old.equals(qc_df):
                raise FileExistsError(
                    f"{qc_path} exists but differs from the newly derived QC table. "
                    "Set OVERWRITE_EXISTING=True only after reviewing the difference."
                )
            qc_action = "reused_existing"
        else:
            if not EXECUTE_WRITES:
                raise RuntimeError("EXECUTE_WRITES=False; refusing to create QC tables")
            write_deterministic_tsv_gz(qc_df, qc_path)
            qc_action = "written"

        tabular_write_rows.extend([
            {
                "day": day,
                "role": "guide quantifications",
                "action": guide_action,
                "output_path": str(guide_path),
                "row_count": len(guide_df),
                "column_count": guide_df.shape[1],
            },
            {
                "day": day,
                "role": "per-cell quality report",
                "action": qc_action,
                "output_path": str(qc_path),
                "row_count": len(qc_df),
                "column_count": qc_df.shape[1],
            },
        ])

        # Build a local column dictionary. No Portal Document accession is invented.
        guide_descriptions = {
            "cell_id": ("string", "yes", "identifier", "NA", "Cell identifier from the day-specific h5ad obs_names.", ""),
            "guide_id": ("string", "yes", "identifier", "NA", "Guide identifier from the h5ad guide-assignment matrix.", ""),
            "assignment_value": ("integer or float", "yes", "as stored", "NA", "Non-zero guide assignment value stored for the cell-guide pair.", ""),
            "guide_umi_count": ("integer or float", "no", "guide UMI count", "NA", "Guide UMI count for the cell-guide pair when a matching matrix is present.", ""),
        }
        for column in guide_df.columns:
            dtype, required, units, missing, description, allowed = guide_descriptions[column]
            column_dictionary_rows.append({
                "table_type": "guide quantifications",
                "column_name": column,
                "data_type": dtype,
                "required": required,
                "units_or_scale": units,
                "allowed_values": allowed,
                "missing_value": missing,
                "description": description,
                "source_analysis_step": "Step 1: CC-Perturb-seq preprocessing and guide assignment",
            })

        qc_descriptions = {
            "cell_id": ("string", "yes", "identifier", "NA", "Cell identifier from h5ad obs_names.", ""),
            "day": ("categorical", "yes", "differentiation day", "NA", "Differentiation timepoint for the cell.", "d0,d1,d2,d3"),
            "biological_sample": ("categorical", "yes", "sample identifier", "NA", "Biological sample or replicate label from h5ad obs.", ""),
            "total_counts": ("numeric", "yes", "RNA UMI counts per cell", "NA", "Total RNA counts for the cell.", "non-negative"),
            "n_genes_by_counts": ("integer", "yes", "genes per cell", "NA", "Number of detected genes for the cell.", "non-negative"),
            "pct_counts_mt": ("float", "yes", "percent", "NA", "Percentage of counts assigned to mitochondrial genes.", "0 to 100"),
            "guides_per_cell": ("integer", "yes", "assigned guides per cell", "NA", "Number of assigned guides for the cell.", "non-negative"),
            "guide_umi_counts": ("numeric", "yes", "guide UMI counts per cell", "NA", "Total guide UMI count for the cell.", "non-negative"),
        }
        for column in qc_df.columns:
            dtype, required, units, missing, description, allowed = qc_descriptions[column]
            column_dictionary_rows.append({
                "table_type": "per-cell quality report",
                "column_name": column,
                "data_type": dtype,
                "required": required,
                "units_or_scale": units,
                "allowed_values": allowed,
                "missing_value": missing,
                "description": description,
                "source_analysis_step": "Step 1: CC-Perturb-seq preprocessing and guide assignment",
            })
    finally:
        del day_adata
        gc.collect()

# The dictionary is shared across days; remove repeated definitions.
column_dictionary = pd.DataFrame(column_dictionary_rows).drop_duplicates().reset_index(drop=True)
column_dictionary_path = DOCUMENTATION_DIR / "column_dictionary_step1.tsv"
column_dictionary.to_csv(column_dictionary_path, sep="\t", index=False, lineterminator="\n")

tabular_write_report = pd.DataFrame(tabular_write_rows)
tabular_write_report.to_csv(WORK_DIR / "tabular_write_report.tsv", sep="\t", index=False)
print(tabular_write_report.to_string(index=False))
print("Column dictionary:", column_dictionary_path)

```

### Text output (stream)

```text
day                    role  action                                                                                                                            output_path  row_count  column_count
 D0   guide quantifications written   /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D0/cc_perturbseq/guide_quantifications.tsv.gz     206613             3
 D0 per-cell quality report written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D0/cc_perturbseq/per_cell_quality_report.tsv.gz      38609             8
 D1   guide quantifications written   /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D1/cc_perturbseq/guide_quantifications.tsv.gz     105523             3
 D1 per-cell quality report written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D1/cc_perturbseq/per_cell_quality_report.tsv.gz      16925             8
 D2   guide quantifications written   /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D2/cc_perturbseq/guide_quantifications.tsv.gz     135227             3
 D2 per-cell quality report written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D2/cc_perturbseq/per_cell_quality_report.tsv.gz      21980             8
 D3   guide quantifications written   /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D3/cc_perturbseq/guide_quantifications.tsv.gz     199129             3
 D3 per-cell quality report written /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/D3/cc_perturbseq/per_cell_quality_report.tsv.gz      44102             8
Column dictionary: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/documentation/column_dictionary_step1.tsv

```

## Markdown cell 14

## 6. Validate all 16 final files and compute MD5 checksums

## Code cell 15

```python

validation_rows = []

# Source cNMF variables for output comparison.
source_anchor_var_names = list(map(str, mdata.mod[ANCHOR_MODALITY].var_names))

for row in plan.to_dict("records"):
    path = Path(row["submitted_file_name"])
    status = "PASS"
    issues = []
    details: dict[str, Any] = {}

    try:
        if not path.is_file():
            raise FileNotFoundError(path)
        if path.stat().st_size == 0:
            raise ValueError("zero-byte file")

        if row["file_format"] == "h5ad":
            validate_hdf5_container(path)
            obj = ad.read_h5ad(path, backed="r")
            try:
                if not obj.obs_names.is_unique:
                    issues.append("duplicate obs_names")
                if not obj.var_names.is_unique:
                    issues.append("duplicate var_names")
                day_column = find_day_column(obj.obs, label=str(path))
                observed = obj.obs[day_column].map(normalize_day)
                if observed.isna().any() or set(observed.unique()) != {row["day"]}:
                    issues.append("wrong or mixed timepoint")
                assignment_key = next((key for key in GUIDE_ASSIGNMENT_KEYS if key in obj.obsm), None)
                if assignment_key is None:
                    issues.append("missing guide-assignment matrix")
                else:
                    matrix = obj.obsm[assignment_key]
                    if matrix.shape[0] != obj.n_obs:
                        issues.append("guide-assignment row count mismatch")
                    guide_names_from_adata(obj, matrix)
                details.update({
                    "n_obs": int(obj.n_obs),
                    "n_vars": int(obj.n_vars),
                    "modalities": "",
                    "row_count": "",
                    "column_count": "",
                    "obs_keys": ",".join(map(str, obj.obs.columns)),
                    "obsm_keys": ",".join(map(str, obj.obsm.keys())),
                    "layers": ",".join(map(str, obj.layers.keys())),
                })
            finally:
                close_backed_anndata(obj)

        elif row["file_format"] == "h5mu":
            validate_hdf5_container(path)
            obj = read_h5mu(path)
            try:
                anchor_name = find_modality(obj, H5MU_ANCHOR_MODALITY)
                anchor_obj = obj.mod[anchor_name]
                if not anchor_obj.obs_names.is_unique:
                    issues.append("duplicate cNMF obs_names")
                observed = anchor_obj.obs["day"].map(normalize_day) if "day" in anchor_obj.obs else pd.Series(dtype=object)
                if observed.empty or observed.isna().any() or set(observed.unique()) != {row["day"]}:
                    issues.append("wrong, mixed, or missing cNMF day")
                if list(map(str, anchor_obj.var_names)) != source_anchor_var_names:
                    issues.append("cNMF program identifiers differ from source")
                expected_h5ad_path = Path(
                    plan.loc[
                        (plan["day"].eq(row["day"])) & (plan["file_format"].eq("h5ad")),
                        "submitted_file_name",
                    ].iloc[0]
                )
                upstream = ad.read_h5ad(expected_h5ad_path, backed="r")
                try:
                    if set(anchor_obj.obs_names.astype(str)) != set(upstream.obs_names.astype(str)):
                        issues.append("h5mu cell IDs do not match day-specific h5ad")
                finally:
                    close_backed_anndata(upstream)
                details.update({
                    "n_obs": int(anchor_obj.n_obs),
                    "n_vars": int(anchor_obj.n_vars),
                    "modalities": ",".join(map(str, obj.mod.keys())),
                    "row_count": "",
                    "column_count": "",
                    "obs_keys": ",".join(map(str, anchor_obj.obs.columns)),
                    "obsm_keys": ",".join(map(str, anchor_obj.obsm.keys())),
                    "layers": ",".join(map(str, anchor_obj.layers.keys())),
                })
            finally:
                close_backed_mudata(obj)

        elif row["file_format"] == "tsv":
            report = validate_tsv_gz(path)
            if report.get("error"):
                issues.append(report["error"])
            table = pd.read_csv(path, sep="\t", compression="gzip")
            if row["role"] == "guide quantifications":
                required = {"cell_id", "guide_id", "assignment_value"}
                if not required.issubset(table.columns):
                    issues.append(f"missing guide columns: {sorted(required - set(table.columns))}")
                if table.duplicated(["cell_id", "guide_id"]).any():
                    issues.append("duplicate cell_id + guide_id rows")
                if (pd.to_numeric(table["assignment_value"], errors="coerce") == 0).any():
                    issues.append("zero assignment values present in non-zero export")
            elif row["role"] == "per-cell quality report":
                required = {"cell_id", "day", *QC_COLUMNS}
                if not required.issubset(table.columns):
                    issues.append(f"missing QC columns: {sorted(required - set(table.columns))}")
                if table["cell_id"].duplicated().any():
                    issues.append("duplicate QC cell_id rows")
                observed = table["day"].map(normalize_day)
                if observed.isna().any() or set(observed.unique()) != {row["day"]}:
                    issues.append("wrong or mixed QC day")
                expected_h5ad_path = Path(
                    plan.loc[
                        (plan["day"].eq(row["day"])) & (plan["file_format"].eq("h5ad")),
                        "submitted_file_name",
                    ].iloc[0]
                )
                upstream = ad.read_h5ad(expected_h5ad_path, backed="r")
                try:
                    if len(table) != upstream.n_obs:
                        issues.append("QC row count does not match h5ad n_obs")
                    if set(table["cell_id"].astype(str)) != set(upstream.obs_names.astype(str)):
                        issues.append("QC cell IDs do not match h5ad")
                finally:
                    close_backed_anndata(upstream)
            details.update({
                "n_obs": "",
                "n_vars": "",
                "modalities": "",
                "row_count": int(report["row_count"]),
                "column_count": int(report["column_count"]),
                "obs_keys": "",
                "obsm_keys": "",
                "layers": "",
            })
        else:
            issues.append(f"unsupported file_format {row['file_format']}")
    except Exception as exc:
        issues.append(f"{type(exc).__name__}: {exc}")

    if issues:
        status = "ERROR"

    validation_rows.append({
        "day": row["day"],
        "profile": row["profile"],
        "role": row["role"],
        "alias": row["alias"],
        "file_set": row["file_set"],
        "submitted_file_name": str(path),
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.is_file() else "",
        "md5sum": md5sum(path) if path.is_file() else "",
        "validation_status": status,
        "issues": "; ".join(issues),
        **details,
    })

validation = pd.DataFrame(validation_rows)
validation_path = WORK_DIR / "file_validation_16.tsv"
validation.to_csv(validation_path, sep="\t", index=False, lineterminator="\n")

print(validation[[
    "day", "role", "size_bytes", "md5sum", "validation_status", "issues"
]].to_string(index=False))

failures = validation[~validation["validation_status"].eq("PASS")]
if not failures.empty:
    raise RuntimeError(
        f"{len(failures)} of 16 files failed validation. Review {validation_path}."
    )
print("All 16 final files passed validation.")

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

### Text output (stream)

```text
day                                    role  size_bytes                           md5sum validation_status issues
 D0 annotated cell by gene and guide matrix   312724458 0d3fe724ad9eca6cb9e811b71e67b2d7              PASS       
 D0             cell by gene program matrix   593116718 eaa9b998f4a15fee771f2c0236c63c5a              PASS       
 D0                   guide quantifications     1701892 e3192ea1c4d02ce875782cacbe1ba99b              PASS       
 D0                 per-cell quality report      853217 ea3ca228602d7f433855184555d855b3              PASS       
 D1 annotated cell by gene and guide matrix   143568531 8ffe7570d640a4f72b2de1d6d68e0a01              PASS       
 D1             cell by gene program matrix   267183838 cdc55d85a86b62355bbba8dd97669d11              PASS       
 D1                   guide quantifications      854776 668a5bc7d6c42d294a8aedae485bf20d              PASS       
 D1                 per-cell quality report      370934 819db5d7e0e81959708efaf4c550fc1e              PASS       
 D2 annotated cell by gene and guide matrix   121143019 f2fe9f4af11a3eda18617ce4b7ec36f6              PASS       
 D2             cell by gene program matrix   235830412 ea85d629fd57646848ba6d31a20cc3b3              PASS       
 D2                   guide quantifications     1096385 591144d0e9db49a972662e727ba9ac64              PASS       
 D2                 per-cell quality report      478069 3f1fc50882a9b6cb469ece4914933f94              PASS       
 D3 annotated cell by gene and guide matrix   229733444 e8c726fb0e634074384bc8b7cf5a4464              PASS       
 D3             cell by gene program matrix   436510560 6920d6bdee235a327d955ff038c6ea82              PASS       
 D3                   guide quantifications     1656396 e0e4a84989f2179d80a91bbd5c45a8c2              PASS       
 D3                 per-cell quality report      956777 ad55846b58d783c7394e6338a104ac2f              PASS       
All 16 final files passed validation.

```

## Markdown cell 16

## 7. Generate upload manifests and IGVF metadata TSVs

## Code cell 17

```python

# Merge final checksums into the plan.
final = plan.merge(
    validation[[
        "alias", "size_bytes", "md5sum", "validation_status", "issues"
    ]],
    on="alias",
    how="left",
    validate="one_to_one",
)
if not final["validation_status"].eq("PASS").all():
    raise RuntimeError("Metadata generation is blocked because a file failed validation")

final["award"] = AWARD
final["lab"] = LAB
final["submitter_comment"] = (
    "Generated as a day-specific, validated submission file by "
    "01_make_igvf_upload_ready_files.ipynb; MD5 calculated from the exact final file."
)

# Master 16-file manifest.
manifest_columns = [
    "day", "profile", "role", "alias", "file_set", "file_format", "content_type",
    "submitted_file_name", "size_bytes", "md5sum", "analysis_step_version",
    "derived_from", "reference_files", "filtered", "controlled_access",
    "derived_manually", "description", "validation_status", "issues",
]
master_manifest = final[manifest_columns].copy()
master_manifest_path = METADATA_DIR / "upload_manifest_16_files.tsv"
master_manifest.to_csv(master_manifest_path, sep="\t", index=False, lineterminator="\n")

# Build profile-specific iu_register inputs. Array properties are represented as
# comma-separated identifiers for igvf_utils.
matrix_columns = [
    "aliases", "award", "lab", "md5sum", "file_format", "file_set",
    "content_type", "reference_files", "filtered", "analysis_step_version",
    "derived_from", "derived_manually", "description", "submitter_comment",
    "submitted_file_name",
]

tabular_columns = [
    "aliases", "award", "lab", "md5sum", "file_format", "file_set",
    "content_type", "controlled_access", "reference_files", "filtered",
    "analysis_step_version", "derived_from", "derived_manually", "description",
    "submitter_comment", "submitted_file_name",
]

meta = final.rename(columns={"alias": "aliases"}).copy()

h5ad_meta = meta[(meta["profile"].eq("matrix_file")) & (meta["file_format"].eq("h5ad"))]
h5mu_meta = meta[(meta["profile"].eq("matrix_file")) & (meta["file_format"].eq("h5mu"))]
tabular_meta = meta[meta["profile"].eq("tabular_file")]

# Enforce submission order within each file.
h5ad_meta = h5ad_meta.sort_values(["day", "aliases"])
h5mu_meta = h5mu_meta.sort_values(["day", "aliases"])
tabular_meta = tabular_meta.assign(
    _role_order=tabular_meta["role"].map({"guide quantifications": 0, "per-cell quality report": 1})
).sort_values(["day", "_role_order", "aliases"])

h5ad_tsv = METADATA_DIR / "01_matrix_file_h5ad.tsv"
h5mu_tsv = METADATA_DIR / "02_matrix_file_h5mu.tsv"
tabular_tsv = METADATA_DIR / "03_tabular_file_step1.tsv"
combined_matrix_tsv = METADATA_DIR / "matrix_file_all_8.tsv"

h5ad_meta[matrix_columns].to_csv(h5ad_tsv, sep="\t", index=False, lineterminator="\n")
h5mu_meta[matrix_columns].to_csv(h5mu_tsv, sep="\t", index=False, lineterminator="\n")
pd.concat([h5ad_meta, h5mu_meta], ignore_index=True)[matrix_columns].to_csv(
    combined_matrix_tsv, sep="\t", index=False, lineterminator="\n"
)
tabular_meta[tabular_columns].to_csv(tabular_tsv, sep="\t", index=False, lineterminator="\n")

# Plain md5sum-compatible file.
md5_path = METADATA_DIR / "md5sums_16_files.txt"
with md5_path.open("w", encoding="utf-8") as handle:
    for row in final.sort_values(["day", "profile", "role"]).to_dict("records"):
        handle.write(f"{row['md5sum']}  {row['submitted_file_name']}\n")

# Readiness summary.
summary = pd.DataFrame([
    {"metric": "planned_files", "value": len(final)},
    {"metric": "validated_files", "value": int(final["validation_status"].eq("PASS").sum())},
    {"metric": "h5ad_files", "value": int(final["file_format"].eq("h5ad").sum())},
    {"metric": "h5mu_files", "value": int(final["file_format"].eq("h5mu").sum())},
    {"metric": "tabular_files", "value": int(final["file_format"].eq("tsv").sum())},
    {"metric": "unique_aliases", "value": int(final["alias"].nunique())},
    {"metric": "unique_output_paths", "value": int(final["submitted_file_name"].nunique())},
    {"metric": "analysis_sets", "value": int(final["file_set"].nunique())},
])
summary_path = METADATA_DIR / "upload_readiness_summary.tsv"
summary.to_csv(summary_path, sep="\t", index=False, lineterminator="\n")

if (
    len(final) != 16
    or not final["validation_status"].eq("PASS").all()
    or final["alias"].nunique() != 16
    or final["submitted_file_name"].nunique() != 16
    or final["file_set"].nunique() != 4
    or len(h5ad_meta) != 4
    or len(h5mu_meta) != 4
    or len(tabular_meta) != 8
):
    raise RuntimeError("Final readiness invariants failed")

ready_marker = OUTPUT_ROOT / "READY_TO_UPLOAD_16_FILES.txt"
ready_marker.write_text(
    "\n".join([
        "CC-Perturb-seq IGVF analysis-file build: READY",
        f"generated_utc: {now_utc()}",
        "validated_files: 16",
        "h5ad_files: 4",
        "h5mu_files: 4",
        "guide_quantification_files: 4",
        "per_cell_quality_report_files: 4",
        f"master_manifest: {master_manifest_path}",
        f"validation_report: {validation_path}",
        f"h5ad_metadata: {h5ad_tsv}",
        f"h5mu_metadata: {h5mu_tsv}",
        f"tabular_metadata: {tabular_tsv}",
        f"column_dictionary: {column_dictionary_path}",
        "",
    ]),
    encoding="utf-8",
)

print(summary.to_string(index=False))
print("\nREADY marker:", ready_marker)
print("Master manifest:", master_manifest_path)
print("Validation report:", validation_path)
print("H5AD metadata:", h5ad_tsv)
print("H5MU metadata:", h5mu_tsv)
print("Tabular metadata:", tabular_tsv)
print("MD5 list:", md5_path)

```

### Text output (stream)

```text
             metric  value
      planned_files     16
    validated_files     16
         h5ad_files      4
         h5mu_files      4
      tabular_files      8
     unique_aliases     16
unique_output_paths     16
      analysis_sets      4

READY marker: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/READY_TO_UPLOAD_16_FILES.txt
Master manifest: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/metadata/upload_manifest_16_files.tsv
Validation report: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/_prep_work/file_validation_16.tsv
H5AD metadata: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/metadata/01_matrix_file_h5ad.tsv
H5MU metadata: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/metadata/02_matrix_file_h5mu.tsv
Tabular metadata: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/metadata/03_tabular_file_step1.tsv
MD5 list: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/metadata/md5sums_16_files.txt

```

## Markdown cell 18

## 8. Generate—but do not execute—the IGVF dry-run commands

## Code cell 19

```python

IGVF_PYTHON = Path("/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python")
IU_REGISTER = Path(
    "/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/"
    "igvf_utils/MetaDataRegistration/iu_register.py"
)

commands_path = METADATA_DIR / "run_igvf_dry_runs.sh"
commands = "\n".join([
    "#!/bin/bash",
    "set -euo pipefail",
    "",
    ': "${IGVF_API_KEY:<REDACTED> IGVF_API_KEY first}"',
    ': "${IGVF_SECRET_KEY:?Export IGVF_SECRET_KEY first}"',
    "",
    f"PY={IGVF_PYTHON}",
    f"IU={IU_REGISTER}",
    f"META={METADATA_DIR}",
    f"LOGS={OUTPUT_ROOT / 'logs'}",
    'mkdir -p "$LOGS"',
    "",
    "# Stage 1: annotated h5ad MatrixFiles. These depend only on existing raw SequenceFiles.",
    '"$PY" "$IU" -m prod -p matrix_file -d -i "$META/01_matrix_file_h5ad.tsv" \\',
    '  2>&1 | tee "$LOGS/01_matrix_file_h5ad_dryrun_$(date +%Y%m%d_%H%M%S).log"',
    "",
    "# Submit the h5ad records for real only after reviewing the dry run. Then run Stage 2.",
    "# The h5mu and tabular files use h5ad aliases in derived_from.",
    '"$PY" "$IU" -m prod -p matrix_file -d -i "$META/02_matrix_file_h5mu.tsv" \\',
    '  2>&1 | tee "$LOGS/02_matrix_file_h5mu_dryrun_$(date +%Y%m%d_%H%M%S).log"',
    "",
    '"$PY" "$IU" -m prod -p tabular_file -d -i "$META/03_tabular_file_step1.tsv" \\',
    '  2>&1 | tee "$LOGS/03_tabular_file_step1_dryrun_$(date +%Y%m%d_%H%M%S).log"',
    "",
])
commands_path.write_text(commands, encoding="utf-8")
commands_path.chmod(0o750)

print(commands)
print("Saved:", commands_path)
print("No Portal command was executed by the notebook.")

```

### Text output (stream)

```text
#!/bin/bash
set -euo pipefail

: "${IGVF_API_KEY:<REDACTED> IGVF_API_KEY first}"
: "${IGVF_SECRET_KEY:?Export IGVF_SECRET_KEY first}"

PY=/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python
IU=/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py
META=/oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/metadata
LOGS=/oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/logs
mkdir -p "$LOGS"

# Stage 1: annotated h5ad MatrixFiles. These depend only on existing raw SequenceFiles.
"$PY" "$IU" -m prod -p matrix_file -d -i "$META/01_matrix_file_h5ad.tsv" \
  2>&1 | tee "$LOGS/01_matrix_file_h5ad_dryrun_$(date +%Y%m%d_%H%M%S).log"

# Submit the h5ad records for real only after reviewing the dry run. Then run Stage 2.
# The h5mu and tabular files use h5ad aliases in derived_from.
"$PY" "$IU" -m prod -p matrix_file -d -i "$META/02_matrix_file_h5mu.tsv" \
  2>&1 | tee "$LOGS/02_matrix_file_h5mu_dryrun_$(date +%Y%m%d_%H%M%S).log"

"$PY" "$IU" -m prod -p tabular_file -d -i "$META/03_tabular_file_step1.tsv" \
  2>&1 | tee "$LOGS/03_tabular_file_step1_dryrun_$(date +%Y%m%d_%H%M%S).log"

Saved: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/metadata/run_igvf_dry_runs.sh
No Portal command was executed by the notebook.

```

## Markdown cell 20


## Final outputs to inspect

The principal outputs are:

```text
${OUTPUT_ROOT}/READY_TO_UPLOAD_16_FILES.txt
${OUTPUT_ROOT}/metadata/upload_manifest_16_files.tsv
${OUTPUT_ROOT}/_prep_work/file_validation_16.tsv
${OUTPUT_ROOT}/metadata/01_matrix_file_h5ad.tsv
${OUTPUT_ROOT}/metadata/02_matrix_file_h5mu.tsv
${OUTPUT_ROOT}/metadata/03_tabular_file_step1.tsv
${OUTPUT_ROOT}/documentation/column_dictionary_step1.tsv
```

Do not run a real upload until the dry runs are reviewed. Submit the four h5ad records first, because the h5mu and tabular metadata use the h5ad aliases in `derived_from`.


