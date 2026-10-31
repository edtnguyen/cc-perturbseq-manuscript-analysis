# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/submit_300gene_analysis_file_format_documents.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/submit_300gene_analysis_file_format_documents.ipynb`
- SHA-256: `6de413323783e2ae8222b58702ca441926199d6d224fb7eee6ef16c068fb6578`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Create and submit file-format specification Documents for the 300-gene CC-Perturb-seq analysis files

This notebook creates **three lab-specific IGVF file-format specification Documents** and attaches each Document to the four previously submitted day-specific files of the matching type:

1. `guide quantifications` — four D0-D3 `TabularFile` records;
2. `per-cell quality report` — four D0-D3 `TabularFile` records;
3. `annotated cell by gene and guide matrix` — four D0-D3 `MatrixFile` records.

The notebook uses the 12 relevant accessions from `igvf_analysis_file_accessions_16.csv`, with an embedded fallback copy so the CSV does not have to be present on Sherlock. It inspects the actual final files on OAK before writing the specifications, creates the three PDFs, posts them as IGVF `Document` records, patches the 12 existing files through `file_format_specifications`, and verifies the final links.

It does **not** alter or re-upload the analysis files.

## Markdown cell 1

## Version 2 — native JSON Document POSTs and corrected PATCH paths

This version incorporates both fixes validated in the earlier
four-Document workflow:

1. IGVF `Document` records are submitted from native `.json` files,
   so `attachment` remains a real JSON object rather than a
   CSV-escaped string inside a TSV.
2. File PATCH operations continue to use TSV and consistently read
   their path from `input_tsv`.

The notebook removes obsolete `POST_document_*.tsv` files, validates
every JSON attachment before dry run, and validates every PATCH TSV
path before dry run or production.

## Markdown cell 2

## 1. Configuration

The output directory sits beside the existing analysis submission. The notebook invokes the established `igvf_utils` environment as an external process, so the Jupyter kernel can remain in the `igvf-upload-ready` analysis environment.

## Code cell 3

```python
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from getpass import getpass
from pathlib import Path
from typing import Any, Sequence

import h5py
import numpy as np
import pandas as pd
import requests

try:
    import anndata as ad
except ImportError as exc:
    raise ImportError(
        "anndata is required to inspect the submitted h5ad files. "
        "Use the same kernel/environment as the upload-ready builder notebook."
    ) from exc

ANALYSIS_ROOT = Path(
    "/oak/stanford/groups/engreitz/Projects/IGVF/"
    "CC-Perturb-seq/300gene/analysis_submission"
)
OUTPUT_ROOT = ANALYSIS_ROOT / "file_format_specifications"
WORK_DIR = OUTPUT_ROOT / "_work"
LOG_DIR = OUTPUT_ROOT / "logs"
ROW_DIR = WORK_DIR / "one_row_tsv"
PDF_DIR = OUTPUT_ROOT / "documents"
SCHEMA_DIR = OUTPUT_ROOT / "schema_tables"

for directory in [OUTPUT_ROOT, WORK_DIR, LOG_DIR, ROW_DIR, PDF_DIR, SCHEMA_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

# Optional live CSV. If absent, the embedded accession table below is used.
ACCESSION_CSV_CANDIDATES = [
    ANALYSIS_ROOT / "metadata/igvf_analysis_file_accessions_16.csv",
    Path.cwd() / "igvf_analysis_file_accessions_16.csv",
    Path.cwd().parent / "igvf_analysis_file_accessions_16.csv",
]

AWARD = "/awards/HG011972/"
LAB = "/labs/jesse-engreitz/"

IGVF_PYTHON = Path(
    "/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python"
)
IU_REGISTER = Path(
    "/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/"
    "igvf_utils/MetaDataRegistration/iu_register.py"
)

DOCUMENT_SPECS = {
    "guide_quantifications": {
        "document_alias": (
            "jesse-engreitz:document_WTC11-EC_300gene_"
            "guide_quantifications_file_format_v1"
        ),
        "pdf_name": (
            "WTC11_EC_300gene_guide_quantifications_"
            "file_format_specification.pdf"
        ),
        "description": (
            "File-format specification for per-cell guide assignment "
            "quantifications from the WTC11-EC 300-gene CC-Perturb-seq "
            "D0-D3 analysis files."
        ),
        "profile": "tabular_file",
        "role": "guide quantifications",
    },
    "per_cell_quality_report": {
        "document_alias": (
            "jesse-engreitz:document_WTC11-EC_300gene_"
            "per_cell_quality_report_file_format_v1"
        ),
        "pdf_name": (
            "WTC11_EC_300gene_per_cell_quality_report_"
            "file_format_specification.pdf"
        ),
        "description": (
            "File-format specification for per-cell quality-control reports "
            "from the WTC11-EC 300-gene CC-Perturb-seq D0-D3 analysis files."
        ),
        "profile": "tabular_file",
        "role": "per-cell quality report",
    },
    "annotated_h5ad": {
        "document_alias": (
            "jesse-engreitz:document_WTC11-EC_300gene_"
            "annotated_cell_gene_guide_h5ad_file_format_v1"
        ),
        "pdf_name": (
            "WTC11_EC_300gene_annotated_cell_gene_guide_h5ad_"
            "file_format_specification.pdf"
        ),
        "description": (
            "File-format specification for day-specific annotated cell-by-gene-"
            "and-guide AnnData matrices from the WTC11-EC 300-gene "
            "CC-Perturb-seq analysis."
        ),
        "profile": "matrix_file",
        "role": "annotated cell by gene and guide matrix",
    },
}

# No production operation runs until these tokens are entered later.
DOCUMENT_PRODUCTION_CONFIRMATION = ""
PATCH_PRODUCTION_CONFIRMATION = ""

print("Analysis root:", ANALYSIS_ROOT)
print("Output root:", OUTPUT_ROOT)
print("igvf_utils Python:", IGVF_PYTHON)
print("iu_register.py:", IU_REGISTER)
```

### Text output (stream)

```text
Analysis root: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission
Output root: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications
igvf_utils Python: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python
iu_register.py: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py

```

## Markdown cell 4

## 2. Check PDF-authoring dependency

The notebook uses ReportLab for stable multipage PDF generation. Install it once in the current Jupyter environment if the import fails.

## Code cell 5

```python
try:
    import reportlab
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER, TA_LEFT
    from reportlab.lib.pagesizes import landscape, letter
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import inch
    from reportlab.platypus import (
        BaseDocTemplate,
        Frame,
        LongTable,
        PageTemplate,
        Paragraph,
        Spacer,
        Table,
        TableStyle,
    )
except ImportError as exc:
    raise ImportError(
        "ReportLab is required. Install it in the active notebook environment, "
        "for example:\n\n"
        "mamba install -c conda-forge reportlab\n\n"
        "Then restart the kernel and rerun this notebook."
    ) from exc

print("reportlab:", reportlab.Version)
print("anndata:", getattr(ad, "__version__", "unknown"))
print("h5py:", h5py.__version__)
print("pandas:", pd.__version__)
```

### Text output (stream)

```text
reportlab: 5.0.1
anndata: 0.12.6
h5py: 3.15.1
pandas: 2.3.3

```

### Text output (stream)

```text
/tmp/ipykernel_22426/418209516.py:27: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('anndata')` instead.
  print("anndata:", getattr(ad, "__version__", "unknown"))

```

## Markdown cell 6

## 3. Load the 12 target file records

The embedded records are taken from the attached `igvf_analysis_file_accessions_16.csv`. If a live copy of that CSV exists on Sherlock, the notebook uses it and verifies that the accessions and paths agree with the embedded records.

## Code cell 7

```python
EMBEDDED_FILE_RECORDS = [
    {
        "day": "D0", "profile": "matrix_file",
        "alias": "jesse-engreitz:matrix_file_WTC11-EC_300gene_CC-Perturbseq_D0_h5ad",
        "accession": "IGVFFI7756RKHT",
        "submitted_file_name": str(ANALYSIS_ROOT / "D0/cc_perturbseq/CCPerturbSeq_D0_annotated_cell_gene_guide.h5ad"),
        "expected_file_set": "IGVFDS0463XTUH",
        "expected_md5sum": "0d3fe724ad9eca6cb9e811b71e67b2d7",
        "format_role": "annotated_h5ad",
    },
    {
        "day": "D1", "profile": "matrix_file",
        "alias": "jesse-engreitz:matrix_file_WTC11-EC_300gene_CC-Perturbseq_D1_h5ad",
        "accession": "IGVFFI4783YGLZ",
        "submitted_file_name": str(ANALYSIS_ROOT / "D1/cc_perturbseq/CCPerturbSeq_D1_annotated_cell_gene_guide.h5ad"),
        "expected_file_set": "IGVFDS0991DBQZ",
        "expected_md5sum": "8ffe7570d640a4f72b2de1d6d68e0a01",
        "format_role": "annotated_h5ad",
    },
    {
        "day": "D2", "profile": "matrix_file",
        "alias": "jesse-engreitz:matrix_file_WTC11-EC_300gene_CC-Perturbseq_D2_h5ad",
        "accession": "IGVFFI5165BSGC",
        "submitted_file_name": str(ANALYSIS_ROOT / "D2/cc_perturbseq/CCPerturbSeq_D2_annotated_cell_gene_guide.h5ad"),
        "expected_file_set": "IGVFDS4402ZGWJ",
        "expected_md5sum": "f2fe9f4af11a3eda18617ce4b7ec36f6",
        "format_role": "annotated_h5ad",
    },
    {
        "day": "D3", "profile": "matrix_file",
        "alias": "jesse-engreitz:matrix_file_WTC11-EC_300gene_CC-Perturbseq_D3_h5ad",
        "accession": "IGVFFI9086SWZC",
        "submitted_file_name": str(ANALYSIS_ROOT / "D3/cc_perturbseq/CCPerturbSeq_D3_annotated_cell_gene_guide.h5ad"),
        "expected_file_set": "IGVFDS2307REVJ",
        "expected_md5sum": "e8c726fb0e634074384bc8b7cf5a4464",
        "format_role": "annotated_h5ad",
    },
    {
        "day": "D0", "profile": "tabular_file",
        "alias": "jesse-engreitz:tabular_file_WTC11-EC_300gene_D0_guide_quantifications",
        "accession": "IGVFFI1608GIRE",
        "submitted_file_name": str(ANALYSIS_ROOT / "D0/cc_perturbseq/guide_quantifications.tsv.gz"),
        "expected_file_set": "IGVFDS0463XTUH",
        "expected_md5sum": "e3192ea1c4d02ce875782cacbe1ba99b",
        "format_role": "guide_quantifications",
    },
    {
        "day": "D1", "profile": "tabular_file",
        "alias": "jesse-engreitz:tabular_file_WTC11-EC_300gene_D1_guide_quantifications",
        "accession": "IGVFFI6431OXJN",
        "submitted_file_name": str(ANALYSIS_ROOT / "D1/cc_perturbseq/guide_quantifications.tsv.gz"),
        "expected_file_set": "IGVFDS0991DBQZ",
        "expected_md5sum": "668a5bc7d6c42d294a8aedae485bf20d",
        "format_role": "guide_quantifications",
    },
    {
        "day": "D2", "profile": "tabular_file",
        "alias": "jesse-engreitz:tabular_file_WTC11-EC_300gene_D2_guide_quantifications",
        "accession": "IGVFFI1924OSEM",
        "submitted_file_name": str(ANALYSIS_ROOT / "D2/cc_perturbseq/guide_quantifications.tsv.gz"),
        "expected_file_set": "IGVFDS4402ZGWJ",
        "expected_md5sum": "591144d0e9db49a972662e727ba9ac64",
        "format_role": "guide_quantifications",
    },
    {
        "day": "D3", "profile": "tabular_file",
        "alias": "jesse-engreitz:tabular_file_WTC11-EC_300gene_D3_guide_quantifications",
        "accession": "IGVFFI7077LERW",
        "submitted_file_name": str(ANALYSIS_ROOT / "D3/cc_perturbseq/guide_quantifications.tsv.gz"),
        "expected_file_set": "IGVFDS2307REVJ",
        "expected_md5sum": "e0e4a84989f2179d80a91bbd5c45a8c2",
        "format_role": "guide_quantifications",
    },
    {
        "day": "D0", "profile": "tabular_file",
        "alias": "jesse-engreitz:tabular_file_WTC11-EC_300gene_D0_per_cell_quality_report",
        "accession": "IGVFFI3905TBLN",
        "submitted_file_name": str(ANALYSIS_ROOT / "D0/cc_perturbseq/per_cell_quality_report.tsv.gz"),
        "expected_file_set": "IGVFDS0463XTUH",
        "expected_md5sum": "ea3ca228602d7f433855184555d855b3",
        "format_role": "per_cell_quality_report",
    },
    {
        "day": "D1", "profile": "tabular_file",
        "alias": "jesse-engreitz:tabular_file_WTC11-EC_300gene_D1_per_cell_quality_report",
        "accession": "IGVFFI0627NEZT",
        "submitted_file_name": str(ANALYSIS_ROOT / "D1/cc_perturbseq/per_cell_quality_report.tsv.gz"),
        "expected_file_set": "IGVFDS0991DBQZ",
        "expected_md5sum": "819db5d7e0e81959708efaf4c550fc1e",
        "format_role": "per_cell_quality_report",
    },
    {
        "day": "D2", "profile": "tabular_file",
        "alias": "jesse-engreitz:tabular_file_WTC11-EC_300gene_D2_per_cell_quality_report",
        "accession": "IGVFFI1451CYGW",
        "submitted_file_name": str(ANALYSIS_ROOT / "D2/cc_perturbseq/per_cell_quality_report.tsv.gz"),
        "expected_file_set": "IGVFDS4402ZGWJ",
        "expected_md5sum": "3f1fc50882a9b6cb469ece4914933f94",
        "format_role": "per_cell_quality_report",
    },
    {
        "day": "D3", "profile": "tabular_file",
        "alias": "jesse-engreitz:tabular_file_WTC11-EC_300gene_D3_per_cell_quality_report",
        "accession": "IGVFFI8282ATZQ",
        "submitted_file_name": str(ANALYSIS_ROOT / "D3/cc_perturbseq/per_cell_quality_report.tsv.gz"),
        "expected_file_set": "IGVFDS2307REVJ",
        "expected_md5sum": "ad55846b58d783c7394e6338a104ac2f",
        "format_role": "per_cell_quality_report",
    },
]

embedded = pd.DataFrame(EMBEDDED_FILE_RECORDS)

live_csv = next((path for path in ACCESSION_CSV_CANDIDATES if path.is_file()), None)
if live_csv is not None:
    live = pd.read_csv(live_csv, dtype=str, keep_default_na=False)
    relevant = live.loc[
        live["accession"].isin(embedded["accession"])
    ].copy()
    if len(relevant) != 12:
        raise RuntimeError(
            f"Live accession CSV contains {len(relevant)} of the expected 12 target files: {live_csv}"
        )
    check = embedded.merge(
        relevant[["accession", "alias", "submitted_file_name"]],
        on="accession",
        how="left",
        suffixes=("_embedded", "_live"),
        validate="one_to_one",
    )
    disagreements = check.loc[
        check["alias_embedded"].ne(check["alias_live"])
        | check["submitted_file_name_embedded"].ne(check["submitted_file_name_live"])
    ]
    if not disagreements.empty:
        raise RuntimeError(
            "The live accession CSV disagrees with the embedded records. Review:\n"
            + disagreements.to_string(index=False)
        )
    target_files = embedded.copy()
    source_label = str(live_csv)
else:
    target_files = embedded.copy()
    source_label = "embedded fallback from attached igvf_analysis_file_accessions_16.csv"

if len(target_files) != 12:
    raise RuntimeError("Expected exactly 12 target files")
if target_files["accession"].nunique() != 12:
    raise RuntimeError("Target accessions are not unique")
if set(target_files.groupby("format_role").size()) != {4}:
    raise RuntimeError("Expected four files for each of the three format roles")

TARGET_FILE_INVENTORY_PATH = WORK_DIR / "target_file_inventory.tsv"
target_files.to_csv(TARGET_FILE_INVENTORY_PATH, sep="\t", index=False, lineterminator="\n")

print("Record source:", source_label)
print(target_files[["day", "format_role", "accession", "alias"]].to_string(index=False))
print("\nSaved:", TARGET_FILE_INVENTORY_PATH)
```

### Text output (stream)

```text
Record source: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/metadata/igvf_analysis_file_accessions_16.csv
day             format_role      accession                                                                   alias
 D0          annotated_h5ad IGVFFI7756RKHT       jesse-engreitz:matrix_file_WTC11-EC_300gene_CC-Perturbseq_D0_h5ad
 D1          annotated_h5ad IGVFFI4783YGLZ       jesse-engreitz:matrix_file_WTC11-EC_300gene_CC-Perturbseq_D1_h5ad
 D2          annotated_h5ad IGVFFI5165BSGC       jesse-engreitz:matrix_file_WTC11-EC_300gene_CC-Perturbseq_D2_h5ad
 D3          annotated_h5ad IGVFFI9086SWZC       jesse-engreitz:matrix_file_WTC11-EC_300gene_CC-Perturbseq_D3_h5ad
 D0   guide_quantifications IGVFFI1608GIRE   jesse-engreitz:tabular_file_WTC11-EC_300gene_D0_guide_quantifications
 D1   guide_quantifications IGVFFI6431OXJN   jesse-engreitz:tabular_file_WTC11-EC_300gene_D1_guide_quantifications
 D2   guide_quantifications IGVFFI1924OSEM   jesse-engreitz:tabular_file_WTC11-EC_300gene_D2_guide_quantifications
 D3   guide_quantifications IGVFFI7077LERW   jesse-engreitz:tabular_file_WTC11-EC_300gene_D3_guide_quantifications
 D0 per_cell_quality_report IGVFFI3905TBLN jesse-engreitz:tabular_file_WTC11-EC_300gene_D0_per_cell_quality_report
 D1 per_cell_quality_report IGVFFI0627NEZT jesse-engreitz:tabular_file_WTC11-EC_300gene_D1_per_cell_quality_report
 D2 per_cell_quality_report IGVFFI1451CYGW jesse-engreitz:tabular_file_WTC11-EC_300gene_D2_per_cell_quality_report
 D3 per_cell_quality_report IGVFFI8282ATZQ jesse-engreitz:tabular_file_WTC11-EC_300gene_D3_per_cell_quality_report

Saved: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/target_file_inventory.tsv

```

## Markdown cell 8

## 4. Validate local files and inspect their schemas

This section verifies every local file against the expected MD5, confirms that the four files within each role share one schema, and records the actual h5ad organization used to generate the matrix-format specification.

## Code cell 9

```python
def md5sum(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv_gz(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        dtype=str,
        keep_default_na=False,
        na_values=[],
        low_memory=False,
    )


def shape_of(value: Any) -> str:
    shape = getattr(value, "shape", None)
    return "x".join(map(str, shape)) if shape is not None else ""


local_audit_rows: list[dict] = []
for record in target_files.to_dict("records"):
    path = Path(record["submitted_file_name"])
    exists = path.is_file()
    observed_md5 = md5sum(path) if exists else ""
    local_audit_rows.append(
        {
            **record,
            "exists": exists,
            "size_bytes": path.stat().st_size if exists else 0,
            "observed_md5sum": observed_md5,
            "md5_matches": observed_md5 == record["expected_md5sum"],
        }
    )

local_audit = pd.DataFrame(local_audit_rows)
LOCAL_FILE_AUDIT_PATH = WORK_DIR / "local_file_audit.tsv"
local_audit.to_csv(LOCAL_FILE_AUDIT_PATH, sep="\t", index=False, lineterminator="\n")

if not local_audit["exists"].all():
    missing = local_audit.loc[~local_audit["exists"], "submitted_file_name"].tolist()
    raise FileNotFoundError("Missing local files:\n" + "\n".join(missing))
if not local_audit["md5_matches"].all():
    raise RuntimeError(
        "At least one local file no longer matches the uploaded MD5. Review "
        f"{LOCAL_FILE_AUDIT_PATH}"
    )

# Tabular schemas.
tabular_schema_rows: list[dict] = []
tabular_headers: dict[str, set[tuple[str, ...]]] = {
    "guide_quantifications": set(),
    "per_cell_quality_report": set(),
}

for role in tabular_headers:
    for record in target_files.loc[target_files["format_role"].eq(role)].to_dict("records"):
        table = read_tsv_gz(Path(record["submitted_file_name"]))
        header = tuple(table.columns.astype(str))
        tabular_headers[role].add(header)
        if table.empty:
            raise RuntimeError(f"Empty table: {record['submitted_file_name']}")
        for order, column in enumerate(header, start=1):
            values = table[column].astype(str).str.strip()
            numeric = pd.to_numeric(values, errors="coerce")
            if values.ne("").all() and numeric.notna().all():
                data_type = (
                    "integer"
                    if np.allclose(numeric.to_numpy(float), np.round(numeric.to_numpy(float)))
                    else "float"
                )
            elif set(values.str.upper().unique()) <= {"TRUE", "FALSE"}:
                data_type = "boolean string"
            else:
                data_type = "string"
            tabular_schema_rows.append(
                {
                    "format_role": role,
                    "day": record["day"],
                    "column_order": order,
                    "column_name": column,
                    "inferred_data_type": data_type,
                    "missing_values": int(values.eq("").sum()),
                    "unique_values": int(values.nunique(dropna=False)),
                    "row_count": len(table),
                }
            )

for role, headers in tabular_headers.items():
    if len(headers) != 1:
        raise RuntimeError(f"{role} does not have one shared D0-D3 schema: {headers}")

TABULAR_SCHEMA_AUDIT_PATH = SCHEMA_DIR / "tabular_schema_audit.tsv"
tabular_schema_audit = pd.DataFrame(tabular_schema_rows)
tabular_schema_audit.to_csv(
    TABULAR_SCHEMA_AUDIT_PATH, sep="\t", index=False, lineterminator="\n"
)

# H5AD structure audit.
h5ad_rows: list[dict] = []
h5ad_component_rows: list[dict] = []
h5ad_schema_signatures: list[dict] = []

for record in target_files.loc[target_files["format_role"].eq("annotated_h5ad")].to_dict("records"):
    path = Path(record["submitted_file_name"])
    with h5py.File(path, "r") as handle:
        top_level_groups = sorted(handle.keys())
    adata = ad.read_h5ad(path, backed="r")
    try:
        obs_columns = list(map(str, adata.obs.columns))
        var_columns = list(map(str, adata.var.columns))
        obsm_keys = sorted(map(str, adata.obsm.keys()))
        layers = sorted(map(str, adata.layers.keys()))
        obsp_keys = sorted(map(str, adata.obsp.keys()))
        varm_keys = sorted(map(str, adata.varm.keys()))
        varp_keys = sorted(map(str, adata.varp.keys()))
        uns_keys = sorted(map(str, adata.uns.keys()))

        assignment_key = next(
            (key for key in ("guide_assignment", "guide_assignments") if key in adata.obsm),
            None,
        )
        guide_names_key = next(
            (key for key in ("guide_names", "guide_ids", "guides", "gRNA_names") if key in adata.uns),
            None,
        )
        if assignment_key is None:
            raise RuntimeError(f"No guide assignment obsm key in {path}")
        if guide_names_key is None:
            raise RuntimeError(f"No guide-name uns key in {path}")

        assignment_shape = getattr(adata.obsm[assignment_key], "shape", None)
        guide_name_count = len(adata.uns[guide_names_key])
        if assignment_shape[0] != adata.n_obs or assignment_shape[1] != guide_name_count:
            raise RuntimeError(f"Guide assignment shape mismatch in {path}")
        if not adata.obs_names.is_unique or not adata.var_names.is_unique:
            raise RuntimeError(f"Non-unique AnnData axis names in {path}")

        h5ad_rows.append(
            {
                "day": record["day"],
                "accession": record["accession"],
                "path": str(path),
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "x_shape": shape_of(adata.X),
                "raw_present": adata.raw is not None,
                "obs_columns": ",".join(obs_columns),
                "var_columns": ",".join(var_columns),
                "obsm_keys": ",".join(obsm_keys),
                "layers": ",".join(layers),
                "obsp_keys": ",".join(obsp_keys),
                "varm_keys": ",".join(varm_keys),
                "varp_keys": ",".join(varp_keys),
                "uns_keys": ",".join(uns_keys),
                "top_level_hdf5_groups": ",".join(top_level_groups),
                "guide_assignment_key": assignment_key,
                "guide_assignment_shape": "x".join(map(str, assignment_shape)),
                "guide_names_key": guide_names_key,
                "guide_name_count": guide_name_count,
            }
        )

        for axis, frame in [("obs", adata.obs), ("var", adata.var)]:
            for order, column in enumerate(frame.columns, start=1):
                h5ad_component_rows.append(
                    {
                        "component": axis,
                        "key": str(column),
                        "order": order,
                        "dtype": str(frame[column].dtype),
                        "shape": str(len(frame)),
                        "day": record["day"],
                    }
                )
        for component, mapping in [
            ("obsm", adata.obsm),
            ("layers", adata.layers),
            ("obsp", adata.obsp),
            ("varm", adata.varm),
            ("varp", adata.varp),
        ]:
            for key in mapping.keys():
                h5ad_component_rows.append(
                    {
                        "component": component,
                        "key": str(key),
                        "order": "",
                        "dtype": str(getattr(mapping[key], "dtype", type(mapping[key]).__name__)),
                        "shape": shape_of(mapping[key]),
                        "day": record["day"],
                    }
                )
        for key in adata.uns.keys():
            value = adata.uns[key]
            h5ad_component_rows.append(
                {
                    "component": "uns",
                    "key": str(key),
                    "order": "",
                    "dtype": type(value).__name__,
                    "shape": shape_of(value),
                    "day": record["day"],
                }
            )

        h5ad_schema_signatures.append(
            {
                "day": record["day"],
                "obs_columns": tuple(obs_columns),
                "var_columns": tuple(var_columns),
                "obsm_keys": tuple(obsm_keys),
                "layers": tuple(layers),
                "obsp_keys": tuple(obsp_keys),
                "varm_keys": tuple(varm_keys),
                "varp_keys": tuple(varp_keys),
                "guide_assignment_key": assignment_key,
                "guide_names_key": guide_names_key,
            }
        )
    finally:
        try:
            adata.file.close()
        except Exception:
            pass

H5AD_FILE_AUDIT_PATH = SCHEMA_DIR / "h5ad_file_audit.tsv"
pd.DataFrame(h5ad_rows).to_csv(
    H5AD_FILE_AUDIT_PATH, sep="\t", index=False, lineterminator="\n"
)
H5AD_COMPONENT_AUDIT_PATH = SCHEMA_DIR / "h5ad_component_audit.tsv"
pd.DataFrame(h5ad_component_rows).to_csv(
    H5AD_COMPONENT_AUDIT_PATH, sep="\t", index=False, lineterminator="\n"
)

signature_frame = pd.DataFrame(h5ad_schema_signatures)
for column in [
    "obs_columns", "var_columns", "obsm_keys", "layers", "obsp_keys",
    "varm_keys", "varp_keys", "guide_assignment_key", "guide_names_key",
]:
    if signature_frame[column].astype(str).nunique() != 1:
        raise RuntimeError(
            f"The four h5ad files differ in {column}. Review {H5AD_FILE_AUDIT_PATH}"
        )

print("All 12 local files match the previously submitted MD5 values.")
print("Guide schema:", next(iter(tabular_headers["guide_quantifications"])))
print("QC schema:", next(iter(tabular_headers["per_cell_quality_report"])))
print("H5AD audit:", H5AD_FILE_AUDIT_PATH)
print("H5AD component audit:", H5AD_COMPONENT_AUDIT_PATH)
```

### Text output (stream)

```text
All 12 local files match the previously submitted MD5 values.
Guide schema: ('cell_id', 'guide_id', 'assignment_value')
QC schema: ('cell_id', 'day', 'biological_sample', 'total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'guides_per_cell', 'guide_umi_counts')
H5AD audit: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/schema_tables/h5ad_file_audit.tsv
H5AD component audit: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/schema_tables/h5ad_component_audit.tsv

```

## Markdown cell 10

## 5. Generate the three PDF specifications and machine-readable schema tables

One Document is shared by all four D0-D3 files of the same type. The h5ad specification is generated from the actual common structure of the four uploaded AnnData files.

## Code cell 11

```python
# Column descriptions for the two tabular formats.
GUIDE_DESCRIPTIONS = {
    "cell_id": "Unique cell identifier matching obs_names in the corresponding day-specific annotated h5ad file.",
    "guide_id": "Guide RNA identifier corresponding to a guide name in the annotated h5ad guide-assignment matrix.",
    "assignment_value": "Nonzero guide-assignment value for the cell-guide pair, exported from the AnnData guide-assignment matrix.",
    "guide_umi_count": "Guide UMI count for the assigned cell-guide pair when a matching guide UMI matrix is present in the AnnData object.",
}

QC_DESCRIPTIONS = {
    "cell_id": "Unique cell identifier matching obs_names in the corresponding day-specific annotated h5ad file.",
    "day": "Differentiation day encoded as d0, d1, d2, or d3; every file contains only its designated day.",
    "biological_sample": "Biological sample or replicate annotation carried from the AnnData obs table.",
    "total_counts": "Total transcript count summary for the cell as stored in AnnData obs.",
    "n_genes_by_counts": "Number of genes with observed counts for the cell as stored in AnnData obs.",
    "pct_counts_mt": "Percentage of the cell's transcript counts assigned to mitochondrial genes.",
    "guides_per_cell": "Number of assigned guides for the cell.",
    "guide_umi_counts": "Total guide UMI count summary for the cell.",
}


def representative_schema(role: str) -> pd.DataFrame:
    rows = tabular_schema_audit.loc[tabular_schema_audit["format_role"].eq(role)].copy()
    first_day = sorted(rows["day"].unique())[0]
    return rows.loc[rows["day"].eq(first_day)].sort_values("column_order")


def tabular_spec_frame(role: str, descriptions: dict[str, str]) -> pd.DataFrame:
    schema = representative_schema(role)
    output = schema[["column_order", "column_name", "inferred_data_type"]].copy()
    output["required"] = "yes"
    output["missing_allowed"] = "no"
    output["description"] = output["column_name"].map(descriptions).fillna(
        "Column observed in the final submitted files; interpretation follows the preprocessing workflow."
    )
    return output


guide_spec = tabular_spec_frame("guide_quantifications", GUIDE_DESCRIPTIONS)
qc_spec = tabular_spec_frame("per_cell_quality_report", QC_DESCRIPTIONS)

GUIDE_SCHEMA_PATH = SCHEMA_DIR / "guide_quantifications_column_specification.tsv"
QC_SCHEMA_PATH = SCHEMA_DIR / "per_cell_quality_report_column_specification.tsv"
guide_spec.to_csv(GUIDE_SCHEMA_PATH, sep="\t", index=False, lineterminator="\n")
qc_spec.to_csv(QC_SCHEMA_PATH, sep="\t", index=False, lineterminator="\n")

# Common h5ad schema from the D0 exemplar after cross-day consistency checks.
common_components = pd.DataFrame(h5ad_component_rows)
common_components = common_components.loc[common_components["day"].eq("D0")].copy()
H5AD_SCHEMA_PATH = SCHEMA_DIR / "annotated_h5ad_component_specification.tsv"
common_components.to_csv(H5AD_SCHEMA_PATH, sep="\t", index=False, lineterminator="\n")

# ReportLab styles.
PAGE = landscape(letter)
LEFT = 0.45 * inch
RIGHT = 0.45 * inch
TOP = 0.50 * inch
BOTTOM = 0.45 * inch
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(
    name="SpecTitle", parent=styles["Title"], fontName="Helvetica-Bold",
    fontSize=18, leading=22, alignment=TA_CENTER,
    textColor=colors.HexColor("#1f2937"), spaceAfter=9,
))
styles.add(ParagraphStyle(
    name="SpecSubtitle", parent=styles["Normal"], fontName="Helvetica",
    fontSize=9, leading=12, alignment=TA_CENTER,
    textColor=colors.HexColor("#4b5563"), spaceAfter=10,
))
styles.add(ParagraphStyle(
    name="SpecH1", parent=styles["Heading1"], fontName="Helvetica-Bold",
    fontSize=12.5, leading=15, textColor=colors.HexColor("#1f4e79"),
    spaceBefore=7, spaceAfter=5,
))
styles.add(ParagraphStyle(
    name="SpecBody", parent=styles["BodyText"], fontName="Helvetica",
    fontSize=8.2, leading=11, spaceAfter=3,
))
styles.add(ParagraphStyle(
    name="Cell", parent=styles["BodyText"], fontName="Helvetica",
    fontSize=6.5, leading=8, splitLongWords=True,
))
styles.add(ParagraphStyle(
    name="CellHead", parent=styles["BodyText"], fontName="Helvetica-Bold",
    fontSize=6.6, leading=8, textColor=colors.white,
))


def esc(value: Any) -> str:
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def para(value: Any, style: str = "SpecBody") -> Paragraph:
    return Paragraph(esc(value), styles[style])


def bullet(value: str) -> Paragraph:
    return Paragraph("&#8226;&nbsp;" + esc(value), styles["SpecBody"])


def footer(canvas, doc) -> None:
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#d1d5db"))
    canvas.setLineWidth(0.4)
    canvas.line(LEFT, 0.33 * inch, PAGE[0] - RIGHT, 0.33 * inch)
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(colors.HexColor("#6b7280"))
    canvas.drawString(LEFT, 0.19 * inch, "Engreitz Lab - IGVF file format specification")
    canvas.drawRightString(PAGE[0] - RIGHT, 0.19 * inch, f"Page {doc.page}")
    canvas.restoreState()


def table_style() -> TableStyle:
    return TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1f4e79")),
        ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#c7d2de")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7fafc")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 3.5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3.5),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ])


def build_document(path: Path, title: str, subtitle: str, story_items: list[Any]) -> None:
    doc = BaseDocTemplate(
        str(path), pagesize=PAGE, leftMargin=LEFT, rightMargin=RIGHT,
        topMargin=TOP, bottomMargin=BOTTOM, title=title,
        author="Engreitz Lab", subject="IGVF file format specification",
    )
    frame = Frame(
        LEFT, BOTTOM, PAGE[0] - LEFT - RIGHT, PAGE[1] - TOP - BOTTOM,
        id="normal",
    )
    doc.addPageTemplates([PageTemplate(id="main", frames=[frame], onPage=footer)])
    story = [para(title, "SpecTitle"), para(subtitle, "SpecSubtitle"), *story_items]
    doc.build(story)


def build_tabular_pdf(role_key: str, spec_frame: pd.DataFrame, title: str, row_unit: str, rules: list[str]) -> Path:
    pdf_path = PDF_DIR / DOCUMENT_SPECS[role_key]["pdf_name"]
    target_rows = target_files.loc[target_files["format_role"].eq(role_key)]
    overview_data = [
        [para("Property", "CellHead"), para("Specification", "CellHead")],
        [para("Applicable files", "Cell"), para(", ".join(target_rows["accession"]), "Cell")],
        [para("IGVF content type", "Cell"), para(DOCUMENT_SPECS[role_key]["role"], "Cell")],
        [para("Physical serialization", "Cell"), para("UTF-8 tab-separated text with one header row, deterministic gzip compression (.tsv.gz), and no comment/preamble lines.", "Cell")],
        [para("Row unit", "Cell"), para(row_unit, "Cell")],
        [para("D0-D3 schema", "Cell"), para("The four day-specific files use the same ordered columns; row counts vary by day.", "Cell")],
    ]
    overview = Table(overview_data, colWidths=[1.75 * inch, 8.45 * inch], repeatRows=1)
    overview.setStyle(table_style())

    headers = ["Order", "Column", "Type", "Required", "Missing", "Definition"]
    rows = [[para(value, "CellHead") for value in headers]]
    for record in spec_frame.to_dict("records"):
        rows.append([
            para(record["column_order"], "Cell"), para(record["column_name"], "Cell"),
            para(record["inferred_data_type"], "Cell"), para(record["required"], "Cell"),
            para(record["missing_allowed"], "Cell"), para(record["description"], "Cell"),
        ])
    schema_table = LongTable(
        rows,
        colWidths=[0.45 * inch, 1.45 * inch, 0.85 * inch, 0.68 * inch, 0.68 * inch, 6.1 * inch],
        repeatRows=1,
    )
    schema_table.setStyle(table_style())

    story = [
        para("1. Scope and serialization", "SpecH1"), overview, Spacer(1, 5),
        para("2. Ordered column specification", "SpecH1"), schema_table,
        para("3. Validation rules", "SpecH1"),
        *[bullet(rule) for rule in rules],
        para("4. IGVF linkage", "SpecH1"),
        bullet("Submit this PDF as an IGVF Document with document_type = file format specification."),
        bullet("Attach the Document to each applicable File record through file_format_specifications."),
        bullet("standardized_file_format is intentionally omitted because these are lab-specific formats and the field is admin-only."),
    ]
    build_document(
        pdf_path,
        title,
        "Lab-specific file format specification, version 1.0 | WTC11-EC 300-gene CC-Perturb-seq",
        story,
    )
    return pdf_path


guide_pdf = build_tabular_pdf(
    "guide_quantifications", guide_spec,
    "Guide quantifications",
    "One nonzero cell-guide assignment per row; (cell_id, guide_id) is unique.",
    [
        "cell_id values are a subset of the corresponding day-specific h5ad obs_names.",
        "Every (cell_id, guide_id) pair is unique; the table contains only assigned/nonzero guide pairs.",
        "assignment_value is numeric; when integer-valued it is serialized as an integer.",
        "guide_umi_count is included only when the source h5ad contains a guide-UMI matrix aligned to the guide-assignment matrix.",
        "No missing values are permitted in columns present in the final schema.",
    ],
)

qc_pdf = build_tabular_pdf(
    "per_cell_quality_report", qc_spec,
    "Per-cell quality report",
    "One cell per row; cell_id is unique.",
    [
        "cell_id exactly matches the corresponding day-specific h5ad obs_names and occurs once.",
        "The day column contains only the designated file day (d0, d1, d2, or d3).",
        "total_counts, n_genes_by_counts, pct_counts_mt, guides_per_cell, and guide_umi_counts are numeric and non-missing.",
        "biological_sample is non-missing.",
        "The number of rows equals n_obs in the corresponding annotated h5ad file.",
    ],
)

# H5AD PDF.
h5ad_pdf_path = PDF_DIR / DOCUMENT_SPECS["annotated_h5ad"]["pdf_name"]
h5ad_targets = target_files.loc[target_files["format_role"].eq("annotated_h5ad")]
h5ad_audit = pd.DataFrame(h5ad_rows).sort_values("day")

h5ad_overview = [
    [para("Property", "CellHead"), para("Specification", "CellHead")],
    [para("Applicable files", "Cell"), para(", ".join(h5ad_targets["accession"]), "Cell")],
    [para("File format", "Cell"), para("h5ad (HDF5 serialization of an AnnData object)", "Cell")],
    [para("IGVF content type", "Cell"), para("annotated cell by gene and guide matrix", "Cell")],
    [para("Observation axis", "Cell"), para("Cells; obs_names are unique cell identifiers.", "Cell")],
    [para("Variable axis", "Cell"), para("Genes/features; var_names are unique identifiers.", "Cell")],
    [para("Primary matrix", "Cell"), para("X is the cell-by-gene numerical matrix stored by the preprocessing workflow. The transformation/normalization semantics are defined by the workflow and are not inferred by this format document.", "Cell")],
    [para("Guide representation", "Cell"), para(f"obsm['{h5ad_audit.iloc[0]['guide_assignment_key']}'] has shape n_obs x n_guides; uns['{h5ad_audit.iloc[0]['guide_names_key']}'] supplies the ordered guide identifiers.", "Cell")],
    [para("Day partitioning", "Cell"), para("Each file contains only one of D0, D1, D2, or D3. Matrix dimensions vary by day; schema and feature annotation structure are shared.", "Cell")],
]
overview_table = Table(h5ad_overview, colWidths=[1.75 * inch, 8.45 * inch], repeatRows=1)
overview_table.setStyle(table_style())

shape_rows = [[para(value, "CellHead") for value in ["Day", "Accession", "n_obs", "n_vars", "X shape", "Guide assignment shape", "Guide count"]]]
for record in h5ad_audit.to_dict("records"):
    shape_rows.append([
        para(record["day"], "Cell"), para(record["accession"], "Cell"),
        para(record["n_obs"], "Cell"), para(record["n_vars"], "Cell"),
        para(record["x_shape"], "Cell"), para(record["guide_assignment_shape"], "Cell"),
        para(record["guide_name_count"], "Cell"),
    ])
shape_table = Table(
    shape_rows,
    colWidths=[0.55 * inch, 1.25 * inch, 0.75 * inch, 0.75 * inch, 1.15 * inch, 1.55 * inch, 0.9 * inch],
    repeatRows=1,
)
shape_table.setStyle(table_style())

component_headers = ["Component", "Key/column", "Order", "dtype/class", "Shape"]
component_rows = [[para(value, "CellHead") for value in component_headers]]
for record in common_components.sort_values(["component", "order", "key"], na_position="last").to_dict("records"):
    component_rows.append([
        para(record["component"], "Cell"), para(record["key"], "Cell"),
        para(record["order"], "Cell"), para(record["dtype"], "Cell"), para(record["shape"], "Cell"),
    ])
component_table = LongTable(
    component_rows,
    colWidths=[0.8 * inch, 3.3 * inch, 0.55 * inch, 1.65 * inch, 1.2 * inch],
    repeatRows=1,
)
component_table.setStyle(table_style())

h5ad_story = [
    para("1. AnnData organization", "SpecH1"), overview_table, Spacer(1, 5),
    para("2. Current day-specific dimensions", "SpecH1"), shape_table,
    para("3. Common stored components", "SpecH1"), component_table,
    para("4. Integrity rules", "SpecH1"),
    bullet("The HDF5 container opens successfully and contains a nonempty AnnData object."),
    bullet("n_obs and n_vars are positive; obs_names and var_names are unique."),
    bullet("The day annotation resolves to exactly one expected day for each file."),
    bullet("The guide-assignment matrix has n_obs rows and one column per ordered guide identifier in uns."),
    bullet("All four day-specific files have the same obs/var schema and the same named multidimensional components; only row counts and day-specific values vary."),
    bullet("Additional AnnData keys may be present only when recorded in the common-component table generated from the submitted files."),
    para("5. IGVF linkage", "SpecH1"),
    bullet("Submit this PDF as an IGVF Document with document_type = file format specification."),
    bullet("Attach the Document to all four MatrixFile records through file_format_specifications."),
    bullet("standardized_file_format is omitted because this is a lab-specific h5ad structure and the field is admin-only."),
]
build_document(
    h5ad_pdf_path,
    "Annotated cell-by-gene-and-guide AnnData matrix",
    "Lab-specific h5ad format specification, version 1.0 | WTC11-EC 300-gene CC-Perturb-seq",
    h5ad_story,
)

pdf_rows = []
for role_key, pdf_path in [
    ("guide_quantifications", guide_pdf),
    ("per_cell_quality_report", qc_pdf),
    ("annotated_h5ad", h5ad_pdf_path),
]:
    if not pdf_path.is_file() or pdf_path.stat().st_size == 0:
        raise RuntimeError(f"PDF was not created: {pdf_path}")
    with pdf_path.open("rb") as handle:
        if handle.read(5) != b"%PDF-":
            raise RuntimeError(f"Invalid PDF header: {pdf_path}")
    pdf_rows.append(
        {
            "format_role": role_key,
            "document_alias": DOCUMENT_SPECS[role_key]["document_alias"],
            "pdf_path": str(pdf_path),
            "size_bytes": pdf_path.stat().st_size,
            "md5sum": md5sum(pdf_path),
        }
    )

PDF_INVENTORY_PATH = WORK_DIR / "pdf_inventory.tsv"
pdf_inventory = pd.DataFrame(pdf_rows)
pdf_inventory.to_csv(PDF_INVENTORY_PATH, sep="\t", index=False, lineterminator="\n")

print(pdf_inventory.to_string(index=False))
print("\nPDF directory:", PDF_DIR)
```

### Text output (stream)

```text
            format_role                                                                         document_alias                                                                                                                                                                                                  pdf_path  size_bytes                           md5sum
  guide_quantifications          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1          /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/documents/WTC11_EC_300gene_guide_quantifications_file_format_specification.pdf        3700 7436e57102fbe13a44798c8c0a2dd845
per_cell_quality_report        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1        /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/documents/WTC11_EC_300gene_per_cell_quality_report_file_format_specification.pdf        4143 9ba3be82cc7db29675b9a017d0b390e1
         annotated_h5ad jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1 /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/documents/WTC11_EC_300gene_annotated_cell_gene_guide_h5ad_file_format_specification.pdf        9919 05473b7f9f7f0148be780ce4be0c2f09

PDF directory: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/documents

```

## Markdown cell 12

## 6. Build the three IGVF Document records as native JSON

Each PDF is represented as a schema-native JSON `Document` payload.
`attachment` is an object of the form
`{"path": "/absolute/path/file.pdf"}`.

A flattened TSV is written only for human review and is never passed
to `iu_register.py`.

## Code cell 13

```python
document_payloads: list[dict] = []

for role_key, spec in DOCUMENT_SPECS.items():
    pdf_path = PDF_DIR / spec["pdf_name"]

    payload = {
        "aliases": [spec["document_alias"]],
        "document_type": "file format specification",
        "description": spec["description"],
        "attachment": {"path": str(pdf_path)},
        "lab": LAB,
        "award": AWARD,
    }

    # Validate the schema-native JSON types immediately.
    payload = json.loads(
        json.dumps(
            payload,
            ensure_ascii=False,
            separators=(",", ":"),
        )
    )

    if not isinstance(payload["aliases"], list):
        raise RuntimeError(
            f"aliases is not an array for {role_key}"
        )
    if not isinstance(payload["attachment"], dict):
        raise RuntimeError(
            f"attachment is not an object for {role_key}"
        )
    if set(payload["attachment"]) != {"path"}:
        raise RuntimeError(
            f"attachment must contain only `path` for {role_key}"
        )

    attachment_path = Path(
        payload["attachment"]["path"]
    )
    if not attachment_path.is_file():
        raise FileNotFoundError(
            f"Document attachment does not exist: {attachment_path}"
        )

    document_payloads.append(payload)

# Aggregate native JSON for audit/recovery.
DOCUMENT_METADATA_PATH = (
    WORK_DIR / "document_file_format_specifications.json"
)
DOCUMENT_METADATA_PATH.write_text(
    json.dumps(
        document_payloads,
        indent=2,
        ensure_ascii=False,
    )
    + "\n",
    encoding="utf-8",
)

reloaded_payloads = json.loads(
    DOCUMENT_METADATA_PATH.read_text(
        encoding="utf-8"
    )
)
if reloaded_payloads != document_payloads:
    raise RuntimeError(
        "Document metadata changed during JSON disk round-trip"
    )

# Human-readable review table only. Do not submit this TSV.
review_rows = [
    {
        "aliases": ",".join(payload["aliases"]),
        "document_type": payload["document_type"],
        "description": payload["description"],
        "attachment_path": payload["attachment"]["path"],
        "lab": payload["lab"],
        "award": payload["award"],
    }
    for payload in document_payloads
]

DOCUMENT_REVIEW_TSV_PATH = (
    WORK_DIR
    / "document_file_format_specifications_review_only.tsv"
)
pd.DataFrame(review_rows).to_csv(
    DOCUMENT_REVIEW_TSV_PATH,
    sep="\t",
    index=False,
    lineterminator="\n",
)

print("Submission JSON:")
print(DOCUMENT_METADATA_PATH)
print(
    json.dumps(
        document_payloads,
        indent=2,
        ensure_ascii=False,
    )
)
print("\nReview-only TSV:")
print(DOCUMENT_REVIEW_TSV_PATH)
```

### Text output (stream)

```text
Submission JSON:
/oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/document_file_format_specifications.json
[
  {
    "aliases": [
      "jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1"
    ],
    "document_type": "file format specification",
    "description": "File-format specification for per-cell guide assignment quantifications from the WTC11-EC 300-gene CC-Perturb-seq D0-D3 analysis files.",
    "attachment": {
      "path": "/oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/documents/WTC11_EC_300gene_guide_quantifications_file_format_specification.pdf"
    },
    "lab": "/labs/jesse-engreitz/",
    "award": "/awards/HG011972/"
  },
  {
    "aliases": [
      "jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1"
    ],
    "document_type": "file format specification",
    "description": "File-format specification for per-cell quality-control reports from the WTC11-EC 300-gene CC-Perturb-seq D0-D3 analysis files.",
    "attachment": {
      "path": "/oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/documents/WTC11_EC_300gene_per_cell_quality_report_file_format_specification.pdf"
    },
    "lab": "/labs/jesse-engreitz/",
    "award": "/awards/HG011972/"
  },
  {
    "aliases": [
      "jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1"
    ],
    "document_type": "file format specification",
    "description": "File-format specification for day-specific annotated cell-by-gene-and-guide AnnData matrices from the WTC11-EC 300-gene CC-Perturb-seq analysis.",
    "attachment": {
      "path": "/oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/documents/WTC11_EC_300gene_annotated_cell_gene_guide_h5ad_file_format_specification.pdf"
    },
    "lab": "/labs/jesse-engreitz/",
    "award": "/awards/HG011972/"
  }
]

Review-only TSV:
/oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/document_file_format_specifications_review_only.tsv

```

## Markdown cell 14

## 7. Portal utilities and credential entry

## Code cell 15

```python
def split_link_items(value: object) -> list[object]:
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]


def terminal_identifier(value: object) -> str:
    if isinstance(value, dict):
        value = value.get("accession") or value.get("@id") or value.get("uuid") or ""
    text = str(value or "").strip().strip("/")
    return text.rsplit("/", 1)[-1] if text else ""


def direct_link_tokens(value: object) -> set[str]:
    tokens: set[str] = set()
    for item in split_link_items(value):
        if isinstance(item, dict):
            for key in ("accession", "uuid", "@id"):
                token = terminal_identifier(item.get(key, ""))
                if token:
                    tokens.add(token)
            for alias in split_link_items(item.get("aliases", [])):
                if str(alias).strip():
                    tokens.add(str(alias).strip())
        else:
            token = terminal_identifier(item)
            if token:
                tokens.add(token)
            if ":" in str(item) and not str(item).startswith("/"):
                tokens.add(str(item).strip())
    return tokens


def portal_urls(identifier: str) -> list[str]:
    token = str(identifier).strip().strip("/")
    base = f"https://api.data.igvf.org/{token}/"
    return [
        base + "?format=json&datastore=database",
        base + "?format=json&datastore=database&frame=object",
        base + "?format=json&datastore=database&frame=edit",
    ]


def portal_get(identifier: str, timeout: int = 45) -> dict | None:
    auth = None
    if os.environ.get("IGVF_API_KEY") and os.environ.get("IGVF_SECRET_KEY"):
        auth = (os.environ["IGVF_API_KEY"], os.environ["IGVF_SECRET_KEY"])
    errors = []
    for url in portal_urls(identifier):
        response = requests.get(
            url, auth=auth, headers={"Accept": "application/json"}, timeout=timeout
        )
        if response.status_code == 404:
            continue
        if response.ok:
            return response.json()
        errors.append(f"{response.status_code} {url}: {response.text[:300]}")
    if errors:
        raise RuntimeError("Portal GET failed:\n" + "\n".join(errors))
    return None


def get_with_retry(identifier: str, attempts: int = 6) -> dict | None:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            record = portal_get(identifier)
            if record is not None:
                return record
        except Exception as error:
            last_error = error
        if attempt < attempts:
            time.sleep(min(2 ** (attempt - 1), 10))
    if last_error is not None:
        raise last_error
    return None


_IDENTITY_CACHE: dict[str, set[str]] = {}


def resolved_link_tokens(value: object) -> set[str]:
    direct = direct_link_tokens(value)
    resolved = set(direct)
    for token in list(direct):
        if token in _IDENTITY_CACHE:
            resolved.update(_IDENTITY_CACHE[token])
            continue
        identities = {token}
        try:
            record = portal_get(token)
        except Exception:
            record = None
        if record:
            identities.update(direct_link_tokens(record))
        _IDENTITY_CACHE[token] = identities
        resolved.update(identities)
    return resolved


def links_equal(left: object, right: object) -> bool:
    left_tokens = resolved_link_tokens(left)
    right_tokens = resolved_link_tokens(right)
    return bool(left_tokens and right_tokens and left_tokens & right_tokens)


def link_set_contains(container: object, required: object) -> bool:
    existing_items = split_link_items(container)
    required_items = split_link_items(required)
    return all(
        any(links_equal(existing, wanted) for existing in existing_items)
        for wanted in required_items
    )


def portal_record_id(record: dict | None) -> str:
    record = record or {}
    return str(record.get("@id") or record.get("accession") or record.get("uuid") or "")


def portal_record_url(record: dict | None) -> str:
    record = record or {}
    if record.get("@id"):
        return "https://data.igvf.org" + str(record["@id"])
    return ""


def row_without_empty_values(row: dict) -> dict:
    return {key: value for key, value in row.items() if str(value).strip() != ""}


def enter_igvf_credentials() -> None:
    api_key = <REDACTED>"IGVF API key: ").strip()
    secret_key = getpass("IGVF secret key: ").strip()
    if not api_key or not secret_key:
        raise ValueError("Both IGVF credentials are required")
    os.environ["IGVF_API_KEY"] = api_key
    os.environ["IGVF_SECRET_KEY"] = secret_key
    del api_key, secret_key
    print("IGVF credentials are set for the current notebook kernel.")


def show_credential_status() -> None:
    for variable in ("IGVF_API_KEY", "IGVF_SECRET_KEY"):
        value = os.environ.get(variable, "")
        print(f"{variable}: {'SET' if value else 'MISSING'}, length={len(value)}")


def run_iu_register(
    *, profile: str, input_tsv: Path, dry_run: bool, patch: bool, label: str
) -> dict:
    if not IGVF_PYTHON.is_file():
        raise FileNotFoundError(IGVF_PYTHON)
    if not IU_REGISTER.is_file():
        raise FileNotFoundError(IU_REGISTER)
    if not input_tsv.is_file():
        raise FileNotFoundError(input_tsv)
    if not os.environ.get("IGVF_API_KEY") or not os.environ.get("IGVF_SECRET_KEY"):
        raise RuntimeError("IGVF credentials are not set")

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    mode = "dryrun" if dry_run else "submit"
    log_path = LOG_DIR / f"{re.sub(r'[^A-Za-z0-9._-]+', '_', label)}_{mode}_{timestamp}.log"
    command = [
        str(IGVF_PYTHON), str(IU_REGISTER), "-m", "prod", "-p", profile,
        "--tries", "3", "--delay", "10", "--backoff", "2",
    ]
    if dry_run:
        command.append("-d")
    if patch:
        command.append("--patch")
    command.extend(["-i", str(input_tsv)])

    print("\nCommand:", " ".join(command))
    print("Log:", log_path)
    completed = subprocess.run(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        env=os.environ.copy(), check=False,
    )
    log_path.write_text(completed.stdout, encoding="utf-8")
    print(completed.stdout)

    fatal = re.search(
        r"Traceback|ValidationFailure|Failed to POST|Failed to PATCH|\b(?:401|403|409|422) Client Error\b",
        completed.stdout,
        flags=re.IGNORECASE,
    )
    if completed.returncode != 0 or fatal:
        raise RuntimeError(f"igvf_utils failed; review {log_path}")
    return {"log_path": str(log_path), "success": True}

print("Portal utilities loaded.")
```

### Text output (stream)

```text
Portal utilities loaded.

```

## Code cell 16

```python
# Run this cell once after every kernel restart.
enter_igvf_credentials()
show_credential_status()
```

### Text output (stream)

```text
IGVF credentials are set for the current notebook kernel.
IGVF_API_KEY: <REDACTED> length=8
IGVF_SECRET_KEY: SET, length=16

```

## Markdown cell 17

## 8. Plan, dry-run, and submit the three Document records

## Code cell 18

```python
if (
    not os.environ.get("IGVF_API_KEY")
    or not os.environ.get("IGVF_SECRET_KEY")
):
    raise RuntimeError(
        "Run the credential-entry cell first"
    )

document_plan_rows = []

for role_key, spec in DOCUMENT_SPECS.items():
    try:
        record = get_with_retry(
            spec["document_alias"]
        )
        error = ""
    except Exception as exc:
        record = None
        error = f"{type(exc).__name__}: {exc}"

    if error:
        action = "BLOCK_LOOKUP_ERROR"
    elif record is None:
        action = "POST_NEW"
    elif (
        str(record.get("document_type", ""))
        != "file format specification"
    ):
        action = "BLOCK_EXISTING_CONFLICT"
    else:
        action = "SKIP_EXISTING"

    document_plan_rows.append(
        {
            "format_role": role_key,
            "document_alias": spec["document_alias"],
            "action": action,
            "portal_id": portal_record_id(record),
            "portal_url": portal_record_url(record),
            "lookup_error": error,
        }
    )

document_plan = pd.DataFrame(
    document_plan_rows
)
DOCUMENT_PLAN_PATH = (
    WORK_DIR / "document_submission_plan.tsv"
)
document_plan.to_csv(
    DOCUMENT_PLAN_PATH,
    sep="\t",
    index=False,
    lineterminator="\n",
)

print(document_plan.to_string(index=False))

if document_plan["action"].str.startswith(
    "BLOCK_"
).any():
    raise RuntimeError(
        f"Document plan is blocked. Review {DOCUMENT_PLAN_PATH}"
    )

metadata_by_alias = {
    payload["aliases"][0]: payload
    for payload in document_payloads
}

# Delete unsafe v1 operation TSVs so they cannot be reused.
removed_legacy_document_tsvs = []
for legacy_path in ROW_DIR.glob(
    "POST_document_*.tsv"
):
    legacy_path.unlink()
    removed_legacy_document_tsvs.append(
        str(legacy_path)
    )

if removed_legacy_document_tsvs:
    print("\nRemoved legacy Document POST TSVs:")
    for legacy_path in (
        removed_legacy_document_tsvs
    ):
        print(" ", legacy_path)

document_operation_rows = []

for plan in document_plan.to_dict(
    "records"
):
    role_key = plan["format_role"]
    alias = plan["document_alias"]

    if plan["action"] == "POST_NEW":
        json_path = (
            ROW_DIR
            / f"POST_document_{role_key}.json"
        )

        payload = metadata_by_alias[alias]
        json_path.write_text(
            json.dumps(
                payload,
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )

        # Validate the exact file that iu_register.py will consume.
        reloaded = json.loads(
            json_path.read_text(
                encoding="utf-8"
            )
        )
        if reloaded != payload:
            raise RuntimeError(
                f"JSON round-trip failed: {json_path}"
            )
        if not isinstance(
            reloaded.get("attachment"),
            dict,
        ):
            raise RuntimeError(
                f"attachment is not an object: {json_path}"
            )
        if set(
            reloaded["attachment"]
        ) != {"path"}:
            raise RuntimeError(
                f"attachment must contain only `path`: {json_path}"
            )
    else:
        json_path = None

    document_operation_rows.append(
        {
            "format_role": role_key,
            "document_alias": alias,
            "action": plan["action"],
            "input_file": (
                str(json_path)
                if json_path is not None
                else ""
            ),
            "input_format": (
                "json"
                if json_path is not None
                else ""
            ),
        }
    )

document_operations = pd.DataFrame(
    document_operation_rows
)
DOCUMENT_OPERATIONS_PATH = (
    WORK_DIR
    / "document_submission_operations.tsv"
)
document_operations.to_csv(
    DOCUMENT_OPERATIONS_PATH,
    sep="\t",
    index=False,
    lineterminator="\n",
)

print("\nOperations:")
print(
    document_operations.to_string(
        index=False
    )
)
```

### Text output (stream)

```text
            format_role                                                                         document_alias   action portal_id portal_url lookup_error
  guide_quantifications          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1 POST_NEW                                  
per_cell_quality_report        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1 POST_NEW                                  
         annotated_h5ad jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1 POST_NEW                                  

Operations:
            format_role                                                                         document_alias   action                                                                                                                                                                     input_file input_format
  guide_quantifications          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1 POST_NEW   /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_guide_quantifications.json         json
per_cell_quality_report        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1 POST_NEW /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_per_cell_quality_report.json         json
         annotated_h5ad jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1 POST_NEW          /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_annotated_h5ad.json         json

```

## Markdown cell 19

### 8A. Validate the exact Document JSON operation files

This cell verifies that every pending Document POST is valid JSON,
that `attachment` is a native object, and that the referenced PDF
exists.

## Code cell 20

```python
document_json_preflight_rows = []

for operation in (
    document_operations.to_dict(
        "records"
    )
):
    if (
        operation["action"]
        == "SKIP_EXISTING"
    ):
        document_json_preflight_rows.append(
            {
                "format_role": (
                    operation[
                        "format_role"
                    ]
                ),
                "action": (
                    operation["action"]
                ),
                "input_file": "",
                "valid_json": True,
                "attachment_is_object": True,
                "attachment_exists": True,
                "status": "SKIP_EXISTING",
            }
        )
        continue

    input_path = Path(
        operation["input_file"]
    )
    payload = json.loads(
        input_path.read_text(
            encoding="utf-8"
        )
    )
    attachment = payload.get(
        "attachment"
    )
    attachment_is_object = (
        isinstance(attachment, dict)
        and set(attachment) == {"path"}
    )
    attachment_path = (
        Path(attachment["path"])
        if attachment_is_object
        else None
    )
    attachment_exists = bool(
        attachment_path
        and attachment_path.is_file()
    )

    document_json_preflight_rows.append(
        {
            "format_role": (
                operation["format_role"]
            ),
            "action": operation["action"],
            "input_file": str(input_path),
            "valid_json": True,
            "attachment_is_object": (
                attachment_is_object
            ),
            "attachment_exists": (
                attachment_exists
            ),
            "status": (
                "PASS"
                if attachment_is_object
                and attachment_exists
                else "FAIL"
            ),
        }
    )

document_json_preflight = (
    pd.DataFrame(
        document_json_preflight_rows
    )
)
DOCUMENT_JSON_PREFLIGHT_PATH = (
    WORK_DIR
    / "document_json_preflight.tsv"
)
document_json_preflight.to_csv(
    DOCUMENT_JSON_PREFLIGHT_PATH,
    sep="\t",
    index=False,
    lineterminator="\n",
)

print(
    document_json_preflight.to_string(
        index=False
    )
)
print(
    "\nSaved:",
    DOCUMENT_JSON_PREFLIGHT_PATH,
)

if not document_json_preflight[
    "status"
].isin(
    {"PASS", "SKIP_EXISTING"}
).all():
    raise RuntimeError(
        "At least one Document JSON "
        "operation failed preflight"
    )
```

### Text output (stream)

```text
            format_role   action                                                                                                                                                                     input_file  valid_json  attachment_is_object  attachment_exists status
  guide_quantifications POST_NEW   /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_guide_quantifications.json        True                  True               True   PASS
per_cell_quality_report POST_NEW /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_per_cell_quality_report.json        True                  True               True   PASS
         annotated_h5ad POST_NEW          /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_annotated_h5ad.json        True                  True               True   PASS

Saved: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/document_json_preflight.tsv

```

## Code cell 21

```python
RUN_DOCUMENT_DRY_RUNS = True
DOCUMENT_DRY_RUN_MARKER = WORK_DIR / "DOCUMENT_DRY_RUN_PASSED.json"

if not RUN_DOCUMENT_DRY_RUNS:
    print("Set RUN_DOCUMENT_DRY_RUNS=True and rerun this cell.")
else:
    results = []
    for operation in document_operations.to_dict("records"):
        if operation["action"] == "SKIP_EXISTING":
            results.append({**operation, "success": True, "log_path": ""})
            continue
        result = run_iu_register(
            profile="document",
            input_tsv=Path(operation["input_file"]),
            dry_run=True,
            patch=False,
            label=f"document_{operation['format_role']}",
        )
        results.append({**operation, **result})

    result_frame = pd.DataFrame(results)
    DOCUMENT_DRY_RUN_RESULTS_PATH = WORK_DIR / "document_dry_run_results.tsv"
    result_frame.to_csv(
        DOCUMENT_DRY_RUN_RESULTS_PATH, sep="\t", index=False, lineterminator="\n"
    )
    if not result_frame["success"].all():
        raise RuntimeError("At least one Document dry run failed")
    marker = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "document_metadata_md5": md5sum(DOCUMENT_METADATA_PATH),
        "pdf_md5s": dict(zip(pdf_inventory["format_role"], pdf_inventory["md5sum"])),
    }
    DOCUMENT_DRY_RUN_MARKER.write_text(json.dumps(marker, indent=2) + "\n")
    print("All Document dry runs passed:", DOCUMENT_DRY_RUN_MARKER)
```

### Text output (stream)

```text

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -d -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_guide_quantifications.json
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/document_guide_quantifications_dryrun_20260930_092553.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:25:53,757:iu_debug:	DRY RUN is enabled.
2026-09-30 09:25:53,758:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:25:53,758:iu_debug:	submission=True: In submission mode.
2026-09-30 09:25:55,122:iu_debug:	
IN post().
2026-09-30 09:25:55,174:iu_debug:	Validating the payload against the schema
2026-09-30 09:25:55,192:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_300gene_guide_quantifications_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOCAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDcgMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL1BhZ2VNb2RlIC9Vc2VOb25lIC9QYWdlcyA3IDAgUiAvVHlwZSAvQ2F0YWxvZwo+PgplbmRvYmoKNiAwIG9iago8PAovQXV0aG9yIChFbmdyZWl0eiBMYWIpIC9DcmVhdGlvbkRhdGUgKEQ6MjAyNjA5MzAwOTI0MTItMDcnMDAnKSAvQ3JlYXRvciAoXCh1bnNwZWNpZmllZFwpKSAvS2V5d29yZHMgKCkgL01vZERhdGUgKEQ6MjAyNjA5MzAwOTI0MTItMDcnMDAnKSAvUHJvZHVjZXIgKFJlcG9ydExhYiBQREYgTGlicmFyeSAtIFwob3BlbnNvdXJjZVwpKSAKICAvU3ViamVjdCAoSUdWRiBmaWxlIGZvcm1hdCBzcGVjaWZpY2F0aW9uKSAvVGl0bGUgKEd1aWRlIHF1YW50aWZpY2F0aW9ucykgL1RyYXBwZWQgL0ZhbHNlCj4+CmVuZG9iago3IDAgb2JqCjw8Ci9Db3VudCAxIC9LaWRzIFsgNCAwIFIgXSAvVHlwZSAvUGFnZXMKPj4KZW5kb2JqCjggMCBvYmoKPDwKL0ZpbHRlciBbIC9BU0NJSTg1RGVjb2RlIC9GbGF0ZURlY29kZSBdIC9MZW5ndGggMjI3NQo+PgpzdHJlYW0KR2F1MEU9XT0/PyY6V2VEb0hRLlxOOS9lX20vdE8zOVsqUD5SQC9Kck1BITBiNERXLVYjPkQpcG5fN2UjXFZMXXQhIVNLMnEsO00nO2c1Ols1Zls1JHMkbmZaaGE/JFxBXFxdazRvaG9HJFxpJElBQFY2QFVOdCgtRjpVNFNaamtwMDpZWUlMJGEuZyFFN1hXQCJPay9AQDs7cUkwdD1acTNScSUyMmBoYllYKVBgVjxqTzhQYkFxYVsjTDFpKUswX1AhMDxlJi4hYlM9bEdWRyhicEkoJD8mYHRiOVEvQklaY0hnczJfRiMuSEdZJDhNa0UmRVcuIXU6ajQ6YCwkJCY2MydwZlQyV0pFV0ZTMCd0OUE/TGI5I1tQN1M/bEopWjwsZG85ODlyJz5WJHRfUkwvOT5wYjZWRkBKY1dMaWAjPjUqZFQsUEwtVjZzX1tWJShKckNqXGcwQDAhb1A3MlE7VC9AYG5DTTZnQkFfOXVSJEJqRS4pWGI9UFNAQ1ZEMmRMVWZkbVszXjo6aEQtaiIuc2IoQTdELihSQkMrJjZEN05lWzEsJjZGcTowNFJVSDtjVlNkb1FDTyFlTWk4Yi87PnAjI1onUUEuTlwuV0xeVV4vPWgpbSxFLWdSc2pNWE91N1YoIz1zInB1WCVqRj4sUGAwZkhGOjItWyZnYjs/XkgsdVhoTCMqPVp1akFHXVhiKUJtPCRyP0sxQFYxcjAmbmc3bVw+RDs2IihqaVZRTVFPb2xGLltVb3NKZWE4JF8nXUJ0Z0w3XjdkYVZTUGpDWTotNTFdXUdKdWJnKEUjTCVRXiNQWm5oRlpJRWlXISMwaCVaWi8qUSMlSHU2azBXZD4/JVhKOCNyK1RaRDtWW05vVDVFP2ZCKScnVURgNStVJ3JpRUAzWlooRWU5RFNfWTBXIylkXTE1TVM1KDpzXks1JU0lblA5RjcpKTBUQGxsaEQtZl5iZz9fP1tOcTxWMzFdRFVqODdLKEdYX24qS0BxLFRfbF0nRmYuXGMxLVoyNGhPUiFsNjpJXFdBaTRaYDMuYl1BRyp1MlA5Y0UyQ0BhcmNJSEhrK2Fpb09ZPkQ3S0RsRmhyZHQlYjdudUg2QjtrQW1JRXEoTjRrOmlNbksiYzdMR1RdbVRQTlAoNzZuWm5JUFtUbFEiITU0MTtgJmVfOWdxUD5wMydOXFpSOVZmRXIwUjNdVURPREdgOCx0Ty09Z2k9O2NmakE2USVhUitrLHJsZjhvaik8TXJfbDNEY3BRQjZyTiZMal8qVTFkQURKWVY0N1c+KCMiQnF0ajxKc3BXcCMoJlYhWV1dZ0xfRnQ6MT1HVi5OYnMmdWJsW3JbXFpJKShZXVNlKGFuOnFZbEJUIVErNSRPJl5XKVApZGE3NSwvIWc0N3NKZGBQRDslXyI8NVQnXEgkQ1BQcmA3I1VaY1A0QjctTVcxQl5OVk5ETUElSTMoYDpzaD0yZCdpalwvOy4uJ1Vaaj1NVWZWaSQ8IzVHTnFoTHFLdTEocjYhZF5dW3JeaERgXltzbTUyS25rMz84aGYiRHBDbDxcO1hda2Y7T1gmaWVjNlVJYHQ4cWNaOkFLJGswPWwhOHJNJGAhXTEoWmcqYkc/NjZSJVVVWStWcjQvZlBlbCEtdTYmLDoyQVJmMT5zZzE2UTc/dSVQJjRMUHJbPXBpZEpZNlhMXC5TZ1JAIjVwa1pAa29AXlFvNDxGSG9MTWJHZWNyazlzPmA2V0g2dUhpR1tyb3I5SEQ1bHNkSWtaUEFvUVlaT0JNT3RPMlYyKlxddFZjP2NaZ0ciUi0ucEtxJlgtb1lcdUZbKzNVQ09hYlpgUTM9VWk3bUNDTEY+OVAnXi1qXUc7PzRNYjo1TCh0SV86O24tPjRZcW1eIyEzVGVNKnE7SHQxZi5BZlZhNmpAQ04uPWNGZlBGT0QyTTlRLVtWXVYkOyk9O2g7Iz9rYXUqSkVfYHJbLE1USzczcF1xPmtkUlhyaE9hTFVrLFk+YnArMFdVXjE6UjsuPko/SHArdEhJN0ReNkdlVk05dUIxO09ILGJNJV4kPyk+amYzTypZZzokXUAsKyI9WDFsTmoxJENhPSx1XkRyN2hPZlxOWzU2ZVAwXDkxTHRpNEY9VmRvciQhdVpAQ2NeJ2NTUzMuaXVaVm9WQ0A5dSpxUG86TE9GKlUjdXBeXGpgcTQ/XGojUmddYUhbJ1JkbnUwV1hQTjJmMkkoXj9JLy9AXTFNUm0wV2BYI2Q5bmNaZDYtXEk+T3FtMWxMbCNUJUVsUlVSLzNsR0kpc0AnPDZBUSpYLmBVdSZ1T01nMTQ4LElfOjdCak8yP0UiI2BvZz0nZGQsPFcsRjE6JFcoXTAwYTRUQGtAbW9JWUEiP1VcJy90JEdvK1ZhLmAzXGFRaVk3J1NqUVE2Qi1OLGxUZSJdVyRaKEU5TXBpbz0qSWcuK1QwKkZtb1krV0M6YV1qXSdmWzw+RWlSbG4zWWk2OyMpUzpOPEphIUVxZjs0J2RkKkhNZlomKyw9UzlDI0RIWG4lJTdeWm00OnIyVWhgN0BnITVfdThXUVNBcV8vW01hPV4+TCtwMytBLFheL20oKlJ1OXM4Pi4ibDc9WUglRCtTaUw/NmhzIyZnXjBob0BJPTdhMShaXz5sZV04cC1PSTxwIks6US5kPiFkbUhBYlY1WzdXSigvTW5LUTA8OE1LZVgxPC5QW1QyRkBzPXUybythPS5QRl09VydhYk5EbF9aJEMjcTlfOygySitqMzhhWz5oN2BaKTglRWRkMl0pQm1mKGMsQnVOSVFqPlwibFp1JmUvJSckM2U8Qmw1UTJMaUtrQU1AYWlWc1o3YnEuTitsZU9ZUzknPDc4QmJQXDgxYzteSi42WEVPOVVfbGx1cDpGNXJdLVouQj5SJU9DRlM8IzFiQ0hMQTNyTzcxaWZvM3BlYGZwTCpvUzg5Ll5uMmwmWy8za3BnNWE5PjVHZ2s/SSxSUT5aK1psSk84Ql8qNk5nbmNgPVhuY011W286PDx+PmVuZHN0cmVhbQplbmRvYmoKeHJlZgowIDkKMDAwMDAwMDAwMCA2NTUzNSBmIAowMDAwMDAwMDYxIDAwMDAwIG4gCjAwMDAwMDAxMDIgMDAwMDAgbiAKMDAwMDAwMDIwOSAwMDAwMCBuIAowMDAwMDAwMzIxIDAwMDAwIG4gCjAwMDAwMDA1MTQgMDAwMDAgbiAKMDAwMDAwMDU4MiAwMDAwMCBuIAowMDAwMDAwODg0IDAwMDAwIG4gCjAwMDAwMDA5NDMgMDAwMDAgbiAKdHJhaWxlcgo8PAovSUQgCls8YTE2YmI2MTYyZjc4NWNhNmYyZDYzYTI5NWU2MzQ1YjM+PGExNmJiNjE2MmY3ODVjYTZmMmQ2M2EyOTVlNjM0NWIzPl0KJSBSZXBvcnRMYWIgZ2VuZXJhdGVkIFBERiBkb2N1bWVudCAtLSBkaWdlc3QgKG9wZW5zb3VyY2UpCgovSW5mbyA2IDAgUgovUm9vdCA1IDAgUgovU2l6ZSA5Cj4+CnN0YXJ0eHJlZgozMzA5CiUlRU9GCg==",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for per-cell guide assignment quantifications from the WTC11-EC 300-gene CC-Perturb-seq D0-D3 analysis files.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 09:25:55,192:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -d -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_per_cell_quality_report.json
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/document_per_cell_quality_report_dryrun_20260930_092555.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:25:55,703:iu_debug:	DRY RUN is enabled.
2026-09-30 09:25:55,704:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:25:55,704:iu_debug:	submission=True: In submission mode.
2026-09-30 09:25:57,004:iu_debug:	
IN post().
2026-09-30 09:25:57,031:iu_debug:	Validating the payload against the schema
2026-09-30 09:25:57,049:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_300gene_per_cell_quality_report_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOCAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDcgMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL1BhZ2VNb2RlIC9Vc2VOb25lIC9QYWdlcyA3IDAgUiAvVHlwZSAvQ2F0YWxvZwo+PgplbmRvYmoKNiAwIG9iago8PAovQXV0aG9yIChFbmdyZWl0eiBMYWIpIC9DcmVhdGlvbkRhdGUgKEQ6MjAyNjA5MzAwOTI0MTMtMDcnMDAnKSAvQ3JlYXRvciAoXCh1bnNwZWNpZmllZFwpKSAvS2V5d29yZHMgKCkgL01vZERhdGUgKEQ6MjAyNjA5MzAwOTI0MTMtMDcnMDAnKSAvUHJvZHVjZXIgKFJlcG9ydExhYiBQREYgTGlicmFyeSAtIFwob3BlbnNvdXJjZVwpKSAKICAvU3ViamVjdCAoSUdWRiBmaWxlIGZvcm1hdCBzcGVjaWZpY2F0aW9uKSAvVGl0bGUgKFBlci1jZWxsIHF1YWxpdHkgcmVwb3J0KSAvVHJhcHBlZCAvRmFsc2UKPj4KZW5kb2JqCjcgMCBvYmoKPDwKL0NvdW50IDEgL0tpZHMgWyA0IDAgUiBdIC9UeXBlIC9QYWdlcwo+PgplbmRvYmoKOCAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAyNzE2Cj4+CnN0cmVhbQpHYXUwRj5CQVEtJnE5U1lmWEVtXk1dbT8uQUJtUE9hLWgldDFgKiFaZSsmQk8rc0RlKyxVZTswPyRaOyxdKmplNU1qLzNsSFllZmVvak4oN0BCbUVFIXUkMFBxQitsI210VT8uUTMiPSIzISYuWSVMYlM8cXIsTm1LYjpROTFCVDw0Z1UwP1E1TVtIZCo7SzZkIW9KPTomN0FeZUBHLUQoNW46VUopSUBLT189cmlsZkksSF4iYCQoLFc4bVpNKWhSJTtmWUJhSWttT1IpOichMCpuLGFIaEU7PGc2LHU6UWhwXzhdciEkPUgnSUxEPnJqIk1vbSpgWiVeSCg3WFgzMyhDbl1SbypMMGA3PiJnSGJDSisxcTUqNWVmXWMlVygxbVlHaCViITBZdVchXTFLbEptJD1qOmZ1a0pJWG1gV2dWRFBVZychMV0/ZiwmJiFuLGYuODg2N2clRl81cV9aWT5KRVkuZ25LRFMqJVY6PzIoVCw0ZUpgPkVGKVhiPkNFYk9RVVdFKVtadEFfPzhRNUIyJ2VHPkJkJSg1a0pXYGBQZ08hbUBmcmEpL2FdX1tdaiMqXU5JJDRFK0hfIllkaUlKKUhDSDIuTTo1ZjNrKjFAKi9ZLmxnIlJkPT1rZjhuYlxTUGdWIlYpYkh1R29ENDxwUVVeRFFsUVFbXFhMOl00TVxqazRUO1hLVmJuMyZEZCokVU0xYV5wP3FJYk5vJ1tLS1c9IlExdUc6QkgpYyw2PDFdLUheaiEzLipGVSkuRGE8SDElL0tjXSZNI0sjOzJSSlNCbFRJQzpAaGhmLThINylzVStVKHFicFpAKmBlQVcjM2FhczM5aFtPP1gwKUQqZVNoY09Vc1xodT1nPU5wc1pJMClELzIpJkFfTmYncCg3dSlYdVRmMl0/LD5NWCVDNUBwZzgpXEIkXVNaKmNvPlJiRTxUQ2liY1NJczJlKGhfRjtZRE40ZipUIXJnUS1Na10jL2BhNSRwVShEck88O0lmSlQ4XWs+RitsaTJHL1VWLjI+Xmo4Kj9vNHQ4VDM7WU45SyIpQCwoRFRvbTNVXWdMWUVzT0Q9JDNGQmMjZTouMTBeJihNVjNpOFRpT2RoV3VvNzY1bF9YRjxeKVJUPDpMKUByY1xPLEJmMU1FXiNiMkMkJ3VjbWEiJ10qYEZZV3VbJmhVWDpTNiZxUlFbOTk/N0xfV0dBMSU8Sy8lZWNIWmMxX1JPSUJGTmwxRSs/XSZtYEFYM1swKVAnbz1eNWFdI01GVWQpb188c2RHMSRLbSU1IjdYSXBVbWhOT2coKiQpQW0/SEdBQ28zblpSLWhPJnA3KXAzRUZvO0hVJkZtdC9CPiFqKT1QJV9EPlslXkdFTEMwWj87OWtrc2E4cikrPlwmYXRQXStWR0dkQHJdI25OLydfS1xXaExKNzFGcnIkdEpwcnNGOiRcM3IzPG9ESDRWO2E9OWBHJykzYXI/RCZRUD8idFYqV2pMRFJHbnVRZUtJYl4/cisuclVvbmpEcTQtIVgmYEtBZy9iayY+QixwQz5PZiI5QVNZRk9nRGRJXSxtOFdTLCgpTmY2b18lTmclTy06TjdROGY+Wl80LltQX01hVVgtVVpfYS0uY3BRUTVBZUZbU28/VFBjSzFNZCNNU3RJMjdIOWJKQG8vb1xdbHVaW0kjO2tcXUY6UDdYS2tWUShVXiopKDshaEJ0JXFEYjs8Mi1tQ2dqP1MkVy9IJipmUVorbjNAPCUjOiRjPD82PnI4Y2hdQ1xNcj9CNSZFdWEhbiVFSHFPcyMlZFc7RzwuSExTOTprZFI6ay8nM3EyV2olRiEhYD0mYyYtQDJXOD4yJSIxNG4sYj9TXnNzM1U1LUdLI2RJVHNrVFBrWVFRTE9oNT8ramIjO14oR1RjI0hqKHQyazYqUkUpLTpsTlNbXmJybFs8Ml1kPi9nTFhHbzBPX3BHVWxbSktJQTtiMVg9Yyg7MVYpWSdHUDNhP1dLLE1vS1o5OSprcG10U2pZTSNrUylVRyNNYCo1NGt1TksxO1BAQTNNdGRHK3FKVykoQztpRC5OcE4kVi1zKDI/IiZCbU0+W1lSIjM1KFEnbUQoIk1jKkRdKDZVNm1gQ01mREFOQ1FgRGIjRShyVCshKCFvKU8pZGZVVmc3OTMpMlRjN0RvNGpRMTs7b0xpMz87QkxOKnFUY09kZl0lJFNLQDQyUSxUSCozNy8sX0ZmTEYjXVFJR3NhPmY+PDc2WEZWT1hrQUU8ZiE2RE1UPzYpPjI4O2ZvSTs3YTJwTHRbLyg/IWVlT2BtU190NiFXQTQ6VCZdPDJJO0xPOTA+b1JRZkZWTUA1WT49MFxdOzh0QV49KjYsKGpXYCMiKWtKMi5lQyRwSC8zOkFKLk9aK04iYSxIKGg4SWFCOzdZTiRkSzo0Ly5bdWQibmFcTVtsMUwzc01HSD1eT0lCTyZmSS0vXWRrbmM5ZF9sJHJhTiMtXidiYUVARGxSZztQMDQpcWUjVUAvKS1FYjsyIlxAXlY/XS9DWCtqclVdVkM0WihvR3AmLksyc05xI2wwMiFWYzFZWmJqcHROKyFVay9nQT4+JyliZXRcWzMhKzc2cWBUS2FZdU5XZWI7VCN1RmQ5VldMTk5FckRSVDtbW11ZXmhmbmdjKytnVlknbWpYPFplNzgwS25jJDtVJCFicnExNl41JjxaXUg0NzUlZipiUiZLLSdjX2MoPDwsOkpZMHJbSkRrQjRZSVBkXkpGcDZba0ZHbC5xTmVzbmEiNVVraiYzJkVNKnFmaFVZcnJBcF91bTI9S3JvWEMjMkZFIklWWWhzMU1gXXFHdUFsaDU7MlVYOEBWVzk0JGYlK1cjKks7ImIuW2tMZVFyWjZmL2hWWjhNLSNDMW1hJyxGanFMJElZLEEpOWxUSDtUYSljKz5KNGtza10qNyJvViU/XWRuaiZxME9aKm1uWWZmaEwtQDpJcCdCLSdoJyhwT1BGNGopM0s6Qic9S0NRI2puUlBmRkkzdXQ5WzNLK091WXI9QWRkQTMpaUUtXU0lSDU0JkZxUDJrIWpGOVs+UHFEPEopdUBGQyRDZSZyWFwzXUFXRVI3MTY7Z2YvbWBKNnFHU1olcyZZQVwmV0hRYSwtPmNncjg6dW9FTjJrXlxiQks7Im9wcl1ySTFjZTlHcE02P0w5RnU8Sz1RVDlFOlhFW1U0W2RoSElxQ0ssV1VaTGNEQE9oZGYtWmxqRTM+cDpHTDBZTygqSVQsaGRjdEwpZlQjWiNgYC5YQidSSzRJM2VJP2okZm1VPiFgKmlwbFgqYEJpTzZjODUsNlVIWG0uVjtjdEFwXCZdRVNuWz43P102XlE3VDg0bklDWm5iMWJCP0orPylhJVRWRSpBRWdKViwsPEVSQHMjR1BCWGlGaS5SM3A7OFdmI2YoKVQuaTVRaVYvQE43T25jREZtT0RFTDhJSFIpcmVsU0xNclZbSE9BO2RDWVBGUDIvQE8qJTxELlM8NkxuLV8/OFRVbk9jTldgKTU0OVUpNz4uV0I/X05haiwjYyUtLHQ6WF9HXDNjKGIvP0ZzUk1hPF1sODNyJ0A2R3EnJFojJ1goIk1na2FNMltFMX4+ZW5kc3RyZWFtCmVuZG9iagp4cmVmCjAgOQowMDAwMDAwMDAwIDY1NTM1IGYgCjAwMDAwMDAwNjEgMDAwMDAgbiAKMDAwMDAwMDEwMiAwMDAwMCBuIAowMDAwMDAwMjA5IDAwMDAwIG4gCjAwMDAwMDAzMjEgMDAwMDAgbiAKMDAwMDAwMDUxNCAwMDAwMCBuIAowMDAwMDAwNTgyIDAwMDAwIG4gCjAwMDAwMDA4ODYgMDAwMDAgbiAKMDAwMDAwMDk0NSAwMDAwMCBuIAp0cmFpbGVyCjw8Ci9JRCAKWzwyM2JhOTEyNDUyZjgzNDNhMmU0YjJlOWVhN2NiMzRjND48MjNiYTkxMjQ1MmY4MzQzYTJlNGIyZTllYTdjYjM0YzQ+XQolIFJlcG9ydExhYiBnZW5lcmF0ZWQgUERGIGRvY3VtZW50IC0tIGRpZ2VzdCAob3BlbnNvdXJjZSkKCi9JbmZvIDYgMCBSCi9Sb290IDUgMCBSCi9TaXplIDkKPj4Kc3RhcnR4cmVmCjM3NTIKJSVFT0YK",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for per-cell quality-control reports from the WTC11-EC 300-gene CC-Perturb-seq D0-D3 analysis files.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 09:25:57,049:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -d -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_annotated_h5ad.json
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/document_annotated_h5ad_dryrun_20260930_092557.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:25:57,563:iu_debug:	DRY RUN is enabled.
2026-09-30 09:25:57,563:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:25:57,563:iu_debug:	submission=True: In submission mode.
2026-09-30 09:25:59,058:iu_debug:	
IN post().
2026-09-30 09:25:59,085:iu_debug:	Validating the payload against the schema
2026-09-30 09:25:59,103:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_300gene_annotated_cell_gene_guide_h5ad_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgMTAgMCBSIC9NZWRpYUJveCBbIDAgMCA3OTIgNjEyIF0gL1BhcmVudCA5IDAgUiAvUmVzb3VyY2VzIDw8Ci9Gb250IDEgMCBSIC9Qcm9jU2V0IFsgL1BERiAvVGV4dCAvSW1hZ2VCIC9JbWFnZUMgL0ltYWdlSSBdCj4+IC9Sb3RhdGUgMCAvVHJhbnMgPDwKCj4+IAogIC9UeXBlIC9QYWdlCj4+CmVuZG9iago1IDAgb2JqCjw8Ci9Db250ZW50cyAxMSAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDkgMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjYgMCBvYmoKPDwKL0NvbnRlbnRzIDEyIDAgUiAvTWVkaWFCb3ggWyAwIDAgNzkyIDYxMiBdIC9QYXJlbnQgOSAwIFIgL1Jlc291cmNlcyA8PAovRm9udCAxIDAgUiAvUHJvY1NldCBbIC9QREYgL1RleHQgL0ltYWdlQiAvSW1hZ2VDIC9JbWFnZUkgXQo+PiAvUm90YXRlIDAgL1RyYW5zIDw8Cgo+PiAKICAvVHlwZSAvUGFnZQo+PgplbmRvYmoKNyAwIG9iago8PAovUGFnZU1vZGUgL1VzZU5vbmUgL1BhZ2VzIDkgMCBSIC9UeXBlIC9DYXRhbG9nCj4+CmVuZG9iago4IDAgb2JqCjw8Ci9BdXRob3IgKEVuZ3JlaXR6IExhYikgL0NyZWF0aW9uRGF0ZSAoRDoyMDI2MDkzMDA5MjQxMy0wNycwMCcpIC9DcmVhdG9yIChcKHVuc3BlY2lmaWVkXCkpIC9LZXl3b3JkcyAoKSAvTW9kRGF0ZSAoRDoyMDI2MDkzMDA5MjQxMy0wNycwMCcpIC9Qcm9kdWNlciAoUmVwb3J0TGFiIFBERiBMaWJyYXJ5IC0gXChvcGVuc291cmNlXCkpIAogIC9TdWJqZWN0IChJR1ZGIGZpbGUgZm9ybWF0IHNwZWNpZmljYXRpb24pIC9UaXRsZSAoQW5ub3RhdGVkIGNlbGwtYnktZ2VuZS1hbmQtZ3VpZGUgQW5uRGF0YSBtYXRyaXgpIC9UcmFwcGVkIC9GYWxzZQo+PgplbmRvYmoKOSAwIG9iago8PAovQ291bnQgMyAvS2lkcyBbIDQgMCBSIDUgMCBSIDYgMCBSIF0gL1R5cGUgL1BhZ2VzCj4+CmVuZG9iagoxMCAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAyOTMxCj4+CnN0cmVhbQpHYXUwRj5CQVEtJnE5U1lmWElgS1FtVFdHb0k3MVhcIzhJSVhZbWU5MWJfQT9iYFBHTVV1Sj8lKzVbNztNYVdzOV1oWGBxIlUhOGBRWicuRSFiPl4hSiM3XDgpaDtqPUUsai1RRUw4Pi9cWCFXMDJ0Mi06cEFXOEM9bGU8XD1KPWRGYU1ncS4lKFFhMWhqaVkoJXFjYj4qI21cPCFaczwzR2BoKz5KJCVdJ1g5Sy5bRFYyZzBEQFBtQkphR202R3JCaCItcyNtSFhCJidzMVdUJlJRdTFkI2xZS1FIQyVXIkwmVS8xKGhxZTNOJUVpamFeYlJRJFU0UkVUYDlCJTo3QCEjZ1FTPGFrWmBnMS8iMnF0XTppcFhWcUdzWzE/akJMLUtKcDpIODwjYms2bC9FajBbO2QmUSJAWk1LXWlhPEJCXSlTW19cNkUoWTlYLy4nO1c6KVUyTy1aSSgkWGJTS14tPEwyYTNXMC4+JmFhL0FhW20rNEc9YFEuI2BWRjg4L2VZcilZIWFROW00LVBGTkknYlgtcWAva0ZhNmZhTWVrK0lfL2UkaXInXUJdbW1MR09YMzJWYGdoVWJaY3NlUXFDM0BxX2tEYV1FaEBrZjs/MmJJIyxYbVAhQCRWOSVfUVU6PWBtLnB1Sm4qUVg1Lko2ZkoiJS8mWSFCTVhJOnEpbVppUkxMaElGQyxBS1lGZiJnNkk8KiNsRUxbNT1nP1JjWzJqPzgvR1FESTBENFFOIz5EUXQsWUtacEpeNUFlTGBYVm4lakMtRHFBUCJVMFlDPSw4KyQyXDJSMnFbWmRAN00ycTw3QTkqOj5eZVMhYWZlK2szJF51UScraz1uKjc2aSpfTFVqQT4nNU0+Ok1oZC5aO1tmPzc+TjlvLHFIUzg9SSNvP0ozZ21lWkU0Q0pvZio5aEIhbnRwYSpdMXNQNyE+WjdEP19IaCUpc1dpdDM9MUZzJzBGaWB1bnFsWmJLSnFHO15vYlpcP1lbaGVJSl0zUVpaU1A0cC4pLCQvMjFhdDNPP0BOYW9UKVxFVmhwLFBbcm5AbEtrUChIKUlsKVklYmtKaV5ddCV0WU5wUC0tclVsNXA7N2pfclU8biVsJDViPWtpTiY0O2Y4RlJOWSEiMXMrcjBFKUJnQkI4T2BnS1xtVGdILElzaj1sQGMxKDIlV1MoPzhRSi1iaF4vT3BCaTttc1xjZDUoOmpGIzMtLih1K25Uc2krSUNQZFwsbGArRU5yPldBLlInYCxhSkIzUjcmLm8uT2YhL2pNK3BtIl5zJystZGVDQEEpaVErVSghXVdcYS9pLF8/JnAzMChQL2g5LWcoVT9CYkc8Rlk0PyVpI1tDbVIyImNKTDovLjlxVipkaWNEaThdXk1TMUg+WyRgYF9eSGwkaUZlOVhqVC0ybSpFLHJzLT8mT0hDYXBYLz4kdEssIltSZGhxKltiTmVBaC8sSU11UEhecTFBaSwqMzBGL1ByYVAqIkM+V0VlUD5aamZhKVwlWElaamxRPitwRDJFNUg3YURtXzlPQ2xcM0JMU10hb21VVjthL0VLS05PdGFEUDE1SltUW0xpVTZrU0hwOVYpMlRKb0xSJnJxbytFSz1PNltYdGpfPihPcz85LiYiOWk+UCxgTysrZlJLVz8vckJnZ3QnZyo+T1kxR2UnJD08TUc+b01hXGcxVlVOY10jYy9BVy0yK1BYNz9wSV1ANFtkNFBUby5xTSFvSi8sKlRIK1pQKWhGXXVPSzo+VnJWTVFaMHE1TFZINytpUUFGYCM5Qy4vJ1lQcj5kQTxBR0hQbypXWjJZYGJwYEA9TDc5TWRlVSFEND03JnNjNWxWQEBEJyR1dStIMSQxV29Cc2cwbDkxS2hFSDBkWEw0Uy4pV2dfNDEyNGBVSERYJ0JQQ3MtOT87ZUB1JXAvLkgxbjVNNENmXz1KO1xGJT1bZHVXOD1WQWwiITA1YlpzPz5NTzU1bUpFUD4nUVVUPzc0KGRgNDtFcWhGWmF1MXRHNi5OQVJXSCdJUSphNCs7TFc2Py5CWyMwTFlwbFtHU05Lai9OLF9tRyJXISI9XmUsM3FSLTMwTilfUUs8UVBWMkBSMlxyaSVSY1oiPForI2hkVlVQUiItJytzY0ZyXHBFSisjRnBtVmV0P0hFdCNWczk2M1JVSTIvZV9CZTNVUUZFRWUsTFFUTDQpOUEyPkw3dTAkYU8lXWRRIjJsTUYycy9vZlFzQ1xRN2RlPGZDP29fVUlFL3QvZlVyZWtraUlQSFRXM1xXNW4iaitYQ2JHU0MiX2k2VixYZT1sOzw2K1gzUFBvPmJpdSpMNHNmTyNUJ2QsMS8sVF8iPm5OJ2NSWkwpWD9qTEpmdDgyLCsnSGxjT2Jbcmo5WyIpJSVTYyJpM1ZwdCtcNzhGS2NtSz1cSmlwNEk4KidsNS9dNDArQUR1T1xANDJHKzpXVGcjWDFXPThLbWpbSHVmUSsqRm08RlVeMiZOOFdZIURgYCg+SSM+dWhXMCFnPkg7NV80VVEzIylHOko2NWFKPWhyJzUlM2VQMD9JNUlYL2hpTVxHNFhHcklWXDFXIzonRmdUInI8XiY9Wy5wXCdNV2g6alRPRmAxSkA0N1giZGIjPTojUzMwU11XcC9QcGN0TlxgV3VeK246XzJiajQ1YihpNSU7bUBjO29eOTlDXUBWSSdkL0FES0FWTWdwNV9PWjk3I1xlMSQ+O0RpRjs0MlVBVjMoZCRjcylWVWcqViY5JGRcYj1zOCwoQUtxZk1qai1pQCg7XDRxUztJKmxJI2RvaEJcNERwVCkmZUNXT1coXFoqNjY0PFkkKXBsNXAwJi5nUF84aHUtQ09WWGZia0lfUWwlXGM3JC9ORz1wIXFvTUxBbFYsXlBNTFNERFE9ZiplZFdbYFcmZFYrcUg9ODM6QGs7b0txRUY/LFVNVSVTcGlrJGc+PTgqQyZzZzQhKWhRSFMrczVhPD9RKlA2UDxLTFVmSiMoY19dSEhqOCE/ZnIvTTgtQCEpSDdIS0swZFtDSVxwblU1YG8mM0hMOkNMZFxFOmxyS0klKGdrR3EzSzlLSTE2Pmc5Oj1KXzU6QXNdLV1GTWQsUlZqTUZYIz9TTVRrLl9wdClhXTFrPitDPS9pOTwnY3NEKVFdQnBBaCVOUSc8J3VmcEImUTxkVnVbODxLZkJxZF8iO1EsWWNAIW05KTU6SzpTQj1NKXEsb0E+RS4hMVFzWk5DYkZVIS4tR2dPTTZpP2xBcEU9Yj5nZCI3XXBiQGI1dTFtT3A2KGkzOEspSzM9KDQ9VDclL2VmO3JzSDkpWFBCSzhIVXAjODosRidhPEQpZU9NcHMsXDYuLUBkVWArVylKUS8hPUpUQXA4WWokN0JJJl1FSGEjT0wuQVROQCRsJSxkVmJwIiYySCoubm8/ZylQJTxlNihWRiNuPjRBcEdSayVMXXVXMDEuZU40ZCdlIl86dS1XTTxiUylRQkxKUT1AMSstaERGZ2kvYFVWaVghTEYlPFdrRVRBdUlGcUAwdXQ1OG4/aFI8KClmTmJlM0FwVSVpKVArIk5QbHBVUFxpIlBnWFRCV187IlQnbjNLU0ZZOEBtJEcuXytzJ08jIVFKWElAYzljNitvcV9GaFgjV1JFcCFfMiwtQExNbCxNUCcnZHNRcGFKNE5fQCtESF1TX3FUKmpYPUNbQ0dVZGo3XjkwSEU9ZkglIWhHTGZcS3VZbGBbbk81SWVXJ2duclxFJXJPOjA4YVpVM0R0UTU8WF1iPmo/WVltdHJFI0chP2lQZ2xnIjQhLFEkclVYPWtWRUFRJ0hkcjQnQVwxajEiLmw8QWNoYG5dOCcvS3JKczQsRUNnOFFjP1ZlTk9XOzk+SDp1UUxxW01VWmk6JEJTMk5Bfj5lbmRzdHJlYW0KZW5kb2JqCjExIDAgb2JqCjw8Ci9GaWx0ZXIgWyAvQVNDSUk4NURlY29kZSAvRmxhdGVEZWNvZGUgXSAvTGVuZ3RoIDIzNTcKPj4Kc3RyZWFtCkdhdTEzPyNTSVUnU2MpSi5zVF9eVzYubmFnaHBlYEQmcCdmXzliOVtFVEMpUUQwYFJyWG83bystJlMjOjg7RFdwLCJUVG1nITw9VU0oImc5LFQoSTc3Ly0ucmE0UUhebS5GVCFLT2hQMi0/USE4bTxlUWhwPyhoNSxTcm5sRi1ibjBhP21AIzlxRyMmcmQ4cT01M2gkcDRdQVknODxtXjdnVGVzOjhhPzthQUtQaDEzNCVEcUBoXDcmPiVjVikvMj1kaCNeMVUiOkNzW00qRV02Y2MuVktFXW0tbDUzKixfSilnZj8yQSQ9UDJuLEBYVjJZRkxmUyolNGJTXzw3KFhtOF1wP0NnQ29ycl1uM19KV2A3WVJyVGc9Vzxva1RSWCRhLk9JWkhybToqKW44S09MQW42Nz1sdTxTN1InYzdcaXE8QVo/LWZuQkhHIi88Q0lISSxrM1NBW20rRCRxQV1xO2JfZS0nRENlZmpFOk0lUC8sM1B0ZnAxSCJDKTRtSy9eO1NMZmRDa1IsYiFobk5nM25qVlQnanNpLV9HUlNpcGlwLjF1WywjPmZCWHBVKDxNcWtzTGstKyVmaSc2Y2hiPU1JNVZjU3Ndc0tSOiMyR1ZNZFg2XF5INzRNYSFyTGRZQU5wNTIqajYzaEZEL0k+SU1mKzJUVComR2YhQGgrNG9nJzNHI1koM1NcPmc7ayVAUXRzYUI3Mj5BcycyU0BWR0N0NmxWPFAuLl1ZITt0SiZWI0dUPlwpPUdMZDwpXig3IT44YypBLkJZLCQpPWtbM1VsbyEndD9Ncj0+JWFEYEFULWtHQG5ZSU5ITSpVPT88YVQ0O1MqZjZSNUR0anNTbic4THBaJFtFUzQrcm4tT048IlQ+KUdJZkJfaGhGSUVaYy1HdG5oYTNVZVIyO2UuL1dAaC1cKnJkRkNfNTdJSE5iaWI6Mmc9alZtVVVxKV9mZzBvbXI3Xy02amxGUyJNKm4uWCpTSD4/Tzs/LT8wKkxnUTQ3RnAqLm8rJkRJJDM/MFgvb2c1RSZoOGNmbVJHYWh1SStua3AsaypsXkByIUsjNltUblQzNEBHSzVZJiMwQEhDZ1ptKDBvMGk8PTdcciFxODl1bFY+T2NXOmsnPEZebjhkL3Bxc1ctQzpjJi45c15aSUBYKyM0KDtYa24+L11WKigxajpVJSolKl0mKXBGUmJEMGg+U11cUkpqZ29gbVRCRVBOO18xVydQWipPUWZmI21dOS1YV086cygkLUNHOVkwb1w+czxnSF82MU5gXDdqXWhPST9maTQtP004PG5xV2VgSVkhUCshWWZpXWxsSTtBTFApVyUqZ0tOMF4lRWwnPEklNkknR0dfIj0sY0sjNWtJVyY2VTVjbi5dJ1UyWFcpbE9MKCNBP2lsYl9cVVFwS0QuLDxcLihtJiFLYT1gWFVGVF8jcSc1MzBQVkxxcU4jLlciZl0pZS1aYy4lRmc5TyIxNSdmJy4mYjtyQyFYP2MqLiUrVGlkUFRRQHFNQGZOXGA5IzctVFdrIjZKN1VATUZEMjNeY2M3VCdVXis/KitEInFRZUUsIysxUm1QQUtnTGYsY3Q7SlNgRTJjKHEtVmo1Q0tjTF9EJ0lPbjpcayIlKVlKb0Y6NShrTEpIQ0YwO3Q6bWY6ayRBJHVWbyFePDAoZ14pajZtK0ViK3BKUj5SKVlRaEpEQScqZVZPdXJUZEApcCJFRlUyRCpgSnU1bmRuLks1cidbZ05xaSdJY2g4cC5yI0ByalBESyU9bEI0LzRmODJCMzBQO0NIYEZoXGpedG5hMmIzLVosYWI3JyJXcmdGMzBKbGFWRU4nbTFFKTQhP1ouaEtyRCJcMT9eTmY2TERWc2wwRWxldDksKk5OMWAxRkBlYkBWLjRCQUFDLWhDaT87aSwnNkMuL0JnTDhrQic5clxtclM9J2ktU1omJiMhSUwpPkw4Kjc+NV0zNWVbK0RcLTRPK1dHakEmXCIldG1cKFFEKDdCKjA7LD4nJ3FLWDE1ZHNDIVAmbC9vKWtbYk0uSilNIlVGK0NXV1taQ1gjK2ZDNDxTakJFZGYtUS5eOUNOWGxQI2Y6bUdGbEJwQ2EpWmo1WWU8NmMjLToqbDQiVzpfLWwpbGFQPDI1c0pDUmVDLkVNOT5YRmQmRnI9J1hiNEtFXk5pWV1yImVlOENKVz00JkoiJiYzPyhZTDhHYDI1XTNNbU4jblhUOCIwLVpna2oyTEpeQ0M3MSxRRF8iKi1iL0BXWWRiKmJYJlZANCo4XnFxWzFLK3VIXldLNCtqPkBXOiUxKDIoTDpAN00oOlxdZltwQlBOZTtfWmo7X2kzWW9kQjxIa3EtUGlacipcOVdWJFRvN2FjbHJaKXE9QzQsKCw+MG8sLjRxNTc7Mjg8XD4mczs9MHJNJjcvRGlQSnVwLW9mcCZEVF0zVGhHZV1nOzlOPVtYMmUiPXBCVyJPYDwvQ0spIUdMUmZXX1ltRlRlZ3A9aUE0cUFoQmVBMVNVRy1rOG9gJF1UWUIwJ0AnXiw2WWZNQ1ZoREtSKC9tUUgndG5XMVc4KyYwWTZyY0A6MWdAJFtSOGBqcGQ/dVBoZV45MTkmYmhla2ZdL1BsO21TMjNnUjprLi0hbERENC9XIzMkbydHIzlgcXJCP1Y2PiFla29NZDNbZiU5aVdyQ0QwYSdLVGNPcSZHQlZQaE1QLGdwT01EOC9kb1NgRltYQ1xxdSYuJCZLOy1pOzliRGxEO00kMiQxMW8lIS8/OWpRPCpXK2liN09LNWZBP2UqJStabjVXYl9xcys1WT1iRTY3cXBiVm5PbzwuZjJsaDQ9MGoncFtYOXFRZFZSKi5jJmNPQV44Y29FZ1tcQS9lKVZwOnFfV2FuWlpqbV9YMSUsanEyZjJtXjcjT0dbY09RJThAJ21UPEM8M2xxNVgoJEVnLm51Ul5mMEgpcighYmdbYG1COEAnbVQ5bnRSW25sN0tBZnA4aSwraklIQzxWZU80OGVYQl8+Z2VLUVhMTWcmNldBTStiKGNTVVAlMDRgMVx0RExeUFxmSGs3P0hgSixKXy5jOGs9PTMuRVJTQixTSVNEZ2s/L05xKDxHay1zdHEocilTVkRUZGgnTXIiOkFHKDxeJX4+ZW5kc3RyZWFtCmVuZG9iagoxMiAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAyNTEyCj4+CnN0cmVhbQpHYXUwRmdRTFUqJlVnOlRedU5GN19eLlIiUVpEMXJQUCRuMW0jKDIkYiM+PD4mLihwNVBHJC5DbSlpXC1ETSJaUFIxTUMiVzI9LzIhaV41WlZdNEhEZlRgbzg2YkMqKTYtLnJUYW5tR1gtbVImIjZASG0lY1RkVF1TZmBFMkUkPzIqREZDIzlqIjBvXzVFQCNOUEApKFNKQGlOWC4sbHJ1MFtDJS07TlIhXDo9LWJ1YWVzbjgzKDJQOEtRWzI2ZUZ1PXNhMkhwaDlNZVhTbkk7aWFYb2tQbScxQ29sK1JkJENxYXRKRis2NXFzLScuLEJuLE11SUckLUg+WSxodVhPU1VBaENtUUhCYVwpNEJvKStuPEl0TThKTG5QVWRGYSpyYFplLnExPy81QWRYRWg9SDEwQkRVSk8xYWZBWCFGLVg2NjJZQCFoLC9jPUApO0hdL1IpXUIjK05VXVFhbyE+LC9IMCNNLnJRXnAqbCFZdUNyTmErNUxwY2k0Lys1Vk5HaG5XYGFxYiJMRFw7QTgwYyFTP1pYcXJxNkEoVGlyaVItcnA6JTA4ZD9BRCVHUyxUKHMvLF9kMDluLSI+UTJcU0dccU5YXWwmWVY6NC8sZiYsSjw3THJecFBtMmRPNytsXFhyOFxEOyZRV008KnFETzg7b1IpNmpqaWMqbVNXamhTNV9hU20lP1ZGPCNeYj9CZCpsRmUtMF5vIUkoPDk8RWVQJWE7QVpjPDs2Wig2Tix1PHVuMTJkUDNLUCw7Y0tYVVtVWSM4Qz9UZmsmTDRCcSRWMz9ANVksU3RUPyxrLShtWUtgdHQzKGdPRWkpMT5HKVljJyEvcUxQOyJzO2ZDISV0YSdFU0h0S18jJCNfalRZazBwOFlWOzE7LUloTlwjXFJMYUUhUDAxP1FMJGY4KWpwUyxIX15YQUxTZ0ZlPUY5OiRWNz5bZHIuUFkwJm5hXjI4IlBORy1cJVVoKEwzNWdpaUYsNEhFWEVUSUBcOTpGO3BVM1RNUEtTQTpATCZFOEZYJ2trbSMzVlQzNV5jI1YpLkpSPSQpJVUpdSpJJWVeRDE5X3BvJmRNLjUtV0pTNEosKk9BQCw8SmYxZjJzZT5IWl1uOFNpcyVVI0QjXmZxZUVzXWNUUUhFNSIlaGZJPHFaKT0xWz5jPzJMKjNkTz5LVlM3QjpYUyRtZENHWEsnRktVMz0mV0tZVVdYOSVGczVFTDBZczlRZiFKSSNYIj5QYyxlNi5MIiFGZz5QZStVZTwwbTt0NzQ+P19jVGgsc11pMC5ZZltXQlNlZkpkQD9QTCdqWDJXWlFyblFpWEAtNSh1OyloPWxOMVVvRzFFJVooWT8oUC03XSVOWTU+I21vKnBPRUAzNik6Qz1kbVRxVCYiUC03YCZMKFtNSXBTZmQyXC5dVlFQKT1nYVxsYi0hOy9xPSc8VjVgWT5xQnRDLkwmTnFFXk91UyxbWz4rJG1ocHBtU21LWUUkbSxpOS5HN2JsPSo4dDdgPzY7PDA/JTpdJUd0OjwvZCxCQT5cWCQ4JXFJMykrbWxVXWkwJ1lfayFgPjRwXm00XHUrYiU4JXFJMykrb0Y4bVNkQlc/blhANE1eT2tEOmczNGE9IiQkcSNSS0M7bEE/Vks3KTVXZDk0WHRORzE4I0ZTMDgiPElcMmQ0IzxOTjpsK0BXRyRBb2tzJDM9Y0YhUHVOSDI2aGVZYjw3XzRHSCxmR0woW3Q5XDNkWUI2M2VUaWgqajdhPT0oPFBAOVpjKmhhTWpHX0k6S3I2RyIsUG1WJDYjcEVJaS1qOCt0MSVCTk86L0BQQCQ5LEErXTxvOWkuNEA/SkMvW3FBNk1CZWdbXycvUXNUQiFOMF5QaGlmOyJtbDA5bW5NIkFcSmxUMlhpWj4kZVA5cWxHOHAmak0lXl1eOzc6cTJFQzV1bkdFJHQkMmEsamxqXEpkTVxDZnNRZGJVUkghT2Q/RFRYRT5EPzFfYzQ0O0s2O2lwNzZYKFVASlQnPitcaClCYUBeMmYiXm5hPyI9NDU6VzVfI1ZnazdnJ1RLTFw3IzFuVGA9TUhtZTNfcCNHUzJjYjBwcis5LSVMZjJjO25PZTBQPjMjYD1YaWI0aGtKbWdQQ2NFaG1wJjRiY14ybUQ3U1M3Tm9hXTNxQlkpYTFOOTVFRERUZDVXJ1NDY0U+R0s3T0olaiIpOylFJj5sOiF1bTk/PVozUmokVDomT2Q1K3VoT010azhCVjo3LjldODVOU3BDKXJfcDluNTVpUS83MnAtP0xTb0hkaUFJUStaIkhAbFw7PCs4LkI9c3NLYF1aaTxHWFksXUhdRlpyYXRNQUchLjQ+VExwWVAiMlZHbF5mWUBwZmVONFYuXnBvKUhMXSFlZjZXQDMwLCN0LGlHZWc4R2g+Uk46J240Yzg3ayVuVlJZW0s6Ok44b1FqUD48Ilc1JXRMOCthUlpLLV5LbHBdTVw2OSVBX2thYmw8UVcvZFprUD00IXBZMEhQTyZMbmY9Vy8yWHM8M2VLOFYmdHJ1LWxQQm5BVERcQzlIYV9VTTlUZDZEWDVScFs9SlllX2ItKSU9X0NKJiwvNDprQ0FqKyJdKjMnaD9HVHBMYiduOWpQK2ElSzQpbjJNbnU4cVpcaCs0Lic/XjRWWEQoOU5eUFxmZkZbOF89KyRpMlxHNShJSk1GYGpMXHVzOihAc1JdaWBjX20mKzN0cmZVXyQ3Y0kzZitxJzxRUlxJYVBhJGRDTCdOKiYyLTJFc2cuWyYuImtgUlZcPixLaFInXXAhUCU9a0hGcFMsIVAnSlc3UCZ1ajVcTUxnPCFJNywkbE8vQCROVSYiJUxDWkxeZDYwIz9jayM4YE1bZVFdcjIkaElMPSkzb0hjRU9pKTU6aFI0PlM5Uyk0aCUlKGo+VTNsb25vMHFebE85M0M4QyU1IiZoSWM6bFJQaWdLQDpUOWA6R2stJTcoMDxlK1NPMFtTYHNKayZxX28pXlNVLT4+V2FSQS4rUkBmOj09cjglVCo4IVtqbHNGXFo8az5iQ2lbS1JyL1NAPyE8ckZoQykwPzglMTs7UkpLWUVEbVssIkQ9ZzxQcGYiWSg9RWRmMFQwR0xPSmFWKjRIbjpKaWdGLUY3ayE3QW9lSFVDIWFDITFkZipjbS5EbVo/KjFfIiUvRlAoN1RqO2Zza1ArVD5BckllNTo6azVvcFc2MDplMVhBT1teNnRBJmpENGBhM3VQWiU7R2o6RCsjO1lXPUA5azFBS1E4J1I+QTs7ZytzVUlWSGRyYD5PXUVzOSdXOWdQR1gmaDVBWT0xOTZvN0lTZWRgP2A5aE5XLU0kWGovNmY0In4+ZW5kc3RyZWFtCmVuZG9iagp4cmVmCjAgMTMKMDAwMDAwMDAwMCA2NTUzNSBmIAowMDAwMDAwMDYxIDAwMDAwIG4gCjAwMDAwMDAxMDIgMDAwMDAgbiAKMDAwMDAwMDIwOSAwMDAwMCBuIAowMDAwMDAwMzIxIDAwMDAwIG4gCjAwMDAwMDA1MTUgMDAwMDAgbiAKMDAwMDAwMDcwOSAwMDAwMCBuIAowMDAwMDAwOTAzIDAwMDAwIG4gCjAwMDAwMDA5NzEgMDAwMDAgbiAKMDAwMDAwMTI5OSAwMDAwMCBuIAowMDAwMDAxMzcwIDAwMDAwIG4gCjAwMDAwMDQzOTMgMDAwMDAgbiAKMDAwMDAwNjg0MiAwMDAwMCBuIAp0cmFpbGVyCjw8Ci9JRCAKWzw3MTE2MThjOTI3ZTRlZGVjYmNhM2VkNDQ5Mjk3NzU4MD48NzExNjE4YzkyN2U0ZWRlY2JjYTNlZDQ0OTI5Nzc1ODA+XQolIFJlcG9ydExhYiBnZW5lcmF0ZWQgUERGIGRvY3VtZW50IC0tIGRpZ2VzdCAob3BlbnNvdXJjZSkKCi9JbmZvIDggMCBSCi9Sb290IDcgMCBSCi9TaXplIDEzCj4+CnN0YXJ0eHJlZgo5NDQ2CiUlRU9GCg==",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for day-specific annotated cell-by-gene-and-guide AnnData matrices from the WTC11-EC 300-gene CC-Perturb-seq analysis.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 09:25:59,104:iu_debug:	DRY RUN is enabled.

All Document dry runs passed: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/DOCUMENT_DRY_RUN_PASSED.json

```

## Code cell 22

```python
REQUIRED_DOCUMENT_CONFIRMATION = "POST_3_300GENE_FILE_FORMAT_DOCUMENTS_TO_IGVF_PROD"

DOCUMENT_PRODUCTION_CONFIRMATION = (
    "POST_3_300GENE_FILE_FORMAT_DOCUMENTS_TO_IGVF_PROD"
)

if DOCUMENT_PRODUCTION_CONFIRMATION != REQUIRED_DOCUMENT_CONFIRMATION:
    print("Document production submission is blocked.")
    print("Set DOCUMENT_PRODUCTION_CONFIRMATION to:")
    print(REQUIRED_DOCUMENT_CONFIRMATION)
else:
    if not DOCUMENT_DRY_RUN_MARKER.is_file():
        raise RuntimeError("Document dry-run marker is missing")
    marker = json.loads(DOCUMENT_DRY_RUN_MARKER.read_text())
    if marker["document_metadata_md5"] != md5sum(DOCUMENT_METADATA_PATH):
        raise RuntimeError("Document metadata changed after dry run")
    for row in pdf_inventory.to_dict("records"):
        if marker["pdf_md5s"][row["format_role"]] != row["md5sum"]:
            raise RuntimeError(f"PDF changed after dry run: {row['pdf_path']}")

    production_results = []
    for operation in document_operations.to_dict("records"):
        if operation["action"] == "SKIP_EXISTING":
            production_results.append({**operation, "success": True, "log_path": ""})
            continue
        result = run_iu_register(
            profile="document",
            input_tsv=Path(operation["input_file"]),
            dry_run=False,
            patch=False,
            label=f"document_{operation['format_role']}",
        )
        production_results.append({**operation, **result})

    pd.DataFrame(production_results).to_csv(
        WORK_DIR / "document_production_results.tsv",
        sep="\t", index=False, lineterminator="\n",
    )

    verified_docs = []
    for role_key, spec in DOCUMENT_SPECS.items():
        record = get_with_retry(spec["document_alias"], attempts=10)
        if record is None or str(record.get("document_type", "")) != "file format specification":
            raise RuntimeError(f"Document verification failed: {spec['document_alias']}")
        verified_docs.append(
            {
                "format_role": role_key,
                "document_alias": spec["document_alias"],
                "document_id": portal_record_id(record),
                "portal_url": portal_record_url(record),
                "attachment_download": str((record.get("attachment") or {}).get("download", "")),
            }
        )

    document_records = pd.DataFrame(verified_docs)
    DOCUMENT_RECORDS_PATH = WORK_DIR / "document_portal_records.csv"
    document_records.to_csv(DOCUMENT_RECORDS_PATH, index=False, lineterminator="\n")
    print(document_records.to_string(index=False))
    print("\nSaved:", DOCUMENT_RECORDS_PATH)
```

### Text output (stream)

```text

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_guide_quantifications.json
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/document_guide_quantifications_submit_20260930_092559.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:25:59,745:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:25:59,746:iu_debug:	submission=True: In submission mode.
2026-09-30 09:26:01,084:iu_debug:	
IN post().
2026-09-30 09:26:01,111:iu_debug:	Validating the payload against the schema
2026-09-30 09:26:01,129:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_300gene_guide_quantifications_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOCAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDcgMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL1BhZ2VNb2RlIC9Vc2VOb25lIC9QYWdlcyA3IDAgUiAvVHlwZSAvQ2F0YWxvZwo+PgplbmRvYmoKNiAwIG9iago8PAovQXV0aG9yIChFbmdyZWl0eiBMYWIpIC9DcmVhdGlvbkRhdGUgKEQ6MjAyNjA5MzAwOTI0MTItMDcnMDAnKSAvQ3JlYXRvciAoXCh1bnNwZWNpZmllZFwpKSAvS2V5d29yZHMgKCkgL01vZERhdGUgKEQ6MjAyNjA5MzAwOTI0MTItMDcnMDAnKSAvUHJvZHVjZXIgKFJlcG9ydExhYiBQREYgTGlicmFyeSAtIFwob3BlbnNvdXJjZVwpKSAKICAvU3ViamVjdCAoSUdWRiBmaWxlIGZvcm1hdCBzcGVjaWZpY2F0aW9uKSAvVGl0bGUgKEd1aWRlIHF1YW50aWZpY2F0aW9ucykgL1RyYXBwZWQgL0ZhbHNlCj4+CmVuZG9iago3IDAgb2JqCjw8Ci9Db3VudCAxIC9LaWRzIFsgNCAwIFIgXSAvVHlwZSAvUGFnZXMKPj4KZW5kb2JqCjggMCBvYmoKPDwKL0ZpbHRlciBbIC9BU0NJSTg1RGVjb2RlIC9GbGF0ZURlY29kZSBdIC9MZW5ndGggMjI3NQo+PgpzdHJlYW0KR2F1MEU9XT0/PyY6V2VEb0hRLlxOOS9lX20vdE8zOVsqUD5SQC9Kck1BITBiNERXLVYjPkQpcG5fN2UjXFZMXXQhIVNLMnEsO00nO2c1Ols1Zls1JHMkbmZaaGE/JFxBXFxdazRvaG9HJFxpJElBQFY2QFVOdCgtRjpVNFNaamtwMDpZWUlMJGEuZyFFN1hXQCJPay9AQDs7cUkwdD1acTNScSUyMmBoYllYKVBgVjxqTzhQYkFxYVsjTDFpKUswX1AhMDxlJi4hYlM9bEdWRyhicEkoJD8mYHRiOVEvQklaY0hnczJfRiMuSEdZJDhNa0UmRVcuIXU6ajQ6YCwkJCY2MydwZlQyV0pFV0ZTMCd0OUE/TGI5I1tQN1M/bEopWjwsZG85ODlyJz5WJHRfUkwvOT5wYjZWRkBKY1dMaWAjPjUqZFQsUEwtVjZzX1tWJShKckNqXGcwQDAhb1A3MlE7VC9AYG5DTTZnQkFfOXVSJEJqRS4pWGI9UFNAQ1ZEMmRMVWZkbVszXjo6aEQtaiIuc2IoQTdELihSQkMrJjZEN05lWzEsJjZGcTowNFJVSDtjVlNkb1FDTyFlTWk4Yi87PnAjI1onUUEuTlwuV0xeVV4vPWgpbSxFLWdSc2pNWE91N1YoIz1zInB1WCVqRj4sUGAwZkhGOjItWyZnYjs/XkgsdVhoTCMqPVp1akFHXVhiKUJtPCRyP0sxQFYxcjAmbmc3bVw+RDs2IihqaVZRTVFPb2xGLltVb3NKZWE4JF8nXUJ0Z0w3XjdkYVZTUGpDWTotNTFdXUdKdWJnKEUjTCVRXiNQWm5oRlpJRWlXISMwaCVaWi8qUSMlSHU2azBXZD4/JVhKOCNyK1RaRDtWW05vVDVFP2ZCKScnVURgNStVJ3JpRUAzWlooRWU5RFNfWTBXIylkXTE1TVM1KDpzXks1JU0lblA5RjcpKTBUQGxsaEQtZl5iZz9fP1tOcTxWMzFdRFVqODdLKEdYX24qS0BxLFRfbF0nRmYuXGMxLVoyNGhPUiFsNjpJXFdBaTRaYDMuYl1BRyp1MlA5Y0UyQ0BhcmNJSEhrK2Fpb09ZPkQ3S0RsRmhyZHQlYjdudUg2QjtrQW1JRXEoTjRrOmlNbksiYzdMR1RdbVRQTlAoNzZuWm5JUFtUbFEiITU0MTtgJmVfOWdxUD5wMydOXFpSOVZmRXIwUjNdVURPREdgOCx0Ty09Z2k9O2NmakE2USVhUitrLHJsZjhvaik8TXJfbDNEY3BRQjZyTiZMal8qVTFkQURKWVY0N1c+KCMiQnF0ajxKc3BXcCMoJlYhWV1dZ0xfRnQ6MT1HVi5OYnMmdWJsW3JbXFpJKShZXVNlKGFuOnFZbEJUIVErNSRPJl5XKVApZGE3NSwvIWc0N3NKZGBQRDslXyI8NVQnXEgkQ1BQcmA3I1VaY1A0QjctTVcxQl5OVk5ETUElSTMoYDpzaD0yZCdpalwvOy4uJ1Vaaj1NVWZWaSQ8IzVHTnFoTHFLdTEocjYhZF5dW3JeaERgXltzbTUyS25rMz84aGYiRHBDbDxcO1hda2Y7T1gmaWVjNlVJYHQ4cWNaOkFLJGswPWwhOHJNJGAhXTEoWmcqYkc/NjZSJVVVWStWcjQvZlBlbCEtdTYmLDoyQVJmMT5zZzE2UTc/dSVQJjRMUHJbPXBpZEpZNlhMXC5TZ1JAIjVwa1pAa29AXlFvNDxGSG9MTWJHZWNyazlzPmA2V0g2dUhpR1tyb3I5SEQ1bHNkSWtaUEFvUVlaT0JNT3RPMlYyKlxddFZjP2NaZ0ciUi0ucEtxJlgtb1lcdUZbKzNVQ09hYlpgUTM9VWk3bUNDTEY+OVAnXi1qXUc7PzRNYjo1TCh0SV86O24tPjRZcW1eIyEzVGVNKnE7SHQxZi5BZlZhNmpAQ04uPWNGZlBGT0QyTTlRLVtWXVYkOyk9O2g7Iz9rYXUqSkVfYHJbLE1USzczcF1xPmtkUlhyaE9hTFVrLFk+YnArMFdVXjE6UjsuPko/SHArdEhJN0ReNkdlVk05dUIxO09ILGJNJV4kPyk+amYzTypZZzokXUAsKyI9WDFsTmoxJENhPSx1XkRyN2hPZlxOWzU2ZVAwXDkxTHRpNEY9VmRvciQhdVpAQ2NeJ2NTUzMuaXVaVm9WQ0A5dSpxUG86TE9GKlUjdXBeXGpgcTQ/XGojUmddYUhbJ1JkbnUwV1hQTjJmMkkoXj9JLy9AXTFNUm0wV2BYI2Q5bmNaZDYtXEk+T3FtMWxMbCNUJUVsUlVSLzNsR0kpc0AnPDZBUSpYLmBVdSZ1T01nMTQ4LElfOjdCak8yP0UiI2BvZz0nZGQsPFcsRjE6JFcoXTAwYTRUQGtAbW9JWUEiP1VcJy90JEdvK1ZhLmAzXGFRaVk3J1NqUVE2Qi1OLGxUZSJdVyRaKEU5TXBpbz0qSWcuK1QwKkZtb1krV0M6YV1qXSdmWzw+RWlSbG4zWWk2OyMpUzpOPEphIUVxZjs0J2RkKkhNZlomKyw9UzlDI0RIWG4lJTdeWm00OnIyVWhgN0BnITVfdThXUVNBcV8vW01hPV4+TCtwMytBLFheL20oKlJ1OXM4Pi4ibDc9WUglRCtTaUw/NmhzIyZnXjBob0BJPTdhMShaXz5sZV04cC1PSTxwIks6US5kPiFkbUhBYlY1WzdXSigvTW5LUTA8OE1LZVgxPC5QW1QyRkBzPXUybythPS5QRl09VydhYk5EbF9aJEMjcTlfOygySitqMzhhWz5oN2BaKTglRWRkMl0pQm1mKGMsQnVOSVFqPlwibFp1JmUvJSckM2U8Qmw1UTJMaUtrQU1AYWlWc1o3YnEuTitsZU9ZUzknPDc4QmJQXDgxYzteSi42WEVPOVVfbGx1cDpGNXJdLVouQj5SJU9DRlM8IzFiQ0hMQTNyTzcxaWZvM3BlYGZwTCpvUzg5Ll5uMmwmWy8za3BnNWE5PjVHZ2s/SSxSUT5aK1psSk84Ql8qNk5nbmNgPVhuY011W286PDx+PmVuZHN0cmVhbQplbmRvYmoKeHJlZgowIDkKMDAwMDAwMDAwMCA2NTUzNSBmIAowMDAwMDAwMDYxIDAwMDAwIG4gCjAwMDAwMDAxMDIgMDAwMDAgbiAKMDAwMDAwMDIwOSAwMDAwMCBuIAowMDAwMDAwMzIxIDAwMDAwIG4gCjAwMDAwMDA1MTQgMDAwMDAgbiAKMDAwMDAwMDU4MiAwMDAwMCBuIAowMDAwMDAwODg0IDAwMDAwIG4gCjAwMDAwMDA5NDMgMDAwMDAgbiAKdHJhaWxlcgo8PAovSUQgCls8YTE2YmI2MTYyZjc4NWNhNmYyZDYzYTI5NWU2MzQ1YjM+PGExNmJiNjE2MmY3ODVjYTZmMmQ2M2EyOTVlNjM0NWIzPl0KJSBSZXBvcnRMYWIgZ2VuZXJhdGVkIFBERiBkb2N1bWVudCAtLSBkaWdlc3QgKG9wZW5zb3VyY2UpCgovSW5mbyA2IDAgUgovUm9vdCA1IDAgUgovU2l6ZSA5Cj4+CnN0YXJ0eHJlZgozMzA5CiUlRU9GCg==",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for per-cell guide assignment quantifications from the WTC11-EC 300-gene CC-Perturb-seq D0-D3 analysis files.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 09:26:01,518:iu_debug:	Success.
2026-09-30 09:26:01,518:iu_debug:	Object posted with identifier: 4ecb21f4-7a3f-49d3-bdbc-3161d349311d


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_per_cell_quality_report.json
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/document_per_cell_quality_report_submit_20260930_092601.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:26:02,079:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:26:02,079:iu_debug:	submission=True: In submission mode.
2026-09-30 09:26:03,597:iu_debug:	
IN post().
2026-09-30 09:26:03,624:iu_debug:	Validating the payload against the schema
2026-09-30 09:26:03,641:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_300gene_per_cell_quality_report_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOCAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDcgMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL1BhZ2VNb2RlIC9Vc2VOb25lIC9QYWdlcyA3IDAgUiAvVHlwZSAvQ2F0YWxvZwo+PgplbmRvYmoKNiAwIG9iago8PAovQXV0aG9yIChFbmdyZWl0eiBMYWIpIC9DcmVhdGlvbkRhdGUgKEQ6MjAyNjA5MzAwOTI0MTMtMDcnMDAnKSAvQ3JlYXRvciAoXCh1bnNwZWNpZmllZFwpKSAvS2V5d29yZHMgKCkgL01vZERhdGUgKEQ6MjAyNjA5MzAwOTI0MTMtMDcnMDAnKSAvUHJvZHVjZXIgKFJlcG9ydExhYiBQREYgTGlicmFyeSAtIFwob3BlbnNvdXJjZVwpKSAKICAvU3ViamVjdCAoSUdWRiBmaWxlIGZvcm1hdCBzcGVjaWZpY2F0aW9uKSAvVGl0bGUgKFBlci1jZWxsIHF1YWxpdHkgcmVwb3J0KSAvVHJhcHBlZCAvRmFsc2UKPj4KZW5kb2JqCjcgMCBvYmoKPDwKL0NvdW50IDEgL0tpZHMgWyA0IDAgUiBdIC9UeXBlIC9QYWdlcwo+PgplbmRvYmoKOCAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAyNzE2Cj4+CnN0cmVhbQpHYXUwRj5CQVEtJnE5U1lmWEVtXk1dbT8uQUJtUE9hLWgldDFgKiFaZSsmQk8rc0RlKyxVZTswPyRaOyxdKmplNU1qLzNsSFllZmVvak4oN0BCbUVFIXUkMFBxQitsI210VT8uUTMiPSIzISYuWSVMYlM8cXIsTm1LYjpROTFCVDw0Z1UwP1E1TVtIZCo7SzZkIW9KPTomN0FeZUBHLUQoNW46VUopSUBLT189cmlsZkksSF4iYCQoLFc4bVpNKWhSJTtmWUJhSWttT1IpOichMCpuLGFIaEU7PGc2LHU6UWhwXzhdciEkPUgnSUxEPnJqIk1vbSpgWiVeSCg3WFgzMyhDbl1SbypMMGA3PiJnSGJDSisxcTUqNWVmXWMlVygxbVlHaCViITBZdVchXTFLbEptJD1qOmZ1a0pJWG1gV2dWRFBVZychMV0/ZiwmJiFuLGYuODg2N2clRl81cV9aWT5KRVkuZ25LRFMqJVY6PzIoVCw0ZUpgPkVGKVhiPkNFYk9RVVdFKVtadEFfPzhRNUIyJ2VHPkJkJSg1a0pXYGBQZ08hbUBmcmEpL2FdX1tdaiMqXU5JJDRFK0hfIllkaUlKKUhDSDIuTTo1ZjNrKjFAKi9ZLmxnIlJkPT1rZjhuYlxTUGdWIlYpYkh1R29ENDxwUVVeRFFsUVFbXFhMOl00TVxqazRUO1hLVmJuMyZEZCokVU0xYV5wP3FJYk5vJ1tLS1c9IlExdUc6QkgpYyw2PDFdLUheaiEzLipGVSkuRGE8SDElL0tjXSZNI0sjOzJSSlNCbFRJQzpAaGhmLThINylzVStVKHFicFpAKmBlQVcjM2FhczM5aFtPP1gwKUQqZVNoY09Vc1xodT1nPU5wc1pJMClELzIpJkFfTmYncCg3dSlYdVRmMl0/LD5NWCVDNUBwZzgpXEIkXVNaKmNvPlJiRTxUQ2liY1NJczJlKGhfRjtZRE40ZipUIXJnUS1Na10jL2BhNSRwVShEck88O0lmSlQ4XWs+RitsaTJHL1VWLjI+Xmo4Kj9vNHQ4VDM7WU45SyIpQCwoRFRvbTNVXWdMWUVzT0Q9JDNGQmMjZTouMTBeJihNVjNpOFRpT2RoV3VvNzY1bF9YRjxeKVJUPDpMKUByY1xPLEJmMU1FXiNiMkMkJ3VjbWEiJ10qYEZZV3VbJmhVWDpTNiZxUlFbOTk/N0xfV0dBMSU8Sy8lZWNIWmMxX1JPSUJGTmwxRSs/XSZtYEFYM1swKVAnbz1eNWFdI01GVWQpb188c2RHMSRLbSU1IjdYSXBVbWhOT2coKiQpQW0/SEdBQ28zblpSLWhPJnA3KXAzRUZvO0hVJkZtdC9CPiFqKT1QJV9EPlslXkdFTEMwWj87OWtrc2E4cikrPlwmYXRQXStWR0dkQHJdI25OLydfS1xXaExKNzFGcnIkdEpwcnNGOiRcM3IzPG9ESDRWO2E9OWBHJykzYXI/RCZRUD8idFYqV2pMRFJHbnVRZUtJYl4/cisuclVvbmpEcTQtIVgmYEtBZy9iayY+QixwQz5PZiI5QVNZRk9nRGRJXSxtOFdTLCgpTmY2b18lTmclTy06TjdROGY+Wl80LltQX01hVVgtVVpfYS0uY3BRUTVBZUZbU28/VFBjSzFNZCNNU3RJMjdIOWJKQG8vb1xdbHVaW0kjO2tcXUY6UDdYS2tWUShVXiopKDshaEJ0JXFEYjs8Mi1tQ2dqP1MkVy9IJipmUVorbjNAPCUjOiRjPD82PnI4Y2hdQ1xNcj9CNSZFdWEhbiVFSHFPcyMlZFc7RzwuSExTOTprZFI6ay8nM3EyV2olRiEhYD0mYyYtQDJXOD4yJSIxNG4sYj9TXnNzM1U1LUdLI2RJVHNrVFBrWVFRTE9oNT8ramIjO14oR1RjI0hqKHQyazYqUkUpLTpsTlNbXmJybFs8Ml1kPi9nTFhHbzBPX3BHVWxbSktJQTtiMVg9Yyg7MVYpWSdHUDNhP1dLLE1vS1o5OSprcG10U2pZTSNrUylVRyNNYCo1NGt1TksxO1BAQTNNdGRHK3FKVykoQztpRC5OcE4kVi1zKDI/IiZCbU0+W1lSIjM1KFEnbUQoIk1jKkRdKDZVNm1gQ01mREFOQ1FgRGIjRShyVCshKCFvKU8pZGZVVmc3OTMpMlRjN0RvNGpRMTs7b0xpMz87QkxOKnFUY09kZl0lJFNLQDQyUSxUSCozNy8sX0ZmTEYjXVFJR3NhPmY+PDc2WEZWT1hrQUU8ZiE2RE1UPzYpPjI4O2ZvSTs3YTJwTHRbLyg/IWVlT2BtU190NiFXQTQ6VCZdPDJJO0xPOTA+b1JRZkZWTUA1WT49MFxdOzh0QV49KjYsKGpXYCMiKWtKMi5lQyRwSC8zOkFKLk9aK04iYSxIKGg4SWFCOzdZTiRkSzo0Ly5bdWQibmFcTVtsMUwzc01HSD1eT0lCTyZmSS0vXWRrbmM5ZF9sJHJhTiMtXidiYUVARGxSZztQMDQpcWUjVUAvKS1FYjsyIlxAXlY/XS9DWCtqclVdVkM0WihvR3AmLksyc05xI2wwMiFWYzFZWmJqcHROKyFVay9nQT4+JyliZXRcWzMhKzc2cWBUS2FZdU5XZWI7VCN1RmQ5VldMTk5FckRSVDtbW11ZXmhmbmdjKytnVlknbWpYPFplNzgwS25jJDtVJCFicnExNl41JjxaXUg0NzUlZipiUiZLLSdjX2MoPDwsOkpZMHJbSkRrQjRZSVBkXkpGcDZba0ZHbC5xTmVzbmEiNVVraiYzJkVNKnFmaFVZcnJBcF91bTI9S3JvWEMjMkZFIklWWWhzMU1gXXFHdUFsaDU7MlVYOEBWVzk0JGYlK1cjKks7ImIuW2tMZVFyWjZmL2hWWjhNLSNDMW1hJyxGanFMJElZLEEpOWxUSDtUYSljKz5KNGtza10qNyJvViU/XWRuaiZxME9aKm1uWWZmaEwtQDpJcCdCLSdoJyhwT1BGNGopM0s6Qic9S0NRI2puUlBmRkkzdXQ5WzNLK091WXI9QWRkQTMpaUUtXU0lSDU0JkZxUDJrIWpGOVs+UHFEPEopdUBGQyRDZSZyWFwzXUFXRVI3MTY7Z2YvbWBKNnFHU1olcyZZQVwmV0hRYSwtPmNncjg6dW9FTjJrXlxiQks7Im9wcl1ySTFjZTlHcE02P0w5RnU8Sz1RVDlFOlhFW1U0W2RoSElxQ0ssV1VaTGNEQE9oZGYtWmxqRTM+cDpHTDBZTygqSVQsaGRjdEwpZlQjWiNgYC5YQidSSzRJM2VJP2okZm1VPiFgKmlwbFgqYEJpTzZjODUsNlVIWG0uVjtjdEFwXCZdRVNuWz43P102XlE3VDg0bklDWm5iMWJCP0orPylhJVRWRSpBRWdKViwsPEVSQHMjR1BCWGlGaS5SM3A7OFdmI2YoKVQuaTVRaVYvQE43T25jREZtT0RFTDhJSFIpcmVsU0xNclZbSE9BO2RDWVBGUDIvQE8qJTxELlM8NkxuLV8/OFRVbk9jTldgKTU0OVUpNz4uV0I/X05haiwjYyUtLHQ6WF9HXDNjKGIvP0ZzUk1hPF1sODNyJ0A2R3EnJFojJ1goIk1na2FNMltFMX4+ZW5kc3RyZWFtCmVuZG9iagp4cmVmCjAgOQowMDAwMDAwMDAwIDY1NTM1IGYgCjAwMDAwMDAwNjEgMDAwMDAgbiAKMDAwMDAwMDEwMiAwMDAwMCBuIAowMDAwMDAwMjA5IDAwMDAwIG4gCjAwMDAwMDAzMjEgMDAwMDAgbiAKMDAwMDAwMDUxNCAwMDAwMCBuIAowMDAwMDAwNTgyIDAwMDAwIG4gCjAwMDAwMDA4ODYgMDAwMDAgbiAKMDAwMDAwMDk0NSAwMDAwMCBuIAp0cmFpbGVyCjw8Ci9JRCAKWzwyM2JhOTEyNDUyZjgzNDNhMmU0YjJlOWVhN2NiMzRjND48MjNiYTkxMjQ1MmY4MzQzYTJlNGIyZTllYTdjYjM0YzQ+XQolIFJlcG9ydExhYiBnZW5lcmF0ZWQgUERGIGRvY3VtZW50IC0tIGRpZ2VzdCAob3BlbnNvdXJjZSkKCi9JbmZvIDYgMCBSCi9Sb290IDUgMCBSCi9TaXplIDkKPj4Kc3RhcnR4cmVmCjM3NTIKJSVFT0YK",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for per-cell quality-control reports from the WTC11-EC 300-gene CC-Perturb-seq D0-D3 analysis files.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 09:26:03,969:iu_debug:	Success.
2026-09-30 09:26:03,969:iu_debug:	Object posted with identifier: d860bb2b-8907-4bf2-88dd-5686e56b102f


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/POST_document_annotated_h5ad.json
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/document_annotated_h5ad_submit_20260930_092604.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:26:04,510:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:26:04,511:iu_debug:	submission=True: In submission mode.
2026-09-30 09:26:05,194:iu_debug:	
IN post().
2026-09-30 09:26:05,221:iu_debug:	Validating the payload against the schema
2026-09-30 09:26:05,239:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_300gene_annotated_cell_gene_guide_h5ad_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgMTAgMCBSIC9NZWRpYUJveCBbIDAgMCA3OTIgNjEyIF0gL1BhcmVudCA5IDAgUiAvUmVzb3VyY2VzIDw8Ci9Gb250IDEgMCBSIC9Qcm9jU2V0IFsgL1BERiAvVGV4dCAvSW1hZ2VCIC9JbWFnZUMgL0ltYWdlSSBdCj4+IC9Sb3RhdGUgMCAvVHJhbnMgPDwKCj4+IAogIC9UeXBlIC9QYWdlCj4+CmVuZG9iago1IDAgb2JqCjw8Ci9Db250ZW50cyAxMSAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDkgMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjYgMCBvYmoKPDwKL0NvbnRlbnRzIDEyIDAgUiAvTWVkaWFCb3ggWyAwIDAgNzkyIDYxMiBdIC9QYXJlbnQgOSAwIFIgL1Jlc291cmNlcyA8PAovRm9udCAxIDAgUiAvUHJvY1NldCBbIC9QREYgL1RleHQgL0ltYWdlQiAvSW1hZ2VDIC9JbWFnZUkgXQo+PiAvUm90YXRlIDAgL1RyYW5zIDw8Cgo+PiAKICAvVHlwZSAvUGFnZQo+PgplbmRvYmoKNyAwIG9iago8PAovUGFnZU1vZGUgL1VzZU5vbmUgL1BhZ2VzIDkgMCBSIC9UeXBlIC9DYXRhbG9nCj4+CmVuZG9iago4IDAgb2JqCjw8Ci9BdXRob3IgKEVuZ3JlaXR6IExhYikgL0NyZWF0aW9uRGF0ZSAoRDoyMDI2MDkzMDA5MjQxMy0wNycwMCcpIC9DcmVhdG9yIChcKHVuc3BlY2lmaWVkXCkpIC9LZXl3b3JkcyAoKSAvTW9kRGF0ZSAoRDoyMDI2MDkzMDA5MjQxMy0wNycwMCcpIC9Qcm9kdWNlciAoUmVwb3J0TGFiIFBERiBMaWJyYXJ5IC0gXChvcGVuc291cmNlXCkpIAogIC9TdWJqZWN0IChJR1ZGIGZpbGUgZm9ybWF0IHNwZWNpZmljYXRpb24pIC9UaXRsZSAoQW5ub3RhdGVkIGNlbGwtYnktZ2VuZS1hbmQtZ3VpZGUgQW5uRGF0YSBtYXRyaXgpIC9UcmFwcGVkIC9GYWxzZQo+PgplbmRvYmoKOSAwIG9iago8PAovQ291bnQgMyAvS2lkcyBbIDQgMCBSIDUgMCBSIDYgMCBSIF0gL1R5cGUgL1BhZ2VzCj4+CmVuZG9iagoxMCAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAyOTMxCj4+CnN0cmVhbQpHYXUwRj5CQVEtJnE5U1lmWElgS1FtVFdHb0k3MVhcIzhJSVhZbWU5MWJfQT9iYFBHTVV1Sj8lKzVbNztNYVdzOV1oWGBxIlUhOGBRWicuRSFiPl4hSiM3XDgpaDtqPUUsai1RRUw4Pi9cWCFXMDJ0Mi06cEFXOEM9bGU8XD1KPWRGYU1ncS4lKFFhMWhqaVkoJXFjYj4qI21cPCFaczwzR2BoKz5KJCVdJ1g5Sy5bRFYyZzBEQFBtQkphR202R3JCaCItcyNtSFhCJidzMVdUJlJRdTFkI2xZS1FIQyVXIkwmVS8xKGhxZTNOJUVpamFeYlJRJFU0UkVUYDlCJTo3QCEjZ1FTPGFrWmBnMS8iMnF0XTppcFhWcUdzWzE/akJMLUtKcDpIODwjYms2bC9FajBbO2QmUSJAWk1LXWlhPEJCXSlTW19cNkUoWTlYLy4nO1c6KVUyTy1aSSgkWGJTS14tPEwyYTNXMC4+JmFhL0FhW20rNEc9YFEuI2BWRjg4L2VZcilZIWFROW00LVBGTkknYlgtcWAva0ZhNmZhTWVrK0lfL2UkaXInXUJdbW1MR09YMzJWYGdoVWJaY3NlUXFDM0BxX2tEYV1FaEBrZjs/MmJJIyxYbVAhQCRWOSVfUVU6PWBtLnB1Sm4qUVg1Lko2ZkoiJS8mWSFCTVhJOnEpbVppUkxMaElGQyxBS1lGZiJnNkk8KiNsRUxbNT1nP1JjWzJqPzgvR1FESTBENFFOIz5EUXQsWUtacEpeNUFlTGBYVm4lakMtRHFBUCJVMFlDPSw4KyQyXDJSMnFbWmRAN00ycTw3QTkqOj5eZVMhYWZlK2szJF51UScraz1uKjc2aSpfTFVqQT4nNU0+Ok1oZC5aO1tmPzc+TjlvLHFIUzg9SSNvP0ozZ21lWkU0Q0pvZio5aEIhbnRwYSpdMXNQNyE+WjdEP19IaCUpc1dpdDM9MUZzJzBGaWB1bnFsWmJLSnFHO15vYlpcP1lbaGVJSl0zUVpaU1A0cC4pLCQvMjFhdDNPP0BOYW9UKVxFVmhwLFBbcm5AbEtrUChIKUlsKVklYmtKaV5ddCV0WU5wUC0tclVsNXA7N2pfclU8biVsJDViPWtpTiY0O2Y4RlJOWSEiMXMrcjBFKUJnQkI4T2BnS1xtVGdILElzaj1sQGMxKDIlV1MoPzhRSi1iaF4vT3BCaTttc1xjZDUoOmpGIzMtLih1K25Uc2krSUNQZFwsbGArRU5yPldBLlInYCxhSkIzUjcmLm8uT2YhL2pNK3BtIl5zJystZGVDQEEpaVErVSghXVdcYS9pLF8/JnAzMChQL2g5LWcoVT9CYkc8Rlk0PyVpI1tDbVIyImNKTDovLjlxVipkaWNEaThdXk1TMUg+WyRgYF9eSGwkaUZlOVhqVC0ybSpFLHJzLT8mT0hDYXBYLz4kdEssIltSZGhxKltiTmVBaC8sSU11UEhecTFBaSwqMzBGL1ByYVAqIkM+V0VlUD5aamZhKVwlWElaamxRPitwRDJFNUg3YURtXzlPQ2xcM0JMU10hb21VVjthL0VLS05PdGFEUDE1SltUW0xpVTZrU0hwOVYpMlRKb0xSJnJxbytFSz1PNltYdGpfPihPcz85LiYiOWk+UCxgTysrZlJLVz8vckJnZ3QnZyo+T1kxR2UnJD08TUc+b01hXGcxVlVOY10jYy9BVy0yK1BYNz9wSV1ANFtkNFBUby5xTSFvSi8sKlRIK1pQKWhGXXVPSzo+VnJWTVFaMHE1TFZINytpUUFGYCM5Qy4vJ1lQcj5kQTxBR0hQbypXWjJZYGJwYEA9TDc5TWRlVSFEND03JnNjNWxWQEBEJyR1dStIMSQxV29Cc2cwbDkxS2hFSDBkWEw0Uy4pV2dfNDEyNGBVSERYJ0JQQ3MtOT87ZUB1JXAvLkgxbjVNNENmXz1KO1xGJT1bZHVXOD1WQWwiITA1YlpzPz5NTzU1bUpFUD4nUVVUPzc0KGRgNDtFcWhGWmF1MXRHNi5OQVJXSCdJUSphNCs7TFc2Py5CWyMwTFlwbFtHU05Lai9OLF9tRyJXISI9XmUsM3FSLTMwTilfUUs8UVBWMkBSMlxyaSVSY1oiPForI2hkVlVQUiItJytzY0ZyXHBFSisjRnBtVmV0P0hFdCNWczk2M1JVSTIvZV9CZTNVUUZFRWUsTFFUTDQpOUEyPkw3dTAkYU8lXWRRIjJsTUYycy9vZlFzQ1xRN2RlPGZDP29fVUlFL3QvZlVyZWtraUlQSFRXM1xXNW4iaitYQ2JHU0MiX2k2VixYZT1sOzw2K1gzUFBvPmJpdSpMNHNmTyNUJ2QsMS8sVF8iPm5OJ2NSWkwpWD9qTEpmdDgyLCsnSGxjT2Jbcmo5WyIpJSVTYyJpM1ZwdCtcNzhGS2NtSz1cSmlwNEk4KidsNS9dNDArQUR1T1xANDJHKzpXVGcjWDFXPThLbWpbSHVmUSsqRm08RlVeMiZOOFdZIURgYCg+SSM+dWhXMCFnPkg7NV80VVEzIylHOko2NWFKPWhyJzUlM2VQMD9JNUlYL2hpTVxHNFhHcklWXDFXIzonRmdUInI8XiY9Wy5wXCdNV2g6alRPRmAxSkA0N1giZGIjPTojUzMwU11XcC9QcGN0TlxgV3VeK246XzJiajQ1YihpNSU7bUBjO29eOTlDXUBWSSdkL0FES0FWTWdwNV9PWjk3I1xlMSQ+O0RpRjs0MlVBVjMoZCRjcylWVWcqViY5JGRcYj1zOCwoQUtxZk1qai1pQCg7XDRxUztJKmxJI2RvaEJcNERwVCkmZUNXT1coXFoqNjY0PFkkKXBsNXAwJi5nUF84aHUtQ09WWGZia0lfUWwlXGM3JC9ORz1wIXFvTUxBbFYsXlBNTFNERFE9ZiplZFdbYFcmZFYrcUg9ODM6QGs7b0txRUY/LFVNVSVTcGlrJGc+PTgqQyZzZzQhKWhRSFMrczVhPD9RKlA2UDxLTFVmSiMoY19dSEhqOCE/ZnIvTTgtQCEpSDdIS0swZFtDSVxwblU1YG8mM0hMOkNMZFxFOmxyS0klKGdrR3EzSzlLSTE2Pmc5Oj1KXzU6QXNdLV1GTWQsUlZqTUZYIz9TTVRrLl9wdClhXTFrPitDPS9pOTwnY3NEKVFdQnBBaCVOUSc8J3VmcEImUTxkVnVbODxLZkJxZF8iO1EsWWNAIW05KTU6SzpTQj1NKXEsb0E+RS4hMVFzWk5DYkZVIS4tR2dPTTZpP2xBcEU9Yj5nZCI3XXBiQGI1dTFtT3A2KGkzOEspSzM9KDQ9VDclL2VmO3JzSDkpWFBCSzhIVXAjODosRidhPEQpZU9NcHMsXDYuLUBkVWArVylKUS8hPUpUQXA4WWokN0JJJl1FSGEjT0wuQVROQCRsJSxkVmJwIiYySCoubm8/ZylQJTxlNihWRiNuPjRBcEdSayVMXXVXMDEuZU40ZCdlIl86dS1XTTxiUylRQkxKUT1AMSstaERGZ2kvYFVWaVghTEYlPFdrRVRBdUlGcUAwdXQ1OG4/aFI8KClmTmJlM0FwVSVpKVArIk5QbHBVUFxpIlBnWFRCV187IlQnbjNLU0ZZOEBtJEcuXytzJ08jIVFKWElAYzljNitvcV9GaFgjV1JFcCFfMiwtQExNbCxNUCcnZHNRcGFKNE5fQCtESF1TX3FUKmpYPUNbQ0dVZGo3XjkwSEU9ZkglIWhHTGZcS3VZbGBbbk81SWVXJ2duclxFJXJPOjA4YVpVM0R0UTU8WF1iPmo/WVltdHJFI0chP2lQZ2xnIjQhLFEkclVYPWtWRUFRJ0hkcjQnQVwxajEiLmw8QWNoYG5dOCcvS3JKczQsRUNnOFFjP1ZlTk9XOzk+SDp1UUxxW01VWmk6JEJTMk5Bfj5lbmRzdHJlYW0KZW5kb2JqCjExIDAgb2JqCjw8Ci9GaWx0ZXIgWyAvQVNDSUk4NURlY29kZSAvRmxhdGVEZWNvZGUgXSAvTGVuZ3RoIDIzNTcKPj4Kc3RyZWFtCkdhdTEzPyNTSVUnU2MpSi5zVF9eVzYubmFnaHBlYEQmcCdmXzliOVtFVEMpUUQwYFJyWG83bystJlMjOjg7RFdwLCJUVG1nITw9VU0oImc5LFQoSTc3Ly0ucmE0UUhebS5GVCFLT2hQMi0/USE4bTxlUWhwPyhoNSxTcm5sRi1ibjBhP21AIzlxRyMmcmQ4cT01M2gkcDRdQVknODxtXjdnVGVzOjhhPzthQUtQaDEzNCVEcUBoXDcmPiVjVikvMj1kaCNeMVUiOkNzW00qRV02Y2MuVktFXW0tbDUzKixfSilnZj8yQSQ9UDJuLEBYVjJZRkxmUyolNGJTXzw3KFhtOF1wP0NnQ29ycl1uM19KV2A3WVJyVGc9Vzxva1RSWCRhLk9JWkhybToqKW44S09MQW42Nz1sdTxTN1InYzdcaXE8QVo/LWZuQkhHIi88Q0lISSxrM1NBW20rRCRxQV1xO2JfZS0nRENlZmpFOk0lUC8sM1B0ZnAxSCJDKTRtSy9eO1NMZmRDa1IsYiFobk5nM25qVlQnanNpLV9HUlNpcGlwLjF1WywjPmZCWHBVKDxNcWtzTGstKyVmaSc2Y2hiPU1JNVZjU3Ndc0tSOiMyR1ZNZFg2XF5INzRNYSFyTGRZQU5wNTIqajYzaEZEL0k+SU1mKzJUVComR2YhQGgrNG9nJzNHI1koM1NcPmc7ayVAUXRzYUI3Mj5BcycyU0BWR0N0NmxWPFAuLl1ZITt0SiZWI0dUPlwpPUdMZDwpXig3IT44YypBLkJZLCQpPWtbM1VsbyEndD9Ncj0+JWFEYEFULWtHQG5ZSU5ITSpVPT88YVQ0O1MqZjZSNUR0anNTbic4THBaJFtFUzQrcm4tT048IlQ+KUdJZkJfaGhGSUVaYy1HdG5oYTNVZVIyO2UuL1dAaC1cKnJkRkNfNTdJSE5iaWI6Mmc9alZtVVVxKV9mZzBvbXI3Xy02amxGUyJNKm4uWCpTSD4/Tzs/LT8wKkxnUTQ3RnAqLm8rJkRJJDM/MFgvb2c1RSZoOGNmbVJHYWh1SStua3AsaypsXkByIUsjNltUblQzNEBHSzVZJiMwQEhDZ1ptKDBvMGk8PTdcciFxODl1bFY+T2NXOmsnPEZebjhkL3Bxc1ctQzpjJi45c15aSUBYKyM0KDtYa24+L11WKigxajpVJSolKl0mKXBGUmJEMGg+U11cUkpqZ29gbVRCRVBOO18xVydQWipPUWZmI21dOS1YV086cygkLUNHOVkwb1w+czxnSF82MU5gXDdqXWhPST9maTQtP004PG5xV2VgSVkhUCshWWZpXWxsSTtBTFApVyUqZ0tOMF4lRWwnPEklNkknR0dfIj0sY0sjNWtJVyY2VTVjbi5dJ1UyWFcpbE9MKCNBP2lsYl9cVVFwS0QuLDxcLihtJiFLYT1gWFVGVF8jcSc1MzBQVkxxcU4jLlciZl0pZS1aYy4lRmc5TyIxNSdmJy4mYjtyQyFYP2MqLiUrVGlkUFRRQHFNQGZOXGA5IzctVFdrIjZKN1VATUZEMjNeY2M3VCdVXis/KitEInFRZUUsIysxUm1QQUtnTGYsY3Q7SlNgRTJjKHEtVmo1Q0tjTF9EJ0lPbjpcayIlKVlKb0Y6NShrTEpIQ0YwO3Q6bWY6ayRBJHVWbyFePDAoZ14pajZtK0ViK3BKUj5SKVlRaEpEQScqZVZPdXJUZEApcCJFRlUyRCpgSnU1bmRuLks1cidbZ05xaSdJY2g4cC5yI0ByalBESyU9bEI0LzRmODJCMzBQO0NIYEZoXGpedG5hMmIzLVosYWI3JyJXcmdGMzBKbGFWRU4nbTFFKTQhP1ouaEtyRCJcMT9eTmY2TERWc2wwRWxldDksKk5OMWAxRkBlYkBWLjRCQUFDLWhDaT87aSwnNkMuL0JnTDhrQic5clxtclM9J2ktU1omJiMhSUwpPkw4Kjc+NV0zNWVbK0RcLTRPK1dHakEmXCIldG1cKFFEKDdCKjA7LD4nJ3FLWDE1ZHNDIVAmbC9vKWtbYk0uSilNIlVGK0NXV1taQ1gjK2ZDNDxTakJFZGYtUS5eOUNOWGxQI2Y6bUdGbEJwQ2EpWmo1WWU8NmMjLToqbDQiVzpfLWwpbGFQPDI1c0pDUmVDLkVNOT5YRmQmRnI9J1hiNEtFXk5pWV1yImVlOENKVz00JkoiJiYzPyhZTDhHYDI1XTNNbU4jblhUOCIwLVpna2oyTEpeQ0M3MSxRRF8iKi1iL0BXWWRiKmJYJlZANCo4XnFxWzFLK3VIXldLNCtqPkBXOiUxKDIoTDpAN00oOlxdZltwQlBOZTtfWmo7X2kzWW9kQjxIa3EtUGlacipcOVdWJFRvN2FjbHJaKXE9QzQsKCw+MG8sLjRxNTc7Mjg8XD4mczs9MHJNJjcvRGlQSnVwLW9mcCZEVF0zVGhHZV1nOzlOPVtYMmUiPXBCVyJPYDwvQ0spIUdMUmZXX1ltRlRlZ3A9aUE0cUFoQmVBMVNVRy1rOG9gJF1UWUIwJ0AnXiw2WWZNQ1ZoREtSKC9tUUgndG5XMVc4KyYwWTZyY0A6MWdAJFtSOGBqcGQ/dVBoZV45MTkmYmhla2ZdL1BsO21TMjNnUjprLi0hbERENC9XIzMkbydHIzlgcXJCP1Y2PiFla29NZDNbZiU5aVdyQ0QwYSdLVGNPcSZHQlZQaE1QLGdwT01EOC9kb1NgRltYQ1xxdSYuJCZLOy1pOzliRGxEO00kMiQxMW8lIS8/OWpRPCpXK2liN09LNWZBP2UqJStabjVXYl9xcys1WT1iRTY3cXBiVm5PbzwuZjJsaDQ9MGoncFtYOXFRZFZSKi5jJmNPQV44Y29FZ1tcQS9lKVZwOnFfV2FuWlpqbV9YMSUsanEyZjJtXjcjT0dbY09RJThAJ21UPEM8M2xxNVgoJEVnLm51Ul5mMEgpcighYmdbYG1COEAnbVQ5bnRSW25sN0tBZnA4aSwraklIQzxWZU80OGVYQl8+Z2VLUVhMTWcmNldBTStiKGNTVVAlMDRgMVx0RExeUFxmSGs3P0hgSixKXy5jOGs9PTMuRVJTQixTSVNEZ2s/L05xKDxHay1zdHEocilTVkRUZGgnTXIiOkFHKDxeJX4+ZW5kc3RyZWFtCmVuZG9iagoxMiAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAyNTEyCj4+CnN0cmVhbQpHYXUwRmdRTFUqJlVnOlRedU5GN19eLlIiUVpEMXJQUCRuMW0jKDIkYiM+PD4mLihwNVBHJC5DbSlpXC1ETSJaUFIxTUMiVzI9LzIhaV41WlZdNEhEZlRgbzg2YkMqKTYtLnJUYW5tR1gtbVImIjZASG0lY1RkVF1TZmBFMkUkPzIqREZDIzlqIjBvXzVFQCNOUEApKFNKQGlOWC4sbHJ1MFtDJS07TlIhXDo9LWJ1YWVzbjgzKDJQOEtRWzI2ZUZ1PXNhMkhwaDlNZVhTbkk7aWFYb2tQbScxQ29sK1JkJENxYXRKRis2NXFzLScuLEJuLE11SUckLUg+WSxodVhPU1VBaENtUUhCYVwpNEJvKStuPEl0TThKTG5QVWRGYSpyYFplLnExPy81QWRYRWg9SDEwQkRVSk8xYWZBWCFGLVg2NjJZQCFoLC9jPUApO0hdL1IpXUIjK05VXVFhbyE+LC9IMCNNLnJRXnAqbCFZdUNyTmErNUxwY2k0Lys1Vk5HaG5XYGFxYiJMRFw7QTgwYyFTP1pYcXJxNkEoVGlyaVItcnA6JTA4ZD9BRCVHUyxUKHMvLF9kMDluLSI+UTJcU0dccU5YXWwmWVY6NC8sZiYsSjw3THJecFBtMmRPNytsXFhyOFxEOyZRV008KnFETzg7b1IpNmpqaWMqbVNXamhTNV9hU20lP1ZGPCNeYj9CZCpsRmUtMF5vIUkoPDk8RWVQJWE7QVpjPDs2Wig2Tix1PHVuMTJkUDNLUCw7Y0tYVVtVWSM4Qz9UZmsmTDRCcSRWMz9ANVksU3RUPyxrLShtWUtgdHQzKGdPRWkpMT5HKVljJyEvcUxQOyJzO2ZDISV0YSdFU0h0S18jJCNfalRZazBwOFlWOzE7LUloTlwjXFJMYUUhUDAxP1FMJGY4KWpwUyxIX15YQUxTZ0ZlPUY5OiRWNz5bZHIuUFkwJm5hXjI4IlBORy1cJVVoKEwzNWdpaUYsNEhFWEVUSUBcOTpGO3BVM1RNUEtTQTpATCZFOEZYJ2trbSMzVlQzNV5jI1YpLkpSPSQpJVUpdSpJJWVeRDE5X3BvJmRNLjUtV0pTNEosKk9BQCw8SmYxZjJzZT5IWl1uOFNpcyVVI0QjXmZxZUVzXWNUUUhFNSIlaGZJPHFaKT0xWz5jPzJMKjNkTz5LVlM3QjpYUyRtZENHWEsnRktVMz0mV0tZVVdYOSVGczVFTDBZczlRZiFKSSNYIj5QYyxlNi5MIiFGZz5QZStVZTwwbTt0NzQ+P19jVGgsc11pMC5ZZltXQlNlZkpkQD9QTCdqWDJXWlFyblFpWEAtNSh1OyloPWxOMVVvRzFFJVooWT8oUC03XSVOWTU+I21vKnBPRUAzNik6Qz1kbVRxVCYiUC03YCZMKFtNSXBTZmQyXC5dVlFQKT1nYVxsYi0hOy9xPSc8VjVgWT5xQnRDLkwmTnFFXk91UyxbWz4rJG1ocHBtU21LWUUkbSxpOS5HN2JsPSo4dDdgPzY7PDA/JTpdJUd0OjwvZCxCQT5cWCQ4JXFJMykrbWxVXWkwJ1lfayFgPjRwXm00XHUrYiU4JXFJMykrb0Y4bVNkQlc/blhANE1eT2tEOmczNGE9IiQkcSNSS0M7bEE/Vks3KTVXZDk0WHRORzE4I0ZTMDgiPElcMmQ0IzxOTjpsK0BXRyRBb2tzJDM9Y0YhUHVOSDI2aGVZYjw3XzRHSCxmR0woW3Q5XDNkWUI2M2VUaWgqajdhPT0oPFBAOVpjKmhhTWpHX0k6S3I2RyIsUG1WJDYjcEVJaS1qOCt0MSVCTk86L0BQQCQ5LEErXTxvOWkuNEA/SkMvW3FBNk1CZWdbXycvUXNUQiFOMF5QaGlmOyJtbDA5bW5NIkFcSmxUMlhpWj4kZVA5cWxHOHAmak0lXl1eOzc6cTJFQzV1bkdFJHQkMmEsamxqXEpkTVxDZnNRZGJVUkghT2Q/RFRYRT5EPzFfYzQ0O0s2O2lwNzZYKFVASlQnPitcaClCYUBeMmYiXm5hPyI9NDU6VzVfI1ZnazdnJ1RLTFw3IzFuVGA9TUhtZTNfcCNHUzJjYjBwcis5LSVMZjJjO25PZTBQPjMjYD1YaWI0aGtKbWdQQ2NFaG1wJjRiY14ybUQ3U1M3Tm9hXTNxQlkpYTFOOTVFRERUZDVXJ1NDY0U+R0s3T0olaiIpOylFJj5sOiF1bTk/PVozUmokVDomT2Q1K3VoT010azhCVjo3LjldODVOU3BDKXJfcDluNTVpUS83MnAtP0xTb0hkaUFJUStaIkhAbFw7PCs4LkI9c3NLYF1aaTxHWFksXUhdRlpyYXRNQUchLjQ+VExwWVAiMlZHbF5mWUBwZmVONFYuXnBvKUhMXSFlZjZXQDMwLCN0LGlHZWc4R2g+Uk46J240Yzg3ayVuVlJZW0s6Ok44b1FqUD48Ilc1JXRMOCthUlpLLV5LbHBdTVw2OSVBX2thYmw8UVcvZFprUD00IXBZMEhQTyZMbmY9Vy8yWHM8M2VLOFYmdHJ1LWxQQm5BVERcQzlIYV9VTTlUZDZEWDVScFs9SlllX2ItKSU9X0NKJiwvNDprQ0FqKyJdKjMnaD9HVHBMYiduOWpQK2ElSzQpbjJNbnU4cVpcaCs0Lic/XjRWWEQoOU5eUFxmZkZbOF89KyRpMlxHNShJSk1GYGpMXHVzOihAc1JdaWBjX20mKzN0cmZVXyQ3Y0kzZitxJzxRUlxJYVBhJGRDTCdOKiYyLTJFc2cuWyYuImtgUlZcPixLaFInXXAhUCU9a0hGcFMsIVAnSlc3UCZ1ajVcTUxnPCFJNywkbE8vQCROVSYiJUxDWkxeZDYwIz9jayM4YE1bZVFdcjIkaElMPSkzb0hjRU9pKTU6aFI0PlM5Uyk0aCUlKGo+VTNsb25vMHFebE85M0M4QyU1IiZoSWM6bFJQaWdLQDpUOWA6R2stJTcoMDxlK1NPMFtTYHNKayZxX28pXlNVLT4+V2FSQS4rUkBmOj09cjglVCo4IVtqbHNGXFo8az5iQ2lbS1JyL1NAPyE8ckZoQykwPzglMTs7UkpLWUVEbVssIkQ9ZzxQcGYiWSg9RWRmMFQwR0xPSmFWKjRIbjpKaWdGLUY3ayE3QW9lSFVDIWFDITFkZipjbS5EbVo/KjFfIiUvRlAoN1RqO2Zza1ArVD5BckllNTo6azVvcFc2MDplMVhBT1teNnRBJmpENGBhM3VQWiU7R2o6RCsjO1lXPUA5azFBS1E4J1I+QTs7ZytzVUlWSGRyYD5PXUVzOSdXOWdQR1gmaDVBWT0xOTZvN0lTZWRgP2A5aE5XLU0kWGovNmY0In4+ZW5kc3RyZWFtCmVuZG9iagp4cmVmCjAgMTMKMDAwMDAwMDAwMCA2NTUzNSBmIAowMDAwMDAwMDYxIDAwMDAwIG4gCjAwMDAwMDAxMDIgMDAwMDAgbiAKMDAwMDAwMDIwOSAwMDAwMCBuIAowMDAwMDAwMzIxIDAwMDAwIG4gCjAwMDAwMDA1MTUgMDAwMDAgbiAKMDAwMDAwMDcwOSAwMDAwMCBuIAowMDAwMDAwOTAzIDAwMDAwIG4gCjAwMDAwMDA5NzEgMDAwMDAgbiAKMDAwMDAwMTI5OSAwMDAwMCBuIAowMDAwMDAxMzcwIDAwMDAwIG4gCjAwMDAwMDQzOTMgMDAwMDAgbiAKMDAwMDAwNjg0MiAwMDAwMCBuIAp0cmFpbGVyCjw8Ci9JRCAKWzw3MTE2MThjOTI3ZTRlZGVjYmNhM2VkNDQ5Mjk3NzU4MD48NzExNjE4YzkyN2U0ZWRlY2JjYTNlZDQ0OTI5Nzc1ODA+XQolIFJlcG9ydExhYiBnZW5lcmF0ZWQgUERGIGRvY3VtZW50IC0tIGRpZ2VzdCAob3BlbnNvdXJjZSkKCi9JbmZvIDggMCBSCi9Sb290IDcgMCBSCi9TaXplIDEzCj4+CnN0YXJ0eHJlZgo5NDQ2CiUlRU9GCg==",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for day-specific annotated cell-by-gene-and-guide AnnData matrices from the WTC11-EC 300-gene CC-Perturb-seq analysis.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 09:26:05,582:iu_debug:	Success.
2026-09-30 09:26:05,583:iu_debug:	Object posted with identifier: e90fd5fc-23b7-4bb7-9444-64db07abc82a

            format_role                                                                         document_alias                                      document_id                                                            portal_url                                                           attachment_download
  guide_quantifications          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1 /documents/4ecb21f4-7a3f-49d3-bdbc-3161d349311d/ https://data.igvf.org/documents/4ecb21f4-7a3f-49d3-bdbc-3161d349311d/          WTC11_EC_300gene_guide_quantifications_file_format_specification.pdf
per_cell_quality_report        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1 /documents/d860bb2b-8907-4bf2-88dd-5686e56b102f/ https://data.igvf.org/documents/d860bb2b-8907-4bf2-88dd-5686e56b102f/        WTC11_EC_300gene_per_cell_quality_report_file_format_specification.pdf
         annotated_h5ad jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1 /documents/e90fd5fc-23b7-4bb7-9444-64db07abc82a/ https://data.igvf.org/documents/e90fd5fc-23b7-4bb7-9444-64db07abc82a/ WTC11_EC_300gene_annotated_cell_gene_guide_h5ad_file_format_specification.pdf

Saved: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/document_portal_records.csv

```

## Markdown cell 23

## 9. Build and apply the 12 file patches

The patch adds the matching Document to `file_format_specifications` on each existing file. Array properties are extended by `iu_register.py --patch`; unrelated existing Document links are retained.

## Code cell 24

```python
# All three Documents must exist before patch planning.
document_record_by_role: dict[str, dict] = {}
for role_key, spec in DOCUMENT_SPECS.items():
    record = get_with_retry(spec["document_alias"], attempts=5)
    if record is None:
        raise RuntimeError(f"Required Document is missing: {spec['document_alias']}")
    document_record_by_role[role_key] = record

patch_plan_rows = []
for record in target_files.to_dict("records"):
    role_key = record["format_role"]
    file_record = get_with_retry(record["accession"], attempts=5)
    if file_record is None:
        action = "BLOCK_FILE_NOT_FOUND"
        error = "File accession not found"
    else:
        expected_doc = portal_record_id(document_record_by_role[role_key])
        existing_docs = file_record.get("file_format_specifications", [])
        action = (
            "SKIP_ALREADY_LINKED"
            if link_set_contains(existing_docs, expected_doc)
            else "PATCH_FILE"
        )
        error = ""

    patch_plan_rows.append(
        {
            **record,
            "document_alias": DOCUMENT_SPECS[role_key]["document_alias"],
            "document_id": portal_record_id(document_record_by_role[role_key]),
            "action": action,
            "portal_url": portal_record_url(file_record),
            "existing_file_format_specifications": ",".join(
                sorted(direct_link_tokens((file_record or {}).get("file_format_specifications", [])))
            ),
            "error": error,
        }
    )

patch_plan = pd.DataFrame(patch_plan_rows)
PATCH_PLAN_PATH = WORK_DIR / "file_patch_plan.tsv"
patch_plan.to_csv(PATCH_PLAN_PATH, sep="\t", index=False, lineterminator="\n")
print(patch_plan[[
    "day", "format_role", "accession", "action",
    "document_alias", "existing_file_format_specifications"
]].to_string(index=False))

if patch_plan["action"].str.startswith("BLOCK_").any():
    raise RuntimeError(f"Patch plan is blocked. Review {PATCH_PLAN_PATH}")

patch_operation_rows = []
for plan in patch_plan.to_dict("records"):
    if plan["action"] == "PATCH_FILE":
        tsv_path = ROW_DIR / f"PATCH_{plan['profile']}_{plan['accession']}.tsv"
        pd.DataFrame([
            {
                "record_id": plan["accession"],
                "file_format_specifications": plan["document_id"],
            }
        ]).to_csv(tsv_path, sep="\t", index=False, lineterminator="\n")
    else:
        tsv_path = None
    patch_operation_rows.append(
        {
            "day": plan["day"],
            "format_role": plan["format_role"],
            "profile": plan["profile"],
            "accession": plan["accession"],
            "action": plan["action"],
            "input_tsv": str(tsv_path) if tsv_path else "",
        }
    )

patch_operations = pd.DataFrame(patch_operation_rows)
PATCH_OPERATIONS_PATH = WORK_DIR / "file_patch_operations.tsv"
patch_operations.to_csv(
    PATCH_OPERATIONS_PATH, sep="\t", index=False, lineterminator="\n"
)
print("\nPatch operations:")
print(patch_operations.to_string(index=False))
```

### Text output (stream)

```text
day             format_role      accession     action                                                                         document_alias existing_file_format_specifications
 D0          annotated_h5ad IGVFFI7756RKHT PATCH_FILE jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1                                    
 D1          annotated_h5ad IGVFFI4783YGLZ PATCH_FILE jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1                                    
 D2          annotated_h5ad IGVFFI5165BSGC PATCH_FILE jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1                                    
 D3          annotated_h5ad IGVFFI9086SWZC PATCH_FILE jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1                                    
 D0   guide_quantifications IGVFFI1608GIRE PATCH_FILE          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1                                    
 D1   guide_quantifications IGVFFI6431OXJN PATCH_FILE          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1                                    
 D2   guide_quantifications IGVFFI1924OSEM PATCH_FILE          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1                                    
 D3   guide_quantifications IGVFFI7077LERW PATCH_FILE          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1                                    
 D0 per_cell_quality_report IGVFFI3905TBLN PATCH_FILE        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1                                    
 D1 per_cell_quality_report IGVFFI0627NEZT PATCH_FILE        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1                                    
 D2 per_cell_quality_report IGVFFI1451CYGW PATCH_FILE        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1                                    
 D3 per_cell_quality_report IGVFFI8282ATZQ PATCH_FILE        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1                                    

Patch operations:
day             format_role      profile      accession     action                                                                                                                                                                 input_tsv
 D0          annotated_h5ad  matrix_file IGVFFI7756RKHT PATCH_FILE  /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI7756RKHT.tsv
 D1          annotated_h5ad  matrix_file IGVFFI4783YGLZ PATCH_FILE  /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI4783YGLZ.tsv
 D2          annotated_h5ad  matrix_file IGVFFI5165BSGC PATCH_FILE  /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI5165BSGC.tsv
 D3          annotated_h5ad  matrix_file IGVFFI9086SWZC PATCH_FILE  /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI9086SWZC.tsv
 D0   guide_quantifications tabular_file IGVFFI1608GIRE PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1608GIRE.tsv
 D1   guide_quantifications tabular_file IGVFFI6431OXJN PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI6431OXJN.tsv
 D2   guide_quantifications tabular_file IGVFFI1924OSEM PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1924OSEM.tsv
 D3   guide_quantifications tabular_file IGVFFI7077LERW PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI7077LERW.tsv
 D0 per_cell_quality_report tabular_file IGVFFI3905TBLN PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI3905TBLN.tsv
 D1 per_cell_quality_report tabular_file IGVFFI0627NEZT PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI0627NEZT.tsv
 D2 per_cell_quality_report tabular_file IGVFFI1451CYGW PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1451CYGW.tsv
 D3 per_cell_quality_report tabular_file IGVFFI8282ATZQ PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI8282ATZQ.tsv

```

## Markdown cell 25

### 9A. Validate the file PATCH operation paths

Document operations use `input_file` because they are JSON. File
PATCH operations use `input_tsv`. This cell verifies that every
`PATCH_FILE` operation has an existing TSV and adds an in-memory
`input_file` compatibility alias for recovery from older sessions.

## Code cell 26

```python
if (
    "input_tsv"
    not in patch_operations.columns
):
    raise RuntimeError(
        "patch_operations is missing "
        "the required input_tsv column"
    )

# Compatibility alias only; input_tsv remains canonical.
patch_operations["input_file"] = (
    patch_operations["input_tsv"]
)

patch_path_rows = []

for operation in (
    patch_operations.to_dict(
        "records"
    )
):
    input_tsv_text = str(
        operation.get(
            "input_tsv",
            "",
        )
    ).strip()

    requires_file = (
        operation["action"]
        == "PATCH_FILE"
    )
    input_exists = (
        Path(input_tsv_text).is_file()
        if input_tsv_text
        else False
    )

    patch_path_rows.append(
        {
            "day": operation["day"],
            "format_role": (
                operation[
                    "format_role"
                ]
            ),
            "accession": (
                operation["accession"]
            ),
            "action": operation["action"],
            "input_tsv": input_tsv_text,
            "requires_file": requires_file,
            "input_exists": input_exists,
            "status": (
                "PASS"
                if (
                    not requires_file
                    or input_exists
                )
                else "FAIL"
            ),
        }
    )

patch_path_preflight = (
    pd.DataFrame(
        patch_path_rows
    )
)
PATCH_PATH_PREFLIGHT_PATH = (
    WORK_DIR
    / "file_patch_path_preflight.tsv"
)
patch_path_preflight.to_csv(
    PATCH_PATH_PREFLIGHT_PATH,
    sep="\t",
    index=False,
    lineterminator="\n",
)

print(
    patch_path_preflight.to_string(
        index=False
    )
)
print(
    "\nSaved:",
    PATCH_PATH_PREFLIGHT_PATH,
)

if not patch_path_preflight[
    "status"
].eq("PASS").all():
    raise RuntimeError(
        "At least one PATCH operation "
        "has a missing TSV path"
    )
```

### Text output (stream)

```text
day             format_role      accession     action                                                                                                                                                                 input_tsv  requires_file  input_exists status
 D0          annotated_h5ad IGVFFI7756RKHT PATCH_FILE  /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI7756RKHT.tsv           True          True   PASS
 D1          annotated_h5ad IGVFFI4783YGLZ PATCH_FILE  /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI4783YGLZ.tsv           True          True   PASS
 D2          annotated_h5ad IGVFFI5165BSGC PATCH_FILE  /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI5165BSGC.tsv           True          True   PASS
 D3          annotated_h5ad IGVFFI9086SWZC PATCH_FILE  /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI9086SWZC.tsv           True          True   PASS
 D0   guide_quantifications IGVFFI1608GIRE PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1608GIRE.tsv           True          True   PASS
 D1   guide_quantifications IGVFFI6431OXJN PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI6431OXJN.tsv           True          True   PASS
 D2   guide_quantifications IGVFFI1924OSEM PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1924OSEM.tsv           True          True   PASS
 D3   guide_quantifications IGVFFI7077LERW PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI7077LERW.tsv           True          True   PASS
 D0 per_cell_quality_report IGVFFI3905TBLN PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI3905TBLN.tsv           True          True   PASS
 D1 per_cell_quality_report IGVFFI0627NEZT PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI0627NEZT.tsv           True          True   PASS
 D2 per_cell_quality_report IGVFFI1451CYGW PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1451CYGW.tsv           True          True   PASS
 D3 per_cell_quality_report IGVFFI8282ATZQ PATCH_FILE /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI8282ATZQ.tsv           True          True   PASS

Saved: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/file_patch_path_preflight.tsv

```

## Code cell 27

```python
RUN_PATCH_DRY_RUNS = True
PATCH_DRY_RUN_MARKER = WORK_DIR / "PATCH_DRY_RUN_PASSED.json"

if not RUN_PATCH_DRY_RUNS:
    print("Set RUN_PATCH_DRY_RUNS=True and rerun this cell.")
else:
    results = []
    for operation in patch_operations.to_dict("records"):
        if operation["action"] == "SKIP_ALREADY_LINKED":
            results.append({**operation, "success": True, "log_path": ""})
            continue
        input_tsv_text = str(
            operation.get("input_tsv", "")
            or operation.get("input_file", "")
        ).strip()
        if not input_tsv_text:
            raise RuntimeError(
                "PATCH operation is missing "
                f"its input TSV path: {operation}"
            )
        input_tsv_path = Path(
            input_tsv_text
        )
        if not input_tsv_path.is_file():
            raise FileNotFoundError(
                input_tsv_path
            )

        result = run_iu_register(
            profile=operation["profile"],
            input_tsv=input_tsv_path,
            dry_run=True,
            patch=True,
            label=f"patch_{operation['format_role']}_{operation['accession']}",
        )
        results.append({**operation, **result})

    result_frame = pd.DataFrame(results)
    PATCH_DRY_RUN_RESULTS_PATH = WORK_DIR / "patch_dry_run_results.tsv"
    result_frame.to_csv(
        PATCH_DRY_RUN_RESULTS_PATH, sep="\t", index=False, lineterminator="\n"
    )
    if not result_frame["success"].all():
        raise RuntimeError("At least one file-patch dry run failed")
    marker = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "patch_plan_md5": md5sum(PATCH_PLAN_PATH),
        "operation_md5s": {
            row["accession"]: md5sum(Path(row["input_tsv"])) if row["input_tsv"] else ""
            for row in patch_operation_rows
        },
    }
    PATCH_DRY_RUN_MARKER.write_text(json.dumps(marker, indent=2) + "\n")
    print("All file-patch dry runs passed:", PATCH_DRY_RUN_MARKER)
```

### Text output (stream)

```text

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p matrix_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI7756RKHT.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_annotated_h5ad_IGVFFI7756RKHT_dryrun_20260930_092654.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:26:55,380:iu_debug:	DRY RUN is enabled.
2026-09-30 09:26:55,380:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:26:55,380:iu_debug:	submission=True: In submission mode.
2026-09-30 09:26:56,668:iu_debug:	
IN patch()
2026-09-30 09:26:56,668:iu_debug:	>>>>>>GET IGVFFI7756RKHT From DACC with URL https://api.data.igvf.org/IGVFFI7756RKHT/?format=json&datastore=database&frame=edit
2026-09-30 09:26:56,838:iu_debug:	<<<<<< PATCHING IGVFFI7756RKHT To IGVF database with URL https://api.data.igvf.org/IGVFFI7756RKHT and this payload:

{
  "file_format_specifications": [
    "e90fd5fc-23b7-4bb7-9444-64db07abc82a"
  ]
}


2026-09-30 09:26:56,838:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p matrix_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI4783YGLZ.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_annotated_h5ad_IGVFFI4783YGLZ_dryrun_20260930_092656.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:26:57,345:iu_debug:	DRY RUN is enabled.
2026-09-30 09:26:57,345:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:26:57,345:iu_debug:	submission=True: In submission mode.
2026-09-30 09:26:58,628:iu_debug:	
IN patch()
2026-09-30 09:26:58,628:iu_debug:	>>>>>>GET IGVFFI4783YGLZ From DACC with URL https://api.data.igvf.org/IGVFFI4783YGLZ/?format=json&datastore=database&frame=edit
2026-09-30 09:26:58,841:iu_debug:	<<<<<< PATCHING IGVFFI4783YGLZ To IGVF database with URL https://api.data.igvf.org/IGVFFI4783YGLZ and this payload:

{
  "file_format_specifications": [
    "e90fd5fc-23b7-4bb7-9444-64db07abc82a"
  ]
}


2026-09-30 09:26:58,841:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p matrix_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI5165BSGC.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_annotated_h5ad_IGVFFI5165BSGC_dryrun_20260930_092658.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:26:59,373:iu_debug:	DRY RUN is enabled.
2026-09-30 09:26:59,373:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:26:59,373:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:00,785:iu_debug:	
IN patch()
2026-09-30 09:27:00,785:iu_debug:	>>>>>>GET IGVFFI5165BSGC From DACC with URL https://api.data.igvf.org/IGVFFI5165BSGC/?format=json&datastore=database&frame=edit
2026-09-30 09:27:00,955:iu_debug:	<<<<<< PATCHING IGVFFI5165BSGC To IGVF database with URL https://api.data.igvf.org/IGVFFI5165BSGC and this payload:

{
  "file_format_specifications": [
    "e90fd5fc-23b7-4bb7-9444-64db07abc82a"
  ]
}


2026-09-30 09:27:00,955:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p matrix_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI9086SWZC.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_annotated_h5ad_IGVFFI9086SWZC_dryrun_20260930_092701.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:01,496:iu_debug:	DRY RUN is enabled.
2026-09-30 09:27:01,496:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:01,496:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:02,203:iu_debug:	
IN patch()
2026-09-30 09:27:02,203:iu_debug:	>>>>>>GET IGVFFI9086SWZC From DACC with URL https://api.data.igvf.org/IGVFFI9086SWZC/?format=json&datastore=database&frame=edit
2026-09-30 09:27:02,396:iu_debug:	<<<<<< PATCHING IGVFFI9086SWZC To IGVF database with URL https://api.data.igvf.org/IGVFFI9086SWZC and this payload:

{
  "file_format_specifications": [
    "e90fd5fc-23b7-4bb7-9444-64db07abc82a"
  ]
}


2026-09-30 09:27:02,396:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1608GIRE.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_guide_quantifications_IGVFFI1608GIRE_dryrun_20260930_092702.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:02,934:iu_debug:	DRY RUN is enabled.
2026-09-30 09:27:02,935:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:02,935:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:04,305:iu_debug:	
IN patch()
2026-09-30 09:27:04,305:iu_debug:	>>>>>>GET IGVFFI1608GIRE From DACC with URL https://api.data.igvf.org/IGVFFI1608GIRE/?format=json&datastore=database&frame=edit
2026-09-30 09:27:04,461:iu_debug:	<<<<<< PATCHING IGVFFI1608GIRE To IGVF database with URL https://api.data.igvf.org/IGVFFI1608GIRE and this payload:

{
  "file_format_specifications": [
    "4ecb21f4-7a3f-49d3-bdbc-3161d349311d"
  ]
}


2026-09-30 09:27:04,461:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI6431OXJN.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_guide_quantifications_IGVFFI6431OXJN_dryrun_20260930_092704.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:04,987:iu_debug:	DRY RUN is enabled.
2026-09-30 09:27:04,987:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:04,987:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:05,657:iu_debug:	
IN patch()
2026-09-30 09:27:05,657:iu_debug:	>>>>>>GET IGVFFI6431OXJN From DACC with URL https://api.data.igvf.org/IGVFFI6431OXJN/?format=json&datastore=database&frame=edit
2026-09-30 09:27:05,802:iu_debug:	<<<<<< PATCHING IGVFFI6431OXJN To IGVF database with URL https://api.data.igvf.org/IGVFFI6431OXJN and this payload:

{
  "file_format_specifications": [
    "4ecb21f4-7a3f-49d3-bdbc-3161d349311d"
  ]
}


2026-09-30 09:27:05,802:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1924OSEM.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_guide_quantifications_IGVFFI1924OSEM_dryrun_20260930_092705.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:06,343:iu_debug:	DRY RUN is enabled.
2026-09-30 09:27:06,343:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:06,343:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:07,858:iu_debug:	
IN patch()
2026-09-30 09:27:07,858:iu_debug:	>>>>>>GET IGVFFI1924OSEM From DACC with URL https://api.data.igvf.org/IGVFFI1924OSEM/?format=json&datastore=database&frame=edit
2026-09-30 09:27:08,019:iu_debug:	<<<<<< PATCHING IGVFFI1924OSEM To IGVF database with URL https://api.data.igvf.org/IGVFFI1924OSEM and this payload:

{
  "file_format_specifications": [
    "4ecb21f4-7a3f-49d3-bdbc-3161d349311d"
  ]
}


2026-09-30 09:27:08,019:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI7077LERW.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_guide_quantifications_IGVFFI7077LERW_dryrun_20260930_092708.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:08,554:iu_debug:	DRY RUN is enabled.
2026-09-30 09:27:08,555:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:08,555:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:09,994:iu_debug:	
IN patch()
2026-09-30 09:27:09,994:iu_debug:	>>>>>>GET IGVFFI7077LERW From DACC with URL https://api.data.igvf.org/IGVFFI7077LERW/?format=json&datastore=database&frame=edit
2026-09-30 09:27:10,155:iu_debug:	<<<<<< PATCHING IGVFFI7077LERW To IGVF database with URL https://api.data.igvf.org/IGVFFI7077LERW and this payload:

{
  "file_format_specifications": [
    "4ecb21f4-7a3f-49d3-bdbc-3161d349311d"
  ]
}


2026-09-30 09:27:10,155:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI3905TBLN.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_per_cell_quality_report_IGVFFI3905TBLN_dryrun_20260930_092710.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:10,670:iu_debug:	DRY RUN is enabled.
2026-09-30 09:27:10,670:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:10,670:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:12,177:iu_debug:	
IN patch()
2026-09-30 09:27:12,177:iu_debug:	>>>>>>GET IGVFFI3905TBLN From DACC with URL https://api.data.igvf.org/IGVFFI3905TBLN/?format=json&datastore=database&frame=edit
2026-09-30 09:27:12,335:iu_debug:	<<<<<< PATCHING IGVFFI3905TBLN To IGVF database with URL https://api.data.igvf.org/IGVFFI3905TBLN and this payload:

{
  "file_format_specifications": [
    "d860bb2b-8907-4bf2-88dd-5686e56b102f"
  ]
}


2026-09-30 09:27:12,335:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI0627NEZT.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_per_cell_quality_report_IGVFFI0627NEZT_dryrun_20260930_092712.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:12,850:iu_debug:	DRY RUN is enabled.
2026-09-30 09:27:12,850:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:12,850:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:14,222:iu_debug:	
IN patch()
2026-09-30 09:27:14,222:iu_debug:	>>>>>>GET IGVFFI0627NEZT From DACC with URL https://api.data.igvf.org/IGVFFI0627NEZT/?format=json&datastore=database&frame=edit
2026-09-30 09:27:14,380:iu_debug:	<<<<<< PATCHING IGVFFI0627NEZT To IGVF database with URL https://api.data.igvf.org/IGVFFI0627NEZT and this payload:

{
  "file_format_specifications": [
    "d860bb2b-8907-4bf2-88dd-5686e56b102f"
  ]
}


2026-09-30 09:27:14,380:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1451CYGW.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_per_cell_quality_report_IGVFFI1451CYGW_dryrun_20260930_092714.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:14,935:iu_debug:	DRY RUN is enabled.
2026-09-30 09:27:14,936:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:14,936:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:16,319:iu_debug:	
IN patch()
2026-09-30 09:27:16,319:iu_debug:	>>>>>>GET IGVFFI1451CYGW From DACC with URL https://api.data.igvf.org/IGVFFI1451CYGW/?format=json&datastore=database&frame=edit
2026-09-30 09:27:16,473:iu_debug:	<<<<<< PATCHING IGVFFI1451CYGW To IGVF database with URL https://api.data.igvf.org/IGVFFI1451CYGW and this payload:

{
  "file_format_specifications": [
    "d860bb2b-8907-4bf2-88dd-5686e56b102f"
  ]
}


2026-09-30 09:27:16,474:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI8282ATZQ.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_per_cell_quality_report_IGVFFI8282ATZQ_dryrun_20260930_092716.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:17,018:iu_debug:	DRY RUN is enabled.
2026-09-30 09:27:17,019:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:17,019:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:18,219:iu_debug:	
IN patch()
2026-09-30 09:27:18,219:iu_debug:	>>>>>>GET IGVFFI8282ATZQ From DACC with URL https://api.data.igvf.org/IGVFFI8282ATZQ/?format=json&datastore=database&frame=edit
2026-09-30 09:27:18,375:iu_debug:	<<<<<< PATCHING IGVFFI8282ATZQ To IGVF database with URL https://api.data.igvf.org/IGVFFI8282ATZQ and this payload:

{
  "file_format_specifications": [
    "d860bb2b-8907-4bf2-88dd-5686e56b102f"
  ]
}


2026-09-30 09:27:18,375:iu_debug:	DRY RUN is enabled.

All file-patch dry runs passed: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/PATCH_DRY_RUN_PASSED.json

```

## Code cell 28

```python
REQUIRED_PATCH_CONFIRMATION = "ATTACH_3_FORMAT_DOCUMENTS_TO_12_300GENE_FILES"

PATCH_PRODUCTION_CONFIRMATION = (
    "ATTACH_3_FORMAT_DOCUMENTS_TO_12_300GENE_FILES"
)

if PATCH_PRODUCTION_CONFIRMATION != REQUIRED_PATCH_CONFIRMATION:
    print("Production PATCH is blocked.")
    print("Set PATCH_PRODUCTION_CONFIRMATION to:")
    print(REQUIRED_PATCH_CONFIRMATION)
else:
    if not PATCH_DRY_RUN_MARKER.is_file():
        raise RuntimeError("Patch dry-run marker is missing")
    marker = json.loads(PATCH_DRY_RUN_MARKER.read_text())
    if marker["patch_plan_md5"] != md5sum(PATCH_PLAN_PATH):
        raise RuntimeError("Patch plan changed after dry run")

    production_results = []
    for operation in patch_operations.to_dict("records"):
        if operation["action"] == "SKIP_ALREADY_LINKED":
            production_results.append({**operation, "success": True, "log_path": ""})
            continue
        input_tsv_text = str(
            operation.get("input_tsv", "")
            or operation.get("input_file", "")
        ).strip()
        if not input_tsv_text:
            raise RuntimeError(
                "PATCH operation is missing "
                f"its input TSV path: {operation}"
            )
        input_tsv = Path(input_tsv_text)
        if not input_tsv.is_file():
            raise FileNotFoundError(
                input_tsv
            )

        if marker["operation_md5s"][operation["accession"]] != md5sum(input_tsv):
            raise RuntimeError(f"Patch TSV changed after dry run: {input_tsv}")
        result = run_iu_register(
            profile=operation["profile"],
            input_tsv=input_tsv,
            dry_run=False,
            patch=True,
            label=f"patch_{operation['format_role']}_{operation['accession']}",
        )
        production_results.append({**operation, **result})

    pd.DataFrame(production_results).to_csv(
        WORK_DIR / "patch_production_results.tsv",
        sep="\t", index=False, lineterminator="\n",
    )
    print("Production PATCH operations completed. Run the reconciliation cell next.")
```

### Text output (stream)

```text

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p matrix_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI7756RKHT.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_annotated_h5ad_IGVFFI7756RKHT_submit_20260930_092718.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:18,949:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:18,950:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:20,243:iu_debug:	
IN patch()
2026-09-30 09:27:20,244:iu_debug:	>>>>>>GET IGVFFI7756RKHT From DACC with URL https://api.data.igvf.org/IGVFFI7756RKHT/?format=json&datastore=database&frame=edit
2026-09-30 09:27:20,438:iu_debug:	<<<<<< PATCHING IGVFFI7756RKHT To IGVF database with URL https://api.data.igvf.org/IGVFFI7756RKHT and this payload:

{
  "file_format_specifications": [
    "e90fd5fc-23b7-4bb7-9444-64db07abc82a"
  ]
}


2026-09-30 09:27:20,989:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p matrix_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI4783YGLZ.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_annotated_h5ad_IGVFFI4783YGLZ_submit_20260930_092721.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:21,512:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:21,512:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:23,032:iu_debug:	
IN patch()
2026-09-30 09:27:23,033:iu_debug:	>>>>>>GET IGVFFI4783YGLZ From DACC with URL https://api.data.igvf.org/IGVFFI4783YGLZ/?format=json&datastore=database&frame=edit
2026-09-30 09:27:23,208:iu_debug:	<<<<<< PATCHING IGVFFI4783YGLZ To IGVF database with URL https://api.data.igvf.org/IGVFFI4783YGLZ and this payload:

{
  "file_format_specifications": [
    "e90fd5fc-23b7-4bb7-9444-64db07abc82a"
  ]
}


2026-09-30 09:27:23,604:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p matrix_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI5165BSGC.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_annotated_h5ad_IGVFFI5165BSGC_submit_20260930_092723.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:24,174:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:24,174:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:25,680:iu_debug:	
IN patch()
2026-09-30 09:27:25,680:iu_debug:	>>>>>>GET IGVFFI5165BSGC From DACC with URL https://api.data.igvf.org/IGVFFI5165BSGC/?format=json&datastore=database&frame=edit
2026-09-30 09:27:25,851:iu_debug:	<<<<<< PATCHING IGVFFI5165BSGC To IGVF database with URL https://api.data.igvf.org/IGVFFI5165BSGC and this payload:

{
  "file_format_specifications": [
    "e90fd5fc-23b7-4bb7-9444-64db07abc82a"
  ]
}


2026-09-30 09:27:26,258:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p matrix_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_matrix_file_IGVFFI9086SWZC.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_annotated_h5ad_IGVFFI9086SWZC_submit_20260930_092726.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:26,831:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:26,831:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:28,317:iu_debug:	
IN patch()
2026-09-30 09:27:28,318:iu_debug:	>>>>>>GET IGVFFI9086SWZC From DACC with URL https://api.data.igvf.org/IGVFFI9086SWZC/?format=json&datastore=database&frame=edit
2026-09-30 09:27:28,519:iu_debug:	<<<<<< PATCHING IGVFFI9086SWZC To IGVF database with URL https://api.data.igvf.org/IGVFFI9086SWZC and this payload:

{
  "file_format_specifications": [
    "e90fd5fc-23b7-4bb7-9444-64db07abc82a"
  ]
}


2026-09-30 09:27:29,066:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1608GIRE.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_guide_quantifications_IGVFFI1608GIRE_submit_20260930_092729.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:29,625:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:29,625:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:31,083:iu_debug:	
IN patch()
2026-09-30 09:27:32,462:iu_debug:	>>>>>>GET IGVFFI1608GIRE From DACC with URL https://api.data.igvf.org/IGVFFI1608GIRE/?format=json&datastore=database&frame=edit
2026-09-30 09:27:32,622:iu_debug:	<<<<<< PATCHING IGVFFI1608GIRE To IGVF database with URL https://api.data.igvf.org/IGVFFI1608GIRE and this payload:

{
  "file_format_specifications": [
    "4ecb21f4-7a3f-49d3-bdbc-3161d349311d"
  ]
}


2026-09-30 09:27:32,970:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI6431OXJN.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_guide_quantifications_IGVFFI6431OXJN_submit_20260930_092733.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:33,507:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:33,507:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:34,196:iu_debug:	
IN patch()
2026-09-30 09:27:34,196:iu_debug:	>>>>>>GET IGVFFI6431OXJN From DACC with URL https://api.data.igvf.org/IGVFFI6431OXJN/?format=json&datastore=database&frame=edit
2026-09-30 09:27:34,357:iu_debug:	<<<<<< PATCHING IGVFFI6431OXJN To IGVF database with URL https://api.data.igvf.org/IGVFFI6431OXJN and this payload:

{
  "file_format_specifications": [
    "4ecb21f4-7a3f-49d3-bdbc-3161d349311d"
  ]
}


2026-09-30 09:27:34,836:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1924OSEM.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_guide_quantifications_IGVFFI1924OSEM_submit_20260930_092734.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:35,398:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:35,398:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:36,821:iu_debug:	
IN patch()
2026-09-30 09:27:37,745:iu_debug:	>>>>>>GET IGVFFI1924OSEM From DACC with URL https://api.data.igvf.org/IGVFFI1924OSEM/?format=json&datastore=database&frame=edit
2026-09-30 09:27:37,903:iu_debug:	<<<<<< PATCHING IGVFFI1924OSEM To IGVF database with URL https://api.data.igvf.org/IGVFFI1924OSEM and this payload:

{
  "file_format_specifications": [
    "4ecb21f4-7a3f-49d3-bdbc-3161d349311d"
  ]
}


2026-09-30 09:27:38,239:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI7077LERW.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_guide_quantifications_IGVFFI7077LERW_submit_20260930_092738.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:38,757:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:38,757:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:40,056:iu_debug:	
IN patch()
2026-09-30 09:27:40,057:iu_debug:	>>>>>>GET IGVFFI7077LERW From DACC with URL https://api.data.igvf.org/IGVFFI7077LERW/?format=json&datastore=database&frame=edit
2026-09-30 09:27:40,216:iu_debug:	<<<<<< PATCHING IGVFFI7077LERW To IGVF database with URL https://api.data.igvf.org/IGVFFI7077LERW and this payload:

{
  "file_format_specifications": [
    "4ecb21f4-7a3f-49d3-bdbc-3161d349311d"
  ]
}


2026-09-30 09:27:40,577:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI3905TBLN.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_per_cell_quality_report_IGVFFI3905TBLN_submit_20260930_092740.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:41,084:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:41,084:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:42,538:iu_debug:	
IN patch()
2026-09-30 09:27:42,539:iu_debug:	>>>>>>GET IGVFFI3905TBLN From DACC with URL https://api.data.igvf.org/IGVFFI3905TBLN/?format=json&datastore=database&frame=edit
2026-09-30 09:27:42,692:iu_debug:	<<<<<< PATCHING IGVFFI3905TBLN To IGVF database with URL https://api.data.igvf.org/IGVFFI3905TBLN and this payload:

{
  "file_format_specifications": [
    "d860bb2b-8907-4bf2-88dd-5686e56b102f"
  ]
}


2026-09-30 09:27:43,140:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI0627NEZT.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_per_cell_quality_report_IGVFFI0627NEZT_submit_20260930_092743.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:43,673:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:43,673:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:44,909:iu_debug:	
IN patch()
2026-09-30 09:27:44,910:iu_debug:	>>>>>>GET IGVFFI0627NEZT From DACC with URL https://api.data.igvf.org/IGVFFI0627NEZT/?format=json&datastore=database&frame=edit
2026-09-30 09:27:45,114:iu_debug:	<<<<<< PATCHING IGVFFI0627NEZT To IGVF database with URL https://api.data.igvf.org/IGVFFI0627NEZT and this payload:

{
  "file_format_specifications": [
    "d860bb2b-8907-4bf2-88dd-5686e56b102f"
  ]
}


2026-09-30 09:27:45,387:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI1451CYGW.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_per_cell_quality_report_IGVFFI1451CYGW_submit_20260930_092745.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:45,895:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:45,896:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:46,581:iu_debug:	
IN patch()
2026-09-30 09:27:46,581:iu_debug:	>>>>>>GET IGVFFI1451CYGW From DACC with URL https://api.data.igvf.org/IGVFFI1451CYGW/?format=json&datastore=database&frame=edit
2026-09-30 09:27:46,737:iu_debug:	<<<<<< PATCHING IGVFFI1451CYGW To IGVF database with URL https://api.data.igvf.org/IGVFFI1451CYGW and this payload:

{
  "file_format_specifications": [
    "d860bb2b-8907-4bf2-88dd-5686e56b102f"
  ]
}


2026-09-30 09:27:47,194:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/one_row_tsv/PATCH_tabular_file_IGVFFI8282ATZQ.tsv
Log: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs/patch_per_cell_quality_report_IGVFFI8282ATZQ_submit_20260930_092747.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 09:27:47,711:iu_debug:	submission=False: In non-submission mode.
2026-09-30 09:27:47,711:iu_debug:	submission=True: In submission mode.
2026-09-30 09:27:48,389:iu_debug:	
IN patch()
2026-09-30 09:27:48,389:iu_debug:	>>>>>>GET IGVFFI8282ATZQ From DACC with URL https://api.data.igvf.org/IGVFFI8282ATZQ/?format=json&datastore=database&frame=edit
2026-09-30 09:27:48,543:iu_debug:	<<<<<< PATCHING IGVFFI8282ATZQ To IGVF database with URL https://api.data.igvf.org/IGVFFI8282ATZQ and this payload:

{
  "file_format_specifications": [
    "d860bb2b-8907-4bf2-88dd-5686e56b102f"
  ]
}


2026-09-30 09:27:49,024:iu_debug:	Success.

Production PATCH operations completed. Run the reconciliation cell next.

```

## Markdown cell 29

## 10. Final read-only reconciliation

## Code cell 30

```python
reconciliation_rows = []
for record in target_files.to_dict("records"):
    role_key = record["format_role"]
    file_record = get_with_retry(record["accession"], attempts=5)
    document_record = get_with_retry(DOCUMENT_SPECS[role_key]["document_alias"], attempts=5)
    expected_doc = portal_record_id(document_record)
    existing_docs = (file_record or {}).get("file_format_specifications", [])
    linked = bool(file_record and document_record and link_set_contains(existing_docs, expected_doc))
    reconciliation_rows.append(
        {
            "day": record["day"],
            "format_role": role_key,
            "profile": record["profile"],
            "accession": record["accession"],
            "file_portal_url": portal_record_url(file_record),
            "document_alias": DOCUMENT_SPECS[role_key]["document_alias"],
            "document_id": expected_doc,
            "document_portal_url": portal_record_url(document_record),
            "portal_file_format_specifications": ",".join(
                sorted(direct_link_tokens(existing_docs))
            ),
            "verification": "PASS" if linked else "FAIL",
        }
    )

reconciliation = pd.DataFrame(reconciliation_rows)
FINAL_RECONCILIATION_PATH = OUTPUT_ROOT / "file_format_specification_reconciliation_12_files.csv"
reconciliation.to_csv(FINAL_RECONCILIATION_PATH, index=False, lineterminator="\n")
print(reconciliation.to_string(index=False))
print("\nSaved:", FINAL_RECONCILIATION_PATH)

passed = int(reconciliation["verification"].eq("PASS").sum())
print(f"\nLinked files: {passed}/12")
if passed == 12:
    marker = OUTPUT_ROOT / "FILE_FORMAT_SPECIFICATIONS_LINKED_12_OF_12.txt"
    marker.write_text(
        "All 12 requested analysis files link to their matching file-format specification Documents.\n"
        f"Results: {FINAL_RECONCILIATION_PATH}\n",
        encoding="utf-8",
    )
    print("Completion marker:", marker)
```

### Text output (stream)

```text
day             format_role      profile      accession                                     file_portal_url                                                                         document_alias                                      document_id                                                   document_portal_url    portal_file_format_specifications verification
 D0          annotated_h5ad  matrix_file IGVFFI7756RKHT  https://data.igvf.org/matrix-files/IGVFFI7756RKHT/ jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1 /documents/e90fd5fc-23b7-4bb7-9444-64db07abc82a/ https://data.igvf.org/documents/e90fd5fc-23b7-4bb7-9444-64db07abc82a/ e90fd5fc-23b7-4bb7-9444-64db07abc82a         PASS
 D1          annotated_h5ad  matrix_file IGVFFI4783YGLZ  https://data.igvf.org/matrix-files/IGVFFI4783YGLZ/ jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1 /documents/e90fd5fc-23b7-4bb7-9444-64db07abc82a/ https://data.igvf.org/documents/e90fd5fc-23b7-4bb7-9444-64db07abc82a/ e90fd5fc-23b7-4bb7-9444-64db07abc82a         PASS
 D2          annotated_h5ad  matrix_file IGVFFI5165BSGC  https://data.igvf.org/matrix-files/IGVFFI5165BSGC/ jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1 /documents/e90fd5fc-23b7-4bb7-9444-64db07abc82a/ https://data.igvf.org/documents/e90fd5fc-23b7-4bb7-9444-64db07abc82a/ e90fd5fc-23b7-4bb7-9444-64db07abc82a         PASS
 D3          annotated_h5ad  matrix_file IGVFFI9086SWZC  https://data.igvf.org/matrix-files/IGVFFI9086SWZC/ jesse-engreitz:document_WTC11-EC_300gene_annotated_cell_gene_guide_h5ad_file_format_v1 /documents/e90fd5fc-23b7-4bb7-9444-64db07abc82a/ https://data.igvf.org/documents/e90fd5fc-23b7-4bb7-9444-64db07abc82a/ e90fd5fc-23b7-4bb7-9444-64db07abc82a         PASS
 D0   guide_quantifications tabular_file IGVFFI1608GIRE https://data.igvf.org/tabular-files/IGVFFI1608GIRE/          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1 /documents/4ecb21f4-7a3f-49d3-bdbc-3161d349311d/ https://data.igvf.org/documents/4ecb21f4-7a3f-49d3-bdbc-3161d349311d/ 4ecb21f4-7a3f-49d3-bdbc-3161d349311d         PASS
 D1   guide_quantifications tabular_file IGVFFI6431OXJN https://data.igvf.org/tabular-files/IGVFFI6431OXJN/          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1 /documents/4ecb21f4-7a3f-49d3-bdbc-3161d349311d/ https://data.igvf.org/documents/4ecb21f4-7a3f-49d3-bdbc-3161d349311d/ 4ecb21f4-7a3f-49d3-bdbc-3161d349311d         PASS
 D2   guide_quantifications tabular_file IGVFFI1924OSEM https://data.igvf.org/tabular-files/IGVFFI1924OSEM/          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1 /documents/4ecb21f4-7a3f-49d3-bdbc-3161d349311d/ https://data.igvf.org/documents/4ecb21f4-7a3f-49d3-bdbc-3161d349311d/ 4ecb21f4-7a3f-49d3-bdbc-3161d349311d         PASS
 D3   guide_quantifications tabular_file IGVFFI7077LERW https://data.igvf.org/tabular-files/IGVFFI7077LERW/          jesse-engreitz:document_WTC11-EC_300gene_guide_quantifications_file_format_v1 /documents/4ecb21f4-7a3f-49d3-bdbc-3161d349311d/ https://data.igvf.org/documents/4ecb21f4-7a3f-49d3-bdbc-3161d349311d/ 4ecb21f4-7a3f-49d3-bdbc-3161d349311d         PASS
 D0 per_cell_quality_report tabular_file IGVFFI3905TBLN https://data.igvf.org/tabular-files/IGVFFI3905TBLN/        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1 /documents/d860bb2b-8907-4bf2-88dd-5686e56b102f/ https://data.igvf.org/documents/d860bb2b-8907-4bf2-88dd-5686e56b102f/ d860bb2b-8907-4bf2-88dd-5686e56b102f         PASS
 D1 per_cell_quality_report tabular_file IGVFFI0627NEZT https://data.igvf.org/tabular-files/IGVFFI0627NEZT/        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1 /documents/d860bb2b-8907-4bf2-88dd-5686e56b102f/ https://data.igvf.org/documents/d860bb2b-8907-4bf2-88dd-5686e56b102f/ d860bb2b-8907-4bf2-88dd-5686e56b102f         PASS
 D2 per_cell_quality_report tabular_file IGVFFI1451CYGW https://data.igvf.org/tabular-files/IGVFFI1451CYGW/        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1 /documents/d860bb2b-8907-4bf2-88dd-5686e56b102f/ https://data.igvf.org/documents/d860bb2b-8907-4bf2-88dd-5686e56b102f/ d860bb2b-8907-4bf2-88dd-5686e56b102f         PASS
 D3 per_cell_quality_report tabular_file IGVFFI8282ATZQ https://data.igvf.org/tabular-files/IGVFFI8282ATZQ/        jesse-engreitz:document_WTC11-EC_300gene_per_cell_quality_report_file_format_v1 /documents/d860bb2b-8907-4bf2-88dd-5686e56b102f/ https://data.igvf.org/documents/d860bb2b-8907-4bf2-88dd-5686e56b102f/ d860bb2b-8907-4bf2-88dd-5686e56b102f         PASS

Saved: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/file_format_specification_reconciliation_12_files.csv

Linked files: 12/12
Completion marker: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/FILE_FORMAT_SPECIFICATIONS_LINKED_12_OF_12.txt

```

## Markdown cell 31

## 11. Output summary

## Code cell 32

```python
summary_paths = [
    ("target inventory", TARGET_FILE_INVENTORY_PATH),
    ("local file audit", LOCAL_FILE_AUDIT_PATH),
    ("tabular schema audit", TABULAR_SCHEMA_AUDIT_PATH),
    ("h5ad file audit", H5AD_FILE_AUDIT_PATH),
    ("h5ad component audit", H5AD_COMPONENT_AUDIT_PATH),
    ("PDF inventory", PDF_INVENTORY_PATH),
    ("Document metadata JSON", DOCUMENT_METADATA_PATH),
    ("Document review-only TSV", DOCUMENT_REVIEW_TSV_PATH),
    ("Document JSON preflight", WORK_DIR / "document_json_preflight.tsv"),
    ("Document plan", DOCUMENT_PLAN_PATH),
    ("file patch plan", PATCH_PLAN_PATH),
    ("final reconciliation", OUTPUT_ROOT / "file_format_specification_reconciliation_12_files.csv"),
]
summary = pd.DataFrame([
    {
        "output": label,
        "path": str(path),
        "exists": path.is_file(),
        "size_bytes": path.stat().st_size if path.is_file() else 0,
    }
    for label, path in summary_paths
])
print(summary.to_string(index=False))
print("\nPDF documents:", PDF_DIR)
print("Logs:", LOG_DIR)
```

### Text output (stream)

```text
                  output                                                                                                                                                                        path  exists  size_bytes
        target inventory                           /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/target_file_inventory.tsv    True        3790
        local file audit                                /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/local_file_audit.tsv    True        4451
    tabular schema audit                    /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/schema_tables/tabular_schema_audit.tsv    True        2791
         h5ad file audit                         /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/schema_tables/h5ad_file_audit.tsv    True        5423
    h5ad component audit                    /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/schema_tables/h5ad_component_audit.tsv    True       10360
           PDF inventory                                   /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/pdf_inventory.tsv    True        1064
  Document metadata JSON            /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/document_file_format_specifications.json    True        1921
Document review-only TSV /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/document_file_format_specifications_review_only.tsv    True        1501
 Document JSON preflight                         /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/document_json_preflight.tsv    True         749
           Document plan                        /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/document_submission_plan.tsv    True         410
         file patch plan                                 /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/_work/file_patch_plan.tsv    True        6221
    final reconciliation     /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/file_format_specification_reconciliation_12_files.csv    True        4288

PDF documents: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/documents
Logs: /oak/stanford/groups/engreitz/Projects/IGVF/CC-Perturb-seq/300gene/analysis_submission/file_format_specifications/logs

```

