# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/prepare_and_submit_both_crispri_screens_igvf_analysis_files.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/prepare_and_submit_both_crispri_screens_igvf_analysis_files.ipynb`
- SHA-256: `c89392481b41b2c5b2d180ad6a69e4e2ec02307195eb793f91839fad7563f697`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Prepare and submit IGVF analysis files for both WTC11-EC CRISPRi screens

This single notebook prepares and submits the four primary analysis files required for the two screens:

| Screen | IGVF content type |
|---|---|
| Focused validation screen | `guide quantifications` |
| Focused validation screen | `differential element quantifications` |
| Genome-wide screen | `guide quantifications` |
| Genome-wide screen | `differential element quantifications` |

Ian's existing AnalysisStepVersion records are used for both screens:

- guide quantifications: `/analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/`
- differential element quantifications: `/analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/`

The notebook first creates and validates all four final `.tsv.gz` files. It then builds one combined four-file submission inventory and handles each Portal record independently as:

- `POST_NEW`
- `PATCH_ASV`
- `SKIP_ALREADY_CORRECT`
- a blocking conflict requiring manual review

Production operations are sequential and explicitly gated. A failure on one record stops the notebook before the next record.

## Markdown cell 1

## Version 2 source-selection fix

The validation-screen source selector now uses the standard MAGeCK RRA file:

```text
/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/
MageckAnalysis/D3pos_v_D3neg/D3pos_D3neg_gene_summary.txt
```

The separate `MLE_D3pos_v_D3neg/D3pos_D3neg.gene_summary.txt` result is intentionally excluded because it is a different MAGeCK analysis and does not match the `neg|...` / `pos|...` schema expected by this submission.

All explicit source overrides now take full precedence over automatic discovery.


## Markdown cell 2

## Version 3 validation guide-count parser fix

`guideCounts.csv` contains at least one unlabeled text column (`Unnamed: 0`).
The earlier parser incorrectly assumed that every column other than the chosen
guide-ID and recognized target column was a numeric sample-count column.

This version instead:

1. identifies the guide-ID column by exact agreement with Jesse's corrected
   guide-metadata IDs;
2. classifies sample-count columns by their actual values;
3. treats other nonnumeric columns as source annotations rather than counts;
4. infers a target-label annotation only when it agrees with the gene label
   parsed from the guide ID;
5. excludes unrecognized source annotations from the upload table but preserves
   them in a QC audit;
6. requires exactly 41 validated numeric count columns.


## Markdown cell 3

## Version 4 validation guide-count reconciliation fix

The attached column profile shows that `guideCounts.csv` has:

- `Guide`: the 1,926-row unique guide identifier;
- `Unnamed: 0`: the 299-level target/gene annotation;
- six aggregate count columns: `D01`, `D02`, `D3neg1`,
  `D3neg2`, `D3pos1`, and `D3pos2`.

This version therefore:

1. matches guide IDs after conservative case/Unicode-dash
   normalization and then writes the exact canonical IDs from
   Jesse's corrected metadata;
2. treats `Unnamed: 0` as the gene annotation rather than an index;
3. expects six aggregate count columns rather than 41 raw PCR
   libraries;
4. cross-checks those six columns against
   `mageck_sgRNA_readcount.txt`;
5. repairs any missing count only when the corresponding MAGeCK
   column agrees exactly at every non-missing row.

The fixed mapping is:

```text
D01     <- D0R1
D02     <- D0R2
D3neg1  <- D3NR1
D3neg2  <- D3NR2
D3pos1  <- D3PR1
D3pos2  <- D3PR2
```

## Markdown cell 4

## Version 5 — explicit NEDD8/NEDD8-MDP1 reconciliation

The attached mismatch audit contains exactly one one-to-one naming
difference:

```text
guideCounts.csv:       NEDD8_+_24701319.23-P1P2
corrected metadata:    NEDD8-MDP1_+_24701319.23-P1P2
```

This version handles that difference with a single explicit curated
alias rather than a broad fuzzy-matching rule. The output guide-count
table uses Jesse's corrected guide ID, while preserving `NEDD8` as the
native MAGeCK analysis group. The guide is counted with the NEDD8
promoter element in the differential-element audit because the native
NEDD8 MAGeCK row contains seven guides.

Every alias application and analysis-group merge is written to a QC
audit. Any additional unmatched guide still stops the notebook.

## Markdown cell 5

## Version 6 — fixed AnalysisSet accessions and canonical genome-wide source

This version hard-codes the production AnalysisSet accessions:

```text
validation screen:  IGVFDS4708PCVX
genome-wide screen: IGVFDS9629VXQK
```

The genome-wide differential-element source is fixed to the documented
canonical results copy:

```text
/oak/stanford/groups/engreitz/Users/dulguun/
230916_Genome-wide_Screen_Analysis/
full_result_files_from_Chad/Pos_vs_Neg/
Pos_vs_Neg.gene_summary.txt
```

The three discovered copies have the same MD5:

```text
257255394c229c685b178e7d91b0c9b2
```

The source-selection helper also now accepts multiple byte-identical
copies and selects the first candidate in the documented priority order.
It still stops when multiple copies differ.

## Markdown cell 6

## Version 7 — explicit credential-entry cell

The credential section now includes a separate executable cell that calls `enter_igvf_credentials()` directly. Running all cells will pause at two hidden prompts before the Portal preflight, preventing the previous `IGVF_API_KEY` / `IGVF_SECRET_KEY` missing-variable error.


## Markdown cell 7

## Version 8 — submitter-safe reference metadata

The prior version used `submitted_assembly=GRCh38`. That field is admin-only on the production Portal. This version removes it and uses the established IGVF reference file `/reference-files/IGVFFI7217ZMJZ/` for all four analysis files. The Portal preflight and final reconciliation verify that link.

## Markdown cell 8

## Version 9 — robust Portal verification and safe resume

This version fixes false post-upload verification failures caused by linked Portal fields being returned as embedded objects rather than plain strings. It queries the editable database frame, compares linked fields using accession/UUID/@id tokens, and saves a partial verification CSV before stopping on a real mismatch.

After an interrupted run, restart the kernel and rerun the notebook through the Portal plan and dry-run sections. Existing matching files are classified as `SKIP_ALREADY_CORRECT`; they are not POSTed again.

## Markdown cell 9

## Version 10 — read-only dependency lookup fix

Version 9 queried every dependency with `frame=edit`. Shared or
DACC-created objects can be readable but not editable by the current
submitter, which can make a valid AnalysisSet, AnalysisStepVersion,
or reference file appear missing during preflight.

Version 10 uses the normal database read endpoint for dependency
lookup and file reconciliation. It still blocks when an object is
genuinely absent. The dependency error now names every missing or
inaccessible record and preserves its lookup error.

Production verification also polls until the newly posted record's
linked metadata becomes visible, reducing false failures caused by
immediate post-submission propagation delay.

## Markdown cell 10

## Version 11 — embedded Portal-link comparison fix

The IGVF API may return linked properties such as `file_set` and
`analysis_step_version` as embedded dictionaries rather than simple
accession/path strings. Version 10 converted those dictionaries to
strings before comparing them, which created a false
`BLOCK_EXISTING_RECORD_CONFLICT`.

Version 11 compares the raw linked objects and writes compact
accession/UUID values to the submission-plan report. An existing
record whose MD5, AnalysisSet, content type, file format,
AnalysisStepVersion, and reference file all match is classified as
`SKIP_ALREADY_CORRECT` and is never POSTed again.

## Markdown cell 11

## 0. Master configuration

Fill in the two AnalysisSet accessions. Source-path overrides may remain blank when the documented Sherlock paths are correct.

Existing-file identifiers should remain blank for new files. Fill one only when that exact Portal file already exists and merely needs Ian's AnalysisStepVersion patched.

## Code cell 12

```python
from pathlib import Path

# -------------------------------------------------------------------------
# REQUIRED production AnalysisSet accessions
# -------------------------------------------------------------------------
VALIDATION_ANALYSIS_SET_INPUT = "IGVFDS4708PCVX"
GENOME_WIDE_ANALYSIS_SET_INPUT = "IGVFDS9629VXQK"

# Established GRCh38 / GENCODE v43 reference file used in the 300-gene metadata.
REFERENCE_FILE_INPUT = "/reference-files/IGVFFI7217ZMJZ/"

# -------------------------------------------------------------------------
# Optional validation-screen source overrides
# -------------------------------------------------------------------------
VALIDATION_GUIDE_COUNTS_OVERRIDE_INPUT = ""
VALIDATION_DIFFERENTIAL_SOURCE_OVERRIDE_INPUT = ""
VALIDATION_CORRECTED_GUIDE_METADATA_OVERRIDE_INPUT = ""

# One documented source-analysis guide-ID alias. Do not generalize this
# mapping: it reconciles one exact row whose spacer/counts are the same but
# whose source-analysis name differs from Jesse's corrected construct name.
VALIDATION_GUIDE_ID_ALIAS_OVERRIDES_INPUT = {
    "NEDD8_+_24701319.23-P1P2": (
        "NEDD8-MDP1_+_24701319.23-P1P2"
    ),
}

# The native validation MAGeCK input grouped that exact guide under NEDD8.
# Preserve that analysis group while using the corrected construct guide ID.
VALIDATION_GUIDE_ANALYSIS_GROUP_OVERRIDES_INPUT = {
    "NEDD8-MDP1_+_24701319.23-P1P2": "NEDD8",
}

# -------------------------------------------------------------------------
# Optional genome-wide source overrides
# -------------------------------------------------------------------------
GENOME_WIDE_GUIDE_COUNTS_OVERRIDE_INPUT = ""
GENOME_WIDE_DIFFERENTIAL_SOURCE_OVERRIDE_INPUT = (
    "/oak/stanford/groups/engreitz/Users/dulguun/"
    "230916_Genome-wide_Screen_Analysis/"
    "full_result_files_from_Chad/Pos_vs_Neg/"
    "Pos_vs_Neg.gene_summary.txt"
)
GENOME_WIDE_DIFFERENTIAL_SOURCE_EXPECTED_MD5_INPUT = (
    "257255394c229c685b178e7d91b0c9b2"
)
GENOME_WIDE_LIBRARY_DESIGN_OVERRIDE_INPUT = ""

# -------------------------------------------------------------------------
# Optional existing Portal records that need only ASV PATCH operations
# -------------------------------------------------------------------------
VALIDATION_EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID_INPUT = ""
VALIDATION_EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID_INPUT = ""
GENOME_WIDE_EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID_INPUT = ""
GENOME_WIDE_EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID_INPUT = ""

# -------------------------------------------------------------------------
# Optional provenance and format-document links
# -------------------------------------------------------------------------
VALIDATION_GUIDE_DERIVED_FROM_INPUT: list[str] = []
VALIDATION_DIFFERENTIAL_DERIVED_FROM_INPUT: list[str] = []
GENOME_WIDE_GUIDE_DERIVED_FROM_INPUT: list[str] = []
GENOME_WIDE_DIFFERENTIAL_DERIVED_FROM_INPUT: list[str] = []

VALIDATION_GUIDE_FORMAT_SPEC_INPUT = ""
VALIDATION_DIFFERENTIAL_FORMAT_SPEC_INPUT = ""
GENOME_WIDE_GUIDE_FORMAT_SPEC_INPUT = ""
GENOME_WIDE_DIFFERENTIAL_FORMAT_SPEC_INPUT = ""

# -------------------------------------------------------------------------
# Output roots
# -------------------------------------------------------------------------
VALIDATION_OUTPUT_ROOT_INPUT = Path(
    "/oak/stanford/groups/engreitz/Users/tri/"
    "IGVF_validation_screen_analysis_submission"
)
GENOME_WIDE_OUTPUT_ROOT_INPUT = Path(
    "/oak/stanford/groups/engreitz/Users/tri/"
    "IGVF_genome_wide_screen_analysis_submission"
)
COMBINED_OUTPUT_ROOT = Path(
    "/oak/stanford/groups/engreitz/Users/tri/"
    "IGVF_both_crispri_screens_analysis_submission"
)

for directory in [
    COMBINED_OUTPUT_ROOT,
    COMBINED_OUTPUT_ROOT / "metadata",
    COMBINED_OUTPUT_ROOT / "logs",
    COMBINED_OUTPUT_ROOT / "work",
    COMBINED_OUTPUT_ROOT / "qc",
]:
    directory.mkdir(parents=True, exist_ok=True)

print("Validation AnalysisSet:", VALIDATION_ANALYSIS_SET_INPUT)
print("Genome-wide AnalysisSet:", GENOME_WIDE_ANALYSIS_SET_INPUT)
print("Reference file:", REFERENCE_FILE_INPUT)
print("Combined submission root:", COMBINED_OUTPUT_ROOT)
```

### Text output (stream)

```text
Validation AnalysisSet: IGVFDS4708PCVX
Genome-wide AnalysisSet: IGVFDS9629VXQK
Reference file: /reference-files/IGVFFI7217ZMJZ/
Combined submission root: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission

```

## Markdown cell 13

# Part A. Focused validation screen

## Markdown cell 14

## 1. Configuration

Before metadata generation, replace `VALIDATION_ANALYSIS_SET` with the production `IGVFDS...` accession for the validation screen.

The notebook deliberately searches only for a file named:

```text
300genes_guide_metadata_v43_corrected.tsv
```

or its `.gz` form. It will not silently use an earlier draft.

## Code cell 15

```python
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from getpass import getpass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import requests

# -------------------------------------------------------------------------
# Documented Sherlock source locations
# -------------------------------------------------------------------------
VALIDATION_ROOT = Path(
    "/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation"
)

GUIDE_COUNTS_SOURCE = VALIDATION_ROOT / "UTestAnalysis/guideCounts.csv"

MAGECK_GUIDE_COUNTS_REFERENCE = (
    VALIDATION_ROOT / "MageckAnalysis/mageck_sgRNA_readcount.txt"
)

# Standard MAGeCK RRA output. Do not substitute the separate MLE result:
# MageckAnalysis/MLE_D3pos_v_D3neg/D3pos_D3neg.gene_summary.txt
DIFFERENTIAL_SOURCE = (
    VALIDATION_ROOT
    / "MageckAnalysis/D3pos_v_D3neg/D3pos_D3neg_gene_summary.txt"
)

VALIDATION_GUIDE_COUNTS_OVERRIDE = VALIDATION_GUIDE_COUNTS_OVERRIDE_INPUT
VALIDATION_DIFFERENTIAL_SOURCE_OVERRIDE = VALIDATION_DIFFERENTIAL_SOURCE_OVERRIDE_INPUT

# Set this when the authoritative corrected guide file is not found by the
# automatic search in Section 3.
CORRECTED_GUIDE_METADATA_OVERRIDE = VALIDATION_CORRECTED_GUIDE_METADATA_OVERRIDE_INPUT

CURATED_GUIDE_ID_ALIASES = dict(
    VALIDATION_GUIDE_ID_ALIAS_OVERRIDES_INPUT
)
CURATED_GUIDE_ANALYSIS_GROUP_OVERRIDES = dict(
    VALIDATION_GUIDE_ANALYSIS_GROUP_OVERRIDES_INPUT
)


# New submission workspace; the source folders are not changed.
OUTPUT_ROOT = Path(VALIDATION_OUTPUT_ROOT_INPUT)

# REQUIRED: replace with the production AnalysisSet accession.
VALIDATION_ANALYSIS_SET = (
    VALIDATION_ANALYSIS_SET_INPUT
    or os.environ.get(
        "IGVF_VALIDATION_ANALYSIS_SET",
        "REPLACE_WITH_VALIDATION_ANALYSIS_SET_IGVFDS",
    )
)

AWARD = "/awards/HG011972/"
LAB = "/labs/jesse-engreitz/"
REFERENCE_FILE = REFERENCE_FILE_INPUT

# Ian-provided AnalysisStepVersion identifiers.
GUIDE_QUANTIFICATION_ASV = (
    "/analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/"
)
DIFFERENTIAL_ELEMENT_ASV = (
    "/analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/"
)

GUIDE_QUANTIFICATION_ALIAS = (
    "jesse-engreitz:tabular_file_WTC11-EC_validation_screen_"
    "guide_quantifications"
)
DIFFERENTIAL_ELEMENT_ALIAS = (
    "jesse-engreitz:tabular_file_WTC11-EC_validation_screen_"
    "D3_CDH5pos_vs_CDH5neg_differential_element_quantifications"
)

# Optional immediate upstream IGVF File identifiers. Leave empty rather than
# guessing. Values are comma-separated in the generated TSV.
GUIDE_QUANTIFICATIONS_DERIVED_FROM = list(VALIDATION_GUIDE_DERIVED_FROM_INPUT)
DIFFERENTIAL_ELEMENTS_DERIVED_FROM = list(VALIDATION_DIFFERENTIAL_DERIVED_FROM_INPUT)

# Optional existing Portal Document identifiers for the table schemas.
GUIDE_QUANTIFICATIONS_FORMAT_SPEC = VALIDATION_GUIDE_FORMAT_SPEC_INPUT
DIFFERENTIAL_ELEMENTS_FORMAT_SPEC = VALIDATION_DIFFERENTIAL_FORMAT_SPEC_INPUT

# If files are already on the Portal and only need their ASV patched, enter
# their accession, UUID, or alias here. Otherwise leave blank.
EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID = VALIDATION_EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID_INPUT
EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID = VALIDATION_EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID_INPUT

# Expected design dimensions.
EXPECTED_GUIDES = 1926
EXPECTED_SAMPLE_COUNT_COLUMNS = 6
EXPECTED_SAMPLE_COUNT_COLUMN_NAMES = [
    "D01",
    "D02",
    "D3neg1",
    "D3neg2",
    "D3pos1",
    "D3pos2",
]
VALIDATION_REFERENCE_COUNT_COLUMN_MAP = {
    "D01": "D0R1",
    "D02": "D0R2",
    "D3neg1": "D3NR1",
    "D3neg2": "D3NR2",
    "D3pos1": "D3PR1",
    "D3pos2": "D3PR2",
}
EXPECTED_MAGECK_ROWS = 299
EXPECTED_BIOLOGICAL_ELEMENT_ROWS = 298
ALLOWED_MULTI_ELEMENT_GENES = {"NOTCH1"}
STRICT_EXPECTED_DIMENSIONS = True

# Existing working igvf_utils environment.
IGVF_PYTHON = Path(
    "/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python"
)
IU_REGISTER = Path(
    "/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/"
    "igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py"
)

# Submission gates.
RUN_IGVF_DRY_RUN = False
REAL_UPLOAD_CONFIRMATION = ""
REQUIRED_REAL_UPLOAD_CONFIRMATION = (
    "UPLOAD_VALIDATION_SCREEN_ANALYSIS_FILES_TO_IGVF_PROD"
)

for directory in [
    OUTPUT_ROOT,
    OUTPUT_ROOT / "files",
    OUTPUT_ROOT / "metadata",
    OUTPUT_ROOT / "documentation",
    OUTPUT_ROOT / "qc",
    OUTPUT_ROOT / "logs",
    OUTPUT_ROOT / "work",
]:
    directory.mkdir(parents=True, exist_ok=True)

print("Output root:", OUTPUT_ROOT)
print("Validation AnalysisSet:", VALIDATION_ANALYSIS_SET)
print("Guide-quantification ASV:", GUIDE_QUANTIFICATION_ASV)
print("Differential-element ASV:", DIFFERENTIAL_ELEMENT_ASV)
print("Reference file:", REFERENCE_FILE)
```

### Text output (stream)

```text
Output root: /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission
Validation AnalysisSet: IGVFDS4708PCVX
Guide-quantification ASV: /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/
Differential-element ASV: /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/
Reference file: /reference-files/IGVFFI7217ZMJZ/

```

## Markdown cell 16

## 2. Utility functions

## Code cell 17

```python
def detect_delimiter(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(
        path,
        "rt",
        encoding="utf-8-sig",
        errors="replace",
    ) as handle:
        for line in handle:
            if line.strip():
                return "\t" if line.count("\t") >= line.count(",") else ","
    raise ValueError(f"No non-empty line found in {path}")


def read_table(path: Path, *, dtype: object = str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(
        path,
        sep=detect_delimiter(path),
        dtype=dtype,
        keep_default_na=False,
        na_values=[],
        low_memory=False,
        compression="infer",
    ).rename(columns=lambda value: str(value).strip())


def clean_empty_columns(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    removed: list[str] = []
    keep: list[str] = []
    for column in frame.columns:
        values = frame[column].astype(str).str.strip()
        if bool(values.eq("").all()):
            removed.append(column)
        else:
            keep.append(column)
    return frame.loc[:, keep].copy(), removed


def drop_sequential_index_artifacts(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    removed: list[str] = []
    output = frame.copy()
    for column in list(output.columns):
        normalized_name = str(column).strip().lower()
        if not normalized_name.startswith("unnamed:"):
            continue
        numeric = pd.to_numeric(output[column], errors="coerce")
        if numeric.notna().all():
            observed = numeric.to_numpy(dtype=int)
            zero_based = np.arange(len(output))
            one_based = np.arange(1, len(output) + 1)
            if np.array_equal(observed, zero_based) or np.array_equal(
                observed, one_based
            ):
                output = output.drop(columns=[column])
                removed.append(column)
    return output, removed


def md5sum(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_tsv_gzip(
    frame: pd.DataFrame,
    plain_path: Path,
    gzip_path: Path,
) -> str:
    plain_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        plain_path,
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    with plain_path.open("rb") as source, gzip_path.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_output,
            compresslevel=9,
            mtime=0,
        ) as compressed:
            shutil.copyfileobj(source, compressed)

    with gzip.open(gzip_path, "rt", encoding="utf-8") as handle:
        roundtrip = pd.read_csv(
            handle,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            na_values=[],
        )
    if list(roundtrip.columns) != list(frame.columns):
        raise RuntimeError(f"Column mismatch after gzip round-trip: {gzip_path}")
    if len(roundtrip) != len(frame):
        raise RuntimeError(f"Row mismatch after gzip round-trip: {gzip_path}")
    return md5sum(gzip_path)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def choose_one_existing(
    candidates: Sequence[Path],
    label: str,
) -> Path:
    """Select one source, accepting only byte-identical duplicate copies."""
    found: list[Path] = []
    seen: set[str] = set()

    for candidate in candidates:
        candidate = candidate.expanduser()
        if not candidate.is_file():
            continue
        key = str(candidate.resolve())
        if key not in seen:
            found.append(candidate)
            seen.add(key)

    if not found:
        attempted = "\n".join(f"  - {path}" for path in candidates)
        raise FileNotFoundError(
            f"No {label} file found. Checked:\n{attempted}"
        )

    if len(found) == 1:
        return found[0]

    checksum_rows = [
        {
            "path": str(path),
            "size_bytes": path.stat().st_size,
            "md5sum": md5sum(path),
        }
        for path in found
    ]
    unique_checksums = {row["md5sum"] for row in checksum_rows}

    if len(unique_checksums) == 1:
        selected = found[0]
        print(
            f"Found {len(found)} byte-identical copies for {label}; "
            f"selecting the first documented candidate:\n  {selected}"
        )
        for row in checksum_rows:
            print(
                f"  equivalent copy: {row['path']} "
                f"(MD5 {row['md5sum']})"
            )
        return selected

    options = "\n".join(
        f"  - {row['path']}  MD5={row['md5sum']}"
        for row in checksum_rows
    )
    raise RuntimeError(
        f"Multiple non-identical {label} files were found. "
        f"Set an explicit source override before continuing:\n{options}"
    )


def exact_find(
    roots: Sequence[Path],
    names: Sequence[str],
    *,
    maxdepth: int,
) -> list[Path]:
    results: set[Path] = set()
    for root in roots:
        if not root.is_dir():
            continue
        for name in names:
            command = [
                "find",
                str(root),
                "-maxdepth",
                str(maxdepth),
                "-type",
                "f",
                "-name",
                name,
                "-print",
            ]
            completed = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            for line in completed.stdout.splitlines():
                candidate = Path(line.strip())
                if candidate.is_file():
                    results.add(candidate)
    return sorted(results)


def infer_guide_id_column(frame: pd.DataFrame) -> str:
    normalized = {
        str(column).strip().lower().replace(" ", "_"): column
        for column in frame.columns
    }
    for candidate in [
        "guide_id",
        "guide",
        "guide_name",
        "guide_names",
        "sgrna",
        "sg_rna",
        "grna",
    ]:
        if candidate in normalized:
            return normalized[candidate]

    # Prefer a unique column whose values look like guide identifiers.
    scores: list[tuple[float, str]] = []
    guide_pattern = re.compile(
        r"(_[+-]_|_SGTN_|SAFE|NON[-_ ]?TARGET)",
        flags=re.IGNORECASE,
    )
    for column in frame.columns:
        values = frame[column].astype(str).str.strip()
        if values.nunique(dropna=False) != len(frame):
            continue
        score = values.map(lambda value: bool(guide_pattern.search(value))).mean()
        scores.append((float(score), column))
    if scores:
        score, column = max(scores)
        if score >= 0.80:
            return column

    first = frame.columns[0]
    if frame[first].astype(str).str.strip().nunique() == len(frame):
        return first
    raise RuntimeError("Could not identify the guide-ID column")


def infer_target_column(
    frame: pd.DataFrame,
    guide_id_column: str,
) -> str | None:
    normalized = {
        str(column).strip().lower().replace(" ", "_"): column
        for column in frame.columns
    }
    for candidate in ["gene", "target", "target_gene", "label"]:
        if candidate in normalized and normalized[candidate] != guide_id_column:
            return normalized[candidate]
    return None


def parse_source_gene_from_guide_id(guide_id: str) -> str:
    match = re.match(r"^(.*?)(?=_[+-]_|_SGTN_)", str(guide_id))
    if not match:
        raise ValueError(f"Cannot parse source gene from guide ID: {guide_id}")
    return match.group(1)


def promoter_label_from_text(value: str) -> str | None:
    match = re.search(r"(?:_|-)(P1P2|P1|P2)$", str(value).strip())
    return match.group(1) if match else None


def normalize_control_label(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value).upper())


def is_control_id(value: str) -> bool:
    return normalize_control_label(value) in {
        "SAFE",
        "SAFETARGETING",
        "NONTARGETING",
        "NEGATIVECONTROL",
    }


def portal_url(identifier: str, *, frame: str | None = None) -> str:
    """Build an authenticated IGVF database-read URL."""
    value = str(identifier).strip().strip("/")
    query = "?format=json&datastore=database"
    if frame:
        query += f"&frame={frame}"
    return f"https://api.data.igvf.org/{value}/{query}"


def portal_get(
    identifier: str,
    *,
    timeout: int = 30,
    prefer_edit: bool = False,
) -> dict | None:
    """Read one Portal object, with a safe fallback across API frames.

    The ordinary database read is tried first unless `prefer_edit=True`.
    This avoids treating shared, released, or DACC-created objects as
    missing merely because the submitter cannot request an edit frame.
    """
    auth = None
    if (
        os.environ.get("IGVF_API_KEY")
        and os.environ.get("IGVF_SECRET_KEY")
    ):
        auth = (
            os.environ["IGVF_API_KEY"],
            os.environ["IGVF_SECRET_KEY"],
        )

    frames: list[str | None] = (
        ["edit", None, "object"]
        if prefer_edit
        else [None, "object", "edit"]
    )

    attempts: list[dict] = []
    saw_non_404_error = False

    for frame in frames:
        url = portal_url(identifier, frame=frame)
        response = requests.get(
            url,
            auth=auth,
            headers={"Accept": "application/json"},
            timeout=timeout,
        )
        attempts.append(
            {
                "frame": frame or "default",
                "status_code": response.status_code,
                "url": url,
                "response_excerpt": response.text[:300],
            }
        )

        if response.status_code == 200:
            return response.json()
        if response.status_code == 404:
            continue

        saw_non_404_error = True

        # Authentication errors should not be obscured by fallback attempts.
        if response.status_code == 401:
            response.raise_for_status()

        # A shared object may reject frame=edit but remain readable through
        # the ordinary object endpoint. Continue to the next frame.
        if response.status_code in {403, 422}:
            continue

        # Retry other frame-dependent failures through the remaining URLs.
        if response.status_code >= 400:
            continue

    if saw_non_404_error:
        raise RuntimeError(
            "IGVF Portal lookup failed for "
            f"{identifier!r}. Attempts: "
            + json.dumps(attempts, indent=2)
        )
    return None


def infer_dictionary_type(series: pd.Series) -> str:
    values = series.astype(str).str.strip()
    nonempty = values.loc[values.ne("")]
    if nonempty.empty:
        return "string"
    normalized = set(nonempty.str.upper().unique())
    if normalized <= {"TRUE", "FALSE"}:
        return "boolean"
    numeric = pd.to_numeric(nonempty, errors="coerce")
    if numeric.notna().all():
        array = numeric.to_numpy(dtype=float)
        if np.allclose(array, np.round(array)):
            return "integer"
        return "float"
    return "string"


print("Utility functions loaded.")
```

### Text output (stream)

```text
Utility functions loaded.

```

## Markdown cell 18

## 3. Locate and audit the authoritative source files

## Code cell 19

```python
# Explicit overrides take precedence. When an override is provided, automatic
# discovery is not run for that source.

if VALIDATION_GUIDE_COUNTS_OVERRIDE:
    guide_counts_path = choose_one_existing(
        [Path(VALIDATION_GUIDE_COUNTS_OVERRIDE)],
        "validation-screen guide-count override",
    )
else:
    validation_guide_count_candidates = [
        GUIDE_COUNTS_SOURCE,
        Path.cwd() / "guideCounts.csv",
        Path.cwd().parent / "guideCounts.csv",
    ]
    guide_counts_path = choose_one_existing(
        validation_guide_count_candidates,
        "validation-screen guide-count",
    )

if VALIDATION_DIFFERENTIAL_SOURCE_OVERRIDE:
    differential_path = choose_one_existing(
        [Path(VALIDATION_DIFFERENTIAL_SOURCE_OVERRIDE)],
        "validation-screen MAGeCK summary override",
    )
else:
    # Use the standard MAGeCK RRA result, which has the expected
    # id/num/neg|*/pos|* schema. The MLE result in
    # MageckAnalysis/MLE_D3pos_v_D3neg is a different analysis and is
    # intentionally excluded.
    differential_candidates = [
        DIFFERENTIAL_SOURCE,
        Path.cwd() / "D3pos_D3neg_gene_summary.txt",
        Path.cwd().parent / "D3pos_D3neg_gene_summary.txt",
    ]
    differential_candidates.extend(
        exact_find(
            [VALIDATION_ROOT / "MageckAnalysis"],
            ["D3pos_D3neg_gene_summary.txt"],
            maxdepth=4,
        )
    )
    differential_candidates = [
        path
        for path in differential_candidates
        if "MLE_D3pos_v_D3neg" not in path.parts
    ]
    differential_path = choose_one_existing(
        differential_candidates,
        "standard MAGeCK RRA D3 CDH5-positive versus CDH5-negative summary",
    )

if CORRECTED_GUIDE_METADATA_OVERRIDE:
    corrected_guide_metadata_path = choose_one_existing(
        [Path(CORRECTED_GUIDE_METADATA_OVERRIDE)],
        "Jesse-corrected guide metadata override",
    )
else:
    guide_metadata_candidates = [
        Path.cwd() / "300genes_guide_metadata_v43_corrected.tsv",
        Path.cwd() / "300genes_guide_metadata_v43_corrected.tsv.gz",
        Path.cwd().parent / "300genes_guide_metadata_v43_corrected.tsv",
        Path.cwd().parent / "300genes_guide_metadata_v43_corrected.tsv.gz",
        VALIDATION_ROOT / "300genes_guide_metadata_v43_corrected.tsv",
        VALIDATION_ROOT / "300genes_guide_metadata_v43_corrected.tsv.gz",
        Path(
            "/oak/stanford/groups/engreitz/Users/tri/"
            "300genes_guide_metadata_v43_corrected.tsv"
        ),
        Path(
            "/oak/stanford/groups/engreitz/Users/tri/"
            "300genes_guide_metadata_v43_corrected.tsv.gz"
        ),
    ]

    if not any(path.is_file() for path in guide_metadata_candidates):
        guide_metadata_candidates.extend(
            exact_find(
                [
                    Path("/oak/stanford/groups/engreitz/Users/tri"),
                    VALIDATION_ROOT,
                    Path("/oak/stanford/groups/engreitz/Projects/IGVF"),
                ],
                [
                    "300genes_guide_metadata_v43_corrected.tsv",
                    "300genes_guide_metadata_v43_corrected.tsv.gz",
                ],
                maxdepth=8,
            )
        )

    corrected_guide_metadata_path = choose_one_existing(
        guide_metadata_candidates,
        "Jesse-corrected guide metadata",
    )

source_inventory = pd.DataFrame(
    [
        {
            "role": "guide quantifications source",
            "path": str(guide_counts_path),
            "size_bytes": guide_counts_path.stat().st_size,
            "md5sum": md5sum(guide_counts_path),
        },
        {
            "role": "differential element source",
            "path": str(differential_path),
            "size_bytes": differential_path.stat().st_size,
            "md5sum": md5sum(differential_path),
        },
        {
            "role": "authoritative corrected guide metadata",
            "path": str(corrected_guide_metadata_path),
            "size_bytes": corrected_guide_metadata_path.stat().st_size,
            "md5sum": md5sum(corrected_guide_metadata_path),
        },
    ]
)
SOURCE_INVENTORY_PATH = OUTPUT_ROOT / "qc/source_inventory.tsv"
write_tsv(source_inventory, SOURCE_INVENTORY_PATH)

print(source_inventory.to_string(index=False))
print("\nSelected standard validation MAGeCK result:", differential_path)
print("Saved:", SOURCE_INVENTORY_PATH)
```

### Text output (stream)

```text
                                  role                                                                                                                       path  size_bytes                           md5sum
          guide quantifications source                             /oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/UTestAnalysis/guideCounts.csv      113791 e02df8df2f5c896eb460b888e89e0e10
           differential element source /oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/MageckAnalysis/D3pos_v_D3neg/D3pos_D3neg_gene_summary.txt       26233 77054081abb533f385731519bb3205a0
authoritative corrected guide metadata             /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/300genes_guide_metadata_v43_corrected.tsv      288643 25950a1cb64d607345840d1755ab6b74

Selected standard validation MAGeCK result: /oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/MageckAnalysis/D3pos_v_D3neg/D3pos_D3neg_gene_summary.txt
Saved: /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/qc/source_inventory.tsv

```

## Markdown cell 20

## 4. Build the final guide-quantification file

The attached source profile identifies `Guide` as the unique guide-ID
field, `Unnamed: 0` as the target/gene annotation, and six aggregate sample
columns.

Guide identifiers are compared after conservative normalization of case,
Unicode dashes, and surrounding whitespace. The final upload table uses
the exact guide IDs from Jesse's corrected metadata. No underscore/hyphen
substitutions are made, because those could merge distinct identifiers.

The six aggregate count columns are cross-checked against the corresponding
columns in `mageck_sgRNA_readcount.txt`. A missing source count is filled
only when the reference column agrees exactly for all non-missing rows.
Otherwise the notebook stops and writes a discrepancy report.

## Code cell 21

```python
import unicodedata

guide_counts_raw = read_table(guide_counts_path)
guide_counts_raw, removed_empty_guide_columns = clean_empty_columns(
    guide_counts_raw
)
guide_counts_raw, removed_index_columns = drop_sequential_index_artifacts(
    guide_counts_raw
)

# Load Jesse's authoritative corrected guide metadata.
guide_metadata = read_table(corrected_guide_metadata_path)
required_metadata_columns = {
    "guide_id",
    "type",
    "intended_target_name",
    "intended_target_chr",
    "intended_target_start",
    "intended_target_end",
}
missing_metadata_columns = sorted(
    required_metadata_columns - set(guide_metadata.columns)
)
if missing_metadata_columns:
    raise RuntimeError(
        "Corrected guide metadata lacks required columns: "
        + ", ".join(missing_metadata_columns)
    )

guide_metadata["guide_id"] = (
    guide_metadata["guide_id"].astype(str).str.strip()
)
if guide_metadata["guide_id"].eq("").any():
    raise RuntimeError(
        "The corrected guide metadata contains blank guide IDs"
    )
if guide_metadata["guide_id"].duplicated().any():
    raise RuntimeError(
        "The corrected guide metadata contains duplicate guide IDs"
    )

DASH_TRANSLATION = str.maketrans(
    {
        "\u2010": "-",
        "\u2011": "-",
        "\u2012": "-",
        "\u2013": "-",
        "\u2014": "-",
        "\u2212": "-",
    }
)

def normalized_guide_id(value: object) -> str:
    # Conservative comparison key. Underscores and dashes remain distinct.
    text = unicodedata.normalize("NFKC", str(value or ""))
    text = text.translate(DASH_TRANSLATION).strip()
    return text.casefold()

guide_metadata["_guide_key"] = guide_metadata["guide_id"].map(
    normalized_guide_id
)
if guide_metadata["_guide_key"].eq("").any():
    raise RuntimeError(
        "The corrected guide metadata contains a blank normalized guide ID"
    )
if guide_metadata["_guide_key"].duplicated().any():
    collision_path = (
        OUTPUT_ROOT
        / "qc/corrected_guide_metadata_normalized_id_collisions.tsv"
    )
    write_tsv(
        guide_metadata.loc[
            guide_metadata["_guide_key"].duplicated(keep=False),
            ["guide_id", "_guide_key", "type"],
        ].sort_values("_guide_key"),
        collision_path,
    )
    raise RuntimeError(
        "Case/Unicode normalization collapses distinct corrected guide IDs. "
        f"Review {collision_path}"
    )

canonical_id_by_key = dict(
    zip(guide_metadata["_guide_key"], guide_metadata["guide_id"])
)
metadata_keys = set(canonical_id_by_key)
metadata_ids = set(guide_metadata["guide_id"])

# Curated source-to-corrected aliases. They are normalized only for lookup;
# the final output always uses the exact corrected guide ID.
curated_alias_key_map = {
    normalized_guide_id(source_id): normalized_guide_id(canonical_id)
    for source_id, canonical_id in CURATED_GUIDE_ID_ALIASES.items()
}
if len(curated_alias_key_map) != len(CURATED_GUIDE_ID_ALIASES):
    raise RuntimeError(
        "Curated guide aliases contain duplicate normalized source IDs"
    )

for source_key, canonical_key in curated_alias_key_map.items():
    if canonical_key not in metadata_keys:
        raise RuntimeError(
            "A curated alias target is absent from Jesse's corrected "
            f"metadata: {source_key!r} -> {canonical_key!r}"
        )
    if source_key in metadata_keys and source_key != canonical_key:
        raise RuntimeError(
            "A curated source alias is itself a distinct corrected guide; "
            f"refusing to collapse {source_key!r}"
        )

canonical_analysis_group_override: dict[str, str] = {}
for configured_id, analysis_group in (
    CURATED_GUIDE_ANALYSIS_GROUP_OVERRIDES.items()
):
    configured_key = normalized_guide_id(configured_id)
    if configured_key not in metadata_keys:
        raise RuntimeError(
            "An analysis-group override references an unknown corrected "
            f"guide ID: {configured_id}"
        )
    canonical_id = canonical_id_by_key[configured_key]
    canonical_analysis_group_override[canonical_id] = str(
        analysis_group
    ).strip()

def apply_curated_guide_alias(normalized_id: str) -> str:
    return curated_alias_key_map.get(normalized_id, normalized_id)

# Profile every source column.
column_profile_rows: list[dict] = []
for column in guide_counts_raw.columns:
    values = guide_counts_raw[column].astype(str).str.strip()
    numeric = pd.to_numeric(values, errors="coerce")
    raw_normalized_values = values.map(normalized_guide_id)
    normalized_values = raw_normalized_values.map(
        apply_curated_guide_alias
    )
    unique_values = int(values.nunique(dropna=False))
    unique_normalized_values = int(
        normalized_values.nunique(dropna=False)
    )
    normalized_set_match = bool(
        len(values) == len(guide_metadata)
        and unique_normalized_values == len(guide_metadata)
        and set(normalized_values) == metadata_keys
    )
    column_profile_rows.append(
        {
            "source_column": column,
            "rows": len(values),
            "unique_values": unique_values,
            "unique_normalized_values": unique_normalized_values,
            "all_nonempty": bool(values.ne("").all()),
            "missing_values": int(values.eq("").sum()),
            "all_numeric": bool(numeric.notna().all()),
            "integer_numeric": bool(
                numeric.notna().all()
                and np.allclose(
                    numeric.to_numpy(dtype=float),
                    np.round(numeric.to_numpy(dtype=float)),
                )
            ),
            "nonnegative_numeric": bool(
                numeric.notna().all() and (numeric >= 0).all()
            ),
            "exact_corrected_guide_id_set_match": bool(
                len(values) == len(guide_metadata)
                and unique_values == len(guide_metadata)
                and set(values) == metadata_ids
            ),
            "normalized_corrected_guide_id_set_match": (
                normalized_set_match
            ),
            "curated_aliases_applied": int(
                (raw_normalized_values != normalized_values).sum()
            ),
            "example_values": " | ".join(values.head(5).tolist()),
        }
    )

GUIDE_COLUMN_PROFILE_PATH = (
    OUTPUT_ROOT / "qc/validation_guide_count_column_profile.tsv"
)
guide_column_profile = pd.DataFrame(column_profile_rows)
write_tsv(guide_column_profile, GUIDE_COLUMN_PROFILE_PATH)

normalized_column_names = {
    str(column).strip().casefold(): column
    for column in guide_counts_raw.columns
}
normalized_id_candidates = guide_column_profile.loc[
    guide_column_profile[
        "normalized_corrected_guide_id_set_match"
    ].astype(bool),
    "source_column",
].tolist()

# The attached source profile explicitly identifies Guide as the unique ID field.
if "guide" in normalized_column_names:
    selected_guide_id_column = normalized_column_names["guide"]
    if selected_guide_id_column not in normalized_id_candidates:
        mismatch_probe = guide_counts_raw[
            selected_guide_id_column
        ].astype(str).str.strip()
        mismatch_path = (
            OUTPUT_ROOT
            / "qc/validation_guide_ids_unmatched_after_normalization.tsv"
        )
        source_keys_probe = set(
            mismatch_probe
            .map(normalized_guide_id)
            .map(apply_curated_guide_alias)
        )
        mismatch_rows = (
            [
                {
                    "normalized_guide_id": value,
                    "issue": "guideCounts.csv_only",
                }
                for value in sorted(source_keys_probe - metadata_keys)
            ]
            + [
                {
                    "normalized_guide_id": value,
                    "issue": "corrected_metadata_only",
                }
                for value in sorted(metadata_keys - source_keys_probe)
            ]
        )
        write_tsv(pd.DataFrame(mismatch_rows), mismatch_path)
        raise RuntimeError(
            "The explicit Guide column still differs from Jesse's corrected "
            "metadata after case/Unicode normalization. Review "
            f"{mismatch_path}"
        )
elif len(normalized_id_candidates) == 1:
    selected_guide_id_column = normalized_id_candidates[0]
elif len(normalized_id_candidates) > 1:
    preferred_names = {
        "guide_id",
        "guide",
        "guide_name",
        "guide_names",
        "sgrna",
        "grna",
    }
    preferred_candidates = [
        column
        for column in normalized_id_candidates
        if str(column).strip().casefold().replace(" ", "_")
        in preferred_names
    ]
    selected_guide_id_column = (
        preferred_candidates[0]
        if len(preferred_candidates) == 1
        else normalized_id_candidates[0]
    )
else:
    raise RuntimeError(
        "No source column matches Jesse's corrected guide identifiers even "
        "after conservative case/Unicode normalization. Review "
        f"{GUIDE_COLUMN_PROFILE_PATH}"
    )

source_guide_ids = (
    guide_counts_raw[selected_guide_id_column]
    .astype(str)
    .str.strip()
)

source_guide_keys_raw = source_guide_ids.map(normalized_guide_id)
source_guide_keys = source_guide_keys_raw.map(
    apply_curated_guide_alias
)
if source_guide_keys.duplicated().any():
    duplicate_path = (
        OUTPUT_ROOT / "qc/validation_normalized_source_guide_collisions.tsv"
    )
    duplicate_rows = pd.DataFrame(
        {
            "source_guide_id": source_guide_ids,
            "normalized_guide_id": source_guide_keys,
        }
    ).loc[source_guide_keys.duplicated(keep=False)]
    write_tsv(
        duplicate_rows.sort_values("normalized_guide_id"),
        duplicate_path,
    )
    raise RuntimeError(
        "Case/Unicode normalization collapses distinct source guide IDs. "
        f"Review {duplicate_path}"
    )

source_keys = set(source_guide_keys)
source_only_keys = sorted(source_keys - metadata_keys)
metadata_only_keys = sorted(metadata_keys - source_keys)
if source_only_keys or metadata_only_keys:
    mismatch_path = (
        OUTPUT_ROOT
        / "qc/validation_guide_id_set_mismatch_after_normalization.tsv"
    )
    mismatch_rows = (
        [
            {
                "normalized_guide_id": value,
                "issue": "guideCounts.csv_only",
            }
            for value in source_only_keys
        ]
        + [
            {
                "normalized_guide_id": value,
                "issue": "corrected_metadata_only",
            }
            for value in metadata_only_keys
        ]
    )
    write_tsv(pd.DataFrame(mismatch_rows), mismatch_path)
    raise RuntimeError(
        "Guide-ID sets still differ after normalization. Review "
        f"{mismatch_path}"
    )

canonical_guide_ids = source_guide_keys.map(canonical_id_by_key)
ID_NORMALIZATION_AUDIT_PATH = (
    OUTPUT_ROOT / "qc/validation_guide_id_normalization_audit.tsv"
)
id_normalization_audit = pd.DataFrame(
    {
        "source_guide_id": source_guide_ids,
        "source_normalized_guide_id": source_guide_keys_raw,
        "curated_alias_normalized_guide_id": source_guide_keys,
        "canonical_corrected_guide_id": canonical_guide_ids,
        "curated_alias_applied": (
            source_guide_keys_raw != source_guide_keys
        ),
        "identifier_changed": (
            source_guide_ids != canonical_guide_ids
        ),
        "case_only_change": (
            source_guide_ids.str.casefold()
            == canonical_guide_ids.str.casefold()
        )
        & (source_guide_ids != canonical_guide_ids),
    }
)
write_tsv(id_normalization_audit, ID_NORMALIZATION_AUDIT_PATH)

CURATED_GUIDE_ALIAS_AUDIT_PATH = (
    OUTPUT_ROOT / "qc/validation_curated_guide_id_alias_audit.tsv"
)
curated_guide_alias_audit = id_normalization_audit.loc[
    id_normalization_audit["curated_alias_applied"]
].copy()
if not curated_guide_alias_audit.empty:
    curated_guide_alias_audit["reason"] = (
        "the source count/MAGeCK input used the NEDD8 analysis-group "
        "name for the guide that Jesse's corrected metadata names "
        "NEDD8-MDP1"
    )
write_tsv(
    curated_guide_alias_audit,
    CURATED_GUIDE_ALIAS_AUDIT_PATH,
)

# Interpret the documented guideCounts.csv schema.
missing_expected_columns = [
    column
    for column in EXPECTED_SAMPLE_COUNT_COLUMN_NAMES
    if column not in guide_counts_raw.columns
]
if missing_expected_columns:
    raise RuntimeError(
        "guideCounts.csv is missing expected aggregate count columns: "
        + ", ".join(missing_expected_columns)
    )

if "Unnamed: 0" in guide_counts_raw.columns:
    source_target_column = "Unnamed: 0"
elif "gene" in normalized_column_names:
    source_target_column = normalized_column_names["gene"]
else:
    raise RuntimeError(
        "Could not identify the target/gene annotation column. The attached "
        "profile showed it as `Unnamed: 0`."
    )

expected_source_columns = {
    selected_guide_id_column,
    source_target_column,
    *EXPECTED_SAMPLE_COUNT_COLUMN_NAMES,
}
duplicate_guide_id_columns: list[str] = []
annotation_columns = [source_target_column]
excluded_annotation_columns = [
    column
    for column in guide_counts_raw.columns
    if column not in expected_source_columns
]

working = pd.DataFrame(
    {
        "guide_id": canonical_guide_ids,
        "gene": (
            guide_counts_raw[source_target_column]
            .astype(str)
            .str.strip()
        ),
    }
)
for column in EXPECTED_SAMPLE_COUNT_COLUMN_NAMES:
    working[column] = (
        guide_counts_raw[column].astype(str).str.strip()
    )

# Validate the target annotation against targeting guide IDs.
type_by_id = dict(
    zip(guide_metadata["guide_id"], guide_metadata["type"])
)
working["_guide_type"] = working["guide_id"].map(type_by_id)

def parsed_source_gene_or_blank(guide_id: str) -> str:
    try:
        return parse_source_gene_from_guide_id(guide_id)
    except Exception:
        return ""

working["_parsed_gene"] = working["guide_id"].map(
    parsed_source_gene_or_blank
)
working["_expected_analysis_group"] = working.apply(
    lambda row: canonical_analysis_group_override.get(
        row["guide_id"],
        row["_parsed_gene"],
    ),
    axis=1,
)
working["_analysis_group_override_applied"] = (
    working["guide_id"].isin(canonical_analysis_group_override)
)
targeting_rows = (
    working["_guide_type"]
    .astype(str)
    .str.strip()
    .str.casefold()
    .eq("targeting")
)
target_match = (
    working["gene"].str.casefold()
    == working["_expected_analysis_group"].str.casefold()
)
ANNOTATION_AUDIT_PATH = (
    OUTPUT_ROOT / "qc/validation_guide_count_annotation_audit.tsv"
)
annotation_audit = working.loc[
    :,
    [
        "guide_id",
        "gene",
        "_parsed_gene",
        "_expected_analysis_group",
        "_analysis_group_override_applied",
        "_guide_type",
    ],
].copy()
annotation_audit["target_label_matches_expected_analysis_group"] = (
    ~targeting_rows | target_match
)
write_tsv(annotation_audit, ANNOTATION_AUDIT_PATH)

ANALYSIS_GROUP_OVERRIDE_AUDIT_PATH = (
    OUTPUT_ROOT
    / "qc/validation_curated_analysis_group_override_audit.tsv"
)
write_tsv(
    annotation_audit.loc[
        annotation_audit["_analysis_group_override_applied"]
    ].copy(),
    ANALYSIS_GROUP_OVERRIDE_AUDIT_PATH,
)

target_mismatch = working.loc[
    targeting_rows & ~target_match,
    [
        "guide_id",
        "gene",
        "_parsed_gene",
        "_expected_analysis_group",
        "_analysis_group_override_applied",
        "_guide_type",
    ],
].copy()
if not target_mismatch.empty:
    target_mismatch_path = (
        OUTPUT_ROOT / "qc/validation_target_label_mismatches.tsv"
    )
    write_tsv(target_mismatch, target_mismatch_path)
    raise RuntimeError(
        "The source gene annotation disagrees with targeting guide IDs. "
        f"Review {target_mismatch_path}"
    )

EXCLUDED_ANNOTATIONS_PATH = (
    OUTPUT_ROOT
    / "qc/validation_guide_count_excluded_source_annotations.tsv"
)
if excluded_annotation_columns:
    write_tsv(
        guide_counts_raw.loc[
            :,
            [selected_guide_id_column, *excluded_annotation_columns],
        ],
        EXCLUDED_ANNOTATIONS_PATH,
    )
else:
    write_tsv(
        pd.DataFrame(
            columns=["guide_id", "excluded_source_column", "value"]
        ),
        EXCLUDED_ANNOTATIONS_PATH,
    )

# Cross-check and repair the six aggregate columns from the MAGeCK count table.
if not MAGECK_GUIDE_COUNTS_REFERENCE.is_file():
    raise FileNotFoundError(
        "The MAGeCK guide-count reference is required because the attached "
        "profile shows an incomplete D3pos1 column: "
        f"{MAGECK_GUIDE_COUNTS_REFERENCE}"
    )

reference_counts = read_table(MAGECK_GUIDE_COUNTS_REFERENCE)
required_reference_columns = {
    "sgRNA",
    "Gene",
    *VALIDATION_REFERENCE_COUNT_COLUMN_MAP.values(),
}
missing_reference_columns = sorted(
    required_reference_columns - set(reference_counts.columns)
)
if missing_reference_columns:
    raise RuntimeError(
        "mageck_sgRNA_readcount.txt lacks required columns: "
        + ", ".join(missing_reference_columns)
    )

reference_counts["_guide_key_raw"] = (
    reference_counts["sgRNA"]
    .astype(str)
    .str.strip()
    .map(normalized_guide_id)
)
reference_counts["_guide_key"] = (
    reference_counts["_guide_key_raw"]
    .map(apply_curated_guide_alias)
)
if reference_counts["_guide_key"].duplicated().any():
    duplicate_reference_path = (
        OUTPUT_ROOT
        / "qc/validation_reference_normalized_guide_collisions.tsv"
    )
    write_tsv(
        reference_counts.loc[
            reference_counts["_guide_key"].duplicated(keep=False)
        ].sort_values("_guide_key"),
        duplicate_reference_path,
    )
    raise RuntimeError(
        "The MAGeCK reference contains duplicate normalized guide IDs. "
        f"Review {duplicate_reference_path}"
    )

reference_keys = set(reference_counts["_guide_key"])
if reference_keys != metadata_keys:
    reference_mismatch_path = (
        OUTPUT_ROOT / "qc/validation_reference_guide_id_mismatch.tsv"
    )
    reference_mismatch_rows = (
        [
            {
                "normalized_guide_id": value,
                "issue": "mageck_reference_only",
            }
            for value in sorted(reference_keys - metadata_keys)
        ]
        + [
            {
                "normalized_guide_id": value,
                "issue": "corrected_metadata_only",
            }
            for value in sorted(metadata_keys - reference_keys)
        ]
    )
    write_tsv(
        pd.DataFrame(reference_mismatch_rows),
        reference_mismatch_path,
    )
    raise RuntimeError(
        "The MAGeCK reference guide set differs from the corrected metadata. "
        f"Review {reference_mismatch_path}"
    )

aligned_reference = (
    reference_counts.set_index("_guide_key")
    .loc[source_guide_keys.tolist()]
    .reset_index(drop=True)
)

count_columns = list(EXPECTED_SAMPLE_COUNT_COLUMN_NAMES)
count_repair_rows: list[dict] = []
count_mapping_rows: list[dict] = []
count_discrepancy_rows: list[dict] = []

for source_column, reference_column in (
    VALIDATION_REFERENCE_COUNT_COLUMN_MAP.items()
):
    source_text = working[source_column].astype(str).str.strip()
    source_numeric = pd.to_numeric(source_text, errors="coerce")
    reference_numeric = pd.to_numeric(
        aligned_reference[reference_column],
        errors="coerce",
    )

    if reference_numeric.isna().any():
        raise RuntimeError(
            f"Reference column {reference_column} has missing/non-numeric "
            "values"
        )
    if (reference_numeric < 0).any():
        raise RuntimeError(
            f"Reference column {reference_column} has negative values"
        )
    if not np.allclose(
        reference_numeric.to_numpy(dtype=float),
        np.round(reference_numeric.to_numpy(dtype=float)),
    ):
        raise RuntimeError(
            f"Reference column {reference_column} has non-integer values"
        )

    source_present = source_text.ne("") & source_numeric.notna()
    source_invalid_nonblank = source_text.ne("") & source_numeric.isna()
    if source_invalid_nonblank.any():
        invalid_path = (
            OUTPUT_ROOT
            / "qc"
            / f"validation_non_numeric_{source_column}.tsv"
        )
        write_tsv(
            pd.DataFrame(
                {
                    "guide_id": working.loc[
                        source_invalid_nonblank, "guide_id"
                    ],
                    "source_value": source_text.loc[
                        source_invalid_nonblank
                    ],
                }
            ),
            invalid_path,
        )
        raise RuntimeError(
            f"Non-numeric non-empty values found in {source_column}; "
            f"review {invalid_path}"
        )

    unequal_present = source_present & (
        source_numeric != reference_numeric
    )
    for row_index in np.flatnonzero(unequal_present.to_numpy()):
        count_discrepancy_rows.append(
            {
                "guide_id": working.iloc[row_index]["guide_id"],
                "source_column": source_column,
                "reference_column": reference_column,
                "source_value": source_numeric.iloc[row_index],
                "reference_value": reference_numeric.iloc[row_index],
            }
        )

    missing_source = ~source_present
    for row_index in np.flatnonzero(missing_source.to_numpy()):
        count_repair_rows.append(
            {
                "guide_id": working.iloc[row_index]["guide_id"],
                "source_column": source_column,
                "reference_column": reference_column,
                "original_value": source_text.iloc[row_index],
                "filled_value": int(reference_numeric.iloc[row_index]),
                "repair_rule": (
                    "filled from MAGeCK reference after exact agreement at "
                    "all non-missing rows"
                ),
            }
        )

    count_mapping_rows.append(
        {
            "source_column": source_column,
            "reference_column": reference_column,
            "nonmissing_source_rows": int(source_present.sum()),
            "missing_source_rows": int(missing_source.sum()),
            "disagreeing_nonmissing_rows": int(
                unequal_present.sum()
            ),
        }
    )

    if not unequal_present.any():
        repaired = source_numeric.copy()
        repaired.loc[missing_source] = reference_numeric.loc[missing_source]
        working[source_column] = repaired.round().astype("int64")

COUNT_COLUMN_MAPPING_PATH = (
    OUTPUT_ROOT
    / "qc/validation_guide_count_reference_column_mapping.tsv"
)
write_tsv(pd.DataFrame(count_mapping_rows), COUNT_COLUMN_MAPPING_PATH)

COUNT_REPAIR_AUDIT_PATH = (
    OUTPUT_ROOT
    / "qc/validation_guide_count_missing_value_repairs.tsv"
)
write_tsv(pd.DataFrame(count_repair_rows), COUNT_REPAIR_AUDIT_PATH)

COUNT_DISCREPANCY_PATH = (
    OUTPUT_ROOT
    / "qc/validation_guide_count_reference_discrepancies.tsv"
)
write_tsv(
    pd.DataFrame(count_discrepancy_rows),
    COUNT_DISCREPANCY_PATH,
)

if count_discrepancy_rows:
    raise RuntimeError(
        "guideCounts.csv disagrees with mageck_sgRNA_readcount.txt at "
        "non-missing rows. No final file was written. Review "
        f"{COUNT_DISCREPANCY_PATH}"
    )

guide_counts = working.loc[
    :,
    ["guide_id", "gene", *count_columns],
].copy()

if guide_counts["guide_id"].duplicated().any():
    raise RuntimeError(
        "Canonicalized validation guide IDs are not unique"
    )
for column in count_columns:
    numeric = pd.to_numeric(guide_counts[column], errors="coerce")
    if numeric.isna().any():
        raise RuntimeError(
            f"Missing or non-numeric values remain in {column}"
        )
    if (numeric < 0).any():
        raise RuntimeError(f"Negative values remain in {column}")
    values = numeric.to_numpy(dtype=float)
    if not np.allclose(values, np.round(values)):
        raise RuntimeError(
            f"Non-integer values remain in {column}"
        )
    guide_counts[column] = numeric.round().astype("int64")

if STRICT_EXPECTED_DIMENSIONS:
    if len(guide_counts) != EXPECTED_GUIDES:
        raise RuntimeError(
            f"Expected {EXPECTED_GUIDES} guides; found "
            f"{len(guide_counts)}"
        )
    if len(count_columns) != EXPECTED_SAMPLE_COUNT_COLUMNS:
        raise RuntimeError(
            f"Expected {EXPECTED_SAMPLE_COUNT_COLUMNS} aggregate count "
            f"columns; found {len(count_columns)}"
        )

count_ids = set(guide_counts["guide_id"])
if count_ids != metadata_ids:
    raise RuntimeError(
        "Internal error: canonicalized output guide IDs do not exactly "
        "match corrected guide metadata"
    )

selected_guide_id_column = str(selected_guide_id_column)
target_column = "gene"

reference_crosscheck = {
    "reference_file_present": True,
    "reference_file": str(MAGECK_GUIDE_COUNTS_REFERENCE),
    "normalized_guide_id_sets_equal": reference_keys == metadata_keys,
    "source_to_reference_column_map": (
        VALIDATION_REFERENCE_COUNT_COLUMN_MAP
    ),
    "missing_source_values_filled": len(count_repair_rows),
    "nonmissing_value_discrepancies": len(
        count_discrepancy_rows
    ),
    "source_ids_changed_to_canonical_form": int(
        id_normalization_audit["identifier_changed"].sum()
    ),
}

GUIDE_PLAIN = (
    OUTPUT_ROOT / "files/validation_screen_guide_quantifications.tsv"
)
GUIDE_GZ = (
    OUTPUT_ROOT
    / "files/validation_screen_guide_quantifications.tsv.gz"
)
GUIDE_MD5 = deterministic_tsv_gzip(
    guide_counts,
    GUIDE_PLAIN,
    GUIDE_GZ,
)

sample_manifest = pd.DataFrame(
    [
        {
            "sample_column": "D01",
            "reference_column": "D0R1",
            "timepoint": "D0",
            "fraction": "baseline",
            "biological_replicate": "1",
            "data_type": "non-negative integer read count",
        },
        {
            "sample_column": "D02",
            "reference_column": "D0R2",
            "timepoint": "D0",
            "fraction": "baseline",
            "biological_replicate": "2",
            "data_type": "non-negative integer read count",
        },
        {
            "sample_column": "D3neg1",
            "reference_column": "D3NR1",
            "timepoint": "D3",
            "fraction": "CDH5-negative",
            "biological_replicate": "1",
            "data_type": "non-negative integer read count",
        },
        {
            "sample_column": "D3neg2",
            "reference_column": "D3NR2",
            "timepoint": "D3",
            "fraction": "CDH5-negative",
            "biological_replicate": "2",
            "data_type": "non-negative integer read count",
        },
        {
            "sample_column": "D3pos1",
            "reference_column": "D3PR1",
            "timepoint": "D3",
            "fraction": "CDH5-positive",
            "biological_replicate": "1",
            "data_type": "non-negative integer read count",
        },
        {
            "sample_column": "D3pos2",
            "reference_column": "D3PR2",
            "timepoint": "D3",
            "fraction": "CDH5-positive",
            "biological_replicate": "2",
            "data_type": "non-negative integer read count",
        },
    ]
)
SAMPLE_MANIFEST_PATH = (
    OUTPUT_ROOT / "documentation/guide_count_sample_columns.tsv"
)
write_tsv(sample_manifest, SAMPLE_MANIFEST_PATH)

print(f"Guide rows: {len(guide_counts):,}")
print(
    "Selected source guide-ID column:",
    selected_guide_id_column,
)
print("Target annotation column:", source_target_column)
print("Count columns:", count_columns)
print(
    "Guide IDs rewritten to exact corrected-metadata form:",
    int(id_normalization_audit["identifier_changed"].sum()),
)
print(
    "Curated NEDD8/NEDD8-MDP1 guide aliases applied:",
    int(id_normalization_audit["curated_alias_applied"].sum()),
)
print("Curated guide alias audit:", CURATED_GUIDE_ALIAS_AUDIT_PATH)
print(
    "Curated analysis-group override audit:",
    ANALYSIS_GROUP_OVERRIDE_AUDIT_PATH,
)
print(
    "Missing count values filled from exact-matching MAGeCK columns:",
    len(count_repair_rows),
)
print("Reference mapping:")
print(pd.DataFrame(count_mapping_rows).to_string(index=False))
print("Output:", GUIDE_GZ)
print("MD5:", GUIDE_MD5)
```

### Text output (stream)

```text
Guide rows: 1,926
Selected source guide-ID column: Guide
Target annotation column: Unnamed: 0
Count columns: ['D01', 'D02', 'D3neg1', 'D3neg2', 'D3pos1', 'D3pos2']
Guide IDs rewritten to exact corrected-metadata form: 91
Curated NEDD8/NEDD8-MDP1 guide aliases applied: 1
Curated guide alias audit: /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/qc/validation_curated_guide_id_alias_audit.tsv
Curated analysis-group override audit: /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/qc/validation_curated_analysis_group_override_audit.tsv
Missing count values filled from exact-matching MAGeCK columns: 2
Reference mapping:
source_column reference_column  nonmissing_source_rows  missing_source_rows  disagreeing_nonmissing_rows
          D01             D0R1                    1926                    0                            0
          D02             D0R2                    1926                    0                            0
       D3neg1            D3NR1                    1926                    0                            0
       D3neg2            D3NR2                    1926                    0                            0
       D3pos1            D3PR1                    1924                    2                            0
       D3pos2            D3PR2                    1926                    0                            0
Output: /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_guide_quantifications.tsv.gz
MD5: 0d4a1e607ef5f0da24f7fffa59e0d00c

```

## Markdown cell 22

## 5. Build the promoter-element reference

The notebook obtains canonical promoter labels from Jesse's corrected file in this order:

1. `description`, when it ends in `P1`, `P2`, or `P1P2`;
2. `putative_target_genes`, when it ends in one of those labels;
3. the guide ID suffix.

`_SGTN_` guides are treated as supplemental guides. They are merged into the gene's sole canonical promoter and are not emitted as separate alt-contig elements.

## Code cell 23

```python
targeting_guides = guide_metadata.loc[
    guide_metadata["type"].astype(str).str.strip().str.lower().eq("targeting")
].copy()

# Current corrected metadata should use versionless ENSG identifiers.
invalid_target_ids = targeting_guides.loc[
    ~targeting_guides["intended_target_name"]
    .astype(str)
    .str.fullmatch(r"ENSG\d+")
]
if not invalid_target_ids.empty:
    invalid_path = (
        OUTPUT_ROOT / "qc/invalid_corrected_intended_target_names.tsv"
    )
    write_tsv(invalid_target_ids, invalid_path)
    raise RuntimeError(
        "The selected corrected guide metadata has non-ENSG "
        f"`intended_target_name` values. Review {invalid_path}"
    )

for coordinate_column in [
    "intended_target_start",
    "intended_target_end",
]:
    numeric = pd.to_numeric(
        targeting_guides[coordinate_column],
        errors="coerce",
    )
    if numeric.isna().any():
        raise RuntimeError(
            f"Non-numeric {coordinate_column} values were found"
        )
    targeting_guides[coordinate_column] = numeric.astype("int64")

if (
    targeting_guides["intended_target_end"]
    < targeting_guides["intended_target_start"]
).any():
    raise RuntimeError(
        "At least one corrected intended-target interval has end < start"
    )

targeting_guides["source_gene"] = targeting_guides["guide_id"].map(
    parse_source_gene_from_guide_id
)
targeting_guides["is_sgtn"] = targeting_guides["guide_id"].str.contains(
    "_SGTN_",
    regex=False,
)
targeting_guides["analysis_group_override"] = (
    targeting_guides["guide_id"]
    .map(canonical_analysis_group_override)
    .fillna("")
)
targeting_guides["is_curated_analysis_group_alias"] = (
    targeting_guides["analysis_group_override"].ne("")
)

def extract_element_label(row: pd.Series) -> str | None:
    for column in ["description", "putative_target_genes", "guide_id"]:
        if column in row.index:
            label = promoter_label_from_text(str(row[column]))
            if label is not None:
                return label
    return None

canonical_guides = targeting_guides.loc[
    ~targeting_guides["is_sgtn"]
    & ~targeting_guides["is_curated_analysis_group_alias"]
].copy()
canonical_guides["promoter_label"] = canonical_guides.apply(
    extract_element_label,
    axis=1,
)

missing_labels = canonical_guides.loc[
    canonical_guides["promoter_label"].isna()
]
if not missing_labels.empty:
    missing_label_path = (
        OUTPUT_ROOT / "qc/corrected_guides_missing_promoter_label.tsv"
    )
    write_tsv(missing_labels, missing_label_path)
    raise RuntimeError(
        "Canonical guides lack a P1/P2/P1P2 element label. "
        f"Review {missing_label_path}"
    )

canonical_guides["gene_symbol"] = canonical_guides["source_gene"]
canonical_guides["element_id"] = (
    canonical_guides["gene_symbol"]
    + "_-"
    + canonical_guides["promoter_label"]
)

consistency_fields = [
    "source_gene",
    "gene_symbol",
    "promoter_label",
    "intended_target_name",
    "intended_target_chr",
    "intended_target_start",
    "intended_target_end",
]
variation = (
    canonical_guides.groupby("element_id")[consistency_fields]
    .nunique(dropna=False)
)
inconsistent_elements = variation.loc[(variation > 1).any(axis=1)]
if not inconsistent_elements.empty:
    inconsistent_path = (
        OUTPUT_ROOT / "qc/inconsistent_corrected_element_annotations.tsv"
    )
    write_tsv(inconsistent_elements.reset_index(), inconsistent_path)
    raise RuntimeError(
        "Element annotations are not internally consistent. "
        f"Review {inconsistent_path}"
    )

element_rows: list[dict] = []
for element_id, group in canonical_guides.groupby(
    "element_id",
    sort=True,
):
    row = group.iloc[0]
    element_rows.append(
        {
            "source_gene": row["source_gene"],
            "gene_symbol": row["gene_symbol"],
            "element_id": element_id,
            "promoter_label": row["promoter_label"],
            "ensembl_gene_id": row["intended_target_name"],
            "element_chr": row["intended_target_chr"],
            "element_start": int(row["intended_target_start"]),
            "element_end": int(row["intended_target_end"]),
            "canonical_guide_count": int(len(group)),
        }
    )

element_reference = pd.DataFrame(element_rows)

elements_per_gene = element_reference.groupby("source_gene")[
    "element_id"
].nunique()
unexpected_multi = elements_per_gene.loc[
    (elements_per_gene > 1)
    & (~elements_per_gene.index.str.upper().isin(ALLOWED_MULTI_ELEMENT_GENES))
]
if not unexpected_multi.empty:
    unexpected_path = (
        OUTPUT_ROOT / "qc/unexpected_multi_promoter_genes.tsv"
    )
    write_tsv(
        unexpected_multi.rename("number_of_elements")
        .reset_index(),
        unexpected_path,
    )
    raise RuntimeError(
        "Genes other than NOTCH1 map to multiple canonical promoter "
        f"elements. Review {unexpected_path}"
    )

# Merge supplemental SGTN guides into the sole canonical promoter for each
# source gene. Stop if the target is ambiguous.
sgtn_counts = (
    targeting_guides.loc[targeting_guides["is_sgtn"]]
    .groupby("source_gene")
    .size()
    .to_dict()
)

element_reference["supplemental_sgtn_guide_count"] = 0
sgtn_audit_rows: list[dict] = []
for source_gene, supplemental_count in sorted(sgtn_counts.items()):
    matching = element_reference.loc[
        element_reference["source_gene"].eq(source_gene)
    ]
    if len(matching) != 1:
        raise RuntimeError(
            f"Cannot merge {supplemental_count} SGTN guides for "
            f"{source_gene}; canonical elements are "
            f"{matching['element_id'].tolist()}"
        )
    index = matching.index[0]
    element_reference.loc[
        index,
        "supplemental_sgtn_guide_count",
    ] = int(supplemental_count)
    sgtn_audit_rows.append(
        {
            "source_gene": source_gene,
            "element_id": matching.iloc[0]["element_id"],
            "supplemental_sgtn_guide_count": int(supplemental_count),
            "rule": "merged into the gene's sole canonical promoter",
        }
    )

element_reference["supplemental_sgtn_guide_count"] = (
    element_reference["supplemental_sgtn_guide_count"].astype(int)
)

# Merge guides with an explicit native-analysis grouping override into the
# destination element used by that analysis. This currently applies only to
# the NEDD8-MDP1 guide that the source MAGeCK input grouped under NEDD8.
element_reference["curated_analysis_alias_guide_count"] = 0
analysis_alias_merge_rows: list[dict] = []
analysis_alias_guides = targeting_guides.loc[
    targeting_guides["is_curated_analysis_group_alias"]
].copy()

for _, alias_guide in analysis_alias_guides.iterrows():
    analysis_group = str(
        alias_guide["analysis_group_override"]
    ).strip()
    matching = element_reference.loc[
        element_reference["source_gene"]
        .astype(str)
        .str.casefold()
        .eq(analysis_group.casefold())
    ]
    if len(matching) != 1:
        raise RuntimeError(
            "Cannot merge curated analysis-group guide "
            f"{alias_guide['guide_id']} into {analysis_group}; "
            f"matching elements are {matching['element_id'].tolist()}"
        )

    element_index = matching.index[0]
    element_reference.loc[
        element_index,
        "curated_analysis_alias_guide_count",
    ] += 1
    analysis_alias_merge_rows.append(
        {
            "corrected_guide_id": alias_guide["guide_id"],
            "corrected_source_gene": alias_guide["source_gene"],
            "analysis_group": analysis_group,
            "destination_element_id": matching.iloc[0]["element_id"],
            "corrected_intended_target_name": (
                alias_guide["intended_target_name"]
            ),
            "reason": (
                "native validation guideCounts/MAGeCK input used the "
                "NEDD8 analysis-group name for this exact guide"
            ),
        }
    )

element_reference["curated_analysis_alias_guide_count"] = (
    element_reference["curated_analysis_alias_guide_count"].astype(int)
)
element_reference["total_guide_count"] = (
    element_reference["canonical_guide_count"]
    + element_reference["supplemental_sgtn_guide_count"]
    + element_reference["curated_analysis_alias_guide_count"]
)

ELEMENT_REFERENCE_PATH = (
    OUTPUT_ROOT
    / "documentation/validation_screen_promoter_element_reference.tsv"
)
write_tsv(element_reference, ELEMENT_REFERENCE_PATH)

SGTN_AUDIT_PATH = OUTPUT_ROOT / "qc/supplemental_sgtn_merge_audit.tsv"
write_tsv(pd.DataFrame(sgtn_audit_rows), SGTN_AUDIT_PATH)

ANALYSIS_ALIAS_MERGE_AUDIT_PATH = (
    OUTPUT_ROOT
    / "qc/validation_curated_analysis_alias_element_merge_audit.tsv"
)
write_tsv(
    pd.DataFrame(analysis_alias_merge_rows),
    ANALYSIS_ALIAS_MERGE_AUDIT_PATH,
)

print(f"Canonical promoter elements: {len(element_reference):,}")
print("Multi-promoter genes:")
print(
    elements_per_gene.loc[elements_per_gene > 1].to_string()
    if (elements_per_gene > 1).any()
    else "none"
)
print(f"Supplemental SGTN guides merged: {sum(sgtn_counts.values())}")
print(
    "Curated analysis-group alias guides merged:",
    len(analysis_alias_merge_rows),
)
print(
    "Curated alias element-merge audit:",
    ANALYSIS_ALIAS_MERGE_AUDIT_PATH,
)
print("Element reference:", ELEMENT_REFERENCE_PATH)
```

### Text output (stream)

```text
Canonical promoter elements: 298
Multi-promoter genes:
source_gene
NOTCH1    2
Supplemental SGTN guides merged: 24
Curated analysis-group alias guides merged: 1
Curated alias element-merge audit: /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/qc/validation_curated_analysis_alias_element_merge_audit.tsv
Element reference: /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/documentation/validation_screen_promoter_element_reference.tsv

```

## Markdown cell 24

## 6. Build the final differential-element-quantification file

Each native MAGeCK biological row is mapped to its promoter element. The one native `NOTCH1` row becomes two rows:

```text
NOTCH1_-P1
NOTCH1_-P2
```

All native MAGeCK statistics are copied unchanged into both rows. The added `element_guide_count` field records the actual number of guides assigned to each promoter.

## Code cell 25

```python
differential_raw = read_table(differential_path)
differential_raw, removed_empty_differential_columns = clean_empty_columns(
    differential_raw
)

required_mageck_columns = [
    "id",
    "num",
    "neg|score",
    "neg|p-value",
    "neg|fdr",
    "neg|rank",
    "neg|goodsgrna",
    "neg|lfc",
    "pos|score",
    "pos|p-value",
    "pos|fdr",
    "pos|rank",
    "pos|goodsgrna",
    "pos|lfc",
]
missing_mageck_columns = [
    column
    for column in required_mageck_columns
    if column not in differential_raw.columns
]
if missing_mageck_columns:
    raise RuntimeError(
        "MAGeCK summary lacks required columns: "
        + ", ".join(missing_mageck_columns)
    )

differential_raw["id"] = differential_raw["id"].astype(str).str.strip()
if differential_raw["id"].duplicated().any():
    raise RuntimeError("The native MAGeCK summary has duplicate IDs")
if (
    STRICT_EXPECTED_DIMENSIONS
    and len(differential_raw) != EXPECTED_MAGECK_ROWS
):
    raise RuntimeError(
        f"Expected {EXPECTED_MAGECK_ROWS} native MAGeCK rows; "
        f"found {len(differential_raw)}"
    )

integer_fields = {
    "num",
    "neg|rank",
    "neg|goodsgrna",
    "pos|rank",
    "pos|goodsgrna",
}
probability_fields = {
    "neg|score",
    "neg|p-value",
    "neg|fdr",
    "pos|score",
    "pos|p-value",
    "pos|fdr",
}

for column in required_mageck_columns[1:]:
    numeric = pd.to_numeric(differential_raw[column], errors="coerce")
    if numeric.isna().any():
        bad_path = (
            OUTPUT_ROOT
            / "qc"
            / f"non_numeric_mageck_{re.sub(r'[^A-Za-z0-9]+', '_', column)}.tsv"
        )
        write_tsv(differential_raw.loc[numeric.isna()], bad_path)
        raise RuntimeError(
            f"Non-numeric values found in {column}; review {bad_path}"
        )
    if column in integer_fields:
        values = numeric.to_numpy(dtype=float)
        if not np.allclose(values, np.round(values)):
            raise RuntimeError(
                f"Non-integer values found in MAGeCK field {column}"
            )
    if column in probability_fields:
        if ((numeric < 0) | (numeric > 1)).any():
            raise RuntimeError(
                f"Values outside [0, 1] found in {column}"
            )

elements_by_source = {
    source_gene.upper(): group.sort_values("element_id").reset_index(
        drop=True
    )
    for source_gene, group in element_reference.groupby(
        "source_gene",
        sort=False,
    )
}

annotation_fields = [
    "id",
    "source_mageck_id",
    "gene_symbol",
    "ensembl_gene_id",
    "element_chr",
    "element_start",
    "element_end",
    "element_guide_count",
    "statistics_duplicated_across_elements",
]
native_stat_fields = [
    column for column in differential_raw.columns if column != "id"
]

output_rows: list[dict] = []
mapping_audit_rows: list[dict] = []
unmapped_ids: list[str] = []

for zero_based_index, source_row in differential_raw.iterrows():
    source_id = str(source_row["id"]).strip()

    if is_control_id(source_id):
        output = {
            "id": source_id,
            "source_mageck_id": source_id,
            "gene_symbol": "",
            "ensembl_gene_id": "",
            "element_chr": "",
            "element_start": "",
            "element_end": "",
            "element_guide_count": "",
            "statistics_duplicated_across_elements": "FALSE",
        }
        output.update(
            {
                column: source_row[column]
                for column in native_stat_fields
            }
        )
        output_rows.append(output)
        mapping_audit_rows.append(
            {
                "source_row_number": int(zero_based_index) + 2,
                "source_mageck_id": source_id,
                "output_element_id": source_id,
                "action": "aggregate control row preserved",
            }
        )
        continue

    elements = elements_by_source.get(source_id.upper())
    if elements is None:
        unmapped_ids.append(source_id)
        continue

    duplicated = len(elements) > 1
    if duplicated and source_id.upper() not in ALLOWED_MULTI_ELEMENT_GENES:
        raise RuntimeError(
            f"Unexpected multi-promoter source target {source_id}: "
            f"{elements['element_id'].tolist()}"
        )

    for _, element in elements.iterrows():
        output = {
            "id": element["element_id"],
            "source_mageck_id": source_id,
            "gene_symbol": element["gene_symbol"],
            "ensembl_gene_id": element["ensembl_gene_id"],
            "element_chr": element["element_chr"],
            "element_start": int(element["element_start"]),
            "element_end": int(element["element_end"]),
            "element_guide_count": int(element["total_guide_count"]),
            "statistics_duplicated_across_elements": (
                "TRUE" if duplicated else "FALSE"
            ),
        }
        output.update(
            {
                column: source_row[column]
                for column in native_stat_fields
            }
        )
        output_rows.append(output)
        mapping_audit_rows.append(
            {
                "source_row_number": int(zero_based_index) + 2,
                "source_mageck_id": source_id,
                "output_element_id": element["element_id"],
                "gene_symbol": element["gene_symbol"],
                "ensembl_gene_id": element["ensembl_gene_id"],
                "element_coordinate": (
                    f"{element['element_chr']}:"
                    f"{int(element['element_start'])}-"
                    f"{int(element['element_end'])}"
                ),
                "element_guide_count": int(element["total_guide_count"]),
                "source_num": source_row["num"],
                "action": (
                    "duplicated unchanged gene-level statistics across "
                    "promoter elements"
                    if duplicated
                    else "mapped one gene-level row to its sole promoter"
                ),
            }
        )

if unmapped_ids:
    unmapped_path = OUTPUT_ROOT / "qc/unmapped_mageck_target_ids.tsv"
    write_tsv(
        pd.DataFrame({"unmapped_mageck_id": unmapped_ids}),
        unmapped_path,
    )
    raise RuntimeError(
        "Native MAGeCK IDs could not be mapped to the corrected guide "
        f"metadata. Review {unmapped_path}"
    )

differential_elements = pd.DataFrame(
    output_rows,
    columns=annotation_fields + native_stat_fields,
)

if differential_elements["id"].duplicated().any():
    duplicate_path = (
        OUTPUT_ROOT / "qc/duplicate_differential_element_ids.tsv"
    )
    write_tsv(
        differential_elements.loc[
            differential_elements["id"].duplicated(keep=False)
        ],
        duplicate_path,
    )
    raise RuntimeError(
        f"Duplicate output element IDs found; review {duplicate_path}"
    )

biological_rows = differential_elements.loc[
    ~differential_elements["source_mageck_id"].map(is_control_id)
].copy()

if (
    STRICT_EXPECTED_DIMENSIONS
    and len(biological_rows) != EXPECTED_BIOLOGICAL_ELEMENT_ROWS
):
    raise RuntimeError(
        f"Expected {EXPECTED_BIOLOGICAL_ELEMENT_ROWS} biological "
        f"promoter-element rows; found {len(biological_rows)}"
    )

notch_rows = differential_elements.loc[
    differential_elements["source_mageck_id"]
    .astype(str)
    .str.upper()
    .eq("NOTCH1")
].copy()
if len(notch_rows) != 2:
    raise RuntimeError(
        f"Expected two NOTCH1 element rows; found {len(notch_rows)}"
    )
if set(notch_rows["id"]) != {"NOTCH1_-P1", "NOTCH1_-P2"}:
    raise RuntimeError(
        "Unexpected NOTCH1 output IDs: "
        + ", ".join(notch_rows["id"].astype(str))
    )
if not (notch_rows[native_stat_fields].nunique(dropna=False) == 1).all():
    raise RuntimeError(
        "NOTCH1 native MAGeCK statistics are not identical across P1/P2"
    )

source_num = pd.to_numeric(
    biological_rows["num"],
    errors="raise",
).astype(int)
actual_element_guides = pd.to_numeric(
    biological_rows["element_guide_count"],
    errors="raise",
).astype(int)
num_mismatch = biological_rows.loc[
    source_num.ne(actual_element_guides)
].copy()

NUM_MISMATCH_PATH = (
    OUTPUT_ROOT / "qc/mageck_num_vs_element_guide_count.tsv"
)
write_tsv(num_mismatch, NUM_MISMATCH_PATH)

DIFF_PLAIN = (
    OUTPUT_ROOT
    / "files/validation_screen_differential_element_quantifications.tsv"
)
DIFF_GZ = (
    OUTPUT_ROOT
    / "files/validation_screen_differential_element_quantifications.tsv.gz"
)
DIFF_MD5 = deterministic_tsv_gzip(
    differential_elements,
    DIFF_PLAIN,
    DIFF_GZ,
)

MAPPING_AUDIT_PATH = (
    OUTPUT_ROOT / "qc/differential_element_mapping_audit.tsv"
)
write_tsv(pd.DataFrame(mapping_audit_rows), MAPPING_AUDIT_PATH)

print(f"Native MAGeCK rows: {len(differential_raw):,}")
print(f"Final total rows: {len(differential_elements):,}")
print(f"Biological promoter-element rows: {len(biological_rows):,}")
print(
    "Control rows:",
    len(differential_elements) - len(biological_rows),
)
print("\nNOTCH1 representation:")
print(
    notch_rows[
        [
            "id",
            "element_chr",
            "element_start",
            "element_end",
            "element_guide_count",
            "num",
        ]
    ].to_string(index=False)
)
print(
    "\nRows where native MAGeCK num differs from element guide count:",
    len(num_mismatch),
)
print("Output:", DIFF_GZ)
print("MD5:", DIFF_MD5)
```

### Text output (stream)

```text
Native MAGeCK rows: 299
Final total rows: 300
Biological promoter-element rows: 298
Control rows: 2

NOTCH1 representation:
        id element_chr element_start element_end element_guide_count num
NOTCH1_-P1        chr9     136545758   136546062                   6  12
NOTCH1_-P2        chr9     136523466   136523919                   6  12

Rows where native MAGeCK num differs from element guide count: 2
Output: /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_differential_element_quantifications.tsv.gz
MD5: 2f0cb5fc62e1ba6d743cbabdef1cf8bf

```

## Markdown cell 26

## 7. Create column dictionaries, QC reports, and checksum manifest

## Code cell 27

```python
guide_dictionary_rows: list[dict] = []
for column in guide_counts.columns:
    if column == "guide_id":
        description = (
            "Unique guide RNA identifier used in the validation-screen "
            "library."
        )
        scale = "not applicable"
        allowed = ""
    elif column == target_column:
        description = (
            "Native target-analysis group from guideCounts.csv. It normally "
            "matches the guide target encoded in the corrected guide ID; the "
            "single documented NEDD8-MDP1 guide remains grouped as NEDD8, "
            "matching the original MAGeCK analysis."
        )
        scale = "not applicable"
        allowed = ""
    elif column in count_columns:
        sample_description = {
            "D01": "Day 0 baseline, biological replicate 1",
            "D02": "Day 0 baseline, biological replicate 2",
            "D3neg1": "Day 3 CDH5-negative fraction, biological replicate 1",
            "D3neg2": "Day 3 CDH5-negative fraction, biological replicate 2",
            "D3pos1": "Day 3 CDH5-positive fraction, biological replicate 1",
            "D3pos2": "Day 3 CDH5-positive fraction, biological replicate 2",
        }[column]
        description = (
            f"Raw guide RNA sequencing read count for {sample_description}."
        )
        scale = "read count"
        allowed = "non-negative integer"
    else:
        raise RuntimeError(
            f"Unexpected guide-quantification output column: {column}"
        )

    guide_dictionary_rows.append(
        {
            "table_type": "guide quantifications",
            "column_name": column,
            "data_type": infer_dictionary_type(guide_counts[column]),
            "required": "yes",
            "units_or_scale": scale,
            "allowed_values": allowed,
            "missing_value": "none expected",
            "description": description,
            "source_analysis_step": GUIDE_QUANTIFICATION_ASV,
        }
    )

guide_dictionary = pd.DataFrame(guide_dictionary_rows)
GUIDE_DICTIONARY_PATH = (
    OUTPUT_ROOT
    / "documentation/"
    "validation_screen_guide_quantifications_column_dictionary.tsv"
)
write_tsv(guide_dictionary, GUIDE_DICTIONARY_PATH)

differential_descriptions = {
    "id": (
        "Unique promoter-element identifier. NOTCH1 has separate P1 and "
        "P2 rows."
    ),
    "source_mageck_id": (
        "Original target identifier in the native gene-level MAGeCK "
        "summary."
    ),
    "gene_symbol": (
        "Gene symbol associated with the perturbed promoter element."
    ),
    "ensembl_gene_id": (
        "Versionless Ensembl gene identifier from Jesse's corrected guide "
        "metadata."
    ),
    "element_chr": "GRCh38 chromosome or contig of the promoter element.",
    "element_start": (
        "Start coordinate of the promoter element in corrected guide "
        "metadata."
    ),
    "element_end": (
        "End coordinate of the promoter element in corrected guide "
        "metadata."
    ),
    "element_guide_count": (
        "Number of validation-library guides assigned to this promoter "
        "element, including SGTN guides merged into a sole canonical "
        "promoter."
    ),
    "statistics_duplicated_across_elements": (
        "TRUE only when one native gene-level MAGeCK result was duplicated "
        "across multiple promoter elements; this occurs for NOTCH1."
    ),
    "num": (
        "Number of sgRNAs in the original MAGeCK target-level test. For "
        "each duplicated NOTCH1 row this remains the original gene-level "
        "value."
    ),
    "neg|score": "MAGeCK negative-selection score.",
    "neg|p-value": "MAGeCK negative-selection nominal p-value.",
    "neg|fdr": "MAGeCK negative-selection false-discovery rate.",
    "neg|rank": "MAGeCK negative-selection rank.",
    "neg|goodsgrna": (
        "Number of MAGeCK good sgRNAs supporting negative selection."
    ),
    "neg|lfc": (
        "MAGeCK log-fold-change value in the negative-selection direction."
    ),
    "pos|score": "MAGeCK positive-selection score.",
    "pos|p-value": "MAGeCK positive-selection nominal p-value.",
    "pos|fdr": "MAGeCK positive-selection false-discovery rate.",
    "pos|rank": "MAGeCK positive-selection rank.",
    "pos|goodsgrna": (
        "Number of MAGeCK good sgRNAs supporting positive selection."
    ),
    "pos|lfc": (
        "MAGeCK log-fold-change value in the positive-selection direction."
    ),
}

probability_columns = {
    "neg|score",
    "neg|p-value",
    "neg|fdr",
    "pos|score",
    "pos|p-value",
    "pos|fdr",
}

differential_dictionary_rows: list[dict] = []
for column in differential_elements.columns:
    if column in {"element_start", "element_end"}:
        scale = "GRCh38 genomic coordinate"
    elif column in probability_columns:
        scale = "unitless, 0 to 1"
    elif column.endswith("|lfc"):
        scale = "MAGeCK log fold change"
    else:
        scale = "not applicable"

    differential_dictionary_rows.append(
        {
            "table_type": "differential element quantifications",
            "column_name": column,
            "data_type": infer_dictionary_type(
                differential_elements[column]
            ),
            "required": "yes",
            "units_or_scale": scale,
            "allowed_values": (
                "TRUE or FALSE"
                if column == "statistics_duplicated_across_elements"
                else ""
            ),
            "missing_value": (
                "blank for aggregate control rows only"
                if column
                in {
                    "gene_symbol",
                    "ensembl_gene_id",
                    "element_chr",
                    "element_start",
                    "element_end",
                    "element_guide_count",
                }
                else "none expected"
            ),
            "description": differential_descriptions.get(
                column,
                "Native MAGeCK output field.",
            ),
            "source_analysis_step": DIFFERENTIAL_ELEMENT_ASV,
        }
    )

differential_dictionary = pd.DataFrame(
    differential_dictionary_rows
)
DIFF_DICTIONARY_PATH = (
    OUTPUT_ROOT
    / "documentation/"
    "validation_screen_differential_element_quantifications_"
    "column_dictionary.tsv"
)
write_tsv(differential_dictionary, DIFF_DICTIONARY_PATH)

FILE_MANIFEST_PATH = (
    OUTPUT_ROOT / "metadata/validation_screen_file_manifest.tsv"
)
file_manifest = pd.DataFrame(
    [
        {
            "role": "guide quantifications",
            "submitted_file_name": str(GUIDE_GZ),
            "md5sum": GUIDE_MD5,
            "size_bytes": GUIDE_GZ.stat().st_size,
            "rows": len(guide_counts),
            "columns": len(guide_counts.columns),
            "content_type": "guide quantifications",
            "analysis_step_version": GUIDE_QUANTIFICATION_ASV,
            "column_dictionary": str(GUIDE_DICTIONARY_PATH),
        },
        {
            "role": "differential element quantifications",
            "submitted_file_name": str(DIFF_GZ),
            "md5sum": DIFF_MD5,
            "size_bytes": DIFF_GZ.stat().st_size,
            "rows": len(differential_elements),
            "columns": len(differential_elements.columns),
            "content_type": "differential element quantifications",
            "analysis_step_version": DIFFERENTIAL_ELEMENT_ASV,
            "column_dictionary": str(DIFF_DICTIONARY_PATH),
        },
    ]
)
write_tsv(file_manifest, FILE_MANIFEST_PATH)

QC_JSON_PATH = (
    OUTPUT_ROOT / "qc/validation_screen_analysis_file_qc.json"
)
qc_summary = {
    "guide_quantifications": {
        "source": str(guide_counts_path),
        "rows": len(guide_counts),
        "sample_count_columns": len(count_columns),
        "unique_guide_ids": int(guide_counts["guide_id"].nunique()),
        "corrected_guide_metadata_id_set_match": count_ids == metadata_ids,
        "selected_source_guide_id_column": selected_guide_id_column,
        "duplicate_guide_id_columns_removed": duplicate_guide_id_columns,
        "count_columns": count_columns,
        "source_annotation_columns": annotation_columns,
        "target_column_retained": target_column is not None,
        "excluded_annotation_columns": excluded_annotation_columns,
        "removed_empty_columns": removed_empty_guide_columns,
        "removed_index_columns": removed_index_columns,
        "column_profile": str(GUIDE_COLUMN_PROFILE_PATH),
        "annotation_audit": str(ANNOTATION_AUDIT_PATH),
        "excluded_annotation_audit": str(EXCLUDED_ANNOTATIONS_PATH),

"guide_id_normalization_audit": str(ID_NORMALIZATION_AUDIT_PATH),
"curated_guide_id_alias_audit": str(
    CURATED_GUIDE_ALIAS_AUDIT_PATH
),
"curated_analysis_group_override_audit": str(
    ANALYSIS_GROUP_OVERRIDE_AUDIT_PATH
),
"curated_analysis_alias_element_merge_audit": str(
    ANALYSIS_ALIAS_MERGE_AUDIT_PATH
),
"curated_guide_id_aliases_applied": int(
    id_normalization_audit["curated_alias_applied"].sum()
),
        "reference_column_mapping": str(COUNT_COLUMN_MAPPING_PATH),
        "missing_value_repair_audit": str(COUNT_REPAIR_AUDIT_PATH),
        "reference_discrepancy_audit": str(COUNT_DISCREPANCY_PATH),
        "missing_count_values_filled": len(count_repair_rows),
        "output": str(GUIDE_GZ),
        "md5sum": GUIDE_MD5,
        "mageck_count_reference_crosscheck": reference_crosscheck,
    },
    "differential_element_quantifications": {
        "source": str(differential_path),
        "source_rows": len(differential_raw),
        "output_rows": len(differential_elements),
        "biological_element_rows": len(biological_rows),
        "control_rows": len(differential_elements) - len(biological_rows),
        "notch1_rows": len(notch_rows),
        "notch1_element_ids": sorted(notch_rows["id"].tolist()),
        "notch1_statistics_duplicated_identically": bool(
            (
                notch_rows[native_stat_fields].nunique(dropna=False)
                == 1
            ).all()
        ),
        "num_vs_element_guide_count_mismatch_rows": len(num_mismatch),
        "removed_empty_columns": removed_empty_differential_columns,
        "output": str(DIFF_GZ),
        "md5sum": DIFF_MD5,
    },
}
QC_JSON_PATH.write_text(
    json.dumps(qc_summary, indent=2) + "\n",
    encoding="utf-8",
)

MD5_MANIFEST_PATH = OUTPUT_ROOT / "metadata/md5sums.txt"
MD5_MANIFEST_PATH.write_text(
    f"{GUIDE_MD5}  {GUIDE_GZ}\n"
    f"{DIFF_MD5}  {DIFF_GZ}\n",
    encoding="utf-8",
)

print(file_manifest.to_string(index=False))
print("\nQC JSON:", QC_JSON_PATH)
print("Column dictionaries:")
print(" ", GUIDE_DICTIONARY_PATH)
print(" ", DIFF_DICTIONARY_PATH)
```

### Text output (stream)

```text
                                role                                                                                                                                    submitted_file_name                           md5sum  size_bytes  rows  columns                         content_type                                         analysis_step_version                                                                                                                                                             column_dictionary
               guide quantifications                /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_guide_quantifications.tsv.gz 0d4a1e607ef5f0da24f7fffa59e0d00c       40170  1926        8                guide quantifications /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/                /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/documentation/validation_screen_guide_quantifications_column_dictionary.tsv
differential element quantifications /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_differential_element_quantifications.tsv.gz 2f0cb5fc62e1ba6d743cbabdef1cf8bf       16449   300       22 differential element quantifications /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/ /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/documentation/validation_screen_differential_element_quantifications_column_dictionary.tsv

QC JSON: /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/qc/validation_screen_analysis_file_qc.json
Column dictionaries:
  /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/documentation/validation_screen_guide_quantifications_column_dictionary.tsv
  /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/documentation/validation_screen_differential_element_quantifications_column_dictionary.tsv

```

## Markdown cell 28

## 8. Generate IGVF `tabular_file` metadata

The metadata uses Ian's AnalysisStepVersion records and the established reference file `/reference-files/IGVFFI7217ZMJZ/`. `submitted_assembly` is deliberately omitted because it is admin-only. Empty optional fields are omitted rather than invented.

## Code cell 29

```python
if not re.fullmatch(r"IGVFDS[A-Z0-9]+", VALIDATION_ANALYSIS_SET):
    raise RuntimeError(
        "VALIDATION_ANALYSIS_SET is not configured. Set it in Section 1 "
        "to the production IGVFDS accession for the validation-screen "
        "AnalysisSet, then rerun Sections 1 and 8."
    )

metadata_rows = [
    {
        "aliases": GUIDE_QUANTIFICATION_ALIAS,
        "award": AWARD,
        "lab": LAB,
        "md5sum": GUIDE_MD5,
        "file_format": "tsv",
        "file_set": VALIDATION_ANALYSIS_SET,
        "content_type": "guide quantifications",
        "reference_files": REFERENCE_FILE,
        "filtered": "false",
        "controlled_access": "false",
        "analysis_step_version": GUIDE_QUANTIFICATION_ASV,
        "derived_manually": "false",
        "description": (
            "Guide-level sequencing read counts for the focused WTC11-EC "
            "validation CRISPRi screen across two biological replicates "
            "and D0, Day 3 CDH5-positive, and Day 3 CDH5-negative PCR "
            "samples."
        ),
        "submitted_file_name": str(GUIDE_GZ),
        "derived_from": ",".join(GUIDE_QUANTIFICATIONS_DERIVED_FROM),
        "file_format_specifications": GUIDE_QUANTIFICATIONS_FORMAT_SPEC,
    },
    {
        "aliases": DIFFERENTIAL_ELEMENT_ALIAS,
        "award": AWARD,
        "lab": LAB,
        "md5sum": DIFF_MD5,
        "file_format": "tsv",
        "file_set": VALIDATION_ANALYSIS_SET,
        "content_type": "differential element quantifications",
        "reference_files": REFERENCE_FILE,
        "filtered": "false",
        "controlled_access": "false",
        "analysis_step_version": DIFFERENTIAL_ELEMENT_ASV,
        "derived_manually": "false",
        "description": (
            "Promoter-element-level MAGeCK results comparing Day 3 "
            "CDH5-positive and CDH5-negative fractions in the WTC11-EC "
            "validation CRISPRi screen. The original NOTCH1 gene-level "
            "statistics are duplicated unchanged across its P1 and P2 "
            "promoter elements following DACC guidance."
        ),
        "submitted_file_name": str(DIFF_GZ),
        "derived_from": ",".join(DIFFERENTIAL_ELEMENTS_DERIVED_FROM),
        "file_format_specifications": DIFFERENTIAL_ELEMENTS_FORMAT_SPEC,
    },
]

metadata = pd.DataFrame(metadata_rows)
if "submitted_assembly" in metadata.columns:
    raise RuntimeError("submitted_assembly is admin-only and must not be submitted")
if "reference_files" not in metadata.columns:
    raise RuntimeError("reference_files is missing from metadata")
if not metadata["reference_files"].astype(str).eq(REFERENCE_FILE).all():
    raise RuntimeError("Unexpected reference_files value in metadata")
empty_optional_columns = [
    column
    for column in metadata.columns
    if metadata[column].astype(str).str.strip().eq("").all()
]
metadata = metadata.drop(columns=empty_optional_columns)

TABULAR_METADATA_PATH = (
    OUTPUT_ROOT / "metadata/tabular_file_validation_screen.tsv"
)
write_tsv(metadata, TABULAR_METADATA_PATH)

patch_rows: list[dict] = []
if EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID:
    patch_rows.append(
        {
            "record_id": EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID,
            "analysis_step_version": GUIDE_QUANTIFICATION_ASV,
        }
    )
if EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID:
    patch_rows.append(
        {
            "record_id": EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID,
            "analysis_step_version": DIFFERENTIAL_ELEMENT_ASV,
        }
    )

PATCH_METADATA_PATH = (
    OUTPUT_ROOT / "metadata/patch_analysis_step_versions.tsv"
)
if patch_rows:
    write_tsv(pd.DataFrame(patch_rows), PATCH_METADATA_PATH)
elif PATCH_METADATA_PATH.exists():
    PATCH_METADATA_PATH.unlink()

print("POST metadata:", TABULAR_METADATA_PATH)
print(metadata.to_string(index=False))
if patch_rows:
    print("\nPATCH metadata:", PATCH_METADATA_PATH)
    print(pd.DataFrame(patch_rows).to_string(index=False))
else:
    print("\nNo existing file IDs configured; no PATCH TSV was created.")
```

### Text output (stream)

```text
POST metadata: /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/metadata/tabular_file_validation_screen.tsv
                                                                                                          aliases             award                   lab                           md5sum file_format       file_set                         content_type                  reference_files filtered controlled_access                                         analysis_step_version derived_manually                                                                                                                                                                                                                                                                   description                                                                                                                                    submitted_file_name
                                     jesse-engreitz:tabular_file_WTC11-EC_validation_screen_guide_quantifications /awards/HG011972/ /labs/jesse-engreitz/ 0d4a1e607ef5f0da24f7fffa59e0d00c         tsv IGVFDS4708PCVX                guide quantifications /reference-files/IGVFFI7217ZMJZ/    false             false /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/            false                                                                                      Guide-level sequencing read counts for the focused WTC11-EC validation CRISPRi screen across two biological replicates and D0, Day 3 CDH5-positive, and Day 3 CDH5-negative PCR samples.                /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_guide_quantifications.tsv.gz
jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications /awards/HG011972/ /labs/jesse-engreitz/ 2f0cb5fc62e1ba6d743cbabdef1cf8bf         tsv IGVFDS4708PCVX differential element quantifications /reference-files/IGVFFI7217ZMJZ/    false             false /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/            false Promoter-element-level MAGeCK results comparing Day 3 CDH5-positive and CDH5-negative fractions in the WTC11-EC validation CRISPRi screen. The original NOTCH1 gene-level statistics are duplicated unchanged across its P1 and P2 promoter elements following DACC guidance. /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_differential_element_quantifications.tsv.gz

No existing file IDs configured; no PATCH TSV was created.

```

## Markdown cell 30

## 9. Local readiness gate and reviewer README

## Code cell 31

```python
for final_file, expected_md5 in [
    (GUIDE_GZ, GUIDE_MD5),
    (DIFF_GZ, DIFF_MD5),
]:
    if not final_file.is_file() or final_file.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty final file: {final_file}")
    with gzip.open(final_file, "rt", encoding="utf-8") as handle:
        first_line = handle.readline()
    if "\t" not in first_line:
        raise RuntimeError(f"Final file is not tab-delimited: {final_file}")
    observed_md5 = md5sum(final_file)
    if observed_md5 != expected_md5:
        raise RuntimeError(
            f"MD5 changed for {final_file}: expected {expected_md5}, "
            f"observed {observed_md5}"
        )

expected_content_types = {
    "guide quantifications",
    "differential element quantifications",
}
if set(metadata["content_type"]) != expected_content_types:
    raise RuntimeError(
        "Unexpected content types in generated Portal metadata"
    )

expected_asvs = {
    GUIDE_QUANTIFICATION_ASV,
    DIFFERENTIAL_ELEMENT_ASV,
}
if set(metadata["analysis_step_version"]) != expected_asvs:
    raise RuntimeError(
        "Generated metadata does not contain Ian's two ASV identifiers"
    )

if set(metadata["file_set"]) != {VALIDATION_ANALYSIS_SET}:
    raise RuntimeError(
        "The two files are not assigned to exactly one validation "
        "AnalysisSet"
    )


if "submitted_assembly" in metadata.columns:
    raise RuntimeError("submitted_assembly must not be present")
if "reference_files" not in metadata.columns:
    raise RuntimeError("reference_files is required")
if not metadata["reference_files"].astype(str).eq(REFERENCE_FILE).all():
    raise RuntimeError("Unexpected reference_files link")

README_PATH = (
    OUTPUT_ROOT / "README_validation_screen_analysis_submission.md"
)
readme_lines = [
    "# WTC11-EC validation CRISPRi screen analysis-file submission",
    "",
    f"AnalysisSet: `{VALIDATION_ANALYSIS_SET}`",
    "",
    "## Files",
    "",
    (
        f"1. `{GUIDE_GZ.name}` — guide quantifications; "
        f"AnalysisStepVersion `{GUIDE_QUANTIFICATION_ASV}`."
    ),
    (
        f"2. `{DIFF_GZ.name}` — differential element quantifications; "
        f"AnalysisStepVersion `{DIFFERENTIAL_ELEMENT_ASV}`."
    ),
    "",
    "## NOTCH1 representation",
    "",
    (
        "The native MAGeCK analysis aggregated all 12 NOTCH1 guides in "
        "one gene-level row. For IGVF per-element harmonization, the same "
        "original statistics are represented in two rows "
        "(`NOTCH1_-P1` and `NOTCH1_-P2`). The statistics were not "
        "independently re-estimated by promoter. `element_guide_count` "
        "reports the actual guide count per promoter, whereas native "
        "`num` remains the original gene-level value in both rows."
    ),
    "",
    "## Supplemental SGTN guides",
    "",
    (
        "Where a gene has one canonical promoter, supplemental SGTN "
        "guides are assigned to that promoter rather than emitted as "
        "separate alternate-contig elements."
    ),
    "",
    "## Supporting files",
    "",
    (
        "Column dictionaries, source checksums, promoter-element mapping, "
        "transformation audits, QC, and Portal metadata are stored under "
        "`documentation/`, `qc/`, and `metadata/`."
    ),
    "",
]
README_PATH.write_text(
    "\n".join(readme_lines),
    encoding="utf-8",
)

READY_MARKER = OUTPUT_ROOT / "READY_FOR_IGVF_DRY_RUN.txt"
READY_MARKER.write_text(
    "Local validation passed for two validation-screen analysis files.\n"
    f"AnalysisSet: {VALIDATION_ANALYSIS_SET}\n"
    f"Guide quantifications ASV: {GUIDE_QUANTIFICATION_ASV}\n"
    f"Differential element ASV: {DIFFERENTIAL_ELEMENT_ASV}\n",
    encoding="utf-8",
)

print(READY_MARKER.read_text())
print("Reviewer README:", README_PATH)
```

### Text output (stream)

```text
Local validation passed for two validation-screen analysis files.
AnalysisSet: IGVFDS4708PCVX
Guide quantifications ASV: /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/
Differential element ASV: /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/

Reviewer README: /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/README_validation_screen_analysis_submission.md

```

## Code cell 32

```python
# Freeze all validation-screen state before the genome-wide section reuses
# the same short variable names.
VALIDATION_PACKAGE = {
    "screen_key": "validation",
    "screen_label": "Focused validation screen",
    "analysis_set": VALIDATION_ANALYSIS_SET,
    "output_root": Path(OUTPUT_ROOT),
    "metadata": metadata.copy(deep=True),
    "metadata_path": Path(TABULAR_METADATA_PATH),
    "ready_marker": Path(READY_MARKER),
    "items": [
        {
            "role": "guide_quantifications",
            "alias": GUIDE_QUANTIFICATION_ALIAS,
            "existing_id": EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID,
            "expected_asv": GUIDE_QUANTIFICATION_ASV,
            "final_file": Path(GUIDE_GZ),
        },
        {
            "role": "differential_element_quantifications",
            "alias": DIFFERENTIAL_ELEMENT_ALIAS,
            "existing_id": EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID,
            "expected_asv": DIFFERENTIAL_ELEMENT_ASV,
            "final_file": Path(DIFF_GZ),
        },
    ],
}

if not VALIDATION_PACKAGE["ready_marker"].is_file():
    raise RuntimeError(
        "Validation-screen local readiness marker was not created"
    )

print("Captured validation-screen package:")
print(" ", VALIDATION_PACKAGE["metadata_path"])
print(" ", VALIDATION_PACKAGE["ready_marker"])
```

### Text output (stream)

```text
Captured validation-screen package:
  /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/metadata/tabular_file_validation_screen.tsv
  /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/READY_FOR_IGVF_DRY_RUN.txt

```

## Markdown cell 33

# Part B. Genome-wide screen

## Markdown cell 34

## 1. Configuration

Before metadata generation, set `GENOME_WIDE_ANALYSIS_SET` to the production `IGVFDS...` accession for the genome-wide screen.

The three source paths can be overridden explicitly when the documented locations differ on Sherlock.

## Code cell 35

```python
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from getpass import getpass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import requests

# -----------------------------------------------------------------------------
# Documented Sherlock source locations
# -----------------------------------------------------------------------------
GW_ROOT_USERS = Path(
    "/oak/stanford/groups/engreitz/Users/dulguun/"
    "230916_Genome-wide_Screen_Analysis"
)
GW_ROOT_USER = Path(
    "/oak/stanford/groups/engreitz/User/dulguun/"
    "230916_Genome-wide_Screen_Analysis"
)

GUIDE_COUNTS_SOURCE = (
    GW_ROOT_USERS
    / "full_result_files_from_Chad/Files_used_for_mageck/"
    "combined_mageck_input_trimmed.txt"
)
DIFFERENTIAL_SOURCE = (
    GW_ROOT_USERS
    / "full_result_files_from_Chad/Pos_vs_Neg/"
    "Pos_vs_Neg.gene_summary.txt"
)
LIBRARY_DESIGN_SOURCE = Path(
    "/oak/stanford/groups/engreitz/Projects/CRISPRDesign/"
    "221205_CRISPRiPool/221205_Genome_Wide.design.txt"
)

# Explicit overrides take priority when non-empty.
GUIDE_COUNTS_OVERRIDE = GENOME_WIDE_GUIDE_COUNTS_OVERRIDE_INPUT
DIFFERENTIAL_SOURCE_OVERRIDE = GENOME_WIDE_DIFFERENTIAL_SOURCE_OVERRIDE_INPUT
DIFFERENTIAL_SOURCE_EXPECTED_MD5 = (
    GENOME_WIDE_DIFFERENTIAL_SOURCE_EXPECTED_MD5_INPUT
)
LIBRARY_DESIGN_OVERRIDE = GENOME_WIDE_LIBRARY_DESIGN_OVERRIDE_INPUT

OUTPUT_ROOT = Path(GENOME_WIDE_OUTPUT_ROOT_INPUT)

# REQUIRED: replace with the production AnalysisSet accession.
GENOME_WIDE_ANALYSIS_SET = (
    GENOME_WIDE_ANALYSIS_SET_INPUT
    or os.environ.get(
        "IGVF_GENOME_WIDE_ANALYSIS_SET",
        "REPLACE_WITH_GENOME_WIDE_ANALYSIS_SET_IGVFDS",
    )
)

AWARD = "/awards/HG011972/"
LAB = "/labs/jesse-engreitz/"
REFERENCE_FILE = REFERENCE_FILE_INPUT

GUIDE_QUANTIFICATION_ASV = (
    "/analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/"
)
DIFFERENTIAL_ELEMENT_ASV = (
    "/analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/"
)

GUIDE_QUANTIFICATION_ALIAS = (
    "jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_"
    "guide_quantifications"
)
DIFFERENTIAL_ELEMENT_ALIAS = (
    "jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_"
    "D3_CDH5pos_vs_CDH5neg_differential_element_quantifications"
)

# Optional immediate-upstream IGVF File identifiers. Leave empty rather than
# guessing; multiple values are comma-separated by the metadata builder.
GUIDE_QUANTIFICATIONS_DERIVED_FROM = list(GENOME_WIDE_GUIDE_DERIVED_FROM_INPUT)
DIFFERENTIAL_ELEMENTS_DERIVED_FROM = list(GENOME_WIDE_DIFFERENTIAL_DERIVED_FROM_INPUT)

# Optional Portal Document identifiers describing the table formats.
GUIDE_QUANTIFICATIONS_FORMAT_SPEC = GENOME_WIDE_GUIDE_FORMAT_SPEC_INPUT
DIFFERENTIAL_ELEMENTS_FORMAT_SPEC = GENOME_WIDE_DIFFERENTIAL_FORMAT_SPEC_INPUT

# Fill these only when the corresponding file records already exist and merely
# need Ian's AnalysisStepVersion patched.
EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID = GENOME_WIDE_EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID_INPUT
EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID = GENOME_WIDE_EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID_INPUT

EXPECTED_GUIDES = 129_249
EXPECTED_COUNT_COLUMNS = ["D0", "Pos", "Neg"]
EXPECTED_DIFFERENTIAL_ROWS = 20_655
EXPECTED_EXCEL_REPAIRS = 198
STRICT_EXPECTED_DIMENSIONS = True

IGVF_PYTHON = Path(
    "/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python"
)
IU_REGISTER = Path(
    "/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/"
    "igvf_utils/MetaDataRegistration/iu_register.py"
)

for directory in [
    OUTPUT_ROOT,
    OUTPUT_ROOT / "files",
    OUTPUT_ROOT / "metadata",
    OUTPUT_ROOT / "documentation",
    OUTPUT_ROOT / "qc",
    OUTPUT_ROOT / "logs",
    OUTPUT_ROOT / "work",
]:
    directory.mkdir(parents=True, exist_ok=True)

print("Output root:", OUTPUT_ROOT)
print("Genome-wide AnalysisSet:", GENOME_WIDE_ANALYSIS_SET)
print("Guide-quantification ASV:", GUIDE_QUANTIFICATION_ASV)
print("Differential-element ASV:", DIFFERENTIAL_ELEMENT_ASV)
print("Reference file:", REFERENCE_FILE)
```

### Text output (stream)

```text
Output root: /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission
Genome-wide AnalysisSet: IGVFDS9629VXQK
Guide-quantification ASV: /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/
Differential-element ASV: /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/
Reference file: /reference-files/IGVFFI7217ZMJZ/

```

## Markdown cell 36

## 2. Utility functions

## Code cell 37

```python
def detect_delimiter(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            if line.strip():
                return "\t" if line.count("\t") >= line.count(",") else ","
    raise ValueError(f"No non-empty line found in {path}")


def read_table(path: Path, *, dtype: object = str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pd.read_csv(
        path,
        sep=detect_delimiter(path),
        dtype=dtype,
        keep_default_na=False,
        na_values=[],
        low_memory=False,
        compression="infer",
    )
    frame.columns = [str(column).strip() for column in frame.columns]
    return frame


def clean_empty_columns(frame: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    removed: list[str] = []
    kept: list[str] = []
    for column in frame.columns:
        values = frame[column].astype(str).str.strip()
        if bool(values.eq("").all()):
            removed.append(column)
        else:
            kept.append(column)
    return frame.loc[:, kept].copy(), removed


def drop_sequential_index_artifacts(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    output = frame.copy()
    removed: list[str] = []
    for column in list(output.columns):
        if not str(column).strip().lower().startswith("unnamed:"):
            continue
        numeric = pd.to_numeric(output[column], errors="coerce")
        if numeric.notna().all():
            observed = numeric.to_numpy(dtype=int)
            if np.array_equal(observed, np.arange(len(output))) or np.array_equal(
                observed, np.arange(1, len(output) + 1)
            ):
                output = output.drop(columns=[column])
                removed.append(column)
    return output, removed


def md5sum(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def deterministic_tsv_gzip(
    frame: pd.DataFrame,
    plain_path: Path,
    gzip_path: Path,
) -> str:
    plain_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(plain_path, sep="\t", index=False, lineterminator="\n")
    with plain_path.open("rb") as source, gzip_path.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_output,
            compresslevel=9,
            mtime=0,
        ) as compressed:
            shutil.copyfileobj(source, compressed)

    with gzip.open(gzip_path, "rt", encoding="utf-8") as handle:
        roundtrip = pd.read_csv(
            handle,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            na_values=[],
        )
    if list(roundtrip.columns) != list(frame.columns) or len(roundtrip) != len(frame):
        raise RuntimeError(f"Gzip round-trip validation failed for {gzip_path}")
    return md5sum(gzip_path)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def choose_one_existing(
    candidates: Sequence[Path],
    label: str,
) -> Path:
    """Select one source, accepting only byte-identical duplicate copies."""
    found: list[Path] = []
    seen: set[str] = set()

    for candidate in candidates:
        candidate = candidate.expanduser()
        if not candidate.is_file():
            continue
        key = str(candidate.resolve())
        if key not in seen:
            found.append(candidate)
            seen.add(key)

    if not found:
        attempted = "\n".join(f"  - {path}" for path in candidates)
        raise FileNotFoundError(
            f"No {label} file found. Checked:\n{attempted}"
        )

    if len(found) == 1:
        return found[0]

    checksum_rows = [
        {
            "path": str(path),
            "size_bytes": path.stat().st_size,
            "md5sum": md5sum(path),
        }
        for path in found
    ]
    unique_checksums = {row["md5sum"] for row in checksum_rows}

    if len(unique_checksums) == 1:
        selected = found[0]
        print(
            f"Found {len(found)} byte-identical copies for {label}; "
            f"selecting the first documented candidate:\n  {selected}"
        )
        for row in checksum_rows:
            print(
                f"  equivalent copy: {row['path']} "
                f"(MD5 {row['md5sum']})"
            )
        return selected

    options = "\n".join(
        f"  - {row['path']}  MD5={row['md5sum']}"
        for row in checksum_rows
    )
    raise RuntimeError(
        f"Multiple non-identical {label} files were found. "
        f"Set an explicit source override before continuing:\n{options}"
    )


def exact_find(
    roots: Sequence[Path],
    names: Sequence[str],
    *,
    maxdepth: int = 7,
) -> list[Path]:
    results: set[Path] = set()
    for root in roots:
        if not root.is_dir():
            continue
        for name in names:
            command = [
                "find", str(root), "-maxdepth", str(maxdepth), "-type", "f",
                "-name", name, "-print",
            ]
            completed = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            for line in completed.stdout.splitlines():
                candidate = Path(line.strip())
                if candidate.is_file():
                    results.add(candidate)
    return sorted(results)


def infer_design_name_column(
    design: pd.DataFrame,
    guide_ids: set[str],
) -> str:
    candidates: list[tuple[int, int, str]] = []
    for column in design.columns:
        values = design[column].astype(str).str.strip()
        overlap = int(values.isin(guide_ids).sum())
        unique = int(values.nunique(dropna=False))
        candidates.append((overlap, unique, str(column)))
    overlap, unique, selected = max(candidates)
    if overlap != len(guide_ids):
        details = sorted(candidates, reverse=True)[:10]
        raise RuntimeError(
            "No design column matched every guide ID. Best column statistics: "
            f"{details}"
        )
    if unique != len(design):
        raise RuntimeError(
            f"The selected design name column {selected!r} is not unique"
        )
    return selected


def target_from_design_name(name: str) -> str | None:
    match = re.match(r"^(.*?)(?=_[+-]_)", str(name).strip())
    return match.group(1) if match else None


EXCEL_DATE_LABEL = re.compile(
    r"^(?:[1-9]|[12][0-9]|3[01])-(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$",
    flags=re.IGNORECASE,
)


def looks_like_excel_date_label(value: str) -> bool:
    return bool(EXCEL_DATE_LABEL.fullmatch(str(value).strip()))


def classify_native_element_id(value: str) -> str:
    text = str(value).strip()
    if re.search(r"(?:_-|-)(?:P1P2|P1|P2)$", text):
        return "P1/P2/P1P2 promoter label"
    if re.search(r"(?:_-|-)(?:ENST\d+(?:\.\d+)?|NM_\d+(?:[._]\d+)?)$", text):
        return "gene plus transcript element label"
    if re.fullmatch(r"(?:ENST\d+(?:\.\d+)?|NM_\d+(?:[._]\d+)?)", text):
        return "transcript-only element label"
    if re.search(r"NM_\d+", text):
        return "other RefSeq-derived element label"
    return "other native element label"


def infer_dictionary_type(series: pd.Series) -> str:
    values = series.astype(str).str.strip()
    nonempty = values.loc[values.ne("")]
    if nonempty.empty:
        return "string"
    normalized = set(nonempty.str.upper().unique())
    if normalized <= {"TRUE", "FALSE"}:
        return "boolean"
    numeric = pd.to_numeric(nonempty, errors="coerce")
    if numeric.notna().all():
        array = numeric.to_numpy(dtype=float)
        if np.allclose(array, np.round(array)):
            return "integer"
        return "float"
    return "string"


def portal_url(identifier: str, *, frame: str | None = None) -> str:
    """Build an authenticated IGVF database-read URL."""
    value = str(identifier).strip().strip("/")
    query = "?format=json&datastore=database"
    if frame:
        query += f"&frame={frame}"
    return f"https://api.data.igvf.org/{value}/{query}"


def portal_get(
    identifier: str,
    *,
    timeout: int = 30,
    prefer_edit: bool = False,
) -> dict | None:
    """Read one Portal object, with a safe fallback across API frames.

    The ordinary database read is tried first unless `prefer_edit=True`.
    This avoids treating shared, released, or DACC-created objects as
    missing merely because the submitter cannot request an edit frame.
    """
    auth = None
    if (
        os.environ.get("IGVF_API_KEY")
        and os.environ.get("IGVF_SECRET_KEY")
    ):
        auth = (
            os.environ["IGVF_API_KEY"],
            os.environ["IGVF_SECRET_KEY"],
        )

    frames: list[str | None] = (
        ["edit", None, "object"]
        if prefer_edit
        else [None, "object", "edit"]
    )

    attempts: list[dict] = []
    saw_non_404_error = False

    for frame in frames:
        url = portal_url(identifier, frame=frame)
        response = requests.get(
            url,
            auth=auth,
            headers={"Accept": "application/json"},
            timeout=timeout,
        )
        attempts.append(
            {
                "frame": frame or "default",
                "status_code": response.status_code,
                "url": url,
                "response_excerpt": response.text[:300],
            }
        )

        if response.status_code == 200:
            return response.json()
        if response.status_code == 404:
            continue

        saw_non_404_error = True

        # Authentication errors should not be obscured by fallback attempts.
        if response.status_code == 401:
            response.raise_for_status()

        # A shared object may reject frame=edit but remain readable through
        # the ordinary object endpoint. Continue to the next frame.
        if response.status_code in {403, 422}:
            continue

        # Retry other frame-dependent failures through the remaining URLs.
        if response.status_code >= 400:
            continue

    if saw_non_404_error:
        raise RuntimeError(
            "IGVF Portal lookup failed for "
            f"{identifier!r}. Attempts: "
            + json.dumps(attempts, indent=2)
        )
    return None


print("Utility functions loaded.")
```

### Text output (stream)

```text
Utility functions loaded.

```

## Markdown cell 38

## 3. Locate and audit the source files

## Code cell 39

```python
search_roots = [GW_ROOT_USERS, GW_ROOT_USER]

if GUIDE_COUNTS_OVERRIDE:
    guide_counts_path = choose_one_existing(
        [Path(GUIDE_COUNTS_OVERRIDE)],
        "genome-wide guide-count override",
    )
else:
    count_candidates = [
        GUIDE_COUNTS_SOURCE,
        GW_ROOT_USER
        / "full_result_files_from_Chad/Files_used_for_mageck/"
        "combined_mageck_input_trimmed.txt",
        Path.cwd() / "combined_mageck_input_trimmed.txt",
        Path.cwd().parent / "combined_mageck_input_trimmed.txt",
    ]
    count_candidates.extend(
        exact_find(search_roots, ["combined_mageck_input_trimmed.txt"])
    )
    guide_counts_path = choose_one_existing(
        count_candidates,
        "genome-wide guide-count",
    )

if DIFFERENTIAL_SOURCE_OVERRIDE:
    differential_path = choose_one_existing(
        [Path(DIFFERENTIAL_SOURCE_OVERRIDE)],
        "genome-wide differential-element source override",
    )
else:
    differential_candidates = [
        DIFFERENTIAL_SOURCE,
        GW_ROOT_USER
        / "full_result_files_from_Chad/Pos_vs_Neg/"
        "Pos_vs_Neg.gene_summary.txt",
        Path.cwd() / "Pos_vs_Neg.gene_summary.txt",
        Path.cwd().parent / "Pos_vs_Neg.gene_summary.txt",
    ]
    differential_candidates.extend(
        exact_find(search_roots, ["Pos_vs_Neg.gene_summary.txt"])
    )
    differential_path = choose_one_existing(
        differential_candidates,
        "genome-wide D3 CDH5-positive versus CDH5-negative MAGeCK summary",
    )

if LIBRARY_DESIGN_OVERRIDE:
    library_design_path = choose_one_existing(
        [Path(LIBRARY_DESIGN_OVERRIDE)],
        "genome-wide library-design override",
    )
else:
    design_candidates = [
        LIBRARY_DESIGN_SOURCE,
        Path.cwd() / "221205_Genome_Wide.design.txt",
        Path.cwd().parent / "221205_Genome_Wide.design.txt",
    ]
    design_candidates.extend(
        exact_find(
            [Path("/oak/stanford/groups/engreitz/Projects/CRISPRDesign")],
            ["221205_Genome_Wide.design.txt"],
            maxdepth=5,
        )
    )
    library_design_path = choose_one_existing(
        design_candidates,
        "authoritative genome-wide library design",
    )

# Confirm the canonical differential source against the checksum supplied after
# comparing all three discovered copies.
observed_differential_md5 = md5sum(differential_path)
if (
    DIFFERENTIAL_SOURCE_EXPECTED_MD5
    and observed_differential_md5
    != DIFFERENTIAL_SOURCE_EXPECTED_MD5
):
    raise RuntimeError(
        "The selected genome-wide differential source MD5 changed. "
        f"Expected {DIFFERENTIAL_SOURCE_EXPECTED_MD5}; observed "
        f"{observed_differential_md5}: {differential_path}"
    )

equivalent_differential_copies = [
    GW_ROOT_USERS
    / "full_result_files_from_Chad/Pos_vs_Neg/"
    "Pos_vs_Neg.gene_summary.txt",
    GW_ROOT_USERS
    / "figures/Pos vs Neg/"
    "Pos_vs_Neg.gene_summary.txt",
    GW_ROOT_USERS
    / "mageck/Pos_vs_Neg/"
    "Pos_vs_Neg.gene_summary.txt",
]
equivalent_copy_rows = []
for candidate in equivalent_differential_copies:
    equivalent_copy_rows.append(
        {
            "path": str(candidate),
            "exists": candidate.is_file(),
            "size_bytes": (
                candidate.stat().st_size if candidate.is_file() else 0
            ),
            "md5sum": md5sum(candidate) if candidate.is_file() else "",
            "selected_canonical_source": (
                candidate.resolve() == differential_path.resolve()
                if candidate.is_file()
                else False
            ),
            "matches_confirmed_md5": (
                md5sum(candidate) == DIFFERENTIAL_SOURCE_EXPECTED_MD5
                if candidate.is_file()
                else False
            ),
        }
    )

EQUIVALENT_COPY_AUDIT_PATH = (
    OUTPUT_ROOT
    / "qc/genome_wide_differential_equivalent_copies.tsv"
)
write_tsv(
    pd.DataFrame(equivalent_copy_rows),
    EQUIVALENT_COPY_AUDIT_PATH,
)

existing_equivalent_rows = [
    row for row in equivalent_copy_rows if row["exists"]
]
if existing_equivalent_rows and not all(
    row["matches_confirmed_md5"]
    for row in existing_equivalent_rows
):
    raise RuntimeError(
        "At least one known genome-wide differential source copy no longer "
        f"matches the confirmed MD5. Review {EQUIVALENT_COPY_AUDIT_PATH}"
    )

source_inventory = pd.DataFrame(
    [
        {
            "role": "guide quantifications source",
            "path": str(guide_counts_path),
            "size_bytes": guide_counts_path.stat().st_size,
            "md5sum": md5sum(guide_counts_path),
        },
        {
            "role": "differential element source",
            "path": str(differential_path),
            "size_bytes": differential_path.stat().st_size,
            "md5sum": md5sum(differential_path),
        },
        {
            "role": "authoritative library design",
            "path": str(library_design_path),
            "size_bytes": library_design_path.stat().st_size,
            "md5sum": md5sum(library_design_path),
        },
    ]
)
SOURCE_INVENTORY_PATH = OUTPUT_ROOT / "qc/source_inventory.tsv"
write_tsv(source_inventory, SOURCE_INVENTORY_PATH)

print(source_inventory.to_string(index=False))
print("\nSelected canonical differential source:", differential_path)
print("Confirmed differential-source MD5:", observed_differential_md5)
print("Equivalent-copy audit:", EQUIVALENT_COPY_AUDIT_PATH)
print("Saved:", SOURCE_INVENTORY_PATH)
```

### Text output (stream)

```text
                        role                                                                                                                                                               path  size_bytes                           md5sum
guide quantifications source /oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Files_used_for_mageck/combined_mageck_input_trimmed.txt     6008551 c8b370516a0acb1f24751a4b65e98013
 differential element source                  /oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Pos_vs_Neg/Pos_vs_Neg.gene_summary.txt     1928870 257255394c229c685b178e7d91b0c9b2
authoritative library design                                                               /oak/stanford/groups/engreitz/Projects/CRISPRDesign/221205_CRISPRiPool/221205_Genome_Wide.design.txt    34316334 80e0e3624b2b6322daafa350f9cc480a

Selected canonical differential source: /oak/stanford/groups/engreitz/Users/dulguun/230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/Pos_vs_Neg/Pos_vs_Neg.gene_summary.txt
Confirmed differential-source MD5: 257255394c229c685b178e7d91b0c9b2
Equivalent-copy audit: /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/qc/genome_wide_differential_equivalent_copies.tsv
Saved: /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/qc/source_inventory.tsv

```

## Markdown cell 40

## 4. Rebuild the final guide-quantification file

The source and design must contain the same 129,249 unique guide IDs. Only target labels that both differ from the design-derived target and match an Excel date pattern are repaired. Any other mismatch blocks the notebook.

## Code cell 41

```python
guide_counts_raw = read_table(guide_counts_path)
guide_counts_raw, removed_empty_count_columns = clean_empty_columns(guide_counts_raw)
guide_counts_raw, removed_index_columns = drop_sequential_index_artifacts(
    guide_counts_raw
)

required_count_columns = ["sgRNA", "gene", *EXPECTED_COUNT_COLUMNS]
missing_count_columns = [
    column for column in required_count_columns if column not in guide_counts_raw.columns
]
if missing_count_columns:
    raise RuntimeError(
        "The guide-count source lacks expected columns: "
        + ", ".join(missing_count_columns)
    )

guide_counts = guide_counts_raw.loc[:, required_count_columns].copy()
for column in ["sgRNA", "gene"]:
    guide_counts[column] = guide_counts[column].astype(str).str.strip()

if guide_counts["sgRNA"].eq("").any():
    raise RuntimeError("Blank sgRNA identifiers were found")
if guide_counts["sgRNA"].duplicated().any():
    duplicate_path = OUTPUT_ROOT / "qc/duplicate_guide_ids.tsv"
    write_tsv(
        guide_counts.loc[guide_counts["sgRNA"].duplicated(keep=False)],
        duplicate_path,
    )
    raise RuntimeError(f"Duplicate guide IDs found; review {duplicate_path}")

for column in EXPECTED_COUNT_COLUMNS:
    numeric = pd.to_numeric(guide_counts[column], errors="coerce")
    if numeric.isna().any():
        bad_path = OUTPUT_ROOT / "qc" / f"non_numeric_counts_{column}.tsv"
        write_tsv(guide_counts.loc[numeric.isna(), ["sgRNA", column]], bad_path)
        raise RuntimeError(f"Non-numeric values in {column}; review {bad_path}")
    if (numeric < 0).any():
        raise RuntimeError(f"Negative counts found in {column}")
    values = numeric.to_numpy(dtype=float)
    if not np.allclose(values, np.round(values)):
        raise RuntimeError(f"Non-integer counts found in {column}")
    guide_counts[column] = numeric.round().astype("int64")

library_design = read_table(library_design_path)
library_design, removed_empty_design_columns = clean_empty_columns(library_design)

guide_ids = set(guide_counts["sgRNA"])
design_name_column = infer_design_name_column(library_design, guide_ids)
library_design[design_name_column] = (
    library_design[design_name_column].astype(str).str.strip()
)

if library_design[design_name_column].duplicated().any():
    raise RuntimeError(
        f"Duplicate names found in design column {design_name_column!r}"
    )

design_ids = set(library_design[design_name_column])
count_only = sorted(guide_ids - design_ids)
design_only = sorted(design_ids - guide_ids)
if count_only or design_only:
    mismatch_rows = [
        {"guide_id": value, "issue": "guide_count_only"} for value in count_only
    ] + [
        {"guide_id": value, "issue": "library_design_only"} for value in design_only
    ]
    mismatch_path = OUTPUT_ROOT / "qc/guide_design_identifier_mismatch.tsv"
    write_tsv(pd.DataFrame(mismatch_rows), mismatch_path)
    raise RuntimeError(
        f"Guide-count and design identifiers do not match; review {mismatch_path}"
    )

name_by_guide = dict(
    zip(
        library_design[design_name_column],
        library_design[design_name_column],
    )
)
guide_counts["design_name"] = guide_counts["sgRNA"].map(name_by_guide)
guide_counts["design_target"] = guide_counts["design_name"].map(
    target_from_design_name
)

biological_mask = guide_counts["design_target"].notna()
mismatch_mask = (
    biological_mask
    & guide_counts["gene"].ne(guide_counts["design_target"])
)
date_mismatch_mask = mismatch_mask & guide_counts["gene"].map(
    looks_like_excel_date_label
)
unexpected_mismatch_mask = mismatch_mask & ~date_mismatch_mask

unexpected_mismatches = guide_counts.loc[
    unexpected_mismatch_mask,
    ["sgRNA", "gene", "design_target", "design_name"],
].copy()
UNEXPECTED_MISMATCH_PATH = (
    OUTPUT_ROOT / "qc/unexpected_non_date_target_mismatches.tsv"
)
write_tsv(unexpected_mismatches, UNEXPECTED_MISMATCH_PATH)
if not unexpected_mismatches.empty:
    raise RuntimeError(
        "Found target-label mismatches that are not Excel-style dates. "
        f"Review {UNEXPECTED_MISMATCH_PATH}; no correction was made."
    )

repair_audit = guide_counts.loc[
    date_mismatch_mask,
    ["sgRNA", "gene", "design_target", "design_name"],
].copy()
repair_audit = repair_audit.rename(
    columns={"gene": "original_gene", "design_target": "corrected_gene"}
)
repair_audit["repair_reason"] = (
    "Excel date conversion repaired from exact-matched design name"
)
REPAIR_AUDIT_PATH = OUTPUT_ROOT / "qc/genome_wide_gene_label_repairs.tsv"
write_tsv(repair_audit, REPAIR_AUDIT_PATH)

repair_count = int(date_mismatch_mask.sum())
if STRICT_EXPECTED_DIMENSIONS and repair_count not in {0, EXPECTED_EXCEL_REPAIRS}:
    raise RuntimeError(
        f"Expected either 0 already-correct rows or {EXPECTED_EXCEL_REPAIRS} "
        f"Excel-date repairs; found {repair_count}. Review {REPAIR_AUDIT_PATH}."
    )

guide_counts.loc[date_mismatch_mask, "gene"] = guide_counts.loc[
    date_mismatch_mask, "design_target"
]
guide_counts = guide_counts.loc[:, required_count_columns]

if STRICT_EXPECTED_DIMENSIONS and len(guide_counts) != EXPECTED_GUIDES:
    raise RuntimeError(
        f"Expected {EXPECTED_GUIDES:,} guide rows; found {len(guide_counts):,}"
    )
if guide_counts.isna().any().any() or guide_counts.astype(str).eq("").any().any():
    raise RuntimeError("Missing values remain in the final guide-count table")

all_zero_count = int(
    guide_counts[EXPECTED_COUNT_COLUMNS].sum(axis=1).eq(0).sum()
)

GUIDE_PLAIN = OUTPUT_ROOT / "files/genome_wide_screen_guide_quantifications.tsv"
GUIDE_GZ = OUTPUT_ROOT / "files/genome_wide_screen_guide_quantifications.tsv.gz"
GUIDE_MD5 = deterministic_tsv_gzip(guide_counts, GUIDE_PLAIN, GUIDE_GZ)

repair_mapping = (
    repair_audit.groupby(["original_gene", "corrected_gene"], dropna=False)
    .size()
    .rename("number_of_guides")
    .reset_index()
)
write_tsv(
    repair_mapping,
    OUTPUT_ROOT / "qc/genome_wide_gene_label_repair_mapping.tsv",
)

print(f"Final guide rows: {len(guide_counts):,}")
print(f"Excel-date labels repaired: {repair_count}")
print(f"All-zero guide rows retained: {all_zero_count}")
print(f"Design name column: {design_name_column}")
print(f"Removed empty source columns: {removed_empty_count_columns}")
print(f"Removed source index columns: {removed_index_columns}")
print(f"Removed empty design columns: {removed_empty_design_columns}")
print("Output:", GUIDE_GZ)
print("MD5:", GUIDE_MD5)
```

### Text output (stream)

```text
Final guide rows: 129,249
Excel-date labels repaired: 198
All-zero guide rows retained: 7
Design name column: name
Removed empty source columns: []
Removed source index columns: []
Removed empty design columns: ['chr', 'start', 'end', 'score', 'strand', 'GuideSequence', 'GuideSequenceMinusG', 'locus', 'GuideSequenceWithPAM', 'guideSet', 'SSC', 'OffTargetScore', 'QualityScore']
Output: /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_guide_quantifications.tsv.gz
MD5: e62c0a100589176793b61af60ed41032

```

## Markdown cell 42

## 5. Rebuild the final differential-element-quantification file

The native MAGeCK table is already grouped by promoter/transcript element. The notebook validates and compresses it without changing its rows, IDs, or statistics. A separate audit classifies identifier formats and reports potential collisions after removing the inconsistent underscore before the promoter separator.

## Code cell 43

```python
differential_raw = read_table(differential_path)
differential_raw, removed_empty_differential_columns = clean_empty_columns(
    differential_raw
)

required_mageck_columns = [
    "id", "num",
    "neg|score", "neg|p-value", "neg|fdr", "neg|rank", "neg|goodsgrna", "neg|lfc",
    "pos|score", "pos|p-value", "pos|fdr", "pos|rank", "pos|goodsgrna", "pos|lfc",
]
missing_mageck_columns = [
    column for column in required_mageck_columns if column not in differential_raw.columns
]
if missing_mageck_columns:
    raise RuntimeError(
        "The MAGeCK summary lacks expected columns: "
        + ", ".join(missing_mageck_columns)
    )

differential_elements = differential_raw.loc[:, required_mageck_columns].copy()
differential_elements["id"] = differential_elements["id"].astype(str).str.strip()
if differential_elements["id"].eq("").any():
    raise RuntimeError("Blank MAGeCK element IDs were found")
if differential_elements["id"].duplicated().any():
    duplicate_path = OUTPUT_ROOT / "qc/duplicate_native_element_ids.tsv"
    write_tsv(
        differential_elements.loc[
            differential_elements["id"].duplicated(keep=False)
        ],
        duplicate_path,
    )
    raise RuntimeError(f"Duplicate native element IDs found; review {duplicate_path}")

if STRICT_EXPECTED_DIMENSIONS and len(differential_elements) != EXPECTED_DIFFERENTIAL_ROWS:
    raise RuntimeError(
        f"Expected {EXPECTED_DIFFERENTIAL_ROWS:,} MAGeCK rows; "
        f"found {len(differential_elements):,}"
    )

integer_fields = {"num", "neg|rank", "neg|goodsgrna", "pos|rank", "pos|goodsgrna"}
probability_fields = {
    "neg|score", "neg|p-value", "neg|fdr",
    "pos|score", "pos|p-value", "pos|fdr",
}
for column in required_mageck_columns[1:]:
    numeric = pd.to_numeric(differential_elements[column], errors="coerce")
    if numeric.isna().any():
        bad_path = OUTPUT_ROOT / "qc" / (
            "non_numeric_mageck_" + re.sub(r"[^A-Za-z0-9]+", "_", column) + ".tsv"
        )
        write_tsv(differential_elements.loc[numeric.isna()], bad_path)
        raise RuntimeError(f"Non-numeric values in {column}; review {bad_path}")
    if column in integer_fields:
        values = numeric.to_numpy(dtype=float)
        if not np.allclose(values, np.round(values)):
            raise RuntimeError(f"Non-integer values found in {column}")
    if column in probability_fields and ((numeric < 0) | (numeric > 1)).any():
        raise RuntimeError(f"Values outside [0,1] found in {column}")

id_audit = pd.DataFrame({
    "id": differential_elements["id"],
    "identifier_class": differential_elements["id"].map(classify_native_element_id),
})
id_audit["underscore_normalized_id_for_audit_only"] = (
    id_audit["id"].str.replace("_-", "-", regex=False)
)
normalized_counts = id_audit.groupby(
    "underscore_normalized_id_for_audit_only"
)["id"].transform("size")
id_audit["potential_collision_after_underscore_normalization"] = (
    normalized_counts > 1
)
ELEMENT_REFERENCE_PATH = (
    OUTPUT_ROOT / "documentation/genome_wide_native_element_id_audit.tsv"
)
write_tsv(id_audit, ELEMENT_REFERENCE_PATH)

collision_audit = id_audit.loc[
    id_audit["potential_collision_after_underscore_normalization"]
].sort_values(["underscore_normalized_id_for_audit_only", "id"])
MAPPING_AUDIT_PATH = (
    OUTPUT_ROOT / "qc/genome_wide_identifier_collision_audit.tsv"
)
write_tsv(collision_audit, MAPPING_AUDIT_PATH)

nedd8_rows = differential_elements.loc[
    differential_elements["id"].str.upper().str.startswith("NEDD8")
]
NEDD8_AUDIT_PATH = OUTPUT_ROOT / "qc/NEDD8_native_rows.tsv"
write_tsv(nedd8_rows, NEDD8_AUDIT_PATH)

DIFF_PLAIN = (
    OUTPUT_ROOT / "files/genome_wide_screen_differential_element_quantifications.tsv"
)
DIFF_GZ = (
    OUTPUT_ROOT / "files/genome_wide_screen_differential_element_quantifications.tsv.gz"
)
DIFF_MD5 = deterministic_tsv_gzip(
    differential_elements,
    DIFF_PLAIN,
    DIFF_GZ,
)

print(f"Final differential-element rows: {len(differential_elements):,}")
print("Identifier classes:")
print(id_audit["identifier_class"].value_counts().to_string())
print(
    "Potential collision groups after underscore-only normalization:",
    collision_audit["underscore_normalized_id_for_audit_only"].nunique(),
)
print("Native NEDD8 rows retained:", len(nedd8_rows))
print("Output:", DIFF_GZ)
print("MD5:", DIFF_MD5)
```

### Text output (stream)

```text
Final differential-element rows: 20,655
Identifier classes:
identifier_class
P1/P2/P1P2 promoter label             18892
gene plus transcript element label     1586
transcript-only element label           177
Potential collision groups after underscore-only normalization: 57
Native NEDD8 rows retained: 2
Output: /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_differential_element_quantifications.tsv.gz
MD5: f8d6aa0564246df55837da95c4888c60

```

## Markdown cell 44

## 6. Create column dictionaries, QC reports, and the file manifest

## Code cell 45

```python
guide_dictionary_descriptions = {
    "sgRNA": "Unique guide RNA identifier in the genome-wide CRISPRi library.",
    "gene": (
        "Target label associated with the guide. Excel-converted date labels were "
        "repaired by exact matching to the authoritative design file."
    ),
    "D0": "Raw guide sequencing read count in the Day 0 baseline sample.",
    "Pos": "Raw guide sequencing read count in the Day 3 CDH5-positive fraction.",
    "Neg": "Raw guide sequencing read count in the Day 3 CDH5-negative fraction.",
}

guide_dictionary = pd.DataFrame([
    {
        "table_type": "guide quantifications",
        "column_name": column,
        "data_type": infer_dictionary_type(guide_counts[column]),
        "required": "yes",
        "units_or_scale": "read count" if column in EXPECTED_COUNT_COLUMNS else "not applicable",
        "allowed_values": "non-negative integer" if column in EXPECTED_COUNT_COLUMNS else "",
        "missing_value": "none expected",
        "description": guide_dictionary_descriptions[column],
        "source_analysis_step": GUIDE_QUANTIFICATION_ASV,
    }
    for column in guide_counts.columns
])
GUIDE_DICTIONARY_PATH = (
    OUTPUT_ROOT / "documentation/genome_wide_guide_quantifications_column_dictionary.tsv"
)
write_tsv(guide_dictionary, GUIDE_DICTIONARY_PATH)

mageck_descriptions = {
    "id": (
        "Native promoter or transcript element identifier used by the genome-wide "
        "MAGeCK analysis; values are preserved without normalization."
    ),
    "num": "Number of sgRNAs assigned to the native MAGeCK element-level test.",
    "neg|score": "MAGeCK negative-selection score.",
    "neg|p-value": "MAGeCK negative-selection nominal p-value.",
    "neg|fdr": "MAGeCK negative-selection false-discovery rate.",
    "neg|rank": "MAGeCK negative-selection rank.",
    "neg|goodsgrna": "Number of MAGeCK good sgRNAs supporting negative selection.",
    "neg|lfc": "MAGeCK log-fold-change value in the negative-selection direction.",
    "pos|score": "MAGeCK positive-selection score.",
    "pos|p-value": "MAGeCK positive-selection nominal p-value.",
    "pos|fdr": "MAGeCK positive-selection false-discovery rate.",
    "pos|rank": "MAGeCK positive-selection rank.",
    "pos|goodsgrna": "Number of MAGeCK good sgRNAs supporting positive selection.",
    "pos|lfc": "MAGeCK log-fold-change value in the positive-selection direction.",
}

differential_dictionary = pd.DataFrame([
    {
        "table_type": "differential element quantifications",
        "column_name": column,
        "data_type": infer_dictionary_type(differential_elements[column]),
        "required": "yes",
        "units_or_scale": (
            "unitless, 0 to 1"
            if column in probability_fields
            else "MAGeCK log fold change"
            if column.endswith("|lfc")
            else "not applicable"
        ),
        "allowed_values": "",
        "missing_value": "none expected",
        "description": mageck_descriptions[column],
        "source_analysis_step": DIFFERENTIAL_ELEMENT_ASV,
    }
    for column in differential_elements.columns
])
DIFF_DICTIONARY_PATH = (
    OUTPUT_ROOT / "documentation/"
    "genome_wide_differential_element_quantifications_column_dictionary.tsv"
)
write_tsv(differential_dictionary, DIFF_DICTIONARY_PATH)

FILE_MANIFEST_PATH = OUTPUT_ROOT / "metadata/genome_wide_file_manifest.tsv"
file_manifest = pd.DataFrame([
    {
        "role": "guide quantifications",
        "submitted_file_name": str(GUIDE_GZ),
        "md5sum": GUIDE_MD5,
        "size_bytes": GUIDE_GZ.stat().st_size,
        "rows": len(guide_counts),
        "columns": len(guide_counts.columns),
        "content_type": "guide quantifications",
        "analysis_step_version": GUIDE_QUANTIFICATION_ASV,
        "column_dictionary": str(GUIDE_DICTIONARY_PATH),
    },
    {
        "role": "differential element quantifications",
        "submitted_file_name": str(DIFF_GZ),
        "md5sum": DIFF_MD5,
        "size_bytes": DIFF_GZ.stat().st_size,
        "rows": len(differential_elements),
        "columns": len(differential_elements.columns),
        "content_type": "differential element quantifications",
        "analysis_step_version": DIFFERENTIAL_ELEMENT_ASV,
        "column_dictionary": str(DIFF_DICTIONARY_PATH),
    },
])
write_tsv(file_manifest, FILE_MANIFEST_PATH)

QC_JSON_PATH = OUTPUT_ROOT / "qc/genome_wide_analysis_file_qc.json"
qc_summary = {
    "guide_quantifications": {
        "source": str(guide_counts_path),
        "library_design": str(library_design_path),
        "rows": len(guide_counts),
        "columns": guide_counts.columns.tolist(),
        "unique_guides": int(guide_counts["sgRNA"].nunique()),
        "excel_date_repairs": repair_count,
        "unexpected_non_date_mismatches": len(unexpected_mismatches),
        "all_zero_rows_retained": all_zero_count,
        "output": str(GUIDE_GZ),
        "md5sum": GUIDE_MD5,
    },
    "differential_element_quantifications": {
        "source": str(differential_path),
        "rows": len(differential_elements),
        "unique_native_ids": int(differential_elements["id"].nunique()),
        "identifier_classes": id_audit["identifier_class"].value_counts().to_dict(),
        "potential_collision_groups_after_underscore_normalization": int(
            collision_audit["underscore_normalized_id_for_audit_only"].nunique()
        ),
        "native_nedd8_rows": nedd8_rows["id"].tolist(),
        "native_ids_modified": False,
        "output": str(DIFF_GZ),
        "md5sum": DIFF_MD5,
    },
}
QC_JSON_PATH.write_text(json.dumps(qc_summary, indent=2) + "\n", encoding="utf-8")

MD5_MANIFEST_PATH = OUTPUT_ROOT / "metadata/md5sums.txt"
MD5_MANIFEST_PATH.write_text(
    f"{GUIDE_MD5}  {GUIDE_GZ}\n{DIFF_MD5}  {DIFF_GZ}\n",
    encoding="utf-8",
)

print(file_manifest.to_string(index=False))
print("\nQC JSON:", QC_JSON_PATH)
```

### Text output (stream)

```text
                                role                                                                                                                                      submitted_file_name                           md5sum  size_bytes   rows  columns                         content_type                                         analysis_step_version                                                                                                                                                        column_dictionary
               guide quantifications                /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_guide_quantifications.tsv.gz e62c0a100589176793b61af60ed41032     1646611 129249        5                guide quantifications /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/                /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/documentation/genome_wide_guide_quantifications_column_dictionary.tsv
differential element quantifications /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_differential_element_quantifications.tsv.gz f8d6aa0564246df55837da95c4888c60      616932  20655       14 differential element quantifications /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/ /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/documentation/genome_wide_differential_element_quantifications_column_dictionary.tsv

QC JSON: /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/qc/genome_wide_analysis_file_qc.json

```

## Markdown cell 46

## 7. Generate IGVF `tabular_file` metadata

The metadata uses Ian's AnalysisStepVersion records and the established reference file `/reference-files/IGVFFI7217ZMJZ/`. `submitted_assembly` is deliberately omitted because it is admin-only. Empty optional fields are omitted rather than invented.

## Code cell 47

```python
if not re.fullmatch(r"IGVFDS[A-Z0-9]+", GENOME_WIDE_ANALYSIS_SET):
    raise RuntimeError(
        "GENOME_WIDE_ANALYSIS_SET is not configured. Set it in Section 1 to the "
        "production IGVFDS accession, then rerun Sections 1 and 7."
    )

metadata_rows = [
    {
        "aliases": GUIDE_QUANTIFICATION_ALIAS,
        "award": AWARD,
        "lab": LAB,
        "md5sum": GUIDE_MD5,
        "file_format": "tsv",
        "file_set": GENOME_WIDE_ANALYSIS_SET,
        "content_type": "guide quantifications",
        "reference_files": REFERENCE_FILE,
        "filtered": "false",
        "controlled_access": "false",
        "analysis_step_version": GUIDE_QUANTIFICATION_ASV,
        "derived_manually": "false",
        "description": (
            "Guide-level sequencing read counts for the genome-wide WTC11-EC "
            "CRISPRi differentiation screen in Day 0 baseline and Day 3 "
            "CDH5-positive and CDH5-negative fractions. Excel-converted target "
            "labels were repaired using the exact-matched authoritative library design."
        ),
        "submitted_file_name": str(GUIDE_GZ),
        "derived_from": ",".join(GUIDE_QUANTIFICATIONS_DERIVED_FROM),
        "file_format_specifications": GUIDE_QUANTIFICATIONS_FORMAT_SPEC,
    },
    {
        "aliases": DIFFERENTIAL_ELEMENT_ALIAS,
        "award": AWARD,
        "lab": LAB,
        "md5sum": DIFF_MD5,
        "file_format": "tsv",
        "file_set": GENOME_WIDE_ANALYSIS_SET,
        "content_type": "differential element quantifications",
        "reference_files": REFERENCE_FILE,
        "filtered": "false",
        "controlled_access": "false",
        "analysis_step_version": DIFFERENTIAL_ELEMENT_ASV,
        "derived_manually": "false",
        "description": (
            "Native element-level MAGeCK results comparing Day 3 CDH5-positive "
            "and CDH5-negative fractions in the genome-wide WTC11-EC CRISPRi "
            "differentiation screen. Native promoter/transcript element identifiers "
            "and MAGeCK statistics are preserved."
        ),
        "submitted_file_name": str(DIFF_GZ),
        "derived_from": ",".join(DIFFERENTIAL_ELEMENTS_DERIVED_FROM),
        "file_format_specifications": DIFFERENTIAL_ELEMENTS_FORMAT_SPEC,
    },
]

metadata = pd.DataFrame(metadata_rows)
if "submitted_assembly" in metadata.columns:
    raise RuntimeError("submitted_assembly is admin-only and must not be submitted")
if "reference_files" not in metadata.columns:
    raise RuntimeError("reference_files is missing from metadata")
if not metadata["reference_files"].astype(str).eq(REFERENCE_FILE).all():
    raise RuntimeError("Unexpected reference_files value in metadata")
empty_optional_columns = [
    column for column in metadata.columns
    if metadata[column].astype(str).str.strip().eq("").all()
]
metadata = metadata.drop(columns=empty_optional_columns)

TABULAR_METADATA_PATH = OUTPUT_ROOT / "metadata/tabular_file_genome_wide_screen.tsv"
write_tsv(metadata, TABULAR_METADATA_PATH)

patch_rows: list[dict] = []
if EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID:
    patch_rows.append({
        "record_id": EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID,
        "analysis_step_version": GUIDE_QUANTIFICATION_ASV,
    })
if EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID:
    patch_rows.append({
        "record_id": EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID,
        "analysis_step_version": DIFFERENTIAL_ELEMENT_ASV,
    })

PATCH_METADATA_PATH = OUTPUT_ROOT / "metadata/patch_analysis_step_versions.tsv"
if patch_rows:
    write_tsv(pd.DataFrame(patch_rows), PATCH_METADATA_PATH)
elif PATCH_METADATA_PATH.exists():
    PATCH_METADATA_PATH.unlink()

print("POST metadata:", TABULAR_METADATA_PATH)
print(metadata.to_string(index=False))
if patch_rows:
    print("\nPATCH metadata:", PATCH_METADATA_PATH)
    print(pd.DataFrame(patch_rows).to_string(index=False))
else:
    print("\nNo existing file IDs configured; no PATCH TSV was created.")
```

### Text output (stream)

```text
POST metadata: /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/metadata/tabular_file_genome_wide_screen.tsv
                                                                                                           aliases             award                   lab                           md5sum file_format       file_set                         content_type                  reference_files filtered controlled_access                                         analysis_step_version derived_manually                                                                                                                                                                                                                                                             description                                                                                                                                      submitted_file_name
                                     jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications /awards/HG011972/ /labs/jesse-engreitz/ e62c0a100589176793b61af60ed41032         tsv IGVFDS9629VXQK                guide quantifications /reference-files/IGVFFI7217ZMJZ/    false             false /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/            false Guide-level sequencing read counts for the genome-wide WTC11-EC CRISPRi differentiation screen in Day 0 baseline and Day 3 CDH5-positive and CDH5-negative fractions. Excel-converted target labels were repaired using the exact-matched authoritative library design.                /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_guide_quantifications.tsv.gz
jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications /awards/HG011972/ /labs/jesse-engreitz/ f8d6aa0564246df55837da95c4888c60         tsv IGVFDS9629VXQK differential element quantifications /reference-files/IGVFFI7217ZMJZ/    false             false /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/            false                           Native element-level MAGeCK results comparing Day 3 CDH5-positive and CDH5-negative fractions in the genome-wide WTC11-EC CRISPRi differentiation screen. Native promoter/transcript element identifiers and MAGeCK statistics are preserved. /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_differential_element_quantifications.tsv.gz

No existing file IDs configured; no PATCH TSV was created.

```

## Markdown cell 48

## 8. Local readiness gate and reviewer README

## Code cell 49

```python
for final_file, expected_md5 in [(GUIDE_GZ, GUIDE_MD5), (DIFF_GZ, DIFF_MD5)]:
    if not final_file.is_file() or final_file.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty final file: {final_file}")
    with gzip.open(final_file, "rt", encoding="utf-8") as handle:
        header = handle.readline()
    if "\t" not in header:
        raise RuntimeError(f"Final file is not tab-delimited: {final_file}")
    if md5sum(final_file) != expected_md5:
        raise RuntimeError(f"MD5 changed unexpectedly for {final_file}")

if set(metadata["content_type"]) != {
    "guide quantifications", "differential element quantifications"
}:
    raise RuntimeError("Unexpected content types in generated metadata")
if set(metadata["analysis_step_version"]) != {
    GUIDE_QUANTIFICATION_ASV, DIFFERENTIAL_ELEMENT_ASV
}:
    raise RuntimeError("Generated metadata does not use Ian's two ASV records")
if set(metadata["file_set"]) != {GENOME_WIDE_ANALYSIS_SET}:
    raise RuntimeError("The two files do not map to exactly one genome-wide AnalysisSet")


if "submitted_assembly" in metadata.columns:
    raise RuntimeError("submitted_assembly must not be present")
if "reference_files" not in metadata.columns:
    raise RuntimeError("reference_files is required")
if not metadata["reference_files"].astype(str).eq(REFERENCE_FILE).all():
    raise RuntimeError("Unexpected reference_files link")

README_PATH = OUTPUT_ROOT / "README_genome_wide_screen_analysis_submission.md"
README_PATH.write_text(
    "\n".join([
        "# WTC11-EC genome-wide CRISPRi screen analysis-file submission",
        "",
        f"AnalysisSet: `{GENOME_WIDE_ANALYSIS_SET}`",
        "",
        "## Submitted analysis files",
        "",
        f"1. `{GUIDE_GZ.name}` — guide quantifications.",
        f"2. `{DIFF_GZ.name}` — differential element quantifications.",
        "",
        "## Guide target-label repair",
        "",
        (
            "The original count table contained Excel-converted labels such as "
            "`1-Mar` and `1-Sep`. Only exact guide-ID matches whose source target "
            "was an Excel-style date and whose design-derived target differed were "
            "repaired. Any other mismatch blocks the notebook."
        ),
        "",
        "## Native element IDs",
        "",
        (
            "The MAGeCK result is already element-level. Its native IDs and statistics "
            "are preserved. An audit reports potential collisions after underscore-only "
            "normalization, including the documented NEDD8 naming inconsistency; the "
            "submission file does not silently merge those rows."
        ),
        "",
    ]) + "\n",
    encoding="utf-8",
)

READY_MARKER = OUTPUT_ROOT / "READY_FOR_IGVF_DRY_RUN.txt"
READY_MARKER.write_text(
    "Local validation passed for two genome-wide screen analysis files.\n"
    f"AnalysisSet: {GENOME_WIDE_ANALYSIS_SET}\n"
    f"Guide quantifications ASV: {GUIDE_QUANTIFICATION_ASV}\n"
    f"Differential element ASV: {DIFFERENTIAL_ELEMENT_ASV}\n",
    encoding="utf-8",
)

print(READY_MARKER.read_text())
print("Reviewer README:", README_PATH)
```

### Text output (stream)

```text
Local validation passed for two genome-wide screen analysis files.
AnalysisSet: IGVFDS9629VXQK
Guide quantifications ASV: /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/
Differential element ASV: /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/

Reviewer README: /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/README_genome_wide_screen_analysis_submission.md

```

## Code cell 50

```python
GENOME_WIDE_PACKAGE = {
    "screen_key": "genome_wide",
    "screen_label": "Genome-wide screen",
    "analysis_set": GENOME_WIDE_ANALYSIS_SET,
    "output_root": Path(OUTPUT_ROOT),
    "metadata": metadata.copy(deep=True),
    "metadata_path": Path(TABULAR_METADATA_PATH),
    "ready_marker": Path(READY_MARKER),
    "items": [
        {
            "role": "guide_quantifications",
            "alias": GUIDE_QUANTIFICATION_ALIAS,
            "existing_id": EXISTING_GUIDE_QUANTIFICATIONS_FILE_ID,
            "expected_asv": GUIDE_QUANTIFICATION_ASV,
            "final_file": Path(GUIDE_GZ),
        },
        {
            "role": "differential_element_quantifications",
            "alias": DIFFERENTIAL_ELEMENT_ALIAS,
            "existing_id": EXISTING_DIFFERENTIAL_ELEMENTS_FILE_ID,
            "expected_asv": DIFFERENTIAL_ELEMENT_ASV,
            "final_file": Path(DIFF_GZ),
        },
    ],
}

if not GENOME_WIDE_PACKAGE["ready_marker"].is_file():
    raise RuntimeError(
        "Genome-wide local readiness marker was not created"
    )

print("Captured genome-wide package:")
print(" ", GENOME_WIDE_PACKAGE["metadata_path"])
print(" ", GENOME_WIDE_PACKAGE["ready_marker"])
```

### Text output (stream)

```text
Captured genome-wide package:
  /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/metadata/tabular_file_genome_wide_screen.tsv
  /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/READY_FOR_IGVF_DRY_RUN.txt

```

## Markdown cell 51

# Part C. Combined four-file submission

The preparation sections above retain separate screen-specific output folders. The cells below create a single combined inventory and perform one coordinated Portal workflow for all four records.

## Code cell 52

```python
SCREEN_PACKAGES = [VALIDATION_PACKAGE, GENOME_WIDE_PACKAGE]

combined_metadata_frames: list[pd.DataFrame] = []
combined_manifest_rows: list[dict] = []
combined_items: list[dict] = []

for package in SCREEN_PACKAGES:
    package_metadata = package["metadata"].copy()
    package_metadata["__screen_key"] = package["screen_key"]
    package_metadata["__screen_label"] = package["screen_label"]
    combined_metadata_frames.append(package_metadata)

    metadata_by_alias = {
        str(row["aliases"]): row
        for row in package["metadata"].to_dict("records")
    }

    for item in package["items"]:
        alias = str(item["alias"])
        expected = metadata_by_alias[alias]
        final_file = Path(item["final_file"])

        if not final_file.is_file() or final_file.stat().st_size == 0:
            raise RuntimeError(
                f"Missing final file for {package['screen_label']}: "
                f"{final_file}"
            )
        observed_md5 = md5sum(final_file)
        if observed_md5 != str(expected["md5sum"]):
            raise RuntimeError(
                f"MD5 mismatch for {final_file}: metadata="
                f"{expected['md5sum']} observed={observed_md5}"
            )

        combined_item = {
            **item,
            "screen_key": package["screen_key"],
            "screen_label": package["screen_label"],
            "analysis_set": package["analysis_set"],
            "output_root": package["output_root"],
            "expected_metadata": expected,
        }
        combined_items.append(combined_item)

        combined_manifest_rows.append(
            {
                "screen_key": package["screen_key"],
                "screen_label": package["screen_label"],
                "role": item["role"],
                "alias": alias,
                "analysis_set": package["analysis_set"],
                "content_type": expected["content_type"],
                "analysis_step_version": expected[
                    "analysis_step_version"
                ],
                "reference_files": expected["reference_files"],
                "submitted_file_name": str(final_file),
                "md5sum": observed_md5,
                "size_bytes": final_file.stat().st_size,
                "screen_metadata_tsv": str(package["metadata_path"]),
                "screen_ready_marker": str(package["ready_marker"]),
            }
        )

combined_metadata = (
    pd.concat(combined_metadata_frames, ignore_index=True)
    .fillna("")
)
combined_manifest = pd.DataFrame(combined_manifest_rows)

if len(combined_items) != 4:
    raise RuntimeError(
        f"Expected four combined submission items; found {len(combined_items)}"
    )
if combined_manifest["alias"].duplicated().any():
    raise RuntimeError("Combined submission aliases are not unique")
if set(combined_manifest["screen_key"]) != {
    "validation",
    "genome_wide",
}:
    raise RuntimeError("Both screen packages were not captured")
if not (
    combined_manifest.groupby("screen_key").size() == 2
).all():
    raise RuntimeError("Expected exactly two files per screen")

COMBINED_MANIFEST_PATH = (
    COMBINED_OUTPUT_ROOT / "metadata/combined_four_file_manifest.tsv"
)
write_tsv(combined_manifest, COMBINED_MANIFEST_PATH)

# This four-row TSV is for review. One-row TSVs are used for actual
# submission so a partial run remains unambiguous.
COMBINED_METADATA_REVIEW_PATH = (
    COMBINED_OUTPUT_ROOT
    / "metadata/tabular_file_both_screens_review.tsv"
)
review_metadata = combined_metadata.copy()
review_metadata = review_metadata[
    ["__screen_key", "__screen_label"]
    + [
        column
        for column in review_metadata.columns
        if not column.startswith("__")
    ]
]
write_tsv(review_metadata, COMBINED_METADATA_REVIEW_PATH)

print(combined_manifest.to_string(index=False))
print("\nCombined manifest:", COMBINED_MANIFEST_PATH)
print("Combined metadata review:", COMBINED_METADATA_REVIEW_PATH)
```

### Text output (stream)

```text
 screen_key              screen_label                                 role                                                                                                              alias   analysis_set                         content_type                                         analysis_step_version                  reference_files                                                                                                                                      submitted_file_name                           md5sum  size_bytes                                                                                                              screen_metadata_tsv                                                                                            screen_ready_marker
 validation Focused validation screen                guide_quantifications                                       jesse-engreitz:tabular_file_WTC11-EC_validation_screen_guide_quantifications IGVFDS4708PCVX                guide quantifications /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/ /reference-files/IGVFFI7217ZMJZ/                  /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_guide_quantifications.tsv.gz 0d4a1e607ef5f0da24f7fffa59e0d00c       40170   /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/metadata/tabular_file_validation_screen.tsv  /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/READY_FOR_IGVF_DRY_RUN.txt
 validation Focused validation screen differential_element_quantifications  jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications IGVFDS4708PCVX differential element quantifications /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/ /reference-files/IGVFFI7217ZMJZ/   /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_differential_element_quantifications.tsv.gz 2f0cb5fc62e1ba6d743cbabdef1cf8bf       16449   /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/metadata/tabular_file_validation_screen.tsv  /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/READY_FOR_IGVF_DRY_RUN.txt
genome_wide        Genome-wide screen                guide_quantifications                                      jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications IGVFDS9629VXQK                guide quantifications /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/ /reference-files/IGVFFI7217ZMJZ/                /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_guide_quantifications.tsv.gz e62c0a100589176793b61af60ed41032     1646611 /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/metadata/tabular_file_genome_wide_screen.tsv /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/READY_FOR_IGVF_DRY_RUN.txt
genome_wide        Genome-wide screen differential_element_quantifications jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications IGVFDS9629VXQK differential element quantifications /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/ /reference-files/IGVFFI7217ZMJZ/ /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_differential_element_quantifications.tsv.gz f8d6aa0564246df55837da95c4888c60      616932 /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/metadata/tabular_file_genome_wide_screen.tsv /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/READY_FOR_IGVF_DRY_RUN.txt

Combined manifest: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/combined_four_file_manifest.tsv
Combined metadata review: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/tabular_file_both_screens_review.tsv

```

## Markdown cell 53

## C1. Enter IGVF credentials

Run the first code cell below to define the secure credential helpers. Then run the **second code cell**, which opens hidden prompts for the API key and secret key.

Paste the raw values without quotation marks. The values are kept only in the current notebook kernel environment and are inherited by the external `igvf_utils` subprocess. They are not written into the notebook or output files.

After a kernel restart, rerun both C1 cells before the Portal preflight in Part C2.


## Code cell 54

```python
import os
from getpass import getpass


def enter_igvf_credentials() -> None:
    """Prompt securely for IGVF production credentials."""
    api_key = <REDACTED>"IGVF API key: ").strip()
    secret_key = getpass("IGVF secret key: ").strip()

    if not api_key:
        <REDACTED> ValueError("IGVF API key was empty")
    if not secret_key:
        raise ValueError("IGVF secret key was empty")

    # Catch a common copy/paste mistake without exposing the credentials.
    for name, value in {
        "IGVF_API_KEY": api_key,
        "IGVF_SECRET_KEY": secret_key,
    }.items():
        if (
            len(value) >= 2
            and value[0] in {"'", '"'}
            and value[-1] == value[0]
        ):
            raise ValueError(
                f"{name} has surrounding quotation marks. "
                "Paste the raw value without quotes."
            )

    # Subprocesses launched by this notebook inherit these variables.
    os.environ["IGVF_API_KEY"] = api_key
    os.environ["IGVF_SECRET_KEY"] = secret_key

    del api_key, secret_key
    print("IGVF credentials are set for the current notebook kernel.")


def show_credential_status() -> None:
    """Show only safe credential diagnostics, never the values."""
    for variable in ("IGVF_API_KEY", "IGVF_SECRET_KEY"):
        value = os.environ.get(variable, "")
        print(
            f"{variable}: {'SET' if value else 'MISSING'}, "
            f"length={len(value)}, "
            f"surrounding_whitespace={value != value.strip()}"
        )

```

## Code cell 55

```python
# Run this cell once after every kernel restart, before Part C2.
enter_igvf_credentials()
show_credential_status()

```

### Text output (stream)

```text
IGVF credentials are set for the current notebook kernel.
IGVF_API_KEY: <REDACTED> length=8, surrounding_whitespace=False
IGVF_SECRET_KEY: SET, length=16, surrounding_whitespace=False

```

## Markdown cell 56

### C1A. Optional direct dependency diagnostic

This read-only cell tests the five dependencies using the safe read
endpoint before the full C2 submission plan is built.

## Code cell 57

```python
dependency_diagnostic_specs = [
    (
        "shared reference file",
        REFERENCE_FILE_INPUT,
    ),
    (
        "validation AnalysisSet",
        VALIDATION_PACKAGE["analysis_set"],
    ),
    (
        "genome-wide AnalysisSet",
        GENOME_WIDE_PACKAGE["analysis_set"],
    ),
    (
        "guide quantifications AnalysisStepVersion",
        GUIDE_QUANTIFICATION_ASV,
    ),
    (
        "differential element AnalysisStepVersion",
        DIFFERENTIAL_ELEMENT_ASV,
    ),
]

dependency_diagnostic_rows = []
for role, identifier in dependency_diagnostic_specs:
    try:
        record = portal_get(str(identifier))
        dependency_diagnostic_rows.append(
            {
                "role": role,
                "identifier": identifier,
                "found": bool(record),
                "accession": str(
                    (record or {}).get("accession", "")
                ),
                "portal_id": str(
                    (record or {}).get("@id", "")
                ),
                "status": str(
                    (record or {}).get("status", "")
                ),
                "object_type": ",".join(
                    (record or {}).get("@type", [])
                ),
                "error": "",
            }
        )
    except Exception as error:
        dependency_diagnostic_rows.append(
            {
                "role": role,
                "identifier": identifier,
                "found": False,
                "accession": "",
                "portal_id": "",
                "status": "",
                "object_type": "",
                "error": (
                    f"{type(error).__name__}: {error}"
                ),
            }
        )

dependency_diagnostic = pd.DataFrame(
    dependency_diagnostic_rows
)
DEPENDENCY_DIAGNOSTIC_PATH = (
    COMBINED_OUTPUT_ROOT
    / "qc/portal_dependency_diagnostic.tsv"
)
write_tsv(
    dependency_diagnostic,
    DEPENDENCY_DIAGNOSTIC_PATH,
)
print(dependency_diagnostic.to_string(index=False))
print("\\nSaved:", DEPENDENCY_DIAGNOSTIC_PATH)
```

### Text output (stream)

```text
                                     role                                                    identifier  found      accession                                                     portal_id      status              object_type error
                    shared reference file                              /reference-files/IGVFFI7217ZMJZ/   True IGVFFI7217ZMJZ                              /reference-files/IGVFFI7217ZMJZ/    archived  ReferenceFile,File,Item      
                   validation AnalysisSet                                                IGVFDS4708PCVX   True IGVFDS4708PCVX                                /analysis-sets/IGVFDS4708PCVX/ in progress AnalysisSet,FileSet,Item      
                  genome-wide AnalysisSet                                                IGVFDS9629VXQK   True IGVFDS9629VXQK                                /analysis-sets/IGVFDS9629VXQK/ in progress AnalysisSet,FileSet,Item      
guide quantifications AnalysisStepVersion /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/   True                /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/ in progress AnalysisStepVersion,Item      
 differential element AnalysisStepVersion /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/   True                /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/ in progress AnalysisStepVersion,Item      
\nSaved: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/qc/portal_dependency_diagnostic.tsv

```

## Markdown cell 58

## C2. Read-only Portal preflight and operation plan

This cell verifies both AnalysisSets and both AnalysisStepVersions, checks all four aliases, and creates an idempotent operation plan.

An existing record is considered compatible only when its MD5, AnalysisSet, content type, and file format match the local file. A conflict blocks all production work.

## Code cell 59

```python
def split_link_items(value: object) -> list[object]:
    if value is None:
        return []
    if isinstance(value, str):
        return [part.strip() for part in value.split(",") if part.strip()]
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]


def direct_link_tokens(value: object) -> set[str]:
    """Collect accession, UUID, @id terminal, and aliases from linked values."""
    tokens: set[str] = set()
    for item in split_link_items(value):
        if isinstance(item, dict):
            for key in ("accession", "uuid", "@id"):
                raw = str(item.get(key, "") or "").strip().strip("/")
                if raw:
                    tokens.add(raw.rsplit("/", 1)[-1])
            for alias in split_link_items(item.get("aliases", [])):
                alias_text = str(alias).strip()
                if alias_text:
                    tokens.add(alias_text)
        else:
            raw = str(item or "").strip().strip("/")
            if raw:
                tokens.add(raw.rsplit("/", 1)[-1])
                if ":" in raw and "/" not in raw:
                    tokens.add(raw)
    return tokens


def terminal_identifier(value: object) -> str:
    tokens = sorted(direct_link_tokens(value))
    return tokens[0] if tokens else ""


def normalized_link_set(value: object) -> set[str]:
    return direct_link_tokens(value)


def compact_link_display(value: object) -> str:
    """Render linked Portal values without serializing large embedded objects."""
    if value is None:
        return ""
    if isinstance(value, dict):
        for key in ("accession", "uuid", "@id"):
            token = terminal_identifier(value.get(key, ""))
            if token:
                return token
        aliases = value.get("aliases", [])
        alias_tokens = sorted(direct_link_tokens(aliases))
        return ",".join(alias_tokens)
    if isinstance(value, (list, tuple, set)):
        displays = [compact_link_display(item) for item in value]
        return ",".join(sorted({item for item in displays if item}))
    return terminal_identifier(value)


def links_equal(left: object, right: object) -> bool:
    left_tokens = direct_link_tokens(left)
    right_tokens = direct_link_tokens(right)
    return bool(left_tokens and right_tokens and left_tokens.intersection(right_tokens))


def link_sets_equal(left: object, right: object) -> bool:
    left_items = split_link_items(left)
    right_items = split_link_items(right)
    if len(left_items) != len(right_items):
        return False
    unmatched = list(right_items)
    for left_item in left_items:
        match = next((idx for idx, right_item in enumerate(unmatched) if links_equal(left_item, right_item)), None)
        if match is None:
            return False
        unmatched.pop(match)
    return not unmatched
def portal_record_url(record: dict | None) -> str:
    record = record or {}
    accession = str(record.get("accession") or "")
    if accession:
        return f"https://data.igvf.org/tabular-files/{accession}/"
    portal_id = str(record.get("@id") or "")
    return f"https://data.igvf.org{portal_id}" if portal_id else ""


def get_with_retry(
    identifier: str,
    *,
    attempts: int = 5,
) -> dict | None:
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


if not (
    os.environ.get("IGVF_API_KEY")
    and os.environ.get("IGVF_SECRET_KEY")
):
    raise RuntimeError(
        "Run enter_igvf_credentials() before the Portal preflight"
    )

for package in SCREEN_PACKAGES:
    if not re.fullmatch(
        r"IGVFDS[A-Z0-9]+",
        str(package["analysis_set"]),
    ):
        raise RuntimeError(
            f"{package['screen_label']} AnalysisSet is not configured: "
            f"{package['analysis_set']}"
        )

dependency_specs = [
    (
        "shared GRCh38 / GENCODE v43 reference file",
        REFERENCE_FILE_INPUT,
    ),
    (
        "validation AnalysisSet",
        VALIDATION_PACKAGE["analysis_set"],
    ),
    (
        "genome-wide AnalysisSet",
        GENOME_WIDE_PACKAGE["analysis_set"],
    ),
    (
        "guide quantifications AnalysisStepVersion",
        GUIDE_QUANTIFICATION_ASV,
    ),
    (
        "differential element AnalysisStepVersion",
        DIFFERENTIAL_ELEMENT_ASV,
    ),
]

dependency_rows: list[dict] = []
for role, identifier in dependency_specs:
    try:
        record = get_with_retry(str(identifier))
        dependency_rows.append(
            {
                "role": role,
                "identifier": identifier,
                "found": bool(record),
                "accession": str((record or {}).get("accession", "")),
                "portal_id": str((record or {}).get("@id", "")),
                "status": str((record or {}).get("status", "")),
                "object_type": ",".join(
                    (record or {}).get("@type", [])
                ),
                "error": "",
            }
        )
    except Exception as error:
        dependency_rows.append(
            {
                "role": role,
                "identifier": identifier,
                "found": False,
                "accession": "",
                "portal_id": "",
                "status": "",
                "object_type": "",
                "error": f"{type(error).__name__}: {error}",
            }
        )

dependency_preflight = pd.DataFrame(dependency_rows)
DEPENDENCY_PREFLIGHT_PATH = (
    COMBINED_OUTPUT_ROOT / "qc/portal_dependency_preflight.tsv"
)
write_tsv(dependency_preflight, DEPENDENCY_PREFLIGHT_PATH)
print(dependency_preflight.to_string(index=False))

missing_dependencies = dependency_preflight.loc[
    ~dependency_preflight["found"].astype(bool)
].copy()

if not missing_dependencies.empty:
    details = []
    for row in missing_dependencies.to_dict("records"):
        error_text = str(row.get("error", "") or "not found")
        details.append(
            f"- {row['role']}: {row['identifier']} :: {error_text}"
        )
    raise RuntimeError(
        "Required Portal dependencies could not be read:\n"
        + "\n".join(details)
        + f"\nFull report: {DEPENDENCY_PREFLIGHT_PATH}"
    )

plan_rows: list[dict] = []

for item in combined_items:
    expected = item["expected_metadata"]
    lookup_identifier = str(
        item["existing_id"] or item["alias"]
    )

    try:
        record = get_with_retry(lookup_identifier)
        lookup_error = ""
    except Exception as error:
        record = None
        lookup_error = f"{type(error).__name__}: {error}"

    portal_md5 = str((record or {}).get("md5sum", ""))
    portal_file_set_raw = (record or {}).get("file_set", "")
    portal_file_set = compact_link_display(portal_file_set_raw)
    portal_content_type = str(
        (record or {}).get("content_type", "")
    )
    portal_file_format = str(
        (record or {}).get("file_format", "")
    )
    portal_asv_raw = (record or {}).get("analysis_step_version", "")
    portal_asv = compact_link_display(portal_asv_raw)
    portal_reference_files = (record or {}).get("reference_files", [])
    expected_reference_files = expected.get("reference_files", "")
    portal_upload_status = str(
        (record or {}).get("upload_status", "")
    )

    conflicts: list[str] = []
    if record:
        if portal_md5 != str(expected["md5sum"]):
            conflicts.append(
                f"md5sum portal={portal_md5!r} "
                f"expected={expected['md5sum']!r}"
            )
        if not links_equal(
            portal_file_set_raw,
            expected["file_set"],
        ):
            conflicts.append(
                "file_set portal="
                f"{compact_link_display(portal_file_set_raw)!r} "
                f"expected={expected['file_set']!r}"
            )
        if portal_content_type != str(
            expected["content_type"]
        ):
            conflicts.append(
                f"content_type portal={portal_content_type!r} "
                f"expected={expected['content_type']!r}"
            )
        if portal_file_format != str(expected["file_format"]):
            conflicts.append(
                f"file_format portal={portal_file_format!r} "
                f"expected={expected['file_format']!r}"
            )
        if not link_sets_equal(
            portal_reference_files,
            expected_reference_files,
        ):
            conflicts.append(
                "reference_files portal="
                f"{sorted(normalized_link_set(portal_reference_files))!r} "
                "expected="
                f"{sorted(normalized_link_set(expected_reference_files))!r}"
            )

    if lookup_error:
        action = "BLOCK_LOOKUP_ERROR"
    elif not record:
        action = "POST_NEW"
    elif conflicts:
        action = "BLOCK_EXISTING_RECORD_CONFLICT"
    elif portal_upload_status in {"file not found", "invalidated"}:
        action = "BLOCK_REUPLOAD_REVIEW_REQUIRED"
    elif links_equal(
        portal_asv_raw,
        expected["analysis_step_version"],
    ):
        action = "SKIP_ALREADY_CORRECT"
    else:
        action = "PATCH_ASV"

    plan_rows.append(
        {
            "screen_key": item["screen_key"],
            "screen_label": item["screen_label"],
            "role": item["role"],
            "alias": item["alias"],
            "lookup_identifier": lookup_identifier,
            "action": action,
            "portal_accession": str(
                (record or {}).get("accession", "")
            ),
            "portal_uuid": str(
                (record or {}).get("uuid", "")
            ),
            "portal_url": portal_record_url(record),
            "portal_upload_status": portal_upload_status,
            "portal_object_status": str(
                (record or {}).get("status", "")
            ),
            "expected_md5sum": str(expected["md5sum"]),
            "portal_md5sum": portal_md5,
            "expected_file_set": str(expected["file_set"]),
            "portal_file_set": compact_link_display(
                portal_file_set_raw
            ),
            "expected_content_type": str(
                expected["content_type"]
            ),
            "portal_content_type": portal_content_type,
            "expected_file_format": str(
                expected["file_format"]
            ),
            "portal_file_format": portal_file_format,
            "expected_analysis_step_version": str(
                expected["analysis_step_version"]
            ),
            "portal_analysis_step_version": compact_link_display(
                portal_asv_raw
            ),
            "expected_reference_files": ",".join(
                sorted(normalized_link_set(expected_reference_files))
            ),
            "portal_reference_files": ",".join(
                sorted(normalized_link_set(portal_reference_files))
            ),
            "conflicts": "; ".join(conflicts),
            "lookup_error": lookup_error,
            "submitted_file_name": str(
                item["final_file"]
            ),
        }
    )

submission_plan = pd.DataFrame(plan_rows)

# A matching existing record must not remain blocked solely because linked
# fields were returned as embedded objects.
false_link_conflicts = submission_plan.loc[
    submission_plan["portal_file_set"]
    .astype(str)
    .str.startswith("{")
]
if not false_link_conflicts.empty:
    raise RuntimeError(
        "Internal comparison error: an embedded file_set dictionary was "
        "serialized into the conflict report. Update the link parser before "
        "submitting."
    )

SUBMISSION_PLAN_PATH = (
    COMBINED_OUTPUT_ROOT / "metadata/igvf_submission_plan.tsv"
)
write_tsv(submission_plan, SUBMISSION_PLAN_PATH)

print("\nCombined submission plan:")
print(submission_plan.to_string(index=False))
print("\nSaved:", SUBMISSION_PLAN_PATH)

blocking_actions = {
    "BLOCK_LOOKUP_ERROR",
    "BLOCK_EXISTING_RECORD_CONFLICT",
    "BLOCK_REUPLOAD_REVIEW_REQUIRED",
}
blocked = submission_plan.loc[
    submission_plan["action"].isin(blocking_actions)
]
if not blocked.empty:
    raise RuntimeError(
        "The plan contains a blocking condition. Do not submit until "
        f"{SUBMISSION_PLAN_PATH} has been reviewed and resolved."
    )
```

### Text output (stream)

```text
                                      role                                                    identifier  found      accession                                                     portal_id      status              object_type error
shared GRCh38 / GENCODE v43 reference file                              /reference-files/IGVFFI7217ZMJZ/   True IGVFFI7217ZMJZ                              /reference-files/IGVFFI7217ZMJZ/    archived  ReferenceFile,File,Item      
                    validation AnalysisSet                                                IGVFDS4708PCVX   True IGVFDS4708PCVX                                /analysis-sets/IGVFDS4708PCVX/ in progress AnalysisSet,FileSet,Item      
                   genome-wide AnalysisSet                                                IGVFDS9629VXQK   True IGVFDS9629VXQK                                /analysis-sets/IGVFDS9629VXQK/ in progress AnalysisSet,FileSet,Item      
 guide quantifications AnalysisStepVersion /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/   True                /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/ in progress AnalysisStepVersion,Item      
  differential element AnalysisStepVersion /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/   True                /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/ in progress AnalysisStepVersion,Item      

Combined submission plan:
 screen_key              screen_label                                 role                                                                                                              alias                                                                                                  lookup_identifier               action portal_accession                          portal_uuid                                          portal_url portal_upload_status portal_object_status                  expected_md5sum                    portal_md5sum expected_file_set portal_file_set                expected_content_type   portal_content_type expected_file_format portal_file_format                                expected_analysis_step_version         portal_analysis_step_version expected_reference_files portal_reference_files conflicts lookup_error                                                                                                                                      submitted_file_name
 validation Focused validation screen                guide_quantifications                                       jesse-engreitz:tabular_file_WTC11-EC_validation_screen_guide_quantifications                                       jesse-engreitz:tabular_file_WTC11-EC_validation_screen_guide_quantifications SKIP_ALREADY_CORRECT   IGVFFI5157STFU 1f433411-ad80-4cbc-a2dd-c97275a7d881 https://data.igvf.org/tabular-files/IGVFFI5157STFU/              pending          in progress 0d4a1e607ef5f0da24f7fffa59e0d00c 0d4a1e607ef5f0da24f7fffa59e0d00c    IGVFDS4708PCVX  IGVFDS4708PCVX                guide quantifications guide quantifications                  tsv                tsv /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/ ce2957c8-747e-415b-91f7-dabd2730816d           IGVFFI7217ZMJZ         IGVFFI7217ZMJZ                                         /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_guide_quantifications.tsv.gz
 validation Focused validation screen differential_element_quantifications  jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications  jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications             POST_NEW                                                                                                                                                     2f0cb5fc62e1ba6d743cbabdef1cf8bf                                     IGVFDS4708PCVX                 differential element quantifications                                        tsv                    /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/                                                IGVFFI7217ZMJZ                                                 /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_differential_element_quantifications.tsv.gz
genome_wide        Genome-wide screen                guide_quantifications                                      jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications                                      jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications             POST_NEW                                                                                                                                                     e62c0a100589176793b61af60ed41032                                     IGVFDS9629VXQK                                guide quantifications                                        tsv                    /analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/                                                IGVFFI7217ZMJZ                                                              /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_guide_quantifications.tsv.gz
genome_wide        Genome-wide screen differential_element_quantifications jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications             POST_NEW                                                                                                                                                     f8d6aa0564246df55837da95c4888c60                                     IGVFDS9629VXQK                 differential element quantifications                                        tsv                    /analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/                                                IGVFFI7217ZMJZ                                               /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_differential_element_quantifications.tsv.gz

Saved: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/igvf_submission_plan.tsv

```

## Markdown cell 60

### C2A. Existing-record interpretation

A record with `pending` upload status is not automatically a conflict.
When all planned metadata fields match, its action should be
`SKIP_ALREADY_CORRECT`; downstream Portal validation can complete later.

## Code cell 61

```python
existing_record_summary = submission_plan.loc[
    submission_plan["portal_accession"].astype(str).ne(""),
    [
        "screen_key",
        "role",
        "portal_accession",
        "action",
        "portal_upload_status",
        "portal_object_status",
        "portal_file_set",
        "portal_analysis_step_version",
        "portal_reference_files",
        "conflicts",
    ],
].copy()

print(
    existing_record_summary.to_string(
        index=False,
        max_colwidth=120,
    )
    if not existing_record_summary.empty
    else "No existing file records were found."
)
```

### Text output (stream)

```text
screen_key                  role portal_accession               action portal_upload_status portal_object_status portal_file_set         portal_analysis_step_version portal_reference_files conflicts
validation guide_quantifications   IGVFFI5157STFU SKIP_ALREADY_CORRECT              pending          in progress  IGVFDS4708PCVX ce2957c8-747e-415b-91f7-dabd2730816d         IGVFFI7217ZMJZ          

```

## Markdown cell 62

## C3. Create one-row operation TSVs

One file is handled per operation. This makes partial success recoverable and keeps logs unambiguous.

## Code cell 63

```python
def safe_stem(value: str) -> str:
    return re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        value,
    ).strip("_")


def row_without_empty_values(row: dict) -> dict:
    return {
        key: value
        for key, value in row.items()
        if str(value).strip() != ""
    }


item_by_key = {
    (item["screen_key"], item["role"]): item
    for item in combined_items
}

ONE_ROW_DIR = (
    COMBINED_OUTPUT_ROOT / "work/igvf_submission_rows"
)
ONE_ROW_DIR.mkdir(parents=True, exist_ok=True)

operation_rows: list[dict] = []

for plan in submission_plan.to_dict("records"):
    key = (plan["screen_key"], plan["role"])
    item = item_by_key[key]
    expected = item["expected_metadata"]
    action = str(plan["action"])
    label = f"{plan['screen_key']}_{plan['role']}"

    if action == "POST_NEW":
        payload = row_without_empty_values(expected)
        tsv_path = ONE_ROW_DIR / f"POST_{safe_stem(label)}.tsv"
        write_tsv(pd.DataFrame([payload]), tsv_path)
    elif action == "PATCH_ASV":
        record_id = (
            str(plan["portal_accession"])
            or str(plan["portal_uuid"])
            or str(plan["lookup_identifier"])
        )
        payload = {
            "record_id": record_id,
            "analysis_step_version": str(
                plan["expected_analysis_step_version"]
            ),
        }
        tsv_path = ONE_ROW_DIR / f"PATCH_{safe_stem(label)}.tsv"
        write_tsv(pd.DataFrame([payload]), tsv_path)
    else:
        tsv_path = None

    operation_rows.append(
        {
            "screen_key": plan["screen_key"],
            "screen_label": plan["screen_label"],
            "role": plan["role"],
            "alias": plan["alias"],
            "lookup_identifier": plan["lookup_identifier"],
            "action": action,
            "input_tsv": str(tsv_path) if tsv_path else "",
        }
    )

submission_operations = pd.DataFrame(operation_rows)
SUBMISSION_OPERATIONS_PATH = (
    COMBINED_OUTPUT_ROOT
    / "metadata/igvf_submission_operations.tsv"
)
write_tsv(
    submission_operations,
    SUBMISSION_OPERATIONS_PATH,
)

print(submission_operations.to_string(index=False))
print("\nSaved:", SUBMISSION_OPERATIONS_PATH)
```

### Text output (stream)

```text
 screen_key              screen_label                                 role                                                                                                              alias                                                                                                  lookup_identifier               action                                                                                                                                                                 input_tsv
 validation Focused validation screen                guide_quantifications                                       jesse-engreitz:tabular_file_WTC11-EC_validation_screen_guide_quantifications                                       jesse-engreitz:tabular_file_WTC11-EC_validation_screen_guide_quantifications SKIP_ALREADY_CORRECT                                                                                                                                                                          
 validation Focused validation screen differential_element_quantifications  jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications  jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications             POST_NEW  /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/work/igvf_submission_rows/POST_validation_differential_element_quantifications.tsv
genome_wide        Genome-wide screen                guide_quantifications                                      jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications                                      jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications             POST_NEW                /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/work/igvf_submission_rows/POST_genome_wide_guide_quantifications.tsv
genome_wide        Genome-wide screen differential_element_quantifications jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications             POST_NEW /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/work/igvf_submission_rows/POST_genome_wide_differential_element_quantifications.tsv

Saved: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/igvf_submission_operations.tsv

```

## Markdown cell 64

## C4. Dry-run all required operations

Change `RUN_DRY_RUNS_NOW` to `True` after reviewing the operation plan.

The successful marker contains a fingerprint of all four final files, the plan, and every one-row TSV. Production submission is blocked if anything changes after the dry run.

## Code cell 65

```python
def sha256sum(
    path: Path,
    block_size: int = 16 * 1024 * 1024,
) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(block_size),
            b"",
        ):
            digest.update(block)
    return digest.hexdigest()


def submission_fingerprint() -> str:
    payload: dict[str, object] = {
        "plan": submission_plan.fillna("").to_dict("records"),
        "operations": [],
        "final_files": {},
    }
    for operation in submission_operations.to_dict("records"):
        input_tsv = str(operation.get("input_tsv", ""))
        if input_tsv:
            tsv_path = Path(input_tsv)
            payload["operations"].append(
                {
                    "screen_key": operation["screen_key"],
                    "role": operation["role"],
                    "action": operation["action"],
                    "input_tsv": str(tsv_path),
                    "sha256": sha256sum(tsv_path),
                }
            )
    for item in combined_items:
        file_path = Path(item["final_file"])
        payload["final_files"][str(file_path)] = (
            sha256sum(file_path)
        )

    serialized = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(
        serialized.encode("utf-8")
    ).hexdigest()


def run_iu_register_operation(
    input_tsv: Path,
    *,
    dry_run: bool,
    patch: bool,
    label: str,
) -> dict:
    if not IGVF_PYTHON.is_file():
        raise FileNotFoundError(IGVF_PYTHON)
    if not IU_REGISTER.is_file():
        raise FileNotFoundError(IU_REGISTER)
    if not input_tsv.is_file():
        raise FileNotFoundError(input_tsv)
    if not (
        os.environ.get("IGVF_API_KEY")
        and os.environ.get("IGVF_SECRET_KEY")
    ):
        raise RuntimeError("IGVF credentials are not set")

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    mode = "dryrun" if dry_run else "submit"
    log_path = (
        COMBINED_OUTPUT_ROOT
        / "logs"
        / f"{label}_{mode}_{timestamp}.log"
    )

    command = [
        str(IGVF_PYTHON),
        str(IU_REGISTER),
        "-m",
        "prod",
        "-p",
        "tabular_file",
    ]
    if dry_run:
        command.append("-d")
    if patch:
        command.append("--patch")
    command.extend(["-i", str(input_tsv)])

    print("\nCommand:", " ".join(command))
    print("Log:", log_path)

    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=os.environ.copy(),
        check=False,
    )
    log_path.write_text(
        completed.stdout,
        encoding="utf-8",
    )
    print(completed.stdout)

    fatal_patterns = [
        r"Traceback",
        r"ValidationFailure",
        r"Failed to POST",
        r"Failed to PATCH",
        r"401 Client Error",
        r"403 Client Error",
        r"409 Client Error",
        r"422 Client Error",
    ]
    detected = [
        pattern
        for pattern in fatal_patterns
        if re.search(
            pattern,
            completed.stdout,
            flags=re.IGNORECASE,
        )
    ]
    success = (
        completed.returncode == 0
        and not detected
    )

    result = {
        "label": label,
        "dry_run": dry_run,
        "patch": patch,
        "input_tsv": str(input_tsv),
        "log_path": str(log_path),
        "returncode": completed.returncode,
        "fatal_patterns": "; ".join(detected),
        "aws_upload_success_in_log": (
            "AWS upload successful" in completed.stdout
        ),
        "success": success,
    }
    if not success:
        raise RuntimeError(
            f"igvf_utils operation failed. Review {log_path}"
        )
    return result


RUN_DRY_RUNS_NOW = True
DRY_RUN_MARKER_PATH = (
    COMBINED_OUTPUT_ROOT / "DRY_RUN_PASSED_4_FILES.json"
)

if not RUN_DRY_RUNS_NOW:
    print(
        "Dry runs are disabled. Review the operation plan, change "
        "RUN_DRY_RUNS_NOW=True in this cell, and rerun it."
    )
else:
    dry_run_results: list[dict] = []

    for operation in submission_operations.to_dict(
        "records"
    ):
        action = str(operation["action"])
        label = (
            f"{operation['screen_key']}_"
            f"{operation['role']}"
        )

        if action == "SKIP_ALREADY_CORRECT":
            dry_run_results.append(
                {
                    "label": label,
                    "dry_run": True,
                    "patch": False,
                    "input_tsv": "",
                    "log_path": "",
                    "returncode": 0,
                    "fatal_patterns": "",
                    "aws_upload_success_in_log": False,
                    "success": True,
                    "note": (
                        "existing record already correct; "
                        "no operation required"
                    ),
                }
            )
            continue

        dry_run_results.append(
            run_iu_register_operation(
                Path(operation["input_tsv"]),
                dry_run=True,
                patch=(action == "PATCH_ASV"),
                label=label,
            )
        )

    dry_run_results_df = pd.DataFrame(
        dry_run_results
    )
    DRY_RUN_RESULTS_PATH = (
        COMBINED_OUTPUT_ROOT
        / "metadata/igvf_dry_run_results.tsv"
    )
    write_tsv(
        dry_run_results_df,
        DRY_RUN_RESULTS_PATH,
    )

    if not dry_run_results_df["success"].all():
        raise RuntimeError(
            "At least one dry run failed. Review "
            f"{DRY_RUN_RESULTS_PATH}"
        )

    marker = {
        "timestamp": time.strftime(
            "%Y-%m-%dT%H:%M:%S%z"
        ),
        "submission_fingerprint": (
            submission_fingerprint()
        ),
        "results_path": str(DRY_RUN_RESULTS_PATH),
    }
    DRY_RUN_MARKER_PATH.write_text(
        json.dumps(marker, indent=2) + "\n",
        encoding="utf-8",
    )

    print("All required dry runs passed.")
    print("Marker:", DRY_RUN_MARKER_PATH)
```

### Text output (stream)

```text

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file -d -i /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/work/igvf_submission_rows/POST_validation_differential_element_quantifications.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/logs/validation_differential_element_quantifications_dryrun_20260904_143031.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-04 14:30:58,248:iu_debug:	DRY RUN is enabled.
2026-09-04 14:30:58,258:iu_debug:	submission=False: In non-submission mode.
2026-09-04 14:30:58,258:iu_debug:	submission=True: In submission mode.
2026-09-04 14:30:59,763:iu_debug:	
IN post().
2026-09-04 14:30:59,763:iu_debug:	Validating the payload against the schema
2026-09-04 14:30:59,825:iu_debug:	<<<<<< POST tabular_file record jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications To IGVF database with URL https://api.data.igvf.org/tabular_file and this payload:

{
  "aliases": [
    "jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications"
  ],
  "analysis_step_version": "/analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/",
  "award": "/awards/HG011972/",
  "content_type": "differential element quantifications",
  "controlled_access": false,
  "derived_manually": false,
  "description": "Promoter-element-level MAGeCK results comparing Day 3 CDH5-positive and CDH5-negative fractions in the WTC11-EC validation CRISPRi screen. The original NOTCH1 gene-level statistics are duplicated unchanged across its P1 and P2 promoter elements following DACC guidance.",
  "file_format": "tsv",
  "file_set": "IGVFDS4708PCVX",
  "file_size": 16449,
  "filtered": false,
  "lab": "/labs/jesse-engreitz/",
  "md5sum": "2f0cb5fc62e1ba6d743cbabdef1cf8bf",
  "reference_files": [
    "/reference-files/IGVFFI7217ZMJZ/"
  ],
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_differential_element_quantifications.tsv.gz"
}


2026-09-04 14:30:59,826:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file -d -i /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/work/igvf_submission_rows/POST_genome_wide_guide_quantifications.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/logs/genome_wide_guide_quantifications_dryrun_20260904_143100.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-04 14:31:00,521:iu_debug:	DRY RUN is enabled.
2026-09-04 14:31:00,521:iu_debug:	submission=False: In non-submission mode.
2026-09-04 14:31:00,521:iu_debug:	submission=True: In submission mode.
2026-09-04 14:31:02,051:iu_debug:	
IN post().
2026-09-04 14:31:02,051:iu_debug:	Validating the payload against the schema
2026-09-04 14:31:02,112:iu_debug:	<<<<<< POST tabular_file record jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications To IGVF database with URL https://api.data.igvf.org/tabular_file and this payload:

{
  "aliases": [
    "jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications"
  ],
  "analysis_step_version": "/analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/",
  "award": "/awards/HG011972/",
  "content_type": "guide quantifications",
  "controlled_access": false,
  "derived_manually": false,
  "description": "Guide-level sequencing read counts for the genome-wide WTC11-EC CRISPRi differentiation screen in Day 0 baseline and Day 3 CDH5-positive and CDH5-negative fractions. Excel-converted target labels were repaired using the exact-matched authoritative library design.",
  "file_format": "tsv",
  "file_set": "IGVFDS9629VXQK",
  "file_size": 1646611,
  "filtered": false,
  "lab": "/labs/jesse-engreitz/",
  "md5sum": "e62c0a100589176793b61af60ed41032",
  "reference_files": [
    "/reference-files/IGVFFI7217ZMJZ/"
  ],
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_guide_quantifications.tsv.gz"
}


2026-09-04 14:31:02,112:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file -d -i /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/work/igvf_submission_rows/POST_genome_wide_differential_element_quantifications.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/logs/genome_wide_differential_element_quantifications_dryrun_20260904_143102.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-04 14:31:02,666:iu_debug:	DRY RUN is enabled.
2026-09-04 14:31:02,667:iu_debug:	submission=False: In non-submission mode.
2026-09-04 14:31:02,667:iu_debug:	submission=True: In submission mode.
2026-09-04 14:31:04,177:iu_debug:	
IN post().
2026-09-04 14:31:04,177:iu_debug:	Validating the payload against the schema
2026-09-04 14:31:04,239:iu_debug:	<<<<<< POST tabular_file record jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications To IGVF database with URL https://api.data.igvf.org/tabular_file and this payload:

{
  "aliases": [
    "jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications"
  ],
  "analysis_step_version": "/analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/",
  "award": "/awards/HG011972/",
  "content_type": "differential element quantifications",
  "controlled_access": false,
  "derived_manually": false,
  "description": "Native element-level MAGeCK results comparing Day 3 CDH5-positive and CDH5-negative fractions in the genome-wide WTC11-EC CRISPRi differentiation screen. Native promoter/transcript element identifiers and MAGeCK statistics are preserved.",
  "file_format": "tsv",
  "file_set": "IGVFDS9629VXQK",
  "file_size": 616932,
  "filtered": false,
  "lab": "/labs/jesse-engreitz/",
  "md5sum": "f8d6aa0564246df55837da95c4888c60",
  "reference_files": [
    "/reference-files/IGVFFI7217ZMJZ/"
  ],
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_differential_element_quantifications.tsv.gz"
}


2026-09-04 14:31:04,239:iu_debug:	DRY RUN is enabled.

All required dry runs passed.
Marker: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/DRY_RUN_PASSED_4_FILES.json

```

## Markdown cell 66

## C5. Submit or patch all four records, then reconcile

After reviewing all dry-run logs, set the confirmation token in this cell. The operations run sequentially in this order:

1. validation guide quantifications
2. validation differential element quantifications
3. genome-wide guide quantifications
4. genome-wide differential element quantifications

Existing matching records may be skipped or patched instead of re-POSTed.

## Code cell 67

```python
def verify_portal_file(
    item: dict,
    *,
    attempts: int = 12,
) -> dict:
    """Poll until the record exists and expected linked fields are visible."""
    expected = item["expected_metadata"]
    identifiers = [
        str(item.get("existing_id") or "").strip(),
        str(item["alias"]).strip(),
    ]
    identifiers = [
        value
        for index, value in enumerate(identifiers)
        if value and value not in identifiers[:index]
    ]

    last_record: dict | None = None
    last_issues: list[str] = ["record not found"]
    last_identifier = ""

    for attempt in range(1, attempts + 1):
        for identifier in identifiers:
            try:
                record = portal_get(identifier)
            except Exception:
                record = None

            if record is None:
                continue

            last_record = record
            last_identifier = identifier
            issues: list[str] = []

            if str(record.get("md5sum", "")) != str(
                expected["md5sum"]
            ):
                issues.append("md5sum mismatch")

            if not links_equal(
                record.get("file_set", ""),
                expected["file_set"],
            ):
                issues.append("file_set mismatch")

            if str(record.get("content_type", "")) != str(
                expected["content_type"]
            ):
                issues.append("content_type mismatch")

            if str(record.get("file_format", "")) != str(
                expected["file_format"]
            ):
                issues.append("file_format mismatch")

            if not links_equal(
                record.get("analysis_step_version", ""),
                expected["analysis_step_version"],
            ):
                issues.append("analysis_step_version mismatch")

            if not link_sets_equal(
                record.get("reference_files", []),
                expected.get("reference_files", ""),
            ):
                issues.append("reference_files mismatch")

            upload_status = str(record.get("upload_status", ""))
            if upload_status in {"file not found", "invalidated"}:
                issues.append(f"upload_status={upload_status}")

            last_issues = issues
            if not issues:
                break

        if last_record is not None and not last_issues:
            break

        if attempt < attempts:
            time.sleep(min(2 + attempt, 10))

    debug_dir = (
        COMBINED_OUTPUT_ROOT / "qc/portal_verification_debug"
    )
    debug_dir.mkdir(parents=True, exist_ok=True)
    debug_path = (
        debug_dir
        / f"{item['screen_key']}_{item['role']}.json"
    )
    debug_path.write_text(
        json.dumps(
            {
                "identifiers_tried": identifiers,
                "identifier_used": last_identifier,
                "expected_metadata": item["expected_metadata"],
                "portal_record": last_record or {},
                "issues": last_issues,
            },
            indent=2,
            default=str,
        )
        + "\\n",
        encoding="utf-8",
    )

    if last_record is None:
        return {
            "screen_key": item["screen_key"],
            "screen_label": item["screen_label"],
            "role": item["role"],
            "alias": item["alias"],
            "found": False,
            "accession": "",
            "portal_url": "",
            "portal_id": "",
            "uuid": "",
            "file_set": "",
            "md5sum": "",
            "content_type": "",
            "file_format": "",
            "analysis_step_version": "",
            "reference_files": "",
            "upload_status": "",
            "object_status": "",
            "file_size": "",
            "verification": "FAIL",
            "issues": "record not found",
            "debug_json": str(debug_path),
        }

    upload_status = str(last_record.get("upload_status", ""))
    return {
        "screen_key": item["screen_key"],
        "screen_label": item["screen_label"],
        "role": item["role"],
        "alias": item["alias"],
        "found": True,
        "accession": str(last_record.get("accession", "")),
        "portal_url": portal_record_url(last_record),
        "portal_id": str(last_record.get("@id", "")),
        "uuid": str(last_record.get("uuid", "")),
        "file_set": compact_link_display(
            last_record.get("file_set", "")
        ),
        "md5sum": str(last_record.get("md5sum", "")),
        "content_type": str(last_record.get("content_type", "")),
        "file_format": str(last_record.get("file_format", "")),
        "analysis_step_version": compact_link_display(
            last_record.get("analysis_step_version", "")
        ),
        "reference_files": ",".join(
            sorted(
                normalized_link_set(
                    last_record.get("reference_files", [])
                )
            )
        ),
        "upload_status": upload_status,
        "object_status": str(last_record.get("status", "")),
        "file_size": last_record.get("file_size", ""),
        "verification": (
            "PASS" if not last_issues else "FAIL"
        ),
        "issues": "; ".join(last_issues),
        "debug_json": str(debug_path),
    }

PRODUCTION_CONFIRMATION = "UPLOAD_BOTH_SCREENS_4_ANALYSIS_FILES_TO_IGVF_PROD"
REQUIRED_PRODUCTION_CONFIRMATION = (
    "UPLOAD_BOTH_SCREENS_4_ANALYSIS_FILES_TO_IGVF_PROD"
)

if (
    PRODUCTION_CONFIRMATION
    != REQUIRED_PRODUCTION_CONFIRMATION
):
    print("Production submission is blocked.")
    print(
        "After reviewing all dry-run logs, set "
        "PRODUCTION_CONFIRMATION to:"
    )
    print(REQUIRED_PRODUCTION_CONFIRMATION)
else:
    if not DRY_RUN_MARKER_PATH.is_file():
        raise RuntimeError(
            "The combined dry-run marker does not exist: "
            f"{DRY_RUN_MARKER_PATH}"
        )

    dry_run_marker = json.loads(
        DRY_RUN_MARKER_PATH.read_text(
            encoding="utf-8"
        )
    )
    if (
        dry_run_marker.get(
            "submission_fingerprint"
        )
        != submission_fingerprint()
    ):
        raise RuntimeError(
            "The files, plan, or operation TSVs changed after "
            "the dry run. Rerun the dry-run cell."
        )

    operation_by_key = {
        (
            operation["screen_key"],
            operation["role"],
        ): operation
        for operation in submission_operations.to_dict(
            "records"
        )
    }

    production_rows: list[dict] = []
    verification_rows: list[dict] = []

    for item in combined_items:
        key = (
            item["screen_key"],
            item["role"],
        )
        operation = operation_by_key[key]
        action = str(operation["action"])
        label = (
            f"{item['screen_key']}_"
            f"{item['role']}"
        )

        if action == "SKIP_ALREADY_CORRECT":
            production_rows.append(
                {
                    "screen_key": item[
                        "screen_key"
                    ],
                    "role": item["role"],
                    "action": action,
                    "operation_success": True,
                    "log_path": "",
                    "note": (
                        "existing record already matched"
                    ),
                }
            )
        elif action in {"POST_NEW", "PATCH_ASV"}:
            result = run_iu_register_operation(
                Path(operation["input_tsv"]),
                dry_run=False,
                patch=(action == "PATCH_ASV"),
                label=label,
            )
            production_rows.append(
                {
                    "screen_key": item[
                        "screen_key"
                    ],
                    "role": item["role"],
                    "action": action,
                    "operation_success": result[
                        "success"
                    ],
                    "log_path": result["log_path"],
                    "aws_upload_success_in_log": result[
                        "aws_upload_success_in_log"
                    ],
                    "note": "",
                }
            )
        else:
            raise RuntimeError(
                f"Unexpected production action: {action}"
            )

        verification = verify_portal_file(item)
        verification_rows.append(verification)
        print("\nPortal verification:")
        print(
            pd.DataFrame(
                [verification]
            ).to_string(index=False)
        )

        PARTIAL_VERIFICATION_PATH = (
            COMBINED_OUTPUT_ROOT
            / "metadata/igvf_verification_results_partial.csv"
        )
        pd.DataFrame(verification_rows).to_csv(
            PARTIAL_VERIFICATION_PATH,
            index=False,
            lineterminator="\n",
        )

        if verification["verification"] != "PASS":
            raise RuntimeError(
                f"Portal verification failed for {label}: "
                f"{verification.get('issues', 'unknown issue')}. "
                f"Review {PARTIAL_VERIFICATION_PATH}. "
                "Do not rerun the production POST blindly; rerun the Portal "
                "plan first so an existing record can be reconciled."
            )

    production_results = pd.DataFrame(
        production_rows
    )
    accession_results = pd.DataFrame(
        verification_rows
    )

    PRODUCTION_RESULTS_PATH = (
        COMBINED_OUTPUT_ROOT
        / "metadata/igvf_production_results.tsv"
    )
    ACCESSION_RESULTS_PATH = (
        COMBINED_OUTPUT_ROOT
        / "metadata/igvf_both_screens_accessions.csv"
    )

    write_tsv(
        production_results,
        PRODUCTION_RESULTS_PATH,
    )
    accession_results.to_csv(
        ACCESSION_RESULTS_PATH,
        index=False,
        lineterminator="\n",
    )

    if (
        len(accession_results) != 4
        or not accession_results[
            "verification"
        ].eq("PASS").all()
        or accession_results[
            "accession"
        ].astype(str).eq("").any()
    ):
        raise RuntimeError(
            "Final four-record reconciliation failed. "
            f"Review {ACCESSION_RESULTS_PATH}"
        )

    COMPLETE_MARKER = (
        COMBINED_OUTPUT_ROOT
        / "IGVF_SUBMISSION_COMPLETE_4_FILES.txt"
    )
    COMPLETE_MARKER.write_text(
        "All four validation/genome-wide screen analysis "
        "file records exist and match the combined local "
        "submission plan.\n"
        f"Accessions: {ACCESSION_RESULTS_PATH}\n",
        encoding="utf-8",
    )

    print("\nCombined production workflow completed.")
    print("Accessions:", ACCESSION_RESULTS_PATH)
    print("Marker:", COMPLETE_MARKER)
```

### Text output (stream)

```text

Portal verification:
screen_key              screen_label                  role                                                                        alias  found      accession                                          portal_url                      portal_id                                 uuid       file_set                           md5sum          content_type file_format                analysis_step_version reference_files upload_status object_status  file_size verification issues                                                                                                                                               debug_json
validation Focused validation screen guide_quantifications jesse-engreitz:tabular_file_WTC11-EC_validation_screen_guide_quantifications   True IGVFFI5157STFU https://data.igvf.org/tabular-files/IGVFFI5157STFU/ /tabular-files/IGVFFI5157STFU/ 1f433411-ad80-4cbc-a2dd-c97275a7d881 IGVFDS4708PCVX 0d4a1e607ef5f0da24f7fffa59e0d00c guide quantifications         tsv ce2957c8-747e-415b-91f7-dabd2730816d  IGVFFI7217ZMJZ       pending   in progress      40170         PASS        /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/qc/portal_verification_debug/validation_guide_quantifications.json

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file -i /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/work/igvf_submission_rows/POST_validation_differential_element_quantifications.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/logs/validation_differential_element_quantifications_submit_20260904_143107.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-04 14:31:07,812:iu_debug:	submission=False: In non-submission mode.
2026-09-04 14:31:07,812:iu_debug:	submission=True: In submission mode.
2026-09-04 14:31:09,166:iu_debug:	
IN post().
2026-09-04 14:31:09,166:iu_debug:	Validating the payload against the schema
2026-09-04 14:31:09,228:iu_debug:	<<<<<< POST tabular_file record jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications To IGVF database with URL https://api.data.igvf.org/tabular_file and this payload:

{
  "aliases": [
    "jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications"
  ],
  "analysis_step_version": "/analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/",
  "award": "/awards/HG011972/",
  "content_type": "differential element quantifications",
  "controlled_access": false,
  "derived_manually": false,
  "description": "Promoter-element-level MAGeCK results comparing Day 3 CDH5-positive and CDH5-negative fractions in the WTC11-EC validation CRISPRi screen. The original NOTCH1 gene-level statistics are duplicated unchanged across its P1 and P2 promoter elements following DACC guidance.",
  "file_format": "tsv",
  "file_set": "IGVFDS4708PCVX",
  "file_size": 16449,
  "filtered": false,
  "lab": "/labs/jesse-engreitz/",
  "md5sum": "2f0cb5fc62e1ba6d743cbabdef1cf8bf",
  "reference_files": [
    "/reference-files/IGVFFI7217ZMJZ/"
  ],
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_differential_element_quantifications.tsv.gz"
}


2026-09-04 14:31:09,839:iu_debug:	Success.
2026-09-04 14:31:09,840:iu_debug:	Object posted with identifier: IGVFFI3975KRIX
2026-09-04 14:31:09,855:iu_debug:	
IN upload_file()

2026-09-04 14:31:09,855:iu_debug:	Attempting to generate new file upload credentials
2026-09-04 14:31:10,245:iu_debug:	Success: upload credentials for 'IGVFFI3975KRIX' regenerated.
2026-09-04 14:31:10,245:iu_debug:	>>>>>>GET IGVFFI3975KRIX From DACC with URL https://api.data.igvf.org/IGVFFI3975KRIX/?format=json&datastore=database
2026-09-04 14:31:13,130:iu_debug:	Running command 'aws s3 cp /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_differential_element_quantifications.tsv.gz s3://igvf-files/2026/09/04/32dd02e0-4512-4e41-8c4c-6e8e2e4c9c4c/IGVFFI3975KRIX.tsv.gz'.

Uploading: [██████████████████████████████████████████████████] 100.00% done (16.1 KiB/16.1 KiB)
2026-09-04 14:31:16,259:iu_debug:	AWS upload successful.


Portal verification:
screen_key              screen_label                                 role                                                                                                             alias  found      accession                                          portal_url                      portal_id                                 uuid       file_set                           md5sum                         content_type file_format                analysis_step_version reference_files upload_status object_status  file_size verification issues                                                                                                                                                              debug_json
validation Focused validation screen differential_element_quantifications jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications   True IGVFFI3975KRIX https://data.igvf.org/tabular-files/IGVFFI3975KRIX/ /tabular-files/IGVFFI3975KRIX/ 32dd02e0-4512-4e41-8c4c-6e8e2e4c9c4c IGVFDS4708PCVX 2f0cb5fc62e1ba6d743cbabdef1cf8bf differential element quantifications         tsv 1ed82e06-e791-4ffd-86f5-eadaae67e8b4  IGVFFI7217ZMJZ       pending   in progress      16449         PASS        /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/qc/portal_verification_debug/validation_differential_element_quantifications.json

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file -i /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/work/igvf_submission_rows/POST_genome_wide_guide_quantifications.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/logs/genome_wide_guide_quantifications_submit_20260904_143118.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-04 14:31:19,461:iu_debug:	submission=False: In non-submission mode.
2026-09-04 14:31:19,461:iu_debug:	submission=True: In submission mode.
2026-09-04 14:31:20,814:iu_debug:	
IN post().
2026-09-04 14:31:20,814:iu_debug:	Validating the payload against the schema
2026-09-04 14:31:20,875:iu_debug:	<<<<<< POST tabular_file record jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications To IGVF database with URL https://api.data.igvf.org/tabular_file and this payload:

{
  "aliases": [
    "jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications"
  ],
  "analysis_step_version": "/analysis-step-versions/ce2957c8-747e-415b-91f7-dabd2730816d/",
  "award": "/awards/HG011972/",
  "content_type": "guide quantifications",
  "controlled_access": false,
  "derived_manually": false,
  "description": "Guide-level sequencing read counts for the genome-wide WTC11-EC CRISPRi differentiation screen in Day 0 baseline and Day 3 CDH5-positive and CDH5-negative fractions. Excel-converted target labels were repaired using the exact-matched authoritative library design.",
  "file_format": "tsv",
  "file_set": "IGVFDS9629VXQK",
  "file_size": 1646611,
  "filtered": false,
  "lab": "/labs/jesse-engreitz/",
  "md5sum": "e62c0a100589176793b61af60ed41032",
  "reference_files": [
    "/reference-files/IGVFFI7217ZMJZ/"
  ],
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_guide_quantifications.tsv.gz"
}


2026-09-04 14:31:21,302:iu_debug:	Success.
2026-09-04 14:31:21,302:iu_debug:	Object posted with identifier: IGVFFI8424MVZE
2026-09-04 14:31:21,302:iu_debug:	
IN upload_file()

2026-09-04 14:31:21,302:iu_debug:	Attempting to generate new file upload credentials
2026-09-04 14:31:21,627:iu_debug:	Success: upload credentials for 'IGVFFI8424MVZE' regenerated.
2026-09-04 14:31:21,627:iu_debug:	>>>>>>GET IGVFFI8424MVZE From DACC with URL https://api.data.igvf.org/IGVFFI8424MVZE/?format=json&datastore=database
2026-09-04 14:31:24,319:iu_debug:	Running command 'aws s3 cp /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_guide_quantifications.tsv.gz s3://igvf-files/2026/09/04/2925b588-dcf5-42d6-acd4-e2d9be677151/IGVFFI8424MVZE.tsv.gz'.

Uploading: [████████                                          ] 16.00% done (256.0 KiB/1.6 MiB)
Uploading: [████████████████                                  ] 32.00% done (512.0 KiB/1.6 MiB)
Uploading: [████████████████████████                          ] 48.00% done (768.0 KiB/1.6 MiB)
Uploading: [███████████████████████████████                   ] 62.50% done (1.0 MiB/1.6 MiB)  
Uploading: [█████████████████████████████████████             ] 75.00% done (1.2 MiB/1.6 MiB)  
Uploading: [██████████████████████████████████████████████    ] 93.75% done (1.5 MiB/1.6 MiB)  
Uploading: [██████████████████████████████████████████████████] 100.00% done (1.6 MiB/1.6 MiB) 
2026-09-04 14:31:25,661:iu_debug:	AWS upload successful.


Portal verification:
 screen_key       screen_label                  role                                                                         alias  found      accession                                          portal_url                      portal_id                                 uuid       file_set                           md5sum          content_type file_format                analysis_step_version reference_files upload_status object_status  file_size verification issues                                                                                                                                                debug_json
genome_wide Genome-wide screen guide_quantifications jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications   True IGVFFI8424MVZE https://data.igvf.org/tabular-files/IGVFFI8424MVZE/ /tabular-files/IGVFFI8424MVZE/ 2925b588-dcf5-42d6-acd4-e2d9be677151 IGVFDS9629VXQK e62c0a100589176793b61af60ed41032 guide quantifications         tsv ce2957c8-747e-415b-91f7-dabd2730816d  IGVFFI7217ZMJZ       pending   in progress    1646611         PASS        /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/qc/portal_verification_debug/genome_wide_guide_quantifications.json

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file -i /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/work/igvf_submission_rows/POST_genome_wide_differential_element_quantifications.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/logs/genome_wide_differential_element_quantifications_submit_20260904_143128.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-04 14:31:28,649:iu_debug:	submission=False: In non-submission mode.
2026-09-04 14:31:28,649:iu_debug:	submission=True: In submission mode.
2026-09-04 14:31:29,877:iu_debug:	
IN post().
2026-09-04 14:31:29,878:iu_debug:	Validating the payload against the schema
2026-09-04 14:31:29,939:iu_debug:	<<<<<< POST tabular_file record jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications To IGVF database with URL https://api.data.igvf.org/tabular_file and this payload:

{
  "aliases": [
    "jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications"
  ],
  "analysis_step_version": "/analysis-step-versions/1ed82e06-e791-4ffd-86f5-eadaae67e8b4/",
  "award": "/awards/HG011972/",
  "content_type": "differential element quantifications",
  "controlled_access": false,
  "derived_manually": false,
  "description": "Native element-level MAGeCK results comparing Day 3 CDH5-positive and CDH5-negative fractions in the genome-wide WTC11-EC CRISPRi differentiation screen. Native promoter/transcript element identifiers and MAGeCK statistics are preserved.",
  "file_format": "tsv",
  "file_set": "IGVFDS9629VXQK",
  "file_size": 616932,
  "filtered": false,
  "lab": "/labs/jesse-engreitz/",
  "md5sum": "f8d6aa0564246df55837da95c4888c60",
  "reference_files": [
    "/reference-files/IGVFFI7217ZMJZ/"
  ],
  "submitted_file_name": "/oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_differential_element_quantifications.tsv.gz"
}


2026-09-04 14:31:30,537:iu_debug:	Success.
2026-09-04 14:31:30,537:iu_debug:	Object posted with identifier: IGVFFI9001HIYI
2026-09-04 14:31:30,537:iu_debug:	
IN upload_file()

2026-09-04 14:31:30,537:iu_debug:	Attempting to generate new file upload credentials
2026-09-04 14:31:30,891:iu_debug:	Success: upload credentials for 'IGVFFI9001HIYI' regenerated.
2026-09-04 14:31:30,892:iu_debug:	>>>>>>GET IGVFFI9001HIYI From DACC with URL https://api.data.igvf.org/IGVFFI9001HIYI/?format=json&datastore=database
2026-09-04 14:31:33,349:iu_debug:	Running command 'aws s3 cp /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_differential_element_quantifications.tsv.gz s3://igvf-files/2026/09/04/34e6d451-d8ad-4089-9a4c-c198c8918f93/IGVFFI9001HIYI.tsv.gz'.

Uploading: [█████████████████████                             ] 42.49% done (256.0 KiB/602.5 KiB)
Uploading: [██████████████████████████████████████████        ] 84.98% done (512.0 KiB/602.5 KiB)
Uploading: [██████████████████████████████████████████████████] 100.00% done (602.5 KiB/602.5 KiB)
2026-09-04 14:31:34,550:iu_debug:	AWS upload successful.


Portal verification:
 screen_key       screen_label                                 role                                                                                                              alias  found      accession                                          portal_url                      portal_id                                 uuid       file_set                           md5sum                         content_type file_format                analysis_step_version reference_files upload_status object_status  file_size verification issues                                                                                                                                                               debug_json
genome_wide Genome-wide screen differential_element_quantifications jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications   True IGVFFI9001HIYI https://data.igvf.org/tabular-files/IGVFFI9001HIYI/ /tabular-files/IGVFFI9001HIYI/ 34e6d451-d8ad-4089-9a4c-c198c8918f93 IGVFDS9629VXQK f8d6aa0564246df55837da95c4888c60 differential element quantifications         tsv 1ed82e06-e791-4ffd-86f5-eadaae67e8b4  IGVFFI7217ZMJZ       pending   in progress     616932         PASS        /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/qc/portal_verification_debug/genome_wide_differential_element_quantifications.json

Combined production workflow completed.
Accessions: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/igvf_both_screens_accessions.csv
Marker: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/IGVF_SUBMISSION_COMPLETE_4_FILES.txt

```

## Markdown cell 68

### C6A. Read-only verification diagnostic

Run this cell after a production attempt stops at verification. It queries the validation guide-quantification alias and prints the exact failed fields without POSTing or PATCHing anything.

## Code cell 69

```python
validation_guide_item = next(
    item
    for item in combined_items
    if item["screen_key"] == "validation"
    and item["role"] == "guide_quantifications"
)
validation_guide_verification = verify_portal_file(
    validation_guide_item,
    attempts=3,
)
VALIDATION_GUIDE_VERIFICATION_PATH = (
    COMBINED_OUTPUT_ROOT
    / "metadata/validation_guide_quantifications_verification.csv"
)
pd.DataFrame([validation_guide_verification]).to_csv(
    VALIDATION_GUIDE_VERIFICATION_PATH,
    index=False,
    lineterminator="\n",
)
print(pd.DataFrame([validation_guide_verification]).T.to_string(header=False))
print("\nSaved:", VALIDATION_GUIDE_VERIFICATION_PATH)
```

### Text output (stream)

```text
screen_key                                                                                                                                                           validation
screen_label                                                                                                                                          Focused validation screen
role                                                                                                                                                      guide_quantifications
alias                                                                                              jesse-engreitz:tabular_file_WTC11-EC_validation_screen_guide_quantifications
found                                                                                                                                                                      True
accession                                                                                                                                                        IGVFFI5157STFU
portal_url                                                                                                                  https://data.igvf.org/tabular-files/IGVFFI5157STFU/
portal_id                                                                                                                                        /tabular-files/IGVFFI5157STFU/
uuid                                                                                                                                       1f433411-ad80-4cbc-a2dd-c97275a7d881
file_set                                                                                                                                                         IGVFDS4708PCVX
md5sum                                                                                                                                         0d4a1e607ef5f0da24f7fffa59e0d00c
content_type                                                                                                                                              guide quantifications
file_format                                                                                                                                                                 tsv
analysis_step_version                                                                                                                      ce2957c8-747e-415b-91f7-dabd2730816d
reference_files                                                                                                                                                  IGVFFI7217ZMJZ
upload_status                                                                                                                                                           pending
object_status                                                                                                                                                       in progress
file_size                                                                                                                                                                 40170
verification                                                                                                                                                               PASS
issues                                                                                                                                                                         
debug_json             /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/qc/portal_verification_debug/validation_guide_quantifications.json

Saved: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/validation_guide_quantifications_verification.csv

```

## Markdown cell 70

## C6. Refresh Portal status later

This cell is read-only. Rerun it after upload to monitor `upload_status` and downstream validation without resubmitting anything.

## Code cell 71

```python
status_rows: list[dict] = []

for item in combined_items:
    try:
        status_rows.append(
            verify_portal_file(
                item,
                attempts=4,
            )
        )
    except Exception as error:
        status_rows.append(
            {
                "screen_key": item["screen_key"],
                "screen_label": item[
                    "screen_label"
                ],
                "role": item["role"],
                "alias": item["alias"],
                "found": False,
                "verification": "FAIL",
                "issues": (
                    f"{type(error).__name__}: "
                    f"{error}"
                ),
            }
        )

portal_status = pd.DataFrame(status_rows)
PORTAL_STATUS_PATH = (
    COMBINED_OUTPUT_ROOT
    / "metadata/igvf_portal_status_latest.csv"
)
portal_status.to_csv(
    PORTAL_STATUS_PATH,
    index=False,
    lineterminator="\n",
)

print(portal_status.to_string(index=False))
print("\nSaved:", PORTAL_STATUS_PATH)

if (
    len(portal_status) == 4
    and portal_status["verification"].eq("PASS").all()
    and "upload_status" in portal_status.columns
    and portal_status["upload_status"].isin(
        {"validated", "validation exempted"}
    ).all()
):
    validated_marker = (
        COMBINED_OUTPUT_ROOT
        / "IGVF_FILES_VALIDATED_4_OF_4.txt"
    )
    validated_marker.write_text(
        "All four screen analysis files are validated "
        "or validation exempted on the IGVF Portal.\n",
        encoding="utf-8",
    )
    print(
        "All four files completed Portal validation:",
        validated_marker,
    )
else:
    print(
        "One or more records are absent, mismatched, or "
        "still awaiting downstream validation. Do not "
        "re-POST existing records; rerun this cell later."
    )
```

### Text output (stream)

```text
 screen_key              screen_label                                 role                                                                                                              alias  found      accession                                          portal_url                      portal_id                                 uuid       file_set                           md5sum                         content_type file_format                analysis_step_version reference_files upload_status object_status  file_size verification issues                                                                                                                                                               debug_json
 validation Focused validation screen                guide_quantifications                                       jesse-engreitz:tabular_file_WTC11-EC_validation_screen_guide_quantifications   True IGVFFI5157STFU https://data.igvf.org/tabular-files/IGVFFI5157STFU/ /tabular-files/IGVFFI5157STFU/ 1f433411-ad80-4cbc-a2dd-c97275a7d881 IGVFDS4708PCVX 0d4a1e607ef5f0da24f7fffa59e0d00c                guide quantifications         tsv ce2957c8-747e-415b-91f7-dabd2730816d  IGVFFI7217ZMJZ       pending   in progress      40170         PASS                        /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/qc/portal_verification_debug/validation_guide_quantifications.json
 validation Focused validation screen differential_element_quantifications  jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications   True IGVFFI3975KRIX https://data.igvf.org/tabular-files/IGVFFI3975KRIX/ /tabular-files/IGVFFI3975KRIX/ 32dd02e0-4512-4e41-8c4c-6e8e2e4c9c4c IGVFDS4708PCVX 2f0cb5fc62e1ba6d743cbabdef1cf8bf differential element quantifications         tsv 1ed82e06-e791-4ffd-86f5-eadaae67e8b4  IGVFFI7217ZMJZ       pending   in progress      16449         PASS         /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/qc/portal_verification_debug/validation_differential_element_quantifications.json
genome_wide        Genome-wide screen                guide_quantifications                                      jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications   True IGVFFI8424MVZE https://data.igvf.org/tabular-files/IGVFFI8424MVZE/ /tabular-files/IGVFFI8424MVZE/ 2925b588-dcf5-42d6-acd4-e2d9be677151 IGVFDS9629VXQK e62c0a100589176793b61af60ed41032                guide quantifications         tsv ce2957c8-747e-415b-91f7-dabd2730816d  IGVFFI7217ZMJZ       pending   in progress    1646611         PASS                       /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/qc/portal_verification_debug/genome_wide_guide_quantifications.json
genome_wide        Genome-wide screen differential_element_quantifications jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications   True IGVFFI9001HIYI https://data.igvf.org/tabular-files/IGVFFI9001HIYI/ /tabular-files/IGVFFI9001HIYI/ 34e6d451-d8ad-4089-9a4c-c198c8918f93 IGVFDS9629VXQK f8d6aa0564246df55837da95c4888c60 differential element quantifications         tsv 1ed82e06-e791-4ffd-86f5-eadaae67e8b4  IGVFFI7217ZMJZ       pending   in progress     616932         PASS        /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/qc/portal_verification_debug/genome_wide_differential_element_quantifications.json

Saved: /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/igvf_portal_status_latest.csv
One or more records are absent, mismatched, or still awaiting downstream validation. Do not re-POST existing records; rerun this cell later.

```

## Markdown cell 72

# Part D. Final combined output summary

## Code cell 73

```python
final_outputs = [
    (
        "validation guide quantifications",
        Path(
            VALIDATION_PACKAGE["items"][0][
                "final_file"
            ]
        ),
    ),
    (
        "validation differential elements",
        Path(
            VALIDATION_PACKAGE["items"][1][
                "final_file"
            ]
        ),
    ),
    (
        "genome-wide guide quantifications",
        Path(
            GENOME_WIDE_PACKAGE["items"][0][
                "final_file"
            ]
        ),
    ),
    (
        "genome-wide differential elements",
        Path(
            GENOME_WIDE_PACKAGE["items"][1][
                "final_file"
            ]
        ),
    ),
    ("combined manifest", COMBINED_MANIFEST_PATH),
    (
        "combined metadata review",
        COMBINED_METADATA_REVIEW_PATH,
    ),
    ("submission plan", SUBMISSION_PLAN_PATH),
    (
        "submission operations",
        SUBMISSION_OPERATIONS_PATH,
    ),
    (
        "dependency preflight",
        DEPENDENCY_PREFLIGHT_PATH,
    ),
    ("dry-run marker", DRY_RUN_MARKER_PATH),
    (
        "combined accession CSV",
        COMBINED_OUTPUT_ROOT
        / "metadata/igvf_both_screens_accessions.csv",
    ),
    (
        "latest Portal status CSV",
        COMBINED_OUTPUT_ROOT
        / "metadata/igvf_portal_status_latest.csv",
    ),
    (
        "completion marker",
        COMBINED_OUTPUT_ROOT
        / "IGVF_SUBMISSION_COMPLETE_4_FILES.txt",
    ),
]

final_summary = pd.DataFrame(
    [
        {
            "output": label,
            "path": str(path),
            "exists": path.is_file(),
            "size_bytes": (
                path.stat().st_size
                if path.is_file()
                else 0
            ),
        }
        for label, path in final_outputs
    ]
)

print(final_summary.to_string(index=False))
print("\nCombined accession output after submission:")
print(
    COMBINED_OUTPUT_ROOT
    / "metadata/igvf_both_screens_accessions.csv"
)
```

### Text output (stream)

```text
                           output                                                                                                                                                     path  exists  size_bytes
 validation guide quantifications                  /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_guide_quantifications.tsv.gz    True       40170
 validation differential elements   /oak/stanford/groups/engreitz/Users/tri/IGVF_validation_screen_analysis_submission/files/validation_screen_differential_element_quantifications.tsv.gz    True       16449
genome-wide guide quantifications                /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_guide_quantifications.tsv.gz    True     1646611
genome-wide differential elements /oak/stanford/groups/engreitz/Users/tri/IGVF_genome_wide_screen_analysis_submission/files/genome_wide_screen_differential_element_quantifications.tsv.gz    True      616932
                combined manifest                           /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/combined_four_file_manifest.tsv    True        3064
         combined metadata review                      /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/tabular_file_both_screens_review.tsv    True        3194
                  submission plan                                  /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/igvf_submission_plan.tsv    True        3012
            submission operations                            /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/igvf_submission_operations.tsv    True        1635
             dependency preflight                                 /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/qc/portal_dependency_preflight.tsv    True         904
                   dry-run marker                                        /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/DRY_RUN_PASSED_4_FILES.json    True         283
           combined accession CSV                          /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/igvf_both_screens_accessions.csv    True        2732
         latest Portal status CSV                             /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/igvf_portal_status_latest.csv    True        2732
                completion marker                               /oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/IGVF_SUBMISSION_COMPLETE_4_FILES.txt    True         253

Combined accession output after submission:
/oak/stanford/groups/engreitz/Users/tri/IGVF_both_crispri_screens_analysis_submission/metadata/igvf_both_screens_accessions.csv

```

## Code cell 74

```python

```

