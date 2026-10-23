# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/patch_crispri_macs_analysis_lineage.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/patch_crispri_macs_analysis_lineage.ipynb`
- SHA-256: `998c576bb0ff3864c5b72c5a742c8835b2dcf211b7a6bf8156c10fbdd2050ad9`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Patch CRISPRi MACS analysis lineage on the IGVF Portal

This notebook repairs two distinct lineage layers for the genome-wide and
focused-validation WTC11-EC CRISPRi MACS screens.

## Correct IGVF properties

An `AnalysisSet` does **not** use `derived_from`. Its upstream lineage is:

```text
AnalysisSet.input_file_sets → MeasurementSet records
```

Individual output files use:

```text
TabularFile.derived_from → input File records
```

The intended lineage is:

```text
GW AnalysisSet IGVFDS9629VXQK
  input_file_sets:
    IGVFDS2967ULMS
    IGVFDS8680ANMB
    IGVFDS8690AKQF

Validation AnalysisSet IGVFDS4708PCVX
  input_file_sets:
    IGVFDS1344YZHW
    IGVFDS1797KANA
    IGVFDS1832ZFCK
    IGVFDS5537WUFJ
    IGVFDS7114EJDN
    IGVFDS7624AWFB

GW guide quantifications IGVFFI8424MVZE
  derived_from: all 362 GW SequenceFiles

Validation guide quantifications IGVFFI5157STFU
  derived_from: all 41 validation SequenceFiles

GW differential quantifications IGVFFI9001HIYI
  derived_from: IGVFFI8424MVZE

Validation differential quantifications IGVFFI3975KRIX
  derived_from: IGVFFI5157STFU
```

The 403 raw SequenceFile accessions are read from the completed seqspec
reconciliation generated previously on Sherlock. The notebook refuses to
proceed unless all 403 rows are unique and have `verification = PASS`.

Existing array values are preserved. Only missing links are patched. By default,
unexpected existing lineage blocks submission.

## Markdown cell 1

## Version 2 — preserve additional AnalysisSet inputs

The initial notebook treated every existing `AnalysisSet.input_file_sets` value
outside the measurement-set list as a blocking conflict. That was too strict.

For an AnalysisSet, the listed screen MeasurementSets are **required inputs**, not
necessarily the complete universe of legitimate upstream file sets. This version
uses subset semantics:

- all required screen MeasurementSets must be present;
- additional existing input file sets are preserved and audited;
- no AnalysisSet PATCH is emitted when the required set is already present;
- `TabularFile.derived_from` remains strict and is checked against its exact
  expected lineage.

This prevents removal or accidental rejection of legitimate existing upstream
sets while retaining strict safeguards for output-file provenance.

## Markdown cell 2

## 1. Configuration

## Code cell 3

```python
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import time
from getpass import getpass
from pathlib import Path
from typing import Any

import pandas as pd
import requests

OAK_ROOT = Path(
    os.environ.get(
        "OAK",
        "/oak/stanford/groups/engreitz",
    )
)

SEQSPEC_OUTPUT_ROOT = (
    OAK_ROOT
    / "Users/tri/tmp/igvf_crispri_macs_seqspec"
)
SEQUENCE_RECONCILIATION_PATH = (
    SEQSPEC_OUTPUT_ROOT
    / "sequencefile_seqspec_reconciliation.csv"
)
IU_REGISTER_WRAPPER = (
    SEQSPEC_OUTPUT_ROOT
    / "_work/iu_register_github_tags_wrapper.py"
)

OUTPUT_ROOT = (
    OAK_ROOT
    / "Users/tri/tmp/igvf_crispri_macs_lineage_patch"
)
WORK_DIR = OUTPUT_ROOT / "_work"
LOG_DIR = OUTPUT_ROOT / "logs"
ROW_DIR = WORK_DIR / "one_row_json"

for directory in (
    OUTPUT_ROOT,
    WORK_DIR,
    LOG_DIR,
    ROW_DIR,
):
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )

IGVF_PYTHON = Path(
    "/oak/stanford/groups/engreitz/Users/tri/"
    "envs/igvf_utils/bin/python"
)
IU_REGISTER = Path(
    "/oak/stanford/groups/engreitz/Users/tri/"
    "envs/igvf_utils/igvf_utils/igvf_utils/"
    "MetaDataRegistration/iu_register.py"
)

ANALYSIS_SETS = {
    "genome_wide": {
        "analysis_set": "IGVFDS9629VXQK",
        "input_file_sets": [
            "IGVFDS2967ULMS",
            "IGVFDS8680ANMB",
            "IGVFDS8690AKQF",
        ],
    },
    "validation": {
        "analysis_set": "IGVFDS4708PCVX",
        "input_file_sets": [
            "IGVFDS1344YZHW",
            "IGVFDS1797KANA",
            "IGVFDS1832ZFCK",
            "IGVFDS5537WUFJ",
            "IGVFDS7114EJDN",
            "IGVFDS7624AWFB",
        ],
    },
}

OUTPUT_FILES = {
    "genome_wide": {
        "guide_quantifications": "IGVFFI8424MVZE",
        "differential_quantifications": "IGVFFI9001HIYI",
    },
    "validation": {
        "guide_quantifications": "IGVFFI5157STFU",
        "differential_quantifications": "IGVFFI3975KRIX",
    },
}

EXPECTED_SEQUENCE_FILE_COUNTS = {
    "genome_wide": 362,
    "validation": 41,
}

PATCH_ANALYSIS_SET_INPUT_FILE_SETS = True
PATCH_GUIDE_QUANTIFICATION_DERIVED_FROM = True
PATCH_DIFFERENTIAL_DERIVED_FROM = True

# AnalysisSet.input_file_sets uses required-subset semantics. Existing
# additional upstream file sets are preserved and audited.
ALLOW_ADDITIONAL_ANALYSIS_SET_INPUT_FILE_SETS = True

# Output-file derived_from lineage remains strict by default.
STRICT_NO_UNEXPECTED_FILE_DERIVED_FROM = True

RUN_DRY_RUNS = False
PRODUCTION_CONFIRMATION = ""
REQUIRED_PRODUCTION_CONFIRMATION = (
    "PATCH_CRISPRI_MACS_ANALYSIS_LINEAGE_TO_IGVF_PROD"
)

print("Seqspec output root:", SEQSPEC_OUTPUT_ROOT)
print("Reconciliation:", SEQUENCE_RECONCILIATION_PATH)
print("Patch output:", OUTPUT_ROOT)
```

### Text output (stream)

```text
Seqspec output root: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_seqspec
Reconciliation: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_seqspec/sequencefile_seqspec_reconciliation.csv
Patch output: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch

```

## Markdown cell 4

## 2. Utility functions

## Code cell 5

```python
def write_tsv(
    frame: pd.DataFrame,
    path: Path,
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
    )


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


def safe_stem(value: object) -> str:
    return re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        str(value),
    ).strip("_")


def split_link_items(
    value: object,
) -> list[object]:
    if value is None:
        return []
    if isinstance(value, str):
        return [
            item.strip()
            for item in value.split(",")
            if item.strip()
        ]
    if isinstance(
        value,
        (list, tuple, set),
    ):
        return list(value)
    return [value]


def terminal_identifier(
    value: object,
) -> str:
    if isinstance(value, dict):
        value = (
            value.get("accession")
            or value.get("@id")
            or value.get("uuid")
            or ""
        )
    text = str(value or "").strip().strip("/")
    return (
        text.rsplit("/", 1)[-1]
        if text
        else ""
    )


def direct_link_tokens(
    value: object,
) -> set[str]:
    tokens: set[str] = set()

    for item in split_link_items(value):
        if isinstance(item, dict):
            for key in (
                "accession",
                "uuid",
                "@id",
            ):
                token = terminal_identifier(
                    item.get(key, "")
                )
                if token:
                    tokens.add(token)

            for alias in split_link_items(
                item.get("aliases", [])
            ):
                alias_text = str(alias).strip()
                if alias_text:
                    tokens.add(alias_text)
        else:
            token = terminal_identifier(item)
            if token:
                tokens.add(token)

            text = str(item).strip()
            if (
                ":" in text
                and not text.startswith("/")
            ):
                tokens.add(text)

    return tokens


def accession_tokens(
    value: object,
    prefix: str,
) -> set[str]:
    return {
        token
        for token in direct_link_tokens(value)
        if token.startswith(prefix)
    }


def record_types(
    record: dict | None,
) -> set[str]:
    if not record:
        return set()

    value = record.get("@type", [])
    if isinstance(value, str):
        return {value}
    if isinstance(value, list):
        return {
            str(item)
            for item in value
        }
    return set()


def portal_auth() -> tuple[str, str] | None:
    if (
        os.environ.get("IGVF_API_KEY")
        and os.environ.get("IGVF_SECRET_KEY")
    ):
        return (
            os.environ["IGVF_API_KEY"],
            os.environ["IGVF_SECRET_KEY"],
        )
    return None


def portal_urls(
    identifier: str,
) -> list[str]:
    token = str(identifier).strip().strip("/")
    base = (
        "https://api.data.igvf.org/"
        f"{token}/"
    )
    return [
        (
            base
            + "?format=json&datastore=database"
            "&frame=object"
        ),
        (
            base
            + "?format=json&datastore=database"
            "&frame=edit"
        ),
        (
            base
            + "?format=json&datastore=database"
        ),
    ]


def portal_get(
    identifier: str,
    *,
    timeout: int = 60,
) -> dict | None:
    errors: list[str] = []

    for url in portal_urls(identifier):
        response = requests.get(
            url,
            auth=portal_auth(),
            headers={
                "Accept": "application/json"
            },
            timeout=timeout,
        )

        if response.status_code == 404:
            continue

        if response.ok:
            return response.json()

        errors.append(
            f"{response.status_code} {url}: "
            f"{response.text[:500]}"
        )

    if errors:
        raise RuntimeError(
            "Portal GET failed:\n"
            + "\n".join(errors)
        )

    return None


def get_with_retry(
    identifier: str,
    *,
    attempts: int = 5,
) -> dict | None:
    last_error: Exception | None = None

    for attempt in range(
        1,
        attempts + 1,
    ):
        try:
            record = portal_get(identifier)
            if record is not None:
                return record
        except Exception as error:
            last_error = error

        if attempt < attempts:
            time.sleep(
                min(
                    2 ** (attempt - 1),
                    10,
                )
            )

    if last_error is not None:
        raise last_error

    return None


def enter_igvf_credentials() -> None:
    api_key = <REDACTED>
        "IGVF API key: "
    ).strip()
    secret_key = getpass(
        "IGVF secret key: "
    ).strip()

    if not api_key or not secret_key:
        raise ValueError(
            "Both IGVF credentials are required"
        )

    os.environ[
        "IGVF_API_KEY"
    ] = api_key
    os.environ[
        "IGVF_SECRET_KEY"
    ] = secret_key

    del api_key, secret_key

    print(
        "IGVF credentials are set for "
        "the current notebook kernel."
    )


def detect_igvf_utils_local_version() -> str:
    probe_code = '''
from importlib.metadata import PackageNotFoundError, version

detected = ""
for distribution_name in ("igvf-utils", "igvf_utils"):
    try:
        detected = version(distribution_name)
        break
    except PackageNotFoundError:
        pass

if not detected:
    try:
        import igvf_utils
        detected = str(
            getattr(igvf_utils, "__version__", "")
        )
    except Exception:
        detected = ""

print(detected)
'''
    completed = subprocess.run(
        [
            str(IGVF_PYTHON),
            "-c",
            probe_code,
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=os.environ.copy(),
        check=False,
    )

    candidates = [
        line.strip().lstrip("v")
        for line in completed.stdout.splitlines()
        if re.fullmatch(
            r"v?\d+(?:\.\d+)+",
            line.strip(),
        )
    ]

    return (
        candidates[-1]
        if candidates
        else "3.1.1"
    )


def validate_submission_tools() -> None:
    for path in (
        IGVF_PYTHON,
        IU_REGISTER,
        IU_REGISTER_WRAPPER,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)

    compile_result = subprocess.run(
        [
            str(IGVF_PYTHON),
            "-m",
            "py_compile",
            str(IU_REGISTER_WRAPPER),
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=os.environ.copy(),
        check=False,
    )
    if compile_result.returncode != 0:
        raise RuntimeError(
            "Existing iu_register wrapper failed syntax validation:\n"
            + compile_result.stdout
        )


def run_iu_register(
    *,
    profile: str,
    input_path: Path,
    dry_run: bool,
    label: str,
) -> dict:
    validate_submission_tools()

    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if portal_auth() is None:
        raise RuntimeError(
            "IGVF credentials are not set"
        )

    local_version = (
        detect_igvf_utils_local_version()
    )
    timestamp = time.strftime(
        "%Y%m%d_%H%M%S"
    )
    mode = (
        "dryrun"
        if dry_run
        else "submit"
    )
    log_path = (
        LOG_DIR
        / (
            f"{safe_stem(label)}_"
            f"{mode}_{timestamp}.log"
        )
    )

    command = [
        str(IGVF_PYTHON),
        str(IU_REGISTER_WRAPPER),
        str(IU_REGISTER),
        "-m",
        "prod",
        "-p",
        profile,
        "--tries",
        "3",
        "--delay",
        "10",
        "--backoff",
        "2",
        "--patch",
    ]

    if dry_run:
        command.append("-d")

    command.extend(
        ["-i", str(input_path)]
    )

    environment = os.environ.copy()
    environment[
        "IGVF_UTILS_LOCAL_VERSION"
    ] = local_version

    print(
        "\nCommand:",
        " ".join(command),
    )
    print("Log:", log_path)

    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=environment,
        check=False,
    )

    log_path.write_text(
        completed.stdout,
        encoding="utf-8",
    )
    print(completed.stdout)

    fatal = re.search(
        (
            r"Traceback|ValidationFailure|"
            r"Failed to PATCH|"
            r"\b(?:401|403|409|422) "
            r"Client Error\b"
        ),
        completed.stdout,
        flags=re.IGNORECASE,
    )

    success = (
        completed.returncode == 0
        and fatal is None
    )

    if not success:
        raise RuntimeError(
            "igvf_utils failed; review "
            f"{log_path}"
        )

    return {
        "success": True,
        "log_path": str(log_path),
    }


print("Utility functions loaded.")
```

### Text output (stream)

```text
Utility functions loaded.

```

## Markdown cell 6

## 3. Load and validate the completed seqspec reconciliation

## Code cell 7

```python
if not SEQUENCE_RECONCILIATION_PATH.is_file():
    raise FileNotFoundError(
        SEQUENCE_RECONCILIATION_PATH
    )

sequence_reconciliation = pd.read_csv(
    SEQUENCE_RECONCILIATION_PATH,
    dtype=str,
    keep_default_na=False,
)

required_columns = {
    "screen",
    "file_set",
    "sequence_file_accession",
    "verification",
}
missing_columns = (
    required_columns
    - set(sequence_reconciliation.columns)
)
if missing_columns:
    raise RuntimeError(
        "Reconciliation is missing required columns: "
        f"{sorted(missing_columns)}"
    )

if not sequence_reconciliation[
    "verification"
].eq("PASS").all():
    bad = sequence_reconciliation.loc[
        ~sequence_reconciliation[
            "verification"
        ].eq("PASS")
    ]
    raise RuntimeError(
        "Not every seqspec reconciliation row passed:\n"
        + bad.to_string(index=False)
    )

if sequence_reconciliation[
    "sequence_file_accession"
].duplicated().any():
    duplicate_rows = sequence_reconciliation.loc[
        sequence_reconciliation[
            "sequence_file_accession"
        ].duplicated(
            keep=False
        )
    ]
    raise RuntimeError(
        "SequenceFile accessions are duplicated:\n"
        + duplicate_rows.to_string(index=False)
    )

sequence_files_by_screen = {
    screen: sorted(
        group[
            "sequence_file_accession"
        ].unique()
    )
    for screen, group in (
        sequence_reconciliation.groupby(
            "screen"
        )
    )
}

expected_screens = {
    "genome_wide",
    "validation",
}
if set(sequence_files_by_screen) != expected_screens:
    raise RuntimeError(
        "Unexpected reconciliation screens: "
        f"{sorted(sequence_files_by_screen)}"
    )

count_rows = []

for screen in sorted(expected_screens):
    accessions = sequence_files_by_screen[
        screen
    ]
    observed_file_sets = set(
        sequence_reconciliation.loc[
            sequence_reconciliation[
                "screen"
            ].eq(screen),
            "file_set",
        ].unique()
    )
    expected_file_sets = set(
        ANALYSIS_SETS[
            screen
        ]["input_file_sets"]
    )

    count_rows.append(
        {
            "screen": screen,
            "sequence_files": len(accessions),
            "expected_sequence_files": (
                EXPECTED_SEQUENCE_FILE_COUNTS[
                    screen
                ]
            ),
            "observed_file_sets": ",".join(
                sorted(observed_file_sets)
            ),
            "expected_file_sets": ",".join(
                sorted(expected_file_sets)
            ),
            "counts_match": (
                len(accessions)
                == EXPECTED_SEQUENCE_FILE_COUNTS[
                    screen
                ]
            ),
            "file_sets_match": (
                observed_file_sets
                == expected_file_sets
            ),
        }
    )

lineage_input_audit = pd.DataFrame(
    count_rows
)
LINEAGE_INPUT_AUDIT_PATH = (
    WORK_DIR
    / "lineage_input_audit.tsv"
)
write_tsv(
    lineage_input_audit,
    LINEAGE_INPUT_AUDIT_PATH,
)

if not (
    lineage_input_audit[
        "counts_match"
    ].all()
    and lineage_input_audit[
        "file_sets_match"
    ].all()
):
    raise RuntimeError(
        "The seqspec reconciliation does not match "
        f"the expected lineage. Review {LINEAGE_INPUT_AUDIT_PATH}"
    )

EXPECTED_SEQUENCE_FILES_PATH = (
    WORK_DIR
    / "expected_sequence_files.tsv"
)
write_tsv(
    sequence_reconciliation[
        [
            "screen",
            "file_set",
            "sequence_file_accession",
            "verification",
        ]
    ].sort_values(
        [
            "screen",
            "file_set",
            "sequence_file_accession",
        ]
    ),
    EXPECTED_SEQUENCE_FILES_PATH,
)

print(
    lineage_input_audit.to_string(
        index=False
    )
)
print(
    "\nVerified SequenceFiles:",
    len(sequence_reconciliation),
)
```

### Text output (stream)

```text
     screen  sequence_files  expected_sequence_files                                                                        observed_file_sets                                                                        expected_file_sets  counts_match  file_sets_match
genome_wide             362                      362                                              IGVFDS2967ULMS,IGVFDS8680ANMB,IGVFDS8690AKQF                                              IGVFDS2967ULMS,IGVFDS8680ANMB,IGVFDS8690AKQF          True             True
 validation              41                       41 IGVFDS1344YZHW,IGVFDS1797KANA,IGVFDS1832ZFCK,IGVFDS5537WUFJ,IGVFDS7114EJDN,IGVFDS7624AWFB IGVFDS1344YZHW,IGVFDS1797KANA,IGVFDS1832ZFCK,IGVFDS5537WUFJ,IGVFDS7114EJDN,IGVFDS7624AWFB          True             True

Verified SequenceFiles: 403

```

## Markdown cell 8

## 4. Define the expected lineage operations

## Code cell 9

```python
operation_specs: list[dict] = []

if PATCH_ANALYSIS_SET_INPUT_FILE_SETS:
    for screen, spec in (
        ANALYSIS_SETS.items()
    ):
        operation_specs.append(
            {
                "operation_key": (
                    f"{screen}_analysis_set_input_file_sets"
                ),
                "screen": screen,
                "profile": "analysis_set",
                "record_id": spec[
                    "analysis_set"
                ],
                "property": (
                    "input_file_sets"
                ),
                "expected_values": sorted(
                    spec[
                        "input_file_sets"
                    ]
                ),
                "expected_prefix": (
                    "IGVFDS"
                ),
                "comparison_mode": (
                    "required_subset"
                ),
            }
        )

if PATCH_GUIDE_QUANTIFICATION_DERIVED_FROM:
    for screen, spec in (
        OUTPUT_FILES.items()
    ):
        operation_specs.append(
            {
                "operation_key": (
                    f"{screen}_guide_quantifications_derived_from"
                ),
                "screen": screen,
                "profile": "tabular_file",
                "record_id": spec[
                    "guide_quantifications"
                ],
                "property": "derived_from",
                "expected_values": (
                    sequence_files_by_screen[
                        screen
                    ]
                ),
                "expected_prefix": (
                    "IGVFFI"
                ),
                "comparison_mode": "exact",
            }
        )

if PATCH_DIFFERENTIAL_DERIVED_FROM:
    for screen, spec in (
        OUTPUT_FILES.items()
    ):
        operation_specs.append(
            {
                "operation_key": (
                    f"{screen}_differential_derived_from"
                ),
                "screen": screen,
                "profile": "tabular_file",
                "record_id": spec[
                    "differential_quantifications"
                ],
                "property": "derived_from",
                "expected_values": [
                    spec[
                        "guide_quantifications"
                    ]
                ],
                "expected_prefix": (
                    "IGVFFI"
                ),
                "comparison_mode": "exact",
            }
        )

expected_lineage = pd.DataFrame(
    [
        {
            "operation_key": spec[
                "operation_key"
            ],
            "screen": spec["screen"],
            "profile": spec["profile"],
            "record_id": spec["record_id"],
            "property": spec["property"],
            "comparison_mode": spec[
                "comparison_mode"
            ],
            "expected_count": len(
                spec["expected_values"]
            ),
            "expected_values_preview": (
                ",".join(
                    spec["expected_values"][
                        :10
                    ]
                )
                + (
                    "..."
                    if len(
                        spec["expected_values"]
                    ) > 10
                    else ""
                )
            ),
        }
        for spec in operation_specs
    ]
)

EXPECTED_LINEAGE_PATH = (
    WORK_DIR
    / "expected_lineage_operations.tsv"
)
write_tsv(
    expected_lineage,
    EXPECTED_LINEAGE_PATH,
)

print(
    expected_lineage.to_string(
        index=False
    )
)
```

### Text output (stream)

```text
                                 operation_key      screen      profile      record_id        property comparison_mode  expected_count                                                                                                                                  expected_values_preview
      genome_wide_analysis_set_input_file_sets genome_wide analysis_set IGVFDS9629VXQK input_file_sets required_subset               3                                                                                                             IGVFDS2967ULMS,IGVFDS8680ANMB,IGVFDS8690AKQF
       validation_analysis_set_input_file_sets  validation analysis_set IGVFDS4708PCVX input_file_sets required_subset               6                                                                IGVFDS1344YZHW,IGVFDS1797KANA,IGVFDS1832ZFCK,IGVFDS5537WUFJ,IGVFDS7114EJDN,IGVFDS7624AWFB
genome_wide_guide_quantifications_derived_from genome_wide tabular_file IGVFFI8424MVZE    derived_from           exact             362 IGVFFI0082TOWX,IGVFFI0093GERJ,IGVFFI0239RKNC,IGVFFI0293IOCZ,IGVFFI0313IASP,IGVFFI0328LHVB,IGVFFI0371ZSBO,IGVFFI0406SABM,IGVFFI0424WQOS,IGVFFI0457GSVL...
 validation_guide_quantifications_derived_from  validation tabular_file IGVFFI5157STFU    derived_from           exact              41 IGVFFI0105AYPZ,IGVFFI0176SATY,IGVFFI0462CBMP,IGVFFI0804VFDP,IGVFFI0861WAON,IGVFFI1258TEQE,IGVFFI1554KQPJ,IGVFFI1760JPCW,IGVFFI2034YLTO,IGVFFI2090ROND...
         genome_wide_differential_derived_from genome_wide tabular_file IGVFFI9001HIYI    derived_from           exact               1                                                                                                                                           IGVFFI8424MVZE
          validation_differential_derived_from  validation tabular_file IGVFFI3975KRIX    derived_from           exact               1                                                                                                                                           IGVFFI5157STFU

```

## Markdown cell 10

## 5. Enter IGVF credentials

## Code cell 11

```python
enter_igvf_credentials()

print(
    "IGVF_API_KEY:",
    "SET"
    if os.environ.get(
        "IGVF_API_KEY"
    )
    else "MISSING",
)
print(
    "IGVF_SECRET_KEY:",
    "SET"
    if os.environ.get(
        "IGVF_SECRET_KEY"
    )
    else "MISSING",
)
```

### Text output (stream)

```text
IGVF credentials are set for the current notebook kernel.
IGVF_API_KEY: <REDACTED>
IGVF_SECRET_KEY: SET

```

## Markdown cell 12

## 6. Portal preflight and patch plan

The preflight checks record types, output-file membership, required inputs,
missing links, and additional existing links. It does not modify the Portal.

For `AnalysisSet.input_file_sets`, required MeasurementSets are evaluated as a
subset. Additional existing file sets are preserved and described in a separate
audit. For `TabularFile.derived_from`, exact lineage remains required.

## Code cell 13

```python
if portal_auth() is None:
    raise RuntimeError(
        "Run the credential-entry cell first"
    )

# -------------------------------------------------------------------------
# Verify the explicitly required screen MeasurementSets exist.
# -------------------------------------------------------------------------
upstream_file_set_rows = []

for screen, config in (
    ANALYSIS_SETS.items()
):
    for accession in (
        config["input_file_sets"]
    ):
        try:
            record = get_with_retry(
                accession
            )
            lookup_error = ""
        except Exception as error:
            record = None
            lookup_error = (
                f"{type(error).__name__}: "
                f"{error}"
            )

        upstream_file_set_rows.append(
            {
                "screen": screen,
                "accession": accession,
                "found": (
                    record is not None
                ),
                "record_types": ",".join(
                    sorted(
                        record_types(
                            record
                        )
                    )
                ),
                "aliases": ",".join(
                    sorted(
                        direct_link_tokens(
                            (record or {}).get(
                                "aliases",
                                [],
                            )
                        )
                    )
                ),
                "status": str(
                    (record or {}).get(
                        "status",
                        "",
                    )
                ),
                "lookup_error": (
                    lookup_error
                ),
            }
        )

upstream_file_set_audit = pd.DataFrame(
    upstream_file_set_rows
)
UPSTREAM_FILE_SET_AUDIT_PATH = (
    WORK_DIR
    / "upstream_file_set_audit.tsv"
)
write_tsv(
    upstream_file_set_audit,
    UPSTREAM_FILE_SET_AUDIT_PATH,
)

if not upstream_file_set_audit[
    "found"
].all():
    raise RuntimeError(
        "At least one required input file set was not found. "
        f"Review {UPSTREAM_FILE_SET_AUDIT_PATH}"
    )

# -------------------------------------------------------------------------
# Build the six-operation plan.
# -------------------------------------------------------------------------
plan_rows = []
additional_analysis_input_accessions: list[
    dict
] = []

for spec in operation_specs:
    record_id = spec["record_id"]

    try:
        record = get_with_retry(
            record_id
        )
        lookup_error = ""
    except Exception as error:
        record = None
        lookup_error = (
            f"{type(error).__name__}: "
            f"{error}"
        )

    action = ""
    existing_values: set[str] = set()

    if lookup_error:
        action = (
            "BLOCK_LOOKUP_ERROR"
        )
    elif record is None:
        action = (
            "BLOCK_RECORD_NOT_FOUND"
        )
    else:
        types = record_types(
            record
        )
        expected_type = (
            "AnalysisSet"
            if spec["profile"]
            == "analysis_set"
            else "TabularFile"
        )

        if expected_type not in types:
            action = (
                "BLOCK_WRONG_RECORD_TYPE"
            )

        existing_values = accession_tokens(
            record.get(
                spec["property"],
                [],
            ),
            spec[
                "expected_prefix"
            ],
        )

        if (
            not action
            and spec["profile"]
            == "tabular_file"
        ):
            expected_analysis_set = (
                ANALYSIS_SETS[
                    spec["screen"]
                ]["analysis_set"]
            )
            observed_file_sets = accession_tokens(
                record.get(
                    "file_set",
                    "",
                ),
                "IGVFDS",
            )
            if (
                observed_file_sets
                != {
                    expected_analysis_set
                }
            ):
                action = (
                    "BLOCK_WRONG_OUTPUT_FILE_SET"
                )

    expected_values = set(
        spec["expected_values"]
    )
    missing_values = (
        expected_values
        - existing_values
    )
    additional_values = (
        existing_values
        - expected_values
    )

    comparison_mode = spec[
        "comparison_mode"
    ]

    if not action:
        if comparison_mode == (
            "required_subset"
        ):
            if missing_values:
                action = (
                    "PATCH_MISSING_REQUIRED"
                )
            elif additional_values:
                action = (
                    "SKIP_REQUIRED_PRESENT_"
                    "WITH_ADDITIONAL_INPUTS"
                )
            else:
                action = (
                    "SKIP_ALREADY_CORRECT"
                )

            for accession in sorted(
                additional_values
            ):
                additional_analysis_input_accessions.append(
                    {
                        "screen": spec["screen"],
                        "analysis_set": (
                            spec["record_id"]
                        ),
                        "accession": accession,
                    }
                )

        elif comparison_mode == "exact":
            if (
                additional_values
                and STRICT_NO_UNEXPECTED_FILE_DERIVED_FROM
            ):
                action = (
                    "BLOCK_UNEXPECTED_EXISTING_LINKS"
                )
            elif missing_values:
                action = "PATCH_MISSING"
            else:
                action = (
                    "SKIP_ALREADY_CORRECT"
                )
        else:
            raise RuntimeError(
                "Unknown comparison mode: "
                f"{comparison_mode}"
            )

    plan_rows.append(
        {
            "operation_key": (
                spec["operation_key"]
            ),
            "screen": spec["screen"],
            "profile": spec["profile"],
            "record_id": record_id,
            "property": spec["property"],
            "comparison_mode": (
                comparison_mode
            ),
            "action": action,
            "expected_count": len(
                expected_values
            ),
            "existing_count": len(
                existing_values
            ),
            "missing_count": len(
                missing_values
            ),
            "additional_count": len(
                additional_values
            ),
            "missing_values": ",".join(
                sorted(missing_values)
            ),
            "additional_values": ",".join(
                sorted(
                    additional_values
                )
            ),
            "lookup_error": lookup_error,
        }
    )

lineage_patch_plan = pd.DataFrame(
    plan_rows
)
LINEAGE_PATCH_PLAN_PATH = (
    WORK_DIR
    / "lineage_patch_plan.tsv"
)
write_tsv(
    lineage_patch_plan,
    LINEAGE_PATCH_PLAN_PATH,
)

# -------------------------------------------------------------------------
# Resolve and describe additional AnalysisSet inputs. These are preserved;
# this audit is descriptive and does not create a PATCH.
# -------------------------------------------------------------------------
additional_input_rows = []

for item in (
    additional_analysis_input_accessions
):
    accession = item["accession"]

    try:
        record = get_with_retry(
            accession
        )
        lookup_error = ""
    except Exception as error:
        record = None
        lookup_error = (
            f"{type(error).__name__}: "
            f"{error}"
        )

    record = record or {}

    lab_tokens = direct_link_tokens(
        record.get("lab", "")
    )

    additional_input_rows.append(
        {
            **item,
            "found": bool(record),
            "record_types": ",".join(
                sorted(
                    record_types(
                        record
                    )
                )
            ),
            "aliases": ",".join(
                sorted(
                    direct_link_tokens(
                        record.get(
                            "aliases",
                            [],
                        )
                    )
                )
            ),
            "file_set_type": str(
                record.get(
                    "file_set_type",
                    "",
                )
            ),
            "title": str(
                record.get(
                    "title",
                    "",
                )
            ),
            "description": str(
                record.get(
                    "description",
                    "",
                )
            ),
            "status": str(
                record.get(
                    "status",
                    "",
                )
            ),
            "same_lab": (
                "jesse-engreitz"
                in lab_tokens
                or "/labs/jesse-engreitz/"
                in lab_tokens
            ),
            "lookup_error": (
                lookup_error
            ),
            "disposition": (
                "PRESERVE_EXISTING_INPUT"
            ),
        }
    )

additional_analysis_inputs_audit = pd.DataFrame(
    additional_input_rows,
    columns=[
        "screen",
        "analysis_set",
        "accession",
        "found",
        "record_types",
        "aliases",
        "file_set_type",
        "title",
        "description",
        "status",
        "same_lab",
        "lookup_error",
        "disposition",
    ],
)
ADDITIONAL_ANALYSIS_INPUTS_AUDIT_PATH = (
    WORK_DIR
    / "additional_analysis_input_file_sets_audit.tsv"
)
write_tsv(
    additional_analysis_inputs_audit,
    ADDITIONAL_ANALYSIS_INPUTS_AUDIT_PATH,
)

print(
    lineage_patch_plan[
        [
            "operation_key",
            "record_id",
            "property",
            "comparison_mode",
            "action",
            "expected_count",
            "existing_count",
            "missing_count",
            "additional_count",
        ]
    ].to_string(index=False)
)
print(
    "\nSaved:",
    LINEAGE_PATCH_PLAN_PATH,
)

if not additional_analysis_inputs_audit.empty:
    print(
        "\nAdditional AnalysisSet inputs "
        "(preserved, not patched):"
    )
    print(
        additional_analysis_inputs_audit[
            [
                "screen",
                "analysis_set",
                "accession",
                "record_types",
                "aliases",
                "file_set_type",
                "status",
                "found",
                "disposition",
            ]
        ].to_string(index=False)
    )
    print(
        "\nSaved:",
        ADDITIONAL_ANALYSIS_INPUTS_AUDIT_PATH,
    )

blocking_actions = {
    "BLOCK_LOOKUP_ERROR",
    "BLOCK_RECORD_NOT_FOUND",
    "BLOCK_WRONG_RECORD_TYPE",
    "BLOCK_WRONG_OUTPUT_FILE_SET",
    "BLOCK_UNEXPECTED_EXISTING_LINKS",
}

blocked = lineage_patch_plan.loc[
    lineage_patch_plan[
        "action"
    ].isin(
        blocking_actions
    )
]
if not blocked.empty:
    raise RuntimeError(
        "The lineage plan contains a blocking condition. "
        f"Review {LINEAGE_PATCH_PLAN_PATH}"
    )
```

### Text output (stream)

```text
                                 operation_key      record_id        property comparison_mode                                       action  expected_count  existing_count  missing_count  additional_count
      genome_wide_analysis_set_input_file_sets IGVFDS9629VXQK input_file_sets required_subset SKIP_REQUIRED_PRESENT_WITH_ADDITIONAL_INPUTS               3               4              0                 1
       validation_analysis_set_input_file_sets IGVFDS4708PCVX input_file_sets required_subset SKIP_REQUIRED_PRESENT_WITH_ADDITIONAL_INPUTS               6               7              0                 1
genome_wide_guide_quantifications_derived_from IGVFFI8424MVZE    derived_from           exact                                PATCH_MISSING             362               0            362                 0
 validation_guide_quantifications_derived_from IGVFFI5157STFU    derived_from           exact                                PATCH_MISSING              41               0             41                 0
         genome_wide_differential_derived_from IGVFFI9001HIYI    derived_from           exact                                PATCH_MISSING               1               0              1                 0
          validation_differential_derived_from IGVFFI3975KRIX    derived_from           exact                                PATCH_MISSING               1               0              1                 0

Saved: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/lineage_patch_plan.tsv

Additional AnalysisSet inputs (preserved, not patched):
     screen   analysis_set      accession            record_types                                                        aliases file_set_type      status  found             disposition
genome_wide IGVFDS9629VXQK IGVFDS4419BVVU CuratedSet,FileSet,Item            jesse-engreitz:curated_set_WTC11-EC_GWscreen_guides    guide RNAs in progress   True PRESERVE_EXISTING_INPUT
 validation IGVFDS4708PCVX IGVFDS8881NOFE CuratedSet,FileSet,Item jesse-engreitz:curated_set_WTC11-EC_GWscreen_validation_guides    guide RNAs in progress   True PRESERVE_EXISTING_INPUT

Saved: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/additional_analysis_input_file_sets_audit.tsv

```

## Markdown cell 14

## 7. Create native JSON PATCH inputs

Only missing links are included. `igvf_utils --patch` extends existing array
values and removes duplicate associations by default.

## Code cell 15

```python
operation_rows = []

for plan in lineage_patch_plan.to_dict(
    "records"
):
    if plan["action"] in {
        "PATCH_MISSING",
        "PATCH_MISSING_REQUIRED",
    }:
        missing_values = [
            value
            for value in str(
                plan["missing_values"]
            ).split(",")
            if value
        ]

        payload = {
            "record_id": (
                plan["record_id"]
            ),
            plan["property"]: (
                missing_values
            ),
        }

        json_path = (
            ROW_DIR
            / (
                "PATCH_"
                + safe_stem(
                    plan[
                        "operation_key"
                    ]
                )
                + ".json"
            )
        )
        json_path.write_text(
            json.dumps(
                payload,
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )

        reloaded = json.loads(
            json_path.read_text(
                encoding="utf-8"
            )
        )
        if reloaded != payload:
            raise RuntimeError(
                "JSON disk round trip changed "
                f"{json_path}"
            )
        if not isinstance(
            reloaded[
                plan["property"]
            ],
            list,
        ):
            raise RuntimeError(
                "PATCH association value must be a JSON array"
            )
    else:
        json_path = None

    operation_rows.append(
        {
            "operation_key": (
                plan["operation_key"]
            ),
            "screen": plan["screen"],
            "profile": plan["profile"],
            "record_id": (
                plan["record_id"]
            ),
            "property": (
                plan["property"]
            ),
            "action": plan["action"],
            "input_json": (
                str(json_path)
                if json_path
                else ""
            ),
            "input_sha256": (
                sha256sum(json_path)
                if json_path
                else ""
            ),
        }
    )

lineage_patch_operations = pd.DataFrame(
    operation_rows
)
LINEAGE_PATCH_OPERATIONS_PATH = (
    WORK_DIR
    / "lineage_patch_operations.tsv"
)
write_tsv(
    lineage_patch_operations,
    LINEAGE_PATCH_OPERATIONS_PATH,
)

print(
    lineage_patch_operations.to_string(
        index=False
    )
)
```

### Text output (stream)

```text
                                 operation_key      screen      profile      record_id        property                                       action                                                                                                                                               input_json                                                     input_sha256
      genome_wide_analysis_set_input_file_sets genome_wide analysis_set IGVFDS9629VXQK input_file_sets SKIP_REQUIRED_PRESENT_WITH_ADDITIONAL_INPUTS                                                                                                                                                                                                                          
       validation_analysis_set_input_file_sets  validation analysis_set IGVFDS4708PCVX input_file_sets SKIP_REQUIRED_PRESENT_WITH_ADDITIONAL_INPUTS                                                                                                                                                                                                                          
genome_wide_guide_quantifications_derived_from genome_wide tabular_file IGVFFI8424MVZE    derived_from                                PATCH_MISSING /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_genome_wide_guide_quantifications_derived_from.json 18a22dc28233e12dcfe559c6b321dd66962f057f84c2e90317c5f0fd141a10f7
 validation_guide_quantifications_derived_from  validation tabular_file IGVFFI5157STFU    derived_from                                PATCH_MISSING  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_validation_guide_quantifications_derived_from.json eca12529fba1e8768c42d7a213601741c9cbe9b6ac07490f1f7f1843bfc07b40
         genome_wide_differential_derived_from genome_wide tabular_file IGVFFI9001HIYI    derived_from                                PATCH_MISSING          /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_genome_wide_differential_derived_from.json f8a73da6c042479b235e782779130dcbf62903fa00802cb7e199df23c9ea510a
          validation_differential_derived_from  validation tabular_file IGVFFI3975KRIX    derived_from                                PATCH_MISSING           /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_validation_differential_derived_from.json dc0dd24f7fda25405850074d190d1ed4dafe2b9add823496f7b9c1b7911fd9d4

```

## Markdown cell 16

## 8. Dry-run all missing lineage PATCHes

## Code cell 17

```python
RUN_DRY_RUNS=True
```

## Code cell 18

```python
DRY_RUN_MARKER = (
    WORK_DIR
    / "LINEAGE_DRY_RUN_PASSED.json"
)

if not RUN_DRY_RUNS:
    print(
        "Dry runs are disabled. Set "
        "RUN_DRY_RUNS=True in Section 1 "
        "and rerun this cell."
    )
else:
    result_rows = []

    for operation in (
        lineage_patch_operations.to_dict(
            "records"
        )
    ):
        if str(
            operation["action"]
        ).startswith("SKIP_"):
            result_rows.append(
                {
                    **operation,
                    "success": True,
                    "log_path": "",
                }
            )
            continue

        if operation["action"] not in {
            "PATCH_MISSING",
            "PATCH_MISSING_REQUIRED",
        }:
            raise RuntimeError(
                "Unexpected operation action: "
                f"{operation}"
            )

        input_path = Path(
            operation["input_json"]
        )

        result = run_iu_register(
            profile=operation[
                "profile"
            ],
            input_path=input_path,
            dry_run=True,
            label=(
                "lineage_"
                + operation[
                    "operation_key"
                ]
            ),
        )
        result_rows.append(
            {
                **operation,
                **result,
            }
        )

    dry_run_results = pd.DataFrame(
        result_rows
    )
    DRY_RUN_RESULTS_PATH = (
        WORK_DIR
        / "lineage_dry_run_results.tsv"
    )
    write_tsv(
        dry_run_results,
        DRY_RUN_RESULTS_PATH,
    )

    if not dry_run_results[
        "success"
    ].all():
        raise RuntimeError(
            "At least one lineage dry run failed"
        )

    marker = {
        "timestamp": time.strftime(
            "%Y-%m-%dT%H:%M:%S%z"
        ),
        "plan_sha256": sha256sum(
            LINEAGE_PATCH_PLAN_PATH
        ),
        "operations_sha256": sha256sum(
            LINEAGE_PATCH_OPERATIONS_PATH
        ),
        "inputs": {
            row["operation_key"]: (
                row["input_sha256"]
            )
            for row in operation_rows
        },
    }

    DRY_RUN_MARKER.write_text(
        json.dumps(
            marker,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(
        dry_run_results[
            [
                "operation_key",
                "action",
                "success",
                "log_path",
            ]
        ].to_string(index=False)
    )
    print(
        "\nDry-run marker:",
        DRY_RUN_MARKER,
    )
```

### Text output (stream)

```text

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_seqspec/_work/iu_register_github_tags_wrapper.py /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -d -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_genome_wide_guide_quantifications_derived_from.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_genome_wide_guide_quantifications_derived_from_dryrun_20261012_182254.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.11, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-10-12 18:23:36,874:iu_debug:	DRY RUN is enabled.
2026-10-12 18:23:36,927:iu_debug:	submission=False: In non-submission mode.
2026-10-12 18:23:36,927:iu_debug:	submission=True: In submission mode.
2026-10-12 18:23:38,564:iu_debug:	
IN patch()
2026-10-12 18:23:38,564:iu_debug:	>>>>>>GET IGVFFI8424MVZE From DACC with URL https://api.data.igvf.org/IGVFFI8424MVZE/?format=json&datastore=database&frame=edit
2026-10-12 18:23:38,800:iu_debug:	<<<<<< PATCHING IGVFFI8424MVZE To IGVF database with URL https://api.data.igvf.org/IGVFFI8424MVZE and this payload:

{
  "derived_from": [
    "IGVFFI8673IIEA",
    "IGVFFI7572XTVQ",
    "IGVFFI1561IYDU",
    "IGVFFI3693RVEA",
    "IGVFFI3802HSST",
    "IGVFFI0518JEMU",
    "IGVFFI7621JZFU",
    "IGVFFI7700DYJX",
    "IGVFFI8682BEND",
    "IGVFFI1492LSCD",
    "IGVFFI7625IOYP",
    "IGVFFI4026IGXM",
    "IGVFFI8868ATZI",
    "IGVFFI6265XXSQ",
    "IGVFFI2429LKDQ",
    "IGVFFI3269KIJB",
    "IGVFFI2804VOWR",
    "IGVFFI6413YENU",
    "IGVFFI1442TUMQ",
    "IGVFFI6604GKMH",
    "IGVFFI7790RPQK",
    "IGVFFI7962EATF",
    "IGVFFI8942VUHE",
    "IGVFFI6421SDER",
    "IGVFFI3109BPDL",
    "IGVFFI3587EAFZ",
    "IGVFFI9234RPTH",
    "IGVFFI7731XWAS",
    "IGVFFI4795YAZT",
    "IGVFFI1979IQSD",
    "IGVFFI0942QWXR",
    "IGVFFI9062HYBX",
    "IGVFFI7382EOKJ",
    "IGVFFI0659HKZW",
    "IGVFFI4777DRZU",
    "IGVFFI8798ULOE",
    "IGVFFI2739CLDH",
    "IGVFFI3063RFAS",
    "IGVFFI4955UYQE",
    "IGVFFI4630NYYD",
    "IGVFFI7138QJSA",
    "IGVFFI5668CVCX",
    "IGVFFI3690GNVJ",
    "IGVFFI2954BRHB",
    "IGVFFI2123DKJW",
    "IGVFFI0789LJCV",
    "IGVFFI5792DDER",
    "IGVFFI7164SHRC",
    "IGVFFI8539TXWA",
    "IGVFFI9419APSX",
    "IGVFFI1842BDSA",
    "IGVFFI1923DCKJ",
    "IGVFFI2920ZKSE",
    "IGVFFI1018FGXW",
    "IGVFFI1390WGFS",
    "IGVFFI6018DRPD",
    "IGVFFI1370SKMJ",
    "IGVFFI4107VWXZ",
    "IGVFFI4950ATNP",
    "IGVFFI5998UCBD",
    "IGVFFI7626MBCD",
    "IGVFFI3899GUJQ",
    "IGVFFI9072KNCM",
    "IGVFFI3675XICB",
    "IGVFFI0627YLRE",
    "IGVFFI6853XNES",
    "IGVFFI7331RYPT",
    "IGVFFI0994SJXC",
    "IGVFFI5456GUKV",
    "IGVFFI1297WNUI",
    "IGVFFI4259MOJY",
    "IGVFFI0082TOWX",
    "IGVFFI8729MUNB",
    "IGVFFI3068MUVM",
    "IGVFFI3203SFBI",
    "IGVFFI1487JIHD",
    "IGVFFI6994WEQW",
    "IGVFFI4475HGVD",
    "IGVFFI2054COBI",
    "IGVFFI3960QQWW",
    "IGVFFI7994UAAL",
    "IGVFFI2750GSCT",
    "IGVFFI0525FJNT",
    "IGVFFI2092RIAE",
    "IGVFFI9581ZRAW",
    "IGVFFI8896JDBM",
    "IGVFFI2756OPOL",
    "IGVFFI7629WNXQ",
    "IGVFFI5641BENB",
    "IGVFFI9587MPBN",
    "IGVFFI0780GVXV",
    "IGVFFI9533SVQC",
    "IGVFFI5459QXVV",
    "IGVFFI2044ZUBM",
    "IGVFFI6225QFPG",
    "IGVFFI1528ELMI",
    "IGVFFI2735MNNL",
    "IGVFFI3371BSMR",
    "IGVFFI6007NBMX",
    "IGVFFI8175CDPB",
    "IGVFFI3001LFJF",
    "IGVFFI9383UQNX",
    "IGVFFI6303KCXL",
    "IGVFFI8836YHEZ",
    "IGVFFI2739XFYW",
    "IGVFFI7283RNDS",
    "IGVFFI3479XOMT",
    "IGVFFI1229GDQS",
    "IGVFFI4961BAEE",
    "IGVFFI2230PAJV",
    "IGVFFI6294EWAU",
    "IGVFFI0093GERJ",
    "IGVFFI1237AHUZ",
    "IGVFFI4686GKSV",
    "IGVFFI9428LBRW",
    "IGVFFI0940ZURT",
    "IGVFFI9295RNWF",
    "IGVFFI8473UFVO",
    "IGVFFI0810JXKW",
    "IGVFFI2601TLWG",
    "IGVFFI7030SFTV",
    "IGVFFI3321LEGN",
    "IGVFFI5459OIKB",
    "IGVFFI5535CFJT",
    "IGVFFI3805EUAD",
    "IGVFFI0935BQGZ",
    "IGVFFI4002JYRL",
    "IGVFFI1243CAAA",
    "IGVFFI8118EEKC",
    "IGVFFI1037CLKB",
    "IGVFFI1439TETW",
    "IGVFFI3904SVMC",
    "IGVFFI5315OJFY",
    "IGVFFI6175UOJJ",
    "IGVFFI6898QPRN",
    "IGVFFI1876KHTG",
    "IGVFFI1485BRHT",
    "IGVFFI8415WIPX",
    "IGVFFI3925DLJR",
    "IGVFFI3564DPKL",
    "IGVFFI6631REXU",
    "IGVFFI4035NLDE",
    "IGVFFI6227NESE",
    "IGVFFI2881ZQMN",
    "IGVFFI1750SKJU",
    "IGVFFI4483QCCF",
    "IGVFFI8009YJQR",
    "IGVFFI6068MBRJ",
    "IGVFFI0769ENCP",
    "IGVFFI4411YZKF",
    "IGVFFI8678UICK",
    "IGVFFI1065BYBZ",
    "IGVFFI3448TENP",
    "IGVFFI1686NAYC",
    "IGVFFI3688GQBS",
    "IGVFFI8921XTUV",
    "IGVFFI4658OPNJ",
    "IGVFFI3602QYHA",
    "IGVFFI4287UVBV",
    "IGVFFI7465KWEO",
    "IGVFFI8720UDCE",
    "IGVFFI2400DVYW",
    "IGVFFI5588VKZU",
    "IGVFFI8118XXKP",
    "IGVFFI3871CNOD",
    "IGVFFI2634CCJZ",
    "IGVFFI2777QCMF",
    "IGVFFI7281LNKI",
    "IGVFFI0654ZSZS",
    "IGVFFI3009OTUU",
    "IGVFFI5065TUPV",
    "IGVFFI9524TIOG",
    "IGVFFI1759JFXL",
    "IGVFFI0371ZSBO",
    "IGVFFI6575EJMW",
    "IGVFFI6764XZRC",
    "IGVFFI8275JTGO",
    "IGVFFI0868UXDF",
    "IGVFFI1216PYFR",
    "IGVFFI8997CAHG",
    "IGVFFI3560SJMG",
    "IGVFFI2127KLUS",
    "IGVFFI3724ODAQ",
    "IGVFFI4351BSDI",
    "IGVFFI1984AWKE",
    "IGVFFI4827RAXU",
    "IGVFFI4764FGFC",
    "IGVFFI4710AQLD",
    "IGVFFI7364QLHK",
    "IGVFFI2393UTVF",
    "IGVFFI1624AESM",
    "IGVFFI2064NOLW",
    "IGVFFI6576YSPO",
    "IGVFFI6804BGIW",
    "IGVFFI1161CZHZ",
    "IGVFFI7125HKYU",
    "IGVFFI3154ORLW",
    "IGVFFI5206YGMY",
    "IGVFFI0989VZAO",
    "IGVFFI1372MCJR",
    "IGVFFI4272ZOKZ",
    "IGVFFI3268AVJT",
    "IGVFFI6417WDPD",
    "IGVFFI1071UNQU",
    "IGVFFI1894RVCX",
    "IGVFFI8402UXJV",
    "IGVFFI3471MDCU",
    "IGVFFI6697OART",
    "IGVFFI5447LVOQ",
    "IGVFFI7488SFZQ",
    "IGVFFI2779BECB",
    "IGVFFI6040OXUF",
    "IGVFFI1229BEWL",
    "IGVFFI3035WMYR",
    "IGVFFI9277KXFI",
    "IGVFFI1857HULC",
    "IGVFFI6048HYGJ",
    "IGVFFI8332GSMS",
    "IGVFFI1510XNOZ",
    "IGVFFI7686ZOES",
    "IGVFFI8467QAYS",
    "IGVFFI1791GVUC",
    "IGVFFI0239RKNC",
    "IGVFFI7845QNBF",
    "IGVFFI9366LIOM",
    "IGVFFI3149SUSU",
    "IGVFFI5868XAWJ",
    "IGVFFI8155VPKE",
    "IGVFFI2555GINI",
    "IGVFFI2169ESWN",
    "IGVFFI8334XKJY",
    "IGVFFI6984HAQH",
    "IGVFFI6034KDPZ",
    "IGVFFI3044HJCV",
    "IGVFFI4358TRRL",
    "IGVFFI2908IDVT",
    "IGVFFI5752ICSY",
    "IGVFFI6492GCUJ",
    "IGVFFI2544MGSP",
    "IGVFFI5547WAII",
    "IGVFFI6569ZAQB",
    "IGVFFI4245UDRO",
    "IGVFFI2310UUSH",
    "IGVFFI1878AGGP",
    "IGVFFI2231XWRD",
    "IGVFFI5191YGON",
    "IGVFFI5565FMRO",
    "IGVFFI4315KHIY",
    "IGVFFI7948HBYB",
    "IGVFFI8279TMNX",
    "IGVFFI8198CNVT",
    "IGVFFI1503PINF",
    "IGVFFI3604YQSX",
    "IGVFFI4740KAZE",
    "IGVFFI8365PMGD",
    "IGVFFI0406SABM",
    "IGVFFI0843EEIV",
    "IGVFFI1237LSIF",
    "IGVFFI1618PRVJ",
    "IGVFFI4501AQQS",
    "IGVFFI2593IMTJ",
    "IGVFFI9638GDPY",
    "IGVFFI8579GPXN",
    "IGVFFI8260BRIF",
    "IGVFFI7536NXTR",
    "IGVFFI0313IASP",
    "IGVFFI3879QMRB",
    "IGVFFI0293IOCZ",
    "IGVFFI8340KYSN",
    "IGVFFI3072WKUI",
    "IGVFFI1322CXHE",
    "IGVFFI4902NFFK",
    "IGVFFI5819TLJZ",
    "IGVFFI6751OXDY",
    "IGVFFI5772VWCJ",
    "IGVFFI0538FVOR",
    "IGVFFI3438IPXE",
    "IGVFFI8676JGFJ",
    "IGVFFI7584SDHV",
    "IGVFFI5344ZWON",
    "IGVFFI1428IGVJ",
    "IGVFFI2002KHKE",
    "IGVFFI9820MPOG",
    "IGVFFI5867UHWK",
    "IGVFFI3136RIWG",
    "IGVFFI8316SDAJ",
    "IGVFFI5217YACS",
    "IGVFFI6460HZCB",
    "IGVFFI5508ETHY",
    "IGVFFI4076KPDZ",
    "IGVFFI8208XDCV",
    "IGVFFI4624BETR",
    "IGVFFI1038FPTO",
    "IGVFFI0701DTPP",
    "IGVFFI5779EZQS",
    "IGVFFI3594WTEB",
    "IGVFFI9467KZBG",
    "IGVFFI2636LFOD",
    "IGVFFI0807FKZN",
    "IGVFFI3180IRWR",
    "IGVFFI6406QIGN",
    "IGVFFI5323RJNF",
    "IGVFFI6951GRID",
    "IGVFFI7632GOZT",
    "IGVFFI1699GKVA",
    "IGVFFI8741VFEY",
    "IGVFFI2996VPKH",
    "IGVFFI4211MUCI",
    "IGVFFI9824HBGZ",
    "IGVFFI3339KFNV",
    "IGVFFI2935ZCQK",
    "IGVFFI1109BLNR",
    "IGVFFI1174AOBM",
    "IGVFFI9389CPTJ",
    "IGVFFI3412KBXJ",
    "IGVFFI4498PVIU",
    "IGVFFI5974AHOF",
    "IGVFFI2692HGWY",
    "IGVFFI7535BSIT",
    "IGVFFI3412UVWL",
    "IGVFFI5158DVZS",
    "IGVFFI1619DBBM",
    "IGVFFI1121TAQA",
    "IGVFFI8636NOVG",
    "IGVFFI4636COCN",
    "IGVFFI1264CUHK",
    "IGVFFI7888YLEY",
    "IGVFFI4441CDOK",
    "IGVFFI5623SIUE",
    "IGVFFI2048BBXB",
    "IGVFFI4799MAEX",
    "IGVFFI5033YLDE",
    "IGVFFI2094NZVJ",
    "IGVFFI6196CBAX",
    "IGVFFI6401RKMU",
    "IGVFFI3543LMYH",
    "IGVFFI4698MLLY",
    "IGVFFI8494CHGZ",
    "IGVFFI2290ASVP",
    "IGVFFI6313FKYW",
    "IGVFFI8046BEES",
    "IGVFFI5404LBVO",
    "IGVFFI9175WEZC",
    "IGVFFI0457GSVL",
    "IGVFFI1489XUIY",
    "IGVFFI0424WQOS",
    "IGVFFI1069BVHD",
    "IGVFFI3659KICP",
    "IGVFFI8014MUBO",
    "IGVFFI2927DJDA",
    "IGVFFI8072CJGA",
    "IGVFFI5472NRIN",
    "IGVFFI0328LHVB",
    "IGVFFI3981ZZOR",
    "IGVFFI8740LVAR",
    "IGVFFI9909DOXG",
    "IGVFFI4353AJJP",
    "IGVFFI5057YKXI",
    "IGVFFI1129PCNY",
    "IGVFFI6006ZKSI",
    "IGVFFI7657ILDG",
    "IGVFFI0580ZUFG"
  ]
}


2026-10-12 18:23:38,800:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_seqspec/_work/iu_register_github_tags_wrapper.py /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -d -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_validation_guide_quantifications_derived_from.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_validation_guide_quantifications_derived_from_dryrun_20261012_182338.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.11, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-10-12 18:23:39,478:iu_debug:	DRY RUN is enabled.
2026-10-12 18:23:39,478:iu_debug:	submission=False: In non-submission mode.
2026-10-12 18:23:39,478:iu_debug:	submission=True: In submission mode.
2026-10-12 18:23:40,994:iu_debug:	
IN patch()
2026-10-12 18:23:40,994:iu_debug:	>>>>>>GET IGVFFI5157STFU From DACC with URL https://api.data.igvf.org/IGVFFI5157STFU/?format=json&datastore=database&frame=edit
2026-10-12 18:23:41,312:iu_debug:	<<<<<< PATCHING IGVFFI5157STFU To IGVF database with URL https://api.data.igvf.org/IGVFFI5157STFU and this payload:

{
  "derived_from": [
    "IGVFFI7496LXLS",
    "IGVFFI5337YCKV",
    "IGVFFI8100DXSO",
    "IGVFFI8773VEGZ",
    "IGVFFI5755XWYU",
    "IGVFFI2034YLTO",
    "IGVFFI0105AYPZ",
    "IGVFFI2654MUEX",
    "IGVFFI9641OBRB",
    "IGVFFI0861WAON",
    "IGVFFI4885IXIT",
    "IGVFFI6713JZVV",
    "IGVFFI0176SATY",
    "IGVFFI6749WUQU",
    "IGVFFI5522LWIQ",
    "IGVFFI7405CGSC",
    "IGVFFI2090ROND",
    "IGVFFI3295AQRD",
    "IGVFFI6714NSWW",
    "IGVFFI1760JPCW",
    "IGVFFI6880HPPP",
    "IGVFFI7447DTDQ",
    "IGVFFI5740WAZW",
    "IGVFFI7579BYIA",
    "IGVFFI7967CUFR",
    "IGVFFI2178KIPK",
    "IGVFFI0804VFDP",
    "IGVFFI5548FGGY",
    "IGVFFI5122UWLR",
    "IGVFFI5586VNDU",
    "IGVFFI8104ORJP",
    "IGVFFI9606SGJV",
    "IGVFFI1258TEQE",
    "IGVFFI4426XNJZ",
    "IGVFFI1554KQPJ",
    "IGVFFI5712QWNG",
    "IGVFFI5050LOJP",
    "IGVFFI6469BUOI",
    "IGVFFI0462CBMP",
    "IGVFFI2240KPCD",
    "IGVFFI6112AOUD"
  ]
}


2026-10-12 18:23:41,312:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_seqspec/_work/iu_register_github_tags_wrapper.py /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -d -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_genome_wide_differential_derived_from.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_genome_wide_differential_derived_from_dryrun_20261012_182341.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.11, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-10-12 18:23:41,944:iu_debug:	DRY RUN is enabled.
2026-10-12 18:23:41,944:iu_debug:	submission=False: In non-submission mode.
2026-10-12 18:23:41,944:iu_debug:	submission=True: In submission mode.
2026-10-12 18:23:43,318:iu_debug:	
IN patch()
2026-10-12 18:23:43,318:iu_debug:	>>>>>>GET IGVFFI9001HIYI From DACC with URL https://api.data.igvf.org/IGVFFI9001HIYI/?format=json&datastore=database&frame=edit
2026-10-12 18:23:43,483:iu_debug:	<<<<<< PATCHING IGVFFI9001HIYI To IGVF database with URL https://api.data.igvf.org/IGVFFI9001HIYI and this payload:

{
  "derived_from": [
    "IGVFFI8424MVZE"
  ]
}


2026-10-12 18:23:43,483:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_seqspec/_work/iu_register_github_tags_wrapper.py /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -d -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_validation_differential_derived_from.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_validation_differential_derived_from_dryrun_20261012_182343.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.11, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-10-12 18:23:44,130:iu_debug:	DRY RUN is enabled.
2026-10-12 18:23:44,130:iu_debug:	submission=False: In non-submission mode.
2026-10-12 18:23:44,130:iu_debug:	submission=True: In submission mode.
2026-10-12 18:23:45,649:iu_debug:	
IN patch()
2026-10-12 18:23:45,649:iu_debug:	>>>>>>GET IGVFFI3975KRIX From DACC with URL https://api.data.igvf.org/IGVFFI3975KRIX/?format=json&datastore=database&frame=edit
2026-10-12 18:23:45,817:iu_debug:	<<<<<< PATCHING IGVFFI3975KRIX To IGVF database with URL https://api.data.igvf.org/IGVFFI3975KRIX and this payload:

{
  "derived_from": [
    "IGVFFI5157STFU"
  ]
}


2026-10-12 18:23:45,817:iu_debug:	DRY RUN is enabled.

                                 operation_key                                       action  success                                                                                                                                                           log_path
      genome_wide_analysis_set_input_file_sets SKIP_REQUIRED_PRESENT_WITH_ADDITIONAL_INPUTS     True                                                                                                                                                                   
       validation_analysis_set_input_file_sets SKIP_REQUIRED_PRESENT_WITH_ADDITIONAL_INPUTS     True                                                                                                                                                                   
genome_wide_guide_quantifications_derived_from                                PATCH_MISSING     True /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_genome_wide_guide_quantifications_derived_from_dryrun_20261012_182254.log
 validation_guide_quantifications_derived_from                                PATCH_MISSING     True  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_validation_guide_quantifications_derived_from_dryrun_20261012_182338.log
         genome_wide_differential_derived_from                                PATCH_MISSING     True          /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_genome_wide_differential_derived_from_dryrun_20261012_182341.log
          validation_differential_derived_from                                PATCH_MISSING     True           /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_validation_differential_derived_from_dryrun_20261012_182343.log

Dry-run marker: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/LINEAGE_DRY_RUN_PASSED.json

```

## Markdown cell 19

## 9. Apply the production lineage PATCHes

## Code cell 20

```python
PRODUCTION_CONFIRMATION = (
    "PATCH_CRISPRI_MACS_ANALYSIS_LINEAGE_TO_IGVF_PROD"
)
```

## Code cell 21

```python
if (
    PRODUCTION_CONFIRMATION
    != REQUIRED_PRODUCTION_CONFIRMATION
):
    print(
        "Production PATCH is blocked."
    )
    print(
        "Set PRODUCTION_CONFIRMATION to:"
    )
    print(
        REQUIRED_PRODUCTION_CONFIRMATION
    )
else:
    if not DRY_RUN_MARKER.is_file():
        raise RuntimeError(
            "Dry-run marker is missing"
        )

    marker = json.loads(
        DRY_RUN_MARKER.read_text(
            encoding="utf-8"
        )
    )

    if marker["plan_sha256"] != sha256sum(
        LINEAGE_PATCH_PLAN_PATH
    ):
        raise RuntimeError(
            "Patch plan changed after dry run"
        )

    if marker[
        "operations_sha256"
    ] != sha256sum(
        LINEAGE_PATCH_OPERATIONS_PATH
    ):
        raise RuntimeError(
            "Patch operations changed after dry run"
        )

    production_rows = []

    for operation in (
        lineage_patch_operations.to_dict(
            "records"
        )
    ):
        if str(
            operation["action"]
        ).startswith("SKIP_"):
            production_rows.append(
                {
                    **operation,
                    "success": True,
                    "log_path": "",
                }
            )
            continue

        input_path = Path(
            operation["input_json"]
        )
        expected_hash = marker[
            "inputs"
        ][operation["operation_key"]]

        if sha256sum(
            input_path
        ) != expected_hash:
            raise RuntimeError(
                "PATCH JSON changed after dry run: "
                f"{input_path}"
            )

        result = run_iu_register(
            profile=operation[
                "profile"
            ],
            input_path=input_path,
            dry_run=False,
            label=(
                "lineage_"
                + operation[
                    "operation_key"
                ]
            ),
        )
        production_rows.append(
            {
                **operation,
                **result,
            }
        )

    production_results = pd.DataFrame(
        production_rows
    )
    PRODUCTION_RESULTS_PATH = (
        WORK_DIR
        / "lineage_production_results.tsv"
    )
    write_tsv(
        production_results,
        PRODUCTION_RESULTS_PATH,
    )

    print(
        production_results[
            [
                "operation_key",
                "action",
                "success",
                "log_path",
            ]
        ].to_string(index=False)
    )
```

### Text output (stream)

```text

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_seqspec/_work/iu_register_github_tags_wrapper.py /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_genome_wide_guide_quantifications_derived_from.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_genome_wide_guide_quantifications_derived_from_submit_20261012_182447.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.11, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-10-12 18:24:48,586:iu_debug:	submission=False: In non-submission mode.
2026-10-12 18:24:48,586:iu_debug:	submission=True: In submission mode.
2026-10-12 18:24:49,965:iu_debug:	
IN patch()
2026-10-12 18:24:49,966:iu_debug:	>>>>>>GET IGVFFI8424MVZE From DACC with URL https://api.data.igvf.org/IGVFFI8424MVZE/?format=json&datastore=database&frame=edit
2026-10-12 18:24:50,149:iu_debug:	<<<<<< PATCHING IGVFFI8424MVZE To IGVF database with URL https://api.data.igvf.org/IGVFFI8424MVZE and this payload:

{
  "derived_from": [
    "IGVFFI6401RKMU",
    "IGVFFI5779EZQS",
    "IGVFFI6196CBAX",
    "IGVFFI1842BDSA",
    "IGVFFI7584SDHV",
    "IGVFFI4698MLLY",
    "IGVFFI0942QWXR",
    "IGVFFI0328LHVB",
    "IGVFFI9524TIOG",
    "IGVFFI2092RIAE",
    "IGVFFI4483QCCF",
    "IGVFFI1618PRVJ",
    "IGVFFI4498PVIU",
    "IGVFFI7962EATF",
    "IGVFFI0627YLRE",
    "IGVFFI7164SHRC",
    "IGVFFI5819TLJZ",
    "IGVFFI2231XWRD",
    "IGVFFI5158DVZS",
    "IGVFFI7700DYJX",
    "IGVFFI3339KFNV",
    "IGVFFI0654ZSZS",
    "IGVFFI3448TENP",
    "IGVFFI8676JGFJ",
    "IGVFFI6413YENU",
    "IGVFFI7632GOZT",
    "IGVFFI0518JEMU",
    "IGVFFI3659KICP",
    "IGVFFI0843EEIV",
    "IGVFFI9581ZRAW",
    "IGVFFI5315OJFY",
    "IGVFFI6492GCUJ",
    "IGVFFI9467KZBG",
    "IGVFFI4777DRZU",
    "IGVFFI8009YJQR",
    "IGVFFI3438IPXE",
    "IGVFFI4245UDRO",
    "IGVFFI8072CJGA",
    "IGVFFI4686GKSV",
    "IGVFFI4026IGXM",
    "IGVFFI2094NZVJ",
    "IGVFFI1237LSIF",
    "IGVFFI0659HKZW",
    "IGVFFI1442TUMQ",
    "IGVFFI3154ORLW",
    "IGVFFI5057YKXI",
    "IGVFFI7283RNDS",
    "IGVFFI3321LEGN",
    "IGVFFI1037CLKB",
    "IGVFFI0994SJXC",
    "IGVFFI3412KBXJ",
    "IGVFFI7948HBYB",
    "IGVFFI6048HYGJ",
    "IGVFFI1129PCNY",
    "IGVFFI2400DVYW",
    "IGVFFI3981ZZOR",
    "IGVFFI3604YQSX",
    "IGVFFI4351BSDI",
    "IGVFFI5641BENB",
    "IGVFFI3805EUAD",
    "IGVFFI8720UDCE",
    "IGVFFI1161CZHZ",
    "IGVFFI1018FGXW",
    "IGVFFI9638GDPY",
    "IGVFFI1492LSCD",
    "IGVFFI1264CUHK",
    "IGVFFI8334XKJY",
    "IGVFFI7030SFTV",
    "IGVFFI6417WDPD",
    "IGVFFI5535CFJT",
    "IGVFFI3693RVEA",
    "IGVFFI2908IDVT",
    "IGVFFI7465KWEO",
    "IGVFFI8473UFVO",
    "IGVFFI3180IRWR",
    "IGVFFI1121TAQA",
    "IGVFFI1984AWKE",
    "IGVFFI2593IMTJ",
    "IGVFFI3724ODAQ",
    "IGVFFI6697OART",
    "IGVFFI2920ZKSE",
    "IGVFFI8798ULOE",
    "IGVFFI3109BPDL",
    "IGVFFI7382EOKJ",
    "IGVFFI6994WEQW",
    "IGVFFI2555GINI",
    "IGVFFI1894RVCX",
    "IGVFFI1624AESM",
    "IGVFFI5998UCBD",
    "IGVFFI6421SDER",
    "IGVFFI3035WMYR",
    "IGVFFI6034KDPZ",
    "IGVFFI2692HGWY",
    "IGVFFI8579GPXN",
    "IGVFFI4961BAEE",
    "IGVFFI8155VPKE",
    "IGVFFI3879QMRB",
    "IGVFFI7731XWAS",
    "IGVFFI6751OXDY",
    "IGVFFI7535BSIT",
    "IGVFFI2290ASVP",
    "IGVFFI5472NRIN",
    "IGVFFI6040OXUF",
    "IGVFFI1297WNUI",
    "IGVFFI3203SFBI",
    "IGVFFI6984HAQH",
    "IGVFFI4002JYRL",
    "IGVFFI1237AHUZ",
    "IGVFFI5456GUKV",
    "IGVFFI5547WAII",
    "IGVFFI7629WNXQ",
    "IGVFFI8365PMGD",
    "IGVFFI9824HBGZ",
    "IGVFFI6175UOJJ",
    "IGVFFI6804BGIW",
    "IGVFFI7331RYPT",
    "IGVFFI4799MAEX",
    "IGVFFI4902NFFK",
    "IGVFFI2634CCJZ",
    "IGVFFI3802HSST",
    "IGVFFI6406QIGN",
    "IGVFFI7125HKYU",
    "IGVFFI8340KYSN",
    "IGVFFI2310UUSH",
    "IGVFFI9072KNCM",
    "IGVFFI3269KIJB",
    "IGVFFI9277KXFI",
    "IGVFFI9820MPOG",
    "IGVFFI6068MBRJ",
    "IGVFFI1038FPTO",
    "IGVFFI0457GSVL",
    "IGVFFI5065TUPV",
    "IGVFFI6303KCXL",
    "IGVFFI2636LFOD",
    "IGVFFI0780GVXV",
    "IGVFFI1619DBBM",
    "IGVFFI0810JXKW",
    "IGVFFI2044ZUBM",
    "IGVFFI4950ATNP",
    "IGVFFI1485BRHT",
    "IGVFFI3925DLJR",
    "IGVFFI8740LVAR",
    "IGVFFI1759JFXL",
    "IGVFFI1109BLNR",
    "IGVFFI9234RPTH",
    "IGVFFI0082TOWX",
    "IGVFFI5974AHOF",
    "IGVFFI6764XZRC",
    "IGVFFI0525FJNT",
    "IGVFFI4955UYQE",
    "IGVFFI8921XTUV",
    "IGVFFI3904SVMC",
    "IGVFFI8275JTGO",
    "IGVFFI2881ZQMN",
    "IGVFFI8678UICK",
    "IGVFFI8046BEES",
    "IGVFFI0538FVOR",
    "IGVFFI8539TXWA",
    "IGVFFI1439TETW",
    "IGVFFI6575EJMW",
    "IGVFFI9383UQNX",
    "IGVFFI1174AOBM",
    "IGVFFI6294EWAU",
    "IGVFFI7621JZFU",
    "IGVFFI3471MDCU",
    "IGVFFI8014MUBO",
    "IGVFFI7626MBCD",
    "IGVFFI3268AVJT",
    "IGVFFI5206YGMY",
    "IGVFFI2954BRHB",
    "IGVFFI6631REXU",
    "IGVFFI5868XAWJ",
    "IGVFFI2169ESWN",
    "IGVFFI1229GDQS",
    "IGVFFI9389CPTJ",
    "IGVFFI6225QFPG",
    "IGVFFI3688GQBS",
    "IGVFFI4358TRRL",
    "IGVFFI4636COCN",
    "IGVFFI5459OIKB",
    "IGVFFI3072WKUI",
    "IGVFFI4475HGVD",
    "IGVFFI3675XICB",
    "IGVFFI6460HZCB",
    "IGVFFI8118EEKC",
    "IGVFFI2002KHKE",
    "IGVFFI7281LNKI",
    "IGVFFI5565FMRO",
    "IGVFFI2429LKDQ",
    "IGVFFI8175CDPB",
    "IGVFFI5792DDER",
    "IGVFFI7364QLHK",
    "IGVFFI8682BEND",
    "IGVFFI5033YLDE",
    "IGVFFI7138QJSA",
    "IGVFFI7572XTVQ",
    "IGVFFI0406SABM",
    "IGVFFI3149SUSU",
    "IGVFFI1923DCKJ",
    "IGVFFI9366LIOM",
    "IGVFFI7536NXTR",
    "IGVFFI9428LBRW",
    "IGVFFI0293IOCZ",
    "IGVFFI1561IYDU",
    "IGVFFI8836YHEZ",
    "IGVFFI1510XNOZ",
    "IGVFFI5459QXVV",
    "IGVFFI8402UXJV",
    "IGVFFI5217YACS",
    "IGVFFI4501AQQS",
    "IGVFFI1487JIHD",
    "IGVFFI2064NOLW",
    "IGVFFI6007NBMX",
    "IGVFFI1503PINF",
    "IGVFFI0239RKNC",
    "IGVFFI2750GSCT",
    "IGVFFI4658OPNJ",
    "IGVFFI5508ETHY",
    "IGVFFI8316SDAJ",
    "IGVFFI4272ZOKZ",
    "IGVFFI8868ATZI",
    "IGVFFI6006ZKSI",
    "IGVFFI9533SVQC",
    "IGVFFI7488SFZQ",
    "IGVFFI6898QPRN",
    "IGVFFI4710AQLD",
    "IGVFFI3899GUJQ",
    "IGVFFI0989VZAO",
    "IGVFFI2756OPOL",
    "IGVFFI2601TLWG",
    "IGVFFI1229BEWL",
    "IGVFFI3044HJCV",
    "IGVFFI2779BECB",
    "IGVFFI2777QCMF",
    "IGVFFI8332GSMS",
    "IGVFFI5668CVCX",
    "IGVFFI3960QQWW",
    "IGVFFI5752ICSY",
    "IGVFFI3371BSMR",
    "IGVFFI7790RPQK",
    "IGVFFI1878AGGP",
    "IGVFFI8636NOVG",
    "IGVFFI0807FKZN",
    "IGVFFI5404LBVO",
    "IGVFFI7845QNBF",
    "IGVFFI8741VFEY",
    "IGVFFI3063RFAS",
    "IGVFFI0789LJCV",
    "IGVFFI3009OTUU",
    "IGVFFI4287UVBV",
    "IGVFFI7994UAAL",
    "IGVFFI5191YGON",
    "IGVFFI8198CNVT",
    "IGVFFI1065BYBZ",
    "IGVFFI2048BBXB",
    "IGVFFI6853XNES",
    "IGVFFI0940ZURT",
    "IGVFFI3136RIWG",
    "IGVFFI5772VWCJ",
    "IGVFFI4076KPDZ",
    "IGVFFI0935BQGZ",
    "IGVFFI3871CNOD",
    "IGVFFI2230PAJV",
    "IGVFFI1876KHTG",
    "IGVFFI1979IQSD",
    "IGVFFI3602QYHA",
    "IGVFFI1243CAAA",
    "IGVFFI5323RJNF",
    "IGVFFI2544MGSP",
    "IGVFFI6951GRID",
    "IGVFFI8208XDCV",
    "IGVFFI9587MPBN",
    "IGVFFI5623SIUE",
    "IGVFFI3690GNVJ",
    "IGVFFI6265XXSQ",
    "IGVFFI5588VKZU",
    "IGVFFI2054COBI",
    "IGVFFI1489XUIY",
    "IGVFFI8729MUNB",
    "IGVFFI9062HYBX",
    "IGVFFI2123DKJW",
    "IGVFFI4624BETR",
    "IGVFFI4441CDOK",
    "IGVFFI6569ZAQB",
    "IGVFFI7686ZOES",
    "IGVFFI1791GVUC",
    "IGVFFI3001LFJF",
    "IGVFFI4259MOJY",
    "IGVFFI0580ZUFG",
    "IGVFFI6313FKYW",
    "IGVFFI9295RNWF",
    "IGVFFI2393UTVF",
    "IGVFFI4211MUCI",
    "IGVFFI1216PYFR",
    "IGVFFI4795YAZT",
    "IGVFFI3594WTEB",
    "IGVFFI4630NYYD",
    "IGVFFI3564DPKL",
    "IGVFFI1750SKJU",
    "IGVFFI8118XXKP",
    "IGVFFI1428IGVJ",
    "IGVFFI1069BVHD",
    "IGVFFI3543LMYH",
    "IGVFFI8415WIPX",
    "IGVFFI4315KHIY",
    "IGVFFI8260BRIF",
    "IGVFFI1528ELMI",
    "IGVFFI0868UXDF",
    "IGVFFI2927DJDA",
    "IGVFFI1071UNQU",
    "IGVFFI3587EAFZ",
    "IGVFFI8494CHGZ",
    "IGVFFI3479XOMT",
    "IGVFFI7888YLEY",
    "IGVFFI5867UHWK",
    "IGVFFI8279TMNX",
    "IGVFFI0093GERJ",
    "IGVFFI2804VOWR",
    "IGVFFI9419APSX",
    "IGVFFI6227NESE",
    "IGVFFI1372MCJR",
    "IGVFFI0769ENCP",
    "IGVFFI5344ZWON",
    "IGVFFI4035NLDE",
    "IGVFFI4353AJJP",
    "IGVFFI8997CAHG",
    "IGVFFI2996VPKH",
    "IGVFFI1857HULC",
    "IGVFFI1686NAYC",
    "IGVFFI1370SKMJ",
    "IGVFFI4411YZKF",
    "IGVFFI4107VWXZ",
    "IGVFFI5447LVOQ",
    "IGVFFI2739CLDH",
    "IGVFFI1699GKVA",
    "IGVFFI9909DOXG",
    "IGVFFI6604GKMH",
    "IGVFFI8896JDBM",
    "IGVFFI6018DRPD",
    "IGVFFI8673IIEA",
    "IGVFFI3068MUVM",
    "IGVFFI8942VUHE",
    "IGVFFI4764FGFC",
    "IGVFFI2935ZCQK",
    "IGVFFI2739XFYW",
    "IGVFFI4740KAZE",
    "IGVFFI4827RAXU",
    "IGVFFI0371ZSBO",
    "IGVFFI0424WQOS",
    "IGVFFI0701DTPP",
    "IGVFFI3560SJMG",
    "IGVFFI1390WGFS",
    "IGVFFI2735MNNL",
    "IGVFFI7625IOYP",
    "IGVFFI3412UVWL",
    "IGVFFI1322CXHE",
    "IGVFFI8467QAYS",
    "IGVFFI0313IASP",
    "IGVFFI9175WEZC",
    "IGVFFI2127KLUS",
    "IGVFFI6576YSPO",
    "IGVFFI7657ILDG"
  ]
}


2026-10-12 18:24:52,340:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_seqspec/_work/iu_register_github_tags_wrapper.py /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_validation_guide_quantifications_derived_from.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_validation_guide_quantifications_derived_from_submit_20261012_182452.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.11, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-10-12 18:24:52,984:iu_debug:	submission=False: In non-submission mode.
2026-10-12 18:24:52,984:iu_debug:	submission=True: In submission mode.
2026-10-12 18:24:53,763:iu_debug:	
IN patch()
2026-10-12 18:24:53,763:iu_debug:	>>>>>>GET IGVFFI5157STFU From DACC with URL https://api.data.igvf.org/IGVFFI5157STFU/?format=json&datastore=database&frame=edit
2026-10-12 18:24:53,983:iu_debug:	<<<<<< PATCHING IGVFFI5157STFU To IGVF database with URL https://api.data.igvf.org/IGVFFI5157STFU and this payload:

{
  "derived_from": [
    "IGVFFI9641OBRB",
    "IGVFFI7579BYIA",
    "IGVFFI5122UWLR",
    "IGVFFI6112AOUD",
    "IGVFFI5522LWIQ",
    "IGVFFI6713JZVV",
    "IGVFFI9606SGJV",
    "IGVFFI8773VEGZ",
    "IGVFFI8104ORJP",
    "IGVFFI7405CGSC",
    "IGVFFI6880HPPP",
    "IGVFFI5548FGGY",
    "IGVFFI8100DXSO",
    "IGVFFI1258TEQE",
    "IGVFFI7496LXLS",
    "IGVFFI2034YLTO",
    "IGVFFI4885IXIT",
    "IGVFFI1554KQPJ",
    "IGVFFI3295AQRD",
    "IGVFFI5050LOJP",
    "IGVFFI7967CUFR",
    "IGVFFI2654MUEX",
    "IGVFFI4426XNJZ",
    "IGVFFI1760JPCW",
    "IGVFFI0804VFDP",
    "IGVFFI5740WAZW",
    "IGVFFI0105AYPZ",
    "IGVFFI5337YCKV",
    "IGVFFI7447DTDQ",
    "IGVFFI5586VNDU",
    "IGVFFI0462CBMP",
    "IGVFFI2240KPCD",
    "IGVFFI6749WUQU",
    "IGVFFI2178KIPK",
    "IGVFFI6714NSWW",
    "IGVFFI0861WAON",
    "IGVFFI5712QWNG",
    "IGVFFI0176SATY",
    "IGVFFI6469BUOI",
    "IGVFFI5755XWYU",
    "IGVFFI2090ROND"
  ]
}


2026-10-12 18:24:54,563:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_seqspec/_work/iu_register_github_tags_wrapper.py /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_genome_wide_differential_derived_from.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_genome_wide_differential_derived_from_submit_20261012_182454.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.11, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-10-12 18:24:55,223:iu_debug:	submission=False: In non-submission mode.
2026-10-12 18:24:55,224:iu_debug:	submission=True: In submission mode.
2026-10-12 18:24:56,536:iu_debug:	
IN patch()
2026-10-12 18:24:56,536:iu_debug:	>>>>>>GET IGVFFI9001HIYI From DACC with URL https://api.data.igvf.org/IGVFFI9001HIYI/?format=json&datastore=database&frame=edit
2026-10-12 18:24:56,710:iu_debug:	<<<<<< PATCHING IGVFFI9001HIYI To IGVF database with URL https://api.data.igvf.org/IGVFFI9001HIYI and this payload:

{
  "derived_from": [
    "IGVFFI8424MVZE"
  ]
}


2026-10-12 18:24:57,360:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_seqspec/_work/iu_register_github_tags_wrapper.py /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/one_row_json/PATCH_validation_differential_derived_from.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_validation_differential_derived_from_submit_20261012_182457.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.11, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-10-12 18:24:58,073:iu_debug:	submission=False: In non-submission mode.
2026-10-12 18:24:58,073:iu_debug:	submission=True: In submission mode.
2026-10-12 18:24:59,336:iu_debug:	
IN patch()
2026-10-12 18:24:59,336:iu_debug:	>>>>>>GET IGVFFI3975KRIX From DACC with URL https://api.data.igvf.org/IGVFFI3975KRIX/?format=json&datastore=database&frame=edit
2026-10-12 18:24:59,506:iu_debug:	<<<<<< PATCHING IGVFFI3975KRIX To IGVF database with URL https://api.data.igvf.org/IGVFFI3975KRIX and this payload:

{
  "derived_from": [
    "IGVFFI5157STFU"
  ]
}


2026-10-12 18:25:00,014:iu_debug:	Success.

                                 operation_key                                       action  success                                                                                                                                                           log_path
      genome_wide_analysis_set_input_file_sets SKIP_REQUIRED_PRESENT_WITH_ADDITIONAL_INPUTS     True                                                                                                                                                                   
       validation_analysis_set_input_file_sets SKIP_REQUIRED_PRESENT_WITH_ADDITIONAL_INPUTS     True                                                                                                                                                                   
genome_wide_guide_quantifications_derived_from                                PATCH_MISSING     True /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_genome_wide_guide_quantifications_derived_from_submit_20261012_182447.log
 validation_guide_quantifications_derived_from                                PATCH_MISSING     True  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_validation_guide_quantifications_derived_from_submit_20261012_182452.log
         genome_wide_differential_derived_from                                PATCH_MISSING     True          /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_genome_wide_differential_derived_from_submit_20261012_182454.log
          validation_differential_derived_from                                PATCH_MISSING     True           /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs/lineage_validation_differential_derived_from_submit_20261012_182457.log

```

## Markdown cell 22

## 10. Final read-only reconciliation

## Code cell 23

```python
final_rows = []

for spec in operation_specs:
    record = get_with_retry(
        spec["record_id"],
        attempts=6,
    )

    if record is None:
        raise RuntimeError(
            "Record missing during final reconciliation: "
            f"{spec['record_id']}"
        )

    observed = accession_tokens(
        record.get(
            spec["property"],
            [],
        ),
        spec["expected_prefix"],
    )
    expected = set(
        spec["expected_values"]
    )

    missing = expected - observed
    additional = observed - expected
    comparison_mode = spec[
        "comparison_mode"
    ]

    if comparison_mode == (
        "required_subset"
    ):
        verification = (
            "PASS"
            if not missing
            else "FAIL"
        )
    elif comparison_mode == "exact":
        verification = (
            "PASS"
            if (
                not missing
                and (
                    not additional
                    or not STRICT_NO_UNEXPECTED_FILE_DERIVED_FROM
                )
            )
            else "FAIL"
        )
    else:
        raise RuntimeError(
            "Unknown comparison mode during final reconciliation: "
            f"{comparison_mode}"
        )

    final_rows.append(
        {
            "operation_key": (
                spec["operation_key"]
            ),
            "screen": spec["screen"],
            "profile": spec["profile"],
            "record_id": (
                spec["record_id"]
            ),
            "property": (
                spec["property"]
            ),
            "comparison_mode": (
                comparison_mode
            ),
            "expected_count": len(
                expected
            ),
            "observed_count": len(
                observed
            ),
            "missing_count": len(
                missing
            ),
            "additional_count": len(
                additional
            ),
            "missing_values": ",".join(
                sorted(missing)
            ),
            "additional_values": ",".join(
                sorted(additional)
            ),
            "verification": (
                verification
            ),
        }
    )

final_reconciliation = pd.DataFrame(
    final_rows
)
FINAL_RECONCILIATION_PATH = (
    OUTPUT_ROOT
    / "analysis_lineage_reconciliation.csv"
)
final_reconciliation.to_csv(
    FINAL_RECONCILIATION_PATH,
    index=False,
    lineterminator="\n",
)

print(
    final_reconciliation[
        [
            "operation_key",
            "record_id",
            "property",
            "comparison_mode",
            "expected_count",
            "observed_count",
            "missing_count",
            "additional_count",
            "verification",
        ]
    ].to_string(index=False)
)
print(
    "\nSaved:",
    FINAL_RECONCILIATION_PATH,
)

if final_reconciliation[
    "verification"
].eq("PASS").all():
    marker = (
        OUTPUT_ROOT
        / "ANALYSIS_LINEAGE_PATCHED_6_OF_6.txt"
    )
    marker.write_text(
        "All six CRISPRi MACS lineage requirements "
        "are satisfied. Additional AnalysisSet input "
        "file sets, when present, were preserved.\n"
        f"Results: {FINAL_RECONCILIATION_PATH}\n",
        encoding="utf-8",
    )
    print(
        "Completion marker:",
        marker,
    )
else:
    raise RuntimeError(
        "Final lineage reconciliation failed"
    )
```

### Text output (stream)

```text
                                 operation_key      record_id        property comparison_mode  expected_count  observed_count  missing_count  additional_count verification
      genome_wide_analysis_set_input_file_sets IGVFDS9629VXQK input_file_sets required_subset               3               4              0                 1         PASS
       validation_analysis_set_input_file_sets IGVFDS4708PCVX input_file_sets required_subset               6               7              0                 1         PASS
genome_wide_guide_quantifications_derived_from IGVFFI8424MVZE    derived_from           exact             362             362              0                 0         PASS
 validation_guide_quantifications_derived_from IGVFFI5157STFU    derived_from           exact              41              41              0                 0         PASS
         genome_wide_differential_derived_from IGVFFI9001HIYI    derived_from           exact               1               1              0                 0         PASS
          validation_differential_derived_from IGVFFI3975KRIX    derived_from           exact               1               1              0                 0         PASS

Saved: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/analysis_lineage_reconciliation.csv
Completion marker: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/ANALYSIS_LINEAGE_PATCHED_6_OF_6.txt

```

## Markdown cell 24

## 11. Output summary

## Code cell 25

```python
output_candidates = [
    (
        "lineage input audit",
        LINEAGE_INPUT_AUDIT_PATH,
    ),
    (
        "expected SequenceFiles",
        EXPECTED_SEQUENCE_FILES_PATH,
    ),
    (
        "expected operations",
        EXPECTED_LINEAGE_PATH,
    ),
    (
        "upstream file-set audit",
        UPSTREAM_FILE_SET_AUDIT_PATH,
    ),
    (
        "additional AnalysisSet inputs audit",
        WORK_DIR
        / "additional_analysis_input_file_sets_audit.tsv",
    ),
    (
        "patch plan",
        LINEAGE_PATCH_PLAN_PATH,
    ),
    (
        "patch operations",
        LINEAGE_PATCH_OPERATIONS_PATH,
    ),
    (
        "dry-run results",
        WORK_DIR
        / "lineage_dry_run_results.tsv",
    ),
    (
        "production results",
        WORK_DIR
        / "lineage_production_results.tsv",
    ),
    (
        "final reconciliation",
        OUTPUT_ROOT
        / "analysis_lineage_reconciliation.csv",
    ),
]

output_summary = pd.DataFrame(
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
        for label, path in (
            output_candidates
        )
    ]
)

print(
    output_summary.to_string(
        index=False
    )
)
print("\nLogs:", LOG_DIR)
```

### Text output (stream)

```text
                             output                                                                                                                            path  exists  size_bytes
                lineage input audit                       /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/lineage_input_audit.tsv    True         440
             expected SequenceFiles                   /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/expected_sequence_files.tsv    True       18953
                expected operations               /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/expected_lineage_operations.tsv    True        1215
            upstream file-set audit                   /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/upstream_file_set_audit.tsv    True        1273
additional AnalysisSet inputs audit /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/additional_analysis_input_file_sets_audit.tsv    True         496
                         patch plan                        /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/lineage_patch_plan.tsv    True        7121
                   patch operations                  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/lineage_patch_operations.tsv    True        1651
                    dry-run results                   /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/lineage_dry_run_results.tsv    True        2332
                 production results                /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/_work/lineage_production_results.tsv    True        2332
               final reconciliation                 /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/analysis_lineage_reconciliation.csv    True         921

Logs: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/logs

```

## Code cell 26

```python

```

## Markdown cell 27

## Patching stale CuratedSet link

## Code cell 28

```python
# Exact AnalysisSet input_file_sets repair — self-contained cell.
#
# First run:
#     APPLY_PATCH = False
#
# After reviewing the plan, rerun with:
#     APPLY_PATCH = True
#     CONFIRMATION = "OVERWRITE_2_CRISPRI_MACS_ANALYSIS_SET_INPUT_FILE_SETS"

from __future__ import annotations

import json
import os
import re
from getpass import getpass
from pathlib import Path
from typing import Any

import pandas as pd
import requests


# =====================================================================
# 1. Safety configuration
# =====================================================================

#APPLY_PATCH = False

#CONFIRMATION = ""

APPLY_PATCH = True
CONFIRMATION = "OVERWRITE_2_CRISPRI_MACS_ANALYSIS_SET_INPUT_FILE_SETS"

REQUIRED_CONFIRMATION = (
    "OVERWRITE_2_CRISPRI_MACS_ANALYSIS_SET_INPUT_FILE_SETS"
)

API_BASE = "https://api.data.igvf.org"

REQUEST_TIMEOUT_SECONDS = 120

OUTPUT_ROOT = Path(
    "/oak/stanford/groups/engreitz/Users/tri/tmp/"
    "igvf_crispri_macs_lineage_patch/"
    "exact_analysis_set_fix"
)

OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

SEQUENCE_RECONCILIATION_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/tri/tmp/"
    "igvf_crispri_macs_seqspec/"
    "sequencefile_seqspec_reconciliation.csv"
)


# These arrays are the exact desired final values.
TARGETS = {
    "IGVFDS9629VXQK": {
        "screen": "genome_wide",
        "expected_input_file_sets": [
            "IGVFDS2967ULMS",
            "IGVFDS8680ANMB",
            "IGVFDS8690AKQF",
        ],
        "known_removable_extras": [
            "IGVFDS4419BVVU",
        ],
    },
    "IGVFDS4708PCVX": {
        "screen": "validation",
        "expected_input_file_sets": [
            "IGVFDS1344YZHW",
            "IGVFDS1797KANA",
            "IGVFDS1832ZFCK",
            "IGVFDS5537WUFJ",
            "IGVFDS7114EJDN",
            "IGVFDS7624AWFB",
        ],
        "known_removable_extras": [
            "IGVFDS8881NOFE",
        ],
    },
}


# These are verified before either AnalysisSet is modified.
OUTPUT_FILES = {
    "genome_wide": {
        "guide_quantifications": "IGVFFI8424MVZE",
        "differential_quantifications": "IGVFFI9001HIYI",
    },
    "validation": {
        "guide_quantifications": "IGVFFI5157STFU",
        "differential_quantifications": "IGVFFI3975KRIX",
    },
}


# =====================================================================
# 2. Credentials and Portal HTTP helpers
# =====================================================================

def _credential(
    name: str,
    prompt: str,
) -> str:
    value = os.environ.get(
        name,
        "",
    ).strip()

    if not value:
        value = getpass(prompt).strip()

    if not value:
        raise RuntimeError(
            f"{name} is required"
        )

    if (
        len(value) >= 2
        and value[0] in {"'", '"'}
        and value[-1] == value[0]
    ):
        raise RuntimeError(
            f"{name} has surrounding quotes. "
            "Paste the raw credential value."
        )

    os.environ[name] = value

    return value


API_KEY = <REDACTED>
    "IGVF_API_KEY",
    "IGVF API key: ",
)

SECRET_KEY = _credential(
    "IGVF_SECRET_KEY",
    "IGVF secret key: ",
)


session = requests.Session()

session.auth = (
    API_KEY,
    SECRET_KEY,
)

session.headers.update(
    {
        "Accept": "application/json",
        "Content-Type": "application/json",
    }
)


def _record_url(
    identifier: str,
) -> str:
    """
    Convert an accession or canonical @id into an API URL.
    """

    token = str(identifier).strip()

    if token.startswith("/"):
        return (
            API_BASE.rstrip("/")
            + token
        )

    return (
        API_BASE.rstrip("/")
        + "/"
        + token.strip("/")
        + "/"
    )


def get_record(
    identifier: str,
    *,
    frame: str = "edit",
    datastore: str | None = "database",
) -> dict:
    """
    Retrieve a record directly from the database.
    """

    params = {
        "format": "json",
        "frame": frame,
    }

    if datastore:
        params["datastore"] = datastore

    response = session.get(
        _record_url(identifier),
        params=params,
        timeout=REQUEST_TIMEOUT_SECONDS,
    )

    if not response.ok:
        raise RuntimeError(
            f"GET failed for {identifier}: "
            f"HTTP {response.status_code}\n"
            f"{response.text[:2000]}"
        )

    payload = response.json()

    if not isinstance(
        payload,
        dict,
    ):
        raise RuntimeError(
            f"Unexpected GET payload for "
            f"{identifier}: "
            f"{type(payload).__name__}"
        )

    return payload


def patch_record(
    record: dict,
    payload: dict,
) -> dict:
    """
    PATCH a record directly.

    Because this sends the complete input_file_sets array directly,
    the Portal receives an exact replacement value rather than an
    array extended by igvf_utils.
    """

    record_id = str(
        record.get("@id")
        or ""
    ).strip()

    if not record_id:
        raise RuntimeError(
            "Record has no @id. "
            "Refusing to PATCH."
        )

    response = session.patch(
        _record_url(record_id),
        json=payload,
        timeout=REQUEST_TIMEOUT_SECONDS,
    )

    if not response.ok:
        raise RuntimeError(
            f"PATCH failed for {record_id}: "
            f"HTTP {response.status_code}\n"
            f"{response.text[:4000]}"
        )

    response_payload = response.json()

    if (
        isinstance(
            response_payload,
            dict,
        )
        and isinstance(
            response_payload.get(
                "@graph"
            ),
            list,
        )
        and response_payload["@graph"]
    ):
        return response_payload[
            "@graph"
        ][0]

    if isinstance(
        response_payload,
        dict,
    ):
        return response_payload

    raise RuntimeError(
        f"Unexpected PATCH response for "
        f"{record_id}: "
        f"{type(response_payload).__name__}"
    )


def linked_accessions(
    value: Any,
    prefix: str,
) -> set[str]:
    """
    Extract accessions from strings, paths, embedded objects, or lists.
    """

    if value is None:
        items: list[Any] = []

    elif isinstance(
        value,
        (
            list,
            tuple,
            set,
        ),
    ):
        items = list(value)

    else:
        items = [value]

    accessions: set[str] = set()

    for item in items:

        candidates: list[str] = []

        if isinstance(
            item,
            dict,
        ):
            for key in (
                "accession",
                "@id",
                "uuid",
            ):
                candidate = item.get(
                    key
                )

                if candidate:
                    candidates.append(
                        str(candidate)
                    )

            aliases = item.get(
                "aliases",
                [],
            )

            if isinstance(
                aliases,
                str,
            ):
                aliases = [aliases]

            candidates.extend(
                map(
                    str,
                    aliases,
                )
            )

        else:
            candidates.append(
                str(item)
            )

        for candidate in candidates:

            token = (
                candidate
                .strip()
                .strip("/")
                .rsplit("/", 1)[-1]
            )

            if token.startswith(
                prefix
            ):
                accessions.add(
                    token
                )

    return accessions


def record_types(
    record: dict,
) -> set[str]:

    value = record.get(
        "@type",
        [],
    )

    if isinstance(
        value,
        str,
    ):
        return {value}

    if isinstance(
        value,
        list,
    ):
        return {
            str(item)
            for item in value
        }

    return set()


def canonical_at_id(
    record: dict,
    expected_accession: str,
) -> str:
    """
    Return a validated canonical Portal @id.
    """

    at_id = str(
        record.get("@id")
        or ""
    ).strip()

    accession = str(
        record.get("accession")
        or ""
    ).strip()

    if (
        accession
        and accession
        != expected_accession
    ):
        raise RuntimeError(
            "Resolved record accession "
            "does not match the requested "
            f"accession: expected "
            f"{expected_accession}, "
            f"observed {accession}"
        )

    if not at_id:
        raise RuntimeError(
            f"Resolved record "
            f"{expected_accession} "
            "has no @id"
        )

    return at_id


# =====================================================================
# 3. Verify that the file-level derived_from lineage is already exact
# =====================================================================

if not SEQUENCE_RECONCILIATION_PATH.is_file():
    raise FileNotFoundError(
        "The completed seqspec reconciliation "
        "is required:\n"
        f"{SEQUENCE_RECONCILIATION_PATH}"
    )


reconciliation = pd.read_csv(
    SEQUENCE_RECONCILIATION_PATH,
    dtype=str,
    keep_default_na=False,
)


required_columns = {
    "screen",
    "sequence_file_accession",
    "verification",
}

missing_columns = (
    required_columns
    - set(
        reconciliation.columns
    )
)

if missing_columns:
    raise RuntimeError(
        "Seqspec reconciliation is "
        "missing columns: "
        f"{sorted(missing_columns)}"
    )


if not reconciliation[
    "verification"
].eq("PASS").all():

    failed = reconciliation.loc[
        ~reconciliation[
            "verification"
        ].eq("PASS")
    ]

    raise RuntimeError(
        "Not every SequenceFile "
        "reconciliation row passed:\n"
        + failed.to_string(
            index=False
        )
    )


if reconciliation[
    "sequence_file_accession"
].duplicated().any():

    duplicates = reconciliation.loc[
        reconciliation[
            "sequence_file_accession"
        ].duplicated(
            keep=False
        )
    ]

    raise RuntimeError(
        "Duplicate SequenceFile accessions "
        "in reconciliation:\n"
        + duplicates.to_string(
            index=False
        )
    )


sequence_files_by_screen = {
    screen: set(
        group[
            "sequence_file_accession"
        ].astype(str)
    )
    for screen, group in (
        reconciliation.groupby(
            "screen"
        )
    )
}


expected_counts = {
    "genome_wide": 362,
    "validation": 41,
}


for (
    screen,
    expected_count,
) in expected_counts.items():

    observed_count = len(
        sequence_files_by_screen.get(
            screen,
            set(),
        )
    )

    if observed_count != expected_count:
        raise RuntimeError(
            f"{screen}: expected "
            f"{expected_count} reconciled "
            "SequenceFiles, observed "
            f"{observed_count}"
        )


derived_audit_rows: list[
    dict
] = []


for (
    screen,
    files,
) in OUTPUT_FILES.items():

    guide_accession = files[
        "guide_quantifications"
    ]

    differential_accession = files[
        "differential_quantifications"
    ]

    guide_record = get_record(
        guide_accession
    )

    differential_record = get_record(
        differential_accession
    )

    guide_derived = linked_accessions(
        guide_record.get(
            "derived_from",
            [],
        ),
        "IGVFFI",
    )

    differential_derived = (
        linked_accessions(
            differential_record.get(
                "derived_from",
                [],
            ),
            "IGVFFI",
        )
    )

    expected_guide_derived = (
        sequence_files_by_screen[
            screen
        ]
    )

    expected_differential_derived = {
        guide_accession
    }

    derived_audit_rows.extend(
        [
            {
                "screen": screen,
                "record_id": (
                    guide_accession
                ),
                "role": (
                    "guide_quantifications"
                ),
                "expected_count": len(
                    expected_guide_derived
                ),
                "observed_count": len(
                    guide_derived
                ),
                "missing_count": len(
                    expected_guide_derived
                    - guide_derived
                ),
                "additional_count": len(
                    guide_derived
                    - expected_guide_derived
                ),
                "status": (
                    "PASS"
                    if guide_derived
                    == expected_guide_derived
                    else "FAIL"
                ),
            },
            {
                "screen": screen,
                "record_id": (
                    differential_accession
                ),
                "role": (
                    "differential_quantifications"
                ),
                "expected_count": 1,
                "observed_count": len(
                    differential_derived
                ),
                "missing_count": len(
                    expected_differential_derived
                    - differential_derived
                ),
                "additional_count": len(
                    differential_derived
                    - expected_differential_derived
                ),
                "status": (
                    "PASS"
                    if differential_derived
                    == expected_differential_derived
                    else "FAIL"
                ),
            },
        ]
    )


derived_audit = pd.DataFrame(
    derived_audit_rows
)

derived_audit_path = (
    OUTPUT_ROOT
    / "derived_from_preflight.tsv"
)

derived_audit.to_csv(
    derived_audit_path,
    sep="\t",
    index=False,
    lineterminator="\n",
)


print(
    "Derived-from preflight:"
)

print(
    derived_audit.to_string(
        index=False
    )
)


if not derived_audit[
    "status"
].eq("PASS").all():

    raise RuntimeError(
        "The four output files do not "
        "currently have the expected "
        "derived_from lineage. Review "
        f"{derived_audit_path}"
    )


# =====================================================================
# 4. Resolve expected MeasurementSets and build an exact replacement plan
# =====================================================================

plan_rows: list[dict] = []

records_before: dict[
    str,
    dict,
] = {}

patch_payloads: dict[
    str,
    dict,
] = {}


for (
    analysis_set_accession,
    spec,
) in TARGETS.items():

    analysis_record = get_record(
        analysis_set_accession,
        frame="edit",
    )

    analysis_object = get_record(
        analysis_set_accession,
        frame="object",
    )

    records_before[
        analysis_set_accession
    ] = analysis_record

    if (
        "AnalysisSet"
        not in record_types(
            analysis_object
        )
    ):
        raise RuntimeError(
            f"{analysis_set_accession} "
            "is not an AnalysisSet: "
            f"{sorted(record_types(analysis_object))}"
        )

    expected_accessions = set(
        spec[
            "expected_input_file_sets"
        ]
    )

    known_extras = set(
        spec[
            "known_removable_extras"
        ]
    )

    existing_accessions = (
        linked_accessions(
            analysis_record.get(
                "input_file_sets",
                [],
            ),
            "IGVFDS",
        )
    )

    unknown_extras = (
        existing_accessions
        - expected_accessions
        - known_extras
    )

    if unknown_extras:
        raise RuntimeError(
            f"{analysis_set_accession} "
            "contains unexpected "
            "input_file_sets outside the "
            "approved removal plan: "
            f"{sorted(unknown_extras)}"
        )

    canonical_expected_ids: list[
        str
    ] = []

    for upstream_accession in sorted(
        expected_accessions
    ):

        upstream_record = get_record(
            upstream_accession,
            frame="object",
        )

        upstream_types = record_types(
            upstream_record
        )

        if (
            "MeasurementSet"
            not in upstream_types
        ):
            raise RuntimeError(
                f"{upstream_accession} "
                "is not a MeasurementSet: "
                f"{sorted(upstream_types)}"
            )

        canonical_expected_ids.append(
            canonical_at_id(
                upstream_record,
                upstream_accession,
            )
        )

    payload = {
        "input_file_sets": (
            canonical_expected_ids
        ),
    }

    patch_payloads[
        analysis_set_accession
    ] = payload

    plan_rows.append(
        {
            "screen": spec["screen"],
            "analysis_set": (
                analysis_set_accession
            ),
            "existing_count": len(
                existing_accessions
            ),
            "expected_count": len(
                expected_accessions
            ),
            "missing_expected": ",".join(
                sorted(
                    expected_accessions
                    - existing_accessions
                )
            ),
            "known_extras_to_remove": (
                ",".join(
                    sorted(
                        existing_accessions
                        & known_extras
                    )
                )
            ),
            "unknown_extras": ",".join(
                sorted(
                    unknown_extras
                )
            ),
            "resulting_input_file_sets": (
                ",".join(
                    sorted(
                        expected_accessions
                    )
                )
            ),
            "action": (
                "SKIP_ALREADY_EXACT"
                if existing_accessions
                == expected_accessions
                else "OVERWRITE_EXACT_ARRAY"
            ),
        }
    )


plan = pd.DataFrame(
    plan_rows
)

plan_path = (
    OUTPUT_ROOT
    / "exact_analysis_set_patch_plan.tsv"
)

plan.to_csv(
    plan_path,
    sep="\t",
    index=False,
    lineterminator="\n",
)


before_path = (
    OUTPUT_ROOT
    / "analysis_sets_before_patch.json"
)

before_path.write_text(
    json.dumps(
        records_before,
        indent=2,
        ensure_ascii=False,
    )
    + "\n",
    encoding="utf-8",
)


payload_path = (
    OUTPUT_ROOT
    / "exact_patch_payloads.json"
)

payload_path.write_text(
    json.dumps(
        patch_payloads,
        indent=2,
        ensure_ascii=False,
    )
    + "\n",
    encoding="utf-8",
)


print(
    "\nExact AnalysisSet "
    "replacement plan:"
)

print(
    plan.to_string(
        index=False
    )
)

print(
    "\nSaved plan:",
    plan_path,
)

print(
    "Saved pre-patch backup:",
    before_path,
)

print(
    "Saved exact payloads:",
    payload_path,
)


# =====================================================================
# 5. Production gate and direct Portal PATCH
# =====================================================================

if not APPLY_PATCH:

    print(
        "\nDRY RUN ONLY — no Portal "
        "records were changed.\n\n"
        "Review the plan above. To "
        "apply the exact array "
        "replacement, change the top "
        "of this cell to:\n\n"
        "APPLY_PATCH = True\n"
        f'CONFIRMATION = '
        f'"{REQUIRED_CONFIRMATION}"\n\n'
        "Then rerun this same cell."
    )

else:

    if (
        CONFIRMATION
        != REQUIRED_CONFIRMATION
    ):
        raise RuntimeError(
            "Production confirmation is "
            "missing or incorrect.\n"
            f'Required:\nCONFIRMATION = '
            f'"{REQUIRED_CONFIRMATION}"'
        )

    patch_results: list[
        dict
    ] = []

    for (
        analysis_set_accession,
        spec,
    ) in TARGETS.items():

        expected_accessions = set(
            spec[
                "expected_input_file_sets"
            ]
        )

        current_record = get_record(
            analysis_set_accession
        )

        current_accessions = (
            linked_accessions(
                current_record.get(
                    "input_file_sets",
                    [],
                ),
                "IGVFDS",
            )
        )

        if (
            current_accessions
            == expected_accessions
        ):

            patch_results.append(
                {
                    "analysis_set": (
                        analysis_set_accession
                    ),
                    "action": (
                        "SKIP_ALREADY_EXACT"
                    ),
                    "http_patch_performed": (
                        False
                    ),
                    "verification": "PASS",
                }
            )

            continue

        allowed_current = (
            expected_accessions
            | set(
                spec[
                    "known_removable_extras"
                ]
            )
        )

        if not (
            current_accessions
            <= allowed_current
        ):
            raise RuntimeError(
                f"{analysis_set_accession} "
                "changed after planning and "
                "now contains unapproved "
                "values: "
                f"{sorted(current_accessions - allowed_current)}"
            )

        patched_response = patch_record(
            current_record,
            patch_payloads[
                analysis_set_accession
            ],
        )

        # Verify against the live database,
        # not a possibly stale search index.
        verified_record = get_record(
            analysis_set_accession
        )

        verified_accessions = (
            linked_accessions(
                verified_record.get(
                    "input_file_sets",
                    [],
                ),
                "IGVFDS",
            )
        )

        if (
            verified_accessions
            != expected_accessions
        ):
            raise RuntimeError(
                f"{analysis_set_accession} "
                "verification failed.\n"
                f"Expected: "
                f"{sorted(expected_accessions)}\n"
                f"Observed: "
                f"{sorted(verified_accessions)}"
            )

        response_path = (
            OUTPUT_ROOT
            / (
                f"{analysis_set_accession}"
                "_patch_response.json"
            )
        )

        response_path.write_text(
            json.dumps(
                patched_response,
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )

        patch_results.append(
            {
                "analysis_set": (
                    analysis_set_accession
                ),
                "action": (
                    "OVERWRITE_EXACT_ARRAY"
                ),
                "http_patch_performed": (
                    True
                ),
                "observed_count_after_patch": (
                    len(
                        verified_accessions
                    )
                ),
                "verification": "PASS",
                "response_json": str(
                    response_path
                ),
            }
        )

    patch_results_frame = (
        pd.DataFrame(
            patch_results
        )
    )

    results_path = (
        OUTPUT_ROOT
        / "exact_analysis_set_patch_results.tsv"
    )

    patch_results_frame.to_csv(
        results_path,
        sep="\t",
        index=False,
        lineterminator="\n",
    )

    completion_marker = (
        OUTPUT_ROOT
        / (
            "ANALYSIS_SET_INPUT_FILE_SETS_"
            "EXACT_2_OF_2.txt"
        )
    )

    completion_marker.write_text(
        "Both CRISPRi MACS "
        "AnalysisSets have the exact "
        "expected MeasurementSet "
        "input_file_sets arrays.\n"
        f"Results: {results_path}\n",
        encoding="utf-8",
    )

    print(
        "\nProduction patch results:"
    )

    print(
        patch_results_frame.to_string(
            index=False
        )
    )

    print(
        "\nResults:",
        results_path,
    )

    print(
        "Completion marker:",
        completion_marker,
    )

    print(
        "\nThe database state is now "
        "exact. Portal audit badges may "
        "take several minutes to reindex "
        "and disappear."
    )
```

### Text output (stream)

```text
Derived-from preflight:
     screen      record_id                         role  expected_count  observed_count  missing_count  additional_count status
genome_wide IGVFFI8424MVZE        guide_quantifications             362             362              0                 0   PASS
genome_wide IGVFFI9001HIYI differential_quantifications               1               1              0                 0   PASS
 validation IGVFFI5157STFU        guide_quantifications              41              41              0                 0   PASS
 validation IGVFFI3975KRIX differential_quantifications               1               1              0                 0   PASS

Exact AnalysisSet replacement plan:
     screen   analysis_set  existing_count  expected_count missing_expected known_extras_to_remove unknown_extras                                                                 resulting_input_file_sets                action
genome_wide IGVFDS9629VXQK               4               3                          IGVFDS4419BVVU                                                             IGVFDS2967ULMS,IGVFDS8680ANMB,IGVFDS8690AKQF OVERWRITE_EXACT_ARRAY
 validation IGVFDS4708PCVX               7               6                          IGVFDS8881NOFE                IGVFDS1344YZHW,IGVFDS1797KANA,IGVFDS1832ZFCK,IGVFDS5537WUFJ,IGVFDS7114EJDN,IGVFDS7624AWFB OVERWRITE_EXACT_ARRAY

Saved plan: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/exact_analysis_set_fix/exact_analysis_set_patch_plan.tsv
Saved pre-patch backup: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/exact_analysis_set_fix/analysis_sets_before_patch.json
Saved exact payloads: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/exact_analysis_set_fix/exact_patch_payloads.json

```

### Text output (error)

```text
[31m---------------------------------------------------------------------------[39m
[31mRuntimeError[39m                              Traceback (most recent call last)
[36mCell[39m[36m [39m[32mIn[21][39m[32m, line 1094[39m
[32m   1082[39m [38;5;28;01mif[39;00m [38;5;129;01mnot[39;00m (
[32m   1083[39m     current_accessions
[32m   1084[39m     <= allowed_current
[32m   1085[39m ):
[32m   1086[39m     [38;5;28;01mraise[39;00m [38;5;167;01mRuntimeError[39;00m(
[32m   1087[39m         [33mf[39m[33m"[39m[38;5;132;01m{[39;00manalysis_set_accession[38;5;132;01m}[39;00m[33m [39m[33m"[39m
[32m   1088[39m         [33m"[39m[33mchanged after planning and [39m[33m"[39m
[32m   (...)[39m[32m   1091[39m         [33mf[39m[33m"[39m[38;5;132;01m{[39;00m[38;5;28msorted[39m(current_accessions[38;5;250m [39m-[38;5;250m [39mallowed_current)[38;5;132;01m}[39;00m[33m"[39m
[32m   1092[39m     )
[32m-> [39m[32m1094[39m patched_response = [43mpatch_record[49m[43m([49m
[32m   1095[39m [43m    [49m[43mcurrent_record[49m[43m,[49m
[32m   1096[39m [43m    [49m[43mpatch_payloads[49m[43m[[49m
[32m   1097[39m [43m        [49m[43manalysis_set_accession[49m
[32m   1098[39m [43m    [49m[43m][49m[43m,[49m
[32m   1099[39m [43m[49m[43m)[49m
[32m   1101[39m [38;5;66;03m# Verify against the live database,[39;00m
[32m   1102[39m [38;5;66;03m# not a possibly stale search index.[39;00m
[32m   1103[39m verified_record = get_record(
[32m   1104[39m     analysis_set_accession
[32m   1105[39m )

[36mCell[39m[36m [39m[32mIn[21][39m[32m, line 252[39m, in [36mpatch_record[39m[34m(record, payload)[39m
[32m    246[39m record_id = [38;5;28mstr[39m(
[32m    247[39m     record.get([33m"[39m[33m@id[39m[33m"[39m)
[32m    248[39m     [38;5;129;01mor[39;00m [33m"[39m[33m"[39m
[32m    249[39m ).strip()
[32m    251[39m [38;5;28;01mif[39;00m [38;5;129;01mnot[39;00m record_id:
[32m--> [39m[32m252[39m     [38;5;28;01mraise[39;00m [38;5;167;01mRuntimeError[39;00m(
[32m    253[39m         [33m"[39m[33mRecord has no @id. [39m[33m"[39m
[32m    254[39m         [33m"[39m[33mRefusing to PATCH.[39m[33m"[39m
[32m    255[39m     )
[32m    257[39m response = session.patch(
[32m    258[39m     _record_url(record_id),
[32m    259[39m     json=payload,
[32m    260[39m     timeout=REQUEST_TIMEOUT_SECONDS,
[32m    261[39m )
[32m    263[39m [38;5;28;01mif[39;00m [38;5;129;01mnot[39;00m response.ok:

[31mRuntimeError[39m: Record has no @id. Refusing to PATCH.
```

## Code cell 29

```python

# Production-only fix for the exact AnalysisSet overwrite.
#
# Paste this cell AFTER the successful dry-run cell. It uses the existing
# TARGETS, patch_payloads, session, get_record, linked_accessions, and OUTPUT_ROOT.
#
# It fixes the previous failure by PATCHing with the known AnalysisSet accession
# instead of requiring @id from a frame=edit response.

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


APPLY_PATCH = True
CONFIRMATION = (
    "OVERWRITE_2_CRISPRI_MACS_ANALYSIS_SET_INPUT_FILE_SETS"
)
REQUIRED_CONFIRMATION = (
    "OVERWRITE_2_CRISPRI_MACS_ANALYSIS_SET_INPUT_FILE_SETS"
)


# ---------------------------------------------------------------------
# 1. Confirm that the successful dry-run state is still present
# ---------------------------------------------------------------------

required_names = [
    "session",
    "_record_url",
    "get_record",
    "linked_accessions",
    "TARGETS",
    "patch_payloads",
    "OUTPUT_ROOT",
    "REQUEST_TIMEOUT_SECONDS",
]

missing_names = [
    name
    for name in required_names
    if name not in globals()
]

if missing_names:
    raise RuntimeError(
        "The current kernel is missing objects created by the dry-run cell: "
        f"{missing_names}. Rerun the dry-run cell first, then run this cell."
    )

if not APPLY_PATCH:
    raise RuntimeError(
        "APPLY_PATCH is False. Set it to True only after reviewing the dry run."
    )

if CONFIRMATION != REQUIRED_CONFIRMATION:
    raise RuntimeError(
        "Production confirmation is missing or incorrect.\n"
        f'Required:\nCONFIRMATION = "{REQUIRED_CONFIRMATION}"'
    )


# ---------------------------------------------------------------------
# 2. Correct PATCH helper
# ---------------------------------------------------------------------

def patch_record_by_identifier(
    identifier: str,
    payload: dict,
) -> dict:
    """
    PATCH by the known accession or canonical path.

    The prior helper incorrectly required @id from frame=edit, but that response
    can omit @id. The accession is already known and is a valid record URL.
    """
    response = session.patch(
        _record_url(identifier),
        json=payload,
        timeout=REQUEST_TIMEOUT_SECONDS,
    )

    if not response.ok:
        raise RuntimeError(
            f"PATCH failed for {identifier}: HTTP {response.status_code}\n"
            f"{response.text[:5000]}"
        )

    response_payload = response.json()

    if (
        isinstance(response_payload, dict)
        and isinstance(response_payload.get("@graph"), list)
        and response_payload["@graph"]
    ):
        return response_payload["@graph"][0]

    if isinstance(response_payload, dict):
        return response_payload

    raise RuntimeError(
        f"Unexpected PATCH response for {identifier}: "
        f"{type(response_payload).__name__}"
    )


# ---------------------------------------------------------------------
# 3. Recheck the current AnalysisSet values immediately before PATCH
# ---------------------------------------------------------------------

pre_patch_rows = []

for analysis_set_accession, spec in TARGETS.items():
    current_record = get_record(
        analysis_set_accession,
        frame="edit",
    )

    current_accessions = linked_accessions(
        current_record.get("input_file_sets", []),
        "IGVFDS",
    )

    expected_accessions = set(
        spec["expected_input_file_sets"]
    )

    known_extras = set(
        spec["known_removable_extras"]
    )

    unapproved_values = (
        current_accessions
        - expected_accessions
        - known_extras
    )

    if unapproved_values:
        raise RuntimeError(
            f"{analysis_set_accession} now contains unapproved values: "
            f"{sorted(unapproved_values)}"
        )

    pre_patch_rows.append(
        {
            "analysis_set": analysis_set_accession,
            "screen": spec["screen"],
            "current_count": len(current_accessions),
            "expected_count": len(expected_accessions),
            "current_values": ",".join(sorted(current_accessions)),
            "expected_values": ",".join(sorted(expected_accessions)),
            "known_extras_present": ",".join(
                sorted(current_accessions & known_extras)
            ),
            "action": (
                "SKIP_ALREADY_EXACT"
                if current_accessions == expected_accessions
                else "OVERWRITE_EXACT_ARRAY"
            ),
        }
    )

pre_patch_plan = pd.DataFrame(pre_patch_rows)

print("Immediate pre-patch plan:")
print(pre_patch_plan.to_string(index=False))


# ---------------------------------------------------------------------
# 4. Apply exact replacement and verify from the live database
# ---------------------------------------------------------------------

results = []

for analysis_set_accession, spec in TARGETS.items():
    expected_accessions = set(
        spec["expected_input_file_sets"]
    )

    current_record = get_record(
        analysis_set_accession,
        frame="edit",
    )

    current_accessions = linked_accessions(
        current_record.get("input_file_sets", []),
        "IGVFDS",
    )

    if current_accessions == expected_accessions:
        results.append(
            {
                "analysis_set": analysis_set_accession,
                "screen": spec["screen"],
                "action": "SKIP_ALREADY_EXACT",
                "http_patch_performed": False,
                "observed_count_after_patch": len(current_accessions),
                "verification": "PASS",
                "response_json": "",
            }
        )
        continue

    allowed_before = (
        expected_accessions
        | set(spec["known_removable_extras"])
    )

    if not current_accessions <= allowed_before:
        raise RuntimeError(
            f"{analysis_set_accession} changed after planning and contains "
            f"unapproved values: {sorted(current_accessions - allowed_before)}"
        )

    # Critical fix: use the accession directly.
    patched_response = patch_record_by_identifier(
        analysis_set_accession,
        patch_payloads[analysis_set_accession],
    )

    verified_record = get_record(
        analysis_set_accession,
        frame="edit",
    )

    verified_accessions = linked_accessions(
        verified_record.get("input_file_sets", []),
        "IGVFDS",
    )

    if verified_accessions != expected_accessions:
        raise RuntimeError(
            f"{analysis_set_accession} verification failed.\n"
            f"Expected: {sorted(expected_accessions)}\n"
            f"Observed: {sorted(verified_accessions)}"
        )

    response_path = (
        Path(OUTPUT_ROOT)
        / f"{analysis_set_accession}_patch_response.json"
    )

    response_path.write_text(
        json.dumps(
            patched_response,
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )

    results.append(
        {
            "analysis_set": analysis_set_accession,
            "screen": spec["screen"],
            "action": "OVERWRITE_EXACT_ARRAY",
            "http_patch_performed": True,
            "observed_count_after_patch": len(verified_accessions),
            "verification": "PASS",
            "response_json": str(response_path),
        }
    )


# ---------------------------------------------------------------------
# 5. Save completion evidence
# ---------------------------------------------------------------------

results_frame = pd.DataFrame(results)

results_path = (
    Path(OUTPUT_ROOT)
    / "exact_analysis_set_patch_results.tsv"
)

results_frame.to_csv(
    results_path,
    sep="\t",
    index=False,
    lineterminator="\n",
)

if not results_frame["verification"].eq("PASS").all():
    raise RuntimeError(
        f"One or more AnalysisSet verifications failed. Review {results_path}"
    )

completion_marker = (
    Path(OUTPUT_ROOT)
    / "ANALYSIS_SET_INPUT_FILE_SETS_EXACT_2_OF_2.txt"
)

completion_marker.write_text(
    "Both CRISPRi MACS AnalysisSets have the exact expected "
    "MeasurementSet input_file_sets arrays.\n"
    f"Results: {results_path}\n",
    encoding="utf-8",
)

print("\nProduction patch results:")
print(results_frame.to_string(index=False))
print("\nResults:", results_path)
print("Completion marker:", completion_marker)
print(
    "\nThe live database values were verified. Portal audit badges may "
    "take several minutes to reindex."
)

```

### Text output (stream)

```text
Immediate pre-patch plan:
  analysis_set      screen  current_count  expected_count                                                                                           current_values                                                                           expected_values known_extras_present                action
IGVFDS9629VXQK genome_wide              4               3                                              IGVFDS2967ULMS,IGVFDS4419BVVU,IGVFDS8680ANMB,IGVFDS8690AKQF                                              IGVFDS2967ULMS,IGVFDS8680ANMB,IGVFDS8690AKQF       IGVFDS4419BVVU OVERWRITE_EXACT_ARRAY
IGVFDS4708PCVX  validation              7               6 IGVFDS1344YZHW,IGVFDS1797KANA,IGVFDS1832ZFCK,IGVFDS5537WUFJ,IGVFDS7114EJDN,IGVFDS7624AWFB,IGVFDS8881NOFE IGVFDS1344YZHW,IGVFDS1797KANA,IGVFDS1832ZFCK,IGVFDS5537WUFJ,IGVFDS7114EJDN,IGVFDS7624AWFB       IGVFDS8881NOFE OVERWRITE_EXACT_ARRAY

Production patch results:
  analysis_set      screen                action  http_patch_performed  observed_count_after_patch verification                                                                                                                         response_json
IGVFDS9629VXQK genome_wide OVERWRITE_EXACT_ARRAY                  True                           3         PASS /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/exact_analysis_set_fix/IGVFDS9629VXQK_patch_response.json
IGVFDS4708PCVX  validation OVERWRITE_EXACT_ARRAY                  True                           6         PASS /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/exact_analysis_set_fix/IGVFDS4708PCVX_patch_response.json

Results: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/exact_analysis_set_fix/exact_analysis_set_patch_results.tsv
Completion marker: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/exact_analysis_set_fix/ANALYSIS_SET_INPUT_FILE_SETS_EXACT_2_OF_2.txt

The live database values were verified. Portal audit badges may take several minutes to reindex.

```

## Markdown cell 30

## Patching Curated Set

## Code cell 31

```python

# Add CuratedSet guide-list provenance to both CRISPRi MACS analyses.
#
# FIRST RUN:
#   APPLY_PATCH = False
#
# Review the plan. Then change only:
#   APPLY_PATCH = True
#   CONFIRMATION = "PATCH_2_CURATED_GUIDE_INPUTS_AND_ANALYSIS_SETS"
#
# and rerun this same cell.
#
# Production order is deliberate:
#   1. Add each guide-list file to guide quantification `derived_from`.
#   2. Add its CuratedSet to AnalysisSet `input_file_sets`.
#
# Differential-result lineage is verified but not modified.

from __future__ import annotations

import json
import os
import re
from getpass import getpass
from pathlib import Path
from typing import Any

import pandas as pd
import requests


# =====================================================================
# 1. Safety configuration
# =====================================================================


APPLY_PATCH = True

CONFIRMATION = (
    "PATCH_2_CURATED_GUIDE_INPUTS_AND_ANALYSIS_SETS"
)

REQUIRED_CONFIRMATION = (
    "PATCH_2_CURATED_GUIDE_INPUTS_AND_ANALYSIS_SETS"
)

API_BASE = "https://api.data.igvf.org"
REQUEST_TIMEOUT_SECONDS = 180

OUTPUT_ROOT = Path(
    "/oak/stanford/groups/engreitz/Users/tri/tmp/"
    "igvf_crispri_macs_lineage_patch/"
    "curated_guide_input_fix"
)
OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

SEQUENCE_RECONCILIATION_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/tri/tmp/"
    "igvf_crispri_macs_seqspec/"
    "sequencefile_seqspec_reconciliation.csv"
)

CONFIG = {
    "genome_wide": {
        "analysis_set": "IGVFDS9629VXQK",
        "measurement_sets": [
            "IGVFDS2967ULMS",
            "IGVFDS8680ANMB",
            "IGVFDS8690AKQF",
        ],
        "curated_set": "IGVFDS4419BVVU",
        "guide_list_file": "IGVFFI9629WVUU",
        "guide_quantifications": "IGVFFI8424MVZE",
        "differential_quantifications": "IGVFFI9001HIYI",
        "expected_fastq_count": 362,
    },
    "validation": {
        "analysis_set": "IGVFDS4708PCVX",
        "measurement_sets": [
            "IGVFDS1344YZHW",
            "IGVFDS1797KANA",
            "IGVFDS1832ZFCK",
            "IGVFDS5537WUFJ",
            "IGVFDS7114EJDN",
            "IGVFDS7624AWFB",
        ],
        "curated_set": "IGVFDS8881NOFE",
        "guide_list_file": "IGVFFI6745HXWN",
        "guide_quantifications": "IGVFFI5157STFU",
        "differential_quantifications": "IGVFFI3975KRIX",
        "expected_fastq_count": 41,
    },
}


# =====================================================================
# 2. Credentials and direct Portal HTTP helpers
# =====================================================================

def _credential(
    name: str,
    prompt: str,
) -> str:
    value = os.environ.get(
        name,
        "",
    ).strip()

    if not value:
        value = getpass(prompt).strip()

    if not value:
        raise RuntimeError(
            f"{name} is required"
        )

    if (
        len(value) >= 2
        and value[0] in {"'", '"'}
        and value[-1] == value[0]
    ):
        raise RuntimeError(
            f"{name} has surrounding quotes. "
            "Paste the raw credential value."
        )

    os.environ[name] = value
    return value


API_KEY = <REDACTED>
    "IGVF_API_KEY",
    "IGVF API key: ",
)
SECRET_KEY = _credential(
    "IGVF_SECRET_KEY",
    "IGVF secret key: ",
)

session = requests.Session()
session.auth = (
    API_KEY,
    SECRET_KEY,
)
session.headers.update(
    {
        "Accept": "application/json",
        "Content-Type": "application/json",
    }
)


def _record_url(
    identifier: str,
) -> str:
    token = str(identifier).strip()

    if token.startswith("/"):
        return (
            API_BASE.rstrip("/")
            + token
        )

    return (
        API_BASE.rstrip("/")
        + "/"
        + token.strip("/")
        + "/"
    )


def get_record(
    identifier: str,
    *,
    frame: str = "edit",
    datastore: str | None = "database",
) -> dict:
    params = {
        "format": "json",
        "frame": frame,
    }

    if datastore:
        params["datastore"] = datastore

    response = session.get(
        _record_url(identifier),
        params=params,
        timeout=REQUEST_TIMEOUT_SECONDS,
    )

    if not response.ok:
        raise RuntimeError(
            f"GET failed for {identifier}: "
            f"HTTP {response.status_code}\n"
            f"{response.text[:3000]}"
        )

    payload = response.json()

    if not isinstance(
        payload,
        dict,
    ):
        raise RuntimeError(
            f"Unexpected GET payload for "
            f"{identifier}: "
            f"{type(payload).__name__}"
        )

    return payload


def patch_record(
    identifier: str,
    payload: dict,
) -> dict:
    response = session.patch(
        _record_url(identifier),
        json=payload,
        timeout=REQUEST_TIMEOUT_SECONDS,
    )

    if not response.ok:
        raise RuntimeError(
            f"PATCH failed for {identifier}: "
            f"HTTP {response.status_code}\n"
            f"{response.text[:5000]}"
        )

    response_payload = response.json()

    if (
        isinstance(
            response_payload,
            dict,
        )
        and isinstance(
            response_payload.get("@graph"),
            list,
        )
        and response_payload["@graph"]
    ):
        return response_payload[
            "@graph"
        ][0]

    if isinstance(
        response_payload,
        dict,
    ):
        return response_payload

    raise RuntimeError(
        f"Unexpected PATCH response for "
        f"{identifier}: "
        f"{type(response_payload).__name__}"
    )


def record_types(
    record: dict,
) -> set[str]:
    value = record.get(
        "@type",
        [],
    )

    if isinstance(
        value,
        str,
    ):
        return {value}

    if isinstance(
        value,
        list,
    ):
        return {
            str(item)
            for item in value
        }

    return set()


def extract_accession(
    value: Any,
    prefix_pattern: str = r"IGVF(?:DS|FI)",
) -> str:
    """
    Extract one IGVF accession from a path, accession, UUID-embedded object,
    or linked object.
    """
    candidates: list[str] = []

    if isinstance(
        value,
        dict,
    ):
        for key in (
            "accession",
            "@id",
            "uuid",
        ):
            if value.get(key):
                candidates.append(
                    str(value[key])
                )
    elif value is not None:
        candidates.append(
            str(value)
        )

    pattern = re.compile(
        rf"\b({prefix_pattern}[A-Z0-9]+)\b"
    )

    for candidate in candidates:
        match = pattern.search(
            candidate
        )
        if match:
            return match.group(1)

    return ""


def link_items(
    value: Any,
) -> list[Any]:
    if value is None:
        return []

    if isinstance(
        value,
        (
            list,
            tuple,
            set,
        ),
    ):
        return list(value)

    return [value]


def link_accessions(
    value: Any,
    *,
    prefix: str,
) -> set[str]:
    accessions: set[str] = set()

    for item in link_items(value):
        accession = extract_accession(
            item,
            re.escape(prefix),
        )
        if accession:
            accessions.add(
                accession
            )

    return accessions


def link_map(
    value: Any,
    *,
    prefix: str,
) -> dict[str, str]:
    """
    Map accession -> existing link representation.

    The existing canonical paths are preserved when building exact replacement
    arrays, avoiding hundreds of unnecessary per-file GET requests.
    """
    output: dict[str, str] = {}

    for item in link_items(value):
        accession = extract_accession(
            item,
            re.escape(prefix),
        )
        if not accession:
            continue

        if isinstance(
            item,
            dict,
        ):
            canonical = str(
                item.get("@id")
                or item.get("accession")
                or accession
            )
        else:
            canonical = str(item)

        output[
            accession
        ] = canonical

    return output


def canonical_at_id(
    record: dict,
    expected_accession: str,
) -> str:
    accession = str(
        record.get("accession")
        or ""
    ).strip()

    if (
        accession
        and accession
        != expected_accession
    ):
        raise RuntimeError(
            f"Resolved accession mismatch: "
            f"expected {expected_accession}, "
            f"observed {accession}"
        )

    at_id = str(
        record.get("@id")
        or ""
    ).strip()

    if not at_id:
        raise RuntimeError(
            f"{expected_accession} has no @id "
            "in its object-frame record"
        )

    return at_id


# =====================================================================
# 3. Load the reconciled FASTQ accession sets
# =====================================================================

if not SEQUENCE_RECONCILIATION_PATH.is_file():
    raise FileNotFoundError(
        SEQUENCE_RECONCILIATION_PATH
    )

reconciliation = pd.read_csv(
    SEQUENCE_RECONCILIATION_PATH,
    dtype=str,
    keep_default_na=False,
)

required_columns = {
    "screen",
    "sequence_file_accession",
    "verification",
}

missing_columns = (
    required_columns
    - set(
        reconciliation.columns
    )
)

if missing_columns:
    raise RuntimeError(
        "Seqspec reconciliation is "
        "missing required columns: "
        f"{sorted(missing_columns)}"
    )

if not reconciliation[
    "verification"
].eq("PASS").all():
    bad = reconciliation.loc[
        ~reconciliation[
            "verification"
        ].eq("PASS")
    ]

    raise RuntimeError(
        "Not every SequenceFile "
        "reconciliation row passed:\n"
        + bad.to_string(
            index=False
        )
    )

if reconciliation[
    "sequence_file_accession"
].duplicated().any():
    duplicates = reconciliation.loc[
        reconciliation[
            "sequence_file_accession"
        ].duplicated(
            keep=False
        )
    ]

    raise RuntimeError(
        "Duplicate SequenceFile accessions "
        "in reconciliation:\n"
        + duplicates.to_string(
            index=False
        )
    )

fastqs_by_screen = {
    screen: set(
        group[
            "sequence_file_accession"
        ].astype(str)
    )
    for screen, group in (
        reconciliation.groupby(
            "screen"
        )
    )
}

for (
    screen,
    config,
) in CONFIG.items():
    observed_count = len(
        fastqs_by_screen.get(
            screen,
            set(),
        )
    )

    if (
        observed_count
        != config[
            "expected_fastq_count"
        ]
    ):
        raise RuntimeError(
            f"{screen}: expected "
            f"{config['expected_fastq_count']} "
            "FASTQ SequenceFiles, observed "
            f"{observed_count}"
        )


# =====================================================================
# 4. Verify CuratedSets, guide-list membership, and current lineage
# =====================================================================

backup_records: dict[
    str,
    dict,
] = {}

plan_rows: list[dict] = []
payloads: dict[
    str,
    dict,
] = {}
expected_final: dict[
    str,
    dict,
] = {}

for (
    screen,
    config,
) in CONFIG.items():

    analysis_set_accession = (
        config["analysis_set"]
    )
    curated_set_accession = (
        config["curated_set"]
    )
    guide_file_accession = (
        config["guide_list_file"]
    )
    guide_quant_accession = (
        config[
            "guide_quantifications"
        ]
    )
    differential_accession = (
        config[
            "differential_quantifications"
        ]
    )

    analysis_edit = get_record(
        analysis_set_accession,
        frame="edit",
    )
    analysis_object = get_record(
        analysis_set_accession,
        frame="object",
    )
    curated_object = get_record(
        curated_set_accession,
        frame="object",
    )
    guide_file_object = get_record(
        guide_file_accession,
        frame="object",
    )
    guide_quant_edit = get_record(
        guide_quant_accession,
        frame="edit",
    )
    guide_quant_object = get_record(
        guide_quant_accession,
        frame="object",
    )
    differential_edit = get_record(
        differential_accession,
        frame="edit",
    )

    backup_records.update(
        {
            analysis_set_accession: (
                analysis_edit
            ),
            curated_set_accession: (
                curated_object
            ),
            guide_file_accession: (
                guide_file_object
            ),
            guide_quant_accession: (
                guide_quant_edit
            ),
            differential_accession: (
                differential_edit
            ),
        }
    )

    if (
        "AnalysisSet"
        not in record_types(
            analysis_object
        )
    ):
        raise RuntimeError(
            f"{analysis_set_accession} "
            "is not an AnalysisSet"
        )

    if (
        "CuratedSet"
        not in record_types(
            curated_object
        )
    ):
        raise RuntimeError(
            f"{curated_set_accession} "
            "is not a CuratedSet: "
            f"{sorted(record_types(curated_object))}"
        )

    guide_file_set = extract_accession(
        guide_file_object.get(
            "file_set",
            "",
        ),
        r"IGVFDS",
    )

    if (
        guide_file_set
        != curated_set_accession
    ):
        raise RuntimeError(
            f"{guide_file_accession} does not "
            "belong to the expected CuratedSet.\n"
            f"Expected: {curated_set_accession}\n"
            f"Observed: {guide_file_set or '<missing>'}"
        )

    guide_file_at_id = canonical_at_id(
        guide_file_object,
        guide_file_accession,
    )
    curated_set_at_id = canonical_at_id(
        curated_object,
        curated_set_accession,
    )

    # Verify differential-result lineage remains exact.
    differential_derived = (
        link_accessions(
            differential_edit.get(
                "derived_from",
                [],
            ),
            prefix="IGVFFI",
        )
    )

    if differential_derived != {
        guide_quant_accession
    }:
        raise RuntimeError(
            f"{differential_accession} has "
            "unexpected derived_from lineage.\n"
            f"Expected: {[guide_quant_accession]}\n"
            f"Observed: {sorted(differential_derived)}"
        )

    # ----------------------------
    # Guide-quantification lineage
    # ----------------------------
    current_derived_items = (
        guide_quant_edit.get(
            "derived_from",
            [],
        )
    )
    current_derived_map = link_map(
        current_derived_items,
        prefix="IGVFFI",
    )
    current_derived = set(
        current_derived_map
    )

    expected_fastqs = (
        fastqs_by_screen[
            screen
        ]
    )
    expected_derived = (
        expected_fastqs
        | {
            guide_file_accession
        }
    )

    missing_fastqs = (
        expected_fastqs
        - current_derived
    )
    unknown_derived = (
        current_derived
        - expected_derived
    )

    if missing_fastqs:
        raise RuntimeError(
            f"{guide_quant_accession} is "
            f"missing {len(missing_fastqs)} "
            "expected FASTQ derived_from links. "
            "This cell only adds guide-list "
            "provenance; rerun the prior FASTQ "
            "lineage repair first."
        )

    if unknown_derived:
        raise RuntimeError(
            f"{guide_quant_accession} contains "
            "unexpected derived_from files: "
            f"{sorted(unknown_derived)}"
        )

    final_derived_links = [
        current_derived_map[
            accession
        ]
        for accession in sorted(
            expected_fastqs
        )
    ] + [
        guide_file_at_id
    ]

    payloads[
        guide_quant_accession
    ] = {
        "derived_from": (
            final_derived_links
        ),
    }

    expected_final[
        guide_quant_accession
    ] = {
        "property": "derived_from",
        "prefix": "IGVFFI",
        "values": expected_derived,
    }

    # ----------------------------
    # AnalysisSet input_file_sets
    # ----------------------------
    current_input_sets = (
        link_accessions(
            analysis_edit.get(
                "input_file_sets",
                [],
            ),
            prefix="IGVFDS",
        )
    )

    expected_measurement_sets = set(
        config[
            "measurement_sets"
        ]
    )
    expected_input_sets = (
        expected_measurement_sets
        | {
            curated_set_accession
        }
    )

    missing_measurement_sets = (
        expected_measurement_sets
        - current_input_sets
    )
    unknown_input_sets = (
        current_input_sets
        - expected_input_sets
    )

    if missing_measurement_sets:
        raise RuntimeError(
            f"{analysis_set_accession} is "
            "missing expected MeasurementSets: "
            f"{sorted(missing_measurement_sets)}"
        )

    if unknown_input_sets:
        raise RuntimeError(
            f"{analysis_set_accession} contains "
            "unexpected input_file_sets: "
            f"{sorted(unknown_input_sets)}"
        )

    canonical_measurement_set_ids = []

    for accession in sorted(
        expected_measurement_sets
    ):
        measurement_object = get_record(
            accession,
            frame="object",
        )

        if (
            "MeasurementSet"
            not in record_types(
                measurement_object
            )
        ):
            raise RuntimeError(
                f"{accession} is not a "
                "MeasurementSet"
            )

        canonical_measurement_set_ids.append(
            canonical_at_id(
                measurement_object,
                accession,
            )
        )

    final_input_set_links = (
        canonical_measurement_set_ids
        + [
            curated_set_at_id
        ]
    )

    payloads[
        analysis_set_accession
    ] = {
        "input_file_sets": (
            final_input_set_links
        ),
    }

    expected_final[
        analysis_set_accession
    ] = {
        "property": "input_file_sets",
        "prefix": "IGVFDS",
        "values": expected_input_sets,
    }

    plan_rows.extend(
        [
            {
                "screen": screen,
                "record_id": (
                    guide_quant_accession
                ),
                "property": "derived_from",
                "current_count": len(
                    current_derived
                ),
                "final_count": len(
                    expected_derived
                ),
                "guide_input": (
                    guide_file_accession
                ),
                "guide_input_file_set": (
                    guide_file_set
                ),
                "action": (
                    "SKIP_ALREADY_CORRECT"
                    if current_derived
                    == expected_derived
                    else "ADD_GUIDE_LIST_FILE"
                ),
            },
            {
                "screen": screen,
                "record_id": (
                    analysis_set_accession
                ),
                "property": (
                    "input_file_sets"
                ),
                "current_count": len(
                    current_input_sets
                ),
                "final_count": len(
                    expected_input_sets
                ),
                "guide_input": (
                    guide_file_accession
                ),
                "guide_input_file_set": (
                    curated_set_accession
                ),
                "action": (
                    "SKIP_ALREADY_CORRECT"
                    if current_input_sets
                    == expected_input_sets
                    else "ADD_CURATED_SET"
                ),
            },
        ]
    )


# =====================================================================
# 5. Save preview evidence
# =====================================================================

plan = pd.DataFrame(
    plan_rows
)

plan_path = (
    OUTPUT_ROOT
    / "curated_guide_input_patch_plan.tsv"
)

plan.to_csv(
    plan_path,
    sep="\t",
    index=False,
    lineterminator="\n",
)

backup_path = (
    OUTPUT_ROOT
    / "records_before_curated_guide_patch.json"
)

backup_path.write_text(
    json.dumps(
        backup_records,
        indent=2,
        ensure_ascii=False,
    )
    + "\n",
    encoding="utf-8",
)

payload_path = (
    OUTPUT_ROOT
    / "curated_guide_patch_payloads.json"
)

payload_path.write_text(
    json.dumps(
        payloads,
        indent=2,
        ensure_ascii=False,
    )
    + "\n",
    encoding="utf-8",
)

print(
    "Curated guide-input patch plan:"
)

print(
    plan.to_string(
        index=False
    )
)

print(
    "\nSaved plan:",
    plan_path,
)

print(
    "Saved backup:",
    backup_path,
)

print(
    "Saved payloads:",
    payload_path,
)


# =====================================================================
# 6. Preview gate or production PATCH
# =====================================================================

if not APPLY_PATCH:

    print(
        "\nPREVIEW ONLY — no Portal records "
        "were changed.\n\n"
        "Expected final counts:\n"
        "  GW guide quantifications derived_from: 363\n"
        "  GW AnalysisSet input_file_sets:          4\n"
        "  Validation guide derived_from:           42\n"
        "  Validation AnalysisSet input_file_sets:  7\n\n"
        "To apply, change the top of this cell to:\n\n"
        "APPLY_PATCH = True\n"
        f'CONFIRMATION = '
        f'"{REQUIRED_CONFIRMATION}"\n\n'
        "and rerun this same cell."
    )

else:

    if (
        CONFIRMATION
        != REQUIRED_CONFIRMATION
    ):
        raise RuntimeError(
            "Production confirmation is "
            "missing or incorrect.\n"
            f'Required:\nCONFIRMATION = '
            f'"{REQUIRED_CONFIRMATION}"'
        )

    results: list[dict] = []

    # -------------------------------------------------------------
    # Stage 1: represent the CuratedSets by adding their guide files
    # to the guide-quantification derived_from arrays.
    # -------------------------------------------------------------
    for (
        screen,
        config,
    ) in CONFIG.items():

        record_id = config[
            "guide_quantifications"
        ]

        expectation = (
            expected_final[
                record_id
            ]
        )

        current_record = get_record(
            record_id,
            frame="edit",
        )

        current_values = (
            link_accessions(
                current_record.get(
                    expectation[
                        "property"
                    ],
                    [],
                ),
                prefix=expectation[
                    "prefix"
                ],
            )
        )

        if (
            current_values
            == expectation["values"]
        ):
            results.append(
                {
                    "stage": 1,
                    "screen": screen,
                    "record_id": (
                        record_id
                    ),
                    "property": (
                        expectation[
                            "property"
                        ]
                    ),
                    "action": (
                        "SKIP_ALREADY_CORRECT"
                    ),
                    "patch_performed": (
                        False
                    ),
                    "observed_count": len(
                        current_values
                    ),
                    "verification": "PASS",
                    "response_json": "",
                }
            )
            continue

        response_payload = patch_record(
            record_id,
            payloads[
                record_id
            ],
        )

        verified_record = get_record(
            record_id,
            frame="edit",
        )

        verified_values = (
            link_accessions(
                verified_record.get(
                    expectation[
                        "property"
                    ],
                    [],
                ),
                prefix=expectation[
                    "prefix"
                ],
            )
        )

        if (
            verified_values
            != expectation["values"]
        ):
            raise RuntimeError(
                f"{record_id} verification "
                "failed after PATCH.\n"
                f"Expected count: "
                f"{len(expectation['values'])}\n"
                f"Observed count: "
                f"{len(verified_values)}\n"
                f"Missing: "
                f"{sorted(expectation['values'] - verified_values)}\n"
                f"Additional: "
                f"{sorted(verified_values - expectation['values'])}"
            )

        response_path = (
            OUTPUT_ROOT
            / (
                f"{record_id}_"
                "curated_guide_patch_response.json"
            )
        )

        response_path.write_text(
            json.dumps(
                response_payload,
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )

        results.append(
            {
                "stage": 1,
                "screen": screen,
                "record_id": record_id,
                "property": (
                    expectation[
                        "property"
                    ]
                ),
                "action": (
                    "ADD_GUIDE_LIST_FILE"
                ),
                "patch_performed": True,
                "observed_count": len(
                    verified_values
                ),
                "verification": "PASS",
                "response_json": str(
                    response_path
                ),
            }
        )

    # -------------------------------------------------------------
    # Stage 2: now that the guide file represents each CuratedSet in
    # derived_from, restore the CuratedSet to AnalysisSet inputs.
    # -------------------------------------------------------------
    for (
        screen,
        config,
    ) in CONFIG.items():

        record_id = config[
            "analysis_set"
        ]

        expectation = (
            expected_final[
                record_id
            ]
        )

        current_record = get_record(
            record_id,
            frame="edit",
        )

        current_values = (
            link_accessions(
                current_record.get(
                    expectation[
                        "property"
                    ],
                    [],
                ),
                prefix=expectation[
                    "prefix"
                ],
            )
        )

        if (
            current_values
            == expectation["values"]
        ):
            results.append(
                {
                    "stage": 2,
                    "screen": screen,
                    "record_id": (
                        record_id
                    ),
                    "property": (
                        expectation[
                            "property"
                        ]
                    ),
                    "action": (
                        "SKIP_ALREADY_CORRECT"
                    ),
                    "patch_performed": (
                        False
                    ),
                    "observed_count": len(
                        current_values
                    ),
                    "verification": "PASS",
                    "response_json": "",
                }
            )
            continue

        response_payload = patch_record(
            record_id,
            payloads[
                record_id
            ],
        )

        verified_record = get_record(
            record_id,
            frame="edit",
        )

        verified_values = (
            link_accessions(
                verified_record.get(
                    expectation[
                        "property"
                    ],
                    [],
                ),
                prefix=expectation[
                    "prefix"
                ],
            )
        )

        if (
            verified_values
            != expectation["values"]
        ):
            raise RuntimeError(
                f"{record_id} verification "
                "failed after PATCH.\n"
                f"Expected: "
                f"{sorted(expectation['values'])}\n"
                f"Observed: "
                f"{sorted(verified_values)}"
            )

        response_path = (
            OUTPUT_ROOT
            / (
                f"{record_id}_"
                "curated_set_patch_response.json"
            )
        )

        response_path.write_text(
            json.dumps(
                response_payload,
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )

        results.append(
            {
                "stage": 2,
                "screen": screen,
                "record_id": record_id,
                "property": (
                    expectation[
                        "property"
                    ]
                ),
                "action": (
                    "ADD_CURATED_SET"
                ),
                "patch_performed": True,
                "observed_count": len(
                    verified_values
                ),
                "verification": "PASS",
                "response_json": str(
                    response_path
                ),
            }
        )

    # -------------------------------------------------------------
    # Final saved reconciliation
    # -------------------------------------------------------------
    results_frame = pd.DataFrame(
        results
    ).sort_values(
        [
            "stage",
            "screen",
        ]
    )

    results_path = (
        OUTPUT_ROOT
        / "curated_guide_input_patch_results.tsv"
    )

    results_frame.to_csv(
        results_path,
        sep="\t",
        index=False,
        lineterminator="\n",
    )

    if not results_frame[
        "verification"
    ].eq("PASS").all():
        raise RuntimeError(
            "One or more final "
            "verifications failed"
        )

    marker_path = (
        OUTPUT_ROOT
        / (
            "CURATED_GUIDE_INPUTS_"
            "PATCHED_4_OF_4.txt"
        )
    )

    marker_path.write_text(
        "Both guide-list files are represented "
        "in guide-quantification derived_from, "
        "and both CuratedSets are represented "
        "in AnalysisSet input_file_sets.\n"
        f"Results: {results_path}\n",
        encoding="utf-8",
    )

    print(
        "\nProduction results:"
    )

    print(
        results_frame.to_string(
            index=False
        )
    )

    print(
        "\nResults:",
        results_path,
    )

    print(
        "Completion marker:",
        marker_path,
    )

    print(
        "\nVerified final counts:\n"
        "  GW guide quantifications derived_from: 363\n"
        "  GW AnalysisSet input_file_sets:          4\n"
        "  Validation guide derived_from:           42\n"
        "  Validation AnalysisSet input_file_sets:  7\n\n"
        "Portal audit badges may take several "
        "minutes to reindex."
    )

```

### Text output (stream)

```text
Curated guide-input patch plan:
     screen      record_id        property  current_count  final_count    guide_input guide_input_file_set              action
genome_wide IGVFFI8424MVZE    derived_from            362          363 IGVFFI9629WVUU       IGVFDS4419BVVU ADD_GUIDE_LIST_FILE
genome_wide IGVFDS9629VXQK input_file_sets              3            4 IGVFFI9629WVUU       IGVFDS4419BVVU     ADD_CURATED_SET
 validation IGVFFI5157STFU    derived_from             41           42 IGVFFI6745HXWN       IGVFDS8881NOFE ADD_GUIDE_LIST_FILE
 validation IGVFDS4708PCVX input_file_sets              6            7 IGVFFI6745HXWN       IGVFDS8881NOFE     ADD_CURATED_SET

Saved plan: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/curated_guide_input_fix/curated_guide_input_patch_plan.tsv
Saved backup: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/curated_guide_input_fix/records_before_curated_guide_patch.json
Saved payloads: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/curated_guide_input_fix/curated_guide_patch_payloads.json

Production results:
 stage      screen      record_id        property              action  patch_performed  observed_count verification                                                                                                                                        response_json
     1 genome_wide IGVFFI8424MVZE    derived_from ADD_GUIDE_LIST_FILE             True             363         PASS /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/curated_guide_input_fix/IGVFFI8424MVZE_curated_guide_patch_response.json
     1  validation IGVFFI5157STFU    derived_from ADD_GUIDE_LIST_FILE             True              42         PASS /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/curated_guide_input_fix/IGVFFI5157STFU_curated_guide_patch_response.json
     2 genome_wide IGVFDS9629VXQK input_file_sets     ADD_CURATED_SET             True               4         PASS   /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/curated_guide_input_fix/IGVFDS9629VXQK_curated_set_patch_response.json
     2  validation IGVFDS4708PCVX input_file_sets     ADD_CURATED_SET             True               7         PASS   /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/curated_guide_input_fix/IGVFDS4708PCVX_curated_set_patch_response.json

Results: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/curated_guide_input_fix/curated_guide_input_patch_results.tsv
Completion marker: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_crispri_macs_lineage_patch/curated_guide_input_fix/CURATED_GUIDE_INPUTS_PATCHED_4_OF_4.txt

Verified final counts:
  GW guide quantifications derived_from: 363
  GW AnalysisSet input_file_sets:          4
  Validation guide derived_from:           42
  Validation AnalysisSet input_file_sets:  7

Portal audit badges may take several minutes to reindex.

```

## Code cell 32

```python

```

