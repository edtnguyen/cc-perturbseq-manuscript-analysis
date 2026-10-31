# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/submit_igvf_file_format_documents_and_attach.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/submit_igvf_file_format_documents_and_attach.ipynb`
- SHA-256: `324015a01c220054a2fa5360ebad84016bf675e653dbcef278efd05fc88d89ba`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Submit IGVF file-format specification Documents and attach them to analysis files

This notebook performs the workflow requested by Ian:

1. locate the four PDF file-format specifications under
   `$OAK/Users/tri/tmp/igvf_submit_validation`;
2. POST each PDF as an IGVF `Document` with
   `document_type = file format specification`;
3. retrieve and verify the four Document records;
4. PATCH the corresponding `TabularFile` records through
   `file_format_specifications`;
5. export a final reconciliation table with Document IDs, target-file
   accessions, Portal links, and verification status.

The notebook is restart-safe:

- an existing matching Document is skipped;
- an analysis file already linked to the correct Document is skipped;
- a target file that does not yet exist is reported as `WAITING_FOR_FILE`;
- production work is gated behind separate confirmation tokens;
- each POST/PATCH uses a one-row TSV and a separate log.

The patch operation **extends** the existing
`file_format_specifications` array; it does not overwrite unrelated
Document links.

## Markdown cell 1

## Version 2 — submit Document payloads as JSON, not TSV

The first version serialized the nested `attachment` object through
`pandas.to_csv()`. That produced CSV-style doubled quotes inside the
TSV field, for example:

```text
"{""path"":""/oak/.../specification.pdf""}"
```

`iu_register.py` reads TSV rows by splitting on tab characters and
then calls `json.loads()` directly on object-valued fields. It does
not undo CSV quote escaping, so the `attachment` value failed with
`JSONDecodeError: Extra data`.

This version submits each Document from a native `.json` file. The
attachment is therefore represented as the actual JSON object:

```json
{"path": "/oak/.../specification.pdf"}
```

TabularFile PATCH operations remain TSV because they do not contain
nested JSON objects.

## Markdown cell 2

## 1. Configuration

## Code cell 3

```python
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import subprocess
import time
from getpass import getpass
from pathlib import Path
from typing import Any, Sequence

import pandas as pd
import requests

# Sherlock directory supplied by the user.
OAK_ROOT = Path(
    os.environ.get(
        "OAK",
        "/oak/stanford/groups/engreitz",
    )
)
BASE_DIR = OAK_ROOT / "Users/tri/tmp/igvf_submit_validation"
WORK_DIR = BASE_DIR / "_igvf_document_submission"
LOG_DIR = WORK_DIR / "logs"
ROW_DIR = WORK_DIR / "one_row_tsv"

for directory in (WORK_DIR, LOG_DIR, ROW_DIR):
    directory.mkdir(parents=True, exist_ok=True)

AWARD = "/awards/HG011972/"
LAB = "/labs/jesse-engreitz/"

IGVF_PYTHON = Path(
    "/oak/stanford/groups/engreitz/Users/tri/"
    "envs/igvf_utils/bin/python"
)
IU_REGISTER = Path(
    "/oak/stanford/groups/engreitz/Users/tri/"
    "envs/igvf_utils/igvf_utils/igvf_utils/"
    "MetaDataRegistration/iu_register.py"
)

# The four analysis files that should receive the format Documents.
# Replace an alias with an IGVFFI accession when one is known.
SPECIFICATIONS = {
    "validation_guide": {
        "pdf_name": (
            "WTC11_EC_validation_screen_"
            "guide_quantifications_file_format_specification.pdf"
        ),
        "document_alias": (
            "jesse-engreitz:document_WTC11-EC_validation_screen_"
            "guide_quantifications_file_format_v1"
        ),
        "description": (
            "File-format specification for guide-level read-count "
            "quantifications from the focused WTC11-EC validation CRISPRi "
            "screen."
        ),
        "target_file_identifier": "IGVFFI5157STFU",
        "target_file_alias": (
            "jesse-engreitz:tabular_file_WTC11-EC_validation_screen_"
            "guide_quantifications"
        ),
    },
    "validation_differential": {
        "pdf_name": (
            "WTC11_EC_validation_screen_"
            "differential_element_quantifications_"
            "file_format_specification.pdf"
        ),
        "document_alias": (
            "jesse-engreitz:document_WTC11-EC_validation_screen_"
            "differential_element_quantifications_file_format_v1"
        ),
        "description": (
            "File-format specification for promoter-element MAGeCK "
            "differential quantifications from the focused WTC11-EC "
            "validation CRISPRi screen."
        ),
        "target_file_identifier": (
            "jesse-engreitz:tabular_file_WTC11-EC_validation_screen_"
            "D3_CDH5pos_vs_CDH5neg_"
            "differential_element_quantifications"
        ),
        "target_file_alias": (
            "jesse-engreitz:tabular_file_WTC11-EC_validation_screen_"
            "D3_CDH5pos_vs_CDH5neg_"
            "differential_element_quantifications"
        ),
    },
    "genome_wide_guide": {
        "pdf_name": (
            "WTC11_EC_genome_wide_screen_"
            "guide_quantifications_file_format_specification.pdf"
        ),
        "document_alias": (
            "jesse-engreitz:document_WTC11-EC_genome_wide_screen_"
            "guide_quantifications_file_format_v1"
        ),
        "description": (
            "File-format specification for guide-level read-count "
            "quantifications from the genome-wide WTC11-EC CRISPRi screen."
        ),
        "target_file_identifier": (
            "jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_"
            "guide_quantifications"
        ),
        "target_file_alias": (
            "jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_"
            "guide_quantifications"
        ),
    },
    "genome_wide_differential": {
        "pdf_name": (
            "WTC11_EC_genome_wide_screen_"
            "differential_element_quantifications_"
            "file_format_specification.pdf"
        ),
        "document_alias": (
            "jesse-engreitz:document_WTC11-EC_genome_wide_screen_"
            "differential_element_quantifications_file_format_v1"
        ),
        "description": (
            "File-format specification for native element-level MAGeCK "
            "differential quantifications from the genome-wide WTC11-EC "
            "CRISPRi screen."
        ),
        "target_file_identifier": (
            "jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_"
            "D3_CDH5pos_vs_CDH5neg_"
            "differential_element_quantifications"
        ),
        "target_file_alias": (
            "jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_"
            "D3_CDH5pos_vs_CDH5neg_"
            "differential_element_quantifications"
        ),
    },
}

print("Base directory:", BASE_DIR)
print("Work directory:", WORK_DIR)
print("igvf_utils Python:", IGVF_PYTHON)
print("iu_register.py:", IU_REGISTER)
```

### Text output (stream)

```text
Base directory: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation
Work directory: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission
igvf_utils Python: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python
iu_register.py: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py

```

## Markdown cell 4

## 2. Utility functions

## Code cell 5

```python
def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        lineterminator="\n",
    )


def md5sum(
    path: Path,
    block_size: int = 16 * 1024 * 1024,
) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(block_size),
            b"",
        ):
            digest.update(block)
    return digest.hexdigest()


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


def safe_stem(value: str) -> str:
    return re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        str(value),
    ).strip("_")


def split_link_items(value: object) -> list[object]:
    if value is None:
        return []
    if isinstance(value, str):
        return [
            item.strip()
            for item in value.split(",")
            if item.strip()
        ]
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]


def terminal_identifier(value: object) -> str:
    if isinstance(value, dict):
        value = (
            value.get("accession")
            or value.get("@id")
            or value.get("uuid")
            or ""
        )
    text = str(value or "").strip().strip("/")
    return text.rsplit("/", 1)[-1] if text else ""


def direct_link_tokens(value: object) -> set[str]:
    tokens: set[str] = set()
    for item in split_link_items(value):
        if isinstance(item, dict):
            for key in ("accession", "uuid", "@id"):
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
            if ":" in text and not text.startswith("/"):
                tokens.add(text)
    return tokens


def portal_urls(identifier: str) -> list[str]:
    token = str(identifier).strip().strip("/")
    base = f"https://api.data.igvf.org/{token}/"
    return [
        base + "?format=json&datastore=database",
        base + "?format=json&datastore=database&frame=object",
        base + "?format=json&datastore=database&frame=edit",
    ]


def portal_get(
    identifier: str,
    *,
    timeout: int = 45,
) -> dict | None:
    auth = None
    if (
        os.environ.get("IGVF_API_KEY")
        and os.environ.get("IGVF_SECRET_KEY")
    ):
        auth = (
            os.environ["IGVF_API_KEY"],
            os.environ["IGVF_SECRET_KEY"],
        )

    errors: list[str] = []
    for url in portal_urls(identifier):
        response = requests.get(
            url,
            auth=auth,
            headers={"Accept": "application/json"},
            timeout=timeout,
        )
        if response.status_code == 404:
            continue
        if response.ok:
            return response.json()
        errors.append(
            f"{response.status_code} {url}: "
            f"{response.text[:300]}"
        )

    if errors:
        raise RuntimeError(
            "Portal GET failed:\n" + "\n".join(errors)
        )
    return None


def get_with_retry(
    identifier: str,
    *,
    attempts: int = 6,
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
            identities.update(
                direct_link_tokens(record)
            )

        _IDENTITY_CACHE[token] = identities
        resolved.update(identities)

    return resolved


def links_equal(
    left: object,
    right: object,
) -> bool:
    left_tokens = resolved_link_tokens(left)
    right_tokens = resolved_link_tokens(right)
    return bool(
        left_tokens
        and right_tokens
        and left_tokens & right_tokens
    )


def link_set_contains(
    container: object,
    required: object,
) -> bool:
    existing_items = split_link_items(container)
    required_items = split_link_items(required)
    return all(
        any(
            links_equal(existing, wanted)
            for existing in existing_items
        )
        for wanted in required_items
    )


def portal_record_id(record: dict | None) -> str:
    record = record or {}
    return str(
        record.get("@id")
        or record.get("accession")
        or record.get("uuid")
        or ""
    )


def portal_record_url(record: dict | None) -> str:
    record = record or {}
    portal_id = str(record.get("@id") or "")
    if portal_id:
        return "https://data.igvf.org" + portal_id

    accession = str(record.get("accession") or "")
    if accession:
        return (
            "https://data.igvf.org/"
            f"tabular-files/{accession}/"
        )
    return ""


def row_without_empty_values(row: dict) -> dict:
    return {
        key: value
        for key, value in row.items()
        if str(value).strip() != ""
    }


def run_iu_register(
    *,
    profile: str,
    input_tsv: Path,
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
        raise RuntimeError(
            "IGVF credentials are not set"
        )

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    mode = "dryrun" if dry_run else "submit"
    log_path = (
        LOG_DIR
        / f"{safe_stem(label)}_{mode}_{timestamp}.log"
    )

    command = [
        str(IGVF_PYTHON),
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
        r"\b401 Client Error\b",
        r"\b403 Client Error\b",
        r"\b409 Client Error\b",
        r"\b422 Client Error\b",
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
    if not success:
        raise RuntimeError(
            f"igvf_utils failed; review {log_path}"
        )

    return {
        "profile": profile,
        "input_tsv": str(input_tsv),
        "dry_run": dry_run,
        "patch": patch,
        "log_path": str(log_path),
        "success": success,
    }


print("Utility functions loaded.")
```

### Text output (stream)

```text
Utility functions loaded.

```

## Markdown cell 6

## 3. Validate the four local PDF attachments and build native JSON Document payloads

## Code cell 7

```python
if not BASE_DIR.is_dir():
    raise FileNotFoundError(
        f"Submission directory does not exist: {BASE_DIR}"
    )

inventory_rows: list[dict] = []
document_payloads: list[dict] = []

for key, spec in SPECIFICATIONS.items():
    pdf_path = BASE_DIR / spec["pdf_name"]

    if not pdf_path.is_file():
        raise FileNotFoundError(
            f"Missing PDF for {key}: {pdf_path}"
        )
    if pdf_path.stat().st_size == 0:
        raise RuntimeError(
            f"PDF is empty: {pdf_path}"
        )

    with pdf_path.open("rb") as handle:
        header = handle.read(5)
    if header != b"%PDF-":
        raise RuntimeError(
            f"File does not have a PDF header: {pdf_path}"
        )

    inventory_rows.append(
        {
            "key": key,
            "pdf_path": str(pdf_path),
            "size_bytes": pdf_path.stat().st_size,
            "md5sum": md5sum(pdf_path),
            "sha256": sha256sum(pdf_path),
            "document_alias": spec["document_alias"],
            "target_file_identifier": (
                spec["target_file_identifier"]
            ),
        }
    )

    # JSON input must use schema-native types:
    # - aliases is an array
    # - attachment is an object, not a JSON string embedded in TSV
    payload = {
        "aliases": [spec["document_alias"]],
        "document_type": "file format specification",
        "description": spec["description"],
        "attachment": {"path": str(pdf_path)},
        "lab": LAB,
        "award": AWARD,
    }

    # Immediate JSON round-trip validation.
    payload = json.loads(
        json.dumps(
            payload,
            ensure_ascii=False,
            separators=(",", ":"),
        )
    )
    if payload["attachment"] != {"path": str(pdf_path)}:
        raise RuntimeError(
            f"Attachment JSON round-trip failed for {key}"
        )

    document_payloads.append(payload)

pdf_inventory = pd.DataFrame(inventory_rows)
PDF_INVENTORY_PATH = (
    WORK_DIR / "pdf_attachment_inventory.tsv"
)
write_tsv(pdf_inventory, PDF_INVENTORY_PATH)

# This JSON file is suitable for direct use with iu_register.py.
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

# A flattened TSV is written only for human review. Do not submit it.
review_rows = []
for payload in document_payloads:
    review_rows.append(
        {
            "aliases": ",".join(payload["aliases"]),
            "document_type": payload["document_type"],
            "description": payload["description"],
            "attachment_path": payload["attachment"]["path"],
            "lab": payload["lab"],
            "award": payload["award"],
        }
    )

DOCUMENT_REVIEW_TSV_PATH = (
    WORK_DIR / "document_file_format_specifications_review_only.tsv"
)
write_tsv(
    pd.DataFrame(review_rows),
    DOCUMENT_REVIEW_TSV_PATH,
)

# Validate the aggregate JSON from disk.
loaded_payloads = json.loads(
    DOCUMENT_METADATA_PATH.read_text(encoding="utf-8")
)
if loaded_payloads != document_payloads:
    raise RuntimeError(
        "Document metadata JSON changed during disk round-trip"
    )

print(pdf_inventory.to_string(index=False))
print("\nSubmission JSON:")
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
                     key                                                                                                                                                          pdf_path  size_bytes                           md5sum                                                           sha256                                                                                          document_alias                                                                                             target_file_identifier
        validation_guide                 /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/WTC11_EC_validation_screen_guide_quantifications_file_format_specification.pdf        5791 0cb375eb11bf64ee66dfe2b388e3c89c 65b47235513478428791fd79a73bac4c8b3fddd3d3a41df730c0d780bf35bcca                 jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1                                                                                                     IGVFFI5157STFU
 validation_differential  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/WTC11_EC_validation_screen_differential_element_quantifications_file_format_specification.pdf        8100 23f4f9f92787f2367355ddb1f53b1ab1 c94d4819c7825f489f076d89a659eabcec7e4d6684bdbca58ffa9c429b1f0283  jesse-engreitz:document_WTC11-EC_validation_screen_differential_element_quantifications_file_format_v1  jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications
       genome_wide_guide                /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/WTC11_EC_genome_wide_screen_guide_quantifications_file_format_specification.pdf        5419 371b8c2433ecb8f6fc6a07d106abcba3 d46bbe16af1e820c5656498d11510536fe780bab4f69e4d6118c9f88c98f2423                jesse-engreitz:document_WTC11-EC_genome_wide_screen_guide_quantifications_file_format_v1                                      jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications
genome_wide_differential /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/WTC11_EC_genome_wide_screen_differential_element_quantifications_file_format_specification.pdf        6275 c27533cfc99a767ca764aa3f064bf56d 9f94de79e6c7abe7311ab5cfe3f0f25ee1e7f878747c4fda44372d2b60a7dbf6 jesse-engreitz:document_WTC11-EC_genome_wide_screen_differential_element_quantifications_file_format_v1 jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications

Submission JSON:
/oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/document_file_format_specifications.json
[
  {
    "aliases": [
      "jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1"
    ],
    "document_type": "file format specification",
    "description": "File-format specification for guide-level read-count quantifications from the focused WTC11-EC validation CRISPRi screen.",
    "attachment": {
      "path": "/oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/WTC11_EC_validation_screen_guide_quantifications_file_format_specification.pdf"
    },
    "lab": "/labs/jesse-engreitz/",
    "award": "/awards/HG011972/"
  },
  {
    "aliases": [
      "jesse-engreitz:document_WTC11-EC_validation_screen_differential_element_quantifications_file_format_v1"
    ],
    "document_type": "file format specification",
    "description": "File-format specification for promoter-element MAGeCK differential quantifications from the focused WTC11-EC validation CRISPRi screen.",
    "attachment": {
      "path": "/oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/WTC11_EC_validation_screen_differential_element_quantifications_file_format_specification.pdf"
    },
    "lab": "/labs/jesse-engreitz/",
    "award": "/awards/HG011972/"
  },
  {
    "aliases": [
      "jesse-engreitz:document_WTC11-EC_genome_wide_screen_guide_quantifications_file_format_v1"
    ],
    "document_type": "file format specification",
    "description": "File-format specification for guide-level read-count quantifications from the genome-wide WTC11-EC CRISPRi screen.",
    "attachment": {
      "path": "/oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/WTC11_EC_genome_wide_screen_guide_quantifications_file_format_specification.pdf"
    },
    "lab": "/labs/jesse-engreitz/",
    "award": "/awards/HG011972/"
  },
  {
    "aliases": [
      "jesse-engreitz:document_WTC11-EC_genome_wide_screen_differential_element_quantifications_file_format_v1"
    ],
    "document_type": "file format specification",
    "description": "File-format specification for native element-level MAGeCK differential quantifications from the genome-wide WTC11-EC CRISPRi screen.",
    "attachment": {
      "path": "/oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/WTC11_EC_genome_wide_screen_differential_element_quantifications_file_format_specification.pdf"
    },
    "lab": "/labs/jesse-engreitz/",
    "award": "/awards/HG011972/"
  }
]

Review-only TSV:
/oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/document_file_format_specifications_review_only.tsv

```

## Markdown cell 8

## 4. Enter IGVF credentials

Run the first cell to define the secure helpers, then run the second cell to
open the hidden API-key and secret-key prompts. Repeat both cells after a
kernel restart.

## Code cell 9

```python
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
                f"{name} has surrounding quotes. "
                "Paste the raw value."
            )

    os.environ["IGVF_API_KEY"] = api_key
    os.environ["IGVF_SECRET_KEY"] = secret_key

    del api_key, secret_key
    print(
        "IGVF credentials are set for "
        "the current notebook kernel."
    )


def show_credential_status() -> None:
    for variable in (
        "IGVF_API_KEY",
        "IGVF_SECRET_KEY",
    ):
        value = os.environ.get(variable, "")
        print(
            f"{variable}: "
            f"{'SET' if value else 'MISSING'}, "
            f"length={len(value)}, "
            f"whitespace={value != value.strip()}"
        )
```

## Code cell 10

```python
# Run this cell before Portal preflight.
enter_igvf_credentials()
show_credential_status()
```

### Text output (stream)

```text
IGVF credentials are set for the current notebook kernel.
IGVF_API_KEY: <REDACTED> length=8, whitespace=False
IGVF_SECRET_KEY: SET, length=16, whitespace=False

```

## Markdown cell 11

## 5. Build the Document submission plan

Existing matching Document aliases are skipped. A pre-existing alias with a
different `document_type` is treated as a blocking conflict.

## Code cell 12

```python
if not (
    os.environ.get("IGVF_API_KEY")
    and os.environ.get("IGVF_SECRET_KEY")
):
    raise RuntimeError(
        "Run the credential-entry cell first"
    )

document_plan_rows: list[dict] = []

for key, spec in SPECIFICATIONS.items():
    alias = spec["document_alias"]

    try:
        record = get_with_retry(alias)
        lookup_error = ""
    except Exception as error:
        record = None
        lookup_error = (
            f"{type(error).__name__}: {error}"
        )

    if lookup_error:
        action = "BLOCK_LOOKUP_ERROR"
    elif record is None:
        action = "POST_NEW"
    elif (
        str(record.get("document_type", ""))
        != "file format specification"
    ):
        action = "BLOCK_EXISTING_DOCUMENT_CONFLICT"
    else:
        action = "SKIP_EXISTING"

    attachment = (
        record.get("attachment", {})
        if record
        else {}
    )

    document_plan_rows.append(
        {
            "key": key,
            "document_alias": alias,
            "action": action,
            "portal_id": str(
                (record or {}).get("@id", "")
            ),
            "portal_uuid": str(
                (record or {}).get("uuid", "")
            ),
            "portal_document_type": str(
                (record or {}).get(
                    "document_type",
                    "",
                )
            ),
            "portal_attachment_download": str(
                attachment.get("download", "")
                if isinstance(attachment, dict)
                else ""
            ),
            "portal_url": portal_record_url(record),
            "lookup_error": lookup_error,
        }
    )

document_plan = pd.DataFrame(
    document_plan_rows
)
DOCUMENT_PLAN_PATH = (
    WORK_DIR / "document_submission_plan.tsv"
)
write_tsv(document_plan, DOCUMENT_PLAN_PATH)

print(document_plan.to_string(index=False))
print("\nSaved:", DOCUMENT_PLAN_PATH)

blocked = document_plan.loc[
    document_plan["action"].str.startswith(
        "BLOCK_"
    )
]
if not blocked.empty:
    raise RuntimeError(
        "The Document plan contains a blocking "
        f"condition. Review {DOCUMENT_PLAN_PATH}"
    )
```

### Text output (stream)

```text
                     key                                                                                          document_alias   action portal_id portal_uuid portal_document_type portal_attachment_download portal_url lookup_error
        validation_guide                 jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1 POST_NEW                                                                                              
 validation_differential  jesse-engreitz:document_WTC11-EC_validation_screen_differential_element_quantifications_file_format_v1 POST_NEW                                                                                              
       genome_wide_guide                jesse-engreitz:document_WTC11-EC_genome_wide_screen_guide_quantifications_file_format_v1 POST_NEW                                                                                              
genome_wide_differential jesse-engreitz:document_WTC11-EC_genome_wide_screen_differential_element_quantifications_file_format_v1 POST_NEW                                                                                              

Saved: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/document_submission_plan.tsv

```

## Markdown cell 13

## 6. Create one-row Document POST JSON files

Each Document is submitted from a native JSON object. The older
`POST_document_*.tsv` files are removed to prevent accidental reuse.

## Code cell 14

```python
metadata_by_alias = {
    payload["aliases"][0]: payload
    for payload in document_payloads
}

# Remove unsafe operation TSVs from version 1. They contain nested JSON
# encoded with CSV quote doubling and must not be submitted.
removed_legacy_files = []
for legacy_path in ROW_DIR.glob("POST_document_*.tsv"):
    legacy_path.unlink()
    removed_legacy_files.append(str(legacy_path))

if removed_legacy_files:
    print("Removed legacy Document POST TSVs:")
    for legacy_path in removed_legacy_files:
        print(" ", legacy_path)

document_operation_rows: list[dict] = []

for plan in document_plan.to_dict("records"):
    key = plan["key"]
    action = plan["action"]
    alias = plan["document_alias"]

    if action == "POST_NEW":
        json_path = (
            ROW_DIR
            / f"POST_document_{safe_stem(key)}.json"
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

        # Validate the exact file that iu_register.py will read.
        reloaded = json.loads(
            json_path.read_text(encoding="utf-8")
        )
        if reloaded != payload:
            raise RuntimeError(
                f"JSON round-trip failed for {json_path}"
            )
        if not isinstance(
            reloaded.get("attachment"),
            dict,
        ):
            raise RuntimeError(
                f"attachment is not a JSON object in {json_path}"
            )
        if set(
            reloaded["attachment"]
        ) != {"path"}:
            raise RuntimeError(
                "attachment must contain exactly one key, 'path': "
                f"{json_path}"
            )
    else:
        json_path = None

    document_operation_rows.append(
        {
            "key": key,
            "document_alias": alias,
            "action": action,
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
    WORK_DIR / "document_submission_operations.tsv"
)
write_tsv(
    document_operations,
    DOCUMENT_OPERATIONS_PATH,
)

print(
    document_operations.to_string(index=False)
)

# Display the first operation payload for visual confirmation.
first_json = next(
    (
        Path(value)
        for value in document_operations["input_file"]
        if str(value).strip()
    ),
    None,
)
if first_json is not None:
    print("\nExample native JSON operation file:")
    print(first_json)
    print(first_json.read_text(encoding="utf-8"))
```

### Text output (stream)

```text
Removed legacy Document POST TSVs:
  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_validation_differential.tsv
  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_validation_guide.tsv
  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_genome_wide_differential.tsv
  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_genome_wide_guide.tsv
                     key                                                                                          document_alias   action                                                                                                                                           input_file input_format
        validation_guide                 jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1 POST_NEW         /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_validation_guide.json         json
 validation_differential  jesse-engreitz:document_WTC11-EC_validation_screen_differential_element_quantifications_file_format_v1 POST_NEW  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_validation_differential.json         json
       genome_wide_guide                jesse-engreitz:document_WTC11-EC_genome_wide_screen_guide_quantifications_file_format_v1 POST_NEW        /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_genome_wide_guide.json         json
genome_wide_differential jesse-engreitz:document_WTC11-EC_genome_wide_screen_differential_element_quantifications_file_format_v1 POST_NEW /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_genome_wide_differential.json         json

Example native JSON operation file:
/oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_validation_guide.json
{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1"
  ],
  "document_type": "file format specification",
  "description": "File-format specification for guide-level read-count quantifications from the focused WTC11-EC validation CRISPRi screen.",
  "attachment": {
    "path": "/oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/WTC11_EC_validation_screen_guide_quantifications_file_format_specification.pdf"
  },
  "lab": "/labs/jesse-engreitz/",
  "award": "/awards/HG011972/"
}


```

## Markdown cell 15

### 6A. Validate the exact JSON operation files

Run this cell before enabling dry runs. It verifies that every pending
Document operation is valid JSON and that `attachment` is a native
one-key object with a readable local PDF path.

## Code cell 16

```python
json_preflight_rows = []

for operation in document_operations.to_dict("records"):
    if operation["action"] == "SKIP_EXISTING":
        json_preflight_rows.append(
            {
                "key": operation["key"],
                "action": operation["action"],
                "input_file": "",
                "valid_json": True,
                "attachment_is_object": True,
                "attachment_exists": True,
                "status": "SKIP_EXISTING",
            }
        )
        continue

    input_path = Path(operation["input_file"])
    payload = json.loads(
        input_path.read_text(encoding="utf-8")
    )
    attachment = payload.get("attachment")
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

    json_preflight_rows.append(
        {
            "key": operation["key"],
            "action": operation["action"],
            "input_file": str(input_path),
            "valid_json": True,
            "attachment_is_object": attachment_is_object,
            "attachment_exists": attachment_exists,
            "status": (
                "PASS"
                if attachment_is_object
                and attachment_exists
                else "FAIL"
            ),
        }
    )

json_preflight = pd.DataFrame(
    json_preflight_rows
)
DOCUMENT_JSON_PREFLIGHT_PATH = (
    WORK_DIR / "document_json_preflight.tsv"
)
write_tsv(
    json_preflight,
    DOCUMENT_JSON_PREFLIGHT_PATH,
)

print(json_preflight.to_string(index=False))
print("\nSaved:", DOCUMENT_JSON_PREFLIGHT_PATH)

if not json_preflight["status"].isin(
    {"PASS", "SKIP_EXISTING"}
).all():
    raise RuntimeError(
        "At least one Document JSON operation failed preflight"
    )
```

### Text output (stream)

```text
                     key   action                                                                                                                                           input_file  valid_json  attachment_is_object  attachment_exists status
        validation_guide POST_NEW         /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_validation_guide.json        True                  True               True   PASS
 validation_differential POST_NEW  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_validation_differential.json        True                  True               True   PASS
       genome_wide_guide POST_NEW        /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_genome_wide_guide.json        True                  True               True   PASS
genome_wide_differential POST_NEW /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_genome_wide_differential.json        True                  True               True   PASS

Saved: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/document_json_preflight.tsv

```

## Markdown cell 17

## 7. Dry-run the missing Document records

Set `RUN_DOCUMENT_DRY_RUNS = True`, rerun the cell, and review every log.

## Code cell 18

```python
RUN_DOCUMENT_DRY_RUNS = True

DOCUMENT_DRY_RUN_MARKER = (
    WORK_DIR / "DOCUMENT_DRY_RUN_PASSED.json"
)

if not RUN_DOCUMENT_DRY_RUNS:
    print(
        "Document dry runs are disabled. "
        "Set RUN_DOCUMENT_DRY_RUNS=True "
        "and rerun this cell."
    )
else:
    results: list[dict] = []

    for operation in (
        document_operations.to_dict("records")
    ):
        if operation["action"] == "SKIP_EXISTING":
            results.append(
                {
                    "key": operation["key"],
                    "action": "SKIP_EXISTING",
                    "success": True,
                    "log_path": "",
                }
            )
            continue

        result = run_iu_register(
            profile="document",
            input_tsv=Path(
                operation["input_file"]
            ),
            dry_run=True,
            patch=False,
            label=(
                f"document_{operation['key']}"
            ),
        )
        results.append(
            {
                "key": operation["key"],
                "action": operation["action"],
                "success": result["success"],
                "log_path": result["log_path"],
            }
        )

    results_frame = pd.DataFrame(results)
    DOCUMENT_DRY_RUN_RESULTS_PATH = (
        WORK_DIR
        / "document_dry_run_results.tsv"
    )
    write_tsv(
        results_frame,
        DOCUMENT_DRY_RUN_RESULTS_PATH,
    )

    if not results_frame["success"].all():
        raise RuntimeError(
            "At least one Document dry run failed"
        )

    marker = {
        "timestamp": time.strftime(
            "%Y-%m-%dT%H:%M:%S%z"
        ),
        "document_metadata_sha256": (
            sha256sum(DOCUMENT_METADATA_PATH)
        ),
        "pdf_sha256": {
            row["key"]: row["sha256"]
            for row in inventory_rows
        },
    }
    DOCUMENT_DRY_RUN_MARKER.write_text(
        json.dumps(marker, indent=2)
        + "\n",
        encoding="utf-8",
    )

    print("All Document dry runs passed.")
    print("Marker:", DOCUMENT_DRY_RUN_MARKER)
```

### Text output (stream)

```text

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -d -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_validation_guide.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/document_validation_guide_dryrun_20260930_083417.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:34:30,826:iu_debug:	DRY RUN is enabled.
2026-09-30 08:34:30,826:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:34:30,826:iu_debug:	submission=True: In submission mode.
2026-09-30 08:34:31,470:iu_debug:	
IN post().
2026-09-30 08:34:32,117:iu_debug:	Validating the payload against the schema
2026-09-30 08:34:32,135:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_validation_screen_guide_quantifications_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOSAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDggMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL0NvbnRlbnRzIDEwIDAgUiAvTWVkaWFCb3ggWyAwIDAgNzkyIDYxMiBdIC9QYXJlbnQgOCAwIFIgL1Jlc291cmNlcyA8PAovRm9udCAxIDAgUiAvUHJvY1NldCBbIC9QREYgL1RleHQgL0ltYWdlQiAvSW1hZ2VDIC9JbWFnZUkgXQo+PiAvUm90YXRlIDAgL1RyYW5zIDw8Cgo+PiAKICAvVHlwZSAvUGFnZQo+PgplbmRvYmoKNiAwIG9iago8PAovUGFnZU1vZGUgL1VzZU5vbmUgL1BhZ2VzIDggMCBSIC9UeXBlIC9DYXRhbG9nCj4+CmVuZG9iago3IDAgb2JqCjw8Ci9BdXRob3IgKEVuZ3JlaXR6IExhYikgL0NyZWF0aW9uRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9DcmVhdG9yIChcKHVuc3BlY2lmaWVkXCkpIC9LZXl3b3JkcyAoKSAvTW9kRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9Qcm9kdWNlciAoUmVwb3J0TGFiIFBERiBMaWJyYXJ5IC0gXChvcGVuc291cmNlXCkpIAogIC9TdWJqZWN0IChJR1ZGIGZpbGUgZm9ybWF0IHNwZWNpZmljYXRpb24pIC9UaXRsZSAoRm9jdXNlZCB2YWxpZGF0aW9uIHNjcmVlbiBHdWlkZSBxdWFudGlmaWNhdGlvbnMgZmlsZSBmb3JtYXQgc3BlY2lmaWNhdGlvbikgL1RyYXBwZWQgL0ZhbHNlCj4+CmVuZG9iago4IDAgb2JqCjw8Ci9Db3VudCAyIC9LaWRzIFsgNCAwIFIgNSAwIFIgXSAvVHlwZSAvUGFnZXMKPj4KZW5kb2JqCjkgMCBvYmoKPDwKL0ZpbHRlciBbIC9BU0NJSTg1RGVjb2RlIC9GbGF0ZURlY29kZSBdIC9MZW5ndGggMjk3OQo+PgpzdHJlYW0KR2IhI14+QkFRMSZVc0pYZk0zT1JTaFhyMkdLKkZdQVolXDZJIWNwMUlIY11HVkk1ZlBdSTNwRWdbck8xKEViMScsIWdzU1IzbixSK0dkX2lfZFhEaSJxTFdEbylHY29wdEhNVjpEYmMlVlVyVjJnUD4jQiopQ2QtbSs9JGM2NzJ0IXAmWzJDbk1jODZmOCM6MTQ4QEBeOidJNUkla2BtIVRuRmdibScoNz9uOShXUioiaCJPK0ReaXA7O1tHKSs5PVFDbi0rR1Y6JzRUR1MwWTFoLjpNdFxSUSE/SSs1byZpIVhvInRxLVBMZ1NKTjlLbF9bP1smLj4pa2Y2TltqSCFVQGBWVmJoWSRoNSxJTjRdRVY3azola2wncWhwSUpjMEE0S3AyWSIkXSpTT1dCaDRKZTFAXlRiR0NZI0FWOTsxUWBdPktfSjhJP2paKWFIREU9ZiFvKj1BITtpWk1mS0g0P1Y7Ty8+OW9GQkdqVFVjKFM2KT9Aby5scVhwcCFNbkcjMVsjIWE4Q2E4OW5RN04wQyJeYHAkcUBOVTRiT3Brb2dLIllWLnBscC9cOitKRDouSi5lPGlbY1xvUi1UJTpKbD49MmpuOV0sdClhUEFzZ09dQ0FILENHNEs9XFhTKWQ2IUczQmswQlBGJ2ciYkgoP0JJOT51NzBdL3UsJzNTWENnZVllIjtycCNtPENDIU4jTEM1P0g0TWduOFluZFBIQyo9UzI2Lk41ZFY4dTwyYF0xMmopW21Xam1MLTMlJVxRbjlaT29CbXNYKGcmO2wtLE1oTy1oNXVbY2dAamcoaC0lbWZGTlNNPEZMMS01UWZjZUBEWkQjOHU7QDk3aj5mMlRPYFo4KiEwNG9gbFxpPEMrTHU2JEFQOyU+OzhyLWlzMjlWZHVjNVgsXmRsc2lqIWNwTWZjXl9vKTI9J3VjKWNRbFVldWQnIWw9MVRZTyhcPDBgOFlhQC4kMjdbQE1iQTouVCo5R2FJTSIiRl5NMCFJNDZuInEkYzssLnJPYk1qIUthK2B0JEApPUYhSVAtNUFSUWJvKHJtNj8laiw1NUozST83aClCXEtvYW1zLjE4ZEFmVnQhREo+NEJWTUtabTshKygyVF03V09rQTRdPUYnZDFPMSRjP2BmYXU7UmIwPThMWjQxc2BhVDM9JUg7V1M5Wkg/QzIxIUNPNTxxNEpHP15yLDwzWT8mal0rKGJDMCJNPDU6P1A+WFZUZjVAXzAnWTxTOzNRMkkhLmJFcWxzLFQ2cmQ7KDAsR1M/QGZUaEtdNE5VVG4nTnIiXDpuUDxWRWwmJSU1MnJVMy9JSzBrSiMwTlpHIjBjUSgsX2Ymcy0yQWhXK2AhanI2LiJTVU1mMm8qIW80I2ZNMXNqaEBjcDFoWlEiSD1Gaz9vKSYrXGNzXic1P0s+W1hROUtnOzYkcE11UjZwMz9TPVchLmojXlNsbm5LY2VAL2NOSFhLNCNqJUJdOkBVXU5TSUlYbENRZlZZR1peUTcyYjhjJnNWKEE3czUpNlUrcVxJXjBPS2NVRTZMYHQ7I1pJS289T1A7TEAqOzorUWJnOGJzIjVNXk9OJ15FUT5WKyYqcGNXSnVIJkZocWQqXC84KWAySGdpYipIcTY1RTlvLlokR2AvOFI6N01gM1dZQjg/QHNsLlEnLmNiPWknJzhYKD8zczE0dDMyVF08bygvUVlra0hcLXFuP1ExWi9dOHNCTilTMEBOWkZWJ1QnciUiVSFccDsmaVZbT2RLbFxgX29BaW9eQT11biFrY19vaj5ALDZjRV8rVDU0c2JiWEtEWklCcUR0LHBEXW1JQCJkQjRnQlZZO1ZyXTxJODxDUHFCIj5SSFtMZEFMRT4uNEhTQT89TEU+K2NVbmFUKSY8N1svN2c1TlxZQmpDOnF0azEhXEVDaEdBXFJOSSJBQmFiRFw5JV4yQjQxKFhcOFkiJC5OLD1EJks4czw7SjwpKzJbWG5YQXM/Lk4rLideJ1sndGxtbVwvM0c8ZWlWZTRVOW0uLyolX1xSYlA8Qk9nZkxeN0JANllnPm9dS1w2MDZiISwzXClUbGAkZj4mQnVqaV5NKi1UNzQpUCQ6VCcmcyhvZ0ZSczkmaio+aC9eUEJIZ2tYTCgiNHBGczBUQF05ZkU1QDdzK2IuPVgvOUIlbmY8NltNXzg1VHVkLVk2WVBTYSNLLTpTPDxmXkUzKWlSXilFXG9ddSFAMm5YI1xCZ3NyZVpNW08sZHM1PTpNU2MsN0huYyY6OmNLU1tFKCs5VGxzaXJbYGZVT25GYCs5P3RFKkdXY15pZlZ1OTxucy9XZDpJcUBjXildUFk1a2VIWnFAaT8pXF5fLms2JSJUTm8waU8vN21jbmoyPW9WKVssa05JU0loO3JcJS8iMlltbD9DY2s8Zl9abmkhPjtJIW1EN0xkZ3RGJWtWI3FAKkwtYSZqIStfO3I4KGJrNEo9OzI1ZFBbQyloKFNROlkrKVxCZF48QDZlcG1pKEtQIkRsal4uTSZEVjBHdUhANCtPP11sP0g6KVZkOTQvWScsI0RLP09mPGxVT1BUcDlIaCM0YytwU3NhLTApMklANlYzY0JnL2koLHRbcD8jJUNeNT9fKl9tPClKP1RwJnVPOWIqc1hfRSUwXDpMWkVxU0JXS1lTRnNWMDZhSmRtZ2xoR1s2VnU0SWZwbVU3Y0s7TEFFJksobDxgYD5AazZmX10rXyQkZGxRQWohUlxVSGJpQDVTNDFAQDw4YjNQZG45OnNzWEgyPzk5P2pEY05bMUhBXl9tO1FVSF8xc15nRSxpaTJFUG4yZk9lTHFNPFYiLFk1dExcamJePiNMTmdOai1iXmxcYm8sX3JCaGhrPzlRU1coIkFGI1RAKE1oWWttKydqUzpMbTNrKi41Nz4iSXM5YzA5Lyc5bGdNYkNtcmkoQGFnaz5HMWNIRkJfKG4qcCMxJVdBdV9PJDJBXEZFYz8qQ0FBIz4lMFRhJSVMUFUvO3NObWpOdFFKYVQ0VEVoQVREKFJgbDhjS0tEWSckS0I6NiJGSU9fSGchZGNMPkZQVUc9XigwXSR1OjRHb1VSUFdTP2NTcWAnOlwwUTteQnBEaTdKVCMua0lLL2JSY3RGTWd1XlwrV1Y7ZkZmNCJZOER1JF9Fam81bywydShNZGZEWGJHb3NNTkhrLGlBdW9YbVYsRD48QmJAO1JpLVEhQEAkMi5wMDdQQkUuYTpVWiM7UjZOMSVhW1kvJE5RWihhXD8lImItJTFoaiR0Vjs/PlpQbzlBTW1UbGwkSHJOT05ucjdlYmwsQkwyXzQhW2QlUycvU1dnWkZjVzwzYUxAQDBvcjV1M1RCaDU7VGhbYjFINGpoajYwYkg8SCpyMUNfSz8jRVdyLjk1IVJKKzIsVlFaLW42OTRYY18tXzxaIWk+OSk7TUdcJ21uJlowZmZxL1BxKzphUVEjYDkjSjoiRU9rUWpbNjJfZVxqImAmKT88JGwpS21eWj5lZlZUQltSR1RqVk5eWUxIRVddazdGQTVdUG03WUFLYkBJZzlJaF5XbSliX0ZLOkQoMFc5RVpKQ2YqQTxhWHVLQStfJm9vQDotLCIxWT5ya24nQDU0TV1QRWYxXDBNY0hOYVlUYUtITW5nJTRoRSo3Ri89MmZkRjM6X0cvJ0MyOGFpO0JAYWFPJXEuY1s0Y1dQQlhdUDVCLU9MYDZmYGFMb0YySFZsPUgxTEREIlRTZ0suPy4tam5oZVdqbHMzRWIsWmFXOSguRnJYYVZYXiw8bnJCKlYwWzdNanAjWTNyZy1KQko7JTdcSnNKOTNoJ0ZdKWYxVHM2K1dEaFUkQT9iOEZUNmQ5JVhnWVZlcUwmR1xKRlpzNDttNVZ1ZEE9UkIrZS1FOWhJWyhrPiZkUSs9YCZFZEx1aEA4dW0jZUZQP0k4LCJtaz02PGo/b2tVPls9WHRMIVY/c2YxJn4+ZW5kc3RyZWFtCmVuZG9iagoxMCAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAxMDAxCj4+CnN0cmVhbQpHYXRtOWdKNktnJjpNbCtiVGo2U1tGXmJ0IVZNIVZKdVs3VjZJX28qZEZvWCckYmNhajliYj1VWmcva15jbWY8KE08cXEkUzJrWSZiZV4udTJwXHU2cGM4TG1iKDJXLjk3X0g2S0ByK2VjM0I9MyhAcSI6JVk5PjptdTBEVXFkSl01RS1qSTg7bGAxalkhVztIbVkiRCVpL2NRKTFYImBwWVZALDJdQV8kZG8vOSckMDxrQzAwXS9gMiVYIlglbi1xP2VcLSNbRmFecG1qLitMK1lNWzUoLG1IWEVbO09KYyo8RiNrNmJuPU5TYiNnJnBRKE1ILSYkc1U5ZjpLZ28lRU5WW0BlMyciI2dWO2M/UE4lVV11S0M4KllHMm5pIiNLTyJQJ247ITxiUlRyNz9hQUU0JnRlbDJaIztQcDc9Pk04cyFKUFw8bCZba0wrNiRMOzhfYW9JMSU7VjJzM0gmS1deMHUkKkllZ0tYZjglcGg3VGNVKENoXjQicW5xUkkkWiRjXTBDVVQ+TFFGSSdlXDxbXmJUR0YhQFw4PHBOLV11VEpLMT1VTyx1ZUhpNHU2aCwvPzhpLS8iK2JvVjw5UDhnLzxKNyZgbUxHKUEmIk8xclVBSTI2UU1ZY3JMb2EjXi5JXTU0MHBkbSpuPUBNTltyXz9UMGdFPi4mL11ibkVWdCRDUFJfKU9KIWpZbVxnMEY4UFFWZlxnN3QscVRwdCwlPGtVSFVtUlBBJXI4WnBvS1snMCQ9SnRrUmouSWQkSDZNZzhMKnQ3UDBYaWFWOUVlYzlSZVovLi5ybVB1Si1jWSsqXGRIQVYuUXMpaVJGWSolTk5ublZeV2pEVUE1Om1HMjkhRVkyakM5cD5WRV1Qb0NuQk1SNjtTW0UsZiNNMEY6OVRBTTIlJG1NP0NKZD8xY09APVRbSlsyJ144WkpNMSE8Ul4pJGw+KF5tP0lNTiZzXmAmb2I0UDx0UiY/Kyo1a0M5UWNCTkdAXUpdISlqLEYmPmxcRGoubSdCNjMpTSw8VExiLiJCMW4oLzxcQk8qKmg5RyslKnFoIiNnVV88RidZZ0xvWlQyVidQUEA3JEshTTAoPzYnR2E6bVJvQE1ROERBIiJpaTk2VFUjRykxaWtnbDM7anBmYz86QG5IOktDRTdTIyxPckEhc2NFaF0uI2NyZGlaMTZjTT9MRzoyRiIwUTJbPCpvKT1hYFpOMmAnMl9FXmxCLHU3PWUqNm5HIUJCTjgpWEQsXlZZb3RnQ007aTBbaD1TNj4zQy1OZmNwWjZXVHFlUTwkOVw0Kz1EN25qLUk4PzBgNGopI09QRlVRQDh+PmVuZHN0cmVhbQplbmRvYmoKeHJlZgowIDExCjAwMDAwMDAwMDAgNjU1MzUgZiAKMDAwMDAwMDA2MSAwMDAwMCBuIAowMDAwMDAwMTAyIDAwMDAwIG4gCjAwMDAwMDAyMDkgMDAwMDAgbiAKMDAwMDAwMDMyMSAwMDAwMCBuIAowMDAwMDAwNTE0IDAwMDAwIG4gCjAwMDAwMDA3MDggMDAwMDAgbiAKMDAwMDAwMDc3NiAwMDAwMCBuIAowMDAwMDAxMTMwIDAwMDAwIG4gCjAwMDAwMDExOTUgMDAwMDAgbiAKMDAwMDAwNDI2NSAwMDAwMCBuIAp0cmFpbGVyCjw8Ci9JRCAKWzxjNTQyYTA0ZjhiMWY2YTU5NjRhOWJiNmMxZjI2MGJhZD48YzU0MmEwNGY4YjFmNmE1OTY0YTliYjZjMWYyNjBiYWQ+XQolIFJlcG9ydExhYiBnZW5lcmF0ZWQgUERGIGRvY3VtZW50IC0tIGRpZ2VzdCAob3BlbnNvdXJjZSkKCi9JbmZvIDcgMCBSCi9Sb290IDYgMCBSCi9TaXplIDExCj4+CnN0YXJ0eHJlZgo1MzU4CiUlRU9GCg==",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for guide-level read-count quantifications from the focused WTC11-EC validation CRISPRi screen.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 08:34:32,135:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -d -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_validation_differential.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/document_validation_differential_dryrun_20260930_083432.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:34:32,803:iu_debug:	DRY RUN is enabled.
2026-09-30 08:34:32,803:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:34:32,803:iu_debug:	submission=True: In submission mode.
2026-09-30 08:34:34,185:iu_debug:	
IN post().
2026-09-30 08:34:34,213:iu_debug:	Validating the payload against the schema
2026-09-30 08:34:34,231:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_validation_screen_differential_element_quantifications_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_validation_screen_differential_element_quantifications_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_validation_screen_differential_element_quantifications_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOSAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDggMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL0NvbnRlbnRzIDEwIDAgUiAvTWVkaWFCb3ggWyAwIDAgNzkyIDYxMiBdIC9QYXJlbnQgOCAwIFIgL1Jlc291cmNlcyA8PAovRm9udCAxIDAgUiAvUHJvY1NldCBbIC9QREYgL1RleHQgL0ltYWdlQiAvSW1hZ2VDIC9JbWFnZUkgXQo+PiAvUm90YXRlIDAgL1RyYW5zIDw8Cgo+PiAKICAvVHlwZSAvUGFnZQo+PgplbmRvYmoKNiAwIG9iago8PAovUGFnZU1vZGUgL1VzZU5vbmUgL1BhZ2VzIDggMCBSIC9UeXBlIC9DYXRhbG9nCj4+CmVuZG9iago3IDAgb2JqCjw8Ci9BdXRob3IgKEVuZ3JlaXR6IExhYikgL0NyZWF0aW9uRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9DcmVhdG9yIChcKHVuc3BlY2lmaWVkXCkpIC9LZXl3b3JkcyAoKSAvTW9kRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9Qcm9kdWNlciAoUmVwb3J0TGFiIFBERiBMaWJyYXJ5IC0gXChvcGVuc291cmNlXCkpIAogIC9TdWJqZWN0IChJR1ZGIGZpbGUgZm9ybWF0IHNwZWNpZmljYXRpb24pIC9UaXRsZSAoRm9jdXNlZCB2YWxpZGF0aW9uIHNjcmVlbiBEaWZmZXJlbnRpYWwgZWxlbWVudCBxdWFudGlmaWNhdGlvbnMgZmlsZSBmb3JtYXQgc3BlY2lmaWNhdGlvbikgL1RyYXBwZWQgL0ZhbHNlCj4+CmVuZG9iago4IDAgb2JqCjw8Ci9Db3VudCAyIC9LaWRzIFsgNCAwIFIgNSAwIFIgXSAvVHlwZSAvUGFnZXMKPj4KZW5kb2JqCjkgMCBvYmoKPDwKL0ZpbHRlciBbIC9BU0NJSTg1RGVjb2RlIC9GbGF0ZURlY29kZSBdIC9MZW5ndGggMzU5Mgo+PgpzdHJlYW0KR2IhO2doL2g9LyZxQlhfWXJERT9HIlJNSUxJZkVcNTAkRkk5OlhML0g3dGk6Ol9aO1JqYkgwYUcrQGU9SWdeZSYuK29qa2tqaV07MWEwQmUwIVArWlorRipNWWpgLzhKJmFSYl5aUCVWUkswVSReJzxTYzAkTjRoZWJuZSJSMk0qKkFkKVFITmI9LCJycCdgPW1EQTdbMVtPZFI8ZStdP0s4bmFrT0chIlVrLGBFMm5FSSNjbHVUbSxjZG0yQlZRVDVoPzNtRUBpPG9dZV9FVDExaHVtT01RbEtGOydoTGpDcUUpP0pRPS5NXzJoLClNM21CKFNiOjhRZj9sTCNebFNPWGVabi1aWHElb0Z1YltoNSNxTXBXInJAUT9WOypvZ0IzcmIzcm1gImVmSEVOOzlpZClfKm5SYWFbMTotbFlvdFsnNUNpQEIxI2ZcTyREJzEhb1pnSTYycWxXTUhqV1NvYEpmISRNZ2I/PGtxXj4qPUBeczU1TVVMSyFdTUAnbz4xJyc+VUM6ZypGKnJZRTluWTEyV05nP05gWFonKXBudSRNUScwJipbM0dmM21ZTG80P1gtJT9IU1R1aG9dLSxcS191RzArRCkxPFZCOlBBMlsqRmdVSnBzXy5mOUFCWnNuRTNgM2VIQHBmUHA8PDhsJjxiRlRTUEJOUW4kJ2U6L2c3PylAb2U7a1JGJDxHM0dAPFpEMmxiTyRDSGN0VlQqNjNrI2h0P1xvTE1IJFNkMUpGLEU/Mlw3LUBDXW8lPkNgXyFNXkJpMClgTUZbcEE2QT9qPlEhU2syb1AkZjljYGJuKz0rblMrKCZJVDMiUUltLFsrRzFBcD4zVz80I0haJUpjLURdU2ppbm9LLGFSOlo3JDk6LE9eRl8vWl1BTUgjMEtrcilOTHUlXkBSIys7Jm46USxSOi9tclQ4P0NmYHNRcmVlVzZJdCtMXTIuVl9AI1hhUFBiZ0s5SlFrW0pQKStvcTtaI0svYDZaJjk0dS5BMTFpQjpyWVg5I3VJWkxVYj9Tc29iVlo3MDg+UXFbI2VZXTYvSl1yXXBYOD87TGk+SDFWRnJLJTRbQWlAXUkhKSRpMzEtVUJZU1txNVpQWzJIYSEyXEVMbl1rKkxgPD0tbmlDJEA7KGovbW5ZREQrXzdIQkIkMWhhWydbMmNyXS8iT2pkQ2BDVUUwPjJCQU1JSTlDbjg+SUBBLFp1cTEyVXJvT20/SkxfTjI5aVc/ZUplOiRxI0YtT2VkQzZzJGkhT3BmTWxFZ2YhT18nJDhiOFc5WnNNMkUzOidYamhqKkNmSix1cEtoNTZCRVNhWistQDBqXWZQTnFbU3JNbChDUFgoNTFsblkscGt1LkdVcVZlMEZqQltYWz82NUM7JjBBXDJcUysmIih1NiJIR1ppIVRmb1JLKjA2U0hBJ0NsSm5JX1AzZGRKWlxXQHVmLC8/JGBxWHNcay9KR0xtMGQpcGVARnIxVyhdYFd1OUpPNGInUSIqJTIiJj5CNVEjMVhLcHNHIkxwVzE5dCVlQ0w0UDRaX0FTbU5FO2VGWnNuJVZUIikoIVkvaGBPMWpKaS9VOFc3KGguYG4jZGcsbkpzMl5CcCNXQDhBOnVALDsvUl5gVTFuVytwUkZrRixrYzs8VFJKSDA5STc5cjEpaFQxTlBsdVJUb1c6YXIkL1s1Y2BJbCd0V1JPOVc5akNCaWZcT2k1JXBpJ3JpWy9LU01XR1MuTUNQWE46akNobWQwWGVqcik6SlwiODVUbFM8cF0wMnFxSDttVGs3Z3JrRmJfMTI9bkkhbkdrMlBKc0ZTIydxMmtUIi0iJ0FpSzRUZyNHcGU8P3JzOD5gY1JIJjVccy1HT2ola0sockxRVDQ5MkVIcU5jVE9mMkxTNWdpMF5wMVVlcTkyPSplYj1xJilEU0prRkskMFlOdF5AJipbPGwxc3RUJmAoJWA+MXN0VCVpPSk5Mk5bLE1gLSpabDhsbVNVVTsuZmAjN2hMYEBpV2EnQjQsUl1uJ19GKkhTZyMqT2tkV0cpJXVyUm83NjhcQF5AR3UiaHF0JmEsOS8obTQ3YW0pSzFfdSRJbUFzJFJyQmctMlxzYignNEE/RCtqbCEsajdEXTUkYikwTEg+ZCYxbydFTjtiblxrNys7Zj1UNHBFJlskSCMmVCdIYyFPTl5DXlFHUjAvKFNRKjZSSnA6L14tQHE2KTNMVm5zKlc0altGRWAiV2A3My0/WXNbV09hV1RBI1M5JkBeb2YrakNocFRvLTJxdD9CQ3NYZy87MyMjJkAmcTFCRmgjYzdvQUNhV2s6VVNzZGdJbHA8LVUqV21QO2c7STkzJk9nc00+S2YuJChMTy4tSGJiPUQuXGRxPSwmXkddZHA9NXQoO1oxKTc/cnA9RVwrXU1LKzlTaigkJmstJEtDViE3cSVjUE1YbUs5MWJKYGgkOyhBNTJsLWgtMiNcNCplRk10KT10Nz9dcSRGQmhOXFFMcFM8NT1MRFZJMCNDTj5qIV0rb1UtTWwkOjNFJ2VGVT5qRCRrQk9bcGJcKmk0KFwuKzQvSyZkRk9MOCxAKy5ZQUNKXFNSaCwsNj9PO046S21WP0A2LC5hQUxyJGtSV2QrRmJIUnQsOjcpNipIZyouXyQrLE5zWT5vZWIlTG8ubDVaR2BzU0Y7Tz9vXVI7cyVSMSFEdF07KGUsWjtuXWtCP2VTck1XbD1HJmZVZVgydXIoZ0xjUHQscT5xLldcTlknbGJrNjcyVzZ1cUMwUlRqcC5cVDN0YkZmZmwoQ0BDIyg9XkZoNGpEcGBgPmEuIUR1IyxGXF1BdVo5dGI3VSJVbHRhMmltMSRAJlZEPzpMaDNebzQ9PT0mWkBBRTVRcD9bZk0tb0dkQVNeM2EjPEgxPXA3Wkw7L15EZjVWPEw2TjBTQiRFTFpmZGVvcTFhcW1qTFZwbFdxTzohMSRTcG1vM2goNiVEIUNhITJiMk0qO0kzWGdUZiNXTjxHLU87YUtKZ1pCTjZYRk0+UFZLXDI5LWsjKGNxcGtXTCchRTU7cjxqPz4iYWJrUVc1PjFWYilOUytpZVlvP1QwRTFsSGFTJGouNlVTcVEuOzVMJW9PSkRYJWU1MEJHSFBsM2ZRYVUmJ0EkN109JiIybiwsdCVRXy5FZ19kPVBWYmtOT05caiZSJk9Fb24iS3E6TkV1KSlEaj5GMidlJVdfZVdOYERLXSMtIWJRWS0oSGw4MUw8OlhXYnRPM2teQz9jJXFFZUpFPGE+V208LUFuTnMnaVdUSXMwcEdjZSI8W2FcKkRtRkhDc3IkQTgkPmhJJC06OFAsN1gmOTAlXU4rIUIxKlY1OzsuNHAlVzNxJkklI3NtQ085REttLnNyKlN0TShwKWBrJiovXG9mI3RMXnE/ZFIyNkdRKFpHVl1tKVFzKyM0TmxSJDBFKCIyYlhHPVdPNj9jNk1TVzlxNDJacFBHLzBCI0k/bE8sYGxuYlN0NGlqdXUlJ01aOzAya3NVZy5OYVgkPnVqJ2orUClWRE89ajY3IkVSa2tsaHU4I2o7Mm82aWlNKj9mJEcvS0AhQzE8MlI1XTcpO1JjXFAsXD1tSlNFO3E6J1FfUy5UQmU6b1ZWNDVqIUtCSlNeJEZWQXJCWTZXbVlqQE1eWWAjW2JcL247OXNANyZQJ0o0J1VOQyFkQVZRdEk7alUjbDtaRkNiYSMwNSc2cWlBNUA0bCMjJ2MtSGUpWFYiPWQuRDVCamU8ZTNQSjYjW25maUc5RyI5ZSsyKitPXTRBNGskWFs4UnFmSzR0NkM6NWY1XWo9YWApPDc1VC9CZFMnbTAsOispJVNGImZOI1tKbkB1NG9wdD9zYE9LV21rdE5qKURKcS4wKFtJPmxiLiNYVWIwNGMzaGpVaylLWFdoNTJpbEM9bltobSNETWRRckplRkErTlooQVcxO21DdCJxUzRESDp0XChFUjVFaT51VyhrMiU9MU9lO11CQ2w6dCpIRls7NSspPltna29rdTQ3PUQpMiVpUzZzcGxhJ1lbcT82MmErS2B1JDc7VGVZYVVpSXM2UENARzpcWFNWIWxVODtULmRLSjxgbz1yZGtuTDw6SClULFdiWD47ZWY3bTU9VXE+clFiSlE2ZCYjJnBOTk5yZVVJKDsybFolJTpAX2VrLDVNcC9razhZc2heND1oNywwcmlSM0NMPWRRRGNGTW9OPzcuWF1OWWVZLmZGPTVjMms2JUMtWGhVQjZpSWxaJiNnVVByJltRKGYvLlZcZ2BNUlhuc0ZMLCMwMDZBLWNVXCRUOlpQaTM+VTtWckRqa2VyXm9KLSxpZ0lKIyNLJ2xET141aV8xKDBOYWMnO1ZhVS9FJGNvQ1xSSi0xNnFXKyk2WGZXKTo7ZTVzKiw+JyMicU91Lz5tLnUiMDRlL1RIcV07R2ZdJm5bcSFTWjo3Km9CcmF1NEdVYChhTz9mOig3WEQ2Xl0oQkM5TkdBU1xjWklKQF1wWWwvOzlyU2knLiQvUEBMW1xgJmpmWFtvOTxfakFkYzk0bi9IQz82bW1WL0VHUWAjYzxqaDxvYVshSi9oO1hTN25lXkhwJ0VfZEBORDZhKXVWWj9GYExuRiYzODpWKlZyNTJxUElsSyI8cFA0SiFyaWY0KydGZjcjZEhDdTJtOXVeJjRYZWVyamhTUFpLJk5aYV80bGdBVXU3ITNqbmQpc2lWRi51LWFvJ3Q8JCxhNFA+cT87VXRZT0RyO2hGcmc3YTRKTjxEczxELyZSXSFIQGE+MCdFJCt+PmVuZHN0cmVhbQplbmRvYmoKMTAgMCBvYmoKPDwKL0ZpbHRlciBbIC9BU0NJSTg1RGVjb2RlIC9GbGF0ZURlY29kZSBdIC9MZW5ndGggMjY4Mgo+PgpzdHJlYW0KR2F1MEY/WlhzWyUiSFUkXmdqZ1JBUCYxPUgrPyk2OF5tWSZcaDxVQVdIVUArVTE/LnBrdShFImwwWkEhLFVIVDJBaSNfLjk6TFNbI185XSQhIiRbSihWYCoyNFhnRDwwNEIxZUpganArIiMiVHA9TUJESGJIb04tOkAnUnBEQEI5JiRpTk1FNUpvSlMiMVdjK0wvYTRIKitoRDZpbz9VSThyMT88MU4sYmVyLWQtT05ZTEJcQkJvW1BkaWFKJUR1aV90aWUuT2w8dVxlUTNSLTIpO284PkxTa0k1azBDWSInP2luKldyVD1BPzxgJGtqLzRnPCRMKmVIXSQpcz9ONUpaaCcvOlVYKSJxNzgwS1RgKV5jT1Y7IWkiKFlEX0BxLkNJZEZ1YUtUWTtpTTNLO0toW0FWR20qJkhWNWluYWBlMmhZPm4rZktxYm5USU0kR1pnRW03NjY3XjA3NyJscFVrQGoiMD9fZ2xLNlJgOXFJI1guUVhdJGJLXGFKMHFbQ05hKnQwRzBdJWRrcS9ELWhqSm1TQi9IMnVqL2hGOWUlIi0yInB1OTVYcDFgWzZhRmA7Km5iVHQxb14xRWNbNFNpci9dUixYL0BtP2JnQV1wTDdWanNyaGlHaGghdD9pcW0yXCgiaW1qXiZrLGFpXiYxJW49OTZVQyVHZj5vX0ojRFxfdT5zTkRPQDZhPzRtMVxgUV5uRHBrKVFdaiNOM3EmPERTTC9PSUs+N2FPJDsoaFk8TyE5I0BQSmxmPypON1pVL3BFUi0rI2llWWYqKCJabzBKRVBrJ0dxOF9IPzBJc2BqPWVhczBwJl4pbTcwZ2E5azRKKW03MmVKImxoIkRAZ3JERDFpJ1tqJnVEJXEzR0BjNSo5cTNBNSx1TjcjbHUhNVZUO0hjWktnJEY6QiJYIzc7MXRCVk9GLiFANV5jUTBYbTdrUm0ibmU4I0RHKzM5MF1GXjhKbzw2ZGNgN2RBX2YwTGNnM05iPyQ7ZVw+Ki5qRktHUysiPzI8RSNkM1wlYGUrUFhgLDkrXHFMW0tULi5dUUhHbVllbTNRRzReclElNV5dUSMocVhLP1lmWVReUzBHN1RwIl5DImkucjkmTD4oN0wmYFhAJihdI00oY1VNRHUxLyVXYSxoTitmS1EiPG1gRk9pTTkhaTQrbDAlOXJKKlgjWFhpZlY3XkQuKkohR2gqNDg9aD1WW2ZSZHFQSi88Y1E0RjkwMmBxLVdqQ2NxJV4rTm07UWQ3TUo3TkJMLDZLX1ZSXiQiWW9ySGE4KCM6OzFMMGZtKSdUQVUmbFQpcEg7cF8jI1wtJiVsUlpUVCdhcSNNT0ZaPGhQRmJdY1VxMSssP1FON0VCIVpQV1teTEVccVA4Vi1KaGBhc0FGZCNOZlxTMD51PWEhSHBTYmtBQTJlR2gyPlFFXl86MSUkVkV1YiY/MEU8X3RYbW9lXU02ZmFfYE5HNXNIXmBAWl8jNTpKTCshX0g9Sy9AL2xuT25lTmJXR0ZdbD9MW3NsU3UpImllXUk4WyxgZXNlTFluV0paJCxjL28vSFVJbC0sIjxeN09XVHEkWyo8V09BX21xcCc8QGY1NWwrQV9ZKlksTHVHUj5GZj4kZTRtUSQ6K3IlYlphWmNbY0grNT5yUnRZVzQuKDMrUU1aVVoqVU5OQGJbO2o+bzd1WiFCU1JdO1NBaSxnMFZdZUJuI05UaV01QGc3YjZ1WVxcIXQ0OkFUWDlKO2Q5akUoI0BYMzxxTiVoZEpfZzs3MVpzRWg2bFdOb1VvRjlTXzdtTCYlV00zS0MsNSYjIi1rTV4iaWxuRCU1VEVEKVBXVTdXUSEzKUN1ZDdcXDUsIjlBMkZFRGw4J0tZSEpWLXFFdUo0STwuMmheN25EOWwudGdIQGMnSENCIkkuZmY1V2NraS5fdGc7VT5ZUFcjJTxhOWxTNE8zMDcndCVcREV0L1k+XSw0UzBwbDpwRCtDP1NnbmQiSXN1NiwxRiVeI1AyR14rajxwWztZPzwrPEVONVFrQ3UlZVJMKVhzLzBQMGYoJGMwZ0dcYjxIOmY6VSVpRVZgZGdHVWlhZUg7cU1CSykxOiNtczlZTWBmbCJYNWBWIzhZJ24xUiEiXVlwNjtaYEc/JUBKWGtVOD4lKXErXDZsVk0/I3ApQ1JQXThlcCphKjpyQiNtWic7UHVTKGU4JDQxRD5KNy9POm1hJScnJ0twLj9xMWRrWCc8YEJrdCIoLiY2cHRpMVdHNyxzWj5uWGpxRnVVazczOlkodSJRZXFoMUFxPkBOSFxwKml1XEBLPiVWKEsiP25SUlJvQnIwKGtTISRDMEJTKGg1PU9IdS07QUlldWBtWUhzcWFxQlgnYlxwYF0lbm4zViFoV3VwOWdKSU47KjNzaV0/S09SVWM6dTUpdC5SXGg2S08oXCVGTWpqYyIqakhPbElmK3RWYFxbSk1AQExlSlhSbU1QcUxHIUpTLj8+byNgIShQWG9UO0I1YileKT5PbmJlSUBKJ2NUUVdJNXUvQ3RfMHBeY2gmI1JGRj1aLnFyVCtsRzNUP1lOMUZealclRT1BNjtlKTNKZic7MVZXL29zU25rSmJMITQzKkBfVURgbFMoVScxUSFROCNIcWJOLjYrKUVZdCsscGNlMklFVmUtOSlIRG0lImBuWC1aKDstR1hDL0xLKmVHUjBvIS5JclNkKDJhQVIhTkFfLXAvY0wvcmhaUUdIaSFzYyFWOk8lPEAzaGQtZEVzK1VgL0xnUlRgTXE6TjY+UU5KcDE7ZlNCYC4mczwsRGIhIT9vUjBhclhbPVRRXnQiVTtTajBuRVBDMSovOlZcMWZKQ1NROGFacXMiXnMkUmAsQVhtJShFbkJpcz1jPDgoND5VQ28kNTVwYlteKjhQLyJHKiRccnFMbCIhbis6LzJpZCRvUldhXCE/N2xIU3AjZVBXQ290U2BbQDw8LClMNWxSbEsoNDJObiphJ2xUYyFTUV9WPHRtV2FyK2BUdHJvdGhycGxVTkBOTyhCc2FFUD1oUjtjUkYsPT1bJVshXHJdY1tQKTloaE5LcSE6Ii8rJG5KSF06Py8lO01WPV1qPUBaakQtb1gjayFZa0c8Nl9EciQ5Q0o5bFcqUz8nJzghVSotS2cxZz1hPVFHKj9tPl9BO3VWWmExMT9ZRmFFMyEwcldLVyRsNWYvViwzWCpMUWdpdEs3ZVZgKlVVWTItPC9rIWxYU0FOWmVcSlNALjNoKiUiT21dXWczTW5QPjBMWl0wJmRMOSVmZWNuZEVwI1pLRER0VWAvWVBLSV8nVk5tW1ArVz5hbWJFJzNfKC1BdSQuYTc0V0lDOjJGMlc3LkVZOzVFKCIjRUxuTDtYJ1o/M1otYCtrYWhMTTgmdUcwOTBTVlRLQG5WYjY9VUV0Imdca2ksNkdYTS1pUjpTVnQ5WC1gJDUkUC5gaUk4S0FZPVokSmUnOlFOayo6STNLRydsOmtyMW4pInFKRSZAOXRiOExFTiFyc3AnP3BzQT42KTs6KiEqcWwiNCsuX0c5TnNPNmlINGMvLlxkJEd0SkdMSDFIJ34+ZW5kc3RyZWFtCmVuZG9iagp4cmVmCjAgMTEKMDAwMDAwMDAwMCA2NTUzNSBmIAowMDAwMDAwMDYxIDAwMDAwIG4gCjAwMDAwMDAxMDIgMDAwMDAgbiAKMDAwMDAwMDIwOSAwMDAwMCBuIAowMDAwMDAwMzIxIDAwMDAwIG4gCjAwMDAwMDA1MTQgMDAwMDAgbiAKMDAwMDAwMDcwOCAwMDAwMCBuIAowMDAwMDAwNzc2IDAwMDAwIG4gCjAwMDAwMDExNDUgMDAwMDAgbiAKMDAwMDAwMTIxMCAwMDAwMCBuIAowMDAwMDA0ODkzIDAwMDAwIG4gCnRyYWlsZXIKPDwKL0lEIApbPDcyOGUxMjMxNmRmNThmNzA5MWY0ZGIwY2FhYzJmZTYwPjw3MjhlMTIzMTZkZjU4ZjcwOTFmNGRiMGNhYWMyZmU2MD5dCiUgUmVwb3J0TGFiIGdlbmVyYXRlZCBQREYgZG9jdW1lbnQgLS0gZGlnZXN0IChvcGVuc291cmNlKQoKL0luZm8gNyAwIFIKL1Jvb3QgNiAwIFIKL1NpemUgMTEKPj4Kc3RhcnR4cmVmCjc2NjcKJSVFT0YK",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for promoter-element MAGeCK differential quantifications from the focused WTC11-EC validation CRISPRi screen.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 08:34:34,231:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -d -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_genome_wide_guide.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/document_genome_wide_guide_dryrun_20260930_083434.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:34:34,854:iu_debug:	DRY RUN is enabled.
2026-09-30 08:34:34,854:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:34:34,854:iu_debug:	submission=True: In submission mode.
2026-09-30 08:34:36,094:iu_debug:	
IN post().
2026-09-30 08:34:36,121:iu_debug:	Validating the payload against the schema
2026-09-30 08:34:36,139:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_genome_wide_screen_guide_quantifications_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_genome_wide_screen_guide_quantifications_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_genome_wide_screen_guide_quantifications_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOSAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDggMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL0NvbnRlbnRzIDEwIDAgUiAvTWVkaWFCb3ggWyAwIDAgNzkyIDYxMiBdIC9QYXJlbnQgOCAwIFIgL1Jlc291cmNlcyA8PAovRm9udCAxIDAgUiAvUHJvY1NldCBbIC9QREYgL1RleHQgL0ltYWdlQiAvSW1hZ2VDIC9JbWFnZUkgXQo+PiAvUm90YXRlIDAgL1RyYW5zIDw8Cgo+PiAKICAvVHlwZSAvUGFnZQo+PgplbmRvYmoKNiAwIG9iago8PAovUGFnZU1vZGUgL1VzZU5vbmUgL1BhZ2VzIDggMCBSIC9UeXBlIC9DYXRhbG9nCj4+CmVuZG9iago3IDAgb2JqCjw8Ci9BdXRob3IgKEVuZ3JlaXR6IExhYikgL0NyZWF0aW9uRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9DcmVhdG9yIChcKHVuc3BlY2lmaWVkXCkpIC9LZXl3b3JkcyAoKSAvTW9kRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9Qcm9kdWNlciAoUmVwb3J0TGFiIFBERiBMaWJyYXJ5IC0gXChvcGVuc291cmNlXCkpIAogIC9TdWJqZWN0IChJR1ZGIGZpbGUgZm9ybWF0IHNwZWNpZmljYXRpb24pIC9UaXRsZSAoR2Vub21lLXdpZGUgc2NyZWVuIEd1aWRlIHF1YW50aWZpY2F0aW9ucyBmaWxlIGZvcm1hdCBzcGVjaWZpY2F0aW9uKSAvVHJhcHBlZCAvRmFsc2UKPj4KZW5kb2JqCjggMCBvYmoKPDwKL0NvdW50IDIgL0tpZHMgWyA0IDAgUiA1IDAgUiBdIC9UeXBlIC9QYWdlcwo+PgplbmRvYmoKOSAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAyODA1Cj4+CnN0cmVhbQpHYXVITmgsRSchJnFCKUprYzknaFtaSmpwXjtFbCU1Vj8raEAyITtKamYwYF4qZlJFJEUoMlYmZUdccGg0WDE4QVlpSDVaLG1EbDBCRFFDbWUnbS0pSlYzdWpxRTwiJnJdLFBpRXNAVEBqOUN1YUhTRSFcWmNVJUUiJDl0M2knXnNDWTBVU2c1USoiaSpkbWQuSkAhZjMiaDVObzMzSilcTFg7K0xUR18+MW9xNV81JyZUV3RDTzZMaGNPMWFYIzI9QGskTzFXOyRXczBvQF1ZajppVUVwTjg8TUk2RTAzaTFKa21KPDVEUF5iJzYkbW00SzUjWE1IUDdoY0I6ND9NZ3JQJzQldClIPVpxTTUpa3VFX09sb2BNJkFXVVVJbV0sVW5ZRUFXTCQwWk0iNjdUOTpmUzhqQD4nKiZAMipPV20hclBtSSpSO15BMkFHT0E7NEUhOjVJdVEqIkcyalU8Sj4hOj9NXFlSMXVCPyM5UiMocVFXWU1caGcncSdzakw1PSkncCMnX2IhR01TJGVEcVVTW0huLllRPm9xKlx1ZzpnVVE1dG0/TG1vaWQuZnJjIls3VG9BbGVbZnE1XlE3O0phOWojPmZhZG9yWSwjL0g7VHNcUjRLdEAuRV1qOCNtRmVRbEMvTkw0b2sjWjBfLiknYi1iP3MsXC8sQlM0Ni9GKjskOXRhcTpoOSddPVNNbVEtZlVNSUU2NyIzXC1OdVAlIm9kT0QwbFMmKTBwLC5cMChRIktwcUJDZVQ6RG8wKEVbJlYxJ0pNWnBVXVUrOl5SWWlCPmJAWSw4QUw+O09fRE9IWVNJWV9zMilYKT0oZFIibk5OaW8uWWlTRD1lZWIlNmNObSFaSTMiTVVtS1w7bnAhaD82TS5UOjVxUDM8UV5EZyIwbUVMNjtSUWdmLjxaJGUpK05CaipuOCRRdVk/XEklalI6Yy5cXU8qOUknSk11QmNCQzpkLWx0MWM+aFZOOG1WPWw+ZkQ5U3IvISsrVkpfW0UvNCskS09BP3Q3SkpYSS4zZDFMQWc8K0ZMSiRpMktVOG0lZ0VjMHNxLl0ja1kyOXFuaWtSIWA7YzNyPG0wWT1mKEkmVT40Z2JbIzNsYztscXBTXDxmZS9vL1dMbzEmI2tIVTllajgpIThqMlcwc1FAXEFqaD9jazctTVI0R3BMOk1pXEtrVEtDUWNvLFJFPEtnXypcYlJyKGFZK2Q/LStwIWo5R1N0QSErbz03XllLMzo8NnAwPHIvNzViLytlZUVaZEY+aSEpKHVQQlYyMk9nMm5WRkxJPkY2IkVRTE01OCJRK2FrMTg/Zms2KE8sbDgsNTpUTGhHI2FWPjp1QSxiMXVcYVZqXEMuO3JxLTNHRENVPGdbY2JvOU5BTUlHRC9zTV86bEc1SSI+RjpmRi0zO0VNUzs5KiZpIl9gUC46ZV8kIWB1ZlFLIStfUiZAb3QzZT43O10hUSdcXjZFYGprRSZaaTcuTHVvIkQwU1FTYlFcPy1UVSFCWyg8T0pdT2xEdUlTSGpNMjZ1WTUlLG0iNmVLJ01mR0ZHNTQxUC5Wbi1tRy5VTW08ZW5EWyM4SCNuInNZR2xjb3JZMVUqTiREKT4yLS9vLSI+TGtjTlNLX2FwMjdfOVRQZkFzdChbLGpxSFtiZWlwYktCQnJxa2ZkP2lMPTEvXWw9aGFXbXNmQmxaOGw/TV04NXMyTF0iI0JncUUvdUZPYkFcOSFdWXBjRFE5R0soRiFqP3M9PVBpK0pIOjRUJVUmPDhkbWxeRWtHcXVSb25NMVdHVjArVVBsPWMtT0BscCt0W0JxcC85PVE6STV0ZCVTXVQ/PFMpUlBgWl4xZiQyY3IqIU1TMylUSFprVWE3R1hrVmA6LTInQDQ2Oi5LLjZARFFIaEdqX110XzdGT2g1TSw6WmNeKjxhTVFTbmZISUZcKURyRGhRRW1JS2lCKVNbNVpjWFdaS1dsWyRlTVdBYTNRNTk2N3JVcFpNPmw0MnJQO0xSPyU+NWBFPjdYNF0sPz8sYXFXMDYkWEtMQUdXMEJIXEpIPWFnPHRjSkZmaDxLJU4qMVc+dFxvS1B0V1ddUEZHRm1bUCYiYTpFMl87UCY7clkpQTZHQjY4Mkdla28/SWg9bkcuQC8vaCFVTWtXZk5eW18+QEA7c1EwVFNSMHBUKi1FX2FTSk1DQmkvIVtwRlFDQygpaE9yIWBZUz08PWlLVmZYXW1hUD5gYissZlxwVm00KyRNTmJaKWhYQThzLEdbdUVmRERlaEYwRisrRVEsa0ZGKFhhW0UsWHEuT2hwW1YkLWRyNjpvOFM4Sk9QNkc+NjZyNHAhMFNjPWthSWs1SGp1cV1YK1BzIz9hc0gmai9sbzd1WT9wcCQ9IWBIbV1nS3JDXTBKLUNucyFHJztPQT0wbC5qWGEtKV9ITmwxL0peP0BbOypWNHJDZ2FmYixIXTdZbyIrIjw8WzBcSGA9VktuSlw4RDlSVystWF1KXHU6JnRhcmZlJXAzXGw6NUU8SUopb1spTjcxWS8+OXFtLDsyS2woVCJhWkhHcSIzSGtWSmsrIzU8MGNaXFNgbk5nMiMjXi1aYVE8UURKU2BWWSdaZkY8MD9fM15gJ09QdCM9QyR0bCZvPHFXKV02PGkmJyswUDo3OjRXYmVmI0kkb2giSSdNJipOcVBzM0FeOitPKixhMzxNWVdoRDFEX19UMFFaKkchJko5JDJvJlxuRTtENkY+KmAwZ2syJlc7SlhIJWdgJWEnTXIzN0FWQldOJlwuZiFqPm1pMCVAWkpiQ1c3QlpZJFdrLGAuVzduJCZkcl1ycklgIUsicEtNaDhIQCImSz48UGJBVyRYYyY7MWkxXmBwYFk1Kl1sOVtdXlBcUF1bcktJLSUmP19LPD9sZnAlcExOZjlhPkBrI3NlRGU1TUN0Tm4tPmRNNG5TJkpMOVsmKmctJmVDJWdiZzR0M0FHOm0nNys3c19FTGAwbUNMXy5PR3QibGwvQDRrYDAoJUhoRlBYQSVEWDZEbiJAV1hmNVFQbW5dSFhbPFpRWi5AbUUsZlY3J11NTlItb1pRYj9wMktkSyg1QF9AK209dENhb0pjbVplazReNC4qPiI+aCJWcVdwZVtcZW5ZJl9Gayg5WDwhQC1mNjw0ITdoUSUnTVtEJWErRWZUQVJUOU5oSEFfa0xWaVJAJVkrVy4/PWFQV29AJlJeMzBVYnArT1k1WTQ+Xys2bixPMlsyYW9PQ1cwTyovNztSUSEwUFskOkkxWmYlcy9APWdcRyY/Z2UyPmMhQD4jOWg4azVyOWsvMT5FOWU3KzVNcHFddURLLkhWa3RgYnUoSmNkXzQyYERCUVFvXXBTR1gsRUYwLi1YYCIyUzpcJF0nZzEkWnE7KklLSzktNDs2UF0jZENxKTNFaStfQWRDYEROOiFwYl5uRyhxRWMoV14yZ0dDRW5Jb148ZTlMYjBSWkNeKC8paF0hUWBdTC5oNE9vTEJyQmBjNDlBME9hQTtbWjYkU20yKzxqNFUpI2giMjNHNkJJVz4iPlVlRkxRTG1wQGgyVmNkS2YjWEhbYjVtLnRgdU1xRjpVI2RjPk82c0QsWExpdS5ZODkvXToxLCtTX0N0IkEsIl1gTnNoTGwyXDJpOjZlSU9lWXVhOlpbaCtsTypvWEEwLTM3NXFOcjI4SVhlWCpITW9HVWxubUttUWF1bGghUCNmaUIoJnJXZEpRQEIwfj5lbmRzdHJlYW0KZW5kb2JqCjEwIDAgb2JqCjw8Ci9GaWx0ZXIgWyAvQVNDSUk4NURlY29kZSAvRmxhdGVEZWNvZGUgXSAvTGVuZ3RoIDgxMQo+PgpzdHJlYW0KR2F0biM/I1EyZCdTYylKLydfXWEuNXA9WDo6b19RVk9IKC0/Jk5rYC41RztIL0tZanQ8cD9nRjBwIjUrMmw1ZjY2bTNEPnBPMT4mQGczLWQ9VDk+LCErLydjOSpXVW4wRWwmJSlvPSw4I1VkQTdgK09YOnBjczVNMGByQXVOcGlHcWdnIVRuZiZOUC1iSlpVbGYoYCduazYvXkRtZz9ScU41RDE7b2xtXlsqXnNyIzFSbkgkQSlJPz4tSVctN0teNDxzUWtIODQ8UnFHVFQmPkpyMFM/UXE1RjdKUzY5bl41LSEvSHI3NF9AVEYhXi87IiwoazBJYWhcTkp0RCh0WmhMUEhicWZDX2hXLGQybl1rSjQzZl0xI2tvPyk0KiNESSg2Ty8vKlQnYTA5MlNkYlBGRDFmX3FDI1kqbl5KXWFEXXUwLlpvZyRpRC8kUSEpLCtDMTFQWDRQO01YUSdeTTNgJzVTMD1eQFwyMWU6TFVjZ0cjOVwnVypIY01daVk+U0pRRHJyVTVsaWA8IWpcQ2FycVlrPGJlclZuRzF0WE9zQjxFMyI+WnUtQUM0VichRCJbVVshOzkjSG1XcG1dc0pJZS5ZamJuKjVFPilmP1MuM0pNQ1Q7YWhibVhSWEstaTZjQWFaSFtFWm8tXkwwX1IsMiUkI15JUCxhLVxhJWhPU2JHOmlOWWolViFsT24+KTYuQj9UVydjVUR0QTZARWYrTkw9XEMwYDU7ODQ2X1A3TzVVYmhZLm08TlsvKFxoRFVNcmdBKkBQKWotIkk8M14rdD1Mb2lhIzBVa3RMQ1dKKUpGPjVXY2AnMDpkJWJOIVdPR0Q7a0clJE1PPVJfVmNSQCVXRGcpZU1sXVVPKVpjMmBmOjFhLVwtOFs0O0xhTSw0W1M4am8lRm0mXCRVTzQ3SzJeaUxdWEM+X102U1ZUcGJeNExOVisnKU1pbHJSZ2BfXz5nZnU0Yz43SElqXCZRXXVpazVUc3F0WyhVSXRQPTJsZ0ZSQ2xxYEskKD4/PzhAWSheUi06ZkVzPy4wRHNtLmQvOWInUV44PjVBOjMtQFoiX0JfOkV+PmVuZHN0cmVhbQplbmRvYmoKeHJlZgowIDExCjAwMDAwMDAwMDAgNjU1MzUgZiAKMDAwMDAwMDA2MSAwMDAwMCBuIAowMDAwMDAwMTAyIDAwMDAwIG4gCjAwMDAwMDAyMDkgMDAwMDAgbiAKMDAwMDAwMDMyMSAwMDAwMCBuIAowMDAwMDAwNTE0IDAwMDAwIG4gCjAwMDAwMDA3MDggMDAwMDAgbiAKMDAwMDAwMDc3NiAwMDAwMCBuIAowMDAwMDAxMTIzIDAwMDAwIG4gCjAwMDAwMDExODggMDAwMDAgbiAKMDAwMDAwNDA4NCAwMDAwMCBuIAp0cmFpbGVyCjw8Ci9JRCAKWzwwZWM0MzZjNjllYmNhZDFkMmM3ZGRhMjI1MDk5N2MxNj48MGVjNDM2YzY5ZWJjYWQxZDJjN2RkYTIyNTA5OTdjMTY+XQolIFJlcG9ydExhYiBnZW5lcmF0ZWQgUERGIGRvY3VtZW50IC0tIGRpZ2VzdCAob3BlbnNvdXJjZSkKCi9JbmZvIDcgMCBSCi9Sb290IDYgMCBSCi9TaXplIDExCj4+CnN0YXJ0eHJlZgo0OTg2CiUlRU9GCg==",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for guide-level read-count quantifications from the genome-wide WTC11-EC CRISPRi screen.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 08:34:36,139:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -d -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_genome_wide_differential.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/document_genome_wide_differential_dryrun_20260930_083436.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:34:36,715:iu_debug:	DRY RUN is enabled.
2026-09-30 08:34:36,715:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:34:36,715:iu_debug:	submission=True: In submission mode.
2026-09-30 08:34:38,048:iu_debug:	
IN post().
2026-09-30 08:34:38,089:iu_debug:	Validating the payload against the schema
2026-09-30 08:34:38,107:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_genome_wide_screen_differential_element_quantifications_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_genome_wide_screen_differential_element_quantifications_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_genome_wide_screen_differential_element_quantifications_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOSAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDggMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL0NvbnRlbnRzIDEwIDAgUiAvTWVkaWFCb3ggWyAwIDAgNzkyIDYxMiBdIC9QYXJlbnQgOCAwIFIgL1Jlc291cmNlcyA8PAovRm9udCAxIDAgUiAvUHJvY1NldCBbIC9QREYgL1RleHQgL0ltYWdlQiAvSW1hZ2VDIC9JbWFnZUkgXQo+PiAvUm90YXRlIDAgL1RyYW5zIDw8Cgo+PiAKICAvVHlwZSAvUGFnZQo+PgplbmRvYmoKNiAwIG9iago8PAovUGFnZU1vZGUgL1VzZU5vbmUgL1BhZ2VzIDggMCBSIC9UeXBlIC9DYXRhbG9nCj4+CmVuZG9iago3IDAgb2JqCjw8Ci9BdXRob3IgKEVuZ3JlaXR6IExhYikgL0NyZWF0aW9uRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9DcmVhdG9yIChcKHVuc3BlY2lmaWVkXCkpIC9LZXl3b3JkcyAoKSAvTW9kRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9Qcm9kdWNlciAoUmVwb3J0TGFiIFBERiBMaWJyYXJ5IC0gXChvcGVuc291cmNlXCkpIAogIC9TdWJqZWN0IChJR1ZGIGZpbGUgZm9ybWF0IHNwZWNpZmljYXRpb24pIC9UaXRsZSAoR2Vub21lLXdpZGUgc2NyZWVuIERpZmZlcmVudGlhbCBlbGVtZW50IHF1YW50aWZpY2F0aW9ucyBmaWxlIGZvcm1hdCBzcGVjaWZpY2F0aW9uKSAvVHJhcHBlZCAvRmFsc2UKPj4KZW5kb2JqCjggMCBvYmoKPDwKL0NvdW50IDIgL0tpZHMgWyA0IDAgUiA1IDAgUiBdIC9UeXBlIC9QYWdlcwo+PgplbmRvYmoKOSAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAzMTMxCj4+CnN0cmVhbQpHYXUwRz8jdU9HJnEwTVhrWmRyRGgtKU1sUyt0K1dOY1REX2lsc05qSFhCWi5LYS8mNmJKNDdKRVU7QyRJcDNPbkNzUyZaWlVsVyQ2SzVAb0Y8a3VyUnQ2OUNSSl9RPiFyMkt1OitkRDwpZU1gTlkzRWdpNGhSPStJNkI+UShBVXU6Qic9Sl9qUGlyZyR1SGZgaF02LitdOyNFSSlfWiNULUEzKjo/PWBuTC4+TFxoL1Q7aE9fXUlCTj4+V29cLzIlKFJaZ1wwTy5NV0p1W1c8J2BGXzFmdERpOyFePW47WF5pdTldMVFabXBjX0Q5OnAvQl9jQzJjLV5SLCRTWjtcQjJuUTghbl5VZmhdSzkhJi8/blpsRDI7SVouMjhDc19zJCotI24nb1tlLWQsWTdGP0QhKFFHPFRbTE9MY2dLJmtXZUM9IlZVVWtkS1VKRDEjJV4qRVlPbUxKPGRHRE5vMEIkPl1KZFY5L3RGWk1ZI2c6PlwpIyhZQlxIKkU9UCtvLyVnUD1SaSQnYyxIQUBeKT9WXEhqLVAyNmt0dWxLJCJGQlhHMnNiQ00sV1NeMERiJSZUJ0ExYzg2aGtsQUBLSXQlJiwwPmc3QlBaN2BZKEtsQDRPP1lBMV1NQW00W19JJyEqKFknOiNWQWMlMnMpNCQyKWpuKC5JVTs7dTtTLDUnbSIzRz1MR2hUUl00TGdPUkQhaTk5TF9qVj1cREQ7NytWIkQiM0thSk9FLGEjcHBwO0kkKjkmL25SVEouUm0sPkxyU20mVGoxQmJiazlXIy4nTSJbdWU3XTkkKThxWnM3OmtkLW1GQ1tsLlhANUhoaWIpWUdEVV5sdFhsL2oycytXIyVrOkc5bzIpci5GL05GYF0/ZGA0PTRXTERbTzMnUnA0Klg9JzFPRSNJQ2FnL1RQSSxfPU9IWl1UaGIrJFcqUlIqX0tobSUyb01bUCojaXInbVVyL2w6OXQ2Yiw3RV43PVpoNkBMIiNoWGksXTcnP0M+RkYoZWwnJ0xrQDNRIWVdT0twMipaIjo1Q11QVVpgMUVHcUpfNWtiOytLTFQsdEk9UF1eNnMjI0AyTUwiQTFAR2U4Z0lDI1lVVipgbURrSE4oXDQzRmcoYG9EJmxOL3REYjVgSUdPQz1qT1JiXUlcOmdNUk5cQEtlQDdhRk05Qk1gQ0I4akBsVTc8SDpEJmM+ITVGLihoKzNEQ0BAYTRmYEBbbFomT3FDSVUxZ0hWVDlZYjNIUEE8J2NbUSNjXilKWy5UXjMhWC05RTwkZmlJODJ1V1IibkBISEJrTjNnTFNBIXRHcnJxMmJUSVRDPyUsLyUtOS1YJm5BIU8xVWw1b3JhXG8hUkksPjFjOFohVEhUNVVzcG5pMmNPJUlGbllCMmY1ImdcKTItOChBb2FCUG45M1grZzVtRk5fQik9Z0YyO2VEPGtpbG9AL2VaO2NRMWM1XWBGRmNyclhmPDg2XWBnZk1oWG08c3MlLFRoPE1YImkmSXJAdG9ANjwqLVdqcU1iZjFfUUJtV3Nrc2tLYkJaIWJNdHRzb2I/N2k2RCEuIT82aVs5ZXUrVG5BOD5bUFkxQWxqKzswS0dCWztLNWUsdCJzLmo9ailYSCRCQUJHNDFuam9PJWg5Q1A7TGZAVjxYKXM9LmhGUVUjQGhwSVRWWy5bUCEqMDJTQC1WLipgJXNQTG1kNldfR28/Rjk4cl84cEU9YkFUQWBmV3BEJT9hcFJmXig0WGJwIyJxZVJPXChZImEvJU4pciIlPShgbV5uXllfbD0rPlZtKFh1LypHUXI3bm5YMmonY0JAQE9tTU91ZnEhWkNnKTIrOmhUb1VRYEQmW287MEclIU9OMGhvYV4pT1MwWDRITjYsblZeRFApJko+cCNNXztVZztTXTJobD5YcWA3bUZsXjNbTXAzIVtqZU87WUNtcDIyMEY+XWJHVmNbIkYkWEpgRm1XOj1sZCssV0xESzkkQi85RWM1J1xbdHAvRy1VaHJicFtTOyxCLW1jUD5KPz0tKUQ/SUAnOWczMidORDA1RTk6ZSU0V2NfPy9KRkYhXEpJSmBpJlw1WjBXZSxkOC5STDpgUF9tNmNLXktwR1xIK1ImNixBNEBDZ2ZFJnUiWFYuQm1sQkA5XyVhS187Z1EtO08/Qj1LMVFxKXIhbVY4SWBbZiZMUlprJlJdQCxvZWk7MkBrMSItJFdXcypFYkFkU0hZJDtlRkZib29rMnBVV2ZGREVzYkUtUS9aQSxOSCdQQ011JW4nbkROJzVNQidZRmdhVmAtYm9UIUBWQVljV1k3R2tJIUgoNCYnX2xKQ0ZkI2ZmZnRZaW0oLVhqczZANV9CQ20xMzBMQz9VXFhcVVdtXzNaJC8jaVxSQG43OXByXmxwXyZNSXE/Ky0rKFdYVEJcTiwlcz0sPi05dTVATS8ibiohRjBtcTNgU20uRTZuOSJYb3BUX11JNVdVUlhmaFoxXWk/OUpmUi0iYGxCazZgPE1mODc+J2U9RnFKLmZgJltcbmtSUlQ0VHRFU0o9UyQ5Yy9hcCItW2JOXycjQGFwbSJVMmNjTUg7N1BhYVVFczhpPGlwMHAjLzFGUWQyP2RWPiYtTnBFJTAjQm9UM1hFNGdwYCZgR2hYXVhpPz0uQGEiZC0pKWVeKy82KGZPYUNAcVdcSmhMTzNvQlRcdTg6OSFgclBkMGVPSXEzTDQ3PyV0Z2Y/OllNcipNKS1QY2ltNUZYcXRZbHBsMypIKXA+bVEmJWxebjc3WUI/K2hAc19CZ2tQNEREc1NHIiI+W283U1BKcS9wO1lrIycqbypvNVMqOl8iUkFcKU0sSDtraFhlMk1bRCNSWzhzSVJbLjlzRihhZ3BjcT9VZXRhakNIUk1TKjpzaW9TQFw7JFVZbzByS1A3RWxudCpbbDowXVFsVFYuOCpaOzlkaGoqL1csLmtBMzdbQXAoPDEnSU5eI2pBJ0Jtai9OL2tTWUMhL0s+XFdMRjRNL0BOXE8hOEBHLVQxJS9oXjJeYi46NEdbM2phdGVObztYbmYnIkkjcUlqc1ksTGVXJitnJ0peZl9hMkdrIVlzbE1idVxKIzxHI3FUblQpZFo6c1BqLlVvTllrXVh0J01KY2pVVyRdbGUnYnFiZGRJQSJkZEo8XzVCIl1mSSkkUCwmSzFFSWgkTG4uMShvSEhdZ0RDXXFWSGJfYWxTJmFHKzZjLlAlLWJpIkwmI2JfZGlbRF1XOlEyX1lMaFY2REtPPSU2UyVPQCdjMEtbbjFsWHE9SSs0Nk1kKEkjMCRBIWJwZUR1LWRuOEk0YSkhM10iUCxwbDxncy5wMmFZK100dEpiY2pZbydvK2p0JWROa0AiMUVVYkk1YCRkTGtpVXByXUByJ1VjTFtKJW9hdV1EKUBTbmVRaylNLGdVRmVPNmQrNjlXUFx0I3JdJy9YMEZqW0olbENwWjc9NGJCUGVGZ20pKi0+Ty1rZzJZMCw5ZyRzSlBkVmRfczlgXGE0Py0sJV9aNFs7a20zRjs7V2xbaDhdUlhaciNHQ0AkR1A7QUNzLitUR0pMZT46Vk5rb1I8S24jWGVkW2dMQCZZJHBdT187QXBYNFAwSnNicnElNUciSF1sWUkhcV9VQ0xiXmRDJl5qODhqK2tOaFsoWCJRNlg4bW1zLmxdW3JfK2UqS2xVXWFNVDFhOTxDbnUzQ2RVZVtgT2pjXFpYKzE8JGN0YGpmbVM0ak5GRzwuUEohT0A/ZFYmOisibGtGZ0NZLHFBOVtJLCN0VkdQRE9GMVNMakdWbVwlUT9nbTlJMk8lSCNEMEw4Wk1hcGpQIip0cSVSM2FBOVspQGsxX0FKLzc/IlI2S0Y2XEFmbkJKPCpqUTBFLWtOQGFfZFUlVTAnZS9nOmUnQSVQTl9dXVFZTzlXMSwhPmc/aUFqZm82QjgicmRBVUY2SzJLNnBePTEuPTY9RFFaP1JzKGc3VCUkM1E8IkgiMGExQGZqV05zZE9OKVckZy1cW1JPODlJVjliXTtTJiRIVCNTW0Ahc2dJVUVbYF5HNEMoWitpK0A3dGgvYztAaDtyI0JqO1EqLk1ST0g/RldFbEdqPDdDbUswJ2ZpXks5PyVrOGFFUio6TDYpS1hfY21lJFpEYzpAaWYtaicraTpEKCQ4LGlmbCplUiN+PmVuZHN0cmVhbQplbmRvYmoKMTAgMCBvYmoKPDwKL0ZpbHRlciBbIC9BU0NJSTg1RGVjb2RlIC9GbGF0ZURlY29kZSBdIC9MZW5ndGggMTMyNQo+PgpzdHJlYW0KR2F0bTpoZiU3LSY6WEBcRUtCPnFsYUpGVWA/KTxFY0FIMSVSTXBNa0AqbT4pQWc8RVY4TkJIQ3FBaCQ8OFwsQl0xO2VpWDs5LTw7XjpYOUAmLURdXW1eSGowSVRfNjokbUNSO1s9cD0pJm0zdSplUz1bLGA2NidvL0JHWk5qIyEyUEZBOz9oOSZULGA2PDRVMSFIQ01XLDIvNUoiJ0lcWkIjQG1gP1JDRWNia2ZTRkRTOVMiUk51JHRjRzwtVT1wUkRvRyFlNj5GLGtUJCVAZ15cOGctLGhFNU9SZktvYWk4LzdIXztLV1MyXzFTQHIncG0pYD9FJTBZW2w3OHJiWkVbXnFiaTRSKkI/cmJbPWJSRiMwYEZgYCsyXzxxQnBzQXQ7KzBfJyohSG8uckYmN0MtbVA0VDI1OyxucEpEWnFyPkI0OG82XmNAamVBb09TYUtZOmlBaUcqR05aYjAkJFVTMFRGN19JJXV0ZTNKXklpUDI1a2NZJHUxZyw2MV1VVGFLby5ZNUNkNSY/RnVdLGohQCM4NCgzbFc7Rzc2cDcpTCdzJFxLMERtQ2dbNGhjVjtiQmxhLiRMTWA7bVVkPSldIWljWEk+NkB1O0s+UEZFNzFdPWpNZWchJyVDRGZcZ0AxaCxoXlstU3AvbzA2MmBvXDFiKG4sRG8tW1VATzokZjRNTWw0I1Ynbk4taUdMKFpANXQ8RWsuOV1ERSJMWVJqVGg+KU0tWD8iVz4sJGNDVHV0aUlvMUQ5SkhgQUktdFVlLyQ3TkdcPSxUS1dKaTM2Pj1jbDJ1IUNjTydtViNHTDhuQUdoW1MvRS5qUCxodThkXSFXZU1hcjFXT2ttQk5POlRaO3JcXURfPGFSNFFJJ1RRL1VIJHE2aCp0XkAoNz4kPloyXCQxLUdnbGFlS29DKEloZ28sPF5oYlgicEldPWQiNGEkTmpkZCkqbyokXXBgKl5COEhqJEFnais5KmppYFVeUzc7VElKNDpzVyg0KTlQI00ldU1YJzYlQGByS19nNWAqYGcqMjlwYE5TTzVRKCM9MWZpJlEscVJDdG9oZEp0XWZoYEsoL11BUWg4VG0vWSU1TT8iL0BRLF1KN1ZRVi9QVWhFa2dQWXVVPXRYSj9yRV02W0kwSzxXQ1A7NEQxS1g5VitFbENeVDU9VS07Zi9CUTpTXmdTPmEhLi1gSjBvKTVPMTFyImg7V3I8VjEyUUdCVF9vaFpEKyN1OUI8UGNXLD5NdCQ9ZisodW1AYEhjbk5qPlkiVjAzMG8lNFdHXkRWRUZvPDQhaExSazgoMjldI1ZGM0wiP1pgQFpQJCdaZiNyQzdVNEUvNi1WPTYhXGVTPWk7IUlHPjhUREM3Ki9FTVQ8TVM4YUZPaVRKcCFOJmhVLjx1Rm0kKFZTYUhbPE1RZVZzWVtMJXVbPk1XJWEpZCkrOWZjQXRNdT5FVllbYFdTQTlNQUFJcyhiN0swLmtfVSpGNUNuX3VQI0wxKjE1cWA3J244ImZUKGNiMUZHbT1mI1onVVE2c0ZLcmA8XW5gLjVYYV1hOSglbGFabk9VLTpDR109JT5bSTk+ZzJCVzY7cDEvLWluaTplKG1iJFNnOillKl4/U2Q8RHI4R09jRmlZa0dgUV1oKVUqdVZxQFZMS2k2Qy88SGM3bGAyOmkwITJTSTAoJigyTlxXTCg6JXJicz83bWdzXz02LUUiPDgwXyQrLWRfWGpoaU00XnBcazkwa0MxMXRycnNUPnJXMlpFSW8kfj5lbmRzdHJlYW0KZW5kb2JqCnhyZWYKMCAxMQowMDAwMDAwMDAwIDY1NTM1IGYgCjAwMDAwMDAwNjEgMDAwMDAgbiAKMDAwMDAwMDEwMiAwMDAwMCBuIAowMDAwMDAwMjA5IDAwMDAwIG4gCjAwMDAwMDAzMjEgMDAwMDAgbiAKMDAwMDAwMDUxNCAwMDAwMCBuIAowMDAwMDAwNzA4IDAwMDAwIG4gCjAwMDAwMDA3NzYgMDAwMDAgbiAKMDAwMDAwMTEzOCAwMDAwMCBuIAowMDAwMDAxMjAzIDAwMDAwIG4gCjAwMDAwMDQ0MjUgMDAwMDAgbiAKdHJhaWxlcgo8PAovSUQgCls8ODZiZGQ0YjZiMGVjZmYzNTE5MzYzM2NjYjc4YmU0ODI+PDg2YmRkNGI2YjBlY2ZmMzUxOTM2MzNjY2I3OGJlNDgyPl0KJSBSZXBvcnRMYWIgZ2VuZXJhdGVkIFBERiBkb2N1bWVudCAtLSBkaWdlc3QgKG9wZW5zb3VyY2UpCgovSW5mbyA3IDAgUgovUm9vdCA2IDAgUgovU2l6ZSAxMQo+PgpzdGFydHhyZWYKNTg0MgolJUVPRgo=",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for native element-level MAGeCK differential quantifications from the genome-wide WTC11-EC CRISPRi screen.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 08:34:38,108:iu_debug:	DRY RUN is enabled.

All Document dry runs passed.
Marker: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/DOCUMENT_DRY_RUN_PASSED.json

```

## Markdown cell 19

## 8. Submit the missing Document records

After reviewing the dry-run logs, enter the confirmation token and rerun
this cell. The notebook verifies all four aliases afterward.

## Code cell 20

```python
DOCUMENT_PRODUCTION_CONFIRMATION = (
    "POST_4_FILE_FORMAT_DOCUMENTS_TO_IGVF_PROD"
)
REQUIRED_DOCUMENT_CONFIRMATION = (
    "POST_4_FILE_FORMAT_DOCUMENTS_TO_IGVF_PROD"
)

if (
    DOCUMENT_PRODUCTION_CONFIRMATION
    != REQUIRED_DOCUMENT_CONFIRMATION
):
    print("Document production POST is blocked.")
    print("Required confirmation token:")
    print(REQUIRED_DOCUMENT_CONFIRMATION)
else:
    if not DOCUMENT_DRY_RUN_MARKER.is_file():
        raise RuntimeError(
            "Document dry-run marker is missing"
        )

    marker = json.loads(
        DOCUMENT_DRY_RUN_MARKER.read_text(
            encoding="utf-8"
        )
    )
    if (
        marker["document_metadata_sha256"]
        != sha256sum(DOCUMENT_METADATA_PATH)
    ):
        raise RuntimeError(
            "Document metadata changed after dry run"
        )

    for inventory_row in inventory_rows:
        if (
            marker["pdf_sha256"][
                inventory_row["key"]
            ]
            != inventory_row["sha256"]
        ):
            raise RuntimeError(
                "A PDF changed after dry run: "
                f"{inventory_row['pdf_path']}"
            )

    production_rows: list[dict] = []

    for operation in (
        document_operations.to_dict("records")
    ):
        if operation["action"] == "SKIP_EXISTING":
            production_rows.append(
                {
                    "key": operation["key"],
                    "action": "SKIP_EXISTING",
                    "success": True,
                    "log_path": "",
                }
            )
            continue

        result = run_iu_register(
            profile="document",
            input_tsv=Path(
                operation["input_file"]
            ),
            dry_run=False,
            patch=False,
            label=(
                f"document_{operation['key']}"
            ),
        )
        production_rows.append(
            {
                "key": operation["key"],
                "action": operation["action"],
                "success": result["success"],
                "log_path": result["log_path"],
            }
        )

    production_frame = pd.DataFrame(
        production_rows
    )
    write_tsv(
        production_frame,
        WORK_DIR
        / "document_production_results.tsv",
    )

    verification_rows: list[dict] = []

    for key, spec in SPECIFICATIONS.items():
        record = get_with_retry(
            spec["document_alias"],
            attempts=10,
        )
        if record is None:
            raise RuntimeError(
                "Document was not found after POST: "
                f"{spec['document_alias']}"
            )
        if (
            str(record.get("document_type", ""))
            != "file format specification"
        ):
            raise RuntimeError(
                "Unexpected document_type for "
                f"{spec['document_alias']}"
            )

        verification_rows.append(
            {
                "key": key,
                "document_alias": (
                    spec["document_alias"]
                ),
                "document_id": portal_record_id(
                    record
                ),
                "document_uuid": str(
                    record.get("uuid", "")
                ),
                "portal_url": portal_record_url(
                    record
                ),
                "document_type": str(
                    record.get(
                        "document_type",
                        "",
                    )
                ),
                "attachment_download": str(
                    (
                        record.get(
                            "attachment",
                            {},
                        )
                        or {}
                    ).get("download", "")
                ),
                "verification": "PASS",
            }
        )

    document_verification = pd.DataFrame(
        verification_rows
    )
    DOCUMENT_RESULTS_PATH = (
        WORK_DIR
        / "document_portal_records.csv"
    )
    document_verification.to_csv(
        DOCUMENT_RESULTS_PATH,
        index=False,
        lineterminator="\n",
    )

    if len(document_verification) != 4:
        raise RuntimeError(
            "Expected four verified Documents"
        )

    DOCUMENT_COMPLETE_MARKER = (
        WORK_DIR
        / "DOCUMENTS_SUBMITTED_4_OF_4.txt"
    )
    DOCUMENT_COMPLETE_MARKER.write_text(
        "Four file-format specification "
        "Documents exist on the IGVF Portal.\n"
        f"Records: {DOCUMENT_RESULTS_PATH}\n",
        encoding="utf-8",
    )

    print(document_verification.to_string(
        index=False
    ))
    print("\nMarker:", DOCUMENT_COMPLETE_MARKER)
```

### Text output (stream)

```text

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_validation_guide.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/document_validation_guide_submit_20260930_083438.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:34:38,729:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:34:38,729:iu_debug:	submission=True: In submission mode.
2026-09-30 08:34:39,978:iu_debug:	
IN post().
2026-09-30 08:34:40,008:iu_debug:	Validating the payload against the schema
2026-09-30 08:34:40,026:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_validation_screen_guide_quantifications_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOSAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDggMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL0NvbnRlbnRzIDEwIDAgUiAvTWVkaWFCb3ggWyAwIDAgNzkyIDYxMiBdIC9QYXJlbnQgOCAwIFIgL1Jlc291cmNlcyA8PAovRm9udCAxIDAgUiAvUHJvY1NldCBbIC9QREYgL1RleHQgL0ltYWdlQiAvSW1hZ2VDIC9JbWFnZUkgXQo+PiAvUm90YXRlIDAgL1RyYW5zIDw8Cgo+PiAKICAvVHlwZSAvUGFnZQo+PgplbmRvYmoKNiAwIG9iago8PAovUGFnZU1vZGUgL1VzZU5vbmUgL1BhZ2VzIDggMCBSIC9UeXBlIC9DYXRhbG9nCj4+CmVuZG9iago3IDAgb2JqCjw8Ci9BdXRob3IgKEVuZ3JlaXR6IExhYikgL0NyZWF0aW9uRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9DcmVhdG9yIChcKHVuc3BlY2lmaWVkXCkpIC9LZXl3b3JkcyAoKSAvTW9kRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9Qcm9kdWNlciAoUmVwb3J0TGFiIFBERiBMaWJyYXJ5IC0gXChvcGVuc291cmNlXCkpIAogIC9TdWJqZWN0IChJR1ZGIGZpbGUgZm9ybWF0IHNwZWNpZmljYXRpb24pIC9UaXRsZSAoRm9jdXNlZCB2YWxpZGF0aW9uIHNjcmVlbiBHdWlkZSBxdWFudGlmaWNhdGlvbnMgZmlsZSBmb3JtYXQgc3BlY2lmaWNhdGlvbikgL1RyYXBwZWQgL0ZhbHNlCj4+CmVuZG9iago4IDAgb2JqCjw8Ci9Db3VudCAyIC9LaWRzIFsgNCAwIFIgNSAwIFIgXSAvVHlwZSAvUGFnZXMKPj4KZW5kb2JqCjkgMCBvYmoKPDwKL0ZpbHRlciBbIC9BU0NJSTg1RGVjb2RlIC9GbGF0ZURlY29kZSBdIC9MZW5ndGggMjk3OQo+PgpzdHJlYW0KR2IhI14+QkFRMSZVc0pYZk0zT1JTaFhyMkdLKkZdQVolXDZJIWNwMUlIY11HVkk1ZlBdSTNwRWdbck8xKEViMScsIWdzU1IzbixSK0dkX2lfZFhEaSJxTFdEbylHY29wdEhNVjpEYmMlVlVyVjJnUD4jQiopQ2QtbSs9JGM2NzJ0IXAmWzJDbk1jODZmOCM6MTQ4QEBeOidJNUkla2BtIVRuRmdibScoNz9uOShXUioiaCJPK0ReaXA7O1tHKSs5PVFDbi0rR1Y6JzRUR1MwWTFoLjpNdFxSUSE/SSs1byZpIVhvInRxLVBMZ1NKTjlLbF9bP1smLj4pa2Y2TltqSCFVQGBWVmJoWSRoNSxJTjRdRVY3azola2wncWhwSUpjMEE0S3AyWSIkXSpTT1dCaDRKZTFAXlRiR0NZI0FWOTsxUWBdPktfSjhJP2paKWFIREU9ZiFvKj1BITtpWk1mS0g0P1Y7Ty8+OW9GQkdqVFVjKFM2KT9Aby5scVhwcCFNbkcjMVsjIWE4Q2E4OW5RN04wQyJeYHAkcUBOVTRiT3Brb2dLIllWLnBscC9cOitKRDouSi5lPGlbY1xvUi1UJTpKbD49MmpuOV0sdClhUEFzZ09dQ0FILENHNEs9XFhTKWQ2IUczQmswQlBGJ2ciYkgoP0JJOT51NzBdL3UsJzNTWENnZVllIjtycCNtPENDIU4jTEM1P0g0TWduOFluZFBIQyo9UzI2Lk41ZFY4dTwyYF0xMmopW21Xam1MLTMlJVxRbjlaT29CbXNYKGcmO2wtLE1oTy1oNXVbY2dAamcoaC0lbWZGTlNNPEZMMS01UWZjZUBEWkQjOHU7QDk3aj5mMlRPYFo4KiEwNG9gbFxpPEMrTHU2JEFQOyU+OzhyLWlzMjlWZHVjNVgsXmRsc2lqIWNwTWZjXl9vKTI9J3VjKWNRbFVldWQnIWw9MVRZTyhcPDBgOFlhQC4kMjdbQE1iQTouVCo5R2FJTSIiRl5NMCFJNDZuInEkYzssLnJPYk1qIUthK2B0JEApPUYhSVAtNUFSUWJvKHJtNj8laiw1NUozST83aClCXEtvYW1zLjE4ZEFmVnQhREo+NEJWTUtabTshKygyVF03V09rQTRdPUYnZDFPMSRjP2BmYXU7UmIwPThMWjQxc2BhVDM9JUg7V1M5Wkg/QzIxIUNPNTxxNEpHP15yLDwzWT8mal0rKGJDMCJNPDU6P1A+WFZUZjVAXzAnWTxTOzNRMkkhLmJFcWxzLFQ2cmQ7KDAsR1M/QGZUaEtdNE5VVG4nTnIiXDpuUDxWRWwmJSU1MnJVMy9JSzBrSiMwTlpHIjBjUSgsX2Ymcy0yQWhXK2AhanI2LiJTVU1mMm8qIW80I2ZNMXNqaEBjcDFoWlEiSD1Gaz9vKSYrXGNzXic1P0s+W1hROUtnOzYkcE11UjZwMz9TPVchLmojXlNsbm5LY2VAL2NOSFhLNCNqJUJdOkBVXU5TSUlYbENRZlZZR1peUTcyYjhjJnNWKEE3czUpNlUrcVxJXjBPS2NVRTZMYHQ7I1pJS289T1A7TEAqOzorUWJnOGJzIjVNXk9OJ15FUT5WKyYqcGNXSnVIJkZocWQqXC84KWAySGdpYipIcTY1RTlvLlokR2AvOFI6N01gM1dZQjg/QHNsLlEnLmNiPWknJzhYKD8zczE0dDMyVF08bygvUVlra0hcLXFuP1ExWi9dOHNCTilTMEBOWkZWJ1QnciUiVSFccDsmaVZbT2RLbFxgX29BaW9eQT11biFrY19vaj5ALDZjRV8rVDU0c2JiWEtEWklCcUR0LHBEXW1JQCJkQjRnQlZZO1ZyXTxJODxDUHFCIj5SSFtMZEFMRT4uNEhTQT89TEU+K2NVbmFUKSY8N1svN2c1TlxZQmpDOnF0azEhXEVDaEdBXFJOSSJBQmFiRFw5JV4yQjQxKFhcOFkiJC5OLD1EJks4czw7SjwpKzJbWG5YQXM/Lk4rLideJ1sndGxtbVwvM0c8ZWlWZTRVOW0uLyolX1xSYlA8Qk9nZkxeN0JANllnPm9dS1w2MDZiISwzXClUbGAkZj4mQnVqaV5NKi1UNzQpUCQ6VCcmcyhvZ0ZSczkmaio+aC9eUEJIZ2tYTCgiNHBGczBUQF05ZkU1QDdzK2IuPVgvOUIlbmY8NltNXzg1VHVkLVk2WVBTYSNLLTpTPDxmXkUzKWlSXilFXG9ddSFAMm5YI1xCZ3NyZVpNW08sZHM1PTpNU2MsN0huYyY6OmNLU1tFKCs5VGxzaXJbYGZVT25GYCs5P3RFKkdXY15pZlZ1OTxucy9XZDpJcUBjXildUFk1a2VIWnFAaT8pXF5fLms2JSJUTm8waU8vN21jbmoyPW9WKVssa05JU0loO3JcJS8iMlltbD9DY2s8Zl9abmkhPjtJIW1EN0xkZ3RGJWtWI3FAKkwtYSZqIStfO3I4KGJrNEo9OzI1ZFBbQyloKFNROlkrKVxCZF48QDZlcG1pKEtQIkRsal4uTSZEVjBHdUhANCtPP11sP0g6KVZkOTQvWScsI0RLP09mPGxVT1BUcDlIaCM0YytwU3NhLTApMklANlYzY0JnL2koLHRbcD8jJUNeNT9fKl9tPClKP1RwJnVPOWIqc1hfRSUwXDpMWkVxU0JXS1lTRnNWMDZhSmRtZ2xoR1s2VnU0SWZwbVU3Y0s7TEFFJksobDxgYD5AazZmX10rXyQkZGxRQWohUlxVSGJpQDVTNDFAQDw4YjNQZG45OnNzWEgyPzk5P2pEY05bMUhBXl9tO1FVSF8xc15nRSxpaTJFUG4yZk9lTHFNPFYiLFk1dExcamJePiNMTmdOai1iXmxcYm8sX3JCaGhrPzlRU1coIkFGI1RAKE1oWWttKydqUzpMbTNrKi41Nz4iSXM5YzA5Lyc5bGdNYkNtcmkoQGFnaz5HMWNIRkJfKG4qcCMxJVdBdV9PJDJBXEZFYz8qQ0FBIz4lMFRhJSVMUFUvO3NObWpOdFFKYVQ0VEVoQVREKFJgbDhjS0tEWSckS0I6NiJGSU9fSGchZGNMPkZQVUc9XigwXSR1OjRHb1VSUFdTP2NTcWAnOlwwUTteQnBEaTdKVCMua0lLL2JSY3RGTWd1XlwrV1Y7ZkZmNCJZOER1JF9Fam81bywydShNZGZEWGJHb3NNTkhrLGlBdW9YbVYsRD48QmJAO1JpLVEhQEAkMi5wMDdQQkUuYTpVWiM7UjZOMSVhW1kvJE5RWihhXD8lImItJTFoaiR0Vjs/PlpQbzlBTW1UbGwkSHJOT05ucjdlYmwsQkwyXzQhW2QlUycvU1dnWkZjVzwzYUxAQDBvcjV1M1RCaDU7VGhbYjFINGpoajYwYkg8SCpyMUNfSz8jRVdyLjk1IVJKKzIsVlFaLW42OTRYY18tXzxaIWk+OSk7TUdcJ21uJlowZmZxL1BxKzphUVEjYDkjSjoiRU9rUWpbNjJfZVxqImAmKT88JGwpS21eWj5lZlZUQltSR1RqVk5eWUxIRVddazdGQTVdUG03WUFLYkBJZzlJaF5XbSliX0ZLOkQoMFc5RVpKQ2YqQTxhWHVLQStfJm9vQDotLCIxWT5ya24nQDU0TV1QRWYxXDBNY0hOYVlUYUtITW5nJTRoRSo3Ri89MmZkRjM6X0cvJ0MyOGFpO0JAYWFPJXEuY1s0Y1dQQlhdUDVCLU9MYDZmYGFMb0YySFZsPUgxTEREIlRTZ0suPy4tam5oZVdqbHMzRWIsWmFXOSguRnJYYVZYXiw8bnJCKlYwWzdNanAjWTNyZy1KQko7JTdcSnNKOTNoJ0ZdKWYxVHM2K1dEaFUkQT9iOEZUNmQ5JVhnWVZlcUwmR1xKRlpzNDttNVZ1ZEE9UkIrZS1FOWhJWyhrPiZkUSs9YCZFZEx1aEA4dW0jZUZQP0k4LCJtaz02PGo/b2tVPls9WHRMIVY/c2YxJn4+ZW5kc3RyZWFtCmVuZG9iagoxMCAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAxMDAxCj4+CnN0cmVhbQpHYXRtOWdKNktnJjpNbCtiVGo2U1tGXmJ0IVZNIVZKdVs3VjZJX28qZEZvWCckYmNhajliYj1VWmcva15jbWY8KE08cXEkUzJrWSZiZV4udTJwXHU2cGM4TG1iKDJXLjk3X0g2S0ByK2VjM0I9MyhAcSI6JVk5PjptdTBEVXFkSl01RS1qSTg7bGAxalkhVztIbVkiRCVpL2NRKTFYImBwWVZALDJdQV8kZG8vOSckMDxrQzAwXS9gMiVYIlglbi1xP2VcLSNbRmFecG1qLitMK1lNWzUoLG1IWEVbO09KYyo8RiNrNmJuPU5TYiNnJnBRKE1ILSYkc1U5ZjpLZ28lRU5WW0BlMyciI2dWO2M/UE4lVV11S0M4KllHMm5pIiNLTyJQJ247ITxiUlRyNz9hQUU0JnRlbDJaIztQcDc9Pk04cyFKUFw8bCZba0wrNiRMOzhfYW9JMSU7VjJzM0gmS1deMHUkKkllZ0tYZjglcGg3VGNVKENoXjQicW5xUkkkWiRjXTBDVVQ+TFFGSSdlXDxbXmJUR0YhQFw4PHBOLV11VEpLMT1VTyx1ZUhpNHU2aCwvPzhpLS8iK2JvVjw5UDhnLzxKNyZgbUxHKUEmIk8xclVBSTI2UU1ZY3JMb2EjXi5JXTU0MHBkbSpuPUBNTltyXz9UMGdFPi4mL11ibkVWdCRDUFJfKU9KIWpZbVxnMEY4UFFWZlxnN3QscVRwdCwlPGtVSFVtUlBBJXI4WnBvS1snMCQ9SnRrUmouSWQkSDZNZzhMKnQ3UDBYaWFWOUVlYzlSZVovLi5ybVB1Si1jWSsqXGRIQVYuUXMpaVJGWSolTk5ublZeV2pEVUE1Om1HMjkhRVkyakM5cD5WRV1Qb0NuQk1SNjtTW0UsZiNNMEY6OVRBTTIlJG1NP0NKZD8xY09APVRbSlsyJ144WkpNMSE8Ul4pJGw+KF5tP0lNTiZzXmAmb2I0UDx0UiY/Kyo1a0M5UWNCTkdAXUpdISlqLEYmPmxcRGoubSdCNjMpTSw8VExiLiJCMW4oLzxcQk8qKmg5RyslKnFoIiNnVV88RidZZ0xvWlQyVidQUEA3JEshTTAoPzYnR2E6bVJvQE1ROERBIiJpaTk2VFUjRykxaWtnbDM7anBmYz86QG5IOktDRTdTIyxPckEhc2NFaF0uI2NyZGlaMTZjTT9MRzoyRiIwUTJbPCpvKT1hYFpOMmAnMl9FXmxCLHU3PWUqNm5HIUJCTjgpWEQsXlZZb3RnQ007aTBbaD1TNj4zQy1OZmNwWjZXVHFlUTwkOVw0Kz1EN25qLUk4PzBgNGopI09QRlVRQDh+PmVuZHN0cmVhbQplbmRvYmoKeHJlZgowIDExCjAwMDAwMDAwMDAgNjU1MzUgZiAKMDAwMDAwMDA2MSAwMDAwMCBuIAowMDAwMDAwMTAyIDAwMDAwIG4gCjAwMDAwMDAyMDkgMDAwMDAgbiAKMDAwMDAwMDMyMSAwMDAwMCBuIAowMDAwMDAwNTE0IDAwMDAwIG4gCjAwMDAwMDA3MDggMDAwMDAgbiAKMDAwMDAwMDc3NiAwMDAwMCBuIAowMDAwMDAxMTMwIDAwMDAwIG4gCjAwMDAwMDExOTUgMDAwMDAgbiAKMDAwMDAwNDI2NSAwMDAwMCBuIAp0cmFpbGVyCjw8Ci9JRCAKWzxjNTQyYTA0ZjhiMWY2YTU5NjRhOWJiNmMxZjI2MGJhZD48YzU0MmEwNGY4YjFmNmE1OTY0YTliYjZjMWYyNjBiYWQ+XQolIFJlcG9ydExhYiBnZW5lcmF0ZWQgUERGIGRvY3VtZW50IC0tIGRpZ2VzdCAob3BlbnNvdXJjZSkKCi9JbmZvIDcgMCBSCi9Sb290IDYgMCBSCi9TaXplIDExCj4+CnN0YXJ0eHJlZgo1MzU4CiUlRU9GCg==",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for guide-level read-count quantifications from the focused WTC11-EC validation CRISPRi screen.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 08:34:40,478:iu_debug:	Success.
2026-09-30 08:34:40,479:iu_debug:	Object posted with identifier: dc8d807a-4038-47ac-9006-16b6e5b64a80


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_validation_differential.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/document_validation_differential_submit_20260930_083440.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:34:41,038:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:34:41,038:iu_debug:	submission=True: In submission mode.
2026-09-30 08:34:42,169:iu_debug:	
IN post().
2026-09-30 08:34:42,196:iu_debug:	Validating the payload against the schema
2026-09-30 08:34:42,214:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_validation_screen_differential_element_quantifications_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_validation_screen_differential_element_quantifications_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_validation_screen_differential_element_quantifications_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOSAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDggMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL0NvbnRlbnRzIDEwIDAgUiAvTWVkaWFCb3ggWyAwIDAgNzkyIDYxMiBdIC9QYXJlbnQgOCAwIFIgL1Jlc291cmNlcyA8PAovRm9udCAxIDAgUiAvUHJvY1NldCBbIC9QREYgL1RleHQgL0ltYWdlQiAvSW1hZ2VDIC9JbWFnZUkgXQo+PiAvUm90YXRlIDAgL1RyYW5zIDw8Cgo+PiAKICAvVHlwZSAvUGFnZQo+PgplbmRvYmoKNiAwIG9iago8PAovUGFnZU1vZGUgL1VzZU5vbmUgL1BhZ2VzIDggMCBSIC9UeXBlIC9DYXRhbG9nCj4+CmVuZG9iago3IDAgb2JqCjw8Ci9BdXRob3IgKEVuZ3JlaXR6IExhYikgL0NyZWF0aW9uRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9DcmVhdG9yIChcKHVuc3BlY2lmaWVkXCkpIC9LZXl3b3JkcyAoKSAvTW9kRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9Qcm9kdWNlciAoUmVwb3J0TGFiIFBERiBMaWJyYXJ5IC0gXChvcGVuc291cmNlXCkpIAogIC9TdWJqZWN0IChJR1ZGIGZpbGUgZm9ybWF0IHNwZWNpZmljYXRpb24pIC9UaXRsZSAoRm9jdXNlZCB2YWxpZGF0aW9uIHNjcmVlbiBEaWZmZXJlbnRpYWwgZWxlbWVudCBxdWFudGlmaWNhdGlvbnMgZmlsZSBmb3JtYXQgc3BlY2lmaWNhdGlvbikgL1RyYXBwZWQgL0ZhbHNlCj4+CmVuZG9iago4IDAgb2JqCjw8Ci9Db3VudCAyIC9LaWRzIFsgNCAwIFIgNSAwIFIgXSAvVHlwZSAvUGFnZXMKPj4KZW5kb2JqCjkgMCBvYmoKPDwKL0ZpbHRlciBbIC9BU0NJSTg1RGVjb2RlIC9GbGF0ZURlY29kZSBdIC9MZW5ndGggMzU5Mgo+PgpzdHJlYW0KR2IhO2doL2g9LyZxQlhfWXJERT9HIlJNSUxJZkVcNTAkRkk5OlhML0g3dGk6Ol9aO1JqYkgwYUcrQGU9SWdeZSYuK29qa2tqaV07MWEwQmUwIVArWlorRipNWWpgLzhKJmFSYl5aUCVWUkswVSReJzxTYzAkTjRoZWJuZSJSMk0qKkFkKVFITmI9LCJycCdgPW1EQTdbMVtPZFI8ZStdP0s4bmFrT0chIlVrLGBFMm5FSSNjbHVUbSxjZG0yQlZRVDVoPzNtRUBpPG9dZV9FVDExaHVtT01RbEtGOydoTGpDcUUpP0pRPS5NXzJoLClNM21CKFNiOjhRZj9sTCNebFNPWGVabi1aWHElb0Z1YltoNSNxTXBXInJAUT9WOypvZ0IzcmIzcm1gImVmSEVOOzlpZClfKm5SYWFbMTotbFlvdFsnNUNpQEIxI2ZcTyREJzEhb1pnSTYycWxXTUhqV1NvYEpmISRNZ2I/PGtxXj4qPUBeczU1TVVMSyFdTUAnbz4xJyc+VUM6ZypGKnJZRTluWTEyV05nP05gWFonKXBudSRNUScwJipbM0dmM21ZTG80P1gtJT9IU1R1aG9dLSxcS191RzArRCkxPFZCOlBBMlsqRmdVSnBzXy5mOUFCWnNuRTNgM2VIQHBmUHA8PDhsJjxiRlRTUEJOUW4kJ2U6L2c3PylAb2U7a1JGJDxHM0dAPFpEMmxiTyRDSGN0VlQqNjNrI2h0P1xvTE1IJFNkMUpGLEU/Mlw3LUBDXW8lPkNgXyFNXkJpMClgTUZbcEE2QT9qPlEhU2syb1AkZjljYGJuKz0rblMrKCZJVDMiUUltLFsrRzFBcD4zVz80I0haJUpjLURdU2ppbm9LLGFSOlo3JDk6LE9eRl8vWl1BTUgjMEtrcilOTHUlXkBSIys7Jm46USxSOi9tclQ4P0NmYHNRcmVlVzZJdCtMXTIuVl9AI1hhUFBiZ0s5SlFrW0pQKStvcTtaI0svYDZaJjk0dS5BMTFpQjpyWVg5I3VJWkxVYj9Tc29iVlo3MDg+UXFbI2VZXTYvSl1yXXBYOD87TGk+SDFWRnJLJTRbQWlAXUkhKSRpMzEtVUJZU1txNVpQWzJIYSEyXEVMbl1rKkxgPD0tbmlDJEA7KGovbW5ZREQrXzdIQkIkMWhhWydbMmNyXS8iT2pkQ2BDVUUwPjJCQU1JSTlDbjg+SUBBLFp1cTEyVXJvT20/SkxfTjI5aVc/ZUplOiRxI0YtT2VkQzZzJGkhT3BmTWxFZ2YhT18nJDhiOFc5WnNNMkUzOidYamhqKkNmSix1cEtoNTZCRVNhWistQDBqXWZQTnFbU3JNbChDUFgoNTFsblkscGt1LkdVcVZlMEZqQltYWz82NUM7JjBBXDJcUysmIih1NiJIR1ppIVRmb1JLKjA2U0hBJ0NsSm5JX1AzZGRKWlxXQHVmLC8/JGBxWHNcay9KR0xtMGQpcGVARnIxVyhdYFd1OUpPNGInUSIqJTIiJj5CNVEjMVhLcHNHIkxwVzE5dCVlQ0w0UDRaX0FTbU5FO2VGWnNuJVZUIikoIVkvaGBPMWpKaS9VOFc3KGguYG4jZGcsbkpzMl5CcCNXQDhBOnVALDsvUl5gVTFuVytwUkZrRixrYzs8VFJKSDA5STc5cjEpaFQxTlBsdVJUb1c6YXIkL1s1Y2BJbCd0V1JPOVc5akNCaWZcT2k1JXBpJ3JpWy9LU01XR1MuTUNQWE46akNobWQwWGVqcik6SlwiODVUbFM8cF0wMnFxSDttVGs3Z3JrRmJfMTI9bkkhbkdrMlBKc0ZTIydxMmtUIi0iJ0FpSzRUZyNHcGU8P3JzOD5gY1JIJjVccy1HT2ola0sockxRVDQ5MkVIcU5jVE9mMkxTNWdpMF5wMVVlcTkyPSplYj1xJilEU0prRkskMFlOdF5AJipbPGwxc3RUJmAoJWA+MXN0VCVpPSk5Mk5bLE1gLSpabDhsbVNVVTsuZmAjN2hMYEBpV2EnQjQsUl1uJ19GKkhTZyMqT2tkV0cpJXVyUm83NjhcQF5AR3UiaHF0JmEsOS8obTQ3YW0pSzFfdSRJbUFzJFJyQmctMlxzYignNEE/RCtqbCEsajdEXTUkYikwTEg+ZCYxbydFTjtiblxrNys7Zj1UNHBFJlskSCMmVCdIYyFPTl5DXlFHUjAvKFNRKjZSSnA6L14tQHE2KTNMVm5zKlc0altGRWAiV2A3My0/WXNbV09hV1RBI1M5JkBeb2YrakNocFRvLTJxdD9CQ3NYZy87MyMjJkAmcTFCRmgjYzdvQUNhV2s6VVNzZGdJbHA8LVUqV21QO2c7STkzJk9nc00+S2YuJChMTy4tSGJiPUQuXGRxPSwmXkddZHA9NXQoO1oxKTc/cnA9RVwrXU1LKzlTaigkJmstJEtDViE3cSVjUE1YbUs5MWJKYGgkOyhBNTJsLWgtMiNcNCplRk10KT10Nz9dcSRGQmhOXFFMcFM8NT1MRFZJMCNDTj5qIV0rb1UtTWwkOjNFJ2VGVT5qRCRrQk9bcGJcKmk0KFwuKzQvSyZkRk9MOCxAKy5ZQUNKXFNSaCwsNj9PO046S21WP0A2LC5hQUxyJGtSV2QrRmJIUnQsOjcpNipIZyouXyQrLE5zWT5vZWIlTG8ubDVaR2BzU0Y7Tz9vXVI7cyVSMSFEdF07KGUsWjtuXWtCP2VTck1XbD1HJmZVZVgydXIoZ0xjUHQscT5xLldcTlknbGJrNjcyVzZ1cUMwUlRqcC5cVDN0YkZmZmwoQ0BDIyg9XkZoNGpEcGBgPmEuIUR1IyxGXF1BdVo5dGI3VSJVbHRhMmltMSRAJlZEPzpMaDNebzQ9PT0mWkBBRTVRcD9bZk0tb0dkQVNeM2EjPEgxPXA3Wkw7L15EZjVWPEw2TjBTQiRFTFpmZGVvcTFhcW1qTFZwbFdxTzohMSRTcG1vM2goNiVEIUNhITJiMk0qO0kzWGdUZiNXTjxHLU87YUtKZ1pCTjZYRk0+UFZLXDI5LWsjKGNxcGtXTCchRTU7cjxqPz4iYWJrUVc1PjFWYilOUytpZVlvP1QwRTFsSGFTJGouNlVTcVEuOzVMJW9PSkRYJWU1MEJHSFBsM2ZRYVUmJ0EkN109JiIybiwsdCVRXy5FZ19kPVBWYmtOT05caiZSJk9Fb24iS3E6TkV1KSlEaj5GMidlJVdfZVdOYERLXSMtIWJRWS0oSGw4MUw8OlhXYnRPM2teQz9jJXFFZUpFPGE+V208LUFuTnMnaVdUSXMwcEdjZSI8W2FcKkRtRkhDc3IkQTgkPmhJJC06OFAsN1gmOTAlXU4rIUIxKlY1OzsuNHAlVzNxJkklI3NtQ085REttLnNyKlN0TShwKWBrJiovXG9mI3RMXnE/ZFIyNkdRKFpHVl1tKVFzKyM0TmxSJDBFKCIyYlhHPVdPNj9jNk1TVzlxNDJacFBHLzBCI0k/bE8sYGxuYlN0NGlqdXUlJ01aOzAya3NVZy5OYVgkPnVqJ2orUClWRE89ajY3IkVSa2tsaHU4I2o7Mm82aWlNKj9mJEcvS0AhQzE8MlI1XTcpO1JjXFAsXD1tSlNFO3E6J1FfUy5UQmU6b1ZWNDVqIUtCSlNeJEZWQXJCWTZXbVlqQE1eWWAjW2JcL247OXNANyZQJ0o0J1VOQyFkQVZRdEk7alUjbDtaRkNiYSMwNSc2cWlBNUA0bCMjJ2MtSGUpWFYiPWQuRDVCamU8ZTNQSjYjW25maUc5RyI5ZSsyKitPXTRBNGskWFs4UnFmSzR0NkM6NWY1XWo9YWApPDc1VC9CZFMnbTAsOispJVNGImZOI1tKbkB1NG9wdD9zYE9LV21rdE5qKURKcS4wKFtJPmxiLiNYVWIwNGMzaGpVaylLWFdoNTJpbEM9bltobSNETWRRckplRkErTlooQVcxO21DdCJxUzRESDp0XChFUjVFaT51VyhrMiU9MU9lO11CQ2w6dCpIRls7NSspPltna29rdTQ3PUQpMiVpUzZzcGxhJ1lbcT82MmErS2B1JDc7VGVZYVVpSXM2UENARzpcWFNWIWxVODtULmRLSjxgbz1yZGtuTDw6SClULFdiWD47ZWY3bTU9VXE+clFiSlE2ZCYjJnBOTk5yZVVJKDsybFolJTpAX2VrLDVNcC9razhZc2heND1oNywwcmlSM0NMPWRRRGNGTW9OPzcuWF1OWWVZLmZGPTVjMms2JUMtWGhVQjZpSWxaJiNnVVByJltRKGYvLlZcZ2BNUlhuc0ZMLCMwMDZBLWNVXCRUOlpQaTM+VTtWckRqa2VyXm9KLSxpZ0lKIyNLJ2xET141aV8xKDBOYWMnO1ZhVS9FJGNvQ1xSSi0xNnFXKyk2WGZXKTo7ZTVzKiw+JyMicU91Lz5tLnUiMDRlL1RIcV07R2ZdJm5bcSFTWjo3Km9CcmF1NEdVYChhTz9mOig3WEQ2Xl0oQkM5TkdBU1xjWklKQF1wWWwvOzlyU2knLiQvUEBMW1xgJmpmWFtvOTxfakFkYzk0bi9IQz82bW1WL0VHUWAjYzxqaDxvYVshSi9oO1hTN25lXkhwJ0VfZEBORDZhKXVWWj9GYExuRiYzODpWKlZyNTJxUElsSyI8cFA0SiFyaWY0KydGZjcjZEhDdTJtOXVeJjRYZWVyamhTUFpLJk5aYV80bGdBVXU3ITNqbmQpc2lWRi51LWFvJ3Q8JCxhNFA+cT87VXRZT0RyO2hGcmc3YTRKTjxEczxELyZSXSFIQGE+MCdFJCt+PmVuZHN0cmVhbQplbmRvYmoKMTAgMCBvYmoKPDwKL0ZpbHRlciBbIC9BU0NJSTg1RGVjb2RlIC9GbGF0ZURlY29kZSBdIC9MZW5ndGggMjY4Mgo+PgpzdHJlYW0KR2F1MEY/WlhzWyUiSFUkXmdqZ1JBUCYxPUgrPyk2OF5tWSZcaDxVQVdIVUArVTE/LnBrdShFImwwWkEhLFVIVDJBaSNfLjk6TFNbI185XSQhIiRbSihWYCoyNFhnRDwwNEIxZUpganArIiMiVHA9TUJESGJIb04tOkAnUnBEQEI5JiRpTk1FNUpvSlMiMVdjK0wvYTRIKitoRDZpbz9VSThyMT88MU4sYmVyLWQtT05ZTEJcQkJvW1BkaWFKJUR1aV90aWUuT2w8dVxlUTNSLTIpO284PkxTa0k1azBDWSInP2luKldyVD1BPzxgJGtqLzRnPCRMKmVIXSQpcz9ONUpaaCcvOlVYKSJxNzgwS1RgKV5jT1Y7IWkiKFlEX0BxLkNJZEZ1YUtUWTtpTTNLO0toW0FWR20qJkhWNWluYWBlMmhZPm4rZktxYm5USU0kR1pnRW03NjY3XjA3NyJscFVrQGoiMD9fZ2xLNlJgOXFJI1guUVhdJGJLXGFKMHFbQ05hKnQwRzBdJWRrcS9ELWhqSm1TQi9IMnVqL2hGOWUlIi0yInB1OTVYcDFgWzZhRmA7Km5iVHQxb14xRWNbNFNpci9dUixYL0BtP2JnQV1wTDdWanNyaGlHaGghdD9pcW0yXCgiaW1qXiZrLGFpXiYxJW49OTZVQyVHZj5vX0ojRFxfdT5zTkRPQDZhPzRtMVxgUV5uRHBrKVFdaiNOM3EmPERTTC9PSUs+N2FPJDsoaFk8TyE5I0BQSmxmPypON1pVL3BFUi0rI2llWWYqKCJabzBKRVBrJ0dxOF9IPzBJc2BqPWVhczBwJl4pbTcwZ2E5azRKKW03MmVKImxoIkRAZ3JERDFpJ1tqJnVEJXEzR0BjNSo5cTNBNSx1TjcjbHUhNVZUO0hjWktnJEY6QiJYIzc7MXRCVk9GLiFANV5jUTBYbTdrUm0ibmU4I0RHKzM5MF1GXjhKbzw2ZGNgN2RBX2YwTGNnM05iPyQ7ZVw+Ki5qRktHUysiPzI8RSNkM1wlYGUrUFhgLDkrXHFMW0tULi5dUUhHbVllbTNRRzReclElNV5dUSMocVhLP1lmWVReUzBHN1RwIl5DImkucjkmTD4oN0wmYFhAJihdI00oY1VNRHUxLyVXYSxoTitmS1EiPG1gRk9pTTkhaTQrbDAlOXJKKlgjWFhpZlY3XkQuKkohR2gqNDg9aD1WW2ZSZHFQSi88Y1E0RjkwMmBxLVdqQ2NxJV4rTm07UWQ3TUo3TkJMLDZLX1ZSXiQiWW9ySGE4KCM6OzFMMGZtKSdUQVUmbFQpcEg7cF8jI1wtJiVsUlpUVCdhcSNNT0ZaPGhQRmJdY1VxMSssP1FON0VCIVpQV1teTEVccVA4Vi1KaGBhc0FGZCNOZlxTMD51PWEhSHBTYmtBQTJlR2gyPlFFXl86MSUkVkV1YiY/MEU8X3RYbW9lXU02ZmFfYE5HNXNIXmBAWl8jNTpKTCshX0g9Sy9AL2xuT25lTmJXR0ZdbD9MW3NsU3UpImllXUk4WyxgZXNlTFluV0paJCxjL28vSFVJbC0sIjxeN09XVHEkWyo8V09BX21xcCc8QGY1NWwrQV9ZKlksTHVHUj5GZj4kZTRtUSQ6K3IlYlphWmNbY0grNT5yUnRZVzQuKDMrUU1aVVoqVU5OQGJbO2o+bzd1WiFCU1JdO1NBaSxnMFZdZUJuI05UaV01QGc3YjZ1WVxcIXQ0OkFUWDlKO2Q5akUoI0BYMzxxTiVoZEpfZzs3MVpzRWg2bFdOb1VvRjlTXzdtTCYlV00zS0MsNSYjIi1rTV4iaWxuRCU1VEVEKVBXVTdXUSEzKUN1ZDdcXDUsIjlBMkZFRGw4J0tZSEpWLXFFdUo0STwuMmheN25EOWwudGdIQGMnSENCIkkuZmY1V2NraS5fdGc7VT5ZUFcjJTxhOWxTNE8zMDcndCVcREV0L1k+XSw0UzBwbDpwRCtDP1NnbmQiSXN1NiwxRiVeI1AyR14rajxwWztZPzwrPEVONVFrQ3UlZVJMKVhzLzBQMGYoJGMwZ0dcYjxIOmY6VSVpRVZgZGdHVWlhZUg7cU1CSykxOiNtczlZTWBmbCJYNWBWIzhZJ24xUiEiXVlwNjtaYEc/JUBKWGtVOD4lKXErXDZsVk0/I3ApQ1JQXThlcCphKjpyQiNtWic7UHVTKGU4JDQxRD5KNy9POm1hJScnJ0twLj9xMWRrWCc8YEJrdCIoLiY2cHRpMVdHNyxzWj5uWGpxRnVVazczOlkodSJRZXFoMUFxPkBOSFxwKml1XEBLPiVWKEsiP25SUlJvQnIwKGtTISRDMEJTKGg1PU9IdS07QUlldWBtWUhzcWFxQlgnYlxwYF0lbm4zViFoV3VwOWdKSU47KjNzaV0/S09SVWM6dTUpdC5SXGg2S08oXCVGTWpqYyIqakhPbElmK3RWYFxbSk1AQExlSlhSbU1QcUxHIUpTLj8+byNgIShQWG9UO0I1YileKT5PbmJlSUBKJ2NUUVdJNXUvQ3RfMHBeY2gmI1JGRj1aLnFyVCtsRzNUP1lOMUZealclRT1BNjtlKTNKZic7MVZXL29zU25rSmJMITQzKkBfVURgbFMoVScxUSFROCNIcWJOLjYrKUVZdCsscGNlMklFVmUtOSlIRG0lImBuWC1aKDstR1hDL0xLKmVHUjBvIS5JclNkKDJhQVIhTkFfLXAvY0wvcmhaUUdIaSFzYyFWOk8lPEAzaGQtZEVzK1VgL0xnUlRgTXE6TjY+UU5KcDE7ZlNCYC4mczwsRGIhIT9vUjBhclhbPVRRXnQiVTtTajBuRVBDMSovOlZcMWZKQ1NROGFacXMiXnMkUmAsQVhtJShFbkJpcz1jPDgoND5VQ28kNTVwYlteKjhQLyJHKiRccnFMbCIhbis6LzJpZCRvUldhXCE/N2xIU3AjZVBXQ290U2BbQDw8LClMNWxSbEsoNDJObiphJ2xUYyFTUV9WPHRtV2FyK2BUdHJvdGhycGxVTkBOTyhCc2FFUD1oUjtjUkYsPT1bJVshXHJdY1tQKTloaE5LcSE6Ii8rJG5KSF06Py8lO01WPV1qPUBaakQtb1gjayFZa0c8Nl9EciQ5Q0o5bFcqUz8nJzghVSotS2cxZz1hPVFHKj9tPl9BO3VWWmExMT9ZRmFFMyEwcldLVyRsNWYvViwzWCpMUWdpdEs3ZVZgKlVVWTItPC9rIWxYU0FOWmVcSlNALjNoKiUiT21dXWczTW5QPjBMWl0wJmRMOSVmZWNuZEVwI1pLRER0VWAvWVBLSV8nVk5tW1ArVz5hbWJFJzNfKC1BdSQuYTc0V0lDOjJGMlc3LkVZOzVFKCIjRUxuTDtYJ1o/M1otYCtrYWhMTTgmdUcwOTBTVlRLQG5WYjY9VUV0Imdca2ksNkdYTS1pUjpTVnQ5WC1gJDUkUC5gaUk4S0FZPVokSmUnOlFOayo6STNLRydsOmtyMW4pInFKRSZAOXRiOExFTiFyc3AnP3BzQT42KTs6KiEqcWwiNCsuX0c5TnNPNmlINGMvLlxkJEd0SkdMSDFIJ34+ZW5kc3RyZWFtCmVuZG9iagp4cmVmCjAgMTEKMDAwMDAwMDAwMCA2NTUzNSBmIAowMDAwMDAwMDYxIDAwMDAwIG4gCjAwMDAwMDAxMDIgMDAwMDAgbiAKMDAwMDAwMDIwOSAwMDAwMCBuIAowMDAwMDAwMzIxIDAwMDAwIG4gCjAwMDAwMDA1MTQgMDAwMDAgbiAKMDAwMDAwMDcwOCAwMDAwMCBuIAowMDAwMDAwNzc2IDAwMDAwIG4gCjAwMDAwMDExNDUgMDAwMDAgbiAKMDAwMDAwMTIxMCAwMDAwMCBuIAowMDAwMDA0ODkzIDAwMDAwIG4gCnRyYWlsZXIKPDwKL0lEIApbPDcyOGUxMjMxNmRmNThmNzA5MWY0ZGIwY2FhYzJmZTYwPjw3MjhlMTIzMTZkZjU4ZjcwOTFmNGRiMGNhYWMyZmU2MD5dCiUgUmVwb3J0TGFiIGdlbmVyYXRlZCBQREYgZG9jdW1lbnQgLS0gZGlnZXN0IChvcGVuc291cmNlKQoKL0luZm8gNyAwIFIKL1Jvb3QgNiAwIFIKL1NpemUgMTEKPj4Kc3RhcnR4cmVmCjc2NjcKJSVFT0YK",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for promoter-element MAGeCK differential quantifications from the focused WTC11-EC validation CRISPRi screen.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 08:34:42,852:iu_debug:	Success.
2026-09-30 08:34:42,852:iu_debug:	Object posted with identifier: babab92f-024f-4734-8a59-55134d3d86b7


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_genome_wide_guide.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/document_genome_wide_guide_submit_20260930_083442.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:34:43,402:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:34:43,402:iu_debug:	submission=True: In submission mode.
2026-09-30 08:34:44,863:iu_debug:	
IN post().
2026-09-30 08:34:44,890:iu_debug:	Validating the payload against the schema
2026-09-30 08:34:44,908:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_genome_wide_screen_guide_quantifications_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_genome_wide_screen_guide_quantifications_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_genome_wide_screen_guide_quantifications_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOSAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDggMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL0NvbnRlbnRzIDEwIDAgUiAvTWVkaWFCb3ggWyAwIDAgNzkyIDYxMiBdIC9QYXJlbnQgOCAwIFIgL1Jlc291cmNlcyA8PAovRm9udCAxIDAgUiAvUHJvY1NldCBbIC9QREYgL1RleHQgL0ltYWdlQiAvSW1hZ2VDIC9JbWFnZUkgXQo+PiAvUm90YXRlIDAgL1RyYW5zIDw8Cgo+PiAKICAvVHlwZSAvUGFnZQo+PgplbmRvYmoKNiAwIG9iago8PAovUGFnZU1vZGUgL1VzZU5vbmUgL1BhZ2VzIDggMCBSIC9UeXBlIC9DYXRhbG9nCj4+CmVuZG9iago3IDAgb2JqCjw8Ci9BdXRob3IgKEVuZ3JlaXR6IExhYikgL0NyZWF0aW9uRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9DcmVhdG9yIChcKHVuc3BlY2lmaWVkXCkpIC9LZXl3b3JkcyAoKSAvTW9kRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9Qcm9kdWNlciAoUmVwb3J0TGFiIFBERiBMaWJyYXJ5IC0gXChvcGVuc291cmNlXCkpIAogIC9TdWJqZWN0IChJR1ZGIGZpbGUgZm9ybWF0IHNwZWNpZmljYXRpb24pIC9UaXRsZSAoR2Vub21lLXdpZGUgc2NyZWVuIEd1aWRlIHF1YW50aWZpY2F0aW9ucyBmaWxlIGZvcm1hdCBzcGVjaWZpY2F0aW9uKSAvVHJhcHBlZCAvRmFsc2UKPj4KZW5kb2JqCjggMCBvYmoKPDwKL0NvdW50IDIgL0tpZHMgWyA0IDAgUiA1IDAgUiBdIC9UeXBlIC9QYWdlcwo+PgplbmRvYmoKOSAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAyODA1Cj4+CnN0cmVhbQpHYXVITmgsRSchJnFCKUprYzknaFtaSmpwXjtFbCU1Vj8raEAyITtKamYwYF4qZlJFJEUoMlYmZUdccGg0WDE4QVlpSDVaLG1EbDBCRFFDbWUnbS0pSlYzdWpxRTwiJnJdLFBpRXNAVEBqOUN1YUhTRSFcWmNVJUUiJDl0M2knXnNDWTBVU2c1USoiaSpkbWQuSkAhZjMiaDVObzMzSilcTFg7K0xUR18+MW9xNV81JyZUV3RDTzZMaGNPMWFYIzI9QGskTzFXOyRXczBvQF1ZajppVUVwTjg8TUk2RTAzaTFKa21KPDVEUF5iJzYkbW00SzUjWE1IUDdoY0I6ND9NZ3JQJzQldClIPVpxTTUpa3VFX09sb2BNJkFXVVVJbV0sVW5ZRUFXTCQwWk0iNjdUOTpmUzhqQD4nKiZAMipPV20hclBtSSpSO15BMkFHT0E7NEUhOjVJdVEqIkcyalU8Sj4hOj9NXFlSMXVCPyM5UiMocVFXWU1caGcncSdzakw1PSkncCMnX2IhR01TJGVEcVVTW0huLllRPm9xKlx1ZzpnVVE1dG0/TG1vaWQuZnJjIls3VG9BbGVbZnE1XlE3O0phOWojPmZhZG9yWSwjL0g7VHNcUjRLdEAuRV1qOCNtRmVRbEMvTkw0b2sjWjBfLiknYi1iP3MsXC8sQlM0Ni9GKjskOXRhcTpoOSddPVNNbVEtZlVNSUU2NyIzXC1OdVAlIm9kT0QwbFMmKTBwLC5cMChRIktwcUJDZVQ6RG8wKEVbJlYxJ0pNWnBVXVUrOl5SWWlCPmJAWSw4QUw+O09fRE9IWVNJWV9zMilYKT0oZFIibk5OaW8uWWlTRD1lZWIlNmNObSFaSTMiTVVtS1w7bnAhaD82TS5UOjVxUDM8UV5EZyIwbUVMNjtSUWdmLjxaJGUpK05CaipuOCRRdVk/XEklalI6Yy5cXU8qOUknSk11QmNCQzpkLWx0MWM+aFZOOG1WPWw+ZkQ5U3IvISsrVkpfW0UvNCskS09BP3Q3SkpYSS4zZDFMQWc8K0ZMSiRpMktVOG0lZ0VjMHNxLl0ja1kyOXFuaWtSIWA7YzNyPG0wWT1mKEkmVT40Z2JbIzNsYztscXBTXDxmZS9vL1dMbzEmI2tIVTllajgpIThqMlcwc1FAXEFqaD9jazctTVI0R3BMOk1pXEtrVEtDUWNvLFJFPEtnXypcYlJyKGFZK2Q/LStwIWo5R1N0QSErbz03XllLMzo8NnAwPHIvNzViLytlZUVaZEY+aSEpKHVQQlYyMk9nMm5WRkxJPkY2IkVRTE01OCJRK2FrMTg/Zms2KE8sbDgsNTpUTGhHI2FWPjp1QSxiMXVcYVZqXEMuO3JxLTNHRENVPGdbY2JvOU5BTUlHRC9zTV86bEc1SSI+RjpmRi0zO0VNUzs5KiZpIl9gUC46ZV8kIWB1ZlFLIStfUiZAb3QzZT43O10hUSdcXjZFYGprRSZaaTcuTHVvIkQwU1FTYlFcPy1UVSFCWyg8T0pdT2xEdUlTSGpNMjZ1WTUlLG0iNmVLJ01mR0ZHNTQxUC5Wbi1tRy5VTW08ZW5EWyM4SCNuInNZR2xjb3JZMVUqTiREKT4yLS9vLSI+TGtjTlNLX2FwMjdfOVRQZkFzdChbLGpxSFtiZWlwYktCQnJxa2ZkP2lMPTEvXWw9aGFXbXNmQmxaOGw/TV04NXMyTF0iI0JncUUvdUZPYkFcOSFdWXBjRFE5R0soRiFqP3M9PVBpK0pIOjRUJVUmPDhkbWxeRWtHcXVSb25NMVdHVjArVVBsPWMtT0BscCt0W0JxcC85PVE6STV0ZCVTXVQ/PFMpUlBgWl4xZiQyY3IqIU1TMylUSFprVWE3R1hrVmA6LTInQDQ2Oi5LLjZARFFIaEdqX110XzdGT2g1TSw6WmNeKjxhTVFTbmZISUZcKURyRGhRRW1JS2lCKVNbNVpjWFdaS1dsWyRlTVdBYTNRNTk2N3JVcFpNPmw0MnJQO0xSPyU+NWBFPjdYNF0sPz8sYXFXMDYkWEtMQUdXMEJIXEpIPWFnPHRjSkZmaDxLJU4qMVc+dFxvS1B0V1ddUEZHRm1bUCYiYTpFMl87UCY7clkpQTZHQjY4Mkdla28/SWg9bkcuQC8vaCFVTWtXZk5eW18+QEA7c1EwVFNSMHBUKi1FX2FTSk1DQmkvIVtwRlFDQygpaE9yIWBZUz08PWlLVmZYXW1hUD5gYissZlxwVm00KyRNTmJaKWhYQThzLEdbdUVmRERlaEYwRisrRVEsa0ZGKFhhW0UsWHEuT2hwW1YkLWRyNjpvOFM4Sk9QNkc+NjZyNHAhMFNjPWthSWs1SGp1cV1YK1BzIz9hc0gmai9sbzd1WT9wcCQ9IWBIbV1nS3JDXTBKLUNucyFHJztPQT0wbC5qWGEtKV9ITmwxL0peP0BbOypWNHJDZ2FmYixIXTdZbyIrIjw8WzBcSGA9VktuSlw4RDlSVystWF1KXHU6JnRhcmZlJXAzXGw6NUU8SUopb1spTjcxWS8+OXFtLDsyS2woVCJhWkhHcSIzSGtWSmsrIzU8MGNaXFNgbk5nMiMjXi1aYVE8UURKU2BWWSdaZkY8MD9fM15gJ09QdCM9QyR0bCZvPHFXKV02PGkmJyswUDo3OjRXYmVmI0kkb2giSSdNJipOcVBzM0FeOitPKixhMzxNWVdoRDFEX19UMFFaKkchJko5JDJvJlxuRTtENkY+KmAwZ2syJlc7SlhIJWdgJWEnTXIzN0FWQldOJlwuZiFqPm1pMCVAWkpiQ1c3QlpZJFdrLGAuVzduJCZkcl1ycklgIUsicEtNaDhIQCImSz48UGJBVyRYYyY7MWkxXmBwYFk1Kl1sOVtdXlBcUF1bcktJLSUmP19LPD9sZnAlcExOZjlhPkBrI3NlRGU1TUN0Tm4tPmRNNG5TJkpMOVsmKmctJmVDJWdiZzR0M0FHOm0nNys3c19FTGAwbUNMXy5PR3QibGwvQDRrYDAoJUhoRlBYQSVEWDZEbiJAV1hmNVFQbW5dSFhbPFpRWi5AbUUsZlY3J11NTlItb1pRYj9wMktkSyg1QF9AK209dENhb0pjbVplazReNC4qPiI+aCJWcVdwZVtcZW5ZJl9Gayg5WDwhQC1mNjw0ITdoUSUnTVtEJWErRWZUQVJUOU5oSEFfa0xWaVJAJVkrVy4/PWFQV29AJlJeMzBVYnArT1k1WTQ+Xys2bixPMlsyYW9PQ1cwTyovNztSUSEwUFskOkkxWmYlcy9APWdcRyY/Z2UyPmMhQD4jOWg4azVyOWsvMT5FOWU3KzVNcHFddURLLkhWa3RgYnUoSmNkXzQyYERCUVFvXXBTR1gsRUYwLi1YYCIyUzpcJF0nZzEkWnE7KklLSzktNDs2UF0jZENxKTNFaStfQWRDYEROOiFwYl5uRyhxRWMoV14yZ0dDRW5Jb148ZTlMYjBSWkNeKC8paF0hUWBdTC5oNE9vTEJyQmBjNDlBME9hQTtbWjYkU20yKzxqNFUpI2giMjNHNkJJVz4iPlVlRkxRTG1wQGgyVmNkS2YjWEhbYjVtLnRgdU1xRjpVI2RjPk82c0QsWExpdS5ZODkvXToxLCtTX0N0IkEsIl1gTnNoTGwyXDJpOjZlSU9lWXVhOlpbaCtsTypvWEEwLTM3NXFOcjI4SVhlWCpITW9HVWxubUttUWF1bGghUCNmaUIoJnJXZEpRQEIwfj5lbmRzdHJlYW0KZW5kb2JqCjEwIDAgb2JqCjw8Ci9GaWx0ZXIgWyAvQVNDSUk4NURlY29kZSAvRmxhdGVEZWNvZGUgXSAvTGVuZ3RoIDgxMQo+PgpzdHJlYW0KR2F0biM/I1EyZCdTYylKLydfXWEuNXA9WDo6b19RVk9IKC0/Jk5rYC41RztIL0tZanQ8cD9nRjBwIjUrMmw1ZjY2bTNEPnBPMT4mQGczLWQ9VDk+LCErLydjOSpXVW4wRWwmJSlvPSw4I1VkQTdgK09YOnBjczVNMGByQXVOcGlHcWdnIVRuZiZOUC1iSlpVbGYoYCduazYvXkRtZz9ScU41RDE7b2xtXlsqXnNyIzFSbkgkQSlJPz4tSVctN0teNDxzUWtIODQ8UnFHVFQmPkpyMFM/UXE1RjdKUzY5bl41LSEvSHI3NF9AVEYhXi87IiwoazBJYWhcTkp0RCh0WmhMUEhicWZDX2hXLGQybl1rSjQzZl0xI2tvPyk0KiNESSg2Ty8vKlQnYTA5MlNkYlBGRDFmX3FDI1kqbl5KXWFEXXUwLlpvZyRpRC8kUSEpLCtDMTFQWDRQO01YUSdeTTNgJzVTMD1eQFwyMWU6TFVjZ0cjOVwnVypIY01daVk+U0pRRHJyVTVsaWA8IWpcQ2FycVlrPGJlclZuRzF0WE9zQjxFMyI+WnUtQUM0VichRCJbVVshOzkjSG1XcG1dc0pJZS5ZamJuKjVFPilmP1MuM0pNQ1Q7YWhibVhSWEstaTZjQWFaSFtFWm8tXkwwX1IsMiUkI15JUCxhLVxhJWhPU2JHOmlOWWolViFsT24+KTYuQj9UVydjVUR0QTZARWYrTkw9XEMwYDU7ODQ2X1A3TzVVYmhZLm08TlsvKFxoRFVNcmdBKkBQKWotIkk8M14rdD1Mb2lhIzBVa3RMQ1dKKUpGPjVXY2AnMDpkJWJOIVdPR0Q7a0clJE1PPVJfVmNSQCVXRGcpZU1sXVVPKVpjMmBmOjFhLVwtOFs0O0xhTSw0W1M4am8lRm0mXCRVTzQ3SzJeaUxdWEM+X102U1ZUcGJeNExOVisnKU1pbHJSZ2BfXz5nZnU0Yz43SElqXCZRXXVpazVUc3F0WyhVSXRQPTJsZ0ZSQ2xxYEskKD4/PzhAWSheUi06ZkVzPy4wRHNtLmQvOWInUV44PjVBOjMtQFoiX0JfOkV+PmVuZHN0cmVhbQplbmRvYmoKeHJlZgowIDExCjAwMDAwMDAwMDAgNjU1MzUgZiAKMDAwMDAwMDA2MSAwMDAwMCBuIAowMDAwMDAwMTAyIDAwMDAwIG4gCjAwMDAwMDAyMDkgMDAwMDAgbiAKMDAwMDAwMDMyMSAwMDAwMCBuIAowMDAwMDAwNTE0IDAwMDAwIG4gCjAwMDAwMDA3MDggMDAwMDAgbiAKMDAwMDAwMDc3NiAwMDAwMCBuIAowMDAwMDAxMTIzIDAwMDAwIG4gCjAwMDAwMDExODggMDAwMDAgbiAKMDAwMDAwNDA4NCAwMDAwMCBuIAp0cmFpbGVyCjw8Ci9JRCAKWzwwZWM0MzZjNjllYmNhZDFkMmM3ZGRhMjI1MDk5N2MxNj48MGVjNDM2YzY5ZWJjYWQxZDJjN2RkYTIyNTA5OTdjMTY+XQolIFJlcG9ydExhYiBnZW5lcmF0ZWQgUERGIGRvY3VtZW50IC0tIGRpZ2VzdCAob3BlbnNvdXJjZSkKCi9JbmZvIDcgMCBSCi9Sb290IDYgMCBSCi9TaXplIDExCj4+CnN0YXJ0eHJlZgo0OTg2CiUlRU9GCg==",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for guide-level read-count quantifications from the genome-wide WTC11-EC CRISPRi screen.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 08:34:45,246:iu_debug:	Success.
2026-09-30 08:34:45,246:iu_debug:	Object posted with identifier: 1672694b-f1ff-4f48-9da1-37899a3fa49b


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p document --tries 3 --delay 10 --backoff 2 -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/POST_document_genome_wide_differential.json
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/document_genome_wide_differential_submit_20260930_083445.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:34:45,785:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:34:45,785:iu_debug:	submission=True: In submission mode.
2026-09-30 08:34:47,015:iu_debug:	
IN post().
2026-09-30 08:34:47,042:iu_debug:	Validating the payload against the schema
2026-09-30 08:34:47,060:iu_debug:	<<<<<< POST document record jesse-engreitz:document_WTC11-EC_genome_wide_screen_differential_element_quantifications_file_format_v1 To IGVF database with URL https://api.data.igvf.org/document and this payload:

{
  "aliases": [
    "jesse-engreitz:document_WTC11-EC_genome_wide_screen_differential_element_quantifications_file_format_v1"
  ],
  "attachment": {
    "download": "WTC11_EC_genome_wide_screen_differential_element_quantifications_file_format_specification.pdf",
    "href": "data:application/pdf;base64,JVBERi0xLjQKJZOMi54gUmVwb3J0TGFiIEdlbmVyYXRlZCBQREYgZG9jdW1lbnQgKG9wZW5zb3VyY2UpCjEgMCBvYmoKPDwKL0YxIDIgMCBSIC9GMiAzIDAgUgo+PgplbmRvYmoKMiAwIG9iago8PAovQmFzZUZvbnQgL0hlbHZldGljYSAvRW5jb2RpbmcgL1dpbkFuc2lFbmNvZGluZyAvTmFtZSAvRjEgL1N1YnR5cGUgL1R5cGUxIC9UeXBlIC9Gb250Cj4+CmVuZG9iagozIDAgb2JqCjw8Ci9CYXNlRm9udCAvSGVsdmV0aWNhLUJvbGQgL0VuY29kaW5nIC9XaW5BbnNpRW5jb2RpbmcgL05hbWUgL0YyIC9TdWJ0eXBlIC9UeXBlMSAvVHlwZSAvRm9udAo+PgplbmRvYmoKNCAwIG9iago8PAovQ29udGVudHMgOSAwIFIgL01lZGlhQm94IFsgMCAwIDc5MiA2MTIgXSAvUGFyZW50IDggMCBSIC9SZXNvdXJjZXMgPDwKL0ZvbnQgMSAwIFIgL1Byb2NTZXQgWyAvUERGIC9UZXh0IC9JbWFnZUIgL0ltYWdlQyAvSW1hZ2VJIF0KPj4gL1JvdGF0ZSAwIC9UcmFucyA8PAoKPj4gCiAgL1R5cGUgL1BhZ2UKPj4KZW5kb2JqCjUgMCBvYmoKPDwKL0NvbnRlbnRzIDEwIDAgUiAvTWVkaWFCb3ggWyAwIDAgNzkyIDYxMiBdIC9QYXJlbnQgOCAwIFIgL1Jlc291cmNlcyA8PAovRm9udCAxIDAgUiAvUHJvY1NldCBbIC9QREYgL1RleHQgL0ltYWdlQiAvSW1hZ2VDIC9JbWFnZUkgXQo+PiAvUm90YXRlIDAgL1RyYW5zIDw8Cgo+PiAKICAvVHlwZSAvUGFnZQo+PgplbmRvYmoKNiAwIG9iago8PAovUGFnZU1vZGUgL1VzZU5vbmUgL1BhZ2VzIDggMCBSIC9UeXBlIC9DYXRhbG9nCj4+CmVuZG9iago3IDAgb2JqCjw8Ci9BdXRob3IgKEVuZ3JlaXR6IExhYikgL0NyZWF0aW9uRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9DcmVhdG9yIChcKHVuc3BlY2lmaWVkXCkpIC9LZXl3b3JkcyAoKSAvTW9kRGF0ZSAoRDoyMDI2MDkwNTAyMDQwNiswMCcwMCcpIC9Qcm9kdWNlciAoUmVwb3J0TGFiIFBERiBMaWJyYXJ5IC0gXChvcGVuc291cmNlXCkpIAogIC9TdWJqZWN0IChJR1ZGIGZpbGUgZm9ybWF0IHNwZWNpZmljYXRpb24pIC9UaXRsZSAoR2Vub21lLXdpZGUgc2NyZWVuIERpZmZlcmVudGlhbCBlbGVtZW50IHF1YW50aWZpY2F0aW9ucyBmaWxlIGZvcm1hdCBzcGVjaWZpY2F0aW9uKSAvVHJhcHBlZCAvRmFsc2UKPj4KZW5kb2JqCjggMCBvYmoKPDwKL0NvdW50IDIgL0tpZHMgWyA0IDAgUiA1IDAgUiBdIC9UeXBlIC9QYWdlcwo+PgplbmRvYmoKOSAwIG9iago8PAovRmlsdGVyIFsgL0FTQ0lJODVEZWNvZGUgL0ZsYXRlRGVjb2RlIF0gL0xlbmd0aCAzMTMxCj4+CnN0cmVhbQpHYXUwRz8jdU9HJnEwTVhrWmRyRGgtKU1sUyt0K1dOY1REX2lsc05qSFhCWi5LYS8mNmJKNDdKRVU7QyRJcDNPbkNzUyZaWlVsVyQ2SzVAb0Y8a3VyUnQ2OUNSSl9RPiFyMkt1OitkRDwpZU1gTlkzRWdpNGhSPStJNkI+UShBVXU6Qic9Sl9qUGlyZyR1SGZgaF02LitdOyNFSSlfWiNULUEzKjo/PWBuTC4+TFxoL1Q7aE9fXUlCTj4+V29cLzIlKFJaZ1wwTy5NV0p1W1c8J2BGXzFmdERpOyFePW47WF5pdTldMVFabXBjX0Q5OnAvQl9jQzJjLV5SLCRTWjtcQjJuUTghbl5VZmhdSzkhJi8/blpsRDI7SVouMjhDc19zJCotI24nb1tlLWQsWTdGP0QhKFFHPFRbTE9MY2dLJmtXZUM9IlZVVWtkS1VKRDEjJV4qRVlPbUxKPGRHRE5vMEIkPl1KZFY5L3RGWk1ZI2c6PlwpIyhZQlxIKkU9UCtvLyVnUD1SaSQnYyxIQUBeKT9WXEhqLVAyNmt0dWxLJCJGQlhHMnNiQ00sV1NeMERiJSZUJ0ExYzg2aGtsQUBLSXQlJiwwPmc3QlBaN2BZKEtsQDRPP1lBMV1NQW00W19JJyEqKFknOiNWQWMlMnMpNCQyKWpuKC5JVTs7dTtTLDUnbSIzRz1MR2hUUl00TGdPUkQhaTk5TF9qVj1cREQ7NytWIkQiM0thSk9FLGEjcHBwO0kkKjkmL25SVEouUm0sPkxyU20mVGoxQmJiazlXIy4nTSJbdWU3XTkkKThxWnM3OmtkLW1GQ1tsLlhANUhoaWIpWUdEVV5sdFhsL2oycytXIyVrOkc5bzIpci5GL05GYF0/ZGA0PTRXTERbTzMnUnA0Klg9JzFPRSNJQ2FnL1RQSSxfPU9IWl1UaGIrJFcqUlIqX0tobSUyb01bUCojaXInbVVyL2w6OXQ2Yiw3RV43PVpoNkBMIiNoWGksXTcnP0M+RkYoZWwnJ0xrQDNRIWVdT0twMipaIjo1Q11QVVpgMUVHcUpfNWtiOytLTFQsdEk9UF1eNnMjI0AyTUwiQTFAR2U4Z0lDI1lVVipgbURrSE4oXDQzRmcoYG9EJmxOL3REYjVgSUdPQz1qT1JiXUlcOmdNUk5cQEtlQDdhRk05Qk1gQ0I4akBsVTc8SDpEJmM+ITVGLihoKzNEQ0BAYTRmYEBbbFomT3FDSVUxZ0hWVDlZYjNIUEE8J2NbUSNjXilKWy5UXjMhWC05RTwkZmlJODJ1V1IibkBISEJrTjNnTFNBIXRHcnJxMmJUSVRDPyUsLyUtOS1YJm5BIU8xVWw1b3JhXG8hUkksPjFjOFohVEhUNVVzcG5pMmNPJUlGbllCMmY1ImdcKTItOChBb2FCUG45M1grZzVtRk5fQik9Z0YyO2VEPGtpbG9AL2VaO2NRMWM1XWBGRmNyclhmPDg2XWBnZk1oWG08c3MlLFRoPE1YImkmSXJAdG9ANjwqLVdqcU1iZjFfUUJtV3Nrc2tLYkJaIWJNdHRzb2I/N2k2RCEuIT82aVs5ZXUrVG5BOD5bUFkxQWxqKzswS0dCWztLNWUsdCJzLmo9ailYSCRCQUJHNDFuam9PJWg5Q1A7TGZAVjxYKXM9LmhGUVUjQGhwSVRWWy5bUCEqMDJTQC1WLipgJXNQTG1kNldfR28/Rjk4cl84cEU9YkFUQWBmV3BEJT9hcFJmXig0WGJwIyJxZVJPXChZImEvJU4pciIlPShgbV5uXllfbD0rPlZtKFh1LypHUXI3bm5YMmonY0JAQE9tTU91ZnEhWkNnKTIrOmhUb1VRYEQmW287MEclIU9OMGhvYV4pT1MwWDRITjYsblZeRFApJko+cCNNXztVZztTXTJobD5YcWA3bUZsXjNbTXAzIVtqZU87WUNtcDIyMEY+XWJHVmNbIkYkWEpgRm1XOj1sZCssV0xESzkkQi85RWM1J1xbdHAvRy1VaHJicFtTOyxCLW1jUD5KPz0tKUQ/SUAnOWczMidORDA1RTk6ZSU0V2NfPy9KRkYhXEpJSmBpJlw1WjBXZSxkOC5STDpgUF9tNmNLXktwR1xIK1ImNixBNEBDZ2ZFJnUiWFYuQm1sQkA5XyVhS187Z1EtO08/Qj1LMVFxKXIhbVY4SWBbZiZMUlprJlJdQCxvZWk7MkBrMSItJFdXcypFYkFkU0hZJDtlRkZib29rMnBVV2ZGREVzYkUtUS9aQSxOSCdQQ011JW4nbkROJzVNQidZRmdhVmAtYm9UIUBWQVljV1k3R2tJIUgoNCYnX2xKQ0ZkI2ZmZnRZaW0oLVhqczZANV9CQ20xMzBMQz9VXFhcVVdtXzNaJC8jaVxSQG43OXByXmxwXyZNSXE/Ky0rKFdYVEJcTiwlcz0sPi05dTVATS8ibiohRjBtcTNgU20uRTZuOSJYb3BUX11JNVdVUlhmaFoxXWk/OUpmUi0iYGxCazZgPE1mODc+J2U9RnFKLmZgJltcbmtSUlQ0VHRFU0o9UyQ5Yy9hcCItW2JOXycjQGFwbSJVMmNjTUg7N1BhYVVFczhpPGlwMHAjLzFGUWQyP2RWPiYtTnBFJTAjQm9UM1hFNGdwYCZgR2hYXVhpPz0uQGEiZC0pKWVeKy82KGZPYUNAcVdcSmhMTzNvQlRcdTg6OSFgclBkMGVPSXEzTDQ3PyV0Z2Y/OllNcipNKS1QY2ltNUZYcXRZbHBsMypIKXA+bVEmJWxebjc3WUI/K2hAc19CZ2tQNEREc1NHIiI+W283U1BKcS9wO1lrIycqbypvNVMqOl8iUkFcKU0sSDtraFhlMk1bRCNSWzhzSVJbLjlzRihhZ3BjcT9VZXRhakNIUk1TKjpzaW9TQFw7JFVZbzByS1A3RWxudCpbbDowXVFsVFYuOCpaOzlkaGoqL1csLmtBMzdbQXAoPDEnSU5eI2pBJ0Jtai9OL2tTWUMhL0s+XFdMRjRNL0BOXE8hOEBHLVQxJS9oXjJeYi46NEdbM2phdGVObztYbmYnIkkjcUlqc1ksTGVXJitnJ0peZl9hMkdrIVlzbE1idVxKIzxHI3FUblQpZFo6c1BqLlVvTllrXVh0J01KY2pVVyRdbGUnYnFiZGRJQSJkZEo8XzVCIl1mSSkkUCwmSzFFSWgkTG4uMShvSEhdZ0RDXXFWSGJfYWxTJmFHKzZjLlAlLWJpIkwmI2JfZGlbRF1XOlEyX1lMaFY2REtPPSU2UyVPQCdjMEtbbjFsWHE9SSs0Nk1kKEkjMCRBIWJwZUR1LWRuOEk0YSkhM10iUCxwbDxncy5wMmFZK100dEpiY2pZbydvK2p0JWROa0AiMUVVYkk1YCRkTGtpVXByXUByJ1VjTFtKJW9hdV1EKUBTbmVRaylNLGdVRmVPNmQrNjlXUFx0I3JdJy9YMEZqW0olbENwWjc9NGJCUGVGZ20pKi0+Ty1rZzJZMCw5ZyRzSlBkVmRfczlgXGE0Py0sJV9aNFs7a20zRjs7V2xbaDhdUlhaciNHQ0AkR1A7QUNzLitUR0pMZT46Vk5rb1I8S24jWGVkW2dMQCZZJHBdT187QXBYNFAwSnNicnElNUciSF1sWUkhcV9VQ0xiXmRDJl5qODhqK2tOaFsoWCJRNlg4bW1zLmxdW3JfK2UqS2xVXWFNVDFhOTxDbnUzQ2RVZVtgT2pjXFpYKzE8JGN0YGpmbVM0ak5GRzwuUEohT0A/ZFYmOisibGtGZ0NZLHFBOVtJLCN0VkdQRE9GMVNMakdWbVwlUT9nbTlJMk8lSCNEMEw4Wk1hcGpQIip0cSVSM2FBOVspQGsxX0FKLzc/IlI2S0Y2XEFmbkJKPCpqUTBFLWtOQGFfZFUlVTAnZS9nOmUnQSVQTl9dXVFZTzlXMSwhPmc/aUFqZm82QjgicmRBVUY2SzJLNnBePTEuPTY9RFFaP1JzKGc3VCUkM1E8IkgiMGExQGZqV05zZE9OKVckZy1cW1JPODlJVjliXTtTJiRIVCNTW0Ahc2dJVUVbYF5HNEMoWitpK0A3dGgvYztAaDtyI0JqO1EqLk1ST0g/RldFbEdqPDdDbUswJ2ZpXks5PyVrOGFFUio6TDYpS1hfY21lJFpEYzpAaWYtaicraTpEKCQ4LGlmbCplUiN+PmVuZHN0cmVhbQplbmRvYmoKMTAgMCBvYmoKPDwKL0ZpbHRlciBbIC9BU0NJSTg1RGVjb2RlIC9GbGF0ZURlY29kZSBdIC9MZW5ndGggMTMyNQo+PgpzdHJlYW0KR2F0bTpoZiU3LSY6WEBcRUtCPnFsYUpGVWA/KTxFY0FIMSVSTXBNa0AqbT4pQWc8RVY4TkJIQ3FBaCQ8OFwsQl0xO2VpWDs5LTw7XjpYOUAmLURdXW1eSGowSVRfNjokbUNSO1s9cD0pJm0zdSplUz1bLGA2NidvL0JHWk5qIyEyUEZBOz9oOSZULGA2PDRVMSFIQ01XLDIvNUoiJ0lcWkIjQG1gP1JDRWNia2ZTRkRTOVMiUk51JHRjRzwtVT1wUkRvRyFlNj5GLGtUJCVAZ15cOGctLGhFNU9SZktvYWk4LzdIXztLV1MyXzFTQHIncG0pYD9FJTBZW2w3OHJiWkVbXnFiaTRSKkI/cmJbPWJSRiMwYEZgYCsyXzxxQnBzQXQ7KzBfJyohSG8uckYmN0MtbVA0VDI1OyxucEpEWnFyPkI0OG82XmNAamVBb09TYUtZOmlBaUcqR05aYjAkJFVTMFRGN19JJXV0ZTNKXklpUDI1a2NZJHUxZyw2MV1VVGFLby5ZNUNkNSY/RnVdLGohQCM4NCgzbFc7Rzc2cDcpTCdzJFxLMERtQ2dbNGhjVjtiQmxhLiRMTWA7bVVkPSldIWljWEk+NkB1O0s+UEZFNzFdPWpNZWchJyVDRGZcZ0AxaCxoXlstU3AvbzA2MmBvXDFiKG4sRG8tW1VATzokZjRNTWw0I1Ynbk4taUdMKFpANXQ8RWsuOV1ERSJMWVJqVGg+KU0tWD8iVz4sJGNDVHV0aUlvMUQ5SkhgQUktdFVlLyQ3TkdcPSxUS1dKaTM2Pj1jbDJ1IUNjTydtViNHTDhuQUdoW1MvRS5qUCxodThkXSFXZU1hcjFXT2ttQk5POlRaO3JcXURfPGFSNFFJJ1RRL1VIJHE2aCp0XkAoNz4kPloyXCQxLUdnbGFlS29DKEloZ28sPF5oYlgicEldPWQiNGEkTmpkZCkqbyokXXBgKl5COEhqJEFnais5KmppYFVeUzc7VElKNDpzVyg0KTlQI00ldU1YJzYlQGByS19nNWAqYGcqMjlwYE5TTzVRKCM9MWZpJlEscVJDdG9oZEp0XWZoYEsoL11BUWg4VG0vWSU1TT8iL0BRLF1KN1ZRVi9QVWhFa2dQWXVVPXRYSj9yRV02W0kwSzxXQ1A7NEQxS1g5VitFbENeVDU9VS07Zi9CUTpTXmdTPmEhLi1gSjBvKTVPMTFyImg7V3I8VjEyUUdCVF9vaFpEKyN1OUI8UGNXLD5NdCQ9ZisodW1AYEhjbk5qPlkiVjAzMG8lNFdHXkRWRUZvPDQhaExSazgoMjldI1ZGM0wiP1pgQFpQJCdaZiNyQzdVNEUvNi1WPTYhXGVTPWk7IUlHPjhUREM3Ki9FTVQ8TVM4YUZPaVRKcCFOJmhVLjx1Rm0kKFZTYUhbPE1RZVZzWVtMJXVbPk1XJWEpZCkrOWZjQXRNdT5FVllbYFdTQTlNQUFJcyhiN0swLmtfVSpGNUNuX3VQI0wxKjE1cWA3J244ImZUKGNiMUZHbT1mI1onVVE2c0ZLcmA8XW5gLjVYYV1hOSglbGFabk9VLTpDR109JT5bSTk+ZzJCVzY7cDEvLWluaTplKG1iJFNnOillKl4/U2Q8RHI4R09jRmlZa0dgUV1oKVUqdVZxQFZMS2k2Qy88SGM3bGAyOmkwITJTSTAoJigyTlxXTCg6JXJicz83bWdzXz02LUUiPDgwXyQrLWRfWGpoaU00XnBcazkwa0MxMXRycnNUPnJXMlpFSW8kfj5lbmRzdHJlYW0KZW5kb2JqCnhyZWYKMCAxMQowMDAwMDAwMDAwIDY1NTM1IGYgCjAwMDAwMDAwNjEgMDAwMDAgbiAKMDAwMDAwMDEwMiAwMDAwMCBuIAowMDAwMDAwMjA5IDAwMDAwIG4gCjAwMDAwMDAzMjEgMDAwMDAgbiAKMDAwMDAwMDUxNCAwMDAwMCBuIAowMDAwMDAwNzA4IDAwMDAwIG4gCjAwMDAwMDA3NzYgMDAwMDAgbiAKMDAwMDAwMTEzOCAwMDAwMCBuIAowMDAwMDAxMjAzIDAwMDAwIG4gCjAwMDAwMDQ0MjUgMDAwMDAgbiAKdHJhaWxlcgo8PAovSUQgCls8ODZiZGQ0YjZiMGVjZmYzNTE5MzYzM2NjYjc4YmU0ODI+PDg2YmRkNGI2YjBlY2ZmMzUxOTM2MzNjY2I3OGJlNDgyPl0KJSBSZXBvcnRMYWIgZ2VuZXJhdGVkIFBERiBkb2N1bWVudCAtLSBkaWdlc3QgKG9wZW5zb3VyY2UpCgovSW5mbyA3IDAgUgovUm9vdCA2IDAgUgovU2l6ZSAxMQo+PgpzdGFydHhyZWYKNTg0MgolJUVPRgo=",
    "type": "application/pdf"
  },
  "award": "/awards/HG011972/",
  "description": "File-format specification for native element-level MAGeCK differential quantifications from the genome-wide WTC11-EC CRISPRi screen.",
  "document_type": "file format specification",
  "lab": "/labs/jesse-engreitz/"
}


2026-09-30 08:34:47,349:iu_debug:	Success.
2026-09-30 08:34:47,350:iu_debug:	Object posted with identifier: c97c6fab-e907-4bc8-bf4b-77321e48b875

                     key                                                                                          document_alias                                      document_id                        document_uuid                                                            portal_url             document_type                                                                            attachment_download verification
        validation_guide                 jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1 /documents/dc8d807a-4038-47ac-9006-16b6e5b64a80/ dc8d807a-4038-47ac-9006-16b6e5b64a80 https://data.igvf.org/documents/dc8d807a-4038-47ac-9006-16b6e5b64a80/ file format specification                 WTC11_EC_validation_screen_guide_quantifications_file_format_specification.pdf         PASS
 validation_differential  jesse-engreitz:document_WTC11-EC_validation_screen_differential_element_quantifications_file_format_v1 /documents/babab92f-024f-4734-8a59-55134d3d86b7/ babab92f-024f-4734-8a59-55134d3d86b7 https://data.igvf.org/documents/babab92f-024f-4734-8a59-55134d3d86b7/ file format specification  WTC11_EC_validation_screen_differential_element_quantifications_file_format_specification.pdf         PASS
       genome_wide_guide                jesse-engreitz:document_WTC11-EC_genome_wide_screen_guide_quantifications_file_format_v1 /documents/1672694b-f1ff-4f48-9da1-37899a3fa49b/ 1672694b-f1ff-4f48-9da1-37899a3fa49b https://data.igvf.org/documents/1672694b-f1ff-4f48-9da1-37899a3fa49b/ file format specification                WTC11_EC_genome_wide_screen_guide_quantifications_file_format_specification.pdf         PASS
genome_wide_differential jesse-engreitz:document_WTC11-EC_genome_wide_screen_differential_element_quantifications_file_format_v1 /documents/c97c6fab-e907-4bc8-bf4b-77321e48b875/ c97c6fab-e907-4bc8-bf4b-77321e48b875 https://data.igvf.org/documents/c97c6fab-e907-4bc8-bf4b-77321e48b875/ file format specification WTC11_EC_genome_wide_screen_differential_element_quantifications_file_format_specification.pdf         PASS

Marker: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/DOCUMENTS_SUBMITTED_4_OF_4.txt

```

## Markdown cell 21

## 9. Refresh Document records and build the target-file patch plan

This read-only cell can be rerun at any time. A target analysis file that
does not yet exist is reported as `WAITING_FOR_FILE` rather than causing the
notebook to fail.

## Code cell 22

```python
document_record_rows: list[dict] = []
document_records: dict[str, dict] = {}

for key, spec in SPECIFICATIONS.items():
    record = get_with_retry(
        spec["document_alias"],
        attempts=4,
    )
    if record is None:
        raise RuntimeError(
            "Required format Document not found: "
            f"{spec['document_alias']}"
        )

    document_records[key] = record
    document_record_rows.append(
        {
            "key": key,
            "document_alias": (
                spec["document_alias"]
            ),
            "document_id": portal_record_id(
                record
            ),
            "portal_url": portal_record_url(
                record
            ),
        }
    )

document_records_frame = pd.DataFrame(
    document_record_rows
)
write_tsv(
    document_records_frame,
    WORK_DIR / "document_records_current.tsv",
)

patch_plan_rows: list[dict] = []

for key, spec in SPECIFICATIONS.items():
    lookup_identifier = (
        spec["target_file_identifier"]
        or spec["target_file_alias"]
    )

    try:
        file_record = get_with_retry(
            lookup_identifier,
            attempts=4,
        )
        lookup_error = ""
    except Exception as error:
        file_record = None
        lookup_error = (
            f"{type(error).__name__}: {error}"
        )

    document_record = document_records[key]
    document_id = portal_record_id(
        document_record
    )
    existing_specs = (
        (file_record or {}).get(
            "file_format_specifications",
            [],
        )
    )

    if lookup_error:
        action = "BLOCK_LOOKUP_ERROR"
    elif file_record is None:
        action = "WAITING_FOR_FILE"
    elif link_set_contains(
        existing_specs,
        document_id,
    ):
        action = "SKIP_ALREADY_LINKED"
    else:
        action = "PATCH_FILE"

    patch_plan_rows.append(
        {
            "key": key,
            "target_file_lookup": (
                lookup_identifier
            ),
            "target_file_alias": (
                spec["target_file_alias"]
            ),
            "action": action,
            "target_accession": str(
                (file_record or {}).get(
                    "accession",
                    "",
                )
            ),
            "target_record_id": (
                portal_record_id(file_record)
            ),
            "target_portal_url": (
                portal_record_url(file_record)
            ),
            "document_alias": (
                spec["document_alias"]
            ),
            "document_id": document_id,
            "existing_file_format_specifications": (
                ",".join(
                    sorted(
                        direct_link_tokens(
                            existing_specs
                        )
                    )
                )
            ),
            "lookup_error": lookup_error,
        }
    )

patch_plan = pd.DataFrame(patch_plan_rows)
PATCH_PLAN_PATH = (
    WORK_DIR / "file_format_patch_plan.tsv"
)
write_tsv(patch_plan, PATCH_PLAN_PATH)

print(patch_plan.to_string(index=False))
print("\nSaved:", PATCH_PLAN_PATH)

blocked = patch_plan.loc[
    patch_plan["action"].eq(
        "BLOCK_LOOKUP_ERROR"
    )
]
if not blocked.empty:
    raise RuntimeError(
        "At least one target-file lookup failed. "
        f"Review {PATCH_PLAN_PATH}"
    )
```

### Text output (stream)

```text
                     key                                                                                                 target_file_lookup                                                                                                  target_file_alias     action target_accession               target_record_id                                   target_portal_url                                                                                          document_alias                                      document_id existing_file_format_specifications lookup_error
        validation_guide                                                                                                     IGVFFI5157STFU                                       jesse-engreitz:tabular_file_WTC11-EC_validation_screen_guide_quantifications PATCH_FILE   IGVFFI5157STFU /tabular-files/IGVFFI5157STFU/ https://data.igvf.org/tabular-files/IGVFFI5157STFU/                 jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1 /documents/dc8d807a-4038-47ac-9006-16b6e5b64a80/                                                 
 validation_differential  jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications  jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications PATCH_FILE   IGVFFI3975KRIX /tabular-files/IGVFFI3975KRIX/ https://data.igvf.org/tabular-files/IGVFFI3975KRIX/  jesse-engreitz:document_WTC11-EC_validation_screen_differential_element_quantifications_file_format_v1 /documents/babab92f-024f-4734-8a59-55134d3d86b7/                                                 
       genome_wide_guide                                      jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications                                      jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications PATCH_FILE   IGVFFI8424MVZE /tabular-files/IGVFFI8424MVZE/ https://data.igvf.org/tabular-files/IGVFFI8424MVZE/                jesse-engreitz:document_WTC11-EC_genome_wide_screen_guide_quantifications_file_format_v1 /documents/1672694b-f1ff-4f48-9da1-37899a3fa49b/                                                 
genome_wide_differential jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications PATCH_FILE   IGVFFI9001HIYI /tabular-files/IGVFFI9001HIYI/ https://data.igvf.org/tabular-files/IGVFFI9001HIYI/ jesse-engreitz:document_WTC11-EC_genome_wide_screen_differential_element_quantifications_file_format_v1 /documents/c97c6fab-e907-4bc8-bf4b-77321e48b875/                                                 

Saved: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/file_format_patch_plan.tsv

```

## Markdown cell 23

## 10. Create one-row `TabularFile` PATCH TSVs

## Code cell 24

```python
patch_operation_rows: list[dict] = []

for plan in patch_plan.to_dict("records"):
    key = plan["key"]
    action = plan["action"]

    if action == "PATCH_FILE":
        record_id = (
            plan["target_accession"]
            or plan["target_record_id"]
            or plan["target_file_lookup"]
        )
        patch_tsv = (
            ROW_DIR
            / f"PATCH_file_format_{safe_stem(key)}.tsv"
        )
        payload = {
            "record_id": record_id,
            "file_format_specifications": (
                plan["document_id"]
            ),
        }
        write_tsv(
            pd.DataFrame([payload]),
            patch_tsv,
        )
    else:
        patch_tsv = None

    patch_operation_rows.append(
        {
            "key": key,
            "action": action,
            "input_tsv": (
                str(patch_tsv)
                if patch_tsv is not None
                else ""
            ),
        }
    )

patch_operations = pd.DataFrame(
    patch_operation_rows
)
PATCH_OPERATIONS_PATH = (
    WORK_DIR / "file_format_patch_operations.tsv"
)
write_tsv(
    patch_operations,
    PATCH_OPERATIONS_PATH,
)

print(patch_operations.to_string(
    index=False
))
```

### Text output (stream)

```text
                     key     action                                                                                                                                               input_tsv
        validation_guide PATCH_FILE         /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_validation_guide.tsv
 validation_differential PATCH_FILE  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_validation_differential.tsv
       genome_wide_guide PATCH_FILE        /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_genome_wide_guide.tsv
genome_wide_differential PATCH_FILE /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_genome_wide_differential.tsv

```

## Markdown cell 25

## 11. Dry-run available file patches

The default `igvf_utils` PATCH behavior extends array values and removes
duplicates. It does not overwrite other `file_format_specifications`
entries.

## Code cell 26

```python
# Compatibility fix for v2:
# Document operations use input_file; File PATCH operations use input_tsv.
if "input_tsv" not in patch_operations.columns:
    raise RuntimeError(
        "patch_operations does not contain the expected input_tsv column. "
        f"Available columns: {patch_operations.columns.tolist()}"
    )

patch_operations["input_file"] = patch_operations["input_tsv"]

# Confirm that every PATCH_FILE row points to an existing TSV.
path_check_rows = []

for operation in patch_operations.to_dict("records"):
    input_tsv = str(operation.get("input_tsv", "")).strip()
    requires_file = operation["action"] == "PATCH_FILE"
    exists = Path(input_tsv).is_file() if input_tsv else False

    path_check_rows.append(
        {
            "key": operation["key"],
            "action": operation["action"],
            "input_tsv": input_tsv,
            "requires_file": requires_file,
            "input_exists": exists,
            "status": (
                "PASS"
                if not requires_file or exists
                else "FAIL"
            ),
        }
    )

patch_path_preflight = pd.DataFrame(path_check_rows)

print(patch_path_preflight.to_string(index=False))

if not patch_path_preflight["status"].eq("PASS").all():
    raise RuntimeError(
        "At least one PATCH operation does not have an existing TSV."
    )
```

### Text output (stream)

```text
                     key     action                                                                                                                                               input_tsv  requires_file  input_exists status
        validation_guide PATCH_FILE         /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_validation_guide.tsv           True          True   PASS
 validation_differential PATCH_FILE  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_validation_differential.tsv           True          True   PASS
       genome_wide_guide PATCH_FILE        /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_genome_wide_guide.tsv           True          True   PASS
genome_wide_differential PATCH_FILE /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_genome_wide_differential.tsv           True          True   PASS

```

## Code cell 27

```python
RUN_PATCH_DRY_RUNS = True

PATCH_DRY_RUN_MARKER = (
    WORK_DIR / "PATCH_DRY_RUN_PASSED.json"
)

if not RUN_PATCH_DRY_RUNS:
    print(
        "Patch dry runs are disabled. "
        "Set RUN_PATCH_DRY_RUNS=True "
        "and rerun this cell."
    )
else:
    results: list[dict] = []

    for operation in (
        patch_operations.to_dict("records")
    ):
        if operation["action"] in {
            "SKIP_ALREADY_LINKED",
            "WAITING_FOR_FILE",
        }:
            results.append(
                {
                    "key": operation["key"],
                    "action": operation["action"],
                    "success": True,
                    "log_path": "",
                }
            )
            continue

        result = run_iu_register(
            profile="tabular_file",
            input_tsv=Path(
                operation["input_file"]
            ),
            dry_run=True,
            patch=True,
            label=(
                f"file_format_patch_"
                f"{operation['key']}"
            ),
        )
        results.append(
            {
                "key": operation["key"],
                "action": operation["action"],
                "success": result["success"],
                "log_path": result["log_path"],
            }
        )

    result_frame = pd.DataFrame(results)
    PATCH_DRY_RUN_RESULTS_PATH = (
        WORK_DIR
        / "file_format_patch_dry_run_results.tsv"
    )
    write_tsv(
        result_frame,
        PATCH_DRY_RUN_RESULTS_PATH,
    )

    if not result_frame["success"].all():
        raise RuntimeError(
            "At least one patch dry run failed"
        )

    marker = {
        "timestamp": time.strftime(
            "%Y-%m-%dT%H:%M:%S%z"
        ),
        "patch_plan_sha256": sha256sum(
            PATCH_PLAN_PATH
        ),
        "operation_files": {
            row["key"]: (
                sha256sum(Path(row["input_tsv"]))
                if row["input_tsv"]
                else ""
            )
            for row in patch_operation_rows
        },
    }
    PATCH_DRY_RUN_MARKER.write_text(
        json.dumps(marker, indent=2)
        + "\n",
        encoding="utf-8",
    )

    print("All available patch dry runs passed.")
    print("Marker:", PATCH_DRY_RUN_MARKER)
```

### Text output (stream)

```text

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_validation_guide.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_validation_guide_dryrun_20260930_084257.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:43:10,365:iu_debug:	DRY RUN is enabled.
2026-09-30 08:43:10,366:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:43:10,366:iu_debug:	submission=True: In submission mode.
2026-09-30 08:43:11,960:iu_debug:	
IN patch()
2026-09-30 08:43:11,960:iu_debug:	>>>>>>GET IGVFFI5157STFU From DACC with URL https://api.data.igvf.org/IGVFFI5157STFU/?format=json&datastore=database&frame=edit
2026-09-30 08:43:12,173:iu_debug:	<<<<<< PATCHING IGVFFI5157STFU To IGVF database with URL https://api.data.igvf.org/IGVFFI5157STFU and this payload:

{
  "file_format_specifications": [
    "dc8d807a-4038-47ac-9006-16b6e5b64a80"
  ]
}


2026-09-30 08:43:12,173:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_validation_differential.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_validation_differential_dryrun_20260930_084312.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:43:12,740:iu_debug:	DRY RUN is enabled.
2026-09-30 08:43:12,740:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:43:12,741:iu_debug:	submission=True: In submission mode.
2026-09-30 08:43:14,188:iu_debug:	
IN patch()
2026-09-30 08:43:14,188:iu_debug:	>>>>>>GET IGVFFI3975KRIX From DACC with URL https://api.data.igvf.org/IGVFFI3975KRIX/?format=json&datastore=database&frame=edit
2026-09-30 08:43:14,353:iu_debug:	<<<<<< PATCHING IGVFFI3975KRIX To IGVF database with URL https://api.data.igvf.org/IGVFFI3975KRIX and this payload:

{
  "file_format_specifications": [
    "babab92f-024f-4734-8a59-55134d3d86b7"
  ]
}


2026-09-30 08:43:14,353:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_genome_wide_guide.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_genome_wide_guide_dryrun_20260930_084314.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:43:14,915:iu_debug:	DRY RUN is enabled.
2026-09-30 08:43:14,915:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:43:14,915:iu_debug:	submission=True: In submission mode.
2026-09-30 08:43:16,314:iu_debug:	
IN patch()
2026-09-30 08:43:16,315:iu_debug:	>>>>>>GET IGVFFI8424MVZE From DACC with URL https://api.data.igvf.org/IGVFFI8424MVZE/?format=json&datastore=database&frame=edit
2026-09-30 08:43:16,464:iu_debug:	<<<<<< PATCHING IGVFFI8424MVZE To IGVF database with URL https://api.data.igvf.org/IGVFFI8424MVZE and this payload:

{
  "file_format_specifications": [
    "1672694b-f1ff-4f48-9da1-37899a3fa49b"
  ]
}


2026-09-30 08:43:16,464:iu_debug:	DRY RUN is enabled.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 -d --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_genome_wide_differential.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_genome_wide_differential_dryrun_20260930_084316.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:43:17,018:iu_debug:	DRY RUN is enabled.
2026-09-30 08:43:17,018:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:43:17,018:iu_debug:	submission=True: In submission mode.
2026-09-30 08:43:18,248:iu_debug:	
IN patch()
2026-09-30 08:43:18,248:iu_debug:	>>>>>>GET IGVFFI9001HIYI From DACC with URL https://api.data.igvf.org/IGVFFI9001HIYI/?format=json&datastore=database&frame=edit
2026-09-30 08:43:18,399:iu_debug:	<<<<<< PATCHING IGVFFI9001HIYI To IGVF database with URL https://api.data.igvf.org/IGVFFI9001HIYI and this payload:

{
  "file_format_specifications": [
    "c97c6fab-e907-4bc8-bf4b-77321e48b875"
  ]
}


2026-09-30 08:43:18,399:iu_debug:	DRY RUN is enabled.

All available patch dry runs passed.
Marker: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/PATCH_DRY_RUN_PASSED.json

```

## Markdown cell 28

## 12. Apply the available file patches

`WAITING_FOR_FILE` entries are not errors. Rerun Sections 9-12 after those
analysis files are posted.

## Code cell 29

```python
PATCH_PRODUCTION_CONFIRMATION = (
    "ATTACH_4_FORMAT_DOCUMENTS_TO_ANALYSIS_FILES"
)
REQUIRED_PATCH_CONFIRMATION = (
    "ATTACH_4_FORMAT_DOCUMENTS_TO_ANALYSIS_FILES"
)

if (
    PATCH_PRODUCTION_CONFIRMATION
    != REQUIRED_PATCH_CONFIRMATION
):
    print("Production PATCH is blocked.")
    print("Required confirmation token:")
    print(REQUIRED_PATCH_CONFIRMATION)
else:
    if not PATCH_DRY_RUN_MARKER.is_file():
        raise RuntimeError(
            "Patch dry-run marker is missing"
        )

    marker = json.loads(
        PATCH_DRY_RUN_MARKER.read_text(
            encoding="utf-8"
        )
    )
    if (
        marker["patch_plan_sha256"]
        != sha256sum(PATCH_PLAN_PATH)
    ):
        raise RuntimeError(
            "The patch plan changed after dry run"
        )

    production_rows: list[dict] = []

    for operation in (
        patch_operations.to_dict("records")
    ):
        if operation["action"] in {
            "SKIP_ALREADY_LINKED",
            "WAITING_FOR_FILE",
        }:
            production_rows.append(
                {
                    "key": operation["key"],
                    "action": operation["action"],
                    "success": True,
                    "log_path": "",
                }
            )
            continue

        input_tsv = Path(
            operation["input_file"]
        )
        if (
            marker["operation_files"][
                operation["key"]
            ]
            != sha256sum(input_tsv)
        ):
            raise RuntimeError(
                "A PATCH TSV changed after dry run: "
                f"{input_tsv}"
            )

        result = run_iu_register(
            profile="tabular_file",
            input_tsv=input_tsv,
            dry_run=False,
            patch=True,
            label=(
                f"file_format_patch_"
                f"{operation['key']}"
            ),
        )
        production_rows.append(
            {
                "key": operation["key"],
                "action": operation["action"],
                "success": result["success"],
                "log_path": result["log_path"],
            }
        )

    production_frame = pd.DataFrame(
        production_rows
    )
    PATCH_PRODUCTION_RESULTS_PATH = (
        WORK_DIR
        / "file_format_patch_production_results.tsv"
    )
    write_tsv(
        production_frame,
        PATCH_PRODUCTION_RESULTS_PATH,
    )

    print(production_frame.to_string(
        index=False
    ))
    print(
        "\nProduction results:",
        PATCH_PRODUCTION_RESULTS_PATH,
    )
```

### Text output (stream)

```text

Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_validation_guide.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_validation_guide_submit_20260930_084318.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:43:18,993:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:43:18,994:iu_debug:	submission=True: In submission mode.
2026-09-30 08:43:20,573:iu_debug:	
IN patch()
2026-09-30 08:43:20,573:iu_debug:	>>>>>>GET IGVFFI5157STFU From DACC with URL https://api.data.igvf.org/IGVFFI5157STFU/?format=json&datastore=database&frame=edit
2026-09-30 08:43:20,731:iu_debug:	<<<<<< PATCHING IGVFFI5157STFU To IGVF database with URL https://api.data.igvf.org/IGVFFI5157STFU and this payload:

{
  "file_format_specifications": [
    "dc8d807a-4038-47ac-9006-16b6e5b64a80"
  ]
}


2026-09-30 08:43:21,086:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_validation_differential.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_validation_differential_submit_20260930_084321.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:43:21,624:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:43:21,624:iu_debug:	submission=True: In submission mode.
2026-09-30 08:43:23,084:iu_debug:	
IN patch()
2026-09-30 08:43:23,084:iu_debug:	>>>>>>GET IGVFFI3975KRIX From DACC with URL https://api.data.igvf.org/IGVFFI3975KRIX/?format=json&datastore=database&frame=edit
2026-09-30 08:43:23,243:iu_debug:	<<<<<< PATCHING IGVFFI3975KRIX To IGVF database with URL https://api.data.igvf.org/IGVFFI3975KRIX and this payload:

{
  "file_format_specifications": [
    "babab92f-024f-4734-8a59-55134d3d86b7"
  ]
}


2026-09-30 08:43:23,632:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_genome_wide_guide.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_genome_wide_guide_submit_20260930_084323.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:43:24,195:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:43:24,195:iu_debug:	submission=True: In submission mode.
2026-09-30 08:43:25,683:iu_debug:	
IN patch()
2026-09-30 08:43:25,683:iu_debug:	>>>>>>GET IGVFFI8424MVZE From DACC with URL https://api.data.igvf.org/IGVFFI8424MVZE/?format=json&datastore=database&frame=edit
2026-09-30 08:43:25,841:iu_debug:	<<<<<< PATCHING IGVFFI8424MVZE To IGVF database with URL https://api.data.igvf.org/IGVFFI8424MVZE and this payload:

{
  "file_format_specifications": [
    "1672694b-f1ff-4f48-9da1-37899a3fa49b"
  ]
}


2026-09-30 08:43:26,252:iu_debug:	Success.


Command: /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/bin/python /oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py -m prod -p tabular_file --tries 3 --delay 10 --backoff 2 --patch -i /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/one_row_tsv/PATCH_file_format_genome_wide_differential.tsv
Log: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_genome_wide_differential_submit_20260930_084326.log
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/auth/__init__.py:54: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/oauth2/__init__.py:40: FutureWarning: You are using a Python version 3.9 past its end of life. Google will update google-auth with critical bug fixes on a best-effort basis, but not with any other fixes or features. Please upgrade your Python version, and then update google-auth.
  warnings.warn(eol_message.format("3.9"), FutureWarning)
/oak/stanford/groups/engreitz/Users/tri/envs/igvf_utils/lib/python3.9/site-packages/google/api_core/_python_version_support.py:242: FutureWarning: You are using a non-supported Python version (3.9.23). Google will not post any further updates to google.api_core supporting this Python version. Please upgrade to the latest Python version, or at least Python 3.10, and then update google.api_core.
  warnings.warn(message, FutureWarning)
Local version:	3.1.1
Remote version:	3.1.1
2026-09-30 08:43:26,825:iu_debug:	submission=False: In non-submission mode.
2026-09-30 08:43:26,825:iu_debug:	submission=True: In submission mode.
2026-09-30 08:43:28,368:iu_debug:	
IN patch()
2026-09-30 08:43:28,369:iu_debug:	>>>>>>GET IGVFFI9001HIYI From DACC with URL https://api.data.igvf.org/IGVFFI9001HIYI/?format=json&datastore=database&frame=edit
2026-09-30 08:43:28,518:iu_debug:	<<<<<< PATCHING IGVFFI9001HIYI To IGVF database with URL https://api.data.igvf.org/IGVFFI9001HIYI and this payload:

{
  "file_format_specifications": [
    "c97c6fab-e907-4bc8-bf4b-77321e48b875"
  ]
}


2026-09-30 08:43:28,807:iu_debug:	Success.

                     key     action  success                                                                                                                                                                log_path
        validation_guide PATCH_FILE     True         /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_validation_guide_submit_20260930_084318.log
 validation_differential PATCH_FILE     True  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_validation_differential_submit_20260930_084321.log
       genome_wide_guide PATCH_FILE     True        /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_genome_wide_guide_submit_20260930_084323.log
genome_wide_differential PATCH_FILE     True /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs/file_format_patch_genome_wide_differential_submit_20260930_084326.log

Production results: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/file_format_patch_production_results.tsv

```

## Markdown cell 30

## 13. Final read-only reconciliation

Rerun this cell after any POST or PATCH attempt. It never changes the
Portal.

## Code cell 31

```python
final_rows: list[dict] = []

for key, spec in SPECIFICATIONS.items():
    document_record = get_with_retry(
        spec["document_alias"],
        attempts=4,
    )
    target_lookup = (
        spec["target_file_identifier"]
        or spec["target_file_alias"]
    )
    target_record = get_with_retry(
        target_lookup,
        attempts=4,
    )

    document_id = portal_record_id(
        document_record
    )

    if document_record is None:
        status = "DOCUMENT_MISSING"
    elif target_record is None:
        status = "WAITING_FOR_FILE"
    elif link_set_contains(
        target_record.get(
            "file_format_specifications",
            [],
        ),
        document_id,
    ):
        status = "LINKED"
    else:
        status = "PATCH_NEEDED"

    final_rows.append(
        {
            "key": key,
            "document_alias": (
                spec["document_alias"]
            ),
            "document_id": document_id,
            "document_portal_url": (
                portal_record_url(
                    document_record
                )
            ),
            "target_file_lookup": target_lookup,
            "target_accession": str(
                (target_record or {}).get(
                    "accession",
                    "",
                )
            ),
            "target_portal_url": (
                portal_record_url(
                    target_record
                )
            ),
            "portal_file_format_specifications": (
                ",".join(
                    sorted(
                        direct_link_tokens(
                            (target_record or {}).get(
                                "file_format_specifications",
                                [],
                            )
                        )
                    )
                )
            ),
            "status": status,
        }
    )

final_reconciliation = pd.DataFrame(
    final_rows
)
FINAL_RESULTS_PATH = (
    WORK_DIR
    / "file_format_document_reconciliation.csv"
)
final_reconciliation.to_csv(
    FINAL_RESULTS_PATH,
    index=False,
    lineterminator="\n",
)

print(final_reconciliation.to_string(
    index=False
))
print("\nSaved:", FINAL_RESULTS_PATH)

linked_count = int(
    final_reconciliation["status"]
    .eq("LINKED")
    .sum()
)
print(
    f"\nLinked analysis files: "
    f"{linked_count}/4"
)

if linked_count == 4:
    complete_marker = (
        WORK_DIR
        / "FILE_FORMAT_DOCUMENTS_LINKED_4_OF_4.txt"
    )
    complete_marker.write_text(
        "All four analysis files link to their "
        "file-format specification Documents.\n"
        f"Results: {FINAL_RESULTS_PATH}\n",
        encoding="utf-8",
    )
    print("Marker:", complete_marker)
```

### Text output (stream)

```text
                     key                                                                                          document_alias                                      document_id                                                   document_portal_url                                                                                                 target_file_lookup target_accession                                   target_portal_url    portal_file_format_specifications status
        validation_guide                 jesse-engreitz:document_WTC11-EC_validation_screen_guide_quantifications_file_format_v1 /documents/dc8d807a-4038-47ac-9006-16b6e5b64a80/ https://data.igvf.org/documents/dc8d807a-4038-47ac-9006-16b6e5b64a80/                                                                                                     IGVFFI5157STFU   IGVFFI5157STFU https://data.igvf.org/tabular-files/IGVFFI5157STFU/ dc8d807a-4038-47ac-9006-16b6e5b64a80 LINKED
 validation_differential  jesse-engreitz:document_WTC11-EC_validation_screen_differential_element_quantifications_file_format_v1 /documents/babab92f-024f-4734-8a59-55134d3d86b7/ https://data.igvf.org/documents/babab92f-024f-4734-8a59-55134d3d86b7/  jesse-engreitz:tabular_file_WTC11-EC_validation_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications   IGVFFI3975KRIX https://data.igvf.org/tabular-files/IGVFFI3975KRIX/ babab92f-024f-4734-8a59-55134d3d86b7 LINKED
       genome_wide_guide                jesse-engreitz:document_WTC11-EC_genome_wide_screen_guide_quantifications_file_format_v1 /documents/1672694b-f1ff-4f48-9da1-37899a3fa49b/ https://data.igvf.org/documents/1672694b-f1ff-4f48-9da1-37899a3fa49b/                                      jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_guide_quantifications   IGVFFI8424MVZE https://data.igvf.org/tabular-files/IGVFFI8424MVZE/ 1672694b-f1ff-4f48-9da1-37899a3fa49b LINKED
genome_wide_differential jesse-engreitz:document_WTC11-EC_genome_wide_screen_differential_element_quantifications_file_format_v1 /documents/c97c6fab-e907-4bc8-bf4b-77321e48b875/ https://data.igvf.org/documents/c97c6fab-e907-4bc8-bf4b-77321e48b875/ jesse-engreitz:tabular_file_WTC11-EC_genome_wide_screen_D3_CDH5pos_vs_CDH5neg_differential_element_quantifications   IGVFFI9001HIYI https://data.igvf.org/tabular-files/IGVFFI9001HIYI/ c97c6fab-e907-4bc8-bf4b-77321e48b875 LINKED

Saved: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/file_format_document_reconciliation.csv

Linked analysis files: 4/4
Marker: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/FILE_FORMAT_DOCUMENTS_LINKED_4_OF_4.txt

```

## Markdown cell 32

## 14. Output summary

## Code cell 33

```python
output_candidates = [
    ("PDF inventory", PDF_INVENTORY_PATH),
    ("Document metadata JSON", DOCUMENT_METADATA_PATH),
        ("Document review-only TSV", DOCUMENT_REVIEW_TSV_PATH),
        ("Document JSON preflight", WORK_DIR / "document_json_preflight.tsv"),
    ("Document plan", DOCUMENT_PLAN_PATH),
    (
        "Document operations",
        DOCUMENT_OPERATIONS_PATH,
    ),
    (
        "Document dry-run marker",
        DOCUMENT_DRY_RUN_MARKER,
    ),
    (
        "Document Portal records",
        WORK_DIR / "document_portal_records.csv",
    ),
    ("File patch plan", PATCH_PLAN_PATH),
    ("File patch operations", PATCH_OPERATIONS_PATH),
    (
        "Patch dry-run marker",
        PATCH_DRY_RUN_MARKER,
    ),
    (
        "Final reconciliation",
        WORK_DIR
        / "file_format_document_reconciliation.csv",
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
        for label, path in output_candidates
    ]
)

print(output_summary.to_string(index=False))
print("\nBase directory:", BASE_DIR)
print("Logs:", LOG_DIR)
```

### Text output (stream)

```text
                  output                                                                                                                                             path  exists  size_bytes
           PDF inventory                        /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/pdf_attachment_inventory.tsv    True        1894
  Document metadata JSON            /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/document_file_format_specifications.json    True        2409
Document review-only TSV /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/document_file_format_specifications_review_only.tsv    True        1830
 Document JSON preflight                         /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/document_json_preflight.tsv    True         859
           Document plan                        /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/document_submission_plan.tsv    True         648
     Document operations                  /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/document_submission_operations.tsv    True        1154
 Document dry-run marker                        /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/DOCUMENT_DRY_RUN_PASSED.json    True         550
 Document Portal records                         /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/document_portal_records.csv    True        1667
         File patch plan                          /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/file_format_patch_plan.tsv    True        1990
   File patch operations                    /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/file_format_patch_operations.tsv    True         741
    Patch dry-run marker                           /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/PATCH_DRY_RUN_PASSED.json    True         548
    Final reconciliation             /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/file_format_document_reconciliation.csv    True        1856

Base directory: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation
Logs: /oak/stanford/groups/engreitz/Users/tri/tmp/igvf_submit_validation/_igvf_document_submission/logs

```

## Code cell 34

```python

```

