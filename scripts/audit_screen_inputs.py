#!/usr/bin/env python3
"""Audit fixed upstream genome-wide and validation-screen inputs on Sherlock.

This script is deliberately read-only. It checks file existence, size, checksum,
text-table headers, selected expected columns, workbook sheet names, and common
identifier issues. It does not rerun MAGeCK/U-tests or modify source files.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import zipfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence
from xml.etree import ElementTree


INPUTS = {
    "gw_mageck_pos_vs_neg": Path(
        "/oak/stanford/groups/engreitz/Users/dulguun/"
        "230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/"
        "Pos_vs_Neg/Pos_vs_Neg.gene_summary.txt"
    ),
    "gw_utest_pos_vs_neg": Path(
        "/oak/stanford/groups/engreitz/Users/dulguun/"
        "230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/"
        "Pos_vs_Neg/combined_utest_P1P2_adj_test.txt"
    ),
    "gw_mageck_d3total_vs_d0": Path(
        "/oak/stanford/groups/engreitz/Users/dulguun/"
        "230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/"
        "D0_vs_Total/D0_vs_EC_total.gene_summary.txt"
    ),
    "gw_utest_d3total_vs_d0": Path(
        "/oak/stanford/groups/engreitz/Users/dulguun/"
        "230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/"
        "D0_vs_Total/combined_utest_P1P2_adj_test.txt"
    ),
    "gw_mageck_columns_readme": Path(
        "/oak/stanford/groups/engreitz/Users/dulguun/"
        "230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/"
        "Mageck_output_README.xlsx"
    ),
    "gw_combined_guide_counts": Path(
        "/oak/stanford/groups/engreitz/Users/dulguun/"
        "230916_Genome-wide_Screen_Analysis/full_result_files_from_Chad/"
        "Files_used_for_mageck/combined_mageck_input.txt"
    ),
    "gw_plot_script": Path(
        "/oak/stanford/groups/engreitz/Users/dulguun/"
        "230916_Genome-wide_Screen_Analysis/generatePlots.R"
    ),
    "validation_main_results": Path(
        "/oak/stanford/groups/engreitz/Users/dulguun/"
        "231023_Validation_Screen_Files/d3pos_d3neg_mageckLFC.csv"
    ),
    "validation_gene_types": Path(
        "/oak/stanford/groups/engreitz/Users/dulguun/"
        "231023_Validation_Screen_Files/validation_genes_type.xlsx"
    ),
    "validation_hit_summary": Path(
        "/oak/stanford/groups/engreitz/Users/dulguun/"
        "231023_Validation_Screen_Files/ec_screen_hits.csv"
    ),
    "validation_plot_script": Path(
        "/oak/stanford/groups/engreitz/Users/dulguun/"
        "231023_Validation_Screen_Files/generatePlots.R"
    ),
    "validation_mageck_counts": Path(
        "/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/"
        "MageckAnalysis/mageck_sgRNA_readcount.txt"
    ),
    "validation_utest_counts": Path(
        "/oak/stanford/groups/engreitz/Users/chloemc/GenomeWideValidation/"
        "UTestAnalysis/guideCounts.csv"
    ),
}


EXPECTED_COLUMNS = {
    "gw_mageck_pos_vs_neg": {
        "id",
        "neg|lfc",
        "neg|p-value",
        "neg|fdr",
        "pos|lfc",
        "pos|p-value",
        "pos|fdr",
    },
    "gw_utest_pos_vs_neg": {"Target", "log2FC", "p_value_corrected"},
    "validation_main_results": {
        "Gene",
        "U.Test.FC.BR1",
        "U.Test.FC.BR2",
        "Average.U.Test.FC",
        "Original.FC.U.Test",
        "Original.FC.Mageck",
    },
}


@dataclass
class AuditRecord:
    key: str
    path: str
    exists: bool
    kind: str
    size_bytes: int | None
    modified_utc: str | None
    sha256: str | None
    hash_status: str
    readable: bool
    notes: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iso_utc(timestamp: float) -> str:
    return datetime.fromtimestamp(timestamp, tz=timezone.utc).isoformat()


def file_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".xlsx":
        return "xlsx"
    if suffix in {".csv", ".tsv", ".txt", ".r", ".py", ".sh"}:
        return "text"
    return suffix.lstrip(".") or "unknown"


def read_sample(path: Path, max_bytes: int = 128 * 1024) -> str:
    with path.open("rb") as handle:
        raw = handle.read(max_bytes)
    return raw.decode("utf-8-sig", errors="replace")


def detect_header(path: Path) -> tuple[str, list[str]]:
    """Return a delimiter label and parsed header fields for a text table."""
    sample = read_sample(path)
    lines = [line for line in sample.splitlines() if line.strip()]
    if not lines:
        return "empty", []
    first = lines[0]

    if "\t" in first:
        return "tab", next(csv.reader([first], delimiter="\t"))
    if "," in first:
        return "comma", next(csv.reader([first], delimiter=","))
    if ";" in first:
        return "semicolon", next(csv.reader([first], delimiter=";"))
    return "whitespace", re.split(r"\s+", first.strip())


def iter_first_column(path: Path, delimiter_label: str, limit: int | None = None) -> Iterable[str]:
    """Yield the first field from a text table, skipping the header."""
    yielded = 0
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        header = handle.readline()
        del header
        for line in handle:
            if not line.strip():
                continue
            if delimiter_label == "tab":
                value = line.rstrip("\n\r").split("\t", 1)[0]
            elif delimiter_label == "comma":
                value = next(csv.reader([line]))[0]
            elif delimiter_label == "semicolon":
                value = next(csv.reader([line], delimiter=";"))[0]
            else:
                value = re.split(r"\s+", line.strip(), maxsplit=1)[0]
            yield value
            yielded += 1
            if limit is not None and yielded >= limit:
                return


def xlsx_sheet_names(path: Path) -> list[str]:
    """Read workbook sheet names with only the Python standard library."""
    namespace = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        xml = archive.read("xl/workbook.xml")
    root = ElementTree.fromstring(xml)
    sheets = root.find("m:sheets", namespace)
    if sheets is None:
        return []
    return [sheet.attrib.get("name", "") for sheet in sheets]


def normalize_gene_identifier(identifier: str) -> str:
    value = identifier.strip()
    value = re.sub(r"_-P\d+(?:P\d+)?$", "", value, flags=re.IGNORECASE)
    value = re.sub(r"-P\d+(?:P\d+)?$", "", value, flags=re.IGNORECASE)
    value = value.rstrip("_")
    return value


def write_tsv(path: Path, rows: Sequence[dict[str, object]], fieldnames: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def audit_one(key: str, path: Path, hash_max_bytes: int) -> AuditRecord:
    if not path.exists():
        return AuditRecord(
            key=key,
            path=str(path),
            exists=False,
            kind=file_kind(path),
            size_bytes=None,
            modified_utc=None,
            sha256=None,
            hash_status="missing",
            readable=False,
            notes="File not found",
        )

    try:
        stat = path.stat()
        readable = os.access(path, os.R_OK)
        checksum: str | None = None
        if not path.is_file():
            hash_status = "not-a-file"
        elif stat.st_size > hash_max_bytes:
            hash_status = f"skipped-size>{hash_max_bytes}"
        elif not readable:
            hash_status = "unreadable"
        else:
            checksum = sha256_file(path)
            hash_status = "computed"
        return AuditRecord(
            key=key,
            path=str(path),
            exists=True,
            kind=file_kind(path),
            size_bytes=stat.st_size,
            modified_utc=iso_utc(stat.st_mtime),
            sha256=checksum,
            hash_status=hash_status,
            readable=readable,
            notes="",
        )
    except OSError as exc:
        return AuditRecord(
            key=key,
            path=str(path),
            exists=True,
            kind=file_kind(path),
            size_bytes=None,
            modified_utc=None,
            sha256=None,
            hash_status="error",
            readable=False,
            notes=f"{type(exc).__name__}: {exc}",
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reproduction_audit/screen_inputs"),
        help="Output directory for audit receipts.",
    )
    parser.add_argument(
        "--hash-max-mib",
        type=float,
        default=250.0,
        help="Skip SHA-256 for files larger than this many MiB (default: 250).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    outdir: Path = args.output
    outdir.mkdir(parents=True, exist_ok=True)
    hash_max_bytes = int(args.hash_max_mib * 1024 * 1024)

    records = [audit_one(key, path, hash_max_bytes) for key, path in INPUTS.items()]
    record_by_key = {record.key: record for record in records}

    header_rows: list[dict[str, object]] = []
    identifier_rows: list[dict[str, object]] = []

    for key, source_path in INPUTS.items():
        record = record_by_key[key]
        if not record.exists or not record.readable or not source_path.is_file():
            continue

        if source_path.suffix.lower() == ".xlsx":
            try:
                sheets = xlsx_sheet_names(source_path)
                header_rows.append(
                    {
                        "key": key,
                        "path": str(source_path),
                        "format": "xlsx",
                        "n_fields": len(sheets),
                        "fields_or_sheets": " | ".join(sheets),
                        "expected_columns_missing": "not-applicable",
                    }
                )
            except (OSError, KeyError, zipfile.BadZipFile, ElementTree.ParseError) as exc:
                header_rows.append(
                    {
                        "key": key,
                        "path": str(source_path),
                        "format": "xlsx-error",
                        "n_fields": 0,
                        "fields_or_sheets": "",
                        "expected_columns_missing": f"{type(exc).__name__}: {exc}",
                    }
                )
            continue

        if source_path.suffix.lower() not in {".csv", ".tsv", ".txt"}:
            continue

        try:
            delimiter_label, fields = detect_header(source_path)
            expected = EXPECTED_COLUMNS.get(key, set())
            missing = sorted(expected.difference(fields))
            header_rows.append(
                {
                    "key": key,
                    "path": str(source_path),
                    "format": delimiter_label,
                    "n_fields": len(fields),
                    "fields_or_sheets": " | ".join(fields),
                    "expected_columns_missing": " | ".join(missing),
                }
            )

            if key == "gw_mageck_pos_vs_neg":
                raw_ids = list(iter_first_column(source_path, delimiter_label))
                normalized = [normalize_gene_identifier(value) for value in raw_ids]
                nedd8_raw = [value for value, gene in zip(raw_ids, normalized) if gene.upper() == "NEDD8"]
                counts: dict[str, int] = {}
                for gene in normalized:
                    counts[gene] = counts.get(gene, 0) + 1
                duplicate_norm = sorted((gene, count) for gene, count in counts.items() if count > 1)
                identifier_rows.extend(
                    [
                        {
                            "check": "row_count",
                            "key": key,
                            "value": len(raw_ids),
                            "details": "Number of data rows",
                        },
                        {
                            "check": "unique_normalized_identifiers",
                            "key": key,
                            "value": len(counts),
                            "details": "After removing promoter suffixes and trailing underscores",
                        },
                        {
                            "check": "nedd8_raw_identifiers",
                            "key": key,
                            "value": len(nedd8_raw),
                            "details": " | ".join(nedd8_raw),
                        },
                        {
                            "check": "normalized_duplicate_count",
                            "key": key,
                            "value": len(duplicate_norm),
                            "details": " | ".join(f"{gene}:{count}" for gene, count in duplicate_norm[:50]),
                        },
                    ]
                )
        except OSError as exc:
            header_rows.append(
                {
                    "key": key,
                    "path": str(source_path),
                    "format": "read-error",
                    "n_fields": 0,
                    "fields_or_sheets": "",
                    "expected_columns_missing": f"{type(exc).__name__}: {exc}",
                }
            )

    input_rows = [asdict(record) for record in records]
    write_tsv(
        outdir / "input_audit.tsv",
        input_rows,
        [
            "key",
            "path",
            "exists",
            "kind",
            "size_bytes",
            "modified_utc",
            "sha256",
            "hash_status",
            "readable",
            "notes",
        ],
    )
    write_tsv(
        outdir / "table_headers.tsv",
        header_rows,
        [
            "key",
            "path",
            "format",
            "n_fields",
            "fields_or_sheets",
            "expected_columns_missing",
        ],
    )
    write_tsv(
        outdir / "identifier_checks.tsv",
        identifier_rows,
        ["check", "key", "value", "details"],
    )

    missing_keys = [record.key for record in records if not record.exists]
    unreadable_keys = [record.key for record in records if record.exists and not record.readable]
    missing_expected = {
        row["key"]: row["expected_columns_missing"]
        for row in header_rows
        if row["expected_columns_missing"] not in {"", "not-applicable"}
    }
    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).resolve()),
        "n_inputs": len(records),
        "n_present": sum(record.exists for record in records),
        "n_readable": sum(record.readable for record in records),
        "missing_keys": missing_keys,
        "unreadable_keys": unreadable_keys,
        "tables_with_missing_expected_columns": missing_expected,
        "hash_max_mib": args.hash_max_mib,
        "output_directory": str(outdir.resolve()),
    }
    with (outdir / "audit_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"\nWrote audit receipts to: {outdir.resolve()}")

    # Missing files are an audit finding, not an execution failure. Return nonzero
    # only when a present file is unreadable or an expected table schema changed.
    if unreadable_keys or missing_expected:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
