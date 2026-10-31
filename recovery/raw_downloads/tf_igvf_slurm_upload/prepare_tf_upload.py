#!/usr/bin/env python3
"""Scan TF CC-Perturb-seq FASTQs and match each file to spreadsheet-derived IGVF file-set metadata.

Uses only the Python standard library. It performs lightweight preflight checks; full gzip/MD5
validation is delegated to the Slurm validation array.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz")
READ_RE = re.compile(r"(?:^|[_\-.])(R[123]|I[12])(?=[_\-.]|$)", re.IGNORECASE)
LANE_RE = re.compile(r"(?:^|[_\-.])L(\d{3})(?=[_\-.]|$)", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", required=True, type=Path)
    p.add_argument("--mapping", required=True, type=Path)
    p.add_argument("--work-dir", required=True, type=Path)
    p.add_argument("--modality", choices=("gex", "grna"), default="gex")
    p.add_argument(
        "--allow-missing-samples",
        action="store_true",
        help="Do not fail if some spreadsheet sample keys have no FASTQ files.",
    )
    p.add_argument(
        "--include-index-reads",
        action="store_true",
        help="Allow I1/I2 FASTQs. By default only R1/R2 are accepted.",
    )
    return p.parse_args()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def strip_fastq_suffix(name: str) -> str:
    lower = name.lower()
    for suffix in FASTQ_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def parse_sample_key(name: str) -> Optional[str]:
    """Extract subN_group_dN_repN from a permissive filename representation."""
    s = strip_fastq_suffix(name).lower()
    sub = re.search(r"sub(?:library)?[_\-.]*0*(\d{1,2})(?!\d)", s)
    group = re.search(r"group[_\-.]*d?[_\-.]*([0-3])(?!\d)", s)
    if group is None:
        group = re.search(r"(?:^|[_\-.])d([0-3])(?=[_\-.]|$)", s)
    rep = re.search(r"rep(?:licate)?[_\-.]*0*([12])(?!\d)", s)
    if not (sub and group and rep):
        return None
    sub_n, group_n, rep_n = int(sub.group(1)), int(group.group(1)), int(rep.group(1))
    return f"sub{sub_n}_group_d{group_n}_rep{rep_n}"


def parse_read(name: str) -> Tuple[Optional[str], Optional[Tuple[int, int]]]:
    stem = strip_fastq_suffix(name)
    matches = list(READ_RE.finditer(stem))
    distinct = sorted({m.group(1).upper() for m in matches})
    if len(distinct) != 1:
        return None, None
    # Use the final exact read token, which is normally the Illumina read field.
    m = [x for x in matches if x.group(1).upper() == distinct[0]][-1]
    return distinct[0], m.span(1)


def pair_key(relative_path: str, span: Tuple[int, int]) -> str:
    p = Path(relative_path)
    stem = strip_fastq_suffix(p.name)
    start, end = span
    paired_name = stem[:start] + "{READ}" + stem[end:]
    return str(p.parent / paired_name) if str(p.parent) != "." else paired_name


def load_mapping(path: Path) -> Dict[str, dict]:
    required = {
        "sample_key",
        "virtual_sample_accession",
        "measurement_set_accession",
        "measurement_set_alias",
        "auxiliary_set_accession",
        "auxiliary_set_alias",
        "award",
        "lab",
    }
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Mapping CSV is missing columns: {sorted(missing)}")
        rows = list(reader)
    mapping: Dict[str, dict] = {}
    for row in rows:
        key = row["sample_key"].strip()
        if key in mapping:
            raise ValueError(f"Duplicate sample_key in mapping: {key}")
        mapping[key] = row
    if not mapping:
        raise ValueError("Mapping CSV has no rows")
    return mapping


def write_csv(path: Path, fieldnames: List[str], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def main() -> int:
    args = parse_args()
    data_dir = args.data_dir.expanduser().resolve()
    mapping_path = args.mapping.expanduser().resolve()
    work_dir = args.work_dir.expanduser().resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    if not data_dir.is_dir():
        print(f"ERROR: data directory does not exist: {data_dir}", file=sys.stderr)
        return 2
    if not mapping_path.is_file():
        print(f"ERROR: mapping CSV does not exist: {mapping_path}", file=sys.stderr)
        return 2

    try:
        mapping = load_mapping(mapping_path)
    except Exception as exc:
        print(f"ERROR: could not load mapping: {exc}", file=sys.stderr)
        return 2

    files = sorted(
        p.resolve()
        for p in data_dir.rglob("*")
        if p.is_file() and any(p.name.lower().endswith(s) for s in FASTQ_SUFFIXES)
    )

    audit: List[dict] = []
    inventory: List[dict] = []
    errors: List[str] = []
    ignored_count = 0
    allowed_reads = {"R1", "R2"} | ({"I1", "I2"} if args.include_index_reads else set())

    if not files:
        errors.append(f"No .fastq.gz or .fq.gz files found under {data_dir}")

    for path in files:
        rel = str(path.relative_to(data_dir))
        if "undetermined" in path.name.lower():
            ignored_count += 1
            audit.append(
                {
                    "path": str(path),
                    "relative_path": rel,
                    "sample_key": "",
                    "read_type": "",
                    "pair_key": "",
                    "file_set_accession": "",
                    "status": "IGNORED",
                    "errors": "undetermined_fastq_has_no_sample_metadata",
                }
            )
            continue
        sample_key = parse_sample_key(path.name)
        read_type, read_span = parse_read(path.name)
        row_errors: List[str] = []

        if sample_key is None:
            row_errors.append("could_not_parse_sample_key")
        elif sample_key not in mapping:
            row_errors.append("sample_key_not_in_spreadsheet_mapping")

        if read_type is None or read_span is None:
            row_errors.append("could_not_parse_unique_read_type")
        elif read_type not in allowed_reads:
            row_errors.append(f"read_type_not_allowed:{read_type}")

        stat = path.stat()
        lane_m = LANE_RE.search(strip_fastq_suffix(path.name))
        lane = int(lane_m.group(1)) if lane_m else ""

        meta = mapping.get(sample_key or "", {})
        if args.modality == "gex":
            file_set_accession = meta.get("measurement_set_accession", "")
            file_set_alias = meta.get("measurement_set_alias", "")
            file_set_type = "measurement_set"
        else:
            file_set_accession = meta.get("auxiliary_set_accession", "")
            file_set_alias = meta.get("auxiliary_set_alias", "")
            file_set_type = "auxiliary_set"

        pkey = pair_key(rel, read_span) if read_span else ""
        status = "ERROR" if row_errors else "OK"
        audit.append(
            {
                "path": str(path),
                "relative_path": rel,
                "sample_key": sample_key or "",
                "read_type": read_type or "",
                "pair_key": pkey,
                "file_set_accession": file_set_accession,
                "status": status,
                "errors": ";".join(row_errors),
            }
        )
        if row_errors:
            errors.append(f"{rel}: {', '.join(row_errors)}")
            continue

        inventory.append(
            {
                "task_index": 0,  # assigned after sorting
                "path": str(path),
                "relative_path": rel,
                "basename": path.name,
                "size_bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "sample_key": sample_key,
                "read_type": read_type,
                "pair_key": pkey,
                "lane": lane,
                "modality": args.modality,
                "file_set_type": file_set_type,
                "file_set_accession": file_set_accession,
                "file_set_alias": file_set_alias,
                "virtual_sample_accession": meta.get("virtual_sample_accession", ""),
                "award": meta.get("award", ""),
                "lab": meta.get("lab", ""),
            }
        )

    # Pair-level checks.
    grouped: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for row in inventory:
        grouped[(row["sample_key"], row["pair_key"])].append(row)
    for (sample_key, pkey), group_rows in sorted(grouped.items()):
        counts = Counter(r["read_type"] for r in group_rows)
        for required_read in ("R1", "R2"):
            if counts[required_read] != 1:
                errors.append(
                    f"Pair {sample_key} / {pkey}: expected exactly one {required_read}, found {counts[required_read]}"
                )

    used_samples = {r["sample_key"] for r in inventory}
    expected_samples = set(mapping)
    missing_samples = sorted(expected_samples - used_samples)
    if missing_samples and not args.allow_missing_samples:
        errors.append(
            f"Missing FASTQs for {len(missing_samples)} spreadsheet sample keys: "
            + ", ".join(missing_samples)
        )

    # Sort deterministically and assign 1-based Slurm indices.
    inventory.sort(key=lambda r: (r["sample_key"], r["pair_key"], r["read_type"], r["relative_path"]))
    for i, row in enumerate(inventory, start=1):
        row["task_index"] = i

    audit_fields = [
        "path",
        "relative_path",
        "sample_key",
        "read_type",
        "pair_key",
        "file_set_accession",
        "status",
        "errors",
    ]
    inventory_fields = [
        "task_index",
        "path",
        "relative_path",
        "basename",
        "size_bytes",
        "mtime_ns",
        "sample_key",
        "read_type",
        "pair_key",
        "lane",
        "modality",
        "file_set_type",
        "file_set_accession",
        "file_set_alias",
        "virtual_sample_accession",
        "award",
        "lab",
    ]
    write_csv(work_dir / "preflight_audit.csv", audit_fields, audit)
    write_csv(work_dir / "inventory.csv", inventory_fields, inventory)

    report_lines = [
        "TF CC-Perturb-seq FASTQ preflight",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        f"Data directory: {data_dir}",
        f"Mapping CSV: {mapping_path}",
        f"Mapping SHA256: {sha256_file(mapping_path)}",
        f"Modality: {args.modality}",
        f"Spreadsheet mapping rows: {len(mapping)}",
        f"FASTQ files found: {len(files)}",
        f"FASTQ files mapped: {len(inventory)}",
        f"FASTQ files ignored (Undetermined): {ignored_count}",
        f"Sample keys represented: {len(used_samples)}",
        f"Sample keys missing: {len(missing_samples)}",
        f"Preflight errors: {len(errors)}",
        "",
    ]
    if missing_samples:
        report_lines += ["Missing sample keys:", *[f"  {x}" for x in missing_samples], ""]
    if errors:
        report_lines += ["Errors:", *[f"  {x}" for x in errors]]
    else:
        report_lines += ["RESULT: PASS"]
    (work_dir / "preflight_report.txt").write_text("\n".join(report_lines) + "\n")

    print("\n".join(report_lines[:10]))
    print(f"Audit: {work_dir / 'preflight_audit.csv'}")
    print(f"Inventory: {work_dir / 'inventory.csv'}")
    print(f"Report: {work_dir / 'preflight_report.txt'}")

    if errors:
        print("Preflight failed. No Slurm jobs should be submitted.", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
