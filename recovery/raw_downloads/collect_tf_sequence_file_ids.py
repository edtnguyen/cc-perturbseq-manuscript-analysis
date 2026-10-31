#!/usr/bin/env python3
"""
Collect and verify IGVF sequence_file accessions for the TF CC-Perturb-seq upload.

This script combines:
  1. Original successful/partially successful upload records from portal_reconciliation.csv
  2. Corrected PATCH metadata for the 55 gRNA records that were initially linked incorrectly
  3. Retry metadata for the 289 records that were missing after the first run

Then it queries the production IGVF Portal to verify the final accession, file_set,
md5sum, aliases, upload_status, and object status for each expected sequence_file.

Run inside the igvf_utils conda environment with IGVF_API_KEY and IGVF_SECRET_KEY exported.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

PORTAL_BASE = "https://data.igvf.org"


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def read_csv(path: Path, delimiter: str = ",") -> List[Dict[str, str]]:
    if not path.exists():
        die(f"Required file not found: {path}")
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def link_id(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        return str(value.get("@id") or value.get("id") or value.get("accession") or "")
    if isinstance(value, list):
        return ";".join(link_id(v) for v in value)
    return str(value)


def portal_url(obj: Dict[str, Any]) -> str:
    object_id = str(obj.get("@id") or "")
    accession = str(obj.get("accession") or "")
    if object_id:
        return PORTAL_BASE + object_id
    if accession:
        return f"{PORTAL_BASE}/sequence-files/{accession}/"
    return ""


def parse_modality(path_or_alias: str) -> str:
    text = path_or_alias.lower()
    if "tf_grna" in text or "_grna_" in text or "grna" in text:
        return "gRNA"
    if "gex" in text or "rna" in text:
        return "GEX"
    return ""


def parse_read(path_or_alias: str) -> str:
    m = re.search(r"(?:_|\b)(R[12])(?:\.|_|\b)", path_or_alias)
    return m.group(1) if m else ""


def aliases_to_string(aliases: Any) -> str:
    if isinstance(aliases, list):
        return ";".join(str(a) for a in aliases)
    if aliases is None:
        return ""
    return str(aliases)


def build_expected(original_work: Path, recovery_dir: Path) -> List[Dict[str, str]]:
    reconciliation_path = original_work / "portal_reconciliation.csv"
    if not reconciliation_path.exists():
        # Some users may have copied reconciliation into the recovery directory.
        alt = recovery_dir / "portal_reconciliation.csv"
        if alt.exists():
            reconciliation_path = alt
        else:
            die(
                "Could not find portal_reconciliation.csv in either "
                f"{original_work} or {recovery_dir}"
            )

    recon = read_csv(reconciliation_path)

    patch_by_accession: Dict[str, Dict[str, str]] = {}
    patch_path = recovery_dir / "patch_55_wrong_grna_file_sets.tsv"
    if patch_path.exists():
        for row in read_csv(patch_path, delimiter="\t"):
            patch_by_accession[row.get("record_id", "")] = row

    retry_manifest_path = recovery_dir / "retry_289_manifest.csv"
    retry_rows = read_csv(retry_manifest_path) if retry_manifest_path.exists() else []

    expected: List[Dict[str, str]] = []

    for row in recon:
        classification = row.get("classification", "")
        if classification not in {"FOUND_CORRECT", "FOUND_WRONG_FILE_SET"}:
            continue

        accession = row.get("portal_accession", "")
        patch = patch_by_accession.get(accession, {})
        expected_file_set = patch.get("file_set") or row.get("expected_file_set") or row.get("file_set", "")
        expected_alias = patch.get("aliases") or row.get("alias", "")
        submitted = row.get("submitted_file_name", "")
        modality = row.get("modality_from_filename") or parse_modality(submitted or expected_alias)

        expected.append(
            {
                "source_group": "original_found_correct" if classification == "FOUND_CORRECT" else "original_found_then_patched",
                "original_task_index": row.get("task_index", ""),
                "recovery_index": "",
                "sample_key": row.get("sample_key", ""),
                "modality": modality,
                "read_type": row.get("read_type") or parse_read(submitted or expected_alias),
                "expected_alias": expected_alias,
                "original_alias": row.get("alias", ""),
                "expected_file_set": expected_file_set,
                "expected_md5sum": row.get("md5sum", ""),
                "submitted_file_name": submitted,
                "known_accession": accession,
                "lookup_identifiers": ";".join(x for x in [accession, expected_alias, row.get("alias", ""), row.get("md5sum", "")] if x),
            }
        )

    for row in retry_rows:
        alias = row.get("aliases", "")
        submitted = row.get("submitted_file_name", "")
        expected.append(
            {
                "source_group": "retry_missing_from_first_run",
                "original_task_index": row.get("original_task_index", ""),
                "recovery_index": row.get("recovery_index", ""),
                "sample_key": row.get("sample_key", ""),
                "modality": row.get("modality") or parse_modality(submitted or alias),
                "read_type": row.get("read_type") or parse_read(submitted or alias),
                "expected_alias": alias,
                "original_alias": "",
                "expected_file_set": row.get("file_set", ""),
                "expected_md5sum": row.get("md5sum", ""),
                "submitted_file_name": submitted,
                "known_accession": "",
                "lookup_identifiers": ";".join(x for x in [alias, row.get("md5sum", "")] if x),
            }
        )

    return expected


def try_get(conn: Any, identifiers: Iterable[str]) -> Tuple[Dict[str, Any], str, str]:
    last_error = ""
    for ident in identifiers:
        ident = (ident or "").strip()
        if not ident:
            continue
        try:
            obj = conn.get(ident, database=True, ignore404=True, frame="object") or {}
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            continue
        if obj:
            return obj, ident, ""
    return {}, "", last_error


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--original-work",
        default="/oak/stanford/groups/engreitz/Users/tri/TF_IGVF_demultiplexed/igvf_gex_upload_20260625",
        type=Path,
        help="Directory containing the first upload/reconciliation outputs.",
    )
    parser.add_argument(
        "--recovery-dir",
        default="/oak/stanford/groups/engreitz/Users/tri/TF_IGVF_demultiplexed/igvf_recovery_20260625",
        type=Path,
        help="Directory containing recovery retry/PATCH outputs.",
    )
    parser.add_argument(
        "--out",
        default="/oak/stanford/groups/engreitz/Users/tri/TF_IGVF_demultiplexed/igvf_recovery_20260625/tf_sequence_file_accessions_all.csv",
        type=Path,
        help="Output CSV path.",
    )
    parser.add_argument(
        "--summary",
        default="/oak/stanford/groups/engreitz/Users/tri/TF_IGVF_demultiplexed/igvf_recovery_20260625/tf_sequence_file_accessions_all.summary.txt",
        type=Path,
        help="Output summary path.",
    )
    args = parser.parse_args()

    if not os.environ.get("IGVF_API_KEY") or not os.environ.get("IGVF_SECRET_KEY"):
        die("IGVF_API_KEY and IGVF_SECRET_KEY must be exported before running this script.")

    expected = build_expected(args.original_work, args.recovery_dir)

    # Import after the credential check so error messages are clearer.
    try:
        from igvf_utils.connection import Connection
    except Exception as exc:
        die(f"Could not import igvf_utils Connection. Activate the igvf_utils env first. Details: {exc}")

    conn = Connection("prod", submission=True, no_log_file=True)

    out_rows: List[Dict[str, str]] = []
    for i, row in enumerate(expected, start=1):
        identifiers = row.get("lookup_identifiers", "").split(";")
        obj, lookup_used, lookup_error = try_get(conn, identifiers)

        accession = str(obj.get("accession") or "") if obj else ""
        portal_id = str(obj.get("@id") or "") if obj else ""
        portal_file_set = link_id(obj.get("file_set")) if obj else ""
        portal_md5sum = str(obj.get("md5sum") or "") if obj else ""
        portal_aliases = aliases_to_string(obj.get("aliases")) if obj else ""
        upload_status = str(obj.get("upload_status") or "") if obj else ""
        object_status = str(obj.get("status") or "") if obj else ""

        checks = []
        if accession:
            checks.append("FOUND")
        else:
            checks.append("NOT_FOUND")
        if row.get("expected_md5sum") and portal_md5sum:
            checks.append("MD5_OK" if row["expected_md5sum"] == portal_md5sum else "MD5_MISMATCH")
        if row.get("expected_file_set") and portal_file_set:
            checks.append("FILE_SET_OK" if row["expected_file_set"] == portal_file_set else "FILE_SET_MISMATCH")
        if row.get("expected_alias") and portal_aliases:
            checks.append("ALIAS_PRESENT" if row["expected_alias"] in portal_aliases.split(";") else "ALIAS_NOT_SEEN")

        verification = ";".join(checks)

        out = {
            "index": str(i),
            "source_group": row.get("source_group", ""),
            "original_task_index": row.get("original_task_index", ""),
            "recovery_index": row.get("recovery_index", ""),
            "sample_key": row.get("sample_key", ""),
            "modality": row.get("modality", ""),
            "read_type": row.get("read_type", ""),
            "accession": accession,
            "portal_url": portal_url(obj) if obj else "",
            "portal_id": portal_id,
            "expected_alias": row.get("expected_alias", ""),
            "original_alias": row.get("original_alias", ""),
            "portal_aliases": portal_aliases,
            "expected_file_set": row.get("expected_file_set", ""),
            "portal_file_set": portal_file_set,
            "expected_md5sum": row.get("expected_md5sum", ""),
            "portal_md5sum": portal_md5sum,
            "submitted_file_name": row.get("submitted_file_name", ""),
            "upload_status": upload_status,
            "object_status": object_status,
            "lookup_used": lookup_used,
            "verification": verification,
            "lookup_error": lookup_error,
        }
        out_rows.append(out)

        if i % 25 == 0 or i == len(expected):
            print(f"Checked {i}/{len(expected)}")

    fields = [
        "index",
        "source_group",
        "original_task_index",
        "recovery_index",
        "sample_key",
        "modality",
        "read_type",
        "accession",
        "portal_url",
        "portal_id",
        "expected_alias",
        "original_alias",
        "portal_aliases",
        "expected_file_set",
        "portal_file_set",
        "expected_md5sum",
        "portal_md5sum",
        "submitted_file_name",
        "upload_status",
        "object_status",
        "lookup_used",
        "verification",
        "lookup_error",
    ]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    total = len(out_rows)
    found = sum(1 for r in out_rows if r["accession"])
    unique_accessions = len({r["accession"] for r in out_rows if r["accession"]})
    md5_bad = [r for r in out_rows if "MD5_MISMATCH" in r["verification"]]
    fs_bad = [r for r in out_rows if "FILE_SET_MISMATCH" in r["verification"]]
    not_found = [r for r in out_rows if not r["accession"]]
    alias_bad = [r for r in out_rows if "ALIAS_NOT_SEEN" in r["verification"]]

    by_group: Dict[str, int] = {}
    for r in out_rows:
        by_group[r["source_group"]] = by_group.get(r["source_group"], 0) + 1

    summary_lines = [
        "TF sequence_file accession collection summary",
        f"Output CSV: {args.out}",
        f"Expected records: {total}",
        f"Found records: {found}",
        f"Unique accessions: {unique_accessions}",
        f"Missing records: {len(not_found)}",
        f"MD5 mismatches: {len(md5_bad)}",
        f"File-set mismatches: {len(fs_bad)}",
        f"Expected alias not seen on Portal: {len(alias_bad)}",
        "Records by source group:",
    ]
    for group, count in sorted(by_group.items()):
        summary_lines.append(f"  {group}: {count}")

    if not_found:
        summary_lines.append("Missing examples:")
        for r in not_found[:20]:
            summary_lines.append(f"  {r['index']} {r['expected_alias']} {r['expected_md5sum']}")
    if fs_bad:
        summary_lines.append("File-set mismatch examples:")
        for r in fs_bad[:20]:
            summary_lines.append(
                f"  {r['accession']} expected={r['expected_file_set']} portal={r['portal_file_set']}"
            )
    if md5_bad:
        summary_lines.append("MD5 mismatch examples:")
        for r in md5_bad[:20]:
            summary_lines.append(
                f"  {r['accession']} expected={r['expected_md5sum']} portal={r['portal_md5sum']}"
            )

    args.summary.write_text("\n".join(summary_lines) + "\n")
    print("\n".join(summary_lines))

    if total != 416 or found != total or unique_accessions != total or md5_bad or fs_bad or not_found:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
