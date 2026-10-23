#!/usr/bin/env python3
"""Combine validation-array results into an IGVF sequence_file TSV and one-row TSV shards."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inventory", required=True, type=Path)
    p.add_argument("--results-dir", required=True, type=Path)
    p.add_argument("--work-dir", required=True, type=Path)
    p.add_argument("--sequencing-platform", required=True)
    p.add_argument("--sequencing-kit", default="")
    p.add_argument("--sequencing-run", type=int, default=1)
    p.add_argument("--alias-lab", default="jesse-engreitz")
    return p.parse_args()


def read_csv(path: Path) -> List[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, fieldnames: List[str], rows: Iterable[dict], delimiter: str = ",") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, delimiter=delimiter, extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def file_set_path(row: dict) -> str:
    accession = row["file_set_accession"]
    if row["file_set_type"] == "measurement_set":
        return f"/measurement-sets/{accession}/"
    if row["file_set_type"] == "auxiliary_set":
        return f"/auxiliary-sets/{accession}/"
    raise ValueError(f"Unsupported file_set_type: {row['file_set_type']}")


def make_alias(row: dict, lab: str, sequencing_run: int) -> str:
    lane = f"_L{int(row['lane']):03d}" if row.get("lane") else ""
    digest = hashlib.sha1(row["relative_path"].encode("utf-8")).hexdigest()[:8]
    modality = row["modality"].upper()
    identifier = (
        f"sequence_file_TF_CC-Perturb-seq_{row['sample_key']}_{modality}"
        f"_run{sequencing_run}{lane}_{row['read_type']}_{digest}"
    )
    identifier = re.sub(r"[^A-Za-z0-9_.-]+", "_", identifier).strip("_")
    return f"{lab}:{identifier}"


def main() -> int:
    args = parse_args()
    work_dir = args.work_dir.resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    platform = args.sequencing_platform.strip()
    if "REPLACE_ME" in platform or not platform.startswith("/platform-terms/") or not platform.endswith("/"):
        print(
            "ERROR: --sequencing-platform must be a real IGVF PlatformTerm path, e.g. /platform-terms/EFO_0022840/",
            file=sys.stderr,
        )
        return 2
    if args.sequencing_run < 1:
        print("ERROR: sequencing run must be >= 1", file=sys.stderr)
        return 2

    inventory = read_csv(args.inventory)
    if not inventory:
        print("ERROR: inventory is empty", file=sys.stderr)
        return 2

    errors: List[str] = []
    validation: Dict[int, dict] = {}
    for row in inventory:
        idx = int(row["task_index"])
        result_path = args.results_dir / f"{idx:06d}.json"
        if not result_path.is_file():
            errors.append(f"Missing validation result: {result_path}")
            continue
        try:
            result = json.loads(result_path.read_text())
        except Exception as exc:
            errors.append(f"Unreadable validation result {result_path}: {exc}")
            continue
        if result.get("status") != "OK":
            errors.append(f"Validation failed for task {idx}: {result.get('error', 'unknown error')}")
            continue
        if result.get("path") != row["path"]:
            errors.append(f"Validation path mismatch for task {idx}")
            continue
        validation[idx] = result

    md5_counts = Counter(r.get("md5sum") for r in validation.values())
    duplicate_md5 = sorted(md5 for md5, count in md5_counts.items() if md5 and count > 1)
    if duplicate_md5:
        errors.append(
            "Duplicate compressed-file MD5 values detected (IGVF identifies files by MD5): "
            + ", ".join(duplicate_md5)
        )

    report_path = work_dir / "finalization_report.txt"
    if errors:
        report_path.write_text(
            "TF sequence_file finalization FAILED\n\n" + "\n".join(f"- {e}" for e in errors) + "\n"
        )
        print(report_path.read_text(), file=sys.stderr)
        return 2

    header = [
        "aliases",
        "award",
        "lab",
        "md5sum",
        "file_format",
        "file_set",
        "content_type",
        "sequencing_platform",
    ]
    if args.sequencing_kit.strip():
        header.append("sequencing_kit")
    header += [
        "sequencing_run",
        "controlled_access",
        "submitted_file_name",
        "illumina_read_type",
        "description",
    ]

    output_rows: List[dict] = []
    manifest_rows: List[dict] = []
    aliases_seen = set()
    for row in sorted(inventory, key=lambda x: int(x["task_index"])):
        idx = int(row["task_index"])
        result = validation[idx]
        alias = make_alias(row, args.alias_lab, args.sequencing_run)
        if alias in aliases_seen:
            errors.append(f"Duplicate generated alias: {alias}")
        aliases_seen.add(alias)
        description = (
            f"TF CC-Perturb-seq {row['modality'].upper()} FASTQ for {row['sample_key']}; "
            f"sequencing run {args.sequencing_run}, {row['read_type']}; source file {row['basename']}."
        )
        out = {
            "aliases": alias,
            "award": row["award"],
            "lab": row["lab"],
            "md5sum": result["md5sum"],
            "file_format": "fastq",
            "file_set": file_set_path(row),
            "content_type": "reads",
            "sequencing_platform": platform,
            "sequencing_run": str(args.sequencing_run),
            "controlled_access": "FALSE",
            "submitted_file_name": row["path"],
            "illumina_read_type": row["read_type"],
            "description": description,
        }
        if args.sequencing_kit.strip():
            out["sequencing_kit"] = args.sequencing_kit.strip()
        output_rows.append(out)
        manifest_rows.append(
            {
                "task_index": idx,
                "alias": alias,
                "sample_key": row["sample_key"],
                "read_type": row["read_type"],
                "file_set": out["file_set"],
                "md5sum": result["md5sum"],
                "submitted_file_name": row["path"],
            }
        )

    if errors:
        report_path.write_text(
            "TF sequence_file finalization FAILED\n\n" + "\n".join(f"- {e}" for e in errors) + "\n"
        )
        print(report_path.read_text(), file=sys.stderr)
        return 2

    sequence_tsv = work_dir / "sequence_file.tsv"
    write_csv(sequence_tsv, header, output_rows, delimiter="\t")

    rows_dir = work_dir / "rows"
    if rows_dir.exists():
        shutil.rmtree(rows_dir)
    rows_dir.mkdir(parents=True)
    for manifest, out in zip(manifest_rows, output_rows):
        write_csv(rows_dir / f"{manifest['task_index']:06d}.tsv", header, [out], delimiter="\t")

    manifest_fields = [
        "task_index",
        "alias",
        "sample_key",
        "read_type",
        "file_set",
        "md5sum",
        "submitted_file_name",
    ]
    write_csv(work_dir / "upload_manifest.csv", manifest_fields, manifest_rows)

    report = [
        "TF sequence_file finalization PASSED",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        f"Rows: {len(output_rows)}",
        f"Unique aliases: {len(aliases_seen)}",
        f"Unique MD5 values: {len(md5_counts)}",
        f"Sequencing platform: {platform}",
        f"Sequencing kit: {args.sequencing_kit.strip() or '(omitted)'}",
        f"Sequencing run: {args.sequencing_run}",
        f"TSV: {sequence_tsv}",
        f"Row shards: {rows_dir}",
    ]
    report_path.write_text("\n".join(report) + "\n")
    print("\n".join(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
