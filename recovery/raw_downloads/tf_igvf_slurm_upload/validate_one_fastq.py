#!/usr/bin/env python3
"""Validate one FASTQ from an inventory row and compute its compressed-file MD5."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inventory", required=True, type=Path)
    p.add_argument("--task-index", required=True, type=int)
    p.add_argument("--results-dir", required=True, type=Path)
    return p.parse_args()


def get_inventory_row(path: Path, task_index: int) -> dict:
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            if int(row["task_index"]) == task_index:
                return row
    raise KeyError(f"task_index {task_index} not found in {path}")


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def check_first_record(path: Path) -> None:
    with gzip.open(path, "rb") as f:
        lines = [f.readline() for _ in range(4)]
    if any(line == b"" for line in lines):
        raise ValueError("FASTQ does not contain a complete first record")
    header, sequence, plus, quality = [x.rstrip(b"\r\n") for x in lines]
    if not header.startswith(b"@"):
        raise ValueError("first FASTQ line does not start with @")
    if not plus.startswith(b"+"):
        raise ValueError("third FASTQ line does not start with +")
    if len(sequence) != len(quality):
        raise ValueError(
            f"first-record sequence/quality lengths differ ({len(sequence)} vs {len(quality)})"
        )


def atomic_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def main() -> int:
    args = parse_args()
    out = args.results_dir / f"{args.task_index:06d}.json"
    started = datetime.now(timezone.utc).isoformat()
    result = {
        "task_index": args.task_index,
        "started_utc": started,
        "status": "ERROR",
    }

    try:
        row = get_inventory_row(args.inventory, args.task_index)
        path = Path(row["path"])
        result.update(
            {
                "path": str(path),
                "relative_path": row["relative_path"],
                "sample_key": row["sample_key"],
                "read_type": row["read_type"],
            }
        )
        if not path.is_file():
            raise FileNotFoundError(path)

        before = path.stat()
        expected_size = int(row["size_bytes"])
        expected_mtime_ns = int(row["mtime_ns"])
        if before.st_size != expected_size:
            raise ValueError(f"file size changed after preflight: {expected_size} -> {before.st_size}")
        if before.st_mtime_ns != expected_mtime_ns:
            raise ValueError("file modification time changed after preflight")

        gzip_proc = subprocess.run(
            ["gzip", "-t", "--", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if gzip_proc.returncode != 0:
            raise ValueError(f"gzip -t failed: {gzip_proc.stderr.strip()}")

        check_first_record(path)
        md5 = md5_file(path)

        after = path.stat()
        if after.st_size != before.st_size or after.st_mtime_ns != before.st_mtime_ns:
            raise ValueError("file changed while it was being validated")

        result.update(
            {
                "status": "OK",
                "md5sum": md5,
                "size_bytes": after.st_size,
                "gzip_test": "PASS",
                "first_fastq_record": "PASS",
            }
        )
        code = 0
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
        code = 2

    result["finished_utc"] = datetime.now(timezone.utc).isoformat()
    atomic_json(out, result)
    print(json.dumps(result, sort_keys=True))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
