#!/usr/bin/env python3
"""Read-only audit of final gRNA demultiplexed FASTQ outputs.

The audit expects gRNA_sub1..gRNA_sub13 and these eight groups:
  d0_rep1, d0_rep2, d1_rep1, d1_rep2,
  d2_rep1, d2_rep2, d3_rep1, d3_rep2

For every expected sample/group pair it requires exactly one R1 and one R2,
then streams both files to EOF to verify:
  * gzip integrity / CRC
  * complete four-line FASTQ records
  * @ header and + separator lines
  * sequence and quality lengths agree
  * R1 and R2 read identifiers agree
  * R1 and R2 record counts agree
  * the pair is non-empty

It never edits or deletes FASTQs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

DEFAULT_ROOT = Path(
    "/oak/stanford/groups/engreitz/Users/tri/TF_IGVF_demux"
)
EXPECTED_SUBS = tuple(range(1, 14))
EXPECTED_GROUPS = (
    "d0_rep1",
    "d0_rep2",
    "d1_rep1",
    "d1_rep2",
    "d2_rep1",
    "d2_rep2",
    "d3_rep1",
    "d3_rep2",
)

SUB_RE = re.compile(r"grna[_-]?sub(\d+)", re.IGNORECASE)
GROUP_RE = re.compile(r"d([0-3])[_\-.]?rep([12])", re.IGNORECASE)
READ_RE = re.compile(r"(?:^|[^A-Za-z0-9])(R[12])(?=$|[^A-Za-z0-9])", re.IGNORECASE)


@dataclass(frozen=True)
class ParsedFastq:
    path: Path
    relpath: str
    size: int
    sub: Optional[int]
    group: Optional[str]
    read: Optional[str]


@dataclass(frozen=True)
class PairResult:
    sub: int
    group: str
    status: str
    r1: str
    r2: str
    r1_size: int
    r2_size: int
    records: int
    elapsed_sec: float
    error_code: str
    error_detail: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit final gRNA_sub1..13 demultiplexed FASTQ pairs."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help=f"Directory containing final gRNA FASTQs (default: {DEFAULT_ROOT})",
    )
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=None,
        help="Report directory (default: ./grna_audit_YYYYmmdd_HHMMSS)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of FASTQ pairs read concurrently (default: 4; keep modest on OAK)",
    )
    parser.add_argument(
        "--inventory-only",
        action="store_true",
        help="Check presence/cardinality only; do not stream FASTQ contents",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=5,
        help="Print progress after this many validated pairs (default: 5)",
    )
    return parser.parse_args()


def parse_fastq(path: Path, root: Path) -> ParsedFastq:
    rel = str(path.relative_to(root))
    searchable = rel.replace(os.sep, "_")

    sub_match = SUB_RE.search(searchable)
    group_match = GROUP_RE.search(searchable)
    read_matches = list(READ_RE.finditer(searchable))

    sub = int(sub_match.group(1)) if sub_match else None
    group = (
        f"d{group_match.group(1)}_rep{group_match.group(2)}"
        if group_match
        else None
    )
    read = read_matches[-1].group(1).upper() if read_matches else None

    try:
        size = path.stat().st_size
    except OSError:
        size = -1

    return ParsedFastq(
        path=path,
        relpath=rel,
        size=size,
        sub=sub,
        group=group,
        read=read,
    )


def normalize_read_id(header: bytes) -> bytes:
    """Normalize common paired FASTQ identifiers.

    Uses the first whitespace-delimited token, removes leading '@', and strips
    terminal /1 or /2. Illumina's second token (1:N:... vs 2:N:...) is ignored.
    """
    token = header.strip().split(None, 1)[0]
    if token.startswith(b"@"):
        token = token[1:]
    if token.endswith(b"/1") or token.endswith(b"/2"):
        token = token[:-2]
    return token


def read_record(handle: gzip.GzipFile, mate: str, record_no: int):
    header = handle.readline()
    if not header:
        return None

    seq = handle.readline()
    plus = handle.readline()
    qual = handle.readline()

    if not seq or not plus or not qual:
        raise ValueError(
            f"TRUNCATED_FASTQ: {mate} ended within record {record_no}"
        )

    header_s = header.rstrip(b"\r\n")
    seq_s = seq.rstrip(b"\r\n")
    plus_s = plus.rstrip(b"\r\n")
    qual_s = qual.rstrip(b"\r\n")

    if not header_s.startswith(b"@"):
        raise ValueError(
            f"BAD_HEADER: {mate} record {record_no} does not start with @"
        )
    if not plus_s.startswith(b"+"):
        raise ValueError(
            f"BAD_PLUS: {mate} record {record_no} separator does not start with +"
        )
    if len(seq_s) != len(qual_s):
        raise ValueError(
            f"SEQ_QUAL_LENGTH_MISMATCH: {mate} record {record_no}: "
            f"sequence={len(seq_s)}, quality={len(qual_s)}"
        )

    return header_s, seq_s, qual_s


def validate_pair(sub: int, group: str, r1: ParsedFastq, r2: ParsedFastq) -> PairResult:
    started = time.monotonic()
    records = 0

    try:
        if r1.size <= 0:
            raise ValueError(f"EMPTY_FILE: R1 compressed size is {r1.size}")
        if r2.size <= 0:
            raise ValueError(f"EMPTY_FILE: R2 compressed size is {r2.size}")

        # Reading to EOF makes Python's gzip reader verify gzip framing and CRC.
        with gzip.open(r1.path, "rb") as f1, gzip.open(r2.path, "rb") as f2:
            while True:
                record_no = records + 1
                rec1 = read_record(f1, "R1", record_no)
                rec2 = read_record(f2, "R2", record_no)

                if rec1 is None and rec2 is None:
                    break
                if rec1 is None:
                    raise ValueError(
                        f"PAIR_COUNT_MISMATCH: R1 ended before R2 at record {record_no}"
                    )
                if rec2 is None:
                    raise ValueError(
                        f"PAIR_COUNT_MISMATCH: R2 ended before R1 at record {record_no}"
                    )

                id1 = normalize_read_id(rec1[0])
                id2 = normalize_read_id(rec2[0])
                if id1 != id2:
                    raise ValueError(
                        f"READ_ID_MISMATCH: record {record_no}: "
                        f"R1={id1[:160]!r}, R2={id2[:160]!r}"
                    )

                records += 1

        if records == 0:
            raise ValueError("ZERO_RECORDS: FASTQ pair contains no reads")

        return PairResult(
            sub=sub,
            group=group,
            status="PASS",
            r1=str(r1.path),
            r2=str(r2.path),
            r1_size=r1.size,
            r2_size=r2.size,
            records=records,
            elapsed_sec=time.monotonic() - started,
            error_code="",
            error_detail="",
        )

    except (gzip.BadGzipFile, EOFError, OSError, ValueError) as exc:
        detail = str(exc)
        code = detail.split(":", 1)[0] if ":" in detail else type(exc).__name__
        if isinstance(exc, (gzip.BadGzipFile, EOFError, OSError)):
            code = "GZIP_OR_IO_ERROR"
        return PairResult(
            sub=sub,
            group=group,
            status="FAIL",
            r1=str(r1.path),
            r2=str(r2.path),
            r1_size=r1.size,
            r2_size=r2.size,
            records=records,
            elapsed_sec=time.monotonic() - started,
            error_code=code,
            error_detail=detail,
        )


def write_tsv(path: Path, header: Iterable[str], rows: Iterable[Iterable[object]]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(list(header))
        for row in rows:
            writer.writerow(list(row))
    tmp.replace(path)


def main() -> int:
    args = parse_args()
    root = args.root.resolve()

    if args.workers < 1:
        print("ERROR: --workers must be >= 1", file=sys.stderr)
        return 2
    if not root.is_dir():
        print(f"ERROR: root directory does not exist: {root}", file=sys.stderr)
        return 2

    stamp = time.strftime("%Y%m%d_%H%M%S")
    report_dir = (
        args.report_dir.resolve()
        if args.report_dir is not None
        else (Path.cwd() / f"grna_audit_{stamp}").resolve()
    )
    report_dir.mkdir(parents=True, exist_ok=False)

    print(f"Audit root:   {root}")
    print(f"Report dir:   {report_dir}")
    print(f"Workers:      {args.workers}")
    print(f"Audit mode:   {'inventory only' if args.inventory_only else 'full FASTQ validation'}")

    paths = sorted(
        (p for p in root.rglob("*.fastq.gz") if p.is_file()),
        key=lambda p: str(p),
    )
    parsed = [parse_fastq(path, root) for path in paths]

    write_tsv(
        report_dir / "inventory.tsv",
        ("relative_path", "bytes", "sublibrary", "group", "read", "classification"),
        (
            (
                item.relpath,
                item.size,
                f"gRNA_sub{item.sub}" if item.sub is not None else "",
                item.group or "",
                item.read or "",
                "PARSED"
                if item.sub is not None and item.group is not None and item.read is not None
                else "UNCLASSIFIED",
            )
            for item in parsed
        ),
    )

    by_key: dict[tuple[int, str, str], list[ParsedFastq]] = defaultdict(list)
    unclassified: list[ParsedFastq] = []
    extras: list[ParsedFastq] = []

    for item in parsed:
        if item.sub is None or item.group is None or item.read is None:
            unclassified.append(item)
            continue
        by_key[(item.sub, item.group, item.read)].append(item)
        if item.sub not in EXPECTED_SUBS:
            extras.append(item)

    issues: list[tuple[str, str, str, str, str, str]] = []

    for item in unclassified:
        issues.append(
            (
                "WARNING",
                "UNCLASSIFIED_FASTQ",
                f"gRNA_sub{item.sub}" if item.sub is not None else "",
                item.group or "",
                item.read or "",
                str(item.path),
            )
        )

    extra_subs = sorted({item.sub for item in extras if item.sub is not None})
    for sub in extra_subs:
        count = sum(1 for item in extras if item.sub == sub)
        issues.append(
            (
                "WARNING",
                "EXTRA_SUBLIBRARY",
                f"gRNA_sub{sub}",
                "",
                "",
                f"Found {count} parsed FASTQs; sublibrary is outside requested sub1..13",
            )
        )

    validatable: list[tuple[int, str, ParsedFastq, ParsedFastq]] = []
    structural_failures: list[PairResult] = []

    for sub in EXPECTED_SUBS:
        for group in EXPECTED_GROUPS:
            r1s = by_key.get((sub, group, "R1"), [])
            r2s = by_key.get((sub, group, "R2"), [])

            if len(r1s) != 1 or len(r2s) != 1:
                if len(r1s) == 0:
                    issues.append(
                        ("ERROR", "MISSING_R1", f"gRNA_sub{sub}", group, "R1", "No matching file")
                    )
                elif len(r1s) > 1:
                    issues.append(
                        (
                            "ERROR",
                            "DUPLICATE_R1",
                            f"gRNA_sub{sub}",
                            group,
                            "R1",
                            " | ".join(str(x.path) for x in r1s),
                        )
                    )

                if len(r2s) == 0:
                    issues.append(
                        ("ERROR", "MISSING_R2", f"gRNA_sub{sub}", group, "R2", "No matching file")
                    )
                elif len(r2s) > 1:
                    issues.append(
                        (
                            "ERROR",
                            "DUPLICATE_R2",
                            f"gRNA_sub{sub}",
                            group,
                            "R2",
                            " | ".join(str(x.path) for x in r2s),
                        )
                    )

                structural_failures.append(
                    PairResult(
                        sub=sub,
                        group=group,
                        status="FAIL",
                        r1=" | ".join(str(x.path) for x in r1s),
                        r2=" | ".join(str(x.path) for x in r2s),
                        r1_size=sum(max(0, x.size) for x in r1s),
                        r2_size=sum(max(0, x.size) for x in r2s),
                        records=0,
                        elapsed_sec=0.0,
                        error_code="PAIR_CARDINALITY",
                        error_detail=f"R1 files={len(r1s)}, R2 files={len(r2s)}",
                    )
                )
            else:
                validatable.append((sub, group, r1s[0], r2s[0]))

    results: list[PairResult] = list(structural_failures)

    if args.inventory_only:
        for sub, group, r1, r2 in validatable:
            results.append(
                PairResult(
                    sub=sub,
                    group=group,
                    status="NOT_RUN",
                    r1=str(r1.path),
                    r2=str(r2.path),
                    r1_size=r1.size,
                    r2_size=r2.size,
                    records=0,
                    elapsed_sec=0.0,
                    error_code="",
                    error_detail="Content validation skipped by --inventory-only",
                )
            )
    else:
        completed = 0
        total = len(validatable)
        print(f"Pairs eligible for full validation: {total}")
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(validate_pair, sub, group, r1, r2): (sub, group)
                for sub, group, r1, r2 in validatable
            }
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                completed += 1
                if result.status != "PASS":
                    print(
                        f"FAIL gRNA_sub{result.sub}/{result.group}: "
                        f"{result.error_code}: {result.error_detail}",
                        file=sys.stderr,
                        flush=True,
                    )
                    issues.append(
                        (
                            "ERROR",
                            result.error_code or "FASTQ_VALIDATION_FAILED",
                            f"gRNA_sub{result.sub}",
                            result.group,
                            "R1+R2",
                            result.error_detail,
                        )
                    )
                elif completed % args.progress_every == 0 or completed == total:
                    print(
                        f"Validated {completed}/{total} pairs; latest PASS: "
                        f"gRNA_sub{result.sub}/{result.group}, {result.records:,} reads",
                        flush=True,
                    )

    results.sort(key=lambda x: (x.sub, EXPECTED_GROUPS.index(x.group)))

    write_tsv(
        report_dir / "pair_results.tsv",
        (
            "sublibrary",
            "group",
            "status",
            "r1_path",
            "r2_path",
            "r1_bytes",
            "r2_bytes",
            "paired_records",
            "elapsed_seconds",
            "error_code",
            "error_detail",
        ),
        (
            (
                f"gRNA_sub{r.sub}",
                r.group,
                r.status,
                r.r1,
                r.r2,
                r.r1_size,
                r.r2_size,
                r.records,
                f"{r.elapsed_sec:.3f}",
                r.error_code,
                r.error_detail,
            )
            for r in results
        ),
    )

    write_tsv(
        report_dir / "issues.tsv",
        ("severity", "code", "sublibrary", "group", "read", "detail"),
        sorted(issues),
    )

    sample_rows = []
    for sub in EXPECTED_SUBS:
        sub_results = [r for r in results if r.sub == sub]
        passed = sum(r.status == "PASS" for r in sub_results)
        failed = sum(r.status == "FAIL" for r in sub_results)
        not_run = sum(r.status == "NOT_RUN" for r in sub_results)
        total_records = sum(r.records for r in sub_results if r.status == "PASS")
        sample_rows.append(
            (
                f"gRNA_sub{sub}",
                len(EXPECTED_GROUPS),
                passed,
                failed,
                not_run,
                total_records,
            )
        )

    write_tsv(
        report_dir / "sample_summary.tsv",
        (
            "sublibrary",
            "expected_group_pairs",
            "pass_pairs",
            "fail_pairs",
            "not_run_pairs",
            "total_paired_records_across_groups",
        ),
        sample_rows,
    )

    expected_pairs = len(EXPECTED_SUBS) * len(EXPECTED_GROUPS)
    pass_pairs = sum(r.status == "PASS" for r in results)
    fail_pairs = sum(r.status == "FAIL" for r in results)
    not_run_pairs = sum(r.status == "NOT_RUN" for r in results)
    error_count = sum(issue[0] == "ERROR" for issue in issues)
    warning_count = sum(issue[0] == "WARNING" for issue in issues)

    if args.inventory_only:
        overall = "INVENTORY_COMPLETE_WITH_ERRORS" if error_count else "INVENTORY_PASS"
    else:
        overall = (
            "PASS"
            if pass_pairs == expected_pairs and fail_pairs == 0 and error_count == 0
            else "FAIL"
        )

    summary_lines = [
        "gRNA FINAL FASTQ AUDIT",
        f"Generated: {time.strftime('%Y-%m-%dT%H:%M:%S%z')}",
        f"Root: {root}",
        f"Mode: {'inventory only' if args.inventory_only else 'full validation'}",
        f"Expected sublibraries: gRNA_sub1..gRNA_sub13",
        f"Expected groups per sublibrary: {len(EXPECTED_GROUPS)}",
        f"Expected R1/R2 pairs: {expected_pairs}",
        f"FASTQ.gz files discovered: {len(parsed)}",
        f"Pairs PASS: {pass_pairs}",
        f"Pairs FAIL: {fail_pairs}",
        f"Pairs NOT_RUN: {not_run_pairs}",
        f"Errors: {error_count}",
        f"Warnings: {warning_count}",
        f"Extra sublibraries detected: {', '.join('gRNA_sub'+str(x) for x in extra_subs) or 'none'}",
        f"OVERALL: {overall}",
        "",
        "Reports:",
        "  inventory.tsv",
        "  pair_results.tsv",
        "  sample_summary.tsv",
        "  issues.tsv",
    ]
    (report_dir / "summary.txt").write_text("\n".join(summary_lines) + "\n")

    marker = report_dir / ("AUDIT_PASS" if overall == "PASS" else "AUDIT_FAIL")
    marker.write_text(overall + "\n")

    print("\n" + "\n".join(summary_lines))
    return 0 if overall in {"PASS", "INVENTORY_PASS"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
