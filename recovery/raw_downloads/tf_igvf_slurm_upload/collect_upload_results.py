#!/usr/bin/env python3
"""Collect per-task upload markers/logs into a CSV and summary."""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work-dir', required=True, type=Path)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    work = args.work_dir.resolve()
    manifest_path = work / 'upload_manifest.csv'
    with manifest_path.open(newline='') as f:
        manifest = list(csv.DictReader(f))

    output = []
    failed = []
    for row in manifest:
        idx = int(row['task_index'])
        marker_path = work / 'upload_done' / f'{idx:06d}.json'
        log_path = work / 'upload_logs' / f'{idx:06d}.log'
        status = 'MISSING'
        accession = ''
        error = ''
        if marker_path.is_file():
            try:
                marker = json.loads(marker_path.read_text())
                status = marker.get('status', 'UNKNOWN')
                accession = marker.get('accession', '')
            except Exception as exc:
                status = 'ERROR'
                error = f'Unreadable marker: {exc}'
        elif log_path.is_file():
            text = log_path.read_text(errors='replace')
            m = re.findall(r'IGVFFI[A-Z0-9]+', text)
            accession = m[-1] if m else ''
            status = 'FAILED'
            tail = [line for line in text.splitlines() if re.search(
                r'Traceback|ValidationFailure|Failed|ERROR|Forbidden|Unprocessable', line, re.I
            )]
            error = ' | '.join(tail[-3:])[:2000]

        out = dict(row)
        out.update({'status': status, 'accession': accession, 'error': error, 'log': str(log_path)})
        output.append(out)
        if status != 'OK' or not accession:
            failed.append(idx)

    fields = list(manifest[0].keys()) + ['status', 'accession', 'error', 'log'] if manifest else []
    out_csv = work / 'sequence_file_upload_results.csv'
    with out_csv.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(output)

    summary = [
        'TF sequence_file Slurm upload summary',
        f'Expected tasks: {len(manifest)}',
        f'Successful tasks: {len(manifest) - len(failed)}',
        f'Failed or missing tasks: {len(failed)}',
        f'Results CSV: {out_csv}',
    ]
    if failed:
        summary.append('Failed task indices: ' + ','.join(map(str, failed)))
    (work / 'upload_summary.txt').write_text('\n'.join(summary) + '\n')
    print('\n'.join(summary))
    return 0 if not failed else 2


if __name__ == '__main__':
    raise SystemExit(main())
