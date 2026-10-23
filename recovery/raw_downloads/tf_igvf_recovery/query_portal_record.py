#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, os, sys
from pathlib import Path
from igvf_utils.connection import Connection

def norm_link(v):
    if isinstance(v, str): return v
    if isinstance(v, dict): return str(v.get('@id') or v.get('id') or '')
    return ''

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--row-tsv', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    a=p.parse_args()
    if not os.environ.get('IGVF_API_KEY') or not os.environ.get('IGVF_SECRET_KEY'):
        raise SystemExit('IGVF credentials are not set')
    with a.row_tsv.open(newline='') as f:
        row=next(csv.DictReader(f, delimiter='\t'))
    conn=Connection('prod', submission=True, no_log_file=True)
    obj={}; method=''
    for ident, name in ((row['aliases'],'alias'), (row['md5sum'],'md5sum')):
        try:
            obj=conn.get(ident, database=True, ignore404=True, frame='object') or {}
        except Exception as exc:
            a.output.parent.mkdir(parents=True, exist_ok=True)
            a.output.write_text(json.dumps({'status':'LOOKUP_ERROR','error':f'{type(exc).__name__}: {exc}'}, indent=2)+'\n')
            return 3
        if obj:
            method=name
            break
    result={
        'status': 'FOUND' if obj else 'NOT_FOUND',
        'lookup_method': method,
        'expected_alias': row['aliases'],
        'expected_md5sum': row['md5sum'],
        'expected_file_set': row['file_set'],
        'accession': str(obj.get('accession') or '') if obj else '',
        'portal_id': str(obj.get('@id') or '') if obj else '',
        'portal_file_set': norm_link(obj.get('file_set')) if obj else '',
        'portal_md5sum': str(obj.get('md5sum') or '') if obj else '',
        'portal_upload_status': str(obj.get('upload_status') or '') if obj else '',
        'portal_object_status': str(obj.get('status') or '') if obj else '',
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    return 0
if __name__=='__main__': raise SystemExit(main())
