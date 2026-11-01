#!/usr/bin/env python3
"""Regenerate plotting tables from a public GEO accession or local submitted H5MU."""
from pathlib import Path
import argparse
import json
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cc_perturbseq.public_data import load_public_paths
from cc_perturbseq.geo_plot_data import prepare_geo_plot_inputs

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=["hippo_timecourse", "ipsc_k562_mix"], required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--geo-accession")
    source.add_argument("--local-file")
    parser.add_argument("--filename")
    parser.add_argument("--sha256")
    args = parser.parse_args()
    _, paths = load_public_paths(ROOT)
    report = prepare_geo_plot_inputs(args.dataset, ROOT, paths,
        mode="geo" if args.geo_accession else "local_file", accession=args.geo_accession or "",
        local_file=args.local_file, filename=args.filename, expected_sha256=args.sha256)
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
