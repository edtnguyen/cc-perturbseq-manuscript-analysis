#!/usr/bin/env python3
"""Download and prepare public figure inputs; no institutional files are read."""
from pathlib import Path
import argparse
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cc_perturbseq.public_data import configure_igvf, load_public_paths, prepare_figure_inputs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", action="append", choices=["Figure_01", "Figure_02", "Figure_03", "Figure_04", "Figure_S02", "Figure_S04", "Figure_S06"])
    parser.add_argument("--no-prompt", action="store_true", help="Use public/environment authentication only, never prompt interactively")
    args = parser.parse_args()
    configure_igvf(allow_prompt=not args.no_prompt)
    _, paths = load_public_paths(ROOT)
    failures = []
    for figure in args.figure or ["Figure_01", "Figure_02", "Figure_03", "Figure_04", "Figure_S02", "Figure_S04", "Figure_S06"]:
        try:
            print(prepare_figure_inputs(figure, ROOT, paths))
        except Exception as exc:
            failures.append(figure)
            print(f"{figure}: {type(exc).__name__}: {exc}", file=sys.stderr)
    raise SystemExit(1 if failures else 0)
