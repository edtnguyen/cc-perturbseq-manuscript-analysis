#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'src'))
from control_defined_state_atlas import run_review, run_finalize

p=argparse.ArgumentParser(description='Build a control/reference-defined coarse developmental-state atlas.')
p.add_argument('--config',required=True)
p.add_argument('--mode',choices=['review','finalize','all'],default='all')
a=p.parse_args()

if a.mode in {'review','all'}:
    out=run_review(a.config)
    print(f'Review outputs: {out}')
if a.mode in {'finalize','all'}:
    out=run_finalize(a.config)
    print(f'Final outputs: {out}')
