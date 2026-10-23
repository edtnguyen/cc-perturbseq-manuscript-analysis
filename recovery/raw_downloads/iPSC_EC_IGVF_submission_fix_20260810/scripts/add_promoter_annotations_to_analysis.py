#!/usr/bin/env python3
"""Add Ensembl ID and promoter coordinates to a MAGeCK gene-summary TSV.

Exact promoter-target IDs (for example NOTCH1_-P1) are matched first.  Gene-only
IDs are matched to the gene-level lookup; genes with multiple promoter targets
retain promoter-labelled semicolon-separated coordinates.
"""
import argparse, csv, re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET_LOOKUP = ROOT / "annotations" / "IGVFFI6745HXWN_analysis_target_annotations.tsv"
GENE_LOOKUP = ROOT / "annotations" / "D3pos_D3neg_gene_summary_annotation_lookup.tsv"
ANN = ["ensembl_gene_id", "promoter_coordinate", "intended_target_chr", "intended_target_start", "intended_target_end"]

def read_lookup(path, key):
    with path.open(newline="", encoding="utf-8") as f:
        return {r[key]: r for r in csv.DictReader(f, delimiter="\t")}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("output", nargs="?")
    a = p.parse_args()
    inp = Path(a.input)
    out = Path(a.output) if a.output else inp.with_name(inp.stem + "_with_promoter_annotations.txt")
    target = read_lookup(TARGET_LOOKUP, "analysis_id")
    gene = read_lookup(GENE_LOOKUP, "gene_symbol")
    with inp.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames or "id" not in r.fieldnames:
            raise SystemExit("Input file has no 'id' column")
        rows = list(r); original = list(r.fieldnames)
    missing = []
    for row in rows:
        raw = row["id"].strip()
        ann = target.get(raw)
        if ann is None:
            symbol = re.sub(r"_-.*$", "", raw)
            ann = gene.get(symbol)
        if ann is None:
            missing.append(raw); continue
        for c in ANN: row[c] = ann.get(c, "")
    if missing:
        raise SystemExit(f"{len(missing)} IDs did not map; examples: {missing[:20]}")
    fields = ["id"] + ANN + [c for c in original if c != "id"]
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"Wrote {out} ({len(rows)} rows)")

if __name__ == "__main__": main()
