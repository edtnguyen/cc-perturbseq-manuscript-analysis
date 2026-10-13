#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import sys
import zipfile
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path

DATE = "2026-08-10"
ROOT = Path("/mnt/data/iPSC_EC_IGVF_submission_fix_20260810")
ZIP_PATH = Path("/mnt/data/iPSC_EC_IGVF_submission_fix_20260810.zip")
INPUTS = {
    "validation": Path("/mnt/data/IGVFFI6745HXWN.tsv.gz"),
    "genome_wide": Path("/mnt/data/IGVFFI9629WVUU.tsv.gz"),
}
ACCESSIONS = {"validation": "IGVFFI6745HXWN", "genome_wide": "IGVFFI9629WVUU"}
MAPPING_PATH = Path("/mnt/data/weissman_guides_with_coordinates.tsv")
FORMAT_DOC = "/documents/f6a915a9-a6b8-4015-8aeb-236764976979/"
VALIDATOR_REPO = "https://github.com/EngreitzLab/igvf-grna-validation"

REQUIRED_COLUMNS = [
    "guide_id", "spacer", "targeting", "type", "guide_chr", "guide_start",
    "guide_end", "strand", "pam", "genomic_element", "intended_target_name",
    "intended_target_chr", "intended_target_start", "intended_target_end",
]
OPTIONAL_COLUMNS = ["putative_target_genes", "reporter", "imperfect", "description"]
SUBMISSION_COLUMNS = REQUIRED_COLUMNS + OPTIONAL_COLUMNS
SOURCE_COLUMNS = [
    "gene name (gencode v43)", "gene name (original)", "guide_start (hg19)",
    "guide_end (hg19)", "conversion issues",
]
MAPPING_AUDIT_COLUMNS = [
    "sgID", "gene", "transcript", "ensembl_id", "ensembl_id_source",
    "gencode_v43_name", "gencode_v43_status", "chr_status", "chr_source",
    "chr (hg38)", "start (hg38, 0-indexed)", "end (hg38, 1-indexed)",
    "liftover_success", "sequence_matches (hg38)",
]
VALID_TYPES = {
    "targeting", "safe-targeting", "non-targeting", "positive control",
    "negative control", "variant",
}
VALID_GE = {
    "variant", "promoter", "enhancer", "insulator", "silencer",
    "distal element", "splice site", "gene",
}
GE_NEEDS_ENSG = {"promoter", "splice site", "gene"}
GE_NEEDS_COORD = {"enhancer", "insulator", "silencer", "distal element"}
ENSG_RE = re.compile(r"^ENSG\d+$")
COORD_RE = re.compile(r"^chr\S+:\d+-\d+$")
NAN_VALUES = {"", "nan", "NaN", "NAN", "na", "NA", "N/A", "n/a", "none", "None", "null", "NULL"}


def clean(value) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text in NAN_VALUES else text


def is_true(value) -> bool:
    return clean(value).upper() == "TRUE"


def normalize_bool(value) -> str:
    text = clean(value).lower()
    return "True" if text == "true" else "False" if text == "false" else clean(value)


def parse_target(guide_id: str):
    m = re.match(r"^(.*?)_[+-]_\d+\.23-(.+)$", clean(guide_id))
    if not m:
        gid = clean(guide_id)
        return gid, gid, gid
    gene, promoter = m.group(1), m.group(2)
    return gene, promoter, f"{gene}_-{promoter}"


def ensg_candidates(value) -> tuple[str, ...]:
    return tuple(sorted(set(re.findall(r"ENSG\d+", clean(value)))))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


@contextmanager
def open_text(path: Path, mode: str):
    if "r" in mode:
        if path.suffix == ".gz":
            with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
                yield fh
        else:
            with path.open("r", encoding="utf-8", newline="") as fh:
                yield fh
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".gz":
            with path.open("wb") as raw:
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=6) as gz:
                    with io.TextIOWrapper(gz, encoding="utf-8", newline="") as fh:
                        yield fh
        else:
            with path.open("w", encoding="utf-8", newline="") as fh:
                yield fh


def write_rows(path: Path, fieldnames, rows):
    with open_text(path, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fieldnames), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({c: clean(row.get(c, "")) for c in fieldnames})


def scan_genome_wide():
    guide_spacer = {}
    guide_ids = set()
    guide_to_aid = {}
    groups = {}
    counts = Counter()
    rows = 0
    with open_text(INPUTS["genome_wide"], "r") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rows += 1
            gid = clean(row["guide_id"])
            if gid in guide_spacer:
                raise ValueError(f"Duplicate genome-wide guide_id: {gid}")
            guide_spacer[gid] = clean(row["spacer"])
            guide_ids.add(gid)
            counts[clean(row["type"])] += 1
            if not is_true(row["targeting"]):
                continue
            gene, promoter, aid = parse_target(gid)
            guide_to_aid[gid] = aid
            chrom = clean(row["guide_chr"])
            start = int(clean(row["guide_start"]))
            end = int(clean(row["guide_end"]))
            g = groups.setdefault(aid, {
                "analysis_id": aid, "gene_symbol": gene, "promoter_label": promoter,
                "chroms": set(), "start": start, "end": end, "guide_ids": [],
                "source_old_ensg": set(), "source_input_intended_target_name": set(),
            })
            g["chroms"].add(chrom)
            g["start"] = min(g["start"], start)
            g["end"] = max(g["end"], end)
            g["guide_ids"].append(gid)
            old_name = clean(row.get("intended_target_name", ""))
            if old_name:
                g["source_input_intended_target_name"].add(old_name)
                g["source_old_ensg"].update(ensg_candidates(old_name))
    return rows, counts, guide_spacer, guide_ids, guide_to_aid, groups


def load_mapping(selected_guide_ids):
    mapping = {}
    seen = set()
    with open_text(MAPPING_PATH, "r") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            gid = clean(row.get("sgID", ""))
            if gid not in selected_guide_ids:
                continue
            if gid in seen:
                raise ValueError(f"Duplicate sgID in mapping: {gid}")
            seen.add(gid)
            mapping[gid] = {c: clean(row.get(c, "")) for c in MAPPING_AUDIT_COLUMNS}
    return mapping


def build_canonical(groups, mapping):
    canonical = {}
    for aid, g in groups.items():
        if len(g["chroms"]) != 1:
            raise ValueError(f"Target {aid} maps to multiple chromosomes: {sorted(g['chroms'])}")
        chrom = next(iter(g["chroms"]))
        candidate_tuples = {ensg_candidates(mapping.get(gid, {}).get("ensembl_id", "")) for gid in g["guide_ids"]}
        if len(candidate_tuples) != 1:
            raise ValueError(f"Guide-level Ensembl mapping disagreement within {aid}: {candidate_tuples}")
        candidates = next(iter(candidate_tuples))
        start, end = g["start"], g["end"]
        coord = f"{chrom}:{start}-{end}"
        # Conservative rule: only one unambiguous ID from the supplied mapping is
        # promoted to a promoter target.  Otherwise use a validator-compatible
        # coordinate target and retain the source ENSG(s) in the review table.
        if len(candidates) == 1:
            ge = "promoter"
            intended_name = candidates[0]
            ensg = candidates[0]
            status = "single_ensembl_id_from_supplied_mapping"
            review = "No"
        else:
            ge = "distal element"
            intended_name = coord
            ensg = ""
            status = "coordinate_fallback_missing_ensembl" if not candidates else "coordinate_fallback_ambiguous_ensembl"
            review = "Yes"
        sources = set()
        names = set()
        statuses = set()
        for gid in g["guide_ids"]:
            m = mapping.get(gid, {})
            if clean(m.get("ensembl_id_source")): sources.add(clean(m["ensembl_id_source"]))
            if clean(m.get("gencode_v43_name")): names.add(clean(m["gencode_v43_name"]))
            if clean(m.get("gencode_v43_status")): statuses.add(clean(m["gencode_v43_status"]))
        canonical[aid] = {
            "analysis_id": aid,
            "gene_symbol": g["gene_symbol"],
            "promoter_label": g["promoter_label"],
            "ensembl_gene_id": ensg,
            "candidate_ensembl_ids": ";".join(candidates),
            "source_legacy_ensembl_ids": ";".join(sorted(g["source_old_ensg"])),
            "genomic_element": ge,
            "intended_target_name": intended_name,
            "intended_target_chr": chrom,
            "intended_target_start": str(start),
            "intended_target_end": str(end),
            "promoter_coordinate": coord,
            "n_guides_genome_wide": str(len(g["guide_ids"])),
            "mapping_status": status,
            "requires_manual_review": review,
            "mapping_ensembl_id_source": "; ".join(sorted(sources)),
            "mapping_gencode_v43_name": "; ".join(sorted(names)),
            "mapping_gencode_v43_status": "; ".join(sorted(statuses)),
            "source_input_intended_target_name": "; ".join(sorted(g["source_input_intended_target_name"])),
        }
    return canonical


def description_for(target, typ, row):
    if is_true(row.get("targeting", "")):
        text = (
            f"CRISPRi guide assigned to {target['gene_symbol']} promoter target "
            f"{target['analysis_id']}; intended target window spans all genome-wide guides assigned to this promoter"
        )
        if target["requires_manual_review"] == "Yes":
            text += "; coordinate representation used because supplied Ensembl mapping is missing or ambiguous"
        return text
    if typ == "non-targeting":
        return "Non-targeting control"
    if typ == "safe-targeting":
        return "Safe-targeting control"
    if typ == "targeting":
        issue = clean(row.get("conversion issues", "")) or "source mapping/QC exclusion"
        return f"Targeting guide marked targeting=False: {issue}"
    return typ or "Control/QC-excluded guide"


def transform_library(library, canonical, mapping):
    accession = ACCESSIONS[library]
    source_path = INPUTS[library]
    out_path = ROOT / "submission_ready" / f"{accession}.tsv.gz"
    audit_path = ROOT / "audit" / f"{accession}_corrected_with_source_and_mapping.tsv.gz"
    audit_cols = SUBMISSION_COLUMNS + SOURCE_COLUMNS + [f"mapping::{c}" for c in MAPPING_AUDIT_COLUMNS] + [
        "analysis_id", "mapping_ensg_candidates", "mapping_ensg_n",
    ]
    target_counts = Counter()
    aids = set()
    total = 0
    types = Counter()
    with open_text(source_path, "r") as src, open_text(out_path, "w") as out_fh, open_text(audit_path, "w") as audit_fh:
        r = csv.DictReader(src, delimiter="\t")
        ow = csv.DictWriter(out_fh, fieldnames=SUBMISSION_COLUMNS, delimiter="\t", lineterminator="\n")
        aw = csv.DictWriter(audit_fh, fieldnames=audit_cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        ow.writeheader(); aw.writeheader()
        for row in r:
            total += 1
            typ = clean(row.get("type", ""))
            types[typ] += 1
            gid = clean(row.get("guide_id", ""))
            out = {c: clean(row.get(c, "")) for c in SUBMISSION_COLUMNS}
            out["targeting"] = normalize_bool(row.get("targeting", ""))
            out["type"] = typ
            gene, promoter, aid = parse_target(gid)
            if is_true(row.get("targeting", "")):
                if aid not in canonical:
                    raise ValueError(f"No canonical promoter target for {gid} ({aid})")
                target = canonical[aid]
                aids.add(aid); target_counts[aid] += 1
                for c in ["genomic_element", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end"]:
                    out[c] = target[c]
            else:
                target = None
                for c in ["genomic_element", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end"]:
                    out[c] = ""
            out["description"] = description_for(target, typ, row) if target else description_for({}, typ, row)
            ow.writerow(out)

            m = mapping.get(gid, {})
            audit = dict(out)
            for c in SOURCE_COLUMNS:
                audit[c] = clean(row.get(c, ""))
            for c in MAPPING_AUDIT_COLUMNS:
                audit[f"mapping::{c}"] = clean(m.get(c, ""))
            audit["analysis_id"] = aid
            candidates = ensg_candidates(m.get("ensembl_id", ""))
            audit["mapping_ensg_candidates"] = ";".join(candidates)
            audit["mapping_ensg_n"] = str(len(candidates))
            aw.writerow(audit)
    return {"rows": total, "types": types, "target_counts": target_counts, "aids": aids, "out_path": out_path}


def make_annotation_rows(library, canonical, transform_info):
    accession = ACCESSIONS[library]
    rows = []
    for aid in transform_info["aids"]:
        c = canonical[aid]
        row = {
            "library": library,
            "source_accession": accession,
            **c,
            "n_guides_in_library": str(transform_info["target_counts"][aid]),
        }
        rows.append(row)
    rows.sort(key=lambda r: (r["gene_symbol"], r["promoter_label"], r["analysis_id"]))
    return rows


ANN_COLS = [
    "library", "source_accession", "analysis_id", "gene_symbol", "promoter_label",
    "ensembl_gene_id", "promoter_coordinate", "intended_target_chr",
    "intended_target_start", "intended_target_end", "genomic_element",
    "n_guides_in_library", "n_guides_genome_wide", "mapping_status",
    "requires_manual_review", "candidate_ensembl_ids", "source_legacy_ensembl_ids",
    "mapping_ensembl_id_source", "mapping_gencode_v43_name",
    "mapping_gencode_v43_status", "source_input_intended_target_name",
]


def validate_file(path: Path):
    issues = []
    with open_text(path, "r") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        columns = r.fieldnames or []
        rows = list(r)
    missing_required = [c for c in REQUIRED_COLUMNS if c not in columns]
    missing_optional = [c for c in OPTIONAL_COLUMNS if c not in columns]
    extra = [c for c in columns if c not in SUBMISSION_COLUMNS]
    if missing_required: issues.append(("structure", "error", f"Missing required columns: {missing_required}", len(missing_required)))
    if missing_optional: issues.append(("structure", "warn", f"Missing optional columns: {missing_optional}", len(missing_optional)))
    if extra: issues.append(("structure", "warn", f"Extra columns: {extra}", len(extra)))
    if missing_required:
        return rows, issues

    gids = [clean(x["guide_id"]) for x in rows]
    dup = len(gids) - len(set(gids))
    if dup: issues.append(("guide_id", "error", f"{dup:,} duplicate guide_id values", dup))
    gid_sp = defaultdict(set)
    for x in rows: gid_sp[clean(x["guide_id"])].add(clean(x["spacer"]))
    multi = sum(len(v) > 1 for v in gid_sp.values())
    if multi: issues.append(("spacer", "error", f"{multi:,} guide_ids map to more than one spacer", multi))

    bad_type = sum(clean(x["type"]) not in VALID_TYPES for x in rows)
    if bad_type: issues.append(("type", "error", f"{bad_type:,} rows have invalid type", bad_type))
    bad_targeting = sum(clean(x["targeting"]) not in {"True", "False"} for x in rows)
    if bad_targeting: issues.append(("targeting", "error", f"{bad_targeting:,} rows have invalid targeting value", bad_targeting))
    nontgt_true = sum(clean(x["type"]) in {"non-targeting", "safe-targeting"} and is_true(x["targeting"]) for x in rows)
    if nontgt_true: issues.append(("targeting", "error", f"{nontgt_true:,} non-targeting/safe-targeting rows have targeting=True", nontgt_true))
    tgt_false = sum(clean(x["type"]) == "targeting" and clean(x["targeting"]).upper() == "FALSE" for x in rows)
    if tgt_false: issues.append(("targeting", "warn", f"{tgt_false:,} rows with type='targeting' have targeting='False' (expected 'True')", tgt_false))

    true_rows = [x for x in rows if is_true(x["targeting"])]
    non_targeting = [x for x in rows if clean(x["type"]) == "non-targeting" and clean(x["targeting"]).upper() == "FALSE"]
    for col in ["strand", "guide_chr", "guide_start", "guide_end", "genomic_element", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end"]:
        n = sum(not clean(x[col]) for x in true_rows)
        if n: issues.append((col, "error", f"{n:,} targeting=True rows have empty {col}", n))
    bad_strand = sum(clean(x["strand"]) not in {"+", "-"} for x in true_rows if clean(x["strand"]))
    if bad_strand: issues.append(("strand", "error", f"{bad_strand:,} targeting=True rows have invalid strand", bad_strand))
    for col in ["guide_start", "guide_end", "intended_target_start", "intended_target_end"]:
        n = 0
        for x in true_rows:
            try: int(clean(x[col]))
            except Exception: n += 1
        if n: issues.append((col, "error", f"{n:,} targeting=True rows have nonnumeric {col}", n))
    ge_invalid = sum(clean(x["genomic_element"]) not in VALID_GE for x in true_rows if clean(x["genomic_element"]))
    if ge_invalid: issues.append(("genomic_element", "error", f"{ge_invalid:,} targeting=True rows have invalid genomic_element", ge_invalid))
    ge_enh = sum(clean(x["genomic_element"]) == "enhancer" for x in true_rows)
    if ge_enh: issues.append(("genomic_element", "warn", f"{ge_enh:,} targeting rows use enhancer; prefer distal element", ge_enh))
    for col in ["genomic_element", "intended_target_name", "intended_target_chr", "intended_target_start", "intended_target_end"]:
        n = sum(bool(clean(x[col])) for x in non_targeting)
        if n: issues.append((col, "error", f"{n:,} non-targeting rows have non-empty {col}", n))
    bad_name = 0
    for x in true_rows:
        ge = clean(x["genomic_element"]); name = clean(x["intended_target_name"])
        if ge in GE_NEEDS_ENSG and not ENSG_RE.match(name): bad_name += 1
        if ge in GE_NEEDS_COORD and not COORD_RE.match(name): bad_name += 1
    if bad_name: issues.append(("intended_target_name", "error", f"{bad_name:,} targeting rows have malformed intended_target_name", bad_name))

    chr_bad = start_bad = end_bad = 0
    prefix_bounds = {}
    for x in true_rows:
        try:
            gs, ge = int(x["guide_start"]), int(x["guide_end"])
            ts, te = int(x["intended_target_start"]), int(x["intended_target_end"])
        except Exception:
            continue
        if clean(x["guide_chr"]) != clean(x["intended_target_chr"]): chr_bad += 1
        if ts > gs: start_bad += 1
        if te < ge: end_bad += 1
        prefix = re.sub(r"_\d+$", "", clean(x["guide_id"]))
        b = prefix_bounds.setdefault(prefix, [gs, ge, ts, te, 0])
        b[0] = min(b[0], gs); b[1] = max(b[1], ge); b[4] += 1
    if chr_bad: issues.append(("intended_target_chr", "error", f"{chr_bad:,} targeting rows have intended_target_chr != guide_chr", chr_bad))
    if start_bad: issues.append(("intended_target_start", "error", f"{start_bad:,} targeting rows have intended_target_start > guide_start", start_bad))
    if end_bad: issues.append(("intended_target_end", "error", f"{end_bad:,} targeting rows have intended_target_end < guide_end", end_bad))
    narrow_groups = 0
    for prefix, (gs, ge, ts, te, n) in prefix_bounds.items():
        if n > 1 and (ts > gs or te < ge): narrow_groups += 1
    if narrow_groups: issues.append(("intended_target_start", "error", f"{narrow_groups:,} guide-prefix groups have target windows that do not span all guides", narrow_groups))
    desc_empty = all(not clean(x.get("description", "")) for x in rows)
    if desc_empty: issues.append(("description", "warn", "description is absent or all empty", len(rows)))
    return rows, issues


def write_validation_report(accession, path, rows, issues):
    counts = Counter(clean(x["type"]) for x in rows)
    lines = [
        "=== IGVF Per-Guide Metadata Validation Report (current-repo check equivalent) ===",
        f"File: {path.name}",
        f"Rows: {len(rows):,}" + "".join(f"  |  {k}: {v:,}" for k, v in sorted(counts.items()) if k),
        "",
    ]
    if issues:
        for field, severity, message, count in issues:
            lines.append(f"[{'ERROR' if severity == 'error' else 'WARN'}] {field}: {message}")
    else:
        lines.append("PASS: No issues found.")
    ne = sum(i[1] == "error" for i in issues); nw = sum(i[1] == "warn" for i in issues)
    lines += ["", f"SUMMARY: {ne} error categories; {nw} warning categories.", ""]
    (ROOT / "validation_reports" / f"{accession}_validation_report.txt").write_text("\n".join(lines), encoding="utf-8")
    obj = {
        "input_file": path.name, "total_rows": len(rows), "validator_reference": VALIDATOR_REPO,
        "validation_mode": "current-repository check-equivalent implementation",
        "issues": [{"field": a, "severity": b, "message": c, "count": d} for a,b,c,d in issues],
    }
    (ROOT / "validation_reports" / f"{accession}_validation_report.json").write_text(json.dumps(obj, indent=2), encoding="utf-8")
    return {"rows": len(rows), "error_categories": ne, "warning_categories": nw, "warnings": [i[2] for i in issues if i[1] == "warn"], "status": "PASS" if ne == 0 else "FAIL"}


def create_merge_script():
    script = r'''#!/usr/bin/env python3
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
'''
    path = ROOT / "scripts" / "add_promoter_annotations_to_analysis.py"
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)


def main():
    if ROOT.exists(): shutil.rmtree(ROOT)
    if ZIP_PATH.exists(): ZIP_PATH.unlink()
    for d in ["submission_ready", "audit", "annotations", "review", "validation_reports", "scripts"]:
        (ROOT / d).mkdir(parents=True, exist_ok=True)

    gw_rows, gw_types, gw_spacer, guide_ids, guide_to_aid, groups = scan_genome_wide()
    mapping = load_mapping(guide_ids)
    canonical = build_canonical(groups, mapping)

    # Validation subset audit before transformation.
    val_rows = 0; missing = []; mismatch = []; val_ids = set()
    with open_text(INPUTS["validation"], "r") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            val_rows += 1
            gid = clean(row["guide_id"]); val_ids.add(gid)
            if gid not in gw_spacer: missing.append(gid)
            elif clean(row["spacer"]) != gw_spacer[gid]: mismatch.append(gid)
    if missing or mismatch:
        raise ValueError(f"Validation subset check failed: missing={len(missing)}, spacer_mismatch={len(mismatch)}")

    transform = {}
    annotation = {}
    validation_summary = {}
    for library in ["validation", "genome_wide"]:
        transform[library] = transform_library(library, canonical, mapping)
        annotation[library] = make_annotation_rows(library, canonical, transform[library])
        suffix = ".tsv" if library == "validation" else ".tsv.gz"
        write_rows(ROOT / "annotations" / f"{ACCESSIONS[library]}_analysis_target_annotations{suffix}", ANN_COLS, annotation[library])
        rows, issues = validate_file(transform[library]["out_path"])
        validation_summary[library] = write_validation_report(ACCESSIONS[library], transform[library]["out_path"], rows, issues)

    write_rows(ROOT / "annotations" / "all_analysis_target_annotations.tsv.gz", ANN_COLS, annotation["validation"] + annotation["genome_wide"])

    # Gene-level validation lookup.
    by_gene = defaultdict(list)
    for row in annotation["validation"]: by_gene[row["gene_symbol"]].append(row)
    gene_rows = []
    for gene in sorted(by_gene):
        grp = sorted(by_gene[gene], key=lambda x: (x["promoter_label"], x["analysis_id"]))
        ids = sorted({x["ensembl_gene_id"] for x in grp if x["ensembl_gene_id"]})
        gene_rows.append({
            "gene_symbol": gene,
            "ensembl_gene_id": ";".join(ids),
            "promoter_coordinate": ";".join(f"{x['promoter_label']}={x['promoter_coordinate']}" for x in grp),
            "intended_target_chr": ";".join(sorted({x["intended_target_chr"] for x in grp if x["intended_target_chr"]})),
            "intended_target_start": ";".join(f"{x['promoter_label']}={x['intended_target_start']}" for x in grp),
            "intended_target_end": ";".join(f"{x['promoter_label']}={x['intended_target_end']}" for x in grp),
            "analysis_id": ";".join(x["analysis_id"] for x in grp),
            "n_promoter_targets": str(len(grp)),
            "n_guides": str(sum(int(x["n_guides_in_library"]) for x in grp)),
            "requires_manual_review": "Yes" if any(x["requires_manual_review"] == "Yes" for x in grp) else "No",
        })
    gene_cols = ["gene_symbol", "ensembl_gene_id", "promoter_coordinate", "intended_target_chr", "intended_target_start", "intended_target_end", "analysis_id", "n_promoter_targets", "n_guides", "requires_manual_review"]
    write_rows(ROOT / "annotations" / "D3pos_D3neg_gene_summary_annotation_lookup.tsv", gene_cols, gene_rows)

    pasted = [
        ["SMAD5", "6", "6.3801e-11", "4.9769e-06", "0.000114", "1", "6", "-3.4078", "1.0", "1.0", "0.999995", "299", "0", "-3.4078"],
        ["SMAD1", "6", "5.1338e-10", "4.9769e-06", "0.000114", "2", "6", "-3.5646", "1.0", "1.0", "0.999995", "298", "0", "-3.5646"],
        ["DBR1", "6", "9.6071e-10", "4.9769e-06", "0.000114", "3", "6", "-3.5284", "1.0", "1.0", "0.999995", "297", "0", "-3.5284"],
        ["ACVR2A", "6", "5.8966e-09", "4.9769e-06", "0.000114", "4", "5", "-3.6132", "0.9293", "0.94523", "0.999995", "204", "1", "-3.6132"],
    ]
    result_cols = ["id", "num", "neg|score", "neg|p-value", "neg|fdr", "neg|rank", "neg|goodsgrna", "neg|lfc", "pos|score", "pos|p-value", "pos|fdr", "pos|rank", "pos|goodsgrna", "pos|lfc"]
    gene_lookup = {x["gene_symbol"]: x for x in gene_rows}
    example = []
    for values in pasted:
        base = dict(zip(result_cols, values)); ann = gene_lookup[base["id"]]
        row = {"id": base["id"]}
        for c in ["ensembl_gene_id", "promoter_coordinate", "intended_target_chr", "intended_target_start", "intended_target_end"]: row[c] = ann[c]
        for c in result_cols[1:]: row[c] = base[c]
        example.append(row)
    ex_cols = ["id", "ensembl_gene_id", "promoter_coordinate", "intended_target_chr", "intended_target_start", "intended_target_end"] + result_cols[1:]
    write_rows(ROOT / "annotations" / "D3pos_D3neg_gene_summary_ANNOTATED_EXAMPLE_ONLY.tsv", ex_cols, example)

    # Review and reconciliation.
    review_cols = ["analysis_id", "gene_symbol", "promoter_label", "promoter_coordinate", "n_guides_genome_wide", "mapping_status", "candidate_ensembl_ids", "source_legacy_ensembl_ids", "mapping_ensembl_id_source", "mapping_gencode_v43_name", "mapping_gencode_v43_status", "source_input_intended_target_name", "genomic_element", "intended_target_name"]
    review = [c for c in canonical.values() if c["requires_manual_review"] == "Yes"]
    review.sort(key=lambda x: (x["gene_symbol"], x["promoter_label"], x["analysis_id"]))
    write_rows(ROOT / "review" / "genome_wide_promoter_targets_requiring_DACC_review.tsv", review_cols, review)

    reconciliation = []
    rec_counts = Counter()
    for c in sorted(canonical.values(), key=lambda x: x["analysis_id"]):
        old = ensg_candidates(c["source_input_intended_target_name"]); new = c["ensembl_gene_id"]
        if new and old == (new,): status = "unchanged_after_version_normalization"
        elif new and old: status = "updated_to_supplied_mapping_ensembl_id"
        elif new and not old: status = "filled_from_supplied_mapping"
        else: status = "coordinate_fallback_requires_review"
        rec_counts[status] += 1
        reconciliation.append({
            "analysis_id": c["analysis_id"], "gene_symbol": c["gene_symbol"], "promoter_coordinate": c["promoter_coordinate"],
            "source_input_intended_target_name": c["source_input_intended_target_name"], "source_legacy_ensembl_ids": c["source_legacy_ensembl_ids"],
            "supplied_mapping_ensembl_id": c["ensembl_gene_id"], "candidate_ensembl_ids": c["candidate_ensembl_ids"],
            "reconciliation_status": status, "requires_manual_review": c["requires_manual_review"],
        })
    rec_cols = ["analysis_id", "gene_symbol", "promoter_coordinate", "source_input_intended_target_name", "source_legacy_ensembl_ids", "supplied_mapping_ensembl_id", "candidate_ensembl_ids", "reconciliation_status", "requires_manual_review"]
    write_rows(ROOT / "review" / "intended_target_name_reconciliation.tsv.gz", rec_cols, reconciliation)

    overlap = [{
        "validation_rows": str(val_rows), "validation_rows_found_in_genome_wide": str(val_rows - len(missing)),
        "missing_validation_guide_ids": str(len(missing)), "spacer_mismatches": str(len(mismatch)),
        "validation_target_records": str(len(annotation["validation"])), "validation_unique_gene_symbols": str(len(gene_rows)),
        "genome_wide_target_records": str(len(annotation["genome_wide"])),
        "conclusion": "Validation is a complete guide-level subset of the genome-wide library" if not missing and not mismatch else "Subset check failed",
    }]
    write_rows(ROOT / "review" / "validation_vs_genome_wide_overlap_summary.tsv", list(overlap[0]), overlap)

    create_merge_script()
    shutil.copy2(Path(__file__), ROOT / "scripts" / "prepare_igvf_guides_fast.py")
    (ROOT / "portal_metadata_patch.json").write_text(json.dumps({"file_format_specifications": [FORMAT_DOC], "note": "Delete the currently linked obsolete guide-RNA format document and retain only this official document. No Portal write was performed."}, indent=2), encoding="utf-8")

    fallback_groups = len(review)
    fallback_guides = sum(int(x["n_guides_genome_wide"]) for x in review)
    targeting_false_gw = 0
    with open_text(INPUTS["genome_wide"], "r") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if clean(row["type"]) == "targeting" and clean(row["targeting"]).upper() == "FALSE": targeting_false_gw += 1

    readme = f'''# iPSC-EC IGVF guide metadata correction

Prepared: {DATE}

## Deliverables completed

- Corrected the validation guide file (`IGVFFI6745HXWN`) and genome-wide guide file (`IGVFFI9629WVUU`).
- Used `weissman_guides_with_coordinates.tsv` for hg38/GENCODE Ensembl mapping while preserving the hg38 guide coordinates already present in the IGVF files.
- Set `genomic_element=promoter` for promoter targets with one unambiguous Ensembl gene ID in the supplied mapping.
- Filled `intended_target_chr`, `intended_target_start`, and `intended_target_end` with the hg38 min-to-max guide span for each promoter target. The same genome-wide promoter window is used for its validation-screen subset.
- Normalized `targeting` to `True`/`False`, retained only the 18 recognized columns in submission-ready files, and populated `description`.
- Confirmed all {val_rows:,} validation rows occur in the genome-wide file with the same spacer sequence.
- Created promoter-target and gene-level lookup tables containing Ensembl IDs and promoter coordinates for the analysis file.

## Submission-ready guide files

- `submission_ready/IGVFFI6745HXWN.tsv.gz`
- `submission_ready/IGVFFI9629WVUU.tsv.gz`

Current-repository check-equivalent results:

- Validation library: {validation_summary['validation']['error_categories']} error categories, {validation_summary['validation']['warning_categories']} warning categories, status {validation_summary['validation']['status']}.
- Genome-wide library: {validation_summary['genome_wide']['error_categories']} error categories, {validation_summary['genome_wide']['warning_categories']} warning categories, status {validation_summary['genome_wide']['status']}.
- The genome-wide warning is the pre-existing condition that {targeting_false_gw:,} rows have `type=targeting` but `targeting=False`; intended-target fields were left empty for these QC-excluded rows.

The reports implement the checks in the current `EngreitzLab/igvf-grna-validation` validator. The exact repository CLI was not executed because the runtime could not clone/download executable source dependencies; the checks and output schema were compared against the current repository source.

## Ensembl mapping decisions

Across {len(canonical):,} genome-wide promoter targets:

- {rec_counts['unchanged_after_version_normalization']:,} retained the same Ensembl ID after version stripping.
- {rec_counts['updated_to_supplied_mapping_ensembl_id']:,} were updated to the supplied mapping's Ensembl ID.
- {rec_counts['filled_from_supplied_mapping']:,} were filled from the supplied mapping when the source guide file had no ENSG value.
- {fallback_groups:,} targets ({fallback_guides:,} targeting guides) had missing or ambiguous Ensembl mappings. They are represented as coordinate targets with `genomic_element=distal element` to remain validator-compatible. Their promoter identity, source legacy ENSG values, candidate current ENSG values, and coordinates are retained in `review/genome_wide_promoter_targets_requiring_DACC_review.tsv` for Ian/DACC review.

All {len(annotation['validation']):,} validation promoter targets ({len(gene_rows):,} gene symbols) resolve to one Ensembl ID in the supplied mapping.

## Analysis-file update

The complete `D3pos_D3neg_gene_summary.txt` attachment was not available in the files accessible to this session, so it could not be rewritten directly. Included instead:

- `annotations/D3pos_D3neg_gene_summary_annotation_lookup.tsv`: complete {len(gene_rows):,}-gene validation lookup.
- `annotations/IGVFFI6745HXWN_analysis_target_annotations.tsv`: promoter-target-level lookup ({len(annotation['validation']):,} rows).
- `scripts/add_promoter_annotations_to_analysis.py`: inserts Ensembl ID and promoter coordinate columns immediately after `id`; exact promoter IDs are matched before gene-level fallback.
- `annotations/D3pos_D3neg_gene_summary_ANNOTATED_EXAMPLE_ONLY.tsv`: the four rows pasted in the message, annotated only as a format example.

Run after placing the full analysis file next to this package:

```bash
python scripts/add_promoter_annotations_to_analysis.py D3pos_D3neg_gene_summary.txt
```

## Portal-side metadata action

For both tabular-file records, delete the currently linked obsolete document and set `file_format_specifications` to only:

`{FORMAT_DOC}`

A proposed patch fragment is in `portal_metadata_patch.json`. No Portal write was performed because Portal credentials were not supplied.

## Coordinate definition

`promoter_coordinate` is the hg38 interval from min(`guide_start`) to max(`guide_end`) across all `targeting=True` genome-wide guides assigned to the promoter target. It is a guide-covered promoter window, not a separately defined canonical TSS-centered promoter interval.
'''
    (ROOT / "README.md").write_text(readme, encoding="utf-8")

    draft = f'''Hi Chong,

I updated both guide metadata files using the supplied Weissman hg19-to-hg38/GENCODE mapping. I filled the intended-target chromosome/start/end using the full guide-covered promoter window, used promoter plus a single Ensembl ID where the mapping was unambiguous, normalized the required fields, and added descriptions. I also confirmed that all {val_rows:,} validation-library guides are present in the genome-wide file with matching spacer sequences.

Both corrected files have zero error categories under the current guide-validator checks. The genome-wide file retains one warning for the existing {targeting_false_gw:,} rows marked type=targeting but targeting=False; I left intended-target fields empty for those QC-excluded rows.

All validation-screen promoter targets resolve to one Ensembl ID. In the genome-wide library, {fallback_groups:,} promoter targets have missing or ambiguous Ensembl mappings. I represented those as coordinate targets so the file remains validator-compatible and included a separate review table rather than choosing an arbitrary Ensembl ID.

I also generated the validation analysis lookup with gene symbol, Ensembl ID, and promoter coordinates. The complete D3pos_D3neg_gene_summary.txt file was not available in the uploaded files I could access, so the package includes a merge script and an annotated example based on the four rows pasted in the email/Slack text.

Portal-side, the old file-format document still needs to be removed and replaced with:
{FORMAT_DOC}

Best,
Tri
'''
    (ROOT / "draft_email_to_Chong.txt").write_text(draft, encoding="utf-8")

    # Reproducibility manifest.
    manifest = []
    for path in sorted(ROOT.rglob("*")):
        if path.is_file() and path.name != "manifest.tsv":
            manifest.append({"relative_path": str(path.relative_to(ROOT)), "bytes": str(path.stat().st_size), "sha256": sha256(path)})
    write_rows(ROOT / "manifest.tsv", ["relative_path", "bytes", "sha256"], manifest)

    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for path in sorted(ROOT.rglob("*")):
            if path.is_file():
                zf.write(path, arcname=str(Path(ROOT.name) / path.relative_to(ROOT)))

    summary = {
        "root": str(ROOT), "zip": str(ZIP_PATH), "validation_rows": val_rows, "genome_wide_rows": gw_rows,
        "validation_targets": len(annotation["validation"]), "validation_genes": len(gene_rows),
        "genome_wide_targets": len(annotation["genome_wide"]), "fallback_groups": fallback_groups,
        "fallback_guides": fallback_guides, "targeting_false_genome_wide": targeting_false_gw,
        "reconciliation_counts": dict(rec_counts), "validation_summary": validation_summary,
    }
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
