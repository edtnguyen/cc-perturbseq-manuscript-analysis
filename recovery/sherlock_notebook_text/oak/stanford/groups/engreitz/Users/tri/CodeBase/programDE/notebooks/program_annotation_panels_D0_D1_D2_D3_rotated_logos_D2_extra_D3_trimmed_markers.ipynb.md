# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/program_annotation_panels_D0_D1_D2_D3_rotated_logos_D2_extra_D3_trimmed_markers.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/program_annotation_panels_D0_D1_D2_D3_rotated_logos_D2_extra_D3_trimmed_markers.ipynb`
- SHA-256: `b0ea563bee36408e54ff793e4915a7258232473a65813bba2bb5d7b5f655599a`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Program annotation panels across D0–D3

This notebook generates:

1. **One annotation panel per time point** for D0, D1, D2, and D3.
2. **One combined D0–D3 figure** with all four time points stacked vertically.

Each time-point panel contains three aligned blocks:

- **Regulator effects:** programs are rows and the union of the top three program-specific regulators is shown as columns.
- **Marker-gene loadings:** programs are rows and representative marker genes are shown as columns.
- **Motif enrichment:** programs remain rows and motifs remain columns. Each motif logo is displayed above its corresponding motif column and rotated by 90 degrees.

The combined figure uses exactly one shared color scale for each data type across all time points:

- regulator effects: shared symmetric approximate log2-fold-change scale;
- marker genes: shared consensus cNMF-loading scale;
- motif enrichment: shared `-log10(q-value)` scale.

Motifs are selected as the top three motifs per program regardless of significance; stars only indicate significance. Heatmap cells are opaque vector rectangles for Illustrator compatibility, and motif logos are vector `logomaker` glyphs.


## Marker-gene selection update

Marker-gene selection is now day-aware:

- **D0 and D1:** retain the curated representative genes from `Programs_annotations.xlsx`.
- **D2:** retain all curated representative genes and append the top three additional genes by consensus cNMF loading for each program, excluding duplicates.
- **D3:** rank the curated representative genes by their loading in the corresponding program and retain only the top three. Programs with only two curated genes retain both.

The number of added/retained genes is configurable in Cell 1.


## Confirmed Sherlock paths

The input paths in Cell 1 were synchronized to the attached reference notebook. In particular, motif inputs now use the `multiome_ipsc_ec_v0` tree:

- enhancer enrichment: `multiome_ipsc_ec_v0/analysis/2.tf_enrichment_in_programs/`
- cluster mapping: `multiome_ipsc_ec_v0/analysis/0.0.motif_id_annotation/tables/`
- compiled MoDISco H5: `multiome_ipsc_ec_v0/archive/_results/compendium/modisco_compiled/modisco_compiled.h5`

The notebook validates all required paths and read permissions before loading any data.


## Markdown cell 1

## Cell 1 — imports, paths, path validation, and settings

## Code cell 2

```python

import os
import re
import json
from pathlib import Path
from collections import OrderedDict

import numpy as np
import pandas as pd
import h5py
import logomaker

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from matplotlib.transforms import Affine2D


# ============================================================
# INPUT PATHS
# Confirmed from:
# program_annotation_panels_D0_D1_D2_D3_motifs_as_rows_
# combined_shared_scales(1).ipynb
# ============================================================

CONSENSUS_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50."
    "spectra.k_50.dt_0_2.consensus.txt"
)

PROGRAM_ANNOT_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/"
    "per_day_umaps/data/Programs_annotations.xlsx"
)

CNMF_WORKBOOK_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/"
    "per_day_umaps/data/cNMF_50_0_2_032526_6co.xlsx"
)

ENHANCER_ENRICHMENT_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec_v0/"
    "analysis/2.tf_enrichment_in_programs/enhancer_enrichment_results.tsv"
)

CLUSTER_MAP_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec_v0/"
    "analysis/0.0.motif_id_annotation/tables/"
    "cluster_to_tf_motif_id_mapping.tsv"
)

MODISCO_H5_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec_v0/archive/"
    "_results/compendium/modisco_compiled/modisco_compiled.h5"
)

OUTDIR = Path(
    "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/"
    "per_day_umaps"
)

SOURCE_DIR = OUTDIR / "_program_annotation_source"

OUTDIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# VALIDATE CONFIRMED INPUT PATHS BEFORE RUNNING ANALYSIS
# ============================================================

INPUT_PATHS = {
    "cNMF consensus matrix": CONSENSUS_PATH,
    "program annotations": PROGRAM_ANNOT_PATH,
    "cNMF workbook": CNMF_WORKBOOK_PATH,
    "enhancer enrichment": ENHANCER_ENRICHMENT_PATH,
    "motif cluster mapping": CLUSTER_MAP_PATH,
    "MoDISco compiled H5": MODISCO_H5_PATH,
}

print("Input-file validation:")

unavailable_inputs = []

for label, path in INPUT_PATHS.items():
    exists = path.is_file()
    readable = os.access(path, os.R_OK) if exists else False

    print(
        f"{label:28s} "
        f"exists={str(exists):5s} "
        f"readable={str(readable):5s}\n"
        f"    {path}"
    )

    if not exists or not readable:
        unavailable_inputs.append(
            (label, path, exists, readable)
        )

if unavailable_inputs:
    details = "\n".join(
        f"{label}: exists={exists}, readable={readable}, path={path}"
        for label, path, exists, readable in unavailable_inputs
    )

    raise FileNotFoundError(
        "One or more confirmed input paths are unavailable:\n"
        + details
    )

print("\nOutput directory:")
print(OUTDIR)


# ============================================================
# FIGURE SETTINGS
# ============================================================

DAYS_TO_PLOT = ["D0", "D1", "D2", "D3"]

N_TOP_REGULATORS = 3
N_TOP_MOTIFS = 3


# ============================================================
# DAY-SPECIFIC MARKER-GENE SELECTION
# ============================================================

# D2: curated genes plus this many additional top-loading genes per program.
D2_N_ADDITIONAL_TOP_LOADING_GENES = 3

# D3: among curated genes, keep at most this many with the highest
# program-specific consensus loading. If only two curated genes exist,
# both are retained.
D3_N_CURATED_MARKERS_TO_KEEP = 3

# D2 output order within each program's marker group:
# True  = curated genes first, then additional top-loading genes
# False = additional top-loading genes first, then curated genes
D2_CURATED_GENES_FIRST = True


REGULATOR_FDR_THRESHOLD = 0.05
MOTIF_QVAL_THRESHOLD = 0.05

PROGRAM_OFFSET = "auto"

VECTOR_LOGOS = True
ROTATE_MOTIF_LOGOS_DEGREES = 90

MANUAL_CLUSTER_ID_COL = None
MANUAL_MOTIF_NAME_COLS = None


# ============================================================
# OUTPUT AND SHARED-SCALE SETTINGS
# ============================================================

MAKE_PER_DAY_PANELS = True
MAKE_COMBINED_ALL_DAYS_FIGURE = True

# None = calculate from every value plotted across D0-D3.
GLOBAL_REGULATOR_VMAX_OVERRIDE = None
GLOBAL_MARKER_VMAX_OVERRIDE = None
GLOBAL_MOTIF_VMAX_OVERRIDE = None

# Per-day geometry
PER_DAY_CELL_W_IN = 0.34
PER_DAY_CELL_H_IN = 0.34
PER_DAY_BLOCK_GAP_IN = 0.65
PER_DAY_LEFT_MARGIN_IN = 3.3
PER_DAY_RIGHT_MARGIN_IN = 0.8
PER_DAY_BOTTOM_MARGIN_IN = 2.0
PER_DAY_TOP_MARGIN_IN = 1.1
PER_DAY_XLABEL_SPACE_IN = 1.25

# Combined figure geometry
COMBINED_CELL_W_IN = 0.30
COMBINED_CELL_H_IN = 0.30
COMBINED_DAY_GAP_IN = 0.75
COMBINED_HEADER_IN = 1.9
COMBINED_DAY_LABEL_SPACE_IN = 0.42
COMBINED_XLABEL_SPACE_IN = 1.30
COMBINED_LEFT_MARGIN_IN = 3.35
COMBINED_RIGHT_MARGIN_IN = 0.9
COMBINED_BLOCK_GAP_IN = 0.65

# Rotated motif logo geometry
ROTATED_LOGO_HEIGHT_IN = 1.05
ROTATED_LOGO_PAD_IN = 0.12

HEATMAP_TILE_LINEWIDTH = 1.0
MISSING_TILE_COLOR = (0.95, 0.95, 0.95, 1.0)

REGULATOR_CMAP = mpl.colormaps["RdBu_r"]
MARKER_CMAP = mpl.colormaps["Greys"]
MOTIF_CMAP = mpl.colormaps["Purples"]


# ============================================================
# ILLUSTRATOR-SAFE MATPLOTLIB DEFAULTS
# ============================================================

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.transparent": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "font.family": "DejaVu Sans",
})


# ============================================================
# VALIDATE INPUTS BEFORE ANALYSIS
# ============================================================

INPUT_PATHS = {
    "cNMF consensus matrix": CONSENSUS_PATH,
    "program annotations": PROGRAM_ANNOT_PATH,
    "cNMF workbook": CNMF_WORKBOOK_PATH,
    "enhancer enrichment": ENHANCER_ENRICHMENT_PATH,
    "motif cluster mapping": CLUSTER_MAP_PATH,
    "MoDISco compiled H5": MODISCO_H5_PATH,
}

print("Input-file validation:")

unavailable = []

for label, path in INPUT_PATHS.items():
    exists = path.is_file()
    readable = os.access(path, os.R_OK) if exists else False

    print(
        f"{label:28s} "
        f"exists={str(exists):5s} "
        f"readable={str(readable):5s}\n"
        f"    {path}"
    )

    if not exists or not readable:
        unavailable.append(
            (label, path, exists, readable)
        )

if unavailable:
    details = "\n".join(
        f"{label}: exists={exists}, readable={readable}, path={path}"
        for label, path, exists, readable in unavailable
    )

    raise FileNotFoundError(
        "One or more required inputs are unavailable:\n"
        + details
    )

print("\nOutput directory:")
print(OUTDIR)

```

### Text output (stream)

```text
Input-file validation:
cNMF consensus matrix        exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50.spectra.k_50.dt_0_2.consensus.txt
program annotations          exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/data/Programs_annotations.xlsx
cNMF workbook                exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/data/cNMF_50_0_2_032526_6co.xlsx
enhancer enrichment          exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec_v0/analysis/2.tf_enrichment_in_programs/enhancer_enrichment_results.tsv
motif cluster mapping        exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec_v0/analysis/0.0.motif_id_annotation/tables/cluster_to_tf_motif_id_mapping.tsv
MoDISco compiled H5          exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec_v0/archive/_results/compendium/modisco_compiled/modisco_compiled.h5

Output directory:
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps
Input-file validation:
cNMF consensus matrix        exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50.spectra.k_50.dt_0_2.consensus.txt
program annotations          exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/data/Programs_annotations.xlsx
cNMF workbook                exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/data/cNMF_50_0_2_032526_6co.xlsx
enhancer enrichment          exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec_v0/analysis/2.tf_enrichment_in_programs/enhancer_enrichment_results.tsv
motif cluster mapping        exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec_v0/analysis/0.0.motif_id_annotation/tables/cluster_to_tf_motif_id_mapping.tsv
MoDISco compiled H5          exists=True  readable=True 
    /oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec_v0/archive/_results/compendium/modisco_compiled/modisco_compiled.h5

Output directory:
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps

```

## Markdown cell 3

## Cell 2 — metadata, regulator, and consensus-loading helper functions

## Code cell 4

```python

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def normalize_token(x):
    if pd.isna(x):
        return None
    x = str(x).strip().upper()
    x = re.sub(r"[^A-Z0-9]+", "", x)
    return x if x else None


def parse_program_num(x):
    m = re.search(r"(\d+)", str(x))
    return int(m.group(1)) if m else None


def split_feature_text(x):
    """
    Split annotation strings like:
      'TAL1/FLI1'
      'MECOM; NRP1; EFNB2'
      'L1TD1, DPPA4, TERF1, SOX2'
    into ordered feature tokens.
    """
    if pd.isna(x):
        return []

    txt = str(x).strip()
    if txt == "":
        return []

    parts = re.split(r"[;,/]", txt)
    cleaned = []
    seen = set()

    for part in parts:
        item = str(part).strip()
        item = re.sub(r"\([^)]*\)", "", item).strip()
        item = item.replace("?", "").strip()
        item = re.sub(r"\s+", " ", item)
        if item == "":
            continue

        # keep first token if there are comments
        item = item.split(" ")[0].strip()
        key = normalize_token(item)

        if key and key not in seen:
            cleaned.append(item)
            seen.add(key)

    return cleaned


def star_text(p):
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def read_matrix_flex(path, index_col=0):
    attempts = [
        {"sep": "\t"},
        {"sep": r"\s+", "engine": "python"},
        {"sep": ","},
    ]
    errors = []

    for kwargs in attempts:
        try:
            df = pd.read_csv(path, index_col=index_col, **kwargs)
            df = df.loc[:, ~df.columns.astype(str).str.startswith("Unnamed")]
            if df.shape[0] >= 2 and df.shape[1] >= 2:
                return df
        except Exception as e:
            errors.append(str(e))

    raise ValueError(
        f"Could not read matrix: {path}\n" + "\n".join(errors)
    )


# ============================================================
# CONSENSUS LOADING MATRIX
# ============================================================

def read_consensus_matrix(path):
    """
    Return consensus matrix oriented as programs x genes.
    """
    mat = read_matrix_flex(path, index_col=0)

    # likely already ~50 x many_genes
    if mat.shape[0] > 200 and mat.shape[1] <= 200:
        mat = mat.T
    elif not (mat.shape[0] <= 200 and mat.shape[1] > 200):
        idx_prog_frac = np.mean([
            bool(re.fullmatch(r"\s*(?:P|Program)?\s*\d+(?:\.0)?\s*", str(x), flags=re.I))
            for x in mat.index
        ])
        col_prog_frac = np.mean([
            bool(re.fullmatch(r"\s*(?:P|Program)?\s*\d+(?:\.0)?\s*", str(x), flags=re.I))
            for x in mat.columns
        ])
        if col_prog_frac > idx_prog_frac:
            mat = mat.T

    mat.columns = mat.columns.astype(str).str.strip()

    if mat.columns.duplicated().any():
        print("[warning] Duplicate gene columns found; collapsing by mean.")
        mat = mat.T.groupby(level=0).mean().T

    return mat.apply(pd.to_numeric, errors="coerce")


def build_program_row_lookup(consensus_df):
    lookup = {}
    for row_label in consensus_df.index:
        p = parse_program_num(row_label)
        if p is not None and p not in lookup:
            lookup[p] = row_label
    return lookup


def choose_program_offset(consensus_df, program_meta):
    """
    Use marker-gene loading consistency to choose offset 0 / +1 / -1.
    """
    if PROGRAM_OFFSET != "auto":
        print("Using manual PROGRAM_OFFSET =", PROGRAM_OFFSET)
        return int(PROGRAM_OFFSET)

    row_lookup = build_program_row_lookup(consensus_df)
    gene_lookup = {normalize_token(g): g for g in consensus_df.columns}

    def score_offset(offset):
        vals = []
        for _, row in program_meta.iterrows():
            p = int(row["program_num"])
            marker_genes = row["marker_genes"]
            raw_p = p + offset
            if raw_p not in row_lookup:
                continue
            raw_row = row_lookup[raw_p]
            for gene in marker_genes:
                g = gene_lookup.get(normalize_token(gene))
                if g is None:
                    continue
                val = consensus_df.loc[raw_row, g]
                if pd.notna(val):
                    vals.append(float(val))
        if len(vals) == 0:
            return -np.inf
        return float(np.mean(vals))

    scores = {off: score_offset(off) for off in [0, 1, -1]}
    print("Consensus program-offset scores:", scores)
    best = max(scores, key=scores.get)
    print("Chosen consensus program offset:", best)
    return best


# ============================================================
# PROGRAM METADATA
# ============================================================

def load_program_metadata():
    """
    Merge:
      Programs_annotations.xlsx: program labels + representative genes
      cNMF Summary: automatic timepoint + top specific regulators columns
    """
    # Programs_annotations.xlsx
    annot = pd.read_excel(PROGRAM_ANNOT_PATH)
    annot.columns = [str(c).strip() for c in annot.columns]

    # Auto-detect expected annotation columns
    prog_col = next((c for c in annot.columns if "program" in c.lower() and "number" in c.lower()), annot.columns[0])
    name_col = next((c for c in annot.columns if c.lower() in ["name", "program name"]), None)
    rep_col = next((c for c in annot.columns if "representative" in c.lower() and "gene" in c.lower()), None)

    if name_col is None:
        name_col = next((c for c in annot.columns if "name" in c.lower()), None)
    if rep_col is None:
        rep_col = next((c for c in annot.columns if "gene" in c.lower()), None)

    if name_col is None or rep_col is None:
        raise ValueError(
            "Could not detect Name / Representative genes columns in Programs_annotations.xlsx. "
            f"Columns: {annot.columns.tolist()}"
        )

    annot = annot.rename(columns={
        prog_col: "program_raw",
        name_col: "program_name",
        rep_col: "representative_genes",
    })
    annot["program_num"] = annot["program_raw"].map(parse_program_num)
    annot = annot.dropna(subset=["program_num"]).copy()
    annot["program_num"] = annot["program_num"].astype(int)
    annot["program_name"] = annot["program_name"].astype(str).str.strip()
    annot["marker_genes"] = annot["representative_genes"].map(split_feature_text)

    # cNMF Summary
    summary = pd.read_excel(CNMF_WORKBOOK_PATH, sheet_name="Summary")
    summary.columns = [str(c).strip() for c in summary.columns]
    summary = summary.rename(columns={
        "program_name": "program_num",
        "Automatic Timepoint": "day",
    })
    summary["program_num"] = summary["program_num"].map(parse_program_num)
    summary = summary.dropna(subset=["program_num"]).copy()
    summary["program_num"] = summary["program_num"].astype(int)
    summary["day"] = summary["day"].astype(str).str.strip().str.upper()

    # merge
    meta = annot.merge(summary, on="program_num", how="left")

    # label
    meta["row_label"] = meta.apply(
        lambda r: f"P{int(r['program_num'])}  {r['program_name']}",
        axis=1
    )

    meta["marker_gene_source"] = "Programs_annotations.xlsx"
    meta = meta.sort_values(["day", "program_num"]).reset_index(drop=True)

    return meta


# ============================================================
# REGULATORS
# ============================================================

def load_perturbation_associations():
    pert = pd.read_excel(CNMF_WORKBOOK_PATH, sheet_name="Perturbation Association 1")
    pert.columns = [str(c).strip() for c in pert.columns]

    rename = {}
    if "target_name" in pert.columns:
        rename["target_name"] = "gene"
    if "program_name" in pert.columns:
        rename["program_name"] = "program"
    if "Sample" in pert.columns:
        rename["Sample"] = "day"
    if "approx_log2FC" in pert.columns:
        rename["approx_log2FC"] = "log2fc"
    if "adj_pval" in pert.columns:
        rename["adj_pval"] = "fdr"

    # Alternate spellings
    if "Perturbed gene" in pert.columns:
        rename["Perturbed gene"] = "gene"
    if "Gene Program" in pert.columns:
        rename["Gene Program"] = "program"
    if "adj.P.Val" in pert.columns:
        rename["adj.P.Val"] = "fdr"

    pert = pert.rename(columns=rename)

    required = ["gene", "program", "day", "log2fc", "fdr"]
    missing = [c for c in required if c not in pert.columns]
    if missing:
        raise ValueError(
            f"Perturbation Association 1 missing columns after rename: {missing}\n"
            f"Available columns: {pert.columns.tolist()}"
        )

    pert["program_num"] = pert["program"].map(parse_program_num)
    pert["day"] = pert["day"].astype(str).str.strip().str.upper()
    pert["gene_display"] = pert["gene"].astype(str).str.strip()
    pert["gene_key"] = pert["gene_display"].map(normalize_token)
    pert["log2fc"] = pd.to_numeric(pert["log2fc"], errors="coerce")
    pert["fdr"] = pd.to_numeric(pert["fdr"], errors="coerce")
    pert["is_sig"] = pert["fdr"] < REGULATOR_FDR_THRESHOLD

    return pert


def select_top_regulators_for_programs(program_meta_day, pert):
    """
    Top 3 specific regulators per program:
      1) use cNMF Summary column: 'Top 5 specific regulators (FDR<0.1) D2/D3'
      2) keep first N_TOP_REGULATORS from that list
      3) retrieve log2FC/FDR from Perturbation Association 1
    """
    day = str(program_meta_day["day"].iloc[0]).upper()
    col = f"Top 5 specific regulators (FDR<0.1) {day}"
    if col not in program_meta_day.columns:
        raise ValueError(f"Missing Summary column: {col}")

    selected_by_program = OrderedDict()
    rows = []

    for _, prow in program_meta_day.iterrows():
        p = int(prow["program_num"])
        regulator_list = split_feature_text(prow.get(col, ""))[:N_TOP_REGULATORS]
        selected_by_program[p] = regulator_list

        for rank, gene in enumerate(regulator_list, start=1):
            sub = pert.loc[
                (pert["program_num"] == p)
                & (pert["day"] == day)
                & (pert["gene_key"] == normalize_token(gene))
            ].copy()

            if len(sub) == 0:
                log2fc = np.nan
                fdr = np.nan
            else:
                # If duplicated, keep lowest FDR
                sub = sub.sort_values("fdr")
                log2fc = sub["log2fc"].iloc[0]
                fdr = sub["fdr"].iloc[0]

            rows.append({
                "day": day,
                "program_num": p,
                "program_label": prow["row_label"],
                "rank_within_program": rank,
                "gene": gene,
                "log2fc": log2fc,
                "fdr": fdr,
                "signif_star": star_text(fdr),
            })

    selected_df = pd.DataFrame(rows)

    # union across programs in program order
    union = []
    seen = set()
    for p in program_meta_day["program_num"]:
        for gene in selected_by_program[int(p)]:
            key = normalize_token(gene)
            if key not in seen:
                union.append(gene)
                seen.add(key)

    return selected_by_program, union, selected_df


def build_regulator_matrices(program_meta_day, regulator_union, pert):
    day = str(program_meta_day["day"].iloc[0]).upper()
    rows = [int(x) for x in program_meta_day["program_num"]]

    sub = pert.loc[pert["day"] == day].copy()

    effect_pivot = sub.pivot_table(
        index="program_num",
        columns="gene_key",
        values="log2fc",
        aggfunc="first"
    )
    fdr_pivot = sub.pivot_table(
        index="program_num",
        columns="gene_key",
        values="fdr",
        aggfunc="min"
    )

    reg_keys = [normalize_token(g) for g in regulator_union]

    effect_mat = pd.DataFrame(index=rows, columns=regulator_union, dtype=float)
    fdr_mat = pd.DataFrame(index=rows, columns=regulator_union, dtype=float)
    star_mat = pd.DataFrame(index=rows, columns=regulator_union, dtype=object)

    for disp, key in zip(regulator_union, reg_keys):
        for p in rows:
            val = np.nan
            q = np.nan
            if p in effect_pivot.index and key in effect_pivot.columns:
                val = effect_pivot.loc[p, key]
            if p in fdr_pivot.index and key in fdr_pivot.columns:
                q = fdr_pivot.loc[p, key]

            effect_mat.loc[p, disp] = val
            fdr_mat.loc[p, disp] = q
            star_mat.loc[p, disp] = star_text(q)

    return effect_mat, fdr_mat, star_mat

```

## Markdown cell 5

## Cell 3 — day-aware marker selection, complete motif selection, and cached vector-logo helpers

## Code cell 6

```python

# ============================================================
# MARKER GENES AND MOTIFS
# ============================================================

def select_marker_genes_for_programs(
    program_meta_day,
    consensus_df,
    consensus_offset,
):
    """
    Day-aware marker-gene selection.

    D0 / D1
        Use curated representative genes from Programs_annotations.xlsx.

    D2
        Use all curated genes plus the top
        D2_N_ADDITIONAL_TOP_LOADING_GENES additional genes ranked by
        consensus cNMF loading in that program. Curated genes and
        duplicate genes are excluded from the added set.

    D3
        Rank only the curated genes by their loading in that program and
        keep at most D3_N_CURATED_MARKERS_TO_KEEP. Thus a program with
        two curated genes retains both, while a program with >3 retains
        the three highest-loading curated genes.

    Returns
    -------
    selected_by_program : OrderedDict[int, list[str]]
    union : list[str]
        Union across programs in program order.
    selected_df : pandas.DataFrame
        One row per selected marker, including its source and loading.
    """
    if program_meta_day.empty:
        return OrderedDict(), [], pd.DataFrame()

    day = str(
        program_meta_day["day"].iloc[0]
    ).strip().upper()

    row_lookup = build_program_row_lookup(
        consensus_df
    )

    gene_lookup = {
        normalize_token(gene): gene
        for gene in consensus_df.columns
    }

    selected_by_program = OrderedDict()
    metadata_rows = []

    def get_program_row(program_num):
        raw_program = (
            int(program_num)
            + int(consensus_offset)
        )

        if raw_program not in row_lookup:
            return None, None

        raw_row_label = row_lookup[raw_program]

        row = pd.to_numeric(
            consensus_df.loc[raw_row_label],
            errors="coerce",
        )

        return raw_row_label, row

    def loading_for_gene(row, requested_gene):
        matrix_gene = gene_lookup.get(
            normalize_token(requested_gene)
        )

        if matrix_gene is None:
            return np.nan, None

        value = row.get(
            matrix_gene,
            np.nan,
        )

        return (
            float(value)
            if pd.notna(value)
            else np.nan,
            matrix_gene,
        )

    for _, program_row in program_meta_day.iterrows():
        program_num = int(
            program_row["program_num"]
        )

        curated_genes = list(
            program_row["marker_genes"]
        )

        raw_row_label, loading_row = get_program_row(
            program_num
        )

        if loading_row is None:
            print(
                f"[warning] No consensus row found for "
                f"P{program_num}; retaining curated markers "
                f"without loading-based modification."
            )

            selected_genes = list(
                curated_genes
            )

            selected_by_program[
                program_num
            ] = selected_genes

            for rank, gene in enumerate(
                selected_genes,
                start=1,
            ):
                metadata_rows.append({
                    "day": day,
                    "program_num": program_num,
                    "program_label": program_row["row_label"],
                    "rank_within_program": rank,
                    "gene": gene,
                    "marker_source": "curated",
                    "selection_rule": "fallback_curated_no_consensus_row",
                    "loading_in_source_program": np.nan,
                    "matrix_gene": None,
                    "raw_consensus_row_label": None,
                })

            continue

        curated_records = []

        for curated_index, gene in enumerate(
            curated_genes,
            start=1,
        ):
            loading, matrix_gene = loading_for_gene(
                loading_row,
                gene,
            )

            curated_records.append({
                "gene": gene,
                "curated_order": curated_index,
                "loading": loading,
                "matrix_gene": matrix_gene,
            })

        # ------------------------------------------------------
        # D2: curated markers + additional top-loading genes
        # ------------------------------------------------------
        if day == "D2":
            curated_keys = {
                normalize_token(gene)
                for gene in curated_genes
            }

            ranked_loading_row = (
                loading_row
                .dropna()
                .sort_values(
                    ascending=False
                )
            )

            additional_genes = []

            for matrix_gene, loading in ranked_loading_row.items():
                display_gene = str(
                    matrix_gene
                ).strip()

                gene_key = normalize_token(
                    display_gene
                )

                if gene_key in curated_keys:
                    continue

                if gene_key in {
                    normalize_token(gene)
                    for gene in additional_genes
                }:
                    continue

                additional_genes.append(
                    display_gene
                )

                if (
                    len(additional_genes)
                    >= D2_N_ADDITIONAL_TOP_LOADING_GENES
                ):
                    break

            if D2_CURATED_GENES_FIRST:
                selected_genes = (
                    curated_genes
                    + additional_genes
                )
            else:
                selected_genes = (
                    additional_genes
                    + curated_genes
                )

            selected_by_program[
                program_num
            ] = selected_genes

            selected_rank = 0

            for gene in selected_genes:
                selected_rank += 1

                curated_match = next(
                    (
                        record
                        for record in curated_records
                        if normalize_token(
                            record["gene"]
                        ) == normalize_token(
                            gene
                        )
                    ),
                    None,
                )

                if curated_match is not None:
                    marker_source = "curated"
                    loading = curated_match["loading"]
                    matrix_gene = curated_match["matrix_gene"]
                    selection_rule = (
                        "D2_curated_plus_top_loading"
                    )
                else:
                    marker_source = (
                        "additional_top_loading"
                    )
                    loading, matrix_gene = loading_for_gene(
                        loading_row,
                        gene,
                    )
                    selection_rule = (
                        "D2_curated_plus_top_loading"
                    )

                metadata_rows.append({
                    "day": day,
                    "program_num": program_num,
                    "program_label": program_row["row_label"],
                    "rank_within_program": selected_rank,
                    "gene": gene,
                    "marker_source": marker_source,
                    "selection_rule": selection_rule,
                    "loading_in_source_program": loading,
                    "matrix_gene": matrix_gene,
                    "raw_consensus_row_label": raw_row_label,
                })

        # ------------------------------------------------------
        # D3: top-loading curated markers only
        # ------------------------------------------------------
        elif day == "D3":
            curated_sort_df = pd.DataFrame(
                curated_records
            )

            if curated_sort_df.empty:
                selected_genes = []

            else:
                # Genes found in the consensus matrix are ranked
                # by loading. Missing genes are placed last while
                # preserving curated order among ties.
                curated_sort_df[
                    "loading_sort"
                ] = curated_sort_df[
                    "loading"
                ].fillna(-np.inf)

                curated_sort_df = (
                    curated_sort_df
                    .sort_values(
                        by=[
                            "loading_sort",
                            "curated_order",
                        ],
                        ascending=[
                            False,
                            True,
                        ],
                    )
                )

                selected_records = (
                    curated_sort_df
                    .head(
                        D3_N_CURATED_MARKERS_TO_KEEP
                    )
                    .copy()
                )

                selected_genes = (
                    selected_records["gene"]
                    .tolist()
                )

            selected_by_program[
                program_num
            ] = selected_genes

            selected_record_lookup = {
                normalize_token(record["gene"]): record
                for record in curated_records
            }

            for rank, gene in enumerate(
                selected_genes,
                start=1,
            ):
                record = selected_record_lookup[
                    normalize_token(gene)
                ]

                metadata_rows.append({
                    "day": day,
                    "program_num": program_num,
                    "program_label": program_row["row_label"],
                    "rank_within_program": rank,
                    "gene": gene,
                    "marker_source": "curated",
                    "selection_rule": "D3_top_loading_curated_only",
                    "loading_in_source_program": record["loading"],
                    "matrix_gene": record["matrix_gene"],
                    "raw_consensus_row_label": raw_row_label,
                })

        # ------------------------------------------------------
        # D0 / D1: curated markers unchanged
        # ------------------------------------------------------
        else:
            selected_genes = list(
                curated_genes
            )

            selected_by_program[
                program_num
            ] = selected_genes

            record_lookup = {
                normalize_token(record["gene"]): record
                for record in curated_records
            }

            for rank, gene in enumerate(
                selected_genes,
                start=1,
            ):
                record = record_lookup[
                    normalize_token(gene)
                ]

                metadata_rows.append({
                    "day": day,
                    "program_num": program_num,
                    "program_label": program_row["row_label"],
                    "rank_within_program": rank,
                    "gene": gene,
                    "marker_source": "curated",
                    "selection_rule": f"{day}_curated_unchanged",
                    "loading_in_source_program": record["loading"],
                    "matrix_gene": record["matrix_gene"],
                    "raw_consensus_row_label": raw_row_label,
                })

    selected_df = pd.DataFrame(
        metadata_rows
    )

    # Union across programs in program order.
    union = []
    seen = set()

    for program_num in program_meta_day[
        "program_num"
    ]:
        for gene in selected_by_program.get(
            int(program_num),
            [],
        ):
            gene_key = normalize_token(
                gene
            )

            if gene_key not in seen:
                union.append(
                    gene
                )
                seen.add(
                    gene_key
                )

    return (
        selected_by_program,
        union,
        selected_df,
    )


def build_marker_loading_matrix(program_meta_day, marker_union, consensus_df, consensus_offset):
    row_lookup = build_program_row_lookup(consensus_df)
    gene_lookup = {normalize_token(g): g for g in consensus_df.columns}

    rows = [int(x) for x in program_meta_day["program_num"]]
    mat = pd.DataFrame(index=rows, columns=marker_union, dtype=float)

    for p in rows:
        raw_p = p + consensus_offset
        if raw_p not in row_lookup:
            continue
        raw_row = row_lookup[raw_p]

        for gene in marker_union:
            g = gene_lookup.get(normalize_token(gene))
            if g is not None:
                mat.loc[p, gene] = consensus_df.loc[raw_row, g]

    return mat


def load_enhancer_enrichment():
    df = pd.read_csv(ENHANCER_ENRICHMENT_PATH, sep="\t")
    df["program_num"] = pd.to_numeric(df["program"], errors="coerce").astype("Int64")
    df["tf_display"] = df["tf"].astype(str).str.strip()
    df["tf_key"] = df["tf_display"].map(normalize_token)
    df["day"] = df["day"].astype(str).str.strip().str.upper()
    df["qval"] = pd.to_numeric(df["qval"], errors="coerce")
    df["neg_log10_qval"] = pd.to_numeric(df["neg_log10_qval"], errors="coerce")
    return df


def select_top_motifs_for_programs(program_meta_day, enrich):
    """
    Complete annotation-panel motif selection.

    For every program:
      - select the top N_TOP_MOTIFS motifs by q-value
      - do NOT require q < MOTIF_QVAL_THRESHOLD
      - keep non-significant motifs in the figure
      - use stars only to mark significance

    This guarantees that programs such as P7/P22 at D0 and P28 at D1
    still contribute motif columns/logos even when no motif passes q < 0.05.
    """
    selected_by_program = OrderedDict()
    rows = []

    for _, prow in program_meta_day.iterrows():
        p = int(prow["program_num"])
        day = str(prow["day"]).strip().upper()

        sub = enrich.loc[
            (enrich["program_num"] == p)
            & (enrich["day"] == day)
        ].copy()

        if sub.empty:
            selected_by_program[p] = []

            rows.append({
                "day": day,
                "program_num": p,
                "program_label": prow["row_label"],
                "rank_within_program": np.nan,
                "motif": np.nan,
                "qval": np.nan,
                "neg_log10_qval": np.nan,
                "AUC": np.nan,
                "significant": False,
                "signif_star": "",
                "note": "no motif enrichment rows found",
            })

            continue

        # Make sure numeric columns are numeric
        sub["qval"] = pd.to_numeric(sub["qval"], errors="coerce")
        sub["neg_log10_qval"] = pd.to_numeric(sub["neg_log10_qval"], errors="coerce")

        if "AUC" in sub.columns:
            sub["AUC"] = pd.to_numeric(sub["AUC"], errors="coerce")
        else:
            sub["AUC"] = np.nan

        # Deduplicate motifs, keeping strongest q-value
        sub = (
            sub.sort_values(
                by=["qval", "neg_log10_qval", "AUC"],
                ascending=[True, False, False],
            )
            .drop_duplicates("tf_key", keep="first")
        )

        # Complete-panel behavior:
        # take top motifs regardless of significance
        top = sub.head(N_TOP_MOTIFS).copy()

        motifs = top["tf_display"].tolist()
        selected_by_program[p] = motifs

        for rank, (_, r) in enumerate(top.iterrows(), start=1):
            q = r["qval"]

            rows.append({
                "day": day,
                "program_num": p,
                "program_label": prow["row_label"],
                "rank_within_program": rank,
                "motif": r["tf_display"],
                "qval": q,
                "neg_log10_qval": r["neg_log10_qval"],
                "AUC": r.get("AUC", np.nan),
                "significant": bool(pd.notna(q) and q < MOTIF_QVAL_THRESHOLD),
                "signif_star": star_text(q),
                "note": "top motif by q-value; significance not required",
            })

    selected_df = pd.DataFrame(rows)

    # Union motifs across programs in program order
    union = []
    seen = set()

    for p in program_meta_day["program_num"]:
        for motif in selected_by_program.get(int(p), []):
            key = normalize_token(motif)

            if key not in seen:
                union.append(motif)
                seen.add(key)

    return selected_by_program, union, selected_df


def build_motif_matrices(program_meta_day, motif_union, enrich):
    rows = [int(x) for x in program_meta_day["program_num"]]
    day = str(program_meta_day["day"].iloc[0]).upper()
    sub = enrich.loc[enrich["day"] == day].copy()

    value_pivot = sub.pivot_table(
        index="program_num",
        columns="tf_key",
        values="neg_log10_qval",
        aggfunc="first"
    )
    qval_pivot = sub.pivot_table(
        index="program_num",
        columns="tf_key",
        values="qval",
        aggfunc="min"
    )

    motif_keys = [normalize_token(m) for m in motif_union]

    value_mat = pd.DataFrame(index=rows, columns=motif_union, dtype=float)
    qval_mat = pd.DataFrame(index=rows, columns=motif_union, dtype=float)
    star_mat = pd.DataFrame(index=rows, columns=motif_union, dtype=object)

    for disp, key in zip(motif_union, motif_keys):
        for p in rows:
            val = np.nan
            qv = np.nan

            if p in value_pivot.index and key in value_pivot.columns:
                val = value_pivot.loc[p, key]
            if p in qval_pivot.index and key in qval_pivot.columns:
                qv = qval_pivot.loc[p, key]

            value_mat.loc[p, disp] = val
            qval_mat.loc[p, disp] = qv
            star_mat.loc[p, disp] = star_text(qv)

    return value_mat, qval_mat, star_mat



# ============================================================
# MOTIF LOGOS: vector logomaker glyphs
# ============================================================

def infer_cluster_map_columns(df):
    if MANUAL_CLUSTER_ID_COL is not None:
        id_col = MANUAL_CLUSTER_ID_COL
    else:
        candidates = []
        for c in df.columns:
            lc = c.lower()
            numeric_frac = pd.to_numeric(df[c], errors="coerce").notna().mean()
            if numeric_frac > 0.80 and re.search(r"(cluster|pattern|cid|id)", lc):
                candidates.append(c)
        if not candidates:
            for c in df.columns:
                numeric_frac = pd.to_numeric(df[c], errors="coerce").notna().mean()
                if numeric_frac > 0.80:
                    candidates.append(c)
        if len(candidates) == 0:
            raise ValueError(
                "Could not auto-detect cluster ID column. "
                f"Available columns: {list(df.columns)}"
            )
        id_col = candidates[0]

    if MANUAL_MOTIF_NAME_COLS is not None:
        name_cols = list(MANUAL_MOTIF_NAME_COLS)
    else:
        name_cols = []
        for c in df.columns:
            if c == id_col:
                continue
            if df[c].dtype == object:
                if not re.search(r"(path|file|url|image)", c.lower()):
                    name_cols.append(c)
        if len(name_cols) == 0:
            raise ValueError(
                "Could not auto-detect motif name columns. "
                f"Available columns: {list(df.columns)}"
            )

    return id_col, name_cols


def motif_alias_variants(name):
    base = str(name).strip()
    variants = set([base])

    cleaned = re.sub(r"_[0-9]+_$", "", base)
    variants.add(cleaned)

    for x in list(variants):
        variants.add(x.replace("_", ""))
        variants.add(x.replace("-", ""))
        variants.add(x.replace("/", ""))

    for part in re.split(r"[;/|]", base):
        part = part.strip()
        if part:
            variants.add(part)

    return [v for v in variants if v]


def build_cluster_lookup(cluster_map_df):
    id_col, name_cols = infer_cluster_map_columns(cluster_map_df)
    print("Cluster-map ID column:", id_col)
    print("Cluster-map motif columns:", name_cols)

    alias_to_cluster = {}

    for _, row in cluster_map_df.iterrows():
        cid = pd.to_numeric(row[id_col], errors="coerce")
        if pd.isna(cid):
            continue
        cid = int(cid)

        aliases = set()
        for c in name_cols:
            val = row[c]
            if pd.isna(val):
                continue
            aliases.add(str(val))
            for tok in re.split(r"[;,/|]", str(val)):
                tok = tok.strip()
                if tok:
                    aliases.add(tok)

        for alias in aliases:
            key = normalize_token(alias)
            if key and key not in alias_to_cluster:
                alias_to_cluster[key] = cid

    return alias_to_cluster


def best_cluster_for_motif(motif_name, alias_to_cluster):
    # exact match on normalized variants
    for variant in motif_alias_variants(motif_name):
        key = normalize_token(variant)
        if key in alias_to_cluster:
            return alias_to_cluster[key]

    # fuzzy prefix/containment match
    motif_keys = [normalize_token(v) for v in motif_alias_variants(motif_name)]
    motif_keys = [k for k in motif_keys if k]

    for mk in motif_keys:
        matches = []
        for alias_key, cid in alias_to_cluster.items():
            if mk.startswith(alias_key) or alias_key.startswith(mk):
                matches.append((len(alias_key), cid))
        if matches:
            matches = sorted(matches, reverse=True)
            return matches[0][1]

    return None


def load_cluster_cwms(h5_path):
    cluster_cwms = {}
    with h5py.File(h5_path, "r") as f:
        for cid_str in f["pos_patterns"].keys():
            cluster_cwms[int(cid_str)] = f["pos_patterns"][cid_str]["contrib_scores"][:]
    return cluster_cwms


def cwm_to_df(cwm):
    arr = np.asarray(cwm)
    if arr.ndim != 2:
        raise ValueError(f"Unexpected CWM ndim: {arr.ndim}")

    if arr.shape[1] == 4:
        mat = arr
    elif arr.shape[0] == 4:
        mat = arr.T
    else:
        raise ValueError(f"Unexpected CWM shape: {arr.shape}")

    return pd.DataFrame(mat, columns=list("ACGT"))


_CLUSTER_ALIAS_LOOKUP_CACHE = None
_CLUSTER_CWMS_CACHE = None


def resolve_logo_paths_for_motifs(motif_names):
    """
    Return vector CWM objects for logomaker.

    The cluster-map lookup and H5 CWMs are cached so D0-D3 do not reload
    the same large input files four times.
    """
    global _CLUSTER_ALIAS_LOOKUP_CACHE
    global _CLUSTER_CWMS_CACHE

    if _CLUSTER_ALIAS_LOOKUP_CACHE is None:
        cluster_map_df = pd.read_csv(
            CLUSTER_MAP_PATH,
            sep="\t",
        )
        _CLUSTER_ALIAS_LOOKUP_CACHE = build_cluster_lookup(
            cluster_map_df
        )

    if _CLUSTER_CWMS_CACHE is None:
        _CLUSTER_CWMS_CACHE = load_cluster_cwms(
            MODISCO_H5_PATH
        )

    alias_to_cluster = _CLUSTER_ALIAS_LOOKUP_CACHE
    cluster_cwms = _CLUSTER_CWMS_CACHE

    motif_to_logo = {}
    motif_rows = []

    for motif in motif_names:
        cluster_id = best_cluster_for_motif(
            motif,
            alias_to_cluster,
        )

        if (
            cluster_id is not None
            and cluster_id in cluster_cwms
        ):
            try:
                cwm_df = cwm_to_df(
                    cluster_cwms[cluster_id]
                )

                motif_to_logo[motif] = {
                    "type": "cwm",
                    "cwm_df": cwm_df,
                    "cluster_id": cluster_id,
                    "label": motif,
                }

                status = "vector_cwm"

            except Exception as exc:
                print(
                    f"[warning] Failed to load CWM for "
                    f"{motif} / cluster {cluster_id}: {exc}"
                )

                motif_to_logo[motif] = {
                    "type": "text",
                    "label": motif,
                    "cluster_id": cluster_id,
                }

                status = "text_fallback_after_cwm_error"

        else:
            motif_to_logo[motif] = {
                "type": "text",
                "label": motif,
                "cluster_id": cluster_id,
            }

            status = "text_fallback_no_cluster"

        motif_rows.append({
            "motif": motif,
            "matched_cluster_id": cluster_id,
            "status": status,
        })

    motif_map_df = pd.DataFrame(motif_rows)

    return motif_to_logo, motif_map_df


```

## Markdown cell 7

## Cell 4 — plotting functions with motif columns, 90° rotated logos, and shared scales

## Code cell 8

```python

# ============================================================
# PLOTTING FUNCTIONS
# Motif heatmaps retain programs as rows and motifs as columns.
# Vector motif logos are rotated by 90 degrees above the columns.
# ============================================================

def _finite_values(dataframes):
    values = []

    for df in dataframes:
        if df is None or df.empty:
            continue

        arr = df.to_numpy(dtype=float)
        arr = arr[np.isfinite(arr)]

        if arr.size:
            values.append(arr)

    if not values:
        return np.array([], dtype=float)

    return np.concatenate(values)


def compute_shared_color_scales(day_packages):
    """
    Calculate one shared scale for each data type across D0-D3.
    """
    regulator_values = _finite_values([
        pkg["reg_effect_mat"]
        for pkg in day_packages.values()
    ])

    marker_values = _finite_values([
        pkg["marker_loading_mat"]
        for pkg in day_packages.values()
    ])

    motif_values = _finite_values([
        pkg["motif_value_mat"]
        for pkg in day_packages.values()
    ])

    regulator_vmax = (
        float(np.max(np.abs(regulator_values)))
        if regulator_values.size
        else 1.0
    )

    marker_vmax = (
        float(np.max(marker_values))
        if marker_values.size
        else 1.0
    )

    motif_vmax = (
        float(np.max(motif_values))
        if motif_values.size
        else 1.0
    )

    if GLOBAL_REGULATOR_VMAX_OVERRIDE is not None:
        regulator_vmax = float(
            GLOBAL_REGULATOR_VMAX_OVERRIDE
        )

    if GLOBAL_MARKER_VMAX_OVERRIDE is not None:
        marker_vmax = float(
            GLOBAL_MARKER_VMAX_OVERRIDE
        )

    if GLOBAL_MOTIF_VMAX_OVERRIDE is not None:
        motif_vmax = float(
            GLOBAL_MOTIF_VMAX_OVERRIDE
        )

    regulator_vmax = max(
        regulator_vmax,
        1e-12,
    )

    marker_vmax = max(
        marker_vmax,
        1e-12,
    )

    motif_vmax = max(
        motif_vmax,
        1e-12,
    )

    scales = {
        "regulator_vmax": regulator_vmax,
        "marker_vmax": marker_vmax,
        "motif_vmax": motif_vmax,
        "regulator_norm": TwoSlopeNorm(
            vmin=-regulator_vmax,
            vcenter=0.0,
            vmax=regulator_vmax,
        ),
        "marker_norm": Normalize(
            vmin=0.0,
            vmax=marker_vmax,
        ),
        "motif_norm": Normalize(
            vmin=0.0,
            vmax=motif_vmax,
        ),
    }

    print("Shared D0-D3 color scales:")
    print(
        "  regulator approximate log2FC:",
        f"{-regulator_vmax:.4g} to "
        f"{regulator_vmax:.4g}",
    )
    print(
        "  marker consensus loading:",
        f"0 to {marker_vmax:.4g}",
    )
    print(
        "  motif -log10(q-value):",
        f"0 to {motif_vmax:.4g}",
    )

    return scales


def text_color_for_value(
    value,
    cmap,
    norm,
    na_color="#000000",
):
    if pd.isna(value):
        return na_color

    rgba = cmap(norm(value))

    luminance = (
        0.2126 * rgba[0]
        + 0.7152 * rgba[1]
        + 0.0722 * rgba[2]
    )

    return (
        "white"
        if luminance < 0.48
        else "black"
    )


def draw_rect_heatmap(
    ax,
    data_df,
    cmap,
    norm,
    xlabels,
    ylabels=None,
    stars_df=None,
    show_ylabels=True,
    xlabels_on_top=False,
    xlab_rotation=60,
    italic_x=False,
    xlab_fontsize=8,
    ylab_fontsize=8,
    star_fontsize=7,
    linewidth=HEATMAP_TILE_LINEWIDTH,
    missing_color=MISSING_TILE_COLOR,
):
    """
    Draw an Illustrator-safe heatmap using opaque vector rectangles.
    """
    n_rows, n_cols = data_df.shape
    values = data_df.to_numpy(dtype=float)

    for row_i in range(n_rows):
        for col_i in range(n_cols):
            value = values[
                row_i,
                col_i,
            ]

            if pd.isna(value):
                facecolor = missing_color
            else:
                rgba = cmap(norm(value))
                facecolor = (
                    rgba[0],
                    rgba[1],
                    rgba[2],
                    1.0,
                )

            ax.add_patch(
                Rectangle(
                    (col_i, row_i),
                    1,
                    1,
                    facecolor=facecolor,
                    edgecolor="white",
                    linewidth=linewidth,
                )
            )

            if stars_df is not None:
                star = stars_df.iloc[
                    row_i,
                    col_i,
                ]

                if isinstance(star, str) and star:
                    ax.text(
                        col_i + 0.5,
                        row_i + 0.5,
                        star,
                        ha="center",
                        va="center",
                        fontsize=star_fontsize,
                        fontweight="bold",
                        color=text_color_for_value(
                            value,
                            cmap,
                            norm,
                        ),
                    )

    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_rows, 0)
    ax.set_aspect("equal")

    ax.set_xticks(
        np.arange(n_cols) + 0.5
    )

    ax.set_xticklabels(
        xlabels,
        rotation=xlab_rotation,
        ha="right",
        rotation_mode="anchor",
        fontsize=xlab_fontsize,
        fontstyle=(
            "italic"
            if italic_x
            else "normal"
        ),
    )

    if xlabels_on_top:
        ax.xaxis.tick_top()
        ax.tick_params(
            axis="x",
            labeltop=True,
            labelbottom=False,
            pad=4,
        )
    else:
        ax.tick_params(
            axis="x",
            labeltop=False,
            labelbottom=True,
            pad=4,
        )

    ax.set_yticks(
        np.arange(n_rows) + 0.5
    )

    if show_ylabels:
        ax.set_yticklabels(
            ylabels,
            fontsize=ylab_fontsize,
        )
    else:
        ax.set_yticklabels(
            [""] * n_rows
        )

    ax.tick_params(
        axis="y",
        length=0,
        pad=4,
    )

    ax.tick_params(
        axis="x",
        length=0,
    )

    for spine in ax.spines.values():
        spine.set_visible(False)


def add_horizontal_colorbar(
    fig,
    left,
    bottom,
    width,
    height,
    cmap,
    norm,
    label,
    ticksize=8,
    labelsize=9,
):
    colorbar_ax = fig.add_axes([
        left,
        bottom,
        width,
        height,
    ])

    scalar = ScalarMappable(
        norm=norm,
        cmap=cmap,
    )

    scalar.set_array([])

    colorbar = fig.colorbar(
        scalar,
        cax=colorbar_ax,
        orientation="horizontal",
    )

    colorbar.outline.set_linewidth(0.8)

    colorbar.ax.tick_params(
        labelsize=ticksize,
        length=2,
    )

    colorbar.set_label(
        label,
        fontsize=labelsize,
        labelpad=2,
    )

    return colorbar_ax, colorbar


def _draw_rotated_vector_logo(
    fig,
    ax,
    motif,
    motif_to_logo,
):
    """
    Draw a vector logomaker logo and rotate its vector patches by 90 degrees.

    The rotation is applied in display coordinates, preserving vector paths
    in SVG/PDF. Text fallbacks are rotated with matplotlib text rotation.
    """
    logo_object = motif_to_logo.get(
        motif,
        {
            "type": "text",
            "label": motif,
        },
    )

    if logo_object.get("type") == "cwm":
        try:
            logomaker.Logo(
                logo_object["cwm_df"],
                ax=ax,
                color_scheme="classic",
                fade_below=0.0,
                shade_below=0.0,
                font_name="DejaVu Sans",
            )

            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_title("")

            for spine in ax.spines.values():
                spine.set_visible(False)

            # Establish the final display-coordinate axis bounding box.
            fig.canvas.draw()

            center_x = (
                ax.bbox.x0
                + ax.bbox.x1
            ) / 2

            center_y = (
                ax.bbox.y0
                + ax.bbox.y1
            ) / 2

            rotation = (
                Affine2D()
                .rotate_deg_around(
                    center_x,
                    center_y,
                    ROTATE_MOTIF_LOGOS_DEGREES,
                )
            )

            for patch in list(ax.patches):
                patch.set_transform(
                    patch.get_transform()
                    + rotation
                )

                patch.set_clip_on(False)
                patch.set_clip_path(None)

            ax.patch.set_alpha(0.0)

            return

        except Exception as exc:
            print(
                f"[warning] Could not draw rotated vector "
                f"logo for {motif}: {exc}"
            )

    ax.text(
        0.5,
        0.5,
        str(
            logo_object.get(
                "label",
                motif,
            )
        ),
        ha="center",
        va="center",
        rotation=ROTATE_MOTIF_LOGOS_DEGREES,
        fontsize=7,
        fontweight="bold",
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def add_rotated_logo_strip(
    fig,
    motif_heatmap_ax,
    motif_names,
    motif_to_logo,
):
    """
    Place one 90-degree rotated vector logo above each motif column.

    Before rotation, each logo axis is wide and short. After rotation, the
    logo is tall and narrow, matching the motif-column geometry.
    """
    if len(motif_names) == 0:
        return []

    figure_width = fig.get_figwidth()
    figure_height = fig.get_figheight()

    heat_pos = motif_heatmap_ax.get_position()

    column_width_fraction = (
        heat_pos.width
        / len(motif_names)
    )

    column_width_in = (
        column_width_fraction
        * figure_width
    )

    # Before rotation:
    # width becomes final vertical logo height;
    # height becomes final horizontal logo width.
    logo_axis_width_fraction = (
        ROTATED_LOGO_HEIGHT_IN
        / figure_width
    )

    logo_axis_height_fraction = (
        max(
            column_width_in * 0.88,
            0.12,
        )
        / figure_height
    )

    logo_center_y = (
        heat_pos.y1
        + ROTATED_LOGO_PAD_IN / figure_height
        + ROTATED_LOGO_HEIGHT_IN
          / (2 * figure_height)
    )

    axes = []

    for col_i, motif in enumerate(motif_names):
        column_center_x = (
            heat_pos.x0
            + (col_i + 0.5)
              * column_width_fraction
        )

        logo_ax = fig.add_axes([
            column_center_x
            - logo_axis_width_fraction / 2,
            logo_center_y
            - logo_axis_height_fraction / 2,
            logo_axis_width_fraction,
            logo_axis_height_fraction,
        ])

        logo_ax.set_zorder(5)

        _draw_rotated_vector_logo(
            fig,
            logo_ax,
            motif,
            motif_to_logo,
        )

        axes.append(logo_ax)

    return axes


def _prepare_day_matrices(day_package):
    program_meta = day_package["meta"]

    program_numbers = [
        int(value)
        for value in program_meta["program_num"]
    ]

    program_labels = (
        program_meta["row_label"]
        .tolist()
    )

    regulator_union = (
        day_package["reg_union"]
    )

    marker_union = (
        day_package["marker_union"]
    )

    motif_union = (
        day_package["motif_union"]
    )

    regulator_matrix = (
        day_package["reg_effect_mat"]
        .loc[
            program_numbers,
            regulator_union,
        ]
    )

    regulator_stars = (
        day_package["reg_star_mat"]
        .loc[
            program_numbers,
            regulator_union,
        ]
    )

    marker_matrix = (
        day_package["marker_loading_mat"]
        .loc[
            program_numbers,
            marker_union,
        ]
    )

    # Keep original orientation:
    # programs are rows and motifs are columns.
    motif_matrix = (
        day_package["motif_value_mat"]
        .loc[
            program_numbers,
            motif_union,
        ]
    )

    motif_stars = (
        day_package["motif_star_mat"]
        .loc[
            program_numbers,
            motif_union,
        ]
    )

    return {
        "program_numbers": program_numbers,
        "program_labels": program_labels,
        "regulator_union": regulator_union,
        "marker_union": marker_union,
        "motif_union": motif_union,
        "regulator_matrix": regulator_matrix,
        "regulator_stars": regulator_stars,
        "marker_matrix": marker_matrix,
        "motif_matrix": motif_matrix,
        "motif_stars": motif_stars,
    }


def make_per_day_panel(
    day,
    day_package,
    shared_scales,
):
    """
    Generate one day-specific panel with the shared D0-D3 scales.
    """
    prepared = _prepare_day_matrices(
        day_package
    )

    n_programs = len(
        prepared["program_numbers"]
    )

    n_regulators = len(
        prepared["regulator_union"]
    )

    n_markers = len(
        prepared["marker_union"]
    )

    n_motifs = len(
        prepared["motif_union"]
    )

    regulator_width_in = max(
        2.2,
        n_regulators * PER_DAY_CELL_W_IN,
    )

    marker_width_in = max(
        2.2,
        n_markers * PER_DAY_CELL_W_IN,
    )

    motif_width_in = max(
        2.2,
        n_motifs * PER_DAY_CELL_W_IN,
    )

    heatmap_height_in = max(
        PER_DAY_CELL_H_IN,
        n_programs * PER_DAY_CELL_H_IN,
    )

    logo_region_height_in = (
        ROTATED_LOGO_HEIGHT_IN
        + ROTATED_LOGO_PAD_IN
    )

    header_height_in = 1.25

    figure_width = (
        PER_DAY_LEFT_MARGIN_IN
        + regulator_width_in
        + PER_DAY_BLOCK_GAP_IN
        + marker_width_in
        + PER_DAY_BLOCK_GAP_IN
        + motif_width_in
        + PER_DAY_RIGHT_MARGIN_IN
    )

    figure_height = (
        PER_DAY_BOTTOM_MARGIN_IN
        + heatmap_height_in
        + logo_region_height_in
        + header_height_in
        + PER_DAY_TOP_MARGIN_IN
    )

    fig = plt.figure(
        figsize=(
            figure_width,
            figure_height,
        ),
        facecolor="white",
    )

    def fx(value):
        return value / figure_width

    def fy(value):
        return value / figure_height

    regulator_left = (
        PER_DAY_LEFT_MARGIN_IN
    )

    marker_left = (
        regulator_left
        + regulator_width_in
        + PER_DAY_BLOCK_GAP_IN
    )

    motif_left = (
        marker_left
        + marker_width_in
        + PER_DAY_BLOCK_GAP_IN
    )

    heatmap_bottom = (
        PER_DAY_BOTTOM_MARGIN_IN
    )

    regulator_ax = fig.add_axes([
        fx(regulator_left),
        fy(heatmap_bottom),
        fx(regulator_width_in),
        fy(heatmap_height_in),
    ])

    marker_ax = fig.add_axes([
        fx(marker_left),
        fy(heatmap_bottom),
        fx(marker_width_in),
        fy(heatmap_height_in),
    ])

    motif_ax = fig.add_axes([
        fx(motif_left),
        fy(heatmap_bottom),
        fx(motif_width_in),
        fy(heatmap_height_in),
    ])

    y_label_fontsize = max(
        7,
        min(
            12,
            95 / max(n_programs, 1),
        ),
    )

    regulator_x_fontsize = max(
        5,
        min(
            9,
            100 / max(n_regulators, 1),
        ),
    )

    marker_x_fontsize = max(
        5,
        min(
            9,
            100 / max(n_markers, 1),
        ),
    )

    motif_x_fontsize = max(
        5,
        min(
            9,
            100 / max(n_motifs, 1),
        ),
    )

    star_fontsize = max(
        4,
        min(
            8,
            80 / max(n_programs, 1),
        ),
    )

    draw_rect_heatmap(
        regulator_ax,
        prepared["regulator_matrix"],
        REGULATOR_CMAP,
        shared_scales["regulator_norm"],
        xlabels=prepared["regulator_union"],
        ylabels=prepared["program_labels"],
        stars_df=prepared["regulator_stars"],
        show_ylabels=True,
        xlabels_on_top=False,
        xlab_rotation=58,
        italic_x=False,
        xlab_fontsize=regulator_x_fontsize,
        ylab_fontsize=y_label_fontsize,
        star_fontsize=star_fontsize,
    )

    draw_rect_heatmap(
        marker_ax,
        prepared["marker_matrix"],
        MARKER_CMAP,
        shared_scales["marker_norm"],
        xlabels=prepared["marker_union"],
        ylabels=prepared["program_labels"],
        stars_df=None,
        show_ylabels=False,
        xlabels_on_top=False,
        xlab_rotation=58,
        italic_x=True,
        xlab_fontsize=marker_x_fontsize,
        ylab_fontsize=y_label_fontsize,
        star_fontsize=star_fontsize,
    )

    draw_rect_heatmap(
        motif_ax,
        prepared["motif_matrix"],
        MOTIF_CMAP,
        shared_scales["motif_norm"],
        xlabels=prepared["motif_union"],
        ylabels=prepared["program_labels"],
        stars_df=prepared["motif_stars"],
        show_ylabels=False,
        xlabels_on_top=False,
        xlab_rotation=58,
        italic_x=False,
        xlab_fontsize=motif_x_fontsize,
        ylab_fontsize=y_label_fontsize,
        star_fontsize=star_fontsize,
    )

    add_rotated_logo_strip(
        fig,
        motif_ax,
        prepared["motif_union"],
        day_package["motif_to_logo"],
    )

    heatmap_top = (
        heatmap_bottom
        + heatmap_height_in
    )

    colorbar_y = (
        heatmap_top
        + logo_region_height_in
        + 0.18
    )

    colorbar_height_in = 0.14

    add_horizontal_colorbar(
        fig,
        fx(regulator_left),
        fy(colorbar_y),
        fx(regulator_width_in),
        fy(colorbar_height_in),
        REGULATOR_CMAP,
        shared_scales["regulator_norm"],
        "approx. log$_2$FC",
    )

    add_horizontal_colorbar(
        fig,
        fx(marker_left),
        fy(colorbar_y),
        fx(marker_width_in),
        fy(colorbar_height_in),
        MARKER_CMAP,
        shared_scales["marker_norm"],
        "consensus cNMF loading",
    )

    add_horizontal_colorbar(
        fig,
        fx(motif_left),
        fy(colorbar_y),
        fx(motif_width_in),
        fy(colorbar_height_in),
        MOTIF_CMAP,
        shared_scales["motif_norm"],
        "$-\\log_{10}$(q-value)",
    )

    title_y = (
        colorbar_y
        + colorbar_height_in
        + 0.28
    )

    fig.text(
        fx(
            regulator_left
            + regulator_width_in / 2
        ),
        fy(title_y),
        "Regulator effects",
        ha="center",
        va="bottom",
        fontsize=15,
        fontweight="bold",
    )

    fig.text(
        fx(
            marker_left
            + marker_width_in / 2
        ),
        fy(title_y),
        "Marker-gene loadings",
        ha="center",
        va="bottom",
        fontsize=15,
        fontweight="bold",
    )

    fig.text(
        fx(
            motif_left
            + motif_width_in / 2
        ),
        fy(title_y),
        "Motif enrichment",
        ha="center",
        va="bottom",
        fontsize=15,
        fontweight="bold",
    )

    day_title = {
        "D0": "Day 0",
        "D1": "Day 1",
        "D2": "Day 2",
        "D3": "Day 3",
    }.get(day, str(day))

    fig.text(
        fx(0.25),
        fy(figure_height - 0.25),
        day_title,
        ha="left",
        va="top",
        fontsize=23,
        fontweight="bold",
    )

    fig.text(
        fx(regulator_left),
        fy(0.30),
        "* < 0.05;  ** < 0.01;  *** < 0.001. "
        "Top three motifs per program are shown regardless of significance.",
        fontsize=8.5,
        color="#444444",
        ha="left",
        va="bottom",
    )

    output_stem = (
        OUTDIR
        / (
            f"Program_annotation_panel_{day}_"
            f"motif_columns_rotated_logos"
        )
    )

    for extension in [
        "svg",
        "pdf",
        "png",
    ]:
        save_kwargs = {
            "bbox_inches": "tight",
            "facecolor": "white",
            "transparent": False,
        }

        if extension == "png":
            save_kwargs["dpi"] = 600

        fig.savefig(
            output_stem.with_suffix(
                f".{extension}"
            ),
            **save_kwargs,
        )

    plt.show()
    plt.close(fig)

    return output_stem


def make_combined_all_days_figure(
    day_packages,
    days,
    shared_scales,
):
    """
    Stack D0-D3 vertically.

    Motif enrichment retains programs as rows and motifs as columns.
    Every motif logo is rotated 90 degrees above its corresponding column.
    Exactly one shared colorbar is used for each data type.
    """
    valid_days = [
        day
        for day in days
        if day in day_packages
    ]

    if not valid_days:
        raise ValueError(
            "No valid day packages were available."
        )

    prepared_by_day = {
        day: _prepare_day_matrices(
            day_packages[day]
        )
        for day in valid_days
    }

    max_regulator_columns = max(
        len(prepared["regulator_union"])
        for prepared in prepared_by_day.values()
    )

    max_marker_columns = max(
        len(prepared["marker_union"])
        for prepared in prepared_by_day.values()
    )

    max_motif_columns = max(
        len(prepared["motif_union"])
        for prepared in prepared_by_day.values()
    )

    regulator_block_width_in = max(
        2.4,
        max_regulator_columns
        * COMBINED_CELL_W_IN,
    )

    marker_block_width_in = max(
        2.4,
        max_marker_columns
        * COMBINED_CELL_W_IN,
    )

    motif_block_width_in = max(
        2.4,
        max_motif_columns
        * COMBINED_CELL_W_IN,
    )

    figure_width = (
        COMBINED_LEFT_MARGIN_IN
        + regulator_block_width_in
        + COMBINED_BLOCK_GAP_IN
        + marker_block_width_in
        + COMBINED_BLOCK_GAP_IN
        + motif_block_width_in
        + COMBINED_RIGHT_MARGIN_IN
    )

    day_layout = {}
    total_day_height = 0.0

    for day in valid_days:
        prepared = prepared_by_day[day]

        n_programs = len(
            prepared["program_numbers"]
        )

        heatmap_height_in = max(
            COMBINED_CELL_H_IN,
            n_programs
            * COMBINED_CELL_H_IN,
        )

        band_height_in = (
            COMBINED_DAY_LABEL_SPACE_IN
            + ROTATED_LOGO_HEIGHT_IN
            + ROTATED_LOGO_PAD_IN
            + heatmap_height_in
            + COMBINED_XLABEL_SPACE_IN
        )

        day_layout[day] = {
            "n_programs": n_programs,
            "heatmap_height_in": heatmap_height_in,
            "band_height_in": band_height_in,
        }

        total_day_height += band_height_in

    total_day_height += (
        COMBINED_DAY_GAP_IN
        * (len(valid_days) - 1)
    )

    bottom_margin_in = 0.75

    figure_height = (
        bottom_margin_in
        + total_day_height
        + COMBINED_HEADER_IN
    )

    fig = plt.figure(
        figsize=(
            figure_width,
            figure_height,
        ),
        facecolor="white",
    )

    def fx(value):
        return value / figure_width

    def fy(value):
        return value / figure_height

    regulator_left = (
        COMBINED_LEFT_MARGIN_IN
    )

    marker_left = (
        regulator_left
        + regulator_block_width_in
        + COMBINED_BLOCK_GAP_IN
    )

    motif_left = (
        marker_left
        + marker_block_width_in
        + COMBINED_BLOCK_GAP_IN
    )

    # ========================================================
    # ONE GLOBAL TITLE AND ONE COLORBAR PER DATA TYPE
    # ========================================================

    fig.text(
        fx(0.25),
        fy(figure_height - 0.18),
        "Program annotations across differentiation",
        ha="left",
        va="top",
        fontsize=24,
        fontweight="bold",
    )

    block_title_y = (
        figure_height
        - 0.62
    )

    colorbar_y = (
        figure_height
        - 1.10
    )

    colorbar_height_in = 0.14

    fig.text(
        fx(
            regulator_left
            + regulator_block_width_in / 2
        ),
        fy(block_title_y),
        "Regulator effects",
        ha="center",
        va="bottom",
        fontsize=16,
        fontweight="bold",
    )

    fig.text(
        fx(
            marker_left
            + marker_block_width_in / 2
        ),
        fy(block_title_y),
        "Marker-gene loadings",
        ha="center",
        va="bottom",
        fontsize=16,
        fontweight="bold",
    )

    fig.text(
        fx(
            motif_left
            + motif_block_width_in / 2
        ),
        fy(block_title_y),
        "Motif enrichment",
        ha="center",
        va="bottom",
        fontsize=16,
        fontweight="bold",
    )

    add_horizontal_colorbar(
        fig,
        fx(regulator_left),
        fy(colorbar_y),
        fx(regulator_block_width_in),
        fy(colorbar_height_in),
        REGULATOR_CMAP,
        shared_scales["regulator_norm"],
        "approx. log$_2$FC",
    )

    add_horizontal_colorbar(
        fig,
        fx(marker_left),
        fy(colorbar_y),
        fx(marker_block_width_in),
        fy(colorbar_height_in),
        MARKER_CMAP,
        shared_scales["marker_norm"],
        "consensus cNMF loading",
    )

    add_horizontal_colorbar(
        fig,
        fx(motif_left),
        fy(colorbar_y),
        fx(motif_block_width_in),
        fy(colorbar_height_in),
        MOTIF_CMAP,
        shared_scales["motif_norm"],
        "$-\\log_{10}$(q-value)",
    )

    fig.text(
        fx(regulator_left),
        fy(figure_height - 1.46),
        "* < 0.05;  ** < 0.01;  *** < 0.001. "
        "Top three motifs per program are shown regardless of significance.",
        fontsize=9,
        color="#444444",
        ha="left",
        va="top",
    )

    # ========================================================
    # DAY BANDS
    # ========================================================

    current_top = (
        figure_height
        - COMBINED_HEADER_IN
    )

    for day_index, day in enumerate(
        valid_days
    ):
        prepared = prepared_by_day[day]
        layout = day_layout[day]

        n_programs = layout["n_programs"]

        n_regulators = len(
            prepared["regulator_union"]
        )

        n_markers = len(
            prepared["marker_union"]
        )

        n_motifs = len(
            prepared["motif_union"]
        )

        heatmap_height_in = (
            layout["heatmap_height_in"]
        )

        day_title_y = (
            current_top
            - 0.08
        )

        heatmap_top = (
            current_top
            - COMBINED_DAY_LABEL_SPACE_IN
            - ROTATED_LOGO_HEIGHT_IN
            - ROTATED_LOGO_PAD_IN
        )

        heatmap_bottom = (
            heatmap_top
            - heatmap_height_in
        )

        regulator_actual_width_in = max(
            COMBINED_CELL_W_IN,
            n_regulators
            * COMBINED_CELL_W_IN,
        )

        marker_actual_width_in = max(
            COMBINED_CELL_W_IN,
            n_markers
            * COMBINED_CELL_W_IN,
        )

        motif_actual_width_in = max(
            COMBINED_CELL_W_IN,
            n_motifs
            * COMBINED_CELL_W_IN,
        )

        regulator_ax = fig.add_axes([
            fx(regulator_left),
            fy(heatmap_bottom),
            fx(regulator_actual_width_in),
            fy(heatmap_height_in),
        ])

        marker_ax = fig.add_axes([
            fx(marker_left),
            fy(heatmap_bottom),
            fx(marker_actual_width_in),
            fy(heatmap_height_in),
        ])

        motif_ax = fig.add_axes([
            fx(motif_left),
            fy(heatmap_bottom),
            fx(motif_actual_width_in),
            fy(heatmap_height_in),
        ])

        y_label_fontsize = max(
            6.0,
            min(
                10.0,
                72 / max(n_programs, 1),
            ),
        )

        regulator_x_fontsize = max(
            4.5,
            min(
                8.0,
                82 / max(n_regulators, 1),
            ),
        )

        marker_x_fontsize = max(
            4.5,
            min(
                8.0,
                82 / max(n_markers, 1),
            ),
        )

        motif_x_fontsize = max(
            4.5,
            min(
                8.0,
                82 / max(n_motifs, 1),
            ),
        )

        star_fontsize = max(
            4.0,
            min(
                7.0,
                55 / max(n_programs, 1),
            ),
        )

        draw_rect_heatmap(
            regulator_ax,
            prepared["regulator_matrix"],
            REGULATOR_CMAP,
            shared_scales["regulator_norm"],
            xlabels=prepared["regulator_union"],
            ylabels=prepared["program_labels"],
            stars_df=prepared["regulator_stars"],
            show_ylabels=True,
            xlabels_on_top=False,
            xlab_rotation=58,
            italic_x=False,
            xlab_fontsize=regulator_x_fontsize,
            ylab_fontsize=y_label_fontsize,
            star_fontsize=star_fontsize,
        )

        draw_rect_heatmap(
            marker_ax,
            prepared["marker_matrix"],
            MARKER_CMAP,
            shared_scales["marker_norm"],
            xlabels=prepared["marker_union"],
            ylabels=prepared["program_labels"],
            stars_df=None,
            show_ylabels=False,
            xlabels_on_top=False,
            xlab_rotation=58,
            italic_x=True,
            xlab_fontsize=marker_x_fontsize,
            ylab_fontsize=y_label_fontsize,
            star_fontsize=star_fontsize,
        )

        draw_rect_heatmap(
            motif_ax,
            prepared["motif_matrix"],
            MOTIF_CMAP,
            shared_scales["motif_norm"],
            xlabels=prepared["motif_union"],
            ylabels=prepared["program_labels"],
            stars_df=prepared["motif_stars"],
            show_ylabels=False,
            xlabels_on_top=False,
            xlab_rotation=58,
            italic_x=False,
            xlab_fontsize=motif_x_fontsize,
            ylab_fontsize=y_label_fontsize,
            star_fontsize=star_fontsize,
        )

        add_rotated_logo_strip(
            fig,
            motif_ax,
            prepared["motif_union"],
            day_packages[day]["motif_to_logo"],
        )

        day_title = {
            "D0": "Day 0",
            "D1": "Day 1",
            "D2": "Day 2",
            "D3": "Day 3",
        }.get(
            day,
            str(day),
        )

        fig.text(
            fx(0.24),
            fy(day_title_y),
            day_title,
            ha="left",
            va="top",
            fontsize=18,
            fontweight="bold",
        )

        band_bottom = (
            heatmap_bottom
            - COMBINED_XLABEL_SPACE_IN
        )

        if day_index < len(valid_days) - 1:
            separator_y = (
                band_bottom
                - COMBINED_DAY_GAP_IN / 2
            )

            separator_ax = fig.add_axes([
                fx(0.20),
                fy(separator_y),
                fx(figure_width - 0.40),
                0.001,
            ])

            separator_ax.axhline(
                0,
                color="#D6D6D6",
                linewidth=0.8,
            )

            separator_ax.axis("off")

        current_top = (
            band_bottom
            - COMBINED_DAY_GAP_IN
        )

    output_stem = (
        OUTDIR
        / (
            "Program_annotation_D0_D1_D2_D3_"
            "combined_motif_columns_rotated_logos"
        )
    )

    for extension in [
        "svg",
        "pdf",
        "png",
    ]:
        save_kwargs = {
            "bbox_inches": "tight",
            "facecolor": "white",
            "transparent": False,
        }

        if extension == "png":
            save_kwargs["dpi"] = 600

        fig.savefig(
            output_stem.with_suffix(
                f".{extension}"
            ),
            **save_kwargs,
        )

    plt.show()
    plt.close(fig)

    return output_stem

```

## Markdown cell 9

## Cell 5 — load data, select features, and save source tables

## Code cell 10

```python

# ============================================================
# LOAD DATA, SELECT FEATURES, SAVE SOURCE TABLES
# ============================================================

program_meta = load_program_metadata()
pert = load_perturbation_associations()
enh = load_enhancer_enrichment()
consensus = read_consensus_matrix(CONSENSUS_PATH)

print("Program metadata:", program_meta.shape)
print("Perturb associations:", pert.shape)
print("Enhancer enrichment:", enh.shape)
print("Consensus matrix:", consensus.shape)

# choose program offset for consensus matrix
consensus_offset = choose_program_offset(consensus, program_meta)

program_meta.to_csv(SOURCE_DIR / "program_metadata_merged.csv", index=False)

day_packages = {}

for day in DAYS_TO_PLOT:
    meta_day = (
        program_meta.loc[program_meta["day"] == day]
        .sort_values("program_num")
        .reset_index(drop=True)
    )

    if meta_day.empty:
        print(f"[warning] no programs found for {day}")
        continue

    # regulators
    reg_by_prog, reg_union, reg_selected_df = select_top_regulators_for_programs(meta_day, pert)
    reg_effect_mat, reg_fdr_mat, reg_star_mat = build_regulator_matrices(meta_day, reg_union, pert)

    # markers
    marker_by_prog, marker_union, marker_selected_df = select_marker_genes_for_programs(
        program_meta_day=meta_day,
        consensus_df=consensus,
        consensus_offset=consensus_offset,
    )
    marker_loading_mat = build_marker_loading_matrix(meta_day, marker_union, consensus, consensus_offset)

    # motifs
    motif_by_prog, motif_union, motif_selected_df = select_top_motifs_for_programs(meta_day, enh)
    motif_value_mat, motif_qval_mat, motif_star_mat = build_motif_matrices(meta_day, motif_union, enh)

    # logos
    motif_to_logo, motif_logo_map_df = resolve_logo_paths_for_motifs(motif_union)

    day_packages[day] = {
        "meta": meta_day,
        "reg_union": reg_union,
        "reg_selected_df": reg_selected_df,
        "reg_effect_mat": reg_effect_mat,
        "reg_fdr_mat": reg_fdr_mat,
        "reg_star_mat": reg_star_mat,
        "marker_union": marker_union,
        "marker_selected_df": marker_selected_df,
        "marker_loading_mat": marker_loading_mat,
        "motif_union": motif_union,
        "motif_selected_df": motif_selected_df,
        "motif_value_mat": motif_value_mat,
        "motif_qval_mat": motif_qval_mat,
        "motif_star_mat": motif_star_mat,
        "motif_to_logo": motif_to_logo,
        "motif_logo_map_df": motif_logo_map_df,
    }

    # Save source tables
    prefix = SOURCE_DIR / f"{day}"

    meta_day.to_csv(f"{prefix}_programs.csv", index=False)

    reg_selected_df.to_csv(f"{prefix}_selected_regulators.csv", index=False)
    reg_effect_mat.to_csv(f"{prefix}_regulator_effect_matrix.csv", index_label="program_num")
    reg_fdr_mat.to_csv(f"{prefix}_regulator_fdr_matrix.csv", index_label="program_num")
    reg_star_mat.to_csv(f"{prefix}_regulator_star_matrix.csv", index_label="program_num")

    marker_selected_df.to_csv(f"{prefix}_selected_marker_genes.csv", index=False)
    marker_loading_mat.to_csv(f"{prefix}_marker_loading_matrix.csv", index_label="program_num")

    motif_selected_df.to_csv(f"{prefix}_selected_motifs.csv", index=False)
    motif_value_mat.to_csv(f"{prefix}_motif_enrichment_matrix.csv", index_label="program_num")
    motif_qval_mat.to_csv(f"{prefix}_motif_qval_matrix.csv", index_label="program_num")
    motif_star_mat.to_csv(f"{prefix}_motif_star_matrix.csv", index_label="program_num")
    motif_logo_map_df.to_csv(f"{prefix}_motif_logo_map.csv", index=False)

    print(f"\n=== {day} ===")
    print("Programs:", meta_day["program_num"].tolist())
    print("Selected regulators:", len(reg_union))
    print("Selected marker genes:", len(marker_union))
    print("Selected motifs:", len(motif_union))


# ============================================================
# SHARED D0-D3 COLOR SCALES
# ============================================================

shared_scales = compute_shared_color_scales(
    day_packages
)

shared_scale_values = {
    "regulator_vmin": -shared_scales["regulator_vmax"],
    "regulator_vmax": shared_scales["regulator_vmax"],
    "marker_vmin": 0.0,
    "marker_vmax": shared_scales["marker_vmax"],
    "motif_vmin": 0.0,
    "motif_vmax": shared_scales["motif_vmax"],
}

with open(
    SOURCE_DIR
    / "shared_color_scales_D0_D3.json",
    "w",
) as handle:
    json.dump(
        shared_scale_values,
        handle,
        indent=2,
    )


# ============================================================
# SAVE CONFIGURATION
# ============================================================

config_used = {
    "CONSENSUS_PATH": str(CONSENSUS_PATH),
    "PROGRAM_ANNOT_PATH": str(PROGRAM_ANNOT_PATH),
    "CNMF_WORKBOOK_PATH": str(CNMF_WORKBOOK_PATH),
    "ENHANCER_ENRICHMENT_PATH": str(
        ENHANCER_ENRICHMENT_PATH
    ),
    "CLUSTER_MAP_PATH": str(CLUSTER_MAP_PATH),
    "MODISCO_H5_PATH": str(MODISCO_H5_PATH),
    "PROGRAM_OFFSET": PROGRAM_OFFSET,
    "consensus_offset_chosen": consensus_offset,
    "N_TOP_REGULATORS": N_TOP_REGULATORS,
    "N_TOP_MOTIFS": N_TOP_MOTIFS,
    "D2_N_ADDITIONAL_TOP_LOADING_GENES": D2_N_ADDITIONAL_TOP_LOADING_GENES,
    "D3_N_CURATED_MARKERS_TO_KEEP": D3_N_CURATED_MARKERS_TO_KEEP,
    "D2_CURATED_GENES_FIRST": D2_CURATED_GENES_FIRST,
    "REGULATOR_FDR_THRESHOLD": (
        REGULATOR_FDR_THRESHOLD
    ),
    "MOTIF_QVAL_THRESHOLD": (
        MOTIF_QVAL_THRESHOLD
    ),
    "VECTOR_LOGOS": VECTOR_LOGOS,
    "ROTATE_MOTIF_LOGOS_DEGREES": (
        ROTATE_MOTIF_LOGOS_DEGREES
    ),
    "MAKE_PER_DAY_PANELS": (
        MAKE_PER_DAY_PANELS
    ),
    "MAKE_COMBINED_ALL_DAYS_FIGURE": (
        MAKE_COMBINED_ALL_DAYS_FIGURE
    ),
    "shared_color_scales": (
        shared_scale_values
    ),
}

with open(
    SOURCE_DIR
    / "figure_config_used.json",
    "w",
) as handle:
    json.dump(
        config_used,
        handle,
        indent=2,
    )

print(
    "\nSource tables written to:",
    SOURCE_DIR,
)

print(
    "Shared color scales written to:",
    SOURCE_DIR
    / "shared_color_scales_D0_D3.json",
)

```

### Text output (stream)

```text
Program metadata: (17, 42)
Perturb associations: (59800, 13)
Enhancer enrichment: (13906, 16)
Consensus matrix: (50, 17538)
Consensus program-offset scores: {0: 5.070201602456141e-05, 1: 0.0014924993410526314, -1: 5.3358765333333336e-05}
Chosen consensus program offset: 1
Cluster-map ID column: cluster_id
Cluster-map motif columns: ['days_present', 'posneg', 'annotation_name0', 'annotation_names', 'TFs', 'tf_motif_id', 'expressed_TFs_d0', 'expressed_TFs_d1', 'expressed_TFs_d2', 'expressed_TFs_d3']

=== D0 ===
Programs: [7, 22, 29]
Selected regulators: 9
Selected marker genes: 9
Selected motifs: 8

=== D1 ===
Programs: [13, 28, 30, 48]
Selected regulators: 10
Selected marker genes: 8
Selected motifs: 7

=== D2 ===
Programs: [14, 15, 17]
Selected regulators: 9
Selected marker genes: 17
Selected motifs: 5

=== D3 ===
Programs: [2, 3, 8, 11, 27, 33, 34]
Selected regulators: 21
Selected marker genes: 20
Selected motifs: 8
Shared D0-D3 color scales:
  regulator approximate log2FC: -9.922 to 9.922
  marker consensus loading: 0 to 0.003387
  motif -log10(q-value): 0 to 56.31

Source tables written to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/_program_annotation_source
Shared color scales written to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/_program_annotation_source/shared_color_scales_D0_D3.json

```

## Markdown cell 11

## Cell 6 — render per-day panels and the combined D0–D3 figure

## Code cell 12

```python

# ============================================================
# RENDER PER-DAY PANELS AND COMBINED D0-D3 FIGURE
# ============================================================

if MAKE_PER_DAY_PANELS:
    print(
        "Rendering per-day panels with programs as rows, "
        "motifs as columns, and 90-degree rotated logos..."
    )

    for day in DAYS_TO_PLOT:
        if day not in day_packages:
            print(
                f"[skip] {day}: no package"
            )
            continue

        output_stem = make_per_day_panel(
            day=day,
            day_package=day_packages[day],
            shared_scales=shared_scales,
        )

        print(f"Saved {day}:")
        print(
            " ",
            output_stem.with_suffix(".svg"),
        )
        print(
            " ",
            output_stem.with_suffix(".pdf"),
        )
        print(
            " ",
            output_stem.with_suffix(".png"),
        )


if MAKE_COMBINED_ALL_DAYS_FIGURE:
    print(
        "\nRendering the combined D0-D3 figure "
        "with one shared colorbar per data type..."
    )

    combined_output_stem = (
        make_combined_all_days_figure(
            day_packages=day_packages,
            days=DAYS_TO_PLOT,
            shared_scales=shared_scales,
        )
    )

    print("Saved combined figure:")
    print(
        " ",
        combined_output_stem.with_suffix(
            ".svg"
        ),
    )
    print(
        " ",
        combined_output_stem.with_suffix(
            ".pdf"
        ),
    )
    print(
        " ",
        combined_output_stem.with_suffix(
            ".png"
        ),
    )

```

### Text output (stream)

```text
Rendering per-day panels with programs as rows, motifs as columns, and 90-degree rotated logos...

```

### Text output (text/plain)

```text
<Figure size 1424x654 with 14 Axes>
```

### Text output (stream)

```text
Saved D0:
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D0_motif_columns_rotated_logos.svg
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D0_motif_columns_rotated_logos.pdf
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D0_motif_columns_rotated_logos.png

```

### Text output (text/plain)

```text
<Figure size 1390x688 with 13 Axes>
```

### Text output (stream)

```text
Saved D1:
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D1_motif_columns_rotated_logos.svg
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D1_motif_columns_rotated_logos.pdf
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D1_motif_columns_rotated_logos.png

```

### Text output (text/plain)

```text
<Figure size 1644x654 with 11 Axes>
```

### Text output (stream)

```text
Saved D2:
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D2_motif_columns_rotated_logos.svg
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D2_motif_columns_rotated_logos.pdf
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D2_motif_columns_rotated_logos.png

```

### Text output (text/plain)

```text
<Figure size 2206x790 with 14 Axes>
```

### Text output (stream)

```text
Saved D3:
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D3_motif_columns_rotated_logos.svg
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D3_motif_columns_rotated_logos.pdf
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_panel_D3_motif_columns_rotated_logos.png

Rendering the combined D0-D3 figure with one shared colorbar per data type...

```

### Text output (text/plain)

```text
<Figure size 2025x2156 with 46 Axes>
```

### Text output (stream)

```text
Saved combined figure:
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_D0_D1_D2_D3_combined_motif_columns_rotated_logos.svg
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_D0_D1_D2_D3_combined_motif_columns_rotated_logos.pdf
  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/per_day_umaps/Program_annotation_D0_D1_D2_D3_combined_motif_columns_rotated_logos.png

```

