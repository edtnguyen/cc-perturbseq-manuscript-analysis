# ============================================================
# Panel D replacement — D3 EC vs non-EC by individual Hippo target
# ============================================================
# This cell assumes an AnnData object named `adata` is already loaded in the notebook.
# It replaces the old D2/D3 three-state pooled-Hippo Panel D section.
#
# What it does:
#   1) D3 only
#   2) Two states only: EC vs non-EC
#   3) Individual Hippo target bars: NF2, LATS2, SAV1, TAOK1, AMOTL2
#   4) Keeps the two control bars: NT excl. Hippo and Strict NT
#   5) Saves Illustrator-editable SVG/PDF and source CSVs

import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import font_manager
from scipy import sparse
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from plotnine import (
    ggplot, aes,
    geom_rect, geom_text,
    scale_fill_manual, scale_x_continuous, scale_y_continuous,
    coord_cartesian, labs,
    theme_classic, theme,
    element_text, element_blank, element_rect
)

warnings.filterwarnings("ignore")

# ============================================================
# Output
# ============================================================

OUTDIR = Path("panelD_D3_EC_vs_nonEC_individual_Hippo")
OUTDIR.mkdir(parents=True, exist_ok=True)

COMPOSITION_SOURCE_OUT = OUTDIR / "panelD_D3_EC_vs_nonEC_individual_Hippo_source.csv"
STATS_SOURCE_OUT = OUTDIR / "panelD_D3_EC_vs_nonEC_individual_Hippo_stats.csv"

# ============================================================
# AnnData fields
# ============================================================

DAY_COL = "day"
CELLTYPE_COL = "cell_type"

GUIDE_NAMES_KEY = "guide_names"
GUIDE_ASSIGNMENT_KEY = "guide_assignment"

# ============================================================
# Guide metadata
# ============================================================

GUIDE_METADATA_TSV = "300genes_guide_metadata_v43.tsv"
GUIDE_ID_COL = "guide_id"
# Robust to both old/new metadata column names.
TARGET_COL_CANDIDATES = ["putative_target_genes", "intended_target_name"]

# ============================================================
# Groups
# ============================================================

HIPPO_TARGETS = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]
CONTROL_TARGET = "non-targeting_control"

# If True, a target bar includes only cells carrying that Hippo target and no other Hippo target.
# If False, each bar is all carriers of that target. This matches the usual target-carrier definition.
EXCLUSIVE_HIPPO_TARGETS = False

CONTROL_PATTERNS = [
    r"non[-_ ]?target",
    r"non[-_ ]?targeting",
    r"\bNTC\b",
    r"negative",
    r"neg[-_ ]?ctrl",
    r"control",
    r"scram",
    r"scrambled",
    r"safe",
    r"aavs1",
]

# ============================================================
# D3 EC state definition
# ============================================================

D3_EC_PATTERNS = [
    r"\bartEC\b",
    # uncomment/add if your EC labels vary:
    # r"arterial.*EC",
    # r"endothelial",
]

D3_STATE_ORDER = ["artEC", "non-EC"]
D3_STATE_COLORS = {
    # artEC in red; non-EC in neutral grey.
    "artEC": "#C73D4C",
    "non-EC": "#B8B8B8",
}

GROUP_ORDER = HIPPO_TARGETS + ["NT excl. Hippo", "Strict NT"]

# ============================================================
# Illustrator-editable vector settings
# ============================================================

logging_level_font = {f.name for f in font_manager.fontManager.ttflist}
FONT_CANDIDATES = ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]
AI_FONT = next((f for f in FONT_CANDIDATES if f in logging_level_font), "DejaVu Sans")

mpl.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": [AI_FONT, "DejaVu Sans"],
    "savefig.dpi": 600,
    "savefig.transparent": False,
    "patch.antialiased": False,
})

print("Using font:", AI_FONT)

# ============================================================
# Helper functions
# ============================================================

def clean_day(x):
    s = str(x).lower().strip().replace("day", "d").replace(" ", "")
    if s in ["0", "0.0"]:
        return "D0"
    if s in ["1", "1.0"]:
        return "D1"
    if s in ["2", "2.0"]:
        return "D2"
    if s in ["3", "3.0"]:
        return "D3"
    return s.upper()


def matches_any(label, patterns):
    s = str(label)
    return any(re.search(p, s, flags=re.IGNORECASE) for p in patterns)


def looks_control(x):
    s = str(x).lower()
    return any(re.search(p, s, flags=re.IGNORECASE) for p in CONTROL_PATTERNS)


def clean_target_string(x):
    if pd.isna(x):
        return []
    s = str(x).strip()
    if s.lower() in ["", "nan", "none", "na"]:
        return []
    s = s.replace("[", "").replace("]", "").replace("'", "").replace('"', "")
    parts = re.split(r"[;,|/]+", s)
    return [p.strip() for p in parts if p.strip()]


def canonical_guide_id(x):
    return str(x).strip().replace(" ", "")


def get_guide_assignment_matrix(adata):
    ga = adata.obsm[GUIDE_ASSIGNMENT_KEY]
    if sparse.issparse(ga):
        return ga.tocsr().astype(bool)
    return np.asarray(ga).astype(bool)


def any_guides_per_cell(guide_assignment, guide_indices):
    if len(guide_indices) == 0:
        return np.zeros(guide_assignment.shape[0], dtype=bool)
    if sparse.issparse(guide_assignment):
        return np.asarray(guide_assignment[:, guide_indices].sum(axis=1)).ravel() > 0
    return guide_assignment[:, guide_indices].any(axis=1)


def classify_D3_EC_vs_nonEC(cell_type):
    if matches_any(cell_type, D3_EC_PATTERNS):
        return "artEC"
    return "non-EC"


def p_to_stars(p):
    if pd.isna(p):
        return "n.s."
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def safe_fisher(table):
    try:
        return fisher_exact(table, alternative="two-sided")[1]
    except Exception:
        return np.nan


def load_guide_target_from_metadata(adata):
    guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))
    guide_name_to_idx = {g: i for i, g in enumerate(guide_names)}

    meta = pd.read_csv(GUIDE_METADATA_TSV, sep="\t")
    if GUIDE_ID_COL not in meta.columns:
        raise ValueError(f"{GUIDE_ID_COL!r} not found. Columns: {meta.columns.tolist()}")

    target_col = None
    for c in TARGET_COL_CANDIDATES:
        if c in meta.columns:
            target_col = c
            break
    if target_col is None:
        raise ValueError(
            f"None of {TARGET_COL_CANDIDATES} found in {GUIDE_METADATA_TSV}. "
            f"Columns: {meta.columns.tolist()}"
        )

    meta = meta[[GUIDE_ID_COL, target_col]].copy()
    meta.columns = ["guide_id", "target_entry"]
    meta["guide_id"] = meta["guide_id"].astype(str)

    # Exact guide matching.
    meta["guide_index"] = meta["guide_id"].map(guide_name_to_idx)
    meta["match_method"] = np.where(meta["guide_index"].notna(), "exact_guide_id_match", None)

    # Fallback canonical matching.
    if meta["guide_index"].isna().any():
        canonical_to_guide = {}
        duplicated = set()
        for g in guide_names:
            cg = canonical_guide_id(g)
            if cg in canonical_to_guide:
                duplicated.add(cg)
            else:
                canonical_to_guide[cg] = g
        for cg in duplicated:
            canonical_to_guide.pop(cg, None)

        for idx in meta.index[meta["guide_index"].isna()]:
            cg = canonical_guide_id(meta.at[idx, "guide_id"])
            if cg in canonical_to_guide:
                matched = canonical_to_guide[cg]
                meta.at[idx, "guide_index"] = guide_name_to_idx[matched]
                meta.at[idx, "match_method"] = "canonical_guide_id_match"

    mapped = meta[meta["guide_index"].notna()].copy()
    mapped["guide_index"] = mapped["guide_index"].astype(int)
    mapped["guide_name"] = [guide_names[i] for i in mapped["guide_index"]]

    records = []
    for _, row in mapped.iterrows():
        guide_index = int(row["guide_index"])
        guide_name = str(row["guide_name"])
        raw_target = row["target_entry"]
        targets = clean_target_string(raw_target)
        if looks_control(guide_name) or looks_control(raw_target):
            targets = [CONTROL_TARGET]
        for target in targets:
            records.append({
                "guide_index": guide_index,
                "guide_name": guide_name,
                "target_name": str(target),
                "raw_target_entry": raw_target,
                "match_method": row["match_method"],
            })

    guide_target = pd.DataFrame(records).drop_duplicates()
    guide_target["target_name"] = guide_target["target_name"].astype(str)
    guide_target.to_csv(OUTDIR / "guide_to_target_map_for_panelD.csv", index=False)
    return guide_target


def count_for(counts, group, state):
    sub = counts[(counts["group"].astype(str) == group) & (counts["state"] == state)]
    if len(sub) == 0:
        return 0
    return int(sub["n_cells"].iloc[0])

# ============================================================
# Build guide/target masks
# ============================================================

guide_target = load_guide_target_from_metadata(adata)
guide_assignment = get_guide_assignment_matrix(adata)
guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))

if guide_assignment.shape[1] != len(guide_names):
    raise ValueError(
        f"guide_assignment has {guide_assignment.shape[1]} columns, "
        f"but adata.uns[{GUIDE_NAMES_KEY!r}] has {len(guide_names)} names."
    )

target_to_guide_indices = (
    guide_target.groupby("target_name")["guide_index"]
    .apply(lambda x: sorted(set(x.astype(int).tolist())))
    .to_dict()
)

target_masks = {
    target: any_guides_per_cell(guide_assignment, idx)
    for target, idx in target_to_guide_indices.items()
}

control_targets = [t for t in target_masks if looks_control(t) or str(t) == CONTROL_TARGET]
if len(control_targets) == 0:
    raise ValueError("No control-like guide targets detected.")

control_mask_raw = np.zeros(adata.n_obs, dtype=bool)
for t in control_targets:
    control_mask_raw |= target_masks[t]

# Hippo masks, by individual target.
hippo_target_masks = {}
missing_hippo = []
for target in HIPPO_TARGETS:
    matches = [t for t in target_masks if t.upper() == target.upper()]
    if len(matches) == 0:
        missing_hippo.append(target)
    else:
        hippo_target_masks[target] = target_masks[matches[0]].copy()

if missing_hippo:
    print("[warning] missing Hippo targets:", missing_hippo)

if len(hippo_target_masks) == 0:
    raise ValueError("No Hippo target masks found.")

hippo_mask_raw = np.zeros(adata.n_obs, dtype=bool)
for m in hippo_target_masks.values():
    hippo_mask_raw |= m

# Count targeting targets per cell, excluding NT controls.
targeting_targets = [
    t for t in target_masks
    if t not in control_targets and not looks_control(t)
]
n_targeting_targets = np.zeros(adata.n_obs, dtype=int)
for t in targeting_targets:
    n_targeting_targets += target_masks[t].astype(int)

# Control groups kept as before.
nt_excl_hippo_mask = control_mask_raw & (~hippo_mask_raw)
strict_nt_mask = control_mask_raw & (n_targeting_targets == 0)

# Individual Hippo groups.
group_masks = {}
for target, mask in hippo_target_masks.items():
    if EXCLUSIVE_HIPPO_TARGETS:
        other_hippo = np.zeros(adata.n_obs, dtype=bool)
        for other_target, other_mask in hippo_target_masks.items():
            if other_target != target:
                other_hippo |= other_mask
        group_masks[target] = mask & (~other_hippo)
    else:
        group_masks[target] = mask

group_masks["NT excl. Hippo"] = nt_excl_hippo_mask
group_masks["Strict NT"] = strict_nt_mask

qc = pd.DataFrame({
    "group": list(group_masks.keys()),
    "n_total_cells": [int(m.sum()) for m in group_masks.values()],
})
qc.to_csv(OUTDIR / "panelD_D3_EC_vs_nonEC_group_mask_qc.csv", index=False)
print("Group mask QC:")
display(qc)

# ============================================================
# Cell metadata and D3 EC/non-EC state
# ============================================================

obs = adata.obs.copy()
obs["cell_id"] = adata.obs_names.astype(str)
obs["day_clean"] = obs[DAY_COL].map(clean_day)
obs["cell_type_clean"] = obs[CELLTYPE_COL].astype(str)
obs["D3_EC_nonEC_state"] = obs["cell_type_clean"].apply(classify_D3_EC_vs_nonEC)

print("D3 cell-type labels:")
print(sorted(obs.loc[obs["day_clean"] == "D3", "cell_type_clean"].unique()))
print("D3 EC/non-EC counts:")
print(obs.loc[obs["day_clean"] == "D3", "D3_EC_nonEC_state"].value_counts())

# ============================================================
# Composition counts: D3 only, EC vs non-EC
# ============================================================

records = []
d3_day_mask = obs["day_clean"].values == "D3"

for group_name, group_mask in group_masks.items():
    use = d3_day_mask & group_mask
    tmp = obs.loc[use, ["cell_id", "cell_type_clean", "D3_EC_nonEC_state"]].copy()
    tmp = tmp.rename(columns={"D3_EC_nonEC_state": "state"})
    tmp["day"] = "D3"
    tmp["group"] = group_name
    records.append(tmp)

composition_cells = pd.concat(records, ignore_index=True)
composition_cells.to_csv(OUTDIR / "panelD_D3_EC_vs_nonEC_cell_membership.csv", index=False)

counts = (
    composition_cells.groupby(["day", "group", "state"], as_index=False)
    .size()
    .rename(columns={"size": "n_cells"})
)

# Complete missing group/state combinations.
complete_records = []
for group in GROUP_ORDER:
    if group not in group_masks:
        continue
    for state in D3_STATE_ORDER:
        sub = counts[(counts["group"] == group) & (counts["state"] == state)]
        if len(sub) == 0:
            complete_records.append({"day": "D3", "group": group, "state": state, "n_cells": 0})
        else:
            complete_records.append(sub.iloc[0].to_dict())

counts = pd.DataFrame(complete_records)
counts["total_cells"] = counts.groupby(["day", "group"])["n_cells"].transform("sum")
counts["fraction"] = counts["n_cells"] / counts["total_cells"].replace(0, np.nan)
counts["fraction"] = counts["fraction"].fillna(0)
counts["group"] = pd.Categorical(counts["group"], categories=[g for g in GROUP_ORDER if g in group_masks], ordered=True)
counts["state"] = pd.Categorical(counts["state"], categories=D3_STATE_ORDER, ordered=True)
counts = counts.sort_values(["group", "state"]).copy()
counts.to_csv(COMPOSITION_SOURCE_OUT, index=False)
print("Saved composition source:", COMPOSITION_SOURCE_OUT)
display(counts)

# ============================================================
# Fisher tests and FDR: each Hippo target vs each control group
# ============================================================

stats_records = []
for hippo_target in HIPPO_TARGETS:
    if hippo_target not in group_masks:
        continue
    for control_group in ["NT excl. Hippo", "Strict NT"]:
        h_ec = count_for(counts, hippo_target, "artEC")
        h_non = count_for(counts, hippo_target, "non-EC")
        c_ec = count_for(counts, control_group, "artEC")
        c_non = count_for(counts, control_group, "non-EC")
        p = safe_fisher([[h_ec, h_non], [c_ec, c_non]])
        stats_records.append({
            "day": "D3",
            "target": hippo_target,
            "comparison": f"{hippo_target} vs {control_group}",
            "test": "Fisher artEC vs non-EC",
            "target_artEC": h_ec,
            "target_nonEC": h_non,
            "control_artEC": c_ec,
            "control_nonEC": c_non,
            "p_value": p,
            "p_stars": p_to_stars(p),
        })

stats = pd.DataFrame(stats_records)
if len(stats):
    stats["fdr"] = multipletests(stats["p_value"].fillna(1), method="fdr_bh")[1]
    stats["fdr_stars"] = stats["fdr"].apply(p_to_stars)
else:
    stats["fdr"] = []
    stats["fdr_stars"] = []

stats.to_csv(STATS_SOURCE_OUT, index=False)
print("Saved stats source:", STATS_SOURCE_OUT)
display(stats)

# ============================================================
# Illustrator-editable stacked bar plot
# ============================================================

def make_panelD_D3_EC_vs_nonEC_individual_Hippo(counts, outdir=OUTDIR):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    plot_df = counts.copy()
    plot_df["group"] = pd.Categorical(plot_df["group"], categories=[g for g in GROUP_ORDER if g in group_masks], ordered=True)
    plot_df["state"] = pd.Categorical(plot_df["state"], categories=D3_STATE_ORDER, ordered=True)
    plot_df = plot_df.sort_values(["group", "state"]).copy()

    # Group labels with n.
    total_by_group = plot_df.groupby("group", observed=True)["n_cells"].sum().to_dict()
    group_order_present = [g for g in GROUP_ORDER if g in group_masks]

    group_label_map = {
        g: f"{g}\nn = {int(total_by_group.get(g, 0)):,}"
        for g in group_order_present
    }

    # Compact y positions: one row per group.
    # Smaller BAR_SPACING and larger BAR_HALF_HEIGHT reduce the white space
    # between horizontal bars.
    BAR_SPACING = 0.78
    BAR_HALF_HEIGHT = 0.34
    y_map = {
        g: (len(group_order_present) - i) * BAR_SPACING
        for i, g in enumerate(group_order_present)
    }
    y_min = min(y_map.values()) - BAR_HALF_HEIGHT - 0.16
    y_max = max(y_map.values()) + BAR_HALF_HEIGHT + 0.16

    plot_df["y"] = plot_df["group"].astype(str).map(y_map)
    plot_df["fraction_pct"] = 100 * plot_df["fraction"]
    plot_df["xmin"] = plot_df.groupby("group", observed=True)["fraction_pct"].cumsum() - plot_df["fraction_pct"]
    plot_df["xmax"] = plot_df.groupby("group", observed=True)["fraction_pct"].cumsum()
    plot_df["xmid"] = (plot_df["xmin"] + plot_df["xmax"]) / 2
    plot_df["ymin"] = plot_df["y"] - BAR_HALF_HEIGHT
    plot_df["ymax"] = plot_df["y"] + BAR_HALF_HEIGHT

    plot_df["segment_label"] = np.where(
        plot_df["fraction_pct"] >= 8,
        plot_df["fraction_pct"].round(0).astype(int).astype(str) + "%",
        ""
    )
    plot_df["label_color"] = np.where(plot_df["state"].astype(str) == "non-EC", "#111111", "white")

    label_white = plot_df[(plot_df["segment_label"] != "") & (plot_df["label_color"] == "white")].copy()
    label_black = plot_df[(plot_df["segment_label"] != "") & (plot_df["label_color"] == "#111111")].copy()

    outline_df = pd.DataFrame({
        "xmin": [0] * len(group_order_present),
        "xmax": [100] * len(group_order_present),
        "ymin": [y_map[g] - BAR_HALF_HEIGHT for g in group_order_present],
        "ymax": [y_map[g] + BAR_HALF_HEIGHT for g in group_order_present],
    })

    out_png = outdir / "panelD_D3_EC_vs_nonEC_individual_Hippo.png"
    out_pdf = outdir / "panelD_D3_EC_vs_nonEC_individual_Hippo_AI_editable.pdf"
    out_svg = outdir / "panelD_D3_EC_vs_nonEC_individual_Hippo_AI_editable.svg"

    p = (
        ggplot()
        + geom_rect(
            data=plot_df,
            mapping=aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="state"),
            inherit_aes=False,
            color=None,
            size=0
        )
        + geom_rect(
            data=outline_df,
            mapping=aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"),
            inherit_aes=False,
            fill=None,
            color="#555555",
            size=0.50
        )
        + geom_text(
            data=label_white,
            mapping=aes(x="xmid", y="y", label="segment_label"),
            inherit_aes=False,
            color="white",
            size=6.5,
            fontweight="bold",
            family=AI_FONT
        )
        + geom_text(
            data=label_black,
            mapping=aes(x="xmid", y="y", label="segment_label"),
            inherit_aes=False,
            color="#111111",
            size=6.5,
            fontweight="bold",
            family=AI_FONT
        )
        + scale_fill_manual(values=D3_STATE_COLORS, limits=D3_STATE_ORDER, name=None)
        + scale_x_continuous(
            limits=(0, 100.5),
            breaks=[0, 25, 50, 75, 100],
            labels=["0", "25", "50", "75", "100"],
            expand=(0, 0)
        )
        + scale_y_continuous(
            breaks=[y_map[g] for g in group_order_present],
            labels=[group_label_map[g] for g in group_order_present],
            limits=(y_min, y_max)
        )
        + coord_cartesian(xlim=(0, 100.5), ylim=(y_min, y_max), expand=False)
        + labs(
            title="D3 EC composition by Hippo target",
            subtitle="artEC vs non-EC among D3 cells",
            x="Cell fraction (%)",
            y=None
        )
        + theme_classic(base_size=11, base_family=AI_FONT)
        + theme(
            figure_size=(6.2, max(2.85, 0.42 * len(group_order_present) + 0.95)),
            panel_border=element_rect(fill=None, color="#6F6F6F", size=0.8),
            panel_grid_major_x=element_blank(),
            panel_grid_major_y=element_blank(),
            panel_grid_minor=element_blank(),
            axis_text_x=element_text(size=10.5, color="#303030", family=AI_FONT),
            axis_text_y=element_text(size=8.3, weight="bold", color="#303030", linespacing=0.95, family=AI_FONT),
            axis_title_x=element_text(size=12.0, weight="bold", color="#202020", family=AI_FONT),
            axis_title_y=element_blank(),
            plot_title=element_text(size=14.2, weight="bold", color="#111111", family=AI_FONT),
            plot_subtitle=element_text(size=10.6, color="#222222", family=AI_FONT),
            legend_position="top",
            legend_direction="horizontal",
            legend_title=element_blank(),
            legend_text=element_text(size=10.0, color="#222222", family=AI_FONT),
            plot_margin=0.03
        )
    )

    width = 6.2
    height = max(2.85, 0.42 * len(group_order_present) + 0.95)
    p.save(str(out_png), width=width, height=height, dpi=600, verbose=False)
    p.save(str(out_pdf), width=width, height=height, verbose=False)
    p.save(str(out_svg), width=width, height=height, verbose=False)

    print("Saved:")
    print("  PNG:", out_png)
    print("  PDF:", out_pdf)
    print("  SVG:", out_svg)

    return p

panelD_plot = make_panelD_D3_EC_vs_nonEC_individual_Hippo(counts, OUTDIR)
