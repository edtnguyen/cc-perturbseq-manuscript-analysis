#!/usr/bin/env python
# export_hippo_celllevel_adjusted_to_user_adata.py
#
# Adjusted to your pasted adata conventions:
# - adata.obsm["guide_assignment"]
# - adata.uns["guide_names"] if guide_assignment is not a DataFrame
# - adata.obsm["cnmf_usage"]
# - guide targets parsed from guide name, with non-targeting controls mapped to "non-targeting"
#
# Use in a notebook where adata is loaded:
#     %run export_hippo_celllevel_adjusted_to_user_adata.py
#
# Or copy cells from Hippo_celllevel_export_adjusted_to_user_adata.ipynb.


import os
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.stats import fisher_exact, ks_2samp, mannwhitneyu, wasserstein_distance, ttest_1samp, wilcoxon

# ============================================================
# CONFIG: adjusted to your pasted adata conventions
# ============================================================
DAY_COL = "day"
CELLTYPE_COL = "cell_type"
PROGRAM_SCORES_KEY = "cnmf_usage"
GUIDE_ASSIGNMENT_KEY = "guide_assignment"

DAY_D1 = "d1"
DAY_D2 = "d2"
DAY_D3 = "d3"

CONTROL_TARGET = "non-targeting"

# Your pasted code labels array-backed cnmf_usage columns as Program 0..Program n-1.
# If your adata.obsm["cnmf_usage"] is a numpy array, keep this as 0.
# If you discover that your array is actually 1-indexed, set this to 1.
ARRAY_PROGRAM_INDEX_BASE = 0

# Set to None to analyze all detected targets.
# Or restrict to a focused list, e.g. ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"].
SELECTED_TARGETS = None

HIPPO_TARGETS = [
    "NF2", "LATS1", "LATS2", "SAV1", "STK3", "STK4",
    "MOB1A", "MOB1B", "TAOK1", "TAOK2", "TAOK3",
    "AMOT", "AMOTL1", "AMOTL2"
]

# Minimum counts
MIN_TARGET_CELLS = 20
MIN_GUIDE_CELLS = 8
MIN_REF_CELLS = 50

# D1 program-state comparison
D1_PROGRAMS = [13, 30]

# D2 state comparison
D2_PRE_EC_LABELS = ["Pre-EC"]
D2_LAT_MESO_LABELS = ["Lat-Meso"]

# D3 state comparison
D3_ENDO_LABELS = ["artEC", "ECM-EC", "Endothelial", "EC", "Angio-EC", "Partial EndMT EC"]
D3_MES_LABELS = ["Fibro-Mesen", "Mesenchymal"]
D3_NEU_MES_LABELS = ["Neu-Mesen", "NeuCrest/Mes"]

# Optional cNMF program summaries
EXPORT_PROGRAM_SUMMARY = True
D2_PROGRAMS = [14, 15, 17, 19, 23, 34]
D3_PROGRAMS = [2, 3, 8, 11, 27, 34, 35]
ALL_PROGRAMS_FOR_AXES = sorted(set(D1_PROGRAMS + D2_PROGRAMS + D3_PROGRAMS))

# Fate-axis definitions. Positive = alternative / failure bias.
FATE_AXIS_DEFS = {
    "D1_alt_P13_minus_P30": {
        "day": DAY_D1,
        "plus": [13],
        "minus": [30],
        "description": "D1 P13-over-P30 alternative/non-EC early bias",
    },
    "D2_alt_P17_minus_P15": {
        "day": DAY_D2,
        "plus": [17],
        "minus": [15],
        "description": "D2 P17-over-P15 Lat-Meso/alternative bias",
    },
    "D2_YAPLat_minus_PreEC": {
        "day": DAY_D2,
        "plus": [17, 34],
        "minus": [15, 23],
        "description": "D2 mean(P17,P34)-mean(P15,P23) YAP/Lat-Meso vs Pre-EC bias",
    },
    "D3_ALT_minus_EC": {
        "day": DAY_D3,
        "plus": [27, 34, 35],
        "minus": [2, 3, 8, 11],
        "description": "D3 mean(P27,P34,P35)-mean(P2,P3,P8,P11) alternative endpoint bias",
    },
}

OUTDIR = "./hippo_adjusted_celllevel_exports"
os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# HELPERS: adjusted from your pasted code
# ============================================================
def get_guide_names(adata):
    ga = adata.obsm[GUIDE_ASSIGNMENT_KEY]
    if isinstance(ga, pd.DataFrame):
        return np.array(ga.columns.astype(str))
    if "guide_names" in adata.uns:
        return np.array(adata.uns["guide_names"]).astype(str)
    raise KeyError("Could not find guide names in adata.uns['guide_names'].")


def parse_guide_target(guide_name):
    g = str(guide_name)
    gl = g.lower()

    if "non-target" in gl or gl.startswith("ntc"):
        return CONTROL_TARGET

    # Your guide names seem to begin with TARGET_
    if "_" in g:
        return g.split("_")[0]

    for sep in ["-", ".", " "]:
        if sep in g:
            return g.split(sep)[0]

    return g


def get_guide_targets(adata):
    guide_names = get_guide_names(adata)
    return np.array([parse_guide_target(g) for g in guide_names], dtype=str)


def get_guide_assignment_df(adata):
    ga = adata.obsm[GUIDE_ASSIGNMENT_KEY]
    guide_names = get_guide_names(adata)

    if isinstance(ga, pd.DataFrame):
        out = ga.copy()
        out.index = adata.obs_names
        out.columns = out.columns.astype(str)
        return out

    if issparse(ga):
        ga = ga.toarray()

    return pd.DataFrame(ga, index=adata.obs_names, columns=guide_names)


def get_target_assignment_df(adata):
    ga = get_guide_assignment_df(adata)
    guide_targets = get_guide_targets(adata)

    out = pd.DataFrame(index=adata.obs_names)
    for target in pd.unique(guide_targets):
        cols = np.where(guide_targets == target)[0]
        out[target] = (ga.iloc[:, cols].sum(axis=1).values > 0)
    return out


def day_mask(adata, day_value):
    return adata.obs[DAY_COL].astype(str).str.lower().values == str(day_value).lower()


def _standardize_program_col(c):
    """
    Convert possible cNMF column names to 'Program N'.
    Supports:
      - P14
      - Program 14
      - numeric column 14
      - string '14'
    """
    c0 = c
    c = str(c)
    if c.startswith("Program "):
        return f"Program {int(c.replace('Program ', ''))}"
    if c.startswith("P"):
        return f"Program {int(c.replace('P', ''))}"
    return f"Program {int(float(c))}"


def get_program_score_df(adata, program_scores_key=PROGRAM_SCORES_KEY, array_program_index_base=ARRAY_PROGRAM_INDEX_BASE):
    X = adata.obsm[program_scores_key]

    if isinstance(X, pd.DataFrame):
        df = X.copy()
        df.columns = [_standardize_program_col(c) for c in df.columns]
        df.index = adata.obs_names
        return df

    if issparse(X):
        X = X.toarray()

    X = np.asarray(X)
    cols = [f"Program {i + array_program_index_base}" for i in range(X.shape[1])]
    return pd.DataFrame(X, index=adata.obs_names, columns=cols)


def zscore_relative_to_reference(score_df, ref_mask, all_mask, programs, min_ref_cells=MIN_REF_CELLS):
    missing = [p for p in programs if p not in score_df.columns]
    if missing:
        raise KeyError(
            f"Missing program columns: {missing}. "
            f"Available program range example: {list(score_df.columns[:5])} ... {list(score_df.columns[-5:])}. "
            f"If adata.obsm['{PROGRAM_SCORES_KEY}'] is an array and indexing is off, try ARRAY_PROGRAM_INDEX_BASE=1."
        )

    if ref_mask.sum() >= min_ref_cells:
        ref = score_df.loc[score_df.index[ref_mask], programs]
        ref_label = f"{CONTROL_TARGET} carriers"
    else:
        ref = score_df.loc[score_df.index[all_mask], programs]
        ref_label = "all cells at day"

    mu = ref.mean(axis=0)
    sd = ref.std(axis=0, ddof=0).replace(0, 1.0)
    z = (score_df[programs] - mu) / sd
    return z, ref_label


def safe_frac(n, d):
    return n / d if d > 0 else np.nan


def log2_or(a, b, c, d, pseudocount=0.5):
    return np.log2(((a + pseudocount) / (b + pseudocount)) / ((c + pseudocount) / (d + pseudocount)))


def map_d1_program_state(row):
    if pd.isna(row["Program 13"]) or pd.isna(row["Program 30"]):
        return pd.NA
    return "P13_high" if row["Program 13"] > row["Program 30"] else "P30_high"


def map_d2_program_state(row):
    if pd.isna(row["Program 17"]) or pd.isna(row["Program 15"]):
        return pd.NA
    return "P17_high" if row["Program 17"] > row["Program 15"] else "P15_high"


def map_d2_state(x):
    x = str(x)
    if x in D2_PRE_EC_LABELS:
        return "Pre-EC"
    if x in D2_LAT_MESO_LABELS:
        return "Lat-Meso"
    return pd.NA


def map_d3_state(x):
    x = str(x)
    if x in D3_ENDO_LABELS:
        return "Endothelial"
    if x in D3_MES_LABELS:
        return "Mesenchymal"
    if x in D3_NEU_MES_LABELS:
        return "Neu/Mesen"
    return pd.NA


def state_counts(df, state_col, positive_states, negative_state):
    vc = df[state_col].value_counts()
    n_pos = int(sum(vc.get(s, 0) for s in positive_states))
    n_neg = int(vc.get(negative_state, 0))
    n_total = n_pos + n_neg
    return n_pos, n_neg, n_total


def cliffs_delta(x, y, max_pairs=2_000_000, random_seed=1):
    x = np.asarray(pd.Series(x).dropna())
    y = np.asarray(pd.Series(y).dropna())
    if len(x) == 0 or len(y) == 0:
        return np.nan

    # Downsample if too large
    if len(x) * len(y) > max_pairs:
        rng = np.random.default_rng(random_seed)
        x = rng.choice(x, size=min(len(x), 1500), replace=False)
        y = rng.choice(y, size=min(len(y), 1500), replace=False)

    gt = sum((xi > y).sum() for xi in x)
    lt = sum((xi < y).sum() for xi in x)
    return (gt - lt) / (len(x) * len(y))


def summarize_distribution_shift(target, mode_name, day_name, score_name, target_values, ctrl_values):
    x = pd.Series(target_values).dropna().values
    y = pd.Series(ctrl_values).dropna().values

    row = {
        "mode": mode_name,
        "target": target,
        "day": day_name,
        "score": score_name,
        "n_target": len(x),
        "n_control": len(y),
        "target_mean": np.mean(x) if len(x) else np.nan,
        "control_mean": np.mean(y) if len(y) else np.nan,
        "delta_mean": (np.mean(x) - np.mean(y)) if len(x) and len(y) else np.nan,
        "target_median": np.median(x) if len(x) else np.nan,
        "control_median": np.median(y) if len(y) else np.nan,
        "delta_median": (np.median(x) - np.median(y)) if len(x) and len(y) else np.nan,
        "wasserstein": wasserstein_distance(x, y) if len(x) and len(y) else np.nan,
        "cliffs_delta": cliffs_delta(x, y) if len(x) and len(y) else np.nan,
        "ks_p": ks_2samp(x, y).pvalue if len(x) >= 5 and len(y) >= 5 else np.nan,
        "mannwhitney_p": mannwhitneyu(x, y, alternative="two-sided").pvalue if len(x) >= 5 and len(y) >= 5 else np.nan,
    }
    return row


def run_binary_enrichment(
    obs,
    guide_names,
    guide_targets,
    guide_assignment,
    target_assignment,
    selected_targets,
    base_mask,
    selector,
    state_col,
    positive_states,
    negative_state,
    positive_name,
    negative_name,
    mode_name,
    day_name,
    min_target_cells=MIN_TARGET_CELLS,
    min_guide_cells=MIN_GUIDE_CELLS,
):
    target_rows = []
    guide_rows = []

    ctrl_mask = (
        (target_assignment[CONTROL_TARGET].values if CONTROL_TARGET in target_assignment.columns
         else np.zeros(obs.shape[0], dtype=bool))
        & base_mask
        & selector
    )
    ctrl_cells = obs.loc[ctrl_mask].copy()
    c_pos, c_neg, c_total = state_counts(ctrl_cells, state_col, positive_states, negative_state)

    for target in selected_targets:
        if target not in target_assignment.columns or target == CONTROL_TARGET:
            continue

        tmask = target_assignment[target].values & base_mask & selector
        sub = obs.loc[tmask].copy()
        n_pos, n_neg, n_total = state_counts(sub, state_col, positive_states, negative_state)

        if n_total >= min_target_cells:
            row = {
                "mode": mode_name,
                "day": day_name,
                "target": target,
                "n_total": n_total,
                f"n_{positive_name}": n_pos,
                f"n_{negative_name}": n_neg,
                f"frac_{positive_name}": safe_frac(n_pos, n_total),
                f"frac_{negative_name}": safe_frac(n_neg, n_total),
                f"ctrl_n_{positive_name}": c_pos,
                f"ctrl_n_{negative_name}": c_neg,
                f"ctrl_frac_{positive_name}": safe_frac(c_pos, c_total),
                f"ctrl_frac_{negative_name}": safe_frac(c_neg, c_total),
                f"delta_frac_{positive_name}_vs_ctrl": safe_frac(n_pos, n_total) - safe_frac(c_pos, c_total),
                f"delta_frac_{negative_name}_vs_ctrl": safe_frac(n_neg, n_total) - safe_frac(c_neg, c_total),
                f"log2OR_{positive_name}_vs_{negative_name}": log2_or(n_pos, n_neg, c_pos, c_neg),
                "median_total_guides": float(sub["n_total_guides"].median()) if n_total > 0 else np.nan,
                "median_total_targets": float(sub["n_total_targets"].median()) if n_total > 0 else np.nan,
                "median_selected_targets": float(sub["n_selected_targets"].median()) if n_total > 0 else np.nan,
            }

            try:
                row[f"fisher_p_{positive_name}_vs_{negative_name}"] = fisher_exact(
                    [[n_pos, n_neg], [c_pos, c_neg]],
                    alternative="two-sided"
                )[1]
            except Exception:
                row[f"fisher_p_{positive_name}_vs_{negative_name}"] = np.nan

            target_rows.append(row)

    target_df = pd.DataFrame(target_rows)

    for j, guide in enumerate(guide_names):
        target = guide_targets[j]
        if target not in selected_targets or target == CONTROL_TARGET:
            continue

        carrier = (guide_assignment.iloc[:, j].values > 0) & base_mask & selector
        sub = obs.loc[carrier].copy()
        n_pos, n_neg, n_total = state_counts(sub, state_col, positive_states, negative_state)

        if n_total >= min_guide_cells:
            row = {
                "mode": mode_name,
                "day": day_name,
                "guide": guide,
                "target": target,
                "n_total": n_total,
                f"n_{positive_name}": n_pos,
                f"n_{negative_name}": n_neg,
                f"frac_{positive_name}": safe_frac(n_pos, n_total),
                f"frac_{negative_name}": safe_frac(n_neg, n_total),
                f"ctrl_n_{positive_name}": c_pos,
                f"ctrl_n_{negative_name}": c_neg,
                f"ctrl_frac_{positive_name}": safe_frac(c_pos, c_total),
                f"ctrl_frac_{negative_name}": safe_frac(c_neg, c_total),
                f"delta_frac_{positive_name}_vs_ctrl": safe_frac(n_pos, n_total) - safe_frac(c_pos, c_total),
                f"log2OR_{positive_name}_vs_{negative_name}": log2_or(n_pos, n_neg, c_pos, c_neg),
            }
            try:
                row[f"fisher_p_{positive_name}_vs_{negative_name}"] = fisher_exact(
                    [[n_pos, n_neg], [c_pos, c_neg]],
                    alternative="two-sided"
                )[1]
            except Exception:
                row[f"fisher_p_{positive_name}_vs_{negative_name}"] = np.nan
            guide_rows.append(row)

    guide_df = pd.DataFrame(guide_rows)
    return target_df, guide_df


# The following block expects an existing variable named `adata`.
if "adata" not in globals():
    raise RuntimeError("This script expects `adata` to already be loaded in the current Python session.")


# ============================================================
# BUILD BASE OBJECTS
# ============================================================
guide_names = get_guide_names(adata)
guide_targets = get_guide_targets(adata)
guide_assignment = get_guide_assignment_df(adata)
target_assignment = get_target_assignment_df(adata)

obs = adata.obs.copy()
obs["cell_id"] = adata.obs_names
obs["cell_type_str"] = obs[CELLTYPE_COL].astype(str)
obs["day_lower"] = obs[DAY_COL].astype(str).str.lower()

obs["d2_state"] = obs["cell_type_str"].map(map_d2_state).astype("object")
obs["d3_state"] = obs["cell_type_str"].map(map_d3_state).astype("object")

mask_d1 = day_mask(adata, DAY_D1)
mask_d2 = day_mask(adata, DAY_D2)
mask_d3 = day_mask(adata, DAY_D3)

mask_d2_relevant = mask_d2 & obs["d2_state"].notna().values
mask_d3_relevant = mask_d3 & obs["d3_state"].notna().values

obs["n_total_guides"] = (guide_assignment.values > 0).sum(axis=1)
obs["n_total_targets"] = target_assignment.sum(axis=1).values

if SELECTED_TARGETS is None:
    selected_targets = sorted(list(pd.unique(guide_targets)))
else:
    selected_targets = list(SELECTED_TARGETS)

selected_targets = [t for t in selected_targets if t in target_assignment.columns]
obs["n_selected_targets"] = target_assignment[selected_targets].sum(axis=1).values

hippo_present = [g for g in HIPPO_TARGETS if g in target_assignment.columns]

# Save base metadata
pd.DataFrame({"guide": guide_names, "target": guide_targets}).to_csv(
    os.path.join(OUTDIR, "guide_to_target_map.csv"), index=False
)
pd.Series(guide_targets, name="target").value_counts().to_csv(
    os.path.join(OUTDIR, "guide_target_counts.csv")
)
obs[[
    "cell_id", DAY_COL, CELLTYPE_COL, "day_lower", "d2_state", "d3_state",
    "n_total_guides", "n_total_targets", "n_selected_targets"
]].to_csv(os.path.join(OUTDIR, "cell_moi_metadata.csv"), index=False)

with open(os.path.join(OUTDIR, "hippo_targets_present.txt"), "w") as f:
    f.write("\n".join(hippo_present) + "\n")

print("n cells:", adata.n_obs)
print("n guides:", len(guide_names))
print("n targets:", len(selected_targets))
print("Hippo targets present:", hippo_present)
print("Control target present:", CONTROL_TARGET in target_assignment.columns)


# ============================================================
# PROGRAM SCORE SUMMARIES
# ============================================================
score_df = get_program_score_df(adata)

# Required columns.
d1_prog_labels = [f"Program {i}" for i in D1_PROGRAMS]
d2_prog_labels = [f"Program {i}" for i in D2_PROGRAMS]
d3_prog_labels = [f"Program {i}" for i in D3_PROGRAMS]
all_axis_program_labels = [f"Program {i}" for i in ALL_PROGRAMS_FOR_AXES]

missing = [p for p in all_axis_program_labels if p not in score_df.columns]
if missing:
    raise KeyError(
        f"Missing required program columns: {missing}. "
        f"First columns: {list(score_df.columns[:8])}; last columns: {list(score_df.columns[-8:])}. "
        f"If your numpy array uses 1-based program numbering, set ARRAY_PROGRAM_INDEX_BASE = 1 and rerun."
    )

ctrl_d1 = (
    (target_assignment[CONTROL_TARGET].values if CONTROL_TARGET in target_assignment.columns
     else np.zeros(adata.n_obs, dtype=bool))
    & mask_d1
)
ctrl_d2 = (
    (target_assignment[CONTROL_TARGET].values if CONTROL_TARGET in target_assignment.columns
     else np.zeros(adata.n_obs, dtype=bool))
    & mask_d2
)
ctrl_d3 = (
    (target_assignment[CONTROL_TARGET].values if CONTROL_TARGET in target_assignment.columns
     else np.zeros(adata.n_obs, dtype=bool))
    & mask_d3
)

z_d1, d1_ref = zscore_relative_to_reference(score_df, ctrl_d1, mask_d1, d1_prog_labels)
z_d2, d2_ref = zscore_relative_to_reference(score_df, ctrl_d2, mask_d2, d2_prog_labels)
z_d3, d3_ref = zscore_relative_to_reference(score_df, ctrl_d3, mask_d3, d3_prog_labels)

# Combine z-scored program columns into obs.
for c in z_d1.columns:
    obs[c] = np.nan
    obs.loc[z_d1.index[mask_d1], c] = z_d1.loc[z_d1.index[mask_d1], c]

for c in z_d2.columns:
    if c not in obs.columns:
        obs[c] = np.nan
    obs.loc[z_d2.index[mask_d2], c] = z_d2.loc[z_d2.index[mask_d2], c]

for c in z_d3.columns:
    if c not in obs.columns:
        obs[c] = np.nan
    obs.loc[z_d3.index[mask_d3], c] = z_d3.loc[z_d3.index[mask_d3], c]

# D1 and D2 program-state labels.
obs["d1_program_state"] = obs.apply(
    lambda row: map_d1_program_state(row) if str(row["day_lower"]) == DAY_D1 else pd.NA,
    axis=1
).astype("object")

obs["d2_program_state"] = obs.apply(
    lambda row: map_d2_program_state(row) if str(row["day_lower"]) == DAY_D2 else pd.NA,
    axis=1
).astype("object")

# Fate-axis scores. Positive = alternative / failure bias.
for axis_name, spec in FATE_AXIS_DEFS.items():
    plus_cols = [f"Program {i}" for i in spec["plus"]]
    minus_cols = [f"Program {i}" for i in spec["minus"]]
    obs[axis_name] = obs[plus_cols].mean(axis=1) - obs[minus_cols].mean(axis=1)

obs.to_csv(os.path.join(OUTDIR, "cell_metadata_with_program_zscores_and_axes.csv"), index=False)

ref_summary = pd.DataFrame([
    {"day": "d1", "reference_used": d1_ref, "n_control_carriers": int(ctrl_d1.sum())},
    {"day": "d2", "reference_used": d2_ref, "n_control_carriers": int(ctrl_d2.sum())},
    {"day": "d3", "reference_used": d3_ref, "n_control_carriers": int(ctrl_d3.sum())},
])
ref_summary.to_csv(os.path.join(OUTDIR, "program_zscore_reference_summary.csv"), index=False)
print(ref_summary)


# ============================================================
# RUN TWO MODES
# ============================================================
modes = {
    "all_carriers": np.ones(adata.n_obs, dtype=bool),
    "single_selected_target": obs["n_selected_targets"].values == 1,
}

all_leadlag = []
axis_distribution_rows = []
composition_rows = []
guide_axis_rows = []

for mode_name, selector in modes.items():

    # D1 P13-high vs P30-high, program-defined
    mask_d1_relevant = mask_d1 & obs["d1_program_state"].notna().values
    d1_target, d1_guide = run_binary_enrichment(
        obs=obs,
        guide_names=guide_names,
        guide_targets=guide_targets,
        guide_assignment=guide_assignment,
        target_assignment=target_assignment,
        selected_targets=selected_targets,
        base_mask=mask_d1_relevant,
        selector=selector,
        state_col="d1_program_state",
        positive_states=["P13_high"],
        negative_state="P30_high",
        positive_name="P13high",
        negative_name="P30high",
        mode_name=mode_name,
        day_name="d1",
    )

    # D2 Lat-Meso vs Pre-EC
    d2_target, d2_guide = run_binary_enrichment(
        obs=obs,
        guide_names=guide_names,
        guide_targets=guide_targets,
        guide_assignment=guide_assignment,
        target_assignment=target_assignment,
        selected_targets=selected_targets,
        base_mask=mask_d2_relevant,
        selector=selector,
        state_col="d2_state",
        positive_states=["Lat-Meso"],
        negative_state="Pre-EC",
        positive_name="latMeso",
        negative_name="preEC",
        mode_name=mode_name,
        day_name="d2",
    )

    # D2 P17-high vs P15-high, program-defined
    mask_d2_program_relevant = mask_d2 & obs["d2_program_state"].notna().values
    d2_prog_target, d2_prog_guide = run_binary_enrichment(
        obs=obs,
        guide_names=guide_names,
        guide_targets=guide_targets,
        guide_assignment=guide_assignment,
        target_assignment=target_assignment,
        selected_targets=selected_targets,
        base_mask=mask_d2_program_relevant,
        selector=selector,
        state_col="d2_program_state",
        positive_states=["P17_high"],
        negative_state="P15_high",
        positive_name="P17high",
        negative_name="P15high",
        mode_name=mode_name,
        day_name="d2",
    )

    # D3 broad alternative vs endothelial
    d3_target_alt, d3_guide_alt = run_binary_enrichment(
        obs=obs,
        guide_names=guide_names,
        guide_targets=guide_targets,
        guide_assignment=guide_assignment,
        target_assignment=target_assignment,
        selected_targets=selected_targets,
        base_mask=mask_d3_relevant,
        selector=selector,
        state_col="d3_state",
        positive_states=["Mesenchymal", "Neu/Mesen"],
        negative_state="Endothelial",
        positive_name="alt",
        negative_name="endo",
        mode_name=mode_name,
        day_name="d3",
    )

    # D3 strict mesenchymal vs endothelial
    d3_target_mes, d3_guide_mes = run_binary_enrichment(
        obs=obs,
        guide_names=guide_names,
        guide_targets=guide_targets,
        guide_assignment=guide_assignment,
        target_assignment=target_assignment,
        selected_targets=selected_targets,
        base_mask=mask_d3_relevant,
        selector=selector,
        state_col="d3_state",
        positive_states=["Mesenchymal"],
        negative_state="Endothelial",
        positive_name="mes",
        negative_name="endo",
        mode_name=mode_name,
        day_name="d3",
    )

    # Program summaries for D1 target table.
    if EXPORT_PROGRAM_SUMMARY and len(d1_target):
        for prog in d1_prog_labels:
            vals = []
            for target in d1_target["target"]:
                tmask = target_assignment[target].values & mask_d1 & selector
                vals.append(obs.loc[tmask, prog].mean())
            d1_target[f"{prog}_mean_z"] = vals

        d1_target["d1_alt_P13_minus_P30"] = d1_target["Program 13_mean_z"] - d1_target["Program 30_mean_z"]

    # Program summaries for D2 target table.
    if EXPORT_PROGRAM_SUMMARY and len(d2_target):
        for prog in d2_prog_labels:
            vals = []
            for target in d2_target["target"]:
                tmask = target_assignment[target].values & mask_d2 & selector
                vals.append(obs.loc[tmask, prog].mean())
            d2_target[f"{prog}_mean_z"] = vals

        d2_target["d2_bifurcation_score"] = (
            d2_target[["Program 17_mean_z", "Program 34_mean_z"]].mean(axis=1)
            - d2_target[["Program 15_mean_z", "Program 23_mean_z"]].mean(axis=1)
        )
        d2_target["d2_alt_P17_minus_P15"] = d2_target["Program 17_mean_z"] - d2_target["Program 15_mean_z"]

    # Program summaries for D3 target table.
    if EXPORT_PROGRAM_SUMMARY and len(d3_target_alt):
        for prog in d3_prog_labels:
            vals = []
            for target in d3_target_alt["target"]:
                tmask = target_assignment[target].values & mask_d3 & selector
                vals.append(obs.loc[tmask, prog].mean())
            d3_target_alt[f"{prog}_mean_z"] = vals

        d3_target_alt["d3_ec_endpoint_score"] = d3_target_alt[
            ["Program 2_mean_z", "Program 3_mean_z", "Program 8_mean_z", "Program 11_mean_z"]
        ].mean(axis=1)

        d3_target_alt["d3_alt_endpoint_score"] = d3_target_alt[
            ["Program 27_mean_z", "Program 34_mean_z", "Program 35_mean_z"]
        ].mean(axis=1)

        d3_target_alt["d3_endpoint_score_EC_minus_ALT"] = (
            d3_target_alt["d3_ec_endpoint_score"] - d3_target_alt["d3_alt_endpoint_score"]
        )
        d3_target_alt["d3_alt_minus_EC"] = (
            d3_target_alt["d3_alt_endpoint_score"] - d3_target_alt["d3_ec_endpoint_score"]
        )

    # Save individual enrichment tables.
    d1_target.to_csv(os.path.join(OUTDIR, f"d1_target_P13high_vs_P30high_{mode_name}.csv"), index=False)
    d1_guide.to_csv(os.path.join(OUTDIR, f"d1_guide_P13high_vs_P30high_{mode_name}.csv"), index=False)

    d2_target.to_csv(os.path.join(OUTDIR, f"d2_target_preEC_vs_latMeso_{mode_name}.csv"), index=False)
    d2_guide.to_csv(os.path.join(OUTDIR, f"d2_guide_preEC_vs_latMeso_{mode_name}.csv"), index=False)

    d2_prog_target.to_csv(os.path.join(OUTDIR, f"d2_target_P17high_vs_P15high_{mode_name}.csv"), index=False)
    d2_prog_guide.to_csv(os.path.join(OUTDIR, f"d2_guide_P17high_vs_P15high_{mode_name}.csv"), index=False)

    d3_target_alt.to_csv(os.path.join(OUTDIR, f"d3_target_alt_vs_endo_{mode_name}.csv"), index=False)
    d3_guide_alt.to_csv(os.path.join(OUTDIR, f"d3_guide_alt_vs_endo_{mode_name}.csv"), index=False)

    d3_target_mes.to_csv(os.path.join(OUTDIR, f"d3_target_mes_vs_endo_{mode_name}.csv"), index=False)
    d3_guide_mes.to_csv(os.path.join(OUTDIR, f"d3_guide_mes_vs_endo_{mode_name}.csv"), index=False)

    # Lead-lag tables: D1 -> D2 -> D3.
    if len(d1_target) and len(d2_target) and len(d3_target_alt):
        leadlag = d1_target.merge(
            d2_target,
            on=["mode", "target"],
            how="inner",
            suffixes=("_d1", "_d2")
        ).merge(
            d3_target_alt,
            on=["mode", "target"],
            how="inner",
            suffixes=("", "_d3_alt")
        )

        # Merge strict mesenchymal log2OR.
        cols_mes = ["mode", "target"]
        for c in ["log2OR_mes_vs_endo", "frac_mes", "frac_endo", "delta_frac_mes_vs_ctrl"]:
            if c in d3_target_mes.columns:
                cols_mes.append(c)
        leadlag = leadlag.merge(d3_target_mes[cols_mes], on=["mode", "target"], how="left")
        leadlag.to_csv(os.path.join(OUTDIR, f"d1_d2_d3_target_leadlag_{mode_name}.csv"), index=False)
        all_leadlag.append(leadlag)

    # Distribution shift tests for all selected targets and fate-axis scores.
    for target in selected_targets:
        if target == CONTROL_TARGET or target not in target_assignment.columns:
            continue

        carrier = target_assignment[target].values & selector
        ctrl = (
            (target_assignment[CONTROL_TARGET].values if CONTROL_TARGET in target_assignment.columns
             else np.zeros(adata.n_obs, dtype=bool))
            & selector
        )

        for axis_name, spec in FATE_AXIS_DEFS.items():
            day = spec["day"]
            day_m = obs["day_lower"].eq(day).values
            target_values = obs.loc[carrier & day_m, axis_name]
            ctrl_values = obs.loc[ctrl & day_m, axis_name]
            if len(pd.Series(target_values).dropna()) >= MIN_TARGET_CELLS:
                axis_distribution_rows.append(
                    summarize_distribution_shift(target, mode_name, day, axis_name, target_values, ctrl_values)
                )

    # Guide-level fate-axis effects.
    for j, guide in enumerate(guide_names):
        target = guide_targets[j]
        if target == CONTROL_TARGET or target not in selected_targets:
            continue

        carrier = (guide_assignment.iloc[:, j].values > 0) & selector
        ctrl = (
            (target_assignment[CONTROL_TARGET].values if CONTROL_TARGET in target_assignment.columns
             else np.zeros(adata.n_obs, dtype=bool))
            & selector
        )

        for axis_name, spec in FATE_AXIS_DEFS.items():
            day = spec["day"]
            day_m = obs["day_lower"].eq(day).values
            x = obs.loc[carrier & day_m, axis_name].dropna()
            y = obs.loc[ctrl & day_m, axis_name].dropna()
            if len(x) >= MIN_GUIDE_CELLS:
                guide_axis_rows.append({
                    "mode": mode_name,
                    "guide": guide,
                    "target": target,
                    "day": day,
                    "score": axis_name,
                    "n_guide_cells": len(x),
                    "n_control_cells": len(y),
                    "guide_mean": x.mean(),
                    "control_mean": y.mean() if len(y) else np.nan,
                    "effect_vs_control": x.mean() - y.mean() if len(y) else np.nan,
                    "guide_median": x.median(),
                    "control_median": y.median() if len(y) else np.nan,
                    "is_hippo_target": target in hippo_present,
                })

# Combined lead-lag table.
if all_leadlag:
    pd.concat(all_leadlag, axis=0, ignore_index=True).to_csv(
        os.path.join(OUTDIR, "d1_d2_d3_target_leadlag_all_modes.csv"),
        index=False
    )

axis_dist = pd.DataFrame(axis_distribution_rows)
axis_dist["is_hippo_target"] = axis_dist["target"].isin(hippo_present) if len(axis_dist) else []
axis_dist.to_csv(os.path.join(OUTDIR, "target_level_fate_axis_distribution_tests.csv"), index=False)

guide_axis = pd.DataFrame(guide_axis_rows)
guide_axis.to_csv(os.path.join(OUTDIR, "guide_level_fate_axis_effects.csv"), index=False)


# ============================================================
# CORRELATION SUMMARIES
# ============================================================
corr_rows = []
for mode_name in modes:
    path = os.path.join(OUTDIR, f"d1_d2_d3_target_leadlag_{mode_name}.csv")
    if not os.path.exists(path):
        continue

    tab = pd.read_csv(path)

    pairs = [
        ("log2OR_P13high_vs_P30high", "log2OR_latMeso_vs_preEC"),
        ("log2OR_P13high_vs_P30high", "log2OR_alt_vs_endo"),
        ("log2OR_latMeso_vs_preEC", "log2OR_alt_vs_endo"),
        ("log2OR_latMeso_vs_preEC", "log2OR_mes_vs_endo"),
        ("d1_alt_P13_minus_P30", "d2_alt_P17_minus_P15"),
        ("d2_alt_P17_minus_P15", "d3_alt_minus_EC"),
        ("d2_bifurcation_score", "d3_alt_minus_EC"),
    ]

    for x, y in pairs:
        if x in tab.columns and y in tab.columns:
            tmp = tab[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
            corr_rows.append({
                "mode": mode_name,
                "x": x,
                "y": y,
                "pearson_r": tmp.corr().iloc[0, 1] if len(tmp) >= 3 else np.nan,
                "spearman_r": tmp.corr(method="spearman").iloc[0, 1] if len(tmp) >= 3 else np.nan,
                "n_targets": len(tmp),
            })

pd.DataFrame(corr_rows).to_csv(os.path.join(OUTDIR, "d1_d2_d3_correlation_summary.csv"), index=False)


# ============================================================
# HIPPO-SPECIFIC EXPORTS
# ============================================================
# Hippo subset of key target-level tables.
for fn in os.listdir(OUTDIR):
    if fn.endswith(".csv") and ("target" in fn or "leadlag" in fn or "distribution" in fn):
        path = os.path.join(OUTDIR, fn)
        try:
            t = pd.read_csv(path)
        except Exception:
            continue
        if "target" in t.columns:
            sub = t[t["target"].isin(hippo_present)].copy()
            if len(sub):
                sub.to_csv(os.path.join(OUTDIR, "hippo_only__" + fn), index=False)

# Hippo carrier cell-level long table for plotting.
hippo_cell_rows = []
rng = np.random.default_rng(1)

for target in hippo_present:
    carrier = target_assignment[target].values
    ctrl = target_assignment[CONTROL_TARGET].values if CONTROL_TARGET in target_assignment.columns else np.zeros(adata.n_obs, dtype=bool)

    # all hippo carriers + sampled controls
    carrier_idx = np.where(carrier)[0]
    ctrl_idx = np.where(ctrl)[0]
    if len(ctrl_idx) > 15000:
        ctrl_idx = rng.choice(ctrl_idx, size=15000, replace=False)

    use_idx = np.concatenate([carrier_idx, ctrl_idx])
    tmp = obs.iloc[use_idx].copy()
    tmp["target"] = target
    tmp["is_target_carrier"] = carrier[use_idx]
    tmp["is_control_carrier"] = ctrl[use_idx]
    hippo_cell_rows.append(tmp)

if hippo_cell_rows:
    hippo_cells = pd.concat(hippo_cell_rows, axis=0, ignore_index=True)
    keep_cols = [
        "cell_id", DAY_COL, CELLTYPE_COL, "day_lower", "cell_type_str",
        "d1_program_state", "d2_state", "d2_program_state", "d3_state",
        "n_total_guides", "n_total_targets", "n_selected_targets",
        "target", "is_target_carrier", "is_control_carrier",
    ]
    keep_cols += [c for c in obs.columns if c.startswith("Program ")]
    keep_cols += list(FATE_AXIS_DEFS.keys())
    keep_cols = [c for c in keep_cols if c in hippo_cells.columns]
    hippo_cells[keep_cols].to_csv(os.path.join(OUTDIR, "hippo_carrier_cell_scores_long.csv"), index=False)

# UMAP/embedding export if available.
embed_key = None
for k in ["X_umap", "X_tsne", "umap", "tsne"]:
    if k in adata.obsm:
        embed_key = k
        break

if embed_key is not None and hippo_present:
    emb = np.asarray(adata.obsm[embed_key])
    embed_df = pd.DataFrame({"embed1": emb[:, 0], "embed2": emb[:, 1]}, index=adata.obs_names)
    umap_rows = []
    for target in hippo_present:
        carrier = target_assignment[target].values
        for day in [DAY_D1, DAY_D2, DAY_D3]:
            day_m = obs["day_lower"].eq(day).values
            carrier_idx = np.where(carrier & day_m)[0]
            bg_idx = np.where((~carrier) & day_m)[0]
            if len(bg_idx) > 5000:
                bg_idx = rng.choice(bg_idx, size=5000, replace=False)
            use_idx = np.concatenate([carrier_idx, bg_idx])
            tmp = obs.iloc[use_idx][["cell_id", DAY_COL, CELLTYPE_COL, "day_lower", "cell_type_str"]].copy()
            tmp["target"] = target
            tmp["is_target_carrier"] = carrier[use_idx]
            tmp["embed1"] = emb[use_idx, 0]
            tmp["embed2"] = emb[use_idx, 1]
            tmp["embedding_key"] = embed_key
            umap_rows.append(tmp)
    pd.concat(umap_rows, ignore_index=True).to_csv(
        os.path.join(OUTDIR, "hippo_umap_overlay_long.csv"),
        index=False
    )

# Hippo guide convergence summary.
if len(guide_axis):
    hippo_guide = guide_axis[guide_axis["is_hippo_target"]].copy()
    conv_rows = []
    for mode_name in sorted(hippo_guide["mode"].unique()):
        for day in sorted(hippo_guide["day"].unique()):
            for score in sorted(hippo_guide["score"].unique()):
                vals = hippo_guide[
                    hippo_guide["mode"].eq(mode_name)
                    & hippo_guide["day"].eq(day)
                    & hippo_guide["score"].eq(score)
                ]["effect_vs_control"].dropna()
                if len(vals) == 0:
                    continue
                conv_rows.append({
                    "mode": mode_name,
                    "day": day,
                    "score": score,
                    "n_guides": len(vals),
                    "mean_effect": vals.mean(),
                    "median_effect": vals.median(),
                    "fraction_positive": (vals > 0).mean(),
                    "fraction_negative": (vals < 0).mean(),
                    "one_sample_t_p": ttest_1samp(vals, 0).pvalue if len(vals) >= 2 else np.nan,
                    "signed_rank_p": wilcoxon(vals).pvalue if len(vals) >= 2 and not np.all(vals == 0) else np.nan,
                })
    pd.DataFrame(conv_rows).to_csv(os.path.join(OUTDIR, "hippo_guide_convergence_summary.csv"), index=False)

# Manifest and README.
pd.DataFrame({
    "file": sorted(os.listdir(OUTDIR)),
    "path": [os.path.abspath(os.path.join(OUTDIR, fn)) for fn in sorted(os.listdir(OUTDIR))]
}).to_csv(os.path.join(OUTDIR, "export_manifest.csv"), index=False)

with open(os.path.join(OUTDIR, "README.txt"), "w") as f:
    f.write(
        "Adjusted Hippo cell-level exports based on user's adata storage conventions.\\n\\n"
        "Key conventions used:\\n"
        "- adata.obsm['guide_assignment'] can be DataFrame/sparse/array.\\n"
        "- guide names come from guide_assignment columns or adata.uns['guide_names'].\\n"
        "- targets are parsed from guide names before '_' and non-targeting guides map to 'non-targeting'.\\n"
        "- adata.obsm['cnmf_usage'] can be DataFrame/sparse/array. Array columns default to Program 0..n-1.\\n"
        "- z-scores are computed relative to non-targeting carriers within day when enough controls exist; otherwise all cells that day.\\n\\n"
        "Positive fate-axis orientation:\\n"
        "- D1_alt_P13_minus_P30: high = P13 over P30.\\n"
        "- D2_alt_P17_minus_P15: high = P17 over P15.\\n"
        "- D2_YAPLat_minus_PreEC: high = mean(P17,P34) over mean(P15,P23).\\n"
        "- D3_ALT_minus_EC: high = mean(P27,P34,P35) over mean(P2,P3,P8,P11).\\n\\n"
        "Most important outputs:\\n"
        "- d1_target_P13high_vs_P30high_*.csv\\n"
        "- d2_target_preEC_vs_latMeso_*.csv\\n"
        "- d3_target_alt_vs_endo_*.csv\\n"
        "- d1_d2_d3_target_leadlag_all_modes.csv\\n"
        "- target_level_fate_axis_distribution_tests.csv\\n"
        "- guide_level_fate_axis_effects.csv\\n"
        "- hippo_carrier_cell_scores_long.csv\\n"
        "- hippo_umap_overlay_long.csv, if embedding exists.\\n"
        "- hippo_guide_convergence_summary.csv\\n"
    )

print("Done. Wrote outputs to:", OUTDIR)
print("Hippo targets present:", hippo_present)

