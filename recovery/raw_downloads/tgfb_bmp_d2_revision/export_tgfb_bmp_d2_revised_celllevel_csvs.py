
"""
Export cell-level CSVs for D2 TGF/BMP-SMAD analysis.

Designed for AnnData with:
- adata.obs["day"] containing d0,d1,d2,d3
- adata.obs["cell_type"]
- adata.uns["guide_names"]
- adata.obsm["guide_assignment"] as numpy array or sparse matrix
- adata.obsm["cnmf_usage"] with columns corresponding to Program 0..Program 49

This is written for high-MOI experiments. Outputs are target-carrier enrichment summaries,
not pure single-gene perturbation fate assignments.
"""

import os
import numpy as np
import pandas as pd
from scipy.sparse import issparse

# =====================
# CONFIG
# =====================
OUTDIR = "./tgfb_bmp_d2_celllevel_exports"
os.makedirs(OUTDIR, exist_ok=True)

DAY_COL = "day"
CELLTYPE_COL = "cell_type"
PROGRAM_SCORES_KEY = "cnmf_usage"
DAY_D2 = "d2"
DAY_D3 = "d3"

def prog(i):
    return f"Program {i}"

CONTROL_TARGET = "non-targeting"
TARGETS_TO_EXPORT = [
    CONTROL_TARGET,
    "ACVR1B", "TGFBR2", "SMAD3",
    "BMPR1A", "SMAD1", "SMAD5",
    "ACVR2A", "ACVR2B", "SMAD4",
    "ETV2", "KDR", "NF2", "LATS2", "AMOTL2"
]

D2_PROGRAMS = [prog(14), prog(17), prog(19), prog(15), prog(23), prog(34)]
D2_EC_PROGRAMS = [prog(15), prog(23)]
D2_P14P17_PROGRAMS = [prog(14), prog(17)]
D3_EC_PROGRAMS = [prog(2), prog(3), prog(8), prog(11)]
D3_ALT_PROGRAMS = [prog(27), prog(34), prog(35)]


# =====================
# HELPERS
# =====================
def get_guide_names(adata):
    ga = adata.obsm["guide_assignment"]
    if isinstance(ga, pd.DataFrame):
        return np.array(ga.columns.astype(str))
    if "guide_names" in adata.uns:
        return np.array(adata.uns["guide_names"]).astype(str)
    raise KeyError("Expected guide names in adata.uns['guide_names'] or guide_assignment columns.")

def parse_guide_target(guide_name):
    g = str(guide_name)
    gl = g.lower()
    if "non-target" in gl or gl.startswith("ntc"):
        return "non-targeting"
    if "_" in g:
        return g.split("_")[0]
    for sep in ["-", ".", " "]:
        if sep in g:
            return g.split(sep)[0]
    return g

def get_guide_targets(adata):
    return np.array([parse_guide_target(g) for g in get_guide_names(adata)], dtype=str)

def get_guide_assignment_df(adata):
    ga = adata.obsm["guide_assignment"]
    guide_names = get_guide_names(adata)
    if isinstance(ga, pd.DataFrame):
        return ga.copy()
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

def get_program_score_df(adata):
    X = adata.obsm[PROGRAM_SCORES_KEY]
    if isinstance(X, pd.DataFrame):
        df = X.copy()
        new_cols = []
        for c in df.columns:
            c = str(c)
            if c.startswith("Program "):
                new_cols.append(c)
            elif c.startswith("P"):
                new_cols.append(f"Program {int(c.replace('P', ''))}")
            else:
                new_cols.append(f"Program {int(c)}")
        df.columns = new_cols
        df.index = adata.obs_names
        return df
    if issparse(X):
        X = X.toarray()
    return pd.DataFrame(X, index=adata.obs_names, columns=[f"Program {i}" for i in range(X.shape[1])])

def day_mask(adata, day):
    return adata.obs[DAY_COL].astype(str).str.lower().values == str(day).lower()

def zscore_within_mask(score_df, mask, programs):
    ref = score_df.loc[score_df.index[mask], programs].copy()
    mu = ref.mean(axis=0)
    sd = ref.std(axis=0, ddof=0).replace(0, 1)
    return (score_df[programs] - mu) / sd

def safe_sem(x):
    x = pd.Series(x).dropna()
    if len(x) <= 1:
        return np.nan
    return x.std(ddof=1) / np.sqrt(len(x))

def summarize_group(df, group_cols, value_cols):
    rows = []
    for key, g in df.groupby(group_cols):
        if not isinstance(key, tuple):
            key = (key,)
        row = dict(zip(group_cols, key))
        row["n_cells"] = len(g)
        for v in value_cols:
            row[f"{v}_mean"] = g[v].mean()
            row[f"{v}_median"] = g[v].median()
            row[f"{v}_sem"] = safe_sem(g[v])
        rows.append(row)
    return pd.DataFrame(rows)


# =====================
# BASE DATA
# =====================
guide_names = get_guide_names(adata)
guide_targets = get_guide_targets(adata)
guide_assignment = get_guide_assignment_df(adata)
target_assignment = get_target_assignment_df(adata)
score_df = get_program_score_df(adata)

n_total_guides = (guide_assignment.values > 0).sum(axis=1)
n_total_targets = target_assignment.sum(axis=1).values
selected_cols = [t for t in TARGETS_TO_EXPORT if t in target_assignment.columns]
n_selected_targets = target_assignment[selected_cols].sum(axis=1).values

cell_meta = pd.DataFrame({
    "cell_id": adata.obs_names,
    "day": adata.obs[DAY_COL].astype(str).values,
    "cell_type": adata.obs[CELLTYPE_COL].astype(str).values,
    "n_total_guides": n_total_guides,
    "n_total_targets": n_total_targets,
    "n_selected_targets": n_selected_targets,
}, index=adata.obs_names)

cell_meta.to_csv(os.path.join(OUTDIR, "cell_moi_metadata.csv"), index=False)
pd.DataFrame({"guide": guide_names, "target": guide_targets}).to_csv(os.path.join(OUTDIR, "guide_to_target_map.csv"), index=False)
pd.Series(guide_targets, name="guide_target").value_counts().to_csv(os.path.join(OUTDIR, "guide_target_counts.csv"))


# =====================
# D2 SCORES
# =====================
mask_d2 = day_mask(adata, DAY_D2)
d2_z = zscore_within_mask(score_df, mask_d2, D2_PROGRAMS).rename(columns={p: f"z_{p}" for p in D2_PROGRAMS})
d2_raw = score_df[D2_PROGRAMS].rename(columns={p: f"raw_{p}" for p in D2_PROGRAMS})

d2_scores = pd.concat([cell_meta, d2_raw, d2_z], axis=1)
d2_scores = d2_scores.loc[d2_scores["day"].str.lower() == DAY_D2].copy()

d2_scores["D2_EC_induction_score"] = d2_scores[[f"z_{p}" for p in D2_EC_PROGRAMS]].mean(axis=1)
d2_scores["D2_P14P17_mesoderm_score"] = d2_scores[[f"z_{p}" for p in D2_P14P17_PROGRAMS]].mean(axis=1)
d2_scores["D2_P19_score"] = d2_scores[f"z_{prog(19)}"]
d2_scores["D2_P19_vs_EC_score"] = d2_scores["D2_P19_score"] - d2_scores["D2_EC_induction_score"]
d2_scores["D2_P19_vs_P14P17_score"] = d2_scores["D2_P19_score"] - d2_scores["D2_P14P17_mesoderm_score"]

# cell x target rows
target_rows = []
for target in TARGETS_TO_EXPORT:
    if target not in target_assignment.columns:
        print(f"Warning: {target} not present in target assignments")
        continue
    cells = target_assignment.index[target_assignment[target].values & mask_d2]
    tmp = d2_scores.loc[cells].copy()
    tmp["target"] = target
    target_rows.append(tmp)

d2_cell_by_target = pd.concat(target_rows, ignore_index=True)
d2_cell_by_target.to_csv(os.path.join(OUTDIR, "d2_cell_scores_by_target.csv"), index=False)

# long program table
long_rows = []
for _, row in d2_cell_by_target.iterrows():
    for p in D2_PROGRAMS:
        long_rows.append({
            "cell_id": row["cell_id"],
            "target": row["target"],
            "program": p,
            "raw_score": row[f"raw_{p}"],
            "z_score": row[f"z_{p}"],
            "D2_EC_induction_score": row["D2_EC_induction_score"],
            "D2_P14P17_mesoderm_score": row["D2_P14P17_mesoderm_score"],
            "D2_P19_score": row["D2_P19_score"],
            "D2_P19_vs_EC_score": row["D2_P19_vs_EC_score"],
            "D2_P19_vs_P14P17_score": row["D2_P19_vs_P14P17_score"],
            "n_total_guides": row["n_total_guides"],
            "n_total_targets": row["n_total_targets"],
            "n_selected_targets": row["n_selected_targets"],
        })
pd.DataFrame(long_rows).to_csv(os.path.join(OUTDIR, "d2_program_scores_long_by_target.csv"), index=False)

d2_summary = summarize_group(
    d2_cell_by_target,
    ["target"],
    ["D2_EC_induction_score", "D2_P14P17_mesoderm_score", "D2_P19_score", "D2_P19_vs_EC_score", "D2_P19_vs_P14P17_score"]
)
d2_summary.to_csv(os.path.join(OUTDIR, "d2_target_scores_summary.csv"), index=False)

# branch fractions
d2_cell_by_target["D2_binary_branch"] = np.where(
    d2_cell_by_target["D2_P19_vs_EC_score"] > 0,
    "P19-high / EC-induction-low",
    "EC-induction-biased"
)
branch_frac = d2_cell_by_target.groupby(["target","D2_binary_branch"]).size().reset_index(name="n_cells")
branch_frac["fraction"] = branch_frac.groupby("target")["n_cells"].transform(lambda x: x / x.sum())
branch_frac.to_csv(os.path.join(OUTDIR, "d2_p19_vs_ec_branch_fraction_by_target.csv"), index=False)

d2_cell_by_target["D2_program_axis"] = d2_cell_by_target[[f"z_{p}" for p in D2_PROGRAMS]].idxmax(axis=1).str.replace("z_", "", regex=False)
prog_frac = d2_cell_by_target.groupby(["target","D2_program_axis"]).size().reset_index(name="n_cells")
prog_frac["fraction"] = prog_frac.groupby("target")["n_cells"].transform(lambda x: x / x.sum())
prog_frac.to_csv(os.path.join(OUTDIR, "d2_dominant_program_axis_fraction_by_target.csv"), index=False)


# =====================
# GUIDE-LEVEL D2 CONCORDANCE
# =====================
guide_rows = []
for j, guide in enumerate(guide_names):
    target = guide_targets[j]
    if target not in TARGETS_TO_EXPORT:
        continue
    carrier = (guide_assignment.iloc[:, j].values > 0) & mask_d2
    if carrier.sum() == 0:
        continue
    sub = d2_scores.loc[adata.obs_names[carrier]]
    guide_rows.append({
        "guide": guide,
        "target": target,
        "n_cells_d2": int(carrier.sum()),
        "D2_EC_induction_score_mean": sub["D2_EC_induction_score"].mean(),
        "D2_P14P17_mesoderm_score_mean": sub["D2_P14P17_mesoderm_score"].mean(),
        "D2_P19_score_mean": sub["D2_P19_score"].mean(),
        "D2_P19_vs_EC_score_mean": sub["D2_P19_vs_EC_score"].mean(),
        "D2_P19_vs_P14P17_score_mean": sub["D2_P19_vs_P14P17_score"].mean(),
    })
guide_d2 = pd.DataFrame(guide_rows)
guide_d2.to_csv(os.path.join(OUTDIR, "d2_guide_level_concordance.csv"), index=False)


# =====================
# D3 ENDPOINT SCORES AND LEAD-LAG
# =====================
mask_d3 = day_mask(adata, DAY_D3)
d3_programs = D3_EC_PROGRAMS + D3_ALT_PROGRAMS
d3_z = zscore_within_mask(score_df, mask_d3, d3_programs).rename(columns={p: f"z_{p}" for p in d3_programs})
d3_raw = score_df[d3_programs].rename(columns={p: f"raw_{p}" for p in d3_programs})
d3_scores = pd.concat([cell_meta, d3_raw, d3_z], axis=1)
d3_scores = d3_scores.loc[d3_scores["day"].str.lower() == DAY_D3].copy()
d3_scores["D3_EC_endpoint_score"] = d3_scores[[f"z_{p}" for p in D3_EC_PROGRAMS]].mean(axis=1)
d3_scores["D3_ALT_endpoint_score"] = d3_scores[[f"z_{p}" for p in D3_ALT_PROGRAMS]].mean(axis=1)
d3_scores["D3_endpoint_score"] = d3_scores["D3_EC_endpoint_score"] - d3_scores["D3_ALT_endpoint_score"]

d3_target_rows = []
for target in TARGETS_TO_EXPORT:
    if target not in target_assignment.columns:
        continue
    cells = target_assignment.index[target_assignment[target].values & mask_d3]
    tmp = d3_scores.loc[cells].copy()
    tmp["target"] = target
    d3_target_rows.append(tmp)
d3_cell_by_target = pd.concat(d3_target_rows, ignore_index=True)
d3_cell_by_target.to_csv(os.path.join(OUTDIR, "d3_cell_endpoint_scores_by_target.csv"), index=False)

d3_summary = summarize_group(
    d3_cell_by_target,
    ["target"],
    ["D3_EC_endpoint_score", "D3_ALT_endpoint_score", "D3_endpoint_score"]
)
d3_summary.to_csv(os.path.join(OUTDIR, "d3_target_endpoint_summary.csv"), index=False)

leadlag = d2_summary.merge(d3_summary, on="target", how="inner")
leadlag.to_csv(os.path.join(OUTDIR, "d2_tgfb_scores_vs_d3_endpoint_leadlag.csv"), index=False)

guide_d3_rows = []
for j, guide in enumerate(guide_names):
    target = guide_targets[j]
    if target not in TARGETS_TO_EXPORT:
        continue
    carrier = (guide_assignment.iloc[:, j].values > 0) & mask_d3
    if carrier.sum() == 0:
        continue
    sub = d3_scores.loc[adata.obs_names[carrier]]
    guide_d3_rows.append({
        "guide": guide,
        "target": target,
        "n_cells_d3": int(carrier.sum()),
        "D3_EC_endpoint_score_mean": sub["D3_EC_endpoint_score"].mean(),
        "D3_ALT_endpoint_score_mean": sub["D3_ALT_endpoint_score"].mean(),
        "D3_endpoint_score_mean": sub["D3_endpoint_score"].mean(),
    })
guide_d3 = pd.DataFrame(guide_d3_rows)
guide_leadlag = guide_d2.merge(guide_d3, on=["guide","target"], how="inner")
guide_leadlag.to_csv(os.path.join(OUTDIR, "guide_level_d2_tgfb_scores_vs_d3_endpoint.csv"), index=False)

print(f"Export complete: {OUTDIR}")
for f in sorted(os.listdir(OUTDIR)):
    print(" -", f)
