
"""
Cell-level D2 TGF/BMP-SMAD diagnostics for high-MOI CC-Perturb-seq.
Run this with `adata` loaded. This exports CSVs that can be plotted in ggplot2.

Interpretation caveat: because MOI is high, target-carrier cells are not pure single-perturbation cells.
Use these as carrier-enrichment analyses, not single-cell fate assignment.
"""
import os
import numpy as np
import pandas as pd
from scipy.sparse import issparse

DAY_COL = "day"
CELLTYPE_COL = "cell_type"
DAY_D2 = "d2"
DAY_D3 = "d3"
PROGRAM_SCORES_KEY = "cnmf_usage"
CONTROL_TARGET = "non-targeting"
TGF_TARGETS = ["ACVR1B", "ACVR2A", "ACVR2B", "BMPR1A", "SMAD1", "SMAD3", "SMAD4", "SMAD5", "TGFBR2"]
TARGETS_TO_EXPORT = [CONTROL_TARGET] + TGF_TARGETS
OUTDIR = "./tgfb_d2_celllevel_csv_exports"
os.makedirs(OUTDIR, exist_ok=True)

def prog(i): return f"Program {i}"
D2_PROGRAMS = [prog(14), prog(15), prog(17), prog(19), prog(23), prog(34)]
D2_EC_PROGRAMS = [prog(15), prog(23)]
D2_FAILURE_PROGRAMS = [prog(14), prog(15), prog(23)]
D3_EC_PROGRAMS = [prog(2), prog(3), prog(8), prog(11)]
D3_ALT_PROGRAMS = [prog(27), prog(34), prog(35)]

def get_guide_names(adata):
    ga = adata.obsm["guide_assignment"]
    if isinstance(ga, pd.DataFrame): return np.array(ga.columns.astype(str))
    if "guide_names" in adata.uns: return np.array(adata.uns["guide_names"]).astype(str)
    raise KeyError("Expected adata.uns['guide_names'] or DataFrame guide_assignment columns.")

def parse_guide_target(g):
    g = str(g); gl = g.lower()
    if "non-target" in gl or gl.startswith("ntc"): return "non-targeting"
    if "_" in g: return g.split("_")[0]
    for sep in ["-", ".", " "]:
        if sep in g: return g.split(sep)[0]
    return g

def get_guide_targets(adata):
    return np.array([parse_guide_target(g) for g in get_guide_names(adata)], dtype=str)

def get_guide_assignment_df(adata):
    ga = adata.obsm["guide_assignment"]
    guide_names = get_guide_names(adata)
    if isinstance(ga, pd.DataFrame): return ga.copy()
    if issparse(ga): ga = ga.toarray()
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
        df = X.copy(); df.index = adata.obs_names
        new_cols = []
        for c in df.columns:
            c = str(c)
            if c.startswith("Program "): new_cols.append(c)
            elif c.startswith("P"): new_cols.append(f"Program {int(c.replace('P',''))}")
            else: new_cols.append(f"Program {int(c)}")
        df.columns = new_cols
        return df
    if issparse(X): X = X.toarray()
    return pd.DataFrame(X, index=adata.obs_names, columns=[f"Program {i}" for i in range(X.shape[1])])

def day_mask(adata, day):
    return adata.obs[DAY_COL].astype(str).str.lower().values == str(day).lower()

def zscore_within_mask(score_df, mask, programs):
    ref = score_df.loc[score_df.index[mask], programs]
    mu = ref.mean(axis=0)
    sd = ref.std(axis=0, ddof=0).replace(0, 1)
    return (score_df[programs] - mu) / sd

def sem(x):
    x = pd.Series(x).dropna()
    if len(x) < 2: return np.nan
    return x.std(ddof=1) / np.sqrt(len(x))

def summarize_group(df, group_cols, value_cols):
    rows = []
    for key, g in df.groupby(group_cols):
        if not isinstance(key, tuple): key = (key,)
        row = dict(zip(group_cols, key)); row["n_cells"] = len(g)
        for v in value_cols:
            row[f"{v}_mean"] = g[v].mean()
            row[f"{v}_median"] = g[v].median()
            row[f"{v}_sem"] = sem(g[v])
        rows.append(row)
    return pd.DataFrame(rows)

# Build matrices
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

# D2 scores
mask_d2 = day_mask(adata, DAY_D2)
d2_z = zscore_within_mask(score_df, mask_d2, D2_PROGRAMS).rename(columns={p:f"z_{p}" for p in D2_PROGRAMS})
d2_raw = score_df[D2_PROGRAMS].rename(columns={p:f"raw_{p}" for p in D2_PROGRAMS})
d2_all = pd.concat([cell_meta, d2_raw, d2_z], axis=1)
d2_all = d2_all.loc[d2_all["day"].str.lower() == DAY_D2].copy()
d2_all["d2_ec_induction_score"] = d2_all[[f"z_{p}" for p in D2_EC_PROGRAMS]].mean(axis=1)
d2_all["d2_p19_minus_p14"] = d2_all[f"z_{prog(19)}"] - d2_all[f"z_{prog(14)}"]
d2_all["d2_p19_failure_score"] = d2_all[f"z_{prog(19)}"] - d2_all[[f"z_{p}" for p in [prog(14), prog(15), prog(23)]]].mean(axis=1)

d2_rows = []
for target in TARGETS_TO_EXPORT:
    if target not in target_assignment.columns: continue
    cells = target_assignment.index[target_assignment[target].values & mask_d2]
    tmp = d2_all.loc[cells].copy(); tmp["target"] = target
    d2_rows.append(tmp)
d2_cell_by_target = pd.concat(d2_rows, ignore_index=True)
d2_cell_by_target.to_csv(os.path.join(OUTDIR, "d2_tgfb_cell_scores_by_target.csv"), index=False)

# long programs
long_rows = []
for _, row in d2_cell_by_target.iterrows():
    for p in D2_PROGRAMS:
        long_rows.append({"cell_id": row["cell_id"], "target": row["target"], "program": p, "raw_score": row[f"raw_{p}"], "z_score": row[f"z_{p}"], "d2_p19_failure_score": row["d2_p19_failure_score"], "n_total_guides": row["n_total_guides"], "n_total_targets": row["n_total_targets"], "n_selected_targets": row["n_selected_targets"]})
pd.DataFrame(long_rows).to_csv(os.path.join(OUTDIR, "d2_tgfb_program_scores_long_by_target.csv"), index=False)

summ = summarize_group(d2_cell_by_target, ["target"], ["d2_ec_induction_score", "d2_p19_minus_p14", "d2_p19_failure_score"] + [f"z_{p}" for p in D2_PROGRAMS])
summ.to_csv(os.path.join(OUTDIR, "d2_tgfb_target_summary.csv"), index=False)

# Binary branch based on P19 failure score > 0
d2_cell_by_target["binary_branch"] = np.where(d2_cell_by_target["d2_p19_failure_score"] > 0, "P19-biased / failed-transition", "EC/lateral-mesoderm biased")
branch = d2_cell_by_target.groupby(["target","binary_branch"]).size().reset_index(name="n_cells")
branch["fraction"] = branch.groupby("target")["n_cells"].transform(lambda x: x / x.sum())
branch.to_csv(os.path.join(OUTDIR, "d2_tgfb_binary_branch_fraction_by_target.csv"), index=False)

# Guide-level concordance
rows = []
for j, guide in enumerate(guide_names):
    target = guide_targets[j]
    if target not in TARGETS_TO_EXPORT: continue
    carrier = (guide_assignment.iloc[:, j].values > 0) & mask_d2
    if carrier.sum() == 0: continue
    sub = d2_all.loc[adata.obs_names[carrier]]
    rows.append({"guide": guide, "target": target, "n_cells": int(carrier.sum()), "d2_p19_failure_score_mean": sub["d2_p19_failure_score"].mean(), "d2_p19_minus_p14_mean": sub["d2_p19_minus_p14"].mean(), "d2_ec_induction_score_mean": sub["d2_ec_induction_score"].mean()})
guide_concordance = pd.DataFrame(rows)
guide_concordance.to_csv(os.path.join(OUTDIR, "d2_tgfb_guide_level_concordance.csv"), index=False)

# D3 endpoint score for lead-lag
mask_d3 = day_mask(adata, DAY_D3)
D3_ALL = D3_EC_PROGRAMS + D3_ALT_PROGRAMS
d3_z = zscore_within_mask(score_df, mask_d3, D3_ALL).rename(columns={p:f"z_{p}" for p in D3_ALL})
d3_all = pd.concat([cell_meta, d3_z], axis=1)
d3_all = d3_all.loc[d3_all["day"].str.lower() == DAY_D3].copy()
d3_all["d3_ec_endpoint_score"] = d3_all[[f"z_{p}" for p in D3_EC_PROGRAMS]].mean(axis=1)
d3_all["d3_alt_endpoint_score"] = d3_all[[f"z_{p}" for p in D3_ALT_PROGRAMS]].mean(axis=1)
d3_all["d3_endpoint_score"] = d3_all["d3_ec_endpoint_score"] - d3_all["d3_alt_endpoint_score"]

d3_rows = []
for target in TARGETS_TO_EXPORT:
    if target not in target_assignment.columns: continue
    cells = target_assignment.index[target_assignment[target].values & mask_d3]
    tmp = d3_all.loc[cells].copy(); tmp["target"] = target
    d3_rows.append(tmp)
d3_cell_by_target = pd.concat(d3_rows, ignore_index=True)
d3_cell_by_target.to_csv(os.path.join(OUTDIR, "d3_tgfb_cell_endpoint_scores_by_target.csv"), index=False)
d3summ = summarize_group(d3_cell_by_target, ["target"], ["d3_endpoint_score", "d3_ec_endpoint_score", "d3_alt_endpoint_score"])
d3summ.to_csv(os.path.join(OUTDIR, "d3_tgfb_target_endpoint_summary.csv"), index=False)
leadlag = summ.merge(d3summ, on="target", how="inner")
leadlag.to_csv(os.path.join(OUTDIR, "d2_tgfb_p19_failure_vs_d3_endpoint_leadlag.csv"), index=False)
print("Done exporting TGF/BMP D2 cell-level CSVs to", OUTDIR)
