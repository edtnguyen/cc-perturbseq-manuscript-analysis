
#!/usr/bin/env python
"""
export_lineage_split_celllevel_csvs.py

Run this in the Python environment where your adata is available.

Notebook use:
    from export_lineage_split_celllevel_csvs import export_lineage_split_celllevel_csvs
    export_lineage_split_celllevel_csvs(
        adata,
        guide_target_csv="combined_TargetKD.csv",
        outdir="lineage_split_celllevel_exports",
        selected_targets=None,  # None = all targets
        highlight_targets=["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"],
        export_neighborhoods=True,
    )

Command-line use:
    python export_lineage_split_celllevel_csvs.py \
        --adata your_adata.h5ad \
        --guide-target-csv combined_TargetKD.csv \
        --outdir lineage_split_celllevel_exports \
        --export-neighborhoods

What this exports:
    1. target/guide cell-type composition and one-vs-rest enrichment
    2. D1/D2/D3 program-axis effects and distribution tests
    3. D1->D2->D3 lead-lag target tables
    4. guide-level concordance summaries
    5. whole-composition Jensen-Shannon distances from non-targeting controls
    6. sample/replicate pseudobulk count tables when replicate columns exist
    7. optional neighborhood enrichment tables as a Milo-like local DA input

Important:
    This is designed for high-MOI screens. These are target-carrier enrichment
    analyses, not pure single-guide single-perturbation fate assignments.
"""

import os
import argparse
import numpy as np
import pandas as pd
from scipy import sparse
from scipy import stats
from scipy.spatial.distance import jensenshannon
from sklearn.neighbors import NearestNeighbors


# ============================================================
# Defaults
# ============================================================
DEFAULT_HIGHLIGHT_TARGETS = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]

DEFAULT_D1_CELL_TYPES = ["PrimStreak", "Lat-Meso", "Neu-Mesen"]
DEFAULT_D2_PRE_EC = ["Pre-EC"]
DEFAULT_D2_LAT_MESO = ["Lat-Meso"]
DEFAULT_D3_EC = ["artEC", "ECM-EC", "Endothelial", "EC", "Angio-EC", "Partial EndMT EC"]
DEFAULT_D3_ALT = ["Fibro-Mesen", "Neu-Mesen", "Mesenchymal", "NeuCrest/Mes"]
DEFAULT_D3_MES = ["Fibro-Mesen", "Mesenchymal"]

PROGRAMS = sorted(set([2, 3, 8, 11, 13, 14, 15, 17, 19, 23, 27, 30, 34, 35]))

AXIS_DEFS = {
    "D1_P13_minus_P30": {
        "day": "d1",
        "plus": [13],
        "minus": [30],
        "description": "D1 program split, positive = P13 over P30",
    },
    "D2_P17_minus_P15": {
        "day": "d2",
        "plus": [17],
        "minus": [15],
        "description": "D2 program split, positive = P17 over P15",
    },
    "D2_YAPLat_minus_PreEC": {
        "day": "d2",
        "plus": [17, 34],
        "minus": [15, 23],
        "description": "D2 bifurcation, positive = mean(P17,P34) over mean(P15,P23)",
    },
    "D3_ALT_minus_EC": {
        "day": "d3",
        "plus": [27, 34, 35],
        "minus": [2, 3, 8, 11],
        "description": "D3 endpoint, positive = mean(P27,P34,P35) over mean(P2,P3,P8,P11)",
    },
}


# ============================================================
# Generic utilities
# ============================================================
def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = np.isfinite(p)
    if mask.sum() == 0:
        return out
    x = p[mask]
    order = np.argsort(x)
    ranks = np.arange(1, len(x) + 1)
    q = np.empty_like(x)
    q[order] = np.minimum.accumulate((x[order][::-1] * len(x) / ranks[::-1]))[::-1]
    out[mask] = np.clip(q, 0, 1)
    return out


def safe_div(a, b):
    return a / b if b else np.nan


def log2_or(a, b, c, d, pseudocount=0.5):
    return np.log2(((a + pseudocount) / (b + pseudocount)) / ((c + pseudocount) / (d + pseudocount)))


def cliffs_delta(x, y, max_pairs=2_000_000, random_seed=1):
    x = np.asarray(pd.Series(x).dropna(), dtype=float)
    y = np.asarray(pd.Series(y).dropna(), dtype=float)
    if len(x) == 0 or len(y) == 0:
        return np.nan
    if len(x) * len(y) > max_pairs:
        rng = np.random.default_rng(random_seed)
        x = rng.choice(x, size=min(len(x), 1500), replace=False)
        y = rng.choice(y, size=min(len(y), 1500), replace=False)
    gt = sum((xi > y).sum() for xi in x)
    lt = sum((xi < y).sum() for xi in x)
    return (gt - lt) / (len(x) * len(y))


def clean_day(x):
    return str(x).strip().lower()


def looks_nt(x):
    s = str(x).lower()
    patterns = [
        "non-target", "nontarget", "non_target", "negative", "control",
        "scramble", "scrambled", "safe", "ntc", "non targeting", "non-targeting"
    ]
    return any(p in s for p in patterns) or s in {"nt", "neg", "ctrl"}


def infer_col(df, candidates):
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    for col in df.columns:
        cl = col.lower()
        for c in candidates:
            if c.lower() in cl:
                return col
    return None


def parse_target_from_guide(guide):
    g = str(guide)
    if looks_nt(g):
        return "non-targeting"
    for sep in ["|", ":", "_", "-", "."]:
        if sep in g:
            return g.split(sep)[0]
    return g


def as_bool_dense(x):
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x) > 0


# ============================================================
# adata extraction helpers
# ============================================================
def get_guide_names(adata, guide_assignment_key="guide_assignment"):
    ga = adata.obsm[guide_assignment_key]
    if isinstance(ga, pd.DataFrame):
        return np.array(ga.columns.astype(str))
    if "guide_names" in adata.uns:
        return np.array(adata.uns["guide_names"]).astype(str)
    raise KeyError("Could not find guide names. Expected adata.uns['guide_names'] or DataFrame columns in adata.obsm['guide_assignment'].")


def get_guide_assignment_df(adata, guide_assignment_key="guide_assignment"):
    ga = adata.obsm[guide_assignment_key]
    guide_names = get_guide_names(adata, guide_assignment_key)
    if isinstance(ga, pd.DataFrame):
        out = ga.copy()
        out.index = adata.obs_names
        out.columns = out.columns.astype(str)
        return out.astype(bool)
    ga = as_bool_dense(ga)
    return pd.DataFrame(ga, index=adata.obs_names, columns=guide_names)


def build_guide_target_table(adata, guide_target_csv=None, guide_assignment_key="guide_assignment"):
    guide_names = pd.Index(get_guide_names(adata, guide_assignment_key), name="guide")
    base = pd.DataFrame({"guide": guide_names.astype(str)})

    if guide_target_csv and os.path.exists(guide_target_csv):
        gt = pd.read_csv(guide_target_csv)
        guide_col = infer_col(gt, ["guide", "guide_name", "sgRNA", "sgrna", "gRNA", "oligo", "spacer"])
        target_col = infer_col(gt, ["target", "target_name", "target_gene", "gene", "symbol", "Target", "TargetGene"])
        if guide_col is None or target_col is None:
            raise ValueError(f"Could not infer guide/target columns in {guide_target_csv}. Columns: {list(gt.columns)}")
        gt = gt[[guide_col, target_col]].copy()
        gt.columns = ["guide", "target"]
        gt["guide"] = gt["guide"].astype(str)
        gt["target"] = gt["target"].astype(str)
        out = base.merge(gt, on="guide", how="left")
        out["target"] = out["target"].fillna(out["guide"].map(parse_target_from_guide))
    else:
        out = base.copy()
        out["target"] = out["guide"].map(parse_target_from_guide)

    out["is_non_targeting"] = out["guide"].map(looks_nt) | out["target"].map(looks_nt)
    out.loc[out["is_non_targeting"], "target"] = "non-targeting"
    return out


def get_target_assignment_df(guide_assignment, guide_target):
    out = pd.DataFrame(index=guide_assignment.index)
    guide_to_target = dict(zip(guide_target["guide"], guide_target["target"]))
    targets = sorted(pd.unique(guide_target["target"].astype(str)))
    for target in targets:
        guides = guide_target.loc[guide_target["target"].astype(str).eq(target), "guide"].astype(str).tolist()
        guides = [g for g in guides if g in guide_assignment.columns]
        if guides:
            out[target] = guide_assignment[guides].any(axis=1).values
        else:
            out[target] = False
    return out.astype(bool)


def standardize_program_col(c):
    s = str(c)
    if s.startswith("Program "):
        return f"Program {int(s.replace('Program ', ''))}"
    if s.startswith("P"):
        return f"Program {int(s.replace('P', ''))}"
    return f"Program {int(float(s))}"


def get_program_score_df(adata, cnmf_key="cnmf_usage", array_program_index_base=0):
    X = adata.obsm[cnmf_key]
    if isinstance(X, pd.DataFrame):
        df = X.copy()
        df.index = adata.obs_names
        df.columns = [standardize_program_col(c) for c in df.columns]
        return df.astype(float)

    if sparse.issparse(X):
        X = X.toarray()
    X = np.asarray(X)
    cols = [f"Program {i + array_program_index_base}" for i in range(X.shape[1])]
    return pd.DataFrame(X, index=adata.obs_names, columns=cols).astype(float)


def zscore_programs_by_day(program_df, obs, day_col="day_clean", nt_mask=None, min_ref_cells=50):
    z = pd.DataFrame(index=program_df.index)
    refs = []
    for day in sorted(obs[day_col].dropna().unique()):
        day_mask = obs[day_col].eq(day)
        ref_mask = day_mask & nt_mask if nt_mask is not None else pd.Series(False, index=program_df.index)
        if ref_mask.sum() < min_ref_cells:
            ref_mask = day_mask
            ref_name = "all cells at day"
        else:
            ref_name = "non-targeting carriers"
        mu = program_df.loc[ref_mask].mean(axis=0)
        sd = program_df.loc[ref_mask].std(axis=0, ddof=0).replace(0, 1.0)
        z.loc[day_mask, program_df.columns] = (program_df.loc[day_mask] - mu) / sd
        refs.append({"day": day, "reference": ref_name, "n_reference": int(ref_mask.sum())})
    return z.astype(float), pd.DataFrame(refs)


def add_axis_scores(cell_df):
    for axis, spec in AXIS_DEFS.items():
        plus = [f"Program {p}" for p in spec["plus"]]
        minus = [f"Program {p}" for p in spec["minus"]]
        cell_df[axis] = cell_df[plus].mean(axis=1) - cell_df[minus].mean(axis=1)
    cell_df["D1_program_state"] = np.where(cell_df["Program 13"] > cell_df["Program 30"], "P13_high", "P30_high")
    cell_df["D2_program_state"] = np.where(cell_df["Program 17"] > cell_df["Program 15"], "P17_high", "P15_high")
    return cell_df


def summarize_distribution(x, y):
    x = pd.Series(x).dropna().astype(float).values
    y = pd.Series(y).dropna().astype(float).values
    out = {
        "n_target": len(x),
        "n_control": len(y),
        "target_mean": np.mean(x) if len(x) else np.nan,
        "control_mean": np.mean(y) if len(y) else np.nan,
        "delta_mean": np.mean(x) - np.mean(y) if len(x) and len(y) else np.nan,
        "target_median": np.median(x) if len(x) else np.nan,
        "control_median": np.median(y) if len(y) else np.nan,
        "delta_median": np.median(x) - np.median(y) if len(x) and len(y) else np.nan,
        "wasserstein": stats.wasserstein_distance(x, y) if len(x) and len(y) else np.nan,
        "cliffs_delta": cliffs_delta(x, y) if len(x) and len(y) else np.nan,
        "ks_p": stats.ks_2samp(x, y).pvalue if len(x) >= 5 and len(y) >= 5 else np.nan,
        "mannwhitney_p": stats.mannwhitneyu(x, y, alternative="two-sided").pvalue if len(x) >= 5 and len(y) >= 5 else np.nan,
    }
    return out


# ============================================================
# Main export function
# ============================================================
def export_lineage_split_celllevel_csvs(
    adata,
    guide_target_csv=None,
    outdir="lineage_split_celllevel_exports",
    selected_targets=None,
    highlight_targets=None,
    day_col="day",
    celltype_col="cell_type",
    cnmf_key="cnmf_usage",
    guide_assignment_key="guide_assignment",
    array_program_index_base=0,
    min_target_cells=20,
    min_guide_cells=8,
    min_ref_cells=50,
    max_cells_per_group_for_plotting=1000,
    export_cell_table=False,
    export_neighborhoods=False,
    neighborhood_k=80,
    neighborhood_n_seeds=800,
    random_seed=1,
):
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(random_seed)
    if highlight_targets is None:
        highlight_targets = DEFAULT_HIGHLIGHT_TARGETS

    obs = adata.obs.copy()
    obs.index = adata.obs_names
    obs["cell_id"] = adata.obs_names
    obs["day_clean"] = obs[day_col].map(clean_day)
    obs["cell_type_clean"] = obs[celltype_col].astype(str)

    guide_assignment = get_guide_assignment_df(adata, guide_assignment_key=guide_assignment_key)
    guide_target = build_guide_target_table(adata, guide_target_csv=guide_target_csv, guide_assignment_key=guide_assignment_key)
    target_assignment = get_target_assignment_df(guide_assignment, guide_target)

    nt_mask = target_assignment["non-targeting"] if "non-targeting" in target_assignment.columns else pd.Series(False, index=obs.index)
    all_targets = [t for t in target_assignment.columns if t != "non-targeting"]
    if selected_targets is None:
        selected_targets = all_targets
    else:
        selected_targets = [str(t) for t in selected_targets if str(t) in target_assignment.columns and str(t) != "non-targeting"]
    highlight_targets = [t for t in highlight_targets if t in target_assignment.columns]

    target_load_selected = target_assignment[selected_targets].sum(axis=1) if selected_targets else pd.Series(0, index=obs.index)
    obs["n_total_guides"] = guide_assignment.sum(axis=1).astype(int)
    obs["n_total_targets"] = target_assignment.sum(axis=1).astype(int)
    obs["n_selected_targets"] = target_load_selected.astype(int)
    obs["is_non_targeting_carrier"] = nt_mask.values

    program_raw = get_program_score_df(adata, cnmf_key=cnmf_key, array_program_index_base=array_program_index_base)
    needed_cols = [f"Program {p}" for p in PROGRAMS]
    missing = [c for c in needed_cols if c not in program_raw.columns]
    if missing:
        raise KeyError(
            f"Missing required cNMF program columns: {missing}. "
            f"First columns: {list(program_raw.columns[:8])}; last columns: {list(program_raw.columns[-8:])}. "
            f"Try array_program_index_base=1 if indexing is off."
        )

    program_z, ref_summary = zscore_programs_by_day(
        program_raw[needed_cols],
        obs,
        day_col="day_clean",
        nt_mask=nt_mask,
        min_ref_cells=min_ref_cells,
    )

    cell_df = pd.concat([obs, program_z], axis=1)
    cell_df = add_axis_scores(cell_df)

    guide_target.to_csv(os.path.join(outdir, "guide_to_target_map.csv"), index=False)
    ref_summary.to_csv(os.path.join(outdir, "program_zscore_reference_summary.csv"), index=False)

    target_counts = []
    for target in target_assignment.columns:
        carrier = target_assignment[target]
        row = {
            "target": target,
            "n_cells_all_days": int(carrier.sum()),
            "n_guides": int((guide_target["target"] == target).sum()),
            "is_highlight_target": target in highlight_targets,
        }
        for d in sorted(cell_df["day_clean"].unique()):
            row[f"n_cells_{d}"] = int((carrier & cell_df["day_clean"].eq(d)).sum())
        target_counts.append(row)
    pd.DataFrame(target_counts).to_csv(os.path.join(outdir, "target_carrier_counts.csv"), index=False)

    # Modes for high-MOI sensitivity.
    modes = {
        "all_carriers": {
            "target_selector": pd.Series(True, index=cell_df.index),
            "control_selector": pd.Series(True, index=cell_df.index),
        },
        "single_selected_target": {
            "target_selector": cell_df["n_selected_targets"].eq(1),
            "control_selector": cell_df["n_selected_targets"].eq(0),
        },
        "low_guide_load_leq2": {
            "target_selector": cell_df["n_total_guides"].le(2),
            "control_selector": cell_df["n_total_guides"].le(2),
        },
    }

    # ----------------------
    # Pseudobulk sample × target × day × cell type counts
    # ----------------------
    sample_cols = [c for c in ["sample", "Sample", "replicate", "Replicate", "batch", "Batch", "library", "donor", "experiment", "well"] if c in cell_df.columns]
    if sample_cols:
        sample_col = sample_cols[0]
        rows = []
        for mode, masks in modes.items():
            for target in selected_targets + ["non-targeting"]:
                if target not in target_assignment.columns:
                    continue
                carrier = target_assignment[target] & masks["target_selector"] if target != "non-targeting" else nt_mask & masks["control_selector"]
                sub = cell_df.loc[carrier, [sample_col, "day_clean", "cell_type_clean"]].copy()
                if len(sub) == 0:
                    continue
                tab = sub.groupby([sample_col, "day_clean", "cell_type_clean"]).size().reset_index(name="n_cells")
                tab["mode"] = mode
                tab["target"] = target
                rows.append(tab)
        if rows:
            pd.concat(rows, ignore_index=True).to_csv(os.path.join(outdir, "pseudobulk_counts_by_sample_target_day_celltype.csv"), index=False)

    # ----------------------
    # Cell-type composition and one-vs-rest enrichment
    # ----------------------
    celltype_rows = []
    guide_celltype_rows = []
    composition_distance_rows = []

    all_cell_types = sorted(cell_df["cell_type_clean"].dropna().astype(str).unique())
    all_days = sorted(cell_df["day_clean"].dropna().unique())

    for mode, masks in modes.items():
        ctrl_base = nt_mask & masks["control_selector"]

        for target in selected_targets:
            carrier_base = target_assignment[target] & masks["target_selector"]
            for day in all_days:
                day_mask = cell_df["day_clean"].eq(day)
                target_cells = cell_df.loc[carrier_base & day_mask]
                ctrl_cells = cell_df.loc[ctrl_base & day_mask]
                if len(target_cells) < min_target_cells or len(ctrl_cells) < min_ref_cells:
                    continue

                target_vc = target_cells["cell_type_clean"].value_counts()
                ctrl_vc = ctrl_cells["cell_type_clean"].value_counts()
                target_total = int(target_vc.sum())
                ctrl_total = int(ctrl_vc.sum())

                p_vec = np.array([target_vc.get(ct, 0) for ct in all_cell_types], dtype=float) + 1e-8
                q_vec = np.array([ctrl_vc.get(ct, 0) for ct in all_cell_types], dtype=float) + 1e-8
                p_vec = p_vec / p_vec.sum()
                q_vec = q_vec / q_vec.sum()

                composition_distance_rows.append({
                    "mode": mode,
                    "target": target,
                    "day": day,
                    "jensen_shannon_distance": float(jensenshannon(p_vec, q_vec, base=2.0)),
                    "total_variation_distance": float(np.abs(p_vec - q_vec).sum() / 2.0),
                    "n_target_cells": len(target_cells),
                    "n_control_cells": len(ctrl_cells),
                    "is_highlight_target": target in highlight_targets,
                })

                for ct in all_cell_types:
                    a = int(target_vc.get(ct, 0))
                    b = target_total - a
                    c = int(ctrl_vc.get(ct, 0))
                    d = ctrl_total - c
                    if (a + c) == 0:
                        lor = np.nan
                        fp = np.nan
                    else:
                        lor = log2_or(a, b, c, d)
                        fp = stats.fisher_exact([[a, b], [c, d]])[1]
                    celltype_rows.append({
                        "mode": mode,
                        "target": target,
                        "day": day,
                        "cell_type": ct,
                        "target_n": a,
                        "target_total": target_total,
                        "control_n": c,
                        "control_total": ctrl_total,
                        "target_fraction": safe_div(a, target_total),
                        "control_fraction": safe_div(c, ctrl_total),
                        "delta_fraction": safe_div(a, target_total) - safe_div(c, ctrl_total),
                        "log2OR_target_vs_control": lor,
                        "fisher_p": fp,
                        "is_highlight_target": target in highlight_targets,
                    })

        # guide-level cell-type one-vs-rest for selected targets
        for j, guide in enumerate(guide_assignment.columns):
            target = guide_target.loc[guide_target["guide"].astype(str).eq(str(guide)), "target"]
            if len(target) == 0:
                continue
            target = str(target.iloc[0])
            if target not in selected_targets:
                continue
            carrier_base = guide_assignment.iloc[:, j] & masks["target_selector"]
            for day in all_days:
                day_mask = cell_df["day_clean"].eq(day)
                guide_cells = cell_df.loc[carrier_base & day_mask]
                ctrl_cells = cell_df.loc[ctrl_base & day_mask]
                if len(guide_cells) < min_guide_cells or len(ctrl_cells) < min_ref_cells:
                    continue
                guide_vc = guide_cells["cell_type_clean"].value_counts()
                ctrl_vc = ctrl_cells["cell_type_clean"].value_counts()
                guide_total = int(guide_vc.sum())
                ctrl_total = int(ctrl_vc.sum())
                for ct in all_cell_types:
                    a = int(guide_vc.get(ct, 0))
                    b = guide_total - a
                    c = int(ctrl_vc.get(ct, 0))
                    d = ctrl_total - c
                    if (a + c) == 0:
                        continue
                    guide_celltype_rows.append({
                        "mode": mode,
                        "guide": guide,
                        "target": target,
                        "day": day,
                        "cell_type": ct,
                        "guide_n": a,
                        "guide_total": guide_total,
                        "control_n": c,
                        "control_total": ctrl_total,
                        "guide_fraction": safe_div(a, guide_total),
                        "control_fraction": safe_div(c, ctrl_total),
                        "delta_fraction": safe_div(a, guide_total) - safe_div(c, ctrl_total),
                        "log2OR_guide_vs_control": log2_or(a, b, c, d),
                        "fisher_p": stats.fisher_exact([[a, b], [c, d]])[1],
                        "is_highlight_target": target in highlight_targets,
                    })

    celltype_enrich = pd.DataFrame(celltype_rows)
    if len(celltype_enrich):
        celltype_enrich["fisher_fdr"] = bh_fdr(celltype_enrich["fisher_p"].values)
        celltype_enrich.to_csv(os.path.join(outdir, "target_celltype_one_vs_rest_enrichment.csv"), index=False)

    guide_celltype_enrich = pd.DataFrame(guide_celltype_rows)
    if len(guide_celltype_enrich):
        guide_celltype_enrich["fisher_fdr"] = bh_fdr(guide_celltype_enrich["fisher_p"].values)
        guide_celltype_enrich.to_csv(os.path.join(outdir, "guide_celltype_one_vs_rest_enrichment.csv"), index=False)

    comp_dist = pd.DataFrame(composition_distance_rows)
    if len(comp_dist):
        comp_dist.to_csv(os.path.join(outdir, "target_composition_distance_from_control.csv"), index=False)

    # ----------------------
    # Binary state contrasts
    # ----------------------
    binary_contrasts = [
        {
            "name": "D1_P13high_vs_P30high",
            "day": "d1",
            "event_label": "P13high",
            "ref_label": "P30high",
            "state_series": cell_df["D1_program_state"],
            "event_values": ["P13_high"],
            "ref_values": ["P30_high"],
        },
        {
            "name": "D2_LatMeso_vs_PreEC",
            "day": "d2",
            "event_label": "latMeso",
            "ref_label": "preEC",
            "state_series": cell_df["cell_type_clean"],
            "event_values": DEFAULT_D2_LAT_MESO,
            "ref_values": DEFAULT_D2_PRE_EC,
        },
        {
            "name": "D2_P17high_vs_P15high",
            "day": "d2",
            "event_label": "P17high",
            "ref_label": "P15high",
            "state_series": cell_df["D2_program_state"],
            "event_values": ["P17_high"],
            "ref_values": ["P15_high"],
        },
        {
            "name": "D3_AltBroad_vs_EC",
            "day": "d3",
            "event_label": "alt",
            "ref_label": "endo",
            "state_series": cell_df["cell_type_clean"],
            "event_values": DEFAULT_D3_ALT,
            "ref_values": DEFAULT_D3_EC,
        },
        {
            "name": "D3_Mesenchymal_vs_EC",
            "day": "d3",
            "event_label": "mes",
            "ref_label": "endo",
            "state_series": cell_df["cell_type_clean"],
            "event_values": DEFAULT_D3_MES,
            "ref_values": DEFAULT_D3_EC,
        },
    ]

    binary_target_rows = []
    binary_guide_rows = []
    axis_target_rows = []
    axis_guide_rows = []
    leadlag_rows = []

    for mode, masks in modes.items():
        ctrl_base = nt_mask & masks["control_selector"]

        for target in selected_targets:
            carrier_base = target_assignment[target] & masks["target_selector"]

            target_lead = {"mode": mode, "target": target, "is_highlight_target": target in highlight_targets}

            for contrast in binary_contrasts:
                day_mask = cell_df["day_clean"].eq(contrast["day"])
                valid = day_mask
                state = contrast["state_series"]
                event = valid & state.isin(contrast["event_values"])
                ref = valid & state.isin(contrast["ref_values"])

                tmask = carrier_base
                cmask = ctrl_base
                a = int((tmask & event).sum())
                b = int((tmask & ref).sum())
                c = int((cmask & event).sum())
                d = int((cmask & ref).sum())
                ntot = a + b
                ctot = c + d

                if ntot >= min_target_cells and ctot >= min_ref_cells:
                    pval = stats.fisher_exact([[a, b], [c, d]])[1] if (a + b + c + d) else np.nan
                    lor = log2_or(a, b, c, d)
                    row = {
                        "mode": mode,
                        "target": target,
                        "contrast": contrast["name"],
                        "day": contrast["day"],
                        "event_label": contrast["event_label"],
                        "ref_label": contrast["ref_label"],
                        "target_event_n": a,
                        "target_ref_n": b,
                        "control_event_n": c,
                        "control_ref_n": d,
                        "target_total_event_ref": ntot,
                        "control_total_event_ref": ctot,
                        "target_event_fraction": safe_div(a, ntot),
                        "control_event_fraction": safe_div(c, ctot),
                        "delta_event_fraction": safe_div(a, ntot) - safe_div(c, ctot),
                        "log2OR_event_vs_ref": lor,
                        "fisher_p": pval,
                        "is_highlight_target": target in highlight_targets,
                    }
                    binary_target_rows.append(row)
                    prefix = contrast["name"]
                    target_lead[f"{prefix}_log2OR"] = lor
                    target_lead[f"{prefix}_delta_fraction"] = row["delta_event_fraction"]
                    target_lead[f"{prefix}_target_fraction"] = row["target_event_fraction"]

            for axis, spec in AXIS_DEFS.items():
                day_mask = cell_df["day_clean"].eq(spec["day"])
                x = cell_df.loc[carrier_base & day_mask, axis]
                y = cell_df.loc[ctrl_base & day_mask, axis]
                if x.dropna().shape[0] >= min_target_cells and y.dropna().shape[0] >= min_ref_cells:
                    dist = summarize_distribution(x, y)
                    row = {
                        "mode": mode,
                        "target": target,
                        "day": spec["day"],
                        "axis": axis,
                        "axis_description": spec["description"],
                        "is_highlight_target": target in highlight_targets,
                    }
                    row.update(dist)
                    axis_target_rows.append(row)
                    target_lead[f"{axis}_delta_mean"] = dist["delta_mean"]
                    target_lead[f"{axis}_target_mean"] = dist["target_mean"]
                    target_lead[f"{axis}_control_mean"] = dist["control_mean"]

            # merge composition distances into lead table
            if len(comp_dist):
                for day in all_days:
                    subd = comp_dist[(comp_dist["mode"].eq(mode)) & (comp_dist["target"].eq(target)) & (comp_dist["day"].eq(day))]
                    if len(subd):
                        target_lead[f"{day}_JS_composition_distance"] = subd["jensen_shannon_distance"].iloc[0]
                        target_lead[f"{day}_TV_composition_distance"] = subd["total_variation_distance"].iloc[0]

            leadlag_rows.append(target_lead)

        # guide-level axes and binary contrasts
        for j, guide in enumerate(guide_assignment.columns):
            target_series = guide_target.loc[guide_target["guide"].astype(str).eq(str(guide)), "target"]
            if len(target_series) == 0:
                continue
            target = str(target_series.iloc[0])
            if target not in selected_targets:
                continue
            carrier_base = guide_assignment.iloc[:, j] & masks["target_selector"]

            for axis, spec in AXIS_DEFS.items():
                day_mask = cell_df["day_clean"].eq(spec["day"])
                x = cell_df.loc[carrier_base & day_mask, axis]
                y = cell_df.loc[ctrl_base & day_mask, axis]
                if x.dropna().shape[0] >= min_guide_cells and y.dropna().shape[0] >= min_ref_cells:
                    row = {
                        "mode": mode,
                        "guide": guide,
                        "target": target,
                        "day": spec["day"],
                        "axis": axis,
                        "n_guide_cells": int(x.dropna().shape[0]),
                        "n_control_cells": int(y.dropna().shape[0]),
                        "guide_mean": x.mean(),
                        "control_mean": y.mean(),
                        "effect_vs_control": x.mean() - y.mean(),
                        "guide_median": x.median(),
                        "control_median": y.median(),
                        "is_highlight_target": target in highlight_targets,
                    }
                    axis_guide_rows.append(row)

            for contrast in binary_contrasts:
                day_mask = cell_df["day_clean"].eq(contrast["day"])
                state = contrast["state_series"]
                event = day_mask & state.isin(contrast["event_values"])
                ref = day_mask & state.isin(contrast["ref_values"])
                a = int((carrier_base & event).sum())
                b = int((carrier_base & ref).sum())
                c = int((ctrl_base & event).sum())
                d = int((ctrl_base & ref).sum())
                ntot = a + b
                ctot = c + d
                if ntot >= min_guide_cells and ctot >= min_ref_cells:
                    binary_guide_rows.append({
                        "mode": mode,
                        "guide": guide,
                        "target": target,
                        "contrast": contrast["name"],
                        "day": contrast["day"],
                        "event_label": contrast["event_label"],
                        "ref_label": contrast["ref_label"],
                        "guide_event_n": a,
                        "guide_ref_n": b,
                        "control_event_n": c,
                        "control_ref_n": d,
                        "guide_total_event_ref": ntot,
                        "control_total_event_ref": ctot,
                        "guide_event_fraction": safe_div(a, ntot),
                        "control_event_fraction": safe_div(c, ctot),
                        "delta_event_fraction": safe_div(a, ntot) - safe_div(c, ctot),
                        "log2OR_event_vs_ref": log2_or(a, b, c, d),
                        "is_highlight_target": target in highlight_targets,
                    })

    binary_target = pd.DataFrame(binary_target_rows)
    if len(binary_target):
        binary_target["fisher_fdr"] = bh_fdr(binary_target["fisher_p"].values)
        binary_target.to_csv(os.path.join(outdir, "target_binary_state_enrichment.csv"), index=False)

    binary_guide = pd.DataFrame(binary_guide_rows)
    if len(binary_guide):
        binary_guide.to_csv(os.path.join(outdir, "guide_binary_state_enrichment.csv"), index=False)

    axis_target = pd.DataFrame(axis_target_rows)
    if len(axis_target):
        for pcol in ["ks_p", "mannwhitney_p"]:
            axis_target[pcol.replace("_p", "_fdr")] = bh_fdr(axis_target[pcol].values)
        axis_target.to_csv(os.path.join(outdir, "target_axis_distribution_tests.csv"), index=False)

    axis_guide = pd.DataFrame(axis_guide_rows)
    if len(axis_guide):
        axis_guide.to_csv(os.path.join(outdir, "guide_axis_effects.csv"), index=False)

    leadlag = pd.DataFrame(leadlag_rows)
    if len(leadlag):
        leadlag.to_csv(os.path.join(outdir, "target_leadlag_metrics.csv"), index=False)

    # Guide concordance summary.
    if len(axis_guide):
        conv_rows = []
        for mode in axis_guide["mode"].unique():
            for target in sorted(axis_guide["target"].unique()):
                for axis in sorted(axis_guide["axis"].unique()):
                    vals = axis_guide[
                        axis_guide["mode"].eq(mode)
                        & axis_guide["target"].eq(target)
                        & axis_guide["axis"].eq(axis)
                    ]["effect_vs_control"].dropna()
                    if len(vals) == 0:
                        continue
                    conv_rows.append({
                        "mode": mode,
                        "target": target,
                        "axis": axis,
                        "n_guides": len(vals),
                        "mean_guide_effect": vals.mean(),
                        "median_guide_effect": vals.median(),
                        "fraction_positive": (vals > 0).mean(),
                        "fraction_negative": (vals < 0).mean(),
                        "one_sample_t_p": stats.ttest_1samp(vals, 0).pvalue if len(vals) >= 2 else np.nan,
                        "signed_rank_p": stats.wilcoxon(vals).pvalue if len(vals) >= 2 and not np.all(vals == 0) else np.nan,
                        "is_highlight_target": target in highlight_targets,
                    })
        guide_conv = pd.DataFrame(conv_rows)
        if len(guide_conv):
            guide_conv.to_csv(os.path.join(outdir, "guide_concordance_summary.csv"), index=False)

    # ----------------------
    # Plotting cell score sample
    # ----------------------
    plot_targets = list(dict.fromkeys(highlight_targets + selected_targets[:0]))
    rows = []
    for target in plot_targets:
        if target not in target_assignment.columns:
            continue
        for day in all_days:
            day_mask = cell_df["day_clean"].eq(day)
            target_idx = cell_df.index[target_assignment[target] & day_mask].to_numpy()
            ctrl_idx = cell_df.index[nt_mask & day_mask].to_numpy()
            if len(target_idx) > max_cells_per_group_for_plotting:
                target_idx = rng.choice(target_idx, size=max_cells_per_group_for_plotting, replace=False)
            if len(ctrl_idx) > max_cells_per_group_for_plotting:
                ctrl_idx = rng.choice(ctrl_idx, size=max_cells_per_group_for_plotting, replace=False)
            for idx, condition in [(target_idx, "target_carrier"), (ctrl_idx, "non_targeting_control")]:
                if len(idx) == 0:
                    continue
                tmp = cell_df.loc[idx, ["cell_id", "day_clean", "cell_type_clean", "n_total_guides", "n_total_targets", "n_selected_targets"] + needed_cols + list(AXIS_DEFS.keys())].copy()
                tmp["target"] = target
                tmp["condition"] = condition
                rows.append(tmp)
    if rows:
        pd.concat(rows, ignore_index=True).to_csv(os.path.join(outdir, "highlight_target_cell_score_sample_for_plotting.csv"), index=False)

    if export_cell_table:
        keep_cols = ["cell_id", "day_clean", "cell_type_clean", "n_total_guides", "n_total_targets", "n_selected_targets", "is_non_targeting_carrier"] + needed_cols + list(AXIS_DEFS.keys())
        cell_df[keep_cols].to_csv(os.path.join(outdir, "all_cells_program_axes_metadata.csv"), index=False)

    # ----------------------
    # Optional neighborhood enrichment / Milo-like local DA table
    # ----------------------
    if export_neighborhoods:
        embed_key = None
        for key in ["X_pca", "X_umap", "X_tsne"]:
            if key in adata.obsm:
                embed_key = key
                break
        if embed_key is not None:
            Xemb = np.asarray(adata.obsm[embed_key])
        else:
            Xemb = program_raw[needed_cols].values
            embed_key = "cnmf_usage_selected_programs"

        n_cells = Xemb.shape[0]
        if n_cells > 0:
            seed_n = min(neighborhood_n_seeds, n_cells)
            seed_idx = rng.choice(np.arange(n_cells), size=seed_n, replace=False)
            nn = NearestNeighbors(n_neighbors=min(neighborhood_k, n_cells), metric="euclidean")
            nn.fit(Xemb)
            neigh_ind = nn.kneighbors(Xemb[seed_idx], return_distance=False)
            nbh_rows = []

            for nid, inds in enumerate(neigh_ind):
                nbh_mask = pd.Series(False, index=cell_df.index)
                nbh_mask.iloc[inds] = True
                center = Xemb[seed_idx[nid]]
                nbh_day = cell_df.iloc[inds]["day_clean"].mode()
                nbh_day = str(nbh_day.iloc[0]) if len(nbh_day) else "mixed"
                for target in highlight_targets:
                    if target not in target_assignment.columns:
                        continue
                    for mode, masks in modes.items():
                        target_sel = masks["target_selector"]
                        control_sel = masks["control_selector"]
                        carrier = target_assignment[target] & target_sel
                        ctrl = nt_mask & control_sel
                        # enrichment of carrier/control in neighborhood
                        a = int((carrier & nbh_mask).sum())
                        b = int((carrier & ~nbh_mask).sum())
                        c = int((ctrl & nbh_mask).sum())
                        d = int((ctrl & ~nbh_mask).sum())
                        if (a + b) < min_target_cells or (c + d) < min_ref_cells:
                            continue
                        row = {
                            "neighborhood_id": nid,
                            "embedding_key": embed_key,
                            "mode": mode,
                            "target": target,
                            "dominant_day": nbh_day,
                            "n_cells_neighborhood": int(nbh_mask.sum()),
                            "target_in_neighborhood": a,
                            "target_outside_neighborhood": b,
                            "control_in_neighborhood": c,
                            "control_outside_neighborhood": d,
                            "log2OR_neighborhood_enrichment": log2_or(a, b, c, d),
                            "fisher_p": stats.fisher_exact([[a, b], [c, d]])[1],
                            "mean_embed1": float(Xemb[inds, 0].mean()) if Xemb.shape[1] >= 1 else np.nan,
                            "mean_embed2": float(Xemb[inds, 1].mean()) if Xemb.shape[1] >= 2 else np.nan,
                        }
                        for axis in AXIS_DEFS:
                            row[f"{axis}_mean"] = cell_df.iloc[inds][axis].mean()
                        nbh_rows.append(row)
            nbh = pd.DataFrame(nbh_rows)
            if len(nbh):
                nbh["fisher_fdr"] = bh_fdr(nbh["fisher_p"].values)
                nbh.to_csv(os.path.join(outdir, "neighborhood_enrichment_selected_targets.csv"), index=False)

    # Manifest / README
    manifest = pd.DataFrame({
        "file": sorted(os.listdir(outdir)),
        "path": [os.path.abspath(os.path.join(outdir, fn)) for fn in sorted(os.listdir(outdir))]
    })
    manifest.to_csv(os.path.join(outdir, "export_manifest.csv"), index=False)

    with open(os.path.join(outdir, "README.txt"), "w") as f:
        f.write(
            "Lineage split cell-level export outputs.\n\n"
            "Interpretation caveat:\n"
            "These are high-MOI target-carrier enrichment analyses, not single-perturbation lineage tracing.\n\n"
            "Core outputs:\n"
            "- target_celltype_one_vs_rest_enrichment.csv\n"
            "- guide_celltype_one_vs_rest_enrichment.csv\n"
            "- target_binary_state_enrichment.csv\n"
            "- guide_binary_state_enrichment.csv\n"
            "- target_axis_distribution_tests.csv\n"
            "- guide_axis_effects.csv\n"
            "- guide_concordance_summary.csv\n"
            "- target_composition_distance_from_control.csv\n"
            "- target_leadlag_metrics.csv\n"
            "- pseudobulk_counts_by_sample_target_day_celltype.csv, if sample metadata exists\n"
            "- neighborhood_enrichment_selected_targets.csv, if export_neighborhoods=True\n"
            "- highlight_target_cell_score_sample_for_plotting.csv\n\n"
            "Axis orientation:\n"
            "- D1_P13_minus_P30: positive = P13 over P30\n"
            "- D2_P17_minus_P15: positive = P17 over P15\n"
            "- D2_YAPLat_minus_PreEC: positive = mean(P17,P34) over mean(P15,P23)\n"
            "- D3_ALT_minus_EC: positive = mean(P27,P34,P35) over mean(P2,P3,P8,P11)\n"
        )

    print(f"Wrote lineage split cell-level exports to: {outdir}")
    print("Selected targets:", len(selected_targets))
    print("Highlight targets present:", highlight_targets)
    print("Files:")
    for fn in sorted(os.listdir(outdir)):
        print(" -", fn)

    return {
        "outdir": outdir,
        "selected_targets": selected_targets,
        "highlight_targets": highlight_targets,
        "manifest": manifest,
    }


# ============================================================
# CLI
# ============================================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata", required=True, help="Path to .h5ad file.")
    parser.add_argument("--guide-target-csv", default=None)
    parser.add_argument("--outdir", default="lineage_split_celllevel_exports")
    parser.add_argument("--array-program-index-base", type=int, default=0)
    parser.add_argument("--selected-targets", default=None, help="Comma-separated target list. Default: all targets.")
    parser.add_argument("--highlight-targets", default="NF2,LATS2,SAV1,TAOK1,AMOTL2")
    parser.add_argument("--export-cell-table", action="store_true")
    parser.add_argument("--export-neighborhoods", action="store_true")
    parser.add_argument("--neighborhood-k", type=int, default=80)
    parser.add_argument("--neighborhood-n-seeds", type=int, default=800)
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.adata)

    selected = None if args.selected_targets is None else [x.strip() for x in args.selected_targets.split(",") if x.strip()]
    highlight = [x.strip() for x in args.highlight_targets.split(",") if x.strip()]

    export_lineage_split_celllevel_csvs(
        adata,
        guide_target_csv=args.guide_target_csv,
        outdir=args.outdir,
        selected_targets=selected,
        highlight_targets=highlight,
        array_program_index_base=args.array_program_index_base,
        export_cell_table=args.export_cell_table,
        export_neighborhoods=args.export_neighborhoods,
        neighborhood_k=args.neighborhood_k,
        neighborhood_n_seeds=args.neighborhood_n_seeds,
    )


if __name__ == "__main__":
    main()
