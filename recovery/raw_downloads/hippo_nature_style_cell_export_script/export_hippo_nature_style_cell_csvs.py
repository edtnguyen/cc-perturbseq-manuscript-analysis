#!/usr/bin/env python
# export_hippo_nature_style_cell_csvs.py
#
# A more Nature-paper-style cell-level export pipeline for Hippo analyses.
#
# Main goals:
#   1. Treat guides as pseudo-replicates, not just pooled cells.
#   2. Export cell-state composition tables analogous to organoid-by-organoid bar plots.
#   3. Export distribution/trajectory-shift tests analogous to pseudotime KS plots.
#   4. Export module-score tables analogous to WGCNA/module-score violin plots.
#   5. Export convergence tables across Hippo targets/guides.
#
# Notebook usage:
#   from export_hippo_nature_style_cell_csvs import export_hippo_nature_style_cell_csvs
#   export_hippo_nature_style_cell_csvs(
#       adata,
#       guide_target_csv="combined_TargetKD.csv",
#       outdir="hippo_nature_style_exports"
#   )
#
# Command line:
#   python export_hippo_nature_style_cell_csvs.py --adata adata.h5ad --guide_target_csv combined_TargetKD.csv --outdir hippo_nature_style_exports

import os
import argparse
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact, ks_2samp, mannwhitneyu, wasserstein_distance


HIPPO_CORE_DEFAULT = [
    "NF2", "LATS1", "LATS2", "SAV1", "STK3", "STK4",
    "MOB1A", "MOB1B", "TAOK1", "TAOK2", "TAOK3",
    "AMOT", "AMOTL1", "AMOTL2",
]

PROGRAMS = [2, 3, 8, 11, 13, 15, 17, 23, 27, 30, 34, 35]

CELLTYPE_GROUPS = {
    "D1_iPSC_like": ["iPSC", "Cyc-iPSC"],
    "D2_PreEC": ["Pre-EC"],
    "D2_LatMeso": ["Lat-Meso"],
    "D2_MesenAlt": ["Neu-Mesen", "Fibro-Mesen"],
    "D3_EC": ["artEC", "ECM-EC"],
    "D3_ALT_broad": ["Fibro-Mesen", "Neu-Mesen"],
    "D3_MES_strict": ["Fibro-Mesen"],
}

GENE_MODULES = {
    "YAP_contractile": ["ANKRD1", "CCN1", "CYR61", "CCN2", "CTGF", "AXL", "TAGLN", "ACTA2", "TGFB2", "COL1A1", "COL1A2"],
    "EC_core": ["KDR", "SOX7", "TAL1", "FLI1", "ERG", "FLT1", "NRP1", "NRP2", "PECAM1", "CDH5", "ESM1"],
    "P17_epithelial_alt": ["CDH1", "KRT8", "KRT19", "DSP", "EPCAM", "TFAP2A", "GRHL2", "PATJ"],
    "P30_mesoderm_entry": ["TBXT", "MSGN1", "DKK1", "DKK4", "NKD1", "FST"],
    "P13_nonvascular_WNT_neural": ["PAX3", "SEMA3C", "EBF3", "GAD1", "WLS", "AXIN2", "RNF43", "NKD1"],
    "PreEC_program_genes": ["KDR", "NRP2", "TAL1", "FLI1", "SOX7", "APLNR", "MEF2C"],
}


def _clean_day(x):
    return str(x).strip().lower()


def _as_bool_array(x):
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x) > 0


def _infer_col(df, candidates):
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    for c in df.columns:
        cl = c.lower()
        for cand in candidates:
            if cand.lower() in cl:
                return c
    return None


def _looks_nt(x):
    s = str(x).lower()
    patterns = [
        "non-target", "non_target", "nontarget", "non targeting",
        "negative", "control", "scramble", "scrambled", "ntc", "safe", "ctrl"
    ]
    return s in {"nt", "neg", "ctrl"} or any(p in s for p in patterns)


def _parse_target_from_guide(g):
    g = str(g)
    for sep in ["|", ":", "_", "-"]:
        if sep in g:
            return g.split(sep)[0]
    return g


def build_guide_target(guide_names, guide_target_csv=None):
    guide_names = pd.Index([str(x) for x in guide_names], name="guide")
    base = pd.DataFrame({"guide": guide_names})
    if guide_target_csv and os.path.exists(guide_target_csv):
        gt = pd.read_csv(guide_target_csv)
        guide_col = _infer_col(gt, ["guide", "guide_name", "sgRNA", "sgrna", "gRNA", "spacer", "oligo"])
        target_col = _infer_col(gt, ["target", "target_name", "target_gene", "gene", "symbol", "TargetGene"])
        if guide_col is None or target_col is None:
            raise ValueError(f"Could not infer guide/target columns from {guide_target_csv}. Columns: {list(gt.columns)}")
        gt = gt[[guide_col, target_col]].copy()
        gt.columns = ["guide", "target"]
        gt["guide"] = gt["guide"].astype(str)
        gt["target"] = gt["target"].astype(str)
        out = base.merge(gt, on="guide", how="left")
        out["target"] = out["target"].fillna(out["guide"].map(_parse_target_from_guide))
    else:
        out = base.copy()
        out["target"] = out["guide"].map(_parse_target_from_guide)
    out["is_nt"] = out["guide"].map(_looks_nt) | out["target"].map(_looks_nt)
    return out


def get_program_df(adata, cnmf_key="cnmf_usage", programs=PROGRAMS):
    usage = adata.obsm[cnmf_key]
    if isinstance(usage, pd.DataFrame):
        u = usage.copy()
        u.index = adata.obs_names
        out = {}
        for p in programs:
            choices = [p, str(p), f"P{p}", f"Program {p}", f"program {p}", f"Program_{p}"]
            found = None
            for c in u.columns:
                if c in choices or str(c) in choices:
                    found = c
                    break
            if found is None:
                found = u.columns[p - 1]
            out[f"Program {p}"] = u[found].astype(float).values
        return pd.DataFrame(out, index=adata.obs_names)
    arr = np.asarray(usage)
    out = {f"Program {p}": arr[:, p - 1].astype(float) for p in programs}
    return pd.DataFrame(out, index=adata.obs_names)


def zscore_by_day(program_df, obs, day_col, background_mask=None, min_background=50):
    z = pd.DataFrame(index=program_df.index, columns=program_df.columns, dtype=float)
    days = obs[day_col].map(_clean_day)
    if background_mask is None:
        background_mask = pd.Series(False, index=program_df.index)
    for d in sorted(days.unique()):
        day_mask = days.eq(d)
        bg = day_mask & background_mask
        if bg.sum() < min_background:
            bg = day_mask
        mu = program_df.loc[bg].mean(axis=0)
        sd = program_df.loc[bg].std(axis=0, ddof=1).replace(0, np.nan)
        z.loc[day_mask] = (program_df.loc[day_mask] - mu) / sd
    return z.astype(float)


def add_fate_axes(cell_df):
    cell_df["D1_alt_axis_P13_minus_P30"] = cell_df["Program 13"] - cell_df["Program 30"]
    cell_df["D2_alt_axis_P17_minus_P15"] = cell_df["Program 17"] - cell_df["Program 15"]
    cell_df["D2_YAPLat_axis_P17P34_minus_P15P23"] = (
        cell_df[["Program 17", "Program 34"]].mean(axis=1)
        - cell_df[["Program 15", "Program 23"]].mean(axis=1)
    )
    cell_df["D3_ALT_axis_P27P34P35_minus_ECcore"] = (
        cell_df[["Program 27", "Program 34", "Program 35"]].mean(axis=1)
        - cell_df[["Program 2", "Program 3", "Program 8", "Program 11"]].mean(axis=1)
    )
    cell_df["D1_program_state"] = np.where(cell_df["Program 13"] > cell_df["Program 30"], "P13_high", "P30_high")
    cell_df["D2_program_state"] = np.where(cell_df["Program 17"] > cell_df["Program 15"], "P17_high", "P15_high")
    return cell_df


def cliff_delta(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    if len(x) == 0 or len(y) == 0:
        return np.nan
    # rank-based efficient approximation:
    # delta = (2*U)/(n_x*n_y) - 1, where U counts x > y + half ties.
    try:
        U = mannwhitneyu(x, y, alternative="two-sided").statistic
        return (2 * U) / (len(x) * len(y)) - 1
    except Exception:
        return np.nan


def compare_distributions(x, y):
    x = np.asarray(pd.Series(x).dropna(), dtype=float)
    y = np.asarray(pd.Series(y).dropna(), dtype=float)
    if len(x) < 2 or len(y) < 2:
        return {"ks_stat": np.nan, "ks_p": np.nan, "mw_p": np.nan, "wasserstein": np.nan, "cliff_delta": np.nan}
    ks = ks_2samp(x, y)
    mw = mannwhitneyu(x, y, alternative="two-sided")
    return {
        "ks_stat": ks.statistic,
        "ks_p": ks.pvalue,
        "mw_p": mw.pvalue,
        "wasserstein": wasserstein_distance(x, y),
        "cliff_delta": cliff_delta(x, y),
    }


def get_expr_matrix(adata, use_raw=False):
    if use_raw and getattr(adata, "raw", None) is not None:
        X = adata.raw.X
        var_names = pd.Index(adata.raw.var_names)
    else:
        X = adata.X
        var_names = pd.Index(adata.var_names)
    if sparse.issparse(X):
        return X.tocsr(), var_names
    return np.asarray(X), var_names


def compute_gene_modules(adata, cell_index, modules=GENE_MODULES, use_raw=False):
    X, var_names = get_expr_matrix(adata, use_raw=use_raw)
    var_upper = pd.Series(np.arange(len(var_names)), index=var_names.astype(str).str.upper())
    out = pd.DataFrame(index=cell_index)
    gene_manifest = []
    for name, genes in modules.items():
        idx = []
        present = []
        for g in genes:
            gu = g.upper()
            if gu in var_upper.index:
                idx.append(int(var_upper.loc[gu]))
                present.append(g)
        gene_manifest.append({"module": name, "genes_requested": ",".join(genes), "genes_present": ",".join(present), "n_present": len(present)})
        if len(idx) == 0:
            out[name] = np.nan
            continue
        vals = X[:, idx]
        if sparse.issparse(vals):
            vals = vals.toarray()
        vals = np.asarray(vals, dtype=float)
        # Gene-wise z, then average. Works for log-normalized or raw-ish matrices.
        mu = np.nanmean(vals, axis=0)
        sd = np.nanstd(vals, axis=0, ddof=1)
        sd[sd == 0] = np.nan
        z = (vals - mu) / sd
        out[name] = np.nanmean(z, axis=1)
    return out, pd.DataFrame(gene_manifest)


def log2_or(a, b, c, d, pc=0.5):
    return np.log2(((a + pc) / (b + pc)) / ((c + pc) / (d + pc)))


def export_hippo_nature_style_cell_csvs(
    adata,
    guide_target_csv="combined_TargetKD.csv",
    outdir="hippo_nature_style_exports",
    hippo_core=None,
    cnmf_key="cnmf_usage",
    guide_assignment_key="guide_assignment",
    day_col="day",
    celltype_col="cell_type",
    use_raw_for_modules=False,
    min_background=50,
    max_background_cells_per_target_day=5000,
    random_seed=7,
):
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(random_seed)

    if hippo_core is None:
        hippo_core = HIPPO_CORE_DEFAULT

    obs = adata.obs.copy()
    obs.index = adata.obs_names
    obs["day_clean"] = obs[day_col].map(_clean_day)
    obs["cell_type_clean"] = obs[celltype_col].astype(str)

    guide_names = pd.Index([str(x) for x in adata.uns["guide_names"]])
    guide_assignment = _as_bool_array(adata.obsm[guide_assignment_key])

    guide_target = build_guide_target(guide_names, guide_target_csv)
    guide_target.to_csv(os.path.join(outdir, "guide_target_mapping_used.csv"), index=False)

    nt_guides = guide_target.loc[guide_target["is_nt"], "guide"].tolist()
    nt_idx = np.where(guide_names.isin(nt_guides))[0]
    nt_carrier = pd.Series(
        guide_assignment[:, nt_idx].any(axis=1) if len(nt_idx) else np.zeros(adata.n_obs, dtype=bool),
        index=adata.obs_names,
        name="nt_carrier",
    )

    program_raw = get_program_df(adata, cnmf_key=cnmf_key, programs=PROGRAMS)
    program_z = zscore_by_day(program_raw, obs, day_col="day_clean", background_mask=nt_carrier, min_background=min_background)

    cell_df = pd.concat([obs[["day_clean", "cell_type_clean"]], program_z], axis=1)
    cell_df = add_fate_axes(cell_df)

    module_scores, module_manifest = compute_gene_modules(adata, cell_df.index, use_raw=use_raw_for_modules)
    module_manifest.to_csv(os.path.join(outdir, "gene_module_manifest.csv"), index=False)
    cell_df = pd.concat([cell_df, module_scores], axis=1)

    embed_key = None
    for key in ["X_umap", "X_tsne", "umap", "tsne"]:
        if key in adata.obsm:
            embed_key = key
            emb = np.asarray(adata.obsm[key])
            cell_df["embed1"] = emb[:, 0]
            cell_df["embed2"] = emb[:, 1]
            cell_df["embedding_key"] = key
            break

    cell_df.insert(0, "cell_id", cell_df.index)
    cell_df.to_csv(os.path.join(outdir, "all_cells_program_axes_and_modules.csv"), index=False)

    non_nt_targets = sorted(guide_target.loc[~guide_target["is_nt"], "target"].dropna().astype(str).unique())
    target_to_guides = {
        t: guide_target.loc[guide_target["target"].astype(str).eq(t) & ~guide_target["is_nt"], "guide"].tolist()
        for t in non_nt_targets
    }
    target_carriers = {}
    for target, guides in target_to_guides.items():
        idx = np.where(guide_names.isin(guides))[0]
        target_carriers[target] = pd.Series(guide_assignment[:, idx].any(axis=1), index=adata.obs_names) if len(idx) else pd.Series(False, index=adata.obs_names)

    target_carrier_df = pd.DataFrame(target_carriers)
    target_count = target_carrier_df.sum(axis=1)

    present_hippo = [g for g in hippo_core if g in target_carriers]
    pd.DataFrame({"target": present_hippo}).to_csv(os.path.join(outdir, "hippo_targets_present.csv"), index=False)

    # Target and guide carrier counts by day/cell type
    count_rows = []
    for target in non_nt_targets:
        carrier = target_carriers[target]
        for day in sorted(obs["day_clean"].unique()):
            for ct in sorted(obs["cell_type_clean"].unique()):
                mask = obs["day_clean"].eq(day) & obs["cell_type_clean"].eq(ct)
                count_rows.append({
                    "level": "target",
                    "target": target,
                    "guide": "",
                    "day": day,
                    "cell_type": ct,
                    "n_carrier": int((carrier & mask).sum()),
                    "n_total_state_day": int(mask.sum()),
                    "n_target_day": int((carrier & obs["day_clean"].eq(day)).sum()),
                    "fraction_within_target_day": (carrier & mask).sum() / max((carrier & obs["day_clean"].eq(day)).sum(), 1),
                    "is_hippo_core": target in present_hippo,
                })
    for gi, guide in enumerate(guide_names):
        row = guide_target.loc[guide_target["guide"].eq(guide)].iloc[0]
        if bool(row["is_nt"]):
            continue
        target = str(row["target"])
        carrier = pd.Series(guide_assignment[:, gi], index=adata.obs_names)
        for day in sorted(obs["day_clean"].unique()):
            for ct in sorted(obs["cell_type_clean"].unique()):
                mask = obs["day_clean"].eq(day) & obs["cell_type_clean"].eq(ct)
                count_rows.append({
                    "level": "guide",
                    "target": target,
                    "guide": guide,
                    "day": day,
                    "cell_type": ct,
                    "n_carrier": int((carrier & mask).sum()),
                    "n_total_state_day": int(mask.sum()),
                    "n_target_day": int((carrier & obs["day_clean"].eq(day)).sum()),
                    "fraction_within_target_day": (carrier & mask).sum() / max((carrier & obs["day_clean"].eq(day)).sum(), 1),
                    "is_hippo_core": target in present_hippo,
                })
    pd.DataFrame(count_rows).to_csv(os.path.join(outdir, "composition_counts_by_target_and_guide.csv"), index=False)

    # Contrasts analogous to cell-type proportion tests.
    contrasts = [
        ("D1_P13high_vs_P30high", "d1", cell_df["D1_program_state"].eq("P13_high"), cell_df["D1_program_state"].eq("P30_high")),
        ("D2_LatMeso_vs_PreEC", "d2", obs["cell_type_clean"].eq("Lat-Meso"), obs["cell_type_clean"].eq("Pre-EC")),
        ("D2_P17high_vs_P15high", "d2", cell_df["D2_program_state"].eq("P17_high"), cell_df["D2_program_state"].eq("P15_high")),
        ("D3_ALTbroad_vs_EC", "d3", obs["cell_type_clean"].isin(["Fibro-Mesen", "Neu-Mesen"]), obs["cell_type_clean"].isin(["artEC", "ECM-EC"])),
        ("D3_FibroMesen_vs_EC", "d3", obs["cell_type_clean"].eq("Fibro-Mesen"), obs["cell_type_clean"].isin(["artEC", "ECM-EC"])),
    ]

    enrich_rows = []
    modes = {
        "all_carriers": {t: target_carriers[t] for t in non_nt_targets},
        "single_selected_target_carriers": {t: target_carriers[t] & target_count.eq(1) for t in non_nt_targets},
    }

    for mode, carrier_dict in modes.items():
        for target, carrier in carrier_dict.items():
            bg_base = (nt_carrier & ~carrier) if nt_carrier.sum() >= min_background else (~carrier)
            for contrast_name, day, event_mask, ref_mask in contrasts:
                day_mask = obs["day_clean"].eq(day)
                event = event_mask & day_mask
                ref = ref_mask & day_mask
                bg = bg_base & day_mask
                if bg.sum() < min_background:
                    bg = (~carrier) & day_mask

                a = int((carrier & event).sum())
                b = int((carrier & ref).sum())
                c = int((bg & event).sum())
                d = int((bg & ref).sum())

                try:
                    _, fp = fisher_exact([[a, b], [c, d]])
                except Exception:
                    fp = np.nan

                enrich_rows.append({
                    "mode": mode,
                    "target": target,
                    "contrast": contrast_name,
                    "day": day,
                    "is_hippo_core": target in present_hippo,
                    "target_event_n": a,
                    "target_ref_n": b,
                    "background_event_n": c,
                    "background_ref_n": d,
                    "target_event_fraction": a / (a + b) if (a + b) else np.nan,
                    "background_event_fraction": c / (c + d) if (c + d) else np.nan,
                    "log2OR_target_vs_background": log2_or(a, b, c, d) if (a + b) and (c + d) else np.nan,
                    "fisher_exact_p": fp,
                })
    pd.DataFrame(enrich_rows).to_csv(os.path.join(outdir, "state_enrichment_target_level.csv"), index=False)

    # Guide-level analogue of organoid replicates.
    guide_effect_rows = []
    guide_level_axes = [
        ("d1", "D1_alt_axis_P13_minus_P30"),
        ("d2", "D2_alt_axis_P17_minus_P15"),
        ("d2", "D2_YAPLat_axis_P17P34_minus_P15P23"),
        ("d3", "D3_ALT_axis_P27P34P35_minus_ECcore"),
    ]
    module_names = list(GENE_MODULES.keys())

    for gi, guide in enumerate(guide_names):
        gmap = guide_target.loc[guide_target["guide"].eq(guide)].iloc[0]
        if bool(gmap["is_nt"]):
            continue
        target = str(gmap["target"])
        carrier = pd.Series(guide_assignment[:, gi], index=adata.obs_names)
        bg_base = (nt_carrier & ~carrier) if nt_carrier.sum() >= min_background else (~carrier)

        for day, score in guide_level_axes:
            day_mask = obs["day_clean"].eq(day)
            cmask = carrier & day_mask
            bmask = bg_base & day_mask
            if bmask.sum() < min_background:
                bmask = (~carrier) & day_mask

            row = {
                "guide": guide,
                "target": target,
                "day": day,
                "score": score,
                "is_hippo_core": target in present_hippo,
                "n_carrier": int(cmask.sum()),
                "n_background": int(bmask.sum()),
                "carrier_mean": cell_df.loc[cmask, score].mean() if cmask.sum() else np.nan,
                "background_mean": cell_df.loc[bmask, score].mean() if bmask.sum() else np.nan,
                "effect_vs_background": (cell_df.loc[cmask, score].mean() - cell_df.loc[bmask, score].mean()) if cmask.sum() and bmask.sum() else np.nan,
            }
            row.update(compare_distributions(cell_df.loc[cmask, score], cell_df.loc[bmask, score]))
            guide_effect_rows.append(row)

        for day in sorted(obs["day_clean"].unique()):
            day_mask = obs["day_clean"].eq(day)
            cmask = carrier & day_mask
            bmask = bg_base & day_mask
            if bmask.sum() < min_background:
                bmask = (~carrier) & day_mask
            for score in module_names:
                row = {
                    "guide": guide,
                    "target": target,
                    "day": day,
                    "score": score,
                    "is_hippo_core": target in present_hippo,
                    "n_carrier": int(cmask.sum()),
                    "n_background": int(bmask.sum()),
                    "carrier_mean": cell_df.loc[cmask, score].mean() if cmask.sum() else np.nan,
                    "background_mean": cell_df.loc[bmask, score].mean() if bmask.sum() else np.nan,
                    "effect_vs_background": (cell_df.loc[cmask, score].mean() - cell_df.loc[bmask, score].mean()) if cmask.sum() and bmask.sum() else np.nan,
                }
                row.update(compare_distributions(cell_df.loc[cmask, score], cell_df.loc[bmask, score]))
                guide_effect_rows.append(row)

    pd.DataFrame(guide_effect_rows).to_csv(os.path.join(outdir, "guide_level_axis_and_module_effects.csv"), index=False)

    # Target-level distribution tests, including multi-day trajectory axes.
    target_effect_rows = []
    for mode, carrier_dict in modes.items():
        for target, carrier in carrier_dict.items():
            bg_base = (nt_carrier & ~carrier) if nt_carrier.sum() >= min_background else (~carrier)
            for day, score in guide_level_axes:
                day_mask = obs["day_clean"].eq(day)
                cmask = carrier & day_mask
                bmask = bg_base & day_mask
                if bmask.sum() < min_background:
                    bmask = (~carrier) & day_mask
                row = {
                    "mode": mode,
                    "target": target,
                    "day": day,
                    "score": score,
                    "is_hippo_core": target in present_hippo,
                    "n_carrier": int(cmask.sum()),
                    "n_background": int(bmask.sum()),
                    "carrier_mean": cell_df.loc[cmask, score].mean() if cmask.sum() else np.nan,
                    "background_mean": cell_df.loc[bmask, score].mean() if bmask.sum() else np.nan,
                    "effect_vs_background": (cell_df.loc[cmask, score].mean() - cell_df.loc[bmask, score].mean()) if cmask.sum() and bmask.sum() else np.nan,
                }
                row.update(compare_distributions(cell_df.loc[cmask, score], cell_df.loc[bmask, score]))
                target_effect_rows.append(row)

    pd.DataFrame(target_effect_rows).to_csv(os.path.join(outdir, "target_level_fate_axis_distribution_tests.csv"), index=False)

    # Cell-level long distributions for plotting only Hippo targets plus sampled background.
    dist_rows = []
    if present_hippo:
        scores_to_export = [
            "D1_alt_axis_P13_minus_P30",
            "D2_alt_axis_P17_minus_P15",
            "D2_YAPLat_axis_P17P34_minus_P15P23",
            "D3_ALT_axis_P27P34P35_minus_ECcore",
        ] + module_names

        for target in present_hippo:
            carrier = target_carriers[target]
            bg_base = (nt_carrier & ~carrier) if nt_carrier.sum() >= min_background else (~carrier)
            for day in sorted(obs["day_clean"].unique()):
                day_mask = obs["day_clean"].eq(day)
                carrier_ids = cell_df.index[carrier & day_mask].to_numpy()
                bg_ids = cell_df.index[bg_base & day_mask].to_numpy()
                if len(bg_ids) > max_background_cells_per_target_day:
                    bg_ids = rng.choice(bg_ids, size=max_background_cells_per_target_day, replace=False)
                for group_name, ids in [("carrier", carrier_ids), ("background", bg_ids)]:
                    tmp = cell_df.loc[ids, ["cell_id", "day_clean", "cell_type_clean"] + scores_to_export].copy()
                    tmp["target"] = target
                    tmp["carrier_group"] = group_name
                    if embed_key is not None:
                        tmp["embed1"] = cell_df.loc[ids, "embed1"].values
                        tmp["embed2"] = cell_df.loc[ids, "embed2"].values
                        tmp["embedding_key"] = embed_key
                    dist_rows.append(tmp)
        if dist_rows:
            pd.concat(dist_rows, axis=0).to_csv(os.path.join(outdir, "hippo_cell_score_distributions_for_plotting.csv"), index=False)

    # UMAP overlay with carrier status for Hippo.
    if embed_key is not None and present_hippo:
        rows = []
        for target in present_hippo:
            carrier = target_carriers[target]
            bg_base = ~carrier
            for day in sorted(obs["day_clean"].unique()):
                day_mask = obs["day_clean"].eq(day)
                carrier_ids = cell_df.index[carrier & day_mask].to_numpy()
                bg_ids = cell_df.index[bg_base & day_mask].to_numpy()
                if len(bg_ids) > max_background_cells_per_target_day:
                    bg_ids = rng.choice(bg_ids, size=max_background_cells_per_target_day, replace=False)
                ids = np.concatenate([carrier_ids, bg_ids])
                tmp = cell_df.loc[ids, [
                    "cell_id", "day_clean", "cell_type_clean", "embed1", "embed2",
                    "D1_alt_axis_P13_minus_P30", "D2_alt_axis_P17_minus_P15",
                    "D2_YAPLat_axis_P17P34_minus_P15P23", "D3_ALT_axis_P27P34P35_minus_ECcore"
                ]].copy()
                tmp["target"] = target
                tmp["is_carrier"] = tmp.index.isin(carrier_ids)
                rows.append(tmp)
        if rows:
            pd.concat(rows, axis=0).to_csv(os.path.join(outdir, "hippo_umap_overlay_for_plotting.csv"), index=False)

    # Convergence across Hippo guides and targets.
    guide_df = pd.DataFrame(guide_effect_rows)
    if not guide_df.empty:
        hippo_guide = guide_df[guide_df["is_hippo_core"]].copy()
        conv = []
        for score in hippo_guide["score"].dropna().unique():
            sub = hippo_guide[hippo_guide["score"].eq(score)].dropna(subset=["effect_vs_background"])
            if len(sub) == 0:
                continue
            conv.append({
                "score": score,
                "n_hippo_guides": len(sub),
                "n_hippo_targets": sub["target"].nunique(),
                "mean_effect": sub["effect_vs_background"].mean(),
                "median_effect": sub["effect_vs_background"].median(),
                "frac_positive": (sub["effect_vs_background"] > 0).mean(),
                "frac_negative": (sub["effect_vs_background"] < 0).mean(),
            })
        pd.DataFrame(conv).to_csv(os.path.join(outdir, "hippo_guide_convergence_summary.csv"), index=False)

    # Manifest
    manifest = pd.DataFrame({
        "file": sorted(os.listdir(outdir))
    })
    manifest.to_csv(os.path.join(outdir, "export_manifest.csv"), index=False)

    print("Wrote outputs to:", outdir)
    print("Hippo targets present:", ", ".join(present_hippo) if present_hippo else "none")
    print("Key files:")
    for fn in [
        "composition_counts_by_target_and_guide.csv",
        "state_enrichment_target_level.csv",
        "guide_level_axis_and_module_effects.csv",
        "target_level_fate_axis_distribution_tests.csv",
        "hippo_cell_score_distributions_for_plotting.csv",
        "hippo_umap_overlay_for_plotting.csv",
        "hippo_guide_convergence_summary.csv",
    ]:
        if os.path.exists(os.path.join(outdir, fn)):
            print(" -", os.path.join(outdir, fn))

    return {
        "outdir": outdir,
        "hippo_targets_present": present_hippo,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata", required=True)
    parser.add_argument("--guide_target_csv", default="combined_TargetKD.csv")
    parser.add_argument("--outdir", default="hippo_nature_style_exports")
    parser.add_argument("--cnmf_key", default="cnmf_usage")
    parser.add_argument("--guide_assignment_key", default="guide_assignment")
    parser.add_argument("--day_col", default="day")
    parser.add_argument("--celltype_col", default="cell_type")
    parser.add_argument("--use_raw_for_modules", action="store_true")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.adata)

    export_hippo_nature_style_cell_csvs(
        adata,
        guide_target_csv=args.guide_target_csv,
        outdir=args.outdir,
        cnmf_key=args.cnmf_key,
        guide_assignment_key=args.guide_assignment_key,
        day_col=args.day_col,
        celltype_col=args.celltype_col,
        use_raw_for_modules=args.use_raw_for_modules,
    )


if __name__ == "__main__":
    main()
