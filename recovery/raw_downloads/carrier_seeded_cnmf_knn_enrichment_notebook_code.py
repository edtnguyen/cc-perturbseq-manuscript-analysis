
# ============================================================
# Carrier-seeded cNMF-kNN local enrichment map
# ============================================================
#
# Goal:
#   For each Hippo gene, use carrier cells of that gene as the neighborhood seeds.
#   Build kNN in cNMF program-usage space.
#   For each carrier seed cell, compute local enrichment of the same Hippo target
#   versus non-targeting controls in that seed's neighborhood.
#   Plot one UMAP per gene, all days pooled, with carrier cells colored by local log2 enrichment.
#
# Important:
#   This is a carrier-seeded local enrichment map, not formal Milo DA.
#   Keep your pertpy Milo output as the formal DA test if needed.
#
# Interpretation:
#   Positive score = the carrier seed lies in a cNMF-neighborhood locally enriched
#                    for that Hippo target relative to non-targeting controls.
#   Negative score = local depletion.
#
# ============================================================

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy import sparse
from scipy.stats import fisher_exact
from sklearn.neighbors import NearestNeighbors


# ============================================================
# Settings
# ============================================================

OUTDIR = Path("carrier_seeded_cnmf_knn_enrichment")
OUTDIR.mkdir(parents=True, exist_ok=True)

HIPPO_TARGETS = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]

DAY_COL = "day"
CELLTYPE_COL = "cell_type"
GUIDE_ASSIGNMENT_KEY = "guide_assignment"
CNMF_KEY = "cnmf_usage"
CONTROL_TARGET = "non-targeting"

# If adata.obsm["cnmf_usage"] is a NumPy array:
#   0 => Program 0, Program 1, ..., Program 49
#   1 => Program 1, Program 2, ..., Program 50
# Use the same value that worked in your other scripts.
CNMF_ARRAY_PROGRAM_INDEX_BASE = 0

# Use all available cNMF programs by default.
# To exclude technical programs, set this to a list like [2,3,8,11,13,15,17,23,27,30,34,35]
PROGRAMS_FOR_KNN = None

# cNMF transform for KNN space:
#   "zscore" = column z-score cNMF usages
#   "clr_zscore" = row centered log-ratio transform, then column z-score
CNMF_TRANSFORM = "zscore"

# Neighborhood parameters.
K = 80
PSEUDOCOUNT = 0.5

# Recommended for developmental time series:
#   True  = a D2 seed uses only D2 neighbors and D2 target/NT background
#   False = a seed can use neighbors from all days, but still uses same-day baseline if SAME_DAY_BASELINE=True
RESTRICT_NEIGHBORS_TO_SAME_DAY = True
SAME_DAY_BASELINE = True

# For same-day neighbor mode, search more than K cells first, then keep the closest same-day K.
K_SEARCH_MULTIPLIER = 8

# High-MOI sensitivity mode:
#   "all_carriers"
#   "low_guide_load_leq2"
#   "single_selected_target"
MODE = "all_carriers"

# UMAP plotting.
UMAP_KEY = "X_umap"
BACKGROUND_POINT_SIZE = 0.8
BACKGROUND_ALPHA = 0.10
CARRIER_POINT_SIZE = 7.0
CARRIER_ALPHA = 0.95
LOG2FC_CLIP = 3.0

# Optional: plot only selected days on background/carriers.
# None = all days.
DAYS_TO_USE = None
# DAYS_TO_USE = ["d1", "d2", "d3"]

RANDOM_SEED = 1


# ============================================================
# Guide and target helpers
# ============================================================

def looks_nt(x):
    s = str(x).lower()
    patterns = [
        "non-target", "nontarget", "non_target", "non targeting",
        "negative", "control", "scramble", "scrambled", "safe", "ntc",
    ]
    return any(p in s for p in patterns) or s in {"nt", "neg", "ctrl"}


def parse_target_from_guide(guide):
    g = str(guide)
    if looks_nt(g):
        return CONTROL_TARGET
    for sep in ["|", ":", "_", "-", ".", " "]:
        if sep in g:
            first = g.split(sep)[0]
            if first:
                return first
    return g


def infer_col(df, candidates):
    lower = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    for col in df.columns:
        cl = str(col).lower()
        for c in candidates:
            if c.lower() in cl:
                return col
    return None


def get_guide_names(adata, guide_assignment_key=GUIDE_ASSIGNMENT_KEY):
    ga = adata.obsm[guide_assignment_key]
    if isinstance(ga, pd.DataFrame):
        return np.asarray(ga.columns.astype(str))
    if "guide_names" in adata.uns:
        return np.asarray(adata.uns["guide_names"]).astype(str)
    raise KeyError("Could not find guide names in guide_assignment columns or adata.uns['guide_names'].")


def get_guide_assignment_df(adata, guide_assignment_key=GUIDE_ASSIGNMENT_KEY):
    ga = adata.obsm[guide_assignment_key]
    guide_names = get_guide_names(adata, guide_assignment_key)

    if isinstance(ga, pd.DataFrame):
        out = ga.copy()
        out.index = adata.obs_names
        out.columns = out.columns.astype(str)
        return out.astype(bool)

    if sparse.issparse(ga):
        ga = ga.toarray()

    return pd.DataFrame(
        np.asarray(ga),
        index=adata.obs_names,
        columns=guide_names,
    ).astype(bool)


def build_guide_target_table(adata, guide_target_csv=None):
    guide_names = pd.Index(get_guide_names(adata), name="guide").astype(str)
    base = pd.DataFrame({"guide": guide_names})

    if guide_target_csv is not None and os.path.exists(guide_target_csv):
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
    out.loc[out["is_non_targeting"], "target"] = CONTROL_TARGET
    return out


def get_target_assignment_df(guide_assignment, guide_target):
    guide_assignment = guide_assignment.astype(bool)

    gt = guide_target[["guide", "target"]].copy()
    gt["guide"] = gt["guide"].astype(str)
    gt["target"] = gt["target"].astype(str)
    gt = gt[gt["guide"].isin(guide_assignment.columns)].copy()

    series_by_target = {}
    for target, sub in gt.groupby("target", sort=True):
        guides = sub["guide"].tolist()
        if len(guides) == 1:
            series_by_target[target] = guide_assignment[guides[0]].astype(bool)
        else:
            series_by_target[target] = guide_assignment.loc[:, guides].any(axis=1)

    return pd.DataFrame(series_by_target, index=guide_assignment.index).astype(bool).copy()


# ============================================================
# cNMF helpers
# ============================================================

def standardize_program_col(c):
    s = str(c)
    if s.startswith("Program "):
        return f"Program {int(s.replace('Program ', ''))}"
    if s.startswith("P"):
        return f"Program {int(s.replace('P', ''))}"
    return f"Program {int(float(s))}"


def get_program_score_df(adata, cnmf_key=CNMF_KEY, array_program_index_base=CNMF_ARRAY_PROGRAM_INDEX_BASE):
    X = adata.obsm[cnmf_key]

    if isinstance(X, pd.DataFrame):
        df = X.copy()
        df.index = adata.obs_names
        df.columns = [standardize_program_col(c) for c in df.columns]

        if pd.Index(df.columns).duplicated().any():
            print("[warning] Duplicate cNMF program columns after standardization; collapsing by mean.")
            df = df.T.groupby(level=0).mean().T

        return df.astype(float).copy()

    if sparse.issparse(X):
        X = X.toarray()

    X = np.asarray(X)
    cols = [f"Program {i + array_program_index_base}" for i in range(X.shape[1])]
    return pd.DataFrame(X, index=adata.obs_names, columns=cols).astype(float).copy()


def make_cnmf_knn_matrix(program_df, programs_for_knn=None, transform="zscore"):
    """
    Returns X_knn and the program columns used.
    """
    if programs_for_knn is None:
        cols = [c for c in program_df.columns if str(c).startswith("Program ")]
    else:
        cols = [f"Program {int(p)}" for p in programs_for_knn]

    missing = [c for c in cols if c not in program_df.columns]
    if missing:
        raise KeyError(f"Missing requested cNMF programs: {missing}")

    X = program_df[cols].copy()
    X = X.apply(pd.to_numeric, errors="coerce")
    X = X.fillna(0.0)

    # Remove zero-variance columns.
    sd = X.std(axis=0, ddof=1)
    cols = sd[sd > 0].index.tolist()
    X = X[cols]

    if transform == "zscore":
        X_mat = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1).replace(0, 1.0)
        X_mat = X_mat.fillna(0.0).values

    elif transform == "clr_zscore":
        # Safe CLR transform for nonnegative program usages.
        # Add a small pseudocount, log, subtract row mean, then column z-score.
        X_nonneg = X.clip(lower=0)
        pc = np.nanpercentile(X_nonneg.values[X_nonneg.values > 0], 1) if np.any(X_nonneg.values > 0) else 1e-6
        pc = max(pc, 1e-6)

        logX = np.log(X_nonneg + pc)
        clr = logX.sub(logX.mean(axis=1), axis=0)
        X_mat = (clr - clr.mean(axis=0)) / clr.std(axis=0, ddof=1).replace(0, 1.0)
        X_mat = X_mat.fillna(0.0).values

    else:
        raise ValueError("transform must be 'zscore' or 'clr_zscore'.")

    return X_mat.astype(np.float32), cols


# ============================================================
# Build base tables
# ============================================================

def prepare_carrier_seeded_inputs(adata, guide_target_csv=None):
    adata_use = adata.copy()

    if not adata_use.obs_names.is_unique:
        adata_use.obs_names_make_unique()

    adata_use.obs["day_clean"] = adata_use.obs[DAY_COL].astype(str).str.lower()

    if CELLTYPE_COL in adata_use.obs.columns:
        adata_use.obs["cell_type_clean"] = adata_use.obs[CELLTYPE_COL].astype(str)
    else:
        adata_use.obs["cell_type_clean"] = "unknown"

    guide_assignment = get_guide_assignment_df(adata_use)
    guide_target = build_guide_target_table(adata_use, guide_target_csv=guide_target_csv)
    target_assignment = get_target_assignment_df(guide_assignment, guide_target)

    if CONTROL_TARGET not in target_assignment.columns:
        raise KeyError(f"Could not find {CONTROL_TARGET!r} in target assignment.")

    targeting_cols = [c for c in target_assignment.columns if c != CONTROL_TARGET]

    adata_use.obs["n_total_guides"] = guide_assignment.sum(axis=1).astype(int).values
    adata_use.obs["n_targeting_targets"] = target_assignment[targeting_cols].sum(axis=1).astype(int).values

    program_df = get_program_score_df(adata_use)

    X_knn, program_cols_used = make_cnmf_knn_matrix(
        program_df=program_df,
        programs_for_knn=PROGRAMS_FOR_KNN,
        transform=CNMF_TRANSFORM,
    )

    if UMAP_KEY not in adata_use.obsm:
        raise KeyError(f"adata.obsm[{UMAP_KEY!r}] not found.")

    umap = np.asarray(adata_use.obsm[UMAP_KEY])
    adata_use.obs["UMAP1"] = umap[:, 0]
    adata_use.obs["UMAP2"] = umap[:, 1]

    # Optional day filtering.
    if DAYS_TO_USE is not None:
        day_mask = adata_use.obs["day_clean"].isin([str(d).lower() for d in DAYS_TO_USE]).values
    else:
        day_mask = np.ones(adata_use.n_obs, dtype=bool)

    print(f"Using {len(program_cols_used)} cNMF programs for KNN.")
    print(f"KNN transform: {CNMF_TRANSFORM}")
    print(f"Mode: {MODE}")
    print(f"Restrict neighbors to same day: {RESTRICT_NEIGHBORS_TO_SAME_DAY}")
    print(f"Same-day baseline: {SAME_DAY_BASELINE}")

    guide_target.to_csv(OUTDIR / "guide_to_target_map.csv", index=False)

    target_counts = pd.DataFrame({
        "target": target_assignment.columns,
        "n_cells": target_assignment.sum(axis=0).values,
    })
    target_counts.to_csv(OUTDIR / "target_carrier_counts.csv", index=False)

    pd.DataFrame({"program_column_used": program_cols_used}).to_csv(
        OUTDIR / "cnmf_programs_used_for_knn.csv",
        index=False,
    )

    return adata_use, target_assignment, X_knn, day_mask


# ============================================================
# Carrier-seeded local enrichment
# ============================================================

def get_masks_for_target(adata, target_assignment, target, mode="all_carriers", day_mask=None):
    if target not in target_assignment.columns:
        raise KeyError(f"Target {target} not found.")

    if day_mask is None:
        day_mask = np.ones(adata.n_obs, dtype=bool)

    target_mask_raw = target_assignment[target].values
    nt_mask_raw = target_assignment[CONTROL_TARGET].values

    if mode == "all_carriers":
        target_mask = target_mask_raw & day_mask
        nt_mask = nt_mask_raw & (~target_mask_raw) & day_mask

    elif mode == "low_guide_load_leq2":
        low_load = adata.obs["n_total_guides"].values <= 2
        target_mask = target_mask_raw & low_load & day_mask
        nt_mask = nt_mask_raw & (~target_mask_raw) & low_load & day_mask

    elif mode == "single_selected_target":
        n_targeting = adata.obs["n_targeting_targets"].values
        target_mask = target_mask_raw & (n_targeting == 1) & day_mask
        nt_mask = nt_mask_raw & (~target_mask_raw) & (n_targeting == 0) & day_mask

    else:
        raise ValueError(f"Unknown mode: {mode}")

    return target_mask, nt_mask


def compute_carrier_seeded_enrichment(
    adata,
    target_assignment,
    X_knn,
    target,
    mode="all_carriers",
    k=K,
    day_mask=None,
    restrict_neighbors_to_same_day=RESTRICT_NEIGHBORS_TO_SAME_DAY,
    same_day_baseline=SAME_DAY_BASELINE,
):
    """
    For each carrier cell of `target`, compute local target-vs-NT enrichment in cNMF-kNN space.
    """
    target_mask, nt_mask = get_masks_for_target(
        adata=adata,
        target_assignment=target_assignment,
        target=target,
        mode=mode,
        day_mask=day_mask,
    )

    seed_idx = np.where(target_mask)[0]

    if len(seed_idx) == 0:
        print(f"[skip] {target}: no carrier seeds for mode={mode}")
        return pd.DataFrame()

    day_vec = adata.obs["day_clean"].astype(str).values
    cell_type_vec = adata.obs["cell_type_clean"].astype(str).values

    if restrict_neighbors_to_same_day:
        k_search = min(max(k * K_SEARCH_MULTIPLIER, k + 1), X_knn.shape[0])
    else:
        k_search = min(k + 1, X_knn.shape[0])

    nbrs = NearestNeighbors(n_neighbors=k_search, metric="euclidean")
    nbrs.fit(X_knn)

    neigh_ind = nbrs.kneighbors(X_knn[seed_idx], return_distance=False)

    rows = []

    for seed_cell_idx, neigh in zip(seed_idx, neigh_ind):
        seed_day = day_vec[seed_cell_idx]

        # Remove the seed itself.
        neigh = neigh[neigh != seed_cell_idx]

        if restrict_neighbors_to_same_day:
            neigh = neigh[day_vec[neigh] == seed_day]

        neigh = neigh[:k]

        if same_day_baseline:
            bg_mask = (day_vec == seed_day)
            bg_target_total = int(np.sum(target_mask & bg_mask)) - 1  # exclude seed from target background
            bg_nt_total = int(np.sum(nt_mask & bg_mask))
        else:
            bg_target_total = int(np.sum(target_mask)) - 1
            bg_nt_total = int(np.sum(nt_mask))

        bg_target_total = max(bg_target_total, 0)

        local_target_n = int(np.sum(target_mask[neigh]))
        local_nt_n = int(np.sum(nt_mask[neigh]))

        # Local ratio relative to same-day/global background ratio.
        local_ratio = (local_target_n + PSEUDOCOUNT) / (local_nt_n + PSEUDOCOUNT)
        bg_ratio = (bg_target_total + PSEUDOCOUNT) / (bg_nt_total + PSEUDOCOUNT)

        log2_local_enrichment = np.log2(local_ratio / bg_ratio)

        # Fisher test: local neighborhood vs same-day/global background outside the neighborhood.
        bg_target_outside = max(bg_target_total - local_target_n, 0)
        bg_nt_outside = max(bg_nt_total - local_nt_n, 0)

        try:
            fisher_p = fisher_exact(
                [
                    [local_target_n, local_nt_n],
                    [bg_target_outside, bg_nt_outside],
                ],
                alternative="two-sided",
            )[1]
        except Exception:
            fisher_p = np.nan

        rows.append({
            "target": target,
            "seed_cell_idx": int(seed_cell_idx),
            "seed_cell_id": adata.obs_names[seed_cell_idx],
            "day": seed_day,
            "cell_type": cell_type_vec[seed_cell_idx],
            "mode": mode,
            "k_requested": k,
            "k_used": int(len(neigh)),
            "restrict_neighbors_to_same_day": bool(restrict_neighbors_to_same_day),
            "same_day_baseline": bool(same_day_baseline),
            "local_target_n": local_target_n,
            "local_nt_n": local_nt_n,
            "background_target_total": int(bg_target_total),
            "background_nt_total": int(bg_nt_total),
            "local_target_fraction_among_target_nt": (
                local_target_n / (local_target_n + local_nt_n)
                if (local_target_n + local_nt_n) > 0
                else np.nan
            ),
            "background_target_fraction_among_target_nt": (
                bg_target_total / (bg_target_total + bg_nt_total)
                if (bg_target_total + bg_nt_total) > 0
                else np.nan
            ),
            "log2_local_enrichment": log2_local_enrichment,
            "fisher_p": fisher_p,
            "UMAP1": float(adata.obs["UMAP1"].iloc[seed_cell_idx]),
            "UMAP2": float(adata.obs["UMAP2"].iloc[seed_cell_idx]),
            "n_total_guides": int(adata.obs["n_total_guides"].iloc[seed_cell_idx]),
            "n_targeting_targets": int(adata.obs["n_targeting_targets"].iloc[seed_cell_idx]),
        })

    out = pd.DataFrame(rows)

    # FDR within each target/mode.
    if len(out):
        p = out["fisher_p"].values.astype(float)
        order = np.argsort(np.nan_to_num(p, nan=1.0))
        q = np.full_like(p, np.nan, dtype=float)
        valid = np.isfinite(p)
        pv = p[valid]
        if len(pv):
            idx = np.where(valid)[0]
            order_valid = idx[np.argsort(pv)]
            ranks = np.arange(1, len(order_valid) + 1)
            q_sorted = p[order_valid] * len(order_valid) / ranks
            q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
            q[order_valid] = np.clip(q_sorted, 0, 1)
        out["fisher_fdr"] = q

    return out


# ============================================================
# Plotting
# ============================================================

def plot_carrier_seeded_umap(
    adata,
    score_df,
    target,
    outdir=OUTDIR,
    log2fc_clip=LOG2FC_CLIP,
):
    """
    One UMAP per gene, all days pooled.
    Background = all cells.
    Colored points = carrier seed cells.
    No outline around points.
    """
    sub = score_df[score_df["target"].eq(target)].copy()

    if len(sub) == 0:
        print(f"No scores for {target}")
        return None

    bg = adata.obs[["UMAP1", "UMAP2"]].copy()

    fig, ax = plt.subplots(figsize=(6.8, 6.2))

    ax.scatter(
        bg["UMAP1"],
        bg["UMAP2"],
        s=BACKGROUND_POINT_SIZE,
        c="lightgrey",
        alpha=BACKGROUND_ALPHA,
        linewidths=0,
        rasterized=True,
    )

    norm = TwoSlopeNorm(
        vmin=-log2fc_clip,
        vcenter=0,
        vmax=log2fc_clip,
    )

    sc = ax.scatter(
        sub["UMAP1"],
        sub["UMAP2"],
        c=sub["log2_local_enrichment"].clip(-log2fc_clip, log2fc_clip),
        cmap="RdBu_r",
        norm=norm,
        s=CARRIER_POINT_SIZE,
        alpha=CARRIER_ALPHA,
        linewidths=0,
        edgecolors="none",
        rasterized=True,
    )

    n_seed = len(sub)
    n_pos = int((sub["log2_local_enrichment"] > 0).sum())
    n_neg = int((sub["log2_local_enrichment"] < 0).sum())
    n_fdr10 = int((sub["fisher_fdr"] < 0.10).sum()) if "fisher_fdr" in sub.columns else 0

    med = sub["log2_local_enrichment"].median()
    mean = sub["log2_local_enrichment"].mean()

    ax.set_title(
        f"{target} carrier-seeded cNMF-kNN enrichment\n"
        f"{MODE}; n={n_seed} seeds; median={med:.2f}; FDR<0.10={n_fdr10}",
        fontsize=12,
        fontweight="bold",
    )

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")

    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("Local enrichment log2FC\nTarget vs NT", fontsize=10)

    fig.tight_layout()

    out_png = outdir / f"carrier_seeded_cnmf_knn_UMAP_{target}_{MODE}.png"
    out_pdf = outdir / f"carrier_seeded_cnmf_knn_UMAP_{target}_{MODE}.pdf"

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.show()

    print("Saved:", out_png)

    return fig


# ============================================================
# Main runner
# ============================================================

def run_carrier_seeded_cnmf_knn_enrichment(
    adata,
    guide_target_csv="combined_TargetKD.csv",
    targets=HIPPO_TARGETS,
):
    adata_use, target_assignment, X_knn, day_mask = prepare_carrier_seeded_inputs(
        adata,
        guide_target_csv=guide_target_csv,
    )

    present_targets = [t for t in targets if t in target_assignment.columns]
    missing_targets = [t for t in targets if t not in target_assignment.columns]

    print("Present targets:", present_targets)
    if missing_targets:
        print("Missing targets:", missing_targets)

    all_scores = []

    for target in present_targets:
        print(f"\nComputing carrier-seeded enrichment for {target}")
        scores = compute_carrier_seeded_enrichment(
            adata=adata_use,
            target_assignment=target_assignment,
            X_knn=X_knn,
            target=target,
            mode=MODE,
            k=K,
            day_mask=day_mask,
            restrict_neighbors_to_same_day=RESTRICT_NEIGHBORS_TO_SAME_DAY,
            same_day_baseline=SAME_DAY_BASELINE,
        )

        if len(scores):
            scores.to_csv(
                OUTDIR / f"carrier_seeded_scores_{target}_{MODE}.csv",
                index=False,
            )
            all_scores.append(scores)

    if not all_scores:
        raise RuntimeError("No carrier-seeded scores were generated.")

    combined = pd.concat(all_scores, ignore_index=True)
    combined.to_csv(
        OUTDIR / f"combined_carrier_seeded_cnmf_knn_scores_{MODE}.csv",
        index=False,
    )

    summary = combined.groupby("target").agg(
        n_seed_cells=("seed_cell_id", "count"),
        median_log2_enrichment=("log2_local_enrichment", "median"),
        mean_log2_enrichment=("log2_local_enrichment", "mean"),
        frac_positive=("log2_local_enrichment", lambda x: np.mean(x > 0)),
        n_fdr10=("fisher_fdr", lambda x: int(np.sum(x < 0.10))),
        n_fdr05=("fisher_fdr", lambda x: int(np.sum(x < 0.05))),
        median_k_used=("k_used", "median"),
    ).reset_index()

    summary.to_csv(
        OUTDIR / f"summary_carrier_seeded_cnmf_knn_{MODE}.csv",
        index=False,
    )

    print("\nSummary:")
    print(summary.to_string(index=False))

    for target in present_targets:
        plot_carrier_seeded_umap(
            adata=adata_use,
            score_df=combined,
            target=target,
        )

    return combined, summary


# ============================================================
# Run this cell
# ============================================================

# Example:
# combined_scores, enrichment_summary = run_carrier_seeded_cnmf_knn_enrichment(
#     adata,
#     guide_target_csv="combined_TargetKD.csv",
#     targets=HIPPO_TARGETS,
# )
