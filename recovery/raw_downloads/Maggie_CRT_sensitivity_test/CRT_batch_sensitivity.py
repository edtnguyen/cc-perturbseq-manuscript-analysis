#!/usr/bin/env python3
"""CRT wrapper for inferred-batch sensitivity analyses.

Differences from the original wrapper:
1. Requested covariates/categorical keys may come from either the cNMF modality
   obs or the guide AnnData obs, aligned by cell barcode.
2. Supports a calibration-only mode that runs NTC pseudo-target CRT tests but
   skips discovery tests over all real genes.
3. Exposes NTC ensemble sizes and seeds as CLI arguments.

The default (without --calibration-only) retains the original full CRT behavior.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
import yaml

from src.sceptre import (
    build_ntc_group_inputs,
    compute_guide_set_null_pvals,
    crt_pvals_for_ntc_groups_ensemble,
    crt_pvals_for_ntc_groups_ensemble_skew,
    make_ntc_groups_ensemble,
    prepare_crt_inputs,
    run_all_genes_union_crt,
)
from src.visualization import qq_plot_ntc_pvals


def _safe_name(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "condition"


def _resolve_obs_column(
    key: str,
    primary_obs: pd.DataFrame,
    guide_obs: pd.DataFrame,
    cell_index: pd.Index,
) -> pd.Series:
    if key in primary_obs.columns:
        s = primary_obs.loc[cell_index, key]
    elif key in guide_obs.columns:
        missing = cell_index.difference(guide_obs.index)
        if len(missing):
            raise KeyError(
                f"{len(missing)} cNMF cells are missing from guide obs while resolving {key!r}"
            )
        s = guide_obs.loc[cell_index, key]
    else:
        raise KeyError(
            f"Requested obs key {key!r} was not found in either cNMF obs or guide obs"
        )
    return pd.Series(s.to_numpy(), index=cell_index, name=key)


def reformat_data_for_crt(
    mdata,
    guide_data,
    *,
    covariates: list[str] | None,
    log_covariates: list[str] | None,
    categorical_key: str,
):
    adata = mdata["cNMF"].copy()
    adata.obsm["cnmf_usage"] = np.asarray(adata.X)

    cells = pd.Index(adata.obs_names)
    g = guide_data[cells].copy()
    guide_obs = g.obs

    # Pull categorical key from the guide object if it was added there.
    if categorical_key not in adata.obs:
        adata.obs[categorical_key] = _resolve_obs_column(
            categorical_key, adata.obs, guide_obs, cells
        ).to_numpy()

    adata.obsm["guide_assignment"] = g.obsm["guide_assignment"].copy()
    adata.obs["guides_per_cell"] = np.asarray(
        adata.obsm["guide_assignment"].sum(axis=1)
    ).ravel()
    adata.uns["program_names"] = list(adata.var_names)

    guide_names = list(guide_data.uns["guide_names"])
    adata.uns["guide_names"] = guide_names
    adata.uns["guide2gene"] = dict(
        zip(guide_names, guide_data.uns["guide_targets"])
    )

    covar_dict: dict[str, pd.Series] = {}
    for key in covariates or []:
        covar_dict[key] = _resolve_obs_column(
            key, adata.obs, guide_obs, cells
        )
    for key in log_covariates or []:
        values = pd.to_numeric(
            _resolve_obs_column(key, adata.obs, guide_obs, cells), errors="coerce"
        )
        if values.isna().any():
            raise ValueError(f"log covariate {key!r} has missing/non-numeric values")
        if (values < 0).any():
            raise ValueError(f"log covariate {key!r} has negative values")
        covar_dict[f"log_{key}"] = np.log1p(values)

    adata.obsm["covar"] = pd.DataFrame(covar_dict, index=cells)
    return adata


def _flatten_ntc_frames(ntc_frames: dict) -> np.ndarray:
    vals = []
    for df in ntc_frames.values():
        arr = df.to_numpy(dtype=float, copy=False).ravel()
        vals.append(arr[np.isfinite(arr)])
    return np.concatenate(vals) if vals else np.empty(0, dtype=float)


def _qq_xy(pvals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    p = np.asarray(pvals, dtype=float)
    p = p[np.isfinite(p)]
    p = np.clip(p, 1e-300, 1.0)
    p.sort()
    m = p.size
    if m == 0:
        return np.empty(0), np.empty(0)
    expected = (np.arange(1, m + 1) - 0.5) / m
    return -np.log10(expected), -np.log10(p)


def plot_calibration_only(
    ntc_raw: dict,
    null_pvals: np.ndarray,
    *,
    title: str,
    output_path: Path,
) -> None:
    ntc = _flatten_ntc_frames(ntc_raw)
    x_ntc, y_ntc = _qq_xy(ntc)
    x_null, y_null = _qq_xy(null_pvals)

    plt.figure(figsize=(7, 6))
    if x_ntc.size:
        plt.plot(x_ntc, y_ntc, label="NTC pseudo-target CRT p-values", linewidth=1.5)
    if x_null.size:
        plt.plot(x_null, y_null, "--", label="CRT-null leave-one-out p-values", linewidth=1.3)
    lim = max(
        1.0,
        float(np.nanmax(x_ntc)) if x_ntc.size else 0.0,
        float(np.nanmax(y_ntc)) if y_ntc.size else 0.0,
        float(np.nanmax(x_null)) if x_null.size else 0.0,
        float(np.nanmax(y_null)) if y_null.size else 0.0,
    )
    plt.plot([0, lim], [0, lim], color="black", linewidth=1, label="y=x")
    plt.xlim(0, lim)
    plt.ylim(0, lim)
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.title(title)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def save_ntc_artifacts(
    *,
    ntc_group_pvals_ens: dict,
    null_pvals: np.ndarray,
    ntc_groups_ens: list,
    output_folder: Path,
    k: int,
    condition: str,
) -> None:
    safe = _safe_name(condition)
    frames = []
    for ensemble, df in ntc_group_pvals_ens.items():
        df2 = df.copy()
        df2.index.name = "group_id"
        df2 = df2.reset_index()
        df2.insert(0, "ensemble", ensemble)
        frames.append(df2)
    pd.concat(frames, ignore_index=True).to_csv(
        output_folder / f"{k}_CRT_{safe}_ntc_null_pvals.txt",
        sep="\t",
        index=False,
    )

    pd.DataFrame({"null_pval": null_pvals}).to_csv(
        output_folder / f"{k}_CRT_{safe}_permutation_null_pvals.txt",
        sep="\t",
        index=False,
    )

    rows = []
    for ensemble, groups in enumerate(ntc_groups_ens):
        for group_id, guides in groups.items():
            for guide in guides:
                rows.append(
                    {
                        "ensemble": ensemble,
                        "group_id": group_id,
                        "guide": guide,
                    }
                )
    pd.DataFrame(rows).to_csv(
        output_folder / f"{k}_CRT_{safe}_ntc_group_guides.txt",
        sep="\t",
        index=False,
    )


def save_full_results(
    out: dict,
    *,
    k: int,
    output_folder: Path,
    condition: str,
    fdr_method: str,
) -> None:
    safe = _safe_name(condition)
    pval_df = out["pvals_raw_df"]
    beta_df = out["betas_df"]

    pval_long = (
        pval_df.reset_index()
        .melt(id_vars="index", var_name="program_name", value_name="p-value")
        .rename(columns={"index": "target_name"})
    )
    beta_long = (
        beta_df.reset_index()
        .melt(id_vars="index", var_name="program_name", value_name="log2FC")
        .rename(columns={"index": "target_name"})
    )
    result = pval_long.merge(beta_long, on=["program_name", "target_name"])
    result = result[["target_name", "program_name", "log2FC", "p-value"]]

    if fdr_method == "BH":
        from statsmodels.stats.multitest import multipletests
        result["adj_pval"] = multipletests(result["p-value"], method="fdr_bh")[1]
    else:
        from multipy.fdr import qvalue
        import builtins
        if not hasattr(builtins, "xrange"):
            builtins.xrange = range
        result["adj_pval"] = qvalue(
            result["p-value"].to_numpy(), threshold=0.05, verbose=False
        )[1]

    result.to_csv(
        output_folder / f"{k}_CRT_{safe}.txt", sep="\t", index=False
    )


def run_condition(
    adata,
    *,
    k: int,
    condition: object,
    output_folder: Path,
    args: argparse.Namespace,
) -> None:
    mask = adata.obs[args.categorical_key].astype(str) == str(condition)
    adata_con = adata[mask].copy()
    safe = _safe_name(condition)
    print(f"Condition {condition!r}: {adata_con.n_obs:,} cells")

    if adata_con.n_obs < args.min_cells_per_condition:
        print(f"Skipping {condition!r}: fewer than {args.min_cells_per_condition} cells")
        return

    inputs = prepare_crt_inputs(
        adata=adata_con,
        usage_key="cnmf_usage",
        covar_key="covar",
        guide_assignment_key="guide_assignment",
        guide2gene_key="guide2gene",
    )

    out = None
    if not args.calibration_only:
        out = run_all_genes_union_crt(
            inputs=inputs,
            B=args.number_permutations,
            n_jobs=args.n_jobs,
            calibrate_skew_normal=args.compute_skew,
            return_raw_pvals=True,
            return_skew_normal=args.compute_skew,
        )

    ntc_labels = args.guide_annotation_key
    ntc_guides, guide_freq, guide_to_bin, real_sigs = build_ntc_group_inputs(
        inputs=inputs,
        ntc_label=ntc_labels,
        group_size=args.number_guide,
        n_bins=args.ntc_frequency_bins,
    )
    ntc_groups_ens = make_ntc_groups_ensemble(
        ntc_guides=ntc_guides,
        ntc_freq=guide_freq,
        real_gene_bin_sigs=real_sigs,
        guide_to_bin=guide_to_bin,
        n_ensemble=args.ntc_ensembles,
        seed0=args.ntc_group_seed,
        group_size=args.number_guide,
        max_groups=args.ntc_max_groups,
        replace=True,
    )
    ntc_raw = crt_pvals_for_ntc_groups_ensemble(
        inputs=inputs,
        ntc_groups_ens=ntc_groups_ens,
        B=args.number_permutations,
        seed0=args.ntc_crt_seed,
    )

    ntc_skew = None
    if args.compute_skew:
        ntc_skew = crt_pvals_for_ntc_groups_ensemble_skew(
            inputs=inputs,
            ntc_groups_ens=ntc_groups_ens,
            B=args.number_permutations,
            seed0=args.ntc_crt_seed,
        )

    guide_to_col = {g: i for i, g in enumerate(inputs.guide_names)}
    null_chunks = []
    for groups in ntc_groups_ens:
        for guides in groups.values():
            null_chunks.append(
                compute_guide_set_null_pvals(
                    guide_idx=[guide_to_col[g] for g in guides],
                    inputs=inputs,
                    B=args.number_permutations,
                ).ravel()
            )
    null_pvals = np.concatenate(null_chunks) if null_chunks else np.empty(0)

    save_ntc_artifacts(
        ntc_group_pvals_ens=ntc_raw,
        null_pvals=null_pvals,
        ntc_groups_ens=ntc_groups_ens,
        output_folder=output_folder,
        k=k,
        condition=str(condition),
    )

    if args.calibration_only:
        plot_calibration_only(
            ntc_raw,
            null_pvals,
            title=f"NTC calibration sensitivity: {condition}",
            output_path=output_folder / f"CRT_{safe}_calibration_only.png",
        )
        return

    assert out is not None
    ax = qq_plot_ntc_pvals(
        pvals_raw_df=out["pvals_raw_df"],
        guide2gene=adata.uns["guide2gene"],
        ntc_genes=ntc_labels,
        pvals_skew_df=out.get("pvals_df") if args.compute_skew else None,
        null_pvals=null_pvals,
        ntc_group_pvals_ens=ntc_raw,
        ntc_group_pvals_skew_ens=ntc_skew,
        show_ntc_ensemble_band=True,
        show_all_pvals=True,
        title=f"QQ: grouped NTC controls vs CRT null for {condition}",
    )
    plt.tight_layout()
    plt.savefig(output_folder / f"CRT_{safe}.png", dpi=150)
    plt.close()

    save_full_results(
        out,
        k=k,
        output_folder=output_folder,
        condition=str(condition),
        fdr_method=args.FDR_method,
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--out_dir", required=True)
    p.add_argument("--run_name", required=True)
    p.add_argument("--components", nargs="*", type=int, default=None)
    p.add_argument("--sel_thresh", nargs="*", type=float, default=None)
    p.add_argument("--mdata_guide_path", required=True)
    p.add_argument("--categorical_key", default="sample")
    p.add_argument("--covariates", nargs="*", default=None)
    p.add_argument("--log_covariates", nargs="*", default=None)
    p.add_argument("--number_guide", type=int, default=6)
    p.add_argument("--number_permutations", type=int, default=1024)
    p.add_argument("--guide_annotation_key", nargs="*", default=["non-targeting"])
    p.add_argument("--FDR_method", choices=["BH", "StoreyQ"], default="BH")
    p.add_argument("--save_dir", default=None)
    p.add_argument("--n_jobs", type=int, default=-1)
    p.add_argument("--calibration-only", action="store_true")
    p.add_argument("--compute-skew", action="store_true")
    p.add_argument("--ntc-ensembles", type=int, default=20)
    p.add_argument("--ntc-max-groups", type=int, default=9)
    p.add_argument("--ntc-frequency-bins", type=int, default=20)
    p.add_argument("--ntc-group-seed", type=int, default=7)
    p.add_argument("--ntc-crt-seed", type=int, default=23)
    p.add_argument("--min-cells-per-condition", type=int, default=500)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.components = args.components or [30]
    args.sel_thresh = args.sel_thresh or [0.4]
    if args.save_dir is None:
        raise ValueError("--save_dir is required for sensitivity runs")

    output_folder = Path(args.save_dir)
    output_folder.mkdir(parents=True, exist_ok=True)
    job_id = os.environ.get("SLURM_JOB_ID", "local")
    with open(output_folder / f"config_{job_id}.yml", "w") as handle:
        yaml.safe_dump(vars(args), handle, sort_keys=False)

    guide_data = mu.read(args.mdata_guide_path)

    for sel_thresh in args.sel_thresh:
        for k in args.components:
            print(f"Processing K={k}, sel_thresh={sel_thresh}")
            mdata_path = Path(args.out_dir) / args.run_name / "adata" / (
                f"cNMF_{k}_{str(sel_thresh).replace('.', '_')}.h5mu"
            )
            mdata = mu.read(mdata_path)
            adata = reformat_data_for_crt(
                mdata,
                guide_data,
                covariates=args.covariates,
                log_covariates=args.log_covariates,
                categorical_key=args.categorical_key,
            )

            U = np.asarray(adata.obsm["cnmf_usage"], dtype=float).copy()
            U = np.maximum(U, 1e-8)
            U /= U.sum(axis=1, keepdims=True)
            adata.obsm["cnmf_usage"] = U

            conditions = pd.unique(adata.obs[args.categorical_key])
            for condition in conditions:
                run_condition(
                    adata,
                    k=k,
                    condition=condition,
                    output_folder=output_folder,
                    args=args,
                )


if __name__ == "__main__":
    main()
