#!/usr/bin/env python3
"""Run NTC calibration with true fixed-count stratified-permutation CRT.

This script keeps the known sample as the outer analysis condition. Within each
sample it fits the propensity model and performs fixed-count permutations inside
propensity-score bins, optionally crossed with an inferred 2- or 4-way technical
stratum. NTC pseudo-target size is configurable (default 15).
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re
from typing import Iterable

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
import scipy.sparse as sp

from src.sceptre import (
    build_ntc_group_inputs,
    compute_guide_set_null_pvals,
    crt_pvals_for_ntc_groups_ensemble,
    make_ntc_groups_ensemble,
    prepare_crt_inputs,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", required=True)
    p.add_argument("--run-name", default="Inference_guideumis")
    p.add_argument("--guide-h5ad", required=True)
    p.add_argument("--save-dir", required=True)
    p.add_argument("--component", type=int, default=30)
    p.add_argument("--sel-thresh", type=float, default=0.4)
    p.add_argument("--outer-key", default="sample")
    p.add_argument(
        "--variant",
        choices=["propensity_only", "inferred2", "inferred4"],
        required=True,
    )
    p.add_argument("--covariates", nargs="*", default=[])
    p.add_argument("--log-covariates", nargs="*", default=[])
    p.add_argument("--ntc-labels", nargs="+", default=["non-targeting"])
    p.add_argument("--guides-per-unit", type=int, default=15)
    p.add_argument("--B", type=int, default=1024)
    p.add_argument("--propensity-bins", type=int, default=20)
    p.add_argument("--min-stratum-size", type=int, default=2)
    p.add_argument("--ntc-frequency-bins", type=int, default=20)
    p.add_argument("--ntc-ensembles", type=int, default=20)
    p.add_argument("--ntc-max-groups", type=int, default=9)
    p.add_argument("--ntc-group-seed", type=int, default=7)
    p.add_argument("--ntc-crt-seed", type=int, default=23)
    p.add_argument("--min-cells", type=int, default=500)
    p.add_argument(
        "--usage-eps",
        type=float,
        default=1e-8,
        help=(
            "Floor cNMF usage entries to this value, then renormalize each "
            "cell to sum to one before CLR transformation (default: 1e-8)."
        ),
    )
    return p.parse_args()


def safe_name(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "condition"


def to_dense(x) -> np.ndarray:
    if sp.issparse(x):
        return x.toarray()
    return np.asarray(x)


def floor_and_renormalize_usage(
    usage,
    *,
    eps: float,
    label: str,
) -> np.ndarray:
    """Return a finite, strictly positive, row-normalized dense usage matrix."""
    U = to_dense(usage).astype(np.float64, copy=True)
    if U.ndim != 2:
        raise ValueError(f"{label}: cNMF usage must be 2D; got shape {U.shape}")
    if eps <= 0:
        raise ValueError(f"{label}: usage eps must be positive; got {eps}")
    if not np.isfinite(U).all():
        bad = int((~np.isfinite(U)).sum())
        raise ValueError(f"{label}: cNMF usage contains {bad} non-finite values")

    n_nonpositive = int((U <= 0).sum())
    n_floored = int((U < eps).sum())
    U = np.maximum(U, eps)
    row_sums = U.sum(axis=1, keepdims=True)
    if not np.isfinite(row_sums).all() or np.any(row_sums <= 0):
        raise ValueError(f"{label}: invalid row sums after flooring")
    U /= row_sums

    if not np.isfinite(U).all() or np.any(U <= 0):
        raise RuntimeError(f"{label}: flooring failed to produce strictly positive usage")

    print(
        f"{label}: eps={eps:g}; nonpositive_before={n_nonpositive:,}; "
        f"floored={n_floored:,}/{U.size:,}; min_after={U.min():.3e}; "
        f"row_sum_range=({U.sum(axis=1).min():.8f}, {U.sum(axis=1).max():.8f})",
        flush=True,
    )
    return U


def resolve_obs_column(
    key: str,
    outcome_obs: pd.DataFrame,
    guide_obs: pd.DataFrame,
    cells: pd.Index,
) -> pd.Series:
    if key in outcome_obs.columns:
        out = outcome_obs.loc[cells, key]
    elif key in guide_obs.columns:
        missing = cells.difference(guide_obs.index)
        if len(missing):
            raise KeyError(f"{len(missing)} outcome cells are missing from guide obs")
        out = guide_obs.loc[cells, key]
    else:
        raise KeyError(f"{key!r} is absent from both cNMF obs and guide obs")
    return pd.Series(out.to_numpy(), index=cells, name=key)


def make_analysis_adata(
    cnmf_adata,
    guide_adata,
    *,
    outer_key: str,
    variant: str,
    covariates: Iterable[str],
    log_covariates: Iterable[str],
    usage_eps: float,
):
    a = cnmf_adata.copy()
    cells = pd.Index(a.obs_names)
    g = guide_adata[cells].copy()

    # Floor and re-close before CLR, exactly as in the production CRT wrapper.
    usage = floor_and_renormalize_usage(
        a.X, eps=usage_eps, label="full analysis matrix"
    )
    a.obsm["cnmf_usage"] = usage
    a.obsm["guide_assignment"] = g.obsm["guide_assignment"].copy()
    a.obs["guides_per_cell"] = np.asarray(
        a.obsm["guide_assignment"].sum(axis=1)
    ).ravel()

    guide_names = list(guide_adata.uns["guide_names"])
    guide_targets = list(guide_adata.uns["guide_targets"])
    a.uns["guide_names"] = guide_names
    a.uns["guide2gene"] = dict(zip(guide_names, guide_targets))
    a.uns["program_names"] = [str(x) for x in a.var_names]

    outer = resolve_obs_column(outer_key, a.obs, g.obs, cells)
    a.obs[outer_key] = outer.to_numpy()

    # Make a same-length raw covariate DataFrame. The special column `batch`
    # is used both as an OLS/propensity covariate and by the stratified sampler.
    covar = pd.DataFrame(index=cells)
    for key in covariates:
        covar[key] = resolve_obs_column(key, a.obs, g.obs, cells).to_numpy()

    for key in log_covariates:
        values = pd.to_numeric(
            resolve_obs_column(key, a.obs, g.obs, cells), errors="coerce"
        )
        if values.isna().any() or (values < 0).any():
            raise ValueError(f"Invalid values in log covariate {key!r}")
        covar[f"log1p_{key}"] = np.log1p(values.to_numpy(float))

    inferred_key = None
    if variant == "inferred2":
        inferred_key = "inferred_batch_2"
    elif variant == "inferred4":
        inferred_key = "inferred_batch_4"

    if inferred_key is not None:
        covar["batch"] = resolve_obs_column(
            inferred_key, a.obs, g.obs, cells
        ).astype(object).to_numpy()

    a.obsm["covar"] = covar
    return a


def qq_coords(pvals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    p = np.asarray(pvals, dtype=float)
    p = np.clip(p[np.isfinite(p)], 1e-300, 1.0)
    p.sort()
    if p.size == 0:
        return np.empty(0), np.empty(0)
    expected = (np.arange(1, p.size + 1) - 0.5) / p.size
    return -np.log10(expected), -np.log10(p)


def save_results(
    *,
    output_dir: Path,
    component: int,
    condition: object,
    ntc_pvals: dict[int, pd.DataFrame],
    null_pvals: np.ndarray,
    group_defs: list[dict[str, list[str]]],
    variant: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    cond = safe_name(condition)

    frames = []
    for ensemble, frame in ntc_pvals.items():
        out = frame.copy()
        out.index.name = "group_id"
        out = out.reset_index()
        out.insert(0, "ensemble", ensemble)
        frames.append(out)
    wide = pd.concat(frames, ignore_index=True)
    wide.to_csv(
        output_dir / f"{component}_CRT_{cond}_ntc_null_pvals.txt",
        sep="\t",
        index=False,
    )

    pd.DataFrame({"null_pval": null_pvals}).to_csv(
        output_dir / f"{component}_CRT_{cond}_permutation_null_pvals.txt",
        sep="\t",
        index=False,
    )

    group_rows = []
    for ensemble, groups in enumerate(group_defs):
        for group_id, guides in groups.items():
            for guide in guides:
                group_rows.append(
                    {"ensemble": ensemble, "group_id": group_id, "guide": guide}
                )
    pd.DataFrame(group_rows).to_csv(
        output_dir / f"{component}_CRT_{cond}_ntc_group_guides.txt",
        sep="\t",
        index=False,
    )

    ntc_flat = wide.drop(columns=["ensemble", "group_id"]).to_numpy(float).ravel()
    x_ntc, y_ntc = qq_coords(ntc_flat)
    x_null, y_null = qq_coords(null_pvals)
    limit = max(
        1.0,
        float(np.nanmax(x_ntc)) if x_ntc.size else 0.0,
        float(np.nanmax(y_ntc)) if y_ntc.size else 0.0,
        float(np.nanmax(x_null)) if x_null.size else 0.0,
        float(np.nanmax(y_null)) if y_null.size else 0.0,
    )
    plt.figure(figsize=(7, 6))
    plt.plot([0, limit], [0, limit], color="black", linewidth=1, label="y=x")
    if x_ntc.size:
        plt.plot(x_ntc, y_ntc, linewidth=1.4, label="NTC calibration pairs")
    if x_null.size:
        plt.plot(x_null, y_null, "--", linewidth=1.2, label="CRT-null LOO")
    plt.xlim(0, limit)
    plt.ylim(0, limit)
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.title(f"True stratified CRT: {variant}; {condition}")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output_dir / f"CRT_{cond}_true_stratified.png", dpi=160)
    plt.close()


def main() -> None:
    args = parse_args()
    out_root = Path(args.out_dir)
    guide_path = Path(args.guide_h5ad)
    save_dir = Path(args.save_dir)

    h5mu = (
        out_root
        / args.run_name
        / "adata"
        / f"cNMF_{args.component}_{str(args.sel_thresh).replace('.', '_')}.h5mu"
    )
    if not h5mu.exists():
        raise FileNotFoundError(h5mu)
    if not guide_path.exists():
        raise FileNotFoundError(guide_path)

    print(f"Reading cNMF data: {h5mu}")
    mdata = mu.read(h5mu)
    print(f"Reading guide data: {guide_path}")
    guide_data = ad.read_h5ad(guide_path)

    analysis = make_analysis_adata(
        mdata["cNMF"],
        guide_data,
        outer_key=args.outer_key,
        variant=args.variant,
        covariates=args.covariates,
        log_covariates=args.log_covariates,
        usage_eps=args.usage_eps,
    )

    for condition in pd.unique(analysis.obs[args.outer_key]):
        cell_mask = analysis.obs[args.outer_key].astype(str) == str(condition)
        a = analysis[cell_mask].copy()
        if a.n_obs < args.min_cells:
            print(f"Skipping {condition!r}: only {a.n_obs} cells")
            continue

        use_batch = args.variant != "propensity_only"
        print(
            f"Condition={condition!r}, cells={a.n_obs:,}, variant={args.variant}, "
            f"stratify_by_batch={use_batch}"
        )

        # Re-floor after subsetting as a defensive check. This is intentionally
        # redundant and guarantees every array task/variant enters CLR with no zeros.
        a.obsm["cnmf_usage"] = floor_and_renormalize_usage(
            a.obsm["cnmf_usage"],
            eps=args.usage_eps,
            label=f"condition {condition!r} / variant {args.variant}",
        )
        assert float(np.min(a.obsm["cnmf_usage"])) > 0.0

        inputs = prepare_crt_inputs(
            a,
            usage_key="cnmf_usage",
            covar_key="covar",
            guide_assignment_key="guide_assignment",
            guide2gene_key="guide2gene",
            batch_key="batch",
        )

        resampling_kwargs = {
            "n_bins": args.propensity_bins,
            "stratify_by_batch": use_batch,
            "batch_key": "batch",
            "min_stratum_size": args.min_stratum_size,
        }

        ntc_guides, guide_freq, guide_to_bin, real_sigs = build_ntc_group_inputs(
            inputs=inputs,
            ntc_label=args.ntc_labels,
            group_size=args.guides_per_unit,
            n_bins=args.ntc_frequency_bins,
        )
        ntc_groups = make_ntc_groups_ensemble(
            ntc_guides=ntc_guides,
            ntc_freq=guide_freq,
            real_gene_bin_sigs=real_sigs,
            guide_to_bin=guide_to_bin,
            n_ensemble=args.ntc_ensembles,
            seed0=args.ntc_group_seed,
            group_size=args.guides_per_unit,
            max_groups=args.ntc_max_groups,
        )

        ntc_pvals = crt_pvals_for_ntc_groups_ensemble(
            inputs=inputs,
            ntc_groups_ens=ntc_groups,
            B=args.B,
            seed0=args.ntc_crt_seed,
            resampling_method="stratified_perm",
            resampling_kwargs=resampling_kwargs,
            expected_group_size=args.guides_per_unit,
            test_stat="ols",
        )

        guide_to_col = {guide: i for i, guide in enumerate(inputs.guide_names)}
        null_chunks = []
        for ensemble_idx, groups in enumerate(ntc_groups):
            for group_idx, guides in groups.items():
                columns = [guide_to_col[g] for g in guides]
                null_chunks.append(
                    compute_guide_set_null_pvals(
                        guide_idx=columns,
                        inputs=inputs,
                        B=args.B,
                        base_seed=(
                            args.ntc_crt_seed
                            + 100_000 * ensemble_idx
                            + (hash(group_idx) & 0xFFFF)
                        ),
                        resampling_method="stratified_perm",
                        resampling_kwargs=resampling_kwargs,
                        test_stat="ols",
                    ).ravel()
                )
        null_pvals = np.concatenate(null_chunks) if null_chunks else np.empty(0)

        save_results(
            output_dir=save_dir,
            component=args.component,
            condition=condition,
            ntc_pvals=ntc_pvals,
            null_pvals=null_pvals,
            group_defs=ntc_groups,
            variant=args.variant,
        )


if __name__ == "__main__":
    main()
