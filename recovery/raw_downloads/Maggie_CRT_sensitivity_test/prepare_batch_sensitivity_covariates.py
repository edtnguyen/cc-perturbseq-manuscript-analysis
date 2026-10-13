#!/usr/bin/env python3
"""Infer 2- and 4-component technical clusters and technical PCs.

The inferred labels are intended for sensitivity analyses only. They are built
from technical QC variables and residualized within known groups (sample,
sublib by default) to reduce the chance of rediscovering known biology/batches.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import RobustScaler


DEFAULT_FEATURES = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "guides_per_cell",
    "guide_umis_per_cell",
]
DEFAULT_LOG1P = [
    "total_counts",
    "n_genes_by_counts",
    "guides_per_cell",
    "guide_umis_per_cell",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input-h5ad", required=True)
    p.add_argument("--output-h5ad", required=True)
    p.add_argument("--diagnostics-dir", default=None)
    p.add_argument("--technical-features", nargs="+", default=DEFAULT_FEATURES)
    p.add_argument("--log1p-features", nargs="*", default=DEFAULT_LOG1P)
    p.add_argument("--known-keys", nargs="*", default=["sample", "sublib"])
    p.add_argument("--n-pcs", type=int, default=4)
    p.add_argument("--fit-max-cells", type=int, default=100_000)
    p.add_argument("--stability-runs", type=int, default=8)
    p.add_argument("--seed", type=int, default=20260619)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def _sanitize_obs_for_h5ad(obs: pd.DataFrame) -> pd.DataFrame:
    """Avoid pandas nullable-string categories that older anndata won't write."""
    out = obs.copy()
    for col in out.columns:
        s = out[col]
        if isinstance(s.dtype, pd.StringDtype):
            out[col] = s.astype(object)
        elif isinstance(s.dtype, pd.CategoricalDtype):
            values = s.astype(object).to_numpy(dtype=object)
            # Preserve missing values while forcing ordinary Python-string categories.
            values = np.array(
                [None if pd.isna(v) else str(v) for v in values], dtype=object
            )
            out[col] = pd.Categorical(values, ordered=s.cat.ordered)
    return out


def _combined_category(obs: pd.DataFrame, keys: Iterable[str]) -> pd.Categorical:
    arrays = []
    for key in keys:
        values = obs[key].astype(object).to_numpy(dtype=object)
        arrays.append(np.array(["<NA>" if pd.isna(v) else str(v) for v in values]))
    joined = np.array(["__".join(parts) for parts in zip(*arrays)], dtype=object)
    return pd.Categorical(joined)


def _orient_labels(
    labels: np.ndarray,
    probs: np.ndarray,
    orientation_values: np.ndarray,
    n_components: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Order labels by median original total_counts for deterministic names."""
    medians = []
    for c in range(n_components):
        medians.append(np.nanmedian(orientation_values[labels == c]))
    order = np.argsort(medians)
    old_to_new = {int(old): int(new) for new, old in enumerate(order)}
    relabeled = np.array([old_to_new[int(x)] for x in labels], dtype=np.int16)
    prob_ordered = probs[:, order]
    return relabeled, prob_ordered


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_h5ad)
    output_path = Path(args.output_h5ad)
    diagnostics_dir = (
        Path(args.diagnostics_dir)
        if args.diagnostics_dir
        else output_path.with_suffix("").with_name(output_path.stem + "_diagnostics")
    )

    if not input_path.exists():
        raise FileNotFoundError(input_path)
    if output_path.exists() and not args.overwrite:
        raise FileExistsError(f"{output_path} exists; pass --overwrite to replace it")

    print(f"Reading {input_path}")
    adata = ad.read_h5ad(input_path)
    features = list(args.technical_features)
    missing = [c for c in features if c not in adata.obs]
    if missing:
        raise KeyError(f"Missing technical features in adata.obs: {missing}")

    X_raw = (
        adata.obs[features]
        .apply(pd.to_numeric, errors="coerce")
        .astype(float)
        .replace([np.inf, -np.inf], np.nan)
    )
    X = X_raw.copy()
    for col in args.log1p_features:
        if col in X:
            X.loc[X[col] < 0, col] = np.nan
            X[col] = np.log1p(X[col])

    X = pd.DataFrame(
        SimpleImputer(strategy="median").fit_transform(X),
        index=adata.obs_names,
        columns=features,
    )

    known_keys = [k for k in args.known_keys if k in adata.obs]
    if known_keys:
        known_df = adata.obs[known_keys].copy()
        known_index = pd.MultiIndex.from_arrays(
            [
                np.array(["<NA>" if pd.isna(v) else str(v) for v in known_df[k]], dtype=object)
                for k in known_keys
            ],
            names=known_keys,
        )
        codes, _ = pd.factorize(known_index, sort=True)
        X = X - X.groupby(codes, sort=False).transform("median")
        print(f"Residualized within known groups: {known_keys}")

    X = X.clip(lower=X.quantile(0.01), upper=X.quantile(0.99), axis=1)
    Z = RobustScaler(quantile_range=(25, 75)).fit_transform(X)
    n_pcs = max(1, min(args.n_pcs, Z.shape[1]))
    pca = PCA(n_components=n_pcs, random_state=args.seed)
    Z_pca = pca.fit_transform(Z)

    for j in range(n_pcs):
        adata.obs[f"technical_pc{j + 1}"] = Z_pca[:, j].astype(np.float32)

    rng = np.random.default_rng(args.seed)
    fit_n = min(args.fit_max_cells, adata.n_obs)
    fit_idx = rng.choice(adata.n_obs, size=fit_n, replace=False)

    models: dict[int, GaussianMixture] = {}
    comparison_rows = []
    for n_components in range(1, 5):
        model = GaussianMixture(
            n_components=n_components,
            covariance_type="full",
            n_init=20,
            reg_covar=1e-6,
            random_state=args.seed,
        ).fit(Z_pca[fit_idx])
        models[n_components] = model
        comparison_rows.append(
            {
                "n_components": n_components,
                "BIC": model.bic(Z_pca[fit_idx]),
                "AIC": model.aic(Z_pca[fit_idx]),
            }
        )

    orientation_values = pd.to_numeric(
        adata.obs["total_counts"], errors="coerce"
    ).to_numpy(float)
    metrics: dict[str, object] = {
        "n_cells": int(adata.n_obs),
        "fit_cells": int(fit_n),
        "features": features,
        "known_keys_residualized": known_keys,
        "pca_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
    }

    eval_n = min(10_000, adata.n_obs)
    eval_idx = rng.choice(adata.n_obs, size=eval_n, replace=False)

    for n_components in (2, 4):
        model = models[n_components]
        raw_labels = model.predict(Z_pca)
        raw_probs = model.predict_proba(Z_pca)
        labels, probs = _orient_labels(
            raw_labels, raw_probs, orientation_values, n_components
        )
        label_name = f"inferred_batch_{n_components}"
        labels_str = np.array([f"batch_{x}" for x in labels], dtype=object)
        adata.obs[label_name] = pd.Categorical(labels_str)
        adata.obs[f"{label_name}_numeric"] = labels
        adata.obs[f"{label_name}_certainty"] = probs.max(axis=1).astype(np.float32)

        if "sample" in adata.obs:
            adata.obs[f"sample_x_{label_name}"] = _combined_category(
                adata.obs.assign(**{label_name: adata.obs[label_name]}),
                ["sample", label_name],
            )

        silhouette = silhouette_score(Z_pca[eval_idx], labels[eval_idx])
        stability_ari = []
        for run in range(args.stability_runs):
            rr = np.random.default_rng(args.seed + 1000 * n_components + run)
            sub_n = min(40_000, fit_n)
            sub_idx = rr.choice(fit_idx, size=sub_n, replace=False)
            m = GaussianMixture(
                n_components=n_components,
                covariance_type="full",
                n_init=10,
                reg_covar=1e-6,
                random_state=args.seed + 1000 * n_components + run,
            ).fit(Z_pca[sub_idx])
            pred = m.predict(Z_pca[eval_idx])
            stability_ari.append(adjusted_rand_score(labels[eval_idx], pred))

        metrics[f"gmm_{n_components}"] = {
            "silhouette": float(silhouette),
            "cluster_fractions": (
                pd.Series(labels).value_counts(normalize=True).sort_index().to_dict()
            ),
            "median_certainty": float(np.median(probs.max(axis=1))),
            "fraction_certainty_ge_0_8": float(np.mean(probs.max(axis=1) >= 0.8)),
            "stability_ari": [float(x) for x in stability_ari],
            "median_stability_ari": float(np.median(stability_ari)),
        }

    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    comparison = pd.DataFrame(comparison_rows).set_index("n_components")
    comparison.to_csv(diagnostics_dir / "gmm_model_comparison.csv")
    pd.DataFrame(
        pca.components_.T,
        index=features,
        columns=[f"technical_pc{i + 1}" for i in range(n_pcs)],
    ).to_csv(diagnostics_dir / "technical_pca_loadings.csv")
    with open(diagnostics_dir / "batch_inference_metrics.json", "w") as handle:
        json.dump(metrics, handle, indent=2)

    plot_n = min(30_000, adata.n_obs)
    plot_idx = rng.choice(adata.n_obs, size=plot_n, replace=False)
    for n_components in (2, 4):
        labels = adata.obs[f"inferred_batch_{n_components}_numeric"].to_numpy()
        plt.figure(figsize=(7, 6))
        sc = plt.scatter(
            Z_pca[plot_idx, 0],
            Z_pca[plot_idx, 1] if n_pcs > 1 else np.zeros(plot_n),
            c=labels[plot_idx],
            s=4,
            alpha=0.35,
        )
        plt.xlabel("technical PC1")
        plt.ylabel("technical PC2" if n_pcs > 1 else "0")
        plt.title(f"Inferred {n_components}-component technical clustering")
        plt.colorbar(sc, label="cluster")
        plt.tight_layout()
        plt.savefig(diagnostics_dir / f"pca_inferred_batch_{n_components}.png", dpi=150)
        plt.close()

    plt.figure(figsize=(6, 4))
    plt.plot(comparison.index, comparison["BIC"], marker="o")
    plt.xticks(comparison.index)
    plt.xlabel("GMM components")
    plt.ylabel("BIC (lower is better)")
    plt.tight_layout()
    plt.savefig(diagnostics_dir / "gmm_bic.png", dpi=150)
    plt.close()

    adata.obs = _sanitize_obs_for_h5ad(adata.obs)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = output_path.with_name(output_path.stem + ".tmp.h5ad")
    if tmp.exists():
        tmp.unlink()
    if output_path.exists() and args.overwrite:
        output_path.unlink()

    # Permit writing any remaining nullable strings from the source object.
    try:
        with ad.settings.override(allow_write_nullable_strings=True):
            adata.write_h5ad(tmp, compression="lzf")
    except AttributeError:
        adata.write_h5ad(tmp, compression="lzf")
    tmp.replace(output_path)
    print(f"Wrote {output_path}")
    print(f"Diagnostics: {diagnostics_dir}")
    print("Added obs columns:")
    print([c for c in adata.obs if c.startswith(("inferred_batch_", "sample_x_inferred_batch_", "technical_pc"))])


if __name__ == "__main__":
    main()
