# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/Maggie.VIC.CRT.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/Maggie.VIC.CRT.ipynb`
- SHA-256: `1ff80f6a602e3217a0dc5e6e86fae291d4d1a424795f8efd66f1e8b50538ae29`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Code cell 0

```python
from pathlib import Path
import json
import warnings

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from IPython.display import display
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    normalized_mutual_info_score,
    silhouette_score,
)
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import RobustScaler


# ============================================================
# 0. Configuration
# ============================================================

INPUT_H5AD = Path(
    "/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/"
    "20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/"
    "20260408-postprocessing/9.cnmf/"
    "20260505_k30_300kvicchd_d9_torch/"
    "CRT_guideorder_allcovars_nullskewcalc/"
    "rna_with_guides.h5ad"
)

OUTPUT_H5AD = INPUT_H5AD.with_name(
    "rna_with_guides_inferred_processing_batch.h5ad"
)

DIAGNOSTICS_DIR = INPUT_H5AD.parent / "inferred_batch_diagnostics"

# Use technical features only. Do not include cNMF usage,
# tested expression outcomes, cell-cycle programs, or Leiden clusters.
TECHNICAL_FEATURES = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "guides_per_cell",
    "guide_umis_per_cell",
]

LOG1P_FEATURES = {
    "total_counts",
    "n_genes_by_counts",
    "guides_per_cell",
    "guide_umis_per_cell",
}

# Remove already-known sample/sublibrary effects before searching for
# an additional hidden processing batch.
#
# Set to [] if you explicitly do not want this residualization.
KNOWN_KEYS_TO_RESIDUALIZE = [
    "sample",
    "sublib",
]

# These are used only to diagnose whether the inferred clusters
# reproduce known sample or biological groups.
DIAGNOSTIC_CATEGORICAL_KEYS = [
    "sample",
    "sublib",
    "leiden_res0_25",
    "cell_type",
    "celltype",
]

ORIENTATION_FEATURE = "total_counts"

SEED = 20260619
N_PCS = 4

# Cap expensive operations for large datasets.
FIT_MAX_CELLS = 100_000
SILHOUETTE_MAX_CELLS = 5_000
PLOT_MAX_CELLS = 30_000

# Stability analysis.
N_STABILITY_RUNS = 8
STABILITY_FIT_MAX_CELLS = 40_000
STABILITY_EVAL_MAX_CELLS = 10_000

# Winsorize technical outliers before PCA/GMM.
WINSOR_LOW = 0.01
WINSOR_HIGH = 0.99

SAVE_UPDATED_H5AD = True
OVERWRITE_OUTPUT = False


# ============================================================
# 1. Load AnnData
# ============================================================

if not INPUT_H5AD.exists():
    raise FileNotFoundError(
        f"Input file does not exist:\n{INPUT_H5AD}"
    )

print(f"Reading:\n{INPUT_H5AD}")
adata = ad.read_h5ad(INPUT_H5AD)

print(f"\nLoaded AnnData:")
print(f"  cells:    {adata.n_obs:,}")
print(f"  features: {adata.n_vars:,}")
print(f"  obs columns: {len(adata.obs.columns):,}")


# ============================================================
# 2. Validate and prepare technical features
# ============================================================

missing_features = [
    c for c in TECHNICAL_FEATURES
    if c not in adata.obs.columns
]

if missing_features:
    raise KeyError(
        "Missing required technical columns in adata.obs:\n"
        f"{missing_features}\n\n"
        f"Available columns include:\n{list(adata.obs.columns)}"
    )

X_numeric = (
    adata.obs[TECHNICAL_FEATURES]
    .apply(pd.to_numeric, errors="coerce")
    .astype(float)
    .replace([np.inf, -np.inf], np.nan)
)

# Check for unusable columns.
all_missing = [
    c for c in TECHNICAL_FEATURES
    if X_numeric[c].notna().sum() == 0
]

if all_missing:
    raise ValueError(
        f"These columns contain no usable numeric values: {all_missing}"
    )

# Log-transform count-like features.
X_transformed = X_numeric.copy()

for col in LOG1P_FEATURES:
    X_transformed.loc[X_transformed[col] < 0, col] = np.nan
    X_transformed[col] = np.log1p(X_transformed[col])

# Median imputation.
imputer = SimpleImputer(strategy="median")

X_imputed = pd.DataFrame(
    imputer.fit_transform(X_transformed),
    index=adata.obs_names,
    columns=TECHNICAL_FEATURES,
)

print("\nTechnical feature missingness before imputation:")
display(X_numeric.isna().mean().sort_values(ascending=False).to_frame("missing_fraction"))


# ============================================================
# 3. Remove known sample/sublibrary shifts
# ============================================================

known_keys = [
    k for k in KNOWN_KEYS_TO_RESIDUALIZE
    if k in adata.obs.columns
]

missing_known_keys = [
    k for k in KNOWN_KEYS_TO_RESIDUALIZE
    if k not in adata.obs.columns
]

if missing_known_keys:
    print(
        "Known keys not found and therefore ignored:",
        missing_known_keys,
    )

if known_keys:
    known_frame = (
        adata.obs[known_keys]
        .astype("string")
        .fillna("<NA>")
    )

    known_index = pd.MultiIndex.from_frame(known_frame)
    known_group_codes, known_group_levels = pd.factorize(
        known_index,
        sort=True,
    )

    known_group_medians = X_imputed.groupby(
        known_group_codes,
        sort=False,
    ).transform("median")

    X_residualized = X_imputed - known_group_medians

    print(
        f"\nResidualized technical features within "
        f"{len(known_group_levels):,} known "
        f"{' × '.join(known_keys)} groups."
    )
else:
    X_residualized = X_imputed.copy()
    print("\nNo known-group residualization was applied.")

# Winsorize residual technical features.
lower = X_residualized.quantile(WINSOR_LOW)
upper = X_residualized.quantile(WINSOR_HIGH)

X_clipped = X_residualized.clip(
    lower=lower,
    upper=upper,
    axis=1,
)

# Robust scaling.
scaler = RobustScaler(quantile_range=(25, 75))
Z = scaler.fit_transform(X_clipped)

if not np.isfinite(Z).all():
    raise ValueError("Non-finite values remain after preprocessing.")


# ============================================================
# 4. PCA on technical features
# ============================================================

n_pcs = max(1, min(N_PCS, Z.shape[1]))

pca = PCA(
    n_components=n_pcs,
    random_state=SEED,
)

Z_pca = pca.fit_transform(Z)

explained = pd.Series(
    pca.explained_variance_ratio_,
    index=[f"PC{i + 1}" for i in range(n_pcs)],
    name="explained_variance_fraction",
)

print("\nPCA variance explained:")
display(explained.to_frame())

loadings = pd.DataFrame(
    pca.components_.T,
    index=TECHNICAL_FEATURES,
    columns=[f"PC{i + 1}" for i in range(n_pcs)],
)

print("\nPCA loadings:")
display(loadings)


# ============================================================
# 5. Fit Gaussian mixture models with 1–4 components
# ============================================================

rng = np.random.default_rng(SEED)

fit_n = min(FIT_MAX_CELLS, adata.n_obs)

fit_idx = rng.choice(
    adata.n_obs,
    size=fit_n,
    replace=False,
)

gmms = {}
bic = {}
aic = {}

for n_components in range(1, 5):
    print(f"Fitting GMM with {n_components} component(s)...")

    model = GaussianMixture(
        n_components=n_components,
        covariance_type="full",
        n_init=20,
        reg_covar=1e-6,
        random_state=SEED,
    )

    model.fit(Z_pca[fit_idx])

    gmms[n_components] = model
    bic[n_components] = model.bic(Z_pca[fit_idx])
    aic[n_components] = model.aic(Z_pca[fit_idx])

model_comparison = pd.DataFrame(
    {
        "n_components": list(bic),
        "BIC": [bic[k] for k in bic],
        "AIC": [aic[k] for k in aic],
    }
).set_index("n_components")

model_comparison["delta_BIC_from_best"] = (
    model_comparison["BIC"]
    - model_comparison["BIC"].min()
)

print("\nModel comparison:")
display(model_comparison)

gmm2 = gmms[2]

raw_labels = gmm2.predict(Z_pca)
raw_probabilities = gmm2.predict_proba(Z_pca)


# ============================================================
# 6. Give cluster labels a deterministic orientation
# ============================================================

orientation_values = X_imputed[ORIENTATION_FEATURE]

orientation_medians = (
    pd.Series(
        orientation_values.to_numpy(),
        index=adata.obs_names,
    )
    .groupby(raw_labels)
    .median()
)

# hidden_1 is defined as the cluster with higher median total_counts.
high_component = int(orientation_medians.idxmax())

if high_component == 1:
    labels = raw_labels.astype(np.int8)
    probabilities = raw_probabilities.copy()
else:
    labels = (1 - raw_labels).astype(np.int8)
    probabilities = raw_probabilities[:, ::-1]

batch_strings = np.where(
    labels == 0,
    "hidden_0",
    "hidden_1",
)

certainty = probabilities.max(axis=1)


# ============================================================
# 7. Add inferred batch labels to adata.obs
# ============================================================

adata.obs["inferred_processing_batch"] = pd.Categorical(
    batch_strings
)

adata.obs["inferred_processing_batch_numeric"] = labels

adata.obs[
    "inferred_processing_batch_prob_hidden_1"
] = probabilities[:, 1].astype(np.float32)

adata.obs[
    "inferred_processing_batch_certainty"
] = certainty.astype(np.float32)

# Diagnostic-only category marking uncertain cells.
high_conf_labels = np.where(
    certainty >= 0.80,
    batch_strings,
    "uncertain",
)

adata.obs[
    "inferred_processing_batch_high_conf"
] = pd.Categorical(high_conf_labels)

# Create the categorical key for your CRT command.
if "sample" in adata.obs.columns:
    combined_key = (
        adata.obs["sample"].astype("string").fillna("<NA>")
        + "__"
        + adata.obs["inferred_processing_batch"].astype("string")
    )
else:
    combined_key = adata.obs[
        "inferred_processing_batch"
    ].astype("string")

adata.obs[
    "sample_x_inferred_processing_batch"
] = pd.Categorical(combined_key)


# ============================================================
# 8. Separation diagnostics
# ============================================================

cluster_counts = pd.Series(
    labels
).value_counts().sort_index()

cluster_fractions = cluster_counts / adata.n_obs

# Use a subset for metrics whose cost grows quickly with n.
evaluation_n = min(
    SILHOUETTE_MAX_CELLS,
    adata.n_obs,
)

evaluation_idx = rng.choice(
    adata.n_obs,
    size=evaluation_n,
    replace=False,
)

silhouette = silhouette_score(
    Z_pca[evaluation_idx],
    labels[evaluation_idx],
)

calinski_harabasz = calinski_harabasz_score(
    Z_pca[evaluation_idx],
    labels[evaluation_idx],
)

davies_bouldin = davies_bouldin_score(
    Z_pca[evaluation_idx],
    labels[evaluation_idx],
)

# Mahalanobis separation of the two GMM components.
mean_difference = (
    gmm2.means_[1] - gmm2.means_[0]
)

pooled_covariance = (
    gmm2.covariances_[0]
    + gmm2.covariances_[1]
) / 2.0

mahalanobis_distance = float(
    np.sqrt(
        mean_difference
        @ np.linalg.pinv(pooled_covariance)
        @ mean_difference
    )
)


# ============================================================
# 9. Stability across repeated GMM fits
# ============================================================

stability_eval_n = min(
    STABILITY_EVAL_MAX_CELLS,
    adata.n_obs,
)

stability_eval_idx = rng.choice(
    adata.n_obs,
    size=stability_eval_n,
    replace=False,
)

base_stability_labels = labels[stability_eval_idx]

stability_ari = []

for run in range(N_STABILITY_RUNS):
    run_rng = np.random.default_rng(SEED + 1_000 + run)

    stability_fit_n = min(
        STABILITY_FIT_MAX_CELLS,
        fit_idx.size,
    )

    stability_fit_idx = run_rng.choice(
        fit_idx,
        size=stability_fit_n,
        replace=False,
    )

    stability_model = GaussianMixture(
        n_components=2,
        covariance_type="full",
        n_init=10,
        reg_covar=1e-6,
        random_state=SEED + 1_000 + run,
    )

    stability_model.fit(
        Z_pca[stability_fit_idx]
    )

    stability_prediction = stability_model.predict(
        Z_pca[stability_eval_idx]
    )

    # ARI is invariant to label switching.
    stability_ari.append(
        adjusted_rand_score(
            base_stability_labels,
            stability_prediction,
        )
    )

stability_ari = np.asarray(stability_ari)


# ============================================================
# 10. Feature differences between inferred clusters
# ============================================================

def standardized_mean_difference(values, group_labels):
    values = np.asarray(values, dtype=float)
    group_labels = np.asarray(group_labels)

    x0 = values[group_labels == 0]
    x1 = values[group_labels == 1]

    pooled_sd = np.sqrt(
        (
            np.nanvar(x0, ddof=1)
            + np.nanvar(x1, ddof=1)
        ) / 2.0
    )

    if not np.isfinite(pooled_sd) or pooled_sd == 0:
        return np.nan

    return (
        np.nanmean(x1) - np.nanmean(x0)
    ) / pooled_sd


smd = pd.Series(
    {
        col: standardized_mean_difference(
            X_imputed[col].to_numpy(),
            labels,
        )
        for col in TECHNICAL_FEATURES
    },
    name="standardized_mean_difference_hidden1_minus_hidden0",
).sort_values(key=np.abs, ascending=False)

raw_feature_summary = (
    X_numeric.assign(
        inferred_processing_batch=batch_strings
    )
    .groupby(
        "inferred_processing_batch",
        observed=True,
    )
    .agg(["count", "median", "mean"])
)

print("\nCluster sizes:")
display(
    pd.DataFrame(
        {
            "n_cells": cluster_counts,
            "fraction": cluster_fractions,
        }
    )
)

print("\nPosterior assignment certainty:")
display(
    pd.Series(
        certainty,
        name="certainty",
    ).describe(
        percentiles=[0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
    ).to_frame()
)

print("\nTechnical feature differences:")
display(smd.to_frame())

print("\nRaw technical-feature summaries by inferred batch:")
display(raw_feature_summary)


# ============================================================
# 11. Does the inferred label merely reproduce known groups?
# ============================================================

diagnostic_keys = [
    key
    for key in DIAGNOSTIC_CATEGORICAL_KEYS
    if key in adata.obs.columns
]

nmi_results = {}

for key in diagnostic_keys:
    key_values = (
        adata.obs[key]
        .astype("string")
        .fillna("<NA>")
    )

    key_codes, _ = pd.factorize(
        key_values,
        sort=True,
    )

    nmi = normalized_mutual_info_score(
        key_codes,
        labels,
    )

    nmi_results[key] = nmi

    print(f"\nInferred-batch proportions within {key}:")
    display(
        pd.crosstab(
            adata.obs[key],
            adata.obs["inferred_processing_batch"],
            normalize="index",
        )
    )

print("\nNormalized mutual information with known labels:")
display(
    pd.Series(
        nmi_results,
        name="NMI_with_inferred_batch",
    ).sort_values(ascending=False).to_frame()
)


# ============================================================
# 12. Summary metrics and cautious interpretation
# ============================================================

bic_1_minus_2 = bic[1] - bic[2]
bic_2_minus_3 = bic[2] - bic[3]
best_component_count = min(bic, key=bic.get)

metrics = {
    "n_cells": int(adata.n_obs),
    "fit_cells": int(fit_n),
    "bic_1_minus_2": float(bic_1_minus_2),
    "bic_2_minus_3": float(bic_2_minus_3),
    "best_component_count_by_bic": int(best_component_count),
    "silhouette_score": float(silhouette),
    "calinski_harabasz_score": float(calinski_harabasz),
    "davies_bouldin_score": float(davies_bouldin),
    "gmm_mahalanobis_distance": float(mahalanobis_distance),
    "median_assignment_certainty": float(np.median(certainty)),
    "fraction_certainty_ge_0_80": float(np.mean(certainty >= 0.80)),
    "fraction_certainty_ge_0_90": float(np.mean(certainty >= 0.90)),
    "smallest_cluster_fraction": float(cluster_fractions.min()),
    "median_stability_ARI": float(np.median(stability_ari)),
    "minimum_stability_ARI": float(np.min(stability_ari)),
}

print("\nCore diagnostics:")
display(
    pd.Series(metrics, name="value").to_frame()
)

checks = {
    "BIC favors 2 over 1": bic_1_minus_2 > 10,
    "both clusters contain at least 5% of cells": (
        cluster_fractions.min() >= 0.05
    ),
    "median posterior certainty >= 0.80": (
        np.median(certainty) >= 0.80
    ),
    "at least 70% of cells have certainty >= 0.80": (
        np.mean(certainty >= 0.80) >= 0.70
    ),
    "silhouette score >= 0.15": silhouette >= 0.15,
    "median stability ARI >= 0.90": (
        np.median(stability_ari) >= 0.90
    ),
}

print("\nHeuristic evidence checks:")
display(
    pd.Series(checks, name="passed").to_frame()
)

n_passed = sum(checks.values())

print("\nInterpretation:")
if n_passed >= 5:
    print(
        "There is strong evidence for a stable two-cluster structure "
        "in the available technical features."
    )
elif n_passed >= 3:
    print(
        "There is moderate evidence for two technical clusters, but "
        "the inferred labels should be treated cautiously."
    )
else:
    print(
        "Evidence for a meaningful two-cluster structure is weak. "
        "Do not treat these labels as recovered processing batches."
    )

if best_component_count != 2:
    print(
        f"\nWARNING: BIC prefers {best_component_count} components rather "
        "than exactly 2. A two-batch model may be an oversimplification."
    )

if nmi_results:
    maximum_nmi_key = max(
        nmi_results,
        key=nmi_results.get,
    )
    maximum_nmi = nmi_results[maximum_nmi_key]

    if maximum_nmi >= 0.80:
        print(
            f"\nWARNING: inferred batch is highly aligned with "
            f"{maximum_nmi_key!r} (NMI={maximum_nmi:.3f}). "
            "It may be rediscovering a known group or biological state."
        )

print(
    "\nImportant: clustering of technical features is not proof that "
    "the groups are the missing processing batches. Use metadata, "
    "library identifiers, sequencing lanes, or hashing/genotype data "
    "for confirmation whenever possible."
)


# ============================================================
# 13. Diagnostic plots
# ============================================================

# PCA scatter
plot_n = min(PLOT_MAX_CELLS, adata.n_obs)
plot_idx = rng.choice(
    adata.n_obs,
    size=plot_n,
    replace=False,
)

plt.figure(figsize=(7, 6))
scatter = plt.scatter(
    Z_pca[plot_idx, 0],
    Z_pca[plot_idx, 1] if n_pcs > 1 else np.zeros(plot_n),
    c=labels[plot_idx],
    s=4,
    alpha=0.35,
)
plt.xlabel(f"Technical PC1 ({pca.explained_variance_ratio_[0]:.1%})")
plt.ylabel(
    f"Technical PC2 ({pca.explained_variance_ratio_[1]:.1%})"
    if n_pcs > 1
    else "0"
)
plt.title("Inferred processing batch from technical QC features")
plt.colorbar(scatter, ticks=[0, 1], label="Inferred batch")
plt.tight_layout()
plt.show()

# BIC model comparison
plt.figure(figsize=(6, 4))
plt.plot(
    model_comparison.index,
    model_comparison["BIC"],
    marker="o",
)
plt.xlabel("Number of Gaussian mixture components")
plt.ylabel("BIC (lower is better)")
plt.title("GMM model comparison")
plt.xticks(model_comparison.index)
plt.tight_layout()
plt.show()

# Certainty histogram
plt.figure(figsize=(6, 4))
plt.hist(
    certainty,
    bins=40,
)
plt.xlabel("Maximum posterior assignment probability")
plt.ylabel("Number of cells")
plt.title("Certainty of inferred batch assignments")
plt.tight_layout()
plt.show()

# Technical-feature SMD
plt.figure(figsize=(7, 4))
smd.sort_values().plot(kind="barh")
plt.axvline(0, linewidth=1)
plt.xlabel("Standardized mean difference: hidden_1 − hidden_0")
plt.title("Technical differences between inferred groups")
plt.tight_layout()
plt.show()

# Stability ARI
plt.figure(figsize=(6, 4))
plt.plot(
    np.arange(1, len(stability_ari) + 1),
    stability_ari,
    marker="o",
)
plt.axhline(0.90, linestyle="--", linewidth=1)
plt.ylim(-0.05, 1.05)
plt.xlabel("Repeated GMM fit")
plt.ylabel("Adjusted Rand index vs primary fit")
plt.title("Stability of inferred two-cluster solution")
plt.tight_layout()
plt.show()


# ============================================================
# 14. Save diagnostics and modified AnnData
# ============================================================

DIAGNOSTICS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

model_comparison.to_csv(
    DIAGNOSTICS_DIR / "gmm_model_comparison.csv"
)

loadings.to_csv(
    DIAGNOSTICS_DIR / "technical_pca_loadings.csv"
)

smd.to_csv(
    DIAGNOSTICS_DIR / "technical_feature_smd.csv"
)

raw_feature_summary.to_csv(
    DIAGNOSTICS_DIR
    / "technical_feature_summary_by_inferred_batch.csv"
)

adata.obs[
    [
        "inferred_processing_batch",
        "inferred_processing_batch_numeric",
        "inferred_processing_batch_prob_hidden_1",
        "inferred_processing_batch_certainty",
        "inferred_processing_batch_high_conf",
        "sample_x_inferred_processing_batch",
    ]
].to_csv(
    DIAGNOSTICS_DIR / "inferred_batch_assignments.csv.gz"
)

with open(
    DIAGNOSTICS_DIR / "inferred_batch_metrics.json",
    "w",
) as handle:
    json.dump(
        {
            **metrics,
            "BIC": {
                str(k): float(v)
                for k, v in bic.items()
            },
            "stability_ARI": stability_ari.tolist(),
            "known_group_NMI": {
                str(k): float(v)
                for k, v in nmi_results.items()
            },
            "features": TECHNICAL_FEATURES,
            "log1p_features": sorted(LOG1P_FEATURES),
            "known_keys_residualized": known_keys,
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad": str(OUTPUT_H5AD),
        },
        handle,
        indent=2,
    )

if SAVE_UPDATED_H5AD:
    if OUTPUT_H5AD.exists() and not OVERWRITE_OUTPUT:
        raise FileExistsError(
            f"Output already exists:\n{OUTPUT_H5AD}\n\n"
            "Set OVERWRITE_OUTPUT = True to replace it."
        )

    print(f"\nWriting updated AnnData:\n{OUTPUT_H5AD}")
    adata.write_h5ad(
        OUTPUT_H5AD,
        compression="lzf",
    )

print(f"\nDiagnostics written to:\n{DIAGNOSTICS_DIR}")

print("\nUse this in the CRT command:")
print(f'  --mdata_guide_path "{OUTPUT_H5AD}"')
print('  --categorical_key "sample_x_inferred_processing_batch"')
```

### Text output (stream)

```text
Matplotlib is building the font cache; this may take a moment.

```

### Text output (stream)

```text
Reading:
/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/9.cnmf/20260505_k30_300kvicchd_d9_torch/CRT_guideorder_allcovars_nullskewcalc/rna_with_guides.h5ad

Loaded AnnData:
  cells:    240,504
  features: 25,788
  obs columns: 29

Technical feature missingness before imputation:

```

### Text output (text/plain)

```text
                     missing_fraction
total_counts                      0.0
n_genes_by_counts                 0.0
pct_counts_mt                     0.0
pct_counts_ribo                   0.0
guides_per_cell                   0.0
guide_umis_per_cell               0.0
```

### Text output (stream)

```text

Residualized technical features within 7 known sample × sublib groups.

PCA variance explained:

```

### Text output (text/plain)

```text
     explained_variance_fraction
PC1                     0.420501
PC2                     0.274432
PC3                     0.156187
PC4                     0.123523
```

### Text output (stream)

```text

PCA loadings:

```

### Text output (text/plain)

```text
                          PC1       PC2       PC3       PC4
total_counts         0.423062 -0.404665  0.272588  0.278249
n_genes_by_counts    0.423890 -0.417010  0.229493  0.297063
pct_counts_mt        0.116155  0.427623  0.844133 -0.299770
pct_counts_ribo      0.170587  0.648146 -0.093623  0.729664
guides_per_cell      0.559996  0.191533 -0.331974 -0.424135
guide_umis_per_cell  0.533993  0.150698 -0.203705 -0.179358
```

### Text output (stream)

```text
Fitting GMM with 1 component(s)...
Fitting GMM with 2 component(s)...
Fitting GMM with 3 component(s)...
Fitting GMM with 4 component(s)...

Model comparison:

```

### Text output (text/plain)

```text
                       BIC           AIC  delta_BIC_from_best
n_components                                                 
1             1.017416e+06  1.017283e+06         28148.628509
2             9.991559e+05  9.988801e+05          9888.071201
3             9.916176e+05  9.911990e+05          2349.719922
4             9.892679e+05  9.887066e+05             0.000000
```

### Text output (stream)

```text

Cluster sizes:

```

### Text output (text/plain)

```text
   n_cells  fraction
0    75369  0.313379
1   165135  0.686621
```

### Text output (stream)

```text

Posterior assignment certainty:

```

### Text output (text/plain)

```text
           certainty
count  240504.000000
mean        0.801360
std         0.142617
min         0.500003
1%          0.509832
10%         0.592193
25%         0.690785
50%         0.813014
75%         0.927507
90%         0.989930
99%         1.000000
max         1.000000
```

### Text output (stream)

```text

Technical feature differences:

```

### Text output (text/plain)

```text
                     standardized_mean_difference_hidden1_minus_hidden0
pct_counts_ribo                                              -2.225867 
guide_umis_per_cell                                          -0.664904 
guides_per_cell                                              -0.513133 
pct_counts_mt                                                -0.277255 
total_counts                                                  0.232673 
n_genes_by_counts                                             0.216065 
```

### Text output (stream)

```text

Raw technical-feature summaries by inferred batch:

```

### Text output (text/plain)

```text
                          total_counts                             \
                                 count       median          mean   
inferred_processing_batch                                           
hidden_0                         75369  8433.740234  10360.984396   
hidden_1                        165135  9915.685547  12376.837124   

                          n_genes_by_counts                       \
                                      count  median         mean   
inferred_processing_batch                                          
hidden_0                              75369  3978.0  4164.300137   
hidden_1                             165135  4365.0  4560.556860   

                          pct_counts_mt                     pct_counts_ribo  \
                                  count    median      mean           count   
inferred_processing_batch                                                     
hidden_0                          75369  7.492577  7.769950           75369   
hidden_1                         165135  6.923040  7.139818          165135   

                                              guides_per_cell         \
                             median      mean           count median   
inferred_processing_batch                                              
hidden_0                   1.236764  1.324757           75369    6.0   
hidden_1                   0.734714  0.728394          165135    4.0   

                                    guide_umis_per_cell                     
                               mean               count median        mean  
inferred_processing_batch                                                   
hidden_0                   7.137046               75369   90.0  156.267564  
hidden_1                   5.104545              165135   39.0   77.579980  
```

### Text output (stream)

```text

Inferred-batch proportions within sample:

```

### Text output (text/plain)

```text
inferred_processing_batch  hidden_0  hidden_1
sample                                       
300k_vic_chd               0.313379  0.686621
```

### Text output (stream)

```text

Inferred-batch proportions within sublib:

```

### Text output (text/plain)

```text
inferred_processing_batch  hidden_0  hidden_1
sublib                                       
300k_sublib1               0.322500  0.677500
300k_sublib2               0.306892  0.693108
300k_sublib3               0.311081  0.688919
300k_sublib4               0.304681  0.695319
300k_sublib5               0.309497  0.690503
300k_sublib6               0.313536  0.686464
300k_sublib7               0.362639  0.637361
```

### Text output (stream)

```text

Inferred-batch proportions within leiden_res0_25:

```

### Text output (text/plain)

```text
inferred_processing_batch  hidden_0  hidden_1
leiden_res0_25                               
0                          0.254943  0.745057
1                          0.207365  0.792635
2                          0.539158  0.460842
3                          0.332701  0.667299
4                          0.543499  0.456501
```

### Text output (stream)

```text

Normalized mutual information with known labels:

```

### Text output (text/plain)

```text
                NMI_with_inferred_batch
leiden_res0_25                 0.027508
sublib                         0.000235
sample                         0.000000
```

### Text output (stream)

```text

Core diagnostics:

```

### Text output (text/plain)

```text
                                     value
n_cells                      240504.000000
fit_cells                    100000.000000
bic_1_minus_2                 18260.557307
bic_2_minus_3                  7538.351279
best_component_count_by_bic       4.000000
silhouette_score                  0.162487
calinski_harabasz_score         845.936490
davies_bouldin_score              2.065398
gmm_mahalanobis_distance          1.446007
median_assignment_certainty       0.813014
fraction_certainty_ge_0_80        0.527833
fraction_certainty_ge_0_90        0.310286
smallest_cluster_fraction         0.313379
median_stability_ARI              0.960854
minimum_stability_ARI             0.439598
```

### Text output (stream)

```text

Heuristic evidence checks:

```

### Text output (text/plain)

```text
                                              passed
BIC favors 2 over 1                             True
both clusters contain at least 5% of cells      True
median posterior certainty >= 0.80              True
at least 70% of cells have certainty >= 0.80   False
silhouette score >= 0.15                        True
median stability ARI >= 0.90                    True
```

### Text output (stream)

```text

Interpretation:
There is strong evidence for a stable two-cluster structure in the available technical features.

WARNING: BIC prefers 4 components rather than exactly 2. A two-batch model may be an oversimplification.

Important: clustering of technical features is not proof that the groups are the missing processing batches. Use metadata, library identifiers, sequencing lanes, or hashing/genotype data for confirmation whenever possible.

```

### Text output (text/plain)

```text
<Figure size 700x600 with 2 Axes>
```

### Text output (text/plain)

```text
<Figure size 600x400 with 1 Axes>
```

### Text output (text/plain)

```text
<Figure size 600x400 with 1 Axes>
```

### Text output (text/plain)

```text
<Figure size 700x400 with 1 Axes>
```

### Text output (text/plain)

```text
<Figure size 600x400 with 1 Axes>
```

### Text output (stream)

```text

Writing updated AnnData:
/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260408-postprocessing/9.cnmf/20260505_k30_300kvicchd_d9_torch/CRT_guideorder_allcovars_nullskewcalc/rna_with_guides_inferred_processing_batch.h5ad

```

### Text output (error)

```text
[31m---------------------------------------------------------------------------[39m
[31mRuntimeError[39m                              Traceback (most recent call last)
[36mCell[39m[36m [39m[32mIn[1][39m[32m, line 906[39m
[32m    900[39m         [38;5;28;01mraise[39;00m [38;5;167;01mFileExistsError[39;00m(
[32m    901[39m             [33mf[39m[33m"[39m[33mOutput already exists:[39m[38;5;130;01m\n[39;00m[38;5;132;01m{[39;00mOUTPUT_H5AD[38;5;132;01m}[39;00m[38;5;130;01m\n[39;00m[38;5;130;01m\n[39;00m[33m"[39m
[32m    902[39m             [33m"[39m[33mSet OVERWRITE_OUTPUT = True to replace it.[39m[33m"[39m
[32m    903[39m         )
[32m    905[39m     [38;5;28mprint[39m([33mf[39m[33m"[39m[38;5;130;01m\n[39;00m[33mWriting updated AnnData:[39m[38;5;130;01m\n[39;00m[38;5;132;01m{[39;00mOUTPUT_H5AD[38;5;132;01m}[39;00m[33m"[39m)
[32m--> [39m[32m906[39m     [43madata[49m[43m.[49m[43mwrite_h5ad[49m[43m([49m
[32m    907[39m [43m        [49m[43mOUTPUT_H5AD[49m[43m,[49m
[32m    908[39m [43m        [49m[43mcompression[49m[43m=[49m[33;43m"[39;49m[33;43mlzf[39;49m[33;43m"[39;49m[43m,[49m
[32m    909[39m [43m    [49m[43m)[49m
[32m    911[39m [38;5;28mprint[39m([33mf[39m[33m"[39m[38;5;130;01m\n[39;00m[33mDiagnostics written to:[39m[38;5;130;01m\n[39;00m[38;5;132;01m{[39;00mDIAGNOSTICS_DIR[38;5;132;01m}[39;00m[33m"[39m)
[32m    913[39m [38;5;28mprint[39m([33m"[39m[38;5;130;01m\n[39;00m[33mUse this in the CRT command:[39m[33m"[39m)

    [31m[... skipping hidden 1 frame][39m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_core/anndata.py:1907[39m, in [36mAnnData.write_h5ad[39m[34m(self, filename, convert_strings_to_categoricals, compression, compression_opts, as_dense)[39m
[32m   1904[39m [38;5;28;01mif[39;00m filename [38;5;129;01mis[39;00m [38;5;28;01mNone[39;00m:
[32m   1905[39m     filename = [38;5;28mself[39m.filename
[32m-> [39m[32m1907[39m [43mwrite_h5ad[49m[43m([49m
[32m   1908[39m [43m    [49m[43mPath[49m[43m([49m[43mfilename[49m[43m)[49m[43m,[49m
[32m   1909[39m [43m    [49m[38;5;28;43mself[39;49m[43m,[49m
[32m   1910[39m [43m    [49m[43mconvert_strings_to_categoricals[49m[43m=[49m[43mconvert_strings_to_categoricals[49m[43m,[49m
[32m   1911[39m [43m    [49m[43mcompression[49m[43m=[49m[43mcompression[49m[43m,[49m
[32m   1912[39m [43m    [49m[43mcompression_opts[49m[43m=[49m[43mcompression_opts[49m[43m,[49m
[32m   1913[39m [43m    [49m[43mas_dense[49m[43m=[49m[43mas_dense[49m[43m,[49m
[32m   1914[39m [43m[49m[43m)[49m
[32m   1915[39m [38;5;66;03m# Only reset the filename if the AnnData object now points to a complete new copy[39;00m
[32m   1916[39m [38;5;28;01mif[39;00m [38;5;28mself[39m.isbacked [38;5;129;01mand[39;00m [38;5;129;01mnot[39;00m [38;5;28mself[39m.is_view:

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/utils.py:320[39m, in [36mno_write_dataset_2d.<locals>.raise_error_if_dataset_2d_present[39m[34m(store, adata, *args, **kwargs)[39m
[32m    313[39m     msg = (
[32m    314[39m         [33m"[39m[33mWriting AnnData objects with a Dataset2D not supported yet. [39m[33m"[39m
[32m    315[39m         [33m"[39m[33mPlease use `ds.to_memory` to bring the dataset into memory. [39m[33m"[39m
[32m    316[39m         [33m"[39m[33mNote that if you have generated this object by concatenating several `AnnData` objects[39m[33m"[39m
[32m    317[39m         [33m"[39m[33mthe original types may be lost.[39m[33m"[39m
[32m    318[39m     )
[32m    319[39m     [38;5;28;01mraise[39;00m [38;5;167;01mNotImplementedError[39;00m(msg)
[32m--> [39m[32m320[39m [38;5;28;01mreturn[39;00m [43mwrite[49m[43m([49m[43mstore[49m[43m,[49m[43m [49m[43madata[49m[43m,[49m[43m [49m[43m*[49m[43margs[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkwargs[49m[43m)[49m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/h5ad.py:97[39m, in [36mwrite_h5ad[39m[34m(filepath, adata, as_dense, convert_strings_to_categoricals, dataset_kwargs, **kwargs)[39m
[32m     89[39m _write_x(
[32m     90[39m     f,
[32m     91[39m     adata,  [38;5;66;03m# accessing adata.X reopens adata.file if it’s backed[39;00m
[32m   (...)[39m[32m     94[39m     dataset_kwargs=dataset_kwargs,
[32m     95[39m )
[32m     96[39m _write_raw(f, adata.raw, as_dense=as_dense, dataset_kwargs=dataset_kwargs)
[32m---> [39m[32m97[39m [43mwrite_elem[49m[43m([49m[43mf[49m[43m,[49m[43m [49m[33;43m"[39;49m[33;43mobs[39;49m[33;43m"[39;49m[43m,[49m[43m [49m[43madata[49m[43m.[49m[43mobs[49m[43m,[49m[43m [49m[43mdataset_kwargs[49m[43m=[49m[43mdataset_kwargs[49m[43m)[49m
[32m     98[39m write_elem(f, [33m"[39m[33mvar[39m[33m"[39m, adata.var, dataset_kwargs=dataset_kwargs)
[32m     99[39m write_elem(f, [33m"[39m[33mobsm[39m[33m"[39m, [38;5;28mdict[39m(adata.obsm), dataset_kwargs=dataset_kwargs)

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/specs/registry.py:518[39m, in [36mwrite_elem[39m[34m(store, k, elem, dataset_kwargs)[39m
[32m    494[39m [38;5;28;01mdef[39;00m[38;5;250m [39m[34mwrite_elem[39m(
[32m    495[39m     store: GroupStorageType,
[32m    496[39m     k: [38;5;28mstr[39m,
[32m   (...)[39m[32m    499[39m     dataset_kwargs: Mapping[[38;5;28mstr[39m, Any] = MappingProxyType({}),
[32m    500[39m ) -> [38;5;28;01mNone[39;00m:
[32m    501[39m [38;5;250m    [39m[33;03m"""[39;00m
[32m    502[39m [33;03m    Write an element to a storage group using anndata encoding.[39;00m
[32m    503[39m 
[32m   (...)[39m[32m    516[39m [33;03m        E.g. for zarr this would be `chunks`, `compressor`.[39;00m
[32m    517[39m [33;03m    """[39;00m
[32m--> [39m[32m518[39m     [43mWriter[49m[43m([49m[43m_REGISTRY[49m[43m)[49m[43m.[49m[43mwrite_elem[49m[43m([49m[43mstore[49m[43m,[49m[43m [49m[43mk[49m[43m,[49m[43m [49m[43melem[49m[43m,[49m[43m [49m[43mdataset_kwargs[49m[43m=[49m[43mdataset_kwargs[49m[43m)[49m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/utils.py:243[39m, in [36mreport_write_key_on_error.<locals>.func_wrapper[39m[34m(*args, **kwargs)[39m
[32m    241[39m     [38;5;28;01mraise[39;00m [38;5;167;01mValueError[39;00m(msg)
[32m    242[39m [38;5;28;01mtry[39;00m:
[32m--> [39m[32m243[39m     [38;5;28;01mreturn[39;00m [43mfunc[49m[43m([49m[43m*[49m[43margs[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkwargs[49m[43m)[49m
[32m    244[39m [38;5;28;01mexcept[39;00m [38;5;167;01mException[39;00m [38;5;28;01mas[39;00m e:
[32m    245[39m     path = _get_display_path(store)

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/specs/registry.py:389[39m, in [36mWriter.write_elem[39m[34m(self, store, k, elem, dataset_kwargs, modifiers)[39m
[32m    386[39m write_func = [38;5;28mself[39m.find_write_func(dest_type, elem, modifiers)
[32m    388[39m [38;5;28;01mif[39;00m [38;5;28mself[39m.callback [38;5;129;01mis[39;00m [38;5;28;01mNone[39;00m:
[32m--> [39m[32m389[39m     [38;5;28;01mreturn[39;00m [43mwrite_func[49m[43m([49m[43mstore[49m[43m,[49m[43m [49m[43mk[49m[43m,[49m[43m [49m[43melem[49m[43m,[49m[43m [49m[43mdataset_kwargs[49m[43m=[49m[43mdataset_kwargs[49m[43m)[49m
[32m    390[39m [38;5;28;01mreturn[39;00m [38;5;28mself[39m.callback(
[32m    391[39m     write_func,
[32m    392[39m     store,
[32m   (...)[39m[32m    396[39m     iospec=[38;5;28mself[39m.registry.get_spec(elem),
[32m    397[39m )

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/specs/registry.py:77[39m, in [36mwrite_spec.<locals>.decorator.<locals>.wrapper[39m[34m(g, k, *args, **kwargs)[39m
[32m     75[39m [38;5;129m@wraps[39m(func)
[32m     76[39m [38;5;28;01mdef[39;00m[38;5;250m [39m[34mwrapper[39m(g: GroupStorageType, k: [38;5;28mstr[39m, *args, **kwargs):
[32m---> [39m[32m77[39m     result = [43mfunc[49m[43m([49m[43mg[49m[43m,[49m[43m [49m[43mk[49m[43m,[49m[43m [49m[43m*[49m[43margs[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkwargs[49m[43m)[49m
[32m     78[39m     g[k].attrs.setdefault([33m"[39m[33mencoding-type[39m[33m"[39m, spec.encoding_type)
[32m     79[39m     g[k].attrs.setdefault([33m"[39m[33mencoding-version[39m[33m"[39m, spec.encoding_version)

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/specs/methods.py:992[39m, in [36mwrite_dataframe[39m[34m(f, key, df, _writer, dataset_kwargs)[39m
[32m    987[39m _writer.write_elem(
[32m    988[39m     group, index_name, df.index._values, dataset_kwargs=dataset_kwargs
[32m    989[39m )
[32m    990[39m [38;5;28;01mfor[39;00m colname, series [38;5;129;01min[39;00m df.items():
[32m    991[39m     [38;5;66;03m# TODO: this should write the "true" representation of the series (i.e. the underlying array or ndarray depending)[39;00m
[32m--> [39m[32m992[39m     [43m_writer[49m[43m.[49m[43mwrite_elem[49m[43m([49m
[32m    993[39m [43m        [49m[43mgroup[49m[43m,[49m[43m [49m[43mcolname[49m[43m,[49m[43m [49m[43mseries[49m[43m.[49m[43m_values[49m[43m,[49m[43m [49m[43mdataset_kwargs[49m[43m=[49m[43mdataset_kwargs[49m
[32m    994[39m [43m    [49m[43m)[49m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/utils.py:243[39m, in [36mreport_write_key_on_error.<locals>.func_wrapper[39m[34m(*args, **kwargs)[39m
[32m    241[39m     [38;5;28;01mraise[39;00m [38;5;167;01mValueError[39;00m(msg)
[32m    242[39m [38;5;28;01mtry[39;00m:
[32m--> [39m[32m243[39m     [38;5;28;01mreturn[39;00m [43mfunc[49m[43m([49m[43m*[49m[43margs[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkwargs[49m[43m)[49m
[32m    244[39m [38;5;28;01mexcept[39;00m [38;5;167;01mException[39;00m [38;5;28;01mas[39;00m e:
[32m    245[39m     path = _get_display_path(store)

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/specs/registry.py:389[39m, in [36mWriter.write_elem[39m[34m(self, store, k, elem, dataset_kwargs, modifiers)[39m
[32m    386[39m write_func = [38;5;28mself[39m.find_write_func(dest_type, elem, modifiers)
[32m    388[39m [38;5;28;01mif[39;00m [38;5;28mself[39m.callback [38;5;129;01mis[39;00m [38;5;28;01mNone[39;00m:
[32m--> [39m[32m389[39m     [38;5;28;01mreturn[39;00m [43mwrite_func[49m[43m([49m[43mstore[49m[43m,[49m[43m [49m[43mk[49m[43m,[49m[43m [49m[43melem[49m[43m,[49m[43m [49m[43mdataset_kwargs[49m[43m=[49m[43mdataset_kwargs[49m[43m)[49m
[32m    390[39m [38;5;28;01mreturn[39;00m [38;5;28mself[39m.callback(
[32m    391[39m     write_func,
[32m    392[39m     store,
[32m   (...)[39m[32m    396[39m     iospec=[38;5;28mself[39m.registry.get_spec(elem),
[32m    397[39m )

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/specs/registry.py:77[39m, in [36mwrite_spec.<locals>.decorator.<locals>.wrapper[39m[34m(g, k, *args, **kwargs)[39m
[32m     75[39m [38;5;129m@wraps[39m(func)
[32m     76[39m [38;5;28;01mdef[39;00m[38;5;250m [39m[34mwrapper[39m(g: GroupStorageType, k: [38;5;28mstr[39m, *args, **kwargs):
[32m---> [39m[32m77[39m     result = [43mfunc[49m[43m([49m[43mg[49m[43m,[49m[43m [49m[43mk[49m[43m,[49m[43m [49m[43m*[49m[43margs[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkwargs[49m[43m)[49m
[32m     78[39m     g[k].attrs.setdefault([33m"[39m[33mencoding-type[39m[33m"[39m, spec.encoding_type)
[32m     79[39m     g[k].attrs.setdefault([33m"[39m[33mencoding-version[39m[33m"[39m, spec.encoding_version)

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/specs/methods.py:1101[39m, in [36mwrite_categorical[39m[34m(f, k, v, _writer, dataset_kwargs)[39m
[32m   1098[39m g.attrs[[33m"[39m[33mordered[39m[33m"[39m] = [38;5;28mbool[39m(v.ordered)
[32m   1100[39m _writer.write_elem(g, [33m"[39m[33mcodes[39m[33m"[39m, v.codes, dataset_kwargs=dataset_kwargs)
[32m-> [39m[32m1101[39m [43m_writer[49m[43m.[49m[43mwrite_elem[49m[43m([49m
[32m   1102[39m [43m    [49m[43mg[49m[43m,[49m[43m [49m[33;43m"[39;49m[33;43mcategories[39;49m[33;43m"[39;49m[43m,[49m[43m [49m[43mv[49m[43m.[49m[43mcategories[49m[43m.[49m[43m_values[49m[43m,[49m[43m [49m[43mdataset_kwargs[49m[43m=[49m[43mdataset_kwargs[49m
[32m   1103[39m [43m[49m[43m)[49m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/utils.py:243[39m, in [36mreport_write_key_on_error.<locals>.func_wrapper[39m[34m(*args, **kwargs)[39m
[32m    241[39m     [38;5;28;01mraise[39;00m [38;5;167;01mValueError[39;00m(msg)
[32m    242[39m [38;5;28;01mtry[39;00m:
[32m--> [39m[32m243[39m     [38;5;28;01mreturn[39;00m [43mfunc[49m[43m([49m[43m*[49m[43margs[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkwargs[49m[43m)[49m
[32m    244[39m [38;5;28;01mexcept[39;00m [38;5;167;01mException[39;00m [38;5;28;01mas[39;00m e:
[32m    245[39m     path = _get_display_path(store)

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/specs/registry.py:389[39m, in [36mWriter.write_elem[39m[34m(self, store, k, elem, dataset_kwargs, modifiers)[39m
[32m    386[39m write_func = [38;5;28mself[39m.find_write_func(dest_type, elem, modifiers)
[32m    388[39m [38;5;28;01mif[39;00m [38;5;28mself[39m.callback [38;5;129;01mis[39;00m [38;5;28;01mNone[39;00m:
[32m--> [39m[32m389[39m     [38;5;28;01mreturn[39;00m [43mwrite_func[49m[43m([49m[43mstore[49m[43m,[49m[43m [49m[43mk[49m[43m,[49m[43m [49m[43melem[49m[43m,[49m[43m [49m[43mdataset_kwargs[49m[43m=[49m[43mdataset_kwargs[49m[43m)[49m
[32m    390[39m [38;5;28;01mreturn[39;00m [38;5;28mself[39m.callback(
[32m    391[39m     write_func,
[32m    392[39m     store,
[32m   (...)[39m[32m    396[39m     iospec=[38;5;28mself[39m.registry.get_spec(elem),
[32m    397[39m )

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/specs/registry.py:77[39m, in [36mwrite_spec.<locals>.decorator.<locals>.wrapper[39m[34m(g, k, *args, **kwargs)[39m
[32m     75[39m [38;5;129m@wraps[39m(func)
[32m     76[39m [38;5;28;01mdef[39;00m[38;5;250m [39m[34mwrapper[39m(g: GroupStorageType, k: [38;5;28mstr[39m, *args, **kwargs):
[32m---> [39m[32m77[39m     result = [43mfunc[49m[43m([49m[43mg[49m[43m,[49m[43m [49m[43mk[49m[43m,[49m[43m [49m[43m*[49m[43margs[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkwargs[49m[43m)[49m
[32m     78[39m     g[k].attrs.setdefault([33m"[39m[33mencoding-type[39m[33m"[39m, spec.encoding_type)
[32m     79[39m     g[k].attrs.setdefault([33m"[39m[33mencoding-version[39m[33m"[39m, spec.encoding_version)

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/anndata/_io/specs/methods.py:1167[39m, in [36mwrite_nullable[39m[34m(f, k, v, _writer, dataset_kwargs)[39m
[32m   1157[39m [38;5;28;01mif[39;00m (
[32m   1158[39m     [38;5;28misinstance[39m(v, pd.arrays.StringArray)
[32m   1159[39m     [38;5;129;01mand[39;00m [38;5;129;01mnot[39;00m settings.allow_write_nullable_strings
[32m   1160[39m ):
[32m   1161[39m     msg = (
[32m   1162[39m         [33m"[39m[33m`anndata.settings.allow_write_nullable_strings` is False, [39m[33m"[39m
[32m   1163[39m         [33m"[39m[33mbecause writing of `pd.arrays.StringArray` is new [39m[33m"[39m
[32m   1164[39m         [33m"[39m[33mand not supported in anndata < 0.11, still use by many people. [39m[33m"[39m
[32m   1165[39m         [33m"[39m[33mOpt-in to writing these arrays by toggling the setting to True.[39m[33m"[39m
[32m   1166[39m     )
[32m-> [39m[32m1167[39m     [38;5;28;01mraise[39;00m [38;5;167;01mRuntimeError[39;00m(msg)
[32m   1168[39m g = f.require_group(k)
[32m   1169[39m values = (
[32m   1170[39m     v.to_numpy(na_value=[33m"[39m[33m"[39m)
[32m   1171[39m     [38;5;28;01mif[39;00m [38;5;28misinstance[39m(v, pd.arrays.StringArray)
[32m   1172[39m     [38;5;28;01melse[39;00m v.to_numpy(na_value=[32m0[39m, dtype=v.dtype.numpy_dtype)
[32m   1173[39m )

[31mRuntimeError[39m: `anndata.settings.allow_write_nullable_strings` is False, because writing of `pd.arrays.StringArray` is new and not supported in anndata < 0.11, still use by many people. Opt-in to writing these arrays by toggling the setting to True.
Error raised while writing key 'categories' of <class 'h5py._hl.group.Group'> to /obs/sample_x_inferred_processing_batch
```

## Code cell 1

```python
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
from IPython.display import display

# ------------------------------------------------------------------
# 1. Load the AnnData object containing inferred 2/4-group labels
# ------------------------------------------------------------------

H5AD_PATH = (
    "/oak/stanford/groups/engreitz/Users/tri/tmp/"
    "Maggie_CRT_sensitivity_test/"
    "rna_with_guides_batch_sensitivity.h5ad"
)

# backed="r" avoids loading the full expression matrix into memory.
adata = ad.read_h5ad(H5AD_PATH, backed="r")

print(adata)
print(f"Cells: {adata.n_obs:,}")
print("\nInferred-group columns:")
print([c for c in adata.obs.columns if "inferred" in c.lower()])

# ------------------------------------------------------------------
# 2. Select the inferred 4-group label
# ------------------------------------------------------------------

inferred_candidates = [
    "inferred_batch_4",
    "inferred_processing_batch_4",
    "sample_x_inferred_batch_4",
]

inferred_col = next(
    (c for c in inferred_candidates if c in adata.obs.columns),
    None,
)

if inferred_col is None:
    raise KeyError(
        "Could not find an inferred 4-group column.\n"
        "Available inferred columns:\n"
        f"{[c for c in adata.obs.columns if 'inferred' in c.lower()]}"
    )

print(f"\nUsing inferred grouping column: {inferred_col}")

# Prefer the inferred-only label over sample × inferred label, if available.
if inferred_col.startswith("sample_x_") and "inferred_batch_4" in adata.obs.columns:
    inferred_col = "inferred_batch_4"
    print(f"Switched to inferred-only label: {inferred_col}")

# ------------------------------------------------------------------
# 3. Find possible biological subtype / Leiden annotations
# ------------------------------------------------------------------

candidate_bio_cols = [
    "cell_type",
    "celltype",
    "cell_subtype",
    "subtype",
    "annotation",
    "broad_cell_type",
    "leiden",
    "leiden_res0_25",
    "leiden_res0_5",
    "leiden_res1",
]

bio_cols = [
    col for col in candidate_bio_cols
    if col in adata.obs.columns
]

if not bio_cols:
    print("\nNo standard annotation names were found.")
    print("Columns containing 'cell', 'type', 'leiden', or 'cluster':")
    print([
        c for c in adata.obs.columns
        if any(term in c.lower() for term in ["cell", "type", "leiden", "cluster"])
    ])
    raise KeyError(
        "Set candidate_bio_cols manually using the appropriate "
        "cell-type or Leiden column."
    )

print("\nBiological annotation columns found:")
print(bio_cols)

# ------------------------------------------------------------------
# 4. Association metrics
# ------------------------------------------------------------------

def bias_corrected_cramers_v(table: pd.DataFrame) -> float:
    """Bias-corrected Cramér's V for two categorical variables."""
    observed = table.to_numpy()
    n = observed.sum()

    if n <= 1 or min(observed.shape) <= 1:
        return np.nan

    chi2 = chi2_contingency(observed)[0]
    rows, cols = observed.shape
    phi2 = chi2 / n

    phi2_corrected = max(
        0.0,
        phi2 - ((cols - 1) * (rows - 1)) / (n - 1),
    )
    rows_corrected = rows - ((rows - 1) ** 2) / (n - 1)
    cols_corrected = cols - ((cols - 1) ** 2) / (n - 1)

    denominator = min(rows_corrected - 1, cols_corrected - 1)

    if denominator <= 0:
        return np.nan

    return float(np.sqrt(phi2_corrected / denominator))


results = []

for bio_col in bio_cols:
    df = (
        adata.obs[[inferred_col, bio_col]]
        .dropna()
        .copy()
    )

    # Convert to ordinary Python strings for compatibility.
    inferred_values = df[inferred_col].astype(object).astype(str)
    biological_values = df[bio_col].astype(object).astype(str)

    count_table = pd.crosstab(
        biological_values,
        inferred_values,
    )

    # Distribution of inferred groups within each biological annotation.
    inferred_within_biology = pd.crosstab(
        biological_values,
        inferred_values,
        normalize="index",
    )

    # Biological composition within each inferred group.
    biology_within_inferred = pd.crosstab(
        biological_values,
        inferred_values,
        normalize="columns",
    )

    nmi = normalized_mutual_info_score(
        biological_values,
        inferred_values,
    )

    ari = adjusted_rand_score(
        biological_values,
        inferred_values,
    )

    cramers_v = bias_corrected_cramers_v(count_table)

    # Proportion captured by the dominant inferred group within each
    # biological annotation.
    weighted_purity = (
        count_table.max(axis=1).sum()
        / count_table.to_numpy().sum()
    )

    # Number of inferred groups representing at least 5% of each subtype.
    groups_present = (inferred_within_biology >= 0.05).sum(axis=1)
    fraction_with_three_or_more = float(
        np.mean(groups_present >= 3)
    )

    results.append({
        "annotation": bio_col,
        "n_annotation_levels": count_table.shape[0],
        "NMI": nmi,
        "ARI": ari,
        "Cramers_V": cramers_v,
        "weighted_purity": weighted_purity,
        "fraction_levels_with_3plus_inferred_groups": (
            fraction_with_three_or_more
        ),
    })

    print("\n" + "=" * 90)
    print(f"Biological annotation: {bio_col}")
    print("=" * 90)
    print(
        f"NMI={nmi:.3f} | "
        f"ARI={ari:.3f} | "
        f"Cramér's V={cramers_v:.3f} | "
        f"weighted purity={weighted_purity:.3f}"
    )

    print("\nInferred-group proportions within each biological group:")
    display(inferred_within_biology.round(3))

    print("\nBiological composition within each inferred group:")
    display(biology_within_inferred.round(3))

# ------------------------------------------------------------------
# 5. Overall summary
# ------------------------------------------------------------------

summary = (
    pd.DataFrame(results)
    .set_index("annotation")
    .sort_values("Cramers_V", ascending=False)
)

print("\n" + "=" * 90)
print("SUMMARY")
print("=" * 90)
display(summary.round(3))

print(
    "\nInterpretation:\n"
    "- Low NMI/Cramér's V and several inferred groups within most cell types "
    "supports a technical effect.\n"
    "- High NMI/Cramér's V, high purity, or one inferred group per subtype "
    "suggests the inferred groups may reproduce biology.\n"
    "- Intermediate association means the inferred groups likely mix "
    "technical and biological variation."
)
```

### Text output (stream)

```text
AnnData object with n_obs × n_vars = 240504 × 25788 backed at '/oak/stanford/groups/engreitz/Users/tri/tmp/Maggie_CRT_sensitivity_test/rna_with_guides_batch_sensitivity.h5ad'
    obs: 'total_counts', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'outlier', 'sublib', 'S_score', 'G2M_score', 'phase', 'leiden_res0_10', 'leiden_res0_25', 'leiden_res0_50', 'leiden_res1_00', 'ipsc_score', 'vic_score', 'ec_score', 'sample', 'guide_umis_per_cell', 'guides_per_cell', 'technical_pc1', 'technical_pc2', 'technical_pc3', 'technical_pc4', 'inferred_batch_2', 'inferred_batch_2_numeric', 'inferred_batch_2_certainty', 'sample_x_inferred_batch_2', 'inferred_batch_4', 'inferred_batch_4_numeric', 'inferred_batch_4_certainty', 'sample_x_inferred_batch_4'
    var: 'gene_id', 'symbol', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'highly_variable', 'highly_variable_rank', 'means', 'variances', 'variances_norm', 'highly_variable_nbatches', 'var_names'
    uns: 'guide_names', 'guide_targets', 'hvg', 'leiden_res0_10', 'leiden_res0_10_colors', 'leiden_res0_25', 'leiden_res0_25_colors', 'leiden_res0_50', 'leiden_res0_50_colors', 'leiden_res1_00', 'leiden_res1_00_colors', 'log1p', 'neighbors', 'pca', 'phase_colors', 'sublib_colors', 'umap'
    obsm: 'X_PCA', 'X_umap', 'guide_assignment'
    varm: 'PCs'
    layers: 'ambiguous', 'counts', 'mature', 'nascent'
    obsp: 'connectivities', 'distances'
Cells: 240,504

Inferred-group columns:
['inferred_batch_2', 'inferred_batch_2_numeric', 'inferred_batch_2_certainty', 'sample_x_inferred_batch_2', 'inferred_batch_4', 'inferred_batch_4_numeric', 'inferred_batch_4_certainty', 'sample_x_inferred_batch_4']

Using inferred grouping column: inferred_batch_4

Biological annotation columns found:
['leiden_res0_25']

==========================================================================================
Biological annotation: leiden_res0_25
==========================================================================================
NMI=0.026 | ARI=0.030 | Cramér's V=0.162 | weighted purity=0.447

Inferred-group proportions within each biological group:

```

### Text output (text/plain)

```text
inferred_batch_4  batch_0  batch_1  batch_2  batch_3
leiden_res0_25                                      
0                   0.200    0.224    0.113    0.463
1                   0.208    0.272    0.086    0.434
2                   0.220    0.216    0.346    0.218
3                   0.150    0.194    0.142    0.515
4                   0.244    0.161    0.315    0.280
```

### Text output (stream)

```text

Biological composition within each inferred group:

```

### Text output (text/plain)

```text
inferred_batch_4  batch_0  batch_1  batch_2  batch_3
leiden_res0_25                                      
0                   0.332    0.322    0.238    0.345
1                   0.261    0.296    0.136    0.244
2                   0.134    0.114    0.266    0.059
3                   0.203    0.228    0.244    0.315
4                   0.070    0.040    0.115    0.036
```

### Text output (stream)

```text

==========================================================================================
SUMMARY
==========================================================================================

```

### Text output (text/plain)

```text
                n_annotation_levels    NMI   ARI  Cramers_V  weighted_purity  \
annotation                                                                     
leiden_res0_25                    5  0.026  0.03      0.162            0.447   

                fraction_levels_with_3plus_inferred_groups  
annotation                                                  
leiden_res0_25                                         1.0  
```

### Text output (stream)

```text

Interpretation:
- Low NMI/Cramér's V and several inferred groups within most cell types supports a technical effect.
- High NMI/Cramér's V, high purity, or one inferred group per subtype suggests the inferred groups may reproduce biology.
- Intermediate association means the inferred groups likely mix technical and biological variation.

```

