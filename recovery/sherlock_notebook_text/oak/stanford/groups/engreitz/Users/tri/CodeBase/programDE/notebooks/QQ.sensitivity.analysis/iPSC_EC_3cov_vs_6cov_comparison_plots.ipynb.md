# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/iPSC_EC_3cov_vs_6cov_comparison_plots.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/iPSC_EC_3cov_vs_6cov_comparison_plots.ipynb`
- SHA-256: `0eb27f02b992df93fc3d9c4c15e173410cd155eea9e8ca83a554a8c4fcdbf9f1`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0


# iPSC-EC: three-covariate versus six-covariate CRT comparison

This notebook compares two completed iPSC-EC `programDE` runs that used the same cells, cNMF programs, guide assignments, CRT settings, and time points but different covariate models.

**Three-covariate model**

- `log1p(guides_per_cell)`
- `log1p(total_counts)`
- `pct_counts_mt`

**Six-source-covariate model**

- `log1p(guides_per_cell)`
- `log1p(total_counts)`
- `pct_counts_mt`
- `log1p(guide_umi_counts)`
- `log1p(n_genes_by_counts)`
- `biological_sample` (dummy encoded; therefore the actual number of design-matrix columns can exceed six)

The comparison is intended to answer two separate questions:

1. **Does the larger model improve null calibration?** This is assessed using NTC pseudo-target p-values: nominal type-I-error rates, QQ error, and tail behavior.
2. **How much do the added covariates change biological conclusions?** This is assessed using beta concordance, effect attenuation/amplification, p-value rank concordance, discovery overlap, sign agreement, Bland–Altman differences, and program-/target-level heterogeneity.

The main metrics produced are:

- Pearson, Spearman, concordance correlation coefficient, regression slope, RMSE, and sign concordance for target–program beta effects;
- median absolute effect change and six-covariate/three-covariate effect-magnitude ratio;
- p-value rank concordance and changes in `-log10(p)`;
- discovery counts, intersection, model-specific discoveries, and Jaccard overlap at `q <= 0.05` and `q <= 0.10`;
- NTC fractions below 0.01, 0.05, and 0.10, uniform QQ RMSE, tail QQ RMSE, and KS distance from uniformity;
- per-program and per-target beta concordance;
- optional comparison of biological-replicate correlations if replicate summary files exist under both model roots.

All plots are saved as **PDF and SVG** for Illustrator editing, plus PNG previews. EPS is deliberately omitted because transparency is not supported reliably by the PostScript backend.


## Code cell 1

```python

from pathlib import Path
import json
import math
import re
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from IPython.display import display

import matplotlib as mpl
from plotnine import (
    ggplot, aes,
    geom_abline, geom_hline, geom_point, geom_line, geom_col,
    geom_tile, geom_boxplot, geom_jitter,
    facet_wrap, coord_equal,
    scale_color_manual, scale_fill_manual, scale_fill_gradient,
    scale_x_continuous, scale_y_continuous, scale_x_log10,
    labs, theme_bw, theme,
    element_text, element_line, element_rect, element_blank,
    position_dodge,
)

# Preserve editable text in Illustrator-compatible vector files.
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"

# ---------------------------------------------------------------------
# Sherlock result paths from the earlier notebooks
# ---------------------------------------------------------------------


THREE_COV_ROOT = Path(
    "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/iPSC_EC_QQ_export_all_days_bundle/results/ipsc_ec_qq_exports/k_50.dt_0_2"
).resolve()

SIX_COV_ROOT = Path(
    "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/iPSC_EC_QQ_export_full_covariates_notebook_bundle/results/ipsc_ec_qq_exports_full_covariates/k_50.dt_0_2"
).resolve()

MODEL_ROOTS = {
    "3cov": THREE_COV_ROOT,
    "6cov": SIX_COV_ROOT,
}

MODEL_LABELS = {
    "3cov": "3 covariates",
    "6cov": "6 covariates",
}

DAYS = ["d0", "d1", "d2", "d3"]
NTC_LABELS = ["SAFE", "non-targeting", "NTC"]
NTC_KEYS = {str(x).strip().lower() for x in NTC_LABELS}

# Compare all p-value types available in both runs. The preferred type is
# used for the main discovery and scatterplot figures.
P_VALUE_TYPES_TO_TRY = ["skew", "empirical"]
PREFERRED_PVALUE_TYPE = "empirical"

Q_CUTOFFS = [0.05, 0.10]
PRIMARY_Q_CUTOFF = 0.10
TOP_N_EFFECTS = [100, 250, 500, 1000, 2000, 5000]

# Set True to stop when run configurations differ in a critical field.
# False still prints and saves the mismatch audit.
STRICT_COMPARABILITY = False

COMPARISON_ROOT = Path(
    "results/ipsc_ec_3cov_vs_6cov_comparison/k_50.dt_0_2"
).resolve()
TABLE_DIR = COMPARISON_ROOT / "tables"
PLOT_DIR = COMPARISON_ROOT / "plots"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

# Plot colors.
MODEL_COLORS = {
    "3 covariates": "#1F5A94",
    "6 covariates": "#B22222",
}
DISCOVERY_COLORS = {
    "Neither": "#BFC5CC",
    "Both": "#17365D",
    "3 covariates only": "#D98E04",
    "6 covariates only": "#7A3E9D",
}
GRID_COLOR = "#E7ECF1"
BORDER_COLOR = "#1D2A36"
DIAGONAL_COLOR = "#AAB3BD"

print("Three-covariate root:", THREE_COV_ROOT)
print("Six-covariate root:  ", SIX_COV_ROOT)
print("Comparison output:   ", COMPARISON_ROOT)

```

### Text output (stream)

```text
Three-covariate root: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/iPSC_EC_QQ_export_all_days_bundle/results/ipsc_ec_qq_exports/k_50.dt_0_2
Six-covariate root:   /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/iPSC_EC_QQ_export_full_covariates_notebook_bundle/results/ipsc_ec_qq_exports_full_covariates/k_50.dt_0_2
Comparison output:    /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/results/ipsc_ec_3cov_vs_6cov_comparison/k_50.dt_0_2

```

## Markdown cell 2

## Helper functions

## Code cell 3

```python

def natural_sort_key(value):
    return [
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", str(value))
    ]


def read_json(path):
    path = Path(path)
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def json_stable(value):
    if isinstance(value, (list, tuple, dict)):
        return json.dumps(value, sort_keys=True)
    if value is None:
        return None
    return str(value)


def bh_qvalues(values):
    p = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p)
    if finite.sum() == 0:
        return q

    pv = np.clip(p[finite], 0.0, 1.0)
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]

    q_finite = np.empty_like(pv)
    q_finite[order] = np.clip(ranked, 0.0, 1.0)
    q[finite] = q_finite
    return q


def safe_pearson(x, y):
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def safe_spearman(x, y):
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    if len(x) < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan
    return float(spearmanr(x, y).statistic)


def concordance_correlation_coefficient(x, y):
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    if len(x) < 3:
        return np.nan

    var_x = np.var(x, ddof=1)
    var_y = np.var(y, ddof=1)
    covariance = np.cov(x, y, ddof=1)[0, 1]
    denominator = var_x + var_y + (np.mean(x) - np.mean(y)) ** 2
    if denominator == 0:
        return np.nan
    return float(2.0 * covariance / denominator)


def regression_slope_intercept(x, y):
    x = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    if len(x) < 3 or np.std(x) == 0:
        return np.nan, np.nan
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(intercept)


def jaccard_index(a, b):
    a = np.asarray(a, dtype=bool)
    b = np.asarray(b, dtype=bool)
    union = np.sum(a | b)
    return 1.0 if union == 0 else float(np.sum(a & b) / union)


def read_target_wide(path, value_name):
    path = Path(path)
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str).str.strip()
    frame.columns = frame.columns.astype(str).str.strip()

    if frame.index.duplicated().any():
        duplicates = frame.index[frame.index.duplicated()].unique().tolist()[:10]
        raise ValueError(f"Duplicate target names in {path}: {duplicates}")
    if frame.columns.duplicated().any():
        duplicates = frame.columns[frame.columns.duplicated()].unique().tolist()[:10]
        raise ValueError(f"Duplicate program names in {path}: {duplicates}")

    frame.index.name = "target_name"
    long = (
        frame.reset_index()
        .melt(
            id_vars="target_name",
            var_name="program_name",
            value_name=value_name,
        )
    )
    long[value_name] = pd.to_numeric(long[value_name], errors="coerce")
    long = long[np.isfinite(long[value_name])].copy()
    long = long[
        ~long["target_name"].str.strip().str.lower().isin(NTC_KEYS)
    ].copy()
    return long


def read_ntc_wide(path):
    path = Path(path)
    frame = pd.read_csv(path, sep="\t")
    id_columns = [c for c in ["ensemble", "group_id"] if c in frame.columns]
    if not id_columns:
        raise ValueError(
            f"Expected ensemble/group_id columns in NTC file {path}; "
            f"found {list(frame.columns)}"
        )

    value_columns = [c for c in frame.columns if c not in id_columns]
    long = frame.melt(
        id_vars=id_columns,
        value_vars=value_columns,
        var_name="program_name",
        value_name="pvalue",
    )
    long["program_name"] = long["program_name"].astype(str).str.strip()
    long["pvalue"] = pd.to_numeric(long["pvalue"], errors="coerce")
    long = long[np.isfinite(long["pvalue"])].copy()
    long["pvalue"] = long["pvalue"].clip(1e-300, 1.0)
    return long


def day_paths(root, day):
    day_dir = Path(root) / day
    return {
        "day_dir": day_dir,
        "config": day_dir / f"{day}_run_config.json",
        "beta": day_dir / f"{day}_target_betas_wide.tsv",
        "target_empirical": day_dir / f"{day}_target_pvalues_empirical_wide.tsv",
        "target_skew": day_dir / f"{day}_target_pvalues_skew_wide.tsv",
        "ntc_empirical": day_dir / f"{day}_ntc_pseudotarget_pvalues_empirical_wide.tsv",
        "ntc_skew": day_dir / f"{day}_ntc_pseudotarget_pvalues_skew_wide.tsv",
    }


def uniform_calibration_metrics(pvalues):
    p = pd.to_numeric(pd.Series(pvalues), errors="coerce").to_numpy(dtype=float)
    p = np.clip(p[np.isfinite(p)], 1e-300, 1.0)
    if len(p) == 0:
        return {
            "n_pvalues": 0,
            "median_p": np.nan,
            "mean_p": np.nan,
            "fraction_p_lt_0_01": np.nan,
            "fraction_p_lt_0_05": np.nan,
            "fraction_p_lt_0_10": np.nan,
            "abs_error_at_0_05": np.nan,
            "ks_distance": np.nan,
            "qq_rmse": np.nan,
            "qq_tail_rmse": np.nan,
            "qq_max_abs_error": np.nan,
        }

    p_sorted = np.sort(p)
    n = len(p_sorted)
    expected = (np.arange(1, n + 1) - 0.5) / n

    observed_log = -np.log10(p_sorted)
    expected_log = -np.log10(expected)
    qq_error = observed_log - expected_log
    tail_n = max(10, int(np.ceil(0.01 * n)))

    # One-sample Kolmogorov distance from Uniform(0,1). The p-value is not
    # reported because NTC pseudo-target tests can be dependent.
    i = np.arange(1, n + 1)
    d_plus = np.max(i / n - p_sorted)
    d_minus = np.max(p_sorted - (i - 1) / n)

    fraction_005 = float((p < 0.05).mean())
    return {
        "n_pvalues": int(n),
        "median_p": float(np.median(p)),
        "mean_p": float(np.mean(p)),
        "fraction_p_lt_0_01": float((p < 0.01).mean()),
        "fraction_p_lt_0_05": fraction_005,
        "fraction_p_lt_0_10": float((p < 0.10).mean()),
        "abs_error_at_0_05": abs(fraction_005 - 0.05),
        "ks_distance": float(max(d_plus, d_minus)),
        "qq_rmse": float(np.sqrt(np.mean(qq_error ** 2))),
        "qq_tail_rmse": float(np.sqrt(np.mean(qq_error[:tail_n] ** 2))),
        "qq_max_abs_error": float(np.max(np.abs(qq_error))),
    }


def qq_points(pvalues):
    p = pd.to_numeric(pd.Series(pvalues), errors="coerce").to_numpy(dtype=float)
    p = np.clip(p[np.isfinite(p)], 1e-300, 1.0)
    p = np.sort(p)
    n = len(p)
    expected = (np.arange(1, n + 1) - 0.5) / n
    return pd.DataFrame({
        "expected_neglog10_p": -np.log10(expected),
        "observed_neglog10_p": -np.log10(p),
    })


def save_plot(p, output_stem, width, height, dpi=300):
    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)

    outputs = {}
    for extension in ["pdf", "svg", "png"]:
        path = output_stem.with_suffix(f".{extension}")
        kwargs = {
            "width": width,
            "height": height,
            "units": "in",
            "verbose": False,
        }
        if extension == "png":
            kwargs["dpi"] = dpi
        p.save(path, **kwargs)
        outputs[extension] = path
    return outputs


def format_metric(value, digits=2):
    if value is None or not np.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"

```

## Markdown cell 4

## 1. Audit files, run settings, and covariate columns

## Code cell 5

```python

for model_code, root in MODEL_ROOTS.items():
    if not root.exists():
        raise FileNotFoundError(
            f"{MODEL_LABELS[model_code]} result root does not exist: {root}\n"
            "Edit THREE_COV_ROOT or SIX_COV_ROOT in the configuration cell."
        )

manifest_lookup = {}
for model_code, root in MODEL_ROOTS.items():
    manifest_path = root / "QQ_export_manifest.csv"
    if manifest_path.exists():
        manifest = pd.read_csv(manifest_path)
        manifest["day"] = manifest["day"].astype(str)
        manifest_lookup[model_code] = manifest.set_index("day").to_dict("index")
    else:
        manifest_lookup[model_code] = {}

file_rows = []
config_rows = []
covariate_rows = []
config_by_model_day = {}
paths_by_model_day = {}

for model_code, root in MODEL_ROOTS.items():
    for day in DAYS:
        paths = day_paths(root, day)
        paths_by_model_day[(model_code, day)] = paths
        config = read_json(paths["config"])
        config_by_model_day[(model_code, day)] = config
        manifest_row = manifest_lookup.get(model_code, {}).get(day, {})

        file_rows.append({
            "model": model_code,
            "model_label": MODEL_LABELS[model_code],
            "day": day,
            "root": str(root),
            "config_exists": paths["config"].exists(),
            "beta_exists": paths["beta"].exists(),
            "target_empirical_exists": paths["target_empirical"].exists(),
            "target_skew_exists": paths["target_skew"].exists(),
            "ntc_empirical_exists": paths["ntc_empirical"].exists(),
            "ntc_skew_exists": paths["ntc_skew"].exists(),
        })

        config_rows.append({
            "model": model_code,
            "model_label": MODEL_LABELS[model_code],
            "day": day,
            "k": config.get("k"),
            "density_threshold_token": config.get("density_threshold_token"),
            "B": config.get("B"),
            "h5mu_path": config.get("h5mu_path"),
            "usage_path": config.get("usage_path"),
            "guide_long_csv": config.get("guide_long_csv"),
            "gex_modality": config.get("gex_modality"),
            "guide_modality": config.get("guide_modality"),
            "ntc_labels": json_stable(config.get("ntc_labels")),
            "ntc_group_size": config.get("ntc_group_size"),
            "ntc_ensembles": config.get("ntc_ensembles"),
            "ntc_bins": config.get("ntc_bins"),
            "ntc_max_groups": config.get("ntc_max_groups"),
            "ntc_group_seed": config.get("ntc_group_seed"),
            "ntc_crt_seed": config.get("ntc_crt_seed"),
            "n_cells": manifest_row.get("n_cells"),
            "n_programs": manifest_row.get("n_programs"),
            "n_guides": manifest_row.get("n_guides"),
            "matched_usage_cells": manifest_row.get("matched_usage_cells"),
            "n_covariate_columns": len(config.get("covariate_columns", [])),
        })

        columns = config.get("covariate_columns", [])
        sources = config.get("covariate_sources", {})
        for position, column in enumerate(columns, start=1):
            source = sources.get(column, {}) if isinstance(sources, dict) else {}
            covariate_rows.append({
                "model": model_code,
                "model_label": MODEL_LABELS[model_code],
                "day": day,
                "position": position,
                "covariate_column": column,
                "source_obs_column": source.get("obs_key"),
                "transform": source.get("transform"),
            })

file_audit = pd.DataFrame(file_rows)
config_audit = pd.DataFrame(config_rows)
covariate_audit = pd.DataFrame(covariate_rows)

print("File audit:")
display(file_audit)

missing_required = file_audit[
    ~file_audit[["config_exists", "beta_exists", "target_empirical_exists", "ntc_empirical_exists"]].all(axis=1)
]
if not missing_required.empty:
    display(missing_required)
    raise FileNotFoundError(
        "One or more required config/beta/empirical-p-value files are missing. "
        "Correct the model roots or complete the corresponding run."
    )

print("\nActual covariate design-matrix columns:")
display(covariate_audit)

# Validate the expected covariate structures.
expected_3cov = {"grna_n_nonzero", "response_n_umis", "response_p_mito"}
expected_6cov_numeric = {
    "pct_counts_mt",
    "log1p_total_counts",
    "log1p_guides_per_cell",
    "log1p_guide_umi_counts",
    "log1p_n_genes_by_counts",
}

covariate_validation_rows = []
for model_code in MODEL_ROOTS:
    for day in DAYS:
        columns = set(config_by_model_day[(model_code, day)].get("covariate_columns", []))
        biological_columns = sorted(c for c in columns if c.startswith("biological_sample"))
        covariate_validation_rows.append({
            "model": model_code,
            "day": day,
            "n_columns": len(columns),
            "has_exact_expected_3cov": columns == expected_3cov if model_code == "3cov" else np.nan,
            "has_all_expected_6cov_numeric": expected_6cov_numeric.issubset(columns) if model_code == "6cov" else np.nan,
            "n_biological_sample_dummy_columns": len(biological_columns),
            "biological_sample_dummy_columns": ", ".join(biological_columns),
            "columns": ", ".join(sorted(columns, key=natural_sort_key)),
        })

covariate_validation = pd.DataFrame(covariate_validation_rows)
print("\nCovariate validation:")
display(covariate_validation)

# Compare all settings that should be held fixed between models.
critical_fields = [
    "k",
    "density_threshold_token",
    "B",
    "h5mu_path",
    "usage_path",
    "guide_long_csv",
    "gex_modality",
    "guide_modality",
    "ntc_labels",
    "ntc_group_size",
    "ntc_ensembles",
    "ntc_bins",
    "ntc_max_groups",
    "ntc_group_seed",
    "ntc_crt_seed",
    "n_cells",
    "n_programs",
    "n_guides",
    "matched_usage_cells",
]

comparability_rows = []
for day in DAYS:
    row3 = config_audit[(config_audit["model"] == "3cov") & (config_audit["day"] == day)].iloc[0]
    row6 = config_audit[(config_audit["model"] == "6cov") & (config_audit["day"] == day)].iloc[0]
    for field in critical_fields:
        value3 = row3[field]
        value6 = row6[field]
        both_missing = pd.isna(value3) and pd.isna(value6)
        same = both_missing or (json_stable(value3) == json_stable(value6))
        comparability_rows.append({
            "day": day,
            "field": field,
            "three_covariates": value3,
            "six_covariates": value6,
            "same": bool(same),
        })

comparability_audit = pd.DataFrame(comparability_rows)
print("\nComparability audit; all rows should be same=True except covariates, which are audited separately:")
display(comparability_audit)

mismatches = comparability_audit[~comparability_audit["same"]]
if not mismatches.empty:
    message = (
        "Critical run settings differ between the three- and six-covariate outputs. "
        "Differences cannot be attributed solely to covariate adjustment."
    )
    display(mismatches)
    if STRICT_COMPARABILITY:
        raise ValueError(message)
    warnings.warn(message)

# Determine p-value types available for every day in both models.
target_pvalue_types = [
    ptype
    for ptype in P_VALUE_TYPES_TO_TRY
    if all(
        paths_by_model_day[(model, day)][f"target_{ptype}"].exists()
        for model in MODEL_ROOTS
        for day in DAYS
    )
]
ntc_pvalue_types = [
    ptype
    for ptype in P_VALUE_TYPES_TO_TRY
    if all(
        paths_by_model_day[(model, day)][f"ntc_{ptype}"].exists()
        for model in MODEL_ROOTS
        for day in DAYS
    )
]

if not target_pvalue_types:
    raise RuntimeError("No target p-value type is available for all days in both model roots.")
if not ntc_pvalue_types:
    raise RuntimeError("No NTC p-value type is available for all days in both model roots.")

PRIMARY_PVALUE_TYPE_USED = (
    PREFERRED_PVALUE_TYPE
    if PREFERRED_PVALUE_TYPE in target_pvalue_types
    else target_pvalue_types[0]
)
PRIMARY_NTC_PVALUE_TYPE_USED = (
    PRIMARY_PVALUE_TYPE_USED
    if PRIMARY_PVALUE_TYPE_USED in ntc_pvalue_types
    else ntc_pvalue_types[0]
)

print("\nTarget p-value types compared:", target_pvalue_types)
print("NTC p-value types compared:   ", ntc_pvalue_types)
print("Primary target type:          ", PRIMARY_PVALUE_TYPE_USED)
print("Primary NTC type:             ", PRIMARY_NTC_PVALUE_TYPE_USED)

file_audit.to_csv(TABLE_DIR / "file_audit.csv", index=False)
config_audit.to_csv(TABLE_DIR / "run_config_audit.csv", index=False)
covariate_audit.to_csv(TABLE_DIR / "covariate_columns_audit.csv", index=False)
covariate_validation.to_csv(TABLE_DIR / "covariate_validation.csv", index=False)
comparability_audit.to_csv(TABLE_DIR / "comparability_audit.csv", index=False)

```

### Text output (stream)

```text
File audit:

```

### Text output (text/plain)

```text
  model   model_label day                                               root  \
0  3cov  3 covariates  d0  /oak/stanford/groups/engreitz/Users/tri/CodeBa...   
1  3cov  3 covariates  d1  /oak/stanford/groups/engreitz/Users/tri/CodeBa...   
2  3cov  3 covariates  d2  /oak/stanford/groups/engreitz/Users/tri/CodeBa...   
3  3cov  3 covariates  d3  /oak/stanford/groups/engreitz/Users/tri/CodeBa...   
4  6cov  6 covariates  d0  /oak/stanford/groups/engreitz/Users/tri/CodeBa...   
5  6cov  6 covariates  d1  /oak/stanford/groups/engreitz/Users/tri/CodeBa...   
6  6cov  6 covariates  d2  /oak/stanford/groups/engreitz/Users/tri/CodeBa...   
7  6cov  6 covariates  d3  /oak/stanford/groups/engreitz/Users/tri/CodeBa...   

   config_exists  beta_exists  target_empirical_exists  target_skew_exists  \
0           True         True                     True                True   
1           True         True                     True                True   
2           True         True                     True                True   
3           True         True                     True                True   
4           True         True                     True                True   
5           True         True                     True                True   
6           True         True                     True                True   
7           True         True                     True                True   

   ntc_empirical_exists  ntc_skew_exists  
0                  True             True  
1                  True             True  
2                  True             True  
3                  True             True  
4                  True             True  
5                  True             True  
6                  True             True  
7                  True             True  
```

### Text output (stream)

```text

Actual covariate design-matrix columns:

```

### Text output (text/plain)

```text
   model   model_label day  position           covariate_column  \
0   3cov  3 covariates  d0         1             grna_n_nonzero   
1   3cov  3 covariates  d0         2            response_n_umis   
2   3cov  3 covariates  d0         3            response_p_mito   
3   3cov  3 covariates  d1         1             grna_n_nonzero   
4   3cov  3 covariates  d1         2            response_n_umis   
5   3cov  3 covariates  d1         3            response_p_mito   
6   3cov  3 covariates  d2         1             grna_n_nonzero   
7   3cov  3 covariates  d2         2            response_n_umis   
8   3cov  3 covariates  d2         3            response_p_mito   
9   3cov  3 covariates  d3         1             grna_n_nonzero   
10  3cov  3 covariates  d3         2            response_n_umis   
11  3cov  3 covariates  d3         3            response_p_mito   
12  6cov  6 covariates  d0         1              pct_counts_mt   
13  6cov  6 covariates  d0         2         log1p_total_counts   
14  6cov  6 covariates  d0         3      log1p_guides_per_cell   
15  6cov  6 covariates  d0         4     log1p_guide_umi_counts   
16  6cov  6 covariates  d0         5    log1p_n_genes_by_counts   
17  6cov  6 covariates  d0         6  biological_sample_d0_rep2   
18  6cov  6 covariates  d1         1              pct_counts_mt   
19  6cov  6 covariates  d1         2         log1p_total_counts   
20  6cov  6 covariates  d1         3      log1p_guides_per_cell   
21  6cov  6 covariates  d1         4     log1p_guide_umi_counts   
22  6cov  6 covariates  d1         5    log1p_n_genes_by_counts   
23  6cov  6 covariates  d1         6  biological_sample_d1_rep2   
24  6cov  6 covariates  d2         1              pct_counts_mt   
25  6cov  6 covariates  d2         2         log1p_total_counts   
26  6cov  6 covariates  d2         3      log1p_guides_per_cell   
27  6cov  6 covariates  d2         4     log1p_guide_umi_counts   
28  6cov  6 covariates  d2         5    log1p_n_genes_by_counts   
29  6cov  6 covariates  d2         6  biological_sample_d2_rep2   
30  6cov  6 covariates  d3         1              pct_counts_mt   
31  6cov  6 covariates  d3         2         log1p_total_counts   
32  6cov  6 covariates  d3         3      log1p_guides_per_cell   
33  6cov  6 covariates  d3         4     log1p_guide_umi_counts   
34  6cov  6 covariates  d3         5    log1p_n_genes_by_counts   
35  6cov  6 covariates  d3         6  biological_sample_d3_rep2   

    source_obs_column           transform  
0                None                None  
1                None                None  
2                None                None  
3                None                None  
4                None                None  
5                None                None  
6                None                None  
7                None                None  
8                None                None  
9                None                None  
10               None                None  
11               None                None  
12      pct_counts_mt            identity  
13       total_counts               log1p  
14    guides_per_cell               log1p  
15   guide_umi_counts               log1p  
16  n_genes_by_counts               log1p  
17  biological_sample  one_hot_drop_first  
18      pct_counts_mt            identity  
19       total_counts               log1p  
20    guides_per_cell               log1p  
21   guide_umi_counts               log1p  
22  n_genes_by_counts               log1p  
23  biological_sample  one_hot_drop_first  
24      pct_counts_mt            identity  
25       total_counts               log1p  
26    guides_per_cell               log1p  
27   guide_umi_counts               log1p  
28  n_genes_by_counts               log1p  
29  biological_sample  one_hot_drop_first  
30      pct_counts_mt            identity  
31       total_counts               log1p  
32    guides_per_cell               log1p  
33   guide_umi_counts               log1p  
34  n_genes_by_counts               log1p  
35  biological_sample  one_hot_drop_first  
```

### Text output (stream)

```text

Covariate validation:

```

### Text output (text/plain)

```text
  model day  n_columns has_exact_expected_3cov has_all_expected_6cov_numeric  \
0  3cov  d0          3                    True                           NaN   
1  3cov  d1          3                    True                           NaN   
2  3cov  d2          3                    True                           NaN   
3  3cov  d3          3                    True                           NaN   
4  6cov  d0          6                     NaN                          True   
5  6cov  d1          6                     NaN                          True   
6  6cov  d2          6                     NaN                          True   
7  6cov  d3          6                     NaN                          True   

   n_biological_sample_dummy_columns biological_sample_dummy_columns  \
0                                  0                                   
1                                  0                                   
2                                  0                                   
3                                  0                                   
4                                  1       biological_sample_d0_rep2   
5                                  1       biological_sample_d1_rep2   
6                                  1       biological_sample_d2_rep2   
7                                  1       biological_sample_d3_rep2   

                                             columns  
0   grna_n_nonzero, response_n_umis, response_p_mito  
1   grna_n_nonzero, response_n_umis, response_p_mito  
2   grna_n_nonzero, response_n_umis, response_p_mito  
3   grna_n_nonzero, response_n_umis, response_p_mito  
4  biological_sample_d0_rep2, log1p_guide_umi_cou...  
5  biological_sample_d1_rep2, log1p_guide_umi_cou...  
6  biological_sample_d2_rep2, log1p_guide_umi_cou...  
7  biological_sample_d3_rep2, log1p_guide_umi_cou...  
```

### Text output (stream)

```text

Comparability audit; all rows should be same=True except covariates, which are audited separately:

```

### Text output (text/plain)

```text
   day                    field  \
0   d0                        k   
1   d0  density_threshold_token   
2   d0                        B   
3   d0                h5mu_path   
4   d0               usage_path   
..  ..                      ...   
71  d3             ntc_crt_seed   
72  d3                  n_cells   
73  d3               n_programs   
74  d3                 n_guides   
75  d3      matched_usage_cells   

                                     three_covariates  \
0                                                  50   
1                                                 0_2   
2                                               10000   
3   /oak/stanford/groups/engreitz/Users/opushkar/c...   
4   /oak/stanford/groups/engreitz/Users/ymo/IGVF_c...   
..                                                ...   
71                                                 23   
72                                              44102   
73                                                 50   
74                                               1926   
75                                              44102   

                                       six_covariates  same  
0                                                  50  True  
1                                                 0_2  True  
2                                               10000  True  
3   /oak/stanford/groups/engreitz/Users/opushkar/c...  True  
4   /oak/stanford/groups/engreitz/Users/ymo/IGVF_c...  True  
..                                                ...   ...  
71                                                 23  True  
72                                              44102  True  
73                                                 50  True  
74                                               1926  True  
75                                              44102  True  

[76 rows x 5 columns]
```

### Text output (stream)

```text

Target p-value types compared: ['skew', 'empirical']
NTC p-value types compared:    ['skew', 'empirical']
Primary target type:           empirical
Primary NTC type:              empirical

```

## Markdown cell 6

## 2. Load and pair target–program beta effects and p-values

## Code cell 7

```python

model_tables = []
matrix_dimension_rows = []

for model_code, root in MODEL_ROOTS.items():
    for day in DAYS:
        paths = paths_by_model_day[(model_code, day)]
        beta = read_target_wide(paths["beta"], "beta")

        matrix_dimension_rows.append({
            "model": model_code,
            "day": day,
            "n_targets": beta["target_name"].nunique(),
            "n_programs": beta["program_name"].nunique(),
            "n_target_program_pairs": len(beta),
        })

        for pvalue_type in target_pvalue_types:
            pvalues = read_target_wide(
                paths[f"target_{pvalue_type}"],
                "pvalue",
            )

            current = beta.merge(
                pvalues,
                on=["target_name", "program_name"],
                how="outer",
                validate="one_to_one",
                indicator=True,
            )

            unmatched = current[current["_merge"] != "both"]
            if not unmatched.empty:
                warnings.warn(
                    f"{model_code} {day} {pvalue_type}: "
                    f"{len(unmatched):,} beta/p-value rows did not match."
                )

            current = current[current["_merge"] == "both"].drop(columns="_merge")
            current["pvalue"] = current["pvalue"].clip(1e-300, 1.0)
            current["qvalue"] = bh_qvalues(current["pvalue"])
            current.insert(0, "pvalue_type", pvalue_type)
            current.insert(0, "day", day)
            current.insert(0, "model_label", MODEL_LABELS[model_code])
            current.insert(0, "model", model_code)
            model_tables.append(current)

model_long = pd.concat(model_tables, ignore_index=True, sort=False)
matrix_dimensions = pd.DataFrame(matrix_dimension_rows)

paired_tables = []
pairing_rows = []

for day in DAYS:
    for pvalue_type in target_pvalue_types:
        left = model_long[
            (model_long["model"] == "3cov")
            & (model_long["day"] == day)
            & (model_long["pvalue_type"] == pvalue_type)
        ][["target_name", "program_name", "beta", "pvalue", "qvalue"]].rename(
            columns={
                "beta": "beta_3cov",
                "pvalue": "pvalue_3cov",
                "qvalue": "qvalue_3cov",
            }
        )

        right = model_long[
            (model_long["model"] == "6cov")
            & (model_long["day"] == day)
            & (model_long["pvalue_type"] == pvalue_type)
        ][["target_name", "program_name", "beta", "pvalue", "qvalue"]].rename(
            columns={
                "beta": "beta_6cov",
                "pvalue": "pvalue_6cov",
                "qvalue": "qvalue_6cov",
            }
        )

        paired_day = left.merge(
            right,
            on=["target_name", "program_name"],
            how="outer",
            validate="one_to_one",
            indicator=True,
        )

        counts = paired_day["_merge"].value_counts(dropna=False)
        pairing_rows.append({
            "day": day,
            "pvalue_type": pvalue_type,
            "three_covariates_only": int(counts.get("left_only", 0)),
            "six_covariates_only": int(counts.get("right_only", 0)),
            "matched": int(counts.get("both", 0)),
            "match_fraction": float((paired_day["_merge"] == "both").mean()),
        })

        paired_day = paired_day[paired_day["_merge"] == "both"].drop(columns="_merge")
        paired_day.insert(0, "pvalue_type", pvalue_type)
        paired_day.insert(0, "day", day)
        paired_tables.append(paired_day)

paired = pd.concat(paired_tables, ignore_index=True, sort=False)
pairing_audit = pd.DataFrame(pairing_rows)

paired["neglog10_p_3cov"] = -np.log10(paired["pvalue_3cov"].clip(1e-300, 1.0))
paired["neglog10_p_6cov"] = -np.log10(paired["pvalue_6cov"].clip(1e-300, 1.0))
paired["delta_neglog10_p"] = paired["neglog10_p_6cov"] - paired["neglog10_p_3cov"]
paired["sign_concordant"] = np.sign(paired["beta_3cov"]) == np.sign(paired["beta_6cov"])

# Betas are identical across p-value-type rows; use one copy for effect-only analyses.
paired_effects = (
    paired[paired["pvalue_type"] == PRIMARY_PVALUE_TYPE_USED]
    [["day", "target_name", "program_name", "beta_3cov", "beta_6cov"]]
    .drop_duplicates(["day", "target_name", "program_name"])
    .copy()
)

paired_effects["beta_difference_6_minus_3"] = paired_effects["beta_6cov"] - paired_effects["beta_3cov"]
paired_effects["beta_mean"] = (paired_effects["beta_6cov"] + paired_effects["beta_3cov"]) / 2.0
paired_effects["abs_beta_3cov"] = paired_effects["beta_3cov"].abs()
paired_effects["abs_beta_6cov"] = paired_effects["beta_6cov"].abs()
paired_effects["max_abs_beta"] = paired_effects[["abs_beta_3cov", "abs_beta_6cov"]].max(axis=1)
paired_effects["sign_concordant"] = np.sign(paired_effects["beta_3cov"]) == np.sign(paired_effects["beta_6cov"])

paired_primary = paired[paired["pvalue_type"] == PRIMARY_PVALUE_TYPE_USED].copy()
paired_primary["sign_concordant"] = np.sign(paired_primary["beta_3cov"]) == np.sign(paired_primary["beta_6cov"])

is_3 = paired_primary["qvalue_3cov"] <= PRIMARY_Q_CUTOFF
is_6 = paired_primary["qvalue_6cov"] <= PRIMARY_Q_CUTOFF
paired_primary["discovery_status"] = np.select(
    [is_3 & is_6, is_3 & ~is_6, ~is_3 & is_6],
    ["Both", "3 covariates only", "6 covariates only"],
    default="Neither",
)
paired_primary["union_discovery"] = is_3 | is_6
paired_primary["intersection_discovery"] = is_3 & is_6

# Add the primary discovery status to the effect-only table.
paired_effects = paired_effects.merge(
    paired_primary[
        ["day", "target_name", "program_name", "discovery_status", "union_discovery", "intersection_discovery"]
    ],
    on=["day", "target_name", "program_name"],
    how="left",
    validate="one_to_one",
)

print("Matrix dimensions after excluding NTC-labeled target rows:")
display(matrix_dimensions)

print("Pairing audit:")
display(pairing_audit)

if (pairing_audit["match_fraction"] < 0.99).any():
    warnings.warn(
        "Less than 99% of target-program pairs matched for at least one day/type. "
        "Inspect target/program naming before interpreting comparisons."
    )

print("Paired target-program rows:", paired.shape)
display(paired.head())

matrix_dimensions.to_csv(TABLE_DIR / "matrix_dimensions.csv", index=False)
pairing_audit.to_csv(TABLE_DIR / "target_program_pairing_audit.csv", index=False)
model_long.to_csv(TABLE_DIR / "model_target_program_results_long.csv", index=False)
paired.to_csv(TABLE_DIR / "paired_target_program_results_all_pvalue_types.csv", index=False)
paired_effects.to_csv(TABLE_DIR / "paired_target_program_beta_effects.csv", index=False)

```

### Text output (stream)

```text
Matrix dimensions after excluding NTC-labeled target rows:

```

### Text output (text/plain)

```text
  model day  n_targets  n_programs  n_target_program_pairs
0  3cov  d0        298          50                   14900
1  3cov  d1        298          50                   14900
2  3cov  d2        298          50                   14900
3  3cov  d3        298          50                   14900
4  6cov  d0        298          50                   14900
5  6cov  d1        298          50                   14900
6  6cov  d2        298          50                   14900
7  6cov  d3        298          50                   14900
```

### Text output (stream)

```text
Pairing audit:

```

### Text output (text/plain)

```text
  day pvalue_type  three_covariates_only  six_covariates_only  matched  \
0  d0        skew                      0                    0    14900   
1  d0   empirical                      0                    0    14900   
2  d1        skew                      0                    0    14900   
3  d1   empirical                      0                    0    14900   
4  d2        skew                      0                    0    14900   
5  d2   empirical                      0                    0    14900   
6  d3        skew                      0                    0    14900   
7  d3   empirical                      0                    0    14900   

   match_fraction  
0             1.0  
1             1.0  
2             1.0  
3             1.0  
4             1.0  
5             1.0  
6             1.0  
7             1.0  
```

### Text output (stream)

```text
Paired target-program rows: (119200, 14)

```

### Text output (text/plain)

```text
  day pvalue_type target_name program_name  beta_3cov  pvalue_3cov  \
0  d0        skew       ACAA1            1  -0.020199     0.796331   
1  d0        skew       ACAA1           10  -0.041764     0.702691   
2  d0        skew       ACAA1           11   0.025892     0.736425   
3  d0        skew       ACAA1           12   0.407773     0.001402   
4  d0        skew       ACAA1           13  -0.211062     0.063878   

   qvalue_3cov  beta_6cov  pvalue_6cov  qvalue_6cov  neglog10_p_3cov  \
0     0.952019  -0.036343     0.767171     0.945219         0.098906   
1     0.926557   0.003820     0.868163     0.976353         0.153236   
2     0.936080   0.029119     0.747228     0.939718         0.132871   
3     0.036776   0.399970     0.001522     0.038892         2.853268   
4     0.388396  -0.214176     0.064687     0.394725         1.194646   

   neglog10_p_6cov  delta_neglog10_p  sign_concordant  
0         0.115108          0.016201             True  
1         0.061399         -0.091837            False  
2         0.126547         -0.006325             True  
3         2.817688         -0.035580             True  
4         1.189182         -0.005464             True  
```

## Markdown cell 8

## 3. Quantitative effect, p-value, and discovery metrics

## Code cell 9

```python

def summarize_effect_group(group, selection):
    x = group["beta_3cov"].to_numpy(dtype=float)
    y = group["beta_6cov"].to_numpy(dtype=float)
    delta = y - x
    slope, intercept = regression_slope_intercept(x, y)
    sd3 = np.std(x, ddof=1) if len(x) > 1 else np.nan
    rmse = float(np.sqrt(np.mean(delta ** 2))) if len(delta) else np.nan
    median_abs_3 = float(np.median(np.abs(x))) if len(x) else np.nan
    median_abs_6 = float(np.median(np.abs(y))) if len(y) else np.nan

    nonzero = (x != 0) & (y != 0)
    sign_concordance = (
        float((np.sign(x[nonzero]) == np.sign(y[nonzero])).mean())
        if nonzero.sum() > 0
        else np.nan
    )

    mean_delta = float(np.mean(delta)) if len(delta) else np.nan
    sd_delta = float(np.std(delta, ddof=1)) if len(delta) > 1 else np.nan

    return {
        "selection": selection,
        "n_pairs": int(len(group)),
        "pearson_r": safe_pearson(x, y),
        "spearman_r": safe_spearman(x, y),
        "concordance_correlation_coefficient": concordance_correlation_coefficient(x, y),
        "regression_slope_6_on_3": slope,
        "regression_intercept_6_on_3": intercept,
        "rmse_beta": rmse,
        "normalized_rmse_by_sd_3cov": rmse / sd3 if np.isfinite(sd3) and sd3 > 0 else np.nan,
        "median_abs_beta_change": float(np.median(np.abs(delta))) if len(delta) else np.nan,
        "p95_abs_beta_change": float(np.quantile(np.abs(delta), 0.95)) if len(delta) else np.nan,
        "mean_beta_difference_6_minus_3": mean_delta,
        "bland_altman_lower_95": mean_delta - 1.96 * sd_delta if np.isfinite(sd_delta) else np.nan,
        "bland_altman_upper_95": mean_delta + 1.96 * sd_delta if np.isfinite(sd_delta) else np.nan,
        "sign_concordance": sign_concordance,
        "median_abs_beta_3cov": median_abs_3,
        "median_abs_beta_6cov": median_abs_6,
        "median_abs_effect_ratio_6_over_3": median_abs_6 / median_abs_3 if median_abs_3 > 0 else np.nan,
    }


effect_summary_rows = []
for day, group in paired_effects.groupby("day", sort=False):
    row = summarize_effect_group(group, "all target-program pairs")
    row["day"] = day
    effect_summary_rows.append(row)

    union = group[group["union_discovery"]].copy()
    if not union.empty:
        row = summarize_effect_group(
            union,
            f"union discoveries at {PRIMARY_PVALUE_TYPE_USED} q <= {PRIMARY_Q_CUTOFF}",
        )
        row["day"] = day
        effect_summary_rows.append(row)

    threshold = group["max_abs_beta"].quantile(0.90)
    top10 = group[group["max_abs_beta"] >= threshold].copy()
    row = summarize_effect_group(top10, "top 10% by max absolute beta")
    row["day"] = day
    effect_summary_rows.append(row)

effect_summary = pd.DataFrame(effect_summary_rows)

pvalue_summary_rows = []
for (day, pvalue_type), group in paired.groupby(["day", "pvalue_type"], sort=False):
    pvalue_summary_rows.append({
        "day": day,
        "pvalue_type": pvalue_type,
        "n_pairs": len(group),
        "pearson_neglog10_p": safe_pearson(group["neglog10_p_3cov"], group["neglog10_p_6cov"]),
        "spearman_neglog10_p": safe_spearman(group["neglog10_p_3cov"], group["neglog10_p_6cov"]),
        "median_abs_change_neglog10_p": float(np.median(np.abs(group["delta_neglog10_p"]))),
        "p95_abs_change_neglog10_p": float(np.quantile(np.abs(group["delta_neglog10_p"]), 0.95)),
        "median_change_neglog10_p_6_minus_3": float(np.median(group["delta_neglog10_p"])),
        "fraction_more_significant_in_6cov": float((group["pvalue_6cov"] < group["pvalue_3cov"]).mean()),
    })
pvalue_summary = pd.DataFrame(pvalue_summary_rows)

discovery_rows = []
for (day, pvalue_type), group in paired.groupby(["day", "pvalue_type"], sort=False):
    for q_cutoff in Q_CUTOFFS:
        sig3 = group["qvalue_3cov"] <= q_cutoff
        sig6 = group["qvalue_6cov"] <= q_cutoff
        intersection = sig3 & sig6
        union = sig3 | sig6

        sign_intersection = (
            float(group.loc[intersection, "sign_concordant"].mean())
            if intersection.sum() > 0
            else np.nan
        )
        sign_union = (
            float(group.loc[union, "sign_concordant"].mean())
            if union.sum() > 0
            else np.nan
        )

        discovery_rows.append({
            "day": day,
            "pvalue_type": pvalue_type,
            "q_cutoff": q_cutoff,
            "n_discoveries_3cov": int(sig3.sum()),
            "n_discoveries_6cov": int(sig6.sum()),
            "n_intersection": int(intersection.sum()),
            "n_3cov_only": int((sig3 & ~sig6).sum()),
            "n_6cov_only": int((~sig3 & sig6).sum()),
            "n_union": int(union.sum()),
            "jaccard": jaccard_index(sig3, sig6),
            "sign_concordance_intersection": sign_intersection,
            "sign_concordance_union": sign_union,
        })
discovery_summary = pd.DataFrame(discovery_rows)

# Agreement as a function of effect strength. Ranking is symmetric between models:
# max(|beta_3cov|, |beta_6cov|), so neither model alone defines the subset.
topn_rows = []
for day, group in paired_effects.groupby("day", sort=False):
    ordered = group.sort_values("max_abs_beta", ascending=False)
    requested = TOP_N_EFFECTS + [len(ordered)]
    for top_n in sorted(set(min(n, len(ordered)) for n in requested)):
        selected = ordered.head(top_n)
        summary = summarize_effect_group(selected, f"top {top_n}")
        topn_rows.append({
            "day": day,
            "top_n": int(top_n),
            "is_all_pairs": bool(top_n == len(ordered)),
            **summary,
        })
topn_summary = pd.DataFrame(topn_rows)

print("Effect concordance summary:")
display(effect_summary)

print("\nP-value concordance summary:")
display(pvalue_summary)

print("\nDiscovery stability summary:")
display(discovery_summary)

print("\nTop-effect stability summary:")
display(topn_summary)

effect_summary.to_csv(TABLE_DIR / "effect_concordance_summary.csv", index=False)
pvalue_summary.to_csv(TABLE_DIR / "pvalue_concordance_summary.csv", index=False)
discovery_summary.to_csv(TABLE_DIR / "discovery_stability_summary.csv", index=False)
topn_summary.to_csv(TABLE_DIR / "effect_stability_by_top_n.csv", index=False)

```

### Text output (stream)

```text
Effect concordance summary:

```

### Text output (text/plain)

```text
                                  selection  n_pairs  pearson_r  spearman_r  \
0                  all target-program pairs    14900   0.986701    0.980407   
1   union discoveries at empirical q <= 0.1      975   0.997663    0.989531   
2              top 10% by max absolute beta     1490   0.995673    0.975688   
3                  all target-program pairs    14900   0.978270    0.972849   
4   union discoveries at empirical q <= 0.1      330   0.996569    0.982504   
5              top 10% by max absolute beta     1490   0.992036    0.952534   
6                  all target-program pairs    14900   0.981104    0.972772   
7   union discoveries at empirical q <= 0.1      595   0.995848    0.979806   
8              top 10% by max absolute beta     1490   0.993342    0.959373   
9                  all target-program pairs    14900   0.989733    0.979115   
10  union discoveries at empirical q <= 0.1     1403   0.997947    0.994264   
11             top 10% by max absolute beta     1490   0.997553    0.985730   

    concordance_correlation_coefficient  regression_slope_6_on_3  \
0                              0.986597                 0.972423   
1                              0.997640                 0.991593   
2                              0.995609                 0.984725   
3                              0.977761                 0.947202   
4                              0.996324                 0.974788   
5                              0.991465                 0.960570   
6                              0.980656                 0.951882   
7                              0.995654                 0.976688   
8                              0.993027                 0.968677   
9                              0.989483                 0.967732   
10                             0.997784                 0.980137   
11                             0.997318                 0.976178   

    regression_intercept_6_on_3  rmse_beta  normalized_rmse_by_sd_3cov  \
0                 -4.888725e-18   0.032966                    0.162542   
1                 -1.269756e-03   0.037355                    0.068461   
2                  1.332349e-03   0.046443                    0.093165   
3                  3.183355e-18   0.047557                    0.207571   
4                  4.524015e-04   0.064130                    0.084689   
5                  5.511779e-03   0.068844                    0.128558   
6                  7.688372e-18   0.046887                    0.193778   
7                 -1.445200e-03   0.070630                    0.092266   
8                  1.077266e-03   0.069588                    0.116595   
9                  2.330671e-18   0.042642                    0.143420   
10                -7.286016e-04   0.051935                    0.065962   
11                -4.996135e-04   0.057508                    0.072430   

    median_abs_beta_change  p95_abs_beta_change  \
0                 0.013039             0.066720   
1                 0.014586             0.072835   
2                 0.018030             0.095398   
3                 0.015546             0.105041   
4                 0.017097             0.141508   
5                 0.020783             0.150411   
6                 0.019533             0.098251   
7                 0.022795             0.153269   
8                 0.026714             0.154740   
9                 0.018212             0.090357   
10                0.024365             0.104323   
11                0.028569             0.119915   

    mean_beta_difference_6_minus_3  bland_altman_lower_95  \
0                     7.758522e-19              -0.064616   
1                    -1.631374e-03              -0.074815   
2                     1.182418e-03              -0.089846   
3                    -1.539596e-18              -0.093214   
4                    -1.411965e-03              -0.127267   
5                     5.604170e-03              -0.128927   
6                     1.050241e-17              -0.091901   
7                    -2.745375e-03              -0.141192   
8                     6.509535e-04              -0.135782   
9                    -7.870289e-19              -0.083581   
10                   -1.144098e-03              -0.102948   
11                   -8.273647e-04              -0.113570   

    bland_altman_upper_95  sign_concordance  median_abs_beta_3cov  \
0                0.064616          0.951409              0.108571   
1                0.071552          1.000000              0.448107   
2                0.092211          0.999329              0.401281   
3                0.093214          0.945906              0.135028   
4                0.124443          1.000000              0.646110   
5                0.140136          1.000000              0.445934   
6                0.091901          0.937114              0.126791   
7                0.135701          1.000000              0.611144   
8                0.137084          1.000000              0.462319   
9                0.083581          0.946510              0.128097   
10               0.100659          1.000000              0.570305   
11               0.111915          1.000000              0.582378   

    median_abs_beta_6cov  median_abs_effect_ratio_6_over_3 day  
0               0.106209                          0.978246  d0  
1               0.446637                          0.996719  d0  
2               0.396186                          0.987303  d0  
3               0.129849                          0.961650  d1  
4               0.631054                          0.976698  d1  
5               0.431153                          0.966854  d1  
6               0.120362                          0.949295  d2  
7               0.589335                          0.964314  d2  
8               0.442564                          0.957270  d2  
9               0.124751                          0.973877  d3  
10              0.563513                          0.988091  d3  
11              0.566994                          0.973584  d3  
```

### Text output (stream)

```text

P-value concordance summary:

```

### Text output (text/plain)

```text
  day pvalue_type  n_pairs  pearson_neglog10_p  spearman_neglog10_p  \
0  d0        skew    14900            0.994638             0.960684   
1  d0   empirical    14900            0.978136             0.943524   
2  d1        skew    14900            0.982882             0.937574   
3  d1   empirical    14900            0.963606             0.924117   
4  d2        skew    14900            0.988217             0.944016   
5  d2   empirical    14900            0.964992             0.920299   
6  d3        skew    14900            0.995685             0.958326   
7  d3   empirical    14900            0.976069             0.936958   

   median_abs_change_neglog10_p  p95_abs_change_neglog10_p  \
0                      0.048336                   0.322253   
1                      0.058524                   0.341106   
2                      0.049900                   0.345445   
3                      0.052318                   0.382203   
4                      0.057206                   0.346061   
5                      0.070247                   0.392084   
6                      0.061884                   0.428857   
7                      0.067219                   0.414973   

   median_change_neglog10_p_6_minus_3  fraction_more_significant_in_6cov  
0                           -0.000704                           0.493960  
1                            0.000000                           0.480067  
2                           -0.001167                           0.491477  
3                           -0.002901                           0.469530  
4                           -0.000335                           0.498121  
5                           -0.004178                           0.468725  
6                           -0.000254                           0.498859  
7                            0.000000                           0.479799  
```

### Text output (stream)

```text

Discovery stability summary:

```

### Text output (text/plain)

```text
   day pvalue_type  q_cutoff  n_discoveries_3cov  n_discoveries_6cov  \
0   d0        skew      0.05                 634                 651   
1   d0        skew      0.10                 911                 904   
2   d0   empirical      0.05                 640                 638   
3   d0   empirical      0.10                 909                 909   
4   d1        skew      0.05                 215                 200   
5   d1        skew      0.10                 285                 281   
6   d1   empirical      0.05                 210                 190   
7   d1   empirical      0.10                 313                 278   
8   d2        skew      0.05                 379                 371   
9   d2        skew      0.10                 521                 505   
10  d2   empirical      0.05                 367                 363   
11  d2   empirical      0.10                 540                 512   
12  d3        skew      0.05                1087                1120   
13  d3        skew      0.10                1462                1500   
14  d3   empirical      0.05                 954                 988   
15  d3   empirical      0.10                1279                1298   

    n_intersection  n_3cov_only  n_6cov_only  n_union   jaccard  \
0              595           39           56      690  0.862319   
1              837           74           67      978  0.855828   
2              586           54           52      692  0.846821   
3              843           66           66      975  0.864615   
4              186           29           14      229  0.812227   
5              257           28           24      309  0.831715   
6              170           40           20      230  0.739130   
7              261           52           17      330  0.790909   
8              336           43           35      414  0.811594   
9              460           61           45      566  0.812721   
10             319           48           44      411  0.776156   
11             457           83           55      595  0.768067   
12            1019           68          101     1188  0.857744   
13            1371           91          129     1591  0.861722   
14             891           63           97     1051  0.847764   
15            1174          105          124     1403  0.836778   

    sign_concordance_intersection  sign_concordance_union  
0                             1.0                     1.0  
1                             1.0                     1.0  
2                             1.0                     1.0  
3                             1.0                     1.0  
4                             1.0                     1.0  
5                             1.0                     1.0  
6                             1.0                     1.0  
7                             1.0                     1.0  
8                             1.0                     1.0  
9                             1.0                     1.0  
10                            1.0                     1.0  
11                            1.0                     1.0  
12                            1.0                     1.0  
13                            1.0                     1.0  
14                            1.0                     1.0  
15                            1.0                     1.0  
```

### Text output (stream)

```text

Top-effect stability summary:

```

### Text output (text/plain)

```text
   day  top_n  is_all_pairs  selection  n_pairs  pearson_r  spearman_r  \
0   d0    100         False    top 100      100   0.998964    0.988263   
1   d0    250         False    top 250      250   0.998363    0.986866   
2   d0    500         False    top 500      500   0.997545    0.979564   
3   d0   1000         False   top 1000     1000   0.996592    0.977508   
4   d0   2000         False   top 2000     2000   0.994921    0.974157   
5   d0   5000         False   top 5000     5000   0.992113    0.971540   
6   d0  14900          True  top 14900    14900   0.986701    0.980407   
7   d1    100         False    top 100      100   0.999121    0.984590   
8   d1    250         False    top 250      250   0.996654    0.972015   
9   d1    500         False    top 500      500   0.995380    0.961323   
10  d1   1000         False   top 1000     1000   0.993040    0.951680   
11  d1   2000         False   top 2000     2000   0.990998    0.951678   
12  d1   5000         False   top 5000     5000   0.986762    0.956202   
13  d1  14900          True  top 14900    14900   0.978270    0.972849   
14  d2    100         False    top 100      100   0.998616    0.981746   
15  d2    250         False    top 250      250   0.997136    0.980350   
16  d2    500         False    top 500      500   0.995392    0.967957   
17  d2   1000         False   top 1000     1000   0.994368    0.960483   
18  d2   2000         False   top 2000     2000   0.992358    0.960618   
19  d2   5000         False   top 5000     5000   0.988285    0.961213   
20  d2  14900          True  top 14900    14900   0.981104    0.972772   
21  d3    100         False    top 100      100   0.999497    0.987063   
22  d3    250         False    top 250      250   0.999350    0.993763   
23  d3    500         False    top 500      500   0.998660    0.989396   
24  d3   1000         False   top 1000     1000   0.998071    0.987871   
25  d3   2000         False   top 2000     2000   0.996830    0.983104   
26  d3   5000         False   top 5000     5000   0.994019    0.973532   
27  d3  14900          True  top 14900    14900   0.989733    0.979115   

    concordance_correlation_coefficient  regression_slope_6_on_3  \
0                              0.998936                 0.991502   
1                              0.998333                 0.990846   
2                              0.997524                 0.991121   
3                              0.996557                 0.988267   
4                              0.994845                 0.983404   
5                              0.991988                 0.977667   
6                              0.986597                 0.972423   
7                              0.999094                 0.995051   
8                              0.996511                 0.980116   
9                              0.995051                 0.971334   
10                             0.992508                 0.963128   
11                             0.990383                 0.958162   
12                             0.986159                 0.953832   
13                             0.977761                 0.947202   
14                             0.998574                 0.994576   
15                             0.997041                 0.986467   
16                             0.995207                 0.977329   
17                             0.994109                 0.971936   
18                             0.991990                 0.965767   
19                             0.987839                 0.959119   
20                             0.980656                 0.951882   
21                             0.999426                 0.987760   
22                             0.999219                 0.983664   
23                             0.998445                 0.978189   
24                             0.997842                 0.976932   
25                             0.996583                 0.974874   
26                             0.993755                 0.971716   
27                             0.989483                 0.967732   

    regression_intercept_6_on_3  ...  normalized_rmse_by_sd_3cov  \
0                  5.572022e-05  ...                    0.045738   
1                  1.938095e-03  ...                    0.057407   
2                  8.488725e-04  ...                    0.070078   
3                  5.725495e-04  ...                    0.082599   
4                  1.921610e-03  ...                    0.100924   
5                  1.926657e-03  ...                    0.125653   
6                 -5.514026e-18  ...                    0.162542   
7                 -6.395880e-03  ...                    0.042276   
8                  3.358854e-03  ...                    0.082677   
9                  6.950737e-03  ...                    0.098203   
10                 6.864893e-03  ...                    0.120528   
11                 4.798995e-03  ...                    0.136381   
12                 2.997728e-03  ...                    0.163611   
13                 6.594094e-18  ...                    0.207571   
14                -9.939245e-03  ...                    0.053043   
15                -7.493796e-03  ...                    0.076372   
16                -3.556168e-03  ...                    0.096930   
17                 8.234050e-04  ...                    0.107276   
18                 1.322668e-03  ...                    0.124854   
19                 9.486218e-04  ...                    0.153657   
20                 1.034591e-17  ...                    0.193778   
21                -4.703110e-03  ...                    0.033509   
22                -4.995308e-03  ...                    0.039135   
23                 1.247605e-03  ...                    0.055150   
24                -7.455571e-04  ...                    0.064977   
25                 4.301237e-04  ...                    0.081743   
26                 2.033197e-03  ...                    0.110501   
27                 2.614899e-18  ...                    0.143420   

    median_abs_beta_change  p95_abs_beta_change  \
0                 0.019167             0.106574   
1                 0.017316             0.088740   
2                 0.018307             0.095865   
3                 0.018337             0.094737   
4                 0.017058             0.094618   
5                 0.015730             0.086086   
6                 0.013039             0.066720   
7                 0.016990             0.103750   
8                 0.019873             0.144902   
9                 0.022279             0.154885   
10                0.020945             0.160642   
11                0.021058             0.147305   
12                0.019863             0.134929   
13                0.015546             0.105041   
14                0.018202             0.172010   
15                0.022344             0.175420   
16                0.027980             0.177239   
17                0.027914             0.158237   
18                0.026472             0.150018   
19                0.025407             0.128908   
20                0.019533             0.098251   
21                0.027856             0.114366   
22                0.030904             0.111062   
23                0.032032             0.129721   
24                0.030424             0.122887   
25                0.028245             0.120561   
26                0.025875             0.117243   
27                0.018212             0.090357   

    mean_beta_difference_6_minus_3  bland_altman_lower_95  \
0                    -6.103965e-04              -0.095068   
1                     1.572774e-03              -0.091766   
2                     4.689684e-04              -0.093726   
3                     2.911106e-04              -0.090782   
4                     1.805173e-03              -0.087860   
5                     1.965633e-03              -0.078985   
6                     8.345301e-19              -0.064616   
7                    -6.511263e-03              -0.095164   
8                     2.120068e-03              -0.136260   
9                     5.665722e-03              -0.131714   
10                    6.903088e-03              -0.133705   
11                    5.008964e-03              -0.126128   
12                    3.023375e-03              -0.114220   
13                   -1.430623e-18              -0.093214   
14                   -1.085003e-02              -0.146590   
15                   -8.888373e-03              -0.161326   
16                   -5.078878e-03              -0.163579   
17                   -4.777643e-05              -0.142227   
18                    1.111314e-03              -0.131750   
19                    9.504707e-04              -0.117423   
20                    1.049124e-17              -0.091901   
21                   -2.720550e-03              -0.122273   
22                   -4.824127e-03              -0.114132   
23                    1.809211e-03              -0.123379   
24                   -7.891788e-04              -0.117666   
25                    6.752029e-05              -0.113879   
26                    2.017802e-03              -0.104674   
27                   -5.960929e-19              -0.083581   

    bland_altman_upper_95  sign_concordance  median_abs_beta_3cov  \
0                0.093848          1.000000              0.918655   
1                0.094912          1.000000              0.702945   
2                0.094664          1.000000              0.566734   
3                0.091364          0.999000              0.459229   
4                0.091471          0.999000              0.360981   
5                0.082917          0.998400              0.241583   
6                0.064616          0.951409              0.108571   
7                0.082141          1.000000              0.959478   
8                0.140500          1.000000              0.727579   
9                0.143045          1.000000              0.601471   
10               0.147511          1.000000              0.502641   
11               0.136146          0.999500              0.408616   
12               0.120267          0.998200              0.288322   
13               0.093214          0.945906              0.135028   
14               0.124890          1.000000              1.150836   
15               0.143550          1.000000              0.854818   
16               0.153421          1.000000              0.674384   
17               0.142131          1.000000              0.526810   
18               0.133973          0.999500              0.417164   
19               0.119324          0.998000              0.288683   
20               0.091901          0.937114              0.126791   
21               0.116832          1.000000              1.555458   
22               0.104484          1.000000              1.201485   
23               0.126998          1.000000              0.931798   
24               0.116088          1.000000              0.699372   
25               0.114014          1.000000              0.496972   
26               0.108710          0.999000              0.307860   
27               0.083581          0.946510              0.128097   

    median_abs_beta_6cov  median_abs_effect_ratio_6_over_3  
0               0.890795                          0.969673  
1               0.691760                          0.984088  
2               0.562941                          0.993309  
3               0.458753                          0.998964  
4               0.355230                          0.984070  
5               0.237836                          0.984490  
6               0.106209                          0.978246  
7               0.947530                          0.987547  
8               0.702339                          0.965310  
9               0.585224                          0.972989  
10              0.484500                          0.963908  
11              0.393160                          0.962174  
12              0.276851                          0.960213  
13              0.129849                          0.961650  
14              1.141628                          0.992000  
15              0.841476                          0.984392  
16              0.658990                          0.977173  
17              0.506944                          0.962289  
18              0.404425                          0.969462  
19              0.276728                          0.958586  
20              0.120362                          0.949295  
21              1.552765                          0.998269  
22              1.166067                          0.970521  
23              0.913002                          0.979829  
24              0.677962                          0.969386  
25              0.485598                          0.977113  
26              0.301484                          0.979288  
27              0.124751                          0.973877  

[28 rows x 21 columns]
```

## Markdown cell 10

## 4. Beta-effect concordance plots

## Code cell 11

```python

all_effect_metrics = effect_summary[
    effect_summary["selection"] == "all target-program pairs"
][["day", "pearson_r", "spearman_r", "concordance_correlation_coefficient", "regression_slope_6_on_3"]].copy()

all_effect_metrics["facet_label"] = all_effect_metrics.apply(
    lambda row: (
        f"{row['day']}\n"
        f"r={format_metric(row['pearson_r'])}; "
        f"rho={format_metric(row['spearman_r'])}\n"
        f"CCC={format_metric(row['concordance_correlation_coefficient'])}; "
        f"slope={format_metric(row['regression_slope_6_on_3'])}"
    ),
    axis=1,
)

beta_plot_df = paired_effects.merge(
    all_effect_metrics[["day", "facet_label"]],
    on="day",
    how="left",
)

background = beta_plot_df[beta_plot_df["discovery_status"] == "Neither"].copy()
foreground = beta_plot_df[beta_plot_df["discovery_status"] != "Neither"].copy()

beta_limit = float(
    np.nanmax(np.abs(beta_plot_df[["beta_3cov", "beta_6cov"]].to_numpy(dtype=float)))
)
beta_limit = beta_limit * 1.04 if beta_limit > 0 else 1.0

p_beta_all = (
    ggplot()
    + geom_abline(intercept=0, slope=1, linetype="dashed", color=DIAGONAL_COLOR, size=0.7)
    + geom_hline(yintercept=0, color=GRID_COLOR, size=0.35)
    + geom_point(
        data=background,
        mapping=aes("beta_3cov", "beta_6cov"),
        color=DISCOVERY_COLORS["Neither"],
        size=0.55,
        alpha=0.28,
        stroke=0,
    )
    + geom_point(
        data=foreground,
        mapping=aes("beta_3cov", "beta_6cov", color="discovery_status"),
        size=1.15,
        alpha=0.82,
        stroke=0,
    )
    + scale_color_manual(
        values=DISCOVERY_COLORS,
        breaks=["Both", "3 covariates only", "6 covariates only"],
        name=f"Discovery status\n{PRIMARY_PVALUE_TYPE_USED} q <= {PRIMARY_Q_CUTOFF}",
    )
    + facet_wrap("~facet_label", ncol=2)
    + coord_equal(xlim=(-beta_limit, beta_limit), ylim=(-beta_limit, beta_limit))
    + labs(
        title="Three- versus six-covariate DE program effects",
        subtitle="All real target-program pairs; colored points are discoveries in either model",
        x="Beta effect: 3-covariate model",
        y="Beta effect: 6-covariate model",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(8.8, 7.4),
        plot_title=element_text(size=13, weight="bold", ha="center"),
        plot_subtitle=element_text(size=9, ha="center"),
        strip_text=element_text(size=8.5),
        legend_position="top",
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
    )
)

display(p_beta_all)
save_plot(
    p_beta_all,
    PLOT_DIR / "beta_concordance_all_pairs_3cov_vs_6cov",
    width=8.8,
    height=7.4,
)

selected_beta_df = beta_plot_df[beta_plot_df["union_discovery"]].copy()
if not selected_beta_df.empty:
    p_beta_selected = (
        ggplot(selected_beta_df, aes("beta_3cov", "beta_6cov", color="discovery_status"))
        + geom_abline(intercept=0, slope=1, linetype="dashed", color=DIAGONAL_COLOR, size=0.7)
        + geom_hline(yintercept=0, color=GRID_COLOR, size=0.35)
        + geom_point(size=1.65, alpha=0.85, stroke=0)
        + scale_color_manual(
            values=DISCOVERY_COLORS,
            breaks=["Both", "3 covariates only", "6 covariates only"],
            name="Discovery status",
        )
        + facet_wrap("~facet_label", ncol=2)
        + coord_equal(xlim=(-beta_limit, beta_limit), ylim=(-beta_limit, beta_limit))
        + labs(
            title="Effect stability among discoveries in either model",
            subtitle=f"Union of {PRIMARY_PVALUE_TYPE_USED} q <= {PRIMARY_Q_CUTOFF} discoveries",
            x="Beta effect: 3-covariate model",
            y="Beta effect: 6-covariate model",
        )
        + theme_bw(base_size=10)
        + theme(
            figure_size=(8.8, 7.4),
            plot_title=element_text(size=13, weight="bold", ha="center"),
            plot_subtitle=element_text(size=9, ha="center"),
            strip_text=element_text(size=8.5),
            legend_position="top",
            panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(color=BORDER_COLOR, size=0.65),
        )
    )
    display(p_beta_selected)
    save_plot(
        p_beta_selected,
        PLOT_DIR / "beta_concordance_union_discoveries_3cov_vs_6cov",
        width=8.8,
        height=7.4,
    )
else:
    print("No union discoveries at the primary q-value cutoff; selected-effect plot skipped.")

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f419e84c190>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f418c74d5d0>
```

## Markdown cell 12

## 5. Bland–Altman effect-difference plot

## Code cell 13

```python

bland_altman_summary = (
    effect_summary[effect_summary["selection"] == "all target-program pairs"]
    [[
        "day",
        "mean_beta_difference_6_minus_3",
        "bland_altman_lower_95",
        "bland_altman_upper_95",
        "median_abs_beta_change",
    ]]
    .copy()
)

bland_altman_summary["facet_label"] = bland_altman_summary.apply(
    lambda row: (
        f"{row['day']}\n"
        f"mean delta={format_metric(row['mean_beta_difference_6_minus_3'], 3)}; "
        f"median |delta|={format_metric(row['median_abs_beta_change'], 3)}"
    ),
    axis=1,
)

ba_df = paired_effects.merge(
    bland_altman_summary[["day", "facet_label"]],
    on="day",
    how="left",
)

ba_background = ba_df[ba_df["discovery_status"] == "Neither"]
ba_foreground = ba_df[ba_df["discovery_status"] != "Neither"]

p_bland_altman = (
    ggplot()
    + geom_hline(yintercept=0, linetype="dashed", color=DIAGONAL_COLOR, size=0.7)
    + geom_point(
        data=ba_background,
        mapping=aes("beta_mean", "beta_difference_6_minus_3"),
        color=DISCOVERY_COLORS["Neither"],
        size=0.55,
        alpha=0.28,
        stroke=0,
    )
    + geom_point(
        data=ba_foreground,
        mapping=aes("beta_mean", "beta_difference_6_minus_3", color="discovery_status"),
        size=1.15,
        alpha=0.82,
        stroke=0,
    )
    + scale_color_manual(
        values=DISCOVERY_COLORS,
        breaks=["Both", "3 covariates only", "6 covariates only"],
        name="Discovery status",
    )
    + facet_wrap("~facet_label", ncol=2, scales="free")
    + labs(
        title="Bland–Altman comparison of beta effects",
        subtitle="Negative values indicate smaller effects in the six-covariate model",
        x="Mean beta across the two models",
        y="Beta difference: six covariates minus three covariates",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(8.8, 7.2),
        plot_title=element_text(size=13, weight="bold", ha="center"),
        plot_subtitle=element_text(size=9, ha="center"),
        strip_text=element_text(size=8.5),
        legend_position="top",
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
    )
)

display(bland_altman_summary)
display(p_bland_altman)
save_plot(
    p_bland_altman,
    PLOT_DIR / "bland_altman_beta_6cov_minus_3cov",
    width=8.8,
    height=7.2,
)

bland_altman_summary.to_csv(TABLE_DIR / "bland_altman_summary.csv", index=False)

```

### Text output (text/plain)

```text
  day  mean_beta_difference_6_minus_3  bland_altman_lower_95  \
0  d0                    7.758522e-19              -0.064616   
3  d1                   -1.539596e-18              -0.093214   
6  d2                    1.050241e-17              -0.091901   
9  d3                   -7.870289e-19              -0.083581   

   bland_altman_upper_95  median_abs_beta_change  \
0               0.064616                0.013039   
3               0.093214                0.015546   
6               0.091901                0.019533   
9               0.083581                0.018212   

                                   facet_label  
0   d0\nmean delta=0.000; median |delta|=0.013  
3  d1\nmean delta=-0.000; median |delta|=0.016  
6   d2\nmean delta=0.000; median |delta|=0.020  
9  d3\nmean delta=-0.000; median |delta|=0.018  
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f418ccdbcd0>
```

## Markdown cell 14

## 6. P-value concordance and discovery-overlap plots

## Code cell 15

```python

primary_p_metrics = pvalue_summary[pvalue_summary["pvalue_type"] == PRIMARY_PVALUE_TYPE_USED].copy()
primary_p_metrics["facet_label"] = primary_p_metrics.apply(
    lambda row: (
        f"{row['day']}\n"
        f"Spearman={format_metric(row['spearman_neglog10_p'])}; "
        f"median |delta|={format_metric(row['median_abs_change_neglog10_p'], 2)}"
    ),
    axis=1,
)

pvalue_plot_df = paired_primary.merge(
    primary_p_metrics[["day", "facet_label"]],
    on="day",
    how="left",
)

pvalue_background = pvalue_plot_df[pvalue_plot_df["discovery_status"] == "Neither"]
pvalue_foreground = pvalue_plot_df[pvalue_plot_df["discovery_status"] != "Neither"]

pvalue_limit = float(
    np.nanmax(
        pvalue_plot_df[["neglog10_p_3cov", "neglog10_p_6cov"]].to_numpy(dtype=float)
    )
)
pvalue_limit = max(1.0, pvalue_limit * 1.03)

p_pvalue = (
    ggplot()
    + geom_abline(intercept=0, slope=1, linetype="dashed", color=DIAGONAL_COLOR, size=0.7)
    + geom_point(
        data=pvalue_background,
        mapping=aes("neglog10_p_3cov", "neglog10_p_6cov"),
        color=DISCOVERY_COLORS["Neither"],
        size=0.55,
        alpha=0.28,
        stroke=0,
    )
    + geom_point(
        data=pvalue_foreground,
        mapping=aes("neglog10_p_3cov", "neglog10_p_6cov", color="discovery_status"),
        size=1.15,
        alpha=0.82,
        stroke=0,
    )
    + scale_color_manual(
        values=DISCOVERY_COLORS,
        breaks=["Both", "3 covariates only", "6 covariates only"],
        name="Discovery status",
    )
    + facet_wrap("~facet_label", ncol=2)
    + coord_equal(xlim=(0, pvalue_limit), ylim=(0, pvalue_limit))
    + labs(
        title=f"Target-program {PRIMARY_PVALUE_TYPE_USED} p-value concordance",
        subtitle="BH q-values were calculated separately within each day and model",
        x="-log10(p): 3-covariate model",
        y="-log10(p): 6-covariate model",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(8.8, 7.4),
        plot_title=element_text(size=13, weight="bold", ha="center"),
        plot_subtitle=element_text(size=9, ha="center"),
        strip_text=element_text(size=8.5),
        legend_position="top",
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
    )
)

display(p_pvalue)
save_plot(
    p_pvalue,
    PLOT_DIR / f"pvalue_concordance_{PRIMARY_PVALUE_TYPE_USED}_3cov_vs_6cov",
    width=8.8,
    height=7.4,
)

primary_discovery = discovery_summary[
    discovery_summary["pvalue_type"] == PRIMARY_PVALUE_TYPE_USED
].copy()
primary_discovery["q_label"] = primary_discovery["q_cutoff"].map(lambda x: f"q <= {x:g}")

discovery_count_plot = primary_discovery.melt(
    id_vars=["day", "q_label"],
    value_vars=["n_intersection", "n_3cov_only", "n_6cov_only"],
    var_name="category_code",
    value_name="n_discoveries",
)
discovery_count_plot["category"] = discovery_count_plot["category_code"].map({
    "n_intersection": "Both",
    "n_3cov_only": "3 covariates only",
    "n_6cov_only": "6 covariates only",
})

p_discovery_counts = (
    ggplot(discovery_count_plot, aes("day", "n_discoveries", fill="category"))
    + geom_col(position="stack", width=0.72)
    + scale_fill_manual(
        values=DISCOVERY_COLORS,
        breaks=["Both", "3 covariates only", "6 covariates only"],
        name="Discovery set",
    )
    + facet_wrap("~q_label", nrow=1, scales="free_y")
    + labs(
        title=f"Discovery overlap using {PRIMARY_PVALUE_TYPE_USED} p-values",
        x="Differentiation day",
        y="Number of target-program discoveries",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(8.2, 4.3),
        plot_title=element_text(size=12.5, weight="bold", ha="center"),
        legend_position="top",
        panel_grid_major_x=element_blank(),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
    )
)

display(p_discovery_counts)
save_plot(
    p_discovery_counts,
    PLOT_DIR / f"discovery_overlap_counts_{PRIMARY_PVALUE_TYPE_USED}",
    width=8.2,
    height=4.3,
)

p_jaccard = (
    ggplot(primary_discovery, aes("q_cutoff", "jaccard", color="day", group="day"))
    + geom_line(size=0.8)
    + geom_point(size=2.2)
    + scale_x_continuous(breaks=Q_CUTOFFS)
    + scale_y_continuous(limits=(0, 1.02), breaks=np.arange(0, 1.01, 0.2))
    + labs(
        title="Discovery-set Jaccard overlap",
        subtitle="1.0 means identical discovery sets",
        x="BH q-value cutoff",
        y="Jaccard overlap",
        color="Day",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(5.8, 4.4),
        plot_title=element_text(size=12.5, weight="bold", ha="center"),
        plot_subtitle=element_text(size=9, ha="center"),
        legend_position="top",
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
    )
)

display(p_jaccard)
save_plot(
    p_jaccard,
    PLOT_DIR / f"discovery_jaccard_{PRIMARY_PVALUE_TYPE_USED}",
    width=5.8,
    height=4.4,
)

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f419d494410>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f419eff73d0>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f41a49e3650>
```

### Text output (text/plain)

```text
{'pdf': PosixPath('/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/results/ipsc_ec_3cov_vs_6cov_comparison/k_50.dt_0_2/plots/discovery_jaccard_empirical.pdf'),
 'svg': PosixPath('/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/results/ipsc_ec_3cov_vs_6cov_comparison/k_50.dt_0_2/plots/discovery_jaccard_empirical.svg'),
 'png': PosixPath('/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/results/ipsc_ec_3cov_vs_6cov_comparison/k_50.dt_0_2/plots/discovery_jaccard_empirical.png')}
```

## Markdown cell 16


## 7. NTC pseudo-target calibration comparison

The NTC metrics assess how close each model's pseudo-target p-values are to Uniform(0,1). Because pseudo-target tests can be dependent, the notebook reports the **KS distance** but not the nominal KS test p-value. Smaller absolute type-I-error deviation, QQ RMSE, tail QQ RMSE, and KS distance indicate better calibration.


## Code cell 17

```python

ntc_tables = []
for model_code, root in MODEL_ROOTS.items():
    for day in DAYS:
        paths = paths_by_model_day[(model_code, day)]
        for pvalue_type in ntc_pvalue_types:
            current = read_ntc_wide(paths[f"ntc_{pvalue_type}"])
            current.insert(0, "pvalue_type", pvalue_type)
            current.insert(0, "day", day)
            current.insert(0, "model_label", MODEL_LABELS[model_code])
            current.insert(0, "model", model_code)
            ntc_tables.append(current)

ntc_long = pd.concat(ntc_tables, ignore_index=True, sort=False)

ntc_metric_rows = []
for (model_code, model_label, day, pvalue_type), group in ntc_long.groupby(
    ["model", "model_label", "day", "pvalue_type"],
    sort=False,
):
    ntc_metric_rows.append({
        "model": model_code,
        "model_label": model_label,
        "day": day,
        "pvalue_type": pvalue_type,
        **uniform_calibration_metrics(group["pvalue"]),
    })
ntc_calibration_summary = pd.DataFrame(ntc_metric_rows)

# Pair the same ensemble/group/program pseudo-target tests across models.
ntc_pair_rows = []
for day in DAYS:
    for pvalue_type in ntc_pvalue_types:
        left = ntc_long[
            (ntc_long["model"] == "3cov")
            & (ntc_long["day"] == day)
            & (ntc_long["pvalue_type"] == pvalue_type)
        ].rename(columns={"pvalue": "pvalue_3cov"})
        right = ntc_long[
            (ntc_long["model"] == "6cov")
            & (ntc_long["day"] == day)
            & (ntc_long["pvalue_type"] == pvalue_type)
        ].rename(columns={"pvalue": "pvalue_6cov"})

        id_columns = [c for c in ["ensemble", "group_id", "program_name"] if c in left.columns and c in right.columns]
        joined = left[id_columns + ["pvalue_3cov"]].merge(
            right[id_columns + ["pvalue_6cov"]],
            on=id_columns,
            how="outer",
            indicator=True,
            validate="one_to_one",
        )
        matched = joined[joined["_merge"] == "both"].copy()
        ntc_pair_rows.append({
            "day": day,
            "pvalue_type": pvalue_type,
            "n_3cov": len(left),
            "n_6cov": len(right),
            "n_matched": len(matched),
            "match_fraction": float((joined["_merge"] == "both").mean()),
            "pearson_neglog10_p": safe_pearson(
                -np.log10(matched["pvalue_3cov"].clip(1e-300, 1.0)),
                -np.log10(matched["pvalue_6cov"].clip(1e-300, 1.0)),
            ),
            "spearman_neglog10_p": safe_spearman(
                -np.log10(matched["pvalue_3cov"].clip(1e-300, 1.0)),
                -np.log10(matched["pvalue_6cov"].clip(1e-300, 1.0)),
            ),
        })
ntc_between_model_summary = pd.DataFrame(ntc_pair_rows)

print("NTC calibration metrics:")
display(ntc_calibration_summary)

print("\nPaired NTC p-value concordance:")
display(ntc_between_model_summary)

# Build QQ coordinates independently for each model/day.
qq_frames = []
for (model_code, model_label, day), group in ntc_long[
    ntc_long["pvalue_type"] == PRIMARY_NTC_PVALUE_TYPE_USED
].groupby(["model", "model_label", "day"], sort=False):
    points = qq_points(group["pvalue"])
    points["model"] = model_code
    points["model_label"] = model_label
    points["day"] = day
    qq_frames.append(points)
ntc_qq = pd.concat(qq_frames, ignore_index=True)

ntc_qq_limit = float(
    np.nanmax(
        ntc_qq[["expected_neglog10_p", "observed_neglog10_p"]].to_numpy(dtype=float)
    )
)
ntc_qq_limit = max(1.0, ntc_qq_limit * 1.03)

p_ntc_qq = (
    ggplot(ntc_qq, aes("expected_neglog10_p", "observed_neglog10_p", color="model_label"))
    + geom_abline(intercept=0, slope=1, linetype="dashed", color=DIAGONAL_COLOR, size=0.75)
    + geom_point(size=0.8, alpha=0.60, stroke=0)
    + scale_color_manual(values=MODEL_COLORS, name="Model")
    + facet_wrap("~day", ncol=2)
    + coord_equal(xlim=(0, ntc_qq_limit), ylim=(0, ntc_qq_limit))
    + labs(
        title=f"NTC pseudo-target calibration: {PRIMARY_NTC_PVALUE_TYPE_USED} p-values",
        x="Expected -log10(p)",
        y="Observed -log10(p)",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(8.2, 7.0),
        plot_title=element_text(size=13, weight="bold", ha="center"),
        legend_position="top",
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
    )
)

display(p_ntc_qq)
save_plot(
    p_ntc_qq,
    PLOT_DIR / f"NTC_QQ_3cov_vs_6cov_{PRIMARY_NTC_PVALUE_TYPE_USED}",
    width=8.2,
    height=7.0,
)

ntc_primary_metrics = ntc_calibration_summary[
    ntc_calibration_summary["pvalue_type"] == PRIMARY_NTC_PVALUE_TYPE_USED
].copy()

p_ntc_type1 = (
    ggplot(ntc_primary_metrics, aes("day", "fraction_p_lt_0_05", color="model_label", group="model_label"))
    + geom_hline(yintercept=0.05, linetype="dashed", color=DIAGONAL_COLOR, size=0.75)
    + geom_line(size=0.8)
    + geom_point(size=2.3)
    + scale_color_manual(values=MODEL_COLORS, name="Model")
    + labs(
        title="NTC nominal type-I-error rate",
        subtitle="Dashed line is the expected 0.05 fraction",
        x="Differentiation day",
        y="Fraction of NTC p-values < 0.05",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(6.2, 4.4),
        plot_title=element_text(size=12.5, weight="bold", ha="center"),
        plot_subtitle=element_text(size=9, ha="center"),
        legend_position="top",
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
    )
)

display(p_ntc_type1)
save_plot(
    p_ntc_type1,
    PLOT_DIR / f"NTC_fraction_p_lt_0_05_{PRIMARY_NTC_PVALUE_TYPE_USED}",
    width=6.2,
    height=4.4,
)

ntc_long.to_csv(TABLE_DIR / "NTC_pvalues_long.csv", index=False)
ntc_calibration_summary.to_csv(TABLE_DIR / "NTC_calibration_summary.csv", index=False)
ntc_between_model_summary.to_csv(TABLE_DIR / "NTC_between_model_concordance.csv", index=False)

```

### Text output (stream)

```text
NTC calibration metrics:

```

### Text output (text/plain)

```text
   model   model_label day pvalue_type  n_pvalues  median_p    mean_p  \
0   3cov  3 covariates  d0        skew       6100  0.481423  0.485957   
1   3cov  3 covariates  d0   empirical       6100  0.488101  0.490848   
2   3cov  3 covariates  d1        skew       4650  0.493024  0.489781   
3   3cov  3 covariates  d1   empirical       4650  0.491201  0.491911   
4   3cov  3 covariates  d2        skew       6350  0.497277  0.498032   
5   3cov  3 covariates  d2   empirical       6350  0.495350  0.499482   
6   3cov  3 covariates  d3        skew       7750  0.454629  0.468736   
7   3cov  3 covariates  d3   empirical       7750  0.504050  0.501486   
8   6cov  6 covariates  d0        skew       6100  0.479589  0.484979   
9   6cov  6 covariates  d0   empirical       6100  0.487051  0.491186   
10  6cov  6 covariates  d1        skew       4650  0.482028  0.488016   
11  6cov  6 covariates  d1   empirical       4650  0.485551  0.491387   
12  6cov  6 covariates  d2        skew       6350  0.492329  0.496010   
13  6cov  6 covariates  d2   empirical       6350  0.498850  0.498082   
14  6cov  6 covariates  d3        skew       7750  0.449628  0.465894   
15  6cov  6 covariates  d3   empirical       7750  0.489351  0.493662   

    fraction_p_lt_0_01  fraction_p_lt_0_05  fraction_p_lt_0_10  \
0             0.014426            0.066557            0.120984   
1             0.012623            0.058689            0.115902   
2             0.010538            0.055484            0.110753   
3             0.012043            0.057634            0.110753   
4             0.008504            0.047874            0.098425   
5             0.010079            0.050551            0.102047   
6             0.024774            0.083097            0.137161   
7             0.011613            0.051355            0.099355   
8             0.013934            0.065574            0.125738   
9             0.012295            0.059672            0.116066   
10            0.009892            0.058280            0.111398   
11            0.011613            0.057419            0.109247   
12            0.010866            0.052913            0.099685   
13            0.010079            0.053543            0.100787   
14            0.027226            0.084645            0.141548   
15            0.013548            0.055097            0.106323   

    abs_error_at_0_05  ks_distance   qq_rmse  qq_tail_rmse  qq_max_abs_error  
0            0.016557     0.028772  0.051865      0.131762          0.614255  
1            0.008689     0.019583  0.034047      0.103845          0.390805  
2            0.005484     0.021607  0.027877      0.097651          0.413581  
3            0.007634     0.014934  0.025003      0.066921          0.253409  
4            0.002126     0.012381  0.015109      0.117914          0.301757  
5            0.000551     0.007165  0.008298      0.071362          0.294180  
6            0.033097     0.054644  0.146971      0.874913          1.868127  
7            0.001355     0.007520  0.015856      0.148521          0.286833  
8            0.015574     0.032366  0.056529      0.200015          0.947160  
9            0.009672     0.020224  0.033298      0.111614          0.457752  
10           0.008280     0.020919  0.029011      0.115838          0.573394  
11           0.007419     0.016684  0.024103      0.111614          0.429500  
12           0.002913     0.011779  0.009793      0.078252          0.410325  
13           0.003543     0.007004  0.010384      0.085282          0.404790  
14           0.034645     0.056991  0.162435      0.996282          2.411291  
15           0.005097     0.016292  0.032975      0.258613          0.550074  
```

### Text output (stream)

```text

Paired NTC p-value concordance:

```

### Text output (text/plain)

```text
  day pvalue_type  n_3cov  n_6cov  n_matched  match_fraction  \
0  d0        skew    6100    6100       6100             1.0   
1  d0   empirical    6100    6100       6100             1.0   
2  d1        skew    4650    4650       4650             1.0   
3  d1   empirical    4650    4650       4650             1.0   
4  d2        skew    6350    6350       6350             1.0   
5  d2   empirical    6350    6350       6350             1.0   
6  d3        skew    7750    7750       7750             1.0   
7  d3   empirical    7750    7750       7750             1.0   

   pearson_neglog10_p  spearman_neglog10_p  
0            0.971321             0.953820  
1            0.948967             0.910579  
2            0.960747             0.936049  
3            0.951259             0.919436  
4            0.973643             0.954965  
5            0.963174             0.934492  
6            0.978109             0.958028  
7            0.958961             0.932483  
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f41a50ef250>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f41a4e83190>
```

## Markdown cell 18

## 8. Program- and target-level heterogeneity

## Code cell 19

```python

program_rows = []
for (day, program_name), group in paired_effects.groupby(["day", "program_name"], sort=False):
    program_rows.append({
        "day": day,
        "program_name": str(program_name),
        "n_targets": len(group),
        "pearson_r": safe_pearson(group["beta_3cov"], group["beta_6cov"]),
        "spearman_r": safe_spearman(group["beta_3cov"], group["beta_6cov"]),
        "ccc": concordance_correlation_coefficient(group["beta_3cov"], group["beta_6cov"]),
        "median_abs_beta_change": float(np.median(np.abs(group["beta_difference_6_minus_3"]))),
        "median_abs_beta_3cov": float(np.median(group["abs_beta_3cov"])),
        "median_abs_beta_6cov": float(np.median(group["abs_beta_6cov"])),
        "sign_concordance": float(group["sign_concordant"].mean()),
    })
program_summary = pd.DataFrame(program_rows)

target_rows = []
for (day, target_name), group in paired_effects.groupby(["day", "target_name"], sort=False):
    target_rows.append({
        "day": day,
        "target_name": str(target_name),
        "n_programs": len(group),
        "pearson_r": safe_pearson(group["beta_3cov"], group["beta_6cov"]),
        "spearman_r": safe_spearman(group["beta_3cov"], group["beta_6cov"]),
        "ccc": concordance_correlation_coefficient(group["beta_3cov"], group["beta_6cov"]),
        "median_abs_beta_change": float(np.median(np.abs(group["beta_difference_6_minus_3"]))),
        "median_abs_beta_3cov": float(np.median(group["abs_beta_3cov"])),
        "median_abs_beta_6cov": float(np.median(group["abs_beta_6cov"])),
        "sign_concordance": float(group["sign_concordant"].mean()),
    })
target_summary = pd.DataFrame(target_rows)

program_order = sorted(program_summary["program_name"].unique(), key=natural_sort_key)
program_summary["program_factor"] = pd.Categorical(
    program_summary["program_name"],
    categories=program_order,
    ordered=True,
)

p_program_heatmap = (
    ggplot(program_summary, aes("day", "program_factor", fill="pearson_r"))
    + geom_tile()
    + scale_fill_gradient(
        low="#F7FBFF",
        high="#08306B",
        limits=(0, 1),
        name="Pearson r",
    )
    + labs(
        title="Per-program beta concordance",
        subtitle="Each correlation is calculated across real perturbation targets",
        x="Differentiation day",
        y="cNMF program",
    )
    + theme_bw(base_size=9)
    + theme(
        figure_size=(5.8, 11.0),
        plot_title=element_text(size=12.5, weight="bold", ha="center"),
        plot_subtitle=element_text(size=8.5, ha="center"),
        axis_text_y=element_text(size=6.5),
        panel_grid=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
    )
)

display(p_program_heatmap)
save_plot(
    p_program_heatmap,
    PLOT_DIR / "program_level_beta_pearson_heatmap",
    width=5.8,
    height=11.0,
)

p_target_distribution = (
    ggplot(target_summary, aes("day", "pearson_r"))
    + geom_boxplot(outlier_alpha=0, width=0.60)
    + geom_jitter(width=0.16, height=0, size=0.65, alpha=0.30, stroke=0)
    + geom_hline(yintercept=0, linetype="dashed", color=DIAGONAL_COLOR, size=0.6)
    + labs(
        title="Per-target beta concordance across programs",
        x="Differentiation day",
        y="Pearson correlation across 50 programs",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(6.2, 4.7),
        plot_title=element_text(size=12.5, weight="bold", ha="center"),
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
    )
)

display(p_target_distribution)
save_plot(
    p_target_distribution,
    PLOT_DIR / "target_level_beta_correlation_distribution",
    width=6.2,
    height=4.7,
)

worst_programs = (
    program_summary
    .sort_values(["day", "pearson_r"], ascending=[True, True])
    .groupby("day", sort=False)
    .head(10)
)
worst_targets = (
    target_summary
    .sort_values(["day", "pearson_r"], ascending=[True, True])
    .groupby("day", sort=False)
    .head(20)
)

print("Worst 10 programs per day by beta Pearson correlation:")
display(worst_programs)

print("\nWorst 20 targets per day by beta Pearson correlation:")
display(worst_targets)

program_summary.to_csv(TABLE_DIR / "program_level_effect_concordance.csv", index=False)
target_summary.to_csv(TABLE_DIR / "target_level_effect_concordance.csv", index=False)
worst_programs.to_csv(TABLE_DIR / "worst_programs_by_effect_concordance.csv", index=False)
worst_targets.to_csv(TABLE_DIR / "worst_targets_by_effect_concordance.csv", index=False)

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f419d2aa590>
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f41a4303e10>
```

### Text output (stream)

```text
Worst 10 programs per day by beta Pearson correlation:

```

### Text output (text/plain)

```text
    day program_name  n_targets  pearson_r  spearman_r       ccc  \
33   d0            4        298   0.885130    0.888732  0.765557   
37   d0           43        298   0.938135    0.929020  0.932756   
48   d0            8        298   0.959877    0.927105  0.941772   
46   d0            6        298   0.967870    0.953882  0.938498   
10   d0           19        298   0.968745    0.961429  0.959920   
32   d0           39        298   0.971243    0.966202  0.932186   
1    d0           10        298   0.973391    0.969991  0.973097   
21   d0           29        298   0.982988    0.979725  0.981854   
42   d0           48        298   0.986295    0.981690  0.980476   
23   d0           30        298   0.986452    0.983431  0.978507   
83   d1            4        298   0.889561    0.881208  0.759822   
66   d1           24        298   0.925746    0.910108  0.851485   
63   d1           21        298   0.926282    0.913422  0.886843   
71   d1           29        298   0.951959    0.946391  0.933875   
87   d1           43        298   0.962253    0.958933  0.947398   
51   d1           10        298   0.962991    0.958344  0.958121   
72   d1            3        298   0.966380    0.951804  0.935064   
98   d1            8        298   0.967588    0.956988  0.948834   
84   d1           40        298   0.968647    0.962294  0.937458   
85   d1           41        298   0.971314    0.972786  0.968638   
133  d2            4        298   0.913099    0.897149  0.803515   
137  d2           43        298   0.930233    0.917738  0.901505   
123  d2           30        298   0.946969    0.943301  0.932623   
101  d2           10        298   0.951994    0.948244  0.934499   
112  d2           20        298   0.961763    0.959746  0.944137   
122  d2            3        298   0.962989    0.950760  0.919248   
147  d2            7        298   0.976111    0.963657  0.949327   
135  d2           41        298   0.976157    0.968151  0.956035   
113  d2           21        298   0.979348    0.976519  0.962617   
110  d2           19        298   0.979964    0.978631  0.967718   
187  d3           43        298   0.922047    0.914819  0.914785   
162  d3           20        298   0.942046    0.939600  0.903898   
183  d3            4        298   0.946554    0.932186  0.859599   
185  d3           41        298   0.960865    0.960651  0.959563   
173  d3           30        298   0.965871    0.965737  0.964749   
160  d3           19        298   0.970650    0.967573  0.967018   
198  d3            8        298   0.976860    0.976339  0.976391   
171  d3           29        298   0.977475    0.972622  0.967542   
151  d3           10        298   0.982834    0.971800  0.982598   
194  d3            5        298   0.984841    0.978252  0.944357   

     median_abs_beta_change  median_abs_beta_3cov  median_abs_beta_6cov  \
33                 0.095897              0.153744              0.095388   
37                 0.028363              0.089048              0.085683   
48                 0.041915              0.096868              0.083893   
46                 0.029582              0.067302              0.057868   
10                 0.035990              0.123363              0.115287   
32                 0.054039              0.106164              0.125840   
1                  0.023962              0.100107              0.097286   
21                 0.021449              0.114264              0.112996   
42                 0.017687              0.080327              0.082807   
23                 0.028141              0.119962              0.112570   
83                 0.122644              0.196918              0.140759   
66                 0.071238              0.177718              0.133014   
63                 0.064297              0.147557              0.122154   
71                 0.048820              0.130970              0.136806   
87                 0.035116              0.108760              0.110428   
51                 0.033480              0.130514              0.127329   
72                 0.070196              0.210164              0.163838   
98                 0.062939              0.164769              0.180805   
84                 0.055303              0.129869              0.121599   
85                 0.032385              0.135620              0.128636   
133                0.118640              0.195804              0.146745   
137                0.027480              0.063759              0.061196   
123                0.043478              0.131604              0.117402   
101                0.050280              0.158496              0.136039   
112                0.043835              0.148605              0.127004   
122                0.071274              0.190076              0.139027   
147                0.066314              0.184249              0.150295   
135                0.031176              0.096788              0.089817   
113                0.041497              0.156993              0.142924   
110                0.028166              0.111620              0.107130   
187                0.030578              0.084076              0.078748   
162                0.070299              0.148673              0.121459   
183                0.109384              0.174972              0.133886   
185                0.024654              0.104975              0.098216   
173                0.027763              0.107076              0.107183   
160                0.024769              0.092251              0.092108   
198                0.023321              0.117122              0.115721   
171                0.034318              0.131269              0.129331   
151                0.040024              0.200399              0.188987   
194                0.067580              0.146751              0.134909   

     sign_concordance program_factor  
33           0.788591              4  
37           0.882550             43  
48           0.852349              8  
46           0.879195              6  
10           0.906040             19  
32           0.855705             39  
1            0.942953             10  
21           0.932886             29  
42           0.926174             48  
23           0.946309             30  
83           0.808725              4  
66           0.875839             24  
63           0.838926             21  
71           0.902685             29  
87           0.895973             43  
51           0.926174             10  
72           0.892617              3  
98           0.885906              8  
84           0.872483             40  
85           0.909396             41  
133          0.835570              4  
137          0.848993             43  
123          0.899329             30  
101          0.889262             10  
112          0.922819             20  
122          0.875839              3  
147          0.869128              7  
135          0.895973             41  
113          0.919463             21  
110          0.932886             19  
187          0.882550             43  
162          0.862416             20  
183          0.808725              4  
185          0.932886             41  
173          0.939597             30  
160          0.922819             19  
198          0.922819              8  
171          0.942953             29  
151          0.932886             10  
194          0.855705              5  
```

### Text output (stream)

```text

Worst 20 targets per day by beta Pearson correlation:

```

### Text output (text/plain)

```text
     day target_name  n_programs  pearson_r  spearman_r       ccc  \
41    d0      CDKN1A          50   0.905701    0.886291  0.894642   
271   d0       TMTC3          50   0.910639    0.898776  0.910459   
283   d0        VWA9          50   0.916379    0.864970  0.916130   
54    d0       DACH1          50   0.922493    0.912989  0.921949   
201   d0        RGP1          50   0.930190    0.946603  0.924131   
...   ..         ...         ...        ...         ...       ...   
1052  d3       PAM16          50   0.949126    0.952365  0.943840   
956   d3        DPM2          50   0.951858    0.942377  0.951694   
924   d3       CAPN5          50   0.952119    0.936807  0.952107   
916   d3        BRD3          50   0.952133    0.950636  0.947952   
1067  d3      POLR3K          50   0.953632    0.932965  0.948885   

      median_abs_beta_change  median_abs_beta_3cov  median_abs_beta_6cov  \
41                  0.021530              0.069079              0.054390   
271                 0.027778              0.105397              0.110070   
283                 0.043972              0.165201              0.161815   
54                  0.019735              0.071352              0.071738   
201                 0.029238              0.116972              0.095833   
...                      ...                   ...                   ...   
1052                0.029040              0.125843              0.106799   
956                 0.029707              0.093090              0.085632   
924                 0.024591              0.091253              0.084388   
916                 0.028349              0.084055              0.097658   
1067                0.041956              0.152323              0.132735   

      sign_concordance  
41                0.82  
271               0.90  
283               0.92  
54                0.92  
201               0.94  
...                ...  
1052              0.96  
956               0.96  
924               0.92  
916               0.90  
1067              0.98  

[80 rows x 10 columns]
```

## Markdown cell 20

## 9. Effect-stability curves for increasingly inclusive effect sets

## Code cell 21

```python

topn_plot = topn_summary.melt(
    id_vars=["day", "top_n", "is_all_pairs"],
    value_vars=["pearson_r", "spearman_r", "sign_concordance"],
    var_name="metric_code",
    value_name="value",
)
topn_plot["metric"] = topn_plot["metric_code"].map({
    "pearson_r": "Pearson r",
    "spearman_r": "Spearman",
    "sign_concordance": "Sign concordance",
})

p_topn = (
    ggplot(topn_plot, aes("top_n", "value", color="metric", group="metric"))
    + geom_line(size=0.8)
    + geom_point(size=1.8)
    + scale_x_log10()
    + facet_wrap("~day", ncol=2)
    + labs(
        title="Effect agreement as progressively weaker effects are included",
        subtitle="Effects are ranked symmetrically by max(|beta 3cov|, |beta 6cov|)",
        x="Number of strongest target-program effects included (log scale)",
        y="Agreement metric",
        color="Metric",
    )
    + theme_bw(base_size=10)
    + theme(
        figure_size=(8.2, 6.5),
        plot_title=element_text(size=12.5, weight="bold", ha="center"),
        plot_subtitle=element_text(size=8.5, ha="center"),
        legend_position="top",
        panel_grid_major=element_line(color=GRID_COLOR, size=0.4),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color=BORDER_COLOR, size=0.65),
    )
)

display(p_topn)
save_plot(
    p_topn,
    PLOT_DIR / "effect_agreement_by_top_n",
    width=8.2,
    height=6.5,
)

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f4195263910>
```

### Text output (text/plain)

```text
{'pdf': PosixPath('/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/results/ipsc_ec_3cov_vs_6cov_comparison/k_50.dt_0_2/plots/effect_agreement_by_top_n.pdf'),
 'svg': PosixPath('/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/results/ipsc_ec_3cov_vs_6cov_comparison/k_50.dt_0_2/plots/effect_agreement_by_top_n.svg'),
 'png': PosixPath('/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/results/ipsc_ec_3cov_vs_6cov_comparison/k_50.dt_0_2/plots/effect_agreement_by_top_n.png')}
```

## Markdown cell 22

## 10. Optional biological-replicate reproducibility comparison

## Code cell 23

```python

def find_replicate_summary(root):
    candidates = [
        Path(root) / "biological_replicate_effects" / "plots" / "biological_replicate_DE_program_effect_correlation_summary.csv",
        Path(root) / "biological_replicate_effects" / "plots" / "DIAGNOSTIC_recomputed_summary_from_current_beta_files.csv",
    ]
    return next((path for path in candidates if path.exists()), None)

replicate_paths = {
    model: find_replicate_summary(root)
    for model, root in MODEL_ROOTS.items()
}

print("Replicate summary paths:")
for model, path in replicate_paths.items():
    print(f"  {MODEL_LABELS[model]}: {path}")

if all(path is not None for path in replicate_paths.values()):
    replicate_tables = []
    for model, path in replicate_paths.items():
        current = pd.read_csv(path)
        current["model"] = model
        current["model_label"] = MODEL_LABELS[model]
        replicate_tables.append(current)

    replicate_comparison = pd.concat(replicate_tables, ignore_index=True, sort=False)
    required = {"day", "pearson_r", "spearman_r"}
    if required.issubset(replicate_comparison.columns):
        replicate_long = replicate_comparison.melt(
            id_vars=["model", "model_label", "day"],
            value_vars=["pearson_r", "spearman_r"],
            var_name="metric_code",
            value_name="correlation",
        )
        replicate_long["metric"] = replicate_long["metric_code"].map({
            "pearson_r": "Pearson r",
            "spearman_r": "Spearman",
        })

        p_replicate = (
            ggplot(replicate_long, aes("day", "correlation", fill="model_label"))
            + geom_col(position=position_dodge(width=0.78), width=0.70)
            + scale_fill_manual(values=MODEL_COLORS, name="Model")
            + facet_wrap("~metric", nrow=1)
            + labs(
                title="Biological-replicate effect correlation by covariate model",
                x="Differentiation day",
                y="Replicate correlation",
            )
            + theme_bw(base_size=10)
            + theme(
                figure_size=(8.0, 4.3),
                plot_title=element_text(size=12.5, weight="bold", ha="center"),
                legend_position="top",
                panel_grid_major_x=element_blank(),
                panel_grid_minor=element_blank(),
                panel_border=element_rect(color=BORDER_COLOR, size=0.65),
            )
        )
        display(replicate_comparison)
        display(p_replicate)
        save_plot(
            p_replicate,
            PLOT_DIR / "biological_replicate_correlation_3cov_vs_6cov",
            width=8.0,
            height=4.3,
        )
        replicate_comparison.to_csv(
            TABLE_DIR / "biological_replicate_correlation_comparison.csv",
            index=False,
        )
    else:
        print("Replicate summary files were found but do not contain the expected correlation columns.")
else:
    print(
        "Replicate comparison skipped. It will run automatically after replicate-specific "
        "summary files exist under both model roots."
    )

```

### Text output (stream)

```text
Replicate summary paths:
  3 covariates: None
  6 covariates: None
Replicate comparison skipped. It will run automatically after replicate-specific summary files exist under both model roots.

```

## Markdown cell 24

## 11. Manuscript-facing summary and changed-result tables

## Code cell 25

```python

# Largest beta changes and model-specific discoveries at the primary threshold.
largest_effect_changes = (
    paired_primary
    .assign(abs_beta_change=lambda x: (x["beta_6cov"] - x["beta_3cov"]).abs())
    .sort_values(["day", "abs_beta_change"], ascending=[True, False])
    .groupby("day", sort=False)
    .head(100)
)

three_only_discoveries = paired_primary[
    paired_primary["discovery_status"] == "3 covariates only"
].sort_values(["day", "qvalue_3cov"])

six_only_discoveries = paired_primary[
    paired_primary["discovery_status"] == "6 covariates only"
].sort_values(["day", "qvalue_6cov"])

both_discoveries = paired_primary[
    paired_primary["discovery_status"] == "Both"
].sort_values(["day", "qvalue_6cov", "qvalue_3cov"])

# One row per day summarizing the central tradeoff: calibration versus result stability.
all_effect_day = effect_summary[
    effect_summary["selection"] == "all target-program pairs"
].copy()

primary_discovery_day = discovery_summary[
    (discovery_summary["pvalue_type"] == PRIMARY_PVALUE_TYPE_USED)
    & np.isclose(discovery_summary["q_cutoff"], PRIMARY_Q_CUTOFF)
].copy()

ntc_primary_wide = (
    ntc_calibration_summary[
        ntc_calibration_summary["pvalue_type"] == PRIMARY_NTC_PVALUE_TYPE_USED
    ]
    .pivot(
        index="day",
        columns="model",
        values=[
            "fraction_p_lt_0_05",
            "abs_error_at_0_05",
            "qq_rmse",
            "qq_tail_rmse",
            "ks_distance",
        ],
    )
)
ntc_primary_wide.columns = [f"{metric}_{model}" for metric, model in ntc_primary_wide.columns]
ntc_primary_wide = ntc_primary_wide.reset_index()

manuscript_summary = (
    all_effect_day.merge(
        primary_discovery_day,
        on="day",
        how="left",
        suffixes=("_effect", "_discovery"),
    )
    .merge(ntc_primary_wide, on="day", how="left")
)

selected_columns = [
    "day",
    "n_pairs",
    "pearson_r",
    "spearman_r",
    "concordance_correlation_coefficient",
    "regression_slope_6_on_3",
    "median_abs_beta_change",
    "median_abs_effect_ratio_6_over_3",
    "sign_concordance",
    "n_discoveries_3cov",
    "n_discoveries_6cov",
    "n_intersection",
    "n_3cov_only",
    "n_6cov_only",
    "jaccard",
    "sign_concordance_union",
    "fraction_p_lt_0_05_3cov",
    "fraction_p_lt_0_05_6cov",
    "abs_error_at_0_05_3cov",
    "abs_error_at_0_05_6cov",
    "qq_rmse_3cov",
    "qq_rmse_6cov",
    "qq_tail_rmse_3cov",
    "qq_tail_rmse_6cov",
    "ks_distance_3cov",
    "ks_distance_6cov",
]
selected_columns = [c for c in selected_columns if c in manuscript_summary.columns]
manuscript_summary = manuscript_summary[selected_columns]

print("Manuscript-facing model comparison summary:")
display(manuscript_summary)

print("\nLargest beta changes, first 20 rows:")
display(largest_effect_changes.head(20))

print("\nThree-covariate-only discoveries:", len(three_only_discoveries))
display(three_only_discoveries.head(20))

print("\nSix-covariate-only discoveries:", len(six_only_discoveries))
display(six_only_discoveries.head(20))

manuscript_summary.to_csv(TABLE_DIR / "manuscript_model_comparison_summary.csv", index=False)
largest_effect_changes.to_csv(TABLE_DIR / "largest_beta_changes_top100_per_day.csv", index=False)
three_only_discoveries.to_csv(TABLE_DIR / "discoveries_3cov_only.csv", index=False)
six_only_discoveries.to_csv(TABLE_DIR / "discoveries_6cov_only.csv", index=False)
both_discoveries.to_csv(TABLE_DIR / "discoveries_in_both_models.csv", index=False)

print("\nInterpretation guide")
print("--------------------")
print("1. Lower NTC abs_error_at_0_05, QQ RMSE, tail QQ RMSE, and KS distance indicate better calibration.")
print("2. Effect Pearson/Spearman quantify rank and linear stability; CCC and slope additionally detect scale attenuation or shift.")
print("3. A six/three median absolute effect ratio below 1 and a regression slope below 1 indicate overall attenuation in the six-covariate model.")
print("4. Discovery Jaccard and model-specific discovery counts quantify inferential stability.")
print("5. Inspect model-specific discoveries and the Bland-Altman plot rather than relying on one global correlation.")
print("6. This comparison identifies the aggregate impact of adding biological_sample, guide_umi_counts, and n_genes_by_counts; it does not isolate each term.")

print("\nOutput directories:")
print("  Tables:", TABLE_DIR)
print("  Plots: ", PLOT_DIR)
print("\nPDF and SVG files are the preferred Illustrator-editable outputs.")

```

### Text output (stream)

```text
Manuscript-facing model comparison summary:

```

### Text output (text/plain)

```text
  day  n_pairs  pearson_r  spearman_r  concordance_correlation_coefficient  \
0  d0    14900   0.986701    0.980407                             0.986597   
1  d1    14900   0.978270    0.972849                             0.977761   
2  d2    14900   0.981104    0.972772                             0.980656   
3  d3    14900   0.989733    0.979115                             0.989483   

   regression_slope_6_on_3  median_abs_beta_change  \
0                 0.972423                0.013039   
1                 0.947202                0.015546   
2                 0.951882                0.019533   
3                 0.967732                0.018212   

   median_abs_effect_ratio_6_over_3  sign_concordance  n_discoveries_3cov  \
0                          0.978246          0.951409                 909   
1                          0.961650          0.945906                 313   
2                          0.949295          0.937114                 540   
3                          0.973877          0.946510                1279   

   ...  fraction_p_lt_0_05_3cov  fraction_p_lt_0_05_6cov  \
0  ...                 0.058689                 0.059672   
1  ...                 0.057634                 0.057419   
2  ...                 0.050551                 0.053543   
3  ...                 0.051355                 0.055097   

   abs_error_at_0_05_3cov  abs_error_at_0_05_6cov  qq_rmse_3cov  qq_rmse_6cov  \
0                0.008689                0.009672      0.034047      0.033298   
1                0.007634                0.007419      0.025003      0.024103   
2                0.000551                0.003543      0.008298      0.010384   
3                0.001355                0.005097      0.015856      0.032975   

   qq_tail_rmse_3cov  qq_tail_rmse_6cov  ks_distance_3cov  ks_distance_6cov  
0           0.103845           0.111614          0.019583          0.020224  
1           0.066921           0.111614          0.014934          0.016684  
2           0.071362           0.085282          0.007165          0.007004  
3           0.148521           0.258613          0.007520          0.016292  

[4 rows x 26 columns]
```

### Text output (stream)

```text

Largest beta changes, first 20 rows:

```

### Text output (text/plain)

```text
      day pvalue_type target_name program_name  beta_3cov  pvalue_3cov  \
29083  d0   empirical        VWA9            4  -0.392296     0.015798   
29533  d0   empirical        XRN1            4  -0.625779     0.000700   
19083  d0   empirical      EXOSC8            4  -0.025502     0.896910   
16283  d0   empirical     C7orf26            4   0.025610     0.892111   
26633  d0   empirical     SNRNP27            4  -0.290747     0.150685   
23333  d0   empirical       PLRG1            4   0.275451     0.290971   
22083  d0   empirical      MRPL32            4  -0.307576     0.123588   
21533  d0   empirical        MAU2            4  -0.420576     0.009599   
24983  d0   empirical        RGP1            4  -0.227842     0.177682   
28883  d0   empirical        UPF2            4  -0.474643     0.011499   
25833  d0   empirical        SKA3            4  -0.320328     0.205679   
17883  d0   empirical     DNAJC17            4  -0.088329     0.608739   
16983  d0   empirical      CDKN1A            4  -0.200134     0.119288   
20533  d0   empirical      INTS10            4  -0.246034     0.078792   
18383  d0   empirical       EIF3M            4   0.354058     0.086291   
28483  d0   empirical       TMTC3            4  -0.092812     0.529747   
25283  d0   empirical       RNF40            4  -0.149772     0.406059   
25483  d0   empirical       SCAF4            4  -0.271584     0.108889   
19633  d0   empirical         GPI            4  -0.076609     0.648535   
17133  d0   empirical     CHORDC1            4  -0.586456     0.000100   

       qvalue_3cov  beta_6cov  pvalue_6cov  qvalue_6cov  neglog10_p_3cov  \
29083     0.184770   0.009678     0.939406     0.986862         1.801386   
29533     0.022237  -0.269018     0.056494     0.370659         3.154945   
19083     0.980349   0.283401     0.083092     0.449715         0.047251   
16283     0.979186   0.314883     0.009499     0.135053         0.049581   
26633     0.572261  -0.002676     0.986301     0.996737         0.821930   
23333     0.717317   0.548078     0.007799     0.118218         0.536150   
22083     0.532797  -0.037526     0.779822     0.952148         0.908025   
21533     0.134930  -0.178912     0.121788     0.531373         2.017772   
24983     0.608054   0.009264     0.949805     0.989348         0.750356   
28883     0.151622  -0.241361     0.094791     0.476833         1.939346   
25833     0.638863  -0.104691     0.563644     0.881867         0.686809   
17883     0.898575   0.126503     0.399460     0.797849         0.215569   
16983     0.525907   0.013750     0.887911     0.976159         0.923403   
20533     0.434816  -0.033156     0.758724     0.947690         1.103517   
18383     0.455614   0.141892     0.379462     0.784404         1.064033   
28483     0.867197   0.118128     0.359964     0.771809         0.275931   
25283     0.797007   0.061135     0.638036     0.910821         0.391410   
25483     0.507173  -0.060911     0.635636     0.910170         0.963016   
19633     0.913098   0.133134     0.222378     0.656310         0.188066   
17133     0.005498  -0.377144     0.000800     0.024275         4.000043   

       neglog10_p_6cov  delta_neglog10_p  sign_concordant   discovery_status  \
29083         0.027147         -1.774240            False            Neither   
29533         1.247995         -1.906950             True  3 covariates only   
19083         1.080442          1.033191            False            Neither   
16283         2.022320          1.972739             True            Neither   
26633         0.005990         -0.815940             True            Neither   
23333         2.107949          1.571798             True            Neither   
22083         0.108005         -0.800020             True            Neither   
21533         0.914396         -1.103376             True            Neither   
24983         0.022366         -0.727990            False            Neither   
28883         1.023235         -0.916110             True            Neither   
25833         0.248995         -0.437814             True            Neither   
17883         0.398527          0.182958            False            Neither   
16983         0.051630         -0.871773            False            Neither   
20533         0.119916         -0.983601             True            Neither   
18383         0.420832         -0.643201             True            Neither   
28483         0.443741          0.167809            False            Neither   
25283         0.195155         -0.196256            False            Neither   
25483         0.196791         -0.766224             True            Neither   
19633         0.652909          0.464842            False            Neither   
17133         3.096953         -0.903090             True               Both   

       union_discovery  intersection_discovery  abs_beta_change  
29083            False                   False         0.401973  
29533             True                   False         0.356761  
19083            False                   False         0.308903  
16283            False                   False         0.289274  
26633            False                   False         0.288071  
23333            False                   False         0.272627  
22083            False                   False         0.270050  
21533            False                   False         0.241664  
24983            False                   False         0.237106  
28883            False                   False         0.233282  
25833            False                   False         0.215637  
17883            False                   False         0.214832  
16983            False                   False         0.213884  
20533            False                   False         0.212878  
18383            False                   False         0.212166  
28483            False                   False         0.210940  
25283            False                   False         0.210907  
25483            False                   False         0.210673  
19633            False                   False         0.209744  
17133             True                    True         0.209312  
```

### Text output (stream)

```text

Three-covariate-only discoveries: 306

```

### Text output (text/plain)

```text
      day pvalue_type target_name program_name  beta_3cov  pvalue_3cov  \
16296  d0   empirical     C7orf26            6   0.250966       0.0003   
29533  d0   empirical        XRN1            4  -0.625779       0.0007   
25297  d0   empirical       RNF40            7   0.481823       0.0009   
17897  d0   empirical     DNAJC17            7   0.522581       0.0013   
26196  d0   empirical     SMARCA4            6   0.173361       0.0015   
28732  d0   empirical       TSTD2           39   0.512921       0.0017   
18447  d0   empirical       EIF4E            7   0.429164       0.0018   
26513  d0   empirical        SMG5           21  -0.425366       0.0018   
22497  d0   empirical         NF2            7   0.408662       0.0019   
27048  d0   empirical       SRP68            8  -0.401059       0.0019   
16598  d0   empirical     CCDC101            8  -0.378153       0.0022   
26983  d0   empirical       SRP54            4  -0.475947       0.0025   
17147  d0   empirical     CHORDC1            7   0.370661       0.0026   
29095  d0   empirical        VWA9           50   0.252623       0.0026   
29535  d0   empirical        XRN1           41  -0.442062       0.0026   
18448  d0   empirical       EIF4E            8  -0.352337       0.0032   
20033  d0   empirical       HDAC2            4  -0.422012       0.0033   
20401  d0   empirical       HUWE1           10   0.398122       0.0033   
20763  d0   empirical       KDM3B           21  -0.381843       0.0033   
19624  d0   empirical         GPI           31   0.317051       0.0035   

       qvalue_3cov  beta_6cov  pvalue_6cov  qvalue_6cov  neglog10_p_3cov  \
16296     0.011951   0.164765     0.006599     0.105054         3.522922   
29533     0.022237  -0.269018     0.056494     0.370659         3.154945   
25297     0.026499   0.355054     0.010699     0.145451         3.045801   
17897     0.034586   0.421795     0.006299     0.102245         2.886100   
26196     0.038597   0.128182     0.009399     0.133887         2.823952   
28732     0.042354   0.452720     0.007899     0.119369         2.769595   
18447     0.043819   0.297747     0.027797     0.259673         2.744771   
26513     0.043819  -0.361217     0.006399     0.103642         2.744771   
22497     0.045364   0.338183     0.006199     0.100841         2.721290   
27048     0.045364  -0.282386     0.010899     0.147096         2.721290   
16598     0.050503  -0.276691     0.007399     0.113894         2.657621   
26983     0.055180  -0.283458     0.027897     0.260118         2.602103   
17147     0.056467   0.259251     0.022498     0.230548         2.585070   
29095     0.056467   0.227049     0.008999     0.130944         2.585070   
29535     0.056467  -0.379997     0.010099     0.141291         2.585070   
18448     0.064864  -0.232393     0.035796     0.297140         2.494893   
20033     0.065817  -0.250051     0.031697     0.278305         2.481529   
20401     0.065817   0.348422     0.011399     0.151106         2.481529   
20763     0.065817  -0.328917     0.009699     0.137372         2.481529   
19624     0.068342   0.279415     0.009099     0.131373         2.455975   

       neglog10_p_6cov  delta_neglog10_p  sign_concordant   discovery_status  \
16296         2.180499         -1.342423             True  3 covariates only   
29533         1.247995         -1.906950             True  3 covariates only   
25297         1.970660         -1.075141             True  3 covariates only   
17897         2.200703         -0.685397             True  3 covariates only   
26196         2.026916         -0.797037             True  3 covariates only   
28732         2.102416         -0.667178             True  3 covariates only   
18447         1.555999         -1.188772             True  3 covariates only   
26513         2.193863         -0.550907             True  3 covariates only   
22497         2.207652         -0.513638             True  3 covariates only   
27048         1.962617         -0.758673             True  3 covariates only   
16598         2.130812         -0.526809             True  3 covariates only   
26983         1.554439         -1.047664             True  3 covariates only   
17147         1.647861         -0.937209             True  3 covariates only   
29095         2.045801         -0.539269             True  3 covariates only   
29535         1.995722         -0.589348             True  3 covariates only   
18448         1.446160         -1.048733             True  3 covariates only   
20033         1.498984         -0.982545             True  3 covariates only   
20401         1.943139         -0.538391             True  3 covariates only   
20763         2.013272         -0.468258             True  3 covariates only   
19624         2.041002         -0.414973             True  3 covariates only   

       union_discovery  intersection_discovery  
16296             True                   False  
29533             True                   False  
25297             True                   False  
17897             True                   False  
26196             True                   False  
28732             True                   False  
18447             True                   False  
26513             True                   False  
22497             True                   False  
27048             True                   False  
16598             True                   False  
26983             True                   False  
17147             True                   False  
29095             True                   False  
29535             True                   False  
18448             True                   False  
20033             True                   False  
20401             True                   False  
20763             True                   False  
19624             True                   False  
```

### Text output (stream)

```text

Six-covariate-only discoveries: 262

```

### Text output (text/plain)

```text
      day pvalue_type target_name program_name  beta_3cov  pvalue_3cov  \
17048  d0   empirical       CENPC            8  -0.481991     0.008899   
17122  d0   empirical     CHORDC1            3  -0.319641     0.009299   
18471  d0   empirical      EIF4G1           29  -0.432926     0.006499   
26560  d0   empirical     SNRNP25           19  -0.387250     0.007599   
29748  d0   empirical       ZC3H8            8  -0.361197     0.016098   
17683  d0   empirical        DBR1            4   0.318938     0.040496   
16096  d0   empirical    C16orf93            6   0.196593     0.010999   
18397  d0   empirical       EIF3M            7   0.396512     0.025997   
16682  d0   empirical        CCNC           39  -0.383787     0.014899   
22231  d0   empirical       NAA25           38  -0.349095     0.006899   
22281  d0   empirical       NCBP1           38  -0.391215     0.007899   
25786  d0   empirical         SF1           42   0.341872     0.006199   
19451  d0   empirical      GNB2L1           10   0.369311     0.007299   
26240  d0   empirical     SMARCA5           46  -0.288807     0.013599   
16299  d0   empirical     C7orf26            9  -0.319533     0.022898   
18353  d0   empirical       EIF3M           12  -0.407007     0.009999   
19896  d0   empirical      GTF2H4            6   0.166620     0.014399   
21860  d0   empirical       MEMO1           19  -0.385631     0.015498   
18172  d0   empirical      DROSHA            3   0.362537     0.008199   
29197  d0   empirical       WAPAL            7  -0.268930     0.059194   

       qvalue_3cov  beta_6cov  pvalue_6cov  qvalue_6cov  neglog10_p_3cov  \
17048     0.127620  -0.568354       0.0007     0.022002         2.050653   
17122     0.131958  -0.367014       0.0013     0.034710         2.031560   
18471     0.104692  -0.503465       0.0013     0.034710         2.187130   
26560     0.115422  -0.431657       0.0013     0.034710         2.119230   
29748     0.186666  -0.417780       0.0016     0.040131         1.793218   
17683     0.315251   0.362632       0.0019     0.045804         1.392588   
16096     0.147113   0.214749       0.0020     0.047675         1.958651   
18397     0.249107   0.504363       0.0022     0.050117         1.585070   
16682     0.177732  -0.477987       0.0024     0.053448         1.826857   
22231     0.108553  -0.381040       0.0024     0.053448         2.161194   
22281     0.118409  -0.420957       0.0024     0.053448         2.102416   
25786     0.101173   0.353037       0.0028     0.060196         2.207652   
19451     0.112471   0.387106       0.0030     0.062511         2.136721   
26240     0.168991  -0.331586       0.0032     0.065488         1.866505   
16299     0.230992  -0.421198       0.0033     0.067165         1.640208   
18353     0.139369  -0.474807       0.0033     0.067165         2.000043   
19896     0.174706   0.179170       0.0037     0.072058         1.841681   
21860     0.182551  -0.425384       0.0037     0.072058         1.809712   
18172     0.121681   0.396316       0.0038     0.073240         2.086230   
29197     0.375316  -0.363206       0.0038     0.073240         1.227722   

       neglog10_p_6cov  delta_neglog10_p  sign_concordant   discovery_status  \
17048         3.154945          1.104292             True  6 covariates only   
17122         2.886100          0.854540             True  6 covariates only   
18471         2.886100          0.698970             True  6 covariates only   
26560         2.886100          0.766870             True  6 covariates only   
29748         2.795923          1.002706             True  6 covariates only   
17683         2.721290          1.328701             True  6 covariates only   
16096         2.699013          0.740363             True  6 covariates only   
18397         2.657621          1.072551             True  6 covariates only   
16682         2.619832          0.792975             True  6 covariates only   
22231         2.619832          0.458638             True  6 covariates only   
22281         2.619832          0.517416             True  6 covariates only   
25786         2.552885          0.345234             True  6 covariates only   
19451         2.522922          0.386202             True  6 covariates only   
26240         2.494893          0.628389             True  6 covariates only   
16299         2.481529          0.841322             True  6 covariates only   
18353         2.481529          0.481486             True  6 covariates only   
19896         2.431842          0.590161             True  6 covariates only   
21860         2.431842          0.622130             True  6 covariates only   
18172         2.420260          0.334030             True  6 covariates only   
29197         2.420260          1.192538             True  6 covariates only   

       union_discovery  intersection_discovery  
17048             True                   False  
17122             True                   False  
18471             True                   False  
26560             True                   False  
29748             True                   False  
17683             True                   False  
16096             True                   False  
18397             True                   False  
16682             True                   False  
22231             True                   False  
22281             True                   False  
25786             True                   False  
19451             True                   False  
26240             True                   False  
16299             True                   False  
18353             True                   False  
19896             True                   False  
21860             True                   False  
18172             True                   False  
29197             True                   False  
```

### Text output (stream)

```text

Interpretation guide
--------------------
1. Lower NTC abs_error_at_0_05, QQ RMSE, tail QQ RMSE, and KS distance indicate better calibration.
2. Effect Pearson/Spearman quantify rank and linear stability; CCC and slope additionally detect scale attenuation or shift.
3. A six/three median absolute effect ratio below 1 and a regression slope below 1 indicate overall attenuation in the six-covariate model.
4. Discovery Jaccard and model-specific discovery counts quantify inferential stability.
5. Inspect model-specific discoveries and the Bland-Altman plot rather than relying on one global correlation.
6. This comparison identifies the aggregate impact of adding biological_sample, guide_umi_counts, and n_genes_by_counts; it does not isolate each term.

Output directories:
  Tables: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/results/ipsc_ec_3cov_vs_6cov_comparison/k_50.dt_0_2/tables
  Plots:  /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/QQ.sensitivity.analysis/results/ipsc_ec_3cov_vs_6cov_comparison/k_50.dt_0_2/plots

PDF and SVG files are the preferred Illustrator-editable outputs.

```

### Text output (error)

```text
[1;31mThe Kernel crashed while executing code in the current cell or a previous cell. 

[1;31mPlease review the code in the cell(s) to identify a possible cause of the failure. 

[1;31mClick <a href='https://aka.ms/vscodeJupyterKernelCrash'>here</a> for more info. 

[1;31mView Jupyter <a href='command:jupyter.viewOutput'>log</a> for further details.
```

