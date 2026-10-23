# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/calculate_missing_numbers_1_to_7_sherlock.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/calculate_missing_numbers_1_to_7_sherlock.ipynb`
- SHA-256: `4321bc2d1cd5799bb2806bf540f4e24dbfa01ec9b463fc080712e8af590bb970`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0


# Calculate missing manuscript numbers on Sherlock — corrected QC section

This replacement notebook calculates the CC-Perturb-seq sequencing and guide-capture QC values without silently skipping fields:

1. average transcriptome UMIs per post-QC cell;
2. average transcriptome reads per post-QC cell, when a true read-count source is supplied;
3. average guide UMIs per assigned guide–cell pair.

The AnnData path and expected QC fields are fixed explicitly from the supplied analysis notebooks. The notebook stops with a clear error if the UMI or guide fields are absent. Transcriptome **read** counts are handled separately because RNA UMIs cannot be substituted for reads.


## Code cell 1

```python

from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from IPython.display import display, Markdown

warnings.filterwarnings("ignore")
pd.set_option("display.max_columns", 200)
pd.set_option("display.max_rows", 200)

class cfg:
    # Primary integrated, post-QC CC-Perturb-seq AnnData.
    ADATA_H5AD = Path(
        "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
        "2.preprocessing_tonys_pipeline/3.integrated_days/"
        "3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad"
    )

    OUTDIR = Path("calculate_missing_numbers_1_to_7_outputs")
    EXPECTED_POST_QC_CELLS = 121_616

    # Explicit fields used in the established programDE/CC-Perturb-seq analysis.
    GEX_UMI_OBS_COL = "total_counts"
    GUIDE_TOTAL_UMI_OBS_COL = "guide_umi_counts"
    N_ASSIGNED_GUIDES_OBS_COL = "guides_per_cell"
    DAY_OBS_COL = "day"
    SAMPLE_OBS_COL = "biological_sample"

    # Optional consistency checks only.
    GUIDE_ASSIGNMENT_OBSM_KEY = "guide_assignment"
    GUIDE_NAMES_UNS_KEY = "guide_names"
    VERIFY_GUIDE_ASSIGNMENT_MATRIX = False

    # A true transcriptome read-count source is required for the read metric.
    # Option A: set this to a numeric per-cell field already present in adata.obs.
    GEX_READ_OBS_COL = None

    # Option B: provide a CSV with one row per library/sample. It may contain:
    #   total_gex_reads + post_qc_cells
    # or:
    #   mean_gex_reads_per_cell + post_qc_cells
    GEX_LIBRARY_METRICS_CSV = None
    LIBRARY_ID_COL = "library"
    LIBRARY_TOTAL_READS_COL = "total_gex_reads"
    LIBRARY_CELL_COUNT_COL = "post_qc_cells"
    LIBRARY_MEAN_READS_PER_CELL_COL = "mean_gex_reads_per_cell"

cfg.OUTDIR.mkdir(parents=True, exist_ok=True)
print("Configured AnnData:", cfg.ADATA_H5AD)

```

### Text output (stream)

```text
Configured AnnData: /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad

```

## Markdown cell 2


## 1. Validate the path and inspect the exact fields

This cell makes path or schema problems visible immediately. It does not continue with a guessed replacement field.


## Code cell 3

```python

if not cfg.ADATA_H5AD.exists():
    raise FileNotFoundError(
        "The integrated AnnData was not found at:\n"
        f"{cfg.ADATA_H5AD}\n\n"
        "Run this notebook on Sherlock with /oak mounted, or update cfg.ADATA_H5AD."
    )

adata = sc.read_h5ad(cfg.ADATA_H5AD, backed="r")

print(adata)
print(f"Post-QC cells: {adata.n_obs:,}")
print(f"Genes/features: {adata.n_vars:,}")

if adata.n_obs != cfg.EXPECTED_POST_QC_CELLS:
    print(
        f"WARNING: expected {cfg.EXPECTED_POST_QC_CELLS:,} post-QC cells from the manuscript, "
        f"but this object contains {adata.n_obs:,}."
    )

expected_fields = [
    cfg.GEX_UMI_OBS_COL,
    cfg.GUIDE_TOTAL_UMI_OBS_COL,
    cfg.N_ASSIGNED_GUIDES_OBS_COL,
    cfg.DAY_OBS_COL,
    cfg.SAMPLE_OBS_COL,
]

field_check = []
for field in expected_fields:
    present = field in adata.obs.columns
    numeric = pd.to_numeric(adata.obs[field], errors="coerce") if present else pd.Series(dtype=float)
    field_check.append(
        {
            "field": field,
            "present": present,
            "dtype": str(adata.obs[field].dtype) if present else "",
            "nonmissing": int(adata.obs[field].notna().sum()) if present else 0,
            "numeric_nonmissing": int(numeric.notna().sum()) if present else 0,
            "example_values": (
                ", ".join(adata.obs[field].astype(str).head(4).tolist()) if present else ""
            ),
        }
    )

field_check_df = pd.DataFrame(field_check)
display(field_check_df)

print("adata.obsm keys:", list(adata.obsm.keys()))
print("adata.uns keys containing 'guide':",
      [str(k) for k in adata.uns.keys() if "guide" in str(k).lower()])
print("adata.layers:", list(adata.layers.keys()))

field_check_df.to_csv(cfg.OUTDIR / "required_field_check.csv", index=False)

missing_required = [
    x for x in [
        cfg.GEX_UMI_OBS_COL,
        cfg.GUIDE_TOTAL_UMI_OBS_COL,
        cfg.N_ASSIGNED_GUIDES_OBS_COL,
    ]
    if x not in adata.obs.columns
]
if missing_required:
    raise KeyError(
        "Required AnnData fields are missing: "
        + ", ".join(missing_required)
        + "\nAvailable count/guide-like fields:\n"
        + "\n".join(
            map(
                str,
                [
                    c for c in adata.obs.columns
                    if any(token in str(c).lower() for token in ["count", "umi", "guide", "grna"])
                ],
            )
        )
    )

```

### Text output (stream)

```text
AnnData object with n_obs × n_vars = 121616 × 17538 backed at '/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad'
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'diff_day', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2'
    var: 'ensembl_id', 'gene_symbol'
    uns: 'guide_names', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2', 'log1p', 'neighbors'
    obsm: 'X_pca', 'X_umap', 'guide_assignment'
    layers: 'raw_counts'
    obsp: 'connectivities', 'distances'
Post-QC cells: 121,616
Genes/features: 17,538

```

### Text output (text/plain)

```text
               field  present     dtype  nonmissing  numeric_nonmissing  \
0       total_counts     True   float32      121616              121616   
1   guide_umi_counts     True   float64      121616              121616   
2    guides_per_cell     True     int64      121616              121616   
3                day     True  category      121616                   0   
4  biological_sample     True  category      121616                   0   

                              example_values  
0  3969.243, 3808.4248, 3441.2764, 805.50006  
1                   195.0, 176.0, 183.0, 3.0  
2                                 5, 5, 4, 0  
3                             d0, d0, d0, d0  
4         d0_rep1, d0_rep2, d0_rep2, d0_rep2  
```

### Text output (stream)

```text
adata.obsm keys: ['X_pca', 'X_umap', 'guide_assignment']
adata.uns keys containing 'guide': ['guide_names']
adata.layers: ['raw_counts']

```

## Markdown cell 4


## 2. Helper functions

The notebook uses arithmetic means for the manuscript placeholders and also reports medians and quantiles for QC.


## Code cell 5

```python

def get_numeric_obs(field: str, *, nonnegative: bool = True) -> pd.Series:
    values = pd.to_numeric(adata.obs[field], errors="coerce").astype(float)
    if values.notna().sum() == 0:
        raise ValueError(f"adata.obs[{field!r}] contains no numeric values.")
    if nonnegative and (values.dropna() < 0).any():
        raise ValueError(f"adata.obs[{field!r}] contains negative values.")
    return values


def summarize_vector(values, metric: str, source: str) -> dict:
    x = pd.Series(values, dtype=float).replace([np.inf, -np.inf], np.nan).dropna()
    if x.empty:
        raise ValueError(f"No usable values for {metric}.")
    return {
        "metric": metric,
        "source": source,
        "n": int(x.size),
        "mean": float(x.mean()),
        "median": float(x.median()),
        "sd": float(x.std(ddof=1)) if x.size > 1 else np.nan,
        "q05": float(x.quantile(0.05)),
        "q25": float(x.quantile(0.25)),
        "q75": float(x.quantile(0.75)),
        "q95": float(x.quantile(0.95)),
        "minimum": float(x.min()),
        "maximum": float(x.max()),
    }


def grouped_summary(frame: pd.DataFrame, group: str, metrics: list[str]) -> pd.DataFrame:
    parts = []
    for metric in metrics:
        if metric not in frame.columns:
            continue
        out = (
            frame.groupby(group, dropna=False)[metric]
            .agg(["count", "mean", "median", "std", "min", "max"])
            .reset_index()
        )
        out.insert(1, "metric", metric)
        parts.append(out)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

```

## Markdown cell 6


## 3. Average transcriptome UMIs per cell

`total_counts` is treated as the post-QC RNA UMI count used by the existing CC-Perturb-seq analyses.


## Code cell 7

```python

summary_rows = []
status_rows = []
cell_qc = pd.DataFrame(index=adata.obs_names.astype(str))

gex_umis = get_numeric_obs(cfg.GEX_UMI_OBS_COL)
cell_qc["transcriptome_umis"] = gex_umis.to_numpy()

gex_umi_summary = summarize_vector(
    gex_umis,
    "Transcriptome UMIs per post-QC cell",
    f"adata.obs[{cfg.GEX_UMI_OBS_COL!r}]",
)
summary_rows.append(gex_umi_summary)
status_rows.append(
    {
        "requested_metric": "Average transcriptome UMIs per cell",
        "status": "CALCULATED",
        "source": gex_umi_summary["source"],
        "value": gex_umi_summary["mean"],
        "note": "Arithmetic mean across all post-QC cells.",
    }
)

display(pd.DataFrame([gex_umi_summary]))

```

### Text output (text/plain)

```text
                                metric                     source       n  \
0  Transcriptome UMIs per post-QC cell  adata.obs['total_counts']  121616   

          mean     median           sd    q05         q25          q75  \
0  2472.064792  2022.8396  1702.234525  667.0  1269.04071  3194.199036   

           q95     minimum       maximum  
0  5790.832031  253.499969  14980.452148  
```

## Markdown cell 8


## 4. Average guide UMIs per guide/cell

For cell \(i\), let \(U_i\) be `guide_umi_counts` and \(G_i\) be `guides_per_cell`.

The main reported metric is

\[
\frac{\sum_i U_i}{\sum_i G_i},
\]

which is the global number of guide UMIs per assigned guide–cell pair. The notebook also reports total guide UMIs per cell and the mean of \(U_i/G_i\) across guide-positive cells.


## Code cell 9

```python

guide_umis = get_numeric_obs(cfg.GUIDE_TOTAL_UMI_OBS_COL)
assigned_guides = get_numeric_obs(cfg.N_ASSIGNED_GUIDES_OBS_COL)

cell_qc["guide_umis_total"] = guide_umis.to_numpy()
cell_qc["assigned_guides"] = assigned_guides.to_numpy()

valid_pair = guide_umis.notna() & assigned_guides.notna()
guide_positive = valid_pair & (assigned_guides > 0)

if guide_positive.sum() == 0:
    raise ValueError("No cells have guides_per_cell > 0.")

ratio_of_totals = float(
    guide_umis.loc[guide_positive].sum()
    / assigned_guides.loc[guide_positive].sum()
)

within_cell_ratio = (
    guide_umis.loc[guide_positive]
    / assigned_guides.loc[guide_positive]
)
cell_qc.loc[guide_positive.to_numpy(), "guide_umis_per_assigned_guide_in_cell"] = (
    within_cell_ratio.to_numpy()
)

guide_umis_per_cell_summary = summarize_vector(
    guide_umis,
    "Total guide UMIs per post-QC cell",
    f"adata.obs[{cfg.GUIDE_TOTAL_UMI_OBS_COL!r}]",
)
assigned_guides_summary = summarize_vector(
    assigned_guides,
    "Assigned guides per post-QC cell",
    f"adata.obs[{cfg.N_ASSIGNED_GUIDES_OBS_COL!r}]",
)
within_cell_summary = summarize_vector(
    within_cell_ratio,
    "Within-cell guide UMIs per assigned guide",
    (
        f"adata.obs[{cfg.GUIDE_TOTAL_UMI_OBS_COL!r}] / "
        f"adata.obs[{cfg.N_ASSIGNED_GUIDES_OBS_COL!r}] among guide-positive cells"
    ),
)
global_guide_metric = {
    "metric": "Guide UMIs per assigned guide-cell pair, ratio of totals",
    "source": (
        f"sum(adata.obs[{cfg.GUIDE_TOTAL_UMI_OBS_COL!r}]) / "
        f"sum(adata.obs[{cfg.N_ASSIGNED_GUIDES_OBS_COL!r}]) among guide-positive cells"
    ),
    "n": int(guide_positive.sum()),
    "mean": ratio_of_totals,
    "median": np.nan,
    "sd": np.nan,
    "q05": np.nan,
    "q25": np.nan,
    "q75": np.nan,
    "q95": np.nan,
    "minimum": np.nan,
    "maximum": np.nan,
}

summary_rows.extend([
    guide_umis_per_cell_summary,
    assigned_guides_summary,
    within_cell_summary,
    global_guide_metric,
])
status_rows.append(
    {
        "requested_metric": "Average guide UMIs per guide/cell",
        "status": "CALCULATED",
        "source": global_guide_metric["source"],
        "value": ratio_of_totals,
        "note": "Recommended metric: ratio of total guide UMIs to total assigned guide-cell pairs.",
    }
)

display(
    pd.DataFrame(
        [
            guide_umis_per_cell_summary,
            assigned_guides_summary,
            within_cell_summary,
            global_guide_metric,
        ]
    )
)

print(f"Guide-positive cells: {guide_positive.sum():,} / {adata.n_obs:,}")
print(f"Recommended guide UMIs per assigned guide-cell pair: {ratio_of_totals:,.3f}")

if cfg.VERIFY_GUIDE_ASSIGNMENT_MATRIX:
    if cfg.GUIDE_ASSIGNMENT_OBSM_KEY not in adata.obsm:
        raise KeyError(
            f"adata.obsm[{cfg.GUIDE_ASSIGNMENT_OBSM_KEY!r}] is absent."
        )
    assignment = adata.obsm[cfg.GUIDE_ASSIGNMENT_OBSM_KEY]
    print(
        "Guide-assignment matrix:",
        assignment.shape,
        assignment.dtype,
        "sparse=" + str(sparse.issparse(assignment)),
    )
    if assignment.shape[0] != adata.n_obs:
        raise ValueError("Guide-assignment matrix row count does not match adata.n_obs.")

```

### Text output (text/plain)

```text
                                              metric  \
0                  Total guide UMIs per post-QC cell   
1                   Assigned guides per post-QC cell   
2          Within-cell guide UMIs per assigned guide   
3  Guide UMIs per assigned guide-cell pair, ratio...   

                                              source       n        mean  \
0                      adata.obs['guide_umi_counts']  121616  187.750485   
1                       adata.obs['guides_per_cell']  121616    4.487600   
2  adata.obs['guide_umi_counts'] / adata.obs['gui...  102392   44.957606   
3  sum(adata.obs['guide_umi_counts']) / sum(adata...  102392   41.451107   

   median          sd   q05   q25         q75    q95  minimum  maximum  
0  101.00  252.415010   6.0  32.0  244.000000  658.0      1.0   5362.0  
1    3.00    4.749326   0.0   1.0    6.000000   14.0      0.0     60.0  
2   31.75   44.986110  13.0  21.0   52.166667  119.0      5.0   2113.0  
3     NaN         NaN   NaN   NaN         NaN    NaN      NaN      NaN  
```

### Text output (stream)

```text
Guide-positive cells: 102,392 / 121,616
Recommended guide UMIs per assigned guide-cell pair: 41.451

```

## Markdown cell 10


## 5. Average transcriptome reads per cell

This metric is calculated only from a real read-count field or sequencing/library metrics file. `total_counts` is RNA UMIs and is not used as a proxy for reads.


## Code cell 11

```python

read_metric_value = None
read_metric_source = ""
read_metric_status = ""

# First preference: explicitly configured per-cell read field.
if cfg.GEX_READ_OBS_COL is not None:
    if cfg.GEX_READ_OBS_COL not in adata.obs.columns:
        raise KeyError(
            f"Configured per-cell read field {cfg.GEX_READ_OBS_COL!r} is not in adata.obs."
        )
    gex_reads = get_numeric_obs(cfg.GEX_READ_OBS_COL)
    cell_qc["transcriptome_reads"] = gex_reads.to_numpy()
    read_summary = summarize_vector(
        gex_reads,
        "Transcriptome reads per post-QC cell",
        f"adata.obs[{cfg.GEX_READ_OBS_COL!r}]",
    )
    summary_rows.append(read_summary)
    read_metric_value = read_summary["mean"]
    read_metric_source = read_summary["source"]
    read_metric_status = "CALCULATED"

# Second preference: library-level read metrics table.
elif cfg.GEX_LIBRARY_METRICS_CSV is not None:
    metrics_path = Path(cfg.GEX_LIBRARY_METRICS_CSV)
    if not metrics_path.exists():
        raise FileNotFoundError(metrics_path)

    library_metrics = pd.read_csv(metrics_path)
    display(library_metrics)

    mean_col = cfg.LIBRARY_MEAN_READS_PER_CELL_COL
    total_col = cfg.LIBRARY_TOTAL_READS_COL
    cells_col = cfg.LIBRARY_CELL_COUNT_COL

    if mean_col in library_metrics.columns:
        means = pd.to_numeric(library_metrics[mean_col], errors="coerce")
        if cells_col in library_metrics.columns:
            weights = pd.to_numeric(library_metrics[cells_col], errors="coerce")
            ok = means.notna() & weights.notna() & (weights > 0)
            if ok.sum() == 0:
                raise ValueError("No valid mean-read/cell-count rows in library metrics.")
            read_metric_value = float(np.average(means[ok], weights=weights[ok]))
            read_metric_source = (
                f"cell-count-weighted {mean_col} from {metrics_path}"
            )
        else:
            read_metric_value = float(means.mean())
            read_metric_source = f"unweighted {mean_col} from {metrics_path}"

    elif (
        total_col in library_metrics.columns
        and cells_col in library_metrics.columns
    ):
        total_reads = pd.to_numeric(library_metrics[total_col], errors="coerce").sum()
        total_cells = pd.to_numeric(library_metrics[cells_col], errors="coerce").sum()
        if total_cells <= 0:
            raise ValueError("Total post-QC cell count in library metrics is zero.")
        read_metric_value = float(total_reads / total_cells)
        read_metric_source = (
            f"sum({total_col}) / sum({cells_col}) from {metrics_path}"
        )
    else:
        raise ValueError(
            "The library metrics file must contain either "
            f"{mean_col!r}, or both {total_col!r} and {cells_col!r}."
        )

    summary_rows.append(
        {
            "metric": "Transcriptome reads per post-QC cell",
            "source": read_metric_source,
            "n": int(len(library_metrics)),
            "mean": read_metric_value,
            "median": np.nan,
            "sd": np.nan,
            "q05": np.nan,
            "q25": np.nan,
            "q75": np.nan,
            "q95": np.nan,
            "minimum": np.nan,
            "maximum": np.nan,
        }
    )
    read_metric_status = "CALCULATED"

else:
    read_metric_status = "NEEDS_READ_METRICS_FILE"
    candidate_read_cols = [
        str(c) for c in adata.obs.columns if "read" in str(c).lower()
    ]
    print(
        "Average transcriptome reads per cell was NOT calculated.\n"
        "Reason: no true read-count field or library-level sequencing metrics file "
        "was configured.\n"
        f"Read-like adata.obs fields found: {candidate_read_cols or 'none'}\n\n"
        "Set cfg.GEX_READ_OBS_COL to a genuine per-cell read-count field, or set "
        "cfg.GEX_LIBRARY_METRICS_CSV to a library-level metrics CSV."
    )

status_rows.append(
    {
        "requested_metric": "Average transcriptome reads per cell",
        "status": read_metric_status,
        "source": read_metric_source,
        "value": read_metric_value,
        "note": (
            "RNA UMIs were not substituted for sequencing reads."
            if read_metric_status != "CALCULATED"
            else "Arithmetic mean or cell-count-weighted library mean."
        ),
    }
)

```

### Text output (stream)

```text
Average transcriptome reads per cell was NOT calculated.
Reason: no true read-count field or library-level sequencing metrics file was configured.
Read-like adata.obs fields found: none

Set cfg.GEX_READ_OBS_COL to a genuine per-cell read-count field, or set cfg.GEX_LIBRARY_METRICS_CSV to a library-level metrics CSV.

```

## Markdown cell 12


## 6. Per-day and per-sample QC summaries


## Code cell 13

```python

if cfg.DAY_OBS_COL in adata.obs.columns:
    cell_qc["day"] = adata.obs[cfg.DAY_OBS_COL].astype(str).to_numpy()
else:
    cell_qc["day"] = "unknown"

if cfg.SAMPLE_OBS_COL in adata.obs.columns:
    cell_qc["sample"] = adata.obs[cfg.SAMPLE_OBS_COL].astype(str).to_numpy()
else:
    cell_qc["sample"] = "unknown"

metric_columns = [
    c for c in [
        "transcriptome_umis",
        "transcriptome_reads",
        "guide_umis_total",
        "assigned_guides",
        "guide_umis_per_assigned_guide_in_cell",
    ]
    if c in cell_qc.columns
]

by_day = grouped_summary(cell_qc, "day", metric_columns)
by_sample = grouped_summary(cell_qc, "sample", metric_columns)

display(Markdown("### By day"))
display(by_day)
display(Markdown("### By biological sample"))
display(by_sample)

by_day.to_csv(cfg.OUTDIR / "sequencing_guide_QC_by_day.csv", index=False)
by_sample.to_csv(cfg.OUTDIR / "sequencing_guide_QC_by_sample.csv", index=False)

```

### Text output (text/plain)

```text
<IPython.core.display.Markdown object>
```

### Text output (text/markdown)

```text
### By day
```

### Text output (text/plain)

```text
   day                                 metric  count         mean  \
0   d0                     transcriptome_umis  38609  3207.869811   
1   d1                     transcriptome_umis  16925  3069.647211   
2   d2                     transcriptome_umis  21980  1904.424294   
3   d3                     transcriptome_umis  44102  1881.478418   
4   d0                       guide_umis_total  38609   183.179336   
5   d1                       guide_umis_total  16925   235.367504   
6   d2                       guide_umis_total  21980   215.386306   
7   d3                       guide_umis_total  44102   159.704911   
8   d0                        assigned_guides  38609     4.472973   
9   d1                        assigned_guides  16925     5.314682   
10  d2                        assigned_guides  21980     5.254231   
11  d3                        assigned_guides  44102     3.800916   
12  d0  guide_umis_per_assigned_guide_in_cell  32963    44.525203   
13  d1  guide_umis_per_assigned_guide_in_cell  14706    47.701963   
14  d2  guide_umis_per_assigned_guide_in_cell  18976    43.161815   
15  d3  guide_umis_per_assigned_guide_in_cell  35747    45.180611   

         median          std         min           max  
0   2802.595947  1840.850553  500.500000  14980.452148  
1   2581.354004  1996.603527  379.099976  14813.840820  
2   1570.741699  1238.457427  259.633331   8197.418945  
3   1530.398926  1268.572158  253.499969   8500.360352  
4    103.000000   239.027533    1.000000   3918.000000  
5    135.000000   295.329954    1.000000   4231.000000  
6    121.000000   269.091871    1.000000   2880.000000  
7     80.000000   232.522081    1.000000   5362.000000  
8      3.000000     4.474284    0.000000     47.000000  
9      4.000000     5.301465    0.000000     57.000000  
10     4.000000     5.374349    0.000000     60.000000  
11     2.000000     4.297464    0.000000     55.000000  
12    31.000000    45.925840    5.000000   2113.000000  
13    33.833333    46.481555    7.000000    947.000000  
14    32.000000    38.657733    5.000000    730.000000  
15    31.500000    46.527648    5.000000    928.000000  
```

### Text output (text/plain)

```text
<IPython.core.display.Markdown object>
```

### Text output (text/markdown)

```text
### By biological sample
```

### Text output (text/plain)

```text
     sample                                 metric  count         mean  \
0   d0_rep1                     transcriptome_umis  17488  2888.669391   
1   d0_rep2                     transcriptome_umis  21121  3472.164919   
2   d1_rep1                     transcriptome_umis   5856  3438.124129   
3   d1_rep2                     transcriptome_umis  11069  2874.706310   
4   d2_rep1                     transcriptome_umis  10895  2071.931235   
5   d2_rep2                     transcriptome_umis  11085  1739.788469   
6   d3_rep1                     transcriptome_umis  21377  1946.089717   
7   d3_rep2                     transcriptome_umis  22725  1820.699727   
8   d0_rep1                       guide_umis_total  17488   182.706656   
9   d0_rep2                       guide_umis_total  21121   183.570712   
10  d1_rep1                       guide_umis_total   5856   268.714993   
11  d1_rep2                       guide_umis_total  11069   217.725178   
12  d2_rep1                       guide_umis_total  10895   253.017164   
13  d2_rep2                       guide_umis_total  11085   178.400451   
14  d3_rep1                       guide_umis_total  21377   145.851055   
15  d3_rep2                       guide_umis_total  22725   172.736986   
16  d0_rep1                        assigned_guides  17488     4.606473   
17  d0_rep2                        assigned_guides  21121     4.362435   
18  d1_rep1                        assigned_guides   5856     5.638661   
19  d1_rep2                        assigned_guides  11069     5.143283   
20  d2_rep1                        assigned_guides  10895     5.820652   
21  d2_rep2                        assigned_guides  11085     4.697519   
22  d3_rep1                        assigned_guides  21377     3.707910   
23  d3_rep2                        assigned_guides  22725     3.888405   
24  d0_rep1  guide_umis_per_assigned_guide_in_cell  15144    42.622708   
25  d0_rep2  guide_umis_per_assigned_guide_in_cell  17819    46.142095   
26  d1_rep1  guide_umis_per_assigned_guide_in_cell   5128    51.970029   
27  d1_rep2  guide_umis_per_assigned_guide_in_cell   9578    45.416867   
28  d2_rep1  guide_umis_per_assigned_guide_in_cell   9770    46.596560   
29  d2_rep2  guide_umis_per_assigned_guide_in_cell   9206    39.516643   
30  d3_rep1  guide_umis_per_assigned_guide_in_cell  17611    41.799422   
31  d3_rep2  guide_umis_per_assigned_guide_in_cell  18136    48.463921   

         median          std         min           max  
0   2433.188599  1789.828474  517.720459  14980.452148  
1   3146.939453  1840.665213  500.500000  14888.178711  
2   2988.112915  2011.528179  383.555481  14735.392578  
3   2349.726807  1960.943090  379.099976  14813.840820  
4   1721.215942  1279.908543  259.916687   8197.418945  
5   1416.250000  1173.297910  259.633331   8047.102539  
6   1595.265869  1287.721373  253.499969   8479.798828  
7   1470.000000  1247.268034  259.500000   8500.360352  
8    102.000000   237.471512    1.000000   3635.000000  
9    104.000000   240.313192    1.000000   3918.000000  
10   159.000000   327.049307    1.000000   3685.000000  
11   122.000000   275.471805    1.000000   4231.000000  
12   153.000000   298.750579    1.000000   2880.000000  
13    95.000000   230.443434    1.000000   2661.000000  
14    75.000000   211.022977    1.000000   5303.000000  
15    86.000000   250.377023    1.000000   5362.000000  
16     3.000000     4.576424    0.000000     37.000000  
17     3.000000     4.384948    0.000000     47.000000  
18     4.000000     5.547919    0.000000     51.000000  
19     4.000000     5.158359    0.000000     57.000000  
20     4.000000     5.688163    0.000000     60.000000  
21     3.000000     4.984863    0.000000     55.000000  
22     2.000000     4.127495    0.000000     43.000000  
23     2.000000     4.449748    0.000000     55.000000  
24    29.825758    41.864703    5.000000    776.000000  
25    32.000000    49.057256    6.000000   2113.000000  
26    37.000000    49.596005    7.000000    858.000000  
27    32.000000    44.560002    7.000000    947.000000  
28    34.733333    40.756551    5.000000    730.000000  
29    29.500000    35.942784    5.500000    646.500000  
30    29.090909    41.977770    5.000000    695.000000  
31    33.857143    50.339744    5.000000    928.000000  
```

## Markdown cell 14


## 7. Final manuscript-ready values and explicit status report


## Code cell 15

```python

summary_df = pd.DataFrame(summary_rows)
status_df = pd.DataFrame(status_rows)

display(Markdown("### Requested-metric status"))
display(status_df)

display(Markdown("### Complete numerical summaries"))
display(summary_df)

recommended = {
    "post_QC_cells": int(adata.n_obs),
    "average_transcriptome_UMIs_per_cell": float(gex_umi_summary["mean"]),
    "average_transcriptome_reads_per_cell": (
        float(read_metric_value) if read_metric_value is not None else None
    ),
    "average_guide_UMIs_per_assigned_guide_cell_pair": float(ratio_of_totals),
    "mean_within_cell_guide_UMIs_per_assigned_guide": float(
        within_cell_summary["mean"]
    ),
    "average_total_guide_UMIs_per_cell": float(
        guide_umis_per_cell_summary["mean"]
    ),
    "average_assigned_guides_per_cell": float(
        assigned_guides_summary["mean"]
    ),
}

summary_df.to_csv(
    cfg.OUTDIR / "CC_Perturbseq_sequencing_guide_QC_summary.csv",
    index=False,
)
status_df.to_csv(
    cfg.OUTDIR / "requested_metric_status.csv",
    index=False,
)
with open(cfg.OUTDIR / "manuscript_QC_values.json", "w") as handle:
    json.dump(recommended, handle, indent=2)

print(json.dumps(recommended, indent=2))

umi_text = f"{recommended['average_transcriptome_UMIs_per_cell']:,.0f}"
guide_text = f"{recommended['average_guide_UMIs_per_assigned_guide_cell_pair']:,.1f}"
read_text = (
    f"{recommended['average_transcriptome_reads_per_cell']:,.0f}"
    if recommended["average_transcriptome_reads_per_cell"] is not None
    else "[READ METRICS REQUIRED]"
)

manuscript_sentence = (
    f"Across {recommended['post_QC_cells']:,} post-QC cells, libraries contained "
    f"an average of {umi_text} transcriptome UMIs per cell, {read_text} "
    f"transcriptome reads per cell, and {guide_text} guide UMIs per assigned "
    "guide-cell pair."
)
print("\nManuscript-ready sentence:\n")
print(manuscript_sentence)

(cfg.OUTDIR / "manuscript_ready_QC_sentence.txt").write_text(
    manuscript_sentence + "\n"
)

```

### Text output (text/plain)

```text
<IPython.core.display.Markdown object>
```

### Text output (text/markdown)

```text
### Requested-metric status
```

### Text output (text/plain)

```text
                       requested_metric                   status  \
0   Average transcriptome UMIs per cell               CALCULATED   
1     Average guide UMIs per guide/cell               CALCULATED   
2  Average transcriptome reads per cell  NEEDS_READ_METRICS_FILE   

                                              source        value  \
0                          adata.obs['total_counts']  2472.064792   
1  sum(adata.obs['guide_umi_counts']) / sum(adata...    41.451107   
2                                                             NaN   

                                                note  
0          Arithmetic mean across all post-QC cells.  
1  Recommended metric: ratio of total guide UMIs ...  
2  RNA UMIs were not substituted for sequencing r...  
```

### Text output (text/plain)

```text
<IPython.core.display.Markdown object>
```

### Text output (text/markdown)

```text
### Complete numerical summaries
```

### Text output (text/plain)

```text
                                              metric  \
0                Transcriptome UMIs per post-QC cell   
1                  Total guide UMIs per post-QC cell   
2                   Assigned guides per post-QC cell   
3          Within-cell guide UMIs per assigned guide   
4  Guide UMIs per assigned guide-cell pair, ratio...   

                                              source       n         mean  \
0                          adata.obs['total_counts']  121616  2472.064792   
1                      adata.obs['guide_umi_counts']  121616   187.750485   
2                       adata.obs['guides_per_cell']  121616     4.487600   
3  adata.obs['guide_umi_counts'] / adata.obs['gui...  102392    44.957606   
4  sum(adata.obs['guide_umi_counts']) / sum(adata...  102392    41.451107   

      median           sd    q05         q25          q75          q95  \
0  2022.8396  1702.234525  667.0  1269.04071  3194.199036  5790.832031   
1   101.0000   252.415010    6.0    32.00000   244.000000   658.000000   
2     3.0000     4.749326    0.0     1.00000     6.000000    14.000000   
3    31.7500    44.986110   13.0    21.00000    52.166667   119.000000   
4        NaN          NaN    NaN         NaN          NaN          NaN   

      minimum       maximum  
0  253.499969  14980.452148  
1    1.000000   5362.000000  
2    0.000000     60.000000  
3    5.000000   2113.000000  
4         NaN           NaN  
```

### Text output (stream)

```text
{
  "post_QC_cells": 121616,
  "average_transcriptome_UMIs_per_cell": 2472.0647921069735,
  "average_transcriptome_reads_per_cell": null,
  "average_guide_UMIs_per_assigned_guide_cell_pair": 41.4511070719212,
  "mean_within_cell_guide_UMIs_per_assigned_guide": 44.95760642459217,
  "average_total_guide_UMIs_per_cell": 187.75048513353505,
  "average_assigned_guides_per_cell": 4.487600315747928
}

Manuscript-ready sentence:

Across 121,616 post-QC cells, libraries contained an average of 2,472 transcriptome UMIs per cell, [READ METRICS REQUIRED] transcriptome reads per cell, and 41.5 guide UMIs per assigned guide-cell pair.

```

### Text output (text/plain)

```text
203
```

## Markdown cell 16


## 8. Close the backed AnnData handle


## Code cell 17

```python

try:
    adata.file.close()
except Exception:
    pass

print("Finished. Outputs:", cfg.OUTDIR.resolve())

```

### Text output (stream)

```text
Finished. Outputs: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/calculate_missing_numbers_1_to_7_outputs

```

