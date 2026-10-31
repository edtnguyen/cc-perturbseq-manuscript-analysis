# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/K562_hiPSC_ParseV3_plotnine_figure_reproduction.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/K562_hiPSC_ParseV3_plotnine_figure_reproduction.ipynb`
- SHA-256: `404daf543b4b333cf0fcbf3e553a30fad6895dc38382a27bd492f72e0d4eadb2`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# K562–hiPSC Parse v3 circular-capture figure reproduction

This notebook reproduces the quantitative plots in **“Testing CC-Perturb-seq with Parse v3 kit”** using **plotnine** and the original Sherlock data paths.

It covers the slide panels for:

- transcriptome QC and sequencing saturation;
- transcriptome yield across RT-primer conditions;
- guide UMI capture across CS1, DC1, dual-primer and no-guide conditions;
- iPSC-guide versus K562-guide mixing/specificity;
- K562/hiPSC UMAPs and marker expression;
- species-specific guide UMIs on the UMAP;
- guide capture stratified by inferred cell type;
- optional paper-style panels based on top-guide UMIs and expected-species guide metrics.

The raw-data workflow follows the attached analysis notebook and low-memory export script. The notebook caches processed tables so subsequent plotting runs do not reread the large guide-count table.

## Markdown cell 1

## 1. Environment

Recommended Sherlock environment packages:

```bash
mamba install -c conda-forge \
  python=3.10 pandas numpy scipy scanpy anndata python-igraph leidenalg \
  plotnine mizani pyarrow jupyterlab
```

SVG output is configured so text remains editable in Adobe Illustrator.

## Code cell 2

```python
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence
import json
import os
import re
import sys
import warnings
import zipfile

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import mannwhitneyu

try:
    import anndata as ad
    import scanpy as sc
except ImportError:
    ad = None
    sc = None
import matplotlib as mpl

from plotnine import (
    ggplot, aes, annotate,
    geom_abline, geom_boxplot, geom_histogram, geom_jitter,
    geom_line, geom_point, geom_text, geom_violin,
    stat_bin_2d,
    facet_grid, facet_wrap,
    coord_equal, coord_flip,
    labs,
    scale_color_cmap, scale_color_manual,
    scale_fill_cmap, scale_fill_manual,
    scale_x_continuous, scale_x_log10,
    scale_y_continuous, scale_y_log10,
    position_dodge, position_jitterdodge,
    theme, theme_bw,
    element_blank, element_line, element_rect, element_text,
)
from mizani.formatters import comma_format, percent_format
from IPython.display import display

# Keep vector text as text in SVG/PDF.
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

warnings.filterwarnings("ignore", category=FutureWarning)

print("Python:", sys.version.split()[0])
print("Scanpy:", sc.__version__ if sc is not None else "not installed (cached-table mode only)")
try:
    import plotnine
    print("plotnine:", plotnine.__version__)
except Exception:
    pass
```

### Text output (stream)

```text
Python: 3.11.14
Scanpy: 1.11.5
plotnine: 0.15.4

```

## Markdown cell 3

## 2. Sherlock paths and analysis settings

## Code cell 4

```python
# Exact project path recovered from the attached analysis notebook and export records.
PROJECT_ROOT = Path(
    "/oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/"
    "250730_eTN10_Parsev3_mixediPSC_K562"
)

DGE_FILTERED_DIR = PROJECT_ROOT / "splitpipe_out/WT/sub1/all-sample/DGE_filtered"
GUIDE_COUNT_TSV = PROJECT_ROOT / "splitpipe_out/CRISPR_combined/grnaxcell_count.tsv"
ALL_GUIDES_CSV = (
    PROJECT_ROOT
    / "splitpipe_out/CRISPR/sub1/all-sample/guide_RNAs_unfiltered/all_guides.csv"
)

SATURATION_FILES = {
    sample: PROJECT_ROOT / f"splitpipe_out/WT/sub1/{sample}/report/tscp_counts_subsampled.csv"
    for sample in ["CS1", "DC1", "Dual", "no_guide"]
}

# Optional gene-saturation files. The notebook uses these only when they exist.
GENE_SATURATION_CANDIDATES = {
    sample: [
        PROJECT_ROOT / f"splitpipe_out/WT/sub1/{sample}/report/gene_counts_subsampled.csv",
        PROJECT_ROOT / f"splitpipe_out/WT/sub1/{sample}/report/genes_counts_subsampled.csv",
    ]
    for sample in ["CS1", "DC1", "Dual", "no_guide"]
}

# On Sherlock, outputs are written under the experiment. Elsewhere, use the notebook directory.
RUN_ROOT = PROJECT_ROOT if PROJECT_ROOT.exists() else Path.cwd() / "k562_parsev3_plotnine_run"
CACHE_DIR = RUN_ROOT / "plotnine_reproduction/cache"
FIGURE_DIR = RUN_ROOT / "plotnine_reproduction/figures"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# The uploaded/precomputed plot-input bundle is an optional fallback.
BUNDLE_CANDIDATES = [
    Path.cwd() / "k562_parsev3_figure_panel_plot_input_bundle.zip",
    Path("/mnt/data/k562_parsev3_figure_panel_plot_input_bundle.zip"),
    Path(
        "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/"
        "figure_1/Figure1_K562_ParseV3_export_v2_lowmem/"
        "k562_parsev3_figure_panel_plot_input_bundle.zip"
    ),
]

SAMPLE_ORDER = ["CS1", "DC1", "Dual", "no_guide"]
PRIMER_SAMPLES = ["CS1", "DC1", "Dual"]
SAMPLE_COLORS = {
    "CS1": "#4E79A7",
    "DC1": "#F28E2B",
    "Dual": "#E15759",
    "no_guide": "#76B7B2",
}
CELLTYPE_COLORS = {"K562": "#4E79A7", "hiPSC": "#F28E2B"}

# QC and analysis settings copied from the original notebook.
MAX_GENES = 5_000
MAX_COUNTS = 5_000
MAX_PCT_MT = 40
N_NEIGHBORS = 20
N_PCS = 30
LEIDEN_RESOLUTION = 0.2
GUIDE_UMI_THRESHOLD = 3  # per-guide threshold used for species-mixing panels
GUIDE_CHUNK_SIZE = 1_000_000
RANDOM_SEED = 0

# SVG keeps text editable; add "pdf" here if a PDF copy is also desired.
EXPORT_FORMATS = ("png", "svg")

# Set True only when a fresh rebuild is required.
FORCE_RECOMPUTE_GEX = False
FORCE_RECOMPUTE_GUIDES = False

print("Raw GEX available:", DGE_FILTERED_DIR.exists())
print("Raw guide table available:", GUIDE_COUNT_TSV.exists())
print("Cache directory:", CACHE_DIR)
print("Figure directory:", FIGURE_DIR)
```

### Text output (stream)

```text
Raw GEX available: True
Raw guide table available: True
Cache directory: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/cache
Figure directory: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures

```

## Markdown cell 5

## 3. General helpers

The workflow first looks for cached outputs. If a cache is absent and the Sherlock raw paths exist, it recomputes the data. If raw data are unavailable, it can load the attached plot-input bundle for all non-UMAP panels.

## Code cell 6

```python
def normalize_cell_id(value: object) -> str:
    """Remove a split-pipe sublibrary suffix such as '__s1'."""
    return re.sub(r"__s\d+$", "", str(value))


def ordered_sample(series: pd.Series) -> pd.Categorical:
    return pd.Categorical(series.astype(str), categories=SAMPLE_ORDER, ordered=True)


def first_existing(paths: Sequence[Path]) -> Path | None:
    return next((p for p in paths if p.exists()), None)


def theme_publication(base_size: int = 10):
    return (
        theme_bw(base_size=base_size)
        + theme(
            panel_grid_minor=element_blank(),
            panel_grid_major=element_line(color="#D8D8D8", size=0.45),
            panel_border=element_rect(color="black", fill=None, size=0.7),
            axis_title=element_text(size=base_size + 1),
            axis_text=element_text(size=base_size),
            legend_title=element_text(size=base_size),
            legend_text=element_text(size=max(base_size - 1, 7)),
            strip_background=element_rect(fill="white", color="black", size=0.7),
            strip_text=element_text(size=base_size, weight="bold"),
            plot_title=element_text(size=base_size + 3, weight="bold"),
            plot_subtitle=element_text(size=base_size + 1),
        )
    )


def save_plot(plot, stem: str, width: float, height: float, dpi: int = 300) -> None:
    """Save PNG plus editable-text SVG/PDF versions."""
    for ext in EXPORT_FORMATS:
        path = FIGURE_DIR / f"{stem}.{ext}"
        plot.save(path, width=width, height=height, units="in", dpi=dpi, verbose=False)
    print(f"Saved: {FIGURE_DIR / stem} in {EXPORT_FORMATS}")


def holm_adjust(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    idx = np.where(valid)[0]
    order = np.argsort(p[idx])
    sorted_idx = idx[order]
    sorted_p = p[sorted_idx]
    m = len(sorted_p)
    adjusted = np.array([(m - k) * pv for k, pv in enumerate(sorted_p)], dtype=float)
    adjusted = np.maximum.accumulate(adjusted)
    adjusted = np.minimum(adjusted, 1.0)
    out[sorted_idx] = adjusted
    return out
```

## Markdown cell 7

## 4. Load and preprocess gene-expression data

## Code cell 8

```python
def infer_dge_files(dge_dir: Path) -> tuple[Path, Path, Path]:
    matrix = first_existing([dge_dir / "count_matrix.mtx", dge_dir / "DGE.mtx"])
    genes = first_existing([dge_dir / "all_genes.csv", dge_dir / "genes.csv"])
    cells = first_existing([dge_dir / "cell_metadata.csv", dge_dir / "cells.csv"])
    if not all([matrix, genes, cells]):
        raise FileNotFoundError(f"Missing matrix/genes/cells in {dge_dir}")
    return matrix, genes, cells


def load_raw_gex(dge_dir: Path):
    if ad is None or sc is None:
        raise ImportError("Raw GEX processing requires scanpy and anndata.")
    matrix_path, genes_path, cells_path = infer_dge_files(dge_dir)
    print("Loading:", matrix_path)
    X = mmread(matrix_path).tocsr()
    genes = pd.read_csv(genes_path)
    cells = pd.read_csv(cells_path)

    gene_col = "gene_name" if "gene_name" in genes.columns else genes.columns[0]
    keep = genes[gene_col].notna() & genes[gene_col].astype(str).ne("nan")

    n_cells, n_genes = len(cells), len(genes)
    if X.shape == (n_genes, n_cells):
        X = X.T.tocsr()
    elif X.shape != (n_cells, n_genes):
        raise ValueError(
            f"Matrix shape {X.shape} does not match cells={n_cells}, genes={n_genes}"
        )

    X = X[:, np.where(keep.to_numpy())[0]].tocsr()
    genes = genes.loc[keep].reset_index(drop=True)

    adata = ad.AnnData(X=X, obs=cells.copy(), var=genes.copy())
    cell_col = "bc_wells" if "bc_wells" in cells.columns else cells.columns[0]
    adata.obs_names = pd.Index([normalize_cell_id(x) for x in cells[cell_col]])
    adata.obs_names_make_unique()
    adata.var_names = pd.Index(genes[gene_col].astype(str))
    adata.var_names_make_unique()
    return adata


def qc_gex(adata):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    before = adata.obs.copy()
    before["cell_id"] = adata.obs_names

    keep = (
        (adata.obs["n_genes_by_counts"] < MAX_GENES)
        & (adata.obs["total_counts"] < MAX_COUNTS)
        & (adata.obs["pct_counts_mt"] < MAX_PCT_MT)
    )
    adata = adata[keep].copy()
    print(f"Cells before QC: {len(before):,}; after QC: {adata.n_obs:,}")
    return adata, before


def extract_gene_vector(adata, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    x = adata[:, gene].X
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def compute_umap_plot_data(adata_counts) -> pd.DataFrame:
    """Notebook-like normalization, HVG selection, PCA, neighbors, UMAP and Leiden."""
    adata = adata_counts.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    marker_genes = ["SOX2", "NANOG", "POU5F1", "LIN28A", "GATA1", "KLF1", "HBB", "HBA1", "HBA2"]
    marker_df = pd.DataFrame(index=adata.obs_names)
    for gene in marker_genes:
        marker_df[gene] = extract_gene_vector(adata, gene)

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.25)
    if int(adata.var["highly_variable"].sum()) < 50:
        raise RuntimeError("Too few highly variable genes were identified.")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.regress_out(adata_hvg, ["total_counts"])
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, svd_solver="arpack", random_state=RANDOM_SEED)
    sc.pp.neighbors(adata_hvg, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=RANDOM_SEED)
    sc.tl.umap(adata_hvg, random_state=RANDOM_SEED)
    sc.tl.leiden(
        adata_hvg,
        resolution=LEIDEN_RESOLUTION,
        random_state=RANDOM_SEED,
        key_added="leiden",
    )

    # Label whole Leiden clusters using lineage-marker scores, matching the slide logic.
    marker_df["hipsc_score"] = marker_df[["SOX2", "NANOG", "POU5F1", "LIN28A"]].mean(axis=1, skipna=True)
    marker_df["k562_score"] = marker_df[["GATA1", "KLF1", "HBB", "HBA1", "HBA2"]].mean(axis=1, skipna=True)
    marker_df["leiden"] = adata_hvg.obs["leiden"].astype(str).to_numpy()

    cluster_scores = marker_df.groupby("leiden")[["hipsc_score", "k562_score"]].mean()
    cluster_to_type = {
        cluster: ("K562" if row["k562_score"] > row["hipsc_score"] else "hiPSC")
        for cluster, row in cluster_scores.iterrows()
    }
    print("Leiden cluster labels:", cluster_to_type)

    out = adata_hvg.obs.copy()
    out["cell_id"] = adata_hvg.obs_names
    out["UMAP1"] = adata_hvg.obsm["X_umap"][:, 0]
    out["UMAP2"] = adata_hvg.obsm["X_umap"][:, 1]
    out["leiden"] = adata_hvg.obs["leiden"].astype(str).to_numpy()
    out["cell_type_umap"] = out["leiden"].map(cluster_to_type)
    for gene in ["SOX2", "NANOG", "GATA1", "KLF1"]:
        out[gene] = marker_df.loc[out["cell_id"], gene].to_numpy()
    return out.reset_index(drop=True)
```

## Code cell 9

```python
UMAP_CACHE = CACHE_DIR / "umap_plot_data.tsv.gz"
GEX_QC_CACHE = CACHE_DIR / "gex_qc_before_filter.tsv.gz"

umap_data: pd.DataFrame | None = None
adata_qc_counts = None

if UMAP_CACHE.exists() and not FORCE_RECOMPUTE_GEX:
    umap_data = pd.read_csv(UMAP_CACHE, sep="\t")
    print("Loaded cached UMAP data:", UMAP_CACHE, umap_data.shape)
elif DGE_FILTERED_DIR.exists():
    adata_raw = load_raw_gex(DGE_FILTERED_DIR)
    adata_qc_counts, gex_before_qc = qc_gex(adata_raw)
    gex_before_qc.to_csv(GEX_QC_CACHE, sep="\t", index=False, compression="gzip")
    umap_data = compute_umap_plot_data(adata_qc_counts)
    umap_data.to_csv(UMAP_CACHE, sep="\t", index=False, compression="gzip")
    print("Cached UMAP data:", UMAP_CACHE)
else:
    print("Raw GEX path is unavailable and no UMAP cache was found. UMAP panels will be skipped.")

if umap_data is not None:
    display(umap_data.head())
    display(umap_data.groupby(["sample", "cell_type_umap"], observed=True).size().rename("n_cells").reset_index())
```

### Text output (stream)

```text
Loading: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/splitpipe_out/WT/sub1/all-sample/DGE_filtered/count_matrix.mtx

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/functools.py:909: ImplicitModificationWarning: Transforming to str index.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/functools.py:909: ImplicitModificationWarning: Transforming to str index.
/tmp/ipykernel_12061/1217581716.py:38: UserWarning: Suffix used (-[0-9]+) to deduplicate index values may make index values difficult to interpret. There values with a similar suffixes in the index. Consider using a different delimiter by passing `join={delimiter}`. Example key collisions generated by the make_index_unique algorithm: ['SNORD115-1', 'SNORD115-2', 'SNORD115-3', 'SNORD115-4', 'SNORD115-5']

```

### Text output (stream)

```text
Cells before QC: 17,225; after QC: 16,818

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/scanpy/preprocessing/_simple.py:666: NumbaPerformanceWarning: [1m[1m'@' is faster on contiguous arrays, called on (Array(float64, 1, 'A', False, aligned=True), Array(float64, 2, 'C', False, aligned=True))[0m[0m
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/scanpy/preprocessing/_simple.py:666: NumbaPerformanceWarning: [1m[1m'@' is faster on contiguous arrays, called on (Array(float64, 1, 'A', False, aligned=True), Array(float64, 2, 'C', False, aligned=True))[0m[0m

```

### Text output (stream)

```text
Leiden cluster labels: {'0': 'K562', '1': 'hiPSC'}
Cached UMAP data: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/cache/umap_plot_data.tsv.gz

```

### Text output (text/plain)

```text
   bc_wells sample species  gene_count  tscp_count  mread_count  bc1_wind  \
0  25_01_18    CS1    hg38        1045        1256         1283        25   
1  25_01_38    CS1    hg38        1960        2731         2787        25   
2  25_01_72    CS1    hg38        1127        1538         1581        25   
3  25_02_67    CS1    hg38        1362        1710         1750        25   
4  25_02_68    CS1    hg38        1159        1574         1616        25   

   bc2_wind  bc3_wind bc1_well  ... pct_counts_mt leiden   cell_id      UMAP1  \
0         1        18       C1  ...      5.175159      1  25_01_18   9.044277   
1         1        38       C1  ...      5.565727      1  25_01_38   9.810464   
2         1        72       C1  ...      7.282185      0  25_01_72  -7.111281   
3         2        67       C1  ...      5.204678      1  25_02_67  10.734505   
4         2        68       C1  ...      6.289708      0  25_02_68  -6.918471   

      UMAP2  cell_type_umap SOX2    NANOG  GATA1  KLF1  
0  4.063761           hiPSC  0.0  0.00000    0.0   0.0  
1  0.536969           hiPSC  0.0  0.00000    0.0   0.0  
2 -2.183275            K562  0.0  0.00000    0.0   0.0  
3  1.336864           hiPSC  0.0  1.92395    0.0   0.0  
4 -1.963791            K562  0.0  0.00000    0.0   0.0  

[5 rows x 25 columns]
```

### Text output (text/plain)

```text
     sample cell_type_umap  n_cells
0       CS1           K562     2178
1       CS1          hiPSC     1800
2       DC1           K562     2321
3       DC1          hiPSC     1969
4      Dual           K562     2186
5      Dual          hiPSC     1789
6  no_guide           K562     2534
7  no_guide          hiPSC     2041
```

## Markdown cell 10

## 5. Process guide-count data with a low-memory cell-by-guide workflow

## Code cell 11

```python
def guide_family(guide_id: str) -> str:
    guide_id = str(guide_id)
    if guide_id.startswith("iPSC_"):
        return "iPSC"
    if guide_id.startswith("K562_"):
        return "K562"
    if re.search(r"(?i)non.?target|ntc|negative|scramble|control", guide_id):
        return "NTC"
    return "other"


def load_filtered_guide_long(
    path: Path,
    valid_cells: Iterable[str],
    chunksize: int = GUIDE_CHUNK_SIZE,
) -> pd.DataFrame:
    valid = set(map(str, valid_cells))
    kept = []
    rows_seen = 0
    rows_kept = 0

    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=["gene", "cell", "count"],
        chunksize=chunksize,
    ):
        rows_seen += len(chunk)
        chunk["cell_id"] = chunk["cell"].map(normalize_cell_id)
        sub = chunk.loc[chunk["cell_id"].isin(valid), ["cell_id", "gene", "count"]].copy()
        if sub.empty:
            continue
        sub["guide_id"] = sub["gene"].astype(str).str.replace(r"^gRNA_", "", regex=True)
        sub["count"] = pd.to_numeric(sub["count"], errors="coerce").fillna(0)
        sub = sub.loc[sub["count"] > 0, ["cell_id", "guide_id", "count"]]
        rows_kept += len(sub)
        kept.append(sub)

    if not kept:
        raise RuntimeError("No guide-count rows matched the QC-passing GEX cell barcodes.")

    long = pd.concat(kept, ignore_index=True)
    long = long.groupby(["cell_id", "guide_id"], as_index=False, observed=True)["count"].sum()
    long["guide_family"] = long["guide_id"].map(guide_family)
    print(f"Guide rows read: {rows_seen:,}; selected before aggregation: {rows_kept:,}; aggregated: {len(long):,}")
    return long


def compute_cell_guide_metrics(
    guide_long: pd.DataFrame,
    cell_annotation: pd.DataFrame,
    threshold: int = GUIDE_UMI_THRESHOLD,
) -> pd.DataFrame:
    base_cols = ["cell_id", "sample"]
    if "cell_type_umap" in cell_annotation.columns:
        base_cols.append("cell_type_umap")
    for col in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
        if col in cell_annotation.columns:
            base_cols.append(col)
    out = cell_annotation[base_cols].drop_duplicates("cell_id").set_index("cell_id")

    total = guide_long.groupby("cell_id", observed=True)["count"].sum().rename("total_gRNA")
    top = guide_long.groupby("cell_id", observed=True)["count"].max().rename("top_guide_umis")

    species_all = (
        guide_long.loc[guide_long["guide_family"].isin(["iPSC", "K562"])]
        .pivot_table(index="cell_id", columns="guide_family", values="count", aggfunc="sum", fill_value=0)
        .rename(columns={"iPSC": "iPSC_total", "K562": "K562_total"})
    )

    thresholded = guide_long.loc[guide_long["count"] >= threshold].copy()
    species_thr = (
        thresholded.loc[thresholded["guide_family"].isin(["iPSC", "K562"])]
        .pivot_table(index="cell_id", columns="guide_family", values="count", aggfunc="sum", fill_value=0)
        .rename(columns={"iPSC": f"iPSC_total_ge{threshold}", "K562": f"K562_total_ge{threshold}"})
    )
    total_thr = thresholded.groupby("cell_id", observed=True)["count"].sum().rename(f"total_gRNA_ge{threshold}")
    n_guides_thr = thresholded.groupby("cell_id", observed=True).size().rename(f"n_guides_ge{threshold}")

    out = out.join(total).join(top).join(species_all).join(species_thr).join(total_thr).join(n_guides_thr)
    numeric_cols = [
        "total_gRNA", "top_guide_umis", "iPSC_total", "K562_total",
        f"iPSC_total_ge{threshold}", f"K562_total_ge{threshold}",
        f"total_gRNA_ge{threshold}", f"n_guides_ge{threshold}",
    ]
    for col in numeric_cols:
        if col not in out.columns:
            out[col] = 0.0
        out[col] = out[col].fillna(0)

    if "cell_type_umap" in out.columns:
        out["expected_species_total_gRNA"] = np.where(
            out["cell_type_umap"].eq("K562"), out["K562_total"], out["iPSC_total"]
        )
        out["off_species_total_gRNA"] = np.where(
            out["cell_type_umap"].eq("K562"), out["iPSC_total"], out["K562_total"]
        )
        out[f"expected_species_total_ge{threshold}"] = np.where(
            out["cell_type_umap"].eq("K562"), out[f"K562_total_ge{threshold}"], out[f"iPSC_total_ge{threshold}"]
        )
        out[f"off_species_total_ge{threshold}"] = np.where(
            out["cell_type_umap"].eq("K562"), out[f"iPSC_total_ge{threshold}"], out[f"K562_total_ge{threshold}"]
        )
        denom = out["expected_species_total_gRNA"] + out["off_species_total_gRNA"]
        out["correct_species_guide_fraction"] = np.divide(
            out["expected_species_total_gRNA"], denom,
            out=np.full(len(out), np.nan), where=denom.to_numpy() > 0,
        )

    out = out.reset_index()
    out["sample"] = ordered_sample(out["sample"])
    return out
```

## Code cell 12

```python
CELL_METRICS_CACHE = CACHE_DIR / "cell_metrics.tsv.gz"
GUIDE_LONG_CACHE = CACHE_DIR / "filtered_cell_guide_long.tsv.gz"

cell_metrics: pd.DataFrame | None = None

if CELL_METRICS_CACHE.exists() and not FORCE_RECOMPUTE_GUIDES:
    cell_metrics = pd.read_csv(CELL_METRICS_CACHE, sep="\t")
    cell_metrics["sample"] = ordered_sample(cell_metrics["sample"])
    print("Loaded cached guide metrics:", CELL_METRICS_CACHE, cell_metrics.shape)

elif GUIDE_COUNT_TSV.exists() and umap_data is not None:
    if GUIDE_LONG_CACHE.exists() and not FORCE_RECOMPUTE_GUIDES:
        guide_long = pd.read_csv(GUIDE_LONG_CACHE, sep="\t")
        print("Loaded cached long guide table:", guide_long.shape)
    else:
        guide_long = load_filtered_guide_long(GUIDE_COUNT_TSV, umap_data["cell_id"])
        guide_long.to_csv(GUIDE_LONG_CACHE, sep="\t", index=False, compression="gzip")
        print("Cached long guide table:", GUIDE_LONG_CACHE)

    cell_metrics = compute_cell_guide_metrics(guide_long, umap_data)
    cell_metrics.to_csv(CELL_METRICS_CACHE, sep="\t", index=False, compression="gzip")
    print("Cached per-cell metrics:", CELL_METRICS_CACHE)

else:
    # Fallback: extract and load the uploaded/precomputed plot-input bundle.
    bundle = first_existing(BUNDLE_CANDIDATES)
    if bundle is None:
        raise FileNotFoundError(
            "Neither Sherlock raw data/caches nor k562_parsev3_figure_panel_plot_input_bundle.zip were found."
        )
    bundle_dir = CACHE_DIR / "uploaded_plot_input_bundle"
    bundle_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(bundle) as zf:
        zf.extractall(bundle_dir)
    cell_metrics = pd.read_csv(bundle_dir / "cell_metrics.tsv.gz", sep="\t")
    cell_metrics["sample"] = ordered_sample(cell_metrics["sample"])
    if "cell_type" in cell_metrics.columns and "cell_type_umap" not in cell_metrics.columns:
        cell_metrics["cell_type_umap"] = cell_metrics["cell_type"].replace({"unassigned": np.nan})
    print("Loaded attached plot-input bundle:", bundle, cell_metrics.shape)

# When UMAP labels are available, use them as the slide-reproduction cell type.
if umap_data is not None:
    label_map = umap_data.set_index("cell_id")["cell_type_umap"]
    cell_metrics["cell_type_slide"] = cell_metrics["cell_id"].map(label_map)
else:
    cell_metrics["cell_type_slide"] = cell_metrics.get("cell_type_umap", cell_metrics.get("cell_type"))

# Ensure species-threshold columns are available in fallback mode.
for family in ["iPSC", "K562"]:
    threshold_col = f"{family}_total_ge{GUIDE_UMI_THRESHOLD}"
    if threshold_col not in cell_metrics.columns:
        # The bundle stores all-count species totals; this is a fallback approximation only.
        cell_metrics[threshold_col] = cell_metrics[f"{family}_total"].where(
            cell_metrics[f"{family}_total"] >= GUIDE_UMI_THRESHOLD, 0
        )

cell_metrics["sample"] = ordered_sample(cell_metrics["sample"])
display(cell_metrics.head())
```

### Text output (stream)

```text
Guide rows read: 2,057,531; selected before aggregation: 860,384; aggregated: 860,384
Cached long guide table: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/cache/filtered_cell_guide_long.tsv.gz
Cached per-cell metrics: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/cache/cell_metrics.tsv.gz

```

### Text output (text/plain)

```text
    cell_id sample cell_type_umap  total_counts  n_genes_by_counts  \
0  25_01_18    CS1          hiPSC          1256               1045   
1  25_01_38    CS1          hiPSC          2731               1960   
2  25_01_72    CS1           K562          1538               1127   
3  25_02_67    CS1          hiPSC          1710               1362   
4  25_02_68    CS1           K562          1574               1159   

   pct_counts_mt  total_gRNA  top_guide_umis  K562_total  iPSC_total  \
0       5.175159       118.0            37.0        52.0        66.0   
1       5.565727       258.0            36.0        82.0       176.0   
2       7.282185       762.0           104.0       614.0       148.0   
3       5.204678       116.0            34.0        55.0        61.0   
4       6.289708       177.0            31.0       138.0        39.0   

   K562_total_ge3  iPSC_total_ge3  total_gRNA_ge3  n_guides_ge3  \
0             0.0            37.0            37.0           1.0   
1             0.0           125.0           125.0           7.0   
2           528.0           113.0           641.0          21.0   
3             0.0            34.0            34.0           1.0   
4            85.0             9.0            94.0           6.0   

   expected_species_total_gRNA  off_species_total_gRNA  \
0                         66.0                    52.0   
1                        176.0                    82.0   
2                        614.0                   148.0   
3                         61.0                    55.0   
4                        138.0                    39.0   

   expected_species_total_ge3  off_species_total_ge3  \
0                        37.0                    0.0   
1                       125.0                    0.0   
2                       528.0                  113.0   
3                        34.0                    0.0   
4                        85.0                    9.0   

   correct_species_guide_fraction cell_type_slide  
0                        0.559322           hiPSC  
1                        0.682171           hiPSC  
2                        0.805774            K562  
3                        0.525862           hiPSC  
4                        0.779661            K562  
```

## Markdown cell 13

## 6. Load transcriptome saturation data

## Code cell 14

```python
def read_saturation_file(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, header=None)
    raw = raw.iloc[1:, :2].copy()
    raw.columns = ["read_count", "transcript_count"]
    raw["read_count"] = pd.to_numeric(raw["read_count"], errors="coerce")
    raw["transcript_count"] = pd.to_numeric(raw["transcript_count"], errors="coerce")
    return raw.dropna().sort_values("read_count")


def load_saturation_files() -> pd.DataFrame:
    rows = []
    for sample, path in SATURATION_FILES.items():
        if path.exists():
            x = read_saturation_file(path)
            x.insert(0, "sample", sample)
            rows.append(x)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


SATURATION_CACHE = CACHE_DIR / "transcriptome_saturation.tsv"
if SATURATION_CACHE.exists():
    saturation = pd.read_csv(SATURATION_CACHE, sep="\t")
else:
    saturation = load_saturation_files()
    if saturation.empty:
        bundle_dir = CACHE_DIR / "uploaded_plot_input_bundle"
        if (bundle_dir / "transcriptome_saturation.tsv").exists():
            saturation = pd.read_csv(bundle_dir / "transcriptome_saturation.tsv", sep="\t")
    if not saturation.empty:
        saturation.to_csv(SATURATION_CACHE, sep="\t", index=False)

if not saturation.empty:
    saturation["sample"] = ordered_sample(saturation["sample"])
    display(saturation.head())
else:
    print("No transcriptome saturation files were found.")
```

### Text output (text/plain)

```text
  sample  read_count  transcript_count
0    CS1         0.0                 0
1    CS1       160.0                42
2    CS1       372.0                98
3    CS1       585.0               155
4    CS1       904.0               237
```

## Markdown cell 15

## 7. Reproducibility checks

The values below should approximately reproduce the slide-level checks:

- 17,225 cells before QC and 16,818 after QC;
- median transcript count near 1,432 per QC-passing cell;
- mean total guide UMIs per cell near 270, 318, 352 and 3 for CS1, DC1, Dual and no-guide, respectively.

## Code cell 16

```python
summary = (
    cell_metrics.groupby("sample", observed=True)
    .agg(
        n_cells=("cell_id", "nunique"),
        median_transcripts=("total_counts", "median"),
        median_genes=("n_genes_by_counts", "median"),
        mean_total_guide_umis=("total_gRNA", "mean"),
        median_total_guide_umis=("total_gRNA", "median"),
        median_top_guide_umis=("top_guide_umis", "median"),
    )
    .reset_index()
)
display(summary)

# Pairwise primer tests on raw total-gRNA counts.
rows = []
for i, a in enumerate(PRIMER_SAMPLES):
    for b in PRIMER_SAMPLES[i + 1:]:
        xa = cell_metrics.loc[cell_metrics["sample"].astype(str).eq(a), "total_gRNA"].to_numpy()
        xb = cell_metrics.loc[cell_metrics["sample"].astype(str).eq(b), "total_gRNA"].to_numpy()
        p = mannwhitneyu(xa, xb, alternative="two-sided", method="asymptotic").pvalue
        rows.append({"sample_a": a, "sample_b": b, "pvalue": p})
pairwise = pd.DataFrame(rows)
pairwise["pvalue_holm"] = holm_adjust(pairwise["pvalue"])
display(pairwise)
```

### Text output (text/plain)

```text
     sample  n_cells  median_transcripts  median_genes  mean_total_guide_umis  \
0       CS1     3978              1358.0        1069.0             269.602313   
1       DC1     4290              1465.0        1147.5             317.973193   
2      Dual     3975              1398.0        1098.0             352.328302   
3  no_guide     4575              1491.0        1165.0               2.823388   

   median_total_guide_umis  median_top_guide_umis  
0                    202.0                   45.0  
1                    207.0                   69.0  
2                    248.0                   70.0  
3                      2.0                    1.0  
```

### Text output (text/plain)

```text
  sample_a sample_b        pvalue   pvalue_holm
0      CS1      DC1  8.653515e-01  8.653515e-01
1      CS1     Dual  1.398706e-18  4.196117e-18
2      DC1     Dual  3.768599e-16  7.537198e-16
```

## Markdown cell 17

## 8. Slide 8 — transcriptome QC distributions

## Code cell 18

```python
qc_long = cell_metrics[["cell_id", "total_counts", "n_genes_by_counts", "pct_counts_mt"]].melt(
    id_vars="cell_id", var_name="metric", value_name="value"
)
qc_long["metric_label"] = qc_long["metric"].map({
    "total_counts": "Transcript counts",
    "n_genes_by_counts": "Gene counts",
    "pct_counts_mt": "% mitochondrial",
})
qc_long["metric_label"] = pd.Categorical(
    qc_long["metric_label"],
    categories=["Transcript counts", "Gene counts", "% mitochondrial"],
    ordered=True,
)

p_qc = (
    ggplot(qc_long, aes(x="metric_label", y="value"))
    + geom_violin(fill="#2C7FB8", alpha=0.8, trim=True, width=0.85)
    + geom_jitter(width=0.13, height=0, size=0.16, alpha=0.20, color="black", raster=True)
    + facet_wrap("~metric_label", scales="free_y", nrow=1)
    + labs(x=None, y="Value", title="Gene-expression quality control")
    + theme_publication(10)
    + theme(
        axis_text_x=element_blank(),
        axis_ticks_major_x=element_blank(),
        axis_title_x=element_blank(),
        strip_text=element_text(size=11, weight="normal"),
    )
)
p_qc
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f9377d15dd0>
```

## Code cell 19

```python
save_plot(p_qc, "slide08_gene_expression_QC", width=10.5, height=3.8)
```

### Text output (stream)

```text
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/slide08_gene_expression_QC in ('png', 'svg')

```

## Markdown cell 20

## 9. Slide 9 — transcriptome yield and saturation across primer conditions

## Code cell 21

```python
yield_df = cell_metrics[["sample", "total_counts"]].dropna().copy()
yield_df["sample"] = ordered_sample(yield_df["sample"])
counts = yield_df.groupby("sample", observed=True).size()
yield_df["sample_label"] = yield_df["sample"].astype(str).map(
    {s: f"{s} (n={int(counts.get(s, 0)):,})" for s in SAMPLE_ORDER}
)
yield_order = [f"{s} (n={int(counts.get(s, 0)):,})" for s in SAMPLE_ORDER]
yield_df["sample_label"] = pd.Categorical(yield_df["sample_label"], categories=yield_order, ordered=True)

p_yield = (
    ggplot(yield_df, aes(x="sample_label", y="total_counts", fill="sample"))
    + geom_violin(trim=True, width=0.9, alpha=0.75)
    + geom_boxplot(width=0.12, outlier_alpha=0, fill="white", color="black")
    + coord_flip()
    + scale_fill_manual(values=SAMPLE_COLORS, guide=None)
    + labs(x=None, y="Total UMI counts per cell", title="Transcriptome yield")
    + theme_publication(10)
    + theme(axis_title_y=element_blank())
)

if not saturation.empty:
    p_saturation = (
        ggplot(saturation, aes(x="read_count", y="transcript_count", color="sample", group="sample"))
        + geom_line(size=0.8)
        + geom_point(size=1.5)
        + scale_color_manual(values=SAMPLE_COLORS)
        + labs(
            x="Mean reads per cell",
            y="Median transcript count",
            color="Primer type",
            title="Sequencing saturation",
        )
        + theme_publication(10)
        + theme(legend_position=(0.18, 0.76))
    )
    p_transcriptome = p_yield | p_saturation
else:
    p_transcriptome = p_yield

p_transcriptome
```

### Text output (text/plain)

```text
<plotnine.composition._beside.Beside object at 0x7f93f4d46750>
```

## Code cell 22

```python
save_plot(p_transcriptome, "slide09_transcriptome_yield_and_saturation", width=11, height=4.2)
```

### Text output (stream)

```text
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/slide09_transcriptome_yield_and_saturation in ('png', 'svg')

```

## Markdown cell 23

## 10. Slide 12 — guide UMI capture across primer conditions

## Code cell 24

```python
guide_df = cell_metrics[["sample", "total_gRNA"]].copy()
guide_df["sample"] = ordered_sample(guide_df["sample"])
guide_df["total_gRNA_plot"] = guide_df["total_gRNA"] + 1

ann = (
    guide_df.groupby("sample", observed=True)
    .agg(mean_raw=("total_gRNA", "mean"), y=("total_gRNA_plot", lambda x: x.quantile(0.995) * 1.35))
    .reset_index()
)
ann["label"] = ann["mean_raw"].round().astype(int).astype(str)
ann["mean_plot"] = ann["mean_raw"] + 1
ann["sample"] = ordered_sample(ann["sample"])

p_guide_primer = (
    ggplot(guide_df, aes(x="sample", y="total_gRNA_plot", fill="sample", color="sample"))
    + geom_jitter(width=0.24, height=0, size=0.55, alpha=0.35, raster=True)
    + geom_boxplot(width=0.48, alpha=0.68, outlier_alpha=0, color="#444444")
    + geom_point(data=ann, mapping=aes(x="sample", y="mean_plot"), inherit_aes=False,
                 color="black", fill="black", shape="D", size=2.2)
    + geom_text(data=ann, mapping=aes(x="sample", y="y", label="label"), inherit_aes=False,
                size=8.5, va="bottom")
    + scale_y_log10(labels=comma_format())
    + scale_fill_manual(values=SAMPLE_COLORS, guide=None)
    + scale_color_manual(values=SAMPLE_COLORS, guide=None)
    + labs(
        x="Primer type",
        y="Total gRNA UMI counts per cell (log scale)",
        title="gRNA UMI counts across RT-primer conditions",
        subtitle="Black diamonds and labels show arithmetic means on the raw-count scale",
    )
    + theme_publication(11)
)
p_guide_primer
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f93872348d0>
```

## Code cell 25

```python
save_plot(p_guide_primer, "slide12_total_gRNA_by_primer", width=7.5, height=5.2)
```

### Text output (stream)

```text
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/slide12_total_gRNA_by_primer in ('png', 'svg')

```

## Markdown cell 26

## 11. Slide 13 — iPSC-guide versus K562-guide mixing

## Code cell 27

```python
i_col = f"iPSC_total_ge{GUIDE_UMI_THRESHOLD}"
k_col = f"K562_total_ge{GUIDE_UMI_THRESHOLD}"

mix_df = cell_metrics[["cell_id", "sample", i_col, k_col]].copy()
mix_df = mix_df.loc[(mix_df[i_col] >= GUIDE_UMI_THRESHOLD) | (mix_df[k_col] >= GUIDE_UMI_THRESHOLD)].copy()
mix_df["iPSC_plot"] = mix_df[i_col] + 1
mix_df["K562_plot"] = mix_df[k_col] + 1

p_mix = (
    ggplot(mix_df, aes(x="iPSC_plot", y="K562_plot"))
    + stat_bin_2d(bins=60)
    + geom_abline(slope=1, intercept=0, linetype="dashed", color="white", alpha=0.7)
    + scale_x_log10(labels=comma_format())
    + scale_y_log10(labels=comma_format())
    + scale_fill_cmap(cmap_name="viridis", trans="log10", name="Cells per bin")
    + coord_equal()
    + labs(
        x="Total iPSC guide UMIs + 1 (log scale)",
        y="Total K562 guide UMIs + 1 (log scale)",
        title=f"iPSC versus K562 guide UMIs per cell (n={len(mix_df):,})",
        subtitle=f"Guide counts were retained when the per-guide UMI count was ≥ {GUIDE_UMI_THRESHOLD}",
    )
    + theme_publication(10)
)

p_mix_x = (
    ggplot(mix_df, aes(x="iPSC_plot"))
    + geom_histogram(bins=60, fill="#4E79A7", alpha=0.8)
    + scale_x_log10(labels=comma_format())
    + labs(x="iPSC guide UMIs + 1", y="Cells", title="iPSC-guide marginal distribution")
    + theme_publication(9)
)

p_mix_y = (
    ggplot(mix_df, aes(x="K562_plot"))
    + geom_histogram(bins=60, fill="#76B7B2", alpha=0.8)
    + scale_x_log10(labels=comma_format())
    + coord_flip()
    + labs(x="K562 guide UMIs + 1", y="Cells", title="K562-guide marginal distribution")
    + theme_publication(9)
)

p_mix
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f9377b7ee10>
```

## Code cell 28

```python
save_plot(p_mix, "slide13_species_mixing_2D_histogram", width=6.2, height=5.4)
save_plot(p_mix_x, "slide13_species_mixing_iPSC_marginal", width=5.2, height=3.2)
save_plot(p_mix_y, "slide13_species_mixing_K562_marginal", width=5.2, height=3.2)
```

### Text output (stream)

```text
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/slide13_species_mixing_2D_histogram in ('png', 'svg')
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/slide13_species_mixing_iPSC_marginal in ('png', 'svg')
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/slide13_species_mixing_K562_marginal in ('png', 'svg')

```

## Markdown cell 29

## 12. Slide 14 — UMAP cell types and marker genes

## Code cell 30

```python
if umap_data is None:
    print("UMAP panel skipped because raw GEX data or a UMAP cache was unavailable.")
    p_umap_celltype = None
    p_umap_markers = None
else:
    cluster_centers = (
        umap_data.groupby("cell_type_umap", observed=True)[["UMAP1", "UMAP2"]]
        .median()
        .reset_index()
    )

    p_umap_celltype = (
        ggplot(umap_data, aes(x="UMAP1", y="UMAP2", color="cell_type_umap"))
        + geom_point(size=0.40, alpha=0.70, raster=True)
        + geom_text(
            data=cluster_centers,
            mapping=aes(x="UMAP1", y="UMAP2", label="cell_type_umap"),
            inherit_aes=False,
            color="black",
            size=11,
            fontweight="bold",
        )
        + scale_color_manual(values=CELLTYPE_COLORS, name="Cell type")
        + coord_equal()
        + labs(title="K562 and hiPSC cell-state separation", x="UMAP1", y="UMAP2")
        + theme_publication(10)
    )

    marker_long = umap_data[["cell_id", "UMAP1", "UMAP2", "SOX2", "NANOG", "GATA1", "KLF1"]].melt(
        id_vars=["cell_id", "UMAP1", "UMAP2"], var_name="marker", value_name="expression"
    )
    marker_long["marker"] = pd.Categorical(
        marker_long["marker"], categories=["SOX2", "NANOG", "GATA1", "KLF1"], ordered=True
    )
    # Order points within each marker so highly expressing cells are drawn last.
    marker_long = marker_long.sort_values(["marker", "expression"], na_position="first")

    p_umap_markers = (
        ggplot(marker_long, aes(x="UMAP1", y="UMAP2", color="expression"))
        + geom_point(size=0.32, alpha=0.85, raster=True)
        + facet_wrap("~marker", ncol=2)
        + scale_color_cmap(cmap_name="Blues", name="log-normalized\nexpression")
        + coord_equal()
        + labs(x="UMAP1", y="UMAP2", title="Lineage-marker expression")
        + theme_publication(9)
    )

    p_umap_slide = p_umap_celltype | p_umap_markers
    display(p_umap_slide)
```

### Text output (text/plain)

```text
<plotnine.composition._beside.Beside object at 0x7f93f31c6e50>
```

## Code cell 31

```python
if umap_data is not None:
    save_plot(p_umap_celltype, "slide14_umap_cell_type", width=5.2, height=4.4)
    save_plot(p_umap_markers, "slide14_umap_marker_genes", width=7.2, height=6.0)
    save_plot(p_umap_slide, "slide14_umap_combined", width=12.2, height=6.0)
```

### Text output (stream)

```text
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/slide14_umap_cell_type in ('png', 'svg')
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/slide14_umap_marker_genes in ('png', 'svg')
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/slide14_umap_combined in ('png', 'svg')

```

## Markdown cell 32

## 13. Slide 15 — species guide totals on the UMAP

## Code cell 33

```python
if umap_data is None:
    print("Guide-on-UMAP panel skipped because UMAP coordinates were unavailable.")
    p_guide_umap = None
else:
    guide_umap = umap_data[["cell_id", "UMAP1", "UMAP2"]].merge(
        cell_metrics[["cell_id", i_col, k_col]], on="cell_id", how="left"
    )
    guide_umap[[i_col, k_col]] = guide_umap[[i_col, k_col]].fillna(0)
    guide_umap["iPSC guide UMIs"] = np.log1p(guide_umap[i_col])
    guide_umap["K562 guide UMIs"] = np.log1p(guide_umap[k_col])
    guide_umap_long = guide_umap.melt(
        id_vars=["cell_id", "UMAP1", "UMAP2"],
        value_vars=["iPSC guide UMIs", "K562 guide UMIs"],
        var_name="guide_source", value_name="log1p_UMIs",
    ).sort_values(["guide_source", "log1p_UMIs"])

    p_guide_umap = (
        ggplot(guide_umap_long, aes(x="UMAP1", y="UMAP2", color="log1p_UMIs"))
        + geom_point(size=0.34, alpha=0.85, raster=True)
        + facet_wrap("~guide_source", nrow=1)
        + scale_color_cmap(cmap_name="Blues", name="log(1 + guide UMIs)")
        + coord_equal()
        + labs(
            x="UMAP1", y="UMAP2",
            title="Species-specific guide UMIs on the transcriptome UMAP",
            subtitle=f"Per-guide counts below {GUIDE_UMI_THRESHOLD} UMIs were excluded before summation",
        )
        + theme_publication(10)
    )
    p_guide_umap
```

## Code cell 34

```python
if p_guide_umap is not None:
    save_plot(p_guide_umap, "slide15_species_guide_UMIs_on_UMAP", width=10.2, height=4.7)
```

### Text output (stream)

```text
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/slide15_species_guide_UMIs_on_UMAP in ('png', 'svg')

```

## Markdown cell 35

## 14. Slide 16 — guide capture by primer condition and inferred cell type

## Code cell 36

```python
celltype_df = cell_metrics.loc[
    cell_metrics["cell_type_slide"].isin(["K562", "hiPSC"]),
    ["sample", "cell_type_slide", "total_gRNA"],
].copy()
celltype_df["sample"] = ordered_sample(celltype_df["sample"])
celltype_df["cell_type_slide"] = pd.Categorical(
    celltype_df["cell_type_slide"], categories=["K562", "hiPSC"], ordered=True
)
celltype_df["total_gRNA_plot"] = celltype_df["total_gRNA"] + 1

p_guide_celltype = (
    ggplot(
        celltype_df,
        aes(x="sample", y="total_gRNA_plot", fill="cell_type_slide", color="cell_type_slide"),
    )
    + geom_jitter(
        position=position_jitterdodge(jitter_width=0.22, dodge_width=0.72),
        size=0.50, alpha=0.32, raster=True,
    )
    + geom_boxplot(
        position=position_dodge(width=0.72), width=0.58,
        alpha=0.70, outlier_alpha=0, color="#444444",
    )
    + scale_y_log10(labels=comma_format())
    + scale_fill_manual(values=CELLTYPE_COLORS, name="Cell type")
    + scale_color_manual(values=CELLTYPE_COLORS, name="Cell type")
    + labs(
        x="Primer type",
        y="Total gRNA UMI counts per cell (log scale)",
        title="gRNA UMI counts for each cell type",
    )
    + theme_publication(11)
)
p_guide_celltype
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f93e6f421d0>
```

## Code cell 37

```python
save_plot(p_guide_celltype, "slide16_total_gRNA_by_primer_and_cell_type", width=8.4, height=5.2)
```

### Text output (stream)

```text
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/slide16_total_gRNA_by_primer_and_cell_type in ('png', 'svg')

```

## Markdown cell 38

## 15. Optional paper-panel versions

The attached export script and publication-ready plot bundle use **top-guide UMIs per cell** for the primary guide-capture panels. The following panels reproduce that paper-oriented metric while retaining the same processed dataset.

## Code cell 39

```python
def make_top_guide_panel(data: pd.DataFrame, title: str, include_no_guide: bool = False):
    samples = SAMPLE_ORDER if include_no_guide else PRIMER_SAMPLES
    x = data.loc[data["sample"].astype(str).isin(samples)].copy()
    x["sample"] = pd.Categorical(x["sample"].astype(str), categories=samples, ordered=True)
    x["top_plot"] = x["top_guide_umis"] + 1
    counts = x.groupby("sample", observed=True).size()
    x["sample_label"] = x["sample"].astype(str).map(
        {s: f"{s}\n(n={int(counts.get(s, 0)):,})" for s in samples}
    )
    x["sample_label"] = pd.Categorical(
        x["sample_label"],
        categories=[f"{s}\n(n={int(counts.get(s, 0)):,})" for s in samples],
        ordered=True,
    )
    ann = (
        x.groupby(["sample", "sample_label"], observed=True)
        .agg(median=("top_guide_umis", "median"), y=("top_plot", lambda z: z.quantile(0.995) * 1.25))
        .reset_index()
    )
    ann["label"] = ann["median"].round().astype(int).astype(str)
    return (
        ggplot(x, aes(x="sample_label", y="top_plot", fill="sample"))
        + geom_violin(trim=True, width=0.90, alpha=0.88)
        + geom_boxplot(width=0.16, fill="white", color="black", outlier_alpha=0)
        + geom_text(data=ann, mapping=aes(x="sample_label", y="y", label="label"),
                    inherit_aes=False, size=8)
        + scale_y_log10(labels=comma_format())
        + scale_fill_manual(values=SAMPLE_COLORS, guide=None)
        + labs(x=None, y="Top guide UMIs per cell", title=title, subtitle="Labels show medians")
        + theme_publication(10)
    )

paper_type_col = "cell_type" if "cell_type" in cell_metrics.columns else "cell_type_slide"

p_k562_top = make_top_guide_panel(
    cell_metrics.loc[cell_metrics[paper_type_col].eq("K562")],
    "Guide capture in K562 cells",
)
p_hipsc_top = make_top_guide_panel(
    cell_metrics.loc[cell_metrics[paper_type_col].eq("hiPSC")],
    "Guide capture in hiPSC cells",
)
p_top_pair = p_k562_top | p_hipsc_top
p_top_pair
```

### Text output (text/plain)

```text
<plotnine.composition._beside.Beside object at 0x7f93e540d2d0>
```

## Code cell 40

```python
save_plot(p_top_pair, "paper_top_guide_capture_K562_hiPSC", width=9.4, height=4.6)
```

### Text output (stream)

```text
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/paper_top_guide_capture_K562_hiPSC in ('png', 'svg')

```

## Code cell 41

```python
# Expected-species guide specificity and expected-species UMIs.
specificity = cell_metrics.loc[
    cell_metrics[paper_type_col].isin(["K562", "hiPSC"])
    & cell_metrics["sample"].astype(str).isin(PRIMER_SAMPLES)
].copy()
specificity["plot_cell_type"] = specificity[paper_type_col]
specificity["sample"] = pd.Categorical(specificity["sample"].astype(str), categories=PRIMER_SAMPLES, ordered=True)

if "correct_species_guide_fraction" not in specificity.columns:
    expected = np.where(
        specificity["plot_cell_type"].eq("K562"), specificity["K562_total"], specificity["iPSC_total"]
    )
    off = np.where(
        specificity["plot_cell_type"].eq("K562"), specificity["iPSC_total"], specificity["K562_total"]
    )
    denom = expected + off
    specificity["correct_species_guide_fraction"] = np.divide(
        expected, denom, out=np.full(len(specificity), np.nan), where=denom > 0
    )
    specificity["expected_species_total_gRNA"] = expected

p_specificity = (
    ggplot(specificity, aes(x="sample", y="correct_species_guide_fraction", fill="sample"))
    + geom_violin(trim=True, alpha=0.85)
    + geom_boxplot(width=0.15, fill="white", color="black", outlier_alpha=0)
    + facet_wrap("~plot_cell_type", nrow=1)
    + scale_fill_manual(values=SAMPLE_COLORS, guide=None)
    + scale_y_continuous(limits=(0, 1.02), labels=percent_format())
    + labs(x=None, y="Expected-species guide fraction", title="Expected-species guide specificity")
    + theme_publication(10)
)

expected_df = specificity.copy()
expected_df["expected_plot"] = expected_df["expected_species_total_gRNA"] + 1
p_expected = (
    ggplot(expected_df, aes(x="sample", y="expected_plot", fill="sample"))
    + geom_violin(trim=True, alpha=0.85)
    + geom_boxplot(width=0.15, fill="white", color="black", outlier_alpha=0)
    + facet_wrap("~plot_cell_type", nrow=1)
    + scale_fill_manual(values=SAMPLE_COLORS, guide=None)
    + scale_y_log10(labels=comma_format())
    + labs(x=None, y="Expected-species guide UMIs per cell", title="Expected-species guide UMIs")
    + theme_publication(10)
)

p_specificity | p_expected
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/layer.py:293: PlotnineWarning: stat_ydensity : Removed 7 rows containing non-finite values.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/layer.py:293: PlotnineWarning: stat_boxplot : Removed 7 rows containing non-finite values.

```

### Text output (text/plain)

```text
<plotnine.composition._beside.Beside object at 0x7f93e52a4090>
```

## Code cell 42

```python
save_plot(p_specificity, "paper_expected_species_specificity", width=7.4, height=4.2)
save_plot(p_expected, "paper_expected_species_guide_UMIs", width=7.4, height=4.2)
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/layer.py:293: PlotnineWarning: stat_ydensity : Removed 7 rows containing non-finite values.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/layer.py:293: PlotnineWarning: stat_boxplot : Removed 7 rows containing non-finite values.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/layer.py:293: PlotnineWarning: stat_ydensity : Removed 7 rows containing non-finite values.
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/plotnine/layer.py:293: PlotnineWarning: stat_boxplot : Removed 7 rows containing non-finite values.

```

### Text output (stream)

```text
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/paper_expected_species_specificity in ('png', 'svg')
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/paper_expected_species_guide_UMIs in ('png', 'svg')

```

## Markdown cell 43

## 16. Guide-detection curves across UMI thresholds

## Code cell 44

```python
thresholds = [1, 3, 5, 10, 20, 50]
rows = []
for cell_type in ["K562", "hiPSC"]:
    for sample in PRIMER_SAMPLES:
        sub = cell_metrics.loc[
            cell_metrics[paper_type_col].eq(cell_type)
            & cell_metrics["sample"].astype(str).eq(sample)
        ]
        for threshold in thresholds:
            rows.append({
                "cell_type": cell_type,
                "sample": sample,
                "threshold": threshold,
                "percent_detected": 100 * (sub["top_guide_umis"] >= threshold).mean() if len(sub) else np.nan,
            })
detection = pd.DataFrame(rows)
detection["sample"] = pd.Categorical(detection["sample"], categories=PRIMER_SAMPLES, ordered=True)

p_detection = (
    ggplot(detection, aes(x="threshold", y="percent_detected", color="sample", group="sample"))
    + geom_line(size=0.9)
    + geom_point(size=2.0)
    + facet_wrap("~cell_type", nrow=1)
    + scale_color_manual(values=SAMPLE_COLORS, name="Primer type")
    + scale_x_continuous(breaks=thresholds)
    + scale_y_continuous(limits=(0, 102), labels=lambda xs: [f"{x:.0f}" for x in xs])
    + labs(
        x="Top-guide UMI threshold",
        y="Cells detected (%)",
        title="Guide detection across UMI thresholds",
    )
    + theme_publication(10)
)
p_detection
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f93e71a46d0>
```

## Code cell 45

```python
save_plot(p_detection, "paper_guide_detection_threshold_curves", width=7.4, height=4.2)
```

### Text output (stream)

```text
Saved: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures/paper_guide_detection_threshold_curves in ('png', 'svg')

```

## Markdown cell 46

## 17. Output manifest

Figures are saved as PNG and editable-text SVG by default. Add `"pdf"` to `EXPORT_FORMATS` to create PDF copies as well. The table below lists the generated files.

## Code cell 47

```python
manifest = pd.DataFrame({
    "file": sorted(str(p.relative_to(FIGURE_DIR)) for p in FIGURE_DIR.glob("*") if p.is_file())
})
display(manifest)

manifest.to_csv(FIGURE_DIR / "figure_manifest.tsv", sep="\t", index=False)
print("Figure directory:", FIGURE_DIR)
```

### Text output (text/plain)

```text
                                              file
0            paper_expected_species_guide_UMIs.png
1            paper_expected_species_guide_UMIs.svg
2           paper_expected_species_specificity.png
3           paper_expected_species_specificity.svg
4       paper_guide_detection_threshold_curves.png
5       paper_guide_detection_threshold_curves.svg
6           paper_top_guide_capture_K562_hiPSC.png
7           paper_top_guide_capture_K562_hiPSC.svg
8                   slide08_gene_expression_QC.png
9                   slide08_gene_expression_QC.svg
10  slide09_transcriptome_yield_and_saturation.png
11  slide09_transcriptome_yield_and_saturation.svg
12                slide12_total_gRNA_by_primer.png
13                slide12_total_gRNA_by_primer.svg
14         slide13_species_mixing_2D_histogram.png
15         slide13_species_mixing_2D_histogram.svg
16        slide13_species_mixing_K562_marginal.png
17        slide13_species_mixing_K562_marginal.svg
18        slide13_species_mixing_iPSC_marginal.png
19        slide13_species_mixing_iPSC_marginal.svg
20                      slide14_umap_cell_type.png
21                      slide14_umap_cell_type.svg
22                       slide14_umap_combined.png
23                       slide14_umap_combined.svg
24                   slide14_umap_marker_genes.png
25                   slide14_umap_marker_genes.svg
26          slide15_species_guide_UMIs_on_UMAP.png
27          slide15_species_guide_UMIs_on_UMAP.svg
28  slide16_total_gRNA_by_primer_and_cell_type.png
29  slide16_total_gRNA_by_primer_and_cell_type.svg
```

### Text output (stream)

```text
Figure directory: /oak/stanford/groups/engreitz/Users/tri/Perturb-Seq/analysis/250730_eTN10_Parsev3_mixediPSC_K562/plotnine_reproduction/figures

```

