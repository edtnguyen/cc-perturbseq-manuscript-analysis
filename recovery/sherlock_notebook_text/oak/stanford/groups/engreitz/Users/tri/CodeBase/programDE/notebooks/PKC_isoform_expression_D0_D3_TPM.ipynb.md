# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/PKC_isoform_expression_D0_D3_TPM.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/PKC_isoform_expression_D0_D3_TPM.ipynb`
- SHA-256: `6e256cc4a97fb9f8d14043ee9f63182affd9ea2b08124a7db5db218d9b26b85d`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0


# PKC isoform expression across D0-D3 — plotnine version

This notebook calculates day-level pseudobulk expression for the canonical human PKC isoforms:

- conventional PKCs: `PRKCA`, `PRKCB`, `PRKCG`;
- novel PKCs: `PRKCD`, `PRKCE`, `PRKCH`, `PRKCQ`;
- atypical PKCs: `PRKCI`, `PRKCZ`.

The primary output is **pseudobulk TPM** calculated from raw UMI counts and gene lengths. The notebook also reports UMI counts and pseudobulk counts per million (CPM).

The newly supplied program-annotation notebook does not itself contain an AnnData path. The integrated AnnData path below is the path consistently used in the previously supplied CC-Perturb-seq UMAP and Hippo notebooks.


## Code cell 1

```python

from pathlib import Path
from collections import OrderedDict, defaultdict
import gzip
import re
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import matplotlib as mpl
import matplotlib.pyplot as plt
from IPython.display import display, Markdown
from plotnine import (
    aes,
    element_blank,
    element_rect,
    element_text,
    facet_grid,
    geom_line,
    geom_point,
    geom_text,
    geom_tile,
    ggplot,
    labs,
    scale_color_brewer,
    scale_fill_cmap,
    scale_x_discrete,
    theme,
    theme_classic,
    theme_minimal,
)

warnings.filterwarnings("ignore")
pd.set_option("display.max_columns", 200)
pd.set_option("display.max_rows", 200)

class cfg:
    ADATA_H5AD = Path(
        "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
        "2.preprocessing_tonys_pipeline/3.integrated_days/"
        "3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad"
    )

    DAY_OBS_COL = "day"
    OUTPUT_DIR = Path("PKC_isoform_expression_D0_D3_outputs")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # True TPM requires gene lengths. This notebook uses the GENCODE v43
    # annotation GTF supplied on Sherlock when a suitable length column is absent.
    GENCODE_V43_GTF = Path("/oak/stanford/groups/engreitz/Users/tri/tmp/gencode.v43.annotation.gtf.gz")

    COUNT_LAYER_CANDIDATES = ["counts", "raw_counts", "count"]
    GENE_SYMBOL_COL_CANDIDATES = [
        "gene_symbol", "gene_name", "symbol", "feature_name"
    ]
    GENE_ID_COL_CANDIDATES = [
        "gene_id", "ensembl_id", "feature_id"
    ]
    GENE_LENGTH_COL_CANDIDATES = [
        "gene_length_bp", "gene_length", "feature_length", "length"
    ]

    CELL_CHUNK_SIZE = 2000

    PKC_GROUPS = OrderedDict(
        [
            ("Conventional PKC", ["PRKCA", "PRKCB", "PRKCG"]),
            ("Novel PKC", ["PRKCD", "PRKCE", "PRKCH", "PRKCQ"]),
            ("Atypical PKC", ["PRKCI", "PRKCZ"]),
        ]
    )

    # These are protein kinase D genes, not PKC isoforms. Set True only
    # when a separate related-kinase comparison is desired.
    INCLUDE_PKD_FAMILY = False
    PKD_GENES = ["PRKD1", "PRKD2", "PRKD3"]

PKC_GENES = [
    gene
    for genes in cfg.PKC_GROUPS.values()
    for gene in genes
]
if cfg.INCLUDE_PKD_FAMILY:
    PKC_GENES += cfg.PKD_GENES

mpl.rcParams.update(
    {
        # Preserve all plotnine labels as native text in SVG and editable
        # TrueType text in PDF for Adobe Illustrator.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "text.usetex": False,
        "font.family": "Liberation Sans",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)

print("AnnData:", cfg.ADATA_H5AD)
print("Genes:", PKC_GENES)

```

### Text output (stream)

```text
AnnData: /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad
Genes: ['PRKCA', 'PRKCB', 'PRKCG', 'PRKCD', 'PRKCE', 'PRKCH', 'PRKCQ', 'PRKCI', 'PRKCZ']

```

## Markdown cell 2


## 1. Load the integrated AnnData and validate the day labels


## Code cell 3

```python

if not cfg.ADATA_H5AD.exists():
    raise FileNotFoundError(
        f"AnnData not found: {cfg.ADATA_H5AD}\n"
        "Run this notebook on Sherlock with /oak mounted, or update cfg.ADATA_H5AD."
    )

adata = sc.read_h5ad(cfg.ADATA_H5AD, backed="r")
print(adata)
print(f"Cells: {adata.n_obs:,}")
print(f"Features: {adata.n_vars:,}")

if cfg.DAY_OBS_COL not in adata.obs.columns:
    raise KeyError(
        f"{cfg.DAY_OBS_COL!r} is not in adata.obs. "
        f"Available columns include: {list(adata.obs.columns[:50])}"
    )

def normalize_day(value):
    text = str(value).strip().upper().replace("DAY", "D")
    match = re.search(r"([0-3])", text)
    return f"D{match.group(1)}" if match else None

day_labels = adata.obs[cfg.DAY_OBS_COL].map(normalize_day)
day_counts = day_labels.value_counts(dropna=False).sort_index()
display(day_counts.rename("cells").to_frame())

missing_days = {"D0", "D1", "D2", "D3"} - set(day_labels.dropna())
if missing_days:
    raise RuntimeError(
        "Missing expected day labels: " + ", ".join(sorted(missing_days))
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
Cells: 121,616
Features: 17,538

```

### Text output (text/plain)

```text
     cells
day       
D0   38609
D1   16925
D2   21980
D3   44102
```

## Markdown cell 4


## 2. Resolve a raw-count matrix and feature metadata

The notebook will not calculate TPM from log-normalized expression. It prefers an explicit raw-count layer, then `adata.raw`, and uses `adata.X` only if sampled values are count-like non-negative integers.


## Code cell 5

```python

def first_existing_column(frame, candidates):
    lowered = {str(c).lower(): c for c in frame.columns}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def nonzero_values(matrix):
    if sparse.issparse(matrix):
        return np.asarray(matrix.data)
    array = np.asarray(matrix)
    return array[np.isfinite(array) & (array != 0)]


def matrix_is_count_like(accessor, n_obs):
    stop = min(n_obs, 1000)
    sample = accessor(0, stop)
    values = nonzero_values(sample)
    if values.size == 0:
        return True
    values = values[: min(values.size, 100000)]
    return (
        np.nanmin(values) >= 0
        and np.mean(np.isclose(values, np.round(values), atol=1e-6)) >= 0.995
    )


count_source = None
matrix_accessor = None
feature_frame = None
feature_names = None

for layer in cfg.COUNT_LAYER_CANDIDATES:
    if layer in adata.layers:
        count_source = f"adata.layers[{layer!r}]"
        matrix_accessor = lambda start, stop, key=layer: adata.layers[key][start:stop, :]
        feature_frame = adata.var.copy()
        feature_names = pd.Index(adata.var_names.astype(str))
        break

if matrix_accessor is None and adata.raw is not None:
    count_source = "adata.raw.X"
    matrix_accessor = lambda start, stop: adata.raw.X[start:stop, :]
    feature_frame = adata.raw.var.copy()
    feature_names = pd.Index(adata.raw.var_names.astype(str))

if matrix_accessor is None:
    candidate_accessor = lambda start, stop: adata.X[start:stop, :]
    if matrix_is_count_like(candidate_accessor, adata.n_obs):
        count_source = "adata.X (verified count-like)"
        matrix_accessor = candidate_accessor
        feature_frame = adata.var.copy()
        feature_names = pd.Index(adata.var_names.astype(str))
    else:
        raise RuntimeError(
            "No trusted raw-count source was found. TPM cannot be calculated "
            "from log-normalized adata.X. Add a raw-count layer or use an h5ad "
            "that preserves counts."
        )

print("Count source:", count_source)
print("Count matrix shape:", (adata.n_obs, len(feature_names)))
display(feature_frame.head())

```

### Text output (stream)

```text
Count source: adata.layers['raw_counts']
Count matrix shape: (121616, 17538)

```

### Text output (text/plain)

```text
                  ensembl_id gene_symbol
TSPAN6    ENSG00000000003.15      TSPAN6
DPM1      ENSG00000000419.14        DPM1
SCYL3     ENSG00000000457.14       SCYL3
C1orf112  ENSG00000000460.17    C1orf112
FUCA2     ENSG00000001036.14       FUCA2
```

## Markdown cell 6


## 3. Resolve gene symbols, Ensembl IDs, and gene lengths

When no suitable length column exists in `adata.var`/`adata.raw.var`, set `cfg.GENCODE_V43_GTF` to the absolute path of `gencode.v43.annotation.gtf` or `.gtf.gz`. The parser calculates union-exon length per gene.


## Code cell 7

```python

symbol_col = first_existing_column(
    feature_frame,
    cfg.GENE_SYMBOL_COL_CANDIDATES,
)
gene_id_col = first_existing_column(
    feature_frame,
    cfg.GENE_ID_COL_CANDIDATES,
)
length_col = first_existing_column(
    feature_frame,
    cfg.GENE_LENGTH_COL_CANDIDATES,
)

gene_symbols = (
    feature_frame[symbol_col].astype(str).to_numpy()
    if symbol_col is not None
    else feature_names.astype(str).to_numpy()
)
gene_ids = (
    feature_frame[gene_id_col].astype(str).str.replace(r"\.\d+$", "", regex=True).to_numpy()
    if gene_id_col is not None
    else feature_names.astype(str).str.replace(r"\.\d+$", "", regex=True).to_numpy()
)

print("Gene-symbol source:", symbol_col or "feature names")
print("Gene-ID source:", gene_id_col or "feature names")
print("Gene-length column:", length_col)


def parse_gtf_attributes(text):
    result = {}
    for field in text.strip().strip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if " " not in field:
            continue
        key, value = field.split(" ", 1)
        result[key] = value.strip().strip('"')
    return result


def merge_interval_length(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    current_start, current_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start + 1
            current_start, current_end = start, end
    total += current_end - current_start + 1
    return total


def gene_lengths_from_gtf(gtf_path, wanted_ids, wanted_symbols):
    path = Path(gtf_path)
    if not path.exists():
        raise FileNotFoundError(path)

    wanted_ids = set(str(x).split(".")[0] for x in wanted_ids)
    wanted_symbols = set(map(str, wanted_symbols))
    exon_intervals = defaultdict(list)
    gene_alias = {}

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "").split(".")[0]
            gene_name = attrs.get("gene_name", "")
            if gene_id not in wanted_ids and gene_name not in wanted_symbols:
                continue
            key = gene_id if gene_id in wanted_ids else gene_name
            exon_intervals[key].append((int(fields[3]), int(fields[4])))
            gene_alias[key] = (gene_id, gene_name)

    by_id = {}
    by_symbol = {}
    for key, intervals in exon_intervals.items():
        length = merge_interval_length(intervals)
        gene_id, gene_name = gene_alias[key]
        if gene_id:
            by_id[gene_id] = length
        if gene_name:
            by_symbol[gene_name] = length
    return by_id, by_symbol


if length_col is not None:
    gene_lengths_bp = pd.to_numeric(
        feature_frame[length_col],
        errors="coerce",
    ).to_numpy()
    length_source = f"feature metadata column {length_col!r}"
else:
    if cfg.GENCODE_V43_GTF is None:
        raise RuntimeError(
            "No gene-length column was found. Set cfg.GENCODE_V43_GTF to "
            "the absolute Sherlock path of gencode.v43.annotation.gtf or .gtf.gz. "
            "The notebook will not label CPM as TPM."
        )

    lengths_by_id, lengths_by_symbol = gene_lengths_from_gtf(
        cfg.GENCODE_V43_GTF,
        gene_ids,
        gene_symbols,
    )
    gene_lengths_bp = np.array(
        [
            lengths_by_id.get(str(gid).split(".")[0],
                              lengths_by_symbol.get(str(symbol), np.nan))
            for gid, symbol in zip(gene_ids, gene_symbols)
        ],
        dtype=float,
    )
    length_source = f"union-exon lengths from {cfg.GENCODE_V43_GTF}"

print("Gene-length source:", length_source)
print(
    f"Features with valid positive lengths: "
    f"{np.sum(np.isfinite(gene_lengths_bp) & (gene_lengths_bp > 0)):,} / "
    f"{len(gene_lengths_bp):,}"
)

if np.sum(np.isfinite(gene_lengths_bp) & (gene_lengths_bp > 0)) < 1000:
    raise RuntimeError(
        "Too few valid gene lengths were resolved for a reliable TPM denominator."
    )

```

### Text output (stream)

```text
Gene-symbol source: gene_symbol
Gene-ID source: ensembl_id
Gene-length column: None
Gene-length source: union-exon lengths from /oak/stanford/groups/engreitz/Users/tri/tmp/gencode.v43.annotation.gtf.gz
Features with valid positive lengths: 17,532 / 17,538

```

## Markdown cell 8


## 4. Sum raw UMI counts by day in chunks


## Code cell 9

```python

days = ["D0", "D1", "D2", "D3"]
day_to_index = {day: i for i, day in enumerate(days)}
pseudobulk_counts = np.zeros(
    (len(days), len(feature_names)),
    dtype=np.float64,
)

day_array = day_labels.to_numpy()

for start in range(0, adata.n_obs, cfg.CELL_CHUNK_SIZE):
    stop = min(start + cfg.CELL_CHUNK_SIZE, adata.n_obs)
    matrix = matrix_accessor(start, stop)
    chunk_days = day_array[start:stop]

    for day in days:
        local_rows = np.where(chunk_days == day)[0]
        if local_rows.size == 0:
            continue
        block = matrix[local_rows, :]
        summed = np.asarray(block.sum(axis=0)).ravel()
        pseudobulk_counts[day_to_index[day], :] += summed

    if start % (cfg.CELL_CHUNK_SIZE * 10) == 0:
        print(f"Processed {stop:,} / {adata.n_obs:,} cells")

print("Completed pseudobulk aggregation.")

```

### Text output (stream)

```text
Processed 2,000 / 121,616 cells
Processed 22,000 / 121,616 cells
Processed 42,000 / 121,616 cells
Processed 62,000 / 121,616 cells
Processed 82,000 / 121,616 cells
Processed 102,000 / 121,616 cells
Processed 121,616 / 121,616 cells
Completed pseudobulk aggregation.

```

## Markdown cell 10


## 5. Calculate pseudobulk CPM and TPM

For day \(d\) and gene \(g\):

\[
\mathrm{RPK}_{dg} = 
rac{\mathrm{UMI}_{dg}}{\mathrm{gene\ length}_{g}/1000},
\]

\[
\mathrm{TPM}_{dg} =
10^6 	imes

rac{\mathrm{RPK}_{dg}}
{\sum_h \mathrm{RPK}_{dh}}.
\]


## Code cell 11

```python

library_totals = pseudobulk_counts.sum(axis=1)
cpm = pseudobulk_counts / library_totals[:, None] * 1e6

valid_length = np.isfinite(gene_lengths_bp) & (gene_lengths_bp > 0)
rpk = np.zeros_like(pseudobulk_counts, dtype=float)
rpk[:, valid_length] = (
    pseudobulk_counts[:, valid_length]
    / (gene_lengths_bp[valid_length][None, :] / 1000.0)
)

rpk_totals = rpk.sum(axis=1)
if np.any(rpk_totals <= 0):
    raise RuntimeError("At least one day has a non-positive RPK total.")

tpm = rpk / rpk_totals[:, None] * 1e6

summary_rows = []
for family, genes in cfg.PKC_GROUPS.items():
    for gene in genes:
        feature_indices = np.where(gene_symbols == gene)[0]
        if feature_indices.size == 0:
            summary_rows.append(
                {
                    "PKC_family": family,
                    "gene": gene,
                    "status": "not found",
                }
            )
            continue

        if feature_indices.size > 1:
            print(
                f"WARNING: {gene} maps to {feature_indices.size} features; "
                "values are summed across matching features."
            )

        for day_index, day in enumerate(days):
            summary_rows.append(
                {
                    "PKC_family": family,
                    "gene": gene,
                    "day": day,
                    "status": "found",
                    "raw_pseudobulk_UMIs": float(
                        pseudobulk_counts[day_index, feature_indices].sum()
                    ),
                    "pseudobulk_CPM": float(
                        cpm[day_index, feature_indices].sum()
                    ),
                    "pseudobulk_TPM": float(
                        tpm[day_index, feature_indices].sum()
                    ),
                    "gene_length_bp": float(
                        np.nanmax(gene_lengths_bp[feature_indices])
                    ),
                    "n_matching_features": int(feature_indices.size),
                    "count_source": count_source,
                    "length_source": length_source,
                    "cells_in_day": int((day_array == day).sum()),
                }
            )

expression_long = pd.DataFrame(summary_rows)
display(expression_long)

expression_long.to_csv(
    cfg.OUTPUT_DIR / "PKC_isoform_expression_D0_D3_long.csv",
    index=False,
)

expression_wide = (
    expression_long.loc[expression_long["status"].eq("found")]
    .pivot(index=["PKC_family", "gene"], columns="day", values="pseudobulk_TPM")
    .reindex(columns=days)
    .reset_index()
)
expression_wide.to_csv(
    cfg.OUTPUT_DIR / "PKC_isoform_TPM_D0_D3_wide.csv",
    index=False,
)
display(expression_wide)

```

### Text output (text/plain)

```text
          PKC_family   gene day status  raw_pseudobulk_UMIs  pseudobulk_CPM  \
0   Conventional PKC  PRKCA  D0  found               5785.0       47.937226   
1   Conventional PKC  PRKCA  D1  found                890.0       17.563262   
2   Conventional PKC  PRKCA  D2  found               1786.0       43.812842   
3   Conventional PKC  PRKCA  D3  found               2717.0       33.590136   
4   Conventional PKC  PRKCB  D0  found               3466.0       28.720903   
5   Conventional PKC  PRKCB  D1  found                  0.0        0.000000   
6   Conventional PKC  PRKCB  D2  found                200.0        4.906253   
7   Conventional PKC  PRKCB  D3  found                271.0        3.350360   
8   Conventional PKC  PRKCG  D0  found               1350.0       11.186734   
9   Conventional PKC  PRKCG  D1  found                246.0        4.854565   
10  Conventional PKC  PRKCG  D2  found                  0.0        0.000000   
11  Conventional PKC  PRKCG  D3  found                257.0        3.177278   
12         Novel PKC  PRKCD  D0  found                972.0        8.054448   
13         Novel PKC  PRKCD  D1  found                601.0       11.860135   
14         Novel PKC  PRKCD  D2  found                369.0        9.052037   
15         Novel PKC  PRKCD  D3  found                606.0        7.491948   
16         Novel PKC  PRKCE  D0  found               2112.0       17.501023   
17         Novel PKC  PRKCE  D1  found                445.0        8.781631   
18         Novel PKC  PRKCE  D2  found                763.0       18.717356   
19         Novel PKC  PRKCE  D3  found              20341.0      251.474770   
20         Novel PKC  PRKCH  D0  found               2482.0       20.567017   
21         Novel PKC  PRKCH  D1  found                615.0       12.136411   
22         Novel PKC  PRKCH  D2  found                410.0       10.057819   
23         Novel PKC  PRKCH  D3  found               6056.0       74.870026   
24         Novel PKC  PRKCQ  D0  found               2142.0       17.749617   
25         Novel PKC  PRKCQ  D1  found               2983.0       58.866529   
26         Novel PKC  PRKCQ  D2  found                599.0       14.694229   
27         Novel PKC  PRKCQ  D3  found                915.0       11.312099   
28      Atypical PKC  PRKCI  D0  found              11071.0       91.739503   
29      Atypical PKC  PRKCI  D1  found               5085.0      100.347401   
30      Atypical PKC  PRKCI  D2  found               2459.0       60.322384   
31      Atypical PKC  PRKCI  D3  found               3557.0       43.975014   
32      Atypical PKC  PRKCZ  D0  found               3729.0       30.900244   
33      Atypical PKC  PRKCZ  D1  found                975.0       19.240652   
34      Atypical PKC  PRKCZ  D2  found                922.0       22.617828   
35      Atypical PKC  PRKCZ  D3  found                598.0        7.393044   

    pseudobulk_TPM  gene_length_bp  n_matching_features  \
0        12.361978         11175.0                    1   
1         4.708046         11175.0                    1   
2        12.103424         11175.0                    1   
3         9.378055         11175.0                    1   
4         7.490964         11049.0                    1   
5         0.000000         11049.0                    1   
6         1.370823         11049.0                    1   
7         0.946056         11049.0                    1   
8         8.132653          3964.0                    1   
9         3.668594          3964.0                    1   
10        0.000000          3964.0                    1   
11        2.500749          3964.0                    1   
12        4.728304          4909.0                    1   
13        7.237351          4909.0                    1   
14        5.692560          4909.0                    1   
15        4.761571          4909.0                    1   
16        5.054550          9978.0                    1   
17        2.636421          9978.0                    1   
18        5.791024          9978.0                    1   
19       78.632024          9978.0                    1   
20        6.624551          8947.0                    1   
21        4.063459          8947.0                    1   
22        3.470410          8947.0                    1   
23       26.108330          8947.0                    1   
24       15.570990          3285.0                    1   
25       53.680442          3285.0                    1   
26       13.809112          3285.0                    1   
27       10.743762          3285.0                    1   
28       45.755303          5778.0                    1   
29       52.024944          5778.0                    1   
30       32.229629          5778.0                    1   
31       23.745263          5778.0                    1   
32       11.572199          7695.0                    1   
33        7.490213          7695.0                    1   
34        9.073955          7695.0                    1   
35        2.997528          7695.0                    1   

                  count_source  \
0   adata.layers['raw_counts']   
1   adata.layers['raw_counts']   
2   adata.layers['raw_counts']   
3   adata.layers['raw_counts']   
4   adata.layers['raw_counts']   
5   adata.layers['raw_counts']   
6   adata.layers['raw_counts']   
7   adata.layers['raw_counts']   
8   adata.layers['raw_counts']   
9   adata.layers['raw_counts']   
10  adata.layers['raw_counts']   
11  adata.layers['raw_counts']   
12  adata.layers['raw_counts']   
13  adata.layers['raw_counts']   
14  adata.layers['raw_counts']   
15  adata.layers['raw_counts']   
16  adata.layers['raw_counts']   
17  adata.layers['raw_counts']   
18  adata.layers['raw_counts']   
19  adata.layers['raw_counts']   
20  adata.layers['raw_counts']   
21  adata.layers['raw_counts']   
22  adata.layers['raw_counts']   
23  adata.layers['raw_counts']   
24  adata.layers['raw_counts']   
25  adata.layers['raw_counts']   
26  adata.layers['raw_counts']   
27  adata.layers['raw_counts']   
28  adata.layers['raw_counts']   
29  adata.layers['raw_counts']   
30  adata.layers['raw_counts']   
31  adata.layers['raw_counts']   
32  adata.layers['raw_counts']   
33  adata.layers['raw_counts']   
34  adata.layers['raw_counts']   
35  adata.layers['raw_counts']   

                                        length_source  cells_in_day  
0   union-exon lengths from /oak/stanford/groups/e...         38609  
1   union-exon lengths from /oak/stanford/groups/e...         16925  
2   union-exon lengths from /oak/stanford/groups/e...         21980  
3   union-exon lengths from /oak/stanford/groups/e...         44102  
4   union-exon lengths from /oak/stanford/groups/e...         38609  
5   union-exon lengths from /oak/stanford/groups/e...         16925  
6   union-exon lengths from /oak/stanford/groups/e...         21980  
7   union-exon lengths from /oak/stanford/groups/e...         44102  
8   union-exon lengths from /oak/stanford/groups/e...         38609  
9   union-exon lengths from /oak/stanford/groups/e...         16925  
10  union-exon lengths from /oak/stanford/groups/e...         21980  
11  union-exon lengths from /oak/stanford/groups/e...         44102  
12  union-exon lengths from /oak/stanford/groups/e...         38609  
13  union-exon lengths from /oak/stanford/groups/e...         16925  
14  union-exon lengths from /oak/stanford/groups/e...         21980  
15  union-exon lengths from /oak/stanford/groups/e...         44102  
16  union-exon lengths from /oak/stanford/groups/e...         38609  
17  union-exon lengths from /oak/stanford/groups/e...         16925  
18  union-exon lengths from /oak/stanford/groups/e...         21980  
19  union-exon lengths from /oak/stanford/groups/e...         44102  
20  union-exon lengths from /oak/stanford/groups/e...         38609  
21  union-exon lengths from /oak/stanford/groups/e...         16925  
22  union-exon lengths from /oak/stanford/groups/e...         21980  
23  union-exon lengths from /oak/stanford/groups/e...         44102  
24  union-exon lengths from /oak/stanford/groups/e...         38609  
25  union-exon lengths from /oak/stanford/groups/e...         16925  
26  union-exon lengths from /oak/stanford/groups/e...         21980  
27  union-exon lengths from /oak/stanford/groups/e...         44102  
28  union-exon lengths from /oak/stanford/groups/e...         38609  
29  union-exon lengths from /oak/stanford/groups/e...         16925  
30  union-exon lengths from /oak/stanford/groups/e...         21980  
31  union-exon lengths from /oak/stanford/groups/e...         44102  
32  union-exon lengths from /oak/stanford/groups/e...         38609  
33  union-exon lengths from /oak/stanford/groups/e...         16925  
34  union-exon lengths from /oak/stanford/groups/e...         21980  
35  union-exon lengths from /oak/stanford/groups/e...         44102  
```

### Text output (text/plain)

```text
day        PKC_family   gene         D0         D1         D2         D3
0        Atypical PKC  PRKCI  45.755303  52.024944  32.229629  23.745263
1        Atypical PKC  PRKCZ  11.572199   7.490213   9.073955   2.997528
2    Conventional PKC  PRKCA  12.361978   4.708046  12.103424   9.378055
3    Conventional PKC  PRKCB   7.490964   0.000000   1.370823   0.946056
4    Conventional PKC  PRKCG   8.132653   3.668594   0.000000   2.500749
5           Novel PKC  PRKCD   4.728304   7.237351   5.692560   4.761571
6           Novel PKC  PRKCE   5.054550   2.636421   5.791024  78.632024
7           Novel PKC  PRKCH   6.624551   4.063459   3.470410  26.108330
8           Novel PKC  PRKCQ  15.570990  53.680442  13.809112  10.743762
```

## Markdown cell 12


## 6. Plot PKC isoform TPM across D0-D3

Both a line plot and a heatmap are displayed and exported as PNG, PDF and Illustrator-editable SVG.


## Markdown cell 13

### Illustrator editability

SVG files are exported with `svg.fonttype = "none"` and no LaTeX/mathtext labels, so all plot labels remain native, editable text objects in Adobe Illustrator. PDF files use TrueType text (`pdf.fonttype = 42`).

## Code cell 14

```python
plot_long = expression_long.loc[
    expression_long["status"].eq("found")
].copy()

plot_long["day"] = pd.Categorical(
    plot_long["day"],
    categories=days,
    ordered=True,
)
plot_long["gene"] = pd.Categorical(
    plot_long["gene"],
    categories=PKC_GENES,
    ordered=True,
)
family_order = list(cfg.PKC_GROUPS.keys())
plot_long["PKC_family"] = pd.Categorical(
    plot_long["PKC_family"],
    categories=family_order,
    ordered=True,
)

line_plot = (
    ggplot(
        plot_long,
        aes(
            x="day",
            y="pseudobulk_TPM",
            group="gene",
            color="gene",
        ),
    )
    + geom_line(size=0.9, alpha=0.9)
    + geom_point(size=2.2, alpha=0.95)
    + facet_grid(
        "PKC_family ~ .",
        scales="free_y",
        space="free_y",
    )
    + scale_x_discrete(limits=days)
    + scale_color_brewer(type="qual", palette="Dark2")
    + labs(
        title="PKC isoform expression during hiPSC-to-EC differentiation",
        subtitle="Day-level pseudobulk expression from raw UMI counts",
        x="Differentiation day",
        y="Pseudobulk TPM",
        color=None,
    )
    + theme_classic(base_size=10.5, base_family="Liberation Sans")
    + theme(
        figure_size=(7.0, 6.8),
        text=element_text(family="Liberation Sans"),
        plot_title=element_text(size=14, weight="bold", ha="left"),
        plot_subtitle=element_text(size=9.5, color="#555555", ha="left"),
        axis_title_x=element_text(size=10.5),
        axis_title_y=element_text(size=10.5),
        axis_text_x=element_text(size=9.5),
        axis_text_y=element_text(size=9),
        strip_background=element_rect(
            fill="#F2F2F2",
            color="#D9D9D9",
            size=0.6,
        ),
        strip_text_y=element_text(
            size=9.5,
            weight="bold",
            angle=0,
        ),
        legend_position="right",
        legend_text=element_text(size=8.3),
        panel_spacing_y=0.10,
    )
)

display(line_plot)
for extension in ["png", "pdf", "svg"]:
    line_plot.save(
        cfg.OUTPUT_DIR / f"PKC_isoform_TPM_D0_D3_plotnine_lineplot_editable_text.{extension}",
        format=extension,
        dpi=350 if extension == "png" else 300,
        verbose=False,
    )

heatmap_long = plot_long[
    ["PKC_family", "gene", "day", "pseudobulk_TPM"]
].copy()
heatmap_long["label"] = heatmap_long["pseudobulk_TPM"].map(
    lambda value: f"{value:.1f}" if pd.notna(value) else ""
)

heatmap_plot = (
    ggplot(
        heatmap_long,
        aes(
            x="day",
            y="gene",
            fill="pseudobulk_TPM",
        ),
    )
    + geom_tile(color="white", size=0.5)
    + geom_text(
        aes(label="label"),
        size=8.0,
        color="#202020",
    )
    + facet_grid(
        "PKC_family ~ .",
        scales="free_y",
        space="free_y",
    )
    + scale_x_discrete(limits=days)
    + scale_fill_cmap(name="TPM", cmap_name="YlGnBu")
    + labs(
        title="PKC isoform expression across D0-D3",
        subtitle="Day-level pseudobulk TPM",
        x="Differentiation day",
        y=None,
    )
    + theme_minimal(base_size=10.5, base_family="Liberation Sans")
    + theme(
        figure_size=(5.4, 6.2),
        text=element_text(family="Liberation Sans"),
        plot_title=element_text(size=14, weight="bold", ha="left"),
        plot_subtitle=element_text(size=9.5, color="#555555", ha="left"),
        axis_title_x=element_text(size=10.5),
        axis_text_x=element_text(size=9.5),
        axis_text_y=element_text(size=10),
        panel_grid=element_blank(),
        strip_background=element_rect(
            fill="#F2F2F2",
            color="#D9D9D9",
            size=0.6,
        ),
        strip_text_y=element_text(
            size=9.5,
            weight="bold",
            angle=0,
        ),
        legend_position="right",
        panel_spacing_y=0.10,
    )
)

display(heatmap_plot)
for extension in ["png", "pdf", "svg"]:
    heatmap_plot.save(
        cfg.OUTPUT_DIR / f"PKC_isoform_TPM_D0_D3_plotnine_heatmap_editable_text.{extension}",
        format=extension,
        dpi=350 if extension == "png" else 300,
        verbose=False,
    )


# Confirm that Illustrator-editable SVG text objects were preserved.
for svg_name in [
    "PKC_isoform_TPM_D0_D3_plotnine_lineplot_editable_text.svg",
    "PKC_isoform_TPM_D0_D3_plotnine_heatmap_editable_text.svg",
]:
    svg_path = cfg.OUTPUT_DIR / svg_name
    svg_source = svg_path.read_text(encoding="utf-8")
    if "<text" not in svg_source:
        raise RuntimeError(f"No native text objects found in {svg_path}")
    print(f"{svg_name}: {svg_source.count('<text')} editable SVG text objects")

```

### Text output (stream)

```text
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7fdcf1c7b650>
```

### Text output (stream)

```text
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7fdcf290b750>
```

### Text output (stream)

```text
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.
findfont: Font family 'Liberation Sans' not found.

```

### Text output (stream)

```text
PKC_isoform_TPM_D0_D3_plotnine_lineplot_editable_text.svg: 35 editable SVG text objects
PKC_isoform_TPM_D0_D3_plotnine_heatmap_editable_text.svg: 61 editable SVG text objects

```

## Markdown cell 15


## 7. Close the backed AnnData handle


## Code cell 16

```python

try:
    adata.file.close()
except Exception:
    pass

print("Outputs:", cfg.OUTPUT_DIR.resolve())

```

### Text output (stream)

```text
Outputs: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/PKC_isoform_expression_D0_D3_outputs

```

