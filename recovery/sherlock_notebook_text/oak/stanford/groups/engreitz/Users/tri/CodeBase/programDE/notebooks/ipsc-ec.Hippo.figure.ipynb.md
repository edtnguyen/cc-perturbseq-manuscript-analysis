# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/ipsc-ec.Hippo.figure.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/ipsc-ec.Hippo.figure.ipynb`
- SHA-256: `f623ec899211436fadd5764e48ecf2255fda5f5ba352b781780549dbb3849774`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Code cell 0

```python
%reset -f
```

## Code cell 1

```python
import sys
print(sys.executable)
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python

```

## Code cell 2

```python
import scanpy as sc
import muon as mu
import pandas as pd
import numpy as np
import os

```

## Code cell 3

```python
adata_int = sc.read_h5ad('/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad')
```

## Code cell 4

```python
adata_int
```

### Text output (text/plain)

```text
AnnData object with n_obs × n_vars = 121616 × 17538
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'diff_day', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2'
    var: 'ensembl_id', 'gene_symbol'
    uns: 'guide_names', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2', 'log1p', 'neighbors'
    obsm: 'X_pca', 'X_umap', 'guide_assignment'
    layers: 'raw_counts'
    obsp: 'connectivities', 'distances'
```

## Code cell 5

```python
adata_int.shape
```

### Text output (text/plain)

```text
(121616, 17538)
```

## Code cell 6

```python
#import scanpy as sc

# 1. Normalize and Log-transform (if starting from raw counts)
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

# 2. Identify Highly Variable Genes (essential for clean clusters)
# sc.pp.highly_variable_genes(adata_int, min_mean=0.0125, max_mean=3, min_disp=0.5)

# 3. Filter to HVGs and Scale
# We store the full data in .raw before filtering
# adata_int.raw = adata_int
# adata_int = adata_int[:, adata_int.var.highly_variable].copy()
# sc.pp.scale(adata_int, max_value=10)
```

## Code cell 7

```python
# sc.tl.pca(adata_int, svd_solver='arpack')

# Optional: Check the elbow plot to see if 50 PCs is enough
# sc.pl.pca_variance_ratio(adata, log=True)
# 1. Compute Neighbors
# Increase n_neighbors (e.g., 20 or 30) for more "global" structure
# sc.pp.neighbors(adata_int, n_neighbors=15, n_pcs=40)

# 2. Compute UMAP
# To spread clusters out more, increase 'min_dist' (default is 0.5)
# sc.tl.umap(adata_int, min_dist=0.5, spread=1.0)
```

## Code cell 8

```python
# Re-run Leiden at your preferred resolution
# sc.tl.leiden(adata_int, resolution=0.5, key_added='leiden_new')
```

## Code cell 9

```python
sc.pl.umap(adata_int, 
           color='leiden_res_0.5', 
           legend_loc='on data', 
           legend_fontsize=9, 
           #add_outline=True, 
           frameon=False, 
           size=2)

cluster_mapping = {
    '0': 'artEC',   
    '1': 'iPSC',
    '2': 'Neuro-iPSC',
    '3': 'PrimStreak', 
    '4': 'Lat-Meso',
    '5': 'Fibro-Mesen',
    '6': 'Pre-EC',
    '7': 'Neu-Mesen',
    '8': 'Others',
}

# 2. Map the labels to a new column in adata.obs
# Replace 'leiden_res_0.8' with whatever your current clustering column is named
adata_int.obs['cell_type'] = adata_int.obs['leiden_res_0.5'].astype(str).map(cluster_mapping)

# Convert the new column to a categorical type (important for Scanpy plotting colors)
adata_int.obs['cell_type'] = adata_int.obs['cell_type'].astype('category')
```

### Text output (text/plain)

```text
<Figure size 640x480 with 1 Axes>
```

## Code cell 10

```python
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import matplotlib.pyplot as plt

CLUSTER_COL = "leiden_res_0.5"
CLUSTER_A = "4"
CLUSTER_B = "6"

OUTDIR = "DE_leiden_res_0.5_cluster4_vs_6"
os.makedirs(OUTDIR, exist_ok=True)

# Use adata_int.X by default.
# This should be log-normalized expression, not raw counts.
# If you have a log-normalized layer, set DE_LAYER to that layer name.
DE_LAYER = None
# Example:
# DE_LAYER = "log1p_norm"

# Make sure cluster labels are strings.
adata_int.obs[CLUSTER_COL] = adata_int.obs[CLUSTER_COL].astype(str)

mask = adata_int.obs[CLUSTER_COL].isin([CLUSTER_A, CLUSTER_B])
adata_de = adata_int[mask].copy()

adata_de.obs[CLUSTER_COL] = adata_de.obs[CLUSTER_COL].astype(str)
adata_de.obs[CLUSTER_COL] = pd.Categorical(
    adata_de.obs[CLUSTER_COL],
    categories=[CLUSTER_A, CLUSTER_B],
    ordered=True
)

print(adata_de)
print(adata_de.obs[CLUSTER_COL].value_counts())
```

### Text output (stream)

```text
AnnData object with n_obs × n_vars = 21983 × 17538
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'diff_day', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2', 'cell_type'
    var: 'ensembl_id', 'gene_symbol'
    uns: 'guide_names', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2', 'log1p', 'neighbors', 'leiden_res_0.5_colors'
    obsm: 'X_pca', 'X_umap', 'guide_assignment'
    layers: 'raw_counts'
    obsp: 'connectivities', 'distances'
leiden_res_0.5
4    16168
6     5815
Name: count, dtype: int64

```

## Markdown cell 11

##### End of Finding DE genes

## Code cell 12

```python
import matplotlib.pyplot as plt

days = sorted(adata_int.obs['day'].unique())

# Adjusted figsize to be more square for a 2x2 layout
fig, axs = plt.subplots(2, 2, figsize=(10, 10)) 

# Flatten the 2x2 matrix into a 1D list of 4 axes
axs = axs.flatten()

for i, day in enumerate(days):
    sc.pl.umap(adata_int, 
               color='day', 
               groups=[day], 
               legend_loc=None,       
               frameon=False, 
               size=2,
               title=f'{day}',
               ax=axs[i],             # Now this works properly again!
               show=False)            

plt.tight_layout()
plt.show()
```

### Text output (text/plain)

```text
<Figure size 1000x1000 with 4 Axes>
```

## Code cell 13

```python
adata_int
```

### Text output (text/plain)

```text
AnnData object with n_obs × n_vars = 121616 × 17538
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'diff_day', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2', 'cell_type'
    var: 'ensembl_id', 'gene_symbol'
    uns: 'guide_names', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2', 'log1p', 'neighbors', 'leiden_res_0.5_colors', 'day_colors'
    obsm: 'X_pca', 'X_umap', 'guide_assignment'
    layers: 'raw_counts'
    obsp: 'connectivities', 'distances'
```

## Code cell 14

```python
adata = adata_int.copy()
```

## Markdown cell 15

#### Import cNMF

## Code cell 16

```python
#1. Define the file path
cnmf_file = "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50.usages.k_50.dt_0_2.consensus.txt"

# 2. Read the dataframe
usage_df = pd.read_csv(cnmf_file, sep='\t', index_col=0)

# 3. Align with adata
# This ensures every cell in adata gets its corresponding row from the file.
# Cells in adata but missing from the file will be filled with 0 
usage_aligned = usage_df.reindex(adata_int.obs_names).fillna(0)

# 4. Save to obsm
adata_int.obsm["cnmf_usage"] = usage_aligned

# 5. Validation
print(f"Loaded usage matrix shape: {adata_int.obsm['cnmf_usage'].shape}")
print(f"Overlap check: {usage_df.index.intersection(adata_int.obs_names).shape[0]} cells matched.")
```

### Text output (stream)

```text
Loaded usage matrix shape: (121616, 50)
Overlap check: 121616 cells matched.

```

## Code cell 17

```python
# Rename columns to program_0, program_1, ..., program_49
adata_int.obsm["cnmf_usage"].columns = [f"Program {i}" for i in range(adata_int.obsm["cnmf_usage"].shape[1])]

# Verify the changes
print(adata_int.obsm["cnmf_usage"].head())
```

### Text output (stream)

```text
                            Program 0  Program 1  Program 2  Program 3  \
barcode                                                                  
CCTAATCCAACAACCATCTATTAC_2   0.062538   0.027402   0.012862   0.000000   
GGTGCGAACCTCCTGAATCCGCGA_2   0.000000   0.019443   0.003818   0.012402   
AACTCACCATCATTCCCTGTCCCG_1   0.047637   0.002082   0.000000   0.000000   
TTCACGCACACTTCGAATTCATGG_2   0.076358   0.091694   0.019964   0.000000   
GTCTGTCAGAGCTGAACCGTTCTA_3   0.021434   0.028805   0.000000   0.022261   

                            Program 4  Program 5  Program 6  Program 7  \
barcode                                                                  
CCTAATCCAACAACCATCTATTAC_2   0.050164   0.127009   0.126219   0.086548   
GGTGCGAACCTCCTGAATCCGCGA_2   0.030527   0.123209   0.065145   0.046808   
AACTCACCATCATTCCCTGTCCCG_1   0.098233   0.006417   0.069543   0.088016   
TTCACGCACACTTCGAATTCATGG_2   0.000000   0.000000   0.550251   0.000000   
GTCTGTCAGAGCTGAACCGTTCTA_3   0.090041   0.038934   0.155332   0.067275   

                            Program 8  Program 9  ...  Program 40  Program 41  \
barcode                                           ...                           
CCTAATCCAACAACCATCTATTAC_2   0.004646   0.039850  ...    0.000000    0.014634   
GGTGCGAACCTCCTGAATCCGCGA_2   0.028964   0.000000  ...    0.002084    0.000000   
AACTCACCATCATTCCCTGTCCCG_1   0.049213   0.000000  ...    0.000000    0.001022   
TTCACGCACACTTCGAATTCATGG_2   0.000000   0.000000  ...    0.000000    0.000000   
GTCTGTCAGAGCTGAACCGTTCTA_3   0.000000   0.023845  ...    0.003669    0.021365   

                            Program 42  Program 43  Program 44  Program 45  \
barcode                                                                      
CCTAATCCAACAACCATCTATTAC_2    0.000000    0.000452    0.018801    0.004714   
GGTGCGAACCTCCTGAATCCGCGA_2    0.000143    0.005316    0.000609    0.006756   
AACTCACCATCATTCCCTGTCCCG_1    0.000000    0.000000    0.000000    0.000000   
TTCACGCACACTTCGAATTCATGG_2    0.000000    0.000000    0.000000    0.006451   
GTCTGTCAGAGCTGAACCGTTCTA_3    0.000062    0.003521    0.000000    0.000006   

                            Program 46  Program 47  Program 48  Program 49  
barcode                                                                     
CCTAATCCAACAACCATCTATTAC_2         0.0    0.000000         0.0    0.001114  
GGTGCGAACCTCCTGAATCCGCGA_2         0.0    0.000000         0.0    0.000000  
AACTCACCATCATTCCCTGTCCCG_1         0.0    0.009507         0.0    0.000000  
TTCACGCACACTTCGAATTCATGG_2         0.0    0.143893         0.0    0.000000  
GTCTGTCAGAGCTGAACCGTTCTA_3         0.0    0.000513         0.0    0.000000  

[5 rows x 50 columns]

```

## Code cell 18

```python
import scanpy as sc

# 1. Loop through all the 'program_X' columns and add them to .obs
for prog in adata_int.obsm["cnmf_usage"].columns:
    # Using .values ensures we bypass any potential index mismatch issues
    adata_int.obs[prog] = adata_int.obsm["cnmf_usage"][prog].values

# Verify they were added
print(adata_int.obs.columns[-5:]) # Look at the last 5 columns added
```

### Text output (stream)

```text
Index(['Program 45', 'Program 46', 'Program 47', 'Program 48', 'Program 49'], dtype='object')

```

## Code cell 19

```python
sc.pl.umap(adata_int, 
           color=['Program 17', 'Program 14','Program 23', 'Program 15', 'Program 19'], 
           legend_loc='on data', 
           legend_fontsize=9, 
           add_outline=True, 
           frameon=False,
           #scale_factor=5,
           cmap="YlGnBu",
           vmin=0, vmax='p99',
           ncols=2,
           size=10)
```

### Text output (text/plain)

```text
<Figure size 1455.6x1440 with 10 Axes>
```

## Code cell 20

```python
sc.pl.umap(
    adata_int,
    color=[
        #'Program 15', 'Program 17', 'Program 14',

        # P15
        #'KDR', 'TAL1', 'FLI1', 'APLNR', 'NRP2',

        # P17
        #'CXCL12', 'TFAP2A', 'KRT19', 'CDH1', 'SLC8A1',

        # P14
        'TBX20', 'HAND1', 'CDH11', 'ALCAM', 'ROR1'
    ],
    legend_loc='on data',
    legend_fontsize=9,
    add_outline=True,
    frameon=False,
    cmap="YlGnBu",
    vmin=0,
    vmax=0.4,
    ncols=3,
    size=10
)
```

### Text output (text/plain)

```text
<Figure size 2183.4x960 with 10 Axes>
```

## Code cell 21

```python
adata_int
```

### Text output (text/plain)

```text
AnnData object with n_obs × n_vars = 121616 × 17538
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'diff_day', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2', 'cell_type', 'Program 0', 'Program 1', 'Program 2', 'Program 3', 'Program 4', 'Program 5', 'Program 6', 'Program 7', 'Program 8', 'Program 9', 'Program 10', 'Program 11', 'Program 12', 'Program 13', 'Program 14', 'Program 15', 'Program 16', 'Program 17', 'Program 18', 'Program 19', 'Program 20', 'Program 21', 'Program 22', 'Program 23', 'Program 24', 'Program 25', 'Program 26', 'Program 27', 'Program 28', 'Program 29', 'Program 30', 'Program 31', 'Program 32', 'Program 33', 'Program 34', 'Program 35', 'Program 36', 'Program 37', 'Program 38', 'Program 39', 'Program 40', 'Program 41', 'Program 42', 'Program 43', 'Program 44', 'Program 45', 'Program 46', 'Program 47', 'Program 48', 'Program 49'
    var: 'ensembl_id', 'gene_symbol'
    uns: 'guide_names', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2', 'log1p', 'neighbors', 'leiden_res_0.5_colors', 'day_colors'
    obsm: 'X_pca', 'X_umap', 'guide_assignment', 'cnmf_usage'
    layers: 'raw_counts'
    obsp: 'connectivities', 'distances'
```

## Code cell 22

```python
adata = adata_int.copy()
```

## Markdown cell 23

#### UMAP neighborhood panel

## Code cell 24

```python
import os
import re
import zipfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from scipy import sparse
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# ============================================================
# Output
# ============================================================

OUTDIR = Path("carrier_seeded_Hippo_all_days_exports")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# adata fields
# ============================================================

DAY_COL = "day"
CELLTYPE_COL = "cell_type"

GUIDE_NAMES_KEY = "guide_names"
GUIDE_ASSIGNMENT_KEY = "guide_assignment"
CNMF_KEY = "cnmf_usage"
UMAP_KEY = "X_umap"

# Guide metadata.
GUIDE_METADATA_TSV = "300genes_guide_metadata_v43.tsv"
GUIDE_ID_COL = "guide_id"
TARGET_COL = "putative_target_genes"

# ============================================================
# Targets
# ============================================================

HIPPO_TARGETS = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]

CONTROL_TARGET = "non-targeting_control"

CONTROL_PATTERNS = [
    r"non[-_ ]?target",
    r"non[-_ ]?targeting",
    r"\bntc\b",
    r"negative",
    r"neg[-_ ]?ctrl",
    r"scram",
    r"scrambled",
    r"safe",
    r"aavs1",
    r"control",
]

# ============================================================
# cNMF / KNN settings
# ============================================================

# If adata.obsm["cnmf_usage"] corresponds to Program 0..49, use 0.
CNMF_ARRAY_PROGRAM_INDEX_BASE = 0

# Use all cNMF programs for KNN by default.
# Example subset: [13, 15, 17, 23, 27, 30, 34, 35]
PROGRAMS_FOR_KNN = None

# "zscore" or "clr_zscore"
CNMF_TRANSFORM = "zscore"

# Neighborhood size.
K = 300
PSEUDOCOUNT = 0.5

# Important for developmental time series:
# True means D2 seeds only compare to D2 neighbors and D2 NT background.
RESTRICT_NEIGHBORS_TO_SAME_DAY = True
SAME_DAY_BASELINE = True

# Carrier modes:
# "all_carriers"
# "low_guide_load_leq2"
# "single_selected_target"
MODE = "all_carriers"

# Days to include on the shared UMAP.
DAYS_TO_USE = ["d0", "d1", "d2", "d3"]

# Plotting/export helper settings.
LOG2FC_CLIP = 3.0
RANDOM_SEED = 1

# Optional: export seed-neighbor membership table.
# This can get large, but is useful for debugging.
EXPORT_NEIGHBOR_MEMBERSHIP = False
```

## Code cell 25

```python
def looks_control(x):
    s = str(x).lower()
    return any(re.search(p, s) for p in CONTROL_PATTERNS)


def clean_target_string(x):
    """
    Convert putative_target_genes entry into a list of target gene names.
    Handles:
      NF2
      NF2;LATS2
      NF2,LATS2
      ['NF2']
    """
    if pd.isna(x):
        return []

    s = str(x).strip()

    if s.lower() in ["", "nan", "none", "na"]:
        return []

    s = (
        s.replace("[", "")
         .replace("]", "")
         .replace("'", "")
         .replace('"', "")
    )

    parts = re.split(r"[;,|/]+", s)
    parts = [p.strip() for p in parts if p.strip() != ""]

    return parts


def canonical_guide_id(x):
    """
    Conservative guide ID normalization for fallback matching.
    """
    return str(x).strip().replace(" ", "")


def get_guide_assignment_matrix(adata):
    ga = adata.obsm[GUIDE_ASSIGNMENT_KEY]

    if sparse.issparse(ga):
        return ga.tocsr().astype(bool)

    return np.asarray(ga).astype(bool)


def any_guides_per_cell(guide_assignment, guide_indices):
    """
    Boolean mask for cells carrying any guide in guide_indices.
    Works with dense or sparse guide_assignment.
    """
    if len(guide_indices) == 0:
        return np.zeros(guide_assignment.shape[0], dtype=bool)

    if sparse.issparse(guide_assignment):
        return np.asarray(guide_assignment[:, guide_indices].sum(axis=1)).ravel() > 0

    return guide_assignment[:, guide_indices].any(axis=1)


def load_guide_target_from_metadata(adata):
    """
    Load guide_name -> target_name from 300genes_guide_metadata_v43.tsv.
    """

    guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))
    guide_name_to_idx = {g: i for i, g in enumerate(guide_names)}

    meta = pd.read_csv(GUIDE_METADATA_TSV, sep="\t")

    if GUIDE_ID_COL not in meta.columns:
        raise ValueError(
            f"{GUIDE_ID_COL!r} not found in {GUIDE_METADATA_TSV}. "
            f"Available columns: {meta.columns.tolist()}"
        )

    if TARGET_COL not in meta.columns:
        raise ValueError(
            f"{TARGET_COL!r} not found in {GUIDE_METADATA_TSV}. "
            f"Available columns: {meta.columns.tolist()}"
        )

    meta = meta[[GUIDE_ID_COL, TARGET_COL]].copy()
    meta.columns = ["guide_id", "putative_target_genes"]
    meta["guide_id"] = meta["guide_id"].astype(str)

    # Exact guide matching.
    meta["guide_index"] = meta["guide_id"].map(guide_name_to_idx)
    meta["match_method"] = np.where(meta["guide_index"].notna(), "exact_guide_id_match", None)

    # Fallback canonical matching.
    if meta["guide_index"].isna().any():
        canonical_to_guide = {}
        duplicated = set()

        for g in guide_names:
            cg = canonical_guide_id(g)
            if cg in canonical_to_guide:
                duplicated.add(cg)
            else:
                canonical_to_guide[cg] = g

        for cg in duplicated:
            canonical_to_guide.pop(cg, None)

        for idx in meta.index[meta["guide_index"].isna()]:
            cg = canonical_guide_id(meta.at[idx, "guide_id"])
            if cg in canonical_to_guide:
                matched_guide = canonical_to_guide[cg]
                meta.at[idx, "guide_index"] = guide_name_to_idx[matched_guide]
                meta.at[idx, "match_method"] = "canonical_guide_id_match"

    mapped = meta[meta["guide_index"].notna()].copy()
    mapped["guide_index"] = mapped["guide_index"].astype(int)
    mapped["guide_name"] = [guide_names[i] for i in mapped["guide_index"]]

    records = []

    for _, row in mapped.iterrows():
        guide_index = int(row["guide_index"])
        guide_name = str(row["guide_name"])
        raw_target = row["putative_target_genes"]

        targets = clean_target_string(raw_target)

        if looks_control(guide_name) or looks_control(raw_target):
            targets = [CONTROL_TARGET]

        for target in targets:
            records.append({
                "guide_index": guide_index,
                "guide_name": guide_name,
                "target_name": str(target),
                "raw_target_entry": raw_target,
                "match_method": row["match_method"],
            })

    guide_target = pd.DataFrame(records).drop_duplicates()

    print("\nGuide metadata mapping QC")
    print(f"adata guide columns:            {len(guide_names):,}")
    print(f"metadata guide rows:            {len(meta):,}")
    print(f"metadata rows matched to adata: {mapped['guide_id'].nunique():,}")
    print(f"expanded guide-target rows:     {len(guide_target):,}")

    if len(guide_target) > 0:
        print(f"unique mapped targets:          {guide_target['target_name'].nunique():,}")
        print("\nMatch methods:")
        print(guide_target["match_method"].value_counts().to_string())

        print("\nExample guide-target rows:")
        print(guide_target.head(20).to_string(index=False))

    matched_adata_guides = set(guide_target["guide_name"].astype(str)) if len(guide_target) else set()
    unmapped = [g for g in guide_names if g not in matched_adata_guides]

    print(f"\nadata guide columns without metadata target: {len(unmapped):,}")
    if len(unmapped) > 0:
        print("First 30 unmapped adata guide names:")
        print(pd.Series(unmapped[:30]).to_string(index=False))

    guide_target.to_csv(OUTDIR / "guide_to_target_map_from_300genes_metadata.csv", index=False)

    return guide_target


def build_target_assignment(adata, guide_target):
    """
    Build target-level assignment DataFrame from guide assignment.
    Columns are target genes, rows are cells.
    """

    guide_assignment = get_guide_assignment_matrix(adata)

    if guide_assignment.shape[1] != len(adata.uns[GUIDE_NAMES_KEY]):
        raise ValueError(
            f"guide_assignment has {guide_assignment.shape[1]} columns but "
            f"adata.uns['{GUIDE_NAMES_KEY}'] has {len(adata.uns[GUIDE_NAMES_KEY])} names."
        )

    series_by_target = {}
    target_to_guide_indices = {}

    for target, sub in guide_target.groupby("target_name"):
        idx = sorted(set(sub["guide_index"].astype(int).tolist()))
        target_to_guide_indices[target] = idx
        series_by_target[target] = any_guides_per_cell(guide_assignment, idx)

    target_assignment = pd.DataFrame(series_by_target, index=adata.obs_names).astype(bool)

    if CONTROL_TARGET not in target_assignment.columns:
        control_like = [c for c in target_assignment.columns if looks_control(c)]
        if len(control_like) == 0:
            raise ValueError("No non-targeting/control guides found.")
        target_assignment[CONTROL_TARGET] = target_assignment[control_like].any(axis=1)

    print("\nTarget assignment QC")
    print(f"Targets with masks: {target_assignment.shape[1]:,}")
    print(f"Control carrier cells: {target_assignment[CONTROL_TARGET].sum():,}")

    print("\nKey target guide/cell counts:")
    for t in HIPPO_TARGETS + ["KDR", "ETV2", "BMPR1A", "SMAD4"]:
        matches = [x for x in target_assignment.columns if x.upper() == t.upper()]
        if matches:
            tt = matches[0]
            print(
                f"  {t:8s}: guides={len(target_to_guide_indices[tt]):3d}, "
                f"carrier_cells={target_assignment[tt].sum():6d}"
            )
        else:
            print(f"  {t:8s}: NOT FOUND")

    target_counts = pd.DataFrame({
        "target": target_assignment.columns,
        "n_carrier_cells": target_assignment.sum(axis=0).values,
        "n_guides": [len(target_to_guide_indices.get(t, [])) for t in target_assignment.columns],
    })
    target_counts.to_csv(OUTDIR / "target_carrier_counts.csv", index=False)

    return target_assignment, target_to_guide_indices
```

## Code cell 26

```python
def standardize_program_col(c):
    s = str(c)

    if s.startswith("Program "):
        return f"Program {int(s.replace('Program ', ''))}"

    if s.startswith("P"):
        return f"Program {int(s.replace('P', ''))}"

    return f"Program {int(float(s))}"


def get_program_score_df(adata):
    """
    Get cNMF program usage as Program 0, Program 1, ...
    from adata.obsm[CNMF_KEY].
    """
    X = adata.obsm[CNMF_KEY]

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
    cols = [f"Program {i + CNMF_ARRAY_PROGRAM_INDEX_BASE}" for i in range(X.shape[1])]

    return pd.DataFrame(X, index=adata.obs_names, columns=cols).astype(float).copy()


def make_cnmf_knn_matrix(program_df, programs_for_knn=None, transform="zscore"):
    """
    Returns X_knn and the program columns used for KNN.
    """
    if programs_for_knn is None:
        cols = [c for c in program_df.columns if str(c).startswith("Program ")]
    else:
        cols = [f"Program {int(p)}" for p in programs_for_knn]

    missing = [c for c in cols if c not in program_df.columns]
    if missing:
        raise KeyError(f"Missing requested cNMF programs: {missing}")

    X = program_df[cols].copy()
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # Remove zero-variance columns.
    sd = X.std(axis=0, ddof=1)
    cols = sd[sd > 0].index.tolist()
    X = X[cols]

    if transform == "zscore":
        X_mat = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1).replace(0, 1.0)
        X_mat = X_mat.fillna(0.0).values

    elif transform == "clr_zscore":
        X_nonneg = X.clip(lower=0)
        positive = X_nonneg.values[X_nonneg.values > 0]
        pc = np.nanpercentile(positive, 1) if len(positive) else 1e-6
        pc = max(pc, 1e-6)

        logX = np.log(X_nonneg + pc)
        clr = logX.sub(logX.mean(axis=1), axis=0)
        X_mat = (clr - clr.mean(axis=0)) / clr.std(axis=0, ddof=1).replace(0, 1.0)
        X_mat = X_mat.fillna(0.0).values

    else:
        raise ValueError("transform must be 'zscore' or 'clr_zscore'.")

    return X_mat.astype(np.float32), cols


def clean_day(x):
    s = str(x).lower().strip()
    s = s.replace("day", "d").replace(" ", "")

    # Numeric days.
    if s in ["0", "0.0"]:
        return "d0"
    if s in ["1", "1.0"]:
        return "d1"
    if s in ["2", "2.0"]:
        return "d2"
    if s in ["3", "3.0"]:
        return "d3"

    return s


def prepare_base_inputs(adata):
    """
    Prepare:
      - obs with UMAP coordinates and guide load
      - target_assignment
      - cNMF KNN matrix
      - day mask
    """

    adata_use = adata.copy()

    if not adata_use.obs_names.is_unique:
        adata_use.obs_names_make_unique()

    adata_use.obs["day_clean"] = adata_use.obs[DAY_COL].map(clean_day)

    if CELLTYPE_COL in adata_use.obs.columns:
        adata_use.obs["cell_type_clean"] = adata_use.obs[CELLTYPE_COL].astype(str)
    else:
        adata_use.obs["cell_type_clean"] = "unknown"

    # Guide metadata and target assignments.
    guide_target = load_guide_target_from_metadata(adata_use)
    target_assignment, target_to_guide_indices = build_target_assignment(adata_use, guide_target)

    # Guide load summary.
    guide_assignment = get_guide_assignment_matrix(adata_use)
    adata_use.obs["n_total_guides"] = (
        np.asarray(guide_assignment.sum(axis=1)).ravel()
        if sparse.issparse(guide_assignment)
        else guide_assignment.sum(axis=1)
    ).astype(int)

    targeting_cols = [c for c in target_assignment.columns if c != CONTROL_TARGET]
    adata_use.obs["n_targeting_targets"] = target_assignment[targeting_cols].sum(axis=1).astype(int).values

    # cNMF KNN matrix.
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

    # Day filtering.
    if DAYS_TO_USE is not None:
        days_clean = [clean_day(d) for d in DAYS_TO_USE]
        day_mask = adata_use.obs["day_clean"].isin(days_clean).values
    else:
        day_mask = np.ones(adata_use.n_obs, dtype=bool)

    # Export background cells for plotting.
    background = adata_use.obs[[
        "UMAP1",
        "UMAP2",
        "day_clean",
        "cell_type_clean",
        "n_total_guides",
        "n_targeting_targets",
    ]].copy()

    background.index.name = "cell_id"
    background.reset_index().to_csv(OUTDIR / "all_days_background_cells_for_plotting.csv", index=False)

    pd.DataFrame({"program_column_used": program_cols_used}).to_csv(
        OUTDIR / "cnmf_programs_used_for_knn.csv",
        index=False,
    )

    print(f"\nUsing {len(program_cols_used)} cNMF programs for KNN.")
    print(f"KNN transform: {CNMF_TRANSFORM}")
    print(f"Mode: {MODE}")
    print(f"Restrict neighbors to same day: {RESTRICT_NEIGHBORS_TO_SAME_DAY}")
    print(f"Same-day baseline: {SAME_DAY_BASELINE}")
    print(f"Days used: {sorted(pd.Series(adata_use.obs.loc[day_mask, 'day_clean']).unique())}")

    return adata_use, target_assignment, X_knn, day_mask
```

## Code cell 27

```python
def get_masks_for_target(adata, target_assignment, target, mode=MODE, day_mask=None):
    if target not in target_assignment.columns:
        raise KeyError(f"Target {target} not found in target_assignment.")

    if day_mask is None:
        day_mask = np.ones(adata.n_obs, dtype=bool)

    target_mask_raw = target_assignment[target].values.astype(bool)
    nt_mask_raw = target_assignment[CONTROL_TARGET].values.astype(bool)

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


def fdr_bh(pvals):
    pvals = np.asarray(pvals, dtype=float)
    qvals = np.full_like(pvals, np.nan, dtype=float)

    valid = np.isfinite(pvals)
    if valid.sum() == 0:
        return qvals

    qvals[valid] = multipletests(pvals[valid], method="fdr_bh")[1]
    return qvals


def compute_carrier_seeded_enrichment(
    adata,
    target_assignment,
    X_knn,
    target,
    mode=MODE,
    k=K,
    day_mask=None,
    restrict_neighbors_to_same_day=RESTRICT_NEIGHBORS_TO_SAME_DAY,
    same_day_baseline=SAME_DAY_BASELINE,
):
    """
    For each carrier cell of `target`, use that cell as seed.
    Find kNN in cNMF space.
    Compute local enrichment of target carriers vs NT controls.

    Positive score = seed is in a local cNMF neighborhood enriched for that target.
    Negative score = local depletion.
    """

    target_mask, nt_mask = get_masks_for_target(
        adata=adata,
        target_assignment=target_assignment,
        target=target,
        mode=mode,
        day_mask=day_mask,
    )

    seed_idx_all = np.where(target_mask)[0]

    if len(seed_idx_all) == 0:
        print(f"[skip] {target}: no carrier seeds.")
        return pd.DataFrame()

    day_vec = adata.obs["day_clean"].astype(str).values
    cell_type_vec = adata.obs["cell_type_clean"].astype(str).values

    rows = []
    neighbor_membership_rows = []

    if restrict_neighbors_to_same_day:
        # Fit separate KNN models per day so every seed gets same-day neighbors.
        days_for_target = sorted(pd.Series(day_vec[seed_idx_all]).unique())

        for day in days_for_target:
            day_cell_idx = np.where((day_vec == day) & day_mask)[0]
            seed_idx = seed_idx_all[day_vec[seed_idx_all] == day]

            if len(day_cell_idx) <= 1 or len(seed_idx) == 0:
                continue

            local_pos = {cell_idx: pos for pos, cell_idx in enumerate(day_cell_idx)}
            seed_local = np.array([local_pos[i] for i in seed_idx])

            X_day = X_knn[day_cell_idx]

            k_query = min(k + 1, len(day_cell_idx))

            nbrs = NearestNeighbors(n_neighbors=k_query, metric="euclidean")
            nbrs.fit(X_day)

            neigh_local = nbrs.kneighbors(X_day[seed_local], return_distance=False)

            for seed_cell_idx, neigh_l in zip(seed_idx, neigh_local):
                # Convert local day-neighbor positions to global indices.
                neigh = day_cell_idx[neigh_l]

                # Remove seed itself.
                neigh = neigh[neigh != seed_cell_idx]
                neigh = neigh[:k]

                rec, members = score_one_seed_neighborhood(
                    adata=adata,
                    target=target,
                    seed_cell_idx=seed_cell_idx,
                    neigh=neigh,
                    target_mask=target_mask,
                    nt_mask=nt_mask,
                    day_vec=day_vec,
                    cell_type_vec=cell_type_vec,
                    same_day_baseline=same_day_baseline,
                )

                rows.append(rec)

                if EXPORT_NEIGHBOR_MEMBERSHIP:
                    neighbor_membership_rows.extend(members)

    else:
        # Global KNN across all selected days.
        selected_idx = np.where(day_mask)[0]
        pos = {cell_idx: i for i, cell_idx in enumerate(selected_idx)}

        seed_idx = seed_idx_all
        seed_local = np.array([pos[i] for i in seed_idx if i in pos])
        seed_idx = np.array([i for i in seed_idx if i in pos])

        X_sel = X_knn[selected_idx]
        k_query = min(k + 1, len(selected_idx))

        nbrs = NearestNeighbors(n_neighbors=k_query, metric="euclidean")
        nbrs.fit(X_sel)

        neigh_local = nbrs.kneighbors(X_sel[seed_local], return_distance=False)

        for seed_cell_idx, neigh_l in zip(seed_idx, neigh_local):
            neigh = selected_idx[neigh_l]
            neigh = neigh[neigh != seed_cell_idx]
            neigh = neigh[:k]

            rec, members = score_one_seed_neighborhood(
                adata=adata,
                target=target,
                seed_cell_idx=seed_cell_idx,
                neigh=neigh,
                target_mask=target_mask,
                nt_mask=nt_mask,
                day_vec=day_vec,
                cell_type_vec=cell_type_vec,
                same_day_baseline=same_day_baseline,
            )

            rows.append(rec)

            if EXPORT_NEIGHBOR_MEMBERSHIP:
                neighbor_membership_rows.extend(members)

    out = pd.DataFrame(rows)

    if len(out) == 0:
        return out

    # FDR per target across all seed cells.
    out["fisher_fdr"] = fdr_bh(out["fisher_p"].values)

    # FDR per target-day.
    out["fisher_fdr_by_day"] = np.nan
    for day, sub in out.groupby("day"):
        idx = sub.index
        out.loc[idx, "fisher_fdr_by_day"] = fdr_bh(sub["fisher_p"].values)

    out["minus_log10_fdr"] = -np.log10(out["fisher_fdr"].clip(lower=1e-300))
    out["minus_log10_fdr_by_day"] = -np.log10(out["fisher_fdr_by_day"].clip(lower=1e-300))
    out["significant_fdr_0_05"] = out["fisher_fdr"] < 0.05
    out["significant_fdr_by_day_0_05"] = out["fisher_fdr_by_day"] < 0.05

    if EXPORT_NEIGHBOR_MEMBERSHIP:
        members_df = pd.DataFrame(neighbor_membership_rows)
        members_df.to_csv(
            OUTDIR / f"neighbor_membership_{target}_{mode}.csv",
            index=False,
        )

    return out


def score_one_seed_neighborhood(
    adata,
    target,
    seed_cell_idx,
    neigh,
    target_mask,
    nt_mask,
    day_vec,
    cell_type_vec,
    same_day_baseline=True,
):
    seed_day = day_vec[seed_cell_idx]

    if same_day_baseline:
        bg_mask = day_vec == seed_day
    else:
        bg_mask = np.ones(len(day_vec), dtype=bool)

    # Exclude seed from target background.
    bg_target_total = int(np.sum(target_mask & bg_mask)) - 1
    bg_target_total = max(bg_target_total, 0)

    bg_nt_total = int(np.sum(nt_mask & bg_mask))

    local_target_n = int(np.sum(target_mask[neigh]))
    local_nt_n = int(np.sum(nt_mask[neigh]))

    local_ratio = (local_target_n + PSEUDOCOUNT) / (local_nt_n + PSEUDOCOUNT)
    bg_ratio = (bg_target_total + PSEUDOCOUNT) / (bg_nt_total + PSEUDOCOUNT)

    log2_local_enrichment = np.log2(local_ratio / bg_ratio)

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

    n_target_nt = local_target_n + local_nt_n

    rec = {
        "target": target,
        "seed_cell_idx": int(seed_cell_idx),
        "seed_cell_id": str(adata.obs_names[seed_cell_idx]),
        "day": seed_day,
        "cell_type": cell_type_vec[seed_cell_idx],
        "mode": MODE,
        "k_requested": int(K),
        "k_used": int(len(neigh)),
        "restrict_neighbors_to_same_day": bool(RESTRICT_NEIGHBORS_TO_SAME_DAY),
        "same_day_baseline": bool(same_day_baseline),
        "local_target_n": local_target_n,
        "local_nt_n": local_nt_n,
        "local_target_plus_nt_n": n_target_nt,
        "background_target_total": int(bg_target_total),
        "background_nt_total": int(bg_nt_total),
        "local_target_fraction_among_target_nt": (
            local_target_n / n_target_nt if n_target_nt > 0 else np.nan
        ),
        "background_target_fraction_among_target_nt": (
            bg_target_total / (bg_target_total + bg_nt_total)
            if (bg_target_total + bg_nt_total) > 0 else np.nan
        ),
        "log2_local_enrichment": float(log2_local_enrichment),
        "fisher_p": float(fisher_p) if np.isfinite(fisher_p) else np.nan,
        "UMAP1": float(adata.obs["UMAP1"].iloc[seed_cell_idx]),
        "UMAP2": float(adata.obs["UMAP2"].iloc[seed_cell_idx]),
        "n_total_guides": int(adata.obs["n_total_guides"].iloc[seed_cell_idx]),
        "n_targeting_targets": int(adata.obs["n_targeting_targets"].iloc[seed_cell_idx]),
    }

    members = []
    if EXPORT_NEIGHBOR_MEMBERSHIP:
        for rank, cell_idx in enumerate(neigh):
            members.append({
                "target": target,
                "seed_cell_id": str(adata.obs_names[seed_cell_idx]),
                "seed_cell_idx": int(seed_cell_idx),
                "neighbor_rank": int(rank),
                "neighbor_cell_idx": int(cell_idx),
                "neighbor_cell_id": str(adata.obs_names[cell_idx]),
                "seed_day": seed_day,
                "neighbor_day": day_vec[cell_idx],
                "neighbor_cell_type": cell_type_vec[cell_idx],
                "neighbor_is_target_carrier": bool(target_mask[cell_idx]),
                "neighbor_is_nt_control": bool(nt_mask[cell_idx]),
            })

    return rec, members
```

## Code cell 28

```python
def run_carrier_seeded_Hippo_all_days(adata, targets=HIPPO_TARGETS):
    """
    Main export function.

    Outputs:
      - all_days_background_cells_for_plotting.csv
      - carrier_seeded_scores_<TARGET>_<MODE>.csv
      - combined_carrier_seeded_scores_<MODE>.csv
      - summary_by_target_day_<MODE>.csv
      - summary_by_target_<MODE>.csv
    """

    adata_use, target_assignment, X_knn, day_mask = prepare_base_inputs(adata)

    present_targets = []
    missing_targets = []

    for t in targets:
        matches = [c for c in target_assignment.columns if c.upper() == t.upper()]
        if len(matches):
            present_targets.append(matches[0])
        else:
            missing_targets.append(t)

    print("\nPresent Hippo targets:", present_targets)
    if missing_targets:
        print("Missing Hippo targets:", missing_targets)

    all_scores = []

    for target in present_targets:
        print("\n" + "=" * 80)
        print(f"Computing carrier-seeded cNMF-kNN enrichment for {target}")

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
            out_file = OUTDIR / f"carrier_seeded_scores_{target}_{MODE}.csv"
            scores.to_csv(out_file, index=False)
            print(f"Saved: {out_file}")
            all_scores.append(scores)

    if not all_scores:
        raise RuntimeError("No carrier-seeded scores were generated.")

    combined = pd.concat(all_scores, ignore_index=True)

    combined_file = OUTDIR / f"combined_carrier_seeded_scores_{MODE}.csv"
    combined.to_csv(combined_file, index=False)

    # Summary by target and day.
    summary_by_target_day = combined.groupby(["target", "day"]).agg(
        n_seed_cells=("seed_cell_id", "count"),
        median_log2_enrichment=("log2_local_enrichment", "median"),
        mean_log2_enrichment=("log2_local_enrichment", "mean"),
        frac_positive=("log2_local_enrichment", lambda x: np.mean(x > 0)),
        n_fdr05=("fisher_fdr", lambda x: int(np.sum(x < 0.05))),
        n_fdr10=("fisher_fdr", lambda x: int(np.sum(x < 0.10))),
        n_fdr_by_day05=("fisher_fdr_by_day", lambda x: int(np.sum(x < 0.05))),
        n_fdr_by_day10=("fisher_fdr_by_day", lambda x: int(np.sum(x < 0.10))),
        median_k_used=("k_used", "median"),
        median_local_target_n=("local_target_n", "median"),
        median_local_nt_n=("local_nt_n", "median"),
    ).reset_index()

    summary_by_target = combined.groupby("target").agg(
        n_seed_cells=("seed_cell_id", "count"),
        median_log2_enrichment=("log2_local_enrichment", "median"),
        mean_log2_enrichment=("log2_local_enrichment", "mean"),
        frac_positive=("log2_local_enrichment", lambda x: np.mean(x > 0)),
        n_fdr05=("fisher_fdr", lambda x: int(np.sum(x < 0.05))),
        n_fdr10=("fisher_fdr", lambda x: int(np.sum(x < 0.10))),
        median_k_used=("k_used", "median"),
    ).reset_index()

    summary_by_target_day_file = OUTDIR / f"summary_by_target_day_{MODE}.csv"
    summary_by_target_file = OUTDIR / f"summary_by_target_{MODE}.csv"

    summary_by_target_day.to_csv(summary_by_target_day_file, index=False)
    summary_by_target.to_csv(summary_by_target_file, index=False)

    print("\nSaved combined output:")
    print(combined_file)

    print("\nSummary by target/day:")
    print(summary_by_target_day.to_string(index=False))

    print("\nSummary by target:")
    print(summary_by_target.to_string(index=False))

    # Zip package.
    zip_path = Path(f"{OUTDIR}.zip")
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for fn in OUTDIR.iterdir():
            z.write(fn, arcname=fn.name)

    print(f"\nZIP package: {zip_path}")

    return {
        "adata_use": adata_use,
        "target_assignment": target_assignment,
        "scores": combined,
        "summary_by_target_day": summary_by_target_day,
        "summary_by_target": summary_by_target,
    }
```

## Code cell 29

```python
adata = adata_int.copy()
```

## Code cell 30

```python
results = run_carrier_seeded_Hippo_all_days(adata)

scores = results["scores"]
summary_by_target_day = results["summary_by_target_day"]

scores.head()
```

### Text output (stream)

```text

Guide metadata mapping QC
adata guide columns:            1,922
metadata guide rows:            1,926
metadata rows matched to adata: 1,922
expanded guide-target rows:     1,922
unique mapped targets:          299

Match methods:
match_method
exact_guide_id_match    1922

Example guide-target rows:
 guide_index                 guide_name target_name raw_target_entry         match_method
           0   ACAA1_+_38178488.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           1   ACAA1_+_38178517.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           2   ACAA1_+_38178559.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           3   ACAA1_+_38178570.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           4   ACAA1_+_38178620.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           5   ACAA1_-_38178575.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           6 ACTL7B_+_111618040.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
           7 ACTL7B_+_111618198.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
           8 ACTL7B_+_111618251.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
           9 ACTL7B_-_111617995.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
          10 ACTL7B_-_111618185.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
          11 ACTL7B_-_111618259.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
          12  ACVR1B_+_52345347.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          13  ACVR1B_+_52345561.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          14  ACVR1B_+_52345762.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          15  ACVR1B_-_52345530.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          16  ACVR1B_-_52345533.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          17  ACVR1B_-_52345539.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          18 ACVR2A_+_148602340.23-P1P2      ACVR2A           ACVR2A exact_guide_id_match
          19 ACVR2A_+_148602343.23-P1P2      ACVR2A           ACVR2A exact_guide_id_match

adata guide columns without metadata target: 0

Target assignment QC
Targets with masks: 299
Control carrier cells: 21,671

Key target guide/cell counts:
  NF2     : guides=  6, carrier_cells=  2354
  LATS2   : guides=  6, carrier_cells=  2695
  SAV1    : guides=  6, carrier_cells=  1595
  TAOK1   : guides=  6, carrier_cells=  2530
  AMOTL2  : guides=  6, carrier_cells=  2381
  KDR     : guides=  6, carrier_cells=  1798
  ETV2    : guides=  6, carrier_cells=  2275
  BMPR1A  : guides=  6, carrier_cells=  1643
  SMAD4   : guides=  6, carrier_cells=  3325

Using 50 cNMF programs for KNN.
KNN transform: zscore
Mode: all_carriers
Restrict neighbors to same day: True
Same-day baseline: True
Days used: ['d0', 'd1', 'd2', 'd3']

Present Hippo targets: ['NF2', 'LATS2', 'SAV1', 'TAOK1', 'AMOTL2']

================================================================================
Computing carrier-seeded cNMF-kNN enrichment for NF2
Saved: carrier_seeded_Hippo_all_days_exports/carrier_seeded_scores_NF2_all_carriers.csv

================================================================================
Computing carrier-seeded cNMF-kNN enrichment for LATS2
Saved: carrier_seeded_Hippo_all_days_exports/carrier_seeded_scores_LATS2_all_carriers.csv

================================================================================
Computing carrier-seeded cNMF-kNN enrichment for SAV1
Saved: carrier_seeded_Hippo_all_days_exports/carrier_seeded_scores_SAV1_all_carriers.csv

================================================================================
Computing carrier-seeded cNMF-kNN enrichment for TAOK1
Saved: carrier_seeded_Hippo_all_days_exports/carrier_seeded_scores_TAOK1_all_carriers.csv

================================================================================
Computing carrier-seeded cNMF-kNN enrichment for AMOTL2

```

### Text output (stream)

```text
Saved: carrier_seeded_Hippo_all_days_exports/carrier_seeded_scores_AMOTL2_all_carriers.csv

Saved combined output:
carrier_seeded_Hippo_all_days_exports/combined_carrier_seeded_scores_all_carriers.csv

Summary by target/day:
target day  n_seed_cells  median_log2_enrichment  mean_log2_enrichment  frac_positive  n_fdr05  n_fdr10  n_fdr_by_day05  n_fdr_by_day10  median_k_used  median_local_target_n  median_local_nt_n
AMOTL2  d0           740                0.116951              0.070449       0.575676        4        9               1               1          300.0                    6.0               58.0
AMOTL2  d1           328                0.084979              0.038545       0.570122        2        3               0               0          300.0                    7.0               69.0
AMOTL2  d2           503                0.256421              0.170287       0.630219       34       58              23              54          300.0                   10.0               75.0
AMOTL2  d3           810                0.324446              0.345291       0.622222      216      240             244             275          300.0                    9.0               63.0
 LATS2  d0           810                0.060711             -0.011305       0.549383       10       16               0               0          300.0                    7.0               58.0
 LATS2  d1           373                0.124966              0.057033       0.584450        4        8               0               0          300.0                    8.0               68.0
 LATS2  d2           516                0.164829              0.061204       0.573643       28       44               5              14          300.0                   10.0               75.0
 LATS2  d3           996                0.616858              0.495177       0.640562      392      446             462             511          300.0                   14.0               59.0
   NF2  d0           718                0.088011              0.024338       0.543175        5       16               0               0          300.0                    6.0               58.0
   NF2  d1           339                0.319628              0.237187       0.696165       10       22               2               8          300.0                    8.0               66.0
   NF2  d2           474                0.298443              0.165061       0.647679       15       31               0               0          300.0                   10.0               75.0
   NF2  d3           823                0.785058              0.650888       0.674362      343      378             388             446          300.0                   14.0               63.0
  SAV1  d0           478                0.052470             -0.055157       0.525105        0        0               0               0          300.0                    4.0               59.0
  SAV1  d1           226                0.025934              0.020881       0.504425        0        0               0               0          300.0                    4.0               67.0
  SAV1  d2           323                0.303603              0.134657       0.634675        1        4               0               0          300.0                    6.0               70.0
  SAV1  d3           568                0.255421              0.156503       0.621479       33       47              55              68          300.0                    5.0               61.5
 TAOK1  d0           788                0.085548              0.010762       0.562183        3        7               0               0          300.0                    7.0               58.0
 TAOK1  d1           384                0.101452              0.009965       0.580729        0        0               0               0          300.0                    8.0               67.0
 TAOK1  d2           496                0.326294              0.180419       0.693548        5       12               0               0          300.0                   10.0               71.5
 TAOK1  d3           862                0.312212              0.203276       0.613689      112      148             167             214          300.0                    8.0               61.0

Summary by target:
target  n_seed_cells  median_log2_enrichment  mean_log2_enrichment  frac_positive  n_fdr05  n_fdr10  median_k_used
AMOTL2          2381                0.195637              0.180645       0.602268      256      310          300.0
 LATS2          2695                0.186223              0.199219       0.592579      434      514          300.0
   NF2          2354                0.298443              0.302379       0.632116      373      447          300.0
  SAV1          1595                0.158926              0.069431       0.578683       34       51          300.0
 TAOK1          2530                0.194211              0.109494       0.608300      120      167          300.0

ZIP package: carrier_seeded_Hippo_all_days_exports.zip

```

### Text output (text/plain)

```text
  target  seed_cell_idx                seed_cell_id day   cell_type  \
0    NF2             55  TGGTGGTAAAGGACACAAACGATA_2  d0        iPSC   
1    NF2             91  CCATCCTCAAACATCGTGTCTGAA_1  d0        iPSC   
2    NF2            150  GAGTTAGCGGAGAACAACTCGTAA_3  d0  Neuro-iPSC   
3    NF2            173  AGAGTCAAAGCCATGCACTCGTAA_1  d0        iPSC   
4    NF2            184  AGTCACTAACAGCAGAATCGCATA_3  d0        iPSC   

           mode  k_requested  k_used  restrict_neighbors_to_same_day  \
0  all_carriers          300     300                            True   
1  all_carriers          300     300                            True   
2  all_carriers          300     300                            True   
3  all_carriers          300     300                            True   
4  all_carriers          300     300                            True   

   same_day_baseline  ...     UMAP1      UMAP2  n_total_guides  \
0               True  ... -4.043896   6.194844               5   
1               True  ... -2.959681   8.991536              18   
2               True  ... -3.995974  12.443013              14   
3               True  ... -0.499485   9.570541               7   
4               True  ... -1.690850   8.355458               8   

   n_targeting_targets  fisher_fdr  fisher_fdr_by_day  minus_log10_fdr  \
0                    5    0.728483                1.0         0.137580   
1                   18    1.000000                1.0        -0.000000   
2                   10    0.815443                1.0         0.088606   
3                    7    0.635627                1.0         0.196797   
4                    8    1.000000                1.0        -0.000000   

   minus_log10_fdr_by_day  significant_fdr_0_05  significant_fdr_by_day_0_05  
0                    -0.0                 False                        False  
1                    -0.0                 False                        False  
2                    -0.0                 False                        False  
3                    -0.0                 False                        False  
4                    -0.0                 False                        False  

[5 rows x 29 columns]
```

## Code cell 31

```python
import os, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

from plotnine import (
    ggplot, aes, geom_point, scale_fill_gradient2, scale_size_continuous,
    coord_equal, labs, theme_bw, theme, element_text, element_blank, element_rect
)

BASE = Path("/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports")
OUTDIR = BASE / "carrier_seeded_single_gene_allDays_UMAP_plotnine"
OUTDIR.mkdir(parents=True, exist_ok=True)

BACKGROUND_CSV = BASE / "all_days_background_cells_for_plotting.csv"
GENES = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]

# Plot settings
LOG2FC_CLIP = 3.0
MAX_BACKGROUND_CELLS = 60000  # stratified subsample for clean rendering; all days are represented
RANDOM_SEED = 17

bg = pd.read_csv(BACKGROUND_CSV)
bg = bg[["UMAP1", "UMAP2", "day_clean"]].copy()

frames = []
for gene in GENES:
    f = BASE / f"carrier_seeded_scores_{gene}_all_carriers.csv"
    df = pd.read_csv(f)
    df["target"] = gene
    frames.append(df)
scores = pd.concat(frames, ignore_index=True)

scores["log2OR_plot"] = scores["log2_local_enrichment"].clip(-LOG2FC_CLIP, LOG2FC_CLIP)
scores["is_sig"] = scores["significant_fdr_0_05"].fillna(False).astype(bool)

# TRUE carrier fraction in neighborhood:
# local target carriers divided by all kNN cells in that seed-cell neighborhood.
scores["carrier_fraction_pct"] = 100 * scores["local_target_n"] / scores["k_used"].replace(0, np.nan)
scores["carrier_fraction_pct"] = scores["carrier_fraction_pct"].clip(lower=0)

# Stratified background sample by day for speed and readability.
rng = np.random.default_rng(RANDOM_SEED)
if len(bg) > MAX_BACKGROUND_CELLS:
    bg_parts = []
    n_days = max(1, bg["day_clean"].nunique())
    n_per_day = MAX_BACKGROUND_CELLS // n_days
    for day, sub in bg.groupby("day_clean"):
        n = min(len(sub), n_per_day)
        idx = rng.choice(sub.index.values, size=n, replace=False) if len(sub) > n else sub.index.values
        bg_parts.append(bg.loc[idx])
    bg_plot = pd.concat(bg_parts, ignore_index=True)
else:
    bg_plot = bg.copy()

# Global limits so each single-gene panel is directly comparable.
xmin = min(bg["UMAP1"].min(), scores["UMAP1"].min())
xmax = max(bg["UMAP1"].max(), scores["UMAP1"].max())
ymin = min(bg["UMAP2"].min(), scores["UMAP2"].min())
ymax = max(bg["UMAP2"].max(), scores["UMAP2"].max())

pad_x = (xmax - xmin) * 0.03
pad_y = (ymax - ymin) * 0.03
XLIM = (xmin - pad_x, xmax + pad_x)
YLIM = (ymin - pad_y, ymax + pad_y)

size_upper = max(30, float(scores["carrier_fraction_pct"].quantile(0.995)))

def make_plot(gene):
    df = scores[scores["target"] == gene].copy()
    nonsig = df[~df["is_sig"]]
    sig = df[df["is_sig"]]

    p = (
        ggplot()
        + geom_point(
            bg_plot,
            aes("UMAP1", "UMAP2"),
            color="#D0D0D0",
            size=0.22,
            alpha=0.18,
            raster=True,
            show_legend=False,
        )
        + geom_point(
            nonsig,
            aes("UMAP1", "UMAP2", fill="log2OR_plot", size="carrier_fraction_pct"),
            shape="o",
            color="#666666",
            stroke=0.30,
            alpha=0.88,
        )
        + geom_point(
            sig,
            aes("UMAP1", "UMAP2", fill="log2OR_plot", size="carrier_fraction_pct"),
            shape="o",
            color="black",
            stroke=0.75,
            alpha=0.96,
        )
        + scale_fill_gradient2(
            low="#2166AC",
            mid="#F7F7F7",
            high="#B2182B",
            midpoint=0,
            limits=(-LOG2FC_CLIP, LOG2FC_CLIP),
            breaks=[-3, -2, -1, 0, 1, 2, 3],
            name="Carrier-seeded\nenrichment\nlog2OR",
        )
        + scale_size_continuous(
            range=(0.4, 3.5),
            breaks=[2, 5, 10, 20, 30],
            limits=(0, size_upper),
            name="Carrier fraction\nin neighborhood (%)",
        )
        + coord_equal(xlim=XLIM, ylim=YLIM, expand=False)
        + labs(
            title=gene,
            subtitle="All D0-D3 cells in one shared UMAP; circles are carrier seed cells; black outline = FDR < 0.05",
            x="UMAP1",
            y="UMAP2",
        )
        + theme_bw(base_size=12)
        + theme(
            figure_size=(7.2, 6.3),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(color="#7f7f7f", fill=None, size=0.8),
            axis_title=element_text(size=13),
            axis_text=element_text(size=10),
            plot_title=element_text(size=22, weight="bold", ha="left"),
            plot_subtitle=element_text(size=10.5),
            legend_title=element_text(size=11, weight="bold"),
            legend_text=element_text(size=10),
            legend_key=element_blank(),
        )
    )
    return p

created = []
for gene in GENES:
    p = make_plot(gene)
    for ext in ["png", "pdf"]:
        fp = OUTDIR / f"carrier_seeded_{gene}_allDays_singleUMAP_plotnine.{ext}"
        if ext == "png":
            p.save(fp, width=7.2, height=6.3, dpi=450, verbose=False)
        else:
            p.save(fp, width=7.2, height=6.3, verbose=False)
        created.append(fp)

# Export plotting sources and summary.
scores.to_csv(OUTDIR / "source_carrier_seeded_scores_for_plotting.csv", index=False)
bg_plot.to_csv(OUTDIR / "source_background_cells_sampled_for_plotting.csv", index=False)

summary = scores.groupby("target", observed=True).agg(
    n_seed_cells=("seed_cell_id", "count"),
    median_log2_enrichment=("log2_local_enrichment", "median"),
    mean_log2_enrichment=("log2_local_enrichment", "mean"),
    median_carrier_fraction_pct=("carrier_fraction_pct", "median"),
    n_fdr05=("is_sig", "sum"),
).reset_index()
summary.to_csv(OUTDIR / "carrier_seeded_plot_summary.csv", index=False)

# Zip package.
zip_path = BASE / "carrier_seeded_single_gene_allDays_UMAP_plotnine_package.zip"
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
    for fp in OUTDIR.iterdir():
        z.write(fp, arcname=fp.name)

print("Created files:")
for fp in created:
    print(fp)
print("ZIP:", zip_path)
print(summary.to_string(index=False))

```

### Text output (stream)

```text
Created files:
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine/carrier_seeded_NF2_allDays_singleUMAP_plotnine.png
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine/carrier_seeded_NF2_allDays_singleUMAP_plotnine.pdf
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine/carrier_seeded_LATS2_allDays_singleUMAP_plotnine.png
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine/carrier_seeded_LATS2_allDays_singleUMAP_plotnine.pdf
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine/carrier_seeded_SAV1_allDays_singleUMAP_plotnine.png
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine/carrier_seeded_SAV1_allDays_singleUMAP_plotnine.pdf
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine/carrier_seeded_TAOK1_allDays_singleUMAP_plotnine.png
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine/carrier_seeded_TAOK1_allDays_singleUMAP_plotnine.pdf
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine/carrier_seeded_AMOTL2_allDays_singleUMAP_plotnine.png
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine/carrier_seeded_AMOTL2_allDays_singleUMAP_plotnine.pdf
ZIP: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_package.zip
target  n_seed_cells  median_log2_enrichment  mean_log2_enrichment  median_carrier_fraction_pct  n_fdr05
AMOTL2          2381                0.195637              0.180645                     2.666667      256
 LATS2          2695                0.186223              0.199219                     2.666667      434
   NF2          2354                0.298443              0.302379                     2.666667      373
  SAV1          1595                0.158926              0.069431                     1.666667       34
 TAOK1          2530                0.194211              0.109494                     2.666667      120

```

## Markdown cell 32

##### Replot with just color

## Code cell 33

```python
# Cell 1
import os
from pathlib import Path

import numpy as np
import pandas as pd

from plotnine import (
    ggplot, aes, geom_point,
    scale_color_gradient2,
    coord_equal,
    labs,
    theme_bw, theme,
    element_text, element_blank, element_rect
)

BASE = Path("carrier_seeded_Hippo_all_days_exports")
OUTDIR = BASE / "carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor"
OUTDIR.mkdir(parents=True, exist_ok=True)

BACKGROUND_CSV = BASE / "all_days_background_cells_for_plotting.csv"
GENES = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]

LOG2OR_CLIP = 2.0
MAX_BACKGROUND_CELLS = 60000
RANDOM_SEED = 17

FOREGROUND_POINT_SIZE = 1.0
BACKGROUND_POINT_SIZE = 0.1
```

## Code cell 34

```python
# Cell 2
# Load background cells
bg = pd.read_csv(BACKGROUND_CSV).copy()

# Keep only needed columns
bg = bg[["UMAP1", "UMAP2"] + ([c for c in ["day_clean", "day"] if c in bg.columns])].copy()

# Load all gene-specific carrier-seeded score files
frames = []
for gene in GENES:
    f = BASE / f"carrier_seeded_scores_{gene}_all_carriers.csv"
    df = pd.read_csv(f).copy()
    df["target"] = gene
    frames.append(df)

scores = pd.concat(frames, ignore_index=True)

# Clip the color scale at [-2, 2]
scores["log2OR_plot"] = scores["log2_local_enrichment"].clip(-LOG2OR_CLIP, LOG2OR_CLIP)

# Background subsampling for speed / clarity
rng = np.random.default_rng(RANDOM_SEED)

if len(bg) > MAX_BACKGROUND_CELLS:
    if "day_clean" in bg.columns:
        bg_parts = []
        n_days = max(1, bg["day_clean"].nunique())
        n_per_day = MAX_BACKGROUND_CELLS // n_days

        for day, sub in bg.groupby("day_clean"):
            n = min(len(sub), n_per_day)
            if len(sub) > n:
                idx = rng.choice(sub.index.values, size=n, replace=False)
                bg_parts.append(sub.loc[idx])
            else:
                bg_parts.append(sub)

        bg_plot = pd.concat(bg_parts, ignore_index=True)

    elif "day" in bg.columns:
        bg_parts = []
        n_days = max(1, bg["day"].nunique())
        n_per_day = MAX_BACKGROUND_CELLS // n_days

        for day, sub in bg.groupby("day"):
            n = min(len(sub), n_per_day)
            if len(sub) > n:
                idx = rng.choice(sub.index.values, size=n, replace=False)
                bg_parts.append(sub.loc[idx])
            else:
                bg_parts.append(sub)

        bg_plot = pd.concat(bg_parts, ignore_index=True)

    else:
        idx = rng.choice(bg.index.values, size=MAX_BACKGROUND_CELLS, replace=False)
        bg_plot = bg.loc[idx].copy()
else:
    bg_plot = bg.copy()

# Common axis limits for all genes
xmin = min(bg["UMAP1"].min(), scores["UMAP1"].min())
xmax = max(bg["UMAP1"].max(), scores["UMAP1"].max())
ymin = min(bg["UMAP2"].min(), scores["UMAP2"].min())
ymax = max(bg["UMAP2"].max(), scores["UMAP2"].max())

pad_x = (xmax - xmin) * 0.03
pad_y = (ymax - ymin) * 0.03

XLIM = (xmin - pad_x, xmax + pad_x)
YLIM = (ymin - pad_y, ymax + pad_y)

print(bg_plot.shape)
print(scores.shape)
scores.head()
```

### Text output (stream)

```text
(60000, 3)
(11555, 30)

```

### Text output (text/plain)

```text
  target  seed_cell_idx                seed_cell_id day   cell_type  \
0    NF2             55  TGGTGGTAAAGGACACAAACGATA_2  d0        iPSC   
1    NF2             91  CCATCCTCAAACATCGTGTCTGAA_1  d0        iPSC   
2    NF2            150  GAGTTAGCGGAGAACAACTCGTAA_3  d0  Neuro-iPSC   
3    NF2            173  AGAGTCAAAGCCATGCACTCGTAA_1  d0        iPSC   
4    NF2            184  AGTCACTAACAGCAGAATCGCATA_3  d0        iPSC   

           mode  k_requested  k_used  restrict_neighbors_to_same_day  \
0  all_carriers          300     300                            True   
1  all_carriers          300     300                            True   
2  all_carriers          300     300                            True   
3  all_carriers          300     300                            True   
4  all_carriers          300     300                            True   

   same_day_baseline  ...      UMAP2  n_total_guides  n_targeting_targets  \
0               True  ...   6.194844               5                    5   
1               True  ...   8.991536              18                   18   
2               True  ...  12.443013              14                   10   
3               True  ...   9.570541               7                    7   
4               True  ...   8.355458               8                    8   

   fisher_fdr  fisher_fdr_by_day  minus_log10_fdr  minus_log10_fdr_by_day  \
0    0.728483                1.0         0.137580                    -0.0   
1    1.000000                1.0        -0.000000                    -0.0   
2    0.815443                1.0         0.088606                    -0.0   
3    0.635627                1.0         0.196797                    -0.0   
4    1.000000                1.0        -0.000000                    -0.0   

   significant_fdr_0_05  significant_fdr_by_day_0_05  log2OR_plot  
0                 False                        False     0.500805  
1                 False                        False    -0.092853  
2                 False                        False     0.363991  
3                 False                        False    -0.689966  
4                 False                        False    -0.115396  

[5 rows x 30 columns]
```

## Code cell 35

```python
# Cell 1 — Illustrator-editable SVG/PDF settings, with automatic font detection

import matplotlib as mpl
from matplotlib import font_manager
from pathlib import Path
import logging

# Silence repeated matplotlib font-manager messages.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

# Detect available fonts on your system.
available_fonts = {f.name for f in font_manager.fontManager.ttflist}

# Prefer Illustrator/Nature-friendly sans-serif fonts if installed.
FONT_CANDIDATES = [
    "Arial",
    "Helvetica",
    "Liberation Sans",
    "Nimbus Sans",
    "DejaVu Sans",
]

AI_FONT = None
for font in FONT_CANDIDATES:
    if font in available_fonts:
        AI_FONT = font
        break

if AI_FONT is None:
    AI_FONT = "DejaVu Sans"

print("Using font:", AI_FONT)

mpl.rcParams.update({
    # Keep SVG text as editable text, not paths.
    "svg.fonttype": "none",

    # Keep PDF text as editable TrueType text.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,

    # Avoid LaTeX rendering, which can convert text into paths.
    "text.usetex": False,

    # Use the detected font.
    "font.family": "sans-serif",
    "font.sans-serif": [AI_FONT, "DejaVu Sans"],

    # Cleaner vector output.
    "savefig.transparent": False,
})

OUTDIR = Path(OUTDIR)
OUTDIR.mkdir(parents=True, exist_ok=True)

# If Illustrator becomes slow, set this to True.
# False = all points are vector/editable.
# True  = gray background points rasterized, but text/axes/foreground remain editable.
RASTERIZE_BACKGROUND_FOR_SPEED = True
```

### Text output (stream)

```text
Using font: DejaVu Sans

```

## Code cell 36

```python
# Cell 2 — Plot function with Illustrator-friendly text settings

from plotnine import (
    ggplot, aes,
    geom_point,
    scale_color_gradient2,
    coord_equal,
    labs,
    theme_bw, theme,
    element_blank, element_rect, element_text
)

def _geom_point_optional_raster(data, mapping, raster=False, **kwargs):
    """
    Helper for plotnine versions that may or may not support raster=True.
    If raster=True is unsupported, it falls back to normal vector points.
    """
    if raster:
        try:
            return geom_point(data, mapping, raster=True, **kwargs)
        except TypeError:
            return geom_point(data, mapping, **kwargs)

    return geom_point(data, mapping, **kwargs)


def make_plot_simple_color(gene):
    df = scores.loc[scores["target"].astype(str) == str(gene)].copy()

    if df.empty:
        print(f"[warning] no rows found for target: {gene}")

    p = (
        ggplot()

        # Background all-day cells.
        # Keep raster=False for fully vector/editable SVG.
        # Use RASTERIZE_BACKGROUND_FOR_SPEED=True only if Illustrator becomes too slow.
        + _geom_point_optional_raster(
            bg_plot,
            aes("UMAP1", "UMAP2"),
            raster=RASTERIZE_BACKGROUND_FOR_SPEED,
            color="#D0D0D0",
            size=BACKGROUND_POINT_SIZE,
            alpha=0.18,
            show_legend=False
        )

        # Carrier-seeded neighborhoods for one gene.
        # These remain vector/editable.
        + geom_point(
            df,
            aes("UMAP1", "UMAP2", color="log2OR_plot"),
            size=FOREGROUND_POINT_SIZE,
            alpha=0.95
        )

        + scale_color_gradient2(
            low="#053061",
            mid="#F7F7F7",
            high="#67001F",
            midpoint=0,
            limits=(-LOG2OR_CLIP, LOG2OR_CLIP),
            breaks=[-2, -1, 0, 1, 2],
            name="Carrier-seeded\nenrichment\nlog2OR"
        )

        + coord_equal(xlim=XLIM, ylim=YLIM, expand=False)

        + labs(
            title=gene,
            subtitle="All D0-D3 cells in one shared UMAP; color = local enrichment vs non-targeting controls",
            x="single-cell state embedding 1",
            y="single-cell state embedding 2"
        )

        + theme_bw(base_size=12, base_family=AI_FONT)

        + theme(
            figure_size=(7.2, 6.3),

            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            panel_border=element_rect(color="#7f7f7f", fill=None, size=0.8),

            axis_title=element_text(size=13, family=AI_FONT),
            axis_text=element_text(size=10, family=AI_FONT),

            plot_title=element_text(
                size=22,
                weight="bold",
                ha="left",
                family=AI_FONT
            ),
            plot_subtitle=element_text(
                size=10.5,
                family=AI_FONT
            ),

            legend_title=element_text(
                size=11,
                weight="bold",
                family=AI_FONT
            ),
            legend_text=element_text(
                size=10,
                family=AI_FONT
            ),
            legend_key=element_blank()
        )
    )

    return p
```

## Code cell 37

```python
# Cell 3 — Save one plot per gene

for gene in GENES:
    p = make_plot_simple_color(gene)

    png_file = OUTDIR / f"carrier_seeded_{gene}_allDays_singleUMAP_simpleColor_plotnine.png"
    pdf_file = OUTDIR / f"carrier_seeded_{gene}_allDays_singleUMAP_simpleColor_plotnine_AI_editable.pdf"
    svg_file = OUTDIR / f"carrier_seeded_{gene}_allDays_singleUMAP_simpleColor_plotnine_AI_editable.svg"

    # PNG preview.
    p.save(
        str(png_file),
        width=7.2,
        height=6.3,
        dpi=450,
        verbose=False
    )

    # Illustrator-friendly vector outputs.
    # SVG: text remains editable because svg.fonttype = "none".
    # PDF: text remains editable TrueType because pdf.fonttype = 42.
    p.save(
        str(pdf_file),
        width=7.2,
        height=6.3,
        verbose=False
    )

    p.save(
        str(svg_file),
        width=7.2,
        height=6.3,
        verbose=False
    )

    print("saved:", gene)
    print("  PNG:", png_file)
    print("  PDF:", pdf_file)
    print("  SVG:", svg_file)
```

### Text output (stream)

```text
saved: NF2
  PNG: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_NF2_allDays_singleUMAP_simpleColor_plotnine.png
  PDF: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_NF2_allDays_singleUMAP_simpleColor_plotnine_AI_editable.pdf
  SVG: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_NF2_allDays_singleUMAP_simpleColor_plotnine_AI_editable.svg
saved: LATS2
  PNG: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_LATS2_allDays_singleUMAP_simpleColor_plotnine.png
  PDF: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_LATS2_allDays_singleUMAP_simpleColor_plotnine_AI_editable.pdf
  SVG: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_LATS2_allDays_singleUMAP_simpleColor_plotnine_AI_editable.svg
saved: SAV1
  PNG: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_SAV1_allDays_singleUMAP_simpleColor_plotnine.png
  PDF: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_SAV1_allDays_singleUMAP_simpleColor_plotnine_AI_editable.pdf
  SVG: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_SAV1_allDays_singleUMAP_simpleColor_plotnine_AI_editable.svg
saved: TAOK1
  PNG: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_TAOK1_allDays_singleUMAP_simpleColor_plotnine.png
  PDF: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_TAOK1_allDays_singleUMAP_simpleColor_plotnine_AI_editable.pdf
  SVG: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_TAOK1_allDays_singleUMAP_simpleColor_plotnine_AI_editable.svg
saved: AMOTL2
  PNG: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_AMOTL2_allDays_singleUMAP_simpleColor_plotnine.png
  PDF: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_AMOTL2_allDays_singleUMAP_simpleColor_plotnine_AI_editable.pdf
  SVG: carrier_seeded_Hippo_all_days_exports/carrier_seeded_single_gene_allDays_UMAP_plotnine_simpleColor/carrier_seeded_AMOTL2_allDays_singleUMAP_simpleColor_plotnine_AI_editable.svg

```

## Code cell 38

```python
# Optional Cell 5: faceted overview
from plotnine import facet_wrap

scores["target"] = pd.Categorical(scores["target"], categories=GENES, ordered=True)

p_facet = (
    ggplot()
    + geom_point(
        bg_plot,
        aes("UMAP1", "UMAP2"),
        color="#D0D0D0",
        size=0.13,
        alpha=0.12,
        show_legend=False
    )
    + geom_point(
        scores,
        aes("UMAP1", "UMAP2", color="log2OR_plot"),
        size=1.7,
        alpha=0.95
    )
    + facet_wrap("~target", ncol=3)
    + scale_color_gradient2(
        low="#053061",
        mid="#F7F7F7",
        high="#67001F",
        midpoint=0,
        limits=(-LOG2OR_CLIP, LOG2OR_CLIP),
        breaks=[-2, -1, 0, 1, 2],
        name="Carrier-seeded\nenrichment\nlog2OR"
    )
    + coord_equal(xlim=XLIM, ylim=YLIM, expand=False)
    + labs(
        title="Carrier-seeded Hippo-gene enrichment across the all-day single-cell state UMAP",
        subtitle="One UMAP per gene; all D0-D3 cells shown in one shared embedding",
        x="single-cell state embedding 1",
        y="single-cell state embedding 2"
    )
    + theme_bw(base_size=12)
    + theme(
        figure_size=(12.0, 8.0),
        panel_grid_major=element_blank(),
        panel_grid_minor=element_blank(),
        panel_border=element_rect(color="#7f7f7f", fill=None, size=0.8),
        axis_title=element_text(size=13),
        axis_text=element_text(size=10),
        plot_title=element_text(size=17, weight="bold"),
        plot_subtitle=element_text(size=10.5),
        strip_background=element_rect(fill="#D9D9D9", color="#7f7f7f", size=0.8),
        strip_text=element_text(size=14, weight="bold"),
        legend_title=element_text(size=11, weight="bold"),
        legend_text=element_text(size=10)
    )
)

p_facet.save(OUTDIR / "carrier_seeded_allHippo_allDays_singleUMAP_simpleColor_faceted_plotnine.png",
             width=12.0, height=8.0, dpi=450, verbose=False)
p_facet.save(OUTDIR / "carrier_seeded_allHippo_allDays_singleUMAP_simpleColor_faceted_plotnine.pdf",
             width=12.0, height=8.0, verbose=False)
p_facet.save(OUTDIR / "carrier_seeded_allHippo_allDays_singleUMAP_simpleColor_faceted_plotnine.svg",
             width=12.0, height=8.0, verbose=False)

p_facet
```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f3ce8dbe610>
```

## Markdown cell 39

#### UMAP showing cells carrying targeting guides

## Code cell 40

```python
adata = adata_int.copy()
```

## Code cell 41

```python
# ============================================================
# Hippo guide-carrier + non-targeting control UMAP location plots
# One giant notebook cell / .py-compatible version
# ============================================================
#
# What this does:
#   - Finds cells carrying guides targeting:
#       ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]
#   - Adds one extra UMAP for non-targeting control cells
#   - Grey background = all cells
#   - Blue points = carrier/control cells
#   - All days pooled
#   - No outlines around dots
#
# This is a location plot only.
# It does NOT calculate neighborhood enrichment or Milo DA.
#
# ============================================================


# ============================================================
# Imports
# ============================================================

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse


# ============================================================
# User settings
# ============================================================

# In notebook, keep this as None because adata is already loaded.
# For .py use, set:
# H5AD_PATH = "your_file.h5ad"
H5AD_PATH = None

OUTDIR = Path("hippo_and_nt_carrier_location_umaps")
OUTDIR.mkdir(parents=True, exist_ok=True)

HIPPO_TARGETS = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]
CONTROL_TARGET = "non-targeting"

# Plot these targets. Includes non-targeting control panel.
PLOT_TARGETS = HIPPO_TARGETS + [CONTROL_TARGET]

GUIDE_ASSIGNMENT_KEY = "guide_assignment"
GUIDE_TARGET_CSV = "combined_TargetKD.csv"  # set to None if unavailable

DAY_COL = "day"
CELLTYPE_COL = "cell_type"

UMAP_KEY = "X_umap"

# Mode options:
#   "all_carriers"
#   "low_guide_load_leq2"
#   "single_selected_target"
MODE = "all_carriers"

# For the non-targeting control panel:
#   True  = plot NT cells with zero targeting genes
#   False = plot all cells carrying an NT guide, even if they also carry targeting guides
PURE_NT_ONLY = True

# Drawing style
BACKGROUND_POINT_SIZE = 1.0
BACKGROUND_ALPHA = 0.10

CARRIER_POINT_SIZE = 7 
CARRIER_ALPHA = 0.9
CARRIER_COLOR = "#2166AC"  # blue

SAVE_PDF = True
DPI = 300


# ============================================================
# Optional adata loading for .py use
# ============================================================

if H5AD_PATH is not None:
    import scanpy as sc
    adata = sc.read_h5ad(H5AD_PATH)

if "adata" not in globals():
    raise RuntimeError(
        "No `adata` found. Load adata first in your notebook, or set H5AD_PATH."
    )


# ============================================================
# Helper functions
# ============================================================

def looks_nt(x):
    s = str(x).lower()

    patterns = [
        "non-target",
        "nontarget",
        "non_target",
        "non targeting",
        "negative",
        "control",
        "scramble",
        "scrambled",
        "safe",
        "ntc",
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

    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]

    for col in df.columns:
        cl = str(col).lower()
        for cand in candidates:
            if cand.lower() in cl:
                return col

    return None


def get_guide_names(adata, guide_assignment_key=GUIDE_ASSIGNMENT_KEY):
    ga = adata.obsm[guide_assignment_key]

    if isinstance(ga, pd.DataFrame):
        return np.asarray(ga.columns.astype(str))

    if "guide_names" in adata.uns:
        return np.asarray(adata.uns["guide_names"]).astype(str)

    raise KeyError(
        "Could not find guide names. Expected adata.uns['guide_names'] "
        "or DataFrame columns in adata.obsm['guide_assignment']."
    )


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

    ga = np.asarray(ga)

    return pd.DataFrame(
        ga,
        index=adata.obs_names,
        columns=guide_names,
    ).astype(bool)


def build_guide_target_table(
    adata,
    guide_target_csv=GUIDE_TARGET_CSV,
    guide_assignment_key=GUIDE_ASSIGNMENT_KEY,
):
    guide_names = pd.Index(
        get_guide_names(adata, guide_assignment_key),
        name="guide",
    ).astype(str)

    base = pd.DataFrame({"guide": guide_names})

    if guide_target_csv is not None and os.path.exists(guide_target_csv):
        gt = pd.read_csv(guide_target_csv)

        guide_col = infer_col(
            gt,
            [
                "guide",
                "guide_name",
                "sgRNA",
                "sgrna",
                "gRNA",
                "oligo",
                "spacer",
            ],
        )

        target_col = infer_col(
            gt,
            [
                "target",
                "target_name",
                "target_gene",
                "gene",
                "symbol",
                "Target",
                "TargetGene",
            ],
        )

        if guide_col is None or target_col is None:
            raise ValueError(
                f"Could not infer guide/target columns in {guide_target_csv}. "
                f"Columns: {list(gt.columns)}"
            )

        gt = gt[[guide_col, target_col]].copy()
        gt.columns = ["guide", "target"]
        gt["guide"] = gt["guide"].astype(str)
        gt["target"] = gt["target"].astype(str)

        out = base.merge(gt, on="guide", how="left")
        out["target"] = out["target"].fillna(
            out["guide"].map(parse_target_from_guide)
        )

    else:
        out = base.copy()
        out["target"] = out["guide"].map(parse_target_from_guide)

    out["is_non_targeting"] = out["guide"].map(looks_nt) | out["target"].map(looks_nt)
    out.loc[out["is_non_targeting"], "target"] = CONTROL_TARGET

    return out


def get_target_assignment_df(guide_assignment, guide_target):
    """
    Build cell x target boolean matrix from cell x guide matrix.
    Efficient version that avoids pandas fragmentation warnings.
    """
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

    out = pd.DataFrame(series_by_target, index=guide_assignment.index)

    return out.astype(bool).copy()


def get_plot_mask_for_target(target, mode=MODE):
    """
    Return boolean mask of cells to highlight for a target.
    Special handling for non-targeting control.
    """
    if target not in target_assignment.columns:
        raise KeyError(f"{target} not found in target_assignment columns.")

    n_total_guides = adata.obs["n_total_guides"].values
    n_targeting_targets = adata.obs["n_targeting_targets"].values

    if target == CONTROL_TARGET:
        nt_mask = target_assignment[CONTROL_TARGET].values

        if PURE_NT_ONLY:
            base_mask = nt_mask & (n_targeting_targets == 0)
        else:
            base_mask = nt_mask

        if mode == "all_carriers":
            mask = base_mask

        elif mode == "low_guide_load_leq2":
            mask = base_mask & (n_total_guides <= 2)

        elif mode == "single_selected_target":
            # For NT, this means pure NT with zero targeting genes.
            mask = nt_mask & (n_targeting_targets == 0)

        else:
            raise ValueError(f"Unknown MODE: {mode}")

        return mask

    # Target-gene panels
    target_mask = target_assignment[target].values

    if mode == "all_carriers":
        mask = target_mask

    elif mode == "low_guide_load_leq2":
        mask = target_mask & (n_total_guides <= 2)

    elif mode == "single_selected_target":
        mask = target_mask & (n_targeting_targets == 1)

    else:
        raise ValueError(f"Unknown MODE: {mode}")

    return mask


def plot_carrier_location_umap(
    adata,
    target,
    mode=MODE,
    carrier_color=CARRIER_COLOR,
    outdir=OUTDIR,
):
    """
    One UMAP per target.
    Grey = all cells.
    Blue = cells carrying guides targeting this gene, or pure non-targeting cells.
    """
    mask = get_plot_mask_for_target(target, mode=mode)

    bg = adata.obs[["UMAP1", "UMAP2"]].copy()

    carriers = adata.obs.loc[
        mask,
        ["UMAP1", "UMAP2", "day_clean", "cell_type_clean"],
    ].copy()

    if len(carriers) == 0:
        print(f"No cells found for {target} in mode={mode}")
        return None

    fig, ax = plt.subplots(figsize=(6.8, 6.2))

    # Background all cells
    ax.scatter(
        bg["UMAP1"],
        bg["UMAP2"],
        s=BACKGROUND_POINT_SIZE,
        c="lightgrey",
        alpha=BACKGROUND_ALPHA,
        linewidths=0,
        rasterized=True,
    )

    # Highlighted cells
    ax.scatter(
        carriers["UMAP1"],
        carriers["UMAP2"],
        s=CARRIER_POINT_SIZE,
        c=carrier_color,
        alpha=CARRIER_ALPHA,
        linewidths=0,
        edgecolors="none",
        rasterized=True,
    )

    n_carriers = len(carriers)
    frac = n_carriers / adata.n_obs

    day_counts = carriers["day_clean"].value_counts().sort_index()
    day_text = ", ".join([f"{d}:{n}" for d, n in day_counts.items()])

    if target == CONTROL_TARGET:
        title_target = "Non-targeting control"
        purity_note = "pure NT only" if PURE_NT_ONLY else "all NT-guide carriers"
    else:
        title_target = target
        purity_note = "guide-carrier cells"

    ax.set_title(
        f"{title_target}\n"
        f"{mode}; {purity_note}; all days pooled",
        fontsize=12,
        fontweight="bold",
    )

    ax.text(
        0.02,
        0.02,
        f"n cells = {n_carriers}\n"
        f"fraction of all cells = {frac:.3f}\n"
        f"day counts: {day_text}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="none",
            alpha=0.75,
        ),
    )

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")

    fig.tight_layout()

    safe_target = str(target).replace("/", "_").replace(" ", "_")

    out_png = outdir / f"Carrier_location_UMAP_{safe_target}_{mode}.png"
    fig.savefig(out_png, dpi=DPI, bbox_inches="tight")

    if SAVE_PDF:
        out_pdf = outdir / f"Carrier_location_UMAP_{safe_target}_{mode}.pdf"
        fig.savefig(out_pdf, bbox_inches="tight")

    plt.show()

    print("Saved:", out_png)

    return fig


def plot_carrier_location_grid(
    adata,
    targets=PLOT_TARGETS,
    mode=MODE,
    carrier_color=CARRIER_COLOR,
    outdir=OUTDIR,
):
    """
    One grid figure with Hippo genes plus non-targeting control.
    """
    targets_present = [t for t in targets if t in target_assignment.columns]

    if len(targets_present) == 0:
        print("No requested targets found.")
        return None

    n = len(targets_present)

    fig, axes = plt.subplots(
        1,
        n,
        figsize=(4.6 * n, 4.8),
        squeeze=False,
    )

    bg = adata.obs[["UMAP1", "UMAP2"]].copy()

    for i, target in enumerate(targets_present):
        ax = axes[0, i]

        mask = get_plot_mask_for_target(target, mode=mode)
        carriers = adata.obs.loc[mask, ["UMAP1", "UMAP2"]].copy()

        ax.scatter(
            bg["UMAP1"],
            bg["UMAP2"],
            s=BACKGROUND_POINT_SIZE,
            c="lightgrey",
            alpha=BACKGROUND_ALPHA,
            linewidths=0,
            rasterized=True,
        )

        ax.scatter(
            carriers["UMAP1"],
            carriers["UMAP2"],
            s=CARRIER_POINT_SIZE,
            c=carrier_color,
            alpha=CARRIER_ALPHA,
            linewidths=0,
            edgecolors="none",
            rasterized=True,
        )

        if target == CONTROL_TARGET:
            label = "NT control"
        else:
            label = target

        ax.set_title(
            f"{label}\n"
            f"n={len(carriers)}",
            fontsize=10,
            fontweight="bold",
        )

        ax.set_xlabel("UMAP1")

        if i == 0:
            ax.set_ylabel("UMAP2")
        else:
            ax.set_ylabel("")

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal", adjustable="box")

    nt_note = "pure NT only" if PURE_NT_ONLY else "all NT-guide carriers"

    fig.suptitle(
        f"Hippo guide-carrier and non-targeting control locations\n"
        f"{mode}; NT panel = {nt_note}; all days pooled",
        fontsize=14,
        fontweight="bold",
    )

    fig.tight_layout()

    out_png = outdir / f"Hippo_plus_NT_carrier_location_grid_{mode}.png"
    fig.savefig(out_png, dpi=DPI, bbox_inches="tight")

    if SAVE_PDF:
        out_pdf = outdir / f"Hippo_plus_NT_carrier_location_grid_{mode}.pdf"
        fig.savefig(out_pdf, bbox_inches="tight")

    plt.show()

    print("Saved:", out_png)

    return fig


# ============================================================
# Prepare adata metadata
# ============================================================

if not adata.obs_names.is_unique:
    adata.obs_names_make_unique()

adata.obs["day_clean"] = adata.obs[DAY_COL].astype(str).str.lower()
adata.obs["cell_type_clean"] = adata.obs[CELLTYPE_COL].astype(str)

if UMAP_KEY not in adata.obsm:
    raise KeyError(f"adata.obsm['{UMAP_KEY}'] not found.")

umap = np.asarray(adata.obsm[UMAP_KEY])
adata.obs["UMAP1"] = umap[:, 0]
adata.obs["UMAP2"] = umap[:, 1]


# ============================================================
# Build guide / target assignments
# ============================================================

guide_assignment = get_guide_assignment_df(adata)
guide_target = build_guide_target_table(adata, guide_target_csv=GUIDE_TARGET_CSV)
target_assignment = get_target_assignment_df(guide_assignment, guide_target)

guide_target.to_csv(OUTDIR / "guide_to_target_map.csv", index=False)

if CONTROL_TARGET not in target_assignment.columns:
    raise KeyError(
        f"Could not find {CONTROL_TARGET!r}. "
        "Check combined_TargetKD.csv or guide naming."
    )

targeting_cols = [c for c in target_assignment.columns if c != CONTROL_TARGET]

adata.obs["n_total_guides"] = guide_assignment.sum(axis=1).astype(int).values
adata.obs["n_targeting_targets"] = target_assignment[targeting_cols].sum(axis=1).astype(int).values


# ============================================================
# Save carrier count tables
# ============================================================

carrier_counts = pd.DataFrame({
    "target": target_assignment.columns,
    "n_cells_all_carriers": target_assignment.sum(axis=0).values,
})

carrier_counts.to_csv(OUTDIR / "target_carrier_counts.csv", index=False)

plot_counts = []

for target in PLOT_TARGETS:
    if target in target_assignment.columns:
        mask = get_plot_mask_for_target(target, mode=MODE)
        plot_counts.append({
            "target": target,
            "mode": MODE,
            "n_cells_plotted": int(mask.sum()),
            "pure_nt_only": PURE_NT_ONLY if target == CONTROL_TARGET else np.nan,
        })

plot_counts = pd.DataFrame(plot_counts)
plot_counts.to_csv(OUTDIR / "plotted_target_counts.csv", index=False)

print("Targets present:")
print([t for t in PLOT_TARGETS if t in target_assignment.columns])

print("\nCells plotted:")
print(plot_counts)


# ============================================================
# Generate plots
# ============================================================

for target in PLOT_TARGETS:
    if target in target_assignment.columns:
        plot_carrier_location_umap(
            adata=adata,
            target=target,
            mode=MODE,
            carrier_color=CARRIER_COLOR,
        )

plot_carrier_location_grid(
    adata=adata,
    targets=PLOT_TARGETS,
    mode=MODE,
    carrier_color=CARRIER_COLOR,
)


print("\nDone.")
print("Output directory:", OUTDIR.resolve())
```

### Text output (stream)

```text
Targets present:
['NF2', 'LATS2', 'SAV1', 'TAOK1', 'AMOTL2', 'non-targeting']

Cells plotted:
          target          mode  n_cells_plotted pure_nt_only
0            NF2  all_carriers             2354          NaN
1          LATS2  all_carriers             2695          NaN
2           SAV1  all_carriers             1595          NaN
3          TAOK1  all_carriers             2530          NaN
4         AMOTL2  all_carriers             2381          NaN
5  non-targeting  all_carriers              748         True

```

### Text output (text/plain)

```text
<Figure size 680x620 with 1 Axes>
```

### Text output (stream)

```text
Saved: hippo_and_nt_carrier_location_umaps/Carrier_location_UMAP_NF2_all_carriers.png

```

### Text output (text/plain)

```text
<Figure size 680x620 with 1 Axes>
```

### Text output (stream)

```text
Saved: hippo_and_nt_carrier_location_umaps/Carrier_location_UMAP_LATS2_all_carriers.png

```

### Text output (text/plain)

```text
<Figure size 680x620 with 1 Axes>
```

### Text output (stream)

```text
Saved: hippo_and_nt_carrier_location_umaps/Carrier_location_UMAP_SAV1_all_carriers.png

```

### Text output (text/plain)

```text
<Figure size 680x620 with 1 Axes>
```

### Text output (stream)

```text
Saved: hippo_and_nt_carrier_location_umaps/Carrier_location_UMAP_TAOK1_all_carriers.png

```

### Text output (text/plain)

```text
<Figure size 680x620 with 1 Axes>
```

### Text output (stream)

```text
Saved: hippo_and_nt_carrier_location_umaps/Carrier_location_UMAP_AMOTL2_all_carriers.png

```

### Text output (text/plain)

```text
<Figure size 680x620 with 1 Axes>
```

### Text output (stream)

```text
Saved: hippo_and_nt_carrier_location_umaps/Carrier_location_UMAP_non-targeting_all_carriers.png

```

### Text output (text/plain)

```text
<Figure size 2760x480 with 6 Axes>
```

### Text output (stream)

```text
Saved: hippo_and_nt_carrier_location_umaps/Hippo_plus_NT_carrier_location_grid_all_carriers.png

Done.
Output directory: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/hippo_and_nt_carrier_location_umaps

```

## Markdown cell 42

#### Panel B dot plot

## Code cell 43

```python
import os
import re
import zipfile
import warnings
import numpy as np
import pandas as pd

from scipy import sparse
from scipy.stats import fisher_exact, pearsonr
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ============================================================
# Input / output
# ============================================================

OUTDIR = "panelB_direct_branch_contrast_from_adata"
os.makedirs(OUTDIR, exist_ok=True)

GUIDE_METADATA_TSV = "300genes_guide_metadata_v43.tsv"

PSEUDOBULK_OUT = os.path.join(
    OUTDIR,
    "pseudobulk_counts_by_sample_target_day_celltype.csv"
)

PANELB_SOURCE_OUT = os.path.join(
    OUTDIR,
    "panelB_direct_branch_contrast_source.csv"
)

# ============================================================
# adata fields
# ============================================================

DAY_COL = "day"
CELLTYPE_COL = "cell_type"

GUIDE_NAMES_KEY = "guide_names"
GUIDE_ASSIGNMENT_KEY = "guide_assignment"

GUIDE_ID_COL = "guide_id"
TARGET_COL = "putative_target_genes"

# Use sample_id by default. Add "well" or "biological_sample" if you want finer pseudobulk groups.
SAMPLE_COLS = ["sample_id"]

# ============================================================
# Analysis mode
# ============================================================
# all_carriers:
#   all cells carrying target are counted for that target.
# low_guide_load_leq2:
#   only cells with <=2 total guides.
# single_selected_target:
#   target cells must carry exactly one targeting target.
MODE = "all_carriers"

# ============================================================
# Targets and controls
# ============================================================

CONTROL_TARGET = "non-targeting_control"

CONTROL_PATTERNS = [
    r"non[-_ ]?target",
    r"non[-_ ]?targeting",
    r"\bNTC\b",
    r"negative",
    r"neg[-_ ]?ctrl",
    r"control",
    r"scram",
    r"scrambled",
    r"safe",
    r"aavs1",
]

HIPPO_TARGETS = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]

HIGHLIGHT_TARGETS = [
    "NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2",
    "KDR", "ETV2", "BMPR1A", "SMAD4"
]

# ============================================================
# Branch-state definitions
# Update these if your cell_type labels differ.
# ============================================================

D2_P15_PATTERNS = [
    r"\bPre-EC\b",
    r"\bP15\b",
    r"aEC progenitor",
]

D2_P17_PATTERNS = [
    r"\bLat-Meso\b",
    r"\bP17\b",
    r"LM.*epithelium",
    r"LM-like",
]

D3_EC_PATTERNS = [
    r"\bartEC\b",
]

D3_ALT_PATTERNS = [
    r"Fibro-Mesen",
]

PSEUDOCOUNT = 0.5

# Minimum number of target-carrier cells in the two-state branch contrast.
MIN_D2_BRANCH_CELLS = 10
MIN_D3_ENDPOINT_CELLS = 10
```

## Code cell 44

```python
def looks_control(x):
    s = str(x).lower()
    return any(re.search(p, s, flags=re.IGNORECASE) for p in CONTROL_PATTERNS)


def clean_target_string(x):
    if pd.isna(x):
        return []

    s = str(x).strip()

    if s.lower() in ["", "nan", "none", "na"]:
        return []

    s = (
        s.replace("[", "")
         .replace("]", "")
         .replace("'", "")
         .replace('"', "")
    )

    parts = re.split(r"[;,|/]+", s)
    parts = [p.strip() for p in parts if p.strip() != ""]

    return parts


def canonical_guide_id(x):
    return str(x).strip().replace(" ", "")


def get_guide_assignment_matrix(adata):
    ga = adata.obsm[GUIDE_ASSIGNMENT_KEY]

    if sparse.issparse(ga):
        return ga.tocsr().astype(bool)

    return np.asarray(ga).astype(bool)


def any_guides_per_cell(guide_assignment, guide_indices):
    if len(guide_indices) == 0:
        return np.zeros(guide_assignment.shape[0], dtype=bool)

    if sparse.issparse(guide_assignment):
        return np.asarray(guide_assignment[:, guide_indices].sum(axis=1)).ravel() > 0

    return guide_assignment[:, guide_indices].any(axis=1)


def load_guide_target_from_metadata(adata):
    guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))
    guide_name_to_idx = {g: i for i, g in enumerate(guide_names)}

    meta = pd.read_csv(GUIDE_METADATA_TSV, sep="\t")

    if GUIDE_ID_COL not in meta.columns:
        raise ValueError(f"{GUIDE_ID_COL} not found. Columns: {meta.columns.tolist()}")

    if TARGET_COL not in meta.columns:
        raise ValueError(f"{TARGET_COL} not found. Columns: {meta.columns.tolist()}")

    meta = meta[[GUIDE_ID_COL, TARGET_COL]].copy()
    meta.columns = ["guide_id", "putative_target_genes"]
    meta["guide_id"] = meta["guide_id"].astype(str)

    # Exact matching.
    meta["guide_index"] = meta["guide_id"].map(guide_name_to_idx)
    meta["match_method"] = np.where(meta["guide_index"].notna(), "exact_guide_id_match", None)

    # Fallback canonical matching.
    if meta["guide_index"].isna().any():
        canonical_to_guide = {}
        duplicated = set()

        for g in guide_names:
            cg = canonical_guide_id(g)
            if cg in canonical_to_guide:
                duplicated.add(cg)
            else:
                canonical_to_guide[cg] = g

        for cg in duplicated:
            canonical_to_guide.pop(cg, None)

        for idx in meta.index[meta["guide_index"].isna()]:
            cg = canonical_guide_id(meta.at[idx, "guide_id"])
            if cg in canonical_to_guide:
                matched_guide = canonical_to_guide[cg]
                meta.at[idx, "guide_index"] = guide_name_to_idx[matched_guide]
                meta.at[idx, "match_method"] = "canonical_guide_id_match"

    mapped = meta[meta["guide_index"].notna()].copy()
    mapped["guide_index"] = mapped["guide_index"].astype(int)
    mapped["guide_name"] = [guide_names[i] for i in mapped["guide_index"]]

    records = []

    for _, row in mapped.iterrows():
        guide_index = int(row["guide_index"])
        guide_name = str(row["guide_name"])
        raw_target = row["putative_target_genes"]

        targets = clean_target_string(raw_target)

        if looks_control(guide_name) or looks_control(raw_target):
            targets = [CONTROL_TARGET]

        for target in targets:
            records.append({
                "guide_index": guide_index,
                "guide_name": guide_name,
                "target_name": str(target),
                "raw_target_entry": raw_target,
                "match_method": row["match_method"],
            })

    guide_target = pd.DataFrame(records).drop_duplicates()

    print("\nGuide metadata mapping QC")
    print(f"adata guide columns:            {len(guide_names):,}")
    print(f"metadata guide rows:            {len(meta):,}")
    print(f"metadata rows matched to adata: {mapped['guide_id'].nunique():,}")
    print(f"expanded guide-target rows:     {len(guide_target):,}")
    print(f"unique mapped targets:          {guide_target['target_name'].nunique():,}")

    print("\nExample guide-target rows:")
    print(guide_target.head(20).to_string(index=False))

    guide_target.to_csv(
        os.path.join(OUTDIR, "guide_to_target_map_from_300genes_metadata.csv"),
        index=False
    )

    return guide_target


def build_target_masks(adata, guide_target):
    guide_assignment = get_guide_assignment_matrix(adata)

    if guide_assignment.shape[1] != len(adata.uns[GUIDE_NAMES_KEY]):
        raise ValueError(
            f"guide_assignment has {guide_assignment.shape[1]} columns but "
            f"adata.uns['{GUIDE_NAMES_KEY}'] has {len(adata.uns[GUIDE_NAMES_KEY])} names."
        )

    target_masks = {}
    target_to_guide_indices = {}

    for target, sub in guide_target.groupby("target_name"):
        idx = sorted(set(sub["guide_index"].astype(int).tolist()))
        target_to_guide_indices[target] = idx
        target_masks[target] = any_guides_per_cell(guide_assignment, idx)

    control_like_targets = [
        t for t in target_masks
        if looks_control(t) or str(t).lower() == CONTROL_TARGET
    ]

    if len(control_like_targets) == 0:
        raise ValueError("No non-targeting/control target found.")

    control_mask = np.zeros(adata.n_obs, dtype=bool)
    for t in control_like_targets:
        control_mask |= target_masks[t]

    # targeting guide count per cell: exclude controls
    targeting_targets = [t for t in target_masks if t not in control_like_targets]
    n_targeting_targets = np.zeros(adata.n_obs, dtype=int)
    for t in targeting_targets:
        n_targeting_targets += target_masks[t].astype(int)

    # total guide count
    ga = get_guide_assignment_matrix(adata)
    if sparse.issparse(ga):
        n_total_guides = np.asarray(ga.sum(axis=1)).ravel().astype(int)
    else:
        n_total_guides = ga.sum(axis=1).astype(int)

    print("\nTarget mask QC")
    print(f"Targets with masks: {len(target_masks):,}")
    print(f"Control-like target labels: {control_like_targets}")
    print(f"Control carrier cells: {control_mask.sum():,}")

    print("\nKey target guide/cell counts:")
    for t in HIPPO_TARGETS + ["KDR", "ETV2", "BMPR1A", "SMAD4"]:
        matches = [x for x in target_masks if x.upper() == t.upper()]
        if matches:
            tt = matches[0]
            print(
                f"  {t:8s}: guides={len(target_to_guide_indices[tt]):3d}, "
                f"carrier_cells={target_masks[tt].sum():6d}"
            )
        else:
            print(f"  {t:8s}: NOT FOUND")

    target_counts = pd.DataFrame({
        "target": list(target_masks.keys()),
        "n_carrier_cells": [target_masks[t].sum() for t in target_masks],
        "n_guides": [len(target_to_guide_indices.get(t, [])) for t in target_masks],
    })

    target_counts.to_csv(
        os.path.join(OUTDIR, "target_carrier_counts.csv"),
        index=False
    )

    return target_masks, target_to_guide_indices, control_mask, n_total_guides, n_targeting_targets
```

## Code cell 45

```python
def clean_day(x):
    s = str(x).lower().strip()
    s = s.replace("day", "d").replace(" ", "")

    if s in ["0", "0.0"]:
        return "D0"
    if s in ["1", "1.0"]:
        return "D1"
    if s in ["2", "2.0"]:
        return "D2"
    if s in ["3", "3.0"]:
        return "D3"

    return s.upper()


def matches_any(label, patterns):
    s = str(label)
    return any(re.search(p, s, flags=re.IGNORECASE) for p in patterns)


def classify_branch_state(day, state_label):
    if day == "D2":
        if matches_any(state_label, D2_P15_PATTERNS):
            return "D2_P15"
        if matches_any(state_label, D2_P17_PATTERNS):
            return "D2_P17"

    if day == "D3":
        if matches_any(state_label, D3_EC_PATTERNS):
            return "D3_EC"
        if matches_any(state_label, D3_ALT_PATTERNS):
            return "D3_ALT"

    return None


def mode_filter_for_target(
    target_mask,
    control_mask,
    n_total_guides,
    n_targeting_targets,
    mode=MODE
):
    if mode == "all_carriers":
        target_use = target_mask.copy()
        control_use = control_mask & (~target_mask)

    elif mode == "low_guide_load_leq2":
        low = n_total_guides <= 2
        target_use = target_mask & low
        control_use = control_mask & (~target_mask) & low

    elif mode == "single_selected_target":
        target_use = target_mask & (n_targeting_targets == 1)
        control_use = control_mask & (~target_mask) & (n_targeting_targets == 0)

    else:
        raise ValueError(f"Unknown mode: {mode}")

    return target_use, control_use
```

## Code cell 46

```python
# ============================================================
# Build guide/target masks
# ============================================================

guide_target = load_guide_target_from_metadata(adata)

target_masks, target_to_guide_indices, control_mask, n_total_guides, n_targeting_targets = (
    build_target_masks(adata, guide_target)
)

# ============================================================
# Build base cell metadata table
# ============================================================

obs = adata.obs.copy()

obs["cell_id"] = adata.obs_names.astype(str)
obs["day_clean"] = obs[DAY_COL].map(clean_day)
obs["cell_type_clean"] = obs[CELLTYPE_COL].astype(str)
obs["branch_state"] = [
    classify_branch_state(day, ct)
    for day, ct in zip(obs["day_clean"], obs["cell_type_clean"])
]

# sample columns
sample_cols_present = [c for c in SAMPLE_COLS if c in obs.columns]

if len(sample_cols_present) == 0:
    obs["sample_id_auto"] = "sample_1"
    sample_cols_present = ["sample_id_auto"]

print("Using sample columns:", sample_cols_present)

# ============================================================
# Export pseudobulk count table
# ============================================================

records = []

# Include all perturbation targets and a collapsed control target.
targets_to_export = [t for t in target_masks.keys() if t != CONTROL_TARGET]
targets_to_export = sorted(targets_to_export)

for target in targets_to_export:
    target_mask = target_masks[target]

    # Use the selected mode.
    target_use, _ = mode_filter_for_target(
        target_mask=target_mask,
        control_mask=control_mask,
        n_total_guides=n_total_guides,
        n_targeting_targets=n_targeting_targets,
        mode=MODE
    )

    if target_use.sum() == 0:
        continue

    sub = obs.loc[target_use, sample_cols_present + ["day_clean", "cell_type_clean", "branch_state"]].copy()
    sub["target"] = target

    grouped = (
        sub
        .groupby(["target"] + sample_cols_present + ["day_clean", "cell_type_clean", "branch_state"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )

    records.append(grouped)

# Add collapsed non-targeting controls.
control_use = control_mask.copy()

if MODE == "low_guide_load_leq2":
    control_use = control_use & (n_total_guides <= 2)

elif MODE == "single_selected_target":
    control_use = control_use & (n_targeting_targets == 0)

sub = obs.loc[control_use, sample_cols_present + ["day_clean", "cell_type_clean", "branch_state"]].copy()
sub["target"] = CONTROL_TARGET

grouped = (
    sub
    .groupby(["target"] + sample_cols_present + ["day_clean", "cell_type_clean", "branch_state"], dropna=False)
    .size()
    .reset_index(name="n_cells")
)

records.append(grouped)

pseudobulk = pd.concat(records, ignore_index=True)

# Rename for easier downstream use.
rename_map = {
    "day_clean": "day",
    "cell_type_clean": "cell_type",
}
pseudobulk = pseudobulk.rename(columns=rename_map)

pseudobulk["mode"] = MODE

pseudobulk.to_csv(PSEUDOBULK_OUT, index=False)

print("Saved:", PSEUDOBULK_OUT)
print("Shape:", pseudobulk.shape)
display(pseudobulk.head())

print("\nBranch-state counts in pseudobulk:")
print(pseudobulk["branch_state"].value_counts(dropna=False))
```

### Text output (stream)

```text

Guide metadata mapping QC
adata guide columns:            1,922
metadata guide rows:            1,926
metadata rows matched to adata: 1,922
expanded guide-target rows:     1,922
unique mapped targets:          299

Example guide-target rows:
 guide_index                 guide_name target_name raw_target_entry         match_method
           0   ACAA1_+_38178488.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           1   ACAA1_+_38178517.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           2   ACAA1_+_38178559.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           3   ACAA1_+_38178570.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           4   ACAA1_+_38178620.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           5   ACAA1_-_38178575.23-P1P2       ACAA1            ACAA1 exact_guide_id_match
           6 ACTL7B_+_111618040.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
           7 ACTL7B_+_111618198.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
           8 ACTL7B_+_111618251.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
           9 ACTL7B_-_111617995.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
          10 ACTL7B_-_111618185.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
          11 ACTL7B_-_111618259.23-P1P2      ACTL7B           ACTL7B exact_guide_id_match
          12  ACVR1B_+_52345347.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          13  ACVR1B_+_52345561.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          14  ACVR1B_+_52345762.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          15  ACVR1B_-_52345530.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          16  ACVR1B_-_52345533.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          17  ACVR1B_-_52345539.23-P1P2      ACVR1B           ACVR1B exact_guide_id_match
          18 ACVR2A_+_148602340.23-P1P2      ACVR2A           ACVR2A exact_guide_id_match
          19 ACVR2A_+_148602343.23-P1P2      ACVR2A           ACVR2A exact_guide_id_match

```

### Text output (stream)

```text

Target mask QC
Targets with masks: 299
Control-like target labels: ['non-targeting_control']
Control carrier cells: 21,671

Key target guide/cell counts:
  NF2     : guides=  6, carrier_cells=  2354
  LATS2   : guides=  6, carrier_cells=  2695
  SAV1    : guides=  6, carrier_cells=  1595
  TAOK1   : guides=  6, carrier_cells=  2530
  AMOTL2  : guides=  6, carrier_cells=  2381
  KDR     : guides=  6, carrier_cells=  1798
  ETV2    : guides=  6, carrier_cells=  2275
  BMPR1A  : guides=  6, carrier_cells=  1643
  SMAD4   : guides=  6, carrier_cells=  3325
Using sample columns: ['sample_id']
Saved: panelB_direct_branch_contrast_from_adata/pseudobulk_counts_by_sample_target_day_celltype.csv
Shape: (149640, 7)

```

### Text output (text/plain)

```text
  target               sample_id day    cell_type branch_state  n_cells  \
0  ACAA1  300genes:eTN5_GEX_sub1  D0  Fibro-Mesen       D2_P15        0   
1  ACAA1  300genes:eTN5_GEX_sub1  D0  Fibro-Mesen       D2_P17        0   
2  ACAA1  300genes:eTN5_GEX_sub1  D0  Fibro-Mesen       D3_ALT        0   
3  ACAA1  300genes:eTN5_GEX_sub1  D0  Fibro-Mesen        D3_EC        0   
4  ACAA1  300genes:eTN5_GEX_sub1  D0  Fibro-Mesen          NaN        0   

           mode  
0  all_carriers  
1  all_carriers  
2  all_carriers  
3  all_carriers  
4  all_carriers  
```

### Text output (stream)

```text

Branch-state counts in pseudobulk:
branch_state
D2_P15    29928
D2_P17    29928
D3_ALT    29928
D3_EC     29928
NaN       29928
Name: count, dtype: int64

```

## Code cell 47

```python
def log2_ratio_of_ratios(a, b, c, d, pseudocount=PSEUDOCOUNT):
    return np.log2(
        ((a + pseudocount) / (b + pseudocount)) /
        ((c + pseudocount) / (d + pseudocount))
    )


def fisher_branch_pval(a, b, c, d):
    try:
        return fisher_exact([[a, b], [c, d]], alternative="two-sided")[1]
    except Exception:
        return np.nan


def count_state_for_mask(mask, branch_state):
    return int(np.sum(mask & (obs["branch_state"].values == branch_state)))


records = []

targets_for_panel = sorted([
    t for t in target_masks.keys()
    if t != CONTROL_TARGET and not looks_control(t)
])

for target in targets_for_panel:
    raw_target_mask = target_masks[target]

    target_use, control_use = mode_filter_for_target(
        target_mask=raw_target_mask,
        control_mask=control_mask,
        n_total_guides=n_total_guides,
        n_targeting_targets=n_targeting_targets,
        mode=MODE
    )

    # D2 branch counts.
    t_d2_p17 = count_state_for_mask(target_use, "D2_P17")
    t_d2_p15 = count_state_for_mask(target_use, "D2_P15")
    c_d2_p17 = count_state_for_mask(control_use, "D2_P17")
    c_d2_p15 = count_state_for_mask(control_use, "D2_P15")

    # D3 endpoint counts.
    t_d3_alt = count_state_for_mask(target_use, "D3_ALT")
    t_d3_ec = count_state_for_mask(target_use, "D3_EC")
    c_d3_alt = count_state_for_mask(control_use, "D3_ALT")
    c_d3_ec = count_state_for_mask(control_use, "D3_EC")

    d2_total = t_d2_p17 + t_d2_p15
    d3_total = t_d3_alt + t_d3_ec

    # Require enough target-carrier cells in both stages.
    if d2_total < MIN_D2_BRANCH_CELLS or d3_total < MIN_D3_ENDPOINT_CELLS:
        continue

    d2_branch_bias = log2_ratio_of_ratios(
        t_d2_p17, t_d2_p15,
        c_d2_p17, c_d2_p15
    )

    d3_endpoint_bias = log2_ratio_of_ratios(
        t_d3_alt, t_d3_ec,
        c_d3_alt, c_d3_ec
    )

    d2_pval = fisher_branch_pval(
        t_d2_p17, t_d2_p15,
        c_d2_p17, c_d2_p15
    )

    d3_pval = fisher_branch_pval(
        t_d3_alt, t_d3_ec,
        c_d3_alt, c_d3_ec
    )

    target_d2_p17_fraction = t_d2_p17 / d2_total if d2_total > 0 else np.nan
    control_d2_p17_fraction = c_d2_p17 / (c_d2_p17 + c_d2_p15) if (c_d2_p17 + c_d2_p15) > 0 else np.nan

    target_d3_alt_fraction = t_d3_alt / d3_total if d3_total > 0 else np.nan
    control_d3_alt_fraction = c_d3_alt / (c_d3_alt + c_d3_ec) if (c_d3_alt + c_d3_ec) > 0 else np.nan

    d2_fraction_shift = target_d2_p17_fraction - control_d2_p17_fraction
    d3_fraction_shift = target_d3_alt_fraction - control_d3_alt_fraction

    records.append({
        "target": target,
        "mode": MODE,

        "d2_branch_bias_log2OR": d2_branch_bias,
        "d3_endpoint_bias_log2OR": d3_endpoint_bias,

        "d2_pval": d2_pval,
        "d3_pval": d3_pval,

        "target_D2_P17_count": t_d2_p17,
        "target_D2_P15_count": t_d2_p15,
        "target_D3_ALT_count": t_d3_alt,
        "target_D3_EC_count": t_d3_ec,

        "control_D2_P17_count": c_d2_p17,
        "control_D2_P15_count": c_d2_p15,
        "control_D3_ALT_count": c_d3_alt,
        "control_D3_EC_count": c_d3_ec,

        "target_D2_P17_fraction": target_d2_p17_fraction,
        "control_D2_P17_fraction": control_d2_p17_fraction,
        "D2_P17_fraction_shift": d2_fraction_shift,

        "target_D3_ALT_fraction": target_d3_alt_fraction,
        "control_D3_ALT_fraction": control_d3_alt_fraction,
        "D3_ALT_fraction_shift": d3_fraction_shift,

        "D3_cell_proportion_shift_abs": abs(d3_fraction_shift),

        "n_D2_branch_cells": d2_total,
        "n_D3_endpoint_cells": d3_total,

        "hippo_core": target in HIPPO_TARGETS,
        "label": target if target in HIGHLIGHT_TARGETS else "",
    })

panel = pd.DataFrame(records)

panel["d2_fdr"] = multipletests(panel["d2_pval"].fillna(1), method="fdr_bh")[1]
panel["d3_fdr"] = multipletests(panel["d3_pval"].fillna(1), method="fdr_bh")[1]

panel["minus_log10_D3_FDR"] = -np.log10(panel["d3_fdr"].clip(lower=1e-300))
panel["min_D2D3_fdr"] = panel[["d2_fdr", "d3_fdr"]].min(axis=1)
panel["minus_log10_min_D2D3_FDR"] = -np.log10(panel["min_D2D3_fdr"].clip(lower=1e-300))

panel = panel.sort_values("minus_log10_D3_FDR", ascending=False).reset_index(drop=True)

panel.to_csv(PANELB_SOURCE_OUT, index=False)

print("Saved:", PANELB_SOURCE_OUT)
print("Shape:", panel.shape)

r = panel[["d2_branch_bias_log2OR", "d3_endpoint_bias_log2OR"]].corr().iloc[0, 1]
print(f"Pearson r = {r:.3f}")

display(panel.head(20))
```

### Text output (stream)

```text
Saved: panelB_direct_branch_contrast_from_adata/panelB_direct_branch_contrast_source.csv
Shape: (298, 30)
Pearson r = 0.678

```

### Text output (text/plain)

```text
    target          mode  d2_branch_bias_log2OR  d3_endpoint_bias_log2OR  \
0    PLCG1  all_carriers               1.823382                 2.852856   
1      NF2  all_carriers               0.717422                 2.271746   
2    LATS2  all_carriers               0.728864                 1.733282   
3      KDR  all_carriers               1.702044                 2.372270   
4     ETV2  all_carriers               0.689091                 2.304550   
5    SMAD4  all_carriers               0.787076                 1.447891   
6    CASP3  all_carriers               0.021217                 1.463692   
7      BAX  all_carriers              -0.058212                 1.324426   
8    MED12  all_carriers               2.100737                 2.280389   
9   AMOTL2  all_carriers               1.077947                 1.490889   
10  SETDB1  all_carriers               1.838170                 2.001938   
11    DBR1  all_carriers               0.800922                 1.450216   
12   PRKCD  all_carriers               0.511057                 1.459106   
13   APAF1  all_carriers              -0.180209                 1.322786   
14   ARIH1  all_carriers               0.426555                 1.749345   
15     VHL  all_carriers               0.502526                 1.434118   
16   HAND1  all_carriers               1.241245                 1.629149   
17   TAOK1  all_carriers               0.345337                 1.191100   
18   RESF1  all_carriers               1.069415                 1.637484   
19   RIBC1  all_carriers               1.103263                 1.372691   

         d2_pval       d3_pval  target_D2_P17_count  target_D2_P15_count  \
0   7.593731e-20  3.605327e-85                  563                   43   
1   1.799903e-04  4.283608e-84                  405                   66   
2   7.206108e-05  3.136838e-60                  437                   71   
3   2.813450e-11  6.533327e-57                  317                   26   
4   1.349379e-04  2.976994e-56                  442                   74   
5   4.116825e-06  5.537012e-54                  529                   82   
6   9.513216e-01  1.358906e-52                  352                   93   
7   7.308542e-01  4.668383e-46                  402                  111   
8   2.375104e-16  2.960366e-39                  370                   23   
9   3.238818e-08  1.090644e-38                  442                   56   
10  1.584422e-12  1.816551e-37                  323                   24   
11  4.018251e-05  1.462574e-34                  395                   61   
12  4.482497e-03  1.344394e-33                  431                   81   
13  3.517833e-01  3.271222e-32                  284                   85   
14  4.124091e-02  1.142835e-28                  297                   59   
15  1.151301e-02  4.314265e-28                  335                   63   
16  1.938005e-08  2.819397e-27                  365                   41   
17  5.266527e-02  4.463139e-27                  406                   86   
18  5.788686e-07  1.312700e-26                  359                   46   
19  1.380586e-08  5.334077e-26                  446                   56   

    target_D3_ALT_count  target_D3_EC_count  ...  \
0                   427                  94  ...   
1                   573                 190  ...   
2                   596                 292  ...   
3                   350                 106  ...   
4                   359                 115  ...   
5                   731                 427  ...   
6                   693                 402  ...   
7                   709                 462  ...   
8                   247                  79  ...   
9                   465                 263  ...   
10                  280                 110  ...   
11                  430                 248  ...   
12                  409                 237  ...   
13                  465                 296  ...   
14                  256                 119  ...   
15                  344                 201  ...   
16                  269                 137  ...   
17                  459                 318  ...   
18                  261                 132  ...   
19                  340                 206  ...   

    D3_cell_proportion_shift_abs  n_D2_branch_cells  n_D3_endpoint_cells  \
0                       0.434528                606                  521   
1                       0.366988                471                  763   
2                       0.291015                508                  888   
3                       0.378941                343                  456   
4                       0.370878                516                  474   
5                       0.245824                611                 1158   
6                       0.248397                445                 1095   
7                       0.225617                513                 1171   
8                       0.367133                393                  326   
9                       0.252799                498                  728   
10                      0.330050                347                  390   
11                      0.246232                456                  678   
12                      0.247725                512                  646   
13                      0.225450                369                  761   
14                      0.293024                356                  375   
15                      0.243691                398                  545   
16                      0.274714                406                  406   
17                      0.203550                492                  777   
18                      0.276006                405                  393   
19                      0.233691                502                  546   

    hippo_core   label        d2_fdr        d3_fdr  minus_log10_D3_FDR  \
0        False          2.262932e-17  1.074388e-82           81.968839   
1         True     NF2  1.449652e-03  6.382576e-82           81.195004   
2         True   LATS2  6.927162e-04  3.115925e-58           57.506413   
3        False     KDR  1.676816e-09  4.867328e-55           54.312709   
4        False    ETV2  1.182691e-03  1.774288e-54           53.750976   
5        False   SMAD4  5.841971e-05  2.750049e-52           51.560659   
6        False          9.843536e-01  5.785059e-51           50.237692   
7        False          8.746770e-01  1.738973e-44           43.759707   
8        False          3.538905e-14  9.802102e-38           37.008681   
9         True  AMOTL2  8.043064e-07  3.250120e-37           36.488101   
10       False          1.573859e-10  4.921201e-36           35.307929   
11       False          4.567371e-04  3.632059e-33           32.439847   
12       False          1.935919e-02  3.081764e-32           31.511201   
13       False          5.403681e-01  6.963030e-31           30.157202   
14       False          1.088362e-01  2.270433e-27           26.643891   
15       False          4.084377e-02  8.035318e-27           26.094997   
16       False          5.775256e-07  4.942238e-26           25.306076   
17        True   TAOK1  1.307854e-01  7.388974e-26           25.131416   
18       False          1.078143e-05  2.058866e-25           24.686372   
19       False          4.571275e-07  7.947775e-25           24.099754   

    min_D2D3_fdr  minus_log10_min_D2D3_FDR  
0   1.074388e-82                 81.968839  
1   6.382576e-82                 81.195004  
2   3.115925e-58                 57.506413  
3   4.867328e-55                 54.312709  
4   1.774288e-54                 53.750976  
5   2.750049e-52                 51.560659  
6   5.785059e-51                 50.237692  
7   1.738973e-44                 43.759707  
8   9.802102e-38                 37.008681  
9   3.250120e-37                 36.488101  
10  4.921201e-36                 35.307929  
11  3.632059e-33                 32.439847  
12  3.081764e-32                 31.511201  
13  6.963030e-31                 30.157202  
14  2.270433e-27                 26.643891  
15  8.035318e-27                 26.094997  
16  4.942238e-26                 25.306076  
17  7.388974e-26                 25.131416  
18  2.058866e-25                 24.686372  
19  7.947775e-25                 24.099754  

[20 rows x 30 columns]
```

## Code cell 48

```python
# ============================================================
# Cell 6 — plot Panel B in plotnine
# Constant dot size, color = -log10(D3 FDR), repel-like labels
# ============================================================

import os
import numpy as np
import pandas as pd

from scipy.stats import pearsonr

from plotnine import (
    ggplot, aes,
    geom_point, geom_text, geom_abline, geom_vline, geom_hline,
    scale_color_cmap,
    labs, theme_bw, theme,
    element_text, element_blank, element_rect,
    guides, guide_colorbar,
    coord_cartesian
)

# Needed for plotnine's adjust_text repel-like behavior.
# Install in notebook if needed:
# !pip install adjustText
try:
    import adjustText
    USE_REPEL = True
except Exception:
    USE_REPEL = False
    print("adjustText is not installed. Labels will use simple nudging instead.")
    print("Install with: pip install adjustText")

# ----------------------------
# Input / output
# ----------------------------

SOURCE_CSV = PANELB_SOURCE_OUT

OUT_PNG = os.path.join(OUTDIR, "panelB_direct_branch_contrast_constant_size_plotnine.png")
OUT_PDF = os.path.join(OUTDIR, "panelB_direct_branch_contrast_constant_size_plotnine.pdf")
OUT_SVG = os.path.join(OUTDIR, "panelB_direct_branch_contrast_constant_size_plotnine.svg")

CMAP_NAME = "plasma"

# Edit this if labels are too crowded.
LABEL_TARGETS = [
    "NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2",
    "KDR", "ETV2", "BMPR1A", "SMAD4"
]

# ----------------------------
# Load source table
# ----------------------------

df = pd.read_csv(SOURCE_CSV).copy()

df["x"] = pd.to_numeric(df["d2_branch_bias_log2OR"], errors="coerce")
df["y"] = pd.to_numeric(df["d3_endpoint_bias_log2OR"], errors="coerce")
df["color_var"] = pd.to_numeric(df["minus_log10_D3_FDR"], errors="coerce")

df = df.dropna(subset=["x", "y"]).copy()
df["color_var"] = df["color_var"].fillna(0)

# Cap color scale so one extreme target does not dominate the legend.
color_cap = max(5, np.nanpercentile(df["color_var"], 99))
df["color_var_plot"] = df["color_var"].clip(0, color_cap)

# ----------------------------
# Fit line and correlation
# ----------------------------

slope, intercept = np.polyfit(df["x"], df["y"], 1)
r, pval_corr = pearsonr(df["x"], df["y"])
n = df.shape[0]

# ----------------------------
# Label data
# ----------------------------

label_df = df[df["target"].isin(LABEL_TARGETS)].copy()
label_df["label"] = label_df["target"]

# ----------------------------
# Axis limits
# ----------------------------

xpad = 0.35
ypad = 0.35

xmin = df["x"].min() - xpad
xmax = df["x"].max() + xpad
ymin = df["y"].min() - ypad
ymax = df["y"].max() + ypad

# ----------------------------
# Repel-like label layer
# ----------------------------

if USE_REPEL:
    label_layer = geom_text(
        data=label_df,
        mapping=aes(label="label"),
        size=8.0,
        color="black",
        ha="center",
        va="center",
        adjust_text={
            "expand": (1.30, 1.55),
            "force_text": (0.80, 1.00),
            "force_static": (0.30, 0.45),
            "force_pull": (0.03, 0.03),
            "max_move": (0.35, 0.35),
            "min_arrow_len": 0.02,
            "arrowprops": {
                "arrowstyle": "-",
                "color": "#666666",
                "lw": 0.5,
                "alpha": 0.8,
            },
        }
    )

    # If your adjustText version errors on "force_static", use this older version instead:
    #
    # label_layer = geom_text(
    #     data=label_df,
    #     mapping=aes(label="label"),
    #     size=8.0,
    #     color="black",
    #     ha="center",
    #     va="center",
    #     adjust_text={
    #         "expand_text": (1.25, 1.45),
    #         "expand_points": (1.25, 1.45),
    #         "force_text": (0.80, 1.00),
    #         "force_points": (0.30, 0.45),
    #         "arrowprops": {
    #             "arrowstyle": "-",
    #             "color": "#666666",
    #             "lw": 0.5,
    #             "alpha": 0.8,
    #         },
    #     }
    # )

else:
    label_layer = geom_text(
        data=label_df,
        mapping=aes(label="label"),
        nudge_x=0.06,
        nudge_y=0.06,
        size=8.0,
        color="black",
        ha="left"
    )

# ----------------------------
# Plot
# ----------------------------

p = (
    ggplot(df, aes("x", "y"))
    + geom_vline(
        xintercept=0,
        linetype="dashed",
        color="#9E9E9E",
        size=0.5
    )
    + geom_hline(
        yintercept=0,
        linetype="dashed",
        color="#9E9E9E",
        size=0.5
    )
    + geom_abline(
        intercept=intercept,
        slope=slope,
        color="black",
        size=0.8
    )
    + geom_point(
        aes(color="color_var_plot"),
        size=3.3,
        alpha=0.90,
        stroke=0.25
    )
    + label_layer
    + scale_color_cmap(
        cmap_name=CMAP_NAME,
        name="-log10(D3 FDR)",
        limits=(0, color_cap)
    )
    + coord_cartesian(
        xlim=(xmin, xmax),
        ylim=(ymin, ymax)
    )
    + labs(
        title="Target-level association between D2 branch bias \nand D3 endpoint bias",
        subtitle=f"Direct branch contrast: n = {n} targets; Pearson r = {r:.2f}; P = {pval_corr:.1e}",
        x="D2 branch bias: \nlog2[(P17/P15)target / (P17/P15)control]",
        y="D3 endpoint bias: \nlog2[(Mesen/artEC)target / (Mesen/artEC)control]"
    )
    + guides(
        color=guide_colorbar(order=1)
    )
    + theme_bw()
    + theme(
        figure_size=(8.4, 6.4),
        panel_grid_minor=element_blank(),
        panel_grid_major=element_blank(),
        panel_background=element_rect(fill="white"),
        plot_background=element_rect(fill="white"),
        axis_title=element_text(size=13),
        axis_text=element_text(size=11),
        plot_title=element_text(size=15, weight="bold"),
        plot_subtitle=element_text(size=11),
        legend_title=element_text(size=11, weight="bold"),
        legend_text=element_text(size=10),
        legend_position="right"
    )
)

# ----------------------------
# Save
# ----------------------------

p.save(OUT_PNG, width=8.4, height=6.4, dpi=300)
p.save(OUT_PDF, width=8.4, height=6.4)
p.save(OUT_SVG, width=8.4, height=6.4)

print("Saved:")
print(OUT_PNG)
print(OUT_PDF)
print(OUT_SVG)
```

### Text output (stream)

```text
Saved:
panelB_direct_branch_contrast_from_adata/panelB_direct_branch_contrast_constant_size_plotnine.png
panelB_direct_branch_contrast_from_adata/panelB_direct_branch_contrast_constant_size_plotnine.pdf
panelB_direct_branch_contrast_from_adata/panelB_direct_branch_contrast_constant_size_plotnine.svg

```

## Markdown cell 49

#### Panel C guide level consistency

## Code cell 50

```python
import os
import re
import zipfile
import warnings
import numpy as np
import pandas as pd

from scipy import sparse
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from plotnine import (
    ggplot, aes,
    geom_point, geom_segment, geom_vline, geom_hline, geom_text,
    scale_fill_manual, scale_y_continuous,
    facet_wrap, coord_cartesian,
    labs, theme_bw, theme,
    element_text, element_blank, element_rect,
    guide_legend
)

warnings.filterwarnings("ignore")

# ============================================================
# Output
# ============================================================

OUTDIR = "panelC1_guide_level_consistency_forest_from_adata"
os.makedirs(OUTDIR, exist_ok=True)

GUIDE_EFFECTS_OUT = os.path.join(OUTDIR, "guide_level_D2_program_D3_cell_effects.csv")
TARGET_EFFECTS_OUT = os.path.join(OUTDIR, "target_level_D2_program_D3_cell_effects.csv")
PLOT_SOURCE_OUT = os.path.join(OUTDIR, "panelC1_forest_plot_source.csv")

# ============================================================
# adata fields
# ============================================================

DAY_COL = "day"
CELLTYPE_COL = "cell_type"

GUIDE_NAMES_KEY = "guide_names"
GUIDE_ASSIGNMENT_KEY = "guide_assignment"

CNMF_KEY = "cnmf_usage"
CNMF_ARRAY_PROGRAM_INDEX_BASE = 0

# ============================================================
# Guide metadata
# ============================================================

GUIDE_METADATA_TSV = "300genes_guide_metadata_v43.tsv"
GUIDE_ID_COL = "guide_id"
TARGET_COL = "putative_target_genes"

# ============================================================
# Targets for main Panel C1
# ============================================================

HIPPO_TARGETS = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]

TARGET_ORDER = ["NT control"] + HIPPO_TARGETS

TARGET_TO_GROUP = {
    "NT control": "NT control",
    "NF2": "Hippo",
    "LATS2": "Hippo",
    "SAV1": "Hippo",
    "TAOK1": "Hippo",
    "AMOTL2": "Hippo",
}

GROUP_COLORS = {
    "NT control": "#BDBDBD",
    "Hippo": "#D73027",
}

# ============================================================
# Controls
# ============================================================

CONTROL_TARGET = "non-targeting_control"

CONTROL_PATTERNS = [
    r"non[-_ ]?target",
    r"non[-_ ]?targeting",
    r"\bNTC\b",
    r"negative",
    r"neg[-_ ]?ctrl",
    r"control",
    r"scram",
    r"scrambled",
    r"safe",
    r"aavs1",
]

# ============================================================
# Branch-state definitions
# ============================================================

D2_P15_PROGRAM = 15
D2_P17_PROGRAM = 17

# D3 endpoint contrast:
# numerator = Fibro-Mesen
# denominator = artAC
#
# If your actual EC label is artEC instead of artAC, change this to:
# D3_EC_PATTERNS = [r"\bartEC\b"]

D3_EC_PATTERNS = [
    r"\bartEC\b",
]

D3_ALT_PATTERNS = [
    r"\bFibro-Mesen\b",
]

# ============================================================
# Analysis mode
# ============================================================

MODE = "all_carriers"
# MODE = "low_guide_load_leq2"
# MODE = "single_selected_target"

PSEUDOCOUNT = 0.5

MIN_GUIDE_D2_CELLS = 5
MIN_GUIDE_D3_CELLS = 5
MIN_CONTROL_D2_CELLS = 20
MIN_CONTROL_D3_CELLS = 20

# Plot
EFFECT_CLIP = 2.5
FDR_CUTOFF = 0.05
```

## Code cell 51

```python
def clean_day(x):
    s = str(x).lower().strip()
    s = s.replace("day", "d").replace(" ", "")

    if s in ["0", "0.0"]:
        return "D0"
    if s in ["1", "1.0"]:
        return "D1"
    if s in ["2", "2.0"]:
        return "D2"
    if s in ["3", "3.0"]:
        return "D3"

    return s.upper()


def matches_any(label, patterns):
    s = str(label)
    return any(re.search(p, s, flags=re.IGNORECASE) for p in patterns)


def classify_D3_endpoint_state(day, cell_type):
    if day != "D3":
        return None

    if matches_any(cell_type, D3_EC_PATTERNS):
        return "D3_EC"

    if matches_any(cell_type, D3_ALT_PATTERNS):
        return "D3_ALT"

    return None


def looks_control(x):
    s = str(x).lower()
    return any(re.search(p, s, flags=re.IGNORECASE) for p in CONTROL_PATTERNS)


def clean_target_string(x):
    if pd.isna(x):
        return []

    s = str(x).strip()

    if s.lower() in ["", "nan", "none", "na"]:
        return []

    s = (
        s.replace("[", "")
         .replace("]", "")
         .replace("'", "")
         .replace('"', "")
    )

    parts = re.split(r"[;,|/]+", s)
    parts = [p.strip() for p in parts if p.strip() != ""]

    return parts


def canonical_guide_id(x):
    return str(x).strip().replace(" ", "")


def get_guide_assignment_matrix(adata):
    ga = adata.obsm[GUIDE_ASSIGNMENT_KEY]

    if sparse.issparse(ga):
        return ga.tocsr().astype(bool)

    return np.asarray(ga).astype(bool)


def any_guides_per_cell(guide_assignment, guide_indices):
    if len(guide_indices) == 0:
        return np.zeros(guide_assignment.shape[0], dtype=bool)

    if sparse.issparse(guide_assignment):
        return np.asarray(guide_assignment[:, guide_indices].sum(axis=1)).ravel() > 0

    return guide_assignment[:, guide_indices].any(axis=1)


def load_guide_target_from_metadata(adata):
    guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))
    guide_name_to_idx = {g: i for i, g in enumerate(guide_names)}

    meta = pd.read_csv(GUIDE_METADATA_TSV, sep="\t")

    if GUIDE_ID_COL not in meta.columns:
        raise ValueError(f"{GUIDE_ID_COL!r} not found. Columns: {meta.columns.tolist()}")

    if TARGET_COL not in meta.columns:
        raise ValueError(f"{TARGET_COL!r} not found. Columns: {meta.columns.tolist()}")

    meta = meta[[GUIDE_ID_COL, TARGET_COL]].copy()
    meta.columns = ["guide_id", "putative_target_genes"]
    meta["guide_id"] = meta["guide_id"].astype(str)

    # Exact match to adata.uns["guide_names"]
    meta["guide_index"] = meta["guide_id"].map(guide_name_to_idx)
    meta["match_method"] = np.where(meta["guide_index"].notna(), "exact_guide_id_match", None)

    # Fallback canonical match
    if meta["guide_index"].isna().any():
        canonical_to_guide = {}
        duplicated = set()

        for g in guide_names:
            cg = canonical_guide_id(g)
            if cg in canonical_to_guide:
                duplicated.add(cg)
            else:
                canonical_to_guide[cg] = g

        for cg in duplicated:
            canonical_to_guide.pop(cg, None)

        for idx in meta.index[meta["guide_index"].isna()]:
            cg = canonical_guide_id(meta.at[idx, "guide_id"])
            if cg in canonical_to_guide:
                matched_guide = canonical_to_guide[cg]
                meta.at[idx, "guide_index"] = guide_name_to_idx[matched_guide]
                meta.at[idx, "match_method"] = "canonical_guide_id_match"

    mapped = meta[meta["guide_index"].notna()].copy()
    mapped["guide_index"] = mapped["guide_index"].astype(int)
    mapped["guide_name"] = [guide_names[i] for i in mapped["guide_index"]]

    records = []

    for _, row in mapped.iterrows():
        guide_index = int(row["guide_index"])
        guide_name = str(row["guide_name"])
        raw_target = row["putative_target_genes"]

        targets = clean_target_string(raw_target)

        if looks_control(guide_name) or looks_control(raw_target):
            targets = [CONTROL_TARGET]

        for target in targets:
            records.append({
                "guide_index": guide_index,
                "guide_name": guide_name,
                "target_name": str(target),
                "raw_target_entry": raw_target,
                "match_method": row["match_method"],
            })

    guide_target = pd.DataFrame(records).drop_duplicates()
    guide_target["target_name"] = guide_target["target_name"].astype(str)

    guide_target.to_csv(
        os.path.join(OUTDIR, "guide_to_target_map_from_300genes_metadata.csv"),
        index=False
    )

    print("\nGuide metadata mapping QC")
    print(f"adata guide columns:            {len(guide_names):,}")
    print(f"metadata guide rows:            {len(meta):,}")
    print(f"metadata rows matched to adata: {mapped['guide_id'].nunique():,}")
    print(f"expanded guide-target rows:     {len(guide_target):,}")
    print(f"unique mapped targets:          {guide_target['target_name'].nunique():,}")

    return guide_target


def get_program_series(adata, program_num):
    """
    Returns cNMF program usage for Program N.
    First tries adata.obs['Program N']; otherwise uses adata.obsm[CNMF_KEY].
    """
    obs_col = f"Program {program_num}"

    if obs_col in adata.obs.columns:
        return pd.to_numeric(adata.obs[obs_col], errors="coerce").values.astype(float)

    X = adata.obsm[CNMF_KEY]

    if sparse.issparse(X):
        X = X.toarray()

    X = np.asarray(X)
    idx = int(program_num) - int(CNMF_ARRAY_PROGRAM_INDEX_BASE)

    if idx < 0 or idx >= X.shape[1]:
        raise ValueError(
            f"Program {program_num} maps to column {idx}, but cnmf_usage has {X.shape[1]} columns. "
            "Check CNMF_ARRAY_PROGRAM_INDEX_BASE."
        )

    return X[:, idx].astype(float)


def log2_ratio_of_ratios(a, b, c, d, pseudocount=PSEUDOCOUNT):
    return np.log2(
        ((a + pseudocount) / (b + pseudocount)) /
        ((c + pseudocount) / (d + pseudocount))
    )


def fisher_branch_pval(a, b, c, d):
    try:
        return fisher_exact([[a, b], [c, d]], alternative="two-sided")[1]
    except Exception:
        return np.nan


def mannwhitney_pval(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]

    if len(x) == 0 or len(y) == 0:
        return np.nan

    try:
        return mannwhitneyu(x, y, alternative="two-sided")[1]
    except Exception:
        return np.nan


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan, dtype=float)

    valid = np.isfinite(pvals)

    if valid.sum() > 0:
        out[valid] = multipletests(pvals[valid], method="fdr_bh")[1]

    return out


def p_to_stars(p):
    if pd.isna(p):
        return "n.s."
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."
```

## Code cell 52

```python
guide_target = load_guide_target_from_metadata(adata)

guide_assignment = get_guide_assignment_matrix(adata)
guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))

if guide_assignment.shape[1] != len(guide_names):
    raise ValueError(
        f"guide_assignment has {guide_assignment.shape[1]} columns, "
        f"but adata.uns['{GUIDE_NAMES_KEY}'] has {len(guide_names)} names."
    )

# ============================================================
# Cell metadata
# ============================================================

obs = adata.obs.copy()
obs["cell_id"] = adata.obs_names.astype(str)
obs["day_clean"] = obs[DAY_COL].map(clean_day)
obs["cell_type_clean"] = obs[CELLTYPE_COL].astype(str)

obs["D3_endpoint_state"] = [
    classify_D3_endpoint_state(day, ct)
    for day, ct in zip(obs["day_clean"], obs["cell_type_clean"])
]

d3_state_vec = obs["D3_endpoint_state"].astype(object).values
d2_mask = obs["day_clean"].values == "D2"

print("\nD2 cell-type labels:")
print(sorted(obs.loc[obs["day_clean"] == "D2", "cell_type_clean"].unique()))

print("\nD3 cell-type labels:")
print(sorted(obs.loc[obs["day_clean"] == "D3", "cell_type_clean"].unique()))

print("\nD3 endpoint-state counts:")
print(obs["D3_endpoint_state"].value_counts(dropna=False))

if np.sum(d3_state_vec == "D3_EC") == 0:
    raise ValueError(
        "No D3_EC cells found. Check D3_EC_PATTERNS. "
        "If your label is artEC, set D3_EC_PATTERNS = [r'\\bartEC\\b']."
    )

if np.sum(d3_state_vec == "D3_ALT") == 0:
    raise ValueError(
        "No D3_ALT cells found. Check D3_ALT_PATTERNS. "
        "If your label differs from Fibro-Mesen, update D3_ALT_PATTERNS."
    )

# ============================================================
# Guide load and control masks
# ============================================================

if sparse.issparse(guide_assignment):
    n_total_guides = np.asarray(guide_assignment.sum(axis=1)).ravel().astype(int)
else:
    n_total_guides = guide_assignment.sum(axis=1).astype(int)

control_guides = guide_target[
    guide_target["target_name"].apply(lambda x: looks_control(x) or str(x) == CONTROL_TARGET)
].copy()

control_guide_indices = sorted(set(control_guides["guide_index"].astype(int).tolist()))
control_mask_raw = any_guides_per_cell(guide_assignment, control_guide_indices)

targeting_guide_target = guide_target[
    ~guide_target["target_name"].apply(lambda x: looks_control(x) or str(x) == CONTROL_TARGET)
].copy()

targeting_target_to_guides = (
    targeting_guide_target
    .groupby("target_name")["guide_index"]
    .apply(lambda x: sorted(set(x.astype(int).tolist())))
    .to_dict()
)

n_targeting_targets = np.zeros(adata.n_obs, dtype=int)

for target, idx in targeting_target_to_guides.items():
    n_targeting_targets += any_guides_per_cell(guide_assignment, idx).astype(int)

print("\nControl guide QC")
print("Number of control guide columns:", len(control_guide_indices))
print("Control carrier cells:", control_mask_raw.sum())

print("\nGuide load summary")
print("Median total guides:", np.median(n_total_guides))
print("Median targeting targets:", np.median(n_targeting_targets))

# ============================================================
# Mode filters
# ============================================================

def apply_mode_to_target_mask(mask):
    if MODE == "all_carriers":
        return mask.copy()

    if MODE == "low_guide_load_leq2":
        return mask & (n_total_guides <= 2)

    if MODE == "single_selected_target":
        return mask & (n_targeting_targets == 1)

    raise ValueError(f"Unknown MODE: {MODE}")


def apply_mode_to_control_guide_mask(mask):
    if MODE == "all_carriers":
        return mask.copy()

    if MODE == "low_guide_load_leq2":
        return mask & (n_total_guides <= 2)

    if MODE == "single_selected_target":
        return mask & (n_targeting_targets == 0)

    raise ValueError(f"Unknown MODE: {MODE}")


def pooled_control_mask(exclude_mask=None):
    if exclude_mask is None:
        exclude_mask = np.zeros(adata.n_obs, dtype=bool)

    if MODE == "all_carriers":
        return control_mask_raw & (~exclude_mask)

    if MODE == "low_guide_load_leq2":
        return control_mask_raw & (~exclude_mask) & (n_total_guides <= 2)

    if MODE == "single_selected_target":
        return control_mask_raw & (~exclude_mask) & (n_targeting_targets == 0)

    raise ValueError(f"Unknown MODE: {MODE}")

# ============================================================
# D2 program-space branch score: P17_z − P15_z
# z-scored using D2 non-targeting controls as reference
# ============================================================

p15 = get_program_series(adata, D2_P15_PROGRAM)
p17 = get_program_series(adata, D2_P17_PROGRAM)

control_ref = pooled_control_mask() & d2_mask & np.isfinite(p15) & np.isfinite(p17)

if control_ref.sum() < 30:
    raise ValueError(f"Too few D2 control cells for program-score z-scoring: {control_ref.sum()}")

p15_mu = np.nanmean(p15[control_ref])
p15_sd = np.nanstd(p15[control_ref])
p17_mu = np.nanmean(p17[control_ref])
p17_sd = np.nanstd(p17[control_ref])

p15_sd = p15_sd if p15_sd > 0 else 1.0
p17_sd = p17_sd if p17_sd > 0 else 1.0

p15_z = (p15 - p15_mu) / p15_sd
p17_z = (p17 - p17_mu) / p17_sd

D2_program_branch_score = p17_z - p15_z
D2_program_branch_score[~d2_mask] = np.nan

print("\nD2 program branch score QC")
print("D2 control reference cells:", control_ref.sum())
print("Mean NT D2 branch score:", np.nanmean(D2_program_branch_score[control_ref]))
print("SD NT D2 branch score:", np.nanstd(D2_program_branch_score[control_ref]))
```

### Text output (stream)

```text

Guide metadata mapping QC
adata guide columns:            1,922
metadata guide rows:            1,926
metadata rows matched to adata: 1,922
expanded guide-target rows:     1,922
unique mapped targets:          299

D2 cell-type labels:
['Fibro-Mesen', 'Lat-Meso', 'Neu-Mesen', 'Neuro-iPSC', 'Others', 'Pre-EC', 'PrimStreak', 'artEC', 'iPSC']

D3 cell-type labels:
['Fibro-Mesen', 'Lat-Meso', 'Neu-Mesen', 'Neuro-iPSC', 'Others', 'Pre-EC', 'PrimStreak', 'artEC', 'iPSC']

D3 endpoint-state counts:
D3_endpoint_state
None      80935
D3_EC     28444
D3_ALT    12237
Name: count, dtype: int64

Control guide QC
Number of control guide columns: 120
Control carrier cells: 21671

Guide load summary
Median total guides: 4.0
Median targeting targets: 4.0

D2 program branch score QC
D2 control reference cells: 4306
Mean NT D2 branch score: -1.592368183949133e-16
SD NT D2 branch score: 1.6413297712483208

```

## Code cell 53

```python
guide_target = load_guide_target_from_metadata(adata)

guide_assignment = get_guide_assignment_matrix(adata)
guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))

if guide_assignment.shape[1] != len(guide_names):
    raise ValueError(
        f"guide_assignment has {guide_assignment.shape[1]} columns, "
        f"but adata.uns['{GUIDE_NAMES_KEY}'] has {len(guide_names)} names."
    )

# ============================================================
# Cell metadata
# ============================================================

obs = adata.obs.copy()
obs["cell_id"] = adata.obs_names.astype(str)
obs["day_clean"] = obs[DAY_COL].map(clean_day)
obs["cell_type_clean"] = obs[CELLTYPE_COL].astype(str)

obs["D3_endpoint_state"] = [
    classify_D3_endpoint_state(day, ct)
    for day, ct in zip(obs["day_clean"], obs["cell_type_clean"])
]

d3_state_vec = obs["D3_endpoint_state"].astype(object).values
d2_mask = obs["day_clean"].values == "D2"

print("\nD2 cell-type labels:")
print(sorted(obs.loc[obs["day_clean"] == "D2", "cell_type_clean"].unique()))

print("\nD3 cell-type labels:")
print(sorted(obs.loc[obs["day_clean"] == "D3", "cell_type_clean"].unique()))

print("\nD3 endpoint-state counts:")
print(obs["D3_endpoint_state"].value_counts(dropna=False))

if np.sum(d3_state_vec == "D3_EC") == 0:
    raise ValueError(
        "No D3_EC cells found. Check D3_EC_PATTERNS. "
        "If your label is artEC, set D3_EC_PATTERNS = [r'\\bartEC\\b']."
    )

if np.sum(d3_state_vec == "D3_ALT") == 0:
    raise ValueError(
        "No D3_ALT cells found. Check D3_ALT_PATTERNS. "
        "If your label differs from Fibro-Mesen, update D3_ALT_PATTERNS."
    )

# ============================================================
# Guide load and control masks
# ============================================================

if sparse.issparse(guide_assignment):
    n_total_guides = np.asarray(guide_assignment.sum(axis=1)).ravel().astype(int)
else:
    n_total_guides = guide_assignment.sum(axis=1).astype(int)

control_guides = guide_target[
    guide_target["target_name"].apply(lambda x: looks_control(x) or str(x) == CONTROL_TARGET)
].copy()

control_guide_indices = sorted(set(control_guides["guide_index"].astype(int).tolist()))
control_mask_raw = any_guides_per_cell(guide_assignment, control_guide_indices)

targeting_guide_target = guide_target[
    ~guide_target["target_name"].apply(lambda x: looks_control(x) or str(x) == CONTROL_TARGET)
].copy()

targeting_target_to_guides = (
    targeting_guide_target
    .groupby("target_name")["guide_index"]
    .apply(lambda x: sorted(set(x.astype(int).tolist())))
    .to_dict()
)

n_targeting_targets = np.zeros(adata.n_obs, dtype=int)

for target, idx in targeting_target_to_guides.items():
    n_targeting_targets += any_guides_per_cell(guide_assignment, idx).astype(int)

print("\nControl guide QC")
print("Number of control guide columns:", len(control_guide_indices))
print("Control carrier cells:", control_mask_raw.sum())

print("\nGuide load summary")
print("Median total guides:", np.median(n_total_guides))
print("Median targeting targets:", np.median(n_targeting_targets))

# ============================================================
# Mode filters
# ============================================================

def apply_mode_to_target_mask(mask):
    if MODE == "all_carriers":
        return mask.copy()

    if MODE == "low_guide_load_leq2":
        return mask & (n_total_guides <= 2)

    if MODE == "single_selected_target":
        return mask & (n_targeting_targets == 1)

    raise ValueError(f"Unknown MODE: {MODE}")


def apply_mode_to_control_guide_mask(mask):
    if MODE == "all_carriers":
        return mask.copy()

    if MODE == "low_guide_load_leq2":
        return mask & (n_total_guides <= 2)

    if MODE == "single_selected_target":
        return mask & (n_targeting_targets == 0)

    raise ValueError(f"Unknown MODE: {MODE}")


def pooled_control_mask(exclude_mask=None):
    if exclude_mask is None:
        exclude_mask = np.zeros(adata.n_obs, dtype=bool)

    if MODE == "all_carriers":
        return control_mask_raw & (~exclude_mask)

    if MODE == "low_guide_load_leq2":
        return control_mask_raw & (~exclude_mask) & (n_total_guides <= 2)

    if MODE == "single_selected_target":
        return control_mask_raw & (~exclude_mask) & (n_targeting_targets == 0)

    raise ValueError(f"Unknown MODE: {MODE}")

# ============================================================
# D2 program-space branch score: P17_z − P15_z
# z-scored using D2 non-targeting controls as reference
# ============================================================

p15 = get_program_series(adata, D2_P15_PROGRAM)
p17 = get_program_series(adata, D2_P17_PROGRAM)

control_ref = pooled_control_mask() & d2_mask & np.isfinite(p15) & np.isfinite(p17)

if control_ref.sum() < 30:
    raise ValueError(f"Too few D2 control cells for program-score z-scoring: {control_ref.sum()}")

p15_mu = np.nanmean(p15[control_ref])
p15_sd = np.nanstd(p15[control_ref])
p17_mu = np.nanmean(p17[control_ref])
p17_sd = np.nanstd(p17[control_ref])

p15_sd = p15_sd if p15_sd > 0 else 1.0
p17_sd = p17_sd if p17_sd > 0 else 1.0

p15_z = (p15 - p15_mu) / p15_sd
p17_z = (p17 - p17_mu) / p17_sd

D2_program_branch_score = p17_z - p15_z
D2_program_branch_score[~d2_mask] = np.nan

print("\nD2 program branch score QC")
print("D2 control reference cells:", control_ref.sum())
print("Mean NT D2 branch score:", np.nanmean(D2_program_branch_score[control_ref]))
print("SD NT D2 branch score:", np.nanstd(D2_program_branch_score[control_ref]))
```

### Text output (stream)

```text

Guide metadata mapping QC
adata guide columns:            1,922
metadata guide rows:            1,926
metadata rows matched to adata: 1,922
expanded guide-target rows:     1,922
unique mapped targets:          299

D2 cell-type labels:
['Fibro-Mesen', 'Lat-Meso', 'Neu-Mesen', 'Neuro-iPSC', 'Others', 'Pre-EC', 'PrimStreak', 'artEC', 'iPSC']

D3 cell-type labels:
['Fibro-Mesen', 'Lat-Meso', 'Neu-Mesen', 'Neuro-iPSC', 'Others', 'Pre-EC', 'PrimStreak', 'artEC', 'iPSC']

D3 endpoint-state counts:
D3_endpoint_state
None      80935
D3_EC     28444
D3_ALT    12237
Name: count, dtype: int64

```

### Text output (stream)

```text

Control guide QC
Number of control guide columns: 120
Control carrier cells: 21671

Guide load summary
Median total guides: 4.0
Median targeting targets: 4.0

D2 program branch score QC
D2 control reference cells: 4306
Mean NT D2 branch score: -1.592368183949133e-16
SD NT D2 branch score: 1.6413297712483208

```

## Code cell 54

```python
guide_effects = pd.read_csv(GUIDE_EFFECTS_OUT)
target_effects = pd.read_csv(TARGET_EFFECTS_OUT)

# Keep only NT control + Hippo targets
guide_effects = guide_effects[guide_effects["requested_target"].isin(TARGET_ORDER)].copy()
target_effects = target_effects[target_effects["requested_target"].isin(TARGET_ORDER)].copy()

# y positions: top to bottom
target_y_map = {
    target: len(TARGET_ORDER) - i
    for i, target in enumerate(TARGET_ORDER)
}

for df in [guide_effects, target_effects]:
    df["y"] = df["requested_target"].map(target_y_map)
    df["target_group"] = df["requested_target"].map(TARGET_TO_GROUP)

# Long tables for two metrics
METRIC_D2 = "D2 program-space branch effect\nΔ mean cNMF score (P17 − P15)"
METRIC_D3 = "D3 cell-level endpoint effect\nlog2OR(Fibro-Mesen / artEC)"

guide_long = pd.concat([
    guide_effects.assign(
        metric=METRIC_D2,
        effect=guide_effects["d2_program_effect_delta_score"],
        pval=guide_effects["d2_pval"],
        fdr=guide_effects["d2_fdr"],
    ),
    guide_effects.assign(
        metric=METRIC_D3,
        effect=guide_effects["d3_cell_effect_log2OR"],
        pval=guide_effects["d3_pval"],
        fdr=guide_effects["d3_fdr"],
    ),
], ignore_index=True)

target_long = pd.concat([
    target_effects.assign(
        metric=METRIC_D2,
        effect=target_effects["d2_program_effect_delta_score"],
        pval=target_effects["d2_pval"],
        fdr=target_effects["d2_fdr"],
    ),
    target_effects.assign(
        metric=METRIC_D3,
        effect=target_effects["d3_cell_effect_log2OR"],
        pval=target_effects["d3_pval"],
        fdr=target_effects["d3_fdr"],
    ),
], ignore_index=True)

# Drop missing effects
guide_long = guide_long.dropna(subset=["effect", "y"]).copy()
target_long = target_long.dropna(subset=["effect", "y"]).copy()

guide_long["effect_plot"] = guide_long["effect"].clip(-EFFECT_CLIP, EFFECT_CLIP)
target_long["effect_plot"] = target_long["effect"].clip(-EFFECT_CLIP, EFFECT_CLIP)

# Target aggregate stars
target_long["star"] = ""
target_long.loc[target_long["fdr"] < 0.05, "star"] = "*"
target_long.loc[target_long["fdr"] < 0.01, "star"] = "**"
target_long.loc[target_long["fdr"] < 0.001, "star"] = "***"
target_long.loc[target_long["requested_target"] == "NT control", "star"] = ""

# Guide-level range by target and metric
range_df = (
    guide_long
    .groupby(["requested_target", "target_group", "metric", "y"], as_index=False)
    .agg(
        guide_min=("effect_plot", "min"),
        guide_max=("effect_plot", "max"),
        n_guides=("guide_name", "nunique"),
    )
)

# Separator between NT control and Hippo rows
separator_df = pd.DataFrame({
    "yintercept": [target_y_map["NT control"] - 0.5]
})

# ============================================================
# Hippo-vs-control statistical tests using guide-level effects
# ============================================================

stats_records = []

for metric in [METRIC_D2, METRIC_D3]:
    sub = guide_long[guide_long["metric"] == metric].copy()

    ctrl = sub[sub["requested_target"] == "NT control"]["effect"].dropna()
    hippo = sub[sub["requested_target"].isin(HIPPO_TARGETS)]["effect"].dropna()

    if len(ctrl) > 0 and len(hippo) > 0:
        pval = mannwhitneyu(hippo, ctrl, alternative="two-sided")[1]
        delta_median = hippo.median() - ctrl.median()
    else:
        pval = np.nan
        delta_median = np.nan

    stats_records.append({
        "metric": metric,
        "n_control_guides": len(ctrl),
        "n_hippo_guides": len(hippo),
        "control_median": ctrl.median() if len(ctrl) else np.nan,
        "hippo_median": hippo.median() if len(hippo) else np.nan,
        "delta_median": delta_median,
        "p_value": pval,
    })

hippo_vs_control = pd.DataFrame(stats_records)
hippo_vs_control["fdr"] = bh_fdr(hippo_vs_control["p_value"].values)
hippo_vs_control["stars"] = hippo_vs_control["fdr"].apply(p_to_stars)

hippo_vs_control["label"] = [
    f"Hippo vs NT guides: Δmedian={d:.2f}; FDR={q:.1e} {s}"
    for d, q, s in zip(
        hippo_vs_control["delta_median"],
        hippo_vs_control["fdr"],
        hippo_vs_control["stars"],
    )
]

hippo_vs_control["x"] = -EFFECT_CLIP + 0.10
hippo_vs_control["y"] = len(TARGET_ORDER) + 0.35

hippo_vs_control.to_csv(os.path.join(OUTDIR, "hippo_vs_control_guide_level_tests.csv"), index=False)

# Save plot sources
guide_long.to_csv(os.path.join(OUTDIR, "panelC1_guide_points_long.csv"), index=False)
target_long.to_csv(os.path.join(OUTDIR, "panelC1_target_points_long.csv"), index=False)
range_df.to_csv(os.path.join(OUTDIR, "panelC1_guide_range.csv"), index=False)

combined_source = pd.concat([
    guide_long.assign(point_type="individual_guide"),
    target_long.assign(point_type="target_aggregate"),
], ignore_index=True)

combined_source.to_csv(PLOT_SOURCE_OUT, index=False)

print("Saved plot source:", PLOT_SOURCE_OUT)
print("\nHippo vs control tests:")
display(hippo_vs_control)
```

### Text output (error)

```text
[31m---------------------------------------------------------------------------[39m
[31mFileNotFoundError[39m                         Traceback (most recent call last)
[36mCell[39m[36m [39m[32mIn[37][39m[32m, line 1[39m
[32m----> [39m[32m1[39m guide_effects = [43mpd[49m[43m.[49m[43mread_csv[49m[43m([49m[43mGUIDE_EFFECTS_OUT[49m[43m)[49m
[32m      2[39m target_effects = pd.read_csv(TARGET_EFFECTS_OUT)
[32m      4[39m [38;5;66;03m# Keep only NT control + Hippo targets[39;00m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/pandas/io/parsers/readers.py:1026[39m, in [36mread_csv[39m[34m(filepath_or_buffer, sep, delimiter, header, names, index_col, usecols, dtype, engine, converters, true_values, false_values, skipinitialspace, skiprows, skipfooter, nrows, na_values, keep_default_na, na_filter, verbose, skip_blank_lines, parse_dates, infer_datetime_format, keep_date_col, date_parser, date_format, dayfirst, cache_dates, iterator, chunksize, compression, thousands, decimal, lineterminator, quotechar, quoting, doublequote, escapechar, comment, encoding, encoding_errors, dialect, on_bad_lines, delim_whitespace, low_memory, memory_map, float_precision, storage_options, dtype_backend)[39m
[32m   1013[39m kwds_defaults = _refine_defaults_read(
[32m   1014[39m     dialect,
[32m   1015[39m     delimiter,
[32m   (...)[39m[32m   1022[39m     dtype_backend=dtype_backend,
[32m   1023[39m )
[32m   1024[39m kwds.update(kwds_defaults)
[32m-> [39m[32m1026[39m [38;5;28;01mreturn[39;00m [43m_read[49m[43m([49m[43mfilepath_or_buffer[49m[43m,[49m[43m [49m[43mkwds[49m[43m)[49m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/pandas/io/parsers/readers.py:620[39m, in [36m_read[39m[34m(filepath_or_buffer, kwds)[39m
[32m    617[39m _validate_names(kwds.get([33m"[39m[33mnames[39m[33m"[39m, [38;5;28;01mNone[39;00m))
[32m    619[39m [38;5;66;03m# Create the parser.[39;00m
[32m--> [39m[32m620[39m parser = [43mTextFileReader[49m[43m([49m[43mfilepath_or_buffer[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkwds[49m[43m)[49m
[32m    622[39m [38;5;28;01mif[39;00m chunksize [38;5;129;01mor[39;00m iterator:
[32m    623[39m     [38;5;28;01mreturn[39;00m parser

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/pandas/io/parsers/readers.py:1620[39m, in [36mTextFileReader.__init__[39m[34m(self, f, engine, **kwds)[39m
[32m   1617[39m     [38;5;28mself[39m.options[[33m"[39m[33mhas_index_names[39m[33m"[39m] = kwds[[33m"[39m[33mhas_index_names[39m[33m"[39m]
[32m   1619[39m [38;5;28mself[39m.handles: IOHandles | [38;5;28;01mNone[39;00m = [38;5;28;01mNone[39;00m
[32m-> [39m[32m1620[39m [38;5;28mself[39m._engine = [38;5;28;43mself[39;49m[43m.[49m[43m_make_engine[49m[43m([49m[43mf[49m[43m,[49m[43m [49m[38;5;28;43mself[39;49m[43m.[49m[43mengine[49m[43m)[49m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/pandas/io/parsers/readers.py:1880[39m, in [36mTextFileReader._make_engine[39m[34m(self, f, engine)[39m
[32m   1878[39m     [38;5;28;01mif[39;00m [33m"[39m[33mb[39m[33m"[39m [38;5;129;01mnot[39;00m [38;5;129;01min[39;00m mode:
[32m   1879[39m         mode += [33m"[39m[33mb[39m[33m"[39m
[32m-> [39m[32m1880[39m [38;5;28mself[39m.handles = [43mget_handle[49m[43m([49m
[32m   1881[39m [43m    [49m[43mf[49m[43m,[49m
[32m   1882[39m [43m    [49m[43mmode[49m[43m,[49m
[32m   1883[39m [43m    [49m[43mencoding[49m[43m=[49m[38;5;28;43mself[39;49m[43m.[49m[43moptions[49m[43m.[49m[43mget[49m[43m([49m[33;43m"[39;49m[33;43mencoding[39;49m[33;43m"[39;49m[43m,[49m[43m [49m[38;5;28;43;01mNone[39;49;00m[43m)[49m[43m,[49m
[32m   1884[39m [43m    [49m[43mcompression[49m[43m=[49m[38;5;28;43mself[39;49m[43m.[49m[43moptions[49m[43m.[49m[43mget[49m[43m([49m[33;43m"[39;49m[33;43mcompression[39;49m[33;43m"[39;49m[43m,[49m[43m [49m[38;5;28;43;01mNone[39;49;00m[43m)[49m[43m,[49m
[32m   1885[39m [43m    [49m[43mmemory_map[49m[43m=[49m[38;5;28;43mself[39;49m[43m.[49m[43moptions[49m[43m.[49m[43mget[49m[43m([49m[33;43m"[39;49m[33;43mmemory_map[39;49m[33;43m"[39;49m[43m,[49m[43m [49m[38;5;28;43;01mFalse[39;49;00m[43m)[49m[43m,[49m
[32m   1886[39m [43m    [49m[43mis_text[49m[43m=[49m[43mis_text[49m[43m,[49m
[32m   1887[39m [43m    [49m[43merrors[49m[43m=[49m[38;5;28;43mself[39;49m[43m.[49m[43moptions[49m[43m.[49m[43mget[49m[43m([49m[33;43m"[39;49m[33;43mencoding_errors[39;49m[33;43m"[39;49m[43m,[49m[43m [49m[33;43m"[39;49m[33;43mstrict[39;49m[33;43m"[39;49m[43m)[49m[43m,[49m
[32m   1888[39m [43m    [49m[43mstorage_options[49m[43m=[49m[38;5;28;43mself[39;49m[43m.[49m[43moptions[49m[43m.[49m[43mget[49m[43m([49m[33;43m"[39;49m[33;43mstorage_options[39;49m[33;43m"[39;49m[43m,[49m[43m [49m[38;5;28;43;01mNone[39;49;00m[43m)[49m[43m,[49m
[32m   1889[39m [43m[49m[43m)[49m
[32m   1890[39m [38;5;28;01massert[39;00m [38;5;28mself[39m.handles [38;5;129;01mis[39;00m [38;5;129;01mnot[39;00m [38;5;28;01mNone[39;00m
[32m   1891[39m f = [38;5;28mself[39m.handles.handle

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/pandas/io/common.py:873[39m, in [36mget_handle[39m[34m(path_or_buf, mode, encoding, compression, memory_map, is_text, errors, storage_options)[39m
[32m    868[39m [38;5;28;01melif[39;00m [38;5;28misinstance[39m(handle, [38;5;28mstr[39m):
[32m    869[39m     [38;5;66;03m# Check whether the filename is to be opened in binary mode.[39;00m
[32m    870[39m     [38;5;66;03m# Binary mode does not support 'encoding' and 'newline'.[39;00m
[32m    871[39m     [38;5;28;01mif[39;00m ioargs.encoding [38;5;129;01mand[39;00m [33m"[39m[33mb[39m[33m"[39m [38;5;129;01mnot[39;00m [38;5;129;01min[39;00m ioargs.mode:
[32m    872[39m         [38;5;66;03m# Encoding[39;00m
[32m--> [39m[32m873[39m         handle = [38;5;28mopen[39m(
[32m    874[39m             handle,
[32m    875[39m             ioargs.mode,
[32m    876[39m             encoding=ioargs.encoding,
[32m    877[39m             errors=errors,
[32m    878[39m             newline=[33m"[39m[33m"[39m,
[32m    879[39m         )
[32m    880[39m     [38;5;28;01melse[39;00m:
[32m    881[39m         [38;5;66;03m# Binary mode[39;00m
[32m    882[39m         handle = [38;5;28mopen[39m(handle, ioargs.mode)

[31mFileNotFoundError[39m: [Errno 2] No such file or directory: 'panelC1_guide_level_consistency_forest_from_adata/guide_level_D2_program_D3_cell_effects.csv'
```

## Code cell 55

```python
OUT_PNG = os.path.join(OUTDIR, "panelC1_guide_level_consistency_forest_plotnine.png")
OUT_PDF = os.path.join(OUTDIR, "panelC1_guide_level_consistency_forest_plotnine.pdf")
OUT_SVG = os.path.join(OUTDIR, "panelC1_guide_level_consistency_forest_plotnine.svg")

target_breaks = [target_y_map[t] for t in TARGET_ORDER]
target_labels = TARGET_ORDER

p = (
    ggplot()

    # separator between control and Hippo rows
    + geom_hline(
        data=separator_df,
        mapping=aes(yintercept="yintercept"),
        inherit_aes=False,
        linetype="dashed",
        color="#9E9E9E",
        size=0.8
    )

    # zero-effect line
    + geom_vline(
        xintercept=0,
        linetype="dashed",
        color="#9E9E9E",
        size=0.7
    )

    # guide range line
    + geom_segment(
        data=range_df,
        mapping=aes(x="guide_min", xend="guide_max", y="y", yend="y"),
        inherit_aes=False,
        color="#9E9E9E",
        size=1.0,
        alpha=0.8
    )

    # individual guide points
    + geom_point(
        data=guide_long,
        mapping=aes(x="effect_plot", y="y", fill="target_group"),
        inherit_aes=False,
        shape="o",
        size=2.9,
        color="#4D4D4D",
        stroke=0.35,
        alpha=0.88
    )

    # target aggregate diamonds
    + geom_point(
        data=target_long,
        mapping=aes(x="effect_plot", y="y", fill="target_group"),
        inherit_aes=False,
        shape="D",
        size=6.0,
        color="black",
        stroke=0.85,
        alpha=0.98
    )

    # target aggregate significance stars
    + geom_text(
        data=target_long[target_long["star"] != ""],
        mapping=aes(x="effect_plot", y="y", label="star"),
        inherit_aes=False,
        nudge_y=0.27,
        size=5.0,
        color="black",
        fontweight="bold"
    )

    # Hippo vs control statistical test annotation
    + geom_text(
        data=hippo_vs_control,
        mapping=aes(x="x", y="y", label="label"),
        inherit_aes=False,
        ha="left",
        size=3.8,
        color="black",
        fontweight="bold"
    )

    + facet_wrap("~ metric", nrow=1)

    + scale_y_continuous(
        breaks=target_breaks,
        labels=target_labels,
        limits=(0.5, len(TARGET_ORDER) + 0.65)
    )

    + scale_fill_manual(
        values=GROUP_COLORS,
        name="Guide group"
    )

    + coord_cartesian(
        xlim=(-EFFECT_CLIP, EFFECT_CLIP),
        expand=False
    )

    + labs(
        title="Guide-level consistency of Hippo branch and endpoint effects",
        subtitle="D2 uses program-space cNMF branch score; D3 uses cell-level Fibro-Mesen/artAC endpoint contrast",
        x="Effect size",
        y=None
    )

    + theme_bw(base_size=12)

    + theme(
        figure_size=(10.8, 4.8),

        panel_grid_major_y=element_blank(),
        panel_grid_minor=element_blank(),
        panel_grid_major_x=element_blank(),

        strip_background=element_rect(fill="#E6E6E6", color="#BDBDBD", size=0.8),
        strip_text=element_text(size=12, weight="bold"),

        axis_text_y=element_text(size=12),
        axis_text_x=element_text(size=11),
        axis_title_x=element_text(size=12, weight="bold"),

        plot_title=element_text(size=16, weight="bold"),
        plot_subtitle=element_text(size=10.5),

        legend_title=element_text(size=11, weight="bold"),
        legend_text=element_text(size=10),
        legend_position="right"
    )
)

p.save(OUT_PNG, width=10.8, height=4.8, dpi=350)
p.save(OUT_PDF, width=10.8, height=4.8)
p.save(OUT_SVG, width=10.8, height=4.8)

print("Saved:")
print(OUT_PNG)
print(OUT_PDF)
print(OUT_SVG)
```

## Markdown cell 56

#### Panel D Proportion

## Code cell 57

```python
import os
import re
import zipfile
import warnings
import numpy as np
import pandas as pd
import matplotlib as mpl

from scipy import sparse
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from plotnine import (
    ggplot, aes,
    geom_col, geom_text,
    scale_fill_manual, scale_color_identity,
    scale_y_continuous, scale_x_discrete,
    coord_flip, labs,
    theme_classic, theme,
    element_text, element_blank, element_rect, element_line,
    position_stack
)

warnings.filterwarnings("ignore")

# Helps remove thin white seams between stacked bar segments in PNG exports.
mpl.rcParams["savefig.dpi"] = 600
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["patch.antialiased"] = False

# ============================================================
# Output
# ============================================================

OUTDIR = "panelD_D2_D3_composition_three_controls_artEC"
os.makedirs(OUTDIR, exist_ok=True)

COMPOSITION_SOURCE_OUT = os.path.join(OUTDIR, "panelD_D2_D3_composition_source.csv")
STATS_SOURCE_OUT = os.path.join(OUTDIR, "panelD_D2_D3_composition_stats.csv")

# ============================================================
# adata fields
# ============================================================

DAY_COL = "day"
CELLTYPE_COL = "cell_type"

GUIDE_NAMES_KEY = "guide_names"
GUIDE_ASSIGNMENT_KEY = "guide_assignment"

# ============================================================
# guide metadata
# ============================================================

GUIDE_METADATA_TSV = "300genes_guide_metadata_v43.tsv"
GUIDE_ID_COL = "guide_id"
TARGET_COL = "putative_target_genes"

# ============================================================
# perturbation groups
# ============================================================

HIPPO_TARGETS = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]

CONTROL_TARGET = "non-targeting_control"

CONTROL_PATTERNS = [
    r"non[-_ ]?target",
    r"non[-_ ]?targeting",
    r"\bNTC\b",
    r"negative",
    r"neg[-_ ]?ctrl",
    r"control",
    r"scram",
    r"scrambled",
    r"safe",
    r"aavs1",
]

# ============================================================
# D2 state definitions
# ============================================================

D2_PRE_EC_PATTERNS = [
    r"\bPre-EC\b",
    r"\bP15\b",
    r"aEC progenitor",
]

D2_LAT_MESO_PATTERNS = [
    r"\bLat-Meso\b",
    r"\bP17\b",
    r"LM.*epithelium",
    r"LM-like",
]

# ============================================================
# D3 state definitions
# Actual label is artEC, not artAC
# ============================================================

D3_EC_PATTERNS = [
    r"\bartEC\b",
]

D3_FIBRO_PATTERNS = [
    r"\bFibro-Mesen\b",
]

# ============================================================
# Plot palette
# ============================================================

D2_STATE_COLORS = {
    "Pre-EC / P15": "#3F6FB3",
    "Lat-Meso / P17": "#CF4332",
    "Other D2": "#B8B8B8",
}

D3_STATE_COLORS = {
    "artEC / EC endpoint": "#3F6FB3",
    "Fibro-Mesen": "#CF4332",
    "Other D3": "#B8B8B8",
}

# Left-to-right order in the stacked bars.
D2_STATE_ORDER = ["Pre-EC / P15", "Lat-Meso / P17", "Other D2"]
D3_STATE_ORDER = ["artEC / EC endpoint", "Fibro-Mesen", "Other D3"]

GROUP_ORDER = ["Hippo targets", "NT excl. Hippo", "Strict NT"]
```

## Code cell 58

```python
def clean_day(x):
    s = str(x).lower().strip()
    s = s.replace("day", "d").replace(" ", "")

    if s in ["0", "0.0"]:
        return "D0"
    if s in ["1", "1.0"]:
        return "D1"
    if s in ["2", "2.0"]:
        return "D2"
    if s in ["3", "3.0"]:
        return "D3"

    return s.upper()


def matches_any(label, patterns):
    s = str(label)
    return any(re.search(p, s, flags=re.IGNORECASE) for p in patterns)


def looks_control(x):
    s = str(x).lower()
    return any(re.search(p, s, flags=re.IGNORECASE) for p in CONTROL_PATTERNS)


def clean_target_string(x):
    if pd.isna(x):
        return []

    s = str(x).strip()

    if s.lower() in ["", "nan", "none", "na"]:
        return []

    s = (
        s.replace("[", "")
         .replace("]", "")
         .replace("'", "")
         .replace('"', "")
    )

    parts = re.split(r"[;,|/]+", s)
    parts = [p.strip() for p in parts if p.strip() != ""]

    return parts


def canonical_guide_id(x):
    return str(x).strip().replace(" ", "")


def get_guide_assignment_matrix(adata):
    ga = adata.obsm[GUIDE_ASSIGNMENT_KEY]

    if sparse.issparse(ga):
        return ga.tocsr().astype(bool)

    return np.asarray(ga).astype(bool)


def any_guides_per_cell(guide_assignment, guide_indices):
    if len(guide_indices) == 0:
        return np.zeros(guide_assignment.shape[0], dtype=bool)

    if sparse.issparse(guide_assignment):
        return np.asarray(guide_assignment[:, guide_indices].sum(axis=1)).ravel() > 0

    return guide_assignment[:, guide_indices].any(axis=1)


def load_guide_target_from_metadata(adata):
    guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))
    guide_name_to_idx = {g: i for i, g in enumerate(guide_names)}

    meta = pd.read_csv(GUIDE_METADATA_TSV, sep="\t")

    if GUIDE_ID_COL not in meta.columns:
        raise ValueError(f"{GUIDE_ID_COL!r} not found. Columns: {meta.columns.tolist()}")

    if TARGET_COL not in meta.columns:
        raise ValueError(f"{TARGET_COL!r} not found. Columns: {meta.columns.tolist()}")

    meta = meta[[GUIDE_ID_COL, TARGET_COL]].copy()
    meta.columns = ["guide_id", "putative_target_genes"]
    meta["guide_id"] = meta["guide_id"].astype(str)

    # Exact matching to adata guide names.
    meta["guide_index"] = meta["guide_id"].map(guide_name_to_idx)
    meta["match_method"] = np.where(meta["guide_index"].notna(), "exact_guide_id_match", None)

    # Fallback canonical matching.
    if meta["guide_index"].isna().any():
        canonical_to_guide = {}
        duplicated = set()

        for g in guide_names:
            cg = canonical_guide_id(g)
            if cg in canonical_to_guide:
                duplicated.add(cg)
            else:
                canonical_to_guide[cg] = g

        for cg in duplicated:
            canonical_to_guide.pop(cg, None)

        for idx in meta.index[meta["guide_index"].isna()]:
            cg = canonical_guide_id(meta.at[idx, "guide_id"])
            if cg in canonical_to_guide:
                matched_guide = canonical_to_guide[cg]
                meta.at[idx, "guide_index"] = guide_name_to_idx[matched_guide]
                meta.at[idx, "match_method"] = "canonical_guide_id_match"

    mapped = meta[meta["guide_index"].notna()].copy()
    mapped["guide_index"] = mapped["guide_index"].astype(int)
    mapped["guide_name"] = [guide_names[i] for i in mapped["guide_index"]]

    records = []

    for _, row in mapped.iterrows():
        guide_index = int(row["guide_index"])
        guide_name = str(row["guide_name"])
        raw_target = row["putative_target_genes"]

        targets = clean_target_string(raw_target)

        if looks_control(guide_name) or looks_control(raw_target):
            targets = [CONTROL_TARGET]

        for target in targets:
            records.append({
                "guide_index": guide_index,
                "guide_name": guide_name,
                "target_name": str(target),
                "raw_target_entry": raw_target,
                "match_method": row["match_method"],
            })

    guide_target = pd.DataFrame(records).drop_duplicates()
    guide_target["target_name"] = guide_target["target_name"].astype(str)

    guide_target.to_csv(
        os.path.join(OUTDIR, "guide_to_target_map_from_300genes_metadata.csv"),
        index=False
    )

    print("\nGuide metadata mapping QC")
    print(f"adata guide columns:            {len(guide_names):,}")
    print(f"metadata guide rows:            {len(meta):,}")
    print(f"metadata rows matched to adata: {mapped['guide_id'].nunique():,}")
    print(f"expanded guide-target rows:     {len(guide_target):,}")
    print(f"unique mapped targets:          {guide_target['target_name'].nunique():,}")

    return guide_target


def classify_D2_state(cell_type):
    if matches_any(cell_type, D2_PRE_EC_PATTERNS):
        return "Pre-EC / P15"
    if matches_any(cell_type, D2_LAT_MESO_PATTERNS):
        return "Lat-Meso / P17"
    return "Other D2"


def classify_D3_state(cell_type):
    if matches_any(cell_type, D3_EC_PATTERNS):
        return "artEC / EC endpoint"
    if matches_any(cell_type, D3_FIBRO_PATTERNS):
        return "Fibro-Mesen"
    return "Other D3"


def p_to_stars(p):
    if pd.isna(p):
        return "n.s."
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def safe_fisher(table):
    try:
        return fisher_exact(table, alternative="two-sided")[1]
    except Exception:
        return np.nan
```

## Code cell 59

```python
# ============================================================
# Guide metadata and guide assignment
# ============================================================

guide_target = load_guide_target_from_metadata(adata)

guide_assignment = get_guide_assignment_matrix(adata)
guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))

if guide_assignment.shape[1] != len(guide_names):
    raise ValueError(
        f"guide_assignment has {guide_assignment.shape[1]} columns, "
        f"but adata.uns['{GUIDE_NAMES_KEY}'] has {len(guide_names)} names."
    )

# ============================================================
# Build target masks
# ============================================================

target_to_guide_indices = (
    guide_target
    .groupby("target_name")["guide_index"]
    .apply(lambda x: sorted(set(x.astype(int).tolist())))
    .to_dict()
)

target_masks = {
    target: any_guides_per_cell(guide_assignment, idx)
    for target, idx in target_to_guide_indices.items()
}

control_targets = [
    t for t in target_masks
    if looks_control(t) or str(t) == CONTROL_TARGET
]

if len(control_targets) == 0:
    raise ValueError("No control-like guide targets detected.")

control_mask_raw = np.zeros(adata.n_obs, dtype=bool)
for t in control_targets:
    control_mask_raw |= target_masks[t]

# Hippo mask: any Hippo target.
hippo_mask_raw = np.zeros(adata.n_obs, dtype=bool)

for target in HIPPO_TARGETS:
    matches = [t for t in target_masks if t.upper() == target.upper()]
    if len(matches) == 0:
        print(f"[warning] missing Hippo target: {target}")
        continue
    hippo_mask_raw |= target_masks[matches[0]]

# Targeting target count per cell, excluding NT controls.
targeting_targets = [
    t for t in target_masks
    if t not in control_targets and not looks_control(t)
]

n_targeting_targets = np.zeros(adata.n_obs, dtype=int)
for t in targeting_targets:
    n_targeting_targets += target_masks[t].astype(int)

# Total guide count per cell.
if sparse.issparse(guide_assignment):
    n_total_guides = np.asarray(guide_assignment.sum(axis=1)).ravel().astype(int)
else:
    n_total_guides = guide_assignment.sum(axis=1).astype(int)

# ============================================================
# Three comparison groups
# ============================================================

# 1. Hippo
hippo_mask = hippo_mask_raw.copy()

# 2. NT guide carriers excluding Hippo carriers.
# May carry other non-Hippo targeting guides.
nt_excl_hippo_mask = control_mask_raw & (~hippo_mask_raw)

# 3. Strict NT cells carrying NT guide and zero targeting target genes.
# Excludes all target genes.
strict_nt_mask = control_mask_raw & (n_targeting_targets == 0)

group_masks = {
    "Hippo targets": hippo_mask,
    "NT excl. Hippo": nt_excl_hippo_mask,
    "Strict NT": strict_nt_mask,
}

qc = pd.DataFrame({
    "group": list(group_masks.keys()),
    "n_cells": [int(m.sum()) for m in group_masks.values()],
})

print("\nGroup mask QC:")
display(qc)

print("\nOverlap QC:")
print("Hippo ∩ NT excl. Hippo:", int(np.sum(hippo_mask & nt_excl_hippo_mask)))
print("Hippo ∩ Strict NT:", int(np.sum(hippo_mask & strict_nt_mask)))

# ============================================================
# Cell metadata
# ============================================================

obs = adata.obs.copy()
obs["cell_id"] = adata.obs_names.astype(str)
obs["day_clean"] = obs[DAY_COL].map(clean_day)
obs["cell_type_clean"] = obs[CELLTYPE_COL].astype(str)

obs["D2_panel_state"] = obs["cell_type_clean"].apply(classify_D2_state)
obs["D3_panel_state"] = obs["cell_type_clean"].apply(classify_D3_state)

print("\nD2 cell-type labels:")
print(sorted(obs.loc[obs["day_clean"] == "D2", "cell_type_clean"].unique()))

print("\nD3 cell-type labels:")
print(sorted(obs.loc[obs["day_clean"] == "D3", "cell_type_clean"].unique()))

print("\nD2 selected state counts:")
print(obs.loc[obs["day_clean"] == "D2", "D2_panel_state"].value_counts())

print("\nD3 selected state counts:")
print(obs.loc[obs["day_clean"] == "D3", "D3_panel_state"].value_counts())
```

### Text output (stream)

```text

Guide metadata mapping QC
adata guide columns:            1,922
metadata guide rows:            1,926
metadata rows matched to adata: 1,922
expanded guide-target rows:     1,922
unique mapped targets:          299

Group mask QC:

```

### Text output (text/plain)

```text
            group  n_cells
0   Hippo targets    10813
1  NT excl. Hippo    18692
2       Strict NT      748
```

### Text output (stream)

```text

Overlap QC:
Hippo ∩ NT excl. Hippo: 0
Hippo ∩ Strict NT: 0

D2 cell-type labels:
['Fibro-Mesen', 'Lat-Meso', 'Neu-Mesen', 'Neuro-iPSC', 'Others', 'Pre-EC', 'PrimStreak', 'artEC', 'iPSC']

D3 cell-type labels:
['Fibro-Mesen', 'Lat-Meso', 'Neu-Mesen', 'Neuro-iPSC', 'Others', 'Pre-EC', 'PrimStreak', 'artEC', 'iPSC']

D2 selected state counts:
D2_panel_state
Lat-Meso / P17    16092
Pre-EC / P15       5765
Other D2            123
Name: count, dtype: int64

D3 selected state counts:
D3_panel_state
artEC / EC endpoint    28444
Fibro-Mesen            12237
Other D3                3421
Name: count, dtype: int64

```

## Code cell 60

```python
composition_records = []

for day_label, state_col in [
    ("D2", "D2_panel_state"),
    ("D3", "D3_panel_state"),
]:
    day_mask = obs["day_clean"].values == day_label

    for group_name, group_mask in group_masks.items():
        use = day_mask & group_mask

        tmp = obs.loc[use, ["cell_id", "cell_type_clean", state_col]].copy()
        tmp = tmp.rename(columns={state_col: "state"})
        tmp["day"] = day_label
        tmp["group"] = group_name

        composition_records.append(tmp)

composition_cells = pd.concat(composition_records, ignore_index=True)

counts = (
    composition_cells
    .groupby(["day", "group", "state"], as_index=False)
    .size()
    .rename(columns={"size": "n_cells"})
)

# Complete missing state/group/day combinations.
state_order_by_day = {
    "D2": D2_STATE_ORDER,
    "D3": D3_STATE_ORDER,
}

complete_records = []

for day in ["D2", "D3"]:
    for group in GROUP_ORDER:
        for state in state_order_by_day[day]:
            sub = counts[
                (counts["day"] == day)
                & (counts["group"] == group)
                & (counts["state"] == state)
            ]

            if len(sub) == 0:
                complete_records.append({
                    "day": day,
                    "group": group,
                    "state": state,
                    "n_cells": 0,
                })
            else:
                complete_records.append(sub.iloc[0].to_dict())

counts = pd.DataFrame(complete_records)

# Fractions
counts["total_cells"] = counts.groupby(["day", "group"])["n_cells"].transform("sum")
counts["fraction"] = counts["n_cells"] / counts["total_cells"].replace(0, np.nan)
counts["fraction"] = counts["fraction"].fillna(0)

counts["group"] = pd.Categorical(counts["group"], categories=GROUP_ORDER, ordered=True)

counts["state_order"] = counts.apply(
    lambda r: state_order_by_day[r["day"]].index(r["state"]),
    axis=1
)

counts = counts.sort_values(["day", "group", "state_order"]).copy()

counts.to_csv(COMPOSITION_SOURCE_OUT, index=False)

print("Saved composition source:", COMPOSITION_SOURCE_OUT)
display(counts)

# ============================================================
# Statistical tests
# ============================================================

def count_for(day, group, state):
    sub = counts[
        (counts["day"] == day)
        & (counts["group"].astype(str) == group)
        & (counts["state"] == state)
    ]

    if len(sub) == 0:
        return 0

    return int(sub["n_cells"].iloc[0])


stats_records = []

# D2 test: Lat-Meso/P17 vs Pre-EC/P15
for control_group in ["NT excl. Hippo", "Strict NT"]:
    h_lat = count_for("D2", "Hippo targets", "Lat-Meso / P17")
    h_pre = count_for("D2", "Hippo targets", "Pre-EC / P15")
    c_lat = count_for("D2", control_group, "Lat-Meso / P17")
    c_pre = count_for("D2", control_group, "Pre-EC / P15")

    p = safe_fisher([[h_lat, h_pre], [c_lat, c_pre]])

    stats_records.append({
        "day": "D2",
        "comparison": f"Hippo vs {control_group}",
        "test": "Fisher Lat-Meso/P17 vs Pre-EC/P15",
        "p_value": p,
        "stars": p_to_stars(p),
    })

# D3 test: Fibro-Mesen vs artEC
for control_group in ["NT excl. Hippo", "Strict NT"]:
    h_fib = count_for("D3", "Hippo targets", "Fibro-Mesen")
    h_ec = count_for("D3", "Hippo targets", "artEC / EC endpoint")
    c_fib = count_for("D3", control_group, "Fibro-Mesen")
    c_ec = count_for("D3", control_group, "artEC / EC endpoint")

    p = safe_fisher([[h_fib, h_ec], [c_fib, c_ec]])

    stats_records.append({
        "day": "D3",
        "comparison": f"Hippo vs {control_group}",
        "test": "Fisher Fibro-Mesen vs artEC",
        "p_value": p,
        "stars": p_to_stars(p),
    })

stats = pd.DataFrame(stats_records)
stats["fdr"] = multipletests(stats["p_value"].fillna(1), method="fdr_bh")[1]
stats["fdr_stars"] = stats["fdr"].apply(p_to_stars)

stats.to_csv(STATS_SOURCE_OUT, index=False)

print("Saved stats source:", STATS_SOURCE_OUT)
display(stats)
```

### Text output (stream)

```text
Saved composition source: panelD_D2_D3_composition_three_controls_artEC/panelD_D2_D3_composition_source.csv

```

### Text output (text/plain)

```text
   day           group                state  n_cells  total_cells  fraction  \
0   D2   Hippo targets         Pre-EC / P15      314         2157  0.145573   
1   D2   Hippo targets       Lat-Meso / P17     1824         2157  0.845619   
2   D2   Hippo targets             Other D2       19         2157  0.008809   
3   D2  NT excl. Hippo         Pre-EC / P15      822         3645  0.225514   
4   D2  NT excl. Hippo       Lat-Meso / P17     2811         3645  0.771193   
5   D2  NT excl. Hippo             Other D2       12         3645  0.003292   
6   D2       Strict NT         Pre-EC / P15       42           93  0.451613   
7   D2       Strict NT       Lat-Meso / P17       51           93  0.548387   
8   D2       Strict NT             Other D2        0           93  0.000000   
9   D3   Hippo targets  artEC / EC endpoint     1246         3753  0.332001   
10  D3   Hippo targets          Fibro-Mesen     2150         3753  0.572875   
11  D3   Hippo targets             Other D3      357         3753  0.095124   
12  D3  NT excl. Hippo  artEC / EC endpoint     3670         6308  0.581801   
13  D3  NT excl. Hippo          Fibro-Mesen     1969         6308  0.312143   
14  D3  NT excl. Hippo             Other D3      669         6308  0.106056   
15  D3       Strict NT  artEC / EC endpoint      339          378  0.896825   
16  D3       Strict NT          Fibro-Mesen       38          378  0.100529   
17  D3       Strict NT             Other D3        1          378  0.002646   

    state_order  
0             0  
1             1  
2             2  
3             0  
4             1  
5             2  
6             0  
7             1  
8             2  
9             0  
10            1  
11            2  
12            0  
13            1  
14            2  
15            0  
16            1  
17            2  
```

### Text output (stream)

```text
Saved stats source: panelD_D2_D3_composition_three_controls_artEC/panelD_D2_D3_composition_stats.csv

```

### Text output (text/plain)

```text
  day               comparison                               test  \
0  D2  Hippo vs NT excl. Hippo  Fisher Lat-Meso/P17 vs Pre-EC/P15   
1  D2       Hippo vs Strict NT  Fisher Lat-Meso/P17 vs Pre-EC/P15   
2  D3  Hippo vs NT excl. Hippo        Fisher Fibro-Mesen vs artEC   
3  D3       Hippo vs Strict NT        Fisher Fibro-Mesen vs artEC   

         p_value stars            fdr fdr_stars  
0   1.133447e-13   ***   1.511263e-13       ***  
1   8.893862e-12   ***   8.893862e-12       ***  
2  6.335505e-153   ***  2.534202e-152       ***  
3   1.737620e-93   ***   3.475241e-93       ***  
```

## Code cell 61

```python
# ============================================================
# Cell 1 — Illustrator-editable SVG/PDF settings
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import font_manager
from pathlib import Path
import logging

from plotnine import (
    ggplot, aes,
    geom_rect, geom_text,
    scale_fill_manual, scale_x_continuous, scale_y_continuous,
    coord_cartesian, labs,
    theme_classic, theme,
    element_text, element_blank, element_rect
)

# ------------------------------------------------------------
# Avoid repeated "Arial not found" warnings
# ------------------------------------------------------------

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

available_fonts = {f.name for f in font_manager.fontManager.ttflist}

FONT_CANDIDATES = [
    "Arial",
    "Helvetica",
    "Liberation Sans",
    "Nimbus Sans",
    "DejaVu Sans",
]

AI_FONT = None
for font in FONT_CANDIDATES:
    if font in available_fonts:
        AI_FONT = font
        break

if AI_FONT is None:
    AI_FONT = "DejaVu Sans"

print("Using font:", AI_FONT)

# ------------------------------------------------------------
# Illustrator-editable vector export settings
# ------------------------------------------------------------

mpl.rcParams.update({
    # SVG: keep text editable instead of converting text to paths.
    "svg.fonttype": "none",

    # PDF/EPS: keep text as TrueType instead of Type 3 paths.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,

    # Avoid LaTeX conversion of text into paths.
    "text.usetex": False,

    # Use installed sans-serif font.
    "font.family": "sans-serif",
    "font.sans-serif": [AI_FONT, "DejaVu Sans"],

    # High-quality export.
    "savefig.dpi": 600,
    "savefig.transparent": False,

    # Cleaner fills in vector output.
    "patch.antialiased": False,
})
```

### Text output (stream)

```text
Using font: DejaVu Sans

```

## Code cell 62

```python
# ============================================================
# Cell 2 — Illustrator-editable stacked-bar plotting function
# ============================================================

def make_publication_stacked_bar_clean(
    counts,
    stats,          # kept so old calls still work, but no FDR is plotted
    outdir,
    day,
    title,
    subtitle,
    state_colors,
    state_order,
    output_prefix,
):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    plot_df = counts[counts["day"] == day].copy()

    # --------------------------------------------------------
    # State and group ordering
    # --------------------------------------------------------
    plot_df["state"] = pd.Categorical(
        plot_df["state"],
        categories=state_order,
        ordered=True
    )

    plot_df["group"] = pd.Categorical(
        plot_df["group"],
        categories=["Hippo targets", "NT excl. Hippo", "Strict NT"],
        ordered=True
    )

    plot_df = plot_df.sort_values(["group", "state"]).copy()

    # --------------------------------------------------------
    # Group labels with n
    # --------------------------------------------------------
    total_by_group = (
        plot_df.groupby("group", observed=True)["n_cells"]
        .sum()
        .to_dict()
    )

    group_label_map = {
        "Hippo targets": f"Hippo targets\nn = {int(total_by_group.get('Hippo targets', 0)):,}",
        "NT excl. Hippo": f"NT excl. Hippo\nn = {int(total_by_group.get('NT excl. Hippo', 0)):,}",
        "Strict NT": f"Strict NT\nn = {int(total_by_group.get('Strict NT', 0)):,}",
    }

    # --------------------------------------------------------
    # Compact y positions
    # --------------------------------------------------------
    y_map = {
        "Hippo targets": 2.05,
        "NT excl. Hippo": 1.35,
        "Strict NT": 0.65,
    }

    BAR_HALF_HEIGHT = 0.30

    plot_df["y"] = plot_df["group"].astype(str).map(y_map)
    plot_df["group_label"] = plot_df["group"].astype(str).map(group_label_map)

    # --------------------------------------------------------
    # Explicit stacked-bar geometry in percent units
    # --------------------------------------------------------
    plot_df["fraction_pct"] = 100 * plot_df["fraction"]

    plot_df = plot_df.sort_values(["group", "state"]).copy()

    plot_df["xmin"] = (
        plot_df.groupby("group", observed=True)["fraction_pct"].cumsum()
        - plot_df["fraction_pct"]
    )
    plot_df["xmax"] = plot_df.groupby("group", observed=True)["fraction_pct"].cumsum()
    plot_df["xmid"] = (plot_df["xmin"] + plot_df["xmax"]) / 2

    plot_df["ymin"] = plot_df["y"] - BAR_HALF_HEIGHT
    plot_df["ymax"] = plot_df["y"] + BAR_HALF_HEIGHT

    # --------------------------------------------------------
    # Percent labels
    # --------------------------------------------------------
    plot_df["segment_label"] = np.where(
        plot_df["fraction_pct"] >= 8,
        plot_df["fraction_pct"].round(0).astype(int).astype(str) + "%",
        ""
    )

    plot_df["label_color"] = np.where(
        plot_df["state"].astype(str).str.contains("Other"),
        "#111111",
        "white"
    )

    label_white = plot_df[
        (plot_df["segment_label"] != "") &
        (plot_df["label_color"] == "white")
    ].copy()

    label_black = plot_df[
        (plot_df["segment_label"] != "") &
        (plot_df["label_color"] == "#111111")
    ].copy()

    # --------------------------------------------------------
    # Full-bar outlines only, no internal segment borders
    # --------------------------------------------------------
    outline_df = pd.DataFrame({
        "xmin": [0, 0, 0],
        "xmax": [100, 100, 100],
        "ymin": [
            y_map["Hippo targets"] - BAR_HALF_HEIGHT,
            y_map["NT excl. Hippo"] - BAR_HALF_HEIGHT,
            y_map["Strict NT"] - BAR_HALF_HEIGHT,
        ],
        "ymax": [
            y_map["Hippo targets"] + BAR_HALF_HEIGHT,
            y_map["NT excl. Hippo"] + BAR_HALF_HEIGHT,
            y_map["Strict NT"] + BAR_HALF_HEIGHT,
        ],
    })

    # --------------------------------------------------------
    # Output paths
    # --------------------------------------------------------
    out_png = outdir / f"{output_prefix}.png"
    out_pdf = outdir / f"{output_prefix}_AI_editable.pdf"
    out_svg = outdir / f"{output_prefix}_AI_editable.svg"

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------
    p = (
        ggplot()

        # Stacked segments, no internal border = no white seams.
        # These are vector rectangles in SVG/PDF.
        + geom_rect(
            data=plot_df,
            mapping=aes(
                xmin="xmin",
                xmax="xmax",
                ymin="ymin",
                ymax="ymax",
                fill="state"
            ),
            inherit_aes=False,
            color=None,
            size=0
        )

        # Full-bar outlines only.
        + geom_rect(
            data=outline_df,
            mapping=aes(
                xmin="xmin",
                xmax="xmax",
                ymin="ymin",
                ymax="ymax"
            ),
            inherit_aes=False,
            fill=None,
            color="#555555",
            size=0.55
        )

        # White labels.
        + geom_text(
            data=label_white,
            mapping=aes(x="xmid", y="y", label="segment_label"),
            inherit_aes=False,
            color="white",
            size=7.0,
            fontweight="bold",
            family=AI_FONT
        )

        # Black labels.
        + geom_text(
            data=label_black,
            mapping=aes(x="xmid", y="y", label="segment_label"),
            inherit_aes=False,
            color="#111111",
            size=7.0,
            fontweight="bold",
            family=AI_FONT
        )

        + scale_fill_manual(
            values=state_colors,
            limits=state_order,
            name=None
        )

        + scale_x_continuous(
            limits=(0, 100.5),
            breaks=[0, 25, 50, 75, 100],
            labels=["0", "25", "50", "75", "100"],
            expand=(0, 0)
        )

        + scale_y_continuous(
            breaks=[
                y_map["Hippo targets"],
                y_map["NT excl. Hippo"],
                y_map["Strict NT"],
            ],
            labels=[
                group_label_map["Hippo targets"],
                group_label_map["NT excl. Hippo"],
                group_label_map["Strict NT"],
            ],
            limits=(0.22, 2.48)
        )

        + coord_cartesian(
            xlim=(0, 100.5),
            ylim=(0.22, 2.48),
            expand=False
        )

        + labs(
            title=title,
            subtitle=subtitle,
            x="Cell fraction (%)",
            y=None
        )

        + theme_classic(base_size=11, base_family=AI_FONT)

        + theme(
            figure_size=(6.0, 2.55),

            panel_border=element_rect(
                fill=None,
                color="#6F6F6F",
                size=0.8
            ),

            panel_grid_major_x=element_blank(),
            panel_grid_major_y=element_blank(),
            panel_grid_minor=element_blank(),

            axis_text_x=element_text(
                size=10.5,
                color="#303030",
                family=AI_FONT
            ),
            axis_text_y=element_text(
                size=8.0,
                weight="bold",
                color="#303030",
                linespacing=0.95,
                family=AI_FONT
            ),

            axis_title_x=element_text(
                size=12.0,
                weight="bold",
                color="#202020",
                family=AI_FONT
            ),
            axis_title_y=element_blank(),

            plot_title=element_text(
                size=14.5,
                weight="bold",
                color="#111111",
                family=AI_FONT
            ),
            plot_subtitle=element_text(
                size=10.8,
                color="#222222",
                family=AI_FONT
            ),

            legend_position="top",
            legend_direction="horizontal",
            legend_title=element_blank(),
            legend_text=element_text(
                size=10.0,
                color="#222222",
                family=AI_FONT
            ),

            plot_margin=0.03
        )
    )

    # --------------------------------------------------------
    # Save outputs
    # --------------------------------------------------------
    p.save(str(out_png), width=6.0, height=2.55, dpi=600, verbose=False)
    p.save(str(out_pdf), width=6.0, height=2.55, verbose=False)
    p.save(str(out_svg), width=6.0, height=2.55, verbose=False)

    print("Saved:")
    print("  PNG:", out_png)
    print("  PDF:", out_pdf)
    print("  SVG:", out_svg)

    return p
```

## Code cell 63

```python
make_publication_stacked_bar_clean(
    counts=counts,
    stats=stats,
    outdir=OUTDIR,
    day="D2",
    title="D2 branch-state composition",
    subtitle="Pre-EC / P15 vs Lat-Meso / P17 among D2 cells",
    state_colors=D2_STATE_COLORS,
    state_order=D2_STATE_ORDER,
    output_prefix="panelD_D2_composition_three_controls_artEC_clean"
)
```

### Text output (stream)

```text
Saved:
  PNG: panelD_D2_D3_composition_three_controls_artEC/panelD_D2_composition_three_controls_artEC_clean.png
  PDF: panelD_D2_D3_composition_three_controls_artEC/panelD_D2_composition_three_controls_artEC_clean_AI_editable.pdf
  SVG: panelD_D2_D3_composition_three_controls_artEC/panelD_D2_composition_three_controls_artEC_clean_AI_editable.svg

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f2610833250>
```

## Code cell 64

```python
make_publication_stacked_bar_clean(
    counts=counts,
    stats=stats,
    outdir=OUTDIR,
    day="D3",
    title="D3 fate composition",
    subtitle="artEC vs Fibro-Mesen among D3 cells",
    state_colors=D3_STATE_COLORS,
    state_order=D3_STATE_ORDER,
    output_prefix="panelD_D3_composition_three_controls_artEC_clean"
)
```

### Text output (stream)

```text
Saved:
  PNG: panelD_D2_D3_composition_three_controls_artEC/panelD_D3_composition_three_controls_artEC_clean.png
  PDF: panelD_D2_D3_composition_three_controls_artEC/panelD_D3_composition_three_controls_artEC_clean_AI_editable.pdf
  SVG: panelD_D2_D3_composition_three_controls_artEC/panelD_D3_composition_three_controls_artEC_clean_AI_editable.svg

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f26009dc410>
```

## Markdown cell 65

#### Panel D - individual genes

## Code cell 66

```python
# ============================================================
# Panel D replacement — D3 EC vs non-EC by individual Hippo target
# ============================================================
# This cell assumes an AnnData object named `adata` is already loaded in the notebook.
# It replaces the old D2/D3 three-state pooled-Hippo Panel D section.
#
# What it does:
#   1) D3 only
#   2) Two states only: EC vs non-EC
#   3) Individual Hippo target bars: NF2, LATS2, SAV1, TAOK1, AMOTL2
#   4) Keeps the two control bars: NT excl. Hippo and Strict NT
#   5) Saves Illustrator-editable SVG/PDF and source CSVs

import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import font_manager
from scipy import sparse
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from plotnine import (
    ggplot, aes,
    geom_rect, geom_text,
    scale_fill_manual, scale_x_continuous, scale_y_continuous,
    coord_cartesian, labs,
    theme_classic, theme,
    element_text, element_blank, element_rect
)

warnings.filterwarnings("ignore")

# ============================================================
# Output
# ============================================================

OUTDIR = Path("panelD_D3_EC_vs_nonEC_individual_Hippo")
OUTDIR.mkdir(parents=True, exist_ok=True)

COMPOSITION_SOURCE_OUT = OUTDIR / "panelD_D3_EC_vs_nonEC_individual_Hippo_source.csv"
STATS_SOURCE_OUT = OUTDIR / "panelD_D3_EC_vs_nonEC_individual_Hippo_stats.csv"

# ============================================================
# AnnData fields
# ============================================================

DAY_COL = "day"
CELLTYPE_COL = "cell_type"

GUIDE_NAMES_KEY = "guide_names"
GUIDE_ASSIGNMENT_KEY = "guide_assignment"

# ============================================================
# Guide metadata
# ============================================================

GUIDE_METADATA_TSV = "300genes_guide_metadata_v43.tsv"
GUIDE_ID_COL = "guide_id"
# Robust to both old/new metadata column names.
TARGET_COL_CANDIDATES = ["putative_target_genes", "intended_target_name"]

# ============================================================
# Groups
# ============================================================

HIPPO_TARGETS = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]
CONTROL_TARGET = "non-targeting_control"

# If True, a target bar includes only cells carrying that Hippo target and no other Hippo target.
# If False, each bar is all carriers of that target. This matches the usual target-carrier definition.
EXCLUSIVE_HIPPO_TARGETS = False

CONTROL_PATTERNS = [
    r"non[-_ ]?target",
    r"non[-_ ]?targeting",
    r"\bNTC\b",
    r"negative",
    r"neg[-_ ]?ctrl",
    r"control",
    r"scram",
    r"scrambled",
    r"safe",
    r"aavs1",
]

# ============================================================
# D3 EC state definition
# ============================================================

D3_EC_PATTERNS = [
    r"\bartEC\b",
    # uncomment/add if your EC labels vary:
    # r"arterial.*EC",
    # r"endothelial",
]

D3_STATE_ORDER = ["artEC", "non-EC"]
D3_STATE_COLORS = {
    # artEC in red; non-EC in neutral grey.
    "artEC": "#C73D4C",
    "non-EC": "#B8B8B8",
}

GROUP_ORDER = HIPPO_TARGETS + ["NT excl. Hippo", "Strict NT"]

# ============================================================
# Illustrator-editable vector settings
# ============================================================

logging_level_font = {f.name for f in font_manager.fontManager.ttflist}
FONT_CANDIDATES = ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]
AI_FONT = next((f for f in FONT_CANDIDATES if f in logging_level_font), "DejaVu Sans")

mpl.rcParams.update({
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": [AI_FONT, "DejaVu Sans"],
    "savefig.dpi": 600,
    "savefig.transparent": False,
    "patch.antialiased": False,
})

print("Using font:", AI_FONT)

# ============================================================
# Helper functions
# ============================================================

def clean_day(x):
    s = str(x).lower().strip().replace("day", "d").replace(" ", "")
    if s in ["0", "0.0"]:
        return "D0"
    if s in ["1", "1.0"]:
        return "D1"
    if s in ["2", "2.0"]:
        return "D2"
    if s in ["3", "3.0"]:
        return "D3"
    return s.upper()


def matches_any(label, patterns):
    s = str(label)
    return any(re.search(p, s, flags=re.IGNORECASE) for p in patterns)


def looks_control(x):
    s = str(x).lower()
    return any(re.search(p, s, flags=re.IGNORECASE) for p in CONTROL_PATTERNS)


def clean_target_string(x):
    if pd.isna(x):
        return []
    s = str(x).strip()
    if s.lower() in ["", "nan", "none", "na"]:
        return []
    s = s.replace("[", "").replace("]", "").replace("'", "").replace('"', "")
    parts = re.split(r"[;,|/]+", s)
    return [p.strip() for p in parts if p.strip()]


def canonical_guide_id(x):
    return str(x).strip().replace(" ", "")


def get_guide_assignment_matrix(adata):
    ga = adata.obsm[GUIDE_ASSIGNMENT_KEY]
    if sparse.issparse(ga):
        return ga.tocsr().astype(bool)
    return np.asarray(ga).astype(bool)


def any_guides_per_cell(guide_assignment, guide_indices):
    if len(guide_indices) == 0:
        return np.zeros(guide_assignment.shape[0], dtype=bool)
    if sparse.issparse(guide_assignment):
        return np.asarray(guide_assignment[:, guide_indices].sum(axis=1)).ravel() > 0
    return guide_assignment[:, guide_indices].any(axis=1)


def classify_D3_EC_vs_nonEC(cell_type):
    if matches_any(cell_type, D3_EC_PATTERNS):
        return "artEC"
    return "non-EC"


def p_to_stars(p):
    if pd.isna(p):
        return "n.s."
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def safe_fisher(table):
    try:
        return fisher_exact(table, alternative="two-sided")[1]
    except Exception:
        return np.nan


def load_guide_target_from_metadata(adata):
    guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))
    guide_name_to_idx = {g: i for i, g in enumerate(guide_names)}

    meta = pd.read_csv(GUIDE_METADATA_TSV, sep="\t")
    if GUIDE_ID_COL not in meta.columns:
        raise ValueError(f"{GUIDE_ID_COL!r} not found. Columns: {meta.columns.tolist()}")

    target_col = None
    for c in TARGET_COL_CANDIDATES:
        if c in meta.columns:
            target_col = c
            break
    if target_col is None:
        raise ValueError(
            f"None of {TARGET_COL_CANDIDATES} found in {GUIDE_METADATA_TSV}. "
            f"Columns: {meta.columns.tolist()}"
        )

    meta = meta[[GUIDE_ID_COL, target_col]].copy()
    meta.columns = ["guide_id", "target_entry"]
    meta["guide_id"] = meta["guide_id"].astype(str)

    # Exact guide matching.
    meta["guide_index"] = meta["guide_id"].map(guide_name_to_idx)
    meta["match_method"] = np.where(meta["guide_index"].notna(), "exact_guide_id_match", None)

    # Fallback canonical matching.
    if meta["guide_index"].isna().any():
        canonical_to_guide = {}
        duplicated = set()
        for g in guide_names:
            cg = canonical_guide_id(g)
            if cg in canonical_to_guide:
                duplicated.add(cg)
            else:
                canonical_to_guide[cg] = g
        for cg in duplicated:
            canonical_to_guide.pop(cg, None)

        for idx in meta.index[meta["guide_index"].isna()]:
            cg = canonical_guide_id(meta.at[idx, "guide_id"])
            if cg in canonical_to_guide:
                matched = canonical_to_guide[cg]
                meta.at[idx, "guide_index"] = guide_name_to_idx[matched]
                meta.at[idx, "match_method"] = "canonical_guide_id_match"

    mapped = meta[meta["guide_index"].notna()].copy()
    mapped["guide_index"] = mapped["guide_index"].astype(int)
    mapped["guide_name"] = [guide_names[i] for i in mapped["guide_index"]]

    records = []
    for _, row in mapped.iterrows():
        guide_index = int(row["guide_index"])
        guide_name = str(row["guide_name"])
        raw_target = row["target_entry"]
        targets = clean_target_string(raw_target)
        if looks_control(guide_name) or looks_control(raw_target):
            targets = [CONTROL_TARGET]
        for target in targets:
            records.append({
                "guide_index": guide_index,
                "guide_name": guide_name,
                "target_name": str(target),
                "raw_target_entry": raw_target,
                "match_method": row["match_method"],
            })

    guide_target = pd.DataFrame(records).drop_duplicates()
    guide_target["target_name"] = guide_target["target_name"].astype(str)
    guide_target.to_csv(OUTDIR / "guide_to_target_map_for_panelD.csv", index=False)
    return guide_target


def count_for(counts, group, state):
    sub = counts[(counts["group"].astype(str) == group) & (counts["state"] == state)]
    if len(sub) == 0:
        return 0
    return int(sub["n_cells"].iloc[0])

# ============================================================
# Build guide/target masks
# ============================================================

guide_target = load_guide_target_from_metadata(adata)
guide_assignment = get_guide_assignment_matrix(adata)
guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))

if guide_assignment.shape[1] != len(guide_names):
    raise ValueError(
        f"guide_assignment has {guide_assignment.shape[1]} columns, "
        f"but adata.uns[{GUIDE_NAMES_KEY!r}] has {len(guide_names)} names."
    )

target_to_guide_indices = (
    guide_target.groupby("target_name")["guide_index"]
    .apply(lambda x: sorted(set(x.astype(int).tolist())))
    .to_dict()
)

target_masks = {
    target: any_guides_per_cell(guide_assignment, idx)
    for target, idx in target_to_guide_indices.items()
}

control_targets = [t for t in target_masks if looks_control(t) or str(t) == CONTROL_TARGET]
if len(control_targets) == 0:
    raise ValueError("No control-like guide targets detected.")

control_mask_raw = np.zeros(adata.n_obs, dtype=bool)
for t in control_targets:
    control_mask_raw |= target_masks[t]

# Hippo masks, by individual target.
hippo_target_masks = {}
missing_hippo = []
for target in HIPPO_TARGETS:
    matches = [t for t in target_masks if t.upper() == target.upper()]
    if len(matches) == 0:
        missing_hippo.append(target)
    else:
        hippo_target_masks[target] = target_masks[matches[0]].copy()

if missing_hippo:
    print("[warning] missing Hippo targets:", missing_hippo)

if len(hippo_target_masks) == 0:
    raise ValueError("No Hippo target masks found.")

hippo_mask_raw = np.zeros(adata.n_obs, dtype=bool)
for m in hippo_target_masks.values():
    hippo_mask_raw |= m

# Count targeting targets per cell, excluding NT controls.
targeting_targets = [
    t for t in target_masks
    if t not in control_targets and not looks_control(t)
]
n_targeting_targets = np.zeros(adata.n_obs, dtype=int)
for t in targeting_targets:
    n_targeting_targets += target_masks[t].astype(int)

# Control groups kept as before.
nt_excl_hippo_mask = control_mask_raw & (~hippo_mask_raw)
strict_nt_mask = control_mask_raw & (n_targeting_targets == 0)

# Individual Hippo groups.
group_masks = {}
for target, mask in hippo_target_masks.items():
    if EXCLUSIVE_HIPPO_TARGETS:
        other_hippo = np.zeros(adata.n_obs, dtype=bool)
        for other_target, other_mask in hippo_target_masks.items():
            if other_target != target:
                other_hippo |= other_mask
        group_masks[target] = mask & (~other_hippo)
    else:
        group_masks[target] = mask

group_masks["NT excl. Hippo"] = nt_excl_hippo_mask
group_masks["Strict NT"] = strict_nt_mask

qc = pd.DataFrame({
    "group": list(group_masks.keys()),
    "n_total_cells": [int(m.sum()) for m in group_masks.values()],
})
qc.to_csv(OUTDIR / "panelD_D3_EC_vs_nonEC_group_mask_qc.csv", index=False)
print("Group mask QC:")
display(qc)

# ============================================================
# Cell metadata and D3 EC/non-EC state
# ============================================================

obs = adata.obs.copy()
obs["cell_id"] = adata.obs_names.astype(str)
obs["day_clean"] = obs[DAY_COL].map(clean_day)
obs["cell_type_clean"] = obs[CELLTYPE_COL].astype(str)
obs["D3_EC_nonEC_state"] = obs["cell_type_clean"].apply(classify_D3_EC_vs_nonEC)

print("D3 cell-type labels:")
print(sorted(obs.loc[obs["day_clean"] == "D3", "cell_type_clean"].unique()))
print("D3 EC/non-EC counts:")
print(obs.loc[obs["day_clean"] == "D3", "D3_EC_nonEC_state"].value_counts())

# ============================================================
# Composition counts: D3 only, EC vs non-EC
# ============================================================

records = []
d3_day_mask = obs["day_clean"].values == "D3"

for group_name, group_mask in group_masks.items():
    use = d3_day_mask & group_mask
    tmp = obs.loc[use, ["cell_id", "cell_type_clean", "D3_EC_nonEC_state"]].copy()
    tmp = tmp.rename(columns={"D3_EC_nonEC_state": "state"})
    tmp["day"] = "D3"
    tmp["group"] = group_name
    records.append(tmp)

composition_cells = pd.concat(records, ignore_index=True)
composition_cells.to_csv(OUTDIR / "panelD_D3_EC_vs_nonEC_cell_membership.csv", index=False)

counts = (
    composition_cells.groupby(["day", "group", "state"], as_index=False)
    .size()
    .rename(columns={"size": "n_cells"})
)

# Complete missing group/state combinations.
complete_records = []
for group in GROUP_ORDER:
    if group not in group_masks:
        continue
    for state in D3_STATE_ORDER:
        sub = counts[(counts["group"] == group) & (counts["state"] == state)]
        if len(sub) == 0:
            complete_records.append({"day": "D3", "group": group, "state": state, "n_cells": 0})
        else:
            complete_records.append(sub.iloc[0].to_dict())

counts = pd.DataFrame(complete_records)
counts["total_cells"] = counts.groupby(["day", "group"])["n_cells"].transform("sum")
counts["fraction"] = counts["n_cells"] / counts["total_cells"].replace(0, np.nan)
counts["fraction"] = counts["fraction"].fillna(0)
counts["group"] = pd.Categorical(counts["group"], categories=[g for g in GROUP_ORDER if g in group_masks], ordered=True)
counts["state"] = pd.Categorical(counts["state"], categories=D3_STATE_ORDER, ordered=True)
counts = counts.sort_values(["group", "state"]).copy()
counts.to_csv(COMPOSITION_SOURCE_OUT, index=False)
print("Saved composition source:", COMPOSITION_SOURCE_OUT)
display(counts)

# ============================================================
# Fisher tests and FDR: each Hippo target vs each control group
# ============================================================

stats_records = []
for hippo_target in HIPPO_TARGETS:
    if hippo_target not in group_masks:
        continue
    for control_group in ["NT excl. Hippo", "Strict NT"]:
        h_ec = count_for(counts, hippo_target, "artEC")
        h_non = count_for(counts, hippo_target, "non-EC")
        c_ec = count_for(counts, control_group, "artEC")
        c_non = count_for(counts, control_group, "non-EC")
        p = safe_fisher([[h_ec, h_non], [c_ec, c_non]])
        stats_records.append({
            "day": "D3",
            "target": hippo_target,
            "comparison": f"{hippo_target} vs {control_group}",
            "test": "Fisher artEC vs non-EC",
            "target_artEC": h_ec,
            "target_nonEC": h_non,
            "control_artEC": c_ec,
            "control_nonEC": c_non,
            "p_value": p,
            "p_stars": p_to_stars(p),
        })

stats = pd.DataFrame(stats_records)
if len(stats):
    stats["fdr"] = multipletests(stats["p_value"].fillna(1), method="fdr_bh")[1]
    stats["fdr_stars"] = stats["fdr"].apply(p_to_stars)
else:
    stats["fdr"] = []
    stats["fdr_stars"] = []

stats.to_csv(STATS_SOURCE_OUT, index=False)
print("Saved stats source:", STATS_SOURCE_OUT)
display(stats)

# ============================================================
# Illustrator-editable stacked bar plot
# ============================================================

def make_panelD_D3_EC_vs_nonEC_individual_Hippo(counts, outdir=OUTDIR):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    plot_df = counts.copy()
    plot_df["group"] = pd.Categorical(plot_df["group"], categories=[g for g in GROUP_ORDER if g in group_masks], ordered=True)
    plot_df["state"] = pd.Categorical(plot_df["state"], categories=D3_STATE_ORDER, ordered=True)
    plot_df = plot_df.sort_values(["group", "state"]).copy()

    # Group labels with n.
    total_by_group = plot_df.groupby("group", observed=True)["n_cells"].sum().to_dict()
    group_order_present = [g for g in GROUP_ORDER if g in group_masks]

    group_label_map = {
        g: f"{g}\nn = {int(total_by_group.get(g, 0)):,}"
        for g in group_order_present
    }

    # Compact y positions: one row per group.
    # Smaller BAR_SPACING and larger BAR_HALF_HEIGHT reduce the white space
    # between horizontal bars.
    BAR_SPACING = 0.78
    BAR_HALF_HEIGHT = 0.34
    y_map = {
        g: (len(group_order_present) - i) * BAR_SPACING
        for i, g in enumerate(group_order_present)
    }
    y_min = min(y_map.values()) - BAR_HALF_HEIGHT - 0.16
    y_max = max(y_map.values()) + BAR_HALF_HEIGHT + 0.16

    plot_df["y"] = plot_df["group"].astype(str).map(y_map)
    plot_df["fraction_pct"] = 100 * plot_df["fraction"]
    plot_df["xmin"] = plot_df.groupby("group", observed=True)["fraction_pct"].cumsum() - plot_df["fraction_pct"]
    plot_df["xmax"] = plot_df.groupby("group", observed=True)["fraction_pct"].cumsum()
    plot_df["xmid"] = (plot_df["xmin"] + plot_df["xmax"]) / 2
    plot_df["ymin"] = plot_df["y"] - BAR_HALF_HEIGHT
    plot_df["ymax"] = plot_df["y"] + BAR_HALF_HEIGHT

    plot_df["segment_label"] = np.where(
        plot_df["fraction_pct"] >= 8,
        plot_df["fraction_pct"].round(0).astype(int).astype(str) + "%",
        ""
    )
    plot_df["label_color"] = np.where(plot_df["state"].astype(str) == "non-EC", "#111111", "white")

    label_white = plot_df[(plot_df["segment_label"] != "") & (plot_df["label_color"] == "white")].copy()
    label_black = plot_df[(plot_df["segment_label"] != "") & (plot_df["label_color"] == "#111111")].copy()

    outline_df = pd.DataFrame({
        "xmin": [0] * len(group_order_present),
        "xmax": [100] * len(group_order_present),
        "ymin": [y_map[g] - BAR_HALF_HEIGHT for g in group_order_present],
        "ymax": [y_map[g] + BAR_HALF_HEIGHT for g in group_order_present],
    })

    out_png = outdir / "panelD_D3_EC_vs_nonEC_individual_Hippo.png"
    out_pdf = outdir / "panelD_D3_EC_vs_nonEC_individual_Hippo_AI_editable.pdf"
    out_svg = outdir / "panelD_D3_EC_vs_nonEC_individual_Hippo_AI_editable.svg"

    p = (
        ggplot()
        + geom_rect(
            data=plot_df,
            mapping=aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="state"),
            inherit_aes=False,
            color=None,
            size=0
        )
        + geom_rect(
            data=outline_df,
            mapping=aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax"),
            inherit_aes=False,
            fill=None,
            color="#555555",
            size=0.50
        )
        + geom_text(
            data=label_white,
            mapping=aes(x="xmid", y="y", label="segment_label"),
            inherit_aes=False,
            color="white",
            size=6.5,
            fontweight="bold",
            family=AI_FONT
        )
        + geom_text(
            data=label_black,
            mapping=aes(x="xmid", y="y", label="segment_label"),
            inherit_aes=False,
            color="#111111",
            size=6.5,
            fontweight="bold",
            family=AI_FONT
        )
        + scale_fill_manual(values=D3_STATE_COLORS, limits=D3_STATE_ORDER, name=None)
        + scale_x_continuous(
            limits=(0, 100.5),
            breaks=[0, 25, 50, 75, 100],
            labels=["0", "25", "50", "75", "100"],
            expand=(0, 0)
        )
        + scale_y_continuous(
            breaks=[y_map[g] for g in group_order_present],
            labels=[group_label_map[g] for g in group_order_present],
            limits=(y_min, y_max)
        )
        + coord_cartesian(xlim=(0, 100.5), ylim=(y_min, y_max), expand=False)
        + labs(
            title="D3 EC composition by Hippo target",
            subtitle="artEC vs non-EC among D3 cells",
            x="Cell fraction (%)",
            y=None
        )
        + theme_classic(base_size=11, base_family=AI_FONT)
        + theme(
            figure_size=(6.2, max(2.85, 0.42 * len(group_order_present) + 0.95)),
            panel_border=element_rect(fill=None, color="#6F6F6F", size=0.8),
            panel_grid_major_x=element_blank(),
            panel_grid_major_y=element_blank(),
            panel_grid_minor=element_blank(),
            axis_text_x=element_text(size=10.5, color="#303030", family=AI_FONT),
            axis_text_y=element_text(size=8.3, weight="bold", color="#303030", linespacing=0.95, family=AI_FONT),
            axis_title_x=element_text(size=12.0, weight="bold", color="#202020", family=AI_FONT),
            axis_title_y=element_blank(),
            plot_title=element_text(size=14.2, weight="bold", color="#111111", family=AI_FONT),
            plot_subtitle=element_text(size=10.6, color="#222222", family=AI_FONT),
            legend_position="top",
            legend_direction="horizontal",
            legend_title=element_blank(),
            legend_text=element_text(size=10.0, color="#222222", family=AI_FONT),
            plot_margin=0.03
        )
    )

    width = 6.2
    height = max(2.85, 0.42 * len(group_order_present) + 0.95)
    p.save(str(out_png), width=width, height=height, dpi=600, verbose=False)
    p.save(str(out_pdf), width=width, height=height, verbose=False)
    p.save(str(out_svg), width=width, height=height, verbose=False)

    print("Saved:")
    print("  PNG:", out_png)
    print("  PDF:", out_pdf)
    print("  SVG:", out_svg)

    return p

panelD_plot = make_panelD_D3_EC_vs_nonEC_individual_Hippo(counts, OUTDIR)

```

### Text output (stream)

```text
Using font: DejaVu Sans
Group mask QC:

```

### Text output (text/plain)

```text
            group  n_total_cells
0             NF2           2354
1           LATS2           2695
2            SAV1           1595
3           TAOK1           2530
4          AMOTL2           2381
5  NT excl. Hippo          18692
6       Strict NT            748
```

### Text output (stream)

```text
D3 cell-type labels:
['Fibro-Mesen', 'Lat-Meso', 'Neu-Mesen', 'Neuro-iPSC', 'Others', 'Pre-EC', 'PrimStreak', 'artEC', 'iPSC']
D3 EC/non-EC counts:
D3_EC_nonEC_state
artEC     28444
non-EC    15658
Name: count, dtype: int64
Saved composition source: panelD_D3_EC_vs_nonEC_individual_Hippo/panelD_D3_EC_vs_nonEC_individual_Hippo_source.csv

```

### Text output (text/plain)

```text
   day           group   state  n_cells  total_cells  fraction
0   D3             NF2   artEC      190          823  0.230863
1   D3             NF2  non-EC      633          823  0.769137
2   D3           LATS2   artEC      292          996  0.293173
3   D3           LATS2  non-EC      704          996  0.706827
4   D3            SAV1   artEC      233          568  0.410211
5   D3            SAV1  non-EC      335          568  0.589789
6   D3           TAOK1   artEC      318          862  0.368910
7   D3           TAOK1  non-EC      544          862  0.631090
8   D3          AMOTL2   artEC      263          810  0.324691
9   D3          AMOTL2  non-EC      547          810  0.675309
10  D3  NT excl. Hippo   artEC     3670         6308  0.581801
11  D3  NT excl. Hippo  non-EC     2638         6308  0.418199
12  D3       Strict NT   artEC      339          378  0.896825
13  D3       Strict NT  non-EC       39          378  0.103175
```

### Text output (stream)

```text
Saved stats source: panelD_D3_EC_vs_nonEC_individual_Hippo/panelD_D3_EC_vs_nonEC_individual_Hippo_stats.csv

```

### Text output (text/plain)

```text
  day  target                comparison                    test  target_artEC  \
0  D3     NF2     NF2 vs NT excl. Hippo  Fisher artEC vs non-EC           190   
1  D3     NF2          NF2 vs Strict NT  Fisher artEC vs non-EC           190   
2  D3   LATS2   LATS2 vs NT excl. Hippo  Fisher artEC vs non-EC           292   
3  D3   LATS2        LATS2 vs Strict NT  Fisher artEC vs non-EC           292   
4  D3    SAV1    SAV1 vs NT excl. Hippo  Fisher artEC vs non-EC           233   
5  D3    SAV1         SAV1 vs Strict NT  Fisher artEC vs non-EC           233   
6  D3   TAOK1   TAOK1 vs NT excl. Hippo  Fisher artEC vs non-EC           318   
7  D3   TAOK1        TAOK1 vs Strict NT  Fisher artEC vs non-EC           318   
8  D3  AMOTL2  AMOTL2 vs NT excl. Hippo  Fisher artEC vs non-EC           263   
9  D3  AMOTL2       AMOTL2 vs Strict NT  Fisher artEC vs non-EC           263   

   target_nonEC  control_artEC  control_nonEC        p_value p_stars  \
0           633           3670           2638   8.715117e-83     ***   
1           633            339             39  6.364833e-112     ***   
2           704           3670           2638   1.766423e-65     ***   
3           704            339             39   3.547358e-97     ***   
4           335           3670           2638   4.472645e-15     ***   
5           335            339             39   6.532900e-56     ***   
6           544           3670           2638   7.372393e-32     ***   
7           544            339             39   5.745777e-73     ***   
8           547           3670           2638   1.317346e-43     ***   
9           547            339             39   6.267596e-83     ***   

             fdr fdr_stars  
0   2.178779e-82       ***  
1  6.364833e-111       ***  
2   2.944039e-65       ***  
3   1.773679e-96       ***  
4   4.472645e-15       ***  
5   9.332715e-56       ***  
6   8.191548e-32       ***  
7   1.149155e-72       ***  
8   1.646682e-43       ***  
9   2.089199e-82       ***  
```

### Text output (stream)

```text
Saved:
  PNG: panelD_D3_EC_vs_nonEC_individual_Hippo/panelD_D3_EC_vs_nonEC_individual_Hippo.png
  PDF: panelD_D3_EC_vs_nonEC_individual_Hippo/panelD_D3_EC_vs_nonEC_individual_Hippo_AI_editable.pdf
  SVG: panelD_D3_EC_vs_nonEC_individual_Hippo/panelD_D3_EC_vs_nonEC_individual_Hippo_AI_editable.svg

```

## Markdown cell 67

#### Scatterplot grid

## Code cell 68

```python
# ============================================================
# Export cell-level target outcome inputs for Jesse-style scatter plots
# ============================================================

import os
import re
import zipfile
import warnings
import numpy as np
import pandas as pd

from scipy import sparse
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ============================================================
# Settings
# ============================================================

OUTDIR = "scatter_plot_celllevel_endpoint_exports"
os.makedirs(OUTDIR, exist_ok=True)

GUIDE_METADATA_TSV = "300genes_guide_metadata_v43.tsv"

DAY_COL = "day"
CELLTYPE_COL = "cell_type"

GUIDE_NAMES_KEY = "guide_names"
GUIDE_ASSIGNMENT_KEY = "guide_assignment"

GUIDE_ID_COL = "guide_id"
TARGET_COL = "putative_target_genes"

CONTROL_TARGET = "non-targeting_control"

CONTROL_PATTERNS = [
    r"non[-_ ]?target",
    r"non[-_ ]?targeting",
    r"\bNTC\b",
    r"negative",
    r"neg[-_ ]?ctrl",
    r"control",
    r"scram",
    r"scrambled",
    r"safe",
    r"aavs1",
]

HIPPO_TARGETS = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]

# D2 branch cell-type definitions
D2_PRE_EC_PATTERNS = [
    r"\bPre-EC\b",
]

D2_LAT_MESO_PATTERNS = [
    r"\bLat-Meso\b",
]

# D3 endpoint definitions
# User requested artEC. If your actual label is artAC, change to r"\bartAC\b".
D3_EC_PATTERNS = [
    r"\bartEC\b",
    # r"\bartAC\b",
]

D3_FIBRO_PATTERNS = [
    r"\bFibro-Mesen\b",
]

# Control modes to export.
# strict_nt = NT-guide cells with zero targeting target genes.
# nt_excl_target = NT-guide cells excluding the target being tested.
# nt_excl_hippo = NT-guide cells excluding Hippo carriers.
CONTROL_MODES = ["strict_nt", "nt_excl_target", "nt_excl_hippo"]

PSEUDOCOUNT = 0.5
MIN_TARGET_ENDPOINT_CELLS = 10
MIN_CONTROL_ENDPOINT_CELLS = 50
```

## Code cell 69

```python
# ============================================================
# Helper functions
# ============================================================

def clean_day(x):
    s = str(x).lower().strip()
    s = s.replace("day", "d").replace(" ", "")

    if s in ["0", "0.0"]:
        return "D0"
    if s in ["1", "1.0"]:
        return "D1"
    if s in ["2", "2.0"]:
        return "D2"
    if s in ["3", "3.0"]:
        return "D3"

    return s.upper()


def matches_any(label, patterns):
    s = str(label)
    return any(re.search(p, s, flags=re.IGNORECASE) for p in patterns)


def looks_control(x):
    s = str(x).lower()
    return any(re.search(p, s, flags=re.IGNORECASE) for p in CONTROL_PATTERNS)


def clean_target_string(x):
    if pd.isna(x):
        return []

    s = str(x).strip()

    if s.lower() in ["", "nan", "none", "na"]:
        return []

    s = (
        s.replace("[", "")
         .replace("]", "")
         .replace("'", "")
         .replace('"', "")
    )

    parts = re.split(r"[;,|/]+", s)
    parts = [p.strip() for p in parts if p.strip() != ""]

    return parts


def canonical_guide_id(x):
    return str(x).strip().replace(" ", "")


def get_guide_assignment_matrix(adata):
    ga = adata.obsm[GUIDE_ASSIGNMENT_KEY]

    if sparse.issparse(ga):
        return ga.tocsr().astype(bool)

    return np.asarray(ga).astype(bool)


def any_guides_per_cell(guide_assignment, guide_indices):
    if len(guide_indices) == 0:
        return np.zeros(guide_assignment.shape[0], dtype=bool)

    if sparse.issparse(guide_assignment):
        return np.asarray(guide_assignment[:, guide_indices].sum(axis=1)).ravel() > 0

    return guide_assignment[:, guide_indices].any(axis=1)


def load_guide_target_from_metadata(adata):
    guide_names = list(map(str, adata.uns[GUIDE_NAMES_KEY]))
    guide_name_to_idx = {g: i for i, g in enumerate(guide_names)}

    meta = pd.read_csv(GUIDE_METADATA_TSV, sep="\t")

    if GUIDE_ID_COL not in meta.columns:
        raise ValueError(f"{GUIDE_ID_COL!r} not found. Columns: {meta.columns.tolist()}")

    if TARGET_COL not in meta.columns:
        raise ValueError(f"{TARGET_COL!r} not found. Columns: {meta.columns.tolist()}")

    meta = meta[[GUIDE_ID_COL, TARGET_COL]].copy()
    meta.columns = ["guide_id", "putative_target_genes"]
    meta["guide_id"] = meta["guide_id"].astype(str)

    # Exact match
    meta["guide_index"] = meta["guide_id"].map(guide_name_to_idx)
    meta["match_method"] = np.where(meta["guide_index"].notna(), "exact_guide_id_match", None)

    # Fallback canonical match
    if meta["guide_index"].isna().any():
        canonical_to_guide = {}
        duplicated = set()

        for g in guide_names:
            cg = canonical_guide_id(g)
            if cg in canonical_to_guide:
                duplicated.add(cg)
            else:
                canonical_to_guide[cg] = g

        for cg in duplicated:
            canonical_to_guide.pop(cg, None)

        for idx in meta.index[meta["guide_index"].isna()]:
            cg = canonical_guide_id(meta.at[idx, "guide_id"])
            if cg in canonical_to_guide:
                matched_guide = canonical_to_guide[cg]
                meta.at[idx, "guide_index"] = guide_name_to_idx[matched_guide]
                meta.at[idx, "match_method"] = "canonical_guide_id_match"

    mapped = meta[meta["guide_index"].notna()].copy()
    mapped["guide_index"] = mapped["guide_index"].astype(int)
    mapped["guide_name"] = [guide_names[i] for i in mapped["guide_index"]]

    records = []

    for _, row in mapped.iterrows():
        guide_index = int(row["guide_index"])
        guide_name = str(row["guide_name"])
        raw_target = row["putative_target_genes"]

        targets = clean_target_string(raw_target)

        if looks_control(guide_name) or looks_control(raw_target):
            targets = [CONTROL_TARGET]

        for target in targets:
            records.append({
                "guide_index": guide_index,
                "guide_name": guide_name,
                "target_name": str(target),
                "raw_target_entry": raw_target,
                "match_method": row["match_method"],
            })

    guide_target = pd.DataFrame(records).drop_duplicates()
    guide_target["target_name"] = guide_target["target_name"].astype(str)

    guide_target.to_csv(
        os.path.join(OUTDIR, "guide_to_target_map_from_metadata.csv"),
        index=False
    )

    print("\nGuide metadata mapping QC")
    print(f"adata guide columns:            {len(guide_names):,}")
    print(f"metadata guide rows:            {len(meta):,}")
    print(f"metadata rows matched to adata: {mapped['guide_id'].nunique():,}")
    print(f"expanded guide-target rows:     {len(guide_target):,}")
    print(f"unique mapped targets:          {guide_target['target_name'].nunique():,}")

    return guide_target


def log2_ratio_of_ratios(a, b, c, d, pseudocount=PSEUDOCOUNT):
    return np.log2(
        ((a + pseudocount) / (b + pseudocount)) /
        ((c + pseudocount) / (d + pseudocount))
    )


def fisher_branch_pval(a, b, c, d):
    try:
        return fisher_exact([[a, b], [c, d]], alternative="two-sided")[1]
    except Exception:
        return np.nan


def safe_fraction(a, b):
    return a / b if b > 0 else np.nan
```

## Code cell 70

```python
# ============================================================
# Build guide/target masks
# ============================================================

guide_target = load_guide_target_from_metadata(adata)

guide_assignment = get_guide_assignment_matrix(adata)

target_to_guide_indices = (
    guide_target
    .groupby("target_name")["guide_index"]
    .apply(lambda x: sorted(set(x.astype(int).tolist())))
    .to_dict()
)

target_masks = {
    target: any_guides_per_cell(guide_assignment, idx)
    for target, idx in target_to_guide_indices.items()
}

control_targets = [
    t for t in target_masks
    if looks_control(t) or str(t) == CONTROL_TARGET
]

if len(control_targets) == 0:
    raise ValueError("No control-like target found in guide metadata.")

control_mask_raw = np.zeros(adata.n_obs, dtype=bool)
for t in control_targets:
    control_mask_raw |= target_masks[t]

targeting_targets = [
    t for t in target_masks
    if t not in control_targets and not looks_control(t)
]

n_targeting_targets = np.zeros(adata.n_obs, dtype=int)
for t in targeting_targets:
    n_targeting_targets += target_masks[t].astype(int)

hippo_mask_raw = np.zeros(adata.n_obs, dtype=bool)
for h in HIPPO_TARGETS:
    matches = [t for t in target_masks if t.upper() == h.upper()]
    for m in matches:
        hippo_mask_raw |= target_masks[m]

qc = pd.DataFrame({
    "metric": [
        "all cells",
        "control guide raw",
        "strict NT: control + no targeting genes",
        "NT excluding Hippo",
        "Hippo carriers",
    ],
    "n_cells": [
        adata.n_obs,
        int(control_mask_raw.sum()),
        int((control_mask_raw & (n_targeting_targets == 0)).sum()),
        int((control_mask_raw & (~hippo_mask_raw)).sum()),
        int(hippo_mask_raw.sum()),
    ]
})
display(qc)
qc.to_csv(os.path.join(OUTDIR, "mask_qc_counts.csv"), index=False)
```

### Text output (stream)

```text

Guide metadata mapping QC
adata guide columns:            1,922
metadata guide rows:            1,926
metadata rows matched to adata: 1,922
expanded guide-target rows:     1,922
unique mapped targets:          299

```

### Text output (text/plain)

```text
                                    metric  n_cells
0                                all cells   121616
1                        control guide raw    21671
2  strict NT: control + no targeting genes      748
3                       NT excluding Hippo    18692
4                           Hippo carriers    10813
```

## Code cell 71

```python
# ============================================================
# Build cell metadata and state labels
# ============================================================

obs = adata.obs.copy()
obs["cell_id"] = adata.obs_names.astype(str)
obs["day_clean"] = obs[DAY_COL].map(clean_day)
obs["cell_type_clean"] = obs[CELLTYPE_COL].astype(str)

obs["D2_branch_state"] = "Other D2"
obs.loc[
    (obs["day_clean"] == "D2") &
    obs["cell_type_clean"].apply(lambda x: matches_any(x, D2_PRE_EC_PATTERNS)),
    "D2_branch_state"
] = "Pre-EC / P15"

obs.loc[
    (obs["day_clean"] == "D2") &
    obs["cell_type_clean"].apply(lambda x: matches_any(x, D2_LAT_MESO_PATTERNS)),
    "D2_branch_state"
] = "Lat-Meso / P17"

obs["D3_endpoint_state"] = "Other D3"
obs.loc[
    (obs["day_clean"] == "D3") &
    obs["cell_type_clean"].apply(lambda x: matches_any(x, D3_EC_PATTERNS)),
    "D3_endpoint_state"
] = "artEC"

obs.loc[
    (obs["day_clean"] == "D3") &
    obs["cell_type_clean"].apply(lambda x: matches_any(x, D3_FIBRO_PATTERNS)),
    "D3_endpoint_state"
] = "Fibro-Mesen"

print("\nD2 cell types:")
print(sorted(obs.loc[obs["day_clean"] == "D2", "cell_type_clean"].unique()))

print("\nD3 cell types:")
print(sorted(obs.loc[obs["day_clean"] == "D3", "cell_type_clean"].unique()))

print("\nD2 branch states:")
print(obs.loc[obs["day_clean"] == "D2", "D2_branch_state"].value_counts())

print("\nD3 endpoint states:")
print(obs.loc[obs["day_clean"] == "D3", "D3_endpoint_state"].value_counts())
```

### Text output (stream)

```text

D2 cell types:
['Fibro-Mesen', 'Lat-Meso', 'Neu-Mesen', 'Neuro-iPSC', 'Others', 'Pre-EC', 'PrimStreak', 'artEC', 'iPSC']

D3 cell types:
['Fibro-Mesen', 'Lat-Meso', 'Neu-Mesen', 'Neuro-iPSC', 'Others', 'Pre-EC', 'PrimStreak', 'artEC', 'iPSC']

D2 branch states:
D2_branch_state
Lat-Meso / P17    16092
Pre-EC / P15       5765
Other D2            123
Name: count, dtype: int64

D3 endpoint states:
D3_endpoint_state
artEC          28444
Fibro-Mesen    12237
Other D3        3421
Name: count, dtype: int64

```

## Code cell 72

```python
# ============================================================
# Helper for choosing controls
# ============================================================

def get_control_mask_for_target(control_mode, target_mask):
    if control_mode == "strict_nt":
        return control_mask_raw & (n_targeting_targets == 0)

    if control_mode == "nt_excl_target":
        return control_mask_raw & (~target_mask)

    if control_mode == "nt_excl_hippo":
        return control_mask_raw & (~hippo_mask_raw)

    raise ValueError(f"Unknown control_mode: {control_mode}")


def count_mask_state(mask, state_vec, state_name):
    return int(np.sum(mask & (state_vec == state_name)))
```

## Code cell 73

```python
# ============================================================
# Export D3 endpoint bias: Fibro-Mesen vs artEC
# ============================================================

d3_state_vec = obs["D3_endpoint_state"].astype(str).values
d3_mask = obs["day_clean"].values == "D3"

d3_records = []

for control_mode in CONTROL_MODES:
    for target, target_mask in target_masks.items():
        if target in control_targets or looks_control(target):
            continue

        target_use = target_mask & d3_mask
        control_use = get_control_mask_for_target(control_mode, target_mask) & d3_mask
        control_use = control_use & (~target_mask)

        t_fibro = count_mask_state(target_use, d3_state_vec, "Fibro-Mesen")
        t_ec = count_mask_state(target_use, d3_state_vec, "artEC")

        c_fibro = count_mask_state(control_use, d3_state_vec, "Fibro-Mesen")
        c_ec = count_mask_state(control_use, d3_state_vec, "artEC")

        t_endpoint_total = t_fibro + t_ec
        c_endpoint_total = c_fibro + c_ec

        if (
            t_endpoint_total < MIN_TARGET_ENDPOINT_CELLS
            or c_endpoint_total < MIN_CONTROL_ENDPOINT_CELLS
        ):
            effect = np.nan
            pval = np.nan
        else:
            effect = log2_ratio_of_ratios(t_fibro, t_ec, c_fibro, c_ec)
            pval = fisher_branch_pval(t_fibro, t_ec, c_fibro, c_ec)

        t_fibro_frac = safe_fraction(t_fibro, t_endpoint_total)
        c_fibro_frac = safe_fraction(c_fibro, c_endpoint_total)

        d3_records.append({
            "target": target,
            "control_mode": control_mode,

            "D3_endpoint_bias_FibroMesen_vs_artEC_log2OR": effect,
            "D3_endpoint_pval": pval,

            "target_D3_FibroMesen_count": t_fibro,
            "target_D3_artEC_count": t_ec,
            "control_D3_FibroMesen_count": c_fibro,
            "control_D3_artEC_count": c_ec,

            "target_D3_endpoint_cells": t_endpoint_total,
            "control_D3_endpoint_cells": c_endpoint_total,

            "target_FibroMesen_fraction": t_fibro_frac,
            "control_FibroMesen_fraction": c_fibro_frac,
            "D3_FibroMesen_fraction_shift": t_fibro_frac - c_fibro_frac,
            "D3_FibroMesen_fraction_shift_abs": abs(t_fibro_frac - c_fibro_frac),
        })

d3_endpoint = pd.DataFrame(d3_records)

d3_endpoint["D3_endpoint_fdr"] = np.nan
for mode, sub in d3_endpoint.groupby("control_mode"):
    valid = sub["D3_endpoint_pval"].notna()
    if valid.sum() > 0:
        d3_endpoint.loc[sub.index[valid], "D3_endpoint_fdr"] = multipletests(
            sub.loc[valid, "D3_endpoint_pval"],
            method="fdr_bh"
        )[1]

d3_endpoint["minus_log10_D3_fdr"] = -np.log10(
    d3_endpoint["D3_endpoint_fdr"].clip(lower=1e-300)
)

d3_endpoint_out = os.path.join(OUTDIR, "target_D3_endpoint_bias_FibroMesen_vs_artEC.csv")
d3_endpoint.to_csv(d3_endpoint_out, index=False)

print("Saved:", d3_endpoint_out)
display(d3_endpoint.head())
```

### Text output (stream)

```text
Saved: scatter_plot_celllevel_endpoint_exports/target_D3_endpoint_bias_FibroMesen_vs_artEC.csv

```

### Text output (text/plain)

```text
   target control_mode  D3_endpoint_bias_FibroMesen_vs_artEC_log2OR  \
0   ACAA1    strict_nt                                     2.414119   
1  ACTL7B    strict_nt                                     2.783547   
2  ACVR1B    strict_nt                                     3.714540   
3  ACVR2A    strict_nt                                     3.423246   
4  ACVR2B    strict_nt                                     3.585079   

   D3_endpoint_pval  target_D3_FibroMesen_count  target_D3_artEC_count  \
0      8.723814e-23                         218                    361   
1      4.001897e-31                         256                    328   
2      3.027991e-59                         362                    243   
3      2.804027e-42                         213                    175   
4      2.498553e-53                         320                    235   

   control_D3_FibroMesen_count  control_D3_artEC_count  \
0                           38                     339   
1                           38                     339   
2                           38                     339   
3                           38                     339   
4                           38                     339   

   target_D3_endpoint_cells  control_D3_endpoint_cells  \
0                       579                        377   
1                       584                        377   
2                       605                        377   
3                       388                        377   
4                       555                        377   

   target_FibroMesen_fraction  control_FibroMesen_fraction  \
0                    0.376511                     0.100796   
1                    0.438356                     0.100796   
2                    0.598347                     0.100796   
3                    0.548969                     0.100796   
4                    0.576577                     0.100796   

   D3_FibroMesen_fraction_shift  D3_FibroMesen_fraction_shift_abs  \
0                      0.275715                          0.275715   
1                      0.337560                          0.337560   
2                      0.497551                          0.497551   
3                      0.448173                          0.448173   
4                      0.475781                          0.475781   

   D3_endpoint_fdr  minus_log10_D3_fdr  
0     9.385186e-23           22.027557  
1     6.276659e-31           30.202271  
2     4.101551e-58           57.387052  
3     1.099474e-41           40.958815  
4     2.326777e-52           51.633245  
```

## Code cell 74

```python
# ============================================================
# Export D3 counts and fractions by target/state
# ============================================================

d3_count_records = []

for control_mode in CONTROL_MODES:
    for target, target_mask in target_masks.items():
        if target in control_targets or looks_control(target):
            continue

        target_use = target_mask & d3_mask
        control_use = get_control_mask_for_target(control_mode, target_mask) & d3_mask
        control_use = control_use & (~target_mask)

        for group_name, mask in [
            ("target", target_use),
            ("control", control_use),
        ]:
            for state_name in ["artEC", "Fibro-Mesen", "Other D3"]:
                n = count_mask_state(mask, d3_state_vec, state_name)
                total = int(mask.sum())
                d3_count_records.append({
                    "target": target,
                    "control_mode": control_mode,
                    "group": group_name,
                    "state": state_name,
                    "n_cells": n,
                    "total_D3_cells": total,
                    "fraction_among_all_D3": safe_fraction(n, total),
                })

d3_counts = pd.DataFrame(d3_count_records)
d3_counts_out = os.path.join(OUTDIR, "target_D3_state_counts_and_fractions.csv")
d3_counts.to_csv(d3_counts_out, index=False)

print("Saved:", d3_counts_out)
display(d3_counts.head())
```

### Text output (stream)

```text
Saved: scatter_plot_celllevel_endpoint_exports/target_D3_state_counts_and_fractions.csv

```

### Text output (text/plain)

```text
  target control_mode    group        state  n_cells  total_D3_cells  \
0  ACAA1    strict_nt   target        artEC      361             634   
1  ACAA1    strict_nt   target  Fibro-Mesen      218             634   
2  ACAA1    strict_nt   target     Other D3       55             634   
3  ACAA1    strict_nt  control        artEC      339             378   
4  ACAA1    strict_nt  control  Fibro-Mesen       38             378   

   fraction_among_all_D3  
0               0.569401  
1               0.343849  
2               0.086751  
3               0.896825  
4               0.100529  
```

## Code cell 75

```python
# ============================================================
# Export D2 cell-type branch bias too, optional for other scatter rows
# Lat-Meso/P17 vs Pre-EC/P15
# ============================================================

d2_state_vec = obs["D2_branch_state"].astype(str).values
d2_mask = obs["day_clean"].values == "D2"

d2_records = []

for control_mode in CONTROL_MODES:
    for target, target_mask in target_masks.items():
        if target in control_targets or looks_control(target):
            continue

        target_use = target_mask & d2_mask
        control_use = get_control_mask_for_target(control_mode, target_mask) & d2_mask
        control_use = control_use & (~target_mask)

        t_lat = count_mask_state(target_use, d2_state_vec, "Lat-Meso / P17")
        t_pre = count_mask_state(target_use, d2_state_vec, "Pre-EC / P15")

        c_lat = count_mask_state(control_use, d2_state_vec, "Lat-Meso / P17")
        c_pre = count_mask_state(control_use, d2_state_vec, "Pre-EC / P15")

        t_total = t_lat + t_pre
        c_total = c_lat + c_pre

        if t_total < MIN_TARGET_ENDPOINT_CELLS or c_total < MIN_CONTROL_ENDPOINT_CELLS:
            effect = np.nan
            pval = np.nan
        else:
            effect = log2_ratio_of_ratios(t_lat, t_pre, c_lat, c_pre)
            pval = fisher_branch_pval(t_lat, t_pre, c_lat, c_pre)

        t_lat_frac = safe_fraction(t_lat, t_total)
        c_lat_frac = safe_fraction(c_lat, c_total)

        d2_records.append({
            "target": target,
            "control_mode": control_mode,

            "D2_branch_bias_LatMeso_vs_PreEC_log2OR": effect,
            "D2_branch_pval": pval,

            "target_D2_LatMeso_count": t_lat,
            "target_D2_PreEC_count": t_pre,
            "control_D2_LatMeso_count": c_lat,
            "control_D2_PreEC_count": c_pre,

            "target_D2_branch_cells": t_total,
            "control_D2_branch_cells": c_total,

            "target_LatMeso_fraction": t_lat_frac,
            "control_LatMeso_fraction": c_lat_frac,
            "D2_LatMeso_fraction_shift": t_lat_frac - c_lat_frac,
            "D2_LatMeso_fraction_shift_abs": abs(t_lat_frac - c_lat_frac),
        })

d2_branch = pd.DataFrame(d2_records)

d2_branch["D2_branch_fdr"] = np.nan
for mode, sub in d2_branch.groupby("control_mode"):
    valid = sub["D2_branch_pval"].notna()
    if valid.sum() > 0:
        d2_branch.loc[sub.index[valid], "D2_branch_fdr"] = multipletests(
            sub.loc[valid, "D2_branch_pval"],
            method="fdr_bh"
        )[1]

d2_branch["minus_log10_D2_branch_fdr"] = -np.log10(
    d2_branch["D2_branch_fdr"].clip(lower=1e-300)
)

d2_branch_out = os.path.join(OUTDIR, "target_D2_branch_celltype_bias_LatMeso_vs_PreEC.csv")
d2_branch.to_csv(d2_branch_out, index=False)

print("Saved:", d2_branch_out)
display(d2_branch.head())
```

### Text output (stream)

```text
Saved: scatter_plot_celllevel_endpoint_exports/target_D2_branch_celltype_bias_LatMeso_vs_PreEC.csv

```

### Text output (text/plain)

```text
   target control_mode  D2_branch_bias_LatMeso_vs_PreEC_log2OR  \
0   ACAA1    strict_nt                                1.425091   
1  ACTL7B    strict_nt                                1.839837   
2  ACVR1B    strict_nt                                2.496546   
3  ACVR2A    strict_nt                                3.120178   
4  ACVR2B    strict_nt                                2.714594   

   D2_branch_pval  target_D2_LatMeso_count  target_D2_PreEC_count  \
0    5.304115e-05                      307                     94   
1    2.953691e-07                      327                     75   
2    3.478162e-11                      338                     49   
3    3.309083e-15                      363                     34   
4    3.609166e-13                      433                     54   

   control_D2_LatMeso_count  control_D2_PreEC_count  target_D2_branch_cells  \
0                        51                      42                     401   
1                        51                      42                     402   
2                        51                      42                     387   
3                        51                      42                     397   
4                        51                      42                     487   

   control_D2_branch_cells  target_LatMeso_fraction  control_LatMeso_fraction  \
0                       93                 0.765586                  0.548387   
1                       93                 0.813433                  0.548387   
2                       93                 0.873385                  0.548387   
3                       93                 0.914358                  0.548387   
4                       93                 0.889117                  0.548387   

   D2_LatMeso_fraction_shift  D2_LatMeso_fraction_shift_abs  D2_branch_fdr  \
0                   0.217199                       0.217199   5.625004e-05   
1                   0.265046                       0.265046   5.715584e-07   
2                   0.324998                       0.324998   4.506488e-10   
3                   0.365971                       0.365971   1.972214e-13   
4                   0.340730                       0.340730   9.777559e-12   

   minus_log10_D2_branch_fdr  
0                   4.249877  
1                   6.242939  
2                   9.346162  
3                  12.705046  
4                  11.009770  
```

## Code cell 76

```python
# ============================================================
# Export D2 counts and fractions by target/state
# ============================================================

d2_count_records = []

for control_mode in CONTROL_MODES:
    for target, target_mask in target_masks.items():
        if target in control_targets or looks_control(target):
            continue

        target_use = target_mask & d2_mask
        control_use = get_control_mask_for_target(control_mode, target_mask) & d2_mask
        control_use = control_use & (~target_mask)

        for group_name, mask in [
            ("target", target_use),
            ("control", control_use),
        ]:
            for state_name in ["Pre-EC / P15", "Lat-Meso / P17", "Other D2"]:
                n = count_mask_state(mask, d2_state_vec, state_name)
                total = int(mask.sum())
                d2_count_records.append({
                    "target": target,
                    "control_mode": control_mode,
                    "group": group_name,
                    "state": state_name,
                    "n_cells": n,
                    "total_D2_cells": total,
                    "fraction_among_all_D2": safe_fraction(n, total),
                })

d2_counts = pd.DataFrame(d2_count_records)
d2_counts_out = os.path.join(OUTDIR, "target_D2_state_counts_and_fractions.csv")
d2_counts.to_csv(d2_counts_out, index=False)

print("Saved:", d2_counts_out)
display(d2_counts.head())
```

### Text output (stream)

```text
Saved: scatter_plot_celllevel_endpoint_exports/target_D2_state_counts_and_fractions.csv

```

### Text output (text/plain)

```text
  target control_mode    group           state  n_cells  total_D2_cells  \
0  ACAA1    strict_nt   target    Pre-EC / P15       94             405   
1  ACAA1    strict_nt   target  Lat-Meso / P17      307             405   
2  ACAA1    strict_nt   target        Other D2        4             405   
3  ACAA1    strict_nt  control    Pre-EC / P15       42              93   
4  ACAA1    strict_nt  control  Lat-Meso / P17       51              93   

   fraction_among_all_D2  
0               0.232099  
1               0.758025  
2               0.009877  
3               0.451613  
4               0.548387  
```

## Code cell 77

```python
# ============================================================
# QC target counts
# ============================================================

target_qc_records = []

for target, mask in target_masks.items():
    if target in control_targets or looks_control(target):
        continue

    target_qc_records.append({
        "target": target,
        "n_total_carrier_cells": int(mask.sum()),
        "n_D2_carrier_cells": int(np.sum(mask & d2_mask)),
        "n_D3_carrier_cells": int(np.sum(mask & d3_mask)),
        "n_D2_PreEC": int(np.sum(mask & d2_mask & (d2_state_vec == "Pre-EC / P15"))),
        "n_D2_LatMeso": int(np.sum(mask & d2_mask & (d2_state_vec == "Lat-Meso / P17"))),
        "n_D3_artEC": int(np.sum(mask & d3_mask & (d3_state_vec == "artEC"))),
        "n_D3_FibroMesen": int(np.sum(mask & d3_mask & (d3_state_vec == "Fibro-Mesen"))),
    })

target_qc = pd.DataFrame(target_qc_records)

target_qc_out = os.path.join(OUTDIR, "target_qc_counts.csv")
target_qc.to_csv(target_qc_out, index=False)

print("Saved:", target_qc_out)
display(target_qc.head())
```

### Text output (stream)

```text
Saved: scatter_plot_celllevel_endpoint_exports/target_qc_counts.csv

```

### Text output (text/plain)

```text
   target  n_total_carrier_cells  n_D2_carrier_cells  n_D3_carrier_cells  \
0   ACAA1                   1966                 405                 634   
1  ACTL7B                   1944                 407                 638   
2  ACVR1B                   1972                 390                 639   
3  ACVR2A                   2012                 400                 585   
4  ACVR2B                   2306                 489                 621   

   n_D2_PreEC  n_D2_LatMeso  n_D3_artEC  n_D3_FibroMesen  
0          94           307         361              218  
1          75           327         328              256  
2          49           338         243              362  
3          34           363         175              213  
4          54           433         235              320  
```

## Markdown cell 78

#### Export %EC for all targets across 2 replicates

## Code cell 79

```python
# ============================================================
# Export D3 %EC by target KD and biological replicate
# ============================================================
# Requires:
#   - adata_int
#   - adata_int.obs["biological_sample"]
#   - adata_int.obs["cell_type"]
#   - adata_int.uns["guide_names"]
#   - adata_int.obsm["guide_assignment"]
#   - 300genes_guide_metadata_v43.tsv
#
# Output:
#   d3_EC_percent_by_target_biorep_wide.csv
#   d3_EC_percent_by_target_biorep_long.csv
#   d3_EC_percent_by_control_biorep.csv
#   guide_to_target_mapping_used_for_d3_biorep_export.csv
#   d3_biorep_export_qc_summary.csv

import os
import re
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse

# ----------------------------
# Settings
# ----------------------------
OUTDIR = Path("d3_EC_percent_by_target_biorep")
OUTDIR.mkdir(parents=True, exist_ok=True)

GUIDE_METADATA_TSV = "300genes_guide_metadata_v43.tsv"

BIO_SAMPLE_COL = "biological_sample"
CELLTYPE_COL = "cell_type"

GUIDE_NAMES_KEY = "guide_names"
GUIDE_ASSIGNMENT_KEY = "guide_assignment"

D3_REPLICATES = ["d3_rep1", "d3_rep2"]

# Default EC definition.
# This currently treats artEC as EC.
# To include ECM-EC as EC too, change to:
# EC_PATTERNS = [r"\bartEC\b", r"\bECM-EC\b"]
EC_PATTERNS = [r"\bartEC\b"]

CONTROL_TARGET_NAME = "non-targeting_control"

CONTROL_PATTERNS = [
    r"non[-_ ]?target",
    r"non[-_ ]?targeting",
    r"\bNTC\b",
    r"negative",
    r"control",
    r"scram",
    r"scrambled",
    r"safe",
    r"aavs1",
]

TARGET_COL_CANDIDATES = [
    "intended_target_name",
    "putative_target_genes",
    "target_name",
    "target",
    "gene",
]

# ----------------------------
# Helper functions
# ----------------------------
def matches_any(x, patterns):
    s = str(x)
    return any(re.search(p, s, flags=re.IGNORECASE) for p in patterns)


def looks_control(x):
    return matches_any(x, CONTROL_PATTERNS)


def clean_target_entry(x):
    """Return one or more target names from metadata target field."""
    if pd.isna(x):
        return []
    s = str(x).strip()
    if s.lower() in ["", "nan", "none", "na"]:
        return []

    s = (
        s.replace("[", "")
         .replace("]", "")
         .replace("'", "")
         .replace('"', "")
    )

    parts = re.split(r"[;,|/]+", s)
    parts = [p.strip() for p in parts if p.strip()]
    return parts


def classify_ec(cell_type):
    return matches_any(cell_type, EC_PATTERNS)


def canonical(x):
    return str(x).strip().replace(" ", "")


def get_guide_assignment_matrix(adata):
    ga = adata.obsm[GUIDE_ASSIGNMENT_KEY]
    if sparse.issparse(ga):
        return ga.tocsr().astype(bool)
    return np.asarray(ga).astype(bool)


def any_guides_per_cell(guide_assignment, guide_indices):
    if len(guide_indices) == 0:
        return np.zeros(guide_assignment.shape[0], dtype=bool)

    if sparse.issparse(guide_assignment):
        return np.asarray(guide_assignment[:, guide_indices].sum(axis=1)).ravel() > 0

    return guide_assignment[:, guide_indices].any(axis=1)


def infer_target_col(meta):
    for c in TARGET_COL_CANDIDATES:
        if c in meta.columns:
            return c
    raise ValueError(
        "Could not find target column in guide metadata. "
        f"Available columns: {meta.columns.tolist()}"
    )


# ----------------------------
# Basic checks
# ----------------------------
for col in [BIO_SAMPLE_COL, CELLTYPE_COL]:
    if col not in adata_int.obs.columns:
        raise ValueError(f"adata_int.obs is missing required column: {col}")

if GUIDE_NAMES_KEY not in adata_int.uns:
    raise ValueError(f"adata_int.uns is missing: {GUIDE_NAMES_KEY}")

if GUIDE_ASSIGNMENT_KEY not in adata_int.obsm:
    raise ValueError(f"adata_int.obsm is missing: {GUIDE_ASSIGNMENT_KEY}")

if not Path(GUIDE_METADATA_TSV).exists():
    raise FileNotFoundError(
        f"Could not find {GUIDE_METADATA_TSV}. "
        "Put it in the notebook working directory or update GUIDE_METADATA_TSV."
    )

# ----------------------------
# Load guide metadata and map guide -> target
# ----------------------------
guide_names = list(map(str, adata_int.uns[GUIDE_NAMES_KEY]))
guide_name_to_idx = {g: i for i, g in enumerate(guide_names)}

meta = pd.read_csv(GUIDE_METADATA_TSV, sep="\t")
target_col = infer_target_col(meta)

if "guide_id" not in meta.columns:
    raise ValueError(f"guide_id column missing from {GUIDE_METADATA_TSV}")

meta = meta.copy()
meta["guide_id"] = meta["guide_id"].astype(str)
meta["raw_target_entry"] = meta[target_col]

# exact match
meta["guide_index"] = meta["guide_id"].map(guide_name_to_idx)
meta["match_method"] = np.where(meta["guide_index"].notna(), "exact", None)

# fallback canonical match
if meta["guide_index"].isna().any():
    canonical_to_guide = {}
    duplicated = set()

    for g in guide_names:
        cg = canonical(g)
        if cg in canonical_to_guide:
            duplicated.add(cg)
        else:
            canonical_to_guide[cg] = g

    for cg in duplicated:
        canonical_to_guide.pop(cg, None)

    for idx in meta.index[meta["guide_index"].isna()]:
        cg = canonical(meta.at[idx, "guide_id"])
        if cg in canonical_to_guide:
            matched_guide = canonical_to_guide[cg]
            meta.at[idx, "guide_index"] = guide_name_to_idx[matched_guide]
            meta.at[idx, "match_method"] = "canonical"

mapped_meta = meta[meta["guide_index"].notna()].copy()
mapped_meta["guide_index"] = mapped_meta["guide_index"].astype(int)
mapped_meta["guide_name"] = mapped_meta["guide_index"].map(lambda i: guide_names[i])

records = []

for _, row in mapped_meta.iterrows():
    guide_id = row["guide_id"]
    guide_name = row["guide_name"]
    guide_index = int(row["guide_index"])
    raw_target = row["raw_target_entry"]

    if looks_control(guide_id) or looks_control(guide_name) or looks_control(raw_target):
        targets = [CONTROL_TARGET_NAME]
    else:
        targets = clean_target_entry(raw_target)

    for target in targets:
        records.append({
            "guide_id": guide_id,
            "guide_name": guide_name,
            "guide_index": guide_index,
            "target": str(target),
            "raw_target_entry": raw_target,
            "match_method": row["match_method"],
        })

guide_target = pd.DataFrame(records).drop_duplicates()
guide_target.to_csv(
    OUTDIR / "guide_to_target_mapping_used_for_d3_biorep_export.csv",
    index=False
)

# ----------------------------
# Build target -> guide indices
# ----------------------------
target_to_guide_indices = (
    guide_target
    .groupby("target")["guide_index"]
    .apply(lambda x: sorted(set(x.astype(int).tolist())))
    .to_dict()
)

control_targets = [
    t for t in target_to_guide_indices
    if t == CONTROL_TARGET_NAME or looks_control(t)
]

targeting_targets = [
    t for t in target_to_guide_indices
    if t not in control_targets and not looks_control(t)
]

guide_assignment = get_guide_assignment_matrix(adata_int)

if guide_assignment.shape[1] != len(guide_names):
    raise ValueError(
        f"guide_assignment has {guide_assignment.shape[1]} columns, "
        f"but adata_int.uns['{GUIDE_NAMES_KEY}'] has {len(guide_names)} guide names."
    )

# Precompute target masks.
target_masks = {
    t: any_guides_per_cell(guide_assignment, idxs)
    for t, idxs in target_to_guide_indices.items()
}

control_mask = np.zeros(adata_int.n_obs, dtype=bool)
for t in control_targets:
    control_mask |= target_masks[t]

any_targeting_mask = np.zeros(adata_int.n_obs, dtype=bool)
for t in targeting_targets:
    any_targeting_mask |= target_masks[t]

strict_nt_mask = control_mask & (~any_targeting_mask)

# ----------------------------
# Cell metadata
# ----------------------------
obs = adata_int.obs.copy()
obs["cell_id"] = adata_int.obs_names.astype(str)
obs["biological_sample_clean"] = obs[BIO_SAMPLE_COL].astype(str)
obs["is_D3_replicate"] = obs["biological_sample_clean"].isin(D3_REPLICATES)
obs["is_EC"] = obs[CELLTYPE_COL].astype(str).map(classify_ec)

print("D3 biological_sample counts:")
display(obs.loc[obs["is_D3_replicate"], "biological_sample_clean"].value_counts())

print("D3 EC/non-EC counts:")
display(
    obs.loc[obs["is_D3_replicate"]]
    .groupby(["biological_sample_clean", "is_EC"])
    .size()
    .rename("n_cells")
    .reset_index()
)

# ----------------------------
# Export target-by-replicate EC percentage
# ----------------------------
rows = []

for target in sorted(targeting_targets):
    carrier_mask = target_masks[target]
    n_guides_for_target = len(target_to_guide_indices[target])

    for rep in D3_REPLICATES:
        rep_mask = obs["biological_sample_clean"].values == rep
        use = carrier_mask & rep_mask

        n_target_cells = int(use.sum())
        n_EC_cells = int((use & obs["is_EC"].values).sum())
        n_nonEC_cells = n_target_cells - n_EC_cells

        pct_EC = np.nan
        if n_target_cells > 0:
            pct_EC = 100 * n_EC_cells / n_target_cells

        rows.append({
            "target": target,
            "replicate": rep,
            "n_guides_mapped_to_target": n_guides_for_target,
            "n_target_carrier_cells": n_target_cells,
            "n_EC_target_carrier_cells": n_EC_cells,
            "n_nonEC_target_carrier_cells": n_nonEC_cells,
            "pct_EC_target_carrier_cells": pct_EC,
        })

long_df = pd.DataFrame(rows)

# Wide table requested by user.
wide_df = long_df.pivot(index="target", columns="replicate")
wide_df.columns = [
    f"{rep}_{metric}"
    for metric, rep in wide_df.columns
]
wide_df = wide_df.reset_index()

# Better column ordering.
ordered_cols = ["target"]

for rep in D3_REPLICATES:
    ordered_cols += [
        f"{rep}_pct_EC_target_carrier_cells",
        f"{rep}_n_target_carrier_cells",
        f"{rep}_n_EC_target_carrier_cells",
        f"{rep}_n_nonEC_target_carrier_cells",
        f"{rep}_n_guides_mapped_to_target",
    ]

wide_df = wide_df[[c for c in ordered_cols if c in wide_df.columns]]

long_out = OUTDIR / "d3_EC_percent_by_target_biorep_long.csv"
wide_out = OUTDIR / "d3_EC_percent_by_target_biorep_wide.csv"

long_df.to_csv(long_out, index=False)
wide_df.to_csv(wide_out, index=False)

print("Saved long table:", long_out)
print("Saved wide table:", wide_out)

display(wide_df.head())

# ----------------------------
# Optional: strict NT control EC percentage by replicate
# ----------------------------
control_rows = []

for rep in D3_REPLICATES:
    rep_mask = obs["biological_sample_clean"].values == rep
    use = strict_nt_mask & rep_mask

    n_cells = int(use.sum())
    n_EC = int((use & obs["is_EC"].values).sum())
    n_nonEC = n_cells - n_EC
    pct_EC = np.nan if n_cells == 0 else 100 * n_EC / n_cells

    control_rows.append({
        "control_group": "Strict NT",
        "replicate": rep,
        "n_control_cells": n_cells,
        "n_EC_control_cells": n_EC,
        "n_nonEC_control_cells": n_nonEC,
        "pct_EC_control_cells": pct_EC,
    })

control_df = pd.DataFrame(control_rows)
control_out = OUTDIR / "d3_EC_percent_by_control_biorep.csv"
control_df.to_csv(control_out, index=False)

print("Saved control table:", control_out)
display(control_df)

# ----------------------------
# QC summary
# ----------------------------
qc_rows = []

qc_rows.append({
    "metric": "n_adata_cells",
    "value": adata_int.n_obs,
})

qc_rows.append({
    "metric": "n_guide_names",
    "value": len(guide_names),
})

qc_rows.append({
    "metric": "n_metadata_guides",
    "value": len(meta),
})

qc_rows.append({
    "metric": "n_metadata_guides_mapped_to_adata_guides",
    "value": len(mapped_meta),
})

qc_rows.append({
    "metric": "n_targeting_targets_exported",
    "value": len(targeting_targets),
})

qc_rows.append({
    "metric": "n_control_targets_detected",
    "value": len(control_targets),
})

qc_rows.append({
    "metric": "n_d3_rep1_cells",
    "value": int((obs["biological_sample_clean"].values == "d3_rep1").sum()),
})

qc_rows.append({
    "metric": "n_d3_rep2_cells",
    "value": int((obs["biological_sample_clean"].values == "d3_rep2").sum()),
})

qc_df = pd.DataFrame(qc_rows)
qc_out = OUTDIR / "d3_biorep_export_qc_summary.csv"
qc_df.to_csv(qc_out, index=False)

print("Saved QC summary:", qc_out)
display(qc_df)
```

### Text output (stream)

```text
D3 biological_sample counts:

```

### Text output (text/plain)

```text
biological_sample_clean
d3_rep2    22725
d3_rep1    21377
Name: count, dtype: int64
```

### Text output (stream)

```text
D3 EC/non-EC counts:

```

### Text output (text/plain)

```text
  biological_sample_clean  is_EC  n_cells
0                 d3_rep1  False     7460
1                 d3_rep1   True    13917
2                 d3_rep2  False     8198
3                 d3_rep2   True    14527
```

### Text output (stream)

```text
Saved long table: d3_EC_percent_by_target_biorep/d3_EC_percent_by_target_biorep_long.csv
Saved wide table: d3_EC_percent_by_target_biorep/d3_EC_percent_by_target_biorep_wide.csv

```

### Text output (text/plain)

```text
   target  d3_rep1_pct_EC_target_carrier_cells  \
0   ACAA1                            58.360656   
1  ACTL7B                            52.786885   
2  ACVR1B                            40.070922   
3  ACVR2A                            29.889299   
4  ACVR2B                            41.516245   

   d3_rep1_n_target_carrier_cells  d3_rep1_n_EC_target_carrier_cells  \
0                             305                                178   
1                             305                                161   
2                             282                                113   
3                             271                                 81   
4                             277                                115   

   d3_rep1_n_nonEC_target_carrier_cells  d3_rep1_n_guides_mapped_to_target  \
0                                   127                                  6   
1                                   144                                  6   
2                                   169                                  6   
3                                   190                                  6   
4                                   162                                  6   

   d3_rep2_pct_EC_target_carrier_cells  d3_rep2_n_target_carrier_cells  \
0                            55.623100                             329   
1                            50.150150                             333   
2                            36.414566                             357   
3                            29.936306                             314   
4                            34.883721                             344   

   d3_rep2_n_EC_target_carrier_cells  d3_rep2_n_nonEC_target_carrier_cells  \
0                                183                                   146   
1                                167                                   166   
2                                130                                   227   
3                                 94                                   220   
4                                120                                   224   

   d3_rep2_n_guides_mapped_to_target  
0                                  6  
1                                  6  
2                                  6  
3                                  6  
4                                  6  
```

### Text output (stream)

```text
Saved control table: d3_EC_percent_by_target_biorep/d3_EC_percent_by_control_biorep.csv

```

### Text output (text/plain)

```text
  control_group replicate  n_control_cells  n_EC_control_cells  \
0     Strict NT   d3_rep1              185                 165   
1     Strict NT   d3_rep2              193                 174   

   n_nonEC_control_cells  pct_EC_control_cells  
0                     20             89.189189  
1                     19             90.155440  
```

### Text output (stream)

```text
Saved QC summary: d3_EC_percent_by_target_biorep/d3_biorep_export_qc_summary.csv

```

### Text output (text/plain)

```text
                                     metric   value
0                             n_adata_cells  121616
1                             n_guide_names    1922
2                         n_metadata_guides    1926
3  n_metadata_guides_mapped_to_adata_guides    1922
4              n_targeting_targets_exported     298
5                n_control_targets_detected       1
6                           n_d3_rep1_cells   21377
7                           n_d3_rep2_cells   22725
```

## Markdown cell 80

#### %EC in individual guides

## Code cell 81

```python
# ============================================================
# Export D3 %EC by individual guide, plus individual NT guides
# ============================================================
# Requires:
#   - adata_int
#   - adata_int.obs["cell_type"]
#   - adata_int.obs["biological_sample"] or adata_int.obs["day"]
#   - adata_int.uns["guide_names"]
#   - adata_int.obsm["guide_assignment"]
#   - 300genes_guide_metadata_v43.tsv
#
# Main output:
#   d3_EC_percent_by_guide_with_individual_NT_excl_Hippo.csv
#
# Rows:
#   - one row per targeting guide
#   - one row per individual non-targeting guide, excluding Hippo carriers
#   - one aggregate row: NT excluding Hippo
#
# EC definition:
#   default EC = artEC only

import os
import re
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse

# ----------------------------
# Settings
# ----------------------------
OUTDIR = Path("d3_EC_percent_by_guide")
OUTDIR.mkdir(parents=True, exist_ok=True)

GUIDE_METADATA_TSV = "300genes_guide_metadata_v43.tsv"

CELLTYPE_COL = "cell_type"
BIO_SAMPLE_COL = "biological_sample"
DAY_COL = "day"

GUIDE_NAMES_KEY = "guide_names"
GUIDE_ASSIGNMENT_KEY = "guide_assignment"

# EC definition.
# Current: EC = artEC only.
# To include ECM-EC too:
# EC_PATTERNS = [r"\bartEC\b", r"\bECM-EC\b"]
EC_PATTERNS = [r"\bartEC\b"]

HIPPO_TARGETS = ["NF2", "LATS2", "SAV1", "TAOK1", "AMOTL2"]
HIPPO_TARGETS_UPPER = {x.upper() for x in HIPPO_TARGETS}

CONTROL_TARGET_NAME = "non-targeting_control"

CONTROL_PATTERNS = [
    r"non[-_ ]?target",
    r"non[-_ ]?targeting",
    r"\bNTC\b",
    r"negative",
    r"control",
    r"scram",
    r"scrambled",
    r"safe",
    r"aavs1",
]

TARGET_COL_CANDIDATES = [
    "intended_target_name",
    "putative_target_genes",
    "target_name",
    "target",
    "gene",
]

# Since strict NT has too few cells, keep this off by default.
EXPORT_STRICT_NT_AGGREGATE = False

EXPORT_CELL_MEMBERSHIP = False

# ----------------------------
# Helper functions
# ----------------------------
def matches_any(x, patterns):
    s = str(x)
    return any(re.search(p, s, flags=re.IGNORECASE) for p in patterns)


def looks_control(x):
    return matches_any(x, CONTROL_PATTERNS)


def clean_target_entry(x):
    if pd.isna(x):
        return []
    s = str(x).strip()
    if s.lower() in ["", "nan", "none", "na"]:
        return []

    s = (
        s.replace("[", "")
         .replace("]", "")
         .replace("'", "")
         .replace('"', "")
    )

    parts = re.split(r"[;,|/]+", s)
    parts = [p.strip() for p in parts if p.strip()]
    return parts


def classify_ec(cell_type):
    return matches_any(cell_type, EC_PATTERNS)


def canonical(x):
    return str(x).strip().replace(" ", "")


def clean_day_from_biosample_or_day(obs):
    if BIO_SAMPLE_COL in obs.columns:
        s = obs[BIO_SAMPLE_COL].astype(str).str.lower()
        return np.where(
            s.str.startswith("d3"),
            "D3",
            np.where(
                s.str.startswith("d2"),
                "D2",
                np.where(
                    s.str.startswith("d1"),
                    "D1",
                    np.where(s.str.startswith("d0"), "D0", s.str.upper())
                )
            )
        )

    if DAY_COL in obs.columns:
        s = (
            obs[DAY_COL]
            .astype(str)
            .str.lower()
            .str.replace("day", "d", regex=False)
            .str.replace(" ", "", regex=False)
        )

        def convert(x):
            if x in ["0", "0.0"]:
                return "D0"
            if x in ["1", "1.0"]:
                return "D1"
            if x in ["2", "2.0"]:
                return "D2"
            if x in ["3", "3.0"]:
                return "D3"
            return x.upper()

        return s.map(convert).values

    raise ValueError(
        f"Could not determine day. Need either obs['{BIO_SAMPLE_COL}'] or obs['{DAY_COL}']."
    )


def get_guide_assignment_matrix(adata):
    ga = adata.obsm[GUIDE_ASSIGNMENT_KEY]
    if sparse.issparse(ga):
        return ga.tocsr().astype(bool)
    return np.asarray(ga).astype(bool)


def any_guides_per_cell(guide_assignment, guide_indices):
    if len(guide_indices) == 0:
        return np.zeros(guide_assignment.shape[0], dtype=bool)

    if sparse.issparse(guide_assignment):
        return np.asarray(guide_assignment[:, guide_indices].sum(axis=1)).ravel() > 0

    return guide_assignment[:, guide_indices].any(axis=1)


def single_guide_mask(guide_assignment, guide_index):
    if sparse.issparse(guide_assignment):
        return np.asarray(guide_assignment[:, guide_index].todense()).ravel().astype(bool)
    return guide_assignment[:, guide_index].astype(bool)


def infer_target_col(meta):
    for c in TARGET_COL_CANDIDATES:
        if c in meta.columns:
            return c
    raise ValueError(
        "Could not find target column in guide metadata. "
        f"Available columns: {meta.columns.tolist()}"
    )


def is_false_like(x):
    return str(x).strip().lower() in ["false", "0", "no", "n"]


def is_hippo_target_string(x):
    targets = [t.strip().upper() for t in str(x).split(";") if t.strip()]
    return any(t in HIPPO_TARGETS_UPPER for t in targets)


def summarize_mask(obs, d3_mask, mask, row_info):
    use = d3_mask & mask

    n_cells = int(use.sum())
    n_EC = int((use & obs["is_EC"].values).sum())
    n_nonEC = n_cells - n_EC
    pct_EC = np.nan if n_cells == 0 else 100 * n_EC / n_cells

    out = dict(row_info)
    out.update({
        "n_D3_cells": n_cells,
        "n_D3_EC_cells": n_EC,
        "n_D3_nonEC_cells": n_nonEC,
        "pct_EC_D3": pct_EC,
    })
    return out


# ----------------------------
# Basic checks
# ----------------------------
if CELLTYPE_COL not in adata_int.obs.columns:
    raise ValueError(f"adata_int.obs is missing required column: {CELLTYPE_COL}")

if GUIDE_NAMES_KEY not in adata_int.uns:
    raise ValueError(f"adata_int.uns is missing: {GUIDE_NAMES_KEY}")

if GUIDE_ASSIGNMENT_KEY not in adata_int.obsm:
    raise ValueError(f"adata_int.obsm is missing: {GUIDE_ASSIGNMENT_KEY}")

if not Path(GUIDE_METADATA_TSV).exists():
    raise FileNotFoundError(
        f"Could not find {GUIDE_METADATA_TSV}. "
        "Put it in the notebook working directory or update GUIDE_METADATA_TSV."
    )

# ----------------------------
# Load guide metadata and map guide -> target
# ----------------------------
guide_names = list(map(str, adata_int.uns[GUIDE_NAMES_KEY]))
guide_name_to_idx = {g: i for i, g in enumerate(guide_names)}

meta = pd.read_csv(GUIDE_METADATA_TSV, sep="\t")
target_col = infer_target_col(meta)

if "guide_id" not in meta.columns:
    raise ValueError(f"guide_id column missing from {GUIDE_METADATA_TSV}")

meta = meta.copy()
meta["guide_id"] = meta["guide_id"].astype(str)
meta["raw_target_entry"] = meta[target_col]

# exact match to adata guide names
meta["guide_index"] = meta["guide_id"].map(guide_name_to_idx)
meta["match_method"] = np.where(meta["guide_index"].notna(), "exact", None)

# fallback canonical match
if meta["guide_index"].isna().any():
    canonical_to_guide = {}
    duplicated = set()

    for g in guide_names:
        cg = canonical(g)
        if cg in canonical_to_guide:
            duplicated.add(cg)
        else:
            canonical_to_guide[cg] = g

    for cg in duplicated:
        canonical_to_guide.pop(cg, None)

    for idx in meta.index[meta["guide_index"].isna()]:
        cg = canonical(meta.at[idx, "guide_id"])
        if cg in canonical_to_guide:
            matched_guide = canonical_to_guide[cg]
            meta.at[idx, "guide_index"] = guide_name_to_idx[matched_guide]
            meta.at[idx, "match_method"] = "canonical"

mapped_meta = meta[meta["guide_index"].notna()].copy()
mapped_meta["guide_index"] = mapped_meta["guide_index"].astype(int)
mapped_meta["guide_name"] = mapped_meta["guide_index"].map(lambda i: guide_names[i])

# Determine whether each metadata row is non-targeting/control.
if "targeting" in mapped_meta.columns:
    mapped_meta["is_control_by_metadata"] = mapped_meta["targeting"].map(is_false_like)
else:
    mapped_meta["is_control_by_metadata"] = False

records = []

for _, row in mapped_meta.iterrows():
    guide_id = row["guide_id"]
    guide_name = row["guide_name"]
    guide_index = int(row["guide_index"])
    raw_target = row["raw_target_entry"]

    is_control = (
        bool(row["is_control_by_metadata"])
        or looks_control(guide_id)
        or looks_control(guide_name)
        or looks_control(raw_target)
    )

    if is_control:
        targets = [CONTROL_TARGET_NAME]
    else:
        targets = clean_target_entry(raw_target)

    if len(targets) == 0:
        targets = ["unknown"]

    target_string = ";".join(sorted(set(targets)))

    records.append({
        "guide_id": guide_id,
        "guide_name": guide_name,
        "guide_index": guide_index,
        "target": target_string,
        "raw_target_entry": raw_target,
        "is_control_guide": is_control,
        "match_method": row["match_method"],
    })

guide_table = pd.DataFrame(records).drop_duplicates(subset=["guide_index"]).copy()
guide_table["is_hippo_guide"] = guide_table["target"].apply(is_hippo_target_string)
guide_table["is_targeting_guide"] = ~guide_table["is_control_guide"]

guide_table.to_csv(
    OUTDIR / "guide_to_target_mapping_used_for_d3_guide_export.csv",
    index=False
)

# ----------------------------
# Build guide assignment and masks
# ----------------------------
guide_assignment = get_guide_assignment_matrix(adata_int)

if guide_assignment.shape[1] != len(guide_names):
    raise ValueError(
        f"guide_assignment has {guide_assignment.shape[1]} columns, "
        f"but adata_int.uns['{GUIDE_NAMES_KEY}'] has {len(guide_names)} guide names."
    )

control_guide_indices = guide_table.loc[
    guide_table["is_control_guide"], "guide_index"
].astype(int).tolist()

targeting_guide_indices = guide_table.loc[
    guide_table["is_targeting_guide"], "guide_index"
].astype(int).tolist()

hippo_guide_indices = guide_table.loc[
    guide_table["is_hippo_guide"], "guide_index"
].astype(int).tolist()

control_mask_raw = any_guides_per_cell(guide_assignment, control_guide_indices)
targeting_mask_any = any_guides_per_cell(guide_assignment, targeting_guide_indices)
hippo_mask_any = any_guides_per_cell(guide_assignment, hippo_guide_indices)

# Main requested control denominator:
# cells with NT guide, excluding Hippo carriers only.
nt_excl_hippo_mask = control_mask_raw & (~hippo_mask_any)

# Optional strict NT, off by default.
strict_nt_mask = control_mask_raw & (~targeting_mask_any)

# ----------------------------
# Cell metadata
# ----------------------------
obs = adata_int.obs.copy()
obs["cell_id"] = adata_int.obs_names.astype(str)
obs["day_clean"] = clean_day_from_biosample_or_day(obs)
obs["cell_type_clean"] = obs[CELLTYPE_COL].astype(str)
obs["is_EC"] = obs["cell_type_clean"].map(classify_ec)

d3_mask = obs["day_clean"].values == "D3"

print("D3 cell-type labels:")
print(sorted(obs.loc[d3_mask, "cell_type_clean"].unique()))

print("\nD3 EC/non-EC counts:")
display(
    obs.loc[d3_mask]
    .assign(EC_status=lambda x: np.where(x["is_EC"], "EC", "non-EC"))
    .groupby("EC_status")
    .size()
    .rename("n_cells")
    .reset_index()
)

# ----------------------------
# Per-guide EC percentage
# ----------------------------
rows = []

# 1. Targeting guide rows.
targeting_guides = guide_table[guide_table["is_targeting_guide"]].copy()

for _, guide in targeting_guides.iterrows():
    guide_idx = int(guide["guide_index"])
    guide_mask = single_guide_mask(guide_assignment, guide_idx)

    rows.append(
        summarize_mask(
            obs=obs,
            d3_mask=d3_mask,
            mask=guide_mask,
            row_info={
                "row_type": "targeting_guide",
                "guide_id": guide["guide_id"],
                "guide_name": guide["guide_name"],
                "guide_index": guide_idx,
                "target": guide["target"],
                "control_definition": "",
                "exclusion_rule": "none",
            },
        )
    )

# 2. Individual NT guide rows, excluding Hippo carriers.
control_guides = guide_table[guide_table["is_control_guide"]].copy()

for _, guide in control_guides.iterrows():
    guide_idx = int(guide["guide_index"])
    guide_mask = single_guide_mask(guide_assignment, guide_idx)

    # This is the key change:
    # individual NT guide denominator excludes Hippo carriers only.
    use_mask = guide_mask & (~hippo_mask_any)

    rows.append(
        summarize_mask(
            obs=obs,
            d3_mask=d3_mask,
            mask=use_mask,
            row_info={
                "row_type": "individual_NT_guide_excl_Hippo",
                "guide_id": guide["guide_id"],
                "guide_name": guide["guide_name"],
                "guide_index": guide_idx,
                "target": "individual_NT_guide",
                "control_definition": "individual NT guide, excluding Hippo carriers",
                "exclusion_rule": "exclude_Hippo_only",
            },
        )
    )

# 3. Aggregate NT excluding Hippo row.
rows.append(
    summarize_mask(
        obs=obs,
        d3_mask=d3_mask,
        mask=nt_excl_hippo_mask,
        row_info={
            "row_type": "control_group",
            "guide_id": "",
            "guide_name": "",
            "guide_index": np.nan,
            "target": "NT excluding Hippo",
            "control_definition": "all NT/control-guide cells, excluding Hippo carriers",
            "exclusion_rule": "exclude_Hippo_only",
        },
    )
)

# 4. Optional strict NT aggregate row.
if EXPORT_STRICT_NT_AGGREGATE:
    rows.append(
        summarize_mask(
            obs=obs,
            d3_mask=d3_mask,
            mask=strict_nt_mask,
            row_info={
                "row_type": "control_group",
                "guide_id": "",
                "guide_name": "",
                "guide_index": np.nan,
                "target": "Strict NT / NT excluding all target",
                "control_definition": "NT/control-guide cells with no targeting guide",
                "exclusion_rule": "exclude_all_targeting_guides",
            },
        )
    )

guide_ec = pd.DataFrame(rows)

# Sort:
# targeting guides first by target,
# then individual NT guides,
# then aggregate controls.
row_order = {
    "targeting_guide": 0,
    "individual_NT_guide_excl_Hippo": 1,
    "control_group": 2,
}
guide_ec["sort_group"] = guide_ec["row_type"].map(row_order).fillna(99)

guide_ec = (
    guide_ec
    .sort_values(["sort_group", "target", "guide_id"], na_position="last")
    .drop(columns=["sort_group"])
)

main_out = OUTDIR / "d3_EC_percent_by_guide_with_individual_NT_excl_Hippo.csv"
guide_ec.to_csv(main_out, index=False)

print("Saved main guide-level table:", main_out)
display(guide_ec.head())

print("\nIndividual NT guide rows:")
display(guide_ec[guide_ec["row_type"] == "individual_NT_guide_excl_Hippo"])

print("\nAggregate control rows:")
display(guide_ec[guide_ec["row_type"] == "control_group"])

# ----------------------------
# Optional per-cell membership table
# ----------------------------
if EXPORT_CELL_MEMBERSHIP:
    membership_rows = []

    for _, guide in targeting_guides.iterrows():
        guide_idx = int(guide["guide_index"])
        guide_mask = single_guide_mask(guide_assignment, guide_idx)
        use = d3_mask & guide_mask

        tmp = obs.loc[use, ["cell_id", "cell_type_clean", "is_EC"]].copy()
        tmp["row_type"] = "targeting_guide"
        tmp["guide_id"] = guide["guide_id"]
        tmp["guide_name"] = guide["guide_name"]
        tmp["guide_index"] = guide_idx
        tmp["target"] = guide["target"]
        tmp["exclusion_rule"] = "none"
        membership_rows.append(tmp)

    for _, guide in control_guides.iterrows():
        guide_idx = int(guide["guide_index"])
        guide_mask = single_guide_mask(guide_assignment, guide_idx)
        use = d3_mask & guide_mask & (~hippo_mask_any)

        tmp = obs.loc[use, ["cell_id", "cell_type_clean", "is_EC"]].copy()
        tmp["row_type"] = "individual_NT_guide_excl_Hippo"
        tmp["guide_id"] = guide["guide_id"]
        tmp["guide_name"] = guide["guide_name"]
        tmp["guide_index"] = guide_idx
        tmp["target"] = "individual_NT_guide"
        tmp["exclusion_rule"] = "exclude_Hippo_only"
        membership_rows.append(tmp)

    use = d3_mask & nt_excl_hippo_mask
    tmp = obs.loc[use, ["cell_id", "cell_type_clean", "is_EC"]].copy()
    tmp["row_type"] = "control_group"
    tmp["guide_id"] = ""
    tmp["guide_name"] = ""
    tmp["guide_index"] = np.nan
    tmp["target"] = "NT excluding Hippo"
    tmp["exclusion_rule"] = "exclude_Hippo_only"
    membership_rows.append(tmp)

    membership = pd.concat(membership_rows, ignore_index=True)
    membership_out = OUTDIR / "d3_EC_percent_by_guide_with_individual_NT_excl_Hippo_cell_membership.csv"
    membership.to_csv(membership_out, index=False)
    print("Saved cell membership table:", membership_out)

# ----------------------------
# QC summary
# ----------------------------
qc_rows = []

qc_rows.append({
    "metric": "n_adata_cells",
    "value": adata_int.n_obs,
})

qc_rows.append({
    "metric": "n_D3_cells",
    "value": int(d3_mask.sum()),
})

qc_rows.append({
    "metric": "n_guide_names",
    "value": len(guide_names),
})

qc_rows.append({
    "metric": "n_metadata_guides",
    "value": len(meta),
})

qc_rows.append({
    "metric": "n_metadata_guides_mapped_to_adata_guides",
    "value": len(mapped_meta),
})

qc_rows.append({
    "metric": "n_targeting_guides_exported",
    "value": len(targeting_guides),
})

qc_rows.append({
    "metric": "n_individual_NT_guides_exported",
    "value": len(control_guides),
})

qc_rows.append({
    "metric": "n_control_guides_detected",
    "value": len(control_guide_indices),
})

qc_rows.append({
    "metric": "n_hippo_guides_detected",
    "value": len(hippo_guide_indices),
})

qc_rows.append({
    "metric": "n_D3_NT_excluding_Hippo_cells",
    "value": int((d3_mask & nt_excl_hippo_mask).sum()),
})

qc_rows.append({
    "metric": "n_D3_Strict_NT_excluding_all_target_cells_optional",
    "value": int((d3_mask & strict_nt_mask).sum()),
})

qc_df = pd.DataFrame(qc_rows)
qc_out = OUTDIR / "d3_EC_percent_by_guide_export_qc_summary.csv"
qc_df.to_csv(qc_out, index=False)

print("Saved QC summary:", qc_out)
display(qc_df)
```

### Text output (stream)

```text
D3 cell-type labels:
['Fibro-Mesen', 'Lat-Meso', 'Neu-Mesen', 'Neuro-iPSC', 'Others', 'Pre-EC', 'PrimStreak', 'artEC', 'iPSC']

D3 EC/non-EC counts:

```

### Text output (text/plain)

```text
  EC_status  n_cells
0        EC    28444
1    non-EC    15658
```

### Text output (stream)

```text
Saved main guide-level table: d3_EC_percent_by_guide/d3_EC_percent_by_guide_with_individual_NT_excl_Hippo.csv

```

### Text output (text/plain)

```text
          row_type                  guide_id                guide_name  \
0  targeting_guide  ACAA1_+_38178488.23-P1P2  ACAA1_+_38178488.23-P1P2   
1  targeting_guide  ACAA1_+_38178517.23-P1P2  ACAA1_+_38178517.23-P1P2   
2  targeting_guide  ACAA1_+_38178559.23-P1P2  ACAA1_+_38178559.23-P1P2   
3  targeting_guide  ACAA1_+_38178570.23-P1P2  ACAA1_+_38178570.23-P1P2   
4  targeting_guide  ACAA1_+_38178620.23-P1P2  ACAA1_+_38178620.23-P1P2   

   guide_index target control_definition exclusion_rule  n_D3_cells  \
0          0.0  ACAA1                              none          87   
1          1.0  ACAA1                              none          99   
2          2.0  ACAA1                              none          72   
3          3.0  ACAA1                              none         110   
4          4.0  ACAA1                              none         136   

   n_D3_EC_cells  n_D3_nonEC_cells  pct_EC_D3  
0             53                34  60.919540  
1             55                44  55.555556  
2             44                28  61.111111  
3             64                46  58.181818  
4             82                54  60.294118  
```

### Text output (stream)

```text

Individual NT guide rows:

```

### Text output (text/plain)

```text
                            row_type             guide_id  \
1862  individual_NT_guide_excl_Hippo  SAFE_TARGETING_1017   
1863  individual_NT_guide_excl_Hippo  SAFE_TARGETING_1044   
1864  individual_NT_guide_excl_Hippo  SAFE_TARGETING_1105   
1865  individual_NT_guide_excl_Hippo  SAFE_TARGETING_1141   
1866  individual_NT_guide_excl_Hippo  SAFE_TARGETING_1187   
...                              ...                  ...   
1857  individual_NT_guide_excl_Hippo  non-targeting_03427   
1858  individual_NT_guide_excl_Hippo  non-targeting_03483   
1859  individual_NT_guide_excl_Hippo  non-targeting_03620   
1860  individual_NT_guide_excl_Hippo  non-targeting_03635   
1861  individual_NT_guide_excl_Hippo  non-targeting_03664   

               guide_name  guide_index               target  \
1862  SAFE_TARGETING_1017       1268.0  individual_NT_guide   
1863  SAFE_TARGETING_1044       1269.0  individual_NT_guide   
1864  SAFE_TARGETING_1105       1270.0  individual_NT_guide   
1865  SAFE_TARGETING_1141       1271.0  individual_NT_guide   
1866  SAFE_TARGETING_1187       1272.0  individual_NT_guide   
...                   ...          ...                  ...   
1857  non-targeting_03427       1917.0  individual_NT_guide   
1858  non-targeting_03483       1918.0  individual_NT_guide   
1859  non-targeting_03620       1919.0  individual_NT_guide   
1860  non-targeting_03635       1920.0  individual_NT_guide   
1861  non-targeting_03664       1921.0  individual_NT_guide   

                                 control_definition      exclusion_rule  \
1862  individual NT guide, excluding Hippo carriers  exclude_Hippo_only   
1863  individual NT guide, excluding Hippo carriers  exclude_Hippo_only   
1864  individual NT guide, excluding Hippo carriers  exclude_Hippo_only   
1865  individual NT guide, excluding Hippo carriers  exclude_Hippo_only   
1866  individual NT guide, excluding Hippo carriers  exclude_Hippo_only   
...                                             ...                 ...   
1857  individual NT guide, excluding Hippo carriers  exclude_Hippo_only   
1858  individual NT guide, excluding Hippo carriers  exclude_Hippo_only   
1859  individual NT guide, excluding Hippo carriers  exclude_Hippo_only   
1860  individual NT guide, excluding Hippo carriers  exclude_Hippo_only   
1861  individual NT guide, excluding Hippo carriers  exclude_Hippo_only   

      n_D3_cells  n_D3_EC_cells  n_D3_nonEC_cells  pct_EC_D3  
1862          32             18                14  56.250000  
1863          28             15                13  53.571429  
1864          44             23                21  52.272727  
1865          46             28                18  60.869565  
1866          50             29                21  58.000000  
...          ...            ...               ...        ...  
1857          39             23                16  58.974359  
1858          28             15                13  53.571429  
1859          96             53                43  55.208333  
1860          95             52                43  54.736842  
1861          86             46                40  53.488372  

[120 rows x 11 columns]
```

### Text output (stream)

```text

Aggregate control rows:

```

### Text output (text/plain)

```text
           row_type guide_id guide_name  guide_index              target  \
1922  control_group                              NaN  NT excluding Hippo   

                                     control_definition      exclusion_rule  \
1922  all NT/control-guide cells, excluding Hippo ca...  exclude_Hippo_only   

      n_D3_cells  n_D3_EC_cells  n_D3_nonEC_cells  pct_EC_D3  
1922        6308           3670              2638  58.180089  
```

### Text output (stream)

```text
Saved QC summary: d3_EC_percent_by_guide/d3_EC_percent_by_guide_export_qc_summary.csv

```

### Text output (text/plain)

```text
                                               metric   value
0                                       n_adata_cells  121616
1                                          n_D3_cells   44102
2                                       n_guide_names    1922
3                                   n_metadata_guides    1926
4            n_metadata_guides_mapped_to_adata_guides    1922
5                         n_targeting_guides_exported    1802
6                     n_individual_NT_guides_exported     120
7                           n_control_guides_detected     120
8                             n_hippo_guides_detected      30
9                       n_D3_NT_excluding_Hippo_cells    6308
10  n_D3_Strict_NT_excluding_all_target_cells_opti...     378
```

