# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/qq_diff_approaches.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/qq_diff_approaches.ipynb`
- SHA-256: `7537b5afb3a07ca2c82ea4951130637ac570250061bf10e78f28c7213f3c0918`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Code cell 0

```python
import sys
print(sys.executable)
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python

```

## Code cell 1

```python
import scanpy as sc
import muon as mu
import pandas as pd
import numpy as np

```

## Code cell 2

```python
import sys
sys.path.append("/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE")
from tests.synthetic_data import make_sceptre_style_synth
from src.sceptre import prepare_crt_inputs
```

## Code cell 3

```python
mdata = mu.read("/oak/stanford/groups/engreitz/Users/ymo/cc-perturb-seq/Results/111025_D0_IGVF_10iter_torch_halsvar_batch_e7_v100s_test/adata/cNMF_100_2_0.h5mu")
mdata_guide = mu.read("/oak/stanford/groups/engreitz/Users/ymo/cc-perturb-seq/Data/IGVF_D0_example.h5mu")
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

## Code cell 4

```python
mdata
```

### Text output (text/plain)

```text
MuData object with n_obs × n_vars = 51979 × 7846
  2 modalities
    rna:	51979 x 7646
      obs:	'batch', 'cov1', 'batch_number', 'n_counts', 'log1p_n_genes_by_counts', 'total_gene_umis', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'percent_mito', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'num_expressed_genes'
      var:	'symbol', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'gene_chr', 'gene_start', 'gene_end'
      obsm:	'guide_assignment'
    cNMF:	51979 x 200
      obs:	'batch', 'cov1', 'batch_number', 'n_counts', 'log1p_n_genes_by_counts', 'total_gene_umis', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'percent_mito', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'num_expressed_genes'
      uns:	'var_names'
      varm:	'loadings'
```

## Code cell 5

```python
mdata['rna'].obs["num_expressed_genes"]
```

### Text output (text/plain)

```text
AAGGTACACATACCAATTACCTCG_0    2393
ACACAGAAGCGAGTAAATAAGCTC_0    1209
CGCTGATCGACTAGTAATTCATGG_0     472
GCCAAGACCACCTTACCCGTTCTA_0     627
ATCATTCCGTCGTAGAATTCATGG_0     893
                              ... 
CGACACACAGAGTCAAATCCGCGA_2     852
TGAAGAGACAGCGTTATGGCGCGC_2    1010
AGGCTAACGAGTTAGCTGTCTGAA_2     335
AATGTTGCAGTCACTATGGTATAC_2    1091
GTCTGTCACAAGGAGCTGGCGCGC_2     125
Name: num_expressed_genes, Length: 51979, dtype: int64
```

## Code cell 6

```python
mdata_guide
```

### Text output (text/plain)

```text
MuData object with n_obs × n_vars = 51979 × 9569
  obs:	'cov1', 'batch'
  uns:	'per_element_results', 'per_guide_results'
  2 modalities
    gene:	51979 x 7646
      obs:	'batch', 'cov1', 'batch_number', 'n_counts', 'log1p_n_genes_by_counts', 'total_gene_umis', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'percent_mito', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'num_expressed_genes'
      var:	'symbol', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'gene_chr', 'gene_start', 'gene_end'
      obsm:	'guide_assignment'
    guide:	51979 x 1923
      obs:	'batch', 'cov1', 'num_expressed_guides', 'total_guide_umis'
      var:	'guide_id', 'spacer', 'targeting', 'type', 'guide_chr', 'guide_start', 'guide_end', 'strand', 'pam', 'intended_target_name', 'intended_target_chr', 'intended_target_start', 'intended_target_end'
      uns:	'capture_method', 'moi'
      layers:	'guide_assignment'
```

## Code cell 7

```python
mdata_guide['guide'].obs['num_expressed_guides']
```

### Text output (text/plain)

```text
AAGGTACACATACCAATTACCTCG_0    188
ACACAGAAGCGAGTAAATAAGCTC_0     22
CGCTGATCGACTAGTAATTCATGG_0      1
GCCAAGACCACCTTACCCGTTCTA_0     33
ATCATTCCGTCGTAGAATTCATGG_0     39
                             ... 
CGACACACAGAGTCAAATCCGCGA_2     75
TGAAGAGACAGCGTTATGGCGCGC_2     34
AGGCTAACGAGTTAGCTGTCTGAA_2     22
AATGTTGCAGTCACTATGGTATAC_2     19
GTCTGTCACAAGGAGCTGGCGCGC_2      3
Name: num_expressed_guides, Length: 51979, dtype: int64
```

## Code cell 8

```python
adata = mdata["cNMF"].copy()
adata.obsm["cnmf_usage"] = np.asarray(adata.X)  # ensure dense float

# Align gene modality to adata cells
g = mdata_guide["gene"][adata.obs_names].copy()

# Guide assignment (cell x guide)
adata.obsm["guide_assignment"] = g.obsm["guide_assignment"].copy()

# Program names
adata.uns["program_names"] = list(adata.var_names)

# Guide names must match columns of guide_assignment
guide_names = list(mdata_guide["guide"].var["guide_id"])
adata.uns["guide_names"] = guide_names

# guide2gene: map guide_name -> gene by taking the first token of guide_name.split('_')
guide2gene = {gn: gn.split("_", 1)[0] for gn in guide_names}
adata.uns["guide2gene"] = guide2gene


```

## Code cell 9

```python
g.obs["batch"]
```

### Text output (text/plain)

```text
AAGGTACACATACCAATTACCTCG_0    1
ACACAGAAGCGAGTAAATAAGCTC_0    1
CGCTGATCGACTAGTAATTCATGG_0    1
GCCAAGACCACCTTACCCGTTCTA_0    1
ATCATTCCGTCGTAGAATTCATGG_0    1
                             ..
CGACACACAGAGTCAAATCCGCGA_2    3
TGAAGAGACAGCGTTATGGCGCGC_2    3
AGGCTAACGAGTTAGCTGTCTGAA_2    3
AATGTTGCAGTCACTATGGTATAC_2    3
GTCTGTCACAAGGAGCTGGCGCGC_2    3
Name: batch, Length: 51979, dtype: category
Categories (3, object): ['1', '2', '3']
```

## Code cell 10

```python
mdata['rna'].obs["percent_mito"]
```

### Text output (text/plain)

```text
AAGGTACACATACCAATTACCTCG_0     3.647181
ACACAGAAGCGAGTAAATAAGCTC_0     7.251676
CGCTGATCGACTAGTAATTCATGG_0     9.743590
GCCAAGACCACCTTACCCGTTCTA_0     8.147175
ATCATTCCGTCGTAGAATTCATGG_0     6.695279
                                ...    
CGACACACAGAGTCAAATCCGCGA_2     6.390978
TGAAGAGACAGCGTTATGGCGCGC_2     9.084355
AGGCTAACGAGTTAGCTGTCTGAA_2    13.793103
AATGTTGCAGTCACTATGGTATAC_2    12.095731
GTCTGTCACAAGGAGCTGGCGCGC_2     3.030303
Name: percent_mito, Length: 51979, dtype: float32
```

## Code cell 11

```python
# Covariates

#(~ log(response_n_nonzero) + log(response_n_umis) +
#         log(grna_n_nonzero) + log(grna_n_umis) + response_p_mito + batch)



covar_df = pd.DataFrame({"response_n_nonzero": mdata['rna'].obs["log1p_n_genes_by_counts"],
            "response_n_umis": mdata['rna'].obs["log1p_total_counts"],
            "grna_n_nonzero": np.log1p(mdata_guide['guide'].obs["num_expressed_guides"]),
            "grna_n_umis": np.log1p(mdata_guide['guide'].obs["total_guide_umis"]),
            "response_p_mito": mdata['rna'].obs["percent_mito"],
            "batch": g.obs["batch"],
           })
covar_df
```

### Text output (text/plain)

```text
                            response_n_nonzero  response_n_umis  \
AAGGTACACATACCAATTACCTCG_0            7.780721         8.260234   
ACACAGAAGCGAGTAAATAAGCTC_0            7.098376         7.403670   
CGCTGATCGACTAGTAATTCATGG_0            6.159095         6.373320   
GCCAAGACCACCTTACCCGTTCTA_0            6.442540         6.635947   
ATCATTCCGTCGTAGAATTCATGG_0            6.795706         7.061335   
...                                        ...              ...   
CGACACACAGAGTCAAATCCGCGA_2            6.748760         6.970730   
TGAAGAGACAGCGTTATGGCGCGC_2            6.918695         7.235619   
AGGCTAACGAGTTAGCTGTCTGAA_2            5.817111         6.008813   
AATGTTGCAGTCACTATGGTATAC_2            6.995766         7.344073   
GTCTGTCACAAGGAGCTGGCGCGC_2            4.836282         4.890349   

                            grna_n_nonzero  grna_n_umis  response_p_mito batch  
AAGGTACACATACCAATTACCTCG_0        5.241747     7.264730         3.647181     1  
ACACAGAAGCGAGTAAATAAGCTC_0        3.135494     4.406719         7.251676     1  
CGCTGATCGACTAGTAATTCATGG_0        0.693147     0.693147         9.743590     1  
GCCAAGACCACCTTACCCGTTCTA_0        3.526361     5.332719         8.147175     1  
ATCATTCCGTCGTAGAATTCATGG_0        3.688879     4.941642         6.695279     1  
...                                    ...          ...              ...   ...  
CGACACACAGAGTCAAATCCGCGA_2        4.330733     6.579251         6.390978     3  
TGAAGAGACAGCGTTATGGCGCGC_2        3.555348     5.356586         9.084355     3  
AGGCTAACGAGTTAGCTGTCTGAA_2        3.135494     4.521789        13.793103     3  
AATGTTGCAGTCACTATGGTATAC_2        2.995732     3.737670        12.095731     3  
GTCTGTCACAAGGAGCTGGCGCGC_2        1.386294     1.386294         3.030303     3  

[51979 rows x 6 columns]
```

## Code cell 12

```python
import scipy.sparse as sp

# Inputs
G = adata.obsm["guide_assignment"]  # cells x guides, sparse/binary
guide2gene = adata.uns["guide2gene"]  # dict: guide_name -> gene_label
ntc_labels = {"SAFE", "non-targeting"}  # set for fast lookup

# Ensure CSR
if not sp.issparse(G):
    raise ValueError("adata.obsm['guide_assignment'] must be a sparse matrix")
G = G.tocsr()

# Guide names must align with columns of G
# Prefer storing explicitly; fallback to adata.var_names
guide_names = adata.uns.get("guide_names", list(adata.var_names))
if len(guide_names) != G.shape[1]:
    raise ValueError("guide_names length must match guide_assignment n_guides")

# Identify NTC guide columns
is_ntc = np.array([guide2gene.get(g, None) in ntc_labels for g in guide_names], dtype=bool)
ntc_cols = np.flatnonzero(is_ntc)
non_ntc_cols = np.flatnonzero(~is_ntc)

# Per-cell counts (works with binary sparse)
ntc_guides_per_cell = np.asarray(G[:, ntc_cols].getnnz(axis=1)).ravel().astype(np.float32)
non_ntc_guides_per_cell = np.asarray(G[:, non_ntc_cols].getnnz(axis=1)).ravel().astype(np.float32)

# log1p transforms
ntc_guides_per_cell_log1p = np.log1p(ntc_guides_per_cell)
non_ntc_guides_per_cell_log1p = np.log1p(non_ntc_guides_per_cell)

covar_df["ntc_guides_per_cell"] = ntc_guides_per_cell_log1p
covar_df["non_ntc_guides_per_cell"] = non_ntc_guides_per_cell_log1p

adata.obsm["covar"] = covar_df

print(
    f"Added covariates. NTC guides: {len(ntc_cols)} / {G.shape[1]} "
    f"(labels={sorted(ntc_labels)})"
)

```

### Text output (stream)

```text
Added covariates. NTC guides: 119 / 1923 (labels=['SAFE', 'non-targeting'])

```

## Code cell 13

```python
adata.obsm["covar"] = covar_df
```

## Code cell 14

```python
adata.obsm["covar"]
```

### Text output (text/plain)

```text
                            response_n_nonzero  response_n_umis  \
AAGGTACACATACCAATTACCTCG_0            7.780721         8.260234   
ACACAGAAGCGAGTAAATAAGCTC_0            7.098376         7.403670   
CGCTGATCGACTAGTAATTCATGG_0            6.159095         6.373320   
GCCAAGACCACCTTACCCGTTCTA_0            6.442540         6.635947   
ATCATTCCGTCGTAGAATTCATGG_0            6.795706         7.061335   
...                                        ...              ...   
CGACACACAGAGTCAAATCCGCGA_2            6.748760         6.970730   
TGAAGAGACAGCGTTATGGCGCGC_2            6.918695         7.235619   
AGGCTAACGAGTTAGCTGTCTGAA_2            5.817111         6.008813   
AATGTTGCAGTCACTATGGTATAC_2            6.995766         7.344073   
GTCTGTCACAAGGAGCTGGCGCGC_2            4.836282         4.890349   

                            grna_n_nonzero  grna_n_umis  response_p_mito  \
AAGGTACACATACCAATTACCTCG_0        5.241747     7.264730         3.647181   
ACACAGAAGCGAGTAAATAAGCTC_0        3.135494     4.406719         7.251676   
CGCTGATCGACTAGTAATTCATGG_0        0.693147     0.693147         9.743590   
GCCAAGACCACCTTACCCGTTCTA_0        3.526361     5.332719         8.147175   
ATCATTCCGTCGTAGAATTCATGG_0        3.688879     4.941642         6.695279   
...                                    ...          ...              ...   
CGACACACAGAGTCAAATCCGCGA_2        4.330733     6.579251         6.390978   
TGAAGAGACAGCGTTATGGCGCGC_2        3.555348     5.356586         9.084355   
AGGCTAACGAGTTAGCTGTCTGAA_2        3.135494     4.521789        13.793103   
AATGTTGCAGTCACTATGGTATAC_2        2.995732     3.737670        12.095731   
GTCTGTCACAAGGAGCTGGCGCGC_2        1.386294     1.386294         3.030303   

                           batch  ntc_guides_per_cell  non_ntc_guides_per_cell  
AAGGTACACATACCAATTACCTCG_0     1             0.693147                 2.708050  
ACACAGAAGCGAGTAAATAAGCTC_0     1             0.000000                 0.693147  
CGCTGATCGACTAGTAATTCATGG_0     1             0.000000                 0.000000  
GCCAAGACCACCTTACCCGTTCTA_0     1             0.000000                 1.945910  
ATCATTCCGTCGTAGAATTCATGG_0     1             0.000000                 1.098612  
...                          ...                  ...                      ...  
CGACACACAGAGTCAAATCCGCGA_2     3             0.000000                 2.564949  
TGAAGAGACAGCGTTATGGCGCGC_2     3             0.000000                 1.386294  
AGGCTAACGAGTTAGCTGTCTGAA_2     3             0.000000                 2.197225  
AATGTTGCAGTCACTATGGTATAC_2     3             0.000000                 1.098612  
GTCTGTCACAAGGAGCTGGCGCGC_2     3             0.000000                 0.000000  

[51979 rows x 8 columns]
```

## Code cell 15

```python
import numpy as np
import scipy.sparse as sp
import pandas as pd

G = adata.obsm["guide_assignment"]
if not sp.issparse(G):
    raise ValueError("guide_assignment must be sparse")
G = G.tocsr()

guide2gene = adata.uns["guide2gene"]
ntc_labels = {"SAFE", "non-targeting"}

guide_names = adata.uns.get("guide_names", list(adata.var_names))
if len(guide_names) != G.shape[1]:
    raise ValueError("guide_names length must match guide_assignment columns")

is_ntc = np.array([guide2gene.get(g) in ntc_labels for g in guide_names], dtype=bool)
ntc_cols = np.flatnonzero(is_ntc)
non_cols = np.flatnonzero(~is_ntc)

# RAW counts (do NOT log yet)
ntc_cnt = np.asarray(G[:, ntc_cols].getnnz(axis=1)).ravel()
non_cnt = np.asarray(G[:, non_cols].getnnz(axis=1)).ravel()
total_cnt = np.asarray(G.getnnz(axis=1)).ravel()

# This must hold if the split is correct
if not np.all(total_cnt == ntc_cnt + non_cnt):
    bad = np.where(total_cnt != ntc_cnt + non_cnt)[0][:10]
    raise RuntimeError(f"Split inconsistent for rows (first 10): {bad}")

# Add BOTH raw and log1p versions (recommended)
covar = adata.obsm["covar"]
covar_df = covar.copy() if isinstance(covar, pd.DataFrame) else pd.DataFrame(covar, index=adata.obs_names)

covar_df["ntc_guides_per_cell"] = ntc_cnt.astype(np.float32)
covar_df["non_ntc_guides_per_cell"] = non_cnt.astype(np.float32)
covar_df["grna_n_nonzero_from_assignment"] = total_cnt.astype(np.float32)

covar_df["log1p_ntc_guides_per_cell"] = np.log1p(ntc_cnt).astype(np.float32)
covar_df["log1p_non_ntc_guides_per_cell"] = np.log1p(non_cnt).astype(np.float32)
covar_df["log1p_grna_n_nonzero_from_assignment"] = np.log1p(total_cnt).astype(np.float32)

adata.obsm["covar"] = covar_df

# Compare to your existing grna_n_nonzero (likely from mdata_guide)
if "grna_n_nonzero" in covar_df.columns:
    a = np.expm1(covar_df["grna_n_nonzero"].to_numpy())          # what it implies in count space
    b = total_cnt.astype(float)
    print("expm1(grna_n_nonzero) vs total_cnt (assignment):")
    print("  corr =", np.corrcoef(a, b)[0, 1])
    print("  median ratio a/b =", np.median((a + 1e-9) / (b + 1e-9)))
    print("  median abs diff =", np.median(np.abs(a - b)))

```

### Text output (stream)

```text
expm1(grna_n_nonzero) vs total_cnt (assignment):
  corr = 0.6655931397223657
  median ratio a/b = 5.799999999040001
  median abs diff = 14.000000000000004

```

## Code cell 16

```python
adata.obsm["covar"]
```

### Text output (text/plain)

```text
                            response_n_nonzero  response_n_umis  \
AAGGTACACATACCAATTACCTCG_0            7.780721         8.260234   
ACACAGAAGCGAGTAAATAAGCTC_0            7.098376         7.403670   
CGCTGATCGACTAGTAATTCATGG_0            6.159095         6.373320   
GCCAAGACCACCTTACCCGTTCTA_0            6.442540         6.635947   
ATCATTCCGTCGTAGAATTCATGG_0            6.795706         7.061335   
...                                        ...              ...   
CGACACACAGAGTCAAATCCGCGA_2            6.748760         6.970730   
TGAAGAGACAGCGTTATGGCGCGC_2            6.918695         7.235619   
AGGCTAACGAGTTAGCTGTCTGAA_2            5.817111         6.008813   
AATGTTGCAGTCACTATGGTATAC_2            6.995766         7.344073   
GTCTGTCACAAGGAGCTGGCGCGC_2            4.836282         4.890349   

                            grna_n_nonzero  grna_n_umis  response_p_mito  \
AAGGTACACATACCAATTACCTCG_0        5.241747     7.264730         3.647181   
ACACAGAAGCGAGTAAATAAGCTC_0        3.135494     4.406719         7.251676   
CGCTGATCGACTAGTAATTCATGG_0        0.693147     0.693147         9.743590   
GCCAAGACCACCTTACCCGTTCTA_0        3.526361     5.332719         8.147175   
ATCATTCCGTCGTAGAATTCATGG_0        3.688879     4.941642         6.695279   
...                                    ...          ...              ...   
CGACACACAGAGTCAAATCCGCGA_2        4.330733     6.579251         6.390978   
TGAAGAGACAGCGTTATGGCGCGC_2        3.555348     5.356586         9.084355   
AGGCTAACGAGTTAGCTGTCTGAA_2        3.135494     4.521789        13.793103   
AATGTTGCAGTCACTATGGTATAC_2        2.995732     3.737670        12.095731   
GTCTGTCACAAGGAGCTGGCGCGC_2        1.386294     1.386294         3.030303   

                           batch  ntc_guides_per_cell  \
AAGGTACACATACCAATTACCTCG_0     1                  1.0   
ACACAGAAGCGAGTAAATAAGCTC_0     1                  0.0   
CGCTGATCGACTAGTAATTCATGG_0     1                  0.0   
GCCAAGACCACCTTACCCGTTCTA_0     1                  0.0   
ATCATTCCGTCGTAGAATTCATGG_0     1                  0.0   
...                          ...                  ...   
CGACACACAGAGTCAAATCCGCGA_2     3                  0.0   
TGAAGAGACAGCGTTATGGCGCGC_2     3                  0.0   
AGGCTAACGAGTTAGCTGTCTGAA_2     3                  0.0   
AATGTTGCAGTCACTATGGTATAC_2     3                  0.0   
GTCTGTCACAAGGAGCTGGCGCGC_2     3                  0.0   

                            non_ntc_guides_per_cell  \
AAGGTACACATACCAATTACCTCG_0                     14.0   
ACACAGAAGCGAGTAAATAAGCTC_0                      1.0   
CGCTGATCGACTAGTAATTCATGG_0                      0.0   
GCCAAGACCACCTTACCCGTTCTA_0                      6.0   
ATCATTCCGTCGTAGAATTCATGG_0                      2.0   
...                                             ...   
CGACACACAGAGTCAAATCCGCGA_2                     12.0   
TGAAGAGACAGCGTTATGGCGCGC_2                      3.0   
AGGCTAACGAGTTAGCTGTCTGAA_2                      8.0   
AATGTTGCAGTCACTATGGTATAC_2                      2.0   
GTCTGTCACAAGGAGCTGGCGCGC_2                      0.0   

                            grna_n_nonzero_from_assignment  \
AAGGTACACATACCAATTACCTCG_0                            15.0   
ACACAGAAGCGAGTAAATAAGCTC_0                             1.0   
CGCTGATCGACTAGTAATTCATGG_0                             0.0   
GCCAAGACCACCTTACCCGTTCTA_0                             6.0   
ATCATTCCGTCGTAGAATTCATGG_0                             2.0   
...                                                    ...   
CGACACACAGAGTCAAATCCGCGA_2                            12.0   
TGAAGAGACAGCGTTATGGCGCGC_2                             3.0   
AGGCTAACGAGTTAGCTGTCTGAA_2                             8.0   
AATGTTGCAGTCACTATGGTATAC_2                             2.0   
GTCTGTCACAAGGAGCTGGCGCGC_2                             0.0   

                            log1p_ntc_guides_per_cell  \
AAGGTACACATACCAATTACCTCG_0                   0.693147   
ACACAGAAGCGAGTAAATAAGCTC_0                   0.000000   
CGCTGATCGACTAGTAATTCATGG_0                   0.000000   
GCCAAGACCACCTTACCCGTTCTA_0                   0.000000   
ATCATTCCGTCGTAGAATTCATGG_0                   0.000000   
...                                               ...   
CGACACACAGAGTCAAATCCGCGA_2                   0.000000   
TGAAGAGACAGCGTTATGGCGCGC_2                   0.000000   
AGGCTAACGAGTTAGCTGTCTGAA_2                   0.000000   
AATGTTGCAGTCACTATGGTATAC_2                   0.000000   
GTCTGTCACAAGGAGCTGGCGCGC_2                   0.000000   

                            log1p_non_ntc_guides_per_cell  \
AAGGTACACATACCAATTACCTCG_0                       2.708050   
ACACAGAAGCGAGTAAATAAGCTC_0                       0.693147   
CGCTGATCGACTAGTAATTCATGG_0                       0.000000   
GCCAAGACCACCTTACCCGTTCTA_0                       1.945910   
ATCATTCCGTCGTAGAATTCATGG_0                       1.098612   
...                                                   ...   
CGACACACAGAGTCAAATCCGCGA_2                       2.564949   
TGAAGAGACAGCGTTATGGCGCGC_2                       1.386294   
AGGCTAACGAGTTAGCTGTCTGAA_2                       2.197225   
AATGTTGCAGTCACTATGGTATAC_2                       1.098612   
GTCTGTCACAAGGAGCTGGCGCGC_2                       0.000000   

                            log1p_grna_n_nonzero_from_assignment  
AAGGTACACATACCAATTACCTCG_0                              2.772589  
ACACAGAAGCGAGTAAATAAGCTC_0                              0.693147  
CGCTGATCGACTAGTAATTCATGG_0                              0.000000  
GCCAAGACCACCTTACCCGTTCTA_0                              1.945910  
ATCATTCCGTCGTAGAATTCATGG_0                              1.098612  
...                                                          ...  
CGACACACAGAGTCAAATCCGCGA_2                              2.564949  
TGAAGAGACAGCGTTATGGCGCGC_2                              1.386294  
AGGCTAACGAGTTAGCTGTCTGAA_2                              2.197225  
AATGTTGCAGTCACTATGGTATAC_2                              1.098612  
GTCTGTCACAAGGAGCTGGCGCGC_2                              0.000000  

[51979 rows x 12 columns]
```

## Code cell 17

```python
covar = adata.obsm["covar"]
if not isinstance(covar, pd.DataFrame):
    covar = pd.DataFrame(covar, index=adata.obs_names)

# Columns to drop (keep log1p versions + keep grna_n_umis, drop raw counts + old grna_n_nonzero)
drop_cols = [
    "grna_n_nonzero",
    "ntc_guides_per_cell",
    "non_ntc_guides_per_cell",
    "grna_n_nonzero_from_assignment",
    "log1p_grna_n_nonzero_from_assignment",
]

missing = [c for c in drop_cols if c not in covar.columns]
if missing:
    print("Warning: these columns were not present and cannot be dropped:", missing)

covar = covar.drop(columns=[c for c in drop_cols if c in covar.columns])

adata.obsm["covar"] = covar
```

## Code cell 18

```python
adata.obsm["covar"]
```

### Text output (text/plain)

```text
                            response_n_nonzero  response_n_umis  grna_n_umis  \
AAGGTACACATACCAATTACCTCG_0            7.780721         8.260234     7.264730   
ACACAGAAGCGAGTAAATAAGCTC_0            7.098376         7.403670     4.406719   
CGCTGATCGACTAGTAATTCATGG_0            6.159095         6.373320     0.693147   
GCCAAGACCACCTTACCCGTTCTA_0            6.442540         6.635947     5.332719   
ATCATTCCGTCGTAGAATTCATGG_0            6.795706         7.061335     4.941642   
...                                        ...              ...          ...   
CGACACACAGAGTCAAATCCGCGA_2            6.748760         6.970730     6.579251   
TGAAGAGACAGCGTTATGGCGCGC_2            6.918695         7.235619     5.356586   
AGGCTAACGAGTTAGCTGTCTGAA_2            5.817111         6.008813     4.521789   
AATGTTGCAGTCACTATGGTATAC_2            6.995766         7.344073     3.737670   
GTCTGTCACAAGGAGCTGGCGCGC_2            4.836282         4.890349     1.386294   

                            response_p_mito batch  log1p_ntc_guides_per_cell  \
AAGGTACACATACCAATTACCTCG_0         3.647181     1                   0.693147   
ACACAGAAGCGAGTAAATAAGCTC_0         7.251676     1                   0.000000   
CGCTGATCGACTAGTAATTCATGG_0         9.743590     1                   0.000000   
GCCAAGACCACCTTACCCGTTCTA_0         8.147175     1                   0.000000   
ATCATTCCGTCGTAGAATTCATGG_0         6.695279     1                   0.000000   
...                                     ...   ...                        ...   
CGACACACAGAGTCAAATCCGCGA_2         6.390978     3                   0.000000   
TGAAGAGACAGCGTTATGGCGCGC_2         9.084355     3                   0.000000   
AGGCTAACGAGTTAGCTGTCTGAA_2        13.793103     3                   0.000000   
AATGTTGCAGTCACTATGGTATAC_2        12.095731     3                   0.000000   
GTCTGTCACAAGGAGCTGGCGCGC_2         3.030303     3                   0.000000   

                            log1p_non_ntc_guides_per_cell  
AAGGTACACATACCAATTACCTCG_0                       2.708050  
ACACAGAAGCGAGTAAATAAGCTC_0                       0.693147  
CGCTGATCGACTAGTAATTCATGG_0                       0.000000  
GCCAAGACCACCTTACCCGTTCTA_0                       1.945910  
ATCATTCCGTCGTAGAATTCATGG_0                       1.098612  
...                                                   ...  
CGACACACAGAGTCAAATCCGCGA_2                       2.564949  
TGAAGAGACAGCGTTATGGCGCGC_2                       1.386294  
AGGCTAACGAGTTAGCTGTCTGAA_2                       2.197225  
AATGTTGCAGTCACTATGGTATAC_2                       1.098612  
GTCTGTCACAAGGAGCTGGCGCGC_2                       0.000000  

[51979 rows x 7 columns]
```

## Code cell 19

```python
adata.obsm["covar"].columns
```

### Text output (text/plain)

```text
Index(['response_n_nonzero', 'response_n_umis', 'grna_n_umis',
       'response_p_mito', 'batch', 'log1p_ntc_guides_per_cell',
       'log1p_non_ntc_guides_per_cell'],
      dtype='object')
```

## Code cell 20

```python
U = adata.obsm["cnmf_usage"]
print("min:", U.min(), "quantile:", np.quantile(U, 1e-4))
print("row sums == 0:", np.sum(U.sum(axis=1) == 0))
```

### Text output (stream)

```text
min: 0.0 quantile: 0.0
row sums == 0: 0

```

## Code cell 21

```python
U0 = adata.obsm["cnmf_usage"]
pos = U0[U0 > 0]
eps = max(1e-8, np.quantile(pos, 0.01) * 0.1)  # e.g. 10% of 1% quantile
Uf = np.where(U0 > 0, U0, eps)
Uf = Uf / Uf.sum(axis=1, keepdims=True)
adata.obsm["cnmf_usage"] = Uf
```

## Code cell 22

```python
# The generator stores usage in adata.obsm["usage"].
# Pass usage_key="usage" when preparing inputs.

inputs = prepare_crt_inputs(
    adata=adata,
    usage_key="cnmf_usage",
    covar_key="covar",
    guide_assignment_key="guide_assignment",
    guide2gene_key="guide2gene"
)
```

## Code cell 23

```python
from src.sceptre import (
    build_ntc_group_inputs,
    compute_ntc_group_null_pvals_parallel,
    crt_pvals_for_ntc_groups_ensemble,
    crt_pvals_for_ntc_groups_ensemble_skew,
    make_ntc_groups_ensemble,
    run_all_genes_union_crt,
)
from src.visualization import qq_plot_ntc_pvals
```

## Code cell 24

```python
import os
cpus_per_task = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
```

## Code cell 25

```python
cpus_per_task
```

### Text output (text/plain)

```text
8
```

## Markdown cell 26

### NTC-Null Empirical p-values

## Code cell 27

```python
from src.sceptre import make_ntc_empirical_qq_plots, run_all_genes_union_crt

out = run_all_genes_union_crt(
    inputs=inputs,
    null_method="ntc_empirical",
    qq_crossfit=True,
    null_kwargs=dict(
        ntc_labels=["SAFE", "non-targeting", "NTC"],
        guides_per_unit=6,
        n_ntc_units=5000,
        batch_mode="meta",
        combine_method="fisher",
        matching=dict(n_n1_bins=10, n_d_bins=10, min_ntc_per_bin=50),
        min_treated=10,
        min_control=10,
        qq_crossfit_seed=11,
        qq_crossfit_frac_A=0.5,
    ),
)
```

## Code cell 28

```python
make_ntc_empirical_qq_plots(
    out=out,
    out_dir="results/qq_ntc",
    programs_to_plot="top_var",
    n_programs=6,
    make_per_batch=True,
    make_meta=True,
    envelope_boot=200,
    seed=0,
)
```

## Markdown cell 29

### Stratified CRT

## Code cell 30

```python
resampling_kwargs = dict(
    n_bins=20,
    stratify_by_batch=True,
    batch_key="batch",
    min_stratum_size=2,
    burden_key="log1p_non_ntc_guides_per_cell",
    n_burden_bins=8,
    burden_bin_method="quantile",
    burden_clip_quantiles=(0.0, 1.0),
)

out = run_all_genes_union_crt(
    inputs=inputs,
    B=1023,
    n_jobs=8,
    resampling_method="stratified_perm",
    resampling_kwargs=resampling_kwargs,
    calibrate_skew_normal=True,
    return_raw_pvals=True,
    return_skew_normal=True,
)
```

## Code cell 31

```python
ntc_labels = ["SAFE","non-targeting"]
# Identify NTC guides and build guide-frequency bins / real-gene signatures
ntc_guides, guide_freq, guide_to_bin, real_sigs = build_ntc_group_inputs(
    inputs=inputs,
    ntc_label=ntc_labels,
    group_size=6,
    n_bins=10,
)
```

## Code cell 32

```python

# Create multiple random partitions (ensembles) of NTC guides into 6-guide groups
ntc_groups_ens = make_ntc_groups_ensemble(
    ntc_guides=ntc_guides,
    ntc_freq=guide_freq,
    real_gene_bin_sigs=real_sigs,
    guide_to_bin=guide_to_bin,
    n_ensemble=10,
    seed0=7,
    group_size=6,
    max_groups=None,
)
```

## Code cell 33

```python
ntc_group_pvals_ens = crt_pvals_for_ntc_groups_ensemble(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=1023,
    seed0=23,
    resampling_method="stratified_perm",
    resampling_kwargs=resampling_kwargs,
)
ntc_group_pvals_skew_ens = crt_pvals_for_ntc_groups_ensemble_skew(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=1023,
    seed0=23,
    resampling_method="stratified_perm",
    resampling_kwargs=resampling_kwargs,
)
null_pvals = compute_ntc_group_null_pvals_parallel(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=1023,
    n_jobs=8,
    backend="threading",
    resampling_method="stratified_perm",
    resampling_kwargs=resampling_kwargs,
)
```

## Code cell 34

```python
ax = qq_plot_ntc_pvals(
    pvals_raw_df=out["pvals_raw_df"],
    guide2gene=adata.uns["guide2gene"],
    ntc_genes=ntc_labels,
    pvals_skew_df=out["pvals_df"],
    null_pvals=null_pvals,
    ntc_group_pvals_ens=ntc_group_pvals_ens,
    ntc_group_pvals_skew_ens=ntc_group_pvals_skew_ens,
    show_ntc_ensemble_band=True,
    show_all_pvals=True,
    title="QQ plot: grouped NTC controls (raw vs skew) vs CRT null",
)

import matplotlib.pyplot as plt
plt.show()
```

## Markdown cell 35

### Test with synthetic data

## Code cell 36

```python
from tests.synthetic_data import make_sceptre_style_synth
from src.sceptre import prepare_crt_inputs

adata = make_sceptre_style_synth(
    N=10000,
    K=20,
    n_target_genes=200,
    guides_per_gene=6,
    ntc_frac_guides=0.15,
    moi_mean=5.5,
    frac_causal_genes=0.10,
    n_effect_programs=3,
    effect_size=1.2,
    confound_strength=0.0,
    seed=0,
)

# The generator stores usage in adata.obsm["usage"].
# Pass usage_key="usage" when preparing inputs.
inputs = prepare_crt_inputs(adata, usage_key="usage", add_intercept=False, standardize=False)
```

## Code cell 37

```python
from src.sceptre import (
    build_ntc_group_inputs,
    compute_ntc_group_null_pvals_parallel,
    crt_pvals_for_ntc_groups_ensemble,
    crt_pvals_for_ntc_groups_ensemble_skew,
    make_ntc_groups_ensemble,
    run_all_genes_union_crt,
)
from src.visualization import qq_plot_ntc_pvals

out = run_all_genes_union_crt(
    inputs=inputs,
    B=1023,
    n_jobs=cpus_per_task,
    calibrate_skew_normal=True,
    return_raw_pvals=True,
    return_skew_normal=True,
)

ntc_labels = ["NTC"]
# Identify NTC guides and build guide-frequency bins / real-gene signatures
ntc_guides, guide_freq, guide_to_bin, real_sigs = build_ntc_group_inputs(
    inputs=inputs,
    ntc_label=ntc_labels,
    group_size=6,
    n_bins=10,
)
# Create multiple random partitions (ensembles) of NTC guides into 6-guide groups
ntc_groups_ens = make_ntc_groups_ensemble(
    ntc_guides=ntc_guides,
    ntc_freq=guide_freq,
    real_gene_bin_sigs=real_sigs,
    guide_to_bin=guide_to_bin,
    n_ensemble=10,
    seed0=7,
    group_size=6,
    max_groups=None,
)
# Compute raw CRT p-values for each NTC group in each ensemble
ntc_group_pvals_ens = crt_pvals_for_ntc_groups_ensemble(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=1023,
    seed0=23,
)
# Compute skew-calibrated CRT p-values for each NTC group in each ensemble
ntc_group_pvals_skew_ens = crt_pvals_for_ntc_groups_ensemble_skew(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=1023,
    seed0=23,
)

# Build CRT-null p-values matched to NTC group units (recommended)
null_pvals = compute_ntc_group_null_pvals_parallel(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=1023,
    n_jobs=cpus_per_task,
    backend="threading",
)

ax = qq_plot_ntc_pvals(
    pvals_raw_df=out["pvals_raw_df"],
    guide2gene=adata.uns["guide2gene"],
    ntc_genes=ntc_labels,
    pvals_skew_df=out["pvals_df"],
    null_pvals=null_pvals,
    ntc_group_pvals_ens=ntc_group_pvals_ens,
    ntc_group_pvals_skew_ens=ntc_group_pvals_skew_ens,
    show_ntc_ensemble_band=True,
    show_all_pvals=True,
    title="QQ plot: grouped NTC controls (raw vs skew) vs CRT null",
)

import matplotlib.pyplot as plt
plt.show()
```

## Markdown cell 38

### CRT without stratification

## Code cell 39

```python

out = run_all_genes_union_crt(
    inputs=inputs,
    B=1023,
    n_jobs=cpus_per_task,
    calibrate_skew_normal=True,
    return_raw_pvals=True,
    return_skew_normal=True,
)

```

## Code cell 40

```python
ntc_labels = ["SAFE","non-targeting"]
# Identify NTC guides and build guide-frequency bins / real-gene signatures
ntc_guides, guide_freq, guide_to_bin, real_sigs = build_ntc_group_inputs(
    inputs=inputs,
    ntc_label=ntc_labels,
    group_size=6,
    n_bins=10,
)
```

## Code cell 41

```python

# Create multiple random partitions (ensembles) of NTC guides into 6-guide groups
ntc_groups_ens = make_ntc_groups_ensemble(
    ntc_guides=ntc_guides,
    ntc_freq=guide_freq,
    real_gene_bin_sigs=real_sigs,
    guide_to_bin=guide_to_bin,
    n_ensemble=10,
    seed0=7,
    group_size=6,
    max_groups=None,
)
# Compute raw CRT p-values for each NTC group in each ensemble
ntc_group_pvals_ens = crt_pvals_for_ntc_groups_ensemble(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=1023,
    seed0=23,
)
# Compute skew-calibrated CRT p-values for each NTC group in each ensemble
ntc_group_pvals_skew_ens = crt_pvals_for_ntc_groups_ensemble_skew(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=1023,
    seed0=23,
) 
```

## Code cell 42

```python
# Build CRT-null p-values matched to NTC group units (recommended)
null_pvals = compute_ntc_group_null_pvals_parallel(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=1023,
    n_jobs=cpus_per_task,
    backend="threading",
)

ax = qq_plot_ntc_pvals(
    pvals_raw_df=out["pvals_raw_df"],
    guide2gene=adata.uns["guide2gene"],
    ntc_genes=ntc_labels,
    pvals_skew_df=out["pvals_df"],
    null_pvals=null_pvals,
    ntc_group_pvals_ens=ntc_group_pvals_ens,
    ntc_group_pvals_skew_ens=ntc_group_pvals_skew_ens,
    show_ntc_ensemble_band=True,
    show_all_pvals=True,
    title="QQ plot: grouped NTC controls (raw vs skew) vs CRT null",
)

import matplotlib.pyplot as plt
plt.show()
```

## Code cell 43

```python

```

