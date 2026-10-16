# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/ipsc_ec_sceptre.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/ipsc_ec_sceptre.ipynb`
- SHA-256: `b0d8e4bf45437928d39442e44256d82357467efd74bf7599eefce9b26de3008b`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Code cell 0

```python
%reset
```

### Text output (stream)

```text
Nothing done.

```

## Code cell 1

```python
import sys
print(sys.executable)
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/bin/python

```

## Code cell 2

```python
import scanpy as sc
import muon as mu
import pandas as pd
import numpy as np
```

### Text output (stream)

```text
Matplotlib is building the font cache; this may take a moment.
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/tqdm/auto.py:21: TqdmWarning: IProgress not found. Please update jupyter and ipywidgets. See https://ipywidgets.readthedocs.io/en/stable/user_install.html
  from .autonotebook import tqdm as notebook_tqdm
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/muon/_core/preproc.py:31: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
  if Version(scanpy.__version__) < Version("1.10"):

```

## Code cell 3

```python
data_dir = "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data"
working_dir = "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks"

```

## Markdown cell 4

### Test

## Code cell 5

```python
d0_mdata = mu.read(f"{data_dir}/d0_normalized_gex_and_guide.h5mu")
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

## Code cell 6

```python
d0_mdata
```

### Text output (text/plain)

```text
MuData object with n_obs × n_vars = 38609 × 18318
  2 modalities
    GEX:	38609 x 16392
      obs:	'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0'
      var:	'gene_id', 'symbol', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'n_counts', 'highly_deviant', 'highly_deviant_2k', 'binomial_deviance', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'orig_highly_variable', 'highly_variable_2k'
      uns:	'hvg', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'neighbors', 'pca', 'umap'
      obsm:	'X_pca', 'X_umap', 'guide_assignment', 'guide_counts'
      varm:	'PCs'
      layers:	'ambiguous', 'analytic_pearson_residuals', 'log1p_norm', 'mature', 'nascent', 'raw_counts', 'scran_normalization'
      obsp:	'connectivities', 'distances'
    guide:	38609 x 1926
      obs:	'bci', 'stype', 'well', 'day', 'cell_pool'
```

## Code cell 7

```python
d0_mdata['GEX']
```

### Text output (text/plain)

```text
AnnData object with n_obs × n_vars = 38609 × 16392
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0'
    var: 'gene_id', 'symbol', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'n_counts', 'highly_deviant', 'highly_deviant_2k', 'binomial_deviance', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'orig_highly_variable', 'highly_variable_2k'
    uns: 'hvg', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'neighbors', 'pca', 'umap'
    obsm: 'X_pca', 'X_umap', 'guide_assignment', 'guide_counts'
    varm: 'PCs'
    layers: 'ambiguous', 'analytic_pearson_residuals', 'log1p_norm', 'mature', 'nascent', 'raw_counts', 'scran_normalization'
    obsp: 'connectivities', 'distances'
```

## Code cell 8

```python
d0_mdata['GEX'].layers["raw_counts"]
```

### Text output (text/plain)

```text
<Compressed Sparse Row sparse matrix of dtype 'float32'
	with 65423017 stored elements and shape (38609, 16392)>
```

## Code cell 9

```python
d0_mdata['GEX'].var_names
```

### Text output (text/plain)

```text
Index(['ENSG00000227232.5', 'ENSG00000238009.6', 'ENSG00000241860.7',
       'ENSG00000279457.4', 'ENSG00000228463.10', 'ENSG00000290385.1',
       'ENSG00000230021.10', 'ENSG00000225972.1', 'ENSG00000225630.1',
       'ENSG00000237973.1',
       ...
       'ENSG00000210196.2', 'ENSG00000271254.7', 'ENSG00000276345.1',
       'ENSG00000276700.1', 'ENSG00000275757.1', 'ENSG00000277196.4',
       'ENSG00000276256.1', 'ENSG00000273748.1', 'ENSG00000274917.1',
       'ENSG00000273730.1'],
      dtype='object', name='gene_id', length=16392)
```

## Code cell 10

```python
d0_mdata['guide'].var_names
```

### Text output (text/plain)

```text
Index(['ACAA1_+_38178488.23-P1P2', 'ACAA1_+_38178517.23-P1P2',
       'ACAA1_+_38178559.23-P1P2', 'ACAA1_+_38178570.23-P1P2',
       'ACAA1_+_38178620.23-P1P2', 'ACAA1_-_38178575.23-P1P2',
       'ACTL7B_+_111618040.23-P1P2', 'ACTL7B_+_111618198.23-P1P2',
       'ACTL7B_+_111618251.23-P1P2', 'ACTL7B_-_111617995.23-P1P2',
       ...
       'SAFE_TARGETING_933', 'SAFE_TARGETING_986', 'WDR61_SGTN_13',
       'WDR61_SGTN_14', 'SUZ12_SGTN_22', 'SUZ12_SGTN_23', 'MESDC1_SGTN_11',
       'MESDC1_SGTN_12', 'TP53_SGTN_01', 'TP53_SGTN_02'],
      dtype='object', name='feature_id', length=1926)
```

## Code cell 11

```python
d0_mdata['guide']
```

### Text output (text/plain)

```text
AnnData object with n_obs × n_vars = 38609 × 1926
    obs: 'bci', 'stype', 'well', 'day', 'cell_pool'
```

## Code cell 12

```python
d0_mdata['guide'].obs["stype"]
```

### Text output (text/plain)

```text
barcode
CCTAATCCAACAACCATCTATTAC_2    T
GGTGCGAACCTCCTGAATCCGCGA_2    T
AACTCACCATCATTCCCTGTCCCG_1    T
TTCACGCACACTTCGAATTCATGG_2    T
GTCTGTCAGAGCTGAACCGTTCTA_3    T
                             ..
CGAACTTAGAGTTAGCATTCATGG_3    T
CCGTGAGAGAATCTGAATCGCATA_3    T
TGAAGAGACCTCTATCTCTATTAC_3    T
CCATCCTCGAGCTGAAAAACGATA_2    T
TCCGTCTAAACCGAGACGTTCGAG_2    T
Name: stype, Length: 38609, dtype: category
Categories (1, object): ['T']
```

## Code cell 13

```python
d0_mdata['GEX'].var['symbol']
```

### Text output (text/plain)

```text
gene_id
ENSG00000227232.5              WASH7P
ENSG00000238009.6     ENSG00000238009
ENSG00000241860.7     ENSG00000241860
ENSG00000279457.4              WASH9P
ENSG00000228463.10    ENSG00000228463
                           ...       
ENSG00000277196.4               PRODH
ENSG00000276256.1     ENSG00000276256
ENSG00000273748.1     ENSG00000273748
ENSG00000274917.1           RNA5-8SN5
ENSG00000273730.1           5_8S_rRNA
Name: symbol, Length: 16392, dtype: category
Categories (16338, object): ['5_8S_rRNA', '7SK', 'A1BG', 'A1BG-AS1', ..., 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']
```

## Code cell 14

```python
guide_meta_df = pd.read_csv(f"{working_dir}/300genes_guide_metadata_v43.tsv", sep='\t')
```

## Code cell 15

```python
guide_meta_df
```

### Text output (text/plain)

```text
                      guide_id                spacer  targeting       type  \
0     ACAA1_+_38178488.23-P1P2   TAGTGCTGGGCCACCTGAG       True  targeting   
1     ACAA1_+_38178517.23-P1P2   GTTGACCTGCGCAATGCAG       True  targeting   
2     ACAA1_+_38178559.23-P1P2   TAACTCCGCGGTCAGTTCC       True  targeting   
3     ACAA1_+_38178570.23-P1P2   CGGCTGTCTGTTAACTCCG       True  targeting   
4     ACAA1_+_38178620.23-P1P2   GGCGGAGTTGAGGCCTTGG       True  targeting   
...                        ...                   ...        ...        ...   
1921             SUZ12_SGTN_23  TTGCAGGCGCTTGCTCTCCG       True  targeting   
1922            MESDC1_SGTN_11  AGAAACTGAAAGGAAGGGCG       True  targeting   
1923            MESDC1_SGTN_12   CCGCGCGCGGAGAGGAAGC       True  targeting   
1924              TP53_SGTN_01   GGACACTTTGCGTTCGGGC       True  targeting   
1925              TP53_SGTN_02  AGCGCCAGTCTTGAGCACAT       True  targeting   

     guide_chr  guide_start   guide_end strand  pam genomic_element  \
0         chr3   38137000.0  38137019.0      -  NGG        promoter   
1         chr3   38137029.0  38137048.0      -  NGG        promoter   
2         chr3   38137071.0  38137090.0      -  NGG        promoter   
3         chr3   38137082.0  38137101.0      -  NGG        promoter   
4         chr3   38137132.0  38137151.0      -  NGG        promoter   
...        ...          ...         ...    ...  ...             ...   
1921     chr17   31937180.0  31937203.0      +  NGG        promoter   
1922     chr15   81001111.0  81001134.0      -  NGG        promoter   
1923     chr15   81001066.0  81001088.0      +  NGG        promoter   
1924     chr17    7687417.0   7687439.0      -  NGG        promoter   
1925     chr17    7687505.0   7687528.0      +  NGG        promoter   

     intended_target_name       ensembl_id intended_target_chr  \
0                   ACAA1  ENSG00000060971                chr3   
1                   ACAA1  ENSG00000060971                chr3   
2                   ACAA1  ENSG00000060971                chr3   
3                   ACAA1  ENSG00000060971                chr3   
4                   ACAA1  ENSG00000060971                chr3   
...                   ...              ...                 ...   
1921                SUZ12  ENSG00000178691               chr17   
1922               TLNRD1  ENSG00000140406               chr15   
1923               TLNRD1  ENSG00000140406               chr15   
1924                 TP53  ENSG00000141510               chr17   
1925                 TP53  ENSG00000141510               chr17   

      intended_target_start  intended_target_end putative_target_genes  \
0                38103129.0           38137242.0                 ACAA1   
1                38103129.0           38137242.0                 ACAA1   
2                38103129.0           38137242.0                 ACAA1   
3                38103129.0           38137242.0                 ACAA1   
4                38103129.0           38137242.0                 ACAA1   
...                     ...                  ...                   ...   
1921             31937007.0           32001038.0                 SUZ12   
1922             81000923.0           81005788.0                TLNRD1   
1923             81000923.0           81005788.0                TLNRD1   
1924              7661779.0            7687538.0                  TP53   
1925              7661779.0            7687538.0                  TP53   

      reporter  imperfect  
0          NaN        NaN  
1          NaN        NaN  
2          NaN        NaN  
3          NaN        NaN  
4          NaN        NaN  
...        ...        ...  
1921       NaN        NaN  
1922       NaN        NaN  
1923       NaN        NaN  
1924       NaN        NaN  
1925       NaN        NaN  

[1926 rows x 18 columns]
```

## Code cell 16

```python
import scipy.io as sio
import pandas as pd
import os
```

## Markdown cell 17

### Day 0 

## Code cell 18

```python
mdata = d0_mdata

# Filter out cells that have no guides
#mdata = mdata[mdata['GEX'].obs['guides_per_cell'] > 0]
```

## Code cell 19

```python
# Define a temporary folder for transfer
temp_dir = working_dir + "/d0_sceptre_transfer"
os.makedirs(temp_dir, exist_ok=True)

# 1. EXTRACT MATRICES (Transpose them!)
# Scanpy is (Cells × Genes), Sceptre needs (Genes × Cells)
print("Writing sparse matrices...")
sio.mmwrite(f"{temp_dir}/response_matrix.mtx", mdata['GEX'].layers['raw_counts'].T)
sio.mmwrite(f"{temp_dir}/grna_matrix.mtx", mdata['guide'].X.T)

# 2. EXTRACT NAMES (Row/Col names are lost in .mtx)
# Save Gene names and Guide names to assign them later in R
pd.Series(mdata['GEX'].var_names).to_csv(f"{temp_dir}/gene_names.csv", index=False, header=False)
pd.Series(mdata['guide'].var_names).to_csv(f"{temp_dir}/guide_names.csv", index=False, header=False)
pd.Series(mdata['GEX'].obs_names).to_csv(f"{temp_dir}/cell_names.csv", index=False, header=False)
```

### Text output (stream)

```text
Writing sparse matrices...

```

## Code cell 20

```python
mdata['GEX']
```

### Text output (text/plain)

```text
AnnData object with n_obs × n_vars = 38609 × 16392
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0'
    var: 'gene_id', 'symbol', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'n_counts', 'highly_deviant', 'highly_deviant_2k', 'binomial_deviance', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'orig_highly_variable', 'highly_variable_2k'
    uns: 'hvg', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'neighbors', 'pca', 'umap'
    obsm: 'X_pca', 'X_umap', 'guide_assignment', 'guide_counts'
    varm: 'PCs'
    layers: 'ambiguous', 'analytic_pearson_residuals', 'log1p_norm', 'mature', 'nascent', 'raw_counts', 'scran_normalization'
    obsp: 'connectivities', 'distances'
```

## Code cell 21

```python
# 3. EXTRACT COVARIATES 
rna = mdata['GEX']
covariates_df = pd.DataFrame({
    "batch": rna.obs['sample_id'],
    "grna_n_nonzero": rna.obs['guides_per_cell'],
    "grna_n_umis": rna.obs['guide_umi_counts'],
    "response_n_nonzero": rna.obs['n_genes_by_counts'],
    "response_n_umis": rna.obs['total_counts'],
    "response_p_mito": rna.obs['pct_counts_mt'],
})
covariates_df.to_csv(f"{temp_dir}/covariates.csv")
covariates_df
```

### Text output (text/plain)

```text
                                             batch  grna_n_nonzero  \
barcode                                                              
CCTAATCCAACAACCATCTATTAC_2  300genes:eTN5_GEX_sub2               5   
GGTGCGAACCTCCTGAATCCGCGA_2  300genes:eTN5_GEX_sub2               5   
AACTCACCATCATTCCCTGTCCCG_1  300genes:eTN5_GEX_sub1               4   
TTCACGCACACTTCGAATTCATGG_2  300genes:eTN5_GEX_sub2               0   
GTCTGTCAGAGCTGAACCGTTCTA_3  300genes:eTN5_GEX_sub3               3   
...                                            ...             ...   
CGAACTTAGAGTTAGCATTCATGG_3  300genes:eTN5_GEX_sub3              12   
CCGTGAGAGAATCTGAATCGCATA_3  300genes:eTN5_GEX_sub3               2   
TGAAGAGACCTCTATCTCTATTAC_3  300genes:eTN5_GEX_sub3               7   
CCATCCTCGAGCTGAAAAACGATA_2  300genes:eTN5_GEX_sub2               4   
TCCGTCTAAACCGAGACGTTCGAG_2  300genes:eTN5_GEX_sub2               5   

                            grna_n_umis  response_n_nonzero  response_n_umis  \
barcode                                                                        
CCTAATCCAACAACCATCTATTAC_2        195.0                2049      3969.242920   
GGTGCGAACCTCCTGAATCCGCGA_2        176.0                2232      3808.424805   
AACTCACCATCATTCCCTGTCCCG_1        183.0                1656      3441.276367   
TTCACGCACACTTCGAATTCATGG_2          3.0                 688       805.500061   
GTCTGTCAGAGCTGAACCGTTCTA_3         46.0                2623      5361.565430   
...                                 ...                 ...              ...   
CGAACTTAGAGTTAGCATTCATGG_3        233.0                1729      4035.552734   
CCGTGAGAGAATCTGAATCGCATA_3         57.0                1205      1955.549194   
TGAAGAGACCTCTATCTCTATTAC_3        344.0                2099      3987.371338   
CCATCCTCGAGCTGAAAAACGATA_2        256.0                1923      3728.120605   
TCCGTCTAAACCGAGACGTTCGAG_2         66.0                 924      1596.008423   

                            response_p_mito  
barcode                                      
CCTAATCCAACAACCATCTATTAC_2        14.708926  
GGTGCGAACCTCCTGAATCCGCGA_2        10.594932  
AACTCACCATCATTCCCTGTCCCG_1         9.923644  
TTCACGCACACTTCGAATTCATGG_2         1.738051  
GTCTGTCAGAGCTGAACCGTTCTA_3        12.086022  
...                                     ...  
CGAACTTAGAGTTAGCATTCATGG_3        10.494226  
CCGTGAGAGAATCTGAATCGCATA_3        13.985841  
TGAAGAGACCTCTATCTCTATTAC_3         9.166440  
CCATCCTCGAGCTGAAAAACGATA_2        14.229690  
TCCGTCTAAACCGAGACGTTCGAG_2        18.796894  

[38609 rows x 6 columns]
```

## Code cell 22

```python
# 4. PREPARE GUIDE MAPPING 
# Sceptre needs a table mapping 'grna_id' to 'grna_target'.

# 1. Create a "skeleton" DataFrame from the actual matrix order
# This ensures the order matches mdata['guide'].X perfectly
# We name the index 'grna_id' to match Sceptre's expectation immediately
guide_mapping = pd.DataFrame(index=mdata['guide'].var_names)
guide_mapping.index.name = 'grna_id'

# 2. Prepare metadata for merging
# Ensure 'guide_id' is the index so we can map efficiently
lookup_df = guide_meta_df.set_index('guide_id')

# 3. Map the targets onto the skeleton
# We use .loc[] to pull the specific target for each guide in the exact order required
guide_mapping['grna_target'] = lookup_df.loc[guide_mapping.index, 'intended_target_name']

# 4. Reset index so 'grna_id' becomes a column (Sceptre requires this)
guide_mapping = guide_mapping.reset_index()

# ==========================================
# VALIDATION CHECKS 
# ==========================================

# Check 1: Are there any NAs? (Guides in mdata that were missing from meta_df)
na_count = guide_mapping['grna_target'].isna().sum()
if na_count > 0:
    print(f"⚠️ WARNING: {na_count} guides have NaN targets! These will crash Sceptre.")
else:
    print("✅ All guides mapped successfully.")

# Check 2: Verify the first few rows manually
print("\nFirst 5 rows of mapping (Must match mdata['guide'].var_names order):")
print(guide_mapping.head())

# 5. Save for Sceptre
guide_mapping.to_csv(f"{temp_dir}/guide_mapping.csv", index=False)
```

### Text output (stream)

```text
✅ All guides mapped successfully.

First 5 rows of mapping (Must match mdata['guide'].var_names order):
                    grna_id grna_target
0  ACAA1_+_38178488.23-P1P2       ACAA1
1  ACAA1_+_38178517.23-P1P2       ACAA1
2  ACAA1_+_38178559.23-P1P2       ACAA1
3  ACAA1_+_38178570.23-P1P2       ACAA1
4  ACAA1_+_38178620.23-P1P2       ACAA1

```

## Code cell 23

```python
guide_mapping
```

### Text output (text/plain)

```text
                       grna_id grna_target
0     ACAA1_+_38178488.23-P1P2       ACAA1
1     ACAA1_+_38178517.23-P1P2       ACAA1
2     ACAA1_+_38178559.23-P1P2       ACAA1
3     ACAA1_+_38178570.23-P1P2       ACAA1
4     ACAA1_+_38178620.23-P1P2       ACAA1
...                        ...         ...
1921             SUZ12_SGTN_23       SUZ12
1922            MESDC1_SGTN_11      TLNRD1
1923            MESDC1_SGTN_12      TLNRD1
1924              TP53_SGTN_01        TP53
1925              TP53_SGTN_02        TP53

[1926 rows x 2 columns]
```

## Code cell 24

```python
# Create a mapping from the Gene Symbol (which matches your guides)
# to the Gene Index (which matches the matrix rownames)

# Check which column in mdata['GEX'].var holds the symbol
symbol_col = 'symbol' 

gene_map = mdata['GEX'].var[[symbol_col]].reset_index()
gene_map.columns = ['response_id', 'grna_target'] # Rename for easy merging in R

# Save it
gene_map.to_csv(f"{temp_dir}/gene_symbol_map.csv", index=False)
print("✅ Gene mapping saved.")
```

### Text output (stream)

```text
✅ Gene mapping saved.

```

## Code cell 25

```python
gene_map
```

### Text output (text/plain)

```text
              response_id      grna_target
0       ENSG00000227232.5           WASH7P
1       ENSG00000238009.6  ENSG00000238009
2       ENSG00000241860.7  ENSG00000241860
3       ENSG00000279457.4           WASH9P
4      ENSG00000228463.10  ENSG00000228463
...                   ...              ...
16387   ENSG00000277196.4            PRODH
16388   ENSG00000276256.1  ENSG00000276256
16389   ENSG00000273748.1  ENSG00000273748
16390   ENSG00000274917.1        RNA5-8SN5
16391   ENSG00000273730.1        5_8S_rRNA

[16392 rows x 2 columns]
```

## Markdown cell 26

#### Import data to SCEPTRE (Rcell)

## Code cell 27

```python
%load_ext rpy2.ipython
```

## Code cell 28

```python
%%R
library(sceptre)
library(Matrix)

# Define path
data_dir <- "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d0_sceptre_transfer"

# 1. LOAD MATRICES
cat("Loading matrices into R...\n")
response_mat <- readMM(file.path(data_dir, "response_matrix.mtx"))
grna_mat <- readMM(file.path(data_dir, "grna_matrix.mtx"))

# 2. ATTACH NAMES
gene_names <- read.csv(file.path(data_dir, "gene_names.csv"), header=FALSE)$V1
guide_names <- read.csv(file.path(data_dir, "guide_names.csv"), header=FALSE)$V1
cell_names <- read.csv(file.path(data_dir, "cell_names.csv"), header=FALSE)$V1

rownames(response_mat) <- gene_names
colnames(response_mat) <- cell_names

rownames(grna_mat) <- guide_names
colnames(grna_mat) <- cell_names
```

### Text output (stream)

```text
Loading matrices into R...

```

## Code cell 29

```python
%%R

# 3. LOAD METADATA
covariates <- read.csv(file.path(data_dir, "covariates.csv"), row.names=1)
guide_mapping <- read.csv(file.path(data_dir, "guide_mapping.csv"))

# ==============================================================================
# FIX: RENAME 'response_p_mito' SO WE KEEP IT
# Sceptre blocks the name 'response_p_mito' but doesn't calculate it.
# We rename it to 'mt_percent' so we can use it as a covariate later.
# ==============================================================================

# 1. Rename response_p_mito if it exists
if ("response_p_mito" %in% colnames(covariates)) {
    cat("Renaming 'response_p_mito' to 'mt_percent'...\n")
    names(covariates)[names(covariates) == "response_p_mito"] <- "mt_percent"
}

# 2. DROP OTHER RESERVED COLUMNS (which Sceptre DOES calculate)
# We still drop n_umis and n_nonzero because Sceptre re-calculates those automatically.
reserved_cols <- c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis")
covariates <- covariates[, !(colnames(covariates) %in% reserved_cols), drop=FALSE]

# 4. CREATE SCEPTRE OBJECT
sceptre_object <- import_data(
  response_matrix = response_mat,
  grna_matrix = grna_mat,
  grna_target_data_frame = guide_mapping,
  extra_covariates = covariates,
  moi = "high"
)

print("✅ Sceptre object created!")
print(head(sceptre_object@covariate_data_frame))
```

### Text output (stream)

```text
Renaming 'response_p_mito' to 'mt_percent'...
[1] "\u2705 Sceptre object created!"
                           response_n_nonzero response_n_umis grna_n_nonzero
CCTAATCCAACAACCATCTATTAC_2               1774        3969.243             13
GGTGCGAACCTCCTGAATCCGCGA_2               1785        3808.425             23
AACTCACCATCATTCCCTGTCCCG_1               1439        3441.277             38
TTCACGCACACTTCGAATTCATGG_2                598         805.500              3
GTCTGTCAGAGCTGAACCGTTCTA_3               2308        5361.566             10
CTGGCATAAACCGAGATGTCTGAA_3               1068        2485.453            137
                           grna_n_umis                  batch mt_percent
CCTAATCCAACAACCATCTATTAC_2         195 300genes:eTN5_GEX_sub2  14.708926
GGTGCGAACCTCCTGAATCCGCGA_2         176 300genes:eTN5_GEX_sub2  10.594932
AACTCACCATCATTCCCTGTCCCG_1         183 300genes:eTN5_GEX_sub1   9.923644
TTCACGCACACTTCGAATTCATGG_2           3 300genes:eTN5_GEX_sub2   1.738051
GTCTGTCAGAGCTGAACCGTTCTA_3          46 300genes:eTN5_GEX_sub3  12.086022
CTGGCATAAACCGAGATGTCTGAA_3        1699 300genes:eTN5_GEX_sub3  13.860650

```

## Code cell 30

```python
%%R
# 1. Load the translation map
gene_map <- read.csv(file.path(data_dir, "gene_symbol_map.csv"))

# 2. Start with guide mapping (which has 'grna_target' as Symbol)
discovery_pairs_highmoi <- guide_mapping

# 3. Merge to attach the correct response_id (Gene ID)
#    This matches 'grna_target' in guides to 'grna_target' in the gene map
discovery_pairs_highmoi <- merge(discovery_pairs_highmoi, gene_map, by="grna_target")

# 4. Select ONLY the two columns Sceptre needs
#    (We filter out non-targeting guides automatically via the merge)
discovery_pairs_highmoi <- discovery_pairs_highmoi[, c("grna_target", "response_id")]

# 5. DEDUPLICATE
#    Sceptre performs analysis at the *Target* level, not the *Guide* level.
#    We only need one row per Target-Gene pair.
discovery_pairs_highmoi <- unique(discovery_pairs_highmoi)

print(paste("Unique pairs to test:", nrow(discovery_pairs_highmoi)))
print(head(discovery_pairs_highmoi))
```

### Text output (stream)

```text
[1] "Unique pairs to test: 281"
   grna_target        response_id
1        ACAA1 ENSG00000060971.19
7       ACVR1B ENSG00000135503.13
13      ACVR2A ENSG00000121989.15
19      ACVR2B ENSG00000114739.14
25       ADCY2 ENSG00000078295.17
31      AMOTL2 ENSG00000114019.15

```

## Code cell 31

```python
%%R 

sceptre_object_highmoi <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs_highmoi,
  side = "both", control_group = "complement")
```

## Code cell 32

```python
%%R

png(file.path(data_dir, "grna_count_distributions.png"), width = 3000, height = 2000, res = 300)

print(plot_grna_count_distributions(
  sceptre_object = sceptre_object_highmoi,
  n_grnas_to_plot = 8
))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "grna_count_distributions.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d0_sceptre_transfer/grna_count_distributions.png

```

## Code cell 33

```python
%%R

plot_covariates(sceptre_object_highmoi, p_mito_threshold = 0.25)
```

### Text output (text/plain)

```text
<IPython.core.display.Image object>
```

## Code cell 34

```python
%%R
sceptre_object_highmoi <- assign_grnas(
  sceptre_object = sceptre_object_highmoi,
  method = "mixture",
  parallel = FALSE
)
```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178488.23-P1P2 (1 of 1926)
Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178620.23-P1P2 (5 of 1926)
Performing gRNA-to-cell assignments for gRNA ACTL7B_-_111617995.23-P1P2 (10 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR1B_+_52345762.23-P1P2 (15 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2A_+_148602343.23-P1P2 (20 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_+_38495339.23-P1P2 (25 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_-_38495819.23-P1P2 (30 of 1926)
Performing gRNA-to-cell assignments for gRNA ADCY2_-_7396429.23-P1P2 (35 of 1926)
Performing gRNA-to-cell assignments for gRNA AMOTL2_+_134093279.23-P1P2 (40 of 1926)
Performing gRNA-to-cell assignments for gRNA ANP32E_+_150208342.23-P1P2 (45 of 1926)
Performing gRNA-to-cell assignments for gRNA APAF1_-_99039071.23-P1P2 (50 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_+_72766725.23-P1P2 (55 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_-_72766755.23-P1P2 (60 of 1926)
Performing gRNA-to-cell assignments for gRNA ARTN_-_44399429.23-P1 (65 of 1926)
Performing gRNA-to-cell assignments for gRNA ASCC3_-_101329193.23-P1P2 (70 of 1926)
Performing gRNA-to-cell assignments for gRNA ASH2L_-_37963059.23-P1P2 (75 of 1926)
Performing gRNA-to-cell assignments for gRNA ATMIN_+_81069661.23-P1P2 (80 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_+_111992112.23-P1P2 (85 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_-_111992129.23-P1P2 (90 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP6AP1_-_153657075.23-P1P2 (95 of 1926)
Performing gRNA-to-cell assignments for gRNA ATXN2L_-_28834668.23-P1P2 (100 of 1926)
Performing gRNA-to-cell assignments for gRNA BAD_+_64052102.23-P1P2 (105 of 1926)
Performing gRNA-to-cell assignments for gRNA BAG6_+_31620398.23-P1P2 (110 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_+_49458268.23-P1P2 (115 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_-_49458323.23-P1P2 (120 of 1926)
Performing gRNA-to-cell assignments for gRNA BCAS2_-_115124045.23-P1P2 (125 of 1926)
Performing gRNA-to-cell assignments for gRNA BMPR1A_-_88516386.23-P1P2 (130 of 1926)
Performing gRNA-to-cell assignments for gRNA BRD3_-_136932911.23-P1P2 (135 of 1926)
Performing gRNA-to-cell assignments for gRNA C17orf49_+_6918114.23-P2 (140 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_+_31532436.23-P1P2 (145 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_-_31532650.23-P1P2 (150 of 1926)
Performing gRNA-to-cell assignments for gRNA CACFD1_-_136325160.23-P1P2 (155 of 1926)

```

### Text output (stream)

```text
Performing gRNA-to-cell assignments for gRNA CACNG1_+_65040756.23-P1P2 (160 of 1926)
Performing gRNA-to-cell assignments for gRNA CAPN5_+_76778225.23-P1P2 (165 of 1926)
Performing gRNA-to-cell assignments for gRNA CASP3_+_185570560.23-P1P2 (170 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_+_107384793.23-P1P2 (175 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_-_107384635.23-P1P2 (180 of 1926)
Performing gRNA-to-cell assignments for gRNA CCM2_-_45039941.23-P1P2 (185 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNC_-_100016267.23-P1P2 (190 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNT2_-_135676416.23-P1P2 (195 of 1926)
Performing gRNA-to-cell assignments for gRNA CCT6B_+_33288381.23-P1P2 (200 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_+_30366413.23-P1P2 (205 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_-_30366630.23-P1P2 (210 of 1926)
Performing gRNA-to-cell assignments for gRNA CD99L2_-_150067120.23-P1P2 (215 of 1926)
Performing gRNA-to-cell assignments for gRNA CDH5_-_66400630.23-P1P2 (220 of 1926)
Performing gRNA-to-cell assignments for gRNA CDKN1A_-_36646538.23-P1P2 (225 of 1926)
Performing gRNA-to-cell assignments for gRNA CENPC_+_68411190.23-P1P2 (230 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_+_30773828.23-P1P2 (235 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_-_30773391.23-P1P2 (240 of 1926)
Performing gRNA-to-cell assignments for gRNA CHMP3_-_86790451.23-P1P2 (245 of 1926)
Performing gRNA-to-cell assignments for gRNA CHORDC1_-_89955894.23-P1P2 (250 of 1926)
Performing gRNA-to-cell assignments for gRNA CNEP1R1_+_50059284.23-P1P2 (255 of 1926)
Performing gRNA-to-cell assignments for gRNA CNOT2_+_70637651.23-P1P2 (260 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_+_83956339.23-P1P2 (265 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_-_83956572.23-P1P2 (270 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS6_-_99686676.23-P1P2 (275 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_-_31634088.23-P1P2 (280 of 1926)
Performing gRNA-to-cell assignments for gRNA CTR9_+_10772881.23-P1P2 (285 of 1926)
Performing gRNA-to-cell assignments for gRNA CUL3_-_225449667.23-P1P2 (290 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_+_145150267.23-P1P2 (295 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_-_145149996.23-P1P2 (300 of 1926)
Performing gRNA-to-cell assignments for gRNA CYCS_+_25164869.23-P1P2 (305 of 1926)
Performing gRNA-to-cell assignments for gRNA DACH1_-_72441101.23-P1P2 (310 of 1926)
Performing gRNA-to-cell assignments for gRNA DBR1_+_137893673.23-P1P2 (315 of 1926)
Performing gRNA-to-cell assignments for gRNA DCTN2_-_57940910.23-P1P2 (320 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_+_61851532.23-P1P2 (325 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_-_61851812.23-P1P2 (330 of 1926)
Performing gRNA-to-cell assignments for gRNA DGCR8_-_20068158.23-P1P2 (335 of 1926)
Performing gRNA-to-cell assignments for gRNA DNAJC17_-_41099215.23-P1P2 (340 of 1926)
Performing gRNA-to-cell assignments for gRNA DNM2_+_10829157.23-P1P2 (345 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM1_+_49575018.23-P1P2 (350 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_+_130700056.23-P1P2 (355 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_-_130700076.23-P1P2 (360 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM3_+_155113020.23-P1P2 (365 of 1926)
Performing gRNA-to-cell assignments for gRNA DROSHA_-_31531982.23-P1P2 (370 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_+_3985081.23-P1 (375 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_-_3985443.23-P1 (380 of 1926)
Performing gRNA-to-cell assignments for gRNA EHMT2_-_31865199.23-P1P2 (385 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3F_+_8009347.23-P1P2 (390 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3M_-_32605405.23-P1P2 (395 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4E_+_99849831.23-P1P2 (400 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184032364.23-P1P2 (405 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184033042.23-P1P2 (410 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G2_+_10830447.23-P1 (415 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF5A_-_7211217.23-P1P2 (420 of 1926)
Performing gRNA-to-cell assignments for gRNA TCEB3_-_24069975.23-P1P2 (425 of 1926)
Performing gRNA-to-cell assignments for gRNA EMC2_+_109456080.23-P1P2 (430 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823220.23-P1P2 (435 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823252.23-P1P2 (440 of 1926)
Performing gRNA-to-cell assignments for gRNA EOMES_-_27763342.23-P1P2 (445 of 1926)
Performing gRNA-to-cell assignments for gRNA ERCC2_+_45873802.23-P1P2 (450 of 1926)
Performing gRNA-to-cell assignments for gRNA ERO1L_+_53162425.23-P1P2 (455 of 1926)
Performing gRNA-to-cell assignments for gRNA ETS1_+_128392044.23-P1P2 (460 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_+_36134592.23-P1P2 (465 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_-_36134789.23-P1P2 (470 of 1926)
Performing gRNA-to-cell assignments for gRNA EWSR1_-_29664347.23-P1P2 (475 of 1926)
Performing gRNA-to-cell assignments for gRNA EXOSC8_-_37574979.23-P1P2 (480 of 1926)
Performing gRNA-to-cell assignments for gRNA FBXW11_-_171433300.23-P1P2 (485 of 1926)
Performing gRNA-to-cell assignments for gRNA FEN1_+_61560782.23-P1P2 (490 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_+_75532080.23-P1P2 (495 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_-_75532461.23-P1P2 (500 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT7_-_139926495.23-P1P2 (505 of 1926)
Performing gRNA-to-cell assignments for gRNA G3BP1_-_151151499.23-P1P2 (510 of 1926)
Performing gRNA-to-cell assignments for gRNA GCK_+_44229021.23-P1P2 (515 of 1926)
Performing gRNA-to-cell assignments for gRNA GOLT1B_+_21654777.23-P1P2 (520 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_+_45000529.23-P1P2 (525 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_-_45000948.23-P1P2 (530 of 1926)
Performing gRNA-to-cell assignments for gRNA GPI_+_34856284.23-P1P2 (535 of 1926)
Performing gRNA-to-cell assignments for gRNA GPR107_-_132816044.23-P1P2 (540 of 1926)
Performing gRNA-to-cell assignments for gRNA GSPT1_-_12009714.23-P1P2 (545 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A1_+_81687215.23-P1P2 (550 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_+_59949355.23-P1P2 (555 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_-_59949657.23-P1P2 (560 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2H4_-_30876185.23-P1P2 (565 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF3C2_+_27579604.23-P1P2 (570 of 1926)
Performing gRNA-to-cell assignments for gRNA HAND1_+_153857732.23-P1P2 (575 of 1926)
Performing gRNA-to-cell assignments for gRNA HDAC2_-_114292253.23-P1P2 (580 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_+_81399312.23-P1P2 (585 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_-_81399351.23-P1P2 (590 of 1926)
Performing gRNA-to-cell assignments for gRNA HINT2_+_35814995.23-P1P2 (595 of 1926)
Performing gRNA-to-cell assignments for gRNA HNF4A_+_43029971.23-P1P2 (600 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPM_-_8509885.23-P1 (605 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPU_+_245027753.23-P1P2 (610 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_+_43702292.23-P1P2 (615 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_-_43702305.23-P1P2 (620 of 1926)
Performing gRNA-to-cell assignments for gRNA HSPA5_-_128003534.23-P1P2 (625 of 1926)
Performing gRNA-to-cell assignments for gRNA HUWE1_-_53710540.23-P1P2 (630 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_-_30007588.23-P1P2 (635 of 1926)
Performing gRNA-to-cell assignments for gRNA INTS10_+_19675000.23-P1P2 (640 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_+_65903046.23-P1P2 (645 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_-_65903523.23-P1P2 (650 of 1926)
Performing gRNA-to-cell assignments for gRNA C7orf26_+_6630358.23-P1P2 (655 of 1926)
Performing gRNA-to-cell assignments for gRNA IPO8_+_30848771.23-P1P2 (660 of 1926)
Performing gRNA-to-cell assignments for gRNA ITGB1BP1_+_9563543.23-P1P2 (665 of 1926)
Performing gRNA-to-cell assignments for gRNA SSFA2_+_182756743.23-P1P2 (670 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_+_109914946.23-P1P2 (675 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_-_109915112.23-P1P2 (680 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM1A_-_23346409.23-P1P2 (685 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM3B_-_137688431.23-P1P2 (690 of 1926)
Performing gRNA-to-cell assignments for gRNA KDR_-_55991659.23-P1P2 (695 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA0232_+_6784651.23-P1P2 (700 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153706791.23-P1P2 (705 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153707248.23-P1P2 (710 of 1926)
Performing gRNA-to-cell assignments for gRNA LAMTOR4_-_99746921.23-P1P2 (715 of 1926)
Performing gRNA-to-cell assignments for gRNA LARP4_-_50794864.23-P1P2 (720 of 1926)
Performing gRNA-to-cell assignments for gRNA LATS2_-_21635336.23-P1P2 (725 of 1926)
Performing gRNA-to-cell assignments for gRNA LEMD3_+_65563441.23-P1P2 (730 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_+_54960428.23-P1P2 (735 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_-_54960616.23-P1P2 (740 of 1926)
Performing gRNA-to-cell assignments for gRNA LHFPL2_-_77944519.23-P1P2 (745 of 1926)
Performing gRNA-to-cell assignments for gRNA LSM12_-_42143603.23-P1P2 (750 of 1926)
Performing gRNA-to-cell assignments for gRNA MANBAL_+_35918136.23-P1P2 (755 of 1926)
Performing gRNA-to-cell assignments for gRNA MAU2_+_19431732.23-P1P2 (760 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_+_47807678.23-P1P2 (765 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_-_47807877.23-P1P2 (770 of 1926)
Performing gRNA-to-cell assignments for gRNA MEAF6_-_37980319.23-P1P2 (775 of 1926)
Performing gRNA-to-cell assignments for gRNA MED1_-_37607140.23-P1P2 (780 of 1926)
Performing gRNA-to-cell assignments for gRNA MED12_+_70338541.23-P1P2 (785 of 1926)
Performing gRNA-to-cell assignments for gRNA MED15_+_20862189.23-P1 (790 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_+_57479335.23-P1P2 (795 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_-_57479437.23-P1P2 (800 of 1926)
Performing gRNA-to-cell assignments for gRNA MEMO1_-_32264882.23-P1 (805 of 1926)
Performing gRNA-to-cell assignments for gRNA MIXL1_+_226411425.23-P1P2 (810 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL28_+_420487.23-P1P2 (815 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL32_+_42972018.23-P1P2 (820 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_+_179322181.23-P1P2 (825 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_-_179322213.23-P1P2 (830 of 1926)
Performing gRNA-to-cell assignments for gRNA LUZP6_-_135612750.23-P1P2 (835 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA20_-_19998023.23-P1P2 (840 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA25_+_112546538.23-P1P2 (845 of 1926)
Performing gRNA-to-cell assignments for gRNA NCBP1_+_100396082.23-P1P2 (850 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_+_172428630.23-P1P2 (855 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_-_172428841.23-P1P2 (860 of 1926)
Performing gRNA-to-cell assignments for gRNA NEDD8_-_24701363.23-P1P2 (865 of 1926)
Performing gRNA-to-cell assignments for gRNA NF2_+_29999635.23-P1P2 (870 of 1926)
Performing gRNA-to-cell assignments for gRNA NFIB_+_14314559.23-P1P2 (875 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_+_139418329.23-P2 (880 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139418349.23-P2 (885 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139440503.23-P1 (890 of 1926)
Performing gRNA-to-cell assignments for gRNA NUFIP2_+_27621073.23-P1P2 (895 of 1926)
Performing gRNA-to-cell assignments for gRNA NUP93_-_56764090.23-P1P2 (900 of 1926)
Performing gRNA-to-cell assignments for gRNA OGT_+_70752948.23-P1P2 (905 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_+_59383170.23-P1P2 (910 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_-_59383348.23-P1P2 (915 of 1926)
Performing gRNA-to-cell assignments for gRNA PAM16_-_4401236.23-P1P2 (920 of 1926)
Performing gRNA-to-cell assignments for gRNA PAXBP1_-_34143614.23-P1P2 (925 of 1926)
Performing gRNA-to-cell assignments for gRNA PBLD_-_70092436.23-P1 (930 of 1926)
Performing gRNA-to-cell assignments for gRNA PCBP2_+_53835614.23-P2 (935 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_+_94292752.23-P2 (940 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_-_94293027.23-P2 (945 of 1926)
Performing gRNA-to-cell assignments for gRNA PGRMC2_+_129208913.23-P1P2 (950 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_-_178866211.23-P1P2 (955 of 1926)
Performing gRNA-to-cell assignments for gRNA PITPNC1_-_65373433.23-P1P2 (960 of 1926)
Performing gRNA-to-cell assignments for gRNA PLCG1_+_39766176.23-P1P2 (965 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_+_151042273.23-P2 (970 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_-_151042422.23-P2 (975 of 1926)
Performing gRNA-to-cell assignments for gRNA PLRG1_-_155471154.23-P1P2 (980 of 1926)
Performing gRNA-to-cell assignments for gRNA PLXNA2_-_208417275.23-P1P2 (985 of 1926)
Performing gRNA-to-cell assignments for gRNA POLE4_-_75185725.23-P1P2 (990 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR2M_+_57998962.23-P1P2 (995 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940045.23-P1P2 (1000 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940414.23-P1P2 (1005 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3K_-_103144.23-P1P2 (1010 of 1926)
Performing gRNA-to-cell assignments for gRNA PPHLN1_-_42720007.23-P1P2 (1015 of 1926)
Performing gRNA-to-cell assignments for gRNA PPIL4_+_149867146.23-P1P2 (1020 of 1926)
Performing gRNA-to-cell assignments for gRNA PPP2R1A_+_52693317.23-P1 (1025 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_+_53195217.23-P1P2 (1030 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_-_53195332.23-P1P2 (1035 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKD1_-_30396996.23-P1P2 (1040 of 1926)
Performing gRNA-to-cell assignments for gRNA PRMT5_+_23398554.23-P1 (1045 of 1926)
Performing gRNA-to-cell assignments for gRNA C2orf53_-_27362230.23-P1P2 (1050 of 1926)
Performing gRNA-to-cell assignments for gRNA PRSS8_+_31147010.23-P1P2 (1055 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_+_36909062.23-P1P2 (1060 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_-_36909230.23-P1P2 (1065 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMD4_-_151227289.23-P1P2 (1070 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMG4_-_3259227.23-P1P2 (1075 of 1926)
Performing gRNA-to-cell assignments for gRNA PTAR1_+_72374835.23-P1P2 (1080 of 1926)
Performing gRNA-to-cell assignments for gRNA PTBP1_+_797463.23-P1P2 (1085 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_+_89623388.23-P1P2 (1090 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_-_89623526.23-P1P2 (1095 of 1926)
Performing gRNA-to-cell assignments for gRNA PTMS_-_6875591.23-P1P2 (1100 of 1926)
Performing gRNA-to-cell assignments for gRNA PYGO2_-_154934135.23-P1P2 (1105 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB18_-_27793292.23-P1P2 (1110 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB6A_+_73471983.23-P1P2 (1115 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_+_76251918.23-P1P2 (1120 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_-_76252008.23-P1P2 (1125 of 1926)
Performing gRNA-to-cell assignments for gRNA GNB2L1_-_180670840.23-P1P2 (1130 of 1926)
Performing gRNA-to-cell assignments for gRNA RAD21_-_117886635.23-P1P2 (1135 of 1926)
Performing gRNA-to-cell assignments for gRNA RANGRF_+_8192257.23-P1P2 (1140 of 1926)
Performing gRNA-to-cell assignments for gRNA RBBP5_+_205091050.23-P1P2 (1145 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM14-RBM4_+_66406280.23-P2 (1150 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM33_+_155437802.23-P1P2 (1155 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM4_+_66406319.23-P2 (1160 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_+_114271352.23-P1P2 (1165 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_-_114271700.23-P1P2 (1170 of 1926)
Performing gRNA-to-cell assignments for gRNA RBX1_-_41347852.23-P1P2 (1175 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1551_+_32112399.23-P1P2 (1180 of 1926)
Performing gRNA-to-cell assignments for gRNA RFX3_-_3525750.23-P1P2 (1185 of 1926)
Performing gRNA-to-cell assignments for gRNA RGP1_+_35749328.23-P1P2 (1190 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_+_95067102.23-P1P2 (1195 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_-_95067164.23-P1P2 (1200 of 1926)
Performing gRNA-to-cell assignments for gRNA RIBC1_-_53449606.23-P1P2 (1205 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1432_+_5629173.23-P1P2 (1210 of 1926)
Performing gRNA-to-cell assignments for gRNA RLF_+_40627522.23-P1P2 (1215 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF115_-_145611148.23-P1P2 (1220 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_+_104296142.23-P1P2 (1225 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_-_104296250.23-P1P2 (1230 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF40_-_30773681.23-P1P2 (1235 of 1926)
Performing gRNA-to-cell assignments for gRNA RPL37A_-_217363600.23-P1P2 (1240 of 1926)
Performing gRNA-to-cell assignments for gRNA RPS15A_-_18801590.23-P1P2 (1245 of 1926)
Performing gRNA-to-cell assignments for gRNA SAV1_+_51134992.23-P1P2 (1250 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_+_33104038.23-P1P2 (1255 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_-_33104343.23-P1P2 (1260 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF8_+_155054976.23-P1P2 (1265 of 1926)
Performing gRNA-to-cell assignments for gRNA EPT1_+_26569050.23-P1P2 (1270 of 1926)
Performing gRNA-to-cell assignments for gRNA SERF2_+_44069659.23-P2 (1275 of 1926)
Performing gRNA-to-cell assignments for gRNA SERPING1_+_57365137.23-P1P2 (1280 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_+_122242041.23-P1P2 (1285 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_-_122242294.23-P1P2 (1290 of 1926)
Performing gRNA-to-cell assignments for gRNA SETDB1_-_150898815.23-P1P2 (1295 of 1926)
Performing gRNA-to-cell assignments for gRNA SF1_-_64545937.23-P1P2 (1300 of 1926)
Performing gRNA-to-cell assignments for gRNA CCDC101_-_28565268.23-P1P2 (1305 of 1926)
Performing gRNA-to-cell assignments for gRNA SKA3_+_21750519.23-P1P2 (1310 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_+_78591497.23-P1P2 (1315 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_-_78591737.23-P1P2 (1320 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC29A1_-_44187438.23-P2 (1325 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC35C1_-_45826711.23-P1P2 (1330 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD1_+_146403428.23-P1P2 (1335 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD3_+_67358229.23-P1P2 (1340 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_+_48556623.23-P1P2 (1345 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_-_48556657.23-P1P2 (1350 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD5_-_135468694.23-P1P2 (1355 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA4_-_11071741.23-P1P2 (1360 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA5_+_144435243.23-P1P2 (1365 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCB1_+_24129228.23-P1 (1370 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_+_47822975.23-P1P2 (1375 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_-_47822893.23-P1P2 (1380 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCD3_-_150974212.23-P1 (1385 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCE1_-_38803669.23-P1P2 (1390 of 1926)
Performing gRNA-to-cell assignments for gRNA SMC2_-_106856631.23-P1P2 (1395 of 1926)
Performing gRNA-to-cell assignments for gRNA SMG5_+_156252537.23-P1P2 (1400 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_+_103852.23-P1P2 (1405 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_-_104242.23-P1P2 (1410 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP27_-_70121181.23-P1P2 (1415 of 1926)
Performing gRNA-to-cell assignments for gRNA SNUPN_-_75917960.23-P1P2 (1420 of 1926)
Performing gRNA-to-cell assignments for gRNA SOX4_+_21594213.23-P1 (1425 of 1926)
Performing gRNA-to-cell assignments for gRNA SPIN1_+_91003737.23-P1P2 (1430 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_+_2328680.23-P1P2 (1435 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_-_2328912.23-P1P2 (1440 of 1926)
Performing gRNA-to-cell assignments for gRNA SPTY2D1_-_18655846.23-P1P2 (1445 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP19_+_112197188.23-P1P2 (1450 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP54_-_35452244.23-P1P2 (1455 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP68_+_74068540.23-P1P2 (1460 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_+_126138730.23-P1P2 (1465 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_-_126138742.23-P1P2 (1470 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPRB_-_133524719.23-P1P2 (1475 of 1926)
Performing gRNA-to-cell assignments for gRNA SRRM2_-_2802714.23-P1P2 (1480 of 1926)
Performing gRNA-to-cell assignments for gRNA SRSF6_+_42086993.23-P1 (1485 of 1926)
Performing gRNA-to-cell assignments for gRNA SSB_+_170655450.23-P1P2 (1490 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_+_57102902.23-P1P2 (1495 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_-_57103355.23-P1P2 (1500 of 1926)
Performing gRNA-to-cell assignments for gRNA STAG2_-_123095658.23-P1 (1505 of 1926)
Performing gRNA-to-cell assignments for gRNA STK40_+_36851507.23-P1P2 (1510 of 1926)
Performing gRNA-to-cell assignments for gRNA SUPT20H_+_37633798.23-P1P2 (1515 of 1926)
Performing gRNA-to-cell assignments for gRNA SUZ12_+_30264157.23-P1P2 (1520 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_+_64901973.23-P1P2 (1525 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_-_64901981.23-P1P2 (1530 of 1926)
Performing gRNA-to-cell assignments for gRNA TADA2B_-_7045121.23-P1P2 (1535 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF10_+_6633406.23-P1P2 (1540 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF2_+_120845032.23-P1P2 (1545 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF3_-_7860540.23-P1P2 (1550 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_+_60640814.23-P1P2 (1555 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_-_60640498.23-P1P2 (1560 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF6L_-_62538954.23-P1P2 (1565 of 1926)
Performing gRNA-to-cell assignments for gRNA TAOK1_-_27717907.23-P1P2 (1570 of 1926)
Performing gRNA-to-cell assignments for gRNA FAM208A_+_56717273.23-P1P2 (1575 of 1926)
Performing gRNA-to-cell assignments for gRNA TBC1D9_+_141677481.23-P1P2 (1580 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_+_170863524.23-P1P2 (1585 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_-_170863481.23-P1P2 (1590 of 1926)
Performing gRNA-to-cell assignments for gRNA TBX3_-_115121998.23-P1P2 (1595 of 1926)
Performing gRNA-to-cell assignments for gRNA T_-_166581333.23-P1P2 (1600 of 1926)
Performing gRNA-to-cell assignments for gRNA TGFBR2_+_30648526.23-P1P2 (1605 of 1926)
Performing gRNA-to-cell assignments for gRNA THBD_+_23029196.23-P2 (1610 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_+_2977592.23-P1P2 (1615 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_-_2977643.23-P1P2 (1620 of 1926)
Performing gRNA-to-cell assignments for gRNA MESDC1_-_81293340.23-P1P2 (1625 of 1926)
Performing gRNA-to-cell assignments for gRNA TMEM41B_+_9336156.23-P1P2 (1630 of 1926)
Performing gRNA-to-cell assignments for gRNA TMTC3_-_88536127.23-P1P2 (1635 of 1926)
Performing gRNA-to-cell assignments for gRNA TOP1_+_39657538.23-P1P2 (1640 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_+_7590531.23-P1P2 (1645 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_-_7590704.23-P1P2 (1650 of 1926)
Performing gRNA-to-cell assignments for gRNA TRIM24_-_138145308.23-P1P2 (1655 of 1926)
Performing gRNA-to-cell assignments for gRNA TRRAP_-_98476465.23-P1P2 (1660 of 1926)
Performing gRNA-to-cell assignments for gRNA TSTD2_+_100395594.23-P1P2 (1665 of 1926)
Performing gRNA-to-cell assignments for gRNA UBE2E1_+_23847541.23-P1P2 (1670 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_+_38924011.23-P1P2 (1675 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_-_38924421.23-P1P2 (1680 of 1926)
Performing gRNA-to-cell assignments for gRNA UPF2_+_12084799.23-P1P2 (1685 of 1926)
Performing gRNA-to-cell assignments for gRNA UXT_-_47518505.23-P1P2 (1690 of 1926)
Performing gRNA-to-cell assignments for gRNA VHL_-_10183525.23-P1P2 (1695 of 1926)
Performing gRNA-to-cell assignments for gRNA VPS54_+_64246104.23-P1P2 (1700 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_+_28821698.23-P1P2 (1705 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_-_28822074.23-P1P2 (1710 of 1926)
Performing gRNA-to-cell assignments for gRNA WAPAL_+_88281533.23-P1P2 (1715 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR12_-_203776383.23-P1P2 (1720 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR77_+_111991860.23-P1P2 (1725 of 1926)
Performing gRNA-to-cell assignments for gRNA WIZ_-_15543016.23-P1P2 (1730 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_+_9983328.23-P1P2 (1735 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_-_9983510.23-P1P2 (1740 of 1926)
Performing gRNA-to-cell assignments for gRNA XPO5_-_43543670.23-P1P2 (1745 of 1926)
Performing gRNA-to-cell assignments for gRNA XRN1_-_142166648.23-P1P2 (1750 of 1926)
Performing gRNA-to-cell assignments for gRNA YARS2_+_32908787.23-P1P2 (1755 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB22_+_33285454.23-P1P2 (1760 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_+_6640093.23-P1P2 (1765 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_-_6640138.23-P1P2 (1770 of 1926)
Performing gRNA-to-cell assignments for gRNA ZC3H8_-_113012435.23-P1P2 (1775 of 1926)
Performing gRNA-to-cell assignments for gRNA ZNF281_+_200379151.23-P1P2 (1780 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00117 (1785 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00233 (1790 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00518 (1795 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00656 (1800 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01226 (1805 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01730 (1810 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02009 (1815 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02288 (1820 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02777 (1825 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03026 (1830 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03252 (1835 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03620 (1840 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_SGTN_04 (1845 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_SGTN_25 (1850 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_SGTN_26 (1855 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1044 (1860 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1414 (1865 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1554 (1870 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1785 (1875 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2009 (1880 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2271 (1885 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2425 (1890 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2614 (1895 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_293 (1900 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_541 (1905 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_725 (1910 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_860 (1915 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_SGTN_14 (1920 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_SGTN_01 (1925 of 1926)

```

## Code cell 35

```python
%%R
library(Matrix)

# 1. Get the assignment matrix
A <- get_grna_assignments(sceptre_object_highmoi)

# 2. Ensure it is a numeric sparse matrix (1/0)
if (inherits(A, "Matrix") || is.matrix(A)) {
    A_num <- if (inherits(A, "Matrix")) drop0(A * 1) else (A * 1)
} else {
    stop("Output is not a matrix.")
}

# --- FIX: ROBUST NAME RECOVERY ---
# We retrieve names using Sceptre's standard metadata accessors
# instead of touching internal slots like @grna_matrix.

# Try to get cell names from the matrix first
current_cell_names <- colnames(A_num)
current_guide_names <- rownames(A_num)

# If missing, grab them from the main object's covariate data frame
if (is.null(current_cell_names)) {
    message("⚠️ Matrix missing colnames. Fetching from sceptre_object...")
    # The covariate data frame ALWAYS has the cell IDs as rownames
    current_cell_names <- rownames(sceptre_object_highmoi@covariate_data_frame)
    
    # Safety check: lengths must match
    if (length(current_cell_names) != ncol(A_num)) {
        stop(paste("Mismatch! Matrix has", ncol(A_num), "cols but object has", length(current_cell_names), "cells."))
    }
    colnames(A_num) <- current_cell_names
}

# If guide names are missing, they should match the unique values in your guide mapping
if (is.null(current_guide_names)) {
    message("⚠️ Matrix missing rownames. Fetching from guide mapping...")
    # Get unique targets from the mapping table you loaded earlier
    # Note: This assumes the matrix rows are sorted/aligned to these unique targets.
    # In 'mixture' method, rows usually correspond to the 'grna_target' labels.
    
    # SAFE FALLBACK: If we can't be 100% sure of order, we warn the user.
    # But usually, they match the 'grna_target' column from the input.
    unique_targets <- sort(unique(guide_mapping$grna_target))
    
    if (length(unique_targets) == nrow(A_num)) {
         rownames(A_num) <- unique_targets
    } else {
         warning("Could not automatically align target names. Using generic IDs.")
         rownames(A_num) <- paste0("target_", seq_len(nrow(A_num)))
    }
}
# ---------------------------------

# 3. Convert to Long Format
s <- summary(A_num)
assign_long <- data.frame(
    cell_id = colnames(A_num)[s$j],
    grna_target = rownames(A_num)[s$i],
    value = s$x
)

# 4. Save
out_path <- file.path(data_dir, "grna_assignments_long.csv")
write.csv(assign_long, out_path, row.names = FALSE)

message("✅ Saved long CSV to: ", out_path)
print(head(assign_long))
```

### Text output (stream)

```text
                     cell_id               grna_target value
1 CCTAATCCAACAACCATCTATTAC_2    HGF_+_81399312.23-P1P2     1
2 CCTAATCCAACAACCATCTATTAC_2  SMG5_+_156252585.23-P1P2     1
3 CCTAATCCAACAACCATCTATTAC_2 WDR77_+_111991472.23-P1P2     1
4 CCTAATCCAACAACCATCTATTAC_2 ZBTB22_-_33285488.23-P1P2     1
5 CCTAATCCAACAACCATCTATTAC_2             WDR61_SGTN_13     1
6 GGTGCGAACCTCCTGAATCCGCGA_2   EIF5A_-_7210911.23-P1P2     1

```

### Text output (text/plain)

```text
⚠️ Matrix missing colnames. Fetching from sceptre_object...
✅ Saved long CSV to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d0_sceptre_transfer/grna_assignments_long.csv

```

## Code cell 36

```python
%%R
png(file.path(data_dir, "MOI.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "MOI.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d0_sceptre_transfer/MOI.png

```

## Code cell 37

```python
%%R
sceptre_object_highmoi <- run_qc(
  sceptre_object = sceptre_object_highmoi,
  p_mito_threshold = 0.25
)

png(file.path(data_dir, "QC.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "QC.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d0_sceptre_transfer/QC.png

```

## Code cell 38

```python
%%R
# Calibration check


sceptre_object_highmoi <- run_calibration_check(
  sceptre_object = sceptre_object_highmoi,
  parallel = FALSE
)

png(file.path(data_dir, "calibration_check.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "calibration_check.png"))

```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Constructing negative control pairs. ✓
Running precomputation on response ENSG00000160741.17 (1 of 259)
Running precomputation on response ENSG00000286221.1 (5 of 259)
Running precomputation on response ENSG00000162971.11 (10 of 259)
Running precomputation on response ENSG00000183087.15 (15 of 259)
Running precomputation on response ENSG00000228566.2 (20 of 259)
Running precomputation on response ENSG00000213235.3 (25 of 259)
Running precomputation on response ENSG00000120645.12 (30 of 259)
Running precomputation on response ENSG00000149100.13 (35 of 259)
Running precomputation on response ENSG00000142046.15 (40 of 259)
Running precomputation on response ENSG00000143748.18 (45 of 259)
Running precomputation on response ENSG00000149311.20 (50 of 259)
Running precomputation on response ENSG00000133980.5 (55 of 259)
Running precomputation on response ENSG00000196312.16 (60 of 259)
Running precomputation on response ENSG00000113522.14 (65 of 259)
Running precomputation on response ENSG00000259673.6 (70 of 259)
Running precomputation on response ENSG00000100580.8 (75 of 259)
Running precomputation on response ENSG00000235655.3 (80 of 259)
Running precomputation on response ENSG00000196323.14 (85 of 259)
Running precomputation on response ENSG00000204520.14 (90 of 259)
Running precomputation on response ENSG00000146834.15 (95 of 259)
Running precomputation on response ENSG00000102978.13 (100 of 259)
Running precomputation on response ENSG00000164758.7 (105 of 259)
Running precomputation on response ENSG00000171097.14 (110 of 259)
Running precomputation on response ENSG00000111704.11 (115 of 259)
Running precomputation on response ENSG00000008086.13 (120 of 259)
Running precomputation on response ENSG00000175390.15 (125 of 259)
Running precomputation on response ENSG00000100300.18 (130 of 259)
Running precomputation on response ENSG00000196678.14 (135 of 259)
Running precomputation on response ENSG00000119638.14 (140 of 259)
Running precomputation on response ENSG00000116120.11 (145 of 259)
Running precomputation on response ENSG00000112667.13 (150 of 259)
Running precomputation on response ENSG00000111203.13 (155 of 259)
Running precomputation on response ENSG00000160606.11 (160 of 259)
Running precomputation on response ENSG00000242612.7 (165 of 259)
Running precomputation on response ENSG00000103404.15 (170 of 259)
Running precomputation on response ENSG00000097046.13 (175 of 259)
Running precomputation on response ENSG00000213722.10 (180 of 259)
Running precomputation on response ENSG00000133104.14 (185 of 259)
Running precomputation on response ENSG00000107771.17 (190 of 259)
Running precomputation on response ENSG00000187720.15 (195 of 259)
Running precomputation on response ENSG00000108344.15 (200 of 259)
Running precomputation on response ENSG00000135424.19 (205 of 259)
Running precomputation on response ENSG00000198408.14 (210 of 259)
Running precomputation on response ENSG00000162374.18 (215 of 259)
Running precomputation on response ENSG00000187987.10 (220 of 259)
Running precomputation on response ENSG00000126464.15 (225 of 259)
Running precomputation on response ENSG00000125149.12 (230 of 259)
Running precomputation on response ENSG00000154079.6 (235 of 259)
Running precomputation on response ENSG00000162769.13 (240 of 259)
Running precomputation on response ENSG00000135763.10 (245 of 259)
Running precomputation on response ENSG00000290018.1 (250 of 259)
Running precomputation on response ENSG00000102241.12 (255 of 259)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00259&non-targeting_02159&non-targeting_02571&non-targeting_03029&non-targeting_03081 (1 of 95)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00355&non-targeting_02025&non-targeting_02089&non-targeting_03081&non-targeting_03384 (5 of 95)
Analyzing pairs containing gRNA group non-targeting_00518&non-targeting_03081&non-targeting_03342&non-targeting_03427&non-targeting_03635&non-targeting_03664 (10 of 95)
Analyzing pairs containing gRNA group non-targeting_01509&non-targeting_01812&non-targeting_02025&non-targeting_02677&non-targeting_02777&non-targeting_03384 (15 of 95)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_00231&non-targeting_01185&non-targeting_02025&non-targeting_03384&non-targeting_03427 (20 of 95)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_02193&non-targeting_02543&non-targeting_02837&non-targeting_03026&non-targeting_03342 (25 of 95)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_00117&non-targeting_00196&non-targeting_00555&non-targeting_02025&non-targeting_02193 (30 of 95)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00355&non-targeting_00518&non-targeting_01601&non-targeting_02677&non-targeting_03427 (35 of 95)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00390&non-targeting_01601&non-targeting_02009&non-targeting_02288&non-targeting_03029 (40 of 95)
Analyzing pairs containing gRNA group non-targeting_00231&non-targeting_00614&non-targeting_00629&non-targeting_00918&non-targeting_02193&non-targeting_02543 (45 of 95)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00390&non-targeting_01732&non-targeting_02193&non-targeting_03081&non-targeting_03384 (50 of 95)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_00196&non-targeting_00390&non-targeting_01730&non-targeting_03081&non-targeting_03252 (55 of 95)
Analyzing pairs containing gRNA group non-targeting_00043&non-targeting_00355&non-targeting_01601&non-targeting_02814&non-targeting_03081&non-targeting_03620 (60 of 95)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_01185&non-targeting_03029&non-targeting_03081&non-targeting_03180&non-targeting_03384 (65 of 95)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00259&non-targeting_00972&non-targeting_01601&non-targeting_02089&non-targeting_03620 (70 of 95)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_01858&non-targeting_02777&non-targeting_03158&non-targeting_03252&non-targeting_03483 (75 of 95)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00614&non-targeting_00625&non-targeting_01601&non-targeting_02837&non-targeting_03180 (80 of 95)
Analyzing pairs containing gRNA group non-targeting_00629&non-targeting_00656&non-targeting_02193&non-targeting_02543&non-targeting_03026&non-targeting_03483 (85 of 95)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_00213&non-targeting_00656&non-targeting_02911&non-targeting_03158&non-targeting_03483 (90 of 95)
Analyzing pairs containing gRNA group non-targeting_00629&non-targeting_00794&non-targeting_02962&non-targeting_03029&non-targeting_03620&non-targeting_03635 (95 of 95)

```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d0_sceptre_transfer/calibration_check.png
In addition: There were 50 or more warnings (use warnings() to see the first 50)

```

## Code cell 39

```python
%%R
calibration_result <- get_result(
  sceptre_object = sceptre_object_highmoi,
  analysis = "run_calibration_check"
)
head(calibration_result)
```

### Text output (stream)

```text
          response_id
               <char>
1:  ENSG00000288053.1
2: ENSG00000102893.18
3:  ENSG00000137970.7
4:  ENSG00000290682.1
5:  ENSG00000248671.9
6:  ENSG00000044574.9
                                                                                                               grna_target
                                                                                                                    <char>
1: non-targeting_00231&non-targeting_00518&non-targeting_00629&non-targeting_00794&non-targeting_02777&non-targeting_02837
2: non-targeting_00518&non-targeting_00656&non-targeting_01730&non-targeting_02574&non-targeting_03029&non-targeting_03158
3: non-targeting_00184&non-targeting_02193&non-targeting_02543&non-targeting_02837&non-targeting_03026&non-targeting_03342
4: non-targeting_00196&non-targeting_00918&non-targeting_01730&non-targeting_02025&non-targeting_02911&non-targeting_03384
5: non-targeting_00231&non-targeting_00518&non-targeting_00629&non-targeting_00794&non-targeting_02777&non-targeting_02837
6: non-targeting_00213&non-targeting_02543&non-targeting_02777&non-targeting_02814&non-targeting_03252&non-targeting_03483
   n_nonzero_trt n_nonzero_cntrl pass_qc     p_value log_2_fold_change
           <int>           <int>  <lgcl>       <num>             <num>
1:            25            1171    TRUE 0.001359857         0.7726258
2:            65            5944    TRUE 0.003678742        -0.4971951
3:           137            8700    TRUE 0.004943804        -0.3413673
4:            62            3031    TRUE 0.007133584         0.4393944
5:            16             845    TRUE 0.007279064         0.7678848
6:           343           24721    TRUE 0.009819351         0.1739685
   significant
        <lgcl>
1:       FALSE
2:       FALSE
3:       FALSE
4:       FALSE
5:       FALSE
6:       FALSE

```

## Code cell 40

```python
%%R
sceptre_object_highmoi <- run_discovery_analysis(
  sceptre_object = sceptre_object_highmoi,
  parallel = FALSE
)

png(file.path(data_dir, "discovery.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "discovery.png"))
```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Running precomputation on response ENSG00000004487.18 (1 of 261)
Running precomputation on response ENSG00000117000.9 (5 of 261)
Running precomputation on response ENSG00000265491.5 (10 of 261)
Running precomputation on response ENSG00000179085.8 (15 of 261)
Running precomputation on response ENSG00000119185.13 (20 of 261)
Running precomputation on response ENSG00000143952.20 (25 of 261)
Running precomputation on response ENSG00000138385.16 (30 of 261)
Running precomputation on response ENSG00000134086.9 (35 of 261)
Running precomputation on response ENSG00000173473.11 (40 of 261)
Running precomputation on response ENSG00000144959.11 (45 of 261)
Running precomputation on response ENSG00000128052.10 (50 of 261)
Running precomputation on response ENSG00000153147.6 (55 of 261)
Running precomputation on response ENSG00000113360.17 (60 of 261)
Running precomputation on response ENSG00000113658.18 (65 of 261)
Running precomputation on response ENSG00000180822.12 (70 of 261)
Running precomputation on response ENSG00000204371.12 (75 of 261)
Running precomputation on response ENSG00000112249.14 (80 of 261)
Running precomputation on response ENSG00000112592.15 (85 of 261)
Running precomputation on response ENSG00000242265.6 (90 of 261)
Running precomputation on response ENSG00000105887.11 (95 of 261)
Running precomputation on response ENSG00000104412.8 (100 of 261)
Running precomputation on response ENSG00000107036.12 (105 of 261)
Running precomputation on response ENSG00000136925.15 (110 of 261)
Running precomputation on response ENSG00000136908.18 (115 of 261)
Running precomputation on response ENSG00000151461.20 (120 of 261)
Running precomputation on response ENSG00000107779.14 (125 of 261)
Running precomputation on response ENSG00000198730.9 (130 of 261)
Running precomputation on response ENSG00000149136.9 (135 of 261)
Running precomputation on response ENSG00000162227.8 (140 of 261)
Running precomputation on response ENSG00000173933.21 (145 of 261)
Running precomputation on response ENSG00000182934.12 (150 of 261)
Running precomputation on response ENSG00000174718.12 (155 of 261)
Running precomputation on response ENSG00000197111.16 (160 of 261)
Running precomputation on response ENSG00000120868.14 (165 of 261)
Running precomputation on response ENSG00000165480.16 (170 of 261)
Running precomputation on response ENSG00000100462.16 (175 of 261)
Running precomputation on response ENSG00000151748.15 (180 of 261)
Running precomputation on response ENSG00000255529.9 (185 of 261)
Running precomputation on response ENSG00000169371.14 (190 of 261)
Running precomputation on response ENSG00000167978.18 (195 of 261)
Running precomputation on response ENSG00000168488.19 (200 of 261)
Running precomputation on response ENSG00000102900.13 (205 of 261)
Running precomputation on response ENSG00000108961.14 (210 of 261)
Running precomputation on response ENSG00000125686.12 (215 of 261)
Running precomputation on response ENSG00000154217.16 (220 of 261)
Running precomputation on response ENSG00000005206.17 (225 of 261)
Running precomputation on response ENSG00000011451.21 (230 of 261)
Running precomputation on response ENSG00000105568.19 (235 of 261)
Running precomputation on response ENSG00000124181.15 (240 of 261)
Running precomputation on response ENSG00000159086.15 (245 of 261)
Running precomputation on response ENSG00000186575.19 (250 of 261)
Running precomputation on response ENSG00000086758.17 (255 of 261)
Running precomputation on response ENSG00000071553.18 (260 of 261)
Analyzing pairs containing gRNA group KDM1A (1 of 261)
Analyzing pairs containing gRNA group RLF (5 of 261)
Analyzing pairs containing gRNA group RNF115 (10 of 261)
Analyzing pairs containing gRNA group DPM3 (15 of 261)
Analyzing pairs containing gRNA group ITGB1BP1 (20 of 261)
Analyzing pairs containing gRNA group VPS54 (25 of 261)
Analyzing pairs containing gRNA group SSB (30 of 261)
Analyzing pairs containing gRNA group VHL (35 of 261)
Analyzing pairs containing gRNA group SMARCC1 (40 of 261)
Analyzing pairs containing gRNA group NCEH1 (45 of 261)
Analyzing pairs containing gRNA group KDR (50 of 261)
Analyzing pairs containing gRNA group SMARCA5 (55 of 261)
Analyzing pairs containing gRNA group DROSHA (60 of 261)
Analyzing pairs containing gRNA group SMAD5 (65 of 261)
Analyzing pairs containing gRNA group PSMG4 (70 of 261)
Analyzing pairs containing gRNA group EHMT2 (75 of 261)
Analyzing pairs containing gRNA group ASCC3 (80 of 261)
Analyzing pairs containing gRNA group TBP (85 of 261)
Analyzing pairs containing gRNA group PEG10 (90 of 261)
Analyzing pairs containing gRNA group MTPN (95 of 261)
Analyzing pairs containing gRNA group EMC2 (100 of 261)
Analyzing pairs containing gRNA group RIC1 (105 of 261)
Analyzing pairs containing gRNA group TSTD2 (110 of 261)
Analyzing pairs containing gRNA group DPM2 (115 of 261)
Analyzing pairs containing gRNA group UPF2 (120 of 261)
Analyzing pairs containing gRNA group BMPR1A (125 of 261)
Analyzing pairs containing gRNA group CTR9 (130 of 261)
Analyzing pairs containing gRNA group SSRP1 (135 of 261)
Analyzing pairs containing gRNA group TAF6L (140 of 261)
Analyzing pairs containing gRNA group RBM4 (145 of 261)
Analyzing pairs containing gRNA group SRPRA (150 of 261)
Analyzing pairs containing gRNA group RESF1 (155 of 261)
Analyzing pairs containing gRNA group PCBP2 (160 of 261)
Analyzing pairs containing gRNA group APAF1 (165 of 261)
Analyzing pairs containing gRNA group SKA3 (170 of 261)
Analyzing pairs containing gRNA group PRMT5 (175 of 261)
Analyzing pairs containing gRNA group SAV1 (180 of 261)
Analyzing pairs containing gRNA group POLR2M (185 of 261)
Analyzing pairs containing gRNA group SNUPN (190 of 261)
Analyzing pairs containing gRNA group SRRM2 (195 of 261)
Analyzing pairs containing gRNA group ATXN2L (200 of 261)
Analyzing pairs containing gRNA group NUP93 (205 of 261)
Analyzing pairs containing gRNA group RANGRF (210 of 261)
Analyzing pairs containing gRNA group MED1 (215 of 261)
Analyzing pairs containing gRNA group PITPNC1 (220 of 261)
Analyzing pairs containing gRNA group SPPL2B (225 of 261)
Analyzing pairs containing gRNA group WIZ (230 of 261)
Analyzing pairs containing gRNA group PPP2R1A (235 of 261)
Analyzing pairs containing gRNA group PLCG1 (240 of 261)
Analyzing pairs containing gRNA group PAXBP1 (245 of 261)
Analyzing pairs containing gRNA group NF2 (250 of 261)
Analyzing pairs containing gRNA group HUWE1 (255 of 261)
Analyzing pairs containing gRNA group ATP6AP1 (260 of 261)

```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d0_sceptre_transfer/discovery.png
In addition: There were 50 or more warnings (use warnings() to see the first 50)

```

## Code cell 41

```python
%%R
discovery_result <- get_result(
  sceptre_object = sceptre_object_highmoi,
  analysis = "run_discovery_analysis"
)

write.csv(discovery_result, file.path(data_dir, "TargetKD.csv"), row.names = FALSE)
message("✅ Saved discovery to: ", file.path(data_dir, "TargetKD.csv"))

```

### Text output (text/plain)

```text
✅ Saved discovery to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d0_sceptre_transfer/TargetKD.csv

```

## Markdown cell 42

### Day 1

## Code cell 43

```python
import scipy.io as sio
import pandas as pd
import os
```

## Code cell 44

```python
d1_mdata = mu.read(f"{data_dir}/d1_normalized_gex_and_guide.h5mu")
mdata = d1_mdata

# Filter out cells that have no guides
#mdata = mdata[mdata['GEX'].obs['guides_per_cell'] > 0]
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

## Code cell 45

```python
# Define a temporary folder for transfer
temp_dir = working_dir + "/d1_sceptre_transfer"
os.makedirs(temp_dir, exist_ok=True)

# 1. EXTRACT MATRICES (Transpose them!)
# Scanpy is (Cells × Genes), Sceptre needs (Genes × Cells)
print("Writing sparse matrices...")
sio.mmwrite(f"{temp_dir}/response_matrix.mtx", mdata['GEX'].layers['raw_counts'].T)
sio.mmwrite(f"{temp_dir}/grna_matrix.mtx", mdata['guide'].X.T)

# 2. EXTRACT NAMES (Row/Col names are lost in .mtx)
# Save Gene names and Guide names to assign them later in R
pd.Series(mdata['GEX'].var_names).to_csv(f"{temp_dir}/gene_names.csv", index=False, header=False)
pd.Series(mdata['guide'].var_names).to_csv(f"{temp_dir}/guide_names.csv", index=False, header=False)
pd.Series(mdata['GEX'].obs_names).to_csv(f"{temp_dir}/cell_names.csv", index=False, header=False)
```

### Text output (stream)

```text
Writing sparse matrices...

```

## Code cell 46

```python
mdata['GEX']
```

### Text output (text/plain)

```text
AnnData object with n_obs × n_vars = 16925 × 13383
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0'
    var: 'gene_id', 'symbol', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'n_counts', 'highly_deviant', 'highly_deviant_2k', 'binomial_deviance', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'orig_highly_variable', 'highly_variable_2k'
    uns: 'hvg', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'neighbors', 'pca', 'umap'
    obsm: 'X_pca', 'X_umap', 'guide_assignment', 'guide_counts'
    varm: 'PCs'
    layers: 'ambiguous', 'analytic_pearson_residuals', 'log1p_norm', 'mature', 'nascent', 'raw_counts', 'scran_normalization'
    obsp: 'connectivities', 'distances'
```

## Code cell 47

```python
# 3. EXTRACT COVARIATES 
rna = mdata['GEX']
covariates_df = pd.DataFrame({
    "batch": rna.obs['sample_id'],
    "grna_n_nonzero": rna.obs['guides_per_cell'],
    "grna_n_umis": rna.obs['guide_umi_counts'],
    "response_n_nonzero": rna.obs['n_genes_by_counts'],
    "response_n_umis": rna.obs['total_counts'],
    "response_p_mito": rna.obs['pct_counts_mt'],
})
covariates_df.to_csv(f"{temp_dir}/covariates.csv")
covariates_df
```

### Text output (text/plain)

```text
                                             batch  grna_n_nonzero  \
barcode                                                              
AACGTGATGTGTTCTATACGGTTA_2  300genes:eTN5_GEX_sub2               3   
CCATCCTCCAACCACATTGGGAGA_3  300genes:eTN5_GEX_sub3              10   
AGTACAAGGCGAGTAATTATTCTG_1  300genes:eTN5_GEX_sub1               6   
GTACGCAAACCACTGTTACGGTTA_2  300genes:eTN5_GEX_sub2               7   
ACAGATTCCCGACAACTACGGTTA_1  300genes:eTN5_GEX_sub1               6   
...                                            ...             ...   
AAACATCGCATACCAATTGTCTTA_2  300genes:eTN5_GEX_sub2               2   
CATCAAGTACTATGCATTGTCTTA_3  300genes:eTN5_GEX_sub3              12   
CACCTTACGCTCGGTATTATTCTG_1  300genes:eTN5_GEX_sub1               0   
GACAGTGCAACAACCACGCGACTA_3  300genes:eTN5_GEX_sub3               0   
ACAAGCTACCATCCTCTACGGTTA_1  300genes:eTN5_GEX_sub1               3   

                            grna_n_umis  response_n_nonzero  response_n_umis  \
barcode                                                                        
AACGTGATGTGTTCTATACGGTTA_2        321.0                1769      2897.497314   
CCATCCTCCAACCACATTGGGAGA_3        463.0                1605      2809.700684   
AGTACAAGGCGAGTAATTATTCTG_1        155.0                1636      2930.134033   
GTACGCAAACCACTGTTACGGTTA_2        102.0                 527       720.762817   
ACAGATTCCCGACAACTACGGTTA_1        332.0                 854      1462.475586   
...                                 ...                 ...              ...   
AAACATCGCATACCAATTGTCTTA_2         43.0                1680      2683.880615   
CATCAAGTACTATGCATTGTCTTA_3        424.0                3617      8417.620117   
CACCTTACGCTCGGTATTATTCTG_1         19.0                1884      3545.153809   
GACAGTGCAACAACCACGCGACTA_3         42.0                2425      3860.632080   
ACAAGCTACCATCCTCTACGGTTA_1         53.0                 577       824.888916   

                            response_p_mito  
barcode                                      
AACGTGATGTGTTCTATACGGTTA_2        12.067884  
CCATCCTCCAACCACATTGGGAGA_3        18.792038  
AGTACAAGGCGAGTAATTATTCTG_1        15.733069  
GTACGCAAACCACTGTTACGGTTA_2        24.765429  
ACAGATTCCCGACAACTACGGTTA_1        18.803732  
...                                     ...  
AAACATCGCATACCAATTGTCTTA_2         9.948281  
CATCAAGTACTATGCATTGTCTTA_3         8.529727  
CACCTTACGCTCGGTATTATTCTG_1         9.999566  
GACAGTGCAACAACCACGCGACTA_3        12.225977  
ACAAGCTACCATCCTCTACGGTTA_1        24.609373  

[16925 rows x 6 columns]
```

## Code cell 48

```python
guide_meta_df = pd.read_csv(f"{working_dir}/300genes_guide_metadata_v43.tsv", sep='\t')
```

## Code cell 49

```python
# 4. PREPARE GUIDE MAPPING 
# Sceptre needs a table mapping 'grna_id' to 'grna_target'.

# 1. Create a "skeleton" DataFrame from the actual matrix order
# This ensures the order matches mdata['guide'].X perfectly
# We name the index 'grna_id' to match Sceptre's expectation immediately
guide_mapping = pd.DataFrame(index=mdata['guide'].var_names)
guide_mapping.index.name = 'grna_id'

# 2. Prepare metadata for merging
# Ensure 'guide_id' is the index so we can map efficiently
lookup_df = guide_meta_df.set_index('guide_id')

# 3. Map the targets onto the skeleton
# We use .loc[] to pull the specific target for each guide in the exact order required
guide_mapping['grna_target'] = lookup_df.loc[guide_mapping.index, 'intended_target_name']

# 4. Reset index so 'grna_id' becomes a column (Sceptre requires this)
guide_mapping = guide_mapping.reset_index()

# ==========================================
# VALIDATION CHECKS 
# ==========================================

# Check 1: Are there any NAs? (Guides in mdata that were missing from meta_df)
na_count = guide_mapping['grna_target'].isna().sum()
if na_count > 0:
    print(f"⚠️ WARNING: {na_count} guides have NaN targets! These will crash Sceptre.")
else:
    print("✅ All guides mapped successfully.")

# Check 2: Verify the first few rows manually
print("\nFirst 5 rows of mapping (Must match mdata['guide'].var_names order):")
print(guide_mapping.head())

# 5. Save for Sceptre
guide_mapping.to_csv(f"{temp_dir}/guide_mapping.csv", index=False)
```

### Text output (stream)

```text
✅ All guides mapped successfully.

First 5 rows of mapping (Must match mdata['guide'].var_names order):
                    grna_id grna_target
0  ACAA1_+_38178488.23-P1P2       ACAA1
1  ACAA1_+_38178517.23-P1P2       ACAA1
2  ACAA1_+_38178559.23-P1P2       ACAA1
3  ACAA1_+_38178570.23-P1P2       ACAA1
4  ACAA1_+_38178620.23-P1P2       ACAA1

```

## Code cell 50

```python
guide_mapping
```

### Text output (text/plain)

```text
                       grna_id grna_target
0     ACAA1_+_38178488.23-P1P2       ACAA1
1     ACAA1_+_38178517.23-P1P2       ACAA1
2     ACAA1_+_38178559.23-P1P2       ACAA1
3     ACAA1_+_38178570.23-P1P2       ACAA1
4     ACAA1_+_38178620.23-P1P2       ACAA1
...                        ...         ...
1921             SUZ12_SGTN_23       SUZ12
1922            MESDC1_SGTN_11      TLNRD1
1923            MESDC1_SGTN_12      TLNRD1
1924              TP53_SGTN_01        TP53
1925              TP53_SGTN_02        TP53

[1926 rows x 2 columns]
```

## Code cell 51

```python
# Create a mapping from the Gene Symbol (which matches your guides)
# to the Gene Index (which matches the matrix rownames)

# Check which column in mdata['GEX'].var holds the symbol
symbol_col = 'symbol' 

gene_map = mdata['GEX'].var[[symbol_col]].reset_index()
gene_map.columns = ['response_id', 'grna_target'] # Rename for easy merging in R

# Save it
gene_map.to_csv(f"{temp_dir}/gene_symbol_map.csv", index=False)
print("✅ Gene mapping saved.")
```

### Text output (stream)

```text
✅ Gene mapping saved.

```

## Code cell 52

```python
gene_map
```

### Text output (text/plain)

```text
             response_id      grna_target
0      ENSG00000227232.5           WASH7P
1      ENSG00000238009.6  ENSG00000238009
2      ENSG00000268903.1  ENSG00000268903
3      ENSG00000241860.7  ENSG00000241860
4      ENSG00000279457.4           WASH9P
...                  ...              ...
13378  ENSG00000275757.1        5_8S_rRNA
13379  ENSG00000277196.4            PRODH
13380  ENSG00000276256.1  ENSG00000276256
13381  ENSG00000274917.1        RNA5-8SN5
13382  ENSG00000273730.1        5_8S_rRNA

[13383 rows x 2 columns]
```

## Markdown cell 53

#### Import data to SCEPTRE (Rcell)

## Code cell 54

```python
%load_ext rpy2.ipython
```

### Text output (stream)

```text
The rpy2.ipython extension is already loaded. To reload it, use:
  %reload_ext rpy2.ipython

```

## Code cell 55

```python
%%R
library(sceptre)
library(Matrix)

# Define path
data_dir <- "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d1_sceptre_transfer"

# 1. LOAD MATRICES
cat("Loading matrices into R...\n")
response_mat <- readMM(file.path(data_dir, "response_matrix.mtx"))
grna_mat <- readMM(file.path(data_dir, "grna_matrix.mtx"))

# 2. ATTACH NAMES
gene_names <- read.csv(file.path(data_dir, "gene_names.csv"), header=FALSE)$V1
guide_names <- read.csv(file.path(data_dir, "guide_names.csv"), header=FALSE)$V1
cell_names <- read.csv(file.path(data_dir, "cell_names.csv"), header=FALSE)$V1

rownames(response_mat) <- gene_names
colnames(response_mat) <- cell_names

rownames(grna_mat) <- guide_names
colnames(grna_mat) <- cell_names
```

### Text output (stream)

```text
Loading matrices into R...

```

## Code cell 56

```python
%%R

# 3. LOAD METADATA
covariates <- read.csv(file.path(data_dir, "covariates.csv"), row.names=1)
guide_mapping <- read.csv(file.path(data_dir, "guide_mapping.csv"))

# ==============================================================================
# FIX: RENAME 'response_p_mito' SO WE KEEP IT
# Sceptre blocks the name 'response_p_mito' but doesn't calculate it.
# We rename it to 'mt_percent' so we can use it as a covariate later.
# ==============================================================================

# 1. Rename response_p_mito if it exists
if ("response_p_mito" %in% colnames(covariates)) {
    cat("Renaming 'response_p_mito' to 'mt_percent'...\n")
    names(covariates)[names(covariates) == "response_p_mito"] <- "mt_percent"
}

# 2. DROP OTHER RESERVED COLUMNS (which Sceptre DOES calculate)
# We still drop n_umis and n_nonzero because Sceptre re-calculates those automatically.
reserved_cols <- c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis")
covariates <- covariates[, !(colnames(covariates) %in% reserved_cols), drop=FALSE]

# 4. CREATE SCEPTRE OBJECT
sceptre_object <- import_data(
  response_matrix = response_mat,
  grna_matrix = grna_mat,
  grna_target_data_frame = guide_mapping,
  extra_covariates = covariates,
  moi = "high"
)

print("✅ Sceptre object created!")
print(head(sceptre_object@covariate_data_frame))
```

### Text output (stream)

```text
Renaming 'response_p_mito' to 'mt_percent'...
[1] "\u2705 Sceptre object created!"
                           response_n_nonzero response_n_umis grna_n_nonzero
AACGTGATGTGTTCTATACGGTTA_2               1487       2897.4972             19
CCATCCTCCAACCACATTGGGAGA_3               1382       2809.7007             48
AGTACAAGGCGAGTAATTATTCTG_1               1415       2930.1340             50
GTACGCAAACCACTGTTACGGTTA_2                448        720.7629             23
ACAGATTCCCGACAACTACGGTTA_1                755       1462.4755             57
AGAGTCAAGGAGAACATTGTCTTA_1               1086       2194.7055             42
                           grna_n_umis                  batch mt_percent
AACGTGATGTGTTCTATACGGTTA_2         321 300genes:eTN5_GEX_sub2   12.06788
CCATCCTCCAACCACATTGGGAGA_3         463 300genes:eTN5_GEX_sub3   18.79204
AGTACAAGGCGAGTAATTATTCTG_1         155 300genes:eTN5_GEX_sub1   15.73307
GTACGCAAACCACTGTTACGGTTA_2         102 300genes:eTN5_GEX_sub2   24.76543
ACAGATTCCCGACAACTACGGTTA_1         332 300genes:eTN5_GEX_sub1   18.80373
AGAGTCAAGGAGAACATTGTCTTA_1         244 300genes:eTN5_GEX_sub1   10.79872

```

## Code cell 57

```python
%%R
# 1. Load the translation map
gene_map <- read.csv(file.path(data_dir, "gene_symbol_map.csv"))

# 2. Start with guide mapping (which has 'grna_target' as Symbol)
discovery_pairs_highmoi <- guide_mapping

# 3. Merge to attach the correct response_id (Gene ID)
#    This matches 'grna_target' in guides to 'grna_target' in the gene map
discovery_pairs_highmoi <- merge(discovery_pairs_highmoi, gene_map, by="grna_target")

# 4. Select ONLY the two columns Sceptre needs
#    (We filter out non-targeting guides automatically via the merge)
discovery_pairs_highmoi <- discovery_pairs_highmoi[, c("grna_target", "response_id")]

# 5. DEDUPLICATE
#    Sceptre performs analysis at the *Target* level, not the *Guide* level.
#    We only need one row per Target-Gene pair.
discovery_pairs_highmoi <- unique(discovery_pairs_highmoi)

print(paste("Unique pairs to test:", nrow(discovery_pairs_highmoi)))
print(head(discovery_pairs_highmoi))
```

### Text output (stream)

```text
[1] "Unique pairs to test: 281"
   grna_target        response_id
1        ACAA1 ENSG00000060971.19
7       ACVR1B ENSG00000135503.13
13      ACVR2A ENSG00000121989.15
19      ACVR2B ENSG00000114739.14
25       ADCY2 ENSG00000078295.17
31      AMOTL2 ENSG00000114019.15

```

## Code cell 58

```python
%%R 

sceptre_object_highmoi <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs_highmoi,
  side = "both", control_group = "complement")
```

## Code cell 59

```python
%%R

png(file.path(data_dir, "grna_count_distributions.png"), width = 3000, height = 2000, res = 300)

print(plot_grna_count_distributions(
  sceptre_object = sceptre_object_highmoi,
  n_grnas_to_plot = 8
))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "grna_count_distributions.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d1_sceptre_transfer/grna_count_distributions.png

```

## Code cell 60

```python
%%R

plot_covariates(sceptre_object_highmoi, p_mito_threshold = 0.25)
```

### Text output (text/plain)

```text
<IPython.core.display.Image object>
```

## Code cell 61

```python
%%R
sceptre_object_highmoi <- assign_grnas(
  sceptre_object = sceptre_object_highmoi,
  method = "mixture",
  parallel = FALSE
)
```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178488.23-P1P2 (1 of 1926)
Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178620.23-P1P2 (5 of 1926)
Performing gRNA-to-cell assignments for gRNA ACTL7B_-_111617995.23-P1P2 (10 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR1B_+_52345762.23-P1P2 (15 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2A_+_148602343.23-P1P2 (20 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_+_38495339.23-P1P2 (25 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_-_38495819.23-P1P2 (30 of 1926)
Performing gRNA-to-cell assignments for gRNA ADCY2_-_7396429.23-P1P2 (35 of 1926)
Performing gRNA-to-cell assignments for gRNA AMOTL2_+_134093279.23-P1P2 (40 of 1926)
Performing gRNA-to-cell assignments for gRNA ANP32E_+_150208342.23-P1P2 (45 of 1926)
Performing gRNA-to-cell assignments for gRNA APAF1_-_99039071.23-P1P2 (50 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_+_72766725.23-P1P2 (55 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_-_72766755.23-P1P2 (60 of 1926)
Performing gRNA-to-cell assignments for gRNA ARTN_-_44399429.23-P1 (65 of 1926)
Performing gRNA-to-cell assignments for gRNA ASCC3_-_101329193.23-P1P2 (70 of 1926)
Performing gRNA-to-cell assignments for gRNA ASH2L_-_37963059.23-P1P2 (75 of 1926)
Performing gRNA-to-cell assignments for gRNA ATMIN_+_81069661.23-P1P2 (80 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_+_111992112.23-P1P2 (85 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_-_111992129.23-P1P2 (90 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP6AP1_-_153657075.23-P1P2 (95 of 1926)
Performing gRNA-to-cell assignments for gRNA ATXN2L_-_28834668.23-P1P2 (100 of 1926)
Performing gRNA-to-cell assignments for gRNA BAD_+_64052102.23-P1P2 (105 of 1926)
Performing gRNA-to-cell assignments for gRNA BAG6_+_31620398.23-P1P2 (110 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_+_49458268.23-P1P2 (115 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_-_49458323.23-P1P2 (120 of 1926)
Performing gRNA-to-cell assignments for gRNA BCAS2_-_115124045.23-P1P2 (125 of 1926)
Performing gRNA-to-cell assignments for gRNA BMPR1A_-_88516386.23-P1P2 (130 of 1926)
Performing gRNA-to-cell assignments for gRNA BRD3_-_136932911.23-P1P2 (135 of 1926)
Performing gRNA-to-cell assignments for gRNA C17orf49_+_6918114.23-P2 (140 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_+_31532436.23-P1P2 (145 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_-_31532650.23-P1P2 (150 of 1926)
Performing gRNA-to-cell assignments for gRNA CACFD1_-_136325160.23-P1P2 (155 of 1926)
Performing gRNA-to-cell assignments for gRNA CACNG1_+_65040756.23-P1P2 (160 of 1926)
Performing gRNA-to-cell assignments for gRNA CAPN5_+_76778225.23-P1P2 (165 of 1926)
Performing gRNA-to-cell assignments for gRNA CASP3_+_185570560.23-P1P2 (170 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_+_107384793.23-P1P2 (175 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_-_107384635.23-P1P2 (180 of 1926)
Performing gRNA-to-cell assignments for gRNA CCM2_-_45039941.23-P1P2 (185 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNC_-_100016267.23-P1P2 (190 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNT2_-_135676416.23-P1P2 (195 of 1926)
Performing gRNA-to-cell assignments for gRNA CCT6B_+_33288381.23-P1P2 (200 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_+_30366413.23-P1P2 (205 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_-_30366630.23-P1P2 (210 of 1926)
Performing gRNA-to-cell assignments for gRNA CD99L2_-_150067120.23-P1P2 (215 of 1926)
Performing gRNA-to-cell assignments for gRNA CDH5_-_66400630.23-P1P2 (220 of 1926)
Performing gRNA-to-cell assignments for gRNA CDKN1A_-_36646538.23-P1P2 (225 of 1926)
Performing gRNA-to-cell assignments for gRNA CENPC_+_68411190.23-P1P2 (230 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_+_30773828.23-P1P2 (235 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_-_30773391.23-P1P2 (240 of 1926)
Performing gRNA-to-cell assignments for gRNA CHMP3_-_86790451.23-P1P2 (245 of 1926)
Performing gRNA-to-cell assignments for gRNA CHORDC1_-_89955894.23-P1P2 (250 of 1926)
Performing gRNA-to-cell assignments for gRNA CNEP1R1_+_50059284.23-P1P2 (255 of 1926)
Performing gRNA-to-cell assignments for gRNA CNOT2_+_70637651.23-P1P2 (260 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_+_83956339.23-P1P2 (265 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_-_83956572.23-P1P2 (270 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS6_-_99686676.23-P1P2 (275 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_-_31634088.23-P1P2 (280 of 1926)
Performing gRNA-to-cell assignments for gRNA CTR9_+_10772881.23-P1P2 (285 of 1926)
Performing gRNA-to-cell assignments for gRNA CUL3_-_225449667.23-P1P2 (290 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_+_145150267.23-P1P2 (295 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_-_145149996.23-P1P2 (300 of 1926)
Performing gRNA-to-cell assignments for gRNA CYCS_+_25164869.23-P1P2 (305 of 1926)
Performing gRNA-to-cell assignments for gRNA DACH1_-_72441101.23-P1P2 (310 of 1926)
Performing gRNA-to-cell assignments for gRNA DBR1_+_137893673.23-P1P2 (315 of 1926)
Performing gRNA-to-cell assignments for gRNA DCTN2_-_57940910.23-P1P2 (320 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_+_61851532.23-P1P2 (325 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_-_61851812.23-P1P2 (330 of 1926)
Performing gRNA-to-cell assignments for gRNA DGCR8_-_20068158.23-P1P2 (335 of 1926)
Performing gRNA-to-cell assignments for gRNA DNAJC17_-_41099215.23-P1P2 (340 of 1926)
Performing gRNA-to-cell assignments for gRNA DNM2_+_10829157.23-P1P2 (345 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM1_+_49575018.23-P1P2 (350 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_+_130700056.23-P1P2 (355 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_-_130700076.23-P1P2 (360 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM3_+_155113020.23-P1P2 (365 of 1926)
Performing gRNA-to-cell assignments for gRNA DROSHA_-_31531982.23-P1P2 (370 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_+_3985081.23-P1 (375 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_-_3985443.23-P1 (380 of 1926)
Performing gRNA-to-cell assignments for gRNA EHMT2_-_31865199.23-P1P2 (385 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3F_+_8009347.23-P1P2 (390 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3M_-_32605405.23-P1P2 (395 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4E_+_99849831.23-P1P2 (400 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184032364.23-P1P2 (405 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184033042.23-P1P2 (410 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G2_+_10830447.23-P1 (415 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF5A_-_7211217.23-P1P2 (420 of 1926)
Performing gRNA-to-cell assignments for gRNA TCEB3_-_24069975.23-P1P2 (425 of 1926)
Performing gRNA-to-cell assignments for gRNA EMC2_+_109456080.23-P1P2 (430 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823220.23-P1P2 (435 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823252.23-P1P2 (440 of 1926)
Performing gRNA-to-cell assignments for gRNA EOMES_-_27763342.23-P1P2 (445 of 1926)
Performing gRNA-to-cell assignments for gRNA ERCC2_+_45873802.23-P1P2 (450 of 1926)
Performing gRNA-to-cell assignments for gRNA ERO1L_+_53162425.23-P1P2 (455 of 1926)
Performing gRNA-to-cell assignments for gRNA ETS1_+_128392044.23-P1P2 (460 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_+_36134592.23-P1P2 (465 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_-_36134789.23-P1P2 (470 of 1926)
Performing gRNA-to-cell assignments for gRNA EWSR1_-_29664347.23-P1P2 (475 of 1926)
Performing gRNA-to-cell assignments for gRNA EXOSC8_-_37574979.23-P1P2 (480 of 1926)
Performing gRNA-to-cell assignments for gRNA FBXW11_-_171433300.23-P1P2 (485 of 1926)
Performing gRNA-to-cell assignments for gRNA FEN1_+_61560782.23-P1P2 (490 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_+_75532080.23-P1P2 (495 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_-_75532461.23-P1P2 (500 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT7_-_139926495.23-P1P2 (505 of 1926)
Performing gRNA-to-cell assignments for gRNA G3BP1_-_151151499.23-P1P2 (510 of 1926)
Performing gRNA-to-cell assignments for gRNA GCK_+_44229021.23-P1P2 (515 of 1926)
Performing gRNA-to-cell assignments for gRNA GOLT1B_+_21654777.23-P1P2 (520 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_+_45000529.23-P1P2 (525 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_-_45000948.23-P1P2 (530 of 1926)
Performing gRNA-to-cell assignments for gRNA GPI_+_34856284.23-P1P2 (535 of 1926)
Performing gRNA-to-cell assignments for gRNA GPR107_-_132816044.23-P1P2 (540 of 1926)
Performing gRNA-to-cell assignments for gRNA GSPT1_-_12009714.23-P1P2 (545 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A1_+_81687215.23-P1P2 (550 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_+_59949355.23-P1P2 (555 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_-_59949657.23-P1P2 (560 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2H4_-_30876185.23-P1P2 (565 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF3C2_+_27579604.23-P1P2 (570 of 1926)
Performing gRNA-to-cell assignments for gRNA HAND1_+_153857732.23-P1P2 (575 of 1926)
Performing gRNA-to-cell assignments for gRNA HDAC2_-_114292253.23-P1P2 (580 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_+_81399312.23-P1P2 (585 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_-_81399351.23-P1P2 (590 of 1926)
Performing gRNA-to-cell assignments for gRNA HINT2_+_35814995.23-P1P2 (595 of 1926)
Performing gRNA-to-cell assignments for gRNA HNF4A_+_43029971.23-P1P2 (600 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPM_-_8509885.23-P1 (605 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPU_+_245027753.23-P1P2 (610 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_+_43702292.23-P1P2 (615 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_-_43702305.23-P1P2 (620 of 1926)
Performing gRNA-to-cell assignments for gRNA HSPA5_-_128003534.23-P1P2 (625 of 1926)
Performing gRNA-to-cell assignments for gRNA HUWE1_-_53710540.23-P1P2 (630 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_-_30007588.23-P1P2 (635 of 1926)
Performing gRNA-to-cell assignments for gRNA INTS10_+_19675000.23-P1P2 (640 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_+_65903046.23-P1P2 (645 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_-_65903523.23-P1P2 (650 of 1926)
Performing gRNA-to-cell assignments for gRNA C7orf26_+_6630358.23-P1P2 (655 of 1926)
Performing gRNA-to-cell assignments for gRNA IPO8_+_30848771.23-P1P2 (660 of 1926)
Performing gRNA-to-cell assignments for gRNA ITGB1BP1_+_9563543.23-P1P2 (665 of 1926)
Performing gRNA-to-cell assignments for gRNA SSFA2_+_182756743.23-P1P2 (670 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_+_109914946.23-P1P2 (675 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_-_109915112.23-P1P2 (680 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM1A_-_23346409.23-P1P2 (685 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM3B_-_137688431.23-P1P2 (690 of 1926)
Performing gRNA-to-cell assignments for gRNA KDR_-_55991659.23-P1P2 (695 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA0232_+_6784651.23-P1P2 (700 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153706791.23-P1P2 (705 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153707248.23-P1P2 (710 of 1926)
Performing gRNA-to-cell assignments for gRNA LAMTOR4_-_99746921.23-P1P2 (715 of 1926)
Performing gRNA-to-cell assignments for gRNA LARP4_-_50794864.23-P1P2 (720 of 1926)
Performing gRNA-to-cell assignments for gRNA LATS2_-_21635336.23-P1P2 (725 of 1926)
Performing gRNA-to-cell assignments for gRNA LEMD3_+_65563441.23-P1P2 (730 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_+_54960428.23-P1P2 (735 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_-_54960616.23-P1P2 (740 of 1926)
Performing gRNA-to-cell assignments for gRNA LHFPL2_-_77944519.23-P1P2 (745 of 1926)
Performing gRNA-to-cell assignments for gRNA LSM12_-_42143603.23-P1P2 (750 of 1926)
Performing gRNA-to-cell assignments for gRNA MANBAL_+_35918136.23-P1P2 (755 of 1926)
Performing gRNA-to-cell assignments for gRNA MAU2_+_19431732.23-P1P2 (760 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_+_47807678.23-P1P2 (765 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_-_47807877.23-P1P2 (770 of 1926)
Performing gRNA-to-cell assignments for gRNA MEAF6_-_37980319.23-P1P2 (775 of 1926)
Performing gRNA-to-cell assignments for gRNA MED1_-_37607140.23-P1P2 (780 of 1926)
Performing gRNA-to-cell assignments for gRNA MED12_+_70338541.23-P1P2 (785 of 1926)
Performing gRNA-to-cell assignments for gRNA MED15_+_20862189.23-P1 (790 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_+_57479335.23-P1P2 (795 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_-_57479437.23-P1P2 (800 of 1926)
Performing gRNA-to-cell assignments for gRNA MEMO1_-_32264882.23-P1 (805 of 1926)
Performing gRNA-to-cell assignments for gRNA MIXL1_+_226411425.23-P1P2 (810 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL28_+_420487.23-P1P2 (815 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL32_+_42972018.23-P1P2 (820 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_+_179322181.23-P1P2 (825 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_-_179322213.23-P1P2 (830 of 1926)
Performing gRNA-to-cell assignments for gRNA LUZP6_-_135612750.23-P1P2 (835 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA20_-_19998023.23-P1P2 (840 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA25_+_112546538.23-P1P2 (845 of 1926)
Performing gRNA-to-cell assignments for gRNA NCBP1_+_100396082.23-P1P2 (850 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_+_172428630.23-P1P2 (855 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_-_172428841.23-P1P2 (860 of 1926)
Performing gRNA-to-cell assignments for gRNA NEDD8_-_24701363.23-P1P2 (865 of 1926)
Performing gRNA-to-cell assignments for gRNA NF2_+_29999635.23-P1P2 (870 of 1926)
Performing gRNA-to-cell assignments for gRNA NFIB_+_14314559.23-P1P2 (875 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_+_139418329.23-P2 (880 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139418349.23-P2 (885 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139440503.23-P1 (890 of 1926)
Performing gRNA-to-cell assignments for gRNA NUFIP2_+_27621073.23-P1P2 (895 of 1926)
Performing gRNA-to-cell assignments for gRNA NUP93_-_56764090.23-P1P2 (900 of 1926)
Performing gRNA-to-cell assignments for gRNA OGT_+_70752948.23-P1P2 (905 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_+_59383170.23-P1P2 (910 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_-_59383348.23-P1P2 (915 of 1926)
Performing gRNA-to-cell assignments for gRNA PAM16_-_4401236.23-P1P2 (920 of 1926)
Performing gRNA-to-cell assignments for gRNA PAXBP1_-_34143614.23-P1P2 (925 of 1926)
Performing gRNA-to-cell assignments for gRNA PBLD_-_70092436.23-P1 (930 of 1926)
Performing gRNA-to-cell assignments for gRNA PCBP2_+_53835614.23-P2 (935 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_+_94292752.23-P2 (940 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_-_94293027.23-P2 (945 of 1926)
Performing gRNA-to-cell assignments for gRNA PGRMC2_+_129208913.23-P1P2 (950 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_-_178866211.23-P1P2 (955 of 1926)
Performing gRNA-to-cell assignments for gRNA PITPNC1_-_65373433.23-P1P2 (960 of 1926)
Performing gRNA-to-cell assignments for gRNA PLCG1_+_39766176.23-P1P2 (965 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_+_151042273.23-P2 (970 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_-_151042422.23-P2 (975 of 1926)
Performing gRNA-to-cell assignments for gRNA PLRG1_-_155471154.23-P1P2 (980 of 1926)
Performing gRNA-to-cell assignments for gRNA PLXNA2_-_208417275.23-P1P2 (985 of 1926)
Performing gRNA-to-cell assignments for gRNA POLE4_-_75185725.23-P1P2 (990 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR2M_+_57998962.23-P1P2 (995 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940045.23-P1P2 (1000 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940414.23-P1P2 (1005 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3K_-_103144.23-P1P2 (1010 of 1926)
Performing gRNA-to-cell assignments for gRNA PPHLN1_-_42720007.23-P1P2 (1015 of 1926)
Performing gRNA-to-cell assignments for gRNA PPIL4_+_149867146.23-P1P2 (1020 of 1926)
Performing gRNA-to-cell assignments for gRNA PPP2R1A_+_52693317.23-P1 (1025 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_+_53195217.23-P1P2 (1030 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_-_53195332.23-P1P2 (1035 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKD1_-_30396996.23-P1P2 (1040 of 1926)
Performing gRNA-to-cell assignments for gRNA PRMT5_+_23398554.23-P1 (1045 of 1926)
Performing gRNA-to-cell assignments for gRNA C2orf53_-_27362230.23-P1P2 (1050 of 1926)
Performing gRNA-to-cell assignments for gRNA PRSS8_+_31147010.23-P1P2 (1055 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_+_36909062.23-P1P2 (1060 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_-_36909230.23-P1P2 (1065 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMD4_-_151227289.23-P1P2 (1070 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMG4_-_3259227.23-P1P2 (1075 of 1926)
Performing gRNA-to-cell assignments for gRNA PTAR1_+_72374835.23-P1P2 (1080 of 1926)
Performing gRNA-to-cell assignments for gRNA PTBP1_+_797463.23-P1P2 (1085 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_+_89623388.23-P1P2 (1090 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_-_89623526.23-P1P2 (1095 of 1926)
Performing gRNA-to-cell assignments for gRNA PTMS_-_6875591.23-P1P2 (1100 of 1926)
Performing gRNA-to-cell assignments for gRNA PYGO2_-_154934135.23-P1P2 (1105 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB18_-_27793292.23-P1P2 (1110 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB6A_+_73471983.23-P1P2 (1115 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_+_76251918.23-P1P2 (1120 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_-_76252008.23-P1P2 (1125 of 1926)
Performing gRNA-to-cell assignments for gRNA GNB2L1_-_180670840.23-P1P2 (1130 of 1926)
Performing gRNA-to-cell assignments for gRNA RAD21_-_117886635.23-P1P2 (1135 of 1926)
Performing gRNA-to-cell assignments for gRNA RANGRF_+_8192257.23-P1P2 (1140 of 1926)
Performing gRNA-to-cell assignments for gRNA RBBP5_+_205091050.23-P1P2 (1145 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM14-RBM4_+_66406280.23-P2 (1150 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM33_+_155437802.23-P1P2 (1155 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM4_+_66406319.23-P2 (1160 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_+_114271352.23-P1P2 (1165 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_-_114271700.23-P1P2 (1170 of 1926)
Performing gRNA-to-cell assignments for gRNA RBX1_-_41347852.23-P1P2 (1175 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1551_+_32112399.23-P1P2 (1180 of 1926)
Performing gRNA-to-cell assignments for gRNA RFX3_-_3525750.23-P1P2 (1185 of 1926)
Performing gRNA-to-cell assignments for gRNA RGP1_+_35749328.23-P1P2 (1190 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_+_95067102.23-P1P2 (1195 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_-_95067164.23-P1P2 (1200 of 1926)
Performing gRNA-to-cell assignments for gRNA RIBC1_-_53449606.23-P1P2 (1205 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1432_+_5629173.23-P1P2 (1210 of 1926)
Performing gRNA-to-cell assignments for gRNA RLF_+_40627522.23-P1P2 (1215 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF115_-_145611148.23-P1P2 (1220 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_+_104296142.23-P1P2 (1225 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_-_104296250.23-P1P2 (1230 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF40_-_30773681.23-P1P2 (1235 of 1926)
Performing gRNA-to-cell assignments for gRNA RPL37A_-_217363600.23-P1P2 (1240 of 1926)
Performing gRNA-to-cell assignments for gRNA RPS15A_-_18801590.23-P1P2 (1245 of 1926)
Performing gRNA-to-cell assignments for gRNA SAV1_+_51134992.23-P1P2 (1250 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_+_33104038.23-P1P2 (1255 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_-_33104343.23-P1P2 (1260 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF8_+_155054976.23-P1P2 (1265 of 1926)
Performing gRNA-to-cell assignments for gRNA EPT1_+_26569050.23-P1P2 (1270 of 1926)
Performing gRNA-to-cell assignments for gRNA SERF2_+_44069659.23-P2 (1275 of 1926)
Performing gRNA-to-cell assignments for gRNA SERPING1_+_57365137.23-P1P2 (1280 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_+_122242041.23-P1P2 (1285 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_-_122242294.23-P1P2 (1290 of 1926)
Performing gRNA-to-cell assignments for gRNA SETDB1_-_150898815.23-P1P2 (1295 of 1926)
Performing gRNA-to-cell assignments for gRNA SF1_-_64545937.23-P1P2 (1300 of 1926)
Performing gRNA-to-cell assignments for gRNA CCDC101_-_28565268.23-P1P2 (1305 of 1926)
Performing gRNA-to-cell assignments for gRNA SKA3_+_21750519.23-P1P2 (1310 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_+_78591497.23-P1P2 (1315 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_-_78591737.23-P1P2 (1320 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC29A1_-_44187438.23-P2 (1325 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC35C1_-_45826711.23-P1P2 (1330 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD1_+_146403428.23-P1P2 (1335 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD3_+_67358229.23-P1P2 (1340 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_+_48556623.23-P1P2 (1345 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_-_48556657.23-P1P2 (1350 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD5_-_135468694.23-P1P2 (1355 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA4_-_11071741.23-P1P2 (1360 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA5_+_144435243.23-P1P2 (1365 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCB1_+_24129228.23-P1 (1370 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_+_47822975.23-P1P2 (1375 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_-_47822893.23-P1P2 (1380 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCD3_-_150974212.23-P1 (1385 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCE1_-_38803669.23-P1P2 (1390 of 1926)
Performing gRNA-to-cell assignments for gRNA SMC2_-_106856631.23-P1P2 (1395 of 1926)
Performing gRNA-to-cell assignments for gRNA SMG5_+_156252537.23-P1P2 (1400 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_+_103852.23-P1P2 (1405 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_-_104242.23-P1P2 (1410 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP27_-_70121181.23-P1P2 (1415 of 1926)
Performing gRNA-to-cell assignments for gRNA SNUPN_-_75917960.23-P1P2 (1420 of 1926)
Performing gRNA-to-cell assignments for gRNA SOX4_+_21594213.23-P1 (1425 of 1926)
Performing gRNA-to-cell assignments for gRNA SPIN1_+_91003737.23-P1P2 (1430 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_+_2328680.23-P1P2 (1435 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_-_2328912.23-P1P2 (1440 of 1926)
Performing gRNA-to-cell assignments for gRNA SPTY2D1_-_18655846.23-P1P2 (1445 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP19_+_112197188.23-P1P2 (1450 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP54_-_35452244.23-P1P2 (1455 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP68_+_74068540.23-P1P2 (1460 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_+_126138730.23-P1P2 (1465 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_-_126138742.23-P1P2 (1470 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPRB_-_133524719.23-P1P2 (1475 of 1926)
Performing gRNA-to-cell assignments for gRNA SRRM2_-_2802714.23-P1P2 (1480 of 1926)
Performing gRNA-to-cell assignments for gRNA SRSF6_+_42086993.23-P1 (1485 of 1926)
Performing gRNA-to-cell assignments for gRNA SSB_+_170655450.23-P1P2 (1490 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_+_57102902.23-P1P2 (1495 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_-_57103355.23-P1P2 (1500 of 1926)
Performing gRNA-to-cell assignments for gRNA STAG2_-_123095658.23-P1 (1505 of 1926)
Performing gRNA-to-cell assignments for gRNA STK40_+_36851507.23-P1P2 (1510 of 1926)
Performing gRNA-to-cell assignments for gRNA SUPT20H_+_37633798.23-P1P2 (1515 of 1926)
Performing gRNA-to-cell assignments for gRNA SUZ12_+_30264157.23-P1P2 (1520 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_+_64901973.23-P1P2 (1525 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_-_64901981.23-P1P2 (1530 of 1926)
Performing gRNA-to-cell assignments for gRNA TADA2B_-_7045121.23-P1P2 (1535 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF10_+_6633406.23-P1P2 (1540 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF2_+_120845032.23-P1P2 (1545 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF3_-_7860540.23-P1P2 (1550 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_+_60640814.23-P1P2 (1555 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_-_60640498.23-P1P2 (1560 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF6L_-_62538954.23-P1P2 (1565 of 1926)
Performing gRNA-to-cell assignments for gRNA TAOK1_-_27717907.23-P1P2 (1570 of 1926)
Performing gRNA-to-cell assignments for gRNA FAM208A_+_56717273.23-P1P2 (1575 of 1926)
Performing gRNA-to-cell assignments for gRNA TBC1D9_+_141677481.23-P1P2 (1580 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_+_170863524.23-P1P2 (1585 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_-_170863481.23-P1P2 (1590 of 1926)
Performing gRNA-to-cell assignments for gRNA TBX3_-_115121998.23-P1P2 (1595 of 1926)
Performing gRNA-to-cell assignments for gRNA T_-_166581333.23-P1P2 (1600 of 1926)
Performing gRNA-to-cell assignments for gRNA TGFBR2_+_30648526.23-P1P2 (1605 of 1926)
Performing gRNA-to-cell assignments for gRNA THBD_+_23029196.23-P2 (1610 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_+_2977592.23-P1P2 (1615 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_-_2977643.23-P1P2 (1620 of 1926)
Performing gRNA-to-cell assignments for gRNA MESDC1_-_81293340.23-P1P2 (1625 of 1926)
Performing gRNA-to-cell assignments for gRNA TMEM41B_+_9336156.23-P1P2 (1630 of 1926)
Performing gRNA-to-cell assignments for gRNA TMTC3_-_88536127.23-P1P2 (1635 of 1926)
Performing gRNA-to-cell assignments for gRNA TOP1_+_39657538.23-P1P2 (1640 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_+_7590531.23-P1P2 (1645 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_-_7590704.23-P1P2 (1650 of 1926)
Performing gRNA-to-cell assignments for gRNA TRIM24_-_138145308.23-P1P2 (1655 of 1926)
Performing gRNA-to-cell assignments for gRNA TRRAP_-_98476465.23-P1P2 (1660 of 1926)
Performing gRNA-to-cell assignments for gRNA TSTD2_+_100395594.23-P1P2 (1665 of 1926)
Performing gRNA-to-cell assignments for gRNA UBE2E1_+_23847541.23-P1P2 (1670 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_+_38924011.23-P1P2 (1675 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_-_38924421.23-P1P2 (1680 of 1926)
Performing gRNA-to-cell assignments for gRNA UPF2_+_12084799.23-P1P2 (1685 of 1926)
Performing gRNA-to-cell assignments for gRNA UXT_-_47518505.23-P1P2 (1690 of 1926)
Performing gRNA-to-cell assignments for gRNA VHL_-_10183525.23-P1P2 (1695 of 1926)
Performing gRNA-to-cell assignments for gRNA VPS54_+_64246104.23-P1P2 (1700 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_+_28821698.23-P1P2 (1705 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_-_28822074.23-P1P2 (1710 of 1926)
Performing gRNA-to-cell assignments for gRNA WAPAL_+_88281533.23-P1P2 (1715 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR12_-_203776383.23-P1P2 (1720 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR77_+_111991860.23-P1P2 (1725 of 1926)
Performing gRNA-to-cell assignments for gRNA WIZ_-_15543016.23-P1P2 (1730 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_+_9983328.23-P1P2 (1735 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_-_9983510.23-P1P2 (1740 of 1926)
Performing gRNA-to-cell assignments for gRNA XPO5_-_43543670.23-P1P2 (1745 of 1926)
Performing gRNA-to-cell assignments for gRNA XRN1_-_142166648.23-P1P2 (1750 of 1926)
Performing gRNA-to-cell assignments for gRNA YARS2_+_32908787.23-P1P2 (1755 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB22_+_33285454.23-P1P2 (1760 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_+_6640093.23-P1P2 (1765 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_-_6640138.23-P1P2 (1770 of 1926)
Performing gRNA-to-cell assignments for gRNA ZC3H8_-_113012435.23-P1P2 (1775 of 1926)
Performing gRNA-to-cell assignments for gRNA ZNF281_+_200379151.23-P1P2 (1780 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00117 (1785 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00233 (1790 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00518 (1795 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00656 (1800 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01226 (1805 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01730 (1810 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02009 (1815 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02288 (1820 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02777 (1825 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03026 (1830 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03252 (1835 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03620 (1840 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_SGTN_04 (1845 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_SGTN_25 (1850 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_SGTN_26 (1855 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1044 (1860 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1414 (1865 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1554 (1870 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1785 (1875 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2009 (1880 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2271 (1885 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2425 (1890 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2614 (1895 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_293 (1900 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_541 (1905 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_725 (1910 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_860 (1915 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_SGTN_14 (1920 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_SGTN_01 (1925 of 1926)

```

## Code cell 62

```python
%%R
library(Matrix)

# 1. Get the assignment matrix
A <- get_grna_assignments(sceptre_object_highmoi)

# 2. Ensure it is a numeric sparse matrix (1/0)
if (inherits(A, "Matrix") || is.matrix(A)) {
    A_num <- if (inherits(A, "Matrix")) drop0(A * 1) else (A * 1)
} else {
    stop("Output is not a matrix.")
}

# --- FIX: ROBUST NAME RECOVERY ---
# We retrieve names using Sceptre's standard metadata accessors
# instead of touching internal slots like @grna_matrix.

# Try to get cell names from the matrix first
current_cell_names <- colnames(A_num)
current_guide_names <- rownames(A_num)

# If missing, grab them from the main object's covariate data frame
if (is.null(current_cell_names)) {
    message("⚠️ Matrix missing colnames. Fetching from sceptre_object...")
    # The covariate data frame ALWAYS has the cell IDs as rownames
    current_cell_names <- rownames(sceptre_object_highmoi@covariate_data_frame)
    
    # Safety check: lengths must match
    if (length(current_cell_names) != ncol(A_num)) {
        stop(paste("Mismatch! Matrix has", ncol(A_num), "cols but object has", length(current_cell_names), "cells."))
    }
    colnames(A_num) <- current_cell_names
}

# If guide names are missing, they should match the unique values in your guide mapping
if (is.null(current_guide_names)) {
    message("⚠️ Matrix missing rownames. Fetching from guide mapping...")
    # Get unique targets from the mapping table you loaded earlier
    # Note: This assumes the matrix rows are sorted/aligned to these unique targets.
    # In 'mixture' method, rows usually correspond to the 'grna_target' labels.
    
    # SAFE FALLBACK: If we can't be 100% sure of order, we warn the user.
    # But usually, they match the 'grna_target' column from your input.
    unique_targets <- sort(unique(guide_mapping$grna_target))
    
    if (length(unique_targets) == nrow(A_num)) {
         rownames(A_num) <- unique_targets
    } else {
         warning("Could not automatically align target names. Using generic IDs.")
         rownames(A_num) <- paste0("target_", seq_len(nrow(A_num)))
    }
}
# ---------------------------------

# 3. Convert to Long Format
s <- summary(A_num)
assign_long <- data.frame(
    cell_id = colnames(A_num)[s$j],
    grna_target = rownames(A_num)[s$i],
    value = s$x
)

# 4. Save
out_path <- file.path(data_dir, "grna_assignments_long.csv")
write.csv(assign_long, out_path, row.names = FALSE)

message("✅ Saved long CSV to: ", out_path)
print(head(assign_long))
```

### Text output (stream)

```text
                     cell_id                grna_target value
1 AACGTGATGTGTTCTATACGGTTA_2 CNEP1R1_+_50059217.23-P1P2     1
2 AACGTGATGTGTTCTATACGGTTA_2 SPTY2D1_-_18655660.23-P1P2     1
3 AACGTGATGTGTTCTATACGGTTA_2   TRRAP_-_98476562.23-P1P2     1
4 CCATCCTCCAACCACATTGGGAGA_3 ATP5F1_+_111992171.23-P1P2     1
5 CCATCCTCCAACCACATTGGGAGA_3  C7orf26_+_6629634.23-P1P2     1
6 CCATCCTCCAACCACATTGGGAGA_3   PLCG1_+_39766234.23-P1P2     1

```

### Text output (text/plain)

```text
⚠️ Matrix missing colnames. Fetching from sceptre_object...
✅ Saved long CSV to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d1_sceptre_transfer/grna_assignments_long.csv

```

## Code cell 63

```python
%%R
png(file.path(data_dir, "MOI.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "MOI.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d1_sceptre_transfer/MOI.png

```

## Code cell 64

```python
%%R
sceptre_object_highmoi <- run_qc(
  sceptre_object = sceptre_object_highmoi,
  p_mito_threshold = 0.25
)

png(file.path(data_dir, "QC.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "QC.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d1_sceptre_transfer/QC.png

```

## Code cell 65

```python
%%R
# Calibration check


sceptre_object_highmoi <- run_calibration_check(
  sceptre_object = sceptre_object_highmoi,
  parallel = FALSE
)

png(file.path(data_dir, "calibration_check.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "calibration_check.png"))

```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Constructing negative control pairs. ✓
Running precomputation on response ENSG00000197747.9 (1 of 242)
Running precomputation on response ENSG00000103035.11 (5 of 242)
Running precomputation on response ENSG00000077150.21 (10 of 242)
Running precomputation on response ENSG00000175198.17 (15 of 242)
Running precomputation on response ENSG00000236859.8 (20 of 242)
Running precomputation on response ENSG00000176171.12 (25 of 242)
Running precomputation on response ENSG00000162733.19 (30 of 242)
Running precomputation on response ENSG00000269547.1 (35 of 242)
Running precomputation on response ENSG00000145996.11 (40 of 242)
Running precomputation on response ENSG00000116141.17 (45 of 242)
Running precomputation on response ENSG00000166025.18 (50 of 242)
Running precomputation on response ENSG00000089775.12 (55 of 242)
Running precomputation on response ENSG00000099204.22 (60 of 242)
Running precomputation on response ENSG00000157224.16 (65 of 242)
Running precomputation on response ENSG00000173207.13 (70 of 242)
Running precomputation on response ENSG00000008869.12 (75 of 242)
Running precomputation on response ENSG00000111145.8 (80 of 242)
Running precomputation on response ENSG00000137815.15 (85 of 242)
Running precomputation on response ENSG00000213614.11 (90 of 242)
Running precomputation on response ENSG00000048342.18 (95 of 242)
Running precomputation on response ENSG00000197302.11 (100 of 242)
Running precomputation on response ENSG00000065243.20 (105 of 242)
Running precomputation on response ENSG00000080546.14 (110 of 242)
Running precomputation on response ENSG00000138835.24 (115 of 242)
Running precomputation on response ENSG00000118971.9 (120 of 242)
Running precomputation on response ENSG00000182287.15 (125 of 242)
Running precomputation on response ENSG00000143786.8 (130 of 242)
Running precomputation on response ENSG00000034713.8 (135 of 242)
Running precomputation on response ENSG00000133874.2 (140 of 242)
Running precomputation on response ENSG00000101413.12 (145 of 242)
Running precomputation on response ENSG00000269609.6 (150 of 242)
Running precomputation on response ENSG00000067829.20 (155 of 242)
Running precomputation on response ENSG00000269955.2 (160 of 242)
Running precomputation on response ENSG00000121897.15 (165 of 242)
Running precomputation on response ENSG00000110651.13 (170 of 242)
Running precomputation on response ENSG00000030419.17 (175 of 242)
Running precomputation on response ENSG00000136866.14 (180 of 242)
Running precomputation on response ENSG00000181830.9 (185 of 242)
Running precomputation on response ENSG00000176681.15 (190 of 242)
Running precomputation on response ENSG00000148730.7 (195 of 242)
Running precomputation on response ENSG00000069869.17 (200 of 242)
Running precomputation on response ENSG00000109670.17 (205 of 242)
Running precomputation on response ENSG00000107521.21 (210 of 242)
Running precomputation on response ENSG00000086015.23 (215 of 242)
Running precomputation on response ENSG00000112343.11 (220 of 242)
Running precomputation on response ENSG00000181027.11 (225 of 242)
Running precomputation on response ENSG00000171314.9 (230 of 242)
Running precomputation on response ENSG00000044090.13 (235 of 242)
Running precomputation on response ENSG00000143727.16 (240 of 242)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00259&non-targeting_02159&non-targeting_02571&non-targeting_03029&non-targeting_03081 (1 of 92)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00355&non-targeting_02025&non-targeting_02089&non-targeting_03081&non-targeting_03384 (5 of 92)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_01730&non-targeting_02193&non-targeting_03180&non-targeting_03620&non-targeting_03664 (10 of 92)
Analyzing pairs containing gRNA group non-targeting_00043&non-targeting_00184&non-targeting_02677&non-targeting_02777&non-targeting_02814&non-targeting_03427 (15 of 92)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00355&non-targeting_00390&non-targeting_00918&non-targeting_01509&non-targeting_01599 (20 of 92)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00972&non-targeting_02193&non-targeting_02574&non-targeting_02777&non-targeting_02962 (25 of 92)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00355&non-targeting_01601&non-targeting_03026&non-targeting_03081&non-targeting_03158 (30 of 92)
Analyzing pairs containing gRNA group non-targeting_00231&non-targeting_00518&non-targeting_00629&non-targeting_00794&non-targeting_02777&non-targeting_02837 (35 of 92)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00390&non-targeting_01601&non-targeting_02009&non-targeting_02288&non-targeting_03029 (40 of 92)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00231&non-targeting_00324&non-targeting_02025&non-targeting_03252&non-targeting_03384 (45 of 92)
Analyzing pairs containing gRNA group non-targeting_00043&non-targeting_00213&non-targeting_00259&non-targeting_00355&non-targeting_00614&non-targeting_01934 (50 of 92)
Analyzing pairs containing gRNA group non-targeting_00629&non-targeting_01185&non-targeting_01732&non-targeting_02089&non-targeting_02677&non-targeting_03483 (55 of 92)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00518&non-targeting_00614&non-targeting_02025&non-targeting_03342&non-targeting_03427 (60 of 92)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00231&non-targeting_00656&non-targeting_02025&non-targeting_02288&non-targeting_02777 (65 of 92)
Analyzing pairs containing gRNA group non-targeting_00614&non-targeting_00972&non-targeting_02911&non-targeting_03081&non-targeting_03180&non-targeting_03384 (70 of 92)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_00555&non-targeting_01934&non-targeting_02193&non-targeting_02677&non-targeting_03252 (75 of 92)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00656&non-targeting_00918&non-targeting_01934&non-targeting_02574&non-targeting_02814 (80 of 92)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_00355&non-targeting_02837&non-targeting_02911&non-targeting_03342&non-targeting_03620 (85 of 92)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00518&non-targeting_00614&non-targeting_01226&non-targeting_02025&non-targeting_02837 (90 of 92)

```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d1_sceptre_transfer/calibration_check.png
In addition: There were 50 or more warnings (use warnings() to see the first 50)

```

## Code cell 66

```python
%%R
calibration_result <- get_result(
  sceptre_object = sceptre_object_highmoi,
  analysis = "run_calibration_check"
)
head(calibration_result)
```

### Text output (stream)

```text
          response_id
               <char>
1: ENSG00000121957.15
2: ENSG00000115414.21
3:  ENSG00000122692.9
4: ENSG00000150995.21
5: ENSG00000048471.14
6: ENSG00000149262.19
                                                                                                               grna_target
                                                                                                                    <char>
1: non-targeting_00259&non-targeting_01934&non-targeting_02677&non-targeting_03081&non-targeting_03252&non-targeting_03664
2: non-targeting_00355&non-targeting_00518&non-targeting_00555&non-targeting_00794&non-targeting_01185&non-targeting_02159
3: non-targeting_00184&non-targeting_00259&non-targeting_02159&non-targeting_02571&non-targeting_03029&non-targeting_03081
4: non-targeting_01504&non-targeting_01601&non-targeting_01732&non-targeting_01812&non-targeting_02962&non-targeting_03483
5: non-targeting_00231&non-targeting_01730&non-targeting_02089&non-targeting_03180&non-targeting_03342&non-targeting_03384
6: non-targeting_00043&non-targeting_00117&non-targeting_00555&non-targeting_00918&non-targeting_03252&non-targeting_03427
   n_nonzero_trt n_nonzero_cntrl pass_qc    p_value log_2_fold_change
           <int>           <int>  <lgcl>      <num>             <num>
1:            31            1118    TRUE 0.02530799         0.5678968
2:            59            4094    TRUE 0.02962683         0.3888507
3:            81            3788    TRUE 0.03600000        -0.3017934
4:            16            1773    TRUE 0.03600000        -0.6002785
5:            41            1598    TRUE 0.04000000         0.3800087
6:            38            1655    TRUE 0.04000000         0.4606734
   significant
        <lgcl>
1:       FALSE
2:       FALSE
3:       FALSE
4:       FALSE
5:       FALSE
6:       FALSE

```

## Code cell 67

```python
%%R
sceptre_object_highmoi <- run_discovery_analysis(
  sceptre_object = sceptre_object_highmoi,
  parallel = FALSE
)

png(file.path(data_dir, "discovery.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "discovery.png"))
```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Running precomputation on response ENSG00000004487.18 (1 of 243)
Running precomputation on response ENSG00000117000.9 (5 of 243)
Running precomputation on response ENSG00000265491.5 (10 of 243)
Running precomputation on response ENSG00000198952.8 (15 of 243)
Running precomputation on response ENSG00000153187.21 (20 of 243)
Running precomputation on response ENSG00000162961.14 (25 of 243)
Running precomputation on response ENSG00000121989.15 (30 of 243)
Running precomputation on response ENSG00000036257.14 (35 of 243)
Running precomputation on response ENSG00000114739.14 (40 of 243)
Running precomputation on response ENSG00000114127.11 (45 of 243)
Running precomputation on response ENSG00000128052.10 (50 of 243)
Running precomputation on response ENSG00000153147.6 (55 of 243)
Running precomputation on response ENSG00000113360.17 (60 of 243)
Running precomputation on response ENSG00000113658.18 (65 of 243)
Running precomputation on response ENSG00000180822.12 (70 of 243)
Running precomputation on response ENSG00000124571.18 (75 of 243)
Running precomputation on response ENSG00000131013.4 (80 of 243)
Running precomputation on response ENSG00000106591.4 (85 of 243)
Running precomputation on response ENSG00000188186.11 (90 of 243)
Running precomputation on response ENSG00000104613.12 (95 of 243)
Running precomputation on response ENSG00000179091.5 (100 of 243)
Running precomputation on response ENSG00000106723.17 (105 of 243)
Running precomputation on response ENSG00000136908.18 (110 of 243)
Running precomputation on response ENSG00000151461.20 (115 of 243)
Running precomputation on response ENSG00000107779.14 (120 of 243)
Running precomputation on response ENSG00000198730.9 (125 of 243)
Running precomputation on response ENSG00000149136.9 (130 of 243)
Running precomputation on response ENSG00000162298.19 (135 of 243)
Running precomputation on response ENSG00000076053.12 (140 of 243)
Running precomputation on response ENSG00000133704.10 (145 of 243)
Running precomputation on response ENSG00000197111.16 (150 of 243)
Running precomputation on response ENSG00000110906.13 (155 of 243)
Running precomputation on response ENSG00000165480.16 (160 of 243)
Running precomputation on response ENSG00000100462.16 (165 of 243)
Running precomputation on response ENSG00000165417.12 (170 of 243)
Running precomputation on response ENSG00000166949.17 (175 of 243)
Running precomputation on response ENSG00000161981.11 (180 of 243)
Running precomputation on response ENSG00000134419.15 (185 of 243)
Running precomputation on response ENSG00000052344.16 (190 of 243)
Running precomputation on response ENSG00000141510.18 (195 of 243)
Running precomputation on response ENSG00000125686.12 (200 of 243)
Running precomputation on response ENSG00000154217.16 (205 of 243)
Running precomputation on response ENSG00000005206.17 (210 of 243)
Running precomputation on response ENSG00000011451.21 (215 of 243)
Running precomputation on response ENSG00000105568.19 (220 of 243)
Running precomputation on response ENSG00000124193.16 (225 of 243)
Running precomputation on response ENSG00000128191.17 (230 of 243)
Running precomputation on response ENSG00000100387.9 (235 of 243)
Running precomputation on response ENSG00000147162.15 (240 of 243)
Analyzing pairs containing gRNA group KDM1A (1 of 243)
Analyzing pairs containing gRNA group RLF (5 of 243)
Analyzing pairs containing gRNA group RNF115 (10 of 243)
Analyzing pairs containing gRNA group SMG5 (15 of 243)
Analyzing pairs containing gRNA group HNRNPU (20 of 243)
Analyzing pairs containing gRNA group DPY30 (25 of 243)
Analyzing pairs containing gRNA group ACVR2A (30 of 243)
Analyzing pairs containing gRNA group CUL3 (35 of 243)
Analyzing pairs containing gRNA group ACVR2B (40 of 243)
Analyzing pairs containing gRNA group XRN1 (45 of 243)
Analyzing pairs containing gRNA group KDR (50 of 243)
Analyzing pairs containing gRNA group SMARCA5 (55 of 243)
Analyzing pairs containing gRNA group DROSHA (60 of 243)
Analyzing pairs containing gRNA group SMAD5 (65 of 243)
Analyzing pairs containing gRNA group PSMG4 (70 of 243)
Analyzing pairs containing gRNA group XPO5 (75 of 243)
Analyzing pairs containing gRNA group PPIL4 (80 of 243)
Analyzing pairs containing gRNA group MRPL32 (85 of 243)
Analyzing pairs containing gRNA group LAMTOR4 (90 of 243)
Analyzing pairs containing gRNA group INTS10 (95 of 243)
Analyzing pairs containing gRNA group CYC1 (100 of 243)
Analyzing pairs containing gRNA group SPIN1 (105 of 243)
Analyzing pairs containing gRNA group DPM2 (110 of 243)
Analyzing pairs containing gRNA group UPF2 (115 of 243)
Analyzing pairs containing gRNA group BMPR1A (120 of 243)
Analyzing pairs containing gRNA group CTR9 (125 of 243)
Analyzing pairs containing gRNA group SSRP1 (130 of 243)
Analyzing pairs containing gRNA group SYVN1 (135 of 243)
Analyzing pairs containing gRNA group RBM7 (140 of 243)
Analyzing pairs containing gRNA group IPO8 (145 of 243)
Analyzing pairs containing gRNA group PCBP2 (150 of 243)
Analyzing pairs containing gRNA group KCTD10 (155 of 243)
Analyzing pairs containing gRNA group SKA3 (160 of 243)
Analyzing pairs containing gRNA group PRMT5 (165 of 243)
Analyzing pairs containing gRNA group GTF2A1 (170 of 243)
Analyzing pairs containing gRNA group SMAD3 (175 of 243)
Analyzing pairs containing gRNA group SNRNP25 (180 of 243)
Analyzing pairs containing gRNA group RPS15A (185 of 243)
Analyzing pairs containing gRNA group PRSS8 (190 of 243)
Analyzing pairs containing gRNA group TP53 (195 of 243)
Analyzing pairs containing gRNA group MED1 (200 of 243)
Analyzing pairs containing gRNA group PITPNC1 (205 of 243)
Analyzing pairs containing gRNA group SPPL2B (210 of 243)
Analyzing pairs containing gRNA group WIZ (215 of 243)
Analyzing pairs containing gRNA group PPP2R1A (220 of 243)
Analyzing pairs containing gRNA group SRSF6 (225 of 243)
Analyzing pairs containing gRNA group DGCR8 (230 of 243)
Analyzing pairs containing gRNA group RBX1 (235 of 243)
Analyzing pairs containing gRNA group OGT (240 of 243)

```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d1_sceptre_transfer/discovery.png
In addition: There were 50 or more warnings (use warnings() to see the first 50)

```

## Code cell 68

```python
%%R
discovery_result <- get_result(
  sceptre_object = sceptre_object_highmoi,
  analysis = "run_discovery_analysis"
)

write.csv(discovery_result, file.path(data_dir, "TargetKD.csv"), row.names = FALSE)
message("✅ Saved discovery to: ", file.path(data_dir, "TargetKD.csv"))

```

### Text output (text/plain)

```text
✅ Saved discovery to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d1_sceptre_transfer/TargetKD.csv

```

## Markdown cell 69

### Day 2

## Code cell 70

```python
d2_mdata = mu.read(f"{data_dir}/d2_normalized_gex_and_guide.h5mu")
mdata = d2_mdata

# Filter out cells that have no guides
#mdata = mdata[mdata['GEX'].obs['guides_per_cell'] > 0]
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

## Code cell 71

```python
# Define a temporary folder for transfer
temp_dir = working_dir + "/d2_sceptre_transfer"
os.makedirs(temp_dir, exist_ok=True)

# 1. EXTRACT MATRICES (Transpose them!)
# Scanpy is (Cells × Genes), Sceptre needs (Genes × Cells)
print("Writing sparse matrices...")
sio.mmwrite(f"{temp_dir}/response_matrix.mtx", mdata['GEX'].layers['raw_counts'].T)
sio.mmwrite(f"{temp_dir}/grna_matrix.mtx", mdata['guide'].X.T)

# 2. EXTRACT NAMES (Row/Col names are lost in .mtx)
# Save Gene names and Guide names to assign them later in R
pd.Series(mdata['GEX'].var_names).to_csv(f"{temp_dir}/gene_names.csv", index=False, header=False)
pd.Series(mdata['guide'].var_names).to_csv(f"{temp_dir}/guide_names.csv", index=False, header=False)
pd.Series(mdata['GEX'].obs_names).to_csv(f"{temp_dir}/cell_names.csv", index=False, header=False)
```

### Text output (stream)

```text
Writing sparse matrices...

```

## Code cell 72

```python
mdata['GEX']
```

### Text output (text/plain)

```text
AnnData object with n_obs × n_vars = 21980 × 12732
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0'
    var: 'gene_id', 'symbol', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'n_counts', 'highly_deviant', 'highly_deviant_2k', 'binomial_deviance', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'orig_highly_variable', 'highly_variable_2k'
    uns: 'hvg', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'neighbors', 'pca', 'umap'
    obsm: 'X_pca', 'X_umap', 'guide_assignment', 'guide_counts'
    varm: 'PCs'
    layers: 'ambiguous', 'analytic_pearson_residuals', 'log1p_norm', 'mature', 'nascent', 'raw_counts', 'scran_normalization'
    obsp: 'connectivities', 'distances'
```

## Code cell 73

```python
# 3. EXTRACT COVARIATES 
rna = mdata['GEX']
covariates_df = pd.DataFrame({
    "batch": rna.obs['sample_id'],
    "grna_n_nonzero": rna.obs['guides_per_cell'],
    "grna_n_umis": rna.obs['guide_umi_counts'],
    "response_n_nonzero": rna.obs['n_genes_by_counts'],
    "response_n_umis": rna.obs['total_counts'],
    "response_p_mito": rna.obs['pct_counts_mt'],
})
covariates_df.to_csv(f"{temp_dir}/covariates.csv")
covariates_df
```

### Text output (text/plain)

```text
                                             batch  grna_n_nonzero  \
barcode                                                              
TAGGATGAATTGAGGATGCTTGGG_2  300genes:eTN5_GEX_sub2               1   
TAGGATGAAGCCATGCCGCCCGGA_2  300genes:eTN5_GEX_sub2               2   
CCAGTTCACAATGGAATGCTTGGG_3  300genes:eTN5_GEX_sub3               5   
ATCATTCCAGTCACTACGCCCGGA_2  300genes:eTN5_GEX_sub2               1   
GCGAGTAAGCTAACGACACAATTG_3  300genes:eTN5_GEX_sub3              10   
...                                            ...             ...   
CCGTGAGAACGCTCGACTTTGGTC_3  300genes:eTN5_GEX_sub3               5   
CACCTTACACAGCAGACACAATTG_1  300genes:eTN5_GEX_sub1              17   
CAAGACTAGCCAAGACTAAATATC_3  300genes:eTN5_GEX_sub3               2   
CAATGGAAGCTCGGTAGTGCTAGC_2  300genes:eTN5_GEX_sub2               6   
CCAGTTCACTGAGCCATAAATATC_2  300genes:eTN5_GEX_sub2               9   

                            grna_n_umis  response_n_nonzero  response_n_umis  \
barcode                                                                        
TAGGATGAATTGAGGATGCTTGGG_2         98.0                1276      2078.232178   
TAGGATGAAGCCATGCCGCCCGGA_2        142.0                 561       849.433411   
CCAGTTCACAATGGAATGCTTGGG_3        245.0                1141      1890.696289   
ATCATTCCAGTCACTACGCCCGGA_2         18.0                 407       549.833374   
GCGAGTAAGCTAACGACACAATTG_3        198.0                 680      1287.249268   
...                                 ...                 ...              ...   
CCGTGAGAACGCTCGACTTTGGTC_3        165.0                1783      3590.150391   
CACCTTACACAGCAGACACAATTG_1       1121.0                1462      2993.500000   
CAAGACTAGCCAAGACTAAATATC_3        134.0                1345      2546.260254   
CAATGGAAGCTCGGTAGTGCTAGC_2        105.0                 714      1144.333252   
CCAGTTCACTGAGCCATAAATATC_2        139.0                 879      1110.250244   

                            response_p_mito  
barcode                                      
TAGGATGAATTGAGGATGCTTGGG_2        13.304576  
TAGGATGAAGCCATGCCGCCCGGA_2        17.541105  
CCAGTTCACAATGGAATGCTTGGG_3        15.073812  
ATCATTCCAGTCACTACGCCCGGA_2        16.914215  
GCGAGTAAGCTAACGACACAATTG_3        22.101393  
...                                     ...  
CCGTGAGAACGCTCGACTTTGGTC_3        14.275167  
CACCTTACACAGCAGACACAATTG_1        23.300484  
CAAGACTAGCCAAGACTAAATATC_3         9.209585  
CAATGGAAGCTCGGTAGTGCTAGC_2        19.181475  
CCAGTTCACTGAGCCATAAATATC_2        14.636339  

[21980 rows x 6 columns]
```

## Code cell 74

```python
guide_meta_df = pd.read_csv(f"{working_dir}/300genes_guide_metadata_v43.tsv", sep='\t')
```

## Code cell 75

```python
# 4. PREPARE GUIDE MAPPING 
# Sceptre needs a table mapping 'grna_id' to 'grna_target'.

# 1. Create a "skeleton" DataFrame from the actual matrix order
# This ensures the order matches mdata['guide'].X perfectly
# We name the index 'grna_id' to match Sceptre's expectation immediately
guide_mapping = pd.DataFrame(index=mdata['guide'].var_names)
guide_mapping.index.name = 'grna_id'

# 2. Prepare metadata for merging
# Ensure 'guide_id' is the index so we can map efficiently
lookup_df = guide_meta_df.set_index('guide_id')

# 3. Map the targets onto the skeleton
# We use .loc[] to pull the specific target for each guide in the exact order required
guide_mapping['grna_target'] = lookup_df.loc[guide_mapping.index, 'intended_target_name']

# 4. Reset index so 'grna_id' becomes a column (Sceptre requires this)
guide_mapping = guide_mapping.reset_index()

# ==========================================
# VALIDATION CHECKS 
# ==========================================

# Check 1: Are there any NAs? (Guides in mdata that were missing from meta_df)
na_count = guide_mapping['grna_target'].isna().sum()
if na_count > 0:
    print(f"⚠️ WARNING: {na_count} guides have NaN targets! These will crash Sceptre.")
else:
    print("✅ All guides mapped successfully.")

# Check 2: Verify the first few rows manually
print("\nFirst 5 rows of mapping (Must match mdata['guide'].var_names order):")
print(guide_mapping.head())

# 5. Save for Sceptre
guide_mapping.to_csv(f"{temp_dir}/guide_mapping.csv", index=False)
```

### Text output (stream)

```text
✅ All guides mapped successfully.

First 5 rows of mapping (Must match mdata['guide'].var_names order):
                    grna_id grna_target
0  ACAA1_+_38178488.23-P1P2       ACAA1
1  ACAA1_+_38178517.23-P1P2       ACAA1
2  ACAA1_+_38178559.23-P1P2       ACAA1
3  ACAA1_+_38178570.23-P1P2       ACAA1
4  ACAA1_+_38178620.23-P1P2       ACAA1

```

## Code cell 76

```python
guide_mapping
```

### Text output (text/plain)

```text
                       grna_id grna_target
0     ACAA1_+_38178488.23-P1P2       ACAA1
1     ACAA1_+_38178517.23-P1P2       ACAA1
2     ACAA1_+_38178559.23-P1P2       ACAA1
3     ACAA1_+_38178570.23-P1P2       ACAA1
4     ACAA1_+_38178620.23-P1P2       ACAA1
...                        ...         ...
1921             SUZ12_SGTN_23       SUZ12
1922            MESDC1_SGTN_11      TLNRD1
1923            MESDC1_SGTN_12      TLNRD1
1924              TP53_SGTN_01        TP53
1925              TP53_SGTN_02        TP53

[1926 rows x 2 columns]
```

## Code cell 77

```python
# Create a mapping from the Gene Symbol (which matches your guides)
# to the Gene Index (which matches the matrix rownames)

# Check which column in mdata['GEX'].var holds the symbol
symbol_col = 'symbol' 

gene_map = mdata['GEX'].var[[symbol_col]].reset_index()
gene_map.columns = ['response_id', 'grna_target'] # Rename for easy merging in R

# Save it
gene_map.to_csv(f"{temp_dir}/gene_symbol_map.csv", index=False)
print("✅ Gene mapping saved.")
```

### Text output (stream)

```text
✅ Gene mapping saved.

```

## Code cell 78

```python
gene_map
```

### Text output (text/plain)

```text
              response_id      grna_target
0       ENSG00000227232.5           WASH7P
1       ENSG00000238009.6  ENSG00000238009
2       ENSG00000241860.7  ENSG00000241860
3       ENSG00000279457.4           WASH9P
4      ENSG00000228463.10  ENSG00000228463
...                   ...              ...
12727   ENSG00000198727.2           MT-CYB
12728   ENSG00000210196.2            MT-TP
12729   ENSG00000271254.7  ENSG00000271254
12730   ENSG00000276345.1  ENSG00000276345
12731   ENSG00000276256.1  ENSG00000276256

[12732 rows x 2 columns]
```

## Markdown cell 79

#### Import data to SCEPTRE (Rcell)

## Code cell 80

```python
%load_ext rpy2.ipython
```

### Text output (stream)

```text
The rpy2.ipython extension is already loaded. To reload it, use:
  %reload_ext rpy2.ipython

```

## Code cell 81

```python
%%R
library(sceptre)
library(Matrix)

# Define path
data_dir <- "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer"

# 1. LOAD MATRICES
cat("Loading matrices into R...\n")
response_mat <- readMM(file.path(data_dir, "response_matrix.mtx"))
grna_mat <- readMM(file.path(data_dir, "grna_matrix.mtx"))

# 2. ATTACH NAMES
gene_names <- read.csv(file.path(data_dir, "gene_names.csv"), header=FALSE)$V1
guide_names <- read.csv(file.path(data_dir, "guide_names.csv"), header=FALSE)$V1
cell_names <- read.csv(file.path(data_dir, "cell_names.csv"), header=FALSE)$V1

rownames(response_mat) <- gene_names
colnames(response_mat) <- cell_names

rownames(grna_mat) <- guide_names
colnames(grna_mat) <- cell_names
```

### Text output (stream)

```text
Loading matrices into R...

```

## Code cell 82

```python
%%R

# 3. LOAD METADATA
covariates <- read.csv(file.path(data_dir, "covariates.csv"), row.names=1)
guide_mapping <- read.csv(file.path(data_dir, "guide_mapping.csv"))

# ==============================================================================
# FIX: RENAME 'response_p_mito' SO WE KEEP IT
# Sceptre blocks the name 'response_p_mito' but doesn't calculate it.
# We rename it to 'mt_percent' so we can use it as a covariate later.
# ==============================================================================

# 1. Rename response_p_mito if it exists
if ("response_p_mito" %in% colnames(covariates)) {
    cat("Renaming 'response_p_mito' to 'mt_percent'...\n")
    names(covariates)[names(covariates) == "response_p_mito"] <- "mt_percent"
}

# 2. DROP OTHER RESERVED COLUMNS (which Sceptre DOES calculate)
# We still drop n_umis and n_nonzero because Sceptre re-calculates those automatically.
reserved_cols <- c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis")
covariates <- covariates[, !(colnames(covariates) %in% reserved_cols), drop=FALSE]

# 4. CREATE SCEPTRE OBJECT
sceptre_object <- import_data(
  response_matrix = response_mat,
  grna_matrix = grna_mat,
  grna_target_data_frame = guide_mapping,
  extra_covariates = covariates,
  moi = "high"
)

print("✅ Sceptre object created!")
print(head(sceptre_object@covariate_data_frame))
```

### Text output (stream)

```text
Renaming 'response_p_mito' to 'mt_percent'...
[1] "\u2705 Sceptre object created!"
                           response_n_nonzero response_n_umis grna_n_nonzero
TAGGATGAATTGAGGATGCTTGGG_2               1101       2078.2321              6
TAGGATGAAGCCATGCCGCCCGGA_2                494        849.4334             28
CCAGTTCACAATGGAATGCTTGGG_3                969       1890.6964             24
ATCATTCCAGTCACTACGCCCGGA_2                355        549.8333              9
GCGAGTAAGCTAACGACACAATTG_3                613       1287.2493             43
TCTTCACACAATGGAACTTTGGTC_1                694       1432.8333             59
                           grna_n_umis                  batch mt_percent
TAGGATGAATTGAGGATGCTTGGG_2          98 300genes:eTN5_GEX_sub2   13.30458
TAGGATGAAGCCATGCCGCCCGGA_2         142 300genes:eTN5_GEX_sub2   17.54111
CCAGTTCACAATGGAATGCTTGGG_3         245 300genes:eTN5_GEX_sub3   15.07381
ATCATTCCAGTCACTACGCCCGGA_2          18 300genes:eTN5_GEX_sub2   16.91421
GCGAGTAAGCTAACGACACAATTG_3         198 300genes:eTN5_GEX_sub3   22.10139
TCTTCACACAATGGAACTTTGGTC_1         417 300genes:eTN5_GEX_sub1   22.12400

```

## Code cell 83

```python
%%R
# 1. Load the translation map
gene_map <- read.csv(file.path(data_dir, "gene_symbol_map.csv"))

# 2. Start with guide mapping (which has 'grna_target' as Symbol)
discovery_pairs_highmoi <- guide_mapping

# 3. Merge to attach the correct response_id (Gene ID)
#    This matches 'grna_target' in guides to 'grna_target' in the gene map
discovery_pairs_highmoi <- merge(discovery_pairs_highmoi, gene_map, by="grna_target")

# 4. Select ONLY the two columns Sceptre needs
#    (We filter out non-targeting guides automatically via the merge)
discovery_pairs_highmoi <- discovery_pairs_highmoi[, c("grna_target", "response_id")]

# 5. DEDUPLICATE
#    Sceptre performs analysis at the *Target* level, not the *Guide* level.
#    We only need one row per Target-Gene pair.
discovery_pairs_highmoi <- unique(discovery_pairs_highmoi)

print(paste("Unique pairs to test:", nrow(discovery_pairs_highmoi)))
print(head(discovery_pairs_highmoi))
```

### Text output (stream)

```text
[1] "Unique pairs to test: 279"
   grna_target        response_id
1        ACAA1 ENSG00000060971.19
7       ACVR1B ENSG00000135503.13
13      ACVR2A ENSG00000121989.15
19      ACVR2B ENSG00000114739.14
25       ADCY2 ENSG00000078295.17
31      AMOTL2 ENSG00000114019.15

```

## Code cell 84

```python
%%R 

sceptre_object_highmoi <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs_highmoi,
  side = "both", control_group = "complement")
```

## Code cell 85

```python
%%R

png(file.path(data_dir, "grna_count_distributions.png"), width = 3000, height = 2000, res = 300)

print(plot_grna_count_distributions(
  sceptre_object = sceptre_object_highmoi,
  n_grnas_to_plot = 8
))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "grna_count_distributions.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/grna_count_distributions.png

```

## Code cell 86

```python
%%R

plot_covariates(sceptre_object_highmoi, p_mito_threshold = 0.25)
```

### Text output (text/plain)

```text
<IPython.core.display.Image object>
```

## Code cell 87

```python
%%R
sceptre_object_highmoi <- assign_grnas(
  sceptre_object = sceptre_object_highmoi,
  method = "mixture",
  parallel = FALSE
)
```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178488.23-P1P2 (1 of 1926)
Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178620.23-P1P2 (5 of 1926)
Performing gRNA-to-cell assignments for gRNA ACTL7B_-_111617995.23-P1P2 (10 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR1B_+_52345762.23-P1P2 (15 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2A_+_148602343.23-P1P2 (20 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_+_38495339.23-P1P2 (25 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_-_38495819.23-P1P2 (30 of 1926)
Performing gRNA-to-cell assignments for gRNA ADCY2_-_7396429.23-P1P2 (35 of 1926)
Performing gRNA-to-cell assignments for gRNA AMOTL2_+_134093279.23-P1P2 (40 of 1926)
Performing gRNA-to-cell assignments for gRNA ANP32E_+_150208342.23-P1P2 (45 of 1926)
Performing gRNA-to-cell assignments for gRNA APAF1_-_99039071.23-P1P2 (50 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_+_72766725.23-P1P2 (55 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_-_72766755.23-P1P2 (60 of 1926)
Performing gRNA-to-cell assignments for gRNA ARTN_-_44399429.23-P1 (65 of 1926)
Performing gRNA-to-cell assignments for gRNA ASCC3_-_101329193.23-P1P2 (70 of 1926)
Performing gRNA-to-cell assignments for gRNA ASH2L_-_37963059.23-P1P2 (75 of 1926)
Performing gRNA-to-cell assignments for gRNA ATMIN_+_81069661.23-P1P2 (80 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_+_111992112.23-P1P2 (85 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_-_111992129.23-P1P2 (90 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP6AP1_-_153657075.23-P1P2 (95 of 1926)
Performing gRNA-to-cell assignments for gRNA ATXN2L_-_28834668.23-P1P2 (100 of 1926)
Performing gRNA-to-cell assignments for gRNA BAD_+_64052102.23-P1P2 (105 of 1926)
Performing gRNA-to-cell assignments for gRNA BAG6_+_31620398.23-P1P2 (110 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_+_49458268.23-P1P2 (115 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_-_49458323.23-P1P2 (120 of 1926)
Performing gRNA-to-cell assignments for gRNA BCAS2_-_115124045.23-P1P2 (125 of 1926)
Performing gRNA-to-cell assignments for gRNA BMPR1A_-_88516386.23-P1P2 (130 of 1926)
Performing gRNA-to-cell assignments for gRNA BRD3_-_136932911.23-P1P2 (135 of 1926)
Performing gRNA-to-cell assignments for gRNA C17orf49_+_6918114.23-P2 (140 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_+_31532436.23-P1P2 (145 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_-_31532650.23-P1P2 (150 of 1926)
Performing gRNA-to-cell assignments for gRNA CACFD1_-_136325160.23-P1P2 (155 of 1926)
Performing gRNA-to-cell assignments for gRNA CACNG1_+_65040756.23-P1P2 (160 of 1926)
Performing gRNA-to-cell assignments for gRNA CAPN5_+_76778225.23-P1P2 (165 of 1926)
Performing gRNA-to-cell assignments for gRNA CASP3_+_185570560.23-P1P2 (170 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_+_107384793.23-P1P2 (175 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_-_107384635.23-P1P2 (180 of 1926)
Performing gRNA-to-cell assignments for gRNA CCM2_-_45039941.23-P1P2 (185 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNC_-_100016267.23-P1P2 (190 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNT2_-_135676416.23-P1P2 (195 of 1926)
Performing gRNA-to-cell assignments for gRNA CCT6B_+_33288381.23-P1P2 (200 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_+_30366413.23-P1P2 (205 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_-_30366630.23-P1P2 (210 of 1926)
Performing gRNA-to-cell assignments for gRNA CD99L2_-_150067120.23-P1P2 (215 of 1926)
Performing gRNA-to-cell assignments for gRNA CDH5_-_66400630.23-P1P2 (220 of 1926)
Performing gRNA-to-cell assignments for gRNA CDKN1A_-_36646538.23-P1P2 (225 of 1926)
Performing gRNA-to-cell assignments for gRNA CENPC_+_68411190.23-P1P2 (230 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_+_30773828.23-P1P2 (235 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_-_30773391.23-P1P2 (240 of 1926)
Performing gRNA-to-cell assignments for gRNA CHMP3_-_86790451.23-P1P2 (245 of 1926)
Performing gRNA-to-cell assignments for gRNA CHORDC1_-_89955894.23-P1P2 (250 of 1926)
Performing gRNA-to-cell assignments for gRNA CNEP1R1_+_50059284.23-P1P2 (255 of 1926)
Performing gRNA-to-cell assignments for gRNA CNOT2_+_70637651.23-P1P2 (260 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_+_83956339.23-P1P2 (265 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_-_83956572.23-P1P2 (270 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS6_-_99686676.23-P1P2 (275 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_-_31634088.23-P1P2 (280 of 1926)
Performing gRNA-to-cell assignments for gRNA CTR9_+_10772881.23-P1P2 (285 of 1926)
Performing gRNA-to-cell assignments for gRNA CUL3_-_225449667.23-P1P2 (290 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_+_145150267.23-P1P2 (295 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_-_145149996.23-P1P2 (300 of 1926)
Performing gRNA-to-cell assignments for gRNA CYCS_+_25164869.23-P1P2 (305 of 1926)
Performing gRNA-to-cell assignments for gRNA DACH1_-_72441101.23-P1P2 (310 of 1926)
Performing gRNA-to-cell assignments for gRNA DBR1_+_137893673.23-P1P2 (315 of 1926)
Performing gRNA-to-cell assignments for gRNA DCTN2_-_57940910.23-P1P2 (320 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_+_61851532.23-P1P2 (325 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_-_61851812.23-P1P2 (330 of 1926)
Performing gRNA-to-cell assignments for gRNA DGCR8_-_20068158.23-P1P2 (335 of 1926)
Performing gRNA-to-cell assignments for gRNA DNAJC17_-_41099215.23-P1P2 (340 of 1926)
Performing gRNA-to-cell assignments for gRNA DNM2_+_10829157.23-P1P2 (345 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM1_+_49575018.23-P1P2 (350 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_+_130700056.23-P1P2 (355 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_-_130700076.23-P1P2 (360 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM3_+_155113020.23-P1P2 (365 of 1926)
Performing gRNA-to-cell assignments for gRNA DROSHA_-_31531982.23-P1P2 (370 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_+_3985081.23-P1 (375 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_-_3985443.23-P1 (380 of 1926)
Performing gRNA-to-cell assignments for gRNA EHMT2_-_31865199.23-P1P2 (385 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3F_+_8009347.23-P1P2 (390 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3M_-_32605405.23-P1P2 (395 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4E_+_99849831.23-P1P2 (400 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184032364.23-P1P2 (405 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184033042.23-P1P2 (410 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G2_+_10830447.23-P1 (415 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF5A_-_7211217.23-P1P2 (420 of 1926)
Performing gRNA-to-cell assignments for gRNA TCEB3_-_24069975.23-P1P2 (425 of 1926)
Performing gRNA-to-cell assignments for gRNA EMC2_+_109456080.23-P1P2 (430 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823220.23-P1P2 (435 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823252.23-P1P2 (440 of 1926)
Performing gRNA-to-cell assignments for gRNA EOMES_-_27763342.23-P1P2 (445 of 1926)
Performing gRNA-to-cell assignments for gRNA ERCC2_+_45873802.23-P1P2 (450 of 1926)
Performing gRNA-to-cell assignments for gRNA ERO1L_+_53162425.23-P1P2 (455 of 1926)
Performing gRNA-to-cell assignments for gRNA ETS1_+_128392044.23-P1P2 (460 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_+_36134592.23-P1P2 (465 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_-_36134789.23-P1P2 (470 of 1926)
Performing gRNA-to-cell assignments for gRNA EWSR1_-_29664347.23-P1P2 (475 of 1926)
Performing gRNA-to-cell assignments for gRNA EXOSC8_-_37574979.23-P1P2 (480 of 1926)
Performing gRNA-to-cell assignments for gRNA FBXW11_-_171433300.23-P1P2 (485 of 1926)
Performing gRNA-to-cell assignments for gRNA FEN1_+_61560782.23-P1P2 (490 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_+_75532080.23-P1P2 (495 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_-_75532461.23-P1P2 (500 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT7_-_139926495.23-P1P2 (505 of 1926)
Performing gRNA-to-cell assignments for gRNA G3BP1_-_151151499.23-P1P2 (510 of 1926)
Performing gRNA-to-cell assignments for gRNA GCK_+_44229021.23-P1P2 (515 of 1926)
Performing gRNA-to-cell assignments for gRNA GOLT1B_+_21654777.23-P1P2 (520 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_+_45000529.23-P1P2 (525 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_-_45000948.23-P1P2 (530 of 1926)
Performing gRNA-to-cell assignments for gRNA GPI_+_34856284.23-P1P2 (535 of 1926)
Performing gRNA-to-cell assignments for gRNA GPR107_-_132816044.23-P1P2 (540 of 1926)
Performing gRNA-to-cell assignments for gRNA GSPT1_-_12009714.23-P1P2 (545 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A1_+_81687215.23-P1P2 (550 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_+_59949355.23-P1P2 (555 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_-_59949657.23-P1P2 (560 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2H4_-_30876185.23-P1P2 (565 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF3C2_+_27579604.23-P1P2 (570 of 1926)
Performing gRNA-to-cell assignments for gRNA HAND1_+_153857732.23-P1P2 (575 of 1926)
Performing gRNA-to-cell assignments for gRNA HDAC2_-_114292253.23-P1P2 (580 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_+_81399312.23-P1P2 (585 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_-_81399351.23-P1P2 (590 of 1926)
Performing gRNA-to-cell assignments for gRNA HINT2_+_35814995.23-P1P2 (595 of 1926)
Performing gRNA-to-cell assignments for gRNA HNF4A_+_43029971.23-P1P2 (600 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPM_-_8509885.23-P1 (605 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPU_+_245027753.23-P1P2 (610 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_+_43702292.23-P1P2 (615 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_-_43702305.23-P1P2 (620 of 1926)
Performing gRNA-to-cell assignments for gRNA HSPA5_-_128003534.23-P1P2 (625 of 1926)
Performing gRNA-to-cell assignments for gRNA HUWE1_-_53710540.23-P1P2 (630 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_-_30007588.23-P1P2 (635 of 1926)
Performing gRNA-to-cell assignments for gRNA INTS10_+_19675000.23-P1P2 (640 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_+_65903046.23-P1P2 (645 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_-_65903523.23-P1P2 (650 of 1926)
Performing gRNA-to-cell assignments for gRNA C7orf26_+_6630358.23-P1P2 (655 of 1926)
Performing gRNA-to-cell assignments for gRNA IPO8_+_30848771.23-P1P2 (660 of 1926)
Performing gRNA-to-cell assignments for gRNA ITGB1BP1_+_9563543.23-P1P2 (665 of 1926)
Performing gRNA-to-cell assignments for gRNA SSFA2_+_182756743.23-P1P2 (670 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_+_109914946.23-P1P2 (675 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_-_109915112.23-P1P2 (680 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM1A_-_23346409.23-P1P2 (685 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM3B_-_137688431.23-P1P2 (690 of 1926)
Performing gRNA-to-cell assignments for gRNA KDR_-_55991659.23-P1P2 (695 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA0232_+_6784651.23-P1P2 (700 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153706791.23-P1P2 (705 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153707248.23-P1P2 (710 of 1926)
Performing gRNA-to-cell assignments for gRNA LAMTOR4_-_99746921.23-P1P2 (715 of 1926)
Performing gRNA-to-cell assignments for gRNA LARP4_-_50794864.23-P1P2 (720 of 1926)
Performing gRNA-to-cell assignments for gRNA LATS2_-_21635336.23-P1P2 (725 of 1926)
Performing gRNA-to-cell assignments for gRNA LEMD3_+_65563441.23-P1P2 (730 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_+_54960428.23-P1P2 (735 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_-_54960616.23-P1P2 (740 of 1926)
Performing gRNA-to-cell assignments for gRNA LHFPL2_-_77944519.23-P1P2 (745 of 1926)
Performing gRNA-to-cell assignments for gRNA LSM12_-_42143603.23-P1P2 (750 of 1926)
Performing gRNA-to-cell assignments for gRNA MANBAL_+_35918136.23-P1P2 (755 of 1926)
Performing gRNA-to-cell assignments for gRNA MAU2_+_19431732.23-P1P2 (760 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_+_47807678.23-P1P2 (765 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_-_47807877.23-P1P2 (770 of 1926)
Performing gRNA-to-cell assignments for gRNA MEAF6_-_37980319.23-P1P2 (775 of 1926)
Performing gRNA-to-cell assignments for gRNA MED1_-_37607140.23-P1P2 (780 of 1926)
Performing gRNA-to-cell assignments for gRNA MED12_+_70338541.23-P1P2 (785 of 1926)
Performing gRNA-to-cell assignments for gRNA MED15_+_20862189.23-P1 (790 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_+_57479335.23-P1P2 (795 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_-_57479437.23-P1P2 (800 of 1926)
Performing gRNA-to-cell assignments for gRNA MEMO1_-_32264882.23-P1 (805 of 1926)
Performing gRNA-to-cell assignments for gRNA MIXL1_+_226411425.23-P1P2 (810 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL28_+_420487.23-P1P2 (815 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL32_+_42972018.23-P1P2 (820 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_+_179322181.23-P1P2 (825 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_-_179322213.23-P1P2 (830 of 1926)
Performing gRNA-to-cell assignments for gRNA LUZP6_-_135612750.23-P1P2 (835 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA20_-_19998023.23-P1P2 (840 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA25_+_112546538.23-P1P2 (845 of 1926)
Performing gRNA-to-cell assignments for gRNA NCBP1_+_100396082.23-P1P2 (850 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_+_172428630.23-P1P2 (855 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_-_172428841.23-P1P2 (860 of 1926)
Performing gRNA-to-cell assignments for gRNA NEDD8_-_24701363.23-P1P2 (865 of 1926)
Performing gRNA-to-cell assignments for gRNA NF2_+_29999635.23-P1P2 (870 of 1926)
Performing gRNA-to-cell assignments for gRNA NFIB_+_14314559.23-P1P2 (875 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_+_139418329.23-P2 (880 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139418349.23-P2 (885 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139440503.23-P1 (890 of 1926)
Performing gRNA-to-cell assignments for gRNA NUFIP2_+_27621073.23-P1P2 (895 of 1926)
Performing gRNA-to-cell assignments for gRNA NUP93_-_56764090.23-P1P2 (900 of 1926)
Performing gRNA-to-cell assignments for gRNA OGT_+_70752948.23-P1P2 (905 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_+_59383170.23-P1P2 (910 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_-_59383348.23-P1P2 (915 of 1926)
Performing gRNA-to-cell assignments for gRNA PAM16_-_4401236.23-P1P2 (920 of 1926)
Performing gRNA-to-cell assignments for gRNA PAXBP1_-_34143614.23-P1P2 (925 of 1926)
Performing gRNA-to-cell assignments for gRNA PBLD_-_70092436.23-P1 (930 of 1926)
Performing gRNA-to-cell assignments for gRNA PCBP2_+_53835614.23-P2 (935 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_+_94292752.23-P2 (940 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_-_94293027.23-P2 (945 of 1926)
Performing gRNA-to-cell assignments for gRNA PGRMC2_+_129208913.23-P1P2 (950 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_-_178866211.23-P1P2 (955 of 1926)
Performing gRNA-to-cell assignments for gRNA PITPNC1_-_65373433.23-P1P2 (960 of 1926)
Performing gRNA-to-cell assignments for gRNA PLCG1_+_39766176.23-P1P2 (965 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_+_151042273.23-P2 (970 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_-_151042422.23-P2 (975 of 1926)
Performing gRNA-to-cell assignments for gRNA PLRG1_-_155471154.23-P1P2 (980 of 1926)
Performing gRNA-to-cell assignments for gRNA PLXNA2_-_208417275.23-P1P2 (985 of 1926)
Performing gRNA-to-cell assignments for gRNA POLE4_-_75185725.23-P1P2 (990 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR2M_+_57998962.23-P1P2 (995 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940045.23-P1P2 (1000 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940414.23-P1P2 (1005 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3K_-_103144.23-P1P2 (1010 of 1926)
Performing gRNA-to-cell assignments for gRNA PPHLN1_-_42720007.23-P1P2 (1015 of 1926)
Performing gRNA-to-cell assignments for gRNA PPIL4_+_149867146.23-P1P2 (1020 of 1926)
Performing gRNA-to-cell assignments for gRNA PPP2R1A_+_52693317.23-P1 (1025 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_+_53195217.23-P1P2 (1030 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_-_53195332.23-P1P2 (1035 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKD1_-_30396996.23-P1P2 (1040 of 1926)
Performing gRNA-to-cell assignments for gRNA PRMT5_+_23398554.23-P1 (1045 of 1926)
Performing gRNA-to-cell assignments for gRNA C2orf53_-_27362230.23-P1P2 (1050 of 1926)
Performing gRNA-to-cell assignments for gRNA PRSS8_+_31147010.23-P1P2 (1055 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_+_36909062.23-P1P2 (1060 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_-_36909230.23-P1P2 (1065 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMD4_-_151227289.23-P1P2 (1070 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMG4_-_3259227.23-P1P2 (1075 of 1926)
Performing gRNA-to-cell assignments for gRNA PTAR1_+_72374835.23-P1P2 (1080 of 1926)
Performing gRNA-to-cell assignments for gRNA PTBP1_+_797463.23-P1P2 (1085 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_+_89623388.23-P1P2 (1090 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_-_89623526.23-P1P2 (1095 of 1926)
Performing gRNA-to-cell assignments for gRNA PTMS_-_6875591.23-P1P2 (1100 of 1926)
Performing gRNA-to-cell assignments for gRNA PYGO2_-_154934135.23-P1P2 (1105 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB18_-_27793292.23-P1P2 (1110 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB6A_+_73471983.23-P1P2 (1115 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_+_76251918.23-P1P2 (1120 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_-_76252008.23-P1P2 (1125 of 1926)
Performing gRNA-to-cell assignments for gRNA GNB2L1_-_180670840.23-P1P2 (1130 of 1926)
Performing gRNA-to-cell assignments for gRNA RAD21_-_117886635.23-P1P2 (1135 of 1926)
Performing gRNA-to-cell assignments for gRNA RANGRF_+_8192257.23-P1P2 (1140 of 1926)
Performing gRNA-to-cell assignments for gRNA RBBP5_+_205091050.23-P1P2 (1145 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM14-RBM4_+_66406280.23-P2 (1150 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM33_+_155437802.23-P1P2 (1155 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM4_+_66406319.23-P2 (1160 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_+_114271352.23-P1P2 (1165 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_-_114271700.23-P1P2 (1170 of 1926)
Performing gRNA-to-cell assignments for gRNA RBX1_-_41347852.23-P1P2 (1175 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1551_+_32112399.23-P1P2 (1180 of 1926)
Performing gRNA-to-cell assignments for gRNA RFX3_-_3525750.23-P1P2 (1185 of 1926)
Performing gRNA-to-cell assignments for gRNA RGP1_+_35749328.23-P1P2 (1190 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_+_95067102.23-P1P2 (1195 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_-_95067164.23-P1P2 (1200 of 1926)
Performing gRNA-to-cell assignments for gRNA RIBC1_-_53449606.23-P1P2 (1205 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1432_+_5629173.23-P1P2 (1210 of 1926)
Performing gRNA-to-cell assignments for gRNA RLF_+_40627522.23-P1P2 (1215 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF115_-_145611148.23-P1P2 (1220 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_+_104296142.23-P1P2 (1225 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_-_104296250.23-P1P2 (1230 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF40_-_30773681.23-P1P2 (1235 of 1926)
Performing gRNA-to-cell assignments for gRNA RPL37A_-_217363600.23-P1P2 (1240 of 1926)
Performing gRNA-to-cell assignments for gRNA RPS15A_-_18801590.23-P1P2 (1245 of 1926)
Performing gRNA-to-cell assignments for gRNA SAV1_+_51134992.23-P1P2 (1250 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_+_33104038.23-P1P2 (1255 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_-_33104343.23-P1P2 (1260 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF8_+_155054976.23-P1P2 (1265 of 1926)
Performing gRNA-to-cell assignments for gRNA EPT1_+_26569050.23-P1P2 (1270 of 1926)
Performing gRNA-to-cell assignments for gRNA SERF2_+_44069659.23-P2 (1275 of 1926)
Performing gRNA-to-cell assignments for gRNA SERPING1_+_57365137.23-P1P2 (1280 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_+_122242041.23-P1P2 (1285 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_-_122242294.23-P1P2 (1290 of 1926)
Performing gRNA-to-cell assignments for gRNA SETDB1_-_150898815.23-P1P2 (1295 of 1926)
Performing gRNA-to-cell assignments for gRNA SF1_-_64545937.23-P1P2 (1300 of 1926)
Performing gRNA-to-cell assignments for gRNA CCDC101_-_28565268.23-P1P2 (1305 of 1926)
Performing gRNA-to-cell assignments for gRNA SKA3_+_21750519.23-P1P2 (1310 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_+_78591497.23-P1P2 (1315 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_-_78591737.23-P1P2 (1320 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC29A1_-_44187438.23-P2 (1325 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC35C1_-_45826711.23-P1P2 (1330 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD1_+_146403428.23-P1P2 (1335 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD3_+_67358229.23-P1P2 (1340 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_+_48556623.23-P1P2 (1345 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_-_48556657.23-P1P2 (1350 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD5_-_135468694.23-P1P2 (1355 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA4_-_11071741.23-P1P2 (1360 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA5_+_144435243.23-P1P2 (1365 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCB1_+_24129228.23-P1 (1370 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_+_47822975.23-P1P2 (1375 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_-_47822893.23-P1P2 (1380 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCD3_-_150974212.23-P1 (1385 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCE1_-_38803669.23-P1P2 (1390 of 1926)
Performing gRNA-to-cell assignments for gRNA SMC2_-_106856631.23-P1P2 (1395 of 1926)
Performing gRNA-to-cell assignments for gRNA SMG5_+_156252537.23-P1P2 (1400 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_+_103852.23-P1P2 (1405 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_-_104242.23-P1P2 (1410 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP27_-_70121181.23-P1P2 (1415 of 1926)
Performing gRNA-to-cell assignments for gRNA SNUPN_-_75917960.23-P1P2 (1420 of 1926)
Performing gRNA-to-cell assignments for gRNA SOX4_+_21594213.23-P1 (1425 of 1926)
Performing gRNA-to-cell assignments for gRNA SPIN1_+_91003737.23-P1P2 (1430 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_+_2328680.23-P1P2 (1435 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_-_2328912.23-P1P2 (1440 of 1926)
Performing gRNA-to-cell assignments for gRNA SPTY2D1_-_18655846.23-P1P2 (1445 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP19_+_112197188.23-P1P2 (1450 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP54_-_35452244.23-P1P2 (1455 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP68_+_74068540.23-P1P2 (1460 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_+_126138730.23-P1P2 (1465 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_-_126138742.23-P1P2 (1470 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPRB_-_133524719.23-P1P2 (1475 of 1926)
Performing gRNA-to-cell assignments for gRNA SRRM2_-_2802714.23-P1P2 (1480 of 1926)
Performing gRNA-to-cell assignments for gRNA SRSF6_+_42086993.23-P1 (1485 of 1926)
Performing gRNA-to-cell assignments for gRNA SSB_+_170655450.23-P1P2 (1490 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_+_57102902.23-P1P2 (1495 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_-_57103355.23-P1P2 (1500 of 1926)
Performing gRNA-to-cell assignments for gRNA STAG2_-_123095658.23-P1 (1505 of 1926)
Performing gRNA-to-cell assignments for gRNA STK40_+_36851507.23-P1P2 (1510 of 1926)
Performing gRNA-to-cell assignments for gRNA SUPT20H_+_37633798.23-P1P2 (1515 of 1926)
Performing gRNA-to-cell assignments for gRNA SUZ12_+_30264157.23-P1P2 (1520 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_+_64901973.23-P1P2 (1525 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_-_64901981.23-P1P2 (1530 of 1926)
Performing gRNA-to-cell assignments for gRNA TADA2B_-_7045121.23-P1P2 (1535 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF10_+_6633406.23-P1P2 (1540 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF2_+_120845032.23-P1P2 (1545 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF3_-_7860540.23-P1P2 (1550 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_+_60640814.23-P1P2 (1555 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_-_60640498.23-P1P2 (1560 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF6L_-_62538954.23-P1P2 (1565 of 1926)
Performing gRNA-to-cell assignments for gRNA TAOK1_-_27717907.23-P1P2 (1570 of 1926)
Performing gRNA-to-cell assignments for gRNA FAM208A_+_56717273.23-P1P2 (1575 of 1926)
Performing gRNA-to-cell assignments for gRNA TBC1D9_+_141677481.23-P1P2 (1580 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_+_170863524.23-P1P2 (1585 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_-_170863481.23-P1P2 (1590 of 1926)
Performing gRNA-to-cell assignments for gRNA TBX3_-_115121998.23-P1P2 (1595 of 1926)
Performing gRNA-to-cell assignments for gRNA T_-_166581333.23-P1P2 (1600 of 1926)
Performing gRNA-to-cell assignments for gRNA TGFBR2_+_30648526.23-P1P2 (1605 of 1926)
Performing gRNA-to-cell assignments for gRNA THBD_+_23029196.23-P2 (1610 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_+_2977592.23-P1P2 (1615 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_-_2977643.23-P1P2 (1620 of 1926)
Performing gRNA-to-cell assignments for gRNA MESDC1_-_81293340.23-P1P2 (1625 of 1926)
Performing gRNA-to-cell assignments for gRNA TMEM41B_+_9336156.23-P1P2 (1630 of 1926)
Performing gRNA-to-cell assignments for gRNA TMTC3_-_88536127.23-P1P2 (1635 of 1926)
Performing gRNA-to-cell assignments for gRNA TOP1_+_39657538.23-P1P2 (1640 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_+_7590531.23-P1P2 (1645 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_-_7590704.23-P1P2 (1650 of 1926)
Performing gRNA-to-cell assignments for gRNA TRIM24_-_138145308.23-P1P2 (1655 of 1926)
Performing gRNA-to-cell assignments for gRNA TRRAP_-_98476465.23-P1P2 (1660 of 1926)
Performing gRNA-to-cell assignments for gRNA TSTD2_+_100395594.23-P1P2 (1665 of 1926)
Performing gRNA-to-cell assignments for gRNA UBE2E1_+_23847541.23-P1P2 (1670 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_+_38924011.23-P1P2 (1675 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_-_38924421.23-P1P2 (1680 of 1926)
Performing gRNA-to-cell assignments for gRNA UPF2_+_12084799.23-P1P2 (1685 of 1926)
Performing gRNA-to-cell assignments for gRNA UXT_-_47518505.23-P1P2 (1690 of 1926)
Performing gRNA-to-cell assignments for gRNA VHL_-_10183525.23-P1P2 (1695 of 1926)
Performing gRNA-to-cell assignments for gRNA VPS54_+_64246104.23-P1P2 (1700 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_+_28821698.23-P1P2 (1705 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_-_28822074.23-P1P2 (1710 of 1926)
Performing gRNA-to-cell assignments for gRNA WAPAL_+_88281533.23-P1P2 (1715 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR12_-_203776383.23-P1P2 (1720 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR77_+_111991860.23-P1P2 (1725 of 1926)
Performing gRNA-to-cell assignments for gRNA WIZ_-_15543016.23-P1P2 (1730 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_+_9983328.23-P1P2 (1735 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_-_9983510.23-P1P2 (1740 of 1926)
Performing gRNA-to-cell assignments for gRNA XPO5_-_43543670.23-P1P2 (1745 of 1926)
Performing gRNA-to-cell assignments for gRNA XRN1_-_142166648.23-P1P2 (1750 of 1926)
Performing gRNA-to-cell assignments for gRNA YARS2_+_32908787.23-P1P2 (1755 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB22_+_33285454.23-P1P2 (1760 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_+_6640093.23-P1P2 (1765 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_-_6640138.23-P1P2 (1770 of 1926)
Performing gRNA-to-cell assignments for gRNA ZC3H8_-_113012435.23-P1P2 (1775 of 1926)
Performing gRNA-to-cell assignments for gRNA ZNF281_+_200379151.23-P1P2 (1780 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00117 (1785 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00233 (1790 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00518 (1795 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00656 (1800 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01226 (1805 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01730 (1810 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02009 (1815 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02288 (1820 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02777 (1825 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03026 (1830 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03252 (1835 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03620 (1840 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_SGTN_04 (1845 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_SGTN_25 (1850 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_SGTN_26 (1855 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1044 (1860 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1414 (1865 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1554 (1870 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1785 (1875 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2009 (1880 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2271 (1885 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2425 (1890 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2614 (1895 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_293 (1900 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_541 (1905 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_725 (1910 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_860 (1915 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_SGTN_14 (1920 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_SGTN_01 (1925 of 1926)

```

## Code cell 88

```python
%%R
library(Matrix)

# 1. Get the assignment matrix
A <- get_grna_assignments(sceptre_object_highmoi)

# 2. Ensure it is a numeric sparse matrix (1/0)
if (inherits(A, "Matrix") || is.matrix(A)) {
    A_num <- if (inherits(A, "Matrix")) drop0(A * 1) else (A * 1)
} else {
    stop("Output is not a matrix.")
}

# --- FIX: ROBUST NAME RECOVERY ---
# We retrieve names using Sceptre's standard metadata accessors
# instead of touching internal slots like @grna_matrix.

# Try to get cell names from the matrix first
current_cell_names <- colnames(A_num)
current_guide_names <- rownames(A_num)

# If missing, grab them from the main object's covariate data frame
if (is.null(current_cell_names)) {
    message("⚠️ Matrix missing colnames. Fetching from sceptre_object...")
    # The covariate data frame ALWAYS has the cell IDs as rownames
    current_cell_names <- rownames(sceptre_object_highmoi@covariate_data_frame)
    
    # Safety check: lengths must match
    if (length(current_cell_names) != ncol(A_num)) {
        stop(paste("Mismatch! Matrix has", ncol(A_num), "cols but object has", length(current_cell_names), "cells."))
    }
    colnames(A_num) <- current_cell_names
}

# If guide names are missing, they should match the unique values in your guide mapping
if (is.null(current_guide_names)) {
    message("⚠️ Matrix missing rownames. Fetching from guide mapping...")
    # Get unique targets from the mapping table you loaded earlier
    # Note: This assumes the matrix rows are sorted/aligned to these unique targets.
    # In 'mixture' method, rows usually correspond to the 'grna_target' labels.
    
    # SAFE FALLBACK: If we can't be 100% sure of order, we warn the user.
    # But usually, they match the 'grna_target' column from your input.
    unique_targets <- sort(unique(guide_mapping$grna_target))
    
    if (length(unique_targets) == nrow(A_num)) {
         rownames(A_num) <- unique_targets
    } else {
         warning("Could not automatically align target names. Using generic IDs.")
         rownames(A_num) <- paste0("target_", seq_len(nrow(A_num)))
    }
}
# ---------------------------------

# 3. Convert to Long Format
s <- summary(A_num)
assign_long <- data.frame(
    cell_id = colnames(A_num)[s$j],
    grna_target = rownames(A_num)[s$i],
    value = s$x
)

# 4. Save
out_path <- file.path(data_dir, "grna_assignments_long.csv")
write.csv(assign_long, out_path, row.names = FALSE)

message("✅ Saved long CSV to: ", out_path)
print(head(assign_long))
```

### Text output (stream)

```text
                     cell_id                grna_target value
1 TAGGATGAATTGAGGATGCTTGGG_2    PTEN_+_89623388.23-P1P2     1
2 TAGGATGAAGCCATGCCGCCCGGA_2  CCNT2_-_135676427.23-P1P2     1
3 TAGGATGAAGCCATGCCGCCCGGA_2      PBLD_-_70092630.23-P1     1
4 TAGGATGAAGCCATGCCGCCCGGA_2 PIK3CA_-_178866211.23-P1P2     1
5 TAGGATGAAGCCATGCCGCCCGGA_2 SMARCC1_+_47822975.23-P1P2     1
6 CCAGTTCACAATGGAATGCTTGGG_3   COPS4_-_83956484.23-P1P2     1

```

### Text output (text/plain)

```text
⚠️ Matrix missing colnames. Fetching from sceptre_object...
✅ Saved long CSV to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/grna_assignments_long.csv

```

## Code cell 89

```python
%%R
png(file.path(data_dir, "MOI.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "MOI.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/MOI.png

```

## Code cell 90

```python
%%R
sceptre_object_highmoi <- run_qc(
  sceptre_object = sceptre_object_highmoi,
  p_mito_threshold = 0.25
)

png(file.path(data_dir, "QC.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "QC.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/QC.png

```

## Code cell 91

```python
%%R
# Calibration check


sceptre_object_highmoi <- run_calibration_check(
  sceptre_object = sceptre_object_highmoi,
  parallel = FALSE
)

png(file.path(data_dir, "calibration_check.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "calibration_check.png"))

```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Constructing negative control pairs. ✓
Running precomputation on response ENSG00000143545.10 (1 of 226)
Running precomputation on response ENSG00000132603.15 (5 of 226)
Running precomputation on response ENSG00000155252.14 (10 of 226)
Running precomputation on response ENSG00000154175.18 (15 of 226)
Running precomputation on response ENSG00000230590.13 (20 of 226)
Running precomputation on response ENSG00000290146.1 (25 of 226)
Running precomputation on response ENSG00000035115.22 (30 of 226)
Running precomputation on response ENSG00000115207.15 (35 of 226)
Running precomputation on response ENSG00000140259.7 (40 of 226)
Running precomputation on response ENSG00000169180.12 (45 of 226)
Running precomputation on response ENSG00000110330.10 (50 of 226)
Running precomputation on response ENSG00000226361.5 (55 of 226)
Running precomputation on response ENSG00000100150.20 (60 of 226)
Running precomputation on response ENSG00000129696.13 (65 of 226)
Running precomputation on response ENSG00000133983.16 (70 of 226)
Running precomputation on response ENSG00000170949.18 (75 of 226)
Running precomputation on response ENSG00000167930.17 (80 of 226)
Running precomputation on response ENSG00000145901.16 (85 of 226)
Running precomputation on response ENSG00000163214.21 (90 of 226)
Running precomputation on response ENSG00000249803.6 (95 of 226)
Running precomputation on response ENSG00000237350.1 (100 of 226)
Running precomputation on response ENSG00000185627.19 (105 of 226)
Running precomputation on response ENSG00000290263.1 (110 of 226)
Running precomputation on response ENSG00000100441.10 (115 of 226)
Running precomputation on response ENSG00000147799.12 (120 of 226)
Running precomputation on response ENSG00000182378.15 (125 of 226)
Running precomputation on response ENSG00000108424.11 (130 of 226)
Running precomputation on response ENSG00000196850.6 (135 of 226)
Running precomputation on response ENSG00000127948.17 (140 of 226)
Running precomputation on response ENSG00000159899.16 (145 of 226)
Running precomputation on response ENSG00000172530.22 (150 of 226)
Running precomputation on response ENSG00000065029.15 (155 of 226)
Running precomputation on response ENSG00000158623.14 (160 of 226)
Running precomputation on response ENSG00000177731.16 (165 of 226)
Running precomputation on response ENSG00000127947.16 (170 of 226)
Running precomputation on response ENSG00000166037.11 (175 of 226)
Running precomputation on response ENSG00000250802.7 (180 of 226)
Running precomputation on response ENSG00000038358.15 (185 of 226)
Running precomputation on response ENSG00000214425.11 (190 of 226)
Running precomputation on response ENSG00000150275.20 (195 of 226)
Running precomputation on response ENSG00000030582.18 (200 of 226)
Running precomputation on response ENSG00000109118.14 (205 of 226)
Running precomputation on response ENSG00000187778.14 (210 of 226)
Running precomputation on response ENSG00000182919.15 (215 of 226)
Running precomputation on response ENSG00000254726.3 (220 of 226)
Running precomputation on response ENSG00000285816.1 (225 of 226)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00259&non-targeting_02159&non-targeting_02571&non-targeting_03029&non-targeting_03081 (1 of 90)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00355&non-targeting_02025&non-targeting_02089&non-targeting_03081&non-targeting_03384 (5 of 90)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_01730&non-targeting_02193&non-targeting_03180&non-targeting_03620&non-targeting_03664 (10 of 90)
Analyzing pairs containing gRNA group non-targeting_00043&non-targeting_00117&non-targeting_01509&non-targeting_02543&non-targeting_03026&non-targeting_03252 (15 of 90)
Analyzing pairs containing gRNA group non-targeting_00259&non-targeting_01732&non-targeting_01934&non-targeting_02193&non-targeting_02837&non-targeting_03081 (20 of 90)
Analyzing pairs containing gRNA group non-targeting_00259&non-targeting_01504&non-targeting_02009&non-targeting_02574&non-targeting_03081&non-targeting_03252 (25 of 90)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_02543&non-targeting_02777&non-targeting_02814&non-targeting_03252&non-targeting_03483 (30 of 90)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00918&non-targeting_01730&non-targeting_02025&non-targeting_02911&non-targeting_03384 (35 of 90)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00390&non-targeting_00518&non-targeting_00625&non-targeting_03081&non-targeting_03252 (40 of 90)
Analyzing pairs containing gRNA group non-targeting_00231&non-targeting_00614&non-targeting_00629&non-targeting_00918&non-targeting_02193&non-targeting_02543 (45 of 90)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00390&non-targeting_01732&non-targeting_02193&non-targeting_03081&non-targeting_03384 (50 of 90)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00518&non-targeting_00614&non-targeting_01226&non-targeting_02025&non-targeting_02837 (55 of 90)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_02193&non-targeting_02571&non-targeting_02911&non-targeting_03081&non-targeting_03180 (60 of 90)
Analyzing pairs containing gRNA group non-targeting_00324&non-targeting_00629&non-targeting_00972&non-targeting_02677&non-targeting_03252&non-targeting_03664 (65 of 90)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00656&non-targeting_02543&non-targeting_02574&non-targeting_03081&non-targeting_03384 (70 of 90)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_00259&non-targeting_00629&non-targeting_01601&non-targeting_01732&non-targeting_02837 (75 of 90)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_00355&non-targeting_02837&non-targeting_02911&non-targeting_03342&non-targeting_03620 (80 of 90)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00614&non-targeting_02009&non-targeting_02574&non-targeting_03026&non-targeting_03180 (85 of 90)
Analyzing pairs containing gRNA group non-targeting_00625&non-targeting_01226&non-targeting_01858&non-targeting_02009&non-targeting_02962&non-targeting_03158 (90 of 90)

```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/calibration_check.png
In addition: There were 50 or more warnings (use warnings() to see the first 50)

```

## Code cell 92

```python
%%R
calibration_result <- get_result(
  sceptre_object = sceptre_object_highmoi,
  analysis = "run_calibration_check"
)
head(calibration_result)
```

### Text output (stream)

```text
          response_id
               <char>
1:  ENSG00000290263.1
2: ENSG00000134371.14
3:  ENSG00000100591.8
4: ENSG00000143742.14
5: ENSG00000163930.10
6:  ENSG00000237350.1
                                                                                                               grna_target
                                                                                                                    <char>
1: non-targeting_00233&non-targeting_00390&non-targeting_01601&non-targeting_02009&non-targeting_02288&non-targeting_03029
2:
```

### Text output (stream)

```text
 non-targeting_00014&non-targeting_00259&non-targeting_00629&non-targeting_01601&non-targeting_01732&non-targeting_02837
3: non-targeting_00614&non-targeting_00972&non-targeting_02911&non-targeting_03081&non-targeting_03180&non-targeting_03384
4: non-targeting_00324&non-targeting_00629&non-targeting_00972&non-targeting_02677&non-targeting_03252&non-targeting_03664
5: non-targeting_00629&non-targeting_01185&non-targeting_01732&non-targeting_02089&non-targeting_02677&non-targeting_03483
6: non-targeting_00014&non-targeting_00117&non-targeting_00196&non-targeting_00555&non-targeting_02025&non-targeting_02193
   n_nonzero_trt n_nonzero_cntrl pass_qc     p_value log_2_fold_change
           <int>           <int>  <lgcl>       <num>             <num>
1:            25             843    TRUE 0.003853224         0.8886069
2:            44            3161    TRUE 0.016145901        -0.4608229
3:            13             997    TRUE 0.018890360        -0.8134084
4:            53            3475    TRUE 0.024000000        -0.2138143
5:            24             776    TRUE 0.024000000         0.6536926
6:            33             964    TRUE 0.032000000         0.6869690
   significant
        <lgcl>
1:       FALSE
2:       FALSE
3:       FALSE
4:       FALSE
5:       FALSE
6:       FALSE

```

## Code cell 93

```python
%%R
sceptre_object_highmoi <- run_discovery_analysis(
  sceptre_object = sceptre_object_highmoi,
  parallel = FALSE
)

png(file.path(data_dir, "discovery.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "discovery.png"))
```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Running precomputation on response ENSG00000004487.18 (1 of 227)
Running precomputation on response ENSG00000116459.11 (5 of 227)
Running precomputation on response ENSG00000163348.4 (10 of 227)
Running precomputation on response ENSG00000138018.19 (15 of 227)
Running precomputation on response ENSG00000124380.11 (20 of 227)
Running precomputation on response ENSG00000138385.16 (25 of 227)
Running precomputation on response ENSG00000170142.12 (30 of 227)
Running precomputation on response ENSG00000114019.15 (35 of 227)
Running precomputation on response ENSG00000128052.10 (40 of 227)
Running precomputation on response ENSG00000109436.8 (45 of 227)
Running precomputation on response ENSG00000113360.17 (50 of 227)
Running precomputation on response ENSG00000120733.14 (55 of 227)
Running precomputation on response ENSG00000124766.7 (60 of 227)
Running precomputation on response ENSG00000112237.13 (65 of 227)
Running precomputation on response ENSG00000213079.10 (70 of 227)
Running precomputation on response ENSG00000136280.17 (75 of 227)
Running precomputation on response ENSG00000188186.11 (80 of 227)
Running precomputation on response ENSG00000104613.12 (85 of 227)
Running precomputation on response ENSG00000179091.5 (90 of 227)
Running precomputation on response ENSG00000106723.17 (95 of 227)
Running precomputation on response ENSG00000136908.18 (100 of 227)
Running precomputation on response ENSG00000151461.20 (105 of 227)
Running precomputation on response ENSG00000171862.14 (110 of 227)
Running precomputation on response ENSG00000179119.15 (115 of 227)
Running precomputation on response ENSG00000156603.20 (120 of 227)
Running precomputation on response ENSG00000248643.5 (125 of 227)
Running precomputation on response ENSG00000182934.12 (130 of 227)
Running precomputation on response ENSG00000134283.18 (135 of 227)
Running precomputation on response ENSG00000174106.3 (140 of 227)
Running precomputation on response ENSG00000111300.10 (145 of 227)
Running precomputation on response ENSG00000120699.14 (150 of 227)
Running precomputation on response ENSG00000129559.13 (155 of 227)
Running precomputation on response ENSG00000165417.12 (160 of 227)
Running precomputation on response ENSG00000166949.17 (165 of 227)
Running precomputation on response ENSG00000167978.18 (170 of 227)
Running precomputation on response ENSG00000103549.22 (175 of 227)
Running precomputation on response ENSG00000141510.18 (180 of 227)
Running precomputation on response ENSG00000125686.12 (185 of 227)
Running precomputation on response ENSG00000154217.16 (190 of 227)
Running precomputation on response ENSG00000005206.17 (195 of 227)
Running precomputation on response ENSG00000011451.21 (200 of 227)
Running precomputation on response ENSG00000105568.19 (205 of 227)
Running precomputation on response ENSG00000124193.16 (210 of 227)
Running precomputation on response ENSG00000128191.17 (215 of 227)
Running precomputation on response ENSG00000100387.9 (220 of 227)
Running precomputation on response ENSG00000147162.15 (225 of 227)
Analyzing pairs containing gRNA group KDM1A (1 of 227)
Analyzing pairs containing gRNA group ATP5PB (5 of 227)
Analyzing pairs containing gRNA group PYGO2 (10 of 227)
Analyzing pairs containing gRNA group SELENOI (15 of 227)
Analyzing pairs containing gRNA group SNRNP27 (20 of 227)
Analyzing pairs containing gRNA group SSB (25 of 227)
Analyzing pairs containing gRNA group UBE2E1 (30 of 227)
Analyzing pairs containing gRNA group AMOTL2 (35 of 227)
Analyzing pairs containing gRNA group KDR (40 of 227)
Analyzing pairs containing gRNA group TBC1D9 (45 of 227)
Analyzing pairs containing gRNA group DROSHA (50 of 227)
Analyzing pairs containing gRNA group KDM3B (55 of 227)
Analyzing pairs containing gRNA group SOX4 (60 of 227)
Analyzing pairs containing gRNA group CCNC (65 of 227)
Analyzing pairs containing gRNA group SCAF8 (70 of 227)
Analyzing pairs containing gRNA group CCM2 (75 of 227)
Analyzing pairs containing gRNA group LAMTOR4 (80 of 227)
Analyzing pairs containing gRNA group INTS10 (85 of 227)
Analyzing pairs containing gRNA group CYC1 (90 of 227)
Analyzing pairs containing gRNA group SPIN1 (95 of 227)
Analyzing pairs containing gRNA group DPM2 (100 of 227)
Analyzing pairs containing gRNA group UPF2 (105 of 227)
Analyzing pairs containing gRNA group PTEN (110 of 227)
Analyzing pairs containing gRNA group SPTY2D1 (115 of 227)
Analyzing pairs containing gRNA group MED19 (120 of 227)
Analyzing pairs containing gRNA group RBM14-RBM4 (125 of 227)
Analyzing pairs containing gRNA group SRPRA (130 of 227)
Analyzing pairs containing gRNA group PPHLN1 (135 of 227)
Analyzing pairs containing gRNA group LEMD3 (140 of 227)
Analyzing pairs containing gRNA group NAA25 (145 of 227)
Analyzing pairs containing gRNA group EXOSC8 (150 of 227)
Analyzing pairs containing gRNA group NEDD8 (155 of 227)
Analyzing pairs containing gRNA group GTF2A1 (160 of 227)
Analyzing pairs containing gRNA group SMAD3 (165 of 227)
Analyzing pairs containing gRNA group SRRM2 (170 of 227)
Analyzing pairs containing gRNA group RNF40 (175 of 227)
Analyzing pairs containing gRNA group TP53 (180 of 227)
Analyzing pairs containing gRNA group MED1 (185 of 227)
Analyzing pairs containing gRNA group PITPNC1 (190 of 227)
Analyzing pairs containing gRNA group SPPL2B (195 of 227)
Analyzing pairs containing gRNA group WIZ (200 of 227)
Analyzing pairs containing gRNA group PPP2R1A (205 of 227)
Analyzing pairs containing gRNA group SRSF6 (210 of 227)
Analyzing pairs containing gRNA group DGCR8 (215 of 227)
Analyzing pairs containing gRNA group RBX1 (220 of 227)
Analyzing pairs containing gRNA group OGT (225 of 227)

```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/discovery.png
In addition: There were 50 or more warnings (use warnings() to see the first 50)

```

## Code cell 94

```python
%%R
discovery_result <- get_result(
  sceptre_object = sceptre_object_highmoi,
  analysis = "run_discovery_analysis"
)

write.csv(discovery_result, file.path(data_dir, "TargetKD.csv"), row.names = FALSE)
message("✅ Saved discovery to: ", file.path(data_dir, "TargetKD.csv"))

```

### Text output (text/plain)

```text
✅ Saved discovery to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/TargetKD.csv

```

## Markdown cell 95

### Day 3

## Code cell 96

```python
import scipy.io as sio
import pandas as pd
import os
```

## Code cell 97

```python
d3_mdata = mu.read(f"{data_dir}/d3_normalized_gex_and_guide.h5mu")
mdata = d3_mdata

# Filter out cells that have no guides
#mdata = mdata[mdata['GEX'].obs['guides_per_cell'] > 0]
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sceptre_env/lib/python3.10/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

## Code cell 98

```python
# Define a temporary folder for transfer
temp_dir = working_dir + "/d3_sceptre_transfer"
os.makedirs(temp_dir, exist_ok=True)

# 1. EXTRACT MATRICES (Transpose them!)
# Scanpy is (Cells × Genes), Sceptre needs (Genes × Cells)
print("Writing sparse matrices...")
sio.mmwrite(f"{temp_dir}/response_matrix.mtx", mdata['GEX'].layers['raw_counts'].T)
sio.mmwrite(f"{temp_dir}/grna_matrix.mtx", mdata['guide'].X.T)

# 2. EXTRACT NAMES (Row/Col names are lost in .mtx)
# Save Gene names and Guide names to assign them later in R
pd.Series(mdata['GEX'].var_names).to_csv(f"{temp_dir}/gene_names.csv", index=False, header=False)
pd.Series(mdata['guide'].var_names).to_csv(f"{temp_dir}/guide_names.csv", index=False, header=False)
pd.Series(mdata['GEX'].obs_names).to_csv(f"{temp_dir}/cell_names.csv", index=False, header=False)
```

### Text output (stream)

```text
Writing sparse matrices...

```

## Code cell 99

```python
mdata['GEX']
```

### Text output (text/plain)

```text
AnnData object with n_obs × n_vars = 44102 × 14712
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0'
    var: 'gene_id', 'symbol', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'n_counts', 'highly_deviant', 'highly_deviant_2k', 'binomial_deviance', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'orig_highly_variable', 'highly_variable_2k'
    uns: 'hvg', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'neighbors', 'pca', 'umap'
    obsm: 'X_pca', 'X_umap', 'guide_assignment', 'guide_counts'
    varm: 'PCs'
    layers: 'ambiguous', 'analytic_pearson_residuals', 'log1p_norm', 'mature', 'nascent', 'raw_counts', 'scran_normalization'
    obsp: 'connectivities', 'distances'
```

## Code cell 100

```python
# 3. EXTRACT COVARIATES 
rna = mdata['GEX']
covariates_df = pd.DataFrame({
    "batch": rna.obs['sample_id'],
    "grna_n_nonzero": rna.obs['guides_per_cell'],
    "grna_n_umis": rna.obs['guide_umi_counts'],
    "response_n_nonzero": rna.obs['n_genes_by_counts'],
    "response_n_umis": rna.obs['total_counts'],
    "response_p_mito": rna.obs['pct_counts_mt'],
})
covariates_df.to_csv(f"{temp_dir}/covariates.csv")
covariates_df
```

### Text output (text/plain)

```text
                                             batch  grna_n_nonzero  \
barcode                                                              
GACAGTGCCCTCTATCTCACTTTA_2  300genes:eTN5_GEX_sub2               4   
ACAGCAGACTGAGCCACTATTTCA_1  300genes:eTN5_GEX_sub1               1   
CCTCCTGAATCCTGTACGCTTAAA_3  300genes:eTN5_GEX_sub3               1   
GTACGCAAATCATTCCCATTTGGA_3  300genes:eTN5_GEX_sub3               2   
CAACCACAAAACATCGGACCTTTC_1  300genes:eTN5_GEX_sub1               3   
...                                            ...             ...   
AATCCGTCTCTTCACAGTTACGTA_1  300genes:eTN5_GEX_sub1               1   
ATCATTCCCGACACACTCACTTTA_2  300genes:eTN5_GEX_sub2               1   
AACAACCAACACAGAAGACCTTTC_1  300genes:eTN5_GEX_sub1               9   
AATGTTGCATCCTGTAGGTGGAGC_1  300genes:eTN5_GEX_sub1               3   
ACCACTGTGGTGCGAACGTCTAGG_1  300genes:eTN5_GEX_sub1               4   

                            grna_n_umis  response_n_nonzero  response_n_umis  \
barcode                                                                        
GACAGTGCCCTCTATCTCACTTTA_2         70.0                 786      1109.333252   
ACAGCAGACTGAGCCACTATTTCA_1        908.0                1086      1671.621948   
CCTCCTGAATCCTGTACGCTTAAA_3        121.0                1493      2472.592041   
GTACGCAAATCATTCCCATTTGGA_3         59.0                1792      3231.275391   
CAACCACAAAACATCGGACCTTTC_1        380.0                 847      1556.620483   
...                                 ...                 ...              ...   
AATCCGTCTCTTCACAGTTACGTA_1         16.0                 870      1732.666626   
ATCATTCCCGACACACTCACTTTA_2         33.0                 607       912.434875   
AACAACCAACACAGAAGACCTTTC_1        249.0                2098      4470.594727   
AATGTTGCATCCTGTAGGTGGAGC_1        223.0                1213      2456.680908   
ACCACTGTGGTGCGAACGTCTAGG_1         75.0                 612      1291.583984   

                            response_p_mito  
barcode                                      
GACAGTGCCCTCTATCTCACTTTA_2         9.375001  
ACAGCAGACTGAGCCACTATTTCA_1         9.092965  
CCTCCTGAATCCTGTACGCTTAAA_3        15.483077  
GTACGCAAATCATTCCCATTTGGA_3        12.533751  
CAACCACAAAACATCGGACCTTTC_1        18.148289  
...                                     ...  
AATCCGTCTCTTCACAGTTACGTA_1         9.320893  
ATCATTCCCGACACACTCACTTTA_2        18.521870  
AACAACCAACACAGAAGACCTTTC_1        14.371690  
AATGTTGCATCCTGTAGGTGGAGC_1        11.784192  
ACCACTGTGGTGCGAACGTCTAGG_1        14.439633  

[44102 rows x 6 columns]
```

## Code cell 101

```python
guide_meta_df = pd.read_csv(f"{working_dir}/300genes_guide_metadata_v43.tsv", sep='\t')
```

## Code cell 102

```python
# 4. PREPARE GUIDE MAPPING 
# Sceptre needs a table mapping 'grna_id' to 'grna_target'.

# 1. Create a "skeleton" DataFrame from the actual matrix order
# This ensures the order matches mdata['guide'].X perfectly
# We name the index 'grna_id' to match Sceptre's expectation immediately
guide_mapping = pd.DataFrame(index=mdata['guide'].var_names)
guide_mapping.index.name = 'grna_id'

# 2. Prepare metadata for merging
# Ensure 'guide_id' is the index so we can map efficiently
lookup_df = guide_meta_df.set_index('guide_id')

# 3. Map the targets onto the skeleton
# We use .loc[] to pull the specific target for each guide in the exact order required
guide_mapping['grna_target'] = lookup_df.loc[guide_mapping.index, 'intended_target_name']

# 4. Reset index so 'grna_id' becomes a column (Sceptre requires this)
guide_mapping = guide_mapping.reset_index()

# ==========================================
# VALIDATION CHECKS 
# ==========================================

# Check 1: Are there any NAs? (Guides in mdata that were missing from meta_df)
na_count = guide_mapping['grna_target'].isna().sum()
if na_count > 0:
    print(f"⚠️ WARNING: {na_count} guides have NaN targets! These will crash Sceptre.")
else:
    print("✅ All guides mapped successfully.")

# Check 2: Verify the first few rows manually
print("\nFirst 5 rows of mapping (Must match mdata['guide'].var_names order):")
print(guide_mapping.head())

# 5. Save for Sceptre
guide_mapping.to_csv(f"{temp_dir}/guide_mapping.csv", index=False)
```

### Text output (stream)

```text
✅ All guides mapped successfully.

First 5 rows of mapping (Must match mdata['guide'].var_names order):
                    grna_id grna_target
0  ACAA1_+_38178488.23-P1P2       ACAA1
1  ACAA1_+_38178517.23-P1P2       ACAA1
2  ACAA1_+_38178559.23-P1P2       ACAA1
3  ACAA1_+_38178570.23-P1P2       ACAA1
4  ACAA1_+_38178620.23-P1P2       ACAA1

```

## Code cell 103

```python
guide_mapping
```

### Text output (text/plain)

```text
                       grna_id grna_target
0     ACAA1_+_38178488.23-P1P2       ACAA1
1     ACAA1_+_38178517.23-P1P2       ACAA1
2     ACAA1_+_38178559.23-P1P2       ACAA1
3     ACAA1_+_38178570.23-P1P2       ACAA1
4     ACAA1_+_38178620.23-P1P2       ACAA1
...                        ...         ...
1921             SUZ12_SGTN_23       SUZ12
1922            MESDC1_SGTN_11      TLNRD1
1923            MESDC1_SGTN_12      TLNRD1
1924              TP53_SGTN_01        TP53
1925              TP53_SGTN_02        TP53

[1926 rows x 2 columns]
```

## Code cell 104

```python
# Create a mapping from the Gene Symbol (which matches your guides)
# to the Gene Index (which matches the matrix rownames)

# Check which column in mdata['GEX'].var holds the symbol
symbol_col = 'symbol' 

gene_map = mdata['GEX'].var[[symbol_col]].reset_index()
gene_map.columns = ['response_id', 'grna_target'] # Rename for easy merging in R

# Save it
gene_map.to_csv(f"{temp_dir}/gene_symbol_map.csv", index=False)
print("✅ Gene mapping saved.")
```

### Text output (stream)

```text
✅ Gene mapping saved.

```

## Code cell 105

```python
gene_map
```

### Text output (text/plain)

```text
              response_id      grna_target
0       ENSG00000227232.5           WASH7P
1       ENSG00000238009.6  ENSG00000238009
2       ENSG00000241860.7  ENSG00000241860
3       ENSG00000279457.4           WASH9P
4      ENSG00000228463.10  ENSG00000228463
...                   ...              ...
14707   ENSG00000275757.1        5_8S_rRNA
14708   ENSG00000276256.1  ENSG00000276256
14709   ENSG00000273748.1  ENSG00000273748
14710   ENSG00000274917.1        RNA5-8SN5
14711   ENSG00000273730.1        5_8S_rRNA

[14712 rows x 2 columns]
```

## Markdown cell 106

#### Import data to SCEPTRE (Rcell)

## Code cell 107

```python
%load_ext rpy2.ipython
```

## Code cell 108

```python
%%R
library(sceptre)
library(Matrix)

# Define path
data_dir <- "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d3_sceptre_transfer"

# 1. LOAD MATRICES
cat("Loading matrices into R...\n")
response_mat <- readMM(file.path(data_dir, "response_matrix.mtx"))
grna_mat <- readMM(file.path(data_dir, "grna_matrix.mtx"))

# 2. ATTACH NAMES
gene_names <- read.csv(file.path(data_dir, "gene_names.csv"), header=FALSE)$V1
guide_names <- read.csv(file.path(data_dir, "guide_names.csv"), header=FALSE)$V1
cell_names <- read.csv(file.path(data_dir, "cell_names.csv"), header=FALSE)$V1

rownames(response_mat) <- gene_names
colnames(response_mat) <- cell_names

rownames(grna_mat) <- guide_names
colnames(grna_mat) <- cell_names
```

### Text output (stream)

```text
Loading matrices into R...

```

## Code cell 109

```python
%%R

# 3. LOAD METADATA
covariates <- read.csv(file.path(data_dir, "covariates.csv"), row.names=1)
guide_mapping <- read.csv(file.path(data_dir, "guide_mapping.csv"))

# ==============================================================================
# FIX: RENAME 'response_p_mito' SO WE KEEP IT
# Sceptre blocks the name 'response_p_mito' but doesn't calculate it.
# We rename it to 'mt_percent' so we can use it as a covariate later.
# ==============================================================================

# 1. Rename response_p_mito if it exists
if ("response_p_mito" %in% colnames(covariates)) {
    cat("Renaming 'response_p_mito' to 'mt_percent'...\n")
    names(covariates)[names(covariates) == "response_p_mito"] <- "mt_percent"
}

# 2. DROP OTHER RESERVED COLUMNS (which Sceptre DOES calculate)
# We still drop n_umis and n_nonzero because Sceptre re-calculates those automatically.
reserved_cols <- c("response_n_nonzero", "response_n_umis", "grna_n_nonzero", "grna_n_umis")
covariates <- covariates[, !(colnames(covariates) %in% reserved_cols), drop=FALSE]

# 4. CREATE SCEPTRE OBJECT
sceptre_object <- import_data(
  response_matrix = response_mat,
  grna_matrix = grna_mat,
  grna_target_data_frame = guide_mapping,
  extra_covariates = covariates,
  moi = "high"
)

print("✅ Sceptre object created!")
print(head(sceptre_object@covariate_data_frame))
```

### Text output (stream)

```text
Renaming 'response_p_mito' to 'mt_percent'...
[1] "\u2705 Sceptre object created!"
                           response_n_nonzero response_n_umis grna_n_nonzero
GACAGTGCCCTCTATCTCACTTTA_2                693       1109.3333             24
ACAGCAGACTGAGCCACTATTTCA_1                920       1671.6219             97
CCTCCTGAATCCTGTACGCTTAAA_3               1281       2472.5923             26
GTACGCAAATCATTCCCATTTGGA_3               1581       3231.2755             16
CAACCACAAAACATCGGACCTTTC_1                759       1556.6205             96
TGGTGGTAAATGTTGCACTATATA_3                515        858.7499             25
                           grna_n_umis                  batch mt_percent
GACAGTGCCCTCTATCTCACTTTA_2          70 300genes:eTN5_GEX_sub2   9.375001
ACAGCAGACTGAGCCACTATTTCA_1         908 300genes:eTN5_GEX_sub1   9.092965
CCTCCTGAATCCTGTACGCTTAAA_3         121 300genes:eTN5_GEX_sub3  15.483077
GTACGCAAATCATTCCCATTTGGA_3          59 300genes:eTN5_GEX_sub3  12.533751
CAACCACAAAACATCGGACCTTTC_1         380 300genes:eTN5_GEX_sub1  18.148289
TGGTGGTAAATGTTGCACTATATA_3         142 300genes:eTN5_GEX_sub3   9.024746

```

## Code cell 110

```python
%%R
# 1. Load the translation map
gene_map <- read.csv(file.path(data_dir, "gene_symbol_map.csv"))

# 2. Start with guide mapping (which has 'grna_target' as Symbol)
discovery_pairs_highmoi <- guide_mapping

# 3. Merge to attach the correct response_id (Gene ID)
#    This matches 'grna_target' in guides to 'grna_target' in the gene map
discovery_pairs_highmoi <- merge(discovery_pairs_highmoi, gene_map, by="grna_target")

# 4. Select ONLY the two columns Sceptre needs
#    (We filter out non-targeting guides automatically via the merge)
discovery_pairs_highmoi <- discovery_pairs_highmoi[, c("grna_target", "response_id")]

# 5. DEDUPLICATE
#    Sceptre performs analysis at the *Target* level, not the *Guide* level.
#    We only need one row per Target-Gene pair.
discovery_pairs_highmoi <- unique(discovery_pairs_highmoi)

print(paste("Unique pairs to test:", nrow(discovery_pairs_highmoi)))
print(head(discovery_pairs_highmoi))
```

### Text output (stream)

```text
[1] "Unique pairs to test: 287"
   grna_target        response_id
1        ACAA1 ENSG00000060971.19
7       ACVR1B ENSG00000135503.13
13      ACVR2A ENSG00000121989.15
19      ACVR2B ENSG00000114739.14
25       ADCY2 ENSG00000078295.17
31      AMOTL2 ENSG00000114019.15

```

## Code cell 111

```python
%%R 

sceptre_object_highmoi <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs_highmoi,
  side = "both", control_group = "complement")
```

## Code cell 112

```python

```

## Code cell 113

```python
%%R

png(file.path(data_dir, "grna_count_distributions.png"), width = 3000, height = 2000, res = 300)

print(plot_grna_count_distributions(
  sceptre_object = sceptre_object_highmoi,
  n_grnas_to_plot = 8
))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "grna_count_distributions.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d3_sceptre_transfer/grna_count_distributions.png

```

## Code cell 114

```python
%%R

plot_covariates(sceptre_object_highmoi, p_mito_threshold = 0.25)
```

### Text output (text/plain)

```text
<IPython.core.display.Image object>
```

## Code cell 115

```python
%%R
sceptre_object_highmoi <- assign_grnas(
  sceptre_object = sceptre_object_highmoi,
  method = "mixture",
  parallel = FALSE
)
```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178488.23-P1P2 (1 of 1926)
Performing gRNA-to-cell assignments for gRNA ACAA1_+_38178620.23-P1P2 (5 of 1926)
Performing gRNA-to-cell assignments for gRNA ACTL7B_-_111617995.23-P1P2 (10 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR1B_+_52345762.23-P1P2 (15 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2A_+_148602343.23-P1P2 (20 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_+_38495339.23-P1P2 (25 of 1926)
Performing gRNA-to-cell assignments for gRNA ACVR2B_-_38495819.23-P1P2 (30 of 1926)
Performing gRNA-to-cell assignments for gRNA ADCY2_-_7396429.23-P1P2 (35 of 1926)
Performing gRNA-to-cell assignments for gRNA AMOTL2_+_134093279.23-P1P2 (40 of 1926)
Performing gRNA-to-cell assignments for gRNA ANP32E_+_150208342.23-P1P2 (45 of 1926)
Performing gRNA-to-cell assignments for gRNA APAF1_-_99039071.23-P1P2 (50 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_+_72766725.23-P1P2 (55 of 1926)
Performing gRNA-to-cell assignments for gRNA ARIH1_-_72766755.23-P1P2 (60 of 1926)
Performing gRNA-to-cell assignments for gRNA ARTN_-_44399429.23-P1 (65 of 1926)
Performing gRNA-to-cell assignments for gRNA ASCC3_-_101329193.23-P1P2 (70 of 1926)
Performing gRNA-to-cell assignments for gRNA ASH2L_-_37963059.23-P1P2 (75 of 1926)
Performing gRNA-to-cell assignments for gRNA ATMIN_+_81069661.23-P1P2 (80 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_+_111992112.23-P1P2 (85 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP5F1_-_111992129.23-P1P2 (90 of 1926)
Performing gRNA-to-cell assignments for gRNA ATP6AP1_-_153657075.23-P1P2 (95 of 1926)
Performing gRNA-to-cell assignments for gRNA ATXN2L_-_28834668.23-P1P2 (100 of 1926)
Performing gRNA-to-cell assignments for gRNA BAD_+_64052102.23-P1P2 (105 of 1926)
Performing gRNA-to-cell assignments for gRNA BAG6_+_31620398.23-P1P2 (110 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_+_49458268.23-P1P2 (115 of 1926)
Performing gRNA-to-cell assignments for gRNA BAX_-_49458323.23-P1P2 (120 of 1926)
Performing gRNA-to-cell assignments for gRNA BCAS2_-_115124045.23-P1P2 (125 of 1926)
Performing gRNA-to-cell assignments for gRNA BMPR1A_-_88516386.23-P1P2 (130 of 1926)
Performing gRNA-to-cell assignments for gRNA BRD3_-_136932911.23-P1P2 (135 of 1926)
Performing gRNA-to-cell assignments for gRNA C17orf49_+_6918114.23-P2 (140 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_+_31532436.23-P1P2 (145 of 1926)
Performing gRNA-to-cell assignments for gRNA C5orf22_-_31532650.23-P1P2 (150 of 1926)
Performing gRNA-to-cell assignments for gRNA CACFD1_-_136325160.23-P1P2 (155 of 1926)
Performing gRNA-to-cell assignments for gRNA CACNG1_+_65040756.23-P1P2 (160 of 1926)
Performing gRNA-to-cell assignments for gRNA CAPN5_+_76778225.23-P1P2 (165 of 1926)
Performing gRNA-to-cell assignments for gRNA CASP3_+_185570560.23-P1P2 (170 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_+_107384793.23-P1P2 (175 of 1926)
Performing gRNA-to-cell assignments for gRNA CBLL1_-_107384635.23-P1P2 (180 of 1926)
Performing gRNA-to-cell assignments for gRNA CCM2_-_45039941.23-P1P2 (185 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNC_-_100016267.23-P1P2 (190 of 1926)
Performing gRNA-to-cell assignments for gRNA CCNT2_-_135676416.23-P1P2 (195 of 1926)
Performing gRNA-to-cell assignments for gRNA CCT6B_+_33288381.23-P1P2 (200 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_+_30366413.23-P1P2 (205 of 1926)
Performing gRNA-to-cell assignments for gRNA CD2BP2_-_30366630.23-P1P2 (210 of 1926)
Performing gRNA-to-cell assignments for gRNA CD99L2_-_150067120.23-P1P2 (215 of 1926)
Performing gRNA-to-cell assignments for gRNA CDH5_-_66400630.23-P1P2 (220 of 1926)
Performing gRNA-to-cell assignments for gRNA CDKN1A_-_36646538.23-P1P2 (225 of 1926)
Performing gRNA-to-cell assignments for gRNA CENPC_+_68411190.23-P1P2 (230 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_+_30773828.23-P1P2 (235 of 1926)
Performing gRNA-to-cell assignments for gRNA C16orf93_-_30773391.23-P1P2 (240 of 1926)
Performing gRNA-to-cell assignments for gRNA CHMP3_-_86790451.23-P1P2 (245 of 1926)
Performing gRNA-to-cell assignments for gRNA CHORDC1_-_89955894.23-P1P2 (250 of 1926)
Performing gRNA-to-cell assignments for gRNA CNEP1R1_+_50059284.23-P1P2 (255 of 1926)
Performing gRNA-to-cell assignments for gRNA CNOT2_+_70637651.23-P1P2 (260 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_+_83956339.23-P1P2 (265 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS4_-_83956572.23-P1P2 (270 of 1926)
Performing gRNA-to-cell assignments for gRNA COPS6_-_99686676.23-P1P2 (275 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_-_31634088.23-P1P2 (280 of 1926)
Performing gRNA-to-cell assignments for gRNA CTR9_+_10772881.23-P1P2 (285 of 1926)
Performing gRNA-to-cell assignments for gRNA CUL3_-_225449667.23-P1P2 (290 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_+_145150267.23-P1P2 (295 of 1926)
Performing gRNA-to-cell assignments for gRNA CYC1_-_145149996.23-P1P2 (300 of 1926)
Performing gRNA-to-cell assignments for gRNA CYCS_+_25164869.23-P1P2 (305 of 1926)
Performing gRNA-to-cell assignments for gRNA DACH1_-_72441101.23-P1P2 (310 of 1926)
Performing gRNA-to-cell assignments for gRNA DBR1_+_137893673.23-P1P2 (315 of 1926)
Performing gRNA-to-cell assignments for gRNA DCTN2_-_57940910.23-P1P2 (320 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_+_61851532.23-P1P2 (325 of 1926)
Performing gRNA-to-cell assignments for gRNA DDX42_-_61851812.23-P1P2 (330 of 1926)
Performing gRNA-to-cell assignments for gRNA DGCR8_-_20068158.23-P1P2 (335 of 1926)
Performing gRNA-to-cell assignments for gRNA DNAJC17_-_41099215.23-P1P2 (340 of 1926)
Performing gRNA-to-cell assignments for gRNA DNM2_+_10829157.23-P1P2 (345 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM1_+_49575018.23-P1P2 (350 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_+_130700056.23-P1P2 (355 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM2_-_130700076.23-P1P2 (360 of 1926)
Performing gRNA-to-cell assignments for gRNA DPM3_+_155113020.23-P1P2 (365 of 1926)
Performing gRNA-to-cell assignments for gRNA DROSHA_-_31531982.23-P1P2 (370 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_+_3985081.23-P1 (375 of 1926)
Performing gRNA-to-cell assignments for gRNA EEF2_-_3985443.23-P1 (380 of 1926)
Performing gRNA-to-cell assignments for gRNA EHMT2_-_31865199.23-P1P2 (385 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3F_+_8009347.23-P1P2 (390 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF3M_-_32605405.23-P1P2 (395 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4E_+_99849831.23-P1P2 (400 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184032364.23-P1P2 (405 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G1_-_184033042.23-P1P2 (410 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF4G2_+_10830447.23-P1 (415 of 1926)
Performing gRNA-to-cell assignments for gRNA EIF5A_-_7211217.23-P1P2 (420 of 1926)
Performing gRNA-to-cell assignments for gRNA TCEB3_-_24069975.23-P1P2 (425 of 1926)
Performing gRNA-to-cell assignments for gRNA EMC2_+_109456080.23-P1P2 (430 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823220.23-P1P2 (435 of 1926)
Performing gRNA-to-cell assignments for gRNA EME2_-_1823252.23-P1P2 (440 of 1926)
Performing gRNA-to-cell assignments for gRNA EOMES_-_27763342.23-P1P2 (445 of 1926)
Performing gRNA-to-cell assignments for gRNA ERCC2_+_45873802.23-P1P2 (450 of 1926)
Performing gRNA-to-cell assignments for gRNA ERO1L_+_53162425.23-P1P2 (455 of 1926)
Performing gRNA-to-cell assignments for gRNA ETS1_+_128392044.23-P1P2 (460 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_+_36134592.23-P1P2 (465 of 1926)
Performing gRNA-to-cell assignments for gRNA ETV2_-_36134789.23-P1P2 (470 of 1926)
Performing gRNA-to-cell assignments for gRNA EWSR1_-_29664347.23-P1P2 (475 of 1926)
Performing gRNA-to-cell assignments for gRNA EXOSC8_-_37574979.23-P1P2 (480 of 1926)
Performing gRNA-to-cell assignments for gRNA FBXW11_-_171433300.23-P1P2 (485 of 1926)
Performing gRNA-to-cell assignments for gRNA FEN1_+_61560782.23-P1P2 (490 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_+_75532080.23-P1P2 (495 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT11_-_75532461.23-P1P2 (500 of 1926)
Performing gRNA-to-cell assignments for gRNA FUT7_-_139926495.23-P1P2 (505 of 1926)
Performing gRNA-to-cell assignments for gRNA G3BP1_-_151151499.23-P1P2 (510 of 1926)
Performing gRNA-to-cell assignments for gRNA GCK_+_44229021.23-P1P2 (515 of 1926)
Performing gRNA-to-cell assignments for gRNA GOLT1B_+_21654777.23-P1P2 (520 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_+_45000529.23-P1P2 (525 of 1926)
Performing gRNA-to-cell assignments for gRNA GOSR2_-_45000948.23-P1P2 (530 of 1926)
Performing gRNA-to-cell assignments for gRNA GPI_+_34856284.23-P1P2 (535 of 1926)
Performing gRNA-to-cell assignments for gRNA GPR107_-_132816044.23-P1P2 (540 of 1926)
Performing gRNA-to-cell assignments for gRNA GSPT1_-_12009714.23-P1P2 (545 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A1_+_81687215.23-P1P2 (550 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_+_59949355.23-P1P2 (555 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2A2_-_59949657.23-P1P2 (560 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF2H4_-_30876185.23-P1P2 (565 of 1926)
Performing gRNA-to-cell assignments for gRNA GTF3C2_+_27579604.23-P1P2 (570 of 1926)
Performing gRNA-to-cell assignments for gRNA HAND1_+_153857732.23-P1P2 (575 of 1926)
Performing gRNA-to-cell assignments for gRNA HDAC2_-_114292253.23-P1P2 (580 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_+_81399312.23-P1P2 (585 of 1926)
Performing gRNA-to-cell assignments for gRNA HGF_-_81399351.23-P1P2 (590 of 1926)
Performing gRNA-to-cell assignments for gRNA HINT2_+_35814995.23-P1P2 (595 of 1926)
Performing gRNA-to-cell assignments for gRNA HNF4A_+_43029971.23-P1P2 (600 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPM_-_8509885.23-P1 (605 of 1926)
Performing gRNA-to-cell assignments for gRNA HNRNPU_+_245027753.23-P1P2 (610 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_+_43702292.23-P1P2 (615 of 1926)
Performing gRNA-to-cell assignments for gRNA HSD17B12_-_43702305.23-P1P2 (620 of 1926)
Performing gRNA-to-cell assignments for gRNA HSPA5_-_128003534.23-P1P2 (625 of 1926)
Performing gRNA-to-cell assignments for gRNA HUWE1_-_53710540.23-P1P2 (630 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_-_30007588.23-P1P2 (635 of 1926)
Performing gRNA-to-cell assignments for gRNA INTS10_+_19675000.23-P1P2 (640 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_+_65903046.23-P1P2 (645 of 1926)
Performing gRNA-to-cell assignments for gRNA VWA9_-_65903523.23-P1P2 (650 of 1926)
Performing gRNA-to-cell assignments for gRNA C7orf26_+_6630358.23-P1P2 (655 of 1926)
Performing gRNA-to-cell assignments for gRNA IPO8_+_30848771.23-P1P2 (660 of 1926)
Performing gRNA-to-cell assignments for gRNA ITGB1BP1_+_9563543.23-P1P2 (665 of 1926)
Performing gRNA-to-cell assignments for gRNA SSFA2_+_182756743.23-P1P2 (670 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_+_109914946.23-P1P2 (675 of 1926)
Performing gRNA-to-cell assignments for gRNA KCTD10_-_109915112.23-P1P2 (680 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM1A_-_23346409.23-P1P2 (685 of 1926)
Performing gRNA-to-cell assignments for gRNA KDM3B_-_137688431.23-P1P2 (690 of 1926)
Performing gRNA-to-cell assignments for gRNA KDR_-_55991659.23-P1P2 (695 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA0232_+_6784651.23-P1P2 (700 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153706791.23-P1P2 (705 of 1926)
Performing gRNA-to-cell assignments for gRNA LAGE3_+_153707248.23-P1P2 (710 of 1926)
Performing gRNA-to-cell assignments for gRNA LAMTOR4_-_99746921.23-P1P2 (715 of 1926)
Performing gRNA-to-cell assignments for gRNA LARP4_-_50794864.23-P1P2 (720 of 1926)
Performing gRNA-to-cell assignments for gRNA LATS2_-_21635336.23-P1P2 (725 of 1926)
Performing gRNA-to-cell assignments for gRNA LEMD3_+_65563441.23-P1P2 (730 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_+_54960428.23-P1P2 (735 of 1926)
Performing gRNA-to-cell assignments for gRNA LENG8_-_54960616.23-P1P2 (740 of 1926)
Performing gRNA-to-cell assignments for gRNA LHFPL2_-_77944519.23-P1P2 (745 of 1926)
Performing gRNA-to-cell assignments for gRNA LSM12_-_42143603.23-P1P2 (750 of 1926)
Performing gRNA-to-cell assignments for gRNA MANBAL_+_35918136.23-P1P2 (755 of 1926)
Performing gRNA-to-cell assignments for gRNA MAU2_+_19431732.23-P1P2 (760 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_+_47807678.23-P1P2 (765 of 1926)
Performing gRNA-to-cell assignments for gRNA MBD1_-_47807877.23-P1P2 (770 of 1926)
Performing gRNA-to-cell assignments for gRNA MEAF6_-_37980319.23-P1P2 (775 of 1926)
Performing gRNA-to-cell assignments for gRNA MED1_-_37607140.23-P1P2 (780 of 1926)
Performing gRNA-to-cell assignments for gRNA MED12_+_70338541.23-P1P2 (785 of 1926)
Performing gRNA-to-cell assignments for gRNA MED15_+_20862189.23-P1 (790 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_+_57479335.23-P1P2 (795 of 1926)
Performing gRNA-to-cell assignments for gRNA MED19_-_57479437.23-P1P2 (800 of 1926)
Performing gRNA-to-cell assignments for gRNA MEMO1_-_32264882.23-P1 (805 of 1926)
Performing gRNA-to-cell assignments for gRNA MIXL1_+_226411425.23-P1P2 (810 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL28_+_420487.23-P1P2 (815 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL32_+_42972018.23-P1P2 (820 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_+_179322181.23-P1P2 (825 of 1926)
Performing gRNA-to-cell assignments for gRNA MRPL47_-_179322213.23-P1P2 (830 of 1926)
Performing gRNA-to-cell assignments for gRNA LUZP6_-_135612750.23-P1P2 (835 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA20_-_19998023.23-P1P2 (840 of 1926)
Performing gRNA-to-cell assignments for gRNA NAA25_+_112546538.23-P1P2 (845 of 1926)
Performing gRNA-to-cell assignments for gRNA NCBP1_+_100396082.23-P1P2 (850 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_+_172428630.23-P1P2 (855 of 1926)
Performing gRNA-to-cell assignments for gRNA NCEH1_-_172428841.23-P1P2 (860 of 1926)
Performing gRNA-to-cell assignments for gRNA NEDD8_-_24701363.23-P1P2 (865 of 1926)
Performing gRNA-to-cell assignments for gRNA NF2_+_29999635.23-P1P2 (870 of 1926)
Performing gRNA-to-cell assignments for gRNA NFIB_+_14314559.23-P1P2 (875 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_+_139418329.23-P2 (880 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139418349.23-P2 (885 of 1926)
Performing gRNA-to-cell assignments for gRNA NOTCH1_-_139440503.23-P1 (890 of 1926)
Performing gRNA-to-cell assignments for gRNA NUFIP2_+_27621073.23-P1P2 (895 of 1926)
Performing gRNA-to-cell assignments for gRNA NUP93_-_56764090.23-P1P2 (900 of 1926)
Performing gRNA-to-cell assignments for gRNA OGT_+_70752948.23-P1P2 (905 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_+_59383170.23-P1P2 (910 of 1926)
Performing gRNA-to-cell assignments for gRNA OSBP_-_59383348.23-P1P2 (915 of 1926)
Performing gRNA-to-cell assignments for gRNA PAM16_-_4401236.23-P1P2 (920 of 1926)
Performing gRNA-to-cell assignments for gRNA PAXBP1_-_34143614.23-P1P2 (925 of 1926)
Performing gRNA-to-cell assignments for gRNA PBLD_-_70092436.23-P1 (930 of 1926)
Performing gRNA-to-cell assignments for gRNA PCBP2_+_53835614.23-P2 (935 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_+_94292752.23-P2 (940 of 1926)
Performing gRNA-to-cell assignments for gRNA PEG10_-_94293027.23-P2 (945 of 1926)
Performing gRNA-to-cell assignments for gRNA PGRMC2_+_129208913.23-P1P2 (950 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_-_178866211.23-P1P2 (955 of 1926)
Performing gRNA-to-cell assignments for gRNA PITPNC1_-_65373433.23-P1P2 (960 of 1926)
Performing gRNA-to-cell assignments for gRNA PLCG1_+_39766176.23-P1P2 (965 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_+_151042273.23-P2 (970 of 1926)
Performing gRNA-to-cell assignments for gRNA PLEKHG1_-_151042422.23-P2 (975 of 1926)
Performing gRNA-to-cell assignments for gRNA PLRG1_-_155471154.23-P1P2 (980 of 1926)
Performing gRNA-to-cell assignments for gRNA PLXNA2_-_208417275.23-P1P2 (985 of 1926)
Performing gRNA-to-cell assignments for gRNA POLE4_-_75185725.23-P1P2 (990 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR2M_+_57998962.23-P1P2 (995 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940045.23-P1P2 (1000 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3H_+_41940414.23-P1P2 (1005 of 1926)
Performing gRNA-to-cell assignments for gRNA POLR3K_-_103144.23-P1P2 (1010 of 1926)
Performing gRNA-to-cell assignments for gRNA PPHLN1_-_42720007.23-P1P2 (1015 of 1926)
Performing gRNA-to-cell assignments for gRNA PPIL4_+_149867146.23-P1P2 (1020 of 1926)
Performing gRNA-to-cell assignments for gRNA PPP2R1A_+_52693317.23-P1 (1025 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_+_53195217.23-P1P2 (1030 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKCD_-_53195332.23-P1P2 (1035 of 1926)
Performing gRNA-to-cell assignments for gRNA PRKD1_-_30396996.23-P1P2 (1040 of 1926)
Performing gRNA-to-cell assignments for gRNA PRMT5_+_23398554.23-P1 (1045 of 1926)
Performing gRNA-to-cell assignments for gRNA C2orf53_-_27362230.23-P1P2 (1050 of 1926)
Performing gRNA-to-cell assignments for gRNA PRSS8_+_31147010.23-P1P2 (1055 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_+_36909062.23-P1P2 (1060 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMB3_-_36909230.23-P1P2 (1065 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMD4_-_151227289.23-P1P2 (1070 of 1926)
Performing gRNA-to-cell assignments for gRNA PSMG4_-_3259227.23-P1P2 (1075 of 1926)
Performing gRNA-to-cell assignments for gRNA PTAR1_+_72374835.23-P1P2 (1080 of 1926)
Performing gRNA-to-cell assignments for gRNA PTBP1_+_797463.23-P1P2 (1085 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_+_89623388.23-P1P2 (1090 of 1926)
Performing gRNA-to-cell assignments for gRNA PTEN_-_89623526.23-P1P2 (1095 of 1926)
Performing gRNA-to-cell assignments for gRNA PTMS_-_6875591.23-P1P2 (1100 of 1926)
Performing gRNA-to-cell assignments for gRNA PYGO2_-_154934135.23-P1P2 (1105 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB18_-_27793292.23-P1P2 (1110 of 1926)
Performing gRNA-to-cell assignments for gRNA RAB6A_+_73471983.23-P1P2 (1115 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_+_76251918.23-P1P2 (1120 of 1926)
Performing gRNA-to-cell assignments for gRNA RABGGTB_-_76252008.23-P1P2 (1125 of 1926)
Performing gRNA-to-cell assignments for gRNA GNB2L1_-_180670840.23-P1P2 (1130 of 1926)
Performing gRNA-to-cell assignments for gRNA RAD21_-_117886635.23-P1P2 (1135 of 1926)
Performing gRNA-to-cell assignments for gRNA RANGRF_+_8192257.23-P1P2 (1140 of 1926)
Performing gRNA-to-cell assignments for gRNA RBBP5_+_205091050.23-P1P2 (1145 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM14-RBM4_+_66406280.23-P2 (1150 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM33_+_155437802.23-P1P2 (1155 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM4_+_66406319.23-P2 (1160 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_+_114271352.23-P1P2 (1165 of 1926)
Performing gRNA-to-cell assignments for gRNA RBM7_-_114271700.23-P1P2 (1170 of 1926)
Performing gRNA-to-cell assignments for gRNA RBX1_-_41347852.23-P1P2 (1175 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1551_+_32112399.23-P1P2 (1180 of 1926)
Performing gRNA-to-cell assignments for gRNA RFX3_-_3525750.23-P1P2 (1185 of 1926)
Performing gRNA-to-cell assignments for gRNA RGP1_+_35749328.23-P1P2 (1190 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_+_95067102.23-P1P2 (1195 of 1926)
Performing gRNA-to-cell assignments for gRNA RHOBTB3_-_95067164.23-P1P2 (1200 of 1926)
Performing gRNA-to-cell assignments for gRNA RIBC1_-_53449606.23-P1P2 (1205 of 1926)
Performing gRNA-to-cell assignments for gRNA KIAA1432_+_5629173.23-P1P2 (1210 of 1926)
Performing gRNA-to-cell assignments for gRNA RLF_+_40627522.23-P1P2 (1215 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF115_-_145611148.23-P1P2 (1220 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_+_104296142.23-P1P2 (1225 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF20_-_104296250.23-P1P2 (1230 of 1926)
Performing gRNA-to-cell assignments for gRNA RNF40_-_30773681.23-P1P2 (1235 of 1926)
Performing gRNA-to-cell assignments for gRNA RPL37A_-_217363600.23-P1P2 (1240 of 1926)
Performing gRNA-to-cell assignments for gRNA RPS15A_-_18801590.23-P1P2 (1245 of 1926)
Performing gRNA-to-cell assignments for gRNA SAV1_+_51134992.23-P1P2 (1250 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_+_33104038.23-P1P2 (1255 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF4_-_33104343.23-P1P2 (1260 of 1926)
Performing gRNA-to-cell assignments for gRNA SCAF8_+_155054976.23-P1P2 (1265 of 1926)
Performing gRNA-to-cell assignments for gRNA EPT1_+_26569050.23-P1P2 (1270 of 1926)
Performing gRNA-to-cell assignments for gRNA SERF2_+_44069659.23-P2 (1275 of 1926)
Performing gRNA-to-cell assignments for gRNA SERPING1_+_57365137.23-P1P2 (1280 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_+_122242041.23-P1P2 (1285 of 1926)
Performing gRNA-to-cell assignments for gRNA SETD1B_-_122242294.23-P1P2 (1290 of 1926)
Performing gRNA-to-cell assignments for gRNA SETDB1_-_150898815.23-P1P2 (1295 of 1926)
Performing gRNA-to-cell assignments for gRNA SF1_-_64545937.23-P1P2 (1300 of 1926)
Performing gRNA-to-cell assignments for gRNA CCDC101_-_28565268.23-P1P2 (1305 of 1926)
Performing gRNA-to-cell assignments for gRNA SKA3_+_21750519.23-P1P2 (1310 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_+_78591497.23-P1P2 (1315 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_-_78591737.23-P1P2 (1320 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC29A1_-_44187438.23-P2 (1325 of 1926)
Performing gRNA-to-cell assignments for gRNA SLC35C1_-_45826711.23-P1P2 (1330 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD1_+_146403428.23-P1P2 (1335 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD3_+_67358229.23-P1P2 (1340 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_+_48556623.23-P1P2 (1345 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD4_-_48556657.23-P1P2 (1350 of 1926)
Performing gRNA-to-cell assignments for gRNA SMAD5_-_135468694.23-P1P2 (1355 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA4_-_11071741.23-P1P2 (1360 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCA5_+_144435243.23-P1P2 (1365 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCB1_+_24129228.23-P1 (1370 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_+_47822975.23-P1P2 (1375 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCC1_-_47822893.23-P1P2 (1380 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCD3_-_150974212.23-P1 (1385 of 1926)
Performing gRNA-to-cell assignments for gRNA SMARCE1_-_38803669.23-P1P2 (1390 of 1926)
Performing gRNA-to-cell assignments for gRNA SMC2_-_106856631.23-P1P2 (1395 of 1926)
Performing gRNA-to-cell assignments for gRNA SMG5_+_156252537.23-P1P2 (1400 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_+_103852.23-P1P2 (1405 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP25_-_104242.23-P1P2 (1410 of 1926)
Performing gRNA-to-cell assignments for gRNA SNRNP27_-_70121181.23-P1P2 (1415 of 1926)
Performing gRNA-to-cell assignments for gRNA SNUPN_-_75917960.23-P1P2 (1420 of 1926)
Performing gRNA-to-cell assignments for gRNA SOX4_+_21594213.23-P1 (1425 of 1926)
Performing gRNA-to-cell assignments for gRNA SPIN1_+_91003737.23-P1P2 (1430 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_+_2328680.23-P1P2 (1435 of 1926)
Performing gRNA-to-cell assignments for gRNA SPPL2B_-_2328912.23-P1P2 (1440 of 1926)
Performing gRNA-to-cell assignments for gRNA SPTY2D1_-_18655846.23-P1P2 (1445 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP19_+_112197188.23-P1P2 (1450 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP54_-_35452244.23-P1P2 (1455 of 1926)
Performing gRNA-to-cell assignments for gRNA SRP68_+_74068540.23-P1P2 (1460 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_+_126138730.23-P1P2 (1465 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPR_-_126138742.23-P1P2 (1470 of 1926)
Performing gRNA-to-cell assignments for gRNA SRPRB_-_133524719.23-P1P2 (1475 of 1926)
Performing gRNA-to-cell assignments for gRNA SRRM2_-_2802714.23-P1P2 (1480 of 1926)
Performing gRNA-to-cell assignments for gRNA SRSF6_+_42086993.23-P1 (1485 of 1926)
Performing gRNA-to-cell assignments for gRNA SSB_+_170655450.23-P1P2 (1490 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_+_57102902.23-P1P2 (1495 of 1926)
Performing gRNA-to-cell assignments for gRNA SSRP1_-_57103355.23-P1P2 (1500 of 1926)
Performing gRNA-to-cell assignments for gRNA STAG2_-_123095658.23-P1 (1505 of 1926)
Performing gRNA-to-cell assignments for gRNA STK40_+_36851507.23-P1P2 (1510 of 1926)
Performing gRNA-to-cell assignments for gRNA SUPT20H_+_37633798.23-P1P2 (1515 of 1926)
Performing gRNA-to-cell assignments for gRNA SUZ12_+_30264157.23-P1P2 (1520 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_+_64901973.23-P1P2 (1525 of 1926)
Performing gRNA-to-cell assignments for gRNA SYVN1_-_64901981.23-P1P2 (1530 of 1926)
Performing gRNA-to-cell assignments for gRNA TADA2B_-_7045121.23-P1P2 (1535 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF10_+_6633406.23-P1P2 (1540 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF2_+_120845032.23-P1P2 (1545 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF3_-_7860540.23-P1P2 (1550 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_+_60640814.23-P1P2 (1555 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF4_-_60640498.23-P1P2 (1560 of 1926)
Performing gRNA-to-cell assignments for gRNA TAF6L_-_62538954.23-P1P2 (1565 of 1926)
Performing gRNA-to-cell assignments for gRNA TAOK1_-_27717907.23-P1P2 (1570 of 1926)
Performing gRNA-to-cell assignments for gRNA FAM208A_+_56717273.23-P1P2 (1575 of 1926)
Performing gRNA-to-cell assignments for gRNA TBC1D9_+_141677481.23-P1P2 (1580 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_+_170863524.23-P1P2 (1585 of 1926)
Performing gRNA-to-cell assignments for gRNA TBP_-_170863481.23-P1P2 (1590 of 1926)
Performing gRNA-to-cell assignments for gRNA TBX3_-_115121998.23-P1P2 (1595 of 1926)
Performing gRNA-to-cell assignments for gRNA T_-_166581333.23-P1P2 (1600 of 1926)
Performing gRNA-to-cell assignments for gRNA TGFBR2_+_30648526.23-P1P2 (1605 of 1926)
Performing gRNA-to-cell assignments for gRNA THBD_+_23029196.23-P2 (1610 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_+_2977592.23-P1P2 (1615 of 1926)
Performing gRNA-to-cell assignments for gRNA TLE6_-_2977643.23-P1P2 (1620 of 1926)
Performing gRNA-to-cell assignments for gRNA MESDC1_-_81293340.23-P1P2 (1625 of 1926)
Performing gRNA-to-cell assignments for gRNA TMEM41B_+_9336156.23-P1P2 (1630 of 1926)
Performing gRNA-to-cell assignments for gRNA TMTC3_-_88536127.23-P1P2 (1635 of 1926)
Performing gRNA-to-cell assignments for gRNA TOP1_+_39657538.23-P1P2 (1640 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_+_7590531.23-P1P2 (1645 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_-_7590704.23-P1P2 (1650 of 1926)
Performing gRNA-to-cell assignments for gRNA TRIM24_-_138145308.23-P1P2 (1655 of 1926)
Performing gRNA-to-cell assignments for gRNA TRRAP_-_98476465.23-P1P2 (1660 of 1926)
Performing gRNA-to-cell assignments for gRNA TSTD2_+_100395594.23-P1P2 (1665 of 1926)
Performing gRNA-to-cell assignments for gRNA UBE2E1_+_23847541.23-P1P2 (1670 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_+_38924011.23-P1P2 (1675 of 1926)
Performing gRNA-to-cell assignments for gRNA UFM1_-_38924421.23-P1P2 (1680 of 1926)
Performing gRNA-to-cell assignments for gRNA UPF2_+_12084799.23-P1P2 (1685 of 1926)
Performing gRNA-to-cell assignments for gRNA UXT_-_47518505.23-P1P2 (1690 of 1926)
Performing gRNA-to-cell assignments for gRNA VHL_-_10183525.23-P1P2 (1695 of 1926)
Performing gRNA-to-cell assignments for gRNA VPS54_+_64246104.23-P1P2 (1700 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_+_28821698.23-P1P2 (1705 of 1926)
Performing gRNA-to-cell assignments for gRNA WAC_-_28822074.23-P1P2 (1710 of 1926)
Performing gRNA-to-cell assignments for gRNA WAPAL_+_88281533.23-P1P2 (1715 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR12_-_203776383.23-P1P2 (1720 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR77_+_111991860.23-P1P2 (1725 of 1926)
Performing gRNA-to-cell assignments for gRNA WIZ_-_15543016.23-P1P2 (1730 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_+_9983328.23-P1P2 (1735 of 1926)
Performing gRNA-to-cell assignments for gRNA WWC3_-_9983510.23-P1P2 (1740 of 1926)
Performing gRNA-to-cell assignments for gRNA XPO5_-_43543670.23-P1P2 (1745 of 1926)
Performing gRNA-to-cell assignments for gRNA XRN1_-_142166648.23-P1P2 (1750 of 1926)
Performing gRNA-to-cell assignments for gRNA YARS2_+_32908787.23-P1P2 (1755 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB22_+_33285454.23-P1P2 (1760 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_+_6640093.23-P1P2 (1765 of 1926)
Performing gRNA-to-cell assignments for gRNA ZBTB48_-_6640138.23-P1P2 (1770 of 1926)
Performing gRNA-to-cell assignments for gRNA ZC3H8_-_113012435.23-P1P2 (1775 of 1926)
Performing gRNA-to-cell assignments for gRNA ZNF281_+_200379151.23-P1P2 (1780 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00117 (1785 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00233 (1790 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00518 (1795 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_00656 (1800 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01226 (1805 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_01730 (1810 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02009 (1815 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02288 (1820 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_02777 (1825 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03026 (1830 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03252 (1835 of 1926)
Performing gRNA-to-cell assignments for gRNA non-targeting_03620 (1840 of 1926)
Performing gRNA-to-cell assignments for gRNA CSNK2B_SGTN_04 (1845 of 1926)
Performing gRNA-to-cell assignments for gRNA INO80E_SGTN_25 (1850 of 1926)
Performing gRNA-to-cell assignments for gRNA PIK3CA_SGTN_26 (1855 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1044 (1860 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1414 (1865 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1554 (1870 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_1785 (1875 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2009 (1880 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2271 (1885 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2425 (1890 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_2614 (1895 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_293 (1900 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_541 (1905 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_725 (1910 of 1926)
Performing gRNA-to-cell assignments for gRNA SAFE_TARGETING_860 (1915 of 1926)
Performing gRNA-to-cell assignments for gRNA WDR61_SGTN_14 (1920 of 1926)
Performing gRNA-to-cell assignments for gRNA TP53_SGTN_01 (1925 of 1926)

```

## Code cell 116

```python
%%R
library(Matrix)

# 1. Get the assignment matrix
A <- get_grna_assignments(sceptre_object_highmoi)

# 2. Ensure it is a numeric sparse matrix (1/0)
if (inherits(A, "Matrix") || is.matrix(A)) {
    A_num <- if (inherits(A, "Matrix")) drop0(A * 1) else (A * 1)
} else {
    stop("Output is not a matrix.")
}

# --- FIX: ROBUST NAME RECOVERY ---
# We retrieve names using Sceptre's standard metadata accessors
# instead of touching internal slots like @grna_matrix.

# Try to get cell names from the matrix first
current_cell_names <- colnames(A_num)
current_guide_names <- rownames(A_num)

# If missing, grab them from the main object's covariate data frame
if (is.null(current_cell_names)) {
    message("⚠️ Matrix missing colnames. Fetching from sceptre_object...")
    # The covariate data frame ALWAYS has the cell IDs as rownames
    current_cell_names <- rownames(sceptre_object_highmoi@covariate_data_frame)
    
    # Safety check: lengths must match
    if (length(current_cell_names) != ncol(A_num)) {
        stop(paste("Mismatch! Matrix has", ncol(A_num), "cols but object has", length(current_cell_names), "cells."))
    }
    colnames(A_num) <- current_cell_names
}

# If guide names are missing, they should match the unique values in your guide mapping
if (is.null(current_guide_names)) {
    message("⚠️ Matrix missing rownames. Fetching from guide mapping...")
    # Get unique targets from the mapping table you loaded earlier
    # Note: This assumes the matrix rows are sorted/aligned to these unique targets.
    # In 'mixture' method, rows usually correspond to the 'grna_target' labels.
    
    # SAFE FALLBACK: If we can't be 100% sure of order, we warn the user.
    # But usually, they match the 'grna_target' column from your input.
    unique_targets <- sort(unique(guide_mapping$grna_target))
    
    if (length(unique_targets) == nrow(A_num)) {
         rownames(A_num) <- unique_targets
    } else {
         warning("Could not automatically align target names. Using generic IDs.")
         rownames(A_num) <- paste0("target_", seq_len(nrow(A_num)))
    }
}
# ---------------------------------

# 3. Convert to Long Format
s <- summary(A_num)
assign_long <- data.frame(
    cell_id = colnames(A_num)[s$j],
    grna_target = rownames(A_num)[s$i],
    value = s$x
)

# 4. Save
out_path <- file.path(data_dir, "grna_assignments_long.csv")
write.csv(assign_long, out_path, row.names = FALSE)

message("✅ Saved long CSV to: ", out_path)
print(head(assign_long))
```

### Text output (stream)

```text
                     cell_id               grna_target value
1 GACAGTGCCCTCTATCTCACTTTA_2    BAD_+_64052144.23-P1P2     1
2 GACAGTGCCCTCTATCTCACTTTA_2   DNM2_-_10828821.23-P1P2     1
3 GACAGTGCCCTCTATCTCACTTTA_2  TCEB3_-_24070127.23-P1P2     1
4 GACAGTGCCCTCTATCTCACTTTA_2    GPI_+_34856284.23-P1P2     1
5 GACAGTGCCCTCTATCTCACTTTA_2            CSNK2B_SGTN_05     1
6 ACAGCAGACTGAGCCACTATTTCA_1 SMAD1_+_146403220.23-P1P2     1

```

### Text output (text/plain)

```text
⚠️ Matrix missing colnames. Fetching from sceptre_object...
✅ Saved long CSV to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d3_sceptre_transfer/grna_assignments_long.csv

```

## Code cell 117

```python
%%R
png(file.path(data_dir, "MOI.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "MOI.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d3_sceptre_transfer/MOI.png

```

## Code cell 118

```python
%%R
sceptre_object_highmoi <- run_qc(
  sceptre_object = sceptre_object_highmoi,
  p_mito_threshold = 0.25
)

png(file.path(data_dir, "QC.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "QC.png"))
```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d3_sceptre_transfer/QC.png

```

## Code cell 119

```python
%%R
# Calibration check


sceptre_object_highmoi <- run_calibration_check(
  sceptre_object = sceptre_object_highmoi,
  parallel = FALSE
)

png(file.path(data_dir, "calibration_check.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "calibration_check.png"))

```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Constructing negative control pairs. ✓
Running precomputation on response ENSG00000163346.17 (1 of 233)
Running precomputation on response ENSG00000253767.3 (5 of 233)
Running precomputation on response ENSG00000234287.1 (10 of 233)
Running precomputation on response ENSG00000115808.12 (15 of 233)
Running precomputation on response ENSG00000136699.20 (20 of 233)
Running precomputation on response ENSG00000117448.14 (25 of 233)
Running precomputation on response ENSG00000264324.1 (30 of 233)
Running precomputation on response ENSG00000113456.20 (35 of 233)
Running precomputation on response ENSG00000178971.17 (40 of 233)
Running precomputation on response ENSG00000177663.15 (45 of 233)
Running precomputation on response ENSG00000171204.13 (50 of 233)
Running precomputation on response ENSG00000131002.15 (55 of 233)
Running precomputation on response ENSG00000165061.15 (60 of 233)
Running precomputation on response ENSG00000181827.16 (65 of 233)
Running precomputation on response ENSG00000291064.1 (70 of 233)
Running precomputation on response ENSG00000290192.1 (75 of 233)
Running precomputation on response ENSG00000171132.14 (80 of 233)
Running precomputation on response ENSG00000067225.21 (85 of 233)
Running precomputation on response ENSG00000007129.18 (90 of 233)
Running precomputation on response ENSG00000177156.11 (95 of 233)
Running precomputation on response ENSG00000156482.11 (100 of 233)
Running precomputation on response ENSG00000258526.8 (105 of 233)
Running precomputation on response ENSG00000230409.3 (110 of 233)
Running precomputation on response ENSG00000047410.14 (115 of 233)
Running precomputation on response ENSG00000196839.14 (120 of 233)
Running precomputation on response ENSG00000106624.11 (125 of 233)
Running precomputation on response ENSG00000163382.13 (130 of 233)
Running precomputation on response ENSG00000126821.8 (135 of 233)
Running precomputation on response ENSG00000164823.11 (140 of 233)
Running precomputation on response ENSG00000273748.1 (145 of 233)
Running precomputation on response ENSG00000181264.9 (150 of 233)
Running precomputation on response ENSG00000128487.19 (155 of 233)
Running precomputation on response ENSG00000172183.16 (160 of 233)
Running precomputation on response ENSG00000186260.17 (165 of 233)
Running precomputation on response ENSG00000137942.18 (170 of 233)
Running precomputation on response ENSG00000175155.11 (175 of 233)
Running precomputation on response ENSG00000136997.22 (180 of 233)
Running precomputation on response ENSG00000137992.15 (185 of 233)
Running precomputation on response ENSG00000189120.5 (190 of 233)
Running precomputation on response ENSG00000255245.4 (195 of 233)
Running precomputation on response ENSG00000086289.13 (200 of 233)
Running precomputation on response ENSG00000182359.15 (205 of 233)
Running precomputation on response ENSG00000189058.9 (210 of 233)
Running precomputation on response ENSG00000186834.4 (215 of 233)
Running precomputation on response ENSG00000116138.13 (220 of 233)
Running precomputation on response ENSG00000101230.6 (225 of 233)
Running precomputation on response ENSG00000171724.3 (230 of 233)
Analyzing pairs containing gRNA group non-targeting_00184&non-targeting_00259&non-targeting_02159&non-targeting_02571&non-targeting_03029&non-targeting_03081 (1 of 91)
Analyzing pairs containing gRNA group non-targeting_01226&non-targeting_01509&non-targeting_01732&non-targeting_02543&non-targeting_02837&non-targeting_03026 (5 of 91)
Analyzing pairs containing gRNA group non-targeting_00518&non-targeting_03081&non-targeting_03342&non-targeting_03427&non-targeting_03635&non-targeting_03664 (10 of 91)
Analyzing pairs containing gRNA group non-targeting_01509&non-targeting_01812&non-targeting_02025&non-targeting_02677&non-targeting_02777&non-targeting_03384 (15 of 91)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_00324&non-targeting_00390&non-targeting_00629&non-targeting_01934&non-targeting_03180 (20 of 91)
Analyzing pairs containing gRNA group non-targeting_00518&non-targeting_00656&non-targeting_01730&non-targeting_02574&non-targeting_03029&non-targeting_03158 (25 of 91)
Analyzing pairs containing gRNA group non-targeting_00213&non-targeting_02543&non-targeting_02777&non-targeting_02814&non-targeting_03252&non-targeting_03483 (30 of 91)
Analyzing pairs containing gRNA group non-targeting_00196&non-targeting_00918&non-targeting_01730&non-targeting_02025&non-targeting_02911&non-targeting_03384 (35 of 91)
Analyzing pairs containing gRNA group non-targeting_00233&non-targeting_00390&non-targeting_01601&non-targeting_02009&non-targeting_02288&non-targeting_03029 (40 of 91)
Analyzing pairs containing gRNA group non-targeting_00014&non-targeting_02159&non-targeting_02193&non-targeting_02677&non-targeting_02962&non-targeting_03029 (45 of 91)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00390&non-targeting_01732&non-targeting_02193&non-targeting_03081&non-targeting_03384 (50 of 91)
Analyzing pairs containing gRNA group non-targeting_00629&non-targeting_01185&non-targeting_01732&non-targeting_02089&non-targeting_02677&non-targeting_03483 (55 of 91)
Analyzing pairs containing gRNA group non-targeting_01185&non-targeting_02543&non-targeting_02677&non-targeting_03029&non-targeting_03342&non-targeting_03664 (60 of 91)
Analyzing pairs containing gRNA group non-targeting_00614&non-targeting_01185&non-targeting_01226&non-targeting_02089&non-targeting_02159&non-targeting_02962 (65 of 91)
Analyzing pairs containing gRNA group non-targeting_00355&non-targeting_00656&non-targeting_01730&non-targeting_02574&non-targeting_02911&non-targeting_03664 (70 of 91)
Analyzing pairs containing gRNA group non-targeting_00355&non-targeting_00518&non-targeting_00555&non-targeting_00794&non-targeting_01185&non-targeting_02159 (75 of 91)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00196&non-targeting_00213&non-targeting_00518&non-targeting_02574&non-targeting_02837 (80 of 91)
Analyzing pairs containing gRNA group non-targeting_00614&non-targeting_01599&non-targeting_02009&non-targeting_02543&non-targeting_02837&non-targeting_03081 (85 of 91)
Analyzing pairs containing gRNA group non-targeting_00117&non-targeting_00555&non-targeting_01226&non-targeting_01812&non-targeting_02911&non-targeting_03180 (90 of 91)

```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d3_sceptre_transfer/calibration_check.png
In addition: There were 50 or more warnings (use warnings() to see the first 50)

```

## Code cell 120

```python
%%R
calibration_result <- get_result(
  sceptre_object = sceptre_object_highmoi,
  analysis = "run_calibration_check"
)
head(calibration_result)
```

### Text output (stream)

```text
          response_id
               <char>
1: ENSG00000175984.15
2: ENSG00000174695.10
3: ENSG00000171570.11
4:  ENSG00000278390.7
5: ENSG00000137992.15
6:  ENSG00000160213.9
                                                                                                               grna_target
                                                                                                                    <char>
1: non-targeting_00233&non-targeting_02193&non-targeting_02571&non-targeting_02911&non-targeting_03081&non-targeting_03180
2: non-targeting_00117&non-targeting_01601&non-targeting_01732&non-targeting_02962&non-targeting_03483&non-targeting_03664
3: non-targeting_00117&non-targeting_00656&non-targeting_01601&non-targeting_03158&non-targeting_03620&non-targeting_03635
4: non-targeting_00233&non-targeting_00972&non-targeting_02193&non-targeting_02574&non-targeting_02777&non-targeting_02962
5: non-targeting_00184&non-targeting_02193&non-targeting_02543&non-targeting_02837&non-targeting_03026&non-targeting_03342
6: non-targeting_00014&non-targeting_00518&non-targeting_00614&non-targeting_03026&non-targeting_03029&non-targeting_03384
   n_nonzero_trt n_nonzero_cntrl pass_qc      p_value log_2_fold_change
           <int>           <int>  <lgcl>        <num>             <num>
1:            20             371    TRUE 0.0005481779         1.2263258
2:            54            3076    TRUE 0.0076809527         0.5482496
3:            32            2086    TRUE 0.0089332634         0.7373303
4:            21             420    TRUE 0.0102692754         0.8778109
5:            17            1575    TRUE 0.0320000000        -0.6913136
6:            64            5269    TRUE 0.0320000000        -0.3982145
   significant
        <lgcl>
1:       FALSE
2:       FALSE
3:       FALSE
4:       FALSE
5:       FALSE
6:       FALSE

```

## Code cell 121

```python
%%R
sceptre_object_highmoi <- run_discovery_analysis(
  sceptre_object = sceptre_object_highmoi,
  parallel = FALSE
)

png(file.path(data_dir, "discovery.png"), width = 2400, height = 2000, res = 300)

print(plot(sceptre_object_highmoi))

dev.off()

message("✅ Saved PNG plot to: ", file.path(data_dir, "discovery.png"))
```

### Text output (stream)

```text
Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.

Running precomputation on response ENSG00000004487.18 (1 of 238)
Running precomputation on response ENSG00000117000.9 (5 of 238)
Running precomputation on response ENSG00000265491.5 (10 of 238)
Running precomputation on response ENSG00000198952.8 (15 of 238)
Running precomputation on response ENSG00000119185.13 (20 of 238)
Running precomputation on response ENSG00000124380.11 (25 of 238)
Running precomputation on response ENSG00000138434.17 (30 of 238)
Running precomputation on response ENSG00000163513.19 (35 of 238)
Running precomputation on response ENSG00000138231.14 (40 of 238)
Running precomputation on response ENSG00000170871.12 (45 of 238)
Running precomputation on response ENSG00000164040.18 (50 of 238)
Running precomputation on response ENSG00000164305.20 (55 of 238)
Running precomputation on response ENSG00000164292.13 (60 of 238)
Running precomputation on response ENSG00000113196.3 (65 of 238)
Running precomputation on response ENSG00000204435.15 (70 of 238)
Running precomputation on response ENSG00000112237.13 (75 of 238)
Running precomputation on response ENSG00000213079.10 (80 of 238)
Running precomputation on response ENSG00000196367.15 (85 of 238)
Running precomputation on response ENSG00000082014.17 (90 of 238)
Running precomputation on response ENSG00000064313.13 (95 of 238)
Running precomputation on response ENSG00000107185.10 (100 of 238)
Running precomputation on response ENSG00000155827.12 (105 of 238)
Running precomputation on response ENSG00000169925.17 (110 of 238)
Running precomputation on response ENSG00000095787.26 (115 of 238)
Running precomputation on response ENSG00000175390.15 (120 of 238)
Running precomputation on response ENSG00000149100.13 (125 of 238)
Running precomputation on response ENSG00000110048.12 (130 of 238)
Running precomputation on response ENSG00000173933.21 (135 of 238)
Running precomputation on response ENSG00000134954.14 (140 of 238)
Running precomputation on response ENSG00000134283.18 (145 of 238)
Running precomputation on response ENSG00000174106.3 (150 of 238)
Running precomputation on response ENSG00000111300.10 (155 of 238)
Running precomputation on response ENSG00000102710.21 (160 of 238)
Running precomputation on response ENSG00000184304.17 (165 of 238)
Running precomputation on response ENSG00000140264.21 (170 of 238)
Running precomputation on response ENSG00000140395.9 (175 of 238)
Running precomputation on response ENSG00000134419.15 (180 of 238)
Running precomputation on response ENSG00000102900.13 (185 of 238)
Running precomputation on response ENSG00000141510.18 (190 of 238)
Running precomputation on response ENSG00000125686.12 (195 of 238)
Running precomputation on response ENSG00000154217.16 (200 of 238)
Running precomputation on response ENSG00000005206.17 (205 of 238)
Running precomputation on response ENSG00000011451.21 (210 of 238)
Running precomputation on response ENSG00000167615.17 (215 of 238)
Running precomputation on response ENSG00000124193.16 (220 of 238)
Running precomputation on response ENSG00000128191.17 (225 of 238)
Running precomputation on response ENSG00000100387.9 (230 of 238)
Running precomputation on response ENSG00000147162.15 (235 of 238)
Analyzing pairs containing gRNA group KDM1A (1 of 238)
Analyzing pairs containing gRNA group RLF (5 of 238)
Analyzing pairs containing gRNA group RNF115 (10 of 238)
Analyzing pairs containing gRNA group SMG5 (15 of 238)
Analyzing pairs containing gRNA group ITGB1BP1 (20 of 238)
Analyzing pairs containing gRNA group SNRNP27 (25 of 238)
Analyzing pairs containing gRNA group ITPRID2 (30 of 238)
Analyzing pairs containing gRNA group TGFBR2 (35 of 238)
Analyzing pairs containing gRNA group DBR1 (40 of 238)
Analyzing pairs containing gRNA group KIAA0232 (45 of 238)
Analyzing pairs containing gRNA group PGRMC2 (50 of 238)
Analyzing pairs containing gRNA group CASP3 (55 of 238)
Analyzing pairs containing gRNA group RHOBTB3 (60 of 238)
Analyzing pairs containing gRNA group HAND1 (65 of 238)
Analyzing pairs containing gRNA group CSNK2B (70 of 238)
Analyzing pairs containing gRNA group CCNC (75 of 238)
Analyzing pairs containing gRNA group SCAF8 (80 of 238)
Analyzing pairs containing gRNA group TRRAP (85 of 238)
Analyzing pairs containing gRNA group SMARCD3 (90 of 238)
Analyzing pairs containing gRNA group TAF2 (95 of 238)
Analyzing pairs containing gRNA group RGP1 (100 of 238)
Analyzing pairs containing gRNA group RNF20 (105 of 238)
Analyzing pairs containing gRNA group BRD3 (110 of 238)
Analyzing pairs containing gRNA group WAC (115 of 238)
Analyzing pairs containing gRNA group EIF3F (120 of 238)
Analyzing pairs containing gRNA group EIF3M (125 of 238)
Analyzing pairs containing gRNA group OSBP (130 of 238)
Analyzing pairs containing gRNA group RBM4 (135 of 238)
Analyzing pairs containing gRNA group ETS1 (140 of 238)
Analyzing pairs containing gRNA group PPHLN1 (145 of 238)
Analyzing pairs containing gRNA group LEMD3 (150 of 238)
Analyzing pairs containing gRNA group NAA25 (155 of 238)
Analyzing pairs containing gRNA group SUPT20H (160 of 238)
Analyzing pairs containing gRNA group PRKD1 (165 of 238)
Analyzing pairs containing gRNA group SERF2 (170 of 238)
Analyzing pairs containing gRNA group SKIC8 (175 of 238)
Analyzing pairs containing gRNA group RPS15A (180 of 238)
Analyzing pairs containing gRNA group NUP93 (185 of 238)
Analyzing pairs containing gRNA group TP53 (190 of 238)
Analyzing pairs containing gRNA group MED1 (195 of 238)
Analyzing pairs containing gRNA group PITPNC1 (200 of 238)
Analyzing pairs containing gRNA group SPPL2B (205 of 238)
Analyzing pairs containing gRNA group WIZ (210 of 238)
Analyzing pairs containing gRNA group LENG8 (215 of 238)
Analyzing pairs containing gRNA group SRSF6 (220 of 238)
Analyzing pairs containing gRNA group DGCR8 (225 of 238)
Analyzing pairs containing gRNA group RBX1 (230 of 238)
Analyzing pairs containing gRNA group OGT (235 of 238)

```

### Text output (text/plain)

```text
✅ Saved PNG plot to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d3_sceptre_transfer/discovery.png
In addition: There were 50 or more warnings (use warnings() to see the first 50)

```

## Code cell 122

```python
%%R
discovery_result <- get_result(
  sceptre_object = sceptre_object_highmoi,
  analysis = "run_discovery_analysis"
)

write.csv(discovery_result, file.path(data_dir, "TargetKD.csv"), row.names = FALSE)
message("✅ Saved discovery to: ", file.path(data_dir, "TargetKD.csv"))

```

### Text output (text/plain)

```text
✅ Saved discovery to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d3_sceptre_transfer/TargetKD.csv

```

## Markdown cell 123

### Plotting

