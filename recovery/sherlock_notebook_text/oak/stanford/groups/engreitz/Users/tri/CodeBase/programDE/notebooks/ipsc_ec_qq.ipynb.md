# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/ipsc_ec_qq.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/ipsc_ec_qq.ipynb`
- SHA-256: `b5662ed051311cfd0f3aee77a558e4d54b24afdff836148e3cfec9b69c10f553`
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

```

## Code cell 3

```python
import sys
sys.path.append("/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE")
from tests.synthetic_data import make_sceptre_style_synth
from src.sceptre import prepare_crt_inputs
```

## Code cell 4

```python
data_dir = "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data"
working_dir = "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks"
```

## Code cell 5

```python
mdata = mu.read(f"{data_dir}/d2_normalized_gex_and_guide.h5mu")
```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)

```

## Code cell 6

```python
adata = mdata['GEX'].copy()
```

## Markdown cell 7

#### Import cNMF matrix

## Code cell 8

```python
# 1. Define the file path
#cnmf_file = "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012326_100k_cells_20iter_torch_halsvar_batch_e7_v100s/012326_100k_cells_20iter_torch_halsvar_batch_e7_v100s.usages.k_50.dt_2_0.consensus.txt"
cnmf_file = "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50.usages.k_50.dt_0_2.consensus.txt"

# 2. Read the dataframe
usage_df = pd.read_csv(cnmf_file, sep='\t', index_col=0)

# 3. Align with adata
# This ensures every cell in adata gets its corresponding row from the file.
# Cells in adata but missing from the file will be filled with 0 
usage_aligned = usage_df.reindex(adata.obs_names).fillna(0)

# 4. Save to obsm
adata.obsm["cnmf_usage"] = usage_aligned

# 5. Validation
print(f"Loaded usage matrix shape: {adata.obsm['cnmf_usage'].shape}")
print(f"Overlap check: {usage_df.index.intersection(adata.obs_names).shape[0]} cells matched.")
```

### Text output (stream)

```text
Loaded usage matrix shape: (21980, 50)
Overlap check: 21980 cells matched.

```

## Code cell 9

```python
# Guide names must match columns of guide_assignment
guide_names = list(mdata["guide"].var_names)
adata.uns["guide_names"] = guide_names
```

## Code cell 10

```python
# Program names
adata.uns["program_names"] = list(adata.obsm["cnmf_usage"].columns)
```

## Markdown cell 11

#### Import guide assignment

## Code cell 12

```python
guide_file = "/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/d2_grna_assignments_long.csv"

# 1. Read
df_long = pd.read_csv(guide_file)

# 2. Pivot
# aggfunc='max' ensures that if a (cell, guide) pair is listed twice, we just keep 1.
guide_dense = df_long.pivot_table(
    index='cell_id', 
    columns='grna_target', 
    values='value', 
    fill_value=0,
    aggfunc='max' 
)

# 3. Align with adata
# This ensures strict ordering matching your AnnData
guide_aligned = guide_dense.reindex(adata.obs_names, fill_value=0)

# 4. Save
adata.obsm["guide_assignment"] = guide_aligned.astype(np.int8)

print(f"✅ Loaded guide assignment: {adata.obsm['guide_assignment'].shape}")
```

### Text output (stream)

```text
✅ Loaded guide assignment: (21980, 1921)

```

## Markdown cell 13

#### Map guide to gene using guide name/ID

## Code cell 14

```python
# guide2gene: map guide_name -> gene by taking the first token of guide_name.split('_')
guide2gene = {gn: gn.split("_", 1)[0] for gn in guide_names}
adata.uns["guide2gene"] = guide2gene
```

## Code cell 15

```python
adata
```

### Text output (text/plain)

```text
AnnData object with n_obs × n_vars = 21980 × 12732
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0'
    var: 'gene_id', 'symbol', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'n_counts', 'highly_deviant', 'highly_deviant_2k', 'binomial_deviance', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'orig_highly_variable', 'highly_variable_2k'
    uns: 'hvg', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'neighbors', 'pca', 'umap', 'guide_names', 'program_names', 'guide2gene'
    obsm: 'X_pca', 'X_umap', 'guide_assignment', 'guide_counts', 'cnmf_usage'
    varm: 'PCs'
    layers: 'ambiguous', 'analytic_pearson_residuals', 'log1p_norm', 'mature', 'nascent', 'raw_counts', 'scran_normalization'
    obsp: 'connectivities', 'distances'
```

## Markdown cell 16

#### Define covariates

## Code cell 17

```python
# Covariates

#(~ log(response_n_nonzero) + log(response_n_umis) +
#         log(grna_n_nonzero) + log(grna_n_umis) + response_p_mito + batch)

covariates_df = pd.DataFrame({
    #"batch": adata.obs['sample_id'],
    "grna_n_nonzero": np.log1p(adata.obs['guides_per_cell']),
    #"grna_n_umis": np.log1p(adata.obs['guide_umi_counts']),
    #"response_n_nonzero": np.log1p(adata.obs['n_genes_by_counts']),
    "response_n_umis": np.log1p(adata.obs['total_counts']),
    "response_p_mito": adata.obs['pct_counts_mt'],
})
```

## Code cell 18

```python
adata
```

### Text output (text/plain)

```text
AnnData object with n_obs × n_vars = 21980 × 12732
    obs: 'guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0'
    var: 'gene_id', 'symbol', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'n_counts', 'highly_deviant', 'highly_deviant_2k', 'binomial_deviance', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'orig_highly_variable', 'highly_variable_2k'
    uns: 'hvg', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'neighbors', 'pca', 'umap', 'guide_names', 'program_names', 'guide2gene'
    obsm: 'X_pca', 'X_umap', 'guide_assignment', 'guide_counts', 'cnmf_usage'
    varm: 'PCs'
    layers: 'ambiguous', 'analytic_pearson_residuals', 'log1p_norm', 'mature', 'nascent', 'raw_counts', 'scran_normalization'
    obsp: 'connectivities', 'distances'
```

## Code cell 19

```python
adata.obsm["covar"] = covariates_df
```

## Code cell 20

```python
adata.obsm["covar"]
```

### Text output (text/plain)

```text
                            grna_n_nonzero  response_n_umis  response_p_mito
barcode                                                                     
TAGGATGAATTGAGGATGCTTGGG_2        0.693147         7.639754        13.304576
TAGGATGAAGCCATGCCGCCCGGA_2        1.098612         6.745746        17.541105
CCAGTTCACAATGGAATGCTTGGG_3        1.791759         7.545229        15.073812
ATCATTCCAGTCACTACGCCCGGA_2        0.693147         6.311432        16.914215
GCGAGTAAGCTAACGACACAATTG_3        2.397895         7.161039        22.101393
...                                    ...              ...              ...
CCGTGAGAACGCTCGACTTTGGTC_3        1.791759         8.186228        14.275167
CACCTTACACAGCAGACACAATTG_1        2.890372         8.004533        23.300484
CAAGACTAGCCAAGACTAAATATC_3        1.098612         7.842773         9.209585
CAATGGAAGCTCGGTAGTGCTAGC_2        1.945910         7.043451        19.181475
CCAGTTCACTGAGCCATAAATATC_2        2.302585         7.013241        14.636339

[21980 rows x 3 columns]
```

## Markdown cell 21

#### Flooring and Renormalizing cNMF matrix

## Code cell 22

```python
U = adata.obsm["cnmf_usage"]
print("min:", U.min(), "quantile:", np.quantile(U, 1e-4))
print("row sums == 0:", np.sum(U.sum(axis=1) == 0))
```

### Text output (stream)

```text
min: 1     0.0
2     0.0
3     0.0
4     0.0
5     0.0
6     0.0
7     0.0
8     0.0
9     0.0
10    0.0
11    0.0
12    0.0
13    0.0
14    0.0
15    0.0
16    0.0
17    0.0
18    0.0
19    0.0
20    0.0
21    0.0
22    0.0
23    0.0
24    0.0
25    0.0
26    0.0
27    0.0
28    0.0
29    0.0
30    0.0
31    0.0
32    0.0
33    0.0
34    0.0
35    0.0
36    0.0
37    0.0
38    0.0
39    0.0
40    0.0
41    0.0
42    0.0
43    0.0
44    0.0
45    0.0
46    0.0
47    0.0
48    0.0
49    0.0
50    0.0
dtype: float64 quantile: 0.0
row sums == 0: 0

```

## Code cell 23

```python
U0 = adata.obsm["cnmf_usage"]
pos = U0[U0 > 0]
eps = max(1e-8, np.quantile(pos, 0.01) * 0.1)  # e.g. 10% of 1% quantile
Uf = np.where(U0 > 0, U0, eps)
Uf = Uf / Uf.sum(axis=1, keepdims=True)
adata.obsm["cnmf_usage"] = Uf
```

## Markdown cell 24

#### Prepare input for programDE pipeline 

## Code cell 25

```python

inputs = prepare_crt_inputs(
    adata=adata,
    usage_key="cnmf_usage",
    covar_key="covar",
    guide_assignment_key="guide_assignment",
    guide2gene_key="guide2gene"
)
```

## Code cell 26

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

## Code cell 27

```python
import os
cpus_per_task = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
```

## Code cell 28

```python
cpus_per_task
```

### Text output (text/plain)

```text
16
```

## Markdown cell 29

### CRT without stratification

## Code cell 30

```python

out = run_all_genes_union_crt(
    inputs=inputs,
    B=10000,
    n_jobs=cpus_per_task,
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
# Compute raw CRT p-values for each NTC group in each ensemble
ntc_group_pvals_ens = crt_pvals_for_ntc_groups_ensemble(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=10000,
    seed0=23,
)
# Compute skew-calibrated CRT p-values for each NTC group in each ensemble
ntc_group_pvals_skew_ens = crt_pvals_for_ntc_groups_ensemble_skew(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=10000,
    seed0=23,
) 
```

## Code cell 33

```python
# Build CRT-null p-values matched to NTC group units (recommended)
null_pvals = compute_ntc_group_null_pvals_parallel(
    inputs=inputs,
    ntc_groups_ens=ntc_groups_ens,
    B=10000,
    n_jobs=cpus_per_task,
    backend="threading",
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

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/IPython/core/pylabtools.py:170: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  fig.canvas.print_figure(bytes_io, **kw)

```

### Text output (error)

```text
[31m---------------------------------------------------------------------------[39m
[31mKeyboardInterrupt[39m                         Traceback (most recent call last)
[36mCell[39m[36m [39m[32mIn[124][39m[32m, line 15[39m
[32m      1[39m ax = qq_plot_ntc_pvals(
[32m      2[39m     pvals_raw_df=out[[33m"[39m[33mpvals_raw_df[39m[33m"[39m],
[32m      3[39m     guide2gene=adata.uns[[33m"[39m[33mguide2gene[39m[33m"[39m],
[32m   (...)[39m[32m     11[39m     title=[33m"[39m[33mQQ plot: grouped NTC controls (raw vs skew) vs CRT null[39m[33m"[39m,
[32m     12[39m )
[32m     14[39m [38;5;28;01mimport[39;00m[38;5;250m [39m[34;01mmatplotlib[39;00m[34;01m.[39;00m[34;01mpyplot[39;00m[38;5;250m [39m[38;5;28;01mas[39;00m[38;5;250m [39m[34;01mplt[39;00m
[32m---> [39m[32m15[39m [43mplt[49m[43m.[49m[43mshow[49m[43m([49m[43m)[49m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/pyplot.py:613[39m, in [36mshow[39m[34m(*args, **kwargs)[39m
[32m    569[39m [38;5;250m[39m[33;03m"""[39;00m
[32m    570[39m [33;03mDisplay all open figures.[39;00m
[32m    571[39m 
[32m   (...)[39m[32m    610[39m [33;03mexplicitly there.[39;00m
[32m    611[39m [33;03m"""[39;00m
[32m    612[39m _warn_if_gui_out_of_main_thread()
[32m--> [39m[32m613[39m [38;5;28;01mreturn[39;00m [43m_get_backend_mod[49m[43m([49m[43m)[49m[43m.[49m[43mshow[49m[43m([49m[43m*[49m[43margs[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkwargs[49m[43m)[49m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib_inline/backend_inline.py:90[39m, in [36mshow[39m[34m(close, block)[39m
[32m     88[39m [38;5;28;01mtry[39;00m:
[32m     89[39m     [38;5;28;01mfor[39;00m figure_manager [38;5;129;01min[39;00m Gcf.get_all_fig_managers():
[32m---> [39m[32m90[39m         [43mdisplay[49m[43m([49m
[32m     91[39m [43m            [49m[43mfigure_manager[49m[43m.[49m[43mcanvas[49m[43m.[49m[43mfigure[49m[43m,[49m
[32m     92[39m [43m            [49m[43mmetadata[49m[43m=[49m[43m_fetch_figure_metadata[49m[43m([49m[43mfigure_manager[49m[43m.[49m[43mcanvas[49m[43m.[49m[43mfigure[49m[43m)[49m[43m,[49m
[32m     93[39m [43m        [49m[43m)[49m
[32m     94[39m [38;5;28;01mfinally[39;00m:
[32m     95[39m     show._to_draw = []

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/IPython/core/display_functions.py:278[39m, in [36mdisplay[39m[34m(include, exclude, metadata, transient, display_id, raw, clear, *objs, **kwargs)[39m
[32m    276[39m     publish_display_data(data=obj, metadata=metadata, **kwargs)
[32m    277[39m [38;5;28;01melse[39;00m:
[32m--> [39m[32m278[39m     format_dict, md_dict = [38;5;28;43mformat[39;49m[43m([49m[43mobj[49m[43m,[49m[43m [49m[43minclude[49m[43m=[49m[43minclude[49m[43m,[49m[43m [49m[43mexclude[49m[43m=[49m[43mexclude[49m[43m)[49m
[32m    279[39m     [38;5;28;01mif[39;00m [38;5;129;01mnot[39;00m format_dict:
[32m    280[39m         [38;5;66;03m# nothing to display (e.g. _ipython_display_ took over)[39;00m
[32m    281[39m         [38;5;28;01mcontinue[39;00m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/IPython/core/formatters.py:238[39m, in [36mDisplayFormatter.format[39m[34m(self, obj, include, exclude)[39m
[32m    236[39m md = [38;5;28;01mNone[39;00m
[32m    237[39m [38;5;28;01mtry[39;00m:
[32m--> [39m[32m238[39m     data = [43mformatter[49m[43m([49m[43mobj[49m[43m)[49m
[32m    239[39m [38;5;28;01mexcept[39;00m:
[32m    240[39m     [38;5;66;03m# FIXME: log the exception[39;00m
[32m    241[39m     [38;5;28;01mraise[39;00m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/decorator.py:235[39m, in [36mdecorate.<locals>.fun[39m[34m(*args, **kw)[39m
[32m    233[39m [38;5;28;01mif[39;00m [38;5;129;01mnot[39;00m kwsyntax:
[32m    234[39m     args, kw = fix(args, kw, sig)
[32m--> [39m[32m235[39m [38;5;28;01mreturn[39;00m [43mcaller[49m[43m([49m[43mfunc[49m[43m,[49m[43m [49m[43m*[49m[43m([49m[43mextras[49m[43m [49m[43m+[49m[43m [49m[43margs[49m[43m)[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkw[49m[43m)[49m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/IPython/core/formatters.py:282[39m, in [36mcatch_format_error[39m[34m(method, self, *args, **kwargs)[39m
[32m    280[39m [38;5;250m[39m[33;03m"""show traceback on failed format call"""[39;00m
[32m    281[39m [38;5;28;01mtry[39;00m:
[32m--> [39m[32m282[39m     r = [43mmethod[49m[43m([49m[38;5;28;43mself[39;49m[43m,[49m[43m [49m[43m*[49m[43margs[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkwargs[49m[43m)[49m
[32m    283[39m [38;5;28;01mexcept[39;00m [38;5;167;01mNotImplementedError[39;00m:
[32m    284[39m     [38;5;66;03m# don't warn on NotImplementedErrors[39;00m
[32m    285[39m     [38;5;28;01mreturn[39;00m [38;5;28mself[39m._check_return([38;5;28;01mNone[39;00m, args[[32m0[39m])

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/IPython/core/formatters.py:402[39m, in [36mBaseFormatter.__call__[39m[34m(self, obj)[39m
[32m    400[39m     [38;5;28;01mpass[39;00m
[32m    401[39m [38;5;28;01melse[39;00m:
[32m--> [39m[32m402[39m     [38;5;28;01mreturn[39;00m [43mprinter[49m[43m([49m[43mobj[49m[43m)[49m
[32m    403[39m [38;5;66;03m# Finally look for special method names[39;00m
[32m    404[39m method = get_real_method(obj, [38;5;28mself[39m.print_method)

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/IPython/core/pylabtools.py:170[39m, in [36mprint_figure[39m[34m(fig, fmt, bbox_inches, base64, **kwargs)[39m
[32m    167[39m     [38;5;28;01mfrom[39;00m[38;5;250m [39m[34;01mmatplotlib[39;00m[34;01m.[39;00m[34;01mbackend_bases[39;00m[38;5;250m [39m[38;5;28;01mimport[39;00m FigureCanvasBase
[32m    168[39m     FigureCanvasBase(fig)
[32m--> [39m[32m170[39m [43mfig[49m[43m.[49m[43mcanvas[49m[43m.[49m[43mprint_figure[49m[43m([49m[43mbytes_io[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkw[49m[43m)[49m
[32m    171[39m data = bytes_io.getvalue()
[32m    172[39m [38;5;28;01mif[39;00m fmt == [33m'[39m[33msvg[39m[33m'[39m:

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/backend_bases.py:2186[39m, in [36mFigureCanvasBase.print_figure[39m[34m(self, filename, dpi, facecolor, edgecolor, orientation, format, bbox_inches, pad_inches, bbox_extra_artists, backend, **kwargs)[39m
[32m   2182[39m [38;5;28;01mtry[39;00m:
[32m   2183[39m     [38;5;66;03m# _get_renderer may change the figure dpi (as vector formats[39;00m
[32m   2184[39m     [38;5;66;03m# force the figure dpi to 72), so we need to set it again here.[39;00m
[32m   2185[39m     [38;5;28;01mwith[39;00m cbook._setattr_cm([38;5;28mself[39m.figure, dpi=dpi):
[32m-> [39m[32m2186[39m         result = [43mprint_method[49m[43m([49m
[32m   2187[39m [43m            [49m[43mfilename[49m[43m,[49m
[32m   2188[39m [43m            [49m[43mfacecolor[49m[43m=[49m[43mfacecolor[49m[43m,[49m
[32m   2189[39m [43m            [49m[43medgecolor[49m[43m=[49m[43medgecolor[49m[43m,[49m
[32m   2190[39m [43m            [49m[43morientation[49m[43m=[49m[43morientation[49m[43m,[49m
[32m   2191[39m [43m            [49m[43mbbox_inches_restore[49m[43m=[49m[43m_bbox_inches_restore[49m[43m,[49m
[32m   2192[39m [43m            [49m[43m*[49m[43m*[49m[43mkwargs[49m[43m)[49m
[32m   2193[39m [38;5;28;01mfinally[39;00m:
[32m   2194[39m     [38;5;28;01mif[39;00m bbox_inches [38;5;129;01mand[39;00m restore_bbox:

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/backend_bases.py:2042[39m, in [36mFigureCanvasBase._switch_canvas_and_return_print_method.<locals>.<lambda>[39m[34m(*args, **kwargs)[39m
[32m   2038[39m     optional_kws = {  [38;5;66;03m# Passed by print_figure for other renderers.[39;00m
[32m   2039[39m         [33m"[39m[33mdpi[39m[33m"[39m, [33m"[39m[33mfacecolor[39m[33m"[39m, [33m"[39m[33medgecolor[39m[33m"[39m, [33m"[39m[33morientation[39m[33m"[39m,
[32m   2040[39m         [33m"[39m[33mbbox_inches_restore[39m[33m"[39m}
[32m   2041[39m     skip = optional_kws - {*inspect.signature(meth).parameters}
[32m-> [39m[32m2042[39m     print_method = functools.wraps(meth)([38;5;28;01mlambda[39;00m *args, **kwargs: [43mmeth[49m[43m([49m
[32m   2043[39m [43m        [49m[43m*[49m[43margs[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43m{[49m[43mk[49m[43m:[49m[43m [49m[43mv[49m[43m [49m[38;5;28;43;01mfor[39;49;00m[43m [49m[43mk[49m[43m,[49m[43m [49m[43mv[49m[43m [49m[38;5;129;43;01min[39;49;00m[43m [49m[43mkwargs[49m[43m.[49m[43mitems[49m[43m([49m[43m)[49m[43m [49m[38;5;28;43;01mif[39;49;00m[43m [49m[43mk[49m[43m [49m[38;5;129;43;01mnot[39;49;00m[43m [49m[38;5;129;43;01min[39;49;00m[43m [49m[43mskip[49m[43m}[49m[43m)[49m)
[32m   2044[39m [38;5;28;01melse[39;00m:  [38;5;66;03m# Let third-parties do as they see fit.[39;00m
[32m   2045[39m     print_method = meth

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/backends/backend_agg.py:481[39m, in [36mFigureCanvasAgg.print_png[39m[34m(self, filename_or_obj, metadata, pil_kwargs)[39m
[32m    434[39m [38;5;28;01mdef[39;00m[38;5;250m [39m[34mprint_png[39m([38;5;28mself[39m, filename_or_obj, *, metadata=[38;5;28;01mNone[39;00m, pil_kwargs=[38;5;28;01mNone[39;00m):
[32m    435[39m [38;5;250m    [39m[33;03m"""[39;00m
[32m    436[39m [33;03m    Write the figure to a PNG file.[39;00m
[32m    437[39m 
[32m   (...)[39m[32m    479[39m [33;03m        *metadata*, including the default 'Software' key.[39;00m
[32m    480[39m [33;03m    """[39;00m
[32m--> [39m[32m481[39m     [38;5;28;43mself[39;49m[43m.[49m[43m_print_pil[49m[43m([49m[43mfilename_or_obj[49m[43m,[49m[43m [49m[33;43m"[39;49m[33;43mpng[39;49m[33;43m"[39;49m[43m,[49m[43m [49m[43mpil_kwargs[49m[43m,[49m[43m [49m[43mmetadata[49m[43m)[49m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/backends/backend_agg.py:429[39m, in [36mFigureCanvasAgg._print_pil[39m[34m(self, filename_or_obj, fmt, pil_kwargs, metadata)[39m
[32m    424[39m [38;5;28;01mdef[39;00m[38;5;250m [39m[34m_print_pil[39m([38;5;28mself[39m, filename_or_obj, fmt, pil_kwargs, metadata=[38;5;28;01mNone[39;00m):
[32m    425[39m [38;5;250m    [39m[33;03m"""[39;00m
[32m    426[39m [33;03m    Draw the canvas, then save it using `.image.imsave` (to which[39;00m
[32m    427[39m [33;03m    *pil_kwargs* and *metadata* are forwarded).[39;00m
[32m    428[39m [33;03m    """[39;00m
[32m--> [39m[32m429[39m     [43mFigureCanvasAgg[49m[43m.[49m[43mdraw[49m[43m([49m[38;5;28;43mself[39;49m[43m)[49m
[32m    430[39m     mpl.image.imsave(
[32m    431[39m         filename_or_obj, [38;5;28mself[39m.buffer_rgba(), [38;5;28mformat[39m=fmt, origin=[33m"[39m[33mupper[39m[33m"[39m,
[32m    432[39m         dpi=[38;5;28mself[39m.figure.dpi, metadata=metadata, pil_kwargs=pil_kwargs)

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/backends/backend_agg.py:382[39m, in [36mFigureCanvasAgg.draw[39m[34m(self)[39m
[32m    379[39m [38;5;66;03m# Acquire a lock on the shared font cache.[39;00m
[32m    380[39m [38;5;28;01mwith[39;00m ([38;5;28mself[39m.toolbar._wait_cursor_for_draw_cm() [38;5;28;01mif[39;00m [38;5;28mself[39m.toolbar
[32m    381[39m       [38;5;28;01melse[39;00m nullcontext()):
[32m--> [39m[32m382[39m     [38;5;28;43mself[39;49m[43m.[49m[43mfigure[49m[43m.[49m[43mdraw[49m[43m([49m[38;5;28;43mself[39;49m[43m.[49m[43mrenderer[49m[43m)[49m
[32m    383[39m     [38;5;66;03m# A GUI class may be need to update a window using this draw, so[39;00m
[32m    384[39m     [38;5;66;03m# don't forget to call the superclass.[39;00m
[32m    385[39m     [38;5;28msuper[39m().draw()

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/artist.py:94[39m, in [36m_finalize_rasterization.<locals>.draw_wrapper[39m[34m(artist, renderer, *args, **kwargs)[39m
[32m     92[39m [38;5;129m@wraps[39m(draw)
[32m     93[39m [38;5;28;01mdef[39;00m[38;5;250m [39m[34mdraw_wrapper[39m(artist, renderer, *args, **kwargs):
[32m---> [39m[32m94[39m     result = [43mdraw[49m[43m([49m[43martist[49m[43m,[49m[43m [49m[43mrenderer[49m[43m,[49m[43m [49m[43m*[49m[43margs[49m[43m,[49m[43m [49m[43m*[49m[43m*[49m[43mkwargs[49m[43m)[49m
[32m     95[39m     [38;5;28;01mif[39;00m renderer._rasterizing:
[32m     96[39m         renderer.stop_rasterizing()

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/artist.py:71[39m, in [36mallow_rasterization.<locals>.draw_wrapper[39m[34m(artist, renderer)[39m
[32m     68[39m     [38;5;28;01mif[39;00m artist.get_agg_filter() [38;5;129;01mis[39;00m [38;5;129;01mnot[39;00m [38;5;28;01mNone[39;00m:
[32m     69[39m         renderer.start_filter()
[32m---> [39m[32m71[39m     [38;5;28;01mreturn[39;00m [43mdraw[49m[43m([49m[43martist[49m[43m,[49m[43m [49m[43mrenderer[49m[43m)[49m
[32m     72[39m [38;5;28;01mfinally[39;00m:
[32m     73[39m     [38;5;28;01mif[39;00m artist.get_agg_filter() [38;5;129;01mis[39;00m [38;5;129;01mnot[39;00m [38;5;28;01mNone[39;00m:

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/figure.py:3257[39m, in [36mFigure.draw[39m[34m(self, renderer)[39m
[32m   3254[39m             [38;5;66;03m# ValueError can occur when resizing a window.[39;00m
[32m   3256[39m     [38;5;28mself[39m.patch.draw(renderer)
[32m-> [39m[32m3257[39m     [43mmimage[49m[43m.[49m[43m_draw_list_compositing_images[49m[43m([49m
[32m   3258[39m [43m        [49m[43mrenderer[49m[43m,[49m[43m [49m[38;5;28;43mself[39;49m[43m,[49m[43m [49m[43martists[49m[43m,[49m[43m [49m[38;5;28;43mself[39;49m[43m.[49m[43msuppressComposite[49m[43m)[49m
[32m   3260[39m     renderer.close_group([33m'[39m[33mfigure[39m[33m'[39m)
[32m   3261[39m [38;5;28;01mfinally[39;00m:

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/image.py:134[39m, in [36m_draw_list_compositing_images[39m[34m(renderer, parent, artists, suppress_composite)[39m
[32m    132[39m [38;5;28;01mif[39;00m not_composite [38;5;129;01mor[39;00m [38;5;129;01mnot[39;00m has_images:
[32m    133[39m     [38;5;28;01mfor[39;00m a [38;5;129;01min[39;00m artists:
[32m--> [39m[32m134[39m         [43ma[49m[43m.[49m[43mdraw[49m[43m([49m[43mrenderer[49m[43m)[49m
[32m    135[39m [38;5;28;01melse[39;00m:
[32m    136[39m     [38;5;66;03m# Composite any adjacent images together[39;00m
[32m    137[39m     image_group = []

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/artist.py:71[39m, in [36mallow_rasterization.<locals>.draw_wrapper[39m[34m(artist, renderer)[39m
[32m     68[39m     [38;5;28;01mif[39;00m artist.get_agg_filter() [38;5;129;01mis[39;00m [38;5;129;01mnot[39;00m [38;5;28;01mNone[39;00m:
[32m     69[39m         renderer.start_filter()
[32m---> [39m[32m71[39m     [38;5;28;01mreturn[39;00m [43mdraw[49m[43m([49m[43martist[49m[43m,[49m[43m [49m[43mrenderer[49m[43m)[49m
[32m     72[39m [38;5;28;01mfinally[39;00m:
[32m     73[39m     [38;5;28;01mif[39;00m artist.get_agg_filter() [38;5;129;01mis[39;00m [38;5;129;01mnot[39;00m [38;5;28;01mNone[39;00m:

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/axes/_base.py:3226[39m, in [36m_AxesBase.draw[39m[34m(self, renderer)[39m
[32m   3223[39m [38;5;28;01mif[39;00m artists_rasterized:
[32m   3224[39m     _draw_rasterized([38;5;28mself[39m.get_figure(root=[38;5;28;01mTrue[39;00m), artists_rasterized, renderer)
[32m-> [39m[32m3226[39m [43mmimage[49m[43m.[49m[43m_draw_list_compositing_images[49m[43m([49m
[32m   3227[39m [43m    [49m[43mrenderer[49m[43m,[49m[43m [49m[38;5;28;43mself[39;49m[43m,[49m[43m [49m[43martists[49m[43m,[49m[43m [49m[38;5;28;43mself[39;49m[43m.[49m[43mget_figure[49m[43m([49m[43mroot[49m[43m=[49m[38;5;28;43;01mTrue[39;49;00m[43m)[49m[43m.[49m[43msuppressComposite[49m[43m)[49m
[32m   3229[39m renderer.close_group([33m'[39m[33maxes[39m[33m'[39m)
[32m   3230[39m [38;5;28mself[39m.stale = [38;5;28;01mFalse[39;00m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/image.py:134[39m, in [36m_draw_list_compositing_images[39m[34m(renderer, parent, artists, suppress_composite)[39m
[32m    132[39m [38;5;28;01mif[39;00m not_composite [38;5;129;01mor[39;00m [38;5;129;01mnot[39;00m has_images:
[32m    133[39m     [38;5;28;01mfor[39;00m a [38;5;129;01min[39;00m artists:
[32m--> [39m[32m134[39m         [43ma[49m[43m.[49m[43mdraw[49m[43m([49m[43mrenderer[49m[43m)[49m
[32m    135[39m [38;5;28;01melse[39;00m:
[32m    136[39m     [38;5;66;03m# Composite any adjacent images together[39;00m
[32m    137[39m     image_group = []

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/artist.py:71[39m, in [36mallow_rasterization.<locals>.draw_wrapper[39m[34m(artist, renderer)[39m
[32m     68[39m     [38;5;28;01mif[39;00m artist.get_agg_filter() [38;5;129;01mis[39;00m [38;5;129;01mnot[39;00m [38;5;28;01mNone[39;00m:
[32m     69[39m         renderer.start_filter()
[32m---> [39m[32m71[39m     [38;5;28;01mreturn[39;00m [43mdraw[49m[43m([49m[43martist[49m[43m,[49m[43m [49m[43mrenderer[49m[43m)[49m
[32m     72[39m [38;5;28;01mfinally[39;00m:
[32m     73[39m     [38;5;28;01mif[39;00m artist.get_agg_filter() [38;5;129;01mis[39;00m [38;5;129;01mnot[39;00m [38;5;28;01mNone[39;00m:

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/collections.py:1008[39m, in [36m_CollectionWithSizes.draw[39m[34m(self, renderer)[39m
[32m   1005[39m [38;5;129m@artist[39m.allow_rasterization
[32m   1006[39m [38;5;28;01mdef[39;00m[38;5;250m [39m[34mdraw[39m([38;5;28mself[39m, renderer):
[32m   1007[39m     [38;5;28mself[39m.set_sizes([38;5;28mself[39m._sizes, [38;5;28mself[39m.get_figure(root=[38;5;28;01mTrue[39;00m).dpi)
[32m-> [39m[32m1008[39m     [38;5;28;43msuper[39;49m[43m([49m[43m)[49m[43m.[49m[43mdraw[49m[43m([49m[43mrenderer[49m[43m)[49m

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/artist.py:71[39m, in [36mallow_rasterization.<locals>.draw_wrapper[39m[34m(artist, renderer)[39m
[32m     68[39m     [38;5;28;01mif[39;00m artist.get_agg_filter() [38;5;129;01mis[39;00m [38;5;129;01mnot[39;00m [38;5;28;01mNone[39;00m:
[32m     69[39m         renderer.start_filter()
[32m---> [39m[32m71[39m     [38;5;28;01mreturn[39;00m [43mdraw[49m[43m([49m[43martist[49m[43m,[49m[43m [49m[43mrenderer[49m[43m)[49m
[32m     72[39m [38;5;28;01mfinally[39;00m:
[32m     73[39m     [38;5;28;01mif[39;00m artist.get_agg_filter() [38;5;129;01mis[39;00m [38;5;129;01mnot[39;00m [38;5;28;01mNone[39;00m:

[36mFile [39m[32m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/matplotlib/collections.py:415[39m, in [36mCollection.draw[39m[34m(self, renderer)[39m
[32m    413[39m     gc.set_antialiased([38;5;28mself[39m._antialiaseds[[32m0[39m])
[32m    414[39m     gc.set_url([38;5;28mself[39m._urls[[32m0[39m])
[32m--> [39m[32m415[39m     [43mrenderer[49m[43m.[49m[43mdraw_markers[49m[43m([49m
[32m    416[39m [43m        [49m[43mgc[49m[43m,[49m[43m [49m[43mpaths[49m[43m[[49m[32;43m0[39;49m[43m][49m[43m,[49m[43m [49m[43mcombined_transform[49m[43m.[49m[43mfrozen[49m[43m([49m[43m)[49m[43m,[49m
[32m    417[39m [43m        [49m[43mmpath[49m[43m.[49m[43mPath[49m[43m([49m[43moffsets[49m[43m)[49m[43m,[49m[43m [49m[43moffset_trf[49m[43m,[49m[43m [49m[38;5;28;43mtuple[39;49m[43m([49m[43mfacecolors[49m[43m[[49m[32;43m0[39;49m[43m][49m[43m)[49m[43m)[49m
[32m    418[39m [38;5;28;01melse[39;00m:
[32m    419[39m     [38;5;28;01mif[39;00m [38;5;28mself[39m._gapcolor [38;5;129;01mis[39;00m [38;5;129;01mnot[39;00m [38;5;28;01mNone[39;00m:
[32m    420[39m         [38;5;66;03m# First draw paths within the gaps.[39;00m

[31mKeyboardInterrupt[39m: 
```

## Code cell 35

```python
out["pvals_raw_df"]
```

### Text output (text/plain)

```text
                      1         2         3         4         5         6  \
ACAA1          0.135486  0.442356  0.048595  0.427557  0.723928  0.606439   
ACTL7B         0.988001  0.691931  0.048695  0.280272  0.225777  0.535646   
ACVR1B         0.622738  0.686131  0.254375  0.040996  0.329067  0.549545   
ACVR2A         0.013399  0.001500  0.126387  0.019998  0.121988  0.690531   
ACVR2B         0.096690  0.000100  0.858214  0.184382  0.055794  0.394561   
...                 ...       ...       ...       ...       ...       ...   
ZBTB22         0.001000  0.698930  0.493751  0.746125  0.186781  0.086991   
ZBTB48         0.087391  0.186881  0.762824  0.601240  0.129587  0.246475   
ZC3H8          0.710629  0.010499  0.012499  0.000500  0.513549  0.660434   
ZNF281         0.054195  0.834817  0.000100  0.674733  0.066193  0.334267   
non-targeting  0.671233  0.188581  0.351365  0.024098  0.040596  0.370663   

                      7         8         9        10  ...        41  \
ACAA1          0.280672  0.590041  0.548145  0.521948  ...  0.236276   
ACTL7B         0.656934  0.936706  0.495950  0.956604  ...  0.280572   
ACVR1B         0.293771  0.955704  0.139886  0.373963  ...  0.070893   
ACVR2A         0.000100  0.012099  0.375062  0.281772  ...  0.229777   
ACVR2B         0.000300  0.152685  0.102690  0.161584  ...  0.181382   
...                 ...       ...       ...       ...  ...       ...   
ZBTB22         0.193181  0.403260  0.124588  0.558844  ...  0.768923   
ZBTB48         0.141586  0.346465  0.484352  0.491451  ...  0.578142   
ZC3H8          0.141586  0.210579  0.097890  0.292571  ...  0.035396   
ZNF281         0.063994  0.254975  0.776522  0.359864  ...  0.747925   
non-targeting  0.250675  0.297570  0.529447  0.266573  ...  0.326767   

                     42        43        44        45        46        47  \
ACAA1          0.337666  0.096590  0.438756  0.196480  0.412259  0.559744   
ACTL7B         0.720228  0.414459  0.532647  0.849015  0.893911  0.027597   
ACVR1B         0.867713  0.307769  0.410359  0.106689  0.073093  0.727727   
ACVR2A         0.150185  0.692331  0.011099  0.956604  0.938006  0.706929   
ACVR2B         0.918408  0.181982  0.094491  0.279572  0.193881  0.200280   
...                 ...       ...       ...       ...       ...       ...   
ZBTB22         0.060194  0.047895  0.283772  0.654035  0.746525  0.384862   
ZBTB48         0.425457  0.917208  0.011499  0.566143  0.789621  0.696330   
ZC3H8          0.436956  0.270573  0.400060  0.911409  0.298170  0.732527   
ZNF281         0.329367  0.089591  0.620238  0.415558  0.782922  0.147885   
non-targeting  0.812219  0.727927  0.014999  0.926207  0.319568  0.442356   

                     48        49        50  
ACAA1          0.955004  0.838916  0.128587  
ACTL7B         0.550945  0.566643  0.930707  
ACVR1B         0.798520  0.187081  0.362964  
ACVR2A         0.171683  0.072793  0.879812  
ACVR2B         0.973803  0.409959  0.641036  
...                 ...       ...       ...  
ZBTB22         0.389761  0.049895  0.711529  
ZBTB48         0.520348  0.987701  0.199980  
ZC3H8          0.254975  0.998000  0.965603  
ZNF281         0.087991  0.014199  0.493151  
non-targeting  0.323268  0.603440  0.350165  

[300 rows x 50 columns]
```

## Code cell 36

```python
out["betas_df"]
```

### Text output (text/plain)

```text
                      1         2         3         4         5         6  \
ACAA1         -0.328642  0.238046  0.629823  0.257959 -0.130893  0.157828   
ACTL7B         0.002802  0.119574  0.613442  0.344181 -0.424045 -0.186709   
ACVR1B         0.106183  0.130888  0.373346 -0.659666 -0.353413  0.190764   
ACVR2A        -0.527418  0.863553 -0.498809 -0.745594 -0.530854 -0.119915   
ACVR2B        -0.308922  0.975739 -0.052427 -0.376760 -0.601753 -0.228451   
...                 ...       ...       ...       ...       ...       ...   
ZBTB22        -0.557315 -0.092820 -0.165844  0.078092 -0.358116 -0.404216   
ZBTB48        -0.365417  0.379368  0.088641 -0.165702 -0.537639 -0.348253   
ZC3H8          0.091746  0.868299 -0.978977 -1.224034  0.274389  0.156210   
ZNF281        -0.424846 -0.064235  1.712410  0.142043 -0.693301 -0.307131   
non-targeting -0.081746  0.250896  0.154683  0.323197 -0.389589 -0.152470   

                      7         8         9        10  ...        41  \
ACAA1          0.323831 -0.169339  0.207468  0.167864  ...  0.375309   
ACTL7B        -0.130359  0.026763  0.234202  0.014184  ... -0.338212   
ACVR1B         0.319974 -0.015995 -0.525219 -0.230868  ...  0.587029   
ACVR2A         1.614839  0.766469 -0.312946 -0.264218  ...  0.380392   
ACVR2B         1.015727  0.410078 -0.517485 -0.320970  ...  0.383880   
...                 ...       ...       ...       ...  ...       ...   
ZBTB22         0.297925  0.198549  0.405326 -0.112814  ...  0.074315   
ZBTB48         0.427763  0.284276 -0.238738 -0.177336  ... -0.170744   
ZC3H8         -0.516155  0.455934  0.674172 -0.300047  ...  0.784291   
ZNF281         0.572575 -0.376259  0.100522  0.246716  ...  0.107233   
non-targeting  0.171549  0.127294  0.086817 -0.135151  ...  0.140337   

                     42        43        44        45        46        47  \
ACAA1          0.229493  0.318265 -0.248040 -0.307044 -0.184281 -0.189916   
ACTL7B         0.085054  0.154688 -0.194134 -0.042407 -0.029233  0.716215   
ACVR1B        -0.038618  0.200475  0.266609  0.390980  0.408582  0.117151   
ACVR2A        -0.351463 -0.077070  0.785341 -0.012772  0.017561 -0.119625   
ACVR2B         0.021768  0.238408  0.478524 -0.231797 -0.263525  0.374646   
...                 ...       ...       ...       ...       ...       ...   
ZBTB22        -0.338259 -0.288896  0.259437 -0.079194 -0.055856  0.216671   
ZBTB48         0.187074 -0.020081  0.801457 -0.133964 -0.059611  0.123193   
ZC3H8          0.216634  0.247687 -0.317015  0.031268 -0.271691 -0.127339   
ZNF281         0.246648  0.339471  0.165096  0.198687  0.066602 -0.487266   
non-targeting -0.022234 -0.027086  0.345972 -0.009136  0.085814  0.100767   

                     48        49        50  
ACAA1         -0.015181  0.043803 -0.259049  
ACTL7B        -0.180391  0.124818  0.014116  
ACVR1B         0.077108 -0.300462 -0.157441  
ACVR2A         0.420215  0.399494  0.026852  
ACVR2B         0.009158  0.167262  0.072806  
...                 ...       ...       ...  
ZBTB22         0.202423  0.331668 -0.048075  
ZBTB48        -0.186538  0.003032  0.211466  
ZC3H8         -0.417202  0.000665  0.009588  
ZNF281         0.559339 -0.569393  0.120430  
non-targeting -0.131727 -0.047372  0.064344  

[300 rows x 50 columns]
```

## Code cell 37

```python
#d2_save_dir = '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/'
#out["pvals_raw_df"].to_csv(d2_save_dir + 'd2_pvals.csv')
#out["betas_df"].to_csv(d2_save_dir + 'd2_betas.csv')
```

## Code cell 38

```python
d2_save_dir = '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/d2_sceptre_transfer/'
out["pvals_raw_df"].to_csv(d2_save_dir + 'd2_pvals.csv')
out["betas_df"].to_csv(d2_save_dir + 'd2_betas.csv')
```

