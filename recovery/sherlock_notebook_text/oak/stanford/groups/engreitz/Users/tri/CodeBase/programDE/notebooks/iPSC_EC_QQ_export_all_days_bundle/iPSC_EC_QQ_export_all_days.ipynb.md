# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/iPSC_EC_QQ_export_all_days_bundle/iPSC_EC_QQ_export_all_days.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/iPSC_EC_QQ_export_all_days_bundle/iPSC_EC_QQ_export_all_days.ipynb`
- SHA-256: `6df06da59640bb0ac71e1ffb92849b9a588ef9ae63fc306bd0bdd00944bb9372`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# iPSC-EC CRT QQ point export for d0–d3

This notebook wrapper runs `iPSC_EC_QQ_export_all_days.py`. It exports per-day observed-target and non-targeting pseudo-target p-values plus QQ coordinates for replotting. The paths and parameter defaults are taken from `qq_sweep_test.ipynb` and `ipsc_ec_qq.ipynb`.

## Code cell 1

```python
from pathlib import Path

PROGRAMDE_ROOT = '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE'
DATA_DIR = '/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data'
CNMF_DIR = '/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all'
GUIDE_LONG_TEMPLATE = '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/{day}_sceptre_transfer/{day}_grna_assignments_long.csv'
OUT_DIR = 'results/ipsc_ec_qq_exports'

DAYS = ['d0', 'd1', 'd2', 'd3']
K = 50
DENSITY_THRESHOLD = '0_2'
B = 10000
N_JOBS = 4

print('Output root:', Path(OUT_DIR) / f'k_{K}.dt_{DENSITY_THRESHOLD}')
```

### Text output (stream)

```text
Output root: results/ipsc_ec_qq_exports/k_50.dt_0_2

```

## Markdown cell 2

## Optional path check

This cell does not run CRT. It prints the H5MU, guide, and cNMF usage candidates that the script will use.

## Code cell 3

```python
!python iPSC_EC_QQ_export_all_days.py \
  --programde-root "{PROGRAMDE_ROOT}" \
  --data-dir "{DATA_DIR}" \
  --cnmf-dir "{CNMF_DIR}" \
  --guide-long-csv-template "{GUIDE_LONG_TEMPLATE}" \
  --out-dir "{OUT_DIR}" \
  --days {' '.join(DAYS)} \
  --k {K} \
  --density-threshold "{DENSITY_THRESHOLD}" \
  --B {B} \
  --n-jobs {N_JOBS} \
  --ntc-label SAFE \
  --ntc-label non-targeting \
  --ntc-label NTC \
  --dry-run
```

### Text output (stream)

```text
iPSC-EC QQ export
days: ['d0', 'd1', 'd2', 'd3']
model: k_50.dt_0_2
data_dir: /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data
cnmf_dir: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all
out_dir: results/ipsc_ec_qq_exports/k_50.dt_0_2
B: 10000 n_jobs: 4

Combined outputs:
  manifest: results/ipsc_ec_qq_exports/k_50.dt_0_2/QQ_export_manifest.csv

```

## Markdown cell 4

## Run export

This performs the CRT runs and writes per-day CSVs. For a quicker smoke test, set `B = 1023` before running this cell.

## Code cell 5

```python
!python iPSC_EC_QQ_export_all_days.py \
  --programde-root "{PROGRAMDE_ROOT}" \
  --data-dir "{DATA_DIR}" \
  --cnmf-dir "{CNMF_DIR}" \
  --guide-long-csv-template "{GUIDE_LONG_TEMPLATE}" \
  --out-dir "{OUT_DIR}" \
  --days {' '.join(DAYS)} \
  --k {K} \
  --density-threshold "{DENSITY_THRESHOLD}" \
  --B {B} \
  --n-jobs {N_JOBS} \
  --ntc-label SAFE \
  --ntc-label non-targeting \
  --ntc-label NTC \
  --skip-if-exists
```

### Text output (stream)

```text
iPSC-EC QQ export
days: ['d0', 'd1', 'd2', 'd3']
model: k_50.dt_0_2
data_dir: /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data
cnmf_dir: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all
out_dir: results/ipsc_ec_qq_exports/k_50.dt_0_2
B: 10000 n_jobs: 4
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/muon/_core/preproc.py:31: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
  if Version(scanpy.__version__) < Version("1.10"):

=== d0: K=50, dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage candidates scored for d0
  overlap=38,609  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt
selected usage file for d0: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/iPSC_EC_QQ_export_all_days_bundle/iPSC_EC_QQ_export_all_days.py:301: UserWarning: Guide metadata file 300genes_guide_metadata_v43.tsv was not found. Falling back to guide_name.split('_', 1)[0].
  warnings.warn(
cells=38,609; programs=50; guides=1,926; usage_matched=38,609; covariates=['grna_n_nonzero', 'response_n_umis', 'response_p_mito']
[0m[0m[0m[0msaved d0: observed target empirical points=14,900; NTC empirical points=6,100; output=results/ipsc_ec_qq_exports/k_50.dt_0_2/d0/d0_QQ_points.csv

=== d1: K=50, dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage candidates scored for d1
  overlap=16,925  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt
selected usage file for d1: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/iPSC_EC_QQ_export_all_days_bundle/iPSC_EC_QQ_export_all_days.py:301: UserWarning: Guide metadata file 300genes_guide_metadata_v43.tsv was not found. Falling back to guide_name.split('_', 1)[0].
  warnings.warn(
cells=16,925; programs=50; guides=1,926; usage_matched=16,925; covariates=['grna_n_nonzero', 'response_n_umis', 'response_p_mito']
[0m[0m[0m[0msaved d1: observed target empirical points=14,900; NTC empirical points=4,650; output=results/ipsc_ec_qq_exports/k_50.dt_0_2/d1/d1_QQ_points.csv

=== d2: K=50, dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage candidates scored for d2
  overlap=21,980  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt
selected usage file for d2: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/iPSC_EC_QQ_export_all_days_bundle/iPSC_EC_QQ_export_all_days.py:301: UserWarning: Guide metadata file 300genes_guide_metadata_v43.tsv was not found. Falling back to guide_name.split('_', 1)[0].
  warnings.warn(
cells=21,980; programs=50; guides=1,926; usage_matched=21,980; covariates=['grna_n_nonzero', 'response_n_umis', 'response_p_mito']
[0m[0m[0m[0msaved d2: observed target empirical points=14,900; NTC empirical points=6,350; output=results/ipsc_ec_qq_exports/k_50.dt_0_2/d2/d2_QQ_points.csv

=== d3: K=50, dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage candidates scored for d3
  overlap=44,102  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.starcat_spectra.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.spectra.k_50.dt_0_2.consensus.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_tpm.k_50.dt_0_2.txt
  overlap=0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.gene_spectra_score.k_50.dt_0_2.txt
selected usage file for d3: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all.usages.k_50.dt_0_2.consensus.txt
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/iPSC_EC_QQ_export_all_days_bundle/iPSC_EC_QQ_export_all_days.py:301: UserWarning: Guide metadata file 300genes_guide_metadata_v43.tsv was not found. Falling back to guide_name.split('_', 1)[0].
  warnings.warn(
cells=44,102; programs=50; guides=1,926; usage_matched=44,102; covariates=['grna_n_nonzero', 'response_n_umis', 'response_p_mito']
[0m[0m[0m[0msaved d3: observed target empirical points=14,900; NTC empirical points=7,750; output=results/ipsc_ec_qq_exports/k_50.dt_0_2/d3/d3_QQ_points.csv

Combined outputs:
  manifest: results/ipsc_ec_qq_exports/k_50.dt_0_2/QQ_export_manifest.csv
  all days QQ points: results/ipsc_ec_qq_exports/k_50.dt_0_2/all_days_QQ_points.csv
[0m
```

## Markdown cell 6

## Read exported points

## Code cell 7

```python
import pandas as pd

model_root = Path(OUT_DIR) / f'k_{K}.dt_{DENSITY_THRESHOLD}'
points_path = model_root / 'all_days_QQ_points.csv'
manifest_path = model_root / 'QQ_export_manifest.csv'

manifest = pd.read_csv(manifest_path)
display(manifest)

points = pd.read_csv(points_path)
display(points.head())
points.groupby(['day', 'pvalue_type', 'series_type']).size().to_frame('n_points')
```

### Text output (text/plain)

```text
  day     status                                 output_dir  \
0  d0  completed  results/ipsc_ec_qq_exports/k_50.dt_0_2/d0   
1  d1  completed  results/ipsc_ec_qq_exports/k_50.dt_0_2/d1   
2  d2  completed  results/ipsc_ec_qq_exports/k_50.dt_0_2/d2   
3  d3  completed  results/ipsc_ec_qq_exports/k_50.dt_0_2/d3   

                                           h5mu_path  \
0  /oak/stanford/groups/engreitz/Users/opushkar/c...   
1  /oak/stanford/groups/engreitz/Users/opushkar/c...   
2  /oak/stanford/groups/engreitz/Users/opushkar/c...   
3  /oak/stanford/groups/engreitz/Users/opushkar/c...   

                                          usage_path  \
0  /oak/stanford/groups/engreitz/Users/ymo/IGVF_c...   
1  /oak/stanford/groups/engreitz/Users/ymo/IGVF_c...   
2  /oak/stanford/groups/engreitz/Users/ymo/IGVF_c...   
3  /oak/stanford/groups/engreitz/Users/ymo/IGVF_c...   

                                      guide_long_csv  n_cells  n_programs  \
0  /oak/stanford/groups/engreitz/Users/tri/CodeBa...    38609          50   
1  /oak/stanford/groups/engreitz/Users/tri/CodeBa...    16925          50   
2  /oak/stanford/groups/engreitz/Users/tri/CodeBa...    21980          50   
3  /oak/stanford/groups/engreitz/Users/tri/CodeBa...    44102          50   

   n_guides  matched_usage_cells  n_observed_target_empirical_points  \
0      1926                38609                               14900   
1      1926                16925                               14900   
2      1926                21980                               14900   
3      1926                44102                               14900   

   n_ntc_pseudotarget_empirical_points  n_observed_target_skew_points  \
0                                 6100                          14900   
1                                 4650                          14900   
2                                 6350                          14900   
3                                 7750                          14900   

   n_ntc_pseudotarget_skew_points  message  \
0                            6100      NaN   
1                            4650      NaN   
2                            6350      NaN   
3                            7750      NaN   

                                       qq_points_csv  qq_points_rows  
0  results/ipsc_ec_qq_exports/k_50.dt_0_2/d0/d0_Q...           42200  
1  results/ipsc_ec_qq_exports/k_50.dt_0_2/d1/d1_Q...           39300  
2  results/ipsc_ec_qq_exports/k_50.dt_0_2/d2/d2_Q...           42700  
3  results/ipsc_ec_qq_exports/k_50.dt_0_2/d3/d3_Q...           45500  
```

### Text output (stream)

```text
/tmp/ipykernel_64318/3847999645.py:10: DtypeWarning: Columns (6) have mixed types. Specify dtype option on import or set low_memory=False.
  points = pd.read_csv(points_path)

```

### Text output (text/plain)

```text
            source_csv day pvalue_type                    series  \
0  d0/d0_QQ_points.csv  d0   empirical  Observed target p-values   
1  d0/d0_QQ_points.csv  d0   empirical  Observed target p-values   
2  d0/d0_QQ_points.csv  d0   empirical  Observed target p-values   
3  d0/d0_QQ_points.csv  d0   empirical  Observed target p-values   
4  d0/d0_QQ_points.csv  d0   empirical  Observed target p-values   

       series_type                                  source target_name  \
0  observed_target  observed target rows from pvals_raw_df      PAXBP1   
1  observed_target  observed target rows from pvals_raw_df       STAG2   
2  observed_target  observed target rows from pvals_raw_df        CTR9   
3  observed_target  observed target rows from pvals_raw_df       EIF4E   
4  observed_target  observed target rows from pvals_raw_df      EIF4G1   

   program_name  pvalue  rank  n_points  expected_pvalue  expected_neglog10_p  \
0             1  0.0001     1     14900         0.000034             4.474216   
1             1  0.0001     2     14900         0.000101             3.997095   
2             2  0.0001     3     14900         0.000168             3.775246   
3             2  0.0001     4     14900         0.000235             3.629118   
4             2  0.0001     5     14900         0.000302             3.519974   

   observed_neglog10_p  ensemble group_id  
0             4.000043       NaN      NaN  
1             4.000043       NaN      NaN  
2             4.000043       NaN      NaN  
3             4.000043       NaN      NaN  
4             4.000043       NaN      NaN  
```

### Text output (text/plain)

```text
                                        n_points
day pvalue_type series_type                     
d0  empirical   ntc_pseudotarget            6100
                observed_control_label       100
                observed_target            14900
    skew        ntc_pseudotarget            6100
                observed_control_label       100
                observed_target            14900
d1  empirical   ntc_pseudotarget            4650
                observed_control_label       100
                observed_target            14900
    skew        ntc_pseudotarget            4650
                observed_control_label       100
                observed_target            14900
d2  empirical   ntc_pseudotarget            6350
                observed_control_label       100
                observed_target            14900
    skew        ntc_pseudotarget            6350
                observed_control_label       100
                observed_target            14900
d3  empirical   ntc_pseudotarget            7750
                observed_control_label       100
                observed_target            14900
    skew        ntc_pseudotarget            7750
                observed_control_label       100
                observed_target            14900
```

## Markdown cell 8

## Minimal replot example

The generated PDFs/PNGs are already saved by the script. This cell shows how to replot from the CSV only.

## Code cell 9

```python
import matplotlib.pyplot as plt

plot_data = points.loc[
    (points['pvalue_type'] == 'empirical')
    & (points['series_type'].isin(['observed_target', 'ntc_pseudotarget']))
].copy()

for day, current in plot_data.groupby('day', sort=True):
    fig, ax = plt.subplots(figsize=(5.5, 5.2))
    for series_type, label in [
        ('ntc_pseudotarget', 'Non-targeting pseudo-targets'),
        ('observed_target', 'Observed targets'),
    ]:
        sub = current.loc[current['series_type'] == series_type]
        ax.scatter(sub['expected_neglog10_p'], sub['observed_neglog10_p'], s=10, alpha=0.65, label=label)
    lim = max(1.0, current['expected_neglog10_p'].max(), current['observed_neglog10_p'].max()) * 1.03
    ax.plot([0, lim], [0, lim], linestyle='--', linewidth=1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('Expected -log10(p)')
    ax.set_ylabel('Observed -log10(p)')
    ax.set_title(f'iPSC-EC {day}: empirical CRT QQ')
    ax.legend(frameon=False)
    plt.show()
```

### Text output (text/plain)

```text
<Figure size 550x520 with 1 Axes>
```

### Text output (text/plain)

```text
<Figure size 550x520 with 1 Axes>
```

### Text output (text/plain)

```text
<Figure size 550x520 with 1 Axes>
```

### Text output (text/plain)

```text
<Figure size 550x520 with 1 Axes>
```

## Code cell 10

```python
from pathlib import Path
from IPython.display import display
from plotnine import theme, guides, guide_legend, element_blank, element_text
import matplotlib as mpl

# Keep text editable in Illustrator vector exports.
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"


def qq_legend_above_panel(p, *, width=4.6, height=4.9):
    """
    Override the old in-panel legend placement.
    Places legend above the panel so it does not cover QQ points.
    """
    return (
        p
        + guides(
            color=guide_legend(
                nrow=1,
                override_aes={"size": 3.6, "alpha": 1},
            )
        )
        + theme(
            figure_size=(width, height),

            # Key change: outside the plotting panel
            legend_position="top",
            legend_direction="horizontal",

            # Remove the boxed legend look
            legend_background=element_blank(),
            legend_key=element_blank(),
            legend_title=element_blank(),
            legend_text=element_text(size=9),

            # Give the title/subtitle/legend enough room
            plot_title=element_text(size=12, weight="bold", ha="center"),
            plot_subtitle=element_text(size=9, weight="bold", ha="center"),
        )
    )


def save_illustrator_formats(
    p,
    output_stem,
    *,
    width=4.6,
    height=4.9,
    dpi=300,
    formats=("pdf", "svg", "eps", "png"),
):
    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)

    saved = {}
    for fmt in formats:
        out = output_stem.with_suffix(f".{fmt}")
        kwargs = dict(width=width, height=height, units="in", verbose=False)
        if fmt == "png":
            kwargs["dpi"] = dpi

        p.save(out, **kwargs)
        saved[fmt] = out

    return saved


qq_csv = Path("results/ipsc_ec_qq_exports/k_50.dt_0_2/all_days_QQ_points.csv")
illustrator_out_dir = qq_csv.parent / "illustrator_editable"

all_saved_files = {}

for d in ["d0", "d1", "d2", "d3"]:

    # Important: output_stem=None here, otherwise the old in-panel legend version
    # gets saved before we apply the legend override.
    p = plot_ipsc_ec_qq(
        qq_csv,
        day=d,
        pvalue_type="empirical",
        title=f"{d}: Real target and NTC empirical p-values",
        subtitle="RNA QC, guide count and guide UMI burden",
        output_stem=None,
    )

    p = qq_legend_above_panel(p, width=4.6, height=4.9)

    saved = save_illustrator_formats(
        p,
        illustrator_out_dir / f"{d}_QQ_plotnine_legend_above",
        width=4.6,
        height=4.9,
        formats=("pdf", "svg", "eps", "png"),
    )

    all_saved_files[d] = saved
    display(p)

print("Saved files:")
for day, files in all_saved_files.items():
    print(f"\n{day}")
    for fmt, path in files.items():
        print(f"  {fmt.upper()}: {path}")
```

### Text output (stream)

```text
/tmp/ipykernel_64318/2500159175.py:59: DtypeWarning: Columns (6) have mixed types. Specify dtype option on import or set low_memory=False.
The PostScript backend does not support transparency; partially transparent artists will be rendered opaque.

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f3dcc100690>
```

### Text output (stream)

```text
/tmp/ipykernel_64318/2500159175.py:59: DtypeWarning: Columns (6) have mixed types. Specify dtype option on import or set low_memory=False.
The PostScript backend does not support transparency; partially transparent artists will be rendered opaque.

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f3dcd525550>
```

### Text output (stream)

```text
/tmp/ipykernel_64318/2500159175.py:59: DtypeWarning: Columns (6) have mixed types. Specify dtype option on import or set low_memory=False.
The PostScript backend does not support transparency; partially transparent artists will be rendered opaque.

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f3dcd871f50>
```

### Text output (stream)

```text
/tmp/ipykernel_64318/2500159175.py:59: DtypeWarning: Columns (6) have mixed types. Specify dtype option on import or set low_memory=False.
The PostScript backend does not support transparency; partially transparent artists will be rendered opaque.

```

### Text output (text/plain)

```text
<plotnine.ggplot.ggplot object at 0x7f3dc57b0fd0>
```

### Text output (stream)

```text
Saved files:

d0
  PDF: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d0_QQ_plotnine_legend_above.pdf
  SVG: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d0_QQ_plotnine_legend_above.svg
  EPS: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d0_QQ_plotnine_legend_above.eps
  PNG: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d0_QQ_plotnine_legend_above.png

d1
  PDF: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d1_QQ_plotnine_legend_above.pdf
  SVG: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d1_QQ_plotnine_legend_above.svg
  EPS: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d1_QQ_plotnine_legend_above.eps
  PNG: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d1_QQ_plotnine_legend_above.png

d2
  PDF: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d2_QQ_plotnine_legend_above.pdf
  SVG: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d2_QQ_plotnine_legend_above.svg
  EPS: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d2_QQ_plotnine_legend_above.eps
  PNG: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d2_QQ_plotnine_legend_above.png

d3
  PDF: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d3_QQ_plotnine_legend_above.pdf
  SVG: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d3_QQ_plotnine_legend_above.svg
  EPS: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d3_QQ_plotnine_legend_above.eps
  PNG: results/ipsc_ec_qq_exports/k_50.dt_0_2/illustrator_editable/d3_QQ_plotnine_legend_above.png

```

