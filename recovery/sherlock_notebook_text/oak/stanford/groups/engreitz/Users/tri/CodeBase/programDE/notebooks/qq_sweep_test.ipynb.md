# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/qq_sweep_test.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/qq_sweep_test.ipynb`
- SHA-256: `ba8097863a745a4effd34f46c50b6d6ac356c01dd187c1eddc5eb4989b6b2b6f`
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

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/muon/_core/preproc.py:31: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
  if Version(scanpy.__version__) < Version("1.10"):

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
!python ../scripts/ipsc_ec_qq_sweep.py \
  --data-dir /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/2.normalized_data/normalized_data \
  --cnmf-dir /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7/012726_100k_cells_20iter_allHVG_torch_halsvar_batch_e7_all \
  --guide-long-csv-template '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/{day}_sceptre_transfer/{day}_grna_assignments_long.csv' \
  --guide-metadata-tsv ./300genes_guide_metadata_v43.tsv \
  --out-dir results/ipsc_ec_qq \
  --ntc-label SAFE --ntc-label non-targeting --ntc-label NTC \
  --B 1023 \
  --n-jobs 4 \
  --skip-if-exists
```

### Text output (stream)

```text
Found 4 dataset files
Found 28 usage files

=== d0_normalized_gex_and_guide.h5mu | k=10 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/muon/_core/preproc.py:31: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
  if Version(scanpy.__version__) < Version("1.10"):
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 10), matched_cells=38609/38609
[saved] results/ipsc_ec_qq/qq_d0.k_10.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=10 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 10), matched_cells=16925/16925
[saved] results/ipsc_ec_qq/qq_d1.k_10.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=10 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 10), matched_cells=21980/21980
[saved] results/ipsc_ec_qq/qq_d2.k_10.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=10 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 10), matched_cells=44102/44102
[saved] results/ipsc_ec_qq/qq_d3.k_10.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=10 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 10), matched_cells=38609/38609
[saved] results/ipsc_ec_qq/qq_d0.k_10.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=10 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 10), matched_cells=16925/16925
[saved] results/ipsc_ec_qq/qq_d1.k_10.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=10 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 10), matched_cells=21980/21980
[saved] results/ipsc_ec_qq/qq_d2.k_10.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=10 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 10), matched_cells=44102/44102
[saved] results/ipsc_ec_qq/qq_d3.k_10.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=20 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 20), matched_cells=38609/38609
[saved] results/ipsc_ec_qq/qq_d0.k_20.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=20 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 20), matched_cells=16925/16925
[saved] results/ipsc_ec_qq/qq_d1.k_20.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=20 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 20), matched_cells=21980/21980
[saved] results/ipsc_ec_qq/qq_d2.k_20.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=20 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 20), matched_cells=44102/44102
[saved] results/ipsc_ec_qq/qq_d3.k_20.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=20 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 20), matched_cells=38609/38609
[saved] results/ipsc_ec_qq/qq_d0.k_20.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=20 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 20), matched_cells=16925/16925
[saved] results/ipsc_ec_qq/qq_d1.k_20.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=20 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 20), matched_cells=21980/21980
[saved] results/ipsc_ec_qq/qq_d2.k_20.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=20 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 20), matched_cells=44102/44102
[saved] results/ipsc_ec_qq/qq_d3.k_20.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=30 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 30), matched_cells=38609/38609
[saved] results/ipsc_ec_qq/qq_d0.k_30.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=30 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 30), matched_cells=16925/16925
[saved] results/ipsc_ec_qq/qq_d1.k_30.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=30 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 30), matched_cells=21980/21980
[saved] results/ipsc_ec_qq/qq_d2.k_30.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=30 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 30), matched_cells=44102/44102
[saved] results/ipsc_ec_qq/qq_d3.k_30.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=30 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 30), matched_cells=38609/38609
[saved] results/ipsc_ec_qq/qq_d0.k_30.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=30 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 30), matched_cells=16925/16925
[saved] results/ipsc_ec_qq/qq_d1.k_30.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=30 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 30), matched_cells=21980/21980
[saved] results/ipsc_ec_qq/qq_d2.k_30.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=30 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 30), matched_cells=44102/44102
[saved] results/ipsc_ec_qq/qq_d3.k_30.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=40 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 40), matched_cells=38609/38609
[saved] results/ipsc_ec_qq/qq_d0.k_40.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=40 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 40), matched_cells=16925/16925
[saved] results/ipsc_ec_qq/qq_d1.k_40.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=40 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 40), matched_cells=21980/21980
[saved] results/ipsc_ec_qq/qq_d2.k_40.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=40 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 40), matched_cells=44102/44102
[saved] results/ipsc_ec_qq/qq_d3.k_40.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=40 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 40), matched_cells=38609/38609
[saved] results/ipsc_ec_qq/qq_d0.k_40.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=40 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 40), matched_cells=16925/16925
[saved] results/ipsc_ec_qq/qq_d1.k_40.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=40 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 40), matched_cells=21980/21980
[saved] results/ipsc_ec_qq/qq_d2.k_40.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=40 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 40), matched_cells=44102/44102
[saved] results/ipsc_ec_qq/qq_d3.k_40.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=50 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 50), matched_cells=38609/38609
[saved] results/ipsc_ec_qq/qq_d0.k_50.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=50 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 50), matched_cells=16925/16925
[saved] results/ipsc_ec_qq/qq_d1.k_50.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=50 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 50), matched_cells=21980/21980
[saved] results/ipsc_ec_qq/qq_d2.k_50.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=50 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 50), matched_cells=44102/44102
[saved] results/ipsc_ec_qq/qq_d3.k_50.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=50 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 50), matched_cells=38609/38609
[saved] results/ipsc_ec_qq/qq_d0.k_50.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=50 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 50), matched_cells=16925/16925
[saved] results/ipsc_ec_qq/qq_d1.k_50.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=50 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 50), matched_cells=21980/21980
[saved] results/ipsc_ec_qq/qq_d2.k_50.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=50 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 50), matched_cells=44102/44102
[saved] results/ipsc_ec_qq/qq_d3.k_50.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=60 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 60), matched_cells=38609/38609
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_60.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=60 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 60), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_60.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=60 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 60), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_60.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=60 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 60), matched_cells=44102/44102
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_60.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=60 dt=2_0 ===
[0m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 60), matched_cells=38609/38609
[0m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/joblib/externals/loky/process_executor.py:782: UserWarning: A worker stopped while some jobs were given to the executor. This can be caused by a too short worker timeout or by a memory leak.
  warnings.warn(
[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_60.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=60 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 60), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_60.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=60 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 60), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_60.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=60 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 60), matched_cells=44102/44102
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_60.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=70 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 70), matched_cells=38609/38609
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_70.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=70 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 70), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_70.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=70 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 70), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_70.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=70 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 70), matched_cells=44102/44102
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_70.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=70 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 70), matched_cells=38609/38609
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_70.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=70 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 70), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_70.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=70 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 70), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_70.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=70 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 70), matched_cells=44102/44102
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_70.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=80 dt=0_2 ===
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 80), matched_cells=38609/38609
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_80.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=80 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 80), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_80.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=80 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 80), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_80.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=80 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 80), matched_cells=44102/44102
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_80.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=80 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 80), matched_cells=38609/38609
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_80.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=80 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 80), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_80.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=80 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 80), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_80.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=80 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 80), matched_cells=44102/44102
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_80.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=90 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 90), matched_cells=38609/38609
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_90.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=90 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 90), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_90.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=90 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 90), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_90.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=90 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 90), matched_cells=44102/44102
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[0m[saved] results/ipsc_ec_qq/qq_d3.k_90.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=90 dt=2_0 ===
[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 90), matched_cells=38609/38609
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_90.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=90 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 90), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_90.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=90 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 90), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_90.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=90 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 90), matched_cells=44102/44102
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[0m[saved] results/ipsc_ec_qq/qq_d3.k_90.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=100 dt=0_2 ===
[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 100), matched_cells=38609/38609
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_100.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=100 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 100), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_100.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=100 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 100), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_100.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=100 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 100), matched_cells=44102/44102
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_100.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=100 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 100), matched_cells=38609/38609
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_100.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=100 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 100), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_100.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=100 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 100), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_100.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=100 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 100), matched_cells=44102/44102
[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[0m[0m[0m[saved] results/ipsc_ec_qq/qq_d3.k_100.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=150 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 150), matched_cells=38609/38609
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_150.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=150 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 150), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_150.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=150 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 150), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_150.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=150 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 150), matched_cells=44102/44102
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_150.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=150 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 150), matched_cells=38609/38609
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_150.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=150 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 150), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_150.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=150 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 150), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_150.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=150 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 150), matched_cells=44102/44102
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_150.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=200 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 200), matched_cells=38609/38609
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_200.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=200 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 200), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_200.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=200 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 200), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[0m[saved] results/ipsc_ec_qq/qq_d2.k_200.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=200 dt=0_2 ===
[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 200), matched_cells=44102/44102
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_200.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=200 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 200), matched_cells=38609/38609
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_200.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=200 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 200), matched_cells=16925/16925
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_200.dt_2_0.png

=== d2_normalized_gex_and_guide.h5mu | k=200 dt=2_0 ===
[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
[0musage: shape=(21980, 200), matched_cells=21980/21980
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[0m[0m[0m[saved] results/ipsc_ec_qq/qq_d2.k_200.dt_2_0.png

=== d3_normalized_gex_and_guide.h5mu | k=200 dt=2_0 ===
[0m/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 200), matched_cells=44102/44102
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_200.dt_2_0.png

=== d0_normalized_gex_and_guide.h5mu | k=250 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 250), matched_cells=38609/38609
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_250.dt_0_2.png

=== d1_normalized_gex_and_guide.h5mu | k=250 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 250), matched_cells=16925/16925
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d1.k_250.dt_0_2.png

=== d2_normalized_gex_and_guide.h5mu | k=250 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(21980, 250), matched_cells=21980/21980
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d2.k_250.dt_0_2.png

=== d3_normalized_gex_and_guide.h5mu | k=250 dt=0_2 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(44102, 250), matched_cells=44102/44102
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d3.k_250.dt_0_2.png

=== d0_normalized_gex_and_guide.h5mu | k=250 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(38609, 250), matched_cells=38609/38609
[0m[0m[0m[0m/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/../scripts/ipsc_ec_qq_sweep.py:566: UserWarning: Creating legend with loc="best" can be slow with large amounts of data.
  ax.figure.savefig(out_path, dpi=300, bbox_inches="tight")
[saved] results/ipsc_ec_qq/qq_d0.k_250.dt_2_0.png

=== d1_normalized_gex_and_guide.h5mu | k=250 dt=2_0 ===
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1598: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("var", axis=0, join_common=join_common)
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/lib/python3.11/site-packages/mudata/_core/mudata.py:1461: FutureWarning: From 0.4 .update() will not pull obs/var columns from individual modalities by default anymore. Set mudata.set_options(pull_on_update=False) to adopt the new behaviour, which will become the default. Use new pull_obs/pull_var and push_obs/push_var methods for more flexibility.
  self._update_attr("obs", axis=1, join_common=join_common)
usage: shape=(16925, 250), matched_cells=16925/16925

```

## Code cell 4

```python

```

