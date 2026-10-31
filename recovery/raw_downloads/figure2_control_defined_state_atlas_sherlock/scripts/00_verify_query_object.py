#!/usr/bin/env python3
from pathlib import Path
import argparse, yaml
p=argparse.ArgumentParser(description="Verify configured CC-Perturb-seq query object.")
p.add_argument('--config', required=True)
a=p.parse_args()
cfg=yaml.safe_load(Path(a.config).read_text())
q=Path(cfg['query_h5ad'])
print(f'query_h5ad: {q}')
if not q.exists():
    raise SystemExit(f'ERROR: query H5AD is not visible: {q}')
print(f'OK: file exists; size={q.stat().st_size/1024**3:.2f} GiB')
try:
    import anndata as ad
    x=ad.read_h5ad(q, backed='r')
    print(f'shape: {x.n_obs:,} cells x {x.n_vars:,} genes')
    print('obs keys:', list(x.obs_keys()))
    print('obsm keys:', list(x.obsm_keys()))
    print('layers:', list(x.layers.keys()))
    print('guide_names in uns:', 'guide_names' in x.uns)
    x.file.close()
except Exception as e:
    print(f'WARNING: metadata inspection failed, but path exists: {type(e).__name__}: {e}')
