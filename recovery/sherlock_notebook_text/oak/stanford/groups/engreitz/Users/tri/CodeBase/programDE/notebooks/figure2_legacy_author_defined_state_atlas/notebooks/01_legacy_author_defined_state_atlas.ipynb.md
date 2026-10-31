# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_legacy_author_defined_state_atlas/notebooks/01_legacy_author_defined_state_atlas.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_legacy_author_defined_state_atlas/notebooks/01_legacy_author_defined_state_atlas.ipynb`
- SHA-256: `328fe4b345ab68955457812e5714a05294ff42c8a0a94438cafc4b2c74e0cb80`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Legacy author-defined coarse-state atlas

This notebook uses the previous `leiden_res_0.5` mapping directly and does not require cluster review.

**Caveat:** the original Leiden boundaries were learned from all cells, so this is not a strictly control-defined atlas.


## Code cell 1

```python
from pathlib import Path
import sys, yaml
PACKAGE_DIR = Path.cwd().resolve()
if not (PACKAGE_DIR / 'src').exists():
    PACKAGE_DIR = Path('..').resolve()
sys.path.insert(0, str(PACKAGE_DIR / 'src'))
CONFIG_LOCAL = PACKAGE_DIR / 'config' / 'legacy_author_state_atlas_config.local.yaml'
print('PACKAGE_DIR:', PACKAGE_DIR)
print('CONFIG_LOCAL:', CONFIG_LOCAL)

```

### Text output (stream)

```text
PACKAGE_DIR: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_legacy_author_defined_state_atlas
CONFIG_LOCAL: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_legacy_author_defined_state_atlas/config/legacy_author_state_atlas_config.local.yaml

```

## Code cell 2

```python
cfg = yaml.safe_load(CONFIG_LOCAL.read_text())
for key in ['query_h5ad','query_cluster_col','query_umap_key','legacy_cluster_state_map_tsv','output_dir']:
    print(f'{key}: {cfg.get(key)}')
query_path = Path(cfg['query_h5ad'])
assert query_path.exists(), f'Query H5AD does not exist: {query_path}'

```

### Text output (stream)

```text
query_h5ad: /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad
query_cluster_col: leiden_res_0.5
query_umap_key: X_umap
legacy_cluster_state_map_tsv: templates/legacy_leiden_res_0.5_state_map.tsv
output_dir: legacy_author_defined_state_atlas_results

```

## Markdown cell 3

## Run the atlas
This step maps the existing Leiden clusters to the previous author-defined names, exports all-cell assignments, identifies pure controls, and creates the figures.


## Code cell 4

```python
from legacy_author_defined_state_atlas import run_legacy_author_atlas
OUTDIR = run_legacy_author_atlas(CONFIG_LOCAL)
print('Outputs written to:', OUTDIR)

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_legacy_author_defined_state_atlas/src/legacy_author_defined_state_atlas.py:336: UserWarning: This figure includes Axes that are not compatible with tight_layout, so results might be incorrect.
  fig.tight_layout()
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_legacy_author_defined_state_atlas/src/legacy_author_defined_state_atlas.py:336: UserWarning: This figure includes Axes that are not compatible with tight_layout, so results might be incorrect.
  fig.tight_layout()
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_legacy_author_defined_state_atlas/src/legacy_author_defined_state_atlas.py:512: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('scanpy')` instead
  "scanpy": sc.__version__,
/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_legacy_author_defined_state_atlas/src/legacy_author_defined_state_atlas.py:513: FutureWarning: `__version__` is deprecated, use `importlib.metadata.version('anndata')` instead.
  "anndata": ad.__version__,

```

### Text output (stream)

```text
Outputs written to: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_legacy_author_defined_state_atlas/legacy_author_defined_state_atlas_results

```

## Code cell 5

```python
from IPython.display import display
import pandas as pd
display(pd.read_csv(Path(OUTDIR) / 'state_atlas_qc.tsv', sep='\t'))
display(pd.read_csv(Path(OUTDIR) / 'state_counts_by_day_and_replicate.tsv', sep='\t').head(30))

```

### Text output (text/plain)

```text
                 population  n_cells  n_states  \
0  all_CC-Perturb-seq_cells   121616         9   
1             pure_controls      748         8   

                     state_definition_source        state_mapping_source  \
0  legacy leiden_res_0.5 all-cell clustering  ipsc-ec.Hippo.figure.ipynb   
1  legacy leiden_res_0.5 all-cell clustering  ipsc-ec.Hippo.figure.ipynb   

   review_bypassed_by_user_request  
0                             True  
1                             True  
```

### Text output (text/plain)

```text
                  population observed_day biological_replicate assigned_state  \
0   all_CC-Perturb-seq_cells           D0              d0_rep1           iPSC   
1   all_CC-Perturb-seq_cells           D0              d0_rep1     Neuro-iPSC   
2   all_CC-Perturb-seq_cells           D0              d0_rep1         Others   
3   all_CC-Perturb-seq_cells           D0              d0_rep1    Fibro-Mesen   
4   all_CC-Perturb-seq_cells           D0              d0_rep1      Neu-Mesen   
5   all_CC-Perturb-seq_cells           D0              d0_rep1     PrimStreak   
6   all_CC-Perturb-seq_cells           D0              d0_rep2     Neuro-iPSC   
7   all_CC-Perturb-seq_cells           D0              d0_rep2           iPSC   
8   all_CC-Perturb-seq_cells           D0              d0_rep2         Others   
9   all_CC-Perturb-seq_cells           D0              d0_rep2    Fibro-Mesen   
10  all_CC-Perturb-seq_cells           D1              d1_rep1     PrimStreak   
11  all_CC-Perturb-seq_cells           D1              d1_rep1         Others   
12  all_CC-Perturb-seq_cells           D1              d1_rep1           iPSC   
13  all_CC-Perturb-seq_cells           D1              d1_rep1       Lat-Meso   
14  all_CC-Perturb-seq_cells           D1              d1_rep1     Neuro-iPSC   
15  all_CC-Perturb-seq_cells           D1              d1_rep1    Fibro-Mesen   
16  all_CC-Perturb-seq_cells           D1              d1_rep2     PrimStreak   
17  all_CC-Perturb-seq_cells           D1              d1_rep2           iPSC   
18  all_CC-Perturb-seq_cells           D1              d1_rep2         Others   
19  all_CC-Perturb-seq_cells           D2              d2_rep1       Lat-Meso   
20  all_CC-Perturb-seq_cells           D2              d2_rep1         Pre-EC   
21  all_CC-Perturb-seq_cells           D2              d2_rep1    Fibro-Mesen   
22  all_CC-Perturb-seq_cells           D2              d2_rep1      Neu-Mesen   
23  all_CC-Perturb-seq_cells           D2              d2_rep1          artEC   
24  all_CC-Perturb-seq_cells           D2              d2_rep1           iPSC   
25  all_CC-Perturb-seq_cells           D2              d2_rep1     Neuro-iPSC   
26  all_CC-Perturb-seq_cells           D2              d2_rep1     PrimStreak   
27  all_CC-Perturb-seq_cells           D2              d2_rep1         Others   
28  all_CC-Perturb-seq_cells           D2              d2_rep2       Lat-Meso   
29  all_CC-Perturb-seq_cells           D2              d2_rep2         Pre-EC   

    n_cells  fraction  row_total  
0      9543  0.545688      17488  
1      7921  0.452939      17488  
2        12  0.000686      17488  
3        10  0.000572      17488  
4         1  0.000057      17488  
5         1  0.000057      17488  
6     10869  0.514606      21121  
7     10191  0.482506      21121  
8        54  0.002557      21121  
9         7  0.000331      21121  
10     5840  0.997268       5856  
11       10  0.001708       5856  
12        3  0.000512       5856  
13        1  0.000171       5856  
14        1  0.000171       5856  
15        1  0.000171       5856  
16    11046  0.997922      11069  
17       13  0.001174      11069  
18       10  0.000903      11069  
19     7740  0.710418      10895  
20     3089  0.283525      10895  
21       20  0.001836      10895  
22       15  0.001377      10895  
23       12  0.001101      10895  
24       10  0.000918      10895  
25        4  0.000367      10895  
26        3  0.000275      10895  
27        2  0.000184      10895  
28     8352  0.753451      11085  
29     2676  0.241407      11085  
```

## Markdown cell 6

## Output files
Open `Fig2_control_defined_state_atlas.pdf` for the selected panel. The downstream occupancy workflow should use `control_state_assignments.tsv`.


