# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/notebooks/01_generate_perturbation_class_occupancy.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/notebooks/01_generate_perturbation_class_occupancy.ipynb`
- SHA-256: `23b534d4d8df124f062b1f3788dd4a5a2b54fdc87a8534d10dbc5ff80ef6d710`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Figure 2: perturbation-class coarse-state occupancy

This notebook creates the requested **perturbation-class occupancy bubble matrix** and the complete target-level supplementary heatmap from the processed CC-Perturb-seq H5AD on Sherlock.

It intentionally does **not** use target→program CRT effects, cNMF differential effects, or perturbation-driven clustering. State assignments are fixed before testing using the existing author-defined mapping of `leiden_res_0.5` recovered from `ipsc-ec.Hippo.figure.ipynb`:

| cluster | state |
|---:|---|
| 0 | artEC |
| 1 | iPSC |
| 2 | Neuro-iPSC |
| 3 | PrimStreak |
| 4 | Lat-Meso |
| 5 | Fibro-Mesen |
| 6 | Pre-EC |
| 7 | Neu-Mesen |
| 8 | Others |

### Primary model

For every day, target, and state, the notebook fits an adjusted one-vs-rest logistic regression:

```text
state membership ~ target-union carrier indicator
                 + biological sample
                 + pct_counts_mt
                 + log1p(total_counts)
                 + log1p(guides_per_cell)
                 + log1p(guide_umi_counts)
                 + log1p(n_genes_by_counts)
```

The main bubble matrix summarizes target-level models within reviewed regulator classes. The analysis never pools all guides in a class into one exposure.

## Markdown cell 1

## 0. Sherlock setup

Launch this notebook from the package root after obtaining an interactive allocation with enough memory and CPUs, for example:

```bash
srun --pty --time=24:00:00 --mem=128G --cpus-per-task=8 bash
source /path/to/your/environment/bin/activate
jupyter lab --no-browser --port=8888
```

The package also includes `scripts/run_full_occupancy_pipeline.sbatch` for a non-interactive run.

## Code cell 2

```python
from pathlib import Path
import os, sys, json, subprocess, zipfile
import yaml
import pandas as pd
from IPython.display import display, Image, Markdown

PACKAGE_DIR = Path.cwd().resolve()
if PACKAGE_DIR.name == "notebooks":
    PACKAGE_DIR = PACKAGE_DIR.parent

CONFIG_PATH = PACKAGE_DIR / "config" / "perturbation_class_occupancy_config.yaml"
EXPORT_DIR = PACKAGE_DIR / "occupancy_compact_export"
RESULT_DIR = PACKAGE_DIR / "occupancy_analysis_results"
N_JOBS = 8
FORCE_REEXPORT = False
FORCE_REFIT = False

print("PACKAGE_DIR:", PACKAGE_DIR)
print("CONFIG_PATH:", CONFIG_PATH)
print("EXPORT_DIR:", EXPORT_DIR)
print("RESULT_DIR:", RESULT_DIR)

```

### Text output (stream)

```text
PACKAGE_DIR: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock
CONFIG_PATH: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/config/perturbation_class_occupancy_config.yaml
EXPORT_DIR: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/occupancy_compact_export
RESULT_DIR: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/occupancy_analysis_results

```

## Markdown cell 3

## 1. Configuration and preflight audit

The verified processed query object from the uploaded notebooks is:

```text
/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad
```

The analysis uses `intended_target_name` from the bundled guide metadata and a target-union carrier indicator.

## Code cell 4

```python
cfg = yaml.safe_load(CONFIG_PATH.read_text())
for key, value in cfg.items():
    print(f"{key}: {value}")

query_path = Path(cfg["query_h5ad"])
assert query_path.exists(), f"Query H5AD does not exist: {query_path}"

for rel_key in ["guide_metadata_tsv", "legacy_state_map_tsv", "class_membership_tsv", "state_definitions_tsv"]:
    p = Path(cfg[rel_key])
    if not p.is_absolute():
        p = PACKAGE_DIR / p
    assert p.exists(), f"Missing required package input: {p}"

print("\nAll configured inputs exist.")

```

### Text output (stream)

```text
query_h5ad: /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad
guide_metadata_tsv: templates/300genes_guide_metadata_v43.tsv
legacy_state_map_tsv: templates/legacy_leiden_res_0.5_state_map.tsv
class_membership_tsv: templates/class_membership_reviewed_screen_hits.tsv
state_definitions_tsv: templates/state_definitions.tsv
output_dir: occupancy_compact_export
cluster_col: leiden_res_0.5
day_col: day
replicate_col: biological_sample
guide_assignment_key: guide_assignment
guide_names_key: guide_names
covariates: ['pct_counts_mt', 'biological_sample', 'total_counts', 'guides_per_cell', 'guide_umi_counts', 'n_genes_by_counts']
minimum_carrier_cells_per_day: 25
minimum_state_cells_per_day: 25
minimum_targets_per_class: 4
fdr_threshold: 0.05
fdr_family: within_day
analysis_population: all_cells_marginal
strict_control_sensitivity: True

All configured inputs exist.

```

## Code cell 5

```python
# Read only metadata/backed structure for a fast audit.
import anndata as ad
adata = ad.read_h5ad(query_path, backed="r")
print(adata)
print("\nobs columns:", adata.obs.columns.tolist())
print("\nobsm keys:", list(adata.obsm.keys()))
print("\nuns keys:", list(adata.uns.keys()))

required_obs = [
    cfg["cluster_col"], cfg["day_col"], cfg["replicate_col"],
    "pct_counts_mt", "total_counts", "guides_per_cell",
    "guide_umi_counts", "n_genes_by_counts",
]
missing_obs = [c for c in required_obs if c not in adata.obs.columns]
print("\nMissing required obs columns:", missing_obs)
assert not missing_obs
assert cfg["guide_assignment_key"] in adata.obsm
assert cfg["guide_names_key"] in adata.uns

display(adata.obs.groupby([cfg["day_col"], cfg["replicate_col"]], observed=True).size().rename("n_cells").reset_index())
display(adata.obs[cfg["cluster_col"]].astype(str).value_counts().sort_index().rename("n_cells").reset_index())

try:
    adata.file.close()
except Exception:
    pass

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

obs columns: ['guides_per_cell', 'guide_umi_counts', 'pct_counts_ribo', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'pool', 'sample_id', 'round1', 'well', 'biological_sample', 'plate', 'bci', 'stype', 'day', 'cell_pool', 'log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'log1p_total_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'S_score', 'G2M_score', 'phase', 'outlier', 'size_factors', 'diff_day', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2']

obsm keys: ['X_pca', 'X_umap', 'guide_assignment']

uns keys: ['guide_names', 'leiden_res_0.0', 'leiden_res_0.1', 'leiden_res_0.2', 'leiden_res_0.3', 'leiden_res_0.4', 'leiden_res_0.5', 'leiden_res_0.6', 'leiden_res_0.7', 'leiden_res_0.8', 'leiden_res_0.9', 'leiden_res_1.0', 'leiden_res_1.1', 'leiden_res_1.2', 'log1p', 'neighbors']

Missing required obs columns: []

```

### Text output (text/plain)

```text
  day biological_sample  n_cells
0  d0           d0_rep1    17488
1  d0           d0_rep2    21121
2  d1           d1_rep1     5856
3  d1           d1_rep2    11069
4  d2           d2_rep1    10895
5  d2           d2_rep2    11085
6  d3           d3_rep1    21377
7  d3           d3_rep2    22725
```

### Text output (text/plain)

```text
  leiden_res_0.5  n_cells
0              0    28458
1              1    19780
2              2    18799
3              3    16895
4              4    16168
5              5    12295
6              6     5815
7              7     3313
8              8       93
```

## Markdown cell 6

## 2. Review the fixed states and regulator classes

These tables are inputs, not results. The state boundaries are not re-estimated. Class membership is the literature-informed reviewed screen-hit annotation assembled earlier in this project; only targets represented in both the guide library and class table contribute to the class-level panel.

## Code cell 7

```python
def resolve_package_path(value):
    p = Path(value)
    return p if p.is_absolute() else PACKAGE_DIR / p

state_map = pd.read_csv(resolve_package_path(cfg["legacy_state_map_tsv"]), sep="\t", dtype={"source_cluster": str})
state_defs = pd.read_csv(resolve_package_path(cfg["state_definitions_tsv"]), sep="\t")
class_membership = pd.read_csv(resolve_package_path(cfg["class_membership_tsv"]), sep="\t")
guide_meta = pd.read_csv(resolve_package_path(cfg["guide_metadata_tsv"]), sep="\t")

print("Fixed cluster-to-state map")
display(state_map)
print("State definitions")
display(state_defs.sort_values("state_order"))
print("Reviewed regulator-class counts")
display(class_membership.groupby("regulator_class")["gene"].nunique().sort_values(ascending=False).rename("n_targets").reset_index())
print("Guide types")
display(guide_meta["type"].value_counts(dropna=False).rename("n_guides").reset_index())

```

### Text output (stream)

```text
Fixed cluster-to-state map

```

### Text output (text/plain)

```text
  source_cluster assigned_state  state_order                    state_label  \
0              0          artEC            6                    Arterial EC   
1              1           iPSC            0                          hiPSC   
2              2     Neuro-iPSC            1        Neuroectoderm-like iPSC   
3              3     PrimStreak            2               Primitive streak   
4              4       Lat-Meso            3               Lateral mesoderm   
5              5    Fibro-Mesen            7  Fibro-mesenchymal alternative   
6              6         Pre-EC            4   Pre-EC / vascular progenitor   
7              7      Neu-Mesen            5             Neural-mesenchymal   
8              8         Others            8             Other / unresolved   

               mapping_source  
0  ipsc-ec.Hippo.figure.ipynb  
1  ipsc-ec.Hippo.figure.ipynb  
2  ipsc-ec.Hippo.figure.ipynb  
3  ipsc-ec.Hippo.figure.ipynb  
4  ipsc-ec.Hippo.figure.ipynb  
5  ipsc-ec.Hippo.figure.ipynb  
6  ipsc-ec.Hippo.figure.ipynb  
7  ipsc-ec.Hippo.figure.ipynb  
8  ipsc-ec.Hippo.figure.ipynb  
```

### Text output (stream)

```text
State definitions

```

### Text output (text/plain)

```text
  assigned_state                    state_label  state_order  \
1           iPSC                          hiPSC            0   
2     Neuro-iPSC        Neuroectoderm-like iPSC            1   
3     PrimStreak               Primitive streak            2   
4       Lat-Meso               Lateral mesoderm            3   
6         Pre-EC   Pre-EC / vascular progenitor            4   
7      Neu-Mesen             Neural-mesenchymal            5   
0          artEC                    Arterial EC            6   
5    Fibro-Mesen  Fibro-mesenchymal alternative            7   
8         Others             Other / unresolved            8   

                                   definition_source  \
1  author-defined mapping from prior notebook; bo...   
2  author-defined mapping from prior notebook; bo...   
3  author-defined mapping from prior notebook; bo...   
4  author-defined mapping from prior notebook; bo...   
6  author-defined mapping from prior notebook; bo...   
7  author-defined mapping from prior notebook; bo...   
0  author-defined mapping from prior notebook; bo...   
5  author-defined mapping from prior notebook; bo...   
8  author-defined mapping from prior notebook; bo...   

   is_final_for_occupancy_panel  
1                          True  
2                          True  
3                          True  
4                          True  
6                          True  
7                          True  
0                          True  
5                          True  
8                          True  
```

### Text output (stream)

```text
Reviewed regulator-class counts

```

### Text output (text/plain)

```text
                       regulator_class  n_targets
0       Transcription / RNA processing         68
1  Metabolism / trafficking / adhesion         41
2           Translation / proteostasis         41
3                 Chromatin regulation         32
4      Apoptosis / cell cycle / stress         19
5            Other reviewed regulators         18
6                      TGF/BMP/Activin          9
7                          Hippo / YAP          6
8            WNT / mesoderm regulators          6
9              VEGF / EC specification          4
```

### Text output (stream)

```text
Guide types

```

### Text output (text/plain)

```text
             type  n_guides
0       targeting      1806
1   non-targeting        60
2  safe-targeting        60
```

## Markdown cell 8

## 3. Export compact cell metadata and target-union assignments

This step reads the H5AD once and exports only the metadata needed by the occupancy models.

Definitions encoded by the exporter:

- **Target carrier:** at least one assigned guide whose metadata maps to that `intended_target_name`.
- Multiple guides for the same target collapse to one target-union indicator; `n_guides_for_target` is retained.
- Cells with multiple biological targets contribute to each corresponding target indicator.
- **Pure control:** at least one safe/non-targeting guide and zero real-targeting guides.
- Cells carrying both control and real-targeting guides remain target carriers in the primary marginal model and are explicitly flagged.

## Code cell 9

```python
export_script = PACKAGE_DIR / "scripts" / "00_export_compact_occupancy_inputs.py"
if FORCE_REEXPORT or not (EXPORT_DIR / "export_audit.json").exists():
    cmd = [sys.executable, str(export_script), "--config", str(CONFIG_PATH)]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=PACKAGE_DIR)
else:
    print("Compact export already exists; set FORCE_REEXPORT=True to rebuild it.")

```

### Text output (stream)

```text
Running: /oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/scripts/00_export_compact_occupancy_inputs.py --config /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/config/perturbation_class_occupancy_config.yaml
Reading H5AD metadata/backed object: /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad
{
  "query_h5ad": "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad",
  "n_cells": 121616,
  "n_genes": 17538,
  "n_guides_h5ad": 1922,
  "n_guides_mapped": 1922,
  "n_targeting_guides_mapped": 1802,
  "n_control_guides_mapped": 120,
  "n_unique_targets_exported": 298,
  "n_cell_target_assignments": 613495,
  "n_pure_control_cells": 748,
  "n_mixed_control_target_cells": 20923,
  "state_source": "existing leiden_res_0.5 mapped with prior author-defined labels",
  "soft_state_probabilities_available": false
}
Wrote compact export: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/occupancy_compact_export.zip

```

## Code cell 10

```python
audit = json.loads((EXPORT_DIR / "export_audit.json").read_text())
print(json.dumps(audit, indent=2))

control_qc = pd.read_csv(EXPORT_DIR / "control_definition_qc.tsv", sep="\t")
state_counts = pd.read_csv(EXPORT_DIR / "state_counts_by_day_and_replicate.tsv", sep="\t")
carrier_counts = pd.read_csv(EXPORT_DIR / "target_carrier_counts_by_day.tsv", sep="\t")

print("Control-definition QC")
display(control_qc)
print("State counts by day and replicate")
display(state_counts.sort_values(["day", "biological_replicate", "assigned_state"]))
print("Carrier-count summary")
display(carrier_counts.groupby("day")["n_carrier_cells"].describe())

# Class coverage before modeling
exp_targets = set(carrier_counts["target_name"].astype(str))
class_cov = (
    class_membership.assign(in_experiment=class_membership["gene"].astype(str).isin(exp_targets))
    .groupby("regulator_class", as_index=False)
    .agg(n_reviewed_targets=("gene", "nunique"), n_targets_in_experiment=("in_experiment", "sum"))
    .sort_values("n_targets_in_experiment", ascending=False)
)
display(class_cov)

```

### Text output (stream)

```text
{
  "query_h5ad": "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad",
  "n_cells": 121616,
  "n_genes": 17538,
  "n_guides_h5ad": 1922,
  "n_guides_mapped": 1922,
  "n_targeting_guides_mapped": 1802,
  "n_control_guides_mapped": 120,
  "n_unique_targets_exported": 298,
  "n_cell_target_assignments": 613495,
  "n_pure_control_cells": 748,
  "n_mixed_control_target_cells": 20923,
  "state_source": "existing leiden_res_0.5 mapped with prior author-defined labels",
  "soft_state_probabilities_available": false
}
Control-definition QC

```

### Text output (text/plain)

```text
  day biological_replicate  n_cells  n_pure_controls  n_mixed_control_target  \
0  D0              d0_rep1    17488               75                    3026   
1  D0              d0_rep2    21121              127                    3547   
2  D1              d1_rep1     5856               27                    1178   
3  D1              d1_rep2    11069               48                    2031   
4  D2              d2_rep1    10895               41                    2284   
5  D2              d2_rep2    11085               52                    1929   
6  D3              d3_rep1    21377              185                    3282   
7  D3              d3_rep2    22725              193                    3646   

   n_any_targeting  
0            16263  
1            19322  
2             5505  
3            10382  
4            10295  
5            10182  
6            19018  
7            19641  
```

### Text output (stream)

```text
State counts by day and replicate

```

### Text output (text/plain)

```text
   day biological_replicate assigned_state  n_cells
0   D0              d0_rep1    Fibro-Mesen       10
1   D0              d0_rep1      Neu-Mesen        1
2   D0              d0_rep1     Neuro-iPSC     7921
3   D0              d0_rep1         Others       12
4   D0              d0_rep1     PrimStreak        1
5   D0              d0_rep1           iPSC     9543
6   D0              d0_rep2    Fibro-Mesen        7
7   D0              d0_rep2     Neuro-iPSC    10869
8   D0              d0_rep2         Others       54
9   D0              d0_rep2           iPSC    10191
10  D1              d1_rep1    Fibro-Mesen        1
11  D1              d1_rep1       Lat-Meso        1
12  D1              d1_rep1     Neuro-iPSC        1
13  D1              d1_rep1         Others       10
14  D1              d1_rep1     PrimStreak     5840
15  D1              d1_rep1           iPSC        3
16  D1              d1_rep2         Others       10
17  D1              d1_rep2     PrimStreak    11046
18  D1              d1_rep2           iPSC       13
19  D2              d2_rep1    Fibro-Mesen       20
20  D2              d2_rep1       Lat-Meso     7740
21  D2              d2_rep1      Neu-Mesen       15
22  D2              d2_rep1     Neuro-iPSC        4
23  D2              d2_rep1         Others        2
24  D2              d2_rep1         Pre-EC     3089
25  D2              d2_rep1     PrimStreak        3
26  D2              d2_rep1          artEC       12
27  D2              d2_rep1           iPSC       10
28  D2              d2_rep2    Fibro-Mesen       20
29  D2              d2_rep2       Lat-Meso     8352
30  D2              d2_rep2      Neu-Mesen       21
31  D2              d2_rep2     Neuro-iPSC        2
32  D2              d2_rep2         Others        2
33  D2              d2_rep2         Pre-EC     2676
34  D2              d2_rep2     PrimStreak        1
35  D2              d2_rep2          artEC        2
36  D2              d2_rep2           iPSC        9
37  D3              d3_rep1    Fibro-Mesen     5806
38  D3              d3_rep1       Lat-Meso       30
39  D3              d3_rep1      Neu-Mesen     1591
40  D3              d3_rep1     Neuro-iPSC        1
41  D3              d3_rep1         Others        3
42  D3              d3_rep1         Pre-EC       23
43  D3              d3_rep1     PrimStreak        1
44  D3              d3_rep1          artEC    13917
45  D3              d3_rep1           iPSC        5
46  D3              d3_rep2    Fibro-Mesen     6431
47  D3              d3_rep2       Lat-Meso       45
48  D3              d3_rep2      Neu-Mesen     1685
49  D3              d3_rep2     Neuro-iPSC        1
50  D3              d3_rep2         Pre-EC       27
51  D3              d3_rep2     PrimStreak        3
52  D3              d3_rep2          artEC    14527
53  D3              d3_rep2           iPSC        6
```

### Text output (stream)

```text
Carrier-count summary

```

### Text output (text/plain)

```text
     count        mean         std   min     25%    50%     75%     max
day                                                                    
D0   298.0  659.258389  178.914545  54.0  553.75  650.0  772.75  1659.0
D1   298.0  336.489933   79.994665  29.0  293.00  334.0  381.75   781.0
D2   298.0  430.291946  109.002264  44.0  367.00  428.5  493.00  1095.0
D3   298.0  632.667785  229.909220  39.0  486.25  614.0  751.00  1680.0
```

### Text output (text/plain)

```text
                       regulator_class  n_reviewed_targets  \
6       Transcription / RNA processing                  68   
3  Metabolism / trafficking / adhesion                  41   
7           Translation / proteostasis                  41   
1                 Chromatin regulation                  32   
0      Apoptosis / cell cycle / stress                  19   
4            Other reviewed regulators                  18   
5                      TGF/BMP/Activin                   9   
2                          Hippo / YAP                   6   
9            WNT / mesoderm regulators                   6   
8              VEGF / EC specification                   4   

   n_targets_in_experiment  
6                       68  
3                       41  
7                       41  
1                       32  
0                       19  
4                       18  
5                        9  
2                        6  
9                        6  
8                        4  
```

## Markdown cell 11

## 4. Fit per-target, per-day, per-state occupancy models

The primary population includes all cells from the relevant day. This is a marginal target-indicator model suited to the high-MOI design; it does not compare raw target fractions against a hand-selected control subset. Biological sample and the approved technical covariates are included in every model.

Benjamini–Hochberg adjustment is applied within each day across the target×state testing family. Replicate-specific models are fitted as a sensitivity check and combined by fixed-effect meta-analysis when estimable.

## Code cell 12

```python
fit_script = PACKAGE_DIR / "scripts" / "01_fit_occupancy_and_plot.py"
if FORCE_REFIT or not (RESULT_DIR / "target_state_occupancy_effects.tsv").exists():
    cmd = [
        sys.executable, str(fit_script),
        "--config", str(CONFIG_PATH),
        "--export-dir", str(EXPORT_DIR),
        "--output-dir", str(RESULT_DIR),
        "--n-jobs", str(N_JOBS),
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=PACKAGE_DIR)
else:
    print("Model results already exist; set FORCE_REFIT=True to rerun all models.")

```

### Text output (stream)

```text
Running: /oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/scripts/01_fit_occupancy_and_plot.py --config /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/config/perturbation_class_occupancy_config.yaml --export-dir /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/occupancy_compact_export --output-dir /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/occupancy_analysis_results --n-jobs 8

```

### Text output (stream)

```text
[Parallel(n_jobs=8)]: Using backend LokyBackend with 8 concurrent workers.
[Parallel(n_jobs=8)]: Done   2 tasks      | elapsed:    3.5s
[Parallel(n_jobs=8)]: Done   9 tasks      | elapsed:    5.2s
[Parallel(n_jobs=8)]: Done  16 tasks      | elapsed:    7.8s
[Parallel(n_jobs=8)]: Done  25 tasks      | elapsed:   10.6s
[Parallel(n_jobs=8)]: Done  34 tasks      | elapsed:   13.2s
[Parallel(n_jobs=8)]: Done  45 tasks      | elapsed:   17.0s
[Parallel(n_jobs=8)]: Done  56 tasks      | elapsed:   20.1s
[Parallel(n_jobs=8)]: Done  69 tasks      | elapsed:   24.5s
[Parallel(n_jobs=8)]: Done  82 tasks      | elapsed:   28.3s
[Parallel(n_jobs=8)]: Done  97 tasks      | elapsed:   33.2s
[Parallel(n_jobs=8)]: Done 112 tasks      | elapsed:   37.6s
[Parallel(n_jobs=8)]: Done 129 tasks      | elapsed:   42.7s
[Parallel(n_jobs=8)]: Done 146 tasks      | elapsed:   47.6s
[Parallel(n_jobs=8)]: Done 165 tasks      | elapsed:   53.9s
[Parallel(n_jobs=8)]: Done 184 tasks      | elapsed:  1.0min
[Parallel(n_jobs=8)]: Done 205 tasks      | elapsed:  1.1min
[Parallel(n_jobs=8)]: Done 226 tasks      | elapsed:  1.2min
[Parallel(n_jobs=8)]: Done 249 tasks      | elapsed:  1.4min
[Parallel(n_jobs=8)]: Done 272 tasks      | elapsed:  1.5min
[Parallel(n_jobs=8)]: Done 297 tasks      | elapsed:  1.6min
[Parallel(n_jobs=8)]: Done 322 tasks      | elapsed:  1.7min
[Parallel(n_jobs=8)]: Done 349 tasks      | elapsed:  1.9min
[Parallel(n_jobs=8)]: Done 376 tasks      | elapsed:  2.0min
[Parallel(n_jobs=8)]: Done 405 tasks      | elapsed:  2.2min
[Parallel(n_jobs=8)]: Done 434 tasks      | elapsed:  2.3min
[Parallel(n_jobs=8)]: Done 465 tasks      | elapsed:  2.5min
[Parallel(n_jobs=8)]: Done 496 tasks      | elapsed:  2.7min
[Parallel(n_jobs=8)]: Done 529 tasks      | elapsed:  2.8min
[Parallel(n_jobs=8)]: Done 562 tasks      | elapsed:  3.0min
[Parallel(n_jobs=8)]: Done 597 tasks      | elapsed:  3.2min
[Parallel(n_jobs=8)]: Done 632 tasks      | elapsed:  3.4min
[Parallel(n_jobs=8)]: Done 669 tasks      | elapsed:  3.6min
[Parallel(n_jobs=8)]: Done 706 tasks      | elapsed:  3.8min
[Parallel(n_jobs=8)]: Done 745 tasks      | elapsed:  4.0min
[Parallel(n_jobs=8)]: Done 784 tasks      | elapsed:  4.3min
[Parallel(n_jobs=8)]: Done 825 tasks      | elapsed:  4.5min
[Parallel(n_jobs=8)]: Done 866 tasks      | elapsed:  4.7min
[Parallel(n_jobs=8)]: Done 909 tasks      | elapsed:  5.0min
[Parallel(n_jobs=8)]: Done 952 tasks      | elapsed:  5.3min
[Parallel(n_jobs=8)]: Done 997 tasks      | elapsed:  5.6min
[Parallel(n_jobs=8)]: Done 1042 tasks      | elapsed:  5.8min
[Parallel(n_jobs=8)]: Done 1089 tasks      | elapsed:  6.2min
[Parallel(n_jobs=8)]: Done 1136 tasks      | elapsed:  6.5min
[Parallel(n_jobs=8)]: Done 1185 tasks      | elapsed:  6.8min
[Parallel(n_jobs=8)]: Done 1234 tasks      | elapsed:  7.0min
[Parallel(n_jobs=8)]: Done 1285 tasks      | elapsed:  7.2min
[Parallel(n_jobs=8)]: Done 1336 tasks      | elapsed:  7.5min
[Parallel(n_jobs=8)]: Done 1389 tasks      | elapsed:  7.7min
[Parallel(n_jobs=8)]: Done 1442 tasks      | elapsed:  7.9min
[Parallel(n_jobs=8)]: Done 1497 tasks      | elapsed:  8.2min
[Parallel(n_jobs=8)]: Done 1552 tasks      | elapsed:  8.4min
[Parallel(n_jobs=8)]: Done 1609 tasks      | elapsed:  8.6min
[Parallel(n_jobs=8)]: Done 1666 tasks      | elapsed:  8.9min
[Parallel(n_jobs=8)]: Done 1725 tasks      | elapsed:  9.1min
[Parallel(n_jobs=8)]: Done 1784 tasks      | elapsed:  9.4min
[Parallel(n_jobs=8)]: Done 1845 tasks      | elapsed:  9.7min
[Parallel(n_jobs=8)]: Done 1906 tasks      | elapsed: 10.0min
[Parallel(n_jobs=8)]: Done 1969 tasks      | elapsed: 10.3min
[Parallel(n_jobs=8)]: Done 2032 tasks      | elapsed: 10.5min
[Parallel(n_jobs=8)]: Done 2097 tasks      | elapsed: 10.8min
[Parallel(n_jobs=8)]: Done 2162 tasks      | elapsed: 11.1min
[Parallel(n_jobs=8)]: Done 2229 tasks      | elapsed: 11.4min
[Parallel(n_jobs=8)]: Done 2296 tasks      | elapsed: 11.7min
[Parallel(n_jobs=8)]: Done 2365 tasks      | elapsed: 12.0min
[Parallel(n_jobs=8)]: Done 2434 tasks      | elapsed: 12.4min
[Parallel(n_jobs=8)]: Done 2505 tasks      | elapsed: 12.7min
[Parallel(n_jobs=8)]: Done 2576 tasks      | elapsed: 13.0min
[Parallel(n_jobs=8)]: Done 2649 tasks      | elapsed: 13.4min
[Parallel(n_jobs=8)]: Done 2722 tasks      | elapsed: 13.7min
[Parallel(n_jobs=8)]: Done 2797 tasks      | elapsed: 14.1min
[Parallel(n_jobs=8)]: Done 2872 tasks      | elapsed: 14.5min
[Parallel(n_jobs=8)]: Done 2949 tasks      | elapsed: 14.8min
[Parallel(n_jobs=8)]: Done 3026 tasks      | elapsed: 15.2min
[Parallel(n_jobs=8)]: Done 3105 tasks      | elapsed: 15.6min
[Parallel(n_jobs=8)]: Done 3184 tasks      | elapsed: 16.0min
[Parallel(n_jobs=8)]: Done 3265 tasks      | elapsed: 16.4min
[Parallel(n_jobs=8)]: Done 3346 tasks      | elapsed: 16.8min
[Parallel(n_jobs=8)]: Done 3429 tasks      | elapsed: 17.2min
[Parallel(n_jobs=8)]: Done 3512 tasks      | elapsed: 17.6min
[Parallel(n_jobs=8)]: Done 3597 tasks      | elapsed: 17.9min
[Parallel(n_jobs=8)]: Done 3682 tasks      | elapsed: 18.3min
[Parallel(n_jobs=8)]: Done 3769 tasks      | elapsed: 18.8min
[Parallel(n_jobs=8)]: Done 3856 tasks      | elapsed: 19.2min
[Parallel(n_jobs=8)]: Done 3874 out of 3874 | elapsed: 19.3min finished

```

### Text output (stream)

```text
{
  "n_cells": 121616,
  "n_targets_with_assignments": 298,
  "n_tests": 3874,
  "n_classes_displayed": 10,
  "soft_sensitivity_run": false,
  "state_source": "existing author-defined Leiden mapping"
}

```

## Markdown cell 13

## 5. Audit target-level and class-level results

## Code cell 14

```python
effects = pd.read_csv(RESULT_DIR / "target_state_occupancy_effects.tsv", sep="\t")
rep_effects = pd.read_csv(RESULT_DIR / "target_state_occupancy_by_replicate.tsv", sep="\t")
class_summary = pd.read_csv(RESULT_DIR / "class_state_occupancy_summary.tsv", sep="\t")
coverage = pd.read_csv(RESULT_DIR / "coverage_and_filtering.tsv", sep="\t")

print("Target-level tests:", len(effects))
print("Converged:", int(effects["converged"].sum()), "/", len(effects))
print("FDR < 0.05:", int(effects["significant"].sum()))
print("Evaluable targets:", effects["target_name"].nunique())
print("Classes displayed:", class_summary["regulator_class"].nunique())

print("Coverage by day")
display(coverage.groupby("day").agg(
    n_targets=("target_name", "nunique"),
    n_passing=("passes_carrier_threshold", "sum"),
    median_carriers=("n_carrier_cells", "median"),
).reset_index())

print("Class-level summary")
display(class_summary.sort_values(["day", "state", "regulator_class"]))

if len(rep_effects):
    print("Replicate-specific models:", len(rep_effects))
    display(effects[[
        "day", "target_name", "state", "log_odds_ratio", "fdr",
        "replicate_sign_concordance", "n_replicates_evaluable"
    ]].sort_values("fdr").head(20))

```

### Text output (stream)

```text
Target-level tests: 3874
Converged: 3874 / 3874
FDR < 0.05: 390
Evaluable targets: 298
Classes displayed: 10
Coverage by day

```

### Text output (text/plain)

```text
  day  n_targets  n_passing  median_carriers
0  D0        298        298            650.0
1  D1        298        298            334.0
2  D2        298        298            428.5
3  D3        298        298            614.0
```

### Text output (stream)

```text
Class-level summary

```

### Text output (text/plain)

```text
                         regulator_class day       state  n_evaluable_targets  \
0        Apoptosis / cell cycle / stress  D0  Neuro-iPSC                   19   
13                  Chromatin regulation  D0  Neuro-iPSC                   32   
26                           Hippo / YAP  D0  Neuro-iPSC                    6   
39   Metabolism / trafficking / adhesion  D0  Neuro-iPSC                   41   
52             Other reviewed regulators  D0  Neuro-iPSC                   18   
..                                   ...  ..         ...                  ...   
77                       TGF/BMP/Activin  D3       artEC                    9   
90        Transcription / RNA processing  D3       artEC                   68   
103           Translation / proteostasis  D3       artEC                   41   
116              VEGF / EC specification  D3       artEC                    4   
129            WNT / mesoderm regulators  D3       artEC                    6   

     median_log_odds_ratio  fraction_significant_in_displayed_direction  \
0                 0.106824                                     0.052632   
13                0.063931                                     0.031250   
26                0.029332                                     0.000000   
39                0.042300                                     0.024390   
52                0.104258                                     0.111111   
..                     ...                                          ...   
77               -0.823696                                     1.000000   
90               -0.101523                                     0.279412   
103              -0.090386                                     0.170732   
116              -1.368988                                     0.750000   
129              -0.367633                                     0.500000   

     n_significant_in_displayed_direction displayed_direction  
0                                       1            enriched  
13                                      1            enriched  
26                                      0            enriched  
39                                      1            enriched  
52                                      2            enriched  
..                                    ...                 ...  
77                                      9            depleted  
90                                     19            depleted  
103                                     7            depleted  
116                                     3            depleted  
129                                     3            depleted  

[130 rows x 8 columns]
```

### Text output (stream)

```text
Replicate-specific models: 8344

```

### Text output (text/plain)

```text
     day target_name        state  log_odds_ratio            fdr  \
3516  D3       SMAD5    Neu-Mesen        2.831828   0.000000e+00   
3501  D3       SMAD1    Neu-Mesen        2.437874   0.000000e+00   
2491  D3      BMPR1A    Neu-Mesen        2.388320  9.045582e-146   
2996  D3       LEMD3    Neu-Mesen        1.404751   3.870738e-89   
3514  D3       SMAD5        artEC       -1.419610   1.964117e-77   
3499  D3       SMAD1        artEC       -1.178888   3.560420e-75   
3115  D3         NF2  Fibro-Mesen        1.360677   3.246243e-62   
3515  D3       SMAD5  Fibro-Mesen       -1.283492   4.627636e-53   
3189  D3       PLCG1        artEC       -1.734196   5.453826e-46   
3114  D3         NF2        artEC       -1.282843   9.631408e-46   
2401  D3      ACVR2A    Neu-Mesen        1.355627   1.015546e-44   
2489  D3      BMPR1A        artEC       -1.696886   2.402224e-44   
3716  D3        TBX3    Neu-Mesen        1.069766   7.398707e-42   
3500  D3       SMAD1  Fibro-Mesen       -0.966329   6.000008e-41   
3190  D3       PLCG1  Fibro-Mesen        1.250508   2.911762e-39   
3509  D3       SMAD4        artEC       -0.823696   4.209028e-37   
2525  D3       CASP3  Fibro-Mesen        0.786982   4.844630e-34   
2529  D3       CBLL1        artEC        0.940854   5.474247e-34   
2990  D3       LATS2  Fibro-Mesen        0.855180   1.566913e-33   
2964  D3         KDR        artEC       -1.412163   1.882574e-32   

      replicate_sign_concordance  n_replicates_evaluable  
3516                         1.0                       2  
3501                         1.0                       2  
2491                         1.0                       2  
2996                         1.0                       2  
3514                         1.0                       2  
3499                         1.0                       2  
3115                         1.0                       2  
3515                         1.0                       2  
3189                         1.0                       2  
3114                         1.0                       2  
2401                         1.0                       2  
2489                         1.0                       2  
3716                         1.0                       2  
3500                         1.0                       2  
3190                         1.0                       2  
3509                         1.0                       2  
2525                         1.0                       2  
2529                         1.0                       2  
2990                         1.0                       2  
2964                         1.0                       2  
```

## Markdown cell 15

## 6. Publication figures

## Code cell 16

```python
main_png = RESULT_DIR / "Fig2_perturbation_class_state_occupancy.png"
heatmap_png = RESULT_DIR / "Fig2_target_state_occupancy_heatmap_supplementary.png"

assert main_png.exists(), main_png
assert heatmap_png.exists(), heatmap_png

print("Main panel")
display(Image(filename=str(main_png)))
print("Supplementary target-level heatmap")
display(Image(filename=str(heatmap_png)))

```

### Text output (stream)

```text
Main panel

```

### Text output (text/plain)

```text
<IPython.core.display.Image object>
```

### Text output (stream)

```text
Supplementary target-level heatmap

```

### Text output (text/plain)

```text
<IPython.core.display.Image object>
```

## Markdown cell 17

## 7. Final output audit and packaging

The SVG exports retain text as text (`svg.fonttype = none`) and PDFs use TrueType fonts (`pdf.fonttype = 42`) for Illustrator compatibility.

## Code cell 18

```python
required_outputs = [
    "Fig2_perturbation_class_state_occupancy.svg",
    "Fig2_perturbation_class_state_occupancy.pdf",
    "Fig2_perturbation_class_state_occupancy.png",
    "Fig2_target_state_occupancy_heatmap_supplementary.svg",
    "Fig2_target_state_occupancy_heatmap_supplementary.pdf",
    "Fig2_target_state_occupancy_heatmap_supplementary.png",
    "target_state_occupancy_effects.tsv",
    "target_state_occupancy_by_replicate.tsv",
    "class_state_occupancy_summary.tsv",
    "class_membership.tsv",
    "state_definitions.tsv",
    "coverage_and_filtering.tsv",
    "draft_panel_legend.md",
    "methods_note.md",
]

missing = [name for name in required_outputs if not (RESULT_DIR / name).exists()]
print("Missing outputs:", missing)
assert not missing

# Check that SVG contains editable text nodes.
svg_text = (RESULT_DIR / "Fig2_perturbation_class_state_occupancy.svg").read_text(errors="ignore")
print("SVG contains <text> elements:", "<text" in svg_text)

zip_path = PACKAGE_DIR / "occupancy_analysis_results.zip"
if zip_path.exists():
    zip_path.unlink()
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for p in RESULT_DIR.rglob("*"):
        if p.is_file():
            zf.write(p, arcname=p.relative_to(RESULT_DIR))
print("Packaged results:", zip_path)

```

### Text output (stream)

```text
Missing outputs: []
SVG contains <text> elements: True
Packaged results: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_perturbation_class_occupancy_notebook_sherlock/occupancy_analysis_results.zip

```

## Markdown cell 19

## Interpretation and reporting caveats

- These are occupancy effects in pre-existing coarse states, not target→program effects.
- State boundaries are the existing author-defined Leiden labels and are not rediscovered from carriers.
- The primary coefficient is the adjusted target-level log odds ratio from a day-specific one-vs-rest logistic model.
- Class dots summarize fitted target-level estimates; class guides are not pooled into a single exposure.
- The legacy state labels are hard assignments, so a genuine soft-probability sensitivity analysis is not available from this object.
- The main figure should be described as an **author-defined coarse-state occupancy analysis**, not as a new control-defined clustering analysis.

