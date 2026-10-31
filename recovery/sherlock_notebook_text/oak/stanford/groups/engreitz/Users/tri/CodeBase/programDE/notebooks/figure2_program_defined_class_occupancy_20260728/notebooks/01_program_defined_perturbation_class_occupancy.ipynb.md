# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_program_defined_class_occupancy_20260728/notebooks/01_program_defined_perturbation_class_occupancy.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_program_defined_class_occupancy_20260728/notebooks/01_program_defined_perturbation_class_occupancy.ipynb`
- SHA-256: `bb4ae094f2f04a65131b82075b193c709a984a21252b9e88949355b07536dbd5`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# Program-defined perturbation-class occupancy analysis

This notebook defines cell states using the **global dominant cNMF program**, not Leiden clusters. It then asks whether biological targets redistribute cells among the curated program-defined states.

The workflow does **not** use target-to-program CRT effects or target-to-program differential-effect results.

Primary state rule:

1. Compute the global argmax across all 50 cell-level cNMF usages.
2. Do not force cells into the curated set.
3. Retain a cell in the primary analysis only when its global dominant program is in the curated primary table and the program's expected day matches the observed collection day.


## Code cell 1

```python
from pathlib import Path
import os, sys, subprocess, json, zipfile
import pandas as pd
import yaml
from IPython.display import display, Image, SVG

NOTEBOOK_PATH = Path.cwd()
PACKAGE_DIR = NOTEBOOK_PATH.parent if NOTEBOOK_PATH.name == "notebooks" else NOTEBOOK_PATH
PACKAGE_DIR = PACKAGE_DIR.resolve()
CONFIG = PACKAGE_DIR / "config" / "program_defined_class_occupancy_config.yaml"

print("PACKAGE_DIR:", PACKAGE_DIR)
print("CONFIG:", CONFIG)
assert CONFIG.exists(), CONFIG

```

### Text output (stream)

```text
PACKAGE_DIR: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_program_defined_class_occupancy_20260728
CONFIG: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_program_defined_class_occupancy_20260728/config/program_defined_class_occupancy_config.yaml

```

## Markdown cell 2

## 1. Inspect the curated program-state table

The primary list was transcribed from the first Hunter/Kyle annotation table. Programs flagged as uncertain or recommended for removal are retained only for sensitivity analysis, and programs marked for elimination are excluded.

## Code cell 3

```python
cfg = yaml.safe_load(CONFIG.read_text())
program_table_path = PACKAGE_DIR / cfg["curated_program_states_tsv"]
program_table = pd.read_csv(program_table_path, sep="\t")

display(program_table)
print("Primary programs:", program_table.query("use_in_primary == True")["program_id"].tolist())
print("Expanded sensitivity programs:", program_table.query("use_in_expanded_sensitivity == True")["program_id"].tolist())

```

### Text output (text/plain)

```text
    program_number program_id expected_day  \
0                7         P7           D0   
1               22        P22           D0   
2               29        P29           D0   
3               28        P28           D1   
4               13        P13           D1   
5               30        P30           D1   
6               48        P48           D1   
7               15        P15           D2   
8               14        P14           D2   
9               17        P17           D2   
10              33        P33           D3   
11              34        P34           D3   
12               2         P2           D3   
13               3         P3           D3   
14               8         P8           D3   
15              11        P11           D3   
16              23        P23           D2   
17              19        P19           D2   
18              27        P27           D3   
19              24        P24           D3   
20               6         P6          NaN   
21              21        P21          NaN   
22              35        P35          NaN   
23              42        P42          NaN   

                          program_name  \
0                         Pluripotency   
1                         Pluripotency   
2                         Pluripotency   
3                     Primitive streak   
4           Posterior primitive streak   
5    Primitive streak, BMP-repressed 1   
6    Primitive streak, BMP-repressed 2   
7                  Vascular progenitor   
8                     Lateral mesoderm   
9    Lateral mesoderm, YAP/TAZ-induced   
10              Anterior-like mesoderm   
11         Mesenchyme, YAP/TAZ-induced   
12                Artery endothelium 1   
13              Endothelial function 1   
14                Artery endothelium 2   
15              Endothelial function 2   
16           Vascular lateral mesoderm   
17                    Unknown mesoderm   
18     Unresolved alternative mesoderm   
19  Unresolved HOX/mesenchymal program   
20                          Eliminated   
21                          Eliminated   
22                          Eliminated   
23                          Eliminated   

                        representative_genes  use_in_primary  \
0                     L1TD1;DPPA4;TERF1;SOX2            True   
1                        POU5F1;TDGF1;DNMT3B            True   
2                           ESRG;DPPA4;NANOG            True   
3                             MIXL1;TBXT;SP5            True   
4                             NKD1;CDX2;MSX1            True   
5                             DKK1;NKD1;TBXT            True   
6                              NKD1;TBXT;WLS            True   
7                             TAL1;FLI1;SOX7            True   
8                                HAND1;HAND2            True   
9                          KRT8;KRT18;ANKRD1            True   
10                    CER1;LHX1;OTX2;CYP26A1            True   
11                         ANKRD1;KRT8;KRT18            True   
12                          MECOM;NRP1;EFNB2            True   
13                            ESM1;ESAM;CDH5            True   
14                 NOTCH1;SOX17;UNC5B;SEMA3G            True   
15                         EDN1;ANGPT2;PDGFB            True   
16                                  KDR;ETV2           False   
17                         WNT8A;WNT5A;PTCH1           False   
18  GATA4;RSPO2;SLIT3;NRP2;COL3A1;COL1A1;FN1           False   
19    HOXD10;LEFTY2;TGFBI;HOTAIR;HMGA2;MLLT3           False   
20                                       NaN           False   
21                                       NaN           False   
22                                       NaN           False   
23                                       NaN           False   

    use_in_expanded_sensitivity curation_status  \
0                          True         curated   
1                          True         curated   
2                          True         curated   
3                          True         curated   
4                          True         curated   
5                          True         curated   
6                          True         curated   
7                          True  reviewed named   
8                          True  reviewed named   
9                          True  reviewed named   
10                         True         curated   
11                         True         curated   
12                         True         curated   
13                         True         curated   
14                         True         curated   
15                         True         curated   
16                         True         flagged   
17                         True     provisional   
18                         True      unresolved   
19                         True      unresolved   
20                        False       eliminate   
21                        False       eliminate   
22                        False       eliminate   
23                        False       eliminate   

                                        curation_note  plot_order  
0                 Included in first-table curated set           1  
1                 Included in first-table curated set           2  
2                 Included in first-table curated set           3  
3                 Included in first-table curated set           4  
4                 Included in first-table curated set           5  
5                 Included in first-table curated set           6  
6                 Included in first-table curated set           7  
7   Named in first table and retained for primary ...           8  
8   Named in first table and retained for primary ...           9  
9   Named in first table and retained for primary ...          10  
10                Included in first-table curated set          11  
11                Included in first-table curated set          12  
12                Included in first-table curated set          13  
13                Included in first-table curated set          14  
14                Included in first-table curated set          15  
15                Included in first-table curated set          16  
16  First table recommends removal from Fig. 2G be...          17  
17              First table suggests possible removal          18  
18                First table requests further review          19  
19                 Identity unresolved in first table          20  
20                 Marked to eliminate in first table          90  
21                 Marked to eliminate in first table          91  
22                 Marked to eliminate in first table          92  
23                 Marked to eliminate in first table          93  
```

### Text output (stream)

```text
Primary programs: ['P7', 'P22', 'P29', 'P28', 'P13', 'P30', 'P48', 'P15', 'P14', 'P17', 'P33', 'P34', 'P2', 'P3', 'P8', 'P11']
Expanded sensitivity programs: ['P7', 'P22', 'P29', 'P28', 'P13', 'P30', 'P48', 'P15', 'P14', 'P17', 'P33', 'P34', 'P2', 'P3', 'P8', 'P11', 'P23', 'P19', 'P27', 'P24']

```

## Markdown cell 4

## 2. Verify the real CC-Perturb-seq H5AD

The path below was recovered from the earlier analysis notebooks. The object should contain D0-D3 metadata, guide assignments, and cell-level cNMF usages.

## Code cell 5

```python
query_path = Path(cfg["query_h5ad"])
print(query_path)
assert query_path.exists(), f"Query H5AD does not exist: {query_path}"
print("Size (GB):", round(query_path.stat().st_size / 1e9, 2))

usage_path = Path(cfg.get("program_usage_file", "")) if cfg.get("program_usage_file") else None
print("External cNMF usage table:", usage_path)
if usage_path is not None:
    assert usage_path.exists(), f"External cNMF usage file does not exist: {usage_path}"

# Preflight: the base H5AD may not persist cnmf_usage/Program columns because
# those were added in memory in ipsc-ec.Hippo.figure.ipynb.
import anndata as ad
_ad = ad.read_h5ad(query_path, backed="r")
program_obs_cols = [c for c in _ad.obs.columns if str(c).startswith("Program ")]
print("Program columns stored in H5AD:", len(program_obs_cols))
print("cnmf_usage stored in H5AD obsm:", "cnmf_usage" in _ad.obsm)
print("Exporter will use external table when neither is persisted.")
_ad.file.close()

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad
Size (GB): 4.22
External cNMF usage table: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50.usages.k_50.dt_0_2.consensus.txt
Program columns stored in H5AD: 0
cnmf_usage stored in H5AD obsm: False
Exporter will use external table when neither is persisted.

```

## Markdown cell 6

## 3. Export program-defined cell states and target-union assignments

This step reads the H5AD once, computes the global dominant program across all 50 programs, maps targeting guides to `intended_target_name`, and writes compact analysis tables.

## Code cell 7

```python
export_script = PACKAGE_DIR / "scripts" / "00_export_program_defined_occupancy_inputs.py"
cmd = [sys.executable, str(export_script), "--config", str(CONFIG)]
print(" ".join(cmd))
proc = subprocess.run(
    cmd, cwd=PACKAGE_DIR, text=True,
    stdout=subprocess.PIPE, stderr=subprocess.STDOUT
)
print(proc.stdout)
if proc.returncode != 0:
    raise RuntimeError(
        f"Exporter failed with return code {proc.returncode}. "
        "The complete child-process traceback is printed immediately above."
    )

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_program_defined_class_occupancy_20260728/scripts/00_export_program_defined_occupancy_inputs.py --config /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_program_defined_class_occupancy_20260728/config/program_defined_class_occupancy_config.yaml
Reading H5AD in backed mode: /oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad
Reading external cNMF usage table: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50.usages.k_50.dt_0_2.consensus.txt
cNMF usage overlap: 121,616/121,616 cells (100.0%)
{
  "query_h5ad": "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/2.preprocessing_tonys_pipeline/3.integrated_days/3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad",
  "n_cells": 121616,
  "n_genes": 17538,
  "program_usage_source": "external file: /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50.usages.k_50.dt_0_2.consensus.txt",
  "n_cells_with_program_usage": 121616,
  "fraction_cells_with_program_usage": 1.0,
  "state_assignment_rule": "global_argmax_then_primary_curated_expected_day",
  "n_primary_curated_programs": 16,
  "n_expanded_curated_programs": 20,
  "n_primary_program_state_cells": 69453,
  "n_target_names": 298,
  "control_definition": "at least one safe/non-targeting guide and zero targeting guides"
}


```

## Markdown cell 8

## 4. Audit program-state coverage

Important checks:

- fraction of each observed day represented by the primary curated program states;
- replicate representation;
- non-curated dominant programs that are excluded from the primary model;
- carrier coverage by day.


## Code cell 9

```python
cfg = yaml.safe_load(CONFIG.read_text())
export_dir = Path(cfg["export_dir"])
if not export_dir.is_absolute():
    export_dir = PACKAGE_DIR / export_dir

for fn in [
    "program_state_assignment_qc.tsv",
    "program_state_counts_by_day_and_replicate.tsv",
    "dominant_program_counts_by_day.tsv",
    "control_definition_qc.tsv",
    "target_carrier_counts_by_day.tsv",
]:
    path = export_dir / fn
    print("\n###", fn)
    display(pd.read_csv(path, sep="\t").head(50))

```

### Text output (stream)

```text

### program_state_assignment_qc.tsv

```

### Text output (text/plain)

```text
  day  n_all_cells  n_has_program_usage  n_primary_curated_dominant  \
0  D0        38609                38609                       16360   
1  D1        16925                16925                       10022   
2  D2        21980                21980                       11614   
3  D3        44102                44102                       31953   

   n_primary_program_state  n_expanded_program_state  n_pure_controls  \
0                    16350                     16350              202   
1                     9977                      9977               75   
2                    11594                     18010               93   
3                    31532                     38002              378   

   n_mixed_control_target  fraction_with_program_usage  \
0                    6573                          1.0   
1                    3209                          1.0   
2                    4213                          1.0   
3                    6928                          1.0   

   fraction_in_primary_program_state  
0                           0.423476  
1                           0.589483  
2                           0.527480  
3                           0.714979  
```

### Text output (stream)

```text

### program_state_counts_by_day_and_replicate.tsv

```

### Text output (text/plain)

```text
   day biological_replicate assigned_state              assigned_program_name  \
0   D0              d0_rep1            P22                       Pluripotency   
1   D0              d0_rep1            P29                       Pluripotency   
2   D0              d0_rep1             P7                       Pluripotency   
3   D0              d0_rep2            P22                       Pluripotency   
4   D0              d0_rep2            P29                       Pluripotency   
5   D0              d0_rep2             P7                       Pluripotency   
6   D1              d1_rep1            P13         Posterior primitive streak   
7   D1              d1_rep1            P28                   Primitive streak   
8   D1              d1_rep1            P30  Primitive streak, BMP-repressed 1   
9   D1              d1_rep1            P48  Primitive streak, BMP-repressed 2   
10  D1              d1_rep2            P13         Posterior primitive streak   
11  D1              d1_rep2            P28                   Primitive streak   
12  D1              d1_rep2            P30  Primitive streak, BMP-repressed 1   
13  D1              d1_rep2            P48  Primitive streak, BMP-repressed 2   
14  D2              d2_rep1            P14                   Lateral mesoderm   
15  D2              d2_rep1            P15                Vascular progenitor   
16  D2              d2_rep1            P17  Lateral mesoderm, YAP/TAZ-induced   
17  D2              d2_rep2            P14                   Lateral mesoderm   
18  D2              d2_rep2            P15                Vascular progenitor   
19  D2              d2_rep2            P17  Lateral mesoderm, YAP/TAZ-induced   
20  D3              d3_rep1            P11             Endothelial function 2   
21  D3              d3_rep1             P2               Artery endothelium 1   
22  D3              d3_rep1             P3             Endothelial function 1   
23  D3              d3_rep1            P33             Anterior-like mesoderm   
24  D3              d3_rep1            P34        Mesenchyme, YAP/TAZ-induced   
25  D3              d3_rep1             P8               Artery endothelium 2   
26  D3              d3_rep2            P11             Endothelial function 2   
27  D3              d3_rep2             P2               Artery endothelium 1   
28  D3              d3_rep2             P3             Endothelial function 1   
29  D3              d3_rep2            P33             Anterior-like mesoderm   
30  D3              d3_rep2            P34        Mesenchyme, YAP/TAZ-induced   
31  D3              d3_rep2             P8               Artery endothelium 2   

    n_cells  median_dominant_fraction  
0       208                  0.157566  
1       213                  0.169667  
2      8746                  0.211443  
3      3743                  0.181114  
4      1049                  0.158684  
5      2391                  0.213371  
6      1650                  0.204289  
7      2376                  0.194283  
8       667                  0.201416  
9        69                  0.199200  
10     2648                  0.228019  
11     1239                  0.181535  
12     1121                  0.222741  
13      207                  0.230355  
14     1656                  0.233971  
15     1836                  0.300382  
16     2014                  0.259524  
17     2468                  0.271814  
18     1740                  0.325332  
19     1880                  0.285227  
20     1047                  0.244876  
21     3833                  0.269916  
22     6760                  0.263218  
23     1284                  0.318333  
24      595                  0.232617  
25     1441                  0.257043  
26     1025                  0.258220  
27     4666                  0.293427  
28     5233                  0.268795  
29     1404                  0.321496  
30      591                  0.231658  
31     3653                  0.273528  
```

### Text output (stream)

```text

### dominant_program_counts_by_day.tsv

```

### Text output (text/plain)

```text
   day global_dominant_program_id       global_dominant_program_name  n_cells  \
0   D0                         P0                Non-curated program      316   
1   D0                         P1                Non-curated program       53   
2   D0                        P10                Non-curated program      158   
3   D0                        P12                Non-curated program       34   
4   D0                        P16                Non-curated program        1   
5   D0                        P18                Non-curated program        7   
6   D0                        P20                Non-curated program        2   
7   D0                        P21                         Eliminated     1719   
8   D0                        P22                       Pluripotency     3951   
9   D0                        P25                Non-curated program     2873   
10  D0                        P26                Non-curated program     2212   
11  D0                        P29                       Pluripotency     1262   
12  D0                        P31                Non-curated program     2001   
13  D0                        P32                Non-curated program       26   
14  D0                        P34        Mesenchyme, YAP/TAZ-induced       10   
15  D0                        P35                         Eliminated        1   
16  D0                        P36                Non-curated program      279   
17  D0                        P37                Non-curated program     1627   
18  D0                        P39                Non-curated program        1   
19  D0                         P4                Non-curated program     1653   
20  D0                        P41                Non-curated program       14   
21  D0                        P44                Non-curated program        3   
22  D0                        P45                Non-curated program      267   
23  D0                        P49                Non-curated program      310   
24  D0                         P5                Non-curated program      601   
25  D0                         P6                         Eliminated     7938   
26  D0                         P7                       Pluripotency    11137   
27  D0                         P9                Non-curated program      153   
28  D1                         P0                Non-curated program      358   
29  D1                         P1                Non-curated program      291   
30  D1                        P10                Non-curated program      139   
31  D1                        P12                Non-curated program       19   
32  D1                        P13         Posterior primitive streak     4298   
33  D1                        P18                Non-curated program     2016   
34  D1                        P20                Non-curated program     2040   
35  D1                        P26                Non-curated program        1   
36  D1                        P28                   Primitive streak     3615   
37  D1                        P30  Primitive streak, BMP-repressed 1     1788   
38  D1                        P32                Non-curated program       12   
39  D1                        P36                Non-curated program        2   
40  D1                         P4                Non-curated program     1303   
41  D1                        P41                Non-curated program        1   
42  D1                        P48  Primitive streak, BMP-repressed 2      276   
43  D1                         P5                Non-curated program      509   
44  D1                         P6                         Eliminated        1   
45  D1                         P7                       Pluripotency       45   
46  D1                         P9                Non-curated program      211   
47  D2                         P0                Non-curated program     1506   
48  D2                         P1                Non-curated program     1647   
49  D2                        P10                Non-curated program       44   

    median_dominant_fraction  
0                   0.154922  
1                   0.183861  
2                   0.133995  
3                   0.136656  
4                   0.366218  
5                   0.114846  
6                   0.147813  
7                   0.174639  
8                   0.180207  
9                   0.276729  
10                  0.190517  
11                  0.160120  
12                  0.200294  
13                  0.138649  
14                  0.218978  
15                  0.155481  
16                  0.222305  
17                  0.225033  
18                  0.158640  
19                  0.156059  
20                  0.165429  
21                  0.185508  
22                  0.194029  
23                  0.199728  
24                  0.146898  
25                  0.231542  
26                  0.211941  
27                  0.155345  
28                  0.174709  
29                  0.228679  
30                  0.161716  
31                  0.138861  
32                  0.218590  
33                  0.235742  
34                  0.242629  
35                  0.129220  
36                  0.189809  
37                  0.214831  
38                  0.165727  
39                  0.221506  
40                  0.175379  
41                  0.271709  
42                  0.218600  
43                  0.162851  
44                  0.178031  
45                  0.141485  
46                  0.169517  
47                  0.194671  
48                  0.265389  
49                  0.151590  
```

### Text output (stream)

```text

### control_definition_qc.tsv

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

### target_carrier_counts_by_day.tsv

```

### Text output (text/plain)

```text
   target_name day  n_carrier_cells  \
0        ACAA1  D0              591   
1        ACAA1  D1              336   
2        ACAA1  D2              405   
3        ACAA1  D3              634   
4       ACTL7B  D0              597   
5       ACTL7B  D1              302   
6       ACTL7B  D2              407   
7       ACTL7B  D3              638   
8       ACVR1B  D0              650   
9       ACVR1B  D1              293   
10      ACVR1B  D2              390   
11      ACVR1B  D3              639   
12      ACVR2A  D0              707   
13      ACVR2A  D1              320   
14      ACVR2A  D2              400   
15      ACVR2A  D3              585   
16      ACVR2B  D0              810   
17      ACVR2B  D1              386   
18      ACVR2B  D2              489   
19      ACVR2B  D3              621   
20       ADCY2  D0              845   
21       ADCY2  D1              422   
22       ADCY2  D2              603   
23       ADCY2  D3              865   
24      AMOTL2  D0              740   
25      AMOTL2  D1              328   
26      AMOTL2  D2              503   
27      AMOTL2  D3              810   
28      ANP32E  D0              499   
29      ANP32E  D1              272   
30      ANP32E  D2              322   
31      ANP32E  D3              552   
32       APAF1  D0              517   
33       APAF1  D1              233   
34       APAF1  D2              372   
35       APAF1  D3              837   
36       ARIH1  D0              500   
37       ARIH1  D1              267   
38       ARIH1  D2              357   
39       ARIH1  D3              436   
40        ARTN  D0              843   
41        ARTN  D1              435   
42        ARTN  D2              527   
43        ARTN  D3              726   
44       ASCC3  D0              769   
45       ASCC3  D1              389   
46       ASCC3  D2              512   
47       ASCC3  D3              533   
48       ASH2L  D0              579   
49       ASH2L  D1              338   

    n_carrier_cells_in_primary_program_states  median_guides_for_target  
0                                         271                       1.0  
1                                         173                       1.0  
2                                         204                       1.0  
3                                         411                       1.0  
4                                         252                       1.0  
5                                         163                       1.0  
6                                         190                       1.0  
7                                         381                       1.0  
8                                         283                       1.0  
9                                         175                       1.0  
10                                        219                       1.0  
11                                        317                       1.0  
12                                        309                       1.0  
13                                        160                       1.0  
14                                        101                       1.0  
15                                        323                       1.0  
16                                        349                       1.0  
17                                        207                       1.0  
18                                        199                       1.0  
19                                        303                       1.0  
20                                        387                       1.0  
21                                        236                       1.0  
22                                        294                       1.0  
23                                        533                       1.0  
24                                        348                       1.0  
25                                        178                       1.0  
26                                        262                       1.0  
27                                        446                       1.0  
28                                        232                       1.0  
29                                        143                       1.0  
30                                        171                       1.0  
31                                        372                       1.0  
32                                        242                       1.0  
33                                        118                       1.0  
34                                        194                       1.0  
35                                        389                       1.0  
36                                        226                       1.0  
37                                        144                       1.0  
38                                        162                       1.0  
39                                        212                       1.0  
40                                        386                       1.0  
41                                        226                       1.0  
42                                        265                       1.0  
43                                        466                       1.0  
44                                        359                       1.0  
45                                        199                       1.0  
46                                        256                       1.0  
47                                        273                       1.0  
48                                        270                       1.0  
49                                        186                       1.0  
```

## Markdown cell 10

## 5. Fit target-level occupancy models and create the class summary

For each observed day, target, and curated dominant-program state, the script fits a one-versus-rest logistic model with the established technical covariates. Target classes are summarized only after the target-level models are fitted.

## Code cell 11

```python
analysis_script = PACKAGE_DIR / "scripts" / "01_fit_program_defined_occupancy_and_plot.py"
cmd = [sys.executable, str(analysis_script), "--config", str(CONFIG)]
print(" ".join(cmd))
subprocess.run(cmd, check=True, cwd=PACKAGE_DIR)

```

### Text output (stream)

```text
/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_program_defined_class_occupancy_20260728/scripts/01_fit_program_defined_occupancy_and_plot.py --config /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_program_defined_class_occupancy_20260728/config/program_defined_class_occupancy_config.yaml

```

### Text output (stream)

```text
[Parallel(n_jobs=8)]: Using backend LokyBackend with 8 concurrent workers.
[Parallel(n_jobs=8)]: Done   2 tasks      | elapsed:    3.1s
[Parallel(n_jobs=8)]: Done   9 tasks      | elapsed:    3.7s
[Parallel(n_jobs=8)]: Done  16 tasks      | elapsed:    5.0s
[Parallel(n_jobs=8)]: Done  25 tasks      | elapsed:    6.7s
[Parallel(n_jobs=8)]: Done  34 tasks      | elapsed:    8.7s
[Parallel(n_jobs=8)]: Done  45 tasks      | elapsed:   10.8s
[Parallel(n_jobs=8)]: Done  56 tasks      | elapsed:   13.0s
[Parallel(n_jobs=8)]: Done  69 tasks      | elapsed:   15.9s
[Parallel(n_jobs=8)]: Done  82 tasks      | elapsed:   18.3s
[Parallel(n_jobs=8)]: Done  97 tasks      | elapsed:   21.4s
[Parallel(n_jobs=8)]: Done 112 tasks      | elapsed:   24.3s
[Parallel(n_jobs=8)]: Done 129 tasks      | elapsed:   27.8s
[Parallel(n_jobs=8)]: Done 146 tasks      | elapsed:   31.1s
[Parallel(n_jobs=8)]: Done 165 tasks      | elapsed:   35.1s
[Parallel(n_jobs=8)]: Done 184 tasks      | elapsed:   38.8s
[Parallel(n_jobs=8)]: Done 205 tasks      | elapsed:   43.2s
[Parallel(n_jobs=8)]: Done 226 tasks      | elapsed:   47.5s
[Parallel(n_jobs=8)]: Done 249 tasks      | elapsed:   51.9s
[Parallel(n_jobs=8)]: Done 272 tasks      | elapsed:   56.4s
[Parallel(n_jobs=8)]: Done 297 tasks      | elapsed:  1.0min
[Parallel(n_jobs=8)]: Done 322 tasks      | elapsed:  1.1min
[Parallel(n_jobs=8)]: Done 349 tasks      | elapsed:  1.2min
[Parallel(n_jobs=8)]: Done 376 tasks      | elapsed:  1.3min
[Parallel(n_jobs=8)]: Done 405 tasks      | elapsed:  1.4min
[Parallel(n_jobs=8)]: Done 434 tasks      | elapsed:  1.5min
[Parallel(n_jobs=8)]: Done 465 tasks      | elapsed:  1.6min
[Parallel(n_jobs=8)]: Done 496 tasks      | elapsed:  1.7min
[Parallel(n_jobs=8)]: Done 529 tasks      | elapsed:  1.8min
[Parallel(n_jobs=8)]: Done 562 tasks      | elapsed:  1.9min
[Parallel(n_jobs=8)]: Done 597 tasks      | elapsed:  2.0min
[Parallel(n_jobs=8)]: Done 632 tasks      | elapsed:  2.1min
[Parallel(n_jobs=8)]: Done 669 tasks      | elapsed:  2.3min
[Parallel(n_jobs=8)]: Done 706 tasks      | elapsed:  2.4min
[Parallel(n_jobs=8)]: Done 745 tasks      | elapsed:  2.5min
[Parallel(n_jobs=8)]: Done 784 tasks      | elapsed:  2.7min
[Parallel(n_jobs=8)]: Done 825 tasks      | elapsed:  2.8min
[Parallel(n_jobs=8)]: Done 866 tasks      | elapsed:  2.9min
[Parallel(n_jobs=8)]: Done 909 tasks      | elapsed:  3.1min
[Parallel(n_jobs=8)]: Done 952 tasks      | elapsed:  3.1min
[Parallel(n_jobs=8)]: Done 997 tasks      | elapsed:  3.2min
[Parallel(n_jobs=8)]: Done 1042 tasks      | elapsed:  3.3min
[Parallel(n_jobs=8)]: Done 1089 tasks      | elapsed:  3.4min
[Parallel(n_jobs=8)]: Done 1136 tasks      | elapsed:  3.5min
[Parallel(n_jobs=8)]: Done 1185 tasks      | elapsed:  3.6min
[Parallel(n_jobs=8)]: Done 1234 tasks      | elapsed:  3.8min
[Parallel(n_jobs=8)]: Done 1285 tasks      | elapsed:  3.9min
[Parallel(n_jobs=8)]: Done 1336 tasks      | elapsed:  4.0min
[Parallel(n_jobs=8)]: Done 1389 tasks      | elapsed:  4.1min
[Parallel(n_jobs=8)]: Done 1442 tasks      | elapsed:  4.2min
[Parallel(n_jobs=8)]: Done 1497 tasks      | elapsed:  4.3min
[Parallel(n_jobs=8)]: Done 1552 tasks      | elapsed:  4.4min
[Parallel(n_jobs=8)]: Done 1609 tasks      | elapsed:  4.5min
[Parallel(n_jobs=8)]: Done 1666 tasks      | elapsed:  4.6min
[Parallel(n_jobs=8)]: Done 1725 tasks      | elapsed:  4.8min
[Parallel(n_jobs=8)]: Done 1784 tasks      | elapsed:  4.9min
[Parallel(n_jobs=8)]: Done 1845 tasks      | elapsed:  5.0min
[Parallel(n_jobs=8)]: Done 1906 tasks      | elapsed:  5.2min
[Parallel(n_jobs=8)]: Done 1969 tasks      | elapsed:  5.3min
[Parallel(n_jobs=8)]: Done 2032 tasks      | elapsed:  5.4min
[Parallel(n_jobs=8)]: Done 2097 tasks      | elapsed:  5.6min
[Parallel(n_jobs=8)]: Done 2162 tasks      | elapsed:  5.8min
[Parallel(n_jobs=8)]: Done 2229 tasks      | elapsed:  6.0min
[Parallel(n_jobs=8)]: Done 2296 tasks      | elapsed:  6.2min
[Parallel(n_jobs=8)]: Done 2365 tasks      | elapsed:  6.4min
[Parallel(n_jobs=8)]: Done 2434 tasks      | elapsed:  6.6min
[Parallel(n_jobs=8)]: Done 2505 tasks      | elapsed:  6.9min
[Parallel(n_jobs=8)]: Done 2576 tasks      | elapsed:  7.1min
[Parallel(n_jobs=8)]: Done 2649 tasks      | elapsed:  7.3min
[Parallel(n_jobs=8)]: Done 2722 tasks      | elapsed:  7.6min
[Parallel(n_jobs=8)]: Done 2797 tasks      | elapsed:  7.8min
[Parallel(n_jobs=8)]: Done 2872 tasks      | elapsed:  8.0min
[Parallel(n_jobs=8)]: Done 2949 tasks      | elapsed:  8.3min
[Parallel(n_jobs=8)]: Done 3026 tasks      | elapsed:  8.6min
[Parallel(n_jobs=8)]: Done 3105 tasks      | elapsed:  8.9min
[Parallel(n_jobs=8)]: Done 3184 tasks      | elapsed:  9.2min
[Parallel(n_jobs=8)]: Done 3265 tasks      | elapsed:  9.4min
[Parallel(n_jobs=8)]: Done 3346 tasks      | elapsed:  9.7min
[Parallel(n_jobs=8)]: Done 3429 tasks      | elapsed: 10.0min
[Parallel(n_jobs=8)]: Done 3512 tasks      | elapsed: 10.4min
[Parallel(n_jobs=8)]: Done 3597 tasks      | elapsed: 10.7min
[Parallel(n_jobs=8)]: Done 3682 tasks      | elapsed: 11.0min
[Parallel(n_jobs=8)]: Done 3769 tasks      | elapsed: 11.3min
[Parallel(n_jobs=8)]: Done 3856 tasks      | elapsed: 11.6min
[Parallel(n_jobs=8)]: Done 3945 tasks      | elapsed: 11.9min
[Parallel(n_jobs=8)]: Done 4034 tasks      | elapsed: 12.3min
[Parallel(n_jobs=8)]: Done 4125 tasks      | elapsed: 12.6min
[Parallel(n_jobs=8)]: Done 4216 tasks      | elapsed: 12.9min
[Parallel(n_jobs=8)]: Done 4309 tasks      | elapsed: 13.3min
[Parallel(n_jobs=8)]: Done 4402 tasks      | elapsed: 13.6min
[Parallel(n_jobs=8)]: Done 4497 tasks      | elapsed: 14.0min
[Parallel(n_jobs=8)]: Done 4592 tasks      | elapsed: 14.3min
[Parallel(n_jobs=8)]: Done 4689 tasks      | elapsed: 14.7min
[Parallel(n_jobs=8)]: Done 4752 out of 4752 | elapsed: 14.9min finished
[Parallel(n_jobs=8)]: Using backend LokyBackend with 8 concurrent workers.
[Parallel(n_jobs=8)]: Done   2 tasks      | elapsed:    0.7s
[Parallel(n_jobs=8)]: Done   9 tasks      | elapsed:    2.7s
[Parallel(n_jobs=8)]: Done  16 tasks      | elapsed:    4.2s
[Parallel(n_jobs=8)]: Done  25 tasks      | elapsed:    6.5s
[Parallel(n_jobs=8)]: Done  34 tasks      | elapsed:    8.7s
[Parallel(n_jobs=8)]: Done  45 tasks      | elapsed:   11.6s
[Parallel(n_jobs=8)]: Done  56 tasks      | elapsed:   14.4s
[Parallel(n_jobs=8)]: Done  69 tasks      | elapsed:   17.8s
[Parallel(n_jobs=8)]: Done  82 tasks      | elapsed:   21.2s
[Parallel(n_jobs=8)]: Done  97 tasks      | elapsed:   24.9s
[Parallel(n_jobs=8)]: Done 112 tasks      | elapsed:   28.7s
[Parallel(n_jobs=8)]: Done 129 tasks      | elapsed:   33.2s
[Parallel(n_jobs=8)]: Done 146 tasks      | elapsed:   37.3s
[Parallel(n_jobs=8)]: Done 165 tasks      | elapsed:   42.0s
[Parallel(n_jobs=8)]: Done 184 tasks      | elapsed:   46.5s
[Parallel(n_jobs=8)]: Done 205 tasks      | elapsed:   51.8s
[Parallel(n_jobs=8)]: Done 226 tasks      | elapsed:   56.8s
[Parallel(n_jobs=8)]: Done 249 tasks      | elapsed:  1.0min
[Parallel(n_jobs=8)]: Done 272 tasks      | elapsed:  1.1min
[Parallel(n_jobs=8)]: Done 297 tasks      | elapsed:  1.3min
[Parallel(n_jobs=8)]: Done 322 tasks      | elapsed:  1.4min
[Parallel(n_jobs=8)]: Done 349 tasks      | elapsed:  1.5min
[Parallel(n_jobs=8)]: Done 376 tasks      | elapsed:  1.6min
[Parallel(n_jobs=8)]: Done 405 tasks      | elapsed:  1.7min
[Parallel(n_jobs=8)]: Done 434 tasks      | elapsed:  1.8min
[Parallel(n_jobs=8)]: Done 465 tasks      | elapsed:  2.0min
[Parallel(n_jobs=8)]: Done 496 tasks      | elapsed:  2.1min
[Parallel(n_jobs=8)]: Done 529 tasks      | elapsed:  2.2min
[Parallel(n_jobs=8)]: Done 562 tasks      | elapsed:  2.4min
[Parallel(n_jobs=8)]: Done 597 tasks      | elapsed:  2.5min
[Parallel(n_jobs=8)]: Done 632 tasks      | elapsed:  2.7min
[Parallel(n_jobs=8)]: Done 669 tasks      | elapsed:  2.8min
[Parallel(n_jobs=8)]: Done 706 tasks      | elapsed:  3.0min
[Parallel(n_jobs=8)]: Done 745 tasks      | elapsed:  3.1min
[Parallel(n_jobs=8)]: Done 784 tasks      | elapsed:  3.3min
[Parallel(n_jobs=8)]: Done 825 tasks      | elapsed:  3.5min
[Parallel(n_jobs=8)]: Done 866 tasks      | elapsed:  3.7min
[Parallel(n_jobs=8)]: Done 909 tasks      | elapsed:  3.8min
[Parallel(n_jobs=8)]: Done 952 tasks      | elapsed:  4.0min
[Parallel(n_jobs=8)]: Done 997 tasks      | elapsed:  4.1min
[Parallel(n_jobs=8)]: Done 1042 tasks      | elapsed:  4.3min
[Parallel(n_jobs=8)]: Done 1089 tasks      | elapsed:  4.4min
[Parallel(n_jobs=8)]: Done 1136 tasks      | elapsed:  4.6min
[Parallel(n_jobs=8)]: Done 1185 tasks      | elapsed:  4.8min
[Parallel(n_jobs=8)]: Done 1234 tasks      | elapsed:  5.0min
[Parallel(n_jobs=8)]: Done 1285 tasks      | elapsed:  5.1min
[Parallel(n_jobs=8)]: Done 1336 tasks      | elapsed:  5.3min
[Parallel(n_jobs=8)]: Done 1389 tasks      | elapsed:  5.5min
[Parallel(n_jobs=8)]: Done 1442 tasks      | elapsed:  5.7min
[Parallel(n_jobs=8)]: Done 1497 tasks      | elapsed:  5.9min
[Parallel(n_jobs=8)]: Done 1552 tasks      | elapsed:  6.1min
[Parallel(n_jobs=8)]: Done 1609 tasks      | elapsed:  6.3min
[Parallel(n_jobs=8)]: Done 1666 tasks      | elapsed:  6.5min
[Parallel(n_jobs=8)]: Done 1725 tasks      | elapsed:  6.7min
[Parallel(n_jobs=8)]: Done 1784 tasks      | elapsed:  6.9min
[Parallel(n_jobs=8)]: Done 1845 tasks      | elapsed:  7.1min
[Parallel(n_jobs=8)]: Done 1906 tasks      | elapsed:  7.3min
[Parallel(n_jobs=8)]: Done 1969 tasks      | elapsed:  7.5min
[Parallel(n_jobs=8)]: Done 2032 tasks      | elapsed:  7.8min
[Parallel(n_jobs=8)]: Done 2097 tasks      | elapsed:  8.0min
[Parallel(n_jobs=8)]: Done 2162 tasks      | elapsed:  8.2min
[Parallel(n_jobs=8)]: Done 2229 tasks      | elapsed:  8.4min
[Parallel(n_jobs=8)]: Done 2296 tasks      | elapsed:  8.7min
[Parallel(n_jobs=8)]: Done 2365 tasks      | elapsed:  8.9min
[Parallel(n_jobs=8)]: Done 2434 tasks      | elapsed:  9.1min
[Parallel(n_jobs=8)]: Done 2505 tasks      | elapsed:  9.4min
[Parallel(n_jobs=8)]: Done 2576 tasks      | elapsed:  9.6min
[Parallel(n_jobs=8)]: Done 2649 tasks      | elapsed:  9.9min
[Parallel(n_jobs=8)]: Done 2722 tasks      | elapsed: 10.2min
[Parallel(n_jobs=8)]: Done 2797 tasks      | elapsed: 10.4min
[Parallel(n_jobs=8)]: Done 2872 tasks      | elapsed: 10.7min
[Parallel(n_jobs=8)]: Done 2949 tasks      | elapsed: 11.0min
[Parallel(n_jobs=8)]: Done 3026 tasks      | elapsed: 11.3min
[Parallel(n_jobs=8)]: Done 3105 tasks      | elapsed: 11.6min
[Parallel(n_jobs=8)]: Done 3184 tasks      | elapsed: 11.9min
[Parallel(n_jobs=8)]: Done 3265 tasks      | elapsed: 12.3min
[Parallel(n_jobs=8)]: Done 3346 tasks      | elapsed: 12.6min
[Parallel(n_jobs=8)]: Done 3429 tasks      | elapsed: 13.0min
[Parallel(n_jobs=8)]: Done 3512 tasks      | elapsed: 13.3min
[Parallel(n_jobs=8)]: Done 3597 tasks      | elapsed: 13.7min
[Parallel(n_jobs=8)]: Done 3682 tasks      | elapsed: 14.0min
[Parallel(n_jobs=8)]: Done 3769 tasks      | elapsed: 14.4min
[Parallel(n_jobs=8)]: Done 3856 tasks      | elapsed: 14.8min
[Parallel(n_jobs=8)]: Done 3945 tasks      | elapsed: 15.1min
[Parallel(n_jobs=8)]: Done 4034 tasks      | elapsed: 15.5min
[Parallel(n_jobs=8)]: Done 4125 tasks      | elapsed: 15.9min
[Parallel(n_jobs=8)]: Done 4216 tasks      | elapsed: 16.3min
[Parallel(n_jobs=8)]: Done 4309 tasks      | elapsed: 16.7min
[Parallel(n_jobs=8)]: Done 4402 tasks      | elapsed: 17.1min
[Parallel(n_jobs=8)]: Done 4497 tasks      | elapsed: 17.5min
[Parallel(n_jobs=8)]: Done 4592 tasks      | elapsed: 17.9min
[Parallel(n_jobs=8)]: Done 4689 tasks      | elapsed: 18.3min
[Parallel(n_jobs=8)]: Done 4768 out of 4768 | elapsed: 18.6min finished

```

### Text output (stream)

```text
{
  "n_all_cells": 121616,
  "n_primary_program_state_cells": 69453,
  "n_targets_with_assignments": 298,
  "n_target_program_state_tests": 4752,
  "n_classes_displayed": 10,
  "state_definition": "global dominant cNMF program; primary curated program list; expected-day match",
  "target_to_program_differential_effects_used": false
}

```

### Text output (text/plain)

```text
CompletedProcess(args=['/oak/stanford/groups/engreitz/Users/tri/envs/sc-dl/bin/python', '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_program_defined_class_occupancy_20260728/scripts/01_fit_program_defined_occupancy_and_plot.py', '--config', '/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_program_defined_class_occupancy_20260728/config/program_defined_class_occupancy_config.yaml'], returncode=0)
```

## Markdown cell 12

## 6. Inspect the main results and coverage

## Code cell 13

```python
cfg = yaml.safe_load(CONFIG.read_text())
result_dir = Path(cfg["output_dir"])
if not result_dir.is_absolute():
    result_dir = PACKAGE_DIR / result_dir

for fn in [
    "class_program_state_occupancy_summary.tsv",
    "class_coverage_summary.tsv",
    "coverage_and_filtering.tsv",
    "target_class_program_state_raw_distribution.tsv",
    "analysis_audit.json",
]:
    path = result_dir / fn
    print("\n###", fn)
    if path.suffix == ".json":
        print(path.read_text())
    else:
        display(pd.read_csv(path, sep="\t").head(80))

```

### Text output (stream)

```text

### class_program_state_occupancy_summary.tsv

```

### Text output (text/plain)

```text
                      analysis_population                  regulator_class  \
0   primary_curated_dominant_expected_day  Apoptosis / cell cycle / stress   
1   primary_curated_dominant_expected_day  Apoptosis / cell cycle / stress   
2   primary_curated_dominant_expected_day  Apoptosis / cell cycle / stress   
3   primary_curated_dominant_expected_day  Apoptosis / cell cycle / stress   
4   primary_curated_dominant_expected_day  Apoptosis / cell cycle / stress   
..                                    ...                              ...   
75  primary_curated_dominant_expected_day        Other reviewed regulators   
76  primary_curated_dominant_expected_day        Other reviewed regulators   
77  primary_curated_dominant_expected_day        Other reviewed regulators   
78  primary_curated_dominant_expected_day        Other reviewed regulators   
79  primary_curated_dominant_expected_day        Other reviewed regulators   

   day state  n_evaluable_targets  median_log_odds_ratio  \
0   D0   P22                   19              -0.086249   
1   D0   P29                   19               0.028386   
2   D0    P7                   19               0.067475   
3   D1   P13                   19               0.094020   
4   D1   P28                   19              -0.099223   
..  ..   ...                  ...                    ...   
75  D3    P2                   18              -0.094316   
76  D3    P3                   18              -0.057749   
77  D3   P33                   18              -0.081740   
78  D3   P34                   18               0.009345   
79  D3    P8                   18              -0.040751   

    fraction_significant_in_displayed_direction  \
0                                      0.000000   
1                                      0.000000   
2                                      0.000000   
3                                      0.000000   
4                                      0.000000   
..                                          ...   
75                                     0.000000   
76                                     0.166667   
77                                     0.000000   
78                                     0.111111   
79                                     0.000000   

    n_significant_in_displayed_direction displayed_direction  
0                                      0            depleted  
1                                      0            enriched  
2                                      0            enriched  
3                                      0            enriched  
4                                      0            depleted  
..                                   ...                 ...  
75                                     0            depleted  
76                                     3            depleted  
77                                     0            depleted  
78                                     2            enriched  
79                                     0            depleted  

[80 rows x 9 columns]
```

### Text output (stream)

```text

### class_coverage_summary.tsv

```

### Text output (text/plain)

```text
                       regulator_class  minimum_evaluable_targets  \
0      Apoptosis / cell cycle / stress                         18   
1                 Chromatin regulation                         28   
2                          Hippo / YAP                          5   
3  Metabolism / trafficking / adhesion                         35   
4            Other reviewed regulators                         17   
5                      TGF/BMP/Activin                          7   
6       Transcription / RNA processing                         61   
7           Translation / proteostasis                         36   
8              VEGF / EC specification                          4   
9            WNT / mesoderm regulators                          5   

   maximum_evaluable_targets  n_day_program_combinations  
0                         19                          16  
1                         31                          16  
2                          6                          16  
3                         41                          16  
4                         18                          16  
5                          9                          16  
6                         68                          16  
7                         41                          16  
8                          4                          15  
9                          6                          16  
```

### Text output (stream)

```text

### coverage_and_filtering.tsv

```

### Text output (text/plain)

```text
                      analysis_population day target_name  n_analysis_cells  \
0   primary_curated_dominant_expected_day  D0       ACAA1             16350   
1   primary_curated_dominant_expected_day  D0      ACTL7B             16350   
2   primary_curated_dominant_expected_day  D0      ACVR1B             16350   
3   primary_curated_dominant_expected_day  D0      ACVR2A             16350   
4   primary_curated_dominant_expected_day  D0      ACVR2B             16350   
..                                    ...  ..         ...               ...   
75  primary_curated_dominant_expected_day  D0       ERCC2             16350   
76  primary_curated_dominant_expected_day  D0       ERO1A             16350   
77  primary_curated_dominant_expected_day  D0        ETS1             16350   
78  primary_curated_dominant_expected_day  D0        ETV2             16350   
79  primary_curated_dominant_expected_day  D0       EWSR1             16350   

    n_carrier_cells  passes_carrier_threshold  n_evaluable_program_states  
0               271                      True                           3  
1               252                      True                           3  
2               283                      True                           3  
3               309                      True                           3  
4               349                      True                           3  
..              ...                       ...                         ...  
75              260                      True                           3  
76              374                      True                           3  
77              389                      True                           3  
78              365                      True                           3  
79              291                      True                           3  

[80 rows x 7 columns]
```

### Text output (stream)

```text

### target_class_program_state_raw_distribution.tsv

```

### Text output (text/plain)

```text
                    regulator_class day assigned_state  \
0   Apoptosis / cell cycle / stress  D0            P22   
1   Apoptosis / cell cycle / stress  D0            P29   
2   Apoptosis / cell cycle / stress  D0             P7   
3   Apoptosis / cell cycle / stress  D1            P13   
4   Apoptosis / cell cycle / stress  D1            P28   
..                              ...  ..            ...   
75        Other reviewed regulators  D3             P2   
76        Other reviewed regulators  D3             P3   
77        Other reviewed regulators  D3            P33   
78        Other reviewed regulators  D3            P34   
79        Other reviewed regulators  D3             P8   

    n_unique_carrier_cells  n_targets_observed  n_state_cells  \
0                     1011                  19           3951   
1                      308                  19           1262   
2                     2946                  19          11137   
3                     1158                  19           4298   
4                     1104                  19           3615   
..                     ...                 ...            ...   
75                    1361                  18           8499   
76                    2058                  18          11993   
77                     786                  18           2688   
78                     416                  18           1186   
79                     797                  18           5094   

    raw_fraction_of_state_cells_with_class_target  
0                                        0.255885  
1                                        0.244057  
2                                        0.264524  
3                                        0.269428  
4                                        0.305394  
..                                            ...  
75                                       0.160136  
76                                       0.171600  
77                                       0.292411  
78                                       0.350759  
79                                       0.156459  

[80 rows x 7 columns]
```

### Text output (stream)

```text

### analysis_audit.json
{
  "n_all_cells": 121616,
  "n_primary_program_state_cells": 69453,
  "n_targets_with_assignments": 298,
  "n_target_program_state_tests": 4752,
  "n_classes_displayed": 10,
  "state_definition": "global dominant cNMF program; primary curated program list; expected-day match",
  "target_to_program_differential_effects_used": false
}

```

## Markdown cell 14

## 7. Display the publication panel and supplementary heatmap

## Code cell 15

```python
main_png = result_dir / "Fig2_perturbation_class_program_state_occupancy.png"
supp_png = result_dir / "Fig2_target_program_state_occupancy_heatmap_supplementary.png"

display(Image(filename=str(main_png)))
display(Image(filename=str(supp_png)))

```

### Text output (text/plain)

```text
<IPython.core.display.Image object>
```

### Text output (text/plain)

```text
<IPython.core.display.Image object>
```

## Markdown cell 16

## 8. Package the completed outputs

## Code cell 17

```python
zip_path = PACKAGE_DIR / "program_defined_occupancy_results.zip"
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for p in result_dir.rglob("*"):
        if p.is_file():
            zf.write(p, arcname=str(p.relative_to(result_dir.parent)))
    for p in export_dir.glob("*.tsv"):
        zf.write(p, arcname=str(p.relative_to(export_dir.parent)))
    for p in export_dir.glob("*.json"):
        zf.write(p, arcname=str(p.relative_to(export_dir.parent)))
print("Wrote:", zip_path)

```

### Text output (stream)

```text
Wrote: /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/figure2_program_defined_class_occupancy_20260728/program_defined_occupancy_results.zip

```

