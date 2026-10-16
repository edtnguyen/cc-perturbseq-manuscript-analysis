# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/remake_D2LM_CLR_heatmap_with_control_and_stats.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/remake_D2LM_CLR_heatmap_with_control_and_stats.ipynb`
- SHA-256: `cf2d754df609269887a6eb70438b19f8656d182cfa472dbdc04cc6049a9dc780`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# D2 LM TDI/IAG selected-program heatmap on the CLR scale

This notebook recomputes small-molecule program effects from the saved zero-based projected AnnData object.

Workflow:

1. apply the same cell QC filters used previously;
2. apply a log2 CLR transform to each cell's 50-program usage composition;
3. average CLR scores within each biological sample;
4. calculate matched replicate effects for D2 LM TDI and IAG relative to D2 LM control;
5. run paired tests across biological replicates;
6. BH-adjust across the 6 treatment × program tests;
7. draw a 3 × 3 square-cell heatmap with a white control/control reference row.

The control row is a constructed zero reference and is not assigned a p-value or significance star.

## Code cell 1

```python
from pathlib import Path
import re

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import ttest_1samp
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

try:
    display
except NameError:
    def display(x):
        print(x if not hasattr(x, "to_string") else x.to_string(index=False))

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "figure.dpi": 160,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",  # editable text in Illustrator
})

BASE = Path("/scratch/users/edtng/20260325_Hippo_timecourse")
H5AD_PATH = (
    BASE / "perturbseq_cNMF_projection" /
    "Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad"
)
OUT_DIR = BASE / "perturbseq_cNMF_projection" / "D2LM_CLR_selected_program_heatmap"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SCORE_KEY = "X_perturbseq_cNMF_usage"
SAMPLE_COL = "sample_name"

# Figure order requested for the D2 LM panel.
PROGRAMS = ["P15", "P14", "P17"]
PROGRAM_LABELS = {
    "P15": "P15\nendothelial / pre-EC",
    "P14": "P14\nlateral mesoderm",
    "P17": "P17\nepithelial-like",
}

CONTROL_CONDITION = "D2_LM"
TREATED_CONDITIONS = {
    "TDI": "D2_LM+TDI",
    "IAG": "D2_LM+IAG",
}
ROW_ORDER = ["Control", "TDI", "IAG"]

# Match the programDE CLR scale as closely as possible.
# Set this to the exact epsilon used in the programDE run if it differed.
CLR_EPS = 1e-6
CLR_LOG_BASE = 2

# Stars default to BH-adjusted q-values across the 8 treatment × program tests.
# Change to "pvalue" only if you explicitly want unadjusted paired-test stars.
STAR_SOURCE = "qvalue"

print("H5AD:", H5AD_PATH, H5AD_PATH.exists())
print("Output:", OUT_DIR)
```

### Text output (stream)

```text
H5AD: /scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_projection/Hippo_timecourse_combined.projected_perturbseq_cNMF_50_dt_0_2.zero_based.h5ad True
Output: /scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_projection/D2LM_CLR_selected_program_heatmap

```

## Code cell 2

```python
# Same condition-specific QC cutoffs used in the preceding time-course analyses.
cond_pct_ribo_cutoffs = {
    'D0_hPSC': 2.4,
    'D1_MPS': 2.3,
    'D1_MPS+TDI': 2.5,
    'D1_MPS+IAG': 2.2,
    'D1.33_LM': 1.9,
    'D1.33_LM+TDI': 2.1,
    'D1.33_LM+IAG': 2.2,
    'D1.67_LM': 2.0,
    'D1.67_LM+TDI': 2.2,
    'D1.67_LM+IAG': 1.9,
    'D2_LM': 2.0,
    'D2_LM+TDI': 2.4,
    'D2_LM+IAG': 1.9,
    'D3_AE': 1.0,
    'D3_AE+TDI': 2.0,
    'D3_AE+IAG': 1.0,
    'D2_LM+TDI->D3_AE+TDI': 2.5,
    'D2_LM+IAG->D3_AE+IAG': 1.0,
}

cond_pct_mt_cutoffs = {
    'D0_hPSC': 5.5,
    'D1_MPS': 3.0,
    'D1_MPS+TDI': 3.0,
    'D1_MPS+IAG': 4,
    'D1.33_LM': 3.5,
    'D1.33_LM+TDI': 3.75,
    'D1.33_LM+IAG': 5,
    'D1.67_LM': 4.0,
    'D1.67_LM+TDI': 5.5,
    'D1.67_LM+IAG': 3.5,
    'D2_LM': 6.5,
    'D2_LM+TDI': 8,
    'D2_LM+IAG': 3.5,
    'D3_AE': 3.25,
    'D3_AE+TDI': 7.0,
    'D3_AE+IAG': 3.75,
    'D2_LM+TDI->D3_AE+TDI': 7,
    'D2_LM+IAG->D3_AE+IAG': 3.75,
}


def parse_sample_name(sample_name):
    sample_name = str(sample_name)
    m = re.match(r"^(.*)_rep(\d+)$", sample_name)
    if m is None:
        raise ValueError(f"Could not parse condition/replicate from {sample_name!r}")
    return m.group(1), f"rep{m.group(2)}"


def clr_transform_rows(usage, eps=1e-6, log_base=2):
    """Floor, renormalize, and CLR-transform rows of a cells × programs matrix."""
    usage = np.asarray(usage, dtype=np.float64)
    if usage.ndim != 2:
        raise ValueError("usage must be 2D")
    if np.any(usage < 0):
        raise ValueError("usage contains negative values")

    row_sums = usage.sum(axis=1, keepdims=True)
    usage = np.divide(usage, row_sums, out=np.zeros_like(usage), where=row_sums > 0)
    floored = np.maximum(usage, eps)
    floored = floored / floored.sum(axis=1, keepdims=True)

    if log_base == 2:
        log_u = np.log2(floored)
    elif log_base == "natural":
        log_u = np.log(floored)
    else:
        raise ValueError("log_base must be 2 or 'natural'")

    clr = log_u - log_u.mean(axis=1, keepdims=True)
    return clr.astype(np.float32)


def bh_adjust(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    qvalues = np.full_like(pvalues, np.nan, dtype=float)
    mask = np.isfinite(pvalues)
    p = pvalues[mask]
    if p.size == 0:
        return qvalues
    order = np.argsort(p)
    ranked = p[order]
    n = ranked.size
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    restored = np.empty_like(adj)
    restored[order] = adj
    qvalues[mask] = restored
    return qvalues


def significance_stars(value):
    if not np.isfinite(value):
        return ""
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return ""
```

## Code cell 3

```python
# Load only obs + obsm from the large object.
adata = sc.read_h5ad(H5AD_PATH, backed="r")
obs = adata.obs.copy()
usage = np.asarray(adata.obsm[SCORE_KEY])

projection_meta = adata.uns.get("perturbseq_cNMF_projection", {})
program_names = [str(x) for x in projection_meta.get("program_names", [])]
if len(program_names) != usage.shape[1]:
    program_names = [f"P{i}" for i in range(usage.shape[1])]

program_to_col = {p: i for i, p in enumerate(program_names)}
missing_programs = [p for p in PROGRAMS if p not in program_to_col]
if missing_programs:
    raise ValueError(f"Missing selected programs: {missing_programs}")

parsed = obs[SAMPLE_COL].astype(str).map(parse_sample_name)
obs["condition_qc"] = [x[0] for x in parsed]
obs["replicate_qc"] = [x[1] for x in parsed]

for col in ["pct_counts_mt", "pct_counts_ribo"]:
    obs[col] = pd.to_numeric(obs[col], errors="coerce")

# Apply condition-specific QC filter.
keep = np.ones(obs.shape[0], dtype=bool)
filter_records = []
for cond in obs["condition_qc"].unique():
    cond = str(cond)
    if cond not in cond_pct_mt_cutoffs or cond not in cond_pct_ribo_cutoffs:
        raise ValueError(f"No QC cutoff specified for {cond!r}")
    cond_mask = obs["condition_qc"].to_numpy() == cond
    fail_mt = cond_mask & (obs["pct_counts_mt"].to_numpy() >= cond_pct_mt_cutoffs[cond])
    fail_ribo = cond_mask & (obs["pct_counts_ribo"].to_numpy() <= cond_pct_ribo_cutoffs[cond])
    fail = fail_mt | fail_ribo
    keep &= ~fail
    filter_records.append({
        "condition": cond,
        "n_before": int(cond_mask.sum()),
        "n_removed": int(fail.sum()),
        "n_after": int(cond_mask.sum() - fail.sum()),
    })

filter_summary = pd.DataFrame(filter_records)
filter_summary.to_csv(OUT_DIR / "QC_filter_summary.tsv", sep="\t", index=False)

obs_f = obs.loc[keep].copy()
usage_f = usage[keep, :]
clr_f = clr_transform_rows(usage_f, eps=CLR_EPS, log_base=CLR_LOG_BASE)

print("Cells before QC:", obs.shape[0])
print("Cells after QC:", obs_f.shape[0])
print("CLR row mean range:", float(clr_f.mean(axis=1).min()), float(clr_f.mean(axis=1).max()))

# Sample-level mean CLR values across all programs.
sample_records = []
for sample_name in obs_f[SAMPLE_COL].astype(str).drop_duplicates():
    mask = obs_f[SAMPLE_COL].astype(str).to_numpy() == sample_name
    condition, replicate = parse_sample_name(sample_name)
    means = clr_f[mask, :].mean(axis=0)
    rec = {
        "sample_name": sample_name,
        "condition": condition,
        "replicate": replicate,
        "n_cells": int(mask.sum()),
    }
    rec.update({p: float(means[j]) for j, p in enumerate(program_names)})
    sample_records.append(rec)

sample_clr = pd.DataFrame(sample_records)
sample_clr.to_csv(OUT_DIR / "sample_mean_cell_level_CLR_usage.tsv", sep="\t", index=False)
display(sample_clr[["sample_name", "condition", "replicate", "n_cells"] + PROGRAMS])

adata.file.close()
```

### Text output (stream)

```text
Cells before QC: 212891
Cells after QC: 99416
CLR row mean range: -3.623962356869015e-07 3.623962356869015e-07

```

### Text output (text/plain)

```text
                  sample_name             condition replicate  n_cells  \
0              D2_LM+IAG_rep2             D2_LM+IAG      rep2     2601   
1                  D2_LM_rep1                 D2_LM      rep1     2909   
2              D2_LM+IAG_rep1             D2_LM+IAG      rep1     2621   
3   D2_LM+TDI->D3_AE+TDI_rep1  D2_LM+TDI->D3_AE+TDI      rep1     3304   
4                D0_hPSC_rep1               D0_hPSC      rep1     2884   
5                  D2_LM_rep3                 D2_LM      rep3     3381   
6               D1.33_LM_rep2              D1.33_LM      rep2     1837   
7                D0_hPSC_rep2               D0_hPSC      rep2     3317   
8           D1.33_LM+TDI_rep2          D1.33_LM+TDI      rep2     1728   
9             D1_MPS+IAG_rep3            D1_MPS+IAG      rep3     2355   
10                D1_MPS_rep2                D1_MPS      rep2     1579   
11              D1.67_LM_rep3              D1.67_LM      rep3     1298   
12             D3_AE+IAG_rep3             D3_AE+IAG      rep3     1706   
13  D2_LM+TDI->D3_AE+TDI_rep2  D2_LM+TDI->D3_AE+TDI      rep2     3152   
14  D2_LM+TDI->D3_AE+TDI_rep3  D2_LM+TDI->D3_AE+TDI      rep3     3097   
15             D3_AE+TDI_rep3             D3_AE+TDI      rep3     2337   
16  D2_LM+IAG->D3_AE+IAG_rep2  D2_LM+IAG->D3_AE+IAG      rep2     1021   
17             D2_LM+TDI_rep3             D2_LM+TDI      rep3     1775   
18             D2_LM+TDI_rep1             D2_LM+TDI      rep1     1745   
19            D1_MPS+TDI_rep3            D1_MPS+TDI      rep3     1803   
20          D1.67_LM+IAG_rep2          D1.67_LM+IAG      rep2     3042   
21              D1.33_LM_rep1              D1.33_LM      rep1     1631   
22                 D3_AE_rep3                 D3_AE      rep3     1317   
23            D1_MPS+IAG_rep1            D1_MPS+IAG      rep1     1785   
24                 D3_AE_rep2                 D3_AE      rep2      999   
25             D2_LM+TDI_rep2             D2_LM+TDI      rep2     2347   
26               D0_hPSC_rep3               D0_hPSC      rep3     2691   
27            D1_MPS+IAG_rep2            D1_MPS+IAG      rep2     2170   
28          D1.33_LM+TDI_rep1          D1.33_LM+TDI      rep1     1364   
29                 D2_LM_rep2                 D2_LM      rep2     3325   
30             D3_AE+IAG_rep1             D3_AE+IAG      rep1      857   
31             D3_AE+IAG_rep2             D3_AE+IAG      rep2     1082   
32          D1.67_LM+TDI_rep2          D1.67_LM+TDI      rep2     1244   
33          D1.33_LM+IAG_rep1          D1.33_LM+IAG      rep1     1725   
34            D1_MPS+TDI_rep2            D1_MPS+TDI      rep2       93   
35          D1.67_LM+IAG_rep3          D1.67_LM+IAG      rep3     3721   
36          D1.67_LM+IAG_rep1          D1.67_LM+IAG      rep1     3356   
37             D3_AE+TDI_rep2             D3_AE+TDI      rep2     1369   
38  D2_LM+IAG->D3_AE+IAG_rep1  D2_LM+IAG->D3_AE+IAG      rep1      586   
39             D2_LM+IAG_rep3             D2_LM+IAG      rep3     1854   
40          D1.33_LM+IAG_rep2          D1.33_LM+IAG      rep2     1507   
41  D2_LM+IAG->D3_AE+IAG_rep3  D2_LM+IAG->D3_AE+IAG      rep3      880   
42             D3_AE+TDI_rep1             D3_AE+TDI      rep1     1715   
43          D1.33_LM+TDI_rep3          D1.33_LM+TDI      rep3     1707   
44            D1_MPS+TDI_rep1            D1_MPS+TDI      rep1     1185   
45              D1.33_LM_rep3              D1.33_LM      rep3     1565   
46                D1_MPS_rep1                D1_MPS      rep1     1743   
47          D1.67_LM+TDI_rep1          D1.67_LM+TDI      rep1      816   
48              D1.67_LM_rep1              D1.67_LM      rep1      334   
49                D1_MPS_rep3                D1_MPS      rep3     1248   
50          D1.33_LM+IAG_rep3          D1.33_LM+IAG      rep3     1418   
51          D1.67_LM+TDI_rep3          D1.67_LM+TDI      rep3     1363   
52              D1.67_LM_rep2              D1.67_LM      rep2      556   
53                 D3_AE_rep1                 D3_AE      rep1      371   

         P15        P14       P17  
0   5.878726  -0.703459 -1.575960  
1  -0.315764  -1.521883 -1.553798  
2   3.593900  -1.556272 -1.676567  
3  -1.507164  -0.400526 -0.627946  
4  -2.438244  -2.438244 -2.438244  
5   7.860542   2.622862 -0.935695  
6  -1.898087  -0.360845 -2.012710  
7  -2.784018  -2.085970 -2.802536  
8  -2.358032  -2.048707 -2.042578  
9  -2.737157  -2.742089 -2.750486  
10 -2.352300  -2.352300 -2.352300  
11 -0.035185   0.994644 -1.431794  
12  5.326709  -2.113669 -2.280143  
13 -1.869558   5.742752  3.450738  
14 -1.591678   6.663464  0.843259  
15  5.321491   4.166621 -1.284109  
16  2.289388   1.449680 -1.314581  
17 -0.538188  13.309165  8.999609  
18 -1.785930  -1.624740 -0.647238  
19 -2.561019  -2.546743 -2.553547  
20  5.813777  -0.110403 -2.084246  
21 -2.071413  -0.174780 -2.197376  
22  3.846312  -2.174086 -2.209179  
23 -3.088985  -3.096180 -3.096180  
24  0.269533  -0.535804 -1.330704  
25 -1.613415   1.210178  5.226427  
26 -2.049848  -1.497250 -2.050397  
27 -3.240085  -3.235714 -3.240085  
28 -1.821136  -0.196045 -1.272537  
29  2.708675  -0.201559 -1.543957  
30  1.324408  -2.096071 -2.096071  
31  5.305177   0.462613 -1.661790  
32 -1.679862  -1.339297 -1.021578  
33 -1.345744   9.478078 -2.489689  
34 -2.513944  -2.513944 -2.513944  
35 -1.030043  -1.607165 -1.721891  
36  0.734543   0.123693 -1.580789  
37  5.019238   1.980059 -1.588721  
38  0.831662   0.293874 -1.338107  
39  7.183838  -1.685515 -2.007866  
40 -1.648330  -0.377174 -1.795291  
41  1.574905  -0.924300 -1.690490  
42 -0.541320  -1.207495 -1.499839  
43 -1.783705  -1.178236 -1.619311  
44 -2.005963  -1.970431 -1.982672  
45 -1.889906   1.448594 -1.844874  
46 -2.224225  -2.224225 -2.224225  
47 -1.539125  -1.136023 -0.746497  
48  4.962663   4.668322 -1.323697  
49 -2.628102  -2.628102 -2.628102  
50 -2.009251  -1.886855 -2.067042  
51 -1.522835  -0.917177  1.776510  
52  1.109866  -0.812171 -1.578023  
53  0.850661   0.580479 -1.358960  
```

## Code cell 4

```python
# Matched replicate CLR effects and paired tests.
sample_index = sample_clr.set_index(["condition", "replicate"])
records = []
paired_records = []

for treatment, treated_condition in TREATED_CONDITIONS.items():
    control_reps = set(sample_clr.loc[sample_clr["condition"] == CONTROL_CONDITION, "replicate"])
    treated_reps = set(sample_clr.loc[sample_clr["condition"] == treated_condition, "replicate"])
    matched_reps = sorted(control_reps & treated_reps)
    if len(matched_reps) < 2:
        raise ValueError(f"Only {len(matched_reps)} matched replicates for {treatment}")

    for program in PROGRAMS:
        control_values = np.array([
            sample_index.loc[(CONTROL_CONDITION, rep), program] for rep in matched_reps
        ], dtype=float)
        treated_values = np.array([
            sample_index.loc[(treated_condition, rep), program] for rep in matched_reps
        ], dtype=float)
        differences = treated_values - control_values

        test = ttest_1samp(differences, popmean=0.0, nan_policy="omit")
        records.append({
            "row": treatment,
            "treated_condition": treated_condition,
            "control_condition": CONTROL_CONDITION,
            "program": program,
            "n_matched_pairs": len(matched_reps),
            "mean_control_clr": float(control_values.mean()),
            "mean_treated_clr": float(treated_values.mean()),
            "mean_diff_clr": float(differences.mean()),
            "sd_diff_clr": float(differences.std(ddof=1)) if len(differences) > 1 else np.nan,
            "paired_t_statistic": float(test.statistic),
            "pvalue": float(test.pvalue),
            "n_positive_diffs": int((differences > 0).sum()),
            "n_negative_diffs": int((differences < 0).sum()),
        })

        for rep, c, t, d in zip(matched_reps, control_values, treated_values, differences):
            paired_records.append({
                "row": treatment,
                "program": program,
                "replicate": rep,
                "control_clr": c,
                "treated_clr": t,
                "paired_diff_clr": d,
            })

effects = pd.DataFrame(records)
effects["qvalue_BH_8tests"] = bh_adjust(effects["pvalue"].to_numpy())
paired_effects = pd.DataFrame(paired_records)

effects.to_csv(OUT_DIR / "D2LM_CLR_selected_program_effects_source.tsv", sep="\t", index=False)
paired_effects.to_csv(OUT_DIR / "D2LM_CLR_selected_program_paired_replicate_values.tsv", sep="\t", index=False)

display(effects)
display(paired_effects)
```

### Text output (text/plain)

```text
   row treated_condition control_condition program  n_matched_pairs  \
0  TDI         D2_LM+TDI             D2_LM     P15                3   
1  TDI         D2_LM+TDI             D2_LM     P14                3   
2  TDI         D2_LM+TDI             D2_LM     P17                3   
3  IAG         D2_LM+IAG             D2_LM     P15                3   
4  IAG         D2_LM+IAG             D2_LM     P14                3   
5  IAG         D2_LM+IAG             D2_LM     P17                3   

   mean_control_clr  mean_treated_clr  mean_diff_clr  sd_diff_clr  \
0          3.417818         -1.312511      -4.730328     3.482276   
1          0.299807          4.298201       3.998394     5.841197   
2         -1.344483          4.526266       5.870749     4.581110   
3          3.417818          5.552155       2.134337     2.462361   
4          0.299807         -1.315082      -1.614889     2.344313   
5         -1.344483         -1.753464      -0.408981     0.576129   

   paired_t_statistic    pvalue  n_positive_diffs  n_negative_diffs  \
0           -2.352820  0.142913                 0                 3   
1            1.185617  0.357546                 2                 1   
2            2.219645  0.156633                 3                 0   
3            1.501315  0.272093                 2                 1   
4           -1.193130  0.355165                 0                 3   
5           -1.229543  0.343884                 0                 3   

   qvalue_BH_8tests  
0          0.357546  
1          0.357546  
2          0.357546  
3          0.357546  
4          0.357546  
5          0.357546  
```

### Text output (text/plain)

```text
    row program replicate  control_clr  treated_clr  paired_diff_clr
0   TDI     P15      rep1    -0.315764    -1.785930        -1.470166
1   TDI     P15      rep2     2.708675    -1.613415        -4.322090
2   TDI     P15      rep3     7.860542    -0.538188        -8.398730
3   TDI     P14      rep1    -1.521883    -1.624740        -0.102856
4   TDI     P14      rep2    -0.201559     1.210178         1.411736
5   TDI     P14      rep3     2.622862    13.309165        10.686303
6   TDI     P17      rep1    -1.553798    -0.647238         0.906559
7   TDI     P17      rep2    -1.543957     5.226427         6.770384
8   TDI     P17      rep3    -0.935695     8.999609         9.935304
9   IAG     P15      rep1    -0.315764     3.593900         3.909665
10  IAG     P15      rep2     2.708675     5.878726         3.170051
11  IAG     P15      rep3     7.860542     7.183838        -0.676704
12  IAG     P14      rep1    -1.521883    -1.556272        -0.034389
13  IAG     P14      rep2    -0.201559    -0.703459        -0.501900
14  IAG     P14      rep3     2.622862    -1.685515        -4.308377
15  IAG     P17      rep1    -1.553798    -1.676567        -0.122769
16  IAG     P17      rep2    -1.543957    -1.575960        -0.032003
17  IAG     P17      rep3    -0.935695    -2.007866        -1.072171
```

## Code cell 5

```python
# Build and plot Control/TDI/IAG × P15/P14/P17 heatmap.
value_matrix = pd.DataFrame(0.0, index=ROW_ORDER, columns=PROGRAMS)
p_matrix = pd.DataFrame(np.nan, index=ROW_ORDER, columns=PROGRAMS)
q_matrix = pd.DataFrame(np.nan, index=ROW_ORDER, columns=PROGRAMS)

for _, row in effects.iterrows():
    value_matrix.loc[row["row"], row["program"]] = row["mean_diff_clr"]
    p_matrix.loc[row["row"], row["program"]] = row["pvalue"]
    q_matrix.loc[row["row"], row["program"]] = row["qvalue_BH_8tests"]

value_matrix.to_csv(OUT_DIR / "D2LM_CLR_selected_program_heatmap_matrix.tsv", sep="\t")
p_matrix.to_csv(OUT_DIR / "D2LM_CLR_selected_program_heatmap_pvalues.tsv", sep="\t")
q_matrix.to_csv(OUT_DIR / "D2LM_CLR_selected_program_heatmap_qvalues.tsv", sep="\t")

values = value_matrix.to_numpy(dtype=float)
vmax = np.nanmax(np.abs(values[1:, :]))
if not np.isfinite(vmax) or vmax == 0:
    vmax = 1.0
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

# 3 columns x 3 rows; equal aspect makes every grid cell square.
fig, ax = plt.subplots(figsize=(6.2, 5.6))
mesh = ax.pcolormesh(
    np.arange(len(PROGRAMS) + 1),
    np.arange(len(ROW_ORDER) + 1),
    values,
    cmap="RdBu_r",
    norm=norm,
    shading="flat",
    linewidth=1.0,
    edgecolors="white",
)
ax.set_aspect("equal", adjustable="box")
ax.set_xlim(0, len(PROGRAMS))
ax.set_ylim(len(ROW_ORDER), 0)
ax.set_xticks(np.arange(len(PROGRAMS)) + 0.5)
ax.set_yticks(np.arange(len(ROW_ORDER)) + 0.5)
ax.set_xticklabels([PROGRAM_LABELS[p] for p in PROGRAMS], rotation=45, ha="right")
ax.set_yticklabels(ROW_ORDER)
ax.set_xlabel("Projected perturb-seq cNMF program")
ax.set_ylabel("D2 LM condition")
ax.set_title("D2 LM treatment effect on selected D2 programs\ncell-level log2 CLR, matched biological replicates")
ax.tick_params(length=0)

star_matrix = q_matrix if STAR_SOURCE == "qvalue" else p_matrix
for i, row_name in enumerate(ROW_ORDER):
    for j, program in enumerate(PROGRAMS):
        val = value_matrix.loc[row_name, program]
        if row_name == "Control":
            annotation = "0.00"
        else:
            sig_value = star_matrix.loc[row_name, program]
            annotation = f"{val:.2f}{significance_stars(sig_value)}"
        # Use white text on darkest cells.
        rgba = mpl.colormaps["RdBu_r"](norm(val))
        luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
        text_color = "white" if luminance < 0.42 else "black"
        ax.text(j + 0.5, i + 0.5, annotation, ha="center", va="center", color=text_color, fontsize=10)

cbar = fig.colorbar(mesh, ax=ax, fraction=0.045, pad=0.04)
cbar.set_label("Δ mean log2 CLR usage\ncondition − matched D2 LM control")

star_label = "BH q-value across 6 tests" if STAR_SOURCE == "qvalue" else "unadjusted paired t-test p-value"
fig.text(0.5, 0.01, f"Stars: * <0.05, ** <0.01, *** <0.001 ({star_label}); Control row is a zero reference.", ha="center", fontsize=7)
fig.subplots_adjust(left=0.18, right=0.88, bottom=0.28, top=0.84)

stem = OUT_DIR / "D2LM_TDI_IAG_selected_program_CLR_heatmap_control_reference_square_cells"
for ext in ["svg", "pdf", "png"]:
    fig.savefig(stem.with_suffix(f".{ext}"), bbox_inches="tight")
plt.show()
plt.close(fig)

print("Saved:", stem.with_suffix(".svg"))
```

### Text output (stream)

```text
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.
findfont: Font family 'Arial' not found.

```

### Text output (text/plain)

```text
<Figure size 992x896 with 2 Axes>
```

### Text output (stream)

```text
Saved: /scratch/users/edtng/20260325_Hippo_timecourse/perturbseq_cNMF_projection/D2LM_CLR_selected_program_heatmap/D2LM_TDI_IAG_selected_program_CLR_heatmap_control_reference_square_cells.svg

```

## Markdown cell 6

## Statistical interpretation

For every treatment × program cell, the test is a one-sample t-test of the three matched replicate CLR differences against zero. This is mathematically equivalent to a paired t-test of treated versus control sample-level CLR values. The displayed default stars use Benjamini–Hochberg adjustment across the six non-control heatmap cells.

With only three biological replicate pairs, the test has low power. A blank star field does not imply absence of a consistent direction; inspect the saved paired-replicate table as well.

