# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/program_level_temporal_map_notebooks/03_program_high_cell_DE_optional_annotation.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/program_level_temporal_map_notebooks/03_program_high_cell_DE_optional_annotation.ipynb`
- SHA-256: `e8cab2b1d15b1a3ebd64f1e585ad83909a49dc21eecf8932d0b5ce7f94f795a0`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Markdown cell 0

# 03. Optional: program-high cell DE for biological annotation

This notebook is optional. It is meant to answer questions like Kyle's in a program-centered way: compare cells high for P17 vs cells high for P14, etc. This should be described as **program-defined cell annotation**, not as the primary perturbation-effect evidence.

It requires an AnnData object named `adata_int` already loaded in the notebook.

## Code cell 1

```python

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy import sparse

try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except Exception:
    HAS_ADJUST_TEXT = False

OUTDIR = "outputs_program_high_cell_DE_annotation_v4"
os.makedirs(OUTDIR, exist_ok=True)

DAY_COL = "day"
DAY_TO_USE = "D2"
PROGRAM_COLS = ["Program 14", "Program 15", "Program 17"]
FDR_CUTOFF = 0.05
LOGFC_CUTOFF = 0.25

```

## Markdown cell 2

## Check AnnData and required columns

## Code cell 3

```python

# The user should load adata_int before running this notebook.
print(adata_int)
print("obs columns containing Program:", [c for c in adata_int.obs.columns if c.startswith("Program ")][:20])
missing = [c for c in PROGRAM_COLS if c not in adata_int.obs.columns]
if missing:
    raise ValueError(f"Missing program usage columns in adata_int.obs: {missing}")

# Clean day column if present.
def clean_day(x):
    s = str(x).strip().lower().replace("day", "d").replace(" ", "")
    if s in ["0", "0.0", "d0"]: return "D0"
    if s in ["1", "1.0", "d1"]: return "D1"
    if s in ["2", "2.0", "d2"]: return "D2"
    if s in ["3", "3.0", "d3"]: return "D3"
    return str(x).strip().upper()

if DAY_COL in adata_int.obs.columns:
    adata_int.obs["day_clean"] = adata_int.obs[DAY_COL].map(clean_day)
    day_mask = adata_int.obs["day_clean"].values == DAY_TO_USE
else:
    print("No day column found. Using all cells.")
    day_mask = np.ones(adata_int.n_obs, dtype=bool)

adata_d2 = adata_int[day_mask].copy()
print(adata_d2)

```

## Markdown cell 4

## Define program-high cell groups

Cells are assigned to the program for which they are high and relatively specific. You can tune quantiles.

## Code cell 5

```python

QUANTILE = 0.80
MIN_MARGIN = 0.05

for col in PROGRAM_COLS:
    adata_d2.obs[col] = pd.to_numeric(adata_d2.obs[col], errors="coerce")

# thresholds per program
thresholds = {col: adata_d2.obs[col].quantile(QUANTILE) for col in PROGRAM_COLS}
thresholds

```

## Code cell 6

```python

# Create a simple program-high label.
prog_values = adata_d2.obs[PROGRAM_COLS].copy()
max_prog = prog_values.idxmax(axis=1)
max_val = prog_values.max(axis=1)
second_val = prog_values.apply(lambda x: np.sort(x.values)[-2], axis=1)

adata_d2.obs["program_high_group"] = "other"
for col in PROGRAM_COLS:
    mask = (max_prog == col) & (max_val >= thresholds[col]) & ((max_val - second_val) >= MIN_MARGIN)
    adata_d2.obs.loc[mask, "program_high_group"] = col.replace("Program ", "P") + "-high"

print(adata_d2.obs["program_high_group"].value_counts())

```

## Markdown cell 7

## Run DE between program-high groups

## Code cell 8

```python

def run_de_between_groups(adata, group_col, group_a, group_b, out_prefix, method="wilcoxon"):
    sub = adata[adata.obs[group_col].isin([group_a, group_b])].copy()
    sub.obs[group_col] = pd.Categorical(sub.obs[group_col].astype(str), categories=[group_a, group_b])
    print(out_prefix, sub.obs[group_col].value_counts())
    sc.tl.rank_genes_groups(sub, groupby=group_col, groups=[group_a], reference=group_b, method=method, use_raw=False, pts=True)
    de = sc.get.rank_genes_groups_df(sub, group=group_a).copy()
    if "gene_symbol" in sub.var.columns:
        # only works if names are var_names; otherwise ignore
        try:
            de["gene_symbol"] = sub.var.loc[de["names"], "gene_symbol"].values
        except Exception:
            de["gene_symbol"] = de["names"]
    else:
        de["gene_symbol"] = de["names"]
    de["comparison"] = f"{group_a}_vs_{group_b}"
    de.to_csv(os.path.join(OUTDIR, f"DE_{out_prefix}.csv"), index=False)
    return de

comparisons = [
    ("P17-high", "P14-high", "P17_high_vs_P14_high"),
    ("P17-high", "P15-high", "P17_high_vs_P15_high"),
    ("P14-high", "P15-high", "P14_high_vs_P15_high"),
]

de_results = {}
for a, b, pref in comparisons:
    if a in set(adata_d2.obs["program_high_group"]) and b in set(adata_d2.obs["program_high_group"]):
        de_results[pref] = run_de_between_groups(adata_d2, "program_high_group", a, b, pref)
    else:
        print("Skipping", pref, "because one group is absent.")

```

## Markdown cell 9

## Volcano plots with label repel

## Code cell 10

```python

def volcano(de, title, out_prefix, label_genes=None, top_n=12):
    df = de.copy().dropna(subset=["logfoldchanges", "pvals_adj"])
    df["neglog10_fdr"] = -np.log10(df["pvals_adj"].clip(lower=1e-300))
    df["sig"] = (df["pvals_adj"] < FDR_CUTOFF) & (df["logfoldchanges"].abs() > LOGFC_CUTOFF)
    fig, ax = plt.subplots(figsize=(5.4, 4.8), dpi=300)
    ax.scatter(df.loc[~df.sig,"logfoldchanges"], df.loc[~df.sig,"neglog10_fdr"], s=8, color="#BDBDBD", alpha=0.35, edgecolor='none')
    ax.scatter(df.loc[df.sig,"logfoldchanges"], df.loc[df.sig,"neglog10_fdr"], s=16, color="#D73027", alpha=0.75, edgecolor='none')
    ax.axvline(0, color="#888888", ls="--", lw=0.8)
    ax.axhline(-np.log10(FDR_CUTOFF), color="#AAAAAA", ls="--", lw=0.7)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("log fold-change")
    ax.set_ylabel("-log10(FDR)")
    # labels: top significant plus custom genes
    lab = df[df.sig].sort_values("pvals_adj").head(top_n).copy()
    if label_genes:
        lab = pd.concat([lab, df[df["gene_symbol"].astype(str).isin(label_genes)]], ignore_index=True).drop_duplicates("gene_symbol")
    texts = []
    for _, row in lab.iterrows():
        texts.append(ax.text(row["logfoldchanges"], row["neglog10_fdr"], str(row["gene_symbol"]), fontsize=8))
    if HAS_ADJUST_TEXT and texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', lw=0.4, color='#666666', alpha=0.7), time_lim=1.5)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, f"volcano_{out_prefix}.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUTDIR, f"volcano_{out_prefix}.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(OUTDIR, f"volcano_{out_prefix}.svg"), bbox_inches="tight")
    plt.show()

custom = ["KRT8", "KRT18", "ANKRD1", "HAND1", "TBX20", "CDH11", "TAL1", "FLI1", "KDR", "SOX7"]
for pref, de in de_results.items():
    volcano(de, pref.replace("_", " "), pref, label_genes=custom)

```

## Markdown cell 11

## Dot plot of top genes for program-defined groups

## Code cell 12

```python

# Choose top genes from each DE comparison to summarize.
top_genes = []
for pref, de in de_results.items():
    top_genes.extend(de.sort_values("pvals_adj").head(8)["names"].tolist())
top_genes = list(dict.fromkeys(top_genes))[:40]
print(top_genes)

if len(top_genes) > 0:
    sc.pl.dotplot(
        adata_d2,
        var_names=top_genes,
        groupby="program_high_group",
        categories_order=["P14-high", "P15-high", "P17-high", "other"],
        standard_scale="var",
        dendrogram=False,
        show=True,
    )

```

## Markdown cell 13

## Package outputs

## Code cell 14

```python

import zipfile
zip_path = OUTDIR + ".zip"
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
    for root, _, files in os.walk(OUTDIR):
        for fn in files:
            fp = os.path.join(root, fn)
            z.write(fp, arcname=os.path.relpath(fp, OUTDIR))
print("Packaged:", zip_path)

```

