# Notebook mirror: `oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/tony_results/splitseq_analysis.ipynb`

- Source archive path: `/oak/stanford/groups/engreitz/Users/tri/CodeBase/cc-perturbseq-figure-reproduction/recovery/sherlock_notebooks/oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE/notebooks/tony_results/splitseq_analysis.ipynb`
- SHA-256: `bb601c110da8b798ec25aa9a0a85aa97715eac5db1e6b84b2781450b3754dc8c`
- Kernel language: `python`
- Full-fidelity notebook storage: Git LFS

## Code cell 0

```python
import scanpy as sc
import numpy as npp
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats as stats
import gffutils
```

## Markdown cell 1

### Process gene expression data (GUIDE CAPTURE)

## Code cell 2

```python
gex = sc.read_h5ad("GEX2.h5ad")

# set the gene expression count matrix to be the sum of all three types of reads
gex.X = gex.layers["ambiguous"]+gex.layers["mature"]+gex.layers["nascent"]
```

## Code cell 3

```python
# count the number of gene expression UMIs per cell
counts = gex.X.sum(axis=1)
counts = np.asarray(counts).flatten()
gex.obs["GEX_UMIS"] = counts
```

## Code cell 4

```python
def plot_knee(umi_counts, ncells):
    # sort UMI counts in descending order
    sorted_counts = np.sort(umi_counts)[::-1]

    plt.figure(figsize=(5,5))
    plt.plot(list(range(len(sorted_counts))), sorted_counts, linewidth=4, color="black")

    plt.xlabel("Barcode rank")
    plt.ylabel("UMI counts")

    plt.yscale("log")
    plt.xscale("log")

    plt.axvline(ncells, linewidth=2, linestyle='--', color="red")
    
    print(f"UMIs at {ncells}th cell: {sorted_counts[ncells]}")
```

## Code cell 5

```python
plt.rcParams["font.size"] = 24
plt.rcParams["font.family"] = "Arial"
plt.rcParams["figure.figsize"] = (5,5)

plot_knee(gex.obs["GEX_UMIS"], ncells=6000)
plt.savefig("DC_GEX_knee.png", dpi=300, bbox_inches="tight")
```

### Text output (stream)

```text
UMIs at 6000th cell: 1151.0000396

```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Code cell 6

```python
# only keep cells with >threshold UMIs
gex = gex[gex.obs["GEX_UMIS"]>1150]
```

## Code cell 7

```python
gex.obs["GEX_UMIS"].mean(), gex.obs["GEX_UMIS"].median()
```

### Text output (text/plain)

```text
(2778.3694621729237, 2474.9992466999997)
```

## Markdown cell 8

### Process guide data (GUIDE CAPTURE)

## Code cell 9

```python
guide = sc.read_h5ad("GUIDE1.h5ad")
guide = guide.to_df()
```

## Code cell 10

```python
# merge the gene expression and guide objects
gex.obs = gex.obs.merge(guide, left_index=True, right_index=True, how="left")
# if a guide is not detected in a cell, set its count to 0
gex.obs = gex.obs.fillna(0)
# count the number of guide UMIs per cell
gex.obs["GUIDE_UMIS"] = gex.obs[guide.columns].sum(axis=1)
```

## Code cell 11

```python
# mixing plot - want points NOT to lie along the x or y axes
#sns.scatterplot(x=gex.obs["GEX_UMIS"], y=gex.obs["GUIDE_UMIS"], alpha=0.1)
```

## Markdown cell 12

### Process gene expression data (CROP)

## Code cell 13

```python
gex_crop = sc.read_h5ad("CROP_GEX1.h5ad")
gex_crop.X = gex_crop.layers["ambiguous"]+gex_crop.layers["mature"]+gex_crop.layers["nascent"]
```

## Code cell 14

```python
counts = gex_crop.X.sum(axis=1)
counts = np.asarray(counts).flatten()
gex_crop.obs["GEX_UMIS"] = counts
```

## Code cell 15

```python
plot_knee(gex_crop.obs["GEX_UMIS"], ncells=6200)
plt.savefig("CROP_GEX_knee.png", dpi=300, bbox_inches="tight")
```

### Text output (stream)

```text
UMIs at 6200th cell: 1179.000992

```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Code cell 16

```python
gex_crop = gex_crop[gex_crop.obs["GEX_UMIS"]>1180]
```

## Code cell 17

```python
gex_crop.obs["GEX_UMIS"].mean(), gex_crop.obs["GEX_UMIS"].median()
```

### Text output (text/plain)

```text
(2844.388430896565, 2521.00037725)
```

## Markdown cell 18

### Process guide data (CROP)

## Code cell 19

```python
guide_crop = sc.read_h5ad("CROP_GUIDE1.h5ad")
guide_crop = guide_crop.to_df()
```

## Code cell 20

```python
gex_crop.obs = gex_crop.obs.merge(guide_crop, left_index=True, right_index=True, how="left")
gex_crop.obs = gex_crop.obs.fillna(0)
gex_crop.obs["GUIDE_UMIS"] = gex_crop.obs[guide_crop.columns].sum(axis=1)
```

## Code cell 21

```python
#sns.scatterplot(x=gex_crop.obs["GEX_UMIS"], y=gex_crop.obs["GUIDE_UMIS"])
```

## Markdown cell 22

### Histogram of number of gene expression UMIs for both libraries

## Code cell 23

```python
counts = pd.DataFrame({"Direct capture": gex.obs["GEX_UMIS"],
                       "CROP-seq": gex_crop.obs["GEX_UMIS"]})
# convert the table from a wide to a long table
counts = counts.melt().dropna()
counts.columns = ["Method", "# GEX UMIs"]
```

## Code cell 24

```python
sns.violinplot(x="# GEX UMIs", y="Method", 
               data=counts, log_scale=True,
               palette=["#ffd479","#76d6ff"])
plt.xlabel("# of gene expression UMIs")
plt.ylabel("")
plt.savefig("GEX_UMIS.png", dpi=300, bbox_inches="tight")
```

### Text output (stream)

```text
/var/folders/hz/kz7qqv9d2_d4sjp7b0zyz0y00000gq/T/ipykernel_66347/1119034430.py:1: FutureWarning: 

Passing `palette` without assigning `hue` is deprecated and will be removed in v0.14.0. Assign the `y` variable to `hue` and set `legend=False` for the same effect.

  sns.violinplot(x="# GEX UMIs", y="Method",

```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Markdown cell 25

### Histogram of number of guide UMIs per guide for both libraries

## Code cell 26

```python
guide_counts = pd.DataFrame({"Direct capture": gex.obs["GUIDE_UMIS"],
                             "CROP-seq": gex_crop.obs["GUIDE_UMIS"]})
guide_counts = guide_counts.melt().dropna()
guide_counts.columns = ["Method", "# GUIDE UMIs"]
```

## Code cell 27

```python
sns.violinplot(x="# GUIDE UMIs", y="Method", 
               data=guide_counts[guide_counts["# GUIDE UMIs"]>0], 
               log_scale=True, palette=["#ffd479","#76d6ff"])
plt.xlabel("# of guide UMIs")
plt.ylabel("")
plt.savefig("GUIDE_UMIS.png", dpi=300, bbox_inches="tight")
```

### Text output (stream)

```text
/var/folders/hz/kz7qqv9d2_d4sjp7b0zyz0y00000gq/T/ipykernel_66347/1850530268.py:1: FutureWarning: 

Passing `palette` without assigning `hue` is deprecated and will be removed in v0.14.0. Assign the `y` variable to `hue` and set `legend=False` for the same effect.

  sns.violinplot(x="# GUIDE UMIs", y="Method",

```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Code cell 28

```python
# calculate fraction of cells with nonzero 
guide_frac_nonzero = np.mean(gex.obs["GUIDE_UMIS"]>0)
crop_guide_frac_nonzero = np.mean(gex_crop.obs["GUIDE_UMIS"]>0)
```

## Code cell 29

```python
sns.barplot(y=["direct", "crop"],
            x=[guide_frac_nonzero, crop_guide_frac_nonzero])
# TODO: place a number on the barplot
```

### Text output (text/plain)

```text
<Axes: >
```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Code cell 30

```python
((gex.obs.iloc[:,1:]>=5).sum(axis=1)>0).mean()
```

### Text output (text/plain)

```text
0.9308448591901349
```

## Code cell 31

```python
((gex_crop.obs.iloc[:,1:]>=5).sum(axis=1)>0).mean()
```

### Text output (text/plain)

```text
0.897258064516129
```

## Markdown cell 32

### Distributions of counts for individual guides (GUIDE)

## Code cell 33

```python
guides = guide.columns.tolist()
```

## Code cell 34

```python
cutoff = 5

for g in ["TP53_+_7590780.23-P1P2"]:
    gcounts = gex.obs[g]
    gcounts = gcounts[gcounts>0]

    gcounts = pd.DataFrame(gcounts)
    gcounts["real"] = False
    gcounts.loc[gcounts[g]>cutoff,"real"] = True

    if gcounts["real"].sum()<100:
        continue
    
    sns.histplot(x=g, hue="real", data=gcounts, log_scale=True, 
                common_norm=False, stat="percent", alpha=0.5,
                palette=["#a9a9a9","#ffd479"])
    
    plt.axvline(cutoff, linestyle='--', color="red")
    plt.legend("", frameon=False)
    plt.savefig("DC_example_guide.png", dpi=300, bbox_inches="tight")
    plt.show()
```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Markdown cell 35

### Distributions of counts for individual guides (CROP)

## Code cell 36

```python
cutoff = 3

for g in ["TP53_+_7590780.23-P1P2"]:
    gcounts = gex_crop.obs[g]
    gcounts = gcounts[gcounts>0]

    gcounts = pd.DataFrame(gcounts)
    gcounts["real"] = False
    gcounts.loc[gcounts[g]>cutoff,"real"] = True

    if gcounts["real"].sum()<100:
        continue
    
    sns.histplot(x=g, hue="real", data=gcounts, log_scale=True, 
                 common_norm=False, stat="percent", alpha=0.5,
                 palette=["#a9a9a9","#76d6ff"])
    
    plt.axvline(cutoff, linestyle='--', color="red")
    plt.legend("", frameon=False)
    plt.savefig("CROP_example_guide.png",  dpi=300, bbox_inches="tight")
    plt.show()
```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Markdown cell 37

### Make saturation curves for guide libraries

## Code cell 38

```python
results = []

for counts in [1000000,3000000,10000000,30000000]:
    adata = sc.read_h5ad(f"GUIDE_S1.{counts}_kb_kite/counts_unfiltered_modified/adata.h5ad")
    cells = list(set(gex.obs.index) & set(adata.obs.index))
    adata = adata[cells]
    med_umis = np.median(np.asarray(adata.X.sum(axis=1)).flatten())
    nreads = counts/len(cells)
    results.append([nreads,med_umis,"Direct capture"])

for counts in [1000000,3000000,10000000,30000000]:
    adata = sc.read_h5ad(f"CROP_GUIDE1_S4.{counts}_kb_kite/counts_unfiltered_modified/adata.h5ad")
    cells = list(set(gex_crop.obs.index) & set(adata.obs.index))
    adata = adata[cells]
    med_umis = np.median(np.asarray(adata.X.sum(axis=1)).flatten())
    nreads = counts/len(cells)
    results.append([nreads,med_umis,"CROP-seq"])
```

## Code cell 39

```python
results = pd.DataFrame(results)
results.columns = ["Reads per cell","UMIs per cell","Method"]
```

## Code cell 40

```python
sns.lineplot(x="Reads per cell", y="UMIs per cell", hue="Method", data=results,
             palette=["#ffd479","#76d6ff"], linewidth=4, marker='o', markersize=15)
plt.xlabel("Number of reads per cell")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.savefig("guide_saturation_curve.png",  dpi=300, bbox_inches="tight")
```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Markdown cell 41

### Make UMAP (guide)

## Code cell 42

```python
db = gffutils.FeatureDB("/Users/tkzeng/Dropbox/engreitz/statin_pilot/gencode.v46.annotation.db")

id_to_name = {}

for gene in db.features_of_type("gene"):
    gene_id = gene["gene_id"][0]
    gene_name = gene["gene_name"][0]
    id_to_name[gene_id] = gene_name
```

## Code cell 43

```python
gex.var.index = gex.var.index.map(id_to_name)
gex.var_names_make_unique()

gex.var["mt"] = gex.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    gex, qc_vars=["mt"], inplace=True, log1p=True
)

# sc.pl.violin(
#     gex,
#     ["GEX_UMIS", "pct_counts_mt"],
#     jitter=0.4,
#     multi_panel=True,
# )

sc.pp.filter_cells(gex, min_genes=100)
sc.pp.filter_genes(gex, min_cells=3)

# Saving count data
gex.layers["counts"] = gex.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(gex)
# Logarithmize the data
sc.pp.log1p(gex)

sc.pp.highly_variable_genes(gex, n_top_genes=2000)

sc.tl.pca(gex)
sc.pp.neighbors(gex)
sc.tl.umap(gex)
```

### Text output (stream)

```text
/Users/tkzeng/miniconda3/envs/scvi-dev/lib/python3.11/site-packages/anndata/utils.py:260: UserWarning: Suffix used (-[0-9]+) to deduplicate index values may make index values difficult to interpret. There values with a similar suffixes in the index. Consider using a different delimiter by passing `join={delimiter}`Example key collisions generated by the make_index_unique algorithm: ['SNORD116-1', 'SNORD116-2', 'SNORD116-3', 'SNORD116-4', 'SNORD116-5']
  warnings.warn(

```

## Code cell 44

```python
# sc.pl.umap(gex)
```

## Markdown cell 45

### Calculate knockdown (guide)

## Code cell 46

```python
gex_copy = gex.copy()
    
guide_counts = gex_copy.obs[guides]
guide_assign = guide_counts[guide_counts>=5]
genes = [x.split('_')[0] for x in guides]
guide_assign.columns = genes

gene_assign = {}
for gene in set(genes):
    gene_assign[gene] = guide_assign[gene].sum(axis=1)>0
gene_assign = pd.DataFrame(gene_assign)
gex_copy.obs = gene_assign
```

## Code cell 47

```python
all_de = []
for gene in set(genes):
    gex_copy.obs["targeted"] = gex_copy.obs[gene].astype(str)
    sc.tl.rank_genes_groups(gex_copy, groupby="targeted", method="t-test_overestim_var")
    de = sc.get.rank_genes_groups_df(gex_copy, group="True")
    de = de[de["names"]==gene]
    de["ncells"] = np.sum(gex_copy.obs["targeted"]=="True")
    all_de.append(de)
all_de = pd.concat(all_de)
```

## Code cell 48

```python
#sns.histplot(all_de["logfoldchanges"][all_de["logfoldchanges"]>-15])
```

## Markdown cell 49

### Make UMAP (crop)

## Code cell 50

```python
gex_crop.var.index = gex_crop.var.index.map(id_to_name)
gex_crop.var_names_make_unique()

gex_crop.var["mt"] = gex_crop.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    gex_crop, qc_vars=["mt"], inplace=True, log1p=True
)

# sc.pl.violin(
#     gex_crop,
#     ["GEX_UMIS", "pct_counts_mt"],
#     jitter=0.4,
#     multi_panel=True,
# )

sc.pp.filter_cells(gex_crop, min_genes=100)
sc.pp.filter_genes(gex_crop, min_cells=3)

# Saving count data
gex_crop.layers["counts"] = gex_crop.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(gex_crop)
# Logarithmize the data
sc.pp.log1p(gex_crop)

sc.pp.highly_variable_genes(gex_crop, n_top_genes=2000)

sc.tl.pca(gex_crop)
sc.pp.neighbors(gex_crop)
sc.tl.umap(gex_crop)
```

### Text output (stream)

```text
/Users/tkzeng/miniconda3/envs/scvi-dev/lib/python3.11/site-packages/anndata/utils.py:260: UserWarning: Suffix used (-[0-9]+) to deduplicate index values may make index values difficult to interpret. There values with a similar suffixes in the index. Consider using a different delimiter by passing `join={delimiter}`Example key collisions generated by the make_index_unique algorithm: ['SNORD116-1', 'SNORD116-2', 'SNORD116-3', 'SNORD116-4', 'SNORD116-5']
  warnings.warn(

```

## Markdown cell 51

### Calculate knockdown (crop)

## Code cell 52

```python
gex_crop_copy = gex_crop.copy()

guide_counts = gex_crop_copy.obs[guides]
guide_assign = guide_counts[guide_counts>=3]
genes = [x.split('_')[0] for x in guides]
guide_assign.columns = genes

gene_assign = {}
for gene in set(genes):
    gene_assign[gene] = guide_assign[gene].sum(axis=1)>0
gene_assign = pd.DataFrame(gene_assign)
gex_crop_copy.obs = gene_assign
```

## Code cell 53

```python
all_de_crop = []
for gene in set(genes):
    gex_crop_copy.obs["targeted"] = gex_crop_copy.obs[gene].astype(str)
    sc.tl.rank_genes_groups(gex_crop_copy, groupby="targeted", method="t-test_overestim_var")
    de = sc.get.rank_genes_groups_df(gex_crop_copy, group="True")
    de = de[de["names"]==gene]
    de["ncells"] = np.sum(gex_crop_copy.obs["targeted"]=="True")
    all_de_crop.append(de)
all_de_crop = pd.concat(all_de_crop)
```

## Code cell 54

```python
all_de["method"] = "Direct capture"
all_de_crop["method"] = "CROP-seq"
de_direct_crop = pd.concat([all_de,all_de_crop])
```

## Code cell 55

```python
de_direct_crop = de_direct_crop[de_direct_crop["names"]!="JUNB"]
de_direct_crop.loc[de_direct_crop["logfoldchanges"]>0,["pvals","pvals_adj"]] = 1
```

## Code cell 56

```python
de_direct_crop[de_direct_crop["method"]=="Direct capture"].sort_values(by="pvals")
```

### Text output (text/plain)

```text
          names    scores  logfoldchanges         pvals     pvals_adj  ncells  \
33251     ITGB1 -9.630595       -3.220815  7.614428e-17  2.531949e-12      86   
33251    CALCRL -7.412618       -3.049972  9.303146e-11  3.093482e-06      48   
33251     ITGA5 -7.007570       -3.693263  2.681223e-09  8.915602e-05      44   
33251      FLI1 -6.254595       -2.725232  1.749041e-08  5.815910e-04      52   
33251      PTK2 -6.178030       -2.438203  4.165887e-08  1.385241e-03      39   
...         ...       ...             ...           ...           ...     ...   
5756      WDR59  0.214375        0.313420  1.000000e+00  1.000000e+00      21   
6779      KDM4B  0.136888        0.151393  1.000000e+00  1.000000e+00      31   
8005   CDC42EP2  0.096128        0.196087  1.000000e+00  1.000000e+00      51   
82         SOX7  1.500038        1.555334  1.000000e+00  1.000000e+00      38   
5954      PDCD7  0.245511        0.551809  1.000000e+00  1.000000e+00      33   

               method  
33251  Direct capture  
33251  Direct capture  
33251  Direct capture  
33251  Direct capture  
33251  Direct capture  
...               ...  
5756   Direct capture  
6779   Direct capture  
8005   Direct capture  
82     Direct capture  
5954   Direct capture  

[242 rows x 7 columns]
```

## Code cell 57

```python
de_direct_crop[de_direct_crop["method"]=="CROP-seq"].sort_values(by="pvals")
```

### Text output (text/plain)

```text
        names    scores  logfoldchanges         pvals  pvals_adj  ncells  \
34000   ITGB1 -6.878323       -1.834114  8.653714e-11   0.000003     100   
34000  CALCRL -7.332406       -3.976672  1.695392e-09   0.000058      37   
34000   ECSCR -5.166112       -5.393016  2.014815e-06   0.046973      71   
34000  IQGAP1 -4.822433       -2.190786  4.344120e-06   0.147704      76   
34000   BMPR2 -4.908681       -2.500965  4.790522e-06   0.162883      54   
...       ...       ...             ...           ...        ...     ...   
8205    SYVN1  0.002484        0.005500  1.000000e+00   1.000000      33   
7813     RAC1  0.049924        0.045391  1.000000e+00   1.000000      32   
1231    USP39  1.019459        1.193625  1.000000e+00   1.000000      78   
6392    BCAR1  0.255573        0.206005  1.000000e+00   1.000000      43   
6065    PDCD7  0.286922        0.633590  1.000000e+00   1.000000      34   

         method  
34000  CROP-seq  
34000  CROP-seq  
34000  CROP-seq  
34000  CROP-seq  
34000  CROP-seq  
...         ...  
8205   CROP-seq  
7813   CROP-seq  
1231   CROP-seq  
6392   CROP-seq  
6065   CROP-seq  

[242 rows x 7 columns]
```

## Code cell 58

```python
p_values_1 = np.sort(de_direct_crop[de_direct_crop["method"]=="Direct capture"]["pvals"])
p_values_2 = np.sort(de_direct_crop[de_direct_crop["method"]=="CROP-seq"]["pvals"])

# Sort the p-values
sorted_p_values_1 = np.sort(p_values_1)
sorted_p_values_2 = np.sort(p_values_2)

# Calculate the theoretical quantiles (expected under the null hypothesis)
theoretical_quantiles = np.linspace(0, 1, max(len(sorted_p_values_1), len(sorted_p_values_2)))

# Transform to uniform distribution
expected_p_values = -np.log10(theoretical_quantiles[1:])  # Avoid log(0)
observed_p_values_1 = -np.log10(sorted_p_values_1[1:])  # Avoid log(0)
observed_p_values_2 = -np.log10(sorted_p_values_2[1:])  # Avoid log(0)

# Create the QQ plot
plt.scatter(expected_p_values, observed_p_values_1, edgecolor='k', linewidth=0.5, facecolor="#ffd479", label="Direct capture", s=100)
plt.scatter(expected_p_values, observed_p_values_2, edgecolor='k', linewidth=0.5, facecolor="#76d6ff", label="CROP-seq", s=100)
plt.plot([min(expected_p_values), max(expected_p_values)],
         [min(expected_p_values), max(expected_p_values)], color="gray", ls='--')
plt.xlabel('Expected -log10(p)')
plt.ylabel('Observed -log10(p)')
plt.title('Target knockdown')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("knockdown_QQ_plot.png",  dpi=300, bbox_inches="tight")
plt.show()
```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Code cell 59

```python
sns.pointplot(y="method", x="logfoldchanges", data=de_direct_crop, join=False)
```

### Text output (stream)

```text
/var/folders/hz/kz7qqv9d2_d4sjp7b0zyz0y00000gq/T/ipykernel_35437/1081230647.py:1: UserWarning: 

The `join` parameter is deprecated and will be removed in v0.15.0. You can remove the line between points with `linestyle='none'`.

  sns.pointplot(y="method", x="logfoldchanges", data=de_direct_crop, join=False)

```

### Text output (text/plain)

```text
<Axes: xlabel='logfoldchanges', ylabel='method'>
```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Code cell 60

```python
"TP53" in de_direct_crop["names"].tolist()
```

### Text output (text/plain)

```text
True
```

## Code cell 61

```python
sns.pointplot(y="method", x="ncells", data=de_direct_crop, join=False)
```

### Text output (stream)

```text
/var/folders/hz/kz7qqv9d2_d4sjp7b0zyz0y00000gq/T/ipykernel_35437/1751036366.py:1: UserWarning: 

The `join` parameter is deprecated and will be removed in v0.15.0. You can remove the line between points with `linestyle='none'`.

  sns.pointplot(y="method", x="ncells", data=de_direct_crop, join=False)

```

### Text output (text/plain)

```text
<Axes: xlabel='ncells', ylabel='method'>
```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Code cell 62

```python

```

## Code cell 63

```python
example_crop = gex_crop_copy.obs.merge(gex_crop_copy[:,"ITGB1"].to_df(),
                                    left_index=True, right_index=True, suffixes=["_assign","_expr"])

sns.violinplot(y="ITGB1_expr",
               x="ITGB1_assign",
               data=example_crop,
               color="#76d6ff")
plt.xlabel("ITGB1 guide")
plt.ylabel("ITGB1 expression")
plt.savefig("example.crop.png", dpi=300, bbox_inches="tight")
```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Code cell 64

```python
example_dc = gex_copy.obs.merge(gex_copy[:,"ITGB1"].to_df(),
                                left_index=True, right_index=True, suffixes=["_assign","_expr"])

sns.violinplot(y="ITGB1_expr",
               x="ITGB1_assign",
               data=example_dc,
               color="#ffd479")
plt.xlabel("ITGB1 guide")
plt.ylabel("ITGB1 expression")
plt.savefig("example.dc.png", dpi=300, bbox_inches="tight")
```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Code cell 65

```python

```

## Code cell 66

```python

```

### Text output (text/plain)

```text
<Axes: xlabel='ITGB1_assign', ylabel='ITGB1_expr'>
```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Code cell 67

```python
pd.concat([gex_crop_copy.obs,
           gex_crop_copy[:,"TP53"].to_df()], axis=1)
```

### Text output (text/plain)

```text
                          PDCD7  MEF2A   UBA5   FDPS    SLK   NEK9  INTS10  \
barcode                                                                      
AAACATCGAACCGAGATTCCGATC  False  False  False  False  False  False   False   
AAACATCGAACTCACCATCCGCGA  False  False  False  False  False  False   False   
AAACATCGAAGAGATCGCTCGCGG  False  False  False  False  False  False   False   
AAACATCGAAGGACACATCCGCGA  False  False  False  False  False  False   False   
AAACATCGAAGGTACACGCTTAAA  False  False  False  False  False  False   False   
...                         ...    ...    ...    ...    ...    ...     ...   
TTCACGCATAGGATGATTGGGAGA   True  False  False  False  False  False   False   
TTCACGCATCCGTCTAGAGCACAA  False  False  False  False  False  False   False   
TTCACGCATCTTCACAGCCGTGTA  False  False  False  False  False  False   False   
TTCACGCATGAAGAGAAGCGAAAC  False   True  False  False  False  False   False   
TTCACGCATGGAACAATCACTTTA  False  False  False  False  False  False   False   

                          GADD45A  PRELID1  METTL3  ...   FLT1  HMGCL  \
barcode                                             ...                 
AAACATCGAACCGAGATTCCGATC    False    False   False  ...  False  False   
AAACATCGAACTCACCATCCGCGA    False    False   False  ...  False  False   
AAACATCGAAGAGATCGCTCGCGG    False    False   False  ...  False  False   
AAACATCGAAGGACACATCCGCGA    False    False   False  ...  False  False   
AAACATCGAAGGTACACGCTTAAA    False    False   False  ...  False  False   
...                           ...      ...     ...  ...    ...    ...   
TTCACGCATAGGATGATTGGGAGA    False    False   False  ...  False  False   
TTCACGCATCCGTCTAGAGCACAA    False    False   False  ...  False  False   
TTCACGCATCTTCACAGCCGTGTA    False    False   False  ...  False  False   
TTCACGCATGAAGAGAAGCGAAAC    False    False   False  ...  False  False   
TTCACGCATGGAACAATCACTTTA    False    False   False  ...  False  False   

                          NSMCE4A  LAMB2    PGF  PECAM1  CDKN1A  YLPM1  \
barcode                                                                  
AAACATCGAACCGAGATTCCGATC    False  False  False   False   False  False   
AAACATCGAACTCACCATCCGCGA    False  False  False   False   False  False   
AAACATCGAAGAGATCGCTCGCGG    False  False  False   False   False  False   
AAACATCGAAGGACACATCCGCGA    False  False  False   False   False  False   
AAACATCGAAGGTACACGCTTAAA    False  False  False   False   False  False   
...                           ...    ...    ...     ...     ...    ...   
TTCACGCATAGGATGATTGGGAGA    False  False  False   False   False  False   
TTCACGCATCCGTCTAGAGCACAA    False  False  False   False   False  False   
TTCACGCATCTTCACAGCCGTGTA    False  False  False   False   False  False   
TTCACGCATGAAGAGAAGCGAAAC    False  False  False   False   False  False   
TTCACGCATGGAACAATCACTTTA    False  False  False   False   False  False   

                          targeted      TP53  
barcode                                       
AAACATCGAACCGAGATTCCGATC     False  0.000000  
AAACATCGAACTCACCATCCGCGA     False  0.000000  
AAACATCGAAGAGATCGCTCGCGG     False  0.000000  
AAACATCGAAGGACACATCCGCGA     False  1.002971  
AAACATCGAAGGTACACGCTTAAA     False  0.000000  
...                            ...       ...  
TTCACGCATAGGATGATTGGGAGA     False  0.000000  
TTCACGCATCCGTCTAGAGCACAA     False  0.000000  
TTCACGCATCTTCACAGCCGTGTA     False  0.000000  
TTCACGCATGAAGAGAAGCGAAAC     False  0.000000  
TTCACGCATGGAACAATCACTTTA     False  0.719439  

[6200 rows x 257 columns]
```

## Code cell 68

```python
gex_crop_copy[:,"TP53"]
```

### Text output (text/plain)

```text
View of AnnData object with n_obs × n_vars = 6200 × 1
    obs: 'PDCD7', 'MEF2A', 'UBA5', 'FDPS', 'SLK', 'NEK9', 'INTS10', 'GADD45A', 'PRELID1', 'METTL3', 'BCAR1', 'SQLE', 'RHOC', 'C18orf21', 'BMP1', 'C17orf58', 'SNAPC2', 'SWAP70', 'ZAK', 'RHOJ', 'MRPS21', 'NFKBIA', 'DTYMK', 'CDC42EP2', 'JUNB', 'NOS3', 'NFKB1', 'MAPK14', 'SRP72', 'IMPDH2', 'WWTR1', 'PRKAR1A', 'RAB24', 'CCDC85C', 'EGFL7', 'LSS', 'MRPL34', 'RHOG', 'CDC42EP5', 'MORF4L1', 'BRD4', 'C7orf26', 'TLN1', 'MYD88', 'INTS1', 'CAMK2G', 'PES1', 'NUP54', 'SCUBE1', 'TMEM2', 'RABGGTA', 'CCL2', 'TAOK1', 'SVIL', 'non-targeting', 'RBM7', 'MAPK1', 'SGK1', 'PMVK', 'SP1', 'ACAT1', 'CDH5', 'FLI1', 'CDC42EP3', 'SMARCC2', 'RAB5C', 'PALLD', 'CLIC1', 'ARFIP2', 'EIF3E', 'RHOB', 'STAT1', 'MAPK11', 'USP39', 'MVD', 'POLR2D', 'CDC42', 'THBS1', 'MEAF6', 'DNTTIP1', 'RAC2', 'HRAS', 'RNF146', 'RPS6KA4', 'CALCRL', 'SLC4A1AP', 'COL4A1', 'ROCK1', 'IBTK', 'MAP2K3', 'VPS39', 'HMGCS1', 'MAPK3', 'CHD4', 'RHOA', 'CDC42EP1', 'PSMD7', 'PRPF3', 'SOX7', 'MAP3K11', 'FOSL1', 'ROCK2', 'ACAT2', 'CDK4', 'DDX6', 'PIP5K1C', 'MICALL2', 'GGT5', 'SMG6', 'HYOU1', 'RAC3', 'CUL3', 'WTAP', 'RHOQ', 'MEF2D', 'PKN1', 'MAPK9', 'SURF4', 'FDFT1', 'WASL', 'GRB2', 'ITGB1BP1', 'USP8', 'MESDC1', 'SRF', 'GDPD5', 'PREX1', 'UPF2', 'KDM2A', 'DNM2', 'KARS', 'NFATC2', 'KDR', 'PDCD10', 'SH2B3', 'BMPR2', 'TP53', 'LDLR', 'EIF3M', 'KIAA1429', 'RPS6KB2', 'PDGFB', 'TSPAN14', 'SRC', 'GNA13', 'ACTG1', 'HSPA5', 'MFF', 'CSNK2B', 'PAK4', 'LRPPRC', 'SMARCB1', 'ZFP36L2', 'RPS6KA2', 'PPAP2B', 'HSPA13', 'KIAA1462', 'RAC1', 'SYVN1', 'XPO7', 'RAP1A', 'KLF2', 'NFE2L2', 'KDM4B', 'DALRD3', 'HIF1A', 'MAP3K3', 'SF3B3', 'HEG1', 'CYFIP1', 'ECSCR', 'CTNNA1', 'TRIP12', 'SERPINH1', 'RCE1', 'HNRNPM', 'ITGB1', 'TEK', 'RPUSD4', 'PNPT1', 'MAP2K7', 'QARS', 'CAMKK2', 'CCM2', 'PHB', 'ITGA5', 'PAK1', 'SRPR', 'CCT8', 'WNK1', 'ARPC2', 'PAK2', 'CLTC', 'HMGCR', 'HMOX1', 'KAT2A', 'TIPARP', 'SMAD3', 'JUND', 'PLRG1', 'FURIN', 'CTNNB1', 'WAC', 'CDC42EP4', 'TFPI', 'PLCG1', 'MAPK7', 'FBN2', 'MVK', 'EDN1', 'DNMT1', 'HDAC7', 'FNTA', 'NUDC', 'AKT1', 'RPS27', 'XPO5', 'NRAS', 'SRCAP', 'RABGGTB', 'MAP2K4', 'SRP68', 'GGPS1', 'BAD', 'GOSR2', 'ID1', 'TERF2IP', 'FOXM1', 'MYC', 'RPS6KA3', 'ELOF1', 'GNA12', 'VEGFB', 'SUZ12', 'PTK2', 'ICAM1', 'KAT5', 'WDR59', 'SMARCAD1', 'COLGALT1', 'HDAC2', 'PIP5K1A', 'JUN', 'IQGAP1', 'CFDP1', 'CCDC92', 'COL4A2', 'FLT1', 'HMGCL', 'NSMCE4A', 'LAMB2', 'PGF', 'PECAM1', 'CDKN1A', 'YLPM1', 'targeted'
    var: 'mt', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
    uns: 'log1p', 'hvg', 'pca', 'neighbors', 'umap', 'rank_genes_groups'
    obsm: 'X_pca', 'X_umap'
    varm: 'PCs'
    layers: 'ambiguous', 'mature', 'nascent', 'counts'
    obsp: 'distances', 'connectivities'
```

## Code cell 69

```python
(gex_crop.obs[guides]>=3).sum(axis=1).mean()
```

### Text output (text/plain)

```text
2.0375806451612903
```

## Code cell 70

```python
sns.histplot((gex.obs[guides]>=5).sum(axis=1).tolist())
```

### Text output (text/plain)

```text
<Axes: ylabel='Count'>
```

### Text output (text/plain)

```text
<Figure size 500x500 with 1 Axes>
```

## Code cell 71

```python
(gex.obs[guides]>=5).sum(axis=1)
```

### Text output (text/plain)

```text
barcode
AAACATCGAACCGAGATTCCGATC    0
AAACATCGAACGCTTAATCGCATA    4
AAACATCGAACGCTTAGGTGGAGC    3
AAACATCGAAGAGATCCGCTTAAA    3
AAACATCGAAGGACACTTGGGAGA    2
                           ..
TTCACGCAGTCTGTCAGAGCACAA    0
TTCACGCATATCAGCACGCGACTA    4
TTCACGCATGAAGAGAACGCCGGC    1
TTCACGCATGAAGAGAGTTCAACA    0
TTCACGCATGAAGAGATGTCTGAA    3
Length: 6001, dtype: int64
```

## Code cell 72

```python

```

