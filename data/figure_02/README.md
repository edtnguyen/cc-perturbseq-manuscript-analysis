# Figure 2 cNMF workbook

Panels 2a and 2b read the original supplied cNMF workbook from this
repository-relative path:

```text
data/figure_02/cNMF_50_0_2_032526_6co.xlsx
```

The expected SHA-256 is:

```text
d2fe4834b1efd3c00e01062b56b907a6178bfbc48a40beb5ea7b347538bfc824
```

The notebook reads `Perturbation Association 1` directly. For panel 2a it
forms a 298 × 200 matrix from D0-D3 and P0-P49, collapses the documented
duplicate alias `NEDD8-MDP1` into canonical `NEDD8`, and retains the pooled
`non-targeting` profile. Panel 2b reads `approx_log2FC` and `adj_pval` for the
displayed genes and programs.

The panel-b EC-output column is recalculated from the integrated AnnData. It
uses all 120 controls (60 non-targeting plus 60 safe-targeting), two-sided
Welch t-tests at the guide level, and BH correction across all 297 target
genes.

Panel 2d also uses all 120 negative-control guides. For each control guide,
cells carrying any KDR, PLCG1, PRKCD, or ETV2 guide are excluded before its D3
EC percentage is calculated. Each displayed target is compared with the
120-guide control distribution by a two-sided Welch t-test, followed by BH
correction across the four displayed genes.

Panels 2a, 2b, 2d, 2e, and 2f are rendered with plotnine. Panel 2f directly
ports the recovered carrier-seeded plotnine notebook:

```text
programDE/notebooks/fig3_carrierseeded_neighborhood_plotnine_style_notebook_package/
fig3_carrierseeded_neighborhood_panels_interactive_plotnine_style.ipynb
```

It retains the same-day 300-neighbor calculation, 0.5 pseudocount, all-carrier
seed definition, pooled 120-guide controls, red/white/blue enrichment scale,
and faceted plotnine UMAP style.

The connected GitHub text-file API cannot ingest a local binary XLSX by path.
Install and push the verified workbook from Sherlock with:

```bash
bash scripts/install_figure2_cnmf_workbook.sh --push
```

The installer verifies the source hash before copying, committing, and pushing
the workbook. At approximately 14 MB, the file is below GitHub's 100 MB
ordinary-blob limit.
