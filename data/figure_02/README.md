# Figure 2 cNMF workbook

Panels 2a and 2b read the original supplied cNMF workbook from this repository-relative path:

```text
data/figure_02/cNMF_50_0_2_032526_6co.xlsx
```

The expected SHA-256 is:

```text
d2fe4834b1efd3c00e01062b56b907a6178bfbc48a40beb5ea7b347538bfc824
```

The notebook reads `Perturbation Association 1` directly. For panel 2a it forms a 298 × 200 matrix from D0–D3 and P0–P49, drops the documented duplicate alias `NEDD8-MDP1`, and retains canonical `NEDD8` plus the pooled `non-targeting` profile. Panel 2b reads `approx_log2FC` and `adj_pval` for the displayed genes and programs.

The panel-b EC-output column is recalculated from the integrated AnnData. It uses all 120 controls (60 non-targeting plus 60 safe-targeting), two-sided Welch t-tests at the guide level, and BH correction across all 297 target genes.

Panels 2a, 2b, 2d and 2e are rendered with plotnine. Panel 2f retains its existing Matplotlib renderer.

The connected GitHub contents API can write UTF-8 source files but cannot ingest the local binary XLSX by path. Install and push the verified workbook from Sherlock with:

```bash
bash scripts/install_figure2_cnmf_workbook.sh --push
```

The installer verifies the source hash before copying, committing and pushing the workbook. At approximately 14 MB, the file is below GitHub's 100 MB ordinary-blob limit.
