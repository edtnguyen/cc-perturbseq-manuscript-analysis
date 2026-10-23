# Figure 2 cNMF workbook

Panels 2a and 2b read the original supplied cNMF workbook from this relative path:

```text
data/figure_02/cNMF_50_0_2_032526_6co.xlsx
```

The expected SHA-256 is:

```text
d2fe4834b1efd3c00e01062b56b907a6178bfbc48a40beb5ea7b347538bfc824
```

The notebook reads `Perturbation Association 1` directly. For panel 2a it forms a 298 × 200 matrix from D0–D3 and P0–P49, drops the documented duplicate alias `NEDD8-MDP1`, and retains canonical `NEDD8` plus the pooled `non-targeting` profile. Panel 2b reads `approx_log2FC` and `adj_pval` for the displayed genes and programs.

Because the connected GitHub API does not accept local binary workbook bytes, install the workbook from Sherlock with:

```bash
bash scripts/install_figure2_cnmf_workbook.sh
```

The installer verifies the source hash before copying. Then commit the workbook normally; at approximately 14 MB it is below GitHub's 100 MB ordinary-blob limit.
