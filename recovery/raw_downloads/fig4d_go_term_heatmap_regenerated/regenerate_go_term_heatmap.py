#!/usr/bin/env python3
"""
Regenerate GO-term x cNMF-program heatmaps for P2, P27, P24, P34.

Input:
  cNMF_50_0_2_032526_6co.xlsx
  Sheet: GO Term Enrichment

Outputs:
  - GO_term_by_program_heatmap_P2_P27_P24_P34.png/pdf
  - GO_term_by_program_heatmap_P2_P27_P24_P34_annotated.png/pdf
  - GO_term_by_program_heatmap_P2_P27_P24_P34_row_normalized.png/pdf
  - source tables

Usage:
  python regenerate_go_term_heatmap.py \
      --cnmf-xlsx /path/to/cNMF_50_0_2_032526_6co.xlsx \
      --outdir /path/to/go_term_program_heatmap_regenerated
"""

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


PROGRAMS = [2, 27, 24, 34]
PROGRAM_LABELS = {
    2: "P2\narterial EC",
    27: "P27\nECM/collagen",
    24: "P24\nHOX/alt context",
    34: "P34\nYAP/TAZ mesenchyme",
}


def clean_term(term):
    """Remove GO ID and shorten long GO term labels for plotting."""
    term = str(term)
    term = re.sub(r"\s*\(GO:\d+\)$", "", term)
    term = term.replace("Regulation Of ", "Reg. of ")
    term = term.replace("Positive Regulation Of ", "Pos. reg. of ")
    term = term.replace("Negative Regulation Of ", "Neg. reg. of ")
    term = term.replace("Extracellular Matrix", "ECM")
    term = term.replace("Extracellular Structure", "Extracellular struct.")
    term = term.replace("External Encapsulating Structure", "External encapsulating struct.")
    term = term.replace("Supramolecular Fiber", "Supramolecular fiber")
    term = term.replace("Cell Migration", "cell migration")
    term = term.replace("Small GTPase Mediated Signal Transduction", "small GTPase signaling")
    term = term.replace("Endothelial Cell", "endothelial cell")
    term = term.replace("Nervous System", "nervous system")
    term = term.replace("Neuron Projection", "neuron projection")
    term = term.replace("Actomyosin Structure", "actomyosin structure")
    term = term.replace("Actin Filament", "actin filament")
    term = term.replace("Collagen Fibril", "collagen fibril")
    term = term.replace("Circulatory System", "circulatory system")
    # Make first character lowercase only if already sentence-like? Keep natural capitalization.
    return term


def load_go(cnmf_xlsx):
    go = pd.read_excel(cnmf_xlsx, sheet_name="GO Term Enrichment")
    required = ["Term", "program_name", "Adjusted P-value", "Genes"]
    missing = [c for c in required if c not in go.columns]
    if missing:
        raise ValueError(f"GO Term Enrichment sheet missing columns: {missing}. Available: {list(go.columns)}")

    go = go.copy()
    go["program_name"] = pd.to_numeric(go["program_name"], errors="coerce").astype("Int64")
    go["Adjusted P-value"] = pd.to_numeric(go["Adjusted P-value"], errors="coerce")
    go = go.dropna(subset=["program_name", "Adjusted P-value", "Term"])
    go["program_name"] = go["program_name"].astype(int)
    go["neglog10_FDR"] = -np.log10(go["Adjusted P-value"].clip(lower=1e-300))
    go["Term_short"] = go["Term"].map(clean_term)
    return go


def select_top_terms(go, programs=PROGRAMS, top_n_per_program=5):
    selected = []
    seen = set()
    for p in programs:
        sub = go[go["program_name"].eq(p)].sort_values("Adjusted P-value", ascending=True).head(top_n_per_program)
        for term in sub["Term"].tolist():
            if term not in seen:
                selected.append(term)
                seen.add(term)
    return selected


def make_matrices(go, selected_terms, programs=PROGRAMS):
    # Actual -log10(FDR)
    mat = pd.DataFrame(index=selected_terms, columns=programs, dtype=float)
    raw_fdr = pd.DataFrame(index=selected_terms, columns=programs, dtype=float)
    genes_mat = pd.DataFrame(index=selected_terms, columns=programs, dtype=object)

    for p in programs:
        sub = go[go["program_name"].eq(p)].copy()
        sub = sub.drop_duplicates("Term", keep="first").set_index("Term")
        for term in selected_terms:
            if term in sub.index:
                mat.loc[term, p] = sub.loc[term, "neglog10_FDR"]
                raw_fdr.loc[term, p] = sub.loc[term, "Adjusted P-value"]
                genes_mat.loc[term, p] = sub.loc[term, "Genes"]
            else:
                mat.loc[term, p] = 0.0
                raw_fdr.loc[term, p] = np.nan
                genes_mat.loc[term, p] = ""

    # Row-normalized specificity matrix
    row_max = mat.max(axis=1).replace(0, np.nan)
    mat_rownorm = mat.div(row_max, axis=0).fillna(0)

    # Plot labels
    term_short = pd.Series({term: clean_term(term) for term in selected_terms})

    return mat, raw_fdr, genes_mat, mat_rownorm, term_short


def plot_heatmap(mat, term_short, outdir, stem, title, colorbar_label, vmax=None, annotate=False, row_norm=False):
    plot_mat = mat.copy()
    labels_y = [term_short.loc[t] for t in plot_mat.index]
    labels_x = [PROGRAM_LABELS.get(int(c), f"P{c}") for c in plot_mat.columns]

    data = plot_mat.values.astype(float)
    if vmax is None:
        vmax = np.nanpercentile(data[data > 0], 95) if np.any(data > 0) else 1.0
        vmax = max(vmax, 1.0)
    vmin = 0

    height = max(6.0, 0.34 * len(labels_y) + 1.5)
    width = max(5.4, 1.05 * len(labels_x) + 3.8)

    fig, ax = plt.subplots(figsize=(width, height))
    im = ax.imshow(data, aspect="auto", cmap="Reds", vmin=vmin, vmax=vmax)

    ax.set_xticks(np.arange(len(labels_x)))
    ax.set_xticklabels(labels_x, fontsize=10)
    ax.set_yticks(np.arange(len(labels_y)))
    ax.set_yticklabels(labels_y, fontsize=8)
    ax.set_title(title, fontsize=12, fontweight="bold")

    if annotate:
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                val = data[i, j]
                if val > 0:
                    if row_norm:
                        txt = f"{val:.2f}"
                    else:
                        txt = f"{val:.1f}"
                    ax.text(j, i, txt, ha="center", va="center", fontsize=6.5, color="black")

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    ax.set_xlabel("cNMF program", fontsize=10, fontweight="bold")
    ax.set_ylabel("GO term", fontsize=10, fontweight="bold")

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label(colorbar_label, fontsize=9)

    fig.tight_layout()
    png = Path(outdir) / f"{stem}.png"
    pdf = Path(outdir) / f"{stem}.pdf"
    fig.savefig(png, dpi=350, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cnmf-xlsx", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--top-n", type=int, default=5, help="Top GO terms per program to include")
    parser.add_argument("--vmax", type=float, default=None, help="Color scale max for -log10(FDR), default 95th percentile")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    go = load_go(args.cnmf_xlsx)
    selected_terms = select_top_terms(go, programs=PROGRAMS, top_n_per_program=args.top_n)
    mat, raw_fdr, genes_mat, mat_rownorm, term_short = make_matrices(go, selected_terms, programs=PROGRAMS)

    # Source outputs
    go[go["program_name"].isin(PROGRAMS)].to_csv(outdir / "source_GO_term_enrichment_selected_programs_all_terms.csv", index=False)
    pd.DataFrame({"Term": selected_terms, "Term_short": [term_short.loc[t] for t in selected_terms]}).to_csv(outdir / "source_selected_GO_terms.csv", index=False)
    mat.to_csv(outdir / "source_GO_term_program_neglog10FDR_matrix.csv")
    raw_fdr.to_csv(outdir / "source_GO_term_program_FDR_matrix.csv")
    genes_mat.to_csv(outdir / "source_GO_term_program_overlap_genes_matrix.csv")
    mat_rownorm.to_csv(outdir / "source_GO_term_program_row_normalized_matrix.csv")

    long = mat.reset_index().rename(columns={"index": "Term"}).melt(id_vars="Term", var_name="program", value_name="neglog10_FDR")
    long = long.merge(raw_fdr.reset_index().rename(columns={"index":"Term"}).melt(id_vars="Term", var_name="program", value_name="FDR"), on=["Term","program"], how="left")
    long = long.merge(genes_mat.reset_index().rename(columns={"index":"Term"}).melt(id_vars="Term", var_name="program", value_name="Genes"), on=["Term","program"], how="left")
    long["Term_short"] = long["Term"].map(term_short)
    long["program_label"] = long["program"].map(lambda x: PROGRAM_LABELS.get(int(x), f"P{x}"))
    long.to_csv(outdir / "source_GO_term_program_support_long.csv", index=False)

    # Plots
    plot_heatmap(
        mat,
        term_short,
        outdir,
        stem="GO_term_by_program_heatmap_P2_P27_P24_P34",
        title="GO term enrichment across selected cNMF programs",
        colorbar_label=r"$-log_{10}$(FDR)",
        vmax=args.vmax,
        annotate=False,
        row_norm=False,
    )

    plot_heatmap(
        mat,
        term_short,
        outdir,
        stem="GO_term_by_program_heatmap_P2_P27_P24_P34_annotated",
        title="GO term enrichment across selected cNMF programs",
        colorbar_label=r"$-log_{10}$(FDR)",
        vmax=args.vmax,
        annotate=True,
        row_norm=False,
    )

    plot_heatmap(
        mat_rownorm,
        term_short,
        outdir,
        stem="GO_term_by_program_heatmap_P2_P27_P24_P34_row_normalized",
        title="GO term enrichment specificity across selected programs",
        colorbar_label="row-normalized\nGO support",
        vmax=1.0,
        annotate=False,
        row_norm=True,
    )

    with open(outdir / "README.txt", "w") as f:
        f.write(
            "GO term x cNMF program heatmap regenerated.\n\n"
            f"Input workbook: {args.cnmf_xlsx}\n"
            "Source sheet: GO Term Enrichment\n"
            "Rows: union of top GO terms per selected program by adjusted p-value.\n"
            "Columns: P2, P27, P24, P34.\n"
            "Main value: -log10(adjusted p-value).\n"
            "Row-normalized version emphasizes program specificity rather than absolute enrichment strength.\n"
        )

    print("Wrote outputs to", outdir)
    print("Selected GO terms:")
    for term in selected_terms:
        print(" -", clean_term(term))


if __name__ == "__main__":
    main()
