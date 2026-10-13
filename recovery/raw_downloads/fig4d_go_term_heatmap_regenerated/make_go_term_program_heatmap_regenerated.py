#!/usr/bin/env python3
"""
Make GO-term x cNMF-program heatmap for selected programs.

Default programs:
    P2, P27, P24, P34

Input:
    cNMF workbook with sheet `GO Term Enrichment`

Outputs:
    GO_term_by_program_heatmap_P2_P27_P24_P34.png/pdf
    GO_term_by_program_heatmap_P2_P27_P24_P34_annotated.png/pdf
    source_GO_term_program_neglog10FDR_matrix.csv
    source_GO_term_program_GO_enrichment_long.csv
    source_GO_term_program_selected_terms.csv

Run:
    python make_go_term_program_heatmap_regenerated.py \
        --cnmf-xlsx /path/to/cNMF_50_0_2_032526_6co.xlsx \
        --outdir /path/to/go_term_program_heatmap_regenerated
"""

import argparse
import os
from pathlib import Path
import textwrap
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PROGRAMS = [2, 27, 24, 34]
PROGRAM_LABELS = {2: "P2", 27: "P27", 24: "P24", 34: "P34"}
PROGRAM_ANNOTATIONS = {
    2: "arterial EC",
    27: "ECM alt.",
    24: "patterning alt.",
    34: "YAP/TAZ mesenchyme",
}


def strip_go_id(term):
    term = str(term)
    term = re.sub(r"\s*\(GO:\d+\)\s*$", "", term)
    return term


def wrap_label(term, width=42):
    term = strip_go_id(term)
    # Keep common abbreviations readable.
    term = term.replace("Regulation Of", "Reg. of")
    term = term.replace("Positive Regulation Of", "Pos. reg. of")
    term = term.replace("Negative Regulation Of", "Neg. reg. of")
    term = term.replace("Extracellular Matrix", "ECM")
    term = term.replace("Extracellular Structure", "ECM")
    term = term.replace("Cell-Cell", "Cell-cell")
    term = term.replace("Cell-Matrix", "Cell-matrix")
    return "\n".join(textwrap.wrap(term, width=width))


def load_go(cnmf_xlsx):
    go = pd.read_excel(cnmf_xlsx, sheet_name="GO Term Enrichment")
    required = ["Term", "program_name", "Adjusted P-value", "Combined Score", "Odds Ratio", "Genes"]
    missing = [c for c in required if c not in go.columns]
    if missing:
        raise ValueError(f"GO Term Enrichment sheet missing columns: {missing}; available={list(go.columns)}")
    go = go.copy()
    go["program_name"] = pd.to_numeric(go["program_name"], errors="coerce").astype("Int64")
    go["Adjusted P-value"] = pd.to_numeric(go["Adjusted P-value"], errors="coerce")
    go["Combined Score"] = pd.to_numeric(go["Combined Score"], errors="coerce")
    go["Odds Ratio"] = pd.to_numeric(go["Odds Ratio"], errors="coerce")
    go["neglog10_FDR"] = -np.log10(go["Adjusted P-value"].clip(lower=1e-300))
    go["Term_noGO"] = go["Term"].map(strip_go_id)
    return go


def select_terms(go, programs, top_n_per_program=5, max_terms=None):
    selected = []
    for p in programs:
        sub = go[go["program_name"].eq(p)].sort_values(["Adjusted P-value", "Combined Score"], ascending=[True, False]).copy()
        selected.extend(sub.head(top_n_per_program)["Term"].tolist())
    # unique preserving order
    terms = []
    for t in selected:
        if t not in terms:
            terms.append(t)
    if max_terms is not None:
        terms = terms[:max_terms]
    return terms


def build_matrix(go, programs, terms):
    sub = go[go["program_name"].isin(programs) & go["Term"].isin(terms)].copy()
    mat = sub.pivot_table(index="Term", columns="program_name", values="neglog10_FDR", aggfunc="max")
    mat = mat.reindex(index=terms, columns=programs)
    mat = mat.fillna(0.0)

    fdr = sub.pivot_table(index="Term", columns="program_name", values="Adjusted P-value", aggfunc="min")
    fdr = fdr.reindex(index=terms, columns=programs)

    genes = sub.pivot_table(index="Term", columns="program_name", values="Genes", aggfunc=lambda x: ";".join(map(str, x)))
    genes = genes.reindex(index=terms, columns=programs)

    # Row order: group by program with maximum score, following program order; then score descending.
    order_df = pd.DataFrame({"Term": mat.index, "max_score": mat.max(axis=1), "argmax_program": mat.idxmax(axis=1)})
    order_df["program_order"] = order_df["argmax_program"].map({p: i for i, p in enumerate(programs)})
    order_df = order_df.sort_values(["program_order", "max_score"], ascending=[True, False])
    mat = mat.loc[order_df["Term"]]
    fdr = fdr.loc[order_df["Term"]]
    genes = genes.loc[order_df["Term"]]

    return mat, fdr, genes, sub


def save_matrix_sources(mat, fdr, genes, go_long, outdir):
    mat_out = mat.copy()
    mat_out.columns = [PROGRAM_LABELS.get(c, str(c)) for c in mat_out.columns]
    mat_out.index.name = "GO_term"
    mat_out.to_csv(outdir / "source_GO_term_program_neglog10FDR_matrix.csv")

    fdr_out = fdr.copy()
    fdr_out.columns = [PROGRAM_LABELS.get(c, str(c)) for c in fdr_out.columns]
    fdr_out.index.name = "GO_term"
    fdr_out.to_csv(outdir / "source_GO_term_program_FDR_matrix.csv")

    genes_out = genes.copy()
    genes_out.columns = [PROGRAM_LABELS.get(c, str(c)) for c in genes_out.columns]
    genes_out.index.name = "GO_term"
    genes_out.to_csv(outdir / "source_GO_term_program_overlap_genes_matrix.csv")

    go_long.to_csv(outdir / "source_GO_term_program_GO_enrichment_long.csv", index=False)

    selected_terms = pd.DataFrame({
        "GO_term": mat.index,
        "GO_term_short": [strip_go_id(x) for x in mat.index],
        "dominant_program": [PROGRAM_LABELS.get(p, str(p)) for p in mat.idxmax(axis=1)],
        "max_neglog10FDR": mat.max(axis=1).values,
    })
    selected_terms.to_csv(outdir / "source_GO_term_program_selected_terms.csv", index=False)


def plot_heatmap(mat, outdir, annotated=False, cap=10.0):
    # Display matrix with rows in reverse so first group is at top with imshow origin upper.
    display = mat.copy()
    vals = display.values.astype(float)
    vals_capped = np.clip(vals, 0, cap)

    row_labels = [wrap_label(t, width=42) for t in display.index]
    col_labels = [f"{PROGRAM_LABELS.get(c, c)}\n{PROGRAM_ANNOTATIONS.get(c, '')}" for c in display.columns]

    n_rows = vals.shape[0]
    n_cols = vals.shape[1]
    fig_h = max(5.2, 0.37 * n_rows + 1.8)
    fig_w = max(4.8, 1.0 * n_cols + 4.6)
    if annotated:
        fig_w += 0.8

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(vals_capped, aspect="auto", cmap="Reds", vmin=0, vmax=cap)

    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(col_labels, rotation=0, ha="center", fontsize=10)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(row_labels, fontsize=8.5)
    ax.set_xlabel("cNMF program", fontsize=11, fontweight="bold")
    ax.set_ylabel("GO term", fontsize=11, fontweight="bold")
    ax.set_title("GO term enrichment across selected D3 programs", fontsize=12, fontweight="bold", pad=16)

    # Gridlines for figure-d-like heatmap
    ax.set_xticks(np.arange(-.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if annotated:
        for i in range(n_rows):
            for j in range(n_cols):
                v = vals[i, j]
                if np.isfinite(v) and v > 0:
                    label = f"{v:.1f}"
                    color = "white" if v >= cap * 0.65 else "black"
                    ax.text(j, i, label, ha="center", va="center", fontsize=8, color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("−log10(FDR)\n(capped)", fontsize=10)

    fig.tight_layout()
    stem = "GO_term_by_program_heatmap_P2_P27_P24_P34_annotated" if annotated else "GO_term_by_program_heatmap_P2_P27_P24_P34"
    png = outdir / f"{stem}.png"
    pdf = outdir / f"{stem}.pdf"
    fig.savefig(png, dpi=350, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cnmf-xlsx", default="/mnt/data/cNMF_50_0_2_032526_6co.xlsx")
    parser.add_argument("--outdir", default="/mnt/data/go_term_program_heatmap_regenerated")
    parser.add_argument("--programs", default="2,27,24,34")
    parser.add_argument("--top-n-per-program", type=int, default=5)
    parser.add_argument("--cap-neglog10-fdr", type=float, default=10.0)
    args = parser.parse_args()

    cnmf = args.cnmf_xlsx
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    programs = [int(x.strip()) for x in args.programs.split(",") if x.strip()]

    if not os.path.exists(cnmf):
        raise FileNotFoundError(cnmf)

    go = load_go(cnmf)
    terms = select_terms(go, programs, top_n_per_program=args.top_n_per_program)
    mat, fdr, genes, go_long = build_matrix(go, programs, terms)
    save_matrix_sources(mat, fdr, genes, go_long, outdir)
    plot_heatmap(mat, outdir, annotated=False, cap=args.cap_neglog10_fdr)
    plot_heatmap(mat, outdir, annotated=True, cap=args.cap_neglog10_fdr)

    with open(outdir / "README.txt", "w") as f:
        f.write(
            "GO term x cNMF program heatmap regenerated from the cNMF workbook.\n\n"
            f"Input workbook: {cnmf}\n"
            "Source sheet: GO Term Enrichment\n"
            f"Programs: {programs}\n"
            f"Top terms per program: {args.top_n_per_program}\n"
            "Matrix value: -log10(adjusted p-value), capped for visualization.\n\n"
            "Outputs:\n"
            "- GO_term_by_program_heatmap_P2_P27_P24_P34.png/pdf\n"
            "- GO_term_by_program_heatmap_P2_P27_P24_P34_annotated.png/pdf\n"
            "- source_GO_term_program_neglog10FDR_matrix.csv\n"
            "- source_GO_term_program_FDR_matrix.csv\n"
            "- source_GO_term_program_overlap_genes_matrix.csv\n"
            "- source_GO_term_program_GO_enrichment_long.csv\n"
        )

    print("Wrote outputs to", outdir)
    print("Selected GO terms:")
    print(pd.read_csv(outdir / "source_GO_term_program_selected_terms.csv").to_string(index=False))


if __name__ == "__main__":
    main()
