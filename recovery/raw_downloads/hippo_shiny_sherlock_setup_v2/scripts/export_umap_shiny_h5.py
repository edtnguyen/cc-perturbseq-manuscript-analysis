#!/usr/bin/env python
"""
Export a compact HDF5 bundle for a Shiny UMAP viewer.

The Shiny app can then display either:
  1. cNMF program score on UMAP
  2. selected gene expression on UMAP

Design:
- UMAP is taken from the integrated AnnData object by default.
- cNMF scores are taken from the cNMF modality of the h5mu file.
- Cells are aligned by obs_names before export.
- Gene expression is stored as gene-by-cell CSR arrays in HDF5, so the Shiny app
  retrieves only the selected gene vector instead of loading the full matrix.

Typical usage on Sherlock:

python export_umap_shiny_h5.py \
  --repo-root /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE \
  --gene-mode all

For a smaller first test:

python export_umap_shiny_h5.py \
  --repo-root /oak/stanford/groups/engreitz/Users/tri/CodeBase/programDE \
  --gene-mode list \
  --gene-list markers_for_shiny.txt
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
from scipy import sparse
import matplotlib as mpl
import matplotlib.pyplot as plt


DEFAULT_MDATA_PATH = (
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/adata/cNMF_50_0_2.h5mu"
)

DEFAULT_ADATA_INT_PATH = (
    "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
    "2.preprocessing_tonys_pipeline/3.integrated_days/"
    "3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad"
)


def find_repo_root(start: Path | None = None) -> Path:
    """Find the project root by looking for data/ and analysis/."""
    if start is None:
        start = Path.cwd()
    start = start.resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "data").exists() and (candidate / "analysis").exists():
            return candidate
    raise FileNotFoundError(
        "Could not auto-detect repo root. Pass --repo-root explicitly."
    )


def parse_program_order(path: Path, deduplicate: bool = True) -> list[int]:
    """Parse program_order.txt. Accepts lines like 15, P15, Program 15, etc."""
    programs: list[int] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        # Strip inline comments.
        line = line.split("#", 1)[0].strip()
        match = re.search(r"(?:Program\s*)?P?(\d+)", line, flags=re.IGNORECASE)
        if match is None:
            raise ValueError(f"Could not parse program from line: {raw!r}")
        programs.append(int(match.group(1)))

    if deduplicate:
        programs = list(dict.fromkeys(programs))
    return programs


def read_gene_list(path: Path) -> list[str]:
    genes: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        genes.append(line.split()[0])
    return list(dict.fromkeys(genes))


def make_label_and_color_maps(prog_meta: pd.DataFrame) -> tuple[dict[int, str], dict[int, str]]:
    if "program_name" not in prog_meta.columns:
        raise KeyError("program_activity_long_v4.csv must contain a `program_name` column.")

    prog_meta = prog_meta.drop_duplicates("program_name").copy()
    program_ids = prog_meta["program_name"].astype(int)

    if "program_label" in prog_meta.columns:
        label_map = dict(zip(program_ids, prog_meta["program_label"].astype(str)))
    else:
        label_map = {int(program): f"P{int(program)}" for program in program_ids}

    preferred_color_columns = [
        "distinct_color",
        "program_color",
        "color",
        "class_color",
    ]
    color_column = next((col for col in preferred_color_columns if col in prog_meta.columns), None)

    prog_color: dict[int, str] = {}
    if color_column is not None:
        for program, color in zip(program_ids, prog_meta[color_column]):
            if pd.notna(color) and str(color).strip():
                prog_color[int(program)] = str(color).strip()
        print(f"Using `{color_column}` from program metadata for program colors")
    else:
        print(
            "No program color column found. Generating fallback colors. "
            f"Available columns: {list(prog_meta.columns)}"
        )

    fallback_cmap = plt.get_cmap("tab20")
    for idx, program in enumerate(program_ids.astype(int)):
        if int(program) not in prog_color:
            prog_color[int(program)] = mpl.colors.to_hex(fallback_cmap(idx % fallback_cmap.N))

    return label_map, prog_color


def get_expr_matrix_and_genes(adata: sc.AnnData, source: str):
    """
    Return a cells x genes matrix and gene names.

    source options:
      - X
      - raw
      - layer:<layer_name>
    """
    source = source.strip()
    if source == "X":
        return adata.X, np.asarray(adata.var_names.astype(str))
    if source == "raw":
        if adata.raw is None:
            raise ValueError("--expression-source raw was requested, but adata.raw is None.")
        return adata.raw.X, np.asarray(adata.raw.var_names.astype(str))
    if source.startswith("layer:"):
        layer_name = source.split(":", 1)[1]
        if layer_name not in adata.layers:
            raise KeyError(
                f"Layer {layer_name!r} not found. Available layers: {list(adata.layers.keys())}"
            )
        return adata.layers[layer_name], np.asarray(adata.var_names.astype(str))
    raise ValueError("--expression-source must be one of: X, raw, layer:<name>")


def as_dense_array(x):
    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


def write_string_dataset(group, name: str, values):
    str_dtype = h5py.string_dtype(encoding="utf-8")
    arr = np.asarray(["" if v is None else str(v) for v in values], dtype=object)
    group.create_dataset(name, data=arr, dtype=str_dtype, compression="gzip")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=None)
    parser.add_argument("--mdata", type=Path, default=Path(DEFAULT_MDATA_PATH))
    parser.add_argument("--adata-int", type=Path, default=Path(DEFAULT_ADATA_INT_PATH))
    parser.add_argument("--program-meta", type=Path, default=None)
    parser.add_argument("--program-order", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument(
        "--gene-mode",
        choices=["all", "list"],
        default="all",
        help="Export all genes or only genes listed in --gene-list.",
    )
    parser.add_argument("--gene-list", type=Path, default=None)
    parser.add_argument(
        "--expression-source",
        default="X",
        help="Gene expression source: X, raw, or layer:<layer_name>.",
    )
    parser.add_argument(
        "--obs-columns",
        default="day,cell_type",
        help="Comma-separated obs columns to export if present.",
    )
    parser.add_argument(
        "--allow-extra-adata-cells",
        action="store_true",
        help="Allow adata_int to contain extra cells; it will be subset/reordered to cNMF cells.",
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve() if args.repo_root is not None else find_repo_root()
    program_meta_path = args.program_meta or repo_root / "data" / "program_activity_long_v4.csv"
    program_order_path = args.program_order or repo_root / "analysis" / "regulator_dotplots" / "program_order.txt"
    out_path = args.out or repo_root / "analysis" / "regulator_dotplots" / "shiny_umap_app" / "umap_gene_program_data.h5"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Repo root:       {repo_root}")
    print(f"cNMF h5mu:       {args.mdata}")
    print(f"Integrated h5ad: {args.adata_int}")
    print(f"Program meta:    {program_meta_path}")
    print(f"Program order:   {program_order_path}")
    print(f"Output HDF5:     {out_path}")

    print("Loading cNMF h5mu...")
    mdata = mu.read_h5mu(args.mdata)
    if "cNMF" not in mdata.mod:
        raise KeyError(f"Expected modality `cNMF`. Available modalities: {list(mdata.mod.keys())}")
    cnmf = mdata.mod["cNMF"]
    cnmf_cells = pd.Index(cnmf.obs_names.astype(str))
    if not cnmf_cells.is_unique:
        raise ValueError("cnmf.obs_names are not unique.")

    print("Loading integrated AnnData into memory...")
    adata_int = sc.read_h5ad(args.adata_int)
    adata_cells = pd.Index(adata_int.obs_names.astype(str))
    if not adata_cells.is_unique:
        raise ValueError("adata_int.obs_names are not unique.")

    missing_from_adata = cnmf_cells.difference(adata_cells)
    extra_in_adata = adata_cells.difference(cnmf_cells)
    if len(missing_from_adata) > 0:
        raise ValueError(
            f"adata_int is missing {len(missing_from_adata)} cNMF cells. "
            f"First few: {missing_from_adata[:5].tolist()}"
        )
    if len(extra_in_adata) > 0 and not args.allow_extra_adata_cells:
        raise ValueError(
            f"adata_int has {len(extra_in_adata)} extra cells not present in cNMF. "
            "If this is expected, rerun with --allow-extra-adata-cells. "
            f"First few extras: {extra_in_adata[:5].tolist()}"
        )

    if not adata_cells.equals(cnmf_cells):
        print("Reordering/subsetting adata_int to match cnmf.obs_names...")
        adata_int = adata_int[cnmf_cells.tolist(), :].copy()
    else:
        print("Cell order already matches between cNMF and adata_int.")

    if "X_umap" not in adata_int.obsm:
        raise KeyError(f"adata_int.obsm does not contain X_umap. Available keys: {list(adata_int.obsm.keys())}")
    umap = np.asarray(adata_int.obsm["X_umap"])
    if umap.shape[0] != cnmf.n_obs or umap.shape[1] < 2:
        raise ValueError(f"Unexpected X_umap shape after alignment: {umap.shape}")

    print("Reading program order and metadata...")
    programs = parse_program_order(program_order_path, deduplicate=True)
    prog_meta = pd.read_csv(program_meta_path)
    label_map, color_map = make_label_and_color_maps(prog_meta)

    score_matrix = cnmf.X
    if sparse.issparse(score_matrix):
        score_matrix = score_matrix.toarray()
    else:
        score_matrix = np.asarray(score_matrix)
    if score_matrix.shape[0] != len(cnmf_cells):
        raise ValueError("cNMF score matrix row count does not match cnmf.obs_names.")
    invalid_programs = [p for p in programs if p < 0 or p >= score_matrix.shape[1]]
    if invalid_programs:
        raise ValueError(
            f"These program IDs are outside the cNMF score matrix column range: {invalid_programs}. "
            f"Score matrix has {score_matrix.shape[1]} columns."
        )
    program_scores = score_matrix[:, programs].T.astype(np.float32)
    program_names = [f"P{p}" for p in programs]
    program_labels = [label_map.get(p, f"P{p}") for p in programs]
    program_colors = [color_map.get(p, "#999999") for p in programs]

    print("Preparing gene expression matrix...")
    X, gene_names_all = get_expr_matrix_and_genes(adata_int, args.expression_source)
    gene_index = pd.Index(gene_names_all.astype(str))
    if not gene_index.is_unique:
        raise ValueError(
            "Gene names are not unique. Make them unique before export, or export a gene list with unique names."
        )

    if args.gene_mode == "list":
        if args.gene_list is None:
            raise ValueError("--gene-mode list requires --gene-list.")
        requested_genes = read_gene_list(args.gene_list)
        missing_genes = [g for g in requested_genes if g not in gene_index]
        if missing_genes:
            print(
                f"Warning: {len(missing_genes)} requested genes are missing and will be skipped. "
                f"First few: {missing_genes[:10]}",
                file=sys.stderr,
            )
        genes_to_export = [g for g in requested_genes if g in gene_index]
        if not genes_to_export:
            raise ValueError("No requested genes were found in adata_int.")
        gene_idx = gene_index.get_indexer(genes_to_export)
        X = X[:, gene_idx]
        gene_names = np.asarray(genes_to_export, dtype=str)
    else:
        gene_names = gene_names_all.astype(str)

    if sparse.issparse(X):
        X_cells_by_genes = X.tocsr()
    else:
        # This can be large for dense expression matrices. Most scRNA-seq matrices should be sparse.
        X_cells_by_genes = sparse.csr_matrix(np.asarray(X))

    print(f"Expression matrix for export: {X_cells_by_genes.shape[0]} cells x {X_cells_by_genes.shape[1]} genes")
    print("Transposing to gene-by-cell CSR for fast per-gene retrieval in Shiny...")
    X_genes_by_cells = X_cells_by_genes.T.tocsr()
    X_genes_by_cells.sort_indices()

    obs_columns_requested = [c.strip() for c in args.obs_columns.split(",") if c.strip()]
    obs_columns = [c for c in obs_columns_requested if c in adata_int.obs.columns]
    missing_obs_cols = [c for c in obs_columns_requested if c not in adata_int.obs.columns]
    if missing_obs_cols:
        print(f"Skipping missing obs columns: {missing_obs_cols}")

    print("Writing HDF5 bundle...")
    with h5py.File(out_path, "w") as h5:
        h5.attrs["description"] = "UMAP + cNMF program scores + sparse gene expression for Shiny viewer"
        h5.attrs["mdata_path"] = str(args.mdata)
        h5.attrs["adata_int_path"] = str(args.adata_int)
        h5.attrs["expression_source"] = args.expression_source

        obs_grp = h5.create_group("obs")
        write_string_dataset(obs_grp, "cell_id", cnmf_cells.tolist())
        obs_grp.create_dataset("umap_1", data=umap[:, 0].astype(np.float32), compression="gzip")
        obs_grp.create_dataset("umap_2", data=umap[:, 1].astype(np.float32), compression="gzip")

        ann_grp = h5.create_group("obs_annotations")
        write_string_dataset(h5, "obs_annotation_names", obs_columns)
        for col in obs_columns:
            write_string_dataset(ann_grp, col, adata_int.obs[col].astype(str).fillna("NA").tolist())

        prog_grp = h5.create_group("programs")
        prog_grp.create_dataset("ids", data=np.asarray(programs, dtype=np.int32), compression="gzip")
        write_string_dataset(prog_grp, "names", program_names)
        write_string_dataset(prog_grp, "labels", program_labels)
        write_string_dataset(prog_grp, "colors", program_colors)
        prog_grp.create_dataset("scores", data=program_scores, compression="gzip", chunks=(1, min(len(cnmf_cells), 4096)))

        gene_grp = h5.create_group("genes")
        write_string_dataset(gene_grp, "names", gene_names)
        gene_grp.create_dataset("indptr", data=X_genes_by_cells.indptr.astype(np.int64), compression="gzip")
        gene_grp.create_dataset("indices", data=X_genes_by_cells.indices.astype(np.int32), compression="gzip")
        gene_grp.create_dataset("data", data=X_genes_by_cells.data.astype(np.float32), compression="gzip")
        gene_grp.create_dataset("shape", data=np.asarray(X_genes_by_cells.shape, dtype=np.int64))

    print("Done.")
    print(f"Wrote: {out_path}")
    print("Use this file path as DATA_H5 in app.R, or set environment variable SHINY_UMAP_H5.")


if __name__ == "__main__":
    main()
