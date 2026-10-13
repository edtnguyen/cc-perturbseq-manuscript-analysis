#!/usr/bin/env python
"""
Export a compact HDF5 bundle for a Shiny UMAP viewer.

The Shiny app can display:
  1. cNMF program score on UMAP
  2. selected gene expression on UMAP
  3. selected cell metadata such as adata_int.obs['day'] on UMAP

UMAP and gene expression are taken from the integrated AnnData object.
cNMF scores are taken from the cNMF modality of the h5mu file.
Cells are aligned by obs_names before export.

This version prints timestamped progress messages and opens the integrated
.h5ad in backed='r' mode by default, so marker-gene tests do not need to load
all genes into memory.
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path
from typing import Any


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


log("Starting export script; importing Python packages...")

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
from scipy import sparse
import matplotlib as mpl
import matplotlib.pyplot as plt

log("Finished importing Python packages.")


DEFAULT_MDATA_PATH = (
    "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/"
    "030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/adata/cNMF_50_0_2.h5mu"
)

DEFAULT_ADATA_INT_PATH = (
    "/oak/stanford/groups/engreitz/Users/opushkar/cc_perturb/"
    "2.preprocessing_tonys_pipeline/3.integrated_days/"
    "3.all_days_concat_pc30_all_genes_sceptre_guides.h5ad"
)


def elapsed(t0: float) -> str:
    return f"{time.time() - t0:,.1f}s"


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
        log(f"Using `{color_column}` from program metadata for program colors")
    else:
        log(
            "No program color column found. Generating fallback colors. "
            f"Available columns: {list(prog_meta.columns)}"
        )

    fallback_cmap = plt.get_cmap("tab20")
    for idx, program in enumerate(program_ids.astype(int)):
        if int(program) not in prog_color:
            prog_color[int(program)] = mpl.colors.to_hex(fallback_cmap(idx % fallback_cmap.N))

    return label_map, prog_color



def infer_all_program_display_ids(cnmf: Any, n_programs: int, label_map: dict[int, str]) -> tuple[list[int], str]:
    """
    Infer display IDs for all cNMF program score columns.

    The score matrix columns are always selected by 0-based column position.
    Display IDs are only used for names/labels/colors shown in Shiny. Metadata
    tables can be either 0-based (0..N-1) or 1-based (1..N), so this function
    tries to infer the most readable display IDs without changing score lookup.
    """
    var_names = getattr(cnmf, "var_names", None)
    if var_names is not None and len(var_names) == n_programs:
        parsed: list[int] = []
        ok = True
        for raw in np.asarray(var_names.astype(str)):
            match = re.search(r"(\d+)$", str(raw))
            if match is None:
                ok = False
                break
            parsed.append(int(match.group(1)))
        if ok and len(set(parsed)) == n_programs:
            return parsed, "cnmf.var_names"

    meta_ids = set(int(x) for x in label_map.keys())
    zero_based_ids = set(range(n_programs))
    one_based_ids = set(range(1, n_programs + 1))

    if zero_based_ids.issubset(meta_ids):
        return list(range(n_programs)), "0-based program metadata"
    if one_based_ids.issubset(meta_ids):
        return list(range(1, n_programs + 1)), "1-based program metadata"

    return list(range(n_programs)), "0-based score-column fallback"


def lookup_program_metadata(mapping: dict[int, str], display_id: int, score_col: int, default: str) -> str:
    """Lookup label/color using display ID first, then common score-column fallbacks."""
    for key in (display_id, score_col, score_col + 1):
        if key in mapping:
            value = mapping[key]
            if value is not None and str(value).strip():
                return str(value)
    return default


def get_expr_matrix_and_genes(adata: Any, source: str):
    """
    Return a cells x genes matrix-like object and gene names.

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


def reorder_rows_to_cnmf(adata_cells: pd.Index, cnmf_cells: pd.Index) -> slice | np.ndarray:
    """Return a row indexer that puts adata rows in cnmf cell order."""
    missing_from_adata = cnmf_cells.difference(adata_cells)
    if len(missing_from_adata) > 0:
        raise ValueError(
            f"adata_int is missing {len(missing_from_adata)} cNMF cells. "
            f"First few: {missing_from_adata[:5].tolist()}"
        )

    if adata_cells.equals(cnmf_cells):
        log("Cell order already matches between cNMF and adata_int.")
        return slice(None)

    log("Computing adata_int row order to match cnmf.obs_names...")
    row_indexer = adata_cells.get_indexer(cnmf_cells)
    if np.any(row_indexer < 0):
        raise ValueError("Internal error: some cNMF cells were not found in adata_int.")
    return row_indexer.astype(np.int64)


def subset_matrix_preserve_order(X: Any, row_indexer: slice | np.ndarray, col_indexer: slice | np.ndarray):
    """
    Subset a cells x genes matrix-like object.

    HDF5-backed arrays often require fancy indices to be increasing. For explicit
    column/row index arrays, read in sorted order then restore requested order.
    """
    # Column subsetting first. For marker-gene export this usually materializes a
    # small cells x marker_genes matrix instead of the full all-gene matrix.
    if isinstance(col_indexer, slice):
        X_sel = X
    else:
        col_indexer = np.asarray(col_indexer, dtype=np.int64)
        col_order = np.argsort(col_indexer)
        col_sorted = col_indexer[col_order]
        log(f"Reading selected expression columns: {len(col_sorted):,} genes")
        X_sel = X[:, col_sorted]
        col_restore = np.argsort(col_order)
        X_sel = X_sel[:, col_restore]

    if not isinstance(row_indexer, slice):
        row_indexer = np.asarray(row_indexer, dtype=np.int64)
        # If X_sel is still backed/HDF5-like, increasing indices are safer.
        row_order = np.argsort(row_indexer)
        row_sorted = row_indexer[row_order]
        if len(row_sorted) and not np.all(row_sorted[:-1] <= row_sorted[1:]):
            raise RuntimeError("Internal sorting failure for row indexer.")
        log(f"Reordering/subsetting expression rows: {len(row_sorted):,} cells")
        X_sel = X_sel[row_sorted, :]
        row_restore = np.argsort(row_order)
        X_sel = X_sel[row_restore, :]

    return X_sel


def materialize_to_csr(X: Any):
    """Convert sparse/dense/backed matrix-like object to scipy CSR."""
    if sparse.issparse(X):
        return X.tocsr()
    if hasattr(X, "to_memory"):
        X_mem = X.to_memory()
        if sparse.issparse(X_mem):
            return X_mem.tocsr()
        return sparse.csr_matrix(np.asarray(X_mem))
    return sparse.csr_matrix(np.asarray(X))


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
    parser.add_argument(
        "--program-source",
        choices=["all", "program_order"],
        default="all",
        help=(
            "Which cNMF programs to export. Default `all` writes every score column "
            "from cnmf.X. Use `program_order` only to intentionally restrict to "
            "--program-order."
        ),
    )
    parser.add_argument(
        "--program-order",
        type=Path,
        default=None,
        help="Optional. Used only when --program-source program_order is selected.",
    )
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
        help="Comma-separated adata_int.obs columns to export if present.",
    )
    parser.add_argument(
        "--required-obs-columns",
        default="day",
        help="Comma-separated adata_int.obs columns that must be present. Use an empty string to disable this check.",
    )
    parser.add_argument(
        "--allow-extra-adata-cells",
        action="store_true",
        help="Allow adata_int to contain extra cells; it will be subset/reordered to cNMF cells.",
    )
    parser.add_argument(
        "--load-adata-into-memory",
        action="store_true",
        help="Use sc.read_h5ad(...), not backed='r'. This is slower/heavier but can help if backed slicing fails.",
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve() if args.repo_root is not None else find_repo_root()
    program_meta_path = args.program_meta or repo_root / "data" / "program_activity_long_v4.csv"
    program_order_path = args.program_order or repo_root / "analysis" / "regulator_dotplots" / "program_order.txt"
    out_path = args.out or repo_root / "analysis" / "regulator_dotplots" / "shiny_umap_app" / "umap_gene_program_data.h5"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    log(f"Repo root:       {repo_root}")
    log(f"cNMF h5mu:       {args.mdata}")
    log(f"Integrated h5ad: {args.adata_int}")
    log(f"Program meta:    {program_meta_path}")
    log(f"Program source:  {args.program_source}")
    if args.program_source == "program_order":
        log(f"Program order:   {program_order_path}")
    else:
        log("Program order:   ignored; exporting every cNMF score column")
    log(f"Output HDF5:     {out_path}")

    t0 = time.time()
    log("Loading cNMF h5mu...")
    mdata = mu.read_h5mu(args.mdata)
    log(f"Loaded cNMF h5mu in {elapsed(t0)}")

    if "cNMF" not in mdata.mod:
        raise KeyError(f"Expected modality `cNMF`. Available modalities: {list(mdata.mod.keys())}")
    cnmf = mdata.mod["cNMF"]
    cnmf_cells = pd.Index(cnmf.obs_names.astype(str))
    if not cnmf_cells.is_unique:
        raise ValueError("cnmf.obs_names are not unique.")

    t0 = time.time()
    if args.load_adata_into_memory:
        log("Loading integrated AnnData into memory...")
        adata_int = sc.read_h5ad(args.adata_int)
    else:
        log("Opening integrated AnnData with backed='r'...")
        adata_int = sc.read_h5ad(args.adata_int, backed="r")
    log(f"Opened integrated AnnData: {adata_int.n_obs:,} cells x {adata_int.n_vars:,} genes in {elapsed(t0)}")

    adata_cells = pd.Index(adata_int.obs_names.astype(str))
    if not adata_cells.is_unique:
        raise ValueError("adata_int.obs_names are not unique.")

    extra_in_adata = adata_cells.difference(cnmf_cells)
    if len(extra_in_adata) > 0 and not args.allow_extra_adata_cells:
        raise ValueError(
            f"adata_int has {len(extra_in_adata)} extra cells not present in cNMF. "
            "If this is expected, rerun with --allow-extra-adata-cells. "
            f"First few extras: {extra_in_adata[:5].tolist()}"
        )

    row_indexer = reorder_rows_to_cnmf(adata_cells, cnmf_cells)
    obs_aligned = adata_int.obs.iloc[row_indexer].copy()

    if "X_umap" not in adata_int.obsm:
        raise KeyError(f"adata_int.obsm does not contain X_umap. Available keys: {list(adata_int.obsm.keys())}")
    log("Reading and aligning X_umap...")
    umap_all = np.asarray(adata_int.obsm["X_umap"])
    umap = umap_all[row_indexer, :]
    if umap.shape[0] != cnmf.n_obs or umap.shape[1] < 2:
        raise ValueError(f"Unexpected X_umap shape after alignment: {umap.shape}")
    umap = umap[:, :2]

    log("Reading program metadata...")
    prog_meta = pd.read_csv(program_meta_path)
    label_map, color_map = make_label_and_color_maps(prog_meta)

    t0 = time.time()
    log("Preparing cNMF program score matrix...")
    score_matrix = cnmf.X
    if sparse.issparse(score_matrix):
        score_matrix = score_matrix.toarray()
    else:
        score_matrix = np.asarray(score_matrix)
    if score_matrix.shape[0] != len(cnmf_cells):
        raise ValueError("cNMF score matrix row count does not match cnmf.obs_names.")

    n_programs_total = int(score_matrix.shape[1])
    if args.program_source == "program_order":
        # Legacy behavior: restrict/order by program_order.txt. The parsed program IDs
        # are assumed to be score-matrix column positions.
        program_col_positions = parse_program_order(program_order_path, deduplicate=True)
        invalid_programs = [p for p in program_col_positions if p < 0 or p >= n_programs_total]
        if invalid_programs:
            raise ValueError(
                f"These program IDs are outside the cNMF score matrix column range: {invalid_programs}. "
                f"Score matrix has {n_programs_total} columns."
            )
        display_ids = list(program_col_positions)
        log(f"Exporting {len(program_col_positions):,} programs from program_order.txt.")
    else:
        # Default behavior: export every cNMF score column. This intentionally
        # ignores program_order.txt.
        program_col_positions = list(range(n_programs_total))
        display_ids, display_id_source = infer_all_program_display_ids(
            cnmf=cnmf,
            n_programs=n_programs_total,
            label_map=label_map,
        )
        log(
            f"Exporting all {n_programs_total:,} cNMF program score columns "
            f"with display IDs inferred from {display_id_source}."
        )

    program_scores = score_matrix[:, program_col_positions].T.astype(np.float32)
    program_names = [f"P{program_id}" for program_id in display_ids]
    program_labels = [
        lookup_program_metadata(label_map, display_id, score_col, f"P{display_id}")
        for display_id, score_col in zip(display_ids, program_col_positions)
    ]
    program_colors = [
        lookup_program_metadata(color_map, display_id, score_col, "#999999")
        for display_id, score_col in zip(display_ids, program_col_positions)
    ]
    log(f"Prepared {len(display_ids):,} programs in {elapsed(t0)}")

    t0 = time.time()
    log(f"Preparing gene expression matrix from source `{args.expression_source}`...")
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
                flush=True,
            )
        genes_to_export = [g for g in requested_genes if g in gene_index]
        if not genes_to_export:
            raise ValueError("No requested genes were found in adata_int.")
        col_indexer = gene_index.get_indexer(genes_to_export).astype(np.int64)
        gene_names = np.asarray(genes_to_export, dtype=str)
    else:
        col_indexer = slice(None)
        gene_names = gene_names_all.astype(str)

    log(f"Genes requested for export: {len(gene_names):,}")
    X = subset_matrix_preserve_order(X, row_indexer, col_indexer)
    X_cells_by_genes = materialize_to_csr(X)
    log(f"Expression matrix for export: {X_cells_by_genes.shape[0]:,} cells x {X_cells_by_genes.shape[1]:,} genes")
    log(f"Expression matrix nonzeros: {X_cells_by_genes.nnz:,}")
    log(f"Prepared gene matrix in {elapsed(t0)}")

    t0 = time.time()
    log("Transposing to gene-by-cell CSR for fast per-gene retrieval in Shiny...")
    X_genes_by_cells = X_cells_by_genes.T.tocsr()
    X_genes_by_cells.sort_indices()
    log(f"Transposed expression matrix in {elapsed(t0)}")

    obs_columns_requested = [c.strip() for c in args.obs_columns.split(",") if c.strip()]
    required_obs_columns = [c.strip() for c in args.required_obs_columns.split(",") if c.strip()]
    available_obs_columns = list(obs_aligned.columns.astype(str))

    missing_required_obs_cols = [c for c in required_obs_columns if c not in obs_aligned.columns]
    if missing_required_obs_cols:
        suggestions = {}
        for col in missing_required_obs_cols:
            suggestions[col] = [c for c in available_obs_columns if c.lower() == col.lower()]
        raise KeyError(
            "Required adata_int.obs columns are missing: "
            f"{missing_required_obs_cols}. Case-insensitive matches: {suggestions}. "
            "If the column has a different name, rerun with --obs-columns and "
            "--required-obs-columns set to the actual column name."
        )

    for col in required_obs_columns:
        if col not in obs_columns_requested:
            obs_columns_requested.append(col)

    obs_columns = [c for c in obs_columns_requested if c in obs_aligned.columns]
    missing_obs_cols = [c for c in obs_columns_requested if c not in obs_aligned.columns]
    if missing_obs_cols:
        log(f"Skipping missing optional obs columns: {missing_obs_cols}")

    log(f"Exporting obs annotations: {obs_columns}")
    if "day" in obs_aligned.columns:
        day_counts = obs_aligned["day"].astype(str).value_counts(dropna=False).to_dict()
        log(f"day counts after cell alignment: {day_counts}")

    t0 = time.time()
    log("Writing HDF5 bundle...")
    with h5py.File(out_path, "w") as h5:
        h5.attrs["description"] = "UMAP + cNMF program scores + sparse gene expression for Shiny viewer"
        h5.attrs["mdata_path"] = str(args.mdata)
        h5.attrs["adata_int_path"] = str(args.adata_int)
        h5.attrs["expression_source"] = args.expression_source
        h5.attrs["program_source"] = args.program_source

        obs_grp = h5.create_group("obs")
        write_string_dataset(obs_grp, "cell_id", cnmf_cells.tolist())
        obs_grp.create_dataset("umap_1", data=umap[:, 0].astype(np.float32), compression="gzip")
        obs_grp.create_dataset("umap_2", data=umap[:, 1].astype(np.float32), compression="gzip")

        ann_grp = h5.create_group("obs_annotations")
        write_string_dataset(h5, "obs_annotation_names", obs_columns)
        for col in obs_columns:
            write_string_dataset(ann_grp, col, obs_aligned[col].astype(str).fillna("NA").tolist())

        prog_grp = h5.create_group("programs")
        prog_grp.create_dataset("ids", data=np.asarray(display_ids, dtype=np.int32), compression="gzip")
        prog_grp.create_dataset("score_columns", data=np.asarray(program_col_positions, dtype=np.int32), compression="gzip")
        write_string_dataset(prog_grp, "names", program_names)
        write_string_dataset(prog_grp, "labels", program_labels)
        write_string_dataset(prog_grp, "colors", program_colors)
        prog_grp.create_dataset(
            "scores",
            data=program_scores,
            compression="gzip",
            chunks=(1, min(len(cnmf_cells), 4096)),
        )

        gene_grp = h5.create_group("genes")
        write_string_dataset(gene_grp, "names", gene_names)
        gene_grp.create_dataset("indptr", data=X_genes_by_cells.indptr.astype(np.int64), compression="gzip")
        gene_grp.create_dataset("indices", data=X_genes_by_cells.indices.astype(np.int32), compression="gzip")
        gene_grp.create_dataset("data", data=X_genes_by_cells.data.astype(np.float32), compression="gzip")
        gene_grp.create_dataset("shape", data=np.asarray(X_genes_by_cells.shape, dtype=np.int64))

    log(f"Wrote HDF5 bundle in {elapsed(t0)}")
    log(f"Done: {out_path}")
    log("Use this file path as DATA_H5 in app.R, or set environment variable SHINY_UMAP_H5.")

    if getattr(adata_int, "isbacked", False):
        try:
            adata_int.file.close()
        except Exception:
            pass


if __name__ == "__main__":
    main()
