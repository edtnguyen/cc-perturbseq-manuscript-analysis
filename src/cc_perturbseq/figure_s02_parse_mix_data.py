"""Source-faithful producers for the mixed hiPSC–K562 Parse v3 panels in Figure S2.

Scientific sources
------------------
* ``recovery/raw_downloads/250904.Parse.v3.iPSC.K562.copy.ipynb``
* Figure S2 caption and Parse-v3 compatibility methods in the manuscript.

This module deliberately implements only panels i, j, l, m and n. Panel k is
part of the earlier hiPSC-only circular-capture/CROP-seq experiment and is not
part of the mixed hiPSC–K562 Parse-v3 experiment.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


SOURCE_NOTEBOOK = "recovery/raw_downloads/250904.Parse.v3.iPSC.K562.copy.ipynb"
SOURCE_NOTEBOOK_GIT_BLOB = "ff940a4b35204cba38db1e453b9968d07372771d"
PRIMER_ORDER = ("CS1", "DC1", "Dual", "no_guide")
CAPTURE_PRIMER_ORDER = ("CS1", "DC1", "Dual")
MANUSCRIPT_TRANSCRIPT_ORDER = ("no_guide", "Dual", "DC1", "CS1")
CELL_TYPE_ORDER = ("K562", "hiPSC")
PRIMER_COLORS = {
    "CS1": "#647FBC",
    "DC1": "#91ADC8",
    "Dual": "#AED6CF",
    "no_guide": "#FAFDD6",
}
CELL_TYPE_COLORS = {"K562": "#4E79A7", "hiPSC": "#F28E2B"}
MARKERS = {
    "K562": ("GATA1", "KLF1", "HBB", "HBA1", "HBA2"),
    "hiPSC": ("SOX2", "NANOG", "POU5F1", "LIN28A"),
}
EXPECTED_MATCHED_POOL_MEDIANS = {
    ("K562", "CS1"): 217.0,
    ("K562", "DC1"): 238.0,
    ("K562", "Dual"): 301.0,
    ("hiPSC", "CS1"): 138.0,
    ("hiPSC", "DC1"): 162.0,
    ("hiPSC", "Dual"): 175.0,
}
EXPECTED_SPECIFICITY_MEDIANS = {
    ("K562", "CS1"): 0.86,
    ("K562", "DC1"): 0.96,
    ("K562", "Dual"): 0.92,
    ("hiPSC", "CS1"): 0.72,
    ("hiPSC", "DC1"): 0.93,
    ("hiPSC", "Dual"): 0.83,
}


@dataclass
class ParseMixResult:
    cell_metrics: pd.DataFrame
    cluster_annotation: pd.DataFrame
    saturation: pd.DataFrame
    median_diagnostics: pd.DataFrame
    provenance: dict[str, Any]


def _file_record(path: Path) -> dict[str, Any]:
    path = Path(path).expanduser().resolve()
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def _canonical_cell_ids(values: Iterable[Any]) -> pd.Index:
    index = pd.Index(values).astype(str)
    return pd.Index(index.str.replace(r"__s1$", "", regex=True), name="cell_id")


def _canonical_guide_names(values: pd.Series) -> pd.Series:
    return values.astype(str).str.replace(r"^gRNA_", "", regex=True)


def _guide_pool(values: pd.Series) -> pd.Series:
    values = values.astype(str)
    pool = np.select(
        [values.str.startswith("iPSC_"), values.str.startswith("K562_")],
        ["hiPSC", "K562"],
        default="other",
    )
    return pd.Series(pool, index=values.index, dtype="object")


def _required_paths(root: Path) -> dict[str, Path]:
    root = Path(root).expanduser().resolve()
    return {
        "root": root,
        "expression_dir": root / "splitpipe_out" / "WT" / "sub1" / "all-sample" / "DGE_filtered",
        "count_matrix": root / "splitpipe_out" / "WT" / "sub1" / "all-sample" / "DGE_filtered" / "count_matrix.mtx",
        "genes": root / "splitpipe_out" / "WT" / "sub1" / "all-sample" / "DGE_filtered" / "all_genes.csv",
        "cell_metadata": root / "splitpipe_out" / "WT" / "sub1" / "all-sample" / "DGE_filtered" / "cell_metadata.csv",
        "guide_counts": root / "splitpipe_out" / "CRISPR_combined" / "grnaxcell_count.tsv",
        "guide_reference": root / "splitpipe_out" / "CRISPR" / "sub1" / "all-sample" / "guide_RNAs_unfiltered" / "all_guides.csv",
    }


def preflight_parse_mix(root: Path) -> dict[str, Path]:
    paths = _required_paths(root)
    missing = [
        str(path)
        for key, path in paths.items()
        if key != "root" and (not path.is_file() or path.stat().st_size == 0)
    ]
    for primer in PRIMER_ORDER:
        saturation = (
            paths["root"]
            / "splitpipe_out"
            / "WT"
            / "sub1"
            / primer
            / "report"
            / "tscp_counts_subsampled.csv"
        )
        paths[f"saturation_{primer}"] = saturation
        if not saturation.is_file() or saturation.stat().st_size == 0:
            missing.append(str(saturation))
    if missing:
        raise FileNotFoundError(
            "Missing or empty Parse-v3 mixed hiPSC–K562 inputs:\n  "
            + "\n  ".join(missing)
        )
    return paths


def read_expression_and_qc(root: Path):
    """Read Split-pipe outputs and apply the source notebook's three QC cutoffs."""
    import anndata as ad
    import scanpy as sc

    paths = preflight_parse_mix(root)
    matrix = sparse.csr_matrix(mmread(paths["count_matrix"]), dtype=np.float32)
    genes = pd.read_csv(paths["genes"])
    cells = pd.read_csv(paths["cell_metadata"])
    if "gene_name" not in genes:
        raise KeyError(f"{paths['genes']} lacks 'gene_name'.")
    if "bc_wells" not in cells or "sample" not in cells:
        raise KeyError(f"{paths['cell_metadata']} requires 'bc_wells' and 'sample'.")

    n_cells, n_genes = len(cells), len(genes)
    if matrix.shape == (n_genes, n_cells):
        matrix = matrix.T.tocsr()
    if matrix.shape != (n_cells, n_genes):
        raise ValueError(
            f"Expression matrix shape {matrix.shape} is incompatible with "
            f"{n_cells} cells and {n_genes} genes."
        )

    keep_gene = genes["gene_name"].notna().to_numpy()
    matrix = matrix[:, keep_gene].tocsr()
    genes = genes.loc[keep_gene].copy()
    obs_names = _canonical_cell_ids(cells["bc_wells"])
    if not obs_names.is_unique:
        duplicates = obs_names[obs_names.duplicated()].unique().tolist()[:10]
        raise ValueError(f"Cell IDs are not unique after removing __s1: {duplicates}")
    cells = cells.copy()
    cells.index = obs_names
    var = genes.copy()
    var.index = pd.Index(var["gene_name"].astype(str), name="gene")
    adata = ad.AnnData(X=matrix, obs=cells, var=var)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    keep = (
        (adata.obs["n_genes_by_counts"] < 5000)
        & (adata.obs["total_counts"] < 5000)
        & (adata.obs["pct_counts_mt"] < 40)
    )
    adata = adata[keep].copy()
    if adata.n_obs == 0:
        raise ValueError("The source QC cutoffs retained no cells.")
    unknown_samples = sorted(set(adata.obs["sample"].astype(str)) - set(PRIMER_ORDER))
    if unknown_samples:
        raise ValueError(f"Unexpected Round-1 primer labels: {unknown_samples}")
    return adata, paths


def assign_transcriptomic_cell_types(adata) -> tuple[pd.Series, pd.DataFrame]:
    """Reproduce the source dimensionality reduction and manuscript marker scoring."""
    import scanpy as sc

    all_markers = [gene for group in CELL_TYPE_ORDER for gene in MARKERS[group]]
    missing = [gene for gene in all_markers if gene not in adata.var_names]
    if missing:
        raise KeyError(f"Expression matrix lacks lineage markers: {missing}")

    normalized = adata.copy()
    sc.pp.normalize_total(normalized, target_sum=1e4)
    sc.pp.log1p(normalized)
    marker_values = normalized[:, all_markers].X
    marker_values = (
        marker_values.toarray()
        if sparse.issparse(marker_values)
        else np.asarray(marker_values)
    )
    marker_frame = pd.DataFrame(
        marker_values, index=normalized.obs_names, columns=all_markers
    )
    per_cell_scores = pd.DataFrame(
        {
            cell_type: marker_frame.loc[:, list(MARKERS[cell_type])].mean(axis=1)
            for cell_type in CELL_TYPE_ORDER
        },
        index=normalized.obs_names,
    )

    sc.pp.highly_variable_genes(
        normalized, min_mean=0.0125, max_mean=3, min_disp=0.25
    )
    if not bool(normalized.var["highly_variable"].any()):
        raise ValueError("No highly variable genes passed the source thresholds.")
    work = normalized[:, normalized.var["highly_variable"]].copy()
    sc.pp.regress_out(work, ["total_counts"])
    sc.pp.scale(work, max_value=10)
    n_comps = min(50, work.n_obs - 1, work.n_vars - 1)
    if n_comps < 30:
        raise ValueError(
            f"Only {n_comps} PCs are available; source workflow requires at least 30."
        )
    sc.tl.pca(work, svd_solver="arpack", n_comps=n_comps, random_state=0)
    sc.pp.neighbors(work, n_neighbors=20, n_pcs=30, random_state=0)
    sc.tl.umap(work, random_state=0)
    sc.tl.leiden(work, resolution=0.2, random_state=0, key_added="leiden")

    clusters = work.obs["leiden"].astype(str)
    score_table = per_cell_scores.join(clusters.rename("leiden"))
    cluster_scores = (
        score_table.groupby("leiden", observed=True)[list(CELL_TYPE_ORDER)]
        .mean()
        .rename(
            columns={
                "K562": "K562_marker_score",
                "hiPSC": "hiPSC_marker_score",
            }
        )
    )
    cluster_scores["n_cells"] = (
        clusters.value_counts().reindex(cluster_scores.index).astype(int)
    )

    def classify(row: pd.Series) -> str:
        k562 = float(row["K562_marker_score"])
        hipsc = float(row["hiPSC_marker_score"])
        if max(k562, hipsc) <= 0 or np.isclose(k562, hipsc):
            return "unassigned"
        return "K562" if k562 > hipsc else "hiPSC"

    cluster_scores["cell_type"] = cluster_scores.apply(classify, axis=1)
    mapping = cluster_scores["cell_type"].to_dict()
    cell_types = clusters.map(mapping).astype("object")
    no_detected_marker = per_cell_scores.max(axis=1).le(0)
    cell_types.loc[no_detected_marker] = "unassigned"
    cell_types.name = "cell_type"
    return cell_types, cluster_scores.reset_index()


def read_guide_pairs(
    guide_counts_path: Path,
    guide_reference_path: Path,
    valid_cells: pd.Index,
    *,
    chunksize: int = 1_000_000,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Read the long UMI table without constructing the source's dense matrix."""
    reference = pd.read_csv(guide_reference_path)
    if "gene_name" not in reference:
        raise KeyError(f"{guide_reference_path} lacks 'gene_name'.")
    guide_names = pd.Index(reference["gene_name"].dropna().astype(str))
    if not guide_names.is_unique:
        duplicated = guide_names[guide_names.duplicated()].unique().tolist()[:10]
        raise ValueError(f"Filtered guide reference is not unique: {duplicated}")
    valid_guides = set(guide_names)
    valid_cell_set = set(_canonical_cell_ids(valid_cells))

    pieces: list[pd.DataFrame] = []
    n_rows_read = 0
    n_rows_retained = 0
    for chunk in pd.read_csv(
        guide_counts_path,
        sep="\t",
        usecols=["gene", "cell", "count"],
        chunksize=chunksize,
    ):
        n_rows_read += len(chunk)
        chunk["cell"] = _canonical_cell_ids(chunk["cell"])
        chunk["guide"] = _canonical_guide_names(chunk["gene"])
        chunk["count"] = pd.to_numeric(chunk["count"], errors="coerce")
        keep = (
            chunk["cell"].isin(valid_cell_set)
            & chunk["guide"].isin(valid_guides)
            & chunk["count"].notna()
            & chunk["count"].gt(0)
        )
        if bool(keep.any()):
            retained = chunk.loc[keep, ["cell", "guide", "count"]].copy()
            pieces.append(retained)
            n_rows_retained += len(retained)
    if not pieces:
        raise ValueError(
            "No guide-count observations matched the QC-passing expression cells."
        )
    pairs = pd.concat(pieces, ignore_index=True)
    pairs = (
        pairs.groupby(["cell", "guide"], observed=True, as_index=False)["count"]
        .sum()
    )
    pairs["pool"] = _guide_pool(pairs["guide"])
    diagnostics = {
        "n_rows_read": int(n_rows_read),
        "n_rows_retained_before_pair_collapse": int(n_rows_retained),
        "n_unique_cell_guide_pairs": int(len(pairs)),
        "n_reference_guides": int(len(guide_names)),
        "n_iPSC_reference_guides": int(
            guide_names.str.startswith("iPSC_").sum()
        ),
        "n_K562_reference_guides": int(
            guide_names.str.startswith("K562_").sum()
        ),
        "n_other_pair_rows": int(pairs["pool"].eq("other").sum()),
    }
    return pairs, diagnostics


def summarize_guide_pairs(
    pairs: pd.DataFrame,
    cell_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Compute the per-cell quantities used by S2i, j, l and n."""
    required_pairs = {"cell", "guide", "count"}
    missing_pairs = required_pairs - set(pairs.columns)
    if missing_pairs:
        raise KeyError(f"Guide-pair table lacks columns: {sorted(missing_pairs)}")
    required_meta = {"sample", "cell_type"}
    missing_meta = required_meta - set(cell_metadata.columns)
    if missing_meta:
        raise KeyError(f"Cell metadata lacks columns: {sorted(missing_meta)}")

    pairs = pairs.copy()
    pairs["cell"] = _canonical_cell_ids(pairs["cell"])
    pairs["guide"] = _canonical_guide_names(pairs["guide"])
    pairs["count"] = pd.to_numeric(pairs["count"], errors="raise")
    if "pool" not in pairs:
        pairs["pool"] = _guide_pool(pairs["guide"])
    pairs = (
        pairs.groupby(
            ["cell", "guide", "pool"], observed=True, as_index=False
        )["count"]
        .sum()
    )
    cells = cell_metadata.copy()
    cells.index = _canonical_cell_ids(cells.index)
    if not cells.index.is_unique:
        raise ValueError("Cell metadata IDs are not unique after canonicalization.")

    metrics = cells.copy()
    metrics["total_guide_umis"] = (
        pairs.groupby("cell", observed=True)["count"]
        .sum()
        .reindex(metrics.index, fill_value=0.0)
    )
    metrics["top_guide_umis"] = (
        pairs.groupby("cell", observed=True)["count"]
        .max()
        .reindex(metrics.index, fill_value=0.0)
    )
    pool_sums = (
        pairs.loc[pairs["pool"].isin(CELL_TYPE_ORDER)]
        .pivot_table(
            index="cell",
            columns="pool",
            values="count",
            aggfunc="sum",
            fill_value=0.0,
        )
        .reindex(
            index=metrics.index,
            columns=list(CELL_TYPE_ORDER),
            fill_value=0.0,
        )
    )
    metrics["K562_guide_umis"] = pool_sums["K562"].to_numpy(float)
    metrics["hiPSC_guide_umis"] = pool_sums["hiPSC"].to_numpy(float)

    thresholded = pairs.loc[
        pairs["count"].ge(3) & pairs["pool"].isin(CELL_TYPE_ORDER)
    ]
    thresholded_sums = (
        thresholded.pivot_table(
            index="cell",
            columns="pool",
            values="count",
            aggfunc="sum",
            fill_value=0.0,
        )
        .reindex(
            index=metrics.index,
            columns=list(CELL_TYPE_ORDER),
            fill_value=0.0,
        )
    )
    metrics["K562_guide_umis_ge3"] = thresholded_sums["K562"].to_numpy(float)
    metrics["hiPSC_guide_umis_ge3"] = thresholded_sums["hiPSC"].to_numpy(float)

    metrics["matched_pool_guide_umis"] = np.select(
        [metrics["cell_type"].eq("K562"), metrics["cell_type"].eq("hiPSC")],
        [metrics["K562_guide_umis"], metrics["hiPSC_guide_umis"]],
        default=np.nan,
    )
    pool_total = metrics["K562_guide_umis"] + metrics["hiPSC_guide_umis"]
    metrics["guide_pool_specificity"] = np.divide(
        metrics["matched_pool_guide_umis"],
        pool_total,
        out=np.full(len(metrics), np.nan, dtype=float),
        where=pool_total.to_numpy(float) > 0,
    )
    metrics.index.name = "cell_id"
    return metrics


def read_saturation(root: Path) -> pd.DataFrame:
    root = Path(root).expanduser().resolve()
    frames = []
    for primer in PRIMER_ORDER:
        path = (
            root
            / "splitpipe_out"
            / "WT"
            / "sub1"
            / primer
            / "report"
            / "tscp_counts_subsampled.csv"
        )
        raw = pd.read_csv(path, header=None)
        if len(raw) < 2 or raw.shape[1] < 2:
            raise ValueError(f"Malformed saturation table: {path}")
        table = raw.iloc[1:, :2].copy()
        table.columns = [
            "mean_assigned_reads_per_cell",
            "median_transcript_umis_per_cell",
        ]
        for column in table:
            table[column] = pd.to_numeric(table[column], errors="coerce")
        table = table.dropna().sort_values("mean_assigned_reads_per_cell")
        if table.empty:
            raise ValueError(f"No numeric saturation points in {path}")
        table["primer"] = primer
        table["source_file"] = str(path)
        frames.append(table)
    result = pd.concat(frames, ignore_index=True)
    result["primer"] = pd.Categorical(
        result["primer"], PRIMER_ORDER, ordered=True
    )
    return result


def manuscript_median_diagnostics(cell_metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric, expected in (
        ("matched_pool_guide_umis", EXPECTED_MATCHED_POOL_MEDIANS),
        ("guide_pool_specificity", EXPECTED_SPECIFICITY_MEDIANS),
    ):
        for (cell_type, primer), recorded in expected.items():
            values = cell_metrics.loc[
                cell_metrics["cell_type"].eq(cell_type)
                & cell_metrics["sample"].astype(str).eq(primer),
                metric,
            ].dropna()
            observed = float(values.median()) if len(values) else np.nan
            rows.append(
                {
                    "metric": metric,
                    "cell_type": cell_type,
                    "primer": primer,
                    "n_cells": int(len(values)),
                    "observed_median": observed,
                    "manuscript_median": float(recorded),
                    "difference": observed - float(recorded),
                    "acceptance_role": "diagnostic_only",
                }
            )
    return pd.DataFrame(rows)


def prepare_parse_mix_benchmark(
    root: Path,
    *,
    guide_chunksize: int = 1_000_000,
) -> ParseMixResult:
    adata, paths = read_expression_and_qc(root)
    cell_types, cluster_annotation = assign_transcriptomic_cell_types(adata)
    cell_metadata = adata.obs[
        ["sample", "total_counts", "n_genes_by_counts", "pct_counts_mt"]
    ].copy()
    cell_metadata["cell_type"] = cell_types.reindex(cell_metadata.index)
    pairs, guide_diagnostics = read_guide_pairs(
        paths["guide_counts"],
        paths["guide_reference"],
        adata.obs_names,
        chunksize=guide_chunksize,
    )
    cell_metrics = summarize_guide_pairs(pairs, cell_metadata)
    saturation = read_saturation(paths["root"])
    diagnostics = manuscript_median_diagnostics(cell_metrics)
    provenance = {
        "scientific_source": SOURCE_NOTEBOOK,
        "source_notebook_git_blob": SOURCE_NOTEBOOK_GIT_BLOB,
        "experiment_root": str(paths["root"]),
        "input_files": {
            key: _file_record(path)
            for key, path in paths.items()
            if key != "root" and path.is_file()
        },
        "source_qc": {
            "n_genes_by_counts": "< 5000",
            "total_counts": "< 5000",
            "pct_counts_mt": "< 40",
        },
        "n_cells_after_qc": int(adata.n_obs),
        "n_genes_after_gene_name_filter": int(adata.n_vars),
        "cell_type_counts": {
            str(key): int(value)
            for key, value in cell_metrics["cell_type"]
            .value_counts(dropna=False)
            .items()
        },
        "guide_processing": {
            **guide_diagnostics,
            "cell_guide_minimum_for_panel_n": 3,
            "guide_prefixes": {"hiPSC": "iPSC_", "K562": "K562_"},
        },
        "historical_medians_are_diagnostic_only": True,
    }
    return ParseMixResult(
        cell_metrics=cell_metrics,
        cluster_annotation=cluster_annotation,
        saturation=saturation,
        median_diagnostics=diagnostics,
        provenance=provenance,
    )


__all__ = [
    "CAPTURE_PRIMER_ORDER",
    "CELL_TYPE_COLORS",
    "CELL_TYPE_ORDER",
    "EXPECTED_MATCHED_POOL_MEDIANS",
    "EXPECTED_SPECIFICITY_MEDIANS",
    "MANUSCRIPT_TRANSCRIPT_ORDER",
    "MARKERS",
    "PRIMER_COLORS",
    "PRIMER_ORDER",
    "ParseMixResult",
    "SOURCE_NOTEBOOK",
    "SOURCE_NOTEBOOK_GIT_BLOB",
    "manuscript_median_diagnostics",
    "prepare_parse_mix_benchmark",
    "preflight_parse_mix",
    "read_guide_pairs",
    "read_saturation",
    "summarize_guide_pairs",
]
