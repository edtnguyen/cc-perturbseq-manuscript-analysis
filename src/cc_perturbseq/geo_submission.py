"""Portable GEO processed objects for Hippo Figure 3f/4c and mixed-cell S2.

No fitting or alternative plotting implementation is introduced. Public reads
use only the H5MU and the existing figure producers in this repository.
"""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import gc
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import tempfile
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

SCHEMA = "ccperturbseq.geo.v1"
HIPPO_FILENAME = "Hippo_timecourse_integrated.h5mu"
MIX_FILENAME = "iPSC_K562_mixed_integrated.h5mu"
BASE_FIGURE_COMMIT = "a69cecdd1ed12625fff5b5fef7b4561d3bc78e18"  # source baseline; implementation hashes below identify the port
METRIC_COLUMNS = (
    "total_guide_umis", "top_guide_umis", "K562_guide_umis",
    "hiPSC_guide_umis", "K562_guide_umis_ge3", "hiPSC_guide_umis_ge3",
    "matched_pool_guide_umis", "guide_pool_specificity",
)


def _json(value: Any) -> str:
    def encode(x):
        if isinstance(x, (Path, datetime)):
            return str(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, np.generic):
            return x.item()
        raise TypeError(f"Not JSON serializable: {type(x).__name__}")
    return json.dumps(value, indent=2, sort_keys=True, default=encode)


def file_checksums(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty file: {path}")
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b""):
            sha.update(block)
            md5.update(block)
    return {"filename": path.name, "size_bytes": path.stat().st_size,
            "sha256": sha.hexdigest(), "md5": md5.hexdigest()}


def _versions() -> dict[str, str]:
    result = {}
    for name in ("anndata", "mudata", "scanpy", "numpy", "pandas", "scipy",
                 "h5py", "matplotlib", "seaborn", "plotnine", "igraph", "leidenalg"):
        try:
            result[name] = version(name)
        except PackageNotFoundError:
            result[name] = "not installed"
    return result


def _code_record() -> dict[str, Any]:
    here = Path(__file__).resolve().parent
    names = ("geo_submission.py", "hippo_projection.py", "figure3_panels_legacy.py",
             "figure4_panels_legacy.py", "figure_s02_parse_mix_data.py",
             "figure_s02_parse_mix_plots.py", "figure_s02_source_plotnine.py", "figure_s02_parse_mix_preflight.py")
    return {name: file_checksums(here / name)["sha256"] for name in names}


def _matrix_hash(matrix) -> str:
    """Hash values, shape and sparse positions without densifying count matrices."""
    digest = hashlib.sha256()
    digest.update(str(tuple(matrix.shape)).encode())
    if sparse.issparse(matrix):
        x = matrix.tocsr(copy=False)
        arrays = (x.indptr, x.indices, x.data)
        digest.update(b"csr")
    else:
        arrays = (np.asarray(matrix),)
        digest.update(b"dense")
    for a in arrays:
        digest.update(str(a.dtype).encode())
        # Read bounded row/entry blocks; do not duplicate the full sparse data.
        for start in range(0, len(a), 100_000):
            digest.update(np.ascontiguousarray(a[start:start + 100_000]).tobytes())
    return digest.hexdigest()


def _count_check(matrix, label: str) -> dict[str, Any]:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).reshape(-1)
    integer_like = True
    for start in range(0, len(values), 1_000_000):
        block = values[start:start + 1_000_000]
        if not np.isfinite(block).all() or (block < 0).any():
            raise ValueError(f"{label}: nonfinite or negative counts")
        integer_like = integer_like and bool(np.allclose(block, np.rint(block), atol=1e-6, rtol=0))
    return {"shape": list(matrix.shape), "nonnegative_finite": True,
            "integer_like": integer_like, "sparse": bool(sparse.issparse(matrix))}


def _safe_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Make analysis tables portable across AnnData/pandas string encodings."""
    out = frame.copy()
    out.index = pd.Index(out.index.astype(str), name=out.index.name)
    out.columns = out.columns.astype(str)
    for column in out:
        if isinstance(out[column].dtype, pd.CategoricalDtype):
            out[column] = out[column].astype(object)
        if not pd.api.types.is_numeric_dtype(out[column]):
            out[column] = out[column].fillna("").astype(str).astype(object)
    return out


def _same_table(observed: pd.DataFrame, expected: pd.DataFrame, label: str):
    left, right = _safe_frame(observed), _safe_frame(expected)
    if list(left.index) != list(right.index) or list(left.columns) != list(right.columns):
        raise AssertionError(f"{label}: row/column identity or order changed")
    try:
        pd.testing.assert_frame_equal(left, right, check_dtype=False,
                                      check_categorical=False, check_exact=False,
                                      rtol=1e-6, atol=1e-7)
    except AssertionError as exc:
        raise AssertionError(f"{label}: values changed after serialization/reconstruction") from exc


def _sample_table(obs: pd.DataFrame, candidates) -> pd.DataFrame:
    columns = [c for c in candidates if c in obs]
    if not columns:
        raise ValueError("No sample-identity columns available")
    return _safe_frame(obs.groupby(columns, observed=True, dropna=False).size()
                       .rename("n_cells").reset_index())


def _assemble(modalities: dict, dataset: str, provenance: dict, references: dict,
              sample_table: pd.DataFrame):
    import mudata as md
    ids = next(iter(modalities.values())).obs_names
    for name, obj in modalities.items():
        if not obj.obs_names.is_unique or not obj.var_names.is_unique:
            raise ValueError(f"{name}: cell/feature IDs must be unique")
        if not ids.equals(obj.obs_names):
            raise ValueError(f"{name}: cell order differs between modalities")
    # Individual modality obs tables are authoritative; shared obs is convenient.
    bundle = md.MuData(modalities)
    bundle.obs = next(iter(modalities.values())).obs.copy()
    bundle.uns["geo_schema"] = SCHEMA
    bundle.uns["dataset"] = dataset
    bundle.uns["provenance_json"] = _json({
        **provenance, "created_utc": datetime.now(timezone.utc).isoformat(),
        "versions": _versions(), "figure_code_baseline_commit": BASE_FIGURE_COMMIT,
        "implementation_sha256": _code_record(),
        "source_paths_are_provenance_only": True,
        "geo_sample_and_raw_read_metadata_still_required": True,
    })
    bundle.uns["sample_summary"] = _safe_frame(sample_table)
    bundle.uns["reference_tables"] = {k: _safe_frame(v) for k, v in references.items()}
    hashes = {}
    for name, obj in modalities.items():
        hashes[f"{name}/X"] = _matrix_hash(obj.X)
        for slot in ("layers", "obsm", "varm"):
            for key, value in getattr(obj, slot).items():
                if not isinstance(value, pd.DataFrame):
                    hashes[f"{name}/{slot}/{key}"] = _matrix_hash(value)
    bundle.uns["matrix_sha256_json"] = _json(hashes)
    return bundle



def _apply_source_style():
    """Same rcParams as figure_notebook.initialize_figure, without path loading."""
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.bbox": "tight",
    })


def _hippo_producers():
    from .figure3_panels import prepare_figure3f_dotplots
    from .figure4_panels import prepare_figure4c_dotplots
    return {"Figure_3f": prepare_figure3f_dotplots,
            "Figure_4c": prepare_figure4c_dotplots}


def _hippo_tables_and_plots(query_path, projection_path, render_dir=None):
    """Call the same producers used in the accepted Figure 3/4 notebooks."""
    import matplotlib.pyplot as plt
    _apply_source_style()
    tables = {}
    panel_metadata = {}
    for panel, prepare in _hippo_producers().items():
        gene_plot, program_plot, gene_table, program_table, metadata = prepare(
            query_path, projection_path)
        tables[f"{panel}_genes"] = gene_table
        tables[f"{panel}_programs"] = program_table
        # Also freeze actual Scanpy dotplot color/size matrices (not just summaries).
        for label, plot in (("genes", gene_plot), ("programs", program_plot)):
            tables[f"{panel}_{label}_dot_color"] = plot.dot_color_df
            tables[f"{panel}_{label}_dot_size"] = plot.dot_size_df
            if render_dir is not None:
                for suffix in ("pdf", "svg", "png"):
                    plot.savefig(Path(render_dir) / f"{panel}_{label}.{suffix}", dpi=220)
        panel_metadata[panel] = metadata
        plt.close("all")
    return tables, panel_metadata


def build_hippo_bundle(query_path, projection_path, spectra_path):
    """Integrate raw query and frozen projection; do not rerun or refit cNMF."""
    import anndata as ad
    from .hippo_projection import read_spectra, validate_projected
    from .figure3_panels import (_hippo_aligned_metadata, _hippo_panel_mask,
                                FIGURE3F_D3_CONDITIONS)
    from .figure4_panels import FIGURE4C_D12_CONDITIONS

    paths = [Path(p).expanduser().resolve() for p in (query_path, projection_path, spectra_path)]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    # Use current invariant-based validation; do not require historical cell counts.
    validate_projected(paths[0], paths[1], strict_historical=False)
    reference_tables, panel_metadata = _hippo_tables_and_plots(paths[0], paths[1])
    query, projection = ad.read_h5ad(paths[0]), ad.read_h5ad(paths[1])
    obs = _hippo_aligned_metadata(query, projection)
    if not query.obs_names.equals(projection.obs_names):
        query = query[projection.obs_names].copy()
    obs["used_in_Figure_3f"] = _hippo_panel_mask(obs, FIGURE3F_D3_CONDITIONS)
    obs["used_in_Figure_4c"] = _hippo_panel_mask(obs, FIGURE4C_D12_CONDITIONS)
    query.obs = obs.copy()
    projection.obs = obs.copy()
    counts = _count_check(query.X, "Hippo RNA X")
    ref = read_spectra(paths[2], strict_historical=False)
    if list(ref.index) != list(projection.var_names):
        raise ValueError("Reference spectra program order differs from projection")
    spectra_digest = file_checksums(paths[2])
    stored_digest = projection.uns.get("reference_spectra_sha256")
    if stored_digest and stored_digest != spectra_digest["sha256"]:
        raise ValueError("The selected reference spectra are not the basis recorded by the projection")
    genes = list(projection.uns["genes_used_for_projection"])
    if any(g not in ref.columns for g in genes):
        raise ValueError("Projection gene list is not contained in the supplied reference")
    projection.varm["reference_gene_spectra_tpm"] = ref.to_numpy()
    projection.uns["reference_gene_names"] = ref.columns.to_numpy(dtype=str)
    H = ref.loc[:, genes].to_numpy(dtype=np.float32).copy()
    H = H / np.maximum(H.sum(axis=1, keepdims=True), np.float32(1e-12))
    projection.varm["fixed_H_selected_normalized"] = H
    projection.uns["fixed_H_gene_names"] = np.asarray(genes, dtype=str)
    query.uns["geo_X_semantics"] = "Original unnormalized RNA counts used by Figure 3f/4c"
    bundle = _assemble(
        {"rna": query, "cnmf": projection}, "hippo_timecourse",
        {"input_files": [file_checksums(p) for p in paths], "RNA_counts": counts,
         "panel_metadata": panel_metadata,
         "scope": ["Figure_3f", "Figure_4c"],
         "other_Figure_3_4_panels_use_separate_experiments": True,
         "normalization": "Original producers normalize each panel subset separately; no substituted global layer",
         "cnmf_refitted": False}, reference_tables,
        _sample_table(obs, ("sample_name", "condition", "replicate", "treatment", "sublibID", "biological_sample")))
    return bundle


def _mix_producers():
    from . import figure_s02_parse_mix as s2
    return s2


def _mix_tables_and_plots(metrics, saturation, render_dir=None):
    import matplotlib.pyplot as plt
    s2 = _mix_producers()
    _apply_source_style()
    tables = {}
    for panel in ("i", "j", "l", "m", "n"):
        args = (metrics, saturation) if panel == "m" else (metrics,)
        result = getattr(s2, f"plot_panel_{panel}")(*args)
        fig, *data = result
        for index, frame in enumerate(data):
            tables[f"Figure_S02{panel}_source_{index}"] = frame
        if render_dir is not None:
            for suffix in ("pdf", "svg", "png"):
                fig.savefig(Path(render_dir) / f"Figure_S02{panel}.{suffix}",
                            dpi=220, bbox_inches="tight", facecolor="white")
        plt.close(fig)
    return tables


def build_mix_bundle(root, *, frozen_metrics_csv=None):
    """Use the existing S2 preprocessing once, then freeze all panel inputs.

    Pass Figure_S02_all_cell_metrics.csv.gz from an accepted run to avoid
    reassigning lineages. Otherwise the existing S2 classifier is run once.
    """
    import anndata as ad
    s2 = _mix_producers()  # Direct repository producers; no notebook/archive dependency.
    from .figure_s02_parse_mix_data import (read_expression_and_qc,
        assign_transcriptomic_cell_types, _canonical_cell_ids, _guide_pool)
    rna, paths = read_expression_and_qc(Path(root))
    _count_check(rna.X, "Mixed-cell RNA X")
    if frozen_metrics_csv is not None:
        frozen = pd.read_csv(frozen_metrics_csv, index_col=0)
        frozen.index = _canonical_cell_ids(frozen.index)
        if not frozen.index.is_unique or set(frozen.index) != set(rna.obs_names):
            raise ValueError("Frozen S2 metrics and QC-passing RNA cells do not match")
        frozen = frozen.loc[rna.obs_names]
        if "cell_type" not in frozen or "sample" not in frozen:
            raise ValueError("Frozen S2 metrics lack cell_type/sample")
        if not np.array_equal(frozen["sample"].astype(str), rna.obs["sample"].astype(str)):
            raise ValueError("Frozen primer labels differ from source metadata")
        cell_types = frozen["cell_type"]
        cluster_annotation = pd.DataFrame({"note": ["Lineages retained from supplied accepted-run metrics"]})
        annotation_source = file_checksums(frozen_metrics_csv)
    else:
        cell_types, cluster_annotation = assign_transcriptomic_cell_types(rna)
        annotation_source = "Existing S2 classifier executed once during export; public plotting uses frozen labels"
    metadata = rna.obs.copy()
    metadata["cell_type"] = cell_types.reindex(metadata.index)
    pairs, guide_diagnostics = s2.read_guide_pairs(
        paths["guide_counts"], paths["guide_reference"], rna.obs_names)
    metrics = s2.summarize_guide_pairs(pairs, metadata)
    if frozen_metrics_csv is not None:
        for key in METRIC_COLUMNS:
            if key not in frozen:
                raise ValueError(f"Frozen metrics lack {key}")
            np.testing.assert_allclose(metrics[key].to_numpy(float), frozen[key].to_numpy(float),
                                       rtol=1e-6, atol=1e-7, equal_nan=True, err_msg=key)
    rna.obs = metrics.copy()
    reference = pd.read_csv(paths["guide_reference"])
    reference = reference.loc[reference["gene_name"].notna()].copy()
    reference.index = pd.Index(reference["gene_name"].astype(str), name="guide")
    if not reference.index.is_unique:
        raise ValueError("Guide-reference IDs must be unique")
    reference["guide_pool"] = _guide_pool(reference["gene_name"]).to_numpy()
    rows = rna.obs_names.get_indexer(pairs["cell"])
    cols = reference.index.get_indexer(pairs["guide"])
    if (rows < 0).any() or (cols < 0).any():
        raise ValueError("Guide pairs do not match RNA cell IDs or guide reference")
    counts = sparse.coo_matrix((pairs["count"].to_numpy(), (rows, cols)),
                               shape=(rna.n_obs, len(reference))).tocsr()
    counts.sum_duplicates()
    counts.eliminate_zeros()
    crispr = ad.AnnData(X=counts, obs=metrics.copy(), var=_safe_frame(reference))
    crispr.uns["geo_X_semantics"] = "Original guide UMIs, summed by cell/guide; not thresholded"
    crispr.uns["panel_n_min_umis_per_cell_guide"] = 3
    rna.uns["geo_X_semantics"] = "QC-passing, unnormalized Split-pipe RNA counts"
    saturation = s2.read_saturation(paths["root"])
    # No live path is needed to use these embedded numeric reports.
    saturation["source_file"] = saturation["source_file"].map(lambda p: Path(p).name)
    tables = _mix_tables_and_plots(metrics, saturation)
    bundle = _assemble(
        {"rna": rna, "crispr": crispr}, "ipsc_k562_mix",
        {"input_files": {k: file_checksums(p) for k, p in paths.items() if p.is_file()},
         "scope": ["Figure_S02i", "Figure_S02j", "Figure_S02l", "Figure_S02m", "Figure_S02n"],
         "annotation_source": annotation_source, "guide_diagnostics": guide_diagnostics,
         "qc": {"n_genes_by_counts": "<5000", "total_counts": "<5000", "pct_counts_mt": "<40"},
         "markers": s2.MARKERS, "plotting_reclusters": False}, tables,
        _sample_table(metrics, ("sample", "cell_type")))
    bundle.uns["saturation"] = _safe_frame(saturation)
    bundle.uns["cluster_annotation"] = _safe_frame(cluster_annotation)
    bundle.uns["manuscript_median_diagnostics"] = _safe_frame(s2.manuscript_median_diagnostics(metrics))
    bundle.uns["raw_saturation_reports_json"] = _json({
        primer: paths[f"saturation_{primer}"].read_text() for primer in s2.PRIMER_ORDER})
    return bundle


def _metrics_from_sparse_guides(bundle) -> pd.DataFrame:
    """Independent sparse-axis reconstruction; never densify cells x guides."""
    rna, guides = bundle.mod["rna"], bundle.mod["crispr"]
    x = sparse.csr_matrix(guides.X)
    metrics = rna.obs.copy()
    metrics["total_guide_umis"] = np.asarray(x.sum(axis=1)).ravel()
    metrics["top_guide_umis"] = x.max(axis=1).toarray().ravel()
    for lineage in ("K562", "hiPSC"):
        columns = guides.var["guide_pool"].astype(str).eq(lineage).to_numpy()
        part = x[:, columns].tocsr()
        metrics[f"{lineage}_guide_umis"] = np.asarray(part.sum(axis=1)).ravel()
        part.data[part.data < 3] = 0
        part.eliminate_zeros()
        metrics[f"{lineage}_guide_umis_ge3"] = np.asarray(part.sum(axis=1)).ravel()
    metrics["matched_pool_guide_umis"] = np.select(
        [metrics["cell_type"].eq("K562"), metrics["cell_type"].eq("hiPSC")],
        [metrics["K562_guide_umis"], metrics["hiPSC_guide_umis"]], default=np.nan)
    total = (metrics["K562_guide_umis"] + metrics["hiPSC_guide_umis"]).to_numpy(float)
    metrics["guide_pool_specificity"] = np.divide(
        metrics["matched_pool_guide_umis"].to_numpy(float), total,
        out=np.full(len(metrics), np.nan), where=total > 0)
    return metrics


@contextmanager
def _local_hippo_inputs(bundle, parent=None):
    """Adapt the deposited object to the existing path-based producers locally."""
    with tempfile.TemporaryDirectory(prefix="hippo_geo_", dir=parent) as directory:
        root = Path(directory)
        query, projected = root / "rna.h5ad", root / "cnmf.h5ad"
        bundle.mod["rna"].write_h5ad(query)
        bundle.mod["cnmf"].write_h5ad(projected)
        yield query, projected


def validate_submission(path, *, render_dir=None, temp_parent=None) -> dict[str, Any]:
    """Read only the submitted H5MU, recompute panel inputs and compare references.

    This function never resolves an Oak path, calls a source-data preflight,
    refits cNMF, or reruns the mixed-cell classifier.
    """
    import mudata as md
    path = Path(path).expanduser().resolve()
    bundle = md.read_h5mu(path)
    if bundle.uns.get("geo_schema") != SCHEMA:
        raise ValueError("Unsupported/missing processed-object schema")
    dataset = bundle.uns["dataset"]
    expected = {"hippo_timecourse": {"rna", "cnmf"},
                "ipsc_k562_mix": {"rna", "crispr"}}
    if dataset not in expected or set(bundle.mod) != expected[dataset]:
        raise ValueError("Incorrect dataset/modality combination")
    rna = bundle.mod["rna"]
    for name, obj in bundle.mod.items():
        if not obj.obs_names.equals(rna.obs_names) or not obj.obs_names.is_unique or not obj.var_names.is_unique:
            raise ValueError(f"{name}: duplicate or misaligned cell/feature IDs")
    counts = {"rna": _count_check(rna.X, "RNA X")}
    hashes = json.loads(bundle.uns["matrix_sha256_json"])
    for name, digest in hashes.items():
        parts = name.split("/")
        obj = bundle.mod[parts[0]]
        value = obj.X if parts[1] == "X" else getattr(obj, parts[1])["/".join(parts[2:])]
        if _matrix_hash(value) != digest:
            raise ValueError(f"Matrix checksum mismatch: {name}")
    if render_dir is not None:
        render_dir = Path(render_dir).expanduser().resolve()
        render_dir.mkdir(parents=True, exist_ok=True)
    if dataset == "hippo_timecourse":
        from .hippo_projection import validate_projected
        from .figure3_panels import _hippo_aligned_metadata, _hippo_panel_mask, FIGURE3F_D3_CONDITIONS
        from .figure4_panels import FIGURE4C_D12_CONDITIONS
        projection = bundle.mod["cnmf"]
        obs = _hippo_aligned_metadata(rna, projection)
        for panel, conditions in (("Figure_3f", FIGURE3F_D3_CONDITIONS), ("Figure_4c", FIGURE4C_D12_CONDITIONS)):
            if not np.array_equal(rna.obs[f"used_in_{panel}"].to_numpy(bool), _hippo_panel_mask(obs, conditions)):
                raise ValueError(f"{panel}: stored inclusion mask disagrees with original producer")
        with _local_hippo_inputs(bundle, temp_parent) as (query, projected):
            validate_projected(query, projected, strict_historical=False)
            tables, _ = _hippo_tables_and_plots(query, projected, render_dir)
    else:
        counts["crispr"] = _count_check(bundle.mod["crispr"].X, "CRISPR X")
        metrics = _metrics_from_sparse_guides(bundle)
        for key in METRIC_COLUMNS:
            np.testing.assert_allclose(metrics[key].to_numpy(float), rna.obs[key].to_numpy(float),
                                       rtol=1e-6, atol=1e-7, equal_nan=True, err_msg=key)
        tables = _mix_tables_and_plots(metrics, bundle.uns["saturation"], render_dir)
    reference = bundle.uns["reference_tables"]
    if set(tables) != set(reference):
        raise ValueError("Panel reference-table set is incomplete")
    for name, table in tables.items():
        _same_table(table, reference[name], name)
        if render_dir is not None:
            table.to_csv(render_dir / f"{name}.csv", index=True)
    report = {"schema": SCHEMA, "dataset": dataset, "status": "PASS",
              "n_cells": rna.n_obs, "modalities": {k: list(v.shape) for k, v in bundle.mod.items()},
              "count_checks": counts, "verified_matrix_hashes": len(hashes),
              "reproduced_reference_tables": sorted(tables),
              "input_file": file_checksums(path),
              "original_source_files_required": False,
              "meaning": "Stored-data roundtrip and numerical panel equivalence; not a GEO submission or independent manuscript audit"}
    if render_dir is not None:
        (render_dir / "validation.json").write_text(_json(report) + "\n")
    return report


def _publish(bundle, output_path, *, overwrite=False, temp_parent=None):
    output = Path(output_path).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not overwrite:
        raise FileExistsError(f"{output} exists; set overwrite=True only to replace a prior export")
    with tempfile.NamedTemporaryFile(dir=output.parent, prefix=f".{output.stem}.", suffix=".h5mu", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        bundle.write_h5mu(temporary, compression="gzip")
        del bundle
        gc.collect()
        # Never publish an unvalidated half-written file as the submission artifact.
        report = validate_submission(temporary, temp_parent=temp_parent)
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)
    report["input_file"] = file_checksums(output)
    output.with_suffix(".validation.json").write_text(_json(report) + "\n")
    for algorithm in ("md5", "sha256"):
        output.with_suffix(output.suffix + f".{algorithm}").write_text(
            f"{report['input_file'][algorithm]}  {output.name}\n")
    return report


def _check_destination(output_path, inputs=(), *, overwrite=False):
    path = Path(output_path).expanduser().resolve()
    if path.suffix != ".h5mu":
        raise ValueError("Use an .h5mu filename for the integrated MuData object")
    if path in {Path(p).expanduser().resolve() for p in inputs}:
        raise ValueError("Output must not overwrite any source input")
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists; validate it or deliberately set overwrite=True")


def export_hippo(query_path, projection_path, spectra_path, output_path, *, overwrite=False, temp_parent=None):
    _check_destination(output_path, (query_path, projection_path, spectra_path), overwrite=overwrite)
    return _publish(build_hippo_bundle(query_path, projection_path, spectra_path),
                    output_path, overwrite=overwrite, temp_parent=temp_parent)


def export_parse_mix(root, output_path, *, frozen_metrics_csv=None, overwrite=False):
    _check_destination(output_path, overwrite=overwrite)
    return _publish(build_mix_bundle(root, frozen_metrics_csv=frozen_metrics_csv),
                    output_path, overwrite=overwrite)
