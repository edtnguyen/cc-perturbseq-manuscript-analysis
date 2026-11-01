"""Maintainer exports of small plotting inputs; never publish single-cell matrices.

Use the accepted Hippo H5AD pair and S2 numerical outputs, or extract only the
saved plotting tables/obs from existing GEO H5MU files. These input routes are
explicit; missing inputs never trigger a different analysis or reclustering.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

from .geo_plot_data import HIPPO_TABLES, load_plot_bundle, write_plot_bundle
from .igvf_download import digest_file

DATASETS = ("hippo_timecourse", "ipsc_k562_mix")
MIX_FILES = {
    "cell_metrics": "Figure_S02_all_cell_metrics.csv.gz",
    "saturation": "panel_m/Figure_S02m_transcript_saturation_source.csv",
    "cluster_annotation": "Figure_S02_transcriptomic_cluster_annotation.csv",
    "median_diagnostics": "Figure_S02_manuscript_median_diagnostics.csv",
}


def resolve_sources(config, *, datasets=DATASETS, hippo_query=None,
                    hippo_projection=None, mixed_figure_dir=None,
                    hippo_h5mu=None, mixed_h5mu=None):
    """Select source files explicitly; validate all files before any computation."""
    defaults = config.get("compact_plot_sources", {})
    sources = {}
    missing = []
    if hippo_h5mu and (hippo_query or hippo_projection):
        raise ValueError("Select the Hippo H5AD pair OR saved H5MU tables, not both.")
    if mixed_h5mu and mixed_figure_dir:
        raise ValueError("Select saved S2 outputs OR saved H5MU tables, not both.")
    for dataset in datasets:
        if dataset not in DATASETS:
            raise ValueError(dataset)
        obj = hippo_h5mu if dataset == "hippo_timecourse" else mixed_h5mu
        if obj:
            spec = {"mode": "saved_h5mu_tables", "files": {"processed": obj}}
        elif dataset == "hippo_timecourse":
            spec = {"mode": "hippo_h5ad_pair", "files": {
                "query": hippo_query or defaults.get("hippo_query"),
                "projection": hippo_projection or defaults.get("hippo_projection")}}
        else:
            directory = mixed_figure_dir or defaults.get("mixed_figure_dir")
            spec = {"mode": "saved_s2_outputs", "files": {
                key: Path(directory) / rel if directory else None
                for key, rel in MIX_FILES.items()}}
            if directory:
                provenance = Path(directory) / "Figure_S02_parse_mix_provenance.json"
                if provenance.is_file():
                    spec["files"]["run_provenance"] = provenance
        for name, raw in spec["files"].items():
            path = Path(raw).expanduser().resolve() if raw else None
            if path is None or not path.is_file() or path.stat().st_size == 0:
                missing.append(f"{dataset} {name}: {path or 'not configured'}")
            else:
                spec["files"][name] = path
        sources[dataset] = spec
    if missing:
        raise FileNotFoundError("Missing inputs for compact plotting tables:\n  "
            + "\n  ".join(missing) + "\nUse --hippo-query/--hippo-projection and "
            "--mixed-figure-dir, or explicitly select existing --hippo-h5mu/--mixed-h5mu. "
            "No integrated H5MU needs to be created. Use --dataset to export one experiment.")
    return sources


def _table_fingerprint(tables):
    """Digest selected numerical tables, not an unread full H5MU file."""
    h = hashlib.sha256()
    for key, frame in sorted(tables.items()):
        h.update(key.encode())
        h.update(frame.to_csv(sep="\t", float_format="%.17g").encode())
    return h.hexdigest()


def _validate_mixed(tables):
    import numpy as np
    from .geo_submission import METRIC_COLUMNS
    cells = tables["cell_metrics"]
    required = {"sample", "cell_type", "total_counts", *METRIC_COLUMNS}
    if required - set(cells) or not cells.index.is_unique or cells.empty:
        raise ValueError("Use the COMPLETE Figure_S02_all_cell_metrics table with unique cell IDs, not a filtered panel source table.")
    if cells[["sample", "cell_type"]].isna().any().any():
        raise ValueError("Frozen sample/lineage labels are missing.")
    for key in (required - {"sample", "cell_type", "guide_pool_specificity", "matched_pool_guide_umis"}):
        a = cells[key].to_numpy(float)
        if not (np.isfinite(a) & (a >= 0)).all():
            raise ValueError(f"Invalid mixed-cell metric: {key}")
    for lineage in ("K562", "hiPSC"):
        if not cells[f"{lineage}_guide_umis_ge3"].le(cells[f"{lineage}_guide_umis"]).all():
            raise ValueError("Thresholded pool counts exceed unthresholded counts.")
    expected = np.select([cells.cell_type.eq("K562"), cells.cell_type.eq("hiPSC")],
        [cells.K562_guide_umis, cells.hiPSC_guide_umis], default=np.nan)
    np.testing.assert_allclose(cells.matched_pool_guide_umis, expected, equal_nan=True)
    total = (cells.K562_guide_umis + cells.hiPSC_guide_umis).to_numpy(float)
    specificity = np.divide(expected, total, out=np.full(len(cells), np.nan), where=total > 0)
    np.testing.assert_allclose(cells.guide_pool_specificity, specificity, equal_nan=True, rtol=1e-6, atol=1e-7)
    sat = tables["saturation"]
    columns = {"primer", "mean_assigned_reads_per_cell", "median_transcript_umis_per_cell"}
    if columns - set(sat) or set(sat.primer) != {"CS1", "DC1", "Dual", "no_guide"}:
        raise ValueError("All four numerical saturation curves are required.")
    for key in columns - {"primer"}:
        if not (np.isfinite(sat[key].to_numpy(float)) & (sat[key].to_numpy(float) >= 0)).all():
            raise ValueError("Invalid saturation values.")


def _saved_h5mu_tables(path, dataset):
    """Read only small named HDF5 groups, never RNA/guide/cNMF count matrices."""
    import h5py
    from anndata.io import read_elem
    from .geo_submission import SCHEMA
    with h5py.File(path, "r") as f:
        def read(key):
            if key not in f:
                raise ValueError(f"Saved plotting information absent: {key}; use source-file export or full GEO regeneration.")
            return read_elem(f[key])
        if read("uns/geo_schema") != SCHEMA or read("uns/dataset") != dataset:
            raise ValueError("Wrong submitted dataset/schema.")
        if dataset == "hippo_timecourse":
            tables = {key: read("uns/reference_tables/" + key) for key in sorted(HIPPO_TABLES)}
            # The submission writer stringifies table row IDs. Restore only the
            # known 0..N-1 summary indices; condition IDs on dot matrices stay intact.
            import pandas as pd
            for key, frame in tables.items():
                if key.endswith(("_genes", "_programs")):
                    if list(map(str, frame.index)) != [str(i) for i in range(len(frame))]:
                        raise ValueError("Hippo saved summary row order is not the submitted order.")
                    frame.index = pd.RangeIndex(len(frame))
            metadata = json.loads(read("uns/provenance_json"))["panel_metadata"]
        else:
            tables = {"cell_metrics": read("mod/rna/obs"),
                      "saturation": read("uns/saturation"),
                      "cluster_annotation": read("uns/cluster_annotation"),
                      "median_diagnostics": read("uns/manuscript_median_diagnostics")}
            metadata = {"n_cells": len(tables["cell_metrics"]), "lineages": "frozen in the submitted object",
                        "cell_guide_threshold_for_mixing": 3}
    return tables, metadata


def export_one(spec, directory, dataset):
    """Write the existing plotting-bundle schema, compatible with all notebooks."""
    import pandas as pd
    mode, files = spec["mode"], spec["files"]
    if mode == "saved_h5mu_tables":
        tables, metadata = _saved_h5mu_tables(files["processed"], dataset)
    elif mode == "hippo_h5ad_pair" and dataset == "hippo_timecourse":
        from .geo_submission import _hippo_tables_and_plots
        tables, metadata = _hippo_tables_and_plots(files["query"], files["projection"])
    elif mode == "saved_s2_outputs" and dataset == "ipsc_k562_mix":
        tables = {}
        for key in MIX_FILES:
            # Cell IDs must remain strings, including possible leading zeroes.
            if key == "cell_metrics":
                index_col = pd.read_csv(files[key], nrows=0).columns[0]
                frame = pd.read_csv(files[key], index_col=0, dtype={index_col: str}, float_precision="round_trip", keep_default_na=False, na_values=[""])
            else:
                frame = pd.read_csv(files[key], float_precision="round_trip", keep_default_na=False, na_values=[""])
            tables[key] = frame
        metadata = {"n_cells": len(tables["cell_metrics"]), "lineages": "frozen in accepted numerical outputs",
                    "cell_guide_threshold_for_mixing": 3}
        if "run_provenance" in files:
            recorded = json.loads(Path(files["run_provenance"]).read_text())
            expected_n = recorded.get("n_cells_after_qc")
            if expected_n is not None and int(expected_n) != len(tables["cell_metrics"]):
                raise ValueError("Complete S2 metrics disagree with the recorded QC cell count.")
    else:
        raise ValueError(f"Unsupported compact export mode: {mode}")
    if dataset == "ipsc_k562_mix":
        _validate_mixed(tables)
    for table in tables.values():
        if "source_file" in table:
            table["source_file"] = table.source_file.map(lambda p: Path(str(p)).name)
    source = {"mode": mode, "filenames": {k: Path(v).name for k, v in files.items()},
        "sha256": _table_fingerprint(tables), "sha256_scope": "selected uncompressed plotting tables",
        "recomputed_from_modalities": mode == "hippo_h5ad_pair",
        "refit_cnmf": False, "reclustered_cells": False,
        "raw_matrix_integrity_verified": False,
        "note": "Saved-table export is not a full verification of source RNA/guide matrices. Optional GEO regeneration performs that verification."}
    return write_plot_bundle(tables, directory, dataset, metadata=metadata, source=source)


def export_all(sources, root, *, apply=False, overwrite=False):
    """Stage all requested experiments before replacing any published bundle."""
    root = Path(root).expanduser().resolve()
    parent = root / "data/reproduction/geo_plots"
    for dataset, spec in sources.items():
        print(f"{dataset}: {spec['mode']}", flush=True)
        for name, path in spec["files"].items():
            print(f"  {name}: {path}", flush=True)
        if (parent / dataset).exists() and not overwrite:
            raise FileExistsError(f"{parent / dataset} exists. Use --overwrite to back it up first.")
    if not apply:
        return {"status": "dry_run", "datasets": list(sources), "single_cell_objects_uploaded": False}
    parent.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".compact-plots-", dir=parent.parent))
    reports = {}
    moved = []
    backups = []
    try:
        for dataset, spec in sources.items():
            manifest = export_one(spec, staging / dataset, dataset)
            reports[dataset] = {"tables": len(manifest["tables"]),
                "size_bytes": sum(p.stat().st_size for p in (staging / dataset).iterdir()),
                "largest_table_bytes": max(v["size_bytes"] for v in manifest["tables"].values())}
        # No .h5mu/.h5ad can enter the publication directory by this route.
        if any(p.suffix in {".h5mu", ".h5ad"} for p in staging.rglob("*")):
            raise RuntimeError("Unexpected single-cell object in compact publication.")
        parent.mkdir(exist_ok=True)
        for dataset in sources:
            dest = parent / dataset
            if dest.exists():
                backup = root / "data" / (dataset + ".backup-" + datetime.now().strftime("%Y%m%d-%H%M%S-%f"))
                dest.rename(backup)
                backups.append((backup, dest))
            (staging / dataset).rename(dest)
            moved.append(dest)
            load_plot_bundle(dest, dataset)
        result = {"status": "exported", "datasets": reports, "single_cell_objects_uploaded": False}
        print(json.dumps(result, indent=2), flush=True)
        print("Commit only data/reproduction/geo_plots/. No files were staged or pushed.", flush=True)
        return result
    except BaseException:
        for dest in reversed(moved):
            shutil.rmtree(dest)
        for backup, dest in reversed(backups):
            backup.rename(dest)
        raise
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    p.add_argument("--source-config", type=Path)
    p.add_argument("--dataset", choices=("all", *DATASETS), default="all")
    p.add_argument("--hippo-query", type=Path)
    p.add_argument("--hippo-projection", type=Path)
    p.add_argument("--mixed-figure-dir", type=Path)
    p.add_argument("--hippo-h5mu", type=Path, help="Extract only saved tables from an existing integrated object")
    p.add_argument("--mixed-h5mu", type=Path, help="Extract only obs/summary tables from an existing integrated object")
    p.add_argument("--apply", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args(argv)
    cfg = args.source_config or args.repo / "config/publication_source_files.json"
    config = json.loads(cfg.read_text())
    sources = resolve_sources(config, datasets=DATASETS if args.dataset == "all" else (args.dataset,),
        hippo_query=args.hippo_query, hippo_projection=args.hippo_projection,
        mixed_figure_dir=args.mixed_figure_dir, hippo_h5mu=args.hippo_h5mu, mixed_h5mu=args.mixed_h5mu)
    return export_all(sources, args.repo, apply=args.apply, overwrite=args.overwrite)
