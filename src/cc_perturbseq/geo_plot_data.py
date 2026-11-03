"""Versioned plotting tables for the two GEO experiments.

Default figure execution reads only these tables. Regeneration reads the
submitted integrated object, invokes the existing numerical producers, and
writes a separate cache without changing the published reference bundle.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import tempfile
import warnings

from .igvf_download import digest_file

SCHEMA = "ccperturbseq.plot_inputs.v1"
DATASETS = ("hippo_timecourse", "ipsc_k562_mix")
HIPPO_PANELS = ("Figure_3f", "Figure_4c")
MIX_TABLES = {"cell_metrics", "saturation", "cluster_annotation", "median_diagnostics"}
HIPPO_TABLES = {f"{p}_{part}" for p in HIPPO_PANELS for part in
                ("genes", "programs", "genes_dot_color", "genes_dot_size", "programs_dot_color", "programs_dot_size")}


def _json(value):
    from .geo_submission import _json as encode
    return encode(value)


def _write_table(frame, path):
    import pandas as pd
    if not frame.index.is_unique or not frame.columns.is_unique or frame.columns.nlevels != 1 or frame.index.nlevels != 1:
        raise ValueError("Precomputed tables require unique one-dimensional row and column identifiers.")
    if not all(isinstance(c, str) for c in frame.columns):
        raise ValueError("Plot table column names must be strings.")
    spec = {"columns": list(frame.columns), "index_name": frame.index.name,
            "index_dtype": str(frame.index.dtype), "dtypes": {}, "categoricals": {}}
    for column in frame:
        if isinstance(frame[column].dtype, pd.CategoricalDtype):
            spec["categoricals"][column] = {"categories": frame[column].cat.categories.tolist(), "ordered": frame[column].cat.ordered}
            spec["dtypes"][column] = "object"
        else:
            spec["dtypes"][column] = str(frame[column].dtype)
    if isinstance(frame.index, pd.CategoricalIndex):
        spec["index_dtype"] = "object"
        spec["index_categories"] = {"categories": frame.index.categories.tolist(), "ordered": frame.index.ordered}
    frame.to_csv(path, sep="\t", index=True, index_label="__row_id__", float_format="%.17g", na_rep="__CCPS_MISSING__",
                 compression={"method": "gzip", "mtime": 0})
    return spec


def _read_table(path, spec):
    import pandas as pd
    # Explicit text types preserve IDs such as 0001 and strings such as NA.
    dtype = {"__row_id__": str, **{k: str for k, v in spec["dtypes"].items() if v in ("object", "string")}}
    table = pd.read_csv(path, sep="\t", index_col=0, dtype=dtype,
                        keep_default_na=False, na_values=["__CCPS_MISSING__"], float_precision="round_trip")
    if table.columns.tolist() != spec["columns"]:
        raise ValueError(f"Precomputed table columns changed: {path}")
    for column, dt in spec["dtypes"].items():
        if column not in spec["categoricals"] and dt not in {"object", "string"}:
            table[column] = table[column].astype(dt)
    for column, cat in spec["categoricals"].items():
        table[column] = pd.Categorical(table[column], **cat)
    if "index_categories" in spec:
        table.index = pd.CategoricalIndex(table.index, **spec["index_categories"])
    elif spec["index_dtype"] not in {"object", "string"}:
        table.index = table.index.astype(spec["index_dtype"])
    table.index.name = spec["index_name"]
    return table


def load_plot_bundle(directory, dataset, *, tables=True):
    """Validate manifest, required table set, sizes and content hashes."""
    directory = Path(directory).resolve()
    manifest_path = directory / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing published plotting-data bundle: {manifest_path}. Retrieve data/reproduction/geo_plots, or explicitly select geo/local_file regeneration. No original-filesystem fallback is used.")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != SCHEMA or manifest.get("dataset") != dataset:
        raise ValueError(f"Wrong plotting-data schema/dataset at {directory}")
    expected = HIPPO_TABLES if dataset == "hippo_timecourse" else MIX_TABLES if dataset == "ipsc_k562_mix" else None
    if expected is None or set(manifest["tables"]) != expected:
        raise ValueError(f"Incomplete precomputed table set for {dataset}")
    loaded = {}
    for key, record in manifest["tables"].items():
        relative = Path(record["file"])
        path = (directory / relative).resolve()
        if relative.is_absolute() or ".." in relative.parts or not path.is_relative_to(directory):
            raise ValueError("Precomputed data path escapes its bundle.")
        if not path.is_file() or path.stat().st_size != record["size_bytes"] or digest_file(path, "sha256") != record["sha256"]:
            raise ValueError(f"Missing or checksum-mismatched plotting table: {key}")
        if tables:
            loaded[key] = _read_table(path, record["table_schema"])
    return loaded, manifest


def write_plot_bundle(tables, directory, dataset, *, metadata, source, overwrite=False):
    """Stage atomically and read back every table before publishing the directory."""
    import pandas as pd
    destination = Path(directory).resolve()
    if destination.exists() and not overwrite:
        raise FileExistsError(f"Plot bundle already exists: {destination}. Use overwrite=True to preserve a backup and replace it.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".plot-data-", dir=destination.parent))
    try:
        records = {}
        for key, frame in tables.items():
            if not key.replace("_", "").isalnum():
                raise ValueError("Unsafe plot-table key.")
            path = staging / f"{key}.tsv.gz"
            spec = _write_table(frame, path)
            if path.stat().st_size >= 95 * 1024 ** 2:
                raise ValueError(f"Plot table {key} exceeds the 95 MiB Git publication guard.")
            records[key] = {"file": path.name, "size_bytes": path.stat().st_size,
                "sha256": digest_file(path, "sha256"), "table_schema": spec}
        manifest = {"schema": SCHEMA, "dataset": dataset, "created_utc": datetime.now(timezone.utc).isoformat(),
                    "source": source, "metadata": metadata, "tables": records}
        (staging / "manifest.json").write_text(_json(manifest) + "\n")
        readback, _ = load_plot_bundle(staging, dataset)
        for key, frame in tables.items():
            pd.testing.assert_frame_equal(readback[key], frame, check_dtype=False, check_index_type=False,
                                          check_categorical=False, rtol=1e-12, atol=1e-12)
        if destination.exists():
            backup = destination.with_name(destination.name + ".backup-" + datetime.now().strftime("%Y%m%d-%H%M%S-%f"))
            destination.rename(backup)
        staging.rename(destination)
        return manifest
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def build_plot_bundle(processed_path, directory, dataset, *, overwrite=False):
    """Recalculate from the actual deposited modalities, not its saved plot tables."""
    import mudata as md
    import numpy as np
    from . import geo_submission as geo
    processed_path = Path(processed_path).expanduser().resolve()
    if not processed_path.is_file():
        raise FileNotFoundError(f"Submitted processed file is missing: {processed_path}")
    # The submitted release has distinct RNA and cNMF/guide feature spaces.
    # An arbitrary expression-only H5AD cannot supply these additional data.
    import h5py
    with h5py.File(processed_path, "r") as handle:
        if "mod" not in handle:
            raise ValueError("This study's submitted object must be the integrated H5MU (RNA plus cnmf or crispr), not an expression-only H5AD. Required guide/saturation/projection information cannot be inferred.")
    bundle = md.read_h5mu(processed_path)
    modalities = {"hippo_timecourse": {"rna", "cnmf"}, "ipsc_k562_mix": {"rna", "crispr"}}
    if dataset not in modalities or bundle.uns.get("geo_schema") != geo.SCHEMA or bundle.uns.get("dataset") != dataset or set(bundle.mod) != modalities[dataset]:
        raise ValueError(f"Not the submitted {dataset} integrated schema ({geo.SCHEMA}).")
    rna = bundle.mod["rna"]
    for name, obj in bundle.mod.items():
        if not obj.obs_names.equals(rna.obs_names) or not obj.obs_names.is_unique or not obj.var_names.is_unique:
            raise ValueError(f"Duplicate or misaligned identifiers in {name}.")
    verified_hashes = 0
    for name, expected in json.loads(bundle.uns.get("matrix_sha256_json", "{}")).items():
        parts = name.split("/")
        obj = bundle.mod[parts[0]]
        matrix = obj.X if parts[1] == "X" else getattr(obj, parts[1])["/".join(parts[2:])]
        if geo._matrix_hash(matrix) != expected:
            raise ValueError(f"Submitted matrix checksum mismatch: {name}")
        verified_hashes += 1
    if dataset == "hippo_timecourse":
        from .figure3_panels import _hippo_aligned_metadata, _hippo_panel_mask, FIGURE3F_D3_CONDITIONS
        from .figure4_panels import FIGURE4C_D12_CONDITIONS
        aligned = _hippo_aligned_metadata(rna, bundle.mod["cnmf"])
        for panel, conditions in (("Figure_3f", FIGURE3F_D3_CONDITIONS), ("Figure_4c", FIGURE4C_D12_CONDITIONS)):
            key = f"used_in_{panel}"
            if key in rna.obs and not np.array_equal(rna.obs[key].to_numpy(bool), _hippo_panel_mask(aligned, conditions)):
                raise ValueError(f"Stored {panel} cell mask differs from the existing QC rule.")
        # Avoid system TMPDIR (which may point to scratch) for regeneration.
        temp_parent = Path(directory).resolve().parent
        temp_parent.mkdir(parents=True, exist_ok=True)
        with geo._local_hippo_inputs(bundle, parent=temp_parent) as (query, projection):
            tables, metadata = geo._hippo_tables_and_plots(query, projection)
        comparison_tables = tables
    else:
        metrics = geo._metrics_from_sparse_guides(bundle)
        for key in geo.METRIC_COLUMNS:
            np.testing.assert_allclose(metrics[key].to_numpy(float), rna.obs[key].to_numpy(float),
                                       rtol=1e-6, atol=1e-7, equal_nan=True, err_msg=key)
        tables = {"cell_metrics": metrics, "saturation": bundle.uns["saturation"].copy(),
                  "cluster_annotation": bundle.uns["cluster_annotation"].copy(),
                  "median_diagnostics": bundle.uns["manuscript_median_diagnostics"].copy()}
        metadata = {"n_cells": len(metrics), "lineages": "frozen in the submitted object", "cell_guide_threshold_for_mixing": 3}
        comparison_tables = geo._mix_tables_and_plots(metrics, tables["saturation"]) if "reference_tables" in bundle.uns else {}
    references = bundle.uns.get("reference_tables", {})
    checked = []
    if references:
        if set(references) != set(comparison_tables):
            raise ValueError("Incomplete embedded reference table set in submitted object.")
        for key in references:
            geo._same_table(comparison_tables[key], references[key], key)
            checked.append(key)
    # Provenance fields are descriptive; never dereference stored source paths.
    for table in tables.values():
        for col in table.columns:
            if str(col) == "source_file":
                table[col] = table[col].map(lambda x: Path(str(x)).name)
    source = {"filename": processed_path.name, "sha256": digest_file(processed_path, "sha256"),
              "geo_schema": geo.SCHEMA, "verified_embedded_matrix_hashes": verified_hashes,
              "matched_embedded_reference_tables": checked, "versions": geo._versions(),
              "recomputed_from_modalities": True, "refit_cnmf": False, "reclustered_cells": False}
    return write_plot_bundle(tables, directory, dataset, metadata=metadata, source=source, overwrite=overwrite)


def load_mixed_plot_data(directory):
    from .figure_s02_parse_mix_data import ParseMixResult
    tables, manifest = load_plot_bundle(directory, "ipsc_k562_mix")
    return ParseMixResult(tables["cell_metrics"], tables["cluster_annotation"], tables["saturation"],
        tables["median_diagnostics"], {"plot_bundle": str(Path(directory).resolve()),
        "source_sha256": manifest["source"]["sha256"],
        "source_sha256_scope": manifest["source"].get("sha256_scope", "submitted processed file"), "lineages": "frozen", "original_source_files_required": False})


def prepare_hippo_dotplots(directory, panel):
    """Render with the same Scanpy class, using explicit published dot matrices."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from .figure3_panels import FIGURE3F_GENES, FIGURE3F_PROGRAMS
    from .figure4_panels import FIGURE4C_GENES, FIGURE4C_PROGRAMS
    tables, manifest = load_plot_bundle(directory, "hippo_timecourse")
    if panel not in HIPPO_PANELS:
        raise ValueError(panel)
    genes = FIGURE3F_GENES if panel == "Figure_3f" else FIGURE4C_GENES
    programs = FIGURE3F_PROGRAMS if panel == "Figure_3f" else FIGURE4C_PROGRAMS
    plots = []
    for kind, features in (("genes", {k: list(v) for k, v in genes.items()}), ("programs", list(programs))):
        colors = tables[f"{panel}_{kind}_dot_color"]
        sizes = tables[f"{panel}_{kind}_dot_size"]
        names = [g for v in features.values() for g in v] if isinstance(features, dict) else features
        if colors.columns.tolist() != names or not sizes.index.equals(colors.index) or not sizes.columns.equals(colors.columns):
            raise ValueError("Hippo dot matrix labels/order do not match the panel specification.")
        if not np.isfinite(colors.to_numpy(float)).all() or not np.isfinite(sizes.to_numpy(float)).all() or not ((sizes.to_numpy(float) >= 0) & (sizes.to_numpy(float) <= 1)).all():
            raise ValueError("Invalid Hippo dot colors/fractions.")
        categories = colors.index.astype(str).tolist()
        # A layout-only container: all scientific values are supplied explicitly
        # below. No cell-level observations or statistics are reconstructed here.
        layout = ad.AnnData(np.zeros((len(categories), len(names)), dtype=np.float32),
            obs=pd.DataFrame({"condition_name": pd.Categorical(categories, categories=categories, ordered=True)}, index=categories),
            var=pd.DataFrame(index=names))
        color = colors.copy(); color.index = pd.Index(categories)
        size = sizes.copy(); size.index = pd.Index(categories)
        plot = sc.pl.DotPlot(layout, features, groupby="condition_name", categories_order=categories,
            dot_color_df=color, dot_size_df=size, standard_scale=None,
            vcenter=0 if kind == "programs" else None)
        plot.style(cmap="Reds" if kind == "genes" else "coolwarm")
        plot.legend(colorbar_title="Normalized expression" if kind == "genes" else "Mean CLR value",
            size_title="Fraction of cells\nexpressing (%)" if kind == "genes" else "Fraction of cells\nwith CLR>0 (%)")
        plots.append(plot)
    return (*plots, tables[f"{panel}_genes"], tables[f"{panel}_programs"],
            {**manifest["metadata"][panel], "plot_input_mode": "precomputed dot matrices", "source_sha256": manifest["source"]["sha256"]})


def prepare_geo_plot_inputs(dataset, root, paths, *, mode="precomputed", accession="", filename=None,
                            local_file=None, expected_sha256=None):
    """Choose explicitly; never fetch/recompute as an implicit fallback."""
    from .public_data import value
    if dataset not in DATASETS:
        raise ValueError(dataset)
    key = "hippo" if dataset == "hippo_timecourse" else "mixed"
    published = value(paths, "precomputed." + key)
    if mode == "precomputed":
        if accession or local_file:
            raise ValueError("Select mode='geo' or 'local_file' to use the supplied GEO accession/processed file.")
        _, manifest = load_plot_bundle(published, dataset, tables=False)
        return {"directory": str(published), "mode": mode, "source_sha256": manifest["source"]["sha256"], "requires_geo_download": False}
    if mode not in {"geo", "local_file"}:
        raise ValueError("GEO experiment mode must be 'precomputed', 'geo', or 'local_file'.")
    if mode == "geo":
        if local_file:
            raise ValueError("Use only a GEO accession in geo mode, or select local_file mode.")
        from .geo_download import GeoClient
        source = GeoClient().download(accession, value(paths, "public_data.geo_cache"), dataset=dataset,
                                      filename=filename, expected_sha256=expected_sha256)
    else:
        if not local_file or accession:
            raise ValueError("local_file mode requires one processed file and no GEO accession.")
        source = Path(local_file).expanduser()
        source = source.resolve() if source.is_absolute() else (Path(root) / source).resolve()
        if not source.is_file():
            raise FileNotFoundError(f"Local submitted object not found: {source}")
    source_hash = digest_file(source, "sha256")
    if expected_sha256 and source_hash != str(expected_sha256).lower():
        raise ValueError("Processed source SHA-256 differs from the supplied expected value.")
    import hashlib
    import importlib.metadata
    module_root = Path(__file__).parent
    implementation_record = {name: digest_file(module_root / name, "sha256") for name in
        ("geo_plot_data.py", "geo_submission.py", "figure3_panels_legacy.py", "figure4_panels_legacy.py", "hippo_projection.py", "figure_s02_parse_mix_plots.py", "figure_s02_source_plotnine.py")}
    implementation_record["versions"] = {name: importlib.metadata.version(name) for name in
        ("numpy", "pandas", "anndata", "mudata", "scanpy", "plotnine")}
    implementation = hashlib.sha256(json.dumps(implementation_record, sort_keys=True).encode()).hexdigest()[:12]
    destination = value(paths, "public_data.derived_cache") / "geo_plots" / dataset / f"{source_hash[:16]}-{implementation}"
    if destination.exists():
        _, manifest = load_plot_bundle(destination, dataset, tables=False)
        if manifest["source"]["sha256"] != source_hash:
            raise ValueError("Cached plot bundle source differs from the requested object.")
    else:
        build_plot_bundle(source, destination, dataset)
    comparison = {"status": "NO_PUBLISHED_REFERENCE_BUNDLE"}
    if (published / "manifest.json").is_file():
        import pandas as pd
        accepted, _ = load_plot_bundle(published, dataset)
        regenerated, _ = load_plot_bundle(destination, dataset)
        differences = []
        for name in accepted:
            try:
                pd.testing.assert_frame_equal(accepted[name], regenerated[name], check_dtype=False,
                    check_index_type=False, check_categorical=False, rtol=1e-6, atol=1e-7)
            except AssertionError:
                differences.append(name)
        comparison = {"status": "MATCH" if not differences else "DIFFERENT", "different_tables": differences}
        if differences:
            warnings.warn("Regenerated GEO data differ from published plotting tables: " + ", ".join(differences))
    (destination / "comparison_to_published.json").write_text(json.dumps(comparison, indent=2) + "\n")
    return {"directory": str(destination), "mode": mode, "accession": accession or None,
        "source_sha256": source_hash, "published_comparison": comparison, "published_data_overwritten": False}
