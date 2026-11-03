"""Public figure inputs: IGVF deposits, versioned companion data, and GEO H5MU.

Nothing in this module searches an institution filesystem or runs an archived
notebook. Numerical representations are preserved; cNMF is not refitted.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import warnings

from .igvf_download import IgvfClient, digest_file

DAYS = ("d0", "d1", "d2", "d3")
_WORKBOOK = "data/figure_02/cNMF_50_0_2_032526_6co.xlsx"
_FIGURES = ("Figure_01", "Figure_02", "Figure_03", "Figure_04", "Figure_S02", "Figure_S04", "Figure_S06")
_CLIENT = None

COMPANION_KEYS = {
    "Figure_01": ["screen.genomewide_utest_pos_vs_neg", "validation.compatibility", "cell_order.expression", "cell_order.cnmf", "ccperturbseq.figure1d_program_metadata", "ccperturbseq.figure1d_program_order"],
    "Figure_02": ["ccperturbseq.cnmf_usage", "ccperturbseq.target_metadata", "cell_order.expression"],
    "Figure_03": ["ccperturbseq.cnmf_usage", "ccperturbseq.target_metadata", "ccperturbseq.cnmf_spectra", "cell_order.expression", "multiome.enhancer_enrichment_results", "multiome.cluster_motif_map", "multiome.modisco_compiled_h5"],
    "Figure_04": ["ccperturbseq.cnmf_usage", "ccperturbseq.cnmf_spectra", "cell_order.expression", "multiome.enhancer_enrichment_results", "multiome.cluster_motif_map", "multiome.modisco_compiled_h5"],
    "Figure_S04": ["s4.*"],
    "Figure_S06": ["ccperturbseq.cnmf_spectra", "ccperturbseq.program_annotations", "multiome.enhancer_enrichment_results", "multiome.cluster_motif_map", "multiome.modisco_compiled_h5"],
}


def configure_igvf(*, allow_prompt=True):
    """Configure optional getpass authentication. Never returns credentials."""
    global _CLIENT
    _CLIENT = IgvfClient(allow_prompt=allow_prompt)


def repo_root():
    return Path(__file__).resolve().parents[2]


def value(paths, key):
    node = paths
    for part in key.split("."):
        node = node[part]
    return Path(node)


def load_public_paths(root):
    """Resolve paths relative to the repository, never the notebook directory."""
    import yaml
    root = Path(root).resolve()
    config = root / "config/paths.example.yaml"
    paths = yaml.safe_load(config.read_text())
    explicit = os.environ.get("CCPERTURBSEQ_PUBLIC_PATHS")
    if explicit:
        config = Path(explicit).expanduser().resolve()
        overrides = yaml.safe_load(config.read_text()) or {}
        def merge(a, b):
            for k, v in b.items():
                if isinstance(v, dict) and isinstance(a.get(k), dict):
                    merge(a[k], v)
                else:
                    a[k] = v
        merge(paths, overrides)
    elif (root / "config/paths.local.yaml").exists() or os.environ.get("CCPERTURBSEQ_PATHS"):
        warnings.warn("Legacy paths.local.yaml/CCPERTURBSEQ_PATHS is not used by public figure notebooks. Use CCPERTURBSEQ_PUBLIC_PATHS for deliberate portable overrides.")
    def resolve(node):
        for k, v in node.items():
            if isinstance(v, dict):
                resolve(v)
            elif isinstance(v, str):
                p = Path(os.path.expandvars(v)).expanduser()
                p = (p if p.is_absolute() else root / p).resolve()
                if str(p).startswith(("/oak/", "/scratch/users/", "/home/users/")) and not p.is_relative_to(root):
                    raise ValueError(f"External institutional path is not allowed in public configuration: {k}. Use downloaded or bundled data.")
                node[k] = str(p)
    resolve(paths)
    # Generated files must stay within this checkout's output/ directory even
    # when a user supplies an override or has SCRATCH/TMPDIR exported.
    output_root = root / "output"
    resolved_output = output_root.resolve()
    if not resolved_output.is_relative_to(root):
        raise ValueError("Repository output/ must not be a symlink outside the checkout.")
    if Path(paths["outputs"]["root"]) != resolved_output:
        raise ValueError("outputs.root must be the repository's output/ directory. Remove the output-root override.")
    generated_keys = (
        "public_data.derived_cache", "screen.validation_guide_counts",
        "ccperturbseq.integrated_anndata", "ccperturbseq.cnmf_h5mu",
        "hippo_timecourse.query_h5ad", "hippo_timecourse.reprojected_cnmf_h5ad",
    )
    for key in generated_keys:
        try:
            path = value(paths, key)
        except KeyError:
            continue
        if not path.is_relative_to(resolved_output):
            raise ValueError(f"Generated file/cache {key} must be beneath repository output/.")
    return config, paths


def _manifest(root):
    return json.loads((Path(root) / "config/public_data_manifest.json").read_text())


def _nonempty(path):
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def _check_lfs(path):
    with Path(path).open("rb") as stream:
        prefix = stream.read(150)
    if prefix.startswith(b"version https://git-lfs.github.com/spec/"):
        raise FileNotFoundError(f"Git LFS pointer, not data: {path}. Run git lfs pull.")


def require_bundle(root, paths, keys):
    """Verify every requested companion file against its publication manifest."""
    root = Path(root).resolve()
    manifest_path = value(paths, "public_data.bundle_manifest")
    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"Companion data has not been published/materialized: {manifest_path}. "
            "The maintainer must run scripts/collect_publication_data.py --apply and publish data/reproduction/. "
            "There is no fallback to original filesystem paths."
        )
    records = json.loads(manifest_path.read_text())["files"]
    lookup = {record["key"]: record for record in records}
    selected = [k for k in lookup if k.startswith("s4.")] if keys == ["s4.*"] else list(keys)
    if not selected:
        raise FileNotFoundError("No S4 numerical files in the companion manifest.")
    resolved = {}
    for key in selected:
        if key not in lookup:
            raise FileNotFoundError(f"Companion manifest lacks required input {key}.")
        record = lookup[key]
        relative = Path(record["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Companion manifest must contain repository-relative paths only.")
        path = (root / relative).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Companion path escapes repository.")
        if not _nonempty(path):
            raise FileNotFoundError(f"Missing companion input: {path}. Retrieve the published data bundle; no institutional fallback is used.")
        _check_lfs(path)
        if path.stat().st_size != record["size_bytes"] or digest_file(path, "sha256") != record["sha256"]:
            raise ValueError(f"Companion checksum mismatch: {path}")
        resolved[key] = path
    return resolved


def _stamp(inputs):
    return hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()


def _is_current(output, signature):
    output = Path(output)
    receipt = output.with_name(output.name + ".inputs.json")
    if not _nonempty(output) or not receipt.is_file():
        return False
    try:
        record = json.loads(receipt.read_text())
        return record["signature"] == signature and record["size_bytes"] == output.stat().st_size and record["mtime_ns"] == output.stat().st_mtime_ns
    except (KeyError, OSError, ValueError):
        return False


def _atomic_write(output, signature, writer):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(dir=output.parent, prefix=output.stem + ".", suffix=output.suffix)
    os.close(fd)
    temporary = Path(name)
    try:
        writer(temporary)
        if not _nonempty(temporary):
            raise ValueError("Derived input writer produced no data.")
        os.replace(temporary, output)
        output.with_name(output.name + ".inputs.json").write_text(json.dumps({
            "signature": signature, "size_bytes": output.stat().st_size,
            "mtime_ns": output.stat().st_mtime_ns, "sha256": digest_file(output, "sha256"),
        }, indent=2) + "\n")
    finally:
        temporary.unlink(missing_ok=True)


def _row_order(path):
    import pandas as pd
    order = pd.read_csv(path, sep="\t", dtype=str)["cell_id"].tolist()
    if len(order) != len(set(order)):
        raise ValueError(f"Duplicate cells in ordering file: {path}")
    return order


def _read_cnmf(path):
    import h5py
    from anndata.io import read_elem
    # Read only the program modality; the much larger RNA modality is unused
    # by Figure 1d. No matrix rescaling, projection or fitting occurs here.
    with h5py.File(path, "r") as handle:
        return read_elem(handle["mod"]["cNMF"])


def assemble_day_objects(day_paths, cell_order_path, output, *, kind="expression", signature=None):
    """Concatenate deposited day partitions, then restore the original row order."""
    if kind not in {"expression", "cnmf"}:
        raise ValueError(kind)
    import anndata as ad
    import numpy as np
    signature = signature or _stamp({"version": 1, "kind": kind,
        "inputs": [(str(p), digest_file(p)) for p in day_paths],
        "order": digest_file(cell_order_path, "sha256")})
    if _is_current(output, signature):
        return Path(output)
    day_paths = list(day_paths)
    if len(day_paths) != 4:
        raise ValueError("Exactly four day deposits are required.")
    objects = []
    reference_genes = reference_guides = None
    for day, path in zip(DAYS, day_paths):
        obj = ad.read_h5ad(path) if kind == "expression" else _read_cnmf(path)
        if not obj.obs_names.is_unique or not obj.var_names.is_unique:
            raise ValueError(f"Duplicate cell/feature IDs in {day} {kind} deposit.")
        genes = obj.var_names.astype(str).tolist()
        guides = [str(g) for g in obj.uns.get("guide_names", [])]
        if reference_genes is None:
            reference_genes, reference_guides = genes, guides
        elif genes != reference_genes or guides != reference_guides:
            raise ValueError(f"Feature or guide-column order differs between day deposits: {day}. Do not concatenate without explicit alignment.")
        if "day" in obj.obs:
            actual = obj.obs["day"].astype(str).str.lower().str.replace("day", "d").str.replace(" ", "")
            if not actual.eq(day).all():
                raise ValueError(f"Day annotation mismatch in {day} deposit.")
        objects.append(obj)
    if len(objects) != 4:
        raise ValueError("All four deposited day objects are required.")
    order = _row_order(cell_order_path)
    all_ids = [str(x) for obj in objects for x in obj.obs_names]
    if len(set(all_ids)) != len(all_ids) or set(order) != set(all_ids):
        raise ValueError("Deposited cells do not exactly match the original ordering manifest.")
    merged = ad.concat(objects, axis=0, join="inner", merge="same", uns_merge="same", index_unique=None)
    merged = merged[order, :].copy()
    if list(merged.var_names.astype(str)) != reference_genes:
        raise ValueError("Feature order changed during concatenation.")
    if reference_guides:
        merged.uns["guide_names"] = np.asarray(reference_guides, dtype=str)
        if merged.obsm["guide_assignment"].shape[1] != len(reference_guides):
            raise ValueError("Guide assignments do not match their explicit identifiers.")
    if kind == "expression":
        if "X_umap" not in merged.obsm:
            raise ValueError("Deposited expression objects lack the required existing UMAP.")
        _atomic_write(output, signature, lambda p: merged.write_h5ad(p, compression="gzip"))
    else:
        import mudata as md
        bundle = md.MuData({"cNMF": merged})
        _atomic_write(output, signature, lambda p: bundle.write_h5mu(p))
    return Path(output)


def adapt_validation_counts(portal_path, compatibility_path, output):
    """Use deposited counts, keeping original guide labels/order and QC totals."""
    import numpy as np
    import pandas as pd
    signature = _stamp({"v": 1, "portal": digest_file(portal_path), "compatibility": digest_file(compatibility_path, "sha256")})
    if _is_current(output, signature):
        return Path(output)
    portal = pd.read_csv(portal_path, sep="\t")
    bridge = pd.read_csv(compatibility_path, sep="\t")
    if portal.guide_id.duplicated().any() or bridge.portal_guide_id.duplicated().any():
        raise ValueError("Validation guide identifiers must map one-to-one.")
    if set(portal.guide_id) != set(bridge.portal_guide_id):
        raise ValueError("Validation guide bridge does not cover the deposited guide set.")
    portal = portal.set_index("guide_id").loc[bridge.portal_guide_id]
    table = bridge[["sgRNA", "Gene", "D3TOTAL1", "D3TOTAL2"]].copy()
    for source, destination in {"D01": "D0R1", "D02": "D0R2", "D3neg1": "D3NR1", "D3neg2": "D3NR2", "D3pos1": "D3PR1", "D3pos2": "D3PR2"}.items():
        values = pd.to_numeric(portal[source], errors="raise").to_numpy()
        if not (np.isfinite(values) & (values >= 0)).all():
            raise ValueError("Invalid deposited validation counts.")
        table[destination] = values
    _atomic_write(output, signature, lambda p: table.to_csv(p, sep="\t", index=False))
    return Path(output)


def require_geo(path, dataset):
    path = Path(path)
    if not _nonempty(path):
        raise FileNotFoundError(
            f"Missing processed GEO object for {dataset}: {path}. "
            "Place the submitted H5MU at this portable path (or configure geo.* in CCPERTURBSEQ_PUBLIC_PATHS). "
            "Use precomputed plotting tables or explicitly select geo mode with a public accession. No source-data fallback is used."
        )
    _check_lfs(path)
    return path


def extract_hippo_geo(path, query_output, projection_output):
    import mudata as md
    from .geo_submission import SCHEMA
    path = require_geo(path, "Hippo time course")
    signature = _stamp({"v": 1, "geo_sha256": digest_file(path, "sha256")})
    if all(_is_current(p, signature) for p in (query_output, projection_output)):
        return
    bundle = md.read_h5mu(path)
    if bundle.uns.get("geo_schema") != SCHEMA or bundle.uns.get("dataset") != "hippo_timecourse":
        raise ValueError("Not the submitted integrated Hippo H5MU schema.")
    if not bundle.mod["rna"].obs_names.equals(bundle.mod["cnmf"].obs_names):
        raise ValueError("Hippo RNA/program cell order differs.")
    for name, output in (("rna", query_output), ("cnmf", projection_output)):
        _atomic_write(output, signature, lambda p, name=name: bundle.mod[name].write_h5ad(p, compression="gzip"))


def load_mixed_geo(path):
    """Use frozen lineage labels; recalculate guide metrics, never recluster."""
    import mudata as md
    import numpy as np
    from .geo_submission import SCHEMA, METRIC_COLUMNS, _metrics_from_sparse_guides
    from .figure_s02_parse_mix_data import ParseMixResult
    path = require_geo(path, "mixed hiPSC–K562")
    bundle = md.read_h5mu(path)
    if bundle.uns.get("geo_schema") != SCHEMA or bundle.uns.get("dataset") != "ipsc_k562_mix":
        raise ValueError("Not the submitted integrated mixed-cell H5MU schema.")
    if not bundle.mod["rna"].obs_names.equals(bundle.mod["crispr"].obs_names):
        raise ValueError("Mixed-cell RNA/guide cell order differs.")
    metrics = _metrics_from_sparse_guides(bundle)
    for key in METRIC_COLUMNS:
        np.testing.assert_allclose(metrics[key].to_numpy(float), bundle.mod["rna"].obs[key].to_numpy(float), rtol=1e-6, atol=1e-7, equal_nan=True, err_msg=key)
    return ParseMixResult(metrics, bundle.uns["cluster_annotation"].copy(),
        bundle.uns["saturation"].copy(), bundle.uns["manuscript_median_diagnostics"].copy(),
        {"processed_file": str(path), "sha256": digest_file(path, "sha256"), "lineages": "frozen in submitted object", "original_source_files_required": False})


def prepare_figure_inputs(figure_id, root, paths, *, geo_options=None):
    """Materialize only the inputs for the requested retained figure."""
    if figure_id not in _FIGURES:
        raise ValueError(f"No public input specification for {figure_id}")
    root = Path(root).resolve()
    bundled = require_bundle(root, paths, COMPANION_KEYS[figure_id]) if figure_id in COMPANION_KEYS else {}
    # Companion overrides may not bypass verification of the actual input.
    for key, verified in bundled.items():
        if key.startswith(("screen.", "ccperturbseq.", "multiome.")):
            if value(paths, key).resolve() != verified:
                raise ValueError(f"Configured {key} differs from its verified companion path.")
    if figure_id in {"Figure_02", "Figure_03", "Figure_04", "Figure_S06"}:
        if not _nonempty(root / _WORKBOOK):
            raise FileNotFoundError(f"Missing repository workbook: {root / _WORKBOOK}")
        _check_lfs(root / _WORKBOOK)
    geo_plotting = None
    if figure_id in {"Figure_03", "Figure_04", "Figure_S02"}:
        from .geo_plot_data import prepare_geo_plot_inputs
        dataset = "ipsc_k562_mix" if figure_id == "Figure_S02" else "hippo_timecourse"
        geo_plotting = prepare_geo_plot_inputs(dataset, root, paths, **(geo_options or {}))
    manifest = _manifest(root)
    def fetch_group(name):
        client = _CLIENT or IgvfClient(allow_prompt=True)
        return [client.download(manifest["files"][acc], value(paths, "public_data.igvf_cache")) for acc in manifest["groups"][name]]
    if figure_id in {"Figure_01", "Figure_02", "Figure_03", "Figure_04"}:
        files = fetch_group("expression")
        assemble_day_objects(files, bundled["cell_order.expression"], value(paths, "ccperturbseq.integrated_anndata"), signature=_stamp({"v": 1, "inputs": [manifest["files"][a]["md5sum"] for a in manifest["groups"]["expression"]], "order": digest_file(bundled["cell_order.expression"], "sha256")}))
    if figure_id == "Figure_01":
        files = fetch_group("cnmf")
        assemble_day_objects(files, bundled["cell_order.cnmf"], value(paths, "ccperturbseq.cnmf_h5mu"), kind="cnmf", signature=_stamp({"v": 1, "inputs": [manifest["files"][a]["md5sum"] for a in manifest["groups"]["cnmf"]], "order": digest_file(bundled["cell_order.cnmf"], "sha256")}))
        mageck, validation = fetch_group("screen")
        paths["screen"]["genomewide_mageck_pos_vs_neg"] = str(mageck)
        adapt_validation_counts(validation, bundled["validation.compatibility"], value(paths, "screen.validation_guide_counts"))
    return {"figure": figure_id, "backend": "IGVF and published plotting companions",
            "geo_plotting": geo_plotting, "institutional_source_files_required": False}
