"""Manuscript motif inputs distributed with the versioned companion data."""
from pathlib import Path
from .public_data import load_public_paths, require_bundle, repo_root

_ROOT = repo_root()
CANONICAL_ENHANCER_ENRICHMENT_PATH = _ROOT / "data/reproduction/motifs/enhancer_enrichment_results_min10.tsv"
CANONICAL_CLUSTER_MOTIF_MAP_PATH = _ROOT / "data/reproduction/motifs/cluster_to_tf_motif_id_mapping.tsv"
CANONICAL_MODISCO_COMPILED_H5_PATH = _ROOT / "data/reproduction/motifs/modisco_compiled.h5"


def _selected(key, configured=None, *, require_exists=True):
    root = repo_root()
    _, paths = load_public_paths(root)
    expected = Path(paths["multiome"][key]).resolve()
    candidate = expected if configured is None else Path(configured).expanduser().resolve()
    if candidate != expected:
        raise ValueError(f"Use the manuscript's bundled motif file for {key}; a different table/compendium was supplied.")
    if require_exists:
        verified = require_bundle(root, paths, ["multiome." + key])["multiome." + key]
        if verified != candidate:
            raise ValueError("Configured motif file differs from the checksum-pinned companion manifest.")
    return candidate


def canonical_enhancer_enrichment_path(configured=None, *, require_exists=True):
    return _selected("enhancer_enrichment_results", configured, require_exists=require_exists)


def canonical_motif_logo_inputs(configured_cluster_map=None, configured_modisco_h5=None, *, require_exists=True):
    return (_selected("cluster_motif_map", configured_cluster_map, require_exists=require_exists),
            _selected("modisco_compiled_h5", configured_modisco_h5, require_exists=require_exists))


def canonical_motif_input_manifest(*, require_exists=True):
    return {key: _selected(key, require_exists=require_exists) for key in
            ("enhancer_enrichment_results", "cluster_motif_map", "modisco_compiled_h5")}
