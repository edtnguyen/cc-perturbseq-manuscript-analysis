"""Canonical motif inputs shared by Figure 3d, Figure 4e, and Figure S6.

The author-selected enrichment result is the active min-10-carrier table.  The
logo inputs are the active non-archive motif map and the standard compiled
MoDISco compendium used by the recovered source notebooks.  Archived
sensitivity tables and the separate ``modisco_compiled_profile`` compendium
must not be substituted silently.
"""
from __future__ import annotations

import os
from pathlib import Path

CANONICAL_ENHANCER_ENRICHMENT_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/"
    "analysis/2.tf_enrichment_in_programs/"
    "enhancer_enrichment_results_min10.tsv"
)
CANONICAL_CLUSTER_MOTIF_MAP_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/"
    "analysis/0.0.motif_id_annotation/tables/"
    "cluster_to_tf_motif_id_mapping.tsv"
)
CANONICAL_MODISCO_COMPILED_H5_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/"
    "results/compendium/modisco_compiled/modisco_compiled.h5"
)


def _normalized(path: str | Path) -> Path:
    return Path(path).expanduser().resolve(strict=False)


def _require_readable(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Canonical {label} is missing: {path}")
    if not os.access(path, os.R_OK):
        raise PermissionError(f"Canonical {label} is not readable: {path}")


def canonical_enhancer_enrichment_path(
    configured: str | Path | None = None,
    *,
    require_exists: bool = True,
) -> Path:
    """Return the one permitted motif-enrichment table.

    ``configured`` is accepted only to audit the merged path configuration.
    Any different path is rejected rather than used as a fallback.
    """
    canonical = _normalized(CANONICAL_ENHANCER_ENRICHMENT_PATH)
    candidate = canonical if configured is None else _normalized(configured)
    if candidate != canonical:
        raise ValueError(
            "Motif enrichment is locked to the active min10 table selected "
            "for the manuscript. "
            f"Configured: {candidate}. Required: {canonical}. "
            "Remove stale paths.local.yaml overrides, including archived "
            "enhancer_enrichment_results.tsv, union-links, pos-only, "
            "zero-filled, carriers-min5, min1, and min20 variants."
        )
    if require_exists:
        _require_readable(canonical, "motif-enrichment table")
    return canonical


def canonical_motif_logo_inputs(
    configured_cluster_map: str | Path | None = None,
    configured_modisco_h5: str | Path | None = None,
    *,
    require_exists: bool = True,
) -> tuple[Path, Path]:
    """Return and audit the motif-name map and compiled MoDISco logo source.

    The standard compiled H5 is retained because it is the file referenced by
    the recovered Figure 3/4/S6 source notebooks.  The later
    ``modisco_compiled_profile`` file is a separate profile compendium and is
    not treated as an interchangeable rename.
    """
    canonical_map = _normalized(CANONICAL_CLUSTER_MOTIF_MAP_PATH)
    canonical_h5 = _normalized(CANONICAL_MODISCO_COMPILED_H5_PATH)
    candidate_map = (
        canonical_map
        if configured_cluster_map is None
        else _normalized(configured_cluster_map)
    )
    candidate_h5 = (
        canonical_h5
        if configured_modisco_h5 is None
        else _normalized(configured_modisco_h5)
    )
    mismatches = []
    if candidate_map != canonical_map:
        mismatches.append(
            f"cluster map configured={candidate_map}, required={canonical_map}"
        )
    if candidate_h5 != canonical_h5:
        mismatches.append(
            f"MoDISco H5 configured={candidate_h5}, required={canonical_h5}"
        )
    if mismatches:
        raise ValueError(
            "Motif-logo inputs are locked to the active non-archive cluster "
            "map and the standard compiled MoDISco H5 used by the recovered "
            "source notebooks. "
            + "; ".join(mismatches)
            + ". Do not substitute archive_v1 mappings or "
            "modisco_compiled_profile/modisco_compiled.h5."
        )
    if require_exists:
        _require_readable(canonical_map, "motif-cluster map")
        _require_readable(canonical_h5, "compiled MoDISco H5")
    return canonical_map, canonical_h5


def canonical_motif_input_manifest(
    *,
    require_exists: bool = True,
) -> dict[str, Path]:
    """Return all motif-specific files used across the retained figures."""
    enrichment = canonical_enhancer_enrichment_path(
        require_exists=require_exists
    )
    cluster_map, modisco_h5 = canonical_motif_logo_inputs(
        require_exists=require_exists
    )
    return {
        "enhancer_enrichment_results": enrichment,
        "cluster_motif_map": cluster_map,
        "modisco_compiled_h5": modisco_h5,
    }


__all__ = [
    "CANONICAL_ENHANCER_ENRICHMENT_PATH",
    "CANONICAL_CLUSTER_MOTIF_MAP_PATH",
    "CANONICAL_MODISCO_COMPILED_H5_PATH",
    "canonical_enhancer_enrichment_path",
    "canonical_motif_input_manifest",
    "canonical_motif_logo_inputs",
]
