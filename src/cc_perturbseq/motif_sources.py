"""Canonical motif-enrichment inputs for every manuscript figure.

The enrichment table is intentionally locked to the file selected by the
manuscript author.  A local paths override must not silently substitute an
older ``_v0`` table or the previously tracked ``*_min10.tsv`` derivative.
"""
from __future__ import annotations

import os
from pathlib import Path

CANONICAL_ENHANCER_ENRICHMENT_PATH = Path(
    "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/"
    "analysis/2.tf_enrichment_in_programs/enhancer_enrichment_results.tsv"
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
            "Motif enrichment is locked to the manuscript-authorized table. "
            f"Configured: {candidate}. Required: {canonical}. "
            "Remove stale paths.local.yaml overrides, including "
            "enhancer_enrichment_results_min10.tsv and multiome_ipsc_ec_v0."
        )
    if require_exists:
        if not canonical.is_file():
            raise FileNotFoundError(
                "Canonical motif-enrichment table is missing: "
                f"{canonical}"
            )
        if not os.access(canonical, os.R_OK):
            raise PermissionError(
                "Canonical motif-enrichment table is not readable: "
                f"{canonical}"
            )
    return canonical


def canonical_motif_logo_inputs(
    *,
    require_exists: bool = True,
) -> tuple[Path, Path]:
    """Return the current non-v0 cluster map and compiled MoDISco H5."""
    paths = (
        _normalized(CANONICAL_CLUSTER_MOTIF_MAP_PATH),
        _normalized(CANONICAL_MODISCO_COMPILED_H5_PATH),
    )
    if require_exists:
        for path in paths:
            if not path.is_file():
                raise FileNotFoundError(f"Canonical motif-logo input is missing: {path}")
            if not os.access(path, os.R_OK):
                raise PermissionError(f"Canonical motif-logo input is not readable: {path}")
    return paths


__all__ = [
    "CANONICAL_ENHANCER_ENRICHMENT_PATH",
    "CANONICAL_CLUSTER_MOTIF_MAP_PATH",
    "CANONICAL_MODISCO_COMPILED_H5_PATH",
    "canonical_enhancer_enrichment_path",
    "canonical_motif_logo_inputs",
]
