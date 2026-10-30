"""Filesystem preflight for the mixed hiPSC–K562 Parse-v3 experiment.

The source experiment stores expression inputs inside a directory named
``DGE_filtered``.  That directory is a required directory, while the matrix,
metadata, guide, and saturation inputs beneath the experiment are required
non-empty files.  Keeping these two path types explicit prevents a valid
``DGE_filtered`` directory from being rejected by ``Path.is_file()``.
"""
from __future__ import annotations

from pathlib import Path

from . import figure_s02_parse_mix_data as _data


def _file_issue(path: Path) -> str | None:
    if not path.exists():
        return f"missing file: {path}"
    if not path.is_file():
        return f"expected file but found another path type: {path}"
    if path.stat().st_size == 0:
        return f"empty file: {path}"
    return None


def _directory_issue(path: Path) -> str | None:
    if not path.exists():
        return f"missing directory: {path}"
    if not path.is_dir():
        return f"expected directory but found another path type: {path}"
    return None


def preflight_parse_mix(root: Path) -> dict[str, Path]:
    """Validate the exact Split-pipe paths used by the source notebook."""
    paths = _data._required_paths(root)
    issues: list[str] = []

    for directory_key in ("root", "expression_dir"):
        issue = _directory_issue(paths[directory_key])
        if issue:
            issues.append(issue)

    for file_key in (
        "count_matrix",
        "genes",
        "cell_metadata",
        "guide_counts",
        "guide_reference",
    ):
        issue = _file_issue(paths[file_key])
        if issue:
            issues.append(issue)

    for primer in _data.PRIMER_ORDER:
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
        issue = _file_issue(saturation)
        if issue:
            issues.append(issue)

    if issues:
        raise FileNotFoundError(
            "Missing, empty, or mistyped Parse-v3 mixed hiPSC–K562 inputs:\n  "
            + "\n  ".join(issues)
        )
    return paths


__all__ = ["preflight_parse_mix"]
