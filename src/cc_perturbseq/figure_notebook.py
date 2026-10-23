from __future__ import annotations

from pathlib import Path
import hashlib
import json
import os
import platform
import subprocess
import sys
from typing import Any

import matplotlib as mpl
import pandas as pd
import yaml


def find_repo_root(start: Path | None = None) -> Path:
    start = (start or Path.cwd()).resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "docs" / "REPRODUCTION_PLAN.md").exists():
            return candidate
    raise FileNotFoundError(
        "Could not locate the repository root. Run the notebook from inside "
        "cc-perturbseq-manuscript-analysis."
    )


def load_paths(root: Path) -> tuple[Path, dict[str, Any]]:
    config = Path(
        os.environ.get(
            "CCPERTURBSEQ_PATHS",
            root / "config" / "paths.local.yaml",
        )
    ).expanduser()
    if not config.exists():
        raise FileNotFoundError(
            f"Missing {config}. Copy config/paths.example.yaml to "
            "config/paths.local.yaml and edit the Sherlock paths."
        )
    return config, (yaml.safe_load(config.read_text()) or {})


def get_path(paths: dict[str, Any], dotted_key: str, *, required: bool = True) -> Path | None:
    node: Any = paths
    for part in dotted_key.split("."):
        if not isinstance(node, dict) or part not in node:
            if required:
                raise KeyError(f"Missing configuration key: {dotted_key}")
            return None
        node = node[part]
    if node in (None, ""):
        if required:
            raise ValueError(f"Empty configuration value: {dotted_key}")
        return None
    return Path(os.path.expandvars(str(node))).expanduser()


def git_commit(root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip()
    except Exception:
        return "unknown"


def audit_inputs(paths: dict[str, Any], keys: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for key in keys:
        try:
            path = get_path(paths, key)
            exists = path.exists()
            rows.append({
                "key": key,
                "path": str(path),
                "exists": exists,
                "is_file": path.is_file() if exists else False,
                "size_bytes": path.stat().st_size if exists and path.is_file() else None,
            })
        except Exception as exc:
            rows.append({
                "key": key,
                "path": None,
                "exists": False,
                "is_file": False,
                "size_bytes": None,
                "error": str(exc),
            })
    return pd.DataFrame(rows)


def load_figure_manifest(figure_id: str, root: Path | None = None) -> dict[str, Any]:
    root = root or find_repo_root()
    manifest_file = root / "config" / "figure_notebooks.yaml"
    if not manifest_file.exists():
        raise FileNotFoundError(manifest_file)
    payload = yaml.safe_load(manifest_file.read_text()) or {}
    try:
        spec = payload["figures"][figure_id]
    except KeyError as exc:
        raise KeyError(f"Figure {figure_id!r} is not defined in {manifest_file}") from exc
    return {"figure_id": figure_id, **spec}


def initialize_figure(manifest: dict[str, Any]) -> dict[str, Any]:
    root = find_repo_root()
    config, paths = load_paths(root)
    output_dir = get_path(paths, "outputs.root") / manifest["figure_id"]
    output_dir.mkdir(parents=True, exist_ok=True)

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.bbox": "tight",
    })

    required = sorted({
        key
        for panel in manifest["panels"]
        for key in panel["required_inputs"]
    })
    panel_table = pd.DataFrame(manifest["panels"])
    panel_table.attrs["figure_id"] = manifest["figure_id"]
    return {
        "root": root,
        "config": config,
        "paths": paths,
        "output_dir": output_dir,
        "git_commit": git_commit(root),
        "panels": panel_table,
        "input_audit": audit_inputs(paths, required),
    }


def panel_directory(context: dict[str, Any], panel: str) -> Path:
    safe = "".join(c if c.isalnum() else "_" for c in panel).strip("_") or "all"
    path = context["output_dir"] / f"panel_{safe}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def code_panel_status(context: dict[str, Any], panel: str, description: str) -> dict[str, Any]:
    return {
        "panel": panel,
        "description": description,
        "status": "CANDIDATES_FOUND",
        "canonical_code_ported": False,
        "output_dir": str(panel_directory(context, panel)),
    }


def source_asset_status(
    context: dict[str, Any],
    panel: str,
    asset_key: str,
) -> dict[str, Any]:
    path = get_path(context["paths"], asset_key, required=False)
    return {
        "panel": panel,
        "status": "SOURCE_ASSET_REQUIRED",
        "asset_key": asset_key,
        "asset_path": str(path) if path else None,
        "exists": bool(path and path.exists()),
    }


def sha256_file(path: Path, *, max_bytes: int = 2_000_000_000) -> str | None:
    if path.stat().st_size > max_bytes:
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_receipt(
    context: dict[str, Any],
    panel: str,
    *,
    status: str,
    inputs: dict[str, Any],
    outputs: dict[str, Any],
    parameters: dict[str, Any],
) -> Path:
    payload = {
        "figure": context["panels"].attrs.get("figure_id"),
        "panel": panel,
        "status": status,
        "git_commit": context["git_commit"],
        "python": sys.version,
        "platform": platform.platform(),
        "inputs": inputs,
        "outputs": outputs,
        "parameters": parameters,
    }
    path = panel_directory(context, panel) / "run_receipt.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return path
