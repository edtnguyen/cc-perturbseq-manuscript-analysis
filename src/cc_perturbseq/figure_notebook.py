from __future__ import annotations

from pathlib import Path
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import warnings
from typing import Any

import matplotlib as mpl
import pandas as pd
import yaml


_PLACEHOLDER_PREFIXES = (
    "REPLACE_WITH_",
    "EMBEDDED_IN_",
)
_UNRESOLVED_ENV_PATTERN = re.compile(r"\$(?:\{[^}]+\}|[A-Za-z_][A-Za-z0-9_]*)")


def find_repo_root(start: Path | None = None) -> Path:
    start = (start or Path.cwd()).resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "docs" / "REPRODUCTION_PLAN.md").exists():
            return candidate
    raise FileNotFoundError(
        "Could not locate the repository root. Run the notebook from inside "
        "cc-perturbseq-manuscript-analysis."
    )


def _is_placeholder(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    text = value.strip()
    if not text:
        return True
    upper = text.upper()
    return upper.startswith(_PLACEHOLDER_PREFIXES)


def _merge_config(defaults: Any, overrides: Any) -> Any:
    """Recursively merge local overrides over the tracked Sherlock defaults.

    Empty or placeholder override values do not replace a concrete tracked value.
    This lets an older paths.local.yaml inherit newly resolved paths from
    paths.example.yaml without manual recopying.
    """
    if isinstance(defaults, dict) and isinstance(overrides, dict):
        merged = dict(defaults)
        for key, value in overrides.items():
            if key in merged:
                merged[key] = _merge_config(merged[key], value)
            else:
                merged[key] = value
        return merged

    if overrides is None or _is_placeholder(overrides):
        return defaults
    return overrides


def _resolve_output_root(paths: dict[str, Any], root: Path) -> None:
    """Ensure outputs.root is usable even when SCRATCH is not exported."""
    outputs = paths.setdefault("outputs", {})
    raw = outputs.get("root")
    if raw in (None, ""):
        outputs["root"] = str(root / "results")
        return

    text = str(raw)
    if "${SCRATCH}" in text or "$SCRATCH" in text:
        scratch = os.environ.get("SCRATCH")
        if scratch:
            text = text.replace("${SCRATCH}", scratch).replace("$SCRATCH", scratch)
        else:
            text = str(root / "results")
            warnings.warn(
                "SCRATCH is not defined; writing generated figures under "
                f"{root / 'results'}.",
                RuntimeWarning,
                stacklevel=2,
            )
    outputs["root"] = text


def load_paths(root: Path) -> tuple[Path, dict[str, Any]]:
    """Load tracked defaults and merge optional local overrides.

    Resolution order:
      1. config/paths.example.yaml supplies tracked Sherlock defaults.
      2. CCPERTURBSEQ_PATHS, when set, supplies an explicit override file.
      3. otherwise config/paths.local.yaml is used when present.
      4. when no local file exists, tracked defaults are used directly.

    The returned configuration path identifies the local override when one was
    used, otherwise the tracked defaults file.
    """
    example = root / "config" / "paths.example.yaml"
    if not example.exists():
        raise FileNotFoundError(f"Missing tracked path defaults: {example}")
    defaults = yaml.safe_load(example.read_text(encoding="utf-8")) or {}

    explicit = os.environ.get("CCPERTURBSEQ_PATHS")
    if explicit:
        config = Path(explicit).expanduser()
        if not config.exists():
            raise FileNotFoundError(
                "CCPERTURBSEQ_PATHS points to a missing file: "
                f"{config}"
            )
        overrides = yaml.safe_load(config.read_text(encoding="utf-8")) or {}
        paths = _merge_config(defaults, overrides)
    else:
        local = root / "config" / "paths.local.yaml"
        if local.exists():
            config = local
            overrides = yaml.safe_load(local.read_text(encoding="utf-8")) or {}
            paths = _merge_config(defaults, overrides)
        else:
            config = example
            paths = defaults
            warnings.warn(
                "config/paths.local.yaml is absent; using the tracked Sherlock "
                "defaults in config/paths.example.yaml. Create paths.local.yaml "
                "only for machine-specific overrides.",
                RuntimeWarning,
                stacklevel=2,
            )

    _resolve_output_root(paths, root)
    return config, paths


def _expand_path_text(value: Any, dotted_key: str) -> str:
    text = str(value).strip()
    if _is_placeholder(text):
        raise ValueError(
            f"Configuration key {dotted_key!r} is still a placeholder: {text}"
        )
    expanded = os.path.expandvars(text)
    unresolved = _UNRESOLVED_ENV_PATTERN.findall(expanded)
    if unresolved:
        raise ValueError(
            f"Configuration key {dotted_key!r} contains unresolved environment "
            f"variables {unresolved}: {text}"
        )
    return expanded


def get_path(
    paths: dict[str, Any],
    dotted_key: str,
    *,
    required: bool = True,
) -> Path | None:
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

    try:
        expanded = _expand_path_text(node, dotted_key)
    except ValueError:
        if required:
            raise
        return None
    return Path(expanded).expanduser()


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
            rows.append(
                {
                    "key": key,
                    "path": str(path),
                    "exists": exists,
                    "is_file": path.is_file() if exists else False,
                    "is_dir": path.is_dir() if exists else False,
                    "size_bytes": (
                        path.stat().st_size
                        if exists and path.is_file()
                        else None
                    ),
                    "error": None if exists else "path does not exist",
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "key": key,
                    "path": None,
                    "exists": False,
                    "is_file": False,
                    "is_dir": False,
                    "size_bytes": None,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
    return pd.DataFrame(rows)


def _format_missing_inputs(
    audit: pd.DataFrame,
    *,
    figure_id: str,
    config: Path,
) -> str:
    missing = audit.loc[~audit["exists"].fillna(False)]
    lines = [
        f"{figure_id} cannot start because {len(missing)} required input(s) "
        "are unavailable.",
        f"Configuration source: {config}",
    ]
    for row in missing.itertuples(index=False):
        detail = row.error or "missing"
        location = row.path or "<unresolved>"
        lines.append(f"  - {row.key}: {location} ({detail})")
    if any(str(key).startswith("hippo_timecourse.") for key in missing["key"]):
        lines.append(
            "Recreate and reconnect the Hippo inputs with Run All in "
            "notebooks/preprocessing/Rebuild_Hippo_Timecourse.ipynb. "
            "Its settings cell identifies the eight processed sublibrary inputs."
        )
    lines.append(
        "Run `python scripts/preflight_main_figure_notebooks.py --check-data` "
        "from the repository root to audit Figures 1–4 together."
    )
    return "\n".join(lines)


def load_figure_manifest(
    figure_id: str,
    root: Path | None = None,
) -> dict[str, Any]:
    root = root or find_repo_root()
    manifest_file = root / "config" / "figure_notebooks.yaml"
    if not manifest_file.exists():
        raise FileNotFoundError(manifest_file)
    payload = yaml.safe_load(manifest_file.read_text(encoding="utf-8")) or {}
    try:
        spec = payload["figures"][figure_id]
    except KeyError as exc:
        raise KeyError(
            f"Figure {figure_id!r} is not defined in {manifest_file}"
        ) from exc
    return {"figure_id": figure_id, **spec}


def initialize_figure(
    manifest: dict[str, Any],
    *,
    require_inputs: bool = True,
) -> dict[str, Any]:
    root = find_repo_root()
    config, paths = load_paths(root)
    output_dir = get_path(paths, "outputs.root") / manifest["figure_id"]
    output_dir.mkdir(parents=True, exist_ok=True)

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.bbox": "tight",
        }
    )

    required = sorted(
        {
            key
            for panel in manifest["panels"]
            for key in panel.get("required_inputs", [])
        }
    )
    input_audit = audit_inputs(paths, required)
    if require_inputs and not input_audit.empty:
        missing_mask = ~input_audit["exists"].fillna(False)
        if bool(missing_mask.any()):
            raise FileNotFoundError(
                _format_missing_inputs(
                    input_audit,
                    figure_id=manifest["figure_id"],
                    config=config,
                )
            )

    panel_table = pd.DataFrame(manifest["panels"])
    panel_table.attrs["figure_id"] = manifest["figure_id"]
    return {
        "root": root,
        "config": config,
        "paths": paths,
        "output_dir": output_dir,
        "git_commit": git_commit(root),
        "panels": panel_table,
        "input_audit": input_audit,
    }


def panel_directory(context: dict[str, Any], panel: str) -> Path:
    safe = "".join(c if c.isalnum() else "_" for c in panel).strip("_") or "all"
    path = context["output_dir"] / f"panel_{safe}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def code_panel_status(
    context: dict[str, Any],
    panel: str,
    description: str,
) -> dict[str, Any]:
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


def sha256_file(
    path: Path,
    *,
    max_bytes: int = 2_000_000_000,
) -> str | None:
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
        "config": str(context["config"]),
        "inputs": inputs,
        "outputs": outputs,
        "parameters": parameters,
    }
    path = panel_directory(context, panel) / "run_receipt.json"
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return path

