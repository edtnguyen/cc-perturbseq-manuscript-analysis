#!/usr/bin/env python3
"""Audit curated main-figure notebooks before expensive Sherlock execution.

Static checks (default):
  * parse and compile every code cell in Figures 1–4;
  * require exactly one repository-path bootstrap cell before support imports;
  * verify every imported repository support module exists;
  * collect every literal get_path(..., "section.key") request;
  * verify those keys have concrete tracked defaults;
  * reject hard-coded /oak or /scratch paths in executable notebook cells;
  * verify the repository-local cNMF workbook and SHA-256.

Sherlock checks (--check-data):
  * merge config/paths.local.yaml over tracked defaults;
  * report all missing required input paths for Figures 1–4 in one pass;
  * verify the output directory can be created.

Runtime checks (--check-runtime):
  * report missing Python packages used by the four notebooks.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
from typing import Any

import yaml


MAIN_NOTEBOOKS = (
    Path("notebooks/main/Figure_01.ipynb"),
    Path("notebooks/main/Figure_02.ipynb"),
    Path("notebooks/main/Figure_03.ipynb"),
    Path("notebooks/main/Figure_04.ipynb"),
)
BOOTSTRAP_TAG = "repo-path-bootstrap"
WORKBOOK = Path("data/figure_02/cNMF_50_0_2_032526_6co.xlsx")
WORKBOOK_SHA256 = "d2fe4834b1efd3c00e01062b56b907a6178bfbc48a40beb5ea7b347538bfc824"
PLACEHOLDER_PREFIXES = ("REPLACE_WITH_", "EMBEDDED_IN_")
ABSOLUTE_SCIENCE_PATH = re.compile(r"""(?P<quote>['"])(?:/oak|/scratch)/[^'"]+(?P=quote)""")
RUNTIME_PACKAGES = (
    "anndata",
    "h5py",
    "logomaker",
    "matplotlib",
    "mudata",
    "numpy",
    "openpyxl",
    "pandas",
    "plotnine",
    "scanpy",
    "scipy",
    "yaml",
)


def repo_root() -> Path:
    root = Path(__file__).resolve().parents[1]
    if not (root / "docs" / "REPRODUCTION_PLAN.md").exists():
        raise FileNotFoundError(f"Could not identify repository root from {__file__}")
    return root


def is_placeholder(value: Any) -> bool:
    if not isinstance(value, str):
        return value is None
    text = value.strip()
    return not text or text.upper().startswith(PLACEHOLDER_PREFIXES)


def dotted_lookup(payload: dict[str, Any], key: str) -> Any:
    node: Any = payload
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node:
            raise KeyError(key)
        node = node[part]
    return node


def deep_merge(defaults: Any, overrides: Any) -> Any:
    if isinstance(defaults, dict) and isinstance(overrides, dict):
        merged = dict(defaults)
        for key, value in overrides.items():
            merged[key] = deep_merge(merged.get(key), value)
        return merged
    if overrides is None or is_placeholder(overrides):
        return defaults
    return overrides


def load_config(root: Path) -> tuple[Path, dict[str, Any]]:
    example = root / "config" / "paths.example.yaml"
    defaults = yaml.safe_load(example.read_text(encoding="utf-8")) or {}
    explicit = os.environ.get("CCPERTURBSEQ_PATHS")
    if explicit:
        local = Path(explicit).expanduser()
        if not local.exists():
            raise FileNotFoundError(f"CCPERTURBSEQ_PATHS does not exist: {local}")
    else:
        local = root / "config" / "paths.local.yaml"

    if local.exists():
        overrides = yaml.safe_load(local.read_text(encoding="utf-8")) or {}
        return local, deep_merge(defaults, overrides)
    return example, defaults


def extract_notebook_facts(
    root: Path,
    relative_path: Path,
) -> tuple[set[str], set[str], list[str]]:
    path = root / relative_path
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("nbformat") != 4:
        raise ValueError(f"{relative_path}: expected nbformat 4")

    code_cells = [
        (index, cell)
        for index, cell in enumerate(payload.get("cells", []))
        if cell.get("cell_type") == "code"
    ]
    bootstrap_cells = [
        index
        for index, cell in code_cells
        if BOOTSTRAP_TAG in cell.get("metadata", {}).get("tags", [])
        or cell.get("id") == BOOTSTRAP_TAG
    ]
    if len(bootstrap_cells) != 1:
        raise ValueError(
            f"{relative_path}: expected one {BOOTSTRAP_TAG!r} cell, "
            f"found {len(bootstrap_cells)}"
        )

    path_keys: set[str] = set()
    support_modules: set[str] = set()
    hard_coded_paths: list[str] = []
    first_support_import: int | None = None

    for index, cell in code_cells:
        source = "".join(cell.get("source", []))
        if not source.strip():
            continue
        tree = ast.parse(source, filename=f"{relative_path}:cell-{index}")
        compile(source, f"{relative_path}:cell-{index}", "exec")

        for match in ABSOLUTE_SCIENCE_PATH.finditer(source):
            hard_coded_paths.append(
                f"cell {index}: {match.group(0)[:180]}"
            )

        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module:
                module = node.module
                if module.startswith("src.cc_perturbseq") or module.startswith(
                    "cc_perturbseq"
                ):
                    support_modules.add(module)
                    if first_support_import is None:
                        first_support_import = index

            if not isinstance(node, ast.Call):
                continue
            func = node.func
            is_get_path = (
                isinstance(func, ast.Name)
                and func.id == "get_path"
                or isinstance(func, ast.Attribute)
                and func.attr == "get_path"
            )
            if not is_get_path or len(node.args) < 2:
                continue
            key_node = node.args[1]
            if isinstance(key_node, ast.Constant) and isinstance(
                key_node.value, str
            ):
                path_keys.add(key_node.value)
            else:
                raise ValueError(
                    f"{relative_path}: get_path key in cell {index} "
                    "must be a string literal for preflight auditing"
                )

    bootstrap_index = bootstrap_cells[0]
    if first_support_import is not None and bootstrap_index > first_support_import:
        raise ValueError(
            f"{relative_path}: bootstrap cell {bootstrap_index} occurs after "
            f"repository support import cell {first_support_import}"
        )

    return path_keys, support_modules, hard_coded_paths


def module_path(root: Path, module: str) -> Path:
    if module == "src.cc_perturbseq":
        return root / "src" / "cc_perturbseq" / "__init__.py"
    if module.startswith("src.cc_perturbseq."):
        suffix = module.removeprefix("src.").replace(".", "/")
        return root / "src" / f"{suffix}.py"
    if module == "cc_perturbseq":
        return root / "src" / "cc_perturbseq" / "__init__.py"
    suffix = module.replace(".", "/")
    return root / "src" / f"{suffix}.py"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expand_path(value: Any, root: Path, key: str) -> Path:
    if is_placeholder(value):
        raise ValueError(f"{key} is unresolved: {value!r}")
    text = str(value)
    if "${SCRATCH}" in text or "$SCRATCH" in text:
        scratch = os.environ.get("SCRATCH")
        text = (
            text.replace("${SCRATCH}", scratch).replace("$SCRATCH", scratch)
            if scratch
            else str(root / "results")
        )
    expanded = os.path.expandvars(text)
    if "$" in expanded:
        raise ValueError(f"{key} has an unresolved environment variable: {text}")
    return Path(expanded).expanduser()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check-data",
        action="store_true",
        help="On Sherlock, verify all configured input paths exist.",
    )
    parser.add_argument(
        "--check-runtime",
        action="store_true",
        help="Verify Python packages used by Figures 1–4 are importable.",
    )
    args = parser.parse_args()

    root = repo_root()
    defaults = yaml.safe_load(
        (root / "config" / "paths.example.yaml").read_text(encoding="utf-8")
    ) or {}

    all_keys: set[str] = set()
    all_modules: set[str] = set()
    failures: list[str] = []

    print("Static audit of main figure notebooks")
    for relative_path in MAIN_NOTEBOOKS:
        try:
            keys, modules, hard_paths = extract_notebook_facts(root, relative_path)
            all_keys.update(keys)
            all_modules.update(modules)
            if hard_paths:
                failures.append(
                    f"{relative_path}: hard-coded Sherlock paths in code: "
                    + "; ".join(hard_paths)
                )
            print(
                f"  PASS {relative_path}: {len(keys)} configured paths, "
                f"{len(modules)} support modules"
            )
        except Exception as exc:
            failures.append(f"{relative_path}: {type(exc).__name__}: {exc}")

    for module in sorted(all_modules):
        path = module_path(root, module)
        if not path.exists():
            failures.append(f"Missing support module for {module}: {path}")

    for key in sorted(all_keys):
        try:
            value = dotted_lookup(defaults, key)
            if is_placeholder(value):
                failures.append(
                    f"Tracked default for main-figure key {key} is unresolved: "
                    f"{value!r}"
                )
        except KeyError:
            failures.append(
                f"Main-figure key {key} is absent from config/paths.example.yaml"
            )

    workbook = root / WORKBOOK
    if not workbook.exists():
        failures.append(f"Missing repository-local workbook: {WORKBOOK}")
    else:
        observed = sha256_file(workbook)
        if observed != WORKBOOK_SHA256:
            failures.append(
                f"Workbook SHA-256 mismatch: observed {observed}, "
                f"expected {WORKBOOK_SHA256}"
            )
        else:
            print(f"  PASS {WORKBOOK}: SHA-256 verified")

    if args.check_runtime:
        missing_packages = [
            package
            for package in RUNTIME_PACKAGES
            if importlib.util.find_spec(package) is None
        ]
        if missing_packages:
            failures.append(
                "Missing Python packages: " + ", ".join(missing_packages)
            )
        else:
            print(
                "  PASS runtime packages: "
                + ", ".join(RUNTIME_PACKAGES)
            )

    if args.check_data:
        config_path, paths = load_config(root)
        print(f"Data audit using {config_path}")
        missing_data: list[str] = []
        for key in sorted(all_keys):
            try:
                path = expand_path(dotted_lookup(paths, key), root, key)
                if not path.exists():
                    missing_data.append(f"{key}: {path}")
                else:
                    kind = "directory" if path.is_dir() else "file"
                    size = path.stat().st_size if path.is_file() else "-"
                    print(f"  PASS {key}: {kind}, size={size}, {path}")
            except Exception as exc:
                missing_data.append(f"{key}: {type(exc).__name__}: {exc}")

        try:
            output_root = expand_path(
                dotted_lookup(paths, "outputs.root"),
                root,
                "outputs.root",
            )
            output_root.mkdir(parents=True, exist_ok=True)
            probe = output_root / ".figure_preflight_write_test"
            probe.write_text("ok\n", encoding="utf-8")
            probe.unlink()
            print(f"  PASS outputs.root is writable: {output_root}")
        except Exception as exc:
            missing_data.append(
                f"outputs.root: {type(exc).__name__}: {exc}"
            )

        if missing_data:
            failures.append(
                "Missing or inaccessible Sherlock inputs:\n    "
                + "\n    ".join(missing_data)
            )

    if failures:
        print("\nPRE-FLIGHT FAILED")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print(
        f"\nPRE-FLIGHT PASSED: {len(MAIN_NOTEBOOKS)} notebooks, "
        f"{len(all_keys)} unique configured inputs, "
        f"{len(all_modules)} support modules."
    )
    if not args.check_data:
        print(
            "On Sherlock, run again with --check-data --check-runtime "
            "before opening the notebooks."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
