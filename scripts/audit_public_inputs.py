#!/usr/bin/env python3
"""Audit active figure code and publication readiness without original data paths."""
from pathlib import Path
import argparse
import ast
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def static_audit(root=ROOT):
    root = Path(root)
    expected = [*(f"notebooks/main/Figure_{i:02d}.ipynb" for i in range(1, 5)),
                *(f"notebooks/supplementary/Figure_S{i:02d}.ipynb" for i in (2, 4, 6))]
    problems, checks = [], []
    for relative in expected:
        notebook = json.loads((root / relative).read_text())
        count = 0
        for i, cell in enumerate(notebook["cells"]):
            if cell["cell_type"] != "code":
                continue
            source = "".join(cell["source"])
            tree = ast.parse(source)
            compile(source, f"{relative}:cell-{i}", "exec")
            for node in ast.walk(tree):
                if isinstance(node, ast.Constant) and isinstance(node.value, str):
                    if node.value.startswith(("/oak/", "/scratch/users/", "/home/users/")):
                        problems.append(f"{relative}:cell-{i}: institutional path literal")
                if isinstance(node, ast.Call):
                    name = node.func.id if isinstance(node.func, ast.Name) else getattr(node.func, "attr", "")
                    if name in {"execute_notebook", "run_path", "patch_source_notebook", "resolve_s4_source_notebook"}:
                        problems.append(f"{relative}:cell-{i}: archived-notebook dependency")
            count += 1
        checks.append({"notebook": relative, "compiled_code_cells": count})
    # These are the executable default-path providers for the retained figures.
    for relative in ["src/cc_perturbseq/motif_sources.py", "src/cc_perturbseq/crt_qq_backend.py", "config/paths.example.yaml"]:
        text = (root / relative).read_text()
        if relative.endswith(".py"):
            tree = ast.parse(text)
            values = [n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)]
            bad = any(v.startswith(("/oak/", "/scratch/users/", "/home/users/")) for v in values)
        else:
            bad = any(v in text for v in ("/oak/", "/scratch/users/", "/home/users/"))
        if bad:
            problems.append(f"{relative}: institutional runtime default")
    return {"status": "PASS" if not problems else "FAIL", "notebooks": checks, "problems": problems,
            "scope": "Active figure source/default paths, not preserved historical output text or the maintainer-only collection manifest."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-data", action="store_true", help="Check companion checksums and GEO presence without downloading")
    args = parser.parse_args()
    result = static_audit()
    if args.check_data:
        from cc_perturbseq.public_data import load_public_paths, require_bundle, require_geo, value
        _, paths = load_public_paths(ROOT)
        readiness = {}
        manifest = value(paths, "public_data.bundle_manifest")
        try:
            keys = [record["key"] for record in json.loads(manifest.read_text())["files"]]
            require_bundle(ROOT, paths, keys)
            readiness["companions"] = f"PASS: {len(keys)} checksummed files"
        except Exception as exc:
            readiness["companions"] = f"MISSING_OR_INVALID: {exc}"
        for dataset, key in (("Hippo", "geo.hippo_h5mu"), ("mixed_cells", "geo.mixed_h5mu")):
            try:
                require_geo(value(paths, key), dataset)
                readiness[dataset] = "PRESENT: processed H5MU; schema checked when read"
            except Exception as exc:
                readiness[dataset] = f"MISSING_OR_INVALID: {exc}"
        result["data_readiness"] = readiness
        if any(v.startswith("MISSING_OR_INVALID") for v in readiness.values()):
            result["status"] = "INCOMPLETE_DATA"
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
