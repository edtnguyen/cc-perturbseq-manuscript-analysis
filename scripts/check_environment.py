#!/usr/bin/env python3
"""Check the figure environment without downloading or processing study data."""
from __future__ import annotations

import argparse
import importlib
from importlib import metadata
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
IMPORT_NAMES = {
    "scikit-learn": "sklearn", "umap-learn": "umap", "pyyaml": "yaml",
    "ipython": "IPython", "jupyter-client": "jupyter_client",
}


def check_environment(root=ROOT, *, check_pip=True):
    root = Path(root).resolve()
    requirements = [s.strip() for s in (root / "config/geo-reproduction-requirements.txt").read_text().splitlines()
                    if s.strip() and not s.lstrip().startswith("#")]
    report = {"python": sys.version.split()[0], "executable": sys.executable,
              "expected_python": "3.11", "packages": [], "local_modules": [], "problems": []}
    if sys.version_info[:2] != (3, 11):
        report["problems"].append("Use Python 3.11 from environment.yml, not the base/system Python.")
    try:
        from packaging.requirements import Requirement
    except ImportError:
        report["problems"].append("packaging is missing; create the environment from environment.yml.")
        report["status"] = "FAIL"
        return report
    for line in requirements:
        req = Requirement(line)
        row = {"requirement": line}
        try:
            installed = metadata.version(req.name)
            row["installed"] = installed
            if installed not in req.specifier:
                raise RuntimeError(f"installed {installed}, expected {req.specifier}")
            module = IMPORT_NAMES.get(req.name.lower(), req.name.replace("-", "_"))
            importlib.import_module(module)
            row["status"] = "PASS"
        except Exception as exc:
            row.update(status="FAIL", error=f"{type(exc).__name__}: {exc}")
            report["problems"].append(f"{req.name}: {row['error']}")
        report["packages"].append(row)
    # Import repository code too, so undeclared direct dependencies are caught.
    sys.path.insert(0, str(root / "src"))
    for path in sorted((root / "src/cc_perturbseq").glob("*.py")):
        if path.name == "__init__.py":
            continue
        name = f"cc_perturbseq.{path.stem}"
        try:
            importlib.import_module(name)
            report["local_modules"].append({"module": name, "status": "PASS"})
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
            report["local_modules"].append({"module": name, "status": "FAIL", "error": message})
            report["problems"].append(f"{name}: {message}")
    if check_pip:
        result = subprocess.run([sys.executable, "-m", "pip", "check"],
                                text=True, capture_output=True, timeout=120)
        report["pip_check"] = {"returncode": result.returncode,
                               "message": (result.stdout + result.stderr).strip()}
        if result.returncode:
            report["problems"].append("pip check found broken/missing package dependencies.")
    report["status"] = "PASS" if not report["problems"] else "FAIL"
    report["scope"] = "Python/package versions and imports; not study-data availability or a numerical reproduction run."
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Print machine-readable results")
    parser.add_argument("--report", type=Path, help="Optional JSON report path")
    args = parser.parse_args()
    report = check_environment()
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("Python:", report["python"])
        print("Executable:", report["executable"])
        for row in report["packages"]:
            print(f"{row['status']:4s}  {row['requirement']:32s} {row.get('installed', 'missing')}")
        print("Repository modules checked:", len(report["local_modules"]))
        for problem in report["problems"]:
            print("ERROR:", problem)
        print("ENVIRONMENT", report["status"])
        print(report["scope"])
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
