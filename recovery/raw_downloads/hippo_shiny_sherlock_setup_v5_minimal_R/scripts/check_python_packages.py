#!/usr/bin/env python
"""Check Python dependencies needed to export the Shiny HDF5 bundle."""
import importlib
import sys

packages = [
    "scanpy",
    "muon",
    "h5py",
    "numpy",
    "pandas",
    "scipy",
    "matplotlib",
]

missing = []
for pkg in packages:
    try:
        mod = importlib.import_module(pkg)
        version = getattr(mod, "__version__", "unknown")
        print(f"OK   {pkg:12s} {version}")
    except Exception as exc:  # noqa: BLE001
        print(f"MISS {pkg:12s} {exc}")
        missing.append(pkg)

if missing:
    print("\nMissing packages:", ", ".join(missing), file=sys.stderr)
    sys.exit(1)
print("\nPython export dependencies are available.")
