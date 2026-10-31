"""Configuration and notebook execution for the Hippo reconstruction workflow.

Scientific reconstruction lives in hippo_rebuild and hippo_projection. This
module connects their validated outputs to the existing figure notebooks.
"""
from __future__ import annotations

from datetime import datetime, timezone
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def _now():
    return datetime.now(timezone.utc).isoformat()


def _atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(text)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def check_runtime(*, include_figures=True):
    """Import actual dependencies early; report all failures in one message."""
    packages = {
        "anndata": "anndata", "h5py": "h5py", "numpy": "numpy",
        "pandas": "pandas", "scipy": "scipy", "sklearn": "scikit-learn",
        "yaml": "PyYAML", "matplotlib": "matplotlib", "plotnine": "plotnine",
    }
    if include_figures:
        packages.update({
            "nbclient": "nbclient", "nbformat": "nbformat",
            "ipykernel": "ipykernel", "scanpy": "scanpy",
            "mudata": "mudata", "logomaker": "logomaker",
            "openpyxl": "openpyxl",
        })
    failures, versions = [], {"python": sys.version, "executable": sys.executable}
    for module, distribution in packages.items():
        try:
            importlib.import_module(module)
            versions[distribution] = importlib.metadata.version(distribution)
        except Exception as exc:
            failures.append(f"{distribution}: {type(exc).__name__}: {exc}")
    if failures:
        raise RuntimeError(
            "The current notebook kernel is missing or cannot import required packages:\n  "
            + "\n  ".join(failures)
            + "\nUse the original Sherlock sc-dl environment and install the missing "
              "packages there before Run All."
        )
    return versions



def validate_hippo_panel_inputs(query_h5ad):
    """Check marker/condition coverage before installing rebuilt figure paths.

    Read only HDF5 observation and feature tables. The potentially large X,
    count layers, and guide matrices are never loaded by this check. The panel
    definitions and source QC are imported from the same code used for plotting.
    """
    from .hippo_projection import _read_tables, add_projection_metadata
    from .figure3_panels import (
        FIGURE3F_GENES, FIGURE3F_D3_CONDITIONS, _hippo_panel_mask,
    )
    from .figure4_panels import FIGURE4C_GENES, FIGURE4C_D12_CONDITIONS

    obs, var = _read_tables(query_h5ad)
    if not obs.index.is_unique or not var.index.is_unique:
        raise ValueError("Hippo query cell IDs and feature IDs must be unique.")
    if "sublibID" not in obs:
        raise ValueError("Hippo query lacks sublibID metadata.")
    allowed_sublibraries = {f"s{i}" for i in range(1, 9)}
    if obs["sublibID"].isna().any():
        raise ValueError("Hippo query has missing sublibrary identifiers.")
    unknown = set(obs["sublibID"].astype(str)) - allowed_sublibraries
    if unknown:
        raise ValueError(f"Hippo query has unrecognized sublibrary identifiers: {sorted(unknown)}")
    if "gene" not in var:
        raise ValueError("Hippo query requires gene symbols in var['gene'].")
    gene_symbols = set(var["gene"].dropna().astype(str))
    metadata = add_projection_metadata(obs)
    panels = {}
    for label, groups, conditions in (
        ("Figure_03f", FIGURE3F_GENES, FIGURE3F_D3_CONDITIONS),
        ("Figure_04c", FIGURE4C_GENES, FIGURE4C_D12_CONDITIONS),
    ):
        markers = [gene for group in groups.values() for gene in group]
        missing = [gene for gene in markers if gene not in gene_symbols]
        if missing:
            raise ValueError(f"Hippo query lacks {label} marker genes: {missing}")
        keep = _hippo_panel_mask(metadata, conditions)
        counts = metadata.loc[keep, "condition"].astype(str).value_counts()
        panels[label] = {
            "markers": markers,
            "n_cells": int(keep.sum()),
            "condition_counts": {condition: int(counts[condition]) for condition in conditions},
        }
    return {
        "query_h5ad": str(Path(query_h5ad).resolve()),
        "n_query_cells": len(obs),
        "n_condition_specific_qc_passing": int(metadata["qc_pass_condition_specific"].sum()),
        "excluded_sublibrary": "s8",
        "panels": panels,
        "validation": "marker availability and all requested conditions after source QC and s8 exclusion",
    }


def write_hippo_overrides(repo_root, query_h5ad, projected_h5ad, *, backup_dir):
    """Update only the two Hippo keys in the active machine-specific config.

The caller must validate the two H5ADs first. Existing unrelated values are
retained, and the previous exact YAML text is backed up before an atomic write.
An explicit CCPERTURBSEQ_PATHS override continues to be honored.
"""
    import yaml
    from .figure_notebook import get_path, load_paths

    root = Path(repo_root).resolve()
    query = Path(query_h5ad).expanduser().resolve()
    projected = Path(projected_h5ad).expanduser().resolve()
    for path in (query, projected):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Cannot configure a missing or empty H5AD: {path}")
    active, _ = load_paths(root)
    example = (root / "config" / "paths.example.yaml").resolve()
    if active.resolve() == example:
        if os.environ.get("CCPERTURBSEQ_PATHS"):
            raise ValueError(
                "CCPERTURBSEQ_PATHS points at the tracked example YAML. Unset it "
                "or select a machine-specific override before rebuilding."
            )
        active = root / "config" / "paths.local.yaml"
    active = active.resolve()
    original = active.read_text(encoding="utf-8") if active.exists() else ""
    overrides = yaml.safe_load(original) or {}
    if not isinstance(overrides, dict):
        raise ValueError(f"Expected a YAML mapping in {active}")
    hippo = overrides.setdefault("hippo_timecourse", {})
    if not isinstance(hippo, dict):
        raise ValueError(f"hippo_timecourse must be a mapping in {active}")
    desired = {"query_h5ad": str(query), "reprojected_cnmf_h5ad": str(projected)}
    backup = None
    if any(hippo.get(key) != value for key, value in desired.items()):
        if original:
            backup = Path(backup_dir) / (
                datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ_") + active.name
            )
            _atomic_text(backup, original)
        hippo.update(desired)
        _atomic_text(active, yaml.safe_dump(overrides, sort_keys=False))
    _, effective = load_paths(root)
    for key, value in desired.items():
        if get_path(effective, f"hippo_timecourse.{key}").resolve() != Path(value):
            raise RuntimeError(f"Effective configuration did not resolve {key} to {value}")
    return {"config": str(active), "backup": str(backup) if backup else None, **desired}


def execute_figure_notebooks(repo_root, run_dir, *,
                             notebook_names=("Figure_03.ipynb", "Figure_04.ipynb"),
                             archive_outputs=True):
    """Run each figure in a fresh kernel using this exact Python interpreter.

Executed copies (including failed cells) and a receipt are saved in run_dir.
Tracked notebooks stay unexecuted. Successful figure outputs are copied to the
run directory so a future scratch cleanup cannot remove that run's figures.
"""
    import nbformat
    from nbclient import NotebookClient
    from jupyter_client import KernelManager
    from jupyter_client.kernelspec import KernelSpecManager
    from .figure_notebook import get_path, git_commit, load_paths

    root, run_dir = Path(repo_root).resolve(), Path(run_dir).resolve()
    active_config, _ = load_paths(root)
    active_config = active_config.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    sources = []
    for name in notebook_names:
        if Path(name).name != name or not name.endswith(".ipynb"):
            raise ValueError(f"Expected a main-notebook filename, received {name!r}")
        source = root / "notebooks" / "main" / name
        if not source.is_file():
            raise FileNotFoundError(source)
        sources.append(source)
    report = {
        "status": "running", "started_utc": _now(), "git_commit": git_commit(root),
        "python_executable": sys.executable, "config": str(active_config), "notebooks": [],
    }
    receipt = run_dir / "figure_execution_receipt.json"

    def save_report():
        _atomic_text(receipt, json.dumps(report, indent=2) + "\n")

    save_report()
    # A temporary kernelspec avoids accidentally using another conda environment.
    # No Jupyter kernel is installed into the user's global kernel registry.
    with tempfile.TemporaryDirectory(prefix="hippo-kernel-") as kernel_root:
        kernel_dir = Path(kernel_root) / "hippo-current-python"
        kernel_dir.mkdir()
        spec = {
            "argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
            "display_name": "Hippo rebuild current Python", "language": "python",
        }
        if os.environ.get("CCPERTURBSEQ_PATHS"):
            # Figure kernels run from root, which can differ from the current
            # notebook directory. Preserve the same explicit file in either cwd.
            spec["env"] = {"CCPERTURBSEQ_PATHS": str(active_config)}
        (kernel_dir / "kernel.json").write_text(json.dumps(spec), encoding="utf-8")
        for source in sources:
            print(f"Executing {source.name} with {sys.executable}", flush=True)
            notebook = nbformat.read(source, as_version=4)
            for cell in notebook.cells:
                if cell.cell_type == "code":
                    cell.outputs, cell.execution_count = [], None
            notebook.metadata.kernelspec = {
                "name": "hippo-current-python", "display_name": spec["display_name"],
                "language": "python",
            }
            destination = run_dir / f"{source.stem}.executed.ipynb"
            record = {"source": str(source), "executed": str(destination),
                      "status": "running", "started_utc": _now()}
            report["notebooks"].append(record)
            save_report()
            manager = KernelManager(
                kernel_name="hippo-current-python",
                kernel_spec_manager=KernelSpecManager(kernel_dirs=[kernel_root]),
            )
            try:
                client = NotebookClient(
                    notebook, km=manager, timeout=None, allow_errors=False,
                    resources={"metadata": {"path": str(root)}},
                )
                client.execute(cleanup_kc=True)
                record["status"] = "executed"
                if archive_outputs:
                    _, paths = load_paths(root)
                    figures = get_path(paths, "outputs.root") / source.stem
                    if not figures.is_dir():
                        raise FileNotFoundError(f"Notebook did not produce its figure directory: {figures}")
                    durable = run_dir / "figures" / source.stem
                    shutil.copytree(figures, durable, dirs_exist_ok=True)
                    record["figure_outputs"] = str(figures)
                    record["archived_figure_outputs"] = str(durable)
            except BaseException as exc:
                record["status"] = "failed"
                record["error"] = f"{type(exc).__name__}: {exc}"
                report["status"] = "failed"
                raise
            finally:
                nbformat.write(notebook, destination)
                record["finished_utc"] = _now()
                report["finished_utc"] = _now()
                save_report()
    report["status"] = "executed"
    report["finished_utc"] = _now()
    report["note"] = "Execution succeeded; comparison with manuscript reference figures is a separate review."
    save_report()
    print(f"Executed notebooks and figure copies: {run_dir}", flush=True)
    return report


def run_main_preflight(repo_root, log_path):
    """Save the full existing main-figure preflight output and fail on an error."""
    root = Path(repo_root).resolve()
    result = subprocess.run(
        [sys.executable, str(root / "scripts" / "preflight_main_figure_notebooks.py"),
         "--check-data", "--check-runtime"],
        cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    _atomic_text(log_path, result.stdout)
    print(result.stdout)
    if result.returncode:
        raise RuntimeError(f"Main-figure preflight failed. See {log_path}")
    return {"status": "passed", "log": str(log_path)}
