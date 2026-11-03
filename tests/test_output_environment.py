"""Output routing and complete installation manifest; no study data required."""
import ast
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import pytest
import yaml
from packaging.requirements import Requirement
from packaging.utils import canonicalize_name

ROOT = Path(__file__).resolve().parents[1]
FIGURES = [f"Figure_{i:02d}" for i in range(1, 5)] + ["Figure_S02", "Figure_S04", "Figure_S06"]


def checkout(tmp_path):
    root = tmp_path / "checkout"
    (root / "config").mkdir(parents=True)
    shutil.copy2(ROOT / "config/paths.example.yaml", root / "config/paths.example.yaml")
    (root / "docs").mkdir()
    (root / "docs/REPRODUCTION_PLAN.md").write_text("Fixture root marker\n")
    return root


@pytest.mark.parametrize("figure", FIGURES)
def test_seven_initializations_write_only_repository_output(tmp_path, monkeypatch, figure):
    from cc_perturbseq import figure_notebook as n, public_data as p
    import scanpy as sc
    root = checkout(tmp_path)
    nested = root / "notebooks" / ("supplementary" if "_S" in figure else "main")
    nested.mkdir(parents=True)
    monkeypatch.chdir(nested)
    monkeypatch.setenv("SCRATCH", str(tmp_path / "must-not-write-scratch"))
    monkeypatch.delenv("CCPERTURBSEQ_PUBLIC_PATHS", raising=False)
    monkeypatch.delenv("CCPERTURBSEQ_PATHS", raising=False)
    monkeypatch.setattr(p, "prepare_figure_inputs", lambda *a, **k: {"fixture": True})
    ctx = n.initialize_figure({"figure_id": figure, "panels": [{"panel": "b", "required_inputs": []}]})
    assert ctx["output_dir"] == root / "output" / figure
    assert Path(sc.settings.figdir) == ctx["output_dir"]
    assert Path(sc.settings.cachedir).is_relative_to(root / "output")
    path = n.panel_directory(ctx, "b") / "source.tsv"
    path.write_text("gene\tvalue\nA\t1\n")
    assert path.is_relative_to(root / "output")
    assert not (root / "results").exists()
    assert not (nested / "output").exists()
    assert not (tmp_path / "must-not-write-scratch").exists()
    for key in ("public_data.derived_cache", "screen.validation_guide_counts",
                "ccperturbseq.integrated_anndata", "ccperturbseq.cnmf_h5mu",
                "hippo_timecourse.query_h5ad", "hippo_timecourse.reprojected_cnmf_h5ad"):
        assert p.value(ctx["paths"], key).is_relative_to(root / "output")


@pytest.mark.parametrize("override", [
    {"outputs": {"root": "results"}},
    {"outputs": {"root": "/tmp/external-results"}},
    {"public_data": {"derived_cache": "data/derived"}},
    {"ccperturbseq": {"integrated_anndata": "data/source.h5ad"}},
])
def test_output_override_cannot_escape_output(tmp_path, monkeypatch, override):
    from cc_perturbseq.public_data import load_public_paths
    root = checkout(tmp_path)
    path = tmp_path / "override.yaml"
    path.write_text(yaml.safe_dump(override))
    monkeypatch.setenv("CCPERTURBSEQ_PUBLIC_PATHS", str(path))
    with pytest.raises(ValueError, match="output"):
        load_public_paths(root)


def test_output_symlink_cannot_redirect_to_scratch(tmp_path, monkeypatch):
    from cc_perturbseq.public_data import load_public_paths
    root = checkout(tmp_path)
    outside = tmp_path / "scratch"; outside.mkdir()
    (root / "output").symlink_to(outside, target_is_directory=True)
    monkeypatch.delenv("CCPERTURBSEQ_PUBLIC_PATHS", raising=False)
    with pytest.raises(ValueError, match="symlink"):
        load_public_paths(root)


def test_yaml_and_requirements_have_identical_python_dependencies():
    env = yaml.safe_load((ROOT / "environment.yml").read_text())
    assert env["name"] == "ccperturbseq-figures"
    assert env["channels"] == ["conda-forge", "nodefaults"]
    assert "python=3.11" in env["dependencies"]
    requirements = [line.strip() for line in (ROOT / "config/geo-reproduction-requirements.txt").read_text().splitlines()
                    if line.strip() and not line.startswith("#")]
    pip = next(x["pip"] for x in env["dependencies"] if isinstance(x, dict))
    assert set(pip) == set(requirements)
    assert len(pip) == len(set(pip))
    names = {canonicalize_name(Requirement(x).name) for x in pip}
    assert {"jupyterlab", "ipython", "ipykernel", "nbclient", "nbformat", "statsmodels", "logomaker"} <= names


def test_all_direct_imports_covered_by_environment():
    env = yaml.safe_load((ROOT / "environment.yml").read_text())
    requirements = next(x["pip"] for x in env["dependencies"] if isinstance(x, dict))
    names = {canonicalize_name(Requirement(x).name) for x in requirements}
    translate = {"yaml": "pyyaml", "sklearn": "scikit-learn", "IPython": "ipython",
                 "umap": "umap-learn", "jupyter_client": "jupyter-client"}
    sources = list((ROOT / "src/cc_perturbseq").glob("*.py")) + list((ROOT / "notebooks").rglob("*.ipynb"))
    for path in sources:
        codes = [path.read_text()] if path.suffix == ".py" else ["".join(c["source"]) for c in json.loads(path.read_text())["cells"] if c["cell_type"] == "code"]
        for code in codes:
            for node in ast.walk(ast.parse(code)):
                imports = ([x.name for x in node.names] if isinstance(node, ast.Import) else
                           [node.module] if isinstance(node, ast.ImportFrom) and not node.level and node.module else [])
                for name in imports:
                    top = name.split(".")[0]
                    if top not in sys.stdlib_module_names and top not in {"cc_perturbseq", "src"}:
                        assert canonicalize_name(translate.get(top, top)) in names, (path, name)


def test_all_seven_notebooks_document_and_print_output():
    import nbformat
    paths = list((ROOT / "notebooks").rglob("*.ipynb"))
    assert len(paths) == 7
    for path in paths:
        nb = nbformat.read(path, as_version=4)
        nbformat.validate(nb)
        text = "\n".join(c.source for c in nb.cells)
        assert f"output/{path.stem}/" in text
        assert 'print("Output directory:", ctx["output_dir"])' in text
        for c in nb.cells:
            if c.cell_type == "code":
                compile(c.source, str(path), "exec")
