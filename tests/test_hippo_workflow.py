"""Bounded integration tests for configuration and fresh-kernel execution."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys

import nbformat
from nbclient.exceptions import CellExecutionError
import pytest
import yaml

from cc_perturbseq.figure_notebook import load_paths
from cc_perturbseq.hippo_workflow import write_hippo_overrides, execute_figure_notebooks


@pytest.fixture
def repo(tmp_path, monkeypatch):
    monkeypatch.delenv("CCPERTURBSEQ_PATHS", raising=False)
    root = tmp_path / "repo"
    (root / "config").mkdir(parents=True)
    (root / "notebooks" / "main").mkdir(parents=True)
    (root / "config" / "paths.example.yaml").write_text(yaml.safe_dump({
        "outputs": {"root": str(tmp_path / "scratch_figures")},
        "hippo_timecourse": {"query_h5ad": "/old/query", "reprojected_cnmf_h5ad": "/old/projected"},
        "reference": {"gtf": "/tracked/annotation.gtf"},
    }))
    return root


def write_inputs(tmp_path):
    query, projected = tmp_path / "query.h5ad", tmp_path / "projected.h5ad"
    query.write_bytes(b"validated-query-fixture")
    projected.write_bytes(b"validated-projection-fixture")
    return query, projected


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_overrides_preserve_unrelated_values_exact_backup_and_idempotence(repo, tmp_path):
    original = ("# Keep a verbatim backup of this comment.\n"
                "outputs:\n  root: /my/existing/figures\n"
                "unrelated:\n  nested: [one, two]\n"
                "hippo_timecourse:\n  query_h5ad: /lost/query\n"
                "  reprojected_cnmf_h5ad: /lost/projected\n  extra: preserved\n")
    active = repo / "config" / "paths.local.yaml"
    active.write_text(original)
    example = repo / "config" / "paths.example.yaml"
    example_before = example.read_bytes()
    query, projected = write_inputs(tmp_path)
    backup_dir = tmp_path / "backups"
    result = write_hippo_overrides(repo, query, projected, backup_dir=backup_dir)
    assert Path(result["backup"]).read_text() == original
    changed = yaml.safe_load(active.read_text())
    assert changed["outputs"] == {"root": "/my/existing/figures"}
    assert changed["unrelated"] == {"nested": ["one", "two"]}
    assert changed["hippo_timecourse"]["extra"] == "preserved"
    assert changed["hippo_timecourse"]["query_h5ad"] == str(query)
    assert example.read_bytes() == example_before
    timestamp = active.stat().st_mtime_ns
    repeated = write_hippo_overrides(repo, query, projected, backup_dir=backup_dir)
    assert repeated["backup"] is None
    assert active.stat().st_mtime_ns == timestamp
    assert len(list(backup_dir.iterdir())) == 1


def test_explicit_config_is_updated_while_local_and_defaults_stay_unchanged(repo, tmp_path, monkeypatch):
    local = repo / "config" / "paths.local.yaml"
    local.write_text("custom_local: intact\n")
    explicit = tmp_path / "other-machine.yaml"
    original = "# Explicit overrides\ncustom_explicit: retained\n"
    explicit.write_text(original)
    monkeypatch.setenv("CCPERTURBSEQ_PATHS", str(explicit))
    query, projected = write_inputs(tmp_path)
    before = digest(repo / "config" / "paths.example.yaml")
    result = write_hippo_overrides(repo, query, projected, backup_dir=tmp_path / "backups")
    assert result["config"] == str(explicit)
    assert Path(result["backup"]).read_text() == original
    assert yaml.safe_load(explicit.read_text())["custom_explicit"] == "retained"
    assert local.read_text() == "custom_local: intact\n"
    assert digest(repo / "config" / "paths.example.yaml") == before
    _, effective = load_paths(repo)
    assert effective["reference"]["gtf"] == "/tracked/annotation.gtf"


def test_new_local_override_is_minimal_and_explicit_example_is_never_rewritten(repo, tmp_path, monkeypatch):
    query, projected = write_inputs(tmp_path)
    example = repo / "config" / "paths.example.yaml"
    before = example.read_bytes()
    result = write_hippo_overrides(repo, query, projected, backup_dir=tmp_path / "backups")
    assert result["backup"] is None
    assert set(yaml.safe_load(Path(result["config"]).read_text())) == {"hippo_timecourse"}
    monkeypatch.setenv("CCPERTURBSEQ_PATHS", str(example))
    with pytest.raises(ValueError, match="tracked example"):
        write_hippo_overrides(repo, query, projected, backup_dir=tmp_path / "backups")
    assert example.read_bytes() == before


def make_notebook(root, name, cells):
    notebook = nbformat.v4.new_notebook(
        cells=[nbformat.v4.new_code_cell(source) for source in cells],
        metadata={"kernelspec": {"name": "must-not-use-this-kernel",
                                 "display_name": "Stale source kernel", "language": "python"}},
    )
    notebook.cells[0].execution_count = 999
    notebook.cells[0].outputs = [nbformat.v4.new_output("stream", name="stdout", text="STALE OUTPUT")]
    path = root / "notebooks" / "main" / name
    nbformat.write(notebook, path)
    return path


def test_execution_uses_current_python_fresh_kernels_and_preserves_sources(repo, tmp_path):
    output_root = tmp_path / "scratch_figures"
    notebooks = []
    for number, stem in enumerate(("Figure_03", "Figure_04")):
        source = f"""from pathlib import Path
import sys, os, json
assert sys.executable == {sys.executable!r}
assert Path.cwd() == Path({str(repo)!r})
assert 'sentinel_from_previous_notebook' not in globals()
sentinel_from_previous_notebook = True
out = Path({str(output_root)!r}) / {stem!r}
out.mkdir(parents=True, exist_ok=True)
(out / 'proof.json').write_text(json.dumps({{'pid': os.getpid(), 'python': sys.executable}}))
print('fresh kernel {number}')
"""
        notebooks.append(make_notebook(repo, stem + ".ipynb", [source]))
    before = [digest(path) for path in notebooks]
    run_dir = tmp_path / "run"
    report = execute_figure_notebooks(repo, run_dir)
    assert report["status"] == "executed"
    assert [digest(path) for path in notebooks] == before
    proofs = []
    for stem in ("Figure_03", "Figure_04"):
        archived = run_dir / "figures" / stem / "proof.json"
        proofs.append(json.loads(archived.read_text()))
        assert archived.read_bytes() == (output_root / stem / "proof.json").read_bytes()
        executed = nbformat.read(run_dir / f"{stem}.executed.ipynb", as_version=4)
        assert executed.cells[0].execution_count == 1
        assert "STALE OUTPUT" not in json.dumps(executed)
    assert proofs[0]["pid"] != proofs[1]["pid"]
    assert all(proof["python"] == sys.executable for proof in proofs)
    saved = json.loads((run_dir / "figure_execution_receipt.json").read_text())
    assert saved["status"] == "executed"
    assert all(record["status"] == "executed" for record in saved["notebooks"])


def test_failure_saves_partial_notebook_and_receipt_without_running_later_sources(repo, tmp_path):
    first = make_notebook(repo, "Figure_03.ipynb", ["print('completed-before-error')",
                                                 "raise ValueError('fixture failure')",
                                                 "raise AssertionError('must not run')"])
    second = make_notebook(repo, "Figure_04.ipynb", ["raise AssertionError('must not start')"])
    before = [digest(first), digest(second)]
    run_dir = tmp_path / "failed-run"
    with pytest.raises(CellExecutionError, match="fixture failure"):
        execute_figure_notebooks(repo, run_dir, archive_outputs=False)
    assert [digest(first), digest(second)] == before
    partial = nbformat.read(run_dir / "Figure_03.executed.ipynb", as_version=4)
    assert partial.cells[0].execution_count == 1
    assert partial.cells[1].outputs[-1].output_type == "error"
    assert partial.cells[2].execution_count is None
    assert not (run_dir / "Figure_04.executed.ipynb").exists()
    saved = json.loads((run_dir / "figure_execution_receipt.json").read_text())
    assert saved["status"] == "failed"
    assert len(saved["notebooks"]) == 1
    assert saved["notebooks"][0]["status"] == "failed"
    assert "fixture failure" in saved["notebooks"][0]["error"]


def test_relative_explicit_config_is_pinned_for_changed_kernel_working_directory(repo, tmp_path, monkeypatch):
    caller = tmp_path / "caller"
    caller.mkdir()
    explicit = caller / "custom.yaml"
    explicit.write_text("custom_explicit: retained\n")
    monkeypatch.chdir(caller)
    monkeypatch.setenv("CCPERTURBSEQ_PATHS", "custom.yaml")
    make_notebook(repo, "Figure_03.ipynb", [f"""from pathlib import Path
import os
configured = Path(os.environ['CCPERTURBSEQ_PATHS'])
assert configured.is_absolute()
assert configured == Path({str(explicit)!r})
assert configured.read_text() == 'custom_explicit: retained\\n'
"""])
    report = execute_figure_notebooks(repo, tmp_path / "relative-config-run",
                                     notebook_names=("Figure_03.ipynb",), archive_outputs=False)
    assert report["config"] == str(explicit)
    assert report["status"] == "executed"
    assert os.environ["CCPERTURBSEQ_PATHS"] == "custom.yaml"


@pytest.fixture
def metadata_only_query(tmp_path):
    import h5py
    import pandas as pd
    from anndata.io import write_elem
    from cc_perturbseq.figure3_panels import FIGURE3F_D3_CONDITIONS, FIGURE3F_GENES
    from cc_perturbseq.figure4_panels import FIGURE4C_D12_CONDITIONS, FIGURE4C_GENES

    conditions = tuple(FIGURE3F_D3_CONDITIONS) + tuple(FIGURE4C_D12_CONDITIONS)
    markers = sorted({gene for groups in (FIGURE3F_GENES, FIGURE4C_GENES)
                      for group in groups.values() for gene in group})
    obs = pd.DataFrame({"sample_name": [condition + "_rep1" for condition in conditions],
                        "sublibID": "s1", "pct_counts_mt": 0.1, "pct_counts_ribo": 8.0},
                       index=[f"cell_{i}" for i in range(len(conditions))])
    var = pd.DataFrame({"gene": markers}, index=[f"ENSG{i}" for i in range(len(markers))])
    path = tmp_path / "metadata-only.h5ad"
    def write(obs_table=obs, var_table=var):
        with h5py.File(path, "w") as handle:
            write_elem(handle, "obs", obs_table)
            write_elem(handle, "var", var_table)
        return path
    return write, obs, var


def test_panel_precheck_reads_only_metadata_and_reports_actual_condition_counts(metadata_only_query):
    import h5py
    from cc_perturbseq.hippo_workflow import validate_hippo_panel_inputs
    write, obs, _ = metadata_only_query
    path = write()
    with h5py.File(path) as handle:
        assert "X" not in handle  # Accessing count matrices would fail this fixture.
    report = validate_hippo_panel_inputs(path)
    assert report["n_query_cells"] == len(obs)
    assert set(report["panels"]) == {"Figure_03f", "Figure_04c"}
    assert all(count == 1 for panel in report["panels"].values()
               for count in panel["condition_counts"].values())
    assert report["panels"]["Figure_03f"]["n_cells"] == 5
    assert report["panels"]["Figure_04c"]["n_cells"] == 10


def test_panel_precheck_rejects_missing_marker_or_empty_post_qc_condition(metadata_only_query):
    from cc_perturbseq.hippo_workflow import validate_hippo_panel_inputs
    write, obs, var = metadata_only_query
    missing_marker = var.loc[var["gene"].ne("CDH5")]
    with pytest.raises(ValueError, match="CDH5"):
        validate_hippo_panel_inputs(write(var_table=missing_marker))
    no_condition = obs.copy()
    no_condition.loc[no_condition.index[0], "sublibID"] = "s8"
    with pytest.raises(ValueError, match="no cells"):
        validate_hippo_panel_inputs(write(obs_table=no_condition))
    no_condition = obs.copy()
    no_condition.loc[no_condition.index[0], "pct_counts_ribo"] = 0.0
    with pytest.raises(ValueError, match="no cells"):
        validate_hippo_panel_inputs(write(obs_table=no_condition))
