import json
from pathlib import Path


NOTEBOOK = (
    Path(__file__).resolve().parents[1]
    / "notebooks"
    / "preprocessing"
    / "Rebuild_Hippo_Timecourse.ipynb"
)


def _sources():
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return notebook, "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook["cells"]
    )


def test_rebuild_notebook_uses_regenerated_sublibraries():
    notebook, source = _sources()

    assert notebook["metadata"]["cc_perturbseq"]["workflow"] == (
        "hippo_timecourse_rebuild_from_regenerated_sublibraries"
    )
    assert 'get_path(paths, "hippo_rebuild.input_root")' in source
    assert "preflight_sublibraries(INPUT_ROOT)" in source
    assert "rebuild_combined" in source
    assert "Hippo:eTN8_GEX_sub{1..8}" in source
    assert "kb_all_main_raw" in source
    assert "counts_filtered" in source

    # The obsolete Oak combined object and its temporary metadata adapter must
    # not be used as the rebuild source.
    assert "SOURCE_QUERY_H5AD" not in source
    assert "REPAIRED_QUERY_H5AD" not in source
    assert "20260325_Hippo_timecourse/Hippo_timecourse_combined.h5ad" not in source


def test_historical_counts_are_diagnostic_but_reference_stays_strict():
    _, source = _sources()

    assert "ENFORCE_EXACT_HISTORICAL_COUNTS = False" in source
    assert "STRICT_REFERENCE = True" in source
    assert 'read_spectra(SPECTRA_TPM, strict_historical=STRICT_REFERENCE)' in source
    assert "historical_count_diagnostics" in source
    assert '"acceptance_role": "diagnostic_only"' in source
    assert 'validate_combined(QUERY_H5AD, strict_historical=False)' in source
    assert "combined_receipt_path" not in source
    assert "except Exception" not in source
    assert "QUERY_H5AD.unlink" not in source


def test_projection_is_non_strict_for_cell_totals_but_fully_validated():
    _, source = _sources()

    assert "project_top_specific_200" in source
    assert "strict_historical=False" in source
    assert "validate_projected" in source
    assert "projection_historical_comparison" in source
    assert "validate_hippo_panel_inputs" in source
    assert "write_hippo_overrides" in source
    assert 'notebook_names=("Figure_03.ipynb", "Figure_04.ipynb")' in source


def test_rebuild_notebook_replaces_stale_artifacts_and_runs_figures():
    _, source = _sources()

    assert "REBUILD_STALE_OUTPUTS = True" in source
    assert "FORCE_REBUILD_ALL = False" in source
    assert "except FileExistsError as exc:" in source
    assert "overwrite=True" in source
    assert "execute_figure_notebooks" in source


def test_rebuild_notebook_code_cells_compile():
    notebook, _ = _sources()

    for index, cell in enumerate(notebook["cells"]):
        if cell.get("cell_type") != "code":
            continue
        compile(
            "".join(cell.get("source", [])),
            f"{NOTEBOOK.name}:cell-{index}",
            "exec",
        )
