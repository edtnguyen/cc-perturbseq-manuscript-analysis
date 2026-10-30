"""Ensure Figure 3/4 preflight differences remain advisory."""
from __future__ import annotations

import json
from pathlib import Path


NOTEBOOK = (
    Path(__file__).resolve().parents[1]
    / "notebooks"
    / "preprocessing"
    / "Rebuild_Hippo_Timecourse.ipynb"
)


def test_figure3_4_preflight_is_advisory_not_a_historical_count_gate():
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    source = "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook["cells"]
    )

    assert 'PANEL_PREFLIGHT_MODE = "advisory"' in source
    assert "RUN_FULL_MAIN_PREFLIGHT = False" in source
    assert "validate_hippo_panel_inputs(QUERY_H5AD)" in source
    assert "except (ValueError, KeyError) as exc:" in source
    assert '"acceptance_role": "advisory_only"' in source
    assert "Continuing because PANEL_PREFLIGHT_MODE='advisory'" in source
    assert "Exact historical " in source
    assert "panel cell totals are not required." in source


def test_advisory_notebook_code_cells_compile():
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    for index, cell in enumerate(notebook["cells"]):
        if cell.get("cell_type") == "code":
            compile(
                "".join(cell.get("source", [])),
                f"{NOTEBOOK.name}:cell-{index}",
                "exec",
            )
