from __future__ import annotations

import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
SUPPLEMENTARY_DIR = ROOT / "notebooks" / "supplementary"
RETAINED_NOTEBOOKS = {
    "Figure_S02.ipynb",
    "Figure_S04.ipynb",
    "Figure_S06.ipynb",
}
RETAINED_MANIFEST_KEYS = {
    "Figure_S02",
    "Figure_S04",
    "Figure_S06",
}


def test_only_requested_supplementary_notebooks_are_tracked():
    observed = {
        path.name
        for path in SUPPLEMENTARY_DIR.glob("Figure_S*.ipynb")
        if path.is_file()
    }
    assert observed == RETAINED_NOTEBOOKS


def test_figure_manifest_contains_only_requested_supplementary_figures():
    manifest = yaml.safe_load(
        (ROOT / "config" / "figure_notebooks.yaml").read_text(encoding="utf-8")
    )
    observed = {
        key for key in manifest["figures"] if key.startswith("Figure_S")
    }
    assert observed == RETAINED_MANIFEST_KEYS


def test_retained_supplementary_notebooks_are_valid_json_and_compile():
    for name in sorted(RETAINED_NOTEBOOKS):
        path = SUPPLEMENTARY_DIR / name
        notebook = json.loads(path.read_text(encoding="utf-8"))

        assert notebook["nbformat"] == 4
        assert isinstance(notebook.get("cells"), list)
        assert notebook["cells"], f"{name} contains no cells"

        for index, cell in enumerate(notebook["cells"]):
            assert cell.get("cell_type") in {"markdown", "code", "raw"}
            assert isinstance(cell.get("source", []), (str, list))
            if cell.get("cell_type") == "code":
                compile(
                    "".join(cell.get("source", []))
                    if isinstance(cell.get("source", []), list)
                    else cell.get("source", ""),
                    f"{name}:cell-{index}",
                    "exec",
                )
