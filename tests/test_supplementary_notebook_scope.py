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
