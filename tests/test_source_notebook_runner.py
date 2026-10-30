from __future__ import annotations

from hashlib import sha1
from pathlib import Path
import sys

import nbformat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cc_perturbseq.source_notebook_runner import (  # noqa: E402
    git_blob_sha1_file,
    patch_source_notebook,
    sha256_file,
)


def _write_notebook(path: Path, cells: list[str]) -> Path:
    notebook = nbformat.v4.new_notebook(
        cells=[nbformat.v4.new_code_cell(source) for source in cells],
        metadata={
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            }
        },
    )
    nbformat.write(notebook, path)
    return path


def test_patch_source_notebook_changes_only_named_path_assignments(tmp_path):
    source = _write_notebook(
        tmp_path / "source.ipynb",
        [
            "from pathlib import Path\n"
            "ENHANCER_ENRICHMENT_PATH = Path(\n"
            "    '/old/'\n"
            "    'enrichment.tsv'\n"
            ")\n"
            "UNCHANGED = 7\n"
        ],
    )
    original_sha = sha256_file(source)
    output = tmp_path / "patched.ipynb"
    notebook, provenance = patch_source_notebook(
        source,
        output,
        expected_sha256=original_sha,
        path_replacements={
            "ENHANCER_ENRICHMENT_PATH": "/new/enhancer_enrichment_results.tsv"
        },
    )
    code = notebook.cells[0].source
    assert 'Path("/new/enhancer_enrichment_results.tsv")' in code
    assert "UNCHANGED = 7" in code
    assert provenance["replacement_counts"] == {
        "ENHANCER_ENRICHMENT_PATH": 1
    }


def test_git_blob_sha1_matches_git_object_formula(tmp_path):
    path = tmp_path / "blob.txt"
    payload = b"source notebook\n"
    path.write_bytes(payload)
    expected = sha1(f"blob {len(payload)}\0".encode() + payload).hexdigest()
    assert git_blob_sha1_file(path) == expected
