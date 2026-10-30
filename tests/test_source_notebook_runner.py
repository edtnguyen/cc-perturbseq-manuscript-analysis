from __future__ import annotations

from pathlib import Path
import sys

import nbformat
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cc_perturbseq.source_notebook_runner import (  # noqa: E402
    identify_s3_bc_cells,
    patch_source_notebook,
    prepare_s3_source_reproduction,
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


def test_s3_source_cells_are_identified_and_only_instrumented(tmp_path):
    source = _write_notebook(
        tmp_path / "source_qc.ipynb",
        [
            "import matplotlib.pyplot as plt\n"
            "total_counts = [1, 2]\n"
            "n_genes_by_counts = [1, 2]\n"
            "round1 = [16, 17]\n"
            "plt.scatter(total_counts, n_genes_by_counts, c=round1)\n"
            "plt.show()\n",
            "pct_counts_mt = [1, 2]\n"
            "pct_counts_ribo = [1, 2]\n"
            "qc_pass = [True, False]\n"
            "fig, axes = plt.subplots(1, 4)\n"
            "axes[0].violinplot(total_counts)\n"
            "axes[1].violinplot(n_genes_by_counts)\n"
            "axes[2].violinplot(pct_counts_mt)\n"
            "axes[3].violinplot(pct_counts_ribo)\n"
            "plt.show()\n",
        ],
    )
    notebook = nbformat.read(source, as_version=4)
    matches = identify_s3_bc_cells(notebook)
    assert [row.cell_index for row in matches["b"]] == [0]
    assert [row.cell_index for row in matches["c"]] == [1]

    prepared_path = tmp_path / "prepared.ipynb"
    prepared, provenance = prepare_s3_source_reproduction(
        source,
        prepared_path,
        capture_output_dir=tmp_path / "captures",
    )
    assert provenance["source_sha256"] == sha256_file(source)
    assert "_cc_capture_s3('b')" in prepared.cells[1].source
    assert "_cc_capture_s3('c')" in prepared.cells[2].source
    assert "Figure_S03_source_capture_manifest.json" in prepared.cells[-1].source


def test_s3_source_without_exact_panels_fails(tmp_path):
    source = _write_notebook(
        tmp_path / "wrong.ipynb",
        ["import matplotlib.pyplot as plt\nplt.plot([1, 2])\n"],
    )
    notebook = nbformat.read(source, as_version=4)
    with pytest.raises(ValueError, match="No generic replacement"):
        identify_s3_bc_cells(notebook)
