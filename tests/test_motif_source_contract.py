from __future__ import annotations

import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CANONICAL = (
    "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/"
    "analysis/2.tf_enrichment_in_programs/enhancer_enrichment_results.tsv"
)


def _notebook_text(path: Path) -> str:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    return "\n".join(
        "".join(cell.get("source", [])) for cell in notebook["cells"]
    )


def test_tracked_motif_path_is_exact_authorized_table():
    config = yaml.safe_load(
        (ROOT / "config" / "paths.example.yaml").read_text(encoding="utf-8")
    )
    assert config["multiome"]["enhancer_enrichment_results"] == CANONICAL


def test_active_motif_consumers_reject_stale_tables():
    paths = [
        ROOT / "src" / "cc_perturbseq" / "motif_sources.py",
        ROOT / "src" / "cc_perturbseq" / "figure3_panels.py",
        ROOT / "src" / "cc_perturbseq" / "figure4_panels.py",
        ROOT / "src" / "cc_perturbseq" / "program_panels.py",
    ]
    combined = "\n".join(path.read_text(encoding="utf-8") for path in paths)
    motif_source = (ROOT / "src" / "cc_perturbseq" / "motif_sources.py").read_text(encoding="utf-8")
    assert "multiome_ipsc_ec" in motif_source
    assert "enhancer_enrichment_results.tsv" in motif_source
    assert "canonical_enhancer_enrichment_path" in combined
    assert "enhancer_enrichment_results_min10.tsv" not in (
        ROOT / "config" / "paths.example.yaml"
    ).read_text(encoding="utf-8")


def test_figure_s6_runs_recovered_notebook_with_path_only_changes():
    text = _notebook_text(
        ROOT / "notebooks" / "supplementary" / "Figure_S06.ipynb"
    )
    assert "program_annotation_panels_D0_D1_D2_D3_rotated_logos_D2_extra_D3_trimmed_markers.ipynb" in text
    assert "b0ea563bee36408e54ff793e4915a7258232473a65813bba2bb5d7b5f655599a" in text
    assert CANONICAL in text
    assert "patch_source_notebook" in text
    assert "source_logic_changed\": False" in text
    assert "plot_multi_day_annotation" not in text


def test_figure_s3_is_panels_b_c_and_not_generic_renderer():
    text = _notebook_text(
        ROOT / "notebooks" / "supplementary" / "Figure_S03.ipynb"
    )
    assert '"panel": "b"' in text
    assert '"panel": "c"' in text
    assert '"panel": "a"' not in text
    assert '"panel": "d"' not in text
    assert "prepare_s3_source_reproduction" in text
    assert "CCPERTURBSEQ_FIGURE_S03_SOURCE_NOTEBOOK" in text
    assert "load_qc_bundle" not in text
    assert "plot_umi_gene_scatter" not in text
    assert "plot_metric_distributions" not in text


def test_notebook_code_cells_compile():
    for name in ("Figure_S03.ipynb", "Figure_S06.ipynb"):
        path = ROOT / "notebooks" / "supplementary" / name
        notebook = json.loads(path.read_text(encoding="utf-8"))
        for index, cell in enumerate(notebook["cells"]):
            if cell.get("cell_type") == "code":
                compile(
                    "".join(cell.get("source", [])),
                    f"{name}:cell-{index}",
                    "exec",
                )
