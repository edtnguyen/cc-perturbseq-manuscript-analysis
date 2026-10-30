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
    motif_source = (
        ROOT / "src" / "cc_perturbseq" / "motif_sources.py"
    ).read_text(encoding="utf-8")
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
    assert (
        "program_annotation_panels_D0_D1_D2_D3_rotated_logos_"
        "D2_extra_D3_trimmed_markers.ipynb"
    ) in text
    assert (
        "b0ea563bee36408e54ff793e4915a7258232473a65813bba2bb5d7b5f655599a"
        in text
    )
    assert CANONICAL in text
    assert "patch_source_notebook" in text
    assert 'source_logic_changed\\": False' in text
    assert "plot_multi_day_annotation" not in text


def test_figure_s4_is_exact_previous_notebook_driver_for_b_c():
    text = _notebook_text(
        ROOT / "notebooks" / "supplementary" / "Figure_S04.ipynb"
    )
    assert "iPSC_EC_QQ_full_covariates_with_replicate_correlation.ipynb" in text
    assert (
        "6e80625abd9b5e71d45da6dd19e888898318d39a159a7ac35600fc1d85489e87"
        in text
    )
    assert '"panels": ["b", "c"]' in text
    assert "collect_s4_source_outputs" in text
    assert 'source_logic_changed\\": False' in text
    assert 'source_style_changed\\": False' in text
    assert "plot_crt_qq" not in text
    assert "plot_replicate_effects" not in text
    assert "plot_cnmf_diagnostics" not in text
    assert "plot_program_gene_overlap" not in text
    assert "plot_usage_correlation" not in text


def test_retained_motif_notebook_code_cells_compile():
    for name in ("Figure_S04.ipynb", "Figure_S06.ipynb"):
        path = ROOT / "notebooks" / "supplementary" / name
        notebook = json.loads(path.read_text(encoding="utf-8"))
        for index, cell in enumerate(notebook["cells"]):
            if cell.get("cell_type") == "code":
                compile(
                    "".join(cell.get("source", [])),
                    f"{name}:cell-{index}",
                    "exec",
                )
