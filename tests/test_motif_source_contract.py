from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cc_perturbseq.motif_sources import (  # noqa: E402
    CANONICAL_CLUSTER_MOTIF_MAP_PATH,
    CANONICAL_ENHANCER_ENRICHMENT_PATH,
    CANONICAL_MODISCO_COMPILED_H5_PATH,
    canonical_enhancer_enrichment_path,
    canonical_motif_input_manifest,
    canonical_motif_logo_inputs,
)

ENRICHMENT = (
    "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/"
    "analysis/2.tf_enrichment_in_programs/"
    "enhancer_enrichment_results_min10.tsv"
)
CLUSTER_MAP = (
    "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/"
    "analysis/0.0.motif_id_annotation/tables/"
    "cluster_to_tf_motif_id_mapping.tsv"
)
MODISCO_H5 = (
    "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/"
    "results/compendium/modisco_compiled/modisco_compiled.h5"
)


def _notebook_text(path: Path) -> str:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    return "\n".join(
        "".join(cell.get("source", [])) for cell in notebook["cells"]
    )


def test_tracked_motif_paths_are_exact_selected_files():
    config = yaml.safe_load(
        (ROOT / "config" / "paths.example.yaml").read_text(encoding="utf-8")
    )
    multiome = config["multiome"]
    assert multiome["enhancer_enrichment_results"] == ENRICHMENT
    assert multiome["cluster_motif_map"] == CLUSTER_MAP
    assert multiome["modisco_compiled_h5"] == MODISCO_H5
    assert str(CANONICAL_ENHANCER_ENRICHMENT_PATH) == ENRICHMENT
    assert str(CANONICAL_CLUSTER_MOTIF_MAP_PATH) == CLUSTER_MAP
    assert str(CANONICAL_MODISCO_COMPILED_H5_PATH) == MODISCO_H5


def test_canonical_manifest_contains_all_three_motif_specific_inputs():
    manifest = canonical_motif_input_manifest(require_exists=False)
    assert {key: str(value) for key, value in manifest.items()} == {
        "enhancer_enrichment_results": ENRICHMENT,
        "cluster_motif_map": CLUSTER_MAP,
        "modisco_compiled_h5": MODISCO_H5,
    }


def test_motif_contract_rejects_archived_and_profile_alternatives():
    with pytest.raises(ValueError, match="active min10"):
        canonical_enhancer_enrichment_path(
            "/oak/stanford/groups/engreitz/Users/opushkar/multiome_ipsc_ec/"
            "analysis/archive_v1/2.tf_enrichment_in_programs/"
            "enhancer_enrichment_results.tsv",
            require_exists=False,
        )
    with pytest.raises(ValueError, match="Motif-logo inputs"):
        canonical_motif_logo_inputs(
            configured_cluster_map=(
                "/oak/stanford/groups/engreitz/Users/opushkar/"
                "multiome_ipsc_ec/analysis/archive_v1/tables/"
                "cluster_to_tf_motif_id_mapping.tsv"
            ),
            configured_modisco_h5=(
                "/oak/stanford/groups/engreitz/Users/opushkar/"
                "multiome_ipsc_ec/results/compendium/"
                "modisco_compiled_profile/modisco_compiled.h5"
            ),
            require_exists=False,
        )


def test_figure3_and_figure4_lock_all_three_motif_inputs():
    for name in ("figure3_panels.py", "figure4_panels.py"):
        text = (
            ROOT / "src" / "cc_perturbseq" / name
        ).read_text(encoding="utf-8")
        assert "canonical_enhancer_enrichment_path" in text
        assert "canonical_motif_logo_inputs" in text
        assert "configured_cluster_map=cluster_map_path" in text
        assert "configured_modisco_h5=modisco_h5_path" in text


def test_figure_s6_runs_recovered_notebook_with_selected_paths():
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
    assert ENRICHMENT in text
    assert CLUSTER_MAP in text
    assert MODISCO_H5 in text
    assert "modisco_compiled_profile" in text
    assert "patch_source_notebook" in text
    assert "'source_logic_changed': False" in text
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
