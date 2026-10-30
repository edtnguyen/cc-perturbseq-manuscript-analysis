from __future__ import annotations

from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import cc_perturbseq.s4_source_runner as s4  # noqa: E402


def test_resolver_accepts_only_frozen_source(monkeypatch, tmp_path):
    source = tmp_path / "source.ipynb"
    source.write_text('{"cells": [], "nbformat": 4, "nbformat_minor": 5}')
    observed = s4.sha256_file(source)
    monkeypatch.setattr(s4, "OAK_SOURCE_SHA256", observed)
    resolved = s4.resolve_s4_source_notebook(tmp_path, source)
    assert resolved.path == str(source.resolve())
    assert resolved.source_kind == "executed_sherlock_source"
    assert resolved.verification_algorithm == "sha256"

    source.write_text("changed")
    with pytest.raises(ValueError, match="does not match either frozen prior version"):
        s4.resolve_s4_source_notebook(tmp_path, source)


def test_collects_exact_source_outputs_without_replotting(tmp_path):
    source_path = tmp_path / "bundle" / "source.ipynb"
    source_path.parent.mkdir(parents=True)
    source_path.write_text("{}")
    source = s4.S4SourceNotebook(
        path=str(source_path),
        source_kind="test",
        verification_algorithm="sha256",
        expected_hash="x",
        observed_hash="x",
    )
    model = source_path.parent / s4.MODEL_OUTPUT_RELATIVE
    qq = model / "illustrator_editable_full_covariates"
    qq.mkdir(parents=True)
    (qq / "d0_QQ_full_covariates.svg").write_text("svg")
    (model / "all_days_QQ_points.csv").write_text("p\n1\n")
    (model / "QQ_export_manifest.csv").write_text("day\nd0\n")

    plots = model / "biological_replicate_effects" / "plots"
    plots.mkdir(parents=True)
    (plots / "biological_replicate_DE_program_effect_correlation_all_days.svg").write_text("svg")
    (plots / "biological_replicate_DE_program_effect_pairs.csv").write_text("beta_x,beta_y\n1,1\n")
    (plots / "biological_replicate_DE_program_effect_correlation_summary.csv").write_text("pearson_r\n1\n")
    (plots.parent / "replicate_beta_manifest.csv").write_text("status\nwritten\n")

    collected = s4.collect_s4_source_outputs(
        source,
        panel_b_directory=tmp_path / "panel_b",
        panel_c_directory=tmp_path / "panel_c",
    )
    assert any(Path(p).suffix == ".svg" for p in collected["panel_b"].values())
    assert any(Path(p).suffix == ".svg" for p in collected["panel_c"].values())
    assert any(Path(p).name == "all_days_QQ_points.csv" for p in collected["panel_b"].values())
    assert any("effect_pairs" in Path(p).name for p in collected["panel_c"].values())
