"""Regression tests for Figure S2 mixed hiPSC–K562 panels."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
import yaml  # noqa: E402

from cc_perturbseq.figure_s02_parse_mix import (  # noqa: E402
    CAPTURE_PRIMER_ORDER,
    CELL_TYPE_ORDER,
    PRIMER_ORDER,
    plot_panel_i,
    plot_panel_j,
    plot_panel_l,
    plot_panel_m,
    plot_panel_n,
    preflight_parse_mix,
    read_saturation,
    summarize_guide_pairs,
)


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "supplementary" / "Figure_S02.ipynb"
SCOPE_CONFIG = ROOT / "config" / "figure_s02_parse_mix.yaml"


def test_scope_is_only_mixed_ipsc_k562_panels():
    config = yaml.safe_load(SCOPE_CONFIG.read_text(encoding="utf-8"))
    assert config["panels"] == ["i", "j", "l", "m", "n"]
    assert "250730_eTN10_Parsev3_mixediPSC_K562" in config["experiment_root"]
    assert "k" in config["excluded_panels"]

    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    source = "\n".join(
        "".join(cell.get("source", [])) for cell in notebook["cells"]
    )
    assert "panels i, j, l, m, and n only" in source
    assert '"panel": "i"' in source
    assert '"panel": "j"' in source
    assert '"panel": "l"' in source
    assert '"panel": "m"' in source
    assert '"panel": "n"' in source
    assert '"panel": "k"' not in source
    assert "hipsc_pilot_anndata" not in source
    assert "telohaec_pilot_anndata" not in source
    assert "standardize_parse_metrics" not in source


def test_notebook_code_cells_compile():
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    for index, cell in enumerate(notebook["cells"]):
        if cell.get("cell_type") == "code":
            source = cell.get("source", "")
            if isinstance(source, list):
                source = "".join(source)
            compile(source, f"{NOTEBOOK.name}:cell-{index}", "exec")


def _write_nonempty(path: Path, text: str = "x\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def test_preflight_accepts_source_directory_and_required_files(tmp_path):
    root = tmp_path / "250730_eTN10_Parsev3_mixediPSC_K562"
    expression_dir = (
        root
        / "splitpipe_out"
        / "WT"
        / "sub1"
        / "all-sample"
        / "DGE_filtered"
    )
    expression_dir.mkdir(parents=True)
    _write_nonempty(expression_dir / "count_matrix.mtx")
    _write_nonempty(expression_dir / "all_genes.csv")
    _write_nonempty(expression_dir / "cell_metadata.csv")
    _write_nonempty(
        root / "splitpipe_out" / "CRISPR_combined" / "grnaxcell_count.tsv"
    )
    _write_nonempty(
        root
        / "splitpipe_out"
        / "CRISPR"
        / "sub1"
        / "all-sample"
        / "guide_RNAs_unfiltered"
        / "all_guides.csv"
    )
    for primer in PRIMER_ORDER:
        _write_nonempty(
            root
            / "splitpipe_out"
            / "WT"
            / "sub1"
            / primer
            / "report"
            / "tscp_counts_subsampled.csv"
        )

    paths = preflight_parse_mix(root)
    assert paths["root"] == root
    assert paths["root"].is_dir()
    assert paths["expression_dir"] == expression_dir
    assert paths["expression_dir"].is_dir()
    assert paths["count_matrix"].is_file()
    assert paths["guide_counts"].is_file()
    assert all(
        paths[f"saturation_{primer}"].is_file() for primer in PRIMER_ORDER
    )


def test_preflight_reports_actual_missing_file_not_valid_directory(tmp_path):
    root = tmp_path / "experiment"
    expression_dir = (
        root
        / "splitpipe_out"
        / "WT"
        / "sub1"
        / "all-sample"
        / "DGE_filtered"
    )
    expression_dir.mkdir(parents=True)

    with pytest.raises(FileNotFoundError) as exc_info:
        preflight_parse_mix(root)

    message = str(exc_info.value)
    assert f"missing directory: {expression_dir}" not in message
    assert f"missing file: {expression_dir / 'count_matrix.mtx'}" in message


def test_preflight_rejects_file_where_expression_directory_is_expected(tmp_path):
    root = tmp_path / "experiment"
    expression_dir = (
        root
        / "splitpipe_out"
        / "WT"
        / "sub1"
        / "all-sample"
        / "DGE_filtered"
    )
    expression_dir.parent.mkdir(parents=True)
    expression_dir.write_text("not a directory", encoding="utf-8")

    with pytest.raises(FileNotFoundError) as exc_info:
        preflight_parse_mix(root)

    assert (
        f"expected directory but found another path type: {expression_dir}"
        in str(exc_info.value)
    )


def test_guide_summaries_match_panel_definitions():
    pairs = pd.DataFrame(
        {
            "cell": ["a__s1", "a__s1", "b", "b", "c"],
            "guide": [
                "gRNA_K562_g1",
                "gRNA_iPSC_g1",
                "iPSC_g2",
                "K562_g2",
                "K562_g3",
            ],
            "count": [5, 2, 4, 4, 2],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample": ["CS1", "DC1", "Dual"],
            "cell_type": ["K562", "hiPSC", "unassigned"],
        },
        index=["a", "b", "c"],
    )
    observed = summarize_guide_pairs(pairs, metadata)

    assert observed.loc["a", "total_guide_umis"] == 7
    assert observed.loc["a", "top_guide_umis"] == 5
    assert observed.loc["a", "matched_pool_guide_umis"] == 5
    assert np.isclose(observed.loc["a", "guide_pool_specificity"], 5 / 7)
    # Panel n removes the 2-UMI iPSC cell-guide observation before pool sums.
    assert observed.loc["a", "K562_guide_umis_ge3"] == 5
    assert observed.loc["a", "hiPSC_guide_umis_ge3"] == 0

    assert observed.loc["b", "matched_pool_guide_umis"] == 4
    assert np.isclose(observed.loc["b", "guide_pool_specificity"], 0.5)
    assert np.isnan(observed.loc["c", "matched_pool_guide_umis"])


def test_saturation_reader_uses_all_four_source_reports(tmp_path):
    for primer in PRIMER_ORDER:
        path = (
            tmp_path
            / "splitpipe_out"
            / "WT"
            / "sub1"
            / primer
            / "report"
            / "tscp_counts_subsampled.csv"
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(",hg38\n10,5\n20,9\n", encoding="utf-8")

    result = read_saturation(tmp_path)
    assert len(result) == 8
    assert set(result["primer"].astype(str)) == set(PRIMER_ORDER)
    assert list(CAPTURE_PRIMER_ORDER) == ["CS1", "DC1", "Dual"]
    assert result["mean_assigned_reads_per_cell"].notna().all()
    assert result["median_transcript_umis_per_cell"].notna().all()


def _synthetic_cell_metrics() -> pd.DataFrame:
    rows = []
    index = []
    counter = 0
    for cell_type in CELL_TYPE_ORDER:
        for sample_index, sample in enumerate(PRIMER_ORDER):
            for replicate in range(6):
                counter += 1
                matched = float(20 + sample_index * 5 + replicate)
                if cell_type == "K562":
                    k562_ge3 = matched
                    hipsc_ge3 = float(replicate % 3)
                else:
                    hipsc_ge3 = matched
                    k562_ge3 = float(replicate % 3)
                rows.append(
                    {
                        "sample": sample,
                        "cell_type": cell_type,
                        "matched_pool_guide_umis": matched,
                        "guide_pool_specificity": matched / (matched + 2.0),
                        "top_guide_umis": matched / 2.0,
                        "total_counts": 1000.0 + sample_index * 100 + replicate,
                        "hiPSC_guide_umis_ge3": hipsc_ge3,
                        "K562_guide_umis_ge3": k562_ge3,
                    }
                )
                index.append(f"cell_{counter}")
    return pd.DataFrame(rows, index=index)


def test_all_s2_plot_functions_smoke_without_seaborn_api_assumptions():
    metrics = _synthetic_cell_metrics()
    saturation = pd.DataFrame(
        [
            {
                "primer": primer,
                "mean_assigned_reads_per_cell": read_count,
                "median_transcript_umis_per_cell": 100 + 0.5 * read_count,
            }
            for primer in PRIMER_ORDER
            for read_count in (10, 20, 30, 40)
        ]
    )
    figures = []
    try:
        for producer in (plot_panel_i, plot_panel_j, plot_panel_l):
            figure, _ = producer(metrics)
            figures.append(figure)
        figure_m, _, _ = plot_panel_m(metrics, saturation)
        figures.append(figure_m)
        figure_n, source_n = plot_panel_n(metrics)
        figures.append(figure_n)
        assert not source_n.empty
        assert all(figure.axes for figure in figures)
    finally:
        for figure in figures:
            plt.close(figure)
