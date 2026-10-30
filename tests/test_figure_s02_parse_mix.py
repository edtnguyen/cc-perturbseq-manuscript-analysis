"""Regression tests for Figure S2 mixed hiPSC–K562 panels."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from cc_perturbseq.figure_s02_parse_mix import (
    CAPTURE_PRIMER_ORDER,
    PRIMER_ORDER,
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
            compile(
                "".join(cell.get("source", [])),
                f"{NOTEBOOK.name}:cell-{index}",
                "exec",
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
