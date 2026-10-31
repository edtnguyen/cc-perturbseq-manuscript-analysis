from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from .figure2_common import FIGURE2_PROGRAMS, FIGURE2_TARGET_GROUPS, sha256_file

def load_figure2_cnmf_workbook(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load the repository-local copy of the supplied cNMF workbook.

    The canonical input is the full workbook supplied by the author. For
    backwards compatibility, the function also accepts the compact plotting
    extract generated during repository curation.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)

    workbook = pd.ExcelFile(path)
    sheets = set(workbook.sheet_names)
    if {"Effect matrix", "Selected associations", "Provenance"}.issubset(sheets):
        matrix = pd.read_excel(workbook, sheet_name="Effect matrix", index_col=0)
        selected = pd.read_excel(workbook, sheet_name="Selected associations")
        provenance = pd.read_excel(workbook, sheet_name="Provenance")
    elif "Perturbation Association 1" in sheets:
        association = pd.read_excel(workbook, sheet_name="Perturbation Association 1")
        association = association.dropna(axis=0, how="all").dropna(axis=1, how="all")
        association = association.loc[:, ~association.columns.astype(str).str.startswith("Unnamed")].copy()
        required_source = {"target_name", "program_name", "Sample", "approx_log2FC", "adj_pval"}
        missing = required_source.difference(association.columns)
        if missing:
            raise ValueError(f"Perturbation Association 1 is missing: {sorted(missing)}")
        association["target_name"] = association["target_name"].astype(str).str.strip()
        association["Sample"] = association["Sample"].astype(str).str.strip().str.upper()
        association["program_name"] = pd.to_numeric(
            association["program_name"], errors="raise"
        ).astype(int)
        association["approx_log2FC"] = pd.to_numeric(
            association["approx_log2FC"], errors="raise"
        )
        association["adj_pval"] = pd.to_numeric(association["adj_pval"], errors="raise")

        # NEDD8-MDP1 is a duplicate target-name alias documented in the screen
        # handoff. Retain canonical NEDD8 and the pooled non-targeting profile.
        canonical = association.loc[association["target_name"].ne("NEDD8-MDP1")].copy()
        feature_order = pd.MultiIndex.from_product(
            [["D0", "D1", "D2", "D3"], range(50)],
            names=["Sample", "program_name"],
        )
        matrix = (
            canonical.pivot(
                index="target_name",
                columns=["Sample", "program_name"],
                values="approx_log2FC",
            )
            .reindex(columns=feature_order)
            .sort_index()
        )
        matrix.columns = [f"{day}_P{program}" for day, program in matrix.columns]

        group_by_target = {
            target: group
            for group, targets in FIGURE2_TARGET_GROUPS.items()
            for target in targets
        }
        label_by_key = {
            (item["day"], int(item["program"])): item["label"]
            for item in FIGURE2_PROGRAMS
        }
        selected_keys = set(label_by_key)
        selected_targets = set(group_by_target)
        selected = canonical.loc[
            canonical["target_name"].isin(selected_targets)
            & canonical.apply(
                lambda row: (str(row["Sample"]), int(row["program_name"])) in selected_keys,
                axis=1,
            )
        ].copy()
        selected.insert(
            0,
            "group",
            selected["target_name"].map(group_by_target),
        )
        selected["program_label"] = [
            label_by_key[(str(day), int(program))]
            for day, program in zip(selected["Sample"], selected["program_name"])
        ]
        provenance = pd.DataFrame(
            {
                "field": [
                    "source_filename",
                    "source_sha256",
                    "source_association_rows",
                    "source_targets",
                    "canonical_filter",
                ],
                "value": [
                    path.name,
                    sha256_file(path),
                    len(association),
                    association["target_name"].nunique(),
                    "drop duplicate alias NEDD8-MDP1; retain canonical NEDD8 and pooled non-targeting",
                ],
            }
        )
    else:
        raise ValueError(
            "Workbook must contain either 'Perturbation Association 1' or the "
            "compact Figure 2 plotting sheets."
        )

    matrix.index = matrix.index.astype(str)
    matrix = matrix.apply(pd.to_numeric, errors="raise")
    expected_features = [
        f"{day}_P{program}"
        for day in ("D0", "D1", "D2", "D3")
        for program in range(50)
    ]
    if matrix.shape != (298, 200):
        raise ValueError(f"Expected a 298 × 200 effect matrix; found {matrix.shape}.")
    if list(matrix.columns) != expected_features:
        raise ValueError("Effect-matrix columns do not match D0–D3 × P0–P49 order.")
    if matrix.isna().any().any():
        raise ValueError("Canonical effect matrix contains missing values.")
    if "non-targeting" not in matrix.index:
        raise ValueError("Effect matrix must retain the pooled non-targeting profile.")
    if "NEDD8-MDP1" in matrix.index:
        raise ValueError("Duplicate alias NEDD8-MDP1 must not be present in the canonical matrix.")

    required = {
        "group", "target_name", "Sample", "program_name", "program_label",
        "approx_log2FC", "adj_pval",
    }
    missing = required.difference(selected.columns)
    if missing:
        raise ValueError(f"Selected-association data are missing: {sorted(missing)}")
    selected = selected.copy()
    selected["target_name"] = selected["target_name"].astype(str)
    selected["Sample"] = selected["Sample"].astype(str)
    selected["program_name"] = pd.to_numeric(
        selected["program_name"], errors="raise"
    ).astype(int)
    selected["approx_log2FC"] = pd.to_numeric(
        selected["approx_log2FC"], errors="raise"
    )
    selected["adj_pval"] = pd.to_numeric(selected["adj_pval"], errors="raise")
    if len(selected) != 182:
        raise ValueError(f"Expected 182 selected target–program rows; found {len(selected)}.")
    return matrix, selected, provenance

def cluster_effect_matrix(matrix: pd.DataFrame) -> tuple[pd.DataFrame, list[str], np.ndarray]:
    values = matrix.to_numpy(float)
    correlation = np.corrcoef(values)
    correlation = np.clip(correlation, -1.0, 1.0)
    np.fill_diagonal(correlation, 1.0)
    distance = np.clip(1.0 - correlation, 0.0, 2.0)
    np.fill_diagonal(distance, 0.0)
    tree = linkage(squareform(distance, checks=False), method="average", optimal_ordering=True)
    order_idx = leaves_list(tree)[::-1]
    order = matrix.index.to_numpy()[order_idx].astype(str).tolist()
    corr_df = pd.DataFrame(correlation, index=matrix.index, columns=matrix.index)
    return corr_df.loc[order, order], order, order_idx
