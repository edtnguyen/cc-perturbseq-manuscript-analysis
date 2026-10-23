from __future__ import annotations

from collections import OrderedDict
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from plotnine import (
    aes,
    coord_fixed,
    element_rect,
    element_text,
    facet_wrap,
    geom_point,
    ggplot,
    scale_color_cmap,
    theme,
    theme_void,
)

from .figure3_panels import (
    _add_horizontal_colorbar,
    _add_logo_strip,
    _draw_rect_heatmap,
    _infer_column,
    _load_logo_objects,
    _normalize_token,
    _parse_program,
    _program_row_lookup,
    _read_spectra,
    _split_features,
)
from .figure2_common import significance_stars
from .program_io import clean_day
from .single_cell_panels import expression_vector, program_vector


FIGURE4C_D12_CONDITIONS = (
    "D1_MPS",
    "D1.33_LM",
    "D1.67_LM",
    "D2_LM",
    "D1.33_LM+IAG",
    "D1.67_LM+IAG",
    "D2_LM+IAG",
    "D1.33_LM+TDI",
    "D1.67_LM+TDI",
    "D2_LM+TDI",
)
FIGURE4C_CONDITION_LABELS = {
    "D1_MPS": "D1 Control",
    "D1.33_LM": "D1.33 Control",
    "D1.67_LM": "D1.67 Control",
    "D2_LM": "D2 Control",
    "D1.33_LM+IAG": "D1-1.33 HIPPOa",
    "D1.67_LM+IAG": "D1-1.67 HIPPOa",
    "D2_LM+IAG": "D1-2 HIPPOa",
    "D1.33_LM+TDI": "D1-1.33 HIPPOi",
    "D1.67_LM+TDI": "D1-1.67 HIPPOi",
    "D2_LM+TDI": "D1-2 HIPPOi",
}
FIGURE4C_CONDITION_ORDER = tuple(
    FIGURE4C_CONDITION_LABELS[condition]
    for condition in FIGURE4C_D12_CONDITIONS
)
FIGURE4C_GENES = OrderedDict(
    {
        "Vascular": ("HAND1", "ETV2", "TAL1", "FLI1"),
        "Mesenchymal": ("ANKRD1", "KRT19", "SLC8A1"),
    }
)
FIGURE4C_PROGRAMS = ("P15", "P17")

FIGURE4D_FEATURES = (
    ("program", 15, "Programs\nP15\nVascular progenitor"),
    ("gene", "FLI1", "Markers\nFLI1"),
    ("program", 17, "Programs\nP17\nLateral mesoderm, YAP/TAZ-induced"),
    ("gene", "SLC8A1", "Markers\nSLC8A1"),
)

FIGURE4E_PROGRAMS = (15, 17)
FIGURE4E_MARKERS = OrderedDict(
    {
        15: ("ETV2", "TAL1", "FLI1"),
        17: ("SLC8A1", "KRT19", "ANKRD1", "HAND1"),
    }
)
FIGURE4E_LABELS = {
    15: "P15  Vascular progenitor",
    17: "P17  Lateral mesoderm, YAP/TAZ-induced",
}


def prepare_figure4c_dotplots(
    query_h5ad_path: str | Path,
    reprojected_h5ad_path: str | Path,
):
    """Prepare Figure 4c from the recovered D1-to-D2 Hippo workflow.

    The function repeats the source notebook's treatment selection, condition
    ordering, condition-specific QC flag, removal of sublibrary s8, library-size
    normalization plus log1p for marker genes, and the top-specific-200 CLR
    projection for P15/P17.
    """
    import anndata as ad
    import scanpy as sc
    from scipy import sparse

    query = sc.read_h5ad(query_h5ad_path)
    projected = sc.read_h5ad(reprojected_h5ad_path)

    projected_ids = pd.Index(projected.obs_names.astype(str))
    query_ids = pd.Index(query.obs_names.astype(str))
    common = projected_ids[projected_ids.isin(query_ids)]
    if not len(common):
        raise ValueError("Query and reprojected AnnData objects share no cells.")
    projected = projected[common].copy()
    query = query[common].copy()

    obs = projected.obs.copy()
    if "condition" not in obs.columns:
        if "sample_name" not in obs.columns:
            raise KeyError("Projected AnnData lacks both condition and sample_name.")
        obs["condition"] = obs["sample_name"].astype(str).str.replace(
            r"_rep\d+$", "", regex=True
        )

    keep = obs["condition"].astype(str).isin(FIGURE4C_D12_CONDITIONS).to_numpy()
    if "qc_pass_condition_specific" in obs.columns:
        keep &= obs["qc_pass_condition_specific"].fillna(False).astype(bool).to_numpy()
    if "sublibID" in obs.columns:
        keep &= ~obs["sublibID"].astype(str).eq("s8").to_numpy()
    if not keep.any():
        raise ValueError("No D1-to-D2 Hippo time-course cells passed the requested filters.")

    projected = projected[keep].copy()
    query = query[keep].copy()
    obs = projected.obs.copy()
    obs["condition_name"] = pd.Categorical(
        obs["condition"].astype(str).map(FIGURE4C_CONDITION_LABELS),
        categories=list(FIGURE4C_CONDITION_ORDER),
        ordered=True,
    )
    if obs["condition_name"].isna().any():
        missing = sorted(
            obs.loc[obs["condition_name"].isna(), "condition"].astype(str).unique()
        )
        raise ValueError(f"Unmapped Figure 4c conditions: {missing}")
    projected.obs["condition_name"] = obs["condition_name"]
    query.obs["condition_name"] = obs["condition_name"]

    if "gene" in query.var.columns:
        query.var_names = query.var["gene"].astype(str)
        query.var_names_make_unique()
    all_genes = [gene for group in FIGURE4C_GENES.values() for gene in group]
    missing_genes = [gene for gene in all_genes if gene not in query.var_names]
    if missing_genes:
        raise KeyError(f"Query AnnData is missing Figure 4c marker genes: {missing_genes}")

    gene_adata = query.copy()
    sc.pp.normalize_total(gene_adata)
    sc.pp.log1p(gene_adata)
    gene_adata = gene_adata[:, all_genes].copy()

    clr_key = "X_cnmf_clr_top_specific_200"
    if clr_key not in projected.obsm:
        raise KeyError(
            f"Reprojected AnnData lacks obsm[{clr_key!r}]. "
            f"Available keys: {list(projected.obsm.keys())}"
        )
    program_names = pd.Index(projected.var_names.astype(str))
    missing_programs = [
        program for program in FIGURE4C_PROGRAMS if program not in program_names
    ]
    if missing_programs:
        raise KeyError(f"Reprojected AnnData is missing programs: {missing_programs}")
    program_positions = [
        int(program_names.get_loc(program)) for program in FIGURE4C_PROGRAMS
    ]
    clr = np.asarray(projected.obsm[clr_key], dtype=float)[:, program_positions]
    program_adata = ad.AnnData(
        X=clr.copy(),
        obs=obs.copy(),
        var=pd.DataFrame(index=pd.Index(FIGURE4C_PROGRAMS, name="program")),
    )
    program_adata.layers["clr_top_specific_200"] = clr.copy()

    gene_dotplot = sc.pl.dotplot(
        gene_adata,
        var_names={key: list(value) for key, value in FIGURE4C_GENES.items()},
        groupby="condition_name",
        standard_scale="var",
        cmap="Reds",
        colorbar_title="Normalized expression",
        size_title="Fraction of cells\nexpressing (%)",
        return_fig=True,
        show=False,
    )
    program_dotplot = sc.pl.dotplot(
        program_adata,
        groupby="condition_name",
        var_names=list(FIGURE4C_PROGRAMS),
        layer="clr_top_specific_200",
        cmap="coolwarm",
        colorbar_title="Mean CLR value",
        size_title="Fraction of cells\nwith CLR>0 (%)",
        vcenter=0,
        return_fig=True,
        show=False,
    )

    gene_matrix = gene_adata.X
    if sparse.issparse(gene_matrix):
        gene_matrix = gene_matrix.toarray()
    gene_matrix = np.asarray(gene_matrix, dtype=float)
    gene_rows: list[dict[str, Any]] = []
    for condition in FIGURE4C_CONDITION_ORDER:
        mask = gene_adata.obs["condition_name"].astype(str).to_numpy() == condition
        for index, gene in enumerate(all_genes):
            values = gene_matrix[mask, index]
            gene_rows.append(
                {
                    "condition": condition,
                    "gene": gene,
                    "n_cells": int(mask.sum()),
                    "mean_log_normalized_expression": float(np.mean(values)),
                    "fraction_expressing": float(np.mean(values > 0)),
                }
            )

    program_rows: list[dict[str, Any]] = []
    for condition in FIGURE4C_CONDITION_ORDER:
        mask = program_adata.obs["condition_name"].astype(str).to_numpy() == condition
        for index, program in enumerate(FIGURE4C_PROGRAMS):
            values = clr[mask, index]
            program_rows.append(
                {
                    "condition": condition,
                    "program": program,
                    "n_cells": int(mask.sum()),
                    "mean_clr": float(np.mean(values)),
                    "fraction_clr_positive": float(np.mean(values > 0)),
                }
            )

    metadata = {
        "n_common_cells": len(common),
        "n_D1_to_D2_QC_cells": int(keep.sum()),
        "conditions": list(FIGURE4C_CONDITION_ORDER),
        "program_projection": clr_key,
        "program_projection_gene_set": "top_specific_200",
        "gene_standard_scale": "var",
        "excluded_sublibrary": "s8 if present",
        "source_workflow": "hippo_tc_analysis_clean.ipynb D1-2 subset",
    }
    return (
        gene_dotplot,
        program_dotplot,
        pd.DataFrame(gene_rows),
        pd.DataFrame(program_rows),
        metadata,
    )


def plot_figure4d_feature_umaps(
    adata,
    usage: pd.DataFrame,
    *,
    point_size: float = 0.42,
    background_size: float = 0.18,
    vmax_quantile: float = 0.99,
):
    """Render P15/P17 and FLI1/SLC8A1 using the recovered plotnine UMAP style."""
    if "X_umap" not in adata.obsm:
        raise KeyError("AnnData lacks obsm['X_umap']")
    coords = np.asarray(adata.obsm["X_umap"])
    feature_order = [label for _, _, label in FIGURE4D_FEATURES]
    background_frames: list[pd.DataFrame] = []
    foreground_frames: list[pd.DataFrame] = []
    source: dict[str, pd.DataFrame] = {}

    for kind, feature, label in FIGURE4D_FEATURES:
        values = (
            program_vector(usage, feature)
            if kind == "program"
            else expression_vector(adata, str(feature))
        )
        values = np.asarray(values, dtype=float)
        finite = values[np.isfinite(values)]
        vmax = float(np.quantile(finite, vmax_quantile)) if finite.size else 1.0
        vmax = max(vmax, np.nextafter(0, 1))
        display = np.clip(np.nan_to_num(values, nan=0.0) / vmax, 0, 1)
        key = f"{kind}_{feature}"
        source[key] = pd.DataFrame(
            {
                "umap_1": coords[:, 0],
                "umap_2": coords[:, 1],
                "value": values,
                "display_vmax_p99": vmax,
                "display_fraction": display,
            },
            index=adata.obs_names,
        )
        background_frames.append(
            pd.DataFrame(
                {
                    "umap_1": coords[:, 0],
                    "umap_2": coords[:, 1],
                    "feature": label,
                }
            )
        )
        positive = np.flatnonzero(np.isfinite(values) & (values > 0))
        positive = positive[np.argsort(values[positive])]
        foreground_frames.append(
            pd.DataFrame(
                {
                    "umap_1": coords[positive, 0],
                    "umap_2": coords[positive, 1],
                    "display_fraction": display[positive],
                    "feature": label,
                }
            )
        )

    background = pd.concat(background_frames, ignore_index=True)
    foreground = pd.concat(foreground_frames, ignore_index=True)
    background["feature"] = pd.Categorical(
        background["feature"], feature_order, ordered=True
    )
    foreground["feature"] = pd.Categorical(
        foreground["feature"], feature_order, ordered=True
    )

    plot = (
        ggplot()
        + geom_point(
            background,
            aes("umap_1", "umap_2"),
            inherit_aes=False,
            color="#D9D9D9",
            size=background_size,
            stroke=0,
            raster=True,
        )
        + geom_point(
            foreground,
            aes("umap_1", "umap_2", color="display_fraction"),
            inherit_aes=False,
            size=point_size,
            stroke=0,
            raster=True,
        )
        + scale_color_cmap(cmap_name="YlGnBu", limits=(0, 1), guide=None)
        + facet_wrap("~feature", ncol=2)
        + coord_fixed()
        + theme_void()
        + theme(
            figure_size=(5.7, 5.2),
            strip_text=element_text(
                family="Liberation Sans", size=8.2, color="black"
            ),
            strip_background=element_rect(fill="white", color="white"),
            panel_spacing=0.10,
            plot_background=element_rect(fill="white", color="white"),
            panel_background=element_rect(fill="white", color="white"),
            plot_margin=0.02,
        )
    )
    return plot, source


def _choose_spectra_offset_figure4(spectra: pd.DataFrame) -> int:
    row_lookup = _program_row_lookup(spectra)
    gene_lookup = {_normalize_token(gene): gene for gene in spectra.columns}
    scores: dict[int, float] = {}
    for offset in (0, 1, -1):
        values: list[float] = []
        for program, genes in FIGURE4E_MARKERS.items():
            raw_program = program + offset
            if raw_program not in row_lookup:
                continue
            row = row_lookup[raw_program]
            for gene in genes:
                matrix_gene = gene_lookup.get(_normalize_token(gene))
                if matrix_gene is None:
                    continue
                value = spectra.loc[row, matrix_gene]
                if pd.notna(value):
                    values.append(float(value))
        scores[offset] = float(np.mean(values)) if values else -np.inf
    return max(scores, key=scores.get)


def _load_figure4e_regulator_tables(
    workbook_path: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str], pd.DataFrame]:
    workbook_path = Path(workbook_path)
    summary = pd.read_excel(workbook_path, sheet_name="Summary")
    summary.columns = [str(column).strip() for column in summary.columns]
    summary_program_col = next(
        (
            column
            for column in summary.columns
            if column.lower() in {"program_name", "program", "program number"}
            or "program" in column.lower()
        ),
        summary.columns[0],
    )
    summary["program_num"] = summary[summary_program_col].map(_parse_program)
    regulator_column = "Top 5 specific regulators (FDR<0.1) D2"
    if regulator_column not in summary.columns:
        raise KeyError(
            f"Workbook Summary is missing {regulator_column!r}; "
            f"available columns: {summary.columns.tolist()}"
        )

    selected_by_program: OrderedDict[int, list[str]] = OrderedDict()
    selection_rows: list[dict[str, Any]] = []
    for program in FIGURE4E_PROGRAMS:
        matched = summary.loc[summary["program_num"].eq(program)]
        if matched.empty:
            raise KeyError(f"Summary has no row for P{program}.")
        regulators = _split_features(matched.iloc[0][regulator_column])[:3]
        selected_by_program[program] = regulators
        for rank, gene in enumerate(regulators, start=1):
            selection_rows.append(
                {
                    "program": program,
                    "rank_within_program": rank,
                    "gene": gene,
                    "selection_source": regulator_column,
                }
            )

    regulator_union: list[str] = []
    seen: set[str] = set()
    for program in FIGURE4E_PROGRAMS:
        for gene in selected_by_program[program]:
            key = _normalize_token(gene)
            if key and key not in seen:
                regulator_union.append(gene)
                seen.add(key)

    pert = pd.read_excel(workbook_path, sheet_name="Perturbation Association 1")
    pert.columns = [str(column).strip() for column in pert.columns]
    rename = {
        "target_name": "gene",
        "program_name": "program",
        "Sample": "day",
        "approx_log2FC": "effect",
        "adj_pval": "fdr",
        "Perturbed gene": "gene",
        "Gene Program": "program",
        "adj.P.Val": "fdr",
    }
    pert = pert.rename(
        columns={key: value for key, value in rename.items() if key in pert.columns}
    )
    required = {"gene", "program", "day", "effect", "fdr"}
    missing = required - set(pert.columns)
    if missing:
        raise ValueError(
            "Perturbation Association 1 is missing columns after rename: "
            f"{sorted(missing)}"
        )
    pert["program_num"] = pert["program"].map(_parse_program)
    pert["day_clean"] = pert["day"].map(clean_day)
    pert["gene_key"] = pert["gene"].map(_normalize_token)
    pert["effect"] = pd.to_numeric(pert["effect"], errors="coerce")
    pert["fdr"] = pd.to_numeric(pert["fdr"], errors="coerce")
    d2 = pert.loc[
        pert["day_clean"].eq("D2")
        & pert["program_num"].isin(FIGURE4E_PROGRAMS)
    ].copy()

    effects = pd.DataFrame(
        index=FIGURE4E_PROGRAMS, columns=regulator_union, dtype=float
    )
    fdrs = effects.copy()
    stars = pd.DataFrame(
        "", index=FIGURE4E_PROGRAMS, columns=regulator_union, dtype=object
    )
    for program in FIGURE4E_PROGRAMS:
        for gene in regulator_union:
            subset = d2.loc[
                d2["program_num"].eq(program)
                & d2["gene_key"].eq(_normalize_token(gene))
            ].sort_values("fdr")
            if subset.empty:
                continue
            effects.loc[program, gene] = subset.iloc[0]["effect"]
            fdrs.loc[program, gene] = subset.iloc[0]["fdr"]
            stars.loc[program, gene] = significance_stars(subset.iloc[0]["fdr"])
    return effects, fdrs, stars, regulator_union, pd.DataFrame(selection_rows)


def _load_figure4e_marker_table(
    spectra_path: str | Path,
) -> tuple[pd.DataFrame, int, pd.DataFrame]:
    spectra = _read_spectra(spectra_path)
    offset = _choose_spectra_offset_figure4(spectra)
    row_lookup = _program_row_lookup(spectra)
    gene_lookup = {_normalize_token(gene): gene for gene in spectra.columns}
    marker_union: list[str] = []
    for program in FIGURE4E_PROGRAMS:
        for gene in FIGURE4E_MARKERS[program]:
            if gene not in marker_union:
                marker_union.append(gene)

    matrix = pd.DataFrame(index=FIGURE4E_PROGRAMS, columns=marker_union, dtype=float)
    source_rows: list[dict[str, Any]] = []
    for program in FIGURE4E_PROGRAMS:
        raw_program = program + offset
        if raw_program not in row_lookup:
            raise KeyError(f"No consensus row found for P{program} with offset {offset}.")
        row_label = row_lookup[raw_program]
        for rank, gene in enumerate(FIGURE4E_MARKERS[program], start=1):
            matrix_gene = gene_lookup.get(_normalize_token(gene))
            value = (
                spectra.loc[row_label, matrix_gene]
                if matrix_gene is not None
                else np.nan
            )
            matrix.loc[program, gene] = value
            source_rows.append(
                {
                    "program": program,
                    "rank_within_program": rank,
                    "gene": gene,
                    "consensus_program_offset": offset,
                    "raw_consensus_row": row_label,
                    "matrix_gene": matrix_gene,
                    "loading": value,
                }
            )
    return matrix, offset, pd.DataFrame(source_rows)


def _load_figure4e_motif_tables(
    motif_path: str | Path,
    *,
    top_n: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str], pd.DataFrame]:
    table = pd.read_csv(motif_path, sep="\t")
    program_col = _infer_column(table, ["program", "program_num", "cnmf"])
    day_col = _infer_column(table, ["day", "timepoint"])
    motif_col = _infer_column(table, ["tf", "motif", "tf_motif_id"])
    q_col = _infer_column(table, ["qval", "q_value", "qvalue", "fdr", "padj"])
    neglog_col = _infer_column(
        table,
        ["neg_log10_qval", "minus_log10_q", "-log10"],
        required=False,
    )
    auc_col = _infer_column(table, ["AUC", "auc"], required=False)

    data = pd.DataFrame(
        {
            "program": table[program_col].map(_parse_program),
            "day": table[day_col].map(clean_day),
            "motif": table[motif_col].astype(str).str.strip(),
            "q_value": pd.to_numeric(table[q_col], errors="coerce"),
            "auc": (
                pd.to_numeric(table[auc_col], errors="coerce")
                if auc_col
                else np.nan
            ),
        }
    )
    data["motif_key"] = data["motif"].map(_normalize_token)
    if neglog_col:
        data["neg_log10_q"] = pd.to_numeric(table[neglog_col], errors="coerce")
    else:
        data["neg_log10_q"] = -np.log10(
            data["q_value"].clip(lower=np.nextafter(0, 1))
        )
    data = data.loc[
        data["day"].eq("D2")
        & data["program"].isin(FIGURE4E_PROGRAMS)
        & data["q_value"].notna()
        & data["motif_key"].notna()
    ].copy()

    selected_by_program: OrderedDict[int, list[str]] = OrderedDict()
    selection_rows: list[dict[str, Any]] = []
    for program in FIGURE4E_PROGRAMS:
        subset = data.loc[data["program"].eq(program)].copy()
        subset = (
            subset.sort_values(
                ["q_value", "neg_log10_q", "auc"],
                ascending=[True, False, False],
                na_position="last",
            )
            .drop_duplicates("motif_key", keep="first")
            .head(top_n)
        )
        motifs = subset["motif"].tolist()
        if len(motifs) < top_n:
            raise ValueError(
                f"Only {len(motifs)} motif rows were available for D2 P{program}; "
                f"{top_n} are required."
            )
        selected_by_program[program] = motifs
        for rank, (_, row) in enumerate(subset.iterrows(), start=1):
            selection_rows.append(
                {
                    "program": program,
                    "rank_within_program": rank,
                    "motif": row["motif"],
                    "q_value": row["q_value"],
                    "neg_log10_q": row["neg_log10_q"],
                    "auc": row["auc"],
                    "significance": significance_stars(row["q_value"]),
                    "selection_rule": (
                        f"top {top_n} motifs by q-value regardless of significance"
                    ),
                }
            )

    motif_union: list[str] = []
    seen: set[str] = set()
    for program in FIGURE4E_PROGRAMS:
        for motif in selected_by_program[program]:
            key = _normalize_token(motif)
            if key and key not in seen:
                motif_union.append(motif)
                seen.add(key)

    values = pd.DataFrame(index=FIGURE4E_PROGRAMS, columns=motif_union, dtype=float)
    qvalues = values.copy()
    stars = pd.DataFrame(
        "", index=FIGURE4E_PROGRAMS, columns=motif_union, dtype=object
    )
    for program in FIGURE4E_PROGRAMS:
        for motif in motif_union:
            subset = data.loc[
                data["program"].eq(program)
                & data["motif_key"].eq(_normalize_token(motif))
            ].sort_values(["q_value", "neg_log10_q"], ascending=[True, False])
            if subset.empty:
                continue
            row = subset.iloc[0]
            values.loc[program, motif] = row["neg_log10_q"]
            qvalues.loc[program, motif] = row["q_value"]
            stars.loc[program, motif] = significance_stars(row["q_value"])
    return values, qvalues, stars, motif_union, pd.DataFrame(selection_rows)


def make_figure4e_annotation(
    *,
    workbook_path: str | Path,
    spectra_path: str | Path,
    motif_path: str | Path,
    cluster_map_path: str | Path,
    modisco_h5_path: str | Path,
    top_motifs_per_program: int = 5,
):
    """Reproduce the focused D2 P15/P17 annotation workflow."""
    effects, effect_fdrs, effect_stars, regulators, regulator_selection = (
        _load_figure4e_regulator_tables(workbook_path)
    )
    markers, spectra_offset, marker_selection = _load_figure4e_marker_table(
        spectra_path
    )
    (
        motif_values,
        motif_qvalues,
        motif_stars,
        motifs,
        motif_selection,
    ) = _load_figure4e_motif_tables(
        motif_path, top_n=top_motifs_per_program
    )
    logo_objects, logo_map = _load_logo_objects(
        motifs, cluster_map_path, modisco_h5_path
    )

    regulator_values = np.abs(effects.to_numpy(float))
    regulator_values = regulator_values[np.isfinite(regulator_values)]
    regulator_vmax = float(np.max(regulator_values)) if regulator_values.size else 1.0
    marker_values = markers.to_numpy(float)
    marker_values = marker_values[np.isfinite(marker_values)]
    marker_vmax = float(np.max(marker_values)) if marker_values.size else 1.0
    motif_finite = motif_values.to_numpy(float)
    motif_finite = motif_finite[np.isfinite(motif_finite)]
    motif_raw_vmax = float(np.max(motif_finite)) if motif_finite.size else 1.0
    motif_vmax = (
        float(np.percentile(motif_finite, 95)) if motif_finite.size else 1.0
    )
    regulator_vmax = max(regulator_vmax, 1e-12)
    marker_vmax = max(marker_vmax, 1e-12)
    motif_vmax = max(motif_vmax, 1e-12)

    regulator_norm = TwoSlopeNorm(
        vmin=-regulator_vmax, vcenter=0, vmax=regulator_vmax
    )
    marker_norm = Normalize(vmin=0, vmax=marker_vmax, clip=True)
    motif_norm = Normalize(vmin=0, vmax=motif_vmax, clip=True)
    regulator_cmap = plt.colormaps["RdBu_r"]
    marker_cmap = plt.colormaps["Greys"]
    motif_cmap = plt.colormaps["Purples"]

    widths = [
        max(2.2, 0.34 * len(regulators)),
        max(2.2, 0.34 * len(markers.columns)),
        max(2.2, 0.34 * len(motifs)),
    ]
    left_margin = 3.25
    right_margin = 0.65
    gap = 0.62
    bottom_margin = 1.85
    heatmap_height = 0.72
    logo_region = 1.20
    header_region = 1.10
    top_margin = 0.35
    figure_width = left_margin + sum(widths) + 2 * gap + right_margin
    figure_height = (
        bottom_margin + heatmap_height + logo_region + header_region + top_margin
    )
    figure = plt.figure(figsize=(figure_width, figure_height), facecolor="white")
    fx = lambda value: value / figure_width
    fy = lambda value: value / figure_height
    lefts = [
        left_margin,
        left_margin + widths[0] + gap,
        left_margin + widths[0] + gap + widths[1] + gap,
    ]
    axes = [
        figure.add_axes(
            [
                fx(lefts[index]),
                fy(bottom_margin),
                fx(widths[index]),
                fy(heatmap_height),
            ]
        )
        for index in range(3)
    ]
    labels = [FIGURE4E_LABELS[program] for program in FIGURE4E_PROGRAMS]
    _draw_rect_heatmap(
        axes[0],
        effects.loc[list(FIGURE4E_PROGRAMS), regulators],
        cmap=regulator_cmap,
        norm=regulator_norm,
        xlabels=regulators,
        ylabels=labels,
        stars=effect_stars.loc[list(FIGURE4E_PROGRAMS), regulators],
        show_ylabels=True,
    )
    _draw_rect_heatmap(
        axes[1],
        markers.loc[list(FIGURE4E_PROGRAMS)],
        cmap=marker_cmap,
        norm=marker_norm,
        xlabels=list(markers.columns),
        ylabels=labels,
        stars=None,
        show_ylabels=False,
        italic_x=True,
    )
    _draw_rect_heatmap(
        axes[2],
        motif_values.loc[list(FIGURE4E_PROGRAMS), motifs],
        cmap=motif_cmap,
        norm=motif_norm,
        xlabels=motifs,
        ylabels=labels,
        stars=motif_stars.loc[list(FIGURE4E_PROGRAMS), motifs],
        show_ylabels=False,
    )
    _add_logo_strip(figure, axes[2], motifs, logo_objects)

    colorbar_bottom = fy(bottom_margin + heatmap_height + logo_region + 0.08)
    colorbar_height = fy(0.12)
    _add_horizontal_colorbar(
        figure,
        left=fx(lefts[0]),
        bottom=colorbar_bottom,
        width=fx(widths[0]),
        height=colorbar_height,
        cmap=regulator_cmap,
        norm=regulator_norm,
        label="approx. log$_2$FC",
    )
    _add_horizontal_colorbar(
        figure,
        left=fx(lefts[1]),
        bottom=colorbar_bottom,
        width=fx(widths[1]),
        height=colorbar_height,
        cmap=marker_cmap,
        norm=marker_norm,
        label="consensus cNMF loading",
    )
    _add_horizontal_colorbar(
        figure,
        left=fx(lefts[2]),
        bottom=colorbar_bottom,
        width=fx(widths[2]),
        height=colorbar_height,
        cmap=motif_cmap,
        norm=motif_norm,
        label="$-\\log_{10}$(q-value), color clipped at 95th percentile",
    )
    title_y = colorbar_bottom + colorbar_height + fy(0.23)
    titles = (
        "Effect of gene CRISPRi\non cNMF program expression",
        "cNMF gene loading",
        "TF motif enrichment",
    )
    for left, width, title in zip(lefts, widths, titles):
        figure.text(
            fx(left + width / 2),
            title_y,
            title,
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )
    figure.text(
        fx(0.20),
        fy(figure_height - 0.15),
        "e",
        ha="left",
        va="top",
        fontsize=16,
        fontweight="bold",
    )
    figure.text(
        fx(left_margin),
        fy(0.22),
        "* q/FDR < 0.05;  ** < 0.01;  *** < 0.001. "
        "The top five motifs per program are displayed regardless of significance.",
        ha="left",
        va="bottom",
        fontsize=7.5,
        color="#444444",
    )

    source_tables = {
        "regulator_effects": effects,
        "regulator_fdr": effect_fdrs,
        "regulator_stars": effect_stars,
        "regulator_selection": regulator_selection,
        "marker_loadings": markers,
        "marker_selection": marker_selection,
        "motif_values_neg_log10_q": motif_values,
        "motif_qvalues": motif_qvalues,
        "motif_stars": motif_stars,
        "motif_selection": motif_selection,
        "motif_logo_mapping": logo_map,
    }
    metadata = {
        "programs": list(FIGURE4E_PROGRAMS),
        "day": "D2",
        "top_regulators_per_program": 3,
        "top_motifs_per_program": top_motifs_per_program,
        "motif_selection_requires_significance": False,
        "consensus_program_offset": spectra_offset,
        "motif_raw_max": motif_raw_vmax,
        "motif_display_vmax_p95": motif_vmax,
        "motif_value": "-log10(BH q-value)",
    }
    return figure, source_tables, metadata
