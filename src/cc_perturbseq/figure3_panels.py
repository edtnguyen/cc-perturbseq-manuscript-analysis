from __future__ import annotations

from collections import OrderedDict
from collections.abc import Sequence
from pathlib import Path
import hashlib
import re
from typing import Any

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.transforms import Affine2D, Bbox
from plotnine import (
    aes,
    annotate,
    coord_fixed,
    element_blank,
    element_rect,
    element_text,
    facet_wrap,
    geom_boxplot,
    geom_hline,
    geom_point,
    geom_text,
    ggplot,
    labs,
    scale_color_cmap,
    scale_fill_manual,
    scale_x_continuous,
    scale_y_continuous,
    theme,
    theme_classic,
    theme_void,
)

from .figure2_carrier_plotnine import (
    compute_figure2f_carrier_seeded,
    plot_figure2f,
)
from .figure2_common import bh_adjust, significance_stars
from .figure2_panel_d_120 import compute_figure2d_guide_ec_120_controls
from .program_io import clean_day
from .single_cell_panels import expression_vector, program_vector


FIGURE3B_FEATURES = (
    ("program", 2, "Programs\nP2\nArtery endothelium 1"),
    ("gene", "CDH5", "Markers\nCDH5"),
    ("program", 34, "Programs\nP34\nMesenchyme, YAP/TAZ-induced"),
    ("gene", "ACTC1", "Markers\nACTC1"),
)
FIGURE3_HIPPO_TARGETS = ("NF2", "LATS2", "AMOTL2", "TAOK1")
FIGURE3D_PROGRAMS = (2, 34)
FIGURE3D_MARKERS = OrderedDict(
    {
        2: ("CDH5", "NRP1", "EFNB2"),
        34: ("TAGLN", "KRT8", "ACTC1", "ANKRD1"),
    }
)
FIGURE3D_LABELS = {
    2: "P2  Artery endothelium 1",
    34: "P34  Mesenchyme, YAP/TAZ-induced",
}
FIGURE3F_D3_CONDITIONS = (
    "D3_AE",
    "D3_AE+IAG",
    "D2_LM+IAG->D3_AE+IAG",
    "D3_AE+TDI",
    "D2_LM+TDI->D3_AE+TDI",
)
FIGURE3F_CONDITION_LABELS = {
    "D3_AE": "D3 Control",
    "D3_AE+IAG": "D2-3 HIPPOa",
    "D2_LM+IAG->D3_AE+IAG": "D1-3 HIPPOa",
    "D3_AE+TDI": "D2-3 HIPPOi",
    "D2_LM+TDI->D3_AE+TDI": "D1-3 HIPPOi",
}
FIGURE3F_CONDITION_ORDER = tuple(
    FIGURE3F_CONDITION_LABELS[condition]
    for condition in FIGURE3F_D3_CONDITIONS
)
FIGURE3F_GENES = OrderedDict(
    {
        "Endothelial": ("CDH5", "NRP1", "EFNB2"),
        "Mesenchymal": ("KRT8", "TAGLN", "ANKRD1"),
    }
)
FIGURE3F_PROGRAMS = ("P2", "P34")


def _normalize_token(value: Any) -> str | None:
    if pd.isna(value):
        return None
    text = re.sub(r"[^A-Z0-9]+", "", str(value).strip().upper())
    return text or None


def _parse_program(value: Any) -> int | None:
    match = re.search(r"(\d+)", str(value))
    return int(match.group(1)) if match else None


def _split_features(value: Any) -> list[str]:
    if pd.isna(value):
        return []
    output: list[str] = []
    seen: set[str] = set()
    for token in re.split(r"[;,/]", str(value)):
        token = re.sub(r"\([^)]*\)", "", token).strip().replace("?", "")
        token = token.split()[0] if token else ""
        key = _normalize_token(token)
        if key and key not in seen:
            output.append(token)
            seen.add(key)
    return output


def plot_figure3b_feature_umaps(
    adata,
    usage: pd.DataFrame,
    *,
    point_size: float = 0.42,
    background_size: float = 0.18,
    vmax_quantile: float = 0.99,
):
    """Render P2/P34 and CDH5/ACTC1 with the recovered Figure 2e plotnine style."""
    if "X_umap" not in adata.obsm:
        raise KeyError("AnnData lacks obsm['X_umap']")
    coords = np.asarray(adata.obsm["X_umap"])
    feature_order = [label for _, _, label in FIGURE3B_FEATURES]
    background_frames: list[pd.DataFrame] = []
    foreground_frames: list[pd.DataFrame] = []
    source: dict[str, pd.DataFrame] = {}

    for kind, feature, label in FIGURE3B_FEATURES:
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


def compute_figure3b_carrier_seeded(
    adata,
    usage: pd.DataFrame,
    guide_metadata_path: str | Path,
    *,
    targets: Sequence[str] = FIGURE3_HIPPO_TARGETS,
    k: int = 300,
    pseudocount: float = 0.5,
):
    """Reuse the recovered same-day carrier-seeded analysis for Hippo targets."""
    return compute_figure2f_carrier_seeded(
        adata,
        usage,
        guide_metadata_path,
        targets=targets,
        k=k,
        pseudocount=pseudocount,
    )


def plot_figure3b_carrier_seeded(
    scores: pd.DataFrame,
    background: pd.DataFrame,
    *,
    targets: Sequence[str] = FIGURE3_HIPPO_TARGETS,
):
    """Reuse the recovered plotnine carrier-neighborhood renderer."""
    return plot_figure2f(
        scores,
        background,
        targets=targets,
    )


def compute_figure3c_guide_ec_120_controls(
    adata,
    guide_metadata_path: str | Path,
    *,
    targets: Sequence[str] = FIGURE3_HIPPO_TARGETS,
):
    """Compute guide-level D3 EC percentages for Figure 3c.

    All 120 negative controls (60 non-targeting plus 60 safe-targeting) are
    included. For each individual control guide, cells carrying any guide
    targeting NF2, LATS2, AMOTL2, or TAOK1 are excluded. Targeting-guide
    percentages use all D3 carriers of the corresponding guide. Welch tests
    compare each target's guide values with the common control-guide values;
    BH correction is across the four displayed genes.
    """
    return compute_figure2d_guide_ec_120_controls(
        adata,
        guide_metadata_path,
        targets=targets,
        control_exclusion_targets=targets,
    )


def plot_figure3c_120_controls(
    guide_table: pd.DataFrame,
    stats: pd.DataFrame,
    *,
    targets: Sequence[str] = FIGURE3_HIPPO_TARGETS,
):
    """Render Figure 3c with plotnine and the manuscript Hippo-pathway color."""
    targets = tuple(targets)
    order = ["Negative controls", *targets]
    positions = {target: index for index, target in enumerate(order)}
    plot_data = (
        guide_table.loc[guide_table["target"].isin(order)]
        .dropna(subset=["pct_EC_D3"])
        .copy()
    )
    plot_data["x_center"] = plot_data["target"].map(positions).astype(float)

    jitter_parts: list[pd.DataFrame] = []
    for target in order:
        frame = plot_data.loc[plot_data["target"].eq(target)].copy()
        seed = int.from_bytes(
            hashlib.sha256(f"Figure3c::{target}".encode()).digest()[:4],
            "little",
        )
        rng = np.random.default_rng(seed)
        spread = 0.10 if target == "Negative controls" else 0.06
        frame["x"] = frame["x_center"] + np.clip(
            rng.normal(0, spread, len(frame)), -0.18, 0.18
        )
        jitter_parts.append(frame)
    plot_data = pd.concat(jitter_parts, ignore_index=True)

    stats_index = stats.set_index("gene")
    star_rows = [
        {
            "x": positions[target],
            "y": 96.0,
            "label": str(stats_index.loc[target, "significance"]),
        }
        for target in targets
        if target in stats_index.index
        and str(stats_index.loc[target, "significance"])
    ]
    star_df = pd.DataFrame(star_rows, columns=["x", "y", "label"])
    control_mean = float(
        plot_data.loc[
            plot_data["target"].eq("Negative controls"), "pct_EC_D3"
        ].mean()
    )
    fill_values = {
        "Negative controls": "white",
        **{target: "#D33682" for target in targets},
    }

    plot = (
        ggplot()
        + geom_hline(
            yintercept=control_mean,
            color="#999999",
            linetype="dashed",
            size=0.8,
        )
        + geom_boxplot(
            plot_data,
            aes(
                x="x_center",
                y="pct_EC_D3",
                group="target",
                fill="target",
            ),
            inherit_aes=False,
            width=0.48,
            outlier_alpha=0,
            color="black",
            size=0.65,
        )
        + geom_point(
            plot_data,
            aes(x="x", y="pct_EC_D3"),
            inherit_aes=False,
            color="#555555",
            size=1.6,
            alpha=0.58,
            stroke=0,
        )
        + geom_text(
            star_df,
            aes(x="x", y="y", label="label"),
            inherit_aes=False,
            family="Liberation Sans",
            size=10,
            fontweight="bold",
        )
        + scale_fill_manual(values=fill_values, guide=None)
        + scale_x_continuous(
            breaks=list(positions.values()),
            labels=[
                f"{target}\nn={int((plot_data['target'] == target).sum())}"
                for target in order
            ],
            limits=(-0.65, len(order) - 0.35),
        )
        + scale_y_continuous(
            breaks=np.arange(0, 101, 20),
            limits=(0, 105),
        )
        + labs(x="gRNAs", y="Day 3 EC cells (%)")
        + theme_classic(base_family="Liberation Sans", base_size=10)
        + theme(
            figure_size=(4.8, 4.2),
            axis_text_x=element_text(rotation=55, ha="right", size=8),
            axis_text_y=element_text(size=8),
            axis_title=element_text(size=10),
            legend_position="none",
            panel_grid=element_blank(),
            plot_background=element_rect(fill="white", color="white"),
            panel_background=element_rect(fill="white", color="white"),
        )
        + annotate(
            "text",
            x=-0.55,
            y=102.5,
            label="c",
            family="Liberation Sans",
            fontweight="bold",
            size=15,
        )
    )
    jitter = plot_data[
        ["target", "guide_id", "control_type", "pct_EC_D3", "x"]
    ].copy()
    return plot, jitter


def _read_spectra(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    attempts = (
        {"sep": "\t"},
        {"sep": r"\s+", "engine": "python"},
        {"sep": ","},
    )
    table = None
    errors: list[str] = []
    for kwargs in attempts:
        try:
            candidate = pd.read_csv(path, index_col=0, **kwargs)
            candidate = candidate.loc[
                :, ~candidate.columns.astype(str).str.startswith("Unnamed")
            ]
            if candidate.shape[0] >= 2 and candidate.shape[1] >= 2:
                table = candidate
                break
        except Exception as exc:
            errors.append(str(exc))
    if table is None:
        raise ValueError(f"Could not read cNMF spectra {path}: {errors}")
    if table.shape[0] > 200 and table.shape[1] <= 200:
        table = table.T
    table.columns = table.columns.astype(str).str.strip()
    if table.columns.duplicated().any():
        table = table.T.groupby(level=0).mean().T
    return table.apply(pd.to_numeric, errors="coerce")


def _program_row_lookup(spectra: pd.DataFrame) -> dict[int, Any]:
    output: dict[int, Any] = {}
    for row in spectra.index:
        program = _parse_program(row)
        if program is not None and program not in output:
            output[program] = row
    return output


def _choose_spectra_offset(spectra: pd.DataFrame) -> int:
    row_lookup = _program_row_lookup(spectra)
    gene_lookup = {_normalize_token(gene): gene for gene in spectra.columns}
    scores: dict[int, float] = {}
    for offset in (0, 1, -1):
        values: list[float] = []
        for program, genes in FIGURE3D_MARKERS.items():
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


def _load_regulator_tables(
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
    regulator_column = "Top 5 specific regulators (FDR<0.1) D3"
    if regulator_column not in summary.columns:
        raise KeyError(
            f"Workbook Summary is missing {regulator_column!r}; "
            f"available columns: {summary.columns.tolist()}"
        )

    selected_by_program: OrderedDict[int, list[str]] = OrderedDict()
    selection_rows: list[dict[str, Any]] = []
    for program in FIGURE3D_PROGRAMS:
        matched = summary.loc[summary["program_num"].eq(program)]
        if matched.empty:
            raise KeyError(f"Summary has no row for P{program}.")
        regulators = _split_features(
            matched.iloc[0][regulator_column]
        )[:3]
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
    for program in FIGURE3D_PROGRAMS:
        for gene in selected_by_program[program]:
            key = _normalize_token(gene)
            if key and key not in seen:
                regulator_union.append(gene)
                seen.add(key)

    pert = pd.read_excel(
        workbook_path, sheet_name="Perturbation Association 1"
    )
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
    d3 = pert.loc[
        pert["day_clean"].eq("D3")
        & pert["program_num"].isin(FIGURE3D_PROGRAMS)
    ].copy()

    effects = pd.DataFrame(
        index=FIGURE3D_PROGRAMS, columns=regulator_union, dtype=float
    )
    fdrs = effects.copy()
    stars = pd.DataFrame(
        "", index=FIGURE3D_PROGRAMS, columns=regulator_union, dtype=object
    )
    for program in FIGURE3D_PROGRAMS:
        for gene in regulator_union:
            subset = d3.loc[
                d3["program_num"].eq(program)
                & d3["gene_key"].eq(_normalize_token(gene))
            ].sort_values("fdr")
            if subset.empty:
                continue
            effects.loc[program, gene] = subset.iloc[0]["effect"]
            fdrs.loc[program, gene] = subset.iloc[0]["fdr"]
            stars.loc[program, gene] = significance_stars(
                subset.iloc[0]["fdr"]
            )
    return (
        effects,
        fdrs,
        stars,
        regulator_union,
        pd.DataFrame(selection_rows),
    )


def _load_marker_table(
    spectra_path: str | Path,
) -> tuple[pd.DataFrame, int, pd.DataFrame]:
    spectra = _read_spectra(spectra_path)
    offset = _choose_spectra_offset(spectra)
    row_lookup = _program_row_lookup(spectra)
    gene_lookup = {_normalize_token(gene): gene for gene in spectra.columns}
    marker_union: list[str] = []
    for program in FIGURE3D_PROGRAMS:
        for gene in FIGURE3D_MARKERS[program]:
            if gene not in marker_union:
                marker_union.append(gene)
    matrix = pd.DataFrame(
        index=FIGURE3D_PROGRAMS, columns=marker_union, dtype=float
    )
    source_rows: list[dict[str, Any]] = []
    for program in FIGURE3D_PROGRAMS:
        raw_program = program + offset
        if raw_program not in row_lookup:
            raise KeyError(
                f"No consensus row found for P{program} with offset {offset}."
            )
        row_label = row_lookup[raw_program]
        for rank, gene in enumerate(FIGURE3D_MARKERS[program], start=1):
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


def _infer_column(
    table: pd.DataFrame,
    candidates: Sequence[str],
    *,
    required: bool = True,
) -> str | None:
    exact = {str(column).strip().lower(): str(column) for column in table.columns}
    for candidate in candidates:
        if candidate.lower() in exact:
            return exact[candidate.lower()]
    for column in table.columns:
        if any(candidate.lower() in str(column).lower() for candidate in candidates):
            return str(column)
    if required:
        raise KeyError(
            f"Could not infer any of {list(candidates)} from "
            f"{table.columns.tolist()}"
        )
    return None


def _load_motif_tables(
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
        data["neg_log10_q"] = pd.to_numeric(
            table[neglog_col], errors="coerce"
        )
    else:
        data["neg_log10_q"] = -np.log10(
            data["q_value"].clip(lower=np.nextafter(0, 1))
        )
    data = data.loc[
        data["day"].eq("D3")
        & data["program"].isin(FIGURE3D_PROGRAMS)
        & data["q_value"].notna()
        & data["motif_key"].notna()
    ].copy()

    selected_by_program: OrderedDict[int, list[str]] = OrderedDict()
    selection_rows: list[dict[str, Any]] = []
    for program in FIGURE3D_PROGRAMS:
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
                f"Only {len(motifs)} motif rows were available for D3 P{program}; "
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
    for program in FIGURE3D_PROGRAMS:
        for motif in selected_by_program[program]:
            key = _normalize_token(motif)
            if key and key not in seen:
                motif_union.append(motif)
                seen.add(key)

    values = pd.DataFrame(
        index=FIGURE3D_PROGRAMS, columns=motif_union, dtype=float
    )
    qvalues = values.copy()
    stars = pd.DataFrame(
        "", index=FIGURE3D_PROGRAMS, columns=motif_union, dtype=object
    )
    for program in FIGURE3D_PROGRAMS:
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
    return (
        values,
        qvalues,
        stars,
        motif_union,
        pd.DataFrame(selection_rows),
    )


def _cluster_map_lookup(path: str | Path) -> dict[str, int]:
    table = pd.read_csv(path, sep="\t")
    id_candidates: list[str] = []
    for column in table.columns:
        numeric_fraction = pd.to_numeric(table[column], errors="coerce").notna().mean()
        if numeric_fraction > 0.80 and re.search(
            r"(cluster|pattern|cid|id)", str(column), flags=re.I
        ):
            id_candidates.append(str(column))
    if not id_candidates:
        for column in table.columns:
            numeric_fraction = pd.to_numeric(
                table[column], errors="coerce"
            ).notna().mean()
            if numeric_fraction > 0.80:
                id_candidates.append(str(column))
    if not id_candidates:
        raise ValueError("Could not infer the motif-cluster identifier column.")
    id_col = id_candidates[0]
    name_cols = [
        str(column)
        for column in table.columns
        if str(column) != id_col
        and table[column].dtype == object
        and not re.search(r"(path|file|url|image)", str(column), flags=re.I)
    ]
    lookup: dict[str, int] = {}
    for _, row in table.iterrows():
        cluster = pd.to_numeric(row[id_col], errors="coerce")
        if pd.isna(cluster):
            continue
        aliases: set[str] = set()
        for column in name_cols:
            value = row[column]
            if pd.isna(value):
                continue
            aliases.add(str(value))
            aliases.update(
                token.strip()
                for token in re.split(r"[;,/|]", str(value))
                if token.strip()
            )
        for alias in aliases:
            key = _normalize_token(alias)
            if key and key not in lookup:
                lookup[key] = int(cluster)
    return lookup


def _motif_variants(name: str) -> list[str]:
    base = str(name).strip()
    variants = {base, re.sub(r"_[0-9]+_$", "", base)}
    for value in list(variants):
        variants.update(
            {
                value.replace("_", ""),
                value.replace("-", ""),
                value.replace("/", ""),
            }
        )
    variants.update(
        token.strip()
        for token in re.split(r"[;/|]", base)
        if token.strip()
    )
    return [value for value in variants if value]


def _best_cluster(name: str, lookup: dict[str, int]) -> int | None:
    for variant in _motif_variants(name):
        key = _normalize_token(variant)
        if key in lookup:
            return lookup[key]
    motif_keys = [
        key
        for key in (_normalize_token(value) for value in _motif_variants(name))
        if key
    ]
    for key in motif_keys:
        matches = [
            (len(alias), cluster)
            for alias, cluster in lookup.items()
            if key.startswith(alias) or alias.startswith(key)
        ]
        if matches:
            return sorted(matches, reverse=True)[0][1]
    return None


def _load_logo_objects(
    motif_names: Sequence[str],
    cluster_map_path: str | Path,
    modisco_h5_path: str | Path,
) -> tuple[dict[str, dict[str, Any]], pd.DataFrame]:
    import h5py

    lookup = _cluster_map_lookup(cluster_map_path)
    cwms: dict[int, np.ndarray] = {}
    with h5py.File(modisco_h5_path, "r") as handle:
        if "pos_patterns" not in handle:
            raise KeyError("MoDISco H5 lacks /pos_patterns.")
        for key in handle["pos_patterns"].keys():
            try:
                cluster = int(key)
            except ValueError:
                parsed = _parse_program(key)
                if parsed is None:
                    continue
                cluster = parsed
            group = handle["pos_patterns"][key]
            if "contrib_scores" in group:
                cwms[cluster] = group["contrib_scores"][:]

    objects: dict[str, dict[str, Any]] = {}
    rows: list[dict[str, Any]] = []
    for motif in motif_names:
        cluster = _best_cluster(motif, lookup)
        status = "text_fallback_no_cluster"
        object_data: dict[str, Any] = {
            "type": "text",
            "label": motif,
            "cluster_id": cluster,
        }
        if cluster is not None and cluster in cwms:
            array = np.asarray(cwms[cluster])
            if array.ndim == 2 and 4 in array.shape:
                matrix = array if array.shape[1] == 4 else array.T
                object_data = {
                    "type": "cwm",
                    "label": motif,
                    "cluster_id": cluster,
                    "cwm_df": pd.DataFrame(matrix, columns=list("ACGT")),
                }
                status = "vector_cwm"
        objects[motif] = object_data
        rows.append(
            {
                "motif": motif,
                "matched_cluster_id": cluster,
                "status": status,
            }
        )
    return objects, pd.DataFrame(rows)


def _text_color(value: float, cmap, norm) -> str:
    if pd.isna(value):
        return "black"
    rgba = cmap(norm(value))
    luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
    return "white" if luminance < 0.48 else "black"


def _draw_rect_heatmap(
    axis,
    values: pd.DataFrame,
    *,
    cmap,
    norm,
    xlabels: Sequence[str],
    ylabels: Sequence[str] | None = None,
    stars: pd.DataFrame | None = None,
    show_ylabels: bool = True,
    italic_x: bool = False,
):
    n_rows, n_cols = values.shape
    array = values.to_numpy(float)
    for row in range(n_rows):
        for column in range(n_cols):
            value = array[row, column]
            color = (0.95, 0.95, 0.95, 1.0) if pd.isna(value) else cmap(norm(value))
            axis.add_patch(
                Rectangle(
                    (column, row),
                    1,
                    1,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=1.0,
                )
            )
            if stars is not None:
                star = stars.iloc[row, column]
                if isinstance(star, str) and star:
                    axis.text(
                        column + 0.5,
                        row + 0.5,
                        star,
                        ha="center",
                        va="center",
                        fontsize=7,
                        fontweight="bold",
                        color=_text_color(value, cmap, norm),
                    )
    axis.set_xlim(0, n_cols)
    axis.set_ylim(n_rows, 0)
    axis.set_aspect("equal")
    axis.set_xticks(np.arange(n_cols) + 0.5)
    axis.set_xticklabels(
        xlabels,
        rotation=58,
        ha="right",
        rotation_mode="anchor",
        fontsize=7.5,
        fontstyle="italic" if italic_x else "normal",
    )
    axis.set_yticks(np.arange(n_rows) + 0.5)
    axis.set_yticklabels(
        ylabels if show_ylabels and ylabels is not None else [""] * n_rows,
        fontsize=8,
    )
    axis.tick_params(length=0)
    for spine in axis.spines.values():
        spine.set_visible(False)


def _add_horizontal_colorbar(
    figure,
    *,
    left: float,
    bottom: float,
    width: float,
    height: float,
    cmap,
    norm,
    label: str,
):
    axis = figure.add_axes([left, bottom, width, height])
    scalar = ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    colorbar = figure.colorbar(scalar, cax=axis, orientation="horizontal")
    colorbar.outline.set_linewidth(0.7)
    colorbar.ax.tick_params(labelsize=7, length=2)
    colorbar.set_label(label, fontsize=8, labelpad=2)


def _draw_rotated_logo(
    figure,
    source_axis,
    motif: str,
    logo_objects: dict[str, dict[str, Any]],
    *,
    center: tuple[float, float],
    target_width: float,
    target_height: float,
):
    object_data = logo_objects.get(motif, {"type": "text", "label": motif})
    if object_data.get("type") != "cwm":
        source_axis.remove()
        figure.text(
            center[0],
            center[1],
            str(object_data.get("label", motif)),
            ha="center",
            va="center",
            rotation=90,
            fontsize=6.2,
            fontweight="bold",
        )
        return
    try:
        import logomaker

        logomaker.Logo(
            object_data["cwm_df"],
            ax=source_axis,
            color_scheme="classic",
            fade_below=0.0,
            shade_below=0.0,
            font_name="DejaVu Sans",
        )
        source_axis.set_axis_off()
        source_axis.patch.set_alpha(0.0)
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        patches = [patch for patch in source_axis.patches if patch.get_visible()]
        boxes = []
        for patch in patches:
            box = patch.get_window_extent(renderer=renderer)
            if (
                np.isfinite([box.x0, box.x1, box.y0, box.y1]).all()
                and box.width > 0
                and box.height > 0
            ):
                boxes.append(box)
        if not boxes:
            raise RuntimeError("No visible logo glyph bounds.")
        content = Bbox.union(boxes)
        source_center = (
            (content.x0 + content.x1) / 2,
            (content.y0 + content.y1) / 2,
        )
        target_center = figure.transFigure.transform(center)
        target_width_px = target_width * figure.bbox.width
        target_height_px = target_height * figure.bbox.height
        scale = min(
            target_width_px / max(content.height, 1e-12),
            target_height_px / max(content.width, 1e-12),
        ) * 0.96
        transform = (
            Affine2D()
            .translate(-source_center[0], -source_center[1])
            .rotate_deg(90)
            .scale(scale)
            .translate(target_center[0], target_center[1])
        )
        for patch in patches:
            patch.set_transform(patch.get_transform() + transform)
            patch.set_clip_on(False)
            patch.set_clip_path(None)
            patch.set_zorder(7)
        source_axis.set_axis_off()
        source_axis.patch.set_visible(False)
    except Exception:
        for patch in list(source_axis.patches):
            patch.remove()
        source_axis.remove()
        figure.text(
            center[0],
            center[1],
            motif,
            ha="center",
            va="center",
            rotation=90,
            fontsize=6.2,
            fontweight="bold",
        )


def _add_logo_strip(
    figure,
    motif_axis,
    motif_names: Sequence[str],
    logo_objects: dict[str, dict[str, Any]],
):
    if not motif_names:
        return
    position = motif_axis.get_position()
    column_width = position.width / len(motif_names)
    logo_height = 0.18
    logo_bottom = position.y1 + 0.025
    for index, motif in enumerate(motif_names):
        center_x = position.x0 + (index + 0.5) * column_width
        center_y = logo_bottom + logo_height / 2
        source_axis = figure.add_axes(
            [
                center_x - logo_height / 2,
                center_y - column_width * 0.42,
                logo_height,
                column_width * 0.84,
            ]
        )
        _draw_rotated_logo(
            figure,
            source_axis,
            motif,
            logo_objects,
            center=(center_x, center_y),
            target_width=column_width * 0.84,
            target_height=logo_height * 0.92,
        )


def make_figure3d_annotation(
    *,
    workbook_path: str | Path,
    spectra_path: str | Path,
    motif_path: str | Path,
    cluster_map_path: str | Path,
    modisco_h5_path: str | Path,
    top_motifs_per_program: int = 5,
):
    """Reproduce the focused D3 P2/P34 annotation workflow."""
    effects, effect_fdrs, effect_stars, regulators, regulator_selection = (
        _load_regulator_tables(workbook_path)
    )
    markers, spectra_offset, marker_selection = _load_marker_table(spectra_path)
    (
        motif_values,
        motif_qvalues,
        motif_stars,
        motifs,
        motif_selection,
    ) = _load_motif_tables(motif_path, top_n=top_motifs_per_program)
    logo_objects, logo_map = _load_logo_objects(
        motifs, cluster_map_path, modisco_h5_path
    )

    regulator_values = np.abs(effects.to_numpy(float))
    regulator_values = regulator_values[np.isfinite(regulator_values)]
    regulator_vmax = (
        float(np.max(regulator_values)) if regulator_values.size else 1.0
    )
    marker_values = markers.to_numpy(float)
    marker_values = marker_values[np.isfinite(marker_values)]
    marker_vmax = float(np.max(marker_values)) if marker_values.size else 1.0
    motif_finite = motif_values.to_numpy(float)
    motif_finite = motif_finite[np.isfinite(motif_finite)]
    motif_raw_vmax = float(np.max(motif_finite)) if motif_finite.size else 1.0
    motif_vmax = (
        float(np.percentile(motif_finite, 95))
        if motif_finite.size
        else 1.0
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
    labels = [FIGURE3D_LABELS[program] for program in FIGURE3D_PROGRAMS]
    _draw_rect_heatmap(
        axes[0],
        effects.loc[list(FIGURE3D_PROGRAMS), regulators],
        cmap=regulator_cmap,
        norm=regulator_norm,
        xlabels=regulators,
        ylabels=labels,
        stars=effect_stars.loc[list(FIGURE3D_PROGRAMS), regulators],
        show_ylabels=True,
    )
    _draw_rect_heatmap(
        axes[1],
        markers.loc[list(FIGURE3D_PROGRAMS)],
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
        motif_values.loc[list(FIGURE3D_PROGRAMS), motifs],
        cmap=motif_cmap,
        norm=motif_norm,
        xlabels=motifs,
        ylabels=labels,
        stars=motif_stars.loc[list(FIGURE3D_PROGRAMS), motifs],
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
        "d",
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
        "programs": list(FIGURE3D_PROGRAMS),
        "top_regulators_per_program": 3,
        "top_motifs_per_program": top_motifs_per_program,
        "motif_selection_requires_significance": False,
        "consensus_program_offset": spectra_offset,
        "motif_raw_max": motif_raw_vmax,
        "motif_display_vmax_p95": motif_vmax,
        "motif_value": "-log10(BH q-value)",
    }
    return figure, source_tables, metadata



def _hippo_aligned_metadata(query, projected) -> pd.DataFrame:
    """Validate the canonical expression/projection pair and align by cell ID.

    The rebuilt expression object deliberately has only the broad producer QC.
    Reconstruct the later condition-specific flag from the recovered projection
    definition, checking any stored annotations instead of silently skipping QC.
    No gene-expression matrices are copied while checking the metadata.
    """
    from .hippo_projection import add_projection_metadata

    query_ids = pd.Index(query.obs_names.astype(str))
    projected_ids = pd.Index(projected.obs_names.astype(str))
    if not query_ids.is_unique or not projected_ids.is_unique:
        raise ValueError("Hippo expression and projection cell IDs must be unique.")
    missing_projection = query_ids.difference(projected_ids)
    missing_query = projected_ids.difference(query_ids)
    if len(missing_projection) or len(missing_query):
        raise ValueError(
            "The canonical Hippo query and projection must cover the same cells; "
            f"{len(missing_projection)} query cells lack projection and "
            f"{len(missing_query)} projected cells lack expression. "
            "Run notebooks/preprocessing/Rebuild_Hippo_Timecourse.ipynb."
        )

    annotated = []
    for label, obj in (("query", query), ("projection", projected)):
        required = ("sample_name", "sublibID", "pct_counts_mt", "pct_counts_ribo")
        missing = [column for column in required if column not in obj.obs]
        if missing:
            raise KeyError(f"Hippo {label} is missing metadata: {missing}")
        obs = add_projection_metadata(obj.obs)
        for column in ("condition", "replicate", "treatment", "qc_pass_condition_specific"):
            if column in obj.obs:
                stored = obj.obs[column]
                rebuilt = obs[column]
                equal = (stored == rebuilt) | (stored.isna() & rebuilt.isna())
                if not bool(equal.fillna(False).all()):
                    raise ValueError(
                        f"Hippo {label} {column!r} disagrees with the recovered "
                        "sample mapping or condition-specific QC definition."
                    )
        if obs["sublibID"].isna().any():
            raise ValueError(f"Hippo {label} has missing sublibrary identifiers.")
        annotated.append(obs.loc[projected_ids].copy())

    query_obs, projected_obs = annotated
    compare = (
        "sample_name", "biological_sample", "sublibID", "plate", "well",
        "condition", "replicate", "treatment", "qc_pass_condition_specific",
        "pct_counts_mt", "pct_counts_ribo", "n_genes_by_counts", "total_counts",
    )
    for column in compare:
        if column not in query_obs or column not in projected_obs:
            continue
        left, right = query_obs[column], projected_obs[column]
        # Categorical storage may differ even though its values are identical.
        if isinstance(left.dtype, pd.CategoricalDtype):
            left = left.astype(object)
        if isinstance(right.dtype, pd.CategoricalDtype):
            right = right.astype(object)
        equal = (left == right) | (left.isna() & right.isna())
        if not bool(equal.fillna(False).all()):
            raise ValueError(f"Hippo cell metadata disagree for {column!r}.")
    return projected_obs


def _hippo_panel_mask(obs: pd.DataFrame, conditions: Sequence[str]) -> np.ndarray:
    """Apply the original panel filters and reject empty requested conditions."""
    keep = (
        obs["condition"].astype(str).isin(conditions).to_numpy()
        & obs["qc_pass_condition_specific"].to_numpy(dtype=bool)
        & ~obs["sublibID"].astype(str).eq("s8").to_numpy()
    )
    retained = set(obs.loc[keep, "condition"].astype(str))
    missing = [condition for condition in conditions if condition not in retained]
    if missing:
        raise ValueError(
            "Hippo conditions have no cells after condition-specific QC and "
            f"excluding sublibrary s8: {missing}"
        )
    return keep


def _hippo_program_values(projected, programs: Sequence[str]) -> np.ndarray:
    """Resolve selected CLR columns using explicit zero-based program labels."""
    from scipy import sparse

    key = "X_cnmf_clr_top_specific_200"
    if key not in projected.obsm:
        raise KeyError(f"Reprojected AnnData lacks obsm[{key!r}].")
    values = projected.obsm[key]
    expected_names = pd.Index([f"P{i}" for i in range(50)])
    if values.shape != (projected.n_obs, 50):
        raise ValueError(f"Hippo {key} must be cells by 50 programs; got {values.shape}.")
    if isinstance(values, pd.DataFrame):
        columns = pd.Index(values.columns.astype(str))
        if not columns.is_unique or set(columns) != set(expected_names):
            raise ValueError("Labelled Hippo CLR columns must be P0 through P49.")
        matrix = values.to_numpy(dtype=float)
        positions = [int(columns.get_loc(program)) for program in programs]
    else:
        names = pd.Index(projected.var_names.astype(str))
        if not names.equals(expected_names):
            raise ValueError("Hippo CLR array requires var_names ordered P0 through P49.")
        if "program" not in projected.var or "program_id_zero_based" not in projected.var:
            raise ValueError(
                "Hippo CLR array lacks explicit program/zero-based column metadata. "
                "Rebuild the canonical projected object."
            )
        labels = projected.var["program"].astype(str).to_numpy()
        ids = pd.to_numeric(projected.var["program_id_zero_based"], errors="raise").to_numpy()
        if not np.array_equal(labels, expected_names.to_numpy()) or not np.array_equal(ids, np.arange(50)):
            raise ValueError("Hippo program-column metadata disagree with P0 through P49.")
        matrix = values.toarray() if sparse.issparse(values) else np.asarray(values, dtype=float)
        positions = [int(names.get_loc(program)) for program in programs]
    if not np.isfinite(matrix).all():
        raise ValueError("Hippo top-specific-200 CLR matrix contains non-finite values.")
    return matrix[:, positions]


def prepare_figure3f_dotplots(
    query_h5ad_path: str | Path,
    reprojected_h5ad_path: str | Path,
):
    """Prepare Figure 3f dotplots using the attached Hippo workflow."""
    import anndata as ad
    import scanpy as sc
    from scipy import sparse

    query = sc.read_h5ad(query_h5ad_path)
    projected = sc.read_h5ad(reprojected_h5ad_path)
    aligned_obs = _hippo_aligned_metadata(query, projected)
    n_common_cells = len(aligned_obs)
    keep = _hippo_panel_mask(aligned_obs, FIGURE3F_D3_CONDITIONS)
    projected = projected[keep].copy()
    query = query[projected.obs_names].copy()
    obs = aligned_obs.loc[projected.obs_names].copy()
    obs["condition_name"] = pd.Categorical(
        obs["condition"].astype(str).map(FIGURE3F_CONDITION_LABELS),
        categories=list(FIGURE3F_CONDITION_ORDER),
        ordered=True,
    )
    projected.obs["condition_name"] = obs["condition_name"]
    query.obs["condition_name"] = obs["condition_name"]

    if "gene" in query.var.columns:
        query.var_names = query.var["gene"].astype(str)
        query.var_names_make_unique()
    all_genes = [
        gene for group in FIGURE3F_GENES.values() for gene in group
    ]
    missing_genes = [gene for gene in all_genes if gene not in query.var_names]
    if missing_genes:
        raise KeyError(f"Query AnnData is missing marker genes: {missing_genes}")
    gene_adata = query.copy()
    sc.pp.normalize_total(gene_adata)
    sc.pp.log1p(gene_adata)
    gene_adata = gene_adata[:, all_genes].copy()

    clr_key = "X_cnmf_clr_top_specific_200"
    clr = _hippo_program_values(projected, FIGURE3F_PROGRAMS)
    program_adata = ad.AnnData(
        X=clr.copy(),
        obs=obs.copy(),
        var=pd.DataFrame(index=pd.Index(FIGURE3F_PROGRAMS, name="program")),
    )
    program_adata.layers["clr_top_specific_200"] = clr.copy()

    gene_dotplot = sc.pl.dotplot(
        gene_adata,
        var_names={key: list(value) for key, value in FIGURE3F_GENES.items()},
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
        var_names=list(FIGURE3F_PROGRAMS),
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
    for condition in FIGURE3F_CONDITION_ORDER:
        mask = (
            gene_adata.obs["condition_name"].astype(str).to_numpy()
            == condition
        )
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
    for condition in FIGURE3F_CONDITION_ORDER:
        mask = (
            program_adata.obs["condition_name"].astype(str).to_numpy()
            == condition
        )
        for index, program in enumerate(FIGURE3F_PROGRAMS):
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
        "n_common_cells": n_common_cells,
        "n_D3_QC_cells": int(keep.sum()),
        "conditions": list(FIGURE3F_CONDITION_ORDER),
        "program_projection": clr_key,
        "program_projection_gene_set": "top_specific_200",
        "gene_standard_scale": "var",
        "excluded_sublibrary": "s8 if present",
    }
    return (
        gene_dotplot,
        program_dotplot,
        pd.DataFrame(gene_rows),
        pd.DataFrame(program_rows),
        metadata,
    )

