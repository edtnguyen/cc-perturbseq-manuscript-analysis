from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path
import gzip
import re
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from .program_io import clean_day, infer_column, read_table
from .screen_panels import PKC_GROUPS


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def _gtf_attributes(text: str) -> dict[str, str]:
    return {
        key: value
        for key, value in re.findall(r'(\S+)\s+"([^"]+)"', text)
    }


def _union_length(intervals: Sequence[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    ordered = sorted(intervals)
    total = 0
    start, end = ordered[0]
    for next_start, next_end in ordered[1:]:
        if next_start <= end + 1:
            end = max(end, next_end)
        else:
            total += end - start + 1
            start, end = next_start, next_end
    return total + end - start + 1


def union_exon_lengths(
    gtf_path: str | Path,
    *,
    wanted_symbols: Sequence[str] | None = None,
    wanted_ids: Sequence[str] | None = None,
) -> tuple[dict[str, int], dict[str, int]]:
    path = Path(gtf_path)
    symbols = set(map(str, wanted_symbols or []))
    identifiers = {str(value).split(".")[0] for value in (wanted_ids or [])}
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    aliases: dict[str, tuple[str, str]] = {}
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attrs = _gtf_attributes(fields[8])
            gene_id = attrs.get("gene_id", "").split(".")[0]
            gene_name = attrs.get("gene_name", "")
            if symbols or identifiers:
                if gene_name not in symbols and gene_id not in identifiers:
                    continue
            key = gene_id or gene_name
            intervals[key].append((int(fields[3]), int(fields[4])))
            aliases[key] = (gene_id, gene_name)
    by_id: dict[str, int] = {}
    by_symbol: dict[str, int] = {}
    for key, values in intervals.items():
        length = _union_length(values)
        gene_id, gene_name = aliases[key]
        if gene_id:
            by_id[gene_id] = length
        if gene_name:
            by_symbol[gene_name] = length
    return by_id, by_symbol


def calculate_pkc_pseudobulk_tpm(
    adata,
    gtf_path: str | Path,
    *,
    day_col: str = "day",
    count_layer: str = "raw_counts",
    pkc_groups: Mapping[str, Sequence[str]] = PKC_GROUPS,
) -> pd.DataFrame:
    genes = [gene for group in pkc_groups.values() for gene in group]
    symbol_column = next(
        (column for column in ("gene_symbol", "gene_name", "symbol") if column in adata.var.columns),
        None,
    )
    symbols = pd.Index(
        adata.var[symbol_column].astype(str) if symbol_column else adata.var_names.astype(str)
    )
    gene_id_column = next(
        (column for column in ("ensembl_id", "gene_id") if column in adata.var.columns),
        None,
    )
    ids = pd.Index(
        adata.var[gene_id_column].astype(str).str.split(".").str[0]
        if gene_id_column
        else adata.var_names.astype(str).str.split(".").str[0]
    )
    by_id, by_symbol = union_exon_lengths(
        gtf_path,
        wanted_symbols=genes,
        wanted_ids=ids[symbols.isin(genes)],
    )
    matrix = adata.layers[count_layer] if count_layer in adata.layers else adata.X
    days = pd.Series(adata.obs[day_col], index=adata.obs_names).map(clean_day)
    rows: list[dict[str, Any]] = []
    for day in ("D0", "D1", "D2", "D3"):
        cell_mask = days.eq(day).to_numpy()
        if not cell_mask.any():
            continue
        for group, group_genes in pkc_groups.items():
            for gene in group_genes:
                positions = np.flatnonzero(symbols == gene)
                if not len(positions):
                    rows.append({"day": day, "group": group, "gene": gene, "counts": np.nan, "tpm": np.nan})
                    continue
                selected = matrix[cell_mask][:, positions]
                counts = float(selected.sum())
                length = by_symbol.get(gene)
                if length is None:
                    for position in positions:
                        length = by_id.get(ids[position])
                        if length:
                            break
                rows.append(
                    {
                        "day": day,
                        "group": group,
                        "gene": gene,
                        "counts": counts,
                        "gene_length_bp": length,
                    }
                )
    output = pd.DataFrame(rows)
    output["rpk"] = output.counts / (output.gene_length_bp / 1000.0)
    for day, index in output.groupby("day").groups.items():
        scale = output.loc[index, "rpk"].sum() / 1e6
        output.loc[index, "tpm"] = output.loc[index, "rpk"] / scale if scale > 0 else np.nan
    return output


def plot_pkc_tpm(table: pd.DataFrame):
    groups = list(PKC_GROUPS)
    fig, axes = plt.subplots(
        len(groups),
        1,
        figsize=(6.2, 5.0),
        sharex=True,
        constrained_layout=True,
    )
    for axis, group in zip(axes, groups):
        data = table[table.group == group]
        for gene, subset in data.groupby("gene", sort=False):
            subset = subset.set_index("day").reindex(["D0", "D1", "D2", "D3"]).reset_index()
            axis.plot(subset.day, subset.tpm, marker="o", markersize=3, linewidth=1, label=gene)
        axis.set_ylabel("Pseudobulk TPM", fontsize=7)
        axis.set_title(group, fontsize=8, loc="right")
        axis.tick_params(labelsize=6)
    axes[-1].set_xlabel("Differentiation day")
    axes[1].legend(frameon=False, fontsize=5.5, ncol=2, bbox_to_anchor=(1.01, 1), loc="upper left")
    return fig


def standardize_knockdown_results(path: str | Path) -> pd.DataFrame:
    table = read_table(path)
    target_col = infer_column(table, ["target", "gene", "target_gene", "response_id"])
    day_col = infer_column(table, ["day", "timepoint"])
    fold_col = infer_column(table, ["fold_change", "relative_expression", "fc"])
    se_col = infer_column(table, ["se_fold_change", "fold_change_se", "se"], required=False)
    q_col = infer_column(table, ["q_value", "qvalue", "FDR", "padj", "adjusted_p_value"], required=False)
    output = pd.DataFrame(
        {
            "target": table[target_col].astype(str),
            "day": table[day_col].map(clean_day),
            "fold_change": pd.to_numeric(table[fold_col], errors="coerce"),
            "se_fold_change": pd.to_numeric(table[se_col], errors="coerce") if se_col else np.nan,
            "q_value": pd.to_numeric(table[q_col], errors="coerce") if q_col else np.nan,
        }
    )
    output = output.dropna(subset=["target", "day", "fold_change"])
    z = 1.96
    valid = (output.fold_change > 0) & output.se_fold_change.notna()
    log2_fc = np.log2(output.loc[valid, "fold_change"])
    se_log2 = output.loc[valid, "se_fold_change"] / (output.loc[valid, "fold_change"] * np.log(2))
    output.loc[valid, "lower"] = 2 ** (log2_fc - z * se_log2)
    output.loc[valid, "upper"] = 2 ** (log2_fc + z * se_log2)
    return output


def plot_knockdown_across_days(
    table: pd.DataFrame,
    *,
    targets: Sequence[str] = ("NF2", "LATS2", "AMOTL2", "TAOK1"),
):
    data = table[table.target.isin(targets)].copy()
    day_position = {day: index for index, day in enumerate(("D0", "D1", "D2", "D3"))}
    fig, ax = plt.subplots(figsize=(7.0, 2.9))
    centers = np.arange(len(targets)) * 5
    rows: list[dict[str, Any]] = []
    for target_index, target in enumerate(targets):
        subset = data[data.target == target].copy()
        subset["x"] = centers[target_index] + subset.day.map(day_position)
        significant = subset.q_value <= 0.05
        ax.errorbar(
            subset.x,
            subset.fold_change,
            yerr=[
                subset.fold_change - subset.lower,
                subset.upper - subset.fold_change,
            ] if subset.lower.notna().any() else None,
            fmt="none",
            ecolor="0.45",
            elinewidth=0.9,
            capsize=2,
        )
        ax.scatter(
            subset.x,
            subset.fold_change,
            s=25,
            color=np.where(significant, "#D95F02", "0.55"),
            zorder=3,
        )
        for _, row in subset.iterrows():
            rows.append(row.to_dict())
        ax.text(centers[target_index] + 1.5, -0.08, target, ha="center", va="top", transform=ax.get_xaxis_transform(), fontsize=7, fontstyle="italic")
    tick_positions = [centers[i] + day_position[day] for i in range(len(targets)) for day in ("D0", "D1", "D2", "D3")]
    tick_labels = [day for _ in targets for day in ("D0", "D1", "D2", "D3")]
    ax.set_xticks(tick_positions, tick_labels, fontsize=6)
    ax.axhline(1.0, color="0.5", linestyle="--", linewidth=0.8)
    ax.set_ylabel("Target-gene expression\nrelative to control")
    ax.set_ylim(bottom=0)
    return fig, pd.DataFrame(rows)
