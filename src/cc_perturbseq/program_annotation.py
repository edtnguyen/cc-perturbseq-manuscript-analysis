from __future__ import annotations

from collections.abc import Mapping, Sequence

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

from .program_io import clean_day, significance_stars


def _top_genes(
    spectra: pd.DataFrame,
    programs: Sequence[int],
    n: int = 4,
) -> list[str]:
    genes: list[str] = []
    for program in programs:
        if program not in spectra.index:
            continue
        for gene in spectra.loc[program].nlargest(n).index:
            if gene not in genes:
                genes.append(str(gene))
    return genes


def _top_motifs(
    motifs: pd.DataFrame,
    day: str,
    programs: Sequence[int],
    n: int = 3,
) -> list[str]:
    data = motifs[
        (motifs.day == day) & motifs.program.isin(programs)
    ].copy()
    data["strength"] = -np.log10(
        data.q_value.clip(lower=np.nextafter(0, 1))
    )
    selected: list[str] = []
    for program in programs:
        for motif in data[data.program == program].nlargest(n, "strength").motif:
            if motif not in selected:
                selected.append(str(motif))
    return selected


def _panel_tables(
    effects: pd.DataFrame,
    spectra: pd.DataFrame,
    motifs: pd.DataFrame,
    day: str,
    programs: Sequence[int],
    targets: Sequence[str],
    genes: Sequence[str],
    motif_names: Sequence[str],
) -> dict[str, pd.DataFrame]:
    differential = effects[
        (effects.day == day) & effects.program.isin(programs)
    ]
    effect_table = differential.pivot_table(
        index="program", columns="target", values="effect", aggfunc="mean"
    ).reindex(index=programs, columns=targets)
    qvalue_table = differential.pivot_table(
        index="program", columns="target", values="q_value", aggfunc="min"
    ).reindex(index=programs, columns=targets)

    loading_table = spectra.reindex(index=programs, columns=genes)
    loading_table = loading_table.div(
        loading_table.max(axis=1).replace(0, np.nan), axis=0
    )

    motif_data = motifs[
        (motifs.day == day)
        & motifs.program.isin(programs)
        & motifs.motif.isin(motif_names)
    ].copy()
    motif_data["display"] = (
        np.sign(motif_data.score.fillna(1))
        * -np.log10(motif_data.q_value.clip(lower=np.nextafter(0, 1)))
    )
    motif_table = motif_data.pivot_table(
        index="program", columns="motif", values="display", aggfunc="max"
    ).reindex(index=programs, columns=motif_names)
    motif_qvalue_table = motif_data.pivot_table(
        index="program", columns="motif", values="q_value", aggfunc="min"
    ).reindex(index=programs, columns=motif_names)
    return {
        "effects": effect_table,
        "qvalues": qvalue_table,
        "loadings": loading_table,
        "motifs": motif_table,
        "motif_qvalues": motif_qvalue_table,
    }


def _draw_three(
    axes,
    tables: Mapping[str, pd.DataFrame],
    programs: Sequence[int],
    labels: Sequence[str],
    targets: Sequence[str],
    genes: Sequence[str],
    motif_names: Sequence[str],
    *,
    titles: bool = True,
):
    effects = tables["effects"]
    qvalues = tables["qvalues"]
    loadings = tables["loadings"]
    motif_values = tables["motifs"]
    motif_qvalues = tables["motif_qvalues"]

    effect_finite = np.abs(effects.to_numpy(float))
    effect_finite = effect_finite[np.isfinite(effect_finite)]
    vmax = (
        max(0.5, float(np.percentile(effect_finite, 98)))
        if effect_finite.size
        else 1.0
    )
    image = axes[0].imshow(
        effects,
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax),
        aspect="auto",
    )
    axes[0].set_xticks(
        np.arange(len(targets)), targets, rotation=60, ha="right", fontsize=6
    )
    axes[0].set_yticks(np.arange(len(programs)), labels, fontsize=7)
    for row, program in enumerate(programs):
        for column, target in enumerate(targets):
            star = significance_stars(qvalues.loc[program, target])
            if star:
                axes[0].text(
                    column, row, star, ha="center", va="center", fontsize=6
                )

    axes[1].imshow(loadings, cmap="Greys", vmin=0, vmax=1, aspect="auto")
    axes[1].set_xticks(
        np.arange(len(genes)), genes, rotation=60, ha="right", fontsize=6
    )
    axes[1].set_yticks([])

    motif_finite = np.abs(motif_values.to_numpy(float))
    motif_finite = motif_finite[np.isfinite(motif_finite)]
    motif_limit = (
        max(1.0, float(np.percentile(motif_finite, 98)))
        if motif_finite.size
        else 1.0
    )
    axes[2].imshow(
        motif_values,
        cmap="PuOr_r",
        norm=TwoSlopeNorm(
            vmin=-motif_limit, vcenter=0, vmax=motif_limit
        ),
        aspect="auto",
    )
    axes[2].set_xticks(
        np.arange(len(motif_names)),
        motif_names,
        rotation=60,
        ha="right",
        fontsize=6,
    )
    axes[2].set_yticks([])
    for row, program in enumerate(programs):
        for column, motif in enumerate(motif_names):
            star = (
                significance_stars(motif_qvalues.loc[program, motif])
                if motif in motif_qvalues.columns
                else ""
            )
            if star:
                axes[2].text(
                    column, row, star, ha="center", va="center", fontsize=6
                )

    if titles:
        axes[0].set_title(
            "Effect of gene CRISPRi\non cNMF program expression", fontsize=8
        )
        axes[1].set_title("cNMF gene loading", fontsize=8)
        axes[2].set_title("TF motif enrichment", fontsize=8)
    return image


def plot_integrated_annotation(
    effects: pd.DataFrame,
    spectra: pd.DataFrame,
    motifs: pd.DataFrame,
    *,
    day: str,
    programs: Sequence[int],
    targets: Sequence[str],
    genes: Sequence[str] | None = None,
    motif_names: Sequence[str] | None = None,
    program_labels: Mapping[int, str] | None = None,
):
    day = clean_day(day) or day
    programs = [int(program) for program in programs if int(program) in spectra.index]
    genes = list(genes) if genes else _top_genes(spectra, programs)
    motif_names = (
        list(motif_names)
        if motif_names
        else _top_motifs(motifs, day, programs)
    )
    labels = [
        (program_labels or {}).get(program, f"P{program}")
        for program in programs
    ]
    tables = _panel_tables(
        effects,
        spectra,
        motifs,
        day,
        programs,
        list(targets),
        genes,
        motif_names,
    )
    widths = [max(1, len(targets)), max(1, len(genes)), max(1, len(motif_names))]
    figure, axes = plt.subplots(
        1,
        3,
        figsize=(max(9.5, 0.35 * sum(widths)), max(2.3, 0.55 * len(programs))),
        gridspec_kw={"width_ratios": widths},
        constrained_layout=True,
    )
    image = _draw_three(
        axes,
        tables,
        programs,
        labels,
        list(targets),
        genes,
        motif_names,
    )
    figure.colorbar(
        image, ax=axes[0], fraction=0.04, pad=0.02, label="log effect"
    )
    return figure, tables


def plot_multi_day_annotation(
    effects: pd.DataFrame,
    spectra: pd.DataFrame,
    motifs: pd.DataFrame,
    *,
    programs_by_day: Mapping[str, Sequence[int]],
    targets_by_day: Mapping[str, Sequence[str]],
    genes_by_day: Mapping[str, Sequence[str]] | None = None,
    motifs_by_day: Mapping[str, Sequence[str]] | None = None,
    program_labels: Mapping[int, str] | None = None,
):
    days = [day for day in ("D0", "D1", "D2", "D3") if day in programs_by_day]
    figure = plt.figure(
        figsize=(13, max(7, 2.2 * len(days))), constrained_layout=True
    )
    outer = figure.add_gridspec(
        len(days),
        3,
        height_ratios=[max(1, len(programs_by_day[day])) for day in days],
    )
    output: dict[str, dict[str, pd.DataFrame]] = {}
    for row, day in enumerate(days):
        programs = [
            int(program)
            for program in programs_by_day[day]
            if int(program) in spectra.index
        ]
        targets = list(targets_by_day.get(day, []))
        genes = (
            list(genes_by_day.get(day, []))
            if genes_by_day and day in genes_by_day
            else _top_genes(spectra, programs)
        )
        motif_names = (
            list(motifs_by_day.get(day, []))
            if motifs_by_day and day in motifs_by_day
            else _top_motifs(motifs, day, programs)
        )
        tables = _panel_tables(
            effects,
            spectra,
            motifs,
            day,
            programs,
            targets,
            genes,
            motif_names,
        )
        output[day] = tables
        labels = [
            (program_labels or {}).get(program, f"P{program}")
            for program in programs
        ]
        axes = [figure.add_subplot(outer[row, column]) for column in range(3)]
        _draw_three(
            axes,
            tables,
            programs,
            labels,
            targets,
            genes,
            motif_names,
            titles=(row == 0),
        )
        axes[0].set_ylabel(
            day, rotation=0, labelpad=20, fontsize=8, fontweight="bold"
        )
    return figure, output
