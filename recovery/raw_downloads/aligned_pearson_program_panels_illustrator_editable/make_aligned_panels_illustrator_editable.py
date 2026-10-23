from __future__ import annotations

from pathlib import Path
import re
import subprocess

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    annotate,
    coord_fixed,
    element_blank,
    element_rect,
    element_text,
    geom_raster,
    geom_rect,
    geom_segment,
    geom_text,
    geom_tile,
    ggplot,
    scale_fill_cmap,
    theme,
    theme_void,
)
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

INPUT_ASSOC = Path('/mnt/data/cNMF_50_0_2_032526_6co(3).xlsx')
ASSOC_SHEET = 'Perturbation Association 1'
INPUT_EC = Path('/mnt/data/d3_EC_percent_by_guide_with_individual_NT_excl_Hippo(3).csv')
CACHED_EC_STATS = Path('/mnt/data/gene_program_EC_heatmap_plotnine_revised/all_298_gene_D3_EC_permutation_results.csv')
OUTDIR = Path('/mnt/data/aligned_pearson_program_panels_illustrator_editable')
OUTDIR.mkdir(parents=True, exist_ok=True)

CMAP_NAME = 'RdBu_r'
PROGRAM_COLOR_LIMIT = 3.0
# Liberation Sans is metrically compatible with Arial for rendering. The final SVG
# declarations are rewritten to Arial so text opens as ordinary editable text in Illustrator.
FONT_FAMILY = 'Liberation Sans'

# Order is shared by panel a insets and panel b rows.
GENE_GROUPS: list[tuple[str, list[str]]] = [
    ('BMP/TGFb pathway', ['BMPR1A', 'SMAD1', 'SMAD5', 'ACVR2A', 'SMAD4']),
    ('MIXL1/EOMES cluster', ['MIXL1', 'EOMES', 'ACVR1B', 'SMAD3']),
    ('EC+ LM regulators', ['KDR', 'KCTD10', 'PRKCD', 'ETV2', 'PLCG1']),
    ('Hippo pathway', ['TAOK1', 'AMOTL2', 'LATS2', 'NF2']),
]

PROGRAMS: list[dict[str, object]] = [
    {'day': 'D1', 'program': 13, 'label': 'P13 Posterior\nprimitive streak', 'x': 13.00},
    {'day': 'D1', 'program': 30, 'label': 'P30 Primitive streak,\nBMP-repressed', 'x': 16.20},
    {'day': 'D1', 'program': 48, 'label': 'P48 Primitive streak', 'x': 19.40},
    {'day': 'D2', 'program': 15, 'label': 'P15 Vascular\nprogenitor', 'x': 27.00},
    {'day': 'D2', 'program': 14, 'label': 'P14 Lateral\nmesoderm', 'x': 30.20},
    {'day': 'D2', 'program': 17, 'label': 'P17 Lateral mesoderm,\nYAP/TAZ-induced', 'x': 33.40},
    {'day': 'D3', 'program': 2, 'label': 'P2 Artery\nendothelium 1', 'x': 41.00},
    {'day': 'D3', 'program': 3, 'label': 'P3 Endothelial\nfunction 1', 'x': 44.20},
    {'day': 'D3', 'program': 8, 'label': 'P8 Artery\nendothelium 2', 'x': 47.40},
    {'day': 'D3', 'program': 11, 'label': 'P11 Endothelial\nfunction 2', 'x': 50.60},
    {'day': 'D3', 'program': 33, 'label': 'P33 Anterior-like\nmesoderm', 'x': 53.80},
    {'day': 'D3', 'program': 34, 'label': 'P34 Mesenchyme,\nYAP/TAZ-induced', 'x': 57.00},
    {'day': 'D3', 'program': 27, 'label': 'P27 Mesenchyme,\nECM', 'x': 60.20},
]

# Shared vertical coordinates ensure that the four gene blocks align across panels.
ROW_STEP = 1.65
BETWEEN_BLOCK_GAP = 1.80
TOP_GENE_Y = 49.00
TILE_WIDTH = 3.00
TILE_HEIGHT = 1.55

PANEL_A_WIDTH = 78.0
PANEL_B_WIDTH = 102.0
PANEL_HEIGHT = 60.0
PANEL_B_COMBINED_OFFSET = PANEL_A_WIDTH

# Panel-a geometry.
A_MAIN_X0 = 2.50
A_MAIN_Y0 = 11.50
A_MAIN_SIZE = 40.00
A_INSET_X0 = 55.00
A_INSET_CELL = ROW_STEP
A_GENE_LABEL_X = 63.80

# Panel-b geometry.
B_GENE_LABEL_X = 9.70
B_EC_X = 70.00


def significance_stars(fdr: float) -> str:
    if fdr < 0.001:
        return '***'
    if fdr < 0.01:
        return '**'
    if fdr < 0.05:
        return '*'
    return ''


def build_gene_y_positions() -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    y_by_gene: dict[str, float] = {}
    center_by_group: dict[str, float] = {}
    title_y_by_group: dict[str, float] = {}
    y = TOP_GENE_Y
    for group_index, (group_name, genes) in enumerate(GENE_GROUPS):
        ys: list[float] = []
        for gene in genes:
            y_by_gene[gene] = y
            ys.append(y)
            y -= ROW_STEP
        center_by_group[group_name] = float(np.mean(ys))
        title_y_by_group[group_name] = max(ys) + 1.70
        if group_index < len(GENE_GROUPS) - 1:
            y -= BETWEEN_BLOCK_GAP
    return y_by_gene, center_by_group, title_y_by_group


def load_association_data() -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray, list[str]]:
    assoc = pd.read_excel(INPUT_ASSOC, sheet_name=ASSOC_SHEET)
    required = {'target_name', 'program_name', 'approx_log2FC', 'adj_pval', 'Sample'}
    missing = required.difference(assoc.columns)
    if missing:
        raise ValueError(f'Missing columns in association sheet: {sorted(missing)}')

    assoc = assoc.loc[assoc['target_name'].ne('non-targeting')].copy()
    feature_order = pd.MultiIndex.from_product(
        [['D0', 'D1', 'D2', 'D3'], list(range(50))],
        names=['Sample', 'program_name'],
    )
    matrix = (
        assoc.pivot(index='target_name', columns=['Sample', 'program_name'], values='approx_log2FC')
        .reindex(columns=feature_order)
        .sort_index()
    )
    if matrix.shape != (298, 200):
        raise ValueError(f'Expected 298 x 200 association matrix; found {matrix.shape}')
    if matrix.isna().any().any():
        raise ValueError(f'Association matrix contains {int(matrix.isna().sum().sum())} missing values')

    corr = np.corrcoef(matrix.to_numpy(dtype=float))
    corr = np.clip(corr, -1.0, 1.0)
    np.fill_diagonal(corr, 1.0)
    distance = np.clip(1.0 - corr, 0.0, 2.0)
    np.fill_diagonal(distance, 0.0)
    tree = linkage(squareform(distance, checks=False), method='average', optimal_ordering=True)
    order = leaves_list(tree)[::-1]
    ordered_genes = matrix.index.to_numpy()[order].tolist()
    return assoc, matrix, corr, order, ordered_genes


def load_ec_statistics() -> pd.DataFrame:
    if not CACHED_EC_STATS.exists():
        raise FileNotFoundError(
            'The cached guide-level permutation results were not found. '
            'Regenerate them with the revised EC heatmap script first.'
        )
    stats = pd.read_csv(CACHED_EC_STATS)
    required = {
        'gene', 'log2_gene_mean_over_NT_mean',
        'permutation_pvalue', 'permutation_fdr_bh_298_genes', 'significance',
    }
    missing = required.difference(stats.columns)
    if missing:
        raise ValueError(f'Cached EC-statistics file lacks columns: {sorted(missing)}')
    return stats


def panel_a_layers(
    corr: np.ndarray,
    order: np.ndarray,
    ordered_genes: list[str],
    matrix_genes: list[str],
    y_by_gene: dict[str, float],
    center_by_group: dict[str, float],
    title_y_by_group: dict[str, float],
    x_offset: float = 0.0,
):
    layers = []
    ordered = corr[np.ix_(order, order)]
    n = ordered.shape[0]
    x = x_offset + A_MAIN_X0 + (np.arange(n) + 0.5) * A_MAIN_SIZE / n
    y = A_MAIN_Y0 + A_MAIN_SIZE - (np.arange(n) + 0.5) * A_MAIN_SIZE / n
    xx, yy = np.meshgrid(x, y)
    main_df = pd.DataFrame({'x': xx.ravel(), 'y': yy.ravel(), 'fill_value': ordered.ravel()})
    layers.append(geom_raster(main_df, aes('x', 'y', fill='fill_value'), inherit_aes=False))
    layers.append(
        geom_rect(
            pd.DataFrame([{
                'xmin': x_offset + A_MAIN_X0,
                'xmax': x_offset + A_MAIN_X0 + A_MAIN_SIZE,
                'ymin': A_MAIN_Y0,
                'ymax': A_MAIN_Y0 + A_MAIN_SIZE,
            }]),
            aes(xmin='xmin', xmax='xmax', ymin='ymin', ymax='ymax'),
            inherit_aes=False, fill=None, color='#B8B8B8', size=0.30,
        )
    )

    gene_to_index = {gene: i for i, gene in enumerate(matrix_genes)}
    rank = {gene: i for i, gene in enumerate(ordered_genes)}
    inset_rows: list[dict[str, object]] = []
    box_rows: list[dict[str, float | str]] = []
    segment_rows: list[dict[str, float | str]] = []

    def add_polyline(group: str, points: list[tuple[float, float]]) -> None:
        for (x1, y1), (x2, y2) in zip(points[:-1], points[1:]):
            segment_rows.append({
                'group': group,
                'x': x1,
                'y': y1,
                'xend': x2,
                'yend': y2,
            })

    elbow_1 = {
        'BMP/TGFb pathway': x_offset + 46.2,
        'MIXL1/EOMES cluster': x_offset + 46.8,
        'EC+ LM regulators': x_offset + 47.4,
        'Hippo pathway': x_offset + 48.0,
    }
    elbow_2 = {
        'BMP/TGFb pathway': x_offset + 52.1,
        'MIXL1/EOMES cluster': x_offset + 52.3,
        'EC+ LM regulators': x_offset + 52.5,
        'Hippo pathway': x_offset + 52.7,
    }

    for group_name, genes in GENE_GROUPS:
        idx = [gene_to_index[g] for g in genes]
        sub = corr[np.ix_(idx, idx)]
        for row_idx, row_gene in enumerate(genes):
            for col_idx, col_gene in enumerate(genes):
                inset_rows.append({
                    'group': group_name,
                    'row_gene': row_gene,
                    'col_gene': col_gene,
                    'x': x_offset + A_INSET_X0 + (col_idx + 0.5) * A_INSET_CELL,
                    'y': y_by_gene[row_gene],
                    'fill_value': float(sub[row_idx, col_idx]),
                })

        ranks = [rank[g] for g in genes]
        start = min(ranks)
        end = max(ranks)
        cell = A_MAIN_SIZE / n
        xmin = x_offset + A_MAIN_X0 + start * cell
        xmax = x_offset + A_MAIN_X0 + (end + 1) * cell
        ymax = A_MAIN_Y0 + A_MAIN_SIZE - start * cell
        ymin = A_MAIN_Y0 + A_MAIN_SIZE - (end + 1) * cell
        box_rows.append({
            'group': group_name,
            'xmin': xmin,
            'xmax': xmax,
            'ymin': ymin,
            'ymax': ymax,
        })
        anchor_y = (ymin + ymax) / 2
        target_y = center_by_group[group_name]
        add_polyline(
            group_name,
            [
                (xmax, anchor_y),
                (elbow_1[group_name], anchor_y),
                (elbow_2[group_name], target_y),
                (x_offset + A_INSET_X0 - 0.75, target_y),
            ],
        )

    inset_df = pd.DataFrame(inset_rows)
    box_df = pd.DataFrame(box_rows)
    segment_df = pd.DataFrame(segment_rows)
    layers.append(
        geom_tile(
            inset_df,
            aes('x', 'y', fill='fill_value'),
            inherit_aes=False,
            width=A_INSET_CELL,
            height=A_INSET_CELL,
            color='white',
            size=0.35,
        )
    )
    layers.append(
        geom_segment(
            segment_df,
            aes(x='x', y='y', xend='xend', yend='yend'),
            inherit_aes=False,
            color='black',
            size=0.48,
            lineend='butt',
        )
    )
    layers.append(
        geom_rect(
            box_df,
            aes(xmin='xmin', xmax='xmax', ymin='ymin', ymax='ymax'),
            inherit_aes=False,
            fill=None,
            color='#555555',
            size=0.56,
        )
    )

    # Pearson color bar, normalized directly to [-1, 1].
    cb_values = np.linspace(-1.0, 1.0, 181)
    cb_x0, cb_x1 = x_offset + 54.0, x_offset + 66.0
    cb_step = (cb_x1 - cb_x0) / len(cb_values)
    cb_df = pd.DataFrame({
        'x': cb_x0 + (np.arange(len(cb_values)) + 0.5) * cb_step,
        'y': np.full(len(cb_values), 8.15),
        'fill_value': cb_values,
    })
    layers.append(
        geom_tile(
            cb_df,
            aes('x', 'y', fill='fill_value'),
            inherit_aes=False,
            width=cb_step + 0.006,
            height=1.35,
            color=None,
        )
    )

    layers.extend([
        annotate('text', x=x_offset + 1.20, y=57.25, label='a', family=FONT_FAMILY, fontweight='bold', size=28, ha='center', va='center'),
        annotate('text', x=x_offset + 22.50, y=54.35, label='298 perturbed genes', family=FONT_FAMILY, size=17.0, ha='center', va='center'),
        annotate('text', x=x_offset + 22.50, y=8.60, label='Pearson r - Effect on program expression', family=FONT_FAMILY, size=10.2, ha='center', va='center'),
        annotate('text', x=cb_x0, y=6.80, label='-1', family=FONT_FAMILY, size=7.8, ha='center', va='center'),
        annotate('text', x=(cb_x0 + cb_x1) / 2, y=6.80, label='0', family=FONT_FAMILY, size=7.8, ha='center', va='center'),
        annotate('text', x=cb_x1, y=6.80, label='1', family=FONT_FAMILY, size=7.8, ha='center', va='center'),
        annotate('text', x=(cb_x0 + cb_x1) / 2, y=4.80, label='Pearson r', family=FONT_FAMILY, size=13.2, ha='center', va='center'),
    ])

    for group_name, genes in GENE_GROUPS:
        title_size = 12.4 if group_name != 'MIXL1/EOMES cluster' else 11.2
        layers.append(
            annotate(
                'text',
                x=x_offset + A_INSET_X0 + len(genes) * A_INSET_CELL / 2,
                y=title_y_by_group[group_name],
                label=group_name,
                family=FONT_FAMILY,
                size=title_size,
                ha='center',
                va='center',
            )
        )
        for gene in genes:
            layers.append(
                annotate(
                    'text',
                    x=x_offset + A_GENE_LABEL_X,
                    y=y_by_gene[gene],
                    label=gene,
                    family=FONT_FAMILY,
                    size=10.2,
                    ha='left',
                    va='center',
                )
            )
    return layers, inset_df, box_df


def panel_b_layers(
    assoc: pd.DataFrame,
    ec_stats: pd.DataFrame,
    y_by_gene: dict[str, float],
    x_offset: float = 0.0,
):
    layers = []
    selected_genes = [gene for _, genes in GENE_GROUPS for gene in genes]
    selected_keys = {(str(p['day']), int(p['program'])) for p in PROGRAMS}
    assoc_selected = assoc.loc[
        assoc['target_name'].isin(selected_genes)
        & assoc.apply(lambda row: (str(row['Sample']), int(row['program_name'])) in selected_keys, axis=1)
    ].copy()
    lookup = {
        (str(row.target_name), str(row.Sample), int(row.program_name)): (float(row.approx_log2FC), float(row.adj_pval))
        for row in assoc_selected.itertuples(index=False)
    }

    tile_rows: list[dict[str, object]] = []
    star_rows: list[dict[str, object]] = []
    for _, genes in GENE_GROUPS:
        for gene in genes:
            for p in PROGRAMS:
                key = (gene, str(p['day']), int(p['program']))
                if key not in lookup:
                    raise KeyError(f'Missing program association: {key}')
                effect, fdr = lookup[key]
                stars = significance_stars(fdr)
                tile_rows.append({
                    'gene': gene,
                    'day': p['day'],
                    'program': p['program'],
                    'x': x_offset + float(p['x']),
                    'y': y_by_gene[gene],
                    # Normalize [-3, +3] onto the shared [-1, +1] color scale.
                    'fill_value': float(np.clip(effect, -PROGRAM_COLOR_LIMIT, PROGRAM_COLOR_LIMIT) / PROGRAM_COLOR_LIMIT),
                    'approx_log2FC': effect,
                    'adj_pval': fdr,
                    'significance': stars,
                })
                if stars:
                    star_rows.append({
                        'x': x_offset + float(p['x']),
                        'y': y_by_gene[gene],
                        'label': stars,
                        'text_color': 'white' if abs(effect) >= 1.45 else 'black',
                    })

    ec_lookup = ec_stats.set_index('gene')
    ec_rows: list[dict[str, object]] = []
    for gene in selected_genes:
        if gene not in ec_lookup.index:
            raise KeyError(f'Missing EC statistics for {gene}')
        row = ec_lookup.loc[gene]
        effect = float(row['log2_gene_mean_over_NT_mean'])
        fdr = float(row['permutation_fdr_bh_298_genes'])
        stars = significance_stars(fdr)
        ec_rows.append({
            'gene': gene,
            'x': x_offset + B_EC_X,
            'y': y_by_gene[gene],
            'fill_value': float(np.clip(effect, -PROGRAM_COLOR_LIMIT, PROGRAM_COLOR_LIMIT) / PROGRAM_COLOR_LIMIT),
            'ec_log2_ratio': effect,
            'permutation_pvalue': float(row['permutation_pvalue']),
            'permutation_fdr': fdr,
            'significance': stars,
        })
        if stars:
            star_rows.append({
                'x': x_offset + B_EC_X,
                'y': y_by_gene[gene],
                'label': stars,
                'text_color': 'white' if abs(effect) >= 1.45 else 'black',
            })

    tile_df = pd.DataFrame(tile_rows)
    ec_df = pd.DataFrame(ec_rows)
    stars_df = pd.DataFrame(star_rows)
    white_stars = stars_df.loc[stars_df['text_color'].eq('white')]
    black_stars = stars_df.loc[stars_df['text_color'].eq('black')]

    layers.extend([
        geom_tile(
            tile_df,
            aes('x', 'y', fill='fill_value'),
            inherit_aes=False,
            width=TILE_WIDTH,
            height=TILE_HEIGHT,
            color='white',
            size=0.35,
        ),
        geom_tile(
            ec_df,
            aes('x', 'y', fill='fill_value'),
            inherit_aes=False,
            width=TILE_WIDTH,
            height=TILE_HEIGHT,
            color='white',
            size=0.35,
        ),
        geom_text(
            white_stars,
            aes('x', 'y', label='label'),
            inherit_aes=False,
            color='white',
            family=FONT_FAMILY,
            size=5.6,
            fontweight='bold',
            va='center',
        ),
        geom_text(
            black_stars,
            aes('x', 'y', label='label'),
            inherit_aes=False,
            color='black',
            family=FONT_FAMILY,
            size=5.6,
            fontweight='bold',
            va='center',
        ),
    ])

    # Program/EC colorbar. Values are normalized to the shared [-1, +1] fill scale.
    cb_values = np.linspace(-1.0, 1.0, 241)
    cb_x0, cb_x1 = x_offset + 77.0, x_offset + 98.0
    cb_step = (cb_x1 - cb_x0) / len(cb_values)
    cb_df = pd.DataFrame({
        'x': cb_x0 + (np.arange(len(cb_values)) + 0.5) * cb_step,
        'y': np.full(len(cb_values), 41.10),
        'fill_value': cb_values,
    })
    layers.append(
        geom_tile(
            cb_df,
            aes('x', 'y', fill='fill_value'),
            inherit_aes=False,
            width=cb_step + 0.004,
            height=1.30,
            color=None,
        )
    )

    # Titles and guide lines.
    lines = pd.DataFrame([
        {'x': x_offset + 11.30, 'xend': x_offset + 21.10, 'y': 52.40, 'yend': 52.40},
        {'x': x_offset + 25.30, 'xend': x_offset + 35.10, 'y': 52.40, 'yend': 52.40},
        {'x': x_offset + 39.30, 'xend': x_offset + 61.90, 'y': 52.40, 'yend': 52.40},
        {'x': x_offset + 68.30, 'xend': x_offset + 71.70, 'y': 52.40, 'yend': 52.40},
    ])
    layers.append(
        geom_segment(
            lines,
            aes(x='x', y='y', xend='xend', yend='yend'),
            inherit_aes=False,
            color='#5A5A5A',
            size=0.48,
        )
    )
    layers.extend([
        annotate('text', x=x_offset + 1.35, y=56.55, label='b', family=FONT_FAMILY, fontweight='bold', size=28, ha='center', va='center'),
        annotate('text', x=x_offset + 36.50, y=57.20, label='Effects on gene programs', family=FONT_FAMILY, fontweight='bold', size=15.5, ha='center', va='center'),
        annotate('text', x=x_offset + B_EC_X, y=57.25, label='Effects on EC\ndifferentiation', family=FONT_FAMILY, fontweight='bold', size=15.0, ha='center', va='center'),
        annotate('text', x=x_offset + 16.20, y=54.00, label='Day 1', family=FONT_FAMILY, fontweight='bold', size=10.5, ha='center', va='center'),
        annotate('text', x=x_offset + 30.20, y=54.00, label='Day 2', family=FONT_FAMILY, fontweight='bold', size=10.5, ha='center', va='center'),
        annotate('text', x=x_offset + 50.60, y=54.00, label='Day 3', family=FONT_FAMILY, fontweight='bold', size=10.5, ha='center', va='center'),
    ])

    for gene in selected_genes:
        layers.append(
            annotate(
                'text',
                x=x_offset + B_GENE_LABEL_X,
                y=y_by_gene[gene],
                label=gene,
                family=FONT_FAMILY,
                size=10.4,
                ha='right',
                va='center',
            )
        )

    for p in PROGRAMS:
        layers.append(
            annotate(
                'text',
                x=x_offset + float(p['x']) + 0.20,
                y=14.25,
                label=str(p['label']),
                family=FONT_FAMILY,
                size=8.0,
                angle=55,
                ha='right',
                va='top',
            )
        )
    layers.append(
        annotate(
            'text',
            x=x_offset + B_EC_X + 0.22,
            y=14.25,
            label='D3 EC\nlog2 ratio vs NT',
            family=FONT_FAMILY,
            size=8.6,
            angle=55,
            ha='right',
            va='top',
        )
    )

    # Right-side legend.
    layers.extend([
        annotate('text', x=x_offset + 77.0, y=45.00, label='Effect / EC ratio', family=FONT_FAMILY, fontweight='bold', size=10.8, ha='left', va='center'),
        annotate('text', x=x_offset + 77.0, y=43.75, label='log2 scale', family=FONT_FAMILY, size=9.7, ha='left', va='center'),
        annotate('text', x=cb_x0, y=39.85, label='-3', family=FONT_FAMILY, size=7.2, ha='center', va='center'),
        annotate('text', x=(cb_x0 + cb_x1) / 2, y=39.85, label='0', family=FONT_FAMILY, size=7.2, ha='center', va='center'),
        annotate('text', x=cb_x1, y=39.85, label='+3', family=FONT_FAMILY, size=7.2, ha='center', va='center'),
        annotate('text', x=x_offset + 77.0, y=36.80, label='Significance', family=FONT_FAMILY, fontweight='bold', size=10.8, ha='left', va='center'),
        annotate('text', x=x_offset + 77.0, y=35.45, label='* FDR < 0.05\n** FDR < 0.01\n*** FDR < 0.001', family=FONT_FAMILY, size=10.2, ha='left', va='top'),
        annotate('text', x=x_offset + 77.0, y=29.10, label='D3 EC ratio', family=FONT_FAMILY, fontweight='bold', size=10.5, ha='left', va='center'),
        annotate('text', x=x_offset + 77.0, y=27.80, label='log2(gene mean / NT mean)', family=FONT_FAMILY, size=9.5, ha='left', va='center'),
    ])
    return layers, tile_df, ec_df


def base_plot(xlim: tuple[float, float], figure_size: tuple[float, float]):
    return (
        ggplot()
        + scale_fill_cmap(cmap_name=CMAP_NAME, limits=(-1.0, 1.0), guide=None)
        + coord_fixed(ratio=1.0, xlim=xlim, ylim=(0.0, PANEL_HEIGHT), expand=False)
        + theme_void()
        + theme(
            figure_size=figure_size,
            plot_background=element_rect(fill='white', color='white'),
            panel_background=element_rect(fill='white', color='white'),
            text=element_text(family=FONT_FAMILY, color='black'),
            axis_title=element_blank(),
            axis_text=element_blank(),
            axis_ticks=element_blank(),
            axis_line=element_blank(),
            plot_margin=0,
        )
    )


def save_plot(plot, path_prefix: Path, figure_size: tuple[float, float]) -> None:
    """Save an Illustrator-oriented SVG/PDF with live text.

    Matplotlib writes text as SVG <text> elements because svg.fonttype is set to
    "none". The SVG font declaration is then changed to Arial (metric-compatible
    with Liberation Sans), and Inkscape converts that SVG to PDF without outlining
    text. The dense 298 x 298 Pearson matrix remains one embedded raster image,
    while labels, stars, titles, legends, inset heatmaps, boxes, and connectors remain
    vector/editable.
    """
    fig = plot.draw(show=False)
    fig.set_size_inches(*figure_size, forward=True)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    svg_path = path_prefix.with_suffix('.svg')
    pdf_path = path_prefix.with_suffix('.pdf')
    png_path = path_prefix.with_suffix('.png')
    temp_svg = path_prefix.with_name(path_prefix.name + '_temp').with_suffix('.svg')

    fig.savefig(temp_svg, format='svg', dpi=300, facecolor='white')
    fig.savefig(png_path, format='png', dpi=300, facecolor='white')
    plt.close(fig)

    svg_text = temp_svg.read_text(encoding='utf-8')
    # Keep all text live and request Arial on the user's system.
    svg_text = svg_text.replace("font-family: 'Liberation Sans'", 'font-family: Arial')
    svg_text = svg_text.replace('font-family: Liberation Sans', 'font-family: Arial')
    svg_path.write_text(svg_text, encoding='utf-8')
    temp_svg.unlink(missing_ok=True)

    # Inkscape preserves the SVG text objects as embedded TrueType text in PDF.
    subprocess.run([
        'inkscape', str(svg_path),
        '--export-type=pdf',
        f'--export-filename={pdf_path}',
    ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def main() -> None:
    mpl.rcParams.update({
        'font.family': FONT_FAMILY,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
        'text.usetex': False,
    })

    y_by_gene, center_by_group, title_y_by_group = build_gene_y_positions()
    assoc, matrix, corr, order, ordered_genes = load_association_data()
    ec_stats = load_ec_statistics()

    # Panel a.
    plot_a = base_plot((0.0, PANEL_A_WIDTH), (7.8, 6.0))
    a_layers, inset_values, anchor_boxes = panel_a_layers(
        corr, order, ordered_genes, matrix.index.tolist(),
        y_by_gene, center_by_group, title_y_by_group,
        x_offset=0.0,
    )
    for layer in a_layers:
        plot_a = plot_a + layer
    save_plot(plot_a, OUTDIR / 'pearson_correlation_clusters_illustrator_editable', (7.8, 6.0))

    # Panel b.
    plot_b = base_plot((0.0, PANEL_B_WIDTH), (10.2, 6.0))
    b_layers, program_values, ec_values = panel_b_layers(assoc, ec_stats, y_by_gene, x_offset=0.0)
    for layer in b_layers:
        plot_b = plot_b + layer
    save_plot(plot_b, OUTDIR / 'gene_program_and_EC_effects_heatmap_illustrator_editable', (10.2, 6.0))

    # Combined aligned figure.
    plot_ab = base_plot((0.0, PANEL_A_WIDTH + PANEL_B_WIDTH), (18.0, 6.0))
    a_layers_ab, _, _ = panel_a_layers(
        corr, order, ordered_genes, matrix.index.tolist(),
        y_by_gene, center_by_group, title_y_by_group,
        x_offset=0.0,
    )
    b_layers_ab, _, _ = panel_b_layers(
        assoc, ec_stats, y_by_gene,
        x_offset=PANEL_B_COMBINED_OFFSET,
    )
    for layer in a_layers_ab + b_layers_ab:
        plot_ab = plot_ab + layer
    save_plot(plot_ab, OUTDIR / 'aligned_panels_a_b_illustrator_editable', (18.0, 6.0))

    # Supporting values.
    pd.DataFrame({
        'cluster_rank': np.arange(1, len(ordered_genes) + 1),
        'target_name': ordered_genes,
    }).to_csv(OUTDIR / 'clustered_gene_order.csv', index=False)
    inset_values.to_csv(OUTDIR / 'pearson_inset_values.csv', index=False)
    anchor_boxes.to_csv(OUTDIR / 'pearson_anchor_boxes.csv', index=False)
    program_values.to_csv(OUTDIR / 'selected_gene_program_effects.csv', index=False)
    ec_values.to_csv(OUTDIR / 'selected_gene_EC_effects.csv', index=False)

    with (OUTDIR / 'METHOD.txt').open('w') as handle:
        handle.write(
            'Panel a: pairwise Pearson correlations among 298 perturbation effect vectors.\n'
            'Each vector contains approx_log2FC across 50 cNMF programs and D0-D3 (200 features).\n'
            'Ordering: average-linkage clustering on distance 1 - Pearson r with optimal leaf ordering.\n'
            'Panel b program effects: approx_log2FC from Perturbation Association 1.\n'
            'Program stars: Benjamini-Hochberg adjusted p-values from that sheet.\n'
            'EC effects: log2 of the mean EC percentage across the six targeting guides divided by the mean across usable individual NT guides.\n'
            'Cells carrying both the query guide and an NTC guide were excluded before guide-level EC percentages were calculated.\n'
            'EC stars: two-sided guide-level permutation test with BH adjustment across all 298 genes, using the previously generated all-gene result table.\n'
            'Shared gene-block order: BMP/TGFb pathway; MIXL1/EOMES cluster; EC+ LM regulators; Hippo pathway.\n'
            'Illustrator export: SVG text is retained as live <text> elements with Arial declarations; PDF is converted from SVG by Inkscape without text-to-path conversion.\n'
            'The dense 298 x 298 Pearson matrix is embedded as a raster image; all text and other figure elements remain editable/vector.\n'
        )

    print('Created outputs in', OUTDIR)


if __name__ == '__main__':
    main()
