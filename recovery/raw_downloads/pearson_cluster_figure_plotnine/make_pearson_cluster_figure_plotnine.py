from __future__ import annotations

from pathlib import Path
from functools import lru_cache

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    annotate,
    coord_fixed,
    geom_raster,
    geom_rect,
    geom_segment,
    geom_tile,
    ggplot,
    scale_fill_cmap,
    theme,
    theme_void,
    element_blank,
    element_rect,
    element_text,
)
from scipy.cluster.hierarchy import linkage, leaves_list, to_tree
from scipy.spatial.distance import squareform

INPUT_XLSX = Path('/mnt/data/cNMF_50_0_2_032526_6co(3).xlsx')
SHEET_NAME = 'Perturbation Association 1'
OUTDIR = Path('/mnt/data/pearson_cluster_figure_plotnine')
OUTDIR.mkdir(parents=True, exist_ok=True)

CMAP_NAME = 'RdBu_r'

BASE_CLUSTERS: dict[str, list[str]] = {
    'BMP/TGFb pathway': ['BMPR1A', 'SMAD1', 'SMAD5', 'ACVR2A', 'SMAD4'],
    'EC+ LM regulators': ['KDR', 'KCTD10', 'PRKCD', 'ETV2', 'PLCG1'],
    'Hippo pathway': ['TAOK1', 'AMOTL2', 'LATS2', 'NF2'],
}

CANVAS_X = (0.0, 100.0)
CANVAS_Y = (0.0, 68.0)
MAIN_X0 = 2.8
MAIN_Y0 = 11.0
MAIN_SIZE = 51.0

INSET_LAYOUT = {
    'BMP/TGFb pathway': dict(x0=69.0, y0=53.0, cell=2.10, title_x=74.25, title_y=65.2, genes_x=88.0),
    'EC+ LM regulators': dict(x0=68.4, y0=39.6, cell=2.10, title_x=73.65, title_y=51.1, genes_x=88.0),
    'Hippo pathway': dict(x0=69.2, y0=28.7, cell=2.10, title_x=73.40, title_y=38.2, genes_x=88.0),
    'MIXL1/EOMES cluster': dict(x0=68.0, y0=18.0, cell=1.95, title_x=74.4, title_y=25.2, genes_x=88.0),
}


def load_effect_matrix() -> pd.DataFrame:
    assoc = pd.read_excel(INPUT_XLSX, sheet_name=SHEET_NAME)
    required = {'target_name', 'program_name', 'approx_log2FC', 'Sample'}
    missing = required.difference(assoc.columns)
    if missing:
        raise ValueError(f'Missing required columns: {sorted(missing)}')

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
        raise ValueError(f'Expected a 298 x 200 effect matrix, found {matrix.shape}')
    if matrix.isna().any().any():
        raise ValueError(f'Effect matrix contains {int(matrix.isna().sum().sum())} missing values')

    requested = {g for genes in BASE_CLUSTERS.values() for g in genes}
    requested.update({'MIXL1', 'EOMES'})
    absent = sorted(requested.difference(matrix.index))
    if absent:
        raise ValueError(f'Requested genes are absent from the worksheet: {absent}')
    return matrix


def pearson_and_order(matrix: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, list[str], np.ndarray]:
    corr = np.corrcoef(matrix.to_numpy(dtype=float))
    corr = np.clip(corr, -1.0, 1.0)
    np.fill_diagonal(corr, 1.0)

    distance = np.clip(1.0 - corr, 0.0, 2.0)
    np.fill_diagonal(distance, 0.0)
    tree = linkage(squareform(distance, checks=False), method='average', optimal_ordering=True)

    order = leaves_list(tree)[::-1]
    ordered_genes = matrix.index.to_numpy()[order].tolist()
    return corr, order, ordered_genes, tree


def smallest_cluster_containing(tree: np.ndarray, gene_names: list[str], required: list[str], ordered_genes: list[str]) -> list[str]:
    root, nodes = to_tree(tree, rd=True)
    name_to_id = {name: i for i, name in enumerate(gene_names)}
    required_ids = {name_to_id[g] for g in required}

    @lru_cache(None)
    def leaf_ids(node_id: int) -> tuple[int, ...]:
        node = nodes[node_id]
        if node.is_leaf():
            return (node.id,)
        return leaf_ids(node.left.id) + leaf_ids(node.right.id)

    candidates = []
    for node in nodes:
        ls = set(leaf_ids(node.id))
        if required_ids.issubset(ls):
            candidates.append((len(ls), node.dist, node.id))
    _, _, best_id = sorted(candidates)[0]
    members = {gene_names[i] for i in leaf_ids(best_id)}
    return [g for g in ordered_genes if g in members]


def get_clusters(matrix: pd.DataFrame, tree: np.ndarray, ordered_genes: list[str]) -> dict[str, list[str]]:
    mixl1_cluster = smallest_cluster_containing(tree, matrix.index.tolist(), ['MIXL1', 'EOMES'], ordered_genes)
    mixl1_cluster = [g for g in mixl1_cluster if g != 'SPIN1']
    clusters = dict(BASE_CLUSTERS)
    clusters['MIXL1/EOMES cluster'] = mixl1_cluster
    return clusters


def main_heatmap_df(corr: np.ndarray, order: np.ndarray) -> pd.DataFrame:
    ordered = corr[np.ix_(order, order)]
    n = ordered.shape[0]
    x = MAIN_X0 + (np.arange(n) + 0.5) * MAIN_SIZE / n
    y = MAIN_Y0 + MAIN_SIZE - (np.arange(n) + 0.5) * MAIN_SIZE / n
    xx, yy = np.meshgrid(x, y)
    return pd.DataFrame({'x': xx.ravel(), 'y': yy.ravel(), 'pearson_r': ordered.ravel()})


def inset_df(title: str, genes: list[str], corr: np.ndarray, gene_to_index: dict[str, int]) -> pd.DataFrame:
    layout = INSET_LAYOUT[title]
    x0, y0, cell = layout['x0'], layout['y0'], layout['cell']
    idx = [gene_to_index[g] for g in genes]
    values = corr[np.ix_(idx, idx)]
    rows = []
    n = len(genes)
    for row in range(n):
        for col in range(n):
            rows.append({'cluster': title, 'x': x0 + (col + 0.5) * cell, 'y': y0 + (n - row - 0.5) * cell, 'pearson_r': float(values[row, col])})
    return pd.DataFrame(rows)


def rank_span_box(ranks: list[int], n: int) -> tuple[float, float, float, float, float, float]:
    start = min(ranks)
    end = max(ranks)
    cell = MAIN_SIZE / n
    xmin = MAIN_X0 + start * cell
    xmax = MAIN_X0 + (end + 1) * cell
    ymax = MAIN_Y0 + MAIN_SIZE - start * cell
    ymin = MAIN_Y0 + MAIN_SIZE - (end + 1) * cell
    xmid = xmax
    ymid = (ymin + ymax) / 2
    return xmin, xmax, ymin, ymax, xmid, ymid


def connector_and_anchor_data(clusters: dict[str, list[str]], ordered_genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rank = {g: i for i, g in enumerate(ordered_genes)}
    segments = []
    boxes = []
    n = len(ordered_genes)

    def add_polyline(title: str, points: list[tuple[float, float]]) -> None:
        for (x1, y1), (x2, y2) in zip(points[:-1], points[1:]):
            segments.append({'cluster': title, 'x': x1, 'y': y1, 'xend': x2, 'yend': y2})

    mids = {}
    for title, genes in clusters.items():
        l = INSET_LAYOUT[title]
        size = len(genes) * l['cell']
        mids[title] = (l['x0'], l['y0'] + size / 2)

    elbow_x = {
        'BMP/TGFb pathway': 59.3,
        'EC+ LM regulators': 57.9,
        'Hippo pathway': 58.7,
        'MIXL1/EOMES cluster': 60.0,
    }
    elbow2_x = {
        'BMP/TGFb pathway': 66.2,
        'EC+ LM regulators': 65.8,
        'Hippo pathway': 66.5,
        'MIXL1/EOMES cluster': 66.8,
    }

    for title, genes in clusters.items():
        ranks = [rank[g] for g in genes]
        xmin, xmax, ymin, ymax, xmid, ymid = rank_span_box(ranks, n)
        boxes.append({'cluster': title, 'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax})
        left, midy = mids[title]
        add_polyline(title, [(xmid, ymid), (elbow_x[title], ymid), (elbow2_x[title], midy), (left - 1.3, midy)])

    return pd.DataFrame(segments), pd.DataFrame(boxes)


def colorbar_df() -> pd.DataFrame:
    x0, x1 = 65.0, 83.0
    values = np.linspace(-1.0, 1.0, 181)
    step = (x1 - x0) / len(values)
    return pd.DataFrame({'x': x0 + (np.arange(len(values)) + 0.5) * step, 'y': np.full(len(values), 12.45), 'pearson_r': values})


def build_plot() -> tuple[ggplot, pd.DataFrame, pd.DataFrame]:
    matrix = load_effect_matrix()
    corr, order, ordered_genes, tree = pearson_and_order(matrix)
    clusters = get_clusters(matrix, tree, ordered_genes)
    main_df = main_heatmap_df(corr, order)
    gene_to_index = {g: i for i, g in enumerate(matrix.index)}
    inset_frames = {title: inset_df(title, genes, corr, gene_to_index) for title, genes in clusters.items()}
    segments, anchors = connector_and_anchor_data(clusters, ordered_genes)
    cb = colorbar_df()

    p = (
        ggplot()
        + geom_raster(main_df, aes('x', 'y', fill='pearson_r'), inherit_aes=False)
        + geom_rect(
            pd.DataFrame([{'xmin': MAIN_X0, 'xmax': MAIN_X0 + MAIN_SIZE, 'ymin': MAIN_Y0, 'ymax': MAIN_Y0 + MAIN_SIZE}]),
            aes(xmin='xmin', xmax='xmax', ymin='ymin', ymax='ymax'),
            inherit_aes=False, fill=None, color='#C7C7C7', size=0.35,
        )
    )

    for title, frame in inset_frames.items():
        cell = INSET_LAYOUT[title]['cell']
        p = p + geom_tile(frame, aes('x', 'y', fill='pearson_r'), inherit_aes=False, width=cell, height=cell, color='white', size=0.45)

    p = (
        p
        + geom_segment(segments, aes(x='x', y='y', xend='xend', yend='yend'), inherit_aes=False, color='black', size=0.52, lineend='butt')
        + geom_rect(anchors, aes(xmin='xmin', xmax='xmax', ymin='ymin', ymax='ymax'), inherit_aes=False, fill=None, color='#575757', size=0.65)
        + geom_tile(cb, aes('x', 'y', fill='pearson_r'), inherit_aes=False, width=(83.0 - 65.0) / len(cb) + 0.01, height=2.10, color=None)
        + scale_fill_cmap(cmap_name=CMAP_NAME, limits=(-1, 1), guide=None)
        + coord_fixed(ratio=1, xlim=CANVAS_X, ylim=CANVAS_Y, expand=False)
        + theme_void()
        + theme(
            figure_size=(12.0, 8.16),
            plot_background=element_rect(fill='white', color='white'),
            panel_background=element_rect(fill='white', color='white'),
            text=element_text(family='Arimo', color='black'),
            axis_title=element_blank(),
            axis_text=element_blank(),
            axis_ticks=element_blank(),
            axis_line=element_blank(),
            panel_grid=element_blank(),
            plot_margin=0,
        )
        + annotate('text', x=1.65, y=66.0, label='a', family='Arimo', fontweight='bold', size=31, ha='center', va='center')
        + annotate('text', x=28.3, y=65.4, label='298 perturbed genes', family='Arimo', size=23, ha='center', va='center')
        + annotate('text', x=28.3, y=8.05, label='Pearson r - Effect on program expression', family='Arimo', size=15.5, ha='center', va='center')
        + annotate('text', x=65.0, y=10.45, label='-1', family='Arimo', size=11.0, ha='center', va='center')
        + annotate('text', x=74.0, y=10.45, label='0', family='Arimo', size=11.0, ha='center', va='center')
        + annotate('text', x=83.0, y=10.45, label='1', family='Arimo', size=11.0, ha='center', va='center')
        + annotate('text', x=74.0, y=7.25, label='Pearson r', family='Arimo', size=20, ha='center', va='center')
    )

    for title, genes in clusters.items():
        l = INSET_LAYOUT[title]
        n = len(genes)
        center_y = l['y0'] + n * l['cell'] / 2
        p = p + annotate('text', x=l['title_x'], y=l['title_y'], label=title, family='Arimo', size=16.5, ha='center', va='center')
        p = p + annotate('text', x=l['genes_x'], y=center_y, label='\n'.join(genes), family='Arimo', size=13.0, ha='left', va='center')

    rank_table = pd.DataFrame({'cluster_rank': np.arange(1, len(ordered_genes) + 1), 'target_name': ordered_genes})
    pair_rows = []
    for title, genes in clusters.items():
        idx = [gene_to_index[g] for g in genes]
        sub = corr[np.ix_(idx, idx)]
        for i, row_gene in enumerate(genes):
            for j, col_gene in enumerate(genes):
                pair_rows.append({'cluster': title, 'row_gene': row_gene, 'column_gene': col_gene, 'pearson_r': sub[i, j]})
    pair_table = pd.DataFrame(pair_rows)
    return p, rank_table, pair_table


def main() -> None:
    mpl.rcParams.update({'font.family': 'Arimo', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
    plot, rank_table, pair_table = build_plot()
    outputs = {
        'pdf': OUTDIR / 'pearson_correlation_clusters_plotnine.pdf',
        'svg': OUTDIR / 'pearson_correlation_clusters_plotnine.svg',
        'png': OUTDIR / 'pearson_correlation_clusters_plotnine.png',
    }
    fig = plot.draw(show=False)
    fig.set_size_inches(12.0, 8.16, forward=True)
    fig.savefig(outputs['pdf'], format='pdf', dpi=300, facecolor='white')
    fig.savefig(outputs['svg'], format='svg', dpi=300, facecolor='white')
    fig.savefig(outputs['png'], format='png', dpi=300, facecolor='white')
    plt.close(fig)
    rank_table.to_csv(OUTDIR / 'clustered_gene_order.csv', index=False)
    pair_table.to_csv(OUTDIR / 'inset_pairwise_pearson_correlations.csv', index=False)
    print(f'Created plotnine outputs in: {OUTDIR}')
    for path in outputs.values():
        print(path)

if __name__ == '__main__':
    main()
